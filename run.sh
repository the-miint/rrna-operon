#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DUCKDB="${DUCKDB:-${SCRIPT_DIR}/../duckdb-miint/build/release/duckdb}"
PARAMS=""
OUTPUT_DIR=""
CLEANUP="none"

# Load the miint extension at the start of every duckdb session. -unsigned
# is required when miint was installed from a non-default repository
# (e.g. `INSTALL miint FROM 'https://ftp.microbio.me/pub/miint'`). Both
# flags are harmless on a static-linked custom duckdb build.
duckdb() { "$DUCKDB" -unsigned -cmd "LOAD miint;" "$@"; }

usage() {
    echo "Usage: $0 [--params FILE] [--output DIR] [--duckdb PATH] [--cleanup MODE] INPUT"
    echo
    echo "  INPUT             fastq(.gz) or parquet file"
    echo "  --params          parameter file (default: params/revio.sql)"
    echo "  --output          output directory (default: ./output)"
    echo "  --duckdb          path to duckdb binary"
    echo "  --cleanup MODE    post-run DB pruning. MODE is one of:"
    echo "                      none        keep all tables (default)"
    echo "                      safe        drop pure-intermediate tables only"
    echo "                      aggressive  also drop per-read tables (reads,"
    echo "                                  reads_unfiltered, bin_reads)"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --params)  PARAMS="$2"; shift 2 ;;
        --output)  OUTPUT_DIR="$2"; shift 2 ;;
        --duckdb)  DUCKDB="$2"; shift 2 ;;
        --cleanup) CLEANUP="$2"; shift 2 ;;
        -h|--help) usage ;;
        -*)        echo "Unknown option: $1"; usage ;;
        *)         INPUT="$1"; shift ;;
    esac
done

case "$CLEANUP" in
    none|safe|aggressive) ;;
    *) echo "Error: --cleanup must be one of: none, safe, aggressive"; exit 1 ;;
esac

if [[ -z "${INPUT:-}" ]]; then
    echo "Error: INPUT file required"
    usage
fi

if [[ ! -f "$INPUT" ]]; then
    echo "Error: input file not found: $INPUT"
    exit 1
fi

INPUT="$(cd "$(dirname "$INPUT")" && pwd)/$(basename "$INPUT")"
PARAMS="${PARAMS:-${SCRIPT_DIR}/params/revio.sql}"
OUTPUT_DIR="${OUTPUT_DIR:-${SCRIPT_DIR}/output}"
mkdir -p "$OUTPUT_DIR"

DB="${OUTPUT_DIR}/pipeline.duckdb"

case "$INPUT" in
    *.fq|*.fastq|*.fq.gz|*.fastq.gz)
        INPUT_SQL="CREATE OR REPLACE VIEW raw_input AS SELECT read_id, sequence1, qual1 FROM read_fastx('${INPUT}');"
        ;;
    *.parquet)
        INPUT_SQL="CREATE OR REPLACE VIEW raw_input AS SELECT read_id, sequence1, qual1 FROM read_parquet('${INPUT}');"
        ;;
    *)
        echo "Error: unrecognized input format: $INPUT"
        echo "Expected .fq, .fastq, .fq.gz, .fastq.gz, or .parquet"
        exit 1
        ;;
esac

run_phase() {
    local sql_file="$1"
    local phase_name
    phase_name="$(basename "$sql_file" .sql)"
    local t0=$SECONDS
    echo "--- ${phase_name} ---"
    duckdb "$DB" <<EOF
${INPUT_SQL}
SET VARIABLE output_dir = '${OUTPUT_DIR}';
.read ${PARAMS}
.read ${sql_file}
EOF
    echo "    (${phase_name}: $(( SECONDS - t0 ))s)"
}

echo "=== rrna-operon pipeline ==="
echo "Input:  $INPUT"
echo "Params: $PARAMS"
echo "Output: $OUTPUT_DIR"
echo "DB:     $DB"
echo

for sql_file in "${SCRIPT_DIR}"/sql/[0-9]*.sql; do
    if [[ "$(basename "$sql_file")" == "99_cleanup.sql" ]]; then
        continue   # cleanup runs after the export step
    fi
    if [[ "$(basename "$sql_file")" == "45_variant_calling.sql" ]]; then
        has_multi=$(duckdb "$DB" -noheader -list "SELECT count(*) FROM (SELECT cluster_id FROM cluster_members GROUP BY 1 HAVING count(*) > 1);" 2>/dev/null)
        if [[ "${has_multi:-0}" -gt 0 ]]; then
            run_phase "$sql_file"
        else
            echo "--- 45_variant_calling (skipped, no multi-member clusters) ---"
        fi
        continue
    fi
    if [[ "$(basename "$sql_file")" == "05_positive_filter.sql" ]]; then
        if grep -q 'positive_ref_path' "$PARAMS" && ! grep -q '^\s*--.*positive_ref_path' "$PARAMS"; then
            run_phase "$sql_file"
        else
            echo "--- 05_positive_filter (skipped, no positive_ref_path) ---"
            duckdb "$DB" <<EOF
${INPUT_SQL}
CREATE OR REPLACE TABLE reads AS SELECT * FROM reads_unfiltered;
EOF
        fi
        continue
    fi
    run_phase "$sql_file"
done

echo "--- export ---"
duckdb "$DB" <<EOF
COPY export_consensus TO '${OUTPUT_DIR}/consensus.parquet' (FORMAT PARQUET, COMPRESSION 'zstd');
COPY export_consensus_fasta TO '${OUTPUT_DIR}/consensus.fa' (FORMAT FASTA);
COPY export_variants TO '${OUTPUT_DIR}/variants.parquet' (FORMAT PARQUET, COMPRESSION 'zstd');
COPY export_variants_fasta TO '${OUTPUT_DIR}/variants.fa' (FORMAT FASTA);
COPY export_unique TO '${OUTPUT_DIR}/unique.parquet' (FORMAT PARQUET, COMPRESSION 'zstd');
EOF

if [[ "$CLEANUP" != "none" ]]; then
    t0=$SECONDS
    echo "--- cleanup (${CLEANUP}) ---"
    drop_per_read=$([[ "$CLEANUP" == "aggressive" ]] && echo true || echo false)
    duckdb "$DB" <<EOF
SET VARIABLE drop_per_read = ${drop_per_read};
.read ${SCRIPT_DIR}/sql/99_cleanup.sql
EOF
    echo "    (cleanup: $(( SECONDS - t0 ))s)"
fi

echo
echo "=== Done ==="
echo "Outputs in: $OUTPUT_DIR"
