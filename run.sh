#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DUCKDB="${DUCKDB:-${SCRIPT_DIR}/../duckdb-miint/build/release/duckdb}"
PARAMS=""
OUTPUT_DIR=""

usage() {
    echo "Usage: $0 [--params FILE] [--output DIR] [--duckdb PATH] INPUT"
    echo
    echo "  INPUT         fastq(.gz) or parquet file"
    echo "  --params      parameter file (default: params/revio.sql)"
    echo "  --output      output directory (default: ./output)"
    echo "  --duckdb      path to duckdb binary"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --params)  PARAMS="$2"; shift 2 ;;
        --output)  OUTPUT_DIR="$2"; shift 2 ;;
        --duckdb)  DUCKDB="$2"; shift 2 ;;
        -h|--help) usage ;;
        -*)        echo "Unknown option: $1"; usage ;;
        *)         INPUT="$1"; shift ;;
    esac
done

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
    echo "--- ${phase_name} ---"
    "$DUCKDB" "$DB" <<EOF
${INPUT_SQL}
SET VARIABLE output_dir = '${OUTPUT_DIR}';
.read ${PARAMS}
.read ${sql_file}
EOF
}

echo "=== rrna-operon pipeline ==="
echo "Input:  $INPUT"
echo "Params: $PARAMS"
echo "Output: $OUTPUT_DIR"
echo "DB:     $DB"
echo

for sql_file in "${SCRIPT_DIR}"/sql/[0-9]*.sql; do
    run_phase "$sql_file"
done

echo
echo "=== Done ==="
echo "Outputs in: $OUTPUT_DIR"
