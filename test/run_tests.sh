#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
DUCKDB="${DUCKDB:-${PROJECT_DIR}/../duckdb-miint/build/release/duckdb}"
WORKDIR="${SCRIPT_DIR}/workdir"
DB="${WORKDIR}/test.duckdb"
FIXTURE="${SCRIPT_DIR}/fixtures/synthetic.fq"
TEST_REF="${SCRIPT_DIR}/fixtures/test_ref.fa"
PARAMS="${PROJECT_DIR}/params/test.sql"

PASS=0
FAIL=0
ERRORS=""

cleanup() { rm -f "$DB" "${DB}.wal" "${WORKDIR}"/*.csv "${WORKDIR}"/*.parquet; }

assert_eq() {
    local desc="$1" expected="$2" actual="$3"
    if [[ "$actual" == "$expected" ]]; then
        PASS=$((PASS + 1)); echo "  PASS: $desc"
    else
        FAIL=$((FAIL + 1)); ERRORS="${ERRORS}\n  FAIL: $desc (expected='$expected', got='$actual')"
        echo "  FAIL: $desc (expected='$expected', got='$actual')"
    fi
}

assert_ge() {
    local desc="$1" min="$2" actual="$3"
    if [[ "$actual" -ge "$min" ]]; then
        PASS=$((PASS + 1)); echo "  PASS: $desc"
    else
        FAIL=$((FAIL + 1)); ERRORS="${ERRORS}\n  FAIL: $desc (expected >= $min, got='$actual')"
        echo "  FAIL: $desc (expected >= $min, got='$actual')"
    fi
}

assert_le() {
    local desc="$1" max="$2" actual="$3"
    if [[ "$actual" -le "$max" ]]; then
        PASS=$((PASS + 1)); echo "  PASS: $desc"
    else
        FAIL=$((FAIL + 1)); ERRORS="${ERRORS}\n  FAIL: $desc (expected <= $max, got='$actual')"
        echo "  FAIL: $desc (expected <= $max, got='$actual')"
    fi
}

assert_file_exists() {
    local desc="$1" path="$2"
    if [[ -f "$path" ]]; then
        PASS=$((PASS + 1)); echo "  PASS: $desc"
    else
        FAIL=$((FAIL + 1)); ERRORS="${ERRORS}\n  FAIL: $desc (file not found: $path)"
        echo "  FAIL: $desc (file not found: $path)"
    fi
}

query() {
    "$DUCKDB" -csv -noheader "$DB" "$1" 2>/dev/null | tr -d '[:space:]'
}

# Session variables that must be set each invocation (not persisted in .duckdb)
VAR_SETUP="SET VARIABLE positive_ref_path = '${TEST_REF}'; SET VARIABLE output_dir = '${WORKDIR}';"

run_phase() {
    local sql_file="$1"
    "$DUCKDB" "$DB" <<EOF
${VAR_SETUP}
.read ${PARAMS}
.read ${sql_file}
EOF
}

# --- Setup ---
echo "=== rrna-operon test suite ==="
echo "DuckDB: $DUCKDB"
echo "Fixture: $FIXTURE"
echo

cleanup
mkdir -p "$WORKDIR"

# Create the persistent raw_input view (persists in .duckdb file)
"$DUCKDB" "$DB" "CREATE OR REPLACE VIEW raw_input AS SELECT read_id, sequence1, qual1 FROM read_fastx('${FIXTURE}');"

# --- Phase 1: Ingest ---
echo "Phase 1: Ingest (00_ingest.sql)"
if [[ -f "${PROJECT_DIR}/sql/00_ingest.sql" ]]; then
    run_phase "${PROJECT_DIR}/sql/00_ingest.sql"
    assert_ge "reads_unfiltered has rows" "1" "$(query "SELECT count(*) FROM reads_unfiltered")"
    assert_eq "reads_unfiltered row count" "11" "$(query "SELECT count(*) FROM reads_unfiltered")"
    assert_eq "all reads pass length filter" "0" \
        "$(query "SELECT count(*) FROM reads_unfiltered WHERE length(seq) < 50 OR length(seq) > 300")"
else
    echo "  SKIP: sql/00_ingest.sql not found"
fi

# --- Phase 1b: Positive filter ---
echo "Phase 1b: Positive filter (05_positive_filter.sql)"
if [[ -f "${PROJECT_DIR}/sql/05_positive_filter.sql" ]]; then
    run_phase "${PROJECT_DIR}/sql/05_positive_filter.sql"
    assert_ge "reads has rows" "1" "$(query "SELECT count(*) FROM reads")"
    assert_le "reads <= reads_unfiltered" \
        "$(query "SELECT count(*) FROM reads_unfiltered")" \
        "$(query "SELECT count(*) FROM reads")"
    assert_eq "all reads come from reads_unfiltered" "0" \
        "$(query "SELECT count(*) FROM reads r WHERE NOT EXISTS (SELECT 1 FROM reads_unfiltered u WHERE u.read_id = r.read_id)")"
else
    echo "  SKIP: sql/05_positive_filter.sql not found"
fi

# --- Phase 2: UMI extraction ---
echo "Phase 2: UMI extraction (10_umi_extract.sql)"
if [[ -f "${PROJECT_DIR}/sql/10_umi_extract.sql" ]]; then
    run_phase "${PROJECT_DIR}/sql/10_umi_extract.sql"
    assert_ge "read_ends has rows" "1" "$(query "SELECT count(*) FROM read_ends")"
    assert_ge "umi_candidates has rows" "1" "$(query "SELECT count(*) FROM umi_candidates")"
    assert_eq "umi_ref has 2 UMI bins" "2" "$(query "SELECT count(*) FROM umi_ref")"
else
    echo "  SKIP: sql/10_umi_extract.sql not found"
fi

# --- Phase 3: UMI binning ---
echo "Phase 3: UMI binning (20_umi_bin.sql)"
if [[ -f "${PROJECT_DIR}/sql/20_umi_bin.sql" ]]; then
    run_phase "${PROJECT_DIR}/sql/20_umi_bin.sql"
    assert_ge "best_bin has rows" "1" "$(query "SELECT count(*) FROM best_bin")"
    assert_eq "bin_pass has 2 passing bins" "2" "$(query "SELECT count(*) FROM bin_pass")"
    assert_ge "bin_reads has rows" "1" "$(query "SELECT count(*) FROM bin_reads")"
else
    echo "  SKIP: not yet implemented"
fi

# --- Phase 4: Consensus ---
echo "Phase 4: Consensus (30_consensus.sql)"
if [[ -f "${PROJECT_DIR}/sql/30_consensus.sql" ]]; then
    run_phase "${PROJECT_DIR}/sql/30_consensus.sql"
    assert_eq "bin_consensus has 2 rows" "2" "$(query "SELECT count(*) FROM bin_consensus")"
    assert_eq "consensus seqs non-empty" "0" \
        "$(query "SELECT count(*) FROM bin_consensus WHERE length(consensus_seq) = 0")"
    assert_ge "high_cov_consensus has rows" "1" "$(query "SELECT count(*) FROM high_cov_consensus")"
else
    echo "  SKIP: not yet implemented"
fi

# --- Phase 5: Variants ---
echo "Phase 5: Variants (40_variants.sql)"
if [[ -f "${PROJECT_DIR}/sql/40_variants.sql" ]]; then
    run_phase "${PROJECT_DIR}/sql/40_variants.sql"
    assert_ge "hp_masked has rows" "1" "$(query "SELECT count(*) FROM hp_masked")"
    # variants may be empty with few UMI bins (no multi-member clusters)
    assert_eq "variants table exists" "0" "$(query "SELECT count(*) FROM variants WHERE false")"
else
    echo "  SKIP: not yet implemented"
fi

# --- Phase 6: Export ---
echo "Phase 6: Export (50_export.sql)"
if [[ -f "${PROJECT_DIR}/sql/50_export.sql" ]]; then
    run_phase "${PROJECT_DIR}/sql/50_export.sql"
    assert_file_exists "consensus.parquet" "${WORKDIR}/consensus.parquet"
    assert_file_exists "variants.parquet" "${WORKDIR}/variants.parquet"
    assert_file_exists "consensus.fa" "${WORKDIR}/consensus.fa"
    assert_file_exists "variants.fa" "${WORKDIR}/variants.fa"
else
    echo "  SKIP: not yet implemented"
fi

# --- Summary ---
echo
echo "=== Results: $PASS passed, $FAIL failed ==="
if [[ $FAIL -gt 0 ]]; then
    echo -e "\nFailures:$ERRORS"
    exit 1
fi
