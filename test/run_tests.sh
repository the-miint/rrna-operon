#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
DUCKDB="${DUCKDB:-${PROJECT_DIR}/../duckdb-miint/build/release/duckdb}"
WORKDIR="${SCRIPT_DIR}/workdir"
FIXTURE="${SCRIPT_DIR}/fixtures/synthetic.fq"
TEST_REF="${SCRIPT_DIR}/fixtures/test_ref.fa"
PARAMS="${PROJECT_DIR}/params/test.sql"

PASS=0
FAIL=0
ERRORS=""

cleanup() { rm -f "${WORKDIR}"/*.csv "${WORKDIR}"/*.parquet "${WORKDIR}"/*.duckdb "${WORKDIR}"/*.wal; }

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

# --- Setup ---
echo "=== rrna-operon test suite ==="
echo "DuckDB: $DUCKDB"
echo "Fixture: $FIXTURE"
echo

cleanup
mkdir -p "$WORKDIR"

# ============================================================================
# SESSION 1: Ingest + SortMeRNA positive filter (in-memory)
# SortMeRNA corrupts WFA2 global state, so it must run in a separate process.
# Exports passing read_ids to CSV; downstream re-reads from read_fastx.
# ============================================================================
echo "Phase 1+1b: Ingest + Positive filter"
if [[ -f "${PROJECT_DIR}/sql/00_ingest.sql" && -f "${PROJECT_DIR}/sql/05_positive_filter.sql" ]]; then
    "$DUCKDB" :memory: <<EOF
CREATE VIEW raw_input AS SELECT read_id, sequence1, qual1 FROM read_fastx('${FIXTURE}');
SET VARIABLE positive_ref_path = '${TEST_REF}';
.read ${PARAMS}
.read ${PROJECT_DIR}/sql/00_ingest.sql
.read ${PROJECT_DIR}/sql/05_positive_filter.sql
COPY (SELECT read_id FROM reads) TO '${WORKDIR}/passing_ids.csv' (HEADER, DELIMITER ',');
EOF
    assert_file_exists "passing_ids checkpoint" "${WORKDIR}/passing_ids.csv"
    reads_count=$("$DUCKDB" -csv -noheader :memory: "SELECT count(*) FROM read_csv('${WORKDIR}/passing_ids.csv')" | tr -d '[:space:]')
    assert_ge "passing reads" "1" "$reads_count"
    assert_le "passing reads <= 11" "11" "$reads_count"
else
    echo "  SKIP: sql/00_ingest.sql or sql/05_positive_filter.sql not found"
fi

# ============================================================================
# SESSIONS 2-5: UMI extraction (1 extract_linked_amplicon call per process)
# WFA2 has a non-deterministic use-after-free: >1 call per process segfaults.
# Each process re-reads from read_fastx (qual can't round-trip through parquet).
# ============================================================================
READS_BOOTSTRAP="CREATE TABLE _ids AS SELECT read_id FROM read_csv('${WORKDIR}/passing_ids.csv'); CREATE TABLE reads AS SELECT read_id, sequence1 AS seq, qual1 AS qual FROM read_fastx('${FIXTURE}') WHERE read_id IN (SELECT read_id FROM _ids); DROP TABLE _ids;"

run_single_extract() {
    local extract_sql="$1" table_name="$2" out_csv="$3"
    local attempt
    for attempt in 1 2 3 4 5; do
        "$DUCKDB" :memory: <<EOF && return 0
.read ${PARAMS}
${READS_BOOTSTRAP}
.read ${PROJECT_DIR}/sql/10a_read_ends.sql
.read ${extract_sql}
COPY ${table_name} TO '${out_csv}' (HEADER, DELIMITER ',');
EOF
        echo "  WARN: extract_linked_amplicon crashed (WFA2 bug), retry $attempt/5"
    done
    echo "  ERROR: extract_linked_amplicon failed after 5 retries"
    return 1
}

echo "Phase 2: UMI extraction (10a-10f)"
if [[ -f "${PROJECT_DIR}/sql/10f_umi_cluster.sql" ]]; then
    run_single_extract "${PROJECT_DIR}/sql/10b_fwd_u1.sql" "_fwd_u1" "${WORKDIR}/fwd_u1.csv"
    run_single_extract "${PROJECT_DIR}/sql/10c_fwd_u2.sql" "_fwd_u2" "${WORKDIR}/fwd_u2.csv"
    run_single_extract "${PROJECT_DIR}/sql/10d_rev_u1.sql" "_rev_u1" "${WORKDIR}/rev_u1.csv"
    run_single_extract "${PROJECT_DIR}/sql/10e_rev_u2.sql" "_rev_u2" "${WORKDIR}/rev_u2.csv"

    # Session 6: combine + cluster (no WFA2 calls)
    "$DUCKDB" :memory: <<EOF
.read ${PARAMS}
CREATE TABLE _fwd_u1 AS SELECT * FROM read_csv('${WORKDIR}/fwd_u1.csv');
CREATE TABLE _fwd_u2 AS SELECT * FROM read_csv('${WORKDIR}/fwd_u2.csv');
CREATE TABLE _rev_u1 AS SELECT * FROM read_csv('${WORKDIR}/rev_u1.csv');
CREATE TABLE _rev_u2 AS SELECT * FROM read_csv('${WORKDIR}/rev_u2.csv');
.read ${PROJECT_DIR}/sql/10f_umi_cluster.sql
COPY (SELECT 'umi_candidates' AS tbl, count(*) AS n FROM umi_candidates
      UNION ALL SELECT 'umi_ref', count(*) FROM umi_ref)
TO '${WORKDIR}/phase2_counts.csv' (HEADER, DELIMITER ',');
EOF
    umi_candidates_n=$("$DUCKDB" -csv -noheader :memory: "SELECT n FROM read_csv('${WORKDIR}/phase2_counts.csv') WHERE tbl='umi_candidates'" | tr -d '[:space:]')
    umi_ref_n=$("$DUCKDB" -csv -noheader :memory: "SELECT n FROM read_csv('${WORKDIR}/phase2_counts.csv') WHERE tbl='umi_ref'" | tr -d '[:space:]')
    assert_ge "umi_candidates has rows" "1" "$umi_candidates_n"
    assert_eq "umi_ref has 2 UMI bins" "2" "$umi_ref_n"
else
    echo "  SKIP: sql/10f_umi_cluster.sql not found"
fi

# ============================================================================
# Remaining phases (TODO — skip for now)
# ============================================================================
echo "Phase 3: UMI binning (20_umi_bin.sql)"
if [[ -f "${PROJECT_DIR}/sql/20_umi_bin.sql" ]]; then
    echo "  TODO: implement"
else
    echo "  SKIP: not yet implemented"
fi

echo "Phase 4: Consensus (30_consensus.sql)"
if [[ -f "${PROJECT_DIR}/sql/30_consensus.sql" ]]; then
    echo "  TODO: implement"
else
    echo "  SKIP: not yet implemented"
fi

echo "Phase 5: Variants (40_variants.sql)"
if [[ -f "${PROJECT_DIR}/sql/40_variants.sql" ]]; then
    echo "  TODO: implement"
else
    echo "  SKIP: not yet implemented"
fi

echo "Phase 6: Export (50_export.sql)"
if [[ -f "${PROJECT_DIR}/sql/50_export.sql" ]]; then
    echo "  TODO: implement"
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
