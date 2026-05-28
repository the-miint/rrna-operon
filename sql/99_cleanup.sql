-- Drop intermediate tables to shrink the pipeline DB. Run after exports.
--
-- Three tiers:
--   1. Always droppable: pure intermediates, no information loss.
--   2. Per-read tables (gated by VARIABLE drop_per_read = true): hold a row
--      per input read; large but rebuildable only by re-running the pipeline.
--   3. Final retained tables: kept for downstream re-analysis (documented in
--      README.md "Database schema" section).
--
-- After dropping, CHECKPOINT compacts the file on disk. Note that DuckDB may
-- still keep some slack; for maximum compaction, EXPORT/IMPORT the database
-- to a new file.

-- ── Drop transient views that reference intermediate tables ──────────────
DROP VIEW IF EXISTS umi_candidates_normalized;
DROP VIEW IF EXISTS umi_unique_for_cluster;
DROP VIEW IF EXISTS mafft_input;

-- ── Tier 1: always drop ───────────────────────────────────────────────────
DROP TABLE IF EXISTS hits1;
DROP TABLE IF EXISTS hits2;
DROP TABLE IF EXISTS pileup;
DROP TABLE IF EXISTS msa_out;
DROP TABLE IF EXISTS cluster_alignments;
DROP TABLE IF EXISTS variant_msa;
DROP TABLE IF EXISTS variant_reads;
DROP TABLE IF EXISTS read_signatures;
DROP TABLE IF EXISTS snp_sites;
DROP TABLE IF EXISTS read_ends;
DROP TABLE IF EXISTS amplicon;
DROP TABLE IF EXISTS hp_masked;
DROP TABLE IF EXISTS bin_consensus;
DROP TABLE IF EXISTS clust1;
DROP TABLE IF EXISTS clust2;
DROP TABLE IF EXISTS umi_unique;
DROP TABLE IF EXISTS umi_clusters;
DROP TABLE IF EXISTS umi_halves;
DROP TABLE IF EXISTS half_owner;
DROP TABLE IF EXISTS umi_candidates;
DROP TABLE IF EXISTS best_bin;

-- ── Tier 2: per-read tables (opt-in) ──────────────────────────────────────
-- Drop these to save the most disk. They are large (rows-per-input-read)
-- and contain the full sequence + quality, so rebuilding them requires
-- re-running 00_ingest, 05_positive_filter, 10_umi_extract, 20_umi_bin.
SET VARIABLE drop_per_read = coalesce(try(getvariable('drop_per_read')), false);

CREATE OR REPLACE TABLE _cleanup_decision AS
SELECT getvariable('drop_per_read')::BOOLEAN AS drop_per_read;

-- Execute the conditional drops by issuing them as separate statements
-- guarded with WHERE; DuckDB does not support IF in DDL, so we materialize
-- a sentinel and use CREATE OR REPLACE TABLE … AS to keep or empty the table.
CREATE OR REPLACE TABLE reads_unfiltered AS
SELECT * FROM reads_unfiltered
WHERE NOT (SELECT drop_per_read FROM _cleanup_decision);

CREATE OR REPLACE TABLE reads AS
SELECT * FROM reads
WHERE NOT (SELECT drop_per_read FROM _cleanup_decision);

CREATE OR REPLACE TABLE bin_reads AS
SELECT * FROM bin_reads
WHERE NOT (SELECT drop_per_read FROM _cleanup_decision);

DROP TABLE _cleanup_decision;

-- ── Compact ───────────────────────────────────────────────────────────────
CHECKPOINT;
