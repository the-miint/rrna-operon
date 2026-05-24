-- Stage 6: UMI→read mapping via Hamming distance
-- Uses canonical umi_candidates from stage 3 (both anchors validated).
-- extract_linked_amplicon doesn't support single-anchor mode, so the SKETCH's
-- "looser re-extraction" is replaced by direct use of validated candidates.

-- Normalize UMI halves: - strand reads have swapped/RC'd halves relative to
-- the + strand form stored in umi_ref. Swap and RC to put all in + orientation.
CREATE OR REPLACE VIEW umi_candidates_normalized AS
SELECT read_id, strand,
       CASE WHEN strand = '+' THEN u1
            ELSE sequence_dna_reverse_complement(u2) END AS u1,
       CASE WHEN strand = '+' THEN u2
            ELSE sequence_dna_reverse_complement(u1) END AS u2
FROM umi_candidates;

-- Adapter views for match_short_barcodes (requires columns: id, sequence)
CREATE OR REPLACE VIEW _query_u1 AS
SELECT read_id || '|' || strand AS id, u1 AS sequence FROM umi_candidates_normalized;

CREATE OR REPLACE VIEW _ref_u1 AS
SELECT umi_id AS id, u1 AS sequence FROM umi_ref;

CREATE OR REPLACE VIEW _query_u2 AS
SELECT read_id || '|' || strand AS id, u2 AS sequence FROM umi_candidates_normalized;

CREATE OR REPLACE VIEW _ref_u2 AS
SELECT umi_id AS id, u2 AS sequence FROM umi_ref;

-- Hamming match: each read terminus UMI half vs every UMI cluster half
CREATE OR REPLACE TABLE hits1 AS
SELECT * FROM match_short_barcodes('_query_u1', '_ref_u1',
    max_nm := getvariable('umi_max_nm_per_half'), report_all := true);

CREATE OR REPLACE TABLE hits2 AS
SELECT * FROM match_short_barcodes('_query_u2', '_ref_u2',
    max_nm := getvariable('umi_max_nm_per_half'), report_all := true);

DROP VIEW IF EXISTS _query_u1;
DROP VIEW IF EXISTS _ref_u1;
DROP VIEW IF EXISTS _query_u2;
DROP VIEW IF EXISTS _ref_u2;

-- Stage 7: pick best UMI per (read, strand)
CREATE OR REPLACE TABLE best_bin AS
WITH pairs AS (
    SELECT split_part(h1.query_id, '|', 1) AS read_id,
           split_part(h1.query_id, '|', 2) AS strand,
           h1.ref_id                        AS umi_id,
           h1.nm + h2.nm                    AS nm_total
    FROM hits1 h1
    JOIN hits2 h2 USING (query_id, ref_id)
    WHERE h1.nm + h2.nm <= 2 * getvariable('umi_max_nm_per_half')
),
ranked AS (
    SELECT *, row_number() OVER (PARTITION BY read_id
                                 ORDER BY nm_total, umi_id) AS rk
    FROM pairs
)
SELECT read_id, strand, umi_id, nm_total
FROM ranked WHERE rk = 1;

-- Stage 8: bin-level filters (RoF, UME, BCR)
CREATE OR REPLACE TABLE bin_pass AS
WITH bin_stats AS (
    SELECT umi_id,
           count(*)                             AS n,
           count(*) FILTER (WHERE strand = '+') AS n_plus,
           count(*) FILTER (WHERE strand = '-') AS n_neg,
           avg(nm_total)                        AS ume_mean,
           stddev_samp(nm_total)                AS ume_sd
    FROM best_bin
    GROUP BY umi_id
),
joined AS (
    SELECT b.*,
           u.size AS cluster_size,
           b.n::DOUBLE / u.size AS bcr
    FROM bin_stats b
    JOIN umi_ref u USING (umi_id)
)
SELECT *
FROM joined
WHERE n        >= getvariable('min_bin_size')
  AND n_plus   >= 2 AND n_neg >= 2
  AND least(n_plus, n_neg)::DOUBLE / n >= getvariable('ro_frac')
  AND ume_mean <= getvariable('ume_mean_max')
  AND (ume_sd IS NULL OR ume_sd <= getvariable('ume_sd_max'))
  AND bcr      <= getvariable('bin_cluster_ratio');

-- Final read→bin assignment (only reads in passing bins)
CREATE OR REPLACE TABLE bin_reads AS
SELECT bb.umi_id AS bin_id, bb.read_id, bb.strand,
       r.seq, r.qual
FROM best_bin bb
JOIN bin_pass bp USING (umi_id)
JOIN reads r USING (read_id);
