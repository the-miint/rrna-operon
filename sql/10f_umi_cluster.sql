-- Stage 3 (combine): merge fwd/rev extractions into umi_candidates
CREATE OR REPLACE TABLE umi_candidates AS
WITH fwd AS (
    SELECT a.read_id, '+' AS strand, a.u1, b.u2
    FROM _fwd_u1 a JOIN _fwd_u2 b USING (read_id)
),
rev AS (
    SELECT a.read_id, '-' AS strand, a.u1, b.u2
    FROM _rev_u1 a JOIN _rev_u2 b USING (read_id)
)
SELECT read_id, strand, u1, u2, u1 || u2 AS umi_pair
FROM (SELECT * FROM fwd UNION ALL SELECT * FROM rev)
WHERE u1 IS NOT NULL AND u1 != ''
  AND u2 IS NOT NULL AND u2 != ''
  AND regexp_matches(u1 || u2, getvariable('umi_pair_pattern'));

DROP TABLE IF EXISTS _fwd_u1;
DROP TABLE IF EXISTS _fwd_u2;
DROP TABLE IF EXISTS _rev_u1;
DROP TABLE IF EXISTS _rev_u2;

-- Stage 4: dereplicate + cluster UMI pairs
CREATE OR REPLACE TABLE umi_unique AS
SELECT umi_pair AS sequence,
       count(*)                                            AS size_in,
       'u_' || row_number() OVER (ORDER BY count(*) DESC) AS umi_id
FROM umi_candidates
GROUP BY umi_pair;

CREATE OR REPLACE VIEW umi_unique_for_cluster AS
SELECT umi_id AS read_id, sequence AS sequence1
FROM umi_unique
ORDER BY size_in DESC;

CREATE OR REPLACE TABLE umi_clusters AS
SELECT read_id, is_centroid, centroid_id, identity
FROM cluster_sequences_vsearch('umi_unique_for_cluster',
     id := getvariable('umi_cluster_id'), strand := 'both');

-- Stage 5: chimera filter (each 18-mer half + RC must be owned by this UMI)
CREATE OR REPLACE TABLE umi_halves AS
SELECT u.umi_id,
       c.centroid_id,
       u.size_in AS size,
       substr(u.sequence,  1, 18) AS u1,
       substr(u.sequence, 19, 18) AS u2,
       sequence_dna_reverse_complement(substr(u.sequence,  1, 18)) AS u1rc,
       sequence_dna_reverse_complement(substr(u.sequence, 19, 18)) AS u2rc
FROM umi_unique u
JOIN umi_clusters c ON c.read_id = u.umi_id
WHERE c.is_centroid;

CREATE OR REPLACE TABLE half_owner AS
WITH all_halves AS (
    SELECT u1   AS half, umi_id, size FROM umi_halves UNION ALL
    SELECT u2  , umi_id, size FROM umi_halves UNION ALL
    SELECT u1rc, umi_id, size FROM umi_halves UNION ALL
    SELECT u2rc, umi_id, size FROM umi_halves
)
SELECT half, arg_max(umi_id, size) AS owner
FROM all_halves
GROUP BY half;

CREATE OR REPLACE TABLE umi_ref AS
SELECT h.umi_id, h.size, h.u1 || h.u2 AS umi_seq, h.u1, h.u2
FROM umi_halves h
JOIN half_owner o1   ON o1.half = h.u1   AND o1.owner = h.umi_id
JOIN half_owner o2   ON o2.half = h.u2   AND o2.owner = h.umi_id
JOIN half_owner o1rc ON o1rc.half = h.u1rc AND o1rc.owner = h.umi_id
JOIN half_owner o2rc ON o2rc.half = h.u2rc AND o2rc.owner = h.umi_id
WHERE h.u1 != h.u1rc;
