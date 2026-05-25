-- Stage 12: hard-mask homopolymers (HP runs >= 3 → length 2)
CREATE OR REPLACE TABLE hp_masked AS
SELECT bin_id, ubs, seq AS original_seq, qual,
       regexp_replace(
         regexp_replace(
           regexp_replace(
             regexp_replace(seq, 'A{3,}', 'AA', 'g'),
                                'T{3,}', 'TT', 'g'),
                                'C{3,}', 'CC', 'g'),
                                'G{3,}', 'GG', 'g') AS seq
FROM high_cov_consensus;

-- Stage 13: two-pass cluster at subcluster_id (default 0.998)
CREATE OR REPLACE VIEW _clust1_input AS
SELECT bin_id AS read_id, seq AS sequence1
FROM hp_masked
ORDER BY ubs DESC;

CREATE OR REPLACE TABLE clust1 AS
SELECT * FROM cluster_sequences_vsearch('_clust1_input',
    id := getvariable('subcluster_id'), strand := 'both');

CREATE OR REPLACE VIEW _clust2_input AS
SELECT c.centroid_id AS read_id, h.seq AS sequence1
FROM clust1 c
JOIN hp_masked h ON c.read_id = h.bin_id
WHERE c.is_centroid;

CREATE OR REPLACE TABLE clust2 AS
SELECT * FROM cluster_sequences_vsearch('_clust2_input',
    id := getvariable('subcluster_id'), strand := 'both');

DROP VIEW IF EXISTS _clust1_input;
DROP VIEW IF EXISTS _clust2_input;

-- Bin reads into multi-member clusters
CREATE OR REPLACE TABLE cluster_members AS
SELECT c2.centroid_id AS cluster_id, c1.read_id AS bin_id
FROM clust1 c1
JOIN clust2 c2 ON c1.centroid_id = c2.read_id;

-- Empty variants table (stages 14-17 populate if multi-member clusters exist)
CREATE OR REPLACE TABLE variants(variant_id VARCHAR, support BIGINT,
                                  seq VARCHAR, qual UTINYINT[]);
