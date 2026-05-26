-- Prepare export views with read_fastx-compatible schema
CREATE OR REPLACE VIEW export_consensus AS
SELECT row_number() OVER () AS sequence_index,
       bin_id AS read_id,
       'ubs=' || ubs AS comment,
       seq AS sequence1,
       NULL::VARCHAR AS sequence2,
       qual AS qual1,
       NULL::UTINYINT[] AS qual2
FROM high_cov_consensus;

CREATE OR REPLACE VIEW export_consensus_fasta AS
SELECT bin_id AS read_id, seq AS sequence1
FROM high_cov_consensus;

CREATE OR REPLACE VIEW export_variants AS
SELECT row_number() OVER () AS sequence_index,
       variant_id AS read_id,
       'support=' || support AS comment,
       seq AS sequence1,
       NULL::VARCHAR AS sequence2,
       qual AS qual1,
       NULL::UTINYINT[] AS qual2
FROM variants;

CREATE OR REPLACE VIEW export_variants_fasta AS
SELECT variant_id AS read_id, seq AS sequence1
FROM variants
WHERE seq IS NOT NULL;

-- Deduplicated unique sequences with UMI-bin abundance.
-- Clusters with variants: one row per phased variant, support = bin count.
-- Clusters without variants: one row for the centroid consensus, support = cluster size.
CREATE OR REPLACE VIEW export_unique AS
WITH clusters_with_variants AS (
    SELECT DISTINCT cluster_id FROM variant_bins
),
unvaried_clusters AS (
    SELECT cm.cluster_id,
           count(*) AS support,
           any_value(h.seq) FILTER (WHERE cm.bin_id = cm.cluster_id) AS seq,
           any_value(h.qual) FILTER (WHERE cm.bin_id = cm.cluster_id) AS qual
    FROM cluster_members cm
    JOIN high_cov_consensus h USING (bin_id)
    WHERE cm.cluster_id NOT IN (SELECT cluster_id FROM clusters_with_variants)
    GROUP BY cm.cluster_id
)
SELECT row_number() OVER () AS sequence_index,
       read_id,
       'support=' || support AS comment,
       seq AS sequence1,
       NULL::VARCHAR AS sequence2,
       qual AS qual1,
       NULL::UTINYINT[] AS qual2
FROM (
    SELECT variant_id AS read_id, support, seq, qual FROM variants
    UNION ALL
    SELECT cluster_id AS read_id, support, seq, qual FROM unvaried_clusters
);

CREATE OR REPLACE VIEW export_unique_fasta AS
SELECT read_id, sequence1 FROM export_unique;
