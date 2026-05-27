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
-- read_fastx-compatible schema + umi_count (molecule count) + total_reads (raw read count).
CREATE OR REPLACE VIEW export_unique AS
WITH clusters_with_variants AS (
    SELECT DISTINCT cluster_id FROM variant_bins
),
variant_reads_agg AS (
    SELECT vb.variant_id,
           vb.support AS umi_count,
           sum(h.ubs)::BIGINT AS total_reads
    FROM variant_bins vb, unnest(vb.read_ids) AS u(rid)
    JOIN high_cov_consensus h ON h.bin_id = u.rid
    GROUP BY vb.variant_id, vb.support
),
unvaried_clusters AS (
    SELECT cm.cluster_id,
           count(*)::BIGINT AS umi_count,
           sum(h.ubs)::BIGINT AS total_reads,
           any_value(h.seq) FILTER (WHERE cm.bin_id = cm.cluster_id) AS seq,
           any_value(h.qual) FILTER (WHERE cm.bin_id = cm.cluster_id) AS qual
    FROM cluster_members cm
    JOIN high_cov_consensus h USING (bin_id)
    WHERE cm.cluster_id NOT IN (SELECT cluster_id FROM clusters_with_variants)
    GROUP BY cm.cluster_id
)
SELECT read_id, seq AS sequence1, qual AS qual1, umi_count, total_reads
FROM (
    SELECT v.variant_id AS read_id, v.seq, v.qual, a.umi_count, a.total_reads
    FROM variants v
    JOIN variant_reads_agg a USING (variant_id)
    UNION ALL
    SELECT cluster_id AS read_id, seq, qual, umi_count, total_reads
    FROM unvaried_clusters
);
