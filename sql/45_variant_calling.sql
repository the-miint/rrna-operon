-- Stages 14-17: variant calling on multi-member clusters.
-- Only run when cluster_members has groups >= variant_min_support.

-- Stage 14: align cluster members to centroid
CREATE OR REPLACE TABLE _variant_queries AS
SELECT cm.bin_id AS read_id, h.seq AS sequence1
FROM cluster_members cm
JOIN high_cov_consensus h USING (bin_id)
WHERE cm.cluster_id IN (
    SELECT cluster_id FROM cluster_members
    GROUP BY 1 HAVING count(*) >= getvariable('variant_min_support')
);

CREATE OR REPLACE TABLE _variant_refs AS
SELECT DISTINCT cm.cluster_id AS read_id, h.seq AS sequence1
FROM cluster_members cm
JOIN high_cov_consensus h ON h.bin_id = cm.cluster_id
WHERE cm.cluster_id = cm.bin_id
AND cm.cluster_id IN (
    SELECT cluster_id FROM cluster_members
    GROUP BY 1 HAVING count(*) >= getvariable('variant_min_support')
);

CREATE OR REPLACE TABLE cluster_alignments AS
SELECT a.read_id, a.reference, a.position, a.cigar,
       q.sequence1 AS sequence, h.qual
FROM align_minimap2('_variant_queries', subject_table='_variant_refs',
     preset := 'asm5', max_secondary := 0, eqx := true) a
JOIN _variant_queries q ON q.read_id = a.read_id
JOIN high_cov_consensus h ON h.bin_id = a.read_id
WHERE a.flags & 4 = 0;

DROP TABLE IF EXISTS _variant_queries;
DROP TABLE IF EXISTS _variant_refs;

-- Stage 15: per-base pileup, find SNP positions
CREATE OR REPLACE TABLE _pileup_refs AS
SELECT DISTINCT cm.cluster_id AS ref_id, h.seq AS sequence
FROM cluster_members cm
JOIN high_cov_consensus h ON h.bin_id = cm.cluster_id
WHERE cm.cluster_id = cm.bin_id;

CREATE OR REPLACE TABLE pileup AS
SELECT * FROM compute_pileup('cluster_alignments', '_pileup_refs')
WHERE insert_pos = 0;

DROP TABLE IF EXISTS _pileup_refs;

CREATE OR REPLACE TABLE snp_sites AS
SELECT ref_id AS cluster_id, ref_pos,
       count(*) FILTER (WHERE query_base IS NOT NULL AND query_base != ref_base) AS alt_depth
FROM pileup
GROUP BY ref_id, ref_pos
HAVING alt_depth >= getvariable('snp_min_alt_depth');

-- Stage 16: per-read SNP signature → phasing
CREATE OR REPLACE TABLE read_signatures AS
SELECT p.ref_id AS cluster_id, p.read_id,
       string_agg(p.query_base, '' ORDER BY p.ref_pos) AS signature
FROM pileup p
JOIN snp_sites s ON s.cluster_id = p.ref_id AND s.ref_pos = p.ref_pos
WHERE p.query_base IS NOT NULL
GROUP BY p.ref_id, p.read_id;

CREATE OR REPLACE TABLE variant_bins AS
SELECT cluster_id, signature,
       count(*) AS support,
       cluster_id || '_var' || row_number() OVER (
           PARTITION BY cluster_id ORDER BY count(*) DESC
       ) AS variant_id,
       list(read_id) AS read_ids
FROM read_signatures
GROUP BY cluster_id, signature
HAVING count(*) >= getvariable('variant_min_support');

-- Stage 17: per-variant phased consensus via MAFFT + Q-aware voting
CREATE OR REPLACE TABLE variant_reads AS
SELECT vb.variant_id, u.read_id_in AS read_id,
       h.seq AS sequence1, h.qual
FROM variant_bins vb, unnest(vb.read_ids) AS u(read_id_in)
JOIN high_cov_consensus h ON h.bin_id = u.read_id_in;

CREATE OR REPLACE TABLE variant_msa AS
SELECT * FROM align_mafft('variant_reads', sample_id := 'variant_id');

CREATE OR REPLACE TABLE variants AS
SELECT m.variant_id,
       any_value(vb.support) AS support,
       (compute_msa_consensus(m.aligned_sequence, vr.qual)).seq AS seq,
       (compute_msa_consensus(m.aligned_sequence, vr.qual)).qual AS qual
FROM variant_msa m
JOIN variant_reads vr USING (variant_id, read_id)
JOIN variant_bins vb USING (variant_id)
GROUP BY m.variant_id;
