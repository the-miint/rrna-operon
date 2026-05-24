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
