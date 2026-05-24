-- Stage 9: per-bin MSA + Q-aware consensus

-- Orient all reads to + strand before MSA (- strand reads are RC'd)
CREATE OR REPLACE VIEW mafft_input AS
SELECT bin_id, read_id,
       CASE WHEN strand = '-' THEN sequence_dna_reverse_complement(seq)
            ELSE seq END AS sequence1
FROM bin_reads;

CREATE OR REPLACE TABLE msa_out AS
SELECT * FROM align_mafft('mafft_input', sample_id := 'bin_id');

-- Join MAFFT aligned_sequence back to oriented qual, aggregate per bin
CREATE OR REPLACE TABLE bin_consensus AS
WITH consensus_raw AS (
    SELECT m.bin_id,
           compute_msa_consensus(m.aligned_sequence,
               CASE WHEN b.strand = '-' THEN list_reverse(b.qual)
                    ELSE b.qual END) AS c,
           count(*) AS ubs
    FROM msa_out m
    JOIN bin_reads b USING (bin_id, read_id)
    GROUP BY m.bin_id
)
SELECT bin_id, (c).seq AS consensus_seq, (c).qual AS consensus_qual, ubs
FROM consensus_raw;

-- Stage 10: amplicon primer trim (both orientations)
CREATE OR REPLACE TABLE amplicon AS
WITH fwd_trim AS (
    SELECT bin_id, ubs,
           extract_linked_amplicon(consensus_seq, consensus_qual,
               getvariable('fw2'),
               getvariable('rv2_rc'),
               getvariable('min_len'), getvariable('max_len'), 0.10) AS t
    FROM bin_consensus
),
rc_trim AS (
    SELECT bin_id, ubs,
           extract_linked_amplicon(
               sequence_dna_reverse_complement(consensus_seq),
               list_reverse(consensus_qual),
               getvariable('fw2'),
               getvariable('rv2_rc'),
               getvariable('min_len'), getvariable('max_len'), 0.10) AS t
    FROM bin_consensus
    WHERE bin_id NOT IN (SELECT bin_id FROM fwd_trim WHERE (t).sequence IS NOT NULL)
)
SELECT bin_id, (t).sequence AS seq, (t).qual AS qual, ubs
FROM (SELECT * FROM fwd_trim UNION ALL SELECT * FROM rc_trim)
WHERE (t).sequence IS NOT NULL;

-- Stage 11: coverage gate
CREATE OR REPLACE TABLE high_cov_consensus AS
SELECT * FROM amplicon WHERE ubs >= getvariable('umi_coverage_min');
