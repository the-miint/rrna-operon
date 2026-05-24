CREATE OR REPLACE TABLE _fwd_u2 AS
SELECT read_id,
       (extract_linked_amplicon(end_seq, end_qual,
           getvariable('rv2_rc'), getvariable('rv1_rc'), 18, 18, 0.10)).sequence AS u2
FROM read_ends;
