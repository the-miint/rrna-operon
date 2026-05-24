CREATE OR REPLACE TABLE _rev_u2 AS
SELECT read_id,
       (extract_linked_amplicon(end_seq, end_qual,
           getvariable('fw2_rc'), getvariable('fw1_rc'), 18, 18, 0.10)).sequence AS u2
FROM read_ends;
