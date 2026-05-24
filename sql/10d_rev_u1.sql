CREATE OR REPLACE TABLE _rev_u1 AS
SELECT read_id,
       (extract_linked_amplicon(start_seq, start_qual,
           getvariable('rv1'), getvariable('rv2'), 18, 18, 0.10)).sequence AS u1
FROM read_ends;
