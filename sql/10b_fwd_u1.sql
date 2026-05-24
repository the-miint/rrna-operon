CREATE OR REPLACE TABLE _fwd_u1 AS
SELECT read_id,
       (extract_linked_amplicon(start_seq, start_qual,
           getvariable('fw1'), getvariable('fw2'), 18, 18, 0.10)).sequence AS u1
FROM read_ends;
