CREATE OR REPLACE TABLE read_ends AS
SELECT
    read_id,
    substr(seq, 1, least(200, length(seq)))              AS start_seq,
    qual[1:least(200, length(qual))]                     AS start_qual,
    substr(seq, greatest(1, length(seq) - 199))          AS end_seq,
    qual[greatest(1, length(qual) - 199):]               AS end_qual
FROM reads;
