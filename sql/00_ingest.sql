-- filter_read(seq, qual, min_length, max_length, qualified_q, max_unqualified_pct, max_n, min_avg_q)
CREATE OR REPLACE TABLE reads_unfiltered AS
SELECT read_id, sequence1 AS seq, qual1 AS qual
FROM (
    SELECT read_id, sequence1, qual1,
           filter_read(sequence1, qual1,
                       getvariable('min_len'),
                       getvariable('max_len'),
                       15, 40, 5,
                       getvariable('min_q')) AS f
    FROM raw_input
)
WHERE (f).fail_reason IS NULL;
