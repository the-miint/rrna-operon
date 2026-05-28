-- Positive 16S filter via minimap2 (map-hifi preset).
-- Keeps reads with >=1kb alignment to the reference on either strand.

CREATE OR REPLACE VIEW reads_for_filter AS
SELECT read_id, seq AS sequence1 FROM reads_unfiltered;

CREATE OR REPLACE TABLE positive_ref AS
SELECT * FROM read_fastx(getvariable('positive_ref_path'));

CREATE OR REPLACE TABLE positive_hits AS
SELECT DISTINCT read_id
FROM align_minimap2('reads_for_filter', subject_table='positive_ref',
     preset := 'map-hifi', max_secondary := 0)
WHERE NOT alignment_is_unmapped(flags)
  AND stop_position - position + 1 >= 1000;

-- Per-read pass/fail manifest: one row per input read with a boolean,
-- exported as positive_filter_status.parquet so callers can recover the
-- set of rejected reads after the working DB is deleted.
CREATE OR REPLACE TABLE positive_filter_status AS
SELECT r.read_id,
       (h.read_id IS NOT NULL) AS passed
FROM reads_unfiltered r
LEFT JOIN positive_hits h USING (read_id);

CREATE OR REPLACE TABLE reads AS
SELECT r.*
FROM reads_unfiltered r
SEMI JOIN positive_hits h USING (read_id);

DROP TABLE IF EXISTS positive_hits;
DROP TABLE IF EXISTS positive_ref;
DROP VIEW IF EXISTS reads_for_filter;
