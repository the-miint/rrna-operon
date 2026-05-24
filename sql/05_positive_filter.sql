CREATE OR REPLACE VIEW reads_for_sortmerna AS
SELECT read_id, seq AS sequence1 FROM reads_unfiltered;

CREATE OR REPLACE TABLE sortmerna_hits AS
SELECT DISTINCT read_id
FROM align_sortmerna_rrna('reads_for_sortmerna',
       ref_paths := [getvariable('positive_ref_path')])
WHERE aligned = 1;

CREATE OR REPLACE TABLE reads AS
SELECT r.*
FROM reads_unfiltered r
SEMI JOIN sortmerna_hits h USING (read_id);
