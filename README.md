# rrna-operon

A SQL-driven pipeline for full-length rRNA operon consensus and intra-genomic
variant calling from PacBio Revio HiFi reads with UMI-tagged amplicons (Karst
et al. protocol). The entire pipeline is bash + DuckDB SQL; all bioinformatic
primitives — read I/O, primer trimming, multiple sequence alignment, pileup,
clustering — are exposed as DuckDB table/scalar functions by the
[`duckdb-miint`](https://github.com/biocore/duckdb-miint) extension. No
external aligners, polishers, or scripting glue.

## Contents

- [What it does](#what-it-does)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Inputs](#inputs)
- [Pipeline architecture](#pipeline-architecture)
- [Parameters](#parameters)
- [Outputs](#outputs)
- [Database schema](#database-schema)
- [Cleanup](#cleanup)
- [Testing](#testing)
- [Approximate runtime](#approximate-runtime)
- [Algorithmic notes](#algorithmic-notes)
- [License](#license)

## What it does

Given a FASTQ of UMI-tagged HiFi reads spanning the bacterial rRNA operon
(16S–ITS–23S, ~4.4 kb), the pipeline produces:

1. **Per-UMI consensus sequences** — one error-corrected sequence per molecule
   (UMI bin).
2. **Phased intra-genomic variants** — for any UMI cluster with multiple
   members, per-bin SNP signatures are phased into distinct operon copies.
3. **Deduplicated unique sequences** — collapsed across UMI bins, with both
   molecule count (UMI bins) and raw read count.

Output formats: Parquet (analysis-ready columnar) and FASTA (for downstream
tools). The intermediate state is preserved in a single DuckDB file so any
stage can be re-queried without re-running.

## Requirements

- Linux or macOS, 8+ cores recommended
- DuckDB CLI v1.5.3 (the version `miint` is built against)
- The `miint` DuckDB extension. Two ways to obtain it:
  1. **Install from the microbio.me extension repository** (recommended):
     ```bash
     duckdb -unsigned -c \
       "INSTALL miint FROM 'https://ftp.microbio.me/pub/miint'; LOAD miint;"
     ```
     `-unsigned` is required because the extension is hosted outside the
     official community repository. After install, every duckdb session
     that uses miint must be started with `-unsigned` (or the equivalent
     `SET allow_unsigned_extensions = true;` at session open). `run.sh`
     and `test/run_tests.sh` handle this automatically.
  2. **Use a custom-built static binary** that already has miint linked in
     (e.g. `duckdb-miint/build/release/duckdb`). Pass the path with
     `--duckdb /path/to/duckdb` or set the `DUCKDB` env var.
- ~30 GB free disk for a 1.5 M-read input (intermediate state; cleanup
  reclaims most of this; see [Cleanup](#cleanup))

The pipeline calls only built-in or `duckdb-miint`-exported functions:

| Function | Source |
|---|---|
| `read_fastx`, `read_parquet` | duckdb-miint / DuckDB core |
| `filter_read`, `extract_linked_amplicon` | duckdb-miint scalar |
| `match_short_barcodes` | duckdb-miint (Hamming, ≤32 bp) |
| `align_minimap2`, `align_abpoa` | duckdb-miint table fns |
| `compute_msa_consensus`, `compute_pileup` | duckdb-miint |
| `cluster_sequences_vsearch`, `search_sequences_vsearch` | duckdb-miint (vsearch lib) |
| `sequence_dna_reverse_complement` | duckdb-miint |

## Quick start

```bash
# 1) Run the pipeline on a FASTQ
./run.sh \
    --params params/revio.sql \
    --output output/ \
    --duckdb /path/to/duckdb-miint/build/release/duckdb \
    reads.fastq.gz

# 2) Shrink the database after exports
./run.sh ... --cleanup safe          # drops pure intermediates
./run.sh ... --cleanup aggressive    # also drops large per-read tables
```

Outputs land in `output/`:

- `consensus.parquet` / `consensus.fa` — per-UMI consensus sequences
- `variants.parquet` / `variants.fa` — phased variant consensuses
- `unique.parquet` — deduplicated unique sequences with abundance
- `pipeline.duckdb` — the working database (queryable post-run)

## Inputs

| Format | Notes |
|---|---|
| `.fastq` / `.fastq.gz` / `.fq` / `.fq.gz` | Raw HiFi reads |
| `.parquet` | Must have `read_id` VARCHAR, `sequence1` VARCHAR, `qual1` UTINYINT[] |

Reads are expected to carry the full Karst-protocol adapter+UMI structure on
both ends. Both strand orientations are accepted; the pipeline canonicalizes
to the + strand internally.

## Pipeline architecture

Stages are SQL scripts under `sql/`, executed in numeric order by `run.sh`:

| Stage | File | Purpose |
|---|---|---|
| 00 | `00_ingest.sql` | Quality and length filtering via `filter_read` (Q ≥ min_q, length range, ≤ N ambiguous bases) |
| 05 | `05_positive_filter.sql` | Optional 16S positive filter — `align_minimap2 map-hifi` against a reference; keep reads with ≥ 1 kb alignment. Skipped if no `positive_ref_path` is set |
| 10 | `10_umi_extract.sql` | Extract 18 bp UMI halves from both read termini using `extract_linked_amplicon` against 4 anchor sequences; build canonical UMI pairs (regex-validated); dereplicate; vsearch-cluster; remove chimeric UMI pairs by requiring all four 18-mer halves (forward + RC for both ends) to be owned by the same UMI cluster |
| 20 | `20_umi_bin.sql` | Map each read to a UMI bin via Hamming distance (`match_short_barcodes`, ≤ `umi_max_nm_per_half` per half); apply RoF (read-orientation fraction), UME (UMI matching error), and BCR (bin-to-cluster ratio) filters |
| 30 | `30_consensus.sql` | Per-bin abPOA MSA + `compute_msa_consensus` (Q-weighted); primer-trim with `extract_linked_amplicon`; gate by coverage ≥ `umi_coverage_min` |
| 40 | `40_variants.sql` | Homopolymer-mask consensuses, two-pass vsearch clustering at `subcluster_id` (symmetric centroid resolution) |
| 45 | `45_variant_calling.sql` | For each multi-member cluster: minimap2 align members to centroid → pileup → find SNP positions → per-read signature → phase reads by signature → abPOA + Q-weighted MSA to produce a phased consensus per variant |
| 50 | `50_export.sql` | Build `export_consensus`, `export_variants`, `export_unique` views with read_fastx-compatible schemas |
| 99 | `99_cleanup.sql` | (Optional) Drop intermediate tables and CHECKPOINT to compact disk |

`run.sh` skips stage 05 when no `positive_ref_path` is set in the params
file, and skips stage 45 when no UMI cluster has ≥ 2 members.

## Parameters

All parameters are set as DuckDB session variables via the params file. The
default file is `params/revio.sql`. To use a different parameter set, copy
the file and pass it as `--params`.

| Variable | Default | Meaning |
|---|---|---|
| `min_q` | 28 | Minimum average Phred quality for `filter_read` |
| `min_len` / `max_len` | 3500 / 6000 | Read length range (bp) |
| `fw1`, `fw2`, `rv1`, `rv2` | Karst protocol | Outer/inner forward and reverse primers; their reverse complements `*_rc` are also set |
| `umi_pair_pattern` | `^([ATCG]{3}[CT][AG]){3}…$` | Regex validation for the canonical 36 bp UMI pair |
| `umi_cluster_id` | 0.90 | vsearch identity for the UMI dereplication cluster step |
| `umi_max_nm_per_half` | 6 | Max Hamming distance per 18 bp UMI half when binning reads to UMI clusters |
| `ume_mean_max`, `ume_sd_max` | 3.0 / 30.0 | UMI matching error filter (mean and SD of summed NM per bin) |
| `ro_frac` | 0.05 | Minimum strand-balance fraction (least(n_plus, n_neg) / n) |
| `bin_cluster_ratio` | 10.0 | Maximum reads-per-bin to UMI-cluster-size ratio (rejects PCR-collapsed bins) |
| `min_bin_size` | 5 | Minimum reads per UMI bin |
| `umi_coverage_min` | 5 | Minimum UMI bin support (`ubs`) for a consensus to be retained |
| `subcluster_id` | 0.995 | vsearch identity for two-pass per-bin consensus sub-clustering |
| `variant_min_support` | 3 | Minimum UMI bins per phased variant |
| `snp_min_alt_depth` | 2 | Minimum alt-base depth at a position to be called a SNP site |
| `positive_ref_path` | (unset) | Path to a 16S reference FASTA. If set, enables the positive filter at stage 05 |

## Outputs

### `consensus.parquet` / `consensus.fa`

Per-UMI-bin consensus sequences (one row per UMI bin that passed all filters).

Parquet schema (read_fastx-compatible):

| Column | Type | Description |
|---|---|---|
| `sequence_index` | BIGINT | Stable index for the row |
| `read_id` | VARCHAR | UMI bin ID (e.g. `u_1234`) |
| `comment` | VARCHAR | `ubs=<n>` — number of reads contributing to this consensus |
| `sequence1` | VARCHAR | Primer-trimmed amplicon sequence (+ strand) |
| `sequence2` | VARCHAR | NULL (single-end) |
| `qual1` | UTINYINT[] | Per-base Phred quality from `compute_msa_consensus` |
| `qual2` | UTINYINT[] | NULL |

### `variants.parquet` / `variants.fa`

Phased intra-genomic variants. Rows exist only for UMI clusters with
`≥ variant_min_support` members and at least one SNP position.

| Column | Type | Description |
|---|---|---|
| `sequence_index` | BIGINT | Stable index |
| `read_id` | VARCHAR | `<cluster_id>_var<N>` — variant identifier |
| `comment` | VARCHAR | `support=<n>` — number of UMI bins supporting the variant |
| `sequence1`, `sequence2`, `qual1`, `qual2` | as above | |

### `unique.parquet`

Deduplicated final set: each row is one unique sequence either from a phased
variant or from a single-member cluster.

| Column | Type | Description |
|---|---|---|
| `read_id` | VARCHAR | Variant ID or single-bin cluster ID |
| `sequence1` | VARCHAR | Consensus sequence |
| `qual1` | UTINYINT[] | Per-base Phred |
| `umi_count` | BIGINT | Number of UMI bins (molecules) supporting this sequence |
| `total_reads` | BIGINT | Sum of `ubs` across supporting bins (raw read count) |

`umi_count` is the abundance metric of choice (PCR-corrected); `total_reads`
is informational for confidence.

## Database schema

After a normal run (no cleanup), `pipeline.duckdb` contains the per-stage
intermediates. After `--cleanup safe` it contains only the tables documented
below; after `--cleanup aggressive` the per-read tables are emptied.

### Final retained tables

| Table | Rows (typical) | Key columns | Purpose |
|---|---|---|---|
| `umi_ref` | 1 per UMI cluster centroid | `umi_id` PK, `size`, `umi_seq`, `u1`, `u2` | UMI cluster definitions (post-chimera-filter) |
| `bin_pass` | 1 per passing UMI bin | `umi_id` PK, `n`, `n_plus`, `n_neg`, `ume_mean`, `ume_sd`, `cluster_size`, `bcr` | Bins that passed RoF/UME/BCR filters; QC stats |
| `high_cov_consensus` | 1 per UMI bin | `bin_id` PK, `seq`, `qual`, `ubs` | Per-bin primer-trimmed consensus with quality |
| `cluster_members` | 1 per (bin, cluster) pair | `cluster_id`, `bin_id` | Bin → sub-cluster mapping from the two-pass vsearch |
| `variants` | 1 per phased variant | `variant_id` PK, `support`, `seq`, `qual` | Phased intra-genomic variant consensus |
| `variant_bins` | 1 per variant | `variant_id` PK, `cluster_id`, `signature`, `support`, `read_ids` | Variant manifest (UMI bin IDs contributing) |

### Per-read tables (retained by `safe` cleanup, dropped by `aggressive`)

| Table | Rows | Key columns | Purpose |
|---|---|---|---|
| `reads_unfiltered` | 1 per input read | `read_id` PK, `seq`, `qual` | Reads after Q/length filter, before 16S filter |
| `reads` | 1 per filtered read | same | Reads after 16S filter (or copy of `reads_unfiltered`) |
| `bin_reads` | 1 per read-to-bin assignment | `bin_id`, `read_id`, `strand`, `seq`, `qual` | Final read→bin mapping (one read may appear in 0 or 1 bin) |

### Export views (always present)

| View | Source | Use |
|---|---|---|
| `export_consensus`, `export_consensus_fasta` | `high_cov_consensus` | Per-bin consensus output |
| `export_variants`, `export_variants_fasta` | `variants` | Phased variant output |
| `export_unique` | `variants` + `variant_bins` + `high_cov_consensus` | Deduplicated unique sequences with abundance |

### Querying the DB after a run

```sql
-- Top 10 most abundant sequences
SELECT read_id, umi_count, total_reads, length(sequence1) AS bp
FROM export_unique
ORDER BY umi_count DESC LIMIT 10;

-- All UMI bins contributing to a specific variant
SELECT vb.variant_id, vb.support, unnest(vb.read_ids) AS bin_id
FROM variant_bins vb
WHERE vb.variant_id = 'u_20_var1';

-- Coverage histogram of UMI bins
SELECT ubs, count(*) AS n_bins
FROM high_cov_consensus
GROUP BY ubs ORDER BY ubs;

-- UMI binning QC: which bins were rejected?
SELECT * FROM bin_pass WHERE bcr > 5;
```

## Cleanup

Three modes, controlled by `--cleanup`:

| Mode | Drops | Use case |
|---|---|---|
| `none` (default) | nothing | Default; keep full intermediate state for debugging |
| `safe` | All purely intermediate tables (Hamming hit tables, pileup, MSA outputs, UMI extraction scratch, clustering scratch) | Long-term storage; preserves per-read provenance |
| `aggressive` | `safe` + per-read tables emptied (`reads`, `reads_unfiltered`, `bin_reads`) | Smallest archive; you cannot re-derive bin-level statistics without re-running stages 00–20 |

The cleanup stage runs `CHECKPOINT` after dropping. For maximum compaction
(reclaim slack that DuckDB still holds), follow up with an
`EXPORT DATABASE`/`IMPORT DATABASE` to a new file.

## Testing

```bash
# Synthetic-fixture test suite (no network, ~5 s)
DUCKDB=/path/to/duckdb-miint/build/release/duckdb ./test/run_tests.sh
```

The fixture under `test/fixtures/synthetic.fq` exercises stages 00 → 50
end-to-end with a small hand-crafted UMI structure. CI runs this same script
on every push.

## Approximate runtime

Measured on Linux, 12 cores, 64 GB RAM, NVMe SSD, ~1.57 M Revio HiFi reads
(~15 GB uncompressed FASTQ, 40 K UMI bins after filtering):

| Stage | Time |
|---|---|
| `00_ingest` | ~2.5 min |
| `05_positive_filter` (minimap2 map-hifi vs 88_otus) | ~15.5 min |
| `10_umi_extract` | ~40 s |
| `20_umi_bin` | 6–15 min (scales with `umi_max_nm_per_half`²) |
| `30_consensus` (abPOA + Q-weighted MSA) | ~6.5 min |
| `40_variants` (vsearch two-pass clustering) | ~60–65 min |
| `45_variant_calling` | ~5 min |
| `50_export` | < 1 s |
| **Total** | **~100–106 min** |

Stage 40 dominates wall-clock; stage 05 is the optional 16S filter and can be
skipped if the input is already a known operon amplicon.

## Algorithmic notes

A few non-obvious choices, with rationale, so the workflow is reproducible
end-to-end:

- **abPOA for MSA, not MAFFT or Racon.** abPOA is a SIMD-vectorized partial
  order alignment library and produces per-position graph alignments suitable
  for direct consensus voting in a single pass. We use it both for per-bin
  consensus (stage 30) and for per-variant phased consensus (stage 45). The
  alternative inside duckdb-miint, `align_mafft`, is preserved but bypassed
  here because abPOA at the bin sizes we encounter (5–50 reads per bin) is
  both faster and producing comparable or better alignments.

- **Q-weighted MSA consensus, not majority vote.** `compute_msa_consensus`
  uses the per-base Phred from the aligned reads to weight each column,
  emitting both a consensus base and a posterior Phred. Revio HiFi quality
  scores are well-calibrated, so weighting by Q materially reduces error in
  low-coverage columns versus an unweighted majority vote.

- **Canonical (+ strand) UMI normalization.** A read sequenced on the
  − strand has its UMI halves swapped and reverse-complemented relative to a
  + strand read. Stage 10 stores both halves in their + strand form so that
  downstream Hamming binning has a single canonical representation per
  molecule. Without this normalization a single UMI bin would split into two
  by strand.

- **Chimeric UMI filter via four-way half ownership.** A chimeric UMI pair
  (two halves from different molecules) is detected by requiring that all
  four 18-mers (u1, u2, RC(u1), RC(u2)) of a UMI cluster centroid be owned
  (in the arg-max-by-size sense) by that same UMI cluster. PCR/sequencing
  chimeras typically fail at least one of these checks.

- **Hamming on packed `uint64_t`, not edit distance.** UMI halves are
  length-uniform 18 bp, so binning uses `match_short_barcodes` (Hamming with
  bit-parallel popcount) instead of an edit-distance aligner. `max_nm` is a
  hard cap; the default of 6 mismatches per half (12 total over the 36 bp
  pair) tolerates per-base sequencing error in Revio HiFi while keeping the
  bin assignment unique.

- **Homopolymer mask before sub-clustering.** Stage 40 collapses runs of
  ≥ 3 identical bases to length 2 (`AAAA…` → `AA`) before vsearch
  clustering. Homopolymer-length errors are the dominant residual error mode
  in HiFi consensus; without this mask the sub-clusterer would over-split
  the same biological sequence based on indel artifacts. The mask is local
  to stage 40 — the final consensus and variant sequences are emitted from
  the unmasked `high_cov_consensus` and `variant_msa` tables.

- **Two-pass symmetric sub-clustering.** vsearch's greedy clustering is
  order-dependent. Pass 1 clusters bins (largest-first); pass 2 clusters the
  resulting centroids again at the same threshold, which collapses cases
  where two near-identical centroids arose because their nearest neighbours
  attached to different seeds. The final `cluster_members` reflects the
  consensus of both passes.

- **Phasing via per-read SNP signatures.** Stage 45 calls pileup-derived
  SNP sites within each multi-member cluster, then groups bins by the exact
  string of their bases at those positions. Each unique signature with
  `≥ variant_min_support` bins becomes a phased variant. This is a single
  pass over the pileup, no iterative EM — adequate because each bin is itself
  a per-molecule consensus.

## License

GPL — see [LICENSE](LICENSE).
