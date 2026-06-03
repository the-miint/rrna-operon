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
- [Install](#install)
- [Quick start](#quick-start)
- [Inputs](#inputs)
- [Pipeline architecture](#pipeline-architecture)
- [Parameters](#parameters)
- [Outputs](#outputs)
- [Operon annotation (optional)](#operon-annotation-optional)
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
- ~30 GB free disk for a 1.5 M-read input (intermediate `pipeline.duckdb`;
  deleted at the end of the run by default — see [Cleanup](#cleanup))
- *(optional)* [barrnap](https://github.com/tseemann/barrnap) ≥ 1.10 for
  [operon annotation](#operon-annotation-optional) (rRNA + ITS tRNA); installed
  by `install.sh` into a dedicated `barrnap` conda env. Not needed to run the
  core pipeline.

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

## Install

The fastest way to a working toolchain is the bundled installer:

```bash
./install.sh
```

It is idempotent (re-running skips anything already present) and sets up three
things:

1. **A conda/mamba package manager.** `mamba` is preferred when present (its
   solver is much faster). If neither `mamba` nor `conda` is found on `PATH`
   or in a standard location,
   [Miniforge](https://github.com/conda-forge/miniforge) is installed to
   `~/miniforge3` (it bundles `mamba` and conda-forge).
2. **barrnap** (with aragorn/infernal/diamond) in a dedicated `barrnap` conda
   env — used for the optional [operon annotation](#operon-annotation-optional).
3. **The DuckDB v1.5.3 CLI + the `miint` extension** — the CLI is placed in
   `./bin/duckdb` and `miint` is registered via
   `INSTALL miint FROM 'https://ftp.microbio.me/pub/miint'`.

Useful flags:

- `./install.sh --no-barrnap` — DuckDB + miint only (skip conda/barrnap); this
  is what CI uses.
- `DUCKDB_DEST=~/.local/bin ./install.sh` — change where the DuckDB CLI lands.

`run.sh` resolves its DuckDB binary from `$DUCKDB`, the `--duckdb` flag, or the
custom static-build path — so point it at the CLI the installer placed:

```bash
DUCKDB="./bin/duckdb" ./run.sh --output out/ reads.fastq.gz
```

## Quick start

```bash
# Run the pipeline on a FASTQ
./run.sh \
    --params params/revio.sql \
    --output output/ \
    --duckdb /path/to/duckdb-miint/build/release/duckdb \
    reads.fastq.gz

# Pass --keep-db to retain pipeline.duckdb for interactive querying.
```

Outputs land in `output/`:

| File | Content |
|---|---|
| `consensus.parquet` / `consensus.fa` | Per-UMI consensus sequences |
| `variants.parquet` / `variants.fa` | Phased variant consensuses |
| `unique.parquet` | Deduplicated unique sequences with abundance |
| `umi_ref.parquet` | UMI cluster definitions |
| `bin_pass.parquet` | Per-bin QC statistics |
| `cluster_members.parquet` | Bin → sub-cluster mapping |
| `variant_bins.parquet` | Variant → contributing-bin manifest |
| `primer_extract_status.parquet` | Per-read pass/fail for UMI primer extraction |

The intermediate `pipeline.duckdb` is deleted by default at the end of the
run. Use `--keep-db` to retain it.

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
| 10 | `10_umi_extract.sql` | Extract 18 bp UMI halves from both read termini using `extract_linked_amplicon` against 4 anchor sequences; build canonical UMI pairs (regex-validated); dereplicate; vsearch-cluster; remove chimeric UMI pairs by requiring all four 18-mer halves (forward + RC for both ends) to be owned by the same UMI cluster. Reads with no recoverable primer pair drop here — see `primer_extract_status.parquet` |
| 20 | `20_umi_bin.sql` | Map each read to a UMI bin via Hamming distance (`match_short_barcodes`, ≤ `umi_max_nm_per_half` per half); apply RoF (read-orientation fraction), UME (UMI matching error), and BCR (bin-to-cluster ratio) filters |
| 30 | `30_consensus.sql` | Per-bin abPOA MSA + `compute_msa_consensus` (Q-weighted); primer-trim with `extract_linked_amplicon`; gate by coverage ≥ `umi_coverage_min` |
| 40 | `40_variants.sql` | Homopolymer-mask consensuses, two-pass vsearch clustering at `subcluster_id` (symmetric centroid resolution) |
| 45 | `45_variant_calling.sql` | For each multi-member cluster: minimap2 align members to centroid → pileup → find SNP positions → per-read signature → phase reads by signature → abPOA + Q-weighted MSA to produce a phased consensus per variant |
| 50 | `50_export.sql` | Build `export_consensus`, `export_variants`, `export_unique` views with read_fastx-compatible schemas |

`run.sh` skips stage 45 when no UMI cluster has ≥ 2 members.

Host filtering is implicit: bacterial UMI primers do not hybridise to
human or other off-target DNA, so reads that lack the flanking primer
sequence are dropped at stage 10 by `extract_linked_amplicon` (which
returns NULL when an anchor cannot be located). An earlier
minimap2-based 16S-mapping filter was removed in favour of relying on
this implicit host-filter, which is both faster and more selective at
the boundary between bacterial and off-target DNA.

After the stages, `run.sh` writes a set of Parquet files (sequence outputs
plus provenance/QC tables — see [Outputs](#outputs)) and, by default,
deletes the working `pipeline.duckdb` file. Pass `--keep-db` to retain it
for interactive querying.

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

### Provenance / QC tables (Parquet, zstd-compressed)

Five additional Parquet files are written alongside the sequence
outputs. Together they let you trace any unique sequence back to the
contributing UMI bins, replay variant calling at different thresholds,
and audit each stage's filters — all without keeping the multi-gigabyte
intermediate database.

#### `umi_ref.parquet` — UMI cluster definitions

| Column | Type | Description |
|---|---|---|
| `umi_id` | VARCHAR | Cluster identifier (`u_<rank>` by size) |
| `size` | BIGINT | Number of dereplicated UMI sequences in the cluster |
| `umi_seq` | VARCHAR | Canonical 36 bp UMI pair (u1 \|\| u2) |
| `u1` | VARCHAR | 18 bp 5′ half |
| `u2` | VARCHAR | 18 bp 3′ half |

#### `bin_pass.parquet` — UMI bin QC

One row per UMI bin that survived the RoF / UME / BCR filters at stage 20.

| Column | Type | Description |
|---|---|---|
| `umi_id` | VARCHAR | FK → `umi_ref.umi_id` and `cluster_members.bin_id` |
| `n` | BIGINT | Total reads assigned to the bin |
| `n_plus`, `n_neg` | BIGINT | Read counts per strand |
| `ume_mean`, `ume_sd` | DOUBLE | UMI matching error (mean/SD of summed NM per bin) |
| `cluster_size` | BIGINT | Original UMI-cluster centroid count from stage 10 |
| `bcr` | DOUBLE | Bin-to-cluster ratio (`n / cluster_size`) |

#### `cluster_members.parquet` — bin → sub-cluster mapping

| Column | Type | Description |
|---|---|---|
| `cluster_id` | VARCHAR | Centroid bin ID from the two-pass sub-clustering |
| `bin_id` | VARCHAR | Member bin ID (== `cluster_id` for the centroid itself) |

#### `primer_extract_status.parquet` — per-read UMI primer extraction outcome

One row per read that entered stage 10 (i.e., passed stage 00).

| Column | Type | Description |
|---|---|---|
| `read_id` | VARCHAR | Input read identifier |
| `umis_extracted` | BOOLEAN | True if `extract_linked_amplicon` located both flanking primers on at least one strand and the resulting UMI pair matched the canonical regex |

`umis_extracted = false` means the read survived quality/length and 16S
filters but had primers too degraded (errors above the `error_rate`
budget), missing, or unrecognizable — common causes include adapter
contamination or read truncation. These reads drop out of UMI binning
silently; this file is the only post-run record of them.

#### `variant_bins.parquet` — variant → contributing bins

One row per phased variant.

| Column | Type | Description |
|---|---|---|
| `variant_id` | VARCHAR | `<cluster_id>_var<N>` |
| `cluster_id` | VARCHAR | FK → `cluster_members.cluster_id` |
| `signature` | VARCHAR | Per-read SNP signature string |
| `support` | BIGINT | Number of UMI bins with this signature |
| `read_ids` | VARCHAR[] | Bin IDs contributing to the variant |

### Querying the outputs

The Parquet files are self-sufficient — no DuckDB database needed.

```sql
-- Top 10 most abundant sequences
SELECT read_id, umi_count, total_reads, length(sequence1) AS bp
FROM 'output/unique.parquet'
ORDER BY umi_count DESC LIMIT 10;

-- All UMI bins contributing to a specific variant
SELECT variant_id, support, unnest(read_ids) AS bin_id
FROM 'output/variant_bins.parquet'
WHERE variant_id = 'u_20_var1';

-- Coverage histogram of UMI bins
SELECT cast(comment[5:] AS BIGINT) AS ubs, count(*) AS n_bins
FROM 'output/consensus.parquet'
GROUP BY ubs ORDER BY ubs;

-- UMI binning QC: bins on the edge of the BCR filter
SELECT * FROM 'output/bin_pass.parquet' WHERE bcr > 5;

-- How many reads lacked recoverable UMI primers (implicit host filter)?
SELECT count(*)                                   AS n_post_ingest,
       count(*) FILTER (WHERE NOT umis_extracted) AS n_no_primers,
       (count(*) FILTER (WHERE NOT umis_extracted))::DOUBLE / count(*) AS frac_no_primers
FROM 'output/primer_extract_status.parquet';

-- Cross-reference: which UMI clusters became multi-variant?
SELECT cluster_id, count(*) AS n_variants, sum(support) AS total_bins
FROM 'output/variant_bins.parquet'
GROUP BY cluster_id HAVING count(*) > 1
ORDER BY n_variants DESC;
```

## Operon annotation (optional)

The per-molecule consensus operons can be annotated for rRNA genes (16S, 23S)
and the intervening ITS tRNAs (tRNA-Ile / tRNA-Ala / tRNA-Glu) with
[barrnap](https://github.com/tseemann/barrnap) — installed by `install.sh`
into the `barrnap` conda env. This is a post-processing step, not part of the
core pipeline.

```bash
conda run -n barrnap barrnap --kingdom bac --trna --fast \
    out/consensus.fa > out/operons.gff
```

Pull the annotation back into Parquet with miint's GFF reader (`read_gff`
parses the attributes column into a DuckDB `MAP`):

```sql
COPY (SELECT seqid AS bin_id, type, position, stop_position, strand,
             attributes['Name'] AS feature
      FROM read_gff('out/operons.gff'))
  TO 'out/operon_annotations.parquet' (FORMAT PARQUET);
```

Notes:

- We run barrnap in **`--fast`** mode for amplicon consensus; see the
  [barrnap docs](https://github.com/tseemann/barrnap) for the speed/accuracy
  trade-off.
- **`--kingdom`** defaults to `bac`; for environmental samples with archaea,
  also run `--kingdom arc` and merge the GFFs.
- **tRNA** is nearly free and resolves the ITS operon type — annotate
  `variants.fa` instead of `consensus.fa` to study per-copy ITS differences.

## Cleanup

By default `run.sh` writes the Parquet outputs (see above), deletes the
working `pipeline.duckdb` and its `.wal` / `.tmp` companions, and exits.
Typical end state on a 1.5 M-read input: ~2 MB of Parquet plus the optional
FASTA files (~170 MB total, dominated by `consensus.fa`).

To keep the intermediate database for interactive querying or debugging,
pass `--keep-db`. The retained DB is large — ~30 GB on a 1.5 M-read run —
and is not pruned.

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
| `10_umi_extract` | ~40 s |
| `20_umi_bin` | 6–15 min (scales with `umi_max_nm_per_half`²) |
| `30_consensus` (abPOA + Q-weighted MSA) | ~6.5 min |
| `40_variants` (vsearch two-pass clustering) | ~60–65 min |
| `45_variant_calling` | ~5 min |
| `50_export` | < 1 s |
| **Total** | **~100–106 min** |

Stage 40 dominates wall-clock.

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
