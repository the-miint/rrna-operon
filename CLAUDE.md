# CLAUDE.md — guidance for AI assistants working on this repo

This file is the orientation for an AI agent (Claude) editing or extending
`rrna-operon`. Read alongside `README.md` (user-facing) — this file covers
the conventions and gotchas that are hard to recover from reading the code
cold.

## What this project is

A bash + DuckDB SQL pipeline that converts PacBio Revio HiFi reads carrying
Karst-protocol UMI structure into per-molecule rRNA operon consensus
sequences plus phased intra-genomic variants. Every primitive (read I/O,
primer trim, MSA, pileup, clustering, alignment) is a function exposed by
the `duckdb-miint` DuckDB extension. There is no Python and no
auto-installed pip dependency tree.

Authoritative dependency: the `miint` extension at
`https://ftp.microbio.me/pub/miint`, built against DuckDB **v1.5.2**.

## Two extension delivery paths — both must keep working

The codebase is invoked in two environments:

1. **Custom static build** (local dev): `../duckdb-miint/build/release/duckdb`
   has miint linked in. `LOAD miint;` is a no-op; `-unsigned` is accepted but
   unused.
2. **Stock DuckDB + community-install** (CI, users): stock v1.5.2 CLI plus
   `INSTALL miint FROM 'https://ftp.microbio.me/pub/miint';`. Every session
   needs `-unsigned -cmd "LOAD miint;"` because the extension is unsigned
   and not auto-loaded.

`run.sh` and `test/run_tests.sh` route every duckdb invocation through a
shell function:

```bash
duckdb() { "$DUCKDB" -unsigned -cmd "LOAD miint;" "$@"; }
```

**If you add a new duckdb invocation anywhere, call `duckdb`, not `$DUCKDB`
directly.** Otherwise it works on the static build and silently fails on
the community-install path.

## Session-variable semantics

DuckDB `SET VARIABLE …` is connection-scoped. The pipeline parameters in
`params/*.sql` are session variables, so **every stage's SQL block must
re-read the params file**. `run_phase` in `run.sh` does this:

```bash
duckdb "$DB" <<EOF
${INPUT_SQL}
SET VARIABLE output_dir = '${OUTPUT_DIR}';
.read ${PARAMS}
.read ${sql_file}
EOF
```

If you add a standalone duckdb invocation that calls `getvariable(…)`, you
must `.read` the params file first or pre-set the variable. The export
block in `run.sh` deliberately doesn't read params because it only references
already-materialized tables/views.

## Stage flow and skip logic

Stages live in `sql/[0-9]*.sql`, executed in lexical order by `run.sh`:

| Stage | File | When skipped |
|---|---|---|
| 00 | `00_ingest.sql` | never |
| 10 | `10_umi_extract.sql` | never |
| 20 | `20_umi_bin.sql` | never |
| 30 | `30_consensus.sql` | never |
| 40 | `40_variants.sql` | never |
| 45 | `45_variant_calling.sql` | when no UMI cluster has ≥ 2 members |
| 50 | `50_export.sql` | never |

The 45 skip is a runtime SELECT against `cluster_members` in `run.sh`.

A previous `05_positive_filter.sql` ran minimap2 against an 88_otus 16S
reference to drop off-target reads. It was removed: bacterial UMI
primers do not hybridise to human/off-target DNA, so the
primer-extraction step at stage 10 implicitly host-filters at far
lower cost. Don't reintroduce a positive filter without evidence
that the implicit one is insufficient.

There is no `99_cleanup.sql` — that file existed historically but was
deleted; the cleanup model is now "export retained tables to parquet, then
delete the working DB" (see below).

## Storage model: parquet outputs, ephemeral DB

`pipeline.duckdb` is treated as a working file. At the end of `run.sh`:

1. Sequence outputs go to `consensus.parquet`, `variants.parquet`,
   `unique.parquet` (+ FASTA mirrors).
2. Provenance/QC tables are COPY'd to:
   - `umi_ref.parquet` — UMI cluster definitions
   - `bin_pass.parquet` — passing UMI bins with QC stats
   - `cluster_members.parquet` — bin → sub-cluster mapping
   - `variant_bins.parquet` — variant → contributing bin IDs
   - `primer_extract_status.parquet` — per-read primer-extraction outcome
3. `pipeline.duckdb`, `pipeline.duckdb.wal`, and `pipeline.duckdb.tmp/`
   are deleted unless `--keep-db` is passed.

**If you add a new "provenance" table — i.e., a per-read or per-bin
manifest with a boolean — follow the same pattern**: build the table in
the appropriate stage SQL file, then COPY it in `run.sh`. The
`primer_extract_status` table is the template; mirror it.

DuckDB `CHECKPOINT` does not shrink the storage file — it only flushes
WAL. To compact, you need EXPORT DATABASE/IMPORT DATABASE. That's why
the cleanup model deletes the DB rather than trying to shrink it.

## CI conventions

`.github/workflows/test.yml`:

- DuckDB version pinned to `v1.5.2` via `env.DUCKDB_VERSION`. This must
  match the version of miint published at
  `https://ftp.microbio.me/pub/miint/v<X.Y.Z>/`. **Bump both together**
  when a new miint release ships.
- Cache key includes `DUCKDB_VERSION` so a bump invalidates the cache.
- The smoke-test step exercises `COPY <retained_table> TO …parquet`
  commands against the post-test fixture DB, not a deleted-stage SQL
  file.

## Karst protocol assumptions (encoded in `params/*.sql`)

- 36 bp UMI = two 18 bp halves, one at each read terminus, with a
  canonical pair regex `^([ATCG]{3}[CT][AG]){3}…$` validating the
  third-base degeneracy.
- 4 anchor primers (`fw1`, `fw2`, `rv1`, `rv2`) flank the UMIs. Reverse
  complements (`*_rc`) are pre-computed in the params file rather than
  derived at runtime.
- UMIs from − strand reads are stored in canonical + strand form by
  swap+RC normalisation in stage 20 — without this, every UMI bin would
  split by strand.

## Known gotchas

- **`extract_linked_amplicon` returns NULL** when an anchor can't be
  located within `ceil(len(anchor) * error_rate)` errors. The pipeline
  silently drops those rows at stage 10. Recover them via
  `primer_extract_status.parquet`.
- **Supplementary minimap2 alignments emit hard-clipped CIGARs**. Stage
  45 must filter with `NOT alignment_is_unmapped(flags) AND
  alignment_is_primary(flags)` — passing secondary/supplementary
  alignments to `compute_pileup` will fail with a CIGAR/sequence-length
  mismatch.
- **vsearch clustering is order-dependent.** Two pipeline runs on
  identical input can produce off-by-one bin/variant counts (e.g.
  39,874 vs 39,876). Don't write tests that assert exact equality on
  large-scale runs.
- **dependencies.sh in upstream `longread_umi` is Latin-1-encoded.** If
  you edit it via a tool that auto-converts to UTF-8, you'll mangle the
  `ø` in "Søren Karst". Use `sed -i` for surgical edits or preserve the
  encoding explicitly.
- **`gh repo create` against the-miint requires explicit visibility
  approval** when run via Claude Code. Ask before defaulting to
  `--public`.

## Testing

```bash
DUCKDB=/path/to/duckdb-miint-build/duckdb ./test/run_tests.sh
```

The fixture (`test/fixtures/synthetic.fq`, 11 hand-crafted reads with 2
distinct UMI bins) exercises stages 00 → 50 end-to-end. Stage 05 is
always skipped under the test params (no 16S reference). The fixture is
too small to exercise stage 45 (no multi-member clusters), so the
`variant_calling` code path is not test-covered — beware when changing
it.

CI runs this same script on every push and PR.

## Related repos and external context

- **`duckdb-miint`** (sibling): the DuckDB extension hosting every
  bioinformatic primitive this pipeline uses (`align_minimap2`,
  `align_abpoa`, `extract_linked_amplicon`, `compute_msa_consensus`,
  `compute_pileup`, `cluster_sequences_vsearch`, etc.). Function docs
  live at `../duckdb-miint/docs/{scalar,table}-functions.md`.
- **`longread_umi`** (Karst et al., upstream reference pipeline): used
  as a benchmarking/validation comparator. PR #59 to
  `SorenKarst/longread_umi` contains portability + vsearch-compat fixes
  derived from this work. The local clone at `~/longread_umi/` has
  additional un-PR'd modifications.
- **Plan document**: `~/.claude/plans/indexed-napping-peach.md` is the
  living design/history doc for this work. Update it when you make a
  load-bearing decision (default change, storage-model change, etc.)
  rather than burying the rationale in commit messages alone.

## When in doubt

- Prefer extending existing tables/files over adding new ones.
- A new provenance file is cheap; a new dependency is expensive.
- The pipeline's invariant is "parquet outputs are self-sufficient" —
  don't add a step that requires `pipeline.duckdb` to be present
  downstream.
- If a stage gets slower or produces different counts on the standard
  E. coli sample, treat that as a regression and find the cause before
  committing.
