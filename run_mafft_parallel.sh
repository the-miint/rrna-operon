#!/usr/bin/env bash
set -euo pipefail

# Parallel MAFFT alignment of per-bin reads.
# Input: TSV with columns (bin_id, read_id, sequence), sorted by bin_id.
# Output: gzipped TSV with columns (bin_id, read_id, aligned_sequence).
# Bins with <2 sequences are passed through unaligned (no MSA needed).

INPUT="${1:?Usage: $0 INPUT.tsv OUTPUT.tsv.gz [JOBS]}"
OUTPUT="${2:?Usage: $0 INPUT.tsv OUTPUT.tsv.gz [JOBS]}"
JOBS="${3:-$(nproc)}"

WORKDIR="$(mktemp -d "${PWD}/mafft_work_XXXXXX")"
trap 'rm -rf "$WORKDIR"' EXIT

BINS_DIR="${WORKDIR}/bins"
ALN_DIR="${WORKDIR}/aln"
mkdir -p "$BINS_DIR" "$ALN_DIR"

echo "Splitting $(wc -l < "$INPUT") reads into per-bin FASTAs..."
awk -F'\t' '{
    f = bins_dir "/" $1 ".fa"
    print ">" $2 >> f
    print $3 >> f
}' bins_dir="$BINS_DIR" "$INPUT"

N_BINS=$(ls "$BINS_DIR" | wc -l)
echo "Running MAFFT on ${N_BINS} bins with ${JOBS} parallel jobs..."

mafft_one() {
    local fa="$1"
    local bin_id
    bin_id="$(basename "$fa" .fa)"
    local out="${ALN_DIR}/${bin_id}.afa"
    local n_seqs
    n_seqs=$(grep -c "^>" "$fa")
    if [[ "$n_seqs" -lt 2 ]]; then
        cp "$fa" "$out"
    else
        mafft --quiet --preservecase --parttree "$fa" > "$out" 2>/dev/null
    fi
}
export -f mafft_one
export ALN_DIR

find "$BINS_DIR" -name '*.fa' | parallel -j "$JOBS" mafft_one {}

echo "Collecting aligned sequences..."
find "$ALN_DIR" -name '*.afa' -size +0c | sort | while read -r afa; do
    bin_id="$(basename "$afa" .afa)"
    awk -v bid="$bin_id" '
        /^>/ { if (seq) print bid "\t" rid "\t" seq; rid=substr($0,2); seq=""; next }
        { seq = seq $0 }
        END { if (seq) print bid "\t" rid "\t" seq }
    ' "$afa"
done | gzip > "$OUTPUT"

N_ALN=$(zcat "$OUTPUT" | wc -l)
echo "Done. ${N_ALN} aligned reads in ${OUTPUT}"
