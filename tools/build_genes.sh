#!/usr/bin/env bash
# Build a tabix-indexed GENCODE transcript-model BED for one genome build.
#
#   build_genes.sh <build> <gencode_gtf_url> [outroot]
#
# Produces  <outroot>/<build>/genes.bed.gz(+.tbi)  and records the source URL in
# <outroot>/<build>/genes.source. Reproduces the gene models behind
# sesameData's genomeInfo$txns (see tools/build_genes.py). bgzip/tabix come from
# $BGZIP/$TABIX (default: on PATH).
set -euo pipefail

build=${1:?usage: build_genes.sh <build> <gencode_gtf_url> [outroot]}
url=${2:?need a GENCODE .gtf.gz URL}
outroot=${3:-..}
here=$(cd "$(dirname "$0")" && pwd)
BGZIP=${BGZIP:-bgzip}
TABIX=${TABIX:-tabix}

# GENCODE release label (v36 / vM25 / ...), for the log line only. The output
# name is unversioned (genes.bed.gz); the build is pinned by git tag and the exact
# source URL is recorded in genes.bed.gz.source.
ver=$(printf '%s\n' "$url" | grep -oE 'gencode\.v[M0-9]+' | head -1 | sed 's/gencode\.//')
name="genes.bed.gz"

dst="$outroot/$build"
mkdir -p "$dst"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

echo "[$build/$ver] fetch $url"
curl -fsSL "$url" -o "$tmp/gencode.gtf.gz"

echo "[$build/$ver] GTF -> BED12+"
python3 "$here/build_genes.py" "$tmp/gencode.gtf.gz" > "$tmp/genes.bed"

echo "[$build/$ver] sort + bgzip + tabix"
LC_ALL=C sort -k1,1 -k2,2n "$tmp/genes.bed" | "$BGZIP" -c > "$dst/$name"
"$TABIX" -p bed -f "$dst/$name"
echo "$url" > "$dst/$name.source"

echo "[$build/$ver] $(gzip -cd "$dst/$name" | wc -l | tr -d ' ') transcripts -> $dst/$name"
