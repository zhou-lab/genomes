# genomes

Genome-level annotation used by [`sesame-cli`](https://github.com/zwdzwd/sesame-cli)
(CNV binning, ideogram plotting, and region/gene tracks) and other zhou-lab
plotting tools. Kept in its own repo, separate from the platform/probe annotation
in [`InfiniumAnnotation`](https://github.com/zhou-lab/InfiniumAnnotation), so it
can be reused independently of any array platform.

## Layout

One folder per genome build, versioned by **git tag** (`v1`, `v2`, …):

```
<build>/SHA256SUMS               coreutils digests of the files below (the trust anchor)
<build>/seqinfo.tsv.gz           chrom <TAB> length              (genome tiling)
<build>/gaps.tsv.gz              chrom <TAB> start <TAB> end     (assembly gaps to skip)
<build>/cytoband.tsv.gz          chrom start end band stain      (ideogram, plot-only)
<build>/genes.bed.gz (+.tbi)     GENCODE transcript models (tabix BED12+)
```

`seqinfo` + `gaps` drive `getBinCoordinates` (tile the genome into bins, carve
out the gaps); `cytoband` is only for the ideogram; these three are exported from
`sesameData` by `tools/export_genomeinfo.R` in the sesame-cli repo. `genes.*` is
built here from GENCODE by `tools/build_genes.sh` (see **Gene models** below).

## Gene models

`genes.bed.gz` is a tabix-indexed transcript-model table — the plot-ready gene
track for `sesame region` and cinderplot's locus browser. One row per transcript,
sorted, `bgzip`ped, and indexed (`tabix -p bed`), so a locus is a single
random-access query:

```sh
tabix hg38/genes.bed.gz chr20:44616522-44655233
```

It is **BED12 + five named columns** (tab-separated; no trailing comma in the
block lists):

```
1 chrom          6 strand         11 blockSizes      16 transcript_name
2 chromStart(0)  7 thickStart     12 blockStarts     17 transcript_type
3 chromEnd       8 thickEnd       13 gene_id
4 transcript_id  9 itemRgb (0)    14 gene_name
5 score (0)     10 blockCount     15 gene_type
```

`thickStart`/`thickEnd` are the CDS span (`[min CDS start, max CDS end]`); a
non-coding transcript has `thickStart == thickEnd == chromEnd` (no thick drawn).
This reproduces the semantics of `sesameData:::build_GENCODE_gtf` (which built
`genomeInfo$txns`) — validated exon-for-exon against the ADA locus — but as a
plain indexable BED, not an R `GRangesList`.

Each build uses the GENCODE release its `sesameData` `txns` came from. The file
name is unversioned (like `seqinfo`/`gaps`/`cytoband`) — the release is pinned by
git tag, and the exact source URL is recorded next to the file in
`genes.bed.gz.source`. Bumping GENCODE means rerunning the builder and cutting a
new tag.

| build | GENCODE release | source GTF |
|---|---|---|
| hg38 | v36  | `gencode.v36.annotation.gtf.gz`  |
| mm10 | vM25 | `gencode.vM25.annotation.gtf.gz` |
| mm39 | vM31 | `gencode.vM31.annotation.gtf.gz` |

```sh
# tools/build_genes.sh <build> <gencode_gtf_url> [outroot]   (needs bgzip/tabix)
cd tools && ./build_genes.sh hg38 \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz ..
```

## Fetch

```sh
sesame fetch genome hg38     # -> <store>/genome/hg38/{seqinfo,gaps,cytoband,genes.*}
```

`sesame` pulls `<build>/SHA256SUMS` first, verifies it against a digest compiled
into the build, then verifies every file it lists against it — a hard chain, and
the file list is exactly `SHA256SUMS`, so the gene models come along automatically.
Files are raw-served over `https://github.com/zhou-lab/genomes/raw/<tag>/<build>/<file>`.

## Builds

| build | seqinfo / gaps / cytoband | gene models |
|---|---|---|
| hg38 | `sesameData_getGenomeInfo("hg38")` | GENCODE v36  |
| mm10 | `sesameData_getGenomeInfo("mm10")` | GENCODE vM25 |
| mm39 | `sesameData_getGenomeInfo("mm39")` | GENCODE vM31 |
