# genomes

Genome-level annotation used by [`sesame-cli`](https://github.com/zwdzwd/sesame-cli)
(CNV binning and ideogram plotting) and other zhou-lab plotting tools. Kept in
its own repo, separate from the platform/probe annotation in
[`InfiniumAnnotation`](https://github.com/zhou-lab/InfiniumAnnotation), so it can
be reused independently of any array platform.

## Layout

One folder per genome build, versioned by **git tag** (`v1`, `v2`, …):

```
<build>/SHA256SUMS        coreutils digests of the files below (the trust anchor)
<build>/seqinfo.tsv.gz    chrom <TAB> length                    (genome tiling)
<build>/gaps.tsv.gz       chrom <TAB> start <TAB> end           (assembly gaps to skip)
<build>/cytoband.tsv.gz   chrom start end band stain            (ideogram, plot-only)
```

`seqinfo` + `gaps` drive `getBinCoordinates` (tile the genome into bins, carve
out the gaps); `cytoband` is only for the ideogram. All are exported from
`sesameData` by `tools/export_genomeinfo.R` in the sesame-cli repo.

## Fetch

```sh
sesame fetch genome hg38     # -> <store>/genome/hg38/{seqinfo,gaps,cytoband}.tsv.gz
```

`sesame` pulls `<build>/SHA256SUMS` first, verifies it against a digest compiled
into the build, then verifies every file against it — a hard chain. Files are
raw-served over `https://github.com/zhou-lab/genomes/raw/<tag>/<build>/<file>`.

## Builds

| build | source |
|---|---|
| hg38 | `sesameData_getGenomeInfo("hg38")` |
| mm10 | `sesameData_getGenomeInfo("mm10")` |
| mm39 | `sesameData_getGenomeInfo("mm39")` |
