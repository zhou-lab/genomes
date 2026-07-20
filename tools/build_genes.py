#!/usr/bin/env python3
"""Convert a GENCODE GTF into a tabix-ready BED12+ transcript-model stream.

Reproduces the semantics of sesameData:::build_GENCODE_gtf (the function that
built genomeInfo$txns) but targets a plain, indexable BED instead of an R
GRangesList: one row per transcript, exon blocks, thick span = the coding
region [min CDS start, max CDS end] (empty for non-coding transcripts).

Reads a GENCODE annotation GTF (plain or .gz, or '-' for stdin) and writes an
UNSORTED BED12+ to stdout; sort + bgzip + tabix are done by the caller.

Columns (tab-separated):
   1 chrom          6 strand              11 blockSizes (CSV, no trailing comma)
   2 chromStart(0)  7 thickStart          12 blockStarts (CSV, relative to col2)
   3 chromEnd       8 thickEnd            13 gene_id
   4 transcript_id  9 itemRgb (0)         14 gene_name
   5 score (0)     10 blockCount          15 gene_type
                                          16 transcript_name
                                          17 transcript_type

Coordinates: GTF is 1-based inclusive; BED is 0-based half-open. Exon blocks are
derived from the 'exon' features (the authoritative bounds), so the BED12 block
invariant chromStart + blockStart[-1] + blockSize[-1] == chromEnd always holds.
"""
import sys
import gzip
import re

_KEYS = ("gene_id", "transcript_id", "gene_type", "gene_name",
         "transcript_type", "transcript_name")
_ATTR = {k: re.compile(k + r' "([^"]*)"') for k in _KEYS}
_TID = _ATTR["transcript_id"]


def _open(path):
    if path == "-":
        return sys.stdin
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _attrs(s):
    return {k: (m.group(1) if (m := rx.search(s)) else "") for k, rx in _ATTR.items()}


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "-"
    tx = {}          # transcript_id -> record
    order = []       # first-seen order (stable; final output is sorted anyway)

    for line in _open(src):
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        feat, start, end, attr = f[2], f[3], f[4], f[8]
        if feat == "transcript":
            a = _attrs(attr)
            tid = a["transcript_id"]
            a.update(chrom=f[0], strand=f[6], exons=[], cds_lo=None, cds_hi=None)
            tx[tid] = a
            order.append(tid)
        elif feat == "exon":
            m = _TID.search(attr)
            if m and (t := tx.get(m.group(1))) is not None:
                t["exons"].append((int(start), int(end)))
        elif feat == "CDS":
            m = _TID.search(attr)
            if m and (t := tx.get(m.group(1))) is not None:
                s, e = int(start), int(end)
                t["cds_lo"] = s if t["cds_lo"] is None else min(t["cds_lo"], s)
                t["cds_hi"] = e if t["cds_hi"] is None else max(t["cds_hi"], e)

    out = sys.stdout
    n = 0
    for tid in order:
        t = tx[tid]
        exons = sorted(t["exons"])
        if not exons:
            continue                     # transcript with no exon feature: skip
        chrom_start = exons[0][0] - 1     # 0-based
        chrom_end = max(e for _, e in exons)
        sizes = [e - s + 1 for s, e in exons]
        starts = [(s - 1) - chrom_start for s, _ in exons]
        if t["cds_lo"] is None:
            thick_start = thick_end = chrom_end            # non-coding: empty thick
        else:
            thick_start = t["cds_lo"] - 1
            thick_end = t["cds_hi"]
        out.write("\t".join(map(str, (
            t["chrom"], chrom_start, chrom_end, tid, 0, t["strand"],
            thick_start, thick_end, 0, len(exons),
            ",".join(map(str, sizes)), ",".join(map(str, starts)),
            t["gene_id"], t["gene_name"], t["gene_type"],
            t["transcript_name"], t["transcript_type"]))) + "\n")
        n += 1
    print("build_genes: %d transcripts" % n, file=sys.stderr)


if __name__ == "__main__":
    main()
