# -*- coding: utf-8 -*-
"""
Extracts every single-base substitution from the MC3 MAF, with no gene filter.

The background mutation rate in dNdScv is estimated from whatever genes it is
given. Supplying only the 964 genes of interest means the background is fitted
on a set that is not a random sample of the genome, which biases the per-gene
dN/dS. The correct design is the whole exome for the background, with the gene
classes applied to the output afterwards.

Downloads nothing: it reads the mc3.v0.2.8.PUBLIC.maf.gz already in place
(get_mc3.py fetches it).

    python extract_all.py

Writes, beside the MAF:
    mutations_all.tsv.gz    sampleID, chr, pos, ref, mut

Expected: 3,284,873 substitutions in 10,294 tumours. Pass the result to
filter_hypermutators.py, which removes the MSI and POLE cases.

Note on chromosomes: MC3 encodes X as 23 and Y as 24, which dNdScv does not
recognise. They are recoded here. Without this, 5,937 mutations are silently
dropped, including those in ATP7A, SLC35A2, C1GALT1C1 and PIGA.
"""
import csv, gzip, os, sys

MAF = "mc3.v0.2.8.PUBLIC.maf.gz"
BASI = {"A", "C", "G", "T"}


def trova():
    for d in (".", "fig1", "..", os.path.expanduser("~")):
        p = os.path.join(d, MAF)
        if os.path.exists(p):
            return p
    sys.exit("Cannot find %s. Put this script in the same folder as the MAF." % MAF)


def main():
    maf = trova()
    out = os.path.join(os.path.dirname(os.path.abspath(maf)), "mutations_all.tsv.gz")
    print("MAF:", os.path.abspath(maf))
    print("%.0f MB\n" % (os.path.getsize(maf) / 1e6))

    n = tenute = 0
    viste = set()
    with gzip.open(maf, "rt", encoding="utf-8", errors="replace") as fh, \
         gzip.open(out, "wt", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(["sampleID", "chr", "pos", "ref", "mut"])
        while True:
            prima = fh.readline()
            if not prima.startswith("#"):
                break
        campi = prima.rstrip("\n").split("\t")
        for m in csv.DictReader(fh, fieldnames=campi, delimiter="\t"):
            n += 1
            if n % 500000 == 0:
                print("  %d rows read, %d kept" % (n, tenute))
            ref = (m.get("Reference_Allele") or "").strip()
            alt = (m.get("Tumor_Seq_Allele2") or "").strip()
            if ref not in BASI or alt not in BASI or ref == alt:
                continue
            ch = str(m.get("Chromosome", "")).replace("chr", "")
            if ch == "23": ch = "X"
            if ch == "24": ch = "Y"
            pos = m.get("Start_Position") or m.get("Start_position")
            sid = (m.get("Tumor_Sample_Barcode") or "")[:15]
            if not (ch and pos and sid):
                continue
            k = (sid, ch, pos, ref, alt)
            if k in viste:
                continue
            viste.add(k)
            w.writerow(k)
            tenute += 1

    print("\n%d rows read from the MAF" % n)
    print("written %s" % out)
    print("%d substitutions, %.0f MB" % (tenute, os.path.getsize(out) / 1e6))
    print("\nnext: python filter_hypermutators.py mutations_all.tsv.gz")


main()
