# -*- coding: utf-8 -*-
"""
filter_hypermutators.py

Removes hypermutated tumours from the exome-wide substitution table, then
writes the dNdScv input.

    python filter_hypermutators.py mutations_all.tsv.gz

Input : mutations_all.tsv.gz  -- from extract_all.py
                                 sampleID, chr, pos, ref, mut  (with header)
Output: exome_filt.tsv        -- same columns, same header
        hypermutators.txt     -- the barcodes removed, for the record

Why
---
MSI and POLE-mutant tumours carry tens of thousands of substitutions each.
dNdScv fits a single background mutation model across all samples, and a
handful of samples contributing a thousand times the median distorts it.
They are also the samples that dominate any per-gene count, so a gene can
look recurrently mutated on the strength of five hypermutators.

The threshold is 3,000 coding substitutions per tumour. The distribution is
strongly bimodal there, so the exact cut is not delicate: anywhere between
2,000 and 5,000 removes the same cluster of samples.

Expected on the MC3 public call set:
    read      3,284,873 substitutions in 10,294 tumours
    removed     126 tumours
    kept      2,096,036 substitutions in 10,168 tumours
"""
import collections
import csv
import gzip
import os
import sys

SOGLIA = 3000


def apri(p, modo):
    return gzip.open(p, modo, encoding="utf-8", newline="") if p.endswith(".gz") \
        else open(p, modo, encoding="utf-8", newline="")


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "mutations_all.tsv.gz"
    if not os.path.exists(src):
        sys.exit("non trovo %s -- esegui prima extract_all.py" % src)

    d = os.path.dirname(os.path.abspath(src))
    out = os.path.join(d, "exome_filt.tsv")
    elenco = os.path.join(d, "hypermutators.txt")

    # pass 1 -- count substitutions per tumour
    conta = collections.Counter()
    n = 0
    with apri(src, "rt") as fh:
        r = csv.reader(fh, delimiter="\t")
        testa = next(r)
        for riga in r:
            n += 1
            conta[riga[0]] += 1
    print("lette %d sostituzioni in %d campioni" % (n, len(conta)))

    scartati = sorted(s for s, c in conta.items() if c > SOGLIA)
    print("scarto %d campioni con piu' di %d sostituzioni" % (len(scartati), SOGLIA))
    if scartati:
        print("  il piu' alto: %s con %d" % (max(conta, key=conta.get),
                                             max(conta.values())))
    scarta = set(scartati)

    with open(elenco, "w") as fh:
        for s in scartati:
            fh.write("%s\t%d\n" % (s, conta[s]))

    # pass 2 -- write what survives
    tenute = 0
    campioni = set()
    with apri(src, "rt") as fh, open(out, "w", newline="") as fo:
        r = csv.reader(fh, delimiter="\t")
        w = csv.writer(fo, delimiter="\t")
        next(r)
        w.writerow(testa)
        for riga in r:
            if riga[0] in scarta:
                continue
            w.writerow(riga)
            campioni.add(riga[0])
            tenute += 1

    print("\nscritto %s" % out)
    print("%d sostituzioni in %d campioni" % (tenute, len(campioni)))
    print("elenco degli scartati in %s" % elenco)


main()
