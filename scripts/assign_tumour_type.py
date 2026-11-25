# -*- coding: utf-8 -*-
"""
assign_tumour_type.py — labels each mutation with its TCGA tumour type.

    python assign_tumour_type.py

Input : exome_filt.tsv   sampleID, chr, pos, ref, mut   (the analysis set)
        tss.tsv          TCGA tissue source site codes -> study name
Output: mut_tipo.tsv     the same five columns plus `tipo`

`per_type.R` reads the output. Run this before it.

How the mapping works
---------------------
A TCGA barcode does not encode the disease directly. It encodes the tissue
source site, the institution that contributed the sample, in characters 6 and
7 — `TCGA-`**`02`**`-0003-01`. Each source site contributed to exactly one
study, so the code determines the tumour type. `tss.tsv` is that lookup, taken
from the GDC code table at
https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables

Expected: all 2,096,036 mutations map, with no unrecognised codes. The script
stops if any code is unknown rather than silently dropping those mutations,
because a missing code would quietly shrink a tumour type.

It also prints the mutation and tumour counts per type. `per_type.R` analyses
the types above 50,000 mutations; there are 11 of them, holding 82.7% of the
analysis set.
"""
import collections
import csv
import os
import sys

MUT = "exome_filt.tsv"
TSS = "tss.tsv"
OUT = "mut_tipo.tsv"
SOGLIA = 50000


def main():
    for f in (MUT, TSS):
        if not os.path.exists(f):
            sys.exit("Cannot find %s in the current directory." % f
                     + ("  Run: gunzip exome_filt.tsv.gz" if f == MUT else ""))

    tss = {}
    with open(TSS, encoding="utf-8") as fh:
        for riga in fh:
            if "\t" not in riga:
                continue
            k, v = riga.rstrip("\n").split("\t", 1)
            tss[k.strip()] = v.strip()
    print("%d tissue source site codes loaded" % len(tss))

    conta = collections.Counter()
    campioni = collections.defaultdict(set)
    ignoti = collections.Counter()
    n = 0

    with open(MUT, encoding="utf-8", newline="") as fh, \
            open(OUT, "w", encoding="utf-8", newline="") as fo:
        r = csv.reader(fh, delimiter="\t")
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        testa = next(r)
        w.writerow(testa + ["tipo"])
        for riga in r:
            n += 1
            codice = riga[0][5:7]
            tipo = tss.get(codice)
            if tipo is None:
                ignoti[codice] += 1
                continue
            w.writerow(riga + [tipo])
            conta[tipo] += 1
            campioni[tipo].add(riga[0])

    if ignoti:
        os.remove(OUT)
        print("\nUNRECOGNISED source site codes, nothing written:")
        for c, k in ignoti.most_common(20):
            print("   %s  %d mutations" % (c, k))
        sys.exit("Check tss.tsv against the current GDC code table.")

    print("\n%d mutations, all mapped, in %d tumour types"
          % (n, len(conta)))
    print("\n%-62s %10s %8s" % ("tumour type", "mutations", "tumours"))
    for tipo, k in conta.most_common():
        segno = "  <- analysed separately" if k >= SOGLIA else ""
        print("%-62s %10d %8d%s" % (tipo[:62], k, len(campioni[tipo]), segno))

    grandi = [t for t, k in conta.items() if k >= SOGLIA]
    quota = sum(conta[t] for t in grandi)
    print("\n%d types above %d mutations, holding %.1f%% of the analysis set"
          % (len(grandi), SOGLIA, 100.0 * quota / n))
    print("\nwritten %s" % OUT)
    print("next: Rscript per_type.R")


main()
