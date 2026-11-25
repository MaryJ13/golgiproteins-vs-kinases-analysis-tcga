# -*- coding: utf-8 -*-
"""
build_source_data.py

Joins the dNdScv output to the four gene classes and to the COSMIC Cancer
Gene Census, producing the two tables the figure is drawn from.

    python build_source_data.py exome_gene.csv \
        Cosmic_CancerGeneCensus_v104_GRCh37.tsv.gz

Input : exome_gene.csv   -- dNdScv sel_cv, all 20,091 reference genes
        the four class .txt files, in the same directory
        the COSMIC Census TSV (gzipped, as downloaded)

Output: Figure1_source_data.csv  -- one row per class gene; Supplementary S1
        cgc_by_class.csv         -- Census counts by class and mechanism,
                                    which is panel C

COSMIC
------
The Census is behind an academic login and cannot be redistributed, so it is
not in this package. Download `Cosmic_CancerGeneCensus_Tsv_v104_GRCh37.tar`
from https://cancer.sanger.ac.uk/cosmic/download/cosmic (free for academic
use), untar it, and pass the .tsv.gz to this script.

The GRCh38 release gives identical answers for every gene reported here,
because only GENE_SYMBOL, TIER, ROLE_IN_CANCER and MUTATION_TYPES are used,
never coordinates. Both were checked.

MUTATION_TYPES is a comma-separated code list. The three groups used in
panel C:

    point mutation   Mis (missense), N (nonsense), F (frameshift),
                     S (splice site), O (other)
    translocation    T
    amplification    A
    large deletion   D

A and D are counted in separate columns. They are both copy-number events but
they point in opposite directions -- A marks a dosage-driven oncogene, D a
deleted tumour suppressor -- so a single "copy number" column is ambiguous.
Among the kinases, for instance, 6 genes carry A (AKT2, AKT3, PAK2, PRKCI,
RAF1, ROCK1) and a different 6 carry D (ACVR1B, ATM, LATS2, MAP2K4, STK11,
TGFBR1). Only the first six mean what "amplified" means in the review.

Panel C of the figure uses only n_census and n_seq; the other columns are
there for the text.

Genes with no coding mutations at all
-------------------------------------
dNdScv reports dN/dS = 0 for a gene in which it saw no substitutions, which
is not an estimate of anything. Those genes are dropped here, as they are in
the figure script. There are three: OST4 and GOLGA7B (both very short) and
SBK3.

These codes were not taken on trust. They were checked against genes whose
mechanism is known independently: MDM2 carries A alone, MYC carries A and T,
EGFR carries A, O and Mis, TP53 carries Mis, N, F, T and D. A gene can carry
several, so the three counts in cgc_by_class.csv are not mutually exclusive
and do not sum to n_census.
"""
import os
import sys

import pandas as pd

CLASSI = ["glycosyltransferase", "Golgi_control_layer",
          "kinase", "receptor_tyrosine_kinase"]

PUNTIFORMI = {"Mis", "N", "F", "S", "O"}


def codici(v):
    if not isinstance(v, str):
        return set()
    return {c.strip() for c in v.split(",") if c.strip()}


def main():
    a = sys.argv[1:] or ["exome_gene.csv",
                         "Cosmic_CancerGeneCensus_v104_GRCh37.tsv.gz"]
    d = os.path.dirname(os.path.abspath(a[0]))

    g = pd.read_csv(a[0])
    cgc = pd.read_csv(a[1], sep="\t", low_memory=False)
    cgc = cgc.drop_duplicates("GENE_SYMBOL").set_index("GENE_SYMBOL")

    # class membership
    mappa = {}
    for cl in CLASSI:
        for s in open(os.path.join(d, cl + ".txt")).read().split():
            assert s not in mappa, "%s is in two classes" % s
            mappa[s] = cl.replace("_", " ")

    g["n_tot"] = g[["n_syn", "n_mis", "n_non", "n_spl"]].sum(axis=1)
    muti = g[g.n_tot > 0]
    senza = sorted(set(g.gene_name[g.n_tot == 0]) & set(mappa))
    if senza:
        print("%d class genes with no coding mutations, dropped: %s"
              % (len(senza), ", ".join(senza)))

    s = muti[muti.gene_name.isin(mappa)].copy()
    s["gene_class"] = s.gene_name.map(mappa)
    mancanti = set(mappa) - set(g.gene_name) - set(senza)
    if mancanti:
        print("%d class genes absent from the dNdScv output: %s"
              % (len(mancanti), ", ".join(sorted(mancanti))))

    out = pd.DataFrame({
        "gene":             s.gene_name,
        "gene_class":       s.gene_class,
        "n_syn":            s.n_syn,
        "n_mis":            s.n_mis,
        "n_non":            s.n_non,
        "n_spl":            s.n_spl,
        "dNdS_missense":    s.wmis_cv.round(3),
        "dNdS_truncating":  s.wnon_cv.round(3),
        "q_global":         s.qallsubs_cv,
    })
    inc = out.gene.isin(cgc.index)
    out["CGC_v104"] = ["yes" if b else "no" for b in inc]
    out["CGC_tier"] = out.gene.map(cgc["TIER"])
    out["CGC_mutation_types"] = out.gene.map(cgc["MUTATION_TYPES"])
    out["CGC_role"] = out.gene.map(cgc["ROLE_IN_CANCER"])

    out = out.sort_values(["gene_class", "q_global"])
    p = os.path.join(d, "Figure1_source_data.csv")
    out.to_csv(p, index=False)
    print("%d genes -> %s" % (len(out), p))

    # panel C
    righe = []
    for cl in CLASSI:
        nome = cl.replace("_", " ")
        sub = out[out.gene_class == nome]
        inn = sub[sub.CGC_v104 == "yes"]
        mt = inn.CGC_mutation_types.map(codici)
        righe.append({
            "gene_class":    nome,
            "n_class":       len(sub),
            "n_census":      len(inn),
            "n_tier1":       int((inn.CGC_tier == 1).sum()),
            "n_tier2":       int((inn.CGC_tier == 2).sum()),
            "n_seq":            int(sum(bool(c & PUNTIFORMI) for c in mt)),
            "n_transloc":       int(sum("T" in c for c in mt)),
            "n_amplification":  int(sum("A" in c for c in mt)),
            "n_large_deletion": int(sum("D" in c for c in mt)),
        })
    c = pd.DataFrame(righe)
    p = os.path.join(d, "cgc_by_class.csv")
    c.to_csv(p, index=False)
    print("\n%s" % p)
    print(c.to_string(index=False))

    # the headline, printed so a rerun can be checked at a glance
    print("\nunder significant selection (q < 0.05), exome-wide family:")
    for cl in CLASSI:
        nome = cl.replace("_", " ")
        sub = out[out.gene_class == nome]
        print("  %-26s %d/%d" % (nome, (sub.q_global < 0.05).sum(), len(sub)))


main()
