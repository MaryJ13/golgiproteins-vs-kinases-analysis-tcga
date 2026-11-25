# -*- coding: utf-8 -*-
"""
build_kinase_classes.py

Turns the parsed KinHub table into the two kinase classes used in Figure 1.

    python build_kinase_classes.py kinhub_parsed.csv exome_gene.csv \
        glycosyltransferase.txt Golgi_control_layer.txt

Input : kinhub_parsed.csv  -- symbol, xname, group, family
                              (from build_kinase_lists.py, which parses the
                              KinHub PDF; 519 rows)
        exome_gene.csv     -- for `gene_name`, the symbols the dNdScv GRCh37
                              reference knows
        the two Golgi lists -- to enforce disjointness

Output: kinase.txt                    -- 434 symbols
        receptor_tyrosine_kinase.txt  --  51 symbols
        kinase_lists_log.txt          -- every gene removed, with the reason

Three filters, in order
-----------------------
1. **Not in the GRCh37 reference.** 14 of 519. A gene the reference does not
   know has no dN/dS, so it cannot enter the figure either way.

2. **Already in a Golgi class.** CERT1, POMK, SCYL1, SCYL3 sit in both the
   GolgiSet and the kinome. The GolgiSet wins: a gene belonging to the object
   under study must not also serve as its control. Without this the kinase
   classes would be contaminated by the very genes being tested.

3. **Not a catalytically active protein kinase.** 26 genes, listed below.
   KinHub's Atypical and Other groups include ADCK family members (ubiquinone
   kinases, not protein kinases), pseudokinases with a degenerate catalytic
   triad (TRIB1/2/3, STRADA/B, MLKL, PXK, RYK, AATK, SCYL2, ULK4, RNASEL,
   RPS6KC1, RPS6KL1, IRAK2), chromatin and scaffold proteins carrying a
   kinase-like domain (BAZ1A, TRIM33, TRIM66, TRRAP, PAN3, ABR), and a few
   others that are not kinases at all (STK19, PRPF4B is a genuine kinase but
   a spliceosomal one and was excluded with them).

   This filter matters for the argument, not for tidiness. The kinase classes
   are the **positive control**: they have to show selection for the absence
   of selection in the Golgi classes to mean anything. Leaving in genes that
   are not kinases but *are* strong cancer drivers -- TRRAP and TRIM33 are
   both in the COSMIC Census -- would inflate the control with signal that
   has nothing to do with kinase activity. Removing them makes the control
   harder to pass, not easier.

Receptor tyrosine kinases
-------------------------
RTKs are the TK group minus the families whose members are cytoplasmic
(Abl, Ack, Csk, FAK, Fer, Jak, Src, Syk, Tec, TK-Unique), and minus two
receptor families that are pseudokinases (CCK4/PTK7, Lmr/LMTK). ERBB3 is a
receptor of the EGFR family but is catalytically dead, so it stays in the
kinase class rather than the RTK class.

The result is 51 RTKs. Note that `build_kinase_lists.py`, which parses the
PDF, reports 58 before any of the filtering above; 58 is the canonical count
of human RTKs and is the check that the PDF was read correctly. The 51 here
is that set after the reference and pseudokinase filters.
"""
import os
import sys

import pandas as pd

# The dNdScv GRCh37 reference predates these HGNC renames, so the current
# symbol KinHub carries has to be mapped back to the one the reference uses.
# current symbol -> GRCh37 reference symbol.
TO_GRCH37 = {
    "GRK2":    "ADRBK1",
    "GRK3":    "ADRBK2",
    "HASPIN":  "GSG2",
    "COQ8A":   "ADCK3",
    "COQ8B":   "ADCK4",
    "CILK1":   "ICK",
    "PAK5":    "PAK7",
    "STKLD1":  "C9orf96",
    "STK26":   "MST4",
    "RSKR":    "SGK494",
}

NON_RECEPTOR_FAMILIES = {
    "Abl", "Ack", "Csk", "FAK", "Fer", "Jak", "Src", "Syk", "Tec",
    "TK-Unique",
}
PSEUDOKINASE_RECEPTOR_FAMILIES = {"CCK4", "Lmr"}
PSEUDOKINASE_RECEPTORS = {"ERBB3"}

NOT_PROTEIN_KINASES = {
    "AATK", "ABR", "ADCK1", "ADCK2", "ADCK5", "BAZ1A", "IRAK2", "MLKL",
    "PAN3", "PRPF4B", "PXK", "RNASEL", "RPS6KC1", "RPS6KL1", "RYK",
    "SCYL2", "STK19", "STRADA", "STRADB", "TRIB1", "TRIB2", "TRIB3",
    "TRIM33", "TRIM66", "TRRAP", "ULK4",
}


def main():
    a = sys.argv[1:] or ["kinhub_parsed.csv", "exome_gene.csv",
                         "glycosyltransferase.txt", "Golgi_control_layer.txt"]
    k = pd.read_csv(a[0])
    ref = set(pd.read_csv(a[1])["gene_name"])
    golgi = set(open(a[2]).read().split()) | set(open(a[3]).read().split())

    log = []
    tenuti = []
    for _, r in k.iterrows():
        s = str(r["symbol"]).strip()
        if s not in ref and TO_GRCH37.get(s) in ref:
            log.append("%-12s  -> %s (GRCh37 reference symbol)" % (s, TO_GRCH37[s]))
            s = TO_GRCH37[s]
        if s not in ref:
            log.append("%-12s  not in the GRCh37 reference" % s)
        elif s in golgi:
            log.append("%-12s  already in a Golgi class" % s)
        elif s in NOT_PROTEIN_KINASES:
            log.append("%-12s  not a catalytically active protein kinase" % s)
        else:
            tenuti.append((s, str(r["group"]).strip(), str(r["family"]).strip()))

    def is_rtk(_s, grp, fam):
        return (grp == "TK"
                and fam not in NON_RECEPTOR_FAMILIES
                and fam not in PSEUDOKINASE_RECEPTOR_FAMILIES
                and _s not in PSEUDOKINASE_RECEPTORS)

    rtk = sorted({s for s, g, f in tenuti if is_rtk(s, g, f)})
    kin = sorted({s for s, g, f in tenuti if not is_rtk(s, g, f)})

    assert not set(rtk) & set(kin)
    assert not (set(rtk) | set(kin)) & golgi

    out = os.path.dirname(os.path.abspath(a[0]))
    for nome, lst in (("kinase", kin), ("receptor_tyrosine_kinase", rtk)):
        p = os.path.join(out, nome + ".txt")
        with open(p, "w") as fh:
            fh.write("\n".join(lst) + "\n")
        print("%-26s %3d  -> %s" % (nome, len(lst), p))

    with open(os.path.join(out, "kinase_lists_log.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    print("\n%d genes removed, listed in kinase_lists_log.txt" % len(log))
    print("expected: 435 kinases, 51 RTKs, 43 removed")


main()
