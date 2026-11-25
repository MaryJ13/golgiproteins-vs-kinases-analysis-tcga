# -*- coding: utf-8 -*-
"""
build_golgi_lists.py

Builds the two Golgi gene lists from the curated GolgiSet.

    python build_golgi_lists.py GolgiSet.xlsx exome_gene.csv

Input : GolgiSet.xlsx   -- 507 rows, curated in-house. The column `Function`
                           carries one of four values:
                             Glycosylation enzyme              234
                             Membrane trafficking component    198
                             Channels and transporters          40
                             lipid trafficking and metabolism   35
        exome_gene.csv  -- used only for its `gene_name` column, i.e. the
                           list of symbols the dNdScv GRCh37 reference knows

Output: glycosyltransferase.txt   -- the glycosylation enzymes
        Golgi_control_layer.txt   -- the other three functions together
        golgi_lists_log.txt       -- every substitution and every drop

Symbols
-------
The dNdScv GRCh37 reference predates several HGNC renames, so a symbol that
is current in the GolgiSet may not exist in the reference. For each row the
script tries, in order, `current_gene_symbol`, then the original
`Gene symbol`, then the explicit alias table below. A row whose gene the
reference does not know under any of those names is dropped and logged;
there is nothing to be done about it, since a gene absent from the reference
has no dN/dS.

The alias table maps current symbol -> the symbol the GRCh37 reference uses.
Each was checked against HGNC.
"""
import os
import sys

import pandas as pd

ALIAS = {
    "GALNT17":  "WBSCR17",
    "LARGE1":   "LARGE",
    "B3GLCT":   "B3GALTL",
    "PXYLP1":   "ACPL2",
    "B4GAT1":   "B3GNT1",
    "DOP1A":    "DOPEY1",
    "DOP1B":    "DOPEY2",
    "VPS26C":   "DSCR3",
    "RIC1":     "KIAA1432",
    "TRAPPC2B": "TRAPPC2P1",
}

GLY = "Glycosylation enzyme"


def main():
    xlsx = sys.argv[1] if len(sys.argv) > 1 else "GolgiSet.xlsx"
    reff = sys.argv[2] if len(sys.argv) > 2 else "exome_gene.csv"
    d = pd.read_excel(xlsx)
    ref = set(pd.read_csv(reff)["gene_name"])

    log = []
    simboli = []
    for _, r in d.iterrows():
        cur = str(r["current_gene_symbol"]).strip()
        orig = str(r["Gene symbol"]).strip()
        scelto = None
        for cand, come in ((cur, "current"), (orig, "original"), (ALIAS.get(cur), "alias")):
            if cand and cand != "nan" and cand in ref:
                scelto = cand
                if come != "current":
                    log.append("SUBSTITUTED  %-12s -> %-12s (%s)" % (cur, cand, come))
                break
        if scelto is None:
            log.append("DROPPED      %-12s  not in the GRCh37 reference" % cur)
        simboli.append(scelto)

    d["sym"] = simboli
    d = d[d.sym.notna()]

    gly = sorted(set(d.sym[d.Function == GLY]))
    oth = sorted(set(d.sym[d.Function != GLY]))

    # A gene cannot be in both classes. There is no overlap by construction --
    # `Function` is a single value per row -- but assert it rather than assume.
    assert not (set(gly) & set(oth)), sorted(set(gly) & set(oth))

    out = os.path.dirname(os.path.abspath(xlsx))
    for nome, lst in (("glycosyltransferase", gly), ("Golgi_control_layer", oth)):
        p = os.path.join(out, nome + ".txt")
        with open(p, "w") as fh:
            fh.write("\n".join(lst) + "\n")
        print("%-24s %3d  -> %s" % (nome, len(lst), p))

    with open(os.path.join(out, "golgi_lists_log.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    print("\n%d substitutions and drops, listed in golgi_lists_log.txt"
          % len(log))
    print("expected: 231 glycosylation enzymes, 270 Golgi control layer,")
    print("          6 genes dropped (ABO, ARFGEF3, LARGE2, MGAT4D, RAB7B, SEC22B)")


main()
