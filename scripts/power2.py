# -*- coding: utf-8 -*-
"""
power2.py — detectable effect size, calibrated on the analysis itself.

power.py computed the Poisson power of the test. It is anti-conservative:
it treats the neutral expectation as known, whereas dNdScv fits it with a
negative binomial whose dispersion adds real uncertainty. Five glycosylation
genes cleared that Poisson threshold and none is significant, which is the
proof that the approximation is too generous.

This version makes no distributional assumption. It reads the detection
boundary off the 20,091 genes that were actually tested: for a gene with a
given neutral expectation, it asks what missense dN/dS genes of comparable
expectation actually needed in order to reach significance.

For each gene, among the genes whose expected missense count lies within a
factor of 1.3, it takes the smallest observed dN/dS that reached
missense-driven significance (BH q < 0.05 on dNdScv's own missense p value,
restricted to genes with dN/dS above 1 so that tumour suppressors detected
through truncating mutations do not set the bar).

That number is the empirical detection threshold: an effect at or above it,
in a gene of that size and mutation burden, was in fact detected somewhere
in this analysis.
"""
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

R = ""   # scripts read and write the current directory
gm = pd.read_csv(R + "genemuts.csv")
sel = pd.read_csv(R + "exome_gene.csv")
d = gm.merge(sel[["gene_name", "wmis_cv", "pmis_cv", "qallsubs_cv"]], on="gene_name")
d = d[d.exp_mis > 0].copy()
d["qmis"] = multipletests(d.pmis_cv.fillna(1), method="fdr_bh")[1]

rilevati = d[(d.qmis < 0.05) & (d.wmis_cv > 1)].sort_values("exp_mis")
print("genes detected on missense (BH q<0.05, dN/dS>1): %d" % len(rilevati))
print("their expected-missense range: %.1f to %.1f\n"
      % (rilevati.exp_mis.min(), rilevati.exp_mis.max()))

e_det = rilevati.exp_mis.values
w_det = rilevati.wmis_cv.values
FATT = 1.3


def soglia(e):
    """smallest dN/dS actually detected among genes of comparable expectation"""
    m = (e_det >= e / FATT) & (e_det <= e * FATT)
    if m.sum() < 3:
        return np.nan
    return w_det[m].min()


d["soglia"] = [soglia(e) for e in d.exp_mis]

classi = {}
for c in ["glycosyltransferase", "Golgi_control_layer", "kinase",
          "receptor_tyrosine_kinase"]:
    for g in open(R + c + ".txt").read().split():
        classi[g] = c.replace("_", " ")
d["classe"] = d.gene_name.map(classi)
d["ntot"] = d[["n_syn", "n_mis", "n_non", "n_spl"]].sum(axis=1)
cl = d[d.classe.notna() & (d.ntot > 0)].copy()

ordine = ["receptor tyrosine kinase", "kinase", "Golgi control layer",
          "glycosyltransferase"]
print("EMPIRICAL DETECTION THRESHOLD (dN/dS that comparable genes needed)\n")
print("%-26s %5s %9s %14s %14s" % ("class", "genes", "median", "observed median",
                                   "observed max"))
for c in ordine:
    s = cl[cl.classe == c]
    t = s.soglia.dropna()
    print("%-26s %5d %9.2f %14.3f %14.2f"
          % (c, len(s), t.median() if len(t) else float("nan"),
             s.wmis_cv.median(), s.wmis_cv.max()))

print("\nGENES WHOSE OBSERVED dN/dS REACHES THEIR OWN EMPIRICAL THRESHOLD\n")
for c in ordine:
    s = cl[(cl.classe == c)].dropna(subset=["soglia"])
    sopra = s[s.wmis_cv >= s.soglia]
    print("%-26s %3d / %3d%s" % (c, len(sopra), len(s),
          ("  -> " + ", ".join(sopra.nlargest(5, "wmis_cv").gene_name)) if len(sopra) else ""))

print("\nSANITY CHECK: do the genes above threshold coincide with the significant ones?")
for c in ordine:
    s = cl[(cl.classe == c)].dropna(subset=["soglia"])
    sopra = set(s[s.wmis_cv >= s.soglia].gene_name)
    sig = set(s[s.qallsubs_cv < 0.05].gene_name)
    print("  %-26s above threshold %2d, actually significant %2d, overlap %2d"
          % (c, len(sopra), len(sig), len(sopra & sig)))

cl[["gene_name", "classe", "n_syn", "n_mis", "exp_mis", "wmis_cv", "soglia",
    "qallsubs_cv"]].rename(columns={
        "gene_name": "gene", "classe": "gene_class",
        "exp_mis": "expected_missense_neutral",
        "wmis_cv": "observed_dNdS_missense",
        "soglia": "empirical_detection_threshold",
        "qallsubs_cv": "q_global"}).sort_values(
    ["gene_class", "observed_dNdS_missense"], ascending=[True, False]).to_csv(
    "power_by_gene.csv", index=False)
print("\nwritten power_by_gene.csv")
