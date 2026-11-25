# -*- coding: utf-8 -*-
"""Supplementary Figure S1: tumour-type stratification and detection power."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from statsmodels.stats.multitest import multipletests

R = ""   # scripts read and write the current directory
PAL = {"glycosyltransferase": ("#2a78d6", "o", "Glycosylation enzymes"),
       "Golgi control layer": ("#1baf7a", "s", "Golgi control layer"),
       "kinase":              ("#eb6834", "^", "Protein kinases"),
       "receptor tyrosine kinase": ("#eda100", "D", "Receptor tyrosine kinases")}
ORD = ["receptor tyrosine kinase", "kinase", "Golgi control layer",
       "glycosyltransferase"]
AB = {"Skin Cutaneous Melanoma": "SKCM", "Lung adenocarcinoma": "LUAD",
      "Lung squamous cell carcinoma": "LUSC",
      "Uterine Corpus Endometrial Carcinoma": "UCEC",
      "Colon adenocarcinoma": "COAD", "Stomach adenocarcinoma": "STAD",
      "Bladder Urothelial Carcinoma": "BLCA", "Breast invasive carcinoma": "BRCA",
      "Head and Neck squamous cell carcinoma": "HNSC",
      "Cervical squamous cell carcinoma and endocervical adenocarcinoma": "CESC",
      "Liver hepatocellular carcinoma": "LIHC"}

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                     "axes.linewidth": .8, "axes.edgecolor": "#444444"})
fig = plt.figure(figsize=(7.2, 8.4))
gs = fig.add_gridspec(2, 1, height_ratios=[1.15, 1], hspace=.34,
                      left=.13, right=.97, top=.945, bottom=.065)

# ---------------------------------------------------------------- A
ax = fig.add_subplot(gs[0])
c = pd.read_csv("per_type_class.csv")
w = c[c.stat == "wall"].copy()
tipi = (w[w.classe == "kinase"].sort_values("dNdS", ascending=False).tipo.tolist())
OFF = {cl: o for cl, o in zip(ORD, [.27, .09, -.09, -.27])}

for cl in ORD:
    col, mk, lab = PAL[cl]
    s = w[w.classe == cl].set_index("tipo")
    for i, t in enumerate(tipi):
        if t not in s.index:
            continue
        r = s.loc[t]
        y = i + OFF[cl]
        sig = r.p < 0.05
        ax.plot([r.lo, r.hi], [y, y], color=col, lw=1.2, alpha=.85, zorder=2)
        ax.plot(r.dNdS, y, mk, color=col, ms=4.6 if sig else 3.8,
                mec="#222222" if sig else col, mew=.7 if sig else 0, zorder=3)
ax.axvline(1, color="#888888", lw=.9, ls="--", zorder=1)
ax.set_yticks(range(len(tipi)))
ax.set_yticklabels([AB.get(t, t) for t in tipi], fontsize=7.5)
ax.set_ylim(-.6, len(tipi) - .4)
ax.invert_yaxis()
ax.set_xlabel("dN/dS pooled over the class, within each tumour type", fontsize=8)
ax.set_xlim(.64, 1.78)
ax.set_title("A   Selection is confined to the kinases in every tumour type",
             loc="left", fontsize=9, fontweight="bold", pad=22)
for s_ in ("top", "right"):
    ax.spines[s_].set_visible(False)
ax.grid(axis="x", color="#eeeeee", lw=.7, zorder=0)
ax.set_axisbelow(True)
h = [Line2D([], [], color=PAL[c_][0], marker=PAL[c_][1], ls="-", ms=4.4,
            label=PAL[c_][2]) for c_ in ORD]
h.append(Line2D([], [], color="none", marker="o", mec="#222222", mew=.7,
                mfc="#cccccc", ms=4.6, ls="none", label="p < 0.05"))
ax.legend(handles=h, fontsize=6.9, ncol=5, loc="lower center",
          bbox_to_anchor=(.5, 1.005), frameon=False, handlelength=1.5,
          columnspacing=1.1, handletextpad=.4)
ax.text(0, -.115,
        "11 tumour types with >50,000 substitutions; each fitted exome-wide in its own run. No Golgi or\n"
        "glycosylation gene reaches significance in any single type (0/11 for both classes, and 0/11 at gene level).",
        transform=ax.transAxes, fontsize=6.6, color="#555555", va="top")

# ---------------------------------------------------------------- B
ax = fig.add_subplot(gs[1])
gm = pd.read_csv(R + "genemuts.csv")
sel = pd.read_csv(R + "exome_gene.csv")
d = gm.merge(sel[["gene_name", "wmis_cv", "pmis_cv", "qallsubs_cv"]], on="gene_name")
d = d[d.exp_mis > 0].copy()
d["qmis"] = multipletests(d.pmis_cv.fillna(1), method="fdr_bh")[1]
det = d[(d.qmis < 0.05) & (d.wmis_cv > 1)].sort_values("exp_mis")

# empirical detection boundary: rolling minimum detected dN/dS
xs = np.logspace(np.log10(det.exp_mis.min()), np.log10(det.exp_mis.max()), 60)
ys = [det.wmis_cv[(det.exp_mis >= x / 1.6) & (det.exp_mis <= x * 1.6)].min()
      for x in xs]
ys = pd.Series(ys).ffill().bfill().rolling(7, center=True, min_periods=1).mean()

classi = {}
for cfile in ["glycosyltransferase", "Golgi_control_layer", "kinase",
              "receptor_tyrosine_kinase"]:
    for g_ in open(R + cfile + ".txt").read().split():
        classi[g_] = cfile.replace("_", " ")
d["classe"] = d.gene_name.map(classi)
d["ntot"] = d[["n_syn", "n_mis", "n_non", "n_spl"]].sum(axis=1)
cl = d[d.classe.notna() & (d.ntot > 0)]

ax.fill_between(xs, ys, 40, color="#f3b0a0", alpha=.30, zorder=0)
ax.plot(xs, ys, color="#c0392b", lw=1.5, zorder=4)
for c_ in ORD:
    col, mk, lab = PAL[c_]
    s = cl[cl.classe == c_]
    ax.scatter(s.exp_mis, s.wmis_cv.clip(.3, 30), s=11, marker=mk,
               facecolor=col, edgecolor="none", alpha=.62, zorder=2, label=lab)
for nome in ["BRAF", "EGFR", "FGFR2", "B3GNT1", "YIPF6"]:
    r = cl[cl.gene_name == nome]
    if len(r):
        r = r.iloc[0]
        ax.annotate(nome, (r.exp_mis, min(r.wmis_cv, 30)),
                    textcoords="offset points", xytext=(5, 3), fontsize=6.6,
                    color="#333333")
ax.axhline(1, color="#888888", lw=.9, ls="--", zorder=1)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("Missense substitutions expected under neutrality (gene size x burden)",
              fontsize=8)
ax.set_ylabel("dN/dS for missense substitutions", fontsize=8)
ax.set_yticks([.5, 1, 2, 5, 10, 20])
ax.set_yticklabels(["0.5", "1", "2", "5", "10", "20"])
ax.set_title("B   What the analysis could have detected", loc="left",
             fontsize=9, fontweight="bold", pad=8)
for s_ in ("top", "right"):
    ax.spines[s_].set_visible(False)
ax.grid(color="#eeeeee", lw=.7, zorder=0); ax.set_axisbelow(True)
ax.legend(fontsize=6.9, loc="upper right", frameon=True, framealpha=.95,
          edgecolor="#dddddd", scatterpoints=1)
ax.text(.012, .022,
        "Red line: the weakest missense dN/dS that genes of comparable expectation actually reached significance with.\n"
        "Shaded: the detectable region. No glycosylation enzyme enters it; the highest is B3GNT1 at 1.95 (q = 0.43).",
        transform=ax.transAxes, fontsize=6.6, color="#555555", va="bottom")

fig.savefig("FigureS1.pdf")
fig.savefig("FigureS1.png", dpi=220)
print("written FigureS1.pdf / .png")
