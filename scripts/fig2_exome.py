import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import beta

PAL={"glycosyltransferase":("#2a78d6","o","Glycosylation enzymes","Glycosylation"),
     "Golgi control layer":("#1baf7a","s","Golgi control layer","Golgi control"),
     "kinase":("#eb6834","^","Protein kinases","Kinases"),
     "receptor tyrosine kinase":("#eda100","D","Receptor tyrosine kinases","RTKs")}
ORD=["receptor tyrosine kinase","kinase","Golgi control layer","glycosyltransferase"]
SURF,INK,INK2,GRID="#fcfcfb","#0b0b0b","#52514e","#dcdbd6"

cls={}
for c in ORD:
    for l in open(c.replace(" ","_")+".txt"):
        s=l.strip()
        if s: cls.setdefault(s,c)
g=pd.read_csv("exome_gene.csv")
g["n"]=g[["n_syn","n_mis","n_non","n_spl"]].sum(axis=1)
g=g[g.n>0].copy()
g["classe"]=g.gene_name.map(cls)
sub=g[g.classe.notna()].copy()
sub["q"]=sub.qallsubs_cv.clip(lower=1e-17); sub["w"]=sub.wmis_cv.clip(0.05,30)
cd=pd.read_csv("exome_class_dnds.csv"); cd=cd[cd.stat=="wall"].set_index("classe")
cg=pd.read_csv("cgc_by_class.csv").set_index("gene_class")
NSIG=int((g.qallsubs_cv<0.05).sum()); NTEST=len(g)

fig=plt.figure(figsize=(7.2,10.2),dpi=300); fig.patch.set_facecolor(SURF)
gs=fig.add_gridspec(3,2,height_ratios=[1.42,0.86,0.74],hspace=0.46,wspace=0.26)
axA=fig.add_subplot(gs[0,:]); axB=fig.add_subplot(gs[1,0])
axC=fig.add_subplot(gs[1,1]); axD=fig.add_subplot(gs[2,:])
for a in (axA,axB,axC,axD): a.set_facecolor(SURF)

# A
axA.axhline(-np.log10(0.05),color="#c01c28",lw=1.0,ls=(0,(4,3)),zorder=1)
axA.text(0.055,-np.log10(0.05)," q = 0.05",color="#c01c28",fontsize=7.5,va="bottom")
axA.axvline(1,color=INK2,lw=0.9,ls=(0,(3,2)),zorder=1)
for c in ["kinase","receptor tyrosine kinase","Golgi control layer","glycosyltransferase"]:
    s=sub[sub.classe==c]; col,m,_,_=PAL[c]
    axA.scatter(s.w,-np.log10(s.q),s=25,c=col,marker=m,alpha=.8,lw=.6,
                edgecolors=SURF,zorder=3)
for nm,dx,dy in [("BRAF",-34,2),("STK11",7,2),("ATM",7,-9),("MAP2K4",7,2),
                 ("MAP3K1",-44,2),("EGFR",7,2),("FGFR2",7,-9),("ERBB2",7,2),
                 ("AKT1",7,-9),("ERBB3",-40,-9)]:
    r=sub[sub.gene_name==nm]
    if len(r):
        r=r.iloc[0]
        axA.annotate(nm,(r.w,-np.log10(r.q)),textcoords="offset points",
                     xytext=(dx,dy),fontsize=7,color=INK2,zorder=5)
axA.set(xscale="log",xlim=(0.045,40),ylim=(-0.6,26))
axA.set_xlabel("dN/dS for missense substitutions",fontsize=9,color=INK)
axA.set_ylabel("$-\\log_{10}$ q",fontsize=9,color=INK)
axA.grid(True,color=GRID,lw=.6); axA.set_axisbelow(True)
for s_ in ("top","right"): axA.spines[s_].set_visible(False)
for s_ in ("left","bottom"): axA.spines[s_].set_color(GRID)
axA.tick_params(colors=INK2,labelsize=8,length=3)
axA.legend(handles=[Line2D([],[],marker=PAL[c][1],color="none",mfc=PAL[c][0],
           mec="none",ms=6,label=PAL[c][2]) for c in ORD],loc="upper left",
           frameon=True,fontsize=8,facecolor=SURF,edgecolor=GRID,labelcolor=INK,
           handletextpad=.5,borderpad=.6)
axA.set_title("A    Selection acts on the kinases, not on the glycosylation machinery",
              fontsize=10.5,color=INK,loc="left",pad=16,fontweight="bold")
axA.text(0,1.015,f"TCGA MC3, 10,168 tumours; dNdScv fitted across {NTEST:,} genes "
         f"exome-wide, {NSIG} significant",transform=axA.transAxes,fontsize=8,color=INK2)

# B
for i,c in enumerate(ORD):
    s=sub[sub.classe==c]; k=int((s.q<0.05).sum()); n=len(s)
    lo,hi=beta.ppf(.025,k+.5,n-k+.5),beta.ppf(.975,k+.5,n-k+.5)
    axB.barh(i,100*k/n,height=.55,color=PAL[c][0],zorder=3)
    axB.plot([100*lo,100*hi],[i,i],color=INK,lw=1.2,zorder=4)
    axB.text(100*hi+0.9,i,f"{k}/{n}",fontsize=8,color=INK2,va="center")
axB.set_xlim(0,27)
axB.set_xlabel("Under positive selection,\nq < 0.05 (% of class)",fontsize=8.5,color=INK)
axB.set_title("B    dN/dS, gene by gene",fontsize=10,color=INK,loc="left",
              pad=10,fontweight="bold")

# C
for i,c in enumerate(ORD):
    r=cg.loc[c]; n=r.n_class
    seq=100*r.n_seq/n; oth=100*(r.n_census-r.n_seq)/n
    axC.barh(i,seq,height=.55,color=PAL[c][0],zorder=3)
    axC.barh(i,oth,left=seq,height=.55,color=PAL[c][0],alpha=.32,zorder=3)
    axC.text(seq+oth+1.4,i,f"{int(r.n_census)}/{int(n)}",fontsize=8,color=INK2,va="center")
axC.set_xlim(0,62)
axC.set_xlabel("In the Cancer Gene Census\n(% of class)",fontsize=8.5,color=INK)
axC.set_title("C    COSMIC Census v104",fontsize=10,color=INK,loc="left",
              pad=10,fontweight="bold")
axC.legend(handles=[Patch(facecolor=INK2,label="point mutation"),
                    Patch(facecolor=INK2,alpha=.32,label="translocation or\ncopy number only")],
           loc="lower right",frameon=True,fontsize=7,facecolor=SURF,edgecolor=GRID,
           labelcolor=INK,handlelength=1.1,handletextpad=.5,borderpad=.5)

for ax in (axB,axC):
    ax.set_ylim(len(ORD)-.5,-.7); ax.set_yticks(range(len(ORD)))
    ax.grid(True,axis="x",color=GRID,lw=.6); ax.set_axisbelow(True)
    for s_ in ("top","right","left"): ax.spines[s_].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=INK2,labelsize=8,length=3); ax.tick_params(axis="y",length=0)
axB.set_yticklabels([PAL[c][3] for c in ORD],fontsize=8.5); axC.set_yticklabels([])

# D
axD.axvline(1,color=INK2,lw=1.0,ls=(0,(4,3)),zorder=1)
for i,c in enumerate(ORD):
    r=cd.loc[c]
    axD.plot([r.lo,r.hi],[i,i],color=PAL[c][0],lw=2.2,solid_capstyle="round",zorder=3)
    axD.plot([r.dNdS],[i],"o",ms=6.5,color=PAL[c][0],mec=SURF,mew=1,zorder=4)
    lab=f"{r.dNdS:.3f}" + ("" if r.p<0.05 else "   n.s.")
    axD.text(r.hi+0.003,i,lab,fontsize=8,color=INK2,va="center")
axD.set_ylim(len(ORD)-.5,-.7); axD.set_yticks(range(len(ORD)))
axD.set_yticklabels([PAL[c][3] for c in ORD],fontsize=8.5)
axD.set_xlim(0.982,1.185)
axD.set_xlabel("dN/dS pooled over the whole class",fontsize=9,color=INK)
axD.text(1,-0.63,"  neutral",color=INK2,fontsize=7,va="bottom",ha="left")
axD.grid(True,axis="x",color=GRID,lw=.6); axD.set_axisbelow(True)
for s_ in ("top","right","left"): axD.spines[s_].set_visible(False)
axD.spines["bottom"].set_color(GRID)
axD.tick_params(colors=INK2,labelsize=8,length=3); axD.tick_params(axis="y",length=0)
axD.set_title("D    Pooled over the class, glycosylation is indistinguishable from neutral",
              fontsize=10,color=INK,loc="left",pad=10,fontweight="bold")

fig.subplots_adjust(left=0.135,right=0.965,top=0.955,bottom=0.058)
fig.savefig("Figure2.pdf",facecolor=SURF); fig.savefig("Figure2.png",facecolor=SURF)
print("significativi esoma:",NSIG,"su",NTEST)
for c in ORD:
    s=sub[sub.classe==c]; print(f"  {c:26s} {int((s.q<0.05).sum())}/{len(s)}")
