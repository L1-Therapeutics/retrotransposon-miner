
import pandas as pd, numpy as np, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
L="/Users/smyan/Desktop/Research/Research_Projects/L1/"; OUT=L+"report/figures/"
apply_figure_style(sizes=(9,8,7))
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Helvetica","Arial","DejaVu Sans"],"savefig.dpi":600})
NAVY="#17365D"; GREY="#8C8C8C"
def bare(ax,sides=("top","right","left")):
    for s in sides: ax.spines[s].set_visible(False)
    ax.tick_params(axis="y",length=0)

# ---------------- FIGURE 2 ----------------
t2=pd.read_csv(L+"artifacts/task2_depletion.csv")
order=[("coding_exon","Coding exon"),("5UTR","5\u2032 UTR"),("3UTR","3\u2032 UTR"),
       ("promoter","Promoter"),("intron","Intron"),("intergenic","Intergenic")]
from scipy.stats import chi2
def pois_ci(k,a=0.05):
    lo = chi2.ppf(a/2,2*k)/2 if k>0 else 0.0
    return lo, chi2.ppf(1-a/2,2*k+2)/2
fig,axes=plt.subplots(1,2,figsize=(6.9,3.3),gridspec_kw=dict(width_ratios=[1.45,1],wspace=0.55))
ax=axes[0]; y=np.arange(len(order))[::-1]
for i,(k,nm) in enumerate(order):
    for null,off,mk,col,fc,msz in [("uniform",0.17,"s",GREY,"white",4.4),("gc_matched",-0.17,"o",NAVY,NAVY,5.0)]:
        r=t2[(t2.null==null)&(t2.element==k)].iloc[0]
        klo,khi=pois_ci(int(r.observed))
        ax.plot([klo/r.expected,khi/r.expected],[y[i]+off]*2,"-",color=col,lw=1.3,
                solid_capstyle="butt",zorder=2)
        ax.plot([r.obs_over_exp],[y[i]+off],marker=mk,ms=msz,mfc=fc,mec=col,mew=1.2,ls="none",zorder=3)
ax.set_yticks(y); ax.set_yticklabels([n for _,n in order]); ax.set_ylim(-0.8,len(order)-0.2)
ax.set_xscale("log"); ax.axvline(1.0,ls=(0,(5,4)),color="#555555",lw=1.0,zorder=1)
ax.set_xlim(0.009,3.2)
ax.set_xticks([0.01,0.1,0.5,1,2]); ax.set_xticklabels(["0.01","0.10","0.50","1.00","2.00"])
ax.set_xlabel("Observed / expected insertions"); bare(ax)
ax.legend(handles=[Line2D([],[],marker="s",ls="none",mfc="white",mec=GREY,mew=1.2,ms=4.4,label="Uniform null"),
                   Line2D([],[],marker="o",ls="none",color=NAVY,ms=5.0,label="GC-matched null")],
          loc="lower left",frameon=False,handletextpad=0.5,borderpad=0.2)
panel_letter(ax,"A",case="upper")
ax=axes[1]
u=t2[(t2.null=="uniform")&(t2.element=="promoter")].iloc[0]
g=t2[(t2.null=="gc_matched")&(t2.element=="promoter")].iloc[0]
ax.plot([u.exp_lo95,u.exp_hi95],[1,1],"-",color=GREY,lw=1.4)
ax.plot([u.expected],[1],marker="s",ms=4.8,mfc="white",mec=GREY,mew=1.3,ls="none")
ax.plot([g.exp_lo95,g.exp_hi95],[0,0],"-",color=NAVY,lw=1.5)
ax.plot([g.expected],[0],marker="o",ms=5.2,color=NAVY,ls="none")
ax.axvline(33,color="#555555",lw=1.0,ls=(0,(5,4)))
ax.set_yticks([1,0]); ax.set_yticklabels(["Uniform null","GC-matched null"]); ax.set_ylim(-0.8,1.8)
ax.set_xlabel("Expected promoter insertions per genome"); ax.set_xlim(18,60)
ax.annotate("observed = 33",xy=(34.2,1.60),ha="left",va="bottom",fontsize=7,color="#555555",annotation_clip=False)
bare(ax); panel_letter(ax,"B",case="upper",dx=-0.30,dy=1.06)
fig.savefig(OUT+"figure2.png",bbox_inches="tight"); plt.close(fig)

# ---------------- FIGURE 3 ----------------
sa=pd.read_csv(L+"artifacts/selfins_sense_antisense.csv").set_index(["family","category"])
hom=pd.read_csv(L+"artifacts/selfins_homology_test.csv").set_index("host_age")
deg=pd.read_csv(L+"artifacts/selfins_en_degenerate_subfamily.csv")
fig,axes=plt.subplots(1,3,figsize=(7.0,3.2),gridspec_kw=dict(wspace=0.62))

ax=axes[0]
rows=[("ALU","nested_sense","Alu, sense"),("ALU","nested_antisense","Alu, antisense"),
      ("LINE1","nested_sense","LINE-1, sense"),("LINE1","nested_antisense","LINE-1, antisense"),
      ("SVA","nested_sense","SVA, sense"),("SVA","nested_antisense","SVA, antisense")]
y=np.arange(len(rows))[::-1]
for i,(f,cat,nm) in enumerate(rows):
    r=sa.loc[(f,cat)]; filled = cat=="nested_sense"
    if r.observed>0:
        ax.plot([1,r.enrichment],[y[i]]*2,"-",color=NAVY,lw=1.0,alpha=.55,zorder=2)
        ax.plot([r.enrichment],[y[i]],marker="o",ms=5.0 if filled else 4.6,
                mfc=NAVY if filled else "white",mec=NAVY,mew=1.3,ls="none",zorder=3)
    else:
        ax.plot([0.016,1],[y[i]]*2,"-",color=GREY,lw=1.0,alpha=.55,zorder=2)
        ax.plot([0.016],[y[i]],marker="x",ms=4.8,color=GREY,mew=1.4,ls="none",zorder=3)
        ax.annotate("none observed (0 of 9)",xy=(1.35,y[i]),ha="left",va="center",
                    fontsize=7,color=GREY)
ax.set_yticks(y); ax.set_yticklabels([n for _,_,n in rows]); ax.set_ylim(-1.0,len(rows)-0.3)
ax.set_xscale("log"); ax.axvline(1,ls=(0,(5,4)),color="#555555",lw=1.0,zorder=1)
ax.set_xlim(0.010,300); ax.set_xticks([0.1,1,10,100]); ax.set_xticklabels(["0.1","1","10","100"])
ax.set_xlabel("Observed / expected"); bare(ax)
panel_letter(ax,"A",case="upper")

ax=axes[1]
hk=[("AluY (youngest)","AluY"),("AluS (intermediate)","AluS"),("AluJ (oldest)","AluJ")]
x=np.arange(3)
for i,(k,nm) in enumerate(hk):
    r=hom.loc[k]
    ax.plot([x[i]]*2,[1,r.enrichment],"-",color=NAVY,lw=1.0,alpha=.55)
    ax.plot([x[i]],[r.enrichment],"o",ms=5.0,color=NAVY,ls="none")
ax.axhline(1,ls=(0,(5,4)),color="#555555",lw=1.0)
ax.set_xticks(x); ax.set_xticklabels([n for _,n in hk]); ax.set_xlim(-0.55,2.55)
ax.set_ylim(0.55,2.05); ax.set_ylabel("Observed / expected"); ax.set_xlabel("Host-element age class")
for s in ("top","right"): ax.spines[s].set_visible(False)
panel_letter(ax,"B",case="upper",dx=-0.30,dy=1.06)

ax=axes[2]
sz=18+(deg.n_insertions/deg.n_insertions.max())*95
ax.scatter(deg.pred_same,deg.observed_same,s=sz,facecolors="none",edgecolors=NAVY,linewidths=1.1,zorder=3)
lims=[0.60,1.04]; ax.plot(lims,lims,ls=(0,(5,4)),color="#555555",lw=1.0,zorder=1)
ax.set_xlim(*lims); ax.set_ylim(0.66,1.045)
ax.set_xlabel("Predicted same-orientation fraction"); ax.set_ylabel("Observed same-orientation fraction")
ax.set_xticks([0.65,0.75,0.85,0.95]); ax.set_yticks([0.7,0.8,0.9,1.0])
for s in ("top","right"): ax.spines[s].set_visible(False)
panel_letter(ax,"C",case="upper",dx=-0.30,dy=1.06)
fig.savefig(OUT+"figure3.png",bbox_inches="tight"); plt.close(fig)
print("fig2, fig3 ok")
