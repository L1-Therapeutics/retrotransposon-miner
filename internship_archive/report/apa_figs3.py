
import pandas as pd, numpy as np, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
L="/Users/smyan/Desktop/Research/Research_Projects/L1/"; OUT=L+"report/figures/"
apply_figure_style(sizes=(9,8,7))
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Helvetica","Arial","DejaVu Sans"],"savefig.dpi":600})
NAVY="#17365D"; GREY="#8C8C8C"
def bare(ax,sides=("top","right","left")):
    for s in sides: ax.spines[s].set_visible(False)
    ax.tick_params(axis="y",length=0)
def wilson(k,n,z=1.959964):
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d
    h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return c-h,c+h

# ---------------- FIGURE 4 ----------------
zy=pd.read_csv(L+"_job9/out/aluY/zones.tsv",sep="\t").set_index("zone")
dq=pd.read_csv(L+"_job9/out/aluY/by_divergence.tsv",sep="\t")
fig,axes=plt.subplots(1,2,figsize=(6.9,3.0),gridspec_kw=dict(width_ratios=[1.25,1],wspace=0.55))
ax=axes[0]
zk=[("within_150_upstream","Within 150 bp upstream\nof an exon"),
    ("within_150_downstream","Within 150 bp downstream"),
    ("overlaps_exon","Overlapping an exon"),
    ("intronic_far","Elsewhere in the intron")]
y=np.arange(len(zk))[::-1]
far=zy.loc["intronic_far","antisense_frac"]
ax.axvline(far,ls=(0,(5,4)),color="#555555",lw=1.0,zorder=1)
for i,(k,nm) in enumerate(zk):
    r=zy.loc[k]; lo,hi=wilson(r.antisense,r.n)
    ax.plot([lo,hi],[y[i]]*2,"-",color=NAVY,lw=1.3,zorder=2)
    ax.plot([r.antisense_frac],[y[i]],"o",ms=5.0,color=NAVY,ls="none",zorder=3)
ax.set_yticks(y); ax.set_yticklabels([n for _,n in zk]); ax.set_ylim(-0.7,len(zk)-0.35)
ax.set_xlabel("Proportion of AluY copies antisense to the host gene")
ax.set_xlim(0.448,0.572); ax.set_xticks([0.46,0.48,0.50,0.52,0.54,0.56])
ax.annotate("far-intronic baseline",xy=(far-0.0015,len(zk)-0.60),ha="right",va="bottom",
            fontsize=7,color="#555555",annotation_clip=False)
bare(ax); panel_letter(ax,"A",case="upper")

ax=axes[1]; x=np.arange(len(dq))
for i,r in dq.iterrows():
    lo,hi=wilson(round(r.antisense_frac_near*r.n_near),int(r.n_near))
    ax.plot([x[i]]*2,[lo,hi],"-",color=NAVY,lw=1.3,zorder=2)
    ax.plot([x[i]],[r.antisense_frac_near],"o",ms=5.0,color=NAVY,ls="none",zorder=3)
    ax.plot([x[i]],[r.antisense_frac_far],marker="_",ms=9,mew=1.6,color=GREY,ls="none",zorder=3)
ax.set_xticks(x)
ax.set_xticklabels([f"Q{i+1}\n(n = {int(r.n_near)})" for i,r in dq.iterrows()],linespacing=1.3)
ax.set_xlim(-0.6,len(dq)-0.4); ax.set_ylabel("Proportion antisense")
ax.set_xlabel("Divergence quartile, youngest to oldest")
for s in ("top","right"): ax.spines[s].set_visible(False)
ax.annotate("far intronic",xy=(len(dq)-1,dq.antisense_frac_far.iloc[-1]),
            xytext=(6,5),textcoords="offset points",ha="right",va="bottom",fontsize=7,color=GREY)
panel_letter(ax,"B",case="upper",dx=-0.26,dy=1.06)
fig.savefig(OUT+"figure4.png",bbox_inches="tight"); plt.close(fig)

# ---------------- FIGURE 5 ----------------
sw=pd.read_csv(L+"_job9/subfamily_summary.tsv",sep="\t").set_index("rep")
rows=["AluY","AluSg","AluSx","AluSz","AluJb","AluJo"]
fig,ax=plt.subplots(figsize=(5.6,3.0))
y=np.arange(len(rows))[::-1]
for i,r in enumerate(rows):
    x=sw.loc[r]
    k_n=round(x.anti_up*x.n_up); lo_n,hi_n=wilson(k_n,int(x.n_up))
    orlo=(lo_n/(1-lo_n))/(x.anti_far/(1-x.anti_far)); orhi=(hi_n/(1-hi_n))/(x.anti_far/(1-x.anti_far))
    ax.plot([orlo,orhi],[y[i]]*2,"-",color=NAVY,lw=1.3,zorder=2)
    ax.plot([x.OR],[y[i]],"o",ms=5.0,color=NAVY,ls="none",zorder=3)
ax.set_yticks(y)
ax.set_yticklabels([f"{r}   ({int(sw.loc[r,'median_div'])} milliDiv, n = {int(sw.loc[r,'n_up']):,})" for r in rows])
ax.set_ylim(-0.7,len(rows)-0.3)
ax.axvline(1.0,ls=(0,(5,4)),color="#555555",lw=1.0,zorder=1)
ax.set_xscale("log"); ax.set_xlim(0.56,1.30)
ax.set_xticks([0.6,0.7,0.8,0.9,1.0,1.2]); ax.set_xticklabels(["0.6","0.7","0.8","0.9","1.0","1.2"])
ax.set_xlabel("Odds ratio for antisense orientation within 150 bp upstream of an exon")
ax.annotate("no depletion",xy=(1.02,len(rows)-0.52),ha="left",va="bottom",fontsize=7,
            color="#555555",annotation_clip=False)
bare(ax)
fig.savefig(OUT+"figure5.png",bbox_inches="tight"); plt.close(fig)
print("fig4, fig5 written")
