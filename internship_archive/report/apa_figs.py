
import pandas as pd, numpy as np, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
L="/Users/smyan/Desktop/Research/Research_Projects/L1/"
OUT=L+"report/figures/"
apply_figure_style(sizes=(9,8,7))
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Helvetica","Arial","DejaVu Sans"],
                     "axes.titlesize":9,"savefig.dpi":600,"figure.dpi":110})
NAVY="#17365D"; GREY="#8C8C8C"; LGREY="#C8C8C8"

def forest(ax, labels, est, lo, hi, filled, color=None, ms=5.2):
    y=np.arange(len(labels))[::-1]
    for i,(e,l,h,f) in enumerate(zip(est,lo,hi,filled)):
        c = (color[i] if color is not None else (NAVY if f else GREY))
        if not (np.isnan(l) or np.isnan(h)):
            ax.plot([l,h],[y[i]]*2,"-",color=c,lw=1.3,solid_capstyle="butt",zorder=2)
        ax.plot([e],[y[i]],marker="o" if f else "D",ms=ms if f else ms-0.7,
                mfc=c if f else "white",mec=c,mew=1.2,zorder=3,ls="none")
    ax.set_yticks(y); ax.set_yticklabels(labels); ax.set_ylim(-0.7,len(labels)-0.3)
    for s in ("top","right","left"): ax.spines[s].set_visible(False)
    ax.tick_params(axis="y",length=0)
    return y

# ---------------- FIGURE 1 ----------------
rows=[("Reference AluY in genes (n = 67,379)",0.4501,0.4463,0.4539,True),
      ("Alu, ten genomes (n = 2,428)",0.4477,0.4280,0.4675,True),
      ("Alu, HG03086 (n = 733)",0.4570,0.4205,0.4939,True),
      ("Alu, Gal-Mark et al. (2008)",0.45,np.nan,np.nan,False),
      ("Alu, Lavi & Carmel (2018)",0.451,np.nan,np.nan,False),
      ("All families, HG03086 (n = 856)",0.4463,0.4126,0.4803,True),
      ("LINE-1, ten genomes (n = 452)",0.3938,0.3498,0.4396,True),
      ("LINE-1, HG03086 (n = 102)",0.3824,0.2879,0.4839,True),
      ("LINE-1, Sela et al. (2007)",0.343,np.nan,np.nan,False),
      ("SVA, ten genomes (n = 121)",0.3719,0.2910,0.4607,True),
      ("SVA, HG03086 (n = 21)",0.3810,0.1811,0.6156,True)]
fig,ax=plt.subplots(figsize=(6.5,3.9))
lab,e,lo,hi,fl=zip(*rows)
forest(ax,lab,e,lo,hi,fl)
ax.axvline(0.50,ls=(0,(5,4)),color="#555555",lw=1.0,zorder=1)
ax.set_xlabel("Proportion of intronic insertions in the same orientation as the host gene")
ax.set_xlim(0.16,0.66); ax.set_xticks(np.arange(0.20,0.66,0.05))
ax.legend(handles=[Line2D([],[],marker="o",ls="none",mfc=NAVY,mec=NAVY,ms=5.2,label="Present study, 95% CI"),
                   Line2D([],[],marker="D",ls="none",mfc="white",mec=GREY,mew=1.2,ms=4.5,label="Published estimate")],
          loc="upper right",frameon=False,handletextpad=0.5,borderpad=0.2,
          bbox_to_anchor=(1.0,1.02))
fig.tight_layout(); fig.savefig(OUT+"figure1.png",bbox_inches="tight"); plt.close(fig)

