"""Figures 4, 6 and 7: origin asymmetry, topology ablations and capacity control."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[2]; HPC=ROOT/"artifacts"/"hpc_results"; OUT=ROOT/"manuscript_rebuild"/"figures"; OUT.mkdir(parents=True,exist_ok=True)
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],"svg.fonttype":"none","pdf.fonttype":42,"font.size":7.3,"axes.titlesize":8.2,"axes.labelsize":7.1,"xtick.labelsize":6.5,"ytick.labelsize":6.5,"axes.linewidth":0.7})
INK="#1F2933"; MUTED="#5B6770"; GRID="#D9E0E5"; NAVY="#28577D"; TEAL="#3FA7A3"; ORANGE="#D77A32"; PURPLE="#7A5A9E"; RED="#B84C4C"; GOLD="#C99724"

def style(ax):
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False); ax.grid(axis="y",color=GRID,lw=.5,alpha=.8); ax.set_axisbelow(True)
def lab(ax,s): ax.text(-.14,1.05,s,transform=ax.transAxes,fontsize=9.3,fontweight="bold",color=INK)
def save(fig,stem):
    for ext,dpi in [("png",600),("tiff",600),("pdf",None),("svg",None)]:
        kw={"dpi":dpi} if dpi else {}; fig.savefig(OUT/f"{stem}.{ext}",bbox_inches="tight",facecolor="white",**kw)
    plt.close(fig)

def figure4():
    ext=pd.read_csv(HPC/"origin_supplement"/"prism_origin_20261110_092718"/"prism_origin_metrics.csv")
    gat=pd.read_csv(HPC/"closure_final"/"origin_paired_20261110_131203"/"origin_paired_metrics.csv")
    boot=pd.read_csv(HPC/"closure_final"/"origin_paired_20261110_131203"/"origin_paired_bootstrap.csv")
    # External tabular data: aggregate source validation vs external target.
    ext2=ext.groupby(["direction","model","evaluation"],as_index=False).roc_auc.agg(["mean","std"]).reset_index()
    fig,axes=plt.subplots(2,2,figsize=(8.8,6.5),gridspec_kw={"width_ratios":[1.12,1]})
    fig.suptitle("PRISM origin labels define a directional transfer estimand",fontsize=10.2,fontweight="bold",y=.99)
    # a paired source/external points
    ax=axes[0,0]; style(ax); lab(ax,"a")
    colors={"logistic":NAVY,"random_forest":TEAL}
    for yi,d in enumerate(["PRISM_discovery_to_validation","PRISM_validation_to_discovery"]):
        for mi,m in enumerate(["logistic","random_forest"]):
            fr=ext[(ext.direction==d)&(ext.model==m)]
            src=fr[fr.evaluation=="source_validation"].roc_auc.to_numpy(); tar=fr[fr.evaluation=="external_target"].roc_auc.to_numpy()
            x=yi*3+mi
            for a,b in zip(src,tar): ax.plot([x-.16,x+.16],[a,b],color=GRID,lw=.7)
            ax.scatter(np.full(len(src),x-.16),src,s=16,color=colors[m],alpha=.55)
            ax.scatter(np.full(len(tar),x+.16),tar,s=16,color=colors[m],alpha=.9)
    ax.set_xticks([.5,3.5],["Discovery → validation","Validation → discovery"]); ax.set_ylabel("AUROC across origin splits"); ax.set_title("External performance falls relative to source validation",loc="left")
    ax.legend([Line2D([0],[0],marker="o",color=colors[m],lw=0,markersize=5) for m in colors],["L2 logistic","Random forest"],frameon=False,fontsize=6,loc="lower left")
    # b origin GAT vs references
    ax=axes[0,1]; style(ax); lab(ax,"b")
    means=gat.groupby(["direction","model"],as_index=False).roc_auc.agg(["mean","std"]).reset_index()
    model_order=["no_prior_gat","elastic_net","random_forest"]; cols={"no_prior_gat":NAVY,"elastic_net":ORANGE,"random_forest":TEAL}; labels={"no_prior_gat":"No-Prior GAT","elastic_net":"Elastic net","random_forest":"Random forest"}
    for di,d in enumerate(["PRISM_discovery_to_validation","PRISM_validation_to_discovery"]):
        for mi,m in enumerate(model_order):
            r=means[(means.direction==d)&(means.model==m)].iloc[0]; x=di*4+mi
            ax.errorbar(x,r["mean"],yerr=r["std"],fmt="o",ms=5,color=cols[m],capsize=2,lw=1.1)
    ax.axhline(.5,color=MUTED,ls="--",lw=.8); ax.set_xticks([1,5],["Discovery → validation","Validation → discovery"]); ax.set_ylabel("AUROC ± split s.d."); ax.set_title("No-Prior GAT is not uniformly superior",loc="left")
    ax.legend([Line2D([0],[0],marker="o",color=cols[m],lw=0,markersize=5) for m in model_order],[labels[m] for m in model_order],frameon=False,fontsize=5.8,loc="lower left")
    # c paired bootstrap forest from origin artifact
    ax=axes[1,0]; style(ax); lab(ax,"c")
    # use per-split delta CIs as a distribution, rather than treating pooled CI as n=1
    boot["comparison_label"]=boot.comparison.str.replace("gat_vs_","GAT − ",regex=False).str.replace("elastic_net","elastic net").str.replace("random_forest","random forest")
    rows=[]
    for d in boot.direction.unique():
        for c in boot.comparison.unique():
            fr=boot[(boot.direction==d)&(boot.comparison==c)]
            rows.append((d,c,fr.delta_auc.mean(),fr.delta_auc.std(ddof=1)/np.sqrt(len(fr))*1.96))
    for yi,(d,c,mu,se) in enumerate(rows):
        col=ORANGE if "elastic" in c else TEAL; ax.errorbar(mu,yi,xerr=se,fmt="o",color=col,capsize=2,lw=1.2,ms=4.5)
    ax.axvline(0,color=MUTED,ls="--",lw=.8); ax.set_yticks(range(len(rows)),[f"{('D→V' if d.endswith('validation') else 'V→D')}\n{c.replace('gat_vs_','GAT − ').replace('elastic_net','elastic net').replace('random_forest','random forest')}" for d,c,_,_ in rows]); ax.set_xlabel("Mean split-specific paired ΔAUROC"); ax.set_title("Origin transfer remains direction dependent",loc="left"); ax.grid(axis="x",color=GRID,lw=.5)
    # d PR-AUC loss on transfer
    ax=axes[1,1]; style(ax); lab(ax,"d")
    for i,d in enumerate(["PRISM_discovery_to_validation","PRISM_validation_to_discovery"]):
        for m in ["logistic","random_forest"]:
            fr=ext[(ext.direction==d)&(ext.model==m)]; src=fr[fr.evaluation=="source_validation"].average_precision.mean(); tar=fr[fr.evaluation=="external_target"].average_precision.mean(); x=i*3+(0 if m=="logistic" else 1)
            ax.plot([x-.15,x+.15],[src,tar],color=colors[m],lw=1.4); ax.scatter([x-.15,x+.15],[src,tar],color=colors[m],s=20,zorder=3)
    ax.set_xticks([.5,3.5],["Discovery → validation","Validation → discovery"]); ax.set_ylabel("Average precision"); ax.set_title("Precision also shifts with origin",loc="left")
    fig.tight_layout(rect=[.02,.02,.99,.95]); save(fig,"Figure4_origin_asymmetry")

def figure6():
    # Keep this legacy entry point safe: the submission figure is generated by
    # the audited four-panel topology script, which excludes the capacity panel
    # and does not treat fixed-k degree correlation as ecological evidence.
    import runpy
    runpy.run_path(str(Path(__file__).with_name('plot_topology_final.py')), run_name='__main__')
    return
    paths={}
    for d in ["identity","no_edge","pearson","shuffled","spearman"]:
        p=list((HPC/"no_prior_gat_fixed").glob(f"no_prior_gat_{d}_*/no_prior_gat_summary.csv")); paths[d]=p[0] if p else None
    rows=[]
    for t,p in paths.items():
        if p:
            df=pd.read_csv(p,header=[0,1]); # two-row header from summary
            # normalize repeated metric columns
            for _,r in df.iterrows(): rows.append({"topology":t,"direction":r.iloc[0],"auc_mean":float(r.iloc[3]),"auc_sd":float(r.iloc[4])})
    topo=pd.DataFrame(rows)
    cap=pd.read_csv(HPC/"capacity_control"/"capacity_20261110_103951"/"capacity_control_metrics.csv")
    gat=pd.read_csv(HPC/"final_supplement"/"no_prior_gat_spearman_20261110_082627"/"no_prior_gat_metrics.csv")
    fig,axes=plt.subplots(2,2,figsize=(8.8,6.4))
    fig.suptitle("Topology ablations localize the relational contribution",fontsize=10.2,fontweight="bold",y=.99)
    # a topology heatmap
    ax=axes[0,0]; lab(ax,"a"); ax.set_title("Mean AUROC across topology controls",loc="left")
    order=["no_edge","identity","shuffled","pearson","spearman"]; dorder=["HMP2_to_PRISM","PRISM_to_HMP2"]
    mat=np.array([[topo[(topo.topology==t)&(topo.direction==d)].auc_mean.iloc[0] for d in dorder] for t in order])
    im=ax.imshow(mat,cmap=mpl.colors.LinearSegmentedColormap.from_list("topo",["#F5F7F8",NAVY,TEAL]),vmin=.45,vmax=.88,aspect="auto")
    ax.set_xticks(range(2),["HMP2 → PRISM","PRISM → HMP2"]); ax.set_yticks(range(len(order)),[x.replace("_"," ").title() for x in order])
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]): ax.text(j,i,f"{mat[i,j]:.3f}",ha="center",va="center",fontsize=6.8,color="white" if mat[i,j]>.68 else INK,fontweight="bold")
    ax.set_xlabel("Transfer direction"); ax.set_ylabel("Topology")
    # b topology deltas relative no-edge
    ax=axes[0,1]; style(ax); lab(ax,"b"); ax.set_title("Spearman topology is the only strong positive contrast",loc="left")
    for d,yi in zip(dorder,[1,0]):
        base=topo[(topo.topology=="no_edge")&(topo.direction==d)].auc_mean.iloc[0]
        for x,t in enumerate(["identity","shuffled","pearson","spearman"]):
            v=topo[(topo.topology==t)&(topo.direction==d)].auc_mean.iloc[0]-base; ax.scatter(v,yi+(x-1.5)*.06,color=[MUTED,ORANGE,TEAL,NAVY][x],s=24)
    ax.axvline(0,color=MUTED,ls="--",lw=.8); ax.set_yticks([0,1],["PRISM → HMP2","HMP2 → PRISM"]); ax.set_xlabel("AUROC change from no-edge"); ax.set_xlim(-.08,.35); ax.grid(axis="x",color=GRID,lw=.5)
    # c split-wise GAT vs capacity MLP
    ax=axes[1,0]; style(ax); lab(ax,"c"); ax.set_title("Capacity alone does not reproduce the favorable direction",loc="left")
    for di,d in enumerate(dorder):
        g=gat[(gat.direction==d)&(gat.model=="no_prior_gat")].sort_values("split").roc_auc.to_numpy(); m=cap[cap.direction==d].sort_values("split").roc_auc.to_numpy(); x=np.arange(len(g))+di*12
        ax.plot(x,g,"o-",color=NAVY,lw=1,ms=3,label="No-Prior GAT" if di==0 else None); ax.plot(x,m,"o-",color=ORANGE,lw=1,ms=3,label="capacity-matched MLP" if di==0 else None)
    ax.axhline(.5,color=MUTED,ls="--",lw=.8); ax.set_xticks([4.5,16.5],["HMP2 → PRISM","PRISM → HMP2"]); ax.set_ylabel("AUROC per source split"); ax.legend(frameon=False,fontsize=5.8,loc="lower left")
    # d edge-set diagnostics from frozen Results values
    ax=axes[1,1]; style(ax); lab(ax,"d"); ax.set_title("Topology overlap is moderate; degree is construction-fixed",loc="left")
    labels=["Pearson vs Spearman\nHMP2 → PRISM","Pearson vs Spearman\nPRISM → HMP2","Spearman vs rewired\nHMP2 → PRISM","Spearman vs rewired\nPRISM → HMP2"]
    vals=[.325,.415,.054,.043]; cols=[PURPLE,PURPLE,ORANGE,ORANGE]; y=np.arange(len(vals))[::-1]
    ax.hlines(y,0,vals,color=cols,lw=4,alpha=.55); ax.scatter(vals,y,color=cols,s=32,zorder=3); ax.set_yticks(y,labels); ax.set_xlim(0,.5); ax.set_xlabel("Edge-set Jaccard similarity"); ax.grid(axis="x",color=GRID,lw=.5); ax.text(.48,.25,"degree corr. = 1.0\n(fixed top-k diagnostic)",ha="right",va="center",fontsize=6.2,color=MUTED)
    fig.tight_layout(rect=[.02,.02,.99,.95]); save(fig,"Figure6_topology_controls")

def figure7():
    cap=pd.read_csv(HPC/"capacity_control"/"capacity_20261110_103951"/"capacity_control_metrics.csv")
    gat=pd.read_csv(HPC/"final_supplement"/"no_prior_gat_spearman_20261110_082627"/"no_prior_gat_metrics.csv")
    fig,axes=plt.subplots(2,2,figsize=(8.5,6.1)); fig.suptitle("Capacity-matched control separates graph message passing from generic nonlinear capacity",fontsize=10.0,fontweight="bold",y=.99)
    dirs=["HMP2_to_PRISM","PRISM_to_HMP2"]; labdir={dirs[0]:"HMP2 → PRISM",dirs[1]:"PRISM → HMP2"}
    # a paired split scatter
    ax=axes[0,0]; style(ax); lab(ax,"a")
    for di,d in enumerate(dirs):
        g=gat[(gat.direction==d)&(gat.model=="no_prior_gat")].sort_values("split").roc_auc.to_numpy(); m=cap[cap.direction==d].sort_values("split").roc_auc.to_numpy(); x=np.arange(len(g))+di*12
        for xx,gg,mm in zip(x,g,m): ax.plot([xx,xx],[mm,gg],color=GRID,lw=.8); ax.scatter(xx,mm,color=ORANGE,s=18); ax.scatter(xx,gg,color=NAVY,s=18)
    ax.axhline(.5,color=MUTED,ls="--",lw=.8); ax.set_xticks([4.5,16.5],[labdir[d] for d in dirs]); ax.set_ylabel("AUROC per split"); ax.set_title("Relational gain is not symmetric",loc="left")
    ax.legend([Line2D([0],[0],marker="o",color=ORANGE,lw=0,markersize=5),Line2D([0],[0],marker="o",color=NAVY,lw=0,markersize=5)],["Capacity-matched MLP","No-Prior GAT"],frameon=False,fontsize=5.8,loc="lower left")
    # b metric means/CI-ish split sd
    ax=axes[0,1]; style(ax); lab(ax,"b"); ax.set_title("Mean AUROC by direction",loc="left")
    for di,d in enumerate(dirs):
        for mi,(name,df,col) in enumerate([("MLP",cap,ORANGE),("GAT",gat,NAVY)]):
            fr=df[df.direction==d].roc_auc; x=di*3+mi; ax.errorbar(x,fr.mean(),yerr=fr.std(ddof=1),fmt="o",color=col,ms=5,capsize=2,lw=1.1)
    ax.axhline(.5,color=MUTED,ls="--",lw=.8); ax.set_xticks([.5,3.5],[labdir[d] for d in dirs]); ax.set_ylabel("AUROC ± split s.d.")
    # c AP scatter
    ax=axes[1,0]; style(ax); lab(ax,"c"); ax.set_title("Precision follows the same directional pattern",loc="left")
    for d in dirs:
        g=gat[gat.direction==d]; m=cap[cap.direction==d]; ax.scatter(g.roc_auc,g.average_precision,color=NAVY,s=24,label="GAT" if d==dirs[0] else None); ax.scatter(m.roc_auc,m.average_precision,color=ORANGE,s=24,label="MLP" if d==dirs[0] else None)
    ax.set_xlabel("AUROC"); ax.set_ylabel("Average precision"); ax.legend(frameon=False,fontsize=6,loc="lower right")
    # d Brier paired differences
    ax=axes[1,1]; style(ax); lab(ax,"d"); ax.set_title("Calibration-related Brier differences are direction dependent",loc="left")
    for di,d in enumerate(dirs):
        g=gat[(gat.direction==d)&(gat.model=="no_prior_gat")].sort_values("split").brier.to_numpy(); m=cap[cap.direction==d].sort_values("split").brier.to_numpy(); x=np.arange(len(g))+di*12
        ax.scatter(x,g-m,color=PURPLE,s=22)
    ax.axhline(0,color=MUTED,ls="--",lw=.8); ax.set_xticks([4.5,16.5],[labdir[d] for d in dirs]); ax.set_ylabel("Brier(GAT) − Brier(MLP)")
    fig.tight_layout(rect=[.02,.02,.99,.95]); save(fig,"Figure7_capacity_control")

if __name__=="__main__":
    figure4(); figure6(); figure7()
