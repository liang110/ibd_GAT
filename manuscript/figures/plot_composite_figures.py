"""Composite figures for the reconstructed manuscript.

The panels are deliberately mixed-modality so that each figure communicates an
evidence chain rather than a single bar chart. All data are read from the
archived HPC result tables.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
HPC = ROOT / "artifacts" / "hpc_results"
OUT = ROOT / "manuscript_rebuild" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 7.4, "axes.titlesize": 8.2, "axes.labelsize": 7.1,
    "axes.linewidth": 0.7, "xtick.labelsize": 6.6, "ytick.labelsize": 6.6,
})

INK, MUTED, GRID = "#1F2933", "#5B6770", "#D9E0E5"
NAVY, TEAL, ORANGE, PURPLE, RED, GOLD = "#28577D", "#3FA7A3", "#D77A32", "#7A5A9E", "#B84C4C", "#C99724"
PALETTE = {"No-Prior Spearman GAT": NAVY, "Elastic net": ORANGE, "Random forest": TEAL}


def style_axes(ax):
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color=GRID, lw=0.5, alpha=0.75); ax.set_axisbelow(True)


def panel_label(ax, letter):
    ax.text(-0.14, 1.05, letter, transform=ax.transAxes, fontsize=9.4, fontweight="bold", color=INK)


def save(fig, stem):
    fig.savefig(OUT / f"{stem}.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.tiff", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def figure2():
    scale = pd.read_csv(HPC / "gate0_20261110_000415" / "input_scale_audit.csv")
    fig, axes = plt.subplots(2, 2, figsize=(8.7, 6.1), gridspec_kw={"height_ratios":[1,1.08]})
    fig.suptitle("Input contract fixes the cross-cohort comparison before model fitting", fontsize=10.2, fontweight="bold", y=0.99)
    # a: scale range as interval + point, with log x axis.
    ax = axes[0,0]; style_axes(ax); panel_label(ax,"a")
    for yi, (_, r) in enumerate(scale.iterrows()):
        color = NAVY if r.cohort == "HMP2" else TEAL
        ax.plot([max(r.sum_min, 1e-3), r.sum_max], [yi, yi], color=color, lw=3, alpha=0.55, solid_capstyle="round")
        ax.scatter(r.sum_median, yi, s=54, color=color, edgecolor="white", linewidth=0.8, zorder=3)
        ax.text(r.sum_median*1.35 if r.sum_median > 1 else 1.15, yi+0.12, f"median {r.sum_median:g}", fontsize=6.4, color=color)
    ax.set_xscale("log"); ax.set_yticks([0,1], ["HMP2", "PRISM"]); ax.set_xlabel("Library sum (log scale; min–max interval)")
    ax.set_title("Raw matrices were on different scales", loc="left")
    ax.grid(axis="x", color=GRID, lw=0.5); ax.grid(axis="y", visible=False)
    # b: feature contract as a flow / set relationship.
    ax = axes[0,1]; style_axes(ax); panel_label(ax,"b"); ax.axis("off")
    ax.set_title("Feature contract after conservative canonicalization", loc="left")
    # scaled horizontal strips
    yvals = [0.75, 0.49, 0.23]
    entries = [("HMP2 raw microbial", 578, NAVY, 0.12, 0.76, 6.25),
               ("shared microbiome", 126, ORANGE, 0.34, 0.36, 6.15),
               ("PRISM raw microbial", 209, TEAL, 0.54, 0.38, 6.2)]
    for y,(lab,n,col,left,width,fs) in zip(yvals, entries):
        ax.add_patch(plt.Rectangle((left,y-0.07), width, 0.13, facecolor=col, alpha=0.16, edgecolor=col, lw=1.1))
        ax.text(left+0.012, y, lab, va="center", ha="left", color=INK, fontsize=fs, fontweight="bold")
        ax.text(left+width-0.012, y, f"{n}", va="center", ha="right", color=col, fontsize=7.2, fontweight="bold")
    ax.annotate("", xy=(0.37,0.64), xytext=(0.78,0.64), arrowprops=dict(arrowstyle="-|>",color=ORANGE,lw=1.0))
    ax.text(0.57,0.68,"shared intersection",ha="center",fontsize=6.4,color=ORANGE)
    ax.annotate("", xy=(0.62,0.38), xytext=(0.37,0.38), arrowprops=dict(arrowstyle="-|>",color=ORANGE,lw=1.0))
    ax.text(0.50,0.42,"same 126 features",ha="center",fontsize=6.4,color=ORANGE)
    ax.text(0.14,0.05,"Exact-name metabolite intersection = 0",fontsize=7.2,color=RED,fontweight="bold")
    # c: source/target preprocessing boundary as a contract matrix.
    ax = axes[1,0]; style_axes(ax); panel_label(ax,"c")
    ax.set_title("Source-only fitting defines the evaluation boundary", loc="left")
    rows = ["Canonicalization", "CLR / scaling", "Feature selection", "Denylist", "Split assignment"]
    vals = np.ones((len(rows),2))
    im = ax.imshow(vals, cmap=mpl.colors.ListedColormap(["#EAF7F5", "#EAF2F8"]), vmin=0, vmax=1, aspect="auto")
    ax.set_xticks([0,1], ["Source training", "Target"]); ax.set_yticks(range(len(rows)), rows)
    for i in range(len(rows)):
        ax.text(0,i,"fit",ha="center",va="center",fontsize=6.4,color=TEAL,fontweight="bold")
        ax.text(1,i,"transform only",ha="center",va="center",fontsize=6.4,color=NAVY)
    ax.set_xlabel("Contract state"); ax.grid(False)
    # d: summary badge / immutable artifacts.
    ax = axes[1,1]; style_axes(ax); panel_label(ax,"d"); ax.set_title("Immutable artifacts carried into every run", loc="left")
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis("off")
    artifacts = [("feature-order hash", "f1a8…23b4", NAVY), ("canonicalizer hash", "5365…31e5", TEAL), ("HMP2 zero-mass exclusions", "3", ORANGE), ("PRISM zero-mass exclusions", "0", PURPLE)]
    for i,(lab,val,col) in enumerate(artifacts):
        y = 0.78-i*0.19
        ax.add_patch(plt.Circle((0.12,y),0.035,facecolor=col,alpha=0.18,edgecolor=col,lw=1.1))
        ax.text(0.12,y,"OK",ha="center",va="center",fontsize=5.7,color=col,fontweight="bold")
        ax.text(0.21,y+0.025,lab,ha="left",va="center",fontsize=6.8,color=MUTED)
        ax.text(0.21,y-0.035,val,ha="left",va="center",fontsize=7.4,color=INK,fontweight="bold")
    fig.tight_layout(rect=[0.02,0.02,0.99,0.95]); save(fig,"Figure2_input_contract")


def figure3():
    leak = pd.read_csv(HPC / "repeated_leakage_20261110_001007" / "repeated_leakage_metrics.csv")
    syn = pd.read_csv(HPC / "synthetic_supplement" / "synthetic_controls_v2_20261110_091610" / "synthetic_control_metrics.csv")
    syn["scenario"] = syn["scenario"].fillna("null").replace({
        "signal_reversing_batch": "reversing_batch",
        "signal_topology_shift": "topology_shift",
    })
    stable = pd.read_csv(HPC / "synthetic_supplement" / "stable_synthetic_v2_20261110_091610" / "domain_stable_sweep_metrics.csv")
    fig, axes = plt.subplots(2,2,figsize=(8.8,6.6))
    fig.suptitle("Participant structure and known-truth controls expose the dominant failure modes", fontsize=10.2, fontweight="bold", y=0.99)
    # a paired leakage AUROC
    ax=axes[0,0]; style_axes(ax); panel_label(ax,"a")
    ax.plot(np.arange(len(leak)), leak.random_split_roc_auc, "o-", color=RED, lw=1, ms=3.5, label="naive sample split")
    ax.plot(np.arange(len(leak)), leak.group_split_roc_auc, "o-", color=NAVY, lw=1, ms=3.5, label="participant grouped")
    ax.axhline(0.5,color=MUTED,ls="--",lw=0.8); ax.set_xlabel("Repeated leakage control seed"); ax.set_ylabel("AUROC")
    ax.set_title("Apparent signal collapses after grouping",loc="left"); ax.legend(frameon=False,fontsize=6.0,loc="lower left")
    # b paired participant overlap
    ax=axes[0,1]; style_axes(ax); panel_label(ax,"b")
    for i,r in leak.iterrows():
        ax.plot([0,1],[r.random_split_subject_overlap,r.group_split_subject_overlap],color=GRID,lw=0.8,zorder=1)
    ax.scatter(np.zeros(len(leak)), leak.random_split_subject_overlap, s=22, color=RED, zorder=2)
    ax.scatter(np.ones(len(leak)), leak.group_split_subject_overlap, s=22, color=NAVY, zorder=2)
    ax.set_xticks([0,1],["naive", "grouped"]); ax.set_ylabel("Participants shared across split")
    ax.set_title("The estimand changes with participant overlap",loc="left")
    # c synthetic controls: dot + mean intervals for scenarios
    ax=axes[1,0]; style_axes(ax); panel_label(ax,"c")
    order=["null","batch_only","signal_only","reversing_batch","topology_shift"]
    labels={"null":"Null","batch_only":"Batch only","signal_only":"Signal only","reversing_batch":"Reversing batch","topology_shift":"Topology shift"}
    colors={"null":MUTED,"batch_only":ORANGE,"signal_only":NAVY,"reversing_batch":TEAL,"topology_shift":PURPLE}
    for yi,sc in enumerate(order):
        fr=syn[syn.scenario==sc]
        if fr.empty: continue
        jitter=np.linspace(-0.13,0.13,len(fr))
        ax.scatter(fr.target_roc_auc, np.full(len(fr),yi)+jitter, s=18, color=colors[sc], alpha=0.55, edgecolor="white",lw=0.3)
        ax.errorbar(fr.target_roc_auc.mean(), yi, xerr=fr.target_roc_auc.std(ddof=1)/np.sqrt(len(fr))*1.96, fmt="o", color=colors[sc], ms=5, capsize=2, lw=1.2)
    ax.axvline(0.5,color=MUTED,ls="--",lw=0.8); ax.set_yticks(range(len(order)),[labels[x] for x in order]); ax.set_xlabel("Target AUROC across seeds")
    ax.set_title("Known-truth scenarios separate signal from batch",loc="left"); ax.set_xlim(0.35,0.9)
    # d stable selector: domain stability as a scatter + pass coding
    ax=axes[1,1]; style_axes(ax); panel_label(ax,"d")
    passed=stable.all_acceptance_pass.astype(bool)
    ax.scatter(stable.n_stable_candidates[~passed], stable.target_roc_auc[~passed], s=26, color=RED, alpha=0.75, label="criteria not met")
    ax.scatter(stable.n_stable_candidates[passed], stable.target_roc_auc[passed], s=26, color=TEAL, alpha=0.8, label="all criteria met")
    ax.axhline(0.7,color=MUTED,ls="--",lw=0.8); ax.set_xlabel("Stable candidates retained"); ax.set_ylabel("Target AUROC")
    ax.set_title("Domain-stable selection improves the operating boundary",loc="left"); ax.legend(frameon=False,fontsize=6.0,loc="lower right")
    fig.tight_layout(rect=[0.02,0.02,0.99,0.95]); save(fig,"Figure3_controls")


def figure5():
    metrics = pd.read_csv(HPC / "final_supplement" / "no_prior_gat_spearman_20261110_082627" / "no_prior_gat_metrics.csv")
    boot = pd.read_csv(HPC / "final_supplement" / "paired_bootstrap_20261110_082627" / "paired_bootstrap_summary.csv")
    models = {"no_prior_gat":"No-Prior Spearman GAT","elastic_net":"Elastic net","random_forest":"Random forest"}
    metrics["model_label"] = metrics.model.map(models)
    dirs=["HMP2_to_PRISM","PRISM_to_HMP2"]; dirlab={dirs[0]:"HMP2 → PRISM", dirs[1]:"PRISM → HMP2"}
    fig, axes = plt.subplots(2,2,figsize=(8.8,6.5), gridspec_kw={"width_ratios":[1.14,1]})
    fig.suptitle("Directional transfer separates a conditional relational gain from a universal GAT claim", fontsize=10.2, fontweight="bold", y=0.99)
    # a per-split AUROC distributions
    ax=axes[0,0]; style_axes(ax); panel_label(ax,"a")
    positions=[]; labels=[]
    for di,d in enumerate(dirs):
        for mi,m in enumerate(models):
            fr=metrics[(metrics.direction==d)&(metrics.model==m)]
            x=di*4+mi
            positions.append(x); labels.append(models[m])
            jitter=np.linspace(-0.08,0.08,len(fr))
            ax.scatter(np.full(len(fr),x)+jitter,fr.roc_auc,s=18,color=PALETTE[models[m]],alpha=0.65,edgecolor="white",lw=0.35,zorder=2)
            ax.errorbar(x,fr.roc_auc.mean(),yerr=fr.roc_auc.std(ddof=1)/np.sqrt(len(fr))*1.96,fmt="o",ms=5,color=PALETTE[models[m]],capsize=2,lw=1.1,zorder=3)
    ax.axhline(0.5,color=MUTED,ls="--",lw=0.8); ax.set_xticks([1,5],["HMP2 → PRISM","PRISM → HMP2"]); ax.set_ylabel("AUROC across source splits"); ax.set_title("Split-level performance",loc="left")
    handles=[Line2D([0],[0],marker="o",color="none",markerfacecolor=PALETTE[v],label=v,markersize=5) for v in models.values()]
    ax.legend(handles=handles,frameon=False,fontsize=5.9,loc="lower left")
    # b paired bootstrap forest
    ax=axes[0,1]; style_axes(ax); panel_label(ax,"b")
    boot["comparison_label"] = boot.comparison.str.replace("no_prior_gat_vs_","GAT − ",regex=False).str.replace("elastic_net","elastic net").str.replace("random_forest","random forest")
    y=np.arange(len(boot))[::-1]
    for yi,(_,r) in zip(y,boot.iterrows()):
        color=ORANGE if "elastic" in r.comparison else TEAL
        ax.plot([r.delta_auc_ci_low,r.delta_auc_ci_high],[yi,yi],color=color,lw=1.5)
        ax.scatter(r.delta_auc,yi,color=color,s=28,zorder=3)
    ax.axvline(0,color=MUTED,ls="--",lw=0.8); ax.set_yticks(y,[f"{dirlab[d]}\n{c}" for d,c in zip(boot.direction,boot.comparison_label)])
    ax.set_xlabel("ΔAUROC (GAT − reference; split-summary interval)"); ax.set_title("Uncertainty is direction-specific",loc="left"); ax.grid(axis="x",color=GRID,lw=0.5)
    # c AUROC vs PR-AUC, per split
    ax=axes[1,0]; style_axes(ax); panel_label(ax,"c")
    for m in models:
        fr=metrics[metrics.model==m]
        ax.scatter(fr.roc_auc, fr.average_precision, s=22, color=PALETTE[models[m]], alpha=0.72, label=models[m])
    ax.plot([0.45,0.95],[0.45,0.95],color=GRID,ls="--",lw=0.8); ax.set_xlabel("AUROC"); ax.set_ylabel("Average precision"); ax.set_title("Ranking and precision agree on the favorable direction",loc="left"); ax.legend(frameon=False,fontsize=6.0,loc="lower right")
    # d compact direction × model heatmap of mean AUROC
    ax=axes[1,1]; panel_label(ax,"d"); ax.set_title("Mean AUROC by direction and model",loc="left")
    means=metrics.pivot_table(index="model_label",columns="direction",values="roc_auc",aggfunc="mean").reindex(list(models.values()))
    im=ax.imshow(means.values,cmap=mpl.colors.LinearSegmentedColormap.from_list("bluegreen",["#EAF2F8",NAVY,TEAL]),vmin=0.5,vmax=0.9,aspect="auto")
    ax.set_xticks(range(2),[dirlab[d] for d in dirs]); ax.set_yticks(range(3),list(models.values()))
    for i in range(3):
        for j in range(2): ax.text(j,i,f"{means.iloc[i,j]:.3f}",ha="center",va="center",fontsize=6.8,color="white" if means.iloc[i,j]>0.68 else INK,fontweight="bold")
    ax.set_xlabel("Transfer direction"); ax.set_ylabel("")
    fig.tight_layout(rect=[0.02,0.02,0.99,0.95]); save(fig,"Figure5_transfer")


if __name__ == "__main__":
    figure2(); figure3(); figure5()
