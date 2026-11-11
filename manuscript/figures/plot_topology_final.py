from pathlib import Path
import numpy as np,pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
ROOT=Path(__file__).resolve().parents[2]; HPC=ROOT/'artifacts'/'hpc_results'; OUT=ROOT/'manuscript_rebuild'/'figures'
mpl.rcParams.update({'font.family':'DejaVu Sans','font.size':7.2,'axes.titlesize':8.2,'axes.labelsize':7.1,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'svg.fonttype':'none','pdf.fonttype':42,'xtick.labelsize':6.5,'ytick.labelsize':6.5})
INK='#1F2933';MUTED='#5B6770';GRID='#D9E0E5';NAVY='#28577D';TEAL='#3FA7A3';ORANGE='#D77A32';PURPLE='#7A5A9E';RED='#B84C4C';

def style(ax):
 ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False); ax.grid(axis='y',color=GRID,lw=.5,alpha=.8); ax.set_axisbelow(True)
def lab(ax,s): ax.text(-.14,1.05,s,transform=ax.transAxes,fontsize=9.3,fontweight='bold',color=INK)
def save(fig):
 for ext,dpi in [('png',600),('tiff',600),('pdf',None),('svg',None)]:
  fig.savefig(OUT/f'Figure6_topology_controls.{ext}',bbox_inches='tight',facecolor='white',**({'dpi':dpi} if dpi else {}))
 plt.close(fig)

def read_metrics():
 rows=[]
 for label,globname in [('No-edge','no_prior_gat_no_edge_*'),('Identity','no_prior_gat_identity_*'),('Rewired','no_prior_gat_shuffled_*'),('Pearson','no_prior_gat_pearson_*'),('Spearman','no_prior_gat_spearman_*')]:
  p=list((HPC/'no_prior_gat_fixed').glob(f'{globname}/no_prior_gat_metrics.csv'))[0]
  q=pd.read_csv(p); q['topology_label']=label; rows.append(q)
 return pd.concat(rows,ignore_index=True)

def main():
 d=read_metrics(); order=['No-edge','Identity','Rewired','Pearson','Spearman']; dirs=['HMP2_to_PRISM','PRISM_to_HMP2'];
 fig,axes=plt.subplots(2,2,figsize=(8.8,6.2)); fig.suptitle('Topology ablations localize the relational contribution',fontsize=10.2,fontweight='bold',y=.99)
 # a split-level both directions
 ax=axes[0,0]; lab(ax,'a'); ax.set_title('Split-level AUROC across topology controls',loc='left'); style(ax)
 cols={'No-edge':'#7C8892','Identity':'#7C8892','Rewired':RED,'Pearson':ORANGE,'Spearman':PURPLE}
 for di,direction in enumerate(dirs):
  for j,t in enumerate(order):
   q=d[(d.direction==direction)&(d.topology_label==t)].sort_values('split'); y=di*6+j+np.linspace(-.10,.10,len(q)); ax.scatter(q.roc_auc,y,s=13,color=cols[t],alpha=.9)
   ax.plot([q.roc_auc.mean()-q.roc_auc.std(ddof=1),q.roc_auc.mean()+q.roc_auc.std(ddof=1)],[di*6+j]*2,color=cols[t],lw=1.7)
 ax.axvline(.5,color=MUTED,ls='--',lw=.8); ax.set_yticks([2,8],['HMP2 → PRISM','PRISM → HMP2']); ax.set_xlabel('AUROC; points are source repeats'); ax.set_xlim(.3,.95)
 ax.legend([Line2D([0],[0],marker='o',color=cols[t],lw=0,markersize=5) for t in order],order,frameon=False,fontsize=5.6,ncol=2,loc='lower right')
 # b heatmap
 ax=axes[0,1]; lab(ax,'b'); ax.set_title('Mean AUROC by topology and direction',loc='left')
 mat=np.array([[d[(d.topology_label==t)&(d.direction==di)].roc_auc.mean() for di in dirs] for t in order]); im=ax.imshow(mat,cmap=mpl.colors.LinearSegmentedColormap.from_list('topo',['#F5F7F8',NAVY,TEAL]),vmin=.45,vmax=.88,aspect='auto')
 ax.set_xticks(range(2),['HMP2 →\nPRISM','PRISM →\nHMP2']); ax.set_yticks(range(5),order); ax.set_xlabel('Transfer direction')
 for i in range(5):
  for j in range(2): ax.text(j,i,f'{mat[i,j]:.3f}',ha='center',va='center',fontsize=6.6,color='white' if mat[i,j]>.68 else INK,fontweight='bold')
 # c edge Jaccard
 ax=axes[1,0]; lab(ax,'c'); ax.set_title('Edge-set overlap is moderate only for Pearson–Spearman',loc='left'); style(ax)
 labels=['Pearson vs Spearman\nHMP2 → PRISM','Pearson vs Spearman\nPRISM → HMP2','Spearman vs rewired\nHMP2 → PRISM','Spearman vs rewired\nPRISM → HMP2']; vals=[.325,.415,.054,.043]; cols2=[PURPLE,PURPLE,ORANGE,ORANGE]; y=np.arange(4)[::-1]
 ax.hlines(y,0,vals,color=cols2,lw=4,alpha=.55); ax.scatter(vals,y,color=cols2,s=32,zorder=3); ax.set_yticks(y,labels); ax.set_xlim(0,.5); ax.set_xlabel('Edge-set Jaccard similarity'); ax.grid(axis='x',color=GRID,lw=.5); ax.text(.49,.15,'J(Ea,Eb) = |intersection| / |union|',ha='right',va='center',fontsize=6.4,color=MUTED)
 # d degree diagnostic schematic
 ax=axes[1,1]; lab(ax,'d'); ax.set_title('Out-degree is fixed by the construction',loc='left'); ax.axis('off')
 ax.add_patch(FancyBboxPatch((.10,.65),.80,.20,boxstyle='round,pad=.02,rounding_size=.03',fc='#F0EBF7',ec=PURPLE,lw=1.2,transform=ax.transAxes)); ax.text(.50,.75,'Each node retains exactly k = 8 outgoing edges',transform=ax.transAxes,ha='center',va='center',fontsize=8,color=INK)
 ax.add_patch(FancyArrowPatch((.50,.62),(.50,.46),transform=ax.transAxes,arrowstyle='-|>',mutation_scale=10,color=MUTED,lw=1.2));
 ax.add_patch(FancyBboxPatch((.10,.25),.80,.20,boxstyle='round,pad=.02,rounding_size=.03',fc='#E7F0F7',ec=NAVY,lw=1.2,transform=ax.transAxes)); ax.text(.50,.35,'conventional degree correlation: undefined',transform=ax.transAxes,ha='center',va='center',fontsize=8,color=INK,fontweight='bold')
 ax.text(.50,.13,'The audit value 1.0 is retained only as a construction-level placeholder; it does not establish conserved ecological topology.',transform=ax.transAxes,ha='center',va='center',fontsize=6.8,color=MUTED,wrap=True)
 fig.tight_layout(rect=[.02,.02,.99,.95]); save(fig)
main()
