from pathlib import Path
import pandas as pd, numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap

ROOT = Path(r'D:\IBD_GAT')
OUT = ROOT / 'manuscript_rebuild' / 'figures'
DATA = ROOT / 'manuscript_rebuild' / 'source_data'
OUT.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 7.2, 'axes.titlesize': 8.2,
    'axes.labelsize': 7.2, 'axes.spines.top': False, 'axes.spines.right': False,
    'axes.linewidth': 0.6, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
    'legend.frameon': False, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
})
NAVY='#2F5D85'; TEAL='#3BA7A0'; ORANGE='#D9792E'; PURPLE='#7E5AA6'; RED='#B54E4E'; GREY='#66737F'; GRID='#D8E1E7'; INK='#18232D'

def panel(ax, label, title):
    ax.text(-0.13, 1.05, label, transform=ax.transAxes, fontweight='bold', fontsize=9.5, color=INK, va='top')
    ax.set_title(title, loc='left', pad=4, fontweight='bold', color=INK)

def save(fig, name):
    base = OUT / name
    fig.savefig(base.with_suffix('.png'), dpi=600, bbox_inches='tight', facecolor='white')
    fig.savefig(base.with_suffix('.pdf'), bbox_inches='tight', facecolor='white')
    fig.savefig(base.with_suffix('.svg'), bbox_inches='tight', facecolor='white')
    plt.close(fig)

leak = pd.read_csv(ROOT / 'artifacts' / 'figure3_leakage_audit_20261116' / 'recomputed_current_contract_metrics.csv')
syn = pd.read_csv(ROOT / 'artifacts' / 'hpc_results' / 'synthetic_supplement' / 'synthetic_controls_v2_20261110_091610' / 'synthetic_control_metrics.csv', keep_default_na=False)
syn['scenario'] = syn['scenario'].fillna('null').replace({'signal_reversing_batch':'reversing_batch','signal_topology_shift':'topology_shift'})
summary = syn.groupby('scenario').agg(auc=('target_roc_auc','mean'), rec=('n_disease_features_recovered','mean'), pass_rate=('acceptance_pass','mean')).reindex(['null','batch_only','signal_only','reversing_batch','topology_shift'])

fig = plt.figure(figsize=(7.3, 4.5))
gs = GridSpec(2, 3, figure=fig, hspace=.68, wspace=.55)

ax = fig.add_subplot(gs[0,0]); panel(ax,'a','Participant leakage')
vals=[leak.random_split_roc_auc.values, leak.group_split_roc_auc.values]
vp=ax.violinplot(vals, positions=[0,1], showmeans=False, showmedians=True, widths=.65)
for body,col in zip(vp['bodies'],[RED,NAVY]): body.set_facecolor(col); body.set_edgecolor(col); body.set_alpha(.28)
vp['cmedians'].set_color(INK); vp['cmedians'].set_linewidth(1.0)
for i,v in enumerate(vals): ax.scatter(np.full(len(v),i),v,s=10,color=[RED,NAVY][i],alpha=.75,zorder=3,edgecolor='white',linewidth=.3)
ax.axhline(.5,color='black',lw=.6,ls='--'); ax.set_xticks([0,1],['Sample-level','Grouped']); ax.set_ylabel('AUROC'); ax.set_ylim(.15,1.05)
ax.text(0, .19, f'mean 0.960', ha='center', va='bottom', color=RED, fontsize=6.5, fontweight='bold')
ax.text(1, .19, f'mean 0.552', ha='center', va='bottom', color=NAVY, fontsize=6.5, fontweight='bold')

ax = fig.add_subplot(gs[0,1]); panel(ax,'b','Overlap audit')
ax.scatter(np.zeros(10), leak.random_split_subject_overlap, color=RED, s=20, label='Sample-level', zorder=2)
ax.scatter(np.ones(10), leak.group_split_subject_overlap, color=NAVY, s=20, label='Grouped', zorder=2)
ax.set_xticks([0,1],['Sample-level','Grouped']); ax.set_ylabel('Shared participants'); ax.set_ylim(-2,50)
ax.text(0, 44.8, 'mean 43', ha='center', va='bottom', color=RED, fontsize=6.5, fontweight='bold')
ax.text(1, 2.0, '0 in every split', ha='center', va='bottom', color=NAVY, fontsize=6.5, fontweight='bold')
ax.legend(fontsize=5.6, loc='upper right')

ax = fig.add_subplot(gs[0,2]); panel(ax,'c','Clinical-field validator')
nrej=int(leak.clinical_leakage_rejected.sum()); ax.barh(['Rejected'],[nrej],color=TEAL,height=.42); ax.set_xlim(0,10); ax.set_xlabel('Seeds (n=10)'); ax.text(nrej-.25,0,f'{nrej}/10',ha='right',va='center',color='white',fontweight='bold')

ax = fig.add_subplot(gs[1,:2]); panel(ax,'d','Known-truth controls')
order=['null','batch_only','signal_only','reversing_batch','topology_shift']; labels=['Null','Batch','Signal','Reversed\nbatch','Topology\nshift']
arr=summary[['auc','rec','pass_rate']].copy(); arr['rec']=arr['rec']/6
im=ax.imshow(arr.values.T, aspect='auto', cmap=LinearSegmentedColormap.from_list('audit',['#edf3f6','#2f5d85']), vmin=0, vmax=1)
ax.set_xticks(range(len(order)), labels, rotation=25, ha='right'); ax.set_yticks([0,1,2],['AUROC','Recovery / 6','All-pass rate'])
for i in range(arr.shape[0]):
    for j in range(arr.shape[1]):
        val=arr.iloc[i,j]; ax.text(i,j,f'{val:.2f}',ha='center',va='center',fontsize=6.1,color='white' if val>.55 else INK)

ax = fig.add_subplot(gs[1,2]); panel(ax,'e','Acceptance by scenario')
accept=syn.groupby('scenario').acceptance_pass.agg(['sum','count']).reindex(order)
colors=[GREY,GREY,TEAL,TEAL,TEAL]
for i,(sc,row) in enumerate(accept.iterrows()):
    ax.barh(i,row['count'],color='#E7EEF2',height=.55)
    ax.barh(i,row['sum'],color=colors[i],height=.55)
    ax.text(row['count']+.15,i,f"{int(row['sum'])}/{int(row['count'])}",va='center',fontsize=6.3)
ax.set_yticks(range(len(order)),labels); ax.set_xlim(0,11); ax.set_xlabel('Accepted runs')
ax.invert_yaxis()
fig.suptitle('Participant structure and known-truth controls define the validity boundary', x=.03, ha='left', fontsize=10.5, fontweight='bold', color=INK)
save(fig,'Figure3_controls_current_contract')
leak.to_csv(DATA/'Figure3_repeated_leakage_current_contract.csv',index=False)
summary.reset_index(names='scenario').to_csv(DATA/'Figure3_known_truth_summary.csv',index=False)
print('wrote', OUT/'Figure3_controls_current_contract.png')
print('random mean', leak.random_split_roc_auc.mean(), 'grouped mean', leak.group_split_roc_auc.mean(), 'overlap', leak.random_split_subject_overlap.mean(), leak.group_split_subject_overlap.max())


