from pathlib import Path
import json
import numpy as np, pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'manuscript_rebuild'/'figures'; DATA=ROOT/'manuscript_rebuild'/'source_data'
OUT.mkdir(parents=True,exist_ok=True); DATA.mkdir(parents=True,exist_ok=True)
mpl.rcParams.update({'font.family':'DejaVu Sans','font.size':7.2,'axes.titlesize':8.2,'axes.labelsize':7.2,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.6,'pdf.fonttype':42,'svg.fonttype':'none','legend.frameon':False,'xtick.major.size':2.5,'ytick.major.size':2.5})
NAVY='#2F5D85'; TEAL='#3BA7A0'; ORANGE='#D9792E'; PURPLE='#7E5AA6'; RED='#B54E4E'; GREY='#66737F'; LIGHT='#F5F7F8'; BLUEPALE='#E7F0F7'; TEALPALE='#E6F5F2'; ORANGEPALE='#FBEFE4'; PURPLEPALE='#F0EBF7';

def panel(ax,label,title):
    ax.text(-0.13,1.05,label,transform=ax.transAxes,fontweight='bold',fontsize=9.5,color='#18232D',va='top')
    ax.set_title(title,loc='left',pad=4,fontweight='bold',color='#18232D')
def save(fig,name):
    base=OUT/name
    fig.savefig(str(base)+'.png',dpi=600,bbox_inches='tight',facecolor='white')
    fig.savefig(str(base)+'.pdf',bbox_inches='tight',facecolor='white')
    fig.savefig(str(base)+'.svg',bbox_inches='tight',facecolor='white')
    plt.close(fig)

def fig2():
    scale=pd.read_csv(ROOT/'artifacts/hpc_results/gate0_20261110_000415/input_scale_audit.csv')
    fig=plt.figure(figsize=(7.2,4.5)); gs=GridSpec(2,3,figure=fig,hspace=.78,wspace=.48)
    ax=fig.add_subplot(gs[0,0]); panel(ax,'a','Input scale range')
    x=np.arange(2); med=scale.sum_median.values; lo=scale.sum_min.values; hi=scale.sum_max.values
    ax.vlines(x,lo,hi,color=[NAVY,TEAL],lw=3,alpha=.75); ax.scatter(x,med,color=[NAVY,TEAL],s=32,zorder=3,edgecolor='white',lw=.7)
    ax.set_xticks(x,scale.cohort); ax.set_ylabel('Library sum'); ax.set_ylim(-5,110); ax.set_yticks([0,50,100]);
    ax.text(0,103,'percentage-like',ha='center',fontsize=6.4,color=NAVY); ax.text(1,6,'proportion-like',ha='center',fontsize=6.4,color=TEAL)
    ax=fig.add_subplot(gs[0,1]); panel(ax,'b','Feature contract')
    mat=np.array([[578,126],[209,126]])
    im=ax.imshow(mat,cmap=LinearSegmentedColormap.from_list('c',['#edf3f6','#2f5d85']),vmin=0,vmax=578,aspect='auto')
    ax.set_xticks([0,1],['Raw','Shared']); ax.set_yticks([0,1],['HMP2','PRISM'])
    for i in range(2):
      for j in range(2): ax.text(j,i,str(mat[i,j]),ha='center',va='center',color='white' if mat[i,j]>300 else '#18232D',fontweight='bold')
    ax.set_xlabel('Microbial features');
    ax.text(.5,-.33,'126-feature primary contract',transform=ax.transAxes,ha='center',fontsize=6.1,color=ORANGE)
    ax=fig.add_subplot(gs[0,2]); panel(ax,'c','Metabolite intersection')
    ax.barh(['Exact shared names'],[0],color=ORANGE,height=.42); ax.set_xlim(-.02,1); ax.set_xticks([0,1]); ax.set_xlabel('Count'); ax.text(.02,0,'0',va='center',ha='left',fontweight='bold',color=ORANGE)
    ax.spines['left'].set_visible(False); ax.spines['bottom'].set_position(('data',-.3));
    ax=fig.add_subplot(gs[1,0]); panel(ax,'d','Samples and exclusions')
    ns=[1638,220]; exc=[3,0]; y=np.arange(2)
    ax.barh(y,ns,color=[NAVY,TEAL],alpha=.85,height=.5); ax.set_yticks(y,['HMP2','PRISM']); ax.set_xlabel('Samples'); ax.set_xlim(0,1750)
    for yy,n,e in zip(y,ns,exc):
      if n > 1000:
        ax.text(n-25,yy,f'n={n}; excluded={e}',va='center',ha='right',fontsize=6.3,color='white',fontweight='bold')
      else:
        ax.text(n+30,yy,f'n={n}; excluded={e}',va='center',fontsize=6.3)
    ax=fig.add_subplot(gs[1,1:]); panel(ax,'e','Analysis boundary')
    ax.axis('off');
    boxes=[(0.02,.48,.22,.26,'HMP2\n1,638 samples\n578 raw features',BLUEPALE,NAVY),(0.02,.12,.22,.26,'PRISM\n220 samples\n209 raw features',TEALPALE,TEAL),(0.38,.30,.24,.32,'Shared microbiome\n126 features',ORANGEPALE,ORANGE),(0.78,.30,.18,.32,'Primary cross-cohort\nmicrobiome only',PURPLEPALE,PURPLE)]
    for x0,y0,w,h,t,fc,ec in boxes:
      p=FancyBboxPatch((x0,y0),w,h,boxstyle='round,pad=.02,rounding_size=.03',fc=fc,ec=ec,lw=1,transform=ax.transAxes);ax.add_patch(p);ax.text(x0+w/2,y0+h/2,t,transform=ax.transAxes,ha='center',va='center',fontsize=7.0,color='#18232D')
    for a,b in [((.24,.61),(.38,.47)),((.24,.25),(.38,.45)),((.62,.46),(.78,.46))]: ax.add_patch(FancyArrowPatch(a,b,transform=ax.transAxes,arrowstyle='-|>',mutation_scale=8,lw=1,color=GREY))
    fig.suptitle('Input contract audit',x=.03,ha='left',fontsize=10.5,fontweight='bold',color='#18232D')
    save(fig,'Figure2_input_contract'); scale.to_csv(DATA/'Figure2_input_scale_audit.csv',index=False)

def fig3():
    leak=pd.read_csv(ROOT/'artifacts/hpc_results/repeated_leakage_20261110_001007/repeated_leakage_metrics.csv')
    synth=pd.read_csv(ROOT/'artifacts/hpc_results/synthetic_supplement/synthetic_controls_v2_20261110_091610/synthetic_control_metrics.csv',keep_default_na=False)
    fig=plt.figure(figsize=(7.3,4.5)); gs=GridSpec(2,3,figure=fig,hspace=.58,wspace=.55)
    ax=fig.add_subplot(gs[0,0]); panel(ax,'a','Participant leakage')
    vals=[leak.random_split_roc_auc.values,leak.group_split_roc_auc.values]
    ax.violinplot(vals,positions=[0,1],showmeans=False,showmedians=True,widths=.65)
    for i,v in enumerate(vals): ax.scatter(np.full(len(v),i),v,s=10,color=[RED,NAVY][i],alpha=.75,zorder=3)
    ax.axhline(.5,color='black',lw=.6,ls='--'); ax.set_xticks([0,1],['Naive split','Grouped split']); ax.set_ylabel('AUROC'); ax.set_ylim(.15,1.05)
    ax=fig.add_subplot(gs[0,1]); panel(ax,'b','Overlap audit')
    ax.scatter(np.zeros(10),leak.random_split_subject_overlap,color=RED,s=20,label='Naive'); ax.scatter(np.ones(10),leak.group_split_subject_overlap,color=NAVY,s=20,label='Grouped');
    ax.set_xticks([0,1],['Naive','Grouped']); ax.set_ylabel('Shared participants'); ax.set_ylim(-2,50); ax.legend(fontsize=6,loc='upper right')
    ax=fig.add_subplot(gs[0,2]); panel(ax,'c','Clinical leakage validator')
    nrej=int(leak.clinical_leakage_rejected.sum()); ax.barh(['Rejected'],[nrej],color=TEAL,height=.42); ax.set_xlim(0,10); ax.set_xlabel('Seeds (n=10)'); ax.text(nrej-.25,0,f'{nrej}/10',ha='right',va='center',color='white',fontweight='bold')
    ax=fig.add_subplot(gs[1,:2]); panel(ax,'d','Known-truth controls')
    summ=synth.groupby('scenario').agg(auc=('target_roc_auc','mean'),rec=('n_disease_features_recovered','mean'),pass_rate=('acceptance_pass','mean')).reindex(['null','batch_only','signal_only','signal_reversing_batch','signal_topology_shift'])
    im=ax.imshow(summ[['auc','rec','pass_rate']].values.T,aspect='auto',cmap='YlGnBu',vmin=0,vmax=1)
    ax.set_xticks(range(len(summ)),['Null','Batch','Signal','Reversed\nbatch','Topology\nshift'],rotation=25,ha='right'); ax.set_yticks([0,1,2],['AUROC','Recovery / 6','All-pass rate'])
    z=summ.copy(); z['rec']=z.rec/6
    arr=z[['auc','rec','pass_rate']].values.T
    for i in range(arr.shape[0]):
      for j in range(arr.shape[1]): ax.text(j,i,f'{arr[i,j]:.2f}',ha='center',va='center',fontsize=6.2,color='white' if arr[i,j]>.55 else '#18232D')
    ax=fig.add_subplot(gs[1,2]); panel(ax,'e','Effect-sweep summary')
    sw=pd.read_csv(ROOT/'artifacts/hpc_results/synthetic_supplement/synthetic_sweep_20261110_091610/synthetic_sweep_metrics.csv',keep_default_na=False)
    ax.scatter(sw.n_disease_features_recovered,sw.target_roc_auc,c=sw.all_acceptance_pass.map({True:TEAL,False:RED}),s=22,edgecolor='white',lw=.3); ax.axhline(.70,color='black',lw=.6,ls='--'); ax.set_xlabel('Recovered disease features'); ax.set_ylabel('AUROC'); ax.set_xlim(0,6.5); ax.set_ylim(.35,.85)
    fig.suptitle('Leakage and known-truth controls',x=.03,ha='left',fontsize=10.5,fontweight='bold',color='#18232D')
    save(fig,'Figure3_controls'); leak.to_csv(DATA/'Figure3_repeated_leakage.csv',index=False); synth.to_csv(DATA/'Figure3_synthetic_controls.csv',index=False)

def fig4():
    d=pd.read_csv(ROOT/'artifacts/hpc_results/closure_final/origin_paired_20261110_131203/origin_paired_metrics.csv')
    bs=pd.read_csv(ROOT/'artifacts/hpc_results/closure_final/origin_paired_20261110_131203/origin_paired_bootstrap.csv')
    fig=plt.figure(figsize=(7.4,4.2)); gs=GridSpec(2,2,figure=fig,hspace=.62,wspace=.52)
    ax=fig.add_subplot(gs[0,0]); panel(ax,'a','Origin-held-out AUROC')
    dirs=['PRISM_discovery_to_validation','PRISM_validation_to_discovery']; methods=['no_prior_gat','elastic_net','random_forest']; colors=[PURPLE,ORANGE,NAVY]
    xpos=np.arange(2); w=.22
    for k,m in enumerate(methods):
      means=d[d.model==m].groupby('direction').roc_auc.mean().reindex(dirs).values; std=d[d.model==m].groupby('direction').roc_auc.std().reindex(dirs).values
      ax.errorbar(xpos+(k-1)*w,means,yerr=std,fmt='o',capsize=2.5,ms=4,color=colors[k],label={'no_prior_gat':'No-Prior GAT','elastic_net':'Elastic net','random_forest':'Random forest'}[m])
    ax.axhline(.5,color='black',lw=.5,ls='--'); ax.set_xticks(xpos,['Discovery→\nvalidation','Validation→\ndiscovery']); ax.set_ylabel('AUROC'); ax.set_ylim(.35,1); ax.legend(fontsize=5.8)
    ax=fig.add_subplot(gs[0,1]); panel(ax,'b','Split-specific GAT−reference')
    for k,comp in enumerate(['gat_vs_elastic_net','gat_vs_random_forest']):
      for j,di in enumerate(dirs):
        q=bs[(bs.comparison==comp)&(bs.direction==di)].sort_values('split'); y=np.arange(len(q))+j*12
        for yy,dl,lo,hi in zip(y,q.delta_auc,q.ci_low,q.ci_high): ax.plot([lo,hi],[yy,yy],color=[ORANGE,NAVY][k],lw=.9); ax.scatter(dl,yy,color=[ORANGE,NAVY][k],s=10)
    ax.axvline(0,color='black',lw=.6); ax.set_xlabel('ΔAUROC (GAT − reference)'); ax.set_yticks([4.5,16.5],['D→V','V→D']); ax.set_ylim(-1,23); ax.set_yticklabels(['D→V','V→D']); ax.tick_params(axis='y',length=0)
    ax=fig.add_subplot(gs[1,0]); panel(ax,'c','Origin transfer matrix')
    meanmat=np.array([[d[(d.direction==di)&(d.model==m)].roc_auc.mean() for m in methods] for di in dirs])
    im=ax.imshow(meanmat,vmin=.4,vmax=.9,cmap='PuBuGn'); ax.set_xticks(range(3),['GAT','Elastic\nnet','Random\nforest']); ax.set_yticks(range(2),['D→V','V→D'])
    for i in range(2):
      for j in range(3): ax.text(j,i,f'{meanmat[i,j]:.3f}',ha='center',va='center',fontsize=6.5,color='white' if meanmat[i,j]>.65 else '#18232D')
    ax.set_xlabel('Model');
    ax=fig.add_subplot(gs[1,1]); panel(ax,'d','Direction contrast')
    contrasts=[]
    for m in methods:
      a=d[(d.direction==dirs[0])&(d.model==m)].roc_auc.mean(); b=d[(d.direction==dirs[1])&(d.model==m)].roc_auc.mean(); contrasts.append(b-a)
    ax.barh(range(3),contrasts,color=colors,height=.5); ax.axvline(0,color='black',lw=.6); ax.set_yticks(range(3),['GAT','Elastic net','Random forest']); ax.set_xlabel('V→D minus D→V AUROC');
    fig.suptitle('PRISM origin-held-out transfer is direction-dependent',x=.03,ha='left',fontsize=10.5,fontweight='bold',color='#18232D')
    save(fig,'Figure4_origin_transfer'); d.to_csv(DATA/'Figure4_origin_metrics.csv',index=False); bs.to_csv(DATA/'Figure4_origin_bootstrap.csv',index=False)

def fig5():
    d=pd.read_csv(ROOT/'artifacts/hpc_results/final_supplement/no_prior_gat_spearman_20261110_082627/no_prior_gat_metrics.csv')
    bs=pd.read_csv(ROOT/'artifacts/hpc_results/final_supplement/paired_bootstrap_20261110_082627/paired_bootstrap_by_split.csv')
    fig=plt.figure(figsize=(7.4,4.3)); gs=GridSpec(2,3,figure=fig,hspace=.62,wspace=.52)
    colors={'no_prior_gat':PURPLE,'elastic_net':ORANGE,'random_forest':NAVY}; names={'no_prior_gat':'No-Prior GAT','elastic_net':'Elastic net','random_forest':'Random forest'}; dirs=['HMP2_to_PRISM','PRISM_to_HMP2']
    ax=fig.add_subplot(gs[0,:2]); panel(ax,'a','Ten-repeat transfer distribution')
    for i,di in enumerate(dirs):
      for j,m in enumerate(['no_prior_gat','elastic_net','random_forest']):
        q=d[(d.direction==di)&(d.model==m)].roc_auc.values; y=np.full(len(q),i*4+j)+np.linspace(-.12,.12,len(q)); ax.scatter(q,y,s=13,color=colors[m],alpha=.8); ax.plot([q.mean()-q.std(),q.mean()+q.std()],[i*4+j]*2,color=colors[m],lw=2)
    ax.axvline(.5,color='black',lw=.5,ls='--'); ax.set_yticks([1,5],['HMP2→PRISM','PRISM→HMP2']); ax.set_ylim(-1,7); ax.set_xlabel('AUROC'); ax.legend([plt.Line2D([0],[0],marker='o',color='w',markerfacecolor=colors[m],label=names[m]) for m in colors],names.values(),fontsize=5.8,ncol=3,loc='lower right')
    ax=fig.add_subplot(gs[0,2]); panel(ax,'b','Mean transfer metrics')
    means=d.groupby(['direction','model']).roc_auc.mean().unstack(); im=ax.imshow(means.values,vmin=.55,vmax=.9,cmap='PuBuGn'); ax.set_xticks(range(3),['GAT','Elastic\nnet','Random\nforest']); ax.set_yticks(range(2),['HMP2→PRISM','PRISM→HMP2'])
    for i in range(2):
      for j in range(3): ax.text(j,i,f'{means.iloc[i,j]:.3f}',ha='center',va='center',fontsize=6.5,color='white' if means.iloc[i,j]>.7 else '#18232D')
    ax=fig.add_subplot(gs[1,:2]); panel(ax,'c','Paired ΔAUROC across source repeats')
    for j,di in enumerate(dirs):
      for k,comp in enumerate(['no_prior_gat_vs_elastic_net','no_prior_gat_vs_random_forest']):
        q=bs[(bs.direction==di)&(bs.comparison==comp)].sort_values('split'); y=np.arange(len(q))+j*12+k*0.25
        ax.scatter(q.delta_auc,y,s=12,color=[ORANGE,NAVY][k],label=comp if j==0 else None)
    ax.axvline(0,color='black',lw=.6); ax.set_yticks([4.5,16.5],['HMP2→PRISM','PRISM→HMP2']); ax.set_xlabel('ΔAUROC (GAT − reference)'); ax.legend(fontsize=5.8,loc='lower right')
    ax=fig.add_subplot(gs[1,2]); panel(ax,'d','PR-AUC direction contrast')
    ap=d.groupby(['direction','model']).average_precision.mean().unstack(); x=np.arange(2); w=.24
    for k,m in enumerate(['no_prior_gat','elastic_net','random_forest']): ax.plot(x+(k-1)*w,ap[m].values,'o-',color=colors[m],ms=4,label=names[m])
    ax.set_xticks(x,['HMP2→PRISM','PRISM→HMP2']); ax.set_ylim(.7,1); ax.set_ylabel('Average precision'); ax.legend(fontsize=5.3)
    fig.suptitle('Directional cross-cohort transfer under a frozen microbiome contract',x=.03,ha='left',fontsize=10.5,fontweight='bold',color='#18232D')
    save(fig,'Figure5_transfer'); d.to_csv(DATA/'Figure5_transfer_metrics.csv',index=False); bs.to_csv(DATA/'Figure5_transfer_bootstrap.csv',index=False)

def fig6():
    # Delegate to the final topology-only renderer.  The previous composite
    # prototype mixed topology and capacity panels and is retained only in the
    # archive; it must not regenerate the submission asset.
    import runpy
    runpy.run_path(str(Path(__file__).with_name('plot_topology_final.py')), run_name='__main__')
    return
    topo_files={'Spearman':'no_prior_gat_fixed/no_prior_gat_spearman_20261110_075752/no_prior_gat_metrics.csv','Pearson':'no_prior_gat_fixed/no_prior_gat_pearson_20261110_075751/no_prior_gat_metrics.csv','Rewired':'no_prior_gat_fixed/no_prior_gat_shuffled_20261110_075752/no_prior_gat_metrics.csv','No-edge':'no_prior_gat_fixed/no_prior_gat_no_edge_20261110_075752/no_prior_gat_metrics.csv','Identity':'no_prior_gat_fixed/no_prior_gat_identity_20261110_075752/no_prior_gat_metrics.csv'}
    dfs=[]
    for t,p in topo_files.items():
      q=pd.read_csv(ROOT/'artifacts/hpc_results'/p); q['topology_label']=t; dfs.append(q)
    d=pd.concat(dfs,ignore_index=True)
    cap=pd.read_csv(ROOT/'artifacts/hpc_results/capacity_control/capacity_20261110_103951/capacity_control_metrics.csv')
    fig=plt.figure(figsize=(7.5,4.7)); gs=GridSpec(2,3,figure=fig,hspace=.78,wspace=.52)
    ax=fig.add_subplot(gs[0,:2]); panel(ax,'a','Topology ablation, AUROC by repeat')
    order=['No-edge','Identity','Rewired','Pearson','Spearman']; cols=[GREY,GREY,RED,ORANGE,PURPLE]
    for j,t in enumerate(order):
      q=d[(d.topology_label==t)&(d.direction=='HMP2_to_PRISM')].sort_values('split'); y=np.full(len(q),j)+np.linspace(-.13,.13,len(q)); ax.scatter(q.roc_auc,y,s=13,color=cols[j]); ax.plot([q.roc_auc.mean()-q.roc_auc.std(),q.roc_auc.mean()+q.roc_auc.std()],[j,j],color=cols[j],lw=2)
    ax.axvline(.5,color='black',lw=.5,ls='--'); ax.set_yticks(range(5),order); ax.set_xlabel('AUROC, HMP2→PRISM'); ax.set_xlim(.3,.95)
    ax=fig.add_subplot(gs[0,2]); panel(ax,'b','Topology × direction')
    mean=d.groupby(['topology_label','direction']).roc_auc.mean().unstack().reindex(order); im=ax.imshow(mean.values,vmin=.45,vmax=.9,cmap='PuBuGn'); ax.set_xticks([0,1],['H→P','P→H']); ax.set_xlabel('Direction'); ax.set_yticks(range(5),order); ax.tick_params(axis='x',pad=4)
    for i in range(5):
      for j in range(2): ax.text(j,i,f'{mean.iloc[i,j]:.3f}',ha='center',va='center',fontsize=6,color='white' if mean.iloc[i,j]>.65 else '#18232D')
    ax=fig.add_subplot(gs[1,0]); panel(ax,'c','Capacity-matched comparator')
    cm=cap.groupby('direction').roc_auc.agg(['mean','std']).reindex(['HMP2_to_PRISM','PRISM_to_HMP2']); gat=d[(d.topology_label=='Spearman')].groupby('direction').roc_auc.agg(['mean','std']).reindex(cm.index)
    x=np.arange(2); ax.errorbar(x-.12,gat['mean'],yerr=gat['std'],fmt='o',color=PURPLE,capsize=2,label='GAT'); ax.errorbar(x+.12,cm['mean'],yerr=cm['std'],fmt='o',color=TEAL,capsize=2,label='MLP'); ax.set_xticks(x,['HMP2→PRISM','PRISM→HMP2']); ax.set_ylim(.45,.9); ax.set_ylabel('AUROC'); ax.legend(fontsize=5.8)
    ax=fig.add_subplot(gs[1,1]); panel(ax,'d','Edge-set Jaccard')
    x=['HMP2→PRISM','PRISM→HMP2']; ps=[.325,.415]; rw=[.054,.043]; xx=np.arange(2); ax.plot(xx,ps,'o-',color=PURPLE,label='Pearson vs Spearman'); ax.plot(xx,rw,'s--',color=RED,label='Spearman vs rewired'); ax.set_xticks(xx,x); ax.set_ylim(0,.5); ax.set_ylabel('Jaccard'); ax.legend(fontsize=5.6,loc='upper left')
    ax=fig.add_subplot(gs[1,2]); panel(ax,'e','Construction diagnostic')
    ax.bar(['Pearson','Spearman','Rewired'],[1,1,1],color=[ORANGE,PURPLE,RED],alpha=.85); ax.set_ylim(0,1.15); ax.set_ylabel('Out-degree corr.'); ax.set_yticks([0,1]); ax.text(1,1.04,'fixed k=8',ha='center',fontsize=6.6,color=GREY)
    fig.suptitle('Topology and capacity controls localize the relational contribution',x=.03,ha='left',fontsize=10.5,fontweight='bold',color='#18232D')
    save(fig,'Figure6_topology_capacity'); d.to_csv(DATA/'Figure6_topology_metrics.csv',index=False); cap.to_csv(DATA/'Figure6_capacity_metrics.csv',index=False)

def main():
    fig2(); fig3(); fig4(); fig5(); fig6()
    captions={
      'Figure1':'Source-only, confound-audited framework for directional transfer. The harmonized HMP2/IBDMDB and PRISM microbiome matrices are reduced to a frozen 126-feature contract before split-local source fitting; selected features enter tabular and No-Prior graph-attention branches, with targets locked until prediction.',
      'Figure2':'Input scale, feature harmonization and the 126-feature cross-cohort contract. Ranges show min–max library sums with median markers; all values are from Gate0.',
      'Figure3':'Participant overlap and known-truth controls. Each point is a seed or scenario replicate; violin medians and heatmap values summarize the audit package.',
      'Figure4':'PRISM origin-held-out transfer. Error bars show across-split standard deviation; paired intervals are plotted separately for each source split and are not pooled intervals.',
      'Figure5':'Directional cross-cohort transfer. Points are ten source repeats; horizontal segments show mean ± s.d.; paired differences are shown per split.',
      'Figure6':'Topology ablations and edge-set diagnostics. (a) Split-level AUROC values across no-edge, identity, rewired, Pearson and Spearman topologies in both transfer directions. (b) Mean AUROC by topology and direction. (c) Edge-set Jaccard overlap relative to the Spearman reference for Pearson and rewired graphs. (d) Fixed top-k construction assigns eight outgoing edges to every retained node; a conventional degree correlation is undefined under this construction and the audit placeholder is not interpreted as ecological conservation.',
      'Figure7':'Capacity-matched no-message comparator. Paired transfer performance of the approximate no-message-passing multilayer perceptron and No-Prior Spearman GAT using the same source-only transformed input contract and similar hidden capacity.'}
    (DATA/'figure_captions.json').write_text(json.dumps(captions,indent=2),encoding='utf-8')
if __name__=='__main__': main()
