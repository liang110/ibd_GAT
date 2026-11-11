#!/usr/bin/env python
"""Source-only shared-space cross-cohort No-Prior GAT ablation."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np, pandas as pd, torch
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from ibd_microbiome_framework.data.matrices import load_canonical_microbiome, load_shared_feature_list
from ibd_microbiome_framework.evaluation.splits import repeated_group_train_validation
from ibd_microbiome_framework.preprocessing.compositional import CLRFeaturePipeline
from ibd_microbiome_framework.models.gat import NoPriorGAT

def edges(x, mode, seed):
    n=x.shape[1]; rng=np.random.default_rng(seed)
    if mode=='identity': return torch.stack([torch.arange(n),torch.arange(n)])
    if mode=='no_edge': return torch.empty((2,0),dtype=torch.long)
    c=np.corrcoef(x.T) if mode=='pearson' else pd.DataFrame(x).corr(method='spearman').to_numpy()
    np.fill_diagonal(c,0); k=min(8,n-1); out=[]
    for i in range(n):
        js=np.argsort(np.abs(c[i]))[-k:]
        out.extend((i,int(j)) for j in js if abs(c[i,j])>0)
    if mode=='shuffled':
        # Null topology: preserve each source node's outgoing degree but
        # randomly rewire destinations. Shuffling edge row order is not a null.
        rewired=[]
        for i,j in out:
            dest=int(rng.integers(0,n-1)); dest += int(dest >= i)
            rewired.append((i,dest))
        out=rewired
    return torch.tensor(out,dtype=torch.long).T.contiguous() if out else torch.empty((2,0),dtype=torch.long)

def make_graph(x, y, ei):
    return Data(x=torch.tensor(x[:,None],dtype=torch.float32), y=torch.tensor([float(y)]), edge_index=ei)

def fit_gat(x,y,ei,seed,epochs=40):
    torch.manual_seed(seed); model=NoPriorGAT(); opt=torch.optim.Adam(model.parameters(),lr=1e-3,weight_decay=1e-3)
    loader=DataLoader([make_graph(a,b,ei) for a,b in zip(x,y)],batch_size=32,shuffle=True)
    model.train()
    for _ in range(epochs):
        for batch in loader:
            opt.zero_grad(); loss=torch.nn.functional.binary_cross_entropy_with_logits(model(batch.x,batch.edge_index,batch.batch),batch.y.view(-1)); loss.backward(); opt.step()
    return model

def predict(model,x,ei):
    model.eval(); loader=DataLoader([make_graph(a,0,ei) for a in x],batch_size=64)
    with torch.no_grad(): return torch.cat([torch.sigmoid(model(b.x,b.edge_index,b.batch)).cpu() for b in loader]).numpy()

def main():
    p=argparse.ArgumentParser(); p.add_argument('--contract-dir',type=Path,required=True); p.add_argument('--legacy-root',type=Path,required=True); p.add_argument('--output-dir',type=Path,required=True); p.add_argument('--n-splits',type=int,default=5); p.add_argument('--topology',default='pearson',choices=['pearson','spearman','no_edge','identity','shuffled']); a=p.parse_args(); a.output_dir.mkdir(parents=True,exist_ok=True)
    shared=load_shared_feature_list(a.contract_dir); hm=pd.read_csv(a.contract_dir/'cohort_manifest_hmp2.csv'); pm=pd.read_csv(a.contract_dir/'cohort_manifest_prism.csv')
    h=load_canonical_microbiome(a.legacy_root/'manual_abundance.csv',hm,shared); q=load_canonical_microbiome(a.legacy_root/'prism_abundance_clean.csv',pm,shared)
    hm=hm.set_index('sample_id', drop=False).loc[h.index].reset_index(drop=True); pm=pm.set_index('sample_id', drop=False).loc[q.index].reset_index(drop=True); rows=[]; predictions=[]
    for direction,source,sm,target,tm in [('HMP2_to_PRISM',h,hm,q,pm),('PRISM_to_HMP2',q,pm,h,hm)]:
        for split,tr,va in repeated_group_train_validation(sm,a.n_splits,.2,20261110):
            ti=sm.iloc[tr].sample_id.tolist(); pipe=CLRFeaturePipeline(top_k=min(100,len(shared)),random_state=20261110+split).fit(source.loc[ti],sm.iloc[tr].outcome_ibd); xs=pipe.transform(source.loc[ti]); xt=pipe.transform(target); ei=edges(xs.to_numpy(),'pearson' if a.topology=='pearson' else a.topology,20261110+split); ytr=sm.iloc[tr].outcome_ibd.to_numpy(); y=tm.outcome_ibd.to_numpy()
            gat=predict(fit_gat(xs.to_numpy(),ytr,ei,20261110+split),xt.to_numpy(),ei)
            en=LogisticRegression(penalty='elasticnet',solver='saga',l1_ratio=.5,C=1.0,class_weight='balanced',max_iter=3000,random_state=20261110+split).fit(xs,ytr).predict_proba(xt)[:,1]
            rf=RandomForestClassifier(n_estimators=300,min_samples_leaf=2,class_weight='balanced',random_state=20261110+split,n_jobs=1).fit(xs,ytr).predict_proba(xt)[:,1]
            for name,prob in [('no_prior_gat',gat),('elastic_net',en),('random_forest',rf)]:
                rows.append({'direction':direction,'split':split,'topology':a.topology,'model':name,'roc_auc':roc_auc_score(y,prob),'average_precision':average_precision_score(y,prob),'brier':brier_score_loss(y,prob),'n_train':len(tr),'n_target':len(y),'n_features':xs.shape[1],'edge_hash':hashlib.sha256(ei.numpy().tobytes()).hexdigest(),'edge_count':ei.shape[1]})
                predictions.extend({'direction':direction,'split':split,'model':name,'sample_id':sid,'participant_id':str(pid),'outcome':int(yy),'probability':float(pp)} for sid,pid,yy,pp in zip(tm.sample_id,tm.participant_id,y,prob))
    frame=pd.DataFrame(rows); frame.to_csv(a.output_dir/'no_prior_gat_metrics.csv',index=False); pd.DataFrame(predictions).to_csv(a.output_dir/'paired_predictions.csv',index=False); summary=frame.groupby(['direction','topology','model'])[['roc_auc','average_precision','brier']].agg(['mean','std']).reset_index(); summary.to_csv(a.output_dir/'no_prior_gat_summary.csv',index=False); (a.output_dir/'manifest.json').write_text(json.dumps({'no_global_prior':True,'shared_features':len(shared),'source_only_graph':True,'n_splits':a.n_splits,'topology':a.topology,'paired_comparators':['elastic_net','random_forest']},indent=2))
if __name__=='__main__': main()
