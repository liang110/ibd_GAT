#!/usr/bin/env python
"""PRISM origin-held-out No-Prior Spearman GAT transfer."""
import argparse, json
from pathlib import Path
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from ibd_microbiome_framework.data.matrices import load_canonical_microbiome, load_shared_feature_list
from ibd_microbiome_framework.evaluation.splits import repeated_group_train_validation
from ibd_microbiome_framework.preprocessing.compositional import CLRFeaturePipeline
from scripts.run_no_prior_gat_transfer import edges, fit_gat, predict

def main():
 p=argparse.ArgumentParser(); p.add_argument('--contract-dir',type=Path,required=True); p.add_argument('--legacy-root',type=Path,required=True); p.add_argument('--output-dir',type=Path,required=True); p.add_argument('--n-splits',type=int,default=5); a=p.parse_args(); a.output_dir.mkdir(parents=True,exist_ok=True)
 shared=load_shared_feature_list(a.contract_dir); m=pd.read_csv(a.contract_dir/'cohort_manifest_prism.csv'); x=load_canonical_microbiome(a.legacy_root/'prism_abundance_clean.csv',m,shared); m=m.set_index('sample_id',drop=False).loc[x.index].reset_index(drop=True); rows=[]
 for so,to in [('discovery','validation'),('validation','discovery')]:
  sm=m[m.origin==so].reset_index(drop=True); tm=m[m.origin==to].reset_index(drop=True); source=x.loc[sm.sample_id]; target=x.loc[tm.sample_id]
  for split,tr,va in repeated_group_train_validation(sm,a.n_splits,.2,20261110):
   ids=sm.iloc[tr].sample_id.tolist(); pipe=CLRFeaturePipeline(top_k=min(100,len(shared)),random_state=20261110+split).fit(source.loc[ids],sm.iloc[tr].outcome_ibd); xs=pipe.transform(source.loc[ids]); xt=pipe.transform(target); ei=edges(xs.to_numpy(),'spearman',20261110+split); model=fit_gat(xs.to_numpy(),sm.iloc[tr].outcome_ibd.to_numpy(),ei,20261110+split); prob=predict(model,xt.to_numpy(),ei); y=tm.outcome_ibd.to_numpy(); rows.append({'direction':f'PRISM_{so}_to_{to}','split':split,'model':'no_prior_gat','roc_auc':roc_auc_score(y,prob),'average_precision':average_precision_score(y,prob),'brier':brier_score_loss(y,prob),'n_train':len(tr),'n_target':len(y),'edge_count':ei.shape[1]})
 frame=pd.DataFrame(rows); frame.to_csv(a.output_dir/'origin_gat_metrics.csv',index=False); frame.groupby(['direction','model'])[['roc_auc','average_precision','brier']].agg(['mean','std']).reset_index().to_csv(a.output_dir/'origin_gat_summary.csv',index=False); (a.output_dir/'manifest.json').write_text(json.dumps({'origin_held_out':True,'no_global_prior':True,'shared_features':len(shared),'source_only_graph':True},indent=2))
if __name__=='__main__': main()
