#!/usr/bin/env python
import argparse,json
from pathlib import Path
import pandas as pd
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_auc_score,average_precision_score,brier_score_loss
from ibd_microbiome_framework.data.matrices import load_canonical_microbiome,load_shared_feature_list
from ibd_microbiome_framework.evaluation.splits import repeated_group_train_validation
from ibd_microbiome_framework.preprocessing.compositional import CLRFeaturePipeline
def main():
 p=argparse.ArgumentParser(); p.add_argument('--contract-dir',type=Path,required=True);p.add_argument('--legacy-root',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--n-splits',type=int,default=10);a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True); shared=load_shared_feature_list(a.contract_dir); hm=pd.read_csv(a.contract_dir/'cohort_manifest_hmp2.csv');pm=pd.read_csv(a.contract_dir/'cohort_manifest_prism.csv');h=load_canonical_microbiome(a.legacy_root/'manual_abundance.csv',hm,shared);q=load_canonical_microbiome(a.legacy_root/'prism_abundance_clean.csv',pm,shared);hm=hm.set_index('sample_id',drop=False).loc[h.index].reset_index(drop=True);pm=pm.set_index('sample_id',drop=False).loc[q.index].reset_index(drop=True);rows=[]
 for direction,x,sm,t,tm in [('HMP2_to_PRISM',h,hm,q,pm),('PRISM_to_HMP2',q,pm,h,hm)]:
  for split,tr,va in repeated_group_train_validation(sm,a.n_splits,.2,20261110):
   ids=sm.iloc[tr].sample_id.tolist();pipe=CLRFeaturePipeline(top_k=min(100,len(shared)),random_state=20261110+split).fit(x.loc[ids],sm.iloc[tr].outcome_ibd);xs=pipe.transform(x.loc[ids]);xt=pipe.transform(t);m=MLPClassifier(hidden_layer_sizes=(32,16),alpha=1e-3,early_stopping=False,max_iter=400,random_state=20261110+split).fit(xs,sm.iloc[tr].outcome_ibd);pr=m.predict_proba(xt)[:,1];y=tm.outcome_ibd.to_numpy();rows.append({'direction':direction,'split':split,'model':'capacity_matched_no_message_passing','roc_auc':roc_auc_score(y,pr),'average_precision':average_precision_score(y,pr),'brier':brier_score_loss(y,pr),'n_features':xs.shape[1]})
 f=pd.DataFrame(rows);f.to_csv(a.output_dir/'capacity_control_metrics.csv',index=False);f.groupby('direction')[['roc_auc','average_precision','brier']].agg(['mean','std']).to_csv(a.output_dir/'capacity_control_summary.csv');(a.output_dir/'manifest.json').write_text(json.dumps({'no_message_passing':True,'source_only':True,'n_splits':a.n_splits},indent=2))
if __name__=='__main__':main()
