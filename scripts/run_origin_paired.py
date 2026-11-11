#!/usr/bin/env python
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score,average_precision_score,brier_score_loss
from ibd_microbiome_framework.data.matrices import load_canonical_microbiome,load_shared_feature_list
from ibd_microbiome_framework.evaluation.splits import repeated_group_train_validation
from ibd_microbiome_framework.preprocessing.compositional import CLRFeaturePipeline
from scripts.run_no_prior_gat_transfer import edges,fit_gat,predict
def main():
 p=argparse.ArgumentParser();p.add_argument('--contract-dir',type=Path,required=True);p.add_argument('--legacy-root',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--n-splits',type=int,default=10);p.add_argument('--n-bootstrap',type=int,default=1000);a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True);shared=load_shared_feature_list(a.contract_dir);m=pd.read_csv(a.contract_dir/'cohort_manifest_prism.csv');x=load_canonical_microbiome(a.legacy_root/'prism_abundance_clean.csv',m,shared);m=m.set_index('sample_id',drop=False).loc[x.index].reset_index(drop=True);rows=[];pred=[];rng=np.random.default_rng(20261110)
 for so,to in [('discovery','validation'),('validation','discovery')]:
  sm=m[m.origin==so].reset_index(drop=True);tm=m[m.origin==to].reset_index(drop=True);source=x.loc[sm.sample_id];target=x.loc[tm.sample_id]
  for split,tr,va in repeated_group_train_validation(sm,a.n_splits,.2,20261110):
   ids=sm.iloc[tr].sample_id.tolist();pipe=CLRFeaturePipeline(top_k=min(100,len(shared)),random_state=20261110+split).fit(source.loc[ids],sm.iloc[tr].outcome_ibd);xs=pipe.transform(source.loc[ids]);xt=pipe.transform(target);ei=edges(xs.to_numpy(),'spearman',20261110+split);ytr=sm.iloc[tr].outcome_ibd.to_numpy();y=tm.outcome_ibd.to_numpy(); models={'no_prior_gat':predict(fit_gat(xs.to_numpy(),ytr,ei,20261110+split),xt.to_numpy(),ei),'elastic_net':LogisticRegression(penalty='elasticnet',solver='saga',l1_ratio=.5,C=1,class_weight='balanced',max_iter=3000,random_state=split).fit(xs,ytr).predict_proba(xt)[:,1],'random_forest':RandomForestClassifier(n_estimators=300,min_samples_leaf=2,class_weight='balanced',random_state=split,n_jobs=1).fit(xs,ytr).predict_proba(xt)[:,1]}
   direction=f'PRISM_{so}_to_{to}'
   for name,pr in models.items():
    rows.append({'direction':direction,'split':split,'model':name,'roc_auc':roc_auc_score(y,pr),'average_precision':average_precision_score(y,pr),'brier':brier_score_loss(y,pr)});pred.extend({'direction':direction,'split':split,'model':name,'sample_id':sid,'participant_id':str(pid),'outcome':int(yy),'probability':float(pp)} for sid,pid,yy,pp in zip(tm.sample_id,tm.participant_id,y,pr))
 d=pd.DataFrame(pred);d.to_csv(a.output_dir/'origin_paired_predictions.csv',index=False);pd.DataFrame(rows).to_csv(a.output_dir/'origin_paired_metrics.csv',index=False);out=[]
 for (direction,split),g in d.groupby(['direction','split']):
  ids=g.participant_id.unique();pan={z:g[g.model==z].set_index('sample_id') for z in g.model.unique()};
  for ref in ['elastic_net','random_forest']:
   vals=[];ga=pan['no_prior_gat'];gb=pan[ref]
   for _ in range(a.n_bootstrap):
    draw=rng.choice(ids,len(ids),replace=True);z=pd.concat([ga[ga.participant_id==q] for q in draw]);w=gb.loc[z.index];vals.append(roc_auc_score(z.outcome,z.probability)-roc_auc_score(w.outcome,w.probability))
   v=np.asarray(vals);out.append({'direction':direction,'split':split,'comparison':'gat_vs_'+ref,'delta_auc':v.mean(),'ci_low':np.quantile(v,.025),'ci_high':np.quantile(v,.975)})
 pd.DataFrame(out).to_csv(a.output_dir/'origin_paired_bootstrap.csv',index=False);(a.output_dir/'manifest.json').write_text(json.dumps({'origin_held_out':True,'paired_unit':'participant','n_splits':a.n_splits,'n_bootstrap':a.n_bootstrap},indent=2))
if __name__=='__main__':main()
