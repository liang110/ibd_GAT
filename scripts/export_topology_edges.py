#!/usr/bin/env python
import argparse,json
from pathlib import Path
import pandas as pd,numpy as np
from ibd_microbiome_framework.data.matrices import load_canonical_microbiome,load_shared_feature_list
from ibd_microbiome_framework.evaluation.splits import repeated_group_train_validation
from ibd_microbiome_framework.preprocessing.compositional import CLRFeaturePipeline
from scripts.run_no_prior_gat_transfer import edges
def main():
 p=argparse.ArgumentParser();p.add_argument('--contract-dir',type=Path,required=True);p.add_argument('--legacy-root',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--n-splits',type=int,default=10);p.add_argument('--topologies',default='pearson,spearman,shuffled');a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True);shared=load_shared_feature_list(a.contract_dir);hm=pd.read_csv(a.contract_dir/'cohort_manifest_hmp2.csv');pm=pd.read_csv(a.contract_dir/'cohort_manifest_prism.csv');h=load_canonical_microbiome(a.legacy_root/'manual_abundance.csv',hm,shared);q=load_canonical_microbiome(a.legacy_root/'prism_abundance_clean.csv',pm,shared);hm=hm.set_index('sample_id',drop=False).loc[h.index].reset_index(drop=True);pm=pm.set_index('sample_id',drop=False).loc[q.index].reset_index(drop=True);rows=[]
 for direction,x,m,t in [('HMP2_to_PRISM',h,hm,q),('PRISM_to_HMP2',q,pm,h)]:
  for split,tr,va in repeated_group_train_validation(m,a.n_splits,.2,20261110):
   ids=m.iloc[tr].sample_id.tolist();pipe=CLRFeaturePipeline(top_k=min(100,len(shared)),random_state=20261110+split).fit(x.loc[ids],m.iloc[tr].outcome_ibd);z=pipe.transform(x.loc[ids]).to_numpy()
   for top in a.topologies.split(','):
    ei=edges(z,'pearson' if top=='pearson' else top,20261110+split); e=ei.numpy().T; pd.DataFrame(e,columns=['source','target']).to_csv(a.output_dir/f'{direction}_split{split}_{top}_edges.csv',index=False);deg=np.bincount(e[:,0],minlength=z.shape[1]) if len(e) else np.zeros(z.shape[1]); rows.append({'direction':direction,'split':split,'topology':top,'edge_count':len(e),'edge_density':len(e)/(z.shape[1]*(z.shape[1]-1)),'degree_mean':deg.mean(),'degree_sd':deg.std()})
 pd.DataFrame(rows).to_csv(a.output_dir/'topology_edge_manifest.csv',index=False);(a.output_dir/'manifest.json').write_text(json.dumps({'edge_sets_exported':True,'degree_vectors_exported':True,'n_splits':a.n_splits},indent=2))
if __name__=='__main__':main()
