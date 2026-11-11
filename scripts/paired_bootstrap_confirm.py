#!/usr/bin/env python
"""Participant-cluster paired bootstrap for frozen target predictions."""
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss

def main():
    p=argparse.ArgumentParser(); p.add_argument('--predictions',type=Path,required=True); p.add_argument('--output-dir',type=Path,required=True); p.add_argument('--n-bootstrap',type=int,default=1000); p.add_argument('--seed',type=int,default=20261110); a=p.parse_args(); a.output_dir.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(a.predictions); rng=np.random.default_rng(a.seed); rows=[]
    for (direction,split),g in d.groupby(['direction','split']):
        ids=g.participant_id.astype(str).unique()
        panels={m:g[g.model==m].set_index('sample_id') for m in g.model.unique()}
        for ref in ['elastic_net','random_forest']:
            if 'no_prior_gat' not in panels or ref not in panels: continue
            a0,b0=panels['no_prior_gat'],panels[ref]; vals=[]
            for _ in range(a.n_bootstrap):
                draw=rng.choice(ids,len(ids),replace=True); pieces=[a0[a0.participant_id==q] for q in draw]; z=pd.concat(pieces); w=b0.loc[z.index]
                vals.append((roc_auc_score(z.outcome,z.probability)-roc_auc_score(w.outcome,w.probability),average_precision_score(z.outcome,z.probability)-average_precision_score(w.outcome,w.probability),brier_score_loss(z.outcome,z.probability)-brier_score_loss(w.outcome,w.probability)))
            v=np.asarray(vals); rows.append({'direction':direction,'split':split,'comparison':'no_prior_gat_vs_'+ref,'delta_auc':float(v[:,0].mean()),'delta_auc_ci_low':float(np.quantile(v[:,0],.025)),'delta_auc_ci_high':float(np.quantile(v[:,0],.975)),'delta_pr_auc':float(v[:,1].mean()),'delta_brier':float(v[:,2].mean()),'n_participants':len(ids),'n_bootstrap':a.n_bootstrap})
    frame=pd.DataFrame(rows); frame.to_csv(a.output_dir/'paired_bootstrap_by_split.csv',index=False); frame.groupby(['direction','comparison'])[['delta_auc','delta_auc_ci_low','delta_auc_ci_high','delta_pr_auc','delta_brier']].mean().reset_index().to_csv(a.output_dir/'paired_bootstrap_summary.csv',index=False); (a.output_dir/'manifest.json').write_text(json.dumps({'unit':'participant_cluster','n_bootstrap':a.n_bootstrap,'target_predictions_frozen':True},indent=2))
if __name__=='__main__': main()
