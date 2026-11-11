#!/usr/bin/env python
"""Summarize topology density and split/hash stability from frozen GAT runs."""
import argparse, json
from pathlib import Path
import pandas as pd
def main():
 p=argparse.ArgumentParser(); p.add_argument('--runs-root',type=Path,required=True); p.add_argument('--output-dir',type=Path,required=True); a=p.parse_args(); a.output_dir.mkdir(parents=True,exist_ok=True); rows=[]
 for f in sorted(a.runs_root.glob('no_prior_gat_*/no_prior_gat_metrics.csv')):
  try: d=pd.read_csv(f)
  except Exception: continue
  for (topology,direction),g in d.groupby(['topology','direction']):
   rows.append({'run':f.parent.name,'topology':topology,'direction':direction,'n_splits':int(g.split.nunique()),'edge_count_mean':float(g.edge_count.mean()),'edge_count_sd':float(g.edge_count.std(ddof=1)),'n_unique_edge_hashes':int(g.edge_hash.nunique()),'edge_hash_all_unique':bool(g.edge_hash.nunique()==len(g))})
 frame=pd.DataFrame(rows); frame.to_csv(a.output_dir/'topology_run_audit.csv',index=False); frame.groupby(['topology','direction'])[['edge_count_mean','edge_count_sd','n_unique_edge_hashes']].mean().reset_index().to_csv(a.output_dir/'topology_stability_summary.csv',index=False); (a.output_dir/'manifest.json').write_text(json.dumps({'purpose':'topology_density_hash_audit','target_labels_not_used':True},indent=2))
if __name__=='__main__': main()
