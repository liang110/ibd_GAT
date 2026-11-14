#!/usr/bin/env python3
"""Recompute known-truth acceptance rates from per-seed records."""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

def parse_bool(value):
    if isinstance(value, bool): return value
    return str(value).strip().lower() in {"true", "1", "yes"}

def main() -> int:
    p=argparse.ArgumentParser(); p.add_argument('--ledger',type=Path,required=True); a=p.parse_args(); d=pd.read_csv(a.ledger)
    required={'experiment','scenario','seed','target_roc_auc','n_disease_features_recovered','acceptance_rule','acceptance_pass'}; miss=required-set(d.columns)
    if miss: raise SystemExit(f'MISSING columns: {sorted(miss)}')
    expected=[]
    for _,r in d.iterrows():
        rule=str(r.acceptance_rule)
        if rule=='positive_auc_0.70_and_recovery_3_of_6': ok=(float(r.target_roc_auc)>=0.70) and (int(r.n_disease_features_recovered)>=3)
        elif rule=='negative_auc_0.35_to_0.65': ok=(0.35<=float(r.target_roc_auc)<=0.65)
        else: raise SystemExit(f'UNKNOWN acceptance rule: {rule}')
        expected.append(ok)
    actual=d.acceptance_pass.map(parse_bool).tolist()
    if expected!=actual:
        mask=[x!=y for x,y in zip(expected,actual)]
        print(d.loc[mask,['experiment','scenario','seed','acceptance_rule','acceptance_pass']].to_string(index=False)); return 2
    if 'truth_denominator' in d.columns:
        positive=d.acceptance_rule.eq('positive_auc_0.70_and_recovery_3_of_6')
        denom=pd.to_numeric(d.loc[positive,'truth_denominator'], errors='coerce')
        if not denom.eq(6).all(): print('FAIL: positive truth denominator is not 6'); return 2
    out=d.groupby(['experiment','scenario','acceptance_rule'],as_index=False).agg(n_seeds=('seed','size'),n_pass=('acceptance_pass','sum'),mean_auc=('target_roc_auc','mean'),mean_recovered=('n_disease_features_recovered','mean')); out['acceptance_rate']=out.n_pass/out.n_seeds; print(out.to_string(index=False)); return 0
if __name__=='__main__': raise SystemExit(main())
