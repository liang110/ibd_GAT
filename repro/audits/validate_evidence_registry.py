#!/usr/bin/env python3
"""Check that compact figure/table source files listed in the registry exist."""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

def main() -> int:
    p=argparse.ArgumentParser(); p.add_argument("--root", type=Path, default=Path(".")); a=p.parse_args()
    reg=a.root/"results/reconciled_evidence/evidence_registry.csv"
    if not reg.exists(): raise SystemExit(f"MISSING registry: {reg}")
    df=pd.read_csv(reg); required={"object","primary_sources","status","limitation"}; miss=required-set(df.columns)
    if miss: raise SystemExit(f"MISSING columns: {sorted(miss)}")
    missing=[]
    for _, row in df.iterrows():
        for raw in str(row.primary_sources).split(";"):
            rel=raw.strip()
            if not rel: continue
            # sources are compact files under results/reconciled_evidence unless already rooted
            path=(a.root/"results/reconciled_evidence"/rel)
            if not path.exists(): missing.append(f"{row['object']}: {rel}")
    if missing:
        print("MISSING SOURCES:"); print("\n".join(missing)); return 2
    print(f"PASS: {len(df)} registered objects and all compact sources exist")
    return 0
if __name__=="__main__": raise SystemExit(main())
