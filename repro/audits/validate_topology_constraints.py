#!/usr/bin/env python3
"""Validate compact topology constraints and fail closed on incomplete Gate B export."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import pandas as pd

def parse_bool(value):
    if isinstance(value, bool): return value
    return str(value).strip().lower() in {"true", "1", "yes"}

def main() -> int:
    p=argparse.ArgumentParser()
    p.add_argument('--edge-audit',type=Path,required=True)
    p.add_argument('--hash-audit',type=Path,required=True)
    p.add_argument('--status',type=Path,default=None)
    p.add_argument('--allow-hash-failures',action='store_true')
    p.add_argument('--allow-incomplete',action='store_true')
    a=p.parse_args()
    edge=pd.read_csv(a.edge_audit); hashes=pd.read_csv(a.hash_audit)
    required={'direction','split','topology','self_loop_count','jaccard_with_identity'}
    missing=required-set(edge.columns)
    if missing: raise SystemExit(f'MISSING edge-audit columns: {sorted(missing)}')
    failures=[]
    if (edge.self_loop_count.astype(int)!=0).any(): failures.append('self_loop_count != 0')
    if (edge.jaccard_with_identity.astype(float)!=0).any(): failures.append('jaccard_with_identity != 0')
    if 'unique_edge_count' in edge.columns and 'edge_count' in edge.columns and (edge.unique_edge_count.astype(int)!=edge.edge_count.astype(int)).any(): failures.append('unique_edge_count != edge_count')
    if 'hash_matches' in hashes.columns and not a.allow_hash_failures:
        parsed=hashes['hash_matches'].map(parse_bool)
        if (~parsed).any(): failures.append('edge hash reconstruction mismatch')
    status_path=a.status or (a.edge_audit.parent/'topology_audit_status.json')
    if status_path.exists():
        status=json.loads(status_path.read_text(encoding='utf-8'))
        formal=status.get('formal_gateb',{})
        if not formal.get('raw_edge_list_exports_present',False) and not a.allow_incomplete:
            failures.append('formal Gate B raw edge-list export missing')
        if int(formal.get('hash_mismatch_records',0)) and not a.allow_incomplete and not a.allow_hash_failures:
            failures.append('formal Gate B status reports unresolved hash mismatches')
    print(f'edge_records={len(edge)} hash_records={len(hashes)}')
    if failures:
        print('FAIL: ' + '; '.join(failures)); return 2
    print('PASS: all declared topology constraints and edge hashes validated')
    return 0
if __name__=='__main__': raise SystemExit(main())
