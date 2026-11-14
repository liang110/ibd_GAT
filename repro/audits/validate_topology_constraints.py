#!/usr/bin/env python3
"""Validate formal Gate B topology constraints and fail closed on mismatches."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import pandas as pd

def parse_bool(value):
    if isinstance(value, bool): return value
    return str(value).strip().lower() in {"true", "1", "yes"}

def main() -> int:
    p=argparse.ArgumentParser(); p.add_argument('--edge-audit',type=Path,required=True); p.add_argument('--hash-audit',type=Path,required=True); p.add_argument('--status',type=Path,default=None); p.add_argument('--allow-hash-failures',action='store_true'); p.add_argument('--allow-incomplete',action='store_true'); a=p.parse_args()
    edge=pd.read_csv(a.edge_audit); hashes=pd.read_csv(a.hash_audit)
    required={'direction','split','condition','self_loop_count','jaccard_with_identity','unique_edge_count','edge_count','n_nodes','out_degree_min','out_degree_max','out_degree_variance'}
    missing=required-set(edge.columns)
    if missing: raise SystemExit(f'MISSING edge-audit columns: {sorted(missing)}')
    failures=[]; non=edge[edge.condition.astype(str).ne('identity')]; ident=edge[edge.condition.astype(str).eq('identity')]
    if (non.self_loop_count.astype(int)!=0).any(): failures.append('non-identity self_loop_count != 0')
    if (non.unique_edge_count.astype(int)!=non.edge_count.astype(int)).any(): failures.append('non-identity duplicate edge detected')
    if (non.edge_count.astype(int)!=8*non.n_nodes.astype(int)).any(): failures.append('non-identity edge_count != 8*n_nodes')
    if (non.out_degree_min.astype(int)!=8).any() or (non.out_degree_max.astype(int)!=8).any() or (non.out_degree_variance.astype(float)!=0).any(): failures.append('non-identity fixed out-degree constraint failed')
    if (non.jaccard_with_identity.astype(float)!=0).any(): failures.append('non-identity identity-Jaccard != 0')
    if len(ident) and ('identity_graph_pass' in ident.columns) and not ident.identity_graph_pass.map(parse_bool).all(): failures.append('identity graph constraint failed')
    if 'hash_matches' in hashes.columns and not a.allow_hash_failures and not hashes.hash_matches.map(parse_bool).all(): failures.append('formal edge hash mismatch')
    if 'selected_feature_hash_matches' in hashes.columns and not a.allow_hash_failures and not hashes.selected_feature_hash_matches.map(parse_bool).all(): failures.append('selected-feature hash mismatch')
    status_path=a.status or (a.edge_audit.parent/'topology_audit_status.json')
    if status_path.exists():
        status=json.loads(status_path.read_text(encoding='utf-8')); formal=status.get('formal_gateb',{})
        if not formal.get('raw_edge_list_exports_present',False) and not a.allow_incomplete: failures.append('formal Gate B raw edge-list export missing')
        if int(formal.get('hash_mismatch_records',0)) and not a.allow_incomplete and not a.allow_hash_failures: failures.append('formal Gate B status reports unresolved hash mismatches')
    print(f'edge_records={len(edge)} hash_records={len(hashes)} nonidentity={len(non)} identity={len(ident)}')
    if failures: print('FAIL: '+'; '.join(failures)); return 2
    print('PASS: formal Gate B edge constraints and hash reconciliation validated'); return 0
if __name__=='__main__': raise SystemExit(main())
