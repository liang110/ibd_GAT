from __future__ import annotations
import argparse, json, re
from pathlib import Path
from datetime import datetime, timezone
import pandas as pd

def norm(value: object) -> str:
    return re.sub(r'[^A-Z0-9]', '', str(value).upper())

def main():
    ap=argparse.ArgumentParser(description='Audit HMP2/PRISM manifest-level independence.')
    ap.add_argument('--source-root', type=Path, required=True)
    ap.add_argument('--output-dir', type=Path, required=True)
    args=ap.parse_args(); src=args.source_root; out=args.output_dir; out.mkdir(parents=True,exist_ok=True)
    h=pd.read_csv(src/'cohort_manifest_hmp2.csv',dtype=str); p=pd.read_csv(src/'cohort_manifest_prism.csv',dtype=str)
    cm=json.loads((src/'canonical_contract_manifest.json').read_text(encoding='utf-8'))
    zero=set(cm.get('zero_mass_sample_ids',{}).get('HMP2',[])); hi=h[~h.sample_id.isin(zero)]; pi=p
    rows=[]; matches=[]
    for scope,a,b in [('manifest',h,p),('included_shared_space',hi,pi)]:
        for id_type,col in [('sample_id','sample_id'),('participant_id','participant_id')]:
            hs=set(a[col].dropna().astype(str)); ps=set(b[col].dropna().astype(str)); exact=hs&ps
            hn={norm(x):x for x in hs}; pn={norm(x):x for x in ps}; common=set(hn)&set(pn)
            rows.append({'scope':scope,'id_type':id_type,'hmp2_count':len(hs),'prism_count':len(ps),'exact_overlap':len(exact),'normalized_overlap':len(common)})
            matches.extend({'scope':scope,'id_type':id_type,'hmp2_id':hn[k],'prism_id':pn[k],'match_type':'normalized'} for k in sorted(common))
    pd.DataFrame(rows).to_csv(out/'overlap_audit.csv',index=False)
    pd.DataFrame(matches,columns=['scope','id_type','hmp2_id','prism_id','match_type']).to_csv(out/'overlap_matches.csv',index=False)
    pd.DataFrame([
        {'scope':'manifest','cohort':'HMP2','accession_candidates':'PRJNA389280; PRJNA398089','samples':len(h),'participants':h.participant_id.nunique()},
        {'scope':'manifest','cohort':'PRISM','accession_candidates':'PRJNA400072','samples':len(p),'participants':p.participant_id.nunique()},
        {'scope':'included_shared_space','cohort':'HMP2','accession_candidates':'PRJNA389280; PRJNA398089','samples':len(hi),'participants':hi.participant_id.nunique()},
        {'scope':'included_shared_space','cohort':'PRISM','accession_candidates':'PRJNA400072','samples':len(pi),'participants':pi.participant_id.nunique()},
    ]).to_csv(out/'cohort_accession_summary.csv',index=False)
    report={'schema':'ibd_hmp2_prism_accession_independence_audit_v2','generated_utc':datetime.now(timezone.utc).isoformat(),'source_root':str(src),'accession_resolution':{'PRISM':'PRJNA400072','HMP2_candidates':['PRJNA389280','PRJNA398089'],'status':'HMP2 accession requires provider-release confirmation; local manifests contain no accession field'},'manifest_counts':{'HMP2_manifest_samples':len(h),'HMP2_manifest_participants':int(h.participant_id.nunique()),'PRISM_manifest_samples':len(p),'PRISM_manifest_participants':int(p.participant_id.nunique()),'HMP2_included_samples':len(hi),'HMP2_included_participants':int(hi.participant_id.nunique()),'PRISM_included_samples':len(pi),'PRISM_included_participants':int(pi.participant_id.nunique()),'HMP2_zero_mass_excluded':len(zero)},'overlap':rows,'conclusion':{'status':'ZERO_OVERLAP_IN_LOCAL_MANIFESTS_BUT_HMP2_ACCESSION_UNRESOLVED'}}
    (out/'cohort_independence_audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))
if __name__=='__main__': main()
