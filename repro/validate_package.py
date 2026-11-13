from pathlib import Path
import csv

ROOT = Path(__file__).resolve().parents[1]
REQUIRED = [
    ROOT / 'results' / 'submission_tables' / 'Table_1_cohort_contract.csv',
    ROOT / 'results' / 'submission_tables' / 'Table_2_primary_benchmark_metrics.csv',
    ROOT / 'results' / 'submission_tables' / 'Supplementary_Table_1_known_truth_controls.csv',
    ROOT / 'manuscript' / 'registry' / 'figure_table_source_registry.csv',
    ROOT / 'results' / 'reconciled_evidence' / 'synthetic_truth_summary.csv',
    ROOT / 'results' / 'reconciled_evidence' / 'absolute_performance_summary.csv',
    ROOT / 'results' / 'reconciled_evidence' / 'paired_bootstrap_reconciled_by_split.csv',
]

def main() -> int:
    missing = [str(p.relative_to(ROOT)) for p in REQUIRED if not p.is_file()]
    if missing:
        for item in missing:
            print(f'MISSING {item}')
        return 1
    for path in REQUIRED:
        with path.open(newline='', encoding='utf-8-sig') as handle:
            rows = list(csv.reader(handle))
        if len(rows) < 2:
            print(f'EMPTY {path.relative_to(ROOT)}')
            return 1
        print(f'OK {path.relative_to(ROOT)} ({len(rows)-1} data rows)')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
