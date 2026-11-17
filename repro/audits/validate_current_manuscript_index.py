#!/usr/bin/env python3
"""Check paths and key scientific identity in the current manuscript index."""
import argparse
import csv
from pathlib import Path

import pandas as pd

EXPECTED = {f'Figure{i}' for i in (1, 2, 3, 4, 6, 7)} | {'Figure5a', 'Figure5b'} | {f'Table{i}' for i in range(1, 5)} | {f'SuppTable{i}' for i in range(1, 11)} | {'SuppFigureS1'}

def require(test, message):
    if not test:
        raise ValueError(message)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, default=Path('.'))
    parser.add_argument('--registry', type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    index = args.registry or root / 'manuscript/registry/figure_table_source_registry_20261117.csv'
    with index.open(encoding='utf8', newline='') as handle:
        rows = list(csv.DictReader(handle))
    ids = [r['display_id'] for r in rows]
    require(len(ids) == len(set(ids)) and set(ids) == EXPECTED, 'display IDs missing, duplicated, or obsolete')
    for row in rows:
        require(bool(row['interpretation']), f"missing interpretation for {row['display_id']}")
        for src in row['source_object'].split(';'):
            candidate = (root / src).resolve()
            require(candidate.is_relative_to(root) and candidate.is_file(), f'missing or unsafe source: {src}')
    s8 = next(r for r in rows if r['display_id'] == 'SuppTable8')
    require(s8['source_object'] == 'results/p0_rechecked_20261117/prism_direct_transfer_paired.csv', 'Table S8 points to a nonfinal summary')
    direct = pd.read_csv(root / s8['source_object'])
    require(set(direct.direction) == {'PRISM_discovery_to_validation', 'PRISM_validation_to_discovery'}, 'wrong direct-transfer directions')
    require(set(direct.reference) == {'elastic_net', 'random_forest', 'no_message_mlp'}, 'wrong reference models')
    require(len(direct) == 6 and (direct.n_bootstrap == 5000).all() and (direct.n_source_splits == 10).all(), 'incorrect transfer dimensions')
    require(set(direct.interval) == {'shared_target_draw_fixed_source_models'}, 'wrong bootstrap estimand')
    expected_s8 = {
        ('PRISM_discovery_to_validation', 'elastic_net'): ('-0.2166', '-0.3093', '-0.1163'),
        ('PRISM_discovery_to_validation', 'random_forest'): ('-0.2742', '-0.3806', '-0.1590'),
        ('PRISM_discovery_to_validation', 'no_message_mlp'): ('-0.1977', '-0.2889', '-0.0962'),
        ('PRISM_validation_to_discovery', 'elastic_net'): ('-0.1519', '-0.2072', '-0.0909'),
        ('PRISM_validation_to_discovery', 'random_forest'): ('-0.1869', '-0.2284', '-0.1405'),
        ('PRISM_validation_to_discovery', 'no_message_mlp'): ('-0.1550', '-0.1975', '-0.1087'),
    }
    for key, expected in expected_s8.items():
        selected = direct[(direct.direction == key[0]) & (direct.reference == key[1])]
        require(len(selected) == 1, f'Table S8 missing or duplicated {key}')
        row = selected.iloc[0]
        actual = tuple(f'{row[c]:.4f}' for c in ('point_mean_delta_auroc', 'ci_low', 'ci_high'))
        require(actual == expected, f'Table S8 manuscript numbers mismatch for {key}: {actual}')
    primary = pd.read_csv(root / 'results/p0_rechecked_20261117/primary_common_draw_summary.csv')
    xgboost = pd.read_csv(root / 'results/p0_rechecked_20261117/xgboost_paired_common_draw.csv')
    # Source tables preserve the historic field name `visit`. Its rows are
    # record-scored for HMP2; participant is the bootstrap sampling unit.
    metrics = {'roc_auc', 'average_precision', 'brier'}
    for label, data in [('primary', primary), ('xgboost', xgboost)]:
        for direction in ('HMP2_to_PRISM', 'PRISM_to_HMP2'):
            pooled = data[(data.direction == direction) & (data.target_stratum == 'pooled') & (data.analysis_level == 'visit')]
            require(set(pooled.metric) == metrics, f'Figure 4 {label} metrics missing for {direction}')
            require(len(pooled) == (9 if label == 'primary' else 3), f'Figure 4 {label} reference coverage mismatch for {direction}')
        for stratum in ('discovery', 'validation'):
            selection = data[(data.direction == 'HMP2_to_PRISM') & (data.target_stratum == stratum) & (data.analysis_level == 'visit') & (data.metric == 'roc_auc')]
            require(len(selection) == (3 if label == 'primary' else 1), f'Figure 5a {label} coverage mismatch for {stratum}')
        for level in ('visit', 'participant'):
            selection = data[(data.direction == 'PRISM_to_HMP2') & (data.target_stratum == 'pooled') & (data.analysis_level == level) & (data.metric == 'roc_auc')]
            require(len(selection) == (3 if label == 'primary' else 1), f'Figure 5b {label} coverage mismatch for {level}')
    row = primary[(primary.direction == 'HMP2_to_PRISM') & (primary.target_stratum == 'pooled') & (primary.analysis_level == 'visit') & (primary.reference == 'random_forest') & (primary.metric == 'roc_auc')]
    require(len(row) == 1 and abs(row.iloc[0].common_draw_bootstrap_mean - .0991) < .00006, 'primary pooled GAT–RF contrast mismatch')
    f3 = pd.read_csv(root / 'manuscript/source_data/Figure3_repeated_leakage.csv')
    require(len(f3) == 10 and abs(f3.random_split_subject_overlap.mean() - 43) < 1e-9, 'Figure 3 seed or overlap mismatch')
    print(f'PASS: {len(rows)} current display entries; paths, estimands, key contrasts and Figure 3 checked')
    print('LIMIT: summaries do not reproduce bootstrap draws or certify access to excluded individual predictions')

if __name__ == '__main__':
    main()
