#!/usr/bin/env python3
"""Check paths and key scientific identity in the current manuscript index."""
import argparse
import csv
import io
import itertools
import json
import zipfile
import subprocess
import sys
from pathlib import Path

import pandas as pd

EXPECTED = {f'Figure{i}' for i in (1, 2, 3, 4, 6)} | {'Figure5a', 'Figure5b', 'Figure5c', 'Figure5d'} | {f'Table{i}' for i in range(1, 5)} | {f'SuppTable{i}' for i in range(1, 16)} | {'SuppFigureS1', 'SuppFigureS2', 'SuppFigureS3'}

def require(test, message):
    if not test:
        raise ValueError(message)


def check_decision_table(root, objects):
    source = 'results/p0_rechecked_20261117/primary_common_draw_summary.csv'
    script = 'scripts/show_target_decisions.py'
    require(objects['Table4'] == [source, script], 'Table 4 must reference the common draws and display script')
    run = subprocess.run([sys.executable, str(root / script), '--root', str(root)],
                         check=True, capture_output=True, text=True)
    rows = list(csv.DictReader(run.stdout.splitlines()))
    expected = [
        ('Combined Franzosa', 'record', 'roc_auc', (.0993, .0476, .1491), 'GAT'),
        ('Combined Franzosa', 'record', 'brier', (-.0279, -.0381, -.0173), 'GAT'),
        ('PRISM discovery', 'record', 'roc_auc', (.1030, .0487, .1590), 'GAT'),
        ('PRISM discovery', 'record', 'brier', (-.0275, -.0382, -.0165), 'GAT'),
        ('Netherlands validation', 'record', 'roc_auc', (.0702, -.0289, .1736), 'unresolved'),
        ('Netherlands validation', 'record', 'brier', (-.0287, -.0511, -.0058), 'GAT'),
        ('HMP2', 'record', 'roc_auc', (-.0835, -.1501, -.0182), 'RF'),
        ('HMP2', 'record', 'brier', (-.0061, -.0320, .0174), 'unresolved'),
        ('HMP2', 'participant', 'roc_auc', (-.0915, -.1644, -.0182), 'RF'),
        ('HMP2', 'participant', 'brier', (.0086, -.0119, .0272), 'unresolved'),
    ]
    require(len(rows) == len(expected), 'Missing or extra Table 4 comparisons')
    for row, (target, unit, metric, triple, choice) in zip(rows, expected):
        key = row['target'], row['scoring_unit'], row['metric']
        numbers = tuple(round(float(row[k]), 4) for k in
                        ('observed_gat_minus_rf', 'conditional_95_low', 'conditional_95_high'))
        require(key == (target, unit, metric) and numbers == triple and row['supported_procedure'] == choice,
                f'Table 4 comparison no longer matches its manuscript entry for {target}/{unit}/{metric}')

def check_taxon_alignment(root, objects):
    archive = 'releases/formal_topology_edge_export_20261114.zip'
    summary_file = 'results/p0_rechecked_20261117/topology_taxon_aligned_summary.csv'
    require(archive in objects['Figure6'] and archive in objects['SuppTable5'], 'Missing formal edge lists')
    summary = pd.read_csv(root / summary_file).set_index('direction', verify_integrity=True)
    prefix = 'formal_topology_edge_export/'
    with zipfile.ZipFile(root / archive) as bundle:
        metadata = json.loads(bundle.read(prefix + 'formal_split_feature_manifest.json'))
        for direction in ['HMP2_to_PRISM', 'PRISM_to_HMP2']:
            fits = sorted((r for r in metadata if r['direction'] == direction), key=lambda r: r['split'])
            require([r['split'] for r in fits] == list(range(1, 11)), 'Incomplete feature order')
            mapped_edges = []
            for r in fits:
                names = r['selected_feature_ids']
                filename = f"{prefix}edge_lists/{direction}/split_{r['split']:02d}/real_spearman.csv"
                with bundle.open(filename) as source:
                    directed = {(int(e['source']), int(e['target']))
                                for e in csv.DictReader(io.TextIOWrapper(source, encoding='utf-8'))}
                require(len(directed) == 8 * len(names) and all(i != j for i, j in directed)
                        and len(names) == len(set(names)), 'Invalid directed Spearman input')
                mapped_edges.append((directed, {(names[i], names[j]) for i, j in directed}))
            for ix, field in enumerate(['local_mean', 'taxon_mean']):
                values = [len(a[ix] & b[ix]) / len(a[ix] | b[ix])
                          for a, b in itertools.combinations(mapped_edges, 2)]
                row = summary.loc[direction]
                require(len(values) == int(row.n_pairs) == 45 and
                        abs(sum(values)/len(values) - row[field]) < 1e-12,
                        f'Edge pairs disagree with {field} in {direction}')
                if field == 'taxon_mean':
                    require(abs(min(values) - row.taxon_min) < 1e-12 and
                            abs(max(values) - row.taxon_max) < 1e-12,
                            'Taxon-aligned overlap range mismatch')


def check_record_scores(root, objects, primary, xgboost):
    filename = 'manuscript/source_data/record_scores_by_fit.csv'
    for display in ['Table3', 'SuppTable6', 'SuppFigureS2']:
        require(filename in objects[display], f'{display} must reference the archived record scores')
    scores = pd.read_csv(root / filename)
    models = {'no_prior_gat','random_forest','elastic_net','capacity_matched_no_message_passing','xgboost'}
    require(len(scores) == 100 and set(scores.model) == models and
            not scores.duplicated(['direction','model','split']).any(), 'Wrong record-score file')
    combined = pd.concat([primary, xgboost])
    metrics = ['roc_auc', 'average_precision', 'brier']
    require(scores[metrics].apply(lambda col: col.between(0, 1).all()).all(), 'Invalid recorded scores')
    for direction, n_records, n_people in [('HMP2_to_PRISM',220,220),('PRISM_to_HMP2',1626,130)]:
        subset = scores[scores.direction == direction]
        require(subset.n_target_records.eq(n_records).all() and subset.n_target_participants.eq(n_people).all()
                and all(set(g.split) == set(range(1,11)) for _,g in subset.groupby('model')),
                'Missing source fits or wrong panel')
        selector = ((combined.direction == direction) & (combined.target_stratum == 'pooled') &
                    (combined.analysis_level == 'visit'))
        common = combined[selector].set_index(['reference','metric'], verify_integrity=True)
        absolute = subset.groupby('model')[metrics].mean()
        for name in models - {'no_prior_gat'}:
            for metric in metrics:
                difference = absolute.loc['no_prior_gat',metric] - absolute.loc[name,metric]
                require(abs(difference - common.loc[(name,metric),'split_mean_point_delta']) < 1e-12,
                        f'Incorrect pooled record score for {direction} / {name} / {metric}')


def check_selector_controls(data):
    sums = {'null_control': 30, 'batch_only': 27, 'signal_only': 20,
            'signal_reversing_batch': 6, 'signal_topology_shift': 20}
    require(len(data) == 50 and set(data.scenario) == set(sums), 'S1 scenario coverage mismatch')
    require(data.candidate_variables.eq(40).all() and data.selected_variables.eq(20).all()
            and data.designated_batch_variables.eq(6).all(), 'S1 candidate set or selector budget mismatch')
    for scenario, rows in data.groupby('scenario'):
        require(set(rows.seed) == set(range(20260815, 20260825)) and len(rows) == 10,
                f'S1 seed coverage mismatch for {scenario}')
        planted = 6 if scenario.startswith('signal_') else 0
        require(rows.implanted_disease_variables.eq(planted).all(), f'S1 disease truth mismatch: {scenario}')
        if planted:
            require(rows.recovered_disease_variables.eq(6).all() and rows.disease_recall.eq(1).all()
                    and rows.expected_disease_hits_uniform_selection.eq(3).all(), 'S1 positive truth mismatch')
        else:
            require(rows[['recovered_disease_variables', 'disease_recall', 'expected_disease_hits_uniform_selection']]
                    .isna().all().all(), 'S1 no-signal disease quantities must be undefined')
        require(rows.selected_batch_variables.isin(range(7)).all()
                and rows.selected_batch_variables.sum() == sums[scenario]
                and (rows.selected_batch_variables / 6 - rows.batch_candidate_fraction).abs().lt(1e-12).all(),
                f'S1 batch selection mismatch: {scenario}')


def check_brier_and_sensitivity(absolute, primary, sensitivity):
    reverse = primary.query("direction == 'PRISM_to_HMP2' and target_stratum == 'pooled'"
                            " and candidate == 'no_prior_gat' and reference == 'random_forest'")
    brier = reverse.query("metric == 'brier'").set_index('analysis_level', verify_integrity=True)
    reported = {'visit': (-.0061, -.0059, -.0320, .0174), 'participant': (.0086, .0085, -.0119, .0272)}
    require(set(brier.index) == set(reported), 'S12 Brier analysis levels mismatch')
    for key, expected in reported.items():
        s12 = absolute[(absolute.target == 'HMP2') & (absolute.scoring_unit == ('record' if key == 'visit' else key))]
        means = s12.set_index('model', verify_integrity=True).brier_mean
        row = brier.loc[key]
        cols = ['split_mean_point_delta', 'common_draw_bootstrap_mean', 'ci_low', 'ci_high']
        require(tuple(round(float(row[c]), 4) for c in cols) == expected and row.n_bootstrap == 5000
                and row.n_source_splits == 10 and row.n_target_participants == 130, 'S12 paired Brier mismatch')
        require(abs(means['Graph attention'] - means['Random forest'] - row.split_mean_point_delta) < 1e-12,
                'S12 observed Brier contrast and absolute scores disagree')
    panels = {'all_eligible': (1626, 1584, 1199), 'published_1595': (1595, 1553, 1169)}
    scores = ('record', 'person_profile_mean', 'person_collection_mean')
    metrics = ('roc_auc', 'average_precision', 'brier')
    s13 = sensitivity.set_index(['panel', 'scoring', 'metric'], verify_integrity=True)
    require(set(s13.index) == {(p,s,m) for p in panels for s in scores for m in metrics},
            'S13 must cover 2 panels, 3 scoring rules and all 3 metrics')
    require(s13.people.eq(130).all() and s13.ibd_people.eq(103).all() and s13.fits.eq(10).all()
            and s13.draws.eq(5000).all(), 'S13 participants, source fits or draw count mismatch')
    require((s13.gat - s13.rf - s13.observed_delta).abs().lt(1e-12).all()
            and s13.gat.between(0,1).all() and s13.rf.between(0,1).all()
            and (s13.ci_lower <= s13.draw_mean).all() and (s13.draw_mean <= s13.ci_upper).all(),
            'S13 absolute scores, contrasts or intervals disagree')
    for panel, counts in panels.items():
        subset = s13.loc[panel]
        require(subset[['profiles','collections','ibd_profiles']].eq(counts).all().all(),
                f'S13 record/collection counts incorrect for {panel}')
    for rule, analysis_level in [('record','visit'),('person_profile_mean','participant')]:
        for metric in metrics:
            source = reverse[(reverse.analysis_level == analysis_level) & (reverse.metric == metric)]
            require(len(source) == 1 and abs(s13.loc[('all_eligible', rule, metric), 'observed_delta']
                    - source.iloc[0].split_mean_point_delta) < 1e-12, 'S13 anchor and primary observed scores differ')
    require(s13.xs('roc_auc',level='metric').ci_upper.lt(0).all(), 'S13 AUROC interval changed sign')
    r = s13.xs('brier',level='metric')
    require(r.ci_lower.lt(0).all() and r.ci_upper.gt(0).all(), 'S13 Brier intervals must include zero')

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
            require(len(selection) == (3 if label == 'primary' else 1), f'Supplementary Figure 3a {label} coverage mismatch for {stratum}')
        for level in ('visit', 'participant'):
            selection = data[(data.direction == 'PRISM_to_HMP2') & (data.target_stratum == 'pooled') & (data.analysis_level == level) & (data.metric == 'roc_auc')]
            require(len(selection) == (3 if label == 'primary' else 1), f'Supplementary Figure 3b {label} coverage mismatch for {level}')
    row = primary[(primary.direction == 'HMP2_to_PRISM') & (primary.target_stratum == 'pooled') & (primary.analysis_level == 'visit') & (primary.reference == 'random_forest') & (primary.metric == 'roc_auc')]
    require(len(row) == 1 and abs(row.iloc[0].common_draw_bootstrap_mean - .0991) < .00006, 'primary pooled GAT–RF contrast mismatch')
    f3 = pd.read_csv(root / 'manuscript/source_data/Figure3_repeated_leakage.csv')
    require(len(f3) == 10 and abs(f3.random_split_subject_overlap.mean() - 43) < 1e-9, 'Figure 3 seed or overlap mismatch')
    objects = {r['display_id']: r['source_object'].split(';') for r in rows}
    check_decision_table(root, objects)
    check_taxon_alignment(root, objects)
    check_record_scores(root, objects, primary, xgboost)
    for display, filename in [('SuppTable14','pooled_auroc_components.csv'),
                              ('SuppTable15','participant_brier_components.csv')]:
        require(objects[display] == ['manuscript/source_data/' + filename,
                'results/p0_rechecked_20261117/primary_common_draw_summary.csv',
                'scripts/decompose_fixed_predictions.py'], f'{display} source mapping mismatch')
    subprocess.run([sys.executable, str(root / 'scripts/decompose_fixed_predictions.py'), '--check',
                    '--table-dir', str(root / 'manuscript/source_data'), '--summary',
                    str(root / 'results/p0_rechecked_20261117/primary_common_draw_summary.csv')], check=True)
    paired_s9 = 'results/reconciled_evidence/identity_aware_topology/paired_split_contrasts.csv'
    require(paired_s9 in objects['SuppTable9'], 'Missing stored taxon-embedding fit contrasts for Table S9')
    pairs = pd.read_csv(root / paired_s9)
    require(len(pairs) == 20 and not pairs.duplicated(['direction', 'split']).any()
            and set(pairs.direction) == {'HMP2_to_PRISM', 'PRISM_to_HMP2'}
            and all(set(group.split) == set(range(1, 11)) for _, group in pairs.groupby('direction'))
            and (pairs.real_auroc - pairs.rewired_auroc - pairs.real_minus_rewired).abs().lt(1e-12).all(),
            'Incorrect Table S9 identity-aware paired source fits')
    s12_sources = ['manuscript/source_data/figure5_model_comparison.csv', 'manuscript/source_data/figure5_gat_rf_by_fit.csv']
    common_draws = 'results/p0_rechecked_20261117/primary_common_draw_summary.csv'
    plot = 'scripts/plot_target_decisions.py'
    figure = 'manuscript/figures/Figure5_target_decisions.png'
    for panel, source in [('a', common_draws), ('b', common_draws),
                          ('c', 'manuscript/source_data/pooled_auroc_components.csv'),
                          ('d', 'manuscript/source_data/participant_brier_components.csv')]:
        require(objects['Figure5' + panel] == [source, plot, figure],
                f'Figure 5{panel} must reference its current source and rendering script')
    require(objects['SuppFigureS3'] == [s12_sources[0]],
            'Supplementary Figure 3 must reference the stratified four-model AUROC comparisons')
    require(objects['SuppTable12'] == s12_sources + [common_draws],
            'S12 requires the stratified AUROC sources and the primary paired Brier summary')
    selector_source = 'manuscript/source_data/supplementary_figure1_selector.csv'
    sensitivity_source = 'manuscript/source_data/hmp2_collection_sensitivity.csv'
    require(objects['SuppFigureS1'] == [selector_source], 'S1 must reference the 40-variable selector experiment')
    require(objects['SuppTable13'] == [sensitivity_source], 'S13 must reference the collection sensitivity analysis')
    check_selector_controls(pd.read_csv(root / selector_source))
    check_brier_and_sensitivity(pd.read_csv(root / s12_sources[0]), primary, pd.read_csv(root / sensitivity_source))
    subprocess.run([sys.executable, str(root / 'scripts/figure5_model_comparison.py'),
                    '--root', str(root), '--check'], check=True)
    summary = pd.read_csv(root / s12_sources[0]).query("model == 'Random forest'")
    fits = pd.read_csv(root / s12_sources[1])
    keys = ['source', 'target', 'scoring_unit']
    require(len(fits) == 40 and not fits.duplicated(keys + ['source_fit']).any(), 'Invalid source-fit table dimensions')
    require(((fits.gat_auroc - fits.rf_auroc - fits.gat_minus_rf_observed).abs() < 1e-12).all(), 'Incorrect fit differences')
    for row in summary.itertuples():
        group = fits[(fits.source == row.source) & (fits.target == row.target) & (fits.scoring_unit == row.scoring_unit)]
        require(set(group.source_fit) == set(range(1, 11)), 'Matched fits 1 through 10 required')
        require(abs(group.gat_minus_rf_observed.mean() - row.gat_minus_model_observed) < 1e-12, 'Fit and summary contrasts disagree')
    s11 = 'manuscript/source_data/siamcat_observed_reference.csv'
    require(objects['SuppTable11'] == [s11], 'S11 must reference the native observed reference')
    reference = pd.read_csv(root / s11)
    expected = {
        ('HMP2','Franzosa'): (220,220,.672147,.854972,.345883),
        ('Franzosa','HMP2'): (385,105,.736230,.872359,.198338),
        ('HMP2','NielsenHB_2014'): (396,318,.667262,.610132,.351597),
        ('Franzosa','NielsenHB_2014'): (396,318,.608827,.522080,.342822),
    }
    require(len(reference) == 4 and set(zip(reference.source, reference.target)) == set(expected), 'Wrong S11 panels')
    for row in reference.itertuples():
        actual = (row.n_target_records, row.n_target_participants, row.auroc_observed_mean,
                  row.average_precision_observed_mean, row.brier_observed_mean)
        require(actual == expected[(row.source,row.target)] and row.n_source_ensembles == 10 and row.scoring_unit == 'record', 'S11 observed summaries or panel mismatch')
    print(f'PASS: {len(rows)} current display entries; Table 4 decisions, pooled fit scores, formal edge overlaps, selector truth, S12–S15 checked')
    print('Intervals refer to shared target-participant draws for the ten fixed source fits.')

if __name__ == '__main__':
    main()
