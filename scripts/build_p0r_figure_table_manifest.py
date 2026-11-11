"""Freeze P0-R result artifacts and emit the BIB figure/table manifest.

The registry is intentionally small and deterministic: it records the exact
analysis artifacts used by the planned figures and tables, their row counts,
and SHA256 digests.  It does not modify the manuscript or regenerate plots.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


ARTIFACTS = [
    ("analysis_manifest", "p0r_analysis_manifest.json", None, "协议、输入和运行参数审计"),
    ("paired_model_deltas", "p0r_paired_model_deltas.csv", 270, "模型配对差值主结果"),
    ("paired_model_delta_draws", "p0r_paired_model_delta_draws.csv", 45000, "模型配对差值 bootstrap 抽样"),
    ("hmp2_participant_by_split", "p0r_hmp2_participant_level_by_split.csv", 360, "HMP2 参与者级敏感性分层"),
    ("hmp2_participant_summary", "p0r_hmp2_participant_level_summary.csv", 36, "HMP2 参与者级摘要"),
    ("hmp2_participant_draws", "p0r_hmp2_participant_level_draws.csv", 6000, "HMP2 参与者级 bootstrap 抽样"),
    ("calibration_by_split", "p0r_direction_calibration_by_split.csv", 360, "来源扰动校准分布"),
    ("calibration_summary", "p0r_direction_calibration_summary.csv", 144, "方向级校准摘要"),
    ("reliability_bins", "p0r_direction_calibration_reliability_bins.csv", 331, "方向级可靠性分箱"),
    ("direction_target_counts", "p0r_direction_target_counts.csv", 6, "六个迁移方向目标结构"),
    ("cohort_contract_summary", "p0r_cohort_contract_summary.csv", 3, "三队列数据合同摘要"),
]

AUDIT_ARTIFACTS = [
    ("source_only_protocol_audit", "../protocol_audit/p0r_source_only_protocol_audit.json", None, "来源端嵌套流程静态审计"),
    ("source_only_protocol_audit_report", "../protocol_audit/P0-R-C_来源端嵌套流程审计.md", None, "来源端嵌套流程审计报告"),
]

MANUSCRIPT_ARTIFACTS = [
    ("standard_model_metrics", "standard_model_benchmark_metrics.csv", 420, "六方向标准模型原始指标"),
    ("standard_model_predictions", "standard_model_benchmark_predictions.csv", 140140, "六方向标准模型逐样本预测"),
    ("standard_model_cluster_bootstrap", "standard_model_aggregate_cluster_bootstrap.csv", 42, "标准模型目标参与者聚簇 bootstrap"),
    ("siamcat_reference_metrics", "siamcat_reference_metrics.csv", 40, "SIAMCAT 独立 source-only 指标"),
    ("siamcat_cluster_bootstrap", "siamcat_aggregate_cluster_bootstrap.csv", 4, "SIAMCAT 独立目标聚簇 bootstrap"),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def row_count(path: Path) -> int | None:
    if path.suffix.lower() != ".csv":
        return None
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return max(sum(1 for _ in csv.reader(handle)) - 1, 0)


def build_yaml(rows: list[dict[str, object]], result_dir: Path) -> str:
    lines = [
        "schema_version: p0r-figure-table-manifest-v2",
        "status: frozen_for_manuscript_reconstruction",
        "scope: three_cohort_ibd_microbiome_transfer_benchmark",
        "manuscript_run_ids_in_body: forbidden",
        "primary_contract:",
        "  shared_species: 126",
        "  directions: 6",
        "  source_perturbations: 10",
        "  trainable_model_families: 6",
        "  target_uncertainty: participant_cluster_bootstrap",
        "  source_uncertainty: source_split_mean_and_sd",
        "  calibration: ensemble_target_predictions_then_one_cluster_bootstrap",
        "figure_numbering:",
        "  main:",
        "    - id: Figure_1",
        "      output_basename: Figure_1_benchmark_architecture",
        "      title: 审计约束下的跨队列 benchmark 架构",
        "      role: design_and_leakage_controls",
        "      source: source_only_protocol_audit + cohort_contract_summary",
        "      uncertainty: source_split_definition_and_target_cluster_definition",
        "      claim_boundary: design_and_audit_only",
        "      panels: [contract, source_only_boundary, six_directions, uncertainty_layers]",
        "    - id: Figure_2",
        "      output_basename: Figure_2_cohort_contract_and_controls",
        "      title: 三队列样本与参与者结构及可迁移性审计",
        "      role: cohort_and_known_truth_controls",
        "      source: cohort_contract_summary + direction_target_counts + known_truth_and_repeated_subject_controls",
        "      uncertainty: descriptive_counts",
        "      claim_boundary: cohort_scope_and_structure_only",
        "      panels: [cohort_contract, target_counts, synthetic_controls, repeated_subject_controls]",
        "    - id: Figure_3",
        "      output_basename: Figure_3_standard_model_transfer_benchmark",
        "      title: 标准模型家族的六方向外部迁移 benchmark",
        "      role: primary_benchmark",
        "      source: frozen_primary_metric_table + standard_model_predictions",
        "      uncertainty: target_participant_cluster_bootstrap_CI_with_source_split_mean_SD",
        "      claim_boundary: direction_and_model_dependence_only",
        "      panels: [auroc, pr_auc, brier, direction_by_model_heatmap]",
        "    - id: Figure_4",
        "      output_basename: Figure_4_paired_deltas_uncertainty",
        "      title: 模型配对差值与两层不确定性",
        "      role: paired_model_differences",
        "      source: paired_model_deltas + paired_model_delta_draws",
        "      uncertainty: paired_target_cluster_bootstrap_and_source_split_variability",
        "      claim_boundary: comparative_descriptive_only",
        "      panels: [delta_forest, source_split_variability, target_cluster_ci]",
        "    - id: Figure_5",
        "      output_basename: Figure_5_hmp2_participant_calibration",
        "      title: HMP2 参与者级敏感性与方向级概率校准",
        "      role: p0r_sensitivity_and_calibration",
        "      source: hmp2_participant_summary + calibration_summary + reliability_bins",
        "      uncertainty: participant_cluster_bootstrap_on_frozen_target_ensemble",
        "      claim_boundary: sensitivity_and_external_calibration_description_only",
        "      panels: [sample_vs_participant, calibration_slope_intercept, reliability_bins]",
        "    - id: Figure_6",
        "      output_basename: Figure_6_native_siamcat_reference",
        "      title: 原生 SIAMCAT source-only 外部标准参考",
        "      role: separate_external_reference",
        "      source: siamcat_reference_metrics + siamcat_cluster_bootstrap",
        "      uncertainty: SIAMCAT_protocol_specific_target_cluster_bootstrap",
        "      claim_boundary: external_reference_not_combined_leaderboard",
        "      panels: [protocol_boundary, direction_metrics, comparison_scope]",
        "  supplementary:",
        "    - id: Supplementary_Figure_1",
        "      title: 来源端嵌套流程与泄漏防护审计",
        "      role: protocol_audit",
        "    - id: Supplementary_Figure_2",
        "      title: PRISM 来源组成与公开元数据边界",
        "      role: cohort_scope_audit",
        "    - id: Supplementary_Figure_3",
        "      title: 亚型范围审计与不作机制外推的依据",
        "      role: scope_audit",
        "    - id: Supplementary_Figure_4",
        "      title: 学习曲线、稳定性和特征/权重敏感性",
        "      role: robustness_secondary",
        "table_numbering:",
        "  main:",
        "    - id: Table_1",
        "      title: 三队列数据合同与目标端有效样本结构",
        "      source: cohort_contract_summary + direction_target_counts",
        "    - id: Table_2",
        "      title: 六方向标准模型 benchmark 指标与两层不确定性",
        "      source: paired_model_deltas + frozen primary metric table",
        "    - id: Table_3",
        "      title: 原生 SIAMCAT source-only 参考结果",
        "      source: separate SIAMCAT source data; never mixed with CLR leaderboard",
        "    - id: Table_4",
        "      title: P0-R HMP2 参与者级敏感性与方向级校准摘要",
        "      source: hmp2_participant_summary + calibration_summary",
        "  supplementary:",
        "    - id: Supplementary_Table_1",
        "      title: 已知真值合成控制",
        "    - id: Supplementary_Table_2",
        "      title: 重复参与者控制与聚簇结构",
        "    - id: Supplementary_Table_3",
        "      title: PRISM 来源组成审计",
        "    - id: Supplementary_Table_4",
        "      title: 亚型范围审计",
        "    - id: Supplementary_Table_5",
        "      title: 特征稳定性",
        "    - id: Supplementary_Table_6",
        "      title: 权重稳定性",
        "    - id: Supplementary_Table_7",
        "      title: 调参稳定性",
        "source_artifacts:",
    ]
    for row in rows:
        lines.extend(
            [
                f"  - key: {row['key']}",
                f"    path: {row['path']}",
                f"    sha256: {row['sha256']}",
                f"    rows: {row['rows'] if row['rows'] is not None else 'na'}",
                f"    role: {row['role']}",
            ]
        )
    lines.extend(
        [
            "display_rules:",
            "  - 主图优先报告 AUROC、PR-AUC、Brier；每个点明确指标和方向。",
            "  - source-split mean/sd 与 target participant-cluster bootstrap CI 必须分栏或分图层呈现。",
            "  - HMP2 主敏感性使用参与者级摘要，同时保留样本级结果作为补充观察。",
            "  - 校准仅作方向级描述，不宣称目标适配、临床风险校准或部署效用。",
            "  - SIAMCAT 作为领域标准参考单列，不能与 CLR 主榜合并排名。",
            "  - 图注和表注不得出现内部 run ID；内部路径与哈希仅保存在本清单和注册表。",
            "  - 不保留图模型、多组学和目标端适配作为核心结果图表。",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--result-dir",
        type=Path,
        default=Path("ibd_framework_scaffold/artifacts/runs/p0r_closure_20260818/analysis"),
    )
    parser.add_argument(
        "--manuscript-dir",
        type=Path,
        default=Path("ibd_framework_scaffold/manuscript"),
    )
    args = parser.parse_args()
    args.manuscript_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    for key, filename, expected, role in ARTIFACTS + AUDIT_ARTIFACTS:
        path = args.result_dir / filename
        if not path.exists():
            raise FileNotFoundError(path)
        observed = row_count(path)
        if expected is not None and observed != expected:
            raise ValueError(f"{filename}: expected {expected} rows, observed {observed}")
        rows.append(
            {
                "key": key,
                "path": str(path).replace("\\", "/"),
                "sha256": sha256(path),
                "rows": observed,
                "role": role,
            }
        )

    for key, filename, expected, role in MANUSCRIPT_ARTIFACTS:
        path = args.manuscript_dir / "source_data" / filename
        if not path.exists():
            raise FileNotFoundError(path)
        observed = row_count(path)
        if expected is not None and observed != expected:
            raise ValueError(f"{filename}: expected {expected} rows, observed {observed}")
        rows.append(
            {
                "key": key,
                "path": str(path).replace("\\", "/"),
                "sha256": sha256(path),
                "rows": observed,
                "role": role,
            }
        )

    registry = args.manuscript_dir / "frozen_result_registry_p0r_v2.csv"
    with registry.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["key", "path", "sha256", "rows", "role"])
        writer.writeheader()
        writer.writerows(rows)

    manifest = args.manuscript_dir / "figure_table_manifest_p0r_v2.yaml"
    manifest.write_text(build_yaml(rows, args.result_dir), encoding="utf-8")
    print(f"wrote {registry}")
    print(f"wrote {manifest}")
    print(f"validated {len(rows)} immutable manuscript result and audit artifacts")


if __name__ == "__main__":
    main()
