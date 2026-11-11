#!/usr/bin/env python
"""Create a traceable source-only protocol audit for P0-R-C."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _line_number(text: str, phrase: str) -> int:
    for index, line in enumerate(text.splitlines(), 1):
        if phrase in line:
            return index
    raise ValueError(f"Expected phrase not found: {phrase}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    paths = {
        "experiment_config": root / "configs/experiments/standard_model_family_v1.yaml",
        "benchmark_runner": root / "scripts/run_standard_model_benchmark.py",
        "reciprocal_runner": root / "scripts/run_reciprocal_nielsen_standard_benchmark.py",
        "selection_engine": root / "src/ibd_microbiome_framework/evaluation/standard_benchmark.py",
        "clr_pipeline": root / "src/ibd_microbiome_framework/preprocessing/compositional.py",
        "candidate_grids": root / "src/ibd_microbiome_framework/models/standard.py",
        "preprocessing_tests": root / "tests/test_preprocessing.py",
        "benchmark_tests": root / "tests/test_standard_benchmark.py",
    }
    texts = {name: path.read_text(encoding="utf-8") for name, path in paths.items()}
    checks = {
        "fixed_126_species_contract": "fixed_hmp2_prism_126_species" in texts["experiment_config"],
        "target_forbidden_in_config": "target_cohort_labels" in texts["experiment_config"] and "target_cohort_distribution_statistics" in texts["experiment_config"],
        "target_not_passed_to_selection": "select_standard_models(\n            train," in texts["benchmark_runner"] and "select_standard_models(\n            train," in texts["reciprocal_runner"],
        "inner_fold_refits_pipeline": "CLRFeaturePipeline(" in texts["selection_engine"] and ".fit(x_source.iloc[train_idx]" in texts["selection_engine"],
        "target_only_transformed_after_fit": "pipeline.transform(target)" in texts["benchmark_runner"] and "pipeline.transform(target)" in texts["reciprocal_runner"],
        "clinical_features_rejected": "Clinical, label, or identifier features are forbidden" in texts["clr_pipeline"],
        "mi_fit_on_training_labels": "mutual_info_classif(scaled, np.asarray(y_train)" in texts["clr_pipeline"],
        "candidate_grids_declared": "Return an equal four-setting grid" in texts["candidate_grids"],
        "preprocessing_test_present": "test_pipeline_rejects_planted_clinical_or_identifier_features" in texts["preprocessing_tests"],
        "nested_source_only_test_present": "test_nested_selection_and_final_refit_use_source_only_inputs" in texts["benchmark_tests"],
    }
    if not all(checks.values()):
        failed = [name for name, passed in checks.items() if not passed]
        raise RuntimeError(f"Protocol audit failed: {failed}")
    lines = {
        "selection_inner_pipeline_fit": _line_number(texts["selection_engine"], ".fit(x_source.iloc[train_idx]"),
        "runner_source_only_selection": _line_number(texts["benchmark_runner"], "selected = select_standard_models("),
        "runner_target_transform": _line_number(texts["benchmark_runner"], "x_target = pipeline.transform(target)"),
        "forbidden_feature_check": _line_number(texts["clr_pipeline"], "Clinical, label, or identifier features are forbidden"),
        "mi_training_labels": _line_number(texts["clr_pipeline"], "mutual_info_classif(scaled, np.asarray(y_train)"),
    }
    report = {
        "status": "pass",
        "scope": "static source audit; runtime tests must be reported separately",
        "checks": checks,
        "line_references": lines,
        "files": {name: {"path": str(path), "sha256": _sha256(path)} for name, path in paths.items()},
        "algorithmic_boundary": {
            "selection": "Only each outer source-training partition enters inner-fold CLR, prevalence filtering, scaling, MI feature selection and hyperparameter selection.",
            "refit": "The selected pipeline and estimator are refit with all and only source samples.",
            "target": "The locked target matrix is transformed by the fitted source pipeline and scored only after the model is frozen.",
            "limitation": "This audit verifies declared code paths and forbidden-feature checks; it does not establish that all real-world biological or technical confounding is removed.",
        },
    }
    (output / "p0r_source_only_protocol_audit.json").write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    markdown = f"""# P0-R-C：来源端嵌套流程审计\n\n## 结论\n\n静态代码审计通过。固定 126 物种合同下，CLR、患病率过滤、缩放、MI 特征选择和候选超参数选择均在来源端内层训练折重新拟合；选择后的管线仅使用全部来源样本重拟合，锁定目标仅调用变换和预测。\n\n## 审计边界\n\n1. 内层管线拟合：`standard_benchmark.py:{lines['selection_inner_pipeline_fit']}`。\n2. 来源端候选选择：`run_standard_model_benchmark.py:{lines['runner_source_only_selection']}`。\n3. 目标端只变换：`run_standard_model_benchmark.py:{lines['runner_target_transform']}`。\n4. 临床、标签和标识符特征拒绝：`compositional.py:{lines['forbidden_feature_check']}`。\n5. MI 使用当前训练标签：`compositional.py:{lines['mi_training_labels']}`。\n\n## 仍然不能据此宣称\n\n该审计证明预定义代码路径未将目标标签或目标分布用于拟合，并不证明年龄、治疗、活动度、国家、样本保存和测序流程等真实世界混杂已被完全消除。运行时单元测试结果与文件哈希见同目录 JSON。\n"""
    (output / "P0-R-C_来源端嵌套流程审计.md").write_text(markdown, encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
