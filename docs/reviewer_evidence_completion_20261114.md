# 审稿人缺失证据补全审计（2026-11-14）

## 结论

当前信息足以完成两项补全：

1. 已知真值正信号的逐种子验收规则和 80%/70%/100% 结果可以从端到端 truth-panel ledger 独立复算；
2. 报告模型的损失函数、类别权重、batch、随机种子、拓扑干预和 bootstrap 单位可以从代码、manifest 与 compact ledger 对齐。

但拓扑证据尚未形成完整的逐边审计闭环。正式 Gate B 的指标文件足以重建 Table 3、Figure 6 以及相应的 Supplementary Figure S3/S4 数值摘要，却没有导出每个 source split × rewiring seed 的正式 edge list，也没有在指标文件中写入 `self_loop_count` 和 `unique_edge_count`。此外，当前重建哈希记录有 24 条不匹配：HMP2→PRISM split 9 的 real+5 rewired，以及 PRISM→HMP2 splits 6–8 的 real+5 rewired。因而在这些记录被重新导出或逐条解释前，不能把整套拓扑干预描述为“已完全可审计”。

这不是要求更换模型或增加第三队列，而是要求把已运行的拓扑结果与其实际边序列重新对齐。

## 1. 拓扑约束与可重建性

### 已核验的证据

- 正式 Gate B `edge_ablation_metrics.csv` 有 140 行：20 个 real、20 个 identity、100 个 rewired 记录。
- 每个方向和 source split 均有 5 个 rewired seeds；指标记录了 `edge_count`、出度范围/方差、`degree_hash`、`edge_hash`、`rewire_swaps_accepted`、模型参数量和 feature/contract 标识。
- 每个 formal rewired 记录的边数等于 `8 × n_nodes`，出度最小值和最大值均为 8，出度方差为 0；同一 split 内的 degree hash 与 real 图一致，edge hash 与 real 图不同。
- 显式 identity 图只有自环，因此在相同节点标识和有向边定义下，real/rewired 与 identity 的 Jaccard 必须为 0。历史 edge audit 的 60 条记录均报告 `self_loop_count=0` 和 `jaccard_with_identity=0`。

### 当前不能宣称的内容

- 正式 Gate B 没有保存每个 rewired graph 的逐边 CSV；仅凭 edge hash 无法让第三方独立核验无自环、无重复边和精确的入度约束。
- `identity_jaccard_hash_audit.csv` 的 140 条记录中有 24 条 `hash_matches=False`，分布为 HMP2→PRISM split 9（6 条）和 PRISM→HMP2 splits 6、7、8（各 6 条）。identity 行匹配，说明问题集中在 real/rewired 的边序列重建，而不是 identity 自环定义。
- `closure_final/topology_edges_20261110_131203` 的旧 shuffled 文件不能替代 formal Gate B 的 rewired 导出；抽查显示存在重复有向边（例如 HMP2 split 1 的文件行数大于唯一边数）。因此它们不满足“simple directed graph”审计要求。

### 对正文图表的影响

`edge_ablation_metrics.csv`、`gate_b_common_draw_summary.csv`、`topology_split_stability.csv`、identity-aware `metrics.csv/summary.csv` 以及 PRISM provenance-stratum 的 `metrics.csv/summary.csv` 已足以重建统计数值和绘图源表。它们不能替代正式 rewired edge list 的逐边核验。投稿前必须二选一：

- 从同一次正式运行导出每个 split × seed 的 edge list 和约束审计表，并解决 24 条哈希不匹配；或
- 将这些记录从“已审计拓扑证据”中隔离，在正文和图注中只保留条件性、实现特定的数值结果，并明确拓扑边界尚未完成逐边复核。

## 2. 已知真值的正信号验收规则

端到端 truth panel 的每个场景有 10 个 seed，疾病真值分母为 6 个已知疾病特征。正信号规则是同一 seed 上同时满足：

\[
A_s = I(\mathrm{AUROC}_s \ge 0.70) \land I(R_s \ge 3/6).
\]

这是 AND，不是 OR。负对照不使用疾病特征恢复条件，而采用包含端点的 AUROC 区间：

\[
0.35 \le \mathrm{AUROC}_s \le 0.65.
\]

接受率为通过 seed 数除以该场景的 10 个 seed。主端到端 truth panel 可复算为：

| 场景 | 通过 seed / 总 seed | 接受率 |
|---|---:|---:|
| signal-only | 8/10 | 80% |
| signal-reversing-batch | 7/10 | 70% |
| signal-topology-shift | 10/10 | 100% |
| batch-only（负对照） | 10/10 | 100% |
| null（负对照） | 10/10 | 100% |

这里的 80%、70% 和 100% 是运行特征验收率，不是 AUROC 置信区间，也不是现实队列中的疾病特征发现率。`selector_audit` 是独立的选择器审计，不得与正文使用的 `end_to_end_truth_panel` 混写。

逐 seed 输入、规则和结果位于：

`results/reconciled_evidence/truth_acceptance/per_seed_acceptance.csv`

规则定义位于：

`results/reconciled_evidence/truth_acceptance/truth_acceptance_rules.json`

## 3. 训练配置与 bootstrap 口径

### 拓扑干预 / identity-aware GAT

- 输入：126 个共享微生物特征；closure/CLR/流行率/MI 参数仅由 source training panel 拟合。
- real 图：source-training Spearman top-k，`k=8`，不自动添加 self-loop。
- identity 图：显式对角 self-loop。
- rewired 图：有向 degree-preserving double-edge swap；运行时拒绝 self-loop 和重复边，并检查图约束。
- 损失：`binary_cross_entropy_with_logits`，不使用 `class_weight` 或 `pos_weight`。
- 优化：Adam，learning rate 0.001，weight decay 0.001，40 epochs。
- Gate B 训练 batch：整个 source panel 一次更新；预测 batch：128。
- identity-aware sensitivity：同一 loss/optimizer/epoch budget，8 维可学习 node-identity 表示；其 real/rewired/identity 结果应标为敏感性分析，不作为原始 GAT 绝对性能复现。
- seed：PyTorch 与 NumPy 均按 source split 设置；模型 seed 为 `base + direction_index × 100000 + split`，rewire seed 为该 split 的显式派生 seed。

### 直接 PRISM provenance-stratum transfer

该 stress test 使用 GAT mini-batch 32、预测 batch 64，并对每个 source split 单独拟合；表格参照模型的 `class_weight`、随机状态和迭代上限保留在对应 benchmark 脚本中。不能把该 stress test 的训练 batch 与 Gate B full-source batch 混称为同一配置。

### bootstrap 单位

Gate B 的共同目标抽样使用 5,000 次 participant-cluster bootstrap；同一方向、目标层级和 draw 在十个 source split 及各拓扑条件间共享，`n_invalid_draws=0`。PRISM provenance-stratum transfer 的正式摘要使用 1,000 次目标参与者抽样。visit-level 和 participant-level 结果即使目标队列相同，也不应强行相等：HMP2 目标面板为 1,626 个 visit、130 名 participant，participant-level 先按参与者聚类抽样并等权汇总；visit-level 直接在观测行上计算。因此小幅差异来自评价单位和聚合顺序，而不是同一抽样条件下的随机漂移。

配置记录：

`repro/audits/training_configuration.json`

## 4. 表图来源注册

为避免手工改数字，复现包新增 `evidence_registry.csv`，明确每个对象的紧凑源表：

| 对象 | 主要源表 | 可重建内容 | 当前限制 |
|---|---|---|---|
| Table 3 | `absolute_performance_summary.csv`、`paired_bootstrap_reconciled_by_split.csv` | 方向、模型、AUROC/AP/Brier 和配对区间 | 不包含原始预测 |
| Figure 6 | `topology_split_stability.csv`、`topology_constraints/` | split-level topology metrics、rewired summary、identity overlap | 逐边 formal rewired audit 未闭环 |
| Supplementary Figure S3 | identity-aware `metrics.csv`、`summary.csv`、`paired_split_contrasts.csv` | real/rewired/identity 的十次 split 比较 | 24 条 hash mismatch 需解释 |
| Supplementary Figure S4 | `prism_discovery_origin_sensitivity/metrics.csv`、`summary.csv` | 两个 PRISM provenance strata 的四模型 transfer 摘要 | 不能解释为地理或来源组成的因果效应 |

## 5. 格式与参考文献验收

最新正文的 Normal 样式为 Times New Roman 10.5 pt、1.25 倍行距；这可用于内部阅读，但正式投稿版应按 BIB 当前指南确认双倍行距要求。正文 46 条参考文献的编号和首次出现顺序已通过结构审计，但目前仍是 paragraph-level scope screen，不是逐条 source-text claim certification。提交前应对每个 Introduction/Methods 论断逐句核对 DOI/摘要或全文，并拆分过宽的群引文；Gate A/B 数字只能由本地 ledger 和 manifest 支撑，不能用一般方法论文替代。

BIB 的 Problem solving protocols 要求真实数据上的方法学应用、足够的重实现信息以及可用的代码/数据入口。近年的 BIB 实践通常采用“核心代码 + 小型测试/汇总数据 + 公开 accession 或下载脚本 + 一条运行命令”，而不是上传受限原始组学、患者标识、HPC 日志或所有 checkpoint。可参照：[BIB Author Guidelines](https://academic.oup.com/bib/pages/author-guidelines)、[Multilevel-GNN](https://pmc.ncbi.nlm.nih.gov/articles/PMC11052635/)、[InterDiff](https://github.com/zephyrdhb/InterDiff)、[PBIP](https://github.com/a1678019300/PBIP) 和 [HyperGraphComplex](https://github.com/LiDlab/HyperGraphComplex)。

## 6. 最终判断

- 真值验收：**信息充分，已可独立复算**。
- 训练配置：**信息充分，已整理到最小复现包**。
- Table 3/Figure 6/S3/S4 的统计源表：**可重建**。
- 正式 rewired 图的逐边约束与 hash：**尚未闭环，当前不能声称完整 auditable**。
- 是否需要新模型或第三队列：**本项核查不需要**。
- 投稿前必做：补出正式 edge list/constraint audit，解决或隔离 24 条 mismatch；随后再进行最终图表和正文数字冻结、双倍行距渲染及逐条引用匹配。

## 7. 复现包已更新内容与验证结果

已在 `D:\GitHub_workplace\ibd_GAT` 仅加入审稿核验所需的紧凑层：

- `results/reconciled_evidence/truth_acceptance/truth_acceptance_rules.json`：机器可读的 AND/区间规则、6 个真值特征分母和 10-seed 汇总；
- `results/reconciled_evidence/truth_acceptance/per_seed_acceptance.csv`：补充 `truth_denominator` 和 `acceptance_logic` 字段；
- `results/reconciled_evidence/topology_constraints/topology_audit_status.json`：明确 formal Gate B 逐边导出缺失和 24 条 hash mismatch；
- `results/reconciled_evidence/identity_aware_topology/`：仅保存 Figure S3 所需的 compact metrics/summary/paired split 表，不含患者标识、原始丰度、checkpoint 或完整预测；
- `results/reconciled_evidence/evidence_registry.csv`：Table 3、Figure 6、S3、S4 的源表映射；
- `repro/audits/validate_evidence_registry.py`、`validate_truth_acceptance.py`、`validate_topology_constraints.py`：分别检查源表存在性、逐种子真值规则和拓扑审计状态；
- `repro/audits/training_configuration.json`：补充 Gate B full-source batch、预测 batch、loss、`pos_weight`、seed 派生规则、线程数与 bootstrap 单位。

验证结果：evidence registry 通过；truth acceptance 通过并复算 80%、70%、100%；topology validator 在默认严格模式下按预期失败，原因是正式 edge-list 缺失和未解决 hash mismatch，在 `--allow-hash-failures --allow-incomplete` 的“仅指标层”模式下才通过。该失败是保护机制，不能改成默认通过。

复现包没有提交 HPC 日志、Slurm 脚本输出、原始队列矩阵、participant identifier、全量 prediction ledger 或 checkpoint，符合 BIB protocol 仓库通常采用的“核心代码 + compact evidence + public accession/download instructions”模式。
