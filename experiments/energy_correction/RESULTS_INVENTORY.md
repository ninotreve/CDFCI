# 实验结果保留与清理清单（2026-10-20）

## 保留

| 数据 | 定位 | 保留原因 |
|---|---|---|
| H2O 1027249，500 万步 | [RESULTS_H2O.md](RESULTS_H2O.md) | 本轮完整同波函数精度实验 |
| C2 1022419，300 万步前缀 | [RESULTS_DEFERRED.md](RESULTS_DEFERRED.md) | 选定的同节点三组性能数据；原任务 TIMEOUT 状态保留 |
| N2 1022420，260 万步 | [RESULTS_DEFERRED.md](RESULTS_DEFERRED.md) | 完整三组性能数据 |
| 988519 / 988520 / 990772 / 990773 | [RESULTS.md](RESULTS.md) | 串行 Olsen 的历史成本证据，供优化前后对照 |
| 1007783 / 1007784 | [RESULTS_PARALLEL.md](RESULTS_PARALLEL.md) | 并行 Olsen 的历史成本证据，供延迟刷新前后对照 |

历史数据不是本轮正式结果，但仍有科学用途，不按日期简单删除。
回归测试输入/参考、正式结果来源、源码冻结快照及 FCIDUMP 不属于临时结果清理对象。
构建目录与积分仍只留本地，不上传二进制或 FCIDUMP。

## 已移出仓库的无用输出

| 原路径（相对仓库根目录） | 原因 |
|---|---|
| `cdfci_h2o_corr_1026536.out` | Slurm spool 路径错误，启动即失败，无有效求解 |
| `experiment_results/energy_correction/h2o_ccpvdz/1027218/`、`cdfci_h2o_corr_1027218.out` | 已取消、短轨迹，旧实验设计已被替代 |
| `experiment_results/energy_correction/h2o_ccpvdz/1027222/`、`cdfci_h2o_corr_1027222.out` | 已取消，九次求解方案未完成，不能作为完整消融实验 |
| `experiment_results/energy_correction/h2o_ccpvdz/analysis_1027222_20261019_141715/` | 对被替代取消任务的重复临时分析 |
| `experiment_results/energy_correction/baseline2_validation/`、`baseline2_validation_ip.json`、`baseline2_validation_ip_olsen.json` | 旧的小规模脚本冒烟输出，非正式论文数据 |
| `experiment_results/energy_correction/three_way_validation/`、`three_way_validation_ip.json`、`three_way_validation_paired.json` | 旧三组流程临时验证输出，已有回归测试替代 |

以上均为原先被 gitignore 忽略的本地结果，所以 PR 不会出现相应的 Git 删除记录。
共 12 个顶层目标，约 6.2 MiB；已从 CDFCI 工作树移除，而非改名后留在正式目录。
为避免误删不可恢复，当前集群暂存恢复副本于仓库外：
`/home/yjzhang/codes/cdfci/experiment_archive/20261020-obsolete-results.nsDSoG/`。
恢复副本不提交、不推送；其中 `energy_correction/` 和 `logs/` 保留相应目录/日志。

`.gitignore` 仅白名单放行本次 `h2o_ccpvdz/1027249/` 及其指定任务/参考日志；
以后生成的 H2O 测试、失败或其他 job 目录仍默认忽略。
