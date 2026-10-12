# C2/N2：并行 Olsen 的 raw / IP / IP+Olsen 三组实验

## 结论与论文表述

在 qiumy 的 `cd303323969b01b5657d7fa1050e8a82b89d7541`（并行 Olsen 报告扫描、内部 row 索引改为 flat map）上，两套生产实验均正常完成。相对 [旧版串行结果](RESULTS.md)，报告点直接 correction 计时缩短约 43 倍，完整求解的 Olsen 增量开销从 C2 24.63% / N2 61.16% 降至 **5.14% / 23.21%**。这验证了优化有效，但 N2 的完整长轨迹仍不能无条件称为 “almost free”。旧/新实现还同时改变了数据结构，且测量不在同一时段，不能把这 43 倍解释为单独 OpenMP 改动的受控加速比。

在 1e-6 Ha 目标精度处，IP+Olsen 相对 raw 的 time-to-target 加速分别为 C2 **5.84×**、N2 **14.11×**；相对 IP-only 分别为 **1.90×**、**2.72×**。论文应区分“达到相同精度所需时间”与“跑满相同步数的成本”，同时报告两者。

建议表述：并行实现大幅降低报告点扫描成本；在 IP/PT2 基础上，Olsen 在 1e-5 至 1e-6 Ha 范围进一步缩短达到目标精度的时间，但完整轨迹的额外维护和扫描成本仍随体系及波函数规模变化。

## 任务完成、数据入口与来源

| 体系 | Slurm ID | 状态 / 退出码 | 作业总时间 | 原始数据 |
|---|---|---|---|---|
| C2/cc-pVDZ | 1007783 | COMPLETED / 0:0 | 03:28:40 | [C2 输出目录](../../experiment_results/energy_correction/1007783/c2_ccpvdz/) |
| N2/cc-pVDZ | 1007784 | COMPLETED / 0:0 | 13:35:55 | [N2 输出目录](../../experiment_results/energy_correction/1007784/n2_ccpvdz/) |

两个作业均运行于 `partition` 分区的 `bigMem5`，所有 benchmark / plotting step 都成功退出。作业总时间包含读取积分、分配和初始化、预热、三组求解和作图；下文性能表取 JSON 中的完整求解计时，不能与 Slurm 总时间混用。Slurm 查询快照保存在各任务根目录的 `slurm_accounting.psv`（时间为 Asia/Shanghai）。根目录日志为 [C2 日志](../../cdfci_three_way_1007783.out) 和 [N2 日志](../../cdfci_three_way_1007784.out)。

运行来源见 [C2 manifest](../../experiment_results/energy_correction/1007783/manifest.json)、[N2 manifest](../../experiment_results/energy_correction/1007784/manifest.json)。两次使用同一个冻结 OpenMP 可执行文件，其 SHA256 为 `5c4569b9b1ceb137e312e519e7f9fbff978d22fbe651114ab8888819535655a8`。manifest 的 `git_commit` 指向求解器运行时的 `cd30332`；当时三组 sbatch 尚未提交，因此同时保存脚本及核心源文件的 SHA256。本次发布保留原始 manifest 和脚本内容，不将运行版本改成发布结果的 commit。

## 参数和公平性

- 两个体系均为 28 个空间轨道 / 56 个自旋轨道；C2 12 个电子、N2 14 个电子，`ms2=0`。
- 64 OpenMP threads、64 coordinates、每 10000 步报告一次；每组求解使用新波函数。
- 申请 64 CPU、96 GiB、18 小时，不独占节点；`OMP_DYNAMIC=false`、`OMP_PROC_BIND=spread`、`OMP_PLACES=cores`，`srun --cpu-bind=cores`。BLAS 等额外线程限制为 1。
- C2：800000 步，`z_threshold=3e-8`，`max_wavefunction_size=1696512081`，80 个报告点。
- N2：2200000 步，`z_threshold=5e-7`，`max_wavefunction_size=848256040`，220 个报告点。
- 顺序为 raw → IP+Olsen → IP-only；前两组使用 `timed` 模式，第三组使用 `ip-only` 模式，均采用 benchmark 的标准 200 步预热。三组都开启相同的 live reporting，计时含求解期间的报告输出。

| 变体 | IP/PT2 | Olsen | 完整计时包括 |
|---|---|---|---|
| raw | 关闭 | 关闭 | 原始 CDFCI 求解 |
| IP-only | 开启 | 关闭 | IP 的迭代更新与报告，不维护 Olsen 专用内部缓存 |
| IP+Olsen | 开启 | 开启 | IP 更新、Olsen 内部缓存维护、并行报告扫描 |

同一体系的三个变体在同一 allocation / 同一节点运行，但每组只有一次完整测量，且节点非独占、顺序固定，仍受节点负载、缓存及 OpenMP 非确定性影响。N2 的 IP-only 比 raw 快约 2.1%，不应解读为开启 IP 必然降低每步计算成本；需要重复测量才能给出置信区间。`06_ip_olsen_summary.json` 的通用 timing note 提到跨作业因素，它适用于旧版独立 baseline2 测量；**本轮同一体系三组是同作业**，以 manifest 和 `runtime_summary.json` 的说明为准。

## 完整长轨迹的成本

| 体系 | raw / min | IP-only / min | IP+Olsen / min | Olsen 额外时间 / min | `(T_IP+Olsen-T_IP)/T_IP` |
|---|---:|---:|---:|---:|---:|
| C2 | 64.2917 | 69.1272 | 72.6817 | 3.5545 | 5.1420% |
| N2 | 255.8049 | 250.3949 | 308.5001 | 58.1052 | 23.2054% |

准确值见 [C2 runtime_summary.json](../../experiment_results/energy_correction/1007783/c2_ccpvdz/runtime_summary.json)、[N2 runtime_summary.json](../../experiment_results/energy_correction/1007784/n2_ccpvdz/runtime_summary.json)，同时提供 CSV。

| 体系 | 旧版报告 correction 累计 / s | 新版报告 correction 累计 / s | 旧/新比值 |
|---|---:|---:|---:|
| C2 | 871.6162 | 20.2206 | 43.11× |
| N2 | 6591.0974 | 152.0733 | 43.34× |

这里的直接报告计时不含迭代内 Olsen 缓存维护，也不等于 IP/PT2 的总成本。IP-only 的报告计时仅 C2 0.000184 s / N2 0.000545 s，因为 IP 主要在迭代更新中计算。不要用接近零的报告计时宣称 IP 免费，也不要把报告扫描的 43 倍加速写成整个求解器加速。

A/B/C 优势图的 Panel C 使用 `(T_corrected-T_raw)/T_corrected`，即 corrected 相对 raw 的累计端到端开销比例，并非上表 Olsen 相对 IP 的开销。JSON 另外保留逐报告点 `correction_seconds` 与直接修正比例。直接 correction time、端到端差值，以及所选分母须在论文图注中明确区分。

## Time-to-target：达到相同精度的时间

所有误差使用保存的长期参考能量。首次达到目标由每 10000 步的离散报告点定义，没有插值；时间含启用修正的成本。raw 列已使用 JSON 的 `baseline_variational_energy` 独立核对，其首次达标步数与现有汇总一致。

| 体系 | 目标误差 / Ha | raw / s | IP-only / s | IP+Olsen / s | raw → IP+Olsen | IP → IP+Olsen |
|---|---:|---:|---:|---:|---:|---:|
| C2 | 1e-4 | 204.018 | 54.969 | 55.521 | 3.675× | 0.990× |
| C2 | 1e-5 | 937.912 | 265.340 | 165.032 | 5.683× | 1.608× |
| C2 | 1e-6 | 3496.073 | 1138.521 | 598.340 | 5.843× | 1.903× |
| N2 | 1e-4 | 724.519 | 76.611 | 80.816 | 8.965× | 0.948× |
| N2 | 1e-5 | 4024.237 | 442.367 | 159.921 | 25.164× | 2.766× |
| N2 | 1e-6 | 14660.607 | 2821.588 | 1038.876 | 14.112× | 2.716× |

1e-4 Ha 时 IP-only / IP+Olsen 都在第一个报告点达标，不能由这一采样得出 Olsen 有额外时间优势。完整时间表的 5.14% / 23.21% 是跑满全部步数的开销，不是达到 1e-6 Ha 时的增量开销。

## 参考能量、图表和精度限制

参考未重跑，也未使用本次轨迹终点：C2 `-75.7319603747` Ha（仓库长期 223400000 步结果），N2 `-109.2821730115` Ha（400000000 步结果）。来源与说明见每个体系的 `reference.json` 和 [历史参考讨论](RESULTS.md#参考能量与图表边界)。参考末尾漂移是收敛诊断，不是精确能量误差的严格上界。

所有优势图只展示误差 **>=1e-7 Ha** 的 reference-resolved 区域；该区域本轮 corrected 误差没有反弹。原始 JSON 保留全部数据，包括低于参考的点：C2 19 个（首次 620000 步），N2 121 个（首次 1000000 步）；最大过冲分别约 1.86e-8 / 3.55e-8 Ha。不能删去这些原始记录，或将低于参考的数值视为已验证的更高精度结论。

- A/B/C：[C2](../../experiment_results/energy_correction/1007783/c2_ccpvdz/05_correction_advantage_reference_resolved.pdf)、[N2](../../experiment_results/energy_correction/1007784/n2_ccpvdz/05_correction_advantage_reference_resolved.pdf)，提供 PDF 和 PNG。
- IP/Olsen 时间消融：[C2](../../experiment_results/energy_correction/1007783/c2_ccpvdz/06_ip_olsen_timing.pdf)、[N2](../../experiment_results/energy_correction/1007784/n2_ccpvdz/06_ip_olsen_timing.pdf)，提供 PDF 和 PNG。
- Panel A 的横轴为非零变分系数数目 `stored_determinants`，Panel C 使用 `stored_wavefunction_entries`，两者不能混称同一个波函数大小。
- 当前 A/B 图沿用现有分析脚本：raw 能量取 joint 运行的未修正能量，raw 时间取独立 raw；JSON 同时保存 `baseline_variational_energy` 以便正式论文使用独立 baseline 能量重画。新 raw/joint 未修正能量最大差为 C2 2.80e-9 / N2 2.59e-8 Ha，目标精度步数已核对一致。

## 重分析与重跑

拉取后无需 FCIDUMP、集群或二进制即可查看全部图表，并由原始 JSON 重新分析（需要 Python、NumPy、Matplotlib）：

```bash
python3 experiments/energy_correction/summarize_large_scale.py \
  experiment_results/energy_correction/1007783/c2_ccpvdz/05_large_scale_with_timings.json \
  experiment_results/energy_correction/reanalysis_parallel_c2
python3 experiments/energy_correction/compare_ip_olsen.py \
  experiment_results/energy_correction/1007783/c2_ccpvdz/05_ip_only.json \
  experiment_results/energy_correction/1007783/c2_ccpvdz/05_large_scale_with_timings.json \
  experiment_results/energy_correction/reanalysis_parallel_c2
```

N2 对应替换 `1007783/c2_ccpvdz` 为 `1007784/n2_ccpvdz` 和输出名。原始结果保持不变；新重分析目录继续被 gitignore 忽略。

重新跑大任务使用 [README](README.md) 中的 `run_three_way.sbatch`，需预先构建 OpenMP benchmark 并本地提供相同 FCIDUMP；可用 `CDFCI_BENCHMARK_APP` 指定冻结可执行文件。节点和分区是本集群配置，其他集群应自行调整 Slurm 资源设置。提交只包含六个明确允许的正式结果目录、对应日志与实验脚本/文档，不上传大积分文件、构建目录或可执行文件。

本轮验证范围：Slurm 全部 step 的成功退出、80/220 个完整报告、raw 无 correction evaluation、IP 的 Olsen 开关关闭且所有 internal correction 为零、联合修正开启、三组参数一致、manifest 核心源文件和脚本哈希匹配、原始 JSON 重分析通过。没有修改求解器或修正算法，也没有新增大任务或重新生成参考能量。
