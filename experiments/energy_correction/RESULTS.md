# C2/N2 能量修正实验：可用于论文的结果与限制

这批结果支持两点：IP/PT2 与 Olsen 联合修正能显著减少达到目标能量精度的时间；在 IP/PT2 基础上加入 Olsen 仍有额外收益。但完整长轨迹不支持 Olsen 无条件地 “almost free”：相对 IP-only，C2 总时间增加 24.6%，N2 增加 61.2%。论文应同时报告精度收益和规模相关的计算成本。

## 数据入口与完成情况

四个正式任务均在 `bigMem5` 上成功完成，Slurm 状态为 `COMPLETED`，退出码为 `0:0`。表中的任务时间包含初始化、预热、求解及结果汇总；后文性能表使用基准 JSON 中的完整求解时间。

| 系统 | 任务 | 运行内容 | Slurm 总时间 | 输出目录 |
|---|---|---|---|---|
| C2/cc-pVDZ | 988519 | raw + IP/PT2+Olsen | 02:58:45 | [988519](../../experiment_results/energy_correction/988519/) |
| N2/cc-pVDZ | 988520 | raw + IP/PT2+Olsen | 11:12:51 | [988520](../../experiment_results/energy_correction/988520/) |
| C2/cc-pVDZ | 990772 | baseline2：仅 IP/PT2 | 01:19:17 | [990772](../../experiment_results/energy_correction/990772/) |
| N2/cc-pVDZ | 990773 | baseline2：仅 IP/PT2 | 04:20:31 | [990773](../../experiment_results/energy_correction/990773/) |

根目录的 `cdfci_pe_<job-id>.out` 保留任务参数、完成日志与打印出的汇总。四个输出目录保留原始 JSON、CSV、PDF/PNG、参考说明和运行 manifest，不需要重新执行大任务即可重分析。

## 参数与比较方式

- 两个系统均为 28 个空间轨道、56 个自旋轨道；C2 为 12 个电子，N2 为 14 个电子，`ms2=0`。
- 64 OpenMP threads、64 coordinates、每 10000 步报告一次；`OMP_DYNAMIC=false`、`OMP_PROC_BIND=spread`、`OMP_PLACES=cores`，`srun --cpu-bind=cores`。
- 每个 sbatch 请求 96 GiB、18 小时，非独占节点。
- C2：800000 步，`z_threshold=3e-8`，`max_wavefunction_size=1696512081`，共 80 个报告点。
- N2：2200000 步，`z_threshold=5e-7`，`max_wavefunction_size=848256040`，共 220 个报告点。

| 变体 | `enabled` | `olsen_enabled` | 实际计算 |
|---|---|---|---|
| raw / baseline | false | 不生效 | 不维护 IP/PT2，不计算 Olsen |
| IP-only / baseline2 | true | false | 维护 IP/PT2，不维护 Olsen 专用内部缓存，不执行 Olsen 扫描 |
| IP+Olsen / corrected | true | true | IP/PT2 + 内部 Olsen 缓存维护与报告点扫描 |

所有变体从新波函数开始，使用相同求解参数。raw 与 IP+Olsen 在同一任务中依次求解；baseline2 在后续独立任务中运行，与保存的 IP+Olsen 数据比较。每个变体只有一次完整测量；跨任务时间差包含节点负载、缓存与 OpenMP 非确定性因素，不能当作带置信区间的性能估计。

## 完整长轨迹的时间开销

| 系统 | raw 时间 / min | IP-only 时间 / min | IP+Olsen 时间 / min | Olsen 相对 IP 的额外时间 |
|---|---:|---:|---:|---:|
| C2 | 79.2873 | 78.2114 | 97.4727 | 24.6272% |
| N2 | 252.4384 | 259.9813 | 418.9891 | 61.1613% |

Olsen 增量开销定义为 `(T_IP+Olsen - T_IP) / T_IP`，包括其内部缓存维护和扫描。原 A/B/C 图的 Panel C 使用另一种分母：`(T_corrected - T_raw) / T_corrected`；这两个指标不能混用。

IP+Olsen 直接修正计时累计为 C2 871.616 s（14.53 min）、N2 6591.097 s（109.85 min）；IP-only 的报告点计时累计分别为 0.000197 s 和 0.000517 s。IP/PT2 的计算发生在迭代更新中，因此其报告点计时不能代表 IP/PT2 的总成本。同理，直接修正计时不包含 Olsen 的迭代内缓存维护，完整求解时间差才是总成本比较。

当前 Olsen 在报告点串行遍历内部非零系数，成本随其数量增长。末次报告的修正耗时约为 C2 21.54 s、N2 58.32 s。因此这批数据不能支持整个长轨迹 “almost free”。

## 达到相同目标精度的收益

表中的首次达到精度均以离散报告点为准，时间已包含启用修正的成本；每 10000 步采样一次，未对真实首次跨越时间插值。

| 系统 | 目标误差 / Ha | raw / s | IP-only / s | IP+Olsen / s | raw → IP+Olsen 加速 | IP → IP+Olsen 加速 |
|---|---:|---:|---:|---:|---:|---:|
| C2 | 1e-4 | 248.914 | 63.205 | 66.516 | 3.742× | 0.950× |
| C2 | 1e-5 | 1156.274 | 306.799 | 195.470 | 5.915× | 1.570× |
| C2 | 1e-6 | 4319.554 | 1311.446 | 713.630 | 6.053× | 1.838× |
| N2 | 1e-4 | 650.021 | 77.545 | 80.017 | 8.124× | 0.969× |
| N2 | 1e-5 | 3749.435 | 446.490 | 159.613 | 23.491× | 2.797× |
| N2 | 1e-6 | 14460.382 | 2904.853 | 1062.303 | 13.612× | 2.734× |

在 IP+Olsen 首次达到 1e-6 Ha 的步数处，加入 Olsen 的同步数累计时间增幅为 C2 7.42%、N2 11.20%。这说明目标精度附近的增量成本较低，但仍应给出数字而非宣称成本为零。在 1e-4 Ha 时，IP-only 和 IP+Olsen 都在首个报告点达标，采样结果没有显示额外的时间优势。

建议论文表述：**相对 IP/PT2 单独修正，加入 Olsen 在 1e-5 至 1e-6 Ha 的目标精度下进一步缩短达到目标的时间；这一收益伴随随内部波函数规模增加的维护和扫描成本。**

## 参考能量与图表边界

使用仓库长期 CDFCI 回归计算的能量，不重新运行参考：

- C2：`-75.7319603747` Ha，223400000 步，来源 `regression_tests/c2/ccpvdz_psi4/output`。
- N2：`-109.2821730115` Ha，400000000 步，来源 `regression_tests/n2/ccpvdz_psi4/output`。

两个参考最后 100 个报告点的能量漂移分别约 4.9e-9 和 5.9e-9 Ha；这是收敛诊断，不是参考真实误差的严格上界。误差图仅保留绝对误差 >=1e-7 Ha 的区域。在该区域 IP+Olsen 误差无反弹、无参考能量以下过冲；更低区域的微小过冲不用于精度论证。不要使用本次有限迭代终点作为误差参考。

不同 OpenMP 运行中的变分能量有微小差异。raw/IP+Olsen 最大差异为 C2 3.33e-9、N2 2.08e-8 Ha；IP-only/IP+Olsen 最大差异为 C2 6.41e-9、N2 3.82e-8 Ha。当前 A/B 图 raw 能量取自 corrected 运行的未修正能量，raw 时间取自独立 baseline；以独立 baseline 能量重新核对后，上表目标精度的首次达到步数相同。正式论文作图可统一使用独立 baseline 能量，避免混合来源。

## 图表入口与重分析

- 原始优势图：[C2 A/B/C](../../experiment_results/energy_correction/988519/c2_ccpvdz/05_correction_advantage_reference_resolved.pdf)、[N2 A/B/C](../../experiment_results/energy_correction/988520/n2_ccpvdz/05_correction_advantage_reference_resolved.pdf)。Panel A 横轴是非零变分系数数目，Panel C 横轴是存储残差/波函数条目数。
- Olsen 时间消融：[C2 IP vs IP+Olsen](../../experiment_results/energy_correction/990772/c2_ccpvdz/06_ip_olsen_timing.pdf)、[N2 IP vs IP+Olsen](../../experiment_results/energy_correction/990773/n2_ccpvdz/06_ip_olsen_timing.pdf)。
- 全精度指标见各目录的 `05_large_scale_summary.json` 和 `06_ip_olsen_summary.json`，逐点数据见 CSV 或原始 JSON。

从仓库根目录重分析，不需要 Slurm 或 FCIDUMP：

```bash
python3 experiments/energy_correction/compare_ip_olsen.py \
  experiment_results/energy_correction/990772/c2_ccpvdz/05_ip_only.json \
  experiment_results/energy_correction/988519/c2_ccpvdz/05_large_scale_with_timings.json \
  experiment_results/energy_correction/reanalysis_c2
python3 experiments/energy_correction/compare_ip_olsen.py \
  experiment_results/energy_correction/990773/n2_ccpvdz/05_ip_only.json \
  experiment_results/energy_correction/988520/n2_ccpvdz/05_large_scale_with_timings.json \
  experiment_results/energy_correction/reanalysis_n2
```

重跑大任务的方法见 [README](README.md) 的 `large` / `baseline2` 提交说明。manifest 中的 `git_commit` 记录运行时 HEAD；baseline2 运行使用了尚未提交的基准辅助代码，其 `source_sha256` 记录了对应源文件的准确内容。求解器与修正算法未因这次消融实验修改。

这批结果覆盖论文第五项大规模实验的性能与精度问题；不替代合成模型阶数、小 CI 初值或失败区间实验。投稿前还需要根据 CSV/JSON 补齐 iteration/Hamiltonian-work 展示、分子几何与参考收敛说明，并明确上述采样和计时限制。
