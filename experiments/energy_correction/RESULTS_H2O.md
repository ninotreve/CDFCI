# H2O/cc-pVDZ：500 万步同波函数三种能量估计

## 完成情况与数据入口

任务 **1027249** 于 2026-10-19 14:38:59 至 23:00:18（Asia/Shanghai）运行，
Slurm 状态 `COMPLETED`、退出码 `0:0`。完成 5,000,000 个 block，
每 10,000 步报告，共 500 个点，没有缩短步数或提前收敛停止。

- [正式结果目录](../../experiment_results/energy_correction/h2o_ccpvdz/1027249/)：
  原始 `runs/joint.json` / `runs/joint.json.jsonl`、逐点 CSV、目标精度表、摘要、图表。
- [任务日志](../../cdfci_h2o_5m_1027249.out)及目录内 `slurm_accounting.psv`。
- [误差 vs 迭代/联合运行时间](../../experiment_results/energy_correction/h2o_ccpvdz/1027249/signed_errors.pdf)、
  [误差 vs Hamiltonian 列数/非零系数](../../experiment_results/energy_correction/h2o_ccpvdz/1027249/errors_vs_work.pdf)，均提供 PNG。
- `manifest.json`、`plan.json`、`command.json` 保留真实运行版本、未提交状态、
  源码/二进制/输入 SHA256、运行参数与完整命令，不改写为发布结果的 commit。

## 实验设计与资源

**本次只有一次 IP+Olsen 求解，不是三组独立 benchmark。**
在每个报告点，对同一个波函数读取变分能量 E、E+IP、E+IP+Olsen；
外部 IP 保留既有增量分母更新实现，Olsen 使用报告时的当前 Rayleigh 能量。
没有独立 raw/IP-only 求解，没有校准、预热或重复计时。

- 24 个空间轨道 / 48 个自旋轨道，10 个电子，`ms2=0`。
- 分区 `partition`，节点 `bigMem3`，64 OpenMP threads / 64 coordinates；
  申请 128 GiB、18 小时，非独占。64 线程不代表已证明占用 64 个独立物理核。
- `OMP_PROC_BIND=spread`、`OMP_PLACES=cores`、`OMP_DYNAMIC=false`，
  `srun --cpu-bind=verbose,cores`；BLAS 等额外线程限制为 1。
- `z_threshold=0`，容量 848256040，`stopping_dx_threshold=0`。
- Slurm 总时间 **08:21:19**；benchmark 完整 solve timer **30066.52569707 s**。
  最后报告时钟 **30053.097429415 s**，不与完整 timer 混用。
- 主求解 step 峰值 RSS **54818956 KiB ≈ 52.28 GiB**。
- 最后非零变分系数数目 **260459958**；存储波函数/残差条目数 **451681060**，
  两者不是同一个“determinant size”；累计 Hamiltonian columns **319999937**。

## 最后报告点的能量

参考值 **−76.2418600697 Ha**，来自用户确认具有严格精度认证的历史参考，
不是本轮轨迹终点。原始参考日志随数据发布于 `reference/h2o_ccpvdz_new.out`，
SHA256 为 `a67420d891f4aeb6a9a7953c471d0d6242831fe3ea0804724ef305e1e150667d`。
认证来源按 manifest 如实记录为用户确认；这里没有重新计算或独立认证参考。

| 能量估计 | 500 万步能量 / Ha | 相对给定参考的有符号误差 / Ha |
|---|---:|---:|
| E | −76.24186006965930 | +4.06999e-11 |
| E+IP | −76.24186006966218 | +3.78151e-11 |
| E+IP+Olsen | −76.24186006966511 | +3.48876e-11 |

全部 500 个点均满足 `E_ref < E+IP+Olsen < E+IP < E`。
联合修正末端有 15 次极微小上跳，最大 **1.42109e-13 Ha**，不是宏观能量反弹。
最后一步不是严格认证的收敛极限。参考仅打印到小数点后 10 位，
末位四舍五入半单位为 **5e-11 Ha**，因此不对末端几位差别作高精度物理结论。
图不任意截掉 1e-7 Ha 以下数据；保留正负号和所有点，symlog 的线性尺度
`5e-11` 仅表示打印分辨率，不是参考真实误差的认证上界。

## 达到目标误差所需工作量

下表是每 10,000 步采样的**首次记录达标点**，不是插值后的精确跨越位置。

| 目标绝对误差 / Ha | E / block | E+IP / block | E+IP+Olsen / block |
|---|---:|---:|---:|
| 1e-5 | 50000 | 20000 | 10000 |
| 1e-6 | 160000 | 60000 | 30000 |
| 1e-7 | 400000 | 210000 | 100000 |
| 1e-8 | 850000 | 540000 | 280000 |
| 1e-9 | 1730000 | 1190000 | 630000 |
| 1e-10 | 3300000 | 2640000 | 1440000 |

在 1e-9 Ha 处，联合修正所需迭代数约为 E 的 **36.4%**、E+IP 的 **52.9%**。
1e-10 Ha 门槛按给定参考值如实导出，但解释时须同时注明上述参考打印分辨率。

## 性能结论边界

500 次报告的直接 correction timer 累计 **1254.459225981 s（20.91 min）**，
约为完整求解时间的 **4.1723%**。此 timer 包括报告时 Olsen 延迟刷新及评估，
但不含全部迭代内 IP/Olsen 维护；它不是 Olsen 相对 IP-only 的端到端增量成本。

三种能量共用联合求解的时间轴。可以展示同波函数的精度改善和达到目标精度
所需迭代/Hamiltonian 工作量的减少；不能将横轴时间解读为三个独立算法的
运行时间，不能从本轮单次任务推导 standalone speedup 或证明 Olsen “almost free”。
独立计时证据仍应使用带各自限制说明的 [C2/N2 数据](RESULTS_DEFERRED.md)。

## 重分析与验证

查看数据和重画图不需要 FCIDUMP、Slurm 或原始绝对路径。
完整命令见 [README_h2o.md](README_h2o.md)；建议先复制正式目录到被忽略的
`experiment_results/energy_correction/reanalysis_h2o` 再运行 `h2o_single_run.py analyze`。

```bash
python3 -m unittest discover -s test -p 'test_h2o_*.py' -v
```

发布检查覆盖 JSON/JSONL/CSV 一致性、500 个报告点、修正开关、能量分解、
首次达标表、参考原始日志 SHA256，以及冻结核心源码与本次发布源码的哈希一致性。
这次整理没有修改求解器或 correction 公式，也没有重新运行大实验。
