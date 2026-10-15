# 最新正式结果：C2 300万步 / N2 260万步（延迟 Olsen 刷新）

## 数据选择与结论

本轮只发布下列两套数据，供合作者直接查看图表和写论文：

| 体系 | 采用步数 / 报告点 | 来源 | 三组是否同节点、同 allocation |
|---|---|---|---|
| C2/cc-pVDZ | 3000000 / 300 | 原任务 1022419 的完整前缀 | 是，bigMem5 |
| N2/cc-pVDZ | 2600000 / 260 | 完整任务 1022420 | 是，bigMem5 |

两套均比较 raw、IP-only、IP+Olsen。**C2 不使用 bigMem3 补跑的数据，也不使用最后60万步**。原 C2 请求360万步，raw / 联合修正完成，IP-only 在347万步报告后超时；前三百万步的三组数据完整。原 Slurm 状态仍是 TIMEOUT，不改写为 COMPLETED；N2 完整完成且退出码0。

本轮结果支持低端到端 Olsen 增量开销：C2 在300万步处为 **+3.2760%**，N2 完整260万步为 **−1.9691%**。N2 负差值只表示这次联合运行略快，不能解释为 Olsen 算法有负计算成本。固定顺序、节点非独占且每组仅一次测量，仍有负载、缓存及 OpenMP 波动；不提供置信区间或无条件“免费”的结论。

## 文件入口

- [C2 正式目录](../../experiment_results/energy_correction/1022419/c2_ccpvdz/prefix_3000000/)：`prefix_summary.json`、三组逐点 JSON/CSV、目标精度 CSV、两类 PDF/PNG，以及仅含前300万步的 `source_paired_prefix.json` / `source_stdout_prefix.out`。
- [N2 正式目录](../../experiment_results/energy_correction/1022420/n2_ccpvdz/)：原始 paired JSON、IP-only JSON、runtime summary、逐点 CSV、目标精度表、PDF/PNG 和参考说明。
- 各 job 根目录保留运行 manifest 和 `slurm_accounting.psv`。N2 [原始任务日志](../../cdfci_three_way_1022420.out) 一并提交。
- C2 原全长 paired JSON/完整 stdout 不上传；发布源片段保持选中报告点的 JSON 精度和 stdout 原行，保存全长源文件原始 SHA256，供本地档案审计。源片段不是伪造的完整360万步 benchmark，省略所有最终全长 solve summary。
- 已取消的补跑 `1024227` 输出、日志、专用脚本已移出仓库，不纳入 PR。之前六套已发布历史结果仍保留用于版本比较，不属于本轮正式结论。

## 版本、参数与可复现性

运行求解器为上游 `5e5cc06b4ac1b836bf06be5b1c185ef9922c1a0f`：Olsen 内部 b 在迭代中按线程排队、借元数据位去重，到报告点再按4096行块并行刷新；随后用原 Olsen 公式和高精度归约。IP 增量更新机制不变。加入 checkpoint 刷新、元数据与并发去重测试，以及刷新统计。

两体系用同一冻结 OpenMP 二进制，SHA256 `48f4dd637b803ebd5ddfaa2c31a296ba054ada3640441af38b843f3edcbd0d70`；GCC 11.4、Release `-O3 -DNDEBUG`。运行 manifest 的 git commit、源码与提交脚本哈希代表真实运行版本，不能用本次发布 commit 或调整后脚本的哈希替换。

- 28空间轨道 / 56自旋轨道；C2 12电子、N2 14电子，ms2=0。
- partition / bigMem5；64 CPU / OpenMP threads、64 coordinates、96 GiB、18小时、非独占。
- `OMP_DYNAMIC=false`、`OMP_PROC_BIND=spread`、`OMP_PLACES=cores`、`srun --cpu-bind=cores`，BLAS等额外线程设为1。
- 每10000步报告；C2 z_threshold=3e-8、capacity=1696512081；N2 z_threshold=5e-7、capacity=848256040。
- 顺序 raw → IP+Olsen → IP-only，每组新波函数，相同 live reporting。上游 timed 模式有200步 raw 预热，IP-only 不预热；保留并明示这一行为。
- raw 全部修正关闭；IP-only 的 IP开启、Olsen关闭；联合变体同时开启 IP 和 Olsen。

后续 `run_three_way.sbatch` 默认值已改为 **C2 300万 / N2 260万步**，其他参数不变；这只影响重跑，不篡改原C2的360万步申请记录。脚本支持冻结二进制和提交时来源快照。

## 时间开销：两套计时基准不能混用

| 体系 | raw / min | IP-only / min | IP+Olsen / min | Olsen 相对 IP 的增量 | 计时基准 |
|---|---:|---:|---:|---:|---|
| C2 300万步 | 351.7753 | 278.5523 | 287.6777 | +3.2760% | 三组第300万步的累计报告时刻时间 |
| N2 260万步 | 238.0844 | 286.8898 | 281.2408 | −1.9691% | 三组完整 benchmark solve timer |

定义为 `(T_IP+Olsen-T_IP)/T_IP`。C2 额外约547.52秒 / 9.13分钟；N2 差值约−338.94秒 / −5.65分钟。C2 不用360万步的最终 raw/joint timer 与300万步 IP clock 相除，也不把报告时刻累计计时伪装成完整进程计时。

直接报告 correction timer：C2 300万步累计629.823秒；N2 260万步累计665.602秒。新版报告点 timer 包含延迟刷新，不能与旧版仅即时维护后的扫描 timer 当作完全相同指标比较。它仍不包含全部迭代内维护成本；端到端差值另行报告。C2 IP report timer 和刷新计数未写入 stdout，保持未知（null），**不补零**。

## 达到同一目标精度的时间

以每10000步的离散报告点定义首次达标，不插值。时间包含当次求解与修正；参考不是本轮末点。

| 体系 | 误差 / Ha | raw / s | IP-only / s | IP+Olsen / s | raw → joint | IP → joint |
|---|---:|---:|---:|---:|---:|---:|
| C2 | 1e-4 | 342.626 | 64.46 | 73.163 | 4.683× | 0.881× |
| C2 | 1e-5 | 1585.223 | 315.66 | 201.839 | 7.854× | 1.564× |
| C2 | 1e-6 | 5753.214 | 1349.51 | 706.207 | 8.147× | 1.911× |
| N2 | 1e-4 | 566.308 | 72.907 | 69.734 | 8.121× | 1.045× |
| N2 | 1e-5 | 3069.133 | 415.030 | 136.374 | 22.505× | 3.043× |
| N2 | 1e-6 | 11597.061 | 2675.124 | 829.166 | 13.986× | 3.226× |

这些是本轮各组的 time-to-target，不是旧/新代码的受控加速比。1e-4 Ha 时 IP 和联合修正都在首个报告点达标，采样太粗，不能过度解释小时间差。

## 精度和图表限制

长期参考不重跑：C2 −75.7319603747 Ha，N2 −109.2821730115 Ha。来源见两体系 reference.json。误差图只展示 **>=1e-7 Ha** 的 reference-resolved 区域；低于参考的原始点仍在数据中，不用于声称更高精度。参考末尾漂移不是误差的严格上界。

C2 IP 从原同一作业的 stdout 恢复，能量打印10位小数（舍入半单位约5e-11 Ha）、时间打印2位小数（约0.005秒），波函数计数为打印整数。解析器检查缺失/重复报告、开关对应的修正分量、raw JSON/log 一致性和时钟单调性；300个报告点完整。这个精度足以支持当前绘图区间，但不能当作原未舍入的 solver JSON。

- [C2 三组 A/B/C 图](../../experiment_results/energy_correction/1022419/c2_ccpvdz/prefix_3000000/05_three_way_reference_resolved.pdf)：raw能量取独立baseline，三组横轴用自身的非零系数数量；Panel C为直接 report correction time / joint time。
- [C2 Olsen 时间图](../../experiment_results/energy_correction/1022419/c2_ccpvdz/prefix_3000000/06_ip_olsen_timing.pdf)：同allocation的IP/joint累计时间及端到端相对增量。
- [N2 A/B/C 图](../../experiment_results/energy_correction/1022420/n2_ccpvdz/05_correction_advantage_reference_resolved.pdf)：沿用现有 summarize 程序，raw能量取joint中的未修正能量、raw时间取独立baseline；目标精度首次跨越已用独立baseline核对一致。Panel C是 `(T_joint-T_raw)/T_joint`，不是C2的直接timer比例，也不是Olsen对IP的比例。
- [N2 Olsen 时间图](../../experiment_results/energy_correction/1022420/n2_ccpvdz/06_ip_olsen_timing.pdf)：提供同allocation的IP/joint对比。泛用06 summary的跨作业说明不表示本轮跨节点；准确来源见manifest/runtime summary。

两套A/B/C的Panel C口径目前不同，**论文不能不改图注就把它们并排当作同一个量**；JSON/CSV都保留重画需要的数据。非零系数数目与存储波函数/残差条目数也不是同一尺度。

## 无需集群或积分文件的重分析

从仓库根目录运行，输出到新的本地目录，避免改动正式结果。只需 Python、NumPy、Matplotlib；无FCIDUMP、Slurm、二进制依赖。

```bash
python3 experiments/energy_correction/analyze_prefix_from_log.py \
  experiment_results/energy_correction/1022419/c2_ccpvdz/prefix_3000000/source_paired_prefix.json \
  experiment_results/energy_correction/1022419/c2_ccpvdz/prefix_3000000/source_stdout_prefix.out \
  experiment_results/energy_correction/reanalysis_c2_3000000 --iterations 3000000
python3 experiments/energy_correction/summarize_large_scale.py \
  experiment_results/energy_correction/1022420/n2_ccpvdz/05_large_scale_with_timings.json \
  experiment_results/energy_correction/reanalysis_n2_2600000
python3 experiments/energy_correction/compare_ip_olsen.py \
  experiment_results/energy_correction/1022420/n2_ccpvdz/05_ip_only.json \
  experiment_results/energy_correction/1022420/n2_ccpvdz/05_large_scale_with_timings.json \
  experiment_results/energy_correction/reanalysis_n2_2600000
```

提交前验证两体系重分析的JSON/CSV与正式结果一致、所有公布源片段不超过选中步数、原始文件SHA256和冻结二进制来源一致，以及无补跑/GPU节点计时混入。本次不修改求解器或 correction 算法，不新增计算作业，不上传积分文件、构建目录或可执行文件。
