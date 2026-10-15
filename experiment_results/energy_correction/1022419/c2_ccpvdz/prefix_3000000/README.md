# C2：原 bigMem5 同一作业的 300 万步完整前缀

按用户要求，本轮 C2 比较只采用任务 `1022419` 的前 3000000 步。
原任务在最后 IP-only 求解时超时，但 stdout 已完整打印至 3470000 步。
因此三组各有 300 个同节点、同 allocation、同参数的有效报告点，不需要跨节点补跑。
原任务仍应标为 TIMEOUT；“完整”只指所选前缀。原始文件未删除或截断。

| 变体 | 300 万步累计报告时间 / min |
|---|---:|
| raw | 351.7753 |
| IP-only | 278.5523 |
| IP+Olsen | 287.6777 |

Olsen 相对 IP 的时间增幅为 **3.2760%**。上述计时都是各组第 300 万步的报告时刻累计时间，不是完整 benchmark 进程的最终 solve timer。
固定执行顺序和共享节点负载仍影响时间差；不是重复测量的置信区间结果。

能量和计时来源：

- raw 使用 paired JSON 的独立 `baseline_variational_energy` 和 `raw_wall_seconds`，非联合运行的未修正能量。
- IP+Olsen 使用原始 paired JSON 的完整精度数据。
- IP-only 从原 stdout 的独立 IP/PT2 段恢复：能量打印 10 位小数、时间打印 2 位小数，分别有约 5e-11 Ha / 0.005 s 的舍入半单位。
- 波函数计数为各变体自身数据；IP report timer 和刷新计数未记录，保持未知，不填零。

`prefix_summary.json` 记录来源 SHA256、精度和计时口径；`prefix_trajectory.json` / CSV 保存三组逐点数据；`target_tolerances.csv` 保存目标精度交叉点。
`05_three_way_reference_resolved.pdf` / PNG 展示三组误差对大小、误差对时间和直接修正报告时间占比；误差仅画 >=1e-7 Ha 区域。
`06_ip_olsen_timing.pdf` / PNG 为同节点 IP/Olsen 时间消融。
直接 correction timer 包含报告点延迟刷新及修正计算，不包含全部迭代内维护成本，不能与端到端差值混用。

`source_paired_prefix.json` 和 `source_stdout_prefix.out` 只含 <=300 万步的选中报告，是可独立重分析的发布源片段；原360万步申请参数/日志段标记作为来源事实保留，最终全长 solve summary 不纳入导出。
原全长源文件的 SHA256 记录在 source_export 元数据，文件本身已移到仓库外可恢复归档，不在本次发布中。

使用 `experiments/energy_correction/analyze_prefix_from_log.py` 将上述两个源片段传入即可重分析；完整命令见 `experiments/energy_correction/RESULTS_DEFERRED.md`。
bigMem3 补跑 `1024227` 已在约 2 分 13 秒时取消，其脚本、输出和日志已从仓库工作树移出，不用于这里的对照。
