# H2O/cc-pVDZ：一次 500 万步的 IP+Olsen 运行

本次目标是判断能量修正是否高估/低估，并得到同一个波函数上的三条完整曲线。
只启动一次 `olsen-only` benchmark（该历史模式名实际同时启用 IP 和 Olsen），
固定 5,000,000 个 block，每 block 64 coordinates，每 10,000 步报告。
不做校准、额外预热、独立 raw/IP 求解、性能对照或重复配对；不按预算缩短步数。

## 参数和输出

- partition / bigMem3，64 OpenMP 线程、128G、18 小时，非独占，spread/cores。
  64 线程不自动等于 64 个独立物理核。
- `z_threshold=0`、capacity=848256040、`stopping_dx_threshold=0`。
- 每个报告点同时保留 `E`、`E+IP`、`E+IP+Olsen`，以及各自有符号误差。
  三者共用本次联合运行的时间，不能解释为三种方法独立的求解耗时。
- 参考为用户确认已有严格精度认证的
  `/home/yjzhang/codes/cdfci/paper-jctc/h2o_ccpvdz/h2o_ccpvdz_new.out`，
  能量 −76.2418600697 Ha。参考文件随源码冻结并校验 SHA256。
  不再使用任意的 1e-7 Ha 未分辨区；图保留正负号及全部点。
- 外部 IP 保留现有增量分母更新实现，Olsen 使用报告时的当前 Rayleigh 能量。
  没有修改求解器、修正公式或坐标更新。

结果默认保存至 `experiment_results/energy_correction/h2o_ccpvdz/<新作业号>/`：

- `manifest.json`、`plan.json`、`command.json`：来源、固定步数、精确执行参数；
- `runs/joint.json.jsonl`：逐报告刷新的完整精度原始记录；
- `three_way_trajectory.csv`：同波函数三种能量、修正项、有符号误差及计数；
- `target_tolerances.csv`：各误差阈值的首次达标报告点（不作时间插值）；
- `summary.json`、`completion.json`：当前进度、低于参考的点数及完成状态；
- `signed_errors.png/pdf`、`errors_vs_work.png/pdf`：每十个报告点及结束时更新。

固定 500 万步是目标，不承诺在 18 小时内一定完成。如果达到作业时限，Slurm 提前
10 分钟发 USR1，驱动器终止自己的 srun 进程组并保留部分曲线和日志，不冒充完整结果。

## 本次提交：由用户执行

已准备的单次实验快照路径为 `build_h2o_experiment/frozen-h2o-single-20261019`。
旧任务 1027222 是九次求解的性能对照方案，应先取消；其已有日志和结果继续保留。
Codex 没有确认旧任务取消，也不代用户提交新任务。

```bash
cd /home/yjzhang/codes/cdfci/perturbation/CDFCI
scancel 1027222
env -u CDFCI_FROZEN_RUN -u CDFCI_EXPERIMENT_OUTPUT \
  sbatch experiments/energy_correction/run_h2o_ccpvdz.sbatch
```

脚本默认选择上述单次实验快照。清除两个旧环境变量，避免误选之前的快照或结果目录。
新日志为 `cdfci_h2o_5m_<新作业号>.out`。不需要使用另一个账户。

## 重新准备快照（本次已有快照时不必执行）

以下是从当前源码重新构建的流程；换用一个尚不存在的新 FROZEN 目录：

```bash
python3 experiments/energy_correction/h2o_single_run.py prepare \
  --root "$PWD" --frozen FROZEN --eigen build_h2o_experiment/_deps/eigen-src \
  --reference-file /home/yjzhang/codes/cdfci/paper-jctc/h2o_ccpvdz/h2o_ccpvdz_new.out
cmake -S FROZEN/src -B FROZEN/build -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON -DCDFCI_BUILD_EXPERIMENTS=ON -DCDFCI_ENABLE_OPENMP=ON \
  -DFETCHCONTENT_SOURCE_DIR_EIGEN="$PWD/FROZEN/eigen"
cmake --build FROZEN/build --target cdfci_energy_correction_benchmark_omp --parallel 4
ctest --test-dir FROZEN/build \
  -R '^(cdfci.h2o_(batch_paths|single_run)|cdfci_energy_correction_benchmark_omp.output)$' \
  --output-on-failure
python3 FROZEN/src/experiments/energy_correction/h2o_single_run.py seal --frozen FROZEN
sbatch --export=ALL,CDFCI_FROZEN_RUN="$PWD/FROZEN" \
  experiments/energy_correction/run_h2o_ccpvdz.sbatch
```

可在结束、运行中或中断后重新生成分析：

```bash
python3 build_h2o_experiment/frozen-h2o-single-20261019/src/experiments/energy_correction/h2o_single_run.py \
  analyze --output experiment_results/energy_correction/h2o_ccpvdz/新作业号
```

旧驱动 `h2o_experiment.py` 保留用于旧任务的恢复与分析，当前 sbatch 不再调用其九次求解流程。
