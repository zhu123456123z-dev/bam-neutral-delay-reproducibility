# BAM 神经网络论文代码与数值数据

本仓库包含论文算例的 MATLAB 主计算、Python 独立求解器复核、参数、已有数值数据以及五幅图。论文正文、老师的论文、个人路径和私人文件不在发布范围内。

## 运行

建议使用 Python 3.10--3.12；本次仓库核查使用 Python 3.12。

```sh
python -m pip install -r requirements.txt
python run_reproduction.py
```

默认命令只核对已有数组和矩阵判据，不重新运行积分，报告生成在 `outputs/`。

```sh
python run_reproduction.py --simulate
```

此命令重新进行 Python 模拟、容差/回溯检查、绘图，并与已保存的 MATLAB 数组比较。新输出放在 `outputs/python/`。

在 MATLAB 中把当前目录切换到仓库根目录，输入 `RUN_ME`，可重新计算并导出五幅图，新结果放在 `outputs/reproduced_matlab/`。归档 MATLAB 数据实际由 R2026b 的 `dde23` 生成；不声称已在 R2021b 上完成运行。本次整理未重新运行 MATLAB。

`04_results/` 保存原始归档数据与补齐的报告，不会被运行入口覆盖。`simulation.mat` 和 `simulation_python.npz` 包含三组完整状态轨迹及导数；CSV 是便于查看的部分字段导出，具体列名和数组维度见 `docs/DATA_DICTIONARY.md`。

注意：论文图4是 `fig5_derivative_errors`，论文图5是 `fig4_regions`。文件名沿用原始导出编号，并非排错。

指数虚线只是斜率参照，不是定理完整误差上界；不同求解器和加严容差之间的差异不是严格认证的真误差界。不同平台或软件版本可能产生微小浮点差异，不保证图片像素完全一致。

英文说明包含完整文件清单及复现指令，见 [README.md](README.md)。尚未由权利人选定开源许可证；公开可访问不等同于已授予任意再利用许可。
