# BAM 神经网络论文代码与数值数据

本仓库提供论文数值算例的 MATLAB 程序、Python 独立求解器、参数、计算结果和五幅图。

## 运行

测试环境为 Python 3.12.14、NumPy 1.26.4、SciPy 1.13.1 和 Matplotlib 3.9.2。

```sh
python -m pip install -r requirements.txt
python run_reproduction.py
python verify_archive.py
```

`run_reproduction.py` 核对矩阵判据并比较已有 MATLAB 与 Python 数组，报告生成在 `outputs/`。`verify_archive.py` 检查归档数据和文件校验值。这两条命令不重新进行数值积分。

```sh
python run_reproduction.py --simulate
```

此命令重新进行 Python 模拟、容差/回溯检查、绘图，并与已保存的 MATLAB 数组比较。新输出放在 `outputs/python/`。

在 MATLAB 中把当前目录切换到仓库根目录，输入 `RUN_ME`，可重新计算并导出五幅图，新结果放在 `outputs/reproduced_matlab/`。归档 MATLAB 数据由 R2026b（26.2.0.3386108）的 `dde23` 生成。本次发布未测试其他 MATLAB 版本。

`04_results/` 保存归档数据和复核报告，运行入口不会覆盖该目录。`simulation.mat` 和 `simulation_python.npz` 包含三组状态轨迹及导数；CSV 导出其中的部分字段，列名和数组维度见 `docs/DATA_DICTIONARY.md`。

文件名沿用调整图序前的编号：论文图4对应 `fig5_derivative_errors`，论文图5对应 `fig4_regions`。

指数虚线仅表示衰减斜率，不是定理中的完整误差上界。求解器对照和容差加严得到的是数值解之间的差异，不是严格认证的真误差界。不同平台或软件版本可能产生微小浮点差异。

完整文件清单及复现指令见 [README.md](README.md)。引用时请注明仓库地址和使用的提交版本。目前尚未指定再利用许可证。

## AI 辅助说明

AI 工具用于辅助代码开发与调试、语言润色与翻译，以及数学表述和数值结果的检查。作者对研究内容承担责任。
