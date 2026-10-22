# BAM neural networks with neutral distributed delays

Numerical reproducibility materials for **Weighted Pseudo-Almost Automorphic Solutions and Stability of High-Order BAM Neural Networks with Neutral Distributed Delays**.

This repository contains the MATLAB simulation and plotting code, an independent Python method-of-steps implementation, model parameters, archived numerical data, and the five figure exports. It does not contain the manuscript or third-party papers. No publication or acceptance is implied.

## Quick start

Use Python 3.10--3.12 (repository verification uses Python 3.12):

```sh
python -m pip install -r requirements.txt
python run_reproduction.py
python verify_archive.py
```

The default command verifies the matrix conditions and compares the **saved** MATLAB and Python arrays. It writes new reports to `outputs/`. It does not rerun either solver.

To recompute the independent Python simulations, convergence checks, and plots:

```sh
python run_reproduction.py --simulate
```

New Python data and plots are written to `outputs/python/`, and are compared with the archived MATLAB data. Runtime depends on the machine. MATLAB is not required for these Python commands.

## MATLAB reproduction

In MATLAB, set the current folder to this repository root and run:

```matlab
RUN_ME
```

This recomputes the `dde23` simulations, tighter-tolerance and longer-history checks, matrix conditions, and all five plots. New results go to `outputs/reproduced_matlab/`; archived data are never overwritten. The saved MATLAB results were generated with **MATLAB 26.2.0.3386108 (R2026b), `dde23`**. The code uses base MATLAB functions and was originally written for the R2021b-era API; a completed R2021b run is not claimed. A new MATLAB run has not been performed during the repository packaging step.

To compare a new MATLAB run with the archived Python run:

```sh
python 03_code/compare_solvers.py --matlab outputs/reproduced_matlab
```

To redraw an existing MATLAB run without repeating integration (this rewrites its figure exports):

```matlab
addpath('03_code');
plot_bam_results(fullfile(pwd, 'outputs', 'reproduced_matlab'));
```

## Contents

| Location | Contents |
| --- | --- |
| `03_code/parameters.json` | All model parameters, constant histories, tolerances and output grid |
| `03_code/run_bam.m` | MATLAB solver, memory-variable formulation and numerical checks |
| `03_code/plot_bam_results.m` | MATLAB PNG (300 dpi), vector PDF and editable FIG exports |
| `03_code/simulate_bam.py` | Independent SciPy DOP853 method-of-steps solver |
| `03_code/check_conditions.py` | Table 1, matrix certificates, decay constant and threshold |
| `03_code/compare_solvers.py` | Saved-array cross-solver comparison and matrix-certificate checks |
| `03_code/make_regions.py` | Alternative Python rendering of the analytic criterion regions |
| `04_results/matlab_derivative/` | Archived MATLAB arrays, CSV, validation records and all five figures |
| `04_results/python/` | Archived Python arrays, CSV and original validation record |
| `04_results/reports/` | Reports recomputed from the packaged parameters and saved arrays |
| `docs/DATA_DICTIONARY.md` | Array shapes, column meanings and figure-number mapping |
| `docs/PROVENANCE.md` | Source provenance and packaging changes |
| `SHA256SUMS.txt` | Checksums of the release files, excluding this checksum file |
| `verify_archive.py` | Integrity, shapes, parameters, figure-presence and paper-value checks |

## Figure numbering

The original export names are preserved to maintain provenance. **The numeric filename prefix is not always the paper figure number.**

| Paper figure | Archived export basename |
| --- | --- |
| 1: reference states | `fig1_states` |
| 2: different constant histories | `fig2_histories` |
| 3: state errors | `fig3_errors` |
| 4: derivative errors | `fig5_derivative_errors` |
| 5: criterion regions | `fig4_regions` |

## Interpretation and limitations

The relative/absolute tolerances are `1e-8`/`1e-10`; the adaptive step-size cap and output-grid spacing are both `0.02`. References start at `-80`, with a longer-history check starting at `-120`. Three constant histories start at zero. Derivatives are evaluated from the model right-hand side, not finite differences. The plotted dashed exponentials are slope guides, **not** complete theoretical error bounds. Cross-solver and refinement differences are empirical diagnostics, not interval-certified bounds on the true error. Figure 5 shows analytic sufficient-condition regions, not an experimentally measured stability region.

The original Python validation record reports SciPy 1.13.1. New runs may differ slightly across platforms, MATLAB/Python versions, and floating-point libraries; exact pixel matching is not promised. See the archived validation files and freshly generated comparison reports rather than comparing screenshots alone.

## Citation and reuse

Until publication, cite the repository URL and the exact commit used. No author names, journal publication details, DOI, or release DOI have been invented. A software/data license has not yet been selected by the rights holders; public availability is not an unrestricted reuse license.

## AI-assisted preparation

AI assistance was used for code development/debugging, language editing/translation, and checking mathematical exposition and numerical results. The repository-packaging changes concern paths, documentation, verification and distribution; they do not constitute independent certification of the mathematical proofs. Authors remain responsible for the research and the manuscript's full AI disclosure.

中文运行说明见 [README_zh.md](README_zh.md)。
