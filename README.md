# BAM neural networks with neutral distributed delays

Numerical reproducibility materials for **Weighted Pseudo-Almost Automorphic Solutions and Stability of High-Order BAM Neural Networks with Neutral Distributed Delays**.

MATLAB and Python implementations of the numerical example, with parameters, saved results, and figure files.

## Quick start

Tested with Python 3.12.14, NumPy 1.26.4, SciPy 1.13.1, and Matplotlib 3.9.2.

```sh
python -m pip install -r requirements.txt
python run_reproduction.py
python verify_archive.py
```

`run_reproduction.py` checks the matrix conditions and compares the saved MATLAB and Python arrays. Reports are written to `outputs/`. `verify_archive.py` checks the archived data and file checksums. Neither command reruns the simulations.

To recompute the independent Python simulations, convergence checks, and plots:

```sh
python run_reproduction.py --simulate
```

The `--simulate` option writes new data and plots to `outputs/python/` and compares them with the saved MATLAB results. The Python workflow does not require MATLAB.

## MATLAB reproduction

In MATLAB, set the current folder to this repository root and run:

```matlab
RUN_ME
```

This runs `dde23`, performs the tolerance and history-horizon checks, and exports all five figures to `outputs/reproduced_matlab/`. It leaves `04_results/` unchanged. The archived MATLAB results were computed with MATLAB R2026b (26.2.0.3386108), using base MATLAB functions. Other MATLAB versions have not been tested for this release.

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
| `03_code/parameters.json` | Model parameters, constant histories, tolerances and output grid |
| `03_code/run_bam.m` | MATLAB solver, memory-variable formulation and numerical checks |
| `03_code/plot_bam_results.m` | MATLAB PNG (300 dpi), vector PDF and editable FIG exports |
| `03_code/simulate_bam.py` | Independent SciPy DOP853 method-of-steps solver |
| `03_code/check_conditions.py` | Table 1, matrix certificates, decay constant and threshold |
| `03_code/compare_solvers.py` | Saved-array cross-solver comparison and matrix-certificate checks |
| `03_code/make_regions.py` | Alternative Python rendering of the analytic criterion regions |
| `04_results/matlab_derivative/` | Archived MATLAB arrays, CSV, validation records and all five figures |
| `04_results/python/` | Archived Python arrays, CSV and original validation record |
| `04_results/reports/` | Matrix and cross-solver reports, plus the Python rerun record |
| `docs/DATA_DICTIONARY.md` | Array shapes, column meanings and figure-number mapping |
| `docs/PROVENANCE.md` | Run history and data provenance |
| `SHA256SUMS.txt` | Checksums of the release files, excluding this checksum file |
| `verify_archive.py` | Integrity, shapes, parameters, figure-presence and paper-value checks |

## Figure numbering

The export names predate the final figure order. Use the following correspondence when matching files to the paper.

| Paper figure | Archived export basename |
| --- | --- |
| 1: reference states | `fig1_states` |
| 2: different constant histories | `fig2_histories` |
| 3: state errors | `fig3_errors` |
| 4: derivative errors | `fig5_derivative_errors` |
| 5: criterion regions | `fig4_regions` |

## Interpretation and limitations

The relative/absolute tolerances are `1e-8`/`1e-10`; the adaptive step-size cap and output-grid spacing are both `0.02`. References start at `-80`, with a longer-history check starting at `-120`. Three constant histories start at zero. Derivatives are evaluated from the model right-hand side, not finite differences. The plotted dashed exponentials are slope guides, **not** complete theoretical error bounds. Cross-solver and refinement differences are empirical diagnostics, not interval-certified bounds on the true error. Figure 5 shows analytic sufficient-condition regions, not an experimentally measured stability region.

The original Python run used SciPy 1.13.1. Small numerical differences can arise across platforms and software versions. The validation records and comparison reports give the relevant numerical differences.

## Citation and reuse

Cite the repository URL and the commit used. A reuse license has not yet been assigned.

## AI assistance

AI tools assisted with code development and debugging, language editing and translation, and checks of mathematical exposition and numerical results. The authors retain responsibility for the work.

