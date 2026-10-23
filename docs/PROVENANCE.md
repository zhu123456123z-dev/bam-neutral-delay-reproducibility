# Run history and data provenance

The MATLAB data were generated on 2026-10-15. The repository was assembled on 2026-10-22 and its documentation revised on 2026-10-23.

The archived MATLAB results were computed in MATLAB Online R2026b (26.2.0.3386108) with `dde23`. The Python results were computed independently with a DOP853 method-of-steps solver using SciPy 1.13.1. The original validation records are stored beside the data.

The numerical archives are unchanged. For distribution, output paths were made relative to the repository and new runs were directed to `outputs/`. The model, histories and solver tolerances were retained. `make_regions.py` provides an alternative Python rendering of the analytic criterion regions; the archived figure was exported from MATLAB.

The matrix and cross-solver reports in `04_results/reports/` were recomputed from the parameters and saved arrays during repository preparation. `packaging_python_validation.json` records a separate Python rerun, which reproduced all archived NPZ arrays exactly in the test environment. The MATLAB simulation was not rerun at that stage.

`SHA256SUMS.txt` records the release file checksums. It is updated when source comments or documentation change; the archived numerical data retain their original checksums.
