# Provenance and packaging changes

Prepared for public repository distribution on 2026-10-22 (Asia/Shanghai).

The archived MATLAB outputs were returned from a MATLAB Online R2026b run on 2026-10-15. The original validation JSON records MATLAB 26.2.0.3386108 and `dde23`. The MATLAB driver, plotting helper and parameters were taken from the accompanying reproduction package. The Python source and archived derivative-inclusive arrays come from the independent DOP853 method-of-steps calculation; its original record reports SciPy 1.13.1.

The numerical archive files are copied without changing their bytes. Source packaging changes redirect newly generated results into `outputs/`, replace machine-specific paths, provide a cross-solver comparison entry point, and add documentation and checksums. The numerical model, history functions, solver tolerances and archived results are not adjusted to improve agreement. `make_regions.py` is an alternative Python plot of analytic inequalities; the archived criterion-region figure is the MATLAB export.

The two missing summary reports are reconstructed from the parameter file and saved MATLAB/Python arrays. They are packaging-time reconstructions, not claimed to be original solver log files. A new Python simulation is a separate verification run; original Python records are preserved. No new MATLAB simulation is claimed during packaging.

Manuscripts, third-party publications, Word-building utilities, browser/authentication information, local absolute paths, chat records and unrelated personal files are deliberately excluded. Research software and numerical data are published, not the entire working directory.
