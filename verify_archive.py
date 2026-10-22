"""Check archived data, metadata, release checksums and selected paper values."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parent

def main():
    matdir = ROOT/'04_results/matlab_derivative'
    mat = loadmat(matdir/'simulation.mat', simplify_cells=True)
    params = json.loads((ROOT/'03_code/parameters.json').read_text())
    for key, expected in params.items():
        np.testing.assert_allclose(np.asarray(mat['p'][key]), np.asarray(expected))
    with np.load(ROOT/'04_results/python/simulation_python.npz', allow_pickle=False) as py:
        np.testing.assert_allclose(mat['t'], py['t'], rtol=0, atol=1e-12)
        assert py['trajectories'].shape == (3,4,2001)
        assert mat['trajectories'].shape == (4,2001,3)
        for data in (mat, py):
            for key in ['t','zr','trajectories','errors','reference_derivatives',
                        'trajectory_derivatives','derivative_errors','joint_errors']:
                assert np.isfinite(data[key]).all(), key
            np.testing.assert_allclose(data['joint_errors'], np.maximum(data['errors'], data['derivative_errors']))
        np.testing.assert_allclose(py['errors'], np.max(np.abs(py['trajectories']-py['zr'][:4]),axis=1))
        np.testing.assert_allclose(mat['errors'],np.max(np.abs(mat['trajectories']-mat['zr'][:4,:,None]),axis=0).T)
    checks = json.loads((ROOT/'04_results/reports/condition_report.json').read_text())
    np.testing.assert_allclose(checks['2.0']['M'],13.121513,rtol=0,atol=5e-7)
    np.testing.assert_allclose(checks['2.0']['lambda_threshold'],0.819437,rtol=0,atol=5e-7)
    for name in ['fig1_states','fig2_histories','fig3_errors','fig4_regions','fig5_derivative_errors']:
        for ext in ['png','pdf','fig']:
            assert (matdir/f'{name}.{ext}').stat().st_size > 0
    count = 0
    for line in (ROOT/'SHA256SUMS.txt').read_text().splitlines():
        digest, name = line.split('  ',1)
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
        count += 1
    print(f'ARCHIVE_VERIFIED: arrays, parameters, figures, paper values and {count} checksums')

if __name__ == '__main__':
    main()
