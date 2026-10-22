"""Check archived or newly computed MATLAB/Python arrays without modifying them."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matlab', type=Path, default=ROOT/'04_results/matlab_derivative')
    parser.add_argument('--python', type=Path, default=ROOT/'04_results/python')
    parser.add_argument('--output', type=Path, default=ROOT/'outputs/cross_solver_validation.json')
    args = parser.parse_args()
    with np.load(args.python/'simulation_python.npz', allow_pickle=False) as py:
        mat = loadmat(args.matlab/'simulation.mat')
        np.testing.assert_allclose(mat['t'].ravel(), py['t'], rtol=0, atol=1e-12)
        histories = np.moveaxis(mat['trajectories'], 2, 0)
        pairs = {
            'reference_state_max_difference': (mat['zr'][:4], py['zr'][:4]),
            'trajectory_state_max_difference': (histories, py['trajectories']),
            'state_error_curve_max_difference': (mat['errors'], py['errors']),
            'reference_derivative_max_difference': (mat['reference_derivatives'], py['reference_derivatives']),
            'trajectory_derivative_max_difference': (np.moveaxis(mat['trajectory_derivatives'],2,0), py['trajectory_derivatives']),
            'derivative_error_curve_max_difference': (mat['derivative_errors'], py['derivative_errors']),
        }
        report = {'grid_points': int(py['t'].size),
                  'comparison': 'dde23 versus DOP853 method of steps; empirical agreement, not a certified true-error bound'}
        for name, (left, right) in pairs.items():
            if left.shape != right.shape or not (np.isfinite(left).all() and np.isfinite(right).all()):
                raise ValueError(f'Invalid arrays for {name}')
            report[name] = float(np.max(np.abs(left-right)))
        if report['reference_state_max_difference'] >= 1e-7 or report['trajectory_state_max_difference'] >= 1e-7:
            raise ValueError('Cross-solver state discrepancy exceeds the verification threshold 1e-7')
    checks = json.loads((args.matlab/'condition_matlab.json').read_text(encoding='utf-8'))
    for c in checks:
        np.testing.assert_allclose(c['rho_J'], np.sqrt(.105125*c['eta']), rtol=0, atol=1e-12)
        np.testing.assert_allclose(c['rho_G'], np.sqrt(.92*c['eta']), rtol=0, atol=1e-12)
        if min(c['residual']) <= 0:
            raise ValueError('Nonpositive matrix-certificate residual')
    report['matlab_matrix_certificates_passed'] = True
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, indent=2))

if __name__ == '__main__':
    main()
