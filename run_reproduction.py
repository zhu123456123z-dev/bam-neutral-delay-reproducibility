"""Verify the archived paper data; optionally rerun the independent Python solver."""
import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--simulate', action='store_true', help='Recompute Python responses and error checks before comparing with archived MATLAB data')
    args = parser.parse_args()
    def run(name, *extra):
        subprocess.run([sys.executable, str(ROOT/'03_code'/name), *map(str,extra)], cwd=ROOT, check=True)
    run('check_conditions.py')
    if args.simulate:
        run('simulate_bam.py')
        run('make_regions.py')
        run('compare_solvers.py', '--python', ROOT/'outputs/python')
    else:
        run('compare_solvers.py')
    print('BAM_REPRODUCTION_CHECKS_PASSED')

if __name__ == '__main__':
    main()
