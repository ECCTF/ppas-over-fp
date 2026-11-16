"""Reproduce the ten-bit experiments, independent audit, and figures.

Run with Python 3 and NumPy, Matplotlib, and ReportLab installed.
The exact restriction computation itself uses only the standard library.
This entry point does not rerun the older exhaustive small-prime census.
"""
from pathlib import Path
import subprocess
import sys

BASE=Path(__file__).resolve().parent
OUT=BASE/'results'
STEPS=[
    ('ten_bit_experiments.py','ten_bit_restriction_run.log'),
    ('ten_bit_jacobian_samples.py','ten_bit_jacobian_run.log'),
    ('verify_ten_bit_results.py','ten_bit_independent_verification.log'),
    ('plot_ten_bit.py','ten_bit_restriction_plot.log'),
    ('plot_ten_bit_jacobian_samples.py','ten_bit_jacobian_plot.log'),
]

if __name__=='__main__':
    OUT.mkdir(exist_ok=True)
    for script,log in STEPS:
        print('Running',script,flush=True)
        with (OUT/log).open('w',encoding='utf8') as stream:
            subprocess.run([sys.executable,str(BASE/script)],cwd=BASE,
                           stdout=stream,stderr=subprocess.STDOUT,check=True)
        print('Passed',script,flush=True)
    print('All ten-bit computations and figures reproduced.',flush=True)
