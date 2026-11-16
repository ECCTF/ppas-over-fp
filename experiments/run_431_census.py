"""Reproduce the complete p=431 census and its independent checks.

Requires Python 3 and NumPy. Every stage uses the same interpreter.
The complete log is retained next to the generated results.
"""
from pathlib import Path
import subprocess
import sys

BASE = Path(__file__).resolve().parent
OUT = BASE / 'results'
STAGES = [
    ['split_gluing_seeds.py', '431'],
    ['restriction_gluing_seeds.py', '431'],
    ['rational_census.py', '431'],
    ['assemble_rational_census.py', '431'],
    ['verify_gluing_seeds.py', '431'],
    ['verify_rational_census.py', 'results/p431_rational_census.json'],
    ['verify_431_arithmetic.py'],
]


def main():
    OUT.mkdir(exist_ok=True)
    with (OUT / 'reproduction_p431.log').open('w', encoding='utf8') as log:
        for stage in STAGES:
            title = 'Running ' + ' '.join(stage)
            print(title, flush=True)
            log.write(title + '\n')
            process = subprocess.Popen(
                [sys.executable, '-X', 'utf8', *stage], cwd=BASE,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, encoding='utf8',
            )
            for line in process.stdout:
                log.write(line)
                log.flush()
                print(line, end='', flush=True)
            if process.wait():
                raise RuntimeError(f'Failed stage: {stage[0]}')
        print('Complete p=431 census and all independent checks passed.', flush=True)
        log.write('Complete p=431 census and all independent checks passed.\n')


if __name__ == '__main__':
    main()
