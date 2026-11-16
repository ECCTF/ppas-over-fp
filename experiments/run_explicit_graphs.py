"""Build the three complete graphs; no class-group connectivity test is used."""
from pathlib import Path
import subprocess,sys
BASE=Path(__file__).resolve().parent
for p in map(int,sys.argv[1:] or [47,71,127]):
    with (BASE/'results'/f'p{p}_explicit_build.log').open('w') as log:
        for stage in [
            ['split_gluing_seeds.py',str(p)],
            ['restriction_gluing_seeds.py',str(p)],
            ['rational_census.py',str(p)],
            ['assemble_rational_census.py',str(p)],
            ['verify_gluing_seeds.py',str(p)],
            ['verify_rational_census.py',f'results/p{p}_rational_census.json'],
        ]:
            print('RUN',*stage,flush=True)
            child=subprocess.Popen([sys.executable,'-u',*stage],cwd=BASE,
                stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
            for line in child.stdout:
                log.write(line);log.flush();print(line,end='',flush=True)
            if child.wait():raise RuntimeError(stage)
