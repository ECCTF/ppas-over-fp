"""Package the published sources and data without build products or Git metadata."""
from pathlib import Path
import zipfile,json,csv
BASE=Path(__file__).resolve().parent;OUT=BASE/'results';ROOT=BASE.parent
files=[ROOT/name for name in ['README.md','requirements.txt','environment.yml','.gitignore','.gitattributes']]
files+=list((ROOT/'docs').glob('*'))+list((ROOT/'.github'/'workflows').glob('*.yml'))
files+=[f for f in BASE.iterdir() if f.suffix in {'.py','.sage','.cpp','.sh'}]
files+=[f for f in OUT.iterdir() if f.suffix in {'.json','.csv','.svg','.png','.html','.txt'} and not f.name.endswith('_maximal_orbits.csv')]
files.append(ROOT/'figures'/'.gitkeep')
for p in [47,71,127]:
    g=json.loads((OUT/f'p{p}_rational_graph.json').read_text());n=len(g['vertices'])
    for name,expected in [('all_vertices',n),('spanning_tree',n-1),('all_edges',len(g['edges']))]:
        with (OUT/f'p{p}_{name}.csv').open(encoding='utf-8-sig') as f:
            assert len(list(csv.DictReader(f)))==expected
build=ROOT/'build';build.mkdir(exist_ok=True)
dest=build/'ppas-over-fp.zip'
with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED) as z:
    for f in sorted(set(files)):
        if f.is_file():z.write(f,Path('ppas-over-fp')/f.relative_to(ROOT))
with zipfile.ZipFile(dest) as z:assert z.testzip() is None
print(dest,'bytes',dest.stat().st_size,flush=True)
