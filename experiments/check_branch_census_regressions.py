"""Check new exhaustive enumeration against the four older full censuses."""
from pathlib import Path
import argparse,csv,json,os,subprocess,tempfile
from rational_census import Census
BASE=Path(__file__).resolve().parent;OUT=BASE/'results'
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--enumerator',type=Path,
                    default=BASE.parent/'build'/('ppas_branch_census.exe' if os.name=='nt' else 'ppas_branch_census'),
                    help='compiled branch enumerator (default: build/ppas_branch_census)')
args=parser.parse_args()
enumerator=args.enumerator.resolve()
if not enumerator.is_file():
    parser.error('branch enumerator missing; first run experiments/reproduce.py or pass --enumerator PATH')
reports=[]
for p in [7,11,13,31]:
    with tempfile.TemporaryDirectory() as td:
        path=Path(td)/'models.csv'
        run=subprocess.run([str(enumerator),str(p),str(path)],capture_output=True,text=True,check=True)
        old=json.loads((OUT/f'p{p}_maximal_graph.json').read_text())
        E=Census(p)
        for v in old['vertices']:
            if v['type']=='Jacobian':E.insert(v['coefficients'])
        n=len(E.vertices);hits=set();rows=0
        for row in csv.DictReader(path.open()):
            f=[int(row[f'c{i}']) for i in range(7)]
            i=E.insert(f);assert i<n
            hits.add(i);rows+=1
        assert len(hits)==n
        reports.append(dict(p=p,old_Jacobian_classes=n,exhaustive_classes=len(hits),model_rows=rows,enumerator_log=run.stdout+run.stderr))
        print(p,'independent exhaustive census matches',n,'Jacobian classes',flush=True)
(OUT/'branch_census_regressions.json').write_text(json.dumps(reports,indent=2))
