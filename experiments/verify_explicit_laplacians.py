"""Independent exact-rational Laplacian ranks (run with Sage Python)."""
from sage.all import matrix,QQ
from pathlib import Path
import json,sys
OUT=Path(__file__).resolve().parent/'results'
for p in map(int,sys.argv[1:] or [47,71,127]):
    g=json.loads((OUT/f'p{p}_rational_graph.json').read_text());n=len(g['vertices'])
    L=matrix(QQ,n,n)
    for e in g['edges']:
        i,j,m=e['source'],e['target'],e['multiplicity'];L[i,i]+=m;L[i,j]-=m
    assert all(sum(row)==0 for row in L.rows())
    rank=L.rank();assert rank==n-1
    report=dict(p=p,vertices=n,rank_Q=rank,nullity=n-rank)
    (OUT/f'p{p}_laplacian_check.json').write_text(json.dumps(report,indent=2))
    print(report,flush=True)
