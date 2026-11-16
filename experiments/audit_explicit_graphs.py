"""Compare an exhaustive branch-model census to graph vertices, then export.

The C++ census never reads graph data. Coverage uses only rational PGL2
normalizations, the maximal branch splitting property, and point counts.
No class-group component theorem is used here.
"""
from pathlib import Path
from collections import Counter,deque
from math import comb
import csv,hashlib,itertools,json,sys,time
from rational_census import Census
from restriction_gluing_seeds import restriction_models

BASE=Path(__file__).resolve().parent
OUT=BASE/'results'

def binary_transform(f,M,p):
    # Separate integer-binomial implementation; no generator polynomial code.
    a,b,c,d=M; result=[0]*7
    for i,co in enumerate(f):
        for j in range(i+1):
            for k in range(7-i):
                result[j+k]=(result[j+k]+co*comb(i,j)*pow(a,j,p)*pow(b,i-j,p)
                    *comb(6-i,k)*pow(c,k,p)*pow(d,6-i-k,p))%p
    return result

def witness(E,f,target):
    F=E.F;p=E.p;g=E.vertices[target]['coefficients']
    cf=E.canonical(E.roots(f));cg=E.canonical_data[target]
    for tg in cg[1]:
        M=E.matmul(E.inverse(tg),cf[1][0])
        scale=F.inv(next(v for v in M if v));M=[F.mul(x,scale) for x in M]
        if any(x>=p for x in M):continue
        assert (M[0]*M[3]-M[1]*M[2])%p
        h=binary_transform(g,M,p);ff=list(f)+[0]*(7-len(f))
        idx=next(i for i,x in enumerate(ff) if x)
        ratio=h[idx]*pow(ff[idx],-1,p)%p
        if ratio and pow(ratio,(p-1)//2,p)==1 and all(x==ratio*y%p for x,y in zip(h,ff)):
            return dict(matrix=M,y_multiplier=next(x for x in range(1,p) if x*x%p==ratio))
    raise AssertionError(('no explicit isomorphism witness',p,f,target))

def dump_csv(path,rows,fields):
    with path.open('w',newline='',encoding='utf-8-sig') as out:
        w=csv.DictWriter(out,fieldnames=fields);w.writeheader();w.writerows(rows)

def equation(f):
    return 'y^2 = '+' + '.join(str(c)+(f'*x^{i}' if i>1 else '*x' if i else '') for i,c in enumerate(f) if c)

def audit(p):
    start=time.time();graphpath=OUT/f'p{p}_rational_graph.json'
    g=json.loads(graphpath.read_text());census=json.loads((OUT/f'p{p}_rational_census.json').read_text())
    E=Census(p)
    for v in census['vertices']:
        assert E.insert(v['coefficients'],v['roots'])==v['id']
    n=len(E.vertices);covered=Counter();cert=[]
    rawpath=OUT/f'p{p}_exhaustive_branch_models.csv'
    with rawpath.open() as stream:
        for row_number,row in enumerate(csv.DictReader(stream),2):
            f=[int(row[f'c{i}']) for i in range(7)];f=list(E.normalize(f))
            idx=E.insert(f)
            assert idx<n,('EXHAUSTIVE CENSUS FOUND MISSING VERTEX',p,idx,f)
            proof=witness(E,f,idx)
            covered[idx]+=1
            cert.append(dict(csv_line=row_number,target_label=f'J{idx}',normalized_coefficients=f,**proof))
    assert set(covered)==set(range(n)),('graph vertex missing in exhaustive census',p)
    # Independent 120-triple isomorphism search, without canonical minimizers.
    for ids in E.buckets.values():
        for i in ids:
            f=E.vertices[i]['coefficients']
            assert E.isomorphisms(f,f)==E.vertices[i]['aut']
        for i,j in itertools.combinations(ids,2):
            assert not E.isomorphisms(E.vertices[i]['coefficients'],E.vertices[j]['coefficients'],True)
    # Restriction forms independently enumerated from ALL supersingular j.
    rmodels=restriction_models(p,E)
    model_by_orbit={tuple(o):(a,b) for a,b,o in rmodels}
    assert set(model_by_orbit)=={tuple(v['j_orbit']) for v in g['vertices'] if v['type']=='restriction'}
    elliptics=g['elliptic_classes'];vertices=g['vertices'];N=len(vertices)
    B={(e['source'],e['target']):e['multiplicity'] for e in g['edges']}
    neigh=[[] for _ in vertices]
    for (i,j),m in B.items():
        assert m>0 and B.get((j,i),0)>0
        assert m*vertices[j]['aut']==B[j,i]*vertices[i]['aut']
        neigh[i].append(j)
    parent={0:None};depth={0:0};queue=deque([0])
    while queue:
        i=queue.popleft()
        for j in sorted(neigh[i]):
            if j not in parent:parent[j]=i;depth[j]=depth[i]+1;queue.append(j)
    assert len(parent)==N
    rows=[]
    for v in vertices:
        i=v['id'];assert sum(B[i,j] for j in neigh[i])==v['degree']
        if v['type']=='Jacobian':model=equation(v['coefficients'])
        elif v['type']=='split':
            factors=[elliptics[k] for k in v['elliptic_ids']]
            model=' x '.join(f"E{e['id']}: y^2=x^3+{e['a']}*x+{e['b']}" for e in factors)
        else:
            a,b=model_by_orbit[tuple(v['j_orbit'])]
            v['elliptic_model_fp2']={'a':a,'b':b}
            model=f'Res(Fp2/Fp)(y^2=x^3+({a%p}+{a//p}*u)*x+({b%p}+{b//p}*u)); u^2={E.F.nu}'
        v['model_description']=model
        rows.append(dict(id=i,label=v['label'],type=v['type'],model=model,aut_order=v['aut'],degree=v['degree'],
            component=0,parent_label='' if parent[i] is None else vertices[parent[i]]['label'],distance_from_S0=depth[i]))
    dump_csv(OUT/f'p{p}_all_vertices.csv',rows,list(rows[0]))
    edge_rows=[dict(source_id=i,source=vertices[i]['label'],target_id=j,target=vertices[j]['label'],multiplicity=m)
               for (i,j),m in sorted(B.items())]
    dump_csv(OUT/f'p{p}_all_edges.csv',edge_rows,list(edge_rows[0]))
    with (OUT/f'p{p}_adjacency.csv').open('w',newline='',encoding='utf-8-sig') as stream:
        w=csv.writer(stream);w.writerow(['source/target']+[v['label'] for v in vertices])
        for i,v in enumerate(vertices):w.writerow([v['label']]+[B.get((i,j),0) for j in range(N)])
    tree=[dict(source=vertices[parent[i]]['label'],target=vertices[i]['label'],multiplicity=B[parent[i],i],depth=depth[i])
          for i in range(1,N)]
    dump_csv(OUT/f'p{p}_spanning_tree.csv',tree,list(tree[0]))
    (OUT/f'p{p}_exhaustive_isomorphism_certificates.json').write_text(json.dumps(cert,indent=2))
    g['verification_scope']='All Jacobian vertices independently covered by exhaustive branch models; graph connectivity by BFS, without class-group criterion.'
    graphpath.write_text(json.dumps(g,indent=2))
    report=dict(p=p,counts=g['counts'],vertices=N,directed_kernels=sum(B.values()),
        adjacency_nonzero_entries=len(B),self_loop_kernels=sum(m for (i,j),m in B.items() if i==j),
        degree_histogram=dict(Counter(v['degree'] for v in vertices)),
        exhaustive_model_rows=len(cert),exhaustive_Jacobian_classes=n,
        every_exhaustive_model_has_explicit_Fp_isomorphism=True,every_graph_Jacobian_reached_by_exhaustive_census=True,
        all_automorphism_orders_rechecked_by_120_triples=True,
        restriction_forms_independently_enumerated=len(rmodels),
        bfs_components=1,spanning_tree_edges=len(tree),max_distance_from_S0=max(depth.values()),
        graph_sha256=hashlib.sha256(graphpath.read_bytes()).hexdigest(),
        exhaustive_csv_sha256=hashlib.sha256(rawpath.read_bytes()).hexdigest(),
        elapsed_seconds=round(time.time()-start,3))
    (OUT/f'p{p}_full_enumeration_audit.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':
    for p in map(int,sys.argv[1:] or [47,71,127]):audit(p)
