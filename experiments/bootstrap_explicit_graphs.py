"""Independent elliptic models for new complete graph experiments (Sage)."""
from sage.all import *
from sage.schemes.elliptic_curves.ell_finite_field import supersingular_j_polynomial
from pathlib import Path
import sys, json, itertools

OUT = Path(__file__).resolve().parent / 'results'

def run(p):
    F = GF(p)
    forms = {}
    for a in range(p):
        for b in range(p):
            if (4*a**3+27*b**2)%p == 0:
                continue
            key = min((a*pow(u,4,p)%p,b*pow(u,6,p)%p) for u in range(1,p))
            if key != (a,b):
                continue
            E = EllipticCurve(F,[a,b])
            if E.cardinality() == p+1:
                forms[key] = E
    rows = [dict(id=i,a=a,b=b,j=int(E.j_invariant()))
            for i,((a,b),E) in enumerate(sorted(forms.items()))]
    ids = {(e['a'],e['b']):e['id'] for e in rows}
    targets = []
    for e in rows:
        ts=[]
        for r in range(p):
            if (r**3+e['a']*r+e['b'])%p == 0:
                a=(-4*e['a']-15*r*r)%p
                b=(-22*r**3-8*e['a']*r)%p
                key=min((a*pow(u,4,p)%p,b*pow(u,6,p)%p) for u in range(1,p))
                ts.append(ids[key])
        targets.append(ts)
    from collections import Counter
    pairs=list(itertools.combinations_with_replacement(range(len(rows)),2))
    pi={v:i for i,v in enumerate(pairs)}
    edges=Counter()
    for i,(a,b) in enumerate(pairs):
        for c in targets[a]:
            for d in targets[b]:edges[i,pi[tuple(sorted((c,d)))]]+=1
    data=dict(p=p,elliptic_classes=rows,ell2=dict(product_vertices=pairs,
              product_edges=[(i,j,n) for (i,j),n in sorted(edges.items())]))
    (OUT/f'p{p}_products.json').write_text(json.dumps(data,indent=2))
    K=GF(p*p,'z'); gen=K.multiplicative_generator()
    js=list(supersingular_j_polynomial(p).change_ring(K).roots(multiplicities=False))
    if p%3==2:js.append(K(0))
    if p%4==3:js.append(K(1728))
    js=sorted(set(js),key=lambda j:tuple(map(int,j.polynomial().list())))
    records=[]
    for j in js:
        if j==0:candidates=[EllipticCurve(K,[0,gen**i]) for i in range(6)]
        elif j==K(1728):candidates=[EllipticCurve(K,[gen**i,0]) for i in range(4)]
        else:
            E=EllipticCurve(j=j);candidates=[E,E.quadratic_twist(gen)]
        maximal=[E for E in candidates if E.cardinality()==(p+1)**2]
        assert len(maximal)==1
        E=maximal[0]
        records.append(dict(a4=list(map(int,E.a4().polynomial().list())),
            a6=list(map(int,E.a6().polynomial().list())),trace=-2*p))
    data=dict(p=p,field_modulus=list(map(int,K.modulus().list())),
        geometric_j=len(js),rational_j=sum(j**p==j for j in js),
        elliptic_curves=records,scope='unique maximal form for every geometric supersingular j')
    assert len(rows)==2*data['rational_j']
    (OUT/f'p{p}_restrictions.json').write_text(json.dumps(data,indent=2))
    print(p,'elliptic Fp forms',len(rows),'geometric j',len(js),flush=True)

if __name__=='__main__':
    for p in map(int,sys.argv[1:]):run(p)
