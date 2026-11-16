"""Independent enumeration of elliptic forms over F_{p^2} and Galois orbits."""
from sage.all import *
from sage.schemes.elliptic_curves.ell_finite_field import supersingular_j_polynomial
from pathlib import Path
import json

OUT=Path(__file__).resolve().parent/'results'
for p in [7,11,13,431]:
    K=GF(p*p,'z'); gen=K.multiplicative_generator()
    js=list(supersingular_j_polynomial(p).change_ring(K).roots(multiplicities=False))
    if p%3==2: js.append(K(0))
    if p%4==3: js.append(K(1728))
    js=sorted(set(js),key=lambda j:tuple(map(int,j.polynomial().list())))
    curves=[]
    for j in js:
        if j==0: curves.extend(EllipticCurve(K,[0,gen**i]) for i in range(6))
        elif j==K(1728): curves.extend(EllipticCurve(K,[gen**i,0]) for i in range(4))
        else:
            E=EllipticCurve(j=j); curves.extend([E,E.quadratic_twist(gen)])
    assert all(not E.is_isomorphic(D) for i,E in enumerate(curves) for D in curves[:i])
    records=[]; sigma=[]
    for i,E in enumerate(curves):
        conjugate=EllipticCurve(K,[a**p for a in E.a_invariants()])
        hit=[j for j,D in enumerate(curves) if conjugate.is_isomorphic(D)]
        assert len(hit)==1
        sigma.append(hit[0])
        records.append({'id':i,'a4':list(map(int,E.a4().polynomial().list())),
                        'a6':list(map(int,E.a6().polynomial().list())),
                        'trace':int(E.trace_of_frobenius()),'conjugate':hit[0]})
    assert all(sigma[sigma[i]]==i for i in range(len(curves)))
    orbits=[sorted(set([i,sigma[i]])) for i in range(len(curves)) if i<=sigma[i]]
    blocks={}
    for orb in orbits:
        b=-records[orb[0]]['trace']//p
        assert all(records[j]['trace']==-b*p for j in orb)
        blocks[str(b)]=blocks.get(str(b),0)+1
    result={'p':p,'field_modulus':list(map(int,K.modulus().list())),
            'geometric_j':len(js),'rational_j':sum(j**p==j for j in js),
            'elliptic_forms_over_fp2':len(curves),'fixed_forms':sum(i==j for i,j in enumerate(sigma)),
            'restriction_classes':len(orbits),'blocks_b_over_p':blocks,
            'elliptic_curves':records,'galois_orbits':orbits,'sage_version':version()}
    expected=json.loads((OUT/f'p{p}_products.json').read_text())['counts']
    assert result['geometric_j']==expected['t']
    assert result['rational_j']==expected['s']
    assert len(orbits)==expected['restriction']
    (OUT/f'p{p}_restrictions.json').write_text(json.dumps(result,indent=2))
    print(p,len(js),len(curves),len(orbits),blocks,flush=True)
