"""Complete p=7, ell=3 graph, reconstructed with a uniqueness proof.

Kani's 3-gluing criterion identifies every decomposed quotient from a split
source: it comes from a degree-2 elliptic map and is isomorphic to the source.
Remaining targets are uniquely determined here by rational 2-torsion rank.
This is not a general genus-two (3,3)-quotient implementation.
"""
from sage.all import *
from pathlib import Path
import json,itertools
OUT=Path(__file__).resolve().parent/'results'
d=json.loads((OUT/'p7_maximal_graph.json').read_text());p=7
F=GF(p);K=GF(p**4,'a');R=PolynomialRing(K,'x');x=R.gen()
es=[EllipticCurve(F,[e['a'],e['b']]) for e in d['elliptic_models']]
eks=[E.change_ring(K) for E in es];f3=GF(3)
bases=[];points=[];frobs=[]
for E in eks:
    roots=E.division_polynomial(3).roots(multiplicities=False);tors=[E(0)]
    for xx in roots:
        for yy in (x*x-(xx**3+E.a4()*xx+E.a6())).roots(multiplicities=False):tors.append(E(xx,yy))
    assert len(tors)==9 and all(3*P==E(0) for P in tors)
    P=next(P for P in tors if P);Q=next(Q for Q in tors if Q not in [E(0),P,2*P])
    coords={(a*P+b*Q):vector(f3,[a,b]) for a,b in itertools.product(range(3),repeat=2)}
    assert len(coords)==9
    def frob(P):return E(0) if not P else E(P[0]**p,P[1]**p)
    T=matrix(f3,[coords[frob(P)],coords[frob(Q)]]).transpose()
    bases.append((P,Q));points.append(coords);frobs.append(T)
pairings=[P.weil_pairing(Q,3) for P,Q in bases]
n=len(d['vertices']);B=matrix(ZZ,n);ranks=[]
for v in d['vertices']:
    if v['type']=='split':
        ranks.append(sum(2 if len(E.division_polynomial(2).roots())==3 else 1 for E in [es[i] for i in v['elliptic_ids']]))
    elif v['type']=='restriction':ranks.append(2)
    else:
        f=PolynomialRing(F,'x')(v['coefficients']);degrees=[g.degree() for g,e in f.factor()]
        if f.degree()==5:degrees.append(1)
        # These two examples both have rational Weierstrass points.
        assert 1 in degrees
        ranks.append(len(degrees)-2)
records=[]
for a,v in enumerate(d['vertices']):
    if v['type']!='split':continue
    i,j=v['elliptic_ids'];P,Q=bases[i];U,V=bases[j]
    anti=[]
    for M in GL(2,3):
        M=matrix(f3,M)
        if M*frobs[i]!=frobs[j]*M:continue
        if pairings[j]**int(M.det())!=pairings[i]**(-1):continue
        anti.append(M)
    assert len(anti)==4
    reducible=set()
    for phi in es[i].isogenies_prime_degree(2):
        for iso in phi.codomain().isomorphisms(es[j]):
            psi=iso*phi
            # Base-extend by evaluating the rational maps in the extension.
            xm,ym=psi.rational_maps()
            def apply(W):return eks[j](xm(W[0],W[1]),ym(W[0],W[1]))
            M=matrix(f3,[points[j][apply(P)],points[j][apply(Q)]]).transpose()
            assert any(M==N for N in anti)
            reducible.add(tuple(map(int,M.list())))
    B[a,a]=len(reducible)
    remaining=4-len(reducible)
    if remaining:
        targets=[b for b,w in enumerate(d['vertices']) if w['type']=='Jacobian' and ranks[b]==ranks[a]]
        assert len(targets)==1
        B[a,targets[0]]=remaining
    records.append({'source':a,'stable_anti_isometries':4,'reducible':len(reducible),'remaining_Jacobian_kernels':remaining})
# All four rational kernels on the restriction are geometric product kernels.
# Their elliptic quotients over F49 are maximal with the sole geometric j=6.
rid=next(i for i,v in enumerate(d['vertices']) if v['type']=='restriction')
E=es[0].change_ring(GF(49,'b'));phis=E.isogenies_prime_degree(3)
assert len(phis)==4 and all(phi.codomain().j_invariant()==E.j_invariant() for phi in phis)
B[rid,rid]=4
autos=[v['aut'] for v in d['vertices']]
for a,v in enumerate(d['vertices']):
    if v['type']!='Jacobian':continue
    for b,w in enumerate(d['vertices']):
        if w['type']=='Jacobian':continue
        k=QQ(B[b,a])*autos[a]/autos[b];assert k in ZZ;B[a,b]=k
    # The other Jacobian has a different rational 2-torsion rank.
    assert all(ranks[b]!=ranks[a] for b,w in enumerate(d['vertices']) if w['type']=='Jacobian' and b!=a)
    B[a,a]=4-sum(B.row(a));assert B[a,a]>=0
for i in range(n):
    assert sum(B.row(i))==4
    for j in range(n):assert B[i,j]*autos[j]==B[j,i]*autos[i]
G=Graph();G.add_vertices(range(n));G.add_edges((i,j) for i in range(n) for j in range(i+1,n) if B[i,j])
d.update({'ell':3,'adjacency':[list(map(int,B.row(i))) for i in range(n)],'outdegrees':[4]*n,
          'directed_kernels':24,'components':G.connected_components(sort=True),'rational_2_torsion_ranks':ranks,
          'construction':'Kani reducibility, explicit elliptic 2-maps, 2-torsion ranks, duality and degree',
          'three_gluing_checks':records})
d.pop('richelot_kernels',None)
(OUT/'p7_maximal_ell3_graph.json').write_text(json.dumps(d,indent=2,default=int))
print('COMPLETE p7 ell3',list(map(len,d['components'])),d['adjacency'],flush=True)
