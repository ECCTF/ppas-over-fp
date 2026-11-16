"""Complete rational (2,2) graphs in the maximal block, retaining F_p forms.

All Jacobian quotients and all decomposed-to-decomposed kernels are computed.
Edges from decomposed vertices to Jacobians are recovered by polarized duality.
Independent symplectic enumeration verifies every completed row degree.
The examples below have all supersingular elliptic j-invariants in F_p.
"""
from sage.all import *
from pathlib import Path
import csv,json,itertools
BASE=Path(__file__).resolve().parent;OUT=BASE/'results'
scope={'__name__':'helper','__file__':str(BASE/'verify_and_graphs.sage')}
exec((BASE/'verify_and_graphs.sage').read_text(),scope)
elliptic_classes=scope['elliptic_classes'];matchings=scope['matchings']

def lagrangians2():
    F=GF(2);J=block_diagonal_matrix(matrix(F,[[0,1],[1,0]]),matrix(F,[[0,1],[1,0]]));out=[]
    for a,b in itertools.combinations(range(4),2):
        slots=[(i,k) for i,c in enumerate([a,b]) for k in range(c+1,4) if k not in [a,b]]
        for vals in itertools.product(range(2),repeat=len(slots)):
            M=matrix(F,2,4);M[0,a]=1;M[1,b]=1
            for (i,k),v in zip(slots,vals):M[i,k]=v
            if M[0]*J*M[1]==0:out.append(M.row_space())
    assert len(out)==15
    return out
LAGS=lagrangians2()
def stable(T):return sum(all(T*v in L for v in L.basis()) for L in LAGS)

def run(p):
    F=GF(p);K=GF(p*p,'a');R=PolynomialRing(K,'x');x=R.gen()
    es=elliptic_classes(p);eks=[E.change_ring(K) for E in es]
    js=sorted({int(E.j_invariant()) for E in es});s=len(js)
    t=p//12+{1:0,5:1,7:1,11:2}[p%12];assert s==t
    assert s==(QuadraticField(-p).class_number() if p%8==7 else 2*QuadraticField(-p).class_number() if p%8==3 else QuadraticField(-p).class_number()/2)
    assert len(es)==2*s
    pairs=list(itertools.combinations_with_replacement(range(len(es)),2))
    vertices=[]
    for i,j in pairs:
        vertices.append({'label':f'S{len(vertices)}','type':'split','elliptic_ids':[i,j],
                         'aut':len(es[i].automorphisms())*len(es[j].automorphisms())*(2 if i==j else 1)})
    pairidx={pair:i for i,pair in enumerate(pairs)};nr=len(vertices);residx={j:nr+i for i,j in enumerate(js)}
    for j in js:
        E=next(E for E in eks if int(E.j_invariant())==j)
        vertices.append({'label':f'R{j}','type':'restriction','j':j,'aut':2*len(E.automorphisms())})
    nj=len(vertices)
    with open(OUT/f'p{p}_maximal_curves.csv') as f:rows=[{k:int(v) for k,v in row.items()} for row in csv.DictReader(f)]
    for row in rows:
        f=R([row[f'c{i}'] for i in range(7)])
        C=HyperellipticCurve(PolynomialRing(F,'x')(f))
        assert C.frobenius_polynomial().list()==[p*p,0,2*p,0,1]
        assert C.hasse_witt_matrix()==zero_matrix(F,2) if hasattr(C,'hasse_witt_matrix') else True
        vertices.append({'label':f'J{row["class_id"]}','type':'Jacobian','coefficients':[row[f'c{i}'] for i in range(7)],'aut':row['aut_order']})
    orbit={int(a):int(b) for a,b in csv.reader(open(OUT/f'p{p}_maximal_orbits.csv'))}
    nu=next(a for a in range(2,p) if not F(a).is_square())
    def curveid(g):
        assert all(c**p==c for c in g)
        d=g.degree();assert d in [5,6]
        lead=1 if F(g[d]).is_square() else nu;g=g*(K(lead)/g[d])
        code=sum(int(F(g[i]))*p**i for i in range(d))+(p**6 if d==5 else 0)+(p**6+p**5 if lead==nu else 0)
        return nj+orbit[code]
    n=len(vertices);B=matrix(ZZ,n);records=[];degrees=[None]*n
    def add(i,j,w=1):B[i,j]+=w
    # Split-to-split product kernels: both elliptic kernels must be rational.
    nbr=[]
    for E in es:
        targets=[]
        for phi in E.isogenies_prime_degree(2):
            hits=[i for i,D in enumerate(es) if phi.codomain().is_isomorphic(D)];assert len(hits)==1
            targets.append(hits[0])
        nbr.append(targets)
    I=identity_matrix(GF(2),2);Z=zero_matrix(GF(2),2);U=matrix(GF(2),[[1,1],[0,1]])
    ef=[]
    for E in es:
        roots=(x**3+K(E.a4())*x+K(E.a6())).roots(multiplicities=False)
        assert len(roots)==3
        ef.append(I if all(r**p==r for r in roots) else U)
    for a,(i,j) in enumerate(pairs):
        degrees[a]=stable(block_diagonal_matrix(ef[i],ef[j]))
        for u in nbr[i]:
            for v in nbr[j]:add(a,pairidx[tuple(sorted((u,v)))])
        # A non-product 2-gluing is decomposed iff induced by an elliptic
        # isomorphism. +/- isomorphisms induce the same graph on 2-torsion.
        used=set()
        for psi in eks[i].isomorphisms(eks[j]):
            u=psi.tuple()[0]
            key=tuple(sorted([str(u),str(-u)]))
            if key in used:continue
            used.add(key)
            if u**p==u:add(a,pairidx[(j,j)])
            elif u**p==-u:add(a,residx[int(es[j].j_invariant())])
    # Restrictions: their geometric product kernels are paired by Frobenius.
    # All E[2] is rational over F_{p^2}; each of the three lines contributes.
    for j in js:
        a=residx[j];E=next(E for E in eks if int(E.j_invariant())==j)
        assert E.trace_of_frobenius()==-2*p
        degrees[a]=stable(block_matrix([[Z,I],[I,Z]]))
        phis=E.isogenies_prime_degree(2);assert len(phis)==3
        for phi in phis:
            jj=phi.codomain().j_invariant();assert jj**p==jj
            add(a,residx[int(F(jj))])
    # Jacobian kernels: enumerate every Frobenius-stable Weierstrass matching.
    def elliptic(g,field):
        g=PolynomialRing(field,'z')(g);a,b,c,d=g[3],g[2],g[1],g[0]
        return EllipticCurve(field,[0,b,0,a*c,a*a*d])
    for r,row in enumerate(rows):
        a=nj+r;f=R([row[f'c{i}'] for i in range(7)])
        roots=f.roots(multiplicities=False)
        if f.degree()==5:roots.append(None)
        assert len(roots)==6
        rho=[roots.index(z**p) if z is not None else roots.index(None) for z in roots]
        degrees[a]=0
        for m in matchings(list(range(6))):
            if {frozenset((rho[i],rho[j])) for i,j in m}!={frozenset(q) for q in m}:continue
            degrees[a]+=1
            qs=[R(prod(x-roots[i] for i in pair if roots[i] is not None)) for pair in m];qs[0]*=f.leading_coefficient()
            A,D,G=qs;delta=matrix(K,[[q[2],q[1],q[0]] for q in qs]).det()
            rec={'source':a,'matching':m,'degenerate':not bool(delta)}
            if delta:
                br=lambda u,v:u.derivative()*v-u*v.derivative()
                g=-br(D,G)*br(G,A)*br(A,D)/delta;b=curveid(g)
                rec['quotient_coefficients']=[int(F(c)) for c in g]
            else:
                disc=(A[1]+x*D[1])**2-4*(A[2]+x*D[2])*(A[0]+x*D[0])
                ts=disc.roots(multiplicities=False)
                fixed=[-(A[1]+v*D[1])/(2*(A[2]+v*D[2])) if A[2]+v*D[2] else None for v in ts]
                if disc.degree()<2:fixed.append(-D[1]/(2*D[2]) if D[2] else None)
                assert len(fixed)==2 and fixed[0]!=fixed[1]
                ls=[x-z if z is not None else R(1) for z in fixed]
                mat=matrix(K,3,2,lambda i,j:(ls[j]**2)[i])
                coeff=[mat.solve_right(vector(K,[q[i] for i in range(3)])) for q in qs]
                g1=prod(u*x+v for u,v in coeff);g2=prod(u+v*x for u,v in coeff)
                rec['fixed_points']=[str(z) for z in fixed]
                if all(z is None or z**p==z for z in fixed):
                    ids=[]
                    for g in [g1,g2]:
                        assert all(c**p==c for c in g)
                        E=elliptic(g,F);hits=[i for i,D in enumerate(es) if E.is_isomorphic(D)];assert len(hits)==1
                        ids.append(hits[0])
                    b=pairidx[tuple(sorted(ids))]
                else:
                    assert fixed[0]**p==fixed[1]
                    E=elliptic(g1,K);assert E.trace_of_frobenius()==-2*p
                    b=residx[int(F(E.j_invariant()))]
            add(a,b);rec['target']=b;records.append(rec)
    # Recover reverse edges only when their forward count was computed above.
    autos=[v['aut'] for v in vertices]
    for i in range(nj):
        for j in range(nj,n):
            count=QQ(B[j,i])*autos[i]/autos[j];assert count in ZZ
            B[i,j]=count
    for i in range(nr,nj):
        for j in range(nr):
            count=QQ(B[j,i])*autos[i]/autos[j];assert count in ZZ
            B[i,j]=count
    for i in range(n):
        assert sum(B.row(i))==degrees[i],(p,vertices[i],list(B.row(i)),degrees[i])
        for j in range(n):assert B[i,j]*autos[j]==B[j,i]*autos[i],(p,i,j)
    # Independent old exhaustive census agrees for the three smallest primes.
    if p in [7,11,13]:
        old=list(csv.DictReader(open(OUT/f'p{p}_curves.csv')))
        selected=[v for v in old if v['weil_a2']==str(2*p)]
        assert len(selected)==len(rows)
        assert sorted(tuple(int(v[f'c{k}']) for k in range(7)) for v in selected)==sorted(tuple(v['coefficients']) for v in vertices[nj:])
    graph=Graph();graph.add_vertices(range(n));graph.add_edges((i,j) for i in range(n) for j in range(i+1,n) if B[i,j])
    data={'p':p,'ell':2,'weil_polynomial':[p*p,0,2*p,0,1],
          'elliptic_models':[{'a':int(E.a4()),'b':int(E.a6()),'j':int(E.j_invariant())} for E in es],
          'vertices':vertices,'adjacency':[list(map(int,B.row(i))) for i in range(n)],
          'outdegrees':degrees,'counts':{'split':nr,'restriction':nj-nr,'Jacobian':n-nj,'total':n},
          'directed_kernels':int(sum(B.list())),'components':sorted(graph.connected_components(sort=True),key=lambda c:(-len(c),c)),
          'richelot_kernels':records,'checks':['Sage Weil polynomials','complete reduced-model census','all 15 geometric kernels partitioned','every rational row degree','weighted duality']}
    (OUT/f'p{p}_maximal_graph.json').write_text(json.dumps(data,indent=2,default=int))
    with open(OUT/f'p{p}_maximal_ell2_adjacency.csv','w',newline='') as f:
        w=csv.writer(f);w.writerow(['source/target']+[v['label'] for v in vertices])
        for i,v in enumerate(vertices):w.writerow([v['label']]+list(B.row(i)))
    print('COMPLETE',p,data['counts'],'edges',data['directed_kernels'],'components',list(map(len,data['components'])),flush=True)

if __name__=='__main__':
    for p in [7,11,13,31]:run(p)
