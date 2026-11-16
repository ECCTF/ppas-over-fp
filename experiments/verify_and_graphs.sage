from sage.all import *
import csv, json, array, itertools, time
from pathlib import Path
BASE=Path(__file__).resolve().parent
OUT=BASE/'results'

def matchings(items):
    if not items:
        yield []
        return
    a=items[0]
    for b in items[1:]:
        for rest in matchings([i for i in items if i not in (a,b)]):
            yield [(a,b)]+rest

def read_curves(p):
    with open(OUT/('p%s_curves.csv'%p)) as stream:
        return [{k:int(v) if v else None for k,v in row.items()} for row in csv.DictReader(stream)]

def census_verify(p):
    rows=read_curves(p); F=GF(p); R=PolynomialRing(F,'x'); x=R.gen()
    counts={}; mass=QQ(0); blocks={}
    for row in rows:
        f=R([row['c%s'%i] for i in range(7)])
        C=HyperellipticCurve(f)
        H=f**((p-1)//2)
        M=matrix(F,[[H[p-1],H[p-2]],[H[2*p-1],H[2*p-2]]])
        prank=(M*M).rank(); anum=2-M.rank()
        assert prank==row['p_rank'] and anum==row['a_number']
        mass+=1/QQ(row['aut_order'])
        key='%s,%s'%(prank,anum);counts[key]=counts.get(key,0)+1
        if prank==0:
            wp=C.frobenius_polynomial()
            assert list(wp)==[p*p,p*row['weil_a1'],row['weil_a2'],row['weil_a1'],1], (p,row,wp)
            key='%s,%s'%(row['weil_a2']//p,anum);blocks[key]=blocks.get(key,0)+1
    assert mass==p**3, (p,mass)
    result={'p':p,'classes':len(rows),'counts_prank_anumber':counts,'ss_blocks_a2_over_p_anumber':blocks,'mass':str(mass),'verified_curves':len(rows)}
    print('CENSUS',result,flush=True)
    (OUT/('p%s_verified.json'%p)).write_text(json.dumps(result,indent=2,default=int))
    return result

def elliptic_classes(p):
    F=GF(p); out=[]
    nu=next(F(a) for a in range(2,p) if not F(a).is_square())
    for j in F:
        E=EllipticCurve(j=j)
        if not E.is_supersingular():continue
        if j==0:
            candidates=[EllipticCurve(F,[0,b]) for b in F if b]
        elif j==F(1728):
            candidates=[EllipticCurve(F,[a,0]) for a in F if a]
        else:candidates=[E,E.quadratic_twist(nu)]
        for E in candidates:
            if not any(E.is_isomorphic(D) for D in out):out.append(E)
    out.sort(key=lambda E:(int(E.j_invariant()),tuple(map(int,E.a_invariants()))))
    return out

def elliptic_product_graphs(p):
    es=elliptic_classes(p); data={'p':p,'elliptic_classes':[{'id':i,'j':int(E.j_invariant()),'a':int(E.a4()),'b':int(E.a6()),'aut':len(E.automorphisms())} for i,E in enumerate(es)]}
    s=len(set(int(E.j_invariant()) for E in es));t=p//12+{1:0,5:1,7:1,11:2}[p%12]
    hminus=QuadraticField(-p).class_number();h3=QuadraticField(-3*p).class_number();D=QuadraticField(p).discriminant()
    zeta=QQ(sum(kronecker(D,a)*a*a for a in range(1,D)))/(24*D)
    # zeta_Q(sqrt(D))(-1) = B_2,chi / 24; B_2,chi = sum chi(a)*a^2 / D.
    chi=kronecker(2,p)
    pp=(QQ(9-2*chi)/2*zeta+QQ(3)/8*hminus+QQ(3+chi)/6*h3) if p%4==1 else (zeta/2+QQ(11-3*chi)/8*hminus+h3/6)
    data['counts']={'s':s,'t':t,'elliptic_forms':len(es),'split':len(es)*(len(es)+1)//2,'restriction':t+s+2*int(p%3==2)+int(p%4==3),'real_block':int(pp),'real_jacobian':int(pp-(t+s)/2),'zeta_minus_one':str(zeta),'h_minus_p':int(hminus),'h_minus_3p':int(h3)}
    assert pp in ZZ
    pairs=list(itertools.combinations_with_replacement(range(len(es)),2)); pairindex={pair:i for i,pair in enumerate(pairs)}
    for ell in [2,3]:
        ee=[]; nbr=[[] for E in es]
        for i,E in enumerate(es):
            for phi in E.isogenies_prime_degree(ell):
                hits=[j for j,D in enumerate(es) if phi.codomain().is_isomorphic(D)]
                assert len(hits)==1
                j=hits[0];nbr[i].append(j);ee.append({'source':i,'target':j,'kernel':[int(c) for c in phi.kernel_polynomial()]})
        edges={}
        for pair in pairs:
            i,j=pair
            for u in nbr[i]:
                for v in nbr[j]:
                    dest=tuple(sorted((u,v)));key=(pairindex[pair],pairindex[dest]);edges[key]=edges.get(key,0)+1
        autos=[len(es[i].automorphisms())*len(es[j].automorphisms())*(2 if i==j else 1) for i,j in pairs]
        for (i,j),w in edges.items():assert w*autos[j]==edges.get((j,i),0)*autos[i]
        data['ell%s'%ell]={'elliptic_edges':ee,'product_vertices':[list(q) for q in pairs],'product_edges':[[i,j,w] for (i,j),w in sorted(edges.items())],'product_only':True}
        with open(OUT/('p%s_product_ell%s_edges.csv'%(p,ell)),'w',newline='') as stream:
            w=csv.writer(stream);w.writerow(['source','target','kernel_multiplicity']);w.writerows(data['ell%s'%ell]['product_edges'])
    (OUT/('p%s_products.json'%p)).write_text(json.dumps(data,indent=2,default=int))
    print('PRODUCTS',p,data['counts'],flush=True)
    return data

def richelot_graph(p):
    rows=[r for r in read_curves(p) if r['p_rank']==0]; byid={r['class_id']:r for r in rows}
    orbit=array.array('i');orbit.frombytes((OUT/('p%s_orbit.bin'%p)).read_bytes())
    F=GF(p);R=PolynomialRing(F,'x');x=R.gen();nu=next(a for a in range(2,p) if not F(a).is_square())
    def classid(g):
        gg=R([F(c) for c in g]);d=gg.degree();assert d in [5,6]
        lead=1 if gg[d].is_square() else nu;gg*=F(lead)/gg[d]
        ix=sum(int(gg[i])*p**i for i in range(d))+(p**6 if d==5 else 0)+(p**6+p**5 if lead==nu else 0)
        return orbit[ix]
    edges=[];frontier=[];degrees={};records=[]
    for row in rows:
        sid=row['class_id'];f=R([row['c%s'%i] for i in range(7)])
        extdeg=lcm([g.degree() for g,e in f.factor()]);K=GF(p**extdeg,'a');S=PolynomialRing(K,'z');z=S.gen()
        roots=list(S(f).roots(multiplicities=False));
        if f.degree()==5:roots.append(None)
        assert len(roots)==6
        rho=[roots.index(t**p) if t is not None else roots.index(None) for t in roots]
        allm=list(matchings(list(range(6))));d=0
        for mi,m in enumerate(allm):
            pairset={frozenset(q) for q in m}
            if {frozenset((rho[a],rho[b])) for a,b in m}!=pairset:continue
            d+=1;qs=[]
            for a,b in m:
                q=S(1)
                for r in (roots[a],roots[b]):
                    if r is not None:q*=z-r
                qs.append(q)
            qs[0]*=K(f.leading_coefficient());A,B,C=qs
            # For columns x^2,x,1 the rational correspondence requires
            # -[B,C][C,A][A,B]/delta, retaining the F_p twist.
            delta=matrix(K,[[q[2],q[1],q[0]] for q in qs]).det()
            record={'source':sid,'matching':mi,'extension_degree':int(extdeg),'pairs':[list(q) for q in m]}
            if not delta:
                frontier.append(record);continue
            bracket=lambda a,b:a.derivative()*b-a*b.derivative()
            g=-bracket(B,C)*bracket(C,A)*bracket(A,B)/delta
            assert all(c**p==c for c in g)
            dest=classid(g);assert dest in byid
            assert byid[dest]['a_number']==row['a_number']
            assert byid[dest]['weil_a1']==row['weil_a1']
            assert byid[dest]['weil_a2']==row['weil_a2']
            edges.append([sid,dest,mi]);record['target']=dest;record['quotient_coefficients']=[int(F(c)) for c in g];records.append(record)
        degrees[str(sid)]=d
    weights={}
    for i,j,m in edges:weights[(i,j)]=weights.get((i,j),0)+1
    for (i,j),w in weights.items():
        assert w*byid[j]['aut_order']==weights.get((j,i),0)*byid[i]['aut_order'],(i,j,w,weights.get((j,i),0))
    result={'p':p,'vertices':rows,'edges':[[i,j,w] for (i,j),w in sorted(weights.items())],'geometrically_decomposed_kernel_frontier':frontier,'degrees':degrees,'quotients':records}
    (OUT/('p%s_richelot.json'%p)).write_text(json.dumps(result,indent=2,default=int))
    print('RICHELOT',p,'vertices',len(rows),'Jacobian kernels',len(edges),'decomposed kernels',len(frontier),'degrees',degrees,flush=True)
    return result

if __name__=='__main__':
    import sys
    job=sys.argv[1]
    if job=='verify':
        for p in [7,11,13]:census_verify(p)
    elif job=='products':
        for p in [7,11,13,431]:elliptic_product_graphs(p)
    elif job=='richelot':richelot_graph(int(sys.argv[2]))
