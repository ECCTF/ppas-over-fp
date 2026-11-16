"""All three rational Richelot kernels of y^2=x^5-x at p=431."""
from sage.all import *
from pathlib import Path
import json
OUT=Path(__file__).resolve().parent/'results';p=431;F=GF(p);R=PolynomialRing(F,'x');x=R.gen()
f=x**5-x;C=HyperellipticCurve(f)
H=f**((p-1)//2);M=matrix(F,[[H[p-1],H[p-2]],[H[2*p-1],H[2*p-2]]])
assert M==0
wp=C.frobenius_polynomial();assert wp.list()==[p*p,0,2*p,0,1]
products=json.loads((OUT/'p431_products.json').read_text())
es=[EllipticCurve(F,[e['a'],e['b']]) for e in products['elliptic_classes']]
bracket=lambda a,b:a.derivative()*b-a*b.derivative()
result={'p':p,'source_coefficients':list(map(int,f.list())),
        'weil_polynomial':list(map(int,wp.list())),'a_number':2,'rational_2_degree':3,'kernels':[]}
for qs in [[x*x+1,x,x*x-1],[x*x+1,x*(x-1),x+1],[x*x+1,x*(x+1),x-1]]:
    assert prod(qs)==f
    delta=matrix(F,[[q[2],q[1],q[0]] for q in qs]).det()
    rec={'quadratics':[list(map(int,q.list())) for q in qs],'delta':int(delta)}
    if delta:
        A,B,D=qs;g=-bracket(B,D)*bracket(D,A)*bracket(A,B)/delta
        # Since x^6*f(1/x)=-f(x), reciprocal x coordinates preserve the
        # rational isomorphism class when the coefficient ratio is nonsquare.
        assert g/f in F
        ratio=F(g/f)
        M=[1,0,0,1] if ratio.is_square() else [0,1,1,0]
        multiplier=ratio if ratio.is_square() else -ratio
        assert multiplier.is_square()
        rec.update({'target':'same_jacobian','quotient_coefficients':list(map(int,g.list())),
                    'coefficient_ratio':int(ratio),'isomorphism_matrix':M,
                    'y_multiplier':int(multiplier.sqrt())})
    else:
        A,B,D=qs
        t=R.gen()
        disc=(A[1]+t*B[1])**2-4*(A[2]+t*B[2])*(A[0]+t*B[0])
        ts=disc.roots(multiplicities=False);assert len(ts)==2
        ss=[-(A[1]+v*B[1])/(2*(A[2]+v*B[2])) for v in ts]
        U=matrix(F,[[s*s for s in ss],[-2*s for s in ss],[1,1]])
        coeff=[U.solve_right(vector(F,q.list()+[0]*(3-len(q.list())))) for q in qs]
        g1=prod(a*x+b for a,b in coeff);g2=prod(a+b*x for a,b in coeff)
        assert all(g.degree()==3 and g.discriminant()!=0 for g in [g1,g2])
        # y^2=a*x^3+b*x^2+c*x+d -> Y^2=X^3+b*X^2+a*c*X+a^2*d.
        def elliptic(g):
            a,b,c,d=g[3],g[2],g[1],g[0]
            return EllipticCurve(F,[0,b,0,a*c,a*a*d])
        ids=[]
        for g in [g1,g2]:
            E=elliptic(g);hits=[i for i,D in enumerate(es) if E.is_isomorphic(D)]
            assert len(hits)==1;ids.append(hits[0])
        ids=sorted(ids);pairid=products['ell2']['product_vertices'].index(ids)
        rec.update({'target':'split_product','elliptic_class_ids':ids,'product_id':pairid,
                    'fixed_points':list(map(int,ss)),
                    'elliptic_cubics':[list(map(int,g.list())) for g in [g1,g2]]})
    result['kernels'].append(rec)
(OUT/'p431_local_richelot.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
