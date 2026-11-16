"""Enumerate Lagrangians; independently verify degree formulas and complete 2-edge totals."""
from sage.all import *
from pathlib import Path
import itertools,json
OUT=Path(__file__).resolve().parent/'results'

def lagrangians(q,J):
    F=GF(q); answer=[]
    for a,b in itertools.combinations(range(4),2):
        slots=[(i,k) for i,c in enumerate([a,b]) for k in range(c+1,4) if k not in [a,b]]
        for vals in itertools.product(range(q),repeat=len(slots)):
            M=matrix(F,2,4); M[0,a]=1;M[1,b]=1
            for (i,k),v in zip(slots,vals): M[i,k]=v
            if M[0]*J*M[1]==0: answer.append(M.row_space())
    assert len(answer)==(q+1)*(q*q+1)
    return answer

def stable_count(Frob,lags):
    return sum(all(Frob*v in L for v in L.basis()) for L in lags)

checks=[]
for p in [7,11,13,431]:
    for ell in [3,5,7,11]:
        if p==ell:continue
        F=GF(ell);j=matrix(F,[[0,1],[-1,0]])
        for sign in [-1,1]:
            if sign==-1:
                A=matrix(F,[[0,-p],[1,0]]); T=block_diagonal_matrix(A,A);J=block_diagonal_matrix(j,j)
            else:
                I=identity_matrix(F,2);Z=zero_matrix(F,2)
                T=block_matrix([[Z,p*I],[I,Z]]);J=block_diagonal_matrix(j,p*j)
            assert T.transpose()*J*T==F(p)*J
            actual=stable_count(T,lagrangians(ell,J));split=F(sign*p).is_square()
            expected=(ell+3 if split else ell+1) if sign==-1 else ((ell+1)**2 if split else ell*ell+1)
            assert actual==expected
            checks.append({'p':p,'ell':ell,'sign':sign,'stable':actual})

F=GF(2);I=identity_matrix(F,2);Z=zero_matrix(F,2);j=matrix(F,[[0,1],[1,0]])
J=block_diagonal_matrix(j,j);lags=lagrangians(2,J)
def elliptic_frob2(E):
    R=PolynomialRing(E.base_field(),'x');x=R.gen()
    n=len((x**3+E.a4()*x+E.a6()).roots())
    if n==3:return I
    if n==1:return matrix(F,[[1,1],[0,1]])
    assert n==0
    return matrix(F,[[0,1],[1,1]])

totals=[]
for p in [7,11,13,431]:
    products=json.loads((OUT/f'p{p}_products.json').read_text())
    es=[EllipticCurve(GF(p),[e['a'],e['b']]) for e in products['elliptic_classes']]
    mats=[elliptic_frob2(E) for E in es]
    split_degrees=[stable_count(block_diagonal_matrix(mats[i],mats[k]),lags)
                   for i,k in products['ell2']['product_vertices']]
    for idx,d in enumerate(split_degrees):
        partial=sum(w for i,k,w in products['ell2']['product_edges'] if i==idx)
        assert partial<=d
    res=json.loads((OUT/f'p{p}_restrictions.json').read_text())
    R=PolynomialRing(GF(p),'z');K=GF(p*p,'z',modulus=R(res['field_modulus']))
    def elt(coeff):return sum(K(c)*K.gen()**i for i,c in enumerate(coeff))
    restrictions=[]
    for orb in res['galois_orbits']:
        row=res['elliptic_curves'][orb[0]];E=EllipticCurve(K,[elt(row['a4']),elt(row['a6'])])
        B=elliptic_frob2(E);T=block_matrix([[Z,B],[I,Z]])
        assert T.transpose()*J*T==J
        restrictions.append(stable_count(T,lags))
    entry={'p':p,'split_2_degrees':split_degrees,'restriction_2_degrees':restrictions,
           'split_2_edges':sum(split_degrees),'restriction_2_edges':sum(restrictions)}
    if p!=431:
        richelot=json.loads((OUT/f'p{p}_richelot.json').read_text())
        entry['jacobian_2_edges']=sum(richelot['degrees'].values())
        entry['all_2_edges']=sum(split_degrees)+sum(restrictions)+entry['jacobian_2_edges']
    totals.append(entry)
    print({k:v for k,v in entry.items() if not k.endswith('_degrees')},flush=True)
(OUT/'degree_checks.json').write_text(json.dumps({'quadratic_checks':checks,'two_edge_totals':totals},indent=2))
