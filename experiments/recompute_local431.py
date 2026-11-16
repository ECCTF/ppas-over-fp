"""Verify the stored local p=431 example with rational twist normalization."""
from pathlib import Path
import json
from rational_census import Census

path=Path(__file__).resolve().parent/'results'/'p431_local_richelot.json'
data=json.loads(path.read_text());E=Census(431);F=E.F;P=E.P
f=data['source_coefficients']
for rec in data['kernels']:
    qs=rec['quadratics']
    assert P.mul(P.mul(qs[0],qs[1]),qs[2])==f
    delta=E.det3([[q[j] if j<len(q) else 0 for j in (2,1,0)] for q in qs])
    if not delta:continue
    def br(a,b):return P.sub(P.mul(E.derivative(a),b),P.mul(a,E.derivative(b)))
    g=E.scale(P.mul(P.mul(br(qs[1],qs[2]),br(qs[2],qs[0])),br(qs[0],qs[1])),F.neg(F.inv(delta)))
    ratio=F.div(g[-1],f[-1]);assert g==E.scale(f,ratio)
    M=[1,0,0,1] if E.chi[ratio]==1 else [0,1,1,0]
    sq=ratio if E.chi[ratio]==1 else (-ratio)%431
    e=pow(sq,(431+1)//4,431);assert e*e%431==sq
    assert E.transform(g,M)==E.scale(f,sq)
    rec.pop('square_ratio',None);rec.pop('square_root',None)
    rec.update(quotient_coefficients=g,coefficient_ratio=ratio,isomorphism_matrix=M,y_multiplier=e)
data['richelot_codomain_normalization']='-[B,C][C,A][A,B]/det(rows(q2,q1,q0))'
path.write_text(json.dumps(data,indent=2)+'\n')
print('Local p=431 smooth quotient remains Fp-isomorphic to the source via reciprocal x coordinates.')
