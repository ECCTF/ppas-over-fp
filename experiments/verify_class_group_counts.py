"""Exact arithmetic consistency checks for the class-group interpretation.

Standard Python 3; also runs under Sage Python. Quadratic ideal class numbers
are independently computed by reduced primitive binary quadratic forms.
Hermitian counts below come from the existing complete geometric computation;
this is NOT an independent implementation of Hermitian lattice enumeration.
"""
from fractions import Fraction
from math import gcd, isqrt
from pathlib import Path
from itertools import permutations
import json

RESULTS = Path(__file__).resolve().parent / 'results'


def class_number(D):
    assert D < 0 and D % 4 in (0, 1)
    forms = []
    for a in range(1, isqrt((-D)//3)+1):
        for b in range(-a, a+1):
            if (b*b-D) % (4*a):
                continue
            c = (b*b-D)//(4*a)
            if a > c or gcd(gcd(a,b),c) != 1:
                continue
            if (abs(b) == a or a == c) and b < 0:
                continue
            forms.append((a,b,c))
    return len(forms)


old = json.loads((RESULTS/'count_formula_verification.json').read_text())
rows = []
for record in old:
    p = record['p']
    hK = class_number(-p if p % 4 == 3 else -4*p)
    hR = class_number(-4*p)
    s, t = record['s'], record['t']
    assert hK == record['h_minus_p']
    assert 2*s == (hR+hK if p % 4 == 3 else hR)
    if p % 4 == 3:
        assert hK % 2 == 1
        assert hR == (1 if p % 8 == 7 else 3)*hK
    base = Fraction(2*record['geometric_Fp_fixed']-record['minimal_Fp_forms']-s*s)
    base += Fraction(s-t,2)
    assert base.denominator == 1
    row = dict(p=p, hK=hK, hR=hR, split=s*(2*s+1), restriction=(t+s)//2,
               Jacobian_trace_term_before_exceptional_correction=int(base))
    graphfile = RESULTS/f'p{p}_maximal_graph.json'
    if graphfile.exists():
        vertices = json.loads(graphfile.read_text())['vertices']
        total = len(vertices)
        jac = sum(v['type'] == 'Jacobian' for v in vertices)
        delta = jac-int(base)
        assert total == 2*record['geometric_Fp_fixed']-record['minimal_Fp_forms']+s*s+2*s+delta
        row.update(total=total, Jacobian=jac, exceptional_correction=delta)
        if p % 4 == 3:
            assert total % hK == 0
        if p % 8 == 7:
            row['quadratic_class_group_orbits_from_freeness'] = total//hK
    if p == 431:
        assert (hK,hR,s,t,row['split'],row['restriction'],int(base)) == (21,21,21,37,903,29,2556)
        row.update(Jacobian_residue_mod_hK=(-row['split']-row['restriction']) % hK,
                   exceptional_correction_residue_mod_hK=(-row['split']-row['restriction']-int(base)) % hK,
                   exact_Jacobian_not_evaluated=True)
    rows.append(row)
assert [r['exceptional_correction'] for r in rows[:4]] == [1,5,0,16]

def verify_p7_descent():
    """Compute constrained twists in the actual 48-element automorphism group.

    F49 = F7[u]/(u^2+1), encoded a+7b. No Sage dependency. The curve is
    the recorded maximal J1: y^2=x^6+x^4+4x^2+1.
    """
    p = 7
    def add(a,b): return (a%p+b%p)%p+p*((a//p+b//p)%p)
    def neg(a): return (-a%p)%p+p*((-(a//p))%p)
    def sub(a,b): return add(a,neg(b))
    def mul(a,b):
        x,y,z,w = a%p,a//p,b%p,b//p
        return (x*z-y*w)%p+p*((x*w+y*z)%p)
    def power(a,n):
        result=1
        while n:
            if n%2: result=mul(result,a)
            a=mul(a,a); n//=2
        return result
    def inv(a):
        assert a
        return power(a,47)
    def mm(a,b):
        return (add(mul(a[0],b[0]),mul(a[1],b[2])),
                add(mul(a[0],b[1]),mul(a[1],b[3])),
                add(mul(a[2],b[0]),mul(a[3],b[2])),
                add(mul(a[2],b[1]),mul(a[3],b[3])))
    def mi(m):
        scale=inv(sub(mul(m[0],m[3]),mul(m[1],m[2])))
        return tuple(mul(v,scale) for v in (m[3],neg(m[1]),neg(m[2]),m[0]))
    def triple(a,b,c):
        return (sub(c,a),neg(mul(b,sub(c,a))),sub(c,b),neg(mul(a,sub(c,b))))
    coeff=[1,0,4,0,1,0,1]
    def binary(a,b):
        value=0
        for i,co in enumerate(coeff):
            value=add(value,mul(co,mul(power(a,i),power(b,6-i))))
        return value
    roots=[a for a in range(49) if binary(a,1)==0]
    assert len(roots)==6
    source=triple(*roots[:3])
    group=set()
    for target in permutations(roots,3):
        m=mm(mi(triple(*target)),source)
        scale=inv(next(v for v in m if v))
        m=tuple(mul(v,scale) for v in m)
        if not all(binary(add(mul(m[0],r),m[1]),add(mul(m[2],r),m[3]))==0 for r in roots):
            continue
        ratio=binary(m[0],m[2])
        for e in range(1,49):
            if mul(e,e)==ratio: group.add(m+(e,))
    assert len(group)==48
    def gm(a,b):
        m=mm(a[:4],b[:4]); e=mul(a[4],b[4])
        scale=inv(next(v for v in m if v))
        return tuple(mul(v,scale) for v in m)+(mul(e,power(scale,3)),)
    def sigma(a): return tuple(power(v,p) for v in a)
    identity=(1,0,0,1,1)
    minus=(1,0,0,1,6)
    inverse={a:next(b for b in group if gm(a,b)==identity) for a in group}
    assert sum(sigma(a)==a for a in group)==8
    counts={}
    for name,sign in [('maximal',identity),('minimal',minus)]:
        allowed={a for a in group if gm(a,sigma(a))==sign}
        remaining=set(allowed); orbits=[]
        while remaining:
            a=next(iter(remaining))
            orbit={gm(gm(b,a),sigma(inverse[b])) for b in group}
            assert orbit <= allowed
            remaining-=orbit; orbits.append(len(orbit))
        burnside=sum(gm(b,a)==gm(a,sigma(b)) for a in allowed for b in group)
        assert burnside==len(group)*len(orbits)
        counts[name]=dict(forms=len(orbits), twist_elements=len(allowed), orbit_sizes=sorted(orbits))
    assert counts['maximal']['forms']==2 and counts['minimal']['forms']==1
    assert counts['maximal']['forms']+counts['minimal']['forms']-2==1
    return dict(p=7, geometric_group_order=48, rational_group_order=8,
                exceptional_Jacobian_correction=1, **counts)


descent=verify_p7_descent()
(RESULTS/'class_group_verification.json').write_text(
    json.dumps(dict(arithmetic=rows, independent_p7_descent=descent), indent=2)+'\n')
print(json.dumps(rows, indent=2))
print(json.dumps(descent, indent=2))
print('All quadratic class-number, orbit-divisibility and exact-correction checks passed; '
      'p=7 correction independently computed from constrained twists.')
