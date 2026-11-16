"""Check the exchanged-factor descent equation for graph kernels.

Standard Python.  Enumerates all anti-isometries of two symplectic
planes and compares psi*v*psi=u with direct Frobenius stability of
their graph bases.  The restriction module has Frobenius (x,y)->(v*y,u*x)
with det(u)=det(v)=p and u*v=v*u=-p.  Every conjugate symplectic basis
choice is covered by this model.  This is a local module check, not an
independent global PPAS census or quotient-isogeny calculation.
"""
from itertools import product
from pathlib import Path
import json


def mul(a, b, q):
    return tuple(sum(a[2*i+k]*b[2*k+j] for k in range(2)) % q
                 for i in range(2) for j in range(2))


def apply(a, x, q):
    return ((a[0]*x[0]+a[1]*x[1]) % q,
            (a[2]*x[0]+a[3]*x[1]) % q)


def verify(p, q):
    assert p != q
    u = (p % q, 0, 0, 1)
    v = ((-1) % q, 0, 0, (-p) % q)
    scalar = ((-p) % q, 0, 0, (-p) % q)
    assert mul(u, v, q) == mul(v, u, q) == scalar
    good = []
    candidates = 0
    for psi in product(range(q), repeat=4):
        if (psi[0]*psi[3]-psi[1]*psi[2]) % q != (-1) % q:
            continue
        candidates += 1
        equation = mul(mul(psi, v, q), psi, q) == u
        # A basis vector (x, psi*x) remains in the graph precisely when
        # its transformed second component is psi of the first component.
        direct = all(apply(u, x, q) ==
                     apply(psi, apply(v, apply(psi, x, q), q), q)
                     for x in ((1, 0), (0, 1)))
        assert equation == direct
        if direct:
            good.append(list(psi))
    if q == 2:
        expected = 4
    else:
        symbol = pow((-p) % q, (q-1)//2, q)
        expected = 2 if symbol == 1 else 0
    assert len(good) == expected
    assert candidates == q*(q*q-1)
    row = dict(p=p, ell=q, anti_isometries=candidates,
               stable_graph_kernels=len(good), product_kernels=q+1,
               total_rational_kernels=q+1+len(good))
    if q == 2:
        # Here u=v=1.  Exchanged descent requires psi^2=1, not the
        # tautological separate-factor commutation psi*1=1*psi.
        row['naive_separate_factor_equivariance_count'] = candidates
        assert candidates == 6 and len(good) == 4
    return row


def main():
    rows = [verify(p, q) for p in (7, 11, 13, 31)
            for q in (2, 3, 5, 7) if p != q]
    result = dict(method='Exact anti-isometry and exchanged graph enumeration',
                  equation='psi*v*psi=u', cases=rows)
    path = Path(__file__).resolve().parent/'results/restriction_descent_verification.json'
    path.write_text(json.dumps(result, indent=2)+'\n', encoding='utf8')
    print(f'{len(rows)} local cases passed; at ell=2, exchanged descent selects 4 of 6 graph kernels.')
    print(f'Wrote {path}')


if __name__ == '__main__':
    main()
