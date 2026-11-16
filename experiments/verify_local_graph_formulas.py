"""Independent finite-linear-algebra audit of the rational kernel formulas.

Enumerates RREF planes directly, without importing graph generators or
quotient data. These checks test local counts, not complete adjacency.
"""
from itertools import combinations, product
from fractions import Fraction
from pathlib import Path
import json


def planes(q):
    for piv in combinations(range(4), 2):
        slots = [(i, j) for i in range(2) for j in range(piv[i]+1, 4)
                 if j not in piv]
        for entries in product(range(q), repeat=len(slots)):
            rows = [[0]*4 for _ in range(2)]
            for i, j in enumerate(piv):
                rows[i][j] = 1
            for (i, j), value in zip(slots, entries):
                rows[i][j] = value
            yield piv, rows


def image(matrix, v, q):
    return [sum(a*b for a, b in zip(row, v)) % q for row in matrix]


def stable(matrix, piv, rows, q):
    for row in rows:
        w = image(matrix, row, q)
        if any((w[j]-sum(w[piv[i]]*rows[i][j] for i in range(2))) % q
               for j in range(4)):
            return False
    return True


def pairing(scale):
    return [[0, 1, 0, 0], [-1, 0, 0, 0],
            [0, 0, 0, scale], [0, 0, -scale, 0]]


def count(q, matrix, form):
    total = rational = products = 0
    for piv, rows in planes(q):
        if sum(a*b for a, b in zip(rows[0], image(form, rows[1], q))) % q:
            continue
        total += 1
        if not stable(matrix, piv, rows, q):
            continue
        rational += 1
        left = (rows[0][0]*rows[1][1]-rows[0][1]*rows[1][0]) % q
        right = (rows[0][2]*rows[1][3]-rows[0][3]*rows[1][2]) % q
        products += left == 0 and right == 0
    assert total == (q+1)*(q*q+1)
    return {'geometric_kernels': total, 'rational_kernels': rational,
            'product_kernels': products, 'graph_kernels': rational-products}


def block(a, b):
    return [row+[0, 0] for row in a]+[[0, 0]+row for row in b]


def main():
    rows = []
    for ell in (3, 5, 7, 11):
        for residue in range(1, ell):
            symbol = 1 if pow((-residue) % ell, (ell-1)//2, ell) == 1 else -1
            a = [[0, -residue], [1, 0]]
            split = count(ell, block(a, a), pairing(1))
            assert split['rational_kernels'] == ell+2+symbol
            assert split['product_kernels'] == (1+symbol)**2
            exchange = [[0, 0, -residue, 0], [0, 0, 0, -residue],
                        [1, 0, 0, 0], [0, 1, 0, 0]]
            restriction = count(ell, exchange, pairing(residue))
            assert restriction['rational_kernels'] == ell+2+symbol
            assert restriction['product_kernels'] == ell+1
            rows.append({'ell': ell, 'p_mod_ell': residue, 'symbol': symbol,
                         'split': split, 'restriction': restriction})
    ident, unipotent = [[1, 0], [0, 1]], [[1, 1], [0, 1]]
    dyadic = {}
    for name, a, b, degree, products in [('I,I', ident, ident, 15, 9),
                                        ('I,U', ident, unipotent, 3, 3),
                                        ('U,U', unipotent, unipotent, 3, 1)]:
        result = count(2, block(a, b), pairing(1))
        assert (result['rational_kernels'], result['product_kernels']) == (degree, products)
        dyadic[name] = result
    exchange = [[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]]
    dyadic['restriction'] = count(2, exchange, pairing(1))
    assert dyadic['restriction']['rational_kernels'] == 7
    assert dyadic['restriction']['product_kernels'] == 3
    # Independent arithmetic and group checks of two printed claims in XOZ v3.
    printed_count = Fraction(7**3+24*7**2+141*7-346, 2880)
    assert printed_count == Fraction(3,4)
    def multiply(a,b):
        return tuple(sum(a[2*i+k]*b[2*k+j] for k in range(2))
                     for i in range(2) for j in range(2))
    sign, swap = (-1,0,0,1), (0,1,1,0)
    assert multiply(sign,swap) != multiply(swap,sign)
    report = {'odd_residue_cases': rows, 'dyadic_types': dyadic,
              'XOZ_v3_theorem_2_4_at_p7': str(printed_count),
              'signed_permutation_group_is_nonabelian': True,
              'scope': 'Exact local kernel counts. No quotient or complete-graph census is inferred.'}
    out = Path(__file__).resolve().parent/'results/local_graph_formula_verification.json'
    out.write_text(json.dumps(report, indent=2), encoding='utf8')
    print('Passed', len(rows), 'odd residue cases and four dyadic cases.')
    print('Verified the geometric count and signed-permutation counterchecks.')


if __name__ == '__main__':
    main()
