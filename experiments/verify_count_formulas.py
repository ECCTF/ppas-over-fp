"""Check field-dependent class numbers and new congruence claims in Sage 10.8.

The asymptotic is proved in the manuscript; finite samples are not a proof of
an asymptotic. This script checks the exact auxiliary formulas, existing
censuses, leading-term evaluation, and the p = 7 quadratic base change.
"""
from sage.all import GF, QQ, QuadraticField, PolynomialRing, matrix, kronecker, HyperellipticCurve
from pathlib import Path
from itertools import permutations
import csv
import json

RESULTS = Path(__file__).resolve().parent / 'results'


def counts(p):
    hp = QuadraticField(-p).class_number()
    h2 = QuadraticField(-2*p).class_number()
    h3 = QuadraticField(-3*p).class_number()
    D = QuadraticField(p).discriminant()
    bernoulli = QQ(sum(kronecker(D, a)*a*a for a in range(1, D))) / D
    c = kronecker(2, p)
    s = hp/2 if p % 4 == 1 else hp if p % 8 == 7 else 2*hp
    t = p//12 + {1: 0, 5: 1, 7: 1, 11: 2}[p % 12]
    geom_j = (QQ((p-1)*(p*p+25*p+166))/2880
              - QQ(1-kronecker(-1, p))/32 + QQ(1-kronecker(-2, p))/8
              + QQ(1-kronecker(-3, p))/18 + QQ(4)/5*(p % 5 == 4))
    if p % 4 == 3:
        tau = (bernoulli/96 + QQ(h2)/8 + QQ(h3)/12
               + (QQ((p-1)*(9-4*c))/48 + QQ(p-c)/16
                  + QQ((1-kronecker(p, 3))*(3-c))/12)*hp)
        nplus = bernoulli/48 + QQ((11-3*c)*hp)/8 + QQ(h3)/6
    else:
        tau = ((9-2*c)*bernoulli/96 + QQ((4*p-1)*hp)/48 + QQ(h2)/8
               + QQ((3+kronecker(-2, p))*h3)/12
               + QQ((1-kronecker(p, 3))*hp)/12)
        nplus = ((9-2*c)*bernoulli/48 + QQ(3*hp)/8 + QQ((3+c)*h3)/6)
    quantities = [s, geom_j, tau, nplus, (t+s)/2]
    assert all(QQ(x).denominator() == 1 for x in quantities)
    return dict(p=p, s=int(s), t=t, h_minus_p=int(hp), B2=str(bernoulli),
                split=int(s*(2*s+1)), restriction=int((t+s)/2),
                maximal_Fp2_products=t*(t+1)//2, maximal_Fp2_Jacobians=int(geom_j),
                geometric_Fp_fixed=int(tau), minimal_Fp_forms=int(nplus),
                Jacobian_leading_term=str(QQ(p*s)/3))


published = {7:1, 11:2, 13:3, 17:5, 19:7, 23:10, 29:18, 31:20,
             37:31, 41:40, 43:45, 47:57, 53:78, 59:104, 61:113}
for p, expected in published.items():
    assert counts(p)['maximal_Fp2_Jacobians'] == expected
print('Geometric formula matches all 15 entries of IKO1986 Remark 3.4.', flush=True)

records = []
for p in [7, 11, 13, 31, 431]:
    row = counts(p)
    file = RESULTS/f'p{p}_maximal_curves.csv'
    if file.exists():
        curves = list(csv.DictReader(file.open()))
        row['maximal_Fp_Jacobians'] = len(curves)
        assert all(int(v['weil_a1']) == 0 and int(v['weil_a2']) == 2*p for v in curves)
        # Geometric invariants identify exactly the image under base change.
        field = GF(p)
        polynomials = PolynomialRing(field, 'x')
        scalings = PolynomialRing(field, 'u')
        u = scalings.gen()
        geometric_representatives = []
        def same_geometric_invariants(left, right):
            common = scalings(0)
            for a, b, weight in zip(left, right, [1, 2, 3, 5]):
                common = common.gcd(a*u**weight-b)
            return common.degree() > 0 and common != u**common.degree()
        for curve in curves:
            f_curve = polynomials([int(curve[f'c{i}']) for i in range(7)])
            inv = HyperellipticCurve(f_curve).igusa_clebsch_invariants()
            if not any(same_geometric_invariants(inv, prev) for prev in geometric_representatives):
                geometric_representatives.append(inv)
        row['geometric_Jacobians_with_maximal_Fp_descent'] = len(geometric_representatives)
        graph = json.loads((RESULTS/f'p{p}_maximal_graph.json').read_text())
        vertices, adjacency = graph['vertices'], graph['adjacency']
        for typ, expected in [('split', row['split']), ('restriction', row['restriction']),
                              ('Jacobian', len(curves))]:
            assert sum(v['type'] == typ for v in vertices) == expected
        for i, v in enumerate(vertices):
            degree = sum(adjacency[i])
            assert degree in ([3, 7] if p % 4 == 1 else [3, 7, 15])
            for j, w in enumerate(vertices):
                assert adjacency[i][j]*w['aut'] == adjacency[j][i]*v['aut']
        row['kernel_count_ell2'] = sum(map(sum, adjacency))
        row['trace_comparison_error'] = (len(vertices) + row['minimal_Fp_forms']
                                         - 2*row['geometric_Fp_fixed'])
    if p in [7, 11, 13]:
        # Independent full enumeration, including all supersingular Weil blocks.
        all_curves = list(csv.DictReader((RESULTS/f'p{p}_curves.csv').open()))
        jminus = sum(v['weil_a1'] == '0' and v['weil_a2'] == str(2*p) for v in all_curves)
        jplus = sum(v['weil_a1'] == '0' and v['weil_a2'] == str(-2*p) for v in all_curves)
        assert jminus == row['maximal_Fp_Jacobians']
        # The plus block has no split products; restrictions use minimal elliptic
        # Fp2 models, again one per conjugation orbit.
        assert jplus + row['restriction'] == row['minimal_Fp_forms']
        row['independent_minimal_Jacobians'] = jplus
    print(row, flush=True)
    records.append(row)

# All four involution types: direct action on perfect matchings of six letters.
def matchings(items):
    if not items:
        yield frozenset()
        return
    a = items[0]
    for b in items[1:]:
        for rest in matchings([x for x in items if x not in (a, b)]):
            yield rest | {frozenset((a, b))}

all_matchings = list(matchings(list(range(6))))
for sigma, expected in [([0,1,2,3,4,5],15), ([1,0,2,3,4,5],3),
                         ([1,0,3,2,4,5],3), ([1,0,3,2,5,4],7)]:
    actual = sum(frozenset(frozenset(sigma[x] for x in pair) for pair in m) == m
                 for m in all_matchings)
    assert actual == expected

# Explicitly check that the two distinct F7 Jacobians coalesce over F49.
F = GF(49, 'u')
R = PolynomialRing(F, 'x')
x = R.gen()
f = x**6+x**4+x**3+4*x**2+5*x
g = x**6+x**4+4*x**2+1
def triple_matrix(a,b,c):
    return matrix(F, [[c-a,-b*(c-a)], [c-b,-a*(c-b)]])
roots_f = f.roots(multiplicities=False)
roots_g = g.roots(multiplicities=False)
source = triple_matrix(*roots_f[:3])
maps = 0
for triple in permutations(roots_g, 3):
    a,b,c,d = (triple_matrix(*triple).inverse()*source).list()
    transformed = sum(g[i]*(a*x+b)**i*(c*x+d)**(6-i) for i in range(7))
    scalar = transformed.leading_coefficient()/f.leading_coefficient()
    if transformed == scalar*f and scalar.is_square():
        maps += 2
assert maps == 48
print('p=7: exactly 48 curve isomorphisms between C0 and C1 over F49.', flush=True)
print('All class-number, congruence, census and weighted-duality checks passed.', flush=True)
(RESULTS/'count_formula_verification.json').write_text(json.dumps(records, indent=2)+'\n')
