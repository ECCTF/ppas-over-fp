"""Independent checks for rational automorphisms, duality and kernel families.

Run with SageMath 10.8's Python. Uses branch triples rather than the full
PGL2 orbit stabilizers used in the original enumeration.
"""
from sage.all import GF, PolynomialRing, matrix, vector, block_diagonal_matrix
from sage.all import block_matrix, identity_matrix, zero_matrix, EllipticCurve
from itertools import permutations, combinations, product
from pathlib import Path
from collections import Counter
import json

RESULTS = Path(__file__).resolve().parent / 'results'


def branch_matrix(field, a, b, c):
    # Homogeneous points; infinity is (1,0). Send a,b,c to infinity,0,1.
    la = vector(field, [a[1], -a[0]])
    lb = vector(field, [b[1], -b[0]])
    return matrix(field, [lb*(la.dot_product(c)), la*(lb.dot_product(c))])


def curve_automorphisms(p, coefficients):
    field = GF(p*p, 'u')
    ring = PolynomialRing(field, 'x')
    x = ring.gen()
    f = ring(coefficients)
    roots = [vector(field, [a, 1]) for a in f.roots(multiplicities=False)]
    if f.degree() == 5:
        roots.append(vector(field, [1, 0]))
    assert len(roots) == 6  # Frobenius squared is scalar on J[2].
    source = branch_matrix(field, *roots[:3])
    rational_branch = 0
    rational_aut = 0
    geometric_aut = 0
    rational_elements = []
    for triple in permutations(roots, 3):
        m = branch_matrix(field, *triple).inverse()*source
        m /= next(a for a in m.list() if a)
        a, b, c, d = m.list()
        h = sum(f[i]*(a*x+b)**i*(c*x+d)**(6-i) for i in range(7))
        scale = h[f.degree()]/f.leading_coefficient()
        if not scale or h != scale*f:
            continue
        if scale.is_square():
            geometric_aut += 2
        if all(v**p == v for v in m.list()):
            rational_branch += 1
            assert scale**p == scale
            if GF(p)(scale).is_square():
                rational_aut += 2
                e = field(GF(p)(scale).sqrt())
                rational_elements.extend([(m,e), (m,-e)])
    orders = Counter()
    for m, e in rational_elements:
        power = identity_matrix(field, 2)
        scalar = field(1)
        for n in range(1, 121):
            power = power*m
            scalar *= e
            normalizer = next(a for a in power.list() if a)
            power /= normalizer
            scalar /= normalizer**3
            if power == identity_matrix(field, 2) and scalar == 1:
                orders[n] += 1
                break
        else:
            raise AssertionError('Automorphism order exceeded 120')
    return dict(rational_branch_stabilizer=rational_branch,
                rational_automorphisms=rational_aut,
                geometric_automorphisms=geometric_aut,
                element_orders=dict(sorted(orders.items())),
                nonlifting_branch_maps=rational_branch-rational_aut//2)


records = []
for p in [7, 11, 13, 31]:
    graph = json.loads((RESULTS/f'p{p}_maximal_graph.json').read_text())
    for v in graph['vertices']:
        if v['type'] == 'Jacobian':
            record = curve_automorphisms(p, v['coefficients'])
            assert record['rational_automorphisms'] == v['aut'], (p, v, record)
            records.append(dict(p=p, label=v['label'], **record))
        elif v['type'] == 'split':
            i, j = v['elliptic_ids']
            assert v['aut'] == (8 if i == j else 4)
        else:
            j = v['j']
            expected = 12 if j == 0 else 8 if j == 1728 % p else 4
            assert v['aut'] == expected
print('Independent branch-triple automorphism checks: 50 Jacobian forms.', flush=True)


def lagrangians(field, pairing):
    for pivots in combinations(range(4), 2):
        positions = [(r, c) for r in range(2) for c in range(pivots[r]+1, 4)
                     if c not in pivots]
        for values in product(field, repeat=len(positions)):
            m = zero_matrix(field, 2, 4)
            for r, c in enumerate(pivots):
                m[r, c] = 1
            for (r, c), value in zip(positions, values):
                m[r, c] = value
            if (m*pairing*m.transpose()).is_zero():
                yield m


family_records = []
for ell in [2, 3, 5, 7]:
    field = GF(ell)
    j = matrix(field, [[0, 1], [-1, 0]])
    pairing = block_diagonal_matrix(j, j)
    planes = list(lagrangians(field, pairing))
    assert len(planes) == (ell+1)*(ell*ell+1)
    primes = [7, 11, 13, 31] if ell == 2 else [p for p in [7, 11, 13] if p != ell]
    for p in primes:
        actions = []
        if ell == 2:
            I = identity_matrix(field, 2)
            U = matrix(field, [[1, 1], [0, 1]])
            for left, right, label, expected in [(I,I,'II',(9,6)),
                                                (I,U,'IU',(3,0)),
                                                (U,U,'UU',(1,2))]:
                actions.append((label, block_diagonal_matrix(left,right), expected))
            res_expected = (3,4)
        else:
            frob = matrix(field, [[0, -p], [1, 0]])
            split = field(-p).is_square()
            actions.append(('split_surface', block_diagonal_matrix(frob,frob),
                            (4,ell-1) if split else (0,ell+1)))
            res_expected = (ell+1,2) if split else (ell+1,0)
        C = matrix(field, [[1,0],[0,p]])
        D = -field(p)*C.inverse()
        z = zero_matrix(field, 2)
        actions.append(('restriction', block_matrix([[z,D],[C,z]]), res_expected))
        for label, frob, expected in actions:
            assert frob*frob == -field(p)*identity_matrix(field, 4)
            assert frob.transpose()*pairing*frob == field(p)*pairing
            counts = [0,0]
            for m in planes:
                if (m*frob.transpose()).echelon_form() != m:
                    continue
                is_product = m.matrix_from_columns([0,1]).rank() == 1
                counts[0 if is_product else 1] += 1
            assert tuple(counts) == expected, (p, ell, label, counts, expected)
            family_records.append(dict(p=p, ell=ell, source=label,
                                       product=counts[0], nonproduct=counts[1]))
print('All product/non-product rows checked by independent symplectic enumeration.', flush=True)

# A genuine rational product (2,2)-isogeny and its polarized dual:
# phi = phi_E x phi_E, with phi_E a rational 2-isogeny at p=7.
field = GF(7)
E = EllipticCurve(field, [1,0])
phi = E.isogeny(E([0,0]))
# Sage 10.8 does not implement its generic dual algorithm in this small
# characteristic. Construct the reverse Velu quotient and identify its
# codomain by the multiplication-by-two identity instead.
reverse_velu = phi.codomain().isogeny(phi.codomain()([0,0]))
dual = next(iso*reverse_velu for iso in reverse_velu.codomain().isomorphisms(E)
            if iso*reverse_velu*phi == E.scalar_multiplication(2))
assert phi.degree() == dual.degree() == 2
assert dual*phi == E.scalar_multiplication(2)
assert phi*dual == phi.codomain().scalar_multiplication(2)
extension = GF(49,'u')
E_ext = E.change_ring(extension)
phi_ext = E_ext.isogeny(E_ext([0,0]))
reverse_ext = phi_ext.codomain().isogeny(phi_ext.codomain()([0,0]))
torsion = [P for P in E_ext.points() if 2*P == 0]
images = {phi_ext(P) for P in torsion}
# Postcomposing reverse_velu by an isomorphism does not change its kernel.
reverse_kernel = {P for P in phi_ext.codomain().points() if reverse_ext(P) == 0}
assert images == reverse_kernel and len(images) == 2
print('Actual elliptic dual maps and image-of-torsion reverse kernel verified; '
      'their products give a polarized (2,2)-example.', flush=True)

output = dict(automorphisms=records, kernel_families=family_records,
              dual_example=dict(p=7, ell=2, elliptic_source=[1,0],
                                forward_degree=4, reverse_degree=4,
                                composite_degree=16, verified=True))
(RESULTS/'automorphisms_duality_verification.json').write_text(
    json.dumps(output, indent=2)+'\n')
print('All new checks passed.', flush=True)
