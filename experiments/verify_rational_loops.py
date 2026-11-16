#!/usr/bin/env python3
"""Independent loop check at p=31, ell=2, for the square with j=1728.

Only labels, elliptic equations, stored automorphism orders, and the final
adjacency row are read from p31_maximal_graph.json.  The counts of geometric
and rational loops are obtained independently by enumerating integral norm
matrices, with no quotient routine from the graph implementation imported.

Mathematical bounds used by this finite enumeration:
* End_Fp(E) lies in O_Q(sqrt(-31)), where an element outside Z has norm at
  least (31+1)/4=8.  Every entry of M^dagger M=2I has norm at most two,
  so a rational loop matrix has integer entries.
* For geometric j=1728 endomorphisms, the standard quaternion order has
  norm (x+z/2)^2+(y+w/2)^2+31(z^2+w^2)/4.  Norm at most two forces z=w=0,
  so the geometric enumeration takes place in Z[i].
Maps having the same kernel differ by a polarized codomain automorphism.
Consequently the norm-two matrix counts are divided by the corresponding
norm-one matrix counts.  Direct enumeration of their kernels on E[2]^2
also confirms this quotient count.

The product quotient and the exchanged descent of the second graph kernel
identify all three codomains.  Thus the existing adjacency row is a final
comparison, not an input to the independently computed loop count.
"""

from itertools import product
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
P = 31
ZERO = (0, 0)
ONE = (1, 0)
I = (0, 1)


def add(x, y):
    return x[0] + y[0], x[1] + y[1]


def mul(x, y):
    return x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0]


def conj(x):
    return x[0], -x[1]


def norm(x):
    return x[0] ** 2 + x[1] ** 2


def matmul(x, y):
    a, b, c, d = x
    e, f, g, h = y
    return add(mul(a, e), mul(b, g)), add(mul(a, f), mul(b, h)), \
        add(mul(c, e), mul(d, g)), add(mul(c, f), mul(d, h))


def norm_matrices(multiplier, gaussian):
    entries = [(a, b) for a, b in product(range(-1, 2), repeat=2)
               if a * a + b * b <= multiplier and (gaussian or b == 0)]
    return [(a, b, c, d) for a, b, c, d in product(entries, repeat=4)
            if norm(a) + norm(c) == multiplier
            and norm(b) + norm(d) == multiplier
            and add(mul(conj(a), b), mul(conj(c), d)) == ZERO]


def two_torsion_action(z, v):
    # E[2] basis P=(i,0), Q=(-i,0).  The automorphism i exchanges P,Q.
    a, b = z
    x, y = v
    return (a * x + b * y) % 2, (b * x + a * y) % 2


def kernel(matrix):
    a, b, c, d = matrix
    result = []
    for v in product(range(2), repeat=4):
        first = add(two_torsion_action(a, v[:2]), two_torsion_action(b, v[2:]))
        second = add(two_torsion_action(c, v[:2]), two_torsion_action(d, v[2:]))
        if all(n % 2 == 0 for n in first + second):
            result.append(v)
    return tuple(sorted(result))


def frobenius(v):
    # p=31 sends i to -i, hence exchanges the two nonrational branch points.
    return v[1], v[0], v[3], v[2]


def points(a, b):
    result = 1
    for x in range(P):
        rhs = (x ** 3 + a * x + b) % P
        result += 1 if rhs == 0 else (2 if pow(rhs, (P - 1) // 2, P) == 1 else 0)
    return result


def pow_mod(z, exponent):
    result = ONE
    for _ in range(exponent):
        result = tuple(n % P for n in mul(result, z))
    return result


def main():
    data = json.loads((ROOT / 'results/p31_maximal_graph.json').read_text())
    assert data['p'] == P and data['ell'] == 2
    vertices = data['vertices']
    labels = {v['label']: i for i, v in enumerate(vertices)}
    s18, s20, r23 = (vertices[labels[s]] for s in ('S18', 'S20', 'R23'))
    assert s18['elliptic_ids'] == [4, 4] and s20['elliptic_ids'] == [5, 5]
    e0, e1 = data['elliptic_models'][4:6]
    assert e0 == {'a': 1, 'b': 0, 'j': 23}
    assert e1 == {'a': 3, 'b': 0, 'j': 23}
    assert r23['type'] == 'restriction' and r23['j'] == 23
    assert 1728 % P == 23 and points(1, 0) == points(3, 0) == P + 1

    # Exact norm bounds justify completeness, rather than an arbitrary search cutoff.
    assert (P + 1) // 4 == 8 > 2 and P > 4 * 2
    rational_aut = norm_matrices(1, False)
    geometric_aut = norm_matrices(1, True)
    rational_maps = norm_matrices(2, False)
    geometric_maps = norm_matrices(2, True)
    assert (len(rational_aut), len(geometric_aut)) == (8, 32)
    assert (len(rational_maps), len(geometric_maps)) == (8, 96)
    rational_kernels = set(map(kernel, rational_maps))
    geometric_kernels = set(map(kernel, geometric_maps))
    assert len(rational_kernels) == len(rational_maps) // len(rational_aut) == 1
    assert len(geometric_kernels) == len(geometric_maps) // len(geometric_aut) == 3
    assert all(len(k) == 4 for k in geometric_kernels)
    assert all(tuple(sorted(map(frobenius, k))) == k for k in geometric_kernels)

    identity_graph = tuple(sorted((x, y, x, y) for x, y in product(range(2), repeat=2)))
    i_graph = tuple(sorted((x, y, y, x) for x, y in product(range(2), repeat=2)))
    product_kernel = tuple(sorted((x, x, y, y) for x, y in product(range(2), repeat=2)))
    assert geometric_kernels == {identity_graph, i_graph, product_kernel}
    assert rational_kernels == {identity_graph}

    # Product kernel: the elliptic quotient is y^2=x^3-4x.
    quotient_a = (-4) % P
    quotient_to_e1 = [u for u in range(1, P) if pow(u, 4, P) * quotient_a % P == 3]
    assert quotient_to_e1 == [14, 17]
    assert not any(pow(u, 4, P) == 3 for u in range(1, P))
    # Over Fp2=Fp[i], i^2=-1, the two elliptic forms are isomorphic.
    extension_root = next(z for z in product(range(P), repeat=2) if pow_mod(z, 4) == (3, 0))

    # i-graph quotient map M=[[1,i],[i,1]] has norm two.
    # Transported target Frobenius is (i*swap)*pi, since M*bar(M)^(-1)=i*swap.
    # Thus the canonical polarized factors are exchanged, identifying R23.
    m = (ONE, I, I, ONE)
    assert m in geometric_maps and kernel(m) == i_graph
    m_squared = matmul(m, m)
    assert m_squared == (ZERO, (0, 2), (0, 2), ZERO)
    swap = (ZERO, ONE, ONE, ZERO)
    restriction_aut = [g for g in geometric_aut
                       if matmul(matmul(swap, tuple(map(conj, g))), swap) == g]
    assert len(restriction_aut) == 8
    assert s18['aut'] == s20['aut'] == r23['aut'] == 8

    predicted_row = {'S18': 1, 'S20': 1, 'R23': 1}
    stored_row = {vertices[j]['label']: n for j, n in enumerate(data['adjacency'][labels['S18']]) if n}
    assert stored_row == predicted_row
    assert data['adjacency'][labels['S20']][labels['S20']] == 1

    result = {
        'p': P, 'ell': 2,
        'reused_data': 'p31_maximal_graph.json: labels, equations, aut orders, and final adjacency comparison',
        'independence': 'Loop counts use complete norm-bounded integer/Gaussian matrix enumeration and direct two-torsion kernels; no graph quotient routine is imported.',
        'norm_bounds': {'rational_noninteger_minimum': 8, 'geometric_outside_gaussian_lower_bound': '31/4', 'entry_norm_maximum': 2},
        'matrix_counts': {'rational_norm_one': 8, 'geometric_norm_one': 32, 'rational_norm_two': 8, 'geometric_norm_two': 96},
        'loop_kernel_counts': {'rational_codomain_same_Fp_form': 1, 'geometric_codomain_same_class': 3, 'Frobenius_stable_geometric_loop_kernels': 3},
        'models': {'S18': 'E0^2, E0: y^2=x^3+x', 'S20': 'E1^2, E1: y^2=x^3+3x', 'R23': 'Res_Fp2/Fp(E0 over Fp2)', 'j': 23, 'Fp_point_counts_of_elliptic_factors': [32, 32], 'Fp_aut_orders': [8, 8, 8]},
        'quotient_identification': {
            'graph_identity': 'S18, via the rational norm-two matrix [[1,1],[-1,1]]',
            'product_kernel': {'codomain': 'S20', 'elliptic_quotient_a': quotient_a, 'Fp_scalings_to_E1': quotient_to_e1},
            'graph_i': {'codomain': 'R23', 'geometric_quotient_map': '[[1,i],[i,1]]', 'target_Frobenius': 'i*swap*pi', 'reason': 'Frobenius exchanges the two canonical polarized factors'},
            'E0_to_E1_over_Fp2': {'field': 'Fp[i], i^2=-1', 'u': extension_root, 'u_fourth_power': [3, 0]}
        },
        'predicted_S18_row': predicted_row, 'stored_S18_row': stored_row,
        'stored_S20_self_loops': 1, 'checks_passed': True,
        'limitation': 'This verifies one complete degree-three row and the rational loop count at two j=1728 squares, not the full p=31 adjacency matrix.'
    }
    output = ROOT / 'results/rational_loop_verification.json'
    output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
