"""Independent arithmetic checks of the saved ten-bit experiment data.

Run with Python 3 and NumPy: python verify_ten_bit_results.py

This script imports neither experiment generator. It enumerates every
rational elliptic j by direct point counting, counts points over Fp^2
using NumPy arrays and the norm character, and enumerates symplectic
planes as triples of nonzero vectors. It also checks the consistency of
the saved kernel records, adjacency matrices, and reverse multiplicities.

The latter checks reuse the saved quotient j values and automorphism
orders. They are not an independent reconstruction of every isogeny,
curve isomorphism test, or complete PPAS census. The auxiliary elliptic
Fp^2 point-count check samples two nonrational j values per prime.
"""
from collections import Counter
from itertools import combinations, combinations_with_replacement
from pathlib import Path
import hashlib
import json
import time

import numpy as np


RESULTS = Path(__file__).resolve().parent / 'results'
PRIMES = (523, 541, 547, 613, 619, 631, 661, 673, 751, 823, 991)
JACOBIAN_PRIMES = (719, 751, 823)
MESSAGES = []


def log(message):
    print(message, flush=True)
    MESSAGES.append(message)


def load(name):
    path = RESULTS / name
    payload = path.read_bytes()
    return json.loads(payload), hashlib.sha256(payload).hexdigest()


class PointCounter:
    """Point enumeration using chi_Fp2(z) = chi_Fp(Norm(z))."""
    def __init__(self, p, nu):
        self.p = p
        self.nu = nu
        assert pow(nu, (p - 1) // 2, p) == p - 1
        self.chi = np.array([0] + [
            1 if pow(i, (p - 1) // 2, p) == 1 else -1
            for i in range(1, p)
        ], dtype=np.int64)
        self.x = np.arange(p, dtype=np.int64)
        self.real_x = np.tile(self.x, p)
        self.imag_x = np.repeat(self.x, p)

    def fp_points(self, coefficients):
        p = self.p
        assert all(0 <= c < p for c in coefficients)
        values = np.zeros(p, dtype=np.int64)
        for c in reversed(coefficients):
            values = (values * self.x + c) % p
        infinity = 1 if len(coefficients) % 2 == 0 else 1 + int(self.chi[coefficients[-1]])
        return p + int(self.chi[values].sum()) + infinity

    def quadratic_points(self, coefficients):
        p, nu = self.p, self.nu
        real = np.zeros(p * p, dtype=np.int64)
        imag = np.zeros(p * p, dtype=np.int64)
        for c in reversed(coefficients):
            real, imag = (
                (real * self.real_x + nu * imag * self.imag_x + c % p) % p,
                (real * self.imag_x + imag * self.real_x + c // p) % p,
            )
        norms = (real * real - nu * imag * imag) % p
        lead = coefficients[-1]
        lead_norm = ((lead % p) ** 2 - nu * (lead // p) ** 2) % p
        infinity = 1 if len(coefficients) % 2 == 0 else 1 + int(self.chi[lead_norm])
        return p * p + int(self.chi[norms].sum()) + infinity


def short_model(p, j):
    if j == 0:
        return 0, 1
    if j == 1728 % p:
        return 1, 0
    ratio = j * pow((1728 - j) % p, -1, p) % p
    return 3 * ratio % p, 2 * ratio % p


def pairing(v, w):
    a = [(v >> i) & 1 for i in range(4)]
    b = [(w >> i) & 1 for i in range(4)]
    return (a[0] * b[1] + a[1] * b[0] + a[2] * b[3] + a[3] * b[2]) % 2


def local_plane_counts(n1, n2):
    # A two-plane has three nonzero vectors a,b,c with a+b+c=0.
    planes = [frozenset((0, a, b, c))
              for a, b, c in combinations(range(1, 16), 3)
              if a ^ b ^ c == 0 and pairing(a, b) == 0]
    assert len(planes) == 15

    def frobenius(v):
        coordinates = [(v >> i) & 1 for i in range(4)]
        if n1 == 1:
            coordinates[0] ^= coordinates[1]
        if n2 == 1:
            coordinates[2] ^= coordinates[3]
        return sum(a << i for i, a in enumerate(coordinates))

    stable = [plane for plane in planes if frozenset(map(frobenius, plane)) == plane]
    product = [plane for plane in stable
               if len(plane & {0, 1, 2, 3}) == 2
               and len(plane & {0, 4, 8, 12}) == 2]
    return len(stable), len(product)


def support_components(adjacency):
    n = len(adjacency)
    remaining, components = set(range(n)), []
    while remaining:
        pending = [min(remaining)]
        component = set()
        while pending:
            i = pending.pop()
            if i in component:
                continue
            component.add(i)
            pending.extend(j for j in remaining - component
                           if adjacency[i][j] or adjacency[j][i])
        remaining -= component
        components.append(sorted(component))
    return sorted(components)


def verify_restrictions(p):
    data, digest = load(f'p{p}_ten_bit_restriction_ell3.json')
    nu = int(data['field_relation'].split('=')[1])
    counter = PointCounter(p, nu)
    rational_j = []
    for j in range(p):
        a, b = short_model(p, j)
        if counter.fp_points([b, a, 0, 1]) == p + 1:
            rational_j.append(j)
    assert rational_j == sorted({e['j'] for e in data['elliptic_fp_forms']})
    assert len(rational_j) == data['s']

    forms = data['elliptic_fp_forms']
    torsion_counts = []
    isomorphism_representatives = set()
    for e in forms:
        a, b = e['a'], e['b']
        assert counter.fp_points([b, a, 0, 1]) == p + 1
        roots = int(np.count_nonzero((counter.x ** 3 + a * counter.x + b) % p == 0))
        assert roots == e['rational_2_kernels']
        torsion_counts.append(roots)
        canonical = min((a * pow(u, 4, p) % p, b * pow(u, 6, p) % p) for u in range(1, p))
        assert canonical not in isomorphism_representatives
        isomorphism_representatives.add(canonical)
    assert len(forms) == 2 * len(rational_j)

    local = {(a, b): local_plane_counts(a, b) for a in (1, 3) for b in (1, 3)}
    histogram, product_count = Counter(), 0
    for i, j in combinations_with_replacement(range(len(forms)), 2):
        degree, product = local[torsion_counts[i], torsion_counts[j]]
        histogram[degree] += 1
        product_count += product
    assert dict(histogram) == {int(k): v for k, v in data['split_ell2_degree_histogram'].items()}
    assert product_count == data['split_ell2_product_kernels']
    assert sum(histogram.values()) == data['split_vertices']

    samples = [e for e in data['elliptic_models'] if e['j'] >= p][:2]
    assert len(samples) == 2
    sample_results = []
    for e in samples:
        count = counter.quadratic_points([e['b'], e['a'], 0, 1])
        assert count in ((p - 1) ** 2, (p + 1) ** 2)
        sample_results.append({'j_encoding': e['j'], 'points_fp2': count,
                               'scalar_frobenius': -p if count == (p + 1) ** 2 else p})

    vertices, adjacency = data['vertices'], data['adjacency']
    n = len(vertices)
    assert n == data['restriction_vertices'] == (data['s'] + data['t']) // 2
    assert len(adjacency) == n and all(len(row) == n for row in adjacency)
    assert all(sum(row) == 4 for row in adjacency)
    orbit_index = {j: i for i, v in enumerate(vertices) for j in v['j_orbit']}
    models = {e['j']: e for e in data['elliptic_models']}
    for i, v in enumerate(vertices):
        kernels = models[v['j_orbit'][0]]['kernels']
        assert len(kernels) == len({e['kernel_x'] for e in kernels}) == 4
        row = Counter(orbit_index[e['target_j']] for e in kernels)
        assert adjacency[i] == [row[j] for j in range(n)]
        for j in range(n):
            assert adjacency[i][j] * vertices[j]['automorphism_order'] == adjacency[j][i] * v['automorphism_order']
    components = support_components(adjacency)
    assert components == sorted(data['components'])
    assert sum(map(sum, adjacency)) == data['directed_kernels']
    log(f'p={p}: all {p} rational j checked, s={len(rational_j)}, '
        f'{len(forms)} distinct elliptic forms, {n} restriction vertices, '
        f'{data["directed_kernels"]} kernels; all checks passed')
    return {'p': p, 'input_sha256': digest, 'all_rational_j_tested': p,
            'supersingular_rational_j': rational_j, 'distinct_fp_elliptic_forms': len(forms),
            'sampled_nonrational_elliptic_models': sample_results,
            'split_degree_histogram': dict(histogram), 'split_product_kernels': product_count,
            'restriction_vertices': n, 'directed_kernels': data['directed_kernels'],
            'component_sizes': sorted(map(len, components), reverse=True),
            'saved_adjacency_matches_saved_kernel_records': True, 'weighted_duality_verified': True}


def verify_jacobians(p):
    data, digest = load(f'p{p}_jacobian_neighborhood.json')
    counter = PointCounter(p, data['extension_nonresidue'])
    results = []
    for v in data['vertices']:
        if v['type'] != 'Jacobian':
            continue
        n1 = counter.fp_points(v['coefficients'])
        n2 = counter.quadratic_points(v['coefficients'])
        assert n1 == v['n1'] == p + 1
        assert n2 == v['n2'] == p * p + 1 + 4 * p
        results.append({'label': v['label'], 'points_fp': n1, 'points_fp2': n2})
    forward = Counter(k['target'] for k in data['kernels'])
    assert sum(forward.values()) == data['source_outdegree'] == 7
    for edge in data['edges']:
        i, j = edge['source'], edge['target']
        assert i == 0 and edge['forward_kernels'] == forward[j]
        a, b = data['vertices'][i]['aut'], data['vertices'][j]['aut']
        assert edge['forward_kernels'] * b == edge['reverse_kernels_by_duality'] * a
        if 'reverse_kernels_directly_verified' in edge:
            assert edge['reverse_kernels_directly_verified'] == edge['reverse_kernels_by_duality']
    log(f'p={p}: independent Fp/Fp2 counts for all {len(results)} saved Jacobian models '
        'and saved reverse-multiplicity consistency passed')
    return {'p': p, 'input_sha256': digest, 'jacobian_point_counts': results,
            'saved_source_kernels': sum(forward.values()),
            'saved_reverse_multiplicity_consistency': True}


def main():
    started = time.time()
    report = {
        'independent_checks': [
            'Complete rational elliptic j census by direct Fp point enumeration',
            'Distinctness of the saved elliptic Fp forms by short-Weierstrass scaling',
            'Actual rational two-torsion roots and all 15 symplectic planes',
            'Direct Fp2 point enumeration for two nonrational auxiliary elliptic models per prime',
            'Direct Fp and Fp2 point enumeration for every saved Jacobian vertex',
        ],
        'consistency_checks_reusing_generator_data': [
            'Adjacency reconstructed from saved quotient-j kernel records',
            'Degrees, weighted duality and components of the saved restriction matrices',
            'Jacobian reverse multiplicities checked against saved automorphism orders',
        ],
        'limitations': [
            'Does not independently recompute every Velu or Richelot quotient',
            'Does not independently repeat the genus-two isomorphism or automorphism computation',
            'Does not provide a full PPAS or Jacobian census at these primes',
            'The I/U Frobenius normal form uses the standard elliptic two-torsion classification',
        ],
        'restriction_experiments': [verify_restrictions(p) for p in PRIMES],
        'jacobian_neighborhoods': [verify_jacobians(p) for p in JACOBIAN_PRIMES],
    }
    report['total_jacobian_models_checked'] = sum(len(r['jacobian_point_counts']) for r in report['jacobian_neighborhoods'])
    assert report['total_jacobian_models_checked'] == 13
    report['seconds'] = round(time.time() - started, 3)
    report['all_checks_passed'] = True
    (RESULTS / 'ten_bit_independent_verification.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf8')
    log(f'ALL CHECKS PASSED in {report["seconds"]} seconds')
    (RESULTS / 'ten_bit_independent_verification.log').write_text('\n'.join(MESSAGES) + '\n', encoding='utf8')


if __name__ == '__main__':
    main()
