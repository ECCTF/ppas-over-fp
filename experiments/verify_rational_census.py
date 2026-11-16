"""Independent point counts and local checks for a rational Richelot census.

The generator is not imported. This audit uses saved equations and edges,
and does not by itself certify global completeness or reconstruct quotients.
"""
from collections import Counter
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import sys
import time

import numpy as np

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE))
from verify_ten_bit_results import PointCounter


def matchings(items):
    if not items:
        yield []
        return
    a = items[0]
    for b in items[1:]:
        for rest in matchings([i for i in items if i not in (a, b)]):
            yield [(a, b)] + rest


MATCHINGS = list(matchings(list(range(6))))


def cartier_manin(f, p):
    n = (p - 1) // 2
    out = np.array([1], dtype=np.int64)
    base = np.array(f, dtype=np.int64)
    while n:
        if n & 1:
            out = np.convolve(out, base)[:2*p] % p
        n >>= 1
        if n:
            base = np.convolve(base, base)[:2*p] % p
    return [int(out[i]) if i < len(out) else 0
            for i in (p-1, p-2, 2*p-1, 2*p-2)]


def verify(path):
    start = time.time()
    raw = path.read_bytes()
    data = json.loads(raw)
    p, nu = data['p'], data['field']['nu']
    counter = PointCounter(p, nu)
    vertices = data['vertices']
    degrees = Counter()
    cm_counts = Counter()
    for k, v in enumerate(vertices):
        f = v['coefficients']
        n1 = counter.fp_points(f)
        n2 = counter.quadratic_points(f)
        assert (n1, n2) == (p+1, p*p+1+4*p), (k, f, n1, n2)
        cm = cartier_manin(f, p)
        assert cm == [0, 0, 0, 0], (k, cm)
        roots = v['roots']
        assert len(set(roots)) == 6
        for z in roots:
            if z is None:
                assert len(f) == 6
                continue
            a, b = z % p, z // p
            real = imag = 0
            for c in reversed(f):
                real, imag = (real*a + nu*imag*b+c) % p, (real*b+imag*a) % p
            assert (real, imag) == (0, 0), (k, z)
        perm = v['frobenius_permutation']
        assert all(perm[perm[i]] == i for i in range(6))
        for i, z in enumerate(roots):
            zp = None if z is None else z % p + p * (-(z // p) % p)
            assert roots[perm[i]] == zp
        stable = sum({frozenset((perm[a], perm[b])) for a, b in m}
                     == {frozenset(pair) for pair in m} for m in MATCHINGS)
        assert v['degree'] == stable
        degrees[stable] += 1
        cm_counts[tuple(cm)] += 1
        if (k+1) % 250 == 0:
            print(f'p={p}: checked {k+1}/{len(vertices)} equations, '
                  f'{time.time()-start:.1f}s', flush=True)

    edge = {(e['source'], e['target']): e['multiplicity'] for e in data['edges']}
    assert len(edge) == len(data['edges'])
    row = Counter()
    for (a, b), n in edge.items():
        assert n > 0
        assert n * vertices[b]['aut'] == edge.get((b, a), 0) * vertices[a]['aut']
        row[a] += n
    for e in data['decomposed_edges']:
        assert e['multiplicity'] > 0
        row[e['source']] += e['multiplicity']
        assert e['multiplicity'] * e['aut'] % vertices[e['source']]['aut'] == 0
    assert all(row[v['id']] == v['degree'] for v in vertices)
    report = {
        'p': p, 'input_file': path.name,
        'input_sha256': hashlib.sha256(raw).hexdigest(),
        'Jacobian_equations_verified': len(vertices),
        'points_fp': p+1, 'points_fp2': p*p+1+4*p,
        'zero_Cartier_Manin_matrices': sum(cm_counts.values()),
        'degree_histogram': dict(sorted(degrees.items())),
        'total_Jacobian_kernel_edges': sum(row.values()),
        'weighted_duality_on_all_Jacobian_pairs': True,
        'integral_reverse_decomposed_multiplicities': True,
        'seconds': round(time.time()-start, 3),
        'scope': 'Independent point counts, Cartier-Manin, stable matchings and saved-edge consistency. '
                 'Completeness and rational-isomorphism distinctness use the separately audited census.'
    }
    dest = path.with_name(f'p{p}_rational_census_independent_check.json')
    dest.write_text(json.dumps(report, indent=2)+'\n', encoding='utf8')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('files', type=Path, nargs='+')
    args = parser.parse_args()
    for path in args.files:
        verify(path.resolve())
