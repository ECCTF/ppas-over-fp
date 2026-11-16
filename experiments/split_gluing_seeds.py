"""Exact rational Jacobian neighbours of every split maximal PPAS.

Census helper with all rational twists retained. Formula: Howe--Leprevost--Poonen (2000),
Proposition 4, https://math.mit.edu/~poonen/papers/large.pdf (p. 8).
The negative sign and the full scalar are retained before normalization
by a square, so the outputs are actual Fp forms, not geometric labels.
"""
import argparse
import itertools
import json
from pathlib import Path
import sys
import time

EXPERIMENTS = Path(__file__).resolve().parent
ROOT = EXPERIMENTS.parent.parent
OUT = EXPERIMENTS / 'results'
sys.path.insert(0, str(EXPERIMENTS))
from ten_bit_experiments import Field, Polys, trace_fp


def elliptic_data(p):
    fn = EXPERIMENTS / 'results' / f'p{p}_products.json'
    if fn.exists():
        rows = json.loads(fn.read_text())['elliptic_classes']
    else:
        fn = EXPERIMENTS / 'results' / f'p{p}_maximal_graph.json'
        rows = json.loads(fn.read_text())['elliptic_models']
        rows = [dict(v, id=i) for i, v in enumerate(rows)]
    assert all(trace_fp(p, e['a'], e['b']) == 0 for e in rows)
    assert len(rows) == 2 * len({e['j'] for e in rows})
    for i, e in enumerate(rows):
        for d in rows[:i]:
            if e['j'] == d['j']:
                assert not any(e['a'] == d['a'] * pow(u, 4, p) % p
                               and e['b'] == d['b'] * pow(u, 6, p) % p
                               for u in range(1, p))
    return rows


def glue(F, P, aa, bb, rational=True):
    """HLP Prop. 4, with beta_i matched to alpha_i. None = decomposed."""
    sub, mul, div = F.sub, F.mul, F.div
    cyc = [(2, 1), (1, 0), (0, 2)]
    da = [sub(aa[i], aa[j]) for i, j in cyc]
    db = [sub(bb[i], bb[j]) for i, j in cyc]
    a1 = b1 = 0
    for a, b in zip(da, db):
        a1 = F.add(a1, div(mul(a, a), b))
        b1 = F.add(b1, div(mul(b, b), a))
    a2 = F.add(F.add(mul(aa[0], sub(bb[2], bb[1])),
                         mul(aa[1], sub(bb[0], bb[2]))),
                         mul(aa[2], sub(bb[1], bb[0])))
    b2 = F.add(F.add(mul(bb[0], sub(aa[2], aa[1])),
                         mul(bb[1], sub(aa[0], aa[2]))),
                         mul(bb[2], sub(aa[1], aa[0])))
    if not a2:
        assert not a1 and not b1 and not b2
        # The root bijection is affine, hence induced by a geometric
        # elliptic isomorphism and its quotient is decomposed.
        return None
    assert a1 and b1 and b2
    df = dg = 1
    for a, b in zip(da, db):
        df = mul(df, mul(a, a))
        dg = mul(dg, mul(b, b))
    A = mul(dg, div(a1, a2))
    B = mul(df, div(b1, b2))
    if rational:
        assert A < F.p and B < F.p
    h = [F.neg(1)]
    # da = (alpha3-alpha2, alpha2-alpha1, alpha1-alpha3).
    # Proposition 4 uses their products (1,2), (0,1), (2,0).
    for i, j in [(1, 2), (0, 1), (2, 0)]:
        h = P.mul(h, [mul(B, mul(db[i], db[j])), 0,
                      mul(A, mul(da[i], da[j]))])
    assert len(h) == 7
    if rational:
        assert all(c < F.p for c in h)
    assert P.gcd(h, [mul(i % F.p, h[i]) for i in range(1, 7)]) == [1]
    if not rational:
        return h
    # Only multiply the equation by a rational square.
    target = 1 if pow(h[-1], (F.p - 1)//2, F.p) == 1 else F.nu
    ratio = target * pow(h[-1], -1, F.p) % F.p
    assert pow(ratio, (F.p - 1)//2, F.p) == 1
    return [c * ratio % F.p for c in h]


def generate_seeds(p, verify_fp=True):
    started = time.time()
    F = Field(p)
    P = Polys(F)
    rows = elliptic_data(p)
    roots = [P.roots([e['b'], e['a'], 0, 1]) for e in rows]
    perms = [[r.index(F.frob(x)) for x in r] for r in roots]
    result = {'p': p, 'extension_nonresidue': F.nu,
              'elliptic_classes': rows, 'kernels': [],
              'decomposed_graph_kernels': [], 'split_rows': []}
    fp_checked = set()
    chi = [0] + [1 if pow(i, (p-1)//2, p) == 1 else -1 for i in range(1, p)]
    for i, j in itertools.combinations_with_replacement(range(len(rows)), 2):
        graph_count = smooth_count = 0
        for perm in itertools.permutations(range(3)):
            if any(perm[perms[i][a]] != perms[j][perm[a]] for a in range(3)):
                continue
            graph_count += 1
            h = glue(F, P, roots[i], [roots[j][a] for a in perm])
            rec = {'split_ids': [i, j], 'root_permutation': list(perm)}
            if h is None:
                result['decomposed_graph_kernels'].append(rec)
                continue
            smooth_count += 1
            if verify_fp and tuple(h) not in fp_checked:
                total = p + 1 + chi[h[-1]]
                for x in range(p):
                    y = 0
                    for c in reversed(h):
                        y = (y*x + c) % p
                    total += chi[y]
                assert total == p + 1, (p, i, j, h, total)
                fp_checked.add(tuple(h))
            rec['coefficients'] = h
            result['kernels'].append(rec)
        dims = [sum(x < p for x in roots[k]) for k in [i,j]]
        expected = 6 if dims == [3,3] else 2 if dims == [1,1] else 0
        assert graph_count == expected
        result['split_rows'].append({'split_ids': [i, j],
                                    'product_kernels': dims[0]*dims[1],
                                    'graph_kernels': graph_count,
                                    'Jacobian_kernels': smooth_count})
    result['distinct_equations'] = len({tuple(k['coefficients']) for k in result['kernels']})
    result['seconds'] = round(time.time()-started, 3)
    result['source'] = 'Howe--Leprevost--Poonen (2000), Proposition 4'
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('primes', type=int, nargs='*', default=[7,31,431])
    args = parser.parse_args()
    for p in args.primes:
        result = generate_seeds(p)
        fn = OUT / f'p{p}_split_gluing_seeds.json'
        fn.write_text(json.dumps(result, indent=2) + '\n')
        print(p, 'Jacobian kernels', len(result['kernels']),
              'decomposed graph kernels', len(result['decomposed_graph_kernels']),
              'distinct equations', result['distinct_equations'],
              'seconds', result['seconds'], flush=True)
