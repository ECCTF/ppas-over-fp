"""Exact arithmetic checks for dyadic connectivity when p is 7 modulo 8.

Only the Python standard library is required. Reduced positive binary
quadratic forms give h(Q(sqrt(-p))). Hensel-lifted ideals above two give
the order of [(2,w)], independently checked by integral norm equations.
The resulting component numbers evaluate the manuscript theorem. Only
p=7 and p=31 are compared with complete, previously computed PPAS graphs.

The finite symplectic check verifies the mod-two step used in the
factorization argument. Neither it nor these examples replaces the proof
of the global connectivity theorem.
"""

from collections import Counter
from itertools import combinations
from math import gcd, isqrt
from pathlib import Path
import hashlib
import json


RESULTS = Path(__file__).resolve().parent / "results"


def is_prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p) + 1))


def reduced_forms(p):
    """Primitive reduced positive forms of fundamental discriminant -p."""
    assert is_prime(p) and p % 8 == 7
    forms = []
    for a in range(1, isqrt(p // 3) + 1):
        for b in range(-a, a + 1):
            numerator = b*b + p
            if numerator % (4*a):
                continue
            c = numerator // (4*a)
            if a > c or gcd(a, gcd(b, c)) != 1:
                continue
            if (abs(b) == a or a == c) and b < 0:
                continue
            forms.append((a, b, c))
    assert len(forms) == len(set(forms)) and len(forms) % 2 == 1
    return forms


def reduce_form(a, b, c):
    """Gauss reduction by determinant-one integral changes of basis."""
    discriminant = b*b - 4*a*c
    assert a > 0 and discriminant < 0
    while True:
        shift = (a-b) // (2*a)
        b, c = b + 2*a*shift, c + b*shift + a*shift*shift
        if a > c:
            a, b, c = c, -b, a
            continue
        if a == c and b < 0:
            b = -b
        assert -a < b <= a <= c and b*b - 4*a*c == discriminant
        return a, b, c


def prime_two_order(p, h):
    """q^n=(2^n,w-r_n), where w^2-w+(p+1)/4=0 and q=(2,w)."""
    constant = (p+1) // 4
    modulus, root = 2, 0
    principal = (1, 1, constant)
    powers = []
    for exponent in range(1, h+1):
        numerator = root*root - root + constant
        assert numerator % modulus == 0
        form = reduce_form(modulus, 1-2*root, numerator // modulus)
        powers.append({"exponent": exponent, "ideal_norm": modulus,
                       "hensel_root": root, "reduced_form": form})
        if form == principal:
            assert h % exponent == 0
            assert len({tuple(row["reduced_form"]) for row in powers}) == exponent
            return exponent, powers
        # The derivative 2*r-1 is odd, so exactly one lift is a root.
        candidates = [r for r in (root, root+modulus)
                      if (r*r-r+constant) % (2*modulus) == 0]
        assert len(candidates) == 1
        root = candidates[0]
        modulus *= 2
    raise AssertionError("The ideal power did not become principal within h steps.")


def pure_prime_two_generators(p, exponent):
    """All a+b*w with norm 2^n and ideal (q^n), via a separate norm search.

    At q, w=0 modulo two; at qbar, w=1. Thus a is even and b is odd.
    The norm equation bounds b through (2*a+b)^2+p*b^2=4*2^n.
    """
    norm = 2**exponent
    bound = isqrt((4*norm) // p)
    generators = []
    for b in range(-bound, bound+1):
        if b % 2 == 0:
            continue
        square = 4*norm - p*b*b
        t0 = isqrt(square)
        if t0*t0 != square:
            continue
        for t in sorted({t0, -t0}):
            if (t-b) % 2:
                continue
            a = (t-b) // 2
            if a % 2:
                continue
            assert a*a+a*b+(p+1)//4*b*b == norm
            assert a % 2 == 0 and (a+b) % 2 == 1
            generators.append([a, b])
    return sorted(generators)


def arithmetic_record(p):
    forms = reduced_forms(p)
    h = len(forms)
    order, powers = prime_two_order(p, h)
    assert all(tuple(row["reduced_form"]) in forms for row in powers)
    independent = {n: pure_prime_two_generators(p, n)
                   for n in range(1, order+1)}
    assert all(not independent[n] for n in range(1, order))
    assert independent[order]
    return {
        "p": p, "fundamental_discriminant": -p,
        "class_number": h, "prime_two_class_order": order,
        "component_count_by_theorem": h//order,
        "connected_by_theorem": h == order,
        "reduced_forms": forms, "prime_two_class_powers": powers,
        "independent_norm_search": independent,
        "norm_search_condition": "a even, b odd, (2*a+b)^2+p*b^2=4*2^n",
        "full_PPAS_graph_enumerated_here": False,
    }


def graph_components(adjacency):
    unused = set(range(len(adjacency)))
    components = []
    while unused:
        start = min(unused)
        component, pending = {start}, [start]
        while pending:
            i = pending.pop()
            for j in range(len(adjacency)):
                if (adjacency[i][j] or adjacency[j][i]) and j not in component:
                    component.add(j)
                    pending.append(j)
        assert component <= unused
        unused -= component
        components.append(sorted(component))
    return components


def compare_existing_graph(record):
    p = record["p"]
    path = RESULTS / f"p{p}_maximal_graph.json"
    data = json.loads(path.read_text(encoding="utf8"))
    assert data["p"] == p and data["ell"] == 2
    components = graph_components(data["adjacency"])
    assert len(components) == record["component_count_by_theorem"]
    assert sorted(map(sorted, components)) == sorted(map(sorted, data["components"]))
    return {"p": p, "vertices": len(data["vertices"]),
            "component_sizes": sorted(map(len, components), reverse=True),
            "theorem_index": record["component_count_by_theorem"],
            "input": path.name, "input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "scope": "Independent reachability on the previously computed complete adjacency matrix."}


def pairing(v, w):
    """Standard symplectic pairing on bit coordinates (e1,e2,f1,f2)."""
    return (((v & 3) & (w >> 2)).bit_count()
            + ((v >> 2) & (w & 3)).bit_count()) % 2


def span(vectors):
    result = {0}
    for vector in vectors:
        result |= {x ^ vector for x in result}
    return frozenset(result)


def apply(columns, vector):
    result = 0
    for i, column in enumerate(columns):
        if vector & (1 << i):
            result ^= column
    return result


def symplectic_factorization_check():
    """All coisotropic spaces and all symplectic involutions, including id."""
    ambient = frozenset(range(16))
    planes = {span(basis) for basis in combinations(range(1,16), 2)}
    hyperplanes = {span(basis) for basis in combinations(range(1,16), 3)}
    hyperplanes = {space for space in hyperplanes if len(space) == 8}
    assert len(planes) == 35 and len(hyperplanes) == 15
    lagrangians = {space for space in planes
                  if all(pairing(v,w) == 0 for v in space for w in space)}
    assert len(lagrangians) == 15
    coisotropic = sorted(lagrangians | hyperplanes | {ambient},
                         key=lambda space: (len(space), sorted(space)))
    contained = {}
    counts = Counter()
    for space in coisotropic:
        perpendicular = {v for v in ambient if all(pairing(v,w) == 0 for w in space)}
        assert perpendicular <= space
        contained[space] = [lag for lag in lagrangians if lag <= space]
        dimension = len(space).bit_length()-1
        assert len(contained[space]) == {2:1,3:3,4:15}[dimension]
        counts[dimension] += 1
    assert dict(counts) == {2:15,3:15,4:1}

    symplectic = []
    for a in range(1,16):
        for c in range(1,16):
            if pairing(a,c) != 1:
                continue
            complement = [x for x in range(1,16)
                          if pairing(a,x) == pairing(c,x) == 0]
            for b in complement:
                for d in complement:
                    if pairing(b,d) == 1:
                        symplectic.append((a,b,c,d))
    assert len(symplectic) == len(set(symplectic)) == 720
    basis = (1,2,4,8)
    involutions = [g for g in symplectic
                   if all(apply(g, apply(g, b)) == b for b in basis)]
    assert len(involutions) == 76
    histogram = Counter()
    for g in involutions:
        for space in coisotropic:
            if frozenset(apply(g,w) for w in space) != space:
                continue
            fixed = [lag for lag in contained[space]
                     if frozenset(apply(g,w) for w in lag) == lag]
            assert len(fixed) >= 1 and len(fixed) % 2 == 1
            dimension = len(space).bit_length()-1
            histogram[(dimension, len(fixed))] += 1
    return {
        "symplectic_group_order": len(symplectic),
        "elements_with_square_identity": len(involutions),
        "coisotropic_spaces_by_dimension": dict(counts),
        "contained_lagrangians_by_dimension": {2:1,3:3,4:15},
        "invariant_pairs_by_dimension_and_fixed_lagrangians": [
            {"dimension": dim,"fixed_lagrangians": fixed,"number_of_pairs": count}
            for (dim,fixed),count in sorted(histogram.items())],
        "all_preserved_coisotropic_spaces_have_fixed_lagrangians": True,
        "scope": "Checks the finite mod-two assertion, not the global Hermitian approximation theorem.",
    }


def main():
    primes = [7,23,31,431,503,1319]
    rows = [arithmetic_record(p) for p in primes]
    expected = [(1,1,1),(3,3,1),(3,3,1),(21,7,3),(21,7,3),(45,9,5)]
    assert [(r["class_number"],r["prime_two_class_order"],r["component_count_by_theorem"])
            for r in rows] == expected
    p431 = next(row for row in rows if row["p"] == 431)
    assert [4,1] in p431["independent_norm_search"][7]
    assert not p431["independent_norm_search"][1]
    result = {
        "method": "Reduced primitive quadratic forms, Hensel ideal powers, independent norm equations, and finite symplectic enumeration; exact integers only.",
        "scope": "Component counts evaluate the proved p=7 mod 8 theorem. No PPAS graph or new curve census is generated. Larger examples are class-arithmetic evaluations, not independent full-graph checks.",
        "arithmetic_examples": rows,
        "p431_certificate": {"w_equation": "w^2-w+108=0", "q": "(2,w)",
                             "alpha": "w+4", "norm_alpha": 128,
                             "residues_at_q_and_qbar": [0,1],
                             "principal_ideal": "(alpha)=q^7", "q_class_order": 7,
                             "class_number": 21, "component_count_by_theorem": 3},
        "complete_graph_comparisons": [compare_existing_graph(r) for r in rows if r["p"] in (7,31)],
        "factorization_local_check": symplectic_factorization_check(),
    }
    target = RESULTS/'dyadic_seven_verification.json'
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    for r in rows:
        print(f"p={r['p']}: h={r['class_number']}, ord(q2)={r['prime_two_class_order']}, "
              f"components by theorem={r['component_count_by_theorem']}")
    print("p=7 and p=31 agree with their complete stored graph components.")
    print("All 76 symplectic involutions and all 31 coisotropic spaces passed the finite factorization check.")
    print(f"Wrote {target}")


if __name__ == '__main__':
    main()
