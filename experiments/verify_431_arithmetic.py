"""Independent arithmetic checks for the complete p = 431 rational census.

Run with ordinary Python 3, without Sage or third-party packages:
    python verify_431_arithmetic.py

Scope: imaginary quadratic field class numbers are computed from primitive
reduced binary quadratic forms. The generalized Bernoulli number is computed
directly from its quadratic character sum. These calculations do not read the
graph or another verification report. They evaluate IK1994 Theorem 2/Remark 3
and XY2024 Theorem 1.1, rather than proving the cited formulas.

The separate graph checks read the corrected census. Curve equations, rather
than the stored type labels, give the number of rational branch points. In the
fixed maximal block all remaining branch points have degree two. Automorphism
orders and edges are inputs from the census, not independently reconstructed
here. Consequently the mass and component comparisons corroborate that census
and do not replace its completeness proof or its explicit quotient checks.

Mass normalization: Hashimoto--Koseki (1989), first paper, Theorem 5.6,
equation (5.5), printed p. 19, gives M_O = h(-p)(p-1)/48 for the rank-two
unimodular maximal-order genus when p == 7 mod 8. At two, the four R = Z[sqrt(-p)]
types have stabilizer indices 1, 3, 3, 1 relative to GL_2(Z_2). They correspond
to O^2, O plus R, nonalternating R^2, and alternating R^2, respectively.
The mod-two groups are GL_2(F_2), a line stabilizer, O_2(F_2, I), and
Sp_2(F_2). Summing their masses gives 8 M_O = h(-p)(p-1)/6. The local types
are those classified in Yu (2012), Theorem 1.1 and Section 4.
"""

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt
from pathlib import Path
import json


P = 431
RESULTS = Path(__file__).resolve().parent / "results"
TYPES = ("1^6", "2 1^4", "2^2 1^2", "2^3")
LOCAL_INDICES = dict(zip(TYPES, (1, 3, 3, 1)))


def reduced_forms(discriminant):
    """Primitive positive definite forms, with the usual boundary convention."""
    assert discriminant < 0 and discriminant % 4 in (0, 1)
    forms = []
    for a in range(1, isqrt(-discriminant // 3) + 1):
        for b in range(-a, a + 1):
            numerator = b*b - discriminant
            if numerator % (4*a):
                continue
            c = numerator // (4*a)
            if c < a or gcd(gcd(a, b), c) != 1:
                continue
            if (abs(b) == a or a == c) and b < 0:
                continue
            forms.append((a, b, c))
    return forms


def jacobi(a, n):
    assert n > 0 and n % 2 == 1
    a %= n
    result = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                result = -result
        a, n = n, a
        if a % 4 == n % 4 == 3:
            result = -result
        a %= n
    return result if n == 1 else 0


def kronecker_positive_denominator(a, n):
    assert n > 0
    result = 1
    while n % 2 == 0:
        n //= 2
        if a % 2 == 0:
            return 0
        if a % 8 in (3, 5):
            result = -result
    return result * jacobi(a, n)


def arithmetic_checks():
    forms = {d: reduced_forms(d) for d in (-P, -8*P, -12*P)}
    hp, h2p, h3p = (len(forms[d]) for d in (-P, -8*P, -12*P))
    assert (hp, h2p, h3p) == (21, 8, 24)
    D = 4*P
    character = [kronecker_positive_denominator(D, a)
                 for a in range(1, D + 1)]
    assert sum(character) == 0
    assert sum(a*chi for a, chi in enumerate(character, 1)) == 0
    square_sum = sum(a*a*chi for a, chi in enumerate(character, 1))
    B = Fraction(square_sum, D)
    c, epsilon3 = jacobi(2, P), jacobi(P, 3)
    tau = (B/96 + Fraction(h2p, 8) + Fraction(h3p, 12)
           + (Fraction((P-1)*(9-4*c), 48) + Fraction(P-c, 16)
              + Fraction((1-epsilon3)*(3-c), 12))*hp)
    nplus = B/48 + Fraction((11-3*c)*hp, 8) + Fraction(h3p, 6)
    s = hp  # p == 7 mod 8, the prime-field geometric elliptic count.
    t = P//12 + {1: 0, 5: 1, 7: 1, 11: 2}[P % 12]
    fixed_decomposed = Fraction(t+s*s, 2)
    fixed_jacobians = tau - fixed_decomposed
    split, restriction = s*(2*s+1), Fraction(t+s, 2)
    jplus = nplus - restriction
    # All generic fixed geometric Jacobians have two forms across the two blocks.
    uncorrected_jminus = 2*fixed_jacobians - jplus
    maximal_order_mass = Fraction((P-1)*hp, 48)
    total_mass = sum(LOCAL_INDICES.values())*maximal_order_mass

    # A primitive norm 2^r element generates q^r or conjugate(q)^r. For even
    # norm, primitivity at two forces odd a,b in (a+b sqrt(-p))/2. This excludes
    # rational powers of two, whose ideals contain both primes above two.
    norm_witnesses = {}
    for r in range(1, 8):
        target = 4*(2**r)
        witnesses = []
        for b in range(1, isqrt(target // P) + 1, 2):
            square = target - P*b*b
            a = isqrt(square)
            if a*a == square and a % 2:
                witnesses.append([a, b])
        norm_witnesses[str(r)] = witnesses
    q_order = next(int(r) for r, witnesses in norm_witnesses.items() if witnesses)
    h_order = hp*2*(1-Fraction(kronecker_positive_denominator(-P, 2), 2))
    assert h_order == hp and hp % 2 == 1 and q_order == 7
    expected_components = hp // q_order
    assert (B, tau, nplus, fixed_decomposed, fixed_jacobians, jplus) == (
        6816, 1586, 167, 239, 1347, 138)
    assert total_mass == 1505 and uncorrected_jminus == 2556
    return {
        "scope": "Computed without reading any census or prior arithmetic report.",
        "field_class_numbers": {str(d): len(f) for d, f in forms.items()},
        "reduced_forms": {str(d): f for d, f in forms.items()},
        "real_discriminant": D,
        "quadratic_character_square_sum": square_sum,
        "generalized_B2": str(B),
        "legendre_2_p": c,
        "legendre_p_3": epsilon3,
        "s": s,
        "t": t,
        "split_count": split,
        "restriction_count": int(restriction),
        "geometric_fixed_PPAS_tau": int(tau),
        "geometric_fixed_decomposed": int(fixed_decomposed),
        "geometric_fixed_Jacobians": int(fixed_jacobians),
        "minimal_block_actual_PPAS_Nplus": int(nplus),
        "minimal_block_actual_Jacobians": int(jplus),
        "maximal_Jacobians_before_exceptional_correction": int(uncorrected_jminus),
        "jacobian_asymptotic_main_term": str(Fraction(P*s, 3)),
        "maximal_order_genus_mass": str(maximal_order_mass),
        "conductor_two_local_mass_indices": LOCAL_INDICES,
        "expected_total_mass": str(total_mass),
        "order_R_class_number": int(h_order),
        "prime_above_two_order": q_order,
        "primitive_norm_power_of_two_witnesses": norm_witnesses,
        "expected_Steinitz_components": expected_components,
        "required_total_count_divisor": hp,
        "required_Jacobian_count_residue_mod_h": (-split-int(restriction)) % hp,
    }


def number_of_roots(coefficients):
    total = 0
    for x in range(P):
        value = 0
        for coefficient in reversed(coefficients):
            value = (value*x + coefficient) % P
        total += value == 0
    return total


def graph_checks(graph, arithmetic):
    assert graph["p"] == P and graph["ell"] == 2
    assert graph["weil_polynomial"] == [P*P, 0, 2*P, 0, 1]
    vertices = {v["id"]: v for v in graph["vertices"]}
    full_two = {}
    for e in graph["elliptic_classes"]:
        roots = number_of_roots([e["b"], e["a"], 0, 1])
        assert roots in (1, 3)
        full_two[e["id"]] = roots == 3
    counts, types = Counter(), Counter()
    masses, type_mass = defaultdict(Fraction), defaultdict(Fraction)
    jacobian_aut = Counter()
    for v in vertices.values():
        category = v["type"]
        counts[category] += 1
        mass = Fraction(1, v["aut"])
        masses[category] += mass
        if category == "split":
            n = sum(full_two[i] for i in v["elliptic_ids"])
            typ = {0: "2^2 1^2", 1: "2 1^4", 2: "1^6"}[n]
        elif category == "restriction":
            typ = "2^3"
        else:
            assert category == "Jacobian"
            coefficients = list(v["coefficients"])
            while coefficients and coefficients[-1] % P == 0:
                coefficients.pop()
            degree = len(coefficients)-1
            assert degree in (5, 6)
            fixed = number_of_roots(coefficients) + (degree == 5)
            typ = {0: "2^3", 2: "2^2 1^2", 4: "2 1^4", 6: "1^6"}[fixed]
            sigma = v["frobenius_permutation"]
            assert sorted(sigma) == list(range(6))
            assert all(sigma[sigma[i]] == i for i in range(6))
            assert sum(sigma[i] == i for i in range(6)) == fixed
            jacobian_aut[v["aut"]] += 1
        types[typ] += 1
        type_mass[typ] += mass
        assert v["degree"] == {"1^6": 15, "2 1^4": 3,
                                "2^2 1^2": 3, "2^3": 7}[typ]

    adjacency = {i: set() for i in vertices}
    degrees = Counter()
    for edge in graph["edges"]:
        i, j, multiplicity = (edge[k] for k in ("source", "target", "multiplicity"))
        assert multiplicity > 0 and i in vertices and j in vertices
        adjacency[i].add(j)
        adjacency[j].add(i)
        degrees[i] += multiplicity
    assert all(degrees[i] == v["degree"] for i, v in vertices.items())
    unseen = set(vertices)
    components = []
    while unseen:
        start = min(unseen)
        unseen.remove(start)
        component, stack = [start], [start]
        while stack:
            for neighbor in adjacency[stack.pop()]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    component.append(neighbor)
                    stack.append(neighbor)
        components.append({
            "vertices": len(component),
            "counts": dict(Counter(vertices[i]["type"] for i in component)),
            "directed_kernels": sum(degrees[i] for i in component),
            "mass": str(sum((Fraction(1, vertices[i]["aut"]) for i in component),
                            Fraction())),
        })

    total = len(vertices)
    total_mass = sum(masses.values(), Fraction())
    h = arithmetic["required_total_count_divisor"]
    assert total % h == 0
    assert counts["Jacobian"] % h == arithmetic["required_Jacobian_count_residue_mod_h"]
    assert len(components) == arithmetic["expected_Steinitz_components"]
    assert len({c["vertices"] for c in components}) == 1
    assert all(c["counts"]["split"] == arithmetic["split_count"]//len(components)
               for c in components)
    assert counts == {"split": 903, "restriction": 29, "Jacobian": 3436}
    assert dict(types) == dict(zip(TYPES, (609, 1575, 1575, 609)))
    for typ, index in LOCAL_INDICES.items():
        assert type_mass[typ] == index*Fraction(arithmetic["maximal_order_genus_mass"])
    assert total_mass == Fraction(arithmetic["expected_total_mass"])
    return {
        "scope": "Types recomputed from equations. Automorphism orders and edges read from census.",
        "counts": dict(counts),
        "total_vertices": total,
        "Frobenius_type_counts": {typ: types[typ] for typ in TYPES},
        "Frobenius_type_masses": {typ: str(type_mass[typ]) for typ in TYPES},
        "type_masses": {typ: str(mass) for typ, mass in masses.items()},
        "Jacobian_automorphism_order_counts": dict(sorted(jacobian_aut.items())),
        "total_mass": str(total_mass),
        "total_divided_by_h": total//h,
        "components_recomputed_from_edges": components,
        "directed_kernels": sum(degrees.values()),
        "Jacobian_exceptional_correction": counts["Jacobian"] - arithmetic[
            "maximal_Jacobians_before_exceptional_correction"],
        "checks_passed": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", type=Path, default=RESULTS/"p431_rational_graph.json")
    parser.add_argument("--output", type=Path,
                        default=RESULTS/"p431_arithmetic_verification.json")
    args = parser.parse_args()
    independent = arithmetic_checks()
    graph_bytes = args.graph.read_bytes()
    census = graph_checks(json.loads(graph_bytes), independent)
    result = {"p": P, "graph_file": args.graph.name,
              "graph_sha256": sha256(graph_bytes).hexdigest(),
              "independent_arithmetic": independent, "census_comparisons": census}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print("p=431: h=21, B2=6816, tau=1586, Nplus=167, fixed geometric Jacobians=1347")
    print("Census: S=903, R=29, J=3436, total=4368; mass=1505; components=3 x 1456")
    print(f"All checks passed. Saved {args.output}")


if __name__ == "__main__":
    main()
