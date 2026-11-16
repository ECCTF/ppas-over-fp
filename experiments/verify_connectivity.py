"""Exact connectivity checks for the supplied maximal-block graphs.

Run with standard Python.  This script reuses the complete adjacency
matrices produced by the curve/isogeny census; it does not regenerate
that census.  It independently computes directed reachability, weak
components, and rational Laplacian rank.  Two-torsion Frobenius types
are recovered from the stored equations, not from component labels.
"""

from fractions import Fraction
from itertools import product
from pathlib import Path
import json
import math


RESULTS = Path(__file__).resolve().parent / "results"


def rational_rank(rows):
    matrix = [[Fraction(x) for x in row] for row in rows]
    rank = 0
    for column in range(len(matrix[0]) if matrix else 0):
        pivot = next((i for i in range(rank, len(matrix))
                      if matrix[i][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        scale = matrix[rank][column]
        matrix[rank] = [x / scale for x in matrix[rank]]
        for i in range(rank + 1, len(matrix)):
            scale = matrix[i][column]
            if scale:
                matrix[i] = [a - scale * b
                             for a, b in zip(matrix[i], matrix[rank])]
        rank += 1
    return rank


def reachable(adjacency, start, undirected=False):
    visited = {start}
    pending = [start]
    while pending:
        i = pending.pop()
        for j in range(len(adjacency)):
            linked = adjacency[i][j] or (undirected and adjacency[j][i])
            if linked and j not in visited:
                visited.add(j)
                pending.append(j)
    return visited


def component_partition(adjacency, strong=False):
    n = len(adjacency)
    closure = [reachable(adjacency, i, undirected=not strong)
               for i in range(n)]
    unused = set(range(n))
    components = []
    while unused:
        i = min(unused)
        component = {j for j in closure[i] if not strong or i in closure[j]}
        assert component <= unused
        components.append(sorted(component))
        unused -= component
    return components


def check_graph(adjacency, labels, automorphism_orders):
    n = len(labels)
    assert len(adjacency) == n and all(len(row) == n for row in adjacency)
    assert all(isinstance(x, int) and x >= 0 for row in adjacency for x in row)
    assert len(automorphism_orders) == n and min(automorphism_orders) > 0
    for i in range(n):
        for j in range(n):
            assert bool(adjacency[i][j]) == bool(adjacency[j][i])
            assert (adjacency[i][j] * automorphism_orders[j]
                    == adjacency[j][i] * automorphism_orders[i])
    weak = component_partition(adjacency)
    strong = component_partition(adjacency, strong=True)
    assert weak == strong
    degrees = [sum(row) for row in adjacency]
    laplacian = [[(degrees[i] if i == j else 0) - adjacency[i][j]
                  for j in range(n)] for i in range(n)]
    rank = rational_rank(laplacian)
    assert n - rank == len(weak)
    # The left weight diag(1/a_i) makes the Laplacian symmetric.
    for i in range(n):
        for j in range(n):
            assert (Fraction(laplacian[i][j], automorphism_orders[i])
                    == Fraction(laplacian[j][i], automorphism_orders[j]))
    return {
        "vertex_count": n,
        "directed_kernel_count": sum(degrees),
        "loop_kernel_count": sum(adjacency[i][i] for i in range(n)),
        "outdegrees": degrees,
        "components": [[labels[i] for i in component] for component in weak],
        "component_sizes": sorted((len(c) for c in weak), reverse=True),
        "weak_component_count": len(weak),
        "strong_component_count": len(strong),
        "rational_laplacian_rank": rank,
        "rational_laplacian_nullity": n - rank,
        "connected": len(weak) == 1,
        "reverse_support_and_weighted_balance_verified": True,
    }


def evaluate(coefficients, x, p):
    value = 0
    for coefficient in reversed(coefficients):
        value = (value * x + coefficient) % p
    return value


def elliptic_two_torsion_action(model, p):
    roots = sum((x**3 + model["a"] * x + model["b"]) % p == 0
                for x in range(p))
    assert roots in (1, 3)
    return "I" if roots == 3 else "U"


def recover_two_torsion_types(data):
    p = data["p"]
    elliptic_actions = [elliptic_two_torsion_action(e, p)
                        for e in data["elliptic_models"]]
    types = {}
    for vertex in data["vertices"]:
        if vertex["type"] == "restriction":
            cycle = "2^3"
        elif vertex["type"] == "split":
            i, j = vertex["elliptic_ids"]
            actions = elliptic_actions[i] + elliptic_actions[j]
            cycle = {"II": "1^6", "IU": "2 1^4", "UI": "2 1^4",
                     "UU": "2^2 1^2"}[actions]
        else:
            coefficients = vertex["coefficients"]
            degree = len(coefficients) - 1
            rational_branch_points = sum(evaluate(coefficients, x, p) == 0
                                         for x in range(p)) + (degree == 5)
            cycle = {6: "1^6", 4: "2 1^4", 2: "2^2 1^2",
                     0: "2^3"}[rational_branch_points]
        types[vertex["label"]] = cycle
    expected = {"2^2 1^2", "2^3"}
    if p % 4 == 3:
        expected |= {"1^6", "2 1^4"}
    assert set(types.values()) == expected
    return types


def local_two_torsion_check():
    vectors = list(product(range(2), repeat=4))
    zero = (0, 0, 0, 0)
    def dot(a, b):
        return sum(x*y for x, y in zip(a, b)) % 2
    def pairing(a, b):
        return (a[0]*b[1] + a[1]*b[0] + a[2]*b[3] + a[3]*b[2]) % 2
    def image(matrix, v):
        return tuple(dot(row, v) for row in matrix)
    def quadratic(b, v):
        return (v[0]*v[1] + v[2]*v[3] + dot(b, v)) % 2
    odd_forms = [b for b in vectors if b[0]*b[1] + b[2]*b[3] == 1]
    assert len(odd_forms) == 6
    planes = set()
    for v in vectors[1:]:
        for w in vectors[1:]:
            if v != w and not pairing(v, w):
                planes.add(frozenset((zero, v, w, tuple(x ^ y for x, y in zip(v, w)))))
    assert len(planes) == 15
    identity = [[int(i == j) for j in range(4)] for i in range(4)]
    iu = [row[:] for row in identity]; iu[2][3] = 1
    uu = [row[:] for row in iu]; uu[0][1] = 1
    exchange = [[0,0,1,0], [0,0,0,1], [1,0,0,0], [0,1,0,0]]
    records = []
    for name, matrix, fixed_expected, kernels_expected in [
        ("I+I", identity, 6, 15), ("I+U", iu, 4, 3),
        ("U+U", uu, 2, 3), ("exchanged factors", exchange, 0, 7),
    ]:
        assert all(pairing(image(matrix, v), image(matrix, w)) == pairing(v, w)
                   for v in vectors for w in vectors)
        assert all(image(matrix, image(matrix, v)) == v for v in vectors)
        fixed_forms = sum(all(quadratic(b, image(matrix, v)) == quadratic(b, v)
                              for v in vectors) for b in odd_forms)
        fixed_planes = sum(frozenset(image(matrix, v) for v in plane) == plane
                           for plane in planes)
        assert (fixed_forms, fixed_planes) == (fixed_expected, kernels_expected)
        records.append({"module": name, "fixed_odd_theta_characteristics": fixed_forms,
                        "fixed_lagrangians": fixed_planes})
    return records


def elliptic_point_count(p, a, b):
    return 1 + sum(1 for x in range(p) for y in range(p)
                   if (y*y - x*x*x - a*x - b) % p == 0)


def elliptic_witnesses():
    records = []
    for p in (7, 11, 19, 23, 31):
        for a, action in ((1, "U"), (p-1, "I")):
            model = {"a": a, "b": 0}
            assert elliptic_point_count(p, a, 0) == p + 1
            assert elliptic_two_torsion_action(model, p) == action
            records.append({"p": p, **model, "trace": 0, "two_torsion_action": action})
    for p in (13, 17, 29):
        model = next({"a": a, "b": b} for a in range(p) for b in range(p)
                     if (4*a**3+27*b*b) % p
                     and elliptic_point_count(p, a, b) == p + 1)
        assert elliptic_two_torsion_action(model, p) == "U"
        records.append({"p": p, **model, "trace": 0, "two_torsion_action": "U"})
    return records


def mixed_star_product_checks():
    """Check the combinatorics of h elliptic stars; not a surface census.

    An input star has center C_i (I on two-torsion) and three leaves
    L_(i,a) (U on two-torsion). A mixed surface is the unordered
    product L_(i,a) x C_j. Its three product kernels lead to
    L_(j,b) x C_i for b=0,1,2. The distinct I/U types make this
    vertex label unambiguous despite the unordered product.
    """
    records = []
    for h in (1, 2, 3):
        vertices = [(i, a, j) for i in range(h) for a in range(3) for j in range(h)]
        labels = [f"L{i},{a} x C{j}" for i, a, j in vertices]
        lookup = {v: i for i, v in enumerate(vertices)}
        adjacency = [[0]*len(vertices) for _ in vertices]
        for n, (i, a, j) in enumerate(vertices):
            for b in range(3):
                adjacency[n][lookup[j, b, i]] += 1
        report = check_graph(adjacency, labels, [4]*len(vertices))
        diagonal = 0
        off_diagonal = 0
        for component in component_partition(adjacency):
            pairs = {tuple(sorted((vertices[n][0], vertices[n][2]))) for n in component}
            assert len(pairs) == 1
            i, j = next(iter(pairs))
            if i == j:
                diagonal += 1
                assert len(component) == 3
                assert all(adjacency[a][b] == 1 for a in component for b in component)
            else:
                off_diagonal += 1
                assert len(component) == 6
                left = [n for n in component if vertices[n][0] == i]
                right = [n for n in component if vertices[n][0] == j]
                assert len(left) == len(right) == 3
                for a in component:
                    for b in component:
                        assert adjacency[a][b] == int((a in left) != (b in left))
        assert diagonal == h and off_diagonal == h*(h-1)//2
        assert report["vertex_count"] == 3*h*h
        assert report["weak_component_count"] == h*(h+1)//2
        assert report["directed_kernel_count"] == 9*h*h
        report.update({"h": h, "source": "Combinatorial construction from h disjoint three-leaf elliptic stars",
                       "scope": "Checks the mixed-product proof only; not a new full PPAS graph or a new prime census.",
                       "J3_components": diagonal, "K3_3_components": off_diagonal,
                       "full_surface_component_lower_bound_if_hypotheses_hold": h*(h+1)//2+1})
        records.append(report)

    data = json.loads((RESULTS / "p11_maximal_graph.json").read_text())
    types = recover_two_torsion_types(data)
    indices = [i for i, vertex in enumerate(data["vertices"])
               if vertex["type"] == "split" and types[vertex["label"]] == "2 1^4"]
    assert len(indices) == 3
    adjacency = data["adjacency"]
    assert all(adjacency[i][j] == 1 for i in indices for j in indices)
    assert all(adjacency[i][j] == adjacency[j][i] == 0
               for i in indices for j in range(len(adjacency)) if j not in indices)
    actual = {
        "p": 11, "ell": 2, "h_Q_sqrt_minus_11": 1,
        "mixed_product_labels": [data["vertices"][i]["label"] for i in indices],
        "mixed_product_adjacency": [[adjacency[i][j] for j in indices] for i in indices],
        "closed_component": True, "matches_J3_prediction": True,
        "complement_vertex_count": len(adjacency)-len(indices),
        "full_component_lower_bound": 2,
    }
    return {"generic_combinatorial_checks": records, "actual_p11_check": actual}


def class_arithmetic_431():
    """Verify the class-number and split-prime witnesses without Sage.

    For D=-431, reduced primitive positive binary quadratic forms
    enumerate the ordinary ideal classes of the maximal order. With
    w=(1+sqrt(-431))/2, norms are a^2+a*b+108*b^2. The norm-128
    element w+4 vanishes at exactly one of the two primes over 2.
    This proves q^7 principal, while the absence of norm-2 elements
    proves q nonprincipal and hence of exact order 7.
    """
    discriminant = -431
    assert discriminant % 4 == 1
    assert all(431 % d for d in range(2, math.isqrt(431)+1))
    forms = []
    for a in range(1, math.isqrt((-discriminant)//3)+1):
        for b in range(-a, a+1):
            numerator = b*b-discriminant
            if numerator % (4*a):
                continue
            c = numerator//(4*a)
            if a > c or math.gcd(a, b, c) != 1:
                continue
            if (abs(b) == a or a == c) and b < 0:
                continue
            assert b*b-4*a*c == discriminant
            forms.append([a, b, c])
    assert len(forms) == 21

    def norm_solutions(n):
        # 4*N(a+b*w)=(2*a+b)^2+431*b^2 bounds b exactly.
        solutions = set()
        bound = math.isqrt((4*n)//431)
        for b in range(-bound, bound+1):
            square = 4*n-431*b*b
            t = math.isqrt(square)
            if t*t != square:
                continue
            for signed_t in {t, -t}:
                if (signed_t-b) % 2 == 0:
                    a = (signed_t-b)//2
                    assert a*a+a*b+108*b*b == n
                    solutions.add((a, b))
        return sorted(solutions)

    norm_two = norm_solutions(2)
    norm_128 = norm_solutions(128)
    assert norm_two == [] and (4, 1) in norm_128
    alpha_a, alpha_b = 4, 1
    alpha_norm = alpha_a**2 + alpha_a*alpha_b + 108*alpha_b**2
    alpha_trace = 2*alpha_a+alpha_b
    assert (alpha_norm, alpha_trace) == (128, 9)
    # w^2-w+108 has the distinct roots 0 and 1 modulo 2.
    roots_mod_two = [x for x in range(2) if (x*x-x+108) % 2 == 0]
    assert roots_mod_two == [0, 1]
    residues = [(alpha_a+alpha_b*x) % 2 for x in roots_mod_two]
    assert residues == [0, 1]
    prime_class_order = 7  # q^7=(w+4), q nonprincipal, and 7 is prime.
    assert 2**prime_class_order == alpha_norm
    steinitz_index = len(forms)//prime_class_order
    assert len(forms) % prime_class_order == 0 and steinitz_index == 3
    return {
        "p": 431,
        "method": "Primitive reduced positive binary quadratic forms and an exact bounded norm-equation search; standard Python.",
        "fundamental_discriminant": discriminant,
        "reduced_form_convention": "|b|<=a<=c, b>=0 when |b|=a or a=c, gcd(a,b,c)=1",
        "reduced_primitive_forms": forms,
        "class_number": len(forms),
        "maximal_order_generator": "w=(1+sqrt(-431))/2",
        "minimal_polynomial_w": "w^2-w+108",
        "norm_formula": "N(a+b*w)=a^2+a*b+108*b^2",
        "norm_search_bound": "431*b^2<=4*n, from (2*a+b)^2+431*b^2=4*n",
        "norm_two_solutions": norm_two,
        "norm_128_solutions": norm_128,
        "split_primes_above_two": ["q=(2,w)", "qbar=(2,w-1)"],
        "alpha": "w+4",
        "alpha_norm": alpha_norm,
        "alpha_trace": alpha_trace,
        "alpha_residues_at_q_qbar": residues,
        "principal_ideal_factorization": "(alpha)=q^7",
        "factorization_reason": "Norm(alpha)=2^7 and alpha is a unit modulo qbar, so q is its only prime-ideal divisor.",
        "q_nonprincipal_reason": "An integral generator of q would have norm 2, and the exact norm equation has no solutions.",
        "order_of_q_in_class_group": prime_class_order,
        "order_of_q_squared_in_class_group": prime_class_order,
        "Steinitz_coset_count": steinitz_index,
        "unoriented_ell_two_component_lower_bound": steinitz_index,
        "exact_component_count_by_dyadic_seven_theorem": steinitz_index,
        "scope": "Verifies the arithmetic inputs for the Steinitz obstruction and the exact dyadic-seven component theorem. The exact count is a theorem evaluation, not a full p=431 PPAS graph enumeration.",
    }


def main():
    records = []
    for p, ell in [(7, 2), (11, 2), (13, 2), (31, 2), (7, 3)]:
        filename = f"p{p}_maximal_graph.json" if ell == 2 else f"p{p}_maximal_ell3_graph.json"
        data = json.loads((RESULTS / filename).read_text())
        labels = [v["label"] for v in data["vertices"]]
        report = check_graph(data["adjacency"], labels, [v["aut"] for v in data["vertices"]])
        stored = sorted(tuple(sorted(labels[i] for i in component)) for component in data["components"])
        computed = sorted(tuple(sorted(component)) for component in report["components"])
        assert stored == computed
        assert report["directed_kernel_count"] == data["directed_kernels"]
        assert report["outdegrees"] == data["outdegrees"]
        types = recover_two_torsion_types(data)
        report.update({"p": p, "ell": ell, "source": filename,
                       "two_torsion_Frobenius_types": types,
                       "nonempty_two_torsion_types": sorted(set(types.values())),
                       "two_torsion_type_counts": {cycle: sum(t == cycle for t in types.values())
                                                   for cycle in sorted(set(types.values()))}})
        if ell == 2:
            degrees_by_type = {"1^6": 15, "2 1^4": 3, "2^2 1^2": 3, "2^3": 7}
            assert report["outdegrees"] == [degrees_by_type[types[label]] for label in labels]
        else:
            for i, label in enumerate(labels):
                for j, target in enumerate(labels):
                    if data["adjacency"][i][j]:
                        assert types[label] == types[target]
            assert all(len({types[label] for label in component}) == 1
                       for component in report["components"])
            assert report["weak_component_count"] == len(set(types.values())) == 4
            for i, vertex in enumerate(data["vertices"]):
                if vertex["type"] == "restriction":
                    assert all(not data["adjacency"][i][j]
                               for j, target in enumerate(data["vertices"])
                               if target["type"] != "restriction")
            ranks = {"1^6": 4, "2 1^4": 3, "2^2 1^2": 2, "2^3": 2}
            assert data["rational_2_torsion_ranks"] == [ranks[types[label]] for label in labels]
            report["odd_prime_component_bound_attained"] = True
            report["restriction_sector_closed"] = True
        records.append(report)
    output = {
        "method": "Standard Python; directed reachability and exact rational Laplacian rank.",
        "scope": "Existing complete adjacency matrices are reused, not recomputed; no new curve census.",
        "local_symplectic_two_torsion_checks": local_two_torsion_check(),
        "elliptic_witnesses": elliptic_witnesses(),
        "mixed_product_star_checks": mixed_star_product_checks(),
        "class_arithmetic_431": class_arithmetic_431(),
        "unoriented_graphs": records,
    }
    (RESULTS / "connectivity_verification.json").write_text(json.dumps(output, indent=2) + "\n", encoding="utf8")
    for record in records:
        print(f"Unoriented p={record['p']}, ell={record['ell']}: "
              f"components={record['component_sizes']}, Laplacian nullity={record['rational_laplacian_nullity']}")
    for record in output["mixed_product_star_checks"]["generic_combinatorial_checks"]:
        print(f"Mixed star-product combinatorics h={record['h']}: "
              f"J3 components={record['J3_components']}, "
              f"K3,3 components={record['K3_3_components']}")
    print("Actual p=11 mixed-product closed J3 component: S3, S6, S8.")
    print("p=431 arithmetic: h(-431)=21, order of [(2,w)]=7, exact component count from the dyadic-seven theorem=3.")
    print("Connectivity, exact Laplacian ranks, Frobenius types, and arithmetic obstructions verified.")


if __name__ == "__main__":
    main()
