"""Independent, standard-Python audit of the manuscript's stored examples.

This does not rerun the complete curve census. It recomputes point counts and
Cartier--Manin matrices for every listed maximal curve, checks stored graph
invariants, constructs Lagrangians as projective-line sets (independently of
the original row-echelon implementation), and directly evaluates the p=7
nondegenerate Richelot quotients including ground-field square classes.
"""

from fractions import Fraction
from itertools import combinations, permutations, product
from pathlib import Path
import csv
import json

RESULTS = Path(__file__).resolve().parent / "results"


def legendre(a, p):
    a %= p
    return 0 if not a else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def curve_checks(f, p):
    """Use Fp[u]/(u^2-nu); quadratic character is that of the norm."""
    nu = next(a for a in range(2, p) if legendre(a, p) == -1)
    degree = 6 if f[6] else 5
    n1 = p + (1 if degree == 5 else 1 + legendre(f[6], p))
    n2 = p*p + (1 if degree == 5 else 2)
    for x in range(p):
        z = 0
        for c in reversed(f):
            z = (z*x + c) % p
        n1 += legendre(z, p)
    for a, b in product(range(p), repeat=2):
        u = v = 0
        for c in reversed(f):
            u, v = (u*a + nu*v*b + c) % p, (u*b + v*a) % p
        n2 += legendre(u*u - nu*v*v, p)
    h = [1]
    for _ in range((p - 1) // 2):
        new = [0] * (len(h) + 6)
        for i, a in enumerate(h):
            for j, b in enumerate(f):
                new[i+j] = (new[i+j] + a*b) % p
        h = new
    cm = [h[p-1], h[p-2], h[2*p-1], h[2*p-2]]
    assert (n1, n2) == (p+1, p*p+1+4*p)
    assert cm == [0, 0, 0, 0]
    return {"points_Fp": n1, "points_Fp2": n2, "Cartier_Manin": cm}


def graph_check(filename, check_curves=True):
    data = json.loads((RESULTS / filename).read_text())
    p, ell = data["p"], data["ell"]
    vertices, adjacency = data["vertices"], data["adjacency"]
    n = len(vertices)
    assert n == sum(data["counts"][t] for t in ["split", "restriction", "Jacobian"])
    assert sum(map(sum, adjacency)) == data["directed_kernels"]
    curves = []
    for i, v in enumerate(vertices):
        assert sum(adjacency[i]) == data["outdegrees"][i]
        for j, w in enumerate(vertices):
            assert adjacency[i][j]*w["aut"] == adjacency[j][i]*v["aut"]
        if ell != 2:
            assert sum(adjacency[i]) == ell + 2 + legendre(-p, ell)
        if check_curves and v["type"] == "Jacobian":
            curves.append({"label": v["label"], **curve_checks(v["coefficients"], p)})
    seen, components = set(), []
    for i in range(n):
        if i in seen:
            continue
        pending, part = [i], set()
        while pending:
            j = pending.pop()
            if j in part:
                continue
            part.add(j)
            pending.extend(k for k in range(n) if adjacency[j][k] and k not in part)
        seen |= part
        components.append(sorted(part))
    assert sorted(components) == sorted(data["components"])
    result = {"p": p, "ell": ell, "vertices": n,
              "directed_kernels": sum(map(sum, adjacency)),
              "component_sizes": sorted(map(len, components)),
              "curves_checked": curves}
    if check_curves and p in [7, 11, 13]:
        with (RESULTS / f"p{p}_curves.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        mass = sum((Fraction(1, int(row["aut_order"])) for row in rows), Fraction())
        assert mass == p**3
        result["stored_full_census_mass"] = str(mass)
    print(f"Graph p={p}, ell={ell}: {n} vertices, "
          f"{result['directed_kernels']} kernels; invariants passed.")
    return result


def normalize(v, q):
    inverse = pow(next(a for a in v if a), -1, q)
    return tuple(inverse*a % q for a in v)


def pairing(v, w, q):
    return (v[0]*w[1]-v[1]*w[0]+v[2]*w[3]-v[3]*w[2]) % q


def image_line(matrix, v, q):
    return normalize(tuple(sum(a*b for a, b in zip(row, v)) % q for row in matrix), q)


def block_diagonal(a, b):
    return [a[0]+[0, 0], a[1]+[0, 0], [0, 0]+b[0], [0, 0]+b[1]]


def kernel_families(q):
    """Build planes from pairs of projective points, not echelon matrices."""
    lines = sorted({normalize(v, q) for v in product(range(q), repeat=4) if any(v)})
    planes = set()
    for a, b in combinations(lines, 2):
        if pairing(a, b, q):
            continue
        plane = frozenset([b] + [normalize(tuple((x+t*y) % q for x, y in zip(a, b)), q)
                                 for t in range(q)])
        planes.add(plane)
    assert len(planes) == (q+1)*(q*q+1)
    records = []
    for p in [7, 11, 13, 31]:
        if p == q:
            continue
        identity, unipotent = [[1, 0], [0, 1]], [[1, 1], [0, 1]]
        if q == 2:
            cases = [("II", block_diagonal(identity, identity), (9, 6)),
                     ("IU", block_diagonal(identity, unipotent), (3, 0)),
                     ("UU", block_diagonal(unipotent, unipotent), (1, 2))]
        else:
            frobenius = [[0, -p], [1, 0]]
            expected = (4, q-1) if legendre(-p, q) == 1 else (0, q+1)
            cases = [("split", block_diagonal(frobenius, frobenius), expected)]
        exchanged = [[0, 0, -p, 0], [0, 0, 0, -1], [1, 0, 0, 0], [0, p, 0, 0]]
        expected = ((3, 4) if q == 2 else
                    ((q+1, 2) if legendre(-p, q) == 1 else (q+1, 0)))
        cases.append(("restriction", exchanged, expected))
        for label, matrix, expected in cases:
            counts = [0, 0]
            for plane in planes:
                if all(image_line(matrix, v, q) in plane for v in plane):
                    is_product = any(v[2:] == (0, 0) for v in plane)
                    counts[not is_product] += 1
            assert tuple(counts) == expected, (p, q, label, counts, expected)
            records.append({"p": p, "ell": q, "type": label,
                            "product": counts[0], "nonproduct": counts[1]})
    print(f"Independent projective-line Lagrangians: ell={q}, {len(planes)} geometric kernels.")
    return records


def matchings(letters):
    if not letters:
        yield []
        return
    for b in letters[1:]:
        for rest in matchings([x for x in letters[1:] if x != b]):
            yield [(letters[0], b)] + rest


class QuadraticFieldArithmetic:
    """Small exact Fp2 arithmetic; encode a+b*u as a+p*b, u^2=nu."""
    def __init__(self, p, nu):
        self.p, self.nu = p, nu

    def add(self, x, y):
        p = self.p
        return (x % p+y % p) % p + p*((x//p+y//p) % p)

    def neg(self, x):
        p = self.p
        return (-(x % p) % p) + p*((-(x//p)) % p)

    def sub(self, x, y):
        return self.add(x, self.neg(y))

    def mul(self, x, y):
        p = self.p
        a, b, c, d = x % p, x//p, y % p, y//p
        return (a*c+self.nu*b*d) % p + p*((a*d+b*c) % p)

    def inv(self, x):
        p = self.p
        a, b = x % p, x//p
        n = pow((a*a-self.nu*b*b) % p, -1, p)
        return a*n % p + p*((-b*n) % p)

    def frobenius(self, x):
        return x % self.p + self.p*((-(x//self.p)) % self.p)

    def pmul(self, a, b):
        result = [0]*(len(a)+len(b)-1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                result[i+j] = self.add(result[i+j], self.mul(x, y))
        return result

    def bracket(self, a, b):
        derivative = lambda f: [self.mul(i % self.p, x) for i, x in enumerate(f)][1:]
        left, right = self.pmul(derivative(a), b), self.pmul(a, derivative(b))
        return [self.sub(left[i] if i < len(left) else 0,
                         right[i] if i < len(right) else 0)
                for i in range(max(len(left), len(right)))]

    def evaluate(self, f, x):
        z = 0
        for c in reversed(f):
            z = self.add(self.mul(z, x), c)
        return z

    def determinant(self, matrix):
        z = 0
        for perm in permutations(range(3)):
            odd = sum(perm[i] > perm[j] for i in range(3) for j in range(i+1, 3)) % 2
            term = 1
            for i in range(3):
                term = self.mul(term, matrix[i][perm[i]])
            z = self.sub(z, term) if odd else self.add(z, term)
        return z

    def binary_orbit(self, f):
        p = self.p
        for a, b, c, d in product(range(p), repeat=4):
            if (a*d-b*c) % p == 0 or next(x for x in [a, b, c, d] if x) != 1:
                continue
            left, right = [[1]], [[1]]
            for _ in range(6):
                left.append(self.pmul(left[-1], [b, a]))
                right.append(self.pmul(right[-1], [d, c]))
            g = [0]*7
            for i in range(7):
                for j, x in enumerate(self.pmul(left[i], right[6-i])):
                    g[j] = self.add(g[j], self.mul(f[i], x))
            yield g

    def same_curve(self, g, f):
        g = (g+[0]*7)[:7]
        for h in self.binary_orbit(f):
            i = next(i for i, x in enumerate(h) if x)
            scalar = self.mul(g[i], self.inv(h[i]))
            if scalar and scalar < self.p and legendre(scalar, self.p) == 1:
                if all(g[j] == self.mul(scalar, h[j]) for j in range(7)):
                    return True
        return False


def richelot_p7():
    field = QuadraticFieldArithmetic(7, 3)
    data = json.loads((RESULTS / "p7_maximal_graph.json").read_text())
    records = []
    for index, vertex in enumerate(data["vertices"]):
        if vertex["type"] != "Jacobian":
            continue
        f = vertex["coefficients"]
        roots = [x for x in range(49) if field.evaluate(f, x) == 0]
        if not f[6]:
            roots.append(None)
        assert len(roots) == 6
        rho = [roots.index(field.frobenius(x)) if x is not None else roots.index(None) for x in roots]
        hits, decomposed, stable = {}, 0, 0
        for matching in matchings(list(range(6))):
            if ({frozenset((rho[a], rho[b])) for a, b in matching}
                    != {frozenset(t) for t in matching}):
                continue
            stable += 1
            factors = []
            for pair in matching:
                factor = [1]
                for i in pair:
                    if roots[i] is not None:
                        factor = field.pmul(factor, [field.neg(roots[i]), 1])
                factors.append(factor+[0]*(3-len(factor)))
            factors[0] = [field.mul(f[6] or f[5], x) for x in factors[0]]
            delta = field.determinant([q[::-1] for q in factors])
            if not delta:
                decomposed += 1
                continue
            a, b, c = factors
            quotient = field.pmul(field.pmul(field.bracket(b, c), field.bracket(c, a)),
                                 field.bracket(a, b))
            quotient = [field.neg(field.mul(x, field.inv(delta))) for x in quotient]
            while quotient and quotient[-1] == 0:
                quotient.pop()
            assert all(x < 7 for x in quotient)
            targets = [j for j, w in enumerate(data["vertices"])
                       if w["type"] == "Jacobian" and field.same_curve(quotient, w["coefficients"])]
            assert len(targets) == 1
            hits[targets[0]] = hits.get(targets[0], 0)+1
        for j, w in enumerate(data["vertices"]):
            if w["type"] == "Jacobian":
                assert hits.get(j, 0) == data["adjacency"][index][j]
        assert decomposed == sum(data["adjacency"][index][j] for j, w in enumerate(data["vertices"])
                                 if w["type"] != "Jacobian")
        records.append({"source": vertex["label"], "stable_matchings": stable,
                        "Jacobian_targets": {data['vertices'][j]['label']: n for j, n in hits.items()},
                        "decomposed_kernels": decomposed})
    print("Direct p=7 Richelot quotient and Fp-isomorphism checks passed, retaining twists.")
    return records


def main():
    result = {
        "method": "Standard Python; point counts, polynomial arithmetic and projective-line sets.",
        "limitations": [
            "Reuses the stored vertex census; does not independently prove exhaustiveness.",
            "Direct Richelot identification is restricted to nondegenerate p=7 quotients.",
            "Decomposed p=7 Richelot kernels are counted but their individual targets are not recomputed.",
            "Stored ell=3 graph is checked for invariants, not reconstructed by a general odd-isogeny algorithm."
        ],
        "graphs": [graph_check(f"p{p}_maximal_graph.json") for p in [7, 11, 13, 31]],
        "ell3_graph": graph_check("p7_maximal_ell3_graph.json", check_curves=False),
        "kernel_families": [record for q in [2, 3, 5, 7] for record in kernel_families(q)],
        "direct_Richelot_p7": richelot_p7(),
    }
    target = RESULTS / "manuscript_audit_verification.json"
    target.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(f"All audit checks passed. Wrote {target}")


if __name__ == "__main__":
    main()
