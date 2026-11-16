"""Fast, read-only integrity checks for the committed graph data.

Uses only the Python standard library. This checks stored graph consistency
and explicit certificates; it does not rerun the exhaustive enumeration or
reconstruct all isogenies. Run from any directory with Python 3.10 or later.
"""
from collections import Counter, deque
import csv
import hashlib
import json
from math import comb
from pathlib import Path
import sys
import time

RESULTS = Path(__file__).resolve().parent / "results"
# split, restriction, Jacobian; degree histogram; directed kernels; components
EXPECTED = {
    7: ((3, 1, 2), {3: 4, 7: 1, 15: 1}, 34, [6]),
    11: ((10, 2, 8), {3: 12, 7: 6, 15: 2}, 108, [3, 17]),
    13: ((3, 1, 4), {3: 6, 7: 2}, 32, [1, 3, 4]),
    31: ((21, 3, 36), {3: 42, 7: 9, 15: 9}, 324, [60]),
    47: ((55, 5, 100), {3: 110, 7: 25, 15: 25}, 880, [160]),
    71: ((105, 7, 210), {3: 224, 7: 49, 15: 49}, 1750, [322]),
    127: ((55, 8, 227), {3: 210, 7: 40, 15: 40}, 1510, [290]),
    431: ((903, 29, 3436), {3: 3150, 7: 609, 15: 609}, 22848,
          [1456, 1456, 1456]),
}
NEW_PRIMES = {47: 976, 71: 2276, 127: 4536}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_json(name):
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def read_csv(name):
    with (RESULTS / name).open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sha256(name):
    return hashlib.sha256((RESULTS / name).read_bytes()).hexdigest()


def integer(value):
    return isinstance(value, int) and not isinstance(value, bool)


def load_edges(graph):
    n = len(graph["vertices"])
    edges = {}
    if "adjacency" in graph:
        matrix = graph["adjacency"]
        require(len(matrix) == n, "wrong adjacency row count")
        for i, row in enumerate(matrix):
            require(len(row) == n, "wrong adjacency column count")
            for j, count in enumerate(row):
                require(integer(count) and count >= 0, "invalid adjacency entry")
                if count:
                    edges[i, j] = count
    else:
        for row in graph["edges"]:
            i, j, count = row["source"], row["target"], row["multiplicity"]
            require(integer(i) and integer(j) and 0 <= i < n and 0 <= j < n,
                    "invalid edge endpoint")
            require(integer(count) and count > 0, "invalid kernel multiplicity")
            require((i, j) not in edges, "duplicate sparse adjacency entry")
            edges[i, j] = count
    return edges


def component_sets(neighbors):
    remaining = set(range(len(neighbors)))
    components = []
    while remaining:
        root = min(remaining)
        remaining.remove(root)
        found, pending = {root}, [root]
        while pending:
            for target in neighbors[pending.pop()]:
                if target in remaining:
                    remaining.remove(target)
                    found.add(target)
                    pending.append(target)
        components.append(frozenset(found))
    return set(components)


def verify_graph(p):
    kind = "maximal" if p in (7, 11, 13, 31) else "rational"
    name = f"p{p}_{kind}_graph.json"
    graph = read_json(name)
    vertices = graph["vertices"]
    n = len(vertices)
    require(graph["p"] == p and graph["ell"] == 2, "wrong graph parameters")
    require(graph["weil_polynomial"] == [p*p, 0, 2*p, 0, 1],
            "wrong ascending Weil polynomial")
    require(len({v["label"] for v in vertices}) == n, "nonunique vertex labels")
    for i, vertex in enumerate(vertices):
        require(vertex.get("id", i) == i, "vertex id/order mismatch")
        require(integer(vertex["aut"]) and vertex["aut"] > 0,
                "invalid rational automorphism order")
    counts, histogram, kernels, sizes = EXPECTED[p]
    expected_counts = dict(zip(("split", "restriction", "Jacobian"), counts))
    require(Counter(v["type"] for v in vertices) == expected_counts,
            "vertex type counts differ from the release summary")
    require(n == sum(counts), "wrong total vertex count")
    edges = load_edges(graph)
    degrees = [0] * n
    neighbors = [set() for _ in vertices]
    for (i, j), count in edges.items():
        degrees[i] += count
        neighbors[i].add(j)
        reverse = edges.get((j, i), 0)
        require(reverse > 0, "directed support is not symmetric")
        require(count*vertices[j]["aut"] == reverse*vertices[i]["aut"],
                "weighted polarized duality failed")
    require(Counter(degrees) == histogram, "wrong row-degree histogram")
    require(sum(degrees) == kernels == graph["directed_kernels"],
            "wrong directed kernel total")
    if "outdegrees" in graph:
        require(graph["outdegrees"] == degrees, "stored outdegrees disagree")
    for i, vertex in enumerate(vertices):
        if "degree" in vertex:
            require(vertex["degree"] == degrees[i], "stored vertex degree disagrees")
    components = component_sets(neighbors)
    require(sorted(map(len, components)) == sizes, "wrong computed component sizes")
    stored = graph["components"]
    require(sum(map(len, stored)) == n and
            set(map(frozenset, stored)) == components,
            "stored component partition disagrees with graph traversal")
    if "total_vertices" in graph:
        require(graph["total_vertices"] == n, "stored total vertices disagree")
    if "counts" in graph:
        require(all(graph["counts"].get(kind) == count
                    for kind, count in expected_counts.items()), "stored counts disagree")
        if "total" in graph["counts"]:
            require(graph["counts"]["total"] == n, "stored counts.total disagrees")
    return name, graph, edges, degrees, neighbors


def transform_binary_sextic(coefficients, matrix, p):
    """F(aX+bZ,cX+dZ), without Sage or the graph-generator routines."""
    powers = [[pow(x, i, p) for i in range(7)] for x in matrix]
    a, b, c, d = powers
    result = [0] * 7
    for i, value in enumerate(coefficients):
        for j in range(i+1):
            for k in range(7-i):
                result[j+k] += (value*comb(i, j)*a[j]*b[i-j]
                                *comb(6-i, k)*c[k]*d[6-i-k])
    return [value % p for value in result]


def verify_new_exports(p, name, graph, edges, degrees, neighbors):
    vertices = graph["vertices"]
    n = len(vertices)
    labels = [v["label"] for v in vertices]
    index = {label: i for i, label in enumerate(labels)}
    rows = read_csv(f"p{p}_all_vertices.csv")
    require(len(rows) == n, "vertex CSV omits rows")
    require([int(v["id"]) for v in rows] == list(range(n)), "vertex CSV id mismatch")
    for i, (row, vertex) in enumerate(zip(rows, vertices)):
        require(row["label"] == vertex["label"] and row["type"] == vertex["type"],
                "vertex CSV label/type mismatch")
        require(row["model"] == vertex["model_description"] and bool(row["model"]),
                "vertex CSV model mismatch")
        require(int(row["aut_order"]) == vertex["aut"] and
                int(row["degree"]) == degrees[i] and int(row["component"]) == 0,
                "vertex CSV numeric fields disagree")
    edge_rows = read_csv(f"p{p}_all_edges.csv")
    csv_edges = {}
    for row in edge_rows:
        i, j = int(row["source_id"]), int(row["target_id"])
        require(row["source"] == labels[i] and row["target"] == labels[j],
                "edge CSV label/id mismatch")
        require((i, j) not in csv_edges, "duplicate edge CSV row")
        csv_edges[i, j] = int(row["multiplicity"])
    require(csv_edges == edges, "edge CSV differs from graph")
    with (RESULTS / f"p{p}_adjacency.csv").open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.reader(stream)
        require(next(reader) == ["source/target"] + labels, "matrix CSV header mismatch")
        matrix = list(reader)
        require(len(matrix) == n, "matrix CSV row count mismatch")
        for i, row in enumerate(matrix):
            require(row[0] == labels[i] and len(row) == n+1, "matrix CSV label mismatch")
            require([int(x) for x in row[1:]] == [edges.get((i, j), 0) for j in range(n)],
                    "matrix CSV differs from sparse graph")
    root = index["S0"]
    distance, queue = {root: 0}, deque([root])
    while queue:
        i = queue.popleft()
        for j in neighbors[i]:
            if j not in distance:
                distance[j] = distance[i]+1
                queue.append(j)
    require(len(distance) == n, "S0 cannot reach all vertices")
    tree = read_csv(f"p{p}_spanning_tree.csv")
    require(len(tree) == n-1, "spanning tree does not have N-1 edges")
    parent = {}
    for row in tree:
        i, j = index[row["source"]], index[row["target"]]
        require(j != root and j not in parent, "spanning tree repeats a child or enters root")
        require(edges.get((i, j), 0) == int(row["multiplicity"]) > 0,
                "spanning tree uses a nonexistent edge")
        require(int(row["depth"]) == distance[j] == distance[i]+1,
                "spanning tree depth mismatch")
        parent[j] = i
    require(set(parent) == set(range(n)) - {root}, "tree does not cover every nonroot")
    for i, row in enumerate(rows):
        require(int(row["distance_from_S0"]) == distance[i], "CSV distance mismatch")
        require(row["parent_label"] == ("" if i == root else labels[parent[i]]),
                "CSV tree parent mismatch")
    audit = read_json(f"p{p}_full_enumeration_audit.json")
    require(audit["p"] == p and audit["vertices"] == n and
            audit["directed_kernels"] == sum(degrees), "enumeration audit totals mismatch")
    require(audit["graph_sha256"] == sha256(name), "graph hash mismatch")
    raw_name = f"p{p}_exhaustive_branch_models.csv"
    require(audit["exhaustive_csv_sha256"] == sha256(raw_name), "exhaustive CSV hash mismatch")
    for flag in ("every_exhaustive_model_has_explicit_Fp_isomorphism",
                 "every_graph_Jacobian_reached_by_exhaustive_census",
                 "all_automorphism_orders_rechecked_by_120_triples"):
        require(audit[flag] is True, "enumeration audit does not report success")
    models = read_csv(raw_name)
    certificates = read_json(f"p{p}_exhaustive_isomorphism_certificates.json")
    require(len(models) == len(certificates) == NEW_PRIMES[p] == audit["exhaustive_model_rows"],
            "exhaustive model/certificate count mismatch")
    reached = set()
    for line, (row, certificate) in enumerate(zip(models, certificates), 2):
        require(certificate["csv_line"] == line, "certificate CSV line mismatch")
        target = index[certificate["target_label"]]
        require(vertices[target]["type"] == "Jacobian", "certificate target is not a Jacobian")
        reached.add(target)
        f = [int(row[f"c{i}"]) for i in range(7)]
        while f[-1] == 0:
            f.pop()
        source = certificate["normalized_coefficients"]
        require(f == source, "certificate source differs from exhaustive CSV")
        m = certificate["matrix"]
        require(len(m) == 4 and all(integer(x) and 0 <= x < p for x in m),
                "invalid rational projective matrix")
        require((m[0]*m[3]-m[1]*m[2]) % p != 0, "singular isomorphism matrix")
        y = certificate["y_multiplier"]
        require(integer(y) and 0 < y < p, "invalid y multiplier")
        transformed = transform_binary_sextic(vertices[target]["coefficients"], m, p)
        require(transformed == [(y*y*x) % p for x in source + [0]*(7-len(source))],
                "explicit binary-sextic isomorphism identity failed")
    jacobians = {i for i, v in enumerate(vertices) if v["type"] == "Jacobian"}
    require(reached == jacobians and len(reached) == audit["exhaustive_Jacobian_classes"],
            "exhaustive certificates do not cover exactly all Jacobian vertices")
    point_audit = read_json(f"p{p}_rational_census_independent_check.json")
    require(point_audit["input_sha256"] == sha256(f"p{p}_rational_census.json"),
            "point-count audit input hash mismatch")
    require(point_audit["Jacobian_equations_verified"] == len(jacobians) and
            point_audit["zero_Cartier_Manin_matrices"] == len(jacobians) and
            point_audit["points_fp"] == p+1 and point_audit["points_fp2"] == p*p+1+4*p,
            "point-count audit summary mismatch")
    require(read_json(f"p{p}_gluing_verification.json")["p"] == p, "missing gluing audit")
    laplacian = read_json(f"p{p}_laplacian_check.json")
    require(laplacian == {"p": p, "vertices": n, "rank_Q": n-1, "nullity": 1},
            "stored exact Laplacian certificate summary mismatch")


def main():
    started = time.monotonic()
    for p in EXPECTED:
        name, graph, edges, degrees, neighbors = verify_graph(p)
        if p in NEW_PRIMES:
            verify_new_exports(p, name, graph, edges, degrees, neighbors)
        print(f"PASS p={p}: {len(graph['vertices'])} vertices, "
              f"{sum(degrees)} directed kernels, components={EXPECTED[p][3]}", flush=True)
    print(f"All eight stored graphs passed in {time.monotonic()-started:.2f}s.")
    print("This is stored-data/certificate validation, not a fresh exhaustive census.")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, IndexError, TypeError, StopIteration) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        sys.exit(1)
