"""Draw every rational kernel of the maximal (2,2) graph at p=13.

The existing exact census is the sole data source. A labelled directed
arrow represents its indicated number of distinct kernels, including loops.
Run with Python, matplotlib and networkx installed.
"""
from pathlib import Path
import json
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.path import Path as DrawPath
import networkx as nx


BASE = Path(__file__).resolve().parent
COLORS = {"split": "#246b96", "restriction": "#965511", "Jacobian": "#2e6949"}
EDGE = "#566773"


def main():
    data = json.loads((BASE / "results/p13_maximal_graph.json").read_text())
    vertices, matrix = data["vertices"], data["adjacency"]
    size = len(vertices)
    assert (data["p"], data["ell"], size) == (13, 2, 8)
    assert [sum(row) for row in matrix] == data["outdegrees"]
    assert sum(map(sum, matrix)) == data["directed_kernels"] == 32
    graph = nx.Graph()
    graph.add_nodes_from(range(size))
    graph.add_edges_from((i, j) for i in range(size) for j in range(i + 1, size)
                         if matrix[i][j] or matrix[j][i])
    components = {frozenset(c) for c in nx.connected_components(graph)}
    assert components == {frozenset(c) for c in data["components"]}
    assert sorted(map(len, components)) == [1, 3, 4]
    for i in range(size):
        for j in range(size):
            assert matrix[i][j] * vertices[j]["aut"] == matrix[j][i] * vertices[i]["aut"]

    # Each component has its own planar layout. The middle component is
    # two triangles sharing one edge, so no distinct edges need to cross.
    position = {
        "S0": (0.0, 1.05), "S2": (-1.10, -0.85), "J0": (1.10, -0.85),
        "S1": (4.25, 1.30), "R5": (3.05, 0.0),
        "J3": (5.45, 0.0), "J2": (4.25, -1.30),
        "J1": (8.20, 0.0),
    }
    loop_angle = {"S0": math.pi / 2, "S2": math.pi, "J0": 0,
                  "S1": math.pi / 2, "R5": math.pi,
                  "J3": 0, "J2": -math.pi / 2, "J1": math.pi / 2}
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42,
                         "font.size": 13, "mathtext.fontset": "dejavusans"})
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    radius = 0.27
    rendered = []

    # Draw both orientations separately. The off-centre labels lie on their
    # own curved arrow, including when the two multiplicities differ.
    for i, j in graph.edges:
        for source, target in ((i, j), (j, i)):
            count = matrix[source][target]
            assert count > 0
            xa, ya = position[vertices[source]["label"]]
            xb, yb = position[vertices[target]["label"]]
            dx, dy = xb - xa, yb - ya
            length = math.hypot(dx, dy)
            rad = 0.19
            # Start and finish exactly at the node boundaries, following
            # the tangent of the corresponding quadratic curve.
            txa, tya = dx + 2 * rad * dy, dy - 2 * rad * dx
            txb, tyb = dx - 2 * rad * dy, dy + 2 * rad * dx
            na, nb = math.hypot(txa, tya), math.hypot(txb, tyb)
            start = (xa + radius * txa / na, ya + radius * tya / na)
            end = (xb - radius * txb / nb, yb - radius * tyb / nb)
            arrow = FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=13,
                                    shrinkA=0, shrinkB=0, linewidth=1.25,
                                    connectionstyle=f"arc3,rad={rad}", color=EDGE,
                                    zorder=1)
            ax.add_patch(arrow)
            offset = rad * length / 2 + 0.025
            ax.text((xa + xb) / 2 + offset * dy / length,
                    (ya + yb) / 2 - offset * dx / length, str(count),
                    ha="center", va="center", fontsize=15, color="#27333b",
                    bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.8},
                    zorder=3)
            rendered.append((source, target, count))

    # Cubic loops point away from their component and have a single
    # numerical label. This also displays the isolated vertex's 3 kernels.
    for i, vertex in enumerate(vertices):
        count = matrix[i][i]
        if not count:
            continue
        x, y = position[vertex["label"]]
        angle = loop_angle[vertex["label"]]
        ux, uy = math.cos(angle), math.sin(angle)
        vx, vy = -uy, ux
        def point(outward, tangent):
            return (x + outward * ux + tangent * vx,
                    y + outward * uy + tangent * vy)
        path = DrawPath([point(0.20, -0.18), point(1.00, -0.64),
                         point(1.00, 0.64), point(0.20, 0.18)],
                        [DrawPath.MOVETO, DrawPath.CURVE4,
                         DrawPath.CURVE4, DrawPath.CURVE4])
        ax.add_patch(FancyArrowPatch(path=path, arrowstyle="-|>", mutation_scale=13,
                                    linewidth=1.25, color=EDGE, zorder=1))
        ax.text(*point(0.97, 0), str(count), ha="center", va="center", fontsize=15,
                color="#27333b", zorder=3)
        rendered.append((i, i, count))

    assert {(i, j, matrix[i][j]) for i in range(size) for j in range(size)
            if matrix[i][j]} == set(rendered)
    assert sum(count for _, _, count in rendered) == 32

    for vertex in vertices:
        x, y = position[vertex["label"]]
        ax.add_patch(Circle((x, y), radius, facecolor=COLORS[vertex["type"]],
                            edgecolor="white", linewidth=1.0, zorder=4))
        ax.text(x, y, vertex["label"], ha="center", va="center", color="white",
                fontsize=15.2, fontweight="normal", zorder=5)

    ax.set_xlim(-2.35, 9.20)
    ax.set_ylim(-2.65, 2.85)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor=c,
                               markersize=10, label=kind)
                        for kind, c in COLORS.items()],
               loc="lower center", ncol=3, frameon=False, fontsize=15.2,
               bbox_to_anchor=(0.5, 0.035))
    fig.subplots_adjust(left=0.015, right=0.985, bottom=0.12, top=0.96)
    output = BASE.parent / "figures/maximal_p13_full"
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output.with_suffix(".png"), bbox_inches="tight", dpi=210)
    plt.close(fig)
    print("Verified 8 vertices, 32 directed kernels and component sizes 4, 3, 1.")
    print(output.with_suffix(".pdf"))


if __name__ == "__main__":
    main()
