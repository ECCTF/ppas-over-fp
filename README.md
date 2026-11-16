# PPAS over prime fields

Computational companion to **Principally Polarized Superspecial Surfaces over Prime Fields: Counts and Isogeny Graphs**.

This repository, `ECCTF/ppas-over-fp`, contains representative equations, rational isogeny graphs, verification reports, and the programs used to construct them. The calculations concern principally polarized abelian surfaces over a prime field with Weil polynomial

$$
P_{A/\mathbb F_p}(T)=(T^2+p)^2,\qquad p\ge 7.
$$

Vertices are **polarized $\mathbb F_p$-isomorphism classes**. Distinct rational twists are retained. An edge is a Frobenius-stable maximal-isotropic subgroup of $A[\ell]$, with its polarized quotient as target. Edge totals count directed kernels, including self-loops and multiplicities; they are not counts of undirected simple edges.

## Quick check

From the repository root, validate the supplied results with:

```sh
python experiments/check_release.py
```

This checks the packaged data using only Python's standard library. It does not replace a fresh curve enumeration or a mathematical proof of the general theorems. For the full workflow and the independence of each check, see [Reproducibility](docs/REPRODUCIBILITY.md); for graph schemas, field encodings, and verification coverage, see [Data format](docs/DATA_FORMAT.md).

## Complete rational (2,2)-graphs

| $p$ | Split products | Weil restrictions | Jacobians | Vertices | Directed kernels | Component sizes |
|---:|---:|---:|---:|---:|---:|---|
| 7 | 3 | 1 | 2 | 6 | 34 | 6 |
| 11 | 10 | 2 | 8 | 20 | 108 | 3, 17 |
| 13 | 3 | 1 | 4 | 8 | 32 | 1, 3, 4 |
| 31 | 21 | 3 | 36 | 60 | 324 | 60 |
| 47 | 55 | 5 | 100 | 160 | 880 | 160 |
| 71 | 105 | 7 | 210 | 322 | 1750 | 322 |
| 127 | 55 | 8 | 227 | 290 | 1510 | 290 |
| 431 | 903 | 29 | 3436 | 4368 | 22848 | 1456, 1456, 1456 |

The completeness arguments differ:

- **7, 11, 13, 31:** the Jacobian classes were exhaustively enumerated independently of the connectivity theorem.
- **47, 71, 127:** a separate exhaustive branch-point enumeration covers all Jacobian classes in the stated Weil block and is matched to the graph by explicit rational isomorphisms. Connectivity is then checked on the actual adjacency matrices.
- **431:** the construction starts from every split-product vertex and closes under rational Richelot quotients. Its global completeness uses the theorem that every component contains a split vertex. Direct traversal finds three components of size 1456 in this closure; it is **not an independent test excluding additional components** without that theorem.

The repository also includes the complete $p=7,\ell=3$ graph, selected complete Weil-restriction sectors, and selected Jacobian neighborhoods. A sector or neighborhood is not a complete surface graph.

## Inspect individual vertices and edges

| $p$ | Vertex catalogue | Vertex CSV | Adjacency CSV | Structured graph |
|---:|---|---|---|---|
| 47 | [HTML](experiments/results/p47_vertex_catalogue.html) | [CSV](experiments/results/p47_all_vertices.csv) | [CSV](experiments/results/p47_adjacency.csv) | [JSON](experiments/results/p47_rational_graph.json) |
| 71 | [HTML](experiments/results/p71_vertex_catalogue.html) | [CSV](experiments/results/p71_all_vertices.csv) | [CSV](experiments/results/p71_adjacency.csv) | [JSON](experiments/results/p71_rational_graph.json) |
| 127 | [HTML](experiments/results/p127_vertex_catalogue.html) | [CSV](experiments/results/p127_all_vertices.csv) | [CSV](experiments/results/p127_adjacency.csv) | [JSON](experiments/results/p127_rational_graph.json) |

Download or clone the repository and open a catalogue locally to view its graph and follow links between neighboring vertices. The diagrams display the undirected support; the tables and machine-readable files retain directions, loops, and kernel multiplicities.

Earlier complete graphs are stored as [p7](experiments/results/p7_maximal_graph.json), [p11](experiments/results/p11_maximal_graph.json), [p13](experiments/results/p13_maximal_graph.json), and [p31](experiments/results/p31_maximal_graph.json). The larger construction is [p431](experiments/results/p431_rational_graph.json). All other records are in [experiments/results](experiments/results).

## Recompute the explicit censuses

Use Python 3.10 or newer, SageMath, a C++17 compiler, NumPy, Matplotlib, and NetworkX. Auxiliary plotting workflows also use ReportLab. The supplied [environment.yml](environment.yml) describes a SageMath environment; [requirements.txt](requirements.txt) lists the Python packages. Installing the Python packages alone does not install SageMath.

```sh
conda env create -f environment.yml
conda activate ppas-over-fp
sage -python experiments/reproduce.py --primes 47 71 127
python experiments/check_release.py
```

On Windows, run the SageMath/C++ workflow inside WSL or another supported Linux environment. The reproduction command regenerates files under `experiments/results`; preserve any local changes before running it.

The workflow constructs the graph, independently enumerates branch models, checks their coverage and rational isomorphisms, verifies the Laplacian ranks, and exports catalogues and drawings. [Reproducibility](docs/REPRODUCIBILITY.md) explains the individual stages, data conventions, and the distinct status of the $p=431$ construction.

## Interpreting the evidence

Finite-field, integer, and rational arithmetic are exact in the mathematical checks. Floating-point arithmetic is used for drawing layouts only. An adjacency-matrix check confirms properties of that matrix; it does not by itself establish that every required vertex or quotient was constructed. Finite numerical examples support the general results but do not prove statements for all primes.

When referring to these computations, cite the accompanying manuscript by its title and record the repository commit used, the primes and degrees checked, and whether results were regenerated or only the stored data were validated.
