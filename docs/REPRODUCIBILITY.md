# Reproducibility and verification scope

This document describes the computational material accompanying **Principally Polarized Superspecial Surfaces over Prime Fields: Counts and Isogeny Graphs**. Commands below are run from the repository root unless stated otherwise.

## Mathematical scope and conventions

The main vertex set consists of principally polarized abelian surfaces $(A,\lambda)/\mathbb F_p$, for primes $p\ge7$, satisfying

$$
P_A(T)=\det(T-\pi_A\mid V_rA)=(T^2+p)^2,
\qquad r\ne p.
$$

Thus $\pi_A^2=[-p]$, and the surfaces become maximal over $\mathbb F_{p^2}$. The three vertex types are split elliptic products, Weil restrictions, and genus-two Jacobians. This is a fixed Weil block, not the set of all supersingular abelian surfaces or all superspecial rational forms.

Two vertices are identified only through a polarization-preserving isomorphism defined over $\mathbb F_p$. Equality of geometric invariants alone is insufficient. The equation comparison includes the rational square multiplier, so nonisomorphic rational twists remain separate.

For $\ell\ne p$, an outgoing edge is a Frobenius-stable maximal-isotropic subgroup $K\subset A[\ell]$. Its quotient has the induced principal polarization. The kernel is required to be defined over $\mathbb F_p$; its individual points need not all be rational. For $\ell=2$, a Jacobian's geometric kernels are the 15 matchings of its six Weierstrass points.

The adjacency entry $B_{ij}$ is the number of these kernels with quotient isomorphic to vertex $j$. Self-loops and different kernels with the same quotient are retained. Although $B_{ij}>0$ if and only if $B_{ji}>0$, the two multiplicities can differ. With $a_i=|\operatorname{Aut}_{\mathbb F_p}(A_i,\lambda_i)|$, the check is

$$
B_{ij}a_j=B_{ji}a_i.
$$

## Environment

The supported reproduction setup uses:

- Python 3.10 or newer;
- SageMath for elliptic-model enumeration, finite-field/class-group calculations, and exact rational matrix ranks;
- a C++17 compiler such as `g++` for exhaustive model enumeration;
- NumPy for independent point counts and polynomial computations;
- Matplotlib and NetworkX for drawings;
- ReportLab for auxiliary plotting workflows.

Use the repository's [environment.yml](../environment.yml), or install the packages in [requirements.txt](../requirements.txt) into an existing SageMath-compatible environment. SageMath is not supplied by `requirements.txt`. The main reproduction command uses SageMath's Python so that the interpreter running the stages can import `sage.all`.

```sh
conda env create -f environment.yml
conda activate ppas-over-fp
```

On Windows, run the SageMath/C++ stages under WSL. No personal installation paths or external accounts are needed.

## Validate the supplied release data

```sh
python experiments/check_release.py
```

This is the entry point for checking the packaged results without rebuilding the entire census. It verifies the available stored graph data and their declared totals. It is a data validation step: reading representative equations and an adjacency matrix does not independently reconstruct the objects or their quotient maps.

The expected complete rational $(2,2)$-graph results are:

| $p$ | Split | Restriction | Jacobian | Total | Directed kernels | Components |
|---:|---:|---:|---:|---:|---:|---|
| 7 | 3 | 1 | 2 | 6 | 34 | 6 |
| 11 | 10 | 2 | 8 | 20 | 108 | 3, 17 |
| 13 | 3 | 1 | 4 | 8 | 32 | 1, 3, 4 |
| 31 | 21 | 3 | 36 | 60 | 324 | 60 |
| 47 | 55 | 5 | 100 | 160 | 880 | 160 |
| 71 | 105 | 7 | 210 | 322 | 1750 | 322 |
| 127 | 55 | 8 | 227 | 290 | 1510 | 290 |
| 431 | 903 | 29 | 3436 | 4368 | 22848 | 1456, 1456, 1456 |

The last column lists component sizes, not their number.

## Rebuild p = 47, 71, 127

```sh
sage -python experiments/reproduce.py --primes 47 71 127
```

The command regenerates outputs in [experiments/results](../experiments/results). The `--primes` argument accepts any subset of 47, 71, and 127; omitting it runs all three. The execution log and completion summary are written to `build/reproduce.log` and `build/reproduction_summary.json`. It carries out these stages:

1. **Elliptic enumeration.** Enumerate the supersingular elliptic $\mathbb F_p$-forms. Enumerate all geometric supersingular $j$-invariants over $\mathbb F_{p^2}$ and select their unique maximal forms. These supply the complete split and restriction vertex families.
2. **Explicit quotient construction.** Compute split-product and restriction gluing seeds, expand every rational Richelot kernel of each discovered Jacobian, and identify quotients over $\mathbb F_p$, retaining twists.
3. **Graph assembly and boundary checks.** Assemble the full integer adjacency matrix. Compare explicit gluing quotients with reversed degenerate Richelot rows and verify row degrees and weighted duality.
4. **Independent global census.** Run the C++ branch-point enumerator, which has no graph or class-group input. Match every surviving model to a graph vertex through an explicit rational coordinate change and square multiplier, and require every graph Jacobian to occur in that census.
5. **Further verification and exports.** Recheck point counts, Cartier-Manin matrices, and branch matchings; compute connectivity by graph traversal; export a spanning tree; compute the exact rational Laplacian rank; produce catalogues, CSV files, and support drawings.

The word *independent* in stage 4 refers to how the candidate models are enumerated. That enumerator uses neither the graph's reachability nor the class-group connectivity criterion. The subsequent isomorphism comparison shares finite-field/model infrastructure with other scripts and is supplemented by explicit polynomial identities and a separate ordered-triple search; this is not a second entirely unrelated software implementation of every operation.

### Why the branch enumeration covers all required curves

For a curve in this block, Frobenius squared is $[-p]$ on its Jacobian. It acts trivially on two-torsion and hence fixes all six branch points. Those points therefore lie in $\mathbb P^1(\mathbb F_{p^2})$, with one of four rationality patterns. The enumerator covers them using rational projective changes of coordinate:

- **Six rational points:** fix three points at infinity, 0, and 1; choose the remaining three.
- **Four rational points:** fix three at infinity, 0, and 1; choose the fourth and every irreducible quadratic pair.
- **Two rational points:** fix them at infinity and 0. Writing a nonrational root as $a+bu$, with $u^2=\nu$ a fixed nonsquare, scaling puts one chosen pair into the form with real part 1, or into $\{u,-u\}$ when $a=0$. Choose the other pair freely.
- **No rational points:** fix one conjugate pair at $\{u,-u\}$. Enumerate the second pair modulo the stabilizer of that pair, then every possible disjoint third pair.

Both leading square classes are retained. The candidates are squarefree by construction; point-count and Cartier-Manin filters select the stated maximal block. The use of these filters relies on the mathematical fact that this maximal block is superspecial, not on any graph connectivity statement.

| $p$ | Candidate branch configurations | Surviving models, with repeats and twists | Distinct Jacobian classes |
|---:|---:|---:|---:|
| 47 | 113572 | 976 | 100 |
| 71 | 400188 | 2276 | 210 |
| 127 | 2333812 | 4536 | 227 |

The complete graphs have rational Laplacian ranks 159, 321, and 289, respectively. The spanning trees have those same numbers of edges and reach all 160, 322, and 290 vertices. The branch enumerator was also compared with the older exhaustive censuses at 7, 11, 13, and 31; see [branch_census_regressions.json](../experiments/results/branch_census_regressions.json).

### How the edges are obtained

For Jacobians, every Frobenius-stable branch matching is evaluated. With descending coefficient rows in the determinant $d$ and brackets $L=B'C-BC'$, $M=C'A-CA'$, $N=A'B-AB'$, the smooth quotient uses the rational normalization $-LMN/d$. Its sign matters over $\mathbb F_p$: changing it can select a different twist.

Split product kernels use elliptic Velu quotients. Nonproduct kernels use explicit Howe-Leprévost-Poonen two-gluing, with affine-map tests for decomposed outputs. Restriction product kernels use all three geometric elliptic two-kernels. Smooth restriction quotients are checked independently using exchanged-factor descent of two-gluing. Some restriction reverse multiplicities in the assembled matrix are obtained by polarized duality from explicitly computed opposite edges; the boundary and row-saturation checks record this dependency. They are not all computed by a separate restriction-to-surface quotient routine.

## Earlier exhaustive examples and the p = 431 construction

The $p=7,11,13,31$ maximal Jacobian censuses use exhaustive normalized curve enumeration. Their representative equations and adjacency matrices are included, along with verification reports. The older `*_curves.csv` and `*_richelot.json` files at small primes also include other supersingular Weil blocks; do not combine their counts with the maximal-block files. Files named `*_maximal_*` or `*_rational_graph.json` identify the graphs tabulated above. Older shell runners and individual programs remain available for the historical computations; the consolidated reproduction command above is the entry point for 47, 71, and 127.

To rebuild the larger construction from its included elliptic data:

```sh
python experiments/run_431_census.py
```

This enumerates gluing seeds from every split product, closes the discovered Jacobians under all rational Richelot kernels, assembles the graph, and runs the associated boundary, point-count, and arithmetic checks.

**Global completeness at 431 depends on the component theorem.** The closure contains 903 split products, 29 restrictions, and 3436 Jacobians. It is closed under all required kernels, and direct graph traversal gives three components of size 1456. To conclude that no additional component lacking a split vertex exists, the construction uses the theorem that every component contains a split vertex. Thus the graph traversal and arithmetic calculation agree, but this experiment does not independently prove that theorem or independently exclude further components. There is no separate full curve census or independent Hermitian-lattice enumeration at 431 in this release.

The arithmetic result $h(\mathbb Q(\sqrt{-431}))=21$ and the order-seven class of a prime above two predict three components. That prediction is a separate calculation of the theorem's arithmetic inputs; it is not an independent geometric completeness certificate.

## Files and encodings

All paths in this section are relative to [experiments/results](../experiments/results). See [Data format](DATA_FORMAT.md) for graph schemas, field encodings, and the checks performed on each file family.

| Pattern | Contents |
|---|---|
| `p*_rational_graph.json` | Sparse directed adjacency, vertex models, automorphism orders, degree and component statistics |
| `p*_maximal_graph.json` | Earlier complete graphs, generally with a dense adjacency array and labelled vertex models |
| `p*_all_vertices.csv` | One row per vertex, its representative model and a path-tree parent |
| `p*_all_edges.csv` | One row per nonzero ordered adjacency entry, with kernel multiplicity |
| `p*_adjacency.csv` | Complete labelled adjacency matrix for 47, 71, 127 |
| `p*_maximal_ell2_adjacency.csv` | Earlier labelled maximal-block adjacency matrices |
| `p*_spanning_tree.csv` | Actual graph edges witnessing reachability from a fixed root |
| `p*_vertex_catalogue.html` | Human-readable model and outgoing-adjacency tables; open locally after downloading the accompanying SVG |
| `p*_complete_support.svg` / `.png` | Drawings of all vertices and the undirected support; loops and multiplicities are omitted only from the drawings |
| `p*_exhaustive_branch_models.csv` | Surviving exhaustive models before rational-isomorphism reduction |
| `p*_exhaustive_isomorphism_certificates.json` | Coordinate-change and square-multiplier certificates matching those models to graph vertices |
| `p*_full_enumeration_audit.json` | Exhaustive coverage, distinctness, automorphism, and reachability checks |
| `p*_laplacian_check.json` | Exact rank and nullity over $\mathbb Q$ |
| `p*_rational_census_independent_check.json` | Separate point counts, Cartier-Manin matrices, branch matching counts, and stored-edge consistency checks |
| `p*_gluing_verification.json` | Explicit split/restriction gluing comparisons and boundary checks |

Polynomial coefficient arrays are in ascending order: `[c0,c1,...]` represents $c_0+c_1x+\cdots$. In the rational-census field encoding, an integer `a+p*b` represents $a+bu$, where `u^2=nu` and `nu` is recorded in the file. A null branch root denotes infinity. Older Sage-generated extension-field files instead give coefficient vectors and a defining modulus; read their field metadata rather than assuming the two encodings coincide.

Labels such as `S0`, `R...`, and `J0` identify vertices within one dataset. Numeric IDs and representative equations need not agree between separately generated historical datasets, even when the represented graph is the same. Compare labelled models or rational isomorphism classes before comparing matrices from different constructions.

## Additional experiments

The supplied $p=7,\ell=3$ dataset is a complete maximal-block graph. The `*_ten_bit_restriction_ell3.json` files describe complete restriction sectors at selected primes, and the `*_jacobian_neighborhood.json` files describe explicit neighborhoods. Neither kind of file claims a complete surface census at those primes. Verification reports distinguish fresh finite-field calculations from checks reusing saved quotients or automorphism orders.

The numerical data do not prove the all-prime degree formulas, asymptotics, or component classification. Likewise, the arithmetic criterion by itself does not enumerate vertices, produce paths, or determine component sizes. The individual experiments are intended to make those different tasks and dependencies inspectable.

For a reproducibility record, retain the repository commit, software versions, command line, regenerated outputs, and the relevant verification reports. Timings, graph-layout coordinates, and file hashes can vary between environments; the mathematical comparisons are made using exact data and the stated isomorphism conventions.
