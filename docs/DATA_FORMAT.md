# Graph data and verification scope

All complete surface graphs in this release use the Weil polynomial
`(T^2+p)^2`. Vertices are actual principally polarized **Fp-isomorphism
classes**, with rational twists retained. An edge counts a Frobenius-stable
maximal-isotropic subgroup, not an isogeny modulo automorphisms.

## Authoritative graph files

Use `experiments/results/p{p}_maximal_graph.json` at `p=7,11,13,31` and
`experiments/results/p{p}_rational_graph.json` at `p=47,71,127,431`.
These are all `(2,2)` graphs. The separately named
`p7_maximal_ell3_graph.json` is a `(3,3)` graph and is not part of the
eight-graph release check.

| p | Split | Restriction | Jacobian | Total vertices | Directed kernel edges | Component sizes |
|---:|---:|---:|---:|---:|---:|---|
| 7 | 3 | 1 | 2 | 6 | 34 | 6 |
| 11 | 10 | 2 | 8 | 20 | 108 | 3, 17 |
| 13 | 3 | 1 | 4 | 8 | 32 | 1, 3, 4 |
| 31 | 21 | 3 | 36 | 60 | 324 | 60 |
| 47 | 55 | 5 | 100 | 160 | 880 | 160 |
| 71 | 105 | 7 | 210 | 322 | 1750 | 322 |
| 127 | 55 | 8 | 227 | 290 | 1510 | 290 |
| 431 | 903 | 29 | 3436 | 4368 | 22848 | 1456, 1456, 1456 |

Vertex indices are zero-based positions in the `vertices` array. When
present, `id` equals that position. Labels such as `S0`, `R6`, and `J12`
are unique **within one graph file**; do not join different files by label
without checking representative models. In particular, regenerated rational
censuses can use different representatives/orderings from older maximal
censuses at the same prime.

Common top-level fields:

- `p`, `ell`: base-field characteristic and kernel prime.
- `weil_polynomial`: coefficients in **ascending** order; for this block,
  `[p*p, 0, 2*p, 0, 1]`.
- `vertices`: all representatives and their rational automorphism orders.
- `adjacency` (older format): full integer matrix `B`, with source rows
  and target columns.
- `edges` (newer format): sparse entries `source`, `target`, `multiplicity`.
  Missing entries are zero; each ordered pair occurs at most once.
- `directed_kernels`: sum of all entries, including diagonal entries.
- `components`: partition of all vertex indices into connected components.

Diagonal entries count self-loop kernels. Parallel kernels have been
aggregated into integer multiplicities. In general `B[i][j] != B[j][i]`.
If `a_i` is the vertex's `aut` value, the checked identity is
`B[i][j] * a_j == B[j][i] * a_i`. The support is symmetric, even when
the weighted matrix is not. The outdegree is the entire row sum, including
loops. Component calculations use this support; repeated edges and loops
do not change components.

## Vertex representatives

Every `aut` value counts automorphisms **over Fp preserving the specified
principal polarization**, not all geometric automorphisms.

- `type: "split"`: `elliptic_ids` indexes the file's `elliptic_classes`
  or `elliptic_models`. An elliptic record with `a,b` means
  `y^2=x^3+a*x+b`. The surface has the product principal polarization.
- `type: "Jacobian"`: `coefficients=[c0,c1,...]` means the smooth
  projective model of `y^2=sum(c_i*x^i)`, with the canonical Jacobian
  polarization. Degree-five models have a branch point at infinity.
- `type: "restriction"`: the surface is the Weil restriction of the
  **maximal** elliptic form over Fp2 with the indicated `j` or `j_orbit`,
  with its descended product polarization. A singleton orbit means a
  rational j-invariant; a two-element orbit is a conjugate pair.
  It is not a choice of an arbitrary elliptic twist with that j-invariant.

In rational census files, `field.nu` specifies
`Fp2 = Fp[u]/(u^2-nu)`. An encoded integer `z=a+p*b`, with
`0 <= a,b < p`, denotes `a+b*u`. `null` in a branch-root list means
infinity. `elliptic_model_fp2` gives encoded short-Weierstrass coefficients
where supplied. The human-readable `model_description` repeats the model
in ordinary text.

Some auxiliary elliptic files instead give coefficient vectors relative
to their own `field_modulus`. Convert field bases before comparing these
with the `a+p*b` encoding; equal integer/vector representations in different
bases do not identify equal field elements.

## Complete exports for p=47,71,127

All CSV files have a header; readers should accept UTF-8 with or without
a byte-order mark. These exports are derived from the corresponding
`p{p}_rational_graph.json`:

- `p{p}_all_vertices.csv`: one row per vertex, with `id,label,type,model,
  aut_order,degree,component,parent_label,distance_from_S0`.
- `p{p}_all_edges.csv`: one row per nonzero ordered adjacency entry, with
  `source_id,source,target_id,target,multiplicity`.
- `p{p}_adjacency.csv`: full labelled matrix; its first column is the
  source label and remaining headers are target labels.
- `p{p}_spanning_tree.csv`: `N-1` actual directed edges, with
  `source,target,multiplicity,depth`. Following the parent relation leads
  to `S0`; depth records shortest-path distance from `S0` in the graph.
- `p{p}_vertex_catalogue.html`: human-readable model and adjacency index.
- `p{p}_complete_support.svg` / `.png`: all vertices and non-loop
  undirected support; these drawings omit multiplicities and loops.

`p{p}_exhaustive_branch_models.csv` is a separate enumeration output,
before deduplication. It includes `rational_branch_points,c0,...,c6`.
Repeated isomorphic models are intentional. There are respectively
976, 2276, and 4536 retained rows at 47, 71, and 127, reducing to
100, 210, and 227 rational Jacobian classes.

`p{p}_exhaustive_isomorphism_certificates.json` has one certificate per
CSV data row. `csv_line` is a one-based physical line number including
the header. `target_label` identifies a graph Jacobian, and
`normalized_coefficients` gives the source polynomial. For
`matrix=[a,b,c,d]` and `y_multiplier=e`, the certificate asserts

```text
F_target(a*X+b*Z, c*X+d*Z) = e^2 * F_source(X,Z),
```

as binary sextics over Fp, with nonzero determinant and `e != 0`.
This explicitly retains the square-class condition needed to distinguish
rational twists. The release checker verifies each polynomial identity
using a separate integer/binomial implementation.

## Audit reports and limitations

- `*_full_enumeration_audit.json` records the model-to-vertex coverage
  checks, graph summary and SHA-256 hashes of the graph and exhaustive CSV.
- `*_rational_census_independent_check.json` records separately computed
  Fp/Fp2 point counts, Cartier-Manin matrices and stable-matching checks,
  together with the input census hash.
- `*_gluing_verification.json` records split/restriction boundary checks.
- `*_laplacian_check.json` records a Sage computation of the rational rank
  of `D-B`, where `D` has the row sums on its diagonal.

Run `python experiments/check_release.py` for a fast read-only check with
the Python standard library. It recomputes row sums, weighted duality and
components for all eight graphs. For the three new primes it also checks
every CSV entry, spanning-tree coverage, audit input hashes and every
stored explicit curve-isomorphism certificate. It verifies the *recorded*
Laplacian and point-count summaries; it does not recompute those expensive
algebraic checks or every quotient isogeny. Syntax compilation is a separate
CI step. Passing CI is not a fresh mathematical proof of census completeness.

Global coverage at p=47,71,127 is established by the standalone exhaustive
branch-configuration enumeration, independently of the class-group
connectivity criterion. To repeat that work, rerun the enumerator and audits
using the reproduction instructions, rather than only this release check.
The older p=7,11,13,31 censuses also have direct model-enumeration inputs.

The p=431 graph was obtained as the Richelot closure of all split-product
vertices. Its three computed components are directly verifiable from the
stored matrix, but the claim that this closure contains every vertex in the
global graph uses the theorem that every component contains a split product.
It must therefore not be described as an entirely independent experimental
verification of that same global component theorem.

Some reverse boundary multiplicities are reconstructed using polarized
duality; independent gluing calculations and row-degree saturation provide
cross-checks. Point counts and duality by themselves do not prove that a
census is globally complete or that two rational twists were distinguished.
