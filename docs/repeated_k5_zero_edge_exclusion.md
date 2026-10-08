# A zero edge and opposite disjoint edge weights in the K5 moment model

**Lemma.** Suppose the ten real edge weights of `K5` have equal sums and
equal cube sums at all five vertices. If one edge is zero and two disjoint
edges have opposite nonzero weights, there must be another zero edge.

This is the repeated-magnitude boundary lemma needed for the multilevel
extension of the [nine-factor moment exclusion](five_row_nine_factor_moment_exclusion.md).
The proof uses two explicit integer polynomial identities, verified by a
[standard-library checker](check_repeated_k5_zero_edge.py). No numerical
feasibility test or unverified symbolic elimination is required.

Relabel the zero edge `01`. Divide every weight by one of the nonzero
opposite weights. There are exactly two configurations under vertex
permutations preserving the zero edge, and possible exchange of the
opposite pair: either both endpoints of `01` belong to the opposite
pair's four endpoints, or just one does. Use representatives

```text
A: x01=0, x02=1, x13=-1;
B: x01=0, x12=1, x34=-1.
```

These normalizations preserve both equal-moment hypotheses. In each case
the four vertex-sum differences and the three fixed-edge constraints have
rank seven. Solving them gives the following complete affine
parametrizations, with three arbitrary real parameters `u,v,w`:

| Edge | Case A | Case B |
| --- | --- | --- |
| 01 | 0 | 0 |
| 02 | 1 | v+w-3 |
| 03 | u+2v+w+1 | u-v+2 |
| 04 | u+w-2 | -u+2v+w-1 |
| 12 | u+v+2w-1 | 1 |
| 13 | -1 | -u+2v+2w-3 |
| 14 | u+v+2 | u |
| 23 | u | v |
| 24 | v | w |
| 34 | w | -1 |

For `i=1,2,3,4`, define the cubic integer polynomial

```text
F_i(u,v,w) = (sum_(e incident to i) x_e^3
                    - sum_(e incident to 0) x_e^3) / 3.
```

The division by three is exact coefficientwise in both parametrizations.
The accompanying [certificate data](repeated_k5_zero_edge_certificates.json)
contain four integer polynomials `L_i` of degree at most four for each
case. Direct expansion verifies

```text
Case A: sum_(i=1)^4 L_i F_i = 3200 u.
Case B: sum_(i=1)^4 L_i F_i = 9600.
```

When all cube sums agree, every `F_i` vanishes. Case B is impossible.
Case A forces `u=x23=0`, giving a second zero edge, distinct from `01`.
This proves the lemma. The identities in fact hold over any field of
characteristic zero.

The JSON uses the edge order `01,02,03,04,12,13,14,23,24,34`.
Every polynomial is an array of monomials `[a,b,c,C]`, meaning
`C u^a v^b w^c`. The variable names retained in the JSON are `x7,x8,x9`
for Case A and `x6,x7,x8` for Case B; they refer respectively to `u,v,w`.
The checker reconstructs the vertex sums and cube sums directly from the
edge forms, checks that the linear parametrizations are complete, and
verifies both identities with integer polynomial arithmetic. The
`denominator` field records the common denominator cleared from the
initial rational multipliers; the stored multipliers are already integer.

Run the complete verification with

```sh
python3 docs/check_repeated_k5_zero_edge.py
```

The identities were found by exact rational elimination on monomial
multiples of the four cubics through multiplier degree four. Their
verification does not depend on that search or on a computer-algebra
package.
