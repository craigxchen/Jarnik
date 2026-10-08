# A bounded Newton scan for degree-two cut orders

This note records an exact bounded check of the simplest possible quadratic
certificates in the eight-point rank-three Plucker coordinates.  It is a
support calculation on the actual phase-torus pullback, so it keeps the row
scaling and Vandermonde factors of the `q_i^2/D_i` Gale presentation.

For a triple `I`, use the polynomial coordinate

```text
R_I(z) = (product_(i in I) z_i^2)
         product_(j<k, j,k not in I) (z_k-z_j).
```

At a cut `S`, set `z_i=epsilon*x_i` for `i in S` and leave the other phases
generic units.  If `a=|I intersect S|`, the direct order of `R_I` is

```text
w_S(I) = 2*a + binomial(|S|-a, 2).
```

The minimum over triples is `m_s=(0,1,3,5)` for minority cut sizes
`s=1,2,3,4`.  Therefore, after primitive local normalization, a nonzero
degree-two pullback `f=F(R)` has generic clean-cut order

```text
nu_S(F) = min_(alpha in supp(f)) sum_(i in S) alpha_i - 2*m_s.
```

The checker expands the five-variable Vandermonde by its 120 permutation
terms and combines integer exponent dictionaries.  It then scans every
signed two- and three-term combination in one representative block of each
nontrivial same-weight type:

* the size-three block has the repeated-label (`r=1`) weight; its one
  three-term conic relation is removed by exact cancellation, leaving the
  expected quotient rank `2`;
* the size-ten block has the disjoint-triple (`r=0`) weight.
  Its exact pullback span has the expected quotient rank `5` (the scan only
  asks about sparse two- and three-term combinations).

The weight histogram is checked independently as `{1:476, 3:280, 10:28}`.
Relabeling carries the two representatives to all blocks of their type, and
size-one blocks are individual monomials.  Every degree-two coordinate
product has total order `106` over the 127 unoriented cuts.

The exact scan gives:

| block size | nonzero signed 2/3-term combinations | total-order values | maximum |
|---:|---:|---|---:|
| 3 | 9 | `90` (3), `106` (6) | `106` |
| 10 | 570 | `74` (240), `82` (225), `86` (45), `106` (60) | `106` |

Thus no nonzero signed two- or three-term combination in these same-weight
blocks has a universal generic clean-cut total above `106`.  The word
“generic” matters: specializing the independent `x_i` can create an extra
zero in a leading coefficient.  Such a specialization is a separate
parameter-dependent cancellation; the support calculation proves that it is
not a universal identity.  A finite-field nonzero evaluation is therefore
an admissible upper-order witness for any proposed specialized cancellation,
while the integer support expansion here already certifies the universal
generic statement.

Run the reproducible check with:

```text
python3 docs/check_eight_point_cut_binomials.py
```
