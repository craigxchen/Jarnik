# Exact weighted clean-cut filtrations in degree two

This note gives an exact support formula for the clean-cut order of every
homogeneous polynomial in the primitive eight-point Gale Plucker
coordinates.  In degree two it converts the proposed weighted certificate
into a finite Newton-support question.  A sparse exact computation settles
that question inside every single row-weight space: no such component has
total clean-cut order greater than `106`.  Thus any quadratic with a gain
must mix row weights.

The last restriction is material.  The calculation below does not prove the
bound for an arbitrary quadratic, and it does not produce a family with no
common zero on the positive actual-circle locus.

## 1. The actual circle pullback

Use complex phase coordinates `z_1,...,z_8` for the eight marked circle
points.  The five moment rows have exponents `-2,-1,0,1,2`.  Complementary
minor duality therefore gives, after one common nonzero factor, the triple
Gale coordinate

```text
R_I=(product_(i in I) z_i^2)
    product_(j<k, j,k not in I)(z_k-z_j),       |I|=3.       (1)
```

This is the complementary Vandermonde formula for the actual moment kernel,
not a substitution into an arbitrary row-scaled conic.  In particular it
retains `sum_i g_i=0`.  It is also equivalent to the `q_i^2/D_i` Gale rows:
both are Plucker presentations of the same kernel of the five circle moment
rows.  The determinant-square and positivity conditions select a real
arithmetic locus inside this presentation; they do not change a generic
formal clean-cut order.

Every `R_I` is homogeneous of total degree `16`, has exponent at most four
in each variable, and satisfies

```text
R_I(1/z_1,...,1/z_8)
 = R_I(z_1,...,z_8)/product_i z_i^4.             (2)
```

The sign in (2) is positive because the complementary Vandermonde has ten
factors.  Hence if `F` is homogeneous of Plucker degree `d` and

```text
f(z)=F(R_I(z)),
```

then the support of `f` is invariant under

```text
alpha -> 4d*(1,...,1)-alpha.                      (3)
```

If `F` is a nonzero class in the degree-`d` coordinate ring, the relevant
case is `f!=0`.  Formula (1), or the earlier boundary-spanning theorem,
checks this injectivity in degree two.  A polynomial which is zero modulo
the Plucker and `sum g_i=0` relations is not a certificate.

## 2. Exact filtration at one clean cut

Represent each cut once by a subset `S` of size `s<=4`; for `s=4`, retain
the representative containing label `1`.  There are

```text
8+28+56+35=127                                      (4)
```

such cuts.  At its generic formal point put

```text
z_i=epsilon*x_i       (i in S),
z_j=u_j               (j outside S),                (5)
```

where all displayed parameters and all differences which should be units
are units.  If `a=|I intersect S|`, (1) gives exactly

```text
ord_S(R_I)=2a+binomial(s-a,2).                      (6)
```

The common primitive-coordinate baselines, obtained by minimizing (6) over
the triples, are

```text
m_1=0,  m_2=1,  m_3=3,  m_4=5.                     (7)
```

Let `Supp(f)` denote the support after all coefficient cancellations have
been performed.  The universal, or generic, clean-cut order is therefore

```text
nu_S(F)=min_(alpha in Supp(f)) alpha(S)-d*m_s.       (8)
```

Thus (8) is an exact filtration, not the minimum of the orders of the
Plucker monomials appearing in a chosen expression for `F`.  It includes
precisely the cancellations relevant to the weighted-order question.
Special boundary parameters may make the order larger, but do not increase
the universally valid order.

Center the reciprocal support by putting

```text
beta=alpha-2d*(1,...,1),
h_S(f)=max_(alpha in Supp(f)) beta(S).               (9)
```

Reciprocity and `sum_i beta_i=0` make `h_S=h_(S^c)`.  Equations (8)-(9)
give

```text
nu_S(F)=2ds-d*m_s-h_S(f).                           (10)
```

Summing the constant term over (4) yields

```text
sum_S (2s-m_s)=373,
sum_S nu_S(F)=373d-W(f),
W(f)=sum_S h_S(f).                                  (11)
```

For `d=2`, this is

```text
sum_S nu_S(F)=746-W(f).                             (12)
```

Directly from (6), every single quadratic Plucker monomial has total order
`106`, so its Newton support has `W=640`.  Consequently the full quadratic
question is exactly

```text
F(R) nonzero  ==>  W(F(R)) >= 640 ?                 (13)
```

A strict weighted gain is the same as a nonzero quadratic pullback with
`W<640`.

## 3. Complete calculation in one row-weight space

Use the 35 coordinates `p_I`, `I subset {1,...,7}`, for `Gr(3,H)`.  Its
degree-two standard monomials split into the familiar one-, two-, and
five-dimensional spaces according as the two triples meet in two or three,
one, or zero labels.  Although the actual phase pullback does not preserve
independent row scaling, it makes sense to ask first whether cancellation
inside one of these spaces can improve (12).

The answer is no.

```text
For every nonzero rational quadratic F contained in one row-weight space,
W(F(R)) >= 640, and hence sum_S nu_S(F) <= 106.       (14)
```

The one-dimensional case is the coordinate-monomial calculation above.
For the repeated-label case, use

```text
m_0=R_123 R_145,   m_1=R_124 R_135,
```

with the third matching represented by `m_0-m_1` up to sign.  The complete
filtration intersection lattice over `F_1009` consists of the whole plane
and three lines.  Its generic width is `656`; all three special lines have
width `640`.

For disjoint triples use the five standard monomials

```text
R_123 R_456, R_124 R_356, R_125 R_346,
R_134 R_256, R_135 R_246.                           (15)
```

Expand (15) using (1).  For each of the 127 cuts, group exponent rows by
their value of `beta(S)`.  On a coefficient subspace `U`, the generic
support height is the largest layer whose coefficient rows do not annihilate
`U`.  The locus where that height drops is exactly the intersection of `U`
with the common kernel of those top-layer rows.  Recursing through these
proper intersections enumerates every possible specialization: any
nongeneric vector must lie in at least one of the children, and dimension
drops strictly at every recursive step.

The resulting finite lattice has 86 subspaces:

| dimension | number of subspaces |
|---:|---:|
| 5 | 1 |
| 3 | 15 |
| 2 | 45 |
| 1 | 25 |

Its minimum width is `640`.  Exactly ten of its finite-field lines attain
that value, namely the ten complementary-triple Plucker monomials expressed
in the basis (15).  The generic width in the five-space is `672`.

This modular enumeration is a characteristic-zero lower-bound certificate.
The five pullbacks in (15) remain independent modulo `1009`.  Given a
nonzero rational coefficient vector, clear denominators and divide its
integer content.  Its reduction modulo `1009` is nonzero, and reduction can
only delete support.  Therefore

```text
W(reduction mod 1009) <= W(characteristic-zero polynomial).
```

Since every nonzero modular vector has width at least `640`, so does every
nonzero rational vector.  Relabeling proves (14) for all spaces of the same
type.

The accompanying checker
[`check_eight_point_weighted_cut_sections.py`](check_eight_point_weighted_cut_sections.py)
expands the sparse alternants exactly, verifies reciprocity, checks the
order `106` for all `1596` coordinate products, verifies independence after
reduction, and performs both complete filtration-lattice enumerations.  It
uses no sampled configurations and no dense `1596`-column rank search.

## 4. The remaining case is now settled by a general theorem

The separate bounds (14) cannot simply be added or minimized across row
weights, because their actual phase supports overlap. The subsequent
[support-multiplicity theorem](polarized_coefficient_support_multiplicity.md)
handles those cancellations directly. It proves `W>=320d` for every
nonzero degree-`d` pullback, and hence (13) for all quadratics as well.

Thus the strict generic weighted-order certificate route is closed in
every degree. The finite calculations above remain independent checks
of its sharp degree-two cases. Actual arithmetic cancellation and
numerical height restrictions are separate questions; no uniform
lattice-circle bound follows from this obstruction.
