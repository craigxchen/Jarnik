# An exact three-row polynomial system with no polynomial extension

This note concerns the explicit construction in
[three_row_quartic_construction.md](three_row_quartic_construction.md).
Its three normalized anchors form a maximal clique in the polynomial
division model, even when a prospective fourth anchor may have arbitrary
degree and arbitrary real coefficients. This is a statement about one
specified polynomial family, not a uniform bound for integer circles.

Write the three real polynomials as

```text
F_1 = (112-152t+121t^2-48t^3+9t^4)/104,
F_2 = (26-421t+8t^2+21t^3-18t^4)/442,
F_3 = (189-944t+162t^2+144t^3-27t^4)/1088.
```

They arise by dividing each real part by the nonzero constant imaginary
part of its Gaussian polynomial numerator. For every pair,

```text
F_i-F_j divides F_i F_j+1 in R[t].                         (1)
```

**Claim.** If a real polynomial `F` satisfies (1) with each of the
three anchors distinct from it, then `F` is one of `F_1,F_2,F_3`.
Thus no fourth distinct normalized polynomial anchor can be added.

## 1. Why this is a finite exhaustive computation

For `F != F_i`,

```text
F-F_i | F F_i+1  iff  F-F_i | F_i^2+1.                    (2)
```

Every difference therefore has degree at most eight. In particular,
every possible `F` has degree at most eight; this degree restriction
is a consequence, not an additional hypothesis.

The roots of `F_i+i`, with multiplicity one, are respectively

```text
R_1 = {i, 2+2i, 8/3-i, 2/3-2i},
R_2 = {i, 2+2i, -8/3-i/3, 11/6-8i/3},
R_3 = {i, 8/3-i, -8/3-i/3, 16/3+i/3}.
```

Each set is disjoint from its conjugate. Consequently `F_i^2+1`
has exactly four distinct irreducible real quadratic factors, up to
scaling:

```text
Q_r(t)=(t-r)(t-bar r),      r in R_i.
```

Its monic real divisors are precisely the sixteen products
`P_(i,S)=product_(r in S) Q_r`, for `S subset R_i`. This includes
the constant direction, the quadratic directions, and the degree-six
and degree-eight directions; no assumption about equal pair degrees
is needed.

Every prospective `F` therefore belongs to one of the sixteen affine
lines

```text
F_i + a P_(i,S),       a in R,                            (3)
```

through each anchor, in the real coefficient space of degree-eight
polynomials. The three anchors are noncollinear. For example, the
two differences from `F_1` are proportional to different monic
quartics, with root pairs `{i,2+2i}` and `{i,8/3-i}` and their
conjugates.

For any point different from the anchors, at least two of its three
anchor lines must be nonparallel: otherwise their common intersection
would make them one line and all anchors collinear. It is therefore
enough to intersect the nonparallel pairs of lines (3), for all three
pairs of anchors. All line data are rational, so every such unique
real intersection has rational coordinates. Irrational coefficients
cannot escape the enumeration.

## 2. The complete intersection certificate

Exactly six distinct coefficient vectors result from the three sets
of `16 x 16` line-pair checks. Besides the three anchors, they are

```text
Q_12 = (-5744+4024t-2277t^2+576t^3-108t^4)/952,
Q_13 = (-13966+11661t-5328t^2+1539t^3-162t^4)/3978,
Q_23 = (37193-4728t-531t^2+1728t^3-324t^4)/19656.
```

Here compatibility with an anchor is interpreted as true when the
candidate is that anchor itself. Direct division gives:

| Candidate | Compatible with `F_1` | With `F_2` | With `F_3` |
|---|---|---|---|
| `F_1` | yes | yes | yes |
| `F_2` | yes | yes | yes |
| `F_3` | yes | yes | yes |
| `Q_12` | yes | yes | no |
| `Q_13` | yes | no | yes |
| `Q_23` | no | yes | yes |

This proves the claim. The accompanying
[exact certificate](check_quartic_extension.py) uses only Python's
standard-library rational arithmetic. It constructs the sixteen
directions, solves each line intersection in all nine coefficients,
asserts that the full candidate set equals the six polynomials above,
and independently tests compatibility by polynomial division. It
does not search over a bounded coefficient box or use numerical
root approximations. Two independent audits checked the exhaustion
argument and exact polynomial arithmetic. This is not a Lean proof.

## 3. What this obstruction does and does not cover

The original three anchors give four actual circle points with
primitive radius comparable to `n^7`, arc length `O(n^3)`, and all
seven intrinsic cut blocks of equal asymptotic weight. The claim
shows that this particular polynomial construction cannot gain a
fifth point by appending another constant-imaginary polynomial
numerator satisfying all three division conditions.

It does not exclude another polynomial family, a modification of the
three anchors, or numerators whose imaginary parts vary with the
parameter. Nor does it exclude isolated integer completions of the
specialized circles: polynomial identities are substantially stronger
than divisibility at selected integer parameters. No transfer from
an arbitrary endpoint sequence to this fixed family has been proved.

The result gives an exact example where all three pairwise candidate
lists have additional intersections, but their simultaneous
intersection contains only the original anchors. Making this
simultaneous obstruction independent of the chosen polynomials is
the further requirement; the requested uniform bound is still open
within this work.

A subsequent [parameter-dependent certificate](quartic_conic_family_maximality.md)
extends this maximality to every nondegenerate real member of the specified
conic family. A [different real family](four_row_quartic_eight_block_family.md)
does have four compatible quartic rows, so the family restriction remains
essential.
