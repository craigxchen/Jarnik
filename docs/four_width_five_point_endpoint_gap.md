# Five points with total primitive width at most four

Let five distinct Gaussian integers of modulus `R` lie on a circular
arc of length `L`. Remove their complete common Gaussian divisor `d`.
At each varying split prime let `e_p` be the range of its allocation
among the five rows. If

```text
E=sum_p e_p <=4,
```

then

```text
L > sqrt(2|d|) sqrt(R)                 arbitrary row units,
L > 2 sqrt(|d|) sqrt(R)                one literal common unit. (1)
```

Thus normalized arcs cannot shrink arbitrarily in the critical
full-prime-rank case `r=4,E=4,M=5`. Indeed (1) excludes that case on
every `C<=sqrt(2)` endpoint arc with arbitrary units, or every `C<=2`
endpoint arc in one common unit class. The statement holds without
full rank, without an asymptotic assumption, and with a moving common
divisor. The constants in (1) are proved lower bounds; sharpness for
actual five-point configurations is not asserted.

The proof uses a weighted four-cube lemma, stronger here than the
quadratic phase-target argument suggested by a simplex determinant.
Earlier four-cube arguments in
[multiplicative_rectangle_separation.md](multiplicative_rectangle_separation.md)
concerned seven-row independent sets in fixed templates. The present
five-row weighted-distance statement applies to literal source rows.
It remains a conditional theorem on primitive width, not a general
uniform point-count bound.

## 1. A weighted four-cube lemma

Give the four coordinates of `{0,1}^4` nonnegative weights `w_j`, with
total `T=sum_j w_j`. Among any five vertices, some pair has weighted
Hamming distance at most `T/2`.

Repeated vertices prove the conclusion immediately. Otherwise suppose
every pair has distance greater than `T/2`.

* No two vertices can have ordinary Hamming distance one. Such a pair
  would force its coordinate weight above `T/2`; a third vertex agrees
  in that coordinate with one member of the pair, giving distance at
  most the weight of the other coordinates, less than `T/2`.
* No two vertices can be complementary. The distances from any third
  vertex to such a pair add to `T`.

Translate a selected vertex to `0000` by coordinatewise XOR. The four
others have ordinary weights two or three. If `k` of them have weight
three, their missing coordinates are distinct. Every weight-two support
must contain all these missing coordinates: otherwise it is contained
in a weight-three support, producing a forbidden Hamming-one pair.

For `k=0`, the six weight-two supports form three complementary pairs,
so at most three may be selected. For `k=2`, at most one weight-two
support remains; for `k=3`, none remains. These cases provide fewer
than five vertices. The only possibilities are therefore:

```text
k=4: 0000 and all four weight-three vertices;
k=1: 0000, one weight-three vertex missing coordinate j,
     and all three weight-two supports containing j.
```

The second is an XOR translate of the first, using its weight-three
vertex as the new origin. In the first configuration the six pairs
among the weight-three vertices differ in the six coordinate pairs.
Their weighted distances sum to `3T`. They cannot all exceed `T/2`.
This proves the lemma.

## 2. Exact Gaussian radius bookkeeping

Write the source rows, with their units retained, as

```text
z_i=d epsilon_i product_p pi_p^a_i(p) bar(pi_p)^(e_p-a_i(p)),
min_i a_i(p)=0,       max_i a_i(p)=e_p,
N_0=product_p p^e_p,  R=|d| sqrt(N_0),  W_0=log N_0.
```

Replace each allocation coordinate by its `e_p` threshold indicators
`1_(a_i(p)>=l)`, each of weight `log p`. This gives `E<=4` binary
coordinates; pad by zero-weight coordinates when needed. The exact
pair norm is

```text
log P_ij=sum_p |a_i(p)-a_j(p)| log p,
```

which is their weighted Hamming distance. The cube lemma supplies a
pair with

```text
P_ij <= exp(W_0/2)=sqrt(N_0)=R/|d|.                    (2)
```

The unit-sensitive primitive chord identity from
[linear_allocation_affine_rigidity.md](linear_allocation_affine_rigidity.md)
gives, for distinct original rows,

```text
|z_i-z_j|^2 >= 2R^2/P_ij >= 2R|d|,
```

or the stronger lower bound `4R|d|` in one literal common unit class.
The positive arc gap is strictly longer than its chord, proving (1).
Coincident allocation vectors cause no exception: distinct source rows
then differ by a unit, and the same exact chord estimate applies.

Combining this with
[minimal_prime_rank_width_bound.md](minimal_prime_rank_width_bound.md)
sharpens its minimal-support consequence: for every fixed `M>=5`,
`r=M-1`, and an endpoint constant in the ranges of (1), the radius is
bounded. If `M>=6`, the total width is already at least five. If `M=5`,
the alternative `E=4` is excluded by (1), and `E>=5` is covered by that
note's Matveev/finite-grid Roth argument. The resulting radius bound
remains generally ineffective and assumes minimal prime support.

## 3. Audit of the proposed determinant/phase route

For five affinely independent vertices of the four-cube, XOR the base
vertex to zero. The resulting `4 x 4` difference matrix `B` is binary.
Augment it to the `5 x 5` sign matrix with first row all ones and other
rows `(1,1-2B_i1,...,1-2B_i4)`. Its determinant has modulus
`16|det B|`. Hadamard's inequality gives

```text
16 |det B| <= 5^(5/2) <64,
1<=|det B|<=3.                                        (3)
```

The exact phase inversion in the minimal-rank note therefore isolates
each prime argument near a member of

```text
T={pi m/(4h) mod 2pi : m in Z, h=1,2,3}.              (4)
```

All slopes in (4) are rational or quadratic. This gives an elementary
uniform lower bound for a split Gaussian prime of norm `p`:

```text
dist(arg(pi_p),T) >= 1/(4p).                           (5)
```

For the odd multiples of `pi/8`, the nonzero integer forms
`b^2+2ab-a^2` and `b^2-2ab-a^2`, with rotations, equal
`sqrt(2)p` times a sine vanishing at the target, with derivative
at most two in angle. For the nonrational `pi/12` targets, use
`p+4ab`, `p-4ab`, `3b^2-a^2`, or `b^2-3a^2`; their angular derivatives
have magnitude at most `4p`. These forms are nonzero at nonzero
integer points in the relevant irrational direction. Axis and diagonal
targets have the stronger integer-coordinate bound of order `p^(-1/2)`.
These facts prove (5) without Roth constants.

A `3 x 3` binary minor has determinant at most two, by the same
augmented Hadamard argument. Adjugate phase inversion thus gives the
coarse error bound `dist(arg(pi_p),T)<=4 delta`, where `delta=L/R`
is the source angular width. From (5), each of the four primes obeys
`p>=1/(16 delta)`. Their product and the **full** common divisor give

```text
R^2=|d|^2 product_p p >= |d|^2/(16 delta)^4,
L >= (sqrt(|d|)/16) sqrt(R).                           (6)
```

So the proposed quadratic route does work and forbids arbitrarily
shrinking normalized arcs. The direct weighted-cube bound (1) is
stronger and also removes its affine-rank requirement.

## Verification scope

The [checker](check_four_width_five_point_endpoint_gap.py) exhausts the
4,368 five-vertex subsets of the four-cube. It verifies the structural
classification used in Section 1, the exact determinant distribution,
and the inverse phase coefficients. It also checks the weighted pair
conclusion for several weight choices and literal Gaussian products
with independent units and common content. The finite checks supplement
the proofs above. No Lean formalization is claimed.
