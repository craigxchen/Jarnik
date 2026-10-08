# Six-factor equal sums: a complete code classification and an ordered cut exclusion

The later [ramified content refinement](sixfactor_moving_content_c_half_exclusion.md)
strengthens the final corollary: every moving rational family covered by
this classification is excluded at `C<=1/2`. This remains a theorem about
the specified exact moment model, with no transfer from arbitrary arcs.
The separate [repeated-magnitude classification](sixfactor_repeated_magnitude_classification.md)
now extends that small-arc exclusion to repeated coefficients by
classifying their aggregate allocation rows.

This note explains the absence of the extreme-pair source cut in the
earlier finite six-factor search. It proves that absence for **all six
real coefficients with distinct nonzero absolute values**, under the
exact equal-first-moment hypothesis. The finite code classification is
computer assisted and exhaustive; the subsequent polynomial sign and
divisibility arguments are explicit. No reduction of arbitrary lattice
circle arcs to this six-factor class is asserted.

The analogous real-moment exclusion already fails with ten factors:
the [certified ten-factor witness](tenfactor_real_outercut_witness.md)
has equal first and third moments and an extreme-pair source column
in fifth-moment order. It does not have a proved rational or Gaussian
integer realization. Thus extending this argument requires additional
information beyond those real moment equations and order.

## Statement on subset sums

Let `a_1,...,a_6` be nonzero real numbers with pairwise distinct absolute
values. Suppose five distinct subsets `S_0,...,S_4` have the same sum
`sum_(j in S_i) a_j`. Put `c_i=sum_(j in S_i) a_j^3`. Then:

1. Up to row permutations, column permutations and independent column
   complements, the five subsets have the single code displayed below.
2. Their five cubic sums `c_i` are distinct.
3. In increasing cubic-sum order, no column has pattern `01110` or
   `10001`: no coefficient occurs in both extreme rows and none of the
   three middle rows, or vice versa.
4. No sixth subset has the same linear sum. In particular every
   subset-sum fiber of six such coefficients has cardinality at most five.

The representative code, with columns numbered one through six, is

```text
S_0={4,6},       S_1={1,3,6},       S_2={1,4,5},
S_3={2,3,5},     S_4={1,2,3,4}.                       (1)
```

This is the code in the existing
[five-point cubic family](five_point_affine_shape_cubic_family.md).
The new classification shows that its order calculation covers every
five-row equal-sum code under the stated coefficient hypotheses.

## Exhaustive finite classification

Write a subset as a six-bit row. Complement columns so that a chosen
row becomes zero, changing the sign of each corresponding coefficient.
If `b` is the chosen reference row, the transformations are

```text
x'_j=x_j xor b_j,       a'_j=(-1)^b_j a_j,
sum_j a'_j x'_j = sum_j a_j x_j - sum_j a_j b_j.        (2)
```

Thus the other four rows have zero weighted sum. No two rows can have
Hamming distance one or two: their equal sums would force a coefficient
to vanish or two coefficients to have equal absolute values.

There are exactly 42 nonzero six-bit masks of weight at least three.
The [exact checker](check_sixfactor_code_classification.py) enumerates
all `binom(42,4)=111930` choices of four masks. Exactly 3150 have all
pairwise distances at least three. For each choice it computes the
rational kernel of the four-row matrix by exact fraction row reduction.
It discards a kernel if it is contained in any of the 36 hyperplanes

```text
a_i=0,                 a_i+a_j=0,      a_i-a_j=0.
```

This is an exact admissibility test over the reals. If the kernel is
contained in none of those hyperplanes, each intersection is a proper
linear subspace; finitely many proper subspaces cannot cover a real
vector space. Conversely a containing hyperplane forbids the required
coefficient vector. The rational kernel basis spans the real kernel,
so checking the defining linear forms on that basis suffices.

Exactly 1800 choices remain, all of rank four. To classify them under
the full symmetry group, the checker tries all 120 row permutations.
Each column is encoded by a five-bit mask, identified with its complement;
sorting the six resulting masks then quotients by column permutations.
Taking the least tuple gives a complete canonical form. Every retained
choice, and (1), has canonical column masks `(3,5,6,9,10,13)`.

This finite exhaustive step proves the unique-orbit assertion. It is
not an inference from coefficients in a bounded numerical range, and
it is not claimed as a Lean formalization.

## Exact order within the unique orbit

The four equal-sum equations for (1) give, directly,

```text
(a_1,...,a_6)=(u,2u,v,u+v,2u+v,3u+v).                (3)
```

For example the first comparison gives `a_4=a_1+a_3`, the second
`a_6=a_1+a_5`, and the last two then give `a_2=2a_1` and
`a_5=2a_1+a_3`. All five subset sums equal `4u+2v`.
Since `u!=0`, put `r=v/u`. Dividing cubic sums by `u^3` and removing
the common term `2r^3` leaves

```text
c_0: 28+30r+12r^2,       c_1: 28+27r+9r^2,
c_2: 10+15r+9r^2,        c_3: 16+12r+6r^2,
c_4: 10+3r+3r^2.                                      (4)
```

If `u<0`, this reverses order, which does not affect the unordered
extreme pair or the prohibited column patterns.

Every real root of every difference in (4) lies in

```text
{-4,-3,-2,-3/2,-1,0,1}.                               (5)
```

At each of these values, a coefficient in (3) is zero or two have equal
absolute values. They are therefore excluded. Between consecutive roots,
all pairwise signs are constant. Factoring the ten differences gives
the following complete increasing-order table:

| Range of r | Row order |
| --- | --- |
| r<-4 | 4,3,1,2,0 |
| -4<r<-3 | 4,1,3,2,0 |
| -3<r<-2 | 1,4,3,0,2 |
| -2<r<-3/2 | 1,2,0,3,4 |
| -3/2<r<-1 | 2,1,0,3,4 |
| -1<r<0 | 2,4,3,0,1 |
| 0<r<1 | 4,2,3,1,0 |
| r>1 | 4,3,2,1,0 |

The possible unordered extreme pairs are consequently
`{0,4},{1,2},{1,4},{2,4}`. Each column of (1) splits the rows two
against three; its two-row side is one of

```text
{0,3},{3,4},{0,2},{1,3},{2,3},{0,1}.
```

The two lists are disjoint. This proves the ordered cut exclusion.

Column complements preserve this argument, not just the linear sums:
the same transformation (2) gives
`sum_j (a'_j)^3 x'_j=sum_j a_j^3 x_j-sum_j a_j^3 b_j`.
Thus it changes all cubic sums by the same additive constant.

For completeness, every one of the other 59 bit rows either never has
sum `4u+2v`, or requires

```text
r in {-5,-4,-3,-5/2,-2,-3/2,-1,-1/2,0,1,2}.            (6)
```

This follows by writing its sum as `Au+Bv`: if `B=2`, only the five
rows in (1) have `A=4`; otherwise equality requires `r=(4-A)/(B-2)`.
Every value in (6) again gives a zero coefficient or equal absolute
values in (3). This proves that a sixth row is impossible. The checker
verifies all 64 equations exactly, along with the difference roots and
all eight sign intervals.

## Consequence for actual Gaussian factor families

Fix integer coefficients satisfying the hypotheses and define

```text
w_i(T)=(-1)^|S_i| product_(j in S_i)(a_j+iT)
                    product_(j notin S_i)(a_j-iT).
```

These are actual equal-norm Gaussian integers for integer `T`. For fixed
coefficients and positive `T` tending to infinity, choose argument lifts
near their common limiting direction. Their differences satisfy

```text
arg w_i(T)-arg w_k(T)
 = -2(sum_(S_i) a_j-sum_(S_k) a_j)/T
   +(2/3)(c_i-c_k)/T^3+O_a(T^(-5)).                    (7)
```

The first term vanishes and all cubic differences are nonzero. Their
eventual actual arc order is the cubic order. Their radius is
asymptotic to `T^6`, and their angular span is `O_a(T^(-3))`, so they
form a fixed endpoint family. Dividing by their common Gaussian gcd
preserves the order and decreases the normalized arc constant.
That gcd also has bounded norm for fixed coefficients: every column
uses both orientations, so a common prime contribution must come from
within-block conjugate collisions or from two distinct blocks. Their
valuations are bounded by the fixed nonzero integers `2a_j` and
`a_j^2-a_k^2`, respectively. Thus the primitive radius still tends to
infinity at order `T^6`, up to coefficient-dependent constants.
The threshold for this order statement depends on the fixed
coefficients; no uniform threshold near the exceptional ratios is used.

Let `K` be the actual forced extreme-pair factor from the
[pair--triple denominator formula](ordered_pair_triple_gap_denominator.md).
It remains bounded in these families even when norm blocks share primes.
More precisely, for every sufficiently large integer `T`,

```text
n_j=T^2+a_j^2,
K | (product_j n_j)/lcm_j n_j
  | product_(j<k) gcd(n_j,n_k)
  | product_(j<k) |a_j^2-a_k^2|.                       (8)
```

To prove the first divisibility, fix an odd split prime. Before common
Gaussian division, its row valuations have the form
`t_i=b+sum_j d_j x_ij`, where `|d_j|<=e_j=v_p(n_j)`.
Choose a column `k` maximizing `e_k`. The one-column allocation
`d_k x_ik` has no extreme-pair gap, by the order theorem, even if
`d_k` is negative. The remaining allocation has oscillation at most
`sum_(j!=k)|d_j|<=sum_j e_j-max_j e_j`. Adding a perturbation can
increase either oriented interval gap by at most its oscillation.
Thus the actual gap multiplicity is at most
`sum_j e_j-max_j e_j`, precisely the prime valuation of the first
quotient in (8). Common Gaussian division shifts all row valuations
equally and leaves the gap unchanged.

For the second divisibility, sort the `e_j`: the sum of all pairwise
minima is at least their sum minus the largest. The last follows from
`gcd(T^2+a_j^2,T^2+a_k^2)` dividing `|a_j^2-a_k^2|`.
The last product is nonzero by hypothesis. Pairwise coprime actual norm
blocks give the sharper conclusion `K=1`.

The denominator bound (8) depends on the coefficients. It is not uniform
when they vary with `T`, and it does not control accidental factors of
the larger index `gcd(Q,N)`.

## A coefficient-uniform small-arc exclusion for the displayed family

A different consequence is uniform even when the coefficients vary.
For any rational `a_1,...,a_6,T` satisfying the coefficient and equal-sum
hypotheses, take the displayed signed products `w_i(T)`, clear a common
Gaussian rational denominator, and divide by their Gaussian gcd. If the
resulting five points are distinct and their primitive radius is `R_0`,
every arc containing them has

```text
arc length / sqrt(R_0) >= 1/sqrt(10).                  (9)
```

This uses the existing
[all-parameter primitive-content bound](six_factor_primitive_content_bound.md)
after the complete code classification. The transfer preserves the actual
points exactly. Complementing column `j` and replacing `a_j` by `-a_j`
multiplies its factor by minus one in every row; it also changes the
row prefactor `(-1)^|S_i|` by minus one. These cancel. Column permutations
do not change the products, and row permutations only relabel them.

Consequently every such family is literally the representative (1) with
coefficients (3). For an integral parameter representative with
`gcd(u,v,T)=1`, that earlier theorem proves
`|G|<=40|P|`, where `P=product_j a_j` and `G` is the full Gaussian gcd.
Its exact chord identity gives normalized span at least
`2 sqrt(|P|/|G|)>=1/sqrt(10)`. Its rational-clearing argument covers the
stated rational parameters as well. No large-`T` limit or cubic-order
threshold is used in (9).

This corollary concerns the specified row prefactors and exact linear
equalities. It does not apply to arbitrary independent row units or to
approximate first moments. The general circle problem has no proved
reduction to this class. In particular (9) is not a proof of a uniform
bound for arbitrary endpoint arcs, and it gives no coefficient-uniform
denominator estimate at larger arc constants.
