# Coupled Gaussian products and the Roth packing threshold

The fixed-target argument in
[the Hadamard Roth note](fixed_hadamard_roth_phase_obstruction.md)
extends from isolated blocks to integral products in the rational span of
the row differences. A fractional packing of total weight greater than
four excludes unbounded endpoint radii. Shared Gaussian primes admit an
exact cancellation refinement. Neither statement supplies a uniform
endpoint bound: the complete balanced-cut model has packing optimum
`2(M-1)/(M-2)`, strictly below four.

## 1. Fixed profile and integral phase certificates

Fix an `m x k` sign matrix `S` and `C>0`. Use nonzero Gaussian integer
blocks `kappa_j`, each coprime to its conjugate, a nonzero Gaussian integer
`d`, and arbitrary row units `u_i`:

```text
z_i = d u_i product_j kappa_j^((1+S_ij)/2)
                         bar(kappa_j)^((1-S_ij)/2),
R = |d| product_j |kappa_j|,
W = sum_j log N(kappa_j) <= 2 log R.                     (1)
```

Assume the points lie in an arc of length at most `C sqrt(R)`. Put
`V=span_Q{S_i-S_0 : i!=0}`. Fix nonzero integer vectors
`v_1,...,v_q in V intersect Z^k`. For each certificate choose once and
for all a rational row vector `lambda_t` such that

```text
lambda_t 1=0,                    lambda_t S=v_t.         (2)
```

Writing `theta_j=arg(kappa_j)`, the lifted row equations yield

```text
v_t dot theta = 2 pi lambda_t l - (pi/2)lambda_t b
               + lambda_t e,
l,b integral,                   |e_i|<=C/(2 sqrt(R)).    (3)
```

Therefore the argument of

```text
A_t = product_j kappa_j^((v_tj)_+)
                bar(kappa_j)^((-v_tj)_+)                (4)
```

is within `O_(S,v_t,C)(R^(-1/2))` of a fixed finite set of rational
multiples of `pi`. Membership in the span of row **differences** matters;
the row span alone need not eliminate the unknown angular center.

## 2. Exact conjugate cancellation, including shared primes

Fix one Gaussian prime `pi_p` above each rational split prime `p` that
occurs in the blocks, with `N(pi_p)=p`. Conjugate primitivity of each
block gives an exact factorization

```text
kappa_j = eta_j product_p pi_p^((e_jp)_+)
                          bar(pi_p)^((-e_jp)_+),
eta_j in {1,-1,i,-i},             e_jp integral,
log N(kappa_j)=sum_p |e_jp| log p.                       (5)
```

No inert prime or `1+i` occurs: either would divide a block and its
conjugate. Different blocks may share primes in either orientation.
Define

```text
E_p(v) = sum_j v_j e_jp,
g_v = product_p pi_p^((E_p(v))_+)
                  bar(pi_p)^((-E_p(v))_+),
P_p(v)=sum_j (v_j e_jp)_+,
Q_p(v)=sum_j (-v_j e_jp)_+.
```

Then (4), with `v=v_t`, satisfies the exact identity

```text
A_v = eta_v [product_p p^min(P_p(v),Q_p(v))] g_v,
eta_v = product_(v_j>=0) eta_j^v_j
        product_(v_j<0) bar(eta_j)^(-v_j).               (6)
```

The bracket is an ordinary positive rational integer. It is a Gaussian
gcd of `A_v` and its conjugate, up to a unit. Thus cancellation changes
the argument only by the recorded unit `eta_v`; no diagonal adjustment
is needed under these hypotheses. Enlarging the target set by the four
unit rotations preserves its finiteness.

The reduced product `g_v` is conjugate-coprime. It is a nonunit exactly
when some `E_p(v)` is nonzero. Its precise logarithmic norm is

```text
h(v) = log N(g_v) = sum_p |E_p(v)| log p
     = sum_j |v_j| log N(kappa_j) - D(v),
D(v)=sum_p [sum_j |v_j e_jp|-|sum_j v_j e_jp|] log p
    =2 sum_p min(P_p(v),Q_p(v)) log p >=0.               (7)
```

This is valid for nested layers sharing the same prime. The assertion
`v!=0` does not imply `g_v` is nonunit in that setting. For example,
`kappa_1=kappa_2` and `v=(1,-1)` give `g_v=1`. Such a certificate must
be omitted from a Roth lower-bound count.

When the blocks have pairwise disjoint rational prime support and are
all nonunits, every nonzero `v` gives `g_v` nonunit and `D(v)=0`.

## 3. Roth lower bound and the fractional packing theorem

For every fixed `epsilon>0`, the phase lemma in
[the preceding note](fixed_hadamard_roth_phase_obstruction.md) gives

```text
h(v_t) >= log R/(2+epsilon) - K_t                       (8)
```

whenever `g_(v_t)` is a nonunit, where `K_t` is independent of the blocks,
their prime supports, and `R`. Indeed the phase error is
`O(R^(-1/2))`, while Roth bounds the error below by a constant times
`|g_(v_t)|^(-2-epsilon)`. The rational-slope and vertical charts, and
the impossibility of an exact root-of-unity direction for a
conjugate-coprime nonunit, are all covered by that lemma.

This invokes the fixed irrational-target theorem of K. F. Roth,
[Rational approximations to algebraic numbers, Mathematika 2 (1955)](https://www.cambridge.org/core/journals/mathematika/article/rational-approximations-to-algebraic-numbers/EFB89E2873019B246F004FAA06E05A7F),
printed page 2. Its constants are generally ineffective.

Choose fixed nonnegative real weights `w_t` satisfying

```text
sum_t w_t |v_tj| <= 1       for every block j,
P = sum_t w_t.                                            (9)
```

Assume every certificate with positive weight has nonunit reduced
product. Multiplying the norm bounds, equivalently summing (7)--(8),
gives

```text
[P/(2+epsilon)] log R - sum_t w_t K_t
 <= sum_t w_t h(v_t)
 <= W <= 2 log R.                                       (10)
```

**Packing theorem.** If `P>4`, then the model (1) has bounded radius
under these hypotheses. Choose `epsilon>0` with `P>4+2 epsilon` in
(10). The threshold depends on the fixed matrix, certificates, weights,
and `C`.

This generalizes the isolated-block theorem: its five certificates are
five standard coordinate vectors with weight one. It is a genuinely
broader hypothesis, but not automatically a stronger result than the
existing elementary product bounds for a given profile.

## 4. What actual cancellation would have to supply

Write `c_j=sum_t w_t|v_tj|` and `D_total=sum_t w_t D(v_t)`. The exact
upper budget is

```text
sum_t w_t h(v_t)
 = sum_j c_j log N(kappa_j) - D_total.                   (11)
```

Under (9), this is at most `W-D_total`. Hence, even if `P<=4`, a family
would be excluded if it had a fixed saving, for some `0<=delta<1`,

```text
D_total >= delta W,
P > 4(1-delta),                                         (12)
```

and all the counted reduced products were nonunits. Choose small positive
`epsilon` so `P>(4+2 epsilon)(1-delta)`; then (8), (11), and `W<=2 log R`
give a strict exponent contradiction. The same proof works if unused
coordinate capacity supplies some or all of this saving.

An alternative exact formulation is a primewise budget: if

```text
sum_t w_t |sum_j v_tj e_jp| <= b sum_j |e_jp|
```

holds at every occurring prime, then `sum_t w_t h(v_t)<=bW`; a fixed
gap `P>4b` suffices. This permits an arithmetic gain from cancellations
that coordinatewise absolute values miss. It does not prove that such
a gain occurs, and cancellation that makes a whole certificate a unit
removes its lower bound as well.

## 5. Fixed power divisibility amplifies the arithmetic bound

There is a separate arithmetic strengthening. Suppose a fixed positive
integer `ell_t` divides every `E_p(v_t)` for a certificate. Then
`g_(v_t)=gamma_t^ell_t` literally with the canonical prime choices in
(5), where `gamma_t` is a conjugate-coprime Gaussian integer. If the
reduced product is nonunit, so is `gamma_t`. Dividing the phase relation
by `ell_t` introduces only the finitely many choices of an `ell_t`-th
root of each target. The same Roth phase lemma now gives

```text
h(v_t) >= [ell_t/(2+epsilon)] log R - K'_t.              (12a)
```

The packing theorem therefore remains valid with its threshold replaced
by `sum_t w_t ell_t>4`. A primewise upper budget `bW` instead requires
`sum_t w_t ell_t>4b`. The powers must be fixed, or chosen from a fixed
finite set: the proof makes no uniform claim for root orders tending to
infinity. A unit multiple of an exact power has the same conclusion
after its four possible unit angles are included.

This is a power-divisibility improvement, not a multivariable Roth
theorem. The importance of exact signed Gaussian powers and retained
units is already explained in
[higher_power_phase_separation.md](higher_power_phase_separation.md).
The present version uses fixed algebraic targets and the exponent two
from Roth rather than the elementary cyclotomic-degree exponent. Its
potential additional scope is when the selected reduced products are
powers without all original blocks being powers. If all blocks are
squares, a global root argument is also available and must not be
mistaken for a new consequence of Roth.

## 6. Complete balanced cuts: the exact packing barrier

Let `M=2k>=4`. Index columns by all `B=binom(M,k)` half-subsets `T`
of the row labels and put `S_i(T)=2 1_(i in T)-1`. Every column sums
to zero, so the span of the row differences equals the row span.

Section 3 of [balanced_pair_products.md](balanced_pair_products.md)
classifies the integral vectors in that span and proves

```text
(1/B) sum_T |v_T| >= a_M := (M-2)/(2(M-1))               (13)
```

for every nonzero integral `v`. Equality is attained by the pair products

```text
v_ij = (S_i+S_j)/2,                 i<j.                 (14)
```

For any fractional packing (9), average the coordinate capacities over
the `B` columns and apply (13):

```text
1 >= sum_t w_t [(1/B)sum_T |(v_t)_T|] >= a_M P.
```

This upper bound is attained. At each column exactly `k(k-1)` unordered
pairs lie on one common side, and these are precisely the pairs for
which `|v_ij(T)|=1`. Giving every pair weight `1/[k(k-1)]` saturates
every coordinate. Therefore

```text
P_max = 1/a_M = 2(M-1)/(M-2) <= 3 < 4.                 (15)
```

Thus **all** fixed integral linear-phase product certificates, not just
pair products selected in advance, fail the Roth packing threshold on
the complete balanced model. For disjoint prime support the exact
cancellation saving in (7) is zero. A saving large enough to repair
the optimum would require

```text
delta > 1-P_max/4 = (M-3)/(2(M-2)),                     (16)
```

which approaches one half. No such saving is available in the disjoint
complete profile. This is a barrier to the displayed certificate method,
not an endpoint realization or counterexample.

The earlier rational-ray argument is stronger on its applicable branches:
integer transverse separation replaces (8) by
`h(v)>=log R-O(1)`, lowering the packing threshold from four to two.
Its central-ray restrictions and complementary-pair exceptions are
handled in [the pair-product note](balanced_pair_products.md). Applying
Roth to arbitrary irrational branches gives wider applicability but a
weaker height exponent; it does not extend that rational-ray exclusion
to all balanced branches.

If every block in this disjoint complete profile is a nonunit Gaussian
square up to a unit, all pair certificates (14) have nonunit reduced
products that are squares. The amplified total mass is
`2P_max=4(M-1)/(M-2)>4`. Thus, for fixed `M,C`, this square-block
subclass has bounded radius even on the irrational central branches.
This relies on the additional even signed exponents; it does not repair
the unrestricted packing obstruction. Related elementary exclusions for
bounded squareclass rank and small arc constants are already recorded in
[squareclass_baseline_and_power_groups.md](squareclass_baseline_and_power_groups.md).
There is substantial further elementary overlap: write
`z_i=d u_i w_i^2` and split the units into their two classes modulo
`{1,-1}`, absorbing a minus sign by replacing `w_i` with `i w_i`.
Within each class coherent square roots lie on an integer circle of
radius `sqrt(R/|d|)` and an arc of length at most
`C/(2 sqrt(|d|))`. Three distinct such roots on a bounded-length arc
are excluded at large root radius by the elementary three-point
Jarnik argument. When `|d|` grows without bound the root arc instead
tends to zero, excluding even two distinct integer roots. Consequently
five or more distinct original square-block points already have bounded
radius by this root argument. The complete disjoint nonunit profile
has distinct points, so for its even `M>=6` the square-block conclusion
is already elementary.

## 7. A profile with no isolated coordinate but a successful packing

Partition nine columns into three triples. Take all `3^3=27` sign rows
having exactly one minus sign in each triple. The row-difference span is

```text
V = {v in Q^9 : the sum on each triple is zero}.          (17)
```

Indeed changing the negative position within one triple gives twice
every coordinate difference in that triple. There is no isolated
coordinate `e_j` in `V`.

Use the three vectors `e_a-e_b` in each triple, with weight one half
on each of the resulting nine vectors. Each column has total load one,
and `P=9/2>4`. If all nine associated primitive pair quotients are
nonunits, the packing theorem excludes unbounded endpoint radii. Pairwise
disjoint nonunit blocks guarantee that condition.

This example verifies that coupled certificates can apply without any
isolated coordinates. It is not claimed as a new exclusion beyond the
elementary methods: its isolating coefficients have denominator two,
so the targets are axis or diagonal directions and rational transverse
separation already gives a stronger bound.

The [rank-one kernel note](rank_one_kernel_relative_roth_obstruction.md)
gives a more substantial application to the eleven-by-eleven core of a
normalized order-twelve Hadamard matrix. It tracks direction collisions
explicitly and obtains a stronger effective bound from the particular
quadratic target directions in that core.

## 8. What has not been obtained

The estimates above apply Roth separately to finitely many scalar slopes
and then use an exact height budget. They do not invoke a simultaneous
approximation theorem. Treating the reduced Gaussian slopes as coordinates
does not automatically give a stronger simultaneous bound: their rational
denominators differ, and the rational height after a common denominator
is introduced must be charged. Their multiplicative identities are already
reflected in the exponent and cancellation formulas (7) and (11).

An improvement through a multivariable theorem would need a new precise
height inequality for these coupled rational points, with hypotheses that
allow the moving prime supports and with its exceptional algebraic sets
controlled. No such theorem is proved or assumed here. The concrete new
output is the fixed coupled-product criterion and its exact cancellation
budget; the complete balanced-cut packing obstruction remains intact.

## Verification

Run the standard-library checker:

```text
python3 docs/check_coupled_block_roth_packing.py
```

The [checker](check_coupled_block_roth_packing.py) verifies 1,491 exact
Gaussian cancellation identities, including 74 unit reductions, and nine
fixed power-root identities. Shared primes, opposite orientations, unit
factors, and a mixed cancellation retaining a nonunit product are included.
It checks the balanced packing primal and dual formulas at `M=4,6,8,10`,
with respective capacities `3,5/2,7/3,9/4`. The general integral-vector
classification used by the dual remains the prose input from
`balanced_pair_products.md`; neither that infinite classification nor
Roth's theorem is being claimed as a numerical test.
