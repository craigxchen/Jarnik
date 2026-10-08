# An ordered conductor-index criterion for a uniform endpoint bound

The uniform circle theorem remains unproved. This note gives a sufficient
arithmetic statement on **five ordered actual lattice points**, using the
conductor part of their affine-circuit index. The statement is not assumed
true. Unlike a sign-profile test alone, its index is calculated from the
actual integer triangle determinants, including their common contents.

The arithmetic input is the proved
[ordered index divisibility theorem](ordered_affine_circuit_conductor_index.md).
The new step here is an ordered sampling argument: adjacent pair separation
forces many cuts to alternate often, which forces a positive average mass
of cuts separating the first and last of five selected rows from the
three middle rows.

## The proposed five-point estimate

For five distinct lattice points `P_1,...,P_5`, in order on an arc, let
`D_ijk` be their signed triangle determinants. Put

```text
g = gcd_(i<j<k) |D_ijk|,
g_L = gcd_(1<=i<j<k<=4) |D_ijk|,
g_R = gcd_(2<=i<j<k<=5) |D_ijk|,
Q = g |D_234|/(g_L g_R).
```

The number `Q` is a positive integer: it is the index of the two primitive
consecutive four-row affine circuits in the saturated affine-relation
lattice. It is invariant under division by a common Gaussian factor.
After that division, let `N_5` be the squared radius of the primitive
five-tuple, and define its conductor-supported index by

```text
J_5 = gcd(Q,N_5).
```

Here is the presently unproved sufficient estimate. For every primitive
five-tuple on an arc of length at most `(1/2)sqrt(R_5)`, suppose

```text
log J_5 <= B,                                           (1)
```

for one absolute `B>=0`. Then every arc of length `C sqrt(R)` has at most

```text
ceil(2C) max(33, 1+8B/log 2)                            (2)
```

lattice points. Nonintegral right sides can be rounded up. A weaker
estimate `log J_5<=delta log N_5+B` also suffices for every fixed
`0<=delta<1/15`; the explicit dependence is given below. This threshold
uses the five-tuple's own primitive norm, including the removal of layers
constant on that tuple.

The bound concerns `gcd(Q,N_5)`, not `Q` itself. The linked index note
exhibits an actual five-point endpoint family with growing `Q` but
`gcd(Q,N_5)=1`. That example does not prove (1) or refute it.

There is now a strictly weaker sufficient target: bound only the forced
cut factor `K`, the reduced denominator of `N_outer N_inner/N_5`, where
the pair and triple are independently Gaussian-primitive. The
[exact denominator formulation](ordered_pair_triple_gap_denominator.md)
proves this dictionary and applies the same `delta<1/15` argument.
Actual and local examples in the
[accidental-conductor classification](ordered_index_accidental_conductor.md)
show why controlling all of `gcd(Q,N_5)` is an additional burden.

## Prime layers and adjacent separation

It suffices first to consider an arc constant `C_0=1/2`. Divide all its
`M` points by their common Gaussian gcd. If its modulus is `|G|>=1`,
the new radius is `R'=R/|G|`, and the allowed arc constant becomes
`C_0/sqrt(|G|)<=C_0`. Thus primitive reduction preserves this class.

Write `N=R^2` for the primitive squared radius and `W=log N`. In a
primitive equal-norm Gaussian tuple, inert and ramified prime factors
are common and have been removed. Each remaining rational prime `p`
splits, and its row allocations have minimum zero and maximum
`e_p=v_p(N)`. Replace each allocation by its threshold bits. A layer
has weight `w=log p`; the sum of all layer weights is `W`.

For a pair of rows put

```text
d_ij = sum_(layers separating i,j) w.
```

Their Gaussian gcd has norm `exp(W-d_ij)`. Dividing their difference
by that gcd leaves a nonzero Gaussian integer, of modulus at least one.
The arc bound consequently gives

```text
exp((W-d_ij)/2) <= |z_i-z_j| <= (1/2)exp(W/4),
d_ij >= W/2+2log 2.                                   (3)
```

This needs no common-unit assumption. If there is more than one point,
`W>0`. Summing (3) over all pairs, while a binary cut with `r` ones
separates at most `M^2/4` pairs, gives

```text
W >= 4(M-1)log 2.                                     (4)
```

Let `T` be the number of transitions of a layer's binary word in the
actual arc order. Summing (3) over adjacent pairs gives

```text
(1/W) sum_layers w T >= (M-1)/2.                      (5)
```

## An exact ordered-pattern bound

For a binary word, count the increasing five-subsets with pattern
`01110` or `10001`. A word with `T` transitions contains an alternating
subsequence of length `T+1`, by taking one entry from every run.
The number of the two patterns in that subsequence is exactly

```text
F(T) = binom(floor((T+4)/2),5)
       + binom(floor((T+3)/2),5),                     (6)
```

where a binomial coefficient is zero when its top entry is less than
five. Thus the original word has at least `F(T)` occurrences. To check
the count, an alternating word of length `2k` has `binom(k+1,5)`
occurrences of each pattern; at length `2k+1` the two counts are
`binom(k+2,5)` and `binom(k+1,5)`. Choosing the three internal entries
and the two enclosing entries gives these binomial counts directly.

Extend `F` linearly between consecutive nonnegative integers. It is
nondecreasing and convex, because Pascal's identity gives

```text
F(T+1)-F(T) = binom(floor((T+3)/2),4).                 (7)
```

Jensen's inequality and (5) show that, if `A_I` is the total layer
weight separating the two extreme rows of a selected ordered five-set
`I` from its three middle rows, then

```text
sum_(|I|=5) A_I >= W F((M-1)/2).                      (8)
```

This already gives a positive asymptotic average, since
`F((M-1)/2)/binom(M,5) -> 1/512`. The full pair inequalities give
the stronger limiting coefficient `1/16`, as follows.

## The full pair inequalities give the fair-pattern limit

Use signs `s_i in {-1,1}` for a layer and the unnormalized weighted
expectation `E f=sum_layers w f`. Thus `E 1=W`. Equation (3) implies

```text
C_ij=E(s_i s_j)=W-2d_ij <= 0,
sum_(i<j) (-C_ij) <= MW/2.                             (9)
```

The second inequality follows from `E(sum_i s_i)^2>=0`. On any
interval of `n` row indices, with sign sum `S`, one likewise has
`E S^2<=nW`. Let `e_4` be the fourth elementary symmetric function
of its signs. Newton's identities give

```text
24 e_4 = S^4-(6n-8)S^2+3n(n-2).
```

For `n>=2`, use `S^4>=S^2` and `E S^2<=nW` for a lower bound.
For `n>=4`, use `S^4<=n^2 S^2` for an upper bound; the cases
`n<=3` have `e_4=0`. Consequently

```text
-n(n-1)W/8 <= E e_4 <= n(n-1)(n-2)W/24.              (10)
```

For a selected ordered five-set, its pattern indicator is

```text
16 I = (1-s_1s_2)(1-s_1s_3)(1-s_1s_4)(1+s_1s_5).
```

After expanding, its ten quadratic coefficients are `+1` or `-1`.
The quartic coefficient is negative when the omitted index is the first
or last of the five, and positive when the omitted index is one of the
three middle indices. Summed over all five-sets, (9) bounds the entire
quadratic contribution from below by
`-(MW/2)binom(M-2,3)`.

For a four-set `a<b<c<d`, its coefficient in the quartic sum is
`2(d-a)-M-2`: adding an index strictly inside its span contributes
positively, and adding an exterior index contributes negatively.
Writing `e_4(full)`, `e_4(prefix r)` and `e_4(suffix M-r)` for the
corresponding sign polynomials, this gives the exact identity

```text
quartic sum = (M-4)e_4(full)
              -2 sum_(r=1)^(M-1) [e_4(prefix r)+e_4(suffix M-r)].
```

Apply (10) and `sum_(n=1)^(M-1)binom(n,3)=binom(M,4)`.
The weighted quartic contribution is at least
`-W[(M-4)M(M-1)/8+binom(M,4)]`. Therefore, for every `M>=5`,

```text
(1/binom(M,5)) sum_(|I|=5) A_I >= beta_M W,
beta_M=(1/16)[1-10/(M-1)-5/(M-4)-15/((M-2)(M-3))].     (11)
```

In particular `beta_M -> 1/16`, and `beta_M>1/32` for `M>=34`.
The latter follows by checking the rational expression at `34`; each
subtracted positive term decreases thereafter. No exchangeability or
independence of the cuts is assumed. The estimate uses their actual
weighted pair correlations and controls the higher correlations directly.

Let `Z_I` be the weight of layers constant on a selected five-set. Its
indicator has the expansion

```text
16 I_constant = 1 + sum_(i<j) s_i s_j
                 + sum_(|A|=4) product_(i in A) s_i.
```

In the sum over all five-sets, every pair occurs `binom(M-2,3)` times
and every four-set occurs `M-4` times. Equations (9)--(10) therefore give

```text
(1/binom(M,5)) sum_(|I|=5) Z_I >= gamma_M W,
gamma_M=(1/16)[1-10/(M-1)-15/((M-2)(M-3))].            (11a)
```

Both `beta_M` and `gamma_M` tend increasingly to `1/16`. Fair pattern
extraction was already proved in
[uniform profile extraction](uniform_profile_extraction.md); sorting
its extracted tuple preserves the assertion about all patterns. The
argument here supplies explicit average estimates for the particular
ordered pattern and for the constant layers.

## Applying the actual index divisibility

For each ordered five-set, remove its own Gaussian gcd. The threshold
layers constant on it disappear; their removal changes neither the
extreme-pair patterns nor the circuit index. At every odd split prime,
the linked index theorem proves that the number `H_p` of remaining
extreme-pair threshold layers satisfies

```text
p^H_p divides gcd(Q_I,N_I).
```

All varying layers arise from odd split primes, so

```text
A_I <= log gcd(Q_I,N_I),       log N_I=W-Z_I.          (12)
```

The norm identity is exact: the layers constant on the selected rows
are precisely those removed by its common Gaussian gcd.

Assume (1). Combining (11)--(12) gives `beta_M W<=B`.
If `M>=34`, then `W<=32B`, and (4) gives
`M<=1+8B/log 2`. If `M<34`, the alternative bound is `33`.
Partitioning an arbitrary allowed arc into `ceil(2C)` subarcs of
length at most `(1/2)sqrt(R)` proves (2). Shared boundary points do
not invalidate the upper bound obtained by summing the subarc counts.

More generally, if the unproved hypothesis is

```text
log gcd(Q_I,N_I) <= delta log N_I+B,     0<=delta<1/15,
```

then averaging `A_I+delta Z_I<=delta W+B` gives
`(beta_M+delta gamma_M-delta)W<=B`. Put
`eta=(1-15delta)/32>0`, and choose an integer `M_delta>=34` such that

```text
beta_(M_delta)+delta gamma_(M_delta)-delta >= eta.
```

Such an integer exists and depends only on `delta`, because the left
side tends to `(1-15delta)/16`. Monotonicity of this displayed lower
bound and (4) yield

```text
M <= max(M_delta-1, 1+B/(4 eta log 2)).                (13)
```

The partition factor for a general arc constant is again `ceil(2C)`.
If all the conductor-supported indices are one, the stronger elementary
observation `F(T)=0 => T<=5` already gives `M<=11` from (5).

## Scope and verification

Neither (1) nor its positive-exponent weakening is proved here. In
particular, the general point-count growth estimate is unchanged.
The available elementary upper bound is only
`J_5<=Q<=|D_234|<=sqrt(R_5)/64=N_5^(1/4)/64` at arc constant
`1/2`, using the circumradius formula and the total arc length.
Its exponent `1/4` is above the sufficient threshold `1/15`.
The new unconditional ingredients are the ordered pattern averaging and
its application to the exact conductor-index divisibility. They reduce
the full goal to a specific five-point arithmetic estimate that retains
the actual common triangle contents and the order of the source points.

The [exact checker](check_ordered_conductor_index_uniformity_criterion.py)
exhausts short binary words, verifies the alternating-pattern count,
discrete convexity and quartic-sum identity, and checks the rational
cutoff, constant-pattern expansion and weighted orthogonal-sign fixtures
in (11)--(11a).
The all-length combinatorial argument and the primewise index theorem
are prose inputs, not conclusions inferred from finite enumeration.
