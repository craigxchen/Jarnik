# Complete fair cuts force near-critical mixtures onto row pairs

This is a combinatorial obstruction to extending the Walsh three-family
height cancellation by a mixture of balanced partition certificates. It
does not assert that a complete fair sign profile has a short-arc Gaussian
realization. The argument applies to the **pointwise half-height covering**
used by such mixtures: the expected absolute coefficient at every physical
column is at most `1/2`, while the total certificate weight must be close
to one so that the common phase lower bound does not lose a fixed fraction
of the full height.

Let `M=2k>=6`. The complete fair profile has one sign column for every
`k`-subset `T` of `[M]`, with `s_i(T)=+1` on `T` and `-1` outside. (Physical
copies of columns do not affect the statements.) Consider any nonzero
certificate `c in {0,+1,-1}^M` with `sum_i c_i=0`; its support has even
size `m`. Its primitive Gaussian log height at a column of log-prime weight
`w_T` has coefficient

```text
v_c(T)=|sum_i c_i s_i(T)|/2.
```

The mean coefficient over all fair columns depends only on `m`:

```text
A_M(m) = E_(|T|=k) v_c(T).
```

For an arbitrary nonnegative mixture `lambda_c`, write
`L=sum_c lambda_c` and `alpha=sum_(|supp c|>=4) lambda_c`. If the mixture
has the pointwise upper baseline

```text
sum_c lambda_c v_c(T) <= 1/2       for every balanced T,       (1)
```

then

```text
L <= (M-1)/M - [(k-2)/(2k-3)] alpha.                         (2)
```

Thus a mixture with `L>=1-O(1/M)` can assign only `O(1/M)` weight to all
certificates supported on four or more rows. A hierarchical partition,
Haar, or martingale mixture consisting mostly of longer balanced signed
blocks cannot have both (1) and near-unit total weight on the full fair
profile. Walsh cosets evade this calculation because the Walsh column set
is a small, structured subset of all fair cuts; its three families have
disjoint label blocks and exact cancellation coefficients.

## Proof of the mean gap

Select one `+1` and one `-1` entry of `c` and call their contribution
`X=s_i(T)-s_j(T)`. The remaining entries sum to zero. Conditional on
the signs at the selected entries, each remaining row has the same mean
sign, so the remaining contribution `Y` has conditional mean zero.
Jensen's inequality gives `E(|X+Y| | s_i,s_j)>=|X|`. Repeating this
after selecting any number of positive-negative pairs proves

```text
A_M(m) >= A_M(m-2) >= ... >= A_M(2).                         (3)
```

Two distinct rows are separated by a random balanced cut with
probability `2k^2/[M(M-1)]=k/(2k-1)`, so

```text
A_M(2)=k/(2k-1)=M/[2(M-1)].                                (4)
```

For `m=4`, take `c=(1,1,-1,-1,0,...)`. If `r_+` and `r_-` denote the
numbers of plus signs of `s(T)` among the two positive and two negative
certificate entries, then `v_c(T)=|r_+-r_-|`. Counting the cases of
difference one and two gives

```text
A_M(4)
 = 4 [binom(M-4,k-2)+2 binom(M-4,k-3)] / binom(M,k)
 = k(3k-5)/[(2k-1)(2k-3)].                                  (5)
```

Consequently `A_M(4)/A_M(2)=1+(k-2)/(2k-3)`, which is at least `4/3`
for `k>=3`. Averaging (1) over all balanced columns and applying
(3)--(5) yields

```text
1/2 >= A_M(2)L+[A_M(4)-A_M(2)]alpha.
```

Dividing by `A_M(2)` gives (2).

## Universal pair mixture and its flip scale

The uniform mixture of row-pair differences `c=e_i-e_j`, with total
weight `(M-1)/M`, saturates (1) for every balanced column. This gives
the exact general balanced-cut analogue of the Walsh baseline:

```text
[(M-1)/M] E_(unordered i,j) v_(e_i-e_j)(T) = 1/2.
```

Suppose one entry of one balanced physical column is flipped. Its plus
set now has size `k+1` or `k-1`. The row-pair separation probability is
`2(k^2-1)/[M(M-1)]`, so that column's mixed coefficient becomes exactly

```text
1/2 - 2/M^2.                                               (6)
```

If distinct physical columns of total log-prime weight `F_tot` receive
one-entry flips, the universal pair mixture has upper height
`W_0/2-2F_tot/M^2`. Its total phase lower bound has the unfilled mass
slack `W_0/(2M)` (plus the explicit arc constants). The flip gain in (6)
cannot close that slack when `F_tot<=W_0`, because
`2F_tot/M^2<=2W_0/M^2`. On a repeated complete fair profile, the number
of balanced columns is exponentially large while one flip per row touches
only `M` physical columns. Pair averaging therefore does not supply a
Walsh-like large gain from the row allocation.

The obstruction is scoped to nonnegative mixtures of zero-sum
`{0,±1}` row certificates with the pointwise upper coefficient (1).
Other certificates, profile-dependent weighted inequalities, or joint
phase constraints can still use information unavailable to this mixture.

The [exact checker](check_complete_fair_partition_mixture_barrier.py)
enumerates all balanced cuts for even `M=4,...,16`, checks every even
certificate support size, and verifies both directions of the pair-mixture
flip decrement. The proof above covers every even `M>=6`.
