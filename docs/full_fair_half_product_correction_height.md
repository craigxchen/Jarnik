# A full fair core forces a large half-product correction

The small-correction extension of
[half_symmetric_star_phase_height_gap.md](half_symmetric_star_phase_height_gap.md)
would start with a half-product identity

```text
(product_A z_i)/(product_B z_i)
 = unit * (Gamma_star/bar(Gamma_star))^k * H/bar(H),
|A|=|B|=k.                                             (1)
```

An actual full fair core does **not** supply `log Norm(H)=o(W)` here.
For eight rows it instead forces

```text
log Norm(H) >= (17/16-o(1)) W,       W=log R^2,          (2)
```

with the original common Gaussian content retained. If that common
content has logarithmic norm `o(W)`, the stronger coefficient is
`136/127`. This is a lower bound on the **reduced** correction, after
all conjugate cancellations and the actual small row factors.

The argument is a valuation obstruction to a proposed perturbative
route. It does not establish a point-count bound or apply Roth with
moving targets. It complements the exact affine-extraction obstruction
in [full_cut_affine_star_extraction_obstruction.md](full_cut_affine_star_extraction_obstruction.md).

## 1. Literal core factorization and the exact bound

Let `z_i` be `2k>=4` distinct Gaussian integers of common modulus `R`,
and let `mathcal C` contain one orientation of each
nonconstant cut class, with the star class oriented as `A`. Suppose
the original rows have the exact factorization

```text
z_i=d_0 epsilon_i K_i product_(S in mathcal C)
               Gamma_S^(1_(i in S)) bar(Gamma_S)^(1_(i notin S)),
epsilon_i in {1,-1,i,-i},       K_i in Z[i]\{0}.         (3)
```

Each `Gamma_S` is conjugate-primitive and supported on odd split
primes. Different classes have disjoint rational prime supports.
All `Gamma_S` are nonunits. Write

```text
w_S=log Norm(Gamma_S),
E_K=sum_i log Norm(K_i),
q(S)=|S intersect A|-|S intersect B|.
```

This is the literal disjoint-core setting supplied by
[central truncation](endpoint_central_truncation.md), whose row factors
are small relative to each retained block. Intrinsic fair cut weights
alone, before such a disjoint-core reduction, are not a substitute for
(3): different nested cuts at one prime can have cancelling imbalances.

The original rows have equal modulus, so the residual quotient `F`
obtained after removing `(Gamma_A/bar(Gamma_A))^k` from the half-product
has norm one in `Q(i)`. Therefore there is a conjugate-primitive Gaussian
integer `H`, unique up to a unit, such that `F=u H/bar(H)`, `u in mu_4`.
Choose the representative with no rational common factor. Its norm is
determined exactly by the signed valuations of `F`.

**Correction-height lemma.** With this normalization,

```text
J:=log Norm(H)
  >= sum_(S != A) |q(S)| w_S - E_K.                    (4)
```

To prove it, fix a rational prime `p` in a nonstar block `Gamma_S`.
Let `rho_p` be its actual oriented Gaussian factor in that block and
let its multiplicity be `m_p`. At this chosen Gaussian prime,

```text
v_(rho_p)(F)=q(S)m_p+c_p,
c_p=sum_(i in A)v_(rho_p)(K_i)
       -sum_(i in B)v_(rho_p)(K_i).                    (5)
```

The star factor has no valuation here, by disjoint support. The
opposite valuation is the negative of (5), since `Norm(F)=1`.
Thus this rational prime contributes exactly
`|q(S)m_p+c_p| log p` to `J`, including any reversal of its orientation
in the reduced numerator. Sum the inequality
`|q m+c|>=|q|m-|c|` over all nonstar core primes. Their correction loss
is at most

```text
sum_p |c_p| log p
 <=sum_i sum_p v_(rho_p)(K_i) log p <= E_K.
```

Other primes can only add to `J`. This proves (4), without assuming
that the correcting factors are coprime to the core, to their own
conjugates, or to each other. The unit in (1) is retained and does not
affect the valuation bound.

## 2. Counting the full-cut imbalance

Put

```text
B=2^(2k-1)-1,
B_k=sum_(S != A) |q(S)|.
```

The last sum runs over unoriented nonconstant cut classes. Complementing
a subset changes the sign of `q`, so the unoriented sum is half the
sum over all subsets. The number of subsets with `q=t-k` equals
`binom(2k,t)`: replace the selected part in `B` by its complement.
The elementary telescoping identity

```text
sum_(t=0)^(2k) |t-k| binom(2k,t) = k binom(2k,k)
```

follows from
`(t-k)binom(2k,t)=k[binom(2k-1,t-1)-binom(2k-1,t)]`
on `t>k`, and symmetry. Empty and full cuts have zero imbalance.
The removed star class contributes `k`. Consequently

```text
B_k=(k/2)[binom(2k,k)-2].                              (6)
```

If every core weight is at least `(1-eta)w`, (4) gives the finite,
one-sided statement

```text
J >= B_k(1-eta)w-E_K.                                 (7)
```

No comparison between `w` and the original radius is implicit in (7).

| Half size k | Rows | B | B_k |
|---:|---:|---:|---:|
| 2 | 4 | 7 | 4 |
| 3 | 6 | 31 | 27 |
| 4 | 8 | 127 | 136 |
| 6 | 12 | 2047 | 2766 |

## 3. Original radius, complete common content, and fairness

Let `d` be the **complete** common Gaussian gcd of the original rows,
and put `D=log Norm(d)`, `W=log R^2`. The displayed factor `d_0` in
(3) need not already be this complete gcd. The core rows themselves
have Gaussian gcd one: at every prime in a nonconstant block, some
row uses each orientation. A common valuation beyond `d_0` is therefore
provided by a correcting factor in a row whose core valuation is zero.
It follows that

```text
0 <= D-log Norm(d_0) <= E_K.
```

All `K_i` have the same norm, because the source norms and core norms
are common. Taking norms in (3) therefore gives

```text
W=D+sum_S w_S+O(E_K).                                 (8)
```

Now keep `k` fixed and assume the genuinely fair regime

```text
w_S=w+o(w) for every S,        E_K=o(w),        w -> infinity.
```

Then `J>=B_k w-o(w)` and `W=D+Bw+o(w)`. The latter identity alone
does not allow replacing `w` by `W/B` when `D` is large.

Actual endpoint pair inequalities bound that common content. Every
pair is separated by exactly `2^(2k-2)` cut classes. At a split prime,
adding the correcting valuations changes the absolute allocation
difference by at most the sum of the two correcting valuations.
Summing over both orientations gives

```text
log P_ij=2^(2k-2)w+o(w),                              (9)
```

with error bounded by the fair core error plus
`log Norm(K_i)+log Norm(K_j)`.
For any fixed endpoint constant `C>0`, the arbitrary-unit chord bound
on distinct original rows gives

```text
log P_ij > W/2+log(2/C^2).
```

Thus

```text
W <= (B+1)w+o(w),       D <= w+o(w),
J/W >= B_k/(B+1)-o(1).                                (10)
```

When `D=o(w)`, including a primitive original normalization, (8)
instead gives the sharper coefficient `B_k/B`. For `k=4`, these
are `136/128=17/16` and `136/127`, respectively. All appearances of
`W,D,w_S,J` are logarithms of **squared moduli**; in particular
`|H|=exp(J/2)` and `R=exp(W/2)`.

The central truncation bounds are stronger than the needed
`E_K=o(w)`. Alternatively, if one initially knows only `E_K=o(W)`,
(9) with its explicit correction error and the endpoint inequality
first imply `W=O(w)`, and hence `E_K=o(w)`. This avoids assuming the
desired radius comparison in advance.

## 4. Consequence for the proposed robustness route

For each fixed `k`, the right side of (10) is positive. The literal
correction in (1) consequently cannot have logarithmic norm `o(W)`
in a full fair extracted core. Even its reduced numerator has positive
radius-scale height. One cannot invoke an estimate for a supposedly
small moving target in this setting.

There is also no fixed-power or bounded-norm condition in the fair-core
factor data. To see the power issue directly, assign one fixed distinct
split prime to each class and let its core multiplicity tend to infinity
through `1 mod k`, choosing the multiplicities so all `w_S=w+O_k(1)`.
Take every `K_i=1`. A singleton nonstar class has `q=+1` or `-1`, so
its prime exponent in the reduced `H` is `1 mod k`. Neither `H` up
to a Gaussian unit nor its ordinary norm is a `k`-th power. These are
literal fair Gaussian core models, with no assertion that their phases
form a short arc. Their role is to delimit what the factorization and
fairness conditions themselves force.

This note makes no claim about arbitrary approximate half-symmetric
tuples outside the full-core setting. Uniform Thue estimates or a valid
moving-target approximation theorem would need their own precise
hypotheses; none is supplied or assumed here.

## Verification scope

The [checker](check_full_fair_half_product_correction_height.py) counts
the imbalance coefficients and pair incidences for `k=2,3,4,6`. For
literal four- and six-row Gaussian examples it verifies the exact
quotient after star removal, reduced orientations, deliberately
overlapping correcting factors, and (4) using integer norm comparisons.
No finite computation is being substituted for the asymptotic proof.
