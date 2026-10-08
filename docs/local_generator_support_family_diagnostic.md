# Support-family diagnostic for localized kernel generators

This note screens support restrictions of the localized determinant
generators against the `(n,d)=(3,4)` target. It records exact support
counts and one-coalition vanishing orders, then separates those facts
from the global integral-lattice statement needed for a height gain.

## Exact support counts

For `m=nd`, let `I` be a `2n`-element support. The localized generator is

```text
L_I(P,Y;v) * H(v),
L_I=det[(P_i v_i^a)_(a=1)^n | (Y_i v_i^a)_(a=1)^n]_(i in I),
```

where `H` is a degree-`d-2` invariant on the complement, expressible as a
linear combination of products of `d-2` rainbow `n`-row determinants.
The complement invariant space has dimension

```text
r_(n,d-2)=(n(d-2))! product_(j=0)^(n-1) j!
             / product_(j=0)^(n-1)(d-2+j)!.
```

Multiplication by nonzero `L_I` is injective, so the fixed-support family
has this exact rank. At `(3,4)`, there are `binom(12,6)=924` supports and
`r_(3,2)=5` directions per support, giving 4620 displayed generators.
The coherent kernel rank is 341.

Fix a coalition `T` of size `r`. If `t=|I intersect T|`, the determinant
factor has exact generic collision order

```text
g(t)=max(0,t-n).
```

The complement factors are disjoint in the vector variables, so they do
not change this common order of the `L_I H` coefficient family. The
number of supports in stratum `t` is exactly

```text
N_t=binom(r,t) binom(m-r,2n-t),
```

with out-of-range binomial coefficients interpreted as zero. For
`(n,d)=(3,4)` and `r=6`, the support counts for `t=0,...,6` are
`1,36,225,400,225,36,1`; each count is multiplied by five for displayed
columns. The corresponding one-coalition orders are
`0,0,0,0,1,2,3` per column. Thus supports contained in `T` (the single
support `I=T`) give five independent directions of order three there;
supports disjoint from `T` have order zero.

For a fixed support `I`, summing `g(|I intersect T|)` over all coalitions
`T subset [m]` gives the exact combinatorial total

```text
2^(m-2n) * sum_(t=n+1)^(2n) binom(2n,t)(t-n).
```

At `n=3,m=12`, this is `64*(15+12+3)=1920` per displayed direction.
This total is a sum of separate collision orders. It is not an asserted
integer divisor of a specialized vector, nor a covolume exponent.

## Exact height and covolume for one fixed support

Write `v_I` for the integer coefficient vector of `L_I`, in the raw
monomial basis. Its entries are signed products of the binary brackets
along perfect matchings of `I`. For any core cut `T`, the minimum number
of matching edges internal to `T intersect I` is `g(|T intersect I|)`.
Each cut has a matching attaining this minimum. Core prime supports for
different cuts are disjoint, so the cut-content divisor

```text
D_I = product_T n_T^g(|T intersect I|),
```

divides the gcd of all matching products. For each core prime, choosing
a matching attaining the minimum at that cut bounds the gcd valuation
above by the exponent in `D_I` plus the residual bracket contribution.
Consequently

```text
D_I divides gcd(matching products),
gcd(matching products) <= D_I * B_I,
log B_I=o(w).
```

under the stated full-profile/core hypotheses. The quotient `B_I` records
the bracket residue and row-correction factors; it need not be bounded
independently of the profile. This gives the exact leading gcd exponent,
without requiring one matching to minimize every cut at once. In the
full profile, every bracket has logarithmic height `2^(m-2)w+o(w)`, so
every matching product has height
`n*2^(m-2)w+o(w)`. Summing the cut exponents gives

```text
log D_I = 2^(m-2n) * sum_(t=n+1)^(2n) binom(2n,t)(t-n) * w + o(w).
```

After division by the gcd, the primitive vector therefore satisfies

```text
log ||v_I|| = B(n,2)*2^(m-2n)w+o(w),
B(n,2)=n*2^(2n-2)-sum_(t=n+1)^(2n) binom(2n,t)(t-n).
```

The complement invariant lattice `W` is saturated in its monomial
lattice. The row supports of `v_I` and `W` are disjoint, so their product
coefficient vector is the tensor `v_I tensor W`. If a rational scalar
combination of its basis vectors has integral coordinates, Bézout on
the coordinates of primitive `v_I` recovers integral coordinates in
`W`; hence this fixed-support product lattice is saturated. In the raw
monomial Euclidean metric,

```text
Gram(v_I tensor basis(W)) = ||v_I||^2 Gram(W),
covol(v_I tensor W)=||v_I||^a covol(W),  a=rank(W).
```

For `(n,d)=(3,4)`, `W` has rank five and fixed covolume. Thus the
fixed-`I` rank-five lattice has exact per-rank logarithmic covolume
exponent `B(3,2)*2^6=18*64=1152`. It is far above the coefficient-floor
target 90, by a factor `12.8`. The full-kernel average
`171600/341 = 503.2258...` is also above 90 by a factor about `5.59`.

## Larger fixed-support factors

More generally, choose a subset `S` of size `n*d0`, with `2<=d0<=d`,
and tensor the entire local coherent kernel `K_S` on `S` with the
constant full invariant lattice on `S^c`, of degree `d-d0`. This is a
different family from multiplying `L_I` by a complement coherent
kernel. Its rank is `h(n,d0) r(n,d-d0)`. The full profile contains
`2^(m-nd0)` choices of outside core cuts for each cut internal to `S`,
so the determinant bound for `K_S` gives the per-rank upper exponent

```text
2^(m-nd0) * B(n,d0)/h(n,d0).
```

At `d0=2`, `K_S` has rank one and this is exactly the fixed-support
theorem above. At `d0>2` it remains only an upper bound, since the
specialized integral covolume of `K_S` is not known exactly. A finite
scan for `2<=n<=6, 2<=d0<=8` finds `B(n,d0)/(h(n,d0) 2^(nd0))`
decreasing with `d0`; it does not make this support-family upper bound
competitive with `M(3,4)=90` on the scanned range.

These exact fixed-support statements do not extend automatically to a
union of different supports: cross-support Gram blocks need not be
orthogonal, and the saturated lattice can acquire an index. A larger
support family that beats 90 would need an exact Gram/minor content
calculation or a globally defined rational subspace with a controlled
saturated integral lattice over the full core profile. The one-cut
support counts alone do not provide that gluing statement.

The subsequent [global index theorem](global_primitive_local_generator_index.md)
now bounds the index of the full primitive generator family by
`exp(o(w))` under these profile hypotheses. It does not give the
missing Gram or shortest-vector estimate: the
[redundant-generator model](redundant_local_generators_height_obstruction.md)
shows why exact generation and the displayed column heights alone
do not imply a vector below 90.

The accompanying exact finite audit is
[check_local_generator_fixed_support.py](check_local_generator_fixed_support.py).
Its profile hypotheses are that each core norm has logarithm
`w+o(w)`, bracket residue/correction products have logarithmic size
`o(w)`, and `(n,d)` stays fixed as `w` grows. The local determinant
estimate `B(n,d0)` is used as an upper bound for `d0>2`; only `d0=2`
has the exact fixed-support covolume calculation here.
