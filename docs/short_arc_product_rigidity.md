# Higher product rigidity on endpoint arcs

This note proves finite arithmetic restrictions for arbitrary lattice
configurations. The full uniform endpoint count remains unproved.
Products here allow repeated points, and equality of multisets includes
their multiplicities.

The later [angular moment argument](angular_moment_product_rigidity.md)
strengthens these conclusions: `k C^2 <= 8` suffices for rigidity
through length `k`, at every radius. In particular, `C=1/2` gives
length thirty-two. The polynomial estimates below remain valid but
are no longer the strongest bounds in this project.

The subsequent [linear-allocation theorem](linear_allocation_affine_rigidity.md)
removes the product-length bound altogether for `C<=sqrt(2)`, with
arbitrary units and no radius threshold. It still does not give a
radius-independent count of the original points.

## Results

Let `S` be the Gaussian integers on an arc of a circle of radius `R>0`,
of length at most `C sqrt(R)`.

1. If `C <= 1/2`, equality of two products of eight elements of `S`
   implies equality of their multisets. Consequently the same holds
   for products of any length from one through eight.
2. If also `R >= 1024`, the conclusion holds through length nine.
3. For every positive integer `k`, `C <= 1/(2k)` suffices for the
   conclusion through length `k`, with no radius threshold.

These are multiplicative `B_k` properties. They do not assert that
`S` has bounded cardinality. Empty sets are harmless; if `S` is
nonempty, its nonzero Gaussian points imply `R>=1`.

## Polynomial coefficient estimate

Suppose `z_1,...,z_k` and `w_1,...,w_k` are on the arc, with the
same product `p`. Form their monic polynomials `F,G`, and let
`e_r(z),e_r(w)` be the elementary symmetric coefficients.
Equal moduli give the exact identities

```text
e_(k-r)(z) = p conjugate(e_r(z)) / R^(2r),
e_(k-r)(w) = p conjugate(e_r(w)) / R^(2r).
```

In particular, equality of the low-index coefficients forces equality
of their reflected coefficients. Suppose `e_j(z)=e_j(w)` for `j<r`,
where `1 <= r <= floor(k/2)`. The constant coefficients already agree.
Then `H=F-G` is divisible by `X^r` and has degree at most `k-r`.
The polynomial `Q=H/X^r` has degree at most `m=k-2r`; its coefficient
of degree `m` is `(-1)^r (e_r(z)-e_r(w))`.

Let `q` be the midpoint of the arc on the circle, let `Delta` be its
angular width, and set `t=Delta/2`, `d=Rt`. Every root is at distance
at most `d` from `q`, and `|q|=R`. Repeated differentiation of a
product of `k` linear factors gives

```text
|F^(j)(q)|/j! <= binom(k,j) d^(k-j),
|G^(j)(q)|/j! <= binom(k,j) d^(k-j).
```

Apply the product rule to `Q(X)=X^(-r) H(X)` at the nonzero point
`q`. The normalized derivative of `X^(-r)` of order `s` has modulus
`binom(r+s-1,s) R^(-r-s)`. Since the normalized derivative of `Q`
of order `m` is its degree-`m` coefficient, this proves

```text
|e_r(z)-e_r(w)|
 <= 2 R^r t^(2r)
      sum_(l=0)^(k-2r) binom(k,2r+l) binom(r+l-1,l) t^l
 <= 2 (C^2/4)^r
      sum_(l=0)^(k-2r) binom(k,2r+l) binom(r+l-1,l)
                         (C/(2 sqrt(R)))^l.                 (*)
```

The second inequality uses `Delta <= C/sqrt(R)` and nonnegative
coefficients. Each coefficient difference is a Gaussian integer.
Whenever the displayed bound is strictly less than one, it is zero.
Induction through `floor(k/2)`, followed by the reflected identities,
therefore makes all coefficients equal. Equality `F=G` proves equality
of the root multisets, including repeated roots.

## A better first-coefficient bound

If `k Delta < 2 pi`, choose arguments of all roots in the centered
interval `[-Delta/2,Delta/2]` relative to `q`. Equality of the products
implies that the two argument sums differ by an integer multiple of
`2 pi`. Their difference has magnitude at most `k Delta < 2 pi`,
so those sums agree exactly.

Using `|exp(i x)-1-i x| <= x^2/2`, cancellation of the constant and
linear terms now gives

```text
|e_1(z)-e_1(w)| <= k R Delta^2/4 <= k C^2/4.          (**)
```

Thus `k C^2/4 < 1` makes the first coefficients equal. This step
requires the stated angular lift condition; estimate (*) does not.

## Eight factors at C=1/2

Take `k=8`. Since `R>=1`, `Delta<=1/2`, and `k Delta<=4<2 pi`.
Estimate (**) is at most `1/2`, so `e_1` agrees. In (*) we may use
`C^2/4=1/16` and `C/(2 sqrt(R))<=1/4`. The exact bounds are

| r | Upper bound in (*) |
|---|---:|
| 2 | 26565/32768 |
| 3 | 275/16384 |
| 4 | 1/32768 |

Every bound is less than one, proving the first assertion. To handle
shorter products, pad both sides with the same point of `S` until
they have length eight, apply the result, and cancel the padding.

## Nine factors for R>=1024

Now `k=9`, `t<=1/128`, and (**) is at most `9/16`. The angular lift
condition also holds. In (*) the exact bounds for `r=2,3,4` are

| r | Upper bound in (*) |
|---|---:|
| 2 | 2198751808323/2199023255552 |
| 3 | 88968581/2147483648 |
| 4 | 289/1048576 |

They are all strictly less than one. Reflected coefficients finish
the proof for nine, and padding gives all smaller lengths.

## Every fixed product length on a sufficiently short arc

For `k=1` the assertion is immediate. For `k>=2`, take
`C<=1/(2k)` and apply (*) from `r=1` onward. The elementary bounds

```text
binom(k,2r+l) <= binom(k,2r) binom(k-2r,l),
binom(r+l-1,l) <= 2^(r+l-1)
```

bound (*) by

```text
(C^2/2)^r binom(k,2r) (1+C/sqrt(R))^(k-2r)
 <= (k^2 C^2/2)^r exp(k C)/(2r)!
 <= (1/8)^r exp(1/2)/(2r)! < 1.
```

Here `R>=1`, `1+x<=exp(x)`, and `binom(k,2r)<=k^(2r)/(2r)!`
were used. This proves all coefficient equalities by induction, with
no assumption about unwrapping the arguments.

## What this contributes to the uniform problem

At `C=1/2`, an `M`-point cluster now has exactly `binom(M+7,8)`
distinct unordered eight-factor products. This excludes all nontrivial
balanced multiplicative relations of product length at most eight,
even with repeated indices. It strengthens the pair-product result
without assumptions on the conductor patterns or a limiting template.

Those products have radius `R^8`, however, and occupy an arc of length
at most `8 C R^(15/2)`. Relative to the new radius this is exponent
`15/16`. Counting them by a purported endpoint theorem would therefore
be invalid. Moreover arbitrarily large abstract multiplicative `B_k`
sets exist for every fixed `k`. The new rigidity does not supply the
missing conductor bonus or a cardinality bound for the original arc.

The proof is in prose. The constants in the two tables were checked
with exact rational arithmetic; no Lean formalization is claimed.
