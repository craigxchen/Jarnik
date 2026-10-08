# Annihilator capacity is endpoint-critical

## Status

This note audits the power-separating annihilator theorem against the elementary pairwise separation already forced by an endpoint cluster.

The conclusion is decisive:

> The annihilator-capacity theorem cannot by itself prove the uniform endpoint bound. Its width threshold is at most one half of the conductor height, while every nontrivial support drawn from an endpoint cluster already has width at least one half of the conductor height up to bounded error.

Thus the framework reproduces the same endpoint equality as the earlier methods, although in a genuinely global language.

## 1. Setup

Let

```text
N = product_j p_j^(e_j),
H = log N,
```

and let `A` be exponent vectors for distinct points in an angular interval of width

```text
delta = C exp(-H/4).
```

For two vectors `a,b`, define the conductor distance

```text
D(a,b) = sum_j |a_j-b_j| log p_j.
```

The quotient of the corresponding norm-one Gaussian rationals is within `O_C(delta)` of `1`. Clearing Gaussian denominators gives the standard endpoint separation

```text
D(a,b) >= H/2 - O_C(1).                         (1.1)
```

## 2. Support width dominates every pair distance

For a signed group-algebra element

```text
P = sum_a c_a [a]
```

with support `S`, define

```text
W(P) = sum_j (max_(a in S) a_j - min_(a in S) a_j) log p_j.
```

For every pair `a,b in S`, coordinatewise

```text
|a_j-b_j| <= max_(x in S) x_j - min_(x in S) x_j.
```

Hence

```text
D(a,b) <= W(P).                                  (2.1)
```

If `P` is nontrivial, its support contains two distinct cluster points. Combining (1.1) and (2.1),

```text
W(P) >= H/2 - O_C(1).                            (2.2)
```

This lower bound is independent of the moment order or coefficient choice.

## 3. The power-separating threshold cannot exceed one half

The power-separating annihilator theorem applies to an order-`q` phase-moment annihilator with support size `s` when asymptotically

```text
W(P)/H < q / (2(s-1)).                           (3.1)
```

For `s` distinct phase values, a nonzero coefficient vector cannot annihilate the scalar moments

```text
1,t,...,t^(s-1)
```

because the corresponding `s x s` Vandermonde matrix is invertible. Therefore

```text
q <= s-1.                                        (3.2)
```

Consequently

```text
q / (2(s-1)) <= 1/2.                             (3.3)
```

The best possible case is `q=s-1`, where the capacity threshold is exactly

```text
W(P) < H/2 - margin.                             (3.4)
```

But (2.2) says that every nontrivial support already satisfies

```text
W(P) >= H/2 - O_C(1).
```

Thus no fixed positive asymptotic margin is available.

## 4. Interpretation

The analytic gain from canceling `q` phase moments is

```text
q H / 4.
```

To eliminate the exact-zero branch by powers, one must test up to power `s-1`. The arithmetic denominator cost at the largest power is

```text
(s-1) W(P) / 2.
```

Since `q <= s-1`, the analytic gain can beat the arithmetic cost only when

```text
W(P) < H/2.
```

Endpoint pairwise separation forbids precisely that.

This is the same critical exponent in a new form:

```text
arc smallness       = exp(-H/4),
quotient height      = H/2,
power separation    = s-1,
moment order         <= s-1.
```

All inequalities meet at equality.

## 5. Consequence for the research program

The annihilator framework remains useful as a language for global signed relations, but it is not a contradiction engine at the endpoint unless it gains something beyond ordinary support width.

A viable repair would need an arithmetic denominator invariant strictly smaller than the coordinate-box width for a signed combination. Examples worth investigating are:

1. cancellation in the common denominator after summing the monomials;
2. ideal-theoretic content shared by the numerator terms;
3. a resultant or norm that uses several power sums simultaneously and has denominator cost below `(s-1)W/2`;
4. a sparse Newton-polytope invariant based on mixed width rather than bounding-box width.

Without such an improvement, searching for annihilators cannot cross the endpoint barrier.
