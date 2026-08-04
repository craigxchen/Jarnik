# Bounded split-prime support is uniformly bounded

## Status

This note proves a genuine global reduction for the endpoint problem.  It does
not prove the final uniform bound, but it eliminates every family whose common
norm has a bounded number of split rational-prime factors.  Consequently any
hypothetical unbounded endpoint family must have genuinely diffuse, unbounded
split support.

The proof combines a half-prime-power pigeonhole with one fixed-order Ramana
subcritical determinant.  It is uniform in the identities and sizes of the
split primes.

## 1. Setup

Let

```text
z_1,...,z_M in Z[i]
```

be distinct points of common Gaussian norm `N`, lying on an arc of length

```text
L <= C N^(1/4).
```

After dividing all points by their common rational Gaussian factor, inert
primes and the ramified prime contribute no moving orientation.  Write the
remaining split part as

```text
N = product_{j=1}^r p_j^(e_j),
```

where `p_j = pi_j conjugate(pi_j)` and `p_j == 1 mod 4`.

Fix `K` and assume `r <= K`.

## 2. Half-power descent at one prime-power block

For a fixed block `p^e`, every point has a unique exponent

```text
0 <= a_i <= e
```

such that its `p`-part is

```text
pi^(a_i) conjugate(pi)^(e-a_i).
```

Put

```text
q = ceil(e/2).
```

For every `a_i`, either

```text
q <= a_i
```

or

```text
q <= e-a_i.
```

Hence at least half the family is divisible by one of

```text
pi^q,  conjugate(pi)^q.
```

This elementary majority statement is formalized in
`GaussianChain/PrimePowerMajority.lean`.

Dividing that subfamily by the chosen Gaussian factor is a similarity.  It
changes

```text
N  -> N' = N / p^q,
L  -> L' = L / p^(q/2),
```

and preserves distinctness and circular order.

## 3. Choose the largest prime-power block

Among at most `K` blocks there is one satisfying

```text
p^e >= N^(1/K).
```

For `q = ceil(e/2)`, put `Q = p^q`.  Then

```text
Q >= (p^e)^(1/2) >= N^(1/(2K)).                 (3.1)
```

The descended squared diameter is at most

```text
B' = C^2 N^(1/2) / Q.                            (3.2)
```

## 4. Fixed-order subcritical inequality

Take a fixed integer `s >= K`.  An injective block of `2s+1` descended points
would contradict Ramana's determinant as soon as

```text
(B')^(s(2s+1)) < (N/Q)^(s^2).                    (4.1)
```

Using (3.2), condition (4.1) is equivalent, up to the harmless natural-floor
rounding in the formalized theorem, to

```text
(s+1) log Q > (1/2) log N + 2(2s+1) log C.       (4.2)
```

By (3.1), the coefficient of `log N` on the left is at least

```text
(s+1)/(2K).
```

Since `s >= K`, this is strictly larger than `1/2`.  Therefore (4.2) holds for
all sufficiently large `N`, with a threshold depending only on `C` and `K`.

## 5. Cardinality consequence

The half-power class has at least `M/2` points.  If

```text
M >= 2(2s+1) = 4s+2,
```

it contains an injective `2s+1`-point block, contradicting (4.1).  Taking
`s=K` gives, for all sufficiently large `N`,

```text
M <= 4K+1.                                             (5.1)
```

The bounded range of smaller `N` contributes another finite constant depending
only on `C,K`.

Thus:

> **Bounded-support theorem.**  For every `C` and `K`, endpoint arcs whose
> common norm has at most `K` split rational-prime factors contain at most
> `B(C,K)` lattice points.

## 6. Consequence for a hypothetical counterexample

Any sequence of endpoint clusters with cardinality tending to infinity must
satisfy

```text
number of distinct split rational-prime factors -> infinity.
```

More quantitatively, applying the argument with `s` comparable to the cluster
size shows that the split support cannot remain substantially smaller than the
cluster.  This isolates the genuinely hard case: many low-weight prime layers
whose orientation cuts can approximate the balanced Plotkin regime.

## 7. Relation to the new prime-layer surplus

For `2s+1` selected points, a split-prime threshold cut of size `k` contributes

```text
s(s+1) - k(2s+1-k) = (k-s)(k-s-1)
```

to the conductor-prime valuation of the Ramana determinant.  This identity and
its integral gap away from `k=s,s+1` are formalized in
`GaussianChain/PrimeLayerSurplus.lean`.

The bounded-support theorem handles the regime where one prime-power block is
large.  The remaining diffuse regime must now be attacked by showing that many
small blocks cannot all maintain balanced cuts while their Gaussian arguments
place the same points in one interval of width `C N^(-1/4)`.
