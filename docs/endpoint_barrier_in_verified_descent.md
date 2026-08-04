# Endpoint barrier in the verified Mertens/descent proof

## Conclusion

The existing certified descent architecture cannot prove a uniform constant bound for endpoint-scale arcs merely by optimizing its parameters.

Its final dichotomy gives

```text
M <= 4 s,
```

but one of its hypotheses is

```text
4 m1 <= s,
```

where `(m0,m1]` is the prime interval used by the Mertens dichotomy. The proof also needs the interval to carry an increasing amount of weighted prime mass. Consequently `m1 -> infinity`, hence `s -> infinity`, and the conclusion `M <= 4s` cannot become uniform.

This identifies the precise theorem that must be strengthened if the verified descent route is to reach the endpoint conjecture.

---

## 1. Where the final cardinality bound comes from

`GaussianChain/MertensDichotomy.lean` proves a saturated dichotomy of the form

```text
card_le_four_mul_s_of_mertens_interval_log:
  ...
  hlarge     : 4*s <= M
  hsmallPrime: 4*m1 <= s
  ...
  -> M <= 4*s.
```

The missing-prime branch gives `M <= 2s`; the prime-divisor descent branch gives `M <= 4s`.

Thus a uniform endpoint bound through this theorem would require choosing `s=O(1)`.

But `hsmallPrime` would then force `m1=O(1)`.

---

## 2. Why a bounded prime interval cannot close the threshold

The Mertens contribution available from `(m0,m1]` is

```text
log m1 - log m0 - O(1).
```

At endpoint scale `L = C sqrt(R)`, the determinant threshold compares a left side of order

```text
s(2s+1) log(L^2) ~ 2s^2 log R
```

with a right side containing

```text
s^2 log N
```

plus the weighted prime mass. Since `N=R^2`, the leading `2s^2 log R` terms are critical and the lower-order prime mass is what supplies the strict inequality after rounding, missing-prime tolls, and fixed constants.

A fixed interval `(m0,m1]` supplies only a fixed amount of mass. It cannot absorb errors uniformly as `R` grows within the present inequality. The certified schedule therefore takes

```text
m0 ~ sqrt(log R),
m1 ~ (log R)^theta,
s  ~ q log R / log log R.
```

In particular both `m1` and `s` diverge.

Even if the threshold arithmetic were sharpened substantially, the structural hypothesis `4m1 <= s` alone prevents fixed `s` once the prime interval grows.

---

## 3. Method-barrier proposition

Any proof that uses the current theorem

```text
card_le_four_mul_s_of_mertens_interval_log
```

with prime intervals whose upper endpoints satisfy `m1(R) -> infinity` necessarily yields a bound tending to infinity:

```text
M <= 4s(R),
4m1(R) <= s(R)
```

implies

```text
M-bound >= 16m1(R) -> infinity.
```

Therefore no choice of `q`, `theta`, `m0`, or other asymptotic schedule inside the current theorem can prove `M=O(1)`.

This is a barrier theorem for the method, not merely a failure of the present constants.

---

## 4. Exact strengthening that would unlock the endpoint

The best application-specific target is a replacement for the single-prime descent lemma in which the admissible prime size is not tied linearly to `s`.

### Large-prime saturated descent target

There should exist an absolute function `C0(s)` independent of `p` such that, for every prime divisor `p|N`, an endpoint-scale cluster satisfying the same circle and monotonicity hypotheses obeys

```text
M <= C0(s)
```

whenever the descended determinant inequality holds, without the hypothesis

```text
4p <= s.
```

Even the weaker condition

```text
p <= exp(O(s))
```

would materially improve the schedule. A condition allowing polynomially growing `p` with fixed `s` would permit a constant endpoint bound.

---

## 5. More plausible multi-prime formulation

Removing `4p<=s` for one prime may be impossible in the existing Vandermonde determinant. A more natural target is simultaneous descent by a product of primes.

Let

```text
Q = product_{p in P} p
```

for a set of prime divisors of `N`. Instead of pigeonholing and descending by one `p`, partition the cluster according to its complete residue/signature data at all primes in `P`, and descend one large class by the corresponding Gaussian divisor of norm approximately `Q`.

The desired tradeoff is

```text
number of classes <= C^|P|
```

while the geometric scale improves by `sqrt(Q)` in the split case or `Q` in the inert case.

The current proof pays for a prime `p` through a condition roughly `p<=s`. A simultaneous argument should pay only for the number of local states, potentially logarithmic in `Q`, rather than for the numerical size of the largest prime.

This is the clearest place where the common-conductor structure can improve the verified proof: every point carries compatible valuation choices at all primes dividing the same norm.

---

## 6. Concrete next lemma

A finite combinatorial version sufficient to test the idea is:

> Let `A` be a cluster of Gaussian integers of common norm `N`. For a finite set `P` of split primes dividing `N`, assign to each point its orientation/valuation state at every `p in P`. Prove that either `|A|` is already bounded in terms of a fixed determinant order `s`, or there is a subset `A'` of size at least `|A|/C^|P|` sharing a common Gaussian divisor `d` with `Norm(d)` comparable to `product_{p in P} p`, such that division by `d` preserves injectivity and sends the arc to one of length reduced by `|d|`.

One then chooses `P` so that the product norm gives a large geometric descent while the state-count loss remains bounded or can be iterated with a contraction.

This target is substantially narrower than a new uniform Subspace Theorem and is directly connected to the already formalized descent proof.

---

## 7. Recommended research path

1. Extract the exact local state partition used in `PrimeDescent` for one split prime.
2. Tensor those partitions over several primes and prove the common-divisor statement.
3. Quantify the class-count loss versus the product descent gain.
4. Determine whether a fixed determinant order `s` can survive one multi-prime descent step.
5. Only if that succeeds, replace the single-prime branch in `MertensDichotomy` and rerun the endpoint threshold.

The Ru--Vojta branch remains useful as a conceptual certificate that a coefficient below `1/2` is plausible, but the most credible constructive route is now a multi-prime upgrade of the verified Gaussian descent machinery.
