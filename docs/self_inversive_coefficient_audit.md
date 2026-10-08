# Self-inversive coefficients on a short equal-radius arc

This note isolates what the polynomial with the arc points as roots gives
beyond the usual product or Vandermonde divisibility. It gives an exact
reduction. The self-inversive identities couple all coefficients, but their
Gaussian-prime content is exactly the complementary product content. The
termwise divisor/reflection identities and the short-arc estimate below do
not by themselves produce a logarithmic margin at a fixed number of roots.

## 1. The exact coefficient identity

Let `z_1,...,z_m` be nonzero Gaussian integers satisfying

```text
|z_j|=R,
```

where `R^2` is an integer. Write

```text
e_k = sum_{|S|=k} product_{j in S} z_j,
P(X)=product_j (X-z_j)=sum_{k=0}^m a_k X^k,
a_k=(-1)^(m-k)e_(m-k).
```

Since `conjugate(z_j)=R^2/z_j`, direct expansion gives, for
`0 <= k <= m`,

```text
e_m conjugate(e_k) = R^(2k)e_(m-k).                 (1)
```

Equivalently, the coefficients obey

```text
a_0 conjugate(a_(m-k)) = R^(2k)a_k.                (2)
```

This is the self-inversive relation for the circle of radius `R`. It is
an identity for every equal-radius root set, including a set with no short
arc. In particular, the high coefficients are determined by the low ones
and the constant coefficient; they are not a second independent source of
small integer values.

Taking absolute values also gives

```text
|e_(m-k)| = R^(m-2k)|e_k|.                            (1a)
```

for every `k` (with the reciprocal interpretation when `k > m/2`). Thus
even an aggregate over all coefficient magnitudes merely repeats the low
coefficient data with explicit powers of `R`; it does not create a second
height budget.

If `g` is a Gaussian gcd of all the `z_j`, put `z_j=g u_j`,
`d=|g|`, and `R'=R/d`. Then `e_k=g^k E_k`, and cancellation of
`|g|^(2k)` in (1) gives

```text
E_m conjugate(E_k) = R'^(2k) E_(m-k).               (3)
```

The primitive points `u_j` therefore satisfy the same identity. Their
angular diameter is unchanged, but the normalized endpoint constant is
improved from `C` to `C/sqrt(d)`, since

```text
Delta <= C/sqrt(R) = (C/sqrt(d))/sqrt(R').          (4)
```

Any coefficient argument should make this gcd reduction before comparing
heights. An inert or ramified common factor is included in `g`.

## 2. Prime valuations: what simultaneous symmetry actually says

Fix a Gaussian prime `pi` and put

```text
E=v_pi(R^2),       alpha_j=v_pi(z_j),
A=sum_j alpha_j=v_pi(e_m).
```

The equal-radius condition gives

```text
v_conjugate(pi)(z_j)=E-alpha_j.
```

Every term of `e_k` has `pi`-valuation at least

```text
mu_k = min_{|S|=k} sum_{j in S} alpha_j
     = sum of the k smallest alpha_j.                 (5)
```

Cancellation can only increase this valuation. Taking valuations in (1)
gives the exact complementary equation

```text
A + v_pi(conjugate(e_k)) = kE + v_pi(e_(m-k)).        (6)
```

Since `v_pi(conjugate(e_k)) = v_conjugate(pi)(e_k)`, (6), together with
its conjugate, simply transfers any extra cancellation in one coefficient
to the complementary coefficient. The baseline valuations (5) are the
valuation content of the individual products. Thus the complete system
(1), for all `k`, imposes compatibility of those contents but gives no
additional divisor from the reflection equations themselves. Formula (6a)
below is an exact statement whenever the coefficients involved are nonzero;
it does not assert that cancellations are absent.

More explicitly, if `delta_k(pi)=v_pi(e_k)-mu_k(pi)` denotes the excess
valuation over the termwise minimum, then the minimum identities imply

```text
delta_(m-k)(pi) = delta_k(conjugate(pi)).             (6a)
```

Thus the all-`k` system only pairs excesses at conjugate primes and
complementary indices; it does not sum them into an additional divisor.
If `e_k=0`, use the convention `v_pi(0)=infinity`. Equation (1) then
forces `e_(m-k)=0` as well (because `e_m` is nonzero), and (6a) is simply
not a finite excess-valuation statement for that pair.

For coefficient divisibility stated without choosing a factorization of the
individual terms, (2) gives the ideal divisor

```text
R^(2k)/gcd(R^(2k),a_0)  divides  conjugate(a_(m-k)),  (7)
```

prime by prime. After the Gaussian gcd is removed, this is the divisor
supplied by the reflection
equation itself and can be weaker than the full termwise minimum profile
(5); it should not be counted once for every coefficient as if the
divisors were independent.

## 3. What the short arc contributes

Let `c=R exp(i phi)` be the midpoint of the containing arc and let
`Delta` be its angular width. Write

```text
z_j=c exp(i delta_j),       |delta_j| <= Delta/2,
eta=exp(i sum_j delta_j).
```

Then

```text
e_k=c^k s_k,
s_k=sum_{|S|=k} exp(i sum_{j in S} delta_j),
e_m=eta c^m.
```

The simultaneous coefficient estimate is only

```text
|e_k - binom(m,k)c^k|
 <= binom(m,k) (k Delta/2) R^k
 <= binom(m,k) k C R^(k-1/2)/2.                    (8)
```

The normalized coefficients satisfy the exact reflection

```text
s_(m-k)=eta conjugate(s_k).                          (9)
```

Here `|eta-1| <= m Delta/2`. Hence (9) gives no second small quantity:
it reflects the same `s_k`, with the same `m`-dependent error. In
particular, for fixed `m`, the coefficient lattice spacing `1` is
overwhelmed by the `R^(k-1/2)` error in (8) as `R` grows. Shrinking
`Delta` to `C/sqrt(R)` therefore does not force a coefficient gap to
vanish.

## 4. Comparison with product rigidity

For two products with the same number `m` of arc points, equality of the
products identifies `e_m`. Applying (1) then shows that equality of the
low coefficients (`e_1,...,e_r`) automatically gives equality of the
reflected coefficients (`e_(m-1),...,e_(m-r)`). This is the coefficient
reflection used in the finite product-rigidity arguments. It can prove
finite `m`-factor uniqueness only after a separate angular-moment
cancellation argument has made the low coefficients equal, as in
[angular_moment_product_rigidity.md](angular_moment_product_rigidity.md).
Estimate (8) by itself grows with `R` and does not provide that equality
for an arbitrary root set. Even finite product rigidity does not bound the
number of original points: products of length `m` have radius `R^m`, so
their arc has normalized exponent `1-1/(2m)`.

For one root set, the reflection calculation above only recycles the
complementary prime weights already present in the termwise products. A
stronger result would require an additional restriction on simultaneous
cancellation in the `e_k`, or a new height bound for a nonzero coefficient
combination. The present audit does not rule out such extra input; it
records only that it is not supplied by self-inversiveness, the termwise
valuation profile, and estimate (8).

This is a scope audit rather than a disproof of the uniform endpoint bound.
