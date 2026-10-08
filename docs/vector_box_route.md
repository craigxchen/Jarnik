# Vector-box route

## 1. Exact exponent-box model

Let

```text
N = prod_i p_i^{e_i}
```

with each split rational prime written

```text
p_i = pi_i * conjugate(pi_i),
pi_i = |pi_i| exp(i theta_i),
0 < theta_i < pi/2.
```

Every Gaussian integer point of norm `N` is, up to a unit and the ramified/inert factors common to every representation,

```text
z(a) = prod_i pi_i^{a_i} conjugate(pi_i)^{e_i-a_i},
0 <= a_i <= e_i.
```

Its modulus is

```text
R = sqrt(N),
```

and its argument is

```text
arg z(a) = const + sum_i (2 a_i-e_i) theta_i mod 2pi.
```

Thus the lattice points are the image of the integer box

```text
B(e) = prod_i {0,...,e_i}
```

under the torus linear form

```text
Phi(a) = 2 sum_i a_i theta_i mod 2pi.
```

An arc of length at most `C sqrt(R)` corresponds to an angular interval of length

```text
delta = C / sqrt(R) = C N^{-1/4}.
```

The uniform arc problem is therefore equivalent to bounding

```text
#{a in B(e) : Phi(a) lies in I}
```

uniformly for every interval `I` of length `C N^{-1/4}`.

## 2. Why an arbitrary box theorem is false

No such statement holds for arbitrary real angles `theta_i`. If all `theta_i` are zero, every vector lies in the same interval. More generally, many very small or rationally related angles can concentrate exponentially many subset sums.

The arithmetic input is essential: each `theta_i` comes from a primitive Gaussian integer `pi_i = x_i+i y_i` with

```text
sin(theta_i) = y_i/sqrt(p_i),
```

and distinct `pi_i` are nonassociate Gaussian primes.

## 3. One-coordinate case

For one split prime `p = pi conjugate(pi)` with exponent `e`, the points are

```text
z_a = pi^a conjugate(pi)^{e-a},
0 <= a <= e,
```

and consecutive arguments differ by `2 theta`.

Since `pi=x+iy` is nonreal and primitive, `|y|>=1`, hence

```text
sin(theta) >= p^{-1/2}.
```

For `0 < theta <= pi/2`, this gives `theta >= p^{-1/2}`. Therefore the angular spacing is at least `2 p^{-1/2}`.

The endpoint interval has length `C p^{-e/4}`. Consequently its occupancy is at most

```text
1 + (C/2) p^{1/2-e/4}.
```

In particular, for `e>=3` this is bounded by a constant depending only on `C`, and for fixed small `e` the box itself has bounded size. Thus the one-prime case satisfies a uniform bound.

This simple argument also identifies the only serious difficulty in higher dimension: cancellation among several prime angles.

## 4. Difference-set reformulation

If `A subset B(e)` maps into an interval of length `delta`, then for every `a,b in A`,

```text
|| 2 sum_i (a_i-b_i) theta_i ||_{R/2piZ} <= delta.
```

Let

```text
D = A-A subset prod_i {-e_i,...,e_i}.
```

Every `d in D` yields a norm-one Gaussian rational

```text
u(d) = prod_i (pi_i/conjugate(pi_i))^{d_i}
```

satisfying

```text
|u(d)-1| << delta.
```

After clearing negative exponents, write

```text
X_d = prod_i pi_i^{max(d_i,0)} conjugate(pi_i)^{max(-d_i,0)},
Y_d = conjugate(X_d).
```

Then

```text
u(d)=X_d/Y_d,
```

and `u(d)` close to one is equivalent to the Gaussian integer

```text
X_d-Y_d = 2i Im(X_d)
```

being unusually small relative to `|X_d|`.

This is the arithmetic form of angular cancellation.

## 5. Precise target theorem

A sufficient vector-box theorem is:

> For every `C>0`, there exists `M(C)` such that for every finite family of pairwise nonassociate split Gaussian primes `pi_i`, every exponent box `B(e)`, and every interval `I` of length `C N^{-1/4}`, the fiber
>
> ```text
> {a in B(e): 2 sum_i a_i arg(pi_i) mod 2pi lies in I}
> ```
>
> has cardinality at most `M(C)`.

This is exactly equivalent to the original uniform arc bound, but it removes heights, exceptional sets, and Subspace-Theorem terminology.

## 6. A genuinely simpler intermediate target

The first nontrivial theorem to prove is a cancellation dichotomy:

> Given a set `A` in one endpoint interval, either
>
> 1. some coordinate projection of `A` has uniformly bounded multiplicity and induction on the remaining coordinates applies; or
> 2. there are two independent difference vectors `d,d' in A-A` whose associated Gaussian integers `X_d-X_d^*` and `X_{d'}-X_{d'}^*` have a nontrivial common divisor carrying a fixed positive proportion of the conductor mass.

The second branch should force a low-complexity multiplicative relation among the vectors. Iterating such relations would place `A` in a bounded-rank affine slice of the box.

This is the box analogue of the generalized-GCD mechanism, but stated entirely in finite combinatorics and Gaussian integer divisibility.

## 7. Immediate research tasks

1. Prove a two-coordinate theorem completely.
2. Quantify how a large fiber forces many repeated coordinate differences.
3. Convert repeated differences into common divisors of `X_d-X_d^*`.
4. Prove that sufficiently many independent such divisibilities force an affine relation among exponent vectors.
5. Use induction on affine rank to bound the fiber.

The determinant and adapted-basis routes are not used in this approach.
