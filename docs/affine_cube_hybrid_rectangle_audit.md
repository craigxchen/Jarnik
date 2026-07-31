# Affine-cube / hybrid-rectangle route audit

## Status

This note studies the affine-cube approach to the uniform endpoint problem. It does not prove the uniform bound. It identifies one genuine mechanism—cube completion inside the conductor box—and pushes it to a sharp entropy/arithmetic split.

## 1. Conductor-box model

After removing common Gaussian factors, write

```text
N = product_j p_j^(e_j)
```

with one Gaussian prime `pi_j` chosen above each split prime. Every lattice point on the circle is represented by an exponent vector

```text
a = (a_j),   0 <= a_j <= e_j,
```

and has angular coordinate

```text
Phi(a) = sum_j (2 a_j - e_j) theta_j mod 2 pi,
theta_j = arg(pi_j).
```

An endpoint cluster is a set `A` whose `Phi`-image lies in an interval of width `delta ~ N^(-1/4)`.

## 2. Cube completion

If `a,b,c` lie in the conductor box and

```text
d = a + b - c
```

also lies in the box, then `d` represents another lattice point on the same circle and

```text
Phi(d) = Phi(a) + Phi(b) - Phi(c) mod 2 pi.
```

Thus if `a,b,c` lie in one interval of width `delta`, the completed point lies in an interval of width at most `3 delta` after choosing compatible lifts.

This is an exact higher-order closure statement. It does not use determinants or pairwise cut metrics.

## 3. Why arbitrary cube abundance does not follow

For three vectors, the condition `a+b-c` lies in the box is coordinatewise. In one coordinate, after sorting three values, one of the three points is the median and subtracting that point gives an in-range completion. But the median point may depend on the coordinate. Hence there need not be one global choice of `c` that works in every coordinate.

In unbounded dimension, coordinate order patterns can be prescribed independently. Therefore a fixed-size set may be affine-cube-free even though all its angular images are concentrated. Pure finite-dimensional Ramsey or density arguments do not apply because the conductor-box dimension is unbounded and the cluster is sparse.

## 4. Hybrid rectangles from two points

For two box points `a,b`, define the coordinate rectangle

```text
R(a,b) = product_j [min(a_j,b_j), max(a_j,b_j)] cap Z.
```

Every point in `R(a,b)` represents a lattice point on the same circle. The cardinality is

```text
|R(a,b)| = product_j (|a_j-b_j| + 1).
```

The endpoints `a,b` are two opposite vertices. The rectangle contains many affine cubes and all coordinatewise hybrids between the two points.

If `M_N(delta)` is the maximum number of circle points in an angular interval of width `delta`, then pigeonholing the angles of all rectangle points around the circle gives

```text
M_N(delta) >= c * delta * product_j (|a_j-b_j| + 1)
```

for an absolute constant `c>0` (with the usual harmless endpoint convention for circular intervals).

Consequently, if `a,b` belong to a cluster realizing `M_N(delta)`, maximality implies

```text
product_j (|a_j-b_j| + 1) <= C * M_N(delta) / delta.
```

At the endpoint `delta ~ N^(-1/4)`, this becomes

```text
sum_j log(|a_j-b_j|+1) <= (1/4) log N + log M_N(delta) + O(1).
```

This is the main new inequality supplied by the affine-cube route.

## 5. Compare with pairwise algebraic separation

For two distinct endpoint points, the norm-one ratio has height

```text
h(a/b) = sum_j |a_j-b_j| log p_j.
```

Elementary algebraic separation and endpoint closeness force, up to constants,

```text
sum_j |a_j-b_j| log p_j >= (1/2) log N - O_C(1).
```

Thus any close pair in a maximal endpoint cluster must satisfy both

```text
sum_j |a_j-b_j| log p_j >= (1/2) log N - O_C(1),
```

and

```text
sum_j log(|a_j-b_j|+1) <= (1/4) log N + log M_N(delta) + O_C(1).
```

The first quantity is arithmetic height. The second is recombination entropy.

## 6. The entropy-rich regime

If a pair satisfies

```text
sum_j log(|a_j-b_j|+1) > (1/4 + epsilon) log N,
```

then its hybrid rectangle alone produces at least `N^epsilon` points in some endpoint-scale interval. In particular, a cluster that is maximal but only bounded in size cannot contain such a pair for large `N`.

So any hypothetical bounded-size extremal cluster must be entropy-poor: its required half-conductor separation is carried by coordinates for which `log p_j` is large compared with the combinatorial choice `log(|a_j-b_j|+1)`.

## 7. Concrete obstruction

There is no universal inequality of the form

```text
sum_j log(|a_j-b_j|+1) >= c * sum_j |a_j-b_j| log p_j
```

with fixed `c>0`. A single exponent difference of one at a very large split prime contributes `log 2` to rectangle entropy but `log p` to arithmetic height.

More generally, many distinct large primes can carry half the conductor height while the rectangle has only `2^r` vertices, with

```text
r << log N.
```

Therefore hybrid-cube proliferation cannot by itself handle the large-information-per-coordinate regime.

This is not the old max-cut obstruction. It is a different and exact failure: combinatorial cube entropy counts choices, while the endpoint scale is governed by logarithmic prime size.

## 8. Relation to inverse Littlewood--Offord theory

The residual regime can be phrased as follows. A family of subset sums of large Gaussian-prime angles has more concentration at scale `N^(-1/4)` than its raw recombination entropy predicts. An inverse Littlewood--Offord theorem would then seek additive structure among those angles.

However, standard inverse results are parameterized by concentration probability relative to `2^r`. Here that probability may be exponentially small in `r`, and the desired structure must be weighted by `sum log p_j`, not by the coordinate count. Standard Freiman and Balog--Szemeredi--Gowers theorems only become useful after one has established large additive energy or small doubling; endpoint angular concentration of a sparse fixed-size set does not supply either automatically.

## 9. Final assessment

The affine-cube approach yields a genuine dichotomy:

```text
recombination-rich pair
    -> many hybrid points
    -> endpoint cluster growth;

recombination-poor pair
    -> half-conductor height stored in large-information coordinates.
```

The first branch is controlled by the hybrid-rectangle inequality. The second branch is not accessible to pure cube combinatorics.

A completion would require a new weighted inverse-concentration theorem for Gaussian-prime angles that converts concentration beyond the hybrid-entropy baseline into either:

1. a bounded set of heavy prime coordinates, enabling descent; or
2. a bounded-rank additive relation among the angles, strong enough to contradict Gaussian unique factorization or a simultaneous logarithmic-form lower bound.

Thus the affine-cube route does not by itself prove the uniform bound, but it produces a concrete, non-determinantal reduction and isolates the exact residual arithmetic regime.
