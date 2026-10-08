# Primitive-chord divisor route

## Status

This note tests a new route to the uniform endpoint bound that does not use the global determinant inequality or Ru--Vojta exceptional sets.  It produces a clean structural lemma for every chord from a fixed lattice point, but its sharp packing inequality is again critical at exponent `1/2`.  The route does not currently prove the uniform bound.

## 1. Primitive chord factorization

Let `z0,z` be distinct Gaussian integers with

```text
|z0| = |z| = R.
```

Write

```text
z-z0 = g w,
```

where `g` is a positive rational integer and `w` is primitive as a lattice vector, meaning that the real and imaginary parts of `w` are coprime.

The circle equation gives

```text
z0 conjugate(w) + conjugate(z0) w = -g w conjugate(w).    (1.1)
```

Reducing modulo `w` shows

```text
w divides z0 conjugate(w).                                (1.2)
```

Let `delta = gcd(w,conjugate(w))` in `Z[i]`.  For primitive `w`, `delta` is a unit or an associate of `1+i`.  Therefore

```text
v = w/delta
```

is a Gaussian divisor of `z0`.  Thus every chord direction is, up to the harmless ramified factor, the direction of a Gaussian divisor of the fixed base point.

This is stronger than the usual exponent-box description because it attaches a divisor directly to each geometric chord rather than to the global point representation.

## 2. Thin-sector geometry

Fix one endpoint `z0` of an arc of length at most

```text
L <= C sqrt(R).
```

For every other point `zj`, let `wj` be the primitive direction of the chord `zj-z0`.

The angle between two chords from the same endpoint is half the difference of the corresponding central angles.  Hence

```text
angle(wi,wj) <= C/(2 sqrt(R)).                             (2.1)
```

Also

```text
|wj| <= |zj-z0| <= C sqrt(R).                             (2.2)
```

Consequently

```text
|det(wi,wj)|
  = |wi| |wj| |sin angle(wi,wj)|
  <= O_C(sqrt(R)).                                        (2.3)
```

Because distinct rays from `z0` meet the circle in at most one further point, distinct points give distinct primitive directions.  Thus the endpoint problem becomes a packing problem for distinct primitive Gaussian divisors of `z0` in a sector of width `O_C(R^(-1/2))`, with pairwise determinants at most `O_C(sqrt(R))`.

## 3. Common-divisor interpretation

For Gaussian integers `u,v`, a Gaussian gcd divides

```text
u conjugate(v) - conjugate(u) v = 2 i det(u,v).
```

Therefore

```text
Norm(gcd(u,v)) <= O(|det(u,v)|^2).                         (3.1)
```

Applied to the primitive chord divisors, (2.3) gives only

```text
Norm(gcd(vi,vj)) <= O_C(R).                               (3.2)
```

The individual norms can also be of order `R`.  Thus the pairwise gcd estimate lies exactly at the half-conductor threshold: two large divisors may share half of the available logarithmic prime mass while their lcm still divides `z0`.

## 4. Why set-system compression does not create a gap

Represent the split-prime factors of `z0` as a weighted universe of total logarithmic mass `2 log R`.  A primitive divisor corresponds to an oriented subset of this universe.  If its norm is of order `R`, its subset has mass about `log R`.  The determinant bound permits pairwise gcd mass up to about `(1/2) log R`.

This is the Plotkin-critical set-system regime.  There are arbitrarily large abstract families of half-mass subsets with pairwise intersection mass one quarter of the universe.  Equivalently, random balanced binary patterns attain

```text
individual mass        = (1/2) total mass,
pairwise intersection  = (1/4) total mass.
```

The primitive-chord inequalities do not distinguish these equality models from genuine Gaussian divisor families.  Hence no argument using only

```text
- divisor norms,
- pairwise Gaussian gcds,
- pairwise chord determinants, or
- the width of the common angular sector
```

can force a uniform cardinality bound.

## 5. A partial range where the route works

If every primitive chord direction satisfies

```text
|wj| <= K R^(1/4),                                       (5.1)
```

then (2.1) gives

```text
|det(wi,wj)| <= O_C(K^2).                                (5.2)
```

A set of distinct primitive integer directions with all pairwise determinants bounded by a fixed number has bounded cardinality.  To see this, fix one direction `u`.  For each nonzero integer `q=det(u,v)` with `|q|<=B`, all vectors with determinant `q` lie in one affine progression parallel to `u`; pairwise determinant boundedness allows only `O(B/|q|)` members in that progression.  Summing over `q` gives `O(B log B)` directions.

Therefore the endpoint cluster is uniformly bounded in the subcase (5.1).  Any counterexample to the full theorem must contain primitive chord directions of size substantially larger than `R^(1/4)`.

This is a genuine but narrow reduction.

## 6. Exact failure point

For unrestricted endpoint chords, `|wj|` may be of order `R^(1/2)`.  The best sector estimate then gives

```text
|det(wi,wj)| <= O_C(R^(1/2)),                             (6.1)
```

not a constant.  The logarithmic divisor/gcd formulation of (6.1) is exactly the balanced-code threshold described in Section 4.  There is no strict exponent gain to iterate.

The route would become viable only with an additional arithmetic statement coupling the determinant

```text
det(wi,wj)
```

to the complementary factors `z0/vi` and `z0/vj` in a way that improves (6.1) by a fixed power of `R`.  No such identity was found.  The circle equation supplies divisibility, but after the ramified gcd is removed its strongest uniform consequence is already recorded in (1.2).

## 7. Conclusion

The primitive-chord factorization is a useful new normal form:

> every chord from a fixed circle point has primitive direction equal, up to `1+i`, to a Gaussian divisor of the base point.

It proves a uniform bound when all primitive chord directions are `O(R^(1/4))`, but the general packing problem reaches the same half-mass equality geometry as the earlier exponent-box models.  Thus the route fails concretely at the endpoint unless one introduces new arithmetic information beyond norms, gcds, and angular spacing.
