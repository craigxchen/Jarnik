# A disjoint-support balanced elliptic six-curve

The mixed translation/reflection case admits a model in which all four
common arguments have different target values.  Consequently its
twenty-five degree-one boundary pullbacks are supported at twenty-five
distinct points.  Repeated boundary support is therefore not forced by
the mixed branch-pair pattern.

## 1. The pencil

Set

```text
P_0=x^2,                  Q_0=x-1,
P_1=-60-73x,              Q_1=x^2-11x-38,
f_lambda=(P_0+lambda P_1)/(Q_0+lambda Q_1).             (1)
```

The determinant is

```text
P_0Q_1-P_1Q_0=(x+1)(x-3)(x-4)(x-5).                   (2)
```

Thus every map in the pencil takes common values at `x=-1,3,4,5`.
Those values are respectively

```text
-1/2,       9/2,       16/3,       25/4,               (3)
```

and are pairwise distinct.

Choose

```text
lambda_a=0,       lambda_b=-240/5329,
lambda_c=8/257.                                          (4)
```

The critical points and branch values are

| map | critical points | branch values |
|---|---|---|
| `f_a` | `0,2` | `u=0,v=4` |
| `f_b` | `-120/73,402794/639337` | `u=0,s=212248804/67144321` |
| `f_c` | `14/5,3244/1069` | `v=4,t=208624/46513` |

All four branch values are distinct.  The branch pairs again have the
mixed path pattern

```text
f_a:{u,v},       f_b:{u,s},       f_c:{v,t}.            (5)
```

None of the common values (3) is a branch value.

## 2. Normalized common arguments

Precompose every map by

```text
x=(16X-1)/(4X+1).                                      (6)
```

This sends

```text
X=0,1,infinity,-3/2       to       x=-1,3,4,5.          (7)
```

Write the resulting maps as `g_a,g_b,g_c`, and let `E` be the
normalization of

```text
g_a(a)=g_b(b)=g_c(c)=T.                                (8)
```

As in the first mixed example, the inertia sign vectors at the four
branch values are

```text
(1,1,0),       (1,0,1),       (0,1,0),       (0,0,1).
```

They span `(Z/2Z)^3`.  Riemann--Hurwitz therefore makes `E` a connected
genus-one curve.  Every fiber over a common value in (3) splits into
eight rational points, so `E` has rational points.  The individual
`a`-involution is free, while the `b`- and `c`-involutions are
reflections.

## 3. All boundary supports are distinct

At the common target corresponding to `X=0`, every nonempty subset of
`{a,b,c}` may choose the root zero.  This gives the seven boundaries

```text
{0} union S,       empty != S subset {a,b,c}.
```

The distinct common targets at `X=1` and `X=infinity` similarly give
seven boundaries each.  At the fourth common target, `X=-3/2` is a
nonanchor, so precisely the subsets of size two or three collide,
giving four boundaries.  Hence

```text
7+7+7+4=25.                                             (9)
```

The four target fibers are distinct, and within each fiber the choice
of inverse root determines the source point.  The twenty-five contacts
in (9) therefore lie at twenty-five distinct points of `E`.

There are no further collisions.  Equation (2) is the exact
cross-product of every pair of pencil members, so moving coordinates
can agree only at the four common arguments.  A moving coordinate can
meet an anchor only at the corresponding three common target values.

Finally, the four roots in (2) are simple, every denominator in (1) is
nonzero there for the parameters (4), and the common values are not
branch values.  Thus each inverse derivative is nonzero, and the
cross-product's simple zeros make the three inverse derivatives
pairwise distinct.  Stable reduction produces one first-level bubble
with distinct marked points at every listed cluster.  All twenty-five
boundary contact orders are exactly one.

This example removes the geometric shared-support obstruction. Its
positive rational rank, proved below, also makes the separate
[integer contact-profile lemma](disjoint_mixed_elliptic_integer_contact_profile.md)
unconditional: an infinite rational orbit subsequence realizes all
twenty-five shared unoriented integer contact cores, with pairwise
coprime supports, equal leading logarithmic sizes, and bounded primitive
row determinant remainders. The Gaussian split-prime restriction and
orientations, six private singleton factors, and first short
determinant-one frame remain outside that realization. It supplies
neither a short-arc counterexample nor the global uniform lattice-arc bound.

Run
[check_disjoint_mixed_elliptic_six_curve.py](check_disjoint_mixed_elliptic_six_curve.py)
for exact rational verification.

## 4. Exact quotient quartic and positive rank

The inverse equation for the normalized coordinate (X) is a binary
quadratic (N_lambda(X)-T D_lambda(X)=0).  Its exact discriminants in
the target coordinate (T), in the order (a,b,c), are

```text
Delta_a(T) = 400 T(T-4),

Delta_b(T) = (400/73^4) T(67144321 T-212248804),

Delta_c(T) = (400/257^2)(T-4)(46513 T-208624).
```

The three sign changes give the multiquadratic cover.  The kernel of the
character `chi_b chi_c` is the full translation subgroup `E[2]`; hence
the quotient by translations has the quartic model

```text
y^2 = F(T)
F(T) = T(67144321 T-212248804)(T-4)(46513 T-208624).       (10)
```

Indeed, `Delta_b Delta_c` is the square of
`400/(73^2*257)` times `F`.  This square extraction is allowed over
`Q`, but the leading coefficient of `F` is
`67144321*46513 = 3123083802673`, which is not a rational square.  Thus
the quartic must be counted with its quadratic twist intact; replacing
`F` by the monic product of its four rational root factors changes the
curve.

For an odd prime of good reduction, the checker counts the affine points
by evaluating the displayed integer polynomial `F` modulo that prime.
It adds two points at infinity exactly when the leading coefficient is a
square modulo the prime.  At (p=23), the four roots are
`0,16,4,2`, the leading coefficient is a square, and the count is
`26+2=28`.  At (p=37), the roots are `0,11,4,23`, the leading
coefficient is a nonsquare, and the count is `40+0=40`.  The roots are
distinct in both fields, so both reductions are good.  Therefore

```text
|E(Q)_tors| = |(E/E[2])(Q)_tors| divides gcd(28,40)=4.              (11)
```

The quartic also gives a small exact witness:

```text
P = (T,Y) = (-1/2, 45299439/4),       Y^2 = F(T).
```

The four rational roots of `F` are the full rational 2-torsion after
choosing `(0,0)` as origin.  Thus the bound (11) makes the rational
torsion exactly those four points, while `P` is not a branch point and
therefore is nontorsion.

The four common target fibers are unramified and split completely: each
of the three discriminants is a nonzero rational square at each of
`-1/2, 9/2, 16/3, 25/4`.  They consequently supply `4*8=32` distinct
rational points on `E`.  Since the quotient by `E[2]` is the curve
`E` again via multiplication by two, (11) rules out rank zero.  The
source therefore has positive rational rank, with no appeal to Mazur's
torsion classification or to floating-point rank computations.

For concrete rational six-point instances from the odd multiples `3R0`
and `5R0`, see the [explicit configuration note](disjoint_mixed_elliptic_explicit_configs.md)
and its exact checker.
