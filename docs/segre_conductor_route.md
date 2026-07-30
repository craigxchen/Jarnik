# Segre common-conductor reduction

## Status

This note gives a fixed geometric compactification for the common-conductor
version of the generalized-GCD route.  It repairs a defect in the naive
four-coordinate construction: the conductor line in `P^3` passes through a
coordinate point and therefore violates the coordinate-position hypothesis in
Yasufuku's codimension-two theorem.

The Segre construction keeps every ambient coordinate monomial in the
`S`-units, makes the moving conductor an honest projective coordinate, and
leaves the point-blowup beta invariant unchanged.  Its numerical coefficient is
strictly below `1/2`.

It does not, by itself, eliminate the dependence of the complete exceptional
set on the finite place set.  The exact remaining theorem is stated in Section
8.

---

## 1. Common-conductor coordinates

Fix a Gaussian conductor `G`, and put

```text
q = G / conjugate(G).
```

For two norm-one divisor ratios from the same conductor, write

```text
x = A / conjugate(A),
y = B / conjugate(B),
```

with the prime exponents of `A,B` lying in the conductor box of `G`.

The ordered pair is represented by

```text
(([1:q],[1:x:y])) in P^1 x P^2.
```

Under the Segre embedding this becomes

```text
[1 : q : x : qx : y : qy] in P^5.                 (1.1)
```

Every displayed coordinate is an `S`-unit for the common support of `G`.

The endpoint center is

```text
Y = P^1 x {[1:1:1]}.
```

Inside the Segre threefold this has codimension two.  Its Segre image is the
line

```text
[1:q:1:q:1:q].                                      (1.2)
```

Unlike the naive line `[1:1:1:q]` in `P^3`, its projective closure contains no
coordinate point of `P^5`: at each endpoint of the `q`-line, three coordinates
remain nonzero.

---

## 2. The blowup is a product

Let

```text
P = [1:1:1] in P^2,
X = P^1 x P^2.
```

Blowing up the product center gives

```text
Bl_{P^1 x {P}}(P^1 x P^2)
    isomorphic to P^1 x Bl_P(P^2).                    (2.1)
```

This follows because the ideal sheaf of the center is pulled back from the
second factor and the Rees algebra commutes with flat base change.

Write

```text
H_q = pullback of O_{P^1}(1),
H   = pullback of a line on P^2,
E   = exceptional divisor.
```

The generalized GCD of `x-1` and `y-1` is the height of `E`.

---

## 3. A small conductor coefficient

Use the rational divisor class

```text
L_tau = tau H_q + 6 H - E,                            (3.1)
```

where `tau>0` is rational.  It is ample because `tau H_q` is ample on the first
factor and `6H-E` is ample on the point blowup.

We take

```text
tau = 1/36.                                           (3.2)
```

An integral representative is

```text
36 L_tau = H_q + 216 H - 36 E.
```

---

## 4. Tensor invariance of the beta calculation

Let `D_i=P^1 x H_i`, where `H_0,H_1,H_2` are the three coordinate lines in
`P^2`.  The divisors `D_i` are pulled back from the second factor.

For every sufficiently divisible integer `N`, Kunneth gives

```text
H^0(N L_tau)
  = H^0(P^1,O(N tau)) tensor
    H^0(Bl_P P^2, N(6H-E)),                            (4.1)
```

and

```text
H^0(N L_tau - m D_i)
  = H^0(P^1,O(N tau)) tensor
    H^0(Bl_P P^2, N(6H-E)-mH_i).                       (4.2)
```

Thus both the finite beta numerator and denominator are multiplied by the same
factor

```text
dim H^0(P^1,O(N tau)).
```

It cancels.  Therefore

```text
beta(L_tau,D_i) = beta(6H-E,H_i).                       (4.3)
```

At the explicit finite level `(N,d)=(1,6)`, the second-factor calculation is

```text
dim H^0(6H-E) = 27,
beta = 50/27,
6 - 3 beta = 4/9.                                      (4.4)
```

These numerical identities are formalized in
`GaussianChain/FiniteBeta.lean`.

---

## 5. Height comparison

For a common-conductor point

```text
Q = ([1:q],[1:x:y]),
```

heights for the product divisor satisfy

```text
h_{L_tau}(Q)
  = tau h(q) + 6 h([1:x:y]) - h_E(Q) + O(1).           (5.1)
```

For ratios from one circle,

```text
h(q) <= W,
h([1:x:y]) <= W,                                      (5.2)
```

where `W=log R` in the circle-radius normalization.

The endpoint arc condition gives

```text
h_E(Q) >= W/2 - O_C(1).                                (5.3)
```

---

## 6. The coefficient supplied by Ru--Vojta

Applying the finite-level Ru--Vojta inequality to the three divisors `D_i`
gives, outside its exceptional set,

```text
3 beta h([1:x:y])
  <= (1+epsilon) h_{L_tau}(Q) + outside-S terms.        (6.1)
```

For the common support `S(G)`, all six Segre coordinates in (1.1) are
`S(G)`-units, so the coordinate outside-`S` contribution vanishes.

Substituting (5.1) and rearranging yields

```text
h_E(Q)
  <= tau h(q) + (6-3 beta) h([1:x:y])
       + O(epsilon W) + O(1).                           (6.2)
```

Using (5.2), (4.4), and `tau=1/36`, the main coefficient is

```text
1/36 + 4/9 = 17/36 < 1/2.                              (6.3)
```

The remaining margin is exactly

```text
1/2 - 17/36 = 1/36.                                    (6.4)
```

These rational identities are formalized in
`GaussianChain/SegreCoefficient.lean`.

Choosing the Ru--Vojta epsilon small enough preserves a strict coefficient
below `1/2`.  Equations (5.3) and (6.2) are then incompatible for sufficiently
large `R`, unless `Q` is exceptional.

---

## 7. Why this is a real improvement

The naive conductor point

```text
[1:x:y:q] in P^3
```

has center

```text
{x=1,y=1},
```

whose projective closure contains the coordinate point `[0:0:0:1]`.  This is
exactly the type of coordinate incidence excluded by the codimension-two
blowup hypotheses: one exceptional component would interact with several
coordinate pullbacks.

The Segre compactification resolves that problem without replacing the
`S`-unit coordinates by additive coordinates.  It also isolates the conductor
height in a line-bundle direction whose coefficient can be made arbitrarily
small while the point-blowup beta invariant remains unchanged.

---

## 8. Exact remaining theorem

For a fixed finite set of places, Ru--Vojta gives a proper exceptional closed
set.  The place set in the circle problem is `S(G)` and moves with the
conductor.  The unrestricted uniform-in-`S` exceptional-degree theorem is
false, as shown by the dense family of unrelated near-one `S`-unit pairs.

The Segre family has an additional coordinate `q` which measures the whole
common conductor.  The counterexample family does not remain critical after
that height is included.

The exact theorem now sufficient for the endpoint problem is:

> **Segre one-scale theorem.**  For the fixed threefold
> `P^1 x Bl_P(P^2)`, the fixed divisors `D_0,D_1,D_2`, and the fixed rational
> line bundle `L_{1/36}`, there are constants `T,Delta,H_0`, independent of the
> Gaussian conductor `G`, such that the conductor-critical points
> `([1:q],[1:x:y])` at the single scale `h(q)` lie in at most `T` proper
> subvarieties of degree at most `Delta` whenever `h(q)>H_0`.

The absolute parametric Subspace Theorem has uniform numerical bounds for one
fixed form-and-weight system.  To prove the displayed theorem one must construct
from `q` one point-independent normalized weight system that covers every
ordered divisor pair from its conductor box.  This is narrower than the false
uniform theorem because the weights are tied to the actual conductor
coordinate.

Once the Segre one-scale theorem is established, intersecting its exceptional
subvarieties with the fiber over `q` gives a uniformly bounded-degree family of
curves in the `(x,y)`-plane.  The existing finite-family grid lemma then bounds
the endpoint cluster uniformly.
