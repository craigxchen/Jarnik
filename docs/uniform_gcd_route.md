# A uniform generalized-GCD route to the square-root arc problem

## Status

This note separates four logically distinct statements.

1. **Circle-side reduction.** A uniform bounded-degree exceptional-curve theorem for the two-variable generalized-GCD inequality implies a uniform bound for lattice points on arcs of length `C * sqrt R`.
2. **Finite-level Ru--Vojta calculation.** On the blowup of `P^2` at `[1:1:1]`, the finite section counts are explicit. The asymptotic beta threshold is already crossed at `(N,d) = (100,347)`.
3. **Canonical-flat uniformity.** At fixed finite level, the local adapted bases produce a fixed finite family of weighted filtrations. Their global filtered degree is supermodular, so maximal-slope subspaces are closed under sums. The unique largest maximal-slope subspace ranges over a fixed finite list independent of the place set. This is developed in `canonical_hn_uniformity.md` and partially formalized in Lean.
4. **Remaining imported input.** One still needs the large-parameter conclusion of the absolute parametric Subspace Theorem uniformly over the normalized local weights: above one threshold depending only on the fixed form family and the approximation margin, every solution lies in the canonical maximal-slope subspace.

This document is therefore a rigorous conditional proof and gap audit. It does not insert the remaining parametric theorem as an axiom.

---

## 1. The uniform exceptional-curve hypothesis

Let `K = Q(i)`. Fix a real number `eta < 1/2`.

### Hypothesis UEC(eta)

There exist constants `D, H0` depending only on `K` and `eta` such that for every finite set of places

```text
S superset M_K^infty
```

there is a nonzero polynomial

```text
F_S(X,Y) in K[X,Y],       deg F_S <= D,
```

with the following property. For all `S`-units `x,y in K^*` with

```text
h([1:x:y]) > H0
```

and `F_S(x,y) != 0`,

```text
log GCD^+(x-1,y-1) < eta * h([1:x:y]).                 (UEC)
```

The polynomial may depend on `S`; its total degree may not.

---

## 2. Conditional uniform bound for square-root arcs

### Theorem 2.1

Assume `UEC(eta)` for some `eta < 1/2`. Then for every `C > 0` there is a constant `B(C)` such that every arc of length at most

```text
C * sqrt R
```

on a centered circle of radius `R`, with `R^2 in N`, contains at most `B(C)` points of `Z[i]`.

### Proof

Let

```text
z_0, z_1, ..., z_{M-1} in Z[i]
```

be distinct, have modulus `R`, and lie on one such arc. Put

```text
W   = log R,
u_j = z_j / z_0        (1 <= j < M).
```

#### Step 1: one common `S`

Let `S` consist of the archimedean place of `Q(i)` and all finite places lying over rational primes dividing `R^2`. Since every `z_j` has norm `R^2`, every quotient `u_j` is an `S`-unit. The same `S` works for the entire cluster.

#### Step 2: archimedean proximity

Chord length is at most arc length, hence

```text
|u_j - 1| = |z_j-z_0| / R <= C * R^(-1/2).            (2.1)
```

For distinct `i,j`, the complex-place contribution to the exceptional-divisor local height at `[1:1:1]` is therefore at least

```text
W/2 - log C.                                           (2.2)
```

At every nonarchimedean place,

```text
|x-1|_v <= max(1,|x|_v),
```

so the corresponding local contribution to `GCD^+` is nonnegative. Consequently

```text
log GCD^+(u_i-1,u_j-1) >= W/2 - log C.                 (2.3)
```

#### Step 3: joint-height upper bound

We claim

```text
h([1:u_i:u_j]) <= W.                                   (2.4)
```

At a split rational prime `p = pi * conjugate(pi)`, write

```text
E = v_p(R^2),
a_k = v_pi(z_k),
v_conj(pi)(z_k) = E-a_k.
```

The combined normalized contribution of the two conjugate places to `h([1:u_i:u_j])` is

```text
(1/2) * (max{a_0,a_i,a_j} - min{a_0,a_i,a_j}) * log p,
```

which is at most `(E/2) log p`. Inert primes contribute nothing to the ratios, and the archimedean contribution is zero because all `u_j` have modulus one. Summing gives (2.4).

#### Step 4: every pair is exceptional for large `R`

Outside `F_S = 0`, `UEC(eta)` and (2.4) give

```text
log GCD^+(u_i-1,u_j-1) < eta W.
```

Together with (2.3),

```text
W/2 - log C < eta W.
```

This is impossible once

```text
log C < (1/2 - eta) W.
```

Hence all ordered off-diagonal pairs `(u_i,u_j)` lie on the same curve `F_S = 0` for sufficiently large `R`.

#### Step 5: grid bound

Let `A = {u_1,...,u_{M-1}}`. A nonzero polynomial of total degree at most `D` has at most `D |A|` zeros on `A x A`. Since all ordered off-diagonal pairs are zeros,

```text
|A|(|A|-1) <= D|A|,
```

so `|A| <= D+1` and the full cluster has at most `D+2` points. The bounded range of smaller radii is absorbed into `B(C)`.

The numerical gap and abstract grid/fiber argument are formalized in `GaussianChain/UniformExceptionalReduction.lean`.

---

## 3. Explicit finite beta level

Let

```text
V = Bl_[1:1:1] P^2,
```

with `H` the pullback of a line and `E` the exceptional divisor. For integer `N > 0` and excess degree `k >= 0`, put `d = N+k`. Then

```text
h^0(V, dH-NE)
  = binom(d+2,2) - binom(N+1,2).
```

The finite beta numerator is

```text
sum_{j=0}^{k-1}
  [ binom(N+j+2,2) - binom(N+1,2) ]
  = k(k+1)(3N+k+2)/6.                                  (3.1)
```

At

```text
N = 100,  k = 247,  d = 347,
```

Lean verifies

```text
h^0(V,347H-100E) = 55676,
beta numerator   = 5604924,
beta              = 45201/44900 > 1.                  (3.2)
```

Thus the section space and the finite family of adapted bases can be fixed independently of `S`.

---

## 4. Canonical maximal-slope subspace

At fixed finite level, order each adapted basis by its local weights. This gives one of finitely many weighted flags in the fixed section space. Summing over the places gives a filtered degree

```text
deg(U) = c dim U + sum_{F in mathcal F} a_F dim(U cap F),
```

where `mathcal F` is fixed and every `a_F >= 0`.

For every fixed `F`, the function `U -> dim(U cap F)` is supermodular. Hence

```text
deg(U) + deg(W)
  <= deg(U cap W) + deg(U + W).                         (4.1)
```

If `mu` is the maximal slope and `U,W` attain it, modularity of dimension and (4.1) force both `U cap W` and `U+W` to attain it. Therefore there is a unique largest maximal-slope subspace.

Its value can only come from a fixed finite list. Indeed, `deg(U)` depends only on the finite incidence signature

```text
(dim U, (dim(U cap F))_{F in mathcal F}).               (4.2)
```

For each signature, take the sum of all subspaces having that signature. If the signature is maximal-slope, this sum is maximal-slope as well. The canonical largest subspace is therefore the sum of a subset of finitely many fixed signature closures.

The abstract slope-closure identity is formalized in `GaussianChain/SlopeClosure.lean`; finite-signature range bounds are formalized in `GaussianChain/FiniteFlatReduction.lean`. A complete proof and the exact scope of the reduction are in `docs/canonical_hn_uniformity.md`.

This removes the naive component-count dependence on `|S|` for the canonical Harder--Narasimhan subspace.

---

## 5. Remaining theorem

The remaining input is now the following uniform asymptotic statement.

> For a fixed finite collection of algebraic linear forms and a fixed positive approximation margin, there is a height/parameter threshold independent of the finite place set and of the normalized local weights, such that every sufficiently large solution of the twisted-height inequality lies in the canonical maximal-slope subspace.

The absolute parametric Subspace Theorem is formulated with a threshold depending on the ambient dimension, the approximation margin, the height of the fixed form family, and the number of distinct forms—not on the number of places. What remains to check line by line is that the finite Ru--Vojta max-over-adapted-bases inequality can be assigned one such parametric system without introducing a point-dependent family outside the fixed filtration framework.

Once that matching is complete, the canonical-flat theorem gives a finite exceptional list independent of `S`; its pullbacks have uniformly bounded total degree; and Section 2 proves the square-root arc theorem.

---

## 6. Formalization boundary

Compiler-checked or represented directly in Lean:

- exact finite beta numerator;
- explicit beta value at `(100,247)`;
- the `1/2 - eta` numerical contradiction;
- the abstract grid/fiber cardinality bound;
- maximal-slope closure under meet and join;
- finite-signature range bounds.

Not yet formalized in mathlib:

- number-field Weil heights and local Weil functions in the needed form;
- generalized `GCD^+`;
- blowups and Ru--Vojta filtrations;
- the absolute parametric Subspace Theorem;
- the exact matching of the finite adapted-basis maximum to one normalized twisted-height system.
