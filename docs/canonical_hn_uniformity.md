# Canonical Harder--Narasimhan flats and uniformity in the place set

## Purpose

The generalized-GCD route reduces the square-root arc theorem to a uniform
exceptional-degree statement. The finite Ru--Vojta level is already fixed in
`GaussianChain/FiniteBeta.lean`. The remaining issue is that, at each place,
the finite-level proof chooses one adapted basis from a fixed finite family.
A naive application of a quantitative Subspace Theorem to every placewise
choice separately can introduce dependence on the number of places.

This note isolates a mechanism that removes that dependence for the canonical
Harder--Narasimhan exceptional subspace.

The argument below is self-contained linear algebra. It does not prove the
parametric Subspace Theorem; instead, it proves that once that theorem assigns
a canonical maximal-slope subspace to a weighted filtration system, the
subspace belongs to one fixed finite list independent of the number of places.

---

## 1. A finite family of weighted filtrations

Let `V` be an `n`-dimensional vector space over a characteristic-zero field.
Fix a finite collection of complete flags

```text
0 = F_{a,0} subset F_{a,1} subset ... subset F_{a,n} = V,
```

indexed by `a` in a finite set `A`.

A weighted filtration of type `a` consists of decreasing real weights

```text
gamma_{a,1} >= gamma_{a,2} >= ... >= gamma_{a,n}.
```

For a subspace `U <= V`, define its local filtered degree by

```text
deg_{a,gamma}(U)
  = sum_{j=1}^n gamma_{a,j}
      (dim(U cap F_{a,j}) - dim(U cap F_{a,j-1})).
```

Telescoping gives

```text
deg_{a,gamma}(U)
  = gamma_{a,n} dim U
    + sum_{j=1}^{n-1}
        (gamma_{a,j} - gamma_{a,j+1}) dim(U cap F_{a,j}).       (1.1)
```

All coefficients in the second sum are nonnegative.

For a finite set of places, choose at each place one flag type and one weight
vector. Summing the local degrees produces a global degree of the form

```text
deg(U) = c dim U + sum_{F in mathcal F} a_F dim(U cap F),       (1.2)
```

where

- `mathcal F` is the finite set of every proper flag step occurring in the
  fixed finite family of flags;
- every `a_F` is nonnegative;
- the coefficients aggregate contributions over all places.

The number of places changes the coefficients `a_F`, but it does not change
`mathcal F`.

---

## 2. Supermodularity

For fixed subspaces `U,W,F <= V`,

```text
(U cap F) + (W cap F) subset (U + W) cap F,
```

and

```text
(U cap F) cap (W cap F) = (U cap W) cap F.
```

Therefore

```text
dim(U cap F) + dim(W cap F)
  <= dim((U cap W) cap F) + dim((U + W) cap F).          (2.1)
```

Using (1.2), while dimension itself is modular, we obtain

```text
deg(U) + deg(W)
  <= deg(U cap W) + deg(U + W).                          (2.2)
```

Thus the filtered degree is supermodular.

Let

```text
mu = max_{0 != U <= V} deg(U) / dim U.
```

If `U` and `W` attain the maximal slope, then (2.2), the modular identity

```text
dim(U cap W) + dim(U + W) = dim U + dim W,
```

and the definition of `mu` force equality throughout. Consequently `U+W`
attains maximal slope, and `U cap W` does as well when it is nonzero.

This closure statement is formalized abstractly in
`GaussianChain/SlopeClosure.lean`.

The sum of all maximal-slope subspaces is therefore maximal slope. It is the
unique largest maximal-slope subspace; call it `T(deg)`.

---

## 3. Only finitely many canonical subspaces can occur

Define the incidence signature of a subspace `U` by

```text
sigma(U) = (dim U, (dim(U cap F))_{F in mathcal F}).       (3.1)
```

There are at most

```text
(n+1)^(|mathcal F|+1)
```

possible signatures.

For each realizable signature `sigma`, let

```text
C_sigma = sum { U <= V : sigma(U) = sigma }.              (3.2)
```

The sum is well-defined: finite dimensionality means it is already generated
by finitely many members of the class.

For every choice of the nonnegative coefficients in (1.2), the value of
`deg(U)` depends only on `sigma(U)`. Hence, if one subspace of signature
`sigma` has maximal slope, every subspace of that signature has maximal slope.
By closure under sums, `C_sigma` then has maximal slope.

Let `Sigma_max(deg)` be the set of signatures attaining maximal slope. The
unique largest maximal-slope subspace is

```text
T(deg) = sum_{sigma in Sigma_max(deg)} C_sigma.           (3.3)
```

The right side ranges over a finite collection of subspaces depending only on
the original finite flag family:

```text
# { T(deg) } <= 2^(# signatures).                         (3.4)
```

It is independent of

- the number of places;
- which adapted flag is selected at each place;
- the numerical local weights.

The finite-signature counting mechanism is formalized abstractly in
`GaussianChain/FiniteFlatReduction.lean`.

---

## 4. Application to the finite Ru--Vojta construction

Fix the explicit level

```text
N = 100,  d = 347,
```

from `FiniteBeta.lean`. The section space

```text
mathcal V = H^0(Bl_P P^2, 347 H - 100 E)
```

is fixed and finite dimensional.

The finite Autissier--Ru--Vojta argument uses only finitely many adapted bases
of `mathcal V`; after ordering each basis by its local weight, these bases give
only finitely many flags. Therefore Sections 1--3 apply to every placewise
choice made in the proof.

For every resulting parametric Subspace-Theorem system, its canonical
Harder--Narasimhan exceptional subspace belongs to one fixed finite list

```text
T_1, ..., T_J subset mathcal V^*,                         (4.1)
```

where `J` is independent of the finite place set `S`.

Pulling the `T_j` back through the fixed finite-level projective map and then
pushing down to `P^2` gives a finite list of plane curves of uniformly bounded
total degree.

This conclusion applies to the canonical maximal-slope subspace. It does not,
by itself, show that every auxiliary subspace appearing in a quantitative
covering theorem belongs to the same list.

---

## 5. Remaining imported theorem

To finish the uniform exceptional-degree theorem, one needs the following
uniform asymptotic statement for the parametric Subspace Theorem.

> For a fixed finite collection of algebraic linear forms and a fixed positive
> approximation margin, there is a parameter threshold independent of the
> finite place set and of the normalized local weights, such that every
> sufficiently large solution lies in the canonical maximal-slope subspace.

The absolute parametric theorem is naturally normalized by

```text
sum_v max_i c_{iv} <= 1,
```

and its large-parameter threshold is stated in terms of the ambient dimension,
the approximation margin, the height of the fixed form family, and the number
of distinct forms. The remaining line-by-line task is to verify that the
finite Ru--Vojta max-over-adapted-bases inequality yields one normalized
parametric system to which this canonical-subspace conclusion applies, without
choosing a different placewise system separately for every point.

Once that matching is established, Sections 1--4 make the exceptional list and
its total pullback degree uniform in `S`.

---

## 6. Status

Proved directly here:

1. filtered degrees have the form (1.2);
2. they are supermodular;
3. maximal-slope subspaces are closed under sums and nonzero intersections;
4. the canonical largest maximal-slope subspace ranges over a fixed finite
   list independent of the place set.

Still to be supplied:

1. the exact matching between the finite Ru--Vojta local maximum and one
   normalized twisted-height system;
2. the uniform canonical-subspace conclusion for that system;
3. pullback-degree bookkeeping for the fixed finite list.

The dependence on the number of places has therefore been removed for the
canonical Harder--Narasimhan object. The unresolved issue is no longer a raw
component count; it is the precise passage from the local maximum to the
canonical parametric system.
