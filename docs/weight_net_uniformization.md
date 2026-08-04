# Weight-net uniformization for the conductor-relative route

## 1. Purpose

The conductor-relative Ru--Vojta strategy has two logically separate problems.

1. For one fixed normalized weight system, prove a strict generalized-GCD or
twisted-height inequality outside a bounded-degree exceptional set.
2. Make the exceptional set and height threshold uniform as the conductor and
its local weights vary.

The second problem is not intrinsically difficult once the conductor dependence
has been compressed into one fixed finite-dimensional weight simplex.  A strict
coefficient margin and a finite epsilon-net give uniformity.

The actual load-bearing issue is therefore the compression step.

---

## 2. Abstract finite-type score

Let `T` be a fixed finite set of local filtration types.  For a normalized
conductor let

```text
omega_t >= 0,              sum_t omega_t = 1.
```

Let `s_t(P)` be the normalized contribution of type `t` at a point `P`, and
assume

```text
|s_t(P)| <= H(P)
```

for every `t`.  Define

```text
S_omega(P) = sum_t omega_t s_t(P).
```

Then

```text
|S_omega(P)-S_omega'(P)|
  <= H(P) sum_t |omega_t-omega'_t|.                     (2.1)
```

This is formalized in `GaussianChain/WeightStability.lean`.

---

## 3. Stability of a strict coefficient

Suppose a fixed-weight theorem gives, outside an exceptional set `Z_omega0`,

```text
S_omega0(P) < eta0 H(P)
```

with

```text
eta0 < eta < 1/2.
```

If

```text
sum_t |omega_t-omega0_t| <= eta-eta0,
```

then (2.1) gives

```text
S_omega(P) < eta H(P)
```

outside the same exceptional set.

Thus one fixed-weight exceptional set controls an entire L1 neighborhood of its
weight vector.  The strict numerical margin is essential.

At the minimal finite Ru--Vojta level `(N,d)=(1,6)`, the zero-error coefficient
is

```text
6 - 3*(50/27) = 4/9.
```

For example one may choose

```text
eta0 = 4/9 + epsilon,
eta  < 1/2,
```

leaving a positive net radius.

---

## 4. Finite-net consequence

The probability simplex on the finite set `T` is compact.  Hence it has a
finite L1 net of any prescribed positive radius.

Apply the fixed-weight theorem at each net center.  Taking the union of the
finitely many exceptional sets and the maximum of their finitely many height
thresholds gives one exceptional family and one threshold valid for every
weight vector.

No quantitative theorem uniform in a continuously varying weight vector is
needed.  Pointwise fixed-weight theorems plus the strict margin are enough.

The abstract finite-net implication is formalized in
`GaussianChain/WeightStability.lean`; only the elementary existence of a finite
net for the finite-dimensional simplex remains topological bookkeeping.

---

## 5. Why this does not yet finish the proof

For the conductor problem, aggregation of the **slope data** is already proved:

```text
sum_v w_v c(type(v),W)
  = sum_t omega_t c(t,W).
```

This is `GaussianChain/ConductorWeights.lean`.

But the local Ru--Vojta expression contains evaluations of sections:

```text
sum_v max_B sum_{s in B} lambda_{s,v}(P).
```

The maximizing basis and the local section values depend on both `v` and `P`.
The slope aggregation alone does not prove an identity

```text
local maximum = sum_t omega_t s_t(P)
```

with a fixed finite collection of global functions `s_t`.

This is the exact missing bridge.

---

## 6. Correct next theorem

A sufficient completion is the following conductor factorization theorem.

> **Finite-type conductor factorization.**  For the fixed 27-dimensional
> section space `H^0(Bl_P P^2,6H-E)`, there are a fixed finite type set `T` and
> fixed normalized score functions `s_t(P)` such that every ordered pair from
> one Gaussian conductor has a probability vector `omega` satisfying
>
> ```text
> Ru--Vojta global max(P) <= sum_t omega_t s_t(P) + O(1),
> ```
>
> with `|s_t(P)| <= H_G`, and with a total error smaller than the available
> margin `1/2-4/9=1/18`.

If this theorem holds, the finite-net argument above upgrades any fixed-weight
parametric Subspace Theorem to a uniform conductor-relative exceptional family.
The existing grid argument then proves the square-root arc bound.

---

## 7. Current logical boundary

Proved and formalized:

1. explicit finite beta coefficient `4/9 < 1/2`;
2. aggregation of linear slope scores by finite type;
3. L1 stability of weighted scores;
4. finite-net uniformization conditional on score factorization;
5. the final off-diagonal grid bound.

Not proved:

1. finite-type factorization of the actual max-over-adapted-bases evaluation
   score;
2. an alternative fixed-family max-type Subspace Theorem that accepts that
   score directly.

Therefore the remaining difficulty is algebraic, not compactness or threshold
uniformity.
