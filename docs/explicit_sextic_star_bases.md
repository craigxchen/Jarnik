# Explicit adapted bases at the minimal sextic level

## 1. The section space

Let

```text
P = [1:1:1] in P^2,
L = 6H-E on Bl_P P^2.
```

The section space is

```text
V = H^0(P^2,I_P(6)).
```

Let `M_6` be the 28 degree-six monomials

```text
X_0^a X_1^b X_2^c,       a+b+c=6.
```

Evaluation at `P` sends every monomial to `1`.  Hence

```text
V = { sum_m c_m m : sum_m c_m = 0 },
```

and `dim V = 27`.

Fix a monomial `m_*`.  Then

```text
B(m_*) = { m-m_* : m in M_6, m != m_* }
```

is a basis of `V`.

---

## 2. The weighted monomial filtration

At one nonarchimedean place, normalize homogeneous coordinates so that

```text
ell_i >= 0,
min(ell_0,ell_1,ell_2)=0,
```

where `ell_i` is the local coordinate height of `X_i`.

Give the monomial `X_0^a X_1^b X_2^c` weight

```text
w(a,b,c)=a ell_0+b ell_1+c ell_2.
```

For a real threshold `t`, the fixed-level Autissier filtration step is

```text
F^t V
 = V cap span {m in M_6 : w(m) >= t}.
```

Indeed, the sum of coordinate ideals with weighted order at least `t` is the
monomial ideal spanned in degree six by precisely those monomials.

Order the monomials

```text
m_1,m_2,...,m_28
```

so that

```text
w(m_1) >= w(m_2) >= ... >= w(m_28).
```

Then

```text
m_1-m_2, m_1-m_3, ..., m_1-m_28
```

is adapted to the filtration.  When the first `r` monomials have entered, their
intersection with `V` has dimension `r-1` and is spanned by the first `r-1`
differences.

Thus one may take the star basis `B(m_1)` rooted at a maximal-weight monomial.

---

## 3. Only three root types occur

Let `i` satisfy

```text
ell_i = max(ell_0,ell_1,ell_2).
```

For every degree-six monomial,

```text
w(a,b,c)
  <= (a+b+c) ell_i
  = 6 ell_i.
```

The pure power `X_i^6` has weight exactly `6 ell_i`.  Hence it is a valid
maximal root.  Ties may be resolved by a fixed rule.

Consequently the complete finite-level construction needs only the three bases

```text
B_0 = B(X_0^6),
B_1 = B(X_1^6),
B_2 = B(X_2^6).
```

The choice of basis is not arbitrary: it is determined by the tropical chamber
in which the coordinate-height vector lies.

---

## 4. Exact filtration weight

Across all 28 degree-six monomials, symmetry gives

```text
sum_m exponent_i(m)=56
```

for each coordinate `i`, because the total sum of all three exponents is

```text
28*6=168.
```

Therefore

```text
sum_{m in M_6} w(m)=56(ell_0+ell_1+ell_2).
```

Removing the pure root `X_i^6` leaves

```text
56(ell_0+ell_1+ell_2)-6 ell_i.
```

Since `ell_i <= ell_0+ell_1+ell_2`, this is at least

```text
50(ell_0+ell_1+ell_2).                         (4.1)
```

The numerical part of (4.1) is formalized in
`GaussianChain/SexticTropicalWeight.lean`.

This is exactly the finite beta numerator

```text
50 = 27*(50/27).
```

No limiting Riemann--Roch calculation, infinite adapted-basis family, or basis
selection lemma is needed at this minimal level.

---

## 5. Evaluation-product form

Define the three fixed product sections

```text
P_i = product_{m in M_6, m != X_i^6} (m-X_i^6).
```

Each `P_i` is a section of `27(6H-E)`.  At a place whose largest coordinate
height is `ell_i`, the valuation estimate

```text
v(m-X_i^6) >= min(v(m),v(X_i^6)) = v(m)
```

and (4.1) give

```text
v(P_i(P)) >= 50(ell_0+ell_1+ell_2).
```

Thus the local Ru--Vojta lower bound is witnessed by exactly one of three fixed
product sections, selected by the largest coordinate height.

---

## 6. What this repairs

This proves, at the level actually needed for the endpoint coefficient, all of
the following.

1. The adapted bases can be made completely explicit.
2. There are only three basis types.
3. Their selection is coupled to a fixed three-chamber tropical partition.
4. The filtration weight is exactly bounded below by `50` times the coordinate
   proximity.

It removes the earlier ambiguity about infinitely many adapted bases.

---

## 7. What still remains

The arithmetic upper bound must handle the point-dependent chamber selection:

```text
sum_v lambda_{P_{i(v,P)},v}(P),
```

where `i(v,P)` is the coordinate of largest local height.

Replacing this selected expression by the unrestricted maximum over all three
`P_i` is too strong; arbitrary `S`-unit pairs give a counterexample to a
uniform exceptional-degree theorem for that maximum.

The remaining viable bridge must therefore exploit the coupling

```text
i(v,P) = argmax_i ell_i(v,P),
```

not discard it.  The explicit star-basis calculation reduces that bridge to a
fixed, three-chamber problem.
