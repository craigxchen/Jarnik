# Finite combinatorial types for fixed-level Autissier filtrations

## Status

This note proves the finite-type statement that was previously used without adequate justification.

At a fixed finite approximation level, the Autissier filtrations arising from a fixed finite collection of ideals have only finitely many distinct chains of subspaces. Consequently one may preselect finitely many adapted bases, including simultaneous adapted bases for any fixed finite number of filtrations. The pointwise maximum occurring in the finite-level Ru--Vojta argument may therefore be taken over a fixed finite family of bases.

This statement does **not** prove the remaining uniform parametric Subspace Theorem threshold. It repairs only the finite-family step.

---

## 1. Abstract setup

Let `K` be a field and let `V` be a finite-dimensional `K`-vector space. Fix a finite index set

```text
B subset N^r
```

and for every `b in B` fix a subspace

```text
V_b subset V.
```

For a weight vector

```text
t = (t_1,...,t_r) in Delta_{r-1}
```

and a real number `x >= 0`, define

```text
F_t(x) = sum_{b in B, <t,b> >= x} V_b.                 (1.1)
```

Here `<t,b> = sum_i t_i b_i`.

The family `x -> F_t(x)` is a decreasing filtration of `V`. We call its ordered chain of distinct subspaces its **filtration type**.

---

## 2. Finite-type theorem

### Theorem 2.1

As `t` varies over the simplex `Delta_{r-1}`, the filtrations `F_t` in (1.1) have only finitely many filtration types.

More precisely, the filtration type depends only on the weak ordering of the finite set of real numbers

```text
{ <t,b> : b in B }.
```

Therefore the number of possible types is at most the number of weak orderings of `B`, and in particular is finite.

### Proof

Fix `t`. The filtration can change only when `x` crosses one of the finitely many values `<t,b>`. If the distinct values are

```text
alpha_1 > alpha_2 > ... > alpha_m,
```

then the distinct filtration steps are among

```text
sum_{b : <t,b> >= alpha_j} V_b,     1 <= j <= m.       (2.1)
```

Thus the entire chain is determined by the nested subsets

```text
{ b in B : <t,b> >= alpha_j }.
```

Those subsets are determined by the weak ordering of the numbers `<t,b>`.

There are finitely many weak orderings of a finite set. Equivalently, the simplex is cut into finitely many relatively open polyhedral cells by the hyperplanes

```text
<t,b-b'> = 0,       b,b' in B.                         (2.2)
```

On each cell all strict comparisons and equalities among the values `<t,b>` are constant, hence the filtration type is constant. This proves the theorem. `square`

---

## 3. Why the effective multiindex set is finite

The Autissier construction is initially written using all multiindices `b in N^r`. At fixed line bundle and fixed approximation level, only finitely many can contribute.

Let `X` be projective, let `L` be a fixed line bundle, and let `D_1,...,D_r` be fixed effective Cartier divisors. Put

```text
V_b = H^0(X, L(-sum_i b_i D_i)) subset H^0(X,L).        (3.1)
```

Assume that for each `i` there is a finite integer `M_i` such that

```text
H^0(X,L(-mD_i)) = 0       for all m > M_i.              (3.2)
```

Then `V_b = 0` whenever some `b_i > M_i`. Hence only the finite box

```text
B_eff = product_i {0,...,M_i}                            (3.3)
```

contributes to any sum of the form (1.1).

Condition (3.2) is automatic in the fixed-level application. For example, after choosing an ample divisor `A`, a nonzero section of `L(-mD_i)` would imply

```text
(L-mD_i).A^(n-1) >= 0,
```

which bounds `m` whenever `D_i.A^(n-1) > 0`. In the explicit projective-plane or blowup calculation the bound is even more elementary: a polynomial of fixed degree cannot vanish to arbitrarily high order along a fixed nonzero ideal.

For the ideal-sheaf formulation used by Autissier, one writes

```text
I(t,x) = sum_{b : <t,b> >= x} I_1^{b_1} ... I_r^{b_r}
```

and then

```text
F_t(x) = H^0(X,L tensor I(t,x)).                         (3.4)
```

Inside the fixed finite-dimensional space `H^0(X,L)`, terms beyond the effective bounds contribute zero. Thus (3.4) reduces to (1.1) with a finite `B_eff`.

---

## 4. Finitely many adapted bases

A filtration type is a finite chain

```text
0 = W_0 subset W_1 subset ... subset W_m = V.           (4.1)
```

Choose once and for all a basis adapted to each filtration type: choose a basis of `W_1`, extend it to one of `W_2`, and continue.

Since Theorem 2.1 gives only finitely many chains, this produces a fixed finite family of bases.

The same argument handles simultaneous adaptation.

### Corollary 4.2

Fix an integer `q >= 1`. For every ordered `q`-tuple of Autissier filtrations at the fixed level, assume a simultaneous adapted basis exists. Then there is a fixed finite family of bases containing a basis simultaneously adapted to every such `q`-tuple.

### Proof

There are finitely many filtration types, so there are finitely many ordered `q`-tuples of types. For each realizable tuple, choose one simultaneous adapted basis. `square`

For the Ru--Vojta use, `q=2` is supplied by the Levin--Autissier simultaneous-adaptation lemma.

---

## 5. Why preselecting a basis does not lose the beta weight

Let `F` be a decreasing filtration with jumps

```text
alpha_1 > ... > alpha_m
```

and steps `F(alpha_j)`. If `s_1,...,s_l` is any adapted basis and

```text
mu_F(s_k) = sup { x : s_k in F(x) },                    (5.1)
```

then the multiset of weights `mu_F(s_k)` depends only on the dimensions of the filtration steps, not on the chosen adapted basis.

Indeed,

```text
# { k : mu_F(s_k) >= x } = dim F(x).                    (5.2)
```

Consequently

```text
sum_k mu_F(s_k)
```

is the same for every adapted basis. This is the finite-level weight sum used in the beta calculation.

The divisor or ideal membership estimate also uses only the statement

```text
s_k in F(mu_F(s_k)),
```

which holds for the preselected adapted basis. Therefore replacing an arbitrary adapted basis by the chosen representative for its filtration type preserves the required lower bound and preserves the total filtration weight.

This is the crucial point: the individual sections may differ, but the Ru--Vojta filtration estimate does not require optimizing over all bases once one adapted representative has been fixed for each type.

---

## 6. Consequence for the local maximum

At the fixed level, let

```text
B_1,...,B_M
```

be the preselected finite family from Section 4. For each place `v` and point `P`, every basis that the Autissier argument needs may be replaced by one of these `B_j`. Hence the local estimate can be written using

```text
max_{1 <= j <= M} sum_{s in B_j} lambda_{s,v}(P).        (6.1)
```

The family of sections occurring in (6.1) is fixed and finite, independent of

- the finite place set `S`;
- the numerical local weight vector `t`;
- the point `P`;
- the conductor support.

The real weights and the maximizing index may vary continuously or pointwise. The underlying linear forms do not.

Thus the finite-form-family assertion in the earlier proof is valid after replacing the unsupported phrase "there are finitely many adapted bases" by the proved statement:

> There are finitely many fixed-level filtration types, and one may preselect one adapted basis for each realizable type (or tuple of types).

---

## 7. What this repairs, and what remains

This proves the missing finite combinatorial reduction. In particular it repairs the first of the three suspect assertions identified in `claimed_complete_proof_salvage_audit.md`.

It does **not** by itself prove either of the following:

1. that the pointwise maximum in (6.1) is represented by one global twisted-height system independent of `P`;
2. that the absolute parametric Subspace Theorem supplies a height threshold and a complete exceptional set uniform in `S` for the resulting finite maximum.

The next step should exploit the now-fixed finite family directly. Since `M` is independent of `S`, one can partition the points by their maximizing basis pattern or encode the maximum by a finite product construction. Any remaining loss must now come from placewise patterns across `S`, not from an infinite family of sections.

---

## 8. Formalization target

The abstract finite-type theorem needs only finite sets, linear orders, and finite sums of subspaces. A Lean module can be organized around:

```text
finite filtration ground set B
score t b
upperSet t x = {b | x <= score t b}
step t x = span/sup of V_b over upperSet t x
```

and prove that `step t x` belongs to the finite image of the powerset of `B`. The stronger hyperplane-cell statement is unnecessary for the basic finiteness result: every step is already the sum of a subset of `B`, so at most `2^|B|` distinct subspaces occur, and every filtration chain is a chain in this finite set.

This gives an even shorter purely combinatorial formalization.