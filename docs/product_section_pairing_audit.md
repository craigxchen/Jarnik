# Product-section and conjugate-place audit

## 1. Exact algebraic reduction

Let `B` be a complete basis of the fixed sextic section space `V`, and define

```text
P_B = product_{s in B} s.
```

For every place `v` and point `x` outside the zero set of the basis sections,

```text
max_B sum_{s in B} log ||s(x)||_v
  = log max_B ||P_B(x)||_v.
```

Thus a local maximum over complete bases can be replaced by a local maximum over
one coordinate for each product section.

This identity is purely pointwise and allows different maximizing bases at
different places.

## 2. Conjugate places do not need the same maximizing basis

At a split rational prime `p = pi * conjugate(pi)`, the two places may choose
bases `B` and `C` independently. The paired local contribution is

```text
log max_B ||P_B(x)||_pi + log max_C ||P_C(x)||_conj(pi)
 = log max_{(B,C)}
     (||P_B(x)||_pi * ||P_C(x)||_conj(pi)).
```

So the correct finite local index set is the Cartesian square of the adapted
basis family. There is no need to force one basis to maximize simultaneously at
both conjugate places.

This is compatible with the higher-degree Subspace Theorem, which permits a
fixed finite list of forms at each place. What must remain uniform is the finite
set of distinct global product sections, not the maximizing index.

## 3. Basis weight sums

For a fixed filtration

```text
V = F^0 superset F^1 superset ...,
```

any basis adapted to the filtration has the same multiset of filtration weights.
Indeed, the number of basis elements of weight at least `m` equals `dim F^m`.
Consequently the total basis weight is intrinsic:

```text
sum_{s in B} wt(s) = sum_{m >= 1} dim F^m.
```

Therefore there is no basis-dependent conductor normalization within one local
filtration. If the conjugate filtration is the reflected filtration, its total
weight is likewise intrinsic. The earlier number `156 rho` should therefore be
understood as a filtration invariant, not as a property of a specially chosen
basis.

The remaining calculation is to verify the exact reflected-filtration formula
for the sextic divisor-pair filtration and to identify the common exponent in
the normalization used by the Diophantine theorem.

## 4. Important limitation: finiteness of adapted bases

The set of every basis adapted to a filtration is infinite. Product sections
`P_B` are therefore useful only if the Ru--Vojta maximum may be restricted to a
fixed finite family of adapted bases or replaced by a finite family up to
uniformly bounded local constants.

A finite flag family alone is not enough: a single flag admits infinitely many
adapted bases, and products of basis vectors vary nontrivially under triangular
basis changes.

Thus the next required lemma is one of the following.

### Finite basis lemma

There is a fixed finite collection `mathcal B` of sextic bases such that every
local filtration occurring in the divisor-pair problem has a basis in
`mathcal B` achieving the finite beta lower bound.

### Finite product domination lemma

There is a fixed finite collection of product sections `P_1,...,P_M` and a
constant `C`, independent of the place, conductor and point, such that

```text
max_{all adapted B} ||P_B(x)||_v
  <= C_v max_{1 <= j <= M} ||P_j(x)||_v,
```

with the global product of the constants contributing only `O(1)` or a loss
below the available endpoint margin.

Without one of these statements, the product-section reduction does not yet
produce a fixed finite projective map.

## 5. Why determinant invariance does not directly solve this

The exterior product of a complete basis is invariant up to determinant, but
evaluating all sections at one point gives a one-dimensional target. The wedge
of their evaluations therefore vanishes as soon as the basis has more than one
member. The useful local quantity is the product of scalar evaluations, not the
exterior determinant.

Hence the infinite triangular-change freedom cannot be removed simply by
passing to the determinant line.

## 6. Revised audit status

The product-section idea successfully removes the requirement that conjugate
places use the same maximizing basis and confirms that total filtration weights
are basis-independent.

It does not yet remove the main finiteness issue. The decisive question is now:

```text
Can the finite beta inequality be witnessed by a fixed finite collection of
explicit sextic bases?
```

If yes, the proof can use a single fixed higher-degree form family and bypass
the varying-weight and canonical-subspace machinery. If no, a finite-product
domination estimate with controlled loss is required.
