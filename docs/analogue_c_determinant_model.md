# Analogue C: the finite-dimensional determinant model

## Goal

Replace the qualitative/quantitative Subspace Theorem machinery by one finite
linear-algebra statement.

The intended conclusion is:

> If a common-conductor cluster contains sufficiently many points, then the
> evaluation vectors of a fixed finite section space are linearly dependent in
> a way that produces one nonzero fixed-degree section vanishing on the entire
> cluster.

The only deep inequality should be a determinant squeeze.

---

## 1. Abstract setup

Let `K = Q(i)`, let `X` be the two-dimensional parameter variety for ordered
pairs `(u_A,u_B)`, and let

```text
V = H^0(X,L)
```

be a fixed `r`-dimensional section space.  The minimal finite Ru--Vojta level is

```text
L = 6H-E,    r = 27.
```

Fix an integral lattice `V_O` in `V` over `O_K`.

For a point `P in X(K)`, evaluation gives a linear functional

```text
ev_P : V -> K.
```

For an ordered `r`-tuple `P_1,...,P_r`, choose an `O_K`-basis
`s_1,...,s_r` of `V_O` and form

```text
Delta(P_1,...,P_r) = det(s_j(P_i))_{1 <= i,j <= r}.
```

After multiplying each row by the fixed conductor denominator dictated by
`L`, this becomes a Gaussian integer.  Changing the integral basis multiplies
`Delta` by a Gaussian unit, so its norm and every finite-place valuation are
basis-independent.

This basis invariance is exactly why different primes may use different adapted
bases without creating a global compatibility problem.

---

## 2. Determinant squeeze

For every `r`-tuple from one conductor `G`, seek bounds

```text
log |Norm Delta| <= A log |G| + O(1),                 (arch)
log Norm(Delta O_K) >= B log |G| - O(1),             (finite)
```

with fixed constants `B > A`.

If `Delta != 0`, the principal ideal identity gives

```text
Norm(Delta O_K) = |Norm_{K/Q}(Delta)|.
```

The two inequalities then contradict each other for sufficiently large `|G|`.
Thus every `r` by `r` evaluation determinant vanishes.

Consequently the evaluation map

```text
V -> K^C,    s |-> (s(P))_{P in C}
```

has rank at most `r-1` for the whole cluster `C`.  Hence its kernel contains a
nonzero section `s`.  The fixed divisor `s=0` contains every point of `C` and
has degree bounded solely in terms of `L`.

The existing grid/fiber lemma then gives the desired uniform cardinality bound.

This route has no exceptional finite remainder and no exceptional-set component
count.

---

## 3. Archimedean estimate

All cluster points satisfy

```text
|u_A-1|, |u_B-1| <= C |G|^(-1/2).
```

Choose local coordinates `x=u_A-1`, `y=u_B-1` at the blown-up point
`[1:1:1]`.  For the sextic space `H^0(6H-E)`, an adapted basis has total
vanishing weight

```text
sum_j ord_E(s_j) = 50.
```

Therefore multilinearity and the standard alternant estimate should give an
archimedean gain of

```text
|G|^(-50/2) = |G|^(-25)
```

relative to the naive sextic denominator size.  The exact normalization must be
written in homogeneous Gaussian-integer coordinates, because this is where the
coefficient

```text
6 - 3*(50/27) = 4/9
```

must reappear.

The crucial point is that this estimate is simultaneous for an entire
`27`-tuple: it is an ordinary determinant estimate, not a pointwise choice of
an exceptional subspace.

---

## 4. Finite-place estimate

At a finite place `v`, the determinant valuation is intrinsic.  We may choose
any `O_v`-basis of the local lattice `V_O tensor O_v`.

For a chosen basis `s_1,...,s_r`, if

```text
v(s_j(P_i)) >= a_{ij},
```

then the determinant expansion gives the elementary tropical bound

```text
v(Delta) >= min_{sigma in S_r} sum_i a_{i,sigma(i)}.   (4.1)
```

Because an `O_v`-unimodular basis change has unit determinant, we may optimize
(4.1) over all local integral bases independently at every place.

Thus the entire mixed-basis issue has been reduced to one finite-dimensional
local statement:

### Local determinant beta lemma

For every `r` conductor-box points `P_1,...,P_r` at one split prime pair
`pi,bar(pi)`, there is an integral basis of `V` for which the assignment minimum
in (4.1), after adding the conjugate-place contribution, is at least the sextic
beta contribution predicted by `50/27`.

Equivalently, the exterior product

```text
ev_{P_1} wedge ... wedge ev_{P_r} in det(V)^*
```

has conductor valuation at least the required amount.

This is now a finite linear-algebra theorem over a DVR.  It contains no global
heights, no Zariski exceptional sets, and no dependence on the number of
places.

---

## 5. Why the infinite adapted-basis family is harmless here

The product-section route required replacing all adapted bases by finitely many
products, which was not justified.

The determinant route does not require a finite basis family.  At each finite
place and for each selected tuple we may choose any integral adapted basis.
Only the determinant valuation matters, and it is unchanged by unimodular basis
change.

This is a genuine simplification, not a reformulation of the same quantifier
problem.

---

## 6. Exact finite theorem to prove next

Let `R` be a DVR with uniformizer `pi`, fraction field `F`, and let `Lambda` be
a free rank-`r` `R`-module.  Let `ell_1,...,ell_r in Lambda^* tensor F` be the
row evaluation functionals supplied by `r` conductor-box points.

Define

```text
w(ell_1,...,ell_r)
  = v(ell_1 wedge ... wedge ell_r).
```

The required theorem is a lower bound for `w` in terms of the three divisor
filtrations at the blown-up point.  At the sextic level the target total is the
integer `50`.

A useful purely combinatorial formulation is:

```text
w = sup over integral bases e_1,...,e_r
      min_{sigma in S_r} sum_i v(ell_i(e_{sigma(i)})).
```

The `<=` direction is immediate from determinant expansion; the reverse
direction is a Smith-normal-form / valuated-matroid statement.  Once this
identity is established, the local beta lower bound becomes a finite assignment
problem.

---

## 7. Proof architecture if the local lemma succeeds

1. Fix `V = H^0(6H-E)`, `dim V = 27`.
2. For every `27` points from one conductor, form the normalized evaluation
   determinant.
3. At each split prime pair, use the local determinant beta lemma.
4. Sum valuations over primes; no factor depending on `|S|` appears.
5. Use the archimedean endpoint estimate and the strict `4/9 < 1/2` margin.
6. Conclude every `27`-tuple determinant is zero.
7. Obtain one nonzero sextic section vanishing on the whole pair grid.
8. Apply the existing grid/fiber bound.

The sole load-bearing new result is now a finite DVR exterior-product inequality.
