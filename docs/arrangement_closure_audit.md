# Audit of the proposed arrangement-closure identity

## Proposed identity

For the canonical Evertse--Ferretti exceptional subspace `T(c)`, we proposed

```text
T(c) = intersection { ker L : L in mathcal L and L|_{T(c)} = 0 }.
```

This identity is false for the original finite form set `mathcal L`.

## Counterexample already implicit in Evertse--Ferretti

In the special system

```text
mathcal L = {X_1, X_2, X_3, X_1 + X_2 + X_3},
```

Evertse--Ferretti, Lemma 15.3, allows canonical spaces of the form

```text
T = {x : sum_{j in I} x_j = 0}
```

for nonempty subsets `I`.  Take `I = {1,2}`:

```text
T = {x in Q^3 : x_1 + x_2 = 0}.
```

No member of the displayed `mathcal L` vanishes identically on `T`:

- `X_1` and `X_2` are nonzero on `(1,-1,0)`;
- `X_3` is nonzero on `(0,0,1)`;
- `X_1+X_2+X_3` is nonzero on `(0,0,1)`.

Hence the intersection over forms in `mathcal L` vanishing on `T` is an empty
intersection, namely all of `Q^3`, not `T`.

Therefore the proposed arrangement-flat identity cannot be used.

## What remains true

Evertse--Ferretti prove a different finiteness statement.  Their canonical
subspace `T(L,c)` satisfies a height bound depending only on the fixed finite
form family `L`:

```text
H_2(T) <= (max H_2(L_i))^(4n).
```

Since `T` is defined over the fixed number field, Northcott finiteness for the
Grassmannian implies that `T(L,c)` belongs to a finite collection depending
only on `L`, even as the local weights `c` vary.

This is the correct source of finite-family uniformity.  It is strictly broader
than the intersection lattice of the original forms.

## Formalization note

`GaussianChain/ArrangementClosure.lean` proves a valid abstract theorem:
if a score is nondecreasing under a chosen closure operator and the canonical
object is the greatest maximizer, then it is closed.  The audit above shows
that the required score-monotonicity hypothesis is not available for the naive
closure generated only by the original forms, so that module must not be cited
as establishing the Evertse--Ferretti identity.
