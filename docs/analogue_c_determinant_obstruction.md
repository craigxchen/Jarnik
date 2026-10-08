# Analogue C determinant obstruction

## Proposed step

The determinant route proposed taking evaluation functionals

```text
ell_1, ..., ell_n in Hom_R(V,R)
```

over a DVR `R`, then relating the valuation of their exterior product to the
sum of pointwise optimal adapted-basis weights.

That implication is false in general.

## Counterexample

Let `R` be a DVR with uniformizer `p`, let `V = R^2`, and define

```text
ell_1(x,y) = x + p^M y,
ell_2(x,y) = p^M x + y.
```

The evaluation matrix in the standard basis is

```text
[ 1    p^M ]
[ p^M  1   ].
```

Its determinant is

```text
1 - p^(2M),
```

which is a unit. Hence

```text
v(ell_1 wedge ell_2) = 0.
```

But each functional separately has an integral kernel vector:

```text
ell_1(-p^M,1) = 0,
ell_2(1,-p^M) = 0.
```

Thus a basis adapted separately to `ell_1` can contain `(-p^M,1)`, while a
basis adapted separately to `ell_2` can contain `(1,-p^M)`. Pointwise
adaptation therefore produces arbitrarily large, even infinite, individual
vanishing weights although the joint exterior product has valuation zero.

## Exact lesson

A determinant permits:

- a different common column basis at each prime;
- row operations that preserve or control determinant valuation.

It does **not** permit a different column basis for each row at the same prime.
The Ru--Vojta local maximum is pointwise, so its adapted basis may vary with the
point. Consequently the beta gain cannot be inserted into a multipoint
determinant without an additional common-basis or clustering theorem.

The safe DVR statement is only

```text
v(det A) >= min_sigma sum_i v(a_{i,sigma(i)}),
```

for one fixed matrix `A` in one fixed common basis. Equality can fail because of
cancellation, but the lower bound is valid. It gives no access to separately
optimized row bases.

## Corrected finite-dimensional target

A determinant proof now requires one of the following genuinely stronger
properties of the conductor evaluation functionals:

1. **Common adapted basis:** for every prime, all rows in a selected tuple admit
   one basis giving the needed total filtration weight.
2. **Primewise clustering:** partition the rows into uniformly few classes,
   each with a common adapted basis, and combine block determinants.
3. **Structured row differences:** after subtracting a root row, the resulting
   functionals lie deeply in a fixed lattice filtration determined only by the
   conductor.
4. **Direct interpolation determinant:** choose explicit sections whose
   determinant factors into pairwise arithmetic differences, as in the existing
   Ramana/Vandermonde determinant.

Without one of these, the naive exterior-product/assignment identity is false.

## Consequence for Analogue C

The abstract finite-dimensional model remains useful, but the right model is
not arbitrary pointwise filtrations. It must include a compatibility condition
forcing many evaluation rows to share one lattice flag. Finding such a
compatibility condition is now the load-bearing step.
