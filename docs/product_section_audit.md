# Product-section audit of the max-over-adapted-bases step

## Main simplification

At the minimal finite level `(N,d)=(1,6)`, let `V` be the fixed 27-dimensional section space.  The Ru--Vojta local term is obtained by maximizing, over a fixed finite family of adapted bases `B`, a sum of logarithmic section evaluations.

For one basis `B`, write

```text
P_B = product_{s in B} s.
```

Then pointwise

```text
max_B sum_{s in B} log |s(P)|_v
  = log max_B |P_B(P)|_v.
```

Thus the max-over-bases expression can be replaced exactly by the ordinary maximum norm of the fixed finite coordinate family `{P_B}`.  The point chooses the maximizing coordinate; it does not choose a new form family or a new weight system.

This suggests using the projective map

```text
Phi : X --> P^{M-1},
Phi(P) = [P_B(P)]_B,
```

where `M` is the number of adapted bases.  All coordinates have the same degree `27 * 6` on `P^2` before removing the common exceptional multiplicity, so this is a fixed higher-degree projective system.

## Why this may remove the earlier obstruction

The problematic quantifier was

```text
for every point P there exists an adapted basis B(P).
```

After multiplying the sections in each basis, this becomes simply

```text
max_B |P_B(P)|_v,
```

which is the standard local sup norm of one fixed coordinate vector.  No pair-dependent choice remains in the data of the Diophantine theorem.

The Evertse--Ferretti projective-variety theorem is natively stated for fixed homogeneous forms and gives proper subvarieties with uniform number and degree bounds once the form family, variety and approximation margin are fixed.  This is therefore a better interface than linearizing every section separately and tracking point-dependent adapted flags.

## Conjugate-place pairing

The earlier calculation showed that for a split Gaussian prime and its conjugate, the total filtration weight of an adapted sextic basis is independent of the valuation ordering.  If that identity applies basis-by-basis, then all product coordinates `P_B` receive the same conductor exponent after pairing `pi` and `bar pi`.

In that case the paired local contribution has the form

```text
log max_B |P_B(P)|_{pi,bar pi} - C_pi log Q_G,
```

with `C_pi` depending only on the conductor exponent, not on the point and not on `B`.  The `Q_G` factor can then be pulled outside the maximum.

This would eliminate the need for a fixed-weight domination estimate and its `1/18` loss: the conversion is exact.

## Items that must be checked before claiming success

1. The local Ru--Vojta expression is exactly a maximum of sums over complete adapted bases, with no basis-dependent additive normalization omitted.
2. Every product `P_B` has the same line-bundle degree and the same vanishing multiplicity at the blown-up point.
3. The conjugate-place total weight is constant for every adapted basis, not merely after summing the multiset of filtration weights abstractly.
4. The product-coordinate system has no base locus on the divisor-pair points under consideration, or the base locus is itself a fixed proper exceptional set.
5. The projective higher-degree theorem accepts the resulting maximum/sup-norm inequality directly; otherwise apply it to the image variety `Phi(X)` with coordinate linear forms.
6. The height comparison between `h(Phi(P))` and the common conductor height has the correct direction and does not consume the strict `1/18` margin.

## Expected shortest route if the checks pass

```text
explicit sextic adapted bases
  -> product sections P_B
  -> one fixed projective map Phi
  -> exact paired-place conductor inequality
  -> quantitative Subspace Theorem on Phi(X)
  -> uniformly bounded exceptional subvarieties
  -> grid/fiber bound.
```

This route would make the Harder--Narasimhan, finite-flat and varying-weight modules non-load-bearing for the final proof, while leaving them as valid side results.
