# Independent phase-coordinate and multiplicity audit

The formulas in [the Newton cut-order note](gale_torus_newton_cut_orders.md)
are consistent in the raw phase coordinates `z_i`. In particular they
must not be recomputed by treating `z_i` as a cotangent without also
performing the corresponding fractional linear substitution and keeping
its factors.

For formal rational phases the squared chord is

```text
delta_ij^2=-(z_i-z_j)^2/(z_i z_j).
```

Thus `z_i=epsilon x_i` on a cut, with generic unit parameters on both
sides, has negative chord valuation one across the cut and zero within
each side. It is exactly the clean conductor normalization used by the
phase formula. The real expression is recovered when `|z_i|=1`;
finite-place orders use its rational-function expression.

The raw coordinates are
`R_I=(product_(i in I)z_i^2) product_(j<k in I^c)(z_k-z_j)`.
They have degree `16` and individual degrees at most `4`. Inversion
introduces ten minus signs, hence no sign, and gives
`R_I(1/z)=R_I(z)/product_i z_i^4`. Homogeneous combinations of degree
`d` inherit the exact coefficient symmetry `alpha -> 4d*1-alpha`,
including after cancellations.

Direct support enumeration gives the minima for cut sizes one through
seven as

```text
0, 1, 3, 5, 7, 9, 12.
```

The larger baselines obey `m_s=m_(8-s)+4s-16`, as required for equal
normalized orders on complementary cuts. On the chosen 127 cuts the
baseline is `373d`, and each coordinate gives width sum `320` and order
sum `53`. Generic minimum support really is the generic cut order:
distinct surviving exponent vectors remain distinct monomials in the
independent cut-unit parameters. Specializations can only raise that
order, subject to retaining arithmetic coefficient denominators.

## A bounded multiplicity result in degree one

There is no nonzero degree-one pullback of diagonal multiplicity greater
than ten. This follows from exact ranks, not a sampled numerical order.

Set `z_i=1+u_i`. The degree-ten leading part of `R_I` is the five-node
Vandermonde `V_(I^c)(u)`. Consider the two linear coefficient maps from
the 56 coordinate coefficients to respectively `R_I(z)` and
`V_(I^c)(u)`.

Both have a common 21-dimensional space of exact relations. To see the
relations, give the complementary Plucker coordinate its usual sign
`p_I=(-1)^(sum I)R_I`, where indices are zero-based. The moment matrix
contains a row of ones, so its kernel lies in the sum-zero hyperplane.
Contracting the Plucker trivector against that row gives, for each pair
`i<j`,

```text
sum_(k different from i,j) p_(k,i,j)=0,
```

where antisymmetry supplies the ordered-coordinate signs. The same
relations hold for the leading Vandermondes, either by taking leading
terms or by the polynomial moment matrix with rows of degrees zero
through four.

The checker verifies all these relations as exact polynomial identities
and verifies their rank is at least 21 using a prime-field minor. It also
finds rank at least 35 for each of the raw and leading coefficient maps
over that same prime. These lower bounds, together with the exact
relations and `56=35+21`, prove that both maps have rational rank exactly
35 and have the identical 21-dimensional kernel. Therefore a nonzero
linear combination of raw coordinates has a nonzero degree-ten leading
part and multiplicity exactly ten. The argument also applies over the
complex numbers.

This rank argument only rules out excess diagonal multiplicity in degree
one. The separate, subsequently proved
[support-multiplicity theorem](polarized_coefficient_support_multiplicity.md)
now establishes `mult_1(f)<=sum_S h_S/32` in every degree, including
all row weights together. Its independent proof audit is recorded in
[the support-multiplicity audit](polarized_support_multiplicity_independent_audit.md).

[check_gale_torus_newton_cut_orders_independent.py](check_gale_torus_newton_cut_orders_independent.py)
performs the exact support, relation, and rank checks. Prime-field ranks
are used only as characteristic-zero lower bounds, alongside the exact
upper bounds from the verified polynomial relations.
