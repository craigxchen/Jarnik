# The first symmetric sum on a full three-block cube

This note records the one positive-sum certificate that survives exact
algebra.  It applies to a fully factorized Boolean cube of eight points;
it does not apply to an arbitrary full-cut eight-point profile, whose row phases
can move after omitted blocks are absorbed.

Let

```text
z_e = d product_{j=1}^3 kappa_j^{(1+e_j)/2}
                     bar(kappa_j)^{(1-e_j)/2},   e in {+1,-1}^3,
rho = |d| product_j |kappa_j|,
```

where `kappa_j=a_j+i b_j` are nonunit, conjugate-primitive Gaussian
integers and `d` is nonzero.  Thus every `z_e` has modulus `rho`.
The first elementary symmetric sum factors exactly:

```text
S_1 = sum_e z_e = d product_j (kappa_j+bar(kappa_j))
    = 8 d product_j a_j.                                  (1)
```

Writing `u_j=b_j^2/(a_j^2+b_j^2)`, the exact deficit is therefore

```text
64 rho^2-|S_1|^2
 =64 rho^2 (1-product_j(1-u_j))
 >=64 rho^2 u_j.                                          (2)
```

The left side is also the equal-modulus identity

```text
64 rho^2-|S_1|^2 = sum_{e<f}|z_e-z_f|^2.                 (3)
```

If these eight points lie in an arc of length `L <= C sqrt(rho)`, (3)
is at most `28 C^2 rho`.  Since every nonunit conjugate-primitive block
has `|b_j|>=1`, (2) gives, for each `j`,

```text
|kappa_j|^2 >= 16 rho/(7 C^2).                            (4)
```

Multiplying (4) for the three blocks and using
`rho=|d| product_j|kappa_j| >= product_j|kappa_j|` gives the exact
radius bound

```text
rho <= (7 C^2/16)^3.                                     (5)
```

The displayed rational comparison has exact denominators: before any
reduction, if `C=A/B` then its right side is
`343 A^6/(4096 B^6)`.  At the algebraic stage the only introduced
denominators are the `|kappa_j|^2` in `u_j`; multiplying through by
their positive integer product gives (2) as an integer inequality.

At `C=1/2` this says `rho <= (7/64)^3 < 1`, impossible for a nonzero
Gaussian integer point.  The constants `28`, `64`, and `16/7` are exact;
there is no asymptotic or numerical step in this certificate.

## Why this does not settle a fair profile

Equation (1) is a product identity for the complete cube.  After a
general full-cut eight-point profile is restricted or regrouped, its retained
row phases have independent row twists, so the first sum is a sum of
monomials with no product factorization.  Its Archimedean consequence

```text
|sum_i z_i|^2 >= m^2 rho^2 - sum_{i<j}|z_i-z_j|^2
```

is only the usual positive first-moment bound.  It gives no primewise
valuation of the sum: Gaussian summands can cancel modulo a selected
prime even when every summand is primitive.  The exact denominator audit
in [sum_kernel_denominator_budget.md](sum_kernel_denominator_budget.md)
shows that natural positive Taylor remainders retain a large real-part
LCD, while [reduced_sum_kernel_cancellation.md](reduced_sum_kernel_cancellation.md)
gives actual circle configurations where the reduced permanent loses a
whole local pole.  Thus positivity cannot be promoted to a uniform
primewise radius bound by itself.

The full-cube certificate is consequently a rederivation of the fixed
profile obstruction already proved in
[semibent_eight_point_phase_obstruction.md](semibent_eight_point_phase_obstruction.md),
where balanced cuts provide an even stronger fourth-power certificate.
It is useful as an exact diagnostic for an actually factorized profile,
but it supplies no extraction theorem for arbitrary fair profiles and
makes no endpoint-family claim.

The same independence gives power-sum factorizations
`sum_e z_e^ell = d^ell product_j(kappa_j^ell+bar(kappa_j)^ell)`.
Newton identities turn these into elementary symmetric sums, but the
factors can vanish or acquire unrelated prime divisors.  No additional
uniform conclusion follows without a new structural relation among the
row phases.
