# Real-order audit for the polynomial division model

This is a bounded check of whether the pointwise order forced by real
pair differences improves the incidence bound in
`polynomial_divisor_incidence_bound.md`.

## What the order gives

For every pair, the divisor identity implies

```text
(F_i-F_j) | F_i^2+1.
```

Thus `F_i-F_j` has no real zero. Its degree is the common even degree
`D`, so its sign is constant on the real line. Consequently the rows
can be totally ordered pointwise (up to reversing the order):
`F_1(x)<...<F_m(x)` for every real `x`.

This order does not produce a new half-plane root count. In fact, for
every real polynomial `F` of even degree `D`, `F+i` has exactly `D/2`
roots in each open half-plane, with multiplicity. Homotope `F` through
real degree-`D` polynomials to its leading term. A root cannot cross
the real axis during this homotopy, since `F(x)+i` is never zero for
real `x`; the leading monomial has `D/2` roots on each side.

There is also a direct ordered-pair version. If `F_i<F_j` on `R`, then

```text
R_ij(z)=(F_i(z)+i)/(F_j(z)+i)
```

has positive imaginary part on the real axis, because
`Im R_ij(x)=(F_j(x)-F_i(x))/(F_j(x)^2+1)>0`. Its values at both ends of
the real axis approach the same real ratio of leading coefficients
(the degree is even). The argument principle therefore gives equal
numbers of upper-half-plane zeros and poles. This merely recovers the
automatic `D/2` count for each numerator and denominator; it does not
couple the rows strongly enough to improve the existing degree-dependent
bound.

## Strict positivity is false

The quotient suggested as a possible positivity obstruction can vanish
on the real line. Set

```text
F_-(t)=-(t^2+1)+t,       F_+(t)=t^2+1+t.
```

Then `F_-<F_+` everywhere,

```text
F_--F_+=-2(t^2+1),
(F_-F_++1)/(F_--F_+)=t^2/2.
```

Hence the divisor condition holds exactly, while the quotient has a
real zero at `t=0`. For this two-row example the lcm condition also
holds: both `F_-+i` and `F_++i` have degree two, so their lcm degree is
at most four (in fact it is three), which is at most `2D`.

The quotient can therefore be nonnegative in special cases, but no
strict positivity or universal no-real-zero assertion is available.

## Consequence for the bound

The four-row quartic construction in `four_row_quartic_eight_block_family.md`
already has constant-sign pair differences, hence admits the same
pointwise total order after sorting its rows, while satisfying
`deg lcm(F_i+i)=8=2D`. It has `m=4`; therefore a degree-independent
bound such as `m<=3` is false even with the real order. The order and
the half-plane count provide no additional height or divisibility
factor. Any stronger bound would need a genuinely new relation beyond
pointwise ordering and the existing lcm/incidence hypotheses.

This conclusion concerns the faithful polynomial model only and makes
no claim about arbitrary integer arc configurations.
