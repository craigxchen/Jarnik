# A differential divisor bound for a real polynomial endpoint jet

This is a bound for a restricted function-field model, not for arbitrary
polynomial endpoint tuples or integer circle arcs. Suppose that distinct
rows `z_j=x_j+i y_j` in `Q(i)[T]` satisfy

```
x_j^2+y_j^2=P,   deg x_0=2d,   d>=1,
deg(x_j-x_0)<=d,   deg y_j<=d.
```

Then there are at most **twelve rows**, independently of `d`. No common-gcd
hypothesis is needed. The same proof works over any characteristic-zero
real field. In particular it applies if a constant rotation makes the
common coefficients above degree `d` real. A general common complex jet
need not admit that rotation.

## Reduction to constant radial shifts

Write `H=x_0`. Comparing degrees in

```
2H(x_j-H)+(x_j-H)^2=y_0^2-y_j^2
```

shows that `x_j-H=c_j` is constant: a nonconstant difference of degree
`r<=d` would give left degree `2d+r>2d`. Set `S=P-H^2`, of degree at most
`2d`. Then

```
y_j^2=S-2c_j H-c_j^2.                         (1)
```

Each constant `c` permits at most two rows, corresponding to the two
polynomial square roots. Choose one root `y_c` for each represented
constant.

## A common differential divisor

Define the actual polynomial

```
K=4(H')^2 S-4HH'S'-(S')^2
 =4(H')^2 P-(P')^2.                          (2)
```

Its degree is at most `6d-2`. Substituting (1) and its derivative gives

```
K=4 y_c^2(H')^2
  -8(H+c)H' y_c y_c'
  -4 y_c^2(y_c')^2.                          (3)
```

Assume first that `K` is nonzero. No `y_c` is identically zero, since that
would make `P=(H+c)^2` and hence `K=0`. We claim the full product divides:

```
product_(represented distinct c) y_c | K.     (4)
```

Work over an algebraic closure. At a root `alpha`, equation (1) is a
quadratic equation in `c`, so at most two represented constants can have
`y_c(alpha)=0`. If only one does, with multiplicity `m>=1`, (3) gives
`ord_alpha K>=2m-1>=m`.

If two distinct constants occur, write their root multiplicities as
`m<=n`. Subtracting their identities (1) and differentiating gives
`ord_alpha H'>=2m-1`. Apply (3) to the root with multiplicity `n`.
Its three terms have respective orders at least

```
2n+4m-2,   2n+2m-2,   4n-2.
```

All are at least `2m+2n-2>=m+n`. This proves (4), including all root
multiplicities.

All but at most one represented constant have `deg y_c=d`: cancellation
of the coefficient of `T^(2d)` in (1) specifies at most one value of `c`.
If `n` distinct constants occur, (4) therefore gives

```
(n-1)d <= sum_c deg y_c <= deg K <= 6d-2.
```

Thus `n<=6`, giving at most twelve rows.

## The identically zero case

If `K=0`, then `P'^2=4(H')^2 P`. Since `H'` is nonzero, `P` is a square
in the rational-function field. A polynomial that is a rational-function
square is a polynomial square, so write `P=Q^2`. Characteristic zero
then gives `Q'^2=H'^2`, hence `Q=H+a` or `Q=-H+a` for a constant `a`.
Changing the sign of `Q` lets us use `Q=H+a`.

Equation (1) becomes

```
y_c^2=(a-c)(2H+a+c).                          (5)
```

The value `c=a` gives just the zero-root row. For two other distinct
values `c,b`, divide their squares by `a-c` and `a-b` and take constant
square roots over the algebraic closure. The resulting polynomial
squares differ by the nonzero constant `c-b`. Factoring that difference
forces both polynomials to be constant, contradicting `deg H=2d`.
Thus there is at most one nontrivial value of `c`, and at most three rows.

## Scope and verification

The existing rational five-row family in
[five_point_affine_shape_cubic_family.md](five_point_affine_shape_cubic_family.md)
has common leading terms `-T^6-i T^5-45T^4` before its constant
normalization. Any constant rotation making its leading coefficient real
leaves a nonreal `T^5` coefficient, so the reduction to (1) is an
additional restriction. This argument supplies no reduction of arbitrary
polynomial rows to that case and no uniform integer-point bound.

The later [complex-jet audit](complex_jet_scalar_discriminant_obstruction.md)
shows that this five-row fixture admits no replacement scalar polynomial
of the natural derivative-square degree: the exact minimum degree is
twenty rather than ten. It also charges the full polynomial content of
a moving rational rotation, which restores the primitive tuple.

The [standard-library checker](check_polynomial_real_jet_differential_divisor_bound.py)
verifies 420 identities in (2)--(3), the root-order inequalities, fourteen
literal four-row polynomial endpoint tuples with their complete product
divisors, and fourteen three-row fixtures in the `K=0` case. Root and Sol
independently audited the degree bounds, full shared-root multiplicities,
and zero-polynomial alternative. The proof of the degree-independent bound
is the multiplicity argument above, not a finite experiment.
