# Polynomial coefficient loss in the moving six-point radius

This note makes the fixed-frame comparisons for both radius and angular
diameter polynomially uniform in the central weights after choosing the
controlled hyperbolic frame.  These comparisons do not by themselves
exclude endpoint-scale configurations.

Let

```text
u=(u_0,...,u_5),       U=max_i |u_i|,
```

be a primitive nonzero integer vector with total sum zero and no vanishing
nonempty proper subsum.  Assume `diag(u)` is hyperbolic over `Q`, as it is
when an actual six-node central circle configuration supplies a maximal
isotropic rowspace.  By
[`rational_hyperbolic_frame_polynomial_height.md`](rational_hyperbolic_frame_polynomial_height.md),
one can choose and clear a rational hyperbolic frame containing the all-one
vector so that the resulting integer pencil has coefficient height

```text
2 <= T <= C_0 U^C.                                  (1)
```

Here and below the exponents and constants are absolute because the
dimension and all polynomial degrees are fixed.

## 1. The fixed-frame constants cost only a power of `T`

The raw cofactor forms and the associated binary forms

```text
Delta=4AC-B^2,             deg Delta=12,
K=C c1^2-B c1 c2+A c2^2-c0 Delta,    deg K=20
```

have integer coefficients bounded by `T^O(1)`.  This follows directly from
their fixed-size determinant formulas.

The twenty minor quadratics have a rank-three coefficient matrix.  Choose
a nonzero three-by-three minor of that matrix.  Its adjugate expresses its
nonzero integer determinant times each of `s^2,st,t^2` as an integer linear
combination of three of the quadratics.  Thus the constant `E_V` in the
minor-content bound may be chosen with

```text
1 <= E_V <= T^O(1).                                  (2)
```

The forms `Delta` and `K` have no common projective root.  The homogeneous
Sylvester resultant and its Bezout cofactors give integers `E_R!=0` and a
fixed exponent `r` for which

```text
E_R s^r, E_R t^r belong to the integer ideal (Delta,K).
```

Every entry of the relevant Sylvester matrices is a coefficient of
`Delta` or `K`.  Their sizes are fixed, so their determinants and adjugate
minors show that `E_R` can be chosen with

```text
1 <= |E_R| <= T^O(1).                                (3)
```

Finally, the maximum of `|K|` on `max(|s|,|t|)=1` is at most `T^O(1)`.
The corresponding positive minimum on the real locus `Delta>=0` is also
at least `T^(-O(1))`.  Here is an elementary proof of the latter statement.
Under the nonresonance hypothesis, `Delta` and `K` are squarefree and
coprime, so the degree-32 binary form

```text
f=Delta K
```

is squarefree.  In each of the two standard bounded charts, dehomogenize
`f`.  The resulting squarefree integer polynomial has degree at most 32 and
coefficient height `T^O(1)`.  Its nonzero integer discriminant has absolute
value at least one.  Cauchy's root bound, followed by the discriminant
product formula, therefore gives pairwise separation

```text
|alpha-beta| >= T^(-O(1))                            (4)
```

for its distinct finite roots.  This argument remains valid if the degree
drops in one chart; the projective root at that chart's infinity occurs in
the other bounded chart.

A nonreal root of `K` stays at least half the separation (4) from the real
axis, because its conjugate is another root.  At a real root of `K`, one has
`Delta<0`: a real triangle factor of `K` vanishes there, so that triple and
its complementary triple lie on two real component lines.  The lines are
nonparallel, whereas a parallel pair is forbidden by nonresonance.  Between
this root and any real parameter with `Delta>=0` lies a real root of `Delta`.
Consequently every allowed real parameter is `T^(-O(1))` away from every
root of `K`.  Factoring the dehomogenized `K` in the bounded chart, whose
nonzero leading coefficient is an integer, gives

```text
T^(-O(1)) <= |K(s,t)| <= T^O(1)
when max(|s|,|t|)=1 and Delta(s,t)>=0.                (5)
```

The argument is homogeneous and includes either parameter infinity through
the other chart.

## 2. Moving-frame radius comparison

Take a primitive integer parameter `(s,t)`, put
`H=max(|s|,|t|)`, and suppose `d^2=Delta(s,t)>0` gives a nonsingular circle
with six distinct nodes.  The exact radius and content formulas in
[`six_point_circle_cover_radius_content.md`](six_point_circle_cover_radius_content.md)
give

```text
|K(s,t)|/(E_V^2 |E_R|^(3/2)) <= N <= |K(s,t)|,       (6)
```

where `N` is the least primitive squared radius.  Combining (2)--(6) and
homogeneity of `K` yields

```text
T^(-C_1) H^20 <= N <= T^C_1 H^20                    (7)
```

after enlarging one absolute exponent `C_1`.  Substitution of (1) gives the
moving-weight consequence

```text
U^(-C_2) H^20 <= N <= U^C_2 H^20.                   (8)
```

The constants refer to the actual controlled frame chosen above.  Equation
(8) is only a polynomial coefficient-height loss and supplies no uniform
rational-point bound on the genus-five cover.  Section 3 controls the
angular comparison separately.

For orientation, the fair six-point scales have
`U=exp(7w+o(w))` and least endpoint cotangent `X=A/L=exp(8w+o(w))`.
A separate theorem, restricted to configurations with
`theta N^(1/4)<=C` for fixed `C`, of the form

```text
X <= U^(c+o(1))              for some c<8/7          (MISSING)
```

would contradict those scales.  Nothing in (8) proves this missing height
comparison.

## 3. Exact chords and polynomial angular constants

The angular comparison can also be made polynomially uniform in `U`.  For
two labels put

```text
dx_ij=x_i-x_j,       dy_ij=y_i-y_j,
P_ij=A dx_ij^2+B dx_ij dy_ij+C dy_ij^2.              (9)
```

Thus `P_ij` is an integer binary form of degree eight and coefficient
height `T^O(1)`.  From the exact lift

```text
Z_i=d(2A x_i+B y_i+c1)+i(Delta y_i+delta)
```

one obtains

```text
|Z_i-Z_j|^2
 =d^2[(2A dx_ij+B dy_ij)^2+d^2 dy_ij^2]
 =4A Delta P_ij.                                    (10)
```

All `Z_i` have norm `4AK`.  Therefore the normalized unit-circle nodes
`q_i=Z_i/Z_0` satisfy the exact identity

```text
|q_i-q_j|^2=Delta P_ij/K.                            (11)
```

On the real locus `Delta>0`, the signs in (11) agree and every `P_ij` is
nonzero because the quadratic form in (9) is definite and the nodes are
distinct.  If

```text
Lambda=h Norm_G(e) Norm_G(f),       N=|K|/Lambda,
```

and `z_i` are the primitive Gaussian integer circle points of squared
radius `N`, (11) equivalently gives the exact integer edge formula

```text
|z_i-z_j|^2=d^2 |P_ij|/Lambda.                       (12)
```

It remains to bound the largest `P_ij` at `Delta=0`, where the quadratic
form has rank one.  Set

```text
S=sum_(i<j) P_ij^2.                                  (13)
```

There is no real projective zero of `S` with `Delta>=0`.  For `Delta>0`
this follows from definiteness and distinctness.  At `Delta=0`, make a
fixed affine coordinate change so that `A!=0`.  Then

```text
4A P_ij=(2A dx_ij+B dy_ij)^2.
```

If every `P_ij` vanished, the affine linear form `2Ax+By` would take the
same value on all six columns, placing all columns on one line.  This
contradicts the rank-three hyperbolic frame.

This nonvanishing has a polynomial quantitative form.  In either bounded
parameter chart, take the primitive squarefree part of `Delta S`.  Its
degree is fixed and its integer coefficient height is `T^O(1)`; this follows
from the fixed-size subresultant matrices for `gcd(Delta S,(Delta S)')`.
Cauchy's bound and the discriminant product formula give
`T^(-O(1))` separation between its distinct roots.  A nonreal root of `S`
is separated from the real axis by its conjugate.  At a real root of `S`,
one has `Delta<0`, so an intervening real root of `Delta` separates it from
every real point with `Delta>=0`.  Factoring `S` and using its nonzero
integer leading coefficient now gives

```text
T^(-O(1)) <= S(s,t) <= T^O(1)
when max(|s|,|t|)=1 and Delta(s,t)>=0.                (14)
```

Together with (5), equations (11), (13), and (14) imply, for
`H=max(|s|,|t|)`,

```text
T^(-C_3) |d|/H^6 <= D <= T^C_3 |d|/H^6,             (15)
D=max_(i<j)|q_i-q_j|.
```

After using the controlled frame (1), this becomes

```text
U^(-C_4) |d|/H^6 <= D <= U^C_4 |d|/H^6.             (16)
```

Let `theta` be the length of the shortest circle arc containing the six
nodes.  Whenever `theta<=pi`, its two endpoint nodes realize the maximum
chord, so

```text
D=2 sin(theta/2),       D<=theta<=(pi/2)D.           (17)
```

In particular, (8), (16), and (17) give the polynomially uniform endpoint
comparison

```text
U^(-C_5) |d|/H <= theta N^(1/4)
                       <= U^C_5 |d|/H                (18)
```

throughout the small-arc regime.  A bounded endpoint-scale sequence with
`N->infinity` is eventually in this regime automatically.  Hence, for an
absolute exponent `C_5` and fixed endpoint bound `C_end`, every such sequence
must give primitive integer solutions of the purely integral conditions

```text
d^2=Delta_u(s,t),       H=max(|s|,|t|),
gcd(s,t)=1,             0<|d| <= C_end U^C_5 H.      (19)
```

Conversely, restrict to the same `Delta>0` distinct-circle locus.  After
enlarging the absolute exponent, the stronger inequality
`|d|<=U^(-C_5)H` makes (16) give `D<=1`.  Relative to any one node, all
other nodes then have principal angular distance at most `pi/3`, so they
lie in an arc of length at most `2pi/3`; the convention in (17) applies.
The upper bound in (18) is consequently absolute, proving a bounded
normalized endpoint scale.  Ruling out the necessary square-value window
(19), or narrowing the polynomial gap between necessity and sufficiency,
remains a uniform arithmetic problem; no frame-dependent compactness
constant is hidden in this formulation.
