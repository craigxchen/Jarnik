# Coefficient height forced by collapse or a small-radius critical image

The [fixed-degree rational-map theorem](fixed_degree_rational_map_inflation.md)
excludes negligible-height non-pure maps on sufficiently many centrally
extracted rows. Its row threshold grows with degree. The elementary
bounds below remain valid when the degree itself grows. They show that
exact interpolation collapse costs the product of the source row
heights, and that retaining even one large source half-angle row on a
target circle of parent-scale radius requires exponentially large map
coefficients once the degree exceeds two. They do not exclude maps of
that height or prove a uniform endpoint bound.

## Setup

Let `P,Q in Z[X,T]` be coprime homogeneous forms of degree `D`, with
`P(1,0)=a>=1` and `Q(1,0)=0`, after reanchoring the output at phase
one. Let `H>=1` bound the absolute value of every coefficient of both
forms, and set `J=log H`. For source pairs `(x_i,t_i)` with
`gcd_Z(x_i,t_i)=1` and
distinct slopes `t_i/x_i`, assume

```text
x_i>0,   t_i!=0,   |log|x_i+i t_i|-W/4|<=E,
|t_i|<=exp(E),     W/4>=2E+log 4.                  (1)
```

Then `x_i>=exp(W/4-E)/2` and `|t_i/x_i|<=1/2`. The output phase is
`(P+iQ)/(P-iQ)`. Divide its coordinates by their integer gcd and, if
both primitive coordinates are odd, divide the Gaussian numerator by
`1+i`, retaining the resulting unit in the phase. Its primitive
Gaussian numerator is `K_i`. Every actual circle carrying the image
relative to the anchor has physical radius at least `|K_i|`.

## Exact collapse requires interpolation height

Suppose `Q(x_i,t_i)=0` at `r` distinct nonanchor rows. If `Q` has a
zero of order `e>=1` at the anchor and is not the zero polynomial,
then Gauss's lemma gives the literal integer-form divisibility

```text
T^e product_(i=1..r)(t_i X-x_i T) | Q(X,T),
D>=e+r.                                             (2)
```

In fact the coefficient bound is exact and needs no Mahler estimate.
Dehomogenize at `T=1` and write
`Q(X,1)=product_i(t_i X-x_i) S(X)` with nonzero `S in Z[X]`.
Let `s_k X^k` be the lowest-degree nonzero term of `S`. The
coefficient of `X^k` in `Q(X,1)` is exactly
`(-1)^r s_k product_i x_i`. Since `|s_k|>=1`, this one coefficient
already proves

```text
H >= product_(i=1..r) x_i,
J >= r(W/4-E-log 2).                             (3)
```

This applies to an attempted map sending all selected rows to its
anchor. Such a map collapses distinct points and is not an endpoint
construction; (3) quantifies the coefficient cost of that exact
interpolation shortcut.

The same coefficient argument covers a map that **retains** selected
directions. Let `g=U/V` be a fixed primitive rational comparison map,
where `U,V` are integer homogeneous forms of degree `e` and coefficient
height at most `B`. Suppose `f=Q/P` is not identically `g` but agrees
with it at the `r` selected slopes (and the comparisons are defined
there). The nonzero homogeneous difference

```text
A=QV-PU
```

has degree `D+e`, vanishes at those slopes, and has coefficient height
at most `2(e+1)BH`. Applying the lowest-nonzero-coefficient argument
to `A` gives the exact finite bound

```text
2(e+1)BH >= product_(i=1..r) x_i,
J >= r(W/4-E-log 2)-log[2(e+1)B].                (4)
```

In particular, for the identity comparison `g=T/X`, one has
`A=XQ-TP` and the sharper coefficient bound `height(A)<=2H`:
every nonidentity rational map fixing these `r` distinct source
directions satisfies `2H>=product_i x_i`. If `A` has an additional
zero at the anchor, its degree also forces `r+1<=D+e`; without that
condition `r<=D+e`. These statements exclude the identical map,
for which `A=0` and there is no coefficient obstruction.

## A degree-growing bound from one retained row

Let `rho=Res(P,Q)!=0`. The coordinate gcd at a primitive pair divides
`|rho|`, even for growing `D`; the binary Sylvester determinant gives

```text
log|rho| <= 2D J+D log(D+1).                       (5)
```

Suppose an image circle has physical radius
`R_new<=exp(W/2+kappa)`, with `kappa>=0`. There are two cases.
If

```text
J <= W/4-2E-log 8,                                (6)
```

then `y=t_i/x_i` satisfies `H|y|<=1/4`. Writing
`P(1,y)=a+sum_(j=1..D)p_j y^j`, the geometric-series bound and
`|y|<=1/2` give `|P(1,y)|>=1/2`. Thus, for any one row,

```text
|K_i| >= x_i^D/(2 sqrt(2)|rho|),
J >= (D-2)W/(8D)
     - [E+log 2+log(D+1)]/2
     - [log(2 sqrt(2))+kappa]/(2D).               (7)
```

If (6) fails, the alternative lower bound is simply

```text
J > W/4-2E-log 8.                                 (8)
```

Equations (7)–(8) are a finite dichotomy with every degree cost
visible. For centrally extracted rows with `E=o(W)`, if
`log(D+1)=o(W)` and `kappa=o(DW)`, every degree `D>=3` therefore
requires

```text
J >= (D-2)W/(8D)-o(W).
```

In particular, if `D` tends to infinity subexponentially in `W`,
then `J>=W/8-o(W)`. An anchor zero of order at least three has
`D>=3`, so the conclusion applies to highly flat critical maps
whose output points are retained on a parent-scale circle. The
degree-two coefficient is zero, consistent with its raw norm
having exactly the parent-radius exponent for one source row.

This is a **one-row** radius test. It does not control the all-row
Gaussian lcm of newly created denominator primes or assert that a
large-coefficient interpolation map preserves distinctness and a
short target arc.

[The exact checker](check_growing_degree_half_angle_height_barrier.py)
constructs degree-growing coprime maps with prescribed collapse roots
and nonidentity maps agreeing with `T/X` at up to seven selected rows.
It verifies the coefficient witness, resultant content, primitive
output pairs, and least Gaussian-lcm radii. These finite fixtures
illustrate the formulas; the inequalities above are proved generally.
