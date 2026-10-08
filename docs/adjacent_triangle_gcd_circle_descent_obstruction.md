# A common adjacent-area gcd does not give a circle descent

The five-point bounded-determinant theorem controls the radius when the
**absolute** adjacent triangle determinants are bounded.  Bounded ratios
between those determinants, even together with an unbounded common gcd, do
not give the same conclusion.  Saturating the difference lattice produces
integral points on an integral conic, but generally changes the Euclidean
circle into an ellipse.  If the saturated image is again an integral
Euclidean circle, the operation is only removal of a common Gaussian factor.

An existing primitive five-point endpoint family realizes the obstruction
exactly: its consecutive triangle determinants have ratios below `4/3` and
common gcd tending to infinity, while its radius tends to infinity.  Thus
this proposed normalization does not improve the growth rate.

## 1. Exact metric after saturating the difference lattice

Let `z_0,...,z_(m-1)` be distinct points of `Z^2` on the centered circle
`|z|=R`, with `m>=5`.  Let

```text
Lambda=sum_i Z(z_i-z_0),
g=[Z^2:Lambda]=gcd_(i<j) |det(z_i-z_0,z_j-z_0)|.
```

The last gcd is also the gcd of all triangle determinants.  Choose an
integer basis matrix `A` for `Lambda`, so `|det A|=g`, and write

```text
z_i=z_0+A u_i,       u_i in Z^2.
```

The saturated coordinates are integral, but their exact conic is

```text
u^T G u+2 ell^T u=0,
G=A^T A,       ell=A^T z_0,       det G=g^2.          (1)
```

Its center is `-A^(-1)z_0`.  Formula (1) records the metric factor that is
lost if one merely divides all triangle areas by `g`.  The primitive conic
through any five of the `u_i` is unique, but uniqueness does not make `G`
a scalar matrix in the standard Euclidean metric.  More explicitly, its
primitive projective coefficient vector is

```text
(G_11,2G_12,G_22,2ell_1,2ell_2)/c,
c=gcd(G_11,2G_12,G_22,2ell_1,2ell_2).               (2)
```

Adding further points only adds equations satisfied by this same vector;
the index identity `det G=g^2` gives no upper bound for its height.

There is a sharp criterion.  The `u_i` lie on a Euclidean circle with
integral center if and only if

```text
A^T A=lambda I,       A^(-1)z_0 in Z^2.              (3)
```

For necessity, the conic (1) and the putative Euclidean circle contain five
distinct points.  Bezout, or the elementary uniqueness of a conic through
these five circle points, makes their equations proportional; comparison of
quadratic and linear terms gives (3).  Sufficiency is immediate.

If (3) holds, the two integer columns of `A` are orthogonal and have equal
length.  Hence, up to reflection, `A` is multiplication by a Gaussian
integer `alpha`, with

```text
lambda=N(alpha),       g=N(alpha).
```

The integral-center condition says `z_0=A c` for an integer vector `c`.
Since every `z_i-z_0` also belongs to `A Z^2`, `alpha` divides every `z_i`
(with conjugation in the reflected case).  Therefore:

```text
For a Gaussian-primitive source tuple, a saturated integral Euclidean
circle is possible only when g=1.                              (4)
```

When `g>1`, at least one of the anisotropic metric or the nonintegral center
is unavoidable.  A further rational map that restores the round metric and
an integral center is a rational similarity; its denominator must divide
the common Gaussian content of the original tuple.  Thus it cannot reduce
the radius after primitive normalization.

This is the exact similarity cost behind the more general affine descent
obstruction.  It uses all five points; normalizing a single triangle alone
can always create integral affine coordinates but cannot preserve the
circle problem.

### Centered affine reembedding corollary

The same argument closes every affine reembedding between centered circles,
without first choosing a difference-lattice basis.  Let `m>=5` distinct
points `z_i` lie on `|z|=R`, and suppose a real affine map

```text
T(x)=M x+b
```

sends them to points on `|w|=R'`, at least five of which are distinct.
The matrix `M` must be invertible: otherwise all image points lie on a line,
which meets a circle
in at most two points.  The inverse image under `T` of the target circle is
therefore a nonsingular ellipse.  It and the source circle have five
distinct common points, so conic uniqueness gives

```text
(Mx+b)^T(Mx+b)-R'^2=lambda(x^T x-R^2).
```

Comparing quadratic and linear terms yields

```text
M^T M=lambda I,       M^T b=0,
```

and hence `b=0`.  Thus `T` is multiplication by a nonzero complex scalar
`c`, possibly followed by complex conjugation, and `R'=|c|R`.

If the `z_i` and their images are Gaussian integers and
`gcd_G(z_0,...,z_(m-1))=1`, then `c` is a Gaussian integer.  In the
orientation-preserving case, write a Gaussian Bezout identity

```text
sum_i a_i z_i=1,       a_i in Z[i].
```

Since every `c z_i` is Gaussian integral,
`c=sum_i a_i(cz_i)` belongs to `Z[i]`.  In the reflected case apply the
conjugate Bezout identity to `bar(z_i)`.  Consequently `|c|>=1`: no affine
reembedding of a primitive Gaussian circle tuple with at least five points
into another centered Gaussian circle can decrease the radius.

Five is essential to this conic argument.  Four points can lie on two
different conics and retain the known Pell-type flexibility.  The statement
does not cover projective, fractional-linear, or other nonlinear maps.

## 2. Bounded ratios do not give a bounded alphabet

Let the consecutive doubled areas be

```text
D_j=|det(z_j-z_(j-1),z_(j+1)-z_j)|,
g_adj=gcd_j D_j.
```

Always `g|g_adj`, but equality need not hold: the adjacent gcd alone need
not even be the index of a common difference lattice.  The example below
has equality, so it tests the strongest possible version of the proposed
normalization.

Even when `g_adj=g`, dividing by this common factor only gives primitive
integers `D_j/g`.  A ratio bound

```text
max_j D_j/min_j D_j <= C_0
```

places them in an interval `[H,C_0H]` with uncontrolled lower endpoint
`H`.  Such an interval contains order `H` possible integers.  Consequently
a growing number of windows does not force a bounded five-point determinant
alphabet unless one separately bounds `H`.  Applying a pigeonhole argument
to ratios or projective area classes does not remove this scale.

The five-point theorem in `five_point_bounded_triangle_determinants.md`
therefore cannot be applied after the formal division by `g_adj`: that
division is an affine determinant normalization, while the theorem's radius
is Euclidean and depends on the metric in (1).

## 3. A primitive endpoint family with bounded adjacent-area ratios

Use the five Gaussian points in
`five_point_affine_shape_cubic_family.md`, with `T=600k`.  They are
Gaussian primitive, lie on an arc shorter than `5sqrt(R)`, and have radius

```text
R(T)=sqrt(product_(a=1)^6(T^2+a^2))/720 asymptotic to T^6/720.
```

Their full triangle gcd is

```text
g=T/15.
```

For the three consecutive triples `012,123,234`, division by `g` gives

```text
A_T=3(T^2+4)/4       =270000k^2+3,
B_T=7T^2/12          =210000k^2,
C_T=25(T^2+36)/36    =250000k^2+25.                 (5)
```

These three integers have gcd one.  Indeed a common divisor divides
`9B_T-7A_T=-21`.  The last value is never divisible by 3, since
`C_T=k^2+1 mod 3`; the first is never divisible by 7, since
`A_T=3(k^2+1) mod 7` and `-1` is not a square modulo 7.  Thus the adjacent
gcd is exactly the full lattice index `g`.

For every `k>=1`, `B_T<C_T<A_T` and

```text
A_T/B_T<4/3.                                        (6)
```

Hence all three raw consecutive areas are uniformly comparable, while
their common gcd tends to infinity.  Their primitive magnitude still grows
as `T^2`, and the radius grows as its third power.

The affine conic in the normalization sending the first three points to
`(0,0),(1,0),(0,1)` makes the metric cost visible.  With `x=T^2`, its
quadratic coefficients are

```text
A=4(x+25)(x+36),
B=10(x^2+43x+360),
C=25(x+9)(x+16),
AC-B^2=(180T)^2,                                    (7)
```

and the physical Gram matrix is

```text
(x+4)/3600 * [[A,B],[B,C]].                         (8)
```

For `T=600k`, its exact projective coefficient content is

```text
gcd(A,2B,C)=3600,
H_conic=C/3600=25(x+9)(x+16)/3600 asymptotic to T^4/144.
```

Indeed, after dividing by `3600` and putting `n=100k^2`, the first two
coefficients are `14400n^2+244n+1` and
`72000n^2+860n+2`.  Their indicated linear combination is
`-3(120n+1)`; combining the first and third similarly gives
`-3(1200n+7)`.  The first coefficient is nonzero modulo three, and the two
parenthesized integers differ after multiplication by ten by three, so the
remaining gcd is one.

Thus the conic has square metric determinant, but its primitive coefficient
height is of order `T^4` and the physical scale in (8) is of order `T^2`.
Neither factor is removed by the common area gcd.  Restoring the original
round metric simply reverses the rational affine normalization and restores
its denominator cost.

## 4. Consequence

Large common adjacent-area content is useful as an index statement when it
equals the full triangle gcd, but it does not create a smaller integral
lattice circle.  Bounded area ratios do not control the remaining primitive
area scale, and the explicit primitive endpoint family (5)--(8) shows that
both phenomena can occur together.

A many-point improvement from the five-point theorem would therefore need
an independent mechanism producing a window with bounded **absolute
primitive** areas, or a conformal divisor that is proved to divide every
Gaussian point.  Adjacent-area gcd normalization alone supplies neither.

The displayed counterfamily has five points, so it does not rule out a
separate global theorem saying that a growing point set must contain a
window with bounded absolute primitive areas.  It rules out deriving that
conclusion from bounded local ratios and a large adjacent gcd alone.  The
similarity obstruction (4), by contrast, applies to every number of points
at least five.

The exact audit `check_adjacent_triangle_gcd_circle_descent.py` reconstructs
the family, checks (5)--(8), verifies the adjacent and full lattice indices,
and tests the saturation metric formula (1).
