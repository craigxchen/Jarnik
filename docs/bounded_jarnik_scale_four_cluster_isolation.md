# Uniform isolation of bounded Jarnik-scale four-point clusters

Four points in an arc of length `K R^(1/3)` are isolated from every additional
point on the same circle at a definite endpoint distance: if the fifth point
is at distance `d` from any of their anchors, then

```text
d+K R^(1/3) >= 4 K^(-3/2) sqrt(R).
```

Section 6 proves this explicit inequality using the five-point determinant,
whose exact pencil calculation is in
[the moving-quadrilateral note](moving_quadrilateral_fifth_point_product.md).
Sections 2--4 give a second proof of qualitative isolation through bounded
triangle determinants and finitely many unimodular affine types. Neither
argument bounds the count in a general endpoint arc.

## 1. Local isolation theorem

Fix `K>0`.  There are constants `c(K)>0` and `R_0(K)>0` with the following
property.  Let four distinct integer points `P_0,P_1,P_2,P_3` lie on a circle
centered at the origin with radius `R>=R_0(K)`.  Suppose they are contained in
one circle arc of length at most

```text
K R^(1/3).                                             (1)
```

Then every integer point `W` on the same circle satisfying

```text
|W-P_0| <= c(K) sqrt(R)                                (2)
```

belongs to `{P_0,P_1,P_2,P_3}`.

The explicit choice in Section 6 works simultaneously for any of the four
anchors.

## 2. The Jarnik-scale span bounds every triangle determinant

Take any three points from the arc and call their lifted angular coordinates
`0<=alpha<=beta<=delta`, after rotating the angular parameter.  Here

```text
delta <= K R^(-2/3).                                  (3)
```

The absolute determinant of the two chords based at the first point is

```text
D=4R^2 sin(alpha/2) sin(beta/2) sin((beta-alpha)/2).
```

Therefore

```text
D <= (R^2/2) alpha beta(beta-alpha)
  <= (R^2/8) delta^3
  <= K^3/8.                                           (4)
```

For the middle inequality, first fix `beta`; the product
`alpha(beta-alpha)` is at most `beta^2/4`, and then `beta<=delta`.
Changing the base vertex only changes the sign or selects the same triangle
area.  Hence every triangle determinant formed from the four points is a
positive integer at most `K^3/8`.  Positivity follows because a line meets a
circle in at most two points.

If `K^3/8<1`, no such four-point configuration exists, so the theorem is
vacuous.  Assume henceforth that `B=floor(K^3/8)>=1`.

## 3. Only finitely many unimodular affine types occur

Translate by `P_0` and write

```text
p_i=P_i-P_0,       1<=i<=3.
```

Choose two of these vectors, say `p_1,p_2`, which are linearly independent,
and put `A=[p_1 p_2]` and `d=|det A|`.  Equation (4) gives `1<=d<=B`.

Choose `S in GL_2(Z)` putting `A` into column Hermite normal form `A'=SA`.
There are only finitely many such matrices with determinant at most `B`.
Moreover

```text
U=d A^(-1)p_3
```

is an integer vector whose two coordinates, up to sign, are
`det(p_3,p_2)` and `det(p_1,p_3)`.  They have modulus at most `B` by (4).
Thus

```text
S p_3=A' U/d
```

also ranges over a finite set.  It follows that the ordered node quadruple

```text
u_0=0, u_1=S p_1, u_2=S p_2, u_3=S p_3              (5)
```

ranges over a finite set depending only on `K`.  Permuting the three nonzero
points to choose an independent pair introduces only finitely many further
possibilities.  In the original coordinates,

```text
P_i=P_0+T u_i,       T=S^(-1) in GL_2(Z).             (6)
```

Thus (1), by itself, extracts a bounded unimodular four-node type.  No
assumption about a fifth point or about saturation of the full difference
lattice is needed.

## 4. Apply fixed-quadrilateral isolation finitely many times

Every type in (5) has no three collinear.  Apply
`fixed_quadrilateral_unimodular_isolation.md` to each of the finitely many
types.  Let `c(K)` be the minimum of their positive isolation constants and
let `R_0(K)` be the maximum of their radius thresholds.  Equation (6) then
shows that every same-circle lattice point satisfying (2) is one of the four
specified points.  This proves the theorem.

If isolation is desired from any of the four anchors, include the finitely
many reanchored and reordered versions of (5) before taking the minimum and
maximum.  This justifies the final sentence of Section 1.

## 5. A five-point gap corollary

Let, for each index `n`, five distinct lattice points lie on a circle of radius
`R_n`, where `R_n` tends to infinity.  Suppose their total arc span is

```text
o(sqrt(R_n)).                                         (7)
```

Then for every choice of four of the five points, their arc span `L_n` obeys

```text
L_n/R_n^(1/3) -> infinity.                            (8)
```

Indeed, if (8) failed for one of the five choices, a subsequence would satisfy
`L_n<=K R_n^(1/3)` for some fixed `K`.  Choose one of those four points as the
anchor.  By (7), the omitted fifth point is at chord distance
`o(sqrt(R_n))` from the anchor, hence eventually at most
`c(K)sqrt(R_n)`.  The local isolation theorem says it must equal one of the
four points, a contradiction.

Because there are only five choices of the omitted point, the conclusion is
simultaneous: the minimum normalized span among all four-point subclusters
tends to infinity.

The hypothesis (7) is essential to this deduction.  The four-point theorem
only excludes an additional point in a fixed endpoint-scale neighborhood of
an anchor; it makes no assertion about a fifth point elsewhere on the circle.

## 6. Exact pencil formula and an explicit separation

There is also a direct quantitative version.  Label the four points in cyclic
order and put

```text
A=D_012, B=D_013, C=D_023, E=D_123,   Pi=A B C E,
```

where `D_ijk` is the positive absolute triangle determinant.  Let
`L_ij(X)=det(P_j-P_i,X-P_i)` with signs chosen consistently, and set
`F=L_01 L_23`, `G=L_02 L_13`.  The circle equation `Q=0`, written in
integer coordinates, lies in this line-pair pencil.  After removing the gcd
of its two pencil coefficients there are coprime positive integers `a,b` and
an integer `s>=1` such that, with consistent signs,

```text
-aF+bG=sQ.                                            (9)
```

For any fifth lattice point `P_4` on the circle, coprimality gives

```text
F(P_4)=k b,       G(P_4)=k a                         (10)
```

for a nonzero integer `k`; it is nonzero because a line meets the circle in
at most two points.  Combining (9)--(10) with the chord formula for each
triangle gives the exact identity

```text
product_(i=0)^3 |P_i-P_4| = 4 |k| s R^2/sqrt(Pi).    (11)
```

Here is a quick check of the normalization in (11).  For four concyclic base
points,

```text
sqrt(Pi)=product_(0<=i<j<=3)|P_i-P_j|/(4R^2),
```

because `D_ijk` is the product of its three side lengths divided by `2R`.
Substituting this and the analogous determinant formulas involving `P_4`
into the two line-pair evaluations in (10) cancels all six base chords and
leaves exactly the factor `4R^2` in (11).

Under (1), equation (4) gives

```text
sqrt(Pi) <= (K^3/8)^2=K^6/64.
```

Write `d=|P_4-P_0|` and `L=K R^(1/3)`.  Every factor on the left of (11) is
at most `d+L`, while `|k|s>=1`.  Consequently

```text
d+K R^(1/3) >= 4 K^(-3/2) sqrt(R).                   (12)
```

Thus one may take, for example,

```text
c(K)=K^(-3/2),
R_0(K)=max(1,K^15/64),                                (13)
```

The small-`K` case in which no four points can exist is vacuous. Indeed (13)
makes `K R^(1/3)<=2K^(-3/2)sqrt(R)`, and (12) forces
`d>=2K^(-3/2)sqrt(R)`, forbidding a distinct fifth point even when
`d<=c(K)sqrt(R)`.  This explicit argument bypasses the finite-type
minimum used in Section 4.
