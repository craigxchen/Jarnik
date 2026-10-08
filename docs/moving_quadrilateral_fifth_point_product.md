# An exact fifth-point obstruction for a moving quadrilateral

Let `P_0,P_1,P_2,P_3` be four distinct integer points in counterclockwise
cyclic order on an origin-centered circle of radius `R`. Write

\[
D_{ijk}=\det(P_j-P_i,P_k-P_i),\qquad
A=D_{012},\quad B=D_{013},\quad C=D_{023},\quad E=D_{123}.
\]

All four numbers are positive, and `E=A-B+C`. Set `Pi=ABCE`.
There is an explicitly computable positive integer `s`, depending on
this quadrilateral, such that every other integer point `W` on the circle
satisfies the exact identity

\[
\boxed{\prod_{i=0}^3|W-P_i|
 =\frac{4ks}{\sqrt{\Pi}}R^2,\qquad
k=\gcd\bigl(|D_{01W}D_{23W}|,|D_{02W}D_{13W}|\bigr)\ge1.}
\tag{1}
\]

In particular, if `|W-P_i|<=L` for all four indices, then

\[
\boxed{L\ge\left(\frac{4s}{\sqrt{\Pi}}\right)^{1/4}\sqrt R
\ge\sqrt2\,\Pi^{-1/8}\sqrt R.}\tag{2}
\]

This retains all triangle-shape dependence and has no asymptotic threshold.
It does not bound `Pi` for unrestricted endpoint clusters.

## 1. Primitive pencil parameters, with no unspecified constants

Use translated physical coordinates `x=P-P_0`, and define integer lines

\[
L_{ij}(x)=\det(P_j-P_i,x-(P_i-P_0)),\qquad
F=L_{01}L_{23},\quad G=L_{02}L_{13}.
\]

The circle equation is

\[
Q(x)=|x|^2+2P_0\cdot x=0.
\]

The two line products form a rational basis of the pencil of conics
through the four points. Indeed those four evaluation conditions on the
six quadratic monomials are independent, by products of lines separating
any selected node from the other three. The two products are independent
because their line factorizations differ.

Consequently there are unique coprime positive integers `a,b` and a
positive integer `s` with

\[
\boxed{-aF+bG=sQ,\qquad a>b>0.}\tag{3}
\]

Here is an explicit construction, also proving all signs and integrality.
Put `U=P_1-P_0`, `V=P_2-P_0`, and

\[
M=\gcd\bigl(B|V|^2,C|U|^2\bigr),\qquad
a=\frac{B|V|^2}{M},\quad b=\frac{C|U|^2}{M},\quad
s=\frac{ABC}{M}.\tag{4}
\]

Initially `s` in (4) is rational. In coordinates `x=U\xi+V\eta`, the
fourth node is `(-C/A,B/A)`. The quadratic coefficients of the pencil
therefore give

\[
s|U|^2=bAB,\qquad s|V|^2=aAC,\qquad
s|U-V|^2=A E(a-b).\tag{5}
\]

The first two determine the ratio in (4), and the pencil representation
then gives (3). Since `-aF+bG` has integer coefficients in physical
coordinates and the coefficient of `x_1^2` in `Q` is one, `s` is an
integer. The third equation in (5), together with `A,E,s>0`, proves
`a>b`. In particular, the integer divisibility `M | ABC` is part of the
conclusion, not an assumption in the construction.

## 2. The invariant determinant form and radius identity

If `H_F,H_G` are the quadratic-part matrices, direct expansion gives

\[
\mathcal D(X,Y)=4\det(XH_F+YH_G)
=-(B-A)^2X^2-2(BC+AE)XY-(A+C)^2Y^2.\tag{6}
\]

Its binary discriminant is exactly

\[
\boxed{\operatorname{disc}\mathcal D=16ABCE=16\Pi.}\tag{7}
\]

Thus the primitive circle parameters obey the explicit square equation

\[
\boxed{4s^2=-(B-A)^2a^2+2(BC+AE)ab-(A+C)^2b^2.}\tag{8}
\]

For example, (7) follows from
`(BC+AE)^2-(B-A)^2(A+C)^2=4ABCE`, using `E=A-B+C`.
The three chord formulas (5) and the circumradius identity for triangle
`P_0P_1P_2` give a second exact invariant equation:

\[
\boxed{4s^3R^2=\Pi\,ab(a-b).}\tag{9}
\]

The triangle determinants, and hence this determinant form and its
discriminant, are unchanged under orientation-preserving unimodular
affine changes. These formulas make explicit the moving-shape parameters
that norm-equivalence constants conceal. They do not assert that `s`
or the triangle determinants remain bounded when the shape moves.

## 3. The fifth-point integer and exact distance product

An additional circle point cannot lie on any line through two of the
four base points. Thus `F(W-P_0)` and `G(W-P_0)` are nonzero integers.
Equation (3) and coprimality give a nonzero integer `t` such that

\[
F(W-P_0)=bt,\qquad G(W-P_0)=at,\qquad |t|=k.\tag{10}
\]

For any three circle points, the absolute determinant is their product
of chord lengths divided by `2R`. Hence

\[
|G(W-P_0)|
=\frac{|P_2-P_0|\,|P_3-P_1|}{4R^2}
  \prod_{i=0}^3|W-P_i|.\tag{11}
\]

Besides `s|P_2-P_0|^2=aAC`, equation (3), or its quadratic part applied
to the direction `P_3-P_1`, gives

\[
s|P_3-P_1|^2=aBE.
\]

Therefore the two diagonal lengths have product `a sqrt(Pi)/s`.
Substitute this and `|G|=ak` in (11) and cancel `a`. This proves (1).

For arbitrary proposed fifth triangle data, (10) and (8) give the exact
necessary square condition

\[
\mathcal D(-G(W-P_0)/k,F(W-P_0)/k)=4s^2.
\]

It is a compatibility condition, not a sufficiency assertion for realizing
arbitrary prescribed determinant data by a fifth lattice point.

## 4. Explicit isolation of every sufficiently tight four-point cluster

Suppose the four base points lie in a minor arc of length `S`. Every one
of their triangle determinants is at most `S^3/(8R)`: if the two successive
arc gaps of a triangle are `u,v`, its chord product is at most
`uv(u+v)<=S^3/4`, and the determinant divides this by `2R`. Thus

\[
\sqrt\Pi\le\frac{S^6}{64R^2}.
\]

If an additional point is at distance `d` from `P_0`, all four distances
in (1) are at most `d+S`. Consequently

\[
\boxed{d+S\ge\frac{4s^{1/4}R}{S^{3/2}}
\ge\frac{4R}{S^{3/2}}.}\tag{12}
\]

In particular, for `S<=K R^(1/3)`,

\[
\boxed{d\ge 4K^{-3/2}\sqrt R-KR^{1/3}.}\tag{13}
\]

If `R>K^15/64`, the right side is strictly larger than
`2K^(-3/2) sqrt(R)`. Therefore no fifth lattice point lies within
`2K^(-3/2) sqrt(R)` of the anchor. This is uniform over all four-point
shapes satisfying the displayed arc-span bound, with explicit constants.
It strengthens the quantitative conclusion available by first extracting
a finite family of unimodular shapes. It does not produce such a tight
quadruple inside an arbitrary endpoint cluster.

More generally, if all five points lie in an arc of length `L`, while
the chosen four have span `S`, equation (1) gives

\[
\boxed{S^3L^2\ge16\sqrt s\,R^2\ge16R^2.}\tag{14}
\]

At `L=C_0 sqrt(R)`, every four-point subarc in a cluster containing a
fifth point must therefore have length at least
`(16/C_0^2)^(1/3) R^(1/3)`. This is a useful separation condition; it
does not turn growing point count into a bounded `K` automatically.

## 5. Gaussian content and an exact nested valuation formula

Write the circle points as Gaussian integers `z_i` and put

\[
\mathsf U=(z_1-z_0)(z_3-z_2),\qquad
\mathsf V=(z_2-z_0)(z_3-z_1).
\]

The two products have the same complex argument: each chord argument is
the mean of its endpoint arguments plus `pi/2`, modulo the orientation
sign, and cyclic ordering makes the product signs agree. Their positive
length ratio is `b/a`. Thus `a U=b V`. Coprimality and Bezout give a
Gaussian integer `gamma` with

\[
\mathsf U=b\gamma,\quad \mathsf V=a\gamma,\quad
(z_3-z_0)(z_2-z_1)=(a-b)\gamma.
\]

It is the Gaussian gcd of the first two products, up to a unit. Their
length formulas imply

\[
\boxed{\operatorname{Norm}\gamma=\Pi/s^2,\quad s^2\mid\Pi,
\qquad 4R^2s=ab(a-b)\operatorname{Norm}\gamma.}\tag{15}
\]

In particular, the shape coefficient in (1) is `s/sqrt(Pi)=1/|gamma|`.
The pencil separates an exact common Gaussian factor; it does not supply
an independent size saving beyond that factorization.

Here is a valuation formula retaining arbitrary nested prime levels.
For a split rational prime `p`, choose `pi` above it, set `e=v_p(R^2)`,
and write

\[
t_i=v_\pi(z_i),\quad \delta_{ij}=v_\pi(z_j-z_i),\quad
u=\delta_{01}+\delta_{23},\quad v=\delta_{02}+\delta_{13}.
\]

Then, interpreting the valuation of a positive square root as half the
valuation of its square,

\[
\boxed{v_p(s/\sqrt\Pi)
=\frac12\sum_{i=0}^3t_i-e-\min(u,v).}\tag{16}
\]

Indeed `v_pi(gamma)=min(u,v)`. The common norm gives
`conj(z_j-z_i)=-R^2(z_j-z_i)/(z_i z_j)`, so the two conjugate product
valuations are respectively `u+2e-sum(t_i)` and
`v+2e-sum(t_i)`. Add the two minima to compute `v_p(Norm(gamma))`,
then use (15). For unequal levels `delta_ij=min(t_i,t_j)`; at equal
levels the additional nonnegative residue is retained in `delta_ij`.

For a two-level quartet with levels `0,e`, formula (16) is exactly zero
for a `2+2` split. It is `-e/2-r` for a `3+1` split, where `r` is the
minimum excess among pairs in the three-point level. The ultrametric
inequality ensures that any two of those three pair excesses include
their common minimum. If all four are at the same extreme level, it is
`-e-min(r_01+r_23,r_02+r_13)`. Thus the simple split budgets are recovered,
with all deeper residue contributions explicit. No improvement over a
many-point prime-cut budget follows merely from this re-expression.

## 6. Example, verification, and relation to earlier notes

Take the five circle points

\[
(5,0),\ (4,3),\ (0,5),\ (-3,4),\ (-5,0).
\]

For the first four, `(A,B,C,E)=(10,20,20,10)`, `Pi=40000`, and
`(a,b,s)=(5,1,20)`. The fifth gives `(F,G)=(300,1500)` and `k=300`.
Equation (8) reads `1600=4s^2`; (9) reads `800000=800000`.
The distance product in (1) equals `3000`.

The checker `check_moving_quadrilateral_fifth_point_product.py` verifies
the exact squared forms of (1), (8), and (9), including the construction
of `s` and both line-product divisibilities, on 438,240 configurations
on integer circles with squared radius at most 2000. It also checks the
evaluation determinant below on a deterministic subsample. These finite
checks supplement the proof.

In fact the product identity is exactly the standard five-point
quadratic-evaluation determinant, with its integer content separated by
the pencil. If `V` is the five-by-five matrix whose rows are
`(1,x_i,y_i,x_i^2,x_i y_i)`, then

\[
\boxed{|\det V|=ks
=\frac{\prod_{0\le i<j\le4}|P_i-P_j|}{16R^4}.}\tag{17}
\]

To verify the constant, write `z_i=x_i+iy_i`, so
`conj(z_i)=R^2/z_i`. The change from the basis
`(1,z,conj(z),z^2,conj(z)^2)` to the five displayed real columns has
absolute determinant `1/16`; the constant contribution to `x^2` does
not affect it. Factor `R^2,R^4` from the conjugate columns, then
multiply each row by `z_i^2`. The remaining columns are the powers
`(z^2,z^3,z,z^4,1)`, whose determinant is the Vandermonde determinant
up to sign. Taking absolute values gives (17), since the total row
factor has modulus `R^10`. Finally the product of the six base chords
is `4R^2 sqrt(Pi)`, so (1) identifies this determinant with `ks`.
In the example, `|det V|=6000=300*20`.

Thus the spatial-isolation corollary is a quantitative application of
the familiar five-point determinant, rather than a new independent
arithmetic obstruction. The pencil calculation makes its integer
factorization and moving quadrilateral discriminant explicit.

The earlier `affine_conic_feature_content_audit.md` treats five-point
metric reconstruction and six-point conic compatibility. The pencil and
fifth-point divisibility here encode the same underlying conic geometry;
no independence from those constraints is claimed. The useful conclusion
is the explicit tight-cluster corollary, without norm-equivalence or
boundary limits. If one simply
uses `S=L` in (14), the result is only `L^5>=16R^2`, the usual five-point
exponent `2/5`. Thus eliminating the shape dependence by the generic
triangle upper bound does not solve the endpoint problem.

For background on the classification of four-point clusters at the
`R^(1/3)` scale, see Cilleruelo and Granville,
[Close Lattice Points on Circles (2009)](https://dms.umontreal.ca/~andrew/PDF/closelattice.pdf).
The local separation proved here should not be read as a claim of a new
general counting exponent.
