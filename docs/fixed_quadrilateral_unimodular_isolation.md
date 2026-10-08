# Fixed rational quadrilaterals are isolated at the endpoint scale

**Theorem.** Fix four distinct integer vectors
\(u_0=0,u_1,u_2,u_3\in\mathbb Z^2\), with no three collinear. There are
constants \(c>0\) and \(R_0>0\), depending only on these vectors, with the
following property. Suppose \(T\in\mathrm{GL}_2(\mathbb Z)\),
\(z_0\in\mathbb Z^2\), and the four points
\[
 z_i=z_0+Tu_i
\]
lie on a circle centered at the origin of radius \(R\ge R_0\). Then every
integer point \(w\) on that circle satisfying
\[
 |w-z_0|\le c\sqrt R
\]
is one of the four \(z_i\).

The proof covers every fixed quadrilateral under the stated hypotheses.
For quadrilaterals whose metric pencil has rational boundary rays, radii
are bounded and the conclusion is vacuous above a fixed threshold. On an
unbounded branch the radius has a uniform cubic relation to metric trace,
which is the key estimate below.

This is a fixed-affine-type theorem. No argument here extracts such a
quadrilateral, with uniformly bounded type, from an arbitrary large
endpoint cluster. Consequently the theorem does not establish a uniform
point-count bound for unrestricted circles.

## 1. A fixed integer conic pencil

For each pair of nodes choose an integer linear polynomial defining their
line, and let
\[
 F=\ell_{01}\ell_{23},\qquad G=\ell_{02}\ell_{13}.
\]
These are degree-two integer polynomials, both with constant coefficient
zero, vanishing at the four nodes. They are linearly independent: their
line factorizations differ, since no three nodes are collinear.

The vector space of degree-at-most-two polynomials with constant zero
vanishing at the other three nodes has dimension two. The three evaluation
conditions are independent; for each nonzero node a product of two lines
through the other nodes vanishes there and at zero, but not at the selected
node. Hence every conic through the four nodes and zero constant lies in
\(\operatorname{span}_{\mathbb Q}(F,G)\).

The map from this pencil to its quadratic parts is injective. A pencil
element with zero quadratic part would be a linear polynomial vanishing
at the four nodes, which is impossible unless it is zero. Thus the
quadratic parts span a rational two-dimensional subspace \(V\) of symmetric
matrices, and the associated linear coefficient is a fixed rational linear
function of the quadratic part.

In the given configuration put
\[
 H=T^TT,\qquad l=T^Tz_0,\qquad h=\operatorname{tr}H.
\]
The equation of the original circle pulled back to node coordinates is
\[
 Q(x)=x^THx+2l\cdot x=0.
\]
It vanishes at all four nodes. Moreover \(H\) is positive definite,
\(\det H=1\), and every coefficient of \(Q\) is integral.
There is a fixed positive integer \(D\) such that
\[
 Q=(A F+B G)/D,\qquad A,B\in\mathbb Z.
\]
For example, choose any nonzero two-by-two coefficient minor of the fixed
pair \((F,G)\) and use its absolute determinant as \(D\).

Let \(g=\gcd(|A|,|B|)\), and write \(A=ga,B=gb\). If \(H_F,H_G\)
are the quadratic-part matrices of \(F,G\), their entries are half-integers,
so \(4\det(AH_F+BH_G)\) is an integer binary quadratic form. Since
\(\det H=1\), it equals \(4D^2\). Consequently \(g^2\mid4D^2\), and
\[
 g\le2D,\qquad \gcd(a,b)=1.
\]
The injectivity above and finite-dimensional norm equivalence give fixed
constants with
\[
 |l|\le C_1h,\qquad \max(|a|,|b|)\ge c_1h. \tag{1}
\]
For every additional integer point, its preimage
\(x=T^{-1}(w-z_0)\) is integral and obeys
\[
 aF(x)+bG(x)=0. \tag{2}
\]

## 2. The two metric boundaries and the rational case

If \(V\) contains no positive definite metric, there are no configurations
to consider. Otherwise \(V\cap\{\operatorname{tr}H=1\}\) is an affine
line whose positive definite portion is a bounded open interval. Its two
endpoints \(E_+,E_-\) are distinct positive semidefinite rank-one matrices.
Indeed its direction is a nonzero traceless symmetric matrix, with negative
determinant, so the determinant on that line is a strictly concave quadratic
with two distinct zeros.

The defining equations have rational coefficients. If either boundary is
rational, both are rational, and the determinant binary form on the pencil
factors over \(\mathbb Q\) into independent rational linear factors.
After fixed denominator clearing, its constant nonzero determinant equation
becomes an equation \(L_1(A,B)L_2(A,B)=C\), with integral linear forms and
fixed nonzero integer \(C\). It has only finitely many integer solutions,
since the product has finitely many divisor pairs and the two linear forms
are independent. Thus \(H,l\), and
\[
 R^2=l^TH^{-1}l
\]
range over finite sets. The theorem follows by taking \(R_0\) larger than
these radii.

There is no repeated-root or parabolic boundary case: the presence of a
positive definite metric made the two determinant roots distinct.

## 3. Irrational boundaries force cubic radius growth

It remains to consider irrational \(E_+,E_-\). Let \(l(E)\) denote the
linear coefficient attached to a normalized metric \(E\), using the fixed
linear map from Section 1. At either boundary, if \(n\) spans its null
space, then
\[
 n\cdot l(E)\ne0. \tag{3}
\]
For otherwise \(E=vv^T\), with \(|v|=1\), and \(l(E)=tv\). The limiting
conic factors as
\[
 (v\cdot x)(v\cdot x+2t)=0.
\]
If \(t=0\), all four nodes lie on one line, impossible. Otherwise the
four nodes lie on two parallel lines. No three are collinear, so each line
contains two nodes. The line through zero contains another nonzero integer
node, forcing its direction and the perpendicular direction \(v\) to be
rational. Then the trace-one matrix \(vv^T\) is rational, contrary to the
irrational boundary assumption. This proves (3).

As \(h\to\infty\), normalized metrics \(H/h\) have determinant \(h^{-2}\)
and approach one of these two boundary matrices. The normalized coefficients
\(l/h\) approach their respective \(l(E)\). By (3), the component of \(l\)
along the small-eigenvalue direction has magnitude at least a fixed positive
multiple of \(h\), once \(h\) is large. The large eigenvalue lies between
\(h/2\) and \(h\); the small eigenvalue is its reciprocal. Therefore
\[
 \boxed{c_2h^3\le R^2=l^TH^{-1}l\le C_2h^3} \tag{4}
\]
for all sufficiently large \(h\). The upper bound also follows directly
from \(|l|\le C_1h\) and \(\|H^{-1}\|\le h\).
Both boundaries are fixed, so the constants and the threshold are uniform
across all configurations of this quadrilateral.

## 4. Pulling back a short chord

Let \(J\) be rotation by ninety degrees. Define
\[
 y_0=T^{-1}z_0=H^{-1}l,\qquad
 v_0=T^{-1}Jz_0=(\det T)^{-1}Jl.
\]
The latter identity follows from
\(T^{-1}JT^{-T}=J/\det T\). By (1),
\[
 |y_0|\le C_3h^2,\qquad |v_0|\le C_3h. \tag{5}
\]
Write a second point on the circle as
\(w=z_0\cos\theta+Jz_0\sin\theta\), and put \(d=|w-z_0|\).
The exact chord identities give
\[
 1-\cos\theta=d^2/(2R^2),\qquad |\sin\theta|\le d/R.
\]
Thus if \(d\le c\sqrt R\), equations (4)--(5) imply
\[
 \begin{aligned}
 |x|&=|T^{-1}(w-z_0)|\\
 &\le C_4(c^2h^{1/2}+ch^{1/4}).
 \end{aligned} \tag{6}
\]
For each fixed \(c>0\), since \(F,G\) are fixed quadratic polynomials
with zero constant term, this yields
\[
 \max(|F(x)|,|G(x)|)\le C_5c^4h+o_c(h). \tag{7}
\]
All constants depend only on the original four nodes. The remainder in (7)
accounts for the mixed term of order \(c^3h^{3/4}\), the quadratic term of
order \(c^2h^{1/2}\), and the linear terms in (6).

## 5. Coprime divisibility isolates the four base points

Choose \(c>0\) small enough that \(C_5c^4<c_1/2\). For sufficiently large
\(h\), (1) and (7) imply
\[
 \max(|F(x)|,|G(x)|)<\max(|a|,|b|).
\]
If \(|a|\) is the larger coefficient, (2) and \(\gcd(a,b)=1\) give
\(a\mid G(x)\). The strict size bound forces \(G(x)=0\), and then
\(F(x)=0\). The case in which \(|b|\) is larger is symmetric. The chosen
larger coefficient is nonzero, so this argument also covers a zero smaller
coefficient.

The simultaneous zero set of \(F,G\) consists exactly of the four nodes:
intersecting a line from \(\{\ell_{01},\ell_{23}\}\) with one from
\(\{\ell_{02},\ell_{13}\}\) yields respectively \(u_0,u_1,u_2,u_3\).
No extra intersection or common line occurs, by the no-three-collinear
hypothesis. Thus \(x\) is one of the \(u_i\), proving the theorem for large
\(h\). Since bounded \(h\) gives bounded \(R\) by (1), enlarging \(R_0\)
completes the proof.

In fact the four original points themselves lie within \(O_u(R^{1/3})\)
of \(z_0\), since \(|Tu_i|\le\sqrt h\,|u_i|\) and (4) holds. They are
therefore all inside the asserted endpoint neighborhood for large \(R\).
The conclusion is isolation of this complete four-point cluster.

## Scope of possible extensions

The same proof works when \(T\) is integral and \(1\le|\det T|\le B\)
for a fixed bound \(B\): there are finitely many determinant values, the
coefficient gcd is bounded in terms of \(B\), and the constants in (4)--(6)
can be made uniform over those values. Here \(x\) need not be integral,
but for \(J=\operatorname{lcm}(1,\ldots,B)\), the vector \(Jx\) is integral and
\(J^2F(x),J^2G(x)\) are integers. Apply the coprime divisibility argument
to these two integers and choose \(c\) smaller by a fixed factor depending
on \(J\). This is still a bounded-affine-type
hypothesis. It does not follow from the existence of many lattice points
on an unrestricted endpoint arc.
