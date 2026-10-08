# Möbius images of a fixed six-node configuration

This note audits the simplest proposed escape from the central-weight height
barrier: apply a parameter-dependent rational Möbius transformation to one fixed
six-node configuration on the rational unit circle. It gives an exact
transformation formula and a scoped obstruction. It does not address moving
six-node families.

Let \(s_1,\ldots,s_6\in\mathbf Q\) be distinct finite cotangent parameters and
write
\[
 q_i=(s_i-i)/(s_i+i).
\]
For a rational Möbius matrix
\[
 M=\begin{pmatrix}a&b\\c&d\end{pmatrix},\qquad ad-bc\ne0,
 \qquad t_i=(as_i+b)/(cs_i+d),
\]
the transformed unit-circle nodes are
\(q_i'=(t_i-i)/(t_i+i)\). Put
\[
 L_i=cs_i+d,\qquad K_i=as_i+b,\qquad
 S_i=K_i^2+L_i^2,\qquad D_i=\prod_{j\ne i}(s_i-s_j).
\]
The central vector for six nodes with finite cotangents is projectively
\[
 w_i=\frac{(s_i^2+1)^2}{D_i}.
\tag{1}
\]
Indeed, the central equations are the conjugate moment system
\(\sum_iw_iq_i^j=\sum_iw_i\overline q_i^{\,j}=0\) for \(j=0,1,2\).
For rational \(w_i\), the conjugate equations follow from the first three;
the resulting real system has rank five. Clearing the fixed factors \(s_i+i\)
in the complex moment equations gives (1).
Since
\[
 t_i-t_j=\frac{(ad-bc)(s_i-s_j)}{L_iL_j},
\]
a direct cancellation gives
\[
 \frac{(t_i^2+1)^2}{\prod_{j\ne i}(t_i-t_j)}
 =
 \frac{\prod_jL_j}{(ad-bc)^5}\,
 \frac{S_i^2}{D_i}.
\tag{2}
\]
Consequently the transformed central vector is exactly the rational vector
\[
 \left(\frac{S_i^2}{D_i}\right)_{i=1}^6
\tag{3}
\]
up to one common scalar. Formula (3) is the quantity that must be primitively
cleared; the Möbius denominator factors do not provide six independent ways
to rescale the rows.

For a fixed base tuple \(s_i\), each \(S_i\) is a quadratic form in
\((a,b,c,d)\); when the matrix entries have maximum degree
\(r\), the maximum degree across the six \(S_i\) is \(2r\), although some
individual evaluations can have lower degree. Thus \(S_i^2/D_i\) has maximum
degree \(4r\). A polynomial or fixed-degree
rational one-parameter matrix family has raw central coefficients of degree
four times the matrix degree. A fixed polynomial common factor can lower this
degree and must be removed explicitly. The content-free hypothesis below
means that, after clearing the fixed denominators \(D_i\) and removing their
fixed gcd, the six specialized polynomials in (3) have no common factor in
\(\mathbf Q[\tau]\); under precisely this hypothesis their primitive
specialization height retains the stated leading degree.

This immediately rules out the intended endpoint escape for such a
content-free fixed-base family. If the transformed cotangents coalesce at a
finite value, the endpoint scale \(X=\min_i t_i\) remains bounded. If they
coalesce at infinity, write the numerator and denominator row degrees as
\(r>s\). Then \(X=O(|\tau|^{r-s})\), whereas (3) has primitive height of
degree \(4r\) in the content-free case (the \(S_i\) have degree \(2r\)).
Hence
\[
 \frac{U}{X}\gg |\tau|^{\,4r-(r-s)}
 =|\tau|^{\,3r+s}\longrightarrow\infty.
\tag{4}
\]
The same conclusion holds after a fixed rational reparametrization of the
parameter. A parameter-dependent gcd could invalidate (4), but producing
such a gcd is precisely a new arithmetic content problem; it cannot be
discarded as a projective row normalization.

## The compressed projective parameter

The six quartics in (3) depend on the matrix only through
\[
 Q=(A,B,C)=(a^2+c^2,\;ab+cd,\;b^2+d^2),
\qquad
 S_i=As_i^2+2Bs_i+C.
\tag{5}
\]
Scaling the matrix scales \(Q\) and all \(S_i^2\) by common factors, so the
primitive object is the primitive integer triple \(Q_0=(A_0,B_0,C_0)\). If
\(Q_0\) is primitive, the coefficient vector of
\[
 (A_0x^2+2B_0x+C_0)^2
\]
has bounded content (at most a fixed power of \(2\)) and height comparable to
\(H(Q_0)^2\). Choose any five of the six distinct \(s_i\); evaluation of
quartics at these five points is an invertible fixed rational matrix. After
clearing its fixed denominators, the gcd of all six evaluated values divides
a fixed determinant times the bounded polynomial content. Therefore there is
a constant \(c_s>0\), depending only on the base tuple, such that the
primitive central vector satisfies
\[
 U\ge c_s H(Q_0)^2.
\tag{6}
\]
The reverse inequality \(U\ll_s H(Q_0)^2\) follows from the same fixed
evaluation matrix and integer clearing, so \(U\asymp_s H(Q_0)^2\). This is
unconditional for the fixed Möbius orbit and does not assume a one-parameter
family. The evaluation matrix here is the fixed weighted matrix
whose \(i\)-th row is \(D_i^{-1}(1,s_i,s_i^2,s_i^3,s_i^4)\), with fixed
denominators cleared once and for all.

The endpoint scale also has an invariant bound once the angular span of the
base nodes is retained. Let \(\delta_s>0\) be the shorter circular angular
diameter of the six \(q_i\). If \(\kappa(M)\) is the condition number of the
real \(2\times2\) matrix \(M\), the angular derivative formula gives
\[
 \left|\frac{d\theta'}{d\theta}\right|
 =\frac{|ad-bc|(1+s^2)}
 {(as+b)^2+(cs+d)^2}
 \ge \kappa(M)^{-1}.
\tag{7}
\]
Consequently an image whose six cotangents are positive and whose minimum is
\(X\) has angular diameter at most \(2\arccot X\), while its diameter is at
least \(\delta_s/\kappa(M)\). To handle the cotangent pole, integrate (7) on
both complementary arcs between a pair attaining \(\delta_s\); the shorter
image arc is still at least \(\delta_s/\kappa(M)\). Hence
\[
 X\le \cot\!\left(\frac{\delta_s}{2\kappa(M)}\right)
 \ll_{\delta_s}\kappa(M).
\tag{8}
\]
The condition number is unchanged by scalar normalization and by output
rotations. In terms of the unscaled Gram matrix \(M^TM\), it is
\(\lambda_{\max}/\sqrt{\det(M^TM)}\); equivalently, use the primitive
\(Q_0\) in this scale-invariant ratio. Combining (6) and (8) gives
\[
 \frac{U}{X}
 \gg_s
 H(Q_0)^2\,
 \frac{\sqrt{\det Q_0}}{\lambda_{\max}(Q_0)}
 \gg_s H(Q_0)\sqrt{\det Q_0}.
\tag{9}
\]
Here \(\det Q_0>0\), and for a rational matrix it is a square after primitive
normalization: after clearing denominators, \(\det(M^TM)=(\det M)^2\), and
dividing the Gram matrix by its entry-content divides \(\det M\). Since
\(\lambda_{\max}(Q_0)\ll H(Q_0)\) and \(\det Q_0\ge1\), (8) also gives
\[
 X\ll_s H(Q_0),\qquad U\gg_s X^2.
\tag{10}
\]
There is also an exact rational version of this estimate. Fix any distinct
base pair and let
\[
 c_{ij}=\frac{4(s_i-s_j)^2}{(1+s_i^2)(1+s_j^2)}>0,\qquad
 J=\det Q_0 .
\]
The squared chord between its images is
\[
 |q_i'-q_j'|^2
 =\frac{4J(s_i-s_j)^2}
 {(A_0s_i^2+2B_0s_i+C_0)(A_0s_j^2+2B_0s_j+C_0)}
 \ge \frac{Jc_{ij}}{(A_0+C_0)^2}.
\]
For positive image cotangents at least \(X\), the same squared chord is at
most \(4/X^2\). Since \(A_0+C_0\le2H(Q_0)\), this proves
\[
 X^2\le\frac{16H(Q_0)^2}{Jc_{ij}},\qquad
 U\ge\frac{c_s c_{ij}}{16}JX^2
\]
whenever \(c_s\) is the constant in (6). This proof avoids any choice of
angular branch and retains the positive integral determinant \(J\).

The angular span in (8) is essential: \(X\) is not invariant under circle
rotations, even though \(Q=M^TM\) is. Thus a fixed Möbius orbit cannot obtain
\(U/X\to0\) through an unbounded common content. A counterfamily must move the
underlying six-node configuration.

## Normalizing three nodes

For three selected distinct base parameters \(p_1,p_2,p_3\), use
\[
 \phi(t)=\frac{(t-p_1)(p_2-p_3)}
 {(t-p_3)(p_2-p_1)}
\]
to send them to \(0,1,\infty\). The other three parameters become the
cross-ratios \(x,y,z\). In the finite chart the relevant nodes are
\(T=(0,1,x,y,z)\), and
\[
 D_t=\prod_{\substack{r\in T\\r\ne t}}(t-r),\qquad
 E_{t,k}=\frac{t^k}{D_t}\quad(0\le k\le4).
\tag{11}
\]
The point at infinity supplies minus the leading-coefficient row separately:
for \(P(x)=A_0x^2+2B_0x+C_0\), its weight is \(-A_0^2\), while the five
finite weights are \(P(t)^2/D_t\). Their sum is zero by the degree-four
Lagrange interpolation identity. For
the finite five rows, let \(L\) be the least positive integer clearing the
denominators of \(E\), and put \(M=LE\). An exact evaluation constant is
\[
 c_{\mathrm{eval}}(x,y,z)
 =\frac{1}{4\,|\det M|\,\|M^{-1}\|_\infty}.
\tag{12}
\]
The factor \(4\) bounds the coefficient content of the square of a primitive
quadratic. Formula (12) is an explicit rational function of \(x,y,z\) after
the denominator-clearing integer \(L\) is specified. Since
\[
 \det E=\pm\left(\prod_{t<r\in T}(t-r)\right)^{-1},
\]
this sufficient lower-bound constant depends explicitly on the cross-ratio
denominators, the Vandermonde of \(0,1,x,y,z\), and the inverse evaluation
matrix. It is not asserted to be optimal. Normalizing three nodes therefore
does not, by this estimate, remove dependence on cross-ratio height.

For example, normalizing the first three entries of
\((1,2,3,4,7,11)\) gives
\((0,1,\infty,-3,-3/2,-5/4)\). Formula (12) gives
\[
 L=315,\qquad
 c_{\mathrm{eval}}=\frac1{9601804800}.
\]
For \((2,3,5,12,17,23)\), the normalized tuple is
\((0,1,\infty,-20/7,-5/2,-7/3)\), with
\[
 L=4950,\qquad
 c_{\mathrm{eval}}=\frac1{4034503242000000}.
\]
These exact values show the normalization constant can become very small
when cross-ratio denominators grow; they do not create a counterexample,
because the estimate remains \(U\gg_s X^2\) for each fixed tuple.

As a finite check, applying all integer matrices with entries in
\([-7,7]\) to the fixed positive tuples
\[
 (2,3,5,12,17,23),\quad (1,2,3,4,7,11),\quad
 (2,4,7,12,19,31)
\]
and retaining six distinct positive images gave minimum primitive
\(U/\min t_i\) values respectively
\[
 27753/2,\qquad 4032,\qquad 11027250/7.
\]
These are diagnostics only, not asymptotic evidence; the reproducible scan is
in [check_mobius_fixed_six_scan.py](check_mobius_fixed_six_scan.py).
The independent exact checker
[check_mobius_gram_compression.py](check_mobius_gram_compression.py)
verifies 1,276 rational matrix images, 540 positive-image chord and height
inequalities, and both displayed normalized evaluation constants, including
the sign of the infinity row.
