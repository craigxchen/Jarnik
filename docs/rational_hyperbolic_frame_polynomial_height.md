# Polynomial-height rational hyperbolic frames

This note records a valid polynomial height bound for the fixed six-dimensional
quadratic space used in the cotangent argument.  It does not claim an optimal
exponent.

Let

\[
 q(x)=\sum_{r=1}^6u_rx_r^2,\qquad 0\ne u_r\in\mathbf Z,
 \qquad U=\max_r|u_r|,
\]

and assume \(\sum_r u_r=0\) and that \(q\) is hyperbolic over \(\mathbf Q\).
Put \(e=(1,\ldots,1)\), and write
\(B(x,y)=\sum_r u_rx_ry_r\), so \(q(x)=B(x,x)\).  The claim is that
there is a rational hyperbolic frame containing \(e\), all of whose reduced
coordinates have height \(U^{O(1)}\), with an implied constant depending only
on the dimension.

## Arithmetic input

We use the following precise form of Cassels' small-zero theorem.  If \(F\) is
an integral quadratic form in \(n\ge 2\) variables, has a nonzero rational
zero, and \(H(F)\) is the maximum absolute value of its integral polynomial
coefficients, then there is a nonzero \(z\in\mathbf Z^n\) with \(F(z)=0\) and
\[
 \|z\|_\infty\le C_n H(F)^{(n-1)/2}.
\tag{C}
\]
Here \(C_n\) depends only on \(n\).  This is Cassels, “Bounds for the least
solutions of homogeneous quadratic equations,” *Proc. Cambridge Philos. Soc.*
**51** (1955), 262–264.  Cambridge's [bibliographic record and 1956
addendum](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/63473834DC9C607E74E60540D3A8115A/S0305004100031625a.pdf/bounds-for-the-least-solutions-of-homogeneous-quadratic-equations.pdf)
record the original volume and pages;
the statement and the exponent \((n-1)/2\) are also reproduced in
[Fukshansky's paper](https://www1.cmc.edu/pages/faculty/lenny/papers/quad.pdf),
§1.  We only use (C) in dimensions four and two.

## The first quotient zero

Let \(E=e^\perp\).  Since \(B(e,e)=\sum u_r=0\), the restriction of \(B\) to
\(E\) has radical exactly \(\mathbf Qe\): the ambient form is nondegenerate,
and \(E\) is the orthogonal complement of an isotropic line.  Hence
\(\bar E=E/\mathbf Qe\) is a regular hyperbolic four-space.

For \(i=1,\ldots,4\), set
\[
 b_i=u_6\mathbf e_i-u_i\mathbf e_6.
\]
These vectors lie in \(E\), and their classes form a basis of \(\bar E\).  In
this basis the integral Gram matrix is
\[
 C_{ij}=B(b_i,b_j)=u_i u_6^2\,\delta_{ij}+u_6u_iu_j.
\tag{1}
\]
Thus \(H(C)\le 2U^3\) (up to the harmless factor 2 used when writing a
quadratic polynomial from a symmetric Gram matrix).  Apply (C) with \(n=4\).
There is a nonzero \(y\in\mathbf Z^4\) with \(y^TCy=0\) and
\[
 \|y\|_\infty\ll U^{9/2}.
\tag{2}
\]
The lift \(Y=\sum_{i=1}^4y_ib_i\) is an integral vector in \(E\), is
\(q\)-isotropic, and has height \(U^{O(1)}\).

## The second quotient zero

The class of \(y\) is isotropic in the regular four-space \(C\).  Choose an
index \(r\) with \((Cy)_r\ne0\).  For every \(j\ne r\), the integral vector
\[
 h_j=(Cy)_r\,\mathbf e_j-(Cy)_j\,\mathbf e_r
\]
lies in \(y^\perp\).  Choose two of these vectors whose classes, together
with \(y\), are independent; call them \(h_a,h_b\).  Such a choice exists
because \(y^\perp\) has dimension three and radical \(\mathbf Qy\).  The
two-by-two Gram matrix
\[
 D=(h_s^TCh_t)_{s,t\in\{a,b\}}
\tag{3}
\]
is regular and isotropic: it represents the regular hyperbolic quotient
\(y^\perp/\mathbf Qy\).  Every entry of \(D\) is bounded by a fixed power of
\(U\), because \(C=U^{O(1)}\), \(y=U^{O(1)}\), and \(h_a,h_b\) are obtained by
one integral matrix multiplication.  Applying (C) with \(n=2\) gives a
nonzero \(z=(z_a,z_b)\in\mathbf Z^2\) with \(z^TDz=0\) and
\(\|z\|_\infty=U^{O(1)}\).  Then
\[
 \widetilde Z=z_ah_a+z_bh_b\in\mathbf Z^4,
 \qquad Z=\sum_{i=1}^4\widetilde Z_i b_i\in\mathbf Z^6
\]
is isotropic, belongs to \(E\), and satisfies \(B(Y,Z)=0\).  Its class is
independent of the class of \(Y\), so \(e,Y,Z\) are independent and span a
maximal totally isotropic three-space.

This recursive quotient step is essential: Cassels' theorem supplies one
small zero at a time; rational linear algebra by itself does not supply a
second independent isotropic vector with controlled height.

## Explicit duals and the complementary half

Let \(V_0=e,V_1=Y,V_2=Z\).  Choose three coordinate positions \(k_0,k_1,k_2\)
for which
\[
 M=(u_{k_s}(V_j)_{k_s})_{0\le j,s\le2}
\]
has nonzero determinant \(\Delta\).  Such a minor exists because the three
linear functionals \(B(V_j,-)\) are independent.  For each \(i=0,1,2\), set
the three selected coordinates of \(G_i\) to the Cramer's-rule solution of
\[
 B(V_j,G_i)=\delta_{ji}\quad(0\le j\le2),
\]
and set the other coordinates to zero.  The numerators are integer minors and
the sole denominator is the nonzero integer \(\Delta\).  Both the numerators
and \(|\Delta|\) are determinants of matrices whose entries already have
height \(U^{O(1)}\), so all \(G_i\) have height \(U^{O(1)}\); no lower bound for
\(|\Delta|\) beyond \(|\Delta|\ge1\) is being hidden.

Put \(A_{ij}=B(G_i,G_j)\), and define
\[
 F_i=G_i-\frac12\sum_{k=0}^2A_{ik}V_k.
\tag{4}
\]
Because the \(V_k\) are mutually orthogonal and isotropic,
\[
 B(V_j,F_i)=\delta_{ji},\qquad B(F_i,F_j)=0.
\]
Therefore
\[
 (V_0,V_1,V_2,F_0,F_1,F_2)
\]
is a rational hyperbolic frame, and every coordinate has height \(U^{O(1)}\).
The entries of \(A\) have denominators dividing \(\Delta^2\), and (4) adds
only the factor 2; thus all denominators are explicitly controlled by
\(2\Delta^2\).

Clearing the frame by the common denominator \(2\Delta^2\) also gives
integer pencil coefficients of polynomial height. The argument proves a
fixed-power bound \(T\le C_0 U^C\) for absolute constants \(C_0,C\) in
dimension six. It does not provide a uniform bound on the
frame height independent of \(U\), and it does not justify obtaining extra
isotropic vectors by linear algebra alone.
