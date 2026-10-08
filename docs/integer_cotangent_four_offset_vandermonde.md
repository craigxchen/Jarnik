# Four finite cotangents: an integral triangle and a Vandermonde alternative

## Status

Four ordered finite cotangents give a genuinely new simultaneous
divisibility that is absent from the three-cotangent Pell example.  After
anchoring at the minimum, the three normalized pair gaps divide three
shared endpoint cofactors.  Their product defines a positive integer
triangle invariant.

The extra row is needed for the invariant \({\cal J}\), but not for the
following conditional exponent-four bound: multiplying the three exact
pair identities among any chosen three finite cotangents gives it directly.
If
\[
 A=X_0<X_1<X_2<X_3,\qquad
 {\cal V}=\prod_{1\le i<j\le3}(X_j-X_i),
\]
then
\[
 \boxed{\displaystyle
 N\ge
 \frac{\bigl(\prod_{i=1}^3(X_i^2+L^2)\bigr)^{2/3}}
      {4L^2{\cal V}^{\,2/3}}
 \ge\frac{A^4}{4L^2{\cal V}^{\,2/3}}.}                  \tag{1}
\]
Consequently \({\cal V}\le K L^3\) implies
\[
 \boxed{\displaystyle
 N\ge\frac1{4K^{2/3}}\left(\frac AL\right)^4.}           \tag{2}
\]

This does not prove the unrestricted integer height target.  It proves
that any countersequence with
\(N=o((A/L)^4)\) must have
\[
 \boxed{\displaystyle
 \frac{\cal V}{L^3}\ge
 \frac18\left(\frac{(A/L)^4}{N}\right)^{3/2}
 \longrightarrow\infty.}                               \tag{3}
\]
Thus the remaining gap is an ordered-divisor separation problem.  This
Vandermonde alternative is a useful reformulation of the pair constraints,
not a new five-point gain.

## 1. Ordered divisors and endpoint cofactors

Let \(L>0\), and suppose four distinct positive integers
\[
 X_0=A<X_1<X_2<X_3
\]
satisfy
\[
 \frac{X_iX_j+L^2}{X_i-X_j}\in\mathbb Z
 \qquad(i\ne j).                                        \tag{4}
\]
Put
\[
 S=A^2+L^2,\qquad d_i=X_i-A,\qquad
 c_i=\frac S{d_i},\qquad
 e_i=2A+d_i+c_i                                         \tag{5}
\]
for \(i=1,2,3\).  The pair consisting of \(A\) and \(X_i\)
shows that \(d_i\mid S\), so \(c_i,e_i\) are positive integers.  The
quadratic norm at \(X_i\) factors exactly as
\[
 X_i^2+L^2=d_ie_i.                                      \tag{6}
\]

For \(i\ne j\), define
\[
 g_{ij}=\gcd(d_i,d_j),\qquad
 r_{ij}=\frac{|d_i-d_j|}{g_{ij}}.                       \tag{7}
\]
Then the remaining pair condition in (4) is equivalent to
\[
 \boxed{\quad r_{ij}\mid e_i
 \quad\Longleftrightarrow\quad r_{ij}\mid e_j.\quad}     \tag{8}
\]
Indeed, write
\[
 d_i=g_{ij}a,\qquad d_j=g_{ij}b,\qquad
 S=g_{ij}ab h,
\]
where \((a,b)=1\).  Then
\[
 c_i=bh,\qquad c_j=ah,\qquad r_{ij}=|a-b|,
\]
and
\[
 e_i-e_j=(a-b)(g_{ij}-h).                               \tag{9}
\]
Moreover \(X_i-X_j=g_{ij}(a-b)\), while
\[
 X_i^2+L^2=d_ie_i.
\]
Removing the coprime factor \(a\) proves (8).

Equation (8) is stronger with three offsets than with two.  At each vertex,
both normalized incident gaps divide the same integer \(e_i\).  Prime by
prime,
\[
 \max(\alpha_{12},\alpha_{13})
 +\max(\alpha_{12},\alpha_{23})
 +\max(\alpha_{13},\alpha_{23})
 \ge\alpha_{12}+\alpha_{13}+\alpha_{23}.
\]
Therefore
\[
 \boxed{\displaystyle
 {\cal J}:=
 \frac{e_1e_2e_3}{r_{12}r_{13}r_{23}}\in\mathbb Z_{>0}.} \tag{10}
\]
This is the positive integer triangle invariant supplied by the extra
finite coordinate.

## 2. Exact edge cofactors

Orient a pair so that \(d_i>d_j\), and put
\[
 Q_{ij}=\frac{X_iX_j+L^2}{X_i-X_j}>0.
\]
Using (6)--(8) gives the two exact identities
\[
\boxed{\begin{aligned}
 \frac{d_i}{g_{ij}}\frac{e_i}{r_{ij}}&=Q_{ij}+X_i,\\
 \frac{d_j}{g_{ij}}\frac{e_j}{r_{ij}}&=Q_{ij}-X_j.
\end{aligned}}                                          \tag{11}
\]
The second identity in particular shows \(Q_{ij}>X_j\).  These formulas
retain the ordinary gcd \(g_{ij}\); replacing it by one would lose the
scale invariance.

The complement reanchoring at \(A\) replaces every offset by
\[
 d_i\longmapsto c_i=S/d_i.
\]
In the notation above,
\[
 \gcd(c_i,c_j)=h,\qquad
 \frac{|c_i-c_j|}{\gcd(c_i,c_j)}=|a-b|=r_{ij}.           \tag{12}
\]
It also fixes \(e_i=2A+d_i+c_i\).  Hence \({\cal J}\) is exactly invariant
under the complement symmetry.  Applying the argument to the complementary
offsets therefore reproduces the same triangle invariant; it does not
supply a second independent height constraint.

## 3. The elementary three-pair radius product

Let \(N\) be the exact all-edge lcm from
[the integer cotangent height formulation](integer_cotangent_lcm_height_target.md).
For any edge cotangent \(Q\), put \(s=\gcd(Q,L)\).  Its primitive edge norm
divides \(N\), so
\[
 Q^2+L^2
 =s^2\epsilon n_Q
 \le2L^2N.
\]
The same statement applies to every anchor edge \(X_i\).  Thus, with
\[
 B=L\sqrt{2N},
\]
one has
\[
 0<X_i\le B,\qquad 0<Q_{ij}\le B.                       \tag{13}
\]

For an oriented pair \(d_i>d_j\), equations (11)--(13) imply
\[
 \frac{e_i}{r_{ij}}\le
 \frac{2Bg_{ij}}{d_i},\qquad
 \frac{e_j}{r_{ij}}\le
 \frac{Bg_{ij}}{d_j}.                                   \tag{14}
\]
Multiplying (14) over the three edges of the offset triangle packages the
calculation in terms of \({\cal J}\).  The left side is \({\cal J}^2\),
while every \(d_i\) occurs twice on the right.  Therefore
\[
 {\cal J}^2\le
 \frac{8B^6(g_{12}g_{13}g_{23})^2}
      {(d_1d_2d_3)^2}.                                  \tag{15}
\]
Substitute (7) and (10), cancel the complete gcd product, and use
\(B^3=2\sqrt2\,L^3N^{3/2}\).  The result is the integral-free inequality
\[
 d_1d_2d_3e_1e_2e_3
 \le8L^3N^{3/2}
      |d_1-d_2||d_1-d_3||d_2-d_3|.                     \tag{16}
\]
Now (6) turns (16) into
\[
 \prod_{i=1}^3(X_i^2+L^2)
 \le8L^3N^{3/2}{\cal V}.                                \tag{17}
\]
Rearranging proves (1), and (2)--(3) follow immediately.

There is also a shorter derivation which shows precisely why (17) is not
an extra-row consequence.  For any three finite cotangents, each pair
satisfies
\[
 (X_i^2+L^2)(X_j^2+L^2)
 =(X_i-X_j)^2(Q_{ij}^2+L^2)
 \le 2L^2N(X_i-X_j)^2.
\]
Multiplying the three pair inequalities and taking a square root gives
(17) immediately.  The invariant \({\cal J}\) records additional integral
compatibility of four finite cotangents, but the present argument extracts
no stronger upper bound on \({\cal V}/L^3\) and no new small cofactor from
that integrality.

## 4. What the new row does and does not force

The three-finite-cotangent Pell family has only two positive offsets after
anchoring at its minimum.  It has one normalized gap divisibility and no
offset triangle or invariant \({\cal J}\).  It does satisfy the elementary
pair-product estimate whenever three finite coordinates are selected, so
it exposes why (1) alone is not a five-point improvement.

For four finite cotangents, (3) is a necessary condition for failure of the
desired exponent-four height.  What is still missing is an upper bound for
\({\cal V}/L^3\), or another radius contribution that grows when this
normalized Vandermonde grows.  The complement symmetry does not give that
contribution because of the exact invariance (12).

The standard family
\[
 L=6,\qquad X_i=6(t+i)\quad(0\le i\le3)
\]
illustrates the controlled-gap branch: \({\cal V}=2L^3\), so (2) applies
uniformly.  Its radius is actually much larger than the lower bound.  This
example is not an endpoint counterexample.

The local inert-prime bounds and the fixed-spectrum identities do not
contain the simultaneous integrality (10): the former treat one valuation
stratum at a time, while the latter do not control entry height.  Yet the
current use of (10) collapses to the elementary three-pair identity above.
Conversely, (17) gives no cardinality bound when the normalized
Vandermonde is unrestricted.  No unconditional improvement to the general
point-count growth rate is claimed.

The identities and the squared form of (17) are checked in
[the persistent exact checker](check_integer_cotangent_four_offset.py).
