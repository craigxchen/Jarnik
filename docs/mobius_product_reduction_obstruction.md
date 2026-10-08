# Product reduction does not construct a Möbius radius descent

## Status

A natural existence strategy is to minimize
\[
 {\cal P}(M)=\prod_i N(MH_i)
\]
over unimodular integral matrices, using reduction of the associated positive
binary forms.  This note gives an exact obstruction to that strategy.  On an
actual primitive four-point Pell endpoint cluster, the global minimizer over
all \(M\in GL_2(\mathbb Z)\) is a shear at the critical scale
\(B\asymp N^{1/4}\).  It decreases \({\cal P}\) by more than four orders of
magnitude, but increases the least squared radius and sends the image far
outside the endpoint range.

Thus small or minimal row-norm product does not imply either \(N'\le N\) or
a bounded target endpoint constant.  This does not exclude a theorem using
eight or more distinct rows; it shows that such a theorem needs a genuinely
many-row link between the product objective, Gaussian lcm, and target angle.
Binary-form reduction alone supplies no such link.

## 1. An actual endpoint source

Use the primitive Gaussian half-angle rows
\[
\begin{split}
 H_0&=(1,0),\\
 H_1&=(-1241,8),\\
 H_2&=(-2008,-1),\\
 H_3&=(-322,1).
\end{split}                                               \tag{1}
\]
These are the index-one instance of the Pell construction in
[the near-real Gaussian divisor checker](check_gaussian_near_real_divisor_reduction.py).
Direct Gaussian lcm reconstruction gives the least squared radius
\[
 N=358853785
   =5\cdot89\cdot233\cdot3461.                            \tag{2}
\]

The four projective directions occur in the order
\[
 -\frac8{1241}<-\frac1{322}<0<\frac1{2008}.
\]
For the two extreme rows,
\[
 \frac{|\det(H_1,H_2)|}{|\langle H_1,H_2\rangle|}
 =\frac{17305}{2491920}=\frac1{144}.
\]
Consequently their phase diameter is
\[
 \Delta=2\arctan(1/144)<2/144.
\]
Since \(N<144^4\), the actual normalized endpoint constant satisfies
\[
 \Delta N^{1/4}<2.                                       \tag{3}
\]
Thus (1) is an actual \(C<2\) source, rather than an unrestricted set of
rational directions.

Its row-norm product is
\[
 {\cal P}(I)
 =643880195044131125
 =5N^2.                                                   \tag{4}
\]

## 2. Exact global product minimization

Let \(M\in GL_2(\mathbb Z)\).  Put
\[
 p=MH_0,\qquad w=-MH_2,\qquad
 A=N(p),\qquad D=\langle p,w\rangle .
\]
The pair \(H_0,-H_2\) is a unimodular basis, so \(p,w\) is also a
unimodular basis.  Lagrange's identity gives
\[
 N(w)=\frac{D^2+1}{A}.
\]
In the source basis,
\[
 H_1=8(-H_2)-17305H_0,\qquad
 H_3=(-H_2)-2330H_0.
\]
Applying Lagrange's identity to these two combinations yields the exact
two-variable formula
\[
 \boxed{\displaystyle
 {\cal P}(M)=
 \frac{(D^2+1)\bigl((8D-17305A)^2+64\bigr)
                  \bigl((D-2330A)^2+1\bigr)}
      {A^2}.}                                             \tag{5}
\]
This retains the determinant-one constraint: \(A\) is a positive integer
and \(A\mid D^2+1\).

Formula (5) makes the global minimum finite and elementary.  If \(A\ge2\),
write \(r=D/A\).  Split the real line into the following regions:

* If \(r\le1080\), then
  \[
  {\cal P}(M)\ge(8665\cdot1250)^2A^2.
  \]
* If \(1080<r<2247\) and
  \(|r-17305/8|\le20\), then
  \[
  {\cal P}(M)\ge
  64(17145/8)^2(1175/8)^2A^2.
  \]
* In the same middle interval but outside that neighborhood,
  \[
  {\cal P}(M)>
  (1080\cdot160\cdot83)^2A^4.
  \]
* If \(r\ge2247\), then
  \[
  {\cal P}(M)\ge(2247\cdot671)^2A^2.
  \]

The weakest of these comparisons is the last one, and already
\[
 4(2247\cdot671)^2
 =9093083444676
 >8481545624500.                                         \tag{6}
\]
Hence a matrix attaining the value on the right of (6) must have \(A=1\).

When \(A=1\), an integral output rotation sends \(p\) to \((1,0)\), and
\(w=(D,\mathord\pm1)\).  The remaining exact objective is
\[
 F(D)=(D^2+1)\bigl((8D-17305)^2+64\bigr)
                    \bigl((D-2330)^2+1\bigr).             \tag{7}
\]
It decreases toward the finite interval from \(D\le0\) and increases away
from it for \(D\ge2330\).  Exhaustion of the 2,331 integers
\(0\le D\le2330\) gives the unique minimum
\[
 D=2163,\qquad F(2163)=8481545624500.                     \tag{8}
\]
Up to integral output rotations and reflections, the minimizing map is
therefore the shear
\[
 M_{155}=\begin{pmatrix}1&155\\0&1\end{pmatrix}.          \tag{9}
\]
Indeed it sends (1) to
\[
 (1,0),\quad(-1,8),\quad(-2163,-1),\quad(-167,1).
                                                               \tag{10}
\]
Here \(137<N^{1/4}<138\), so \(B=155\) is precisely at the critical
\(N^{1/4}\) shear scale.  Its intrinsic singular-value ratio satisfies
\[
 \kappa+\kappa^{-1}=\operatorname {tr}(M_{155}^TM_{155})
 =155^2+2=24027,
\]
which is consistent with the critical \(N^{1/2}\) geometric scale.
The necessary condition-number theorem for eight or more points does not
apply to this four-point example.

## 3. The minimizing image is neither a radius descent nor endpoint

Gaussian lcm reconstruction of (10) gives
\[
 \begin{split}
 N'&=84815456245\\
   &=5\cdot13^2\cdot17\cdot29\cdot73\cdot2789,
 \end{split}                                               \tag{11}
\]
and
\[
 {\cal P}(M_{155})=100N'.                                 \tag{12}
\]
Thus the globally minimized product is less than the source product, but
\[
 \frac{N'}N
 =\frac{16963091249}{71770757}>236.                       \tag{13}
\]

The angular failure is also exact.  The first two image rows in (10) have
unit-circle phase chord
\[
 |q'_0-q'_1|=\frac{16}{\sqrt{65}}.
\]
If the target lay in an endpoint arc with constant at most two, every pair
chord would be at most \(2(N')^{-1/4}\).  The required inequality fails
because
\[
 4096N'>65^2.                                             \tag{14}
\]

The obstruction is the distinction between a product and an lcm.  The
shear makes one row exceptionally short, which dominates the product
objective, while the other primitive image norms acquire incompatible
Gaussian factors.  Those factors enter the least-radius lcm, and the short
row simultaneously creates a macroscopic projective angle from the anchor.

## 4. Consequence for the existence program

Minimizing \(\prod_iN(MH_i)\), reducing an associated binary form, or choosing
the corresponding Julia/Hermite representative cannot by itself construct
the descent required by the uniform endpoint problem.  Even on an actual
endpoint source, its exact global optimum can have \(N'>N\) and a normalized
target angle far above two.  Unlike the negligible-height maps excluded in
[the earlier radius-inflation audit](projective_radius_inflation.md), the
minimizer here has the critical height \(B\asymp N^{1/4}\).

The example has four distinct points, so it does not refute an additional
theorem specific to \(m\ge8\).  A successful eight-row argument would have
to prove that the extra rows prevent both failures in Section 3.  In
particular, it must control the Gaussian lcm after primitive row reduction;
a bound for the row-norm product, trace, or condition number is
insufficient.

The exact arithmetic and the global-minimum certificate are checked in
[the persistent checker](check_mobius_product_reduction_obstruction.py).
