# Exact conformal quotient normalization of rational Möbius matrices

This note removes artificial matrix height coming from rational output
rotations, and classifies the larger quotient by input and output rotations.
It does not prove a uniform endpoint bound. In particular, the larger quotient
cannot automatically be substituted into a theorem requiring a fixed input
anchor.

For a rational real matrix of nonzero determinant write its real-linear action
on the complex plane as
\[
 Mz=\alpha z+\beta\bar z,
 \qquad
 \alpha=\frac{a+d+i(c-b)}2,
 \quad \beta=\frac{a-d+i(c+b)}2.
\]
Then
\[
 T=\operatorname{tr}(M^TM)=2(N\alpha+N\beta),\quad
 \delta=N\alpha-N\beta,\quad
 K=T^2-4\delta^2=16N\alpha N\beta.
 \tag{1}
\]
Assume throughout that \(K>0\), so both coefficients are nonzero. A rational
orientation-preserving conformal map is multiplication by an element of
\(\mathbb Q(i)^*\). Its induced action on phases is a rational circle
rotation, and hence preserves the least squared realization radius of the
whole configuration.

## 1. Exact minimum in the output-conformal class

Write
\[
 \frac\beta\alpha=\frac vu,
 \qquad u,v\in\mathbb Z[i]\setminus\{0\},\qquad (u,v)=1.
\]
The coprime pair is unique up to a common Gaussian unit. Set
\[
 U=Nu,\qquad V=Nv,\qquad
 \epsilon=\frac4{N\gcd(2,u-v)}\in\{1,2,4\}.
 \tag{2}
\]
Here the gcd is Gaussian. The value of \(\epsilon\) is independent of the
common unit. Since the determinant is nonzero, \(U\ne V\).

**Normalization theorem.** Among all integral matrices obtained from \(M\)
by a rational output conformal map, the simultaneous minima are
\[
 \boxed{T_* =\frac\epsilon2(U+V),\qquad
 K_*=\epsilon^2UV.}
 \tag{3}
\]
They are attained by an integer-primitive matrix. Consequently that matrix
minimizes every product \(T^rK^s\) with \(r,s>0\) in the class.

To prove this, write the coefficients of any matrix in the class as
\[
 2\alpha'=\eta u,\qquad 2\beta'=\eta v,
 \qquad \eta\in\mathbb Q(i)^*.
\]
If the real matrix is integral, both displayed products are Gaussian integers.
Gaussian Bezout for \(u,v\) then implies \(\eta\in\mathbb Z[i]\).
Conversely, for a Gaussian integer \(\eta\), reconstructing the four real
entries shows that integrality is equivalent to
\[
 \eta u\equiv\eta v\pmod{2\mathbb Z[i]}.
 \tag{4}
\]
Thus allowable multipliers form the principal ideal
\[
 \eta\in\frac{2}{\gcd(2,u-v)}\mathbb Z[i].
\]
The least possible nonzero norm is precisely \(\epsilon\). Substituting
\(N\alpha'=N\eta\,U/4\) and \(N\beta'=N\eta\,V/4\) into (1) proves
(3). If the minimizing real matrix had an ordinary common divisor larger
than one, division by it would give an integral matrix in the same class
with smaller trace, a contradiction.

There is an explicit construction: take any Gaussian generator
\(\gamma\) of the ideal in (4), and use
\[
 M_*z=\frac\gamma2(uz+v\bar z).
 \tag{5}
\]
Its four possible minimal choices differ by output quarter-turns. This is
an exact arithmetic normalization, with no constants depending on a fixed
source configuration.

## 2. Safe consequence for the anchored radius bound

Suppose the hypotheses of the determinant-free anchored radius theorem hold:
a primitive fixed-unit \(m\)-cluster of squared radius \(N\), with anchor
\(H_0=1\), endpoint constant \(C\le2\), and all transformed phases retained.
Write
\[
 E=\binom m2,\qquad A_{m,C}=\frac12(2/C)^{2(m-1)+E}.
\]
The theorem proved using
[determinant-free transfer](mobius_determinant_free_transfer_audit.md) gives
\[
 N'\ge A_{m,C}\frac{N^{(3m-4)/8}}{T^mK^{m-2}}.
\]
Output conformal normalization changes neither the source nor \(N'\).
Therefore it gives the sharper, invariant version
\[
 \boxed{
 N'\ge
 A_{m,C}\frac{2^mN^{(3m-4)/8}}
 {\epsilon^{3m-4}(U+V)^m(UV)^{m-2}}.}
 \tag{6}
\]
Let \(B=U+V\). Since \(\epsilon\le4\) and \(UV\le B^2/4\),
\[
 T_*^mK_*^{m-2}\le 2^{3m-4}B^{3m-4}.
\]
Thus for \(m\ge5\), a radius-nonincreasing image requires
\[
 \boxed{B\gg_{m,C}
 N^{3(m-4)/(8(3m-4))}.}
 \tag{7}
\]
The exponent tends to \(1/8\); \(B\) is comparable to the square of the
least matrix height in this output-conformal class. This is stronger in
interpretation than a restriction on an arbitrarily chosen representative:
large output-rotation denominators alone cannot satisfy (7).

The exact bound also retains the improvement for unequal coefficient norms.
Writing \(\theta=4UV/B^2\in(0,1)\), equation (6) is
\[
 N'\ge A_{m,C}(2/\epsilon)^{3m-4}
 \theta^{-(m-2)}\frac{N^{(3m-4)/8}}{B^{3m-4}}.
 \tag{10}
\]
Discarding \(\theta\) can lose an arbitrarily large factor when the map is
close to conformal. All constants in (6) and (10) are explicit.

The invariant \(B\) is nevertheless unbounded over rational maps. Equations
(6)--(7) therefore remain necessary conditions, not a uniform endpoint
obstruction.

## 3. Complete classification under input and output conformal maps

Under output multiplication by \(\ell\) and input multiplication by \(r\),
the coefficient ratio transforms by
\[
 w=\beta/\alpha\longmapsto w\,\bar r/r.
\]
Consequently
\[
 t=Nw=V/U\in\mathbb Q_{>0}\setminus\{1\}
 \tag{8}
\]
is invariant. It is a complete invariant for the double quotient.

Indeed, if \(Nw_1=Nw_2\), then \(h=w_2/w_1\) has norm one. There exists
\(r\in\mathbb Q(i)^*\) with \(\bar r/r=h\): for \(h\ne-1\), use
\(r=1+\bar h\), and for \(h=-1\), use \(r=i\). An output multiplier
then matches the two \(\alpha\)-coefficients. Allowing reflections identifies
\(t\) with \(1/t\). The conformal and anticonformal cases correspond to
the omitted values zero and infinity.

Write \(t=a/b\) in coprime positive integers. Both \(a\) and \(b\) are
sums of two integer squares. For each prime congruent to three modulo four,
the valuation of the rational Gaussian norm \(t\) is even; coprimality
implies the exponents in \(a\) and \(b\) separately are even. The usual
integer sum-of-two-squares criterion proves the assertion.

Choose \(\alpha_0,\beta_0\in\mathbb Z[i]\) with norms \(b,a\).
The integral real matrix of \(\alpha_0z+\beta_0\bar z\) lies in the double
class and has maximum-entry height at most \(\sqrt a+\sqrt b\).
Dividing its ordinary content only improves this height.

Conversely, any integral real matrix in the class with maximum-entry height
\(H\) satisfies
\[
 a/b=N(2\beta)/N(2\alpha),\qquad
 N(2\alpha),N(2\beta)\le8H^2.
\]
Reduction of the fraction yields
\[
 \boxed{
 \sqrt{\max(a,b)/8}\le H_{\mathrm{double},\min}
 \le\sqrt a+\sqrt b\le2\sqrt{\max(a,b)}.}
 \tag{9}
\]
Thus even the full double quotient has unbounded intrinsic arithmetic height,
with uniform absolute comparison constants.

## 4. Two obstructions to using the larger quotient

First, bounded geometric distortion does not bound this arithmetic height.
For
\[
 M_n=\operatorname{diag}(n+1,n),\qquad n\ge1,
\]
one has
\[
 t_n=(2n+1)^{-2},\qquad
 H_{\mathrm{double},\min}\ge(2n+1)/\sqrt8,
\]
while the singular-value ratio \((n+1)/n\) tends to one. Even maps
arbitrarily close to a conformal map can retain unbounded height after both
quotients. The relation between distortion and the invariant is
\[
 K/\delta^2=16t/(1-t)^2.
\]

Second, the input quotient is not free in the anchored estimate. Replacing
\(M\) by \(LMR\) requires replacing source columns by \(R^{-1}H_i\)
to preserve the corresponding output configuration. Its old anchor is now
\(R^{-1}1\), not \(1\). Returning this particular anchor to \(1\) by an
input conformal normalization undoes \(R\), projectively. Alternatively,
keeping the original source columns and using \(LMR\) generally changes
the image configuration and its least radius.

This is material: the anchored radius argument uses \(f_0=1\) and bounds
each nonanchor primitive numerator through its chord to that anchor. A
general input rotation preserves relative phases and \(N\), but it need not
preserve those absolute primitive numerator estimates. Thus (9) cannot be
inserted into (6) without a further theorem handling a moving anchor.

Output normalization (3)--(7) avoids this problem completely. It removes
coordinate inflation that the current proof can legitimately discard; the
remaining invariant still requires new control from the source if this
route is to imply a uniform bound.

## 5. Explicit removal of scalar and conformal padding

Scalar padding \(nM\) fixes every projective phase and makes raw matrix height
arbitrarily large; passing to a primitive integral representative removes it.
Conformal padding persists even among primitive integral matrices. For
\[
 M_0=\begin{pmatrix}2&0\\0&1\end{pmatrix},\qquad
 L_n=\begin{pmatrix}n&-1\\1&n\end{pmatrix},\qquad
 M_n=L_nM_0=\begin{pmatrix}2n&-1\\2&n\end{pmatrix},
 \quad n\ge1,
\]
every \(M_n\) is integer primitive and has height \(2n\). For any source
configuration, its output phases satisfy
\[
 q_i^{(n)}=\frac{n+i}{n-i}\,q_i^{(0)}.
\]
Thus all relative phases and the least realization radius are exactly
unchanged. Nevertheless
\[
 T_n=5(n^2+1),\qquad K_n=9(n^2+1)^2.
\]
The raw denominator \(T_n^mK_n^{m-2}\) is inflated by
\((n^2+1)^{3m-4}\). Our normalization has
\(u=3,v=1,\epsilon=1\), hence \(T_*=5,K_*=9\), independently of \(n\).
This shows why primitiveness alone does not turn raw matrix height into an
intrinsic obstruction, and verifies the role of the post-conformal quotient
on an explicit family.

## 6. Exact finite verification

The root agent's persistent
[checker](check_mobius_conformal_normalization.py) verifies 1,472 coprime
coefficient pairs, including the integrality condition, primitive minimizing
representative, and simultaneous minimum for every tested integral
multiplier. It also verifies 120 conformal-padding cases by reconstructing
the least radius from all pair phases. The proofs above establish the
general statements independently of these finite checks.

## 7. The same normalization optimizes the joint determinant gain

Every integral representative of this output-conformal class has the form
`M_lambda=L_lambda M_*`, with lambda a nonzero Gaussian integer. Indeed
the allowed multiplier ideal is generated by gamma, so eta=gamma lambda.
Writing s=Norm(lambda), its invariants are

`T_lambda=s T_*`, `K_lambda=s^2 K_*`, `delta_lambda=s delta_*`.

Let D_delta denote the part of the absolute determinant coprime to K.
If a prime divides s it also divides K_lambda and contributes nothing to
D_delta(M_lambda). At every other prime the determinant and K valuations
are those of M_*. Thus `D_delta(M_lambda)` divides `D_delta(M_*)`.
Consequently M_* also maximizes the lower bounds with numerator
`D_delta^(2m-3)` and denominators `T^m K^E` or
`T^(m+2) K^(m-2)`, as well as the maximum of the latter with the
original compressed bound. Artificial output rotations cannot improve
these determinant gains. This observation uses the joint-content theorem
in [mobius_joint_content_conductor.md](mobius_joint_content_conductor.md).
