# All-source-anchor normalization for rational Möbius maps

This note carries out the reanchoring proposed at the end of
[`endpoint_continuation_status.md`](endpoint_continuation_status.md).  The
allowed anchors give an exact geometric-mean form of the radius-content
bound.  Its new denominator is governed by a one-dimensional distance
formula at each split prime.  That formula also gives a sharp obstruction:
prime powers occurring in opposite orientations in the two Möbius
coefficients, but absent from the source conductor, survive with the same
size at every actual source anchor.  Consequently the source hypotheses
alone cannot force one anchor to have bounded normalized loss.  This is
shown on an actual four-point Pell cluster with endpoint constant at most
two.

The countermodel below does not assert that the transformed cluster has
least squared radius at most the source radius.  A result which also uses
that descent hypothesis remains possible; the present result identifies
exactly the extra information such a proof must use.

## 1. Exact normalization at every actual anchor

Let
\[
 H_0=1,H_1,\ldots,H_{m-1}\in\mathbb Z[i]
\]
be the primitive, conjugate-coprime, odd-norm half-angle columns of a
fixed-unit source cluster.  Put
\[
 f_j=N(H_j),\qquad
 e_{ij}=N(\gcd_{\mathbb Z[i]}(H_i,H_j)),\qquad
 b_{ij}=\frac{|\det(H_i,H_j)|}{e_{ij}},
 \quad B_*=\prod_{i<j}b_{ij}.
\]
Changing the source anchor from \(0\) to \(j\) replaces the column for row
\(i\) by the full primitive transition
\[
 H_i^{(j)}=
 \frac{H_i\overline{H_j}}{e_{ij}},\qquad
 f_{ij}:=N(H_i^{(j)})=rac{f_if_j}{e_{ij}^2},             \tag{1}
\]
with \(H_j^{(j)}=1\).  This retains the actual Gaussian common content;
omitting \(e_{ij}\) would give the wrong anchor heights.  The intrinsic
pair residues \(b_{i\ell}\), and hence \(B_*\), are unchanged.

Write a nonconformal real-linear map as
\[
 Mz=\alpha z+\beta\bar z,
 \qquad w=\frac\beta\alpha=\frac vu,\qquad (u,v)=1
 \quad\hbox{in }\mathbb Z[i].
\]
With the convention in the status note, anchor \(j\) changes the coefficient
ratio to
\[
 w_j=w\frac{\overline{H_j}}{H_j}
     =\frac{v\overline{H_j}}{uH_j}.                        \tag{2}
\]
Let
\[
 G_j=\gcd_{\mathbb Z[i]}(v\overline{H_j},uH_j),
 \quad c_j=N(G_j),
 \quad u_j=uH_j/G_j,\quad v_j=v\overline{H_j}/G_j,
 \quad U_j=\frac{N(u)f_j}{c_j},\quad
 V_j=\frac{N(v)f_j}{c_j}.                                 \tag{3}
\]
These are the norms of the reduced denominator and numerator in (2).
Because both \((u,v)=1\) and \((H_j,\overline{H_j})=1\), the complete
cross-content factorization is
\[
 G_j=\gcd(v,H_j)\gcd(u,\overline{H_j})                    \tag{4}
\]
up to a Gaussian unit.  The two factors on the right are coprime.  Formula
(4) includes all prime powers and is often the most convenient way to
compute the permitted cancellation.

Reduce the positive rational norm ratio as
\[
 \frac{N(v)}{N(u)}=\frac ab,\qquad (a,b)=1,\quad a\ne b.
\]
There is therefore a positive integer \(k_j\) such that
\[
 U_j=bk_j,\qquad V_j=ak_j.                                \tag{5}
\]
The parity factor
\[
 \epsilon_j=\frac4{N(\gcd(2,u_j-v_j))}
\]
is independent of \(j\); denote it by \(\epsilon\).  Indeed, \(H_j\) and
\(G_j\) have odd norm and are units modulo \(2\mathbb Z[i]\), while
\(H_j\equiv\overline{H_j}\pmod{2\mathbb Z[i]}\).  Multiplication by these
units preserves the exact divisor of \(2\) in the numerator-denominator
difference.

The simultaneous output-conformal minima at anchor \(j\) are consequently
\[
 \boxed{
 T_j=\tau k_j,\qquad K_j=\kappa k_j^2,
 \qquad
 \tau=\frac\epsilon2(a+b),\quad \kappa=\epsilon^2ab.}     \tag{6}
\]
Also
\[
 |\delta_j|=\frac\epsilon4|a-b|k_j.                       \tag{7}
\]
Let \(D_j\) be the full positive part of \(|\delta_j|\) supported on
rational primes not dividing \(K_j\), including its two-part when
appropriate.  This is the same good-determinant quantity as in
[`mobius_joint_content_conductor.md`](mobius_joint_content_conductor.md);
no replacement by a radical or by a squarefree part is intended.

## 2. The split-prime distance formula

Fix an odd split rational prime \(p=\pi\bar\pi\), and define
\[
 d=v_\pi(w),\qquad e=v_{\bar\pi}(w),
 \qquad t_j=v_\pi(H_j)-v_{\bar\pi}(H_j).
\]
Coordinate primitivity of \(H_j\) means that at most one of the last two
valuations is positive.  From (2), the two valuations of \(w_j\) are
\[
 d-t_j,\qquad e+t_j.
\]
The common rational norm factor in the reduced numerator and denominator is
therefore exactly
\[
 \boxed{
 v_p(k_j)=
 \operatorname {dist}\bigl(t_j,
 [\min(d,-e),\max(d,-e)]\bigr).}                          \tag{8}
\]
To see this directly, if \(t_j\) lies between \(d\) and \(-e\), the two
displayed valuations have the same sign and the reduced numerator and
denominator norms have no common \(p\)-factor.  Outside that interval they
have opposite signs, and the smaller absolute valuation is precisely the
distance to the nearer endpoint.  Odd inert primes are unaffected by
reanchoring, and the prime above two is already accounted for by the common
factor \(\epsilon\).

Thus permitted anchors sample the signed source valuations \(t_j\); they do
not sample arbitrary input rotations.  A useful anchor can remove a
prime-power loss only when one of the actual \(t_j\)'s approaches the
interval in (8).

## 3. The aggregate radius inequality

Set
\[
 P_j=\prod_{i\ne j}f_{ij},
 \qquad
 \Gamma_j=\max\left\{1,\frac{D_j^{\,2m-3}}{T_j^2}\right\}.
\]
Apply the determinant-free compressed content theorem and its joint-content
improvement with anchor \(j\).  Since the transformed anchor is retained,
each of the \(m\) valid inequalities is
\[
 \boxed{
 N'\ge
 \frac{P_j\Gamma_j}
 {2B_*T_j^mK_j^{m-2}}.}                                  \tag{9}
\]
Here \(N'\) is the least squared radius of the same full transformed phase
configuration in every row of (9).  Multiplying all \(m\) inequalities is
therefore legitimate.  Each unordered transition occurs twice, so
\[
 \prod_jP_j=\prod_{i<j}f_{ij}^2
 =\frac{(\prod_jf_j)^{2(m-1)}}{(\prod_{i<j}e_{ij})^4}.     \tag{10}
\]
Taking the geometric mean and using (6) gives the all-anchor form
\[
 \boxed{
 N'\ge
 \frac{\displaystyle
       \prod_{i<j}f_{ij}^{\,2/m}
       \left(\prod_j\Gamma_j\right)^{1/m}}
 {\displaystyle
       2B_*\tau^m\kappa^{m-2}
       \left(\prod_jk_j\right)^{(3m-4)/m}}.}              \tag{11}
\]
Equation (10) is the version that explicitly displays every old Gaussian
common content.  Equation (11) is an actionable reduction of the proposed
route: after all actual anchors are combined, the remaining task is an
upper bound for the product of the distances in (8), possibly aided by the
good-determinant factors \(\Gamma_j\).  Merely proving that the \(t_j\)'s
are spread or averaging their absolute values does not supply such an upper
bound.

## 4. Exact external-prime countermodel on a Pell cluster

Take any sufficiently large member of the primitive common-unit four-point
Pell family in
[`four_point_bonus_counterexample.md`](four_point_bonus_counterexample.md).
Its endpoint constants tend to zero, so fix one member whose four points lie
on an arc of length at most \(2\sqrt R\).  Its squared radius is \(N=R^2\),
and every half-angle norm \(f_j\) divides \(N\).

Choose a rational prime \(p\equiv1\pmod4\) with \(p\nmid10N\), choose
\(p=\pi\bar\pi\), and put \(F=1+2i\).  For every \(L\ge1\), take the
coprime Gaussian coefficient pair
\[
 u_L=\pi^L,\qquad v_L=F\bar\pi^L.                         \tag{12}
\]
It defines a nonconformal rational projective map because
\[
 N(u_L)=p^L,\qquad N(v_L)=5p^L.
\]
The exact conformal-quotient construction produces a primitive integral
real matrix for every such coprime pair, so these are genuine rational
Möbius maps rather than formal coefficient ratios.
Thus \(a=5,b=1\).  Since \(p\) divides none of the \(f_j\), every source
signed valuation at \(p\) is \(t_j=0\).  At the chosen orientation,
\[
 d=-L,\qquad e=L,\qquad [\min(d,-e),\max(d,-e)]=\{-L\}.
\]
Formula (8) now gives the exact simultaneous statement
\[
 p^L\mid k_j\qquad(0\le j<4).                             \tag{13}
\]
No actual source anchor cancels even one copy of this external prime.

Equations (6)--(7) specialize to
\[
 T_j=3\epsilon k_j,\qquad
 K_j=5\epsilon^2k_j^2,\qquad
 |\delta_j|=\epsilon k_j.                                \tag{14}
\]
Every prime dividing \(\delta_j\) therefore divides \(K_j\), so in fact
\[
 D_j=1,\qquad \Gamma_j=1.                                \tag{15}
\]
The complete remaining denominator at every anchor satisfies
\[
 T_j^mK_j^{m-2}
 =3^m5^{m-2}\epsilon^{3m-4}k_j^{3m-4}
 \ge 3^m5^{m-2}p^{L(3m-4)}.                              \tag{16}
\]
For the four-point Pell cluster this grows at least as \(p^{8L}\) at all
four anchors.  The source circle, its full anchor points, its Gaussian
contents, and its endpoint constant remain fixed while \(L\) grows.

This falsifies the target that one actual source anchor must have bounded,
or source-controlled, normalized Möbius loss.  It also shows directly why
the product in (11) has no source-only upper bound: the prime in (12) adds
\(mL\) to \(v_p(\prod_jk_j)\) and contributes nothing to any source factor
in (10) or to any \(D_j\).

The example is not a radius descent: no claim \(N'\le N\) is made.  Such a
hypothesis relates the external Möbius primes to the arithmetic of the
target configuration, and could exclude (12), for example through
three-point projective-height reconstruction.  Hence the viable
continuation is narrower than an all-anchor averaging argument: it must use
the target radius (or equivalent target arithmetic) to bound the external
opposite-orientation part of \(\prod_jk_j\).  The exact quantity that needs
that new input is displayed in (8) and (11).

## 5. The external part is charged to every target pair

There is a useful two-ended repair of the obstruction.  At an odd split
prime retain the interval
\[
 I_p=[\min(d,-e),\max(d,-e)]
\]
from (8), and let
\[
 J_p=[\min_jt_j,\max_jt_j],\qquad
 g_p=\operatorname {dist}(I_p,J_p),\qquad
 G_{\rm ext}=\prod_{p\ {\rm split}}p^{g_p}.               \tag{17}
\]
Thus \(G_{\rm ext}\) is the common part of all the \(k_j\)'s caused by
coefficient intervals lying completely outside the source valuation range.
Primes whose intervals meet \(J_p\) contribute zero.  Then
\[
 \boxed{G_{\rm ext}\mid b'_{i\ell}\quad\hbox{for every }i<\ell,} \tag{18}
\]
where \(b'_{i\ell}\) is the full primitive pair residue of the target
configuration.

Here is the local proof, including row contents and Gaussian pair contents.
Work over \(\mathbb Z_p\) in isotropic coordinates, and write the four
coefficient valuations as
\[
 A=v_\pi(\alpha),\quad B=v_\pi(\beta),\quad
 \bar B=v_\pi(\bar\beta),\quad \bar A=v_\pi(\bar\alpha).
\]
Thus \(d=B-A\) and \(e=\bar B-\bar A\).  Remove the common rational
matrix content, so the minimum of these four nonnegative valuations is
zero.

Suppose first that \(J_p\) lies strictly to the right of \(I_p\).  For a
source row of signed valuation \(t_i\), the two raw output-coordinate
valuations are
\[
 B+(-t_i)_+,\qquad \bar A+(-t_i)_+.
\]
There is no cancellation here because \(t_i>d\) and \(t_i>-e\).  If
\(s=\min(B,\bar A)\), the ordinary row content has exact valuation
\[
 v_p(r_i)=s+(-t_i)_+,
\]
and every primitive output row has the same exact signed Gaussian valuation
\(B-\bar A\).  Hence the Gaussian gcd of any two primitive output rows has
norm valuation \(|B-\bar A|\), with no unrecorded common factor.

The determinant identity
\[
 \det(MH_i,MH_\ell)=\delta\det(H_i,H_\ell)
\]
now gives
\[
\begin{aligned}
 v_p(b'_{i\ell})
 &\ge v_p(\delta)-2s-|B-\bar A|\\
 &\quad
   +v_p(e_{i\ell})-(-t_i)_+-(-t_\ell)_+\\
 &=v_p(\delta)-(B+\bar A)+\min(t_i,t_\ell).       \tag{19}
\end{aligned}
\]
The last identity uses
\[
 v_p(e_{i\ell})-(-t_i)_+-(-t_\ell)_+
 =\min(t_i,t_\ell);
\]
it holds for equal or opposite source orientations.  Finally,
\[
 v_p(\delta)\ge\min(A+\bar A,B+\bar B),
\]
so
\[
 v_p(\delta)-(B+\bar A)
 \ge\min(-d,e)=-\max(d,-e).
\]
Equation (19) is therefore at least
\[
 \min_jt_j-\max I_p=g_p.
\]
The case in which \(J_p\) lies to the left of \(I_p\) is the conjugate
calculation and gives \(\min I_p-\max_jt_j=g_p\).  This proves (18).
The proof used the actual ordinary row contents \(r_i\), the full old
contents \(e_{i\ell}\), and the full target Gaussian gcds.

Suppose now that the target itself is a primitive \(m\)-point configuration
of least squared radius \(N'\), contained in an arc of length at most
\(C'(N')^{1/4}\).  No common-unit hypothesis is needed.  Let
\[
 d'_{i\ell}=N\left(
 \frac{B_i\overline{B_\ell}}
 {N(\gcd_{\mathbb Z[i]}(B_i,B_\ell))}
 \right)
\]
be the primitive target pair norm, and let \(r\) of the \(m\) primitive
target rows have two odd coordinates.  Prime-power layer cake at every odd
split prime, followed by the binary parity cut at two, gives
\[
 \prod_{i<\ell}d'_{i\ell}
 \le 2^{r(m-r)}(N')^{m^2/4}.                              \tag{20}
\]
At two, the pair norm has valuation one exactly when one row is odd-odd and
the other is not, explaining the exact exponent \(r(m-r)\).

The target chord estimate, pair by pair, is
\[
 b'_{i\ell}\le
 \frac{C'}2\sqrt{d'_{i\ell}}\,(N')^{-1/4}.
\]
Writing \(E=\binom m2\), multiplying this estimate and using (20) yields
\[
 B'_*:=\prod_{i<\ell}b'_{i\ell}
 \le
 \left(\frac{C'}2\right)^E
 2^{r(m-r)/2}(N')^{m/8}.                                 \tag{21}
\]
Since (18) gives \(G_{\rm ext}^E\mid B'_*\), one obtains the explicit
two-ended bound
\[
 \boxed{
 G_{\rm ext}\le
 \frac{C'}2\,
 2^{\,r(m-r)/(2E)}
 (N')^{1/(4(m-1))}.}                                    \tag{22}
\]
In particular \(r(m-r)/(2E)\le m/(4(m-1))\).

Thus the external opposite-orientation factor in the Pell countermodel
cannot remain arbitrary if its image is also required to be an endpoint
cluster of controlled radius.  Formula (22) handles exactly the part of
\(\prod_jk_j\) coming from intervals disjoint from the source valuation
range.  The unresolved part comes from primes for which \(I_p\) intersects
\(J_p\): their individual distances in (8) can still have a large product,
and neither (20) nor source-only averaging presently bounds it sharply
enough to close radius minimality or the uniform point count.
