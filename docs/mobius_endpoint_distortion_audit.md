# Möbius endpoint distortion and the two-ended content bound

This note isolates the archimedean information supplied by a real rational
projective map and tests it against the determinant-free content theorem.  The
angular calculation is exact, but it does not remove the arithmetic matrix
loss.  A symmetric source/target estimate is valid under symmetric anchored
fixed-unit hypotheses; in general the target anchor and unit class prevent
applying it to the original full configuration.

## 1. Exact angular and chord formulas

Let
\[
 M=\begin{pmatrix}a&b\\c&d\end{pmatrix}\in M_2(\mathbb Z),\qquad
 \delta=\det M\ne0,
\]
be an integer-primitive representative of a rational projective map.  Put
\[
 Q=M^TM,\qquad T=\operatorname {tr}Q,\qquad
 K=T^2-4\delta^2>0.
\]
Write a real direction as
\(u(\theta)=(\cos\theta,\sin\theta)^T\), and let
\(\phi=2\theta\) be the argument of its circle phase
\(u/\bar u=e^{i\phi}\).  If \(\phi'\) is a continuous lift of the
transformed phase argument, direct differentiation gives
\[
 \boxed{\left|\frac{d\phi'}{d\phi}\right|
 =\frac{|\delta|}{u(\theta)^TQ u(\theta)}.}                 \tag{1}
\]
The factor two in the phase arguments cancels.  If \(\sigma_1\ge\sigma_2>0\)
are the singular values of \(M\), and
\(\kappa=\sigma_1/\sigma_2\), then
\[
 \frac1\kappa\le \left|\frac{d\phi'}{d\phi}\right|\le\kappa,
 \qquad
 \kappa+\kappa^{-1}=\frac T{|\delta|},\qquad
 \kappa-\kappa^{-1}=\frac{\sqrt K}{|\delta|}.              \tag{2}
\]
More precisely, on a lifted interval \(I\),
\[
 \frac{|\delta|}{\max_I u^TQu}
 \le \frac{|\phi'(I)|}{|\phi(I)|}
 \le \frac{|\delta|}{\min_I u^TQu}.                       \tag{3}
\]

There is also an exact two-point formula.  For primitive integer columns
\(H_i,H_j\), set
\[
 f_i=N(H_i),\qquad \widetilde F_i=N(MH_i).
\]
For the associated unit-circle phases \(q_i=H_i/\bar H_i\) and their
images \(q_i'\), determinant invariance gives
\[
 \boxed{
 |q_i'-q_j'|
 =\frac{|\delta|\sqrt{f_if_j}}
 {\sqrt{\widetilde F_i\widetilde F_j}}\,|q_i-q_j|.}       \tag{4}
\]
Ordinary coordinate contents cancel from (4), as they must: they change a
chosen integral column but not its projective direction.

## 2. Forward and backward endpoint constants

Let \(N=R^2\) and \(N'=(R')^2\) be the least squared radii of two primitive
integral realizations.  Let \(\Theta,\Theta'<\pi\) be the lengths of their
shortest containing angular arcs, and define their **actual** normalized arc
constants
\[
 \mathcal C=\Theta N^{1/4},\qquad
 \mathcal C'=\Theta'(N')^{1/4}.                            \tag{5}
\]
The restriction below \(\pi\) holds automatically in the endpoint regime
with constant at most two.  The projective circle map and its inverse are
both \(\kappa\)-Lipschitz in angular distance.  Taking the angular diameter
of the finite sets therefore gives
\[
 \boxed{
 \frac1\kappa\left(\frac{N'}N\right)^{1/4}\mathcal C
 \le \mathcal C'
 \le
 \kappa\left(\frac{N'}N\right)^{1/4}\mathcal C.}           \tag{6}
\]
Equivalently, if
\[
 \rho=\frac{\Theta'}{\Theta},
\]
then
\[
 \mathcal C'=\rho\left(\frac{N'}N\right)^{1/4}\mathcal C,
 \qquad \kappa^{-1}\le\rho\le\kappa.                    \tag{7}
\]
Here \(\rho\) is an exact ratio of the two finite-set diameters, not
necessarily the integral in (3) over one common extremal pair, since the
extremal labels can change under the map.

The distinction between an actual constant and an allowed upper bound is
essential.  From hypotheses \(\mathcal C\le C\) and
\(\mathcal C'\le C'\), (6) gives no comparison between \(N\) and \(N'\):
there is no lower bound comparable to \(C\) for the actual source constant.
A configuration can occupy a much shorter subarc than the one permitted by
the statement.  Thus one cannot replace \(\mathcal C,\mathcal C'\) by
\(C,C'\) in (6) in the direction needed for a radius bound.

If the content theorem is stated using chord rather than arc length, the
actual chord constant associated with (5) is
\[
 \widehat{\mathcal C}
 =2N^{1/4}\sin\!\left(\frac{\mathcal C}{2N^{1/4}}\right)
 \le\mathcal C,                                          \tag{8}
\]
and similarly at the target.  Formula (4), rather than a formal replacement
of arc length by chord length, is the exact transfer law.

## 3. What substitution into the content theorem gives

Let
\[
 E=\binom m2,\qquad
 s=2(m-1)+E=\frac{m^2+3m-4}{2},\qquad
 \alpha=\frac{3m-4}{8},
\]
and abbreviate the arithmetic loss by
\[
 D(M)=T^mK^{m-2}.
\]
For the anchored primitive fixed-unit source in the determinant-free theorem,
using an arc bound in place of the slightly sharper chord bound gives
\[
 N'\ge 2^{s-1}\mathcal C^{-s}
 \frac{N^\alpha}{D(M)}.                                  \tag{9}
\]
Substitution of the exact identity (7) merely rewrites this as
\[
 D(M)\ge
 2^{s-1}\rho^s(\mathcal C')^{-s}
 (N')^{s/4-1}N^{-m^2/8}.                                 \tag{10}
\]
The negative source exponent in (10) is the explicit failure of the proposed
archimedean cancellation.  The endpoint relation does not absorb either
\(T\) or \(K\); it has exchanged the useful source endpoint factor for a
generally weaker mixture of the two radii.

There is a legitimate symmetric statement if the target independently
satisfies the same primitive odd fixed-unit hypotheses with the same \(m\)
labels.  Let \(M_{\rm f}\) be the primitive integral representative used in
the forward anchored theorem and let \(M_{\rm rev}\) be the primitive integral
representative used after the target's own anchor normalization in the reverse
theorem.  Put
\[
 D_{\rm f}=T_{\rm f}^mK_{\rm f}^{m-2},\qquad
 D_{\rm rev}=T_{\rm rev}^mK_{\rm rev}^{m-2}.
\]
Applying (9) in both directions gives
\[
 \begin{aligned}
 N'D_{\rm f}&\ge2^{s-1}\mathcal C^{-s}N^\alpha,\\
 ND_{\rm rev}&\ge2^{s-1}(\mathcal C')^{-s}(N')^\alpha,
 \end{aligned}                                          \tag{11}
\]
and hence
\[
 \boxed{
 \sqrt{D_{\rm f}D_{\rm rev}}\ge
 2^{s-1}(\mathcal C\mathcal C')^{-s/2}
 (NN')^{3(m-4)/16}.}                                    \tag{12}
\]
This is a two-ended arithmetic-height obstruction.  It does not improve the
radius exponent.  When \(N'=N\) and the endpoint constants are comparable,
(12) has exactly the same \(N^{\alpha-1}\) power as either one-sided bound.
For a genuine descent \(N'<N\), the forward inequality in (11) is stronger
than their geometric mean, apart from a possible advantage from a much
smaller target endpoint constant.  One may replace the geometric mean in
(12) by one common \(D\) only when the chosen projective representative sends
the distinguished source anchor to the distinguished target anchor, so that
the reverse representative is its primitive adjugate, or when equality of the
two pairs \((T,K)\) is otherwise established explicitly.

## 4. An anchor-sensitive local refinement

Before bounding every transformed norm by \(f_i/T\), the compressed
signed-valuation estimate is
\[
 N'\ge
 \frac{\prod_i F_i^*}
 {2(F_0^*)^2K^{m-2}\prod_{i<j}|b_{ij}|}
 =\frac{\prod_{i>0}F_i^*}
 {2F_0^*K^{m-2}\prod_{i<j}|b_{ij}|},                    \tag{13}
\]
where \(F_i^*\) is the norm after removing the ordinary coordinate content
of \(MH_i\).  Thus even without any local angular hypothesis, using
\(F_i^*\ge f_i/T\) only for \(i>0\) gives the modest improvement
\[
 N'\ge
 \frac{\prod_{i>0}f_i}
 {2F_0^*T^{m-1}K^{m-2}\prod_{i<j}|b_{ij}|}.              \tag{14}
\]
The published \(T^m\) form follows by \(F_0^*\le T\); retaining (14) is
strictly sharper when the transformed anchor has small primitive norm.

For a unimodular map all ordinary contents equal one, and there is a cleaner
short-interval form.  Write
\[
 A=N(Me_1)=F_0,\qquad D=N(Me_2),\qquad
 \vartheta=\max_i|\theta_i|,\qquad
 \eta=\sqrt{D/A}\tan\vartheta.
\]
Assume the source anchor is \(e_1\), all source directions have
\(|\theta_i|<\pi/2\), and \(\eta<1\).  If \(\Theta\) is the phase-arc
width from (5), then \(\vartheta\le\Theta/2\).  The reverse triangle
inequality gives
\[
 \begin{aligned}
 \|Mu(\theta_i)\|
 &\ge \sqrt A\cos\theta_i-\sqrt D\,|\sin\theta_i|\\
 &\ge \sqrt A\cos\vartheta(1-\eta).
 \end{aligned}
\]
Consequently, with
\[
 L_0=\cos^2\vartheta(1-\eta)^2,
\]
one has \(F_i\ge A L_0f_i\) for every nonanchor row.  Substitution into
(13), followed by the same endpoint product estimate as in (9), proves
\[
 \boxed{
 N'\ge 2^{s-1}\mathcal C^{-s}L_0^{m-1}
 N^\alpha\left(\frac A K\right)^{m-2}.}                 \tag{15}
\]
The corresponding local angular statement keeps the anchor stretch exactly:
\[
 \frac1{A(1+\eta)^2}
 \le\left|\frac{d\phi'}{d\phi}\right|
 \le\frac1{A L_0}
 \quad (|\theta|\le\vartheta),\qquad
 \left|\frac{d\phi'}{d\phi}(0)\right|=\frac1A.           \tag{15a}
\]
Indeed the missing upper norm estimate is
\(\|Mu(\theta)\|\le\sqrt A(1+\eta)\).  Thus the ratio of endpoint
widths can determine \(A\), up to the displayed local factors, when the same
lifted interval supplies both widths.

This is the useful local version of the suggested cancellation: when the
stretch is nearly constant on the source arc, the coarse factor
\(T^mK^{m-2}\) is replaced by
\(L_0^{-(m-1)}(K/A)^{m-2}\).

It still does not close the moving-map problem.  If \(\eta\ge1\), the
transverse column already satisfies
\[
 \frac DA\ge\cot^2\vartheta,
\]
which is a matrix-anisotropy lower bound but not a contradiction.  In the
\(\eta<1\) branch, angular distortion near the anchor controls the raw
quantity \(A\), while (15) retains the independent arithmetic discriminant
\(K\).  Endpoint upper bounds do not bound \(K/A\).

For a general nonsingular integral matrix the gap is larger.  If
\[
 r_i=\gcd_{\mathbb Z}((MH_i)_1,(MH_i)_2),
\]
then the geometric stretch is
\[
 \frac{N(MH_i)}{f_i}=\frac{r_i^2F_i^*}{f_i}.
\]
Equations (1)--(7) see the left side but (13) contains \(F_i^*\).  Although
\(r_i\mid\delta\), the individual \(r_i\) can remove exactly the large raw
stretch detected by the angular formula.  The distortion identities alone do
not control these ordinary contents.  The
[joint row-content theorem](mobius_joint_content_conductor.md) does give
partial arithmetic control at determinant primes away from \(2K\), through
the source pair contents, but that information is not contained in either
endpoint width.  Therefore (15) is valid as stated for unimodular maps;
extending this particular local cancellation determinant-free requires
combining it with arithmetic control of the \(r_i/r_0\), not another
distortion inequality.

## 5. Why angle control cannot produce a canonical small matrix

Equations (1)--(7) see only the scale-invariant condition number.  The content
bound sees the absolute arithmetic discriminant of the primitive integral
representative.  These quantities are independent in precisely the direction
needed to defeat the argument.  Take
\[
 M_n=\begin{pmatrix}n+1&0\\0&n\end{pmatrix}.
\]
Then
\[
 \kappa(M_n)=1+\frac1n\longrightarrow1,
\]
while
\[
 T_n=2n^2+2n+1,\qquad K_n=(2n+1)^2,
\]
and therefore
\[
 D(M_n)\asymp_m n^{4m-4}.                               \tag{16}
\]
Thus every angular distance, in both directions and on every arc, is changed
by a factor between \(n/(n+1)\) and \((n+1)/n\), even though the arithmetic
quantity in the determinant-free compressed estimate is unbounded.

This example already retains the distinguished anchor correctly:
\(M_n(1,0)^T=(n+1,0)^T\), whose primitive projective column is again
\((1,0)^T\).  It also survives the allowed output-conformal normalization.
In the notation \(M z=\alpha z+\beta\bar z\), the reduced ratio is
\(\beta/\alpha=1/(2n+1)\), so the minimal output-conformal representative
has
\[
 U=(2n+1)^2,\qquad V=1,
\]
and recovers exactly the displayed \(T_n,K_n\).  Hence neither scalar removal
nor a rational output rotation supplies a bounded-height representative.

The newer joint row-content estimate reduces the effective loss for this
example, so (16) should not be read as the best available content bound.  Here
\(\delta_n=n(n+1)\) and \(\gcd(\delta_n,K_n)=1\).  The ramified extension of the joint-content theorem includes two, so
the full good determinant part equals \(\delta_n\) for every n.  Its compressed gain
\(D_{\rm good}^{\,2m-3}/T_n^2\) changes the effective power from
\(n^{4m-4}\) to \(n^6\).  This is a substantial arithmetic improvement, but
the residual loss is still unbounded while \(\kappa(M_n)\to1\).  The example
therefore establishes only the claimed separation between direct angular
distortion and arithmetic loss; it does not claim optimality of (16).

An input rotation is not available in the anchored theorem.  Replacing
\(M\) by \(LMR\) while keeping the same output configuration requires
replacing the source columns by \(R^{-1}H_i\); the source anchor is then
\(R^{-1}(1,0)^T\), not the unit column used to obtain the nonanchor norm
lower bounds.  Rotating it back undoes \(R\) projectively.  This is why the
smaller representative in a double conformal class cannot be inserted into
(9).

## 6. The inverse-anchor qualification

A rational output rotation \(S\) can send the image of the source anchor to
phase one without changing \(N'\), angular distances, or the image
configuration.  That observation alone does not justify using the original
forward loss \(D(M)\) in the second line of (11).  The reverse map is
\(M^{-1}S^{-1}\), not \(M^{-1}\), and a primitive integral representative of
it need not have the same \(T,K\) as \(M\).  Equivalently, if one replaces the
forward map by \(SM\), then the common-invariant adjugate statement concerns
\(D(SM)\), not \(D(M)\).

There is also a separate unit-class issue.  The target rows may occupy four
different literal Gaussian-unit classes.  The anchored endpoint theorem
requires the retained nonanchor rows and its anchor to come from one suitable
odd fixed-unit class, and a rational output rotation need not preserve that
chosen arithmetic presentation.

One may partition the labels simultaneously by source and target unit class,
choose a joint class, and then choose its anchor.  This costs up to sixteen
classes.  More seriously, discarding the other labels can expose a new common
Gaussian divisor, so the intrinsic squared radii of the two retained
subconfigurations need not be the parent values \(N,N'\).  For a fixed
\(m\)-tuple there is no uniform way to restore those parent radii.  Thus
(12), with the geometric mean
\(\sqrt{D_{\rm f}D_{\rm rev}}\), is a theorem only under the explicitly
symmetric anchored hypotheses.  The stronger same-\(D\) version additionally
requires anchor compatibility or separately verified equality of invariants.
Neither version is an automatic corollary for an arbitrary image of the
original source cluster.

For a large parent configuration, the existing common-unit cohort and
deletion-stability estimates do control this loss.  Combining them with
cross-ratio invariance gives the stronger matrix-independent result already
recorded in `projective_orbit_radius_rigidity.md`:
\[
 \frac{\log R'}{\log R}=1+O_{C,C'}(1/M).
\]
That theorem uses both endpoint scales and permits arbitrary projective matrix
height, but still allows an \(\exp(O(\log R/M))\) radius change and therefore
does not prove a radius-independent bound on the number of lattice points.

## Conclusion

The exact archimedean transfer is (1), equivalently (4) pair by pair, and the
sharp global endpoint comparison is (6).  Combined with determinant-free
content transfer it yields the conditional symmetric obstruction (12) for the
geometric mean of the separately normalized forward and reverse losses, but
no stronger descent exponent.  Identifying those two losses requires the map
to preserve the chosen anchor projectively or an independent invariant check.
The local unimodular refinement (15) is stronger when the map has nearly
constant stretch on the source arc.  The family (16), with the source anchor
retained, shows why endpoint distortion cannot control the integral
invariants \(T,K\) in general, even after the available joint-content gain is
included.  Inverse use additionally requires a target fixed-unit anchor or a
cohort extraction whose radius loss must be charged.  Consequently this route
does not close the uniform `C sqrt(R)` arc problem.
