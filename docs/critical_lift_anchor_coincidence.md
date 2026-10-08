# All actual-anchor discriminant lifts have the same rays and center

Actual reanchoring does not create independent near-square circles or new
radial coordinates. After removing their common ordinary integer content,
all discriminant lifts give one canonical oriented integer tuple and one
rational near-square center. Their radial coordinates are exactly the
original matrix's quadratic stretches. The minimal dilation making the
center integral can improve a triangular lift by at most a small diagonal
factor.

This is an exact covariance result. As in
[critical_shear_near_square_lift.md](critical_shear_near_square_lift.md), its
use for a critical shear is conditional on the admissible-map hypothesis;
no such map is constructed for an arbitrary endpoint cluster.

## 1. Covariance of the discriminant factor

For a nonsingular rational real-linear map write
\[
 Mz=\alpha z+\beta\bar z,\qquad
 \Gamma=-4\alpha\bar\beta,\quad
 T=2(\operatorname{Norm}(\alpha)+\operatorname{Norm}(\beta)),\quad\delta=\operatorname{Norm}(\alpha)-\operatorname{Norm}(\beta).
\]
Assume \(\alpha\beta\ne0\), so \(\operatorname{Norm}(\Gamma)=K=T^2-4\delta^2>0\).
Let the source phases be \(q_i=H_i/\bar H_i\), with \(q_0=1\).
Actual reanchoring at \(j\) replaces them by \(q_i/q_j\), and precomposes
\(M\) with multiplication by \(H_j\). After any allowed nonzero rational
Gaussian output multiplier \(\lambda_j\), the coefficients are
\[
 \alpha_j=\lambda_j\alpha H_j,\qquad
 \beta_j=\lambda_j\beta\bar H_j.
\]
Set \(s_j=\operatorname{Norm}(\lambda_j)\,\operatorname{Norm}(H_j)>0\). Directly,
\[
 \boxed{\Gamma_j=\operatorname{Norm}(\lambda_j)\,\Gamma H_j^2,\qquad
 \Gamma_j(q_i/q_j)=s_j\Gamma q_i.} \tag{1}
\]
The other invariants satisfy
\[
 T_j=s_jT,\qquad\delta_j=s_j\delta,\qquad K_j=s_j^2K. \tag{2}
\]
Ordinary rational scaling is included in \(\lambda_j\), and output
conjugation leaves \(M^TM\) and \(\Gamma\) unchanged. Thus minimal
output-conformal normalization and subsequent triangularization introduce
no missing angular factor into (1).

## 2. One canonical oriented tuple

There is a unique positive rational \(r\) such that
\[
 C_i=r\Gamma q_i\in\mathbb Z[i]
\]
and the gcd of every real and imaginary coordinate of the tuple is one.
Existence follows by clearing rational denominators and then their ordinary
integer gcd; uniqueness follows because a rational multiplier between two
such primitive integer coordinate tuples must equal one.

For any anchor let \(h_j\) be an integer clearing its discriminant lift.
By (1), its tuple is a positive rational multiple of \((C_i)\), hence
\[
 h_j\Gamma_j(q_i/q_j)=n_j C_i,\qquad n_j\in\mathbb Z_{>0}.
\]
Its ordinary coordinate gcd is exactly \(n_j\). Therefore all actual
anchors give the same \((C_i)\) after that removal. In addition (2) gives
\[
 \boxed{h_jT_j/n_j=rT=:R_{\rm can},\qquad
 h_j\delta_j/n_j=r\delta.} \tag{3}
\]
Both the oriented tuple and the rational center are independent of anchor.
The canonical squared radius is the integer
\[
 N_{\rm can}=r^2K=R_{\rm can}^2-(2r\delta)^2.
\]
The tuple can still have a nonunit common Gaussian divisor with relatively
prime ordinary coordinates. Dividing by that divisor rotates the tuple,
and does not preserve this canonical oriented near-axis presentation.

If \(v\) is the positive denominator of \(R_{\rm can}\) in lowest terms,
then \((vC_i)\) is the smallest positive integer dilation of the canonical
tuple having an integer near-square center. It also makes \(2vr\delta\)
an integer: its square is the difference of two integers, and a rational
number whose square is an integer is integral. Every integer-center anchor
lift is an integer dilation of this same minimal tuple, because \(v\mid n_j\).
There are not many different integer centers to intersect; their apparent
differences are exactly these common dilations.

## 3. Size of the removable factor in a triangular lift

For an ordinary-primitive triangular matrix
\(M=\left(\begin{smallmatrix}A&B\\0&D\end{smallmatrix}\right)\),
let \(h\) be the least positive integer clearing all \(h\Gamma q_i\), and
let \(d\) be the ordinary coordinate gcd of that lifted tuple. The
shared-content audit proves \(\gcd(d,h)=1\). Thus the minimal integer-center
lift just described is
\[
 \boxed{(h\Gamma q_i/d_0)_i,\qquad
 d_0=\gcd(d,hT)=\gcd(d,T).} \tag{4}
\]
Its center is \(hT/d_0\); its squared radius is \(h^2K/d_0^2\).

Moreover
\[
 \boxed{d_0\mid2|A|.} \tag{5}
\]
To prove this, \(d\mid\Gamma\) coordinatewise because \(\gcd(d,h)=1\).
Since \(T-\Re\Gamma=2A^2\) and \(\Im\Gamma=2AB\), one has
\(d_0\mid2A^2\) and \(d_0\mid2AB\). Also \(d_0^2\mid K\) and
\(d_0\mid T\), whence \(d_0^2\mid4A^2D^2\), or \(d_0\mid2AD\).
Taking the gcd of these three multiples and using \(\gcd(A,B,D)=1\)
proves (5).

Thus common-content optimization while retaining the integer near-square
center gains at most a small diagonal factor. It does not automatically
remove the \(N^{O(1/m)}\) factors appearing in the strip estimate.

## 4. The common radial levels are exactly quadratic stretches

Write
\(M^TM=\left(\begin{smallmatrix}x&z\\z&y\end{smallmatrix}\right)\).
Then
\[
 \Gamma=y-x+2iz.
\]
For every nonzero Gaussian column \(H\), direct expansion gives
\[
 T-\Re(\Gamma H/\bar H)=2\|MH\|^2/\operatorname{Norm}(H).
\]
Hence the canonical tuple satisfies the exact radial identity
\[
 \boxed{R_{\rm can}-\Re C_i
 =2r\,\frac{\|MH_i\|^2}{\operatorname{Norm}(H_i)}.} \tag{6}
\]
Choosing many anchors simply samples the same quadratic stretches already
present in (6). It cannot create additional radial levels on the canonical
circle. Although \(R_{\rm can}\) can be rational, differences between the
radial levels are integers because the real coordinates of \(C_i\) are
integers.

For a source endpoint cluster with constant at most two, the existing
pair-product rigidity additionally implies that at most one unordered pair
of distinct source points can share a radial coordinate in this oriented
circle. Each such pair has a common product fixed by \(\Gamma\), and two
different pairs would violate that rigidity. Thus there are at least
\(m-1\) distinct radial levels. If their total range is \(H\), this gives
\(m\le H+2\). With the currently available range
\(H=N^{O(1/m)}\), even this observation gives only the familiar
\(m\log m=O(\log N)\), with the input parameter dependence retained.

The exact covariance, shared center, and stretch identity therefore rule
out the proposed gain from treating different anchor lifts as independent
near-square-circle information. Improving the growth rate would require a
new bound on the common stretch range or a different arithmetic input.

The persistent exact verification is
[check_mobius_canonical_discriminant_lift.py](check_mobius_canonical_discriminant_lift.py),
which checks covariance, canonical tuples, rational centers, and (5).


## 5. A necessary strip condition in this aligned family

For the minimal integer-center triangular lift (4), the anchor alone gives
\[
 \mathcal H_{\min}\ge
 \frac{2|AB|\sqrt{h/d_0}}{K^{1/4}}
 \ge\frac{\sqrt{2h|A|}\,|B|}{K^{1/4}}.
\]
The last inequality uses \(d_0\le2|A|\). In the critical regime
\(|B|/K^{1/4}\to1\), a polylogarithmic normalized strip therefore
requires \(h|A|\) to be polylogarithmic. More exactly, Temur's strip
threshold \(6(\log N_{\min})^{\kappa/4}\) requires
\[
 h|A|\le18\frac{K^{1/2}}{B^2}
 (\log N_{\min})^{\kappa/2},
 \qquad N_{\min}=h^2K/d_0^2,
\]
when \(B\ne0\). This is a necessary condition, not a sufficient one:
other points can lie farther from the real axis than the anchor.

This lower bound applies only to the canonical aligned family produced by
this map and its actual-anchor normalizations. An additional rational
rotation, including a different Gaussian multiplier of equal norm, changes
the orientation and requires a separate arithmetic analysis. No global
optimality over all near-square realizations is asserted.

The canonical checker passed 2,205 anchor/output-factor covariance cases
and 3,840 primitive triangular cases, including 1,680 nontrivial removable
factors. The separate
[check_critical_shear_near_square_lift.py](check_critical_shear_near_square_lift.py)
checks the source divisor and shared-content relations on 1,720 cases.
