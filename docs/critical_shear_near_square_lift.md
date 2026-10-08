# An exact near-square lift for a critical triangular shear

The proposed lift exists: all source points can be placed on one explicit
near-square circle, with a positive integer multiplier bounded by the
missing clipped conductor. Choosing an integer multiplier preserves the
explicit near-real argument. The resulting near-axis strip is still generally too
wide for Temur's Theorem 6. This construction remains conditional on the
existence of the admissible Möbius map used to obtain the critical triangular
normal form; it is not a reduction available for every arbitrary cluster.

## 1. Data and the near-real discriminant factor

Let the source phases be \(q_i=H_i/\bar H_i\), with actual anchor \(q_0=1\),
and let their least squared realization radius be \(N\). Let
\[
 M=\begin{pmatrix}A&B\\0&D\end{pmatrix},\qquad AD\ne0,
\]
be an integral triangular representative of the given map. Set
\[
 T=A^2+B^2+D^2,\quad K=T^2-4A^2D^2>0,
\]
\[
 U=A+D-iB,\qquad V=A-D+iB.
\]
Its coefficient ratio is \(w=V/U\). Define
\[
 \boxed{\Gamma=-U\bar V=B^2-A^2+D^2+2iAB.} \tag{1}
\]
Then \(\operatorname{Norm}(\Gamma)=K\). In the critical-shear regime
\(A,D=N^{O(1/m)}\), \(|B|=N^{1/4+O(1/m)}\), its real part is positive
and dominates its imaginary part for sufficiently large \(N\) and suitable
fixed large \(m\).

The same-anchor normal form and its conditional hypotheses are established
in [mobius_finite_sample_critical_form.md](mobius_finite_sample_critical_form.md).
The discussion below uses its supplied matrix; it does not assert the
existence of such a map for an arbitrary source.

## 2. The oriented interval contains the clipped source data

At an odd split prime \(p=\pi\bar\pi\), write
\[
 a=v_\pi U,\quad b=v_\pi V,\quad
 c=v_{\bar\pi}U,\quad d=v_{\bar\pi}V.
\]
These are nonnegative integers. The coefficient interval is
\[
 I=[\min(b-a,c-d),\max(b-a,c-d)].
\]
The two valuations of \(\Gamma\) are \(a+d\) and \(b+c\). Therefore the
interval of signed phase valuations that \(\Gamma\) clears is
\[
 J=[-(a+d),b+c].
\]
Both endpoints of \(I\) lie in \(J\), and \(0\in J\): explicitly,
\(b-a\ge-a-d\), \(c-d\ge-a-d\),
\(b-a\le b+c\), and \(c-d\le b+c\). Thus
\[
 \boxed{I\cup\{0\}\subseteq J.} \tag{2}
\]
This is the orientation information missing from norm divisibility alone.
It is independent of cancellation common to \(U,V\); such common factors
only enlarge \(J\).

Write \(t_i=v_\pi(H_i)-v_{\bar\pi}(H_i)\), and
\(s_- =\min_i t_i\), \(s_+=\max_i t_i\). The source range contains zero
because \(q_0=1\). The maximum denominator deficits of the products
\(\Gamma q_i\) at the two Gaussian primes are
\[
 r_\pi=\max(0,-s_- -a-d),\qquad
 r_{\bar\pi}=\max(0,s_+-b-c). \tag{3}
\]
Their sum is the source width outside \(J\). By (2), this is at most the
source width outside \(I\), namely the local exponent of
\(N/N_{\rm clip}\). This remains true if the coefficient interval is
entirely outside the source range.

## 3. Exact common multiplier and one circle containing every point

Let \(h\) be the least positive integer such that all \(h\Gamma q_i\)
are Gaussian integers. Formula (3) gives
\[
 v_p(h)=\max(r_\pi,r_{\bar\pi})
 \le r_\pi+r_{\bar\pi}
 \le v_p(N/N_{\rm clip}).
\]
No other primes occur. At inert primes every rational unit-circle phase
has zero Gaussian valuation. At the ramified prime the same follows from
\(q_i\bar q_i=1\), since conjugation fixes that prime. Since \(\Gamma\)
is integral, it introduces no denominator there. Gaussian unit factors in
any source row are likewise harmless. Consequently
\[
 \boxed{h\mid N/N_{\rm clip}.} \tag{4}
\]

For an explicitly Gaussian version of the clearing operation, take the
Gaussian denominator \(\eta\) with valuations \(r_\pi,r_{\bar\pi}\).
Then all \(\eta\Gamma q_i\) are integral, with squared norm
\(K\operatorname{Norm}(\eta)\), and
\(\operatorname{Norm}(\eta)\mid N/N_{\rm clip}\). Multiplication by
\(\bar\eta\) yields the real scalar \(\operatorname{Norm}(\eta)\), eliminating its argument.
The least integer \(h\) in (4) may improve this scalar, since its local
exponent is the maximum rather than the sum of the two deficits.

The promised lifted tuple is therefore
\[
 \boxed{Z_i=h\Gamma q_i\in\mathbb Z[i],\qquad
 \mathcal N=N(Z_i)=h^2K=(hT)^2-4(hAD)^2.} \tag{5}
\]
Every source point is retained, and distinct phases remain distinct. The
anchor is exactly \(h\Gamma\), whose argument has the controlled near-real
form in (1). If \(z_i\) is the original primitive realization with anchor
\(z_0\), the common multiplier is \(h\Gamma/z_0\). It is a Gaussian
integer by Bezout, because it maps every \(z_i\) to the integral tuple
(5). This explicitly checks that the construction embeds the entire
original tuple, rather than unrelated pairwise lifts. In particular,
\[
 N\mid\mathcal N,\qquad
 \operatorname{Norm}\bigl(\gcd_{\mathbb Z[i]}(Z_0,\ldots,Z_{m-1})\bigr)
 =\mathcal N/N. \tag{5a}
\]
Indeed the source tuple has Gaussian gcd one, so the gcd of its image is
the common multiplier \(h\Gamma/z_0\), up to a unit. This also identifies
the complete content cost of the lift.

The lifted tuple need not be primitive. Removing its common Gaussian gcd
would forfeit the displayed near-square presentation and is neither needed
nor justified for this application: the circle theorem counts all integral
points at the displayed norm \(\mathcal N\).

No assumption \(N\mid K\) is made. The precise remaining factor is the
square \(h^2\), bounded by \((N/N_{\rm clip})^2\).

## 4. Exact angular and radial bounds

Suppose every anchored source phase has an argument \(\phi_i\) with
\(|\phi_i|\le\Delta\). Since \(|\Im\Gamma|=2|AB|\), (5) gives
\[
 |\Im Z_i|\le h(2|AB|+\sqrt K\,\Delta).
\]
Normalize this strip by the fourth root of the lifted norm:
\[
 \boxed{\frac{|\Im Z_i|}{\mathcal N^{1/4}}
 \le\sqrt h\left(\frac{2|AB|}{K^{1/4}}+K^{1/4}\Delta\right)
 =:\mathcal H.} \tag{6}
\]
This bound does not suppress the possibly large scalar \(h\).
For \(\Delta\le C N^{-1/4}\), critical shear sizes, and
\(N/N_{\rm clip}=N^{O(1/m)}\), it gives
\[
 \mathcal N=N^{1+O(1/m)},\quad
 4(hAD)^2=N^{O(1/m)},\quad
 \mathcal H=N^{O(1/m)}.
\]
Any constants depending on \(m\) in the input bounds must also be retained
when \(m\) varies.

For a radial-coordinate version, when the lifted real coordinates are
positive, write \(R_0=hT\). Since
\(R_0^2-\mathcal N=4(hAD)^2\),
\[
 0\le R_0-\Re Z_i
 =\frac{4(hAD)^2+(\Im Z_i)^2}{R_0+\Re Z_i}.
\]
The critical estimates make this radial width \(N^{O(1/m)}\), but do not
make it bounded or polylogarithmic.

## 5. Both hypotheses needed for Temur's result

[Temur, Theorem 6](https://arxiv.org/pdf/2012.10784) applies for sufficiently
large \(\mathcal N\), with \(\kappa=1/2-\varepsilon\) and
\(0<\varepsilon<10^{-2}\), if
\[
 4(hAD)^2\le16\exp(2(\log\mathcal N)^\kappa)
 \quad\text{and}\quad
 \mathcal H\le6(\log\mathcal N)^{\kappa/4}. \tag{7}
\]
Under both conditions it bounds the lifted set by twenty points. These are
two separate requirements. The first concerns the near-square defect; the
second concerns the allowed near-axis strip.

Ignoring displayed parameter-dependent constants only for the purpose of
comparing scales, a bound \(\log(hAD)=O(\log N/m)\) reaches the first
condition around \(m\gg(\log N)^{1-\kappa}\). A strip bound
\(\mathcal H=N^{O(1/m)}\) generally needs the stronger scale
\(m\gtrsim\log N/\log\log N\). Satisfying the defect condition therefore
does not justify invoking the theorem on all lifted points.

The weighted radial estimate in
[Temur, Theorem 3](https://arxiv.org/pdf/2012.10784) instead gives
\(m\ll_\lambda H^\lambda+2\) for radial width \(H\), \(\lambda>1/2\).
With \(H=N^{O(1/m)}\), this yields only
\(m\log m=O(\log N)\), subject to the input constants. It supplies no
uniform point-count conclusion from the present lift. The cited statements
were checked against the primary paper on 2026-09-13.

Thus the exact orientation and common-multiplier obstruction is resolved,
but the available strip width and the conditional map-existence premise
remain substantive limitations. Neither may be discarded in claiming an
endpoint bound.

Independent exact checks covered 1,120 full-tuple lifts built from nested
allocations at the split primes 5, 13, and 17 and triangular matrices with
\(1\le A,D\le4\), \(-4\le B\le4\). They computed the least integer \(h\)
directly from the rational coordinates, verified \(h\mid N/N_{\rm clip}\),
checked every lifted coordinate and its common norm, and checked the single
Gaussian multiplier on the original primitive tuple. All checks passed.

The persistent [exact lift checker](check_critical_shear_near_square_lift.py)
independently passes 1,720 matrix/source cases, including 1,718 nontrivial
integer denominators and 620 cases with missing conductor on both
orientations. It also verifies the exact common-content identity (5a).
