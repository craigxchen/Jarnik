# Shared content of the critical lift: exact identities and limits

The common Gaussian content of the near-square lift supplies an exact
small-cofactor factorization of the source conductor. It does not by itself
shrink the normalized near-axis strip. Removing that content can destroy
the displayed integer near-square center, and Gaussian lattice index does
not give a corresponding radial spacing gain. The statements below remain
conditional on the admissible map and critical triangular normal form from
[mobius_finite_sample_critical_form.md](mobius_finite_sample_critical_form.md).

Retain the notation of
[critical_shear_near_square_lift.md](critical_shear_near_square_lift.md):
\[
 M=\begin{pmatrix}A&B\\0&D\end{pmatrix},\quad
 T=A^2+B^2+D^2,\quad K=T^2-4A^2D^2,
\]
\[
 \Gamma=B^2-A^2+D^2+2iAB,\qquad \operatorname{Norm}(\Gamma)=K.
\]
Let the original primitive Gaussian tuple be \(z_i\), with common norm
\(N\), anchor \(z_0\), and \(q_i=z_i/z_0\). The least positive integer
\(h\) making all \(h\Gamma q_i\) integral satisfies
\(h\mid N/N_{\rm clip}\).

## 1. Exact common content and a coprimality consequence of minimality

The common multiplier
\[
 L=h\Gamma/z_0\in\mathbb Z[i]
\]
exists by primitive-source Bezout. The lifted tuple is precisely \(Lz_i\),
so its Gaussian gcd is \(L\), up to a unit. Hence
\[
 \boxed{g=\operatorname{Norm}(L)=h^2K/N\in\mathbb Z_{>0},\qquad N\mid h^2K.} \tag{1}
\]
Writing
\(e=N/N_{\rm clip}\) and \(k=K/N_{\rm clip}\), this becomes
\[
 \boxed{g=h^2k/e,\qquad h\mid e,\qquad e\mid h^2k.} \tag{2}
\]
In particular \(g\le ek\), so the critical estimates make this Gaussian
content norm \(N^{O(1/m)}\).

There is a further constraint from the least-integer choice of \(h\).
Let \(d=\gcd(|\Re L|,|\Im L|)\). If a rational prime divided both \(d\)
and \(h\), dividing the entire lifted tuple by that prime would contradict
minimality of \(h\). Therefore
\[
 \boxed{\gcd(d,h)=1.} \tag{3}
\]
Since \(Lz_0=h\Gamma\), the integer \(d\) divides both coordinates of
\(h\Gamma\); (3) gives
\[
 \boxed{d\mid\gcd(|B^2-A^2+D^2|,2|AB|).} \tag{4}
\]
At primes dividing \(h\), the common Gaussian content can occupy only one
of the two conjugate orientations. It is not necessarily coprime to \(h\)
as a Gaussian integer, and (3) must not be strengthened to that assertion.

## 2. A small-cofactor divisor relation inside the source conductor

Define
\[
 K_+=B^2+(A+D)^2,\qquad K_-=B^2+(A-D)^2,
\]
so \(K=K_+K_-\) and \(K_+-K_-=4AD\). Put
\[
 n_+=\gcd(N_{\rm clip},K_+),\qquad
 n_-=N_{\rm clip}/n_+.
\]
Since \(N_{\rm clip}\mid K_+K_-\), primewise comparison proves
\(n_-\mid K_-\). Thus the positive integer cofactors
\[
 u_+=K_+/n_+,\qquad u_-=K_-/n_-
\]
satisfy the exact identities
\[
 \boxed{N=e n_+n_-,\qquad u_+u_-=k,\qquad
 u_+n_+-u_-n_-=4AD.} \tag{5}
\]
Moreover
\[
 \boxed{\gcd(n_+,n_-)\mid4AD.} \tag{6}
\]
No coprimality between \(K_+,K_-\) is assumed. In the critical regime,
\(K_\pm=N^{1/2+O(1/m)}\) and \(e,k=N^{O(1/m)}\); hence
\[
 n_\pm=N^{1/2+O(1/m)},\qquad e,u_+,u_-=N^{O(1/m)}.
\]
This is a precise source-arithmetic consequence: almost all of its conductor
splits into two nearly balanced divisors satisfying a linear relation with
small integer coefficients and small nonzero residual \(4AD\).

Its parameters still have polynomial size \(N^{O(1/m)}\). The relation
therefore does not turn the existing near-square or strip error into a
bounded or polylogarithmic quantity. Multiplying (5) back together recovers
\(kN_{\rm clip}=T^2-4A^2D^2\), so an application based solely on that
quadratic identity is the same near-square condition already present in the
lift, with all cofactors retained.

## 3. Dividing out content does not preserve the integer near-square center

The lifted norm has integer center \(hT\):
\[
 gN=(hT)^2-(2hAD)^2.
\]
Removing the full Gaussian content \(L\) restores exactly the original
source tuple and its original angular placement. Even dividing just by the
positive integer content \(d\) replaces the center by \(hT/d\), which need
not be an integer. Relation (4) does not imply \(d\mid T\).

For a concrete exact example, take
\[
 M=\begin{pmatrix}1&3\\0&2\end{pmatrix},\quad
 T=14,\quad K=180,\quad\Gamma=12+6i,
\]
and the primitive source pair \(z_0=2-i\), \(z_1=2+i\), of norm five.
Then the least integer denominator is \(h=5\), and
\[
 L=h\Gamma/z_0=18+24i=6(3+4i).
\]
The lifted points are \(60+30i\) and \(12+66i\), of norm \(4500\).
Their integer content is six, coprime to \(h\), but \(6\nmid hT=70\).
After dividing by six the circle norm is 125, and the inherited center is
\(70/6\), not an integer. This example is an algebraic obstruction to the
content-removal step, not an assertion of a large critical endpoint tuple.

## 4. Gaussian content does not force radial gaps of its modulus or norm

There is also no general lattice-spacing gain from the index
\(g=\operatorname{Norm}(L)\).
For integers \(n,t\ge1\), set
\[
 L_n=n-(n+1)i,\quad
 z_t=(t+1)+ti,\quad z_t'=t+(t+1)i.
\]
The source pair has equal norm \(2t^2+2t+1\) and Gaussian gcd one: a common
divisor would divide \(-1+i\), but both source norms are odd. The two
lifted points have common Gaussian gcd exactly \(L_n\), whose norm
\(2n^2+2n+1\) is unbounded. Nevertheless
\[
 \Re(L_nz_t')-\Re(L_nz_t)=1.
\]
Explicitly the lifted coordinates are
\[
 ((2n+1)t+n,\ -t-n-1),\qquad
 ((2n+1)t+n+1,\ n-t).
\]
Thus arbitrarily large primitive Gaussian common content is compatible with
a radial-coordinate gap of one. For \(n\to\infty\) and \(t\gg n\), the
lifted direction approaches the real axis while its strip factor relative
to the fourth root of its squared norm has order \(\sqrt{t/n}\), which
can grow. This example addresses the proposed spacing inference only; it
does not claim all the critical-shear hypotheses simultaneously.

## 5. Changing the actual anchor gives the same canonical lift

For a general nonsingular nonconformal real matrix write
\(M(z)=\alpha z+\beta\bar z\), and set
\[
 \Gamma_M=-4\alpha\bar\beta,\quad
 T_M=\operatorname{tr}(M^TM),\quad \delta_M=\det M.
\]
For integral M, the factor \(\Gamma_M\) is a Gaussian integer of norm
\(T_M^2-4\delta_M^2\). If
\(M^TM=\left(\begin{smallmatrix}x&z\\z&y\end{smallmatrix}\right)\),
then \(\Gamma_M=y-x+2iz\).

Reanchor the source at the actual phase \(q_j=H_j/\bar H_j\), so its
new phases are \(q_i/q_j\). The associated input precomposition is
multiplication by \(H_j\). If an output conformal factor \(\lambda\)
and a rational scalar normalization are also used, the resulting
coefficients have the form
\[
 \alpha_j=\lambda\alpha H_j,\qquad
 \beta_j=\lambda\beta\bar H_j.
\]
The rational scalar can be absorbed in \(\lambda\). Consequently,
with the same positive rational factor \(s_j\),
\[
 \Gamma_{M_j}(q_i/q_j)=s_j\Gamma_M q_i,\qquad
 T_{M_j}=s_jT_M,\qquad |\delta_{M_j}|=s_j|\delta_M|,
 \quad s_j=\operatorname{Norm}(\lambda)\operatorname{Norm}(H_j).
 \tag{7}
\]
An output reflection preserves \(\Gamma_M\) and T and changes only the
sign of the determinant, so it gives the same unsigned statement.
Primitive reduction of individual half-angle rows does not change their
phases and therefore does not affect (7).

Let d be the common ordinary coordinate gcd of all \(h\Gamma_M q_i\).
It equals the ordinary coordinate content of L used above: each divisibility
direction follows by multiplying the primitive source tuple and then using
its Gaussian Bezout identity. Define
\[
 C_i=(h/d)\Gamma_M q_i\in\mathbb Z[i],\quad
 R_{\rm can}=hT_M/d,\quad E_{\rm can}=2h|\delta_M|/d.
 \tag{8}
\]
The tuple \((C_i)\) has common ordinary coordinate gcd one. Every positive
rational multiple of it with integral coordinates is a positive integer
multiple, by ordinary Bezout applied to all its coordinates. Thus (7)
shows that all actual-anchor lifts give the same ordered canonical tuple,
the same rational center \(R_{\rm can}\), and the same rational defect
root \(E_{\rm can}\), after removing their ordinary content. In particular
\[
 \operatorname{Norm}(C_i)=R_{\rm can}^2-E_{\rm can}^2.
 \tag{9}
\]
These are not independent near-square circles produced by different anchors.

There is an exact interpretation of their radial coordinates. For every
nonzero real direction represented by H,
\[
 T_M-\Re(\Gamma_M H/\bar H)
 =2\frac{\|MH\|^2}{\|H\|^2}.
\]
Hence the canonical radial offsets are exactly the row stretches:
\[
 R_{\rm can}-\Re C_i
 =\frac{2h}{d}\frac{\|MH_i\|^2}{\|H_i\|^2}. \tag{10}
\]
Reanchoring does not provide an additional collection of offsets beyond
these same coordinates.

## 6. The least aligned lift retaining the displayed integer center

The least integer dilation of \((C_i)\) for which its inherited center
is an integer is
\[
 \widetilde C_i=(h/d_0)\Gamma q_i,\qquad
 d_0=\gcd(d,hT)=\gcd(d,T), \tag{11}
\]
where the second equality uses \(\gcd(d,h)=1\). Indeed the denominator
of \(hT/d\) in lowest terms is \(d/d_0\).

If the triangular matrix is ordinary primitive,
\(\gcd(A,B,D)=1\), then
\[
 \boxed{d_0\mid 2|A|.} \tag{12}
\]
To prove this, d divides both coordinates of \(\Gamma\), and d0 divides T.
Since \(T-\Re\Gamma=2A^2\) and \(\Im\Gamma=2AB\), d0 divides
both \(2A^2\) and \(2AB\). Also \(d_0^2\mid K\) and
\(d_0^2\mid T^2\), so \(d_0\mid2AD\). The gcd of these three
integers is \(2|A|\gcd(A,B,D)=2|A|\). This argument includes the
prime two. It also proves that the defect root in (11) is an integer.

The anchor in (11) has imaginary coordinate \(2hAB/d_0\). Consequently
the maximum absolute imaginary coordinate divided by the fourth root of
the squared radius is at least
\[
 \frac{2|AB|}{K^{1/4}}\sqrt{h/d_0}
 \ \ge\ \sqrt{2h|A|}\,\frac{|B|}{K^{1/4}}. \tag{13}
\]
In the critical regime with \(|A|,|D|=o(|B|)\), the last ratio tends to
one. Thus even the smallest aligned lift with the displayed integer center
requires h|A| to be polylogarithmic if its normalized strip is to be in
Temur's stated range. The current power bounds do not guarantee this.

This is an optimization only within the aligned family in (7)--(11).
A further rational rotation or a different Gaussian multiplier of the
same norm could change the imaginary coordinates; no optimality over
those additional operations is claimed.

## Consequence for the growth-rate question

Equations (1)--(6) are useful exact arithmetic constraints, but neither
content division nor a lattice-index spacing argument improves the strip
estimate established in the lift note. The two separate requirements for
Temur's theorem remain: a sufficiently small near-square defect and a strip
of polylogarithmic width after fourth-root normalization. Polynomial-small
common content alone does not establish the second requirement.

No new contradiction for growing cluster size follows here. Any stronger
argument using (5) would need an additional estimate on its small
coefficients, on the allowed Gaussian orientations, or on their simultaneous
occurrence across the actual source rows. The existence of the admissible
map producing this critical form also remains an independent hypothesis.

The [exact lift checker](check_critical_shear_near_square_lift.py) verifies
the shared-content and complementary-divisor identities on 1,720 cases.
The independent [canonical-lift checker](check_mobius_canonical_discriminant_lift.py)
passes 2,205 actual-anchor/output-factor covariance cases and 3,840 primitive
triangular cases, including 1,680 nontrivial integer-center content factors.
It checks the rational centers, defect roots, radial/stretch identity, and
the divisibility (12) exactly.
