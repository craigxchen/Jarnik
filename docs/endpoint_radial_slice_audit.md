# Independent audit of the effective radial slice

This audits Sections 1--4 of endpoint_independent_attack.md. The
primitive-chord identities, the finite simultaneous-Pell threshold, and
the resulting bounds of three and five points are correct. They remain
a restricted-height result, not a proof of the unrestricted endpoint
bound.

## 1. Algebra, hypotheses, and constants

The source used in the note is indeed
[Bugeaud, Theorem 1.2, page 3](https://arxiv.org/pdf/1907.10253).
Its constant is absolute and effective, its variables are positive
integers, its two residues are nonzero, and its nonsquare conditions
apply to all three of \(a,b,ab\). The coefficient dependence quoted in
the note is retained.

The primitive chord calculation also covers zero coordinates:
\(\gcd(a,a^2+b^2)=1\) remains valid when \(a=0\), because then
\(b=\pm1\). The Gaussian cofactor has imaginary part
\(2(ya-xb)/(a^2+b^2)\), with the displayed sign. The deduction
\(q\mid2N\) uses both \(q\mid4N\) and \(v_2(q)\leq1\).
The parity conclusions and the converse construction follow directly
from \(2P=v(-d+it)\).

For the radial slice, the transformed radius
\(S=R^2/g\) is integral and satisfies \(R\leq S\leq R^2\).
The offset bound uses the original angle:
\[
h=S(1-\cos\theta)\leq C^2Q_{\rm rad}/2.
\]
Replacing the transformed arc constant by the original one would be
incorrect; the note does not make that replacement.

For two offsets with the same squarefree part,
\(h_i=st_i^2,\ h_j=st_j^2\), one has
\[
(t_iY_j)^2\geq S s t_i^2t_j^2\geq S,
\qquad
(t_jY_i)^2\geq S.
\]
They are unequal integers, because eliminating \(S\) gives the nonzero
quantity \(h_ih_j(h_i-h_j)\). Hence the absolute difference of their
squares is at least their sum, and therefore at least \(2\sqrt S\).
Multiplication by \(s\geq1\) preserves that lower bound. This verifies
the use of \(H^3<2\sqrt S\) to exclude repeated squarefree parts.

For three distinct squarefree parts the proposed substitutions give
\[
\begin{aligned}
(h_3Y_1)^2-(h_1h_3)Y_3^2
  &=h_1h_3^2(h_3-h_1),\\
(h_3Y_2)^2-(h_2h_3)Y_3^2
  &=h_2h_3^2(h_3-h_2).
\end{aligned}
\]
Both residues are positive. The two coefficients are nonsquares,
and their product is nonsquare because its squarefree part is that of
\(h_1h_2\). All theorem hypotheses are satisfied.

The bounds
\[
\sqrt{ab}\leq H^2,\qquad
(\log a)(\log b)\leq4(\log H)^2,\qquad
\log\max(|u|,|v|,2)\leq4\log H
\]
give exactly \(16K H^2(\log H)^3\). Since \(Y_3\geq\sqrt S\), the
strict condition
\[
\log S>32K H^2(\log H)^3
\]
is sufficient for contradiction. No coefficient-independent
Diophantine approximation constant has been inserted.

If the reference point is an endpoint, at most two positive offsets
remain, giving three points including the reference point. If it is
interior, each side has at most two other points; the total is five.
The angle assumption \(C/\sqrt R<\pi/2\) ensures positive imaginary
coordinates after reflection and distinct increasing offsets on each
side. It holds throughout the asserted sufficiently large radius range.

## 2. A slightly stronger growth range from the same finite threshold

The stated powers \((\log R)^\sigma\), \(\sigma<1/2\), can be enlarged
without changing the proof. Let \(K>0\) be the constant in the cited
theorem. Fix \(C>0\) and
\[
0<A<\frac{1}{C^2\sqrt K}.
\tag{1}
\]
For all sufficiently large \(R\), the same three-point and five-point
bounds hold whenever
\[
\boxed{\quad
Q_{\rm rad}\leq
A\frac{\sqrt{\log R}}{(\log\log R)^{3/2}}.
\quad}
\tag{2}
\]

To check the constant, write \(L=\log R\), \(c=C^2/2\), and
\[
H=\max(2,cQ_{\rm rad}),\qquad S=Q_{\rm rad}R.
\]
The upper bound in (2) gives
\[
\limsup_{R\to\infty}
\frac{32K H^2(\log H)^3}{\log R}
\leq4Kc^2A^2
=KC^4A^2<1.
\]
Here the factor \(1/8\) comes from
\((\log H)^3/(\log L)^3\leq1/8+o(1)\).
Since \(\log S\geq\log R\), the strict logarithmic threshold follows.
The other conditions \(H\leq S\) and \(H^3<2\sqrt S\), and the
angle restriction, follow immediately for large \(R\). All thresholds
can be chosen effectively because the input constant is effective and
the comparisons involve elementary monotone bounds.

This extension still leaves \(Q_{\rm rad}=R\), including squarefree
circle norms, completely outside the proved range.
