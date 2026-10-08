# A primitive-chord model and an effective growing-height radial slice

## Status

This note does **not** prove the uniform endpoint bound or the missing
inequality (D11) in Appendix D of multipoint_continuation.md.

It gives two independent pieces of progress:

1. An exact divisor/cofactor description of primitive chords from a lattice
   point, including the parity conditions.
2. A uniform bound of three points on a one-sided endpoint arc, or five
   points on an arc containing the reference point, when that point's
   primitive radial height is bounded by any fixed power
   \((\log R)^\sigma\) with \(\sigma<1/2\). This uses a verified effective simultaneous-Pell
   theorem, with its coefficient dependence retained explicitly.

The second statement strengthens the bounded-height rational-normal slice
in multipoint_continuation.md: the permitted radial height now tends to
infinity. Squarefree common-conductor obstructions generally have radial
height comparable to \(R\), so this still misses the hardest configurations.

## 1. Exact primitive chord and cofactor

Let \(P=x+iy\in\mathbf Z[i]\), \(N=x^2+y^2>0\), and let
\(Q\ne P\) be another Gaussian integer of norm \(N\). Write uniquely

\[
Q-P=dv,\qquad d=\gcd(|\Re(Q-P)|,|\Im(Q-P)|)>0,\qquad
v=a+ib,\qquad \gcd(a,b)=1.
\]

Put
\[
q=a^2+b^2,\qquad A=xa+yb,\qquad B=ya-xb.
\]

Equality of the endpoint norms is exactly
\[
2A=-dq. \tag{1}
\]

Because \(\gcd(a,q)=\gcd(a,b^2)=1\), the identity
\[
a(2B)=2yq-b(2A)
\]
shows that \(q\mid2B\). Consequently
\[
\boxed{\quad
\frac{2P}{v}=-d+it\in\mathbf Z[i],\qquad
t=\frac{2B}{q}\in\mathbf Z,\qquad
2P=v(-d+it),\quad 2Q=v(d+it).
\quad} \tag{2}
\]

In particular \(v\mid2P\). Taking norms,
\[
q(d^2+t^2)=4N. \tag{3}
\]

A primitive sum of two squares has \(v_2(q)\leq1\), since two odd
coordinates give \(q\equiv2\pmod4\). Combining this with (3) gives
\[
q\mid2N. \tag{4}
\]

There is useful parity information:

- If \(q\) is odd, then \(a,b\) have opposite parity. Reducing (1)
  modulo two makes \(d\) even; reducing \(v(-d+it)=2P\) modulo two
  then makes \(t\) even. In this case \(v\mid P\).
- If \(q\) is even, then \(a,b\) are odd and the two coordinates of
  \(v(-d+it)\) show \(d\equiv t\pmod2\).

Conversely, a primitive \(v\in\mathbf Z[i]\) with \(v\mid2P\) and
\(\Re(2P/v)=-d<0\) for an integer \(d\) produces
\(Q=P+dv\in\mathbf Z[i]\), and (1) gives \(|Q|=|P|\). Thus this is
an exact parametrization, not merely a necessary divisor condition.
For a fixed oriented primitive \(v\), there is at most one such \(Q\).

If the containing arc has length at most \(C\sqrt R\), where
\(R=\sqrt N\), then
\[
d\sqrt q=|Q-P|\leq C\sqrt R,
\qquad
\frac{d}{\sqrt{d^2+t^2}}\leq\frac{C}{2\sqrt R}. \tag{5}
\]

Thus the cofactor in (2) is a Gaussian divisor of \(2P\) lying close
to the imaginary axis. The original primitive chord has norm
\(q\leq C^2R\) and divides \(2N\) at the level of rational norms.
These facts alone do not give a bounded divisor count.

## 2. Integral rotation to a perfect-square circle

For a nonzero lattice point \(P=x+iy\) of modulus \(R\), set
\[
g=\gcd(|x|,|y|),\qquad
Q_{\rm rad}=\frac{R}{g},\qquad
\alpha=\frac{\overline P}{g}\in\mathbf Z[i].
\]

Here \(Q_{\rm rad}\) is the length of the primitive integer radial
vector. Multiplication by \(\alpha\) injects every lattice point on the
original circle into the integer lattice on the circle of radius
\[
S=Q_{\rm rad}R=\frac{R^2}{g}\in\mathbf Z_{>0}. \tag{6}
\]

The integrality of \(S\) follows from \(g^2\mid R^2=x^2+y^2\).
The selected point \(P\) maps to \(\alpha P=S\) on the positive real
axis. All angles and their order are preserved.

Suppose \(P\) is at an endpoint of an arc of length \(C\sqrt R\),
and \(C/\sqrt R<\pi/2\). Reflect if necessary so the other points
have positive imaginary part after rotation. Write them as
\[
X_i+iY_i,\qquad Y_i>0,\qquad h_i=S-X_i\in\mathbf Z_{>0}.
\]

Along this one-sided arc the \(h_i\)'s are strictly increasing, and
\[
Y_i^2=h_i(2S-h_i),\qquad
h_i=S(1-\cos\theta_i)
 \leq\frac{C^2Q_{\rm rad}}2. \tag{7}
\]

The angle bound, rather than an estimate using only the transformed
radius, is essential. Integral rotation preserves the endpoint angle,
but multiplies the effective arc constant by \(\sqrt{Q_{\rm rad}}\).

## 3. A quantitative three-offset obstruction

We use Bugeaud, *Effective simultaneous rational approximation to pairs
of real quadratic numbers*, Theorem 1.2
([primary preprint, page 3](https://arxiv.org/pdf/1907.10253)).
It supplies an absolute effectively computable \(K>0\) such that
\[
x^2-ay^2=u,\qquad z^2-by^2=v
\quad\Longrightarrow\quad
\max(x,y,z)\leq
\max(|u|,|v|,2)^{K\sqrt{ab}\log a\log b}, \tag{8}
\]
for positive integer \(x,y,z,a,b\), provided \(u,v\ne0\) and
none of \(a,b,ab\) is a square. The dependence on \(a,b\) is part of
the theorem and is not suppressed.

**Lemma.** Let \(S\geq2\) be an integer and \(H\geq2\) be real.
Suppose
\[
H\leq S,\qquad H^3<2\sqrt S,\qquad
\log S>32K H^2(\log H)^3. \tag{9}
\]
Then there do not exist three distinct positive integers
\(h_1<h_2<h_3\leq H\) and positive integers \(Y_1,Y_2,Y_3\) with
\[
Y_i^2=h_i(2S-h_i). \tag{10}
\]

*Proof.* For \(i<j\), eliminating \(S\) gives
\[
h_iY_j^2-h_jY_i^2=h_ih_j(h_i-h_j)\ne0. \tag{11}
\]

First, the three \(h_i\)'s have distinct squarefree parts. Otherwise,
write \(h_i=st_i^2\), \(h_j=st_j^2\), with \(s\) squarefree. Since
\(h_i,h_j\leq S\), both \(t_iY_j\) and \(t_jY_i\) are at least
\(\sqrt S\). They are unequal integers by (11), so
\[
\begin{aligned}
|h_iY_j^2-h_jY_i^2|
 &=s\,|(t_iY_j)^2-(t_jY_i)^2|\\
 &\geq 2\sqrt S>H^3.
\end{aligned}
\]
This contradicts the upper bound \(H^3\) for the right side of (11).

Now take
\[
a=h_1h_3,\quad b=h_2h_3,\quad
x=h_3Y_1,\quad y=Y_3,\quad z=h_3Y_2.
\]
Then
\[
x^2-ay^2=h_1h_3^2(h_3-h_1),\qquad
z^2-by^2=h_2h_3^2(h_3-h_2). \tag{12}
\]
The residues are nonzero. Distinctness of the squarefree parts makes
each of \(a,b,ab\) nonsquare. Moreover
\[
a,b\leq H^2,\quad
\sqrt{ab}\leq H^2,\quad
U:=\max(|u|,|v|,2)\leq H^4.
\]
Applying (8),
\[
\log\max(x,y,z)
\leq16K H^2(\log H)^3. \tag{13}
\]
But \(Y_3^2\geq S h_3\geq S\), so the left side is at least
\(\tfrac12\log S\). This contradicts (9). \(\square\)

This argument uses the common denominator in the two Pell equations
and controls their coefficient growth. It does not apply a theorem
uniformly over unspecified quadratic fields.

## 4. The growing-height uniform slice

Fix \(C,A>0\) and \(0<\sigma<1/2\). There is an effectively computable
\(R_0(C,A,\sigma)\) such that the following holds for \(R>R_0\).

> If an arc of length at most \(C\sqrt R\) begins at a lattice point
> \(P\) whose primitive radial height satisfies
> \[
> Q_{\rm rad}\leq A(\log R)^\sigma,
> \tag{14}
> \]
> the arc contains at most **three** lattice points.
>
> If \(P\) can lie anywhere on the arc, the bound is **five**.

Indeed, apply the lemma with
\[
H=\max\left(2,\frac{C^2Q_{\rm rad}}2\right),
\qquad S=Q_{\rm rad}R.
\]
Since \(1\leq Q_{\rm rad}\leq R\),
\[
\log R\leq\log S\leq2\log R.
\]
Condition (14) gives
\[
H=O_{C,A}((\log R)^\sigma),\qquad
H^2(\log H)^3
 =O_{C,A,\sigma}((\log R)^{2\sigma}(\log\log R)^3)
 =o(\log R).
\]
It also gives \(H^3=o(\sqrt S)\) and \(H\leq S\) for large \(R\).
Thus (9) eventually holds. Three points other than \(P\) would give
the forbidden three offsets. If \(P\) is interior, apply the one-sided
bound to each side and count \(P\) once, giving \(3+3-1=5\).

The finite inequalities (9), together with (6)--(7), are the explicit
version of this statement; the asymptotic corollary only makes the
permitted growth easy to read.

There is also a directly applicable weaker published slice:
Temur, *Discrete fractional integrals, lattice points on short arcs,
and diophantine approximation*, Theorem 6
([primary preprint](https://arxiv.org/html/2012.10784)),
bounds by 20 the points in an axis window of width
\(6S^{1/2}(\log S^2)^{\kappa/4}\), for \(\kappa<1/2\) in its stated
range and sufficiently large integer-radius \(S\). The perfect-square
case of its radius hypothesis holds automatically in (6). This would
allow \(Q_{\rm rad}\leq36C^{-2}(\log S^2)^{\kappa/2}\).
The direct proof above reaches larger radial heights by exploiting the
zero residual in the perfect-square circle.

## 5. Why this does not supply (D11)

The target (D11) asks for fixed \(\delta>0\), \(B\geq0\) such that
every eight-point endpoint subcluster at constant \(1/2\), with its
inherited conductor layers, satisfies
\[
D_T+\delta W_{T,4}\leq2W+B.
\]

Sections 2--4 exclude such subclusters only when at least one selected
point has the low primitive radial height (14). For a squarefree
circle norm, every lattice position is primitive over the rational
integers, so \(g=1\) and \(Q_{\rm rad}=R\). This far exceeds the
proved range. The exactly balanced squarefree conductor models can
therefore survive this necessary condition.

Applying (8) after integral rotation without bounding \(Q_{\rm rad}\)
only gives (13) with \(H\) as large as a constant times \(R\). Its
right side can then have size \(R^2(\log R)^3\); this yields no
contradiction against a lower bound of order \(\log R\), and it
supplies no positive term proportional to \(W_{T,4}\).

Similarly, the primitive-chord identities (2)--(5) retain exact
Gaussian factors but give no support-independent count of their
allowed divisors. No uniform \(S\)-unit result or bounded residual
clique number has been inferred from these reformulations.

The new rigorous restriction is therefore: any hypothetical unbounded
endpoint cluster must eventually avoid every primitive radial height
\(A(\log R)^\sigma\), for every fixed \(A>0\) and \(\sigma<1/2\).
This is compatible with the known arithmetic frontier and does not
resolve the uniform endpoint theorem.
