# Positive even/odd Bezoutians can miss balanced conductor depth

The even/odd decomposition of a positive-root polynomial supplies a
positive definite integral Bezoutian. Its determinant, however, records
opposite roots, rather than repeated roots at one isotropic orientation.
An arbitrarily deep, perfectly balanced conductor cut can be invisible to
the entire local matrix. The examples below are actual positive integer
cotangent cliques, with all pair quotients integral.

This does **not** disprove an endpoint-restricted refinement. Their least
squared radius is at least a constant times the fifth power of the
normalized minimum, and grows still faster as their size increases.
No fourth-power endpoint inequality is proved or refuted here. Compare
the [height target](integer_cotangent_lcm_height_target.md) and the
[derivative concentration audit](integer_cotangent_derivative_concentration.md).

## 1. The positive matrix and its exact determinant

Let \(k\ge2\), \(X_i>0\), and
\[
 P(z)=\prod_{i=1}^k(z+X_i)=E(z^2)+zO(z^2),\qquad
 m=\lfloor k/2\rfloor.
\]
Thus \(P(z)=(-1)^kF(-z)\) for the positive-root polynomial
\(F(t)=\prod_i(t-X_i)\) in the height formulation.
Here \(\deg E=m\); when \(k=2m\), \(E\) is monic and
\(\deg O=m-1\), while for \(k=2m+1\), \(O\) is monic and both degrees
are \(m\). Define the symmetric \(m\times m\) coefficient matrix \(B\) by
\[
 \frac{E(s)O(t)-E(t)O(s)}{s-t}
 =(1,s,\ldots,s^{m-1})B(1,t,\ldots,t^{m-1})^{\mathsf T}.
 \tag{1}
\]
It is integral when the \(X_i\) are integers. Repeated \(X_i\) are
permitted in this section. Then
\[
 \boxed{B>0,\qquad \det B=\prod_{i<j}(X_i+X_j).}       \tag{2}
\]

For positivity, write, for \(y>0\),
\[
 P(iy)=\mathcal R(y)e^{i\phi(y)},\qquad
 \phi(y)=\sum_i\arctan(y/X_i),\qquad
 \phi'(y)=\sum_i\frac{X_i}{X_i^2+y^2}>0.
\]
The phase increases from zero to \(k\pi/2\). Consequently \(E\) has
exactly \(m\) distinct negative roots \(r_j=-y_j^2\), corresponding to
the \(m\) finite levels \((j+1/2)\pi\), \(0\le j<m\), of this phase.
At each such root,
\[
 E'(r_j)O(r_j)=
 \frac{\mathcal R(y_j)^2\phi'(y_j)}{2y_j^2}>0.
\]
Evaluation of (1) at pairs of these roots gives
\[
 V B V^{\mathsf T}=\operatorname{diag}(E'(r_j)O(r_j)),
 \qquad V_{ja}=r_j^a.
\]
The Vandermonde matrix \(V\) is invertible, proving positivity.

Here is a normalization check for the determinant in both parities.
The same evaluation identity gives
\[
 \det B=(-1)^{m(m-1)/2}\operatorname{Res}(E,O).          \tag{3}
\]
In general the right side has the additional factor
\(\operatorname{lc}(E)^{m-\deg O}\). It equals one here: \(E\) is
monic in even degree, and the exponent is zero in odd degree.

Direct evaluation at the roots of \(P\) gives
\[
 \operatorname{Res}(P(z),P(-z))
 =2^k P(0)\prod_{i<j}(X_i+X_j)^2.
\]
On the other hand, subtracting \(P(z)\) from \(P(-z)\) inside the
resultant gives
\[
 \operatorname{Res}(P(z),P(-z))
 =2^k P(0)\operatorname{Res}(P(z),O(z^2))
 =2^k P(0)\operatorname{Res}(E,O)^2.                   \tag{4}
\]
For the last equality, evaluate at both square roots of every root of
\(O\). Their two \(P\)-values multiply to \(E^2\). The leading
coefficient factors agree because \(k=2\deg E\) in even degree, and
\(O\) is monic in odd degree. This also covers constant \(O\) at \(k=2\).
Equations (3)--(4) determine the square of the determinant; positivity
chooses the positive sign in (2).

For example, at \(k=5\), writing
\(P=z^5+e_1z^4+e_2z^3+e_3z^2+e_4z+e_5\),
\[
 B=\begin{pmatrix}
 e_3e_4-e_2e_5&e_1e_4-e_5\\
 e_1e_4-e_5&e_1e_2-e_3
 \end{pmatrix}.
\]
Each coefficient \(B_{ab}\), indexed from zero, is homogeneous of
degree \(2k-3-2a-2b\) in the \(X_i\). Thus a common scaling by an integer
\(L\) prime to \(p\) changes \(B\bmod p\) by invertible diagonal factors
and an invertible scalar; local invertibility is unchanged.

## 2. Repeated isotropic roots need not affect this matrix

Fix an odd split prime \(p\nmid L\), and identify \(i\) with one chosen
square root of \(-1\) in \(\mathbb Z_p\). The conductor records the maximum
depths of \(X_i-iL\) and \(X_i+iL\), after the common \(L\)-content is
removed. A block of many \(X_i\equiv iL\pmod p\) can therefore contribute
an arbitrarily large conductor exponent.

By (2), \(B\) is nevertheless invertible over \(\mathbb Z_p\) whenever
\[
 X_i+X_j\not\equiv0\pmod p\qquad(i<j).                \tag{5}
\]
Repeated roots at \(iL\) do not violate (5). An opposite root at
\(-iL\), or another opposite pair of root residues, is required.
In particular, no positive conductor valuation is forced in any Smith
invariant factor of \(B\) by a balanced, one-orientation cluster alone.

For a six-row configuration, including the infinity anchor, take
\(p=17\) and five normalized finite root residues
\[
 (X_i/L)_i\equiv(4,4,4,1,2)\pmod {17}.
\]
Since \(4^2=-1\pmod {17}\), the first three rows can carry a balanced
\(3\mid3\) cut. Yet
\[
 B\equiv\begin{pmatrix}14&10\\10&4\end{pmatrix}\pmod {17},
 \qquad\det B\equiv7\pmod {17}                       \tag{6}
\]
for the normalized roots. Even all five leading Hurwitz determinants
are units in this example: their residues are \(15,4,13,7,12\).
Thus neither the full even/odd Bezoutian nor these leading Hurwitz
pivots acquire this conductor prime.

## 3. Actual positive cliques with a deep balanced cut

The following construction realizes this local phenomenon for every
even row count \(2s\ge6\). There are \(k=2s-1\) finite cotangents.
Choose a fixed split prime
\[
 p>s^2+1,
\]
and put \(P=p^e\), with \(e\ge1\). Choose \(a\) satisfying
\[
 P<a<(p+1)P,\qquad
 v_p((a+jP)^2+1)=e\quad(0\le j<s).                  \tag{7}
\]
For completeness, lift a root \(a_0^2=-1\pmod P\), with
\(0<a_0<P\), and choose \(a=a_0+cP\), \(1\le c\le p\).
After division by \(P\), each expression in (7) excludes one residue
of \(c\bmod p\); its coefficient is \(2a_0\not\equiv0\pmod p\).
There are only \(s<p\) excluded residues, so a choice exists.

Define distinct positive integer slopes
\[
 y=(a,a+P,\ldots,a+(s-1)P,\ P+1,P+2,\ldots,P+s-1).
 \tag{8}
\]
Let \(L\) be the lcm of the reduced denominators of
\((y_i y_j+1)/(y_i-y_j)\), and set \(X_i=Ly_i\).
All \(X_i\) and all required pair cotangents are then integers.
Moreover
\[
 p\nmid L,\qquad A/L=P+1.                            \tag{9}
\]
Indeed the within-cluster differences have valuation exactly \(e\),
and their numerators have valuation at least \(e\). Every other
difference is a \(p\)-unit. To see the latter, neither square root of
\(-1\pmod p\) belongs to \(\{1,\ldots,s\}\): otherwise
\(p\le s^2+1\). The outside residues \(1,\ldots,s-1\) are distinct and
are neither isotropic root.

The three statements below are exact statements about this integer
clique and its least primitive squared radius \(N\):
\[
 \boxed{v_p(N)=e;\quad p\nmid\det B;\quad
 N\ge\frac{P^{2s-1}}{C_s},\qquad
 C_s=2^{s-1}\prod_{1\le b<c\le s-1}(c-b)^2.}          \tag{10}
\]
For the first statement, the \(s\) inside slopes have the same signed
phase valuation of absolute value \(e\), while the anchor and the
\(s-1\) outside slopes have valuation zero. The all-edge valuation
width is therefore exactly \(e\). Every one of its \(e\) layers is a
balanced \(s\mid s\) cut.

For the second statement, let \(r=a\bmod p\). Inside-inside sums
equal \(2r\not\equiv0\). An inside-outside sum could vanish only if
\(p-r\le s-1\), which was excluded above. Outside-outside sums lie
between \(2\) and \(2s-2<p\). Equation (5), and then \(p\nmid L\),
prove the claim for the actual integer roots \(X_i\).

For the radius lower bound, the reduced Gaussian denominator associated
with \(P+b+i\), \(1\le b<s\), has norm at least
\(((P+b)^2+1)/2\ge P^2/2\). Its \(p\)-valuation is zero.
A Gaussian gcd of the denominators for \(b,c\) divides the ordinary
integer \(b-c\), so its norm is at most \((b-c)^2\). Primewise,
\[
 \operatorname{Norm}\operatorname{lcm}_{\mathbb Z[i]}(D_b)
 \ge\frac{\prod_b\operatorname{Norm}(D_b)}
 {\prod_{b<c}\operatorname{Norm}\gcd(D_b,D_c)}.
\]
The additional inside denominator contributes the coprime norm factor
\(P\). The exact Gaussian denominator lcm formula now gives (10).
For \(s=3,p=17\), it reads \(N\ge P^5/4\).

Since \(A/L=P+1\), (10) gives
\[
 \frac{N}{(A/L)^4}\ge\frac{P^{2s-5}}{16C_s}
 \longrightarrow\infty.                              \tag{11}
\]
These are consequently not endpoint counterexamples. The common
denominator also has only the upper bound
\[
 L\le \operatorname{lcm}(1,\ldots,s-1)
             ((p+s)P)^{s(s-1)};                      \tag{12}
\]
the within-group denominators divide bounded index differences, and
there are \(s(s-1)\) cross-group pairs. In particular this construction
does not establish subpower residue content or a full fair profile.

## 4. The remaining arithmetic

The positive Bezoutian detects opposite-root sums. The ordinary
real-root Hankel/derivative determinant instead detects repeated-root
differences, leading back to the existing Vandermonde and derivative
budgets. Positivity alone does not identify the two kinds of content.
The examples prove that even actual pair integrality and a balanced
deep cut do not force the first matrix to record that cut.

An endpoint argument using these matrices would still need a *global*
inequality linking their actual archimedean sizes or their opposite-root
sum divisors to the total conductor. Equations (1)--(12) do not give
that inequality, and do not rule out a refinement using the low-conductor
hypothesis itself.

The [exact checker](check_integer_cotangent_positive_bezoutian.py) checks
both degree parities, positive definiteness by all leading principal
minors, the determinant and scaling identities, the displayed modular
matrix and Hurwitz pivots, and the complete all-edge lcm of the literal
cliques. No spectral approximation or unfactored radius substitute is
used.
