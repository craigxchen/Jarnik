# Quantitative audit of the common-conductor Ru--Vojta route

## Verdict: Outcome B

The desired uniform endpoint arc bound does **not** follow from the cited
Ru--Vojta / Evertse--Ferretti results or from the relevant 2025--2026
refinements checked below.

The fixed-geometric part of the proposed argument is sound and is slightly
stronger than some of the earlier audit documents suggested:

- the level-\(1\) section ratio is \(50/27\);
- the resulting GCD coefficient can be made strictly smaller than \(1/2\);
- the adapted bases actually used in the Ru--Vojta proof form a fixed finite
  family;
- in the common-conductor specialization, that family can be replaced at all
  finite places by two explicit integral bases, chosen only according to the
  orientation of the conductor prime and not according to the point;
- the same paired construction works in every degree, and degree \(2\)
  already has strict endpoint margin: with five sections and
  \(\epsilon=1/4\), its conditional GCD coefficient is \(10/21<1/2\);
- their coefficients, their heights, and the summed Weil-function error are
  independent of the conductor primes and of \(|S|\).

The first failure depends on which substitute is used after the fixed
finite max-over-bases inequality:

1. the intrinsic qualitative theorem applies, but leaves a finite residual
   set for which the cited theorem supplies neither a uniform cardinality nor
   a uniform height bound; the moving family below rules out a uniform height
   bound;
2. the support-independent fixed-twisted theorem still cannot be applied to
   the whole conductor grid, because the Mahler reduction supplies a weight
   tuple \(\mathbf c(P)\) after the point is chosen; an explicit pair of
   points from one conductor, proved below, shows that no single weight tuple
   works even for those two points in the natural fixed two-basis system; and
3. the fixed-product theorem now applies without any partition by local
   bases, but the 2013 exceptional count itself contains \(n|S|\) in its
   exponent.  Already in the minimal quadratic system the exact factor is
   \(900^{5|S|}\).  The 2008 higher-degree proof similarly discretizes
   \(q=(n+1)|S|\) local terms; and
4. Corvaja--Zannier Proposition 2 makes the positive-dimensional exceptional
   directions support-independent in this specialization, but leaves a
   finite residual with no uniform cardinality or height bound.  An explicit
   multiplicatively independent common-conductor endpoint pair shows that
   this residual moves to arbitrarily large height as the conductor primes
   vary.

Ngo--Quang 2025 makes the previously implicit large-height cutoff explicit,
but its exceptional count is still exponential in \((n+1)|S|\); Grieve
2025 strengthens the qualitative description without giving a uniform
count, degree, or threshold.  The direct 2026 Ru--Vojta and GCD applications
identified in the current-literature check remain qualitative.

The common conductor fixes the global form family, the local basis at each
oriented place, and the height scale.  It does not fix one twisted weight
tuple valid for every divisor pair in the conductor grid, and the published
fixed-product bounds do not suppress their explicit support-size dependence.

No PR description was modified.

This audit was performed on branch `agent/uniform-gcd-formalization` at
commit `94b5dbd`.  All files listed in the handoff were read with later audit
documents controlling wherever they corrected earlier optimistic claims.

## 1. Fixed geometry and the exact finite-level inequality

Work over \(K=\mathbf Q(i)\), and put

\[
X=\operatorname{Bl}_{P_0}\mathbf P^2,\qquad
P_0=[1:1:1],\qquad L=6H-E.
\]

Here \(H\) is the pullback of a line and \(E\) is the exceptional divisor.
Let \(D_0,D_1,D_2\) be the pullbacks of the three coordinate lines.  The
blow-up center is on none of those lines, so the \(D_i\) are fixed effective
Cartier divisors in proper position.

There are two height conventions in the primary sources.  Ru--Vojta,
Section 2.1, uses heights relative to \(K\), whereas Evertse--Ferretti and
Corvaja--Zannier use absolute normalized heights.  Since
\([K:\mathbf Q]=2\), the former are twice the latter here.  All global
height and proximity formulas below use the absolute normalization: when
Ru--Vojta Theorem 2.10 is invoked in Section 1.4, its inequality is divided
by two.  This common factor changes no beta coefficient, exceptional set, or
numerical margin.

Concretely, if \(p=\pi\overline\pi\) is a split rational prime, then

\[
\lVert\pi\rVert_\pi=p^{-1/2},
\qquad
-\log\lVert\pi\rVert_\pi=\frac12\log p.
\]

### 1.1 The number \(50/27\)

At level \(N=1\),

\[
h^0(X,L)=\binom{8}{2}-1=27
\]

and, for each \(i\),

\[
\begin{aligned}
\sum_{m\geq 1}h^0(X,L-mD_i)
 &=20+14+9+5+2\\
 &=50.
\end{aligned}
\]

Thus the finite-level ratio is

\[
\beta_0=\frac{50}{27}.
\]

This is the quantity proved in
`GaussianChain/FiniteBeta.lean`.  It should not be identified with the
asymptotic beta constant of Ru--Vojta, Definition 1.9.  Directly taking the
limit is also elementary here.  For every \(N\geq1\),

\[
h^0(X,NL)=\binom{6N+2}{2}-\binom{N+1}{2},
\]

and, since \(NL-mD_i=(6N-m)H-NE\),

\[
\sum_{m\geq1}h^0(X,NL-mD_i)
=\binom{6N+2}{3}-\binom{N+2}{3}
-5N\binom{N+1}{2}.
\]

The leading coefficients are respectively \(35/2\) and \(100/3\), so

\[
\beta(L,D_i)=\frac{40}{21}>\frac{50}{27}.
\]

Accordingly, \(50/27\) is a valid rational number below the asymptotic beta
constant and is the finite-level value used in the direct proof.

In Section 6 of Ru--Vojta, equation (49), the choices \(N=1\),
\(\beta_i=50/27\), and \(\dim X=2\) require

\[
\left(1+\frac{2}{b}\right)
\left(1+\frac{\epsilon_1}{27}\right)<1+\rho.
\tag{1}
\]

For each fixed \(\rho>0\), compatible \(b\) and \(\epsilon_1\) can therefore
be chosen once and for all.

### 1.2 The adapted bases are fixed, not point-dependent forms

The infinite set of all bases adapted to all possible filtrations is not the
family used in the published proof.  In Ru--Vojta, Section 6:

- the finite sets \(\Sigma\) and \(\Delta_\sigma\) are chosen from the fixed
  geometric data;
- one adapted basis \(\mathcal B_{\sigma;\mathbf a}\) is chosen for each
  ordered pair \((\sigma,\mathbf a)\);
- immediately after equation (40), the paper explicitly notes that only
  finitely many such ordered pairs occur;
- equations (49)--(51) take the union of those bases,
  \(\{s_1,\ldots,s_{T_2}\}\), and invoke Theorem 2.10 once using these same
  sections at all \(v\in S\).

For the present three divisors there are three singleton and three pairwise
nonempty members of \(\Sigma\).  Writing
\(a_i=(27/50)m_i\), the defining equation for
\(\Delta_\sigma\) is \(\sum_{i\in\sigma}m_i=b\).  Hence

\[
\#\Delta_{\{i\}}=1,\qquad
\#\Delta_{\{i,j\}}=b+1.
\]

There are exactly \(3+3(b+1)=3b+6\) filtration pairs.  Thus one may bound
the number of selected bases by \(T_1\leq3b+6\), and the number of distinct
selected sections by

\[
T_2\leq 27(3b+6).
\tag{2}
\]

The locally maximizing basis may vary with \((P,v)\), but its index ranges
over this fixed list.  No section coefficient is chosen from a conductor
prime, and no new algebraic form is created from \(P\).

The \(O_v(1)\) in Ru--Vojta equation (50) is likewise an
\(M_K\)-constant arising from comparisons among fixed Weil functions.  It
has finite support, so summing it over a varying \(S\) gives one fixed
\(O(1)\), not an \(O(|S|)\) error.

These facts resolve the finite-basis uncertainty retained in the later
`docs/product_section_pairing_audit.md`,
`docs/three_track_consensus_audit.md`, and the handoff: although all bases
adapted to a given filtration form an infinite collection, the direct
Ru--Vojta proof selects one basis for each member of a fixed finite
filtration list before the rational point and the set \(S\) vary.

### 1.3 A point-independent two-basis system for a common conductor

The preceding finite library can be reduced further in the
common-conductor geometry.  Write

\[
M_{a,b}=X_0^{6-a-b}X_1^aX_2^b,
\qquad R_a=X_1^aX_2^{6-a},
\]

and put

\[
\begin{aligned}
\mathcal B_{12}
 &=\{R_a-R_{a+1}:0\leq a\leq5\}
   \cup\{M_{a,b}-R_a:a+b\leq5\},\\
\mathcal B_0
 &=\{M-X_0^6:M\text{ is a sextic monomial},\ M\ne X_0^6\}.
\end{aligned}
\tag{2a}
\]

Both are integral bases of

\[
V=\{f\in K[X_0,X_1,X_2]_6:f(1,1,1)=0\}.
\]

Indeed, the 28 monomials are vertices and the displayed differences are
the 27 oriented edges of a spanning tree; edge differences along any
spanning tree form a \(\mathbf Z\)-basis of the coefficient-sum-zero
lattice.  In particular the transition determinant between the two bases
is \(\pm1\).

Let \(\nu_i(f)\) be the largest power of \(X_i\) dividing \(f\).  Directly
from (2a), or from the five filtration dimensions \(20,14,9,5,2\),

\[
\sum_{f\in\mathcal B_{12}}\nu_1(f)
=\sum_{f\in\mathcal B_{12}}\nu_2(f)=50.
\tag{2b}
\]

For a monomial \(M=X_0^{e_0}X_1^{e_1}X_2^{e_2}\), symmetry of the 28
sextic monomials gives

\[
\sum_{M\ne X_0^6}e_0=50,
\qquad
\sum_{M\ne X_0^6}e_1
=\sum_{M\ne X_0^6}e_2=56.
\tag{2c}
\]

Choose an orientation of the split part of a conductor \(G\), with
\((G,\overline G)=1\), and use \(\mathcal B_{12}\) at every
\(\pi\mid G\) and \(\mathcal B_0\) at every conjugate place
\(\overline\pi\).  This choice depends on \(G\), but not on the divisors
\(A,B\mid G\) or the point
\(P=[1:A/\overline A:B/\overline B]\).

Fix \(\pi\mid G\), write

\[
a=v_\pi(A/\overline A),\qquad
b=v_\pi(B/\overline B),\qquad m=\max(a,b),
\]

and let \(\ell_\pi=-\log\lVert\pi\rVert_\pi\).  At \(\pi\) the primitive
coordinate valuation vector is \((0,a,b)\).  The 21 interior edges of
\(\mathcal B_{12}\) contribute at least \(35(a+b)\ell_\pi\).  The six
boundary edges contain an additional factor \(X_2-X_1\), and contribute
at least \((15(a+b)+6\min(a,b))\ell_\pi\).  Consequently

\[
\sum_{f\in\mathcal B_{12}}\lambda_{f,\pi}(P)
\geq\bigl(56(a+b)-6m\bigr)\ell_\pi.
\tag{2d}
\]

At \(\overline\pi\) it is \((m,m-a,m-b)\).  Equation (2c) gives

\[
\sum_{f\in\mathcal B_0}\lambda_{f,\overline\pi}(P)
\geq\bigl(162m-56(a+b)\bigr)\ell_\pi.
\tag{2e}
\]

The sum of (2d) and (2e) is therefore at least

\[
156m\ell_\pi,
\tag{2f}
\]

where cancellation in a boundary difference can only increase the left
side.  The sum of the three coordinate-divisor proximities over the
conjugate pair is exactly \(3m\ell_\pi\).  Hence (2f) contains the local
factor \(50\) needed for the finite beta value \(50/27\), plus a uniform
surplus \(6m\ell_\pi\).  Places outside the conductor and the complex
place contribute zero to the coordinate-divisor side for norm-one
coordinates.  Thus the local maximizing-basis index
\(J(P,v)\) is not an obstruction in this application: the two bases in
(2a) give one point-independent placewise form tuple for the whole
conductor grid.

### 1.4 Specialization to the GCD inequality

Let \(x,y\) be \(S\)-units and \(P=[1:x:y]\).  With fixed height
normalizations,

\[
m_S(D_i,P)=h_H(P)+O(1)
\]

for each coordinate divisor, while

\[
h_L(P)=6h_H(P)-h_E(P)+O(1),
\]

where

\[
h_E(P)=\log\operatorname{GCD}^{+}(x-1,y-1)+O(1).
\]

Ru--Vojta Theorem 2.10, inserted into the proof at equations (49)--(51),
therefore gives, outside its exceptional set,

\[
\frac{50}{27}\sum_{i=0}^2m_S(D_i,P)
\leq (1+\rho)h_L(P)+O(1).
\]

After substitution,

\[
h_E(P)\leq
\eta(\rho)h_H(P)+O(1),
\qquad
\eta(\rho)=6-\frac{50}{9(1+\rho)}.
\tag{3}
\]

The limiting coefficient is

\[
\eta(0)=\frac49,
\]

and the exact error and margin are

\[
\eta(\rho)-\frac49
=\frac{50}{9}\frac{\rho}{1+\rho},
\qquad
\frac12-\eta(\rho)
=\frac{1-99\rho}{18(1+\rho)}.
\tag{4}
\]

Consequently,

\[
\eta(\rho)<\frac12
\quad\Longleftrightarrow\quad
0<\rho<\frac1{99}.
\tag{5}
\]

If the uniform additive constant in (3) were \(C_{\mathrm{RV}}\), and the
endpoint comparison lost \(\log C_{\mathrm{arc}}\), the numerical part of
the argument would close once

\[
\log|G|>
\frac{18(1+\rho)}
     {1-99\rho}
\bigl(C_{\mathrm{RV}}+\log C_{\mathrm{arc}}\bigr),
\tag{6}
\]

in addition to any theorem-specific height threshold.  The obstruction
below is precisely that the available exceptional-set results do not supply
the required \(S\)-uniform exceptional complexity and threshold.

### 1.5 The conductor pairing improves the conditional coefficient to \(2/9\)

The surplus in (2f) gives a stronger fixed-product inequality than the
intrinsic beta calculation alone.  Summing over conductor-prime pairs gives

\[
\sum_v\sum_{i=1}^{27}\lambda_{L_i^{(v)},v}(P)
\geq156h_H(P)+O(1).
\tag{6a}
\]

Therefore the usual product-theorem conclusion outside its exceptional
subspaces would imply

\[
156h_H(P)\leq(27+\epsilon)h_L(P)+O(1),
\]

and hence

\[
h_E(P)\leq
\eta_{\mathrm{prod}}(\epsilon)h_H(P)+O(1),
\qquad
\eta_{\mathrm{prod}}(\epsilon)
=6-\frac{156}{27+\epsilon}.
\tag{6b}
\]

In particular

\[
\eta_{\mathrm{prod}}(0)=\frac29,
\qquad
\eta_{\mathrm{prod}}(\epsilon)<\frac12
\Longleftrightarrow
0<\epsilon<\frac{15}{11}.
\tag{6c}
\]

Evertse--Ferretti 2013, Corollary 3.2 assumes
\(0<\epsilon\leq1\), so even its maximal allowed value gives
\(\eta_{\mathrm{prod}}(1)=3/7<1/2\).  Equivalently, an endpoint point has
\(h_L\leq(11/2)h_H+O_C(1)\), while

\[
\frac{156}{11/2}=\frac{312}{11}=27+\frac{15}{11}.
\]

Thus the analytic coefficient is no longer tight on the fixed-product
route.  Its failure below is wholly in uniform exceptional complexity, not
in the available numerical margin.

### 1.6 Degree \(2\) is the smallest-dimensional level with endpoint margin

For every \(d\geq1\), let

\[
V_d=\{f\in K[X_0,X_1,X_2]_d:f(1,1,1)=0\},
\qquad
n_d=\binom{d+2}{2}-1,\qquad
S_d=\binom{d+2}{3}.
\]

The boundary-tree and root-star bases obtained from (2a) in degree \(d\)
satisfy, for \(m=\max(a,b)\),

\[
\begin{aligned}
W_{12}^{(d)}(a,b)&\geq S_d(a+b)-dm,\\
W_0^{(d)}(a,b)&\geq n_d\,d\,m-S_d(a+b).
\end{aligned}
\]

The middle terms cancel over a conjugate pair, giving

\[
W_{12}^{(d)}+W_0^{(d)}\geq d(n_d-1)m.
\tag{6d}
\]

Since \(h_{dH-E}=dh_H-h_E+O(1)\), the conditional product-theorem
coefficient is

\[
\eta_d(\epsilon)
=d-\frac{d(n_d-1)}{n_d+\epsilon}.
\tag{6e}
\]

Its exact endpoint condition is

\[
\eta_d(\epsilon)<\frac12
\quad\Longleftrightarrow\quad
\epsilon<\frac{d(d-1)}{4d-2}.
\tag{6f}
\]

Degree \(1\) has no positive margin.  Degree \(2\) is therefore minimal:

\[
d=2,\qquad n_2=5,\qquad
W_{12}^{(2)}+W_0^{(2)}\geq8m.
\]

Taking \(\epsilon=1/4<1/3\) gives

\[
\eta_2(1/4)=2-\frac8{5+1/4}=\frac{10}{21}<\frac12.
\tag{6g}
\]

Thus, conditionally on a support-independent exceptional theorem at this
level,

\[
h_E(P)\leq \frac{10}{21}h_H(P)+O(1),
\qquad
\frac12-\frac{10}{21}=\frac1{42}.
\tag{6g'}
\]

The large-height interface would also be uniform.  The elementary bound
\(h_E\leq h_H+O(1)\), together with the endpoint lower bound
\(h_E\geq(1/2)\log|G|-O_C(1)\), gives

\[
h_{2H-E}=2h_H-h_E+O(1)
\geq \frac12\log|G|-O_C(1).
\tag{6g''}
\]

Consequently every fixed theorem threshold is exceeded once \(|G|\) is
large in terms of the arc constant.  What is missing is the uniform
exceptional family, not this conditional threshold passage.

The two quadratic bases contain at most ten distinct integral forms, their
transition determinant is \(\pm1\), and hence \(D=1\) and
\(\Delta_{\mathcal L}=1\).  Substitution in Evertse--Ferretti 2013,
Corollary 3.2 gives the exact exceptional count

\[
\boxed{
900^{5|S|}
10^{10}2^{16}5^{15}
\log12\,\log(4\log3).
}
\tag{6h}
\]

Thus lowering the section degree substantially improves every fixed
parameter but does not remove the exponential support dependence.

## 2. What the qualitative theorem actually supplies

Ru--Vojta Theorem 2.10 is a fixed-section version of Schmidt's Subspace
Theorem.  Its inputs include \(k,S,X,L,V,s_1,\ldots,s_q,\epsilon,c\), and it
produces a proper closed set \(Z\).  Thus the theorem itself permits \(Z\)
to depend on \(S\).

Ru--Vojta Remark 2.12 gives the strongest relevant refinement.  When the
linear-system map is generically finite, it writes the exceptional locus as

\[
Z_1\cup F^{\mathrm{RV}}_{k,S,\epsilon,c},
\tag{7}
\]

where \(Z_1\) depends only on \(L,V,s_1,\ldots,s_q\), but
\(F^{\mathrm{RV}}_{k,S,\epsilon,c}\) is a finite union of points.  The system
\(|6H-E|\) is base-point-free and its map is generically finite, so this
remark applies.

Rousseau--Turchet--Wang, Theorems 3.4 and 3.5 in the linked arXiv version,
make the same distinction explicit: the proper closed set can
be made independent of the number field and \(S\), but the inequality is
asserted only for “all but finitely many” remaining rational points.

Neither result bounds

\[
\#F^{\mathrm{RV}}_{k,S,\epsilon,c}
\quad\text{or}\quad
\max_{P\in F^{\mathrm{RV}}_{k,S,\epsilon,c}}h(P)
\tag{8}
\]

independently of \(S\).  This is already fatal for the proposed application.
For each \(G\), its divisor-pair grid is finite, so “all but finitely many”
allows that entire grid to be contained in (8).  Equivalently, the scale
\(\log|G|\) may lie below an uncontrolled \(S_G\)-dependent cutoff.

This is the first failure if one follows the intrinsic qualitative
Ru--Vojta proof.

## 3. Audit of the quantitative substitutes

There are three different dimensions denoted by \(n\) in the source papers.
In Ru--Vojta equation (49), \(n=\dim X=2\).  In
Evertse--Ferretti 2013, Theorem 2.1 and Corollary 3.2, \(n\) is the vector
dimension, here
\(\dim H^0(X,L)=27\) (projective dimension \(26\)).  In
Evertse--Ferretti 2008, Theorem 1.3, \(n\) is again the dimension of its
input projective variety.  The bounds below use each paper's own convention.

The quantitative papers offer strong uniformity for **one fixed ordered
local form system and one fixed local weight system**.  A maximizing
coordinate that varies with the point inside that one tuple is harmless:
the maximum is built into the definition of a twisted height.  Section 1.3
supplies the fixed ordered form system in the present application.  What is
not supplied is a point-independent normalized weight tuple valid for the
whole conductor grid.

### 3.1 Evertse--Ferretti 2013: fixed twisted system

Evertse--Ferretti 2013, equations (2.4)--(2.10) and Theorem 2.1, fix the
entire tuple

\[
\mathcal L=(L_i^{(v)}),\qquad
\mathbf c=(c_{iv})
\]

before the rational point varies.  The weights satisfy

\[
\sum_i c_{iv}=0,\qquad
\sum_v\max_i c_{iv}\leq1,
\]

and have finite support.

The form parameter is

\[
r=\#\bigcup_v\{L_1^{(v)},\ldots,L_n^{(v)}\},
\qquad R\geq r.
\tag{9}
\]

This is the number of **literally distinct forms**, not the number of
indexed occurrences and not, as stated, the number modulo proportionality.
Identical repetitions at arbitrarily many places do not enlarge \(r\);
proportional but unequally normalized forms do enlarge it unless a separate
normalization argument is supplied.

The determinant parameters are

\[
\Delta_{\mathcal L}
=\prod_v\|\det(L_1^{(v)},\ldots,L_n^{(v)})\|_v
\tag{10}
\]

and

\[
H_{\mathcal L}
=\prod_v\max_{1\leq i_1<\cdots<i_n\leq r}
\|\det(L_{i_1},\ldots,L_{i_n})\|_v.
\tag{11}
\]

Theorem 2.1 gives

\[
t_0\leq
10^6\,2^{2n}n^{10}\delta^{-3}
\log(3\delta^{-1}R)
\log(\delta^{-1}\log 3R)
\tag{12}
\]

and applies for

\[
Q\geq C_0
=\max\left(H_{\mathcal L}^{1/R},n^{1/\delta}\right).
\tag{13}
\]

For the two integral bases in (2a), \(R\leq54\), all coefficient heights
are fixed, and every transition determinant is \(\pm1\).  Taking
\(\mathcal B_0\) as the reference system outside the conductor gives

\[
\Delta_{\mathcal L}=1
\tag{13a}
\]

for every \(G\).  Thus \(t_0\) and \(C_0\) are independent of \(|S|\), the
support of \(\mathbf c\), and the conductor-prime identities.  In
particular neither coefficient heights nor determinant products hide a raw
\(|S|\)-loss.

The exceptional subspaces nevertheless depend on the full ordered
placewise form tuple and weight tuple
\((\mathcal L,\mathbf c)\).  Both are quantified before the points.

The product-to-component reduction, however, has the opposite quantifier
order.  Evertse--Ferretti 2013, equations (1.1)--(1.2), and the explicit
Mahler reduction in Evertse--Ferretti 2008, Lemma 5.2 and equation (5.12),
give

\[
\forall P\text{ in the conductor grid}\;
\exists\mathbf c(P)
\tag{14}
\]

making that point small for a twisted height.  Theorem 2.1 instead needs

\[
\exists\mathbf c_G\;
\forall P\text{ in the conductor grid}
\tag{15}
\]

with the already fixed form tuple \(\mathcal L_G\).  This is not just a gap
in the present construction.  The following explicit family disproves (15)
for the natural two-basis tuple.

Let \(T=30q+10\), with \(q\to\infty\), and set

\[
H_1=T+i,\qquad H_2=T+3i,\qquad G_T=H_1H_2,
\]

\[
u_1=H_1/\overline H_1,\qquad
u_2=H_2/\overline H_2,
\qquad
P_T=[1:u_1:u_2],\quad P'_T=[1:u_2:u_1].
\tag{15a}
\]

The four Gaussian integers \(H_1,H_2,\overline H_1,\overline H_2\) are
pairwise coprime.  The same-factor gcds divide \(2i\) or \(6i\), the
cross-factor gcds divide \(2i\) or \(4i\), and the congruence on \(T\)
excludes the primes above \(2\) and \(3\).  At primes dividing
\(H_1,\overline H_1,H_2,\overline H_2\), respectively, the passive
coordinate residues are

\[
-\frac12,\quad-2,\quad\frac12,\quad2.
\]

The only rational primes that can make one of their first six powers equal
to \(1\) are \(3,5,7,11,31\).  The four primes other than \(5\) are inert
and cannot divide a primitive \(H_j\) (with \(3\mid H_2\) separately
excluded by \(T\equiv1\pmod3\)); and \(T\equiv0\pmod5\) gives
\(N(H_1)\equiv1\), \(N(H_2)\equiv4\pmod5\).  Therefore the following
valuation profiles are exact, including prime multiplicity.

For \(\mathcal B_{12}\),

\[
\sum\nu_1=\sum\nu_2=50,
\qquad
\sum_{f\in\mathcal B_{12}}|\nu_1(f)-\nu_2(f)|=62.
\tag{15b}
\]

For \(\mathcal B_0\), if \(e_i(M)\) is the exponent of \(X_i\) in \(M\),

\[
\sum(6-e_1)=\sum(6-e_2)=106,
\qquad
\sum_{M\ne X_0^6}|e_1(M)-e_2(M)|=68.
\tag{15c}
\]

Use \(\mathcal B_{12}\) at primes dividing \(G_T\) and \(\mathcal B_0\)
at conjugate primes and everywhere else.  Put

\[
\mathfrak h_T=\log|G_T|
=\frac12\bigl(\log N(H_1)+\log N(H_2)\bigr),
\tag{15d}
\]

where \(|\cdot|\) is the ordinary complex modulus, consistently with the
absolute values in Evertse--Ferretti.  Then

\[
h_H(P_T)=h_H(P'_T)=\mathfrak h_T+O(1).
\]

Moreover,

\[
u_1-1=\frac{2i}{T-i},\qquad
u_2-1=\frac{6i}{T-3i},
\]

so \(\lambda_{E,\infty}=\mathfrak h_T/2+O(1)\) and

\[
h_L(P_T),h_L(P'_T)\leq\frac{11}{2}\mathfrak h_T+O(1).
\tag{15e}
\]

Equations (15b)--(15c) give, for both points,

\[
\sum_v\sum_{i=1}^{27}\lambda_{L_i^{(v)},v}(P)
=156\mathfrak h_T+O(1).
\tag{15f}
\]

Since \(156>28\cdot11/2=154\), both points satisfy the fixed-form product
inequality of Evertse--Ferretti 2013, equation (3.11), with
\(\epsilon=1\), and hence with every smaller positive \(\epsilon\).

It remains to test whether one weight tuple can make both points twisted
small.  For \(a,b\in\mathbf R^d\), put
\(\phi_c(a)=\min_i(a_i+c_i)\).  The exact minimax identity

\[
\sup_{\sum_i c_i=0}\bigl(\phi_c(a)+\phi_c(b)\bigr)
=\overline a+\overline b
-\frac1d\min_{r\in\mathbf R}\sum_i|a_i-b_i-r|
\tag{15g}
\]

follows because prescribed minima \(s,t\) are feasible exactly when
\(\sum_i\max(s-a_i,t-b_i)\leq0\); write \(u=s+t\), \(r=s-t\), and use
\(\max(x,y)=(x+y+|x-y|)/2\).  The median in (15g) is \(0\) for both
profile pairs in (15b)--(15c), so their combined local optimum over a
conjugate prime pair is

\[
\frac{100-62}{27}+\frac{212-68}{27}=\frac{182}{27}.
\tag{15h}
\]

For a finite-support zero-sum weight system define

\[
F_{\mathbf c,Q}(z)
=\sum_v\min_i\bigl(\lambda_{L_i^{(v)},v}(z)+c_{iv}\log Q\bigr).
\]

At all other places a minimum is at most the coordinatewise mean, while
the \(\mathcal B_0\)-coordinate First Main identity is

\[
\sum_v\sum_i\lambda_{\mathcal B_{0,i},v}(z)=27h(z).
\]

Accounting for the change from \(\mathcal B_0\) to
\(\mathcal B_{12}\) at the oriented places, (15h) gives

\[
F_{\mathbf c,Q}(P_T)+F_{\mathbf c,Q}(P'_T)
\leq2h(z_T)-\frac{30}{27}\mathfrak h_T+O(1).
\tag{15i}
\]

Coordinate exchange preserves the \(\mathcal B_0\)-height.  Hence at
least one point \(P_T^*\) satisfies

\[
F_{\mathbf c,Q}(P_T^*)
\leq h(z_T)-\frac{15}{27}\mathfrak h_T+O(1).
\]

By Evertse--Ferretti equation (2.13),

\[
\log H_{\mathcal L,\mathbf c,Q}(z)=h(z)-F_{\mathbf c,Q}(z),
\]

and consequently

\[
\log H_{\mathcal L,\mathbf c,Q}(P_T^*)
\geq\frac59\mathfrak h_T+O(1)\longrightarrow+\infty.
\tag{15j}
\]

The two swapped points have the same \(\mathcal B_0\)-height.  Mahler's
product-to-twisted reduction therefore chooses the same forced scale
\(Q=H(z)^{1+\epsilon/27}\) for both.  The minimax estimate excludes a
common weight at that common scale; in fact it proves the stronger assertion
for every \(Q\geq1\).

This holds for every \(Q\geq1\) and every finite-support weight system
satisfying equation (2.8), even without the normalization (2.9).  Since
\(\Delta_{\mathcal L}=1\), the selected point is not in the small locus
of equation (2.16) for any \(\delta>0\).  Rescaling the individual local
forms does not evade the conclusion: subtracting each place's geometric
mean scale is absorbed into another zero-sum weight, and the remaining
mean cancels against \(\Delta_{\mathcal L}^{1/27}\) in (2.16).

Thus (15), for the direct \(\mathcal B_{12}/\mathcal B_0\) system, is
false.  This does not rule out a genuinely new conductor-dependent
embedding or form system.  It proves that the first point-dependent datum
remaining in the direct route is exactly

\[
P\longmapsto\mathbf c(P),
\tag{15k}
\]

not a local basis index.  Applying Theorem 2.1 separately to these weights
does not help: although (12) is uniform for each application, the union of
the output subspaces as \(P\) varies is not bounded by \(t_0\).

The same obstruction already occurs in the minimal quadratic system of
Section 1.6.  For the two points (15a), the swapped
\(\mathcal B_{12}^{(2)}\)-profiles and reflected
\(\mathcal B_0^{(2)}\)-profiles satisfy

\[
\begin{array}{c|cc}
&\text{two profile sums}&\text{absolute-difference sum}\\ \hline
\mathcal B_{12}^{(2)}&2,\ 2&4\\
\mathcal B_0^{(2)}&6,\ 6&6.
\end{array}
\tag{15l}
\]

Equation (15g) therefore makes the best shared local minimum

\[
\frac{4-4}{5}+\frac{12-6}{5}=\frac65,
\tag{15m}
\]

whereas point-dependent weights attain \(16/5\).  Here

\[
h_{2H-E}(P_T)=h_{2H-E}(P'_T)
=\frac32\mathfrak h_T+O(1),
\]

and each point has product proximity \(8\mathfrak h_T+O(1)\).  Thus both
are solutions of Evertse--Ferretti equation (3.11) for
\(\epsilon=1/4\), since
\(8>(5+1/4)(3/2)\).  Repeating the First Main calculation with five
coordinates gives

\[
F_{\mathbf c,Q}(P_T)+F_{\mathbf c,Q}(P'_T)
\leq2h(z_T)-\frac65\mathfrak h_T+O(1),
\]

so one of them satisfies

\[
\log H_{\mathcal L,\mathbf c,Q}(P_T^*)
\geq\frac35\mathfrak h_T+O(1)\longrightarrow+\infty.
\tag{15n}
\]

Hence passing from 27 sections to the quantitatively optimal five sections
does not produce a common twisted weight.

The canonical-subspace refinement does not repair this.  Evertse--Ferretti
2013, Theorem 2.3, defines \(T(\mathcal L,\mathbf c)\) and gives

\[
m_0=
\left[
10^5\,2^{2n}n^{10}\delta^{-2}
\log(3\delta^{-1}R)
\right],
\qquad
\omega_0=\delta^{-1}\log 3R.
\tag{16}
\]

It also introduces exceptional intervals

\[
[Q_j,Q_j^{\omega_0}),\qquad 1\leq j\leq m_0.
\tag{17}
\]

The cited theorem states no upper bound for the endpoints
\(Q_j(\mathcal L,\mathbf c)\).  Thus it supplies neither a uniform
large-\(Q\) threshold as \(\mathbf c\) varies nor control at the single
conductor scale \(Q_G\).  Proposition 17.5 bounds the height of each
canonical subspace by the fixed form height raised to \(4^n\); as the paper
notes after equation (2.21), the canonical space therefore belongs to a
finite collection depending only on \(\mathcal L\).  That observation does
not control the small loci at the parameter intervals (17), or the spaces
introduced when one derives a covering theorem as
\((\mathcal L,\mathbf c)\) varies.

This preserves the positive conclusion of
`docs/canonical_hn_uniformity.md`: for a fixed finite flag/form library, the
canonical HN space ranges over a fixed finite list independently of \(S\)
and of the numerical weights.  It is the interval endpoints and the
additional covering spaces, not the canonical HN space, that remain
uncontrolled.  The exponent in Proposition 17.5 is \(4^n\), as the primary
TeX and its proof state; the display \(4n\) in
`docs/arrangement_closure_audit.md` is a transcription error.

### 3.2 Evertse--Ferretti 2013: the fixed product system

Evertse--Ferretti 2013, Corollary 3.2 begins after one ordered independent
\(n\)-tuple of forms has already been fixed at each \(v\).  For that one
placewise tuple, with \(s=|S|\), it bounds the number of exceptional
subspaces by

\[
\boxed{
(9n^2\epsilon^{-1})^{ns}
\;10^{10}2^{2n}n^{15}\epsilon^{-3}
\log(3\epsilon^{-1}D)
\log(\epsilon^{-1}\log 3D)
}.
\tag{18}
\]

Here \(D\) bounds the coefficient-field degrees
\([K(L_i^{(v)}):K]\).  It is fixed for the present \(K\)-defined section
library, and in fact \(D=1\).  By (6c), the direct conductor estimate has
enough slack to take any \(0<\epsilon\leq1\), including \(\epsilon=1\).
The fixed Weil-function and form-normalization errors can be absorbed above
a threshold depending on the fixed form library and the gap
\(15/11-\epsilon\), but not on \(S\).  This does not remove the dependence
on \(s\).
The corollary states a lower condition \(H(x)\geq H_0\), but
\(H_0\) is not defined in the corollary or in the immediately preceding
notation; its proof is omitted.  This audit therefore does not infer any
threshold uniformity from that symbol.

For the full level-\(1\) section space, \(n=27\), the exact bound for the
point-independent two-basis tuple (with \(D=1\)) is

\[
(6561\epsilon^{-1})^{27|S|}
10^{10}2^{54}27^{15}\epsilon^{-3}
\log(3\epsilon^{-1})
\log(\epsilon^{-1}\log3).
\tag{19}
\]

Write \(T_{\mathrm{EF13}}(S)\) for the right-hand side of (19).
This is an exact, published \(|S|\)-dependent exceptional-count parameter.
It is present although the set of distinct global forms and the tuple at
each place are fixed.  No additional partition by adapted bases is needed;
the earlier factor \((3b+6)^{|S|}\) was an artefact of a nonoptimal
interface and is not part of the corrected obstruction.  The
fixed-componentwise-system Theorem 3.1 has no corresponding raw \(s\)-factor,
but Section 3.1 proves that its one-weight hypothesis fails for the direct
two-basis system.

At the minimal quadratic level, \(n=5\), \(D=1\), and
\(\epsilon=1/4\).  The exact same corollary gives

\[
\boxed{
T_{\mathrm{EF13}}^{(2)}(S)
=900^{5|S|}
10^{10}2^{16}5^{15}
\log12\,\log(4\log3).
}
\tag{19a}
\]

This is the lowest-dimensional explicit specialization used in the audit.
It still grows
exponentially with \(|S|\).

### 3.3 Evertse--Ferretti 2008: direct higher-degree theorem

The same obstruction appears even more transparently in
Evertse--Ferretti 2008, Theorem 1.3.  With

- \(s=|S|\);
- \(n=\dim X\), \(d=\deg X\);
- \(C=\max[K(f_i^{(v)}):K]\);
- \(\Delta=\operatorname{lcm}\deg f_i^{(v)}\); and
- \(H=\log(2N)+h(X)+\max_{v,i}h(1,f_i^{(v)})\),

and under the additional hypothesis

\[
X(\overline{\mathbf Q})\cap
\{f_0^{(v)}=\cdots=f_n^{(v)}=0\}=\varnothing
\quad(v\in S),
\]

the theorem defines

\[
\begin{aligned}
A_1={}&(20n\delta^{-1})^{(n+1)s}\\
&\cdot
\exp\left(
2^{12n+16}n^{4n}\delta^{-2n}d^{2n+2}
\Delta^{n(2n+2)}
\right)
\log(4C)\log\log(4C),\\
A_2={}&
(8n+6)(n+2)^2d\Delta^{n+1}\delta^{-1},\\
A_3={}&
\exp\left(
2^{6n+20}n^{2n+3}\delta^{-n-1}d^{n+2}
\Delta^{n(n+2)}
\log(2Cs)
\right).
\end{aligned}
\tag{21}
\]

It gives at most \(A_1\) exceptional hypersurfaces, each of degree at most
\(A_2\), only for

\[
h(P)\geq A_3H.
\tag{22}
\]

Thus both the total exceptional complexity \(A_1A_2\) and the threshold
\(A_3H\) depend on \(|S|\), although the individual degree \(A_2\) does
not.

The source of the dependence is not merely the number of distinct
polynomials.  In the proof, Lemma 5.2 is applied to

\[
q=(n+1)s
\tag{23}
\]

separately indexed local logarithmic terms.  Lemma 5.3 then gives at most

\[
(17n\delta^{-1})^{(n+1)s-1}
\tag{24}
\]

weight systems, over whose exceptional families the proof takes a union.
Repeated use of the same \(K\)-defined form can keep the projective
coordinate variety fixed, but it does not reduce \(q\): evaluations at
different places remain different local terms.

For comparison, Evertse--Ferretti 2008, Theorem 2.1, gives genuinely
support-independent bounds for one fixed twisted-height system:

\[
\begin{aligned}
B_1&=
\exp(2^{10n+4}\delta^{-2n}D^{2n+2})
\log(4R)\log\log(4R),\\
B_2&=(4n+3)D\delta^{-1},\\
B_3&=
\exp(2^{5n+4}\delta^{-n-1}D^{n+2}\log(4R)).
\end{aligned}
\tag{25}
\]

It yields at most \(B_1\) hypersurfaces of degree at most \(B_2\) when
\(\log Q\geq B_3(h(Y)+1)\).  But its hypersurfaces depend on the fixed
variety \(Y\) and weight vector \(\mathbf c\), so it again cannot be applied
after the pointwise existential choice in (14).  The formulas above audit
the available generic higher-degree theorem; they are not presented as a
completed construction of its auxiliary \(Y\) for the Ru--Vojta maximum.

There are two natural interfaces, and neither is complete.  Applying
Theorem 1.3 to the embedded surface makes its parameter \(n=2\), so it takes
only three polynomials per place and does not directly encode a product of
all \(27\) members of a basis.  Applying it instead to
\(\mathbf P^{26}\) takes \(n=N=26\) and can use the \(27\) basis
hyperplanes, but an output hypersurface that is proper in
\(\mathbf P^{26}\) may contain the entire embedded surface; its restriction
would then give no exceptional curve.  An additional bridge is required in
either formulation.

For reference, the ambient-\(\mathbf P^{26}\) numerical specialization
\((n,N,d,C,\Delta)=(26,26,1,1,1)\) would already give

\[
\begin{aligned}
A_1&=(520\delta^{-1})^{27|S|}
\exp(2^{328}26^{104}\delta^{-52})\log4\log\log4,\\
A_2&=167776\,\delta^{-1},\\
A_3&=\exp(2^{176}26^{55}\delta^{-27}\log(2|S|)),
\end{aligned}
\]

while Lemma 5.3 contributes
\((442\delta^{-1})^{27|S|-1}\) weight systems.  These are illustrative
dependency checks, not a completed pushdown theorem.

### 3.4 The 2025 refinements still retain the support-size dependence

Ngo--Quang 2025 makes explicit the height bound left undefined in
Evertse--Ferretti 2013, Corollary 3.2, but does not remove the obstruction.
Their Theorem 1.8 fixes \(s=|S|\), one independent
\((n+1)\)-tuple at every place, and the number \(R\) of distinct forms.  Its
equation (1.24) gives

\[
t_2=
\left(4e(n+1)^2\delta^{-1}\right)^{(n+1)s}
10^6\,4^{\,n+7/2}(n+1)^{10}\delta^{-3}
\log^2\!\left(\frac{24R^3D}{\delta}\right).
\]

The now-explicit lower bound in equation (1.25) is

\[
H(x)\geq H_2(n,H^*,\delta)
=
\max\left\{
2(H^*)^{1/[2(n+3)]},
(n+1)^{2(n+1)/\delta}
\right\}.
\]

For fixed form height \(H^*\), this multiplicative-height threshold is
independent of \(s\); thus Theorem 1.8 repairs that particular ambiguity in
the 2013 statement.  But its exceptional count still contains
\((n+1)s\) in the exponent.  In the present
\(\mathbf P^{26}\) section embedding, that exponent is \(27|S|\).

The non-absolute Corollary 1.11 has the same count (with \(D\) absent from
the logarithm) and equations (1.30)--(1.31) give

\[
\begin{aligned}
K_{\mathrm{basic}}
&=s(n+1)h^*
  +\left(\frac{s}{2}+1\right)(n+1)\log(n+1),\\
h_{\mathrm{basic}}
&=\max\left\{
\frac{h^*+\log2}{2(n+3)},
\frac{2(n+1)}{\delta}\log(n+1)
\right\}.
\end{aligned}
\]

Although \(h_{\mathrm{basic}}\) itself is independent of \(s\), absorbing
the additive \(K_{\mathrm{basic}}\) into a fixed coefficient slack requires
a height threshold proportional to \(K_{\mathrm{basic}}\), hence depending
on \(s\).  The more general Corollary 1.14, equation (1.42), states the
dependence directly:

\[
h_{\mathrm{gen}}(n,q,s,h^*,\delta)
=
\frac{(q-n+1)(n+2)(s+1)}{\delta}
\bigl(\log(n+1)+h^*\bigr).
\]

That corollary additionally assumes that the entire local family is in
general position.  No such hypothesis is proved for the union of the
Ru--Vojta adapted bases.  It also assumes
\(\|L_i^{(v)}\|_v\geq1\) in equation (1.37).  This normalization condition
is harmless here: because \(K=\mathbf Q(i)\) and the form library is finite,
normalize each form once to a primitive Gaussian-integer coefficient
vector.  Its coefficient sup-norm is then \(1\) at every nonarchimedean
place and at least \(1\) at the complex place, without using a
conductor-prime-dependent scalar.  Even granting general position would
leave both the exceptional count and the height threshold dependent on
\(s\).

Ngo--Quang's Theorem 1.20 does not provide a hidden support-independent
escape.  Its visible exceptional set is the union of a fixed conjunction
\(\widetilde{\mathcal L}\) of hyperplanes, but its conclusion has an
unspecified \(O(1)\).  There is also an internal discrepancy in the source:
the statement's equation (1.49) prints the coefficient \(q-1+\delta\) for
the sum of the whole \(q+1\)-indexed family, whereas the proof after
equation (6.19) concludes with \(q+\delta\).  The argument here relies on
neither value; neither version is the required dimension-level
max-over-independent-subsets statement.  Moreover, Theorem 1.20 assumes
that the whole form family is non-subdegenerate and is supplied with a
conjunction; those hypotheses have not been proved for the union of the
Ru--Vojta adapted bases.  More decisively, Lemma 6.2, on which the proof
rests, gives in equations (6.4) and (6.7)

\[
\begin{aligned}
t_{\mathrm{mindep}}(t,s,\delta)
&=
\left(4et^2\delta^{-1}\right)^{ts}
10^6\,4^{\,t+5/2}t^{10}\delta^{-3}
\log^2\!\left(\frac{24(t+1)^3}{\delta}\right),\\
h_{\mathrm{mindep}}(t,s,h^*,\delta)
&=
\frac{2(t+1)(s+1)}{\delta}
\bigl(\log t+h^*\bigr).
\end{aligned}
\]

The proof of Theorem 1.20 iterates those finite exceptional families and,
after equation (6.19), absorbs the heights of the remaining finitely many
points into the final \(O(1)\).  Thus the fixed visible conjunction does
not make either that additive constant or the associated effective
large-height threshold uniform in \(s\).

Grieve 2025 gives a later result directly in the Ru--Vojta linear-system
framework, but it is qualitative.  Theorem 1.4 and Corollary 1.5 fix \(S\)
and produce finitely many scattering classes and proper linear sections
without bounding their number uniformly in \(S\).  Theorem 1.6 applies the
Ru--Vojta filtration construction and describes the exceptional set using
finitely many linear sections, again with no bound for their number or
total degree as \(S\) varies.  In the proof, the integer \(m\) is chosen so
that the stable base locus of \(L\) is the base locus of
\(|L_{\overline K}^{\otimes m}|\); hence \(m\) may be fixed from \((X,L)\),
independently of \(S\).  This does not control the number of sections or
scattering classes, their total union complexity, or the finite set left
inside those sections.  Lemma 7.1 is applied to index sets \(S\) and
\(S\times\{0,\ldots,n\}\), while the Northcott step adds finitely many
low-height solutions.  Neither finiteness assertion is quantitative or
uniform in \(S\).  These later results therefore confirm rather than close
the gap above.

### 3.5 Direct 2026 candidates remain qualitative or inapplicable

A current-paper and citation search through 4 September 2026 found no
primary result that changes this conclusion.  In particular:

- Wang--Xiao 2026, Theorem 1.1, restates the Ru--Vojta arithmetic theorem
  with a proper \(Z\) independent of \(k\) and \(S\), but the inequality is
  asserted only for **all but finitely many** points outside \(Z\).
  Theorem 1.6 is especially close to the present interface: its left side
  contains a maximum over a fixed finite collection
  \((I,Y)\in\mathcal M\).  Its conclusion still has the same uncounted
  finite residual set.  The GCD-specific Theorem 1.13 likewise holds for
  all but finitely many points of the integral set outside a proper
  \(Z\).  None of the three theorems bounds the residual cardinality,
  residual heights, or total exceptional degree uniformly in \(S\).
- Yasufuku 2026, Theorem 1.2 (PDF pp. 2--3), is directly about GCD
  inequalities arising from codimension-\(2\) blowups, but its exceptional
  set is explicitly \(Z=Z(k,S,F_1,F_2)\).  Theorem 3.2 and equation (8)
  (PDF p. 6) are the qualitative Ru--Vojta input after \(S\) has been
  fixed.  In the proof of Theorem 1.2 (PDF p. 9), bounded-height points are
  added to the exceptional set without a count or a bound for the cutoff.
  The result therefore supplies no \(S\)-uniform exceptional complexity or
  large-height threshold.
- Adiceam--Shirandami 2026, Theorems 1.3 and 1.6, give probabilistic
  effectivity for varying algebraic/archimedean systems.  Section 1.1
  explicitly says that the “pseudo-deterministic regime,” in which the
  rational-solution height dominates the form height, is not treated.
  Those density estimates are therefore not a deterministic exceptional
  theorem for every conductor system.

This literature search cannot prove that no unpublished or differently
formulated theorem exists.  It does show that the closest identified
post-2025 primary results do not supply the missing quantifier exchange or
uniform constants.

### 3.6 Corvaja--Zannier: uniform subgroup directions, nonuniform finite residue

There is a useful specialization outside the Ru--Vojta product formalism.
Corvaja--Zannier, Proposition 2, states that for fixed \(K,S,\epsilon\),
all but finitely many \(S\)-unit solutions of

\[
-\sum_{\mu\in M_K}\log^-
\max\{|u-1|_\mu,|v-1|_\mu\}
>\epsilon\max\{h(u),h(v)\}
\tag{25a}
\]

lie on subgroups

\[
u^p=v^q,\qquad (p,q)=1,\qquad
\max(|p|,|q|)\leq\epsilon^{-1}.
\tag{25b}
\]

For an endpoint common-conductor pair, the left side of (25a) is at least
\((1/2)\log|G|-O_C(1)\), while each coordinate height is at most
\(\log|G|\).  Taking \(\epsilon=4/9\), every sufficiently large endpoint
pair satisfies (25a), and the integer bound in (25b) becomes
\(\max(|p|,|q|)\leq2\), independently of \(S\).

For non-root off-diagonal pairs, this positive-dimensional part is not
merely of uniformly bounded degree; its incidence with the endpoint grid is
uniformly bounded.  Relative to an arc-end point, write

\[
g_j=\pi_j/\overline\pi_j,\qquad
u=\zeta\prod_jg_j^{d_j},\qquad
v=\xi\prod_jg_j^{d'_j},\qquad |d_j|,|d'_j|\leq e_j.
\]

Valuation at \(\pi_j\) in (25b) gives \(pd_j=qd'_j\).  Since the reference
point is at one end of the arc, all nontrivial ratios have positive lifted
argument.  Opposite-sign relations and the diagonal give no non-root
off-diagonal pair; coordinate-axis relations are precisely the root pairs.
The only remaining non-root relations are

\[
u^2=v\qquad\text{and}\qquad u=v^2.
\tag{25c}
\]

Consider \(u^2=v\), so \(d'_j=2d_j\).  Orient each prime pair for this
ordered pair by setting \(\rho_j=\pi_j\) if \(d_j\geq0\), and
\(\rho_j=\overline\pi_j\) otherwise, and put

\[
W=\prod_j\rho_j^{|d_j|},\qquad
G_d=\prod_j\rho_j^{e_j}.
\]

Then \(u=\zeta W/\overline W\), \(W^2\mid G_d\), and
\(|G_d|=|G|\).  Notice that \(G_d\) may depend on the ordered pair; no
global orthant assumption is being made.  The endpoint estimate

\[
\left|\frac{\zeta W}{\overline W}-1\right|
\leq\kappa_C|G|^{-1/2}
\]

and the fact that a nonzero twisted conjugate difference has modulus at
least \(\sqrt2\) give

\[
|W|\geq\sqrt2\,\kappa_C^{-1}|G|^{1/2},
\qquad
|D|:=|G_d/W^2|\leq\kappa_C^2/2.
\tag{25d}
\]

There are only \(O_C(1)\) possible Gaussian integers \(D\).  Moreover,
\(W^2\mid G_d\) gives \(|W|\leq|G|^{1/2}\), so the numerator in the
preceding endpoint estimate -- one of
\(2a,2b,(1+i)(a+b),(1+i)(a-b)\), up to a unit when
\(W=a+bi\) -- is \(O_C(1)\).  For fixed \(D\), the radius
\(|W|^2=|G|/|D|\) is fixed; a bounded real or imaginary coordinate on that
circle leaves only \(O_C(1)\) Gaussian integers \(W\).  Units add only a
factor four, and \(v=u^2\).  The swapped relation is identical.  Thus all
non-root subgroup directions together contribute only \(O_C(1)\) ordered
pairs; the root/axis directions contribute the separate \(2m\) term in
(25f).

The failure is exactly the proposition's remaining finite set
\(F^{\mathrm{CZ}}_{K,S,\epsilon}\).  It has no stated cardinality or height
bound.  This
is not merely an omission that a uniform height cutoff could repair.  Let

\[
T=30q+10,\qquad
G_T=(T+i)(T+3i),\qquad
u_T=\frac{T+i}{T-i},\qquad
v_T=\frac{T+3i}{T-3i}.
\tag{25e}
\]

The four Gaussian factors are pairwise coprime, so \(u_T,v_T\) are
multiplicatively independent.  Both are
\(O(|G_T|^{-1/2})\)-close to \(1\).  Hence (25a) holds for every fixed
\(\epsilon<1/2\) when \(q\) is large, but no relation (25b) is possible.
The pair (25e) must therefore belong to
\(F^{\mathrm{CZ}}_{K,S_{G_T},\epsilon}\), and
its height tends to infinity.

In the proof of Proposition 2, the first loss appears in the invocation of
Ridout's theorem for “all but finitely many” solutions; the subsequent
partition by \(S^+\) and \(S^-\) has up to \(2^{|S|}\) profiles, and the
Subspace-Theorem/Lang step remains qualitative.  For a cluster of \(m\)
ratios, the theorem yields only

\[
m(m-1)\leq \#F^{\mathrm{CZ}}_{K,S,4/9}+2m+O_C(1).
\tag{25f}
\]

Thus this route reduces the missing input to a uniform bound for the
intersection of the finite residual with one endpoint conductor grid (or
only its off-diagonal clique number), but it does not prove such a bound.

## 4. Dependency table

| Result | \(T\), number | \(\Delta\), individual degree | \(H_0\) / \(Q_0\) | Decisive dependence |
|---|---:|---:|---:|---|
| Ru--Vojta Theorem 2.10 | No quantitative bound | No quantitative bound | No quantitative bound | \(Z\) may depend on \(S\) |
| Ru--Vojta Remark 2.12 / RTW Theorems 3.4--3.5 | Fixed proper closed part; residual finite set uncounted | Fixed part's degree not quantified; residual points have degree \(1\) | No bound stated for residual-point heights | The theorem supplies no uniform count for \(F^{\mathrm{RV}}_{K,S,\rho,c}\), and no uniform height bound |
| EF 2013 Theorem 2.1, one fixed \((\mathcal L,\mathbf c)\) | (12), independent of support size | Linear subspaces | (13), uniform for fixed forms | Output subspaces depend on the full point-independent pair \((\mathcal L,\mathbf c)\) |
| EF 2013 Theorem 2.3 | One canonical subspace; at most (16) exceptional parameter intervals | Canonical subspace linear | No upper bound stated for endpoints (17) | \(Q_j\) depend on \((\mathcal L,\mathbf c)\) |
| EF 2013 Corollary 3.2, one fixed product system | Contains (18) | Linear subspaces | Symbol \(H_0\) is not defined in the statement | Explicit \((9n^2\epsilon^{-1})^{n|S|}\) for the corrected two-basis tuple |
| EF 2008 Theorem 1.3, generic higher-degree bound | \(A_1\), nonuniform | \(A_2\), uniform individually | \(A_3H\), nonuniform | Illustrative only here: \(A_1,A_3\) contain \(s=|S|\), but no proper-surface specialization of the 27-section system is proved |
| EF 2008 Theorem 2.1, one fixed \(\mathbf c\) | \(B_1\), support-independent | \(B_2\), support-independent | \(B_3(h(Y)+1)\) | Hypersurfaces vary with \(\mathbf c\) |
| Ngo--Quang 2025 Theorem 1.8 / Corollary 1.11 | \(t_2\), exponential in \((n+1)|S|\) | Linear subspaces | Theorem 1.8 threshold is explicit and support-independent for fixed \(H^*\); absorbing Corollary 1.11's \(K_{\mathrm{basic}}\) is not | The exact count, and \(K_{\mathrm{basic}}\), retain \(s=|S|\) |
| Ngo--Quang 2025 Corollary 1.14 | \(t_{\mathrm{gen}}\), exponential in \((n+1)|S|\) | Linear subspaces | \(h_{\mathrm{gen}}\), linear in \(|S|+1\) | Also requires the whole local form family to be in general position |
| Ngo--Quang 2025 Theorem 1.20 / Lemma 6.2 | Fixed visible conjunction, but intermediate \(t_{\mathrm{mindep}}\) is exponential in \(t|S|\) | Hyperplanes in the final qualitative statement | Final \(O(1)\) unspecified; \(h_{\mathrm{mindep}}\) is linear in \(|S|+1\) | Non-subdegeneracy/conjunction is unproved here; terminal finite points are absorbed into \(O(1)\), not uniformly bounded |
| Grieve 2025 Theorems 1.4 and 1.6 | Finite, not quantified | The embedding level \(m\) can be fixed from \((X,L)\); the number of classes/sections and total union degree are not quantified uniformly | “Sufficiently large”; low-height points added by Northcott | \(S\) is fixed before the unquantified class count and finite residue are obtained |
| Wang--Xiao 2026 Theorems 1.1, 1.6, and 1.13 | Proper \(Z\), plus an uncounted finite residual set | Not quantitatively bounded | Not quantitatively bounded | The finite max in Theorem 1.6 does not bound the residual set |
| Yasufuku 2026 Theorems 1.2 and 3.2 | Not quantified | Not quantified | Not quantified | \(Z=Z(k,S,F_1,F_2)\); qualitative Ru--Vojta input |
| Adiceam--Shirandami 2026 Theorems 1.3 and 1.6 | Probabilistic density estimates, not a deterministic cover | Not applicable | Wrong height regime for this application | Section 1.1 excludes the pseudo-deterministic regime |
| Corvaja--Zannier Proposition 2 | At most eight subgroup directions when \(\epsilon=4/9\), plus a finite residual with no stated uniform count | Subgroup equations have degree at most \(3\); non-root off-diagonal incidence is \(O_C(1)\), with \(2m\) root/axis pairs | No residual height bound; (25e) shows no uniform one exists | The precise uncontrolled object is \(F^{\mathrm{CZ}}_{K,S,4/9}\) |

For the application, the field degree, variety degree and height, line bundle,
distinct section family, approximation margin, coefficient heights, and
Weil-function constants can all be fixed.  Equations (10)--(11) show that
the determinant products are also harmless once one placewise tuple is
fixed; over the finite basis library, (10) ranges through a fixed finite set.
The remaining nonuniform datum is not a conductor-prime coefficient height.
It is the pointwise Mahler weight \(\mathbf c(P)\), indexed across the
\(|S|\) places, together with the explicit \(|S|\)-dependence of the
fixed-product covering theorem.

## 5. Why a common conductor does not remove the failure

For a fixed Gaussian conductor \(G\), the support \(S_G\) and the scale
\(Q_G\) are common to all ratios

\[
u_A=A/\overline A,\qquad A\mid G.
\]

That gives

\[
h([1:u_A:u_B])\leq\log|G|
\]

and the common endpoint lower bound.  It also means that one Ru--Vojta
application with \(S=S_G\) and the fixed sections covers the whole grid
**qualitatively**.

It does not make the local evaluations common.  At a prime \(v\mid G\), the
valuation of \(u_A\), and hence the componentwise Mahler weight selected from
the local evaluations, changes with the divisor \(A\).  The two bases (2a)
already give a point-independent \(\mathcal L_G\); nevertheless the two
points (15a) prove that no point-independent \(\mathbf c_G\) works for that
tuple.  Fixing the common scale therefore does not interchange the
quantifiers in (14)--(15).

Conjugate-place pairing is what produces both the fixed tuple and the exact
surplus in (2f).  It removes the former local-basis assignment count, but it
does not correlate the numerical component weights selected at different
split primes.  Nor can the repeated-form observation remove (18) or (24):
those factors count locally indexed logarithmic occurrences or their
discretized weight patterns, not distinct algebraic equations.  The
conductor condition imposes no identity known to collapse those separate
local values.

Therefore:

- the fixed-form heights are uniform;
- the local form tuple is point-independent and uses only two integral
  bases;
- normalization does not introduce conductor-prime coefficients in the
  direct Ru--Vojta construction;
- the summed \(M_K\)-error is uniform;
- the paired calculation leaves the stronger coefficient \(2/9\) in the
  zero-margin limit;
- the minimal quadratic specialization gives \(10/21\) at
  \(\epsilon=1/4\), with only five sections;
- but no cited product route supplies both uniformly bounded
  exceptional complexity and a usable uniform threshold.  The exceptional
  count is the decisive nonuniform parameter; some qualitative and
  non-absolute formulations additionally leave the threshold or low-height
  points uncontrolled.

## 6. Conditional pushdown and grid endgame

This section records what would follow if the missing uniform exceptional
theorem were supplied; it is not a claim that it has been supplied.

The complete system \(|6H-E|\) is base-point-free.  A hyperplane section of
its image pulls back to a divisor linearly equivalent to \(L\), whose
pushdown is a plane sextic through \(P_0\).  Hence:

- one exceptional linear subspace contributes a plane curve of degree at
  most \(6\);
- \(T\) exceptional linear subspaces contribute total plane degree at most
  \(6T\);
- more generally, \(T\) auxiliary hypersurfaces of degree at most
  \(\Delta\), **provided each restricts nontrivially to the embedded
  surface**, contribute total plane degree at most \(6T\Delta\).

The proviso in the last item is essential.  A proper ambient hypersurface
may contain the whole embedded surface, as in the unresolved
\(\mathbf P^{26}\) interface in Section 3.3, and then supplies no
exceptional curve on \(X\).  Proper linear subspaces do not have this
problem because the image of a complete linear system is linearly
nondegenerate.

Here \(T\) is understood to include every exceptional component.  If one
instead starts from the qualitative decomposition \(Z_1\cup F\), choose a
fixed plane curve of degree \(D_{Z_1}\) containing the pushdown of the fixed
proper set \(Z_1\); the total curve degree becomes
\(D_{Z_1}+6T\Delta\).  This adds a fixed constant but does not control the
uncounted residual set \(F\).  In that convention the final cluster bound
below is \(D_{Z_1}+6T\Delta+2\), not merely \(6T\Delta+2\).

At the minimal level, \(|2H-E|\) embeds \(X\) as the cubic scroll in
\(\mathbf P^4\).  A proper exceptional linear subspace lies in a hyperplane,
whose pullback pushes down to a conic through \(P_0\).  Consequently \(T\)
linear subspaces at this level have total plane degree at most \(2T\), and
the same formal grid argument would give the sharper conditional cluster
bound

\[
\#\{\text{endpoint-cluster points}\}\leq2T+2.
\tag{26a}
\]

The hypotheses formalized in
`GaussianChain/FiniteExceptionalFamily.lean` would then give

\[
\#\{\text{ratio points}\}\leq 6T\Delta+1
\]

and, after reinstating the distinguished root,

\[
\#\{\text{lattice points in the endpoint cluster}\}
\leq 6T\Delta+2.
\tag{26}
\]

For a single degree-\(D\) curve,
`GaussianChain/UniformExceptionalReduction.lean` gives the corresponding
\(D+2\) bound.  Equations (4)--(6) verify that the analytic error can fit
inside the exact \(1/18\) margin once \(\rho<1/99\) and a genuinely uniform
additive constant and height threshold are available.

The Lean files prove these arithmetic and conditional combinatorial steps.
They do not prove the missing quantitative exceptional-family theorem.

## 7. Earliest precise failed implication and the missing theorem

The earliest failed implication is:

\[
\begin{gathered}
\text{fixed variety and the point-independent two-basis product system}\\
+\ \text{the uniform local gain }8h_H\text{ already in degree }2
\end{gathered}
\quad\not\Longrightarrow_{\text{cited results}}\quad
\begin{gathered}
\text{\(O(1)\) exceptional components}\\
\text{of \(O(1)\) total degree, above an \(O(1)\) threshold,}\\
\text{uniformly in }S.
\end{gathered}
\tag{27}
\]

More precisely, the route-specific earliest failures are:

- **Intrinsic qualitative route:** The theorem applies, but gives no
  uniform bound for the residual finite set in (7)--(8).
- **Fixed-twisted route:** Hypothesis matching fails first.  The available
  reduction produces

  \[
  P\longmapsto\mathbf c(P),
  \]

  and (15a)--(15j) prove that no one \(\mathbf c_G\) can replace even the
  two displayed weights for the direct two-basis system.  This failure
  precedes any appeal to the support-independent numerical bounds
  (12)--(13).
- **Fixed-product route:** Corollary 3.2 applies directly to the two-basis
  tuple, but its lowest-dimensional explicit specialization has the
  explicitly \(|S|\)-dependent exceptional count (19a).
- **Corvaja--Zannier route:** Proposition 2 reduces all non-root
  off-diagonal positive-dimensional subgroup directions to \(O_C(1)\)
  endpoint-grid incidences (with \(2m\) root/axis pairs), but leaves the
  finite residual \(F^{\mathrm{CZ}}_{K,S,4/9}\) uncounted.
  The moving family (25e) shows that its heights cannot be bounded uniformly.

Thus the exact uncontrolled objects or theorem parameters are:

\[
\begin{gathered}
\#F^{\mathrm{RV}}_{K,S,\rho,c},\quad
\max_{P\in F^{\mathrm{RV}}_{K,S,\rho,c}}h(P),\quad
\#F^{\mathrm{CZ}}_{K,S,4/9},\quad
\max_{P\in F^{\mathrm{CZ}}_{K,S,4/9}}h(P),
\quad P\mapsto\mathbf c(P),\\
\text{and, on the fixed-product route,}\quad
T_{\mathrm{EF13}}^{(2)}(S)\text{ in (19a)}
\end{gathered}
\tag{28}
\]

Here (19a) is the lowest-dimensional primary exact quantitative obstruction
recorded in this audit: it is \(|S|\)-dependent for the one fixed product
system already constructed.
Equation (19) records the corresponding sextic count.
The Evertse--Ferretti 2008 weight count (24) and the resulting
\(A_1,A_3\) in (21) corroborate that dependence, but, as Section 3.3
explains, they are not a completed specialization of the 27-section system
to proper exceptional curves on the surface.

The narrowest new input that would close the route is one of the following
alternative sufficient statements, preferably specialized to the fixed
five-dimensional quadratic section space and its two fixed integral bases:

1. **Uniform fixed-product theorem.**  For the quadratic versions of the
   two bases (2a) and \(\epsilon=1/4\), the solutions of the product
   inequality above a threshold
   independent of \(S\) lie in at most \(T\) proper subspaces, with \(T\)
   independent of \(S\).  Linear subspaces already have degree one, and
   (6g) would give the coefficient \(10/21\).
2. **Bounded conductor-weight cover.**  Replace the false single-weight
   assertion (15) by at most \(M\) twisted systems, where \(M\) is independent
   of \(G\) and \(|S_G|\).  This must use more than the direct
   quadratic \(\mathcal B_{12}^{(2)}/\mathcal B_0^{(2)}\) tuple, or otherwise
   evade the explicit two-point obstruction (15l)--(15n).
3. **Uniform finite-residual incidence theorem.**  The intersection of
   \(F^{\mathrm{CZ}}_{K,S,4/9}\) with an endpoint common-conductor grid, or
   just the clique
   number of its off-diagonal-pair graph, is \(O_C(1)\).  Equations
   (25c)--(25d) already handle every subgroup direction.
4. **Conductor-relative one-scale theorem.**  At \(Q=Q_G\), the union of
   small loci over all conductor-realized local weight patterns is contained
   in \(O(1)\) hypersurfaces of uniformly bounded total degree.

The fourth formulation is not an independent shortcut: by
`docs/conductor_relative_equivalence.md`, in the present application it is
essentially the endpoint assertion rewritten in algebraic-geometric terms.

None of these statements is a theorem in the primary sources audited here.
Until one is proved, the common-conductor Ru--Vojta route does not establish
the uniform endpoint arc bound.

## Primary sources

1. M. Ru and P. Vojta, “A birational Nevanlinna constant and its
   consequences,” *American Journal of Mathematics* 142 (2020), 957--991:
   [Theorem 2.10, Remark 2.12, and Section 6 equations
   (40), (49)--(51)](https://doi.org/10.1353/ajm.2020.0022);
   [public author manuscript](https://par.nsf.gov/servlets/purl/10155671).
2. J.-H. Evertse and R. G. Ferretti, “A further improvement of the
   Quantitative Subspace Theorem,” *Annals of Mathematics* 177 (2013),
   513--590: [Theorems 2.1 and 2.3, Corollary 3.2, and Proposition
   17.5](https://doi.org/10.4007/annals.2013.177.2.4);
   [journal page and PDF](https://annals.math.princeton.edu/2013/177-2/p04).
3. J.-H. Evertse and R. G. Ferretti, “A generalization of the Subspace
   Theorem with polynomials of higher degree,” *Diophantine Approximation*,
   Dev. Math. 16 (2008), 175--198:
   [Theorems 1.3 and 2.1, Lemmas 5.2--5.3](https://doi.org/10.1007/978-3-211-74280-8_9);
   [primary preprint](https://arxiv.org/abs/math/0408381).
4. E. Rousseau, A. Turchet, and J. T.-Y. Wang, “Divisibility of polynomials
   and degeneracy of integral points,” *Mathematische Annalen* 388 (2024),
   1969--1999: [arXiv Theorems 3.4--3.5](https://arxiv.org/abs/2106.11337);
   [version of record](https://doi.org/10.1007/s00208-023-02564-3).
5. H. T. Ngo and S. D. Quang, “On absolute and quantitative subspace
   theorems,” *Forum Mathematicum* 37 (2025), 821--849:
   [Theorems 1.8 and 1.20, Corollaries 1.11 and 1.14, and Lemma
   6.2](https://doi.org/10.1515/forum-2023-0247);
   [primary author PDF](https://viasm.edu.vn/Cms_Data/Contents/viasm/Media/2023/tienanpham/Quang-Viasm2023-03.pdf).
6. N. Grieve, “On qualitative aspects of the quantitative subspace
   theorem,” *Rocky Mountain Journal of Mathematics* 55 (2025), 121--146:
   [Theorems 1.4 and 1.6 and Lemma
   7.1](https://doi.org/10.1216/rmj.2025.55.121);
   [author preprint](https://arxiv.org/abs/2306.16583).
7. J. T.-Y. Wang and Z. Xiao, “Families of unit equations and exponential
   Diophantine problems via integral points,” arXiv:2604.26497 (2026):
   [Theorems 1.1, 1.6, and
   1.13](https://arxiv.org/abs/2604.26497).
8. Y. Yasufuku, “GCD inequalities arising from codimension-2 blowups,”
   *Bulletin of the London Mathematical Society* 58 (2026), no. 4,
   e70349: [Theorems 1.2 and
   3.2](https://doi.org/10.1112/blms.70349).
9. F. Adiceam and V. Shirandami, “Probabilistic effectivity in the
   Subspace Theorem,” *Research in Number Theory* 12 (2026), article 8:
   [Theorems 1.3 and 1.6 and Section
   1.1](https://doi.org/10.1007/s40993-025-00692-0).
10. P. Corvaja and U. Zannier, “A lower bound for the height of a rational
    function at S-unit points,” *Monatshefte für Mathematik* 144 (2005),
    203--224:
    [Proposition 2](https://arxiv.org/abs/math/0311030).
