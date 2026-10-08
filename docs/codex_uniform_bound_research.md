# Research continuation toward the uniform endpoint bound

## Status

The endpoint theorem is **not proved** here. This note records the progress
obtained after the quantitative audit and isolates the narrower remaining
problem.

The subsequent [verified continuation](endpoint_continuation_status.md)
adds a direct uniform elimination of all multiplicatively dependent endpoint
pairs and a rational-ray obstruction for the complete `m = 7` hypersimplex
model with arbitrary Gaussian block values. Its unit and subset-extraction
audits explain why this does not yet give a general endpoint bound.

All heights and local absolute values in this note use the absolute
normalization from the audit. Thus, if
\(p=\pi\overline\pi\) is a split rational prime, then

\[
\lVert\pi\rVert_\pi=p^{-1/2},
\qquad
\ell_\pi=-\log\lVert\pi\rVert_\pi=\frac12\log p.
\]

The main positive advance is that two explicit integral bases give one
point-independent product system at every finite place. Conjugate-place
pairing raises the local coefficient from the beta-required value \(150\) to
\(156\). Conditional on a support-independent exceptional theorem, the GCD
coefficient would therefore be \(2/9\) in the zero-margin limit, rather than
\(4/9\).

The main negative advance is equally precise. The remaining Mahler weight
cannot be made common even for two explicit endpoint points from the same
conductor in this natural system. A simple determinant cannot average the
point-dependent weights: it acquires a nonnegative assignment defect, with
equality only when all row weights agree.

A second positive advance comes from the circle determinant. Its local
valuation is an exact nonnegative cut defect.  A seven-row
discriminant-code argument eliminates every squarefree simplex realization
once \(R>(C/2)^6\).  The broader combinatorial obstruction is now sharp:
complete hypersimplex incidence matrices give genuine common-conductor
families of arbitrarily many valuation profiles with exact zero defect and
Plotkin equality, while no bounded row restriction supplies the more than
two disjoint coordinate collisions required by the phase-rigid lemma.
Direct Gaussian separation has unweighted optimum
\((m+1)/(2m)>1/2\), and with the actual near-equal prime weights every
relation remains on the wrong side of \(1/2\). The unresolved input is
therefore a simultaneous arithmetic estimate for this phase kernel, not a stability
theorem for the valuation code alone.

The exact value-level reduction is a positive integral Plücker system with
small cofactors. In the first case \(m=5\), its Gaussian origin supplies a
rank-one Hermitian lift (equivalently a rank-two integral real Gram lift) and
ten compatible near-square identities.  The trivial one of its three cubic
phase branches is impossible for large radius; the other two carry explicit
coprime Pell residuals.  Their individual square-root gap is still short by
one full private-prime factor. Fixed-degree polynomial block families cannot
realize the system at the endpoint, while no available theorem controls
arbitrary conductor-prime values. Standard CRT decoding misses the projective
message class; even if that gap is granted away, its one-sided Johnson
denominator is \(1/(4m^2)\) and its list bound is \(m^2\). Using both
conjugate prime-ideal coordinates improves the formal ceiling to \(m+1\), but
the complete fair weak-flip family has exactly \(m+1\) words and attains
equality. The 2026 rational-code extension handles bad denominator primes,
but its largest possible deterministic radius misses the endpoint. Artificial
interleaving also fails by an exact joint-support calculation, and the
beyond-radius theorems are probabilistic rather than worst-case.

## 1. One fixed two-basis product system

Let

\[
V=\{f\in\mathbf Q(i)[X_0,X_1,X_2]_6:f(1,1,1)=0\},
\qquad \dim V=27.
\]

For

\[
M_{a,b}=X_0^{6-a-b}X_1^aX_2^b,
\qquad R_a=X_1^aX_2^{6-a},
\]

define

\[
\begin{aligned}
\mathcal B_{12}
 &=\{R_a-R_{a+1}:0\leq a\leq5\}
   \cup\{M_{a,b}-R_a:a+b\leq5\},\\
\mathcal B_0
 &=\{M-X_0^6:M\text{ a sextic monomial},\ M\ne X_0^6\}.
\end{aligned}
\tag{1.1}
\]

Each list is a \(\mathbf Z\)-basis of the coefficient-sum-zero lattice in
\(V\). Regard the 28 sextic monomials as vertices. The differences in either
list are the oriented edges of a spanning tree, and edge differences in a
spanning tree form a basis of the augmentation lattice. The transition
determinant is consequently \(\pm1\).

Orient the split part of a Gaussian conductor \(G\), with
\((G,\overline G)=1\). Use \(\mathcal B_{12}\) at every
\(\pi\mid G\), and \(\mathcal B_0\) at every
\(\overline\pi\mid\overline G\) and at all other places. This placewise
tuple depends on the orientation of \(G\), but not on the point

\[
P=[1:A/\overline A:B/\overline B],\qquad A,B\mid G.
\]

Fix \(\pi\mid G\), and write

\[
\alpha=v_\pi(A/\overline A),\qquad
\beta=v_\pi(B/\overline B),\qquad
m=\max(\alpha,\beta),\qquad
\ell_\pi=-\log\lVert\pi\rVert_\pi.
\]

At \(\pi\), the 21 interior edges of \(\mathcal B_{12}\) contribute at
least \(35(\alpha+\beta)\ell_\pi\). The six boundary edges have the form

\[
X_1^aX_2^{5-a}(X_2-X_1),\qquad 0\leq a\leq5,
\]

and contribute at least

\[
\bigl(15(\alpha+\beta)+6\min(\alpha,\beta)\bigr)\ell_\pi.
\]

Thus

\[
\sum_{f\in\mathcal B_{12}}\lambda_{f,\pi}(P)
\geq\bigl(56(\alpha+\beta)-6m\bigr)\ell_\pi.
\tag{1.2}
\]

At \(\overline\pi\), the primitive coordinate valuation vector is
\((m,m-\alpha,m-\beta)\). Across the 27 nonroot monomials, the monomial
exponent totals in \(X_0,X_1,X_2\) are \(50,56,56\), respectively.
Therefore

\[
\sum_{f\in\mathcal B_0}\lambda_{f,\overline\pi}(P)
\geq\bigl(162m-56(\alpha+\beta)\bigr)\ell_\pi.
\tag{1.3}
\]

Adding (1.2) and (1.3) gives the point-independent paired estimate

\[
\boxed{
\sum_{f\in\mathcal B_{12}}\lambda_{f,\pi}(P)
+\sum_{f\in\mathcal B_0}\lambda_{f,\overline\pi}(P)
\geq156m\ell_\pi.
}
\tag{1.4}
\]

Possible cancellation in a boundary factor only increases the left side.
The sum of the three coordinate-divisor proximities over the same conjugate
pair is \(3m\ell_\pi\). Hence (1.4) contains the
\(50\cdot3m\ell_\pi\) required by the finite beta computation, with an
unconditional surplus \(6m\ell_\pi\).

Globally,

\[
\sum_v\sum_{i=1}^{27}\lambda_{L_i^{(v)},v}(P)
\geq156h_H(P)+O(1).
\tag{1.5}
\]

All forms have fixed integral coefficients, there are at most 54 distinct
forms, and

\[
\Delta_{\mathcal L}=1
\]

in Evertse--Ferretti 2013, equation (2.11). Thus neither a local basis index
nor a determinant-height parameter depends on the conductor primes.

### 1.1 The same identity in every degree

The sextic level is not minimal once conjugate places are paired. For

\[
V_d=\{f\in K[X_0,X_1,X_2]_d:f(1,1,1)=0\},
\qquad n_d=\binom{d+2}{2}-1,
\]

define \(\mathcal B_{12}^{(d)}\) and \(\mathcal B_0^{(d)}\) by the same
boundary-tree and root-star constructions as (1.1). Put

\[
S_d=\binom{d+2}{3}.
\]

The identical exponent count gives

\[
W_{12}^{(d)}(\alpha,\beta)
\geq S_d(\alpha+\beta)-d\max(\alpha,\beta),
\tag{1.6}
\]

and

\[
W_0^{(d)}(\alpha,\beta)
\geq n_d\,d\max(\alpha,\beta)-S_d(\alpha+\beta).
\tag{1.7}
\]

Consequently

\[
\boxed{
W_{12}^{(d)}+W_0^{(d)}
\geq d(n_d-1)\max(\alpha,\beta).
}
\tag{1.8}
\]

The cancellation of \(S_d(\alpha+\beta)\) is exact. The sextic value is
\(6(27-1)=156\), but already

\[
d=2,\qquad n_2=5,\qquad d(n_2-1)=8
\tag{1.9}
\]

has strict endpoint margin. Degree \(1\) gives equality at the endpoint and
no positive Subspace-Theorem margin, so \(d=2\) is minimal for this paired
construction.

## 2. The improved conditional GCD coefficient

The section embedding has

\[
h_L(P)=6h_H(P)-h_E(P)+O(1).
\]

If the fixed product theorem gave, outside uniformly bounded exceptional
subspaces,

\[
156h_H(P)\leq(27+\epsilon)h_L(P)+O(1),
\]

then

\[
h_E(P)\leq
\left(6-\frac{156}{27+\epsilon}\right)h_H(P)+O(1).
\tag{2.1}
\]

The exact limiting and endpoint margins are

\[
6-\frac{156}{27}=\frac29,
\]

and

\[
6-\frac{156}{27+\epsilon}<\frac12
\quad\Longleftrightarrow\quad
\epsilon<\frac{15}{11}.
\tag{2.2}
\]

Evertse--Ferretti 2013, Corollary 3.2, assumes
\(0<\epsilon\leq1\). Its largest permitted value already gives

\[
6-\frac{156}{28}=\frac37<\frac12.
\]

Equivalently, the endpoint condition gives
\(h_L\leq(11/2)h_H+O_C(1)\), while

\[
\frac{156}{11/2}=27+\frac{15}{11}.
\]

The product inequality therefore has comfortable numerical slack. The
published Corollary 3.2 still covers it by

\[
(6561\epsilon^{-1})^{27|S|}
10^{10}2^{54}27^{15}\epsilon^{-3}
\log(3\epsilon^{-1})
\log(\epsilon^{-1}\log3)
\tag{2.3}
\]

subspaces. The exponent \(27|S|\), not the analytic coefficient, is the
remaining quantitative obstruction on this route.

### 2.1 The smallest-dimensional quantitative specialization is quadratic

In general degree, (1.8) would give

\[
\eta_d(\epsilon)
=d-\frac{d(n_d-1)}{n_d+\epsilon}.
\tag{2.4}
\]

The exact endpoint condition is

\[
\eta_d(\epsilon)<\frac12
\quad\Longleftrightarrow\quad
\epsilon<
\frac{d(d-1)}{4d-2}.
\tag{2.5}
\]

For \(d=2\), take \(\epsilon=1/4<1/3\). Then

\[
\eta_2(1/4)
=2-\frac8{5+1/4}
=\frac{10}{21}<\frac12.
\tag{2.6}
\]

The exact Evertse--Ferretti 2013, Corollary 3.2 count becomes

\[
\boxed{
900^{5|S|}
10^{10}2^{16}5^{15}
\log 12\,\log(4\log3).
}
\tag{2.7}
\]

This is exponentially smaller than (2.3), but still exponentially dependent
on \(|S|\). The complete system \(|2H-E|\) embeds the one-point blowup as
the cubic scroll in \(\mathbf P^4\). A proper exceptional linear subspace is
contained in a hyperplane, whose pullback pushes down to a conic through
\([1:1:1]\). Thus a hypothetical support-independent bound of \(T\)
subspaces at this quadratic level would give total plane degree at most
\(2T\), and the formal grid endgame would give at most \(2T+2\) points in
the original endpoint cluster.

## 3. Even the quadratic system has no common twisted weight

The weight obstruction persists at the minimal degree. Let

\[
T=30q+10,\qquad H_1=T+i,\qquad H_2=T+3i,
\qquad G_T=H_1H_2,
\]

and

\[
u_1=H_1/\overline H_1,\qquad
u_2=H_2/\overline H_2,\qquad
P_T=[1:u_1:u_2],\quad P'_T=[1:u_2:u_1].
\tag{3.1}
\]

The four factors \(H_1,H_2,\overline H_1,\overline H_2\) are pairwise
coprime. At primes dividing them, the passive coordinate residues are
\(-1/2,-2,1/2,2\), respectively. Parity and
\(T\equiv1\pmod 3\) exclude the only residue characteristics that can
create cancellation in the quadratic basis evaluations. In particular the
profiles below are exact, including prime multiplicity.

For \(\mathcal B_{12}^{(2)}\), the two swapped profiles have sums and
absolute difference

\[
\sum a_i=\sum b_i=2,
\qquad
\sum_i|a_i-b_i|=4.
\tag{3.2}
\]

For \(\mathcal B_0^{(2)}\), the reflected profiles have

\[
\sum a'_i=\sum b'_i=6,
\qquad
\sum_i|a'_i-b'_i|=6.
\tag{3.3}
\]

For arbitrary \(a,b\in\mathbf R^n\), put
\(\phi_c(a)=\min_i(a_i+c_i)\). The exact identity

\[
\sup_{\sum_i c_i=0}
\bigl(\phi_c(a)+\phi_c(b)\bigr)
=\overline a+\overline b
-\frac1n\min_{r\in\mathbf R}\sum_i|a_i-b_i-r|
\tag{3.4}
\]

follows by writing the two minima as \(s,t\). Feasibility is equivalent to
\(\sum_i\max(s-a_i,t-b_i)\leq0\); now set \(u=s+t\), \(r=s-t\), and use the
formula for the maximum of two real numbers.

The medians in (3.2)--(3.3) are zero. Hence a weight shared by \(P_T\) and
\(P'_T\) has combined local optimum

\[
\frac{4-4}{5}+\frac{12-6}{5}=\frac65
\tag{3.5}
\]

per logarithmic conductor unit, whereas point-dependent weights attain
\(16/5\).

Use \(\mathcal B_{12}^{(2)}\) at oriented conductor primes and
\(\mathcal B_0^{(2)}\) elsewhere. Put

\[
\mathfrak h_T=\log|G_T|
=\frac12\bigl(\log N(H_1)+\log N(H_2)\bigr).
\]

Then

\[
h_H(P_T)=h_H(P'_T)=\mathfrak h_T+O(1),
\qquad
h_{2H-E}(P_T)=h_{2H-E}(P'_T)
=\frac32\mathfrak h_T+O(1).
\tag{3.6}
\]

Each point has product proximity \(8\mathfrak h_T+O(1)\), so both satisfy
Evertse--Ferretti equation (3.11) for every
\(0<\epsilon<1/3\), including \(\epsilon=1/4\).

For a finite-support zero-sum weight tuple, define

\[
F_{\mathbf c,Q}(z)
=\sum_v\min_i\bigl(
\lambda_{L_i^{(v)},v}(z)+c_{iv}\log Q
\bigr).
\]

The \(\mathcal B_0^{(2)}\)-coordinate identity is

\[
\sum_v\sum_{i=1}^5
\lambda_{\mathcal B_{0,i}^{(2)},v}(z)=5h(z).
\]

Using (3.5) at the conductor places and the bound
\(\min_i x_i\leq5^{-1}\sum_i x_i\) elsewhere gives

\[
F_{\mathbf c,Q}(P_T)+F_{\mathbf c,Q}(P'_T)
\leq2h(z_T)-\frac65\mathfrak h_T+O(1).
\tag{3.7}
\]

Coordinate exchange preserves \(h(z_T)\), so at least one point \(P_T^*\)
satisfies

\[
F_{\mathbf c,Q}(P_T^*)
\leq h(z_T)-\frac35\mathfrak h_T+O(1).
\]

Evertse--Ferretti equation (2.13) now yields

\[
\log H_{\mathcal L,\mathbf c,Q}(P_T^*)
\geq\frac35\mathfrak h_T+O(1)\longrightarrow+\infty.
\tag{3.8}
\]

Thus no one weight tuple makes both endpoint points twisted-small, for any
\(Q\geq1\) and any \(\delta>0\). This conclusion concerns the direct fixed
\(\mathcal B_{12}^{(2)}/\mathcal B_0^{(2)}\) tuple. It does not rule out a
new residue-adapted or conductor-dependent embedding.

## 4. Why two natural weight-averaging repairs fail

### 4.1 A bounded weight net cannot cover the full valuation box

Take \(q\) equal-logarithmic conductor blocks and, in each block, choose one
of the two swapped quadratic profiles in (3.2)--(3.3). This gives a binary
cube of \(2^q\) valuation profiles.

For one profile, the point-dependent local mean is \(8/5\). For two
profiles that disagree in one block, (3.5) shows that the largest shared
minimum drops from \(16/5\) to \(6/5\), a loss of \(10/5\). If both points
are to remain twisted-small even at exponent zero at the endpoint height
\(3/2\), their combined minimum must be at least \(3=15/5\) per block on
average. Thus two binary profiles captured by the same weight tuple must
have Hamming distance at most

\[
\frac q{10}.
\tag{4.1}
\]

Any class with diameter at most \(q/10\) lies in a Hamming ball of that
radius around any one of its members. Consequently a weight cover of the
entire valuation cube requires at least

\[
\frac{2^q}{\displaystyle
  \sum_{j\leq\lfloor q/10\rfloor}\binom qj}
=2^{(1-H_2(1/10)+o(1))q}
\tag{4.2}
\]

systems, where \(1-H_2(1/10)=0.5310\ldots\). This is an exponential lower
bound for valuation-only robustification.

It is **not** a counterexample to the endpoint theorem. The \(2^q\) formal
divisor profiles need not all satisfy the simultaneous archimedean
short-arc condition. Equation (4.2) instead says exactly what any positive
argument must use: the archimedean phase compatibility has to eliminate
almost all of the valuation cube.

### 4.2 A mixed evaluation determinant pays a wrong-sign assignment cost

Suppose \(n\) points \(x_r\) have row-dependent zero-sum weights
\(c_{ri}\), and put at one place

\[
M_r=\max_i
\bigl|L_i(x_r)\bigr|_vQ^{-c_{ri}}.
\]

Expanding the evaluation determinant gives, at a nonarchimedean place,

\[
\left|\det(L_i(x_r))\right|_v
\leq Q^{\Gamma_v}\prod_{r=1}^n M_r,
\qquad
\Gamma_v=\max_{\sigma\in S_n}
\sum_{r=1}^n c_{r,\sigma(r)}.
\tag{4.3}
\]

At an archimedean place there is only the additional fixed factor \(n!\).
The average of the permutation sums in (4.3) is zero because every row of
\(c\) sums to zero. Hence

\[
\Gamma_v\geq0.
\tag{4.4}
\]

Moreover, equality holds exactly when all row weight vectors are identical.
Indeed, equality forces every permutation sum to be zero. Comparing two
permutations that differ by a transposition gives

\[
c_{ri}+c_{sj}=c_{rj}+c_{si}.
\]

Thus \(c_{ri}=\alpha_r+\beta_i\); the zero row sums force all
\(\alpha_r\) to be equal, so the rows coincide. The estimate is sharp on a
permutation evaluation matrix.

Conjugate places do not cancel the defect. Replacing \(c\) by \(-c\) gives

\[
\Gamma_v(c)+\Gamma_{\overline v}(-c)
=\max_\sigma C(\sigma)-\min_\sigma C(\sigma)\geq0,
\tag{4.5}
\]

again with equality only in the aligned case. Therefore the direct
multipoint determinant cannot exchange the quantifiers
\(\forall P\,\exists\mathbf c(P)\) into one useful determinant estimate:
it pays a nonnegative error precisely when the weights differ.

In the quadratic basis, after fixed column operations, an evaluation row is

\[
(X,Y,X^2,XY,Y^2),\qquad X=x-1,\quad Y=y-1.
\]

Every five-row determinant has total tangent degree eight, exactly the paired
conductor coefficient in (1.9). Each of its five columns can be used only
once in a determinant term, so taking more points creates no hidden copy of
the local gain.

## 5. The Ramana determinant gives an exact conductor-defect statistic

Let \(M=2s+1\) distinct Gaussian integers \(z_i\) have common norm \(N\).
Ramana's Lemma 2.1 and Definition 2.1, specialized to \(k=s\), give a
Gaussian integer determinant \(J_s\) satisfying

\[
N_{\mathbf Q(i)/\mathbf Q}(J_s)
=\frac{\prod_{i<j}N(z_j-z_i)}{N^{s^2}}.
\tag{5.1}
\]

Fix a split rational prime \(p=\pi\overline\pi\) with \(p^e\Vert N\), and
write

\[
a_i=v_\pi(z_i),\qquad 0\leq a_i\leq e.
\]

Since \(v_{\overline\pi}(z_i)=e-a_i\),

\[
v_p\bigl(N(z_j-z_i)\bigr)\geq e-|a_i-a_j|.
\]

Consequently \(p^{D_p}\mid N(J_s)\), where

\[
D_p=e\,s(s+1)-\sum_{i<j}|a_i-a_j|.
\tag{5.2}
\]

More precisely, if

\[
\kappa_\pi=\sum_{i<j}
\bigl(v_\pi(z_j-z_i)-\min(a_i,a_j)\bigr)
\]

and \(\kappa_{\overline\pi}\) is defined conjugately, then the determinant
identity gives the exact valuation

\[
v_\pi(J_s)+v_{\overline\pi}(J_s)
=D_p+\kappa_\pi+\kappa_{\overline\pi},
\qquad \kappa_\pi,\kappa_{\overline\pi}\geq0.
\tag{5.2a}
\]

For the threshold cuts

\[
r_t=\#\{i:a_i\geq t\},\qquad 1\leq t\leq e,
\]

the identity

\[
\sum_{i<j}|a_i-a_j|=\sum_{t=1}^e r_t(2s+1-r_t)
\]

turns (5.2) into

\[
\boxed{
D_p=\sum_{t=1}^e(r_t-s)(r_t-s-1).
}
\tag{5.3}
\]

Every summand is a nonnegative even integer. It vanishes precisely for
\(r_t=s\) or \(r_t=s+1\). If \(D_p=0\), the multiset of valuations is

\[
\{\underbrace{0,\ldots,0}_{s},\ k,
  \underbrace{e,\ldots,e}_{s}\}
\quad\text{for some }0\leq k\leq e,
\tag{5.4}
\]

with the evident multiplicity merger when \(k=0\) or \(e\).

If all points lie on an arc of length at most \(C N^{1/4}\), every chord
norm is at most \(C^2N^{1/2}\). Equation (5.1) therefore gives

\[
N(J_s)\leq C^{2s(2s+1)}N^{s/2}.
\]

Combining the upper and lower bounds first gives the residue-sensitive
inequality

\[
\sum_{p\ {\rm split}}
\bigl(D_p+\kappa_{\pi,p}+\kappa_{\overline\pi,p}\bigr)\log p
\leq \frac{s}{2}\log N+2s(2s+1)\log C.
\tag{5.5a}
\]

Any determinant valuations away from the conductor only strengthen the
left side. Dropping the nonnegative \(\kappa\)-terms yields the cut-only
stability inequality

\[
\boxed{
\sum_{p\ {\rm split}}\sum_{t=1}^{e_p}
  (r_{p,t}-s)(r_{p,t}-s-1)\log p
\leq \frac{s}{2}\log N+2s(2s+1)\log C.
}
\tag{5.5}
\]

For \(M=2s\), Ramana's Lemma 2.1 with its parameters
\(m=2s,\ k=s-1\) gives

\[
N^{s(s-1)}
N_{\mathbf Q(i)/\mathbf Q}
\bigl(\det_{s-1}(\mathbf z,\overline{\mathbf z})\bigr)
=\prod_{i<j}N(z_j-z_i).
\]

The resulting layer defect is
\(s^2-r_t(2s-r_t)=(r_t-s)^2\), and the residual exponent on the right of
(5.5) remains \(s/2\). Thus an endpoint cluster constrains the aggregate
weighted squared imbalance of its conductor cuts. It does not force every
cut, or even most unweighted cuts, to be close to half-and-half.

### 5.1 The determinant bound is exactly the weighted Plotkin bound

This limitation is algebraic, not an artifact of estimating the Ramana
determinant. Index all prime-exponent layers by \(\ell\), give layer \(\ell\)
weight \(w_\ell\), and put

\[
W=\sum_\ell w_\ell,
\qquad
d(i,j)=\sum_\ell w_\ell
\mathbf 1_{\{\ell\text{ separates }i,j\}}.
\]

For \(2s+1\) rows, a layer whose cut has size \(r_\ell\) separates
\(r_\ell(2s+1-r_\ell)\) unordered pairs. Therefore

\[
\boxed{
\sum_{i<j}d(i,j)
=s(s+1)W-D,
\qquad
D=\sum_\ell w_\ell(r_\ell-s)(r_\ell-s-1).
}
\tag{5.6}
\]

After removing factors common to every point, the endpoint chord estimate
for each pair is

\[
d(i,j)\geq\frac W2-2\log C.
\tag{5.7}
\]

After first removing every factor common to all points, and every
self-conjugate factor (which only rescales the common circle), one has
\(W=\log N\).  Summing (5.7) over all pairs and substituting (5.6) then
gives exactly

\[
D\leq\frac s2W+2s(2s+1)\log C,
\tag{5.8}
\]

which is (5.5). Thus the multipoint determinant supplies no amplification
beyond the sum of the pairwise endpoint inequalities.

The equality structure is nevertheless exact. At a prime of exponent
\(e\), its threshold cuts are nested. If their total defect vanishes, every
cut has size \(s\) or \(s+1\); nested equal-cardinality cuts coincide, so
there is at most one strict transition. Up to relabeling, the valuation
multiset is

\[
\{\underbrace{0,\ldots,0}_{s},a,
  \underbrace{e,\ldots,e}_{s}\},
\qquad 0\leq a\leq e.
\tag{5.9}
\]

Since every nonzero layer defect is an even integer at least two, local
defect \(K\) also permits at most \(K/2\) nonbalanced thresholds. This
laminar rigidity does not globalize for squarefree conductors: each prime
then contributes only one layer, and those balanced cuts may be unrelated.

The algebraic identities behind (5.3), its equality case and parity, as well
as the seven-row matrix identity, mod-\(4\) kernel parity, and support-cardinality
classification used below, are compiler-checked in
[GaussianChain/ConductorDefect.lean](../GaussianChain/ConductorDefect.lean).

## 6. Valuation balance is insufficient; the modular phase gap is exact

For \(M=2^r-1\), index rows and columns by the nonzero vectors of
\(\mathbf F_2^r\), and put

\[
b_{x,a}=a\mathbin\cdot x\pmod2.
\]

Every column has \(2^{r-1}\) ones and \(2^{r-1}-1\) zeroes, so every
squarefree conductor layer has zero defect in (5.3). Any two rows differ in
exactly \(2^{r-1}\) columns. This constructs arbitrarily large valuation
patterns that evade every argument based only on the nonnegative defects
\(D_p\).

These patterns are realized by genuine Gaussian points on a common circle.
Choose distinct oriented split Gaussian primes \(\pi_a\), put

\[
G=\prod_{a\ne0}\pi_a,\qquad
A_x=\prod_{a\ne0}\pi_a^{b_{x,a}},\qquad
z_x=A_x\,\overline{G/A_x}.
\tag{6.1}
\]

Then \(N(z_x)=N(G)\), the points are distinct, and their conductor
valuations are exactly the simplex columns. This valuation realization alone
does not address the endpoint arc condition; Section 6.1 below shows that the
same patterns cannot occupy one such arc at unbounded radius.

There is no automatic residue surplus hidden in (5.5a). For example, take

\[
(\pi_a)_{a=1}^7
=(10+7i,\,3+2i,\,5+4i,\,8+3i,\,6+5i,\,11+6i,\,9+4i).
\tag{6.1a}
\]

Their norms
\(149,13,41,73,61,157,97\) are distinct rational primes. Exact Gaussian
division on the \(21\) pair differences in the seven-row simplex gives

\[
v_{\pi_a}(z_y-z_x)=\min(b_{x,a},b_{y,a}),\qquad
v_{\overline\pi_a}(z_y-z_x)
=\min(1-b_{x,a},1-b_{y,a})
\]

for every \(a\) and every distinct \(x,y\). Thus every
\(\kappa_{\pi_a}=\kappa_{\overline\pi_a}=0\) in this genuine configuration.

The archimedean condition can be written exactly. Let

\[
\omega_a=\pi_a/\overline\pi_a,\qquad
\theta_a=\arg\omega_a,
\]

and let \(B=(b_{x,a})\). If all \(u_x=\prod_a\omega_a^{b_{x,a}}\) lie in
an angular interval of width
\(\Delta=C\,N(G)^{-1/4}\), there are \(c\in\mathbf R\),
\(m\in\mathbf Z^M\), and \(\varepsilon\in\mathbf R^M\), with
\(\|\varepsilon\|_\infty\leq\Delta\), such that

\[
B\theta=c\mathbf1+2\pi m+\varepsilon.
\tag{6.2}
\]

Put \(H_{x,a}=(-1)^{a\cdot x}\). The punctured Walsh matrices satisfy

\[
B^{-1}=-\frac{2}{M+1}H.
\]

Since \(H\mathbf1=-\mathbf1\), equation (6.2) implies, for
\(k=(M+1)/2\),

\[
k\theta_a
=c-2\pi(Hm)_a-(H\varepsilon)_a.
\]

Therefore

\[
\left|k\theta_a-c\pmod{2\pi}\right|
\leq M\Delta.
\tag{6.3}
\]

This is the correct Fourier consequence. The integers \(m_x\) cannot be
dropped: the original prime angles need not be close; only the growing
powers \(\omega_a^k\) must align. The crude one-ratio algebraic separation

\[
\left|\left(\frac{\omega_a}{\omega_b}\right)^k-1\right|
\gg (\lvert\pi_a\rvert\lvert\pi_b\rvert)^{-k}
\tag{6.4}
\]

is too weak at \(k=(M+1)/2\) to contradict the endpoint width in general.
It does not, however, use the full congruence structure in (6.2). That
structure eliminates the exact simplex pattern uniformly.

### 6.1 A phase-rigid matrix lemma

There is a useful abstract criterion. Choose \(t\) points and suppose the
points have first been restricted to one of the four Gaussian-unit classes,
which costs at most a factor four in the eventual cardinality bound. Suppose
also that the nonconstant conductor factors can be grouped into pairwise coprime nonunit
Gaussian blocks \(H_a\), with
\(\gcd(H_a,\overline H_b)=1\) for every \(a,b\), including \(a=b\). This is
automatic for an oriented split squarefree conductor satisfying
\((G,\overline G)=1\), and more generally for binary endpoint
profiles in which the whole prime power belongs to one block; it is not
automatic after splitting the threshold layers of a higher-power prime.
Suppose the resulting \(t\) cuts are the columns of an
invertible matrix \(A\in\{0,1\}^{t\times t}\), and

\[
A^{-1}\mathbf1=\lambda\mathbf1.
\]

Put

\[
\Lambda_A=A^{-1}\mathbf Z^t/\mathbf Z^t\subset(\mathbf Q/\mathbf Z)^t
\]

and define \(\mu(A)\) as the minimum, over
\(\gamma\in\Lambda_A\), of the largest number of disjoint coordinate pairs
\((a,b)\) satisfying \(\gamma_a=\gamma_b\). Also put

\[
K(A)=\frac12\max_{a,b}
\|\operatorname{row}_a(A^{-1})-
  \operatorname{row}_b(A^{-1})\|_1.
\tag{6.5}
\]

If the point phases lie in an interval of width \(\Delta\), midpoint lifts
give

\[
A\theta=c\mathbf1+2\pi m+\varepsilon,
\qquad \|\varepsilon\|_\infty\leq\Delta/2.
\]

For \(\gamma=A^{-1}m\pmod{\mathbf Z^t}\), every equal-coordinate pair in
\(\gamma\) then satisfies

\[
|\theta_a-\theta_b|_{\mathbf R/2\pi\mathbf Z}
\leq K(A)\Delta.
\tag{6.6}
\]

Since the blocks are pairwise coprime together with their conjugates,

\[
\left|\frac{H_a}{\overline H_a}
      -\frac{H_b}{\overline H_b}\right|
=\frac{|H_a\overline H_b-\overline H_aH_b|}{|H_aH_b|}
\geq\frac2{|H_aH_b|}.
\tag{6.7}
\]

For \(\mu=\mu(A)\) disjoint pairs, (6.6)--(6.7) imply

\[
R_{\rm var}:=\prod_a|H_a|
\geq\left(\frac2{K(A)\Delta}\right)^\mu.
\tag{6.8}
\]

Blocks constant on the selected points contribute a common factor \(U\),
so the full circle radius is
\(R=|U|R_{\rm var}\geq R_{\rm var}\). Hence

\[
\Delta\geq\frac2{K(A)}R^{-1/\mu}.
\]

Combining this with the endpoint upper bound
\(\Delta\leq C R^{-1/2}\) proves, whenever \(\mu(A)>2\),

\[
\boxed{
R\leq\left(\frac{C K(A)}2\right)^{2\mu(A)/(\mu(A)-2)}.
}
\tag{6.9}
\]

The relevant invariant is the *embedded* discriminant code
\(\Lambda_A\), not just the Smith invariant factors of \(A\): the proof
needs coordinate collisions in every codeword.

### 6.2 Seven rows rule out every large squarefree simplex realization

For the seven-point simplex incidence matrix \(B\), put

\[
\mathcal H_{x,a}=(-1)^{a\cdot x},
\qquad \mathcal H=J-2B.
\]

Then

\[
B^{-1}=-\frac14\mathcal H,\qquad
\operatorname{SNF}(B)=\operatorname{diag}(1,1,1,2,2,2,4),
\qquad K(B)=1.
\tag{6.10}
\]

The \(32\) elements of \(\Lambda_B\) can be classified without enumeration.
Writing \(\gamma=-\mathcal Hm/4\), the congruence
\(B(-\mathcal Hm)=0\pmod4\) holds. Since every entry of \(\mathcal H\)
is odd, the seven coordinates of \(-\mathcal Hm\) all have parity
\(\sum_xm_x\pmod2\). Write this numerator as
\(c_0\mathbf1+2e\pmod4\). Since \(B\mathbf1=4\mathbf1\), the remaining
condition is \(Be=0\pmod2\). Its row indexed by \(x\) is
\(x\cdot\bigoplus_{a:e_a=1}a\), so the Boolean support has xor zero in
\(\mathbf F_2^3\), and its size is
\(0,3,4,\) or \(7\). Every codeword therefore has coordinate multiplicities
\(7\), or \(4+3\), and contains three disjoint equal-coordinate pairs. Thus

\[
\mu(B)=3.
\]

Equation (6.9) now gives the sharp output of this argument:

\[
\boxed{R\leq(C/2)^6.}
\tag{6.11}
\]

The same bound holds for every squarefree block-simplex with
\(M=2^r-1\), \(r\geq3\). Restrict its rows
to the seven nonzero vectors in any three-dimensional subspace. Aggregate
the original conductor columns by their restriction to that subspace;
the seven nonzero restriction classes give the blocks \(H_a\), while the
zero class is common to all selected rows. The seven-row proof applies
unchanged.

Thus the exact squarefree simplex patterns demonstrate why valuation defect
alone is insufficient, but they cannot form an unbounded endpoint family.

### 6.3 A second phase-rigid family: cyclic interval cuts

The phenomenon is not special to Walsh matrices. For a squarefree
block system with \(n=2s+1\), index
rows and columns cyclically and put

\[
b_{x,a}=1
\quad\Longleftrightarrow\quad
a\in\{x,x+1,\ldots,x+s-1\}\pmod n.
\tag{6.12}
\]

Every cut has size \(s\), hence zero layer defect. If
\(\phi_x=\sum_{j=0}^{s-1}\theta_{x+j}\) are the point phases and all lie in
an interval of width \(\Delta\), then

\[
|\theta_{x+s}-\theta_x|_{\mathbf R/2\pi\mathbf Z}
=|\phi_{x+1}-\phi_x|_{\mathbf R/2\pi\mathbf Z}
\leq\Delta.
\]

Because \(\gcd(s,2s+1)=1\), the edges \(x\leftrightarrow x+s\) form one
odd cycle and contain \(s\) disjoint pairs. Gaussian separation as in
(6.7) gives

\[
R_{\rm var}\geq(2/\Delta)^s.
\]

For \(s>2\), comparison with \(\Delta\leq C R^{-1/2}\) yields

\[
\boxed{
R\leq(C/2)^{2s/(s-2)}.
}
\tag{6.13}
\]

Thus two extremal zero-defect designs with very different combinatorics are
phase-rigid at endpoint scale.

### 6.4 The rectangular obstruction is an unbounded Plotkin-equality family

The \(7\times35\) example is not isolated.  Let

\[
m=2s+1\geq5,\qquad k=s+1=\frac{m+1}{2},\qquad
L=\binom mk,
\]

index rows by \([m]\), columns by the \(k\)-subsets \(S\subset[m]\), and
put

\[
A_{iS}=\mathbf1_{\{i\in S\}}.
\tag{6.14}
\]

Every column has size \(s+1\), hence zero odd-cluster defect.  If

\[
r=\binom{m-1}{k-1}=L\frac{m+1}{2m},
\qquad
\lambda=\binom{m-2}{k-2}=\frac r2,
\]

then every row has weight \(r\), two rows meet in \(\lambda\) columns, and

\[
d(i,j)=2(r-\lambda)=r>\frac L2,
\qquad
AA^{\mathsf T}=\frac r2(I+J).
\tag{6.15}
\]

Append one all-zero row.  The resulting \(M=m+1\) codewords are
equidistant, every column has exactly \(M/2\) ones, and

\[
d=\frac{LM}{2(M-1)},
\tag{6.16}
\]

which is equality in the even binary Plotkin bound.  In coding terminology
this is the fair weak-flip code: Lin--Moser--Chen, Definition 16, uses every
weak-flip column equally often; Lemma 26, equation (41), is the
parity-sensitive Plotkin bound; Corollary 40, equation (63), makes every
\(q\)-wise distance equal; Theorem 43, equation (71), is the generalized
\(q\)-wise Plotkin bound; and Theorem 45, equation (79), with Corollary 46
propagates equality through those bounds. Their Example 20 is
exactly \((M,L,d)=(8,35,20)\), the former \(7\times35\) matrix with its zero
row.  Thus neither exact pairwise Plotkin equality nor equality in all of
their higher-distance bounds forces a Walsh/simplex block.

This is a genuine common-conductor valuation family.  Choose distinct
rational primes \(p_S\equiv1\pmod4\), oriented factors
\(p_S=\pi_S\overline\pi_S\), and put

\[
G=\prod_S\pi_S,\qquad A_0=1,\qquad
A_i=\prod_{S\ni i}\pi_S,\qquad
z_i=A_i\overline{G/A_i}.
\tag{6.17}
\]

Then \(N(z_i)=N(G)\), the \(M=m+1\) points are distinct, and their exact
prime-valuation cuts are the columns above.  For fixed \(m\), the prime
number theorem in the progression \(1\pmod4\) lets the \(L\) primes be
chosen in \([X,X^{1+\eta}]\).  If \(\eta\leq1/m\), with
\(w_S=\log p_S\) and \(W=\sum_Sw_S\), then every pair satisfies

\[
d_w(i,j)\geq r\log X
\geq\frac{L(1+\eta)}2\log X
\geq\frac W2.
\tag{6.18}
\]

So the family obeys the actual conductor identities, has zero layer defect,
and satisfies the endpoint necessary pair inequalities with no
\(O_C(1)\) loss.  Equation (6.18) is not an assertion that its phases lie in
one endpoint arc.

### 6.5 No bounded row restriction restores coordinate phase rigidity

Return here to the original \(m\) nonzero rows; the appended zero row in
Section 6.4 was used to identify the Plotkin-equality code.  The integral
failure also scales.  Let

\[
Y=\{y\in\mathbf Z^m:\sum_i y_i=0\},\qquad
\Phi(y)_S=\sum_{i\in S}y_i.
\tag{6.19}
\]

Up to a unimodular choice of independent row differences, \(\Phi\) is
\(D_A^{\mathsf T}\).  Its image is saturated.  Indeed, if
\(y\in Y\otimes\mathbf Q\) and every \(k\)-subset sum is integral, swapping
\(i\) and \(j\) in a common \((k-1)\)-subset shows
\(y_i-y_j\in\mathbf Z\).  Write \(y_i=c+n_i\).  A \(k\)-subset sum gives
\(kc\in\mathbf Z\), while \(\sum_i y_i=0\) gives
\(mc\in\mathbf Z\).  Since \(\gcd(m,k)=1\), one has
\(c\in\mathbf Z\).

No coordinate difference \(e_U-e_V\) lies even in the rational image of
\(\Phi\).  If it did, then for fixed \(i\ne j\), choose a
\((k-1)\)-subset \(T\) disjoint from \(i,j\) for which neither
\(T\cup\{i\}\) nor \(T\cup\{j\}\) is \(U\) or \(V\).  There are
\(\binom{m-2}{k-1}>2\) choices and at most two are excluded.  Subtracting
the two unaffected coordinates gives \(y_i=y_j\).  Hence all \(y_i\) are
equal, and the zero-sum condition gives the contradiction \(y=0\).
Therefore

\[
\operatorname{coker}(D_A^{\mathsf T})
\quad\text{is torsion-free on every coordinate-difference class.}
\tag{6.20}
\]

If the zero row is included in the row-difference matrix, the full cokernel
has torsion subgroup exactly \(\mathbf Z/k\mathbf Z\). Indeed, if
\(y\in\mathbf Q^m\) has all \(k\)-subset sums integral, then
\(y_i=c+n_i\) with \(n_i\in\mathbf Z\) and \(kc\in\mathbf Z\); conversely
\(A^{\mathsf T}(k^{-1}\mathbf1)=\mathbf1\), and
\(t\mathbf1\in A^{\mathsf T}\mathbf Z^m\) only when \(k\mid t\).
Nevertheless, the unaffected-swap argument still proves
that no coordinate difference lies in its rational row span.  Thus that
global torsion supplies no pair to which (6.6) applies.

Nor can one select a bounded number of rows and then use a hidden rigid
minor.  Select \(t\) rows.  If \(m\geq2t+3\), every restriction profile
\(R=S\cap[t]\) occurs at least twice, since its multiplicity is

\[
\binom{m-t}{k-|R|}\geq2.
\]

Every selected-row combination is constant on a profile fibre, so no
original coordinate difference belongs to its rational row span.  The only
legitimate compression is to multiply all prime blocks with the same
profile.  After deleting the two constant profiles, the aggregate map is

\[
\Psi(y)_R=\sum_{i\in R}y_i,\qquad
\varnothing\ne R\subsetneq[t],\qquad \sum_i y_i=0.
\tag{6.21}
\]

It is saturated because singleton coordinates recover every \(y_i\), and
\(\Psi(y)_{R^c}=-\Psi(y)_R\).  A coordinate difference in its image would
therefore have to be \(e_R-e_{R^c}\); singleton coordinates rule this out
for \(t\geq3\).  For \(t=2\) there is exactly one such pair, whereas the
phase-rigid bound (6.9) needs more than two disjoint pairs. Thus no bounded
restriction among the nonzero rows satisfies the present coordinate-torsion
criterion.

The same conclusion holds when the selected rows include the appended zero
row. If there are also \(t\) nonzero rows and \(m\geq2t+3\), compression by
restriction profile gives

\[
\Psi_0(y)_R=\sum_{i\in R}y_i,
\qquad \varnothing\ne R\subseteq[t],
\qquad y\in\mathbf Z^t.
\tag{6.21a}
\]

Every profile fibre again has multiplicity at least two. The map is
saturated because its singleton coordinates recover the \(y_i\). If a
coordinate difference belonged to its rational image, those singleton
coordinates would force \(y\) to have at most two nonzero entries. For
\(t\geq3\), either one entry is nonzero, producing at least
\(2^{t-1}\) nonzero profile coordinates, or
\(y=e_i-e_j\), producing \(2^{t-2}\) positive and \(2^{t-2}\) negative
coordinates. Neither is a coordinate difference. For \(t=2\), precisely
\(e_{\{1\}}-e_{\{2\}}\) occurs, again only one pair. Thus the assertion
really covers every bounded restriction of the full \(m+1\)-row family.

The columns in (6.14) are the bases of the uniform matroid \(U_{k,m}\), so
the phase kernel is the hypersimplex toric kernel.  Lasoń--Michałek,
Theorem 2, proves that the toric ideal is generated by quadratic symmetric
exchange binomials for strongly base-orderable matroids, including uniform
matroids; their Theorem 3 proves White's conjecture up to saturation for
every matroid.  These abundant small *kernel* trades explain the large real
phase kernel.  They do not create a two-supported *row-space* vector of the
kind required by (6.6).

### 6.6 The exact limit of row-relation Gaussian separation

There is a general positive criterion, but the family above lies on its
wrong side.  First restrict to one of the four Gaussian-unit classes; this
costs at most a factor four in a cardinality bound and makes the unit factor
common to all selected points.  Suppose their phases \(\phi_i\) occupy an
interval of width \(\Delta\), the oriented blocks satisfy the coprimality
hypothesis of Section 6.1, and normalize by one selected point.  A
zero-total-exponent relation on all selected points is then encoded by an
arbitrary \(y\in\mathbf Z^m\) on the remaining points, with coefficient
\(-\sum_i y_i\) on the normalized point. Put

\[
f_y(a)=\sum_i y_iA_{ia},\qquad
H_y=\prod_{f_y(a)>0}H_a^{f_y(a)}
    \prod_{f_y(a)<0}\overline H_a^{-f_y(a)}.
\]

Assume \(f_y\ne0\), and write

\[
\widetilde y=(-\textstyle\sum_i y_i,y_1,\ldots,y_m).
\]

The common unit cancels because the entries of \(\widetilde y\) sum to
zero. The coprimality hypothesis and \(f_y\ne0\) imply
\(\alpha_y=H_y/\overline H_y\ne1\), and midpoint lifts give the exact
two-sided estimate

\[
\frac2{|H_y|}
\leq|\alpha_y-1|
\leq\frac{\lVert\widetilde y\rVert_1}{2}\Delta.
\tag{6.22}
\]

Consequently, if

\[
|H_y|\leq R^{1/2-\delta}
\]

for some \(\delta>0\), the endpoint bound
\(\Delta\leq C R^{-1/2}\) forces

\[
\boxed{
R\leq\left(\frac{C\lVert\widetilde y\rVert_1}{4}\right)^{1/\delta}.
}
\tag{6.23}
\]

More generally, multiply (6.22) for relations
\(y^{(1)},\ldots,y^{(h)}\).  If

\[
\sum_{\nu=1}^h|f_{y^{(\nu)}}(a)|\leq L_0
\quad\text{for every block }a,
\qquad h>2L_0,
\]

and

\[
Y_0=\left(\prod_{\nu=1}^h
  \frac{\lVert\widetilde y^{(\nu)}\rVert_1}{2}\right)^{1/h},
\]

then

\[
\boxed{
R\leq\left(\frac{CY_0}{2}\right)^{2h/(h-2L_0)}.
}
\tag{6.24}
\]

This extends the phase-rigidity strategy to rectangular matrices whenever
one can exhibit enough low-load integral row relations.

For the complete \(k\)-subset matrix, the obstruction includes the zero-row
relations just introduced. If \(f=A^{\mathsf T}y\), then (6.15) gives the
exact inverse identity

\[
r y_i=\sum_{|S|=k}
  (2\mathbf1_{\{i\in S\}}-1)f(S).
\]

Therefore every nonzero integral \(y\) satisfies

\[
\frac1L\sum_{|S|=k}|f(S)|
\geq\frac rL
=\frac{m+1}{2m}>\frac12.
\tag{6.25}
\]

Equality is attained by a single row and also by an ordinary row pair. For
the near-equal prime weights used in (6.18), division by \(1+\eta\)
preserves the strict inequality when \(\eta<1/m\). Thus every relation has
weighted load greater than \(W/2\), where
\(w_a=\log|H_a|\) and \(W=\sum_a w_a=\log R\). For a collection of \(h\)
relations,

\[
\frac1W\sum_a w_a\sum_{\nu=1}^h|f_{y^{(\nu)}}(a)|>\frac h2.
\]

Consequently any per-block load bound \(L_0\) in (6.24) satisfies
\(L_0>h/2\), the exact negation of its hypothesis \(h>2L_0\). Hence no
single row relation, nor any product of such relations, reaches a Gaussian
denominator exponent below \(1/2\).

Finally, averaging the weighted loads of \(e_i-e_j\) over all row pairs
shows that some pair has load at most

\[
\frac{k(m-k)}{\binom m2}W
=\frac{m+1}{2m}W.
\]

Even the best direct consequence is therefore only

\[
2R^{-(m+1)/(2m)}\leq C R^{-1/2},
\qquad\text{or}\qquad
\frac2C\leq R^{1/(2m)},
\tag{6.26}
\]

which becomes automatic as \(R\) grows.

Exact phase equality remains impossible by Gaussian unique factorization,
because \(A\) has full row rank and the blocks are mutually
conjugate-coprime.  What is now isolated is strictly simultaneous:
endpoint-scale approximation to the high-dimensional hypersimplex phase
kernel.  Coordinate-torsion separation, arbitrary products of integral row
relations, Plotkin stability, and defect-only arguments all fail on this
same genuine family.  A positive proof needs an arithmetic estimate coupling
the kernel trades, with constants independent of
\(L=\binom m{(m+1)/2}\) and of the conductor support.

### 6.7 Fixed-degree polynomial templates cannot realize the obstruction

There is nevertheless a rigorous obstruction stronger than the fixed-shift
calculation in Section 7. It rules out every fixed polynomial block
template, even when its block degrees are unequal.

Write the columns by their complements

\[
\mathcal T=\{T\subset[m]:|T|=s\},
\]

and, for every \(T\in\mathcal T\), let
\(H_T(X)\in\mathbf Z[i][X]\) be nonconstant, with positive real leading
coefficient and

\[
\gcd_{\mathbf Q(i)[X]}(H_T,\overline H_U)=1
\qquad(T,U\in\mathcal T).
\tag{6.27}
\]

Put

\[
G=\prod_TH_T,
\qquad
A_i=\prod_{i\notin T}H_T,
\qquad
B_i=\prod_{i\in T}H_T=G/A_i,
\tag{6.28}
\]

and set \(Q=\deg G=\sum_Tq_T\), where \(q_T=\deg H_T\). At a positive
integer parameter \(t\), the corresponding common-circle points are

\[
z_0(t)=\overline{G(t)},
\qquad
z_i(t)=A_i(t)\overline{B_i(t)},
\qquad
R(t)=|G(t)|.
\tag{6.29}
\]

Suppose that, on an unbounded sequence of \(t\), these points lie on an arc
of length at most \(C R(t)^\alpha\). Since every
\(H_T(t)/\overline{H_T(t)}\) tends to one, there is an analytic logarithm
\(\ell_T(t)\to0\). If

\[
a_i=\sum_{i\notin T}\ell_T,
\qquad
F=\sum_T\ell_T,
\]

the arc condition and \(R(t)\asymp t^Q\) give

\[
a_i=O_C\bigl(t^{-Q(1-\alpha)}\bigr).
\tag{6.30}
\]

Every \(s\)-set \(T\) is omitted by \(k=s+1\) indices, so the exact identity

\[
\sum_i a_i=kF
\]

gives \(F=O_C(t^{-Q(1-\alpha)})\). Hence the complementary logarithm

\[
\log\frac{B_i(t)}{\overline{B_i(t)}}=F-a_i
\]

has the same bound. Let

\[
D_i=\deg B_i=\sum_{i\in T}q_T.
\]

Double counting gives \(\sum_iD_i=sQ\), so some \(i_0\) satisfies

\[
D_{i_0}\leq\frac{s}{m}Q.
\tag{6.31}
\]

If

\[
\alpha<1-\frac{s}{m}=\frac{m+1}{2m},
\tag{6.32}
\]

then

\[
B_{i_0}(t)-\overline{B_{i_0}(t)}
=O_C\bigl(t^{D_{i_0}-Q(1-\alpha)}\bigr)=o(1).
\tag{6.33}
\]

The left side is a fixed Gaussian polynomial evaluated on an unbounded real
sequence, so (6.33) forces
\(B_{i_0}=\overline{B_{i_0}}\) identically. But (6.27) makes
\(B_{i_0}\) coprime to its conjugate, forcing it to be constant, contrary
to (6.28). Thus no such polynomial template exists. In particular, the
endpoint \(\alpha=1/2\) has the strict degree margin \(Q/(2m)\). The same
argument covers multiscale Prouhet templates and substitutions such as
\(t\mapsto a^n\) for fixed \(a>1\): the degrees \(q_T\) were arbitrary.

Without (6.27), the proof identifies the exact failure instead: one
low-degree complementary star is conjugation-invariant, so
\(\operatorname{Res}(H_T,\overline H_U)=0\) for some \(T,U\). A successful
fixed template must therefore reuse a nonconstant conjugate factor and no
longer represents independent oriented conductor columns. What the proof
does not cover is a sequence whose blocks or algebraic degrees change with
the radius.

### 6.8 The surviving value-level system is a small-cofactor Plücker system

The arbitrary value-level problem can be made equally explicit. Let the
complete \(k\)-subset columns carry pairwise coprime nonunit blocks
\(H_S\in\mathbf Z[i]\), and assume
\(\gcd(H_S,\overline H_U)=1\) for every \(S,U\), including \(S=U\). Put

\[
R=\prod_S|H_S|,
\qquad
A_i=\prod_{S\ni i}H_S.
\tag{6.34}
\]

Assume that \(1,A_1/\overline A_1,\ldots,A_m/\overline A_m\) lie in an
angular interval of width
\(\Delta\leq C R^{-1/2}\). Multiplication of each \(A_i\) by a sign lets
us write

\[
A_i=x_i+iy_i,
\qquad
|y_i|\leq\frac{\Delta}{2}|A_i|.
\tag{6.35}
\]

For distinct \(i,j\), define

\[
K_{ij}=\prod_{S\supset\{i,j\}}H_S,
\qquad
P_{ij}=N(K_{ij}),
\qquad
D_{ij}=x_i y_j-x_j y_i.
\tag{6.36}
\]

The factor \(\overline K_{ij}K_{ij}=P_{ij}\) divides
\(\overline A_iA_j\), and hence divides its imaginary part \(D_{ij}\).
The coprimality hypotheses and the fact that two incidence rows differ show
that \(D_{ij}\ne0\). Therefore

\[
t_{ij}:=\frac{D_{ij}}{P_{ij}}\in\mathbf Z\setminus\{0\},
\qquad
|t_{ij}|\leq
\frac C2\frac{|A_iA_j|}{P_{ij}}R^{-1/2}.
\tag{6.37}
\]

For three distinct indices, put

\[
P_{ijk}=N\!\left(\prod_{S\supset\{i,j,k\}}H_S\right),
\qquad
E_{jk;i}=P_{jk}/P_{ijk}.
\]

The three integers \(E_{jk;i},E_{ik;j},E_{ij;k}\) use the mutually
disjoint block profiles that contain exactly the displayed pair among
\(\{i,j,k\}\), so they are pairwise coprime. Dividing the elementary
minor identity by \(P_{ijk}\) gives

\[
y_iE_{jk;i}t_{jk}-y_jE_{ik;j}t_{ik}
  +y_kE_{ij;k}t_{ij}=0.
\tag{6.38}
\]

There is also a coefficient-only four-index form. For
\(I=\{i,j,k,l\}\), let

\[
Q_{ij\mid kl}=
\prod_{S:\,S\cap I\in\{\{i,j\},\{k,l\}\}}N(H_S).
\]

After dividing the ordinary Plücker relation by the common factors from
the three- and four-element intersection profiles, one obtains

\[
Q_{ij\mid kl}t_{ij}t_{kl}
-Q_{ik\mid jl}t_{ik}t_{jl}
+Q_{il\mid jk}t_{il}t_{jk}=0.
\tag{6.39}
\]

The three \(Q\)'s in (6.39) are pairwise coprime: every remaining block has
an intersection of size two with \(I\), and the three perfect matchings
partition those six profiles.

The signs carry additional information. The full \(M=m+1=2k\) fair code
contains one representative of every complementary pair of \(k\)-subsets
of its \(M\) rows. We may therefore choose the geometric arc-end point as
the zero row and reorient every cut away from it; the remaining incidence
matrix is again the complete \(k\)-subset matrix on \(m\) rows. Order those
rows by increasing phase. The signs in (6.35) can then be chosen so that
\(y_i>0\) and \(t_{ij}>0\) for \(i<j\). Consequently (6.39) is the
totally positive relation

\[
Q_{ik\mid jl}t_{ik}t_{jl}
=Q_{ij\mid kl}t_{ij}t_{kl}
+Q_{il\mid jk}t_{il}t_{jk}
\qquad(i<j<k<l).
\tag{6.39a}
\]

For \(m=5\), indexing a block by its two-element complement gives
\(Q_{ij\mid kl}=p_{ij}p_{kl}\), where
\(p_{ij}=N(H_{[5]\setminus\{i,j\}})\). Hence
\(w_{ij}=p_{ij}t_{ij}\) is a positive integral Plücker vector on
\(\operatorname{Gr}(2,5)\). This is a useful finite-dimensional
reformulation, but positivity alone does not contradict the size bounds
below.

In fact, the three-index equations do not provide independent extra rank
conditions in this case. Put

\[
c_i=\prod_{a\ne i}p_{\min\{a,i\},\max\{a,i\}},
\qquad z_i=c_i y_i.
\tag{6.39b}
\]

For \(i<j<k\), multiply (6.38) by
\(p_{ij}p_{ik}p_{jk}\). Using
\(E_{jk;i}=p_{il}p_{im}\), where \(\{l,m\}\) is the complement of
\(\{i,j,k\}\), gives exactly

\[
z_iw_{jk}-z_jw_{ik}+z_kw_{ij}=0.
\tag{6.39c}
\]

Because \(w\) is a nonzero decomposable bivector, the ten equations
(6.39c) say precisely that \(z\) lies in the same rank-two plane represented
by \(w\). Thus the surviving arithmetic information is not additional
Grassmannian rank: it is the simultaneous divisibility
\(p_{ij}\mid w_{ij}\), \(c_i\mid z_i\), together with the small quotient
bounds on \(t_{ij}\) and \(y_i\).

The Gaussian origin does impose a further exact *metric lift* that an
arbitrary positive Plücker vector need not possess. Put

\[
\Pi=\prod_{a<b}p_{ab}=R^2.
\]

For \(m=5\), direct partition of the ten edges of \(K_5\) gives

\[
N(A_i)=\frac{\Pi}{c_i},
\qquad
P_{ij}=\frac{p_{ij}\Pi}{c_ic_j}.
\tag{6.39d}
\]

The common Gaussian factor \(K_{ij}\) shows that \(P_{ij}\) divides both
the real and imaginary parts of \(A_i\overline{A_j}\). Hence

\[
u_{ij}:=\frac{x_ix_j+y_iy_j}{P_{ij}}\in\mathbf Z,
\qquad
b_{ij}:=u_{ij}-\mathbf i t_{ij}
=\frac{A_i\overline{A_j}}{P_{ij}}\in\mathbf Z[\mathbf i].
\tag{6.39e}
\]

Taking norms in (6.39e) yields ten simultaneous near-square identities

\[
u_{ij}^2+t_{ij}^2=\frac{c_ic_j}{p_{ij}^2}.
\tag{6.39f}
\]

They can also be expressed as an integral rank-two Gram constraint. Define

\[
X_i=c_ix_i,
\qquad
s_{ij}=p_{ij}u_{ij}.
\]

Then (6.39d)--(6.39e) give

\[
X_i^2+z_i^2=\Pi c_i,
\qquad
X_iz_j-X_jz_i=\Pi w_{ij},
\qquad
X_iX_j+z_iz_j=\Pi s_{ij}.
\tag{6.39g}
\]

Thus the symmetric integer matrix with diagonal \(c_i\) and off-diagonal
entries \(s_{ij}\) is positive semidefinite of rank two, and

\[
c_ic_j-s_{ij}^2=w_{ij}^2.
\tag{6.39h}
\]

The ten Gaussian numbers in (6.39e) also obey an exact cocycle. If
\(i,j,k\) are distinct, \(\{r,s\}=[5]\setminus\{i,j,k\}\), and
\(b_{ji}=\overline{b_{ij}}\), then

\[
b_{ij}b_{jk}=p_{ik}p_{jr}p_{js}b_{ik}.
\tag{6.39i}
\]

Indeed, before substituting (6.39d), the real coefficient on the right is
\(N(A_j)P_{ik}/(P_{ij}P_{jk})\). Equations (6.39f)--(6.39i) are necessary
Gaussian reconstruction data beyond the abstract plane-incidence statement
(6.39c).

They still stop at a sharp individual gap barrier. If all \(p_{ab}\) have
common scale \(P\), put \(d_{ij}=c_ic_j/p_{ij}^2\asymp P^6\). For large
\(R\), the short arc makes \(u_{ij}>0\), and (6.40a) below gives

\[
0<\sqrt{d_{ij}}-u_{ij}
=\frac{t_{ij}^2}{\sqrt{d_{ij}}+u_{ij}}
\ll_C P^{-2}.
\tag{6.39j}
\]

If \(d_{ij}\) were a square, (6.39j) would be impossible for large \(P\).
If it is nonsquare, the upper bound in (6.39j) eventually makes
\(u_{ij}=\lfloor\sqrt{d_{ij}}\rfloor\), but the elementary lower bound is
only \((2\sqrt{d_{ij}})^{-1}\asymp P^{-3}\). Thus the individual
near-square estimate is short by one factor of \(P\); any gain must use the
ten identities and the private prime supports simultaneously.

The metric lift can be sharpened from real rank two to Hermitian rank one.
Let

\[
B_i=G/A_i,
\qquad N(B_i)=c_i.
\]

Then (6.39d)--(6.39e) give

\[
p_{ij}b_{ij}=\overline{B_i}B_j.
\tag{6.39k}
\]

Consequently the Hermitian matrix with diagonal entries \(c_i\) and
off-diagonal entries \(p_{ij}b_{ij}\) is exactly
\((\overline{B_i}B_j)_{i,j}\), hence has rank one.  Equivalently, if
\(v_i=c_iA_i=G\overline{B_i}\), that matrix is
\(\Pi^{-1}(v_i\overline{v_j})_{i,j}\).  Thus the cocycles in (6.39i) are
Gaussian reconstruction identities, not independent inequalities; an
abstract positive Plücker vector need not possess this lift.

There is nevertheless a new phase restriction.  For the original products,
before choosing representatives of their arguments modulo \(\pi\),

\[
\prod_{i=1}^5A_i=G^3,
\qquad
\prod_{i=1}^5B_i=G^2,
\qquad
\prod_{i=1}^5|B_i|=R^2.
\tag{6.39l}
\]

Use the geometric endpoint as phase zero and choose the signs in (6.35) so
that

\[
\theta_i=\arg A_i\in[0,\Delta/2],
\qquad
\mu=\frac13\sum_i\theta_i.
\]

Taking arguments in (6.39l), modulo the signs, gives one common
\(n\in\{0,1,2\}\) such that

\[
\arg B_i\equiv\frac{n\pi}{3}+\varepsilon_i\pmod\pi,
\qquad
\varepsilon_i=\mu-\theta_i,
\qquad
|\varepsilon_i|\leq\frac{2\Delta}{3}.
\tag{6.39m}
\]

The branch \(n=0\) is impossible once the radius is large.  Indeed,
(6.39l) supplies an index with \(|B_i|\leq R^{2/5}\), and then

\[
|\operatorname{Im}B_i|
 \leq |B_i|\,|\varepsilon_i|
 \leq\frac{2C}{3}R^{-1/10}.
\tag{6.39n}
\]

If \(R>(2C/3)^{10}\), the integral imaginary part must vanish.  This makes
the nonunit \(B_i\) equal to its conjugate, contradicting
\(\gcd(B_i,\overline{B_i})=1\).  Every sufficiently large survivor is
therefore in one of the two nontrivial cubic branches \(n=1,2\).

Those two branches have an exact Pell-type residue.  Write
\(B_i=a_i+\mathbf i b_i\), and define

\[
F_i=3a_i^2-b_i^2
   =B_i^2+B_i\overline{B_i}+\overline{B_i}^{\,2},
\qquad c_i=a_i^2+b_i^2.
\tag{6.39o}
\]

The function \(3\cos^2\phi-\sin^2\phi\) vanishes at both nontrivial branch
angles and has derivative of absolute value at most \(4\).  Equations
(6.39m) and \(\Delta\leq CR^{-1/2}\) therefore give

\[
0<|F_i|\leq\frac{8}{3}\Delta c_i
          \leq\frac{8C}{3}R^{-1/2}c_i.
\tag{6.39p}
\]

Moreover \(\gcd(F_i,c_i)=1\).  For if a Gaussian prime
\(\varpi\mid B_i\), conjugate-coprimality and (6.39o) give
\(F_i\equiv\overline{B_i}^{\,2}\not\equiv0\pmod\varpi\); the same argument
applies to a prime dividing \(\overline{B_i}\).  Since
\(\prod_i c_i=R^4\), at least one index has \(c_i\leq R^{4/5}\), and hence

\[
0<|F_i|\leq\frac{8C}{3}R^{3/10}.
\tag{6.39q}
\]

Finally, for

\[
L_{ij}=3a_ia_j-b_ib_j,
\]

the two-by-two determinant identity and
\(\operatorname{Im}(\overline{B_i}B_j)=-p_{ij}t_{ij}\) yield

\[
L_{ij}^2-F_iF_j=3p_{ij}^2t_{ij}^2,
\qquad
L_{ij}^2\equiv F_iF_j\pmod{p_{ij}^2}.
\tag{6.39r}
\]

Each \(F_i\) is coprime to its four incident \(p_{ij}\).  This is a
strictly stronger arithmetic target than the bare Plücker equations, but
it still lands at the same square-root barrier: the residues \(F_iF_j\) and
the pair-dependent corrections \(p_{ij}^2t_{ij}^2\) are not fixed, and
the available size bound does not make either side smaller than the private
modulus.

The scale separation in these equations is substantial. In the formal
equal-log model \(|H_S|=X\), so \(R=X^L\), equations (6.35) and (6.37)
give

\[
|y_i|,|t_{ij}|\leq\frac C2R^{1/(2m)}.
\tag{6.40}
\]

For \(m=5\), a block norm has scale \(p_{ij}=R^{1/5}\), so (6.40) is
exactly

\[
|y_i|,|t_{ij}|\leq\frac C2\sqrt{p_{ab}}
\tag{6.40a}
\]

in the equal-log model. Reducing a triple or Plücker equation modulo one of
its private primes therefore leaves products such as \(y_it_{jk}\) or
\(t_{ij}t_{kl}\) of size only \(O_C(p_{ab})\), not \(o(p_{ab})\).
Consequently size alone cannot turn those congruences into integral
equalities. This is the exact square-root boundary at which the elementary
rational-reconstruction uniqueness inequality \(2PQ<p\) ceases to give a
contradiction.

Each \(E\) in (6.38), and each \(Q\) in (6.39), has size

\[
R^{\gamma_m},
\qquad
\gamma_m=\frac{m^2-1}{4m(m-2)},\qquad
\gamma_m-\frac1m=\frac{(m-2)^2+3}{4m(m-2)}>0.
\tag{6.41}
\]

Thus the coefficients in every three-term relation are larger, by a fixed
power of \(R\), than the small multipliers
\(y_it_{jk}\) or \(t_{ij}t_{kl}\). One equation alone still gives no
contradiction: three large pairwise coprime integers can satisfy a linear
relation with small coefficients. The new target is the *simultaneous*
arithmetic of (6.38)--(6.39), with its shared minors and disjoint
conductor-prime supports. In the first surviving case \(m=5\), equations
(6.39b)--(6.40a) show that the abstract rank relations themselves collapse
to one plane-incidence condition; the prime divisibilities and their
square-root-sized quotients must do the remaining work. Standard
unit-equation bounds treat the union of those supports as a multiplicative
group of rank growing with \(L\); the bounds cited in Section 9 therefore do
not control this system uniformly.

For prime blocks there is an equivalent projective CRT interpretation. If
\(p_S=H_S\overline H_S\) and \(\iota_S\) is the image of the Gaussian
unit \(\mathbf i\) in
\(\mathbf Z[\mathbf i]/(H_S)\cong\mathbf F_{p_S}\), then, for a row
index \(a\),

\[
\bigl[x_a:y_a\bigr]=\bigl[-\iota_S:1\bigr]
\quad\text{in }\mathbf P^1(\mathbf F_{p_S})
\quad\Longleftrightarrow\quad a\in S.
\tag{6.42}
\]

Indeed, equality in (6.42) is exactly divisibility of \(A_a\) by \(H_S\);
the projective coordinate also permits the primes that divide \(y_a\) when
\(a\notin S\). Consequently the \(m\) rational slopes form a projective
CRT list in which every candidate agrees in
\(r=L(m+1)/(2m)\) positions and every pair co-agrees in
\(\lambda=r/2\) positions. Explicitly,

\[
mr=Lk,
\qquad
\binom m2\lambda=L\binom k2.
\tag{6.43}
\]

Every coordinate has the same agreement multiplicity \(k\), so the raw
pair count and the Cauchy step in the usual agreement second-moment/Johnson
calculation are equalities. This does not say that the standard CRT-list
theorem applies to the projective rational messages in (6.42); the exact
message-class mismatch is recorded in Section 9. Like the Plotkin equality in
Section 6.4, list decoding based only on these agreement moments cannot
supply the missing strict margin; an improvement must use the small
archimedean representatives in (6.35)--(6.40) simultaneously.

## 7. The standard fixed-shift construction stops at six at the endpoint

Cilleruelo--Granville, Section 12, take balanced signs and form

\[
v_\sigma(a)=\prod_{j=1}^{2m}(a+j+i\sigma_j),
\qquad \sum_j\sigma_j=0.
\tag{7.1}
\]

For \(m=2\), the six choices lie on one circle. The two extreme conjugate
points have real and imaginary parts

\[
X=a^4+10a^3+37a^2+60a+32,
\qquad
Y=4a^2+20a+26.
\]

Writing

\[
R_a=\prod_{j=1}^4\sqrt{(a+j)^2+1},
\]

their containing arc has length \(2R_a\arctan(Y/X)\), and

\[
\frac{2R_a\arctan(Y/X)}{\sqrt{R_a}}\longrightarrow8.
\tag{7.2}
\]

Thus any eventual bound satisfies \(B(C)\geq6\) for every \(C>8\).

There is also a rigidity statement for the whole fixed-shift template.
Let \(t_1,\ldots,t_n\) be distinct fixed real numbers and

\[
V_\sigma(a)=\prod_{j=1}^n(a+t_j+i\sigma_j).
\]

Suppose a family of sign vectors has pairwise differences
\(V_\sigma(a)-V_\tau(a)=O(a^{n/2})\). Comparing coefficients and using
Newton identities gives:

- for \(n=2m\), all signed moments
  \(\sum_j\sigma_jt_j^r\), \(0\leq r\leq m-2\), are common;
- for \(n=2m+1\), the same holds for \(0\leq r\leq m-1\).

For even \(n=2m\) with \(m\geq2\), the common zeroth moment implies that
the sets \(S_\sigma=\{j:\sigma_j=1\}\) have a common cardinality.
Simultaneously replace all sets by their complements if that cardinality
exceeds half the ground-set size. Then \(|S_\sigma|\leq m\). Sets of
size at most \(m-2\) are determined by their moments. Sets of size \(m-1\)
have monic root polynomials differing only in the constant term and hence
are pairwise disjoint. Sets of size \(m\) have root polynomials differing
by a linear polynomial, hence intersect pairwise in at most one point; pair
counting gives

\[
|\mathcal F|\binom m2\leq\binom{2m}{2}.
\]

Therefore an even fixed-shift family has at most six members (with the
trivial \(n=2\) case bounded by four). In the odd case the analogous
constant-term argument gives at most three members, and at most two once
\(n\geq5\). The six-point family (7.1) is consequently extremal within
this entire fixed-shift mechanism. This does not classify arbitrary
Gaussian conductors.

### 7.1 A withdrawn super-endpoint claim cannot be used

Oganesyan, arXiv:2107.09991v1, Theorem 1 claimed unbounded occupancy for
arcs of length \(R^\alpha\) when \(\alpha>1/2\). Its Lemma 3 and Corollary
3 claimed more: for

\[
\phi(N)\geq(2\sqrt2+\epsilon)\sqrt N,\qquad
\phi(N)=o(N),
\]

the paper asserted

\[
\sum_{j\geq3}j^2|D_j(\phi,N)|
\asymp\frac{\phi(N)^3}{N}\log\frac{\phi(N)^2}{N},
\qquad
\sum_{j\geq3}|D_j(\phi,N)|
\ll\frac{\phi(N)^3}{N}.
\tag{7.3}
\]

If (7.3) were valid with \(\epsilon\) fixed, substituting
\(\phi(x)=C\sqrt x\) would formally give
\(\max j\gg_\epsilon\sqrt{\log C}\) for
\(C\geq2\sqrt2+\epsilon\). It would therefore provide useful information
even at fixed endpoint constants.

It is not an available theorem. The arXiv v2 record says exactly that the
paper “is withdrawn because of a mistake in Lemma 2.” The proof of Lemma 3
says that its assertions follow immediately from the proofs of Lemmas 1
and 2, and the second estimate in (7.3) is the generalization of Lemma 2.
Thus neither (7.3), Corollary 3, nor the displayed consequence for the
endpoint counting function may be used without a new proof. This source
status is important because search indexes still display the v1 abstract
as though Theorem 1 were established.

## 8. Corvaja--Zannier reduces the positive-dimensional part to bounded incidence

Corvaja--Zannier, Proposition 2, states that, for fixed \(K,S,\epsilon\),
all but finitely many \(S\)-unit pairs satisfying

\[
-\sum_{\mu\in M_K}\log^-
\max\{|u-1|_\mu,|v-1|_\mu\}
>\epsilon\max\{h(u),h(v)\}
\tag{8.1}
\]

lie on subgroups

\[
u^p=v^q,\qquad (p,q)=1,\qquad
\max(|p|,|q|)\leq\epsilon^{-1}.
\tag{8.2}
\]

For an endpoint conductor grid, the left side of (8.1) is at least
\((1/2)\log|G|-O_C(1)\), while both coordinate heights are at most
\(\log|G|\). Taking \(\epsilon=4/9\), every sufficiently large endpoint
pair satisfies (8.1), and (8.2) has \(\max(|p|,|q|)\leq2\).

Relative to an arc-end point, write

\[
g_j=\pi_j/\overline\pi_j,\qquad
u=\zeta\prod_jg_j^{d_j},\qquad
v=\xi\prod_jg_j^{d'_j},\qquad |d_j|,|d'_j|\leq e_j.
\]

At \(\pi_j\), equation (8.2) gives \(pd_j=qd'_j\). Choose the reference
point at one end of the short arc, so every nontrivial ratio has a small
positive lifted argument. Relations with \(pq<0\) and the diagonal give no
non-root off-diagonal pair; coordinate-axis relations are the root pairs.
The only genuine non-root possibilities are

\[
u^2=v\qquad\text{or}\qquad u=v^2.
\tag{8.3}
\]

In the first case \(d'_j=2d_j\). Orient each prime pair for this ordered
pair by taking \(\rho_j=\pi_j\) when \(d_j\geq0\), and
\(\rho_j=\overline\pi_j\) otherwise, and put

\[
W=\prod_j\rho_j^{|d_j|},\qquad G_d=\prod_j\rho_j^{e_j}.
\]

Then \(u=\zeta W/\overline W\), \(W^2\mid G_d\), and
\(|G_d|=|G|\); the orientation of \(G_d\) may depend on this ordered pair.
Endpoint proximity supplies a fixed
\(\kappa_C\) with

\[
\left|\frac{\zeta W}{\overline W}-1\right|
\leq\kappa_C|G|^{-1/2}.
\]

Unless the ratio is exactly one, the numerator on the left is a nonzero
Gaussian integer of modulus at least \(\sqrt2\). Hence

\[
|W|\geq\frac{\sqrt2}{\kappa_C}|G|^{1/2},
\qquad
\left|\frac{G_d}{W^2}\right|\leq\frac{\kappa_C^2}{2}.
\tag{8.4}
\]

There are only \(O_C(1)\) possible Gaussian integers
\(D=G_d/W^2\). Also \(|W|\leq|G|^{1/2}\), so the numerator of the preceding
endpoint estimate -- one of
\(2a,2b,(1+i)(a+b),(1+i)(a-b)\), up to a unit for
\(W=a+bi\) -- is \(O_C(1)\). For fixed \(D\), the radius
\(|W|^2=|G|/|D|\) is fixed; a bounded coordinate on this circle leaves only
\(O_C(1)\) possibilities for \(W\). Units add a fixed factor and \(v=u^2\).
The swapped relation is identical. Thus all subgroup curves in (8.2), taken
together, contain only \(O_C(1)\) non-root ordered pairs from an endpoint
conductor grid. This is stronger than merely bounding their degrees.

### 8.1 The exact remaining object is the finite residual

Let \(F^{\mathrm{CZ}}_{K,S,4/9}\) denote the finite residual left by
Proposition 2. For a
cluster of \(m\) ratios, the preceding calculation gives only

\[
m(m-1)\leq \#F^{\mathrm{CZ}}_{K,S,4/9}+2m+O_C(1).
\tag{8.5}
\]

The theorem supplies no support-independent bound for the cardinality or
height of this residual. The height failure is real even inside the
common-conductor endpoint model. With

\[
T=30q+10,\quad G_T=(T+i)(T+3i),\quad
u_T=\frac{T+i}{T-i},\quad
v_T=\frac{T+3i}{T-3i},
\tag{8.6}
\]

the four Gaussian factors are pairwise coprime, so \(u_T,v_T\) are
multiplicatively independent. Yet both are \(O(|G_T|^{-1/2})\)-close to
one. For every fixed \(\epsilon<1/2\), (8.1) holds for large \(q\), while
no relation (8.2) holds. Thus (8.6) belongs to the finite residual for
\(S=S_{G_T}\), and its height tends to infinity.

The first nonuniform step in the proof of Proposition 2 is its invocation of
Ridout's theorem for “all but finitely many” points; the next step partitions
according to \(S^+\) and \(S^-\), giving up to \(2^{|S|}\) sign profiles.
The final Subspace-Theorem/Lang step is also qualitative. The narrow missing
statement on this route is therefore a uniform bound for the intersection
of \(F^{\mathrm{CZ}}_{K,S,4/9}\) with one endpoint conductor grid, or merely a uniform
bound for the clique number of its off-diagonal pair graph.

## 9. Other rank-independent slices do not yet match the conductor grid

Evertse's 2002 survey, Section 2.4, equation (2.3.8), makes the standard
weight issue explicit: the local weights depend on the point, and their
approximation uses a set of size at most \(c_9(N)^r\), where \(r\) is the
rank of the multiplicative group. For

\[
g_j=\pi_j/\overline\pi_j,
\]

valuation at \(\pi_j\) proves that the \(g_j\) are multiplicatively
independent modulo the four Gaussian units. A conductor with \(s\) split
prime factors therefore has norm-one rank \(s\), and the pair group has
rank \(2s\). Norm-one pairing halves place bookkeeping but does not bound
the rank.

Nasserden--Xiao, Theorem 1.4, gives an \(S\)-independent bound for an
\(S\)-unit equation only after fixing a cyclic Kummer extension \(L/K\), an
auxiliary inert place, and a residue class. The paragraph following their
theorem permits up to \(C_0 8^{|S|}\) possible Kummer fields before that
slice is fixed. The independent conductor exponent vectors do not select
one common \(L\), and the bound also retains dependence on the auxiliary
place.

Similarly, quantitative unit-equation bounds remain rank-dependent:
Beukers--Schlickewei, Theorem 1.1, gives \(2^{8r+8}\) for a two-term unit
equation, while Evertse--Schlickewei--Schmidt, Theorem 1.1, gives a bound
exponential in \(r+1\). No cited theorem turns the norm-one conductor box
into bounded rank.

That support dependence is not merely an artifact of known upper bounds in
the unrestricted rational setting.  Ha--Soundararajan, Theorem 1, constructs
sets \(S\) of \(s\) primes for which \(a+1=c\) has
\(\gg\exp(s^{1/4}/\log s)\) solutions supported on \(S\).  Their proof uses
varying products of prescribed numbers of primes from two disjoint
intervals, so it already has a Boolean exponent core.  It is not a
counterexample here: a common multiplier remains in the construction, and
the theorem does not say that both consecutive integers are squarefree or
that the solutions lie in the Gaussian norm-one conductor grid.

The projective CRT formulation (6.42) does not evade this obstruction.
Guruswami--Sahai--Sudan, Theorem 1 and Corollary 1, prove a weighted
Johnson bound for an arbitrary finite-alphabet code, so that combinatorial
theorem does apply to the code with alphabets
\(\mathbf P^1(\mathbf F_{p_S})\). Section 3, however, specializes their
*algorithm* to the one-dimensional CRT message set

\[
0\leq M<K=\prod_{j=1}^{k_0}p_j,
\qquad M\longmapsto(M\bmod p_j)_j.
\tag{9.1}
\]

The objects in (6.42) are instead projective rational messages
\([x_a:y_a]\). Some conductor primes can divide \(y_a\), so there is no
single affine dehomogenization to which (9.1) applies; even when all the
denominators are invertible, their least CRT representatives have no height
bound of the form \(M<K\) furnished by (6.35). Thus the published CRT
decoding algorithm does not contain the required message class. Its general
Johnson theorem still yields the following applicable, but nonuniform,
bound.

Applying that general combinatorial theorem to the projective code still
has the wrong dependence. In the equal-log prime-block model put

\[
W=\sum_S\log p_S=2\log R,
\qquad
\rho=\frac{m+1}{2m}.
\]

Each of the \(m\) candidates has received agreement \(\rho W\). If two
candidates agree projectively modulo \(p_S\), then \(p_S\mid D_{ab}\).
Equations (6.35)--(6.36) give

\[
\sum_{p_S\mid D_{ab}}\log p_S
 \leq \log|D_{ab}|
 \leq \frac{m+2}{2m}\log R+O_C(1)
 =b_mW+O_C(1),
\qquad b_m=\frac{m+2}{4m}.
\tag{9.2}
\]

For a list of \(M_0\) candidates, the weighted agreement second moment
therefore gives, after letting the common block size tend to infinity with
\(m\) fixed,

\[
M_0^2\rho^2-M_0\rho
 \leq M_0(M_0-1)b_m,
\qquad
M_0\leq\frac{\rho-b_m}{\rho^2-b_m}=m^2.
\tag{9.3}
\]

Here \(\rho-b_m=1/4\), but the decisive denominator is
\(\rho^2-b_m=1/(4m^2)\). The hypersimplex obstruction has only
\(M_0=m\) candidates, so it is fully compatible with (9.3). The exact
Johnson margin degenerates like \(m^{-2}\), and this route supplies no
uniform list bound.

Keeping both conjugate prime-ideal coordinates improves the calculation,
but only to exact equality.  Still in the formal equal-log model, restore
the appended zero row, so the complete fair code has \(M=m+1\) words.  Use
the Gaussian integer words \(z_i=A_i\overline{G/A_i}\), and for every block
retain the two residue coordinates

\[
z\longmapsto z\bmod(H_S),
\qquad
z\longmapsto z\bmod(\overline H_S),
\]

with received symbol zero and unit weight per block-coordinate.  The exact
coprimality assumptions give

\[
H_S\mid z_i\Longleftrightarrow i\in S,
\qquad
\overline H_S\mid z_i\Longleftrightarrow i\notin S.
\]

Thus every word agrees with the doubled received word in exactly one of the
two coordinates, so its relative agreement is

\[
\rho'=\frac12.
\tag{9.3a}
\]

Two words agree together in one coordinate precisely when they lie on the
same side of the corresponding balanced cut.  Equation (6.16) shows that
the fraction of such cuts is \((m-1)/(2m)\); relative to the doubled total
weight their co-agreement is therefore

\[
b'_m=\frac{m-1}{4m}.
\tag{9.3b}
\]

The same second-moment inequality now gives

\[
M_0\leq
\frac{\rho'-b'_m}{(\rho')^2-b'_m}
=m+1,
\qquad
(\rho')^2-b'_m=\frac1{4m}.
\tag{9.3c}
\]

This is sharper than (9.3), but the fair weak-flip construction itself has
exactly \(m+1\) words.  It saturates (9.3c), including both conjugate
coordinates.  Thus squarefreeness and the full norm-one symmetry remove one
factor of \(m\) from the formal equal-log list ceiling but create no strict
margin.  With unequal block weights \(w_S\), \(\rho'=1/2\) remains exact,
but the individual pair co-agreement is

\[
\frac{\sum_{S:\,a,b\ \mathrm{on\ the\ same\ side}}w_S}
     {2\sum_Sw_S};
\]

only its average over word pairs equals \(b'_m\).  Thus unequal weights do
not turn (9.3c) into a support-independent worst-case improvement.

The recent rational-code refinement closes the literal bad-denominator gap,
but not the radius or quantifier gaps. Abbondati--Guerrini--Lebreton,
Definition 4.1,
introduces a multi-precision encoding of \(f/g\) that records
the truncated valuation
\(\nu_{p_j}(g)=\min\{\operatorname{Val}_{p_j}(g),\lambda_j\}\). For
squarefree moduli \(\lambda_j=1\), so this records exactly whether
\(p_j\mid g\), and the symbol is exactly the reduction of the reduced
projective point \([f:g]\) in \(\mathbf P^1(\mathbf F_{p_j})\).
Proposition 4.14 proves the key-lattice conclusion
\(S_{R,2^d}\subset v_C\mathbf Z\), under the assumption that the actual
codeword is within distance \(d\) of the reduced received representative.
In the surrounding exact-SVP algorithm (\(\beta_{\rm SVP}=1\)), this gives
deterministic unique decoding under the sufficient condition

\[
d\leq\log_2\sqrt{\frac{N}{2FG}}.
\tag{9.4}
\]

For (6.42), the equal-log endpoint scales are

\[
N=R^2,
\qquad F\asymp R^\rho,
\qquad 1<G\ll_C R^{1/(2m)},
\]

whereas the disagreement distance from the received projective word is

\[
d=(1-\rho)\log_2N
  =\frac{m-1}{m}\log_2R.
\]

Even with the smallest admissible denominator bound, the right side of
(9.4) is at most

\[
\frac{3m-1}{4m}\log_2R+O_C(1),
\]

so the deterministic radius is missed by at least

\[
\frac{m-3}{4m}\log_2R-O_C(1),
\qquad(m\geq5),
\tag{9.5}
\]

which is positive for all sufficiently large \(R\).

The beyond-uniqueness gain of their Theorems 4.16 and 4.17 comes from
interleaving \(\ell>1\) fractions with one common denominator and a shared
error support. The \(m\) candidates in (6.42) are separate fractions whose
agreement supports differ, so the direct specialization is \(\ell=1\).
Writing \(d=d_v+d_e\), and denoting their shortest-vector approximation
factor by \(1\leq\beta_{\rm SVP}<3^\ell\), their equation (16) then
requires

\[
d_e\leq\frac12\left(
 \log_2\frac{N}{2FG}-\log_2(3\beta_{\rm SVP})-2d_v
 \right),
\]

and hence

\[
d\leq\frac12\left(
 \log_2\frac{N}{2FG}-\log_2(3\beta_{\rm SVP})
 \right),
\tag{9.5a}
\]

The same theorems also impose the separate hypothesis

\[
d_v\leq\log_2\sqrt{\frac{N}{6FG\beta_{\rm SVP}}}.
\tag{9.5b}
\]

Definition 1.1 and Proposition 4.14 print \(\log_2\), whereas Definition
4.3, equation (16), and Theorems 4.16--4.17 print \(\log\) without a
subscript.  Their use of \(\Lambda\leq2^d\) and of powers of \(2\) fixes
the consistent interpretation as base \(2\); the displays here use
\(\log_2\) throughout.

Retaining (9.5b) only narrows their applicable range. Thus (9.5a) is below
(9.4), not beyond it.

Artificially clearing denominators and treating several candidates as one
interleaved word does not restore the gain. Choose \(\ell\) of the
candidates and repeat the received symbol in all \(\ell\) components. A
prime block \(p_S\) is then correct as a vector coordinate exactly when all
\(\ell\) chosen indices lie in \(S\). In the equal-log model the joint
agreement fraction and distance are therefore

\[
a_\ell=\frac{\binom{m-\ell}{k-\ell}}{\binom mk}
=\frac{(k)_\ell}{(m)_\ell},
\qquad
d_\ell^{\rm joint}=2(1-a_\ell)\log_2R,
\tag{9.5c}
\]

Here \((q)_\ell=q(q-1)\cdots(q-\ell+1)\), and \(a_\ell=0\) when
\(\ell>k\). Any reduced common-denominator representation still has
\(F\gg_C R^\rho\): conjugate-coprimality makes each \(x_i/y_i\) reduced,
and its numerator has size \((1+o(1))R^\rho\). Even granting the
optimistic bound \(G=O(1)\), equation (16) is largest when \(d_v=0\) and
gives at most

\[
d_\ell^{\rm AGL}
\leq\frac{\ell}{\ell+1}(2-\rho)\log_2R+O_C(1)
=\frac{\ell}{\ell+1}\frac{3m-1}{2m}\log_2R+O_C(1).
\tag{9.5d}
\]

For \(\ell\geq2\),

\[
a_\ell\leq a_2=\frac{m+1}{4m},
\qquad
d_\ell^{\rm joint}\geq\frac{3m-1}{2m}\log_2R,
\]

which is strictly beyond (9.5d), already before accounting for denominator
growth, \(\beta_{\rm SVP}\), or valuation errors. Thus neither the direct
\(\ell=1\) use nor any interleaved grouping of the hypersimplex candidates
reaches the endpoint received word. Moreover, those theorems choose the
received word *uniformly at random* from specified error balls and bound a
probability of failure. They neither bound a worst-case list nor apply to
the deterministic hypersimplex received word. Thus the most directly
relevant 2026 extension still leaves (9.3) untouched.

Number-field CRT decoding does not add a hidden uniform list theorem.
Biasse--Quintin, Theorem 1, returns messages having at least
\(t\geq\sqrt{k(n+\varepsilon)}\) agreeing places.  Its weighted
Proposition 1 requires

\[
\sum_i a_i z_i\log N\mathfrak p_i
\geq
\sqrt{
 \log(2^{d_K^2}B^{d_K})
 \left(\sum_i z_i^2\log N\mathfrak p_i
       +\varepsilon z_{\max}^2\right)},
\tag{9.5e}
\]

where \(d_K=[K:\mathbf Q]\).  Corollary 1 gets the same Johnson radius only
under the additional separation hypothesis

\[
\log N\mathfrak p_{k+1}
\geq\max\{2d_Kk\log N\mathfrak p_k,,2d_K^2\}.
\tag{9.5f}
\]

The near-equal conductor primes violate (9.5f), and conjugate prime ideals
have equal norms rather than the strict norm ordering required in that
corollary.  Their Coppersmith specialization, Theorem 3, only reaches

\[
e<n-\sqrt{kn\,
  \frac{\log N\mathfrak p_n}{\log N\mathfrak p_1}}.
\tag{9.5g}
\]

It is again an affine integral-message theorem, not the projective rational
message class in (6.42).  Finally, the output list is obtained as the root
set of an auxiliary polynomial whose degree is proved after equation (7) to
be polynomial in \(\log N\), \(1/\varepsilon\), \(d_K\), and
\(\log|\Delta_K|\); this supplies no height-independent list cardinality.
Boneh, Theorem 2.1 and Corollary 2.2, give the rational-integer predecessor
of (9.5g), with pairwise-coprime affine moduli, and have the same mismatch.

Huang--Levin, Theorem 4.3, gives another useful qualitative slice. Applied
to the point \(Y=[1:\cdots:1]\subset\mathbf P^m\), its dimension parameter
is \(n=m\) and its exceptional-dimension parameter is
\(\ell=0\): for a **fixed** finite set \(S\) and \(\varepsilon>0\), the
large-proximity \(S\)-unit points are confined to a finite union of
zero-dimensional subgroup translates. In the prime-block conductor model,
the local conductor contributions give exactly

\[
h([1:u_1:\cdots:u_m])
 =\frac12\sum_S\log p_S=\log R,
\tag{9.6}
\]

while the common endpoint interval gives, for the fixed Weil function of
\(Y\),

\[
m_{Y,S}([1:u_1:\cdots:u_m])
 \geq\frac12\log R-O_{C,m}(1).
\tag{9.7}
\]

Hence any \(\varepsilon<1/2\) forces a sufficiently high endpoint point
into \(Y\cup Z\) once \((1/2-\varepsilon)\log R\) exceeds the additive
constant in Theorem 4.3. But the theorem quantifies both that finite union
\(Z\) and its \(O(1)\) term after \(S\); it gives no bound for the number,
degree, or translation height of the exceptional points, or for this height
threshold, as \(|S|\) varies. Consequently it does not bound the endpoint
conductor grid uniformly; moreover, the theorem's ambient dimension
\(n=m\) grows with the putative cluster.

A very recent toric gap theorem has a suggestive shape but a different
hypothesis. Dujella, arXiv:2609.03448v1, Theorem 1.1, bounds every positive
reduced \(D(\pm p^r)\)-tuple by \(2^{121}\), uniformly in the prime and its
exponent. Its Lemma 6.1 constructs a five-variable monomial-saturated toric
eliminant, and Proposition 7.2 turns that certificate into a uniform gap
principle with \(\ell<2^{35}\) and \(\gamma>2^{-26}\).

That theorem does not apply to (6.39h). A \(D(\pm p^r)\)-tuple has one
fixed correction and one common full modulus,

\[
a_i a_j\pm p^r=x_{ij}^2
\qquad(i\ne j).
\]

Our metric lift instead has

\[
c_ic_j-w_{ij}^2=s_{ij}^2,
\qquad w_{ij}^2=p_{ij}^2t_{ij}^2,
\]

where the ten corrections depend on the pair and the \(p_{ij}\) have
disjoint prime supports.  The Pell version (6.39r) still has the moving
coprime residues \(F_iF_j\) and the same pair-dependent corrections.
Neither Theorem 1.1 nor Lemma 6.1 supplies an eliminant for those moving
quantities. Adapting its toric-saturation method to
(6.39f)--(6.39r) is a concrete possible next step, but it would require a
new certificate that uses the Gaussian cocycle and the cubic-branch
constraints; it is not an invocation of the published \(D(n)\)-tuple
theorem.

## 10. Research frontier

The work above leaves two sharply formulated routes.

1. **Finite-residual route.** Prove that the Corvaja--Zannier residual
   \(F^{\mathrm{CZ}}_{K,S,4/9}\) has uniformly bounded intersection, or uniformly bounded
   clique number, on endpoint common-conductor grids. Equations (8.3)--(8.4)
   already dispose of every positive-dimensional subgroup direction.
2. **Phase-stability route.** Combine the weighted near-balance theorem
   (5.5) with the modular phase constraints. Sections 6.4--6.8 show that
   pure valuation stability, bounded row restriction, coordinate-torsion
   separation, and arbitrary products of one-row-space-relation estimates
   all fail on the exact hypersimplex family. Fixed polynomial templates
   cannot realize that family at the endpoint, but arbitrary value-level
   blocks leave the simultaneous positive Plücker system (6.38)--(6.41).
   At \(m=5\), its triple equations collapse to the same rank-two incidence
   condition; the Gaussian metric lift (6.39d)--(6.39r) is Hermitian rank
   one, eliminates the trivial cubic phase branch for
   \(R>(2C/3)^{10}\), and leaves two branches with ten simultaneous
   near-square/Pell defects.  Their individual gap estimates still miss by
   one private-prime factor. A positive proof must control that coupled
   prime-divisibility system, with constants uniform in
   \(\binom m{(m+1)/2}\), and remain stable under the positive defect allowed
   by (5.5).

Either result would be genuinely new. The currently published theorems and
the calculations in this repository do not yet prove the uniform endpoint
arc bound.

## Primary sources

1. J. Cilleruelo and A. Granville, “Close lattice points on circles,”
   *Canadian Journal of Mathematics* 61 (2009), 1214--1238:
   [Theorem 1.1, Section 12 Problem 2, and Conjecture 12.1](https://doi.org/10.4153/CJM-2009-057-2).
2. D. S. Ramana, “Arithmetical applications of an identity for the
   Vandermonde determinant,”
   *Acta Arithmetica* 130 (2007), 351--359:
   [Lemma 2.1 and Definition 2.1](https://doi.org/10.4064/aa130-4-4).
3. P. Corvaja and U. Zannier, “A lower bound for the height of a rational
   function at S-unit points,” *Monatshefte für Mathematik* 144 (2005),
   203--224: [Proposition 2](https://arxiv.org/html/math/0311030).
4. J.-H. Evertse and R. G. Ferretti, “A further improvement of the
   Quantitative Subspace Theorem,” *Annals of Mathematics* 177 (2013),
   513--590: [Theorems 2.1 and 3.1 and Corollary 3.2](https://doi.org/10.4007/annals.2013.177.2.4).
5. J.-H. Evertse, “Points on subvarieties of tori,” in *A Panorama in
   Number Theory* (2002), 214--230:
   [Section 2.4](https://pub.math.leidenuniv.nl/~evertsejh/00-torus.pdf).
6. A. Nasserden and A. Xiao, “Uniformity of fibres of period mappings and
   the S-unit equation”:
   [Theorem 1.4](https://arxiv.org/html/1910.14122).
7. F. Beukers and H. P. Schlickewei, “The equation \(x+y=1\) in finitely
   generated groups,” *Acta Arithmetica* 78 (1996), 189--199:
   [Theorem 1.1](https://eudml.org/doc/206709).
8. J.-H. Evertse, H. P. Schlickewei, and W. M. Schmidt, “Linear equations
   in variables which lie in a multiplicative group,” *Annals of
   Mathematics* 155 (2002), 807--836:
   [Theorem 1.1](https://doi.org/10.2307/3062133).
9. H.-Y. Lin, S. M. Moser, and P.-N. Chen, “Weak Flip Codes and their
   Optimality on the Binary Erasure Channel,” *IEEE Transactions on
   Information Theory* 64 (2018), 5191--5218:
   [Definition 16, Lemma 26, Theorem 43, Corollaries 40 and 46, Theorem 45,
   and Example 20](https://doi.org/10.1109/TIT.2018.2834924).
10. M. Lasoń and M. Michałek, “On the toric ideal of a matroid,”
    *Advances in Mathematics* 259 (2014), 1--12:
    [Theorems 2 and 3](https://doi.org/10.1016/j.aim.2014.03.004).
11. M. Plotkin, “Binary Codes with Specified Minimum Distance,”
    *IRE Transactions on Information Theory* 6 (1960), 445--450:
    [Theorem 1](https://doi.org/10.1109/TIT.1960.1057584).
12. K. Oganesyan, “Lattice points on small arcs,”
    arXiv:2107.09991v1 (2021):
    [Theorem 1, Lemma 3, and Corollary 3; withdrawn in v2 because of a
    mistake in Lemma 2](https://export.arxiv.org/api/query?id_list=2107.09991).
13. V. Guruswami, A. Sahai, and M. Sudan, “Soft-decision decoding of
    Chinese remainder codes,” *Proceedings of the 41st Annual Symposium on
    Foundations of Computer Science* (2000), 159--168:
    [Theorem 1, Corollary 1, and Section 3](https://doi.org/10.1109/SFCS.2000.892076).
14. K. Huang and A. Levin, “Greatest Common Divisors on the Complement of
    Numerically Parallel Divisors,” arXiv:2207.14432v1 (2022):
    [Theorem 4.3](https://arxiv.org/html/2207.14432#S4.Thmtheorem3).
15. M. Abbondati, E. Guerrini, and R. Lebreton, “Simultaneous Rational
    Number Codes: Decoding Beyond Half the Minimum Distance with
    Multiplicities and Bad Primes,” *Journal of Symbolic Computation* 132
    (2026), article 102481:
    [Definition 4.1, Proposition 4.14, and Theorems 4.16--4.17](https://arxiv.org/html/2504.08472).
16. A. Dujella, “Prime-power Diophantine tuples,” arXiv:2609.03448v1
    (2026): [Theorem 1.1, Lemma 6.1, and Proposition
    7.2](https://arxiv.org/html/2609.03448).
17. J.-F. Biasse and G. Quintin, “An algorithm for list decoding number
    field codes,” arXiv:1107.2321v2 (2012): [Theorems 1 and 3, Proposition
    1, and Corollary 1](https://arxiv.org/abs/1107.2321).
18. D. Boneh, “Finding smooth integers in short intervals using CRT
    decoding,” *Journal of Computer and System Sciences* 64 (2002),
    768--784: [Theorem 2.1 and Corollary
    2.2](https://doi.org/10.1006/jcss.2002.1827).
19. J. Ha and K. Soundararajan, “Many solutions to the \(S\)-unit equation
    \(a+1=c\),” arXiv:1902.07397v1 (2019):
    [Theorem 1](https://arxiv.org/abs/1902.07397).
