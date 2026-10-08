# Exact nested-layer extraction with \(1/M\) residual height

## Scope

This is a sharper arithmetic extraction from an actual origin-centered short-arc cluster. It retains every prime threshold layer and therefore has **no correcting factors**. Its core norms can share rational primes along nested chains. The prime-disjoint four-witness and random-sign theorems cannot be applied to these blocks without a new argument. The uniform lattice-point theorem remains unproved.

Fix \(0<C_0<\sqrt2\). From a cluster of lattice points on an arc of length at most \(C_0\sqrt R\), retain one Gaussian-unit class and divide that class by its complete Gaussian gcd. Call its size \(M\), its squared radius \(N=R^2>1\), and \(W=\log N\). This changes the original cardinality by at most a factor four and can only decrease the normalized arc length. All remaining rational primes in \(N\) are split and odd. For each choose a Gaussian prime \(\pi_p\) above \(p\), and write
\[
 z_a=\prod_p\pi_p^{a_a(p)}\bar\pi_p^{e_p-a_a(p)},\qquad
 N=\prod_p p^{e_p},\qquad 0\le a_a(p)\le e_p.
\]
The common unit has been suppressed. Put
\[
d_{ab}=\sum_p|a_a(p)-a_b(p)|\log p,\qquad
\delta_{ab}=d_{ab}-W/2,\qquad
\lambda=\log(2/C_0)>0.
\]
The primitive chord inequality gives \(\delta_{ab}\ge2\lambda\) for distinct rows. The layer-cake identity gives
\[
\sum_{a<b}d_{ab}
 =\sum_{p,t}n_{p,t}(M-n_{p,t})\log p
 \le WM^2/4,\qquad
\sum_{a<b}\delta_{ab}\le WM/4.                         \tag{1}
\]
Consequently \(W\ge4(M-1)\lambda\).

## Simultaneous balanced and low-excess selection

Fix \(2\le k\le M\), put \(m=k-1\) and \(b=2^m\), and choose an ordered \(k\)-tuple of distinct rows uniformly. On each weighted threshold layer \((p,t)\), let \(f_j(p,t)=2\mathbf1_{a_j(p)\ge t}-1\). For nonempty \(S\subseteq[k]\), let
\[
\mu_S=W^{-1}\sum_{p,t}(\log p)\prod_{j\in S}f_j(p,t).
\]
The obtuse-vector Bessel estimate from [uniform_profile_extraction.md](uniform_profile_extraction.md) gives
\[
\mathbb E A\le a:=\frac{2(2^k-1)}{M-k+1},\qquad
A=\sum_{\varnothing\ne S\subseteq[k]}\mu_S^2.
\]
For \(B=\sum_{a<b\text{ in the tuple}}\delta_{ab}\), (1) gives
\[
\mathbb E B\le b_0:=\frac{k(k-1)W}{4(M-1)}.
\]
Both \(a,b_0\) are positive. Since \(\mathbb E(A/a+B/b_0)\le2\), one tuple satisfies
\[
A\le\frac{4(2^k-1)}{M-k+1},\qquad
B\le\frac{k(k-1)W}{2(M-1)}.                             \tag{2}
\]
In particular **every** pair in it has
\[
0<\delta_{ij}\le\frac{k(k-1)W}{2(M-1)}.                \tag{3}
\]
Fourier inversion on the sign cube shows that each complete sign pattern has relative weight error at most
\[
\eta=\frac{2(2^k-1)}{\sqrt{M-k+1}}.                    \tag{4}
\]
This is simultaneous with (3), rather than a separate tuple selection.

## Exact nested core blocks

Relabel the selected rows \(0,1,\ldots,m\), with \(0\) as anchor. Define the actual primitive anchor numerator
\[
P_i=\prod_p\pi_p^{(a_i-a_0)_+}
             \bar\pi_p^{(a_0-a_i)_+}
    =x_i+iy_i,\qquad z_i/z_0=P_i/\bar P_i.              \tag{5}
\]
For a threshold \((p,t)\), set \(T(p,t)=\{i\in[m]:f_i(p,t)\ne f_0(p,t)\}\). For each nonempty \(T\subseteq[m]\), form \(H_T\) by multiplying \(\pi_p\) for every layer with \(T(p,t)=T\) and \(t>a_0(p)\), and \(\bar\pi_p\) for every layer with \(T(p,t)=T\) and \(t\le a_0(p)\). Write \(n_T=N(H_T)\) and \(w=W/b\).

Below the anchor threshold, \(T(p,t)=\{i:a_i(p)<t\}\) runs through a nested increasing chain contained in the lower-row set \(\{i:a_i<a_0\}\). Above the anchor threshold, \(T(p,t)=\{i:a_i(p)\ge t\}\) runs through a nested decreasing chain contained in the disjoint upper-row set \(\{i:a_i>a_0\}\). Thus one nonempty \(T\) cannot occur on both sides for the same \(p\); each \(H_T\) is conjugate-primitive. Distinct \(H_T\), however, can contain powers of the **same** rational prime.
For a fixed row \(i\) and prime \(p\), all its contributing layers lie on one side of \(a_0(p)\), so \(P_i\) is also conjugate-primitive.
For example, \(e_p=4\), \(a_0=2\), and nonanchor allocations
\((0,1,3,4)\) put one \(p\)-layer in each of
\(H_{\{1\}},H_{\{1,2\}},H_{\{3,4\}},H_{\{4\}}\). The first two
blocks use \(\bar\pi_p\), the latter two use \(\pi_p\), and all four
norms share \(p\).

The layer construction gives, without correcting factors,
\[
P_i=\prod_{T\ni i}H_T,\qquad
|\log n_T-w|\le\eta w
\quad(\varnothing\ne T\subseteq[m]).                  \tag{6}
\]
Indeed a fixed nonempty \(T\) merges the two opposite complete sign patterns relative to the anchor, each of weight \((1\pm\eta)W/2^k\).

More strongly, for every nonempty \(S\subseteq[m]\),
\[
N\bigl(\gcd_{\mathbb Z[i]}(P_i:i\in S)\bigr)
 =\prod_{T\supseteq S} n_T.                           \tag{7}
\]
At one \(p\), both sides count exactly the threshold layers on which all \(i\in S\) differ from the anchor in the same orientation. In particular \(N(P_i)=\prod_{T\ni i}n_T\) and, for \(i\ne j\),
\[
G_{ij}:=N(\gcd_{\mathbb Z[i]}(P_i,P_j))
       =\prod_{T\ni i,j}n_T.                           \tag{8}
\]
The Gaussian gcd of the selected **circle points** has norm \(n_\varnothing\), the product of rational norms on the empty \(T\)-layers. The selected tuple's primitive squared radius is \(N/n_\varnothing\). The formal oriented empty-layer product is generally the conjugate of that common Gaussian gcd; only its norm is used here. All displayed \(W,w,d_{ij}\) bounds refer to the pre-selection normalization \(N\).

## Exact residuals and aggregate gain

Distinct points give \(y_i\ne0\) and \(\Delta_{ij}=\operatorname{Im}(\bar P_iP_j)\ne0\). By (8),
\[
t_{ij}=\Delta_{ij}/G_{ij}\in\mathbb Z\setminus\{0\}.
\]
The same chord computation as in [endpoint_full_profile_quantifiers.md](endpoint_full_profile_quantifiers.md) now has **zero core-gcd loss**:
\[
|y_i|\le (C_0/2)e^{\delta_{0i}/2},\qquad
|t_{ij}|\le(C_0/2)e^{\delta_{ij}/2}.                   \tag{9}
\]
Because \(C_0<2\) and all displayed integers are nonzero, (2) gives the aggregate estimate
\[
\sum_{i=1}^m\log|y_i|
 +\sum_{1\le i<j\le m}\log|t_{ij}|
 \le\tfrac12 B
 \le\frac{k(k-1)W}{4(M-1)}.                            \tag{10}
\]
Each individual log residual obeys the same bound. Relative to \(w=W/2^m\), the aggregate height is \(O(k^2 2^m/M)\), while the pattern error remains \(O(2^k/\sqrt M)\). For fixed \(k\) and \(M\to\infty\), (1) forces \(w\to\infty\), and (6), (10) are a balanced full profile with subpower residuals and \(K_i=1\). The same holds in growing dimension when \(2^k/\sqrt M\to0\).

The aggregate height is independent of the selected anchor. For every
edge let \(A_{ij}\) be the primitive numerator formed from the allocation
differences \(a_i-a_j\). Exact cancellation gives
\(A_{ij}=P_i\bar P_j/G_{ij}\), up to reversing the edge, so its
imaginary part has absolute value \(|t_{ij}|\); for anchor edges it
has absolute value \(|y_i|\). Thus (10) is the sum of the logarithms
of the primitive imaginary residues on all \(\binom{k}{2}\) edges.

This gives a smaller residual scale than the disjoint-core extraction by retaining a different class of blocks; it does not prove an arithmetic height gap. Pairwise rational-prime disjointness was used in the random signed-block lemma and in coprime-norm four-witness rank arguments. Independent sign flips of the \(H_T\) can even place both Gaussian orientations above one rational prime in a product, so those arguments cannot simply be relabelled.

## Conditional growth consequences, with the new class retained

Put \(H_\Sigma=\sum_i\log|y_i|+\sum_{i<j}\log|t_{ij}|\).
Suppose that for one fixed \(m\ge4\), constants \(\varepsilon_0>0\)
and \(w_0\ge0\), and a positive function \(f\), one proves the following
**new arithmetic assertion** for the exact nested class above:

> Every profile with \(|\log n_T-w|\le\varepsilon_0 w\),
> \(w\ge w_0\), and \(H_\Sigma\le\varepsilon_0 w\) satisfies
> \(H_\Sigma\ge f(w)\).

The small-height hypothesis can be omitted if it is not needed by the
proposed theorem. The block factorization, two-chain prime structure,
conjugate-primitivity and all exact gcd/residual identities must be
allowed. A theorem only for prime-disjoint blocks cannot be substituted.
The neighborhood has fixed positive width; exact balance alone is
insufficient. No such growing lower bound is proved in this note.
By restricting the neighborhood, assume without loss that
\(0<\varepsilon_0\le1/2\). The selected block norms then exceed one
whenever \(w>0\).

Let \(k=m+1\), \(b=2^m\), and define the fixed threshold

\[
Q=\max\left\{
k,\quad k-1+\frac{4(2^k-1)^2}{\varepsilon_0^2},\quad
1+\frac{k(k-1)b}{4\varepsilon_0},\quad
1+\frac{b w_0}{4\lambda}
\right\}.                                             \tag{11}
\]

For every common-unit cluster with \(M\ge Q\), (1), (4), and (10)
give all the hypotheses of the proposed assertion. Combining its lower
bound with (10) gives

\[
f(W/b)\le\frac{k(k-1)W}{4(M-1)},\qquad
M\le\max\left\{Q,\ 1+\frac{k(k-1)W}{4f(W/b)}\right\}.  \tag{12}
\]

For \(M<Q\) the first term suffices. Formula (12) is used when
\(W/b\ge w_0\); smaller \(W\) already bounds \(M\) through (1).
To return to an arbitrary original arc, subdivide into
\(s=\max(1,\lceil C/C_0\rceil)\) pieces and retain a largest unit
class in a largest piece. Its cardinality is at least the original
cardinality divided by \(4s\). Its primitive radius is no larger than
the original radius. This transfers any eventually nondecreasing bound
in (12); an arbitrary nonmonotone expression requires its monotone
envelope over smaller radii.

For fixed dimension and fixed positive implicit constants, the resulting
conditional comparisons are:

| Eventual logarithmic-height lower bound | Consequence for the original count |
| --- | --- |
| \(f(w)=c\log w\) | \(O_C(\log R/\log\log R)\), the existing order |
| \(f(w)=c\log w\log\log w\) | \(O_C(\log R/(\log\log R\,\log\log\log R))\) |
| \(f(w)=c w^\alpha\), \(0<\alpha<1\) | \(O_C((\log R)^{1-\alpha})\) |
| \(f(w)=c w/\log w\) | \(O_C(\log\log R)\) |
| \(f(w)=c w\) | A constant independent of \(R\) |

In particular, an eventual \(f(w)/\log w\to\infty\) would improve
the general growth order. For a nonmonotone \(f\), this follows by first
minorizing \(f(w)/\log w\) by its increasing tail infimum; the resulting
upper envelope is still \(o(W/\log W)\). The old disjoint-core comparison
required a lower bound exceeding \(\sqrt{w\log w}\) by an unbounded
factor. This is a weaker prospective input for a class that permits
nested prime overlap, not an unconditional growth improvement.

As before, merely sublinear bounds do not force uniformity through
this comparison. For any family \(f_m(w)=o(w)\), choose \(M_n=n\).
At this cardinality only dimensions with \(m+1\le n\) and
\(Q_m\le n\) are available; increasing \(W_n\) does not relax
these selection requirements. Choose \(W_n\) sufficiently large to
pass their finitely many height thresholds and make
\(f_m(W_n/2^m)/(W_n/2^m)\le n^{-2}\) for every available dimension.
Every fixed \(m\) eventually becomes available, but these inequalities
are still compatible with (10) throughout that finite set. This numerical
observation constructs no Gaussian profile; it records why even the
improved extraction still needs a linear gap at some fixed dimension,
or a different arithmetic mechanism, to establish a uniform count.

## A valuation tool that survives core overlap

The positive Hankel divisor of [positive_hankel_full_profile_budget.md](positive_hankel_full_profile_budget.md) extends to these nested blocks. Let \(L>0\) be divisible by every \(|y_i|\) and \(|t_{ij}|\), put \(Q_i=Lx_i/y_i\) and \(C_i=Q_i^2+L^2=(L/y_i)^2N(P_i)\). For \(i,j\in T\), (8) and Gaussian primitivity imply, at every rational \(p\),
\[
v_p(C_i)\ge\sum_{T\ni i}v_p(n_T),\qquad
v_p(Q_i-Q_j)\ge\sum_{T\ni i,j}v_p(n_T).
\]
If \(p\mid G_{ij}\), then \(p\nmid y_i y_j\). Define the positive Hankel numerator
\[
D_r=\sum_{\substack{S\subseteq[m]\\|S|=r}}
 \left(\prod_{\substack{i<j\\i,j\in S}}(Q_i-Q_j)^2\right)
 \left(\prod_{i\notin S}C_i\right).
\]
Its Cauchy--Binet summand indexed by a size-\(r\) set \(S\) therefore has \(p\)-valuation at least
\[
\sum_T f_T(S)v_p(n_T),\qquad
f_T(S)=|T\setminus S|+|T\cap S|(|T\cap S|-1).
\]
Hence its positive determinant \(D_r\) has the guaranteed divisor
\[
\prod_p p^{c_p(r)}\mid D_r,\qquad
c_p(r)=\min_{|S|=r}\sum_T f_T(S)v_p(n_T)
\ge\sum_T\bigl(\min_{|S|=r}f_T(S)\bigr)v_p(n_T).       \tag{13}
\]
This proves the previous product-of-blocks divisor even when their norms overlap; the left-hand side can be stronger when a prime spans several threshold blocks. It supplies no universal new gap: a squarefree prime occupies one threshold layer, where (13) reduces to the old budget with its strict slack for proper minors.

The [exact checker](check_exact_nested_profile_residual_extraction.py) tests local allocation identities, nested-chain orientation, multirow Gaussian gcds, the first-moment layer identity, and overlap-aware Hankel divisibilities. These finite tests supplement the arguments above; they do not check the endpoint selection existence or establish a uniform circle bound.
