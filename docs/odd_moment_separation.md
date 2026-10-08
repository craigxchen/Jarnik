# Odd phase moments force large Hamming distance

This note proves an exact separation lemma for positive norm-one
coefficients, allowing the differing subsets to have unequal sizes.
It extends the quadratic phase regime of
[the shared-ray classification](pell_shared_ray_classification.md).
The resulting distance estimates imply further asymptotic
coefficient-one bonuses in specified Pell degree ranges, but do
not establish a uniform additive bound for arbitrary circles.

## 1. Exact algebraic separation

Let \(K\) be a real quadratic field with nontrivial automorphism
\(\tau\). Let \(r_1,\ldots,r_m\) be distinct positive elements of
\(K\) satisfying \(\tau(r_j)=r_j^{-1}\).
Let \(s,v\in\{-1,1\}^m\) be distinct rows whose numbers of plus
signs have the same parity. Suppose, for an integer \(a\geq0\),
\[
\sum_js_jr_j^\ell=\sum_jv_jr_j^\ell
\qquad(\ell=1,3,\ldots,2a-1).
\tag{1}
\]
For \(a=0\) this list is empty. Then
\[
\boxed{d_H(s,v)\geq4a+2.}
\tag{2}
\]
Equivalently, if \(\nu\) is positive and odd and all positive
odd moments below \(\nu\) agree, then \(d_H(s,v)\geq2\nu\).

### Proof, including unequal subset sizes

Cancel the common coordinates and set
\[
A=\{j:s_j=1,\ v_j=-1\},\qquad
B=\{j:s_j=-1,\ v_j=1\}.
\]
These sets are disjoint, their total size is
\(d=|A|+|B|=d_H(s,v)\), and \(d\) is even.
Their cardinalities need not be equal. Equation (1) and its
quadratic conjugate give
\[
\sum_A r_j^\ell=\sum_B r_j^\ell,\qquad
\sum_A r_j^{-\ell}=\sum_B r_j^{-\ell}
\quad(\ell=1,3,\ldots,2a-1).
\tag{3}
\]
Form
\[
P(X)=\prod_A(1+r_jX)\prod_B(1-r_jX),\qquad
F(X)=P(X)-P(-X).
\tag{4}
\]
The positive moment equalities imply
\[
X^{2a+1}\mid F(X).
\tag{5}
\]
Indeed, in the formal power-series ring \(K[[X]]\),
\[
\log P(X)-\log P(-X)
=2\sum_{\ell\ {\rm positive\ odd}}
 \frac{\sum_A r_j^\ell-\sum_B r_j^\ell}{\ell}X^\ell.
\]
Since \(P(0)=1\), the claimed vanishing of the difference
follows by exponentiation.

Let \(p_d\ne0\) be the leading coefficient of \(P\). Its normalized
reversal is
\[
Q(X)=\frac{X^dP(1/X)}{p_d}
=\prod_A(1+X/r_j)\prod_B(1-X/r_j).
\tag{6}
\]
The inverse moment equalities imply
\(X^{2a+1}\mid Q(X)-Q(-X)\). Because \(d\) is even,
\[
Q(X)-Q(-X)=\frac{X^dF(1/X)}{p_d}.
\]
Consequently
\[
\deg F\leq d-(2a+1).
\tag{7}
\]
Finally \(F\ne0\). If \(P\) were even, every negative root
\(-1/r_j\) from \(A\) would have the opposite positive root
\(1/r_j\) from \(B\), and conversely. Positivity and distinctness
would force \(A=B\), contrary to disjointness and the rows
being different. Thus
\[
2a+1\leq\deg F\leq d-(2a+1),
\]
which proves (2).

Only evenness of \(|A|+|B|\) is used in the reversal step.
An assumption \(|A|=|B|\) is unnecessary.

## 2. Equality and sharpness information

At \(d=4a+2\), the proof gives the exact necessary identity
\[
\boxed{P(X)-P(-X)=cX^{2a+1}\qquad(c\ne0).}
\tag{8}
\]
Conversely, for a polynomial (4) of degree \(4a+2\), identity
(8) implies all the positive and inverse moment equalities
in (3). Thus equality is characterized by one central odd monomial.

Over arbitrary positive real coefficients, with both sets of
moment equalities imposed explicitly, the bound is sharp for
every \(a\). Choose distinct positive \(b_1,\ldots,b_{2a+1}\) and
perturb
\[
P_\epsilon(X)=
\prod_{j=1}^{2a+1}(1-b_j^2X^2)+\epsilon X^{2a+1}.
\tag{9}
\]
For sufficiently small nonzero real \(\epsilon\), all roots remain
real, distinct, and nonzero, with \(2a+1\) of each sign.
There are no opposite roots, by (8). Their reciprocal magnitudes
produce disjoint positive sets satisfying (3).
This argument does not place those coefficients in one quadratic
norm-one group.

There is an exact sharp example in \(\mathbf Q(\sqrt5)\) for
\(a=1\):
\[
\begin{aligned}
\alpha&=(423+29\sqrt5)/418,\\
\beta&=(1202+504\sqrt5)/418,\\
\gamma&=(1263-533\sqrt5)/418.
\end{aligned}
\tag{10}
\]
Each is positive with norm \(1\); the six elements consisting of
these three and their inverses are distinct. Moreover
\[
\alpha+\beta+\gamma
=\alpha^{-1}+\beta^{-1}+\gamma^{-1}
=76/11.
\tag{11}
\]
The complementary three-element subsets attain distance \(6\).
The exact norm checks are
\[
423^2-5\cdot29^2
=1202^2-5\cdot504^2
=1263^2-5\cdot533^2=418^2.
\]
No sharpness claim within a fixed quadratic norm-one group for
every higher \(a\) is made. Its arithmetic could still further
restrict the equality identity (8).

## 3. Moments forced by a Pell endpoint family

Use the shared-leading-ray model
\[
H_j(n)=h q_jt+\tau(h)\tau(q_j)t^{-1},\qquad
t=\lambda^n,\quad\lambda=9+4\sqrt5.
\]
Assume the quadratic phase regime
\(\tau(h/\overline h)=-h/\overline h\).
Then \(\kappa=\tau(h)/h\) is nonzero and purely imaginary.
Conjugate primitivity and disjoint norm supports exclude repeated
or opposite \(\tau(q_j)/q_j\); coordinate reversals make these
coefficients positive and distinct. The rows have one common
weight parity, preserved by these reversals.

The full phase expansion is
\[
\arg H_j
=\arg(hq_j)+
\operatorname{Im}\sum_{\ell\geq1}
 \frac{(-1)^{\ell+1}}{\ell}
 \kappa^\ell r_j^\ell t^{-2\ell}.
\tag{12}
\]
Every even term has zero imaginary part; every odd coefficient
is nonzero. With \(m\) degree-one blocks and fixed private factors,
the radius satisfies \(R_n\asymp t^m\). Angular width
\(O(R_n^{-1/2})=O(t^{-m/2})\) forces cancellation of every odd
phase moment with
\[
4\ell<m.
\tag{13}
\]
This inequality is strict: a term at \(4\ell=m\) can remain
compatible with a fixed endpoint constant.

The number of cancelled positive odd moments is
\[
a=a(m)=\left\lfloor\frac{m+3}{8}\right\rfloor.
\tag{14}
\]
Thus every two distinct endpoint rows have distance at least
\[
d_0=4a+2.
\tag{15}
\]

## 4. Eight-row defect counting

For eight rows, let \(c_j\in\{0,\ldots,8\}\) be the size of
column \(j\), and write
\[
\mathcal D=\sum_{j=1}^m(c_j-4)^2,\qquad
b=\#\{j:c_j=4\}.
\]
The exact total separation count is
\[
\sum_{s<v}d_H(s,v)
=\sum_jc_j(8-c_j)=16m-\mathcal D.
\]
Therefore
\[
\boxed{\mathcal D\leq16m-28(4a+2).}
\tag{16}
\]
Eight rows cannot exist when
\[
m<7a+\frac72,
\quad\text{equivalently}\quad m\leq7a+3.
\tag{17}
\]

The coefficient-one conductor bonus concerns
\(\mathcal D+b\), not just \(\mathcal D\).

### 4.1. Retaining inherited constant columns

With all columns retained, squared defects are at most \(16\).
Thus \(\mathcal D\leq16(m-b)\), and
\[
\mathcal D+b
\leq m+\frac{15}{16}\mathcal D
\leq16m-\frac{105}{4}(4a+2).
\tag{18}
\]
This is at most \(2m\) whenever
\[
\boxed{m\leq\frac{15}{8}(4a+2).}
\tag{19}
\]
The radius and \(m\) in (13)--(19) are those before removing
constant growing columns.

### 4.2. Removing constant columns and using the primitive radius

If a block occurs in one orientation in all eight rows, divide
all points by that common Gaussian factor. For a product \(G_n\)
of such factors, the new radius and arc satisfy
\[
R_n'=R_n/|G_n|,\qquad
\frac{\text{new arc length}}{\sqrt{R_n'}}
\leq\frac C{\sqrt{|G_n|}}\leq C
\]
for all sufficiently large indices. The normalized family
therefore still satisfies the endpoint condition.
Recompute \(m\) using only its active blocks and use
\(R_n'\asymp t^m\) in (13)--(14).
The stronger phase decay from the original radius may contain
additional information, but is unnecessary for the next bound.

Every active column has size in \(\{1,\ldots,7\}\), so an
unbalanced squared defect is at most \(9\). Consequently
\[
\mathcal D+b
\leq m+\frac89\mathcal D
\leq\frac{137m-224(4a+2)}9.
\tag{20}
\]
The primitive count is forced when
\[
\boxed{m\leq\frac{32}{17}(4a+2).}
\tag{21}
\]
Here the primitive radius \(R_n'\) is used throughout.
Removing constant columns is not asserted to preserve their
original inherited defect.

## 5. Degree ranges and the remaining critical sequence

For \(a\geq1\), the moment window is
\[
8a-3\leq m\leq8a+4.
\tag{22}
\]
The distance bound excludes eight rows below \(7a+4\).
Within this window, (19) forces the inherited count through
\(\lfloor(30a+15)/4\rfloor\), while (21) forces the primitive
count through \(\lfloor(128a+64)/17\rfloor\).

| \(a\) | Moment window | Eight rows excluded through | Inherited count through | Primitive count through |
|---:|---:|---:|---:|---:|
| 1 | 5–12 | 10 | 11 | 11 |
| 2 | 13–20 | 17 | 18 | 18 |
| 3 | 21–28 | 24 | 26 | 26 |
| 4 | 29–36 | 31 | 33 | 33 |
| 5 | 37–44 | 38 | 41 | 41 |
| 6 | 45–52 | 45 | 48 | 48 |
| 7 | 53–60 | none | 56 | 56 |
| 8 | 61–68 | none | 63 | 64 |
| 9 | 69–76 | none | 71 | 71 |
| 10 | 77–84 | none | 78 | 79 |
| 11 | 85–92 | none | 86 | 86 |
| 12 | 93–100 | none | 93 | 94 |
| 13 | 101–108 | none | 101 | 101 |
| 14 | 109–116 | none | none | 109 |
| \(a\geq15\) | \(8a-3\)–\(8a+4\) | none | none | none |

These are sufficient ranges from separation and cut counting.
The first unresolved degree is \(12\). Failure of a sufficient
inequality does not construct a binary code, a Pell family, or a
lattice-point configuration. The table identifies where this
distance argument alone is insufficient.

Asymptotically, \(a\sim m/8\), so (15) supplies distance only
about \(m/2\). The inherited estimate (19) would require distance
at least \(8m/15\), and the primitive estimate (21) requires
at least \(17m/32\).

## 6. Precise arithmetic consequence

For a fixed Pell family with disjoint block supports, every
variable block has log norm \(2\log t+O_{\rm family}(1)\).
With fixed private factors, its conductor therefore satisfies
\[
D+W_4-2W
=2(\mathcal D+b-2m)\log t+O_{\rm family}(1).
\tag{23}
\]
Thus \(\mathcal D+b<2m\) proves
\(D+W_4\leq2W\) for all sufficiently large valid indices of
that fixed family. Equality in the leading count gives only an
additive bound depending on the family. Neither assertion provides
a uniform radius threshold or additive constant over unrestricted
block coefficients and private heights.

The algebraic separation lemma itself is exact. Its application
to lattice points remains scoped to the shared-leading-ray
quadratic phase model, including its coprimality and positivity
reductions. General mixed eight-point configurations have not
been reduced to that model.

## 7. A stronger corollary for pure Pell shifts

For the subclass consisting of distinct fixed integer shifts of
one block, there is no degree ceiling. More precisely, suppose
\(q_j=q_0\lambda^{e_j}\) for distinct fixed integers \(e_j\).
Then
\[
r_j=r_0\lambda^{-2e_j},\qquad
r_0=\tau(q_0)/q_0\ne0.
\]
These coefficients admit no nontrivial relation with coefficients
in \(\{-1,0,1\}\). After dividing by \(r_0\) and selecting the
largest participating magnitude \(M\), the sum of every possible
smaller magnitude is at most
\[
M\sum_{k\geq1}\lambda^{-2k}
=\frac{M}{\lambda^2-1}<M.
\]
Thus even the first-moment equality cannot hold between two
distinct rows.

For every \(m\geq5\), an endpoint family forces that equality,
so a pure-shift family has at most one persistent endpoint row.
For \(m\leq4\), the earlier leading-phase classification applies:
eight rows require exactly four blocks and a complete parity
class, with every growing cut balanced. Consequently no pure-shift
family of any degree supplies a mixed critical eight-point
counterexample under these fixed-private-factor hypotheses.

This includes the block types in the canonical four-point
construction, which has only three growing blocks. It does not
extend to arbitrary fixed \(q_j\in K\); (10)--(11) explicitly
demonstrate a first-moment relation available for more general
norm-one coefficients.

## 8. A rational signed linear-block collision lemma

There is a separate elementary obstruction for literal linear blocks,
which does not use positivity, norm-one coefficients, or a common row
unit. Let

\[
 H_j(T)=a_j+iT,\qquad a_j\in\mathbf Q^\times,
\]

and let rows choose \(H_j\) or \(\overline{H_j}\), with an arbitrary
row unit in \(\mu_4\). Suppose two rows differ in \(d\leq 2s\) block
positions and their signed moments of orders
\(1,3,\ldots,2s-1\) agree. Put
\(b_j=d_j a_j\) on the changed positions, where
\(d_j\in\{-1,1\}\) records the direction of the flip. Since the
orders are odd,

\[
 p_\ell:=\sum_j b_j^\ell=0,
 \qquad \ell=1,3,\ldots,2s-1.                 \tag{24}
\]

Newton's identities give \(e_1=e_3=\cdots=0\) through the available
degree. If \(d\) is odd, the root polynomial has zero constant term
and is odd, so one of the nonzero \(b_j\) would have to be zero. If
\(d\) is even, including \(d=2s\), the root polynomial is even and
its roots occur in opposite pairs. For each such pair,
\(a_u=a_v\) with opposite flip directions, or \(a_u=-a_v\) with the
same flip direction. In the first case the two blocks coincide; in
the second,

\[
 H_{-a_u}(T)=-\overline{H_{a_u}(T)},
\]

so in either case the two factors in the actual row ratio cancel
identically. All changed factors cancel in pairs. The two complete row
values consequently differ only by their fixed row-unit ratio. Once
the leading unit directions have aligned—as they must for a shrinking
arc—that ratio is \(1\), and the two actual Gaussian points coincide.
Thus distinct rows in such a literal family have

\[
 d_H(s^{(r)},s^{(r')})\geq 2s+1.             \tag{25}
\]

This argument also covers coincident \(a_j\)'s; pairwise
nonassociateness is only needed if one wants to exclude coincident
blocks before applying the lemma. It does not require a global common
content hypothesis. The statement concerns single-use linear blocks;
an allowed block multiplicity must be expanded with its multiplicity
in the signed moments.

For \(n\) rows and \(m\) blocks, the Plotkin count gives

\[
 \binom n2(2s+1)
 \leq m\left\lfloor\frac{n^2}{4}\right\rfloor.       \tag{26}
\]

Indeed, a column with \(c\) plus signs separates \(c(n-c)\) pairs,
which is at most \(\lfloor n^2/4\rfloor\). For \(n=8,s=3\), (26)
requires

\[
 m\geq \frac{7\cdot28}{16}=12.25,
\]

and hence \(m\geq13\). A parity-bit refinement is sometimes stronger:
append to each binary row its Hamming-weight parity. Every new distance
is even and at least \(2s+2\), so

\[
 (m+1)\left\lfloor\frac{n^2}{4}\right\rfloor
 \geq (2s+2)\binom n2.                         \tag{27}
\]

For \(n=8,s=2\), this excludes \(m\leq9\) (it gives \(m+1\geq10.5\),
hence \(m\geq10\)); \(m=10\) is the first open length for moments
\(1,3\). For \(s=3\), (27) again gives \(m\geq13\).

The endpoint scale must be kept separate from this collision lemma. If
there are \(m=10\) blocks and the primitive radius retains degree ten,
then \(R(T)\asymp T^{10}\), so an arc of
length \(O(\sqrt R)\) has angular width \(O(T^{-5})\). The first and
third phase terms must vanish, while the fifth term is exactly at the
allowed threshold and need not vanish. Sharing the fifth moment would
give the stronger \(O(T^{-7})\) cancellation, but it is not forced by
the endpoint hypothesis.

## 9. Equal-weight degree windows and the first possible bonus failure

Here is a sufficient literal arithmetic hypothesis for this application.
Fix nonzero rationals \(a_j=u_j/v_j\), in lowest terms with \(v_j>0\),
whose absolute values are pairwise distinct. Use the integral blocks
\(q_j(t)=u_j+i v_jt\) at positive integer \(t\), and suppose every
sign column is active: each orientation occurs in at least one row.
Clearing these fixed denominators multiplies every row by one fixed
constant. Arbitrary row units are permitted.

All cross-block content losses are bounded independently of \(t\).
The identities
\[
v_kq_j-v_jq_k=v_ku_j-v_ju_k,\qquad
v_kq_j+v_j\overline{q_k}=v_ku_j+v_ju_k
\]
bound cross-orientation Gaussian gcds by fixed nonzero integers.
The self-conjugate gcd divides \(2u_j\). Moreover
\[
\gcd(N(q_j),N(q_k))\mid
|v_k^2u_j^2-v_j^2u_k^2|\ne0.
\]
The row polynomials have common gcd one over \(\mathbf Q(i)[t]\):
all roots \(\pm i a_j\) are distinct and each orientation is absent
from at least one row. A polynomial Bezout identity with cleared
denominators bounds their specialized common Gaussian gcd by a fixed
nonzero Gaussian integer. Therefore the primitive radius satisfies
\(R(t)\asymp t^m\). Each active cut has total prime-layer weight
\(2\log t+O_{\mathrm{family}}(1)\); all overlaps and nonsplit content
contribute only a bounded error. If \(\mathcal D\) denotes the actual
weighted defect, then
\[
\mathcal D+W_4-2W
=2(D+b-2m)\log t+O_{\mathrm{family}}(1).
\]
These error bounds depend on the fixed coefficients and sign template.
They are not uniform over moving coefficients. A strict negative
coefficient below proves the intrinsic bonus eventually for each fixed
template, not uniformly for all Gaussian endpoint clusters.

With this primitive degree and content normalization, the endpoint scale forces the odd
moments through \(2s-1\), where
\[
 s=s(m):=\left\lfloor\frac{m+1}{4}\right\rfloor.       \tag{28}
\]
Indeed \(2s-1<m/2\), while the next odd phase order may be critical.
The collision lemma gives minimum Hamming distance \(2s+1\), and the
parity-bit Plotkin bound gives
\[
 B_m:=16m-56s-40\geq0.                                  \tag{29}
\]
Here \(B_m\) is an upper bound for the original squared defect
\[
 D=\sum_j(c_j-4)^2,
 \qquad c_j=\#\{\text{plus signs in column }j\},
\]
because the parity extension has length \(m+1\), minimum distance
\(2s+2\), and \(D\) is no larger than its defect.

For active columns, the possible defects are \(0,1,4,9\). If
\(n_9,n_4,n_1\) count the columns with those nonzero defects, then the
equal-weight bonus quantity obeys the exact finite optimization bound
\[
 D+b-2m
 \leq -m+8n_9+3n_4,\qquad
 9n_9+4n_4+n_1\leq B_m,\quad
 n_9+n_4+n_1\leq m,                                   \tag{30}
\]
where \(b\) is the number of balanced columns. This follows columnwise
from \(q+\mathbf1_{q=0}=1,1,4,9\) for \(q=0,1,4,9\), respectively.
It is a combinatorial upper bound; realizing equality also requires the
moment equations and rational block coefficients.

The first allowed lengths and the resulting largest possible bonus
excess are
\[
\begin{array}{c|rrrrrr}
m&6&10&13&14&17&18\\
s(m)&1&2&3&3&4&4\\
B_m&0&8&0&16&8&24\\
\max(D+b-2m)&-6&-4&-13&-2&-11&1
\end{array}                                             \tag{31}
\]

Lengths \(7,8,9,11,12,15,16\) have \(B_m<0\) and are already excluded
by the parity-extended Plotkin inequality. Thus all equal-weight
single-use degree-\(m\) families through \(m=17\) have strict
\(D+b<2m\) at the level of this defect budget; \(m=18\) is the first
length at which the budget permits a positive excess. One maximizing
profile at \(m=18\) is two defect-\(9\) columns, one defect-\(4\)
column, and fifteen balanced columns, giving \(D=22\), \(b=15\), and
\(D+b-2m=1\). This profile is permitted only by the scalar budget;
Section 10 excludes its realization even as a binary code.

## 10. Padded incidence check at the first positive budget

The scalar budget at \(m=18\) can be sharpened without classifying the
eight-row code. Pad the parity extension by

\[
 h:=4s+3-m
\]

constant columns, so its length is \(L=4s+4\). Write \(n_r\) for the
number of original columns whose minority side has size \(r\),
\(1\le r\le4\), and let \(\beta\in\{0,1,2,3,4\}\) be the parity
column's deficiency from balance. Set

\[
 D_0=9n_1+4n_2+n_3,\qquad B_0=3n_1+2n_2+n_3,
\]

\[
 D_*=D_0+\beta^2+16h,\qquad
 B_*=B_0+\beta+4h,
\]

\[
 K=L-D_*/2,\qquad q=(2B_*-L)/4.
\]

The padded sign matrix has

\[
 TT^{\mathsf T}=L I-4A,
\]

where, writing \(d_{ij}\) for the padded Hamming distance,
\(A_{ij}:=(2d_{ij}-L)/4\) for \(i\ne j\), and \(A_{ii}=0\), the
matrix \(A\) is symmetric with nonnegative integral off-diagonal
entries. Indeed, all padded distances are even and at least
\(L/2=2s+2\). Orient each column so its minority set \(S_j\) has size
at most four, and put \(b_j:=4-|S_j|\). Its column sum is \(2b_j\),
while its pair-separation count is \(16-b_j^2\). Consequently, with
\(D_*:=\sum_j b_j^2\),

\[
 \sum_{i<j}d_{ij}=16L-D_*,\qquad
 K:=\sum_{i<j}A_{ij}=L-D_*/2.                  \tag{32}
\]

If \(W_i\) is the sum of the deficiencies of the minority sets
containing row \(i\), and \(B_*:=\sum_j b_j\), then
\[
 \sum_{j\ne i}d_{ij}=4L-B_*+2W_i,
 \qquad
 \deg_A(i)=W_i-\frac{2B_*-L}{4}.
\]
Thus \(q:=(2B_*-L)/4\) must be integral and

\[
 W_i=q+\deg_A(i)\le q+K.                         \tag{33}
\]

Thus the presence of a singleton minority column (deficiency \(3\))
requires \(q+K\ge3\). This is a necessary row-incidence condition,
separate from the scalar defect budget.

For \(m=18\), \(s=4,h=1,L=20\). Exact enumeration of the integer
profiles \((n_1,n_2,n_3,n_4,\beta)\), with \(K\ge0\), integral \(q\),
and positive original excess \(D_0+n_4-2m\), gives only

\[
\begin{array}{c|c|c|c|c|c}
(n_1,n_2,n_3,n_4)&\beta&D_*&K&q&q+K\\ \hline
(2,1,0,15)&0&38&1&1&2\\
(2,1,1,14)&1&40&0&2&2\\
(2,1,2,13)&0&40&0&2&2
\end{array}
\]

Every one contains singleton columns and is therefore excluded by
(33). Consequently the first positive scalar-budget profiles that
survive this incidence test occur at \(m=22\), so the \(m=18\)
positive entry in (31) is only a budget-level possibility and not a
code-level possibility. For example,
\((n_1,n_2,n_3,n_4,\beta)=(0,8,0,14,0)\) has original excess \(2\),
\(D_*=48\), \(K=0\), and \(q=4\); the incidence test alone does not
exclude it. No signed-moment realization is asserted by this profile
calculation.

For completeness, the only zero-excess profile at \(m=18\) surviving
the same integrality and \(K\ge0\) checks is
\((n_1,n_2,n_3,n_4,\beta)=(0,6,0,12,0)\). It has \(D_*=40\),
\(K=0\), and \(q=3\). Every nonconstant minority set is a pair, so
each \(W_i\) is even, whereas (33) gives \(W_i=3\); it is impossible.
Therefore the actual defect margin at \(m=18\) is already strict.
Together with the adjacent finite-window values
\[
\begin{array}{c|rrrr}
m&18&19&20&21\\ \hline
B_m&24&-16&0&16\\
\text{scalar excess bound}&1&\text{excluded}&-20&-9
\end{array}
\]
this proves strictness through degree \(21\) in the stated fixed-template model.
The first scalar profile surviving the incidence test is at \(m=22\),
as stated above. The tensor restriction of four rows of an order-twelve
Hadamard matrix with both rows of an order-two matrix realizes the
displayed degree-22 sign profile after deleting its constant column and
one balanced column. Thus the code-level boundary is nonempty. Its
odd-moment lift is separately excluded by the integer-certificate argument
in [the even-column continuation](critical_pell_moment_codes.md).

The [exact degree-window checker](check_odd_moment_degree_windows.py)
checks bounded Newton fixtures, the scalar budgets, the degree-18
profile exclusions, and the padded identities on actual sign matrices.
None of these statements establishes the uniform endpoint theorem.
