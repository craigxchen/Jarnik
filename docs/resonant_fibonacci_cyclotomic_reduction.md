# Resonant Fibonacci factors: exact layers and effective endpoint height

This note treats fixed products of the Gaussian Fibonacci factors
\[
G_d=f_{(d+1)/2}+i f_{(d-1)/2}\qquad(d\geq1\text{ odd}),
\]
where \(f_0=0,f_1=1\). It computes the primitive radius after
growing norm overlaps have been removed. A one-class specialization
turns the endpoint condition into a finite cyclotomic-polynomial
fibre problem. The resulting height inequality does not yet give a
uniform count for arbitrary affine rates. The earlier factors
\(H_r=f_{6r+1}+if_{6r}\) are \(G_{12r+1}\).

Within one proportional class, Sections 7--8 prove at most four rows
when all active layer indices are odd prime powers, and at most
\(2^{\lfloor4/c_r\rfloor}\) when every active index has at most
\(r\) distinct prime divisors. In particular, two distinct prime
divisors per index permit at most 128 rows. These are fixed-template
theorems; the uniform bound for arbitrary circle configurations remains
unproved.

## 1. Oriented Gaussian strong divisibility

The Fibonacci identity \(f_{m+1}^2+f_m^2=f_{2m+1}\) gives
\[
N(G_d)=f_d.                                                   \tag{1}
\]
Every common Gaussian divisor of \(G_d\) and \(\overline{G_d}\)
divides \(2\), because consecutive Fibonacci coordinates are
coprime. For odd \(e\mid d\), write \(e=2m+1\) and
\(d=(2t+1)e\). The Gaussian sequence
\(J_r=f_{r+1}+if_r\) obeys
\[
J_{r+e}=f_{e+1}J_r+f_e J_{r-1}.
\]
Since \(G_e=J_m\mid f_e=N(G_e)\), induction gives
\(G_e\mid J_{m+te}=G_d\). Combining this with
\(\gcd(f_d,f_e)=f_{\gcd(d,e)}\) proves the exact oriented identity
\[
\gcd_{\mathbf Z[i]}(G_d,G_e)\doteq G_{\gcd(d,e)}
                 \qquad(d,e\text{ odd}).                     \tag{2}
\]
Here \(\doteq\) means association by a Gaussian unit. In
particular, the shared growing factor has the *same* orientation.
For \(e=1\), the divisibility is immediate from \(G_1=1\);
the recurrence argument above is used only for \(e>1\).

For the concrete resonance in
[the fixed Pell-template note](pell_product_endpoint_obstruction.md#6-a-nonresonant-extension-with-different-pell-rates),
put \(d=12n+1\), \(A=H_n=G_d\), and
\(B=H_{13n+1}=G_{13d}\). Then \(B=AQ_n\) exactly, where
\[
N(Q_n)=f_{13d}/f_d,\qquad |Q_n|\asymp\lambda^{12n+1},
\quad \lambda=9+4\sqrt5.                                    \tag{3}
\]
The two equal-norm rows \(A\overline B,\overline A B\) equal
\(N(A)\overline{Q_n},N(A)Q_n\). Their primitive normalization
removes \(N(A)\), leaving radius \(\asymp\lambda^{12n}\).
Their angle difference is \(\asymp\lambda^{-2n}\), so they do
*not* form a fixed-endpoint pair. In contrast, the two rows
\(AB,\overline A B\) share the whole factor \(B\) and normalize
to \(A,\overline A\), with radius \(\asymp\lambda^n\).
Thus the raw radius exponent \(14\) has no invariant endpoint
meaning before the Gaussian gcd is removed.

## 2. Cyclotomic layers in a resonant affine class

Set \(\varphi=(1+\sqrt5)/2\). Choose algebraic integers
\(\alpha>0,\beta=i/\alpha\) with \(\alpha^2=\varphi\).
Then \(\alpha\beta=i\), \(\alpha^2+\beta^2=1\), both
\(\alpha,\beta\) are units, and the recurrence with initial
values \(G_1=1,G_3=1+i\) gives
\[
G_d=\frac{\alpha^d-\beta^d}{\alpha-\beta}
                       \qquad(d\text{ odd}).                 \tag{4}
\]

Suppose a finite set of odd affine indices has one proportional
class:
\[
d_s(n)=m_s d_0(n),\qquad d_0(n)=An+B>0,
\quad A>0\text{ even},\ B\text{ odd},\ m_s\text{ odd}.       \tag{5}
\]
The \(m_s\) are fixed. For each odd \(e\) dividing at least one
\(m_s\), define
\[
C_1(n)=G_{d_0(n)},\qquad
C_e(n)=\Phi_e(\alpha^{d_0(n)},\beta^{d_0(n)})
       \quad(e>1),                                           \tag{6}
\]
where \(\Phi_e(X,Y)\) is the homogeneous cyclotomic polynomial.
These are Gaussian integers. For \(e>1\) odd,
\(\Phi_e(X,Y)\) is symmetric of even degree \(\phi(e)\);
it is a polynomial in \((X+Y)^2\) and \(XY\). At the arguments
in (6),
\[
(\alpha^{d_0}+\beta^{d_0})^2
 =L_{d_0}+2i^{d_0}\in\mathbf Z[i],\qquad
\alpha^{d_0}\beta^{d_0}=i^{d_0},
\]
where \(L_k=\varphi^k+(-\varphi^{-1})^k\) is the integer Lucas
number. The cyclotomic identity and (4) give the exact factorization
\[
G_{m d_0}=\prod_{e\mid m} C_e
                         \qquad(m\text{ odd}).                \tag{7}
\]
Each layer has logarithmic norm
\[
\log N(C_e)=\phi(e)d_0\log\varphi+O_e(1),
\qquad \phi(1)=1.                                           \tag{8}
\]

Different fixed layers have only bounded Gaussian gcd as
\(n\to\infty\), including \(C_1\). To see this for \(e,f>1\),
put \(z=(\alpha/\beta)^{d_0}\). Since \(\beta\) is a unit,
a prime ideal over a common Gaussian divisor of \(C_e,C_f\)
annuls both \(\Phi_e(z)\) and \(\Phi_f(z)\), hence divides the
fixed nonzero integer resultant
\(\operatorname{Res}(\Phi_e,\Phi_f)\). The resultant Bezout
identity bounds its valuation, not merely its prime support.
For \(C_1,C_e\), a common prime annuls \(z-1\) and
\(\Phi_e(z)\), so the fixed resultant is \(\Phi_e(1)\).
Every \(C_e\) divides some \(G_{m d_0}\). Taking a common
multiple \(m\) of two layer indices and using
\(\gcd(G_{m d_0},\overline{G_{m d_0}})\mid2\) also bounds
all cross-conjugate layer gcds. Distinct proportional classes
have bounded mutual norm gcd: for nonproportional affine indices
\(a_sn+b_s,a_tn+b_t\), their integer gcd divides the fixed
nonzero determinant \(|a_sb_t-a_tb_s|\), and (1)--(2) apply.

## 3. Effective primitive radius for arbitrary fixed affine rates

Assume all affine indices are positive for sufficiently large \(n\),
with positive even slopes and odd intercepts. Partition them into classes
\(c\) of proportional coefficient pairs. In class \(c\), choose
the primitive coefficient pair with \(A_c>0\)
\((A_c,B_c)\) so that
\(d_c(n)=A_cn+B_c\) is odd and every member is
\(m_s d_c(n)\) with odd fixed \(m_s\). The intercept \(B_c\)
may be negative; \(d_c(n)>0\) eventually. For each factor \(s\),
fix a multiplicity \(r_s\geq0\). Let the equal-norm rows be
\[
Z_i(n)=u_i\prod_s
G_{d_s(n)}^{a_{is}}\overline{G_{d_s(n)}}^{\,r_s-a_{is}},
\quad 0\leq a_{is}\leq r_s,                                \tag{9}
\]
where \(u_i\in\mathbf Z[i]\setminus\{0\}\) have one common
norm and are fixed as \(n\) grows. For each class and layer put
\[
M_{c,e}=\sum_{\substack{s\in c\\e\mid m_s}}r_s,\quad
A_{i,c,e}=\sum_{\substack{s\in c\\e\mid m_s}}a_{is},\quad
W_{c,e}=\max_i A_{i,c,e}-\min_i A_{i,c,e}.                   \tag{10}
\]
The explicit product
\[
D_0(n)=\prod_{c,e}
C_{c,e}^{\,\min_i A_{i,c,e}}\,
\overline{C_{c,e}}^{\,M_{c,e}-\max_i A_{i,c,e}}             \tag{11}
\]
divides every row. After dividing by it, each row has exact
squared modulus
\[
N(u_i)\prod_{c,e}N(C_{c,e})^{W_{c,e}}.
\]
The true global Gaussian gcd can exceed \(D_0\) only by a
bounded-norm factor: away from the finitely many fixed prefactor
primes and the bounded layer overlaps, some row has zero
valuation at every orientation. Therefore the primitive radius
has the rigorous height formula
\[
\boxed{\displaystyle
\log R_n=
\frac{\log\varphi}{2}
\sum_c d_c(n)\sum_e\phi(e)W_{c,e}+O(1).}                  \tag{12}
\]
This is the effective exponent after *all* growing Gaussian
overlaps, rather than the raw sum of factor rates.

For a distinct pair of endpoint rows, restrict to one parity
subsequence of \(n\). The exact angle expansion of
\(J_r=G_{2r+1}\) has \(q=-\varphi^{-2}\), so after the limiting
phases cancel, their angle difference is a convergent sum of
terms with decay exponents \(\ell a_s\), where
\(\ell\geq1\) is odd and \(a_s\) is the slope of \(d_s(n)\).
Let \(\tau_{ij}\) be the least exponent with a nonzero combined
coefficient on that subsequence. It is finite for distinct
eventual points. If
\[
\mathcal L=\sum_c A_c\sum_e\phi(e)W_{c,e},
\]
then (12) and the fixed-endpoint arc condition force
\[
\boxed{\mathcal L\leq4\tau_{ij}
\quad\text{for every distinct pair }i,j.}                 \tag{13}
\]
One must compute \(\tau_{ij}\) after collecting coincident odd
harmonics from all affine rates; the least raw rate need not survive.
Equation (13) is a reusable necessary height inequality. It does
not alone bound the number of rows.

## 4. One class becomes a reciprocal cyclotomic fibre

Assume there is one class (5), and use
\(a_{i,e}=A_{i,e}-\min_j A_{j,e}\), so
\(0\leq a_{i,e}\leq W_e\). Put
\[
z=(\beta/\alpha)^{d_0}=i^{d_0}\varphi^{-d_0},
\qquad L=\sum_e\phi(e)W_e.                               \tag{14}
\]
Since \(d_0\) is odd, \(z\) is purely imaginary. Equations
(4),(6) become
\[
C_1=\frac{\alpha^{d_0}(1-z)}{\alpha-\beta},\quad
\overline{C_1}=\frac{\alpha^{d_0}(1+z)}{\alpha+\beta},
\quad
C_e=\alpha^{d_0\phi(e)}\Phi_e(1,z),\quad
\overline{C_e}=\alpha^{d_0\phi(e)}\Phi_e(1,-z)\ (e>1).
                                                               \tag{15}
\]
The ratio \((\alpha+\beta)/(\alpha-\beta)
=(1+2i)/\sqrt5\) has even powers in \(\mathbf Q(i)\)
and odd powers outside it. Equality of the rows' limiting
arguments, together with the equal norms of the fixed \(u_i\),
therefore forces all \(A_{i,1}\) to have the same parity.
Because \(\min_i a_{i,1}=0\), every \(a_{i,1}\) and \(W_1\)
is even. The fixed prefactors cancel the differing constant
phases in (15), leaving one common nonzero multiplier times
\[
\alpha^{d_0 L}\mathcal A_i(z),\qquad
\mathcal A_i(z)=
(1-z)^{a_{i,1}}(1+z)^{W_1-a_{i,1}}
\prod_{e>1}
\Phi_e(1,z)^{a_{i,e}}\Phi_e(1,-z)^{W_e-a_{i,e}}.
                                                               \tag{16}
\]
Each \(\mathcal A_i\) is monic reciprocal in
\(\mathbf Z[z]\), with constant coefficient one and even
degree \(L\). Its roots are on the unit circle, and
\[
\mathcal A_i(z)\mathcal A_i(-z)
=\prod_e\bigl(\Psi_e(z)\Psi_e(-z)\bigr)^{W_e}
                                                               \tag{17}
\]
is independent of \(i\), where
\(\Psi_1(z)=1-z\) and \(\Psi_e(z)=\Phi_e(z)\) for
\(e>1\). The explicit common factor (11), rather than the
possibly larger true gcd, is used in (16); the latter differs
only by bounded norm and has no effect on relative angles.

For distinct rows, let
\(h_{ij}=\operatorname{ord}_{z=0}(\mathcal A_i-\mathcal A_j)\).
Equation (17) forces \(h_{ij}\) to be odd: an even first
difference would change the equal norm at that order. Thus
their angle difference is \(\asymp|z|^{h_{ij}}\).
The radius in (12) is
\(\asymp\varphi^{d_0L/2}\); therefore the endpoint condition
becomes
\[
h_{ij}\geq\left\lceil L/4\right\rceil.                    \tag{18}
\]
Precisely the coefficients of orders \(k<h_{ij}\) agree.
Writing \(D=L/2\), there are unique monic
\(P_i\in\mathbf Z[u]\) with
\(\mathcal A_i(z)=z^D P_i(z+z^{-1})\).
They are real-rooted in \([-2,2]\), have fixed
\(P_i(u)P_i(-u)\) up to the common sign \((-1)^D\),
and agree in their top coefficient positions
\(k<\lceil D/2\rceil\). Section 6 shows that every finite selection
of these layer orientations with even base-layer differences is
incidence-compatible after adding finitely many rates and fixed
equal-norm prefactors. Thus the remaining one-class problem is a
genuine cyclotomic fibre question: can its cardinality be bounded
independently of \(D\)? No such bound is proved here.

## 5. Exact odd-moment constraints and their limit

For every odd \(e\), including \(e=1\), let \(c_e(k)\) be the
Ramanujan sum. The formal logarithm about zero satisfies
\[
\log\Psi_e(z)=-\sum_{k\geq1}\frac{c_e(k)}{k}z^k.
\]
For two rows set \(b_e=a_{i,e}-a_{j,e}\).
Taking odd parts of the logarithms of (16) gives
\[
\log\mathcal A_i(z)-\log\mathcal A_j(z)
=-2\sum_{\substack{k\geq1\\k\text{ odd}}}
\frac{z^k}{k}\sum_e b_e c_e(k).                         \tag{19}
\]
The exact divisor formula for Ramanujan sums is
\[
c_e(k)=\sum_{d\mid(e,k)}d\,\mu(e/d).
\]
Consequently, if \(S_k=\sum_e b_ec_e(k)\), then
\[
S_k=\sum_{d\mid k}d B_d,\qquad
B_d=\sum_{\substack{e\\d\mid e}}b_e\mu(e/d).             \tag{20}
\]
Möbius inversion in the odd variable \(k\) shows:
if \(S_k=0\) for all odd \(k<h\), then \(B_d=0\) for all
odd \(d<h\). In particular a nontrivial pair with contact
order \(h\) must differ on some layer index \(e\geq h\).
This is a partial rigidity statement; \(\phi(e)\) can be much
smaller than \(e\), so it gives no uniform fibre bound by itself.

High contact need not bound the degree. For every odd \(r\),
\[
\mathcal A_+(z)=1+z^r+z^{2r},\qquad
\mathcal A_-(z)=1-z^r+z^{2r}                            \tag{21}
\]
are reciprocal cyclotomic products of degree \(2r\) with
contact order \(r\). In the full Fibonacci orbit they arise
from the exact identity
\[
\frac{G_{3k}}{G_k}
=\alpha^{2k}+(\alpha\beta)^k+\beta^{2k}
=L_k+i^k\qquad(k\text{ odd}).                            \tag{22}
\]
Taking \(k=r d_0(n)\), the two cross-oriented equal-norm rows
\(G_k\overline{G_{3k}},\overline{G_k}G_{3k}\) normalize
to \(L_k-i^k,L_k+i^k\) if \(3\mid k\). If \(3\nmid k\),
the Lucas number \(L_k\) is odd, and one must divide both
quotients by their extra common factor \(1+i\). Thus the
primitive chord is respectively \(2\) or \(\sqrt2\), while
the radius grows like \(\varphi^k\). The example does not make
the fibre cardinality grow. It shows why (18) alone cannot
bound degrees or rates.

## 6. Converse: arbitrary finite layer orientations lift to a template

The incidence map in (10) imposes no further restriction if one may
choose any finite collection of rates in one proportional class.
Start with any finite set of odd layer indices and orientation counts
\(A_{i,e}\in\mathbf Z\) for finitely many rows. Extend the index set
to a finite divisor-closed set \(E\) of odd positive integers, setting
unused-layer counts to zero. Möbius inversion on the divisibility
poset gives unique integral raw exponents
\[
a^{\rm raw}_{i,m}
=\sum_{\substack{e\in E\\m\mid e}}\mu(e/m)A_{i,e},
\qquad
A_{i,e}=\sum_{\substack{m\in E\\e\mid m}}a^{\rm raw}_{i,m}.
                                                               \tag{23}
\]
For each \(m\in E\), put
\(c_m=-\min_i a^{\rm raw}_{i,m}\),
\(r_m=\max_i a^{\rm raw}_{i,m}-\min_i a^{\rm raw}_{i,m}\), and
\(a_{i,m}=a^{\rm raw}_{i,m}+c_m\in[0,r_m]\). The resulting
actual layer counts are
\[
\widetilde A_{i,e}=A_{i,e}+\sum_{\substack{m\in E\\e\mid m}}c_m.
                                                               \tag{24}
\]
The offset is independent of the row. Removing the explicit common
factor (11) therefore gives exactly the prescribed normalized
orientations \(A_{i,e}-\min_jA_{j,e}\) and widths
\(W_e=\max_iA_{i,e}-\min_iA_{i,e}\). Use the fixed factors
\(G_{m d_0(n)}\), each with multiplicity \(r_m\), to realize (23).

Suppose in addition that all base counts \(A_{i,1}\) have the same
parity, as is necessary for coincident limiting arguments. Write
\(S_i=\widetilde A_{i,1}\), choose a reference row \(0\), and set
\(k_i=(S_i-S_0)/2\in\mathbf Z\). With \(F=1+2i\), choose a fixed
integer \(K\geq\max_i|k_i|\) and prefactors
\[
u_i=F^{K-k_i}\overline F^{\,K+k_i}.
                                                               \tag{25}
\]
They are Gaussian integers of the same norm \(5^{2K}\), and
\(u_i/u_0=(\overline F/F)^{k_i}\). Since the base-layer
limiting phase ratio in (15) is
\((\alpha+\beta)/(\alpha-\beta)=F/\sqrt5\) and
\((F/\sqrt5)^2=F/\overline F\), (25) makes every row's
limiting argument identical. Multiplying all \(u_i\) by a common
fixed Gaussian integer is harmless.

Thus any finite family of orientation polynomials (16) with a common
norm product (17) and even base-layer differences has an exact
one-class Fibonacci template realization. Its true global gcd differs
from the explicit layer factor only by bounded norm, by Section 2.
For that fixed template, contact orders satisfying (18) put the
distinct rows in an arc of length \(C\sqrt{R_n}\) for some finite
template-dependent \(C\) along a fixed parity subsequence. A
uniform point-count theorem at one prescribed \(C\) would require
control of these leading constants as the template varies; the
converse alone supplies no such control.

## 7. A sharp four-row bound for odd prime-power layers

Suppose every active nonbase index \(e\) in (16) is a power
of an odd prime. Multiplicities are arbitrary. Let
\(\mathcal F\) be a family of distinct orientation polynomials
with even base-layer counts and pairwise endpoint contact (18).
Then
\[
|\mathcal F|\leq4.                                         \tag{26}
\]
Families with at most one polynomial are immediate; below assume
there are at least two, so a minimum pair contact is defined.
The bound is attained: \(W_1=2,W_3=1\) gives four distinct
polynomials of degree \(L=4\), all with contact at least one.
Prime powers of different primes may occur together; no
disjoint-support assumption is needed.

For a pair of rows let \(b_e\) be its difference of layer
orientation counts, with missing layers assigned zero.
The Euler-transform coordinates (20) simplify to
\[
B_1=b_1-\sum_p b_p,\qquad
B_{p^r}=b_{p^r}-b_{p^{r+1}}\quad(r\geq1).                \tag{27}
\]
Their first nonzero odd index is exactly the contact order:
if all earlier \(B_d\) vanish, then \(S_d=dB_d\) by (20).
Let \(h\) be the minimum pair contact in \(\mathcal F\);
it is odd and \(L\leq4h\). If \(L\leq4\), only the
nonbase layers \(3,5\) are possible. The boxes
\((W_1,W_3,W_5)=(\leq4,0,0),(0,0,1),(0,2,0)\), or
\((\leq2,1,0)\) contain respectively at most \(3,2,3,4\)
even-base orientations. Thus assume \(L>4\), so \(h\geq3\).

For each odd prime \(p\), let \(q_p=p^{k_p}\) be the first
power of \(p\) with \(q_p\geq h\). Since all pair differences
in (27) vanish at \(p^r<h\), the counts on
\(p,p^2,\ldots,q_p\) differ by the same amount between any
two rows. Denote their common width by \(v_p\), and write
\(f_{i,p}=a_{i,p}\) for the prefix coordinate of row \(i\).
Every layer above \(q_p\) is a *tail*. Let \(U\) be the sum
of tail widths and \(T=\sum_pv_p\). Summing Euler totients
over each prefix gives the exact degree budget
\[
L=W_1+\sum_p(q_p-1)v_p+
  \sum_{\substack{p,a\\p^a>q_p}}\phi(p^a)W_{p^a},
\qquad \phi(p^a)\geq(p-1)q_p\geq2h
\quad\text{in the tail}.                                 \tag{28}
\]
The low coordinate \(B_1\) also vanishes for every pair, so
\[
S=a_{i,1}-\sum_p f_{i,p}
\quad\text{is independent of }i.                         \tag{29}
\]
Since \(a_{i,1}\) is even, the prefix tuple belongs to one
parity class of \(\sum_p f_{i,p}\), and (29) determines the
base count from that tuple. Thus \(U\leq2\).

If \(U=2\), the tails alone cost at least \(4h\), leaving
no base or prefix width. One tail of width two has three
orientations, and two tails of width one have four. If
\(U=1\), the tail has two orientations and leaves at most
\(2h\) for base and prefixes. For \(h\geq5\), three
prefix units would cost \(3(h-1)>2h\), so \(T\leq2\).
For \(h=3\), \(T\leq3\); equality forces all three units
onto the sole \(q_p=3\), the tail to have cost six, and
\(W_1=0\). Then (29) fixes the one prefix coordinate.
For \(T\leq2\), a fixed parity class has at most two
prefix tuples. Hence \(U=1\) also gives at most four rows.

It remains to consider \(U=0\). If all \(v_p=0\), only
the base count can vary, but two distinct base orientations
have contact one, contrary to \(h\geq3\). Otherwise let
\(q=\min\{q_p:v_p>0\}\). Equation (27) shows that
\(h=q\): all lower coordinates agree, and a pair varying
the \(q\)-prefix first differs at \(q\). The distinct
labels \(q_p\) are odd integers at least \(q\), and
\[
L=W_1+\sum_p(q_p-1)v_p\leq4q.                            \tag{30}
\]
This is the same integer box as for distinct prime layers.
If \(T\leq3\), one parity class has at most four tuples.
If \(T=4\), four distinct labels cost at least \(4q+8\),
and three labels of widths \((2,1,1)\) cost at least
\(4q+2\), so neither occurs. One label of width four
has at most three parity-compatible counts; two labels
of widths \((3,1)\) have at most four. Two labels of
widths \((2,2)\) fit (30) only if they are \(q,q+2\)
and \(W_1=0\); then (29) fixes their sum and leaves at
most three tuples.

For \(T\geq5\), (30) implies \(q\leq5\). At \(q=5\),
only \(T=5\), a sole label \(5\) of width five, and
\(W_1=0\) fit; (29) fixes its count. At \(q=3\),
\(T\leq6\). The \(T=6\) case is a sole label \(3\)
of width six and \(W_1=0\). For \(T=5\), either
only label \(3\) varies, with \(W_1\leq2\) and at
most two rows, or the labels \(3,5\) have widths
\((4,1)\), \(W_1=0\), and a fixed sum leaves at most
two tuples. This proves (26). Layers with indices
divisible by two different odd primes are not covered.

## 8. A uniform bound when the totient ratio is bounded below

The Euler transform also gives a count for broader classes of layer
indices. Let \(\mathcal F\) be a finite family in (16)--(18), with
the same layer widths \(W_e\), common norm product, and even base
counts. For the active indices, including \(e=1\) when active, put
\[
\eta=\min_{e:W_e>0}\frac{\phi(e)}e>0,
\qquad \phi(1)=1.
\]
If there are at least two distinct polynomials, then
\[
\boxed{\displaystyle
|\mathcal F|\leq2^{\lfloor4/\eta\rfloor}.}                 \tag{31}
\]
Families with at most one polynomial need no estimate. In particular,
the empty active set causes no undefined minimum in the assertion.

To prove (31), let \(h\) be the minimum pair contact. Thus
\(h\geq1\) is odd and \(L\leq4h\). Project every row to just
its orientation counts at indices \(e\geq h\). This projection
is injective. Indeed, if two rows agree there, their difference
\(b_e\) is supported at indices below \(h\). Common jets and
(20) give \(B_d=0\) for every odd \(d<h\). If the difference
were nonzero, choose its largest nonzero index \(d\). Then
\[
B_d=b_d+\sum_{\substack{e>d\\d\mid e}}\mu(e/d)b_e=b_d\ne0,
\]
a contradiction. For \(h=1\), the projection simply keeps all
coordinates and is already injective; there are no low equations
to impose. Thus the base layer and the small-degree cases require
no exception to this proof.

Let \(M=\sum_{e\geq h}W_e\) be the total width of the retained
coordinates. Their degree cost satisfies
\[
\eta h M
\leq\sum_{e\geq h}\phi(e)W_e
\leq L\leq4h,
\qquad M\leq\lfloor4/\eta\rfloor.                         \tag{32}
\]
There are at most \(\prod_{e\geq h}(W_e+1)\leq2^M\)
integer orientation tuples in this projected box. The evenness
restriction at the base layer only reduces this count. Injectivity
therefore proves (31), with arbitrary multiplicities allowed.

For a useful uniform specialization, suppose every active odd
index has at most \(r\) distinct prime divisors. If
\(p_1=3,p_2=5,p_3=7,\ldots\) are the odd primes in increasing
order, set
\[
c_r=\prod_{j=1}^r\left(1-\frac1{p_j}\right),
\qquad c_0=1.
\]
The formula \(\phi(e)/e=\prod_{p\mid e}(1-1/p)\) gives
\(\eta\geq c_r\), and hence
\[
|\mathcal F|\leq2^{\lfloor4/c_r\rfloor}.                    \tag{33}
\]
For \(r=2\), \(c_2=8/15\), so (33) gives at most \(128\)
rows. For \(r=3\), \(c_3=16/35\), it gives at most \(256\).
This includes arbitrary prime powers and mixed-prime composites;
for prime-power layers the sharper bound four in Section 7 remains
available.

For fixed \(r\), (33) is independent of the degree, the layer
indices and their multiplicities, the fixed rates, and the radius.
It is a bound within the one-class fixed-template model proved
above. In the unrestricted model \(\phi(e)/e\) can decrease as
the number of prime divisors grows, so (31) does not supply an
absolute uniform bound there. Neither conclusion is a point-count
theorem for arbitrary lattice circles or for templates varying
with the radius.

There is also an elementary degree-dependent estimate without a bound
on the number of prime divisors. If an odd \(e\geq3\) has \(s\)
distinct prime divisors, their increasing list satisfies
\(p_j\geq2j+1\). Therefore
\[
\log\frac e{\phi(e)}
\leq\sum_{j=1}^s\frac1{p_j-1}
\leq\frac{1+\log s}{2},\qquad s\leq\log_3 e,
\]
and consequently
\[
\phi(e)\geq\frac e{\sqrt{\exp(1)\log_3 e}}.
\]
The function on the right is increasing for real \(e\geq3\).
When \(h\geq3\), every retained coordinate thus costs at least
\(h/\sqrt{\exp(1)\log_3h}\) per unit of width. The same injectivity
argument gives
\[
|\mathcal F|\leq
2^{\left\lfloor4\sqrt{\exp(1)\log_3h}\right\rfloor}.
                                                               \tag{34}
\]
A nonzero difference of reciprocal degree-\(L\) polynomials has
matching lowest and highest nonzero positions \(h\) and \(L-h\),
so \(h\leq L/2\). Hence (34) is also
\(\exp(O(\sqrt{\log L}))\). If \(h=1\), then \(L\leq4\)
and the small-degree count is at most four. This improves the crude
degree-dependent count inside this model, but remains nonuniform in
the degree and gives no new bound for arbitrary circle radii.

## 9. Pairwise degree cost and a complete finite audit

For two orientation rows, cancel their polynomial gcd. Their remaining
polynomials have the form \(U(z),U(-z)\), with
\[
\deg U=d=\sum_e\phi(e)|a_{i,e}-a_{j,e}|.
\]
The base exponents remaining after cancellation are even, so both
polynomials are monic reciprocal. The canceled factor has constant
coefficient one and does not change the contact order \(h_{ij}\).
The nonzero reciprocal difference \(U(z)-U(-z)\) consequently has
lowest and highest nonzero coefficient positions \(h_{ij}\) and
\(d-h_{ij}\). Thus
\[
d\geq2h_{ij}\geq2\lceil L/4\rceil.                      \tag{35}
\]
This distance restriction alone is not a uniform count theorem.

The [complete degree-40 checker](check_resonant_cyclotomic_fibres.py)
enumerates every odd layer index with \(\phi(e)\leq40\), every
multiplicity profile of total degree at most 40, and every allowed
orientation. Completeness in the layer index uses
\(\phi(e)^2\geq e\) for odd \(e\): this holds on each odd prime
power and then follows by multiplicativity. Thus it suffices to
search \(e\leq40^2\); the 26 nonbase indices found reach 75.
The common-prefix signatures use the exact odd moments (19), retaining
all orders \(k<\lceil L/4\rceil\).

The exhaustive run contains 20,886 profiles, 1,506,841 orientation
vectors and 1,437,928 separate profile/signature pairs. The largest
fiber has four rows, and all 70,528 pairs within fibers satisfy (35).
The checker asserts these finite counts. This is evidence only for
degrees at most 40, not a proof of four rows at arbitrary degree.

One regression retains a common prefix that is not constant:
take \(W_3=W_5=W_9=1\), \(L=12\), and
\[
(a_3,a_5,a_9)=(1,0,0),(1,0,1),(0,1,0),(0,1,1).
\]
All four polynomials start with coefficients \((1,0,1)\), and their
pairwise weighted distances are six or twelve. This example prevents
replacing agreement of the jets by an assumption that every row's
lower coefficients vanish.

The checker also verifies an explicitly moving center. Take
\(W_1=2\), \(W_{11}=W_{13}=W_{17}=1\), \(L=40\), and rows
\[
(a_1,a_{11},a_{13},a_{17})
=(0,0,0,0),(2,0,1,1),(2,1,0,1),(2,1,1,0).
\]
Their coefficients at every position below ten agree, but their
common coefficient at position one is \(-1\). On the imaginary
parameter axis the common angular drift therefore starts at order
one, much earlier than the pairwise separations. Agreement of odd
jets cannot be replaced by their individual vanishing.

The [exact checker](check_resonant_fibonacci_cyclotomic_reduction.py)
tests strong Gaussian divisibility, small layer factorizations,
bounded overlap fixtures, the normalized radius formula, (16)--(20),
the high-contact quotient (21), and the converse incidence map (23).

## 10. Later rank bounds and a five-row counterexample

The general whole-fiber estimates in Section 8 have been improved in
[the affine-rank note](cyclotomic_affine_rank_reduction.md). If `rho` is
the nullity of the low Möbius matrix, the pairwise degree inequality
(35) gives a weighted affine embedding with pairwise nonpositive inner
products. Consequently every fiber has at most `2rho+1` rows, or
`rho+1` when `L<4h`. This gives `O(sqrt(log L))` without any restriction
on the indices, and 15 rows when each active index has at most two
distinct prime divisors. These improve the corresponding estimates
128 and `exp(O(sqrt(log L)))` above; the prime-power bound four remains
stronger in its scope.

The subsequent [general prime-power compression](cyclotomic_uniform_rank.md)
removes all restrictions on the active odd indices. It sends the support
to a divisor downset without increasing its charged degree or decreasing
nullity. A roots-of-unity intersection graph then gives `rho<=1125`.
Thus the **entire one-class fiber has at most 2251 rows**, independently
of its degree, prime-power indices and arbitrary multiplicities. The
squarefree proof was an intermediate step. The constants are deliberately
not optimized. Nonproportional affine classes and arbitrary circles
have not been reduced to this model.

The tempting universal bound of four suggested by the degree-40 search
is **false**. The [explicit five-row family](cyclotomic_five_row_subendpoint_family.md)
has `h=8211`, `L=32812<4h=32844`, even base differences, and an exact
Gaussian realization with normalized arc length tending to zero.
Separately, the proposed bound `rho<=3` is refuted by
[exact nullity-four supports](cyclotomic_nullity_four_counterexample.md).
These examples do not supply an unbounded family or refute the main
uniform circle goal. The finite degree-40 results remain correct in
their stated range.
