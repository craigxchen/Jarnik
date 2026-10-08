# One-shot multiprime packets and the cost of further refinement

Several excluded primes can be handled algebraically by an integral
packet split with a uniform factor-three width-cost bound. This is a
conditional projection lemma: the individual class's moments must
already vanish at the low harmonics coprime to the excluded primes.
It does not prove that coupled affine classes satisfy those hypotheses.

Two exact examples delimit the conclusion. A two-prime split can
necessarily increase cost, even after optimizing its integer gauge.
Refining completely into base-layer packets can have an unbounded
cost ratio. Neither example is asserted to be an endpoint family.

## 1. The one-shot projection

Fix distinct odd primes \(p_1,\ldots,p_r\), put
\(P=\prod_jp_j\), and fix a finite integer orientation family
\(a_{i,e}\) on odd cyclotomic indices, with actual coordinate widths
\(W_e\). Its cost is

\[
D=\sum_e\phi(e)W_e.
\]

Let \(h\) be a positive odd integer. Assume that every pair difference
\(b_e=a_{i,e}-a_{j,e}\) satisfies

\[
S_k=\sum_e b_e c_e(k)=0
\quad (k<h\text{ odd},\ \gcd(k,P)=1).                 \tag{1}
\]

There is a linear integral projection, followed if desired by
row-independent coordinate translations, into packets whose slopes
are multiplied by one of the \(p_j\). Its total weighted width cost
is at most \(3D\), and it preserves all pair phase coefficients at
odd harmonics below \(h\). Multiplying every cost by a common affine
slope gives the same statement in physical-frequency coordinates.

First handle every index \(e\) divisible by \(p_j^2\) for at least
one \(j\). Assign it to any one such prime. The exact identity

\[
\Psi_e(z)=\Psi_{e/p_j}(z^{p_j})
\]

copies this coordinate to a \(p_j\)-packet. Its cost is unchanged,
since \(p_j\phi(e/p_j)=\phi(e)\). Each such coordinate is assigned
only once. These indices have zero Ramanujan sum at harmonics coprime
to \(P\), so removing them does not change (1).

Write each remaining index uniquely as

\[
e=s\prod_{j\in S}p_j,
\quad \gcd(s,P)=1,\quad S\subseteq\{1,\ldots,r\}.
\]

Use zero coordinates where an index is absent. For a pair difference
define the full alternating difference

\[
u_s=\sum_S(-1)^{|S|}b_{S,s}.
\]

When \(\gcd(k,P)=1\), multiplicativity gives
\(S_k=\sum_su_s c_s(k)\). For \(\gcd(s,P)=1\), removing factors
of \(P\) from \(k\) leaves \(c_s(k)\) unchanged. Therefore (1)
implies

\[
\sum_su_s c_s(k)=0\quad\text{for every odd }k<h.        \tag{2}
\]

Replace only the empty-subset coordinate by

\[
b'_{\varnothing,s}=b_{\varnothing,s}-u_s
 =-\sum_{S\ne\varnothing}(-1)^{|S|}b_{S,s}.
\]

This kills the full alternating difference. By (2), the replacement
preserves every low phase coefficient. Apply this fixed integral map
to all rows; negative intermediate exponents cause no problem, because
pair differences are the relevant data and each final coordinate can
be translated by its minimum. The cost proof below does not charge
the temporarily modified empty coordinate.

Now split the zero-alternating-difference cube recursively. Take its
level-one slice in the \(p_1\) direction as a \(p_1\)-packet, using

\[
\Psi_s(z)\Psi_{p_1s}(z)=\Psi_s(z^{p_1}).
\]

Subtracting its lift leaves the level-zero slice minus the level-one
slice. That residual has zero full alternating difference in the
remaining primes. Continue in the fixed prime order. Equivalently,
the packet at prime \(p_j\), with residual subset
\(R\subseteq\{j+1,\ldots,r\}\), has coefficient

\[
c_{j,R,s}=\sum_{T\subseteq\{1,\ldots,j-1\}}
                 (-1)^{|T|}b_{T\cup\{j\}\cup R,s}.    \tag{3}
\]

The decomposition after killing \(u\) is an exact cyclotomic-product
identity. Thus its only change to pair phase coefficients was the
high-contact term (2).

To bound widths, use the triangle inequality in each coordinate of
(3). An original nonempty subset \(S\), whose active positions are
\(j_1<\cdots<j_m\), is charged relative to its original totient
weight by at most

\[
\sum_{a=1}^m
 \frac{p_{j_a}}{p_{j_a}-1}
 \prod_{b<a}\frac1{p_{j_b}-1}
 \le\frac32\sum_{a=1}^m2^{-(a-1)}<3.                  \tag{4}
\]

The empty-subset coordinate is not used in (3). Summing (4), and
adding the cost-preserving higher-depth packets, proves the bound.
Packets with identical final parameters may be combined by addition;
their actual width is at most the sum of the separate widths.

This is a simultaneous family bound, using fixed linear combinations
of the original rows. It is not obtained by choosing an unrelated
gauge for each pair.

## 2. A necessary two-prime cost increase

For two distinct odd primes \(p,q\), write a cube as
\((b_{00},b_{10},b_{01},b_{11})\), where the first bit is the
\(p\)-coordinate. The difference

\[
b=(1,0,0,-1)                                         \tag{5}
\]

has zero mixed difference. Its original width cost is
\(1+(p-1)(q-1)\). Any decomposition into one \(p\)-packet and
one \(q\)-packet on this cube has the form
\(b_{ab}=x_b+y_a\). Writing

\[
y_0=g,\quad x_0=1-g,\quad x_1=-g,\quad y_1=g-1,
\]

its weighted cost is

\[
C(g)=p|1-g|+p(q-1)|g|+q|g|+q(p-1)|g-1|.
\]

This is the two-point weighted absolute-value objective

\[
C(g)=(pq+p-q)|1-g|+(pq-p+q)|g|.
\]

Its minimum, even over real \(g\), is

\[
\min_g C(g)=pq-|p-q|
 =1+(p-1)(q-1)+2(\min(p,q)-1).                        \tag{6}
\]

For \(p=5,q=7\), the cost must increase from \(25\) to at least
\(33\). Multiplying (5) by two makes its base difference even and
changes these costs to \(50\) and \(66\).

For a two-row family, allowing a different gauge in each row merely
allows an arbitrary difference \(g\), already covered by (6).
Thus a row-dependent gauge cannot restore cost preservation in this
example. The lower bound concerns the indicated cube packet split;
it is not a theorem excluding every larger auxiliary construction.

## 3. Complete base-packet refinement has unbounded cost ratio

Let \(P\) be any squarefree product of odd primes and set

\[
b_P=1,\qquad b_1=-\mu(P),\qquad b_e=0\text{ otherwise}.
\]

Its moments vanish at every harmonic coprime to \(P\), and its
original cost is \(\phi(P)+1\). The exact Euler-product identity is

\[
\Psi_P(z)\Psi_1(z)^{-\mu(P)}
 =\prod_{\substack{d\mid P\\d>1}}
       \Psi_1(z^d)^{\mu(P/d)}.                        \tag{7}
\]

The exponents in a finite expansion into base packets \(\Psi_1(z^d)\)
are unique: choose the largest index where two proposed expansions
differ, or invert the cyclotomic divisor-incidence matrix. Consequently
the cost of complete base refinement in (7) is necessarily

\[
\sum_{\substack{d\mid P\\d>1}}d
 =\sigma(P)-1=\prod_{p\mid P}(p+1)-1.                 \tag{8}
\]

The ratio of (8) to \(\phi(P)+1\) is unbounded. Indeed,

\[
\prod_{p\mid P}\frac{p+1}{p-1}
 =\prod_{p\mid P}(1-p^{-2})
  \prod_{p\mid P}(1-p^{-1})^{-2}
\]

diverges as \(P\) runs through products of all odd primes up to a
growing cutoff. The first product stays bounded away from zero; the
second diverges by the elementary Euler-product bound for the harmonic
series. The added or subtracted ones in (8) and \(\phi(P)+1\) do
not change this conclusion. Scaling all exponents by two again
enforces even base difference without changing the ratio.

## 4. Scope of the remaining problem

Section 1 gives a genuine one-shot integral projection with a uniform
cost factor and low-jet preservation. It does not justify applying
that factor afresh an unbounded number of times. Section 3 shows that
one cannot simply claim a universal cost factor for arbitrary complete
base refinement. A stopping rule may behave differently, and these
examples do not refute every possible stopping argument.

There is also an earlier arithmetic issue. For many excluded primes,
vanishing of a coupled class-group coefficient at the surviving
harmonics has not been converted uniformly into the individual-class
hypothesis (1). The four-prime argument used in the divisibility-chain
theorem avoids only one prescribed prime. Neither this issue nor the
general incomparable-slope count is resolved here. No improvement to
the general circle count is claimed.

The [exact checker](check_multiprime_packet_cost.py) verifies the
discarded-moment identity and width bound on 160 five-row families,
including higher prime powers and nonzero discarded terms whose low
moments vanish. It also checks the two-prime gauge formula and finite
divisor expansions. Unbounded cost growth rests on the Euler-product
argument in Section 3, not the finite samples.
