# Balanced pair products, rational directions, and a height classification

This note uses integral products of the original Gaussian factors that
are not generally integer combinations of the zero-sum valuation rows.
It gives a finite exclusion for every exactly balanced common-unit
cluster on a rational Gaussian direction branch, including complementary
rows. It also classifies the small integral monomials in the linear phase
space of the complete equally weighted half-cut model. Neither result
establishes balance or a suitable direction for arbitrary clusters.

## 1. Pair products in an exactly balanced cluster

Let M = 2k >= 4 distinct primitive circle points have the same Gaussian
unit epsilon, radius R > 1, and representation

\[
z_i=\epsilon\prod_p\pi_p^{a_i(p)}
                 \overline{\pi_p}^{e_p-a_i(p)}.
\]

Suppose every nonconstant prime threshold layer has exactly k rows.
After primitive normalization, the nested threshold sets at a fixed
prime must all coincide. Hence every a_i(p) is either zero or e_p, and
exactly k rows have the latter value. Put W = log(R^2).

For an unordered pair i,j, define

\[
\gamma_{ij}(p)=a_i(p)+a_j(p)-e_p,
\]

\[
h_{ij}=\prod_{\gamma_{ij}>0}\pi_p^{\gamma_{ij}}
        \prod_{\gamma_{ij}<0}\overline{\pi_p}^{-\gamma_{ij}},\qquad
H_{ij}=\log N(h_{ij}).
\tag{1}
\]

These are pair products rather than the pair quotients used for chord
residues. Exactly,

\[
\frac{z_i z_j}{R^2}=\epsilon^2
 \frac{h_{ij}}{\overline{h_{ij}}},\qquad
\sum_{i<j}H_{ij}=k(k-1)W.
\tag{2}
\]

For the sum, each balanced prime separates k rows from k rows. There
are 2 binomial(k,2) = k(k-1) pairs on the same side, each contributing
e_p log p to (1); the other pairs contribute zero.

The factor h_ij is a unit precisely when the two valuation rows are
complementary. Every row has at most one complementary row, because
the points are distinct and their units coincide. Therefore, if L is
the number of nonunit pair products,

\[
L\geq\binom M2-\frac M2=\frac{M(M-2)}2.
\tag{3}
\]

Every nonunit h_ij is conjugate-coprime. It lies on neither an axis nor
a diagonal: an axis would force a common nonunit conjugate factor, and
a diagonal would force divisibility by 1+i.

## 2. A finite rational-branch exclusion

Choose argument lifts phi_i in the arc, of angular width Delta. Since
each prime is balanced,

\[
\prod_i z_i=\epsilon^M R^M,
\qquad
\overline\phi:=\frac1M\sum_i\phi_i
 =\arg\epsilon+\frac{2\pi l}{M}
\tag{4}
\]

for some integer l. Equation (2) shows

\[
\arg h_{ij}\equiv\frac{2\pi l}{M}
 +\frac{\phi_i+\phi_j}{2}-\overline\phi\pmod\pi.
\]

The final error is at most (M-2) Delta/M: it is (M-2)/M times the
difference between the mean of the pair and the mean of the other rows.

Suppose

\[
\boxed{M\mid8l.}
\tag{5}
\]

The central ray is then an axis or diagonal. For every nonunit h_ij,
its appropriate integer transverse form has operator norm at most
sqrt(2), and its value is nonzero. Consequently

\[
1\leq\sqrt2\frac{M-2}{M}\Delta\,|h_{ij}|.
\]

Multiplying over all L nonunit pairs, and using (2), yields

\[
1\leq
 \left(\sqrt2\frac{M-2}{M}\Delta\right)^L
 R^{M(M-2)/4}.
\tag{6}
\]

If Delta <= C/sqrt(R), the power of R on the right after substitution
is M(M-2)/4 - L/2 <= 0 by (3). Thus such a cluster is impossible whenever

\[
\boxed{C<\frac{M}{\sqrt2(M-2)}.}
\tag{7}
\]

In particular C = 1/2 excludes this branch for every even M >= 4, with
no radius threshold. Complementary valuation rows were included in (3).
For M = 8, condition (5) always holds, recovering the balanced eight-row
exclusion. For general M it is a real restriction: the order of the root
exp(2 pi i l/M) must divide eight. Changing all arc lifts by a common
full turn, or changing the argument chosen for epsilon, does not change
that condition.

### A complementary pair is excluded without a branch condition

If rows i,j are complementary, then h_ij = 1. For every other row t,

\[
\arg h_{it}\equiv(\phi_t-\phi_j)/2\pmod\pi,
\qquad
1\leq|\operatorname{Im}h_{it}|\leq(\Delta/2)|h_{it}|.
\]

There are M-2 such nonunit products. At every prime, exactly k-1 other
rows are on the same side as row i. Consequently their moduli multiply
to R^(k-1), and

\[
1\leq(\Delta/2)^{M-2}R^{(M-2)/2}.
\]

Thus Delta sqrt(R) >= 2. In particular C < 2 excludes complementary
rows in every exactly balanced common-unit cluster, regardless of l.
This is a finite statement, not an asymptotic branch argument.

## 3. All integral linear-phase monomials in the complete half-cut model

The remaining statements concern the complete half-cut profile system
with equal formal weights; this is not an assertion about equal norms
of distinct Gaussian primes.

Index columns by the B = binomial(M,M/2) half-subsets S. Write
sigma_i(S) = 2 indicator(i in S)-1. Suppose an integer exponent vector
v lies in the real row span of this sign matrix. Then, and only then,

\[
\boxed{v_S=\frac12\sum_i n_i\sigma_i(S)
       =\sum_{i\in S}n_i-\frac12\sum_i n_i,}
\tag{8}
\]

where all n_i are integers and their total sum is even. Adding one
common integer to all n_i leaves v unchanged.

To prove necessity, initially write v_S = sum_i b_i sigma_i(S) with
real b_i. Swapping i and j between two half-subsets shows
2(b_i-b_j) is an integer. Set n_i = 2(b_i-b_1). Since every sign column
has sum zero, this gives (8). Its integrality forces sum n_i to be even.
The displayed formula also proves sufficiency.

Give each column weight W/B and define its monomial height by

\[
\mathcal H(v)=\frac WB\sum_S|v_S|.
\tag{9}
\]

For every even M >= 4, the nonzero v satisfying
mathcal H(v) <= W/2 are precisely

\[
\boxed{v=\pm\frac{\sigma_i+\sigma_j}{2}\quad(i\ne j).}
\tag{10}
\]

Their common height is

\[
\boxed{\mathcal H(v)=\frac{M-2}{2(M-1)}W<\frac W2.}
\tag{11}
\]

### Proof of the height classification

Let D = max n_i - min n_i. The same conditional two-index argument
used in the integer-combination note gives

\[
\mathcal H(v)/W\geq\frac{MD}{4(M-1)}.
\tag{12}
\]

If D >= 2, this is strictly greater than 1/2. If D = 0, v vanishes.
Thus a nonzero monomial of height at most W/2 must have D = 1. After
subtracting a constant, n consists of zeros and ones. Its number of
ones is even, say 2j, with 1 <= j <= k-1.

For a uniform half-subset, its intersection size X with these 2j
positions is hypergeometric, and

\[
F_j:=\mathcal H(v)/W=\mathbb E|X-j|
 =\frac{j(k-j)}{k}
   \frac{\binom{k}{j}^2}{\binom{2k}{2j}}.
\tag{13}
\]

For completeness, the probabilities are proportional to
binomial(k,x) binomial(k,2j-x). Put
lambda_x = (k-x)(2j-x) and mu_x = x(k-2j+x). They satisfy
lambda_(x-1) p_(x-1) = mu_x p_x and mu_x-lambda_x = 2k(x-j).
Summing this telescoping identity over x > j, then using symmetry,
gives E|X-j| = j(k-j) p_j/k, which is (13).

The sequence is symmetric under j -> k-j, and

\[
\frac{F_{j+1}}{F_j}
 =\frac{(k-j-1)(2j+1)}{j(2k-2j-1)}\geq1
 \quad\text{when }2j\leq k-1.
\]

For k >= 4 and 2 <= j <= k-2, therefore,

\[
F_j\geq F_2
 =\frac{3(k-1)(k-2)}{(2k-1)(2k-3)}>\frac12.
\]

There are no interior j when k = 2 or 3. Only j = 1 or k-1 remains,
giving (10) and F_1 = (k-1)/(2k-1), as asserted.

## 4. Consequence and limitation of the classification

An integral factor product of height H has modulus exp(H/2). A nonzero
integer projection onto a rational ray, combined with an endpoint angular
error of order exp(-W/4), can force a large-height contradiction when
H < W/2. The classification shows that all products in this height range
whose phases are linear consequences of the balanced row phases are
the pair products (10).

Their central rays have arguments plus or minus 2 pi l/M modulo pi.
If the order of exp(2 pi i l/M) does not divide eight, these rays are
neither axes nor diagonals.
Thus there is no alternative small monomial in this linear phase space
that repairs the rational-ray argument in the complete equal-weight
model. This does not exclude polynomial projections onto irrational
rays, other factor identities, unequal-weight methods, or a different
extraction argument.

The first theorem is a finite Gaussian-integer result. The height
classification is an exact statement about a formal equal-weight
profile model. Neither is the unrestricted endpoint theorem.
