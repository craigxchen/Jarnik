# Uniform counts with at most three affine classes per frequency bucket

Three nonproportional affine Fibonacci classes can have incomparable
slopes and genuinely overlapping frequencies. Nevertheless the same
absolute template count holds when each 2-adic slope bucket contains
at most three original classes. At a frequency where at most two
classes contribute, quadratic irrationality separates their moments
without any coefficient-size threshold. This supplies a packet
promotion that is unavailable from slope divisibility alone.

The count concerns fixed growing affine Fibonacci templates. The
general mixed-class model and the arbitrary-circle theorem are not
proved here.

## 1. Hypotheses and statement

Use the original rows and normalization of
[the divisibility-chain packet theorem](divisibility_chain_affine_uniform_count.md).
There are finitely many nonproportional affine class representatives

\[
d_c(n)=A_cn+B_c,\qquad A_c>0\text{ even},\quad B_c\text{ odd},
\]

on a fixed parity subsequence. Their integer layer counts have actual
widths \(W_{c,e}\), with odd indices \(e\), and positive weighted
degree

\[
L=\sum_{c,e}A_c\phi(e)W_{c,e}>0.
\]

Assume fixed equal-norm Gaussian prefactors, aligned limiting phases,
and distinct eventual rows on a fixed \(C\sqrt{R_n}\) arc after
division by the entire Gaussian gcd. A family of at most one row is
trivial. Otherwise let \(\tau\) be the least nonzero pair contact
frequency. The original geometry gives \(L\le4\tau\), common total
base parity, and the global reciprocal/obtuse count in terms of the
original affine dimension.

Assume each bucket of classes with the same \(v_2(A_c)\) contains
at most three **original classes**, not merely three distinct slope
values. The classes in such a bucket may have arbitrary incomparable
slopes and arbitrary distinct ratios \(B_c/A_c\).

**Theorem.** With \(D_0=2^{32}\), the original row affine dimension
is at most \(D_0+1125\), and the number of rows is at most

\[
2(D_0+1125)+1=8589936843.                               \tag{1}
\]

One may mix these buckets with those of the preceding theorem: a
bucket may instead have arbitrarily many classes if its distinct
slopes form a divisibility chain or if it has at most two distinct
slope values. The proof processes the buckets separately and charges
their total lost dimension only once.

## 2. Two contributing classes separate exactly

For a pair difference write

\[
S_{c,k}=\sum_e b_{c,e}c_e(k),\qquad
\gamma_c=i^{B_c}(-1)^{A_cn/2}\varphi^{-B_c},
\quad \varphi=(1+\sqrt5)/2.
\]

At a positive physical frequency \(t\), only classes with
\(A_c\mid t\) and \(k_c=t/A_c\) odd contribute. After clearing
the logarithmic denominator, their equation is

\[
\sum_c A_c S_{c,k_c}\gamma_c^{k_c}=0.                  \tag{2}
\]

For contributing classes,

\[
\gamma_c^{k_c}
 =(-1)^{tn/2}i^{m_c}\varphi^{-m_c},
\qquad m_c=B_ct/A_c\text{ odd}.                        \tag{3}
\]

Distinct nonproportional classes have distinct \(m_c\). Thus the
ratio of the phase constants of two contributing classes is

\[
\pm\varphi^{m_{c'}-m_c},
\]

where the exponent is a nonzero even integer. This ratio is
irrational. For example, its quadratic conjugate is its reciprocal;
if it were rational it would equal that reciprocal, forcing modulus
one, contrary to the nonzero exponent.

Consequently if at most two classes contribute to (2), every integer
coefficient \(A_cS_{c,k_c}\) is zero. This assertion has no size
condition on the coefficients, slopes, or frequency. It also holds
with negative odd intercepts. The fixed-parity sign in (3) is common
to the contributing classes and does not affect the argument.

## 3. A promotion that excludes one of three colors

Retain original class labels throughout all packet operations. A class
is rescaled only as \((A_c,B_c)\mapsto(pA_c,pB_c)\), so its ratio
\(B_c/A_c\) is invariant. Distinct original colors remain
nonproportional. A class with no active coordinates may be discarded.

Suppose a bucket has more than one current slope. Let \(A\) be the
numerically smallest one and choose any larger current slope \(B\).
There is an odd prime \(p\) with

\[
v_p(A)<v_p(B).                                         \tag{4}
\]

Otherwise \(B\mid A\), contradicting \(A<B\); the common 2-adic
valuation excludes the prime two from this argument.

At a frequency \(Ak\) with \(k\) odd and \(p\nmid k\), the class
chosen at slope \(B\) cannot contribute. At most two of the original
colors in this bucket remain. Other 2-adic buckets have disjoint
odd-harmonic frequencies. For \(Ak<\tau\), the low coefficient is
zero, and Section 2 separates it. Therefore each class at slope \(A\)
satisfies

\[
S_{c,k}=0\quad(k<h_A\text{ odd},\ p\nmid k),
\qquad h_A=\min\{h\ge1:h\text{ odd},\ Ah\ge\tau\}.
                                                               \tag{5}
\]

Apply the integral single-prime packet projection to every class at
\(A\), as proved in Section 3 of the preceding note. Specifically,
the difference of its entire level-zero and level-one \(p\)-slices
has every low moment zero. Equalize those slices by copying the one
with smaller aggregate totient width, using one choice for the whole
class. The identities

\[
\Psi_s(z)\Psi_{ps}(z)=\Psi_s(z^p),\qquad
\Psi_{p^js}(z)=\Psi_{p^{j-1}s}(z^p)\quad(j\ge2)
\]

then rescale the class to \((pA,pB_c)\). Its integer width cost does
not increase, and each original color's physical pair coefficients
below \(\tau\) are preserved individually.

Unlike the arithmetic separation used in the preceding theorem,
this step does **not** require \(L/A\ge D_0\). The exact
two-term irrationality argument supplied (5) directly. Thus it is
legitimate to equalize these buckets even when an intermediate slope
already exceeds the usual budget threshold.

## 4. Finite equalization and the final cutoff

Let \(H\) be the least common multiple of the bucket's original
slopes. During the preceding process every current slope divides
\(H\). Indeed, in (4) the new \(p\)-valuation of \(pA\) is no
larger than the valuation of the already present slope \(B\).
Every processed original color has its slope strictly increased by
a factor of at least three and bounded by \(H\). There are finitely
many original colors. Hence after finitely many steps the bucket has
at most one slope value, unless all its coordinates have disappeared.
This termination argument remains valid when a projected color
disappears.

Now all surviving classes in the bucket share a slope \(A\). If
\(A>L/D_0\), leave them unchanged. Otherwise use the one-excluded-prime
arithmetic lemma and the single-prime projection from the preceding
theorem, with \(p=3\), upper budget \(\lceil L/A\rceil\), and cutoff
\(h_A\). Repeat until the surviving slope exceeds \(L/D_0\) or all
coordinates disappear. No new class labels are created.

The same final cutoff is reached in divisibility-chain buckets and
two-slope buckets by their existing procedures. Across all buckets,
the integer width cost remains at most the fixed original \(L\).
Thus fewer than \(D_0\) final active coordinates survive.

The intermediate rows may lose total base parity or actual Gaussian
circle realization. Neither is used here. The preserved data are
integrality, width cost, and each original color's low physical pair
coefficients. The original geometric embedding is used only at the end.

## 5. One composed-kernel bound

Let \(V\) be the span of the original row differences, and \(T\) the
composition of the fixed packet coordinate maps. All statements about
preserved coefficients hold linearly on \(V\). If \(Tv=0\), every
original color's original moments vanish below \(\tau\). Hence

\[
\ker(T|_V)\subseteq\bigoplus_c\ker M_c,
\qquad
M_c[d,e]=\mu(e/d)\mathbf1_{d\mid e},\quad
d<h_{A_c}\text{ odd},                                 \tag{6}
\]

with original slopes and original active supports on the right.
Ramanujan inversion gives this implication without a coprimality
assumption on the layer indices.

The original support budget is at most \(L\le4\tau\). The weighted
colored root-group compression therefore bounds the sum of the
nullities in (6) by \(1125\), independently of the number of buckets
or promotions. The image dimension is less than \(D_0\). This proves
the affine bound in Section 1. The original global total-base parity
and sigma-reciprocal degree inequality give pairwise nonpositive inner
products in the original weighted embedding, so its affine dimension
bounds the number of distinct rows by \(2(D_0+1125)+1\).

More generally, Section 3 is valid whenever at most two original
colors contribute at every relevant off-\(p\) frequency for the
selected slope, even if the bucket has more colors in total. It is
the bound on contributors, not merely a bound on slope values, that
permits exact separation without an integer-size hypothesis.

The formerly obstructing three slopes \(6,10,14\) are now covered
when there is at most one original class at each slope. Arbitrarily
many classes at those three slopes are not covered by this argument.
General incomparable mixed templates and the full radius-independent
circle theorem remain unresolved.

The [exact checker](check_three_class_affine_promotion.py) verifies
quadratic-field independence of two contributing phases, denominator
clearing, and all promotion paths for 1,680 selected small slope
triples. The general statements rest on the proofs above.
