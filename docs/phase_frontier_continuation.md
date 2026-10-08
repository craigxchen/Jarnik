# Rational-ray obstruction for exactly balanced eight-point systems

## Status

The uniform endpoint theorem is **not proved**. The result below is new
progress for arbitrary Gaussian block values: the complete seven-row,
35-column hypersimplex model requires arc length at least a constant times
\(R^{4/7}\). More generally, an exactly balanced eight-point model with
total Gaussian unit one cannot occupy an arc of length
\(\frac12\sqrt R\). Extraction from a general cluster remains unresolved.

The additional observation is that all four rays \(n\pi/4\) are rational
Gaussian directions. Besides the axes, the diagonal rays are detected by
the integer forms \(a-b\) and \(a+b\).

## 1. Exact balanced model and phase branches

Let \(\mathcal F\) be a nonempty family of four-element subsets of \([7]\),
carrying nonunit Gaussian blocks \(H_S\). Assume

\[
\gcd(H_S,\overline H_T)=1\qquad(S,T\in\mathcal F).
\]

Pairwise coprimality between distinct blocks gives the usual valuation
interpretation, but only conjugate-coprimality is needed in the proof.
Set

\[
G=\prod_SH_S,\quad R=|G|>1,\quad
A_i=\prod_{S\ni i}H_S,\quad B_i=G/A_i.
\]

The points \(z_0=\overline G\) and \(z_i=A_i\overline{B_i}\) all have
absolute value \(R\). Every column has four ones among these eight rows.
Missing four-subset profiles and arbitrary block sizes are allowed. We have

\[
\prod_{i=1}^7A_i=G^4,\qquad
\prod_{i=1}^7B_i=G^3,\qquad
\prod_{i=1}^7|B_i|=R^3. \tag{1}
\]

Suppose the eight distinct points occupy an arc of angular width
\(\Delta<\pi/2\). Keep the original zero row, so the condition that its
complementary products are nonunit does not change. The normalized phases
\(1,A_1/\overline A_1,\ldots,A_7/\overline A_7\) have half-angle lifts
\(\theta_i\in[a,a+d]\), where \(d=\Delta/2\) and
\(-d\leq a\leq0\), since the interval contains the zero-row phase \(0\).
Signs chosen for these lifts are handled modulo \(\pi\); the signed
products need not satisfy (1) literally.

Put \(\mu=\frac14\sum_i\theta_i\). A common \(n\in\{0,1,2,3\}\) satisfies

\[
\arg B_i\equiv n\pi/4+\mu-\theta_i\pmod\pi,\qquad
|\mu-\theta_i|\leq\frac{3\Delta}{4}. \tag{2}
\]

Indeed,
\(\mu-\theta_i=(\sum_{j\ne i}\theta_j-3\theta_i)/4\), so

\[
\frac{3a-3d}{4}\leq\mu-\theta_i\leq\frac{3a+6d}{4}.
\]

Using \(-d\leq a\leq0\) puts both endpoints between
\(-3\Delta/4\) and \(3\Delta/4\). Thus the same constants hold for an
arbitrary original zero row, without re-rooting. When the zero row is an
arc endpoint, one can take \(a=0\) and recover the sharper lower bound
\(-3\Delta/8\), which is not needed below.

For \(B=a+ib\), respectively use the integer forms

\[
\lambda_0(B)=b,\quad\lambda_1(B)=a-b,\quad
\lambda_2(B)=a,\quad\lambda_3(B)=a+b.
\]

Each has operator norm at most \(\sqrt2\) and vanishes on its indicated
ray. The inequality \(|\sin t|\leq|t|\) gives

\[
|\lambda_n(B_i)|\leq\frac{3\sqrt2}{4}\Delta|B_i|. \tag{3}
\]

For nonunit \(B_i\), this integer is nonzero: a Gaussian integer on an
axis is associated to its conjugate; a Gaussian integer on a diagonal is
divisible by \(1+i\), as is its conjugate. Both contradict
\(\gcd(B_i,\overline B_i)=1\).

## 2. Arbitrary-value exponent \(4/7\) obstruction

If all \(B_i\) are nonunit, multiply (3) and use (1):

\[
1\leq\prod_i|\lambda_n(B_i)|
 \leq\left(\frac{3\sqrt2}{4}\Delta\right)^7R^3.
\]

Therefore

\[
\boxed{\Delta\geq\frac{4}{3\sqrt2}R^{-3/7},\qquad
R\Delta\geq\frac{4}{3\sqrt2}R^{4/7}.} \tag{4}
\]

In the complete 35-block model, every \(B_i\) contains 15 nonunit blocks,
so this hypothesis is automatic. For \(\Delta\leq CR^{-1/2}\),

\[
R\leq\left(\frac{3\sqrt2 C}{4}\right)^{14}. \tag{5}
\]

For \(C<4/(3\sqrt2)\), (5) is impossible at every \(R>1\).
Unlike the fixed-polynomial obstruction of research-note §6.7, this uses
no parameterization or degree hypothesis.

## 3. A unit complementary product

If \(B_i=1\), then \(A_i=G\). Distinctness permits at most one such row in
the common-unit model. For every other \(j\),

\[
\arg B_j\equiv\theta_i-\theta_j\pmod\pi,\qquad
1\leq|\operatorname{Im}B_j|\leq\frac{\Delta}{2}|B_j|.
\]

There are six remaining rows, and their \(|B_j|\) still have product
\(R^3\). Multiplication yields

\[
1\leq(\Delta/2)^6R^3,\qquad\Delta\geq2R^{-1/2}. \tag{6}
\]

Thus \(C=1/2\) excludes any exactly four-balanced eight-point common-unit
model, including models with missing profiles and unequal block sizes.

## 4. Independent Gaussian units: audit and extension

For independent normalized units, write the ratios as
\(\eta_iA_i/\overline A_i\), with \(\eta_i=i^{e_i}\) and
\(e_i\in\{0,1,2,3\}\). The half-angle lifts satisfy

\[
\arg A_i\equiv\theta_i-\pi e_i/4\pmod\pi.
\]

For \(E=\sum_i e_i\), the exact generalization of (2) is

\[
\arg B_i\equiv\mu-\theta_i+\frac{\pi e_i}{4}
 -\frac{\pi E}{16}+\frac{n\pi}{4}\pmod\pi. \tag{7}
\]

Arbitrary units can produce rays at multiples of \(\pi/16\), and sign
units can still produce multiples of \(\pi/8\). Those need not be rational
Gaussian directions. Silently absorbing independent units into \(A_i\)
would invalidate an unrestricted rational-ray proof.

However, the condition

\[
\prod_i\eta_i=1,\qquad\text{equivalently }E\equiv0\pmod4 \tag{8}
\]

makes every ray in (7) an integer multiple of \(\pi/4\). The ray normal
may depend on \(i\), but its norm is at most \(\sqrt2\). The same product
argument proves (4)--(5). Thus total unit product one suffices; individual
units need not all coincide.

If some \(B_i=1\), no total-unit hypothesis is needed. Each remaining
\(B_j\) lies within \(\Delta/2\) of the rational ray
\((e_j-e_i)\pi/4\). At \(\Delta<\pi/2\), there is at most one unit
complementary product: identical valuation rows with different units give
angular separation at least \(\pi/2\). Consequently

\[
1\leq(\Delta/\sqrt2)^6R^3,\qquad
\Delta\geq\sqrt2R^{-1/2}. \tag{9}
\]

So \(C=1/2\) excludes all exactly balanced eight-point systems satisfying
(8), and any such system with a unit complementary product regardless of
total unit. The cases with unit product \(-1,i,-i\) and every \(B_i\)
nonunit remain unsettled.

The product of the eight relative units is invariant under changing the
zero row, because the eighth power of a Gaussian unit is one. Reorienting
cuts conjugates block representatives without introducing independent
point units.

## 5. General branch statement and the remaining gap

For \(k\geq2\), an exactly \(k\)-balanced common-unit model on \(m+1=2k\) points satisfies the following. The
same argument gives

\[
\prod_i|B_i|=R^{k-1},\quad
\arg B_i\equiv n\pi/k+\varepsilon_i\pmod\pi,\quad
|\varepsilon_i|\leq\frac{k-1}{k}\Delta. \tag{10}
\]

Whenever \(k\mid4n\), this is a rational Gaussian ray. If all
complementary products are nonunit,

\[
\Delta\geq\frac{k}{\sqrt2(k-1)}
 R^{-(k-1)/(2k-1)}. \tag{11}
\]

For \(k=2\), both rays are axes and the \(\sqrt2\) can be removed.
For \(k=4\), every ray is rational, giving the new arbitrary-value result.
For \(k=3\), the cubic/Pell branches remain. Some branches remain for
\(k=5\) and \(k=8\) too.

Restricting a larger fair weak-flip family to eight rows creates profiles
of sizes other than four. Then \(\prod A_i=G^4\) acquires a nontrivial
Gaussian correction, moving the ray. Neither the present argument nor
the existing weighted defect estimate supplies the needed extraction or
stability theorem. The full endpoint bound remains open in this work.

## 6. Formalization boundary

[GaussianChain/RationalPhaseObstruction.lean](../GaussianChain/RationalPhaseObstruction.lean)
formalizes the arithmetic product step with distinct integer ray normals.
Its theorem seven_factor_radius_bound derives

\[
R\leq(9C^2/8)^7=(3\sqrt2 C/4)^{14}
\]

from seven nonzero integer transverse forms, their squared upper bounds,
and the product of their norm squares \(R^6\). This arithmetic theorem has
been checked by the Lean compiler. The Gaussian construction, phase
branches, conjugate-coprimality nonvanishing, trigonometric estimates, and
unit-complement case are proved in prose above; they are not claimed to
have been formalized in that file.


## 7. Why direct fourth-power clearing does not give stability

Allow general column sizes \(r_S=|S|\) on the seven nonzero rows, rather
than imposing \(r_S=4\). The exact product identity becomes

\[
\prod_{i=1}^7A_i=G^4E,\qquad
E=\prod_SH_S^{\,r_S-4}\in\mathbf Q(i)^\times. \tag{12}
\]

The correction \(E\) changes the common phase ray. One direct attempt to
clear its fourth root is to use the Gaussian rational numbers

\[
EB_i^4.
\]

Their angles are close to integer multiples of \(\pi\), since
\(\arg(EB_i^4)\equiv\sum_j\theta_j-4\theta_i\pmod\pi\) in the
common-unit model. This error has absolute value at most \(3\Delta\).

After cancelling common Gaussian factors in the norm-one ratio of
\(EB_i^4\), the signed exponent at block \(S\) is

\[
f_i(S)=(r_S-4)+4(1-\mathbf1_{\{i\in S\}})
      =r_S-4\mathbf1_{\{i\in S\}}. \tag{13}
\]

Thus its Gaussian denominator size is exactly

\[
\prod_S|H_S|^{\,|f_i(S)|},
\]

under the pairwise and conjugate-coprimality assumptions. Summing the
absolute exponent loads over the seven rows gives

\[
\sum_{i=1}^7|f_i(S)|
=r|r-4|+(7-r)r
=\begin{cases}
11r-2r^2,&0\leq r\leq4,\\
3r,&4\leq r\leq7.
\end{cases} \tag{14}
\]

Already for balanced columns \(r=4\), this sum is \(12\), rather than the
\(3\) obtained from the unpowered complementary products. Multiplication
of seven elementary nonzero-imaginary-part bounds therefore pays
denominator \(R^{12}\), while the seven endpoint errors contribute only
\(R^{-7/2}\). Even the exactly balanced case gives no contradiction from
that estimate. If one of these rational numbers is real, its corresponding
nonzero lower bound is unavailable, making this repair still weaker.

This calculation rules out only the direct fourth-power/root-clearing
repair just described. It is not a proof that every stability argument
fails. A different method would need to retain the arithmetic benefit of
the rational rays without multiplying all complementary exponents by four.
