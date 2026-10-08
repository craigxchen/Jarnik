# Endpoint maps with controlled radius must approach rank one

The erased-conductor budget and the complete determinant-core bound together
give a geometric restriction on the map itself. For at least eight source
points, a nonconformal rational map between two endpoint configurations with
\(N'\le N\) must have condition number growing by an explicit positive
power of \(N\). The constants are uniform in source configurations, matrix
height, determinant, prime powers, and target Gaussian units. Section 7
extends the conclusion to arbitrary source units as well, retaining the full
point count with an explicit absolute constant.

This is a conditional classification of admissible maps. It is not a uniform
bound on the number of lattice points in an endpoint arc.

## 1. Hypotheses and the invariant ratio

Assume a primitive fixed-unit source of \(m\ge5\) points, of least squared
radius \(N\), satisfies the endpoint arc bound with \(C\le2\). Its full
image under a nonsingular rational real-linear map has least squared radius
\(N'\) and satisfies the endpoint bound with \(C'\le2\). The target may
have arbitrary Gaussian units; no target subset is extracted.

Write
\[
 Mz=\alpha z+\beta\bar z,\qquad \alpha\beta\ne0,
 \qquad N(\beta/\alpha)=a/b,
\]
with coprime positive integers \(a\ne b\). Retain the exact normalization
quantities
\[
 \tau=\epsilon(a+b)/2,\qquad c=\epsilon|a-b|/4>0,
 \qquad\epsilon\in\{1,2,4\}.
\]
At every actual source anchor, \(T=\tau k\) and \(|\delta|=ck\).
Consequently
\[
 \boxed{\mathcal R:=\frac{T}{|\delta|}=\frac\tau c
 =\frac{2(a+b)}{|a-b|}.}
 \tag{1}
\]
This ratio is unchanged by arbitrary nonzero scalar normalization and by
input or output conformal maps, including reflections. It is therefore an
intrinsic geometric quantity, not the height of a chosen integer matrix.

Conformal and anticonformal maps are excluded from the hypotheses above;
they preserve least radius and have condition number one. Their existence
does not conflict with the nonconformal conclusion.

## 2. Two established budgets, including all parity factors

Use the constants from
[the erased-conductor budget](mobius_erased_conductor_budget.md):
\[
 k=\lfloor m/4\rfloor,\quad s=m-2k,\quad
 A_m=\frac{k}{s-1},\quad
 c_m=\binom k2+A_m(k-m/2)^2,
\]
and put
\[
 F_m=\left\lfloor\frac{(m-1)^2}{4}\right\rfloor,
 \qquad W=\log N,\quad W'=\log N'.
\]
The subscripted cut constant \(c_m\) is distinct from the determinant core
\(c\). If \(E_0=W-\log N_{\rm clip}\), the erased-conductor theorem gives
\[
 E_0\le\frac{mW'/8+A_mmW/4+\eta}{c_m}.
 \tag{2}
\]
For arbitrary target units and \(C'\le2\), the explicit parity cut permits
\(\eta\le m^2\log2/8\).

The full-core theorem in
[joint source-target residue control](joint_source_target_residue_core.md),
including its ramified-prime calculation, gives
\[
 \log c\le
 \frac{m}{8F_m}(W+W')+
 \frac{m(m-1)}{2F_m}\log2.
 \tag{3}
\]
Using uniform parity bounds separately in (2)--(3) avoids any need to choose
a common target parity cohort or to match parity counts across different
output normalizations. No two-part of \(c\) is omitted.

At every odd split prime, the clipped source width is at most the length
of the coefficient interval, which equals \(|v_p(a/b)|=v_p(ab)\). Hence
\[
 \boxed{N_{\rm clip}\mid ab.}
 \tag{4}
\]
This uses coprimality of \(a,b\); inert and ramified source primes are
absent under the fixed-unit odd-norm source normalization. It follows that
\[
 \tau\ge\epsilon\sqrt{ab}\ge\epsilon\sqrt{N_{\rm clip}}.
 \tag{5}
\]

## 3. A two-ended bound without a radius-descent assumption

Define explicitly
\[
 \alpha_m=\frac12-\frac{A_mm}{8c_m}-\frac{m}{8F_m},\qquad
 \beta_m=\frac{m}{16c_m}+\frac{m}{8F_m},
\]
\[
 H_m=\frac{m^2}{16c_m}+\frac{m(m-1)}{2F_m}.
\]
Combining (2)--(5) gives
\[
 \boxed{
 \mathcal R\ge
 \epsilon\,2^{-H_m}N^{\alpha_m}(N')^{-\beta_m}.}
 \tag{6}
\]
Indeed, take one half of the lower bound on \(\log N_{\rm clip}=W-E_0\)
and subtract the upper bound (3) on \(\log c\). The coefficient of \(W\)
is \(\alpha_m\), the coefficient of \(W'\) is \(-\beta_m\), and the
two parity charges sum to \(H_m\log2\).

Thus even without \(N'\le N\), a relation \(N'\le N^\theta\) forces
\(\mathcal R\) to grow whenever \(\alpha_m-\theta\beta_m>0\).
For example, at \(m=8\),
\[
 \alpha_8=31/132,\qquad\beta_8=29/132,
\]
so any fixed \(\theta<31/29\) gives a positive exponent. Formula (6)
retains the full two-radius information rather than assuming comparable
radii in advance.

## 4. Explicit positive exponent for every m at least eight

If \(N'\le N\), define
\[
 q_m=\frac{m(1+2A_m)}{8c_m}
 =\frac{m(m-1)}{2k[(m-1)(m-2k-1)+1]},
\]
\[
 z_m=\alpha_m-\beta_m
 =\frac{1-q_m}{2}-\frac{m}{4F_m}.
\]
Then
\[
 \boxed{\mathcal R\ge\epsilon\,2^{-H_m}N^{z_m}.}
 \tag{7}
\]
Some exact exponents are
\[
\begin{array}{c|rrrrrr}
 m&5&6&7&8&9&10\\\hline
 z_m&-53/144&-7/32&-103/900&1/66&61/704&3/23.
\end{array}
\]
For this stated cut choice, eight is the first positive threshold, and the
exponent is positive for every \(m\ge8\). Here is a uniform verification.
For fixed \(k\), the displayed rational function \(q_m\) decreases on
\(4k\le m\le4k+3\). At the block start,
\[
 q_{4k}=\frac{4k-1}{4k^2-3k+1},
\]
which decreases for \(k\ge2\). Hence \(q_m\le7/11\) for \(m\ge8\).
Also \(m/(4F_m)\le1/6\) in this range. Therefore
\[
 \boxed{z_m\ge\frac12-\frac7{22}-\frac16=\frac1{66}>0.}
 \tag{8}
\]
As \(m\to\infty\), \(q_m=4/m+O(m^{-2})\) and
\(m/(4F_m)=1/m+O(m^{-2})\), so
\[
 z_m=\frac12-\frac3m+O(m^{-2}),\qquad H_m=3+O(m^{-1}).
\]

There is also a simple absolute constant valid simultaneously for all
\(m\ge8\). The erased-budget estimate
\(c_m+m/8\ge m^2/24\) gives
\(c_m\ge m(m-3)/24\), so
\[
 H_m\le\frac{3m}{2(m-3)}+\frac{m(m-1)}{2F_m}
 \le\frac{12}{5}+\frac73<5.
\]
Since \(\epsilon\ge1\), (7)--(8) imply
\[
 \boxed{\mathcal R\ge\frac1{32}N^{z_m}
 \ge\frac1{32}N^{1/66}.}
 \tag{9}
\]
The exact constants in (7) are preferable when a particular \(m\) is fixed.

## 5. Condition number and closeness of the coefficient norms

Let \(\kappa=\sigma_{\max}(M)/\sigma_{\min}(M)\ge1\) be the Euclidean
condition number. The identity
\[
 \mathcal R=\kappa+\kappa^{-1}
\]
implies \(\kappa\ge\mathcal R/2\). Thus under the descent hypotheses and
\(m\ge8\),
\[
 \boxed{\kappa\ge\frac1{64}N^{z_m},\qquad
 \frac{\sigma_{\min}}{\sigma_{\max}}\le64N^{-z_m}.}
 \tag{10}
\]
After normalizing the largest singular value to one, the distance in
operator norm to the set of rank-one matrices is exactly
\(\sigma_{\min}/\sigma_{\max}\). This gives a precise meaning to the
statement that such maps approach rank one as the source radius grows.

The same conclusion can be written solely in terms of the rational norm
ratio:
\[
 \boxed{\frac{|a-b|}{a+b}=\frac2{\mathcal R}
 \le64N^{-z_m}.}
 \tag{11}
\]
For \(r_0=\min(a,b)/\max(a,b)\), this gives
\(1-r_0\le128N^{-z_m}\). If one keeps the unsymmetrized ratio \(a/b\),
then whenever \(\mathcal R\ge4\),
\[
 |1-a/b|\le8/\mathcal R\le256N^{-z_m}.
\]
These bounds concern a near-singular norm balance \(a/b\to1\), not
closeness to a conformal map.

## 6. Relation to earlier projective-orbit restrictions

[Projective-orbit radius rigidity](projective_orbit_radius_rigidity.md)
already gives a near-unit relation between the logarithms of endpoint radii,
without restricting matrix height or determinant. The present consequence
does not claim a new qualitative radius-comparison principle. Instead it
combines the newer erased-conductor and complete-core budgets to constrain
the invariant condition number of a nonconformal map, uniformly while the
source configuration moves.

The earlier [fixed-six obstruction](mobius_fixed_six_obstruction.md) uses
the condition number to compare the image angular scale with the diameter
of a fixed source configuration. Its constants depend on that source.
Equations (6)--(10) have no such fixed-source constants.

Near-rank-one rational maps still exist at arbitrary arithmetic height.
Nothing here proves that they cannot map large endpoint configurations to
other endpoint configurations, or excludes \(N'=N\). A uniform point-count
conclusion remains unproved.

## 7. Arbitrary source units without losing any points

The [arbitrary-unit audit](mobius_arbitrary_unit_retention_addendum.md)
extends both budgets to a full primitive source configuration. Here are the
normalizations and constants needed to extend the condition-number theorem.
Keep all \(m\) source points and all \(m\) target points. Choose one source
phase anchor \(H_0=1\), and fix its primitive integral map normalization.
Let \(\sigma\) source half-angle rows and \(r\) target half-angle rows have
even norms. The source half-angle rows need not have odd norms.

The least primitive source and target circle squared radii remain odd:
an even squared radius would force every Gaussian circle point to be
divisible by \(1+i\), contradicting primitivity of the whole realization.
Mixed half-angle parities therefore contribute a factor two to certain
relative pair norms, not to the least circle radius. All of \(\log N\)
is still accounted for by the odd split-prime source widths.

The full source pair product has the exact parity factor
\(2^{\sigma(m-\sigma)}\). Hence its odd-prime Plotkin deficit satisfies
\[
 D_{\rm odd}\le mW/4+\sigma(m-\sigma)\log2.
\]
The same odd-prime erased-cut proof consequently gives
\[
 E_0\le\frac{mW'/8+A_mmW/4+
 [r(m-r)/2+A_m\sigma(m-\sigma)]\log2}{c_m}.
 \tag{12}
\]
The full-core audit, including both source and target Gaussian gcd parities
at two, proves
\[
 c^{F_m}\mid
 2^{\binom\sigma2+\binom r2}BB'.
\]
Combining with the two endpoint residue products gives
\[
 \log c\le\frac{m}{8F_m}(W+W')+
 \frac{(\sigma+r)(m-1)}{2F_m}\log2
 \le\frac{m}{8F_m}(W+W')+
 \frac{m(m-1)}{F_m}\log2.
 \tag{13}
\]
This retains the full two-part of the determinant core. Even-norm source
reanchoring can change \(\epsilon\); no constant-\(\epsilon\) all-anchor
claim is needed here. The proof uses one fixed phase anchor, the intrinsic
ratio \(\mathcal R\), and the odd coefficient intervals only.

The divisibility \(N_{\rm clip}\mid ab\) and inequality (5) remain valid.
Combining them with (12)--(13) gives exactly the same radius exponents as
before, with the enlarged explicit parity constant
\[
 H_m^{\rm all}=
 \frac{m^2(1+2A_m)}{16c_m}+\frac{m(m-1)}{F_m}.
\]
Thus, without assuming radius descent,
\[
 \boxed{\mathcal R\ge
 \epsilon\,2^{-H_m^{\rm all}}N^{\alpha_m}(N')^{-\beta_m}.}
 \tag{14}
\]
If \(N'\le N\), replace the radius factor by \(N^{z_m}\), with the same
positive threshold \(m\ge8\) and the same \(z_m\ge1/66\).

For completeness, \(H_m^{\rm all}<8\) for every \(m\ge8\). Put
\(m=4k+j\), \(0\le j\le3\), and \(d=m-2k-1\). The exact formula
\[
 c_m=\frac{k(m-1)}4+\frac{k}{4d}
\]
gives
\[
 16\left(c_m-\frac{m(m-4)}{16}\right)
 =4k(3-j)+j(4-j)+4k/d>0.
\]
Also
\[
 m+2k(m-2k-2)=m^2/4+j(4-j)/4\ge m^2/4.
\]
Using the exact expression for \(G_m=c_m+m/8\), the latter inequality
implies
\[
 \frac{m^2(1+2A_m)}{16c_m}
 \le\frac{2G_m}{c_m}
 =2+\frac{m}{4c_m}
 \le2+\frac4{m-4}\le3.
\]
The remaining term is at most \(14/3\) for \(m\ge8\), so
\(H_m^{\rm all}\le23/3<8\).

Consequently, for the full arbitrary-unit source and target configurations,
\[
 \boxed{\begin{aligned}
 \mathcal R&\ge N^{z_m}/256,\\
 \kappa(M)&\ge N^{z_m}/512,\\
 \frac{\sigma_{\min}(M)}{\sigma_{\max}(M)}
 &=\frac1{\kappa(M)}\le512N^{-z_m},\\
 \frac{|a-b|}{a+b}&\le512N^{-z_m}.
 \end{aligned}}
 \tag{15}
\]
There is no unit-class reduction in the value of \(m\). This extends the
geometric classification to the full point count while preserving the
limitations stated in Section 6.

Root independently checked the algebra and the same exponent and parity
constants through \(m=100\) in
[the persistent joint checker](check_mobius_joint_anchor_residues.py),
whose 2,187 configurations include mixed source and target half-angle
parities. Astra also checked the exact constants through \(m=10{,}000\).
These finite checks supplement the unrestricted proofs above.
