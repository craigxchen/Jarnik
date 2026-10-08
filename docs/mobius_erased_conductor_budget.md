# Erased conductor is controlled by source cut deficit and target residues

The clipped-conductor and same-side-overlap theorems combine to give a
matrix-height-independent radius comparison. For a fixed-unit source of
\(m\ge5\) points with endpoint constant at most two, and an arbitrary-unit
target with endpoint constant at most two,
\[
 \boxed{N'\ge\frac18 N^{1-4/m}.} \tag{1}
\]
If the target pair-residue budget has no extra unit factor, the stronger
pure-power bound \(N'\ge N^{\rho_m}\), with \(\rho_m\ge1-4/m\), holds.
This is another proof of a near-unit radius exponent, not a uniform bound
on cluster cardinality; the earlier projective-orbit work already obtains
that qualitative conclusion by a different route.

Here \(N\) is the least squared radius of the full primitive fixed-unit
source tuple, and \(N'\) is the least squared radius of its full image.
No extraction of target units or redefinition of either tuple is used.
The map is any nonsingular nonconformal rational real-linear map with both
complex coefficients nonzero. Conformal maps preserve least radius and need
no argument. A purely antilinear map is a conjugate rotation and likewise
preserves it.

## 1. Three exact budgets

Write \(W=\log N\), \(W'=\log N'\),
\(P=\prod_{i<j}f_{ij}\), and \(B'=\prod_{i<j}b'_{ij}\). The source
endpoint lower bound at constant \(C\le2\) gives
\[
 \log P\ge\frac{m(m-1)}4W.
\]
At an odd split prime, each unit-depth cut of the source valuation interval
has some \(n\in\{1,\ldots,m-1\}\) rows on one side. Its contribution to
\(\log P\) is \(n(m-n)\log p\). Thus the exact Plotkin deficit is
\[
 D=\frac{m^2}4W-\log P
 =\sum_{\text{source cuts}}(n-m/2)^2\log p
 \le\frac m4W. \tag{2}
\]
This identity uses arbitrary nested valuations, not binary allocations.

Let \(I_p\) be the coefficient interval, and let \(N_{\rm clip}\) be the
conductor from clipping the source valuations to those intervals. By
[mobius_capped_anchor_clipping.md](mobius_capped_anchor_clipping.md),
\[
 N_{\rm clip}\mid\gcd(N,N').
\]
Define the erased logarithmic width
\[
 E_0=W-\log N_{\rm clip}\ge0.
\]
Consequently
\[
 W'\ge W-E_0. \tag{3}
\]

Every erased source cut lies strictly outside its coefficient interval on
one side. If \(n\) source rows lie beyond that cut on the same side, the
same-side overlap theorem contributes \(\binom n2\log p\) to the target
residue budget. Integrating the pairwise overlap over these integer cuts
is precisely the minimum-distance formula from
[mobius_internal_overlap_residues.md](mobius_internal_overlap_residues.md).
Hence
\[
 \log B'\ge
 \sum_{\text{erased source cuts}}\binom n2\log p. \tag{4}
\]
Distance lying outside the full source interval contributes further
nonnegative overlap and can safely be omitted. On erased cuts the
orientation determining \(n\) may differ from the orientation used in
(2), but the square \((n-m/2)^2\) is invariant under \(n\mapsto m-n\).

## 2. An optimized elementary cut inequality

For \(A\ge0\), put
\[
 c(A)=\min_{1\le n\le m-1}
 \left\{\binom n2+A(n-m/2)^2\right\}.
\]
Equations (2) and (4) give
\[
 \boxed{c(A)E_0\le\log B'+AD.} \tag{5}
\]
Suppose for the moment that the target estimate is
\[
 \log B'\le\frac m8W'+\eta. \tag{6}
\]
Combining (2), (3), (5), and (6), whenever \(c(A)>0\), yields
\[
 W'\ge
 \frac{c(A)-Am/4}{c(A)+m/8}W
 -\frac{\eta}{c(A)+m/8}. \tag{7}
\]

The optimal coefficient within this one-parameter charge can be stated
exactly. Let
\[
 k=\lfloor m/4\rfloor,\quad s=m-2k,\quad
 A_m=\frac{k}{s-1},
\]
\[
 c_m=\binom k2+A_m(k-m/2)^2,
 \qquad G_m=c_m+m/8.
\]
At this value the minimum defining \(c_m\) is attained at the two adjacent
integers \(k,k+1\). The resulting ratio is
\[
 \boxed{\rho_m=
 \frac{2k(m-2k-2)}{m+2k(m-2k-2)}.} \tag{8}
\]
Equivalently,
\[
 \rho_m=
 \begin{cases}
 1-4/m,&m\equiv0\pmod4,\\
 1-4m/(m^2+4),&m\equiv2\pmod4,\\
 1-4m/(m^2+3),&m\text{ odd}.
 \end{cases}
\]
In particular \(\rho_m\ge1-4/m\), and \(\rho_m>0\) for \(m\ge5\).

For completeness, the minimum of this family of affine functions of \(A\)
is piecewise affine. On each piece the ratio in (7) has derivative of
constant sign, so a maximum occurs at an adjacent-line breakpoint. The
breakpoint between \(n=j\) and \(n=j+1\) is
\(A=j/(m-2j-1)\), for positive denominator. Substitution gives the ratio
\(2j(m-2j-2)/(m+2j(m-2j-2))\). It is maximized by maximizing
\(j(m-2-2j)\), whose maximizing integer is \(k=\lfloor m/4\rfloor\)
(with a tie when \(m\equiv0\pmod4\)). The boundary \(A=0\) gives zero,
and the limit at infinity is not larger. This proves the claimed optimality
within (5), without asserting optimality among all possible cut estimates.

The exact conclusion with an arbitrary additive target budget is therefore
\[
 \boxed{W'\ge\rho_mW-\eta/G_m.} \tag{9}
\]
Also (5) directly bounds the amount erased:
\[
 E_0\le\frac{mW'/8+A_mmW/4+\eta}{c_m}.
\]
Its coefficients are \(O(1/m)\). Thus erased conductor is bounded by
\(O((W+W')/m+\eta/m^2)\), uniformly in matrix height and prime depths.

## 3. Target units and the clean radius statement

If the target belongs to one suitable unit-parity class, the usual endpoint
product bound is
\[
 B'\le(C'/2)^{\binom m2}(N')^{m/8}.
\]
For \(C'\le2\), one may take \(\eta=0\), proving
\(N'\ge N^{\rho_m}\).

For arbitrary target units, let \(r\) of its primitive half-angle columns
have two odd coordinates. Exactly \(r(m-r)\) pair half-angle norms acquire
the ramified factor two. The odd prime norms still obey Plotkin at the
least target radius. Therefore
\[
 B'\le(C'/2)^{\binom m2}
 2^{r(m-r)/2}(N')^{m/8}.
\]
Thus (9) holds, without target extraction, with
\[
 \eta=\binom m2\log(C'/2)+\frac{r(m-r)}2\log2. \tag{10}
\]
For \(C'\le2\), this is at most \(m^2\log2/8\). Meanwhile
\[
 G_m=\frac{(m-1)[m+2k(s-2)]}{8(s-1)}
 \ge\frac{m^2}{24}\qquad(m\ge5).
\]
Indeed the bracket is at least \(m^2/4\), while
\(s-1\le(m+1)/2\), giving
\(G_m\ge m^2(m-1)/(16(m+1))\ge m^2/24\).
It follows that \(\eta/G_m\le3\log2\), proving (1). Retaining (9)--(10)
gives better explicit constants for any prescribed \(m,r,C'\).

## Scope

The proof couples actual source near-balanced cut budgets to the target
residue budget and the conductor that survives clipping. It addresses
arbitrary coefficient interval locations and lengths; external intervals
and outside-source primes cause only omitted nonnegative terms. It also
accounts for arbitrary target units with a fixed multiplicative loss.

The exponent comparison does not exclude \(N'=N\), nor does it bound
\(m\) independently of radius. Its qualitative near-unit exponent was
already available by the distinct cross-ratio/deletion argument in
[projective_orbit_radius_rigidity.md](projective_orbit_radius_rigidity.md).
The contribution here is the exact cut-deficit/erased-conductor proof and
its explicit constants.

Independent exact checks tested all 4,944 cut charges for \(5\le m\le100\)
and \(1\le n\le m-1\), together with 9,600 deterministic sampled integer
profiles and coefficient intervals, including disjoint intervals and point
intervals. They verified (5), the optimized ratio, and \(G_m\ge m^2/24\).
All checks passed. The proofs above establish the unrestricted statements.
