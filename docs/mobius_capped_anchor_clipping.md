# Capped anchor normalization and the clipped source conductor

This note proves three further consequences of the exact anchor formulas:
output-conformal minimization still optimizes the source-capped bounds;
clipped source pair norms divide the corresponding target pair norms; and
the capped all-anchor inequality combines with overlap residues to eliminate
the entire product of anchor scale factors. Full interval coverage therefore
rules out strict radius descent. The remaining cases do not yield a uniform
endpoint bound here.

Retain a primitive odd-norm fixed-unit source tuple with anchor \(H_0=1\)
and least squared radius \(N\). Let \(M\) be a nonsingular nonconformal
rational real-linear map. In the exact actual-anchor notation, set
\[
 \tau=\epsilon(a+b)/2,\qquad K_c=\epsilon^2ab,\qquad
 c=\epsilon(b-a)/4\ne0,
\]
so the minimal integral representative at actual anchor \(j\) has
\[
 T_j=\tau k_j,\qquad K_j=K_ck_j^2,\qquad\delta_j=ck_j.
\]
All these quantities are as defined in
[actual anchor normalization](mobius_actual_anchor_normalization.md).

## 1. Output minimization also optimizes the source cap

For an integral matrix with invariants \((T,K,\delta)\), define
\(J=\gcd(K,|\delta|N)\). Output multiplication by a Gaussian integer
\(\lambda\ne0\), of norm \(s\), gives invariants
\((sT,s^2K,s\delta)\). Its new source cap is exactly
\[
 \boxed{J(\lambda M)=s\gcd(sK,|\delta|N),\qquad sJ\mid J(\lambda M).}
 \tag{1}
\]
Every integral representative of a fixed output-conformal class is such a
Gaussian padding of its minimal representative, by the normalization ideal
theorem. Its good determinant part divides that of the minimal representative:
primes of \(s\) enter the discriminant and cease to be good, while all
other determinant valuations are unchanged. Therefore both
\[
 \frac1{T^mJ^{m-2}},\qquad
 \frac{D^{2m-3}}{T^{m+2}J^{m-2}}
\]
are maximized at the same minimal representative. Passing to the maximum of
these bounds does not change the conclusion.

At actual anchor \(j\), the exact source cap consequently is
\[
 \boxed{J_j=k_jA_j,\qquad A_j=\gcd(K_ck_j,|c|N).}
 \tag{2}
\]
In particular \(A_j\mid |c|N\), rather than requiring an unrestricted
power of \(K_c\).

## 2. Clipped pair distances survive in the target

At an odd split prime \(p=\pi\bar\pi\), write
\[
 d=v_\pi(w),\quad e=v_{\bar\pi}(w),\quad
 t_i=v_\pi(H_i)-v_{\bar\pi}(H_i),\quad
 I=[\min(d,-e),\max(d,-e)],
\]
where \(w=\beta/\alpha\). Let \(z_i=\operatorname{clip}_I(t_i)\),
and let \(t_i'\) be the signed target phase valuations in any common output
normalization. Then
\[
 \boxed{|t_i'-t_j'|\ge|z_i-z_j|\quad\text{for every pair }i,j.}
 \tag{3}
\]

Here is a complete endpoint-cancellation audit. Put \(r_i=H_i/\bar H_i\).
After a common output rotation, the target phase is
\[
 r_i\longmapsto\frac{r_i+w}{\bar w r_i+1}.
\]
Ignoring the common signed output shift, its valuation is the difference
between the numerator and denominator valuations.

If \(d+e>0\), the interval is \([-e,d]\). Away from its endpoints, the
output valuation is exactly \(\operatorname{clip}_{[-e,d]}(t_i)\).
At the lower endpoint \(t_i=-e\), only the denominator can cancel, giving
\(t_i'\le-e\). At the upper endpoint \(t_i=d\), only the numerator can
cancel, giving \(t_i'\ge d\). Thus endpoint exceptions move outward.

If \(d+e<0\), the interval is \([d,-e]\). Away from its endpoints, the
output valuation is
\[
 d-e-\operatorname{clip}_{[d,-e]}(t_i).
\]
At the lower endpoint \(t_i=d\), numerator cancellation gives
\(t_i'\ge-e\); at the upper endpoint \(t_i=-e\), denominator cancellation
gives \(t_i'\le d\). These are again outward deviations, now after a
reflection.

For two distinct clipped values, their order is therefore preserved or
reversed with their separation at least unchanged. If the clipped values
coincide, the right side of (3) is zero and no cancellation restriction is
needed. Finally, if \(d+e=0\), the interval is a point and (3) is trivial.
This proves the assertion in all cases, including source valuations exactly
at interval endpoints. A common output shift does not affect pair distances.

Define the clipped pair norms
\[
 C_{ij}=\prod_{p\equiv1\pmod4}p^{|z_{i,p}-z_{j,p}|},
\]
and let \(f_{ij}'\) be the primitive target relative half-angle norm.
At odd split primes, its exponent is \(|t_i'-t_j'|\). Consequently
\[
 \boxed{C_{ij}\mid\operatorname{oddpart}(f_{ij}').}
 \tag{4}
\]
No target unit-class assumption is required. Possible target parity factors
cannot interfere with this odd divisor.

## 3. A source conductor divides the target radius

Let
\[
 S_p=[\min_i t_{i,p},\max_i t_{i,p}],\qquad
 V_p=\max_i z_{i,p}-\min_i z_{i,p},\qquad
 N_{\rm clip}=\prod_p p^{V_p}.
\]
The integer \(V_p\) is the length of \(I_p\cap S_p\) if the intervals
meet, and is zero otherwise. Since clipped widths cannot exceed source
widths, \(N_{\rm clip}\mid N\). By (3), the target valuation width is at
least \(V_p\). Therefore
\[
 \boxed{N_{\rm clip}\mid\gcd(N,N').}
 \tag{5}
\]

In particular, if the coefficient interval contains the entire source
valuation range at every prime dividing \(N\), then
\[
 \boxed{N\mid N'.}
 \tag{6}
\]
Thus this coverage regime admits no strict radius descent. If \(N'\le N\),
it forces \(N'=N\). This conclusion is stronger than an archimedean height
restriction and does not require the determinant core \(c\) to be small.
Coefficient behavior at primes outside \(N\) is irrelevant to (6).

Full coverage implies \(N\mid ab\), since the interval length at a source
prime is \(|v_p(a/b)|=v_p(ab)\). The divisibility \(N\mid ab\) alone is
not a sufficient coverage condition: it specifies lengths but not the
locations of the coefficient intervals.

If coverage also holds at every prime outside \(N\), then \(k_j=1\) for
all anchors. In this case
\[
 A_j=\gcd(K_c,|c|N)=N\gcd(K_c,|c|)\le4N.
 \tag{7}
\]
For the equality, \(N\mid ab\) and \((a,b)=1\) imply \((N,c)=1\), since
\(N\) is odd. At odd primes, \(K_c\) and \(c\) are coprime. At two,
the common factor is at most four: if \(\epsilon=1\), then \(a,b\) and
\(K_c\) are odd; if \(\epsilon=2\), then \(a,b\) are odd and
\(v_2K_c=2\); if \(\epsilon=4\), the reduced norms have opposite parity
and \(c=b-a\) is odd. These parity alternatives follow directly from the
Gaussian gcd defining \(\epsilon\).

## 4. Exact elimination of the anchor product

Write
\[
 P=\prod_{i<j}f_{ij},\quad Q=\prod_{i<j}C_{ij},\quad
 B=\prod_{i<j}b_{ij},\quad B'=\prod_{i<j}b_{ij}',
\]
and
\[
 \Gamma_j=\max\{1,D_j^{2m-3}/T_j^2\}.
\]
Applying the source-capped bound at every actual anchor and taking the
geometric mean gives
\[
 N'\ge
 \frac{P^{2/m}(\prod_j\Gamma_j)^{1/m}}
 {2B\tau^m(\prod_jk_j)^{(2m-2)/m}
             (\prod_jA_j)^{(m-2)/m}}.
 \tag{8}
\]
The same-side overlap theorem gives
\[
 (\prod_jk_j)^{m-1}Q\le(B')^2P.
\]
Raising this inequality to power \(2/m\) and substituting into (8) cancels
the complete factor \(P^{2/m}\). The resulting bound is
\[
 \boxed{
 N'\ge
 \frac{Q^{2/m}(\prod_j\Gamma_j)^{1/m}}
 {2B\tau^m(B')^{4/m}(\prod_jA_j)^{(m-2)/m}}.}
 \tag{9}
\]
Unlike the uncapped all-anchor formula, (9) contains no product of anchor
scales \(k_j\). Those scales remain only inside the exact gcd caps \(A_j\)
and the beneficial \(\Gamma_j\). In particular one can bound
\((\prod_jA_j)^{(m-2)/m}\le(|c|N)^{m-2}\).

If both endpoint product estimates are available, write
\[
 B\le L N^{m/8},\qquad B'\le L'(N')^{m/8},
\]
where \(L=(C/2)^{\binom m2}\). For a target with \(r\) odd-odd primitive
columns one may take
\(L'=(C'/2)^{\binom m2}2^{r(m-r)/2}\), without a target unit extraction.
Then (9) becomes
\[
 \boxed{
 (N')^{3/2}\ge
 \frac{Q^{2/m}(\prod_j\Gamma_j)^{1/m}}
 {2L(L')^{4/m}\tau^m N^{m/8}(\prod_jA_j)^{(m-2)/m}}.}
 \tag{10}
\]
All quantities in this inequality are explicit source, target, or normalized
map invariants. No unproved averaging step is used.

## 5. What this closes and what remains

Equation (6) closes the proposed strict-descent obstruction when coefficient
intervals cover the entire source profiles. Equation (5) handles partial
coverage prime by prime, and (9) is an exact stronger way to combine the
capped content theorem with the overlap theorem in the remaining cases.

These results still permit \(N'=N\), and do not control how much source
valuation width can lie beyond the interval endpoints in arbitrary profiles.
Although (9) removes the full anchor product, it retains
\(\tau=\epsilon(a+b)/2\); a bound on \(c=\epsilon(b-a)/4\) alone does
not bound that sum. No deduction of a uniform endpoint count from these
remaining inequalities is asserted.
