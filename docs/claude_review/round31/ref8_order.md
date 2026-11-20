# Referee report, block "order": Part VIII.2, the height–order inequality (VIII.6)

Scope: `astra_free_VIII/answer.md` lines 1709–1739 (VIII.2 and its function-field test), the statement of (VIII.6) on line 1675, and item 9 of `final8.txt`. VIII.1 and VIII.3–4 are refereed elsewhere. I use them only where (VIII.6) is compared with them.

My scripts are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref8_order/`:
* `check10.py` → `check10_out.txt`: the 10-point cluster, in exact integer and rational arithmetic, with mpmath at 60 digits for the angles.
* `walsh.py` → `walsh_out.txt`: the Walsh-block lattices.

I did not import Astra's scripts.

## 0. Results

**Notation.** The normalised cluster is z_i = ε_i ∏_p π_p^{a_ip} π̄_p^{e_p−a_ip} with i = 1..M, and N' = ∏ p^{e_p}, W = log N'.
* B is the (M−1)×r integer matrix with rows B_i = 2(a_i − a_1), i = 2..M.
* L = row_Z(B) ⊆ Z^r, and s = rank B.
* For v ∈ Z^r, w(v) = Σ|v_p| log p.
* D_L(v) is the least D ≥ 1 with Dv ∈ L. It exists if and only if v ∈ L⊗Q; such v are said to have finite order.
* t_p = 2 arg(π_p)/π, the angle in quarter-turns.

**Theorem (VIII.6), referee-checked. Proved.**
* *Constants.* There are an absolute effective constant c_* > 0 and a threshold M_0(C) depending only on the arc constant C.
* *Statement.* Take any normalised cluster of M ≥ M_0(C) distinct lattice points on an arc of length ≤ C√R of x²+y² = N. Take any nonzero v ∈ Z^r of finite order modulo L. Then (1 + w(v)) log(2 D_L(v)) ≥ c_* W.
* *Uniformity.* Nothing depends on the radius, the position of the arc, the prime set or the profile.
* *Value of the constant.* One can take c_* = 1/(13(2K+3)), where K is the constant of the two-logarithm Baker–Wüstholz bound in §1(6). That value is my computation.

**Corollary (low-height consequence). Proved.**
* *Statement.* Fix A > 0. Let M ≥ M_0(C, A) and let p | N' be a prime with p ≤ M^A whose coordinate vector e_p has finite order modulo L. Then log(2 D_L(e_p)) ≥ (3c_*/A − o_C(1)) M.
* *The window.* In P_+(N') ≤ M^{2−δ} this gives D_L(e_p) ≥ exp(c M) with c = 3c_*/2. The constant can be chosen independent of δ. The window hypothesis is used only for the single prime p, not for the others.

**Corollary, in plain terms.** Suppose an actual cluster has a prime p ≤ M^A whose coordinate is torsion in Z^r/L. Equivalently, D e_p is an integer combination of the exponent differences for some D. Then the smallest such D is exponentially large in M. A second way to say it is that arg π_p lies within exp(−W/12 + o(W)) of a rational multiple (π/2)(k/D) of π/2 whose denominator D is at least exp(cM).

**Proved, but only conditionally, and it gives no growth.** Let f(M) = exp(o(M)) be a fixed function. The size of clusters is bounded in terms of C and f alone among clusters that have one finite-order coordinate at a prime ≤ M^A with order ≤ f(M). No known mechanism produces such coordinates. The trivial upper bound for any finite order is H ≤ exp((1/6 + o(1))W), and the lower bound exp(cM) is far below it. So (VIII.6) is consistent with every growth rate that Theorem 1.1 allows.

**Computed.** The data of the 10-point record cluster are confirmed exactly:
* the Smith invariants of B are 2, 2, 2, 2, 24;
* the coordinate orders are D_L(e_p) = 24, 8, 6, 4, 24 for p = 5, 13, 17, 29, 113;
* the exact phase is t* = (7/24, 3/8, 1/6, 1/4, 11/24), and it satisfies B t* = k on all nine rows.

At M = 10, however, the lemma's own sufficient condition (s+1)H‖u‖_∞ < 1 **fails**: its value is 65.2. The record cluster therefore illustrates the identities. It does not illustrate the asymptotic mechanism.

**Function-field failure. Proved after an arithmetic check.** The quantities in VI.5 are as Astra states:
* w(v) = W/(M−1) for the oriented indicator of one Walsh block;
* D_L(v) = M exactly, computed for M = 4, 8, 16, 32.

So (1 + w) log(2D)/W → 0. Her explanation of the failure (depth amplification under p-th powers) is a correct heuristic. It is not a proof that every two-logarithm estimate fails in characteristic p, and it is not needed: the explicit family already settles the matter.

**Novelty. Low.**
* (VIII.6) is the record's `growth/auxiliary.md` Proposition B2 ("Baker-improved perfect-power certificates"; Matveev in place of Baker–Wüstholz) with the certificate size ‖y‖_1 bounded by Cramer and Hadamard. That file already states the moral: "the general Baker loss is the missing small exponent of Z^s/Λ_B". `growth/auxiliary.md` itself is unrefereed.
* The analytic lemma is identical to the round-31 two-logarithm lemma already refereed in `round31/ref_baker.md`.
* What is new is the explicit uniform packaging: the 1/6 minor bound, the 1/12 margin, the t* correction that removes the dependence on ‖y‖, and the one-line window corollary.

**Significance.** Its only content beyond the record is a structural fact about the pinned and partly pinned case in the window. That case is the record's acknowledged blind spot (round 21, J3). (VIII.6) says that in the pinned case every coordinate order is at least exp(cM). Liouville alone gives only about M/4. It does not bound M, and it does not improve 2/3. **It is a curiosity with modest map value.** It has not been tested against Astra's own window fakes (VIII.3), which is the obvious next test (§3).

## 1. Step-by-step verdicts

### (1) The set-up B t = k + u. CORRECT
* *Angle formula.* arg z_i = arg ε_i + Σ_p (2a_ip − e_p) φ_p with φ_p = arg π_p. Hence arg z_i − arg z_1 ≡ (π/2)(B t)_i mod π/2, because unit ratios contribute multiples of π/2.
* *Size of the error.* Put u_i = (2/π)(θ_i − θ_1), where θ_i are the lifted arguments on the arc. Then B t = k + u with k ∈ Z^{M−1}.
* *Effect of normalisation.* Division by the common gcd g maps the arc to an arc of length ≤ C√R/|g| ≤ C√(R/|g|) on radius R' = √N'. So every angular difference is ≤ C/√R' = C e^{−W/4}, and ‖u‖_∞ ≤ (2C/π)e^{−W/4}.
* *Common phase.* It is eliminated by differencing, not set to zero: θ_1 is arbitrary and never used. The Gaussian gcd is divided out exactly, and its cost is the passage from R to R', which only lowers C. No standing caution is violated.
* *Computed.* On the record cluster: gcd = i (a unit), C' = 8.2393269, and ‖u‖_∞ = 0.0283199. This equals (2C'/π)e^{−W/4} to all digits, because the end pair attains C'. The identity B t − k = u holds to 10^{−61}.

### (2) The minor bound (VIII.7): log H ≤ (1/6 + o_C(1)) W. CORRECT
* *Hadamard.* The entries of column p have size ≤ 2e_p. Hadamard gives |det| ≤ (2√j)^j ∏_{p∈J} e_p for a j-minor.
* *AM–GM.* Every p is ≥ 5, so Σ_J e_p ≤ W/log 5. AM–GM then gives the displayed bound j^{j/2}(2W/(j log 5))^j.
* *First factor.* j ≤ s ≤ M−1, so the first factor is ≤ (M/2) log M. By Theorem 1.1 applied to the normalised cluster (on |z|² = N', arc constant ≤ C), (M/2) log M ≤ (W + O_C(M))/6 = (1/6 + o_C(1))W.
* *Second factor.* j log(A/j) with A = 2W/log 5 increases for j < A/e. Here A/e ≈ 0.46W ≫ M. So the maximum is at most M log(2W/(M log 5)) = W · log(cx)/x with x = W/M ≥ 3 log M − O_C(1), and this is o_C(W).
* *Constants.* The o_C(1) depends on C only.
* *Computed.* At M = 10, log H/W = 0.285 (H = 384). This exceeds 1/6, as expected at small M. The bound is asymptotic.

### (3) The correction t* and B t* = k on all rows. CORRECT
* *Definition.* T = B_{I,J} is nonsingular s×s with d = |det T|, and t* = t − (T^{−1}u_I on J). Then B_I t* = k_I.
* *Other rows.* Every row b is c·B_I with c = b_J T^{−1}, since the rows I span the row space over Q. By Cramer, the entries of dc are, up to sign, s-minors of B obtained by replacing one row of T by b_J. So |(dc)_m| ≤ H.
* *Integrality.* d(b t* − k_b) = (dc)·k_I − d k_b is an integer.
* *Smallness.* b t* − k_b = u_b − c·u_I. So the same integer has absolute value ≤ d|u_b| + sH‖u‖_∞ ≤ (s+1)H‖u‖_∞, using d ≤ H.
* *Conclusion.* (s+1)H‖u‖_∞ ≤ M (2C/π) exp(−(1/12 − o_C(1))W) → 0, because W ≥ 3M log M − O_C(M). So the integer is 0 once M ≥ M_0(C). The threshold depends on C only.
* *Computed.* At M = 10 the left side is 65.2, so the hypothesis is not met. B t* = k nevertheless holds on all nine rows: B has full column rank 5 here, t* = T^{−1}k_I is the unique solution, and I checked the four other rows exactly.

### (4) Finite order, D | d, and the 1/12. CORRECT
* *Finite order and the row space.* v ∈ L⊗Q gives v = (v_J T^{−1})B_I, and d·v_J T^{−1} is integral (Cramer again). So dv ∈ L. The set {D : Dv ∈ L} is the ideal D_L(v)Z, hence D_L(v) | d.
* *Integrality of D v·t*.* Dv = λB with λ integral, so Dv·t* = λ·k ∈ Z.
* *The estimate.* ‖Dv·t‖_{R/Z} ≤ D|v_J T^{−1}u_I| ≤ (D/d) s‖v‖_1 H‖u‖_∞. The adjugate entries are (s−1)-minors, which are ≤ H.
* *The 1/12.* For w(v) ≤ W, ‖v‖_1 ≤ W/log 5. The bound is then exp(log M + log W + (1/6 − 1/4 + o_C(1))W) = exp(−(1/12 − o_C(1))W).
* *Large weight.* For w(v) > W the inequality is trivial, with any c_* ≤ log 2.

### (5) α, its height, non-torsion, α^{2D}. CORRECT
* *Form of α.* α = γ/γ̄ with γ = ∏_{v_p>0} π_p^{v_p} ∏_{v_p<0} π̄_p^{|v_p|}. Here γ and γ̄ are coprime, since π_p and π̄_p are non-associate for split odd p. N(γ) = e^{w(v)}.
* *Height.* The only places with |α|_v > 1 are those dividing γ̄, so h(α) = (1/2) log N(γ̄) = w(v)/2 in the absolute logarithmic normalisation.
* *Non-torsion.* By unique factorisation, α is not a root of unity when v ≠ 0.
* *The power.* arg α = 2Σ v_p φ_p = π v·t, so α^{2D} = exp(2πi D v·t). The factor 2 in t_p = 2 arg π_p/π is right.
* *Non-vanishing.* Dv·t ∉ Z, because α^{2D} ≠ 1.

### (6) The two-logarithm bound. CORRECT, with the source checked only secondarily
* *Statement used (Baker–Wüstholz, Crelle 442 (1993) 19–62, main theorem).* For Λ = b_1 log α_1 + … + b_n log α_n ≠ 0 with α_i in a number field of degree d: log|Λ| > −c(n,d) h'(α_1)⋯h'(α_n) log B. Here:
  * h'(α) = max(h(α), |log α|/d, 1/d);
  * B = max(|b_i|, e);
  * c(n,d) = 18(n+1)! n^{n+1}(32d)^{n+2} log(2nd).
* *Specialisation.* n = 2, d = 2, α_1 = α, α_2 = −1 = e^{iπ}, b_1 = 2D and b_2 = −2j, where j is the nearest integer to D v·t. Then:
  * h'(α) ≤ (π/2)(1 + w(v));
  * h'(−1) = π/2;
  * B ≤ 2D + 2;
  * |α^{2D} − 1| ≥ (2/π)|Λ|, since |Im Λ| ≤ π.
* *Result.* −log|α^n − 1| ≤ K(1 + h(α)) log(2n) for every n ≥ 1, with K absolute. The bound is linear in h and logarithmic in n, as claimed.
* *Constants.* Nothing depends on the number of primes: the whole product is one algebraic number of degree 2.
* *Case j = 0.* This is Liouville, or the same theorem with b_2 = 0.
* *Source.* I could not open the Crelle paper; the network blocks it. I did not read it. I read:
  * the referee-checked statement in `round31/ref_baker.md` §4 (same specialisation; it cites arXiv:1611.04664, which states −log|α^b − 1| ≪ D^4 log D · h'(α) log b);
  * `round8/lw.md` line 126, which gives the same product-of-heights shape.

  A web search returned only snippets (the de Gruyter landing page and the arXiv:1611.04664 entry). The snippets confirm the shape "h'(α) times log b" but not the constants. The constant c(n,d) above is quoted from memory, as in ref_baker.md.
* *Alternatives (from memory, not read).*
  * Laurent–Mignotte–Nesterenko (1995) gives a bound with (log b')² in place of log B. That would yield only (1 + w)(log 2D)² ≥ cW, and so D ≥ exp(c√M) in the window.
  * Gouillon (2006) gives an explicit two-logarithm bound that is linear in log b', and it would serve equally well.

  For the statement as written, BW (or Matveev's degree-two form, as used in `growth/auxiliary.md`) is the right citation.

### (7) Combination, the constant c_*, and the window corollary. CORRECT
* *Combining the two bounds.* |α^{2D} − 1| = 2|sin(π Dv·t)| ≤ 2π‖Dv·t‖. So (1/12 − o_C(1))W ≤ K(1 + w/2) log(4D) + log 2π ≤ (2K + 3)(1 + w) log(2D).
* *The constant.* For M ≥ M_0(C) this gives (VIII.6) with c_* = 1/(13(2K + 3)).
* *Window corollary.* Put v = e_p, w = log p ≤ A log M, and W ≥ 3M log M − O_C(M). Then log(2D) ≥ c_*(3M log M − O_C(M))/(1 + A log M) = (3c_*/A − o_C(1))M.

**Do finite-order coordinates exist?**
* *Full column rank.* If s = r (the pinned case, m = 0), every coordinate has finite order dividing d. Then the exponent of Z^r/L is ≥ exp(cM) in the window. It cannot be larger than the index, which is ≤ H ≤ exp((1/6 + o)W). No contradiction arises.
* *Rank deficit.* If s < r, the generic coordinate has infinite order, and (VIII.6) says nothing about it. The record found no coordinate certificate in 484 tested rank-deficient small clusters (`growth/auxiliary.md`, aux_checks C).
* *Relation to VIII.1.* VIII.1 gives s ≥ M − 1 − K_C, so r ≥ M − 1 − K_C. In the window there are π(M^{2−δ}; 4, 1) ≫ M primes available, so r may far exceed s. Pinned window clusters then require r ≤ M − 1. VIII.1 and (VIII.6) do not combine into any bound on M.
* *Liouville comparison.* The monomial gap for Dv, which angle models satisfy by definition, gives only D w(v) ≳ W/6, that is D ≳ M/4 in the window. So (VIII.6) is the Baker strengthening, from D to log D, of an inequality that fakes already satisfy. Whether window fakes violate it is untested (§3).

### (8) The 10-point cluster. COMPUTED
* *Cluster.* N = 1176852625 = 5³·13²·17·29·113 has 384 lattice points. The best 10-window has C = 8.2393269 (< 33/4).
* *Normalisation.* The common gcd is the unit i, so N' = N. With π_p = a + bi, a > b > 0, the exponent rows are as in `check10_out.txt`.
* *Lattice.* rank B = 5 = r, so the cluster is pinned. The Smith invariants are 2, 2, 2, 2, 24 and the index is 384.
* *Orders.* D_L(e_p) = 24, 8, 6, 4, 24.
* *Exact phase.* t* = (7/24, 3/8, 1/6, 1/4, 11/24). B t* = k holds on all rows, the denominators of t* equal the coordinate orders, and ‖t − t*‖_∞ = 0.0107.

All of Astra's numbers are reproduced.
* At this size, (1 + log p) log(2D)/W ranges over 0.43–1.06. This is consistent with (VIII.6) for any c_* ≤ 0.43.
* ‖D t_p‖ ranges over 0.005–0.084. These values are not small: u is not small at M = 10.

### (9) The function-field failure on VI.5. CORRECT
* *Profile.* Marks are indexed by functionals ℓ on E_0 (M = 2^{d−1}). The block of v ∈ E_0∖0 has column ℓ ↦ ℓ(v) ∈ {0, 1}, the same for all primes of the factor T + θ + v.
* *The relation.* Σ_ℓ (−1)^{ℓ(v)} B_ℓ = −M·1_{block v}, by orthogonality and zero mean. So D_L(v) | M.
* *Exact order.* walsh.py gives D_L = M exactly for M = 4, 8, 16, 32, and rank M − 1.
* *Weight.* A prime π of F_4[t] of degree e lies over a split base prime of degree 2e, and |N π| = 2^{2e}. Hence w(v) = 2L log 2. Since deg N = 2(M−1)L, W = 2(M−1)L log 2, and w(v) = W/(M−1) as stated.
* *Conclusion.* (1 + w) log(2D)/W ≤ (1/W + 1/(M−1)) log(2M) → 0. The failure is genuine.
* *Labels.* The sentence about principal units in characteristic p is a correct heuristic explanation (Heuristic). It is not needed.

### (10) Novelty and significance

| record item | relation to (VIII.6) |
| --- | --- |
| `growth/auxiliary.md` Prop B2 (unrefereed) | Same certificate: Σy_i a_i = Dβ, Matveev two-log, W/4 < 2^40 V log(e(4D+1)) + log(C‖y‖_1). It already says that the Baker loss is "the missing small exponent of Z^s/Λ_B" and that Cramer bounds it by s^{s/2}H^s. (VIII.6) = B2 + Cramer/Hadamard bound on ‖y‖ (or the t* device) + Theorem 1.1. |
| `growth/walsh-extraction.md` Thm B, item 276 | The Hadamard-core instance (D = s/2, β = e_k) is already used there to bound pure Hadamard cores. |
| `round31/ref_baker.md` (Part I §2.2) | The identical two-logarithm lemma, with the same specialisation, refereed CORRECT. There it is applied to the dual lattice of five-tuples. |
| overview §5 rows 10–11; expert note §5 item 4 | Baker/Matveev and single-form LW are capped: the skew double core has w(v) ≥ W/4 for every character, so it satisfies (VIII.6) trivially. Those fakes lie outside the window, though. |
| round 21, J3 and Rem 4.4 | The pinned case is untouched by the measure statements; there, concentration onto (π/2)Z/D_j is forced. (VIII.6) adds D_j ≥ exp(cM) in the window. This is the only new content. |

Verdict on novelty: low. It is a clean uniform restatement of record material, with one new corollary. Verdict on significance:
* it is not strong enough to say anything about growth, alone or combined with VIII.1;
* it gives no restriction at all on unpinned coordinates;
* it is a curiosity, apart from the pinned-case remark.

## 2. Corrections

The statements are correct as written. Two clarifications:
1. **The window corollary needs no window.** It holds for any single prime p ≤ M^A, with constant 3c_*/A. The δ-dependence of "c_δ" is spurious.
2. **The "COMPUTED" sentence needs a caveat.** The 10-point example lies outside the asymptotic range of the lemma: (s+1)H‖u‖ = 65 > 1. The exact phase there exists because B has full column rank, not because of the lemma.

## 3. What would test it

Two tests would decide whether (VIII.6) has any bite.
* *Against Astra's fakes.* Does (VIII.6) separate actual clusters from the VIII.3 window fakes? That happens if some finite-order v in those profiles has (1 + w(v)) log(2D_L(v)) = o(W). Astra did not test this. Because VIII.3 uses fewer than q/2 Paley rows on q columns (rank < r), coordinates are probably of infinite order there. Then the question is the minimal weight of finite-order characters, which is the record's row-11 mechanism.
* *Against pinned window fakes.* Do pinned window fakes exist with polynomial-size coordinate orders? (VIII.6) would exclude them.

Neither test was run here.

Sources:
* [Baker–Wüstholz, Crelle 442 (1993)](https://www.degruyterbrill.com/document/doi/10.1515/crll.1993.442.19/html): landing page only.
* [arXiv:1611.04664](https://arxiv.org/pdf/1611.04664): search snippet only.
