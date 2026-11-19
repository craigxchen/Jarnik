# Referee report, block "conseq": Part IV, "Where this construction stops on the integer problem" (assuming (IV.2))

Scripts are mine, in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref4_conseq/`:
* `consts.py`: exact-rational recomputation of every constant in (b) and (c), the optimal sample sizes, the Turán fractions, the hypergeometric second moment checked by exact enumeration, and the separation identity.
* `ffstar.py`: an independent check that a genuine C(t)/Q(i)(t)-cluster (z zbar = N) maps to Part IV's C[t] data, with E_ij = 0 and delta = max deg(z_i - z_j) - D/2. It covers squarefree and non-squarefree parents, and evaluates (IV.2) on the result.
* `remainder.py`: sanity checks of the one-parameter remainder lemma and the Bezout constant, plus the two-parameter failure of the same argument.

Throughout I **assume (IV.2)**, as the block instructs. Its proof (and Langer's Thm 0.1) belongs to another block (`ref4_bmy.md`). I have not verified Langer's statement myself.

## Verdicts in one table

| claim | verdict |
|---|---|
| (a) specialisation lemma (one parameter, nonconstant parent, fixed C1) | CORRECT WITH CORRECTIONS (scope) |
| (a') several parameters / parameter-dependent C1 | GAP (not covered by the proof; repairable only under extra hypotheses, stated below) |
| (b) 1/136 on any 48 marks; Turán 1/47; (IV.4) | CORRECT. 48 is not the optimal sample size (49 is, at 75/10192); the global maximum excess is in fact >= (19/1480 - o(1)) L |
| (c) mean missing >= (1/140 - 1/(8(m-1))) L; 34/35 | CORRECT (needs pair separation >= L/2 for 34/35; without it the bound is 139/140) |
| (d) "no polynomial template gives >= 25 points at fixed C" | CORRECT WITH CORRECTIONS. It holds for rational-branch Q(i)(t)-templates, uniformly in the degree, and by a **more direct and stronger route than Astra's**: it implies M(D,0) <= 24 for all D, which closes a problem the record lists as OPEN. Unit (Pell/golden) templates, multi-parameter families and isolated parents are untouched. |

## (a) The specialisation lemma

**Setting.** Polynomials N, n_i, P_i = (X_i, Y_i) in Z[t] (or integer-valued Q[t]) supply the data of an integral fake at integers T_k, with |T_k| -> infinity. The proof has five steps, and each is correct.

1. **Degrees.** Take logarithms of 0 < |D_ij(T)| <= C1 sqrt(n_i n_j) N^(-1/4) and divide by log|T|. Since |F(T)| = |lc F| |T|^deg F (1+o(1)) for a single variable, this gives deg D_ij <= (alpha_i+alpha_j)/2 - L/4 + o(1). The degrees are fixed quarter-integers, so the o(1) disappears and delta = 0. The same holds if C1 = C1(T) with log C1(T) = o(log|T|).
2. **Primitive normalisation.** First remove a common polynomial factor h of X_i, Y_i. The primitive part of h(T)(X',Y') equals the primitive part of (X',Y'). Write gcd(X'(T),Y'(T)) = c_i(T). It divides the fixed integer Res(X',Y') (Bezout), so it is bounded; `remainder.py` checks this on an example (gcds 1, 2, 4, 41, 82, 164 all divide Res = 164). The integer minor is D'_ij(T)/(c_i c_j), so degrees are unchanged.
3. **Divisibility goes the right way.** Let G = gcd(n_i, n_j) in Z[t], taken primitive. Then G(T) | g_ij(T). Also g_ij(T) | D_int(T), and D_int(T) divides D'_ij(T). So G(T) | D'_ij(T) on the sequence.
4. **Remainder argument.** Write c D' = Q'G + R' in Z[t], with deg R' < deg G. If deg G >= 1, then |R'(T)| < |G(T)| for large |T|, and G(T) | R'(T) forces R'(T) = 0 infinitely often, so R' = 0. (If deg G = 0, G = ±1 and there is nothing to prove.) Divisibility over Q[t] is the same as over C[t], so E_ij = 0. The same argument gives n_i | N.
5. **Robustness.** The argument works on any infinite subsequence (progressions, sparse sets). It also works if the integer data equal the polynomial values only up to bounded factors (for example the fixed gcd g(t0) of `bridge.md` Thm 1.2 on its optimal progression). This is because G(T) | B·F(T) with B fixed still forces G | F.

**Is "nonconstant parent" enough? Yes, for one parameter.** L >= 1 is all that (IV.2) needs. Pair separation, norms, the cone and n_i being nonconstant are all unnecessary. Finitely many sign or unit choices that vary with T are handled by passing to a subsequence. Rational parameters with bounded denominators reduce to integers by T = r + M s.

**Parameter-dependent C1.** Let log C1(T) <= kappa log|T|. Then delta <= kappa, and (IV.2) excludes m marks (E = 0) only when (3/4)(m-24)L > m(m+3)kappa. In particular 25 marks are excluded iff kappa < 3L/2800. Any polynomial growth of C1 that is not small compared with L/m^2 therefore escapes. In the arc problem C is fixed and C1 = C/2 (the star normalisation, `ptolemy.md` Thm B(i)), so this case does not arise there.

**Several parameters: GAP.** The proof uses one variable twice:
* |F(T)| is about |T|^deg F;
* division with remainder.

Both fail in several variables. Take G = T1 and F = T2^2 + 1. Then G(T) | F(T) on the Zariski-dense set {(a, b0 + ka) : b0^2 = -1 mod a} (`remainder.py` lists points), but T1 does not divide T2^2 + 1. Sizes also depend on the direction of approach: log|F(T)|/log|T| is not deg F near the zero set of the top form.

The repair works only under an extra hypothesis. Suppose infinitely many of the parameter points lie on one polynomial curve T = phi(s), along which the parent is nonconstant and the data are given by the composed polynomials. Then the one-parameter lemma applies to the composite. Without such a hypothesis, Part IV says nothing. Astra's wording ("fixed integer polynomials ... at infinitely many unbounded integer parameters") is one-parameter, so the text is correct as written. Any reading that extends it to families in several parameters is unjustified.

## (b) Any 48 marks

* **Excess per sample.** With E = 0, apply (IV.2) to a k-subset, keeping the whole parent, and drop the sum of squares. The maximum edge excess in the subset is then >= 3(k-24)L/(4k(k+3)). At k = 48 this is 72/9792 = **1/136**, as claimed (exact, `consts.py`).
* **Optimal sample size: 48 is not optimal.** The real maximiser is k = 24 + 18 sqrt 2 = 49.46. The best integer is **k = 49**, giving 75/10192 = 1/135.89; k = 50 gives 39/5300. The difference is negligible, and Astra does not claim that 48 is optimal.
* **Turán.** Call a pair good if its excess is < L/136. The good graph is K_48-free, so it has at most (1 - 1/47)m^2/2 edges. The bad pairs therefore number at least m(m-47)/94, a fraction >= (m-47)/(47(m-1)) of all pairs, which is 1/47 - o(1). This is correct. Sample size k trades the excess 3(k-24)/(4k(k+3)) against the fraction 1/(k-1). Their product is maximised at k = 36.
* **A strict improvement Astra missed (global maximum).** Average (IV.2) over all k-subsets and keep the square term. By the hypergeometric second moment from (c), E sum_h (s'_h - k/2)^2 >= L k(m-k)/(4(m-1)). This gives

  max_{i<j} excess >= L[(k-18) - k(k-1)/(4(m-1))]/(k(k+3)).

  This is maximised at k = 37, where it equals **(19/1480 - o(1)) L, about L/77.9**. The improvement is for the maximum over all pairs. It does not improve the "every 48-subset" statement: a balanced 48-subset can have zero square term, so 1/136 (or 75/10192) is the right per-subset constant by this method.
* **(IV.4).** Let W = lambda L + o(W) and log n_i = lambda alpha_i + o(W), uniformly. On a bad pair the genuine bound gives log|D_ij| <= lambda((alpha_i+alpha_j)/2 - L/4) + O(1) + o(W). Since lambda deg D_ij >= that + lambda L/136, the inequality (IV.4) follows. It is correct, but it is only the degree statement restated: every integer instance has a trivial constant presentation (L = 0), to which it does not apply. At a single specialisation the "shortfall" can also come from the bounded-but-large Bezout contents c_i. Astra's own caveat ("the definition supplies neither that accuracy nor a polynomial presentation") is accurate.

## (c) Partial divisibility

Hypotheses: delta = 0 for every pair, and m >= 36. Applying (IV.2) to a uniform k-subset and taking expectations gives

  (3/4)(k-24)L + sum_h E(s'_h - k/2)^2 <= 4 C(k,2) Ebar.

* **The square term.** E(s'-k/2)^2 = k x(1-x)(m-k)/(m-1) + k^2(x-1/2)^2, where x is the layer's global marked fraction. In y = x - 1/2 the coefficient of y^2 is k^2 - k(m-k)/(m-1) = k m(k-1)/(m-1) >= 0. So the minimum is at x = 1/2, with value k(m-k)/(4(m-1)): 9(m-36)/(m-1) at k = 36. I also checked this by exact enumeration at m = 50, k = 36 (minimum at s = 25, value 18/7).
* **The mean missing degree.** The bound is Ebar >= L[(k-18)/(2k(k-1)) - 1/(8(m-1))]. This was verified symbolically, and at k = 36 it is exactly (1/140 - 1/(8(m-1)))L. The constant g(k) = (k-18)/(2k(k-1)) equals 1/140 at both 35 and 36 (and 4/561, 19/2664 at 34, 37), as claimed. Sanity check at m = 36: (IV.2) directly gives sum E >= 9L/4, so Ebar >= L/280, which agrees.
* **The Cauchy–Schwarz step.** Pair separation >= L/2 together with sum_h s_h(m - s_h) = sum_{i<j} sep_ij (randomised identity check) gives sum_h (s_h - m/2)^2 <= mL/4. Cauchy–Schwarz then gives sum_h |s_h - m/2| <= L sqrt(m)/2. Hence the mean prescribed common degree is <= L/4 + L/(2 sqrt m) + O(L/m).
* **The fraction.** The missing fraction is >= (1/140 - O(1/m))/(1/4 + O(1/sqrt m)), which tends to **1/35**. So the retained fraction is <= **34/35 + O(m^(-1/2))**, which is correct. Explicit lower bounds on the missing fraction are 0.0192 at m = 100, 0.0280 at 10^4 and 0.02851 at 10^6, against 1/35 = 0.02857.
* **What the hypotheses do.** Without pair separation the mean shared degree can approach L, and the same argument gives only <= 139/140. The bound transfers to fixed one-parameter families: the integer gcd exceeds the polynomial gcd only by resultant-bounded factors, and C1 fixed gives delta = 0.
* **What it does not say.** Nothing on the method side matches it. III.6's "half via CRT" is an integer construction, not a polynomial family with partial divisibility. For kill criterion (iii), the best attainable fraction is open between whatever a polynomial construction achieves (none was exhibited) and 34/35.

## (d) What is and is not excluded for the integer problem

**A direct route that Astra did not use, and that gives more.** (IV.2) applies **verbatim to genuine function-field clusters**. Let z_1..z_m be monic of degree D in C[t] (or Q(i)[t]), with z zbar = N. Take the cluster primitive with no real roots: a real root divides every z_j. Then:
* choose one orientation f of each conjugate root pair, and put P_i = prod f^{a_i,f} and F = prod f^{e_f};
* then z_i = Fbar P_i / Pbar_i. The phase is 1 because everything is monic, so no drifting-phase issue arises;
* X_i = Re P_i and Y_i = Im P_i are coprime, since the orientation is one-sided;
* n_i = P_i Pbar_i divides N;
* gcd(n_i, n_j) divides D_ij = Im(Pbar_i P_j). The valuations are min(a_i, a_j) at both f and fbar;
* deg D_ij = deg(z_i - z_j) - D + (alpha_i + alpha_j)/2, with L = 2D.

So E_ij = 0 and delta = K - D/2. `ffstar.py` checks all of this exactly over Q(i)[t]: D = 4 (6 points, delta = 0; 16 points, delta = 1), D = 6 (20 points, delta = 1) and a non-squarefree D = 6 parent (10 points). (IV.2) holds in every case.

Hence, assuming (IV.2):
* **M(D,0) <= 24 for every D**. For odd D, delta = -1/2, which is also fine.
* **sup_D M(D,c) < infinity for every fixed c.** For m = 25 the bound fails once 1.5D > 700c, and the record's 2^(2c+2)(D+2) covers D <= 467c.

The record lists a bound uniform in D as **OPEN** ("truncated function-field Vojta", rounds 29–30). So this consequence, not the ones in Part IV's text, is the significant one. The bmy referee reports the same observation independently.

**The corollary in (d) is correct with these hypotheses.** Let the template be:
* a one-parameter family z_j(u) in Q(i)[u] of fixed degree D;
* specialised along integers, or rationals with bounded denominators, tending to infinity;
* with the integer points being the gcd-reduced values w_j = z_j(t0)/g(t0);
* with arc constant C(t0) bounded along the sequence.

By `families.md` Lemma P / `bridge.md` Thm 1.2, bounded C forces K <= D/2. The direct route then gives at most 24 template points, uniformly in D, C and the template. Astra's specialisation lemma gives the same for integral-fake families. Consequences for the record:
* rational-branch templates can never give M_infinity(C) >= 25 for any C;
* bridge §0.6 asks for "unbounded Q(i)(t)-clusters" to get sup_C M_infinity(C) = infinity; under (IV.2) such clusters do not exist;
* Lemma 3.1 still allows the 10-point cluster to be a c = 0 specialisation (10 <= 24).

**Not excluded.**
1. **Unit (Pell/golden) templates.** These hold the record (C = 5.8977 at 8 points). They are function-field clusters at level c = D/4, so delta = L/8. Then (IV.2) is vacuous: m(m+3)/8 >= (3/4)(m-24) for every m, because m^2 - 3m + 144 > 0.
2. **Families in several parameters** whose good parameter set contains infinitely many points of no polynomial curve (see (a)).
3. **Parameter-dependent constants** C1(T) = T^kappa with kappa >= 3L/2800.
4. **Templates producing only part of a cluster.** A template caps only its own points. Extra sporadic points at a given t0 are not controlled.
5. **Isolated parents, the worst-case parent, and every finite-R statement.** The threshold beyond which a template's specialisations have C(t0) close to C_infinity is template-dependent, and nothing here is uniform in R. The growth rate (2/3 + eps) log R/loglog R and the uniform bound for integer circles are unchanged.

## Novelty and significance

* **(a)** is the standard "polynomial divides on infinitely many integers implies polynomial divisibility" lemma. It is correct and routine.
* **(b) and (c)** are routine sampling and Turán corollaries of (IV.2). Their constants are right, and are slightly improvable as noted.
* The orbifold-BMY input itself is new against the record. Its only precursor is round4/L7.md §7, which found slack at seven marks with Kummer weights.
* **On the integer goal the significance of this section is nil.** It excludes only the weakest construction mechanism (rational-branch polynomial templates, per bridge §3, "the weakest mechanism"), and Astra says so.
* The one significant consequence is the function-field corollary M(D,0) <= 24 uniformly in D. That is the first exceptional-set-free uniform bound in any model of the problem. It is not stated in Part IV and it has no arithmetic transfer, so it depends entirely on (IV.2) surviving the bmy block.
