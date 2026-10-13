# Round 9: L4, L5, L1/L2 after refereeing

Every write-up of this round was attacked by an adversarial referee whose report sits next to it
(`referee_*.md`). **§§1-3 give the referee-corrected statements.** §§4-6 add remarks by this README (marked
*README remark*, *sketch* or *Prop. H*). Nobody but the coordinator has checked those remarks, so they are
unrefereed. Conventions and labels are those of `round7/README.md`:
* `N = R^2` and `W = log N`;
* `mu_xy = d_xy - W/2`;
* near-extremal means `W <= c_0 M log M`;
* `X = 3 M^A`;
* the labels are proved / GRH / model / exact / evidence / heuristic / unrefereed.

The question is the growth rate only.

## 0. Bottom line

**No growth improvement for actual circles was proved, either unconditionally or from a standard
conjecture.** `M <= (2/3 + eps) log R/loglog R` (`two_thirds/proof.md`) is still the best bound. **No new
barrier (cap) was proved either.** P3-all\* is still open, and its fixed-`A` bracket
`3 M log M - O(M) <= W_min^all <= (66A + o(1)) M log^2 M/loglog M` is unchanged. All three tasks survive
refereeing as partial results, each with several major corrections.

**What the round established (refereed):**
1. **The regime that matters (L4 §1, proved).** By `stability.md` Cor. 3.6, growth needs P3-all\* only for
   primes `p_j <= M^(2 c_0 - 3 + o(1))` and moduli up to `3 M^A`, with `A` chosen after `c_0`. Every
   construction of rounds 5-7 has `p_j > X`, so it lies outside this regime.
2. **L6 is proved in that regime, unconditionally** (Linnik and Heath-Brown's quadratic large sieve).
   Columns cannot be certified exactly, Hamming-approximately (for almost all columns) or with few pieces by
   Gaussian integers of norm `<= M^h` and sub-genuine height. Exact coherence with `Norm zeta <= M^h` forces
   `p_j Norm zeta` to be a square or twice a square, once `A >= L(h_0 + h) + 2`. The GRH hypothesis of round 7 is gone where it matters.
3. **L4 is undecided.**
   * Construction D (norm `<= 2`) is valid only for `p_j > X`, where it is certified in a generalised
     sense. In regime (1.4) it violates the height condition.
   * The headline "L4 is false up to norm 2; inverse theorems must use norm `>= 3`" is withdrawn.
   * (I1\*) is not a usable inverse step.
4. **L5 is open.** Conditionally on a Hall condition only, it is equivalent to the same statement on the
   small moduli `q < 2cM + 1`. Genuine labels satisfy the other hypothesis (IFL(a)). L5+ is the weaker lemma
   and suffices. Prop. 5.3 and the Hadamard-core and dense-interface claims are corrected or withdrawn.
5. **L1/L2 for actual clusters.**
   * L1 is the statement that all reduced chords are polynomial in `M`.
   * L1 and L2 are mutually exclusive on any near-extremal cluster. They are not each other's negation.
   * Route 1 is an exact split of the growth statement into an archimedean half and an arithmetic half,
     not a reduction.
   * The bisector lemma is new: `M <= sqrt 2 C^2 |g_xy| + 2`.
   * The spectral form of L1 is vacuous by `p3lower.md` Thm A.

**What this README adds (unrefereed; §5):**
* **R1.** Every residue-class cap (rows 2-4 of round 7 §2.1) has `p_j > X`, so its constant grows at least
  like `8A`. Hence these caps block only arguments that fix `A` before `c_0`. No `A`-uniform cap of any
  residue class is proved.
* **Sketch R2.** For `A_X` and P3-char such caps very probably exist, through small auxiliary primes with
  repair, or injective classes at `q >= 8M` together with dense character budgets. For P3-all\* nothing
  similar is available.
* **Prop. H (elementary, with an exact check).** Fix `K`. Take genuine labels below
  `Q_K = 8K(2r+1)^K + 1` and greedy parity-respecting labels above it. The result meets every
  norm-`<= K` height condition. So bounded norms force nothing above `M^(K + o(1))`. An inverse theorem
  that feeds Linnik's Thm A must therefore use norms `K >~ L(h_0 + h)` at moduli up to `M^(K + o(1))`.

**Next step (§6).** Decide the scale-dependent inverse step L4(K) in regime (1.4) at the threshold that
Prop. H fixes. It is the only purely combinatorial open step left in route 3. A proof reduces route 3 to
L7 (Poisson tails of genuine-type data). A counterexample family is the first step to capping P3-all\*
where growth needs it. A finite exact first test is specified.

## 1. L4 (the inverse step (I1)): referee-corrected statement

Write-up: `L4.md`. Report: `referee_L4.md`, verdict "the proved theorems survive; two of the headline
consequences do not". Scripts: `L4_checks/` (author) and `referee_L4_checks/` (referee, independent code).

**Status: L4 is undecided.** Nothing here bounds actual circles beyond `(2/3 + eps) log R/loglog R`.
Every statement is about label systems (designs) in P3-all\*. Actual clusters satisfy them with
`zeta_j = pi_j`, the genuine Gaussian primes.

**Setting.** A label system is `ell_j(q) in Z/m_Q` with the forced parity
`ell_j(q) = [(p_j/q) = -1] mod 2` (`sharp.md` Cor. 2.2). It is *coherent* with Gaussian integers `zeta_j`
if `ell_j(q) = log(zeta_j/conj zeta_j)`, up to `mu_4` and an odd unit at each `q`. `Q_7` is the set of
primes `q = 7 mod 8`. Such `q` never divide `N`, and at them the parity is a Legendre symbol computable by
reciprocity. The **regime (1.4)** is `p_j <= M^(h_0)` with `h_0 = 2 c_0 - 3 + o(1)` and `X = 3 M^A`, where
`A` may be chosen as large as needed in terms of `c_0`.

1. **Reduction (T') (proved; actual circles).** By `round8/stability.md` Cor. 3.6(1), growth follows from
   the following statement:

   > **(T')** For every `c_0` there is `A(c_0)` such that, for large `M`, no member of
   > `F^all_(3M^A)(M, C)` with all `p_j <= M^(2 c_0 - 3 + o(1))` has `W <= c_0 M log M`.

   This is round 7's (T) with the quantifier order `forall c_0 exists A`. In this regime `X` exceeds any
   fixed power of the `p_j` and of `Norm zeta` for `Norm zeta <= M^h` with `h` fixed.
   **All earlier constructions (rounds 5-7) have `p_j > X`, so they lie outside regime (1.4).** §5 below
   draws the consequence for the barrier map.

2. **Exact coherence (Theorem A and Cor. A.1; proved, given Linnik's theorem on the least prime in a
   progression).** Work in regime (1.4) with `A >= 2 L h_0 + 2`, where `L` is Linnik's exponent. Suppose
   a column agrees with `u_q log(zeta/conj zeta)` modulo `mu_4` at every `q in Q_7`, `q <= X`, where `u_q`
   is a unit that may depend on the column.
   * Then `|zeta_j|^2 >= p_j`.
   * If `Norm zeta <= M^h` and `A >= L(h_0 + h) + 2`, then `v_(p_j)(Norm zeta)` is odd.

   So coherent certificates have zero slack. This is the unconditional half of L6 for exact coherence,
   valid in the regime growth needs. It replaces the GRH Theorem 3.2 of `p3Dlower.md` there.

3. **Hamming coherence (Theorem B; proved, given Heath-Brown's quadratic large sieve).** For every
   parity-respecting system, all but `O(M^0.45 log^2 M)` columns have the following property. For every
   `zeta` with `Norm zeta <= M^h` and `v_(p_j)(Norm zeta)` even, the column agrees in parity with `zeta` on
   at most `(1 + delta)/2` of the `Q_7` log-mass.
   * The norm bound is necessary. A Gaussian prime over `p' = p mod 4 prod_(q <= X, q in Q_7) q` is
     sub-genuine and parity-exact on all of `Q_7 ∩ [1, X]` (referee, explicit for `X = 200, 500, 1000`).
   * The bold wording of `L4.md` §0/§1 ("no column is certified by Gaussian integers with `p_j` not
     dividing `Norm zeta_j`") omits this bound and is false without it.
   * *Genuine height* should mean `v_(p_j)(Norm zeta)` odd (or `|zeta|^2 >= p_j`), not `p_j | Norm zeta`.
   * The proof's own accounting needs `delta X >= 60 log D`, which is harmless. The GRH remark needs
     `X >> (log D)^4/delta^2`.

4. **Piecewise coherence (Theorems C, C'; proved, given Linnik, resp. GRH).** Covering one column by `P`
   sub-genuine coherent pieces needs `P > (A - 1)/(L(h_0 + h))` unconditionally, and
   `P >= (1/2 - o(1)) log_2 X` under GRH.
   * The duality Lemma 4.2 (relations among the characters have even size and sign product `+1`) was
     checked exactly on 3000 `F_2` instances.
   * The twisted-coherent labels of `p3Dlower.md` Def. 4.1 have norms `> X`. For them the unconditional
     Theorem C is vacuous, and only the GRH Theorem C' applies.

5. **Proposition D (proved, any regime).** There is an explicit parity-respecting Construction D with
   `m_v | Im(Z_v^16)` for every monomial with `||v||_1 <= 2`. Its ingredients are unit shifts at
   `q = +-3 mod 8`, matched `zeta_j` at `{17, 31, 47}`, a direct parity factor at `b_q = 3`, and an odd
   injective orbit encoding (Lemma 5.1) at `b_q >= 4`. The two regimes behave differently.
   * *Construction regime* (`p_j >= M^(16h+3) > X`). D meets the height condition up to norm 2 with
     slack `(3/2) log M` per unit of norm, while independent labels already fail at norm 2 (Cor. D.1).
     But here D **is certified** in the generalised sense, by the 16th-power Gaussian certificate of
     height `M^h`, which is round 7's plain-language reading of L4. It is not shown to be uncertified in
     any other sense without GRH, since Theorems A-C are vacuous when `p_j > X`. **Cor. D.2's clause "or
     in the setting of Corollary D.1" is wrong.**
   * *Regime (1.4).* D is uncertified in the sense of items 2-4, but **it violates the height condition
     already at norm `<= 2`**. On exact instances (`X = 2*10^4` and `5*10^4`, primes `<= 400`), 93-134 of
     100-144 norm-2 monomials and 8-10 of 10-12 norm-1 monomials violate it, against 0 for the genuine
     labels on the same profile.
   * **Withdrawn** as general claims: "L4 as stated is false up to norm 2" and "any inverse theorem must
     use norm `>= 3` and only the moduli `q = +-1 mod 2^B`". They hold only with the narrow notion
     "certified = F1-coherent on most of `Q_7`", only in the construction regime, and only for inverse
     theorems that must also cover `p_j > X`. **In regime (1.4) both questions are open.**

6. **(I1\*) (formal reformulation only).** (I1\*) asks whether bounded demand forces Hamming-closeness
   to genuine-height certificates. It is (I1) in the per-column Hamming sense, combined with items 2-3.
   It is **not** a usable inverse step, for two reasons.
   * The threshold `1/2 + delta` is mis-calibrated. Against the genuine column, a uniformly random
     parity-respecting column agrees on 0.540-0.555 of the `Q_7` log-mass in expectation (exact per
     modulus, 40 columns, `X = 2*10^5`). Construction D agrees on 0.56 to 0.70.
   * Per-column agreement does not transfer the collision sets of multi-column monomials, which is what
     (I3)/L7 needs.

   A usable version needs system-level agreement (one unit per modulus) on clearly more than half of
   `Q_7`.

7. *Minor corrections.*
   * The density of `{q = +-1 mod 2^(B+1)}` is `2^(1-B)`.
   * The independence in §6.3 holds modulo squares and 2.
   * The additive-label divisor count of §5 and the whole of §6.3 are heuristics.
   * Lemma 5.1 also holds at `N = 12`.
   * Linnik's exponent and Heath-Brown's large sieve are cited from memory. They affect only how large
     `A` must be.

The bracket `3 M log M - O(M) <= W_min^all <= (66 A + o(1)) M log^2 M/loglog M` is unchanged. There is no
growth improvement and no new cap.

## 2. L5 (twin channel): referee-corrected statement

Write-up: `L5.md`. Report: `referee_L5.md`, verdict "survives as a partial reduction, with six major
corrections and no fatal error". Scripts: `L5_checks/` and `referee_L5_checks/`.

**Status: L5 is open.** No construction, no obstruction and no growth statement is proved.

**Setting** (`L5.md` §1). The profile is a twin profile with `e = 1`: each row `x` owns two columns
`f(x), s(x)` that differ only in row `x`. The other variables are:
* the free labels `f = (ell_J, U)`;
* the twin labels `Lambda_x = eps_x ell_(f(x))`;
* the point classes `kappa_y = <phi_y, f> + Lambda_y`.

Demands are `log q`-weighted: `delta_s(v) = log m_(v,s)`. L5 is the following statement.

> **L5.** For all `A, beta, K, rho` and large `M`, no twin data with `r <= rho M` satisfy all three of
> * **(F)** `delta_s(w, 0) <= beta ||w||_1 log M` for the free monomials;
> * **(T)** `delta_s(w, eps e_y) <= beta (||w||_1 + 1) log M` for the twin monomials;
> * **(P)** unit-1 pair demand `delta(x, y) <= K log M` for all pairs.

Write `Q_large(c) = {q : n_q >= cM}` and `Q_med(c) = {q : n_q < cM}`, where `n_q = m_q/2`. Every
`q in Q_med(c)` is `< 2cM + 1`.

1. **Formulation (proved; actual circles).**
   * (F), (T) and (P) hold for near-extremal actual twin-profile clusters with margins `<= K' log M`, with
     `beta = 4/eps + o(1)` and `K = K'/2 + O(1)`. So L5 would give growth for that subclass.
   * Only unit-1 collisions matter in (P). The other units have a budget of order `W` (Remark 1.1). This
     clarifies the round-7 wording rather than correcting it.
   * The unweighted wording ("all but `O(||w||_1 + 1)` moduli") is **not known to be valid for actual
     circles**. Genuine labels can collide at about `log p/(2 loglog p)` small moduli: there are exact
     examples with 12 moduli at `log p = 70.6`.
   * **Withdrawn:** "the unweighted wording is false for every system because every monomial collides at
     `q = 3, 5`". Two moduli are absorbed by `O(||w||_1 + 1)`.

2. **Identities (proved; actual clusters).** Props. 2.1-2.3, Cor. 2.4 and Lemma 2.5.
   * (T) says that each `kappa_y` avoids a dense translate `<phi_y, f>` of the sparse sumset of the free
     labels.
   * *Class-mate transfer:* where `kappa_x = kappa_y` at `q^a`, the twin monomials of `y` collide exactly
     when those of `x`, translated by the dense vector `a_x - a_y`, collide.
   * Cor. 2.4 needs no weight-ratio hypothesis in the near-extremal regime.

3. **Per-modulus freedom (Lemma 3.1; proved).** At any modulus with `n_q > M - 1 + 4 max_y |V_y|` (in
   effect `Q_large`), any prescribed twin monomials can be avoided with `kappa` injective.
   * That a single modulus cannot prove L5 is trivial, since one modulus adds at most `log X` to any
     demand.
   * **Heuristic only:** "no counting argument proves L5 on `Q_med`".

4. **Localisation (Thm 3.2, Cor. 3.3; proved, conditional on the Hall condition only).**
   * IFL(a) is **not** the hypothetical object of L4. Genuine labels of any distinct split primes
     `p_j <= M^(O(1))` satisfy it, by the height bound. So does any multiplicatively independent coherent
     design of height `O(log M)`.
   * Assume the Hall condition IFL(b) at every `q in Q_large(c)`. Then switching each twin label among `t`
     candidates gives data on `Q_large(c)` with zero pair demand, (F) with `beta_0`, and (T) with
     `t beta_0 (1 + 2 c_0)`.
   * Hence, along profile sequences with `M`-independent constants, **L5 <=> L5_med(c)**: the same
     statement restricted to `q < 2cM + 1`.
   * *Status of Hall.*
     * For genuine candidates it is open. Evidence: it holds at all 194 moduli tested (`M = 12`, `t = 5`,
       `c = 2`).
     * In the random-candidate model it holds simultaneously at all `q in Q_large` with probability
       `1 - O_t(M^(4-2t) log(A+1))`. This uses the sharper count; the per-modulus Lemma 3.4 bound
       `1.5/n_q` does not sum to `< 1`.
   * *Correction:* pair collisions are forced only where some parity class has more than `n_q` rows,
     i.e. `n_q < M/2` for balanced classes. Nothing is forced on `M/2 <= n_q < cM`.

5. **L5+ (proved implication).**
   * L5 implies L5+ (all `||c'||_1 <= k`). So L5+ is the **weaker** lemma, not the stronger one, and it is
     what the route needs.
   * Its localisation to `n_q < M^(2+eps)` needs a common-index hypothesis. It also needs the ideal bound
     on the sparse elements of `L`; Cor. 3.6 has a gap when the switched monomial lies in `L`.

6. **Certification barrier (Thm 4.1, Cor. 4.2; proved; construction method only).** Take coherent designs
   with multiplicatively independent columns, and suppose (F) and (T) are certified by height. Then the
   pair heights have geometric mean `>= M^(1/2)/(C log M)^(b-1)` on the Paley twin profile.
   * This extends `p3Dconstruct.md` Thm 2 beyond the Paley profile, in a weaker form.
   * It beats `K log M` only for `M >~ e^450` at `b = 33`.

7. **Non-twin profiles (proved for rank-`M` profiles only).**
   * Thm 5.1: the profile decouples wherever `gcd(Delta_M(A), m_Q/2) = 1`.
   * **Prop. 5.3 is false as written.** It makes `kappa'` injective instead of `kappa` (0/1 counterexample)
     and uses the top level instead of level 1 (`q = 3` counterexample). The corrected level-1 statement
     holds by a greedy proof: if `n_q/gcd(d_M, n_q) >= M`, a collision-free parity-compatible `kappa`
     exists.
   * The Hadamard-core examples `(J + H)/2` are not profiles, since they have a constant column. The valid
     Hadamard-core profiles have rank `M - 1` and are not covered. `d_M = M/2` is a theorem about the
     non-profile.
   * **Withdrawn:** "generic profiles have dense interfaces". Exact or constructive integral sections of
     random 0/1 profiles have row norm 2-3 (`M = 8`, MILP-optimal), `<= 7` (`M = 12`) and `<= 11`
     (`M = 16`) (evidence). So at these sizes the twin localisation transfers to generic profiles.
     Bounded section norm as `M -> infinity` is open.

## 3. L1/L2 (uniform margins; heavy pair): referee-corrected statement

Write-up: `L1L2.md`. Report: `referee_L1L2.md`, verdict "the mathematics survives; the headline does not".
Scripts: `L1L2_checks/` and `referee_L1L2_checks/`. The referee re-derived every proof and re-ran the scan
by an independent method, obtaining exactly the same 22,907 clusters.

**Status: partial.** There is no growth improvement and no new barrier. Every item concerns normalised
actual clusters unless labelled otherwise. L1 and L2 are stated for actual clusters, with the reduced
chord `beta_xy = (z_x - z_y)/gcd(z_x, z_y)` and the normalised gap `t_xy in (0, C]`.

1. **Reduced chord (Thm 1.1; proved).**
   * `beta_xy = i^(k_x) gamma_s B_xy`, with `gamma_s in {2i, 1-i, 2, 1+i}` and `B_xy = D_s(A_(a_x - a_y))`
     a nonzero **rational** integer.
   * For odd `q` not dividing `N`, `z_x = z_y mod q^a` iff `q^a | B_xy`.
   * Exactly, `mu_xy = 2 log|B_xy| + log nu_s - 2 log t_xy`.
   * Window collisions are agreements of the discrete-log code `w_x(q) = <a_x, ell(q)> + k_x u(q)` built
     from the genuine labels (Prop. 3.1).
   * This is the pair case, with equality, of `sharp.md` Prop. 1.2 and the `p3lower.md` Def. 1.1 bullets.
     *Notation clash:* `p3lower.md`'s `B_xy` is the Gaussian integer `zhat_x/g_xy`.

2. **L1 is a polynomial-size statement (Thm 2.1; proved).** For `M >= M_0(C)`, L1 on near-extremal
   clusters is equivalent to the following: for some fixed `b`, every `|B_xy| <= M^b` and every
   `t_xy >= M^(-b)`. Equivalently, every pair gcd has norm `>= N^(1/2) M^(-O(1))`. This is a one-line
   consequence of the margin identity and `mubar = Theta(log M)`.
   * At `M = 4` there are 26 near-extremal clusters with mean margin 0. On them L1 fails for every `K`,
     which is why `M >= M_0` is needed.

3. **L1 and L2 cannot coexist (Thm 3.2; proved).** On a near-extremal cluster with `K`-uniform margins,
   every `B_xy` has fewer than `K c_0/(4(1-eta)) + o(1)` prime factors above `M^(1-eta)`. So
   `K`-uniformity and an L2 pair never occur on the same near-extremal cluster with `M >= M_1`.
   * This is the contrapositive of the proof of Theorem C.
   * **Withdrawn:** "L1 and L2 are each other's negation" and "they cannot be pursued as independent
     lemmas". As universal lemmas both are true, vacuously, if growth holds at `c_0`. Their conjunction is
     equivalent to there being no near-extremal clusters for large `M`.

4. **Heavy pair without the window (Thm 3.3; proved).** `K`-uniform margins, together with one pair with
   `|B|/t >= M^f`, give `W >= (2(M-1)/K)(2 f log M + log 2)`. This is Prop. B.1 of `p3lower.md` with
   `2 log(|B|/t)` in place of the pair demand. The heavy pair is assumed, not derived; "without `H`" is
   misleading.

5. **Route 1 for actual clusters (Cor. 3.4; proved).** Let `(P_b)` say that every `|B_xy| <= M^b` and
   every `t_xy >= M^(-b)`. Route 1 splits into two contrary statements:
   * **(R1-arch):** some `b` has `(P_b)` on every near-extremal cluster;
   * **(R1-arith):** for every `b`, every near-extremal cluster with `M >= M_1(b)` violates `(P_b)`.

   Their conjunction is the growth statement at `c_0`. The same is true of any route phrased about actual
   near-extremal clusters. Each half can be pursued separately, but neither can be tested on data.
   * The actual-cluster reformulation drops round 7's intent, which was to derive L1 from valid inputs in a
     larger class (`A_X`, P3-all\*).

6. **Bisector lemma (Lemma 1.3; proved, new).** `M <= sqrt 2 C^2 |g_xy| + 2` for every pair. Hence
   `mu_xy <= W/2 - 2 log M + O(log C)`, and complementary pairs need `M <= sqrt 2 C^2 + 2`.
   * It is the only known **`M`-dependent** upper bound on an individual margin. The triple identity
     already gives the `M`-free bound `mu_xy <= W/2 - 2 log 2 + 4 log C`.
   * It is an integrality argument giving an upper bound. So the claim "integrality gives only lower
     bounds" must be restricted to divisibility of the reduced chords.

7. **Spectral form of L1 (Prop. 2.2, corrected).** `max_x D_x <= (1/2) sqrt(c_0(c_0 - 2)) M^(3/2) log M (1 + o(1))`
   on near-extremal clusters.
   * For `lambda_max(D)` this is **superseded**: `p3lower.md` Thm A, which is refereed and valid on `A_X`,
     gives `lambda_max(D) <= (c_0/2) M log M + (M - 1) c_C^+`.
   * So the alternative form of L1 ("`lambda_max(D) >> M log M`") also holds only vacuously. Open point
     6.3 of `L1L2.md` is answered.

8. **Other proved items.**
   * Markov-type uniformity is automatic and useless (Prop. 2.3(a,b)).
   * The Poisson rarity of L2 pairs is a model statement (Prop. 2.3(c)).
   * Exact relations `sum c_x z_x = 0` of height `M^(c_0/2 + 1/2 + o(1))` exist, and their content is
     lattice-ness (Prop. 5.2).
   * The cocycle, area and 4-point identities reduce to profile + `B_xy` + residue partitions (Props. 5.1,
     5.3). Prop. 5.1(a) must keep the qualifier "at a single odd `q`".

9. **Heuristic.** Genuineness of the labels alone does not force L2 for rows adapted to the labels.
   Whether genuineness together with near-Plotkin profile structure forces L2 is open. The bold "L2 must
   use the short-arc condition" overstates even the heuristic.

10. **Evidence (exact, `M <= 10`).**
    * Thm 1.1 and Lemma 1.3 hold without exception on 22,907 clusters (`N <= 10^9`, `C <= 4`) and on the
      7-, 8-, 9- and 10-point records.
    * Small actual clusters are non-uniform even without complementary pairs. Among clusters with no
      complementary pair, the median `K` is 2.67 (`M = 4`) and 3.54 (`M = 5`). On the records `K` is
      6.6-13.7.
    * In the golden family, 24 of 28 margins stay bounded while 4 grow with `W`. This is the
      adapted-margin pattern, not uniformity.
    * **Withdrawn:** the figures 1.17/3.28 (an artefact of the statistic), "every margin grows linearly in
      `W`", and the reading "non-uniformity is a symmetry effect confined to bounded `M`".

## 4. Updated status of L1-L8

The lemmas are those of `round7/README.md` §3. "Regime (1.4)" means `p_j <= M^(h_0)` and `X = 3 M^A` with
`A` large in terms of `c_0`; by (T') (§1 item 1) this is the only regime growth needs. **Apart from L6 in
regime (1.4), no lemma was proved this round, and none was refuted in the regime that matters.**

| Lemma | Status after round 9 | What changed |
|---|---|---|
| **L1** (uniform margins) | **open**. For actual clusters it is (R1-arch): every `|B_xy|` and `1/t_xy` is `M^(O(1))` (Thm 2.1). It is true vacuously if growth holds, and not derived from any valid input. | Exact polynomial-size form (proved). The spectral form `lambda_max(D) >> M log M` is vacuous by `p3lower.md` Thm A. The only `M`-dependent individual-margin bound is the bisector Lemma 1.3, far from L1. Small actual clusters are non-uniform (evidence). Round 7's intent (L1 from a valid input in a larger class) is untouched, and Thm A shows that the pair inputs allow adapted margins. |
| **L2** (`H_1`, heavy pair) | **open**. On `K`-uniform near-extremal clusters it cannot hold (Thm 3.2), because L1 and L2 together contradict near-extremality. | L2 is one mechanism for (R1-arith), "some `|B_xy|/t_xy` is super-polynomial". Thm 3.3 shows that any super-polynomial pair suffices once L1 holds. Route 1 = (R1-arch) + (R1-arith), whose conjunction is exactly the growth statement. Route 1 is therefore a split of the problem, not a reduction of it. |
| **L3** (prime floor, (Q)) | **open**; not attempted | none |
| **L4** ((I1), inverse step) | **undecided**. The certified side is settled in regime (1.4): no sub-genuine certification, exact (Thm A), Hamming for almost all columns (Thm B, with the norm bound) or piecewise with few pieces (Thm C). The inverse side is open in regime (1.4). By README Prop. H (§5 item 4), norms `<= K` force nothing above `Q_K = M^(K + o(1))`, so an inverse theorem feeding Thm A must use norms `K >~ L(h_0 + h)`. | Construction D meets the height condition up to norm 2 only in the construction regime `p_j > X`, where it is certified in the generalised 16th-power sense. In regime (1.4) it fails at norms 1 and 2. "L4 false up to norm 2" and "inverse theorems must use norm `>= 3`" are withdrawn for regime (1.4); Prop. H gives the correct regime-(1.4) replacement. (I1\*) is a formal reformulation, mis-calibrated (random baseline 0.54-0.55 > 1/2), and not a usable inverse step. |
| **L5** (twin channel) | **open**. Conditional on the Hall condition only, it is equivalent to L5_med: the same statement on `q < 2cM + 1` (Thm 3.2, Cor. 3.3). | IFL(a) holds for genuine labels, so it is not a hypothesis. Hall holds at 194/194 genuine moduli (evidence) and in the random-candidate model. L5 implies L5+ (the weaker lemma, which suffices). Per-modulus freedom holds on `Q_large` only. Certification barrier on all profiles (Thm 4.1), with an astronomical threshold. Non-twin rank-`M` profiles decouple where `gcd(Delta_M, n_Q) = 1`. Prop. 5.3 is corrected to level 1. Rank-deficient profiles are not covered. |
| **L6** ((I2), certified => genuine size) | **proved in regime (1.4)** for exact coherence (Thm A, Cor. A.1, given Linnik; `A >= 2 L h_0 + 2`). Also proved there for Hamming coherence on almost all columns and for piecewise coherence with `P <= (A-1)/(L(h_0+h))` pieces, both with the norm bound `Norm zeta <= M^h`. | The GRH hypothesis of round 7 is removed in the regime growth needs. Twisted coherence (`Norm zeta > X`) is out of scope there, and its genuine size is trivial (`|zeta|^2 > X > p_j`). The open part of round 7's L6 (approximate certification with `Norm zeta << p_j`) lives only in the construction regime `p_j > X`, which growth does not need. What remains is the interface with L4: a **system-level** agreement notion (one unit per modulus, clearly above the random baseline) under which L4's conclusion and L6's hypothesis match. |
| **L7** ((I3), Poisson tails) | **open**, even for genuine primes; not attempted | It is now the common bottleneck. L5_med needs an L7/`H`-type input on dense translates of the free system. (R1-arith) is an `H`-type smooth-part statement about the `B_xy`. In regime (1.4) the coherent side of L4 is genuine-type labels (Thm A), whose collision statistics are exactly L7. |
| **L8** (slab covering) | **open**; not attempted | none. The README remark of round 7 stands: a covering lemma cannot be a threshold on the sum of slab measures. |

Construction-side items C1-C3 and route-4 items L9-L11 were not attempted. C2 (minimise `D_max` on the skew
double core) is still the undone round-8 "decisive first test".

**Net effect on the routes.**
* *Route 1* (L1 + L2) is, for actual clusters, an exact split of the growth statement (Cor. 3.4). It is not a
  weaker target.
* *Route 3* (L4 + L6 + L7, with L5 on twin profiles). In regime (1.4), L6 is done for exact, Hamming and
  few-piece coherence. The small-height coherent alternative of the L4 dichotomy is genuine-height labels:
  if `Norm zeta_j <= M^h`, then `p_j Norm zeta_j` is a square or twice a square (Thm A). That alternative is the territory of
  route 4 (P3-G) and of L7. So routes 3 and 4 meet in regime (1.4). The open steps are the inverse step L4
  there, and L7.
* *Twin channel.* Given Hall, everything in L5 happens at the small moduli `q < 2cM + 1`.

## 5. Status of P3-all\* with respect to the growth rate

**P3-all\* is still open** (round 7 §2.1 row 5). It is not capped, and no growth is proved for it. At fixed
`X = 3 M^A` the bracket is unchanged:

```text
3 M log M - O(M)  <=  W_min^all  <=  (66 A + o(1)) M log^2 M / loglog M.
```

Its upper end is `sharp.md` Thm 6.4, whose members have all `p_j > X`. Round 9 adds the following.

1. **The regime that matters (proved; (T') of §1).** Growth needs P3-all\* only in regime (1.4): primes
   `p_j <= M^(2 c_0 - 3 + o(1))`, and moduli up to `3 M^A` with `A` chosen **after** `c_0`. There `X` exceeds
   every fixed power of the `p_j`, so Legendre-symbol information about the `p_j` is available
   unconditionally. Linnik and Heath-Brown give:
   * no exact, Hamming or few-piece certification of a column by sub-genuine Gaussian integers of norm
     `<= M^h`;
   * every column exactly coherent with some `zeta_j` of norm `<= M^h` (and `A >= L(h_0 + h) + 2`) has
     `p_j Norm zeta_j` equal to a square or twice a square.

   **So in regime (1.4) a small-height coherent design (D3 of round 7 §2.2, `Norm zeta_j <= M^h`) is
   either impossible or of genuine type.** Round 7 had this only at the certificate level and under GRH.
   Coherence with `Norm zeta_j > X` (twisted designs, D4) is not excluded; there `|zeta_j|^2 > X > p_j`
   holds trivially. Such designs have collision integers of size `>= X^(||v||_1/2)` against budgets
   `<= M^(h_0 ||v||_1/2)`, so their certificates are useless once `A > h_0` (heuristic for their actual
   demands).

2. **README remark R1 (elementary, unrefereed): the residue caps are not uniform in `A`.** Round 7 §2.1
   rows 2, 3, 4 and 4a list caps for `A_X` and P3-char. Every one of them is built with all `p_j > X`:
   * `sharp.md` Def. 3.2 and Thm 5.1 (`r = b(M-1)`, `b >= 33`);
   * `construct.md` Thm 2 (`r = 9M - 8`, `p_j > 160 M^2 L/C` with `L > Q_0 = X`);
   * `lower.md` Thm 5.1, whose auxiliary primes `ell_x` and `p_j` both exceed `X`.

   By `p3lower.md` Lemma D.0, `W >= r log X` when all `p_j > X`. So these caps have
   `W/(M log M) >= 8A (1 - o(1))` at `X = M^A`, matching their quoted constants `9(2 + max(K,1))`,
   `16 + 32A` and `b max(A, 2)`.
   * Each cap is correct for its fixed `A`. The dictionary entry "capped => no argument using only these
     inputs proves growth" holds for arguments that fix `A` before `c_0`.
   * A proof by (T') chooses `A = A(c_0)`. It is blocked only by a cap whose constant stays bounded as
     `A -> infinity`, with `p_j <= M^(O(1))`.
   * **No such `A`-uniform cap is known for any residue class.** That covers row 2 (`A_X`), row 3,
     row 4 (P3-char), row 5 (P3-all\*) and row 6 (P3-G).
   * Row 1 (`A_0`: profile, chords, cut, Liouville gaps, rational characters) involves no `X`, and its
     skew-double-core cap stands.
   * The growth-status sentence of `docs/claude_review/README.md` ("every local input ... residues at the
     point level, and prime-level residues used through characters" is capped) is therefore proved only at
     each fixed polynomial level `X = M^A`.

3. **README sketch R2 (heuristic, unrefereed): rows 2-4 probably extend uniformly in `A`; P3-all\*
   does not.**
   * *Where the `A`-dependence comes from.*
     * In `lower.md` it comes from Lemma 4.1. There the auxiliary primes satisfy `ell_x > X`, so that
       `varpi_x` is a unit modulo every `q <= X`, and the bound is `m_c^gen <= Y^(2||c||_1)` with `Y > X`.
     * In `sharp.md` it comes from the choice `p_j > X`, which makes `N` prime to every modulus.
   * *Mechanism (a), `A_X` only.* Take distinct `ell_x in (M + 1, M^(O(1))]` and repair at the single bad
     modulus `q = ell_x` of point `x`. There choose `rho_x(q)`:
     * distinct from every `rho_y(q)`, so no pair collides there;
     * generic otherwise, at a union-bound cost of `O(||c||_1 log M)` per character.

     The demands become independent of `X`.
   * *Mechanism (b), P3-char and hence `A_X`.* Keep the parity-class design `D_A` at `q < 8M`. At
     `q >= 8M`, make the point classes injective modulo `mu_4` inside each parity class, which is possible
     since `m_q/8 >= M`, and otherwise independent.
     * Pairs then have no demand at `q >= 8M`.
     * Every non-pair character on the `sharp.md` profile has slack `s_c >= (tau/2) kappa_b M - 1/2`
       (Prop. 4.3), i.e. a budget `Theta(M log M)`. Independent collisions at `q <= M^A` add `log X` on
       average.
     * A union bound over the `<= (2t+1)^M` characters of amplitude `t` closes for each fixed `A` once
       `M >= M_0(A)`. Exceeding the slack needs `j` collisions with `j` of order `t kappa_b M/A`, at
       probability `~ e^(-j log j)`.
     * The small-coherent mechanism (a) is **not** available here. A class `c_q + log(varpi_x/conj varpi_x)`
       has the parity of `(ell_x/q)`, not the forced parity of `(N_x/q)`, and the two agree for all `x` only
       at a sparse set of `q in Q_1`.
     * Level data at `q = p_j <= X` (same-level pairs only) can be chosen injective.

     This would give P3-char caps with `p_j = M^(2 + o(1))` and a constant independent of `A`, valid for
     `M >= M_0(A)`. That suffices to block (T') for P3-char. **Neither mechanism is written out; both are to
     be checked.**
   * *Why P3-all\* is different.*
     * Its sparse monomials of norm `k` have budgets `<= (k h_0/2) log M`, independent of `A`, and there
       are `e^(Theta(k log M))` of them.
     * Per-modulus avoidance handles them only at `q >~ (2r)^k`.
     * For independent labels the union bound fails for every `k` from about `2A/h_0` up to about
       `(2A/h_0) M^(2A/h_0)`. This range **grows** with `A`: it is the D1 phenomenon (`sharp.md`
       Prop. 6.5) in regime (1.4).
     * The single columns have Legendre patterns of conductor `p_j << X`, so Thm A excludes small-height
       coherence outright.
   * *Reading (heuristic).* If R2 is confirmed, R1 is only a quantifier repair for rows 2-4. Then
     P3-all\* (with P3-G inside it) stays the frontier, as round 7 said. What separates it from P3-char in
     regime (1.4) is the sparse monomials together with the parity law at small conductor.

4. **README Proposition H (elementary, unrefereed): bounded norms see only bounded moduli.** Fix genuine
   primes `p_1, ..., p_r` and `K >= 1`, and put `Q_K := 8K(2r+1)^K + 1`. Build a label system as follows.
   * At `q <= Q_K` (`q` not dividing `N`), take the genuine labels `log(pi_j/conj pi_j)`.
   * At each `q > Q_K`, choose `ell_1(q), ..., ell_r(q)` in turn, each in its forced parity class. Avoid
     every value for which some `v` with `v_j != 0`, `supp v ⊂ {1, ..., j}` and `||v||_1 <= K` collides
     at level 1 with some unit.

   Then the system is parity-respecting. No monomial of norm `<= K` collides at any level at any `q > Q_K`.
   So every such monomial has per-unit demand at most its genuine demand, and the system meets the per-unit
   height condition, and `beta`-goodness with the genuine `beta`, for all `||v||_1 <= K`, at every `X`.
   * *Proof.* There are at most `(2r+1)^K` such `v`. For each `v` and each of the 4 units, the equation
     `v_j ell_j = c mod m_q` has at most `gcd(v_j, m_q) <= K` solutions. So at most `4K(2r+1)^K < m_q/2`
     values are forbidden, while the parity class has `m_q/2` elements. QED
   * At every greedy step at least half of the parity class is admissible. So above `Q_K` the labels are
     free up to that constraint; for example they can be random. No coherence is forced there.
   * *Exact check* (`README_checks/propH.py`, output `propH_out.txt`).
     * The greedy systems at the first three primes above `Q_K` were verified exhaustively (every
       `||v||_1 <= K`, every unit) for `(r, K) = (4, 2), (4, 3), (6, 2)`.
     * *Evidence:* greedy already succeeds at every prime above 41, 103 and 47 respectively. So the true
       threshold is far below `Q_K`, about the number of monomials, but it is still `M^(K + o(1))` when
       `r ~ M`.

   *Consequences in regime (1.4).*
   * **Restated L4 questions.** The referee left two questions open in regime (1.4): whether L4 fails up to
     norm 2, and which norms an inverse theorem must use. Here are the answers.
     * For every `K`, norm-`<= K` goodness does not force coherence at moduli above `Q_K = M^(K + o(1))`
       (`r ~ M`). Every such counterexample is genuine below `Q_K`.
     * An inverse theorem that certifies coherence up to `Y` must use norms `K` with `Q_K >= Y`. To reach
       the range of Thm A, i.e. `Y >= C_L (8 p_j Norm zeta)^L`, it must use
       `K >~ L (h_0 + h) log M/log(2r+1)`, which is about `L (h_0 + h)` for `r ~ M`.
     * This replaces the withdrawn "norm `>= 3` and deep 2-adic moduli".
   * **Why D fails and H does not.** D's collision moduli are bounded only by the 16th-power certificate
     `Im(Z_v^16)`, of size about `|Z_v|^16`. When `p_j <= M^(h_0)` this exceeds the budget `e^(w(v)/2)`,
     which is the referee's regime mismatch. Prop. H is genuine wherever collisions are forced and
     collision-free above `Q_K`.
   * **Model remark.** For labels independent across moduli, a single monomial of budget `B` meets it with
     probability `<~ B/log X` (Mertens: one collision at any `q > e^B` is too many). So random labels fail
     already at norm 1 in regime (1.4). Prop. H shows that this is a defect of randomness, not a
     constraint.

**Status by class and regime** (`A`-uniform means: for every `A`, with `p_j <= M^(O(1))`).

| class | fixed `A` (`p_j > X` allowed) | `A`-uniform, regime (1.4) (what growth needs) |
|---|---|---|
| `A_0` (row 1) | capped (skew double core, `8 + o(1)`) | capped (no `X`) |
| `A_X` (rows 2-3) | capped, constant `Theta(A)` | not proved; cap expected by sketch R2 (a) or (b) |
| P3-char (rows 4, 4a) | capped, constant `b max(A,2)` | not proved; cap expected by sketch R2 (b), since its non-pair characters have budgets `Theta(M log M)`; the small-coherent route (a) is blocked by parity |
| P3-all\* (rows 5, 5a) | open, bracket above | **open**; small-height coherent designs are genuine-type (Thm A); bounded norms force nothing above `M^(K+o(1))` (Prop. H); the inverse step at the scale of Thm A is open |
| P3-G (row 6) | open | open; it is the coherent alternative of L4 in regime (1.4) |

**Bottom line for P3-all\*.** It remains the strongest open local class. Round 9 locates the problem in
regime (1.4), where:
* the certification side is settled unconditionally (L6 done);
* the coherent alternative is genuine-type, i.e. route 4 and L7;
* the inverse step L4 is open at the only scale that can feed Thm A: norms `K >~ L(h_0 + h)` and moduli
  up to `Q_K` (Prop. H).

Remark R1 is a quantifier correction to the barrier map. Sketch R2 says that it is probably harmless for
rows 2-4. Until R2 is written out and refereed, the caps of rows 2-4 should be quoted as "capped at each
fixed `A`".

## 6. Next most promising step

**Recommended: the inverse step L4 in regime (1.4), at the scale fixed by Prop. H.**

*Statement to decide.*

> **L4(K).** Let `p_j <= M^(h_0)` and `r <= rho M`. Suppose a parity-respecting label system is `beta`-good
> for every `v` with `||v||_1 <= K` at the moduli `q <= Q_K` (`Q_K = M^(K + o(1))`). Is it then
> **system-level** coherent with Gaussian integers of norm `<= M^h` on all but a small part of the
> `Q_7` log-mass below `Q_K`? System-level means one unit `u_q` per modulus for all columns, plus the F1
> operations.

Take `K` with `Q_K >= C_L (8 M^(h_0 + h))^L`, i.e. `K` about `L(h_0 + h)`. Then L4(K) together with Thm A gives
genuine height (`p_j Norm zeta_j` a square or twice a square) on all but a few columns.

*Why this step.*
* After round 9 it is the only open step of route 3 that is purely combinatorial. In regime (1.4):
  * L6 is proved unconditionally for exact, Hamming and few-piece coherence (§1 items 2-4);
  * the coherent alternative is genuine-type, i.e. P3-G and L7;
  * L5 is localised to `q < 2cM + 1`.
* Prop. H says exactly where L4 must act: norms up to `K` and moduli up to `M^(K + o(1))`. Bounded norms
  cannot see further, so any inverse theorem lives at that threshold. The question there is a
  load-balancing problem.
  * Pigeonhole forces about `#{||v||_1 <= K}/q` collisions at each `q <= Q_K`.
  * The budgets allow each monomial only about `h_0/2` collisions near the top of the range.
  * Genuine labels meet this through coherence: `v` collides at `q` iff `q | D_s(A_v)`. Independent labels
    fail the union bound (§5 item 4, model remark).
* Either answer is decisive for route 3.
  * *A proof* is the first inverse theorem. With Thm A it reduces route 3 on near-uniform-margin members
    of P3-all\* to L7 (Poisson tails of genuine-type data), plus L5 on twin profiles.
  * *A counterexample family for all `K`* is a non-coherent "multiplicative code with prescribed Legendre
    parities" at the threshold scale. It refutes route 3 as formulated. If it can be completed to the slab
    conditions, it caps P3-all\* in regime (1.4), i.e. C1 in the regime growth needs.

*First finite test (exact; evidence).* Take `r = 4-8` split primes `<= 200`, `K = 2, 3`, and moduli up to
`Q_K` (`Q_2 = 1297` and `Q_3 = 17497` at `r = 4`).
* Search for parity-respecting systems that meet every norm-`<= K` per-unit budget but are system-level far
  from every coherent system on the top half of the log-mass. Use greedy or ILP load-balancing of the
  forced collisions across moduli.
* Compare their maximal normalised demand with that of the genuine labels on the same primes.

This replaces the referee's regime-mismatched test of Construction D. The referee's own regime-(1.4)
instances (`referee_L4_checks/`, primes `<= 400`, `X <= 5*10^4`) give the harness.

**Also to do (cheap, required for the record).**
1. Write out and referee README remark R1, sketch R2 and Prop. H (§5).
   * R1 corrects the quantifier reading of the barrier map: the residue caps are proved only at fixed `A`.
   * R2 should restore rows 2-4 uniformly in `A`. Until then they are "capped at each fixed `A`".
   * Prop. H replaces the withdrawn norm statements of `L4.md` §6.2.
2. **L5 as a finite problem.** `L5_med(c) => L5` holds unconditionally, by restriction (Lemma 2.5). The
   converse holds given Hall. So an obstruction at the moduli `q < 2cM + 1` proves L5, and with it a restricted
   growth theorem: no near-extremal twin-profile actual clusters with `O(log M)` margins for large `M`. The first data would be exact
   feasibility tests (ILP/SAT) of (F), (T), (P) on `Q_med` for `M <= 16`.

**Lower priority, or not recommended.**
* Route 1 as stated *is* the growth statement for actual clusters (Cor. 3.4). Further work on L1 alone
  (bisector-type bounds towards (R1-arch)) cannot give growth without an arithmetic half.
* The spectral form of L1 is vacuous (Thm A of `p3lower.md`).
* L7/`H` is the common analytic bottleneck of every route: L5_med, (R1-arith) and (I3). No new handle on it
  appeared this round.
* C2 (skew double core with P3-all data) is still undone. By sketch R2(b), P3-char data alone probably
  leave `D_max = O(log M)` there, so C2 is informative only with the full P3-all monomial conditions.

## 7. Files

All paths are relative to this folder (`round9/`). Earlier rounds are under `docs/claude_review/` in the
repository.

| file | content | status |
|---|---|---|
| `L4.md`, `L4_notes.md`, `L4_checks/` (`check_D.py`, `check_AC.py`, `tq.py`) | L4 write-up: Theorems A-C, Prop. D, (I1\*) | refereed: theorems survive; regime claims withdrawn |
| `referee_L4.md`, `referee_L4_checks/` (`c1`-`c7`, `consD.py`, `glib.py`, `notes.md`) | report; regime-(1.4) failures of D, random baseline 0.54-0.55, no-norm-bound counterexample | report |
| `L5.md`, `L5_notes.md`, `L5_checks/` (`check_identities.py`, `check_volume.py`, `check_nontwin.py`) | L5 write-up: identities, localisation, volume inequality, non-twin profiles | refereed: survives with six major corrections |
| `referee_L5.md`, `referee_L5_checks/` (`ref_hall.py`, `ref_identities.py`, `ref_nontwin.py`, ...) | report; Hall evidence, Prop. 5.3 counterexamples, MILP section norms, unweighted-wording examples | report |
| `L1L2.md`, `L1L2_notes.md`, `L1L2_checks/` (`scan.py`, `stats.py`, `known.py`, `margins_detail.py`) | L1/L2 write-up: Thm 1.1, Lemma 1.3, Thms 2.1, 3.2, 3.3, Cor. 3.4, data | refereed: mathematics survives; headline and two data readings withdrawn |
| `referee_L1L2.md`, `referee_L1L2_checks/` (`scan_ref*.py`, `records*.py`, `golden7.py`, `neg_mubar.py`, ...) | report; independent scan (22,907 clusters), records including the 10-point record | report |
| `README.md` (this file), `README_notes.md` | synthesis; R1, R2, Prop. H (unrefereed) | coordinator |
| `README_checks/propH.py`, `propH_out.txt` | exact check of Prop. H | exact (small cases) / evidence (threshold) |

*Citations to check* (from `L4.md`, flagged by both author and referee):
* Linnik's theorem with an explicit exponent `L` (Xylouris: `L = 5`);
* Heath-Brown, "A mean value estimate for real character sums", Acta Arith. 72 (1995), for the quadratic
  large sieve.

They affect only how large `A` must be.

*Caveats.*
* §§4-6 contain coordinator remarks (R1, R2, Prop. H) that no referee has seen. R1 and Prop. H are
  elementary. R2 is a sketch.
* Before any of them is quoted as proved, R1 and Prop. H should be refereed, and R2 written out.

---------------------------------------------------------------------------------------------

**Repository note.** Three pickled data files over 1 MB are not committed
(`L1L2_checks/s9.pkl`, `L1L2_checks/s9_stats.pkl`, `referee_L1L2_checks/ref_1e9.pkl`). The scripts in
the same directories regenerate them.
