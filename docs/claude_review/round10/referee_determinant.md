# Referee report on round10/determinant.md

Referee: adversarial review of `round10/determinant.md` (2026-10-13). I re-derived every proof and tested the
claims with my own exact computations. Scripts and outputs are in `round10/referee_determinant_checks/`.

## 0. Verdict (short)

**The result survives as stated (partial: no growth improvement).** I found no error that affects a theorem.

* Theorem 2.1, Theorem 2.3, Corollary 2.4, Theorem 4.1, Proposition 4.2, Corollary 4.3, Lemmas 3.1-3.2 and the
  conditional Theorem 5.7 are correct.
* The exact thresholds (`n = 2137, 2679, 4278`) reproduce with exact big-integer binomials.
* The Liouville/forced-divisibility chain behind Theorem 4.1 passes an end-to-end exact test on 240 kernel
  vectors over actual lattice clusters.

The corrections are of interpretation and wording.

1. **(Major, interpretive.) The fake separation is "almost every", generic, and not by itself progress.**
   * Corollary 4.3 is true, but only because a nonzero trigonometric polynomial has a null zero set. The trivial
     property "the points are Gaussian integers" separates in exactly the same way.
   * "*Every* fake violates (4.2)" is not proved. It is a special case of the unproved Hypothesis ZE.
   * The separating input is integrality plus forced divisibility of multi-term forms, already listed as outside
     the barrier in `outside.md` §6 item 1.
   * The genuinely new content is that this input is LP-supercritical in the growth regime (Theorem 2.1).
2. **(Minor.) Theorem 4.1 may be vacuous.** It needs actual clusters of at least `2137` points on arcs of length
   `C sqrt R`, and none is known. The write-up acknowledges that it gives no bound.
3. **(Minor.) ZE evidence and status.**
   * At the author's evidence sizes, `V_H = V(A,m)` exactly (`ref_VH_small.py`), so its caveat is moot there.
   * In the regime, `V_H` may be strictly smaller and `Z_P` larger.
   * I add a rigorous positive case (Proposition R1): ZE holds when the slabs do not cut. This shows the whole
     difficulty lies in the slab cuts.
4. **(Minor.)** Wording and constant slips:
   * Theorem 4.1's proof applies `|c|_1 <= UH` to every integer vector in the span;
   * Lemma 3.2 cites the `(pi/2)Z` gap where the `(pi/4)Z` gap is needed;
   * Prop 5.1's "65% coincide" must say "with the base point `n`" (demonstrated numerically);
   * Prop 5.2's application silently takes `kappa = 0` (justified here);
   * Prop 5.3(iii) has the wrong constant;
   * Prop 5.6's "at best critical" is heuristic, and its trailing remark is garbled.

The corrected statement is in §7.

## 1. What the write-up claims

* (1) **Theorem 2.1.** For every profile with `r` columns, the slab-box support gives
  `ord/sum omega_p h(beta_p) >= 0.12 sqrt((n-1)/ln 4r)`. Theorem 2.3 bounds the Hadamard profile by `sqrt n`.
  Corollary 2.4 says LP-critical profiles need `r >= e^{0.0144(n-1)}/4`.
* (2) **Theorem 4.1 (unconditional rigidity).** Every actual cluster with `n - 1 >= 300 ln 4r` and
  `W >= 16 ln n + 8 log^+ C + 14` satisfies `sum_{k in A} c_k e^{i<k,delta>} = 0` for all `c in V_H`, where
  `dim V_H >= |A|/2`.
* (3) **Corollary 4.3.** (4.2) is false for almost every round-5 fake and skew-double-core fake.
* (4) **Conditional.** Lemma 3.1 (rank `>= n/4 - 2` on squarefree norms) and Lemma 3.2 (a nondegenerate full-rank
  subcluster) together with Hypothesis ZE give `n = O(log W)` on squarefree norms (Theorem 5.7).
* §5 shows that the natural nonvanishing certificates fail or are capped.

The claimed fake-violated input is the rigidity relation (4.2), equivalently the vanishing of every
forced-divisibility auxiliary form in `V_H`.

## 2. Theorem 2.1 (LP supercriticality), Theorem 2.3, Corollary 2.4

**Verdict: correct.** I re-derived every step.

* (a) Hoeffding. Each summand `k_x(beta_{p,x} - beta_{p,n})` has range length `2T|beta_{p,x}-beta_{p,n}| <= 4T`, so
  `P(|X_p| > w) <= 2 exp(-2w^2/(16(n-1)T^2)) = 2 exp(-ln 4r) = 1/(2r)`. The union bound over `r` columns gives
  `|A| >= |B_T|/2`. Correct.
* (b) The closed form `m_0` follows from `binom(a+b,b) <= (e(a+b)/b)^b`; `e 8^{1/63} = 2.8095 < 2.8097`;
  `0.7118n - 0.644 >= 0.7n` needs `n >= 54.6`. Correct.
* (c) Counting Siegel. The pigeonhole count, the bound "a `d`-dimensional subspace meets `[-H,H]^U` in at most
  `(2H+1)^d` integer points", and the final inequality `d >= U(1.25L - 0.275)/(2L + 1.1) >= U/2` (needs only
  `L >= 3.3`; here `L > 300`) are all correct.
* (d) `4 sqrt(8/300) = 0.653 < 0.7`. The general ratio is `0.7/(2 sqrt 8) = 0.1237`. Correct.

**Exact check** (`referee_determinant_checks/ref_params.py`, output `ref_params_out.txt`). I computed `m` with exact
big-integer binomials (`8 binom(m+n-2,n-1) <= (2n+1)^{n-1}`), independently of the author's `lgamma` code.
* `m >= 0.7 n(n-1)` holds at every tested `n` from 55 to 4278.
* The thresholds for `m >= 4w` at `T = n` are exactly the author's: `n = 2137` for `r = n`, `2679` for
  `r = floor(n ln n)+1`, and `4278` for `r = n^2`. The failure rows are `n = 2136, 2678, 4277`, with exact `m` and
  `4w` printed.
* At the boundary `r = floor(e^{(n-1)/300}/4)`, `m >= 4w` holds for `n = 2100, ..., 6000`, as Theorem 2.1(d)
  predicts.
* The ratio `m/(2w)` is at least `0.12 sqrt((n-1)/ln 4r)` in every row.

**Theorem 2.3** (Hadamard profile, `ord <= sqrt n Lambda_H`) is correct. It uses `H H^T = nI`, that `h_A` is a
seminorm with `h_A(1) = 0` on a slice, `E|<eps,h_p>| <= sqrt n`, and (I2). **Corollary 2.4** is correct as stated.

**Relation to the growth table.** The correction to `round8/README.md` §2.4 is fair. The `threshold/`
criticality concerns a form in `M` variables that is fixed and is evaluated at an *uncontrolled* tuple through
the grid lemma. In this write-up the form lives in all `n` cluster variables and is chosen after the profile.
The price is that nonvanishing is lost. `round8/README.md` §2.5 item 5 already records that "no proved cap
covers" single-place auxiliary polynomials in the growth regime, so Theorem 2.1 contradicts no refereed
statement. The table entry should read "exactly critical for a form fixed before the profile (or used through the
grid lemma); supercritical at the LP level when the form is chosen after the profile, with nonvanishing open".

**Size of the regime.** The theorem needs `n >= 2137` cluster points even at `r = n`, and `two_thirds` then forces
`W = log N >= 3 n log n - (14 + 4 log^+ C) n`, about `1.9 * 10^4` at `n = 2137`. No actual cluster of that size is
known, so the regime is purely asymptotic (see §3).

## 3. Theorem 4.1 (rigidity): re-derivation

**Verdict: correct, with one wording slip.** This is a Thue-Siegel dichotomy, and each link checks.

1. *Monomials.* `D = 2(n-1)n` is even and `|k|_1` is even on `Z^n_0`, so `G_k(z) in Z[i]`. On the circle,
   `G_k(z) = N^{(D-|k|_1)/2} R^{|k|_1} e^{i<k,theta>} = N^{D/2} e^{i<k,delta>}`, because `sum k = 0`.
2. *Taylor.* `<k,delta> = sum_{x<n} k_x(delta_x - delta_n)` is linear in `k'`, so `<k,delta>^t` has degree `t`
   and `c in V(A,m)` annihilates it for `t < m`. Then
   `|e^{ix} - sum_{t<m}(ix)^t/t!| <= |x|^m/m!`, and `|<k,delta>| <= (|k|_1/2) Theta <= D Theta/2`.
3. *Valuations.* `v_pi(G_k) = e_p(D-|k|_1)/2 + sum k^+ a + sum k^-(e_p - a) = e_p D/2 + <k,a_p>`. The inert part
   is `q^{fD}`, and `(1+i)^{sD}` likewise. So `Norm F_c >= N^D prod p^{-h_{supp c}(a_p)}` whenever `F_c != 0`.
   Also `<k,a_p> = (e_p/2)<k,beta_p>` on the slice.
4. *Comparison.* `(m/2 - w) W <= 2 ln|c|_1 + 2m ln(DC/2) - 2 ln m!`. With `m >= 4w`, `ln|c|_1 <= 3 ln U + 2(m-1) ln n`
   and `ln m! >= m ln(m/e)`, this gives `W <= 16 ln n + 8 ln C + 13.46`, contradicting (4.1).
   * Exact check (`ref_params.py`): the right side minus `16 ln n` is at most `12.79` (at `n = 64`) and decreases
     to `10.5` as `n` grows.
   * `eD/(2m) <= 3.69` at the exact `m`, against the `3.884` used.

**Slip (minor).** The proof says "Let `c in V_H` be an integer vector ... `|c|_1 <= UH`". An integer vector in the
*span* `V_H` need not have `|c|_inf <= H`. The argument should be run for the generators, the integer
`c in V(A,m)` with `|c|_inf <= H`, and then extended to `V_H` by linearity. The proof's last sentence ("By
linearity the same holds on the span") does exactly that, so the conclusion is unaffected. Only the first
sentence of the proof needs rewording.

**End-to-end exact test on actual lattice clusters** (`ref_chain.py`, output `ref_chain_out.txt`). This is
independent of the author's `exact_checks.py`, which tested only the valuation formula and the face congruence.
* Data: 5 norms (two with `e_p = 2`), 4-5 consecutive lattice points, random supports `A` in the `T`-box, and
  240 integer kernel vectors of the order-`m` moment matrix (`m = 3, 4`).
* Results, all exact or at 60 digits:
  * (P) the phase identity `F_c(z) = N^{D/2} sum_k c_k e^{i<k,delta>}` holds, with `delta` measured from the arc
    midpoint;
  * (U) `|F_c(z)| <= N^{D/2} sum|c_k||<k,delta>|^m/m! <= N^{D/2}|c|_1(D Theta/2)^m/m!` holds;
  * (L) `Norm F_c(z) prod_p p^{h_{supp c}(a_p)} >= N^D` holds, as an exact integer inequality.
* All 240 pass. No `F_c` vanished: at `n <= 5` the LP ratio is below 1 (§2), so no vanishing is forced, as
  expected.

**Is the theorem vacuous?** It is a correct implication, but its hypothesis needs an actual cluster with
`n >= 2137` points (exact threshold at `r = n`; the proof's sufficient condition `n - 1 >= 300 ln 4r` needs
`n >~ 2700`) on an arc of length `C sqrt R`, hence `log N >= 3n log n - (14 + 4 log^+ C)n`, about `1.9 * 10^4`. No such cluster is known.
Either no actual cluster lies in the regime, which is the growth statement itself, or such clusters satisfy
(4.2). The theorem cannot tell these apart. The write-up says this ("It does not produce a contradiction"),
but the §0 phrase "a rigidity theorem for actual clusters" should carry the qualifier "if any exist in the
regime".

**What the theorem really uses.** Forced divisibility is essential. With integrality alone (`|F_c| >= 1`), the
same comparison needs `m > 2D`, which `threshold/upper.md` (`tau* < 2`) forbids. The slab-box lowers the
requirement from `2D` to `4w = o(m)`. Neither residues, prime sizes nor the pair inequality are used.

## 4. Fake test: does (4.2) fail for the fakes?

**Proposition 4.2 and Corollary 4.3: correct as stated ("almost every fake").** I checked them against the source
constructions.

* *Round-5 fakes* (`round5/construct.md` Thm 2, steps 3-5). `delta` is drawn uniformly from `[0,Delta]^M`, the good
  set has measure `> 1/2`, and the profile has rank `M` (twins), so every `delta` is realisable. `phi` is chosen
  afterwards, on the realisation coset (Theorem A.1(ii)).
* *Skew double core* (`round5/exact.md` Thm 5.1, steps 3-5). The same mechanism applies: Lemma 5.0, a union bound
  `<= 0.414`, and full rank, because `1_{S'_a} - 1_{S_a} = e_a - e_0`.
* `Z_P` is invariant under `delta -> delta + t 1`. "Positive Lebesgue measure" should therefore be read in
  `[0,Delta]^M` (the fakes' parametrisation) or in `R^n/R1`. With the write-up's convention of centring at the
  arc midpoint, `delta` lies in the null set `{max + min = 0}`. This is cosmetic: Fubini along `1` fixes it.
* *Regime.* (4.1) needs `M - 1 >= 300 ln(4r)`. That is `M >= 3527` for the round-5 profile (`r = 9M - 8`; least
  Paley order `3528`) and `M >= 3030` for the skew double core (`r = 2(M-1)`; least Paley order `3068`). Both
  families exist there.
* *Fakes are nondegenerate.* For a full-rank profile and a nonzero zero-sum `c`, `v(c) != 0`. If `c.delta = 0`,
  then `<v(c),phi> = (pi/4) c.k` lies in `(pi/4)Z`, which the fakes' Lemma 0.1 gaps forbid. So the degenerate
  strata of Propositions 5.1-5.2, which do lie in `Z_P`, contain no fake.

**Small-scale illustration** (`ref_fakes.py`, output `ref_fakes_out.txt`; evidence, far below the regime).
* Profiles: the skew double core at `q = 7` (`M = 8`, `r = 14`) and the round-5 Paley-I profile with `b = 8` at
  `q = 7` (`r = 64`).
* Rank: both have exponent rank `7`, which is full.
* Support: a slab-cut box `T = 1`, `m = 5`, `|A| = 1765` and `1781`, `dim V(A,m) >= |A|/2`, exact integer basis.
* Positions: 200 draws `delta ~ U[0, 10^-6]^8`.
* Result: every one of the 400 positions is *certified* off `Z_P`. Some `c` has
  `|sum_k c_k <k,delta>^m/m!| > sum_k |c_k||<k,delta>|^{m+1}/(m+1)!` (the leading Taylor term beats the
  remainder), so `sum_k c_k e^{i<k,delta>} != 0`.

**But the separation is generic, and the write-up should say how weak that is.**

1. *"Almost every", not "every".* Corollary 4.3 does not exclude a fake whose `delta` lies on `Z_P`. Fakes are
   nondegenerate and have `max|delta_x| <= C e^{-W/4} << 1/(4D)`. So **"no fake lies on `Z_P`" is a special case of
   Hypothesis ZE**, and a fake on `Z_P` would refute ZE. This matters for the barrier logic.
   * The missing nonvanishing step could be proved from angle-model inputs alone (the classes (I)-(IV) of
     `outside.md` Cor A.4) only if no fake lies on `Z_P`.
   * If a barrier construction could be steered onto a nondegenerate point of `Z_P`, the barrier would transfer
     to the nonvanishing step.
   * Nothing is known either way. The good sets are closed sets of positive measure, and the nondegenerate part of
     `Z_P` near `0` is unknown.
2. *The violation has no quantitative content.* Proposition 4.2 says only that a nonzero trigonometric polynomial
   has a null zero set. Every exact analytic identity forced on actual clusters is violated by almost every fake in
   the same way. An example is "`R e^{i theta_x} in Z[i]`", which almost every fake also violates. So Corollary 4.3
   certifies that the mechanism is *not blocked* by the barrier. It does not certify progress.
   * The input that does the separating is the one `outside.md` §6 item 1 already lists as outside the barrier:
     integrality of multi-term forms `F_c(z)`, together with their forced divisibility.
   * For a fake, `F_c := N^{D/2} sum c_k e^{i<k,delta>}` is defined, but the dichotomy "`F_c = 0` or
     `Norm F_c >= N^D prod p^{-h}`" fails whenever `g_c(delta) != 0`. In the regime, the computation of Theorem 4.1
     shows it fails for almost every fake.
   * What is new is Theorem 2.1. In the growth regime this input is *supercritical* at the LP level, so it yields
     an exact relation rather than nothing. That is a genuine and correct observation.
3. *The escape from `outside.md` Cor A.2 is correctly explained.* Class (IV) consists of properties holding off a
   Haar-null subset of each realisation torus. (4.2) depends only on `(P, delta)`, which is constant on each
   realisation torus, and it holds only on a null set of `delta`. So it is not in class (IV).

**Answer to the mandatory question.** The input is the polynomial forced-divisibility dichotomy for `F_c`,
`c in V_H`. Equivalently, it is its consequence (4.2).
* It is true for every actual cluster in the regime (Theorem 4.1, which I verified).
* It is provably false for almost every round-5 fake and almost every skew-double-core fake at the same
  `(n, r, W)`.
* It is **not** proved false for *every* fake. That statement is part of ZE.

## 5. Lemmas 3.1, 3.2, Theorem 5.7 and Hypothesis ZE

**Lemma 3.1 (rank `>= n/4 - 2`): correct.**
* `G = S diag(w) S^T` has `G_xx = 1` and `G_xy = 1 - 2d_xy/W'`.
* (I4) is the gcd argument `Norm gcd(z_x,z_y) = N' e^{-d_xy} <= |z_x - z_y|^2 <= C'^2 N'^{1/2}`, which I re-derived.
* `P <= n(n-1)eta <= n`, `Q <= n + P`, `||G||_F^2 <= 4n`, and `rank >= (tr G)^2/||G||_F^2`.
* `rank S <= rank(a) + 1` and `k >= rank(a) - 1`.

**Lemma 3.2: correct, with a citation fix.** The proof needs `<2u,phi> notin (pi/2)Z`, that is
`<u,phi> notin (pi/4)Z`. That is the *second* inequality of `outside.md` Lemma 0.1 (the distance to
`pi/4 + (pi/2)Z`), not the `(pi/2)Z` gap. Both inequalities are in Lemma 0.1, so nothing breaks. The same remark
gives multiplicative independence of `w_x = z_x conj(z_n)/N`, which §7 item 1 uses.

**Theorem 5.7: a correct conditional.**
* Lemma 3.1 gives rank `k >= n/4 - 2`, Lemma 3.2 a nondegenerate subcluster `J` with `n' = k+1` points, and
  Theorem 4.1 applies to `J` (same `N`, same `r`, `n' <= n`).
* For ZE one needs `max|delta_x| <= C' e^{-W'/4}/2 <= 1/(8n'(n'-1))`. This is `W' >= 4 ln(4C'n'(n'-1))`, which is
  implied.
* The conclusion `n < 4 max(n_1, 300 ln 4r + 1) + 4` follows, with `r <= W'/log 5`.
* Nothing is hidden. The whole weight is on ZE.

**Hypothesis ZE: status and evidence.**
* *The author's evidence is about the right object at its own sizes.* The write-up warns that its scripts use the
  full `V(A,m)`, not `V_H`. I computed an exact integer basis of `V(A,m)` for the three evidence instances
  (`ref_VH_small.py`). Its sup-norms are `2835`, `87584` and `1536`, against `H = 2.7 * 10^8`, `1.4 * 10^11` and
  `3.4 * 10^7`. So `V_H = V(A,m)` there, and the caveat is moot at those sizes.
* *The regime differs.* There, `V_H` can be a proper subspace: Siegel only guarantees half the dimension, and the
  slab cuts rule out the short "product finite-difference" vectors (see below). That makes `Z_P` larger. The
  evidence does not reach this situation.
* *A rigorous positive case (new; Proposition R1).* Suppose the slabs do not cut (`A = B_T`; for example, every
  column is constant on the slice) and `2^{2T(n-1)} <= H` (true for `T = n >= 64`). Then
  `Z_P = {delta : 2T #{x < n : delta_x != delta_n mod 2pi} < m}`. Every point of this set with `n >= 3` is
  degenerate, so ZE holds for such profiles.
  * *Proof.* In the Newton basis, `f_delta(k) = prod_x zeta_x^{-T} sum_{alpha <= 2T} binom(k+T, alpha) prod_x (zeta_x-1)^{alpha_x}`,
    with `zeta_x = e^{i(delta_x - delta_n)}`.
  * `V(B_T,m)` is spanned by the dual functionals of `binom(k+T,alpha)` with `|alpha| >= m`. These are the
    Newton forward differences `Delta^alpha`, with integer entries `<= 2^{|alpha|} <= 2^{2T(n-1)} <= H`.
  * Hence `V_H = V(B_T,m)`. Then `f_delta in V_H^perp = Poly_{<m}|_B` iff every `prod_x(zeta_x-1)^{alpha_x}` with
    `|alpha| >= m`, `alpha <= 2T` vanishes. QED
  * Numerical consistency: `ref_fullbox.py`.
* *Where the difficulty sits.* `V(A,m)` is contained in `V(B_T,m)` (extend by zero), so cutting the box only
  *enlarges* the base locus. Any sub-box of `A` with `|alpha| >= m` would have width about `2|alpha| > 2w` in a
  generic `beta_p` direction, so the short `Delta^alpha` vectors do not fit in `A`. ZE is exactly the claim that
  the cut-out half-box, at low height, still has no non-subtorus components near `1`. I could neither prove nor
  refute this. The write-up's description of it as an unlikely-intersection statement is fair. Note also that ZE
  must exclude translates `gH` (with `g` near `1`, `g notin H`) of proper subtori, not only subtori through `1`.
  The anomalous-intersection structure theorems (Bombieri-Masser-Zannier) do not exclude such translates
  effectively.
* *ZE contains the fake question.* As shown in §4, ZE implies that no fake lies on `Z_P`.

## 6. Section 5 obstructions (negative results)

These are statements that nothing closes the gap. I checked each.

* **Prop 5.1: correct, but the "65%" gloss is wrong as stated.** The box `B_T` bounds `|k_x| <= T` only for
  `x < n`, so the base point `n` is special.
  * Moving point `n` alone makes all `zeta_x` equal to some `zeta != 1`. Then `f_delta(k) = zeta^{-k_n}` with
    `|k_n| <= (n-1)T`. This is not in `Poly_{<m}` unless `2T(n-1) < m`, which fails.
  * The correct gloss is "every `delta` in which more than 65% of the points coincide *with point `n`*".
  * Exact demonstration (`ref_basepoint.py`): at `(n,T,m) = (4,1,3), (5,1,3), (5,2,5), (6,1,4)`, moving one point
    `x < n` gives residual `<= 2 * 10^-15`, while moving point `n` gives residual `0.36-0.74`. The author's own
    `n = 5` data agree: 85 zeros, all on coordinate axes, none on the diagonal.
  * This changes no theorem.
* **Prop 5.2: correct, with one unstated step in its application.** For an actual cluster of rank `k`, the
  write-up writes `delta = sum_{p in S} t_p a_p + (pi/2) kappa + t_0 1` and then drops `h_A(kappa)` when it says
  "with bounded exponents this means `k <~ 0.25 sqrt(n/ln 4r)`". This can be justified.
  * Cramer's rule gives a spanning set of `(span{a_p} + R1)^perp cap Z^n` with `|v|_1 <= n(sqrt(k+1) e_max)^{k+1}`.
  * For `k <~ sqrt n` and `W >= 3n log n - O(n)`, this is `< pi/Theta`.
  * Since `<v,delta> = -(pi/2)<v,kappa>` and `|<v,delta>| < pi/2`, we get `<v,delta> = 0`.
  * So `delta in span{a_p} + R1`, and `kappa` can be taken to be `0`. The write-up should say so.
* **Prop 5.3: (i) and (ii) correct; (iii) has a wrong constant.** By the write-up's own convention,
  `<v,delta> = 2<u,phi> - (pi/2)<v,kappa>`. When `<v,kappa>` is odd, the relevant gap is the distance of `<u,phi>`
  to `(pi/4)Z`.
  * Lemma 0.1 gives `|<v,delta>| >= 2 arcsin(e^{-X_v/2}/sqrt 2)`, not `2 arcsin(e^{-X_v/2})`.
  * Together with `|<v,delta>| <= |v|_1 Theta/2`, the correct conclusion is
    `X_v >= W/2 - 2 log(C|v|_1/(2 sqrt 2))`.
  * The stated `W/2 - log(C|v|_1/2) - 1` is stronger by about `log(C|v|_1)` when `C|v|_1` is large. This does not
    matter for the conclusion ("up to `O(log(C|v|_1))`").
* **Lemma 5.4, Prop 5.5: correct.**
  * I re-derived the face congruence:
    `G_k = pi^{eD/2+<k,a>} conj(pi)^{eD/2-<k,a>} (N/p^e)^{D/2} prod y_x^{k_x}` in the localisation, and
    `conj(pi) == sigma (mod pi)`.
  * I re-derived the "no unique minimiser" argument, which uses the level-sum lemma with `h + 1 <= m` levels.
  * (iii) is a tautology under Theorem 4.1, as the author says.
* **Prop 5.6: the expectation identity is correct.** For `t` uniform in `[n]^M` and a balanced sign column,
  `beta_p o t` is uniform on `{+-1}^M`. The conclusion "a grid-lemma argument is at best critical" is a heuristic
  reading of an averaging statement, not a theorem, and should be labelled as such.
  * The trailing remark on "`M = 1` ... `h(eps) = 2(T_s(v) - B_s(v)) >= ...`" uses undefined notation (`T_s`,
    `B_s`), and "`M = 1`" clashes with `M` variables. I could not verify it as written.
* **§5.6 table: accurate.**

## 7. Corrected statement

**No growth improvement is proved. What is proved is the following.**

1. **Theorem 2.1 (unchanged; verified with exact binomials).**
   * Hypotheses: `n >= 64`, any `r` columns `beta_p in [-1,1]^n`, weights with `sum omega_p <= 1`.
   * Construction: the slab-box `A` (`T = n`, `w = n sqrt(8(n-1) ln 4r)`) and `V_H`.
   * Conclusion: `dim V_H >= |A|/2`, `ord >= m >= 0.7 n(n-1)`, `h_A(beta_p) <= 2w`, and
     `m/(2w) >= 0.12 sqrt((n-1)/ln 4r)`.
   * `m >= 4w` holds iff `n >= 2137` at `r = n`, `2679` at `r = n ln n`, and `4278` at `r = n^2`. It always holds
     if `n - 1 >= 300 ln 4r`.
   * Theorem 2.3 (`ord <= sqrt n Lambda_H`) and Corollary 2.4 are unchanged.
2. **Theorem 4.1 (unchanged; proof wording fixed).**
   * Hypotheses: an actual cluster with `n >= 64`, `n - 1 >= 300 ln 4r` and `W >= 16 ln n + 8 log^+ C + 14`.
   * Conclusion: `F_c(z) = 0` for every integer `c in V(A,m)` with `|c|_inf <= H`, hence for all of `V_H`.
   * It needs `n >= 2137` cluster points. No actual cluster in this regime is known, so the theorem may be
     vacuous. It gives no bound.
3. **Corollary 4.3 (sharpened wording).**
   * For a fixed profile `P` satisfying (4.1), (4.2) fails for Lebesgue-almost every `delta in [0,Delta]^M`
     (equivalently, almost every `delta` modulo `R 1`). Hence it fails for almost every round-5 fake
     (`M >= 3527`) and almost every skew-double-core fake (`M >= 3030`), at Paley orders.
   * It is **not** proved to fail for every fake. "No fake lies on `Z_P`" is a special case of Hypothesis ZE.
   * The separation is generic (a nonzero trigonometric polynomial has a null zero set). The input behind it is
     integrality plus forced divisibility of multi-term forms, which `outside.md` §6 item 1 already places
     outside the barrier.
   * What is new and correct: in the growth regime this input is supercritical at the LP level (Theorem 2.1). So it
     yields an exact algebraic constraint on the positions, where before it yielded nothing.
4. **Conditional result (unchanged).**
   * Lemma 3.1: rank `>= n/4 - 2` on squarefree norms.
   * Lemma 3.2, using the `(pi/4)Z` gap of Lemma 0.1.
   * Theorem 5.7: ZE implies `n = O(log W)` on squarefree norms.
   * ZE is unproved and non-standard. It holds when the slabs do not cut (Proposition R1 of this report). The
     small-case evidence is about `V_H = V(A,m)` at `n = 4, 5` and does not reach the regime, where `V_H` may be
     smaller and the base locus larger.
5. **Minor corrections to §5.**
   * Prop 5.1: "65% of the points coincide *with the base point `n`*".
   * Prop 5.2: in the application, `kappa = 0` can be taken by a Cramer-rule argument (§6).
   * Prop 5.3(iii): read `X_v >= W/2 - 2 log(C|v|_1/(2 sqrt 2))`.
   * Prop 5.6: its "at best critical" conclusion is heuristic.

## 8. Scripts (`round10/referee_determinant_checks/`)

| script | what it does | status |
| --- | --- | --- |
| `ref_params.py` | exact big-integer `m`, `0.7n(n-1)`, thresholds for `m >= 4w`, boundary `r`, Theorem 4.1 constant | exact |
| `ref_chain.py` | end-to-end Theorem 4.1 chain (phase identity, Taylor upper bound, forced-divisibility lower bound) on 240 kernel vectors over actual clusters, including `e_p = 2` | exact integer test plus 60-digit archimedean side |
| `ref_VH_small.py` | `V_H = V(A,m)` at the author's evidence sizes (integer basis sup-norm `<< H`) | exact |
| `ref_fakes.py` | full rank of both fake profiles at `q = 7`; 400 fake positions certified off `Z_P` (leading Taylor term > remainder) | rank exact; certification at 50 digits; illustration only (below regime) |
| `ref_fullbox.py` | consistency of Proposition R1 (uncut box: base locus = sparse strata) | floating point, evidence |
| `ref_basepoint.py` | Prop 5.1 base-point asymmetry | floating point, evidence (residuals `1e-15` vs `0.36-0.74`) |

Each script writes its output to the `*_out.txt` file of the same name.
