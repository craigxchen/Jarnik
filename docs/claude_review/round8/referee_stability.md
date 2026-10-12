# Referee report on `round8/stability.md` (stability of near-extremal clusters)

Referee: adversarial, round 8. Scripts: `round8/referee_stability_checks/`. Each runs from that directory with
`python3 <script> [args]`; outputs are stored next to it as `*_out.txt`. The repository, the earlier round
folders and `round8/stability_checks/` were not modified.

## 0. Verdict

**The write-up survives as a restricted theorem plus a barrier. It claims no growth improvement for general
circles, and it proves none.** I re-derived every proof in Sections 1-7 and found no mathematical error. The
exact computations below confirm every checkable statement on actual clusters:
* about 127,000 arbitrary subsets (not only windows) for Prop 2.2 and Lemma 2.1;
* 435,196 windows for Theorem 3.1, exhaustively for `N <= 5*10^6` and `C <= 8`;
* 117,892 further windows, with `C` up to 40 and `M` up to 17, for Theorem 3.3, Props 3.5, 4.1 and 5.1,
  Theorem 5.2, and the algebra of Prop 6.1, Lemma 8.1 and Prop 8.3.

**Checked against the fakes.** The only growth statement is the restricted corollary of Theorem 3.1: if
`log P_+(N)/loglog R -> infinity`, then `M = o(log R/loglog R)`. Its input is a large prime factor of `N`. The
round-5 fakes violate that input: their primes have `log p ~ 3 log M ~ 3 loglog R`. So this corollary is
compatible with them. On smooth norms the transfer (Prop 7.1) is correct, so no conclusion of Sections 2-5
separates actual clusters from the fakes. This is the barrier the write-up states.

What does not survive as stated:

1. **(major, attribution and novelty)** The headline "new" results are, up to the second-order term
   `-3 log(M/2) + K_0`, a two-line consequence of the **Cilleruelo-Cordoba theorem** (Proc. AMS 115 (1992)
   899-905), applied after the same division by a shared Gaussian factor. These results are:
   * the crude gcd lemma;
   * the smoothness theorem `P_+(N) < N^(2/M) C_+^4`;
   * the restricted growth statement;
   * the restricted uniform statement `M < 2/lambda`.

   For odd `s`, Cilleruelo-Cordoba even gives the exponent `1/(s+1)` instead of `1/s` (§2.4). The write-up's
   own proof (through `two_thirds`) is correct and self-contained. Only the "new" labels and the summary's
   framing ("the key new tool") are overstated. Genuinely new: the `F*` refinement of the second-order term,
   the systematic sub-cluster packaging, and Sections 4-5.
2. **(major, Section 8.6 only; no theorem affected)** The pointer "by `round8/weakest.md` (draft), five
   simultaneous forms do [suffice]" rests on a **false** hypothesis. Vertical templates
   `u_l = t/b_l + i` with moment conditions give five independent forms in logarithms of Gaussian-prime
   ratios with `max |Lambda| = (prod H)^(-gamma)`, `gamma -> (2s+1)/(2(5+s))`. That is `1/2` for `s = 4` and
   `13/22` for `s = 6`, and it tends to `1` as `s` grows. So `SLW(5, g, kappa)` is false for every `g < 1`, in
   particular for the `g < 1/4` that `weakest.md` Theorem 2.3 assumes (§6.2; my `r4_slw_templates.py` re-derives
   the construction of `referee_weakest_checks/`). Section 8.6 should not present this as a route to growth.

Minor issues (§7) include:
* the Matveev constant is "quoted from memory" and the standard statement is not given;
* the task summary drops several `C <= 1` and `log^+ C` qualifiers;
* in Section 7, "integrality of the reduced chords is equivalent to the genuineness of the prime angles"
  is only one implication;
* the fake check of Prop 2.2 is profile-level, with `C_B = 1` and not the model's chords;
* "complete content of recursion on sub-clusters" (Remark 2.3(3), §8.5) is an unproved assertion;
* `lw.md` is called "draft" in Remark 6.2(3) and "refereed" in §8.6, and the barrier holds only for the max
  form with `Log i` among the logarithms.

## 1. What was checked, and how

Every check below uses my own code (`rlib.py`), written independently of `stability_checks/stab_lib.py`. It has
its own Gaussian gcd, its own lattice-point generator and its own evaluator of `F*(s)`, a shortest-path DP
cross-checked against brute force. "Exact" means Python integers: every inequality of the form
`gamma + 8F*/s^2 <= ...` is exponentiated, because `8F*(s)` is an integer combination of logarithms of primes.

| item | method | result |
| --- | --- | --- |
| inputs (T1)-(T3) | re-read `two_thirds/proof.md` Thms 5.1, 5.2, 7.5 and §6. (T1) is (3.6)+(2.5) for a normalised configuration, and only chords are used. (T2) is step 2 of Thm 7.5, valid for every `s >= 2` (Lemma 7.2 needs only integers `n >= 2`), with `F* >= D` termwise | correct as quoted |
| `F*(s)` evaluator | `r1_fstar.py`: DP equals brute force for `Psi_p` (`s <= 15`); my `F*` equals `stab_lib.Fstar` exactly for `s <= 80`; `F*(16) = 108.77`, `F*(48) = 1686.65` (two_thirds table) | agrees |
| (T2), (3.0) | `r1_fstar.py`: `8F*(s)/s > 3 s log s - 13.9092 s + 3 + 2 log 2` for `2 <= s <= 250` (min margin `7.87 s`); `8F*(s)/s^2 > 3 log(M/2) - K_0` for `ceil(M/2) <= s <= M <= 250` | holds |
| Lemma 2.1, Prop 2.2 | `r2_gcd_subsets.py`: **127,410 arbitrary subsets**: random subsets of the whole circle (huge `C_B`), random subsets of windows, all subsets of 8-point windows. Ten norms carry content (`2`, `3^2`, `7^2`, `11^2`, repeated split primes). Exact integer forms of (2.1) and (2.2), with levels recomputed on `u = z/G_B` | 0 failures; equality in (2.2) occurs exactly when `s = 2`, `\|beta\|^2 = 2` (4,623 cases) |
| Theorem 3.1 | `r5_smooth_exhaustive.py`: **every** `N <= 5*10^6`, every window with `M >= 4`, `C <= 8`, every `q \| N`. Exact refined form (3.1), and exactly the stronger crude form `log q < 2W/M + 2 log^+ C` (which implies the stated `+ 4 log^+ C`) | 435,196 windows, 0 failures; tightest ratio `log q / max_s b_C(s) = 0.994` (`N = 3301381 = 61*54121`, `M = 4`, `C = 7.98`) |
| Thm 3.3, Props 3.5, 4.1, 5.1, Thm 5.2, algebra of 6.1, 8.1, 8.3 | `r3_window_scan.py` on four ranges (`N <= 3*10^6, C <= 6`; `N <= 10^9, C <= 10`; `N <= 10^11, C <= 40`; `N <= 10^12, C <= 14`), circles with many points and with content. Exact for T3.1, the cut identity, the `Pi` identity of 6.1, the triangle identity 8.1 and the parallel-chord identity 8.3; 40-digit `mpmath` otherwise | 117,892 windows (`M` up to 17), 0 failures |
| Prop 5.1 | the three excesses computed from scratch (`v_q(beta)` by factoring `\|beta\|^2`, collision counts `kappa_(q^a)`, level-wise `E`), identity and `0 <= X <= Xi` | holds in every window |
| `M = 12` fake (Section 7 table) | re-ran `check_fake_structure.py`; re-read its construction against `construct.md` §2 | reproduced; profile-level only (§5.3) |
| Cilleruelo-Cordoba derivation (§2.4) | `r6_cc_gcd.py`: the gcd bounds derived from CC, with the actual arc length, on all subsets of short windows, `N <= 20000` | 271,110 subsets, 0 violations; CC bound smaller at `s = 3` (72,864 cases) |
| `SLW(5, g)` of `weakest.md` | `r4_slw_templates.py`: exact moment conditions and coprimality, 2000-digit angles | `gamma_N -> 1/2` (s = 4), `-> 13/22` (s = 6) |

## 2. Sections 1-3: the gcd lemma and the smoothness theorem

### 2.1 Inputs and Lemma 2.1

* **(T1).** `two_thirds` (3.6) gives `sum log|c| >= (ramified + collision terms)`. Adding `Q/2` gives the left side
  of (5.1), which is `>= F*(M)` by Thm 5.2 step 1. The upper bound is (2.3) plus the chord bound, (2.5).
  * Only `gcd = 1` and the chord bound `|z_x - z_y| <= c n^(1/4)` are used.
  * `gcd = 1` already forces `n` odd with no inert factor (Lemma 1.1(a),(b)), and it forces `n_0, n_e >= 1` in
    `Psi_p`.
  * So (T1) applies verbatim to `u^B = z/G_B`, whose points are distinct and have gcd 1.
* **(T2).** Step 2 of Theorem 7.5 is `8D(M)/M > 3M log M - K_0 M + 3 + 2 log 2` for every integer `M >= 2`.
  `F*(M) >= D(M)` holds termwise: each inert term is `>= E(M, l+1) log l`, and each `G_p >= (c_p M^2/4 - M/2)`
  on both branches (Lemma 4.3). So (T2) is correct as quoted.
* **Lemma 2.1.** Re-derived. The coefficient of `gamma_B` is `K_s/2 - s(s-2)/8 = s^2/8`, and the reduced chords of
  `u^B` are those of `z`. Exact on 127,410 subsets, including non-normalised `z` (content `2, 3^2, 7^2, 11^2`).

### 2.2 Proposition 2.2

The rescaling is right:
* `|u_x - u_y| <= C_B N^(1/4)/|G_B| = (C_B |G_B|^(-1/2)) N_B^(1/4)`;
* the `gamma_B` coefficient is `s/8 + K_s/4 = s^2/8`;
* (2.2) follows.

The consequences use `(4(s-1)/s) log C_B <= 4 log^+ C` and (T2). The equality case of Remark 2.3(1) is exactly
`|beta|^2 = 2`, which I confirmed both ways on 4,623 exact equality instances.

Remark 2.3(2) is consistent: for `C <= sqrt 2`, `F*(s) >= K_s (log 2)/2` gives
`(4(s-1)/s) log C - 8F*/s^2 <= 0`, so `b'(s) <= W'/s <= 2W'/M` for `s >= M/2`. This recovers the notes' obtuse
Bessel bound.

### 2.3 Theorem 3.1, Corollary 3.2, Theorem 3.3, Propositions 3.4-3.5, Corollary 3.6

All re-derived; no gap.

* **Theorem 3.1, case 2.** Every `x` lies in `U` (`a_xp >= 1`) or in `L` (`a_xp <= e_p - 1`). Each `u_x` with
  `x in U` is divisible by `pi_p`, and each with `x in L` by `conj pi_p`. The larger side has
  `M/2 <= s <= M - 1`. Prop 2.2 applies to the normalised configuration with `C'_B <= C`, and `b'_C <= b_C`.
* **Theorem 3.1, case 1** (content) uses `B = [M]`.
* **Corollary 3.2** is algebra. For `M >= 4` and `log P_+ >= lambda W + 4 log^+ C`, (3.1) gives `M < 2/lambda`.
* **Theorem 3.3(2).**
  * The minimality of `t*` gives `|U*| > M/2`, and `|L*| >= M/2`; both are `>= 2` when `M >= 4`.
  * `L* = [M]` forces `t* = e`, where the bound is trivial.
* **Theorem 3.3(4)** uses `m_xy >= log 2 - 2 log C'`, `C' <= C`, and the fact that each cut weighs at most its
  class.
* **Proposition 3.4.**
  * `Q < MW'/4` for `C <= 1` (from (T1) and `|beta| >= sqrt 2`).
  * The unbalanced weight is `< W' M^(-1/3)`.
  * The balanced cuts at `p` form an interval `[tau_1, tau_2]`, and `conj pi^(tau_2 - tau_1 + 1)` divides
    every `u_x` with `a_xp < tau_1`.
  * `|L| >= M/2 - h >= 2` for `M >= 8` (`M = 8` gives exactly 2).
  * Summing gives `omega > (M/2)(1 - M^(-1/3))^2`.
* **Proposition 3.5.**
  * The constant: `4 * 1.25506 = 5.0202 <= 5.03`.
  * Rosser-Schoenfeld (Cor. 1 to Thm 2, `pi(x) < 1.25506 x/log x`, `x > 1`) is quoted correctly.
* **Corollary 3.6.** The `o(1)` bookkeeping is right: `log M >= loglog R - logloglog R + log eps` from
  `H(eps)`, and `log M <= loglog R` from (T3).

The checks of `r3`/`r5` cover:
* items (1)-(4) of Theorem 3.3, both with the stated `beta'` (built from `C`) and with the stronger `C'`;
* the strict `C <= 1` forms whenever such windows occur;
* Prop 3.5 at `Y = 5, 13, 30, 100, M`.

There were no failures. The summary's near-extremal list drops the hypothesis `C <= 1` from ">M/2 distinct
partitions", "rows differ in >M/4 cuts" and "`omega > M/4`". For general `C` these read
`(M/2)(1 - O(log^+ C/log M))` and similar (minor, §7).

### 2.4 Attribution: the crude forms follow from Cilleruelo-Cordoba

**Quoted theorem (CC).** J. Cilleruelo and A. Cordoba, *Trigonometric polynomials and lattice points*, Proc. Amer.
Math. Soc. 115 (1992) 899-905: on the circle `x^2 + y^2 = rho^2` (`rho^2 in N`), an arc of length
`sqrt2 rho^(1/2 - 1/(4 floor(m/2) + 2))` contains at most `m` lattice points (`m >= 1`).

Provenance of the quote:
* I quote it from memory.
* A web search in this session returned the same statement as Theorem 1.1 of Cilleruelo-Granville, *Close
  lattice points on circles*, Canad. J. Math.
* The full text could not be fetched (DNS blocked).
* The project already cites the qualitative form (`growth/README.md`, `conditional.md`).

**Derivation.**
1. *Division.* Let `s >= 2` points `z_x` of an arc of length `L <= C sqrt R` of `|z| = R` (`R^2 = N`) all be
   divisible by `g in Z[i]`, with `n = Norm g`. The `z_x/g` lie on an arc of length `L/|g|` of the circle of
   radius `rho = R/|g|`, with `rho^2 = N/n`, because division is a similarity.
2. *Apply CC with `m = s - 1`.* An arc holding `s` points must be longer than the CC length, so
   `C rho^(1/2) |g|^(-1/2) >= L/|g| > sqrt2 rho^(1/2 - 1/(4 floor((s-1)/2) + 2))`.
3. *Take logarithms.* With `W = log N`:

   ```text
   s even:  log n < W/s + (4(s-1)/s) log C - (2(s-1)/s) log 2,
   s odd:   log n < (W + 4 s log C - 2 s log 2)/(s + 1).
   ```

**Consequences.**
* For even `s`, the leading term `W/s + (4(s-1)/s) log C` is exactly that of Prop 2.2. For odd `s`, the
  exponent `1/(s+1)` beats Prop 2.2's `1/s`.
* Feeding this into the halving step of Theorem 3.1 (the side `B` with `|B| >= M/2` shares `pi_p` or
  `conj pi_p`, and the content is shared by all points) gives `log q < 2W/M + 4 log^+ C`. That is the crude
  form of Theorem 3.1.
* It also gives Corollary 3.2, the restricted growth statement, and `M < max(4, 2/lambda)`.

**What Prop 2.2 adds over CC.** The `F*` gain `8F*(s)/s^2 = 3 log s - O(1)` in place of
`(2(s-1)/s) log 2`. This yields the `-3 log(M/2) + K_0` refinement, hence
`P_+ <= (log R)^(4/eps - 3 + o(1))` instead of `(log R)^(4/eps + o(1))` under `H(eps)`. It also yields the
self-contained derivation and the systematic use in Sections 3-5.

`r6_cc_gcd.py` checks the two displayed bounds on 271,110 subsets of short arcs, using the actual arc length, and finds no
violation. This checks my derivation, not CC.

So the proofs stand. The labels "proved, new" (table rows Thm 3.1, Cor 3.2) and "the key new tool is the gcd
lemma" should be replaced by an attribution: the crude forms follow from CC via the halving step, and the new
content is the `F*`-refined second-order term and its use.

## 3. Sections 4-5: Plotkin balance, margins, collision stability

**Proposition 4.1.**
* **(a)** `|u_x - u_y|^2 = |g_xy/G|^2 |beta|^2` and `log|g_xy/G|^2 = W' - d_xy` give `m = 2 log|beta| - 2 log C_xy`.
  The lower bound is `two_thirds` Lemma 3.1 with Lemma 3.3(R).
* **(b)** The Gram identities follow from `<f_x, f_y> = W' - 2d_xy` and `|sum f|^2 = sum w (M - 2S)^2 = 4Q`.
  The per-prime cut identity was checked in integers in every window.
* **(c)** Re-derived as follows.
  * `sum log|beta| >= (K/2) log 2 + sum_l E(M, l+1) log l` (depth one, `l = 3 mod 4`).
  * `2E >= M^2/(l+1) - M`, which may be negative and so may be added for all `l <= M`.
  * `A(M) >= S_3(M) - a_3` and `S_3 = (O - X)/2 >= (log M - beta - xi + 1/M)/2`.
  * `theta_3(M) < M log 4`.

  The constant is `(2.1019419 + 2.70189)/2 + 0.1743774 + log 4 = 3.9625877 < 3.9626`.
* **(d)** `mbar <= MW'/(4K) = W'/(2(M-1))`. The lower bound (4.2) holds for all `M >= 2`.
  * The table's simplified form `log M - 7.3 - 2 log^+ C` needs
    `(M/(M-1))(log M - 7.9252) + log 2 > log M - 7.3`, which holds exactly for `M >= 58` (checked), so the
    stated `M >= 120` is safe.
  * The task summary drops the `-2 log^+ C` term and the threshold on `M` (minor).
* **(e)** Upper: (T1). Lower: as in (d). The Markov count uses `log|beta|^2 >= log 2 > 0`. The reading
  `|beta|^2_geo in [M^(1-o(1)), C^2 M^(1/eps + o(1))]` under `H(eps)` is correct.

All of (a)-(e) hold in every one of the 117,892 windows.

**Proposition 5.1.** The decomposition is an identity once two facts are used:
* `v_p(beta) = 0` for pairs on different `p`-levels (`two_thirds` Lemma 3.3), so `sum_pairs v_q(beta) = sum_a kappa_(q^a)`;
* `(K/2) log 2 + sum_q log q sum_a E_(q^a) + Q/2` is literally the left side of `two_thirds` (5.1).

Each excess is nonnegative: Lemma 3.1, pigeonhole, and Theorem 5.2 step 1 there. The bound
`Xi <= (M^2/8)(W'/M - 3 log M + K')` uses `F*(M) >= D(M)`, (T2) and `K log^+ C <= (M^2/2) log^+ C`. I computed
`X_chord`, `X_coll` and `X_prof` from scratch in every window:
* `v_q(beta)` from the factorisation of `|beta|^2`;
* `kappa_(q^a)` by direct count;
* level-wise `E` at `p | N'`.

The identity holds to 30 digits, and `0 <= X <= Xi` holds throughout. The phrase "`(A-3)/3 + O(1/log M)` times
the main term" uses `F*(M) >= (3/8) M^2 log M - O(M^2)`, which is (T2), so it is fine.

**Theorem 5.2.** I re-derived the termwise comparison with `D(M)`.
* For split `p <= M` with `p` not dividing `N'`, the actual term is `>= E(M, p-1) log p`.
* `D` uses `(c_p M^2/4 - M/2) log p` there.
* With `E(M, p-1) >= M^2/(2(p-1)) - M/2` and `c_p <= 1/(p-1)`, the difference is `>= M^2 log p/(4(p-1))`.

The final constant `(W' - 3M log M + K'M)/(2M)` absorbs `4K log^+ C/M^2 < 2 log^+ C`. The Mertens input in the
class `1 mod 4` is quoted (Williams 1974), not proved; it enters only the interpretation, not the inequality.

**The model reading `A >= 4 - o(1)` (Prop 7.1)** is correct and is a small new fact about the class
`A_1^>(Q_0)`, `Q_0 >= 5M/4`. In that class every split `p <= M` is free, so the split gain doubles from
`(M^2/8) log M` (chain lemma) to `(M^2/4) log M` (pigeonhole). This is consistent with every known fake
(`A >= 8`). It concerns the constant, not the growth rate.

## 4. Section 6: the Matveev statement

**Algebra.**
* `u/conj u = eps^2 prod (pi/conj pi)^(2a - e)`, and `sum c = 0` removes `e`.
* `s = prod eps_x^(2c_x)` lies in `{+-1}`.
* `|Pi - 1| <= ||c||_1 Delta`, after centring the angular interval.
* `Pi != 1` when `v != 0`, by unique factorisation, since `pi_p` and `conj pi_p` are non-associate.

`r3` checked `Pi = +- prod (pi_p/conj pi_p)^(2v_p)` in cross-multiplied Gaussian integers for 352,358 random
zero-sum `c`, and `|Pi - 1| <= ||c||_1 (angular width)` for each. Both held.

**Heights.** The polynomial `p X^2 - 2(a^2 - b^2) X + p` is primitive, because `p` does not divide `a^2 - b^2`.
Hence `h(pi/conj pi) = (1/2) log p` and `A_p = max(log p, pi)`.

**Final inequality.** With `n - 1 = |supp v| <= (1 - eta) L/log L`, the right side is
`(1 - eta/2) L + O(L/log L)` against `L - O(1)`. The contradiction is correct. Only the exponential dependence
`c(n) <= c_0^n` is used.

**The quoted input.** It is not stated precisely, and the write-up says "quoted from memory" (minor).

*Standard form (also from memory):* E. M. Matveev, Izv. Math. 64 (2000) 1217-1269, Corollary 2.3, as restated in
Bugeaud, *Linear Forms in Logarithms and Applications* (EMS 2018), Thm 2.2. Let:
* `K` have degree `D`, with `kappa = 1` if `K` is real and `2` otherwise;
* `Lambda = alpha_1^(b_1) ... alpha_n^(b_n) - 1 != 0`;
* `A_j >= max(D h(alpha_j), |log alpha_j|, 0.16)`;
* `B >= max |b_j|`.

Then

```text
log|Lambda| > - min{ kappa^(-1) (e n/2)^kappa 30^(n+3) n^(3.5), 2^(6n+20) } D^2 (1 + log D)(1 + log(nB)) A_1...A_n.
```

The factor `1 + log(nB) <= (1 + log n)(1 + log B)` is absorbed into `c_0^n`. So Prop 6.1 is correct with any
version of this form, and the write-up should state one with its source.

**Scope.** Prop 6.1 is genuine, since it uses the actual `arg pi_p`. Its conclusion is weak:
* the support must exceed `(1-eta) log M/loglog M`;
* the `construct.md` profiles have supports `>= 2 sqrt((b-1)M) - 1`, which is 17.3 at `M = 12` against
  measured minima of 53 and 20;
* the skew double core has supports `>= M/2`.

So it separates nothing from the fakes, as Remark 6.2(3) says. I re-derived the uncertainty bound of Remark
6.2(3):
* unflipped copies of labels in `T`;
* flipped columns of rows `x in supp c` with `a(x)` not in `T`;
* `a` injective except for `a*`;
* `|supp c| |T| >= M` from `H^T H = M I`, then AM-GM.

## 5. Section 7: transfer to the angle models, and the growth check against the round-5 fakes

### 5.1 Proposition 7.1 (transfer)

I re-derived the four inputs.

* **Pair lower bound.** Write the point angles as `psi_x = <2a_x - e, phi> - (pi/2) k_x`. Then
  `|beta_xy| = 2 e^(d_xy/2) |sin((psi_x - psi_y)/2)|` and
  `(psi_x - psi_y)/2 = <a_x - a_y, phi> - psi_1(e_x - e_y)`.
  So the pair case of (R) with `u = 1` is `|beta| >= sqrt2 D_xy`, ramified factor included.
* **Depth compatibility.** Summing over depths needs `rho^(q^a) = rho^(q^(a-1)) mod q^(a-1)`, which is Definition
  1.2 of `construct.md`. Collisions mod `q^a` then imply collisions mod `q^(a-1)`, and
  `sum_pairs log D_xy = sum_q log q sum_a #(collisions mod q^a)`.
* **Pigeonhole.** It uses the conic fibre `(q - chi(q)) q^(a-1)` (Lemma 3.2 of `two_thirds`, `q` not dividing
  `N`).
* **Moduli.** The moduli that `F*(s)` uses have fibre `< s`. Hence `q^a < s` (inert) or
  `q^a < s q/(q-1) <= 5s/4` (split), so they are `<= Q_0` when `Q_0 >= 5M/4`.
* **Sub-models.** These are the sub-profile on `B` with the cuts constant on `B` removed. The residue data and the
  pair data are unchanged, because a fake common factor cancels in collisions and the reduced chords of
  `B`-pairs are unchanged. The chord constant is `C_B e^(-gamma_B/4)`.

So (T1) holds for every sub-model with `F_free >= F*`, and every derivation in Sections 2-5 goes through.

**Gaps (minor).**
1. Prop 5.1 needs `v_q(beta)` and Prop 7.1 defines only `|beta|`. In the model, read `v_q(beta)` as the exponent
   of `q` in `D_xy`; then the decomposition holds with `X_chord = sum (log|beta| - (1/2) log 2 - log D_xy)`.
2. Cor 7.2 at `M` not of the form `q + 1` uses `construct.md` Lemma 1.3 (row deletion). Deleting rows can make
   columns constant, so the profile is no longer normalised.
   * Deleting those columns lowers `W` and keeps every axiom: the monomial gaps concern `phi` only, and the
     residue data can be divided, for each modulus `m`, by one fixed class of `(Z[i]/m)^*` of norm `p^e`. Such
     a class exists because every `p > Q_0`, so `p` does not divide `m`, and dividing by it preserves
     collisions, the conic condition (new norm `N/p^e`) and depth compatibility.
   * The bound `(27 + o(1)) M log M` survives because `q + 1 = (1 + o(1)) M` (primes in progressions).
   * This step should be written out.

### 5.2 Growth check against the fakes (the task's main test)

The write-up claims no growth improvement for general circles. Its only growth statements are the following.

| claimed growth statement | input used | do the `construct.md` / `lower.md` / `exact.md` fakes satisfy the input? |
| --- | --- | --- |
| Cor 3.2: `M = o(log R/loglog R)` if `log P_+(N)/loglog R -> infinity` | a prime factor of `N` larger than every power of `log R` | **no**. In `construct.md` the primes lie in one window above `160 M^2 L/C`, with `log p = (3 + o(1)) log M = O(loglog R)` (the `M = 12` instance has `log p = 16.1`). In general, Thm 3.1 transferred by Prop 7.1 gives `log P_+ <= 2W/M + O(1) = O(log M)` for **every** model in `A_1^>(5M/4)` with `W = O(M log M)`, so no fake of that size can have this input |
| Cor 3.2: `M < max(4, 2/lambda)` if `P_+(N) >= N^lambda C_+^4` | a prime factor `>= N^lambda` | **no**: `log P_+/W <= 2/M + O(1/W)` in every such model |

Both are correct for actual circles. Both use an input that the fakes violate, so they do not contradict the
barrier. As §2.4 shows, both are also consequences of Cilleruelo-Cordoba. Conversely, Prop 7.1 shows that on
`M^(O(1))`-smooth norms no conclusion of Sections 2-5 separates actual clusters from the fakes. The claimed barrier
(Cor 7.2) is therefore correct relative to its stated input class.

That class is the angle-model class `A_1(Q_0)`: point-level residues at moduli `<= Q_0` used through pairs, all
monomial gaps, and Haar-generic angles. It does not cover:
* the prime-level parity law of `sharp.md`;
* residues at all moduli;
* additive relations.

The write-up says this in §7 and §9.1.

### 5.3 The exact illustration on the `M = 12` fake

I re-ran `check_fake_structure.py`; the output is reproduced byte for byte. The profile is that of `construct.md`
§2: Paley-I with `q = 11`, normalised, `b = 8` unflipped copies, and one flipped copy per row with `a* = 1`. The
check is **profile-level only**:
* the angles and residue data of `construct.md` Theorem 2 are not built;
* Prop 2.2 is tested with `C_B := 1` (the arc bound), not the model's sub-cluster chords.

`C_B < 1` makes `b(s)` smaller, and for `s = 2` Prop 2.2 is equivalent to the pair axiom `|beta|^2 >= 2`. So the
table's "max `gamma_B - b(s) = -47.7`" certifies only the weight part. The rest is guaranteed by Prop 7.1's
proof, not by the computation. The table should say so (minor).

### 5.4 "Integrality of the reduced chords is equivalent to the genuineness of the prime angles"

Only one direction is shown. Genuine angles make the model an actual configuration, so its `beta` are integral.

The converse is not proved, and the statement is not even well-defined, since Prop 7.1 defines only `|beta|`. An
example shows the converse fails pointwise. On `|z|^2 = 13` take the two points `sqrt12 +- i`, a fake
`pi = sqrt12 + i` with `e = 1`. The fake gcd is 1 and `beta = 2i` is integral, while the points are not Gaussian
integers. I did not check whether this two-point example satisfies every Lemma 0.1 gap. It only shows that
"equivalent" needs an argument (minor).

## 6. Section 8: the exploitation routes, and the pointer to `weakest.md`

### 6.1 Sections 8.1-8.5

* **Lemma 8.1.**
  * *Coprimality.* At each Gaussian prime only the pair avoiding the minimiser can exceed `G_3`.
  * *Identity.* (8.1) is the triangle identity divided by `G_3`.
  * *Norm bounds.* `log|A_w|^2 = W'/2 - m_xy - gamma_3`; Prop 2.2 at `s = 3` with `F*(3) = (3/2) log 2`;
    `m >= log 2 - 2 log^+ C`.
  * *Checks.* Exact for the identity and the coprimality, and to 40 digits for the norm identity and both
    bounds, on 605,328 random triangles.
* **Prop 8.2.** It is a correct "this route cannot work" remark. The references are quoted as in the standard
  literature, from my memory as well:
  * Beukers-Schlickewei, Acta Arith. 78 (1996): `<= 2^(8(r+1))` solutions of `x + y = 1` in a group of rank `r`;
  * Erdos-Stewart-Tijdeman, Compositio 66 (1988): `exp((4 - o(1)) (s/log s)^(1/2))` solutions for suitable
    `S`.

  The deduction (bounds depending only on `|S| = Theta(M)` cannot see `binom(M,3)` solutions) is right.
* **Prop 8.3.** `a_x + a_y = a_w + a_v` forces `u_x u_y = (unit) u_w u_v`, and the unit is 1 because the angle
  difference is `< pi/2`. Exact on all 14,928 such quadruples found in the scans: `u_x u_y = u_w u_v` exactly.
  Otherwise Prop 6.1 with `||c||_1 = 4` applies.
* **§8.2 (abc).** The computation assumes the three terms coprime "since common factors only lower the height".
  Dividing out a common factor also lowers the radical, so the reduction needs a line. Since the paragraph only
  argues that abc is vacuous, this is cosmetic.
* **§8.3 (smooth numbers near squares).** The algebra `N' - X = |u_x - u_y|^2/2` is correct. The conclusion
  ("no reduction was found") is a report, not a theorem.
* **§8.4** is correctly labelled heuristic. GRH for Hecke L-functions of `Q(i)` is invoked without a precise
  statement. That is acceptable for evidence only.
* **Remark 2.3(3) and §8.5.** The claims that Prop 2.2 "is the complete content of recursion on sub-clusters"
  and that iterating "only re-proves Theorem 3.1" are unproved assertions. Other recursions exist, e.g. applying
  (T1) to `u^B` and keeping its collision terms, as Prop 5.1 does for `B = [M]`. Rephrase as "the natural
  recursion gives exactly Prop 2.2" (minor).

### 6.2 The pointer "by `round8/weakest.md` (draft), five simultaneous forms do"

**The hypothesis.** `weakest.md` Theorem 2.3 assumes `SLW(k, g, kappa)` with `g < 1/4`:

```text
max_i |Lambda_i| >= c^(n+1) prod_j (1 + B_j)^(-kappa) prod_(j in U) A(alpha_j)^(-g)
```

for `k` forms with independent coefficient rows in logarithms of `alpha_j in Q(i)`, `|alpha_j| = 1`,
multiplicatively independent mod roots of unity.

**The counterexample family (vertical templates).** Fix:
* distinct `b_1, ..., b_D`;
* `s` and `k = D - s`;
* independent `sigma_i in Z^D` with `sum_l sigma_il b_l^(2j-1) = 0` for `j <= s`.

For `t` divisible by `lcm(b) P^2` (`P` the product of the primes `<= 2 max b`), put `u_l = t/b_l + i`. The
`u_l` are pairwise coprime and coprime to every `conj u_m`. A common prime would divide the integer
`t/b_l -+ t/b_m = (lcm/b_l -+ lcm/b_m) P^2 T`, but `u_l = i (mod p)` for every `p` dividing `t/b_l`, and the
remaining factor has only primes `<= max b`. `r4` checks the coprimality exactly for the `t` used. Define
`beta_i = prod u_l^(sigma_il)` (conjugates for negative exponents). Then:
* `Lambda_i = Log(beta_i/conj beta_i) = 2i sum_l sigma_il arctan(b_l/t) = O(t^(-(2s+1)))`, because the first `s`
  odd moments vanish;
* `Lambda_i` is an integral linear form in `Log i` and the `Log(pi_p/conj pi_p)` over the Gaussian primes of the
  `u_l`;
* the rows are independent because the `sigma_i` are and the prime sets are disjoint;
* `prod A(alpha_p) <= 2^n prod_l N(u_l) ~ t^(2D)`, with `n = sum omega(N(u_l)) <= (2D + o(1)) log t/loglog t`.

For the coefficients, the bound `B_j <= max|sigma| v_j` holds, where `v_j` is the exponent of the `j`-th prime
and `sum v_j = O(log t)`. By concavity,

```text
sum_j log(1 + B_j) <= n log(1 + max|sigma|) + n log(1 + sum v_j/n) = O(log t * logloglog t/loglog t) = o(log t).
```

The coefficient of `Log i` is `O(n max|sigma| log t)`. So the factors `c^(n+1) prod (1+B_j)^(-kappa) 2^(-gn)`
are `t^(-o(1))`, and SLW fails for every `g < (2s+1)/(2D) = (2s+1)/(2(k+s))`.

**Numbers.**

| `b`, `s` | forms `k` | `gamma` limit | `gamma_N` (at `log t`) along increasing `t` |
| --- | --- | --- | --- |
| `1..9`, `s = 4` | 5 | `1/2` | `0.4929 (41), 0.4953 (62), 0.4968 (89), 0.4980 (145)` |
| `1..11`, `s = 6` | 5 | `13/22 = 0.591` | `0.5835 (49), 0.5857 (70), 0.5872 (98), 0.5886 (153)` |

Here `gamma_N := log(1/max|Lambda|)/sum_l log N(u_l)`.

**Consequence.**
* For fixed `k`, `SLW(k, g, kappa)` is false for every `g < 1`; let `s -> infinity`.
* In particular `SLW(5, g)` with `g < 1/4` is false.
* So `weakest.md` Theorem 2.3 is vacuous, and the Dirichlet heuristic of its Remark 2.1(3) is beaten by
  polynomial families.

This is the construction of `referee_weakest_checks/slw_template_counterexample.py`, re-derived independently
(`r4_slw_templates.py`, no factoring needed). `weakest.md`'s own referee report is not yet written, so I flag
this here.

**Required correction.** In Section 8.6, delete "five simultaneous forms do". Instead say that the simultaneous
hypotheses `SLW(k, g)` of `weakest.md` are false for fixed `k`, and that no simultaneous-approximation hypothesis
known to imply growth is currently on the table.

The rest of §8.6 is a correct reformulation of Theorem 3.1:
* restriction to Gaussian primes of norm `<= C_+^4 (log R)^(4/eps)`;
* coefficients `|v_p| <= ||c||_1 e_p`;
* `e_p log p < (8/eps) loglog R + 8 log^+ C`;
* supports `> (1-eta) log M/loglog M` by Prop 6.1.

**Inconsistency (minor).** The same paragraph calls `lw.md` "refereed ... survives as a barrier", while Remark
6.2(3) calls it "(draft)". `referee_lw.md` says it survives for the **max form with `Log i` in the list** only,
not for the product form of Lang/Waldschmidt. §8.6's "in its max form" should add "with `i` among the numbers".

## 7. List of issues

**Major (neither affects the correctness of Sections 2-7):**

1. **Attribution of the crude smoothness results (§2.4).** The crude gcd lemma, the crude Theorem 3.1
   (`P_+(N) < N^(2/M) C_+^4`), Corollary 3.2, the restricted growth corollary and `M < max(4, 2/lambda)` follow in
   two lines from the Cilleruelo-Cordoba theorem (1992) applied after dividing by the shared Gaussian factor.
   * For odd `s`, CC gives the better exponent `1/(s+1)`.
   * On 271,110 window subsets, the CC-derived bound is the smaller one in 72,864 cases (all at `s = 3`).
   * The table labels "proved, new" and "the key new tool is the gcd lemma" overstate the novelty.
   * New are: the `F*` second-order gain `-3 log(M/2) + K_0` and the resulting exponent `4/eps - 3`; the
     systematic sub-cluster use in Thm 3.3 and Props 3.4-3.5; Sections 4-5; Prop 6.1; and the transfer
     (Prop 7.1), including the model fact `A >= 4 - o(1)` on `A_1^>(5M/4)`.
2. **Section 8.6 cites a false hypothesis as a route to growth (§6.2).** `SLW(k, g)` of `weakest.md` is false for
   every fixed `k` and every `g < 1` (vertical templates). The sentence "five simultaneous forms do" must go.

**Minor:**

3. Matveev's theorem is "quoted from memory" without its exact constant or form. It should be stated (§4). The
   argument is unaffected.
4. The task summary drops qualifiers that the write-up has:
   * "`omega(N') > M/4`", "`> M/2` distinct partitions" and "rows differ in `> M/4` cuts" need `C <= 1`;
     otherwise read `(1 - O(log^+ C/log M))` times these;
   * the average-margin lower bound needs `-2 log^+ C` and `M >= 58`.
5. §7: "integrality of the reduced chords is equivalent to genuineness of the prime angles" proves only one
   direction, and in the models `beta` is defined only by its modulus (§5.4).
6. The `M = 12` illustration is profile-level. Prop 2.2 is tested with `C_B = 1`, not with model chords. The
   table should say "weight part of Prop 2.2" (§5.3).
7. Prop 7.1 needs a definition of `v_q(beta)` in models (use `D_xy`) for Prop 5.1. Cor 7.2 for `M != q+1`
   needs the row-deletion and renormalisation step written out (§5.1).
8. Remark 2.3(3) and §8.5 ("complete content of recursion on sub-clusters", "only re-proves Theorem 3.1") are
   unproved assertions.
9. `lw.md` is labelled both "draft" (Remark 6.2(3)) and "refereed" (§8.6). The barrier is for the max form with
   `Log i` included (`referee_lw.md`).
10. §8.2 (abc): the reduction to coprime terms needs a line, since dividing by a common factor lowers the radical
    too.
11. Cosmetic: the threshold `M >= 120` for the simplified margin bound can be `M >= 58`.

**Checked and found correct:**
* (T1)-(T3) as quoted;
* Lemma 2.1, Prop 2.2 and Remark 2.3(1)-(2);
* Thm 3.1, Cor 3.2, Thm 3.3, Props 3.4-3.5 and Cor 3.6;
* Prop 4.1, including the constants `3.9626` and `7.9252`;
* Prop 5.1, (5.3) and Thm 5.2;
* Prop 6.1 (given Matveev) and Remark 6.2(3);
* Prop 7.1 and Cor 7.2 (given `construct.md`);
* Lemma 8.1, Props 8.2-8.3;
* all numbers in the Section 7 table, which I reproduced.

## 8. Corrected statement

**No growth improvement for general circles is proved, and none is claimed.** Let `z_1, ..., z_M` be distinct
lattice points with `|z_x - z_y| <= C sqrt R` on `x^2 + y^2 = R^2 = N`, with `W = log N`, normalised as in
`two_thirds` (`G, N', W', C' = C|G|^(-1/2)`). Proved unconditionally:

1. **gcd lemma (Prop 2.2).** For every subset `B` of size `s >= 2`,
   `log Norm(gcd_B z) <= W/s + (4(s-1)/s) log C_B - 8F*(s)/s^2 < W/s - 3 log s + K_0 + 4 log^+ C`, where
   `K_0 < 13.9092`.
   * Equality holds iff `s = 2` and `|beta|^2 = 2`.
   * *Attribution:* the weaker bound `W/s + (4(s-1)/s) log C - (2(s-1)/s) log 2` for even `s`, and
     `(W + 4s log C - 2s log 2)/(s+1)` for odd `s`, follows from Cilleruelo-Cordoba (1992) by dividing by the
     gcd. The `F*` gain is what is new.
2. **Smoothness (Thm 3.1, Cor 3.2).** For `M >= 4`, every prime `q | N` satisfies
   `log q <= max_(ceil(M/2) <= s <= M) b_C(s) <= 2W/M - 3 log(M/2) + K_0 + 4 log^+ C`.
   * The crude form `q < N^(2/M) C_+^4` and its corollaries follow equally from Cilleruelo-Cordoba plus the
     halving step:
     * `M < 2W/(log P_+ - 4 log^+ C)`;
     * `M = o(log R/loglog R)` whenever `log P_+(N)/loglog R -> infinity` (`C` fixed);
     * `M < max(4, 2/lambda)` if `P_+(N) >= N^lambda C_+^4`.
   * These are restricted theorems. Their input (a prime factor of `N` larger than every power of `log R`) is
     violated by every fake of size `W = O(M log M)`. In `construct.md`, `log p = (3 + o(1)) log M`. In general,
     every model in `A_1^>(5M/4)` has `log P_+ <= 2W/M + O(1)` by Prop 7.1.
3. **Near-extremal structure (Cor 3.6, Thm 3.3, Props 3.4-3.5).** Fix `eps in (0,1]` and `C`, and assume
   `M >= eps log R/loglog R` with `R >= R_0(eps, C)`.
   * `P_+(N) <= M^(4/eps - 3 + o(1))`, and `p^e || N'` implies `p^e < C_+^8 (log R)^(8/eps)`.
   * `(M/4)(1 - o(1)) < omega(N') <= Omega(N') <= (2/eps + o(1)) M`. For `C <= 1`: `omega(N') > M/4`,
     `omega(N') >= (M/2)(1 - M^(-1/3))^2` (`M >= 8`), more than `M/2` distinct cut partitions, and any two rows
     separated by more than `M/4` cuts.
   * The primes below `M^(1-delta)` carry an `O(M^(-delta)/log M)` fraction of `W'`.
4. **Profile and collision stability (Props 4.1, 5.1, Thm 5.2).** Let `A = W'/(M log M)`.
   * `(M/(M-1))(log M - 7.9252) + log 2 - 2 log^+ C < mbar <= W'/(2(M-1))`.
   * `Q < MW'/4 - (M^2/2) log M + 3.9626 M^2 + 2K log^+ C`.
   * Under `H(eps)`, the geometric mean of `|beta|^2` lies in `[M^(1-o(1)), C^2 M^(1/eps + o(1))]`.
   * Each of `X_chord`, `X_coll`, `X_prof` is `<= (M^2/8)(W'/M - 3 log M + K')`, i.e. at most
     `(A - 3)/3 + O(1/log M)` times `F*(M)`.
   * `sum_(p <= M, p = 1 (4), p not dividing N') log p/(p-1) <= ((A-3) log M + K')/2`.
5. **Genuine input (Prop 6.1, given Matveev's theorem in a precisely stated form).** In the near-extremal regime,
   for `R >= R_1(eps, eta, C)`, every zero-sum `c` with `||c||_1 <= exp((log R)^(eta/2))` and `v = sum c_x a_x != 0`
   has `|supp v| > (1 - eta) loglog R/logloglog R`.
6. **Barrier (Prop 7.1, Cor 7.2).** All of 1-4 (with `F*` replaced by `F_free`, and with `A >= 4 - o(1)` in
   place of Thm 5.2) hold in the angle-model class `A_1^>(Q_0)`, `Q_0 >= 5M/4`.
   * The refereed `construct.md` models (`W <= (27 + o(1)) M log M + 9M log(1/C)`, `C <= 1`) lie in this class.
     Row-deleted sub-profiles must be renormalised for `M != q+1`.
   * Their profiles satisfy the conclusion of 5.
   * Hence no argument whose inputs are the `A_1(Q_0)` axioms (point-level residues at moduli `<= Q_0` through
     pairs, all monomial gaps, Haar-generic angles) plus 1-5 can prove `M = o(log R/loglog R)`.
   * The class does not include the prime-level parity law, residues at all moduli, or additive relations.
7. **Section 8.**
   * The exploitation routes (S-unit counting, small-`S` relations = parallel chords, abc, smooth numbers near
     squares) give no growth, for the stated reasons. GRH counting is heuristic.
   * The sentence "five simultaneous forms do" is withdrawn. The `SLW(k, g)` hypotheses of `weakest.md` are
     false for every fixed `k` and `g < 1`.

## 9. Files

All files are in `round8/referee_stability_checks/`. They are pure Python 3 with integers, plus `mpmath`; `r1` imports `stab_lib` only for the
cross-check, and `r4` uses `sympy` for the nullspace. Each output ends with `ALL CHECKS PASSED: True`.

| script | what it checks | output |
| --- | --- | --- |
| `rlib.py` | independent toolkit: Gaussian arithmetic, lattice points, levels, exact `8F*(s)` exponent vectors (DP and brute force) | (library) |
| `r1_fstar.py 250` | `Psi_p` DP = brute force (`s <= 15`); `F*` = `stab_lib.Fstar` (`s <= 80`); (T2) with `+3+2log2` and (3.0) for `s, M <= 250` | `r1_fstar_out.txt` |
| `r2_gcd_subsets.py 6000 7` | Lemma 2.1 and Prop 2.2 in exact integer form on 127,410 arbitrary subsets of ten circles with content; the equality case | `r2_gcd_subsets_out.txt` |
| `r3_window_scan.py NMAX CMAX R2MIN MMAX` | Thm 3.1 (exact) and crude forms; Thm 3.3(1)-(4) (stated and `C'` forms); Prop 3.5; Prop 4.1(a)-(e); Prop 5.1 decomposition and bounds; Thm 5.2; algebra of Prop 6.1; Lemma 8.1; Prop 8.3. Runs: `3e6 6 96 10`, `1e9 10 384 14`, `1e11 40 1024 22`, `1e12 14 1536 18` (17,648 + 42,200 + 49,184 + 8,860 windows) | `r3_window_scan_*_out.txt` |
| `r4_slw_templates.py` | vertical templates refuting `SLW(5, g)` for `g < 1/2` (`s = 4`) and `g < 13/22` (`s = 6`); exact moment conditions and coprimality | `r4_slw_templates_out.txt` |
| `r5_smooth_exhaustive.py 5000000 8` | Thm 3.1, exact refined and crude forms, on **every** window with `M >= 4`, `C <= 8`, `N <= 5*10^6` (435,196 windows) | `r5_smooth_exhaustive_5e6_C8_out.txt` |
| `r6_cc_gcd.py 20000 6` | the gcd bounds derived from Cilleruelo-Cordoba (§2.4) on 271,110 subsets of short arcs: 0 violations. The CC bound is the smaller one in 72,864 cases, all at `s = 3`; Prop 2.2 is smaller at every even `s` and at `s = 5, 7` | `r6_cc_gcd_out.txt` |

The tightest instance of Theorem 3.1 found (`r5`) has ratio `log q / max_s b_C(s) = 0.994`: `N = 3301381 = 61 * 54121`,
`M = 4`, `C = 7.98`, `q = 54121`. So (3.1) is essentially sharp at `M = 4`, as the write-up says.
