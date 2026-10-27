# Referee report: `rigSLtheory.md` (round 15)

Referee: adversarial, independent code (`referee_rigSLtheory_checks/`, nothing imported from the author or the
repository; the author's `classify.py` was *read* once to understand a table convention, see I5). Labels as in the
write-up: *proved*, *exact* (integer/rational arithmetic, complete enumeration), *evidence*.

## 0. Verdict

**Survives as a restricted theorem, after corrections. No growth improvement, no uniform bound; the author says so
and is right.** `M <= (2/3 + eps) log R/loglog R` stands; P3-all\* stays uncapped in regime (1.4); Rig^SL is
undecided.

* **Re-derived and correct:** Lemma 2.1 (J0 closed form OPP/ADJ/MAX), Cor. 2.2 (four-unit laws), Theorem A (C/A0/A1),
  Lemma B0, Lemma B1, Theorem B (pinning), Cor. B2, inequality (2.1) of Theorem C, Theorem P (pinning closure), the
  transfer Theorem T as a conditional statement. Theorem S(a) is `stability.md` Thm 3.3(3) (round 8, refereed);
  S(b) is a plausible input check (not independently re-derived down to `two_thirds`).
* **Every finite claim I tested reproduces exactly** with my own code: forced sets (`K = 4`, all four prime sets,
  identical), one-modulus PU/ANG/SL counts, the A0/A1 profile of 83/163, triple generation `[Z^4 : Lambda_T] = 1` at
  every forced `q_0 >= 7` (`K = 4, 5`, four sets), the one-modulus SL segments and first flexible moduli, lifted
  all-level witnesses with exact `(1.2)` at `phi_0`, strictness data, and TC certificates for `r = 16, 32, 64`.
* **Two claims are wrong as stated (major, but neither propagates to a headline number):**
  1. **Theorem C's reading** "an extra beyond the cofactor slack on a monomial with saturated opposite unit forces a
     dropped collision at a *second* modulus" / "deviations do not come alone" is false: the forced drop can be at
     `q_0` itself. Explicit counterexample in §3 (108 instances for `first4`, `K = 4`). Inequality (2.1) itself is
     correct.
  2. **Cor. S1 / summary** "every `q <= 2M` is forced already at norm 2" is proved only for `C <= 1`. For actual
     clusters that case is empty (an arc of length `sqrt 2 R^(1/2)` holds at most 2 lattice points, so `M >= 3`
     needs `C > sqrt 2`). The clause "(or `M >= M_0(C)`)" is false as derived. The statement that matters,
     `Q_eta ⊆ F_2^ex` for `M >= M_0(C)`, survives.
* Minor issues are listed in §8.

## 1. Class and quantifiers

* **Class.** The SL tests use genuine labels at all levels on `F' = F_K^ex \ {q_0}`, a parity-respecting level-1
  label at `q_0`, and one angle vector for all `v in V_K` and all four units. That is a relaxation of Rig^SL(K)
  (Lemma R(2)), so the rigidity statements transfer. Every flexibility witness is a full system: genuine on `F'`,
  a level-2 lift at `q_0` with no level-2 collisions, collision-free labels on the free moduli (these exist by the
  definition of "free"), and the explicit angle `phi = phi_0`, checked exactly as `m_(v,s) <= |D_s(Z_v)|`. This is
  the SL/P3-all\* truncation, not PU. **Confirmed.**
* **Symmetry group.** "Outside every `eps`-class mod `n`" (with `u in (Z/n)^*`) implies outside every
  `(eps, t, u)`-class: the associate twist `(m/2) t` vanishes mod `n = m/4`. **Confirmed.**
* **Quantifiers.** Theorem B is stated for every label `c'` and every `phi`. It is vacuous when no `phi` exists, and
  the proof handles that case correctly. Theorem P quantifies over members of `F^all_X` that are strict on `S` at
  **all levels `q^a <= X`**, and TC is defined with `q^a <= X`. **Correct.**

## 2. Re-derivation of §2.1-2.4 (theory)

* **Lemma 2.1.** `B` is the image of the four squared conditions. `B` is convex, so it meets the unit circle iff
  `B != ∅` and `min_B |z| <= 1 <= max_B |z|`. Here `d(0,[l,h]) = max(0, l, -h)` on a non-empty interval, which gives
  (ADJ), and `max |z|^2` is attained at a corner, which gives (MAX). (OPP) is non-emptiness. **Correct.**
* **Cor. 2.2.** I re-did the identities `m_0 m_2 <= (N/2)|sin 2x|`, `m_1 m_3 <= N|cos 2x|` (using
  `sin(x - 3pi/4) = -cos(x - pi/4)`), the product law `<= (N^2/4)|sin 4x|`, the PU product `2N^2`, and
  `Im(Z^4) = 4 D_0 D_1 D_2 D_3`. **Correct.**
* **Theorem A.** (1) holds by (1.1). (2) A0 gives `a_s' + a_(s'+2) > 1`, so OPP fails. (3) The OPP and PU
  inequalities are as stated; (3) claims only OPP, and the ADJ/MAX conditions are not part of the claim. (4) The
  chain SL-flex ⇒ ANG-flex ⇒ (no A1) only-C ⇒ SL at `phi_0` is valid. **Correct.**
* **Lemma B0.** This is unique factorisation in `Z[i]`: the `pi_j` and `conj pi_j` are pairwise non-associate, since
  the `p_j` are distinct and `p_j = 1 mod 4`. It holds verbatim for the twisted `pi'_j = i^(t_j) pi_j^(eps_j)`, whose
  arguments differ from `eps_j arg pi_j` by multiples of `pi/2`. **Correct.**
* **Lemma B1 and Theorem B.** In (1), a label collision at an odd unit on `v in P0` is A0. Otherwise the odd demands
  are genuine, and (1.1) puts the odd pair at full budget. (2) Write `<y, eta> = tau_y (pi - 4x_0(y)) mod 2pi`. The
  two branches cannot coincide for `y != 0`: that would need `4x_0(y) = pi mod 2pi`, which Lemma B0 excludes. Take a
  triple `(v, w, v + w)`. Lemma B0 gives `tau_v v + tau_w w = tau_u (v + w)` and `tau_u - tau_v - tau_w` even. I
  checked all 8 patterns for `w = lambda v` (`lambda != 0, -1`): only the constant patterns survive, and parity then
  forces `tau = 0`. (3) A shift by `pi Z^r` changes each `<v, phi>` by a multiple of `pi`. **Correct.**
* **Cor. B2.** If every odd collision is retained on `U^F_1 = U^F_3 = 1`, then `m' >= |D|`. Cor. 2.2(1) then forces
  equality, which pins `phi`. **Correct.**

## 3. Theorem C: the inequality is right, the "second modulus" reading is false (major)

(2.1) is Cor. 2.2(1) with `m'_t = g_t E_t/delta_t` and `|D_t| = U^F_t g_t`. Its special case is also correct: if
`U_(s+2) = 1` and `E_s > U_s delta_s`, then `delta_(s+2) > 1`, i.e. a genuine `(v, s+2)` collision is dropped **at some
`q_1 in F`**.

The text then says more, in §0 ("forces a dropped collision at a second modulus") and in §2.5 ("forces a deviation
at some `q_1 in F \ {q_0}`", "deviations do not come alone"). **This does not follow.** `U^F` is the cofactor
relative to **all** of `F`, `q_0` included. Suppose `q_0 | D_(s+2)(Z_v)` genuinely. A label at `q_0` that moves this
collision from unit `s+2` to unit `s` creates the extra `E_s = q_0` and the drop `delta_(s+2) = q_0^a` at the same
modulus. (2.1), and even the full J0 test for `v`, can then hold with no deviation anywhere else.

**Counterexample (exact, `thmC_search.py`).** `first4`, `K = 4`, `F = F_4^ex`, `v = (-3, 1, 0, 0)`,
`Z_v = (1-2i)^3 (2+3i) = -28 - 29i`, `N_v = 1625`, `(D_0, D_1, D_2, D_3) = (-29, 1, -28, -57)`, `q_0 = 19`, `s = 1`.
* The opposite unit is saturated: `D_3 = -3 * 19` with `3, 19 in F`, so `U^F_3 = 1`.
* The extra exceeds the cofactor slack: `U^F_1 = 1 < 19 = E_1`.
* Move the 19-collision from unit 3 to unit 1. Then `m_1 = 19`, `m_3 = 3`, and `m_0 = 1`, `m_2 = 7` (29 divides `N`;
  `7 in F`). OPP holds: `361 + 9 = 370 <= 2N = 3250`. ADJ and MAX hold as well.

So (2.1) holds and the only drop is at `q_0 = 19`. The script finds 108 such `(v, s, q_0)` for `first4` and 80 for
`rand4` (`K = 4`). The A1 collisions at `rand4`/31 and `bal4`/41 are not of this kind (`thmC_check.py`: 0
instances), so the pinning data do not depend on this point.

*Corrected statement.* The drop is at some `q_1 in F` dividing `D_(s+2)(Z_v)`. It is a **second** modulus
(`q_1 != q_0`) whenever `q_0` does not divide `D_(s+2)(Z_v)`, i.e. when `U_(s+2)` is computed relative to
`F \ {q_0}` (as in Theorem A). The "transport/propagation" heuristic of §2.5 and §6.1 must carry this hypothesis.

## 4. Theorem P (pinning closure): proved; scope as stated

**Proof re-derived.**
* Let `v in P0'(S)`. Strictness up to `(eps, t)` with `u = 1 mod 4` gives a collision of `(v, 1)` at every
  `q^a || D_1(Z'_v)` and of `(v, 3)` at every `q^a || D_3(Z'_v)`. These are allowed: `D_1` and `D_3` are odd and
  coprime, since `gcd(X, Y) = 1` (`Z_v` contains no `pi_j conj pi_j`). Hence `m_(v,1) >= |D_1|` and
  `m_(v,3) >= |D_3|`.
* (e-all) with Cor. 2.2(1) forces equality, and Lemma B1 pins `sin 2<v, phi>`.
* Theorem B(2), with Lemma B0 for the `pi'_j`, gives `phi = phi'_0 + pi n`. The shift contributes the sign
  `(-1)^(<e, n>)`, which is **independent of `x`**, so the points are `+- i^(-k_x) w_x`.
* Distinct rows give non-associate `w_x`, and the arc has length `R Delta = C sqrt R`.
* The scaling must be `u = 1 mod 4`. With `u = 3 mod 4` the odd units swap and the pin changes. This is consistent
  with the round-14 symmetry group.

TC symmetry: `Z'_v = i^(<t,v>) Z_(eps∘v)`, so `Re(Z'^2_v) = +- Re(Z^2_(eps∘v))` and `P0'_(eps,t) = eps∘P0'`.
**Correct.**

**Computation (independent, `tc.py`, exact).** My own incremental integer echelon form gives a generating set of
additive triples, and hence an exact certificate of TC:
* all moduli charged (`gcd(Re Z_v^2, N) = 1`, triples inside `V_3`): first `r = 4, 8, 16, 32, 64` split primes;
* `Re(Z_v^2)` `Y`-smooth, triples `v, w in V_2`, `v + w in V_4`: `r = 16` at `Y = 64`, `r = 32` at `Y = 128`, and
  `r = 64` at `Y = 256` and at `Y = 128 = 2r`.

All have index 1. The author's TC claims are confirmed; I did not check the negative claims ("fails at `Y <= 2r`").

**Caveats (minor).**
1. The all-moduli computation silently needs `X >= max |Re Z_v^2|` over the triples used (`X >= p_max^6` for
   `V_3`), because `P0'` requires `q^a <= X`. In regime (1.4) this holds once `A >= 6 h_0`, but it should be stated.
2. TC is a property of the prime set, and the fake chooses its primes. "P3-G collapses" is proved **on prime sets
   with TC**. TC is verified for the first `r <= 64` split primes and the three loose 4-sets; genericity is
   evidence.
3. Theorem P caps nothing. It identifies the freedom a fake must use: non-strictness on every TC family. Correctly
   stated in the write-up.

## 5. Theorem S (scale): correct, but the `C <= 1` forms are empty for actual clusters (minor), and Cor. S1 overreaches (major)

* **(a)** is `stability.md` Thm 3.3(3) and Prop. 3.4 verbatim. I re-read the proof: it uses no near-extremality.
  However, **by Cilleruelo-Córdoba an arc of length `sqrt 2 R^(1/2)` holds at most two lattice points**, so for
  actual clusters with `M >= 3` one has `C > sqrt 2`. The statements "`r > M/4` for `C <= 1`, any `M >= 4`" and
  "`r >= (M/2)(1 - M^(-1/3))^2`" are therefore empty for actual clusters. They have content only for members (fakes)
  with `C <= 1`. For actual clusters the operative bound is
  `r > (M/4)(1 + 2M log^+ C/W')^(-1) = (M/4)(1 - O_C(1/log M))`. That bound is correct and answers round 14's
  "`r >= M^gamma`?" with `gamma = 1 - o(1)`. For actual clusters this is a **record correction** (round 8 already
  had it), not a new theorem.
* **(b) Members.** The four inputs of `stability.md` Prop. 7.1 are re-checked for sub-models. I verified two points.
  First, the level residues of the sub-model differ from the member's by an `x`-independent factor
  `r^alpha conj(r)^(e - beta)`, so pair collisions are unchanged. Second, a column constant on `B` turns `p_j` into a
  modulus not dividing `N_B` whose pair collisions are exactly the same-level collisions charged in (e) of `A_X`. The
  argument then rests on `lower.md` Prop. 1.3 and `two_thirds` (T1) using **only pair collisions** at such primes. I
  did not re-derive that from `two_thirds`. **Plausible; unrefereed at that depth.**
* **Cor. S1 (major as a statement, no downstream damage).** The pigeonhole step (`2r + 1` points of `S_1` in
  `Z/n_q`) is correct. Applying it to all `q <= 2M` needs `n_q <= (2M+1)/4 < 2r + 1`, i.e. `r > M/4 - 3/8`. Theorem S
  gives this only for `C <= 1`. For `C > 1`, `r` is only `> (M/4)(1 - O_C(1/log M)) = M/4 - Theta_C(M/log M)`. So
  "**(or `M >= M_0(C)`)**, every `q <= 2M` not dividing `2N` is forced for norm 2" is **not proved**. The pigeonhole
  covers only `q < 2M (1 + 2M log^+ C/W')^(-1) + 2`.
  - What survives (proved): every `q <= M`, so `Q_eta`, lies in `F_2^ex` once `M >= M_0(C)`, since
    `(M+1)/4 < 2r + 1`.
  - The summary line "every `q <= 2M` is forced already at norm 2" must be restricted accordingly (`C <= 1`, or
    `q <= (2 - c_C/log M) M`).
* **Cor. S2.** (1) is fine for the window. (2) ("dense regime `n_q <= r(1 + o(1))` on `Q_eta`") is correct. (3) is a
  heuristic and is labelled as one.

## 6. Independent computations (exact)

All code is my own (`referee_rigSLtheory_checks/`):
* `rlib.py`: Gaussian arithmetic; labels by discrete log in `T_q`, normalised so that `g^(m/4) = i`;
  `pi_j = a + bi` with `0 < a < b`.
* `enum.c`: complete column-wise DFS over parity-respecting level-1 labels, with per-`v` forbidden-unit masks and
  class membership mod `n`.
* `onemod.py`: PU masks `q_0^2 g_s^2 nu_s^2 > N`, the cofactor mask `q_0 g_s > |D_s|`, A0/A1 typing by `U`, exact J0
  via Lemma 2.1 in rationals, and the triple lattice by integer echelon form.
* `lift.py`, `strict.py`, `strict_lift.py`, `tgonly.py`, `tc.py`, `thmC_*.py`.

**Sanity check.** Genuine label collisions agree with `q | D_s(Z_v)` in 23,040 checks (`first4`, 18 moduli,
`K = 4`). The parity `L_j = [(p_j/q) = -1] mod 2` is asserted.

**Results.**
1. **Forced sets, `K = 4`** (`forced.py`, exhaustive DFS per `q <= 1500`): identical to the author's/round-14 files
   for `first4` (40, max 199), `rand4` (40, max 199), `bal4` (40, max 271) and `last4` (44, max 263). For `K = 5`, see
   item 7.
2. **83/163.** `first4`, `K = 4`, `q_0 = 83`: PU 36, of which 24 are outside every `eps`-class. All 24 carry 2 or 4 A0
   collisions and none carries an A1 collision; ANG and SL outside are 0. `K = 5`, `q_0 = 163`: PU 160/80, every
   outside label has exactly 2 A0 collisions and no A1; SL 0. **Confirmed.** The `q_0 = 107` row also matches
   (216/144 PU, 108 SL outside, 36 die by A0).
3. **Triple generation.** `[Z^4 : Lambda_T] = 1` at every forced `q_0 >= 7` for all four sets at `K = 4` and `K = 5`.
   The exceptions are `q_0 = 3` (index 2), `q_0 = 5` (index 2, or rank `< 4` for `last4` at `K = 4`). **Confirmed.**
   The write-up's "`|P0|` between 18 and 81" is wrong: the range is 18-164 at `K = 4` and 28-338 at `K = 5`
   (harmless).
4. **One-modulus SL segments** (by Theorem B(3), the exact count is the cofactor count):

   | set, `K` | SL-rigid (outside every `eps`-class) | first SL-flexible `q_0` (count) | status |
   |---|---|---|---|
   | `first4`, 4 | `[3, 103]` | 107 (108) | confirmed |
   | `first4`, 5 | `[3, 163]` | 167 (24); also 173 (42); 179 rigid again | confirmed |
   | `rand4`, 4 | `[3, 47]` | 59 (752) | confirmed |
   | `bal4`, 4 | `[3, 41]` | 43 (40) | confirmed |
   | `last4`, 4 | `[3, 29]` | 31 (24) | confirmed |
   | `rand4`, 5 | `[3, 59]` | 67 (16) | **new (referee)** |
   | `bal4`, 5 | `[3, 61]` | 67 (16) | **new (referee)** |
   | `last4`, 5 | `[3, 41]` | 43 (60) | **new (referee)**; matches round 14's `<= 41` |

   `rand4`/31 and `bal4`/41 have SL outside 0 with index 1, so they are SL-rigid in integer arithmetic. **Confirmed.**
5. **Witnesses (`lift.py`).** I took 20 (or all, if fewer) SL labels outside every `eps`-class at each first flexible
   modulus: 107, 59, 43, 31 (`K = 4`) and 167, 67, 67, 43 (`K = 5`). Each was re-checked directly as
   `m_(v,s) = g_s q_0^[collision] <= |D_s(Z_v)|` for all `v in V_K` and all `s`, re-checked for class membership by
   brute force over `(u, eps)`, and lifted to level 2 with no level-2 collision. **All pass.** This gives:
   * `Q_core^SL(rand4, 4) <= 47`, `Q_core^SL(bal4, 4) <= 41`, `Q_core^SL(last4, 4) <= 29`;
   * `Q_core^SL(first4, 5) <= 163` (with 167 forced at `K = 5`, which I checked by exhaustive DFS);
   * new: `Q_core^SL(rand4, 5) <= 59`, `Q_core^SL(bal4, 5) <= 61`, `Q_core^SL(last4, 5) <= 41`.
6. **Strictness (`strict.py`, `strict_lift.py`).** `K = 4`: every class survivor is `(1,0)`-strict at every forced
   `q_0 in [7, 101]` except 79 (32 survivors / 32 strict for some `(eps, t)` / 16 `(1,0)`-strict), 89 (80 / 40 / 20)
   and 103 (48 / 48 / 24). `K = 5`: all `(1,0)`-strict for `q_0 <= 151`; 157 (48 / 48 / 24) and 163 (80 / 80 / 40).
   **Confirmed exactly.** Three non-strict 89-labels were lifted to level 2 without level-2 collisions, so these are
   genuine all-level SL systems at `phi_0`.
   - Wording: at 79, 103, 157 and 163 the non-`(1,0)` survivors are strict for *some* `(eps, t)`; that may be a
     local conjugation as well as a local twist.
7. **`K = 5` forced sets.** The `K = 5` scans used the round-14 forced files (copied by the author). Independently,
   167 (`first4`) is forced at `K = 5` (exhaustive DFS). The full recomputation is reported in §6.1.

## 7. Implications for growth and the barrier

* **Growth: none, and the write-up does not claim any. Agreed.**
  * Theorem T is conditional (RigS + H^gen + L1). For actual clusters RigS holds trivially with `(eps, t) = (1, 0)`,
    so its "in particular" for actual clusters is conditional on H^gen + L1 alone. Those are open statements about
    genuine Gaussian integers (L7-type).
  * Theorem P turns strict members into actual clusters; it bounds nothing.
  * Theorem S is a size statement.
* **Barrier map.**
  * Correct:
    - the class of all new finite statements is SL (truncated P3-all\*), not PU;
    - the loose-set one-modulus segments are upper bounds on the joint SL cores, not lower bounds (the write-up
      says so);
    - no Rig^SL joint core is decided.
  * The real gains for the barrier map:
    - Theorem B: pinning makes the one-modulus SL problem angle-free;
    - Theorem P: the exact freedom a P3-all\* fake must use, conditional on TC of its prime set;
    - Theorem S: the "`r >= M^gamma`" worry of round 14 disappears, for actual clusters through a record correction
      from round 8.
  * The "propagation" mechanism of §2.5/§6.1 is weaker than stated (§3). A single-modulus swap of a genuine
    collision to the opposite unit is allowed by (2.1). In joint SL it is limited only by pinning (Cor. B2), not by
    transport.
* **Outcome table (§4.2).** The logic is right, given the corrections:
  * "strict coherence on `[3, Y]` with TC ⇒ P3-all\* exact" is valid under Theorem P **for prime sets with TC**;
  * "stall ⇒ C1 seeds, no cap" is right;
  * none of the rows moves the growth rate.

## 8. Minor issues

* **I1.** Theorem S: the `C <= 1` forms are empty for actual clusters (§5). For actual clusters it is round 8's
  theorem, not new; the "no near-extremality" point was already in `stability.md`.
* **I2.** Cor. S1: the "`(or M >= M_0(C))`" range claim. Listed as major in §0 because it is in the result summary;
  it does not propagate.
* **I3.** TC "all moduli" needs `X >= max |Re Z_v^2|` on the triples (§4). Generic TC is evidence only.
* **I4.** The `|P0|` range is misreported (§6.3).
* **I5.** Table §2.3: the "C" counts exclude collisions with the genuine unit (`classify.py` counts them as GEN),
  which is not stated. My counts include them; A0/A1 counts agree.
* **I6.** Strictness wording: "local twist" should read "local `(eps, t)` symmetry" (§6.6).
* **I7.** The `K = 5` claims rest on the round-14 `F_5^ex` files; see §6.1 for my independent recomputation.
* **I8.** The interval B&B certificates (`slab_iv.py`) are superseded by Theorem B(3) plus exact index-1 data. I did
  not re-run them; nothing depends on them any more.

## 9. Corrected statement

**No growth improvement and no uniform bound.** `M <= (2/3 + eps) log R/loglog R` stands. P3-all\* is uncapped in
regime (1.4), and Rig^SL(K) is undecided (no joint SL core is determined).

1. **Proved (refereed here).**
   * Lemma 2.1 (J0 = OPP ∧ ADJ ∧ MAX) and the four-unit laws.
   * Theorem A (C/A0/A1). A0 collisions fit at no `phi`.
   * Lemmas B0, B1; Theorem B (one-modulus pinning: `phi in phi_0 + pi Lambda_T^*`; if `Lambda_T = Z^r`, one-modulus
     SL = the integer cofactor criterion); Cor. B2.
   * Theorem C's inequality (2.1), with the corrected reading: an extra beyond the cofactor slack, on a monomial whose
     opposite unit is saturated, forces a dropped genuine opposite-unit collision at some `q_1 in F` dividing
     `D_(s+2)(Z_v)`. This is a modulus other than `q_0` only if `q_0` does not divide `D_(s+2)(Z_v)`. Counterexample to
     "second modulus": `first4`, `v = (-3,1,0,0)`, `q_0 = 19`.
   * Theorem P: a P3-all\* member that is `(eps, t)`-strict (`u = 1 mod 4`, all levels `<= X`) on a set `S` with
     TC(`S`) is an actual cluster. So P3-G has no fake on prime sets with TC.
2. **Scale.**
   * Every actual cluster has `r > (M/4)(1 + 2M log^+ C/W')^(-1) = (M/4)(1 - O_C(1/log M))` (`stability.md`
     Thm 3.3(3); round 14's "`r >= M^gamma`?" is answered with `gamma = 1 - o(1)`).
   * The same holds for members of `A_X ⊇ P3-char ⊇ P3-all*` with `X >= 2M`. That proof is a check of the inputs of
     `stability.md` Prop. 7.1: plausible, not re-derived to `two_thirds`.
   * The forms `r > M/4` and `r >= (M/2)(1 - M^(-1/3))^2` (`C <= 1`) are empty for actual clusters.
   * Proved: `Q_eta ⊆ F_2^ex` for `M >= M_0(C)`. Not proved for `C > 1`: "every `q <= 2M` is forced at norm 2"; only
     `q < 2M (1 + 2M log^+ C/W')^(-1) + 2` is.
3. **Exact, independently reproduced.**
   * The 83/163 A0 explanation.
   * `[Z^4 : Lambda_T] = 1` at all forced `q_0 >= 7` (four sets, `K = 4, 5`).
   * One-modulus SL segments `first4` `[3,103]` (`K = 4`) and `[3,163]` (`K = 5`), `rand4` `[3,47]`, `bal4` `[3,41]`,
     `last4` `[3,29]` (`K = 4`), with lifted witnesses verified at `phi_0`. Hence
     `Q_core^SL(rand4,4) <= 47`, `Q_core^SL(bal4,4) <= 41`, `Q_core^SL(last4,4) <= 29`,
     `109 <= Q_core^SL(first4,5) <= 163`.
   * New (referee): `Q_core^SL(rand4,5) <= 59`, `Q_core^SL(bal4,5) <= 61`, `Q_core^SL(last4,5) <= 41`.
   * Strictness data (79/89/103 at `K = 4`, 157/163 at `K = 5`).
   * TC certificates for the first `r = 4..64` split primes (all moduli) and at `Y = 4r` (`r = 16, 32, 64`), and
     `Y = 2r` (`r = 64`).
4. **Evidence only.**
   * Genericity of TC and of covering.
   * Growth of the one-modulus SL segments with `K` (loose sets: `rand4` first flexible `59 -> 67`, `bal4`
     `43 -> 67`, `last4` `31 -> 43`). One-modulus segments give no lower bound on the joint SL core.
   * Strictness growing with `K`.
