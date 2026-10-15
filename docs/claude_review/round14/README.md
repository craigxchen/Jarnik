# Round 14 -- judge's report: Rig(K), C1 in regime (1.4), pal11

## 0. Bottom line

**No growth improvement and no uniform bound.** `M <= (2/3 + eps) log R/loglog R` stands. P3-all\* is still open
and uncapped in regime (1.4) (`p_j <= M^(h_0)`, `X = 3M^A`, `A` chosen after `c_0`), with
`3 M log M - O(M) <= W_min^all <= O(M log^2 M/loglog M)`. No fake of any kind is known there.

**Main finding of the round: Rig(K), as round 13 posed it, used the wrong class.** "Good = per-unit budget" (PU)
drops the angle coupling of (1.2): one `phi` must serve all four units of every monomial (found by the rigcompute
referee; the judge re-proves its exact per-monomial form as Lemma J0, §5, and settles the referees'
disagreement). Rigidity certificates in PU transfer to P3-all\*; flexibility witnesses do not. Every joint PU
witness of the round fails the coupling: 10 rigcompute witnesses (referee) and 24 rigtheory witnesses (judge,
exact).

**Rig(K): undecided.**
* *Growth in `K`* is supported for the smallest primes `{5,13,17,29}`. Certified lower bounds on the coherent
  joint core are `61, 109, 193, 239` for `K = 4..7`, and these transfer to P3-all\*. In PU the core is pinned
  down exactly: 37 at `K = 3` (per-modulus conjugation; not coherent there), then 61, 109, 193 for `K = 4, 5, 6`
  (coherent; the `K >= 5` certificates are single-implementation). In SL `61 <= Q <= 103` at `K = 4`, so growth
  `K = 4 -> 5` is strict in SL too.
* *Growth in `r` saturates* at `K = 4`: core `<= 173` for `r = 6, 7, 8` (survives in SL).
* *A proved ceiling* `q_core(K) = O_K(r^floor(K/2))` means bounded `K` reaches the scale `M` only when
  `r >= M^gamma`.
* *The loose-prime "stall"* is proved only in PU.
* *The transfer theorem* (rigidity + H^gen + L1 ⇒ `W >~ M log^2 M/loglog M`) is valid only after a correction
  for associate twists.

**C1 in regime (1.4): open.** Proved obstructions kill the random fibre (norm 1), the torsion-avoiding fibre
(norm 2), the L4K systems (norm 4) and exactly coherent designs. No exhibited flexible top admits angles. The
problem is localised to DeepBal + PairCompat. Norm-1 pinning holds only on prime sets satisfying a condition that
small-prime sets usually violate.

**pal11: a proved negative for route-5 certificates** (no product certificate of any kind; exact, re-proved
independently). The tangent-cone question is undecided, and the result holds equally for fakes and actual clusters.
It also corrects a round-13 record: `zeta(M8') = 10`.

**Scores** (rubric of `round10/README.md` §1; details in §7):

| write-up | (a) | (b) | (c) | **overall** |
| --- | --- | --- | --- | --- |
| rigtheory | 6 | 2 | 1 | **4** |
| rigcompute | 6 | 2 | 1 | **4** |
| c1construct | 5 | 3 | 1 | **4** |
| pal11cert | 8 | 0 | 0 | **3** |

**Recommendation (§8).** Decide **Rig^SL**, the norm-`<= K` truncation of P3-all\* with labels and one angle vector,
starting with the open branch: do loose primes still stall once angles are imposed? Certify rigidity with the exact
box test of Lemma J0, and require a verified `phi` for every witness. Alongside it, decide whether near-extremal
clusters in regime (1.4) must have `r >= M^gamma`. Without that, bounded-`K` rigidity cannot reach the scale `M`.

## 1. Referee-corrected statements

Four write-ups, four adversarial referees with independent code. Labels: *proved*; *exact* (finite computation in
exact integer/rational arithmetic, infeasibility certified by a complete search); *solver* (floating-point MILP or
interval branch-and-bound, no certificate); *evidence*. Throughout, the prime set `first4 = {5,13,17,29}` is used
unless said otherwise; `F_K^ex` is the exact norm-`K` forced set (`#F = 7, 17, 40, 65` for `K = 2..5`).

**Two label classes must be kept apart in everything below** (this is the round's main correction, from
referee_rigcompute I1, confirmed by the judge in §5):
* **PU** ("per-unit"): parity-respecting labels at all levels with `m_(v,s)^2 <= Norm Z_v * 2^[s odd]` for every
  `1 <= ||v||_1 <= K` and every unit `s`. This is Rig(K) exactly as round 13 §7.2 worded it, and the class of both
  rigcompute and rigtheory.
* **SL** ("slab"): PU plus the slab system (1.2) of `sharp.md` Def. 1.1 for those `v`, with **one** angle vector
  `phi in R^r`. SL is necessary for a label system to be the norm-`<= K` part of a P3-all\* member. PU is a strict
  relaxation of SL. A cheap necessary condition for SL is per-`v` angle compatibility **ANG** (§5, Lemma J0).

Consequently: **lower bounds on rigid cores (certificates) proved in PU transfer to SL and to P3-all\*; flexibility
witnesses found in PU do not**. The rigtheory referee's statement that Rig(K)'s class is "exactly the P3-all\* budget
data restricted to norm `<= K`" is right about budgets, parity, units and moduli, but overlooks the angle coupling;
on this point the rigcompute referee is correct (judge, §5).

### 1.1 rigcompute (first finite test of Rig(K)): computations survive, class-dependent readings do not

Write-up `rigcompute.md`, report `referee_rigcompute.md`. Referee verdict: computations correct, result as headlined
does not survive.

*Corrected statement (evidence about Rig(K); finite results exact in the class named).*
1. **Lower bounds (exact in PU, hence valid in SL and P3-all\*).** The coherent joint core (every good system is
   coherent with one column-sign conjugate `eps∘L` of the genuine labels on `F_K^ex ∩ [3,Q]`, all moduli free) is
   * `Q_core(first4, K) >= 61, 109, 193, 239` for `K = 4, 5, 6, 7`;
   * `Q_core(first5, K) >= 101, 163` for `K = 4, 5`.

   These come from complete C searches over the exact relaxation of Lemma R. The `first4, K = 4` bound is
   independently confirmed by the referee's MILP (infeasible at `Q = 61`, feasible at `Q = 67`; solver). The
   `K >= 5` certificates are single-implementation. Lemma M (monotonicity, side condition verified) carries them
   to the all-norm core.
2. **Upper bounds, exact in PU only.** `Q_core^PU = 61, 109, 193` (`first4`, `K = 4, 5, 6`), `101` (`first5`,
   `K = 4`), `7` (`rand4`, `bal4`, `last4`, `K = 4, 5`), `<= 173` (`first6..8`, `K = 4`). All 10 joint witnesses
   are verified at all levels, but **none lies in SL**: 8 fail ANG (exact; e.g. `rand4`, `K = 4`, `v = (1,0,0,1)`:
   `65^2 + 87^2 = 11794 > N = 9169`), the two `last4` witnesses pass ANG but have no `phi` for their slab system
   (complete interval branch-and-bound with a `10^-9` margin, solver; minimal infeasible core of 10 monomials).
3. **In SL** (referee, verified witnesses at `phi = arg pi_j` plus transferred lower bounds):
   * `61 <= Q_core^SL(first4, 4) <= 103` (SL witnesses at 107); hence `Q_core^SL(first4, 5) >= 109 >
     Q_core^SL(first4, 4)`: **strict growth from `K = 4` to `K = 5` holds in SL**;
   * `Q_core^SL <= 173` for `first6..8`, `K = 4` (the first-FLEX witnesses at 179 satisfy the slab system at
     `phi = arg pi_j`);
   * `7 <= Q_core^SL(last4, 4) <= 29`, `7 <= Q_core^SL(last4, 5) <= 41`;
   * one-modulus: at `first4`, `q = 83` (`K = 4`) and `q = 163` (`K = 5`) are rigid under ANG (exact), so the PU
     values `q_flex1 = 83, 163` are artefacts of the relaxation; `rand4` at 31 and `bal4` at 41 are rigid in SL at
     `K = 4` (solver).
4. **Not established:** the loose-prime stall in SL or P3-all\*; any C1 seed; the scale law
   "core `~ 0.4-0.7 |S_floor(K/2)|`" beyond PU at `r = 4, 5`; anything in regime (1.4).
5. **Minor corrections.** "local `<=` joint `<=` one-modulus" is an observed comparison, not an implication;
   Lemma K proves `K`-independence only for `q <= p_min^((K+1)/2)`, and "local rigidity ends at 31" is checked for
   `q <= 113`; local rigidity is not vacuous for `rand4` (`q = 17`); the `first8` witness at 179 needed a lift
   (supplied and verified by the referee); "tight" sets have `log p/log r` up to 2.43, not near 1; depth
   certificates need `X >= 27, 49, 121, 169`.

### 1.2 rigtheory (theory of Rig(K)): survives as "undecided", with a corrected transfer theorem

Write-up `rigtheory.md`, report `referee_rigtheory.md`. Referee verdict: survives after corrections.

*Corrected statement.*
1. **Proved.** Theorem 2.5 (only `F_K^ex` matters; cores monotone in `K`; `q_core(K) <= 4A_K + 5 =
   O_K(r^floor(K/2))`); Lemma 3.1 (genuine slack identity: the slack of `(v,s)` is
   `log U(D_s(Z_v)) + log(1/|sin(arg Z_v - s pi/4)|)`, `U` the uncharged cofactor) with Cor. 3.2 (read "not a genuine
   collision with unit `s`"); Lemma 3.3; **Theorem 3.4 (anticode rigidity)**: if the genuine label maps an
   `l1`-anticode of diameter `<= K-1` onto `Z/n_q`, every label outside the genuine class has a non-genuine collision
   of norm `<= K`; Cor. 3.5 (no F2 deformation at covered moduli; the joint transport inequality);
   Prop. 3.6(1-2); Prop. 3.7.
2. **Theorem 4.1, corrected.** The class has the symmetry group `{+-1}^r x (Z/2)^r` (conjugations and associate
   twists `pi_j -> i pi_j`; referee Prop. T, proved). "Strict and coherent" must be read up to `(eps, t)`:
   `ell(q^a) = u_q (eps∘L(q^a) + (m_(q^a)/2) t)`, `u_q = 1 mod 4`, one `(eps, t)` on `Q_eta`. Then
   RigS + H^gen (for the primes `i^(t_j) pi_j^(eps_j)`) + L1 gives
   `W >= ((8(1-eta) - eps)/K_u) M log^2 M/loglog M` for P3-all\* members in regime (1.4). **As written, RigS implies
   the class is empty in the regime** (Lemma T': unconditionally if some `p_j <= (log M)^C`, under GRH always), and
   the twisted genuine system is non-strict at 14/17, 37/40, 61/65 forced moduli (`K = 3, 4, 5`; exact). The size
   claim "`q_core >= M` at bounded `K`" needs the extra hypothesis `r >= M^gamma`, which regime (1.4) does not give.
   For actual circles alone the theorem is vacuous: it only routes P3-all\* to statements about genuine labels.
3. **Exact, independently reproduced.** Forced sets `K = 2..5`; one-modulus relaxed cores
   `q_core^1 = 7, 41, 79, 157` (`K = 2..5`) with every round-10 count compared (167: 120/96, 199: 120/40, all F1);
   `K = 6`: 263 rigid, 269 flexible; all 22 joint PU witnesses (all levels).
   `Core^eps_J(3) = F_3^ex ∩ [3,37]` exactly (37 already from the subset `{3, 37}`, i.e. from the 3-adic levels up to
   81: a small-modulus effect); conditional on genuine labels at 3 and 7, `Core^eps_J(4) ⊇ F ∩ ([11,61] ∪ {73})`
   exactly. Solver-only: `K = 4` coherence given 3, 7; the `K = 5` claim given 3, 7 and every modulus above 163.
4. **Evidence only.** The "anticode threshold explains round 10" law (Theorem 3.4 excludes only F2 at covered moduli;
   137-157 are rigid but uncovered at `K = 5`, the onset at 163 is F1-only; the law fails at `K = 7` and `r = 5`);
   growth of the one-modulus core at `r <= 6`; "the all-norm reference is more rigid".
5. **Judge's addendum (§5, exact, unrefereed).** All 24 stored joint witnesses of rigtheory are PU-good but fail ANG
   (14-28 monomials each), and none of their 16 associate-twist variants is PU- and ANG-good. So every rigtheory
   upper bound on a joint core is a PU statement only. Combined with rigcompute's certificate,
   `Core^coh_J(first4, 4) ⊇ F ∩ [3, 61]` **unconditionally** (no genuine labels fixed at 3 and 7), which removes the
   referee's "only `37 <= q_core^(J,eps)(4)`": in PU the joint core grows `37 -> 61 -> 109 -> 193` for
   `K = 3, 4, 5, 6`.

### 1.3 c1construct (C1 in regime (1.4)): restricted theorems survive; pinning consequences do not

Write-up `c1construct.md`, report `referee_c1construct.md`. Referee verdict: survives as a restricted theorem; no C1
fake, no lower bound, growth unchanged.

*Corrected statement (regime (1.4): `p_j <= M^(h_0)`, `X = 3M^A`, `A` after `c_0`).*
1. **Break theorems (proved).**
   * B1: for the random fibre of `sharp.md` Thm 6.4, each non-twin column meets its norm-1 per-unit budget with
     probability `<= (1+o(1)) log Y max(log Y, log M)/(log X)^2`, `Y = sqrt(2 p_j)` (the product skips `q | N`);
     `-> 0` as `A` grows, for `A > h_0/2`. A statement about a random law, not about every member of its support.
   * B2: the torsion-avoiding fibre fails at norm 2 for `A > h_0` (same caveat).
   * B3: `K in {2,3}`, `h = 1`, `beta_j = j`: in the L4K Thm 2 systems at least `r^3/36 - O(r^2)` norm-4 monomials
     vanish identically at every `q > q_T` for every choice of permutations (norm 4 = `K+1` for `K = 3`, `K+2` for
     `K = 2`).
2. **Norm-1 pinning (proved, narrow).** Lemma P0, Prop. P, Cor. P1-P2: a column satisfying condition (i) (no `p_l`
   divides the pinning numbers `D_s(pi_j)`) whose labels are genuine per unit at its pinning prime powers is pinned
   at norm 1; a P3-G member satisfying (i) in every column is an actual cluster. **But (i) fails for 80-91% of the
   columns of the first `r` split primes, and for all but one column of selected sets** (exact); even with genuine
   labels everywhere, 50-77% of such columns are unpinned at norm 1 (witness: `p = 13`, `5 | D_3(3+2i)`, free arc
   of `phi_13` of length 0.307). So "P3-G collapses at norm 1", "flexible tops cannot seed a fake", "Rig_pin makes
   P3-all\* exact" and "the label-angle step is settled" are **not** established; they hold only on prime sets
   satisfying (i), which the fake controls.
3. **Coherent designs (proved, Linnik/Siegel quoted).** Lemma L (matching the parity on a deep class forces
   `p_j Norm zeta` to be a square or twice a square); Prop. G and Thm C (exactly coherent members with
   `Norm xi_j <= M^h` and the cofactor conditions are actual clusters).
4. **Pairs (proved).** Thm R: on the Paley profile, G-coherent designs with `X >= 3M^4` have `C + 2Hbar > (M-1)/8`
   (round 7: `M^(1/2-o(1))`).
5. **Shallow classes (proved).** Thm S, Cor. S1: the classes `k(q) <= E_0` are served `X`-uniformly by
   root-repaired coherence, for monomials outside the character lattice, at cost `(2^E_0 - 2) sum |v_j| log l_j`.
6. **Open for C1:** DeepBal and PairCompat. Root-repaired coherence on the deep classes is uncertified, **not
   excluded** (one-instance evidence against it).

### 1.4 pal11cert (route-5 barrier data point): survives; out of this round's scope

Write-up `pal11cert.md`, report `referee_pal11cert.md`. Referee verdict: survives; major issue is scope only.

*Corrected statement.*
1. **Theorem 4.3 (exact, computer-assisted, independently re-proved by the referee by a simpler route).** On the
   pal11 box (`n = 12`, `T = 1`, `L = 2`, `|A| = 12093`), with `V` the 1519 test directions of width `tau_A(v) < 8`,
   there are no nonzero real `z'_1..z'_8 in R^11` with `#{i : <z'_i, v> != 0} <= tau_A(v)` for all `v in V`. The
   width counts are necessary for real generators (fatpoint Prop. 2.6(c), whose proof holds verbatim for real
   `z_i`), so `J(A)_8^perp` contains no product of 8 real linear forms: pal11 has no product certificate at all.
2. **Not decided:** whether `TC(pal11, 8)` has a nondegenerate real point (numerical evidence only: 200 minimisations
   end on certified coordinate subspaces or within `3.4e-4` of collision strata). Prop. 6.2: `TC(pal11,8)` contains
   the 12 five-dimensional twin-pattern coordinate spaces (exact).
3. **Symmetry.** `Aut(A) = Z_5 x {+-1}` (exhaustive); the "Z_11" of the task belongs to Paley-12.
4. **Round-13 correction (confirmed).** fatpoint's `prodsearch.py` capped base incidence at 2, while
   `tau_A(1) = 6, 8, 6, 4` (M8, M8', pal7 `L = 2`, pal11). Correct exact tables: M8' (syl8 coordinates, `L = 3`) at
   `m = 10` has 1092 count-feasible products, all zonotopes, so `zeta(M8') = 10` (round 13 said none); M8 at
   `m = 8`: 42 count-feasible, 0 certificates (rigorous); pal7 (`L = 2`) at `m = 7`: 2134 count-feasible, exactly the
   300 zonotopes are certificates.
5. **Scope.** The statement concerns one fixed box and holds identically for fakes and actual clusters with this
   pattern set: it caps nothing, uncaps nothing, and has no bearing on Rig(K) or P3-all\*.

## 2. Status of Rig(K)

**Verdict: undecided. Neither the growing nor the stalling alternative is proved; the finite evidence is mixed and
depends on the prime set; and the question as worded in round 13 used the wrong class.**

1. **The question must be restated (proved, §5 Lemma J0).** Rig(K) with the per-unit budget (PU) drops the angle
   coupling of (1.2). Certificates (rigid cores) in PU are valid for P3-all\*; flexibility in PU says nothing about
   P3-all\*. The correct finite question is **Rig^SL(K)**: the norm-`<= K` truncation of P3-all\* itself, i.e.
   parity-respecting labels **and** one `phi in R^r` satisfying (1.2) for every `1 <= ||v||_1 <= K`. Every witness
   must carry a verified `phi`.
2. **Evidence for growth in `K` at fixed small `r` (both classes).** `first4`:

   | `K` | core, PU (exact both sides unless marked) | core, SL |
   |---|---|---|
   | 3 | `Core^eps_J = [3, 37]` (rigtheory + referee, exact) | `>= 37`; no SL witness known (all 12 PU witnesses fail ANG, §5); `F_3^ex` ends at 103 |
   | 4 | `Q_core^coh = 61` (rigcompute; MILP-confirmed) | `61 <= Q <= 103` (referee; judge re-verified, §5) |
   | 5 | `109` (single-implementation certificate) | `>= 109` |
   | 6 | `193` (single-implementation certificate) | `>= 193` |
   | 7 | `>= 239` | `>= 239` |

   In SL the only proved strict increase is `K = 4 -> 5` (`<= 103` vs `>= 109`). In PU the increase is strict at
   every step `K = 3..6`.
3. **Evidence against polynomial growth in `M`.**
   * *Saturation in `r` (survives in SL).* At `K = 4`, `Q_core <= 173` for `first6, first7, first8`: adding primes
     41, 53, 61 does not move the first flexible modulus 179. Growth in `r` stops once added primes are no longer
     the smallest.
   * *A proved ceiling.* `q_core(K) <= 4A_K + 5 = O_K(r^floor(K/2))` (rigtheory Thm 2.5). So at bounded `K` the
     core can reach the scale `M` that route 3 needs only if `r >= M^gamma` (`K >= 2/gamma + 1`). Regime (1.4)
     only gives `r <= rho M`; near-extremality forces `sum e_j >~ 3M/h_0`, not many columns. For clusters with
     `r = M^(o(1))` a Rig(K) theorem at bounded `K` cannot give RigS at all; `K` must grow with `M`, where none of
     this round's theory applies.
   * *Loose primes (PU only).* For `rand4`, `bal4`, `last4` (all `p_j >= 53`) the PU core is 7 (vacuous) at
     `K = 4, 5`. In SL this is unproved: the stored witnesses all fail SL, and SL rigidity appears at 31 (`rand4`)
     and 41 (`bal4`) one modulus at a time (solver). In SL only `7 <= Q_core^SL(last4, 4) <= 29`,
     `<= 41` at `K = 5`, are known.
4. **Theory.** Proved: reductions to `F_K^ex`, monotonicity in `K`, the genuine slack identity, the anticode
   rigidity theorem (no F2 deformation on covered moduli), the transport inequality, and the corrected transfer
   Theorem 4.1 (RigS up to `(eps, t)` + H^gen + L1 ⇒ `W >~ M log^2 M/loglog M` for P3-all\* members). Missing:
   **F1 exclusion** (a counting lemma: no index-`n` lattice other than the genuine one has all its short vectors in
   the affordable set), covering (an L7/H-type equidistribution statement that fakes share), and the joint,
   coherent, strict upgrades on the window `(M^(1-eta), M]`.
5. **What either answer would buy (unchanged logic, corrected).** A growing Rig^SL with the `(eps, t)`-strict
   conclusion on `Q_eta` routes P3-all\* to genuine-label statements (H^gen, L1, i.e. L7-type inputs) in regime
   (1.4). It does **not** by itself give growth, and for actual circles it is vacuous. A stalling Rig^SL would not
   cap P3-all\* unless the flexible tops complete to all norms, all moduli up to `X`, and the pair conditions; no
   exhibited top completes even at norm `<= K` jointly (§3).

## 3. Status of the C1 construction in regime (1.4)

**Verdict: open. No C1 fake was built, none was obstructed, and no fake of any kind is known in regime (1.4).
The round adds proved obstructions to every concrete design tried, and no surviving seed.**

* **Refuted as seeds (proved / exact).**
  * The random fibre of `sharp.md` Thm 6.4 fails at norm 1 with probability `-> 1` as `A` grows (B1), the
    torsion-avoiding fibre at norm 2 (B2), and the deterministic L4K Thm 2 systems at norm 4 for every choice of
    permutations (B3). B1-B2 are statements about random laws.
  * Exactly coherent designs with small height are actual clusters (Thm C, given cofactor conditions); G-coherent
    designs on the Paley profile pay `C + 2Hbar > (M-1)/8` (Thm R).
  * **Every joint flexible top exhibited this round fails SL**: the 10 rigcompute witnesses (referee: 8 fail ANG
    exactly, 2 have an infeasible slab system, solver) and the 24 rigtheory witnesses (judge, §5: all fail ANG
    exactly, also after any associate twist). So the round-13 plan "take the flexible tops as seeds and add the
    slab conditions with angles chosen after the labels" fails at its first step for every top found.
* **What survives as possible raw material (exact, one modulus, norm `<= K` only).** SL-valid non-coherent
  one-modulus deviations exist: `first4`, `K = 4`, at `q = 107` (108 labels outside every `eps`-class at level 1;
  three completed to all levels and checked at `phi = arg pi_j`, referee and judge independently); `first6..8`,
  `K = 4`, at 179; `last4` at 31 (`K = 4`) and 43 (`K = 5`). They are single-modulus deviations at bounded norm,
  paid for by large uncharged cofactors `U(D_s(Z_v))` at `phi = arg pi_j` (Lemma 3.1). A C1 member needs the
  opposite: deviations at a positive proportion of moduli up to `X`, all norms, and pairs.
* **Localisation (proved).** Thm S serves the shallow classes `k(q) <= E_0` `X`-uniformly by root-repaired
  coherence for monomials outside the character lattice. What remains for C1 is **DeepBal** (balanced,
  non-coherent labels on the deep classes) and **PairCompat** (pair demands `O(log M)` uniformly in `A`). No
  mechanism is known for either; root-repaired coherence on the deep classes is uncertified, not excluded.
* **Pinning (narrowed).** Norm-1 pinning (Prop. P) forces actual clusters only on prime sets satisfying (i). Sets
  with small split primes violate (i) in 80-91% of columns, so the fake keeps the freedom to choose primes where
  norm-1 pinning fails. Higher pinning closures (norm `>= 2`, `p3G.md`) are open; on `{5,13,17,29}` norm-2
  relations pin columns 13 and 17 numerically (evidence).

## 4. Status of the pal11 barrier data point

**Verdict: proved negative for the certificate method; the tangent-cone question that matters is undecided.**

* **Proved (exact, computer-assisted, two independent routes).** The fake pattern box pal11 (`T = 1`, `L = 2`,
  `m = 8`) admits no product certificate of any kind: no 8 real linear forms satisfy the width counts. Product
  certificates therefore cannot decide route 5 on this first overdetermined, faithful, fake box.
* **Undecided.** Whether `TC(pal11, 8)` has a nondegenerate real point, which is what MIN^pair/ZE-TC needs on this
  fake (round 12-13 barrier test). Numerical evidence (200 minimisations) leans towards "no", but it is evidence
  only, and Theorem 4.3 holds identically for fakes and actual clusters with this pattern set.
* **Record correction (confirmed).** Round 13's product searches were not exhaustive (`prodsearch.py` capped base
  incidence at 2): M8' (syl8, `L = 3`) has 1092 zonotope certificates at `m = 10`, `zeta(M8') = 10`.
* **Barrier reading.** Neither for nor against the fake barrier for route 5: the result neither shows a fake
  violates a route-5 input nor that one satisfies it. It supports keeping route 5 paused (round 13 §7.1): the
  cheapest remaining certificate test is closed with a negative.

## 5. Judge's own checks (unrefereed)

All scripts are in `judge_checks/` (own code; the helper `jglib.py` is a copy of the rigtheory referee's
`glib.py`, residue-field discrete logs; nothing is imported from the repository or from an author's directory;
author pickles were opened only after checking that they contain data opcodes only). **None of this has been
refereed.**

**Lemma J0 (per-`v` angle compatibility; proved).** Let `N = Norm Z_v`, `x = <v, phi>`, `y = 2x`, and
`m_s = m_(v,s)`. Condition (1.2) of `sharp.md` for `s = 0..3` (`nu_0 = nu_2 = 1`, `nu_1 = nu_3 = 2^(-1/2)`) is
equivalent to

```text
2 m_2^2/N - 1 <= cos y <= 1 - 2 m_0^2/N,        m_3^2/N - 1 <= sin y <= 1 - m_1^2/N.
```

*Proof.* `sin^2 x = (1 - cos y)/2`, `cos^2 x = (1 + cos y)/2`, `sin^2(x - pi/4) = (1 - sin y)/2`,
`sin^2(x - 3pi/4) = (1 + sin y)/2`; square the four conditions. QED. Hence (1.2) for one `v` holds for some `x`
iff the rational box `B_v` above meets the unit circle, iff `B_v` is non-empty and
`min_(B_v) |z|^2 <= 1 <= max_(B_v) |z|^2` (`B_v` is convex). In particular (1.2) implies
`m_0^2 + m_2^2 <= N` and `m_1^2 + m_3^2 <= 2N`, while the per-unit budget allows `2N` and `4N`. Genuine labels
satisfy `B_v` at `y = 2 arg Z_v`. This settles the disagreement between the two Rig(K) referees in favour of the
rigcompute referee: PU is a strict relaxation of the norm-`<= K` part of P3-all\*, and the condition that is lost
is exactly this per-`v` coupling plus the linear coupling `x_v = <v, phi>` across `v`. At `phi = arg pi_j`,
(1.2) for `(v, s)` reads exactly `m_(v,s) <= |D_s(Z_v)|` (an integer test).

**J1 (exact): every rigtheory joint witness fails ANG** (`jang.py`, output `jang_out.txt`).
* Input: the 24 stored joint witnesses of `rigtheory_checks/` (`first4`; `K = 3`: genuine at 3, 7, deviation at
  `q_0 = 19, 41, 43, 47, 53, 59, 61, 67, 71, 83, 103`; `K = 4`: genuine at 3, 7, 11, deviation at `q_0 = 67, 71, 79,
  83, 89, 97, 101, 103, 107, 127, 151, 199`). Demand: genuine moduli at all levels, others at level 1 (level-2
  collision-freeness re-checked).
* Result: all 24 are PU-good (0 violations, agreeing with the referee's `vwit.py`), and **all 24 fail ANG**, on
  14-28 monomials each (e.g. `K = 3`: `v = (-2,0,0,-1)`, `N = 725`, `m = (7, 33, 23, 37)`,
  `33^2 + 37^2 = 2458 > 2N = 1450`). Applying any of the 16 associate twists `ell -> ell + (m/2) t` to the free
  labels never gives a system that is both PU-good and ANG-good (only `t = 0` is PU-good).
* Sanity (`jang_sanity.py`, `jang_sanity_out.txt`): with genuine labels substituted at the same moduli, label
  collisions agree with `q | D_s(Z_v)` in 660,480 checks (0 mismatches) and every `B_v` test passes. So the unit
  conventions of labels, `D_s` and (1.2) match.
* Consequence: in SL there is no known upper bound on the joint core at `K = 3`, and the rigtheory bounds
  `q_core^(J,eps)(3) <= 37`, `(4) <= 61` are PU statements only.

**J2 (exact): independent confirmation of the SL one-modulus data** (`jonemod.py`, `jlift.py`; outputs
`jonemod_out.txt`, `jlift_out.txt`). Own exhaustive enumeration of all parity-respecting level-1 labels at `q_0`,
genuine labels at all levels on `F_K^ex \ {q_0}`, `first4`:

| `K` | `q_0` | PU: all / non-genuine / outside every `eps`-class | ANG | SL at `phi = arg pi_j` |
|---|---|---|---|---|
| 4 | 83 | 36 / 24 / 24 | 12 / 0 / 0 | 12 / 0 / 0 |
| 4 | 107 | 216 / 144 / 144 | 180 / 108 / 108 | 180 / 108 / 108 |
| 5 | 163 | 160 / 80 / 80 | 80 / 0 / 0 | 80 / 0 / 0 |

The PU and ANG columns agree exactly with referee_rigcompute §4.6 (36/24 -> 12/0 at 83; 216/144 -> 180/108 at 107;
160/80 -> 80/0 at 163). So 83 (`K = 4`) and 163 (`K = 5`) are one-modulus rigid under ANG, confirmed independently
(`K = 5` run: `jonemod_fast.py`, which merges `v` and `-v`; it reproduces the 83 row of the unmerged script).
New: at 107 **every** ANG survivor satisfies the whole norm-`<= 4` slab system at the genuine angles. Three of them,
`ell = (1,4,61,92), (1,88,43,50), (1,88,97,104)` (judge's conventions: `pi = 1+2i, 2+3i, 1+4i, 2+5i`), were lifted
to level 2 without level-2 collisions (lifts in `jlift_out.txt`) and re-checked with the all-level demand: 0 PU
violations, 0 ANG failures, 0 slab violations at `phi = arg pi_j`, label outside every `eps`-class. This
independently confirms `Q_core^SL(first4, 4) <= 103`.

**J3 (synthesis of refereed results; no new computation).** rigcompute's coherent core and rigtheory's
`Core^eps_J` use the same class notion mod `n` (`c' = u eps∘L`, `u` a unit of `Z/n`), and coherence on a set
implies membership of each modulus in `Core^eps_J`. Hence rigcompute's certificate `Q_core(first4, 4) >= 61`
(all moduli free, MILP-confirmed) gives `Core^eps_J(first4, 4) ⊇ F ∩ [3, 61]` without the condition "genuine at 3
and 7" that the rigtheory referee had to impose. With rigtheory's verified witness at 67, `q_core^(J,eps)(4) = 61`
in PU, and the PU joint core grows strictly `37 -> 61` from `K = 3` to `K = 4`. (PU only; in SL see §2.)

## 6. Growth rate and uniform bound

**No change.** No growth improvement over `log R/loglog R` and no uniform bound is proved, unconditionally or from
a standard conjecture. `M <= (2/3 + eps) log R/loglog R` (`two_thirds`) stands, and so does the P3-all\* bracket
`3 M log M - O(M) <= W_min^all <= O(M log^2 M/loglog M)`. P3-all\* remains open and uncapped in regime (1.4); no
fake of any kind is known there.

What the round changes in the map (barrier-map corrections only, which do not count as growth progress):
* **The route-3 chain is re-specified.** Its first link must be Rig^SL (labels plus one angle vector), not
  per-unit Rig(K). Its conclusion must be strictness up to `(eps, t)` (conjugation and associate twist). Its size
  claim holds at bounded `K` only when `r >= M^gamma`. With these corrections, Theorem 4.1 (rigtheory) is a valid
  conditional statement: Rig^SL-type rigidity on `(M^(1-eta), M]` + H^gen + L1 ⇒
  `W >= ((8(1-eta) - eps)/K_u) M log^2 M/loglog M` for P3-all\* members, i.e. growth. For actual circles alone the
  rigidity hypothesis is vacuous and the inputs that matter are H^gen and L1 (L7-type, round 9), which are
  unchanged and open.
* **The cheapest C1 seeds are gone.** Every PU flexible top fails SL; the random and L4K designs fail at norms 1, 2
  and 4; exactly coherent small-height designs are actual clusters. C1 is localised to DeepBal + PairCompat.
* **Route 5** stays paused; pal11 closes its cheapest certificate test with a negative.

## 7. Scores

Rubric of `round10/README.md` §1: (a) correctness after refereeing; (b) fake-violated input; (c) movement of growth
or uniformity. The overall score is a judgement weighing (b) and (c) above (a). No write-up moves the growth rate.

| write-up | (a) | (b) | (c) | **overall** | reason |
| --- | --- | --- | --- | --- | --- |
| **rigtheory** | 6 | 2 | 1 | **4** | Real theorems: reductions to `F_K^ex` with the ceiling `q_core = O_K(r^floor(K/2))`, the genuine slack identity, anticode rigidity (no F2 deformation), the transport inequality, and a transfer theorem that becomes valid after the `(eps, t)` correction. Exact one-modulus data reproduced throughout. Deductions: RigS as written contains its conclusion (missed associate-twist symmetry); the size claim silently needs `r >= M^gamma`; joint rigidity claims were solver outputs labelled exact and were conditional; the "anticode threshold explains round 10" headline is a coincidence. |
| **rigcompute** | 6 | 2 | 1 | **4** | The round's most solid data: exact forced sets, complete-search certificates `Q_core(first4, K) >= 61, 109, 193, 239` (transfer to P3-all\*), and the `r`-saturation at 179 (survives in SL). Deductions: the class drops the angle coupling, so "exact on both sides", the loose stall and the C1-seed reading fail; all 10 witnesses lie outside SL; `q_flex1 = 83, 163` are relaxation artefacts. |
| **c1construct** | 5 | 3 | 1 | **4** | Correct break theorems B1-B3, Lemma L (unconditional Linnik version of the GRH parity obstruction), Thm C, Thm R (linear-in-`M` pair cost for coherent designs on Paley), and Thm S (shallow classes served `X`-uniformly), which localise C1 to DeepBal + PairCompat. Deductions: condition (i) of Prop. P was called rare on vacuous evidence (it fails for 80-91% of small-prime columns), which voids four headline consequences; "every coordination mechanism excluded" overclaims. |
| **pal11cert** | 8 | 0 | 0 | **3** | Theorem 4.3 is exact and independently re-proved by a simpler route; it also finds and fixes a real round-13 record error (non-exhaustive `prodsearch.py`; `zeta(M8') = 10`). Deductions: it answers the paused route-5 item, not this round's question; it holds identically for fakes and actual clusters, so it has no barrier content; the TC question it was meant to inform stays open. |

Partial (c) credit is for corrections to the barrier map only: the restated Rig^SL question, the `(eps, t)` and
`r >= M^gamma` corrections to the transfer theorem, and the removal of the known C1 seeds. No (c) score exceeds 1
because no statement about actual circles moved.

## 8. Single most promising next step

**Recommended: decide Rig^SL, the norm-`<= K` truncation of P3-all\* itself (labels and one angle vector), on the
branch that this round left open: whether loose primes stall once the angles are imposed.**

**Why this one, from the whole record.**
1. P3-all\* in regime (1.4) is still the only local input class that is valid for actual circles, needs no
   conjecture and is not known to be capped (README growth-rate status; rounds 7-9, 13). Route 5 is paused (round
   13 §7.1; pal11 closes its cheapest certificate test). Vojta stays conditional.
2. This round shows that the decisive finite question was posed in the wrong class. In PU it looked answered
   ("growth for tight primes, stall for loose primes"). In SL the growth half survives (`K = 4 -> 5` strict for
   `first4`), while the stall half and every C1 seed fall away. The slab coupling itself acts as a rigidity
   mechanism: it turns the PU-flexible moduli 83 and 163 rigid and moves the loose one-modulus segments up (rand4
   to 31, bal4 to 41). So the branch that decides "C1 or rigidity" is open again, and it is now finite and
   well-posed.
3. Either answer is decisive for the barrier map, as in round 13, but now with the right class:
   * **Loose primes stall in SL** (non-coherent systems with a verified `phi` on most forced moduli, at every `K`
     tested): these are the first genuine C1 seeds, since labels and angles are both already consistent at norm
     `<= K`. What remains is DeepBal at larger norms, the pair conditions, and the deep classes (Thm S handles the
     shallow ones).
   * **Loose primes become rigid in SL**: rigidity holds across prime types, the SL coupling is the mechanism, and
     the target becomes a proof of Rig^SL on `(M^(1-eta), M]` for `r >= M^gamma`, feeding the corrected Theorem
     4.1.

**Concrete test (exact where it matters).**
* Class: parity-respecting labels at all levels on `F_K^ex` (free moduli collision-free) and one `phi in R^r`,
  satisfying (1.2) for every `1 <= ||v||_1 <= K` and every unit. Report cores in this class only.
* Certificates (rigidity): keep the Lemma R relaxation and add the exact rational box test of Lemma J0 per `v`.
  That is still a finite signature problem, and it already decides 83/163 for `first4`.
* Witnesses (flexibility): every witness must carry an explicit `phi` checked against (1.2) in interval arithmetic
  with a rigorous margin, or exactly at `phi = arg pi_j` (where (1.2) is `m_(v,s) <= |D_s(Z_v)|`). Search near
  `phi = arg pi_j` first: there the only slack is the uncharged `(2N)`-part of `D_s(Z_v)` (rigtheory Lemma 3.1),
  so candidate deviations can be enumerated from the factorisations of `D_s(Z_v)`.
* Instances: `last4`, `rand4`, `bal4` at `K = 4, 5, 6` first (the open branch); then `first4` at `K = 3, 4` to pin
  `Q_core^SL` (now `>= 37` and `[61, 103]`); then `r = 6..8` to test whether the `r`-saturation at 179 persists.
* State with every result: prime set, `r`, `K`, the forced set, the level depths used (`X >= 27, 49, 121, 169`),
  and whether it is exact, solver or evidence.

**The theory question that goes with it.** The ceiling `q_core(K) <= O_K(r^floor(K/2))` means that bounded-`K`
rigidity reaches the scale `M` only when `r >= M^gamma`. Decide whether near-extremal clusters in regime (1.4)
(`W <= c_0 M log M`, `p_j <= M^(h_0)`) must have `r >= M^gamma`. If they must not, the route-3 chain needs a
separate mechanism for few primes with large exponents, and Rig at bounded `K` cannot supply it.

**Not recommended.** Further per-unit (PU) Rig(K) computations; any flexible-top "seed" that has no verified `phi`;
P3-G pinning closures (round 7: they decide row 6 only, not growth); more route-5 certificate searches.

## 9. Records to correct and files

**Corrections to earlier records.**
* `round13/README.md` §7.2: Rig(K) as worded ("meet the per-unit budget") is a strict relaxation of P3-all\* at
  norm `<= K`. Replace it by Rig^SL (§8). The plan "take the flexible tops as C1 seeds and add the slab conditions
  with angles chosen after the labels" fails for every top found (§3).
* `round13/README.md` §1.4 item 7 and `round13/fatpoint.md` §5.2: `prodsearch.py` was not exhaustive. M8' (syl8,
  `L = 3`) at `m = 10` has 1092 zonotope certificates, so `zeta(M8') = 10`. The M8 `m = 8` negative becomes "42 pass
  the counts, 0 are certificates". The `m = 7` certificate counts are 300 (pal7), 1746 (M8), 24 (sdc7).
* `round10/README.md` §2 (rigid core `{3..47}` at `K = 4`): with all moduli free the coherent PU core is `[3, 61]`.
* The class symmetry group is `{+-1}^r x (Z/2)^r` (conjugation and associate twist), not `{+-1}^r`. Any "strict"
  notion must be read up to `(eps, t)`.

**Files** (all under `round14/`).

| file | content |
| --- | --- |
| `rigcompute.md`, `rigcompute_checks/` | Rig(K) finite test (PU class) |
| `referee_rigcompute.md`, `referee_rigcompute_checks/` | class mismatch I1, SL results, MILP check |
| `rigtheory.md`, `rigtheory_checks/`, `rigtheory_notes.md` | Rig(K) theory, one-modulus and joint data |
| `referee_rigtheory.md`, `referee_rigtheory_checks/` | associate twist (Prop. T, Lemma T'), exact joint certificates |
| `c1construct.md`, `c1construct_checks/` | C1 attempt: B1-B3, pinning, L/G/C, R, S |
| `referee_c1construct.md`, `referee_c1construct_checks/` | condition (i) failure, corrected B1 |
| `pal11cert.md`, `pal11cert_checks/` | pal11 product-certificate theorem, round-13 correction |
| `referee_pal11cert.md`, `referee_pal11cert_checks/` | independent re-proof of Thm 4.3 |
| `judge_checks/` | J1 `jang.py` (+ `jang_sanity.py`), J2 `jonemod.py`, `jonemod_fast.py`, `jlift.py`; helper `jglib.py`; outputs `*_out.txt` |

**Repository note.** Generated data files over 1 MB, pickles (ignored by `.gitignore`) and compiled
binaries (`*.so` and executables built from the `.c` sources next to them) are not committed. The scripts
regenerate them; compile the C helpers with `cc -O2` (shared objects with `-shared -fPIC`) as the Python
drivers expect. Two stray agent notes found in the repository root during the round are kept as
`rigtheory_checks/stray_rigtheory_notes.md` and `rigtheory_checks/jointG_K3_rest.txt`.
