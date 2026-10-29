# Round 27 — judge's report (new-mathematics round: dp2, stab, split, interp)

Judge: reads `dp2.md`, `stab.md`, `split.md`, `interp.md` and the four kill reports `kill_*.md`; own checks in
`judge_checks/` (all **unrefereed**). Labels: **proved** (proof read and re-derived by the killer; not refereed by a
third party unless said), **conditional**, **exact** (finite exact computation), **numerical**, **evidence**,
**heuristic**, **unrefereed**, **open**. "Killer" = the adversarial checker of that lane.

## 0. Verdict

**No lane moved the growth rate or the uniform bound.** `M <= (2/3+eps) log R/loglog R` (two referees) stands, and no
`M = o(log R/loglog R)` statement was proved, conditionally on anything weaker than what the record already had, or
unconditionally. Measured by the round's own bar ("rediscovering a mapped barrier is a failed outcome"), three of
four lanes failed: the killers rated dp2, stab and split **repackaged known barriers**, and interp **genuinely new,
modestly**. I agree with all four verdicts after reading the proofs and checks.

What is new **and** correct (none of it is a mechanism that beats a mapped cap on the goal problem):

| item | lane | label | goal relevance |
|---|---|---|---|
| Thm 3.3: a mean-width-extremal measure is unique up to scalar on its support (`dim A_t(S) <= 1` for integer `t >= MW_n(S)`) | stab | proved (killer re-derived; 82 703 sets; judge J1: 2 130 more sets, 913 tight, 0 violations) | closes the two-section part of the mapped loophole "equality case of Theorem B" (construct.md §3 item 3) and the uniqueness half of Conjecture E. Constant-level only. |
| Thm 3.1: multi-shift zero criterion (`J(A)_e = cap_l J(A cap (A-l))_e` for `e < m` ⇒ no `Y` with all `Y^l != 1` in `Bs V(A,m)`) | interp | proved, elementary; new to the record; literature status unclear (possibly folklore) | sufficient criterion for route 5's ZE target (round 11's ZE-braid); Cor 4.1: `(SHIFT)_A` ⇒ `M = O_C(log log R)` (conditional). Hypothesis is at least as hard as the open lemma ZE^all_H. |
| Thm 3.1 + killer completion: the mnn31 (`C = 7/2`) five-divisor curve sits on a quadric cone at `alpha = 7/2`; every admissible split `Z` contains 15 cones, 15 cubic threefolds `iota(Q_sigma)`, and `W` (+ `iota W`) | split | proved, exact | calibration of the split exceptional set (and of absolute Erdős #887). Split-only; does not port to `T_5`. |
| Prop K (killer): at the even-set vertex `x_[5] = 1[|B| even]` (and its inversion) no rational curve has `mult_e >= 4` | split (killer) | proved by the killer; judge J2 confirms the linear algebra on 8 instances | closes the `N = 1` curve route to the split cap WTB there. Split-only. |
| Transfer lemma (killer): on `Bl_e Lbar_k`, `(K + sum Zt_ab + E).lt = -1` on lines through `e`, every `k >= 3`; `tau_e(2(-K) - B_{1,4}) = 5` exactly | dp2 (killer) | proved | integrality w.r.t. the diagonals never carries a Lang–Vojta prediction at `e`, for any `k`. Barrier (negative). |
| certified dP2 instance `(S,P)`: `Pic = ZK`, `tau = 3/2`, `beta = 17/18`, real-cubic branch, no `Q`-curve of degree `<= 6` with `alpha < 2`; exact two-point effective cone of `Bl_{P,P'}S` | dp2 | proved / exact | calibration only; every unconditional statement is in the T4 (section) class. |
| Castravet–Tevelev link: `Bl_e LM_5` is a contraction of `Mbar_{0,7}` | split (killer) | derived by the killer; literature seen via search results only | `tau_e = 10/3`, Theorem B, WTB and Conjecture E are statements about `Eff(Mbar_{0,7})`; new to the record. |

**Honest summary.** The round produced correct new mathematics *around* the barrier — one genuinely new lemma
(stab Thm 3.3), one genuinely new elementary criterion (interp Thm 3.1), one new calibration family (split cones), one
new obstruction (Prop K), one uniform negative transfer (killer's line lemma), and a literature link (Castravet–Tevelev)
— but no mechanism that gets past any mapped cap on the goal problem. No lane produced `alpha^ess > tau` at any
rational point with `tau < dim`, on any variety, under any hypothesis weaker than Vojta.

## 1. What the round means for the goal

**Growth rate.** Unchanged. The only growth-relevant statement is interp's Corollary 4.1:

> **(conditional; implication proved)** If the cluster's `T = 1` slab box `A` satisfies `(SHIFT_e)` w.r.t.
> `L_pair = {±e_x, ±(e_x - e_y)}` for all `e <= ceil(4w) - 2`, plus the top-degree height condition (3.3) at
> `H = U^2`, then `n - 1 < 316 ln 4r`, hence `M = O_C(log log R)`.

The hypothesis implies the height-cut ZE-braid at `ceil(4w)`, which is strictly stronger than the record's open
lemma ZE^all_H and at least as hard as the exact barrier question (killer §1). It is proved only for `e <= j - 1`,
and it is *false* at Sylvester orders `> 2j log_2 n` (Prop 3.4 = the round-11 referee's ZE-braid refutation). Model
evidence is mixed: in 27 random multi-cut slab boxes the criterion's reach never exceeds the zonotope (Liouville-
capped) order (below it in 20), but on M8' (syl8, `L = 3`) it reaches the truth `rho = 11` while the zonotope order is
10 (killer §5.2(d); mod `p`, two primes). If "reach <= zonotope order" persisted at regime size, Corollary 4.1 would be
vacuous. This is a new *form* of route 5's missing lemma (linear algebra on one explicit set per profile), not a
weaker one. The hypothesis is not comparable with Vojta; it is not "a conditional result under a hypothesis strictly
weaker than Vojta" in any useful sense.

**Uniform bound.** Unchanged. Question B is open in every form. What the round adds:
* *Stability reading.* Near-tight balanced stability `Stab(eta, eps)` is equivalent to `(L_5°)`; `A_5(1/30 - eta)`
  only implies it (stab §1.3, killer-confirmed).
* *Pencils.* "Second-order integrality" at Liouville-tight points, in its only available form (two tight sections
  with point-dependent coefficients), is void: tight pencils do not exist (stab Thm 3.3 / Cor 3.4). The killer is
  right that this is the section cap at its boundary face, not a mechanism outside it.
* *Integral points.* Integrality w.r.t. the diagonals (bounded pair cofactors) has no Lang–Vojta content at `e` for
  any `k` (killer's line lemma), and amplifies nothing beyond sections wherever the weighted cap WTB holds (split
  Lemma 3.5). On the dP2 the analogous pair *is* of log general type (dp2 Thm C). So balance must enter any proof on
  `X`; integrality and proximity alone carry nothing.
* *Profiles.* Non-simple (chain-carrying) balanced profiles beat 15/4 by their chain deficit (stab Prop 3.5,
  proved, small novelty); asymmetric split profiles do not help symmetric sections (split Thm 3.2). Both locate the
  hard core at **simple, symmetric, balanced profiles**, which is exactly where the twisted problem sits and where
  all recorded fakes live.
* *Split side.* The split exceptional set must contain at least `W`, `iota W`, 15 quadric cones and 15 cubic
  threefolds; any proof of split Question B (hence of absolute #887) must allow all of them.

## 2. Lane dp2 (`dp2.md`; `kill_dp2.md`) — survives; novelty: repackaged known barrier

**Killer-corrected statement.** `S : w^2 = f` (the write-up's explicit quartic), `P = (0:0:1:1)`.
* **[proved / exact]** `S` smooth (over `Q`, `F_3`, `F_5`); `P` off the 56 lines and the ramification curve;
  `Pic(S) = ZK` (Frobenius at 3 on `NS` has char. poly `(X-1)(X+1)Phi_18`, a Coxeter class of `W(E_7)`, from
  exhaustive counts over `F_{3^n}`, `n <= 8`); `eps_P(-K) = 4/3`, `tau_P(-K) = 3/2`, `beta_P(-K) = 17/18`; the
  `tau`-computing member `D in |-2K|` (triple point at `P`) is unique, irreducible, rational, with branch field the
  cubic field of `19t^3 - 10t^2 + 20t + 8` (disc `-2^10·1667`, one real place), so `alpha_{P,D} = 2`; no `Q`-curve of
  `-K`-degree `<= 6` through `P` has `alpha < 2`. Hence `3/2 <= alpha^ess_P(-K)`, Vojta predicts `>= 2`.
* **(A) [proved; T4 class]** The `Q`-effective cone of `Bl_{P,P'}S` (`P'` the Geiser image) is spanned by the fibre
  `F` and `Gamma_k = k sigma'` (killer closed the cone step: `NE ∩ F^perp = R_+F` via Zariski decomposition). Liouville
  with forced divisibility at `P, P'` gives exactly `x <= x_sec(nu_+, nu_-)`; this is `1/tau = 2/3` when
  `gcd(a,b) = 1`, and better only if `nu_+ > 0` or `nu_- > 2/9` (domain heuristically no richer than the window).
* **(A3) [proved; mechanism = McKinnon–Roth MR-L Thm 4.5, uncited]** On `S_0`, `alpha_{P_0} = tau = 3/2` is attained
  on a cuspidal tangent-plane cubic whose points are fibre-locally generic; no bound uniform in fibre data alone beats
  `3/2`. (`S_0` is special: `Q`-lines, rational branch of `D_0`, `f_1 = f_2 = 0`.)
* **(C) [conditional on Lang–Vojta for `S \ D`]** For each `K`, `{X^{3/2}H <= K}` lies in a proper closed `Z_K`
  (depending on `K`). Weaker than `alpha^ess > 3/2` and weaker than ERR Prop 0.1 / LMS Prop 1.4 (Vojta ⇒ `>= 2`).
  Subspace route collapses to McKinnon–Roth's `beta = 17/18` (threshold is `3/2`, not 1).
* **(Transfer) [proved, every `k >= 3`, by the killer]** see §0. `(UE)` is a restatement of the window, not a
  reduction; §8.1 of the write-up is a sketch.

**What survived as new.** A certified instance and an exact two-point cone (map value); the killer's line lemma.
**Mechanism:** none outside T4/Vojta. **Growth / uniform:** none.

## 3. Lane stab (`stab.md`; `kill_stab.md`) — survives; novelty: repackaged barrier closed by a new lemma

**Killer-corrected statement.**
1. **[proved; new to the record]** Thm 3.3: for finite `S ⊂ Z^n` and every integer `t >= MW_n(S)`,
   `dim A_t(S) <= 1`. Proof: slice peeling (Theorem B's), keeping the `P' = ∅` term, which rules out two critical
   slice indices. Through the Cayley transform (referee_construct.md §2) it is uniqueness of extremal measures for the
   multiplicity-Schwartz–Zippel / mean-width bound. Literature: no equality characterisation found (killer, WebSearch).
2. **[proved]** Cor 3.4: on `Bl_e` of either form of the Losev–Manin 4-fold, `h^0(pi^*D - tE) <= 1` when
   `8t >= D.beta`. Vacuous for multiples of `-K` (`8t >= 30m` exceeds `floor(10m/3)`); content on the face
   `gamma_hat^perp`, whose effective part in all samples is the cone of the 10 diagonal divisors (Conjecture E's
   structure half, **open**). Minimal pencil defect `s = 4` is **evidence**; `s >= 1` is proved.
3. **[proved; barrier]** Pencil principle (Thm 3.2) is the section cap at its boundary face (its inputs are
   single-section Liouville cofactors, gcd-invariant by Lemma K1). It closes, unconditionally, the two-section part of
   the mapped loophole "constant-level effects in the equality case of Theorem B". The claim that second-order
   integrality "can be produced in two ways only" is **heuristic**.
4. **[proved; restricted; small novelty]** Prop 3.5: near-tight balanced one-class points with chain deficit
   `>= (s_S + c)w` (some `S` in a fixed finite list) lie on `{f_S = 0}`. Satisfiable, but needs chain primes carrying
   `≳ 2w` of log-weight; void on squarefree profiles (where fakes live). Makes explicit round25/intra.md 523–524.
5. **[proved]** Elementary forms: Cox parametrisation is a bijection exactly on the chain locus (Ch);
   `Q_ab Q_bc = G^(b)_abc Q_ac` for pair cofactors (not reduced chords); `x_b t_a - x_a t_b = G_ab t_ab`;
   `Stab ⇔ (L_5°)`, `A_5(1/30 - eta) ⇒ Stab`; arccot form exact only for `t = 1`; sign fix
   `x_a ≡ iota(R_a) t_a (mod N(R_a))`. (These verify, with two corrections, the coordinator's elementary form.)
6. **[proved; no gain]** Fibrewise section arguments = direction-sensitive one-section Liouville (uncapped class;
   every computed instance exactly tight); four singleton conditions decouple at fixed shared blocks (tautology).
7. Numerical corrections: split chain deficit reaches 11.15 over all four Pell orbits (not `<= 6.8`); `k = 4`
   `|c_ab|^2` reaches 230 400 at random `s` (true bound 705 600). Conclusions unchanged.

**What survived as new.** Thm 3.3 (the round's best single theorem), Cor 3.4, the corrected elementary forms.
**Mechanism:** the pencil principle would have been one, but it is provably void. **Growth / uniform:** none.

## 4. Lane split (`split.md`; `kill_split.md`) — survives with corrections; novelty: repackaged barrier + one new calibration

**Killer-corrected statement.**
1. **[proved; exact]** Thm 3.1: the mnn31 `C = 7/2` curve has `B.C = 7` (7 distinct cuts once each) and two conjugate
   real-quadratic branches of order 1, so `alpha = 7/2`; its 3-parameter conic family is dense in the quadric cone
   `(u2-u1)(u3-u1) = (u4-u1)(u5-u1)`. Every admissible split `Z` in `A_5^split(delta)`, `delta < 1/14`, contains the
   15 cones, their 15 inversion images `iota(Q_sigma)` (cubic threefolds; the co-divisor tuples lie on no cone), and
   `W`, `iota(W)`. In twisted `T_5` only `Q ∩ sigma Q` can carry rational points.
2. **[proved]** Cost identity and Thm 3.2: `Pi_nt` (11-dimensional) is the polytope of cut-balanced curve classes; for
   centrally symmetric supports the forced-divisibility cost is profile-independent, and `cost >= 4 reg`
   (`cost_{1/2} = (1/4) sum_P w_P`; the write-up's `1/8` is a slip).
3. **[standard]** Prop 3.7 = easy half of BDPP: WTB at `x` ⇔ `pi^*gamma_x - 4l` movable on `Bl_e LM_5`, a contraction
   of `Mbar_{0,7}` (Castravet–Tevelev, Duke 2015, arXiv 1311.7673 Thm 1.1; search results only).
4. **[numerical, reproduced independently at 60 digits]** Dominant four-branch rational families at 14 of 16
   one-sided classes, spanning `Pi_nt`: WTB holds numerically on an 11-dimensional neighbourhood of `x = 1/2`. Not
   interval-certified.
5. **[proved by the killer]** Prop K (see §0); WTB at `x_[5]` is **open** (needs `N >= 2`, positive genus, or a
   divisor-side proof). If WTB fails there, a section beats 15/4 near the even profile — split-only.
6. **[proved]** Lemma 3.5 (Delta-poles gain nothing where WTB holds); the claimed reduction of genuine fixed-`S`
   Subspace on `LM_5 \ ∪ Delta_ab` to Lemma 3.5 is **withdrawn** (the S-unit equation is a counterexample to the
   principle). Thm 3.6 (`x^tau × y^tau = (t_bc, t_ca, t_ab)`, exact).
7. **[derived]** Phase 0: forum reduction = complete intersections of quadrics (K3 at 4 divisors, general type from
   5); this is Vojta's Büchi paradigm. Known families are Pell conics.

**What survived as new.** The cone calibration (with its inversion completion) and Prop K. **Mechanism:** none; the
section cap carried to asymmetric profiles. **Growth / uniform:** none; nothing ports to `T_5`.

## 5. Lane interp (`interp.md`; `kill_interp.md`) — survives; novelty: genuinely new, modest

**Killer-corrected statement.**
1. **[proved; new to the record]** Thm 3.1 (multi-shift criterion; see §0) with tangent-cone and height versions. With
   `L = L_pair` its conclusion is round 11's ZE-braid(A,m) (uncited in the write-up); Prop 3.4 is the round-11
   referee's ZE-braid refutation at Sylvester orders, transported.
2. **[proved, implication only]** Cor 4.1 (see §1). Hypothesis **open**: proved for `e <= j - 1` (Prop 3.5), false at
   Sylvester orders `> 2j log_2 n`.
3. **[proved; repackaging]** Thm 2.2: under the gap hypothesis every nonzero form of `V(A,m')` has a Zariski-dense set
   of nondegenerate zeros, so no single form proves ZE^all (the record's zonotope cap plus a standard classification).
4. **[evidence, mod p unless exact]** Reach w.r.t. `L_pair`: sdc7 `L = 2`: 6 (zonotope order 7, exact); syl8 `L = 2`:
   6, or 7 with all cube shifts (zonotope 7; ZE^all true at 8); M8' = syl8 `L = 3`: 11 = rho = truth (zonotope 10,
   exact), reproving round 11 Thm 5.2(a) by linear algebra; 27 random slab boxes: never above zonotope order. The
   write-up's "reaches exactly the zonotope order" and "the failure is genuine" are **false** as stated.

**What survived as new.** The criterion itself and the M8' data point (killer's). **Mechanism:** yes — a system-level
(non-single-form) certification that passes T4 formally and is non-zonotopal on M8'. It is the only round-27 item that
is both outside every mapped cap and attached to a growth statement. **Growth / uniform:** none proved.

## 6. Pre-registered tests T1–T5, consolidated (killer reruns; judge agrees)

| test | dp2 | stab | split | interp |
|---|---|---|---|---|
| T1 fakes (growth claims) | vacuous (no growth claim); `S_0` analogue passes, "statistically identical" overstated | vacuous; Prop 3.5 void on squarefree profiles | vacuous | passes: rigidity is fake-violated; `(SHIFT)_A` depends on the profile only, and Thm 3.1(3) forces it to fail on any exact fake |
| T2 `k = 4` control | passes: line lemma ⇒ integrality pair predicts nothing at `k = 4` (`tau = 2 < 3`) | passes vacuously (pencil route void at every `k`); residues bounded on the Thm 4.1 family | passes: `k = 4` cap 7/2; conditional exclusions do not fire | not applicable (no Question-B claim; regime needs `n >= 2966`) |
| T3 split / #887 | vacuous | no #887 claim; five-divisor curve not excluded | m-czech on `W`, mnn31 on the cone; both in `Z`; cubics added | not applicable; flagged: any transplant would be Galois-blind and would have to prove absolute #887 |
| T4 criticality | every unconditional statement is in the T4 class (Thm A *is* `tau`) | Thm 3.2 overstated as outside the cap: it is the cap's boundary face | inside the section class throughout | passes formally; mostly at or below zonotope order, one exception (M8') |
| T5 calibration | `D`-points, `E_0`, `C_0` all respected | round-24 family, Pell curves, 8/10-point families untouched | `W`-points protected; cones added to `Z` | criterion stops before every known nondegenerate base point (sdc7 at 8, 1-D at 5, Sylvester subtori) |
| owner constraints | respected (`Z_K` explicitly `K`-dependent) | respected | respected after withdrawing the fixed-`S` Subspace principle | respected |

## 7. Is the dP2 calibration the right 2-dimensional analogue? What would it take there?

**Yes for the literature form, no for the form the problem needs.**

*Where it is right.* The literature question behind Question B is "is any `alpha^ess > tau` provable at a rational
point with `tau < dim`?" The dP2 point is the cleanest such instance: dimension 2, `Pic = ZK`, `tau = 3/2 < 2`,
generic (`P` off lines and ramification, real-cubic branch), with the whole section/two-point story computed exactly
(Theorem A). If an expert knows an unconditional `alpha^ess_P(-K) > 3/2` at such a point, that is a mechanism of a kind
the record has never seen, and it should be studied for transfer. The papers check (§0, §3 item 5) and the dp2
literature pass both found none.

*Where it is wrong.* Question B does not need `alpha^ess > tau` on `X`; `tau_e = 10/3` and the cap that must be beaten
is the **balanced** Liouville value 15/4. Two facts make the dP2 a poor *mechanism* analogue:
1. `Pic(S) = ZK` has no notion of balance, while three lanes this round (dp2 transfer, stab Prop 3.5, split Thm 3.2)
   and the record (Lemma K1) all point to balance — simple, symmetric, balanced profiles — as the hard core.
2. The structure differs at the top of the window: on the dP2, near-tight points are integral points of a
   **log general type** pair (Thm C), so Lang–Vojta for one surface pair already gives stability; on `X` the analogous
   pair is log-uniruled through `e` for every `k` (killer's line lemma). The dP2 is conjecturally *easier*.

*There is no lower-dimensional balanced analogue in the Losev–Manin tower* (judge J3, exact arithmetic; it is the
record's sign of `kappa` restated). On `Lbar_k^tw` the balanced Liouville constant at `e` is
`alpha_bal(k) = 4(1 - 2^{1-k}) = (2^{k-1}-1)/2^{k-3}`: `3, 7/2, 15/4, 31/8, …` for `k = 3, 4, 5, 6`, against
`dim = k - 1 = 2, 3, 4, 5`. It is below the dimension exactly for `k >= 5`. So `k = 5` is the smallest case of the
window in this tower, and for `k = 3, 4` balanced Liouville is attained densely (gauss Thm 4.1), consistent with Vojta.
A faithful 2-dimensional analogue would need a surface with Picard rank `>= 2`, a rational point with
`tau < 2` *and* a "balanced" sub-cone of height directions whose own Liouville constant is `< 2`. I do not know one;
finding one is a well-posed small question (not attempted; toric surfaces are excluded because Huang's results give
`tau >= dim` there).

*What it would take on the dP2* (all **open**):
* an unconditional non-density theorem for integral points of `(Bl_P S, Dt + E)`, i.e. on `S \ D` with one
  irreducible rational boundary curve `D in |-2K|` — beyond Corvaja–Zannier/Levin-type results, which need several
  boundary components (memory-quote; the killer agrees, kill_dp2 §5.2); this gives only the `K`-dependent stability
  `Z_K`, not `alpha^ess > 3/2`;
* or `(UE)` for one explicit pencil: `log c <= (3 - eta) log A'` for `gcd(a,b) = 1` points off a fixed curve — which
  is the window statement itself;
* or a refutation: a dense supply of `Q`-curves through `P` with `alpha < 2` (e.g. infinitely many with `alpha < 7/4`),
  which would contradict Vojta.

Honest probability that anyone proves `alpha^ess_P(-K) > 3/2` at this dP2 point within a year: **about 1%**. Value of
the calibration to the goal: map value (a sharp expert question), not a route.

## 8. Cross-lane ideas

1. **One cone-theoretic dictionary (stab Thm 3.3 + split BDPP + Castravet–Tevelev).** `tau_e(-K) = 10/3`, Theorem B,
   WTB and Conjecture E are all statements about `Eff(Bl_e LM_5)`, hence about `Eff(Mbar_{0,7})` pushed forward under
   the Castravet–Tevelev contraction. Under it the 10 diagonals `Z_ab` are the light–light boundary divisors and `E` is
   `D_{0∞}`. Conjecture E's open structure half becomes: *the face of `Eff(Mbar_{0,7})` orthogonal to the balanced
   class `gamma_hat = pi^*beta - 8l` (nonnegative on effective divisors by Theorem B, so in the closed movable cone of
   curves) is spanned by boundary divisors.* That is a question about a known, studied cone
   (hypertree and non-hypertree extremal divisors on `Mbar_{0,7}`; memory-quote, not checked) and can go to
   Castravet or Tevelev. **Goal value: none** (constant-level loophole only); map value: real.
2. **One integrality barrier (dp2 Thm C + split Lemma 3.5 + killer line lemma).** Three lanes, three proofs, one
   conclusion: integrality w.r.t. the diagonals / bounded pair cofactors carries no information at `e` beyond
   sections — log-uniruled complement (Lang–Vojta predicts nothing, every `k`), and Vandermonde ratio exactly 15/4 on
   `Pi_nt` (Delta-poles amplify nothing where WTB holds). Only genuine fixed-`S` Subspace on `LM_5 \ ∪ Delta_ab`
   remains formally open, and by the line lemma it cannot use proximity alone. Proposed as one barrier-map item
   (§12).
3. **The hard core is the same in three lanes.** stab Prop 3.5 (chains help), split Thm 3.2 (asymmetry is invisible
   to symmetric forms) and the record's simple-profile caps all locate the obstruction at *simple, symmetric, balanced*
   profiles, where the twisted problem sits and where every fake lives. Any future Question-B mechanism should be
   tested there first; nothing else is informative.
4. **stab Thm 3.3 → interp (trivial, judge's observation, proved).** If `t = MW(A)` is an integer, `V(A,t)` has
   dimension `<= 1`, so its base locus is the zero set of one form; under interp Thm 2.2's hypothesis that set has
   nondegenerate points, so ZE^all fails at order exactly `MW(A)`. Route 5 must work at `m' <= MW(A) - 1`. I do not
   know of a regime box where this binds (the recommended `m' = ceil(4w)` is far below `MW(A)` heuristically); recorded
   only for completeness.
5. **interp → Question B: no feed.** One-point Liouville at `k = 5` needs non-vanishing off a fixed `Z`, which is
   automatic; base-locus criteria matter only for route 5 (growth). Conversely the dp2 pencil and stab's fibrations
   land in the one class named uncapped by two lanes, *direction-sensitive one-section Liouville*; at
   `eps`-balanced points the approach direction `(t_a5 e^{-d_a5/2})_a` varies by factors `e^{±O(eps w)}`, so I expect
   (heuristic) at most an `O(eps)` gain there. Not pursued.

## 9. Judge's own checks (all **unrefereed**; scripts and outputs in `judge_checks/`)

* **J1 (stab Thm 3.3; exact, own code `j1_rigidity.py` → `j1_rigidity.out`).** Exact ranks (python-flint `fmpq_mat`)
  of moment matrices on 2 130 sets in `Z^1..Z^4`: Minkowski sums of root binomials (`e_i`, `e_i - e_j`), the same with
  one or two extra points or one point swapped, and random small sets. `dim A_{ceil(MW)}(S) <= 1` in all cases;
  913 sets are tight (`A_MW != 0`, `MW` integral); 0 violations of Thm 3.3 or of Theorem B. Independent of the
  author's and the killer's code. (In this sample `A_{ceil MW} != 0` occurred only when `MW` is an integer.)
* **J2 (killer's Prop K, steps 4–6; exact, `j2_propK.py` → `j2_propK.out`).** For 8 random rational configurations of
  five lines `w = rho_a y + kappa_a` (distinct slopes, no three concurrent, none through `O`): the quartics through the
  10 nodes form a 5-dimensional space spanned by the five products of four lines; those with a triple point at `O`
  form the zero space. With four concurrent lines the triple-point space is 4-dimensional (the proof's only escape,
  excluded by non-degeneracy). Steps 1–3 (normal form) were checked by hand, not by computer.
* **J3 (balanced Liouville vs dimension; exact arithmetic, `j3_balanced_liouville.py` → `.out`).**
  `alpha_bal(k) = 4/(1 + 2 delta_L(k)) = (2^{k-1}-1)/2^{k-3}` with `delta_L = 1/(2^k - 2)`; `< k - 1` iff `k >= 5`;
  `kappa(delta_L) = -1/3, -1/7, 1/15, 9/31` for `k = 3..6`. Used in §7.
* **Reading checks.** I re-derived Thm 3.3's critical-index step (the two (∗) inequalities sum to `j - 1 < j - 1`),
  the killer's line lemma (`K.lt = -k + (k-2)`, `Zt_ab.lt = 0`, `E.lt = 1`), and Cor 3.4's vacuity for multiples of
  `-K` (`30m/8 > 10m/3`). No discrepancy with the kill reports.

## 10. Scores

Rubric of `round10/README.md` §1 (as in rounds 22–26): (a) correctness after the kill, (b) fake / barrier status of
the input, (c) movement of growth or uniformity; the overall score weighs (b) and (c) above (a). As in rounds 22–26,
no movement of growth or uniformity caps the overall score at 3.

| lane | (a) | (b) | (c) | overall | reasons |
|---|---|---|---|---|---|
| **dp2** | 6 | 2 | 1 | **2** | All exact claims reproduce; Thm A, A1, A3, B1, C correct. Every unconditional statement is a section/`tau` statement; A3's mechanism is MR-L Thm 4.5 (uncited); Thm C is a Vojta consequence weaker than ERR/LMS. Value: a certified calibration and a precise picture of what the 2-dimensional window needs. |
| **stab** | 7 | 3 | 1 | **3** | Thm 3.3 is correct, general (all `n`) and genuinely new to the record — the round's best theorem. The "new barrier" is the mapped equality-case loophole; T4 overstated; small numerical slips. |
| **split** | 5 | 2 | 1 | **2** | New exact calibration (cones; completed by the killer with the cubics), correct cost identity; but a false display, an overclaimed Subspace principle, an incomplete `Z`, a forward claim refuted by the killer (Prop K), and missed literature (Castravet–Tevelev, BDPP). Split-only. |
| **interp** | 6 | 3 | 2 | **3** | New elementary criterion, outside every mapped cap, non-zonotopal on M8', attached to a growth corollary; but two false model claims, its target is round 11's uncited ZE-braid, and its hypothesis is at least as hard as ZE^all_H. |

Round as a whole: **2.5 / 10** on the goal; the round's own bar ("new mathematics, not a rediscovered barrier") was met
only by stab Thm 3.3 and interp Thm 3.1, and neither reaches the goal.

## 11. Single most promising next step (with honest probabilities)

**Decide, by computation, whether the multi-shift criterion outruns the Liouville-capped zonotope order as slab boxes
grow.** It is the only round-27 mechanism that is simultaneously (i) outside every mapped cap (T4 passes; it is not a
single form), (ii) non-zonotopal on at least one slab box (M8': reach 11 = truth > zonotope order 10), and
(iii) attached to a growth statement (Cor 4.1: `M = O_C(log log R)`). Concretely:

* On the round-11/12 slab-box families (Sylvester and sdc profiles, `T = 1`) at `d = 7, 8` (and `9` if memory allows,
  under `prlimit --as=3500000000`), for each admissible `L`, compute mod two primes the reach of `(SHIFT_e)` w.r.t.
  `L_pair` and w.r.t. all `{-1,0,1}^d` shifts, the top degree `rho`, and the exact zonotope order `zeta` (DFS as in
  `kill_interp_checks/k6_zono.py`, certified with an exact `V(A,zeta) \ V(A,zeta+1)` check). Tabulate `reach - zeta`
  and `reach/rho` against `rho/L`, and add the multi-cut random boxes at `d = 8`.
* Decision rule fixed in advance. *Favourable:* `reach > zeta` on a growing share of multi-cut boxes with
  `reach - zeta` increasing in `rho`; then attack `(SHIFT)_A` by proof via nested shift chains `S_e ⊂ S_{e+1}` (interp
  open point 3) and the linkage form `H_A(e) = T_d(e) - H_{cube\A}(2d - e)` (open point 4). *Unfavourable:*
  `reach <= zeta` on every multi-cut box beyond M8'; then Cor 4.1 is presumptively vacuous, the interpolation route to
  route 5 is closed, and internal rounds should stop with Questions A, B and the lemma ZE^all_H sent to experts (the
  standing recommendation since round 22), now accompanied by the certified dP2 instance (§7) and the
  Castravet–Tevelev reformulation (§8 item 1).

Honest probabilities (judgement):
* the computation is feasible at `d = 8` within the memory limit: about 70%;
* the trend is favourable in the sense above: **about 15%**;
* conditional on that, `(SHIFT)_A` is proved on regime boxes within a few rounds: about 5–10%;
* hence a proved growth improvement `M = o(log R/log log R)` along this route: **about 1%**;
* Question B (uniform bound) proved by any route in the next year: **< 1%**; an expert pointer to an unconditional
  `alpha^ess > tau` at some rational point with `tau < dim`: about 3–5% (unchanged from papers_check §4).

## 12. Proposed barrier-map updates (text only; `expert_note.md` not edited)

1. **§5 item 3, new bullet (stab Thm 3.3, Cor 3.4; proved, killer-checked, unrefereed).** Mean-width-extremal measures
   are unique up to scalar on their support (`dim A_t(S) <= 1` for integer `t >= MW_n(S)`, all `n`). Hence no pencil of
   tight sections exists on `Bl_e X` (both forms, every `k`), and "second-order integrality" at Liouville-tight points
   — two tight sections with point-dependent coefficients — is void. This closes the two-section part of the
   equality-case loophole of Theorem B; the single-form part still rests on Conjecture E's structure half.
2. **New item, "integrality at the root arrangement" (dp2 Thm C and transfer, split Lemma 3.5, killer line lemma;
   proved).** Bounded pair cofactors make near-tight tuples integral w.r.t. the diagonals, but on `Bl_e Lbar_k` the
   pair `(sum Zt_ab + E)` has log canonical degree `-1` on lines through `e` for every `k >= 3`, so Lang–Vojta predicts
   nothing at `e`; and Delta-poles amplify nothing beyond sections wherever WTB holds (in particular at the symmetric
   twisted profile, where WTB is Theorem B). Genuine fixed-`S` Subspace on `LM_5 \ ∪ Delta_ab` is not covered.
3. **§3, split calibration (split Thm 3.1 + killer; proved, exact).** Every admissible split exceptional set contains,
   besides `W` and `iota W`, the 15 quadric cones `(u2-u1)(u3-u1) = (u4-u1)(u5-u1)` (mnn31 curves, `alpha = 7/2`) and
   their 15 inversion images (cubic threefolds). Only `Q ∩ sigma Q` survives in the twisted form.
4. **§3, Question B literature form, add (dp2; proved/exact).** A certified dP2 point (`Pic = ZK`, `tau = 3/2 < 2`,
   real-cubic `tau`-curve with `alpha = 2`, no `Q`-curve of degree `<= 6` with `alpha < 2`) is the cleanest test of
   "any `alpha^ess > tau` with `tau < dim`"; its two-point section cap is exactly computed and never beats 3/2 at
   `gcd(a,b) = 1`. Caveat: it has no balance, its top-of-window problem is a Lang–Vojta problem for a log general type
   pair (unlike `X`), and the Losev–Manin tower has no lower-dimensional balanced analogue
   (`alpha_bal(k) = (2^{k-1}-1)/2^{k-3} < k - 1` iff `k >= 5`).
5. **§5 item 3, route 5 bullet, add (interp Thm 3.1; proved).** A linear-algebra sufficient condition for route 5's
   target: if the degree-`e` least space of the slab box is the sum of those of its shift-intersections
   `A ∩ (A - l)`, `l in L_pair`, for all `e < m`, then `Bs V(A,m)` lies in the braid-and-coordinate arrangement
   (round 11's ZE-braid). Evidence on slab boxes is mixed (non-zonotopal success on M8'; below the zonotope order on
   20 of 27 random boxes).
6. **Pointer (split killer; literature via search results only).** `Bl_e LM_5` is a contraction of `Mbar_{0,7}`
   (Castravet–Tevelev, Duke 2015), so `tau_e(-K)`, Theorem B, WTB and Conjecture E are statements about
   `Eff(Mbar_{0,7})`.

## Files

* This report: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/README.md`
* Write-ups: `dp2.md`, `stab.md`, `split.md`, `interp.md` (same directory), with `*_checks/` and `*_notes.md`.
* Kill reports: `kill_dp2.md`, `kill_stab.md`, `kill_split.md`, `kill_interp.md`, with `kill_*_checks/`.
* Judge checks: `judge_checks/j1_rigidity.py` (+ `.out`), `judge_checks/j2_propK.py` (+ `.out`),
  `judge_checks/j3_balanced_liouville.py` (+ `.out`).
