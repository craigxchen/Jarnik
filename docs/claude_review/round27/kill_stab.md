# Adversarial review of round27/stab.md (lane stab)

Status: complete (round 27, adversarial killer for lane stab). Scripts and outputs:
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/kill_stab_checks/`
(written from scratch; none of the author's scripts is imported). Notes: `kill_stab_checks/notes.md`.

## 0. Verdict summary

**The mathematics survives; the novelty claim for the barrier does not.**

* **Correctness.** I found no mathematical error that changes a conclusion. I re-derived each of the following
  by hand and recomputed it with my own exact code:
  * Theorem 3.3 (uniqueness of mean-width-extremal measures, all `n`);
  * Corollary 3.4;
  * Theorem 3.2 (pencil principle);
  * Lemma 3.1 (cofactor identity);
  * Prop 3.5 (chain deficit);
  * Props 1.1–1.3 (Cox form, triangle relations, `Stab ⟺ (L_5^∘)`);
  * the T2 and T3 data.

  The only errors are numerical over-statements and a sign slip, all harmless (§7, issues m1–m4).
* **Novelty.**
  * *Theorem 3.3 is a genuinely new proof in the record.* The record had only evidence for it:
    `threshold/referee_construct.md` §4 verified exhaustively that 2097 tight kernels are one-dimensional.
  * The write-up's "**new barrier**" is not new. "Second-order integrality" means a pencil of tight sections. That
    is the **equality case of the record's mean-width/section cap (Theorem B)**, already listed as an open
    constant-level loophole in `threshold/construct.md` §3 item 3 ("Constant-level effects in the equality case of
    Theorem B. These are closed under Conjecture E or C_bal").
  * Theorem B is itself the Cayley transform of the multiplicity Schwartz–Zippel lemma (DKSS), according to
    `threshold/referee_construct.md` §2. stab.md does not mention this. So Theorem 3.3 is an equality-case uniqueness
    statement for a known inequality. My searches found no such statement in the literature, so its literature
    status is **unclear**.
  * The pencil principle uses only single-section Liouville data, namely the heights `h_L(y)` of divisors
    `L = π^*D − tE` on `Bl_e X`. It is Northcott for the Iitaka map of `L`, and it can only act on the face
    `L·γ̂ = 0` of the section cap. The T4 claim that it lies "outside Lemma K1 and the section cap" is therefore
    overstated: it is the section cap pushed to its boundary.
* **Prop 3.5 (restricted theorem) is correct, but its novelty is small.** The record already flags that width is not
  additive along chains and that the caps are unproved on non-simple profiles (`round25/intra.md` lines 523–524).
  Prop 3.5 turns that flag into an explicit gain. Its hypotheses are satisfiable by balanced profiles, and it touches
  none of the fakes (T1), the k = 4 family (T2), the split construction (T3) or any calibration family (T5).
* **Props 3.6–3.8 are correct reductions that gain nothing.** Prop 3.8's "decoupling" holds because each singleton
  variable occurs only in its own `Π_a`.
* **Tests T1–T5.** I reran all five independently, and all pass. T2 and T3 pass vacuously for the pencil principle:
  Theorem 3.3 holds at every `k` and for both forms, so tight pencils never exist.
* **Overall classification:** `obstruction` stands. Novelty verdict: **repackaged known barrier** (the
  equality/constant-level case of the Theorem B section cap = multiplicity Schwartz–Zippel). It is closed in its
  pencil part by a new and correct lemma (Theorem 3.3).

## 1. What was checked and how

| item | method | script → output | result |
|---|---|---|---|
| Thm 3.3 proof | line-by-line re-derivation (notes.md) | — | correct, incl. edge cases `t−i+1 ≤ 0`, `j = i+1` |
| Thm 3.3 + Thm B | mod-p ranks (certify `dim_Q ≤ dim_p`), exact `fmpz_mat` rank when `dim_p ≥ 2` | `thm33.py` → `thm33.out` | 82 703 sets, 1 638 tight, **0 failures** |
| Cor 3.4 on `P_D ∩ M` of the 4-fold | 30-ray permutohedral fan, `h^0(π^*D − tE) = dim A_t(P_D∩M)` | `cor34.py`, `cor34b.py`, `cor34c.py` | 0 failures; min pencil defect found = 4 |
| Props 1.1–1.2 | own Gaussian arithmetic, **opposite** prime normalisation `π = a+bi, b>a>0` | `cox.py` → `cox.out` | 400 profiles, all identities exact; (Ch) ⟺ primitive & canonical 3000/3000 |
| Lemma 3.1 | exact `Z[i]` division + 60-digit `|f(y)|` | `lemma31.py` → `lemma31.out` | 300 non-simple cases, 0 failures |
| T2 (k = 4) | gauss Thm 4.1 instance, 38 values of `s` incl. random up to `10^9` | `k4_check.py` → `k4_check.out` | residues bounded, chain deficit O(1), nothing fires |
| T3 (split) | all **four** Pell orbits of `z²−56x²=65`, 25 points `x ≤ 4·10^10` | `split_check.py` → `split_check.out` | author's orbit reproduced exactly; deficit O(1) |
| literature | WebSearch (3 queries) + record grep | notes.md | §8 |

## 2. Cox form verification (stab.md §1)

**Re-derived by hand:**
* Prop 1.1(i)–(ii): layer regrouping, and primitivity from (Ch) by choosing `a` in the smallest `π̄`-set.
* Prop 1.2(i): conjugate-primitivity of `Q_ab`, via the four cases (1)&(3), (2)&(4), (1)&(4), (2)&(3).
* Prop 1.2(ii): comparison of cut types by `A ∩ {a,b,c}`. The 4 cuts separating `b` from `{a,c}` contribute
  `γ_Aγ̄_A = n_A`.
* Prop 1.2(iii): `Im(Π_aΠ̄_b) = x_bt_a − x_at_b`. The discriminant of the quadratic form is
  `4y² − 4M_{a∖b}M_{b∖a} = −4t_ab²`, using `y² + t² = N(Q_ab)`.
* Prop 1.3(i)–(iv): `G_ab` is a product over the 4 cuts `A ⊇ {a,b}`, so `|t_ab| ≤ 2e^{ηw}X^{8ε}`. The equivalence
  `Stab ⟺ (L_5^∘)` is immediate from `|t| ≥ 1`. I agree that it is an equivalence with `(L_5^∘)` and only an
  implication from `A_5`.

**Exact check** (`cox.out`). I used 400 random genuine one-class profiles (2–7 primes ≤ 400, exponents ≤ 3, random
unit), with the Gaussian prime normalised **oppositely** to the author's. All of the following held exactly:
* the parametrisation and primitivity;
* `∏ N(γ_A) = N'`;
* conjugate-primitivity of `Q_ab`;
* `ζ_a = G'Q_ab`;
* `N(ζ_a − ζ_b)N(Q_ab) = 4t_ab²N'`;
* the triangle relation for all 60 ordered triples;
* `Π_aΠ̄_b = G_abQ_ab`, the `x/t` identity and the quadratic identity.

181 pairs have `t_ab = 0`; these are coincident points, as in the author's 24 flagged rows. Over 3000 random
`γ`-vectors, (Ch) ⟺ "primitive and canonical" held in 3000/3000. The Prop 1.1(iii) counterexample reproduces: the
input blocks `π̄_5` in `γ_13, γ_3, γ_34` give canonical blocks `π̄_5²` in `γ_3` and `π̄_5` in `γ_134`.

**Assessment.** Both "corrections to the handoff" are right, but mild.
* Correction (i) concerns the converse direction (formula ⇒ canonical blocks), which the handoff never asserted.
* Correction (ii) concerns a triangle relation for the reduced chords that does not appear in the computed task text.

The three "new elementary forms" of §1.5 are correct restatements.
* The arccot form is exact only for the sub-case `t ≡ 1`. In general `cot(θ_a−θ_b) = y_ab/t_ab`, which the write-up
  says.
* **Sign slip (m3).** In §1.5(3) and Prop 3.8, `R | ι+i` gives `x ≡ −it ≡ ιt (mod R)`, not `−ιt`. This is harmless:
  the existence criterion uses the absolute least residue, which is symmetric under `x ↦ −x`.

## 3. Theorem 3.3 (uniqueness of mean-width-extremal measures)

**Proof re-derived, step by step.**
* *Step 1 (dimension peeling).* For `ord μ ≥ t`, multiplying by `(x_1 − c)` lowers the order by 1 and kills the slice
  `x_1 = c`. The kernel is the set of measures on that slice of order `≥ t`. Order in `Z^n` equals order in
  `Z^{n−1}`, because degree-`<t` polynomials restrict onto degree-`<t` polynomials on the hyperplane. This gives
  `dim A_t(S) ≤ Σ_i dim A_{t−i+1}(S_(i))`. Correct.
* *Step 2.* `A_{t−i+1}(S_(i)) ≠ 0` forces `t−i+1 ≤ reg ≤ g_(i) ≤ K−i+1 ≤ t−i+1`, hence `t = K` and `g_(i) = K−i+1`.
  The edge case `t−i+1 ≤ 0` forces `g_(i) = 0`, i.e. a singleton slice, so `dim A_0 = 1`. Correct.
* *Step 3.* (∗) is correct. The `P' = ∅` term has weight `2^{1−n}` and exact value `c_max − c_min`, while (2.5)
  replaces it by `c' − c`.
  * A critical `j` forces `C_v` to be `j` consecutive integers with both ends at `g_(j)` (needs `j ≥ 2`, automatic
    when `i < j`).
  * For two critical indices `i < j`, the two applications of (∗) are strict and sum to `j−1 < j−1`.
  * For `j = i+1` the containment `[b,b+i−1] ⊂ (a, a+i)` is already impossible.

  Correct.

**Exact check** (`thm33.out`; independent code; certified upper bounds via `rank_p ≤ rank_Q`). Test sets:
* all subsets of the boxes `[0,4]×[0,1]`, `{0,1}^3`, `[0,2]×{0,1}^2`, `[0,5]×[0,1]` and `{0,1}^4` (65 535
  subsets);
* random sets in `Z^3` (2000), in `Z^5` (1500) and in `Z^6` (800). The author went only to `n ≤ 4`, so `n = 5, 6` is
  new evidence;
* zonotope unions and perturbations in `n = 3, 4, 5` (2800);
* **cut closures** in `[0,3]^3` and `[0,2]^4` (600). These are the maximal sets with given widths, where a second
  tight measure would have to live, by the construct.md referee's reduction.

Total: 82 703 sets and 1 638 tight sets, with **0** violations of Theorem B and **0** of Theorem 3.3.

**What it is.** Theorem B (`reg ≤ MW`) is, by `threshold/referee_construct.md` §2, the Cayley transform
`w = (i−t)/(i+t)` of the multiplicity Schwartz–Zippel inequality on the grid `{i,−i}^n` (Dvir–Kopparty–Saraf–Sudan;
recorded as inequality (1) of `hermitian_polynomial_barrier.md`). Theorem 3.3 is therefore an **equality-case
statement for multiplicity Schwartz–Zippel on a 2-point grid**: among measures with a given support set, the extremal
one is unique up to scalar. stab.md cites construct.md's Theorem B but not this identification, and it did no
literature search.

My searches found no characterization of equality in multiplicity Schwartz–Zippel (§8), so I rate Theorem 3.3 as
**new to the record, literature status unclear**. In the record, its content had been checked only as evidence:
`referee_construct.md` §4 enumerated every cut-closed set in four boxes and found all 2097 tight kernels
one-dimensional. That referee also observed that `dim K_m(T) ≥ 2` would refute Conjecture E. So "uniqueness half of
Conjecture E" is an accurate description. The structure half (tight ⇒ root-binomial product) remains open.

**A shortcut that does not work.** If some point `s ∈ S` were the unique optimiser of some `±x(P)`, then removing it
would lower `MW`, and `dim ≤ 1` would follow in one line. But vertices can tie in all `2^n` directions (in `{0,1}^2`,
the vertex `(1,0)` ties in every direction). So the slicing argument is genuinely needed. This is a point in favour
of the theorem's non-triviality.

## 4. Corollary 3.4 and the pencil principle (Theorem 3.2)

### 4.1 Correctness

**Corollary 3.4.** The proof is correct:
* With `M = 5`, Cor B1 gives `8·MW_4 = Σ_{15 cuts} width`.
* `width_{P_D}(1_A) ≤ d_A + d_{A^c}`.
* `D·β = Σ_ρ d_ρ`, and `β` is a well-defined numerical class because `Σ_ρ⟨m,v_ρ⟩ = 0`.
* `h^0(π^*D − tE)` is the number of Laurent polynomials on `P_D ∩ M` of order `≥ t` at `e`.
* The `σ_X`-invariant part over `Q` has dimension `≤ 1` as well.

Exact tests:
* (A) 120 root-product classes with `D·β = 8t` have `h^0` **exactly** 1.
* (B′) For sums of 2–3 root-product divisors, `P_D` equals the Minkowski sum (generalized permutohedra), and
  `h^0 = 1`.
* (B″) 8794 random, mostly non-nef classes on the face `D·β = 8t` all have `h^0 = 0`.

So in every sample, the effective part of the face `γ̂^⊥` consists of the root-product classes, i.e. the cone spanned
by the 10 diagonal divisors `Δ̃_ij`. Proving that would be the structure half of Conjecture E, which is still open.

Note: for `D = m(−K)` Cor 3.4 is vacuous. `8t ≥ 30m` forces `t ≥ 15m/4 > ⌊10m/3⌋ = τ`-order, so `h^0` is already 0.
The corollary has content only on the tight face.

**Pencil defect.** Search (C) over 3000 supports in `M` independently gives minimum defect `s = 4`, attained by the
triangle `{0, e_1−e_2, e_1−e_3}`. The defects found were `4, 6, 7, 8, …`. Search (D) over 1736 random classes with
`D·β − 8t ∈ {1,2,3}` found no pencil. So the evidence label "`s_min = 4`" stands.

**Theorem 3.2** is correct as stated:
* `κ_S = 4t + ½Σwidth` dominates `4tε + tη + ½εΣwidth`.
* (a) needs `(f_1,f_2)(y) ≠ (0,0)`.
* (b) counts the points of `P^1(Q(i))` of height `O((η+ε)w)`.
* (c) has `s/30 = (s/2)w / 15w`.

Cor 3.4(iii) is also correct. On a Néron–Severi orbit, `c_2χ_1(λ)F_1 − c_1χ_2(λ)F_2` has a zero whenever
`F_1F_2 ≠ 0` and `χ_1 ≠ χ_2`, so the relation does not descend.

### 4.2 What Theorem 3.2 really is (T4 re-examined), and why the "new barrier" is a known one

At a balanced near-tight point `y`, Lemma 3.1 gives, for every section `f` of `L = π^*D − tE` with support `S`:

`log|c_f(y)| = log|f(y)| + log H_S(y) ≤ h_L(y) + O_f(1)`, with `h_L(y) = h_D(y) − tλ_E(y) = (L·γ̂/2)·w + O((η+ε)w)`.

On simple profiles and on generalized-permutohedral supports, `log H_S` is exactly `h_D` at the finite places. So:

1. **The only arithmetic input of Theorem 3.2 is the size of `c_f(y)` for single sections.** This is the same quantity
   that Liouville with a section bounds below by 1. Theorem 3.2 adds Northcott on `P^1`:
   `h(φ_L(y)) ≤ h_L(y) + O(1)`. This is the textbook companion of "Liouville with a section", applied to the Iitaka
   map of `L`. It is not an inequality, as the author says, but it sees `y` only through the numbers `h_L(y)`, and
   at balanced points those depend only on `L·γ̂`.
2. Consequently it can say something only on the face `L·γ̂ = 0` of the cone that the section cap controls.
   * Theorem B says `L·γ̂ ≥ 0` for every effective `L`: this is the cap `15/4`.
   * The pencil principle needs `κ(L) ≥ 1` on the face `L·γ̂ = 0`.
   * Theorem 3.3 / Cor 3.4 says `κ = 0` there.

   So the pencil route is **the section cap pushed to its boundary**, not a mechanism outside it.
3. **Lemma K1 applies to its inputs.** `c_f(y)` is unchanged by gcd division, since `F(gζ) = g^aḡ^bF(ζ)` and the forced
   divisor scales the same way. So Theorem 3.2 is blind to the gcd reduction, exactly as K1 predicts for the record's
   methods. The T4 sentence "Lemma K1 … do[es] not apply to it" is true only in the narrow sense that Theorem 3.2 is
   not an inequality.
4. **The record already lists this loophole.** `threshold/construct.md` §3 item 3 lists among "inputs not covered by
   the barrier": "Constant-level effects in the equality case of Theorem B. These are closed under Conjecture E or
   C_bal". `referee_construct.md` §6 repeats: "Route (C), constant level … For other tight forms it is conditional on
   Conjecture E".
   * A pair of tight sections with bounded `c`-values is precisely a constant-level effect in the equality case of
     Theorem B.
   * Under Conjecture E it is excluded, because a zonotope determines its root product.
   * The referee's finite reduction ("if some `T` has `dim K_m(T) ≥ 2`, Conjecture E fails") is literally the
     statement that a tight pencil would refute Conjecture E.

So "candidate (a) = pencil principle, void by Theorem 3.3" **closes, unconditionally and only in its two-section
part, a loophole the record had mapped and conditionally closed**. It is a worthwhile clean-up, and Theorem 3.3 is
the new ingredient. But it is not a new barrier, and by the round's own standard it is not new mechanism-level
mathematics.

**Residual content of the record's loophole.** A *single* tight section that is not a root-binomial product would
still give a constant-level inequality not implied by the pair inequalities. Theorem 3.3 does not exclude this, and
stab.md does not claim it does. It remains conditional on the structure half of Conjecture E. My samples (§4.1)
support that half: every effective class on the face is a root-product class.

### 4.3 Is "second-order integrality" exhausted by pencils?

The author's argument that (II) is the only other way to produce a further integer `< 1` is informal ("can be
produced in two ways only"). I checked the natural extensions, and none gives anything new:
* **Different supports.** The ratio `f_1(y)/f_2(y)` carries `D_{S_1}/D_{S_2}`, whose height is of order `w`.
  Translated supports carry `y^{m'}`, of height order `w`. Taking the union support `S_1 ∪ S_2` reduces to Theorem
  3.3.
* **`N`-dependent coefficients.** Theorem 3.3 is coefficient-free.
* **Linear systems of dimension `r ≥ 2` with `s ≥ 1`.** These have `e^{(r+1)(s/2)w}` candidate images, which is no
  finiteness.
* **Pencils with defect `s` on profiles whose chain deficit is about `s·w`.** This is the pencil version of Prop 3.5,
  and stab.md does not mention it. It lives only on the large-chain-deficit class, which Prop 3.5 already handles by
  a single section once `Δ ≥ (s+c)w`. Note also that the minimal pencil (the triangle) has a generalized-permutohedral
  support, hence chain deficit identically 0.

So the informal "only two ways" is right in practice. It is not a theorem; it should be labelled a heuristic
classification.

## 5. Prop 3.5 (chain deficit restricted theorem)

**Correctness.** Lemma 3.1 holds exactly on 300 independent non-simple cases (`lemma31.out`: divisibility,
`N(F)/N(D_S) = |f(y)|²∏p^{width}` to 30 digits, and subadditivity, all with 0 failures). The inequality of Prop 3.5
then follows in one line, and I re-derived it: `2κ_S(η+ε)` dominates `εΣwidth + 8tε + 2tη`.

The example `S = {0, e_1−e_2+e_3−e_4}` has `Σwidth = 12`, `t = 1` and `s = 4`. For `a(p) = (2,1,0,1,0)` (chain
`{1,2,4} ⊃ {1}`) the deficit is `2 log p`; both are confirmed.

**Non-vacuity (my addition).** The two `±2` cuts `{1,3}` and `{2,4}` can never be cancelled along a chain: subsets
of `{1,3}` have `⟨m,·⟩ ≥ 0`, and no superset has a negative value. The eight `±1` cuts pair into the nested
opposite-sign chains `{1}⊂{1,2,4}`, `{3}⊂{2,3,4}`, `{2}⊂{1,2,3}` and `{4}⊂{1,3,4}`. So balanced profiles reach
`Δ_S` up to about `8w > s·w = 4w`. The hypothesis of Prop 3.5 is satisfiable, and the proposition is a genuine (if
narrow) restricted theorem.

**Scope corrections.**
* The headline "on balanced profiles that are not simple, one-section Liouville beats 15/4" overstates. The gain
  needs `Δ_S ≥ (s_S + c)w` with `s_S ≥ 1`. In every example `s_S ≥ 4`, so chain primes must carry at least
  `(s_S + c)/2 · w ≥ 2w` of log-weight. A non-simple profile with `o(w)` chain weight gains nothing. The correct
  wording is "on balanced profiles whose chain deficit for some fixed `S` exceeds `(s_S + c)w`".
* **Novelty is small.** `round25/intra.md` lines 523–524 (open item 3) already record: "Non-simple balanced profiles
  (`e_p ≥ 2`): … for chains `Σ_t 1_{A_t}` the width is not additive and Theorem 5.1 is not proved (as in
  `round24/theory.md`, which also assumes simple profiles)". Prop 3.5 converts this flag into an explicit statement
  that the cap does fail there, with the gain equal to the deficit.

  In the language of heights, the chain deficit is extra non-archimedean proximity of `y` to the base locus of the
  section. Newton-polygon forced divisibility `min_S⟨m,a(p)⟩` exceeds the toric-divisor bound `−ψ_{D_S}(a(p))`,
  which is additive on the cones of the permutohedral fan. The deficit vanishes identically for
  generalized-permutohedral supports, in particular for all root-binomial products and for every pair inequality.

  This is the record's class `𝓛` (round 24: "toric forced divisibility at every place"), which round 24 capped only
  on simple profiles. So Prop 3.5 is new as a stated inequality, but it sits in a class and a gap the record had
  already identified.
* **The heuristic weight is low.** Profiles in which chain primes cover `≥ 2` cut-weights lose configuration freedom
  (one prime set serves two cuts), so the heuristic count is even more negative there. Prop 3.5 proves rigidity where
  rigidity was most expected, and nothing on the simple (squarefree) core. The write-up says this itself ("removes
  no part of the hard core").
* **Calibration** (independent). The k = 4 dense family has maximum chain deficit 7.2 against `s_S·w ≈ 80` at
  `s ~ 10^9` (`k4_check.out`). The split Pell points have maximum deficit **11.15** (not `≤ 6.8`, see m1) against
  `4w ≈ 36`. Neither fires, as it must not.

## 6. Props 3.6–3.8 (candidates b, c)

* **Prop 3.6.** This is Theorem 3.2 for maps to `P^n`. It is correct, and it inherits the `s ≥ 1` of Cor 3.4.
* **Prop 3.7.** The algebra is correct: on `F_c`, `H = f_2^n g_c`, and the forced divisor of `H` is `D_S^n D_{S'}`.
  The projection-from-`e` example is exact: `c_aσ = Π_a/Π̄_a − 1 = 2it_a/Π̄_a`, and `σ = 2i/Γ̄`.

  The record does list direction-sensitive bounds as uncovered: `round25/README.md` line 73, "Direction-sensitive
  bounds (route 5's leading form `Φ_c(u)`) are **not** covered", and `referee_intra.md` line 36. So the placement is
  right. **Citation slip (m4):** stab.md attributes this to "round 25, intra item 3", but that item concerns
  non-simple profiles. The record also notes (`intra.md` lines 279–282) that direction refinements "do not touch the
  non-archimedean side, where the gap is", and Prop 3.7 confirms this in its one computed instance (exactly tight).
  The genus `≥ 1` discussion is a reasonable pointer, not a result. The uniform-Mordell attribution (Dimitrov–Gao–
  Habegger, Kühne) is from memory and labelled so; I did not re-check it.
* **Prop 3.8 ("decoupling").** It is true, but it is a tautology. The singleton `γ_{\{a\}}` occurs only in `Π_a`, so at
  fixed non-singletons the four conditions are on four different variables. The one-point-per-line count
  (`|R_a| ≈ X^{7/2} > 2X^{(1+ε)/2}`) and the heuristic count `X^{11}·X^{−12} = X^{−1}` are correct. They are the
  coordinator's heuristic in residue form. "No transference" adds no information: all transference would have to come
  from the joint variation of the `R_a` through the 11 shared blocks, which is the original problem.

  The verdict "(c) reduces to counting below one point per scale" is right, and it is the record's existing
  assessment ("density mechanism excluded only heuristically, count `e^{−w}` per scale"). Nothing new as a barrier.

## 7. Pre-registered tests T1–T5 (rerun independently) and the issue list

* **T1 (fakes).** No growth claim is made.
  * Theorem 3.3 and Cor 3.4 concern measures and classes, not points.
  * Prop 3.5 is identically void on simple profiles: one layer per prime gives `Δ_S = 0` exactly. Every recorded fake
    is squarefree, hence simple.
  * Theorem 3.2 would act on fakes as well, since its inputs are single-section Liouville quantities and fakes satisfy
    every monomial gap. But it never applies (Cor 3.4).

  **Pass.**
* **T2 (k = 4).** On the gauss Thm 4.1 instance (`k4_check.out`):
  * `Im P_x(s) ≡ −420, −180, −240` (sympy, symbolic);
  * over 38 values of `s` (8 structured, 30 random in `[10^5, 10^9]`): gcd norm `≤ 10368000`, balance deviation
    `→ 0.09`, residues bounded, chain deficit `≤ 7.2 ≪ s_S w`, and Prop 3.5 never fires.

  Theorem 3.3 holds for all `n` (also checked at `n = 3`, inside `thm33.out`), so no tight pencil exists at k = 4 and
  the mechanism proves nothing there. **Pass, but vacuously:** the pencil principle is void at every `k`, so T2 does not
  discriminate. As the author notes, even a hypothetical tight pencil would only give a `T`-dependent `Z_{T,K}`, which
  the dense family does not contradict.

  **Numerical over-statement (m2).** The author's "`|c_ab|² ≤ 14400`" is a sample maximum over `s = 10^k + 7`, where
  11 of 12 values give 162. Random `s` give up to **230400** (`|t| = 240`). The true bound is `4·max|Im P_x|² = 705600`.
  Boundedness, and hence T2, is unaffected.
* **T3 (split).** I used all four Pell orbits of `z² − 56x² = 65` (`split_check.out`; 25 points, `x ≤ 4·10^10`; the
  author used one orbit, 7 points).
  * The identity `F_1F̄_1 = n` holds, and so does the twin representation.
  * The author's orbit is reproduced exactly: `α = 6.5143, 4.6651, 4.3679, 4.1969, 4.0970, 4.0071, 3.9913`, and
    residues `[8,5,9,5], [1,5,1,1], [8,1,1,1], [1,1,9,5]`.
  * On the other orbits `α` dips to `3.948` at `x ≈ 2.9·10^9`, still decreasing towards `15/4`.
  * Residues stay in `{1,2,4,5,8,9,14,28}`.

  Nothing in stab.md proves split Question B, so there is no #887 claim. The construction is not in Prop 3.5's
  exceptional set: across all orbits the maximum chain deficit is **11.15** (`x = 98 382 962`, `w = 9.0`), against
  `s·w = 36`. **Pass.**

  **Numerical over-statement (m1).** "Deficits stay `≤ 6.8`" is true on the author's orbit only. Over all four orbits
  the maximum is 11.15. It is still O(1), so the conclusion is unchanged.
* **T4 (criticality).** Partially fails as a *novelty* claim (§4.2).
  * Prop 3.5 is correctly placed: class `𝓛`, on profiles where the class-`𝓛` cap is unproved.
  * Prop 3.7 is correctly placed in the direction-sensitive class, which is uncapped. The citation is wrong (m4).
  * Theorem 3.2 is placed "outside K1, the section cap, class `𝓛` and the multipoint reduction". That is overstated. Its
    only inputs are single-section quantities `c_f(y)`, which are gcd-invariant (K1) and controlled by `L·γ̂`. Its only
    possible domain is the boundary face of the section cap. Its failure is the equality case of Theorem B, which the
    record had already mapped (`construct.md` §3 item 3). It is not a new barrier.
* **T5 (calibration).** None of the following is excluded:
  * the round-24 real family `M`: no `Q`-points; Theorem 3.2 never applies; Prop 3.5 is void on simple profiles;
  * the split Pell curves: checked above;
  * the k = 4 family: checked above;
  * Huang's `Y_4`, `P^2`, and the 8- and 10-point families. stab.md's results are a measure-theoretic statement plus a
    proved Liouville inequality, which cannot exclude existing points, and there is no growth claim.

  **Pass.**
* **Handoff constraints.** I checked all four, and none is violated:
  * no fixed-`S` Subspace;
  * the Gaussian gcd is kept exactly (Lemma 3.1, exact check);
  * no drifting phase is set to zero;
  * template-dependent thresholds (`Z_{T,K}`, Prop 3.5's family) are stated as such.

**Issue list.**

| id | severity | issue |
|---|---|---|
| M1 | major (novelty) | The "new barrier" (pencil principle void) is the equality / constant-level case of Theorem B, which the record already lists as an open loophole closed conditionally under Conjecture E (`construct.md` §3 item 3; `referee_construct.md` §§4, 6). Theorem 3.3 closes its two-section part unconditionally. It is a repackaged known barrier with a new lemma, not new mechanism-level mathematics. |
| M2 | major (T4 claim) | The claim that Theorem 3.2 is outside Lemma K1 and the section cap is overstated. Its inputs are single-section Liouville quantities (heights `h_L`, `L = π^*D − tE`), gcd-invariant and determined by `L·γ̂` at balanced points. It is Northcott for the Iitaka map on the boundary face of the cap. |
| m0 | minor (novelty bookkeeping) | Theorem B is the Cayley transform of multiplicity Schwartz–Zippel (DKSS), per `referee_construct.md` §2, which stab.md omits. So Theorem 3.3 is an equality-case uniqueness for a known inequality. Its literature status is unclear (no hit in my searches); the author did no search. The record already had exhaustive evidence (2097 tight kernels, all of dimension 1). |
| m1 | minor | "Chain deficits `≤ 6.8`" holds on one Pell orbit only; over all four orbits the maximum is 11.15. Still O(1). |
| m2 | minor | "`|c_ab|² ≤ 14400`" (k = 4) is a sample maximum over structured `s`; random `s` give 230400; true bound 705600. Still bounded. |
| m3 | minor | Sign slip `x ≡ −ιt` vs `ιt` in §1.5(3) and Prop 3.8. Harmless. |
| m4 | minor | The direction-sensitive class is cited as "intra item 3". It should be `round25/README.md` line 73 / `referee_intra.md` line 36. Intra item 3 is non-simple profiles. |
| m5 | minor | The Prop 3.5 headline "beats 15/4 on non-simple profiles" needs "with chain deficit `≥ (s_S+c)w`", which requires chain primes carrying `≥ 2w`. The non-additivity is already flagged in `round25/intra.md` 523–524. |
| m6 | minor | "Second-order integrality can be produced in two ways only" is an informal classification. Label it heuristic; §4.3 checks the obvious extensions. |
| m7 | minor | Cor 3.4 is vacuous for multiples of `−K` (`15/4 > 10/3`). Its content is on the face `γ̂^⊥`, whose effective part in all my samples is the root-product cone (Conjecture E's open structure half). |

## 8. Literature check

* **Record (grep).** As cited above:
  * `threshold/construct.md` §2.4 Remarks and §3 item 3;
  * `threshold/referee_construct.md` §2 (Theorem B = Cayley transform of multiplicity Schwartz–Zippel), §4 (finite
    reduction; 2097 tight kernels of dimension 1), §6;
  * `round25/intra.md` 279–282, 523–524;
  * `round25/README.md` 73;
  * `round24/theory.md` (simple profiles only).
* **WebSearch** (snippets only; labelled as such).
  * Multiplicity Schwartz–Zippel: Dvir–Kopparty–Saraf–Sudan, with an algorithmic version by Bhandari–Harsha et al.
    (arXiv 2111.11072). Related hypercube work: Sauermann–Wigderson (arXiv 2010.00077), Ghosh et al. (arXiv
    2307.16881), and weighted multiplicities (arXiv 1707.01354). **No characterization of equality cases found.**
  * Regularity index and lattice ideals (arXiv 1301.4249, 1206.1892): nothing on mean width or uniqueness.
  * Castravet–Tevelev, "`M̄_{0,n}` is not a Mori dream space" (arXiv 1311.7673), and Castravet–Laface–Tevelev–Ugaglia
    (arXiv 2009.14298): their framework is exactly `Bl_e` of toric varieties with sections given by Laurent
    polynomials vanishing at `e`. The snippet states rational contractions `Bl_e LM_{n+1} ⇢ M̄_{0,n} ⇢ Bl_e LM_n`, so
    Cor 3.4 is a statement about a face of the effective cone of a birational model of `M̄_{0,7}`. I did **not** check
    whether that face is described there. This is an open pointer for a referee, not a priority claim.
* **Not checked:** full texts of the above; `hermitian_polynomial_barrier.md` (outside `docs/claude_review/`).

## 9. Novelty decision

**Overall: repackaged known barrier**, with one genuinely new lemma.

* **Barrier (candidate (a)).** It is the equality-case / constant-level loophole of the record's section cap
  (Theorem B, mean width; equivalently multiplicity Schwartz–Zippel through the Cayley transform). That loophole is
  listed in `threshold/construct.md` §3 item 3 and was conditionally closed under Conjecture E. The pencil principle
  is Northcott for the Iitaka map of a divisor on the boundary face `γ̂^⊥` of that cap.
* **Theorem 3.3.** It is new to the record as a proof: the record previously had exhaustive evidence only. Its
  literature status is **unclear**. It is an equality-case uniqueness statement for multiplicity Schwartz–Zippel on a
  2-point grid, and I found no prior statement of it. It is correct for all `n`.
* **Prop 3.5.** Correct. It makes explicit a gap already flagged in `round25/intra.md` (non-additivity on chains).
  Its novelty is small, and it gives no gain on the simple core.
* **Props 3.6–3.8 and §1.5.** Correct reformulations and reductions. Prop 3.8 is a tautology at fixed non-singletons.

Measured against the round's stated goal (new mathematics, without rediscovering a mapped barrier), the lane:
* produced one correct new lemma (Theorem 3.3);
* used it to close a mapped loophole;
* did **not** produce a new mechanism or a new barrier.

## 10. Corrected statement

*Status: obstruction. No proof of Question B and no growth improvement; `M ≤ (2/3+ε) log R/log log R` stands.*

1. **Theorem 3.3** (proved; independently re-derived and checked on 82 703 sets in dimensions 1–6, with 0 failures).
   * For finite `S ⊂ Z^n` and integers `t ≥ MW_n(S)`: `dim A_t(S) ≤ 1`.
   * Equivalently, via the Cayley transform of `referee_construct.md` §2: an extremal measure for the
     multiplicity-Schwartz–Zippel / mean-width bound is unique up to scalar on its support.
   * This is the uniqueness half of `construct.md`'s Conjecture E, which the record had only verified exhaustively in
     small boxes. The structure half stays open.
   * New to the record; literature status unclear.
2. **Corollary 3.4** (proved).
   * Statement: on `Bl_e` of either form of the Losev–Manin 4-fold, `h^0(π^*D − tE) ≤ 1` when `8t ≥ D·β`. Every
     effective divisor on the face `γ̂^⊥` therefore has Iitaka dimension 0.
   * Scope: it is vacuous for multiples of `−K`. In samples, the effective part of the face is the cone of the 10
     diagonal divisors.
   * Consequence: there is no pencil of tight sections. Every pencil has defect `s ≥ 1`; `s_min = 4` is evidence.
3. **What this closes** (this replaces "new barrier"). It closes, unconditionally, the *two-section* part of the
   record's mapped loophole "constant-level effects in the equality case of Theorem B" (`construct.md` §3 item 3).
   * Bounded-residue stability can never come from a pencil of tight sections. Such a pencil, through Northcott for its
     map to `P^1`, is the only form "second-order integrality" takes in the practical classification of §4.3, which
     is heuristic, not proved.
   * The pencil principle is not a mechanism outside the section cap. Its inputs are single-section Liouville
     quantities. They are gcd-invariant, consistent with Lemma K1, and determined by `L·γ̂` at balanced points.
   * The single-tight-form part of the loophole remains conditional on Conjecture E's structure half.
4. **Prop 3.5** (proved). It is a restricted statement:
   * On balanced near-tight one-class points whose chain deficit for some `S` in a fixed finite list exceeds
     `(s_S + c)w`, one-section Liouville forces `f_S(y) = 0`.
   * The hypothesis is satisfiable (deficit up to about `8w` versus `s = 4` for `S = {0, e_1−e_2+e_3−e_4}`). It
     requires chain primes carrying at least about `2w` of log-weight.
   * It is void on simple (squarefree) profiles. It makes explicit the non-additivity already flagged in
     `round25/intra.md` lines 523–524.
5. **Elementary forms** (proved).
   * The Cox parametrisation is a bijection exactly on (Ch).
   * The triangle relation `Q_abQ_bc = G^{(b)}_{abc}Q_ac` holds for pair cofactors.
   * `x_bt_a − x_at_b = G_abt_ab`, and the quadratic identity holds.
   * `Stab ⟺ (L_5^∘)`, while `A_5(1/30−η)` only implies `Stab`.
   * The arccot form is exact for `t ≡ 1`.
   * Sign fix: `x_a ≡ ι(R_a)t_a (mod N(R_a))` when `R_a | ι + i`.
6. **Reductions** (proved, no gain).
   * Fibrewise section arguments are global direction-sensitive section arguments. This class is uncovered in the
     record (`round25/README.md` line 73), and the one computed instance is exactly tight.
   * At fixed non-singleton blocks the four singleton conditions are trivially independent. Candidate (c) is
     counting below one point per scale.
7. **Tests.**
   * T1, T3 and T5 pass.
   * T2 passes vacuously; the pencil route is void at every `k`.
   * T4 passes for Prop 3.5 and Prop 3.7 and is overstated for Theorem 3.2 (M2).
   * Numerical corrections: the split chain deficit reaches 11.15 over all Pell orbits, not "`≤ 6.8`". The k = 4
     residues `|c_ab|²` reach 230400 at random `s`, not "`≤ 14400`"; the bound is `≤ 705600`.
