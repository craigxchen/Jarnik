# Referee report, block "reconcile": Part V (the fixed-field F_2 cluster) against the rest of the record

Scripts and outputs: `round31/ref5_reconcile/`. All code is mine.
* `f24.py`: F_2[t] arithmetic (integers as bit vectors, windowed carry-less multiplication) and F_4[t] arithmetic (pairs A + Bθ, θ² = θ + 1). Unit test: `test_f24.py`.
* `cluster.py`: builds the Part V cluster. It checks the norms and the pair differences, computes the small-prime valuations and the gcd degrees (spot-checked against Euclidean gcds), computes the cut slack Q in two independent ways, and evaluates both sides of (IV.2) and the transferred cut identity. Output: `cluster_m*_L*.json`.
* `dict_check.py`: the explicit star dictionary in characteristic 2 (U, P_i, columns, n_i, D_ij), with every hypothesis of (IV.2) asserted. It evaluates (IV.2) raw and after normalisation. Output: `dict_m*_L*.json`, `dict_m64.out`.
* `fstar_f2.py`: the transfer of `two_thirds/proof.md` Thm 5.2 to F_2[t], computed exactly. Output: `fstar_f2.out`.
* `five.py`: all 5-point subconfigurations (m = 16), or random samples (30000 at m = 32, 8000 at m = 64). It reports the reduced exponent δ_S, α = 4/(1+2δ_S), the 15 cut weights, H_Σ, and an exact check of the Liouville identity.
* `separable.py`: tests separability of the classifying point (cross-ratios and ratios).

## 0. Verdict

* **(a)** The Part V cluster gives exact data for (IV.2) in characteristic 2, with δ = 0 and E_ij = 0, and it violates (IV.2) for m = 32 and m = 64. The only step of the proof of (IV.2) that fails is Langer's Theorem 0.1. Every other step works in any characteristic: the surface, the blow-ups, effectivity, the lc/SNC conditions and the ledger. The resulting rational surface over F̄_2 has 3e_orb − (K+B)² ≤ −13781 (m = 32) and ≤ −305692 (m = 64). The classifying point is **separable**, so this is not the Frobenius failure.
* **(b)** The record's 3M log M bound transfers to F_2[t] exactly, as W ≥ 3M log₂M − O(M). Part V satisfies it, with a factor of about (2/3)m to spare.
* **(c)** The Part V cluster refutes the F_2(t)-analogues of A_5(δ) for every δ < 1/30, of (L_5°), and asymptotically of (L_8). Its independent 5-subtuples are almost exactly balanced, with exponent α → 15/4 from above and H_Σ/w ≤ 40/m. So in characteristic 2 the balanced Liouville exponent 15/4 is attained and cannot be confined to any proper exceptional set. In this model the answer to the analogue of Question B is **no**.
* **(d)** Part V is not the inseparable/Frobenius failure of Vojta in char p. It is the additive (Artin–Schreier) mechanism: Galois acts on Y = t^L + θ by translation by 1, and subspace polynomials are additive. It belongs to the same family as the known characteristic-p counterexamples to BMY that come from configurations existing only in characteristic p.
* **(e)** F3+ holds for uniform statements and nothing weaker. It kills every recorded unconditional route to the uniform bound, and every char-free reduction step (J4, (L_8) ⇒ uniform, Vojta ⇒ A_5). It does not touch any growth-rate argument.

## (a) From the F_2 cluster to star data, and the violation of (IV.2)

**Construction re-checked independently.** I rebuilt z_H = ∏_{v∈H∖0}(Y+v) with Y = t^L + θ and L = (d−2)2^{d−1} for m = 4, 8, 16, 32, 64. For every mark, z_H·z̄_H = N ∈ F_2[t]. The points are distinct, all of degree D = (m−1)L, and the maximum pair-difference degree is exactly D/2: 6, 56, 360, 1984, 10080. This extends Astra's table by m = 64 (2016 pairs). With L = 129 (m = 32, D = 3999 odd) the maximum is 1999 = ⌊D/2⌋, so odd D with δ = −1/2 also works.

**Dictionary in characteristic 2 (CORRECT; the proof in ref4_bridge §1 needs no characteristic).**
* No self-conjugate prime divides any Y+v. A common root would give θ ∈ F_2(α), so every prime of N splits. This is the analogue of "no real roots".
* Y+v and Y+v+1 = conj(Y+v) are coprime.
* An F_4-prime f can divide Y+v and Y+v′ only if P = f f̄ divides v+v′, which has degree < d. So the "big parts" of distinct Y+v are pairwise coprime, and conj(big(v)) = big(v+1).
* Orientation: U = ∏_{v even} big(v) · ∏ f^{v_f(N)}, with one f chosen from each small conjugate pair. Then P_H = gcd(z_H, U) is explicit.
* In characteristic 2, θ̄ − θ = 1, so P_iP̄_j − P̄_iP_j = D_ij = X_iY_j + X_jY_i. With P_i = GA and P_j = GB this gives D_ij = GḠ(AB̄ + ĀB), so g_ij = GḠ divides D_ij and E_ij = 0.
* Computed:
  * m = 8, 16, 32, all pairs (28 / 120 / 496): all hypotheses hold. That is, UŪ = N, gcd(U,Ū) = 1, z_i P̄_i = P_iŪ, primitive columns, n_i | N, D_ij ≠ 0, g_ij | D_ij, the exact identity deg(z_i−z_j) = D + deg D_ij − (α_i+α_j)/2, and deg D_ij ≤ (α_i+α_j)/2 − W/4 (the maximum excess is exactly 0).
  * m = 64, 300 random pairs: the same results. gcd(U,Ū) = 1 was not recomputed there; it holds by the structural argument above.

**Both sides of (IV.2).** L_IV = W = deg N. Q = Σ_h (s_h − m/2)² is computed (i) from the layers and (ii) as m²W/4 − Σ_{i<j}(W − 2 deg gcd(z_i,z_j)). The two agree.

| m | W | Q | δ | LHS = ¾(m−24)W + Q | RHS | 3e_orb − (K+B)² ≤ −(2/9)(LHS−RHS) |
|---|---|---|---|---|---|---|
| 32 (raw) | 7936 | 14400 | 0 | **62016** | 0 | −13781.3 |
| 32 (normalised, global gcd of degree 8 removed) | 7920 | 10304 | −4 | 57824 | −4480 | — |
| 32, L = 129 | 7998 | 15264 | −1/2 | 63252 | −560 | −14180.4 |
| 64 (raw) | 40320 | 166016 | 0 | **1375616** | 0 | −305692.4 |
| 64 (normalised, gcd of degree 30 removed) | 40260 | 104576 | −15 | 1312376 | −64320 | — |

For m = 16 the LHS is −7408 and no violation occurs, as expected for m ≤ 24. Note that the big-part layers all have s_h = m/2. Only primes of degree < d unbalance the layers, and at those primes there are layers with s_h ∈ {0, m}: the cluster has a small common gcd. Removing it only makes the violation larger.

**Which step fails.** The surface is built as in Part IV: sections of P¹×P¹ → P¹ over F̄_2, elementary transformations at ∞, blow-ups of the collision trees, weights 2/3, 1/2 and 0. This is a blown-up ruled (rational) surface over F̄_2. In any characteristic the following hold:
* the contact orders are valuations of minors;
* SNC, klt, and Q-effectivity of K+B (by adjunction and the discrepancy count);
* the ledger, which is pure intersection theory together with ℓ-adic Euler numbers.

The ledger gives 3e_orb − (K+B)² ≤ −(2/9)(LHS − RHS), which is **negative** (table). So the one non-characteristic-free input, Langer's orbifold BMY (Thm 0.1), is what fails in characteristic 2. Its proof in characteristic 0 rests on Bogomolov–Miyaoka semipositivity of the orbifold cotangent sheaf (the Kodaira–Spencer input); that property is false in characteristic p. BMY is known to fail in positive characteristic (Lang 1983; Easton, PAMS 136 (2008), using covers of P² branched over line configurations that exist only in characteristic p). Here the counterexample is a **rational surface with an SNC pair of standard-type weights**. It comes from a genuine anisotropic circle cluster over a **fixed** constant field, and its classifying map is separable (see (d)). I found no such example in the literature (search not exhaustive).

**Verdict (a): CORRECT.** The F_2 cluster is an instance of the hypotheses of (IV.2) with δ = 0, E = 0, m = 32 and m = 64 (> 24). Characteristic-2 failure of (IV.2) was already known from the trace seed (ref4_arith), but only over F_{2^d} with d growing and with the split form. The new content is the fixed field together with the anisotropic genuine norm.

## (b) The transferred growth bound

Here is the transfer of `two_thirds/proof.md` §§1–5, checked step by step.
* Primes of F_2[t] of odd degree e are inert. Their norm fibre mod P has 2^e + 1 elements.
* Primes of even degree e split. Their fibre has 2^e − 1 elements.
* Hensel lifting works in characteristic 2. The polar form of A² + AB + B² is (y, x), which is nonzero mod P.
* There is no ramified prime, so the term K log2/2 is absent.
* log p becomes deg P, and log R becomes D.
* The cut identity, multiplicativity of collisions, the chain lemma and normalisation are unchanged.

> **Theorem (transfer; proof as in `two_thirds`).** For M distinct z_j ∈ F_4[t] with z_j z̄_j = N and deg(z_i − z_j) ≤ D/2 + c:
> Σ_{i<j} deg c_ij + Q/2 = Σ deg(z_i−z_j) − M(M−2)D/4 ≤ MD/4 + C(M,2)c.
> Hence F₂*(M) ≤ MW/8 + C(M,2)c, i.e. **W ≥ 8F₂*(M)/M − 4(M−1)c = 3M log₂M − O(M) − 4(M−1)c**.
> Here F₂*(M) = Σ_{e odd} e·I(e)·Σ_a E(M,(2^e+1)2^{e(a−1)}) + Σ_{e even} e·I(e)·G_e(M), with I(e) the number of irreducibles of degree e and G_e as in Thm 5.2 with p replaced by 2^e. Equivalently M ≤ (1/3 + o(1)) W/log₂W.

The leading constant is as for the integers. The inert primes give (M²/2)·(log₂M)/2 and the split primes give (M²/4)·(log₂M)/2.

**Exact values against Part V** (`fstar_f2.out`):

| m | W_min = 8F₂*/m | 3m log₂m | Part V deg N = 2m(m−1)(log₂m − 1) | ratio |
|---|---|---|---|---|
| 32 | 223 | 480 | 7936 | 35.6 |
| 64 | 674 | 1152 | 40320 | 59.8 |
| 128 | 1611 | 2688 | 195072 | 121 |
| 256 | 4210 | 6144 | 913920 | 217 |

The construction is **consistent**. I also verified the exact cut identity on the cluster: m = 32 gives 30960 = 30960 with budget mD/4 = 31744, and m = 64 gives 318432 = 318432 with budget 322560. The construction is a factor of about (2/3)m above the transferred lower bound: it has m ≈ √(W/(2 log₂W)) points, against the permitted ≈ W/(3 log₂W). Part V says nothing about the growth constant over F_2.

**Verdict (b): CORRECT** (the transfer is routine). It is new only in its fixed-field form; ref4_arith did the F_q version with q growing.

## (c) Questions A and B in the F_2 model

**Five-point subtuples, exact formula.** Take five marks with linearly independent forms. Then:
* |∩H| = m/16, so g_S ≥ (m/16 − 1)L;
* the reduced degree is D′ ≤ (15m/16)L;
* each pair satisfies deg(z_a − z_b) = (m/2 − 1)L + p_ab with p_ab = deg P_J(1) ≤ L/2.

Hence δ_S = (L/2 − p + g_S/2)/(D − g_S), which gives

  **1/30 − 8/(15m) ≤ δ_S ≤ 1/30**, i.e. α_S = 4/(1+2δ_S) → 15/4 from above.

The upper bound is char-free Liouville at near-exact balance. Each of the 15 cuts receives m/16 conjugate pairs {v, v+1}, with weight ≈ 2L each. So w ≈ mL/8, the balance ratio tends to 1, and H_Σ = Σ deg c_ab ≤ 5L, giving H_Σ/w ≤ 40/m. All points are monic, so they lie in one unit class.

Computed (`five.py`, rank-5 tuples; the exact identity deg(ζ_a−ζ_b) = D″ − d_ab/2 + deg c_ab held on every tuple checked):

| m | δ_S range | α range | max balance ratio | max H_Σ/w̄ |
|---|---|---|---|---|
| 16 | [0, 0.0106] | [3.92, 4.00] | 1.25 | 2.31 |
| 32 | [0.0176, 0.0221] | [3.83, 3.86] | 1.081 | 1.11 |
| 64 | [0.0261, 0.0280] | [3.79, 3.80] | 1.036 | 0.53 |

At m = 64 every independent 5-subtuple is an A-point for δ = 1/40. Dependent (rank-4) tuples have 7 effective cuts and δ_S → 1/14, which is the Liouville value for that profile and above 1/30. Those may lie in a special subvariety; the argument does not need them.

**Escape from every exceptional set.** I ran J4's argument (round 23) over F_4(t). Steps 1–7 hold in any characteristic: the factorisation, the averaging, the A-point inequality (which at c = 0 needs no condition on R), and D″ ≥ D/2. Schwartz–Zippel works over any field. So for every proper closed Z ⊂ T_5, including Z defined over F_2(t), a homogeneous F of degree D_Z has at most D_Z m⁴ zeros on Σ⁵. The number of ordered independent 5-tuples is m(m−1)(m−2)(m−4)(m−8), which is larger once m ≫ D_Z.

**Hence (PROVED, given the transfer above):**
* The F_2(t)/F_4(t) analogue of **A_5(δ) is false for every δ < 1/30**. It suffices to take m ≥ 8/(15(1/30 − δ)).
* So is the analogue of **(L_5°)** for every c > 0 and η_0 > 0.
* So is the analogue of **V_E(5)**, via gauss Thm 4.2, whose deduction is algebraic.
* In McKinnon–Roth terms, the essential exponent at e along balanced sequences is exactly the Liouville value 15/4. The open window (15/4, 4) collapses to its lower end in this model.

**Question A.** Take eight independent forms (m ≥ 128). All 127 cuts get equal weight w = mL/64, up to small primes. The reduced chords have degree ≤ L/2, so H_Σ ≤ 14L and H_Σ/w ≤ 896/m. The same Schwartz–Zippel escape works on (P¹)⁸ → M̄_{0,8}. So the FF-(L_8) and (L_8^c)^all analogues fail once m ≫ 896/c.

Caveat: I identify the (L_8) quantity H_Σ with reduced-chord degrees through the star/Plücker dictionary (ptolemy Thm B). I did not re-derive that dictionary in characteristic 2, and I did not compute m = 128 (D ≈ 9.75·10⁴). Status: **CORRECT_WITH_CORRECTIONS**. The statement holds asymptotically, not at m = 128, where H_Σ/w can be as large as 7.

So Part V is a genuine char-2 cluster that violates both Question A and Question B in the function-field model. The record's remark that (L_8) is "fake-blind" is irrelevant here, because this is not a fake.

## (d) Relation to the known failures of Vojta in characteristic p

* **Not Frobenius / inseparable.** All cross-ratios of 4 points and all ratios z_a/z_b tested have nonzero derivative (60/60 at m = 8, 16 and 32, and at L = 128 and L = 129). The classifying map P¹ → M_{0,m} or T_5 is therefore separable, and it is not a Frobenius twist. The record already noted that (IV.2) is invariant under Frobenius base change, so pure Frobenius could never violate it.
* **The mechanism is additive (Artin–Schreier).** F_4(t)/F_2(t) is the Artin–Schreier extension θ² + θ = 1. Its Galois action is translation by 1 on Y = t^L + θ. The subspace polynomial P_J is F_2-linear, which gives P_J(Y + a) = P_J(Y) + P_J(a) and makes every pair difference factor completely: P_J(1)·P_J(Y)/Y. This is the same additive structure as in the trace seed. It is also behind BMY/Hirzebruch failures for configurations that exist only in characteristic p, and behind wild ramification. It is not the "p-th power" failure in which one point is replaced by its Frobenius iterate.
* **Precise status.** The record's conditional routes are Vojta for number fields, so they are not contradicted. Part V shows that the **function-field (char 2, fixed constant field, separable) analogues** are false:
  * Vojta D = E on the toric 4-fold, via A_5;
  * D = 0 on M̄_{0,8}, via (L_8);
  * D = 0 on Bl_Y X_7, via the uniform bound, assuming the record's deduction transfers (not checked).

  Characteristic-p Vojta is not expected to hold without extra hypotheses, so this matches known expectations. But Part V is a separable, non-isotrivial instance. It therefore rules out the most common repair, "assume the point separable/non-isotrivial", for these varieties. Whether a published characteristic-p formulation (for example one with Hasse–Wronskian non-degeneracy hypotheses) excludes it, I could not determine: NOT_CHECKED.

## (e) Filter F3+

> **F3+.** Call an argument *char-free* if every input it uses holds verbatim for F_2[t] with the anisotropic norm of F_4(t)/F_2(t), where deg plays the role of log and the inert degree-one place ∞ plays the role of the complex place. Admissible inputs are: unique factorisation, exact norms, residue counts, pigeonhole at small primes, cut and chord identities, the chain lemma, Liouville/toric forced divisibility, parity/anisotropy of Galois, Plücker relations, Schwartz–Zippel, and intersection theory on surfaces short of BMY. **No char-free argument can prove the uniform bound, or any bound M ≤ f(W) with f(W) = o(√(W/log W)).** Arguments for growth-rate improvements (M ≤ εW/log W) are **not** touched: the fixed-field growth question (b) is open.

Tests against the record:

| Argument | Fate |
|---|---|
| Theorem 1.1, 3M log M | transfers; it is a growth bound and survives |
| (IV.2)/BMY | uses characteristic 0, so it passes, but it has no integer analogue (arithmetic BMY is at least Szpiro strength) |
| J4 (A_5 ⇒ uniform), gauss reduction A_5 ⇒ (L_5°) | char-free, so A_5 cannot be proved char-free |
| (L_8) ⇒ uniform | the same |
| Liouville/K1-class methods, anisotropy/parity lemma | sharp in characteristic 2 too (15/4 attained); killed |
| route 5's zero-estimate step | killed if algebraic and char-free (NOT_CHECKED) |
| Vojta routes | survive only because number-field Vojta is not char-free; any "function-field style" proof must use characteristic 0 (a Kodaira–Spencer/semipositivity input) |
| F3 (round 31 IV) | its "or small-residue collisions" escape clause is now dead for uniform statements: residue fields here are F_{2^e} with e small, yet clusters are unbounded |

What survives: inputs genuinely special to number fields, namely the archimedean real structure of the angles, characteristic 0, and Baker-type estimates (their char-p status was not checked).

## Novelty and significance

* **Novel against the record:**
  * a fixed-field anisotropic unbounded cluster (round 31 IV had only growing q and the split form);
  * an explicit characteristic-2 rational-surface BMY violation coming from genuine clusters;
  * FF-A_5 and FF-(L_5°) refuted at the Liouville exponent;
  * F3+.

  A grep of the record finds no subspace-polynomial construction outside round 31.
* **Significance for the goal: low (negative calibration).** Nothing changes for integers: the (2/3+ε) bound stands, and Questions A and B are open over Z. What Part V does is make "Questions A and B need characteristic 0 or genuinely archimedean input" a theorem-level calibration. It is the sharpest no-go in the record for char-free uniform proofs. Score: 3/10.

## Verdicts

1. (a) Dictionary → (IV.2) instance with δ = E = 0 and m > 24; violation at m = 32 and m = 64; Langer fails; 3e_orb − (K+B)² < 0: **CORRECT**.
2. (b) Transferred bound W ≥ 3M log₂M − O(M); Part V consistent at about (2/3)m times the bound: **CORRECT**.
3. (c) The FF analogues of A_5(δ < 1/30), (L_5°) and V_E(5) fail; α → 15/4 from above: **CORRECT** (via the J4 + Schwartz–Zippel transfer).
4. (c) The FF-(L_8) analogue fails: **CORRECT_WITH_CORRECTIONS** (asymptotic only, m ≫ 896/c; dictionary not re-derived in characteristic 2).
5. (d) Not a Frobenius/inseparable instance; separable additive (Artin–Schreier) mechanism: **CORRECT** (separability computed). Whether it is covered by a published characteristic-p Vojta formulation: NOT_CHECKED.
6. Astra's V.2 bounds (V.3)/(V.4): NOT_CHECKED (outside this block).

Sources: [Easton, Surfaces violating BMY in positive characteristic, arXiv math/0511455](https://arxiv.org/pdf/math/0511455); [PAMS 136 (2008)](https://www.ams.org/proc/2008-136-07/S0002-9939-08-09466-5/); [Algebraic surfaces in positive characteristic, arXiv 0912.4291](https://arxiv.org/pdf/0912.4291).
