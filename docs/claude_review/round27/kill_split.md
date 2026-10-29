# Adversarial review (killer): round 27, lane split

Target: /tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/split.md
Scripts: /tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/kill_split_checks/

## 0. Verdict

**The obstruction survives with corrections. Novelty: a repackaged known barrier**, namely the section/Liouville
`τ`-cap (Theorem B, Lemma K1, `τ_e = 10/3`) carried to asymmetric split profiles. Its curve form is the standard BDPP
duality on `Bl_e LM_5`, a contraction of `M̄_{0,7}` (Castravet–Tevelev). There is one genuinely new calibration fact:
the `7/2` quadric cones, which must be completed by their inversions. No mechanism beats `15/4` anywhere, and nothing
ports to `T_5`, as the write-up itself says.

What I re-derived or recomputed with my own code (`kill_split_checks/`):
* *Exactly:* Thm 3.1 (symbolic and exact); `dim Π_nt = 11`, the 278 one-sided vertices, the code hull; Lemma 1.1, the
  cost identity and Thm 3.2; Lemma 3.5 (exact `V(y)` formula on 200 tuples); Thm 3.6 (200 tuples); the phase-0
  `k`-tuples.
* *Numerically and independently:* Cor 3.8, refined to 60 digits.

Defects, most serious first:
1. **Proposition K (new, proved here) refutes the write-up's forward claims at the even-set vertex.**
   * At the `S_5`-invariant vertex `x_{[5]} = 1[|B| even]` (and its inversion) there is **no rational curve of class
     `γ_{x_{[5]}}` with `mult_e ≥ 4`**.
   * So Prop 3.7's remark that covering families are "expected at every rational near-tight profile" (expected
     dimension `3N+1 ≥ 4`) is false at `N = 1`.
   * The proposed next step Q-split-1′ ("first: the even-set vertex `x_{[5]}` … a Pell-type construction with profile
     `x_{[5]}`, or a degeneration") is impossible at `N = 1`.
   * The author's diagnosis "a solver limitation" is right for 14 classes and wrong for this pair. Their own failure
     at the `W`-class came from their solver, while for this pair it is a theorem.
2. **Defect 3: an overclaim on fixed-`S` Subspace.** "Subspace with rational forms and no special `p`-adic behaviour
   is Liouville" is false in general (`S`-unit equation). Lemma 3.5 covers only single Δ-pole functions, so novelty
   item (b) must be cut back.
3. **Missed literature, affecting novelty.** WTB is a moving-curve statement on `Bl_e LM_5`, a contraction of
   `M̄_{0,7}`. Prop 3.7 is the easy half of BDPP. The phase-0 reading is Vojta's Büchi paradigm.
4. **Defect 1: incomplete exceptional set.** Split inversion `d ↦ n/d` preserves A-points. The mnn31 co-divisor
   tuples lie on none of the 15 cones (exact), so `Z` must also contain the 15 cubics `ι(𝒬_σ)` (and `ι(W)` if it
   differs from `W`).
5. **Defect 2: false display.** In the proof of Thm 3.2, `cost_{1/2} = (1/8)Σ_P w_P` should read `(1/4)Σ_P w_P`; the
   conclusion is unaffected.

**Strengthening found along the way (numerical).**
* Dominant four-branch families exist at **14 of the 16** one-sided classes (complex-start damped Gauss–Newton;
  the author's random-start solver, whose start distribution I did not inspect, failed even at the `W`-class).
* Those 276 vertices span `Π_nt`, so WTB holds numerically on an **11-dimensional** neighbourhood of `x = 1/2`, not
  just on the 10-dimensional slab `Q_code`.
* The remaining open spot is a neighbourhood of `x_{[5]}` and of `1 − x_{[5]}` (and the non-0/1 vertices). There
  WTB is genuinely undecided, and the `N = 1` curve route is provably closed.

## 1. Theorem 3.1 (mnn31 curve, cone)

**Re-derived and recomputed (own code: `k_thm31.py`, `k_thm31b.py`). The theorem survives, but its conclusion is incomplete.**

* *Symbolic identities (exact).* I parametrised the conic `u(u+a+b) = k w(w+c+d)` (`k = ab/cd`) by the lines `w = tu`
  through `(0,0)`. Then `u(u+a+b) − k w(w+c+d) ≡ 0`, `(u+a)(u+b) ≡ k(w+c)(w+d)` and `Y_2Y_3 − Y_4Y_5 ≡ 0`, with `a,b,c,d,t`
  symbolic. The cone identity is proved.
* *Profile (exact).* I factored every `D_j(t)` over `Q` at `(a,b,c,d) = (2,3,1,2)` (mnn31, `k = 3`) and at
  `(3/7, 5/2, 2/3, 11/5)`, and computed the valuation vector of `(D_1..D_5)` modulo the diagonal at every place.
  There are seven rational places with 0/1 vectors on seven distinct cuts, each hit once (`B·C = 7`). The degree-2 place
  `1 − kt² = 0` has vector `0` and maps to `e` (all `Y_j → 0`). For mnn31 the sides are exactly the author's
  `{4,5},{2,3},{2,3,4,5},{1,2,4},{1,2,5},{1,3,4},{1,3,5}`.
* *Exponent (exact integers).* I generated 58 Pell points of `(2u+5)² − 3(2w+3)² = −2`. At every one of them the five
  numbers divide `n`, are distinct, satisfy the conic, and lie on cone no. 0 (exact). `log(dist^{7/2}·H)` (primitive
  tuple, `H = (max/min)·lcm`) is periodic with values `{4.064, 4.267, 4.960}` up to `u ≈ 10^{14.7}`, and
  `C_actual → 3.4641 = 2√3`. So `α = 7/2` along the curve. (The author's "`lcm/n → 7/8`" means `log lcm/log n → 7/8`;
  I observe `0.861…0.900`, oscillating around `0.875`.)
* *Family (exact).* The Jacobian of `(a,b,c,d,t) ↦ (Y_2..Y_5)` has exact rank 3 at 4 random rational points. With the
  cone irreducible of dimension 3, the union of the curves is dense in it.
* *A-points.* Dirichlet approximants to the real quadratic parameter `t_0 = ±1/√k` give `dist ≍ H^{-2/7}`. Since
  `2/7 > (1+2δ)/4` iff `δ < 1/14`, each curve with `k > 0` non-square carries infinitely many A-points at every
  `δ < 1/14`. The "Consequence" (every admissible `Z` contains the cone and its 15 `S_5`-images) is **correct**.
* *m-czech.* The generic W-curve lies on none of the 15 cones, checked symbolically by reducing the cone polynomial
  modulo the conic at `(a,b,c) = (1,2,4)` and `(3/2,5/3,7/5)`. The integral m-czech points with `x ≥ 7` lie on none of
  them. (The points with `x ≤ 4` lie on several cones only because they are degenerate: some `F_j` vanishes.)

**Defect 1 (minor, but it corrects a stated result): the exceptional set must also contain the inverted cones.**
Split inversion `u_a ↦ lcm(u)/u_a` (divisor `d ↦ n/d`) maps primitive tuples to primitive tuples. It preserves lcm,
`H` and `dist` up to `O(1)`, so it maps A-points to A-points, and an admissible `Z` can always be replaced by
`Z ∩ ι(Z)`. The co-divisor tuples of mnn31 have `gcd = 1` and the same lcm (exact, all 25 points tested), and they lie on
**none** of the 15 cones (exact). So every admissible split `Z` contains, besides the 15 cones `𝒬_σ`, the 15 cubic
threefolds `ι(𝒬_σ)`: `u_1u_2u_3 − u_1u_4u_5 − u_2u_3u_4 − u_2u_3u_5 + u_2u_4u_5 + u_3u_4u_5 = 0` and its images
(`sympy.factor` of the inverted quartic gives `−u_1·(cubic)`). If `ι(W) ≠ W`, as the author's Prop 3.3(e) asserts
numerically, then `Z` must also contain `ι(W)`; round 23 states only `Z ⊇ W`.

The author's twisted remark is still right. Over `Q(i)` the real structure is `u ↦ 1/ū`, so `σ(𝒬) = ι(𝒬)` as
equations. `𝒬 ∪ ι𝒬` is therefore a `Q`-subvariety of the twisted `T_5` whose rational points lie in `𝒬 ∩ ι𝒬`, a
surface. In the split form both threefolds are forced into `Z`; in the twisted form at most this surface is. This is the
same phenomenon as the parity lemma. It belongs to the exceptional set and is not a mechanism.

**Status.** The theorem is correct, and the exceptional-set statement should read "contains the 15 cones and the 15
cubics `ι(𝒬_σ)`". The construction is mnn31's (Lean, 2026). The cone and the exponent `7/2` are a new observation for
the record. It is a calibration fact (it enlarges `Z`) and not a mechanism. Neither the cones nor the cubics are
cut-balanced (`B·C = 7`), so they do not touch `(L_5^split)` or the twisted Question B.

## 2. Profile polytope, cost identity, Theorem 3.2

**Recomputed (`k_profiles.py`, exact).** The 15 cut-balance equations and the 4 equations `X_a = X_5` have rank 19, so
`dim Π_nt = 11`. There are exactly 278 one-sided 0/1 near-tight profiles in 16 `S_5`-classes. The 20 code vertices
`x_R`, `1 − x_R` (`|R| = 3`) span an affine 10-space. Their uniform average is exactly `1/2`, so `1/2` is in the
relative interior of the hull. Adding `x_{[5]} = 1[|B| even]` and its complement gives affine dimension 11. All of this
agrees with the author.

**Is the "gap in the record" real? Yes.** Round 23 `split.md` §2 proves that the pair (binomial) bound *attains* `15/4`
at every cut-balanced profile, by averaging over the 20 ordered pairs. That is a lower bound on `α`. For the *cap* it
says only "the section cap is the same, because the cap is a statement about Laurent polynomials and widths and holds
over any field". At a split point of profile `x`, forced divisibility gives the one-sided cost
`Σ_B x_B π_B(S)`, not the width sum `(1/2)Σ_cuts w_B(S)`. The two agree at `x = 1/2`, and also for every `x` when `S`
is symmetric. So the record's transfer of the cap to cut-balanced split profiles was unproved off `x = 1/2`, and the
author has located a real, if small, gap. I checked that the translation invariance needed for this to be well posed
holds: on `Π_nt`, `cost_x(S + c') = cost_x(S) − Σ_a c'_a X_a = cost_x(S)`.

**Lemma 1.1, cost identity, Theorem 3.2: correct, with one false display.** I re-derived all three.
* `P_ba = 8 − P_ab` and `P_ab − P_ba = X_a − X_b` hold.
* `cost_x = cost_{1/2} − Σ_B(x_B − 1/2)mid_B` holds.
* For `S = 2c − S`, `Σ_B x_B c(B) = Σ_a c_a X_a = 0` on `Π_nt`.

**Defect 2 (minor; false display, true conclusion).** The proof of Thm 3.2 writes
`cost_{1/2}(S) = (1/8)Σ_{P⊆[5]} w_P(S) ≥ 4·reg(S)`. The correct identity is `cost_{1/2} = (1/2)Σ_cuts w_B = (1/4)Σ_P w_P`.
Theorem B (`threshold/construct.md` §2.4: `reg(S) ≤ 2^{1−n}Σ_{P⊆[n]} w_P(S)`, here `n = 5`) gives
`4 reg ≤ (1/4)Σ_P w_P`. With `(1/8)` the display would need `reg ≤ (1/32)Σ_P w_P`, which is false: a root binomial has
`reg = 1` and `Σ_P w_P = 16`. The binomial bullet, `cost = (1/4)Σ_{B}|m(B)|` over 30 sets, `= 4` for roots, uses the
right constant. So the theorem stands once the display is corrected.

**Assessment of Theorem 3.2.** It is a two-line consequence of linearity. The near-tight conditions `X_a = const` are
exactly the balancing conditions that make `x` a curve class (see §3), and a symmetric support has midpoint function
`mid_B = c(B)`, which is linear in `B` and pairs to zero with every curve class. It is correct. As mathematics it is
book-keeping that closes the corner of the record's barrier map for symmetric supports. It is not a new mechanism.

## 3. Prop 3.7 curve criterion, Prop 3.3, Cor 3.8 (W-vertex, Q_code)

**3.1 Prop 3.7 is correct, and it is the standard duality between moving curves and effective divisors.** Its proof
is right: zeros over `e` are at least `ord_e(f)·mult_e C`, and poles are at most `Σ_P h_S(−v_P) ≤ Σ_B N x_B π_B(f)` by
sublinearity of the support function on the smooth fan.

Two points the write-up does not make:
* The condition `X_1 = … = X_5` defining `Π_nt` is exactly the balancing condition `Σ_B x_B 1_B ∈ R·1`. So `Π_nt` is
  the set of cut-balanced curve classes `γ_x` on `LM_5`, with `E_B·γ_x = x_B`.
* Write `D_f = Σ_B π_B(f)E_B − ord_e(f)E` on `Bl_e LM_5`. Then `D_f·(π^*γ_x − 4ℓ) = cost_x(f) − 4ord_e(f)`, where
  `ℓ` is a line in the exceptional divisor. So WTB at `x` says exactly that `π^*γ_x − 4ℓ` is nonnegative on every
  effective divisor, i.e. lies in the closed movable cone of the smooth projective variety `Bl_e LM_5`
  (Boucksom–Demailly–Păun–Peternell). Prop 3.7 is the easy half of that duality. By BDPP the criterion is also
  necessary up to closure and multiples: WTB at `x` holds iff `π^*γ_x − 4ℓ` is a limit of classes of covering
  families. These covering families may have any genus; Prop 3.7 as stated uses only rational curves.
* *Literature link missed by the author and by the record (grep: no hit in `docs/claude_review`).* My own derivation:
  the Hassett reduction `M̄_{0,7} → LM_5` (points `0,∞` heavy) has fibre over `e` equal to the Cartier divisor
  `D_{0∞} = D_{12345}`, so it factors through `Bl_e LM_5`. It contracts exactly the 15 divisors `D_T`, `T ⊆ [5]`,
  `|T| ∈ {3,4}`; the Picard ranks agree, `42 − 15 = 27 = 26 + 1`. This is Castravet–Tevelev's morphism
  `M̄_{0,n} → Bl_e \overline{LM}_n` (Duke 2015, arXiv 1311.7673, Thm 1.1, seen through search results only; their `LM_n`
  has dimension `n − 3`). Consequently `τ_e(−K) = 10/3`, Theorem B and WTB are statements about the effective cone of
  `M̄_{0,7}`. The known extremal effective divisors there (boundary, hypertree divisors, and the non-hypertree examples
  in the literature, recalled from memory and not re-checked) are the natural test set for WTB. The cone itself is not
  known to me to be determined for `n = 7`.

**3.2 Prop 3.3(a,b): correct** (same zero/pole count on the `W`-curves; Newton polytopes add). (c)–(e) are numerical or
derived and are not needed below.

**3.3 Cor 3.8: independently reproduced, and extended (numerical).** My formulation is my own (`k_cover.py`).
* *Set-up.* For a one-sided 0/1 vertex `F`, put `u_a(s) = c_a Π_{B∈F, B∋a}(s − r_B)` with `c_5 = 1`. The four branch
  parameters are `s = 0, 1, −1, s_4` (this fixes the Möbius gauge). The equations are
  `c_a Π_{B∈F: a∈B∌5}(s_i − r_B) = Π_{B∈F: 5∈B∌a}(s_i − r_B)`: 16 equations in 20 unknowns. On a one-sided near-tight
  vertex `P_{a5} = 4`, so every equation has four factors on each side.
* *Method.* Complex random starts, damped Gauss–Newton (minimum-norm steps), at most 40 s per class.
* *Result.* Non-degenerate solutions were found at **14 of the 16** one-sided classes, including the `W`-class and its
  `ι`-image. Non-degenerate means root separation `≥ 3·10^{-3}`, `c_a ∉ {0,1}`, and the 16×20 Jacobian of rank 16
  (smallest singular value `≥ 4·10^{-3}`). In every case the dominance matrix (the tangent of the 4-dimensional
  solution set together with `∂_s`, mapped to `log(u_a/u_5)`) has four nonzero singular values. The smallest ranges
  from `3·10^{-4}` to `0.38`.
* *60-digit refinement at the `W`-class* (`k_refine.py`: augmented square system, mpmath).
  * Newton residuals go `1.3·10^{-14} → 5.5·10^{-27} → 6.4·10^{-53}`, so convergence is quadratic.
  * The 16×20 Jacobian has smallest singular value `0.109`.
  * The dominance singular values are `(6.28, 1.09, 0.242, 0.0289)`.

  This is far from degenerate. It is still not interval-certified, so the label stays "numerical", but it is
  independent of the author's code and much sharper than "8/8 targets".
* So the author's statement that random-start Newton fails at every class "even at the `W`-class" describes their own solver
  (I did not inspect its start distribution); complex random starts succeed at once.
* **Extension (exact hull computation, numerical input; `k_hull14.py`).** The 276 one-sided vertices in those 14
  classes span the full 11-dimensional affine hull of `Π_nt`, and their average is exactly `1/2`. By linearity and
  `S_5`, WTB therefore holds (numerically) on an **11-dimensional** neighbourhood of the symmetric profile, not only on
  the author's 10-dimensional slab `Q_code`.

**3.4 New obstruction (proved here): the even-set vertex has no `N = 1` rational curves at all.** The two classes where
both searches fail are exactly the `S_5`-fixed pair `x_{[5]} = 1[|B| even]` and `ι x_{[5]} = 1[|B| odd]`. In
`k_cover2.py`, 600 s per class, every converged start was degenerate, with coincident roots or `c_a ∈ {0,1}`. This is
not a solver failure.

> **Proposition K.** There is no irreducible rational curve `C ⊂ LM_5` not contained in `∂ ∪ ⋃Δ̄_ab` with
> `E_B·C = 1` for `|B|` even, `E_B·C = 0` for `|B|` odd, and `mult_e C ≥ 4`. By `ι`, the same holds for the odd
> profile.

*Why the class determines the shape.*
* *Distinct places.* Places with non-nested types cannot merge without changing the class. A merged place
  `r_{ab} = r_{ac}` or `r_{ab} = r_{cd}` has a valuation vector outside every chain cone, so it meets different
  divisors. Likewise `ρ_a = ρ_b` cannot occur. So all `r_{ab}` are distinct, and all `ρ_a` are distinct.
* *Nested merges are harmless.* The only merge that keeps the class is `r_{ab} = ρ_e`
  (`{a,b} ⊂ [5]∖e`). It is harmless below.
* *Birationality.* The normalisation map is birational, since `γ_{x_{[5]}}` is primitive.

So the hypotheses used in the proof are automatic. Exact random-instance checks follow the proof.

*Proof.*
1. *Normal form.* Write `u_a ∝ φ_a(s) = c_a Q_a(s)/(s − ρ_a)`, with `Q_a = Π_{j≠a}(s − r_{aj})`. Here `r_{aj}` is the
   place on `E_{\{a,j\}}`, `ρ_a` the place on `E_{[5]∖a}`, and `c_a ≠ c_b` (else `C ⊂ Δ̄_ab`).
2. *A common quartic `q`.* `φ_a − φ_b` has numerator of degree 5, divisible by `s − r_{ab}`, with leading coefficient
   `c_a − c_b`. It vanishes to total order `≥ 4` over `e`. Hence
   `φ_a − φ_b = (c_a − c_b) q(s)(s − r_{ab})/((s − ρ_a)(s − ρ_b))` with one monic quartic `q` for all pairs. This
   includes the cuspidal case.
3. *Normal form for the `φ_a`.* Put `κ_a = Res_{ρ_a}φ_a/q(ρ_a)`. Then `φ_a − κ_a q/(s − ρ_a)` is a polynomial of degree
   `≤ 3`, and it does not depend on `a`: the differences have no poles and vanish at `∞`. So
   `φ_a = Ψ + κ_a q/(s − ρ_a)` with all `κ_a ≠ 0`.
4. *Passing to a plane configuration.* `r_{ab}` is a common root of `N_a = Ψ·(s − ρ_a) + κ_a q` and of `N_b`. In the
   plane `(y, w = ys)` this says that the node `P_{ab}` of the lines `L_a: w = ρ_a y + κ_a` and `L_b` lies on the
   quartic `G = y⁴q(w/y) + y³Ψ(w/y)`, which has order `≥ 3` at `O = (0,0)` and is nonzero. Non-degeneracy (all `r`
   distinct) means: distinct slopes, no three lines concurrent, no line through `O`.
5. *Quartics through the ten nodes.* Restricting to `L_5` and inducting (cubics through the 6 nodes of 4 lines, conics
   through a triangle), the quartics through the 10 nodes are exactly `span{Λ_c = Π_{b≠c}L_b}`.
6. *The contradiction.* Then `G/ΠL_b = Σ_c λ_c/L_c` must vanish to order 3 at `O`. With `v_c = (ρ_c, −1)/κ_c` this
   says `Σ_c μ_c(1, v_c, v_c⊗v_c) = 0`, i.e. the five points `v_c` impose dependent conditions on conics. That forces
   four of them to be collinear, i.e. four lines `L_c` concurrent or parallel, a contradiction. `□`

Exact random-instance checks (`k_even_thm.py`, 5 instances): quartics through the nodes form a 5-dimensional space,
those with a triple point at `O` a 0-dimensional one, and the Veronese rank is 5. With four concurrent lines the
triple-point space becomes 4-dimensional, as the proof predicts.

**Consequences for the write-up.**
1. *The "expected at every profile" claim is false.* Prop 3.7's remark ("covering families are expected at every
   rational near-tight profile", expected dimension `3N+1`) fails at `N = 1` on the vertex `x_{[5]}`: the expected
   4-dimensional family is empty.
2. *The suggested route to WTB is closed at `N = 1`.* Q-split-1′ proposes `x_{[5]}` first, via "a Pell-type
   construction with profile `x_{[5]}`, or a degeneration". Proposition K shows no such rational curve (Pell or not,
   smooth branches or cusps) exists at `N = 1`.
3. *What remains open.* WTB at `x_{[5]}` needs `N ≥ 2`, higher-genus covering families (an elliptic family of the same
   class has expected dimension 3, which is enough to cover), or a divisor-side proof. Otherwise it fails, and a
   section beats `15/4` at the even profile.
   * *The split into pairs and quadruples fails too.* The obvious decomposition
     `γ_{x_{[5]}} − 4ℓ = (γ_{pairs} − 3ℓ) + (γ_{quads} − ℓ)` is not realised by rational curves either. A rational curve
     of class `γ_{pairs}` (`B·C = 10`) with `mult_e ≥ 3` leads, by the same argument, to a quartic `y·q(s) + Ψ(s)`
     through the 10 nodes of the lines `y = c_a s − κ_a`, and to the same Veronese contradiction (derived; not
     separately checked by computer).
   * *Searches at `N = 2`.* My `N = 2` searches (`k_coverN2.py` and `k_coverN2b.py`, 39 unknowns, 32 equations, about
     1300–2300 starts per run; the second variant also allows tangency to a single `E_B`) found no non-degenerate
     solution, either at `x_{[5]}` or at the `W`-class **control**. The `N = 2` negative is therefore uninformative.
   * *Sections at `x_{[5]}`.* `S_5`-symmetric supports `{m_a+m_b ≥ −α, m_c ≤ β}`, `α ≤ 6`, give
     `max 4·reg/cost_{x_{[5]}} = 1` exactly (at `(α,β) = (3,2)`, where the Vandermonde is tight;
     `k_even_reg.py`, two primes). No violation was found.
4. *The even profile has no analogue of `W`.* There is no `N = 1` Pell curve through `e` with two cusps at profile
   `x_{[5]}`.

## 4. Bounded residues: Lemma 3.5, Theorem 3.6

**Recomputed (`k_misc.py`, exact on 200 random primitive 5-tuples, `n_B` computed by factorisation).**
* The exact Vandermonde formula `V(y) = Π_{a<b} t_ab · Π_B n_B^{−|B|(5−|B|)/2}` holds 200/200. This is the identity
  behind Lemma 3.5.
* `Σ_B x_B |B|(5−|B|)/2 = 40` holds on every cut-balanced profile (it is `(1/2)Σ_{a≠b}P_ab`); near-tightness is not
  needed.
* Theorem 3.6, `x^τ × y^τ = (t_bc, t_ca, t_ab)`, holds for all 10 triples of all 200 tuples. It is a one-line identity
  (`α_ab = x_a y_b`).

**Lemma 3.5: correct as stated.**
* (i) With `m ∈ (k·cost/40, k·ord/10)` the forced numerator `Π n_B^{m f(|B|) − kπ_B(g)}` makes `|f(y)| ≥ X^{−cost_x(f)}`
  with `cost_x(f) < 0`. Together with the archimedean bound this forces `g(y) = 0` for large `X`.
* (ii) is the mediant inequality.
* The general case of Δ-poles (`g/Π Δ_ab^{m_ab}` with unequal exponents) also follows: on `Π_nt` each `Δ_ab` has
  order 1 and cost `P_ba = 4`.

**Defect 3 (major for the novelty claim, not for any theorem): overclaims about fixed-`S` Subspace.**
* *The false principle.* The write-up says "Corvaja–Zannier-type arguments with `S = S_K` and rational linear forms
  reduce to Lemma 3.5 (the Subspace Theorem with rational forms and no special `p`-adic behaviour at `p ∈ S_K` is
  Liouville)". Novelty item (b) builds on it: "Integrality at the root arrangement … the natural fixed-`S`
  Siegel/Corvaja–Zannier input, adds exactly nothing beyond sections". The principle in parentheses is false as a
  general statement. Finiteness for the `S`-unit equation `x + y = 1` is the Subspace Theorem applied with the
  *rational* forms `X, Y, X + Y`, and it is not a Liouville statement.
* *What Lemma 3.5 actually proves.* A **single** regular function on `X°` with Δ-poles, used through the product
  formula, gains nothing over sections where WTB holds. Corvaja–Zannier and Levin arguments apply the Subspace Theorem
  to many such functions at once, at all places of `S`. Lemma 3.5 does not touch them.
* *Corrected scope of (b).* Δ-pole sections add nothing, because `Δ_ab` has ratio exactly 4 on `Π_nt`. Whether a
  genuine fixed-`S` Subspace argument on `X°` adds anything is open. For instance, Levin's general-position theorems do
  not apply as stated: the `Δ̄_ab` are not ample and are not in general position, since
  `Δ̄_ab ∩ Δ̄_bc ⊂ Δ̄_ac`.
* The record's papers check marks the analogous claim for `S = {∞}` as an inference to be checked (`papers_check`
  §3 item 6).
* There is no harm to the proved statements; the defect is in the barrier-map item (b).

**Prop 3.4: correct.** One small point: "`[u]` is `S_K`-integral, i.e. every `z_ab ∈ Z[1/S_K]`" identifies integrality
with integrality of the `z_ab`. That is fine at the level of heights (finite part of `h_{D_ab} = log|t_ab|`). For a
literal integral-model statement one would have to name generators of `O(X°)`.

## 5. Phase 0

* *The identity.* Re-derived: `z_j² − z_1² = 2k_js_1 + k_j²` follows from `s_j² − z_j² = 4n`.
* *The two families* (`k_misc.py` and by hand). The m-czech family gives `k = (0,224,225,324,900)` at all 13 Pell points
  with `x ≥ 37`, `x < 10⁶` (the identity holds exactly, and `C → 30`). For mnn31, `s_j − s_1 = 0, 4, 9, 3, 12` by
  direct expansion on the conic, so `k = (0,3,4,9,12)`.
* *The window for mnn31.* By hand, `D_1 − √n > 3` and the largest offset is `≈ 6w ≈ 2√3·n^{1/4}`. All five divisors lie
  in `(√n, √n + (7/2)n^{1/4})` (and `C_act → 2√3`, see §1).
* *The canonical class.* `K_{S_k} = O(M − 4)` for a complete intersection of `M − 1` quadrics in `P^{M+1}` is standard
  adjunction. The author correctly flags smoothness as unchecked.
* *Novelty.* The general-type-surface plus Lang reading is the paradigm of Vojta's conditional treatment of Büchi's
  problem: Vojta (1999/2000, Contemp. Math. 270) proved that a weak Bombieri–Lang conjecture implies a positive answer
  for systems of diagonal quadrics, i.e. "linear forms simultaneously squares". This comes from search results only;
  Vojta's chapter was not read. The author labels this material "rediscovered, not claimed" in the summary, which is
  appropriate. I would add the Vojta–Büchi pointer.
* *Literature.* The erdosproblems.com snippets I saw still say "infinitely many `n` with 4 divisors … whether 4 is best
  possible" (search index, date unknown). The five-divisor constructions rest on the record and the fetched Lean file.
  I verified the mnn31 construction mathematically and exactly (§1), so the mathematics does not depend on the source.
  The search snippet shows Chan's square case with a width `n^{1/4}(log n)^{1/7}` interval, slightly different from the
  author's paraphrase.

## 6. Pre-registered tests T1-T5 (independent)

**T1 (fakes, growth): vacuous pass.** No growth claim is made. Theorem 3.2, Prop 3.3, Prop 3.7, Lemma 3.5 and WTB
say which auxiliary forms *fail*, so fakes satisfy them trivially. Theorem 3.1 and Proposition K are about curves
over `Q`/`C` and say nothing about fakes. Nothing in the write-up is an input that a fake could violate.

**T2 (`k = 4` control): pass.**
* *No exclusion fires.* The only exclusion mechanisms in the write-up are conditional on a section beating the cap:
  Lemma 3.5(i) and Prop 3.3(b). The dense `k = 4` family of gauss Thm 4.1 is twisted, so it sits at the symmetric
  profile. There Theorem B with `n = 4` gives `2·reg ≤ cost_{1/2}` (cap `7/2`), so neither mechanism fires.
* *The `k = 4` analogues hold.* I re-derived the `k = 4` analogue of Thm 3.2 (same proof). For the `k = 4` analogue of
  Lemma 3.5:
  * each pair is separated by 4 sets `B`, so `P_ab + P_ba = 4` on cut-balanced profiles;
  * hence `Σ_B x_B|B|(4−|B|)/2 = (1/2)Σ_{a<b}(P_ab+P_ba) = 12`, against `ord_e V_4 = 6`;
  * the cost/order ratio is `2`, i.e. exponent `7/2`, which is exactly the balanced `k = 4` cap. Consistent.
* *Phase 0 says nothing at `k = 4`.* For four divisors `S_k` is K3, so Lang's conjecture is silent.

**T3 (split / #887): pass, with one correction.**
* *Nothing that would imply #887 is claimed.* Split Question B is not claimed, so absolute #887 is not claimed.
* *The 2026 constructions land in the exceptional set.* The m-czech curve lies on `W` and on none of the 15 cones
  (exact, symbolic). The mnn31 curve lies on cone no. 0 (exact, 58 Pell points).
* *Correction (Defect 1).* The co-divisor tuples of mnn31 are A-points with the same lcm and lie on none of the cones.
  So `Z` must also contain the 15 cubics `ι(𝒬_σ)`, and `ι(W)` if it differs from `W`.
* *No exclusion of these points.* Nothing in the write-up excludes any of them.

**T4 (criticality): the write-up is inside the mapped section/Liouville class, as it says itself.**
* *Thm 3.2* is Lemma K1/Theorem B book-keeping at asymmetric profiles. For symmetric supports the profile drops out,
  because the near-tight conditions are the balancing conditions of curve classes.
* *Prop 3.7* is the elementary half of BDPP duality on `Bl_e LM_5`, which is a contraction of `M̄_{0,7}`
  (Castravet–Tevelev).
* *Lemma 3.5* is class L with Δ-poles. The Δ-poles are toric again by the exact formula
  `V(y) = Πt_ab·Π n_B^{−f(|B|)}`.
* *Multi-point arguments* are not used.
* *Fixed-`S` Subspace* is not used in any proof. The write-up's assertion that it reduces to Lemma 3.5 is an unproved
  and over-general claim (Defect 3), and novelty item (b) must be cut back accordingly.
* *Where something could escape the cap.* The one place where a section *could* escape the mapped caps is an
  asymmetric split profile where WTB fails. My Proposition K shows that the author's proposed curve route cannot close
  the even-set vertex at `N = 1`, so this possibility is more alive there than the write-up suggests.
  * It cannot port to `T_5`: the twisted near-tight profile is `1/2`.
  * So even a failure of WTB would be split-only, and it would bear on absolute #887 only through split Question B
    restricted to that profile.

**T5 (calibration): pass.**
* *Points already known.* Nothing excludes existing points. Round 24's real twisted family sits at `x = 1/2`, where
  Theorem B gives no exclusion. The split Pell `W`-points (profile `x_R`) are protected: I reproduced the dominant
  four-branch family at `x_R`, so no section beats `15/4` there. The cone points (`α = 7/2`) are put into `Z`.
* *Other calibrations.* Huang's `Y_4` and `P²` are not touched (no exclusion is claimed on any variety). The record's
  8- and 10-point families are not affected (no growth claim).
* *Proposition K.* It is a non-existence statement for rational curves of one class at mult `≥ 4`. No known point or
  curve has that class (the known split curves have `B·C = 7` or profile `x_R`), so it excludes nothing that exists.

**Handoff constraints.** Respected by the proved parts. The fixed-`S` sentence (Defect 3) does not assume a uniform
exceptional set; it overstates what Subspace reduces to, which is a different error.

## 7. Literature claims

The author's literature table was checked against what WebSearch shows today (arXiv itself does not resolve in this
session; snippets only).

| claim | check | verdict |
|---|---|---|
| mnn31 five-divisor construction, `C = 7/2` | not visible in search snippets. The erdosproblems.com 887 snippets still show the 4-divisor status (index date unknown). The author's source is the fetched Lean file | the mathematics was verified exactly here (§1, §5), so the claim is safe regardless of the source |
| m-czech `C = 30`/`31` | record (round 23); re-verified exactly here (`C → 30.0000`, `k`-tuple) | fine |
| Chan, squares, at most five divisors | snippet: "at most 5 divisors in an interval of width about `n^{1/4}(log n)^{1/7}`" | close to the author's paraphrase; the width differs slightly |
| "Subspace with rational forms … is Liouville" | no source given | **false as a general principle** (the `S`-unit equation); see Defect 3 |
| curve criterion (Prop 3.7) as new | standard (BDPP duality / covering families; Castravet–Tevelev for `M̄_{0,n} → Bl_e LM_n`) | not new; missing citations |
| phase-0 Lang reading | Vojta's conditional Büchi theorem (Contemp. Math. 270, 2000; seen via Wikipedia/arXiv 2412.16740 snippets) is the same paradigm | the author says "rediscovered", which is fine; add the pointer |

Sources consulted (search results only): arXiv 1311.7673 (Castravet–Tevelev, "M̄_{0,n} is not a Mori Dream Space",
Thm 1.1 as summarised in search results), arXiv 2009.14298 (Castravet–Laface–Tevelev–Ugaglia), erdosproblems.com/887
and its forum thread, en.wikipedia.org/wiki/Büchi's_problem, arXiv 2412.16740.

## 8. Novelty decision

**Verdict: repackaged known barrier, plus one new calibration fact.** The barrier being repackaged is the
section/Liouville `τ`-cap: Theorem B, Lemma K1 and the `τ_e = 10/3` cap, transported to asymmetric split profiles,
together with BDPP duality on `Bl_e LM_5`.

The author's novelty items, one by one:
* **(a) 11-dimensional profile space; symmetric supports are profile-independent; the curve criterion; the slab
  `Q_code`.**
  * The profile space and the gap in round 23's transfer are genuine observations about the record.
  * Theorem 3.2 is two lines of linear algebra.
  * Prop 3.7 is the standard covering-family half of BDPP duality, applied on `Bl_e LM_5`. By Castravet–Tevelev that
    variety is a birational contraction of `M̄_{0,7}`, so WTB is a statement about the moving cone of `M̄_{0,7}`, a
    classical object the write-up does not mention.
  * The numerical slab is correct, and I extend it to an 11-dimensional neighbourhood. The proposed completion route
    is provably blocked at the even-set vertex for `N = 1` (Proposition K).
  * New as book-keeping for the record; not new mathematics.
* **(b) Integrality at the root arrangement adds nothing beyond sections.** Proved only for single Δ-pole functions
  (mediant lemma, trivial once `Δ_ab` has ratio 4). The fixed-`S` Subspace part is unproved, and its stated principle
  is false in general (Defect 3).
* **(c) The new picture of `W` as cuspidal degenerations of a dominant four-branch family.** Confirmed numerically and
  independently (§3.3), and it is a nice structural picture. It is a statement about where the section cap is
  attained, i.e. inside the mapped barrier.
* **(d) Quadric cones at `α = 7/2`.** New and correct; incomplete without the cubics `ι(𝒬_σ)`. Calibration of `Z`,
  not a mechanism.
* **(e) The `SL_2` structure of triples.** The identity is correct and trivial. Consequences (1)–(2) were not
  re-derived beyond the identity. They have no bearing on `k = 5`, as the author says.
* **Phase 0.** The known paradigm (Vojta/Büchi, forum reduction); labelled "rediscovered" by the author.

No item is a new mechanism that beats `15/4` anywhere, and nothing transfers to `T_5`. The write-up says the same.

## 9. Corrected statement

> **Round 27, lane split (corrected).** There is no proof of split Question B, of Erdős #887 or of any growth bound.
> `M ≤ (2/3+ε) log R/loglog R` stands. The proved content is the following.
>
> 1. **(Thm 3.1, exact.)** The mnn31 `C = 7/2` five-divisor curve has `B·C = 7` (seven cuts, each once) and two
>    conjugate real-quadratic branches of order 1, so `α = 7/2`. Its 3-parameter conic family is dense in the quadric
>    cone `(u_2−u_1)(u_3−u_1) = (u_4−u_1)(u_5−u_1)`. Hence every admissible split `Z` in `A^split_5(δ)`, `δ < 1/14`,
>    contains the 15 cones `𝒬_σ` **and their 15 inversion images `ι(𝒬_σ)`** (cubic threefolds), besides `W` and
>    `ι(W)`. In the twisted `T_5` only `𝒬 ∩ σ𝒬` (a surface) can carry rational points.
> 2. **(Cost identity, Thm 3.2, proved.)** The forced-divisibility cost `Σ_B x_B π_B(S)` of a centrally symmetric
>    support is the same at every near-tight split profile (`Π_nt` = the 11-dimensional polytope of cut-balanced curve
>    classes). By Theorem B it is `≥ 4·reg(S)` (with `cost_{1/2} = (1/4)Σ_P w_P`).
> 3. **(Prop 3.7 = BDPP, proved.)** WTB at `x` holds iff `π^*γ_x − 4ℓ` is movable on `Bl_e LM_5`, a birational
>    contraction of `M̄_{0,7}` (Castravet–Tevelev). A dominating family of curves of class `π^*γ_{Nx} − 4Nℓ`
>    suffices.
> 4. **(Numerical, reproduced independently.)** Dominant 4-dimensional families of rational curves with four smooth
>    branches at `e` exist at 14 of the 16 one-sided classes, including the `W`-class (60-digit refinement). So WTB holds
>    numerically on an **11-dimensional** neighbourhood of `x = 1/2` in `Π_nt`, not only on a 10-dimensional slab.
> 5. **(Proposition K, proved here.)** At the `S_5`-invariant even-set vertex `x_{[5]}` (and, by `ι`, the odd one),
>    no rational curve of class `γ_{x_{[5]}}` has `mult_e ≥ 4`. The `N = 1` curve route to WTB is closed there, and the
>    "expected dimension `3N+1`" heuristic fails. WTB at `x_{[5]}` remains open. It needs `N ≥ 2`, positive genus, or a
>    divisor-side argument; otherwise a section beats `15/4` near the even profile, a split-only restricted statement
>    that cannot port to `T_5`.
> 6. **(Lemma 3.5, proved.)** Sections with poles along the root arrangement, paid for by bounded pair cofactors, gain
>    nothing where WTB holds (each `Δ_ab` has ratio 4 on `Π_nt`). Whether a genuine fixed-`S` Subspace argument on
>    `LM_5 ∖ ⋃Δ̄_ab` gains anything is **open**. The write-up's claim that it reduces to Lemma 3.5 is withdrawn.
> 7. **(Thm 3.6.)** The exact identity `x^τ × y^τ = (t_bc, t_ca, t_ab)` holds.
> 8. **(Phase 0, derived.)** The forum reduction gives complete intersections of quadrics (K3 for 4 divisors, general
>    type for `≥ 5`; this is the Vojta–Büchi paradigm). Both known families are Pell conics with
>    `k = (0,224,225,324,900)` and `(0,3,4,9,12)`.
>
> Nothing ports to `T_5`.

## 10. Scripts (all in `kill_split_checks/`, own code; heavy runs under `prlimit`)

| script | what it checks | label |
|---|---|---|
| `k_thm31.py` | conic/cone identities (symbolic), valuation vectors at all places for two parameter points (`B·C = 7`, 7 distinct cuts, quadratic place → `e`), exact Jacobian rank 3 at 4 points | symbolic / exact |
| `k_thm31b.py` | 58 mnn31 Pell points: divisibility, cone no. 0, `log(dist^{7/2}H)` periodic, `C → 2√3`; co-divisor tuples off all cones; m-czech points and generic W-curves on no cone; `ι(cone)` factorisation | exact |
| `k_profiles.py` | `dim Π_nt = 11`, 278 one-sided vertices / 16 classes, code hull dim 10, with `x_{[5]}` dim 11 | exact |
| `k_misc.py` | exact `V(y)` formula (200 tuples), Thm 3.6 identity (200 tuples × 10 triples), `Σx_B f(|B|) = 40`, m-czech phase-0 `k`-tuple | exact |
| `k_cover.py` → `k_cover_out_1.json` | four-branch covering families at all 16 one-sided classes: found at 14, dominance rank 4 | numerical |
| `k_cover2.py` → `k_cover2_c13.out`, `k_cover2_c4.out` | 600 s each at the `x_{[5]}` pair: only degenerate convergence | numerical |
| `k_refine.py` | 60-digit Newton refinement and dominance SVD at the `W`-class | numerical (high precision) |
| `k_hull14.py` | the 276 vertices of the 14 classes span `Π_nt` (dim 11), average `1/2` | exact |
| `k_even_thm.py` | Proposition K ingredients on random exact instances; four-concurrent-lines converse | exact |
| `k_even_reg.py` → `k_even_reg.out` | `reg` vs `cost_{x_{[5]}}` for symmetric supports `P(α,β)`: max ratio 1 | exact mod two primes |
| `k_coverN.py`, `k_coverN2.py`, `k_coverN2b.py` → `k_N2*_*.out` | `N = 2` searches (control failed, so inconclusive) | numerical, negative |

Working notes: `../kill_split_notes.md` (N1–N17).
