# Referee report on round24/theory.md (anisotropy candidates reduce to split caps)

Status: complete. Adversarial referee of `round24/theory.md`. Verdict: **survives with corrections** (2 major gaps repaired, 6 minor issues).

## 0. Verdict

**The result survives, with corrections; status `reduces_to_split_cap` is confirmed.** I found no error that changes
the conclusion: within the class `𝓛` that theory.md defines, no argument beats the split caps (`15/4` on balanced simple
profiles, `δ = 1/10` for arbitrary profiles). `M ≤ (2/3+ε) log R/loglog R` stands and nothing here bears on growth or
uniformity. Every number in theory.md was recomputed with independent code and agrees. I also finished the two
computations theory.md left open: the isotypic orders for `m = 6, 7`, and index additivity on `V_1 ⊗ V_2`, for which
I give a three-line general proof, so no computation is needed.

I found **two proof gaps**, both repaired below (§2.5, §2.6), and **several overstatements**. None of them changes the verdict:

* **(M1, major, repaired) Group words in Theorem 7.1(b).** Item 2 of Theorem 7.1 is proved only "with pairwise disjoint
  primes". A tuple that contains group words in A-points, such as `(y_1, y_1y_2)` or `(y, ιy)`, never has
  pairwise disjoint primes. The cited Prop 4.3 therefore does not cover the group words that the theorem and the
  summary ("the whole class 𝓛 ... group words ... is capped") claim. Repair: Lemma R1 (pullback along the
  homomorphism to the generating points), §2.6.
* **(M2, major, repaired; touches an owner constraint) Prop 6.2, non-homogeneous `G`.** The step "otherwise fix the
  representative with `P_1 = 1`" covers only arguments that evaluate `G` at that representative. A `G` that is not
  homogeneous in the common factor is not a function on `T_5`. Evaluated at any other representative, such as one tied to
  the position of the arc, it depends on a drifting common phase. Replacing that representative by `P_1 = 1` is
  exactly "replacing a drifting common phase by zero". The identity `G = Π P̄_j^{deg_j} g(x)` also needs `G`
  to be multihomogeneous in each pair `(U_j, V_j)`, which the statement does not say. Repair: decompose into homogeneous
  components (Lemma R2, §2.5). The conclusion of Cor 6.4 is unchanged, but "identical" becomes "implied by".
* **(m1-m6, minor)** These are: uniformity in `F` of the `O_F(1)` terms (Prop 4.3(iii) and Thm 7.1 need
  `log deg F = o(w)` and `log H(F) = o(deg F · w)`); a list of anisotropy-sensitive inputs that claims to be exhaustive
  but misses congruence divisibility at the inert and ramified places; the prose of Prop 5.3 ("weaker than or equal
  to the split ones" holds only against the better of `f` and `ι^*f`); "simple" is missing from Cor 3.3 and from the
  summary of Prop 4.3; "distance `= |det|/√(Φ_0Φ_0)`" holds only up to a factor 2 (twisted) or `1+O(dist)` (split); and
  the phrase "every argument equals its average".

The corrected statement is in §5.

## 1. What was checked and how

* I read theory.md §§0-8 in full, together with the cited record: expert_note §3, round22 `toric.md` §1 and §3
  (Lemma 1.3, Lemma 3.1, Thms 3.2-3.3), round23 `width.md` Thm 2.5 and README §1.1 (`split.md`, parity lemma),
  and `gauss.md` §1.2 (`t_ab`).
* I re-derived every proof in §§2-7, with all quantifiers. The derivations are in `round24/referee_notes.md`.
* I wrote new scripts from scratch, without reusing the author's code, in `round24/referee_theory_checks/`:
  * `iso_ref.py` → `iso_ref.out`, `iso_ref67.out`: isotypic maximal orders. The rows are split by parity, and ranks are
    computed mod `2^61-1` with python-flint. This certifies the upper bounds, since full rank mod p implies full
    rank over Q, for `m = 1..7`.
  * `sec_ref.py` → `sec_ref.out`: the building blocks `V, V_abc, D_ab` and the table of Prop 3.1(b). Support,
    `ι`-parity and order at `e` are computed directly from `Σ c_s⟨s,y⟩^j`, the orders for `m ≤ 4` without assuming additivity.
  * `h90_ref.py` → `h90_ref.out`: Prop 2.2(i)-(iv) and Prop 6.1. These are exact Gaussian-integer computations on 304
    random primitive 5-tuples and 3040 pairs, with an independent construction of `P_j` (primitive part of `N(1+x_j)`).
  * `misc_ref.py` → `misc_ref.out`: the Picard ranks, the identity of Prop 5.3 together with a count of the cases it
    contradicts (prose), and the width half of Lemma R1.
* Resource use: all runs were under `prlimit --as=3500000000`, with at most one heavy process at a time.

## 2. Section-by-section audit
### 2.1 §2 Galois structure

* **Lemma 2.1: correct.** (a) `v_{π̄}(x_j) = v_π(x̄_j) = -v_π(x_j)`. (b) A conjugation-stable place `w` has
  `v_w = -v_w`. (c) `x x̄ = 1`. The descent convention `σ_X = ι∘c` agrees with toric.md §1.1 (`σ` acts by `-1` on
  `N` and `M`, and `1_A ↦ 1_{A^c}`).
* **The "two totally real slices" paragraph: correct.** A homogeneous polynomial that vanishes on `R^4` or on `iR^4`
  is zero, so orders at `e` can be read on either slice.
* **Prop 2.2: correct, with one imprecision in the "Reading" paragraph (m5).**
  * Proofs re-derived:
    * (i) Hilbert 90 per coordinate. `P/P̄ = P'/P̄'` iff `P'/P ∈ Q^*`, so a primitive `P` is unique up to sign.
    * (ii) `P_aP̄_b - P̄_aP_b = 2i[P_b,P_a]`.
    * (iii) A primitive tuple has `min_j a_j(p) = 0` and `max_j a_j(p) = k_p`, so `Π p^{range ε} = R''^2`.
    * (iv) The case check gives `c_ab - (|ε_a|+|ε_b|)/2 = -|ε_a-ε_b|/2` in both sign cases.
  * Exact check (`h90_ref.out`): 304 tuples and 3040 pairs, zero failures for all of the following:
    * `P_j/P̄_j = ζ_j/ζ_1`, with `P_j` built independently as the primitive part of `ζ_jζ̄_1 + N`;
    * `N(P_j) = 2^{e_j}Π p^{|ε_j|}`, `ε_j - ε_1 = a_j - a_1`, and `Π p^{range ε} = R''^2`;
    * `Π p^{c_ab} | [P_a,P_b] ≠ 0`;
    * `t_ab^2` computed from the record's definition `|ζ_a-ζ_b| = 2|t_ab|R''e^{-d_ab/2}` equals
      `[P_a,P_b]^2/(Π p^{2c_ab}2^{e_a+e_b})`;
    * `t_ab ∈ Z` when `e_a = e_b`, and `2t_ab^2` is an odd integer when `e_a ≠ e_b`.
  * New, exact: **`e_a = e_b` iff the units of `ζ_a, ζ_b` agree up to sign** (`u_a/u_b = ±1`). So the record's
    integrality of `t_ab` within one unit class extends to unit classes modulo `±1`. This is consistent with the
    record, which only claims sufficiency.
  * Imprecision (m5): the Reading says "distance `= |det|/√(Φ_0Φ_0)`". For the twisted form the exact identity has
    a factor 2 (`2|[P_a,P_b]|/(|P_a||P_b|)`). For the split form the exact identity is
    `|u_a/u_b - 1| = |det|/|a_b b_a|`, and `√|Φ_0(P_a)Φ_0(P_b)| = |a_b b_a|·(u_a/u_b)^{1/2}`, so the two agree only
    up to `1 + O(dist)`. The statement should read "`≍`". Nothing downstream uses the equality.
* **Lemma 2.3: ranks correct** (`misc_ref.out`: 30 rays, 15 orbits under `A ↦ A^c`, `M^Gal = 0`, so ranks 15 and
  26). `Pic X = (Pic X̄)^G` holds exactly because `X(Q) ≠ ∅`. The counting statement is correctly labelled heuristic.
### 2.2 §3 Candidate 1 (isotypic sections)

* **Prop 3.1: correct and independently confirmed.**
  * (a): `ι^*` acts on the degree-`j` Taylor part by `(-1)^j`. `σ_X` multiplies it by `(-1)^j` after conjugation.
    `ι` and `σ_X` commute (`σ ι^* f = ι^* σ f = f̄`), so Speiser's lemma (descent for a semilinear involution
    over `Q(i)/Q`) applies inside `V_m^±`.
  * (b): the upper bounds follow from width Thm 2.5 (refereed) and the parity of (a). The lower bounds follow from the
    table. I checked the parities of `⌊10m/3⌋` for `m = 3s, 3s+1, 3s+2` (even, odd, even).
  * `sec_ref.out`, computed directly:
    * `V`: least `m` = 3, even, order 10. `V_abc`: least `m` = 1, odd, order 3. `D_ab`: least `m` = 1, even, order 2.
    * The table entries for `m = 1..7` lie in `S_m` with the stated parity.
    * Their orders for `m ≤ 4` are `(2,3), (6,5), (10,9), (12,13)`, computed from `Σ c_s⟨s,y⟩^j` and not from additivity.
  * `iso_ref.out` and `iso_ref67.out` (new rank computation, rows split by parity): for every `m ≤ 7` the claimed
    isotypic order `t` has full rank at degree `≤ t`, which certifies the upper bound, and a deficient rank at degree
    `≤ t-1`, as the explicit sections require:
    `(t_+,t_-) = (2,3), (6,5), (10,9), (12,13), (16,15), (20,19), (22,23)`.
    The `m = 6` value, which the author stopped for time, and `m = 7` are new. They agree with (b).
* **Prop 3.2: correct.**
  * (a) Over `K`, `2log|α| = -Σ_{w∤∞} log|α|_w`. The bounds at `π, π̄` give `-width·log p`. Elsewhere the bound is `≤ 0`.
  * (b) `f^*(y) = \overline{f(y)}` on norm-one points. The inequality for `ff^*` is twice that for `f`.
  * (c) A `Q`-rational `f` has `c_{-s} = c̄_s`, so its support is automatically centrally symmetric. The "reality"
    bound `max(min⟨s,v⟩, -max⟨s,v⟩)` then equals `-½width`. Reality adds nothing. I confirm this.
  * Remark (not an error): the word "exactly" in (a) means the product-formula form with `width_{Q_f}(v_π)` at each
    prime. On primes with several layers this form is *sharper* than toric Thm 3.2's cut-sum form. The support function
    is only subadditive along a chain `v = Σ_t 1_{A_t}`, with equality on every chain when `Q_f` is a generalised permutohedron (normal fan coarser than the braid fan) and possibly strict otherwise.
    The caps (MW), (τ) are proved against the cut-sum form, so they hold only on **simple** profiles. Theorem 7.1
    says this, but Cor 3.3 ("on `η`-balanced profiles the best is `μ ≤ 4+O(η)`") omits "simple" (m4). For the barrier
    reading it is harmless: simple balanced configurations already have to be excluded by any proof of `(L_5^∘)`.
* **Cor 3.3: correct, with "simple" added.** `δ = 2/τ_e - 1/2 = 3/5 - 1/2 = 1/10`.
### 2.3 §4 Candidate 2 (conjugate products, Thue)

* **Prop 4.1: correct.**
  * `x^s⊗x^{s'} ↦ x^{s-s'}` and `S_{d_1} - S_{d_2} ⊆ S_{d_1+d_2}`.
  * Each Taylor part of `g` of degree `j` comes from parts `F_J` with `j_1+j_2 = j`.
  * At `π` the lower bound for `F` uses the valuation vector `(v,-v)`, and `Q_g ⊆ {s-s'}`.
  * Hence `t_F λ ≤ ord_e(g) λ ≤ RHS_g ≤ RHS_F`. This uses `λ ≥ 0`, which holds for A-points and should be stated.
* **Prop 4.2: correct.** The lex-minimal proof is valid. A shorter proof, which I used in place of the unfinished
  `V_1⊗V_2` computation:
  * Let `W_i^{≤k} := span{(s^α)_{s∈S_i} : |α| ≤ k} ⊆ Q^{S_i}`. By Lemma 3.1, `W_i^{≤R_i} = Q^{S_i}`.
  * "Weighted index `> Σ_i w_iR_i`" forces `F_J = 0` for every `J ≤ (R_1,…,R_n)` componentwise. So the
    coefficient tensor is orthogonal to `⊗_i W_i^{≤R_i} = ⊗_i Q^{S_i}`, and hence is zero.
  * The product of extremal sections attains `Σ_i w_iR_i`.
  * So the exact checks on `V_1⊗V_1` (6 and 9) and the expected `V_1⊗V_2` value 9 are special cases of a
    general statement. No computation is needed.
* **Prop 4.3: (i) and (ii) correct. (iii) correct up to the uniformity qualifier m1.**
  * (iii) needs simple profiles, disjoint primes and `λ_i ≥ 0`.
  * (MW) applies to `S_i = Q_i ∩ M` because `conv(S_i) = Q_i`: the projection of a lattice polytope has lattice vertices.
  * (m1) "`(*)` holds whenever every `μ(y_i) ≤ 4(1-η)`" ignores the `O_F(1)`, which can have either sign. It
    contains `log` of the Taylor constant, which is `≥ -t log t` because the leading Taylor coefficients of an integral
    `F` have denominators dividing `t!`. In a Thue–Roth argument, `F` and its degrees are chosen after the points, so a
    threshold `w_0(F)` is template-dependent. The correct statement is the following.
  * **Corrected (iii):** `(*)` holds at the configuration whenever `μ(y_i) ≤ 4(1-η) - C(F)/w_min` for all `i`, where
    `C(F) ≪ log H(F) + Σ_i d_i log(2+d_i)`. So the cap applies uniformly to every family with
    `log H(F) + Σ d_i log d_i = o(Σ_i d_i w_i)`. Families from Siegel's lemma whose degrees are at most polynomial in
    `e^{w}` satisfy this. Without such a qualifier, "holds" is not proved.
* **Remark (1)** (shared primes, open point 1): the logic is correct. A proof of `A_5` must also exclude hypothetical
  disjoint-prime configurations, so the barrier does not need the shared-prime cap. That remark covers only
  *configurations of A-points*, though, not *group words* in them (M1, §2.6).
* **Cor 4.4: correct**, with m1. The arbitrary-profile computation `R_i ≤ (5/3)·max_A width_{Q_i}(1_A)` follows
  from `ff^*` and Thm 2.5. It gives `δ = 1/10`. The witnessing profile depends on `F` and is labelled as such.
### 2.4 §5 Candidate 3 (point parity)

* **Prop 5.1: correct.** The finite local data of a twisted point are `(v_π, -v_π, 0, …)`, read off from the Gaussian
  profile.
* **Cor 5.2: correctly labelled record-conditional.** It applies only to inequalities proved from local data **and
  monomial gaps**, which is what the fake statement covers. The summary sentence "fakes satisfy it" refers to the
  parity itself and is true.
* **Prop 5.3: the identity is correct, the prose is not (m3).**
  * Identity: `r_tw = ½(r_sp(v) + r_sp(-v))` holds on 20000 random `(Q,v)` (`misc_ref.out`), and (a), (b) are correct.
  * The sentence "parity can only make the twisted Liouville inequalities weaker than or equal to the split ones at the
    same cut heights" is false section by section. For `f = x_1/x_2 - 1` and `v = 1_{\{2\}}`, `r_sp(v) = 1` but
    `r_tw = ½`, so the twisted inequality for this `f` is *stronger*. Twisted is strictly smaller in 8285 of 20000
    random cases.
  * The true statement compares with the **better of `f` and `ι^*f`**. Since `r_sp^{ι^*f}(v) = r_sp^f(-v)`, the
    best split bound is `min(r_sp(v), r_sp(-v)) ≤ r_tw`. This is what (b) states for the pair (`min(X,Y)` against
    `(X+Y)/2`) and what the barrier needs.
  * The slogan "the anisotropy changes which points exist, not which section inequalities hold" should read
    "...not the caps of section inequalities". The two families of inequalities differ (symmetrised against
    one-sided local terms); only their caps coincide.
* **The curve paragraph and Cor 5.4: correct.** `B·C ∈ 30Z`, contact `≥ 8`. The concentration step is correctly called open.
### 2.5 §6 Candidate 4 (Hilbert 90)

* **Prop 6.1: correct** (exact, `h90_ref.out`; see §2.1).
* **Prop 6.2: the statement is incomplete and one proof step is wrong (M2). The conclusion holds after repair.**
  * The identity `G(P,P̄) = Π_j P̄_j^{deg_j} g(x)` requires every monomial `U^cV^d` of `G` to have `c_j + d_j = deg_j`
    for each `j`, i.e. multihomogeneity in each pair. "Degree `deg_j` in `(U_j,V_j)` jointly", read as total degree
    `≤ deg_j`, is not enough.
  * If `G` is not homogeneous in the common factor (`Σ_j c_j` not constant), then `G(P,P̄)` depends on the representative
    `x = ρx^{(1)}` of `y` (`ρ ∈ K`, `N(ρ) = 1`; `x^{(1)}` is the `P_1 = 1` representative): `g(ρx^{(1)}) = Σ_k ρ^k g_k(x^{(1)})`.
    Fixing `P_1 = 1` handles only the argument that evaluates `G` at `x^{(1)}`. An argument whose `G` is evaluated at the
    lattice points as they sit on the arc carries the phase `ρ`, and that phase drifts with `R`. Saying "fix `P_1 = 1`"
    then replaces a drifting common phase by zero, which the owner's constraint forbids. The constraint audit of §8
    asserts the opposite.
  * Example: `G = U_1V_2 + U_1U_2` gives `g = x_1 + x_1x_2`. At `x^{(1)}` this is `1 + x_2`; at `ρx^{(1)}` it is
    `ρ + ρ^2x_2`. These are different functions of `y`.
* **Lemma R2 (repair; proved).** Let `0 ≠ G ∈ Z[i][U,V]` be arbitrary, evaluated at the primitive `P_j` of any
  representative `x = ρx^{(1)}` of `y ∈ T_5(Q)`, with `G(P,P̄) ≠ 0`.
  * *Decomposition.* Write `G = Σ_D G_D`, with `G_D` multihomogeneous of degree `D_j` in `(U_j,V_j)`. Write
    `G_D(P,P̄) = Π_j P̄_j^{D_j} Σ_k g_{D,k}(x)`, splitting by total degree `k`. Put
    `h_{D,k} := x_1^{-k} g_{D,k} ∈ K[M]`, a section. Then `G(P,P̄) = Σ_{D,k} Π P̄_j^{D_j} ρ^k (x^{(1)}_1)^k h_{D,k}(y)`.
  * *Statement.* Any inequality of type (c) in `𝓛` for `G` at `y` is implied, up to `O_G(1)`, by (TL) for one
    section `h_{D*,k*}` with `h_{D*,k*}(y) ≠ 0`.
  * *Proof.*
    1. Upper bound: `|G(P,P̄)| ≤ Σ_{(D,k) ∈ 𝒩} |P|^D|h_{D,k}(y)| ≤ C max_{𝒩} |P|^D dist^{t_{D,k}}`, where `𝒩` is
       the set of components with `h_{D,k}(y) ≠ 0` (nonempty) and `t_{D,k} = ord_e h_{D,k}`. This uses `|ρ| = |x_1| = 1`.
       A coarser upper bound, such as the order of `G` along `{x_j}` equal, which is `≤ min t_{D,k}`, only weakens the
       inequality. Let `(D*,k*)` attain the maximum.
    2. Lower bound: forced divisibility at `π, π̄` for `G` uses the minimum over **all** monomials of `G`. This is at
       most the minimum over the component `(D*,k*)`.
    3. On that component, `Σ_j (c_jα_j + d_jβ_j) = Σ_j D*_jβ_j + ⟨c, ε^{(1)}⟩ + k*·v_π(ρ)`, where the last term is
       constant on the slice `Σc = k*`. So the `π` and `π̄` bounds add to `Σ_j D*_j(α_j+β_j) - width_{Q_h}(v_π(y))`.
    4. With `log|P|^{D*} = ½Σ_pΣ_j D*_j(α_j+β_j) log p + O(1)`, the inequality for `G` reads
       `t*λ ≤ RHS_G` with `RHS_G ≥ ½Σ_p width_{Q_{h_{D*,k*}}}(v_π(y)) log p + O(1) = RHS_{TL}(h_{D*,k*})`.
    5. The latter is implied by (TL) for `h_{D*,k*}` at `y`, which is nonzero there. `□`
  * Consequence: Cor 6.4 holds as stated, with "identical" (Prop 6.2) replaced by "implied by" for non-homogeneous `G`.
    The no-half-integer remark is unaffected: every `h_{D,k}` has integer exponents.
* **Prop 6.3: correct.** `G_{d,c} = ±\overline{G_{c,d}}` makes the support symmetric about `deg/2`. The same holds
  per component after Lemma R2.
* **Plücker/Ptolemy paragraph:** correct as a pointer to the record's abc barrier (refereed elsewhere). Not re-derived here.
### 2.6 §7 Theorem 7.1

* **Item 1: correct for (a) and (c), the latter after Lemma R2.** For (b) it is correct as stated: an `n`-point
  inequality at the tuple of words.
* **Item 2 for (b) with group words: the proof does not apply (M1).**
  * The proof cites Prop 4.3 ("group words are points of `T_5(Q)`"). But Prop 4.3(iii) needs the **entries** to have
    pairwise disjoint primes. A tuple containing a word together with one of its letters, such as
    `(y_1, y_1y_2)`, `(y, y^2)` or `(y, ιy)`, never does.
  * The "(with pairwise disjoint primes in case (b))" clause therefore makes item 2 vacuous for exactly the
    group-word arguments that the summary lists as capped. Prop 4.1 rescues only the special word `ιy` at one point.
* **Lemma R1 (repair; proved).** Let `y_1,…,y_r ∈ T_5(Q)`, let `φ: T^r → T^n` be the homomorphism with integer
  matrix `(n_{ki})`, `w_k = Π_i y_i^{n_{ki}}`, and `F ∈ Z[i][M^n]` with `F(w) ≠ 0`. Put `G := φ^*F ∈ Z[i][M^r]` and
  `λ'_k := min{λ(y_i) : n_{ki} ≠ 0}`. Then:
  * (a) `G(y) = F(w) ≠ 0`.
  * (b) `RHS_F(w) ≥ RHS_G(y)` in `(*)`.
  * (c) `LHS_G(y) ≥ min_{F_J ≠ 0} Σ_k j_kλ'_k`.
  * (d) If `λ(w_k) ≤ λ'_k + c_0` for all `k`, then `(*)` for `G` at `y` implies `(*)` for `F` at `w`, up to `c_0·deg F`.
  * *Proof.*
    * (b) `v_π(w) = φ_*(v_π(y))` and `⟨s, φ_*v⟩ = ⟨φ^*s, v⟩`, so `width_{Q_F}(φ_*v) = width_{φ^*Q_F}(v)`. Also
      `Q_G ⊆ φ^*Q_F`, since cancellation only removes points. Checked on 5000 random cases (`misc_ref.out`).
    * (c) `G(exp Y) = Σ_J F_J(LY)` with `(LY)_k = Σ_i n_{ki}Y_i`. Every monomial of `F_J(LY)` has multidegree `a`
      with `Σ_i a_iλ_i ≥ Σ_k j_kλ'_k`, and every nonzero Taylor term of `G` comes from some nonzero `F_J`.
    * (d) `LHS_F ≤ min_J Σ_k j_k(λ'_k + c_0) ≤ LHS_G + c_0 deg F`. `□`
  * *Use for the barrier.* Take a hypothetical configuration of `η`-balanced simple A-points `y_i` with pairwise
    disjoint primes and lacunary scales, `λ_{i+1} ≥ λ_i + log(2Σ_{k,i}|n_{ki}|) + 1`. Then
    `|Σ_i n_{ki}Y_i| ≥ ½|Y_{i_0}|` (with `Y = log y ∈ iR^4`, `i_0` the letter of smallest `λ`), so
    `λ(w_k) ≤ λ'_k + O(1)`. By (d) and Prop 4.3(iii) applied to `G` at `y`, every `(*)` for every `F` at every finite
    tuple of words holds once `μ(y_i) ≤ 4(1-η) - o(1)`. Item 2 for group words then holds with "pairwise disjoint
    primes **of the generating A-points**". Lacunarity is a choice of hypothetical configuration, which is enough for
    a barrier; it is not a claim about actual A-points.
* **Item 3 (split): correct.** On split-balanced profiles that put half of each cut's height on each side, the split
  right side equals the twisted one, so (MW) gives the same cap. On one-sided profiles, choosing the side with the
  larger `a_A(Q_f)` gives `≥ ½Σwidth`, again by minimax over a convex set.
* **"Mechanism" paragraph:** "every argument equals its average" (m6) is exact for (a) (Prop 3.2(b)) and for the
  pair `(y, ιy)`. For (c) and for words it is "is implied by its average" (Lemmas R1, R2).
* **"What is not in `𝓛`" list: not exhaustive (m2).** It omits **congruence (residue-field) divisibility at the
  inert places `q ≡ 3 (4)` and at `(1+i)`**.
  * Twisted points always reduce into the finite group `T_5(F_q) = μ_{q+1}^5/μ_{q+1}`; split points at a prime may have
    contact there instead. So a section that vanishes mod `q^c` on `T_5(Z/q^c)` gives `v_q(f(y)) ≥ c` for every
    twisted point. This is twist-sensitive (the group is `μ_{q+1}`, not `μ_{ℓ-1}`), and it strengthens (TL) by
    `Σ_q c_q(f) log q`.
  * *Proved bound:* for `q`-primitive `f ∈ Z[i][S_m]` and `q ≥ 2m`, `c_q(f) = 0`. Exponents of `S_m` have `|s_j| ≤ m`, so
    distinct `s ∈ S_m` give distinct characters of `μ_{q+1}^5/μ_{q+1}`. Hence `f mod q ≠ 0` on `T_5(F_q)`. `T_5` is
    `Q`-rational, so weak approximation lifts a non-vanishing residue class to `T_5(Q)`.
  * So for a fixed finite family the gain is `O_F(1)`. The bound `Σ_q c_q log q = O(m log m)` for growing `m`, by the
    integer-valued-polynomial (`v_q((Cm)!)`) argument, is a **sketch** and is not needed for the verdict.
  * This is the same uniformity issue as m1: it is harmless when `log m = o(w)`. The list should name it.
### 2.7 §8 Constraint audit

* *Fixed-`S` Subspace:* none is used. Confirmed.
* *Gaussian gcd and denominators:* (TL) comes from the full product formula over `K`. Tuples are primitive, and
  `2^{e_j}` is kept. Confirmed. The omitted inert-place divisibility (m2) is a gain that was left out, not a cost
  that was discarded, so it does not violate this constraint.
* *Drifting phase:* **not fully respected** in Prop 6.2 (M2). Lemma R2 removes the dependence on the representative.
* *Template-dependent thresholds:* respected for the witnessing profiles. **Not respected** for the `O_F(1)`
  thresholds in Prop 4.3(iii) and Thm 7.1 item 2 (m1). Corrected in §5.
* *Evidence against proof:* labels are accurate. Lemma 2.3's count is heuristic, Cor 5.2 is record-conditional, and
  the exact checks are used only as checks. Lemma K1 is cited as killer-proved, which is the record's own label.

## 3. Independent computations

All computations are exact or certified mod a prime in the valid direction. Scripts are in `round24/referee_theory_checks/`.

| claim | script / output | result |
|---|---|---|
| Prop 3.1(c): isotypic orders, `m ≤ 5` | `iso_ref.py` / `iso_ref.out` | `(2,3),(6,5),(10,9),(12,13),(16,15)`; full rank at the claimed degree, deficient one below |
| Prop 3.1(b) at `m = 6, 7` (author: not finished) | `iso_ref.py 6 7` / `iso_ref67.out` | `(20,19)`, `(22,23)`, as predicted |
| Prop 3.1(b) building blocks and table | `sec_ref.py` / `sec_ref.out` | orders 10, 3, 2; table orders for `m ≤ 4` computed directly; supports and parities for `m ≤ 7` |
| Prop 2.2(i)-(iv), Prop 6.1 | `h90_ref.py` / `h90_ref.out` | 304 tuples, 3040 pairs, 0 failures; new: `e_a = e_b` iff `u_a/u_b = ±1` |
| Lemma 2.3 ranks | `misc_ref.py` / `misc_ref.out` | 15 orbits, so ranks 15 and 26 |
| Prop 5.3 identity; prose | `misc_ref.out` | identity holds on 20000 cases; twisted RHS `<` split RHS at `v` in 8285 of them (m3) |
| Prop 4.2 (`V_1⊗V_2`, unfinished) | proof in §2.3 | general three-line proof; no computation needed |
| Lemma R1(b) | `misc_ref.out` | 0 failures in 5000 random cases |

## 4. Issues list

| id | severity | place | issue | status |
|---|---|---|---|---|
| M1 | major | Thm 7.1 item 2, case (b); summary item 6 | group-word tuples never have pairwise disjoint primes, so Prop 4.3 does not apply | repaired by Lemma R1 (§2.6); conclusion unchanged |
| M2 | major (constraint) | Prop 6.2, §8 audit | non-homogeneous `G`: "fix `P_1 = 1`" sets a drifting phase to zero; per-pair multihomogeneity is not stated | repaired by Lemma R2 (§2.5); "identical" becomes "implied by" |
| m1 | minor | Prop 4.3(iii), Cor 4.4, Thm 7.1(2) | `O_F(1)` thresholds are template-dependent; "holds whenever `μ ≤ 4(1-η)`" needs `w ≥ w_0(F)`, i.e. `log H(F) + Σd_i log d_i = o(Σd_iw_i)` | corrected in §5 |
| m2 | minor | §7 "not in `𝓛`" list | omits congruence divisibility at inert places and at `(1+i)` (twist-sensitive: `T_5(F_q) = μ_{q+1}^4`); effect `O_F(1)` (proved: `c_q = 0` for `q ≥ 2m`) | add to the list |
| m3 | minor | Prop 5.3 prose; §7 slogan | "weaker than or equal to the split ones" is false section by section; true against `min(f, ι^*f)` | reword |
| m4 | minor | Cor 3.3; summary of Prop 4.3 | "simple" is missing: the caps are proved for the cut-sum form, which equals the product-formula form only on simple profiles | add "simple" |
| m5 | minor | Prop 2.2 Reading | "distance `= |det|/√(Φ_0Φ_0)`" holds only up to a factor 2 (twisted) or `1+O(dist)` (split) | replace by `≍` |
| m6 | minor | §7 mechanism | "every argument equals its average" is exact only for (a) and `(y, ιy)`; otherwise "is implied by" | reword |

Also noted (not an error): `λ ≥ 0` should be stated in Props 4.1 and 4.3, and Lemma K1 remains killer-proved,
not independently refereed.

## 5. Corrected statement

**Theorem (round 24 theory, referee-corrected; proved except where labelled).** Let `𝓛` be the class of arguments
that derive, at points where the auxiliary object does not vanish, "archimedean Taylor upper bound `≥` product-formula
lower bound". The lower bound may use integrality, toric forced divisibility at every place of `K = Q(i)`, point
parity and reality of values. The auxiliary objects are of three kinds:
(a) `f ∈ Z[i][M]` of any `σ_X`/`ι`-isotypic type;
(b) `F ∈ Z[i][M^n]` evaluated at a tuple of group words (including `ι`) in A-points `y_1,…,y_r`;
(c) any `G ∈ Z[i][P,P̄]` evaluated at the Hilbert-90 vectors of any representative of `y`.
Then:
1. *(reduction)*
   * Each inequality from (a) is the (TL) inequality for the `Q`-rational symmetric section `ff^*` (exact).
   * Each inequality from (c) is implied by (TL) for one homogeneous component `h_{D,k} ∈ K[S_m]` with
     `h_{D,k}(y) ≠ 0` (Lemma R2).
   * Each inequality from (b) is implied, at configurations with `λ(w_k) ≤ λ'_k + O(1)`, by the `n`-point inequality
     `(*)` for `φ^*F` at the generating points (Lemma R1, Prop 4.1).
2. *(balanced cap)* Let `y_i` be `η`-balanced A-points with **simple** profiles and pairwise disjoint primes, with
   lacunary scales in case (b). Every inequality in `𝓛` holds at them whenever `μ(y_i) ≤ 4(1-η) - C(F)/w_min`, with
   `C(F) ≪ log H(F) + Σ_i d_i log(2+d_i)`. Hence no argument in `𝓛` whose auxiliary family satisfies
   `log H(F) + Σ_i d_i log d_i = o(Σ_i d_i w_i)` proves `(L_5^∘)` or `A_5(δ)` for any `δ < 1/30`. In particular, no
   argument with a fixed finite family does.
3. *(arbitrary profiles)* Under the same proviso, no argument in `𝓛` proves `A_5(δ)` for `δ < 1/10`. The witnessing
   profile depends on the family.
4. *(form-independence)* The same numbers hold for the split form. The isotypic maximal orders are `t_max(m) = ⌊10m/3⌋`
   for the part of the same parity and `t_max(m) - 1` for the other part, for every `m`. This is proved, and exact for
   `m ≤ 7`.
5. *(not covered)*
   * `(*)` on shared-prime configurations of genuine A-points (open; not needed for 2).
   * Arguments outside `𝓛`: global inputs (curve concentration plus the parity lemma, conditional on McKinnon/Vojta)
     and congruence divisibility at the inert and ramified places. The latter is twist-sensitive, but its effect is
     `O_F(1)` for fixed families (`c_q(f) = 0` for `q ≥ 2m`).

**Status: `reduces_to_split_cap` confirmed.** No growth improvement and no uniform bound:
`M ≤ (2/3+ε) log R/loglog R` stands. Question B is unchanged.
