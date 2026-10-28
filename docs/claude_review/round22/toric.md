# Round 22 — toric.md: Seshadri and approximation constants of the twisted Losev–Manin 4-fold at e

Status: complete draft, unrefereed (round 22, lens: toric Seshadri / approximation constants at the identity).

## 0. Summary and verdict

**Verdict: `restricted_theorem` + barrier.** No growth improvement and no uniform bound. `A_5(δ)` for `δ < 1/30` is
neither proved nor refuted. What is new is an exact placement of `A_5` / `(L_5^∘)` among the approximation constants of
`(X, B, e) = (L̄_5^tw, -K, identity)`, with proofs and exact computations:

* **(i) Seshadri constant (Thm 2.1, proved).** `ε_e(-K) = 2`, for the split and the twisted form alike (descent is
  trivial for `ε`; the extremal curves descend). The polytope of `-K` is the root polytope
  `R_5 = conv(e_a - e_b) = {Σx = 0, ‖x‖_1 ≤ 2}`, `-K = H + H'` (two blow-downs to `P^4`), `-K` nef, big, not ample,
  `(-K)^4 = 70` (Prop 1.1, exact). The minimum is attained by the 15 `Q`-conics `C_{1_A}` (closures of norm-one
  subtori), which all lie in the diagonals `Z_diag` (repeated points). So `α_e(-K) = ε_e(-K) = 2` on the whole variety,
  in accordance with McKinnon's picture, and **irrelevant to `A_5`**.
* **(ii) Approximation constants (§4).** Along 1-PS conics `α = 2(max v - min v)` (`≥ 8` off `Z_diag`); along a
  rational curve `α = (B·C)/μ_0`; on balanced curves off `Z_diag`, `α ≥ 15/4` with equality iff the curve is null for
  the pair-product divisor (Prop 2.3, Cor 4.2); for balanced sequences in general `α ≥ 15/4` (Liouville). The `k = 4`
  family of `gauss.md` Thm 4.1 is exactly a covering family of such null curves at `k = 4` (class `2γ_4`, `B·C = 14`,
  contact 4, `α = 7/2`; Thm 4.4, exact). At `k = 5` covering rational families through `e` have `α ≥ 4` (Prop 4.3),
  above the balanced Liouville constant `15/4`; the switch happens exactly between `k = 4` and `k = 5`
  (`k - 1` vs `4(1 - 2^{1-k})`).
* **New exact computation: the Liouville threshold `τ_e(-K)`.** Toric Liouville with any auxiliary polynomials, forced
  divisibility included, is governed by `τ_e(-K) = lim reg(mR_5 ∩ M)/m` (Lemma 3.1, Thm 3.3, proved: minimax over
  profiles). Proved `10/3 ≤ τ_e(-K) ≤ 15/4`; exact: **no section of `m(-K)`, `m ≤ 7`, vanishes at `e` to order
  `> ⌊10m/3⌋`** (Thm 3.4; certified by full rank mod a prime). So toric Liouville cannot pass `δ = 1/30` (proved), and on
  the data it does not even improve Cilleruelo–Córdoba's `δ > 1/10` off a closed set (the open "intermediate exponent"
  item of `gauss.md` §7 is out of reach of toric Liouville with polynomials of degree `≤ 7`). At `k = 4`: `3m` for
  `m ≤ 9`, against the balanced `7/2` — the all-profile and balanced constants genuinely differ.
* **(iii) Form of the needed inequality (§5).** `A_5(δ<1/30)` ⟺ `α^ess_e(-K) > 15/4` off an exceptional set (Lemma 1.3);
  `(L_5^∘)` is the same over balanced sequences. This is **strictly beyond Liouville** at a **rational** point (beyond
  `τ_e`, beyond the balanced Liouville constant), far above anything Seshadri-proportional (`≤ 2`), so no
  McKinnon–Roth-type theorem reaches it. It is **not McKinnon's conjecture** (which concerns `α_e = 2`, attained in
  `Z_diag`); it is an essential (off-`Z`) statement of Vojta type: Vojta for `(Bl_e X, E, {∞})` gives `α^ess ≥ 4 = dim X`
  (Prop 5.1), which is `gauss.md` Thm 4.2. lit18's "drifting phase" dismissal of Huang / McKinnon–Roth no longer applies
  in `T_5` (the phase is quotiented out, the target is fixed), but the tools fail for the sharper reason above.
* **(iv) Toric arguments (§6): no, for proved reasons.** Sections: capped at `δ = 1/30` (mean width; Thm 3.3).
  Archimedean toric data: the points lie on the compact real torus where `Log_∞` and every moment map are constant, so
  proximity to `e` is invisible to every toric invariant (Prop 6.1). Cocharacter/tropical data: fake-blind (record's
  balanced 5-row fakes at `R''^{7/15}`). Curves: covering rational families stop at `α = 4`; the minimal-degree
  candidates for extremal balanced curves form a 0-dimensional (expected) PTE system, and a numerical search with
  calibration found **no** nondegenerate real solution (§7; evidence, inconclusive). So toric means neither prove nor
  refute `(L_5^∘)`.

Constraints of the owner: respected (§8).

## 1. Setup: the toric model, the height class, the twist

**Conventions.** *Theorem/Proposition/Lemma* = proved here with all quantifiers, or cited with the place of the proof.
*Exact* = finite exact computation (integers, `Fraction`, sympy, or ranks over a prime field used only in the direction
"full rank mod p ⇒ full rank over Q"), scripts in `toric_checks/`. *Numerical* = floating point, used only to locate
objects that are then certified, or labelled as numerical. *Heuristic* = no proof. *Conditional* = proved from a named
hypothesis. Nothing here is refereed.

### 1.1 The variety

* `N = Z^5/Z·1` (cocharacters of `T_5` over `Q(i)`), `M = {m ∈ Z^5 : Σ m_j = 0}` (characters), `x^m = Π u_j^{m_j}`.
  Over `Q(i)`, `T_5 ≅ G_m^5/G_m` via `u_j = x_j + i y_j`; complex conjugation `σ` sends `u_j` to `v_j = u_j^{-1}` on the
  norm-one torus, so `σ` acts on `N` and `M` by `-1`.
* `Σ` = the braid (permutohedral) fan: rays `ρ_A = 1_A` for `∅ ≠ A ⊊ [5]` (30 rays), cones = chains
  `A_1 ⊊ … ⊊ A_r`; 120 maximal cones (permutations). `σ(1_A) = -1_A = 1_{A^c}`, so `σ(Σ) = Σ`.
  `X := L̄_5^tw` is the Galois form over `Q` of the smooth projective toric 4-fold `X_Σ` (the Losev–Manin space
  `LM_5`), descended along this action (the fan is `σ`-stable and carries the `σ`-stable ample class `Σ_A a E_A`
  with any strictly convex symmetric support function, so descent is effective). `e ∈ T_5(Q)` is the identity.
  `E_A` = the toric divisor of `1_A` (defined over `Q(i)`); `E_A + E_{A^c}` is defined over `Q`.
* *Height class.* `B := -K_X = Σ_A E_A`. By `gauss.md` Lemma 1.1 (re-derived by the killer), for the class
  `y = [ζ]` of a primitive tuple on `|ζ| = R''`: `h_{E_A + E_{A^c}}(y) = log n_A` exactly, and `h_B(y) = W' = log R''^2`.
  By `σ`-symmetry `h_{E_A}(y) = h_{E_{A^c}}(y) = ½ log n_A` (heights over `Q(i)`, absolute normalisation).
* *Distance.* `dist(y, e) := max_j |ζ_j/ζ_1 - 1| = max_j|ζ_j - ζ_1|/R''` (archimedean; comparable near `e` with the
  distance of any smooth metric). Since `max_j|ζ_j - ζ_1| ≤ diam ≤ 2 max_j |ζ_j-ζ_1|`,
  `diam/(2R'') ≤ dist(y,e) ≤ diam/R''`.
* *Diagonals.* `Δ_{ab} = {u_a = u_b}` (codimension-one subtori through `e`); `Z_diag := ∪_{a<b} Δ_{ab}`. The tuple is
  distinct iff `y ∉ Z_diag`.

**Proposition 1.1 (the polytope of `B`; proved, exact check `basics.py`).** For `k ≥ 3` points (so `M ⊂ Z^k`):
(a) `P_B = R_k := {x ∈ M_R : ⟨x, 1_A⟩ ≤ 1 for all A} = {x : Σx_j = 0, Σ_j |x_j| ≤ 2} = conv{e_a - e_b : a ≠ b}`,
the root polytope of type `A_{k-1}`;
(b) `R_k = Δ + (-Δ)` with `Δ = conv(e_1,…,e_k)`, i.e. `B = H + H'`, where `H = π_1^*O_{P^{k-1}}(1)`,
`H' = π_2^*O_{P^{k-1}}(1)` for the two toric blow-downs `π_1, π_2 : X_Σ → P^{k-1}` (fans with rays `{1_{\{j\}}}`,
resp. `{1_{[k]∖j}}`);
(c) `B` is nef and big; for `k ≥ 4` it is not ample; `B^{k-1} = binom(2k-2, k-1)` (`= 6, 20, 70` for `k = 3,4,5`);
(d) for every `v ∈ N`, `width_v(R_k) := max_{R_k}⟨·,v⟩ - min_{R_k}⟨·,v⟩ = 2(max_j v_j - min_j v_j)`.

*Proof.* (a) The vertex of `P_B` attached to the maximal cone of the permutation `s` solves `⟨m, 1_{\{s_1..s_t\}}⟩ = -1`
(`t < k`), `Σ m = 0`, so `m = e_{s_k} - e_{s_1}`. Hence `P_B = conv{e_a - e_b}`; and `⟨x, 1_A⟩ ≤ 1 ∀A` is
`max_A Σ_{j∈A}x_j = Σ_j (x_j)_+ ≤ 1`, i.e. `‖x‖_1 ≤ 2` on `Σx = 0`, whose vertices are the `e_a - e_b`. (b)
`conv(e_a) + conv(-e_b) = conv(e_a - e_b)` (the points `e_a - e_a = 0` are interior). Each braid cone of the chain
`A_1 ⊊ … ⊊ A_{k-1}` lies in the cone spanned by `{1_{\{j\}} : j ∈ A_{k-1}}`, resp. by `{1_{[k]∖j} : j ∉ A_1}`, so both
`π_i` are toric morphisms, birational, isomorphisms over the torus; `π_1^*O(1)` has polytope `Δ` (translated),
`π_2^*O(1)` has `-Δ`. (c) Nef and big: `P_B` is full-dimensional and is the polytope of the support function. For
`k ≥ 4` the vertex `e_{s_k} - e_{s_1}` is shared by `(k-2)!` maximal cones, so the support function is not strictly
convex and `B` is not ample. `B^{k-1} = (k-1)!·vol(R_k)`; the Ehrhart polynomial was computed exactly
(`35m^4/12 + 35m^3/6 + 85m^2/12 + 25m/6 + 1` at `k = 5`), leading coefficient `35/12`, so `B^4 = 70`. (d) The maximum
of `⟨e_a - e_b, v⟩` is `max v - min v`, the minimum its negative. `□`

**Lemma 1.2 (degrees of one-parameter orbits over `Q`; proved).** For primitive `v ∈ N` let `C_v` be the closure in
`X` of the image of the `Q`-homomorphism `φ_v : R^1_{Q(i)/Q}G_m → T_5`, `t ↦ (t^{v_j})_j`. Then `C_v` is a smooth
`Q`-conic through `e` (over `Q(i)` it is the closure of the 1-PS `λ_v` and of `λ_{-v}`), and `B·C_v = 2(max v - min v)`.
Its rational points near `e` are `y = φ_v(P/P̄)`, `P = a + bi` primitive, with `h_B(y) = (max v - min v)·log N(P)` and
`dist(y,e) ≍_v |b|/|P|`.
*Proof.* `φ_v` is defined over `Q` because each power map of `R^1G_m` is. Over `Q(i)` the image is `λ_v(G_m)`, and
`σ∘λ_v = λ_{-v}∘σ` has the same image, so the closure is `σ`-stable. The degree of `L` on the closure of a 1-PS orbit
through a torus point is `width_v(P_L)` (pull back the toric sections `x^m`, `m ∈ P_L`, to `t ↦ t^{⟨m,v⟩}`); with
Prop 1.1(d) this is `2(max v - min v)`. The orbit map is an immersion, so `C_v` is smooth at `e`. For the points: write
`v` with `min v = 0`, `max v = m'`; then `ζ_j = P^{v_j}P̄^{m'-v_j}` is a primitive tuple (gcd `= P^0 P̄^0`) on
`|ζ| = |P|^{m'}`, so `h_B = log|P|^{2m'}`, and `ζ_j/ζ_1 - 1 = (P/P̄)^{v_j - v_1} - 1`, `|P/P̄ - 1| = 2|b|/|P|`. `□`

### 1.2 The approximation constants

Following McKinnon–Roth (Invent. Math. 200 (2015) 513–583; definitions checked against the abstract and the `P^1`
example only, arXiv being blocked; nothing below depends on their theorems), for a sequence `y_i ∈ X(Q)` with
`dist(y_i, e) → 0` put `α({y_i}) := inf{γ ≥ 0 : dist(y_i,e)^γ H_B(y_i) is bounded}`, and for a set `Ω ⊆ X(Q)`:
`α_e(B; Ω) := inf` over such sequences in `Ω`; the **essential** constant is
`α^ess_e(B; Ω) := sup_{Z ⊊ X closed} α_e(B; Ω ∖ Z)`. Small `α` = good approximation. On `P^1` at a rational point with
`O(1)`, `α = 1 = ε` (Liouville is sharp).

**Lemma 1.3 (dictionary; proved).** Let `Ω_dist = T_5(Q) ∖ Z_diag`.
(a) If `A_5(δ)` holds with exceptional set `Z`, then every sequence in `Ω_dist ∖ Z` tending to `e` has
`α ≥ 4/(1+2δ)`; so `α^ess_e(B; Ω_dist) ≥ 4/(1+2δ)`.
(b) Conversely, if `α_e(B; Ω_dist ∖ Z) > 4/(1+2δ)` for some proper closed `Z`, then `A_5(δ')` holds for every `δ' > δ`,
with exceptional set `Z ∪ Z_diag`.
(c) (Balanced form.) For a point in one unit class, `η`-balanced at scale `w`, off `Z_diag`, put
`μ(y) := λ(y)/w` with `λ(y) = -log dist(y,e)`. Then
`(4 - μ)w - O(ηw) - O(1) ≤ H_Σ(y) ≤ 10(4 - μ)w + O(ηw) + O(1)`, and every `log|t_ab| ≥ 0`. Hence
(Liouville) `μ ≤ 4 + O(η)`, i.e. `dist ≥ H_B^{-4/15 - O(η)}`, and `(L_5^∘)` is equivalent to: *there are `c', η_0` and `Z`
such that balanced points off `Z` with `w ≥ w_0` have `μ ≤ 4 - c'`*, i.e. the essential approximation constant over
`η_0`-balanced sequences is `≥ 15/(4-c') > 15/4`.
*Proof.* (a), (b): `dist ≤ diam/R''` and `H_B = R''^2` turn `diam ≤ R''^{1/2-δ}` into
`dist ≤ H_B^{-(1+2δ)/4}`, and `dist ≥ diam/(2R'')` gives the converse with the constant `2` absorbed by `δ' > δ`;
distinctness is `y ∉ Z_diag`. (c) `|ζ_a - ζ_b| = 2|t_ab| R'' e^{-d_ab/2}` and each pair is separated by 8 of the 15
cuts, so under `η`-balance `d_ab/2 = 4w + O(ηw)`. Now `log|t_ab| = log|ζ_a - ζ_b| - log R'' + d_ab/2 - log 2`, and
`log|ζ_a - ζ_b| - log R'' ≤ log(2 dist) = -λ + log 2` for every pair, with equality up to `O(1)` for a pair realising
`max_j |ζ_j - ζ_1|`. So every `log|t_ab| ≤ (4-μ)w + O(ηw) + O(1)` and one of them is `≥ (4-μ)w - O(ηw) - O(1)`;
within one unit class every `log|t_ab| ≥ 0`. Summing gives the two bounds, and `0 ≤ (4-μ)w + O(ηw) + O(1)` is the
Liouville bound. If `H_Σ ≥ cw` then `μ ≤ 4 - c/10 + O(η)`; if `μ ≤ 4 - c'` then `H_Σ ≥ (c' - O(η))w - O(1)`. This is
the equivalence (with `η_0` small against `c`, `c'`). `□`

So, in this language, **`A_5(δ)` for some `δ < 1/30` is the statement `α^ess_e(B; Ω_dist) > 15/4`, and `(L_5^∘)` is the
same statement restricted to balanced sequences.**

## 2. (i) Seshadri constant at e

`ε_e(L) := inf_{C ∋ e} (L·C)/mult_e C` over irreducible curves of `X_{\bar Q}`; equivalently
`sup{t : π^*L - tE nef}` on `Y = Bl_e X`.

**Theorem 2.1 (proved).** `ε_e(B) = 2`, for `B = -K_X` on `X = L̄_5^tw` (and on the split form `LM_5`). More generally
`ε_e(-K_{LM_k}) = 2` for every `k ≥ 2`. The infimum is attained by the 30 (over `Q`: 15) conics `C_{1_A}` of Lemma 1.2,
all of which are defined over `Q` and all of which lie in `Z_diag` (for `k = 5`).
*Proof.* Upper bound: `B·C_{1_A} = 2(1 - 0) = 2` and `C_{1_A}` is smooth at `e` (Lemma 1.2). Lower bound: by
Prop 1.1(b), `B = H + H'`. Let `C ∋ e` be an irreducible curve over `\bar Q`. The exceptional locus of `π_1` lies in
the toric boundary, which misses `e`, so `C` is not contracted, and `π_1` is an isomorphism near `e`; hence
`H·C = deg(C → π_1 C)·deg π_1(C) ≥ mult_{π_1 e}(π_1 C) = mult_e C` (a hyperplane through `π_1(e)` not containing
`π_1(C)` meets it with multiplicity `≥ mult`). The same for `H'`. So `B·C ≥ 2 mult_e C`. Distinct tuples: `C_{1_A}` is
the set of classes of `(ζ_j) = (P` for `j ∈ A`, `P̄` for `j ∉ A)`, which has a repeated point as soon as `|A| ≥ 2` or
`|A^c| ≥ 2`, always the case at `k = 5`. `□`

*Galois descent.* Nefness of `π^*B - tE` and intersection numbers are invariant under field extension, so `ε_e` is the
same over `Q`, `Q(i)` and `\bar Q`; the twist enters only through which curves are defined over `Q`, and the extremal
`C_{1_A}` are (images of `Q`-homomorphisms `R^1G_m → T_5`; `C_{1_A} = C_{1_{A^c}}`, so 15 of them). In general a
`σ`-anti-invariant cocharacter `v` always gives a `Q`-curve here because `σ` acts by `-1` on all of `N`: **every**
1-PS orbit through `e` descends to a `Q`-conic, although `T_5` is anisotropic (`N^σ = 0`).

**Corollary 2.2 (what Seshadri-type theorems can give; proved).** Any inequality of the form
`α_e(B) ≥ c·ε_e(B)` with `c ≤ 1` (the McKinnon–Roth Roth-type theorem has `c = 1/2`; any theorem valid on `P^1` has
`c ≤ 1`, because Liouville is sharp at rational points of `P^1`) gives at most `α ≥ 2`, i.e. `dist ≥ H_B^{-1/2}`,
i.e. `diam ≥ c R''^0`: nothing. And `α_e(B; T_5(Q)) ≤ 2` is attained on `C_{1_A}` (Lemma 1.2 with `b = 1`), so
`α_e(B) = ε_e(B) = 2` *on the whole variety*, computed in `Z_diag`, in accordance with McKinnon's principle that best
approximations lie on rational curves (here: on minimal-degree 1-PS conics, the analogue for this twisted form of
Huang's "minimal centred primitive collections"; Huang's theorem is for split toric varieties and is not used).

So the Seshadri constant of the height class at `e` is **irrelevant to `A_5`**: it is computed inside the exceptional
set `Z_diag`, and the needed value `15/4 < α` is far above it. What matters is the **essential** constant, i.e. the
Seshadri / Liouville constants of the variety in which the diagonals are made exceptional, and the restriction to
balanced points. Two "restricted Seshadri constants" are relevant:

**Proposition 2.3 (Seshadri constant off the diagonals; proved bounds).** Let
`ε^diag_e(B) := inf{(B·C)/mult_e C : C ∋ e irreducible, C ⊄ Z_diag}`. Then `10/3 ≤ ε^diag_e(B) ≤ 8`, and for every
irreducible curve `C ∋ e`, `C ⊄ Z_diag`, whose toric intersection profile is balanced
(`E_A·C` independent of `A`), `(B·C)/mult_e C ≥ 15/4`.
*Proof.* Upper: a 1-PS `C_v` with `v = (0,1,2,3,4)` is not in `Z_diag` and has `B·C_v = 8`, `mult = 1`. Lower: the
Vandermonde `V = Π_{a<b}(x_a - x_b)·(x_1⋯x_5)^{-2}` is a section of `3B` (its Newton polytope is
`Perm(-2,…,2) ⊂ 3R_5`, `basics.py`) with `ord_e V = 10` and zero set exactly `Z_diag` inside `T_5`. For `C ⊄ Z_diag`
through `e`, `V|_C ≠ 0` and `3B·C ≥ ord_e(V|_C) ≥ 10·mult_e C`. Balanced: write `D_V = Σ_A a_A E_A` (the divisor of `V`
as a section of the toric bundle with polytope `Q_V`), `a_A + a_{A^c} = width_{Q_V}(1_A) = |A||A^c|`. If
`E_A·C = κ` for all `A`, then (V is a section of the nef toric bundle `O(D_V)`, nonzero on `C`)
`D_V·C = κ Σ_A a_A = κ Σ_{cuts}|A||A^c| = κ(5·4 + 10·6) = 80κ ≥ ord_e(V|_C) ≥ 10 mult` and `B·C = 30κ`, so
`mult ≤ 8κ`, i.e. `(B·C)/mult ≥ 30κ/(8κ) = 15/4`. `□`

The value `15/4` reappears in §§3–5 as the balanced Liouville constant; Proposition 2.3 is its curve-wise form.

## 3. The Liouville (pseudo-effective) threshold tau_e and toric Liouville

`τ_e(L) := sup{t : π^*L - tE is big on Y = Bl_e X} = lim_m t_max(m)/m`, `t_max(m) := max{ord_e f : 0 ≠ f ∈ H^0(X, mL)}`.
Always `ε_e(L) ≤ (L^4)^{1/4} ≤ τ_e(L)`. For `B`: `H^0(mB)` = Laurent polynomials with support in `S_m := mR_5 ∩ M`
(over `Q(i)`; a `σ`-stable condition, so by Galois descent for vector spaces the extremal `f` can be taken over `Q`).

**Lemma 3.1 (vanishing order = interpolation regularity; proved).** For finite `S ⊂ M` let `reg(S)` be the least `d`
with `Poly_{≤d}|_S = Q^S` (polynomial functions on `M ⊗ R` restricted to `S`). Then
`max{ord_e f : 0 ≠ f, supp f ⊆ S} = reg(S)`. In particular `t_max(m) = reg(S_m)`.
*Proof.* Put `x = exp(y)`, `y ∈ M^∨⊗C` near 0. `f = Σ_{s∈S} c_s x^s` has `ord_e f ≥ t` iff `Σ_s c_s ⟨s,y⟩^j = 0` for
`j < t`, i.e. iff `c ⊥ Poly_{<t}|_S`. A nonzero such `c` exists iff `Poly_{≤ t-1}|_S ≠ Q^S`, i.e. iff `t ≤ reg(S)`. `□`
(This is the regularity `reg` of `threshold/construct.md` §2.4, here applied to dilates of the root polytope.)

**Theorem 3.2 (toric Liouville; proved).** Let `0 ≠ f ∈ H^0(mB)` have coefficients in `Z[i]`. Then every
`y = [ζ] ∈ T_5(Q)` with `f(y) ≠ 0` satisfies `dist(y,e)^{ord_e f} ≥ c_f H_B(y)^{-m}`. More precisely (forced
divisibility), with `Q_f` = Newton polytope of `f` and `D_{Q_f} = Σ_A a_A(Q_f)E_A` its toric divisor,
`ord_e f · λ(y) ≤ h_{D_{Q_f}}(y) + O_f(1) = ½ Σ_{cuts {A,A^c}} width_{Q_f}(1_A) · log n_A + O_f(1)`.
*Proof.* First form: with `x_j = ζ_j/ζ_1` and `s ∈ S_m`, `Σ_j s_j = 0` and `Σ_j (s_j)_- ≤ m`, so
`x^s = Π_j ζ_j^{(s_j)_+} ζ̄_j^{(s_j)_-} / N'^{Σ(s_j)_-}` and `f(y) = G(ζ)/N'^m` with
`G = Σ_s c_s Π ζ_j^{(s_j)_+}ζ̄_j^{(s_j)_-} N'^{m - Σ(s_j)_-} ∈ Z[i]`. `G ≠ 0` gives `|f(y)| ≥ N'^{-m} = H_B(y)^{-m}`,
while `|x_j| = 1` and Taylor's formula give `|f(y)| ≤ C_f dist^{ord_e f}`. Second form: the product formula over `Q(i)`
for the section `f` of the toric bundle `O(D_{Q_f})` with its toric (Batyrev–Tschinkel) metrics, whose finite local
heights are `log max_{s ∈ Q_f} |x^s|_v ≥ log|f(y)|_v`; the archimedean toric local height vanishes (`|x_j|_∞ = 1`);
and `h_{E_A}(y) = ½ log n_A` (§1.1), `a_A + a_{A^c} = width_{Q_f}(1_A)`. `□`

**Theorem 3.3 (all toric Liouville arguments are governed by `τ_e(B)`; proved).** Let `F` be any finite set of
nonzero Laurent polynomials and `Z_F = ∪_{f∈F} Z(f) ∩ T_5`. The best exponent that Theorem 3.2 (second form, any
`f ∈ F` chosen per point) gives for points off `Z_F` of arbitrary profile is: such points satisfy
`diam ≥ R''^{1/2 - δ_F - o(1)}` with `δ_F := inf over the multiplicative semigroup generated by F of
max_A width_f(1_A)/ord_e f - 1/2` (a minimum over the convex hull of the points `width_f/ord_e f`, `f ∈ F`), and there
is a profile at which no `f` in that semigroup does better. Moreover
`inf_F δ_F = 2/τ_e(B) - 1/2`.
*Proof.* For a point with cut weights `ℓ_A = log n_A / W'` (a probability vector on the 15 cuts), Theorem 3.2 gives
`λ ≤ (W'/2) Σ_A width_f(1_A) ℓ_A / ord_e f`; `A_5(δ)`-points have `λ ≥ (¼ + δ/2)W' - O(1)`. The set
`K = {width_f/ord_e f}` over the semigroup is convex (widths and orders add under products: Newton polytopes add) and
`ℓ` ranges over a compact simplex, so by the minimax theorem
`max_ℓ min_{w∈K} ⟨w,ℓ⟩ = min_{w∈K} max_ℓ⟨w,ℓ⟩ = min_{w∈ K} max_A w_A`, which is the claimed `δ_F` after
`¼ + δ/2 ≤ ½⟨w,ℓ⟩`. For the last claim replace `f` by `f·f^*`, `f^*(x) := \bar f(x^{-1})`: on norm-one points
`|f^*(y)| = |f(y)|`, the order doubles and `Q_{ff^*} = Q_f + (-Q_f)` is centrally symmetric, so it lies in `m'R_5`
with `m' = max_A width_{Q_f}(1_A)` (`R_5` is cut out by `⟨·,1_A⟩ ≤ 1`). Hence
`sup_f 2 ord_e f / max_A width_f(1_A) = sup_m t_max(m)/m = τ_e(B)`. `□`

So **the strongest conclusion any toric Liouville argument (any auxiliary Laurent polynomial, forced divisibility
included, finitely many of them, off their zero sets) can give for arbitrary profiles is `A_5(δ)` for
`δ > 2/τ_e(B) - 1/2`**, and for balanced profiles `δ > 1/30` (§4). Bounds on `τ_e(B)`:

**Theorem 3.4 (`τ_e(B)`; proved bounds + exact finite computation).**
(a) `10/3 ≤ τ_e(B) ≤ 15/4`.
(b) (exact) `t_max(m) = ⌊10m/3⌋` for `m = 1, …, 7`; i.e. no section of `mB`, `m ≤ 7`, vanishes at
`e` to order `> 10m/3`.
(c) (exact, comparison) For `k = 4` (the 3-fold `LM_4`), `t_max(m) = 3m` for `m = 1, …, 9`, although the balanced
Liouville constant there is `7/2`; for `k = 3` (`dP_6`), `τ_e = 3` exactly (attained by the Vandermonde and equal to the
mean-width bound).
*Proof.* (a) Lower: `V^a·Π(3-point Vandermondes)` give `ord ⌊10m/3⌋` in `H^0(mB)` (`V ∈ H^0(3B)` with order 10; a
3-point Vandermonde `V_{abc}` has polytope `Perm(-1,0,1)` on coordinates `a,b,c`, inside `R_5`, order 3). Upper:
`threshold/construct.md` Thm B (homogeneous form B1, refereed) gives `ord_e f ≤ 2^{-4}Σ_{P⊆[5]} width_{Q_f}(1_P)
≤ 2^{-4}·30·2m = 15m/4` for `Q_f ⊆ mR_5`. (b) Lower bound as in (a); upper bound: by Lemma 3.1 it suffices that
`Poly_{≤⌊10m/3⌋}|_{S_m}` has full rank `|S_m|` over `Q`, and full rank modulo the prime `2^{61}-1` implies full rank
over `Q`. Computed (`reg.py`, `reg5.out`, m ≤ 5, by increasing `d` until full rank, which also shows rank deficiency
at `d = ⌊10m/3⌋ - 1`, consistent with the explicit lower bound; `reg_cert.py` with all monomial columns: m = 6, rank 5321 = |S_6|; m = 7, rank 9381 = |S_7|). (c) Same
programme for `k = 4` (`reg.py 4 9`); `k = 3`: `V_{123} ∈ H^0(B)` has order 3 and the mean-width bound is
`2^{-2}·6·2 = 3`. `□`

**Remarks.** (1) The data say `τ_e(B) = 10/3` is very likely (**evidence**, not proof): no section of `m(-K)`,
`m ≤ 7`, beats the Vandermonde products, so on the data **toric Liouville does not even improve Cilleruelo–Córdoba (`δ > 1/10`) off a
closed set**, let alone reach `δ < 1/30`. A proof of `τ_e(B) ≤ 10/3` would follow from the weighted mean-width
inequality `ord_e f ≤ (1/6) Σ_{|A| = 2} width_{Q_f}(1_A)` for all `f` (the 2|3 cuts; equality for every root binomial
and for `V`); I did not prove it. It cannot come from rational curves: a covering family of rational curves through
`e` has `(B·C)/mult ≥ 4` (Prop 4.3 below). (2) At `k = 4` the gap between `τ_e(B) = 3` (data) and the balanced
constant `7/2` shows that the two constants really differ: the all-profile Liouville threshold is set by unbalanced
(2|2-heavy) profiles, the balanced one by the mean width.

## 4. (ii) Approximation constants: 1-PS curves, rational curves, balanced locus, the k = 4 curves

**Proposition 4.1 (approximation along a rational curve; proved).** Let `ν : P^1_Q → X` be non-constant, `C = ν(P^1)`,
`ν(O) = e` for some `O ∈ P^1(Q)`, and `μ_0 := mult_O ν^*E` (contact order of the branch at `O`). Then the rational points
`ν(t)`, `t → O`, have `α = (B·C)/μ_0` (with `B·C := deg ν^*B`), and no sequence of points of `C` approaching `e`
through this branch does better. If `E_A·C = κ` for all `A` (a **balanced curve**), these points are `η`-balanced
with `η = O(1/w)`, and `(4-μ)w - O(1) ≤ H_Σ ≤ 10(4 - μ)w + O(1)` with `μ = μ_0/(2κ) = 15μ_0/(B·C)` (for points in one
unit class; across unit classes read `t_ab := |c_ab|/2`, which changes `H_Σ` by `O(1)`, `kill_gauss.md` §1.2 item 2).
*Proof.* `B` is base-point free (nef toric), so `h_B(ν(t)) = (B·C)h(t) + O(1)`; likewise
`h_{E_A}(ν(t)) = h_{ν^*E_A}(t) + O(1) = (E_A·C)h(t) + O(1)` for `C ⊄ E_A` (Weil height machine on `P^1`). Near `O`,
`dist(ν(t),e) ≍ |t - O|^{μ_0}`; for `t ∈ Q`, `|t - O| ≥ c/H(t)` with equality up to constants on infinitely many `t`
(Liouville and its converse on `P^1` at a rational point). So `dist ≍ H(t)^{-μ_0} = H_B^{-μ_0/(B·C)}`. In the balanced
case all `log n_A = 2κ h(t) + O(1)`, so `η = O(1/h)`; then use Lemma 1.3(c) with `w = 2κ h(t) + O(1)`,
`λ = μ_0 h(t) + O(1)`, and `B·C = 30κ`. `□`

**Corollary 4.2 (the values; proved).**
(a) *One-parameter subgroups.* Along `C_v`: `α = B·C_v = 2(max v - min v)` (Lemma 1.2: `μ_0 = 1`). Minimum `2` at
`v = 1_A` (inside `Z_diag`); off `Z_diag` (all `v_j` distinct) `α ≥ 8`, i.e. diameter `≍ R''^{3/4}`. Their profile is
maximally unbalanced: at most 4 of the 15 cut norms are nontrivial (one chain).
(b) *Balanced curves.* For a balanced curve `C ⊄ Z_diag` through `e`, `α_C = 30κ/μ_0 ≥ 15/4`, i.e. `μ_0 ≤ 8κ`
(Prop 2.3), with equality iff `C̃` is `V`-null: `(π^*D_V - 10E)·C̃ = 0`. Equality = a rational curve of balanced points
with bounded residues at the exponent `7/15` (Prop 4.1 with `μ = 4`).
(c) *Restricted to the balanced locus.* For every sequence of balanced points (`η_i → 0`, one unit class, off `Z_diag`):
`α ≥ 15/4` (Lemma 1.3(c), Liouville). On curves this is Prop 2.3; on arbitrary sequences it is the pair bound. So
`α_e(B; balanced, off Z_diag) ≥ 15/4`, and it equals `15/4` iff balanced bounded-residue (`H_Σ = o(w)`) sequences exist,
on curves iff a `V`-null balanced rational curve exists (§7).

**Proposition 4.3 (covering rational families cannot reach the Liouville constant at `k = 5`; proved).** Let
`{ν_s}` be a family of rational curves through `e` (contact `≥ μ_0` at a fixed point `O` of `P^1`, `O` mapping to
`e`), whose members cover `X`. Then `(B·C)/μ_0 ≥ k - 1` (= 4 at `k = 5`, 3 at `k = 4`).
*Proof.* This is the inequality `#fixed·μ_0 ≤ -K_Y·C̃ = B·C - (k-2)E·C̃` of `gauss.md` Thm 4.3(ii) (free-curve
nefness of the subsheaf `F_ν ⊂ ν^*T_Y` of fields with the prescribed contact; re-derived in `kill_gauss.md` §4) with
one fixed point and `E·C̃ ≥ μ_0`: `μ_0 + (k-2)μ_0 ≤ B·C`. `□`

So at `k = 5` a covering family of rational curves through `e` gives at best `α = 4` (the counting / Vojta value,
`δ = 0`), strictly above the balanced Liouville value `15/4`. **The extremal objects for `(L_5^∘)` can therefore not
move in a covering family of rational curves**; whatever attains `15/4` densely must be rigid curves of unbounded
degree, non-rational curves, or sporadic points (cf. gauss Thm 4.3).

**Theorem 4.4 (the `k = 4` family is the `k = 4` instance of these extremal curves; proved, exact).** Each member
(`(X,m)` fixed) of `gauss.md` Theorem 4.1's family is a rational curve `C ⊂ L̄_4^tw` through `e` with:
`E_A·C = 1` for all 14 rays (balanced, class `2γ_4`), `B·C = 14`, a single branch at `e` (`s = ∞`) of multiplicity 4,
so `α_C = 7/2` = the balanced Liouville constant at `k = 4` (exponent `3/7`), and `C̃` is `V_4`-null:
`D_{V_4}·C = Σ_cuts |A||A^c| = 24 = ord_e V_4 · mult = 6·4`. These curves cover `T_4` (gauss, exact Jacobian rank 3),
which is allowed at `k = 4` since `7/2 ≥ k - 1 = 3` (Prop 4.3).
*Proof/check.* `u_x = P_x/P̄_x` with `P_x = Π_{T∋x}(s - r_T)`: at `s = r_T` the valuation vector is `1_T`, at
`s = r̄_T` it is `-1_T = 1_{T^c}`; the 14 roots are distinct, so every ray is met once. At `s = ∞`,
`u_x - 1 = 2i c_x/P̄_x(s)` with `c_x = Im P_x ≠ 0` constant and `deg P_x = 4`, so the branch has multiplicity 4.
`k4class.py` verifies this exactly for the instance `(X,m) = (-3,3)` and four other rational parameters (constant
`Im P_x`, 14 rays each hit once, `B·C = 14`, `mult = 4`, ratio `7/2`). `□`

The comparison with `k = 5` is the whole story:

| | `k = 4` | `k = 5` |
|---|---|---|
| Seshadri `ε_e(B)` (computed in `Z_diag`) | 2 | 2 |
| `τ_e(B)` (all-profile Liouville) | `3` (data `m ≤ 9`) | `10/3` (data `m ≤ 7`); proved `∈ [10/3, 15/4]` |
| balanced Liouville constant (mean width) | `7/2` | `15/4` |
| covering rational families: `(B·C)/μ_0 ≥` | `3` | `4` |
| counting / Vojta value `dim X` | `3` | `4` |
| balanced extremal curves (expected real moduli dim, block degree `d`) | `1 + 2d` (`d = 1`: covering family, Thm 4.1) | `2 - 2d` (`d = 1`: finitely many expected; §7) |

At `k = 4` the balanced Liouville constant `7/2` lies **above** the covering threshold `3`, and is attained by a
covering family; at `k = 5` it lies **below** the covering threshold `4` and below the Vojta value `4`. The interval
`(15/4, 4)` of `α` (i.e. `δ ∈ (0, 1/30)`) is the whole content of `(L_5^∘)` / `A_5`.

## 5. (iii) Where A_5 sits: beyond Liouville, not McKinnon–Roth

All constants for `(X, B, e)`, in the normalisation `α = W'/λ` (`δ = 2/α - 1/2`):

| constant | value | status | meaning for arcs |
|---|---|---|---|
| `α_e(B)` over all of `X(Q)` | `2` = `ε_e(B)` | proved (Thm 2.1, Lemma 1.2) | computed by the conics `C_{1_A} ⊂ Z_diag`: repeated points |
| McKinnon–Roth Roth-type bound `ε/2` | `1` | their theorem (statement from abstract) | vacuous |
| any bound `α ≥ c·ε_e(B)`, `c ≤ 1` | `≤ 2` | proved (Cor 2.2) | vacuous |
| `α^ess_e(B; Ω_dist)` lower bound by toric Liouville | `≥ τ_e(B) ≥ 10/3` | proved (Thms 3.2–3.4) | Cilleruelo–Córdoba `δ > 1/10` |
| best possible toric Liouville, all profiles | `τ_e(B) ≤ 15/4`; `= 10/3` on data | proved / exact `m ≤ 7` | cannot pass `δ = 1/30`; likely not even `1/10` |
| balanced Liouville constant | `15/4` | proved (Lemma 1.3(c), Prop 2.3) | `δ > 1/30` excluded on balanced points |
| covering rational families through `e` | `≥ 4` | proved (Prop 4.3) | `δ ≤ 0` on such families |
| **needed**: `α^ess` on balanced sequences | **`> 15/4`** | open | `(L_5^∘)` |
| **needed**: `α^ess_e(B; Ω_dist)` | **`> 15/4`** | open | `A_5(δ)`, `δ < 1/30` |
| Vojta `V_E(5)` prediction | `≥ 4` | conditional (Prop 5.1) | `A_5(δ)` for all `δ > 0` |
| counting heuristic | `= 4 = dim X` | heuristic | `δ = 0` critical |

**Proposition 5.1 (Vojta at a rational point gives `α^ess ≥ dim`; conditional, proved implication).** Let `X` be
smooth projective of dimension `n` over `Q`, `x ∈ X(Q)`, `-K_X` big, `Y = Bl_x X`. If Vojta's Main Conjecture holds for
`(Y, D = E, S = {∞})`, then `α^ess_x(-K_X) ≥ n`. For `X = L̄_5^tw`, `x = e` this is `gauss.md` Thm 4.2 (`κ = 2δ`).
*Proof.* `K_Y = π^*K_X + (n-1)E`. Off `Z_ε`, `λ_{E,∞} + h_{K_Y} ≤ ε h_A + C_ε`; finite local heights of `E` are
bounded below, so `h_E ≥ λ_{E,∞} - O(1)` and `n λ_{E,∞} ≤ h_{-K_X} + ε h_A + O(1)`. With `-K_X` big, `h_A ≤ c h_{-K_X}`
off a proper closed set, so `λ_{E,∞} ≤ (1/n + εc/n) h_{-K_X} + O(1)`, i.e. `dist ≥ H_{-K_X}^{-1/n - ε'}`. `□`

**Answer to (iii).** The inequality `A_5` needs is
`α^ess_e(B; Ω_dist) > 15/4` (or the same over balanced sequences), where `15/4` is simultaneously the balanced Liouville
constant (exactly sharp for the pair product, Prop 2.3 / mean width) and an upper bound for every section-based
constant `τ_e(B)`. So it is **strictly beyond Liouville** in the precise sense: beyond the pseudo-effective threshold
of every auxiliary section, at a **rational** point. It is not of the form `α ≥ c·ε_e` (Seshadri-proportional bounds
are `≤ 2` here, computed in `Z_diag`), and not of the form "`α ≥ τ` off the stable base locus" (that gives at most
`τ_e(B) ≤ 15/4`; with the data, `10/3`). McKinnon–Roth's Roth-type theorem gains over Liouville only through the
`1/2`-loss for **irrational** points; at a rational point it is weaker than Liouville. So no McKinnon–Roth-type theorem
reaches `A_5`.

**Is it McKinnon's conjecture for this variety?** No. McKinnon's conjecture (as summarised in Huang's abstract: best
approximations are achieved on rational curves through the point, of minimal degree) is about the infimum
`α_e(B) = 2`, which is attained on the conics `C_{1_A}` and is compatible with it; that statement says nothing about the
essential constant. The needed statement is an **essential** approximation statement off an exceptional set, which is
of Vojta type (`D = E` on `Bl_e`, Prop 5.1). A hypothetical strengthening, *(McK^ess): Zariski-dense sets of
approximations of exponent `α < dim X` arise only from covering families of rational curves through `x` with
`(L·C)/mult ≤ α`*, would give `A_5(δ)` for all `δ > 0` by Prop 4.3; but (McK^ess) is not in the literature that I could
check (arXiv blocked; abstracts only) and is itself a Vojta-strength statement. **Re-examination of lit18's dismissal.**
In the `T_5` formulation the target is the fixed rational point `e`, the common phase is quotiented out exactly
(`[ζ]` is similarity-invariant, no phase is set to zero), so the "drifting phase" objection no longer applies. The
fixed-target tools still fail, now for a different, sharper reason: their constant is `α_e = ε_e = 2`, computed on
`Z_diag`, and the essential constant they would need lies beyond every Liouville/Seshadri-type bound.

## 6. (iv) Can a toric argument prove A_5(delta), delta < 1/30? Exactly why not

**Answer: no.** A toric argument can enter in three ways, and each is excluded at `δ < 1/30` for a proved reason.

**6.1 Sections (Newton polytopes, the fan, Liouville with forced divisibility).** Every argument that uses an auxiliary
Laurent polynomial (equivalently a section of a toric line bundle on `X`, or of `π^*L - tE` on `Bl_e X`, or of any
toric blow-up of `X` in the boundary) is Theorem 3.2, and Theorem 3.3 shows the best it can give for all profiles is
`δ > 2/τ_e(B) - 1/2 ≥ 2/(15/4) - 1/2 = 1/30`, and for balanced profiles exactly `δ > 1/30` (mean width,
`construct.md` Thm B; the pair product attains it). The fan, the polytope and the twist only change *which* sections
exist; Theorem B bounds all of them. *Proved.* Moreover, by Theorem 3.4(b) no section of `mB`, `m ≤ 7`, beats the
Vandermonde, so on the data even the intermediate range `δ ∈ (1/30, 1/10]` is out of reach of toric Liouville.

**6.2 The moment map / archimedean toric data.** **Proposition 6.1 (proved).** Every `y ∈ T_5(Q)` lies, at the
archimedean place, in the compact real torus `T_5(R) = {|u_j| = 1}/S^1 ≅ (S^1)^4`, on which the logarithm map
`Log_∞ = (log|x_j|)_j` and the moment map of any toric Kähler (or Batyrev–Tschinkel) metric are **constant**. All
archimedean toric local heights vanish on `T_5(Q)`, and the compact torus `T_5(R)` acts on `X(C)` by metric-preserving
automorphisms, transitively on `T_5(R)`. Hence no `T_5(R)`-invariant archimedean quantity distinguishes `y` from `e`:
the condition "`y` close to `e`" is invisible to every toric invariant at `∞`, and enters only through non-invariant
objects (characters `x^m - 1` and sections, i.e. 6.1).
*Proof.* `|ζ_j/ζ_1| = 1`. The toric metric of `O(D_P)` at `∞` is `‖x^m‖ = |x^m|/max_{m'∈P}|x^{m'}|`, which is `1` at
`|x_j| = 1`; the moment map is a function of `Log_∞`. `T_5(R)` acts by `u ↦ tu`, preserving `|x^m|`. `□`

So "convexity of the moment polytope" carries no information about the short arc: on the norm-one torus the moment
map is a single point (the barycentre). This is the toric form of the fact that the record's problem is an *angle*
problem.

**6.3 Cocharacter-lattice / tropical data (geometry of numbers on `N`).** The finite-place toric data of `y` are its
tropicalisations `trop_π(y) ∈ N` (one chain of cuts per split prime, `gauss.md` Lemma 1.1) and the norms `p`. An
argument that uses the primes only through `p` and `trop_π(y)` (any lattice-point counting, Minkowski or Siegel
argument on `N`, any convexity statement on the profile) cannot distinguish `y` from a **fake** point
`y_φ = Π_p e^{iφ_p·trop_π(y)}` with arbitrary real angles `φ_p` in place of `arg π_p`. For balanced 5-row profiles there
are fakes with diameter `≍ R''^{7/15}` satisfying every monomial gap (`gauss.md` §5 item 1, confirmed in
`kill_gauss.md` §5 item 1; the angle models of `overview/overview.md` §5.1 at `k = 5`). These violate `(L_5^∘)`. So no
argument of this type proves `(L_5^∘)` — this is the record's "local inputs" barrier restated torically, and it applies
verbatim. (Status: the fake statement is the record's, killer-checked, not refereed.)

**6.4 Curves (the refutation direction).** A toric *refutation* would be a Zariski-dense supply of balanced points at
`μ ≥ 4 - o(1)`; the toric curves through `e` are the 1-PS conics (`α ≥ 8` off `Z_diag`, maximally unbalanced, Cor 4.2)
and their deformations. Covering rational families stop at `α ≥ 4` (Prop 4.3); bounded-degree covering families with
infinitely many A-points are excluded (`gauss.md` Thm 4.3(ii)); and the natural minimal-degree candidates for
`V`-null balanced curves (Gaussian linear blocks) have expected real moduli dimension `2 - 2d`, i.e. at most finitely
many at `d = 1` and none expected for `d ≥ 2` (§7). The `k = 4` picture (a covering family *at* the Liouville
exponent, Thm 4.4) cannot repeat at `k = 5` because the covering threshold `k - 1 = 4` overtakes the balanced Liouville
constant exactly between `k = 4` and `k = 5`: the balanced Liouville constant at `k` points is
`α_bal(k) = 2/(1/2 + 1/(2^k-2)) = 4(1 - 2^{1-k})` (`7/2`, `15/4`, `31/8`, …), the covering threshold is `k - 1`, and
`k - 1 < α_bal(k)` iff `k ≤ 4` (`3 < 7/2`; `4 > 15/4`; for `k ≥ 5`, `k - 1 ≥ 4 > α_bal(k)`).

**Conclusion of (iv).** The needed statement is not provable by toric means (6.1–6.3, proved caps plus the record's
fakes) and not refutable by toric curve families (6.4). It is not McKinnon's conjecture (§5). It is a statement about
the essential approximation constant at a rational point **strictly between** the best Liouville constant (`15/4`,
balanced) and the Vojta/counting value (`4`), i.e. exactly the window `δ ∈ (0, 1/30)` that `V_E(5)` predicts and no
unconditional method in the record, toric or not, can enter.

## 7. The minimal-degree extremal curves at k = 5 (computation)

**Setting (proved reformulation).** A balanced rational curve through `e` of class `E_A·C = 1` for all 30 rays, with
its contact at a rational point (normalised to `s = ∞`), is given by Gaussian linear blocks `H_T = s - r_T`,
`T` nonempty `⊆ [4]` (cut `{T, [5]∖T}`), `z_5 = Π_T H̄_T`, `z_j = Π_{T∋j}H_T Π_{T∌j}H̄_T`, so
`u_j = z_j/z_5 = P_j/P̄_j`, `P_j = Π_{T∋j}(s - r_T)` (degree 8). Its contact at `s = ∞` is
`μ_0 = min_j ord_∞(u_j - 1) = min_j ord_∞(2i Im P_j/P̄_j)`, which equals `8 = 4·2κ` (the `V`-null value of Cor 4.2(b),
`κ = 1`) iff every `Im P_j` is a nonzero **constant**; otherwise `μ_0 ≤ 7` and `α_C ≥ 30/7 > 4`. By Newton's identities
this is: for `j = 1..4`, the 8-element multisets `R_j = {r_T : T ∋ j}` and `\bar R_j` have equal power sums
`p_1, …, p_7` — four coupled complex ideal Prouhet–Tarry–Escott problems of size 8 sharing 15 roots. Unknowns 30 real,
minus the real affine reparametrisations `s ↦ λs + ν` (2); equations 28: **expected dimension 0** (the `2 - 2d` of
`gauss.md` §4.1 at `d = 1`). Nondegeneracy: all `Im r_T ≠ 0` (else the block is real and the cut is trivial), the
`c_j = Im P_j` nonzero and distinct (else points coincide). Two easy necessary conditions (proved): (i) by the argument
principle, `Im P_j ≡ c_j ≠ 0` keeps `arg P_j(s)` in `(0,π)` or `(-π,0)` on `R`, so **each `R_j` has exactly 4 roots in
each half-plane**; (ii) any such curve has bounded residues and exponent exactly `7/15` (Prop 4.1), and the same curve
cannot occur in a positive-dimensional family (Prop 4.3 would be violated by a covering family; a non-covering
positive-dimensional family is not excluded by this argument).

**Numerical search (evidence, inconclusive).**

| method | script | starts | converged | nondegenerate |
|---|---|---|---|---|
| Levenberg–Marquardt on the 28 equations | `kpte.py 5`, `k5curves.py` | 300 + 200 | 297 + 200 | 0 |
| LM with Rabinowitsch variables `y_T Im r_T = 1`, `w_j c_j = 1` (47 × 47) | `k5pte.py`, `k5curves2.py` | 300 + ~100 | 0 | 0 |
| symmetric ansatz `r_{gT} = -\bar r_T`, `g = (12)(34)` (14 × 14) | `k5sym.py` | 400 | 372 | 0 |
| complexified system, Newton from random starts | `k5complex.py` | 400 | 0 | 0 |
| complexified system, Newton homotopy `F(z) = (1-t)F(z_0)` | `k5homotopy.py`, `k5hom_run.py` | 300 paths | all paths end at singular points | 0 |

Every converged real solution is degenerate, and in the 60 inspected (`k5inspect.py`) all four `c_j` vanish: each
`R_j` is closed under conjugation, i.e. the five points coincide. *Calibration:* the same LM code on the `k = 4`
system cut down to expected dimension 0 by three random affine constraints finds nondegenerate isolated real solutions
from 11–45% of random starts (`kpte_calib.py`, five instances); on the uncut 3-dimensional `k = 4` family 55/200
(`kpte.py 4`). So the absence at `k = 5` is mild evidence that **no balanced `V`-null rational curve of minimal degree
exists**, i.e. that the balanced Liouville constant `15/4` is not attained on curves of the lowest possible degree. It is
not a proof (the `k = 5` system has degree 7 and 28 unknowns; basins of isolated solutions can be tiny), and nothing
downstream depends on it: if such curves exist they are finitely many in each degree class and can be put into `Z`.

**What a positive answer would and would not mean.** A `Q`-rational solution would be an explicit rational curve of
balanced bounded-residue 5-tuples at exponent `7/15` (the `k = 5` analogue of one member of Theorem 4.1's family): it
would show that `Z ≠ ∅` is necessary in `(L_5^∘)` and that `α` restricted to balanced curves equals `15/4`, but it
would not refute `(L_5^∘)` (one curve, or finitely many, lie in a proper closed set).

## 8. Constraints, labels, computations, open points

**Owner's constraints.** (1) No uniform exceptional set from a fixed-`S` Subspace theorem is assumed: the only
exceptional sets are base loci of explicit sections (`Z(V) = Z_diag`, `Z_F`), the hypothetical `Z` of `A_5`/`(L_5^∘)`,
and Vojta's `Z_ε` for a fixed variety with `S = {∞}` (conditional, Prop 5.1). (2) No Gaussian gcd or denominator cost is
discarded: points are classes of **primitive** tuples, `h_B = log R''^2` exactly, and Theorem 3.2's proof keeps every
`N'`-power; forced divisibility enters through the exact toric heights. (3) No drifting phase is replaced by zero: the
common phase is removed by the similarity quotient (exact invariance of `[ζ]`), not set to zero; the fixed target `e`
is the genuine identity class. (4) No template-dependent threshold is treated as uniform: the computed `t_max(m)` are
statements for each fixed `m` and are used only as caps on a method; the curve computations are refutation attempts.

**Labels.** Proved: Prop 1.1, Lemmas 1.2, 1.3, Thm 2.1, Cor 2.2, Prop 2.3, Lemma 3.1, Thms 3.2, 3.3, 3.4(a), Prop 4.1,
Cor 4.2, Prop 4.3 (from `gauss.md` Thm 4.3(ii)'s inequality, killer-checked), Prop 6.1, the reformulation in §7.
Exact: Prop 1.1(c) numbers, Thm 3.4(b)(c), Thm 4.4 checks. Conditional: Prop 5.1. Evidence: `τ_e(-K) = 10/3`
(finite `m`), §7 search. Cited from the record (killer-checked, unrefereed): mean-width Thm B (refereed in
`threshold/`), balanced 5-row fakes (§6.3), `gauss.md` Thm 4.3(ii). Literature: McKinnon–Roth (Invent. Math. 200
(2015)), their Liouville companion (Eur. J. Math. 2016), Huang (Algebra & Number Theory 15 (2021) 461–512): definitions
and main statements known to me only from abstracts/search snippets (arXiv, ar5iv and the authors' pages are blocked
here); no argument above uses their theorems.

**Computations** (`toric_checks/`, run under `prlimit --as=3500000000`):

| script | content | result |
|---|---|---|
| `basics.py` | vertices of `P_{-K}` = roots; Ehrhart polynomials; `(-K)^{k-1}`; widths; mean width `Ω(R_k)`; Vandermonde ratio | `k = 3,4,5`: `(-K)^{k-1} = 6, 20, 70`; `Ω = 3, 7/2, 15/4`; Vandermonde `3, 3, 10/3` |
| `reg.py` → `reg5.out` | `reg(mR_k ∩ M)` by ranks mod `2^{61}-1` | `k = 5`: `3, 6, 10, 13, 16` (`m = 1..5`); `k = 4`: `3m` (`m = 1..9`) |
| `reg_cert.py` → `reg_cert_m6.out`, `reg_cert_m7.out` | full rank of `Poly_{≤⌊10m/3⌋}` on `S_m` | `m = 6`: rank 5321 = `|S_6|`; `m = 7`: rank 9381 = `|S_7|` (all 17550 monomials) |
| `k4class.py` | Theorem 4.1 curves: rays met, `B·C`, contact, ratio, `V`-null | instance + 4 parameters: 14 rays once each, `B·C = 14`, contact 4, ratio `7/2`, `D_V·C = 24 = 6·4` |
| `kpte.py`, `k5curves.py`, `k5curves2.py`, `k5pte.py`, `k5sym.py`, `k5complex.py`, `k5homotopy.py`, `k5hom_run.py`, `k5inspect.py`, `kpte_calib.py` | §7 search and calibration | no nondegenerate `k = 5` solution; calibration finds `k = 4` isolated solutions at 11–45% |

(`reg_cert_m67.out` records a failed variant: a random subset of `|S|+40` columns had rank 5035 at `m = 6`; this says
nothing, and the full-column run is the certificate.)

**Open points.**
1. Prove `τ_e(-K_{LM_5}) = 10/3` (equivalently the weighted mean-width inequality
   `ord_e f ≤ (1/6)Σ_{|A|=2} width_{Q_f}(1_A)`), or find a section beating the Vandermonde for some `m ≥ 8`. A section
   with `ord/m > 10/3` would give the first improvement of Cilleruelo–Córdoba at five points, off an explicit closed set
   (`δ > 2/τ - 1/2`); it would not touch the uniform theorem (needs `δ < 1/30`).
2. Decide whether the 0-dimensional `k = 5` PTE system of §7 has nondegenerate real (and then `Q(i)`-rational)
   solutions, e.g. by a certified homotopy or monodromy solver. Calibration value only.
3. The window `α^ess ∈ (15/4, 4)` at the rational point `e` — `(L_5^∘)` itself — is the expert question; the toric
   analysis shows it is not McKinnon's conjecture but an essential, Vojta-type statement strictly beyond every Liouville
   constant.
