# Round 27, lane stab: Question B at the top of the window as stability for balanced Liouville

Status: complete draft, unrefereed (round 27, lane stab). Scripts in `stab_checks/`, notes in `stab_notes.md`.

## 0. Summary and verdict

**Verdict: `obstruction` with one new theorem and one restricted theorem. No growth improvement, no uniform bound,
no proof of Question B.** `M ≤ (2/3+ε) log R/loglog R` stands.

1. **Elementary form (§1; proved, exact checks).** The coordinator's Gaussian Cox form is correct on the
   chain-condition locus (Prop 1.1), with two corrections: (i) off (Ch) the formula also produces primitive tuples
   whose block norms are *not* the `N(γ_A)`, so balance and level varieties must be taken on (Ch); (ii) the triangle
   relation is `Q_ab Q_bc = G^{(b)}_{abc} Q_ac` for the **pair cofactors** `Q_ab` (`Im Q_ab = t_ab`), not for the reduced
   chords `c_ab = 2i t_ab·unit`. New consequences: `x_b t_a - x_a t_b = G_ab t_ab`, a binary quadratic form of
   discriminant `-4t_ab^2` representing `G_ab t_ab^2`, and three equivalent elementary forms (Venn form, arccotangent
   cliques, singleton-decoupled residue form). `Stab(η,ε)` (near-tight balanced ⇒ `Z`) is **equivalent** to
   `(L_5^∘)` (Prop 1.3); `A_5(1/30-η)` implies it. The split twin is verified on the round-23 Pell curve.
2. **Main (§3.1).** Second-order integrality, made precise, is a **pencil principle** (Theorem 3.2, proved): two
   independent sections `f_1, f_2` with the same support and zero mean-width defect would map every bounded-residue
   balanced point to a point of bounded height on `P^1`, hence into finitely many fibres — exact algebraic relations,
   not a Liouville inequality. **Theorem 3.3 (new, proved for all `n`): measures attaining the mean-width bound are
   unique on their support set** (the uniqueness half of construct.md's Conjecture E). Corollary 3.4: every effective
   divisor on `Bl_e X` orthogonal to the balanced class `γ̂` has Iitaka dimension 0, for every `k` and both forms; so
   defect-0 pencils do not exist and mechanism (a) is void. Any pencil has defect `s ≥ 1` (search minimum 4), so its
   fibres through near-tight points have height a fixed positive proportion of the point's height.
3. **Restricted theorem (§3.2, Prop 3.5, proved).** On balanced profiles that are *not* simple, one-section Liouville
   beats `15/4` by the **chain deficit** `Δ_S`; near-tight points with `Δ_S ≥ (s_S + c)w` lie on `{f_S = 0}`. The
   record's caps are proved only for simple profiles, so this is consistent and new; it localises the hard core of
   Question B to small-chain-deficit (e.g. squarefree) profiles, where all fakes live, so it gives no growth.
4. **(b), (c) reduced (§3.3–3.4).** Fibrewise section arguments are global single-section arguments with a
   direction-sensitive archimedean bound (Prop 3.7); the projection from `e` is exactly tight; genus `≥ 1` fibres put
   the point in the small-point regime of its fibre, where no uniform Siegel theorem is known. The four Pell-type
   conditions decouple completely once the 11 non-singleton blocks are fixed (Prop 3.8): no transference exists, and
   (c) is counting below one point per scale.
5. **Tests.** T1–T5 are reported in §4; all pass (the mechanism proves nothing at `k = 4`, nothing for the split form,
   and excludes none of the calibration families).

## 1. The elementary Gaussian Cox form (verification of the coordinator's derivation)

**Conventions.** *Proved* = proved here with all quantifiers. *Exact* = finite exact integer computation in
`stab_checks/` (Python integers, `fractions`, python-flint `fmpz_mat.rank`, sympy). *Evidence* = computation that no
proof uses. *Heuristic* = no proof. Nothing here is refereed. Gaussian integers are written `z = x + iy`; `N(z) = |z|^2`.
Fix once and for all a Gaussian prime `π_p` above each `p ≡ 1 (4)` (in the scripts: `π_p = a + bi`, `a > b > 0`).

**Index set.** A cut of `[5]` is named by its side `A` not containing `5`; so the 15 cuts are the nonempty
`A ⊆ {1,2,3,4}`. For `a ∈ [4]`, the cuts separating `a` from `5` are exactly the `A ∋ a` (8 of them); the cuts
separating `a, b ∈ [4]` are the `A` containing exactly one of them (8 of them).

### 1.1 The parametrisation

Let `ζ_1..ζ_5` be a primitive tuple of lattice points of `|ζ|^2 = N'` **in one unit class**, i.e.
`ζ_j = ε ∏_p π_p^{a_j(p)} π̄_p^{e_p - a_j(p)}` with one unit `ε` for all `j` (`gauss.md` §1.2; `N'` odd, all
`p | N'` split, `min_j a_j(p) = 0`, `max_j a_j(p) = e_p` by primitivity). For a layer `(p,t)`, `1 ≤ t ≤ e_p`, let
`S = {j : a_j(p) ≥ t}`. If `5 ∉ S`, put the factor `π_p` into `γ_S`; if `5 ∈ S`, put `π̄_p` into `γ_{[5]∖S}`.

**Chain condition (Ch).** For every `p`: the sets `A` with `π_p | γ_A` form a chain (with multiplicity), the sets `A`
with `π̄_p | γ_A` form a chain, and the union of the first chain is disjoint from the union of the second.

> **Proposition 1.1 (proved; exact check `cox_form.py`).**
> (i) Every primitive one-class 5-tuple is
> `ζ_a = ε ∏_{A∋a} γ_A ∏_{A∌a} γ̄_A` (`a ≤ 4`), `ζ_5 = ε ∏_A γ̄_A`, with `γ_A` built as above; the `γ_A` satisfy
> (Ch), are products of the fixed `π_p, π̄_p` only, and `N(γ_A) = n_A` is the block norm of the cut (so
> `h_{E_A+E_{A^c}}([ζ]) = log N(γ_A)`, `R''^2 = N' = ∏_A N(γ_A)`).
> (ii) Conversely, if `γ_A` are products of the `π_p, π̄_p` satisfying (Ch), the tuple so defined is primitive, its
> layer decomposition returns exactly these `γ_A`, and its block norms are `N(γ_A)`. So, for fixed `ε`, the map
> `(γ_A) ↦ ζ` is a **bijection** from (Ch)-tuples onto primitive one-class tuples with unit `ε`.
> (iii) **Correction to the handoff form.** Without (Ch) the formula can still produce a *primitive* tuple, but then
> its true block norms are not the `N(γ_A)` (example `γ_{13}, γ_3, γ_{34}` each carrying one `π̄_5`: the tuple is
> primitive, its canonical blocks are `π̄_5^2` in `γ'_3` and `π̄_5` in `γ'_{134}`). Every statement about block norms,
> balance or level varieties must therefore be made on the (Ch)-locus of `Z[i]^{15}`.
>
> *Proof.* (i) `ζ_j = ε ∏_{(p,t)} (π_p if a_j(p) ≥ t else π̄_p)`; regroup the layers by cut. For a layer with
> `5 ∉ S = A`, the points `j ∈ A` receive `π_p` and the others (including `5`) `π̄_p`, which is the contribution of
> `γ_A ∋ π_p`; for a layer with `5 ∈ S`, `A = S^c`, the points in `S` (including 5) receive `π_p = conj(π̄_p)` and
> those in `A` receive `π̄_p`, which is the contribution of `γ_A ∋ π̄_p`. The sets `S_t = {a_j ≥ t}` decrease in `t`;
> those not containing 5 form the `π`-chain, the complements of those containing 5 form the `π̄`-chain, and since
> `S_{t'} ⊆ S_{t_0}` for the layers above the last one containing 5, the two unions are disjoint. Norms: each layer
> contributes `p` to exactly one `N(γ_A)`. (ii) At `π = π_p`: `v_π(ζ_a) = #{π-chain sets ∋ a} + #{π̄-chain sets ∌ a}`
> (`a ≤ 4`) and `v_π(ζ_5) = #π̄-chain`. If the `π̄`-chain is empty, `v_π(ζ_5) = 0`; otherwise pick `a` in its
> smallest member: `a` lies in every `π̄`-set and, by disjointness, in no `π`-set, so `v_π(ζ_a) = 0`. Symmetrically
> at `π̄`. Rational primes `q ≡ 3 (4)` and `1+i` do not occur. The layer sets of the resulting profile are, by the
> same count, the `S_t` read off from the two chains, so the construction of (i) returns the `γ_A`. (iii) is the
> displayed example (it occurs in `dbg_conv.py`). `□`

*Exact check.* `cox_form.out`: 300 random genuine profiles (3–9 primes up to 5000, exponents up to 3) satisfy every
identity of §§1.1–1.3 (the 24 flagged rows are random coincidences `ζ_a = ζ_b`, where every identity still holds with
`t_ab = 0`); 3000/3000 random `γ`-vectors satisfy "(Ch) ⟺ primitive and canonical".

### 1.2 Pair cofactors, residues, and the corrected triangle relation

For `a ≠ b` in `[5]` put
`Q_ab := ∏_{A ∋ a, A ∌ b} γ_A · ∏_{A ∋ b, A ∌ a} γ̄_A` (8 factors; `N(Q_ab) = e^{d_ab}`), and `Π_a := Q_{a5} = ∏_{A∋a} γ_A`.

> **Proposition 1.2 (proved; exact).** Under (Ch):
> (i) `ζ_a = G'_ab Q_ab`, `ζ_b = G'_ab Q̄_ab` with `G'_ab = ε ∏_{A ⊇ \{a,b\}} γ_A ∏_{A ∩ \{a,b\} = ∅} γ̄_A`, and `Q_ab` is
> conjugate-primitive; hence `ζ_a/ζ_b = Q_ab/Q̄_ab`, `gcd(ζ_a,ζ_b) = G'_ab`, the reduced chord is
> `c_ab = Q_ab - Q̄_ab = 2i·Im Q_ab`, and the record's residue is `t_ab = ±Im Q_ab ∈ Z`, with
> `|ζ_a - ζ_b| = 2|t_ab| R'' e^{-d_ab/2}`. In particular `t_{a5} = Im Π_a` and `|ζ_a - ζ_5| = 2R''|Im Π_a|/|Π_a|`.
> (ii) (**Triangle relations.**) For distinct `a, b, c ∈ [5]`:
> `Q_ab Q_bc = G^{(b)}_{abc} Q_ac`, `G^{(b)}_{abc} := ∏_{A separates b from \{a,c\}} n_A ∈ Z_{>0}` (4 cuts).
> For the triple `(a, 5, b)` (`Q_{5b} = Π̄_b`): `Π_a Π̄_b = G_ab Q_ab`, `G_ab := ∏_{A ⊇ \{a,b\}} n_A`.
> (iii) Writing `Π_a = x_a + i t_a`, `Q_ab = y_ab + i t_ab`, `M_{a∖b} := ∏_{A∋a, A∌b} n_A`:
> `x_b t_a - x_a t_b = G_ab t_ab`, `x_a x_b + t_a t_b = G_ab y_ab`, and the **binary quadratic form identity**
> `M_{b∖a} t_a^2 - 2 y_ab t_a t_b + M_{a∖b} t_b^2 = G_ab t_ab^2` (a definite form of discriminant `-4t_ab^2`).
>
> *Correction.* The handoff's "`c_ab c_bc = G_abc c_ac`" is false for the reduced chords themselves
> (`c_ab = unit·2i t_ab` are `2i` times small integers); the relation holds for the **pair cofactors** `Q_ab`, whose
> imaginary parts are the residues.
>
> *Proof.* (i) Split the 15 cuts by their intersection with `{a,b}`. Conjugate-primitivity of `Q_ab`: `π | Q_ab` needs
> (1) a `π`-chain set containing `a` but not `b`, or (2) a `π̄`-chain set containing `b` but not `a`; `π̄ | Q_ab` needs
> (3) a `π̄`-chain set containing `a` but not `b`, or (4) a `π`-chain set containing `b` but not `a`. The pairs
> (1)&(3) and (2)&(4) give intersecting members of the two chains; (1)&(4) and (2)&(3) give two non-nested members of
> one chain. (Ch) excludes all four, so `π` and `π̄` never both divide `Q_ab` (for `b = 5` only (1)&(3) is possible). Then `gcd(Q_ab, Q̄_ab) = 1` and the rest is
> `ζ_a - ζ_b = G'_ab(Q_ab - Q̄_ab)`. (ii) Compare, cut by cut, the exponents of `γ_A, γ̄_A` on both sides according to
> `A ∩ {a,b,c}` (`5` lies in no `A`): for `{a}` and `{a,b}` both sides carry `γ_A`; for `{c}` and `{b,c}` both carry
> `γ̄_A`; for `∅` and `{a,b,c}` neither does; for `{b}` (`γ̄_A` from `Q_ab`, `γ_A` from `Q_bc`) and for `{a,c}` (`γ_A` from
> `Q_ab`, `γ̄_A` from `Q_bc`) the left side carries `γ_Aγ̄_A = n_A` and the right side nothing; these are the 4 cuts
> separating `b` from `{a,c}`. (iii) Imaginary and real parts
> of `Π_aΠ̄_b = G_ab Q_ab`; the quadratic identity follows by eliminating `x_a, x_b` with `x^2 = N - t^2`,
> `y_ab^2 = M_{a∖b}M_{b∖a} - t_ab^2`, `N_a = G_ab M_{a∖b}`. `□`

### 1.3 Balanced scaling and the exact place of "stability" among the record's statements

Balanced at scale `X = e^w`: `n_A ∈ [X^{1-ε}, X^{1+ε}]`; then `|γ_A| ≈ X^{1/2}`, `|Π_a| ≈ X^4`, `R'' ≈ X^{15/2}`, and
`|ζ_a - ζ_5| = 2|t_a|X^{7/2 ± O(ε)}`. An A-point at `δ` (`diam ≤ R''^{1/2-δ}`) has `|t_a| ≤ X^{1/4 - 15δ/2 + O(ε)}`;
Liouville is `|t_a| ≥ 1`. Define

**Stab(η, ε).** There are a proper Zariski-closed `Z ⊂ T_5` and `w_0` such that every primitive one-class 5-tuple that
is `ε`-balanced at scale `w ≥ w_0` and has `max_a |t_{a5}| ≤ e^{ηw}` has `[ζ] ∈ Z`.

> **Proposition 1.3 (proved).** (i) `|t_ab| ≤ 2e^{(η + 8ε)w}` for all pairs, under the hypotheses of Stab.
> (ii) `(L_5^∘)(c, η_0, Z) ⇒ Stab(η, ε)` (same `Z`) whenever `ε ≤ η_0` and `10η + 48ε < c`.
> (iii) `Stab(η, ε) ⇒ (L_5^∘)` with `c = η`, `η_0 = ε`. Hence "`Stab(η,ε)` for some `η, ε > 0`" is **equivalent** to
> `(L_5^∘)`, and it implies the uniform bound (gauss Thm 3.4′).
> (iv) `A_5(1/30 - η) ⇒ Stab(15η/2 - Cε, ε)` for a constant `C`; the converse is not claimed (`A_5` also covers
> unbalanced profiles and all unit classes). So the handoff's "A_5(1/30-η) says: near-tight ⇒ Z" is correct as an
> implication; the equivalence is with `(L_5^∘)`, not with `A_5`.
>
> *Proof.* (i) `t_ab = (x_b t_a - x_a t_b)/G_ab`, `|x_a| ≤ X^{4(1+ε)}`, `G_ab ≥ X^{4(1-ε)}`. (ii) `H_Σ ≤ 4ηw +
> 6(η+8ε)w + O(1) < cw`. (iii) If `H_Σ < ηw` then, as every `|t_ab| ≥ 1`, each `log|t_{a5}| ≤ H_Σ < ηw`. (iv) Insert
> the A-point bound into the displayed formula for `|ζ_a - ζ_5|`. `□`

### 1.4 The split twin (verified on the Pell curve)

For a primitive positive integer 5-tuple `u` (base point `u_5`), the same layer construction with `v_p` in place of
`v_π` gives coprime-at-each-prime integer pairs `(u_A, v_A)` with `u_a = ∏_{A∋a}u_A ∏_{A∌a}v_A`, `u_5 = ∏_A v_A`,
`u_a/u_5 = U_a/V_a` (`U_a = ∏_{A∋a}u_A`, `V_a = ∏_{A∋a}v_A`) and `dist = |U_a - V_a|/V_a` (`γ_A ↔ u_A`,
`γ̄_A ↔ v_A`, `Im Π_a ↔ (U_a - V_a)/2i`). Exact check `split_pell.out` on the round-23 curve `(1,2,4)`,
`z^2 = 56x^2 + 65`, at the first seven Pell points (`x` up to `2.7·10^{10}`): the representation holds and the split
residues `U_a - V_a` are bounded and periodic (`[8,5,9,5]`, `[1,5,1,1]`, `[8,1,1,1]`, `[1,1,9,5]`), with
`α = log H/(-log dist)` decreasing `6.51, 4.67, 4.37, 4.20, 4.10, 4.01, 3.99` towards `15/4`.

### 1.5 Three equivalent elementary forms (proved reformulations; new to the record)

1. **Venn form.** If the blocks are pairwise coprime, `γ_A` is the part of `∏ P_a` (with `P_a := Π_a = x_a + i t_a`)
   that divides exactly the `P_a` with `a ∈ A`; so Question B is about **four conjugate-primitive Gaussian integers
   `x_a + i t_a` with bounded imaginary parts whose 15-region "Venn diagram" of prime factors is balanced**. In general
   (chains) the regions are read from the signed profiles `ε_a(p) = v_π(P_a) - v_{π̄}(P_a)` together with `ε_5 = 0`.
2. **Arccot form** (`t_a = t_ab = 1`). The near-tight condition reads `(x_b - x_a) | (x_a^2 + 1)` for all pairs, i.e.
   `cot(θ_a - θ_b) ∈ Z` for the five angles `θ_a = arccot x_a`, `θ_5 = 0`: a 4-term Machin/Lehmer-type arccotangent
   clique with balanced factor structure. (Proof: `(x_a+i)(x_b-i) = (x_ax_b+1) + i(x_b - x_a)`, and `d | x_ax_b+1` with
   `d = x_b - x_a` iff `d | x_a^2+1`.) In general `cot(θ_a - θ_b) = y_ab/t_ab` with all ten `t`'s small.
3. **Residue (singleton-decoupled) form.** Put `R_a := ∏_{A ∋ a, |A| ≥ 2} γ_A` (7 factors, `N(R_a) ≈ X^7`). Given the 11
   non-singleton `γ_A`, the condition `Im(γ_a R_a) = t_a` is a linear equation for `γ_a ∈ Z^2` whose solutions are
   spaced `|R_a| ≈ X^{7/2}` apart on a line at distance `|t_a|/|R_a|` from `0`; so **each singleton is determined by
   `(t_a, R_a)`**, and exists iff `x_a :≡ -ι(R_a) t_a (mod N(R_a))` (`ι(R_a)` the square root of `-1` with
   `R_a | ι(R_a) + i`) has a representative with `|x_a| ≲ X^4`; equivalently `‖t_a ι(R_a)/N(R_a)‖ ≲ X^{-3}`.
   The four conditions **decouple** completely once the 11 non-singletons are fixed (they are conditions on four
   different lattices `N(R_a)Z`); all coupling sits in the joint distribution of the four "Kloosterman fractions"
   `ι(R_a)/N(R_a) mod 1` as the 11 shared `γ_A` vary. Heuristic count: `X^{11}·(X^{-3})^4 = X^{-1}` per dyadic scale
   (sequentially `Π_1: X^4`, `Π_2: X^0`, `Π_3: X^{-2}`, `Π_4: X^{-3}`). Trivial bound for the number of
   near-tight (Ch)-points at scale `X` with `|t| ≤ T`: `≪ X^{11+o(1)}T^4`.

## 2. Level varieties on the Gaussian Cox space

Identify `Z[i]^{15} ⊂ A^{30}` (coordinates `Re γ_A, Im γ_A`). For `t ∈ Z^4` put
`V_t := {Im Π_a(γ) = t_a, a = 1..4}`.

* **Structure (proved).** `V_t` is cut out by four real forms of degree 8 (each `Im Π_a` involves 8 of the 15
  Gaussian variables), so it is a complete intersection, smooth of dimension 26 at every point where all
  `R_a = Π_a/γ_{\{a\}}` are nonzero: the singleton variable `γ_{\{a\}}` occurs only in `Π_a`, and the restriction of
  `d Im Π_a` to the `γ_{\{a\}}`-plane is the nonzero real form `Im(R_a·dγ_{\{a\}})`, so the Jacobian is block-triangular
  of rank 4. It is invariant under the anisotropic 11-dimensional torus
  `{γ_A ↦ u_A γ_A : |u_A| = 1, ∏_{A∋a} u_A = 1 ∀a}` (the norm-one part of the Néron–Severi torus fixing the four
  `Π_a`), and under the real scalings only up to `t ↦ (∏_{A∋a}λ_A) t_a`. The map `V_t ∩ (Ch) → T_5`,
  `γ ↦ [ζ]`, is injective on integral points up to the unit `ε` (Prop. 1.1(ii)).
* **Count (heuristic).** Integral points of `V_t` in the box `|γ_A| ≤ B = X^{1/2}`: `B^{30 - 4·8} = B^{-2} = X^{-1}`
  per dyadic scale, matching the record's `e^{-w}`. The projective closure of `V_t` (`t ≠ 0`) in `P^{30}` is a
  complete intersection of four octics, `K = O(32 - 31) = O(1)`: "general type" if its singularities are mild
  (**not checked**; label heuristic). Bombieri–Lang for `V_t` would not suffice anyway: the box condition, not
  rationality, is what is rare (`V_t(Q)` is Zariski dense, by the real scalings).
* **One `Z`, not one `Z_t`.** With `|t_a| ≤ X^{η}` there are `X^{4η}` level varieties at scale `X`; a statement of
  the form "`V_t ∩ box ∩ (Ch)` maps into `Z_t`" proves `Stab` only if `∪_t Z_t` is contained in one proper closed set.
  This is exactly the point where the mechanism of §3.1 would have needed a second step (Theorem 3.2(c)).

## 3. Main: stability mechanisms

### 3.1 Candidate (a), second-order integrality: it is a pencil argument, and pencils of tight sections do not exist

**The cofactor.** Let `S ⊂ M = {m ∈ Z^5 : Σm = 0}` be finite, `f = Σ_{m∈S} c_m x^m` with `c_m ∈ Z[i]`, and `y = [ζ]`
a primitive one-class tuple with profile `a(p) ∈ Z^5`. Put `m_0 = (min_{m∈S} m_j)_j`, `F(ζ) = Σ_m c_m ζ^{m - m_0} ∈ Z[i]`
(a polynomial in the `ζ_j`, hence in the Cox coordinates `γ_A, γ̄_A`), and let `D_S(ζ)` be the Gaussian integer
`∏_p π_p^{min_S⟨m,a(p)⟩ - ⟨m_0,a(p)⟩} π̄_p^{-max_S⟨m,a(p)⟩ - ⟨m_0, e_p 1 - a(p)⟩}` (all exponents `≥ 0`).

> **Lemma 3.1 (cofactor; proved).** `D_S(ζ) | F(ζ)`, and `c_f(y) := F(ζ)/D_S(ζ) ∈ Z[i]` is `Z[i]`-linear in `f` (for
> fixed `S, y`) and satisfies `|c_f(y)| = |f(y)|·H_S(y)`, `log H_S(y) := ½ Σ_p width_S(a(p)) log p`, where
> `width_S(v) = max_S⟨m,v⟩ - min_S⟨m,v⟩`. If `ord_e f ≥ t` then `|c_f(y)| ≤ C_f dist(y,e)^t H_S(y)` for `dist ≤ 1`.
> On a simple profile `log H_S(y) = ½Σ_{cuts} width_S(1_A) log n_A`; on any profile it is at most that
> (`width_S` is subadditive along the layer chain of each prime).
>
> *Proof.* `v_π(ζ^{m-m_0}) = ⟨m - m_0, a(p)⟩` and `v_{π̄}(ζ^{m-m_0}) = ⟨m - m_0, e_p1 - a(p)⟩ = -⟨m,a(p)⟩ - ⟨m_0, e_p 1 - a(p)⟩`
> (as `Σm = 0`); take minima over `m`. Since `|ζ_j| = R''` for all `j`, `|F(ζ)| = |f(y)| R''^{-Σ m_0}`, and
> `|D_S(ζ)| = R''^{-Σm_0} ∏_p p^{-width_S(a(p))/2}`. The Taylor bound is the record's (TL) archimedean side. `□`

This is the whole content of "second-order integrality": Liouville is `|c_f(y)| ≥ 1`; at a near-tight balanced point
every MW-tight `f` (`8t = Σ_{cuts} width_S`) has `|c_f(y)| ≤ C_f e^{O((η+ε)w)}`, i.e. `c_f(y)` is a Gaussian integer of
size `X^{O(η+ε)}` — the "O(1) integers" of the handoff. A *further* integer forced to be `< 1` can be produced in two
ways only: (I) a single section whose bound is `< 1` (impossible beyond the cap, by MW on simple profiles; see
Prop. 3.5 for non-simple ones), or (II) a combination whose coefficients depend on `y` through these small integers:
`g_y = c_{f_2}(y) f_1 - c_{f_1}(y) f_2`, which vanishes at `y` identically. (II) is a genuinely different mechanism
(not a Liouville inequality), so it is worth stating exactly what it would give.

> **Theorem 3.2 (pencil principle; proved).** Let `f_1, f_2 ∈ Z[i][M]` be linearly independent, supported in `S`,
> with `ord_e f_i ≥ t`, and let `s := Σ_{cuts} width_S(1_A) - 8t` (the *defect*; an integer). Let `y` be a primitive
> one-class 5-tuple, `ε`-balanced at scale `w`, with `max_a |t_{a5}| ≤ e^{ηw}`, and `(f_1(y), f_2(y)) ≠ (0,0)`.
> Then `φ(y) := [f_1(y) : f_2(y)] = [c_{f_1}(y) : c_{f_2}(y)] ∈ P^1(Q(i))` has
> `log H(φ(y)) ≤ (s/2)w + κ_S (η + ε) w + O_{f}(1)`, `κ_S = 4t + ½Σ width_S`.
> (a) If `s = 0`, then for every `T, K` the points with all `|t_ab| ≤ T` and `|log n_A - w| ≤ K` lie on the union of
> the finitely many hypersurfaces `{β f_1 = α f_2}`, `α, β ∈ Z[i]`, `|α|,|β| ≤ C(f,T,K)`, together with
> `{f_1 = f_2 = 0}`: a proper Zariski-closed `Z_{T,K}` — a "bounded-residue stability" statement.
> (b) If `s = 0` and only `max|t_{a5}| ≤ e^{ηw}`, `y` lies on one of `e^{O((η+ε)w)}` fibres of `φ`, each of height
> `O((η+ε)w)`, i.e. tiny compared with `h(y) = 15w`; a second step on these low-height fibres (where Roth/Schmidt-type
> theorems with polynomial dependence on the fibre do apply) would then be needed.
> (c) If `s ≥ 1`, the fibre of `φ` through `y` has height up to `(s/2)w`, a fixed positive proportion `s/30` of
> `h(y)`.
>
> *Proof.* `c_{f_i}(y) ∈ Z[i]` and `|c_{f_i}(y)| ≤ C dist^t H_S(y)`. Near-tight and balanced:
> `-log dist = λ(y) ≥ d_{a5}/2 - ηw - log 2 ≥ (4 - 4ε - η)w - 1` and `log H_S(y) ≤ ½(1+ε)w Σ width_S`, so
> `log|c_{f_i}| ≤ (s/2)w + κ_S(η+ε)w + O(1)`; the height of `[c_1 : c_2]` is at most `max log|c_i|`. (a) Under the
> bounded hypotheses the same computation gives `|c_{f_i}(y)| ≤ C(f,T,K)`, and `y ∈ {c_{f_2}(y) f_1 = c_{f_1}(y) f_2}`.
> Each such hypersurface is proper because `f_1, f_2` are independent. (b), (c) are the same bound. `□`

So the pencil principle would give a genuinely new kind of statement at `s = 0` — exactly the "exact algebraic
relation from second-order integrality" asked for. The obstruction is that `s = 0` pencils do not exist:

> **Theorem 3.3 (rigidity of mean-width extremal measures; proved, all `n`).** For finite `S ⊂ Z^n` let
> `MW_n(S) = 2^{1-n}Σ_{P⊆[n]} w_P(S)` (construct.md Thm B) and `A_t(S)` the space of measures on `S` (equivalently
> Laurent polynomials supported on `S`) of order `≥ t` at the identity. Then `dim A_t(S) ≤ 1` for every integer
> `t ≥ MW_n(S)`. In particular a measure attaining Theorem B (`ord = MW_n(S)`) is unique up to scalar on its support set.
>
> *Proof.* Induction on `n`; `n = 0` is trivial. For `n ≥ 1` slice by `x_1`: slices `S_c` (`c` the occurring values of
> `x_1`), `g_c := MW_{n-1}(S_c)`, `K := MW_n(S)`. Order the slices `S_{(1)}, …, S_{(r)}` by decreasing `g` (ties
> arbitrarily).
> *Step 1 (peeling for dimensions).* For any ordering and any `t`,
> `dim A_t(S) ≤ Σ_{i=1}^r dim A_{t-i+1}(S_{(i)})`. Indeed `μ ↦ (x_1 - c_{(1)})μ` maps `A_t(S)` into
> `A_{t-1}(S ∖ S_{(1)})` with kernel the order-`≥t` measures supported on the hyperplane `x_1 = c_{(1)}`, which is
> `A_t(S_{(1)})` (order computed in `Z^{n-1}`); iterate.
> *Step 2 (which terms can be nonzero).* Theorem B gives `reg(S_{(i)}) ≤ g_{(i)}`, and Theorem B's counting step gives
> `g_{(i)} + i - 1 ≤ K`. Hence `A_{t-i+1}(S_{(i)}) ≠ 0` forces `t - i + 1 ≤ g_{(i)} ≤ K - i + 1 ≤ t - i + 1`, i.e.
> `t = K` and `g_{(i)} = K - i + 1` (call such `i` *critical*); then by induction (`t - i + 1 = MW_{n-1}(S_{(i)})`)
> that term is `≤ 1`.
> *Step 3 (at most one critical index).* Theorem B's key inequality, keeping the `P' = ∅` term, gives for slices
> `c < c'`:
> `K ≥ (c' - c) + (g_c + g_{c'})/2 + 2^{1-n}[(c_max - c_min) - (c' - c)]`. (∗)
> (For `P' = ∅` the left term is `w_{\{1\}}(S) = c_max - c_min` and the right term `c' - c`; all other terms are as in
> Theorem B; `E_{P'} l_c(P') = g_c/2`.) Let `C_v = {c : g_c ≥ v}`. By (∗) without the last term, `C_v` lies in an
> interval of length `K - v`. If `j` is critical, `v = g_{(j)} = K - j + 1` and `|C_v| ≥ j`, so `C_v` consists of
> exactly `j` consecutive integers `[a, a+j-1]` (the first `j` slices), and (∗) at its two ends gives
> `g_a + g_{a+j-1} ≤ 2v`, so both ends have `g = g_{(j)}`. Suppose `i < j` are both critical. Then `g_{(i)} > g_{(j)}`,
> so `C_{g_{(i)}} = [b, b+i-1]` (the first `i` slices) lies strictly inside `(a, a+j-1)`. Apply (∗) to
> `(c,c') = (a, b+i-1)` and to `(b, a+j-1)`; in both cases `c' - c < (a+j-1) - a ≤ c_max - c_min`, so the last term
> is positive and
> `b - a < (j-i)/2` and `(a+j-1) - b < (i+j)/2 - 1`.
> Adding: `j - 1 < j - 1`, a contradiction. Hence at most one critical index, and `dim A_K(S) ≤ 1` by Steps 1–2. For
> `t > K`, `A_t ⊆ A_{⌈K⌉}`. `□`

*Exact checks (`rigidity.out`).* `dim A_{⌈MW⌉}(S) ≤ 1`, the critical-index count `≤ 1`, and Theorem B were checked on
all `4095` nonempty subsets of `{0..3}×{0..2}` and of `{0,1,2}×{0,1}^2`, all subsets of `{0,1,2}^2`, `{0,1}^3`, `{0..5}`,
4000 random sets in `Z^3`, 3000 in `Z^4`, 1500 zonotope sumsets of root binomials in the `k = 5` slice and 1500
random subsets of them: no violation; every tight set (`reg = MW`, 2 477 cases in total) has `dim A_MW = 1`.

*Relation to the record and novelty.* Theorem 3.3 is new to the record; I did not search the literature for it
(label: unrefereed, literature unchecked). construct.md's **Conjecture E** ("equality in Theorem B holds only for root-binomial
products") was supported by 1662 examples. Theorem 3.3 proves its **uniqueness half** for all `n` (if Conjecture E
holds, uniqueness follows from it, because the zonotope `Σ m_ij[e_i,e_j]` determines `m`; Theorem 3.3 does not need
it). The structure half (that the unique extremal measure is a root-binomial product) remains open.

> **Corollary 3.4 (the tight face of `Eff(Bl_e X)` is rigid; proved).** Let `X` be either form of the Losev–Manin
> 4-fold, `D` a torus-invariant divisor over `Q(i)`, `β` the balanced class (`β·E_ρ = 1`), and `t` an integer with
> `8t ≥ D·β`. Then `h^0(Bl_e X, π^*D - tE) ≤ 1`, and the same for every multiple. Hence:
> (i) every effective `L` on `Bl_e X` with `L·γ̂ = 0` (`γ̂` = balanced class of multiplicity 8 at `e`) has Iitaka
> dimension 0;
> (ii) **no pencil has defect `s = 0`**: Theorem 3.2(a),(b) never apply; every pencil of sections of order `≥ t`
> has `s ≥ 1`;
> (iii) the same holds with `f_1, f_2` of different supports but the same divisor class (a monomial correction
> `x^{m_0}` puts both in `H^0(π^*D_1 - tE)`), and a non-homogeneous Cox relation `c_2F_1 - c_1F_2 = 0` between
> sections of different classes does not descend to a Zariski condition on `T_5` (its zero set in the universal torsor is
> not invariant under the Néron–Severi torus, whose orbits are the fibres of the torsor, so it maps dominantly to `X`).
>
> *Proof.* Sections of `π^*D - tE` are Laurent polynomials supported in `P_D ∩ M` of order `≥ t`; and
> `8 MW(P_D ∩ M) = Σ_{cuts} width_{P_D ∩ M}(1_A) ≤ Σ_{cuts}(d_A + d_{A^c}) = D·β ≤ 8t`, so Theorem 3.3 (with `n = 4`
> after dropping the fifth coordinate, construct.md Cor B1) gives `dim ≤ 1`. (i) `L·γ̂ = D·β - 8t`. (ii) A pencil with
> `s = 0` would give two independent elements of `A_t(S)` with `t = MW(S)`. (iii) as stated. Over `Q` the
> `σ_X`-invariant subspace has dimension `≤ 1` as well. `□`

*Minimal defect (evidence, `defect.out`).* In a search over 2500 supports (zonotope sumsets plus 1–2 extra points) the
smallest defect of an actual pencil is `s = 4`, attained by the triangle `{0, e_1-e_2, e_1-e_3}` (`t = 1`,
`Σwidth = 12`; its pencil `φ = (x_1/x_2 - 1)/(x_1/x_3 - 1)` is the similarity shape of the sub-triangle `ζ_1ζ_2ζ_3`).
Its fibres are `pr_{123}^{-1}` of conics through `e` in `T_3`, on which `e` is approximated to exponent 2: useless.
`s ≥ 1` is proved (Cor. 3.4); `s_min = 4` is only evidence.

**Verdict on (a).** Second-order integrality at Liouville-tight points, in the only form available (Lemma 3.1 plus
point-dependent coefficients), is the pencil principle; by Theorem 3.3/Cor. 3.4 it is void at the top of the window,
for every `k`, for both forms. This is a **new barrier** (it is not a section bound: Theorem 3.2 is not an inequality,
and it would have produced exact relations), and it is closed by a **new theorem** (Theorem 3.3).

### 3.2 A restricted theorem: Liouville beats 15/4 on chain-cancelling profiles

Lemma 3.1 uses `width_S(a(p))` per prime, which is subadditive along the chain of layers of `p`. The record's caps
(round 24 theory item 2, round 25 intra §5) are proved only for **simple** profiles (one layer per prime), and
`round25/intra.md` lists non-simple balanced profiles as not covered. The cap indeed fails there:

> **Proposition 3.5 (chain cancellation; proved).** For finite `S ⊂ M` and an integer `t ≥ 1` with
> `A_t(S) ∋ f ≠ 0`, define the *chain deficit* of `y`:
> `Δ_S(y) := Σ_{cuts} width_S(1_A) log n_A - Σ_p width_S(a(p)) log p ≥ 0`.
> If `y` is `ε`-balanced at scale `w` with `max_a|t_{a5}| ≤ e^{ηw}` and `f(y) ≠ 0`, then
> `Δ_S(y) ≤ s·w + 2κ_S(η+ε)w + O_f(1)`, `s = Σ width_S - 8t`.
> Hence for any fixed finite family `{(S_i, f_i)}` and `c > 0`: near-tight balanced points with
> `Δ_{S_i}(y) ≥ (s_i + c)w` for some `i` lie in the fixed proper closed set `∪_i {f_i = 0}` once `η, ε` are small and
> `w` is large. *Example.* `S = {0, m}`, `m = e_1 - e_2 + e_3 - e_4`, `f = x_1x_3/(x_2x_4) - 1`, `s = 4`: a prime whose
> layer chain contains a cut with `⟨m,1_A⟩ = +1` nested with one with `⟨m,1_{A'}⟩ = -1` (e.g. `a(p) = (2,1,0,1,0)`,
> chain `{1,2,4} ⊃ {1}`) contributes `0` instead of `2 log p`; balanced profiles carrying more than `≈ 2w` of such
> cancellation are near-tight only on the subtorus `x_1x_3 = x_2x_4`.
>
> *Proof.* Lemma 3.1: `1 ≤ |c_f(y)| ≤ C_f e^{-tλ(y)} H_S(y)` with `2 log H_S(y) = Σ_{cuts} width_S(1_A) log n_A - Δ_S(y)
> ≤ (1+ε)w Σwidth_S - Δ_S(y)` and `λ(y) ≥ (4 - 4ε - η)w - 1`. `□`

*Scope (honest).* Non-simple 5-point profiles come only from primes at which the five points take at least three
distinct exponent values, so they need `e_p ≥ 2` in the cluster. Clusters with squarefree reduced norm produce only
simple 5-subtuple profiles, and all recorded fakes are squarefree. So Prop. 3.5 removes no part of the hard core; it
**locates** it: the 15/4 cap is a simple-profile (more precisely, small-chain-deficit) phenomenon, and Question B may
be restricted to profiles with `Δ_S ≤ (s_S + o(1))w` for every `S` in any fixed finite list.

*Exact checks (`cofactor.out`).* Lemma 3.1's divisibility and norm identity hold exactly in 316/316 random cases
(genuine profiles with up to 4 layers per prime, random supports and coefficients; chain deficits up to 42 seen).
On the split Pell points the deficits of all 120 segments `{0, e_a - e_b + e_c - e_d}` stay `≤ 6.8` while `w` grows to
12: the calibration family is not touched by Prop. 3.5 (T5).

### 3.3 Candidate (b): fibrations through `e`

Let `φ = [f_0 : … : f_n]`, `n ≥ 1`, be a rational map given by sections of order `≥ t` at `e` with common support `S`
(every rational map whose fibre closures pass through `e` is of this form after clearing denominators).

> **Proposition 3.6 (fibre heights; proved).** At an `ε`-balanced one-class point with `max|t_{a5}| ≤ e^{ηw}`,
> off the base locus, `h(φ(y)) ≤ (s/2)w + κ_S(η+ε)w + O(1)` with `s = Σwidth_S - 8t ≥ 1` (Cor. 3.4). So
> section-defined fibrations can only certify fibres of height a positive proportion `≥ 1/30` (per unit of `s`) of
> `h_B(y) = 15w`; they never certify fibres of height `o(w)`.

> **Proposition 3.7 (fibrewise sections are global direction-sensitive sections; proved).** Let `φ = f_1/f_2`, and
> suppose an auxiliary function on the fibre `F_c = {f_1 = c f_2}` is written `g_c = Σ_{j≤n} c^j g_j` with `g_j`
> supported on `S'`. Put `H := Σ_j f_1^j f_2^{n-j} g_j` (support `nS + S'`). Then `H|_{F_c} = f_2^n g_c|_{F_c}`, and
> for `y ∈ F_c(Q)`, `c = c_1/c_2` in lowest terms, the fibrewise Liouville inequality
> `1 ≤ |c_2^n g_c(y)|·H_{S'}(y)·(forced part)` is the single-section Liouville inequality for `H` at `y`, with the
> archimedean bound `|H(y)| ≤ |f_2(y)|^n·|g_c(y)|` in place of `C·dist^{ord_e H}`.
>
> *Proof.* Substitute `f_1 = c f_2` on `F_c`; `c_H(y) = c_{f_2}(y)^n·(c_2^n g_c(y))/D_{S'}(y)` up to units by Lemma 3.1
> (the forced divisor of `H` is `D_S^n D_{S'}`). `□`

So fibrewise **section** arguments are exactly global Liouville arguments with a **direction-sensitive** Taylor
bound — `|H(y)|` is small because `y` approaches `e` inside the special fibre `F_{φ(y)}` on which `H` vanishes to
higher order than on `X`. This class is explicitly *not* covered by the record's caps (round 25, intra item 3), so it
is a legitimate place to look; but every instance computed is exactly tight. *Example (projection from `e`, proved):*
the map `y ↦ [t_a ∏_{A∌a} γ̄_A]_{a ≤ 4} ∈ P^3(Q(i))` (the direction of `y - e` in the `x_a/x_5` chart, height
`≈ (7/2)w`) has fibres the lines `x_a/x_5 = 1 + c_a σ` through `e`, and on the fibre through `y` the parameter is
exactly `σ = 2i/Γ̄`, `Γ = ∏_A γ_A` (`|Γ| = R''`; proof: `c_aσ = Π_a/Π̄_a - 1 = 2it_a/Π̄_a`). So `dist = |c||σ| ≈
X^{7/2}/X^{15/2}`: the fibre (height `≈ (7/2)w`) plus the trivial bound `|σ| ≳ 1/|Γ|` on it reproduce exactly the pair
bound `|t_{a5}| ≥ 1`, with nothing to spare. Genus `≥ 1` fibres replace sections by Siegel/Faltings, which are not section arguments, but
Prop. 3.6 says the point `y` has height at most `30/s` times that of its fibre, i.e. `y` is a **small point** of its
fibre. In that regime no theorem with constants uniform in the fibre is known to me (not literature-checked this round; Siegel/Roth constants are ineffective
in the height of the target, and the elliptic-logarithm bounds depend on `h(E)` and on `e^{rank}`: the record's "Lim");
genus `≥ 2` fibres give finitely many points per fibre (Faltings) or `c(g)^{1+rank}` (uniform Mordell, recalled from
memory: Dimitrov–Gao–Habegger, Kühne), neither of which bounds the number of fibres. **Verdict on (b):** reduced to
either the direction-sensitive single-section class (open, every computed instance exactly tight) or a uniform
small-point Siegel theorem (not available). No gain.

### 3.4 Candidate (c): geometry of numbers on the four Pell-type conditions

> **Proposition 3.8 (decoupling; proved).** Fix the 11 non-singleton `γ_A` (with (Ch)), and `t ∈ (Z∖0)^4`. The set of
> singleton 4-tuples `(γ_{\{1\}},…,γ_{\{4\}})` completing them to a point of `V_t` is the product of four sets, the
> `a`-th being the lattice line `{γ : Im(γ R_a) = t_a}` (direction `R̄_a`, spacing `|R_a|`, at distance `|t_a|/|R_a|`
> from `0`); it meets the disc `|γ| ≤ X^{(1+ε)/2}` in at most one point once `X^{1/2+ε} < |R_a|/2`, and in one point
> iff the least residue of `-ι(R_a)t_a mod N(R_a)` is `≲ X^4`.

So there is **no transference between the four conditions at fixed non-singletons**: they are four independent
one-dimensional problems on four different lattices. All coupling is in how the four moduli `N(R_a)` and residues
`ι(R_a)` vary jointly with the 11 shared `γ_A` (`R_a ∩ R_b` = 4 shared factors, triple intersections 2, all four 1).
A geometry-of-numbers or equidistribution argument for `(ι(R_a)/N(R_a))_a mod 1` at scale `X^{-3}` would be a counting
statement whose main term is `X^{-1}` per scale (`< 1`); it can give upper bounds of the form `X^{11-θ}` (bilinear
Kloosterman-fraction estimates; not attempted) but cannot by itself place the survivors in one closed set: that
requires a structure theorem for the configurations where equidistribution fails, which is the statement itself.
**Verdict on (c):** the elementary content (one-point-per-line, the residue form) is verified; as a stability
mechanism it reduces to counting below one point per scale, where no counting method can conclude.

### 3.5 What the four candidates leave

| candidate | outcome | type of result |
| --- | --- | --- |
| (a) second-order integrality | equals the pencil principle (Thm 3.2); void because MW-extremal sections are unique (Thm 3.3, Cor 3.4) | **new theorem + new barrier** |
| (a′) chain cancellation | Liouville beats 15/4 on balanced profiles with large chain deficit (Prop 3.5) | **restricted theorem** (non-simple profiles only) |
| (b) fibrations | fibres have height `≥ (s/2)w`, `s ≥ 1`; fibrewise sections = global direction-sensitive sections (Prop 3.6–3.7) | reduction; no gain |
| (c) geometry of numbers | complete decoupling at fixed non-singletons (Prop 3.8); coupling only through the shared factors | reduction to counting `< 1` per scale |
| (d) new forms | Venn form, arccot cliques, residue form (§1.5) | reformulations for experts |

## 4. Pre-registered tests T1-T5

* **T1 (fakes; growth claims).** No growth or uniformity claim is made. Theorem 3.3 is a statement about measures, not
  about points. Prop 3.5 is a one-section Liouville inequality; on simple profiles (all recorded fakes are squarefree,
  hence simple) it is the MW-capped inequality and excludes nothing, so fakes are untouched and no growth follows.
* **T2 (k = 4 negative control; `k4_check.out`).** On gauss Thm 4.1's instance (`X = -3, m = 3`): `Im Π_x(s) ≡ -420,
  -180, -240` exactly; along `s` up to `1.2·10^{20}` the reduced 4-tuples have Gaussian gcd norm `≤ 10368000`,
  `|c_ab|^2 ≤ 14400` (bounded residues), `log N/log s → 14`; the singleton decoupling (Prop 3.8 at `k = 4`) holds.
  (a) Theorem 3.3 holds for all `n`, in particular `n = 3` (`k = 4`): the tight supports of `V_4`, `V_4^2`,
  `(12)(34)`, `(12)(13)(23)` (all lattice points of their Newton polytopes) have `dim A_MW = 1`. So the pencil principle
  is void at `k = 4` as at `k = 5`: the mechanism proves nothing at `k = 4`, as it must. (Theorem 3.2(a) itself is not
  contradicted by the dense family in any case: its `Z_{T,K}` depends on `T`, and the family's residue bound grows
  with the height of `(X, m)`.) (a′) Prop 3.5 at `k = 4`: the family's blocks are values of distinct linear
  polynomials, so repeated primes come only from the bounded resultants and the chain deficit is `O(1)`: nothing is
  excluded. (b), (c) are reductions and exclude nothing.
* **T3 (split test).** Theorem 3.3 and Cor 3.4 are form-independent (both forms have the same Laurent sections and
  the same MW), so the pencil route is void for split Question B as well; nothing here proves split Question B, so no
  claim about Erdős #887 arises. Prop 3.5 has a split analogue (same proof with `v_p`); it is a restricted statement
  (large chain deficit) and does not touch #887's construction: on the `(1,2,4)` Pell points the chain deficits stay
  `≤ 6.8` (`cofactor.out`), and the construction's residues are bounded and periodic (`split_pell.out`). It lands in
  the exceptional hypersurface `W` of the record, which nothing here shrinks.
* **T4 (criticality).** Theorem 3.2 is *not* a Liouville inequality and not a section bound: it uses two sections and
  point-dependent integer coefficients and outputs an exact relation, so Lemma K1, the `τ_e = 10/3` cap, the class-`L`
  cap and the multi-point reduction do not apply to it. It fails for a different, new reason: Theorem 3.3 (no two
  independent tight sections on one support). Prop 3.5 *is* a one-section Liouville argument; it is consistent with
  the class-`L` cap because that cap is proved for simple profiles only, and it gains exactly the chain deficit that
  the cap's proof discards. Prop 3.7 places fibrewise section arguments in the direction-sensitive class, which none
  of the four caps covers.
* **T5 (calibration).** Nothing excludes the round-24 real family `M` (no `Q`-points known; Theorem 3.2 never
  applies; Prop 3.5 needs chain deficit, absent on simple profiles), the split Pell curves (bounded residues, chain
  deficit `≤ 6.8`), or the record's 8- and 10-point families (no growth claim; Prop 3.5 concerns near-tight balanced
  5-tuples only).
* **Handoff constraints.** No fixed-`S` Subspace exceptional set is used. Lemma 3.1 keeps the full Gaussian gcd
  (`R''^{-Σm_0}` and the forced divisor are tracked exactly; exact check). No phase is replaced by zero: the pencil
  value is a function on `T_5`, and Cor 3.4(iii) records that non-homogeneous Cox relations do not descend. Template
  dependence: Theorem 3.2(a)'s `Z_{T,K}` and Prop 3.5's thresholds depend on `T` resp. on the fixed family, and are
  stated as such; nothing is claimed uniform in the radius that is not.

## 5. What this round adds to the barrier map

1. **Elementary form verified, with two corrections** (Prop 1.1–1.2): the Cox parametrisation is a bijection only on
   the chain-condition locus (off it, primitive tuples arise with wrong block norms); the triangle relation holds for
   the pair cofactors `Q_ab` (imaginary part `t_ab`), not for the reduced chords. Equivalence: `Stab ⟺ (L_5^∘)`
   (Prop 1.3); `A_5(1/30-η) ⇒ Stab` only.
2. **New theorem (Thm 3.3):** mean-width-extremal measures are unique on their support set, in every dimension (the
   uniqueness half of construct.md's Conjecture E). Corollary: the face `γ̂^⊥` of `Eff(Bl_e X)` consists of divisors
   of Iitaka dimension 0, for every `k` and both forms.
3. **New barrier:** "second-order integrality at Liouville-tight points" = the pencil principle (Thm 3.2), which
   would give exact relations (bounded-residue stability) if a defect-0 pencil existed; Thm 3.3 shows none exists. Any
   section-defined pencil has defect `s ≥ 1` (search minimum 4), so its fibres through near-tight points have height a
   positive proportion of the point's.
4. **Restricted theorem (Prop 3.5):** Liouville beats `15/4` on balanced profiles with large chain deficit; the hard
   core of Question B is the small-chain-deficit (in particular simple, e.g. squarefree) profiles.
5. **Reductions:** fibrewise section arguments are direction-sensitive global sections (Prop 3.7; every computed
   instance exactly tight); genus `≥ 1` fibres need a uniform small-point Siegel theorem; the four Pell-type conditions
   decouple completely given the 11 non-singletons (Prop 3.8), so (c) is pure counting below one point per scale.
6. **New elementary forms** for experts (§1.5): four Gaussian integers near the real axis with a balanced 15-region
   Venn diagram; 4-term arccotangent cliques `(x_b - x_a) | x_a^2 + 1` with balanced factor structure.

## Appendix: scripts

All in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/stab_checks/`:
* `gauss_lib.py` — exact Gaussian integer arithmetic (gcd, primes, divisibility).
* `cox_form.py` → `cox_form.out` — Prop 1.1–1.2 on 300 genuine profiles; converse (Ch) test (3000/3000); unit mixing
  (1020/1999 pairs non-integral `t_ab`, reproducing kill_gauss minor 2). `dbg_conv.py` — the (Ch) counterexample of
  Prop 1.1(iii).
* `rigidity.py` → `rigidity.out` — Theorem 3.3, Theorem B and the critical-index lemma, exhaustive and random.
* `defect.py` → `defect.out` — minimal pencil defect in a search (evidence: 4).
* `cofactor.py` → `cofactor.out` — Lemma 3.1 exact (316/316) and chain deficits; split Pell deficits.
* `split_pell.py` → `split_pell.out` — split twin on the `(1,2,4)` Pell curve.
* `k4_check.py` → `k4_check.out` — T2 on gauss Thm 4.1's instance.
Notes: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/stab_notes.md`.
