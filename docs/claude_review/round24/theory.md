# Round 24 (theory): does the anisotropy of T_5 give alpha^ess_e(-K) > 15/4?

*Round 24, theory lens. Status: complete draft, unrefereed. Scripts and outputs in `theory_checks/`.*

## 0. Summary and verdict

**Verdict: `reduces_to_split_cap`. No gain, no growth improvement, no uniform bound.** `M ≤ (2/3+ε) log R/loglog R`
stands. None of the four anisotropy candidates yields an inequality for twisted A-points that split points violate
and that beats the Liouville constant: each one reduces, by a proved statement, to toric Liouville for a single
`Q`-rational, centrally symmetric section (or to its `n`-point version), whose caps are the split caps
(`15/4` balanced, `δ = 1/10` arbitrary). The common reason (Theorem 7.1) is that Galois averaging is lossless:
`σ_X` fixes `e`, swaps `π ↔ π̄`, and preserves orders at `e`, so every argument equals its Galois average, which is a
symmetric toric Liouville inequality governed by the form-independent invariants `t_max(m) = ⌊10m/3⌋` and the mean
width. **The anisotropy changes which points exist, not which section inequalities hold.**

Results (all proved here unless labelled; nothing refereed):
1. **Galois structure (§2).** Lemma 2.1: twisted points have `v_{π̄} = -v_π`, no contact at `q ≡ 3 (4)` or `2`, and
   `|x| = 1` at `∞`; split points have `v_{π̄} = v_π`. Prop 2.2 (Hilbert 90): both forms are "five primitive vectors
   in `Z^2` modulo `SO(Φ_0)`" for the binary form `Φ_0 = ab` (split) or `a^2 + b^2` (twisted); height = product of the
   ranges of the signed contacts with the roots of `Φ_0`; distance = `|det|/√(Φ_0Φ_0)`. Exact check `h90.out`.
   Lemma 2.3: Picard ranks `15` (twisted) and `26` (split): only log powers in counts change (heuristic as input).
2. **Candidate 1, isotypic sections (§3).** Prop 3.1: for **every** `m`, the `ι`-isotypic parts of `H^0(m(-K))` have
   maximal orders `t_max(m)` and `t_max(m) - 1` (parity of `t_max` decides which), with explicit extremal products of
   `V`, `V_{abc}`, `D_{ab}`; the `σ_X`-rational parts have the same orders (Galois descent). Exact ranks for
   `m ≤ 5` agree (`isotypic.out`). Prop 3.2: reality of values of `Q`-rational sections and the pair `f, f^*` give
   exactly the width bound. Cor 3.3: reduces.
3. **Candidate 2, conjugate products and Thue (§4).** Prop 4.1: `(y, ιy)` is the graph of `ι`, so two-point arguments
   at `(y, ιy)` are single sections of `(d_1+d_2)(-K)`. Prop 4.2: at the rational target `(e,…,e)` the weighted index
   is additive over factors (filtered tensor products; exact check `index_add.out` on `V_1 ⊗ V_1` with weights `(1,1)`, `(1,2)`: maximal index `6`, `9`). Prop 4.3: the `n`-point Liouville
   inequality is implied by the one-point caps on balanced data with disjoint primes. Cor 4.4: reduces.
4. **Candidate 3, parity of point profiles (§5).** Prop 5.1: point parity *is* the Gaussian profile, so fakes satisfy
   it (Cor 5.2, record-conditional). Prop 5.3: the twisted local right side of toric Liouville is the **mean** of the
   split right sides at `v` and `-v`, so parity can only weaken Liouville (pair bound: `½(X+Y)` against `min(X,Y)`).
   The anisotropy enters only through curves (the parity lemma), never pointwise.
5. **Candidate 4, Hilbert 90 (§6).** Prop 6.1: the bracket divisibility `Π p^{c_ab} | [P_a,P_b]` is exactly the record's
   pair bound `|t_ab| ≥ 1`, with the exact formula `t_ab^2 = [P_a,P_b]^2/(Π p^{2c_ab}2^{e_a+e_b})` (also explaining
   the record's `√(odd/2)` across unit classes). Prop 6.2: every polynomial in `(P, P̄)` is a section times a unit
   (`G = Π P̄_j^{deg_j}·g(x)`, integer exponents, so no half-integer loophole in the mean-width bound). Prop 6.3:
   real-valued brackets give nothing. Plücker/Ptolemy relations fall under the abc barrier. Cor 6.4: reduces.
6. **Theorem 7.1.** The whole class `𝓛` (any isotypic sections, conjugates, group words, several points,
   Hilbert-90 polynomials, parity, reality) is capped at the split values; Lemma K1 applies to each member.

**What would be needed.** An input outside `𝓛` that is also not local: one that uses where A-points lie globally
(on which curves, how many), combined with the curve parity lemma. All such inputs in the record are conditional
(`V_E(5)`, McKinnon) or concern excluded families (bounded-degree covering curves). Question B is unchanged.

## 1. Setting, notation, what is used

**Labels.** *Proved* = proved here with all quantifiers (or cited with the place of a refereed/killer-checked proof).
*Exact* = finite exact computation in `theory_checks/` (integers, `Fraction`, sympy, or ranks modulo a prime used only
in the direction "full rank mod p implies full rank over Q"). *Record* = proved earlier in the record, cited.
*Heuristic* = no proof. Nothing in this file is refereed.

**1.1 The variety and the two forms.** `N = Z^5/Z·1`, `M = {s ∈ Z^5 : Σ s_j = 0}`, `x^s = Π x_j^{s_j}`. `LM_5` is the
smooth projective toric 4-fold of the permutohedral fan (rays `1_A`, `∅ ≠ A ⊊ [5]`, cones = chains). Its torus over a
field `F` is `G_m^5/G_m`. The **split form** is `LM_5` over `Q`, torus `T`, `T(Q) = (Q^*)^5/Q^*`. The **twisted form**
`X = Lbar_5^tw` is the `Q`-form of `LM_5 ⊗ Q(i)` whose Galois action is

    σ_X := ι ∘ c,     c = complex conjugation of coordinates,   ι(x) = x^{-1} (inversion, a toric automorphism),

so that `σ_X` acts by `-1` on `N` and `M` and sends the ray `1_A` to `1_{A^c}` (toric §1.1). Its torus is the
anisotropic norm-one torus `T_5`; `T_5(Q) = {x ∈ (Q(i)^*)^5/Q(i)^* : x̄ = x^{-1}}` (coordinates of norm one, modulo a
common norm-one factor). `e = [1:…:1]` is rational in both forms. `B = -K = Σ_A E_A` (nef, big, not ample; polytope
`R_5 = {s ∈ M_R : Σ|s_j| ≤ 2} = conv(e_a - e_b)`; toric Prop 1.1). `H^0(X_{Q(i)}, mB)` = Laurent polynomials with
support in `S_m = mR_5 ∩ M`.

**1.2 Heights, distance, A-points (record).** For `y = [ζ] ∈ T_5(Q)`, `ζ` a primitive tuple on `|ζ| = R''`:
`h_B(y) = log R''^2 = W'` and `h_{E_A + E_{A^c}}(y) = log n_A` exactly (`gauss.md` Lemma 1.1); `dist(y,e) =
max_j |ζ_j/ζ_1 - 1|`, `λ(y) = -log dist(y,e)`. An **A-point at exponent `γ`** is `y ∉ Z_diag` with
`dist(y,e)^γ H_B(y) ≤ 1`; `γ = 4/(1+2δ)`. `A_5(δ)` for some `δ < 1/30` ⟺ `α^ess_e(B; T_5(Q)∖Z_diag) > 15/4`
(toric Lemma 1.3). A point is **`η`-balanced at scale `w`** if `|log n_A - w| ≤ ηw` for all 15 cuts; then
`μ(y) := λ(y)/w ≤ 4 + O(η)` (balanced Liouville, toric Lemma 1.3(c)), and `(L_5^∘)` ⟺ `μ ≤ 4 - c'` off `Z` on balanced
one-unit-class points. For the split form: `H_B(u) = (max u/min u)·lcm(u)` for primitive positive `u`, ray heights
`h_{E_A}(u) = log n_A + O(1)` with `n_A`, `n_{A^c}` independent (split §1.4).

**1.3 The section caps used throughout (record, refereed).**
* **(TL) Toric Liouville** (toric Thm 3.2). If `0 ≠ f ∈ Z[i][M]` has Newton polytope `Q_f` and `f(y) ≠ 0`, then
  `ord_e(f)·λ(y) ≤ h_{D_{Q_f}}(y) + O_f(1) = ½ Σ_{cuts} width_{Q_f}(1_A)·log n_A + O_f(1)`.
* **(MW) Mean width** (`threshold/construct.md` Thm B, refereed): for finite `S ⊂ Z^5`,
  `reg(S) ≤ 2^{-4} Σ_{P ⊆ [5]} w_P(S)`; with toric Lemma 3.1, `ord_e f ≤ 2^{-3} Σ_{cuts} width_{Q_f}(1_A)`.
  Consequence (toric Thm 3.3, §4): on `η`-balanced points (TL) never gives better than `μ ≤ 4 + O(η)`, i.e. the
  balanced cap is the Liouville constant `15/4`.
* **(τ) Width theorem** (round 23 `width.md` Thm 2.5, refereed): `t_max(m) = ⌊10m/3⌋` for every `m`, split and twisted;
  hence for arbitrary profiles (TL) gives at best `δ > 1/10` (toric Thm 3.3 with `τ_e = 10/3`).
* **(K1)** Lemma K1 (`kill_gauss.md` §6.1, proved by the killer): bihomogeneous forms, forced divisors, angular width,
  pair data and residues are unchanged by dividing a tuple by its Gaussian gcd.
* **(Fake)** For balanced 5-row Gaussian profiles there are fake circles (real angles `φ_p` in place of `arg π_p`,
  genuine primes and genuine exponent profile) with diameter `≍ R''^{7/15}` that satisfy every monomial gap
  (`gauss.md` §5 item 1, `kill_gauss.md` §5 item 1, killer-checked; restated in toric §6.3). They violate `(L_5^∘)`.

**1.4 What "a gain" and "reduces" mean here.** A *gain* is a proved inequality, valid for all `y ∈ T_5(Q)` off a
proper closed `Z` fixed in advance, which on `η`-balanced points gives `μ ≤ 4 - c` for a constant `c > 0`
independent of `y`, `R''`, the primes and the profile (or which gives `A_5(δ)` for some `δ < 1/30`). A candidate
*reduces to the split section cap* if every inequality it produces is implied, up to `O(1)` or `O(η w)`, by (TL) for a
single section on `X` (or by its `n`-point version, Prop 4.3), so that (MW)/(τ) cap it at `15/4` (balanced) and
`10/3` (arbitrary profiles), the same values as on the split form; or if it is satisfied by the fakes of (Fake).

## 2. Galois structure of T_5: what the twist does to points, sections and contacts

Throughout, `K = Q(i)`, `p` runs over primes `≡ 1 (4)`, `p = π π̄` with a fixed choice of `π`, and `q` over primes
`≡ 3 (4)`. For `x ∈ T(K) = (K^*)^5/K^*` and a finite place `w` of `K` put `v_w(x) := (v_w(x_j))_j ∈ N`.

**Lemma 2.1 (local structure of twisted points; proved).** Let `y ∈ T_5(Q)`, with `K`-coordinates `x` (so `x̄ = x^{-1}`).
(a) `v_{π̄}(x) = -v_π(x)` for every split `p`.
(b) `v_w(x) = 0` for `w = q` inert and for `w = (1+i)`.
(c) `|x_j|_∞ = 1` for all `j`.
For a split point `u ∈ T(Q) ⊂ T(K)` instead: `v_{π̄}(u) = v_π(u) = v_p(u)`, `v_q(u)` is arbitrary, and `u` is real.
*Proof.* (a) `v_{π̄}(x_j) = v_π(x̄_j) = v_π(x_j^{-1})`. (b) For `w` stable under conjugation, `v_w(x_j) = v_w(x̄_j) =
-v_w(x_j)`. (c) `x_j x̄_j = 1`. The split statements: `u_j ∈ Q`, so `ū_j = u_j`. `□`

So at the places of `Q`: `T_5` is split at `p ≡ 1 (4)` and anisotropic (compact, no boundary contact) at `q`, `2`
and `∞`; the split torus is split everywhere. Lemma 2.1(a) is the point-level form of the curve parity lemma
(split §1.4): over `K`, `h_{E_A}(y) = h_{E_{A^c}}(y) = ½ log n_A` for every `A` (toric §1.1), whereas for split points
the ray heights `h_{E_A}`, `h_{E_{A^c}}` are independent.

**Two real structures on one complex picture.** Both forms sit inside `T(K)`: the split points as the fixed points of
`c` (real `x`), the twisted points as the fixed points of `ι∘c` (unitary `x`). Near `e`, writing `x = exp(y)`,
split points have `y ∈ R^4` and twisted points have `y ∈ iR^4`, two totally real slices of `C^4`. A holomorphic
function vanishes to order `t` at `0` iff its restriction to either slice does (a totally real `R`-subspace of maximal
dimension is a uniqueness set for Taylor coefficients). So **orders of vanishing at `e` do not see the form.**

**Proposition 2.2 (Hilbert 90: the binary-form model; proved, exact check `theory_checks/h90.py`).** Let
`Φ_0(a,b)` be a binary quadratic form over `Q` with nonzero discriminant, `Z_0 = {Φ_0 = 0} ⊂ P^1` (two geometric
points), and `SO(Φ_0)` the one-dimensional torus fixing `Z_0` pointwise. Then the Losev–Manin torus of the pair
`Z_0` is `T^{Φ_0} = (P^1∖Z_0)^5/SO(Φ_0)`; `Φ_0 = ab` gives the split form (`Z_0 = {0,∞}`, `u_j = a_j/b_j`) and
`Φ_0 = a^2 + b^2` gives `T_5` (`Z_0 = {±i}`, `x_j = P_j/P̄_j` with `P_j = a_j + i b_j`). Concretely, for `T_5`:
(i) every `y ∈ T_5(Q)` is `y = [P_j/P̄_j]` with `P_j ∈ Z[i]` **primitive** (no rational prime divides `P_j`); each
`P_j` is unique up to sign, and the tuple up to `P_j ↦ γP_j` (`γ ∈ K^*`, then re-primitivised);
(ii) with `[P,P'] := Im(P̄ P') = a b' - b a' ∈ Z`: `x_a/x_b - 1 = 2i[P_b,P_a]/(P̄_a P_b)`, so
`|x_a/x_b - 1| = 2|[P_a,P_b]|/(|P_a||P_b|) = 2|sin(arg P_a - arg P_b)|`;
(iii) with `ε_j(p) := v_π(P_j) - v_{π̄}(P_j)`: `N(P_j) = 2^{e_j} Π_p p^{|ε_j(p)|}` (`e_j ∈ {0,1}`) and
`H_B(y) = Π_p p^{max_j ε_j(p) - min_j ε_j(p)} = R''^2`;
(iv) (pair divisibility) `v_p([P_a,P_b]) ≥ c_ab(p) := max(0, min(ε_a,ε_b), min(-ε_a,-ε_b))`, and
`c_ab(p) - (|ε_a(p)| + |ε_b(p)|)/2 = -|ε_a(p) - ε_b(p)|/2`; hence
`|x_a/x_b - 1| ≥ 2^{1-(e_a+e_b)/2} e^{-d_ab/2}`, `d_ab = Σ_p |ε_a(p) - ε_b(p)| log p = Σ_{A sep a,b} log n_A`.
For the split form the same statements hold with `P_j = (a_j,b_j) ∈ Z^2` primitive, `[P,P'] = det`,
`|a_j b_j| = Π_p p^{|v_p(u_j)|}` over **all** primes, and the extra archimedean factor `max u/min u` in `H_B`.
*Proof.* (i) Hilbert 90: `x_j ∈ K` with `N(x_j) = 1` is `P/P̄` for some `P ∈ K^*`, unique up to `Q^*`
(`P/P̄ = P'/P̄'` iff `P'/P` is real); clearing denominators and rational content makes `P` primitive and unique up to
sign. The common factor of `T_5` is `γ/γ̄`. (ii) `P_aP̄_b - P̄_aP_b = 2i Im(P_aP̄_b) = 2i[P_b,P_a]`, divide by
`P̄_aP_b`, and `|P̄_a P_b| = |P_a||P_b|`. (iii) Primitivity excludes `q | P_j`, `2 | P_j` and `π π̄ | P_j`, so at each
`p` at most one of `v_π(P_j)`, `v_{π̄}(P_j)` is positive. `v_π(x_j) = v_π(P_j) - v_π(P̄_j) = ε_j(p)`, and
`v_{π̄}(x_j) = -ε_j(p)`. The toric height of `R_5` at `w` is `max_{a,b}|x_a/x_b|_w`; at `π` and at `π̄` it is
`p^{range_j ε_j}`, each with weight `½` in the absolute normalisation; it is `1` at `q`, `(1+i)` and `∞` (Lemma 2.1).
Since `ζ_j/ζ_1 = x_j/x_1` and `v_π(ζ_j/ζ_1) = a_j(p) - a_1(p)` (`gauss.md` §1.2), `range ε = range a`, so the product
is `R''^2` (`gauss.md` Lemma 1.1). (iv) Put `α_j = v_π(P_j)`, `β_j = v_{π̄}(P_j)` (one of them is `0`). Then
`v_π(P_aP̄_b) = α_a + β_b` and `v_π(P̄_aP_b) = β_a + α_b`, so `v_π([P_a,P_b]) ≥ min(α_a+β_b, β_a+α_b)`, and the
same number bounds `v_{π̄}` (swap `α`, `β`). Since `[P_a,P_b] ∈ Z`, `v_p = v_π`. A case check on the signs of
`ε_a, ε_b` gives the displayed value and the identity. Multiply over `p` and use `|P_j| = N(P_j)^{1/2}`. The
equality `Σ_p|ε_a - ε_b| log p = d_ab` is `gauss.md` §1.2 (the layers `(p,t)` separating `a,b` are the `t` between
`a_a(p)` and `a_b(p)`). The split statements are the same computation with `v_p(a_j)`, `v_p(b_j)` in place of
`α_j, β_j`. `□`

*Reading.* Hilbert 90 turns both problems into the same one: five primitive integer vectors in `Z^2` with nearly
equal directions, height = product of the ranges of the signed contacts with the two roots of `Φ_0`, distance =
`|det|/√|Φ_0Φ_0|`. Proposition 2.2(iv) is exactly the record's pair bound `|t_ab| ≥ 1` (`gauss.md` §1.2): the
"extra divisibility" of `[P_a,P_b]` is the pair gcd and nothing more. The only differences between `Φ_0 = ab` and
`Φ_0 = a^2+b^2` are (1) which primes can carry contact (all `p`, versus `p ≡ 1 (4)`), (2) the sign rule of
Lemma 2.1(a) (the two roots `0, ∞` are rational, so contact with `0` at `p` and contact with `0` at `p'` are recorded
in one rational integer `a_j`; the roots `±i` are conjugate, so contact with `i` at `π` and with `-i` at `π̄` are
recorded in one Gaussian integer `P_j` and its conjugate), (3) definiteness at `∞`. Sections 3–6 show that none of
(1)–(3) enters any inequality of the four candidate types beyond the split section cap.

**Lemma 2.3 (counting: only log powers change; proved for the ranks, the counting itself is heuristic).** The Picard
rank of the twisted form over `Q` is `15` (Galois orbits of the 30 rays, minus the rank `0` of the `Q`-character
group of the anisotropic `T_5`); that of the split form is `30 - 4 = 26`. In Manin–Batyrev–Tschinkel counts this changes
the power of `log B` (`14` against `25`), not the exponent of `B`. So the expected number of A-points at exponent `γ`
in a dyadic height range is `B^{1-4/γ+o(1)}` in both forms (**heuristic**), and the counting threshold stays at
`γ = 4`.
*Proof of the ranks.* For a smooth projective toric variety over `Q` with split torus over `K`, `Pic(X_K)` is the
cokernel of `M → Z^{rays}`; Galois-invariants give `rank Pic(X) = #(ray orbits) - rank M^{Gal}` (the `H^1` term is
finite). Rays: `σ(1_A) = 1_{A^c}` and `A ≠ A^c`, so `15` orbits; `M^{Gal} = {m : -m = m} = 0`. Split: `30` orbits,
`rank M = 4`. `□`

## 3. Candidate 1: σ-equivariant auxiliary sections (isotypic parts, exact orders)

`V_m := H^0(X_K, mB) = K[S_m]`. Two involutions act and commute: the **linear** involution `ι^*f(x) = f(x^{-1})` (`ι` is
an automorphism of `X` over `Q` fixing `e`, and `S_m = -S_m`), and the **semilinear** Galois action
`σ_X f := f^*`, `f^*(x) = f̄(x^{-1})` (coefficients conjugated). `V_m^{σ} = H^0(X, mB)` (the `Q`-rational sections:
`c_{-s} = c̄_s`), and `V_m = V_m^σ ⊕ iV_m^σ`; `V_m = V_m^+ ⊕ V_m^-` (`ι`-even/odd). On `V_m^σ` the two "isotypic"
refinements are `V_m^σ ∩ V_m^+` (coefficients `c_s = c_{-s} ∈ Q`) and `V_m^σ ∩ V_m^-` (`c_s = -c_{-s} ∈ iQ`). Put
`t_±(m) := max{ord_e f : 0 ≠ f ∈ V_m^±}` and `t(m) := t_max(m) = ⌊10m/3⌋` (width Thm 2.5).

**Proposition 3.1 (isotypic orders; proved, exact check `theory_checks/isotypic.py`, `sections.py`).**
(a) For `f ∈ V_m`, `ord_e f = min(ord_e f^+, ord_e f^-)`, where `f^± = (f ± ι^*f)/2`; `ord_e f^+` is even and
`ord_e f^-` is odd (or `∞`). Also `ord_e(σ_X f) = ord_e f`, and the maximal order on `V_m^σ` and on each
`V_m^σ ∩ V_m^±` equals the maximal order on `V_m` resp. `V_m^±`.
(b) For **every** `m ≥ 1`: `t_ε(m) = t(m)` for the parity `ε` of `t(m)`, and `t_{-ε}(m) = t(m) - 1`. Explicitly, with
`m = 3s + r`:

| `m` | `t(m)` | even part `t_+` | odd part `t_-` |
|---|---|---|---|
| `3s` (`s ≥ 1`) | `10s` | `10s`: `V^s` | `10s - 1`: `V^{s-1}V_{123}V_{145}V_{245}` |
| `3s+1` | `10s+3` | `10s+2`: `V^s D_{12}` | `10s+3`: `V^s V_{123}` |
| `3s+2` | `10s+6` | `10s+6`: `V^s V_{123}V_{145}` | `10s+5`: `V^s V_{123}D_{45}` |

Here `V = Π_{a<b}(x_a-x_b)/(x_1⋯x_5)^2 ∈ V_3^+` (order 10), `V_{abc} = (x_a-x_b)(x_b-x_c)(x_a-x_c)/(x_ax_bx_c) ∈ V_1^-`
(order 3), `D_{ab} = x_a/x_b + x_b/x_a - 2 ∈ V_1^+` (order 2).
(c) Exact: `t_+(m), t_-(m)` computed by certified ranks for `m = 1..5` (`isotypic.out`): `(2,3), (6,5), (10,9),
(12,13), (16,15)`, in agreement with (b).
*Proof.* (a) Write `F(y) = f(exp y) = Σ_j F_j(y)` (homogeneous parts, `y ∈ N_K`, `y_5 = 0`). Then
`(ι^*f)(exp y) = F(-y) = Σ (-1)^jF_j`, so `f^+ ↔ Σ_{j even}F_j` and `f^- ↔ Σ_{j odd}F_j`; the claims on orders
follow. `(σ_Xf)(exp y) = \overline{F(-ȳ)}`, whose homogeneous parts are `(-1)^j \overline{F_j(ȳ)}`, nonzero iff
`F_j ≠ 0`. The subspaces `{ord_e ≥ t}` are `σ_X`-stable (and `ι`-stable), so by Galois descent for `K`-vector
spaces with a semilinear involution they are spanned by their `σ_X`-fixed vectors; hence the maximal orders on
`V_m^σ`, `V_m^σ ∩ V_m^±` equal those on `V_m`, `V_m^±` (`ι` and `σ_X` commute, so descent applies inside `V_m^±`).
(b) Upper bounds: `t_± ≤ t(m)` by width Thm 2.5; by (a) the part whose parity differs from that of `t(m)` cannot
reach `t(m)`, so it is `≤ t(m) - 1`. Lower bounds: the table. Products: `supp(fg) ⊆ supp f + supp g`, and
`S_{m_1} + S_{m_2} ⊆ S_{m_1+m_2}` (`R_5` is convex, `m_1R_5 + m_2R_5 = (m_1+m_2)R_5`); `ord_e(fg) = ord_e f + ord_e g`
(initial forms multiply in a domain); parities multiply. The building blocks (support, parity, order) are checked
exactly in `sections.out`: `V ∈ S_3`, even, order 10; `V_{123} ∈ S_1`, odd, order 3; `D_{12} ∈ S_1`, even, order 2;
`V_{123}V_{145}` (even, 6, `S_2`), `V_{123}V_{145}V_{245}` (odd, 9, `S_3`), `V_{123}D_{45}` (odd, 5, `S_2`). `□`

So **no isotypic part has a higher vanishing threshold**: the "Galois-twisted vanishing at `e`" is never forced
higher, because `e` is fixed by both `ι` and `σ_X` and orders are read off initial forms, which the involutions only
multiply by signs. The best each part can do is `t(m)` or `t(m) - 1`.

**Proposition 3.2 (Liouville with the Galois structure gives nothing beyond widths; proved).** Let `0 ≠ f ∈ Z[i][M]`
with Newton polytope `Q_f`, and `y ∈ T_5(Q)` with `f(y) ≠ 0`. Write `v = v_π(y)` at each split `p`.
(a) The product formula over `K`, with the ultrametric bounds at all finite places, gives exactly
`ord_e f · λ(y) ≤ ½ Σ_p width_{Q_f}(v_π(y)) log p + O_f(1)`; the right side depends on `Q_f` only through
`Q_f + (-Q_f)`.
(b) `f^*(y) = \overline{f(y)}`; hence `f`, `f^*` and `f f^*` (which is `Q`-rational, has polytope `Q_f - Q_f` and order
`2 ord_e f`) give the **same** inequality.
(c) If `f` is `Q`-rational then `f(y) ∈ Q`, so `v_π(f(y)) = v_{π̄}(f(y))`. Using this "reality" instead of the
separate bounds at `π`, `π̄` gives `v_p(f(y)) ≥ max(min_{Q_f}⟨s,v⟩, -max_{Q_f}⟨s,v⟩)`, which equals
`-½ width_{Q_f}(v)`, the bound already used in (a), because `Q_f = -Q_f`.
*Proof.* (a) Normalise `|π|_π = p^{-1}`, `|·|_∞` the usual absolute value; the product formula over `K` reads
`2 log|α| = -Σ_{w∤∞} log|α|_w` for `α ∈ K^*`. For `α = f(y)` and coefficients in `Z[i]`:
`log|α|_π ≤ max_{s∈Q_f} log|x^s|_π = -min_{s∈Q_f}⟨s,v⟩ log p`, `log|α|_{π̄} ≤ -min_{s∈Q_f}⟨s,-v⟩ log p =
max_{s∈Q_f}⟨s,v⟩ log p` (Lemma 2.1(a)), and `log|α|_w ≤ 0` elsewhere (Lemma 2.1(b)). Summing:
`2 log|f(y)| ≥ -Σ_p width_{Q_f}(v) log p`. Taylor's formula at `e` (with `|x_j| = 1`) gives
`log|f(y)| ≤ -ord_e f·λ(y) + O_f(1)`. Finally `width_Q(v) = h_Q(v) + h_Q(-v) = h_{Q+(-Q)}(v)`.
(b) `f^*(y) = f̄(y^{-1}) = f̄(ȳ) = \overline{f(y)}` since `x^{-1} = x̄` (Lemma 2.1(c)); `Q_{f^*} = -Q_f`; the Newton
polytope of a product is the Minkowski sum; `ord_e f^* = ord_e f` by Prop 3.1(a); `width_{Q-Q} = 2 width_Q`.
(c) For `Q`-rational `f`, `c_{-s} = c̄_s`, so `\overline{f(y)} = Σ c̄_s x^{-s} = Σ c_{-s}x^{-s} = f(y)`; and a rational
number has equal valuations at `π` and `π̄`. If `Q_f = -Q_f`, then `min_{Q_f}⟨s,v⟩ = -max_{Q_f}⟨s,v⟩ = -½width`. `□`

**Corollary 3.3 (Candidate 1 reduces to the split section cap; proved).** For any finite family `𝔉 ⊂ K[M]` of any
isotypic types, used through Liouville with forced divisibility and with the reality of values on `T_5(Q)`, the
conclusions are exactly those of the `Q`-rational, centrally symmetric family `{ff^* : f ∈ 𝔉}`. Hence (toric
Thm 3.3 with (τ)) for arbitrary profiles the best exponent is `δ = 2/τ_e - 1/2 = 1/10` (Cilleruelo–Córdoba), and on
`η`-balanced profiles the best is `μ ≤ 4 + O(η)` (MW), i.e. the Liouville constant `15/4`. The split form has the same
caps (`t_max` is the same, width Thm 2.5; the split balanced cap `15/4` is attained on the split hypersurface `W`,
split Thm 2). *Where the twist entered and why it is harmless:* on twisted points `|f| = |f^*|`, which is what allows
the symmetrisation in toric Thm 3.3; the symmetrisation does not raise `t_max`.
*Proof.* Prop 3.2(b),(c) and Prop 3.1; then toric Thm 3.3 and (MW) as recalled in §1.3. `□`

## 4. Candidate 2: a point times its conjugate; Thue-type use of several approximations

**4.1 The conjugate of a twisted point.** As a point of `X`, `y ∈ T_5(Q)` is its own Galois conjugate. In the
`K`-coordinates of the split model the coordinate-conjugate point is `c(y) = ι(y) = y^{-1}`, the class of the
conjugate tuple `ζ̄`. It is again an A-point with the same height and distance (`|x̄ - 1| = |x - 1|`, `h_B∘ι = h_B`
since `R_5 = -R_5`), and `v_π(ι y) = -v_π(y) = v_{π̄}(y)`: the conjugate point carries the same contacts, exchanged
between `π` and `π̄`. The group-theoretic combinations are trivial: `y·ι(y) = e`, `y·ι(y)^{-1} = y^2`.

**Proposition 4.1 (conjugate pairs are single points; proved).** Let `F ∈ Z[i][M × M]` with Newton polytope
`Q_F ⊆ d_1R_5 × d_2R_5`, and `y ∈ T_5(Q)` with `F(y, ιy) ≠ 0`. Put `g := Δ_ι^*F`, `g(x) = F(x, x^{-1})`. Then
`g ∈ K[S_{d_1+d_2}]`, `g(y) = F(y,ιy)`, `ord_e g ≥ min{j_1 + j_2 : F_{(j_1,j_2)} ≠ 0}`, and the toric Liouville
inequality (TL) for `g` at `y` implies the two-point Liouville inequality for `F` at `(y, ιy)`. Conversely every
`g ∈ K[S_{d}]` is `Δ_ι^*(g ⊗ 1)`. So the conjugate pair gives exactly the single-point class (TL), capped as in
Cor 3.3.
*Proof.* `x^s ⊗ x^{s'} ↦ x^{s - s'}` and `S_{d_1} - S_{d_2} = S_{d_1} + S_{d_2} ⊆ S_{d_1+d_2}`. Taylor: with
`x = exp y`, `F(exp y, exp(-y)) = Σ_J F_J(y, -y)`, so every homogeneous part of `g` of degree `j` comes from
`F_J`, `j_1 + j_2 = j`; hence the order bound. Two-point Liouville for `F` at `(y, ιy)`: upper
`log|F| ≤ -min_J (j_1+j_2) λ(y) + O(1)` (`λ(ιy) = λ(y)`); lower, at `π`: `-min_{(s,s')∈Q_F}(⟨s,v⟩ + ⟨s',-v⟩)
= -min_{s - s' : (s,s') ∈ Q_F}⟨s - s', v⟩ ≥ -min_{Q_g}⟨·,v⟩` because `Q_g ⊆ {s - s'}` (cancellation can only shrink
`Q_g`); the same at `π̄`. So the (TL) inequality for `g` has a larger left side and a smaller right side. `□`

**4.2 Several distinct approximations (Thue/Roth on `X^n`).** Thue's method combines several *distinct* good
approximations through one auxiliary polynomial on `X^n` vanishing to high weighted index at `(e,…,e)`. The gain in
Roth comes from targets of degree `≥ 2`; at a rational target the index is additive.

**Proposition 4.2 (index additivity at a rational target; proved).** Let `S_1,…,S_n ⊂ M` be finite, `V_i = K[S_i]`,
`R_i := max{ord_e f : 0 ≠ f ∈ V_i} = reg(S_i)` (toric Lemma 3.1). Let `0 ≠ F ∈ V_1 ⊗ ⋯ ⊗ V_n ⊂ K[M^n]` and
`F(exp y_1, …, exp y_n) = Σ_J F_J` its multihomogeneous Taylor expansion at `(e,…,e)`. Then some `J` with
`F_J ≠ 0` has `j_i ≤ R_i` for every `i`. Consequently, for all weights `w_i > 0`,
`min{Σ_i w_i j_i : F_J ≠ 0} ≤ Σ_i w_i R_i`, with equality for `F = ⊗_i f_i`, `ord_e f_i = R_i`.
*Proof.* For each `i` choose a basis `(f_{i,k})_k` of `V_i` adapted to the order filtration
`V_i^{≥t} = {ord_e ≥ t}`: for each `t`, the basis vectors of order exactly `t` map to a basis of
`V_i^{≥t}/V_i^{≥t+1}`. The initial-form map `V_i^{≥t}/V_i^{≥t+1} → Sym^t` is injective, so the initial forms
`in(f_{i,k})` of the basis vectors of a given order are linearly independent. Expand
`F = Σ_k c_k f_{1,k_1} ⊗ ⋯ ⊗ f_{n,k_n}` and put `J(k) = (ord f_{i,k_i})_i`. Choose `k^*` with `c_{k^*} ≠ 0` and
`J^* := J(k^*)` lexicographically minimal among such `k`. A term `k` contributes to the multidegree-`J^*` part only if
`J(k) ≤ J^*` componentwise; if moreover `J(k) ≠ J^*` then `J(k) <_lex J^*`, so `c_k = 0`. Hence
`F_{J^*} = Σ_{k : J(k) = J^*} c_k ⊗_i in(f_{i,k_i})`. Products of linearly independent families of polynomials in
disjoint sets of variables are linearly independent, so `F_{J^*} ≠ 0`, and `j^*_i = ord f_{i,k^*_i} ≤ R_i`. `□`

**Proposition 4.3 (the `n`-point Liouville inequality is capped by the one-point caps; proved).** Let
`y_1,…,y_n ∈ T_5(Q)` and `0 ≠ F ∈ Z[i][M^n]` with `F(y_1,…,y_n) ≠ 0`, Newton polytope `Q_F`, projections
`Q_i = pr_i(Q_F)`, `S_i = Q_i ∩ M`. The `n`-point Liouville inequality with toric forced divisibility is
`(*)   min_{F_J ≠ 0} Σ_i j_i λ(y_i) ≤ ½ Σ_p width_{Q_F}((v_π(y_i))_i) log p + O_F(1).`
(i) Its left side is `≤ Σ_i R_i λ(y_i)`, `R_i = reg(S_i)` (Prop 4.2).
(ii) If the `y_i` have pairwise disjoint sets of primes, its right side equals `Σ_i ½Σ_p width_{Q_i}(v_π(y_i)) log p
+ O(1)`, the sum of the one-point right sides for polytopes `Q_i`.
(iii) Consequently, for data with pairwise disjoint primes, simple profiles (one layer per prime) and
`η`-balanced cuts, `(*)` holds whenever every `μ(y_i) ≤ 4(1-η)`. So the `n`-point method excludes no configuration
of balanced points in the window `μ ∈ (4 - c, 4(1-η)]`: its balanced cap is `15/4`, as for one point. The same
holds with `ιy_i` among the points (Prop 4.1) and for the split form.
*Proof.* (*) is derived exactly as Prop 3.2(a), with Taylor's formula on `X^n` giving
`log|F(y)| ≤ log Σ_{J ∈ Min} e^{-Σ j_iλ_i} + O_F(1)` over the finite set `Min` of minimal multidegrees of the Taylor
support (Dickson), and with `v_π(y_1,…,y_n) = (v_π(y_i))_i ∈ N^n`. (i) is Prop 4.2 with `F ∈ ⊗ K[S_i]`. (ii) If only
`y_i` has contact at `p`, then `width_{Q_F}(0,…,v,…,0) = width_{Q_i}(v)`. (iii) For a simple profile,
`width_{Q_i}(v_π) = width_{Q_i}(1_A)` for the one cut `A` at `p`, so the `i`-th right side is
`½ Σ_{cuts} width_{Q_i}(1_A) log n_A^{(i)} ≥ ½(1-η)w_i Σ_{cuts} width_{Q_i}(1_A) ≥ 4(1-η) w_i R_i` by (MW)
(`R_i = reg(S_i) ≤ 2^{-3}Σ_{cuts} width_{Q_i}(1_A)`), while the left side is `≤ Σ_i R_i μ(y_i) w_i`. `□`

*Remarks.* (1) Without disjointness the right side of (*) can be smaller than the sum of one-point sides (e.g.
`F(x,x') = g(x/x')` sees only `v_1 - v_2`); but then (*) is (TL) for `g` at the **point** `y_1y_2^{-1} ∈ T_5(Q)`, with
`λ(y_1y_2^{-1}) ≥ min(λ_1,λ_2) - log 2`: a one-point inequality for another twisted point, capped by Cor 3.3. A
general `F` interpolates between these two cases; I do not claim a cap for every overlapping configuration, and none
is needed: a proof of `A_5(δ)` must handle A-points with disjoint prime sets, and on those (iii) applies.
(2) This is the standard reason Thue–Siegel–Roth cannot beat Liouville at a rational point (on `P^1` the same
argument gives index `≤ Σ_i 1`), now with the toric forced divisibility included. The twist enters nowhere: (i)–(iii)
hold verbatim for split points (`v_π = v_{π̄}`), and the conjugate point `ιy` is Prop 4.1.
(3) Products of A-points in the group `T_5(Q)` are points; any inequality for them is a one-point inequality. The
round-23 referee's correction stands: products give no lower bound on `α` (they are upper-bound constructions).

**Corollary 4.4 (Candidate 2 reduces; proved).** Conjugate pairs reduce to single sections of `(d_1+d_2)B`
(Prop 4.1); several distinct approximations at the rational target `e` have additive index (Prop 4.2) and, on
balanced data with disjoint primes, the `n`-point Liouville inequality is implied by the one-point caps (Prop 4.3).
No gain; cap `15/4` (balanced). For arbitrary profiles the same argument with (τ) in place of (MW) gives the cap
`δ = 1/10` for each fixed `F`: put the whole height of `y_i` on a cut maximising `width_{Q_i}(1_A)`; then the `i`-th right
side is `½ max_A width_{Q_i}(1_A)·W'_i ≥ R_i W'_i/τ_e` (toric Thm 3.3, `τ_e = 10/3`). The witnessing profile depends on
`F` (as in toric Thm 3.3) and is not claimed uniform in `F`.

## 5. Candidate 3: parity of the contact profiles of points

The curve parity lemma (split §1.4) says: a curve over `Q` meets `E_A` and `E_{A^c}` equally. Its point form is
Lemma 2.1(a): `v_{π̄}(y) = -v_π(y)`, so over `K`, `h_{E_A}(y) = h_{E_{A^c}}(y)` for every `A`.

**Proposition 5.1 (point parity is the Gaussian profile itself; proved).** For `y = [ζ] ∈ T_5(Q)`, the complete
family of finite local data `(v_w(y))_w` is determined by the Gaussian exponent profile `(a_j(p))`: `v_π(y) =
(a_j(p))_j mod Z·1`, `v_{π̄}(y) = -v_π(y)`, all other `v_w = 0`. Conversely the profile is `v_π(y)` up to the
normalisation of `π`. Hence point parity is **not** an additional constraint on a point beyond its Gaussian profile:
any object built on a Gaussian profile (in particular every fake circle of the record) satisfies it automatically.
*Proof.* `gauss.md` Lemma 1.1 computes `v_π(ζ_j/ζ_1) = a_j(p) - a_1(p)`; Lemma 2.1 gives the rest. A fake (angle
model) keeps "the genuine primes and the genuine exponent profile" (`overview/overview.md` §5.1) and replaces only
the angles `arg π_p` by reals; its contact data at `π`, `π̄` are therefore `±v_π`. `□`

**Corollary 5.2 (parity cannot prove `(L_5^∘)` without global input; record-conditional).** Any inequality whose
hypotheses are the primes, the Gaussian profile (equivalently the local data at all finite places, with parity),
and the archimedean distance, and which holds for every actual point, is a "local input" in the sense of the barrier
map. By (Fake), balanced fakes at diameter `≍ R''^{7/15}` satisfy every monomial gap. So no inequality whose proof
uses only these data **and monomial gaps** can give `μ ≤ 4 - c`. *Status:* this uses the record's fake statement
(killer-checked, not refereed); for inequalities proved via sections the cap is unconditional (Cor 3.3).

**Proposition 5.3 (parity weakens, never strengthens, toric Liouville; proved).** For `f ∈ Z[i][M]` write
`r^{sp}_f(v) := -min_{s∈Q_f}⟨s, v⟩` (the split local right side at a prime with valuation vector `v`, weight 1) and
`r^{tw}_f(v) := ½ width_{Q_f}(v)` (the twisted local right side at `p`, both places). Then
`r^{tw}_f(v) = ½(r^{sp}_f(v) + r^{sp}_f(-v))`. Consequently: (a) if the split data `(λ, (v_p))` and `(λ, (-v_p))`
both satisfy (TL) for `f`, then the twisted data `(λ, (v_π))` satisfy it; (b) for the pair inequality (`f = x_a/x_b -
1`), the split bound is `λ ≤ min(X_{ab}, Y_{ab})` with `X_{ab} = Σ_{A ∋ a, A ∌ b} log n_A^{sp}`,
`Y_{ab} = Σ_{A ∋ b, A ∌ a} log n_A^{sp}`, and the twisted bound is `λ ≤ ½ d_ab`, the mean of the two split sums for the
symmetrised ray heights. So the parity of points can only make the twisted Liouville inequalities weaker than or
equal to the split ones at the same cut heights.
*Proof.* `width_Q(v) = max_Q⟨s,v⟩ - min_Q⟨s,v⟩ = r^{sp}(-v) + r^{sp}(v)`. (a) Sum over primes; (TL) is linear in the
local right sides. (b) For `f = x_a x_b^{-1} - 1`, `Q_f = {0, e_a - e_b}`, `r^{sp}(v) = max(0, v_b - v_a)`; summing over
the layers of a split point gives `ord·λ ≤ Σ_{A∋b, A∌a} log n_A` and, from `x_b/x_a - 1`, the other sum; split §2
(pair bound) is the same statement. For twisted points `½ width = ½|v_a - v_b|`, summing to `½ d_ab`
(Prop 2.2(iv)). `□`

*Reading.* In the split form a point may put all of a cut's height on one side (`n_A` large, `n_{A^c} = 1`); the
split hypersurface `W` does exactly this (weight `1` on one side of each cut) and attains `15/4`. In the twisted form
the height of every cut is shared equally between `E_A` and `E_{A^c}`. Prop 5.3 shows that this sharing is a
*symmetrisation* of the split local data, and symmetrisation can only weaken Liouville. So the anisotropy that
kills the `W`-curves (a statement about curves) does not translate into a pointwise inequality.

**What point parity does give (for the record; proved, but not an inequality on A-points).** Combined with the
curve parity lemma: if a curve `C` over `Q` carries infinitely many points of `T_5(Q)` with given cut proportions
`(ℓ_A)`, then `E_A·C = E_{A^c}·C` and `ℓ_A ∝ 2E_A·C`. For balanced points this forces `B·C ∈ 30Z_{>0}`, hence, at
`α_C ≤ 15/4`, total branch contact `≥ 8` at `e` (round 23). This constrains *where Liouville-extremal points could
lie if they lie on curves of bounded degree*; turning it into a statement about points off a fixed closed set would
need that A-points concentrate on such curves, which is McKinnon/Vojta-type and open (and bounded-degree covering
families are already excluded, `gauss.md` Thm 4.3(ii)). So the parity enters the problem only through curves.

**Corollary 5.4 (Candidate 3 reduces).** Point parity is the Gaussian profile (Prop 5.1); with sections it gives at
most the split section cap (Prop 5.3 and Cor 3.3); with local data and monomial gaps it is fake-satisfied (Cor 5.2).

## 6. Candidate 4: the Hilbert-90 parametrisation `x = P/P̄` and its divisibility

Prop 2.2 writes every twisted point as `x_j = P_j/P̄_j` with primitive `P_j ∈ Z[i]`. Candidate "extra divisibility"
comes in three forms: divisibility of the brackets `[P_a,P_b]`; polynomial identities among the `P_j` (Plücker /
Ptolemy, norm relations); and auxiliary polynomials in `(P, P̄)` used with forced divisibility.

**Proposition 6.1 (the bracket divisibility is the pair bound; proved, exact `h90.out`).** For all `a ≠ b`:
`v_p([P_a,P_b]) ≥ c_ab(p)` (Prop 2.2(iv)), and with `t_ab` as in `gauss.md` §1.2,
`t_ab^2 = [P_a,P_b]^2 / (Π_p p^{2c_ab(p)} · 2^{e_a+e_b})`.
So `|t_ab| ≥ 1` (the record's pair inequality) is exactly the statement `Π_p p^{c_ab(p)} | [P_a,P_b]` plus
`[P_a,P_b] ≠ 0`, when `e_a = e_b`; when `e_a ≠ e_b`, `t_ab^2 ∈ Z/2` (the record's "`√(odd/2)`" across unit classes,
`kill_gauss.md` §7 item 1). Hilbert 90 adds no divisibility beyond the pair gcd.
*Proof.* Prop 2.2(ii),(iv): `|x_a/x_b - 1| = 2|[P_a,P_b]|/(|P_a||P_b|) = 2|t_ab| e^{-d_ab/2}`, `|P_j|^2 = 2^{e_j}
Π_p p^{|ε_j(p)|}`, and `c_ab - (|ε_a|+|ε_b|)/2 = -|ε_a - ε_b|/2`. Exact check on 300 random primitive 5-tuples
(3000 pairs; 1744 with `e_a = e_b`, all with `t_ab ∈ Z`; 1256 with `e_a ≠ e_b`, all with `2t_ab^2 ∈ Z`). `□`

**Proposition 6.2 (polynomials in `(P,P̄)` are sections; proved).** Let `G ∈ Z[i][U_1,…,U_5,V_1,…,V_5]` have degree
`deg_j` in `(U_j,V_j)` jointly, for each `j`, and evaluate at `U_j = P_j`, `V_j = P̄_j`. Then
`G(P,P̄) = Π_j P̄_j^{deg_j} · g(x)`, `g(x) := Σ_{(c,d)} G_{c,d} x^c` (`x_j = P_j/P̄_j`, `0 ≤ c_j ≤ deg_j`). The
Liouville inequality for `G` with forced divisibility (`G(P,P̄) ∈ Z[i]`, `v_π(G) ≥ min_{(c,d)} Σ_j(c_j v_π(P_j) +
d_j v_π(P̄_j))`, same at `π̄`, and the archimedean bound from the order `t` of `G` along `{x_1 = ⋯ = x_5}`) is
**identical** to (TL) for `g` at `y`: `t·λ(y) ≤ ½Σ_p width_{Q_g}(v_π(y)) log p + O(1)`. If `G` is homogeneous in the
common factor (`G(γP, γ̄P̄) = γ^aγ̄^bG`), then `g` is homogeneous and `g/x_1^{Σc}` lies in `K[M]`; otherwise fix the
representative with `P_1 = 1` (allowed, Prop 2.2(i)), which makes `x` the coordinates of `y` and `g` a section.
*Proof.* `U^cV^d/V^{c+d} = (U/V)^c`. Valuations: with `α_j = v_π(P_j)`, `β_j = v_π(P̄_j) = v_{π̄}(P_j)`,
`Σ_j(c_jα_j + d_jβ_j) = Σ_j deg_jβ_j + ⟨c, ε(p)⟩` (`d_j = deg_j - c_j`, `ε = α - β`). So `v_π(G) ≥ Σ deg_jβ_j +
min_c⟨c,ε⟩` and `v_{π̄}(G) ≥ Σ deg_jα_j - max_c⟨c,ε⟩`. Adding and using `log|Π P̄_j^{deg_j}| = ½Σ_p Σ_j deg_j
(α_j+β_j) log p + O(1)` (`N(P_j) = 2^{e_j}Π_p p^{α_j+β_j}`, Prop 2.2(iii)), the Liouville inequality for `G` becomes
`t λ ≤ ½Σ_p width_{c}(ε(p)) log p + O(1)` with `ε(p) = v_π(y)` (Prop 2.2(iii)). The order of `G` along
`{x_j equal}` equals `ord_e g` (the factor `Π P̄_j^{deg_j}` is a unit near that locus). `□`

*Remark (why no half-integer gain).* Writing `U^cV^d = x^{(c-d)/2}N^{(c+d)/2}` would suggest Laurent polynomials in
`x^{1/2}`, for which the mean-width bound (MW) is false (e.g. `S = {0, ½} ⊂ ½Z`: `reg = 1 > MW = ½`). Prop 6.2 avoids
this: all exponents of a given `G` lie in one coset `deg/2 + Z^5`, and dividing by `Π P̄_j^{deg_j}` instead of
`Π N(P_j)^{deg_j/2}` lands in integer exponents `c`. So (MW) and (τ) apply to `g` as they stand.

**Proposition 6.3 (reality of brackets gives nothing; proved).** If `G` is real on all points by a polynomial identity
(`\overline{G(P,P̄)} = ±G(P,P̄)` identically, e.g. `G = [P_a,P_b]` or any polynomial in brackets and norms), then
`G(P,P̄) ∈ Z ∪ iZ` and `v_π(G) = v_{π̄}(G)`, but the forced bounds at `π` and `π̄` already coincide, so the "reality"
bound `max(f_π, f_π̄)` equals the bound used in Prop 6.2.
*Proof.* The identity means `G_{d,c} = ±\overline{G_{c,d}}`, i.e. the support `{c}` of `g` satisfies
`c ↦ deg - c`; then `Q_g` is symmetric about `deg/2`, and `min_c⟨c,ε⟩ - ⟨deg/2,ε⟩ = -(max_c⟨c,ε⟩ - ⟨deg/2,ε⟩)`,
which makes the two forced bounds equal (compare Prop 3.2(c)). `□`

**Identities (Plücker = Ptolemy, norms).** `[P_a,P_b][P_c,P_d] - [P_a,P_c][P_b,P_d] + [P_a,P_d][P_b,P_c] = 0` and
`N(P_j) = 2^{e_j}Π_p p^{|ε_j(p)|}` are the record's Gaussian Ptolemy and norm relations (`growth/ptolemy.md` Thm A(3);
`gauss.md` §4.4). As auxiliary polynomials they are covered by Prop 6.2; as three-term relations fed to abc they fall
under barrier item 6 (`round8/weakest.md` Prop 4.2: vacuous on squarefree profiles, which contain every recorded
fake). All of them are bihomogeneous in `(ζ, ζ̄)` after `x_j = P_j/P̄_j`, so Lemma K1 applies: they are unchanged by the
gcd reduction and cannot see the move from `1/2` to `7/15`.

**Corollary 6.4 (Candidate 4 reduces; proved).** Under Hilbert 90 every divisibility statement is either the pair
gcd (Prop 6.1), or (TL) for a section `g` (Props 6.2, 6.3), or a three-term identity covered by the abc barrier.
Caps: Cor 3.3. The parametrisation is form-blind: the same Props hold for the split form with `P_j = (a_j, b_j)`,
`[·,·] = det`, `P̄_j` replaced by the second coordinate, and all primes allowed.

## 7. The common reason: Galois averaging is lossless (a Lemma K1-type invariance)

**Theorem 7.1 (Galois-blind toric Liouville; proved, from §§3–6 and the record's caps).** Let `𝓛` be the class of
arguments of the following form: choose finitely many auxiliary objects, each of one of the types
(a) `f ∈ Z[i][M]` of any `σ_X`/`ι`-isotypic type (§3);
(b) `F ∈ Z[i][M^n]`, `n ≥ 1`, evaluated at a tuple whose entries are A-points, their conjugates `ιy`, or group words
in them (§4);
(c) `G ∈ Z[i][P, P̄]` in the Hilbert-90 variables (§6);
and derive, for points where the object does not vanish, the inequality "archimedean Taylor upper bound `≥` product
formula lower bound", with the lower bound taken from integrality, toric forced divisibility at every finite place of
`K`, point parity (Lemma 2.1) and reality of values (Props 3.2(c), 6.3). Then:
1. every inequality produced is, up to `O(1)`, either (TL) for one section `g ∈ K[S_m]` at one point of `T_5(Q)`, or
   the `n`-point inequality (*) of Prop 4.3;
2. on `η`-balanced simple profiles (with pairwise disjoint primes in case (b)), every inequality produced holds for
   all data with `μ ≤ 4(1-η)`; so no argument in `𝓛` proves `(L_5^∘)` or `A_5(δ)` for any `δ < 1/30`; for arbitrary
   profiles no argument in `𝓛` gives `δ < 1/10`, the witnessing profile depending on the family;
3. the same statements, with the same numbers, hold for the split form (where the balanced value `15/4` is moreover
   attained on the hypersurface `W`, split Thm 2).
*Proof.* Item 1: (a) Prop 3.2; (b) Props 4.1 and 4.3 (group words are points of `T_5(Q)`); (c) Props 6.2, 6.3;
parity is Lemma 2.1, already used in (TL) (Prop 3.2(a)). Item 2: Cor 3.3 for one point (via (MW), (τ) and toric
Thm 3.3), Prop 4.3(iii) and Cor 4.4 for several. Item 3: the proofs never used `v_{π̄} = -v_π` except to evaluate
(TL); for split points the local right side is `r^{sp}` (Prop 5.3), and (MW), (τ), Prop 4.2 are statements about
Laurent polynomials, valid for both forms (width Thm 2.5 is stated for both). `□`

**The mechanism in one line.** The Galois involution `σ_X` fixes `e`, permutes the places above each split `p`
(`π ↔ π̄`) and acts on sections preserving orders at `e`. Toric Liouville is a sum of local terms over places, so
averaging an argument over its Galois orbit (`f ↦ f f^*`, `(y, ιy) ↦ Δ_ι^*`, `G ↦ g`) loses nothing, and the averaged
argument is a `Q`-rational, centrally symmetric toric Liouville inequality. Its strength is governed by `t_max(m)`
and the mean width, which are invariants of the lattice polytope `R_5`, i.e. of the geometric (form-independent)
variety. So **the anisotropy changes which points exist, not which section inequalities hold.**

**Relation to Lemma K1.** Lemma K1 says the gcd reduction is invisible to bihomogeneous forms; Theorem 7.1 says
the Galois twist is invisible to them too. Both have the same proof shape: the method's inequality is invariant
under an operation (rescaling by `g`; Galois averaging) that changes the arithmetic of the points. Every object in
`𝓛` is bihomogeneous in `(ζ, ζ̄)` after `x_j = ζ_j/ζ_1`, `x̄_j = ζ̄_j/ζ̄_1`, so Lemma K1 applies to each one:
nothing in `𝓛` sees either the reduction `1/2 → 7/15` or the twist.

**What is not in `𝓛` and is anisotropy-sensitive (exhaustive for the record, to my knowledge).**
* *The curve parity lemma* (split §1.4): a statement about curves over `Q`, not about points (§5, last paragraph).
* *The set of primes* (`p ≡ 1 (4)` only) and *the Picard rank* `15` versus `26` (Lemma 2.3): these change constants
  and powers of `log B` in counts, not exponents (**heuristic** as an input to `α`).
* *Definiteness at `∞`* (`T_5(R)` compact): invisible to toric invariants at `∞` (toric Prop 6.1), and near `e` both
  forms are totally real slices of the same complex picture (§2), so orders and Taylor bounds coincide.
A gain must therefore come from an input outside `𝓛` that is also not "local" in the barrier sense: an input that
uses the *global* position of the A-points (which curves they lie on, how many there are), combined with the parity
of those curves. The only such statements available are conditional (Vojta `V_E(5)`, McKinnon's curve conjecture)
or excluded families (bounded-degree covering curves, `gauss.md` Thm 4.3(ii)).

## 8. What survives; constraints; open points

**Survives (proved here, unrefereed).** Lemma 2.1, Prop 2.2 (Hilbert-90 binary-form model, with the exact `t_ab`
formula), Lemma 2.3 (ranks), Props 3.1–3.2 (isotypic orders for all `m`; reality is lossless), Props 4.1–4.3
(conjugate pairs; index additivity at a rational target; `n`-point cap on disjoint-prime balanced data), Props 5.1,
5.3 (point parity is the profile; it symmetrises and weakens Liouville), Props 6.1–6.3 (Hilbert-90 polynomials are
sections), Theorem 7.1. Exact checks: `isotypic.out` (`m ≤ 5`; `m = 6` stopped for time), `sections.out`, `h90.out`,
`index_add.out`. **Record-conditional:** Cor 5.2 (uses the killer-checked fakes).

**Constraint audit.**
* *No uniform exceptional set from a fixed-`S` Subspace theorem:* none is used; the only exceptional sets are zero sets
  of auxiliary sections (as in toric Thm 3.3) and the hypothetical `Z` of `A_5`.
* *No Gaussian gcd or denominator cost discarded:* (TL) is derived from the full product formula over `K`, with every
  finite place (Prop 3.2(a)); Hilbert-90 polynomials carry the explicit factor `Π P̄_j^{deg_j}` and the `2`-adic
  `2^{e_j}` (Props 2.2, 6.1, 6.2); tuples are primitive and `R''` is exact (`h90.py` divides by the Gaussian gcd).
* *No drifting common phase replaced by zero:* the phase is quotiented out in `T_5`; in Prop 6.2 a non-homogeneous `G`
  is evaluated at the representative `P_1 = 1`, which is a choice of coordinates on `T_5`, not a phase set to zero.
* *No template-dependent threshold treated as uniform:* the caps are statements about inequalities (they hold for every
  family); the arbitrary-profile witnessing profile depends on the family and is said so (Cor 4.4, Thm 7.1).

**Open points.**
1. A cap for the `n`-point inequality (*) on configurations with **shared** primes (Prop 4.3 Remark (1)). Not
   needed for the negative conclusion, but it would close `𝓛` completely. *Assessment:* the interpolation between
   `g(x/x')` (a point inequality for `y_1y_2^{-1}`) and tensor products suggests the same cap.
2. An anisotropy-sensitive **global** input. The only candidate I can name is curve-based: show that balanced A-points
   in the window concentrate (off a fixed `Z`) on curves over `Q`, then use parity (`B·C ∈ 30Z`, branch contact `≥ 8`).
   The concentration step is McKinnon/Vojta-type and open; with it, `gauss.md` Thm 4.3(ii) would already finish for
   bounded degree, so parity would only sharpen degrees.
3. The twisted search ansatz of round 23 (two conjugate real-quadratic branches of order 4, `B·C = 30`) is the
   natural place where the anisotropy would show up in a construction; it is a computation, not part of this theory
   file.
