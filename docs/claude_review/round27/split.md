# Round 27, lane split: split Question B with integer-specific tools

Status: complete draft, unrefereed (round 27, lane split). Scripts in `split_checks/`, notes in `split_notes.md`.

## 0. Summary and verdict

**Verdict: `obstruction` (new, partly proved, partly conjectural with exact evidence), plus a new piece of the
exceptional set. No growth improvement, no uniform bound, no proof of split Question B or of Erdős #887.**
`M ≤ (2/3+ε) log R/loglog R` stands.

1. **New exceptional hypersurfaces (Theorem 3.1, proved, exact checks).** The second 2026 five-divisor construction
   (mnn31, Lean, `C = 7/2`: `n = u(u+2)(u+3)(u+5)`, `(u+2)(u+3) = 3(w+1)(w+2)`) is a rational curve with two conjugate
   real-quadratic branches of order 1, `B·C = 7`, hence **`α = 7/2 < 15/4`** (unbalanced: 7 of 15 cuts). It belongs to a
   3-parameter family of conics whose union is dense in the **quadric cone**
   `(u_2−u_1)(u_3−u_1) = (u_4−u_1)(u_5−u_1)`. So every admissible `Z` in split `A_5(δ)` (any `δ < 1/14`) contains these
   15 cones besides round 23's `W`. The cone is not defined over `Q` in the twisted `T_5`.
2. **Split profiles are 11-dimensional, and sections are profile-dependent (§1.3, Lemma 3.0).** The record transferred
   the balanced cap `15/4` to the split form only at the symmetric profile `x ≡ 1/2`. **Theorem 3.2 (proved):** for every
   auxiliary form with centrally symmetric support (all binomial and Vandermonde products, all known equality cases)
   the forced-divisibility cost is *independent* of the near-tight profile, so the cap `15/4` holds on the whole
   near-tight polytope `Π_nt` (exact identity + the refereed Theorem B). **Prop 3.3 (proved):** at the `W`-vertex (a
   "code" profile `1[|B ∩ R| even]`, `|R| = 3`) the cap holds for every form not divisible by the defining polynomial
   `g_W` (via the `W`-curves). **Prop 3.7 (proved)** turns the cap into a curve question: a *dominating* family of
   rational curves through `e` with profile `N x` and multiplicity `4N` caps every form at `x`; such families have
   expected dimension `3N + 1 ≥ 4` at every near-tight profile. At the `W`-vertex such a family exists
   (**numerical**: splitting the two cusps of a `W`-curve gives a 4-dimensional family of curves with four smooth branches
   at `e`, which hits 8/8 random targets while the `W`-family hits 0/8), so the cap holds there for all forms including
   `g_W`, and by `S_5`, inversion and linearity **on a 10-dimensional polytope `Q_code ∋ 1/2`** (Cor 3.8). Picture: the
   Liouville-extremal `W`-curves are the codimension-2 *cuspidal* degenerations of a dominant family of slope-`15/4`
   curves whose general member has four smooth branches at `e` and carries no A-points.
   **Conjecture WTB** (cap at every near-tight split profile): 28 712 random sets, a local search, one-sided polytopes
   and a `k = 4` control, all with maximum ratio exactly 1.
3. **Bounded pair cofactors (Prop 3.4, Lemma 3.5, Theorem 3.6; proved).** Five divisors with `|t_ab| ≤ K` are
   `S_K`-integral points of the affine variety `LM_5 ∖ ⋃Δ̄_ab` (fixed `S`, so Siegel/Corvaja–Zannier-type arguments
   are legitimate). But poles along the root arrangement, paid for by integrality, multiply only what plain sections
   already give: the Vandermonde has ratio exactly `15/4` on all of `Π_nt` (mediant lemma). So integrality is worthless
   wherever WTB holds, and decisive (no distance condition at all) wherever it fails. Each triple satisfies the exact
   identity `x^τ × y^τ = (t_bc, t_ca, t_ab)`, so bounded-residue triples are Möbius (`SL_2(Z)`/Farey) images of fixed
   rationals; this explains `A_3` false and gives no constraint at `k = 5`.
4. **Phase 0 (#887 proper; §3.6).** The forum's reduction (`z_j² − z_1² = 2k_j s + k_j²`, `k_j ≤ C² + 1`) is a family of
   complete intersections of quadrics: K3 for four divisors, general type from five on. The two known families are
   Pell conics on these surfaces (`k = (0,224,225,324,900)` and `(0,3,4,9,12)`, exact). Under Lang + Siegel, absolute
   #887 reduces to the absence of such curves for `K+1` divisors. Phase-0 only; not used for the torus problem.
5. **Transfer to `T_5` (§5).** Nothing ports as a proof. The twisted problem sits exactly at the symmetric profile,
   where the section cap is Theorem B; the split-only freedoms are profile asymmetry (worthless to symmetric forms,
   Theorem 3.2) and weight-1 curves (`W`, the cones), which the twist removes (parity lemma; Galois moves the cone).

**Status labels.** Theorem 3.1, Lemma 3.0, the cost identity, Theorem 3.2, Prop 3.3(a,c), Prop 3.4, Lemma 3.5,
Theorem 3.6, Prop 3.7: **proved** (unrefereed). Prop 3.3(d,e): **derived**. Existence of the dominant family at the
`W`-vertex, hence Cor 3.8: **numerical** (Newton continuation from an exact `W`-curve, not interval-certified).
Conjecture WTB outside `Q_code` and symmetric supports: **open** (evidence). §3.6: **derived**
(conditional reading only). §3.7: **heuristic**.

## 1. Setting and normalisations (verified)

Labels: **proved** (proof here), **exact** (exact computation in `split_checks/`), **record** (refereed earlier, cited),
**evidence**, **heuristic**, **snippet** (literature seen through search results or a raw file only).

**1.1 Cox data.** A primitive 5-tuple of distinct positive integers `u = (u_1..u_5)` determines, for every proper
nonempty `B ⊂ [5]`, the integer `n_B = Π_{layers of type B} p` (round 23 split §1.4). Then `u_a = Π_{B∋a} n_B`,
`g_ab := gcd(u_a,u_b) = Π_{B⊇{a,b}} n_B`, `lcm(u) = Π_B n_B`, and `H(u) = (max u/min u)·lcm(u)` (record). Put
`α_ab = Π_{B∋a, b∉B} n_B` (8 factors), so `u_a = g_ab α_ab`, `u_b = g_ab α_ba`, and the **pair cofactor** is
`t_ab = α_ab − α_ba = (u_a − u_b)/g_ab`.

**1.2 The coordinator's split twin (verified).** Index the 15 cuts by nonempty `A ⊆ [4]` and put `u_A = n_A`,
`v_A = n_{[5]∖A}`. Then `u_a/u_5 = U_a/V_a` with `U_a = Π_{A∋a} u_A`, `V_a = Π_{A∋a} v_A` (8 factors each), and
`U_a − V_a = t_{a5}` (proved: the sets `B ∋ a, 5 ∉ B` are the `A ∋ a`, and `B ∋ 5, a ∉ B` are the `[5]∖A`, `A ∋ a`). `gcd(u_A, v_A) = 1`
because the layers of one prime form a chain. For `a,b ∈ [4]`, `t_ab = (t_{a5}V_b − t_{b5}V_a)/Π_{A∋a,b} u_A v_A`
(proved by expanding `U_aV_b − U_bV_a`), so near-tightness of the four pairs `(a,5)` controls all ten.

**1.3 Profiles.** Write `n_B = X^{x_B + o(1)}`. *Cut-balanced* means `x_B + x_{B^c} = 1` for all 15 cuts (`log H ~ 15 log X`).
Put `X_a = Σ_{B∋a} x_B` (so `u_a = X^{X_a}`) and `P_ab = Σ_{B∋a,b∉B} x_B` (so `α_ab = X^{P_ab}`).

*Lemma 1.1 (proved).* On cut-balanced profiles `P_ab + P_ba = 8` and `P_ab − P_ba = X_a − X_b`.
*Proof.* `P_ba = Σ_{B∋b,a∉B} x_B = Σ_{C∋a,b∉C} x_{C^c} = Σ_{C∋a,b∉C}(1−x_C) = 8 − P_ab`; and
`X_a − X_b = Σ_{B∋a,b∉B}x_B − Σ_{B∋b,a∉B}x_B`. `□`

Hence for cut-balanced points with `u_a ≈ u_b` (all `X_a` equal) every `P_ab = 4`, so `dist ≍ max_ab |t_ab| X^{−4}`
and `H = X^{15}`: an A-point at `δ` has `max|t_ab| ≤ X^{1/4 − 15δ/2 + o(1)}`, i.e. the coordinator's elementary form.
The **near-tight profile polytope**

    Π_nt = { x ∈ [0,1]^30 : x_B + x_{B^c} = 1 ∀B,  X_1 = … = X_5 }

has **dimension 11** (exact, `profiles.py`: rank 19 of the 35 defining equations; by Lemma 1.1 only 4 conditions
beyond cut balance). It contains the symmetric profile `x ≡ 1/2` (the only one the twisted problem can produce, since
`|γ_A| = |γ̄_A|`) and exactly 278 one-sided 0/1 vertices (exhaustive over `2^15` side choices, `onesided_enum.py`) in 16 `S_5`-classes, one of which is the profile `x_W` of the round-23
Pell family `W` (exact). The record treated only `x ≡ 1/2` when it transferred the section cap to the split form.

**1.4 Moduli picture (proved, standard).** `T = G_m^5/G_m ⊂ M_{0,7}` (points `0, ∞, u_1..u_5`), and
`M̄_{0,7} → LM_5` contracts the 16 divisors `D_T`, `T ⊆ [5]`, `|T| ≥ 3`, and maps `D_{ab}` onto the subtorus
closure `Δ̄_ab = {u_a = u_b}`. The 30 toric divisors are `D_{0∪B} = E_B`. The finite part of `h_{D_ab}` is
`log|t_ab|` up to the deeper residues. So **bounded pair cofactors** (`|t_ab| ≤ K`) imply that the point is
`S`-integral, `S = {∞} ∪ {p ≤ K}`, with respect to `D_int = Σ_{T⊆[5],|T|≥2} D_T` (26 components); conversely
`S`-integrality says exactly that all `t_ab` are `S`-units, i.e. that the stable reduction at every `p ∉ S` is a
Losev–Manin chain (every node separates `0` from `∞`).

## 2. Literature: Erdős #887/#886, divisors in short intervals (2026-10-29)

Access: WebSearch; raw GitHub files (fetched, treated as data, copy in `split_checks/dl/`). `erdosproblems.com`
does not resolve in this session (forum read only through search snippets). No paper full text was read.

| item | what is proved / claimed | label |
|---|---|---|
| formal-conjectures `887.lean` (raw, fetched today) | `parts.i`, `parts.ii` (`∃K ∀C`, eventually `≤ K` divisors in `(√n, √n + C n^{1/4})`): **research open**. `rosenfeld_infinite` (4 divisors, some `C`): solved. `rosenfeld_4` (4 is the largest such count): **solved in the negative**, formal proof link `mnn31/erdos-lean …/Erdos887.lean#L538` | fetched |
| formal-conjectures `886.lean` (raw) | `erdos_886` (Ruzsa, `n^{1/2-ε}` window, `O_ε(1)`): **open**; `rosenfeld_bound` (`≤ 1 + C²` divisors in `[√n, √n + C n^{1/4}]`) and `rosenfeld_infinite` (4 divisors in `16 n^{1/4}`): solved | fetched |
| m-czech (forum, July 2026), via PR #7218 | 5 divisors infinitely often, `C = 31` (Lean `C = 40`): `n = (x²−1)(x²−4)(x²−16)(x²−49)`, `z² = 56x² + 65`; balanced Pell curve at `α = 15/4` (round 23 Thm 2) | record + snippet |
| **mnn31 `Erdos887.lean` (fetched, read)** | `five_divisors_infinite`: 5 divisors infinitely often with **`C = 7/2`**: `n = u(u+2)(u+3)(u+5)`, `(u+2)(u+3) = 3(w+1)(w+2)` (`(2u+5)² − 3(2w+3)² = −2`), divisors `(u+2)(u+3), 3(w+1)(w+3), (u+2)(u+5), (u+3)(u+5), 3(w+2)(w+3)`. Also `card_window_le_one` (`C = 1`: at most 1) and `card_window_le_four` (`C = 2`: at most 4) | fetched; **re-verified exactly here** (§3.1, `czech72.py`) |
| erdosproblems forum #887 (snippet) | the problem "reduces to showing that if `Q_1(x,y),…,Q_K(x,y)` are certain polynomials of degree 1, then there exist only finitely many `(x,y)` such that they are simultaneously perfect squares"; points to Letendre (2025) | snippet |
| T. H. Chan, arXiv 1303.2069 | a perfect square has at most **five** divisors in `[√n − c n^{1/4}, √n + c n^{1/4}]`; three-divisor examples for squares come from Pell equations (`X² − 2Y² = 2`) | snippet |
| Chan 1303.2249, 2507.07094; Chan–Holmes–Liu–Villarreal 2508.02818 | three / four close factorisations, Pell structure, phase 0 only | record (round 23) |
| Cilleruelo–Jiménez-Urroz, JTNB 12 (2000) 87–92; Granville–Jiménez-Urroz, QJM 51 (2000) 343–352 | arcs of `xy = N` with `k` lattice points; lower bounds for the lcm of `k` integers in a short interval (theorem statements not seen) | snippet; this is the Liouville regime |
| Letendre, arXiv 2503.12146 | divisors in `n^{1/4−δ}` windows (`O(1/δ)`), below the critical scale | record |
| Ford (math/0607473 and the short-interval paper), Tenenbaum | distribution in `n` of integers with divisors in short/dyadic intervals; averaged over `n`, not pointwise at the critical scale | titles only (record) |

**Exactly what is proved for 4 and 5 divisors (as far as found).**
* *Phase 0 (`X = √n`).* At most `1 + C²` divisors (Erdős–Rosenfeld); `≤ 1` for `C = 1` and `≤ 4` for `C = 2` (mnn31,
  Lean). Four divisors infinitely often (Erdős–Rosenfeld, `C = 16`); five divisors infinitely often (m-czech `C = 31`;
  mnn31 `C = 7/2`). No six-divisor construction and no absolute bound `K` were found (search-negative). Perfect squares:
  at most five (Chan).
* *General phase / torus form.* Nothing beyond Liouville: split `A_5(δ < 1/30)` is open and implies #887 `parts.ii`
  phase-uniformly (round 23 Thm 1, refereed).
* *The forum reduction* is the phase-0 Erdős–Rosenfeld substitution `s_j = d_j + n/d_j`, `z_j = d_j − n/d_j`:
  `z_j² − z_1² = 2k_j s_1 + k_j²` with `k_j = s_j − s_1 ∈ [0, C² + 1]` (re-derived here, §3.6). It is a statement at
  phase 0 only and is not used for the torus problem.

## 3. Main: mechanism(s) and results

**3.0 What was attempted, and the outcome.** I looked for an integer-specific mechanism in four places, each chosen
because the split form has structure the twisted form lacks: (i) the 11-dimensional freedom of split profiles
(`x_B ≠ x_{B^c}`; the twisted form only has `x ≡ 1/2`); (ii) integrality at the root arrangement, i.e. bounded pair
cofactors (`S`-integral points of an affine variety, so Siegel/Corvaja–Zannier type arguments become available);
(iii) the exact `SL_2(Z)` structure of bounded-residue triples; (iv) phase-0 integrality (for #887 itself). None
proves split Question B. Instead:
* a new piece of the **exceptional set**: the second 2026 five-divisor construction lies on a **quadric cone**, at
  `α = 7/2 < 15/4` (Theorem 3.1);
* a new **obstruction**: freedoms (i) and (ii) are invisible to every Liouville/section argument with centrally
  symmetric support (Theorem 3.2, Lemma 3.5), for all `f` at the symmetric profile (Theorem B, record), and — via a
  curve criterion (Prop 3.7) and a numerically constructed dominant family — for all `f` on a 10-dimensional polytope
  of profiles through the symmetric one, including the `W`-vertex (Prop 3.3, Cor 3.8). In general this is a precise
  conjecture (WTB, §3.4) with exact evidence;
* structural theorems for bounded residues (Theorem 3.6) and a geometric reading of the phase-0 reduction (§3.6).

### 3.1 Theorem 3.1: the `C = 7/2` construction lies on a quadric cone at `α = 7/2` (proved; exact checks)

Let `a, b, c, d ∈ Q^*`, `k = ab/(cd)`, and let `Q_{abcd}` be the conic `u(u+a+b) = k·w(w+c+d)`. Put

    D_1 = (u+a)(u+b),  D_2 = (u+a)(u+a+b),  D_3 = (u+b)(u+a+b),  D_4 = k(w+c)(w+c+d),  D_5 = k(w+d)(w+c+d),

and `φ_{abcd} = [D_1 : … : D_5] : Q_{abcd} → T`. (mnn31's family is `(a,b,c,d,k) = (2,3,1,2,3)`, `u ↦ u`, `w ↦ w`.)

(i) *(proved, symbolic: `fam72.py`)* On `Q_{abcd}`, `(u+a)(u+b) = k(w+c)(w+d)` identically (subtract the conic
equation: `ab = kcd`). The conic passes through `(0,0)`, so `Q_{abcd} ≅ P^1_Q`.

(ii) *(proved)* With `Y_j = D_j/D_1 − 1`: `Y_2 = a/(u+b)`, `Y_3 = b/(u+a)`, `Y_4 = c/(w+d)`, `Y_5 = d/(w+c)`, and (i) gives

    Y_2 Y_3 = Y_4 Y_5,   i.e.   (D_2 − D_1)(D_3 − D_1) = (D_4 − D_1)(D_5 − D_1).

So every curve lies on the **quadric cone** `𝒬 = {(u_2−u_1)(u_3−u_1) = (u_4−u_1)(u_5−u_1)}` (vertex `e` in the
chart `Y`). *(exact, `cone_checks.py`: the identity holds exactly at six integral Pell points of mnn31's family.)*

(iii) *(proved)* Profile. The eight rational places `(u,w) ∈ {(0,0),(0,−c−d),(−a−b,0),(−a−b,−c−d),(−a,−c),(−a,−d),(−b,−c),(−b,−d)}`
carry the valuation vectors `0` (constant place `(0,0)`) and `1_B` for
`B = {4,5}, {2,3}, {2,3,4,5}, {1,2,4}, {1,2,5}, {1,3,4}, {1,3,5}`: seven distinct cuts, each once (each `D_j` is a
product of two linear forms, each vanishing at two of these places; all `D_j` have degree 2 at the two points at
infinity). So `B·C = 7`; the curve is **not** cut-balanced (8 of 15 cuts are missed).

(iv) *(proved)* At infinity `u/w → ±√k`. If `k > 0` is not a square the two branches are conjugate real-quadratic,
and `D_j − D_1 = O(u)` against `D_j ≍ u²`, so each branch has order `μ_0 = 1`. As in round 23 Thm 2(iv) (Liouville
on `P^1` at a quadratic point, `h_{−K}∘φ = 7h + O(1)`), `dist ≍ H_{−K}^{−2/7}` infinitely often and never better:
`α_C = 7/(2·1) = 7/2` (McKinnon–Roth Thm 2.16 with `r_q = 2`, `m_q = 1`). *(exact, `czech72.py`: along 40 Pell
points of mnn31's curve `log(dist^{7/2}·H)` is periodic with values `4.064, 4.267, 4.96`; `γ → 7/2`; `lcm/n → 7/8`;
`C_actual → 2√3`.)*

(v) *(proved, exact rank: `fam72.py`)* The map `(a,b,c,d,t) ↦ (D_j/D_1)_{j=2..5}` (conic parametrised by `w = tu`)
has exact Jacobian rank 3 at 4 random rational points, and `𝒬` is an irreducible 3-fold. So the union of the curves
is Zariski dense in `𝒬`.

**Consequence (proved).** For every `δ < 1/14` (in particular every `δ < 1/30`) each curve with `k > 0` non-square
carries infinitely many A-points at `δ`, with distinct positive coordinates; such `(a,b,c,d)` are Zariski dense. Hence
**every admissible `Z` in split `A_5(δ)` contains `𝒬` and its 15 `S_5`-images** (choice of the base index and of the
pairing of the other four). This is in addition to round 23's `W`; `W ⊄ 𝒬_σ` for all 15 cones (exact, `cone_checks.py`:
m-czech points lie on none). The cones are unbalanced, so they do not enter `(L_5^split)`; they do enter `A_5^split`
and the absolute #887 reduction (round 23 Thm 1 uses `A_5`).

**Twisted form (proved, `cone_checks.py` (d)).** The real structure of `T_5` acts by `u ↦ 1/ū`; it maps `𝒬` to the
cubic hypersurface `u_1u_2u_3 − u_1u_4u_5 − u_2u_3u_4 − u_2u_3u_5 + u_2u_4u_5 + u_3u_4u_5 = 0` (base index `u_1`), which does not contain
`𝒬`. So `𝒬` is not defined over `Q` in `T_5`; its twisted trace is at most the surface `𝒬 ∩ σ𝒬`. This is a second
instance (after the parity lemma) of the twisted exceptional set being forced to be smaller.

*Note on the mechanism.* Both known five-divisor families are "two factorisations of one quartic on a Pell conic":
m-czech uses the norm identity `F_1 F̄_1 = n(x)`, mnn31 uses `u(u+a+b)·(u+a)(u+b) = k²·w(w+c+d)·(w+c)(w+d)` with a
common difference of the two factor pairs. In the torus they give rational curves with two conjugate real-quadratic
branches, the only way (McKinnon–Roth Thm 2.16, `r_q = 2`) a single rational curve can carry infinitely many points at
`α = B·C/(2μ_0)`.

### 3.2 Theorem 3.2: split profile asymmetry is invisible to symmetric sections (proved)

**Setup.** For a finite `S ⊂ M = {m ∈ Z^5 : Σm_a = 0}` (the support of a Laurent polynomial `f` in `y_a = u_a/u_5`)
put `m(B) = Σ_{a∈B} m_a`, `π_B(S) = max_{m∈S}(−m(B))` (toric pole order along `E_B`), `w_B(S) = π_B + π_{B^c}`
(width), `mid_B(S) = (max_S m(B) + min_S m(B))/2`, and `reg(S)` = the regularity index = the largest order at `e`
of a nonzero Laurent polynomial supported on `S` (record: `threshold/construct.md` §2.4, Cor B1). For a profile `x`
put `cost_x(S) = Σ_B x_B π_B(S)`.

**Lemma 3.0 (the section inequality at a profile; proved, standard).** Let `f ≠ 0` be supported on `S`, with
integer coefficients. If `y` is a primitive tuple with simple profile `x` (`n_B = X^{x_B+o(1)}`, `Σx_B = 15`) and
`f(y) ≠ 0`, then `reg`-order Taylor at `e` and forced divisibility (`v_p(f(y)) ≥ −π_B` per layer of type `B`) give
`ord_e(f)·(−log dist(y,e)) ≤ cost_x(S)·log X + O_f(1)`. Hence A-points with profile `x` at exponent
`γ < 15·ord_e(f)/cost_x(S)` lie on `{f = 0}` once `X` is large. *(This is the record's "Liouville with toric forced
divisibility", written per profile.)* So the **section exponent at `x`** is `γ_sec(x) = 15·sup_S reg(S)/cost_x(S)`.

**Identity (proved).** For every `x ∈ Π_nt` and every finite `S ⊂ M`:

    cost_x(S) = (1/2) Σ_cuts w_B(S) − Σ_B x_B·mid_B(S)  =  cost_{1/2}(S) − Σ_B (x_B − 1/2)·mid_B(S).

*Proof.* `π_B = w_B/2 − mid_B`; summing `x_B w_B/2` over the two sides of a cut gives `(x_B + x_{B^c})w_B/2 = w_B/2`;
and `Σ_B mid_B = 0` since `mid_{B^c} = −mid_B`. `□`

**Theorem 3.2 (proved).** If `S` is centrally symmetric (`S = 2c − S` for some `c`) then `cost_x(S)` is the same for
all `x ∈ Π_nt`, and `4·reg(S) ≤ cost_x(S)` for every `x ∈ Π_nt`. In particular:
* (binomials) for every `m ∈ M ∖ 0` and every `x ∈ Π_nt`, `cost_x({0,m}) = (1/4)Σ_{B}|m(B)| ≥ 4`, with equality for
  roots `m = e_a − e_b`;
* (zonotopes) every product of binomials `χ^{m_i} − 1`, in particular every Vandermonde-type product, has section
  exponent `≤ 15/4` at **every** near-tight split profile.

*Proof.* If `S = 2c − S` then `mid_B(S) = c(B)`, and `Σ_B x_B c(B) = Σ_a c_a X_a = X_c Σ_a c_a = 0` on `Π_nt` (all
`X_a` equal, `c ∈ M_R`). So `cost_x(S) = cost_{1/2}(S) = (1/8)·Σ_{P⊆[5]} w_P(S) ≥ 4·reg(S)` by the refereed Theorem B
(`construct.md` §2.4 with `n = 5`, Cor B1: `reg ≤ 2^{1-5}Σ_P w_P`, and `Σ_P w_P = 2Σ_cuts w_B`). For binomials,
`(m(B))^+ = (|m(B)| + m(B))/2` and `Σ_B x_B m(B) = 0` give `cost = (1/2)Σ_cuts|m(B)|`; pairing `B ↔ B △ {a}` for one
`a` with `m_a ≠ 0` gives `Σ_{P⊆[5]}|m(P)| ≥ 16`. `□`

**Exact data (`profiles.py`).** `Π_nt` has dimension 11; random-objective LPs found 787 vertices (coordinates in
`{0, 1/3, 1/2, 2/3, 1}`), among them all 278 one-sided 0/1 vertices (16 `S_5`-classes; exhaustive check `onesided_enum.py`); the minimum binomial cost over all of them, for
`|m_a| ≤ 3`, is exactly 4. The one-sided classes include the **code vertices** `x_R = 1[|B ∩ R| even]` for
`|R| = 3` (the class of `W`: in round 23's labels `R = {2,3,4}`, i.e. the three even polynomials `F_2, F_3, F_4`)
and `|R| = 5` (even sets). *(proved: `|B∩R| + |B^c∩R| = |R|` is odd, so exactly one side of each cut is even, and
`X_a = 8` for all `a`; exhaustively, the near-tight code families are exactly the ten `F_R`, `|R| = 3`, and `F_{[5]}`.)*

**Meaning.** The record transferred the balanced section cap `15/4` to the split form at the symmetric profile only.
A split near-tight point has an 11-dimensional profile, and the section inequality is profile-dependent
(Lemma 3.0). Theorem 3.2 shows that for every auxiliary form with centrally symmetric Newton polytope — all binomial
products, all Vandermonde products, all equality cases of Theorem B found by the record (Conjecture E) — the cost
does not depend on the profile at all. The split asymmetry is a direction along which these forms are exactly flat.

### 3.3 Proposition 3.3: the cap at the `W`-vertex (reduction proved; resolution numerical)

Let `x_R` (`|R| = 3`) be the code vertex of `W`, and let `g_W` be an irreducible Laurent polynomial defining the
hypersurface `W` (round 23 Thm 2(v), exact rank). 

**Proposition 3.3.** (a) Every Laurent polynomial `f` with `g_W ∤ f` satisfies `4·ord_e(f) ≤ cost_{x_R}(f)`.
(b) Hence the section cap `γ_sec(x_R) = 15/4` holds iff `4·ord_e(g_W) ≤ cost_{x_R}(g_W)`; if this fails, a power of
`g_W` beats `15/4` at `x_R`, and every near-tight split point with profile near `x_R` lies on `W ∪ {f_i = 0}` for a fixed
finite family `f_i` — a restricted theorem.
(c) `W` has degree `> 14` in `P^4` (exact over `F_{1000003}`: no form of degree `≤ 14` vanishes on 3100 random points of
`W`, `w_equation.py`), and its tangent cone at `e` contains the cone over the irreducible octic surface

    Σ_W :  4Y_1²Y_2²Y_3²Y_4² + Y_2⁴Y_3⁴ + Y_2⁴Y_4⁴ + Y_3⁴Y_4⁴ − 2Y_2²Y_3²Y_4²(Y_2² + Y_3² + Y_4²) = 0

of branch directions `(±2√(abcd), (a+b)(a+c), (a+b)(b+c), (a+c)(b+c))` (exact, `tangent_W.py`), so `ord_e(g_W) ≥ 8`.
(d) *(derived)* Along a branch, `F_1, …, F_4` are even in `s = 1/x` and `F_0` has odd part `−e_3 x`
(`e_3 = abc − d(ab+ac+bc)`), so each branch germ is `Y(s) = ξ s² + e_3(1,1,1,1) s³ + O(s⁴)`: a cusp, and `W` has two
sheets (`±s`) over each direction of `Σ_W`; so `ord_e(g_W) ≥ 16` if no degenerate member of the family adds sheets
through `e`.
(e) *(derived; numerical support `w_member2.py`)* `ι(W) ≠ W` for the inversion `ι: u ↦ 1/u`: `ι` sends the germ with
data `(ξ, η)` to one with data `(ξ, −iη)` (substitute `s = −iσ` in `−ξs² − ηs³`), while `W`'s sheets over `ξ` carry
`±η` modulo `ξ`; equality would force `(1,1,1,1) ∥ ξ(θ)` for all `θ`, which fails. Numerically, the inversions of two
`W`-points have membership residuals `1.4·10^{−1}` and `5·10^{−2}` against `≤ 1.3·10^{−15}` for `W`-points and for the
swap `F_2 ↔ F_3` (a symmetry of `W`). Hence, by (a) applied to the family `ι(C_{abc})` (profile `ι x_R`), `g_W`
**does** satisfy the cap at the complementary vertex: `4·ord_e(g_W) ≤ cost_{ι x_R}(g_W)`. If `ι(W)` had equalled `W`,
the identity `cost_{x_R} + cost_{ιx_R} = 2cost_{1/2}` and Theorem B would have closed (b); it does not.

*Proof.* (a) The curves `C_{abc}` (over `C`) form a 2-parameter family sweeping `W`, each with `E_B·C = 1_{B∈F_R}`
and `mult_e C = 4` (two branches of order 2; round 23 Thm 2). If `g_W ∤ f` then `f` does not vanish on a general
`C = C_{abc}`; pulling back to the normalisation `P^1`, the zeros of `f` over `e` number at least `ord_e(f)·mult_e(C)`
and all zeros number at most the poles, `Σ_B π_B(f)(E_B·C) = cost_{x_R}(f)`. (b) `ord_e` and `π_B` are additive
under products (Newton polytopes add), so write `f = g_W^k h` with `g_W ∤ h` and apply (a) to `h`. The restricted
theorem follows from Lemma 3.0 and continuity of `cost_x` in `x`. (c) as stated. `□`

So at the `W`-vertex the cap is **proved for every auxiliary form except multiples of `g_W`**, and it is decided by a
single inequality `cost_{x_R}(g_W) ≥ 4·ord_e(g_W)` (`≥ 64` if `ord_e g_W = 16`, as (d) suggests). Since `deg W > 14`, every small-support
search (§3.4) is blind to this case. I could not compute the Newton polytope of `g_W` (a rough intersection count of `O(4)^3` on the sextic
`Y ⊂ P(1,1,1,1,3)` suggests `deg W ≤ 128`, not checked; interpolation beyond degree 14 exceeds the memory limit).
Random-start numerical intersection counts (`w_count.py`, lower bounds only, visibly incomplete) gave
`#(W ∩ generic line) ≥ 26`, `#(W ∩ line through e, off e) ≥ 13` and `#(W ∩ λC_R) ≥ 79`; they decide nothing.

**Resolution at the `W`-vertex (numerical; §3.4, Prop 3.7).** The question is settled the other way: the
cusps of a `W`-curve can be split. Starting from the exact `W`-curve with `(a,b,c) = (3/2, 5/3, 7/5)` in product form
(14 finite places, one at infinity, cusps at `t = ±√m`), Newton continuation in the 4-point system
`u_j(s_i) = u_5(s_i)` (`i = 1..4`) finds curves of the **same profile** `x_R` passing through `e` at four distinct
parameters with smooth branches; the solution set has dimension 6 = 4 + 2 (affine reparametrisations), the expected
value, along the whole path `ε ∈ [0.01, 0.8]` (`w_deform2.py`); and this 4-dimensional family hits 8/8 random targets
of `T` near it (residuals `≤ 10^{−12}`), while the cusp-constrained (`W`) family hits 0/8 (`w_deform3.py`). So the family
dominates `T`, its general member is not in `W`, and Prop 3.7 gives `4·ord_e(f) ≤ cost_{x_R}(f)` for **every** `f`,
including `g_W`. **The section cap holds at the `W`-vertex; there is no restricted theorem there.** (Certification of
the Newton solutions by interval arithmetic was not done; label: numerical.)

### 3.4 Conjecture WTB (weighted Theorem B on `Π_nt`) and evidence

> **Conjecture WTB.** For every finite `S ⊂ M` and every near-tight split profile `x ∈ Π_nt`:
> `4·reg(S) ≤ Σ_B x_B π_B(S)`. Equivalently: no Liouville/section argument beats `15/4` at any near-tight split
> profile.

Proved: at `x = 1/2` (Theorem B); for centrally symmetric `S` (Theorem 3.2); at `x_R` for `S` not containing the
support of a multiple of `g_W` in the sense of Prop 3.3(a); and trivially whenever `Σ_B (x_B − 1/2)mid_B(S) ≤ 0`.
*Evidence (exact mod `2^31 − 1`; mod-`p` rank is a lower bound for the rational rank, so mod-`p` `reg` is an upper
bound for the true `reg` and cannot produce false violations):*
* random polytopes `Q(p) = {⟨m,1_B⟩ ≥ −p_B}` (`p ∈ {0..3}^30`, and lopsided `p`): 28 712 sets with `|S| ≤ 1500`, LP
  minimum of `cost_x(S)` over `Π_nt` (floating-point LP): maximum of `4reg/mincost` is `1` (to `2·10^{-16}`), attained
  many times (`search_beat.py`, outputs `sb_rand.out`, `sb_lop.out`);
* local search from root zonotopes with random deletions/insertions (`local_search.py`, 2 runs × 600 s): maximum of
  `4reg − mincost` is `0` (to `4·10^{-15}`), attained by sets of size up to 77 (`ls7.out`, `ls8.out`);
* one-sided free polytopes `Q_F(n) = {⟨m,1_B⟩ ≥ −n, B ∈ F}`: for the `W`-class, `reg/n = 3, 3.5, 3.33, 3.5, ≤ 3.4`
  (`n = 1..5`), for the class `1^5 3^{10}` `reg/n ≤ 3.33`; normalised volumes `4!·vol(Q_F(1)) ∈ [89.2, 118.3]` for all 16
  classes, below the `(15/4)^4 = 197.75` that a dimension count would need (`qf_reg.py`, `qf_reg_big.py`, `onesided_vol.py`);
* the quadric-cone form `(u_2−u_1)(u_3−u_1) − (u_4−u_1)(u_5−u_1)` (ord 2): minimum cost over `Π_nt` = `25/3 > 8`
  (`cone_section.py`);
* **k = 4 control** (`k4_wtb.py`): the same statement with `2·reg ≤ cost_x` on the `k = 4` near-tight polytope,
  18 875 random sets, maximum ratio exactly 1.

*Why a proof is not immediate (derived).* Theorem B's slicing step pairs `w_{P'}(S) + w_{P'∪1}(S)` to recover the
widths of two slices; one-sided costs `π_B` give only two of the four extremes, so the induction does not close for
non-symmetric weights. A proof needs a genuinely asymmetric slicing lemma, or the curve criterion below.

**Proposition 3.7 (curve criterion for the cap; proved).** Let `x ∈ Π_nt ∩ Q^{30}` and `N x ∈ Z^{30}`. Suppose there is
an irreducible family of rational curves `C ⊂ LM_5` through `e`, with `E_B·C = N x_B` for all `B` and
`mult_e C ≥ 4N`, whose members cover a dense subset of `LM_5`. Then `4·ord_e(f) ≤ cost_x(f)` for **every** Laurent
polynomial `f ≠ 0`; i.e. WTB holds at `x`. *Proof.* A general member is not contained in `{f = 0}`; on its
normalisation `f` has at least `ord_e(f)·mult_e C` zeros over `e`, and at most `Σ_P h_{supp f}(−v_P) ≤ Σ_B N x_B π_B(f)`
poles (boundary points with valuation vectors `v_P`; sublinearity of the support function). `□`
*Expected dimension (derived).* On `Bl_e LM_5` such curves have `−K·C̃ = 15N − 3·4N = 3N`, so the family of rational
curves of this class has expected dimension `3N + 1 ≥ 4 = dim`: covering families are expected at **every** rational
near-tight profile, and the record's covering inequality (`3E·C̃ ≤ B·C`) is satisfied. The `W`-curves are the
codimension-2 sublocus of the `N = 1`, `x = x_R` family where the four smooth branches through `e` merge into two
cusps with conjugate real-quadratic tangents — exactly the members that carry A-points (`α = 15/(2·2)`); a general
member has four smooth branches and `α ≥ 15/2`. Since WTB is linear in `x`, it suffices to exhibit such a family at
each vertex of `Π_nt` (values in `{0, 1/3, 1/2, 2/3, 1}`, so `N ∈ {1, 2, 3, 6}`).

**Corollary 3.8 (WTB on a 10-dimensional slab through the symmetric profile; numerical at one input).** The dominant
`N = 1` family at `x_R` (§3.3, resolution) and its images under `S_5` and the inversion `ι` (both fix `e` and preserve
multiplicities) give WTB at all 20 code vertices `x_R`, `ι x_R` (`|R| = 3`). Their convex hull `Q_code` has affine
dimension **10** inside the 11-dimensional `Π_nt` and contains `1/2` in its relative interior (uniform weights `1/20`;
exact, `code_hull.py`). By linearity, **every auxiliary Laurent polynomial is capped at `15/4` at every profile in
`Q_code`.** Adding the even-set vertex `x_{[5]} = 1[|B| even]` and its complement would make the hull 11-dimensional
(`code_hull.py`), i.e. a full neighbourhood of the symmetric profile; no starting curve is known there, and random-start
Newton (`cover_family.py`) found no solutions at any one-sided class (a solver limitation: it also failed at the
`W`-class, where solutions exist).

### 3.5 Bounded pair cofactors: integral points, and why integrality adds nothing to sections

**Proposition 3.4 (integral-point model; proved).** Let `S_K = {∞} ∪ {p ≤ K}`. A primitive tuple has all
`|t_ab| ≤ K` only if `[u]` is an `S_K`-integral point of the **affine** 4-fold
`X° = LM_5 ∖ ⋃_{a<b} Δ̄_ab` (`Δ̄_ab` = closure of the subtorus `u_a = u_b`), i.e. every regular function
`z_ab = 1/(1 − u_a/u_b) = −α_ba/t_ab` takes values in `Z[1/S_K]`. For `K = 1` the condition is exactly that every ratio
`u_a/u_b`, in lowest terms, is a ratio of consecutive integers (`(u_a − u_b) | u_a`).
*Proof.* `1 − u_a/u_b = −t_ab/α_ba` in lowest terms, since `gcd(t_ab, α_ba) = gcd(α_ab, α_ba) = 1`. `X°` is affine because
`Σ_{a<b} Δ̄_ab ~ (1/2)Σ_B |B|(5−|B|) E_B` (from `div(χ^{e_a−e_b} − 1)`), the permutohedral divisor
(`f(k) = k(5−k)/2` strictly concave), which is ample on `LM_5`. `□`

So "five divisors with bounded pair cofactors" is an integral-point problem on a fixed affine variety with a fixed
finite `S` — the setting where Siegel/Corvaja–Zannier/Levin arguments (Subspace with **fixed** `S`, so with a
legitimately fixed exceptional set) apply. Lemma 3.5 shows what they can and cannot add.

**Lemma 3.5 (Δ-pole amplification; proved).** Fix a profile `x ∈ Π_nt`. Let `V = Π_{a<b}(u_a − u_b)/Π u_a^2` (the
Vandermonde section: `ord_e V = 10`, `cost_x(V) = 40` for all `x ∈ Π_nt` by Theorem 3.2). For a Laurent polynomial `g`
and `m ≥ 0`, the rational function `f = g/V^m` is regular on `X°` up to toric poles, and at a bounded-residue point
with profile `x` and comparable pair distances (`|u_a − u_b| ≍ max|u_c − u_d|`, automatic when near-tight) the
product formula gives
`(ord_e g − 10m)(−log dist) ≤ (cost_x(g) − 40m) log X + O_{g,m,K}(1)` whenever `f(y) ≠ 0`.
(i) If `4·ord_e(g) > cost_x(g)` for some `g`, then for suitable `k, m` the function `f = g^k/V^m` has
`ord_e f > 0 ≥ cost_x(f)`, so **every** bounded-residue point with profile near `x` and small enough distance to `e`
(no A-point condition at all) lies on `{g = 0}`.
(ii) If `4·ord_e(g) ≤ cost_x(g)` for all `g` (WTB at `x`), then every `f = g/V^m` with `ord_e f > 0` has
`cost_x(f) ≥ 4·ord_e f`: the Δ-poles give nothing.
*Proof.* At `p ∤ K`, `v_p(V(y))` equals its toric value exactly (the cofactors are prime to `p`), so `f`'s effective
pole along `E_B` is `π_B(g) − m·f(|B|)`; at `p | K` the loss is `O(m log K)`; at `∞`,
`|f(y)| ≪ dist^{ord_e g}/Π|y_a − y_b|^m ≍ dist^{ord_e g − 10m}`. For (i) choose `k` large and `m` in
`(k·cost_x(g)/40, k·ord_e(g)/10)`. For (ii), `cost − 40m ≥ 4ord_e g − 40m = 4(ord_e g − 10m)`. `□`

**Consequence.** Bounded residues (integrality at the root arrangement) and Liouville give exactly the same exponent
as Liouville alone, at every profile; and by Theorem 3.2, for symmetric supports this exponent is `15/4` at every
near-tight profile. Integrality would help precisely at profiles where WTB fails — then it helps a lot (no distance
condition). Corvaja–Zannier-type arguments with `S = S_K` and rational linear forms reduce to Lemma 3.5 (the
Subspace Theorem with rational forms and no special `p`-adic behaviour at `p ∈ S_K` is Liouville).

**Theorem 3.6 (cross-product structure of triples; proved, exact check).** For a triple `τ = {a,b,c}` with complement
`{d,e}` put `N_τ(T) = Π_{R⊆{d,e}} n_{T∪R}` (`T ⊆ τ`), `x^τ = (N_τ(a), N_τ(b), N_τ(c))`,
`y^τ = (N_τ(bc), N_τ(ca), N_τ(ab)) ∈ Z^3_{>0}`. Then

    x^τ × y^τ = (t_bc, t_ca, t_ab).

*Proof.* `α_ab = Π_{B∋a,b∉B} n_B = N_τ(a)N_τ(ac) = x^τ_a y^τ_b`, `α_ba = x^τ_b y^τ_a`, and `t_ab = α_ab − α_ba`. `□`
*(exact, `cone_checks.py` (c): on the mnn31 point and on 200 random 5-tuples, computing the `n_B` by factorisation.)*

*Consequences (proved).* (1) For fixed cofactors `t^τ` with `g = gcd(t^τ)`, `x^τ, y^τ` span a sublattice of index
`g` of the rank-2 lattice `Λ = (t^τ)^⊥ ∩ Z^3`; so `(x^τ, y^τ) = (e_1, e_2)·γ` for a fixed basis and
`γ ∈ M_2(Z)`, `det γ = g`, and `[u_a : u_b : u_c] = [γξ_a : γξ_b : γξ_c]` with `ξ_i ∈ P^1(Q)` the fixed slopes of the
rows of `(e_1 e_2)`: **every bounded-residue triple projection is a Möbius image of three fixed rationals**, and the
operations `u ↦ u + kΠ y^τ` (translation by a multiple of the diagonal) and its inverse-dual generate `SL_2(Z)`.
(2) The **Möbius k-tuples** `u_a = P_a Π_{b≠a} Q_b`, `(P_a, Q_a) = (⟨E_a, γ_1⟩, ⟨E_a, γ_2⟩)`, have
`t_ab | det γ·(E_a ∧ E_b)` (bounded), only singleton cuts, and `α = k`; for fixed `E` their Zariski closure is a
**surface** in `T_k` (the diagonal torus of `PGL_2` acts by scaling). For `k = 3` this is all of `T_3`, which explains
`A_3` false; for `k = 5` all bounded-residue Möbius tuples lie on finitely many surfaces and are not A-points.
(3) Near-tightness says that for each of the ten triples the two integer vectors `x^τ, y^τ` (size `X^2` at the
symmetric profile) are parallel to within `X^{−4+O(η)}`. Each triple condition alone is satisfied by the convergents
of an arbitrary real direction; the content of Question B is the coupling of the ten triples through the 30 shared
`n_B`. I found no way to exploit this coupling beyond Lemma 3.0.

### 3.6 Phase 0 (Erdős #887 proper): the forum reduction as a family of general-type surfaces (derived)

At `X = √n` put `s_j = d_j + n/d_j`, `z_j = d_j − n/d_j` (integers), `k_j = s_j − s_1 ∈ [0, C² + 1]`. Then
`s_j² − z_j² = 4n` for all `j`, i.e.

    z_j² − z_1² = 2k_j s_1 + k_j²      (j = 2..M),

the forum's "linear polynomials simultaneously squares" (`x = s_1`, `y = z_1²`). For fixed `k = (k_1..k_M)` the
projective closure `S_k ⊂ P^{M+1}` is a complete intersection of `M − 1` quadrics, a surface with
`K_{S_k} = O(2(M−1) − M − 2) = O(M − 4)`: a quartic del Pezzo for `M = 3`, a **K3 surface for `M = 4`**, and of
**general type for `M ≥ 5`** (when smooth enough; not checked). Since `k_j ≤ C² + 1`, for fixed `C` only finitely
many surfaces occur, so **#887 `parts.ii` with constant `K` ⟸ for every `(K+1)`-tuple `k`, `S_k` has only finitely
many integral points with `n > 0`** (this per-template finiteness is legitimate for #887 because `n_0(C)` may depend on
`C`; it is *not* usable for the torus problem, whose phase drifts). The two known five-divisor families are Pell
conics on general-type surfaces: m-czech `k = (0, 224, 225, 324, 900)` (`= offset²`; offsets `0, 2√56, 15, 18, 30`),
mnn31 `k = (0, 3, 4, 9, 12)` (offsets `0, √3, 2, 3, 2√3`) (exact, `phase0_k.py`: the identity holds and `k` is
constant along 4 Pell points of each family). Lang's conjecture (rational points of a general-type
surface lie on finitely many curves) plus Siegel would reduce absolute #887 to: *for some `K`, no `S_k` with `K+1`
entries contains a rational curve with at most two (conjugate real-quadratic) points at infinity.* This is the
phase-0 shadow of the record's "`W`-question"; the 4-divisor case is a K3 density question (Erdős–Rosenfeld's
infinite family). I make no claim beyond this conditional reading.

### 3.7 The 4+1 decomposition (heuristic, with exact bookkeeping)

Removing index 5: the projection to `{1,2,3,4}` of a near-tight balanced 5-tuple is a near-tight balanced 4-tuple at
the `k = 4` Liouville exponent `7/2` (`H_4 = X^{14}`, `dist = X^{−4}`; proved by the bookkeeping of §1), and those are
Zariski dense (gauss Thm 4.1 in the twisted form). Given the 4-tuple (cut numbers `M_S = n_S n_{S∪5}`), the fifth
point needs a split `M_S = n_S·n_{S∪5}` and integers `n_5, n_{[4]} ≍ X^{1/2}` with
`n_5/n_{[4]} = ρ(1 + O(X^{−4}))`, `ρ = Π_{S∋1}M_S/Π_S n_{S∪5}` a rational of height `≍ X^7`: one Diophantine
approximation of a rational by a rational of height `X^{1/2}` at precision `X^{−4}`, probability `X^{−3}`, against a
supply of `X^2` four-tuples per scale: `X^{−1}` in total, the coordinator's count. In the Pell families the large
partial quotient of `ρ` is produced by the unit of the real-quadratic field.

## 4. Pre-registered tests T1–T5

**T1 (fake test; growth claims).** No growth claim is made: `M ≤ (2/3+ε) log R/loglog R` stands. The barrier
statements (Thm 3.2, Prop 3.3, Lemma 3.5, WTB) are statements about which auxiliary forms *fail*; they are consistent
with fakes, which satisfy every section inequality at the cube profile (record Cor C3). Theorem 3.1 concerns curves
over `Q` and is invisible to fakes. Nothing here is an input that fakes could violate, so T1 is vacuous and passed.

**T2 (`k = 4` negative control).** 
* *Barrier statements.* Their `k = 4` analogues hold and are consistent with `A_4(δ < 1/14)` being false: Theorem B
  with `n = 4` caps sections at `7/2`; the cost identity and Theorem 3.2 hold verbatim on the `k = 4` near-tight
  polytope (same proof); and the `k = 4` WTB was tested on 18 875 random sets with maximum ratio exactly 1
  (`k4_wtb.py`). The `k = 4` analogue of Prop 3.7 (curves of profile `Nx` through `e` with multiplicity `2N`,
  `−K_{Bl}·C̃ = 7N − 4N = 3N`, expected dimension `3N ≥ 3 = dim`) likewise predicts the cap `7/2`. A capped section
  method cannot exclude the dense `k = 4` families, as it must not.
* *The only exclusion mechanism in this write-up* is conditional: Lemma 3.5(i)/Prop 3.3(b) ("if WTB fails at `x`,
  then points with profile near `x` are degenerate"). At `k = 4` the dense families of gauss Thm 4.1 are twisted, hence
  at the symmetric profile, where WTB is Theorem B; so the conditional mechanism does not fire on them. Correct
  behaviour.
* *Phase 0.* The surface `S_k` is K3 for four divisors (density allowed: Erdős–Rosenfeld's infinite family) and of
  general type only from five divisors on; the conditional Lang reading (§3.6) therefore says nothing at `k = 4`.
  Correct behaviour.

**T3 (split test).** Nothing here proves split Question B, so absolute #887 is not claimed. Every statement is
compatible with the 2026 constructions: the m-czech curve (`α = 15/4`, profile `x_R`) lies on `W`, and Prop 3.3(a)
shows that no auxiliary form outside `g_W·Z[χ]` can exclude its points (and by Cor 3.8 no form at all beats `15/4`
at its profile); the mnn31 curve (`α = 7/2`) lies on the cone
`𝒬` (exact), which Theorem 3.1 adds to every admissible exceptional set. Both land in the exceptional set; neither is
excluded by anything claimed.

**T4 (criticality test).** My results live inside the section class and say where it stops:
* *Lemma K1 / `τ_e = 10/3` / balanced `15/4`.* The record's balanced cap was proved at the symmetric profile only
  (Theorem B). For split points the profile is 11-dimensional and the section inequality depends on it (Lemma 3.0).
  So a section mechanism at an asymmetric split profile would **not** be one of the mapped caps. Theorem 3.2 shows
  that for all forms with centrally symmetric support (all binomial and Vandermonde products) the cap is nevertheless
  `15/4` at every near-tight profile; Prop 3.3 and Prop 3.7 with the numerical family show the same at the `W`-vertex
  for every form, hence on the 10-dimensional slab `Q_code` through `1/2` (Cor 3.8); WTB conjectures it in general.
  So this lane does not escape the cap; it closes the corner (asymmetric split profiles) that the record had not
  examined, except for profiles off `Q_code` approached by non-symmetric forms.
* *Class `L` (round 24).* Lemma 3.5 adds integrality at the root arrangement (poles along `Δ̄_ab` paid by bounded
  residues) and shows it multiplies only what single sections already give.
* *Multi-point / Thue / Dyson.* Not used. The 4+1 decomposition (§3.7) is a heuristic count, not an argument.
* *Fixed-`S` Subspace.* Used only with `S = {∞} ∪ {p ≤ K}` genuinely fixed (Prop 3.4), where it reduces to
  Lemma 3.5; no exceptional set is taken uniform over varying prime sets.

**T5 (calibration).** Nothing excludes existing points: round-24's real twisted family is at the symmetric profile
(Theorem B, consistent); split Pell curves at `15/4` are protected by Prop 3.3(a); the cone points at `7/2` are put
into `Z`, not excluded; the record's 8- and 10-point families are untouched (no growth claim).

**Handoff constraints.** No uniform exceptional set from a fixed-`S` Subspace theorem (only `S = S_K` fixed with `K`
fixed); no Gaussian gcd or denominator cost discarded (tuples primitive, exact toric heights, `t_ab` exact); no drifting
phase set to zero (the phase-0 material of §3.6 is stated for #887 at `X = √n` only, where the phase is genuinely 0,
and is explicitly not used for the torus problem); no template-dependent threshold treated as uniform (in §3.6 the
per-template finiteness is legitimate for #887 because `n_0(C)` may depend on `C`, and nothing is claimed for the
torus problem from it).

## 5. Which steps use splitness; Gaussian analogue for `T_5`

| result | uses splitness? | twisted (`T_5`) analogue |
|---|---|---|
| Thm 3.1 (cones at `7/2`) | yes: the construction lives on the split torus; the cone is a split `Q`-subvariety | the real structure `u ↦ 1/ū` moves `𝒬` to a cubic; only `𝒬 ∩ σ𝒬` (a surface) could be twisted-rational. Second instance (after parity) of the twisted exceptional set being smaller |
| Lemma 3.0, Thm 3.2 | uses the 11-dim split profile space | twisted near-tight profiles are only `x ≡ 1/2`; there the statement is Theorem B (record). The split profile freedom is exactly what the twist removes, and Thm 3.2 shows it gives symmetric forms nothing |
| Prop 3.3 | yes: the `W`-vertex has weight-1 cuts | none: parity lemma (twisted curves have even weight on every cut) |
| Prop 3.4, Lemma 3.5 | no | identical over `Z[i]`: bounded residues = integrality at the Gaussian root arrangement; the Vandermonde has ratio `15/4` at `x = 1/2`, so Δ-poles add nothing |
| Thm 3.6 (triples) | integer factorisation of `α_ab` | over `Z[i]` the same identity holds with `x^τ, y^τ ∈ Z[i]^3`, and conjugation relates `y^τ` to `x̄^τ` (cut sides swap), giving `GL_2(Z[i])`-structure with a reality condition; not pursued |
| §3.6 phase 0 | #887 only | the near-axis bound (record lit18 §2.3) |

**Decision.** Nothing in this lane ports to a proof of Question B for `T_5`. What ports is negative and precise: the
two split-only freedoms (profile asymmetry; integrality of residues) are worthless to section methods for symmetric
supports, and the twisted problem sits exactly at the symmetric profile where the cap is already a theorem. The split-only
positive possibility this lane looked for — sections beating `15/4` at an asymmetric profile, e.g. via `g_W` at the
`W`-vertex — is closed at the `W`-vertex (Cor 3.8, numerical) and could only survive off `Q_code`; even there it could
not port to `T_5`, whose only near-tight profile is the symmetric one.

## 6. Scripts and certificates

All under `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/split_checks/`; runs under
`prlimit --as=3500000000`. The repository was not edited.

| file | checks | label |
|---|---|---|
| `dl/Erdos887.lean` | mnn31 Lean file (downloaded as data, read, not executed) | source |
| `czech72.py` | `C = 7/2` family along 40 Pell points: 5 divisors, `C_act → 2√3`, `γ → 7/2`, `log(dist^{7/2}H)` periodic | exact integers |
| `fam72.py` | second factorisation identity; exact Jacobian rank 3 of the 4-parameter family | symbolic / exact |
| `cone_checks.py` | (a) cone identity at 6 Pell points; (b) m-czech points on none of the 15 cones; (c) triple identity `x×y = t` on 201 tuples; (d) Galois image of the cone | exact |
| `phase0_k.py` | phase-0 `k`-tuples and the identity `z_j² − z_1² = 2k_js + k_j²` on both families | exact |
| `profiles.py` | `dim Π_nt = 11`, 787 vertices found (278 one-sided), binomial cost `≥ 4` at all of them (`|m_a| ≤ 3`; the general proof is Thm 3.2) | exact rank; LP vertices (float) |
| `onesided_enum.py` | exactly 278 one-sided near-tight profiles (`X_c ∈ {5,…,10}`); near-tight code families = `F_R`, `|R| ∈ {3,5}` | exact |
| `reglib.py` | regularity index mod `2^31−1` (checked: `reg(nR_5) = ⌊10n/3⌋`, `n ≤ 4`) | exact mod `p` |
| `search_beat.py` → `sb_rand.out`, `sb_lop.out` | 28 712 sets: `max 4reg/min_{Π_nt} cost = 1` | mod-`p` reg, float LP |
| `local_search.py` → `ls7.out`, `ls8.out` | local search from zonotopes: `max(4reg − mincost) = 0` | same |
| `onesided_vol.py`, `qf_reg.py`, `qf_reg_big.py` → `qf_reg_big.out` | one-sided free polytopes: volumes, `reg(nQ_F)/n ≤ 3.5` (`n ≤ 5`) | exact mod `p`; float volumes |
| `cone_section.py` | quadric-cone form: min cost over `Π_nt` = `25/3 > 8` | mod-`p`, LP |
| `k4_wtb.py` → `k4_wtb.out` | `k = 4` control: 18 875 sets, max ratio 1 | mod-`p`, LP |
| `w_equation.py` | no form of degree `≤ 14` vanishes on `W` (3100 points over `F_{1000003}`) | exact mod `p` |
| `tangent_W.py` | branch directions of `W` at `e` span the irreducible octic `Σ_W` | symbolic |
| `w_member2.py` → `w_member2.out` | numerical membership in `W`: controls `≤ 1.3·10^{-15}`; inversion images `1.4·10^{-1}`, `5·10^{-2}`; swap `F_2↔F_3` `2.8·10^{-16}` | numerical |
| `w_count.py` → `wc_line.out`, `wc_linee.out`, `wc_curve.out` | random-start intersection counts (lower bounds 26, 13, 79; incomplete) | numerical, weak |
| `w_deform.py`, `w_deform2.py` → `w_deform2.out` | product form of an exact `W`-curve; continuation to four-branch curves of profile `x_R` through `e`; solution-set dimension 6 = 4 + 2 along `ε ∈ [0.01, 0.8]` | numerical |
| `w_deform3.py` → `w_deform3.out` | dominance: four-branch family hits 8/8 random targets (residual `≤ 10^{-12}`), cusp-constrained `W`-family 0/8 | numerical |
| `cover_family.py` | random-start search for such families at the 16 one-sided classes: no solutions found (solver limitation; it also fails at the `W`-class) | numerical, negative |
| `cover_vertices.py` | the same for all vertex classes (`N ∈ {1,2,3,6}`); written, **not run** | — |
| `code_hull.py` | the 20 code vertices span an affine 10-space with `1/2` in the relative interior; with `x_{[5]}`, `ι x_{[5]}`: 11 | exact rank |
| `../split_notes.md` | working notes N1–N7 | — |

## 7. Open questions produced by this lane

* **Q-split-1 (answered numerically: yes).** Does `g_W` satisfy the cap at its own vertex `x_R`? By Prop 3.7 and the
  dominant four-branch family obtained by splitting the cusps of `W`-curves (§3.3), yes. *To do:* certify one family
  member and the rank-4 differential by interval Newton (Krawczyk) so that Cor 3.8 becomes a theorem.
* **Q-split-1′ (the remaining vertices).** Exhibit dominating families of rational curves through `e` with profile
  `N v` and multiplicity `4N` at the other vertices of `Π_nt` (first: the even-set vertex `x_{[5]}`, which together with
  `Q_code` spans all of `Π_nt`). Positive answers at all vertex classes prove WTB, i.e. close the section method on
  split near-tight points completely. Random-start Newton is not good enough; a start curve (a Pell-type construction
  with profile `x_{[5]}`, or a degeneration) is needed.
* **Q-split-2 (WTB).** Prove `4·reg(S) ≤ Σ_B x_B π_B(S)` on `Π_nt`, by curves (Q-split-1′) or by an asymmetric
  slicing lemma (§3.4).
* **Q-split-3.** Is `K_{M̄_{0,7}} + D_int` big (log general type of `LM_5 ∖ ⋃Δ̄_ab`)? If not, Lang–Vojta does not even
  predict degeneracy of bounded-residue points without the A-point condition.
* **Q-split-4 (phase 0).** Does any complete intersection `S_k` (six entries) contain a rational curve with two conjugate
  real-quadratic points at infinity (a six-divisor Pell family)?
