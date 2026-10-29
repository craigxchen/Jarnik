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
  symmetric support (Theorem 3.2, Lemma 3.5), for all `f` at the symmetric profile (Theorem B, record), and for all
  `f` not divisible by `g_W` at the `W`-vertex (Proposition 3.3). In general this is a precise conjecture (WTB, §3.4)
  with exact evidence; at the `W`-vertex it reduces to **one inequality** for the defining polynomial of `W`;
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
