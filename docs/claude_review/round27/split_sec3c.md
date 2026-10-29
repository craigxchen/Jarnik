### 3.5 Bounded pair cofactors: integral points, and why integrality adds nothing to sections

**Proposition 3.4 (integral-point model; proved).** Let `S_K = {∞} ∪ {p ≤ K}`. A primitive tuple has all
`|t_ab| ≤ K` only if `[u]` is an `S_K`-integral point of the **affine** 4-fold
`X° = LM_5 ∖ ⋃_{a<b} Δ̄_ab` (`Δ̄_ab` = closure of the subtorus `u_a = u_b`), i.e. every regular function
`z_ab = 1/(1 − u_a/u_b) = −α_ba/t_ab` takes values in `Z[1/S_K]`. For `K = 1` the condition is exactly that every ratio
`u_a/u_b`, in lowest terms, is a ratio of consecutive integers (`(u_a − u_b) | u_a`).
*Proof.* `1 − u_a/u_b = t_ab/α_ba` in lowest terms, since `gcd(t_ab, α_ba) = gcd(α_ab, α_ba) = 1`. `X°` is affine because
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
dividing `g` of the rank-2 lattice `Λ = (t^τ)^⊥ ∩ Z^3`; so `(x^τ, y^τ) = (e_1, e_2)·γ` for a fixed basis and
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
mnn31 `k = (0, 3, 4, 9, 12)` (offsets `0, √3, 2, 3, 2√3`). Lang's conjecture (rational points of a general-type
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
