# Round 20 — gauss.md: circle-specific Gaussian arithmetic at small k

Status: complete draft, unrefereed (round 20, lens: circle-specific Gaussian arithmetic at small k).

## 0. Summary and verdict

**Verdict: `reduced_to_cap` with one new reduction and one new refutation.** No growth improvement and no
uniform bound. The lens "circle-specific Gaussian arithmetic at small `k`" gives the following.

* **(Reduction; proved, Theorems 3.4, 3.4′.)** Dividing a random `k`-tuple of a short-arc cluster by its
  own Gaussian gcd (norm `≈ N^{2^{1-k}}`) puts it on an arc `R''^{1/2-1/(2^k-2)}` of its reduced circle,
  strictly below the critical exponent. Hence the uniform theorem follows from **`A_k(δ)`**: *primitive
  `k`-tuples of lattice points on arcs of length `R^{1/2-δ}` lie, for large `R`, in a fixed proper
  Zariski-closed set of similarity classes*, for one `k` and one `δ < 1/(2^k-2)` (`δ ≤ 2^{-k-1}` with a
  completely elementary proof). It also follows from the balanced-residue version `(L_k^circ)` through
  Steps 0–5 of `ptolemy.md` §6.2 (checked step by step, §3.2; Step 2 must be run in the similarity torus
  `T_k`, not in `X_k`).
* **(a) (§2.)** The `(L_5)`/`(L_6)` families are not circle data: read with isotropic points `±i` they
  satisfy essentially none of the isotropy/norm divisibilities (`φ ≈ 0`, exact), and their three
  mechanisms (density, one rational curve, covering genus-one family with positive rank) provably do not
  transfer to circles at `k ≥ 5` (Theorem 4.3).
* **(b) (§3.)** The weakest circle statements that still imply the uniform theorem: `(L_5^circ)`
  (`A_5(δ<1/30)`), `(L_6^circ)` (`A_6(δ<1/62)`), `(L_7^circ)` (`A_7(δ<1/126)`); they are logically
  incomparable, with increasing heuristic margin.
* **(c) (§4.)** `k = 4` is **refuted** (Theorem 4.1: an explicit rational 2-parameter family of Gaussian
  linear blocks, dense in `T_4`, of balanced bounded-residue 4-tuples on arcs `≍ R^{3/7}`), so five points is
  the minimum. `k = 5, 6, 7`: **open**. They are predicted by a non-vacuous instance of Vojta on a 4-fold
  (`Bl_e` of the twisted Losev–Manin space, `D = E`, `S = {∞}`; Theorem 4.2, conditional), are immune to
  the curve-family refutation mechanisms (Theorem 4.3: Siegel for genus `≥ 1`, Roth plus log-free
  deformations for bounded-degree rational families), and are out of reach of every Gaussian tool in the
  record. Those tools (norm equations, isotropy, units, item-532 nested-chain identities, Gaussian Ptolemy)
  are auxiliary polynomials with forced divisibility, and `construct.md`'s mean-width theorem makes the
  circle class `γ̂` movable for every `k`. So Liouville is tight **exactly** at the needed exponent `7/15`.
  The target point is rational, so Roth/Schmidt cannot help. Status: same as `(L_8)`, at five points
  instead of eight.
* **What changed in the barrier map.** Item 3 of the checklist extends to the full Gaussian model
  (`γ̂` movable for all `k`, previously checked only at `k = 4`). Item 7's analogue is proved for circles at
  `k = 5`. The Ptolemy threshold `8` drops to `5` once the Gaussian structure and the tuple gcd are used. The
  circle structure buys three points of `k`, and no method.


## 1. Setup and what the Ptolemy system forgets

**Conventions.** *Theorem/Proposition/Lemma* = proved here with all quantifiers, or cited with the
exact place of the proof. *Exact* = finite exact computation (Python integers, `Fraction`, sympy
rational arithmetic) in `gauss_checks/`. *Evidence* = a computation no proof uses. *Heuristic* = no
proof. *Conditional* = proved from a named hypothesis. Nothing here has been refereed.

### 1.1 Objects

* `X_k = {x_1^2+y_1^2 = ... = x_k^2+y_k^2} ⊂ P^{2k-1}` (the equal-norm variety of
  `growth/conditional.md`), `u_j = x_j+iy_j`, `v_j = x_j-iy_j`.
* `T_k` := the **similarity torus**: `T_k(Q) = {(w_1,...,w_k) : w_j ∈ Q(i), N(w_j) = 1}/diagonal`, i.e.
  `(R^1_{Q(i)/Q} G_m)^k / R^1_{Q(i)/Q} G_m`, a `(k-1)`-dimensional torus over `Q` (non-split; it splits
  over `Q(i)` as `G_m^k/G_m` via `w ↦ u`). A `k`-tuple of distinct points `z_1..z_k` of one circle
  `|z|^2 = N` defines `[z] := (z_j/z_1)_j ∈ T_k(Q)`; `[z]` is invariant under every similarity
  `z ↦ γz`, `γ ∈ Q(i)^*`. Over `C`, the map `C_N^k → T_k`, `(z_j) ↦ (u_j/u_1)` is surjective onto the
  torus `(C^*)^k/C^*`.
* `T_k` is the open torus of the **twisted Losev–Manin space** `L̄_k^tw` (the toric variety of the
  permutohedral fan in `Z^k/Z·1`, whose rays are the `1_A`, `∅ ≠ A ⊊ [k]`; complex conjugation acts on
  cocharacters by `-1` and maps the ray `1_A` to `1_{A^c}`, so the fan is Galois-stable and `L̄_k^tw` is a
  smooth projective toric variety over `Q`). Blowing up `L̄_k^tw` at the identity `e` and along the
  diagonal subtori gives the twisted form `M̄^tw_{0,k+2}` of `M̄_{0,k+2}` (marks `1..k` plus the two
  isotropic points `P = -i`, `Q = +i`, swapped by Galois): this is `threshold/L8.md` §6.1's setting.

### 1.2 Dictionary (this is the "full Gaussian model")

Let `ζ_1..ζ_k` be a **primitive** tuple (`gcd_{Z[i]} = 1`) of distinct lattice points of `|ζ|^2 = N'`.
By `two_thirds/proof.md` Lemma 1.1, `N'` is odd with only split prime factors and
`ζ_j = ε_j ∏_p π_p^{a_j(p)} π̄_p^{e_p - a_j(p)}`. Every threshold layer `(p,t)`, `1 ≤ t ≤ e_p`, gives a
nontrivial cut `{A, A^c}` of `[k]` (`A = {j: a_j(p) ≥ t}`; nontrivial because the tuple is primitive).
Put `n_{A} = ∏_{layers with cut {A,A^c}} p` (the Gaussian **block norm** of the cut) and
`W' = log N' = Σ_cuts log n_A`.

| circle datum | divisor on `Bl_e L̄_k^tw` / `M̄^tw_{0,k+2}` | height at `[ζ]` |
|---|---|---|
| Gaussian block of the cut `{A,A^c}` (norm `n_A`) | toric pair `E_A + E_{A^c}` = `D_{P∪A} + D_{P∪A^c}` | **exactly** `log n_A` (finite), archimedean part `0` (Lemma 1.1 below) |
| singleton cuts `{j}|rest` and, in anchored form, the top block `T=[m]` | `D_{P j}, D_{Q j}` (a `D_2`!) | `log n_{{j}}` |
| short arc of length `ℓ` | `E_e` = `D_{PQ}` (all `k` marks collide away from `±i`) | archimedean `λ_{E_e,∞} = log(R'/ℓ) + O(1)`; finite part = common `p`-adic collisions not through `±ρ_p` (≤ residues) |
| residues `t_ij`, `Y_j` (`ptolemy.md` §1) | diagonal divisors `D_B`, `B ⊆ [k]`, `2 ≤ |B| ≤ k-1`, and finite part of `D_{PQ}` | `O(H_Σ)` |

**Lemma 1.1 (toric heights are block norms; proved).** For the standard toric Weil functions on
`L̄_k^tw`, `h_{E_A + E_{A^c}}([ζ]) = log n_A` exactly, and `h_{-K_{L̄}}([ζ]) = Σ_A h_{E_A} = W'`.
*Proof.* `[ζ] = (ζ_j/ζ_1)_j`. At the place `π | p`, `v_π(ζ_j/ζ_1) = a_j(p) - a_1(p)`, so the valuation
vector is `(a_j(p))_j mod Z·1 = Σ_{t=1}^{e_p} 1_{A_{p,t}}`, a sum over a **nested chain** of cuts, i.e. a
point of one cone of the permutohedral fan. The toric local height of the ray divisor `E_ρ` at a vector of
the cone spanned by a chain is `log p` times the coefficient of `v_ρ` (Courant function), i.e. the number
of layers with `A_{p,t} = A`. At `π̄`, `v_π̄(ζ_j/ζ_1) = -(a_j - a_1)`, giving the complementary chain
`Σ 1_{A^c_{p,t}}`. Summing over `π, π̄` gives `log n_A` for `E_A + E_{A^c}`. At `∞`,
`|ζ_j/ζ_1| = 1`, so every archimedean toric local height is `0`. `□`

So the **Ptolemy/`M̄_{0,k}` data forget** (cf. `growth/ptolemy.md` Theorem B(iii)): the singleton and
top blocks (the `D_{P j}`), the isotropy `ρ_T^2 = -1` (the existence of the conjugate pair `±i` in every
`p`-adic tree, which forces every `p`-adic cluster of the `k` points to be a cluster **around `±ρ_p`**),
and the archimedean contact with `D_{PQ}` (the short arc itself, which in `M̄_{0,k}` is invisible: `k`
points close together have bounded cross-ratios). These are exactly the three extra rows of the table.
`threshold/L8.md` §6.1 introduced the balanced class `γ̂` of this model (`γ̂.D_{PQ} = 1/4`,
`γ̂.D_{P∪A} = 2^{-k}`); everything below is consistent with it (`gauss_checks/kclass.py`).

### 1.3 The one structural fact used repeatedly: the tuple gcd shortens the arc

**Lemma 1.2 (gcd reduction; proved).** Let `z_1..z_k` be distinct lattice points of `|z| = R` with all
pairwise distances `≤ L`, let `g = gcd_{Z[i]}(z_1..z_k)`, `ζ_j = z_j/g`, `R'' = R/|g|`. Then the `ζ_j` are
distinct, primitive, on `|ζ| = R''`, with pairwise distances `≤ L/|g|`; `[ζ] = [z]` in `T_k`; and
`R'' ≥ R/L` (so `R'' ≥ √R` if `L ≤ √R`). If `L ≤ √R`, `log|g|^2 = W_∅` and `W = log R^2`, then the reduced
tuple has all pairwise distances `≤ R''^{1/2-δ}` whenever `δ ≤ W_∅/(2(W - W_∅))`; in particular whenever
`W_∅ ≥ θW` and `δ ≤ θ/2`.
*Proof.* Distinctness, primitivity and `[ζ]=[z]` are clear. `|g| ≤ |z_a - z_b| ≤ L` for one pair, as
`z_a - z_b = g(ζ_a - ζ_b)` with `ζ_a ≠ ζ_b` Gaussian integers; so `R'' = R/|g| ≥ R/L`. Finally
`L/|g| ≤ e^{W/4 - W_∅/2}` and `R''^{1/2-δ} = e^{(W-W_∅)(1/4-δ/2)}`; the inequality between exponents is
`(W-W_∅)δ/2 ≤ W_∅/4`. `□`

For a balanced tuple taken from a large cluster, `W_∅ ≈ 2^{1-k}W` (two of the `2^k` sign patterns), so the
reduced tuple sits on an arc `R''^{1/2 - 1/(2^k-2)}`: **strictly below the critical exponent `1/2`**. At
`k = 3` this is `R''^{1/3}`, exactly Cilleruelo's sharp three-point scale `(16R)^{1/3}`; at `k = 5` it is
`R''^{7/15}`, against the Cilleruelo–Córdoba five-point threshold `R''^{2/5}`.

## 2. (a) Do the (L_5)/(L_6) counterexample families come from circles?

**Answer: no, in three increasingly strong senses.**

**2.1 As given (isotropic points at `±i`), they violate isotropy and the norm condition almost
completely (exact, `gauss_checks/familycheck.py`).** Circle data require `n_T | N(P_j) = X_j^2+Y_j^2` for
every visible block `T ∋ j` (this *is* isotropy: `X_j ≡ ρ_T Y_j mod n_T`, `ρ_T^2 ≡ -1`). Measure
`φ = min_j log gcd(N(P_j), ∏_{T∋j} n_T)/log ∏_{T∋j} n_T` (circle data need `φ = 1`).
* `ptolemy.md` C1 (`k=5`, `Y_i = 1`): every resultant `Res_s(n_T(s), X_i(s)^2+1)` is nonzero (max
  `|Res| = 41408180101`), so `n_T(s) | X_i(s)^2+1` fails as soon as `|n_T(s)|` exceeds it; on the explicit
  class `s ≡ 52346140003 mod 85900454640`, `φ = 0.0000`.
* `threshold/L8.md` `E_t` (`k=6`, `t_0 = (2,3,-5,1)`, the 13 anchored points `N = 2..14`, `w` up to `495`):
  `max_j φ ≤ 0.027`, and `≤ 0.002` for `N ≥ 4`.
* `L8.md` Prop C (`k=5`, integers with pairwise coprime differences): circle data would need
  `x_i - x_j | x_i^2+1` for all pairs, a density-zero condition (it forces every prime of every
  difference to be `≡ 1 mod 4`).

**2.2 Any circle model.** A form-blind configuration `Φ ∈ M_{0,k}(Q)` has a 2-parameter family of circle
models (choose the conjugate pair `ρ, ρ̄ ∈ P^1(Q(i))`; with `P_0 = (1,0)`, `ρ = b + ai`, `P_j ↦ (aX_j+bY_j, Y_j)`
up to scaling). A circle model of a balanced bounded-residue member, if balanced over all cuts, is a
counterexample to `(L_k^circ)`; so "can the families come from circles" is exactly "do they contain
`(L_k^circ)`-counterexamples", and the forgetful image of the circle class is the form-blind one
(`π_*γ̂ = 2^{1-k}β`, `L8.md` §6.1). The single curve C1 cannot matter at all: its preimage in `T_5` is a
proper closed subset (dimension `3 < 4`) and can be put into `Z`.

**2.3 The refutation mechanisms provably do not transfer (Theorem 4.3 below).** The `(L_5)`/`(L_6)`
refutations work by *density* (Prop C), by *one rational curve with infinitely many points* plus
re-choosing the exceptional set (C1), and by *a covering family of genus-one curves each carrying
infinitely many points* (`E_t`). For circle data at `k ≥ 5`:
* a fixed curve of genus `≥ 1` carries only finitely many circle points with subcritical arcs (Siegel;
  Theorem 4.3(i)) — the `E_t` mechanism is dead;
* no family of rational curves of bounded degree whose members each carry infinitely many such points
  covers `T_k` (log-freeness plus Roth; Theorem 4.3(ii)) — the C1 mechanism, run over a covering family,
  is dead;
* the density count of balanced bounded-residue circle configurations is `e^{-w}` per scale at `k = 5`
  (against `e^{+6w}` form-blind), §4.2 (**heuristic**).

What separates the circle from the `(L_5)`/`(L_6)` examples is precisely the archimedean contact with
`D_{PQ}` (the short arc) together with the toric (Gaussian-block) heights; in the language of §4.2,
`(K + D_{PQ}).γ̂` changes sign between `k = 4` and `k = 5`, while the form-blind `K.β` changes sign
between `7` and `8`.

## 3. (b) The circle-specific statements (L_k^circ), k = 5, 6, 7, and the implication through Steps 0-5

### 3.1 Two statements

**Definition (`A_k(δ)`, few points on a subcritical arc, off a closed set).** `A_k(δ)` asserts: there are
a proper Zariski-closed `Z ⊂ T_k` and `R_0` such that every **primitive** `k`-tuple of distinct lattice
points of a circle of radius `R ≥ R_0`, with all pairwise distances `≤ R^{1/2-δ}`, has `[ζ] ∈ Z`.

**Definition (`(L_k^circ)`).** For a primitive tuple `ζ_1..ζ_k` with cut norms `n_A` (§1.2) put
`d_ab = Σ_{cuts separating a,b} log n_A` and define `t_ab` by `|ζ_a - ζ_b| = 2|t_ab| R'' e^{-d_ab/2}`.
Then `t_ab` is a nonzero integer: it is the primitive imaginary residue of the edge `ab`
(`growth/ptolemy.md` Theorem B(i); item 532, eq. (9) and the anchor-independence paragraph), equal to `Y_j`
on anchor edges and to `t_ij` otherwise. Put `H_Σ = Σ_{a<b} log|t_ab|`. `(L_k^circ)` asserts: there are
`c > 0`, `η_0 > 0`, `w_0` and a proper Zariski-closed `Z ⊂ T_k` such that every primitive `k`-tuple with
`[ζ] ∉ Z` and `|log n_A - w| ≤ η_0 w` for all `2^{k-1}-1` cuts, `w ≥ w_0`, has `H_Σ ≥ c w`.

`(L_k^circ)` is `ptolemy.md`'s `(L_k)` restricted to actual circle tuples, with three changes, each in the
direction of a *weaker* statement: (i) the data are actual Gaussian tuples (norms, isotropy, singleton and
top blocks, nested chains, conjugate-primitivity, units are all present, nothing is relaxed); (ii) the
exceptional set may be any proper closed subset of the `(k-1)`-dimensional `T_k`, not only of the
`(k-3)`-dimensional `M_{0,k}`; (iii) it is required only at balanced profiles over **all** `2^{k-1}-1` cuts
(including the singletons, which Ptolemy cannot see). Every hypothesis Steps 0–5 supply is kept; none is
added that they do not supply (the constraints of the task: no uniform exceptional set from fixed-`S`
Subspace, no discarded gcd, no drifting phase replaced by zero, no template-dependent threshold — `Z`,
`c`, `η_0`, `w_0` are absolute).

**Proposition 3.1 (`A_k ⇒ (L_k^circ)`; proved).** If `A_k(δ)` holds for some `δ < 1/(2^k-2)`, then
`(L_k^circ)` holds with the same `Z`.
*Proof.* Let `m = k-1`. Each pair is separated by exactly `2^{m-1}` of the `2^m-1` cuts, so under balance
`d_ab ≥ 2^{m-1}(1-η_0)w` and `W' = log N' ≤ (2^m-1)(1+η_0)w`. If `H_Σ < cw` then every `|t_ab| < e^{cw}`
and `|ζ_a-ζ_b| < 2 e^{cw} R'' e^{-2^{m-2}(1-η_0)w}`. Since `R'' = e^{W'/2}`, the exponent of `R''` is at most
`1 - 2^{m-1}(1-η_0)/((2^m-1)(1+η_0)) + 2c/((2^m-1)(1-η_0))`, which is `1/2 - 1/(2^k-2) + O(η_0 + c)`. Choose
`η_0, c` so small that it is `≤ 1/2 - δ'` with `δ < δ' < 1/(2^k-2)`; then for `w ≥ w_0` the factor `2` is
absorbed and all distances are `≤ R''^{1/2-δ}`. `A_k(δ)` gives `[ζ] ∈ Z` or `R'' < R_0`; the latter bounds
`w`. Enlarging `w_0` finishes. `□`

### 3.2 The implication into the uniform theorem

**Lemma 3.2 (averaging the tuple gcd; proved).** Let `z_1..z_M` (`M ≥ 2k`) be a normalised cluster
(`gcd_{Z[i]} = 1`), `W = log N`, and for an ordered `k`-tuple `x` of distinct indices let
`W_∅(x) = log N(gcd_{Z[i]}(z_{x_1},..,z_{x_k}))`. Then
`W_∅(x) = Σ_{p,t} log p · 1[all x_j on one side of layer (p,t)]`, and over uniform ordered tuples
`E[W_∅] ≥ ρ_{M,k} W`, `ρ_{M,k} = min_{1≤S≤M-1} (C(S,k)+C(M-S,k))/C(M,k) ≥ 2^{1-k}(1 - k(k-1)/(2(M-k)))`.
Consequently, if `M ≥ k + 2k(k-1)`, at least a fraction `2^{-k-1}` of ordered tuples have
`W_∅ ≥ 2^{-k} W`.
*Proof.* `gcd = ∏_p π^{min_j a}π̄^{e_p - max_j a}`, and `min a + (e - max a)` counts the layers on which all
`a_{x_j}` lie on one side. Normalisation gives `min_j a_j(p) = 0`, `max = e_p`, so `1 ≤ S_{p,t} ≤ M-1` for
every layer, and the probability that `k` distinct random rows lie on one side is the displayed ratio. With
`φ(x) = C(x,k)` for `x ≥ k-1` and `0` below (convex), `C(S,k)+C(M-S,k) ≥ 2φ(M/2)`, and
`φ(M/2)/C(M,k) = 2^{-k}∏_{j<k}(M-2j)/(M-j) ≥ 2^{-k}(1 - Σ_j j/(M-k))`. Since `0 ≤ W_∅ ≤ W`,
`P(W_∅ ≥ θW) ≥ (E W_∅/W - θ)/(1-θ)`; with `θ = 2^{-k}` and `k(k-1)/(M-k) ≤ 1/2` this is `≥ 2^{-k-1}`.
`□` (Exact check of `ρ_{M,k}` and its bound, and of `E[W_∅] ≥ ρW` on the record's 10-point cluster:
`gauss_checks/gcdred.py`.)

**Lemma 3.3 (counting; proved).** For a proper Zariski-closed `Z ⊂ T_k` there is `c(Z)` such that for every
circle and every set `Σ` of `M` of its lattice points, at most `c(Z)M^{k-1}` ordered `k`-tuples from `Σ`
have `[z] ∈ Z`.
*Proof.* Over `C`, `u = x+iy` identifies the circle with a subset of `C^*`, and `(z_j) ↦ [z]` is the
quotient map `(C^*)^k → (C^*)^k/C^*`, which is surjective; so the preimage `Z'` of `Z` is a proper closed
subset of `(C^*)^k ⊂ (P^1)^k`, independent of the circle. The bound `#(Σ'^k ∩ Z') ≤ c(Z')M^{k-1}` for
`Σ' ⊂ P^1`, `|Σ'| = M`, is the induction of `ptolemy.md` Lemma 6.1 (project to the first `k-1` factors;
over the proper closed set where the fibre is all of `P^1` use induction, elsewhere fibres have bounded
size). `□`

**Theorem 3.4 (reduction to five points on a subcritical arc; proved).** If `A_k(δ)` holds for some
`k ≥ 3` and some `0 < δ ≤ 2^{-k-1}`, then for every `C` there is `M(C)` bounding the number of lattice points
on any arc of length `C√R` of any circle `x^2+y^2 = R^2`.
*Proof.* Cut the arc into `s = ⌈C⌉` arcs of length `≤ √R` and keep one with `≥ M_orig/s` points; divide its
points by their Gaussian gcd `g_0`: the new radius `R_1 = R/|g_0|` and the points have pairwise distances
`≤ √R/|g_0| ≤ √R_1`. Call the number of points `M`. If `M ≥ M_k(Z) := ` the least integer `≥ k+2k(k-1)`
with `2^{-k-1}M(M-1)···(M-k+1) > c(Z)M^{k-1}`, Lemmas 3.2–3.3 give an ordered tuple with
`W_∅ ≥ 2^{-k}W` and `[z] ∉ Z`. By Lemma 1.2 (`θ = 2^{-k}`, `δ ≤ θ/2`) its reduction is a primitive tuple
on a circle of radius `R'' ≥ √R_1` with distances `≤ R''^{1/2-δ}` and the same class `[z] ∉ Z`. `A_k(δ)`
forces `R'' < R_0`, so `R_1 < R_0^2` and `M ≤ R_0 + 1` (an arc of length `√R_1` has `≤ √R_1 + 1` lattice
points). Hence `M_orig ≤ ⌈C⌉ max(M_k(Z), R_0 + 1)`. `□`

No balance, no residues, no unit class and no item-532 input is used: only the gcd averaging.
Using the refereed balanced selection of item 532 / `ptolemy.md` Step 1 (each of the `2^k` sign patterns
carries `(1±η)W/2^k`, `η = 2√2(2^k-1)/√(M-k+1)`, for at least half of the tuples) the same proof gives:

**Theorem 3.4′ (sharp exponent).** `A_k(δ)` for some `δ < 1/(2^k-2)` implies the uniform theorem; so does
`(L_k^circ)` (Steps 0–5 of `ptolemy.md` §6.2 with Step 2 run in `T_k` via Lemma 3.3 and Step 4 replaced by
`(L_k^circ)`: Step 1 supplies `η → 0` and `H_Σ ≤ k(k-1)W/(2(M-1)) = O(k^2 2^k w/M)`, `w = W/2^{k-1} → ∞`).
*Check of the implication (b).* Step 0 (subdivision, unit class, gcd) only shrinks radii and arcs. Step 1
holds verbatim for any `k` (item 532 eqs. (2)–(4), (10)). Step 2 needs the exceptional set to be avoided by
a positive proportion of tuples of **one** circle: Lemma 3.3 does this for `Z ⊂ T_k` because `[z]` is
similarity-invariant (it would **fail** for a non-invariant `Z ⊂ X_k`: the reductions `z/g_x` of different
tuples are rotated by different `g_x`, so they are not tuples of one circle in `X_k`'s coordinates; this is
why `A_k` and `(L_k^circ)` are stated on `T_k`). Step 4 is the lemma itself. Step 5: `H_Σ ≥ cw` against
`H_Σ ≤ O(k^2 2^k w/M)` bounds `M`. No step uses a fixed-`S` exceptional set, discards a gcd, or replaces
a phase by zero; `Z, c, η_0, w_0, R_0` are absolute.

### 3.3 The three candidates and their margins

| `k` | needed `A_k(δ)` | crude (Thm 3.4) | proved range (Cilleruelo–Córdoba) | `(K+D_{PQ}).γ̂` | balanced-count exponent per `w` |
|---|---|---|---|---|---|
| 3 | `δ < 1/6` | `δ ≤ 1/16` | `δ > 1/6` (three points need `(16R)^{1/3}`) | `-1/4` | `+1` |
| 4 | `δ < 1/14` | `δ ≤ 1/32` | `δ > 1/6` | `-1/8` | `+1` |
| **5** | **`δ < 1/30`** | `δ ≤ 1/64` | `δ > 1/10` (`R^{2/5}`) | **`+1/16`** | **`-1`** |
| 6 | `δ < 1/62` | `δ ≤ 1/128` | `δ > 1/10` | `+9/32` | `-9` |
| 7 | `δ < 1/126` | `δ ≤ 1/256` | `δ > 1/14` | `+33/64` | `-33` |

(`gauss_checks/kclass.out`.) The "proved range" column is the Cilleruelo–Córdoba theorem with `Z = ∅`:
`k` points need an arc `≥ √2 R^{1/2 - 1/(4⌊(k-1)/2⌋+2)}`. It reaches the needed range only at `k = 3`, and
only at its boundary (where constants decide, and the averaging loses `e^{O(W/M)}`; this is the
`log R/log log R` mechanism again). The last two columns are explained in §4.2: they change sign
exactly between `k = 4` and `k = 5`.

**Which `k` is "strongest".** `A_5(δ) ⇒ A_6(δ) ⇒ A_7(δ)` for the same `δ` (forget points; `Z` pulled back;
sub-tuples only get shorter after their own gcd reduction), but the needed `δ` shrinks with `k`, so
`A_5(<1/30)`, `A_6(<1/62)`, `A_7(<1/126)` are logically incomparable. The heuristic margin
`(k-1)(1+2δ_k)/4 - 1` (with `δ_k = 1/(2^k-2)`; this is `(K+E_e).β/W'` of §4.2) grows with `k` (`1/15 ≈ 0.067`, `0.29`, `0.52`), so `k = 7` has the most room and `k = 5` the
fewest points. In the balanced form, `(L_k^circ)` is weaker (easier) than `A_k` for each `k`
(Proposition 3.1). The **weakest single statement in this record that implies the uniform theorem and
involves at most five points** is `(L_5^circ)`; in plain words: *five primitive lattice points of a circle
of radius `R` whose Gaussian profile is balanced over all 15 cuts cannot all lie within
`R^{7/15 + c}` of each other, off a fixed closed set of similarity classes*.

## 4. (c) Can (L_k^circ) be proved or refuted with Gaussian structure?

**Summary.** `k = 4`: refuted (Theorem 4.1), so the reduction cannot work below five points. `k = 5, 6, 7`:
neither proved nor refuted. They are predicted by Vojta's conjecture in a form that is non-vacuous already
at `k = 5` (Theorem 4.2, conditional). The curve-family refutation mechanisms are excluded
(Theorem 4.3). Every unconditional tool that uses the Gaussian structure is a Liouville argument, and
those are tight **exactly** at the exponent `1/2 - 1/(2^k-2)` that is needed (§4.4).

### 4.1 `k = 4`: an explicit dense family; `(L_4^circ)` and `A_4(δ)`, `δ < 1/14`, are false

**Construction.** For a nonempty `T ⊆ {1,2,3}` let `H_T(s) = q_T(s - r_T)` with roots `r_T ∈ Q(i)` and
`q_T ∈ Z` clearing denominators, `z_0 = ∏_T H̄_T`, `z_x = ∏_{T∋x}H_T ∏_{T∌x}H̄_T`, `P_x = ∏_{T∋x}(s-r_T)`.
Then `z_x/z_0 = P_x/P̄_x`, so the four points lie on one circle and are within `O(|z|/|P_x|)` of each other
as soon as `Im P_x(s)` is a **constant**. For a monic quartic this means `e_1, e_2, e_3` of its roots are
real. Given three roots with `Im σ_1 ≠ 0`, the fourth root
`r = (Re σ_1 Im σ_1 - Im σ_2)/Im σ_1 - i Im σ_1` makes `e_1, e_2` real, and then `e_3` is real iff
(`gauss_checks/k4phi.py`, exact)

```text
Φ(r_1,r_2,r_3) = Σ_{a<b} v_a v_b (u_a-u_b)^2 + (v_1+v_2)(v_1+v_3)(v_2+v_3)(v_1+v_2+v_3) = 0,   r_a = u_a + i v_a.
```

So the family is: `r_123 = i` (normalised by `s ↦ λs+μ`), `r_12, r_13, r_23` with `Φ(i, ·, ·) = 0` on all
three pairs, singleton roots from the formula. Its moduli are 3-dimensional; a **rational 2-parameter
subfamily** is given by a common real part: `r_12 = u + i(x_a-1)`, `r_13 = u + i(x_b-1)`,
`r_23 = u + i(x_c-1)` with `x_a + x_b + x_c = 1`, `u^2 = x_a x_b x_c` (then `Φ = (v_a+v_b)(u^2 - x_ax_bx_c) = 0`),
parametrised by `x_a = X`, `x_b = X(1-X)/(m^2+X)`, `u = m x_b`, `(X, m) ∈ Q^2`.

**Exact checks.** (i) Instance `X = -3, m = 3`: `r_12 = -6-4i`, `r_13 = -6-3i`,
`r_23 = -6+5i`, `r_1 = -7+6i`, `r_2 = -3-2i`, `r_3 = -4-3i` gives `Im P_1 = -420`, `Im P_2 = -180`,
`Im P_3 = -240` identically; for integer `s` up to `10^30` the four lattice points have bounded common
Gaussian gcd (norm `≤ 10368000`), bounded residues (`|t_ij| ≤ 45`, `Y = -420,-180,-240`), balance
`η → 0`, and reduced diameter `/ R''^{3/7}` between `8.3` and `43.4` (`k4verify.out`). (ii) The map
`(X, m, s) ↦ [z] ∈ T_4` has Jacobian rank `3` at four rational points (exact forward-mode differentiation
over `Q`, `k4subfam.out`), so it is dominant; the full 4-dimensional family `(S, s) → P(X_4)` has rank 7
(submersive) at four rational points (`k4dominance.out`).

**Theorem 4.1 (proved, with the exact rank computation as input).** For every `δ < 1/14`, every proper
closed `Z ⊂ T_4` and every `R_0`, there are primitive 4-tuples of distinct lattice points of circles of
radius `≥ R_0`, with all distances `≤ R^{1/2-δ}`, balanced over the 7 cuts with `η → 0`, `H_Σ = O(1)`, and
`[ζ] ∉ Z`. Hence `A_4(δ)` is false for all `δ < 1/14` and `(L_4^circ)` is false for every `c, η_0, Z`.
*Proof.* `Ψ: (X,m,s) ↦ [z]` is dominant, so `Ψ^{-1}(Z)` is contained in a proper closed subset of `A^3`
and `B = {(X,m): {(X,m)} × A^1 ⊆ closure Ψ^{-1}(Z)}` is a proper closed subset of `A^2`. Also proper
closed is the locus where two of the seven roots coincide or are conjugate, or some `Im σ_1 = 0`. Take
`(X,m) ∈ Q^2` outside both. Then the `z_x(s) ∈ Q(i)[s]` have no common factor (a common root would be a
root of some `H̄_T`, hence a coincidence of roots), so `Σ_x a_x(s) z_x(s) = c` with `a_x ∈ Z[i][s]`,
`0 ≠ c ∈ Z[i]`, and the Gaussian gcd `g(s)` of the four lattice points divides `c`. The reduced tuple has
radius `R'' = |z_0(s)|/|g(s)| ≍ s^7` and diameter `≤ 2 max_x |Im P_x|·|z_0|/(|g||P_x|) ≍ s^3`, i.e.
`≤ C R''^{3/7} ≤ R''^{1/2-δ}` for large `s`; only finitely many `s` give `[z] ∈ Z`. The cut norms differ from
`n_T(s) = q_T^2|s-r_T|^2` by bounded factors (divisors of `c` and of resultants), so `η → 0`, and
`|t_ab| = |ζ_a-ζ_b|e^{d_ab/2}/(2R'') = O(s^3 s^4/s^7) = O(1)`. `□`

So the threshold `k = 5` is sharp: Theorem 3.4′ at `k = 4` needs exactly the statement that this family
refutes (the family sits at exponent `3/7 = 1/2 - 1/14`, the balanced reduced exponent). The dimension
count agrees: polynomial (degree-`d` block) families have real moduli dimension `1+2d` at `k=4` and
`2-2d` at `k=5`.

### 4.2 `k = 5`: what Vojta's conjecture says (conditional)

Let `Y_k = Bl_e L̄_k^tw` (smooth projective over `Q`, dimension `k-1`), `E` the exceptional divisor over
the identity `e`. `K_{Y_k} = -π^*Σ_ρ E_ρ + (k-2)E`. Vojta's Main Conjecture for `(Y_k, D = E, S = {∞})`
(Vojta, LNM 1239, Conj. 3.4.3) asserts: for every `ε > 0` there are a proper closed `Z_ε` and `C_ε` with
`λ_{E,∞}(y) + h_{K_{Y_k}}(y) ≤ ε h_A(y) + C_ε` for `y ∈ Y_k(Q) \ Z_ε` (`A` ample). Call this `V_E(k, ε)`.
`D = E` is smooth and `S` is the fixed set `{∞}`: no uniformity in an unbounded `S` is involved.

**Theorem 4.2 (conditional on `V_E(k,ε)`).** Let `k ≥ 5`, `δ > 0` with `κ := (k-1)(1+2δ)/4 - 1 > 0`, and
`ε < κ/c_A`. Then `A_k(δ)` holds with `Z = Z_ε ∩ T_k`. In particular `V_E(5, ε)` for one
`ε < (1/32)/c_A` implies the uniform theorem (Theorem 3.4 with `δ = 1/64`).
*Proof.* Let `ζ` be primitive on `|ζ| = R''` with diameter `≤ R''^{1/2-δ}`, `y = [ζ] ∈ T_k(Q)`,
`W' = log R''^2`. By Lemma 1.1, `h_{π^*E_ρ}(y)` are the toric heights and `Σ_ρ h_{E_ρ} = W'`. The
archimedean Weil function of `E` near `e` is `-log max_j |ζ_j/ζ_1 - 1| + O(1) ≥ (1/2+δ) log R'' - O(1)`,
and `h_E ≥ λ_{E,∞} - O(1)` (finite local heights of the effective `E` are `≥ 0`). Hence
`λ_{E,∞} + h_K ≥ (k-1)(1/2+δ)W'/2 - W' - O(1) = κW' - O(1)`. Take `A = π^*L - ε'E` with `L = Σ a_ρE_ρ`
ample toric, `a_ρ > 0`; then `h_A ≤ (max a_ρ) W' + O(1) =: c_A W' + O(1)`. Off `Z_ε`,
`κW' ≤ εc_AW' + O(1)`, so `W'` is bounded. `□`

*Remarks.* (1) With `D = 0` the same computation needs `(k-2)(1+2δ)/4 > 1`, i.e. `k ≥ 6` — this is the
pointer of `threshold/L8.md` §6.1 (`K.γ̂ > 0` first at `k = 6`). With the record's `X_k` and no gcd
reduction, `D = 0` needs `k ≥ 7` (`growth/conditional.md`). The two gains are: the gcd reduction (§1.3)
and the archimedean divisor `D = E` at the fixed place `∞`. (2) This is a third, logically independent
Vojta instance implying the uniform theorem (4-dimensional, against 7 for `X_7` and 5 for `M̄_{0,8}`); its
value here is that it **predicts `(L_5^circ)`**, i.e. the five-point statement is not expected to be false.
(3) `(K + D_{PQ}).γ̂` on `M̄^tw_{0,k+2}` equals `-1/4, -1/8, +1/16, +9/32, +33/64` for `k = 3..7`
(`kclass.out`); the balanced lattice count (blocks of norm `≈ e^w`, `k-1` angle conditions of precision
`e^{-2^{k-3}w}`) is `exp((2^{k-3}(5-k) - 1)w) = exp(-(K+D_{PQ}).γ̂ · 2^{k-1} w)` (**heuristic**):
`+1, +1, -1, -9, -33`.

### 4.3 The refutation mechanisms fail at `k ≥ 5`

Call `y ∈ Y_k(Q)` an *A-point* (for `δ > 0`) if it is the class of a primitive tuple with diameter
`≤ R''^{1/2-δ}`; then `λ_{E,∞}(y) ≥ (1/4+δ/2) h_B(y) - O(1)` with `B = π^*Σ_ρE_ρ` (`h_B = W'`). Every
`(L_k^circ)`-counterexample with small `η_0, c` is an A-point (proof of Prop. 3.1).

**Theorem 4.3 (proved, modulo the standard facts cited).** Let `k ≥ 5`, `δ > 0`.
(i) An irreducible curve `C ⊂ Y_k` over `Q` with infinitely many A-points has geometric genus `0`; its
normalisation `ν: P^1 → Y_k` has, for the point `O` approached, either `O ∈ P^1(Q)` with
`μ_0 := mult_O ν^*E ≥ (1/4+δ/2)(B.C)/(1+ε)`, or `O` real quadratic with conjugate `O'`,
`mult_O = mult_{O'} = μ_0 ≥ (1/4+δ/2)(B.C)/(2+ε)`, for every `ε > 0` and infinitely many of the points.
(ii) For every `d`, the union of the curves `C` with `B.C ≤ d` carrying infinitely many A-points is
contained in a proper closed subset of `Y_k`. In particular no bounded-degree covering family of curves
whose members carry infinitely many A-points exists.
*Proof.* (i) If `B.C = 0` all toric heights are bounded on `C` and there are finitely many tuples. Otherwise
`h_B|_C` is a height of positive degree. Genus `≥ 2`: `C(Q)` is finite (Faltings). Genus `1`: for
`O ∈ C ∩ E` and rational points `P → ∞`, `λ_{O,∞}(P) = o(h(P))` (Siegel's theorem in the form of
Silverman, *AEC*, Thm IX.3.1), so `λ_{E,∞} = o(h_B)` along `C(Q)`, contradicting
`λ_{E,∞} ≥ (1/4+δ/2)h_B - O(1)`. Genus `0`: `λ_{E,∞}(ν(t)) = Σ_O mult_O λ_{O,∞}(t) + O(1)`; for real `t`
approaching one real point `O`, the other terms are bounded; Liouville (`O` rational, `λ_O ≤ h + O(1)`) or
Roth (`O` irrational algebraic, `λ_O ≤ (2+ε)h + C_ε`) with `h_B(ν(t)) = (B.C)h(t) + O(1)` gives the bounds;
non-real `O` cannot be approached, and an irrational `O` of degree `≥ 3` is excluded by the computation in
(ii) (it only strengthens the contradiction).
(ii) Curves of bounded degree form a finite-type family. If the union in question were dense, some
irreducible component `F` of the space of maps `P^1 → Y_k` of bounded degree, with a dense set of members
of the type in (i), would dominate `Y_k`; the condition "`ν^*E` has multiplicity `≥ μ_0` at one (resp. two)
points" is closed, so the general member of `F` has it; reparametrise so that the point(s) are `∞`
(resp. `0, ∞`). The tangent space to maps with this fixed contact is contained in
`H^0(P^1, F_ν)`, `F_ν ⊂ ν^*T_{Y_k}` the subsheaf of fields `v` with `dh(v)` vanishing to order `μ_0` at the
fixed point(s) (`h` a local equation of `E`), a subsheaf of colength `μ_0` per point. Dominance and generic
smoothness make `F_ν` generically globally generated, hence nef: `deg F_ν ≥ 0`, i.e.
`(#fixed points)·μ_0 ≤ -K_{Y_k}.C = B.C - (k-2)E.C`. With `E.C ≥ (#orbit)μ_0`:
rational `O`: `(k-1)μ_0 ≤ B.C`, against `μ_0 ≥ (1/4+δ/2)B.C/(1+ε)`, impossible for `k ≥ 5`, `ε < 2δ`;
quadratic: `2(k-1)μ_0 ≤ B.C` against `μ_0 ≥ (1/4+δ/2)B.C/(2+ε)`, impossible for `k ≥ 5`, `ε < 4δ`;
degree `d ≥ 3` orbit (fix two points): `(2+(k-2)d)μ_0 ≤ B.C` against the Roth bound, impossible. `□`

At `k = 4` the rational case gives only `3(1/4+δ/2) ≤ 1`, i.e. `δ ≤ 1/6`, and Theorem 4.1's family
(`δ = 1/14`, `μ_0 = 4`, `B.C = 14`, `3·4 ≤ 14`) is exactly such a covering family. What Theorem 4.3 does
**not** exclude, as in round 18 (correction M1): infinitely many non-covering rational curves of unbounded
degree, covering families of higher-genus curves with finitely many A-points per member, and sporadic
points. A refutation of `(L_5^circ)` would refute `V_E(5,ε)`.

### 4.4 Why the Gaussian structure does not prove `(L_5^circ)`: every such tool is Liouville, and Liouville is tight

* **Pair bound = residue integrality.** `|ζ_a - ζ_b| = 2|t_ab|R''e^{-d_ab/2}` with `|t_ab| ≥ 1` is the only
  inequality supplied by one pair, and for a balanced tuple it gives diameter `≥ 2R''^{7/15}` at `k = 5`:
  **exactly** the needed exponent, with nothing to spare. Averaging it over pairs is Cilleruelo–Córdoba:
  `Σ_{a<b} d_ab = Σ_A |A||A^c| log n_A ≤ 6W'` forces `δ ≤ 1/10`; balanced profiles force `δ ≤ 1/30`.
* **Every auxiliary polynomial with forced Gaussian divisibility** (`threshold/construct.md` Theorems B,
  C2, C3, refereed): for a Laurent polynomial `f` on `T_k` with Newton polytope `Q`, the forced-divisibility
  Liouville inequality excludes balanced bounded-residue points iff `2^{k-2} ord_e f > Σ_cuts width_Q(1_A)`,
  and `ord_e f ≤ 2^{1-k}Σ_{A⊆[k]} width_Q(1_A)` (construct.md Thm B, "mean width"). In the language of
  §1: `Z_f.γ̂ ≥ 0` for the divisor `Z_f` of every `f`, i.e. **`γ̂` is movable on the twisted `M̄_{0,k+2}` for
  every `k`** (the statement `L8.md` §6.1 verified only at `k = 4`). So no single-place auxiliary-polynomial
  argument on the full Gaussian model proves `(L_k^circ)` for any `k`, and products of pair differences
  (Vandermonde, Cilleruelo–Córdoba) reach ratio exactly `1`.
* **Roth/Schmidt/Ru–Vojta at `S = {∞}`.** The target `e` is a **rational** point and `E` is geometrically
  irreducible; Ru–Vojta's bound for one irreducible divisor is `m_{E,S} ≤ (β(L,E)^{-1}+ε)h_L` with
  `β(L,E) = ∫_0^∞ vol(L-tE)dt/vol(L) ≤ τ(L,E) := sup{t : L - tE effective}`, never better than Liouville.
  Roth-type gains need non-rational (Galois-conjugate) components in the target; there are none.
* **`S`-unit / Subspace with `S` = primes of the tuple:** checklist item 5 (thresholds through `|S|` only).
* **Norm equations, units, nested-chain identities (items 515–533).** The Gaussian Ptolemy identity
  (`ptolemy.md` Thm A(3)), the anchored identity (Thm B(i)), item 532's Hankel divisor (13), and the
  norm relations `N(P_i) = ∏_{T∋i} n_T` are polynomial identities with forced divisibility; used as
  auxiliary polynomials they fall under the previous bullet. Units only change bounded factors. Isotropy
  is the toric structure itself (Lemma 1.1). I found no Gaussian identity that is not of this form.

**Status of `(L_5^circ)`, `(L_6^circ)`, `(L_7^circ)`: open.** Predicted by `V_E(k)` (Theorem 4.2);
not refutable by density or covering curve families (Theorem 4.3, §4.2 count); not provable by any
single-place auxiliary polynomial (construct.md) or by Roth/Schmidt at `∞` (rational target). This is the
same position as `(L_8)` in round 18, three points earlier.

## 5. Barrier checklist, item by item

The "proposal" of this round is the reduction **uniform theorem ⟸ `(L_5^circ)` (or `A_5(δ)`, `δ < 1/30`)**
together with the claim that this is the smallest fixed-size circle statement that can work. The
reduction is proved (Theorems 3.4, 3.4′); the question is whether the target statement is capped.

1. *No fake-satisfied input.* `(L_5^circ)` is not an input but the target. It is **violated by fakes**:
   a balanced angle model on five rows with angles placed at the pair-gap boundary satisfies every
   monomial gap (Lemma 5.1 of the overview is exactly `|t_ab| ≥ 1`) and has reduced diameter
   `≍ R''^{7/15}`. So any proof must use an input fakes violate; the reduction itself uses only the
   gcd identity of Lemma 3.2, which fakes also satisfy, and therefore proves nothing on its own. Not capped,
   not a mechanism.
2. *No robust/averaged analytic statement.* Not used; `(L_5^circ)` is pointwise for one 5-tuple.
3. *No single-place auxiliary polynomial on `X_M`, none on `M̄_{0,k}`.* The circle analogue is **also capped**:
   `γ̂` is movable on the twisted `M̄_{0,k+2}` for every `k` (§4.4, from `construct.md` Thm B, refereed), so no
   Liouville argument on the full Gaussian model proves `(L_k^circ)`. Ru–Vojta at `S = {∞}` is no better
   than Liouville here (one irreducible divisor, rational target).
4. *No single linear form in logarithms; SLW false; `SLW_w` Vojta-strength.* Not used. The only conditional
   input discussed is `V_E(5)` (Theorem 4.2), which is Vojta strength and is offered as a prediction, not a
   route.
5. *No `S`-integer argument with `|S|`-thresholds.* Not used; Theorem 4.2 has `S = {∞}` fixed and `D = E`.
6. *Route 5 (`R_H`, CIRCLE(q)).* Not used.
7. *No `(L_5)/(L_6)`-type family mechanism at `k = 8`; `(L_8)` untestable.* The analogue for circles at
   `k = 5` is proved here (Theorem 4.3): genus `≥ 1` curves carry finitely many A-points, bounded-degree
   covering rational families are impossible. A census is uninformative for the same reason as round 18's Thm E: balanced 5-tuples need 15
   comparable blocks (several distinct split primes, since the blocks at one prime form two nested chains),
   and the predicted count is `e^{-w}` per scale. What exists as evidence:
   for `N ≤ 10^{12}` no five lattice points lie on an arc shorter than `1.9962√R` (`growth/families.md`,
   reproduced to `10^{10}` by its referee), so `A_5(δ)` holds with `Z = ∅` in that range for every `δ > 0`
   (**evidence**, small radii, says nothing asymptotically).
8. *abc unit-pencil lane closed.* Not used.
* *Dirichlet-exponent fact (`round4/hardness.md`).* The reduction moves the problem from the critical
  exponent `1/2` to `7/15`, where the balanced count is subcritical (`e^{-w}`). But the gain is exactly
  eaten by Liouville's tightness: at the balanced profile the pair bound **is** the exponent `7/15`, so the
  needed statement sits at the Liouville exponent of the restricted class, and since the target point `e`
  is rational, Roth/Schmidt-type methods (which beat Liouville only for irrational targets) do not apply.
  This is the same obstruction in a sharper form, not a way around it.

**Constraints of the task.** No uniform exceptional set from a fixed-`S` Subspace theorem is assumed (the
only exceptional sets are a hypothesis `Z` in `A_k`/`(L_k^circ)` or Vojta's for a fixed variety with
`S = {∞}`). No Gaussian gcd or denominator is discarded (the tuple gcd is the whole point of §1.3;
residues are exact integers). No drifting phase is replaced by zero (no phases appear; the family of §4.1 is
exact). No template-dependent threshold is used as uniform (Theorem 4.1 is a refutation and may depend on
the template; Theorems 3.4, 3.4′, 4.2 have absolute constants).

## 6. Computations (scripts in gauss_checks/)

All under `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round20/gauss_checks/`;
outputs in the `.out` files. Exact unless marked.

| script | content | result |
|---|---|---|
| `kclass.py` | `K.γ̂`, `(K+D_{PQ}).γ̂` on `M̄_{0,k+2}` from `K = Σ(s(n-s)/(n-1)-2)D_S`; `X_k`/`Y_k` margins; direct block count | `K.γ̂` matches `L8.md` §6.1 (`-3/8,-3/16,1/32,17/64,...`); `(K+D_{PQ}).γ̂ = -1/4,-1/8,+1/16,+9/32,+33/64` (`k=3..7`); count exponent `= -(K+D_{PQ}).γ̂·2^{k-1}` for `k = 3..9` |
| `gcdred.py` | Lemma 3.2: `ρ_{M,k}` and its lower bound; on the record's 10-point cluster (`n = 1176852625`) all 5- and 3-subsets: `R'' ≥ √R`, `E[W_∅] ≥ ρW`, reduced exponents | bound valid at all tested `(M,k)`; `R'' ≥ √R` always; `E[W_∅/W] = 0.0235 ≥ ρ = 0.0079` (`k=5`); reduced exponents `0.61–0.70` (**evidence**: small `R`, `C ≈ 8`; irrelevant asymptotically) |
| `k4phi.py` | symbolic derivation of `Φ` | `Φ = Σ v_av_b(u_a-u_b)^2 + (v_1+v_2)(v_1+v_3)(v_2+v_3)(v_1+v_2+v_3)` |
| `k4search.py`, `k4search2.py` | rational solutions with small height | 30 nondegenerate rational triples at height `≤ 3`; the common-real-part pattern found and then proved (§4.1) |
| `k4verify.py` | lattice points of the instance for `s` up to `10^{30}` | `Im P = -420,-180,-240`; bounded gcd and residues; diameter/`R''^{3/7}` in `[8.3, 43.4]` |
| `k4dominance.py`, `k4subfam.py` | exact Jacobian ranks by forward-mode differentiation over `Q` | rank 3 onto `T_4` for the rational subfamily (4 points); rank 7 onto `P(X_4)` for the full family (4 points) |
| `familycheck.py` | isotropy/norm test of C1 and `E_t` with isotropic points `±i` | resultants nonzero (max `4.1·10^{10}`); `φ = 0` on C1's class; `φ ≤ 0.027` on `E_t` (`N = 2..14`) |


## 7. Open problems / what is left

1. **`(L_5^circ)`** (equivalently, for the uniform theorem, `A_5(δ)` for one `δ < 1/30`): five primitive
   lattice points of a circle of radius `R`, balanced over all 15 Gaussian cuts, cannot all lie within
   `R^{7/15+c}` of each other, off a fixed proper closed set of similarity classes. Open. Predicted by
   `V_E(5)`. Its unrestricted form `A_5(δ)` is a Cilleruelo–Córdoba-type statement (five points on an arc
   `R^{1/2-δ}`) needed for `δ` just below `1/30`, against the proved `δ > 1/10`.
2. **Intermediate exponents.** Prove `A_5(δ)` for some `δ ∈ (1/30, 1/10]` (off a closed set). This would be
   the first improvement of Cilleruelo–Córdoba at five points since 1992, but by Theorem 3.4′ it implies
   nothing about the critical scale. Pair-bound averaging shows the configurations to be excluded there are
   those with singleton-cut weight `W'_1 ≤ (1-10δ)W'/2`; Liouville may or may not reach them (open; the
   cube-profile barrier only covers `δ ≤ 1/30`).
3. **Rotation-sensitive exceptional sets.** `A_4(δ)` with `Z ⊂ X_4` not similarity-invariant is not refuted
   by Theorem 4.1 (the primitive tuples of a family form a 3-dimensional subset of the 4-dimensional `X_4`);
   the full 4-dimensional family is submersive onto `P(X_4)` at rational points, so it would be refuted if
   the rational points of its 3-dimensional parameter variety are Zariski dense. This version does not feed
   Theorem 3.4 (Lemma 3.3 needs invariant `Z`), so it is of calibration interest only.
4. **For experts.** Is there any unconditional method for approximating a *rational* point of a variety at
   a single place with heights restricted to a fixed ray (here: toric heights balanced on the permutohedral
   rays, archimedean contact with the exceptional divisor over `e`) that beats Liouville? Every such
   method would apply to `(L_5^circ)`. The movability of `γ̂` (§4.4) says no auxiliary section does.

