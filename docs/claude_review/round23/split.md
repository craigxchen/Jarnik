# Round 23: the split calibration of Question B (split Losev-Manin LM_5 at e)

## 0. Verdict

**No growth improvement and no uniform bound.** `M ≤ (2/3 + ε) log R/loglog R` stands. The split calibration is
**not easier** than Question B. It is a dead end as a source of a mechanism and useful as a calibration. Status: `partial`.

1. **Precise split question (§2).** In McKinnon–Roth normalisation, on the split `LM_5` at `e` with `-K` (polytope
   `R_5`), it asks whether `α^ess_e(-K) > 15/4`. The window is the same, `(15/4, 4)`, and so is the section cap
   `τ_e ≤ 15/4`.
2. **Exact elementary form (§3).** The height of a primitive positive 5-tuple is `H = (max u/min u)·lcm(u)`, and
   `dist ≍ (max u - min u)/min u`. The split `A_5(δ)` says: five distinct integers with gcd 1, lying in
   `[U, U(1 + lcm^{-(1+2δ)/4})]`, have ratio class in a fixed proper closed `Z` once their lcm is large. Equivalently,
   five factorisations `N = d_j·(N/d_j)` within **hyperbolic angle** `N^{-(1+2δ)/4}` on `xy = N`. At `δ = 0` this is
   the Lorentzian Cilleruelo–Granville scale `radius^{1/2}`, and at the vertex it is the Erdős–Rosenfeld window
   `√n + C n^{1/4}`.
3. **Theorem 1 (proved).** `A^split_5(δ)` for one `δ < 1/30` implies Erdős #887 in its absolute-`K` form, uniformly in
   the position on the hyperbola. It also gives Ruzsa's #886 for `ε > (15/64)(1+2δ)`, which is below `1/4`. So the split
   statement is at least as hard as these **open** Erdős problems, and it is **not known** (a). (Canonical
   boundedness of all split `LM_k` would give #886 in full.)
4. **Theorem 2 (proved; exact and numerical checks).** The July 2026 five-divisor #887 construction (m-czech, via
   formal-conjectures PR #7218) is, in torus language, a **cut-balanced rational curve at exactly the Liouville constant
   `15/4`**: `B·C = 15`, each of the 15 cuts hit once, and two real-quadratic branches of order 2. It belongs to an
   explicit 2-parameter family `C_{abc}` (`d = a+b+c`, conic `z² = abcd·x² + B`) that sweeps a **hypersurface
   `W ⊂ T`**. So every admissible `Z` in the split `A_5` / `(L_5)` must contain `W`. The twisted form provably cannot
   carry such curves: Galois forces even weight on every cut (parity lemma 1.4). The record's Gaussian Pell and golden
   families sit at `α = 4`. Whether its cyclotomic five-row family (`C → 0`) has `α < 4` was not checked.
5. **True or false (b): undecided.** One-parameter subgroups (`α ≥ 8`), coordinate products, covering families
   (`α ≥ 4`, one or two branches) and polynomial blocks do not refute it. Refutation would need infinitely many
   Pell-type hypersurfaces `W_i`, a K3-type density question I could not decide. Vojta predicts that the split statement
   is true.
6. **Transfer (c).** There is no unconditional split mechanism. What transfers is the parity lemma, a cleaner reduction
   proof, and a new twisted search ansatz (two real-quadratic branches, `μ_0 = 4`), which has map value only. Logical
   remark (6.3): a proof of Question B that works for every form of `LM_5` would solve #887 (absolute `K`).

## 1. Setting: split LM_5, the identity, -K, heights and distance

Labels: **proved** (proof here), **exact** (exact computation in `split_checks/`), **record** (proved earlier in the
record and killer-checked, cited), **snippet** (literature seen only through search snippets or a GitHub page),
**evidence**, **heuristic**.

**1.1 The variety.** `N = Z^5/Z·1`, `M = {m ∈ Z^5 : Σ m_j = 0}`. The permutohedral fan in `N_R` (cones = chains of
proper nonempty subsets, rays `1_A`, `A ⊊ [5]`, 30 rays) defines the **split Losev–Manin 4-fold** `LM_5` over `Q`. Its
open torus is `T = G_m^5/G_m`; `T(Q) = (Q^*)^5/Q^*`. Every class has a unique representative
`u = (u_1,…,u_5) ∈ Z^5`, `gcd(u) = 1`, up to the sign `±u`; near `e` all `u_j` have one sign, so we take `u_j > 0`.
`e = [1:1:1:1:1]`. The twisted form `Lbar_5^tw` of the record is the same fan with Galois acting by `-1` on `N`;
its torus is the anisotropic `T_5`. All intersection-theoretic data are the same over `Qbar` (record, round 22):
`-K = Σ_A E_A`, nef and big, not ample; polytope `P_{-K} = R_5 = conv(e_a - e_b)`; `(-K)^4 = 70`;
`ε_e(-K) = 2`; `10/3 ≤ τ_e(-K) ≤ 15/4`; balanced Liouville constant `15/4`; covering rational families through `e`
have `(B·C)/μ_0 ≥ 4`. Only the **arithmetic** (which points and curves are rational) differs.

**1.2 Height (proved; judge J2 checked it exactly on 300 tuples).** For primitive positive `u`,

    H(u) := H_{-K}([u]) = Π_v max_{a,b} |u_a/u_b|_v = (max u / min u) · lcm(u).

*Proof.* The toric (polytope) height of `R_5` has local factor `max_{m ∈ R_5} |χ^m|_v = max_{a,b}|u_a/u_b|_v`. At `p`
it is `p^{max_j v_p(u_j) - min_j v_p(u_j)}`; the product over `p` is `lcm(u)/gcd(u) = lcm(u)`. At `∞` it is
`max u/min u`. `□` Near `e`, `H(u) = lcm(u)(1 + o(1))`. We write `L(u) = lcm(u)`.

**1.3 Distance.** In the affine chart `y_j = u_j/u_1 - 1` (`j = 2..5`), with the sup norm, put
`dist([u], e) = max_j |u_j/u_1 - 1|`. With `D = max u - min u` and `U = min u`,
`D/(U + D) ≤ dist ≤ D/U`, so `dist ≍ D/U` near `e`. Any other choice of local coordinates changes `dist` by a bounded
factor, which never matters below.

**1.4 Cut heights.** For a prime `p` and a threshold `s ≥ 1` put `A_{p,s} = {j : v_p(u_j) ≥ s}`. When
`∅ ≠ A_{p,s} ≠ [5]` the pair `(p,s)` is a *layer* of type `A_{p,s}`. Put `n_A = Π_{layers of type A} p` (30 numbers, one
per ray) and `N_{A|A^c} = n_A n_{A^c}` (15 numbers, one per cut). Then (proved, same bookkeeping as 1.2)

    lcm(u) = Π_A n_A,  gcd(u_a,u_b) = Π_{A ∋ a,b} n_A,  u_a = Π_{A ∋ a} n_A,  h_{E_A}([u]) = log n_A + O(1).

**Difference from the twisted form (proved).** In the split form `n_A` and `n_{A^c}` are independent. In the twisted
form, complex conjugation sends a Gaussian layer of type `A` to one of type `A^c`, so the record's cut norm `n_A`
(`gauss.md` Lemma 1.1) is attached to the *cut*, and for a curve `C` over `Q`, `E_A·C = E_{A^c}·C`. Indeed, a place of
the curve fixed by conjugation has valuation vector `v` with `v = -v` in `N`, so `v = 0`. Places with `v ≠ 0` come in
conjugate pairs with vectors `1_A`, `1_{A^c}`. Consequence: **a twisted curve over `Q` has even weight on every cut;
a split curve can have weight 1 on a cut.** This is the one geometric-arithmetic asymmetry that matters in §5.

## 2. The split analogue of Question B, precisely

Normalisation as in round 22: an **A-point at exponent `γ`** is `y ∈ T(Q)` with `dist(y,e)^γ · H(y) ≤ 1`; the record's
`δ` corresponds to `γ = 4/(1+2δ)` (`15/4 ↔ δ = 1/30`, `10/3 ↔ δ = 1/10`, `4 ↔ δ = 0`).

* **Essential constant.** For a set `Ω ⊂ T(Q)`, `α^ess_e(-K; Ω)` is the supremum of `γ` for which some proper
  Zariski-closed `Z ⊂ T` leaves only finitely many A-points at exponent `γ` in `Ω \ Z`. (For a sequence tending to `e`,
  this matches McKinnon–Roth's `α`, as in round 22 toric Lemma 1.3; only the "finitely many outside `Z`" form is used.)
* **`A^split_5(δ)`.** There are a proper closed `Z ⊂ T` and `H_0` such that every primitive 5-tuple of distinct
  positive integers with `H(u) ≥ H_0` and `dist([u],e) ≤ H(u)^{-(1+2δ)/4}` has `[u] ∈ Z`.
  For `δ ∈ (0, 1/30)`, `A^split_5(δ)` is implied by `α^ess_e > 4/(1+2δ)` and implies `α^ess_e ≥ 4/(1+2δ)`. So
  "**`A^split_5(δ)` for some `δ < 1/30`**" is the split form of `α^ess_e(-K) > 15/4`.
* **`(L^split_5)` (balanced form).** The same, restricted to points whose 15 cut heights satisfy
  `|log N_{A|A^c} - w| ≤ η_0 w`. The split balanced Liouville constant is `15/4` (proved): the pair bound
  `|u_a - u_b| ≥ gcd(u_a,u_b)` gives `dist ≥ 1/Π_{A ∋ a, A ∌ b} n_A` for each ordered pair; the sum over the 20 ordered
  pairs of the logarithms is `Σ_A |A||A^c| log n_A = Σ_cuts |A||A^c| log N_{A|A^c} = 80w + O(η_0 w)`, so the best pair has
  `-log dist ≤ 4w(1+O(η_0))`, while `log H = 15w(1+O(η_0))`. Hence `α ≥ 15/4 - O(η_0)` on balanced points, and §5
  shows `15/4` is attained.
* **Strongest form.** `CB^split_e`: `LM_5` (or its anticanonical model `X'_split`) is canonically bounded at `e` in
  McKinnon–Satriano's sense, i.e. `α^ess_e(-K) ≥ 4`. It implies `A^split_5(δ)` for every `δ > 0` (round 22 §2 Step 1,
  verbatim: the argument uses only Northcott and that A-points avoid `Z_diag`).

> **Question B^split.** Is `α^ess_e(-K_{LM_5}) > 15/4` on `T(Q)` (the split torus), at least on cut-balanced
> sequences? Equivalently: is `A^split_5(δ)` true for some `δ < 1/30`? Vojta's conjecture for `(Bl_e LM_5, E, {∞})`
> predicts `≥ 4` (gauss Thm 4.2 / toric Prop 5.1; the derivation is geometric plus the height dictionary, which holds
> verbatim for the split form).

The window is the same as in the twisted case: `(15/4, 4)`. The section cap (`τ_e ≤ 15/4`) is the same, because the
cap is a statement about Laurent polynomials and widths and holds over any field (round 22, Prop 1.3 of lit).

## 3. Elementary translation: five integers with small lcm in a short interval

**3.1 Tuple form (exact restatement of §2, proved by 1.2-1.3).** Let `u_1 < … < u_5` be positive integers with
`gcd = 1`, `U = u_1`, `D = u_5 - u_1`, `L = lcm(u)`. Then `[u]` is an A-point at exponent `γ` iff, up to the factor
`(1 + D/U)` from 1.3,

    (D/U)^γ · (u_5/u_1) · L ≤ 1,    i.e.    D ≤ U · L^{-1/γ} (1 + o(1)).

So `A^split_5(δ)`, `δ < 1/30`, says: **five distinct integers with gcd 1 lying in an interval `[U, U + U·L^{-(1+2δ)/4}]`,
where `L` is their lcm, have their ratio class in a fixed proper closed `Z`, once `L` is large.** Liouville, the
balanced bound and the counting value read `L^{-3/10}`, `L^{-4/15}` and `L^{-1/4}` respectively in place of
`L^{-(1+2δ)/4}`.

**3.2 Divisor (hyperbola) form (proved).** Primitive `u` with `lcm(u) = N` is the same thing as five divisors `d_j = u_j`
of `N` such that `gcd(d_j) = 1` and `gcd(N/d_j) = 1` (the second condition is `lcm(d) = N`). The points `(d_j, N/d_j)`
lie on the hyperbola `xy = N`. The boosts `(x,y) ↦ (λx, y/λ)` preserve it, and they are the split analogue of rotations.
The invariant (Lorentzian) arc between `(d, N/d)` and `(d', N/d')` is `√N·|log(d'/d)|`; it is the arc length in the
Minkowski metric `dx dy`. So:

* `dist ≍ D/U` is the hyperbolic angle;
* an A-point at `γ = 4` is a set of five factorisations within hyperbolic angle `N^{-1/4}`, i.e. within invariant arc
  length `N^{1/4} = (√N)^{1/2}`. **This is exactly the Cilleruelo–Granville scale `C·radius^{1/2}`**, with radius `√N`.
* At the vertex `X = √N` the window `[X, X(1 + C N^{-1/4})]` is `[√N, √N + C N^{1/4}]`. That is the Erdős–Rosenfeld window of
  Erdős problems #886/#887 (§4).
* At `X = N^θ`, `θ ≠ 1/2`, the invariant window `X·C N^{-1/4}` differs from the Euclidean-curvature window `X^{1/2}`
  that lit18 `direct.md` §2.2 uses. The torus question is the **Lorentzian** one.

**3.3 The phase.** The position `log(X/√N)` along the hyperbola is the split analogue of the drifting rotation phase.
At phase `0` (the vertex) there is an archimedean integrality: `s = d + N/d` is an integer and
`s - 2√N ≍ (d - √N)^2/√N`. So divisors in `[√N, √N + C N^{1/4}]` give distinct integers `s` in an interval of length about
`C^2`. This is Erdős–Rosenfeld's bound `1 + C^2`. At `X^2/N = p/q` (lowest terms) the same trick with `s = qd + p(N/d)`
gives `O(1 + C^2 √(pq))`. For a general phase, Dirichlet only gives `pq ≲ N^{1/4}` and the bound is useless. This is the
exact parallel of lit18 §2.3 for circles (near-axis versus general direction). The phase is quotiented out in `T`, so
`A^split_5` sees no phase, exactly as `A_5` does in `T_5` (owner constraint: no phase is set to zero).

**3.4 Theorem 1 (proved). `A^split_5(δ)` for one `δ ∈ (0, 1/30)` implies an absolute, phase-uniform Erdős #887.**
Assume `A^split_5(δ)` with `Z`, `H_0`. There is `K = K(δ, Z)` such that for every `C > 0` there is `n_0(C)` with

    #{d | n : X ≤ d ≤ X(1 + C n^{-1/4})} ≤ K    for all n ≥ n_0(C) and all X > 0.

In particular (`X = √n`) Erdős #887 holds in its absolute-`K` form (formal-conjectures `erdos_887.parts.ii`). Also,
Ruzsa's #886 holds for every `ε > (15/64)(1 + 2δ)`: for those `ε`, `n` has `O_{δ,Z}(1)` divisors in
`(√n, √n + n^{1/2-ε})`. Since `(15/64)(1+2δ) < 1/4`, this reaches into the open range `ε ≤ 1/4`.

*Proof.* Let `d_1 < … < d_M` be the divisors of `n` in the window, with `M ≥ 5`. Choose a 5-subset `S` uniformly at
random. Let `g_S = gcd(d_S)`, `u_S = d_S/g_S` and `L_S = lcm(u_S) = lcm(d_S)/g_S`; `u_S` is primitive.

*(i) Distance.* `dist([u_S], e) ≤ (max d_S - min d_S)/min d_S ≤ C n^{-1/4}`.

*(ii) Averaging the height.* For `p^t ‖ n` and `1 ≤ s ≤ t`, let `q_s` be the fraction of the `d_i` with `v_p(d_i) ≥ s`.
Then `v_p(L_S) = max_S v_p - min_S v_p = Σ_{s=1}^t 1[S meets both {v_p ≥ s} and {v_p < s}]`. With `κ = q_s M`,
the probability that `S` lies on one side is `p(M,κ) = (C(κ,5) + C(M-κ,5))/C(M,5)`. Since `κ ↦ C(κ,5)` has second
difference `C(κ,3) ≥ 0` on integers and `C(x,5)` increases for `x ≥ 4`, for `M ≥ 9` we get
`p(M,κ) ≥ 2C((M-1)/2,5)/C(M,5) = 2^{-4} Π_{i=0}^{4} (M-1-2i)/(M-i) ≥ 2^{-4}(1 - 15/(M-4))`
(exact check `hypergeom.py`: `min_κ p ≥ 1/16 - 2/M` for `5 ≤ M ≤ 3000`). Summing over `p, s`:

    E[log L_S] ≤ (15/16 + 15/(16(M-4))) log n,   and   0 ≤ log L_S ≤ log n.

*(iii) Many good subsets.* Put `Y = 1 - log L_S/log n ∈ [0,1]`, so `E Y ≥ 1/16 - 1/(M-4)`. For `η ∈ (0, 1/16)`,
`P(Y ≥ 1/16 - η) ≥ E Y - (1/16 - η) ≥ η - 1/(M-4)` (because `Y ≤ 1`). Assume `M ≥ 9 + 2/η`. Then at least
`(η/2)·C(M,5)` subsets
have `L_S ≤ n^{15/16 + η}`, hence `H(u_S) ≤ (1 + o(1)) n^{15/16+η}`.

*(iv) They are A-points at `δ`.* Choose `η = η(δ) > 0` with `(15/16 + η)(1 + 2δ) < 1`; this is possible since `δ < 1/30`.
Then `C n^{-1/4} ≤ (2 n^{15/16+η})^{-(1+2δ)/4}` for `n ≥ n_0(C)`. Also `H(u_S) ≥ H_0` for `n ≥ n_0(C)`: by
Northcott, only finitely many primitive tuples have `H < H_0`, they all have distinct coordinates, so their
distances to `e` are bounded below, while `dist([u_S], e) → 0`. By `A^split_5(δ)`, every good `[u_S]` lies in `Z`.

*(v) Schwartz–Zippel.* `Z ≠ T` is proper, so some nonzero homogeneous `F ∈ Q[u_1..u_5]` vanishes on the cone over
`Z`. The number of ordered `(s_1,…,s_5) ∈ {d_1..d_M}^5` with `F(s) = 0` is at most `deg F · M^4`. Each good subset,
listed in increasing order (the ordering used to define `[u_S]`), is one such tuple. So `(η/2) C(M,5) ≤ deg F · M^4`,
which forces `M ≤ M_1(deg F, η)`. Take
`K = max(9 + 2/η, M_1)`. This depends only on `δ` and `Z`; `C` entered only through `n_0(C)`.

*(vi) #886.* Use the window `(√n, √n + n^{1/2-ε})`: `dist ≤ n^{-ε}(1+o(1))`. Step (iv) needs
`ε > (15/16 + η)(1+2δ)/4`. For `ε > (15/64)(1+2δ)`, choose `η` small enough. `□`

*Remarks.*
(1) Nothing is discarded. The gcd `g_S` is divided out exactly, `H` is the exact toric height, the phase `X` is
arbitrary, and `Z` is the single exceptional set of the hypothesis. It is not a set produced by a fixed-`S` Subspace
theorem.
(2) The same proof with `k`-subsets gives the following (proved, conditional on its hypothesis): **if split `LM_k` is
canonically bounded at `e` for every `k`, then Ruzsa's conjecture #886 holds for every `ε > 0`**. Here
`E log L_S ≤ (1 - 2^{1-k}) log n` and the condition is `ε > (1 - 2^{1-k})/(k-1) → 0`.
(3) Theorem 1 is the split twin of gauss Thm 3.4 / 3.4′, the twisted reduction. The split proof is shorter: the
uniform 5-subset already averages the profile, and no unit class or balanced selection is needed.
(4) No converse is claimed. #887 (absolute `K`) does not obviously imply `A^split_5(δ)`.

## 4. Literature: is the split statement known?

Sources were read through WebSearch snippets and two GitHub pages (WebFetch). `erdosproblems.com`, `numdam.org`,
`math.colgate.edu`, `arxiv.org`, `themoonlight.io` and `althofer.de` do not resolve in this session. **No full text of a
paper was read.**

| item | content (as seen) | status | relation to `A^split_5` |
|---|---|---|---|
| Erdős #886 (Ruzsa), [ErRo97] = Erdős–Rosenfeld, *The factor-difference set of integers*, Acta Arith. (1997) 353–359 | Is the number of divisors of `n` in `(n^{1/2}, n^{1/2} + n^{1/2-ε})` `O_ε(1)`? ErRo: `≤ 1 + C^2` divisors in `[√n, √n + C n^{1/4}]`; infinitely many `n` with 4 divisors in `(√n, √n + 16 n^{1/4})` | open; trivially true for `ε > 1/4` (window `< n^{1/4}`) | Thm 1: `A^split_5(δ) ⇒` #886 for `ε > (15/64)(1+2δ)` |
| Erdős #887 (formal-conjectures `ErdosProblems/887.lean`, fetched) | `parts.ii`: `∃K ∀C>0`, for large `n`, at most `K` divisors in `(√n, √n + C n^{1/4})` | open (the PR below says the main question is unchanged) | Thm 1: `A^split_5(δ) ⇒` `parts.ii`, uniformly in the phase |
| formal-conjectures PR #7218 (fetched; cites an erdosproblems forum post by m-czech, July 2026) | `rosenfeld_4` ("4 is the largest `K` with infinitely many `n`") is false: with `n = (x²-1)(x²-4)(x²-16)(x²-49)`, `z² = 56x² + 65`, there are 5 divisors in the window with `C = 31` (Lean proof with `C = 40`) | claimed solved, with a Lean proof; **re-verified exactly here** (§5, `pell5.py`: `C → 30`) | it is a balanced 15/4 curve in split `LM_5` (§5) |
| Cilleruelo–Jiménez-Urroz, *Divisors in a Dedekind domain*, Acta Arith. 85 (1998) 229–233 | arcs of `x² - dy² = N` of length `N^α`, `α ≤ 1/4 - 1/(8⌊k/2⌋+4)`, contain at most `k` lattice points | snippet | split/real-quadratic Cilleruelo–Córdoba: the Liouville regime (`γ ≤ τ`), not the window |
| Cilleruelo–Jiménez-Urroz, *The hyperbola xy = N*, JTNB 12 (2000) 87–92; Granville–Jiménez-Urroz, *The least common multiple and lattice points on hyperbolas*, Q. J. Math. 51 (2000) 343–352 | bounds for intervals of `xy = N` with `k` lattice points, via lower bounds for the lcm of `k` integers in a short interval | snippet (theorem statements not seen) | the split Liouville / lcm lower bounds; the Euclidean-arc version of §3.2 |
| T. H. Chan, *Numbers with three close factorizations* (arXiv 1303.2249; Integers 18 (2018) A13); Chan–Holmes–Liu–Villarreal, *Numbers with four close factorizations*, arXiv 2508.02818 (2025); arXiv 2507.07094 (2025) | 3 and 4 close factorisations at the vertex, Pell equations, sharp constants (`0.04742` for four) | snippet | vertex (phase 0) only; `k ≤ 4` |
| P. Letendre, *Divisors of an integer in a short interval*, arXiv 2503.12146 (2025); lit18 `direct.md` R3 | `D_n(n^θ, n^η) ≪ τ(n)^{1-ξ(θ,η)}…`; specialised to `θ = 1/2`, `O(1/δ)` divisors in `n^{1/4-δ}` (forum); Conjecture 1 at `√n` | snippet | below the critical scale only (Liouville / CC regime) |
| Ford, *Integers with a divisor in (y,2y]* (arXiv math/0607473); *The distribution of integers with at least two divisors in a short interval* (hal-00091331; authorship not confirmed in the snippet) | distribution in `n` of integers having divisors in dyadic or short intervals | titles only (search results) | statements about how many `n` have such divisors, averaged over `n`; not pointwise, not at the critical scale (assessment) |
| McKinnon–Satriano TAMS 2021; Huang | canonical boundedness is a **hypothesis**; Huang's lower bounds are for split toric varieties but section-type (round 22, row 21) | snippet | give `≤ τ_e` unconditionally, as in the twisted case |

**Verdict (a), search-negative but with a hard implication.** No theorem found states `α^ess_e(-K_{LM_5}) > 15/4`, or
any of its elementary forms in 3.1–3.2. More strongly, by Theorem 1 **`A^split_5(δ)` for any `δ < 1/30` would settle
Erdős #887 (absolute `K`, even uniformly in the phase) and a new range of Ruzsa's #886**. Both are listed as open as of
the July–October 2026 snippets. So the split statement is **not known**, unless #887 has been solved in a source the
snippets did not show. What *is* known in the split case is the Liouville side (Cilleruelo–Jiménez-Urroz, Granville–
Jiménez-Urroz, Letendre's `O(1/δ)`), plus the phase-0 integrality bound `1 + C^2`. That bound is invisible in `T`.

## 5. True or false? Split families at and below the Liouville exponent

**5.1 Theorem 2 (the Pell-balanced hypersurface; proved, with exact and numerical checks as labelled).** Let
`a, b, c ∈ Q`, `d = a + b + c`, and let `S, e_2, e_3` be the power sum and the elementary symmetric functions of
`a², b², c², d²`. Put `m = abcd` and `B = (e_3 - S m)/4`, and let `Q_{abc}: z² = m x² + B` be the conic. Define

    F_0 = (x-a)(x-b)(x-c)(x+d),  F_1 = x⁴ - (S/2)x² + m + 2xz,
    F_2 = (x²-b²)(x²-c²),  F_3 = (x²-a²)(x²-c²),  F_4 = (x²-a²)(x²-b²),

and `φ_{abc} = [F_0 : … : F_4] : Q_{abc} ⇢ LM_5`. Then the following hold.

(i) *(proved, symbolic: `pell_family.py`)* `F_1 · F̄_1 = Π_{r ∈ {a,b,c,d}} (x² - r²) =: n(x)` on `Q_{abc}`, and the
coefficient of `x²` in the conic equation is `(S²/4 + 2m - e_2)/4 = abc(a+b+c) = m`. For each root `r`, `m r² + B` is an
explicit square (printed by the script), so `Q_{abc}` has rational points, e.g. `(a, z_a)`, and `Q_{abc} ≅ P^1_Q`.
All 16 places of `Q_{abc}` over `x = ±a, ±b, ±c, ±d` are rational.

(ii) *(proved symbolically for the vanishing rule; exact at `(a,b,c) = (1,2,4), (1,3,5), (2,3,7), (1,5,6), (1,3/2,5/3)`)*
The 5-vectors of valuations at these 16 places are 0/1. Exactly one place has a constant vector, and the other 15 hit
the 15 cuts of `[5]` **once each**. So `B·C = 15`, the curve is cut-balanced with weight 1, and `E_A·C ∈ {0,1}` with
`E_A·C + E_{A^c}·C = 1`. At infinity all `F_j` have degree 4, so the infinite places contribute 0. The rule
"`F_1` vanishes at `(r, z_r)` iff `z_r = -(r⁴ - (S/2)r² + m)/(2r)`" is a polynomial identity, so the pattern is the same
on the Zariski-open set of parameters where the 16 places are distinct.

(iii) *(proved)* The two points at infinity of `Q_{abc}`, with `z/x → ±√m`, are conjugate over `Q(√m)`. There
`F_j - F_0` has order `x²` against `x⁴`. Symbolically, the `x²`-coefficients of `F_1 - F_0, …, F_4 - F_0` are
`2(z/x) → ±2√m`, `(a+b)(a+c)`, `(a+b)(b+c)` and `(a+c)(b+c)`, which are generically nonzero; the other terms are `O(x)`.
So `max_j |F_j/F_0 - 1| ≍ x^{-2}`. So `C = φ(Q_{abc})` passes
through `e` with two branches of order `μ_0 = 2`, and `mult_e C = 4`.

(iv) *(proved)* If `m` is not a square, then `α_C = B·C/(2μ_0) = 15/4`, attained. Let `t` be a rational parameter
on `Q_{abc} ≅ P^1`, and let `t_0` be the real quadratic parameter of a branch. Then `h_{-K}(φ(t)) = 15 h(t) + O(1)`
(B is base-point free) and `dist(φ(t), e) ≍ |t - t_0|²`. By Liouville on `P^1` at a quadratic point,
`|t - t_0| ≫ H(t)^{-2}`, with equality up to constants along the convergents of `t_0`. So
`dist ≍ H_{-K}^{-4/15}` infinitely often, and never better. Exact check (`w_points.py`): along 60 convergents of
`√m` on `C_{1,2,4}`, `C_{1,5,6}` and `C_{2,3,7}`, `log(dist^{15/4} H)` stays in `[6.8, 10.9]`, `[9.3, 13.7]` and
`[10.6, 25.2]` (periodic), up to `log H ≈ 1500`.

(v) *(numerical, 30-digit evaluation)* The Jacobian of `(a,b,c,x) ↦ [F]` (with `z = +√(m x² + B)`) has rank 3 at 4 random
rational points: the singular values are about `10^{-3}`, plus one of about `10^{-20}`. The `10^{-20}` value is the
scaling kernel, since `F_j` is homogeneous of degree 4 in `(a,b,c,x)`, with `z` of degree 3. So the closure of the
union of the curves is a **hypersurface `W ⊂ T`**, and `W ⊄ Z_diag` (the offsets are distinct).

(vi) *(proved from (iv), (v))* For every `δ < 1/30`, the A-points at `δ` are Zariski dense in `W`. Each curve with
`m` non-square contributes infinitely many, and such `(a,b,c) ∈ Q^3` are Zariski dense in the parameter space. So
**every admissible `Z` in `A^split_5(δ)` or in `(L^split_5)` must contain the hypersurface `W`**, and
`α^ess_e(-K; W(Q)) ≤ 15/4`, with equality along each curve `C_{abc}` (its points are cut-balanced, so the
balanced Liouville bound of §2 applies to them; points of `W` off these curves are not controlled).

**The #887 example is the member `(a,b,c) = (1,2,4)`** (`pell5.py`, `pell5_curve.py`, exact). The five divisors of
`n = (x²-1)(x²-4)(x²-16)(x²-49)` in `(√n, √n + 30.000 n^{1/4}]` are `F_0, F_1, F_2, F_3, F_4`, with offsets
`0, 2√56 - 35, 15, 18, 30` times `x²`. They lie in a window of the critical scale `n^{1/4}`. Their lcm is about `n^{15/16}`,
because one of the 16 places divides none of them. In torus language this is a subcritical (`α = 15/4`) point. This is the
split instance of the phenomenon of the record's Lemma 1.2: the reduction moves a critical cluster to the balanced
Liouville exponent (`α = 15/4`, the twisted `7/15`).

**5.2 What `W` says.**
* In the record's terms, `C_{abc}` is a **balanced, bounded-residue (V-null) rational curve at the Liouville constant**.
  `D_V·C = Σ_cuts |A||A^c|/2 = 40 = ord_e V · mult_e C = 10·4`. The record searched for twisted ones
  (`toric.md` §7; single rational branch, `μ_0 = 8`, 30 roots) and found no nondegenerate real solution. **In the split
  form they exist in a 2-parameter family.** The reason is the parity in 1.4: a split curve can hit each cut once
  (`B·C = 15`). A twisted curve over `Q` hits each cut an even number of times, so cut-balance forces `B·C ≥ 30` and
  then `mult_e C ≥ 8` at `15/4`.
* The family does **not** cover `T`. That is forced: the record's inequality (gauss Thm 4.3(ii); toric Prop 4.3),
  `#fixed·μ_0 + (k-2) E·C̃ ≤ B·C`, applied with the two branch points fixed and `E·C̃ ≥ 2μ_0`, gives `8μ_0 ≤ B·C`,
  i.e. `α ≥ 4` for covering families with two branches. *Label:* this is the record's inequality with
  `#fixed = 2`. I did not re-derive its two-point form; it is the standard free-curve argument with two marked points.
* The record's Gaussian Pell / golden families sit at `α = 4`, not below. For example, every 5-subset of the golden
  8-point family separates all 4 blocks, so `H ≍ φ^{8n}` and `dist ≍ φ^{-2n}` (derived from `growth/families.md` §2).
  The split family sits strictly below the counting value, at exactly `15/4`. I did not compute `α` for the record's
  cyclotomic five-row family (`C → 0` at radii `~10^{3428}`), so I make no priority claim.

**5.3 Could `A^split_5` be false?** For a given `δ`, `A^split_5(δ)` fails iff the A-points at `δ` are Zariski dense
in `T`. The following was checked:
1. *Bounded covering families of rational curves* (one or two branches): excluded at every `α < 4` (Prop 4.3, which is
   geometric and holds for the split form; two-branch form as above).
2. *1-PS orbits and products of coordinate curves.* A 1-PS orbit through `e` with distinct coordinates has
   `α = 2(max v - min v) ≥ 8` (toric Cor 4.2(a), geometric). Group products `y·y'` satisfy
   `h(yy') ≤ h(y) + h(y') + O(1)` and `dist(yy',e) ≲ max(dist(y,e), dist(y',e))`, so they never improve `α`. The
   split-specific objects the task names give nothing below 4.
3. *Polynomial-block curves* `u_j = Π_{A∋j} q_A(x)` with a rational branch at `x = ∞`. These need `μ_0 > B·C/4`. For
   block degree `s` and residual degree `r` this requires `s > r` (2|3 profile) or `s > 2r` (balanced). The coupled
   coefficient conditions are then about `10(6s - r)` or `10(8s - r)` against `20s` or `30s` unknowns:
   overdetermined (heuristic count, consistent with toric §7: expected dimension `2 - 2d`).
4. *Infinitely many non-covering families* (hypersurfaces `W_i`, rigid curves of unbounded degree). This is the only
   way left. It is a K3-type density question. Generalise Thm 2 to `n = Π_{i ≤ m}(x² - r_i²)`, with five divisors of
   degree `m` (in the ring of a conic) agreeing to order `μ_0` at two quadratic branches. Then `α < 4` needs
   `B·C < 8μ_0`. With up to `4m` places, a profile using most of them needs `μ_0 > m/2`. A naive count gives about
   `4(μ_0 - 1) ≈ 2m` coefficient conditions on about `m` roots. But symmetric patterns remove most of them: at `m = 4`
   the even polynomials `F_2, F_3, F_4` and the conic norm `F_1` leave the single condition `d = a+b+c`. So the count is
   **not reliable** and the question is open. Vojta for `(Bl_e LM_5, E, {∞})` predicts that the union of all curves
   with `B·C ≤ (4 - ε)·mult_e C` carrying infinitely many rational points near `e` lies in a fixed proper closed set.
5. *Sporadic points.* Manin–Batyrev–Tschinkel counting for the split torus is `~ c B (log B)^{25}` (Picard rank
   `30 - 4 = 26`). The expected number of points of height `~ B` in the ball of radius `B^{-1/γ}` is
   `B^{1-4/γ+o(1)}`, which is summable for `γ < 4` (**heuristic**).
6. *Data (exact, low height).* `scan.c` covers all primitive 5-tuples with `lcm = N ≤ 2·10^8` inside the invariant
   window `N^{-1/5}`; it finds **no** point with `γ ≤ 4.0`, and the minimum is `γ = 4.277`. This says nothing about `15/4`:
   the `W`-family needs `log H ≳ 60` before `γ < 4` (`pell5.py`). This is the same height-scale obstruction as in the record.

**Verdict (b).** `A^split_5(δ < 1/30)` is **neither proved nor refuted**. The split form is *not* refuted by the extra
split structure: one-parameter subgroups, coordinate products and covering families all stay at `α ≥ 4`. But the split
form carries a **proper, explicit, Zariski-dense-in-a-hypersurface** family of Liouville-extremal points that the
twisted form provably cannot carry with the same profile. So the split exceptional set is strictly larger in a
precise sense. If infinitely many such hypersurfaces exist, the split statement is false while the twisted one could
still be true. That would make the anisotropy of `T_5` essential. I could not decide this.

## 6. Mechanism and transfer to the anisotropic T_5

**6.1 There is no unconditional split mechanism to transfer.** Everything unconditional in the split case stops at
`τ_e(-K) ≤ 15/4`, as in the twisted case. That covers the Liouville and lcm bounds (Cilleruelo–Jiménez-Urroz,
Granville–Jiménez-Urroz, Letendre's `O(1/δ)` below `n^{1/4}`), sections and Seshadri bounds, and Huang's split torsor
bounds (the last by assessment, round 22 row 21). The section cap is a statement about Laurent polynomials and widths, so it holds for every form. The only
genuinely split tools are:
* **Phase-0 integrality** (`d + N/d ∈ Z`, Erdős–Rosenfeld `1 + C²`; Chan's three/four close factorisations). This is
  archimedean information at one value of the phase. It is quotiented out in `T`. Its Gaussian twin is the near-axis
  bound (lit18 §2.3), which the record already has. Using it for a drifting phase would set the phase to zero
  (owner constraint).
* **Fixed-template Pell/Schmidt arguments** (Chan; Temur in the Gaussian case). These give template-dependent,
  ineffective thresholds (owner constraint; lit18 §2.4).

So (c) does not arise: the split statement is not known to be provable, and no split proof exists to transfer.

**6.2 What does transfer (all proved above; map value).**
1. *Reduction.* The split proof of Theorem 1 is a uniform random 5-subset, profile averaging by thresholds, and
   Schwartz–Zippel against one equation of `Z`. It gives the absolute-`K` form with no balanced selection.
   *Assessment, not checked here:* the same three steps should work in the Gaussian setting within one unit class, and
   would show that `A_5(δ)` gives a bound `M(C) ≤ K` with `K` independent of `C`. That would be a mild sharpening of
   gauss Thm 3.4's statement.
2. *Parity lemma (1.4).* A twisted curve over `Q` has even weight on every cut. This is a proved structural
   difference: the cut-balanced weight-1 curves of Theorem 2 have **no** twisted counterpart. A twisted balanced curve
   at `15/4` needs `B·C ≥ 30` and `mult_e C ≥ 8`.
3. *A new twisted search ansatz (suggestion only).* The record searched one rational branch with `μ_0 = 8`
   (`toric.md` §7). Theorem 2 suggests trying **two conjugate real-quadratic branches with `μ_0 = 4`** on a Pell
   conic, with `B·C = 30`: Gaussian blocks on `z² = m x² + B`. By the two-branch covering bound (5.2), such curves
   cannot cover `T_5`. At best they enlarge `Z`, like `W`. Map value only.

**6.3 The calibration (a logical remark, not a meta-theorem).** Suppose a proof of `α^ess_e(-K) > 15/4` works for
every `Q`-form of `LM_5` at once, as Vojta's conjecture would. Then it proves `A^split_5(δ)` for some `δ < 1/30`, and
by Theorem 1 it settles Erdős #887 in its absolute form. So a proof of Question B either is at least as strong as a
solution of #887 (absolute `K`), or uses the anisotropy of `T_5`. Theorem 2 shows that the anisotropy is visible:
the parity lemma excludes in `T_5` the weight-1 balanced curves that exist in the split form. This does not show that
any particular method needs the anisotropy.

## 7. Honest assessment: easier, equal, or dead end

**Not easier. As a source of a mechanism it is a dead end; as a calibration it is useful.**
* *Same geometry, same caps.* `ε_e = 2`, `τ_e ∈ [10/3, 15/4]`, balanced Liouville `15/4`, covering threshold `4`,
  Vojta `4`. Every unconditional tool stops at `τ_e` in both forms.
* *At least as hard as an open Erdős problem.* `A^split_5(δ < 1/30)` implies #887 in its absolute form, uniformly in the
  phase, and a new range of #886 (Theorem 1). Canonical boundedness of all split `LM_k` would imply #886 in full.
* *Strictly more exceptional structure.* The split form has an explicit hypersurface `W` of balanced Liouville-extremal
  points (Theorem 2), which every admissible `Z` must contain. The twisted form cannot have the same curves (parity).
  The split statement could be false, through infinitely many such `W_i`, while the twisted one is true. A proof that
  works for all forms at once would prove the split statement too (6.3).
* *No refutation.* One-parameter subgroups, coordinate products, covering families and polynomial blocks all stay at
  `α ≥ 4` or are overdetermined. The split statement is open in the same window `(15/4, 4)`.

**Probabilities (judgement).**
* The split statement `A^split_5(δ < 1/30)` is known in the literature: < 1% (it would solve #887).
* It is false through infinitely many Pell-type hypersurfaces: 5–10%. (Theorem 2 makes this more plausible than
  before; Vojta says no.)
* Unconditional progress on Question B from this calibration: < 0.5%.
* An Erdős-problems or close-factorisation expert recognising "is `W` the only obstruction?" as tractable: small, but
  this is a new and well-posed question in their language.

**Questions for experts (new, elementary).**
> **B^split (elementary).** Are there a proper Zariski-closed `Z ⊂ T = G_m^5/G_m` and `γ > 15/4` such that only finitely many
> primitive 5-tuples `u` of distinct positive integers with `[u] ∉ Z` satisfy
> `((max u - min u)/min u)^γ · (max u/min u) · lcm(u) ≤ 1`?
> Necessarily `Z ⊇ W` (Theorem 2). A yes implies Erdős #887 (absolute `K`, at every position on the hyperbola).
>
> **W-question.** Is the set of Liouville-extremal (cut-balanced, `α = 15/4`) rational curves through `e` in split
> `LM_5` contained in a proper closed subset? In particular, are there only finitely many "Pell-balanced" patterns
> like Theorem 2's?

## 8. Files and checks

All files are under `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round23/`.
The repository was not edited. Runs used `prlimit --as=3500000000`.

| file | what it checks | label |
|---|---|---|
| `split_checks/pell5.py` | the m-czech #887 family along 40 Pell solutions: 5 divisors in the window, `C_actual = 30.000`, gcd, `log L/log n → ≈ 0.93`, `γ` decreasing to `3.915` at `x ≈ 1.4·10^14` | exact integers |
| `split_checks/pell5_curve.py` | valuation vectors at the 16 places for `(1,2,4)`: one constant place, 15 cuts each once, `B·C = 15`; differences `F_j - F_0` of order `x²` | exact |
| `split_checks/pell_family.py` | Theorem 2 symbolically in `(a,b,c)`: `N(F_1) = n`, `A = m`, `m r² + B` square; profile at 5 parameter points; Jacobian rank 3 at 4 random points | symbolic / exact / numerical (30 digits) |
| `split_checks/w_points.py` → `w_points.out` | rational points along convergents of `√m` on `C_{1,2,4}`, `C_{1,3,5}`, `C_{2,3,7}`, `C_{1,5,6}`: `γ → 3.78…3.79` at `log H ≈ 1000`; `log(dist^{15/4}H)` bounded and periodic over 60 convergents | exact rationals |
| `split_checks/hypergeom.py` | `min_κ p(M,κ) ≥ 1/16 - 2/M` for `5 ≤ M ≤ 3000` (Theorem 1 step (ii)) | exact fractions |
| `split_checks/scan.c` → `scan_3e6.out`, `scan_2e8.out` | all primitive 5-tuples with `lcm = N ≤ 2·10^8` in the invariant window `N^{-1/5}`: min `γ = 4.277`, none `≤ 4` | exact integers, floating `γ` |
| `split_notes.md` | working notes (derivations N1–N8) | — |

Literature access: WebSearch snippets; WebFetch of
`github.com/google-deepmind/formal-conjectures` (`ErdosProblems/886.lean`, `887.lean`, PR #7218). Not resolvable here:
erdosproblems.com, arxiv.org, numdam.org, math.colgate.edu, themoonlight.io, althofer.de.

**Owner-constraint audit.**
* *No uniform exceptional set from a fixed-`S` Subspace theorem.* The only exceptional sets are the hypothetical `Z`
  of `A^split_5` and the explicit `W`. Chan/Temur-type Pell–Schmidt arguments are named and not used (6.1).
* *No gcd or denominator cost discarded.* Tuples are made primitive by exact division. `H` is the exact toric
  height, including the archimedean factor. The rational points of Theorem 2 have their denominators cleared exactly
  (`w_points.py`).
* *No drifting phase replaced by zero.* The phase `log(X/√n)` is arbitrary in Theorem 1. The phase-0 integrality is
  identified and explicitly not used (3.3, 6.1).
* *No template-dependent threshold treated as uniform.* Theorem 1's `n_0(C)` depends on `C`, and nothing claims
  uniformity in `C` beyond what the hypothesis gives. Theorem 2 and the scans make no uniformity claim.
