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

**Exact data (`profiles.py`).** `Π_nt` has dimension 11; 787 vertices were found (coordinates in
`{0, 1/3, 1/2, 2/3, 1}`), 278 of them 0/1 (16 `S_5`-classes); the minimum binomial cost over all of them, for
`|m_a| ≤ 3`, is exactly 4. The one-sided classes include the **code vertices** `x_R = 1[|B ∩ R| even]` for
`|R| = 3` (the class of `W`: in round 23's labels `R = {2,3,4}`, i.e. the three even polynomials `F_2, F_3, F_4`)
and `|R| = 5` (even sets). *(proved: `|B∩R| + |B^c∩R| = |R|` is odd, so exactly one side of each cut is even, and
`X_a = 8` for all `a`.)*

**Meaning.** The record transferred the balanced section cap `15/4` to the split form at the symmetric profile only.
A split near-tight point has an 11-dimensional profile, and the section inequality is profile-dependent
(Lemma 3.0). Theorem 3.2 shows that for every auxiliary form with centrally symmetric Newton polytope — all binomial
products, all Vandermonde products, all equality cases of Theorem B found by the record (Conjecture E) — the cost
does not depend on the profile at all. The split asymmetry is a direction along which these forms are exactly flat.

### 3.3 Proposition 3.3: the cap at the `W`-vertex reduces to one polynomial (proved)

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

*Proof.* (a) The curves `C_{abc}` (over `C`) form a 2-parameter family sweeping `W`, each with `E_B·C = 1_{B∈F_R}`
and `mult_e C = 4` (two branches of order 2; round 23 Thm 2). If `g_W ∤ f` then `f` does not vanish on a general
`C = C_{abc}`; pulling back to the normalisation `P^1`, the zeros of `f` over `e` number at least `ord_e(f)·mult_e(C)`
and all zeros number at most the poles, `Σ_B π_B(f)(E_B·C) = cost_{x_R}(f)`. (b) `ord_e` and `π_B` are additive
under products (Newton polytopes add), so write `f = g_W^k h` with `g_W ∤ h` and apply (a) to `h`. The restricted
theorem follows from Lemma 3.0 and continuity of `cost_x` in `x`. (c) as stated. `□`

So at the `W`-vertex the cap is **proved for every auxiliary form except multiples of `g_W`**, and it is decided by a
single inequality `cost_{x_R}(g_W) ≥ 4·ord_e(g_W)` (`≥ 32` if `ord_e g_W = 8`). Since `deg W > 14`, every small-support
search (§3.4) is blind to this case. I could not compute the Newton polytope of `g_W` (`deg W ≤ 128` by the
intersection count on `Y ⊂ P(1,1,1,1,3)`, minus base-locus corrections; interpolation beyond degree 14 exceeds the
memory limit). This is the one place where sections could still beat `15/4` on split near-tight points, and it is a
finite algebraic computation (§7, Q-split-1).

### 3.4 Conjecture WTB (weighted Theorem B on `Π_nt`) and evidence

> **Conjecture WTB.** For every finite `S ⊂ M` and every near-tight split profile `x ∈ Π_nt`:
> `4·reg(S) ≤ Σ_B x_B π_B(S)`. Equivalently: no Liouville/section argument beats `15/4` at any near-tight split
> profile.

Proved: at `x = 1/2` (Theorem B); for centrally symmetric `S` (Theorem 3.2); at `x_R` for `S` not containing the
support of a multiple of `g_W` in the sense of Prop 3.3(a); and trivially whenever `Σ_B (x_B − 1/2)mid_B(S) ≤ 0`.
*Evidence (exact mod `2^31 − 1`; mod-`p` rank is a lower bound for the rational rank, so mod-`p` `reg` is an upper
bound for the true `reg` and cannot produce false violations):*
* random polytopes `Q(p) = {⟨m,1_B⟩ ≥ −p_B}` (`p ∈ {0..3}^30`, and lopsided `p`), 3 790 sets with `|S| ≤ 1500`, LP
  minimum of `cost_x(S)` over `Π_nt`: maximum of `4reg/mincost` **exactly 1.0000** (`search_beat.py`);
* local search from root zonotopes with random deletions/insertions (`local_search.py`, 2×600 s): maximum of
  `4reg − mincost` exactly `0`;
* one-sided free polytopes `Q_F(n) = {⟨m,1_B⟩ ≥ −n, B ∈ F}`: for the `W`-class, `reg/n = 3, 3.5, 3.33, 3.5, ≤ 3.4`
  (`n = 1..5`), for the class `1^5 3^{10}` `reg/n ≤ 3.33`; normalised volumes `4!·vol(Q_F(1)) ∈ [89.2, 118.3]` for all 16
  classes, below the `(15/4)^4 = 197.75` that a dimension count would need (`qf_reg.py`, `qf_reg_big.py`, `onesided_vol.py`);
* the quadric-cone form `(u_2−u_1)(u_3−u_1) − (u_4−u_1)(u_5−u_1)` (ord 2): minimum cost over `Π_nt` = `25/3 > 8`
  (`cone_section.py`);
* **k = 4 control** (`k4_wtb.py`): the same statement with `2·reg ≤ cost_x` on the `k = 4` near-tight polytope,
  18 875 random sets, maximum ratio exactly 1.

*Why a proof is not immediate (derived).* Theorem B's slicing step pairs `w_{P'}(S) + w_{P'∪1}(S)` to recover the
widths of two slices; one-sided costs `π_B` give only two of the four extremes, so the induction does not close for
non-symmetric weights. A proof needs a genuinely asymmetric slicing lemma.
