# Kill report: round22/toric.md

Status: complete (killer report, round 22). Nothing here is refereed.

## 0. Verdict

**SURVIVES as `restricted_theorem`, with corrections. No mathematical error was found in any statement labelled
proved.** I recomputed every number with my own code, and all of them check: Ehrhart polynomials, `(-K)^4 = 70`, the
root-polytope vertices and widths, `Σ|A||A^c| = 80`, `α_bal(k)`, the Vandermonde orders, `t_max(m)` for `m ≤ 7` at
`k = 5`, `m ≤ 9` at `k = 4` and `m ≤ 3` at `k = 3` (ranks mod `2^31-1`, against the author's `2^61-1`), and the
`k = 4` curve data at five parameter points. Nothing in the paper proves the uniform bound, a growth improvement, or
`A_5(δ)` for `δ < 1/30`, and the paper claims none of these.

**What is wrong or overstated** (details in §§2–6):
* **F1 (major overclaim, outside the proved core).** The finite data `t_max(m) = ⌊10m/3⌋`, `m ≤ 7`, do **not** show
  that "toric Liouville with polynomials of degree `≤ 7` cannot improve Cilleruelo–Córdoba's `δ > 1/10`".
  * They cap only a *single* auxiliary polynomial whose maximal cut width is `≤ 7`, not all of `H^0(7B)`. My new
    `m = 8` certificate, `t_max(8) = 26`, extends this to width `≤ 8`.
  * They say nothing about finite *families* used point by point. By Thm 3.3's own minimax, those are governed by
    products of large width.
  * The proved cap `δ ≥ 1/30` is unaffected.
* **F4 (literature/novelty).** McKinnon–Satriano (Trans. AMS 2021, plus a 2023 erratum) define *canonical
  boundedness*, `α(-K) ≥ dim X` for all Zariski-dense sequences, and state that Vojta implies it when `-K` is big. That
  is toric.md's Prop 5.1. Canonical boundedness at `e` implies `A_5(δ)` for all `δ > 0` (§5). So "(McK^ess) is not in
  the literature" is inaccurate, and "not McKinnon's conjecture" needs the qualifier "in its 2007 infimum form".
  * The Huang journal reference (ANT 15 (2021) 461–512) could not be confirmed. The arXiv paper exists and its
    content matches toric.md's description.
* **F2 (minor).** `α_e(B) = 2` is asserted but only `≤ 2` is proved. I supply the missing `≥ 2` (§2).
* **F3 (minor).** "No toric argument can prove `A_5`" is a meta-claim, not a theorem.
  * Prop 6.1 is a tautology about invariants of a transitive action.
  * §6.3 restates the record's fake barrier.
  * Only §6.1 (sections) is a genuine cap.
* **F5, F6 (cosmetic).**
  * "Smooth conic": 1-PS orbit closures can be cuspidal at the boundary.
  * The width form of the forced-divisibility bound is an upper bound for the exact toric height, with equality on
    squarefree profiles, which is enough for the minimax.
  * §7's reformulation should say that coinciding roots are allowed only for nested blocks.

**Novelty assessment.** Several items restate the record in `α`-language:
* the balanced constant `15/4` (= `construct.md` Cor C3 / gauss §4.4);
* the covering-family bound (= gauss Thm 4.3(ii));
* the fake barrier;
* Prop 5.1 (= gauss Thm 4.2, and McKinnon–Satriano's erratum).

Genuinely new and correct:
* `ε_e(-K) = 2`, and it is computed in `Z_diag`;
* the all-profile threshold `τ_e(-K) ∈ [10/3, 15/4]`, its minimax characterisation, and the exact `t_max(m)` data
  (reproduced);
* the identification of gauss Thm 4.1's family as balanced `V`-null curves (reproduced);
* the minimal-degree PTE reformulation at `k = 5`.

**Owner's constraints: respected** (§7). The re-examination of lit18's "drifting phase" dismissal is correct: in `T_5`
the phase is quotiented out exactly.

## 1. Normalisations checked against gauss.md / kill_gauss.md

All checked by hand. Scripts are in `kill_toric_checks/`; they are my own code and import nothing from `toric_checks/`.

* **Height.** `h_B([ζ]) = W' = log R''^2` for the primitive tuple. This is gauss Lemma 1.1, re-derived in kill_gauss
  §4. The symmetric split `h_{E_A} = h_{E_{A^c}} = ½ log n_A` follows from `σ(1_A) = 1_{A^c}`. **Correct.**
* **Dictionary (Lemma 1.3).**
  * Write `dist ≤ diam/R''` and `H_B = R''^2`. Then `diam ≤ R''^{1/2-δ}` gives `dist ≤ H_B^{-(1+2δ)/4}`.
  * So `A_5(δ)` gives `α ≥ 4/(1+2δ)` off `Z`.
  * At `δ = 1/30` this is `4/(16/15) = 15/4`, and `δ = 2/α - 1/2` matches gauss Prop 3.1 (`1/2 - 1/(2^k-2)`).
  * The converse (b): if `A_5(δ')` fails for every `Z`, the A-points off `Z ∪ Z_diag` form a sequence with
    `R'' → ∞`, so `α ≤ 4/(1+2δ')`. **Correct.**
  * The balanced part (c) uses `d_ab/2 = 4w + O(ηw)`: each pair is separated by 8 of the 15 cuts. It also uses
    `log|t_ab| ≥ 0` within one unit class, as kill_gauss §1.2 item 2 corrected. The equivalence of `(L_5^∘)` with
    "`μ ≤ 4 - c'` off `Z`" holds with `c' ≍ c/10` and `η_0` small. **Correct.**
* **No phase is replaced by zero.** `[ζ] = (ζ_j/ζ_1)_j` is exactly invariant under `ζ ↦ γζ`, and `dist(y,e)` is
  computed from the genuine classes. So the lit18 objection is answered correctly in `T_5`; see §7.
* **gcd and denominators.** Theorem 3.2 keeps every power of `N'`: `f(y) = G(ζ)/N'^m` with `G ∈ Z[i]`. Points are
  classes of primitive tuples. **Correct.**
* **Covering inequality (Prop 4.3).** This is gauss Thm 4.3(ii)'s `#fixed·μ_0 ≤ B·C - (k-2)E·C̃`, which kill_gauss
  §4 checked. With one fixed point and `E·C̃ ≥ μ_0` it reads `-(K_Y+E)·C̃ ≥ 0`, i.e. `B·C ≥ (k-1)μ_0`. At `k = 4`
  it is the calibration `12 ≤ 14` of gauss §4.3. **Correct.**

## 2. Claim (i): Seshadri constant, root polytope, (-K)^4

**Recomputed** (`basics_own.py` → `basics_own.out`). Lattice points of `mR_k`, enumerated by brute force:

| `k` | point counts, `m = 0, 1, …, k` | Ehrhart polynomial (exact interpolation) | `B^{k-1}` | `binom(2k-2, k-1)` |
|---|---|---|---|---|
| 3 | 1, 7, 19, 37 | `3m^2+3m+1` | 6 | 6 |
| 4 | 1, 13, 55, 147, 309 | `10m^3/3+5m^2+11m/3+1` | 20 | 20 |
| 5 | 1, 21, 131, 471, 1251, 2751 | `35m^4/12+35m^3/6+85m^2/12+25m/6+1` | **70** | 70 |

* The interpolations were checked against brute force at two further values of `m`.
* For all 120 permutations, the vertex of the max cone is the root `e_{s_k} - e_{s_1}`.
* `{⟨x,1_A⟩ ≤ 1 ∀A} = {‖x‖_1 ≤ 2}` holds on a box.
* `width_v(R_k) = 2(max v - min v)` was checked for all `v ∈ {-2..2}^k`.
* Independently, for `k = 4`: `Bl` of `P^3` at 4 points and then 6 lines gives `(-K)^3 = 64 - 32 - 6·2 = 20`.

**Theorem 2.1 (`ε_e(B) = 2`): correct.**
* Upper bound: `B·C_{1_A} = width = 2`, and the curve is smooth at `e`.
* Lower bound: `π_1` is an isomorphism over the torus, so `C → π_1(C)` is birational. A curve in `P^4` has degree
  `≥ mult` at any point, so `H·C ≥ mult_e C`, and likewise `H'·C ≥ mult_e C`.
* `C_{1_A} ⊂ Z_diag` at `k = 5` (since `|A| ≥ 2` or `|A^c| ≥ 2`).
* Galois invariance of `ε` holds.

The statement that every 1-PS orbit descends to a `Q`-curve is correct, because `σ = -1` on `N`.

**Minor gap (F2).** Cor 2.2 asserts `α_e(B) = ε_e(B) = 2` "on the whole variety". The text proves only `≤ 2`, via
`C_{1_A}`. The lower bound `α_e(B) ≥ 2` is not proved there, and McKinnon–Roth's Liouville theorem is not invoked
for it. It is true. My proof:
* `f_{ab} := (x_a - x_b)(x_a^{-1} - x_b^{-1}) = 2 - x_a/x_b - x_b/x_a` has support `{0, ±(e_a - e_b)} ⊂ R_5`, so
  `f_{ab} ∈ H^0(B)`. It has integer coefficients and `ord_e f_{ab} = 2`.
* On norm-one points, `|f_{ab}(y)| = |u_a - u_b|^2`.
* For `y ≠ e`, choose `a, b` with `|u_a - u_b|` maximal. Then `f_{ab}(y) ≠ 0` and `|u_a - u_b| ≤ 2 dist(y,e)`.
* Theorem 3.2 (first form, `m = 1`) gives `|u_a - u_b|^2 ≥ H_B^{-1}`, so `dist ≥ ½H_B^{-1/2}`, i.e. `α ≥ 2` on all of
  `T_5(Q)`.

**Minor (F5).** Lemma 1.2 calls `C_v` a "smooth `Q`-conic". The closure of a 1-PS orbit is smooth at `e`, but at
its two boundary points it is a monomial curve `(t^{a_1}, …)`, where the `a_i` are the coefficients of `v` in its braid
cone. For the braid fan the `a_i` are the consecutive gaps of the sorted entries of `v`. The curve is singular there when
every gap is `≥ 2`. Example: `v = (0,2,4,6,9)` is primitive in `Z^5/Z1` (its gaps `2,2,2,3` have gcd 1), and the
branch is `(t^2,t^2,t^2,t^3)`, a cusp. Only smoothness at `e` is used, so nothing changes.

## 3. Claim (ii): approximation constants, 15/4, k=4 family, k=5 covering families

* **Prop 4.1** (`α = B·C/μ_0` along a rational curve with contact at a rational `O`): **correct.** It needs
  Liouville on `P^1` at a rational point, together with `h_B(ν(t)) = (B·C)h(t) + O(1)` (`B` is base-point free). The
  balanced bookkeeping `w = 2κh`, `λ = μ_0 h`, `B·C = 30κ` gives `μ = 15μ_0/(B·C)`.
* **Cor 4.2(a): correct.** If all `v_j` are distinct then `max v - min v ≥ 4`, so `α ≥ 8`.
* **Prop 2.3 / Cor 4.2(b): correct.**
  * `V = Π(x_a-x_b)/(x_1⋯x_5)^2` has support `Perm(2,1,0,-1,-2)`: 120 monomials, max `ℓ_1 = 6`, so `V ∈ H^0(3B)`.
  * The widths are `4` (singleton cuts) and `6` (2|3 cuts), i.e. `|A||A^c|`, so
    `Σ_cuts |A||A^c| = 80 = k(k-1)2^{k-3}`. Hence `mult ≤ 8κ` and `α ≥ 30/8 = 15/4`.
  * For general `k`, `α_bal(k) = (2^k-2)·C(k,2)/(k(k-1)2^{k-3}) = 4(1-2^{1-k})`, recomputed exactly for `k = 3..8`.
    Its `δ`-value `2/α - 1/2` equals `1/(2^k-2)` exactly, matching gauss Prop 3.1. **Confirmed.**
* **Thm 4.4: confirmed exactly** (`k4curve_own.py` → `.out`, sympy). The singleton roots were rebuilt from gauss
  §4.1's fourth-root formula at five parameter points `(X,m) = (-3,3), (2,5), (1/2,3/2), (-5,1), (7,2)`. At each:
  * `Im P_x` is constant: `(-420, -180, -240)` at `(-3,3)`, matching gauss and kill_gauss;
  * the 14 roots are distinct, and each of the 14 rays is met exactly once;
  * `B·C = 14`, the contact at `∞` is `4`, so the ratio is `7/2`;
  * `D_{V_4}·C = 24 = 6·4`.
* **Prop 4.3 and the `k = 4 → 5` switch: correct** (see §1). Covering families give `α ≥ k-1`. The balanced
  Liouville constant `4(1-2^{1-k})` is `> k-1` only for `k ≤ 4`. This is a restatement of gauss Thm 4.3(ii)
  (killer-checked) in `α`-language, not a new theorem.
* **§7 reformulation: correct, with one omission.**
  * Each balanced class-`κ = 1` curve with its contact at `s = ∞` has `u_j = c_j P_j/P̄_j`.
  * Here `c_j = 1` because `u_j(∞) = 1`.
  * A point `s_0` meeting the boundary has valuation in a braid cone. Balance forces every ray to occur once with
    coefficient 1, so coinciding roots are allowed **only for nested blocks**. A coincidence `r_T = r_{T'}` with
    `T, T'` not nested would give `1_{T∪T'} + 1_{T∩T'}` and unbalance the class.
  * The text does not mention this. It is harmless, because nested coincidences lie in the closure of the
    parametrisation that was searched.
  * The count is `30 - 2 = 28` unknowns against `4·7 = 28` equations. This gives `2 - 2d` in general and `1 + 2d` at
    `k = 4`.
  * The argument principle gives 4 roots of each `R_j` in each half-plane.
  * **Correct.**
* **§7 numerics (evidence): independently reproduced in kind.** `pte_own.py` is my own formulation: least squares on
  the 28 equations with gauge `r_{[4]} = i`.
  * *Calibration at `k = 4`:* 100/100 converged, and 24 were nondegenerate.
  * *`k = 5`* (`pte_own2.py 5 80 11 3000` → `pte_k5b.out`): 44/80 converged and **0 were nondegenerate**. In every
    converged solution some roots coincide: the largest minimal root distance is `2·10^{-12}`. This **matches** the
    author's finding in kind. A longer run (`pte_own.py 5 400 7`, `max_nfev = 20000`) was stopped after about 30 CPU
    minutes without output and is not counted. Evidence only. Nothing depends on it.

## 4. Liouville threshold tau_e(-K) and the section computation

**Lemma 3.1, Thm 3.2: correct.**
* Theorem 3.2 (first form) is the integrality `G(ζ) ∈ Z[i]`, keeping `N'^m`.
* In the second form, for a polytope `Q_f` whose normal fan does not coarsen the braid fan, `D_{Q_f}` is not a
  Cartier divisor on `X`. The correct exact finite local height is `-ψ_{Q_f}(trop_v y)`. It is `≤ Σ_A a_A·`(coefficient
  of `1_A` in `trop_v y`) by concavity of `ψ_Q`, so the stated inequality holds as an *upper bound*.
* It is an equality on squarefree profiles (one layer per prime). There it is all that a toric height can give,
  which is what the minimax needs (F6, wording only).

**Thm 3.3 (minimax): correct.**
* `K = {width_f/ord_e f}` over the semigroup consists of mediants, so its closure is convex.
* Sion's theorem applies because the simplex of profiles is compact.
* The identity `inf_F δ_F = 2/τ_e - 1/2` uses three facts:
  * `Q_f - Q_f ⊂ wR_5` with `w = max_A width_f(1_A)`;
  * `ord(ff^*) = 2 ord f`;
  * Fekete, since `t_max` is superadditive, so `τ = sup_m t_max(m)/m`.
* All are checked.

**Thm 3.4(a): correct.**
* Lower bound: `V^a·V_{abc}^b` gives `⌊10m/3⌋` for every `m`. Recomputed for `m ≤ 9`.
* Upper bound: `construct.md` Cor B1 (homogeneous form, refereed) with `M = 5`, `s = 0`:
  `ord ≤ 2^{-4}·30·2m = 15m/4`.
* So `2/τ - 1/2 ≥ 1/30` for all profiles. This is a mild sharpening (all profiles, any finite family) of what
  `construct.md` Cor C3 / gauss §4.4 already gave for balanced profiles. **Survives** as stated.

**Thm 3.4(b)(c): reproduced exactly with independent code** (`regrank.c`, own C, ranks mod `p = 2^31-1`; the author used
`2^61-1`):

| `k` | `m` | `|S_m|` | `d` | rank at `d` | rank at `d-1` |
|---|---|---|---|---|---|
| 5 | 1, 2, 3, 4, 5 | 21, 131, 471, 1251, 2751 | `⌊10m/3⌋` = 3, 6, 10, 13, 16 | full | 15, 111, 470, 1245, 2731 (deficient mod p) |
| 5 | 6 | 5321 | 20 | **5321 full** (all 10626 monomials) | n/c |
| 5 | 7 | 9381 | 23 | **9381 full** (all 17550 monomials) | n/c |
| 5 | 8 | 15421 | 26 | **15421 full** (all 27405 monomials; **new**, beyond the author's range) | n/c |
| 4 | 1..9 | 13 … 2869 | `3m` | full | deficient |
| 3 | 1..3 | 7, 19, 37 | `3m` | full | (m=3) deficient |

Over `Q`, full rank mod `p` implies full rank. Deficiency at `d-1` mod `p` is not a proof over `Q`, but the explicit
Vandermonde products give the lower bound anyway. **Thm 3.4(b)(c) confirmed.**

**F1 (major as an overclaim; the proved statements are unaffected).** In §0 ("the open intermediate exponent item of
gauss.md §7 is out of reach of toric Liouville with polynomials of degree `≤ 7`") and in Remark (1) and §6.1 ("by
Theorem 3.4(b) no section of `mB`, `m ≤ 7`, beats the Vandermonde, so on the data even the intermediate range
`δ ∈ (1/30, 1/10]` is out of reach of toric Liouville"), the task summary says: "For `m ≤ 7` it does not even improve
Cilleruelo–Córdoba's `δ > 1/10` off a closed set". The data do not support this as stated.

1. **For a single `f`,** what matters in Thm 3.3 is `max_A width_f(1_A)/ord_e f`, not `ord_e f/m`.
   * A section `f ∈ H^0(7B)` has `ord_e f ≤ 23`, but its maximal cut width `w` can be anything up to `14`.
   * The data cap `w/ord ≥ 3/5` only through `ff^* ∈ H^0(wB)`, i.e. only for `w ≤ 7`. That class contains all of
     `H^0(mB)` only for `m ≤ 3`.
   * An `f ∈ H^0(7B)` with `w = 11` and `ord ≥ 19` would beat `3/5`. It is not excluded by `t_max(7) = 23`; it would
     need `t_max(11) ≥ 38 > ⌊110/3⌋ = 36`, which is not computed.
2. **For a finite family `F`** with each `f` used separately, point by point, the guaranteed exponent is
   `min_{w∈conv K_F} max_A w_A` (Thm 3.3's own minimax).
   * Interior points of `conv K_F` correspond to products `Π f_i^{a_i}`, whose widths add.
   * So even if every `f_i` has width `≤ 7`, the hull is not controlled by `t_max(m ≤ 7)`.
   * Two polynomials of small width that peak on different cuts are capped only by `t_max` at larger `m`.

What the data prove (the author's `m ≤ 7`, extended here to `m = 8`: `t_max(8) = 26`) is this: *every single Laurent
polynomial `f` with `max_A width_f(1_A) ≤ 8` (equivalently, `ff^* ∈ H^0(8B)`) gives at best `δ > 1/10`.* For finite families and larger widths, "toric Liouville does not
improve Cilleruelo–Córdoba" is exactly the conjecture `τ_e = 10/3`, which the paper itself labels evidence. The proved
cap `δ ≥ 1/30` (Thm 3.4(a)) is unaffected.

## 5. Claim (iii): equivalence A_5 <-> alpha^ess > 15/4; Vojta; McKinnon-Roth; citations

**Mathematics: correct.**
* The equivalence `A_5(δ<1/30) ⟺ α^ess_e(B; Ω_dist) > 15/4` is Lemma 1.3; see §1.
* `α^ess ≥ τ_e`: each `f` contributes its own closed set.
* `τ_e ≤ 15/4` (Thm 3.4(a)), hence "beyond every section constant".
* Prop 5.1, `K_Y = π^*K_X + (n-1)E`, gives `nλ_E ≤ h_{-K} + εh_A`, and `h_A ≤ c h_{-K}` off a closed set because `-K`
  is big. So `α^ess ≥ n`. **Correct.**

**Citations (re-verified by my own searches; queries and snippets are in `kill_toric_notes.md`).**

| citation in toric.md | my query | finding |
|---|---|---|
| McKinnon–Roth, Invent. Math. 200 (2015) 513–583 | "McKinnon Roth Seshadri constants, diophantine approximation, and Roth's theorem for arbitrary varieties Inventiones" | **confirmed**: vol. 200, 513–583, DOI 10.1007/s00222-014-0540-1, arXiv 1306.2976 |
| their Liouville companion, Eur. J. Math. 2016 | "An analogue of Liouville's theorem and an application to cubic surfaces" | **confirmed**: Eur. J. Math. 2 (2016) 929–959, arXiv 1306.2977. Snippet: "when the locus being approximated is defined over the ground field, Liouville's Theorem is strictly better" than Roth; they also give a version with the asymptotic base locus (a τ-type, section-based bound). This is consistent with toric.md §5 |
| Huang, Algebra & Number Theory 15 (2021) 461–512 | "Zhizhong Huang Rational approximations on toric varieties Algebra & Number Theory"; "msp.org ant 2021 15-2 Huang …" | **journal reference NOT confirmed** by any snippet. arXiv 1809.02001 exists (v5, 30 Apr 2021): universal torsor method, "general fixed point", **split** toric varieties, best approximations on minimal-degree rational curves (minimal centred primitive collections), "under certain hypotheses". The content matches toric.md's description; the ANT volume and pages are **unverified** |
| Vojta, LNM 1239, Conj. 3.4.3 | (inherited from gauss.md; not re-searched) | not re-verified |

**F4 (minor–major: literature and novelty).** toric.md does not cite **McKinnon–Satriano, "Approximating rational points
on toric varieties", Trans. AMS 374 (2021) 3557–3577 (arXiv 2004.05212), with a 2023 erratum**. By the search
snippets (my queries: "McKinnon Satriano Approximating rational points on toric varieties abstract …"; "canonically
bounded OR canonical boundedness McKinnon Satriano …"):
* They define *`X` is canonically bounded at `P`* as `α_{P,{x_i}}(-K_X) ≥ dim X` **for every Zariski-dense sequence**.
* The erratum states that, **when `-K_X` is big, Vojta's Main Conjecture implies canonical boundedness at `P`**.
* They prove, for `Q`-factorial terminal **split** toric varieties, that `P` is better approximated on a rational
  curve than by any Zariski-dense sequence. For surfaces they prove McKinnon's conjecture conditional on Vojta.
* Related: Monahan–Satriano, arXiv 2308.11847 (horospherical, Vojta-conditional); Lehmann–McKinnon, arXiv 2403.02480
  ("essentially bounded").

These are all snippet-level; I could not read the papers (arXiv, ar5iv and math.uwaterloo.ca do not resolve here).

Consequences:
* toric.md **Prop 5.1 is, at snippet level, the corrected McKinnon–Satriano Prop 1.3.** It should be cited, not
  presented as new.
* The sentence "(McK^ess) … is not in the literature that I could check" is **inaccurate**. The essential-constant
  framework comparing Zariski-dense sequences with rational curves is exactly the McKinnon–Satriano/Lehmann–McKinnon
  setting.
* The claim "it is not McKinnon's conjecture" is true only for the 2007 infimum form. In the
  McKinnon–Satriano form ("any sequence that approximates `P` sufficiently well must lie on a rational curve"), the
  relevant notion is canonical boundedness.
* My own short argument shows that canonical boundedness of `L̄_5^tw` at `e`, even restricted to the `Q`-points
  `T_5(Q)`, implies `A_5(δ)` for **every** `δ > 0`, hence the uniform bound:
  * if the A-points for `δ` are not Zariski dense, their closure is the `Z` of `A_5(δ)`;
  * if they are dense, an enumeration is a Zariski-dense sequence tending to `e` with `α ≤ 4/(1+2δ) < 4 = dim X`.
* So the expert question can be stated in the literature's own terms: *is the twisted Losev–Manin 4-fold canonically
  bounded at the identity (or at least `α ≥ 15/4 + c` for Zariski-dense sequences of norm-one points)?*
* None of this changes a mathematical conclusion of toric.md. The literature results are conditional (Vojta), are for
  split varieties, or concern the infimum. But the novelty and literature claims of §5 must be corrected.

**Seshadri-proportional statement: correct with qualification.** "Any bound `α ≥ c·ε_e`, `c ≤ 1`, gives `≤ 2`" is
right for the unrestricted `ε_e`. Prop 2.3's restricted constant satisfies `10/3 ≤ ε^diag ≤ 8`, so a hypothetical
"essential Seshadri" theorem `α^ess ≥ ε^ess` is **not** capped at 2 by this argument. toric.md does not claim it is,
but §0's "far above anything Seshadri-proportional (`≤ 2`)" should name the unrestricted constant.

## 6. Claim (iv): "no toric argument" meta-claims

The task summary states (iv) as "No toric argument can prove `A_5(δ)`, `δ < 1/30`". **This is not a theorem**, because
"toric argument" is not defined. What is proved is narrower:
* **6.1 (sections):** proved. Thm 3.3 plus Thm 3.4(a) cap every finite family of auxiliary Laurent polynomials used
  through Liouville with toric forced divisibility at `δ ≥ 1/30` for arbitrary profiles. For balanced profiles this is
  exactly the cube-profile cap `construct.md` Cor C3 / gauss §4.4 / kill_gauss §5 item 3, already in the record.
  The extra sentence about `δ ∈ (1/30, 1/10]` "on the data" is F1.
* **6.2 (Prop 6.1):** correct but **tautological**. A function invariant under a group acting transitively on
  `T_5(R)` is constant there. The conclusion that proximity "enters only through non-invariant objects (characters
  and sections, i.e. 6.1)" is an unproved identification of all non-invariant inputs with sections (F3).
* **6.3 (cocharacter/tropical data):** a restatement of the record's fake barrier (gauss §5 item 1, kill_gauss §5
  item 1; killer-checked, not refereed). It is correct as a restatement, and its status is inherited.
* **6.4 (curves):** covering rational families have `α ≥ 4` (proved, gauss Thm 4.3(ii)). Non-covering families,
  unbounded degree, higher genus with finitely many points per member, and sporadic points are not excluded, and
  toric.md says so.

So the correct label for (iv) is: *barrier for three named classes (sections, `T(R)`-invariant archimedean data,
fake-blind tropical data), plus a restatement of the curve-family exclusion*. It is not "no toric argument".

## 7. Owner's constraints

1. *Uniform exceptional sets from fixed-`S` Subspace:* none used. The exceptional sets are `Z(f)`, the hypothetical
   `Z` and Vojta's `Z_ε` with `S = {∞}` fixed, conditionally. **Respected.**
2. *Gaussian gcd or denominator costs:* none discarded. Primitive tuples are used, `N'^m` is kept in Thm 3.2, and
   `h_B = log R''^2` is exact. **Respected.**
3. *Drifting phase replaced by zero:* no. The quotient by `R^1G_m` removes the common phase exactly. The re-examination
   of lit18's dismissal is correct: in `T_5` the target `e` is genuinely fixed. The reason fixed-target tools still
   fail (the infimum `α_e = ε_e = 2` is computed in `Z_diag`, and section-based constants are `≤ τ_e ≤ 15/4`) is
   correct. **Respected.** One caveat: Huang's and McKinnon–Satriano's theorems are for **split** toric varieties. Over
   `Q(i)` the variety splits, but applying them to `Q(i)`-points changes the height normalisation and the set of
   points (non-norm-one points). toric.md uses none of them, so nothing is violated.
4. *Template-dependent threshold treated as uniform:* none. `t_max(m)` is used only as a cap on a method, and the
   PTE search is a refutation attempt. The only problem is F1: a finite-`m` cap is presented as a cap on a broader
   class than it covers. That is an overclaim about a method, not a uniformity violation.

## 8. What survives (corrected statement)

Let `X = L̄_5^tw`, `B = -K_X`, `e` the identity, and `Ω_dist = T_5(Q) ∖ Z_diag`. Heights are `h_B([ζ]) = log R''^2` for
primitive tuples, and `dist(y,e) = max_j|ζ_j/ζ_1 - 1|`.

1. **(proved)** `B` has polytope `R_5`; `B = H + H'`; `B` is nef and big but not ample; `B^4 = 70`. Also
   `ε_e(B) = 2`, attained only by the conics `C_{1_A} ⊂ Z_diag`, and `α_e(B) = 2` (the lower bound by the argument of
   §2 here).
2. **(proved)** `A_5(δ)` for some `δ < 1/30` ⟺ `α^ess_e(B; Ω_dist) > 15/4`. `(L_5^∘)` is the same over `η_0`-balanced
   one-unit-class sequences.
3. **(proved)** Along rational curves, `α = B·C/μ_0`.
   * Balanced curves off `Z_diag` have `α ≥ 15/4`, with equality iff they are `V`-null.
   * Covering rational families through `e` have `α ≥ 4` (gauss Thm 4.3(ii)).
   * Gauss Thm 4.1's family consists of `V`-null balanced curves at `k = 4` (`B·C = 14`, contact 4, `α = 7/2`).
4. **(proved)** Toric Liouville is governed by `τ_e`.
   * Finitely many auxiliary Laurent polynomials, used through Liouville with toric forced divisibility off their zero
     sets, give for arbitrary profiles exactly `A_5(δ)`, `δ > 2/τ_e - 1/2`.
   * `10/3 ≤ τ_e ≤ 15/4`, so this never passes `δ = 1/30`.
   * **(exact)** `t_max(m) = ⌊10m/3⌋` for `m ≤ 7`, reproduced here; `m = 8` is new here.
   * **Corrected consequence:** a single auxiliary polynomial with maximal cut width `≤ 8` gives at best `δ > 1/10`.
     The same for finite families or larger widths is equivalent to `τ_e = 10/3`, which is open (evidence only).
5. **(conditional on Vojta, `-K` big)** `α^ess ≥ 4`, i.e. canonical boundedness in the sense of McKinnon–Satriano.
   This implies `A_5(δ)` for every `δ > 0`, and hence the uniform bound. This is gauss Thm 4.2 and (snippet-level)
   McKinnon–Satriano's corrected Prop 1.3. It is not new.
6. **(barrier, three named classes; not "all toric arguments")**
   * sections: proved;
   * `T(R)`-invariant archimedean data: trivial;
   * fake-blind tropical data: the record's barrier.
7. **(evidence)** No nondegenerate real solution of the minimal-degree `k = 5` PTE system has been found, in the
   author's search or in mine (§3); the `k = 4` calibration finds solutions readily.

Not proved, and not claimed: the uniform bound, any growth improvement, `A_5(δ)` for `δ < 1/30`, `(L_5^∘)`, and
`τ_e = 10/3`.

## Scripts

All under `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round22/kill_toric_checks/`,
run under `prlimit --as=3500000000`. All are my own code.

| script | output | content |
|---|---|---|
| `basics_own.py` | `basics_own.out` | Ehrhart polynomials, `B^{k-1}`, vertices, widths, `Σ|A||A^c|`, `α_bal(k)`, Vandermonde support/orders |
| `regrank.c` (gcc -O3) | `reg_small.out`, `reg_m6.out`, `reg_m7.out`, `reg_m8.out` | exact ranks mod `2^31-1` of monomial evaluation on `mR_k ∩ M` |
| `k4curve_own.py` | `k4curve_own.out` | sympy exact check of Thm 4.4 at five parameters |
| `pte_own.py`, `pte_own2.py` (same, flushed, `max_nfev` argument) | `pte_k5b.out` (k=5: 44/80 converged, 0 nondegenerate); k=4 calibration printed in the session (100/100 converged, 24 nondegenerate) | least-squares search for minimal-degree `V`-null balanced curves |

Literature queries and snippets are in `../kill_toric_notes.md`. arXiv, ar5iv and math.uwaterloo.ca do not resolve in
this session, so all literature findings are at abstract/snippet level.
