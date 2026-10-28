# Referee report on round23/split.md (adversarial)

Referee scripts: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round23/referee_split_checks/`
(all written independently of `split_checks/`; outputs `*.out`). Notes: `referee_split_notes.md`. Repository read only.

## 0. Verdict summary

**The result survives.** Status `partial` is right. No growth improvement and no uniform bound follow, and split.md does
not claim either.

* **Theorem 1** (split `A_5(δ<1/30)` ⇒ Erdős #887 in the absolute-`K` form, uniform in the phase; ⇒ #886 for
  `ε > (15/64)(1+2δ)`). I re-derived every step with all quantifiers and found it **correct**. There are two cosmetic
  gaps: the Northcott wording and the integer window of the Lean statement (§2).
* **Theorem 2** (the m-czech construction is a cut-balanced curve at exactly 15/4, in a 2-parameter family sweeping a
  hypersurface `W`). Every symbolic identity was **re-verified independently**. The one step that was only numerical,
  Jacobian rank 3 in (v), I have now **proved exactly**: exact rational rank 3 at 4 rational points, plus an exact
  weighted-homogeneity invariance that gives rank ≤ 3. So (vi), "every admissible `Z` contains `W`", is now fully proved.
* **Parity lemma.** Correct. I re-derived it through Galois orbits. It holds in the stronger form `E_A·C = E_{A^c}·C`
  for every `A`.
* **One genuine error (major, but outside the main results).** §5.3 item 2 says that group products "never improve α".
  The inequalities given for it point the wrong way. The claim is also **false, and it contradicts Theorem 2**: every
  point of `T(Q)` is a product of points on the coordinate one-parameter subgroups, and these have `α ≥ 8`. The words
  "coordinate products do not refute it" must go from the result.
* Minor issues: a wrong offset `2√56 − 35` (it should be `2√56 ≈ 14.97`). The claim "`log H ≳ 60` before `γ < 4`" is not
  reproduced: I find `log H ≈ 106` at the earliest. That is harmless because it is a weaker lower bound. Small
  citation omissions are listed in §6.

## 1. Translation: geometric ↔ elementary (heights, distances, balance)

**Checked and correct.**
* *Polytope.* `P_{-K} = {m ∈ M : ⟨m,1_A⟩ ≥ -1 ∀A}` is the polar of `conv(1_A)`. Also
  `conv(1_A) = {n ∈ N_R : max n − min n ≤ 1}`, which is the polar of `conv(e_a − e_b)`. So `P_{-K} = R_5`, as stated.
* *Height (1.2).* The local factor at `p` is `p^{max v_p − min v_p}`, and the product over `p` is `lcm/gcd = lcm`. At
  `∞` the factor is `max u/min u`. Re-derived. Round 22 J2 checked it exactly on 300 tuples.
* *Cut bookkeeping (1.4).* For each `p` with `min_j v_p = 0` (forced by `gcd = 1`), the layers `s = 1..max v_p` are
  proper and nonempty. So `lcm = Π n_A`, `u_a = Π_{A∋a} n_A` and `gcd(u_a,u_b) = Π_{A⊇{a,b}} n_A`.
  - Small caveat: `h_{E_A} = log n_A + O(1)` holds near `e`, where the archimedean proximity to `E_A` is bounded. It does
    not hold globally. Only points near `e` are used, so nothing breaks.
* *Balanced Liouville 15/4 (§2).* Sum over the 20 ordered pairs: `Σ_A |A||A^c| log n_A`. The cut weights are
  `5·4 + 10·6 = 80`, so the best pair gives `−log dist ≤ 4w`, against `log H = 15w`. Correct. It agrees with J2's
  `30/(15−7)`.
* *A^split_5 versus α^ess.* `A`-points at `γ` increase with `γ` because `dist < 1`. Hence "α^ess > γ ⇒ A(δ) ⇒ α^ess ≥ γ"
  with `γ = 4/(1+2δ)`. Correct.
* *3.1 and 3.2.* "Exact" means exact up to the factors `(1 + D/U)` and `max u/min u = 1 + o(1)`. These are absorbed
  because `δ < 1/30` is an open condition. The hyperbola dictionary is correct:
  - `gcd(N/d_j) = N/lcm(d)`;
  - with the Minkowski form `dx dy` and `x = √N e^s`, the arc length is `√N |Δs|`;
  - `γ = 4` corresponds to arc length `N^{1/4}`.
  In the vertex statement, `N` is the lcm of the five divisors, not the `n` of #887. Theorem 1 handles that gap
  correctly, so this matters only for the heuristic sentence.
* *3.3, phase-`p/q` bound.* With `s = qd + pN/d` and `X² = pN/q`, one gets `s − 2qX ≈ qX h² ≤ C²√(pq)`. The values of `s`
  are distinct, because `dd' = X²` is impossible for two distinct divisors `≥ X`. So the bound is `O(1 + C²√(pq))`, as
  stated.

## 2. Theorem 1 (split A_5 ⇒ Erdős #887 / #886)

**Correct.** I re-derived each step.
* **(ii)** `v_p(L_S) = Σ_{s≤t} 1[S meets both sides]`, using `L_S = lcm(d_S)/g_S`, which divides `n`.
  - The one-side probability is `p(M,κ)`. Its second difference is `C(κ,3) ≥ 0`, and `C(x,5)` is increasing for `x ≥ 4`.
    So `min_κ p ≥ 2C((M−1)/2,5)/C(M,5) = 2^{-4}Π(M−1−2i)/(M−i) ≥ 2^{-4}(1 − Σ(i+1)/(M−i)) ≥ 2^{-4}(1 − 15/(M−4))`.
  - Independent exact check (`hyper_indep.py`): **0 violations** of either `2^{-4}(1−15/(M−4))` (for `9 ≤ M ≤ 3000`) or
    `1/16 − 2/M` (for `5 ≤ M ≤ 3000`).
* **(iii)** Since `Y ≤ 1`, `P(Y ≥ a) ≥ EY − a`. With `M ≥ 9 + 2/η` we have `1/(M−4) < η/2`. Correct.
* **(iv)** `(15/16 + η)(1+2δ) < 1 ⇔` (for small `η`) `δ < 1/30`. The inequality
  `C n^{-1/4} ≤ (2n^{15/16+η})^{-(1+2δ)/4}` holds for `n ≥ n_0(C)` **independently of `X`**. Monotonicity in `H` is
  fine, since a smaller `H` only loosens the condition.
  - *Wording gap (minor).* "Only finitely many primitive tuples have `H < H_0`; they all have distinct coordinates" is
    false as written: `(1,1,1,1,1)` has `H = 1` and `dist = 0`. The fix is immediate. Apply Northcott to the finite set
    of **distinct-coordinate** primitive tuples with `H < H_0`. Their distances are bounded below by some `c_0 > 0`,
    while `dist(u_S) ≤ C n^{-1/4} < c_0` for `n ≥ n_0(C)`.
* **(v)** Schwartz–Zippel on the grid `{d_i}^5` against one nonzero form `F` that vanishes on the cone over `Z`.
  Distinct subsets give distinct sorted tuples. `K = max(9 + 2/η, M_1(deg F, η))` is independent of `C` and `X`.
  Correct.
* **Match with the formal statement.** I fetched `raw.githubusercontent.com/.../ErdosProblems/887.lean` myself.
  `parts.ii` is `∃K ∀C>0 ∀ᶠ n, #{d ∈ Ioo ⌊√n⌋₊ ⌈√n + C n^{1/4}⌉₊ | d ∣ n} ≤ K`. That integer window lies inside
  `[√n, √n + (C+1)n^{1/4}]`, so apply Theorem 1 with `C + 1` (minor; this should be said). I also fetched `886.lean`:
  `∀ε>0 ∃K ∀ᶠ n, … ≤ K`. Step (vi) gives exactly this for `ε > (15/64)(1+2δ)`.
  - Arithmetic check: `(15/64)(1 + 2/30) = 1/4`, so the range really goes below `1/4`. The range `ε > 1/4` is already
    known from `rosenfeld_bound` (`1 + C²`).
* **Remark (2), `k`-subsets.** The inequality `E log L_S ≤ (1 − 2^{1-k}) log n` holds only up to `+O(1/M)` (checked:
  `k = 7`, `M = 50` gives `min p = 0.0096 < 2^{-6}`). The `η`-slack absorbs this. The threshold
  `ε > (1 − 2^{1-k})/(k−1)` is correct, given canonical boundedness (`α ≥ k−1` on Zariski-dense sequences).
* **Is it the main logical point?** Yes. Theorem 1 is a clean, correct proof that split `A_5(δ<1/30)` is at least as
  strong as the open #887 (absolute `K`). So "the split statement is not known" is justified, modulo the open status,
  which I verified (§6).

## 3. Theorem 2 (m-czech construction; family C_abc; hypersurface W)

* **(i)** Exact, symbolic in `a, b, c` (`thm2_sym.py`):
  - `F_1 F̄_1 ≡ Π(x² − r²)` modulo the conic;
  - the `x²`-coefficient of the conic equals `m`;
  - `m r² + B` is a perfect square for every `r ∈ {a, b, c, d}`. For example,
    `m a² + B = ((a²b + a²c + ab² + ac² − b²c − bc²)/2)²`.
* **(ii) Profile.** I reproved it by hand, and the hand proof does not depend on which sign of `z` makes `F_1` vanish.
  Over each of the 8 values `x = ±a, ±b, ±c, ±d` lie two places. Their valuation sets are `S_x` and `S_x ∪ {1}`, where
  `S_a = {0,3,4}`, `S_{−a} = {3,4}`, `S_b = {0,2,4}`, `S_{−b} = {2,4}`, `S_c = {0,2,3}`, `S_{−c} = {2,3}`, `S_d = ∅`
  and `S_{−d} = {0}`. Reducing modulo complement gives `∅` once and every 1- and 2-subset of `{0..4}` exactly once.
  `F_1` has no other zeros because `F_1 F̄_1 = n(x)`. All `F_j` have degree 4 at the two points at infinity.
  - Exact, with series in the local parameter (`thm2_sym.py`), at the 6 points `(1,2,4)`, `(1,3,5)`, `(2,3,7)`, `(1,5,6)`,
    `(1,3/2,5/3)` and `(3,5,11)`: one constant place, all 15 cuts hit exactly once, multiplicities ≤ 1. **`B·C = 15`
    confirmed.**
* **(iii)** Symbolically, the `x⁴` and `x³` coefficients of `F_j − F_0` vanish. The `x³` term cancels because
  `d = a+b+c`. The `x²` coefficients are `(a+b)(a+c)`, `(a+b)(b+c)` and `(a+c)(b+c)`.
  - For `F_1`: `F_1 − F_0 = 2xz + 0·x² + O(x)`. The pure `x²` term vanishes because `−S/2 − e_2(a,b,c,−d) = 0`, so the
    order-`x²` term is `2xz ~ 2√m x²`.
  - Hence the branch order is 2 at each conjugate point at infinity, and `mult_e = 4`. Correct.
  - **Typo:** split.md §5.1 lists the offsets as `0, 2√56 − 35, 15, 18, 30`. The second offset is **`2√56 ≈ 14.967`**
    (`2√56 − 35 ≈ −20.03` would put `F_1` outside the window). `split_notes.md` line 48 carries the same slip.
* **(iv)** `h_{-K}∘φ = 15h + O(1)` (bpf), `dist ≍ |t−t_0|²` and Liouville at a real-quadratic `t_0` give `α_C = 15/4`.
  This is correct. Any branch type gives `mult_e ≥ 8` for a twisted analogue (§4).
  - Independent exact check (`conv_indep.py`): 70 convergents of `√m` with the slope parametrisation through
    `(a, z_a)`. The ranges of `log(dist^{15/4}H)` are `[6.81, 10.82]` for `(1,2,4)`, `[9.27, 13.61]` for `(1,5,6)` and
    `[10.56, 25.16]` for `(2,3,7)`. These are **identical** to split.md's ranges. The heights reach `log H` up to
    `1916`, and the last values are `γ = 3.771…3.791`.
  - Integral Pell points of `z² = 56x² + 65` (`pell_indep.py`, 118 points, `log H ≤ 766`):
    `log(dist^{15/4}H) ∈ [6.00, 13.95]`, `C_actual → 30.0000` and `log L/log n → 0.934…`, which tends to `15/16`.
* **(v) Upgraded to exact (`jac_exact.py`).** Parametrise the conic rationally by the slope `t` through `(a, z_a)`:
  `x = (a(t²+m) − 2z_a t)/(t² − m)`.
  - The map `(a,b,c,t) ↦ (F_j/F_0)_{j=1..4}` is **exactly invariant** under `(a,b,c,t) ↦ (λa,λb,λc,λ²t)`, verified
    symbolically. So its rank is ≤ 3 everywhere.
  - Its Jacobian has **exact rank 3** at the 4 rational points `(1,2,4,3/7)`, `(1,3,5,−2/5)`, `(2,3,7,11/3)` and
    `(1/2,5,6,1/9)`.
  - Hence the closure `W` of the union of the `C_{abc}` is an irreducible hypersurface of `T`: the total space of the
    conic bundle is irreducible because the generic conic is nondegenerate. This is now **proved**, not numerical.
* **(vi)** Each `C_{abc}` with `m` a non-square carries infinitely many A-points at every `γ > 15/4`. A closed `Z` that
  contains all but finitely many of them contains `C_{abc}`. The parameters with `m` a non-square are Zariski dense, so
  `Z ⊇ W`. Correct, and with (v) upgraded it is fully proved. The points lie on the positive real branch and have
  distinct coordinates, so they satisfy the "distinct positive integers" clause of `A^split_5`.
* **Scope.** Theorem 2 is not a refutation. `W` is proper, and A-points on `W` off the curves are not controlled, as
  split.md says.
* **Not mentioned in split.md.** PR #7218 also states that the same Lean file proves "infinitely many `n` with five
  divisors already for `C = 7/2`" (no published source). That is a second, different five-divisor family. It may give
  another Pell-balanced surface `W'`. It is unexamined and bears on the "only `W`?" question in §7.

## 4. Parity lemma

**Correct.** Let `K = Q(i)`. `Gal(Q̄/K)` acts trivially on `N`, and `σ` acts by `−1`. A `Gal(Q̄/Q)`-orbit `O` of boundary
places splits into at most two `Gal(Q̄/K)`-orbits, `O_1` and `σO_1`, carrying the vectors `v` and `−v`.
* If `O = O_1`, then `v = −v`, so `2v = 0`, and `v = 0` because `N` is torsion-free.
* Otherwise the vectors `v = Σ c_B 1_B` (along a chain) and `−v = Σ c_B 1_{B^c}` (along the complementary chain) occur
  equally often.

Hence `E_A·C = E_{A^c}·C` for every `A`, even for non-ray vectors. So the weight on every cut is even.

Consequences:
* A cut-balanced twisted curve has `B·C ≥ 30`.
* At `α_C = 15/4`, `mult_e C ≥ 8` for every branch type:
  - a rational branch gives `α = B·C/μ_0`, so `μ_0 ≥ 8`;
  - an algebraic branch of degree `k ≥ 2` gives, by Roth, `α = B·C/(2μ_0)` and `mult_e ≥ kμ_0 ≥ 8`.

split.md states `mult_e ≥ 8` without the case split, but the claim is correct.

## 5. Truth-value discussion, scan data, heuristics

* **5.3 item 1 (covering families, `α ≥ 4`).** This rests on the record's inequality
  `#fixed·μ_0 + (k−2)E·C̃ ≤ B·C`. In its two-point form it gives `2μ_0 + 3·2μ_0 = 8μ_0 ≤ B·C`, and in its one-point form
  `4μ_0 ≤ B·C`. split.md says openly that the two-point form was not re-derived. I did not re-derive it either.
  *Label: record, plus an unverified extension.* That is acceptable as stated.
* **5.3 item 2, 1-PS orbits: correct.** For `t ↦ t^v` with `v_j` distinct, `B·C = 2(max v − min v) ≥ 8`. The branch is
  rational of order 1, so `α = B·C ≥ 8`.
* **5.3 item 2, group products: ERROR (major for the claim, no effect on Theorems 1 and 2).** The text says
  `h(yy') ≤ h(y) + h(y') + O(1)` and `dist(yy') ≲ max(dist(y), dist(y'))`, "so they never improve α".
  - These are **upper** bounds on both `H(yy')` and `dist(yy')`. They bound the ratio `log H/(−log dist)` of `yy'` from
    **above**, roughly by twice the worse ratio of `y` and `y'`. They give no lower bound, because cancellation
    (`y' ≈ y^{-1}w`) can bring `yy'` arbitrarily close to `e`.
  - The claim is in fact **false**: `[u] = Π_j ε_j(u_j)`, where `ε_j` is the `j`-th coordinate 1-PS. Each factor lies on
    a 1-PS orbit through `e` with `α ≥ 8`. Yet Theorem 2's points are such products, and they have `α = 15/4`. So the
    principle contradicts split.md's own Theorem 2.
  - "Coordinate products do not refute it" is therefore **vacuous**: every point is a coordinate product. It should be
    deleted from §5.3, §7 and the result summary. The verdict "undecided" is unaffected.
* **5.3 item 3 (polynomial blocks).** A heuristic count, labelled as such in the body. split.md itself shows that the
  count is unreliable (item 4: symmetry collapses it at `m = 4`). The result summary says "do not refute it (a heuristic
  count says overdetermined)". That is acceptable only with the word *heuristic* kept.
* **5.3 item 5.** The Manin-type count `B(log B)^{25}` (`ρ = 26`) gives `B^{1−4/γ}` per ball. Heuristic, and labelled so.
* **5.3 item 6 (scan).**
  - I re-read `scan.c`. The window `N^{-1/5}` is a valid necessary condition because `H ≥ N`. Divisors are enumerated
    with a smallest-prime-factor sieve, and the gcd and lcm conditions are imposed exactly.
  - Independent Python enumeration (`scan_indep.py`, different method, `N < 2^20`): `2^18`: count 2, min `γ = 4.9399`
    at `(336,340,350,357,360)`; `2^19`: count 4, min `γ = 4.7170`. This **matches `scan_2e8.out` exactly**. I did not
    rerun the full `2·10^8` range. Min `γ = 4.277` is therefore *reported*, cross-checked only on the low ranges.
  - "The `W`-family needs `log H ≳ 60` before `γ < 4` (`pell5.py`)" was not reproduced. On the `(1,2,4)` integral Pell
    points the first `γ < 4` is at `log H = 105.8`, and on convergents of the four members at
    `log H = 139…225`. As a lower bound, "≳ 60" is still true, and the conclusion (the scan is uninformative) is
    reinforced.
* **Vojta ⇒ split true.** Consistent with the McKinnon–Satriano erratum, which I saw as a snippet: Vojta implies
  canonical boundedness for toric varieties (big `−K`). Canonical boundedness gives `α^ess ≥ 4` and hence
  `A^split_5(δ)` for every `δ > 0`. The probabilities in §7 are judgement and are labelled as such.

## 6. Citations / literature

These were my own checks, made on 2026-10-28. `erdosproblems.com`, arxiv/ar5iv, `jtnb.centre-mersenne.org` and
`galoisrepresentations.org` do not resolve here. `raw.githubusercontent.com` does resolve.

| item in split.md | my check | result |
|---|---|---|
| formal-conjectures `887.lean` | fetched raw: `parts.i`/`parts.ii` are `research open`; `rosenfeld_4` is now the negation and marked solved | **confirmed** (open status as of the fetch) |
| `886.lean` | fetched raw: `erdos_886` open; `rosenfeld_infinite` (4 divisors, `16n^{1/4}`) and `rosenfeld_bound` (`1+C²`) solved | **confirmed** |
| PR #7218 | GitHub search: closed 2026-10-22. Body: m-czech forum post (July 2026), `n = (x²−1)(x²−4)(x²−16)(x²−49)`, `z² = 56x²+65`, `C = 31`, Lean `C = 40`; "main question `parts.i`/`parts.ii` still open"; also claims five divisors for `C = 7/2` | **confirmed**. I re-verified the construction exactly (`pell_indep.py`). The `C = 7/2` claim was omitted by split.md |
| Erdős–Rosenfeld [ErRo97] | Acta Arith. **79** (1997) 353–359, doi 10.4064/aa-79-4-353-359 | confirmed (split.md omits the volume number) |
| Cilleruelo–Jiménez-Urroz, Acta Arith. 85 (1998) 229–233 | bibliographic record confirmed; content not seen | the content claim (`α ≤ 1/4 − 1/(8⌊k/2⌋+4)`) is **snippet**, as split.md labels it |
| Cilleruelo–Jiménez-Urroz, JTNB 12 (2000) 87–92 | confirmed: arcs of `xy = N` containing `k` lattice points | confirmed |
| Granville–Jiménez-Urroz, QJM 51 (2000) 343–352 | confirmed bibliographically | confirmed, content unseen |
| Chan 1303.2249; 2507.07094; Chan–Holmes–Liu–Villarreal 2508.02818 | confirmed. The 0.04742 figure comes from a KSU symposium abstract, not the arXiv paper | confirmed, secondary for 0.04742 |
| Letendre 2503.12146 | confirmed; Conjecture 1 at `√n` | confirmed |
| McKinnon–Satriano TAMS 374 (2021) 3557–3577, erratum | canonical boundedness is a hypothesis; erratum: Vojta ⇒ canonical boundedness for toric (big `−K`) | confirmed (snippet) |
| **missing** | Chan, *Factors of almost squares and lattice points on circles*, IJNT 11 (2015), arXiv 1406.2230 (already in the record, lit18 R5); Chan, *Perfect squares have at most five divisors close to its square root*, arXiv 1303.2069 | these should be cited. They are phase-0 or special-`n` results and do not change the verdict |

No search turned up a construction with ≥ 6 divisors in the `C n^{1/4}` window, or a resolution of #887 (search-negative).
That point matters: a refutation of #887's absolute-`K` form would, by Theorem 1, **refute** split `A_5(δ<1/30)`.

## 7. Owner constraints

* *Fixed-`S` Subspace theorem.* Not used. The only exceptional sets are the hypothetical `Z` and the explicit `W`.
* *Gaussian gcd or denominator costs.* `g_S` is divided out exactly and `H` includes the archimedean factor. My
  convergent check clears denominators exactly.
* *Drifting phase.* Theorem 1's `n_0(C)` is independent of `X` (verified in §2). Phase-0 integrality is identified and not
  used.
* *Template-dependent thresholds.* `K = K(δ,Z)` and `n_0(C)`. Nothing is claimed uniform beyond the hypothesis.

**Passes.**

## 8. Corrected statement

See the structured output (`corrected_statement`). It equals split.md's result with the following changes:
1. (v) and (vi) of Theorem 2 are now fully proved (exact Jacobian rank).
2. The phrase "coordinate products do not refute it" is deleted. The group-product principle is false.
3. The offset of `F_1` is `2√56`.
4. The Northcott step in Theorem 1 is restricted to distinct-coordinate tuples, and the Lean window is matched with `C+1`.
5. The `C = 7/2` five-divisor claim of PR #7218 and Chan 1406.2230 / 1303.2069 are added to the literature list.

## 9. Scripts

All are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round23/referee_split_checks/`:
* `thm2_sym.py` — Thm 2 (i)–(iii) symbolic; exact profile at 6 parameter points.
* `jac_exact.py` — exact homogeneity plus exact Jacobian rank 3 (Thm 2 (v)).
* `pell_indep.py` → `pell_indep.out` — m-czech integral Pell points, exact.
* `conv_indep.py` → `conv_indep.out` — convergent points on 4 members, exact.
* `scan_indep.py` → `scan_indep.out` — independent low-height scan.
* `hyper_indep.py` — hypergeometric bound for Thm 1 (ii).

Runs used `prlimit --as=3500000000`.
