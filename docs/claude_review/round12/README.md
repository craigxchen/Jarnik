# Round 12 — Judge's README (route 5 at box size T = 1)

## 0. Bottom line

**No growth improvement and no uniform bound were proved this round**, either unconditionally or from a standard
conjecture. `M <= (2/3 + eps) log R/loglog R` (`two_thirds/proof.md`) is still the best bound for actual circles.

**Route 5 at box size `T = 1` is now rigorous up to a single, explicitly stated, open hypothesis about the real
tangent cone of the cluster's own slab box. That hypothesis is not proved for any regime box, and no computation
bears on it, in either direction.**

What this round established. Everything is refereed unless it is marked J (this README's own unrefereed findings,
§4).
1. **The `T = 1` machinery is proved, for every profile and every norm.** This covers:
   * the slab box: `m_max >= 0.6417(n-1) - 2.2137`, and `m_max >= 4w` once `n - 1 >= 316 ln 4r`;
   * rigidity for every order in `[ceil(4w), m_max]` with explicit `W_1`;
   * Step 0;
   * transversal vanishing of order `m' - rho_A(F)` along every projection-cheap span;
   * the structure of all cheap families: at `m' = ceil(4w)` at most 4 long directions and 16 classes, so cheap
     spans lie in the braid arrangement;
   * the **arithmetic separation lemma**: `dist(hat u, L) >= (1/2) exp(-binom(n,2) abar/(2P(L)))` for every
     subspace `L` with `P(L)` collision pairs, for any cluster on any norm.

   This proves the judge's round-11 sketches J3 and §5.3, and it uses no squarefree, rank or nondegeneracy
   hypothesis.
2. **The R11 target "ZE-TC" (cheap-span form) is false**, including at the recommended order:
   * in an admissible, zonotope-free, overdetermined, LP-supercritical regime profile (Sylvester `2^12` plus one
     light column, `m' = ceil(4w) = 2256`);
   * in the faithful models `M8` (`m = 8`) and `M8 x M8` (`m = 16`).

   All failing points are degenerate. The correct cost of a rational subspace is the projection cost / tangent
   degree `nu_A`, not a width sum.
3. **The final hypothesis has been located, and it is quantitative.** Step 0, the separation lemma and one
   Lojasiewicz-type inequality on the profile's own box give `n = O(loglog R)` **on all norms** (Theorems 4.1/6.8,
   proved implications).
   * The weakest form that suffices is the referee's profile-specific **MIN^pair**.
   * Every form that depends only on dimension, degree, height and the qualitative zero set provably cannot supply
     the needed size (nearstrata Prop 2.7).
4. **Obstructions at finite `n`.**
   * MIN\* is false on the regime Sylvester box `n = 4096` at `g > 328` (referee).
   * **J5:** the uniform-`Q_0` collision hypothesis that the t1rigour referee recommended (and the author's
     (Q1)+(Q2)) is false on the Sylvester boxes `n = 2^12, ..., 2^17` for **every** `Q_0`, at every level
     `abar >= 5.42`. That includes the whole near-extremal range. With per-subspace radii, the coordinate-span
     test points refute it at `n = 2^12` (`W >~ 2.7 W_min`), but not at `n = 2^13` or `2^14`.
   * These are finite-`n` failures, so they push `n_1` up but do not touch the asymptotic statements.
5. **The fake test is now a theorem-level reduction.** A nondegenerate real tangent-cone point of a (Mar) fake pattern
   box is the direction of gap-admissible fakes (fake Thm 2.1). At the barrier scale this also needs a Diophantine
   condition. So "route 5 beats the barrier" is equivalent to (quantitative) ZE-TC on the fake pattern boxes
   themselves.
   * The round's computational evidence that overdetermination kills such points did not survive refereeing. Three
     of the four "overdetermined" models are zonotopal. The sdc11 cone did not shrink: a certified zero has
     `sep >= 0.106`, 64% of the maximum. The landscape optimiser misreports minima.
6. **J6 (unrefereed, exact):** the tangent cone has an algebraic description.
   * The leading-form ideal of a one-cut box is exactly `(k_x^3, sigma_S^(2L+1))`. That is the strong Lefschetz
     property of `k[x]/(x^3)`, and it gives a short proof of zetc Thm 2.1's `rho`.
   * In the faithful model `M8` the leading-form ideal is generated, in all degrees `<= m`, by width powers and
     pullbacks from **3-column projections**.
   * By apolarity, `TC(A,m)` always contains the real base locus of an explicit fat-point linear system. Cheap spans
     are exactly its linear base cycles.

**Scores** (rubric of `round10/README.md` §1; details in §6):

| write-up | (a) | (b) | (c) | **overall** |
| --- | --- | --- | --- | --- |
| t1rigour | 8 | 3 | 3 | **5** |
| nearstrata | 7 | 3 | 3 | **4** |
| zetc | 7 | 2 | 2 | **4** |
| fake | 4 | 3 | 2 | **3** |

**Recommended next step** (§7): establish the **local structure of the leading-form ideal** of regime `T = 1` boxes
in degrees `<= ceil(4w)`. Then use it to decide the strata half of the per-subspace collision hypothesis: every
subspace in the tangent cone must satisfy `o_L binom(n,2)/P(L) < d`.

## 1. Setting and notation

Round 12 had four write-ups, each attacked by an adversarial referee whose report sits next to it:

| write-up | task |
| --- | --- |
| `t1rigour.md` | make R11's J1-J3 and the `T = 1` reduction rigorous; cheap families; separation |
| `zetc.md` | the qualitative tangent-cone statement ZE-TC: formulation, few cuts, products, small models |
| `nearstrata.md` | the quantitative (Lojasiewicz) half, separation, assembly |
| `fake.md` | the adversarial fake test for the `T = 1` programme |

Conventions are those of `round10/README.md` and `round11/README.md` §1. *Proved* means proved with all quantifiers
in a refereed write-up. *Exact* means integer or exact-rational computation, or mod-`p` computation used only in the
valid direction. *Evidence* means floating point, or a finite experiment that no proof uses. **J**-items are this
README's own findings: nobody but the judge has checked them, so they are **unrefereed**. Scripts are in
`judge_checks/` and notes in `judge_notes.md`.

**`T = 1` objects.** For a cluster of `n` points (`d = n - 1`) with `r` split primes and profile columns
`beta_p in [-1,1]^n`:
* `w = sqrt(8 d ln 4r)` and `j = floor(w/2)`;
* `A = {k in Z^n_0 : |k_x| <= 1 (x < n), |<k,beta_p>| <= w}`;
* `m_max` is the largest `m` with `8 binom(m+d-1, d) <= 3^d`, and orders `m' in [ceil(4w), m_max]`;
* `H = U^2`; `V_H` is the span of the integer `c in V(A,m')` with `|c|_inf <= H`;
* `Phi_c(u) = sum_k c_k <k',u>^(m')`; `kappa_c(u) = |Phi_c(u)|/(|c|_1 d^(m'))`; and
  `TC_H = {u : Phi_c(u) = 0 for all c in V_H}`.

Further notation:
* `J(A)` is the leading-form ideal of the vanishing ideal of `A`. The tangent degree is
  `nu_A(u) = max{i : <k,u>^i not in J(A)}`, and `TC(A,m) = {nu_A < m}` (zetc Lemma 1.1).
* `rho_A(F)` is the interpolation degree of the image `F_*(A)` (projection cost), and `o_F = m' - rho_A(F)`.
* `d_F = sum h_A(v_i)` is the width cost; `F` is *cheap* if `d_F < m'`.
* Pair excesses: `eps_xy = d_xy - W/2 + 2 ln C' >= 0` and `abar = W/(2d) + 2 ln C'`, with
  `sum eps_xy <= binom(n,2) abar`.
* `Gamma(u) = (2/N_2) sum ln^+(1/|u_x - u_y|)`, with `N_2 = binom(n,2)`.
* `P(L)` is the number of collision pairs of a subspace `L`.
* `eps_0 = C d e^(-W/4)/(m'+1)` is the Step-0 level.

## 2. Referee-corrected statements

### 2.1 t1rigour (`t1rigour.md`, `referee_t1rigour.md`): restricted theorem, survives

**Setting.** `T = 1`, `w = sqrt(8(n-1) ln 4r)`, `j = floor(w/2)`, `H = U^2`. All statements hold for every profile,
every norm and every base point. No squarefree, rank or nondegeneracy hypothesis is used.

1. **Slab box (proved; constants interval-certified; thresholds exact).**
   * `|A| >= 3^d/2` and `m_max > 0.6417 d - 2.2137`.
   * `dim V_H >= |A|/2` for every `m' <= m_max` (`n >= 5`).
   * `m_max >= 4w` whenever `d >= 316 ln 4r`. The exact thresholds are `2900 / 3623 / 5809` for
     `r = n, n ln n, n^2`.
2. **Rigidity (proved).** Rigidity holds for every `m' in [ceil(4w), m_max]` once
   `W > W_1 = 24 d ln 3/m' + 8 ln(e d C'/m')`. This is at most `52.71 + 8 ln C'` at `m_max`, and
   `O(sqrt(n/log r))` at `ceil(4w)`.
3. **Step 0 (proved).** `kappa_c(hat u) <= eps_0 = C d e^(-W/4)/(m'+1)` for every `c in V_H`.
4. **Transversal vanishing (proved).**
   * Along the span of any family with `rho_A(F) < m'`, every `Phi_c` vanishes to order `>= m' - rho_A(F)`, with
     ceiling `kappa <= min(2^(1-m'), (2+eta)^(m') eta^(m'-rho_A(F)))`.
   * A local converse holds, with ineffective constants.
   * New (referee, Lemma R1): `rho_A(F) >= 2 min(j, n - largest F-class)`.
5. **Cheap families (proved).** For any order:
   * total window length `<= floor((m'-1)/(2j))`;
   * at most `(m'-1)/2` exceptional points;
   * at most `2^b` core classes.

   At `ceil(4w)` (`w >= 10`): at most 4 long directions and at most 16 classes (both attained), and cheap spans lie
   in the braid arrangement. At Sylvester orders this fails exactly when `m' > 2j log_2 n`, which happens inside the
   window iff `n = 2^k` with `k >= 16`. The tail `k >= 23` is the referee's estimate.
6. **Separation (proved).** `dist(hat u, L) >= (1/2) exp(-binom(n,2) abar/(2P(L)))` for every subspace `L` with
   `P(L) >= 1`. In particular the distance to every cheap span at `ceil(4w)` is at least
   `(1/2) exp(-16(1+O(w/n)) abar/2)`, which is `n^(-O(1))`.
7. **Refutation (proved).** The cheap-span form of ZE-TC is false.
   * Witness: Sylvester `2^12` plus one light column, at `m' = ceil(4w) = 2256 <= m_max = 2631`. The model is
     admissible, zonotope-free (Parseval: `zeta <= 2w(n-1)/n < m'`), overdetermined and LP-supercritical.
   * The `n = 7` model of Thm 6.11(a) is only a toy: it contains a zonotope of length 8 and lies outside the regime.
8. **Conditional reduction (proved implication; referee's reformulation).** The hypothesis
   ZE-TC^coll_eps`(A, ceil(4w); Q_0)` implies `n = O(log W) = O(loglog R)` on all norms. It says: for all profiles
   with `n >= n_1` and `d >= 316 ln 4r` and all `W > W_1`, `kappa(hat u) > eps_0` whenever `hat u` is at distance
   `>= (1/2) e^(-Q_0 abar/2)` from every collision subspace with `>= binom(n,2)/Q_0` coincident pairs. The author's
   (Q1)+(Q2) imply it.
   * Two remarks on its form. ZE-TC^rat is vacuous at rational zeros. Definition 5.4 needs a finite family `𝓛`.
   * **J5 below shows that this uniform-`Q_0` form is false on Sylvester boxes `n = 2^12..2^17` for every `Q_0`.**
     The per-subspace radii of Lemma 6.6 should be used instead (§4.1).
9. **Withdrawn / corrected.**
   * Withdrawn: "scales compatible from `n >= 5185`" (Prop 5.5). Its premise `rho_A(F) >= 2(js+e)` is false for
     near-twin columns: exactly, `rho_A = 2j + 2 delta`. With the proved Lemma R1, compatibility holds only from
     `n ~ 2.55e5` (`r = n`, `C = 1`; double-precision evidence), asymptotically `n > (48+o(1))w`. No theorem depends
     on it.
   * Corrected: J1's `p`-adic point becomes cheap at `T = 1` inside the window from `n = 2^13` (`lambda = 1`). It is
     not cheap at `ceil(4w)` for genuine widenings.
10. **J1, J2, J3 (R11) are confirmed.** J1 has exact `m(2^13) = 49388999`.

### 2.2 zetc (`zetc.md`, `referee_zetc.md`): survives with corrections

**Setting.** `T = 1` slab boxes in the ternary cube, cut by `|k(S_p)| <= L_p`, at orders `1 <= m <= rho(A)`.
*Degenerate* means lying on a rational hyperplane.

1. **The tangent degree (proved).** `TC(A,m) = {nu_A < m}`. `nu_A` is:
   * subadditive, and additive on products;
   * monotone under inclusion (`A ⊂ A'` gives `nu_A <= nu_(A')`);
   * bounded above by `rho` of projections, and below by zonotopes.

   Along every rational `M` with `nu_A(M) < m`, every `Phi_c` vanishes transversally to order `>= m - nu_A(M)`. For
   `m <= rho(A)`, ZE-TC^cheap implies ZE-TC^rat, which implies ZE-TC^deg. (For `m > rho`, `TC = C^d`.)
2. **One cut (proved).** `nu_S(u) = min(2 sigma(u), 2L + 2t(u), s + L)`, and `rho = s + L`. The slice step needs the
   referee's corrected indicator `prod_Z (1 - k_x)`. ZE-TC^rat holds for one cut, and for disjoint cuts, at every
   `m <= rho`. *J6 gives an independent algebraic proof of the `rho` statement and identifies `J(A)` (§4.2).*
3. **ZE-TC^cheap is false (proved; strengthened by the referee).** It fails:
   * in zonotopal one-cut models (`|S| + L < m <= rho`);
   * in the zonotope-free model `M8` at `m = 8`. Here `nu <= |S| + L = 6` on every column span, while the rational
     hull is 4-dimensional, so it costs `>= 8`;
   * in the LP-supercritical, zonotope-free, overdetermined, twin-free product `M8 x M8` at `m = 16`.

   All failing points are degenerate. R11 §8.1's "ZE-TC holds in `M8` (J4)" is true only in the **deg** form.
4. **Few cuts are zonotopal (computer-assisted, exact in the direction used, independently re-verified).**
   `zeta = rho` holds, hence ZE-TC^deg for every `m <= rho`, in:
   * all two-cut blocks with `|S_1 u S_2| <= 7`;
   * all 843 three-cut types with block size `<= 6` and `L in {1,2}`.

   So an induction on cuts that propagates zonotopes cannot reach the regime.
5. **Products, sections and the Sylvester tower (proved; deg form only).** ZE-TC^deg passes to products and from
   lattice sections. Conditional on J4 (unrefereed, soundly derived), ZE-TC^deg holds for:
   * `Syl(2^k, 2)` at `m = 8`;
   * `Syl(2^k, 3)` at `m = 11`.

   Both are LP-subcritical. **Deg-type statements control neither the height of the normal vectors nor the
   transversal order, so they do not feed the near-strata step** (referee C2).
6. **Small models (rigorous upper bounds on `rho`).**
   * No Hadamard `T = 1` model with `n <= 12` is LP-supercritical with `V != 0`.
   * `Syl(16, L)`: `rho <= 11` at `L = 2` (referee; the author has `<= 15`), `<= 22` at `L = 3`, and the one-cut
     bound covers `L >= 4`.
   * **The smallest LP-supercritical Sylvester model has `n >= 32`.**
   * Products such as `M8 x M8` meet every stated criterion and carry no new information, so test models must also
     be *indecomposable*.
7. **Paley-12, `L = 2`, `m = 8` (exact).**
   * `rho = 8` and zonotope-free.
   * All 330 eight-point sections have an irreducible, indefinite top octic, hence nondegenerate real zeros. So no
     single section certifies anything.
   * That the full system has no nondegenerate real zero is **evidence only**.
8. **Not proved:** ZE-TC in any form for any indecomposable, LP-supercritical, zonotope-free model.

### 2.3 nearstrata (`nearstrata.md`, `referee_nearstrata.md`): proofs survive, framing corrected

**Setting.** `T = 1`, `w = sqrt(8 d ln 4r)`, `H = U^2`, any profile, any order `m'' in [ceil(4w), m_max]`.

1. **Duality (proved).** `d^m kappa(u)` is the best uniform approximation error of `<k',u>^m` on `A` by polynomials of
   degree `< m`.
2. **Effective Lojasiewicz in zonotopal models (proved).** On the uncut ternary cube and on every one-cut box,
   `kappa_H >= (m!/(2d)^m) (dist_inf(u, TC_H)/(n-1))^m`, with `TC_H` an explicit union of rational subspaces. The
   gadget count was exhaustively verified for `s <= 9`. These boxes are zonotopal, so the mechanism never applies in
   the regime.
3. **What qualitative statements can and cannot give (proved).**
   * `Lambda_A(g) > 0` for all `g` iff ZE-TC^braid. That gives a collision Lojasiewicz inequality only
     ineffectively.
   * **No bound that depends only on dimension, degree, height and "zero set in the braid arrangement" can reach the
     needed exponent `d/2`** (Prop 2.7, explicit forms).
4. **Separation (proved, any cluster, any norm).**
   * Each pair has its own gap: `|hat u_x - hat u_y| >= e^(-eps_xy/2)`.
   * `sum eps_xy <= N_2 abar`, so `Gamma(hat u) <= abar`.
   * `dist(hat u, L) >= (1/2) e^(-N_2 abar/(2P(L)))`.
   * Every width-cheap span at `ceil(4w)` has `<= 16` classes outside `<= (m'-1)/2` points.
5. **Assembly (proved implications).** Step 0 plus separation plus **(MIN)**,
   `Lambda_A(abar) > eps_0` with `Lambda_A(g) = min{kappa_H(u) : |u|_inf = 1, Gamma(u) <= g}`, contradict rigidity
   (Thm 4.1). If (MIN) holds for all regime profiles at all `W > W_1`, then `n = O(loglog R)` on all norms (Cor 4.2).
6. **Strengthened (referee).** Theorem 4.1 holds with the weaker, profile-specific
   **(MIN^pair)**: `min{kappa_H(u) : |u|_inf = 1, |u_x - u_y| >= e^(-eps_xy/2) for all pairs} > eps_0`.
   * The empty minimum counts as `+infinity`.
   * MIN^pair, not MIN\*, is the minimal statement relative to the pair data.
   * Its qualitative core is weaker than ZE-TC^braid. *See J7 (§4.3) on the referee's witness.*
7. **Necessary conditions (corrected wording).**
   * MIN\* at **all** levels forces `o_L N_2/P(L) <= d` for every linear subspace `L` of `TC_H`. Nothing is proved
     for (MIN) at a single level, or at the bounded levels needed for growth.
   * **MIN\* is false on the regime Sylvester box `n = 4096` (`m' = 2256`) for `g > 328.4`** (3-column span,
     `o >= 570`, `o N_2/P = 4568 > 4095`). So `n_1 > 4096` in Cor 4.2.
   * (MIN) is false on Sylvester boxes at orders in `(2jk, m_max]`, `k >= 16`. This obstructs the averaged route only.
   * Width-cheap spans meet the MIN\* necessary condition:
     * whatever their order, for `n >= 473,716 / 564,084 / 947,350`;
     * with the proved orders, for `n >= 5185 / 6398 / 10280`;
     * for `r = n, n ln n, n^2` only. These are scans, and they are not monotone in `r`.
8. **Corrected sufficient condition.** A collision Lojasiewicz inequality with exponent `theta < d` suffices if
   `ln(1/gamma) <= (1 - theta/d)(3/4 - o(1)) n ln n`. This equals `(3/8) n ln n` only when `theta <= d/2`.
9. **Not proved:** MIN\*, MIN, MIN^pair or ZE-TC^braid for any regime box.

### 2.4 fake (`fake.md`, `referee_fake.md`): theorems survive, computational reading does not

1. **Thm 2.1, fake directions (proved, corrected constant `theta' = theta/(8 sqrt n)`).** Let the pattern data
   satisfy (Mar), and let `u*` be a nondegenerate real point of `TC_H(A,m)`.
   * For `P -> infinity` there are gap-admissible fakes whose direction is within `theta e^(-W/4)` of `u*`, with
     `kappa_c(hat u) <= kappa_c(u*) + m theta e^(-W/4)/sqrt d`. These fakes are **not** barrier-scale:
     `W/(n log n)` is unbounded.
   * At the barrier scale `W = O((a/rho_1 + 1) r log n)` the same holds if `|lambda.u*| >= n^(-a||lambda||_1)` for
     all relevant `lambda` with `||lambda||_1 <= r/rho_1 + 1`. That Diophantine condition (2.2) is needed for all such
     `lambda`, not only for pairs.
   * So a nondegenerate real `TC_H` point of a (Mar) box makes ZE-TC_eps fail for that box, with gap-admissible
     witnesses. At the barrier scale this needs (2.2).
   * "Neither easier nor harder than ZE-TC" is overstated: the two directions need different hypotheses.
2. **Thm 2.2, lifting (proved, modulo a cited uniform sublevel lemma).** In underdetermined injective systems, a
   regular, non-flat real common zero lifts to fakes lying exactly on `Z_P`. Gap-admissible fakes on `Z_P` exist for:
   * sdc7: a certified zero with `sep >= 0.2496`, at 99.8% of the true maximum `1/4`;
   * sdc11: `dim V(A,8) = 5`, a certified zero with `sep >= 0.1062`, at 64% of the maximum `1/6` (referee).

   Both models are `T = 1`, `L = 2`, `m = 8`, for every base point. They are LP-subcritical and outside the regime.
3. **Exact negatives (referee).**
   * sdc7, pal7 and sdc11 at `L = 2`, `m = 7` contain braid zonotopes of length `7 = m`. Their "collapse onto
     collision strata" is a single character, and carries no regime information.
   * pal7 at `L = 3` is `M8'` only at base `inf`, where the twins are invisible.
   * The landscape optimiser misreported minima: it found positive minima where the certified value is 0, and those
     rows were omitted from the table.
4. **The only faithful overdetermined fake model computed is pal11** (`L = 2`, `m = 8`, `zeta <= 22/3`,
   `dim V <= 436` mod `p`). For it there is only unreliable floating-point evidence.
5. **Not proved:** whether regime-size fake `T = 1` boxes have nondegenerate real tangent-cone points. No computation
   so far bears on it.

## 3. Status table: ZE-TC / ZE-TC_eps / separation / ZE at T = 1

All entries are at `T = 1`, `m' = ceil(4w)` unless stated otherwise. "Regime box" means the slab box of an admissible
profile with `n >= n_R(r)` (`m_max >= 4w`). Status values: **proved**, **false**, **open**, **evidence**.

### 3.1 Inputs of the reduction (all proved)

| statement | status | source |
| --- | --- | --- |
| Slab box, Siegel count, `m_max >= 4w` from `d >= 316 ln 4r` | **proved** (thresholds exact) | t1rigour §2 |
| Rigidity at every `m' in [ceil(4w), m_max]`, `W > W_1` | **proved** | t1rigour Thm 3.1 |
| Step 0: `kappa(hat u) <= eps_0` | **proved** | t1rigour Thm 4.2; nearstrata Lemma 1.1 |
| Transversal order `>= m' - rho_A(F)`; Lemma R1 `rho_A(F) >= 2 min(j, n - max class)` | **proved** | t1rigour §5; referee |
| Cheap families: `<= 4` long directions, `<= 16` classes; cheap spans in the braid arrangement at `ceil(4w)` | **proved** | t1rigour Thm 6.3, Cor 6.4 |
| **Arithmetic separation**: `dist(hat u, L) >= (1/2) e^(-N_2 abar/(2P(L)))` and per-pair gaps, any cluster, any norm | **proved** | t1rigour Lemma 6.6; nearstrata Lemmas 3.1-3.3 |
| Scale compatibility of separation radii with transversal ceilings (cheap spans) | per-subspace: necessary check passes for `n >= 5185` (`r = n`; proved orders; scan). Uniform `Q_0`: **fails** for `n = 2^12..2^17` (J5) and is not established below `~2.55e5` | nearstrata Cor 3.7, Prop 5.3; referee_t1rigour F1; J5 |

### 3.2 Qualitative tangent-cone statements

| statement | status | source |
| --- | --- | --- |
| ZE-TC^cheap (R11 §8.2: real zeros lie in cheap spans) | **false**: in the regime-type Sylvester `2^12`+1 profile at `ceil(4w)`, in `M8` (`m = 8`) and in `M8 x M8` (`m = 16`); all failing points degenerate | t1rigour Thm 6.11(b); referee_zetc C1 |
| ZE-TC^rat (real zeros lie in rational subspaces `M` with `nu_A(M) < m`) | **proved** for one cut and disjoint cuts; vacuous at rational zeros; **open** in the regime. Not sufficient by itself: a collision-free rational stratum is fatal | zetc Cor 2.2; referee_t1rigour F4; nearstrata §2.4 |
| ZE-TC^deg (real zeros on some rational hyperplane) | **holds** in zonotopal models, all small two/three-cut blocks, `Syl(2^k, 2/3)` at `m = 8/11` (conditional on J4, LP-subcritical). **Insufficient**: no height or order control | zetc §3-4; referee_zetc C2 |
| ZE-TC^braid (`Lambda_A(g) > 0` for all `g`) | **open** in the regime. Cheap spans satisfy it at `ceil(4w)`. False at Sylvester orders `(2jk, m_max]`, `k >= 16` | nearstrata Prop 2.6; t1rigour Prop 6.5 |
| ZE-TC^coll (every real zero has `>= N_2/Q_0` coincidences) | **open**; consistent with everything proved | referee_t1rigour §6.2 |
| Tangent cone of zonotope-free LP-supercritical **indecomposable** models | **no example computed**. Sylvester needs `n >= 32`, beyond exact computation | zetc §5 |
| Paley-12 `L = 2`, `m = 8`, full system | **evidence** (numerical) for no nondegenerate real zero; every section has one (exact) | zetc §5.5; referee_zetc |
| pal11 (faithful, overdetermined fake model) | **evidence only**, from an unreliable optimiser | referee_fake |

### 3.3 Quantitative statements and the final hypothesis

| statement | status | source |
| --- | --- | --- |
| Effective ZE-TC_eps (Lojasiewicz with explicit constants) | **proved** on the uncut cube and every one-cut box (zonotopal only) | nearstrata Thms 2.4-2.5 |
| Effective Lojasiewicz from dimension, degree, height and qualitative zero set alone | **impossible** (cannot reach exponent `d/2`) | nearstrata Prop 2.7 |
| ZE-TC^coll_eps with uniform `Q_0` (referee's single target; author's (Q1)+(Q2)) | **false** on Sylvester `n = 2^12..2^17` for every `Q_0`, at every `abar >= 5.42` (J5, unrefereed). Asymptotically open; needs `n >~ 48w` | §4.1 |
| ZE-TC^coll_sigma, per-subspace radii `sigma_pi = (1/2) e^(-N_2 abar/(2P(pi)))` (J5) | **false** at `n = 4096` for `abar >= 14.9` (`W >= 2.7 W_min`); not refuted by coordinate spans for `n = 2^13, 2^14`; **open** | §4.1 |
| MIN\* (`Lambda_A(g) > (d/(m'+1)) e^(-dg/2)` for all `g >= 1.5 ln n - 7`) | **false** at `n = 4096` for `g > 328.4`; **open** for large `n` | referee_nearstrata §2.5 |
| (MIN) at the single level `abar` | **open** | nearstrata Thm 4.1 |
| **MIN^pair** (gap region of the actual pair data) | **open**; the weakest form that suffices | referee_nearstrata §4 |
| Implications | ZE-TC^coll_eps(`Q_0`) ⇒ ZE-TC^coll_sigma ⇒ (MIN) ⇒ MIN^pair; MIN\* ⇒ (MIN) (J7) | §4.3 |

### 3.4 ZE and the fake test

| statement | status | source |
| --- | --- | --- |
| **ZE at `T = 1`** (no nondegenerate `delta` near `e` in `Z_P`) | **open**. Reduced (proved) to MIN^pair for the cluster's own box, which suffices on all norms | §2.3 |
| Fake test: nondegenerate real `TC_H` points of (Mar) fake boxes give gap-admissible fakes (barrier scale needs (2.2)) | **proved** (fake Thm 2.1, corrected) | §2.4 |
| Fakes exactly on `Z_P` | **exist** in the underdetermined sdc7 and sdc11 models (`T = 1`, `L = 2`, `m = 8`; LP-subcritical) | fake Cor 2.3, 3.3; referee |
| Overdetermination removes the nondegenerate points of fake boxes | **no evidence either way**. The claimed evidence was zonotopal, optimiser-based or refuted | referee_fake |

**Overall.** Separation is **proved**. The cheap-span form of ZE-TC is **refuted**. The qualitative core (braid or
coll) and every sufficient quantitative form (MIN^pair and stronger) are **open** for regime boxes, with **no
evidence for or against** in the regime. ZE at `T = 1` is **open**. It is reduced by proved implications to one
inequality about the profile's own box, and that inequality is exactly the fake test (fake Thm 2.1).

## 4. Judge's own findings (unrefereed)

Three findings: J5, J6 and J7 (R11 used J1-J4). J5 and J6 are backed by exact computation (scripts in
`judge_checks/`). J7 is logic only. **None of them has been refereed.**

### 4.1 J5: the uniform-`Q_0` collision hypothesis is false on Sylvester regime boxes; use per-subspace radii

The referee of t1rigour proposed one hypothesis, ZE-TC^coll_eps`(A, ceil(4w); Q_0)`, as the target. It is the
author's (Q1)+(Q2) with the same single `Q_0`. Write `sigma = (1/2) e^(-Q_0 abar/2)`. The hypothesis says:
`kappa(hat u) > eps_0` whenever `|hat u|_inf = 1` and `dist(hat u, Delta_pi) >= sigma` for every partition `pi` with
`P(pi) >= N_2/Q_0`. Here `Delta_pi = {u : u_x = u_y for x ~_pi y, u_base = 0}`.

**Proposition J5.** Let `n = 2^k` with `12 <= k <= 17`. Take the Sylvester profile:
* points `x in F_2^k`;
* sign columns `(-1)^(a.x)`, `a != 0`, so `r = d`;
* any primes, any base point, `C <= 1`.

Then `m' = ceil(4w) <= m_max` (exact). For **every** `Q_0 >= 1` the hypothesis fails at every level
`abar = W/(2d) >= 5.42`, which includes all `W >= 3n ln n - 14n`.

*Proof.*
1. **Setup.** Let `S_i = {x : x_i = 0}`. On `A`, `<k, beta> = 2k(S_i)`, so `|k(S_i)| <= j`. The family
   `F_b = (1_(S_1), ..., 1_(S_b))` has image in `{-j..j}^b`, so `rho_A(F_b) <= 2jb`. For `b <= 4` the computed values
   give `8j <= m' - 1`. Lemma 5.0 of nearstrata therefore gives transversal order `o_b >= m' - 2jb` and
   `kappa(hat(y + eta v)) <= 4^(m') eta^(o_b)`, for `y in L_b = span F_b` with `|y|_inf = 1`, `|v|_inf <= 1` and
   `eta <= 1/2`.
2. **The test point `y`.** Take `y = sum_(i<=b) 2^(i-1)(1_(S_i) - 1_(S_i)(base))/(2^b - 1)`, normalised. Then
   `y in L_b`, and it takes distinct values on the `2^b` classes of `(x_1, ..., x_b)`, with gaps `>= 1/(2^b - 1)`.
   Also `P(L_b) = 2^b binom(n/2^b, 2)`; write `Q_b = N_2/P(L_b)`.
3. **Case A: `Q_0 < Q_4`** (`Q_4 ~ 16.06`).
   * Every partition with `P(pi) >= N_2/Q_0 > P(L_4)` must merge two classes of `L_4`. So `dist(y, Delta_pi) >= 1/30`.
   * Also `1/30 >= sigma` once `Q_0 abar >= 2 ln 15`, in particular for `abar >= 5.42`.
   * But `kappa(y) = 0`.
4. **Case B: `Q_0 >= Q_4`.**
   * Take `u = y + eta v`, where `y in L_1` and `v` takes distinct values `rank/(n/2)` inside each class, with `0` at
     the base point.
   * A partition with `P(pi) >= N_2/Q_0` either merges the two classes, giving distance `>= (1-eta)/2`, or refines
     them. In the second case it has a block of size `>= 1 + (n-1)/Q_0`, giving distance
     `>= eta(n-1)/(Q_0 n)`.
   * With `eta = 2 sigma Q_0 n/(n-1)`, the point `hat u` is `sigma`-far from the whole family.
   * Step 1 gives `kappa <= 4^(m') eta^(o_1)`, with `o_1 >= m' - 2j`. This is `<= eps_0` as soon as
     `abar >= g*(Q_0) = 2[m' ln 4 + o_1 ln(Q_0 n/(n-1)) - ln(d/(m'+1))]/(o_1 Q_0 - d)`.
   * `g*` decreases for `Q_0 >= e`: the analytic sign of `dg*/dQ_0` agrees with a grid check. So `g*` is at most
     `g*(Q_4)`, which is `0.68, 0.73, 0.80, 0.94, 1.24, 2.18` for `k = 12..17`.

QED. Margins: `judge_checks/uniformQ0_sylvester_out.txt` (exact `m_max` test, Fraction collision counts, 60-digit
logs). For `k >= 18` we get `o_1 Q_4 < d`, and this mechanism stops.

**Consequences.**
* The single target recommended by referee_t1rigour, and the author's (Q1)+(Q2), force `n_1 > 2^17` in Theorem 6.8.
  They do so at every near-extremal `W`, not only at absurd `W`.
* This turns the referee's withdrawal of Prop 5.5 ("compatibility not established below `~2.55e5`") into an actual
  failure on admissible boxes.
* **The fix is free.** Lemma 6.6 already supplies **per-subspace radii** `sigma_pi = (1/2) e^(-N_2 abar/(2P(pi)))`,
  and the first part of Theorem 6.8 accepts any family with these radii. Call the corresponding statement
  **ZE-TC^coll_sigma**: `kappa(hat u) > eps_0` whenever `dist(hat u, Delta_pi) >= sigma_pi` for every partition `pi`
  with `P(pi) >= 1`. This is a finite family of subspaces.
* With per-subspace radii the same test points (`judge_checks/persubspace_sylvester_out.txt`) give:
  * failure at `n = 4096` only, for `abar >= 14.90` (3-column span; `W >= 1.22e5 ~ 2.7 W_min`) and for
    `abar >= 17.93` (2-column span; `o_2 Q_2 = 4531 > 4095`). The referee found only the 3-column span, and only for
    MIN\* at `g > 328`.
  * no failure from coordinate spans at `n = 8192, 16384`.

### 4.2 J6: algebraic structure of the tangent cone (apolarity, Lefschetz, locality)

**(a) Fat-point lower bound (proved here; a short argument).** For every integer direction `v` with
`h = h_A(v) < m`, the polynomial `prod_t (<k',v> - t)` over the `h + 1` values vanishes on `A`. So
`<k',v>^(h+1) in J(A)`. Let `J_w` be the ideal these generate (it contains `k_x^3`).
* By compute R2, `span{Phi_c : c in V(A,m)} = J(A)_m^perp ⊂ (J_w)_m^perp`.
* By apolarity, `(J_w)_m^perp` is the space of degree-`m` forms `Phi` with `(v.∂)^(h_A(v)+1) Phi = 0` for all `v`.
  That is: `Phi(tv + w)` has `t`-degree `<= h_A(v)`, i.e. `Phi` vanishes to order `>= m - h_A(v)` at `[v]`.

Hence

```text
TC(A, m)  ⊇  Bs_R( L_m ),     L_m = { degree-m forms with multiplicity >= m - h_A(v) at every [v] in P^(d-1) } .
```

A linear span of independent points `[v_1], ..., [v_s]` lies in `Bs(L_m)` when `sum_i (m - mult_i) < m`, i.e.
`sum h_A(v_i) < m`. **The cheap spans are exactly the linear base cycles of this fat-point system.** Non-linear base
components of fat-point systems (rational normal curves, secant varieties, Cremona-special loci) are classical. By
the display above, any of them would also lie in `TC`. This gives a source of tangent-cone points beyond cheap spans
that is computable from widths alone.

**(b) One cut is the Lefschetz property (exact checks; proof uses a quoted theorem).** For a one-cut box
`A_S(L) ⊂ {-1,0,1}^d`, `J(A_S(L)) = (k_x^3 (all x), sigma_S^(2L+1))` in **all** degrees.

*Proof.*
1. The inclusion `⊇` is (a).
2. `k[x_S]/(x^3)` has the strong Lefschetz property with Lefschetz element `sigma_S` in characteristic 0. This is
   Stanley 1980 for monomial complete intersections; it is quoted, not reproved.
3. Hence `S/(x^3, sigma_S^(2L+1))` has Hilbert function `max(h_i - h_(i-2L-1), 0)`, where `h` is the coefficient
   sequence of `(1+t+t^2)^s`, tensored with `(1+t+t^2)^(d-s)`.
4. Its total is `3^(d-s) sum_(|i-s|<=L) h_i = |A_S(L)|`, which is also `dim S/J(A)`. So the two ideals are equal.

QED. `judge_checks/onecut_slp_out.txt` gives exact agreement of `H_A`, `H(S/J_0)` and the Lefschetz formula on 6
boxes, all degrees (mod-`p` sandwich, see below). This gives zetc's `rho = s + L` in one line. It also makes
Theorem 2.1 a statement about which powers `l_u^m` lie in an explicit ideal.

**(c) Locality in the faithful model `M8` (exact).** In `M8` (`Syl(8)`, `L = 2`, `m = 8`; `judge_checks/jw_test_M8_out.txt`,
`jw_M8_exact6_out.txt`, `jproj3_test_M8_out.txt`):
* `J(A) != J_w`. Exactly over `Q`, `dim (J_w)_6 = 616 < 637 = dim J(A)_6`. The 21 missing degree-6 generators are not
  width powers.
* Adding pullbacks `pi_F^* J(pi_F(A))` for all **pairs** of columns adds nothing. This is mod `p`, i.e. evidence.
* Adding them for all **triples** of columns gives `J(A)` **exactly in every degree `<= 8`**: dims
  `[0,0,0,7,49,203,637,1541,2961]`.
* In the one-cut boxes, two two-cut boxes and `Syl(8), L = 1`, `J(A) = J_w` holds in all degrees `<= m`
  (`jw_test_small_out.txt`, `jw_test_M8L1_out.txt`).

*How the equality is certified.* Pulled-back vanishing polynomials of a projection vanish on `A`, and their leading
forms are the pulled-back leading forms. So `J_w + sum_F pi_F^* J_F ⊂ J(A)` over `Q`. The mod-`p` dimension of the
generated ideal is `<=` its rational dimension, and the mod-`p` rank of the evaluation matrix is `<=` its rational
rank. So equality of the cumulative counts mod `p` proves equality over `Q`. The 2D/3D nullspaces are computed over
`Z`.

*Reading.*
* In the only faithful model computed, the leading-form ideal is **local** up to degree `m`: it is generated by
  pullbacks from boxes cut by at most 3 columns.
* If the same holds in the regime with a bounded number of columns, the tangent forms are
  `J(A)_m^perp = ∩_F (pi_F^* J_F)_m^perp`. That is an intersection of local inverse systems, each a few-cut problem
  of the kind zetc already solves (one cut exactly, few cuts zonotopally).
* This is a hypothesis only. The regime's `m' ~ 8j` lets far more columns interact than `M8`'s `m = 8 = 2 * 2L`.

### 4.3 J7: the hierarchy of final hypotheses, and the MIN^pair witness

**(a) Hierarchy (proved here; one line each).** Fix a cluster's data, and let
`G = {u : |u_x - u_y| >= e^(-eps_xy/2)}` be the gap region.
* `G ⊂ {Gamma <= abar}`, because `Gamma(u) <= (1/N_2) sum eps_xy <= abar` on `G`.
* `{Gamma <= abar} ⊂ R_per`, the per-subspace far region (nearstrata Lemma 3.3).
* `R_per ⊂ R_unif(Q_0)`, because for `P(pi) >= N_2/Q_0` the per-subspace radius is the larger one.

Hence

```text
ZE-TC^coll_eps(Q_0)  =>  ZE-TC^coll_sigma  =>  (MIN) at level abar  =>  MIN^pair ,      MIN* => (MIN) for abar >= 1.5 ln n - 7 .
```

All of these suffice for Theorem 4.1 / 6.8. J5 refutes the first at `n <= 2^17`. The other three are open.

**(b) The referee's MIN^pair witness has an empty gap region.** Its data satisfy `(n-1) e^(-max eps/2) > 2`. So the
pair inequalities **alone** forbid any configuration with that profile, those primes and that `W`: no actual cluster,
and no gap-admissible fake, carries these data. On data carried by an actual cluster, `G` contains `hat u` and is
non-empty. So the example shows that MIN^pair is weaker than (MIN) **as a hypothesis**. Whether it is weaker on
non-empty gap regions is **open**, and that is the case that matters. For (Mar) fake data, fake Thm 2.1 shows that
MIN^pair fails as soon as `TC_H` contains a nondegenerate point satisfying (2.2): the fakes it produces lie in `G`.

## 5. What it means for route 5 and the growth rate

**Growth rate: unchanged.** `M <= (2/3 + eps) log R/loglog R` is still the best bound for actual circles. Nothing
this round improves it, unconditionally or from a standard conjecture. The uniform bound is not touched.

**Conditional theorems on record after round 12.**

| hypothesis | conclusion | implication | hypothesis status |
| --- | --- | --- | --- |
| ZE(`n' >= n_1`), `T = n` box | `n = O(loglog R)`, squarefree norms | proved ([D] Thm 5.7) | open |
| ZE-col\*(univ) | `M = O_C((loglog R)^2)`, all norms | proved (R11 subcluster Thm 5.4) | open |
| **MIN^pair** for the cluster's own `T = 1` box at `ceil(4w)`, all regime profiles, `n >= n_1`, `W > W_1` | `n = O(loglog R)`, **all norms**, no rank lemma, no nondegenerate subcluster | proved (nearstrata Thm 4.1 + referee; t1rigour Thm 6.8) | open; weakest known sufficient form |
| (MIN), ZE-TC^coll_sigma | same | proved | open (`n_1 > 4096`) |
| ZE-TC^coll_eps with uniform `Q_0` | same | proved | **false for `n <= 2^17`** (J5); asymptotically open |
| SL, SL\*, ZE-braid, ZE-TC^cheap | — | — | **false** |

**Meaning for route 5.**
1. **The reduction is finished, and its residue is one inequality.**
   * Round 12 proves every arithmetic ingredient of the `T = 1` route: rigidity, Step 0, separation and the cheap
     structure.
   * It does so on **all norms**, with no rank or squarefree input.
   * What is left is a single inequality about the profile's own slab box: `kappa > eps_0` on the gap region, or on
     the per-subspace far region. That is an improvement on R11, where ZE at `T = 1` was tied to the squarefree
     rank lemma.
2. **The residue is genuinely quantitative and genuinely global.**
   * Qualitative ZE-TC (braid or coll) is not enough, because general effective Lojasiewicz theorems cannot supply the
     size (Prop 2.7).
   * Model computations cannot see the regime:
     * the first LP-supercritical Sylvester model has `n >= 32`;
     * products are decomposable;
     * few-cut boxes are zonotopal;
     * every "overdetermined" fake model computed this round except pal11 is zonotopal at the tested order.
3. **The residue is the barrier itself.** By fake Thm 2.1, MIN^pair for the (Mar) fake pattern boxes is equivalent
   (up to the Diophantine condition (2.2)) to the absence of nondegenerate real tangent-cone points there.
   * If such points exist, gap-admissible fakes exist whose directions satisfy Step 0, and route 5 is capped as
     formulated.
   * If they do not exist, with the needed size, route 5 gives `n = O(loglog R)`.
   * Either outcome is decisive. Neither is in sight.
4. **Fake-violated input (rubric b).** This is unchanged from R10-R11.
   * The input fakes violate is rigidity, the exact relation (4.2) at the actual point. Fakes satisfy it only on
     `Z_P`, and generic fakes are off `Z_P`.
   * The separation lemma uses only pair inequalities, which fakes satisfy.
   * The open step is exactly "no fake direction in the near tangent cone".
5. **Finite-`n` obstructions are real but not asymptotic.**
   * MIN\* fails at `n = 4096`.
   * The uniform-`Q_0` form fails up to `2^17` (J5).
   * The per-subspace form fails at `4096` from `W ~ 2.7 W_min`.

   So any proof must be asymptotic, and must use per-subspace radii, or better the actual pair data.

## 6. Scores

Rubric of `round10/README.md` §1: (a) correctness after refereeing; (b) fake-violated input; (c) movement of growth
or uniformity. The overall score is a judgement and weighs (b) and (c) above (a).

| write-up | (a) | (b) | (c) | **overall** | reason |
| --- | --- | --- | --- | --- | --- |
| **t1rigour** | 8 | 3 | 3 | **5** | Every theorem the reduction uses survives independent re-derivation. It turns R11's J1-J3 and §5.3 sketches into proofs: closed-form `T = 1` box, `W_1`, Step 0, the cheap-family structure theorem, and the separation lemma on all norms. It refutes cheap-span ZE-TC on an admissible, zonotope-free, LP-supercritical regime profile at `ceil(4w)`. Deductions: the scale-compatibility claim (Prop 5.5) rested on a false premise; the `n = 7` toy is zonotopal; J1 at `T = 1` was mis-costed; and (Q1)+(Q2) was badly formed (infinite family; uniform `Q_0`, which J5 shows is false at `n <= 2^17`). |
| **nearstrata** | 7 | 3 | 3 | **4** | All proofs survive. Duality, effective bounds on one-cut boxes, the impossibility result Prop 2.7, the separation lemma, and the clean assembly theorem (MIN implies growth on all norms) are useful. The necessary conditions (Props 5.1-5.2) are new. Deductions: the "minimal statement" claim was false (MIN^pair is weaker); the necessary condition was mis-attributed to single-level (MIN); the CL budget was misstated; constants were swapped in Prop 5.3; the thresholds are not monotone in `r`. MIN\* is false at `n = 4096`, which the write-up missed. |
| **zetc** | 7 | 2 | 2 | **4** | The tangent degree `nu_A` is the right object, with transversal order `m - nu_A(M)`. The one-cut formula is exact. The few-cut zonotopality theorems (computer-assisted) and the small-model bounds (Sylvester needs `n >= 32`) are informative. Deductions: the headline refutation used an unfaithful witness (the referee supplied faithful ones); the deg-type results (Sylvester tower, products) cannot feed the quantitative step; a slice-lemma degree error; "every order" should read `m <= rho`; and the Paley claims were overstated. |
| **fake** | 4 | 3 | 2 | **3** | Thm 2.1 is the right formal reduction of the fake test to the tangent cone of fake boxes. Thm 2.2 and its certified underdetermined realisations are correct. But the computational half fails: three of four "overdetermined" models are zonotopal; the claimed cone shrinkage is refuted by a certified zero at 64% of the maximum; the optimiser misreports minima, and those rows were omitted; and the normalisation was wrong. |

No write-up moves the growth rate. Partial credit in (c) is for:
* the proved all-norms reduction (t1rigour, nearstrata);
* the refutation of the R11 target (t1rigour, zetc);
* the proved reduction of the fake test (fake).

## 7. Single most promising next step

**Recommended: determine the local structure of the leading-form ideal `J(A)` of regime `T = 1` boxes in degrees
`<= m' = ceil(4w)`, then use it to decide the strata half of ZE-TC^coll_sigma.**

### 7.1 Statement to decide

> **Locality `Loc(b_0)`.** For every admissible profile with `n >= n_R(r)` and every `i <= ceil(4w)`,
>
> ```text
> J(A)_i  =  sum over families F of <= b_0 integer directions with rho(pi_F(A)) < ceil(4w)  of  ( pi_F^* J(pi_F(A)) )_i ,
> ```
>
> for a constant `b_0` independent of `n`.

By apolarity (J6a), `Loc(b_0)` gives the tangent forms as an intersection of local inverse systems:

```text
span{Phi_c} = ∩_F (pi_F^* J_F)^perp_(m') .
```

So `TC(A, m')` becomes the real base locus of an explicit linear system, built from boxes cut by at most `b_0`
directions.

### 7.2 Why this step

* **It is the only structural handle on the regime tangent cone left standing.**
  * Model computations cannot reach the regime (§5 item 2).
  * General effective Lojasiewicz theorems provably cannot give the size (nearstrata Prop 2.7).
  * The arithmetic is spent: separation is proved, and only the box's algebra remains.
* **It holds in every model computed** (J6):
  * exactly for one cut, by the Lefschetz property;
  * exactly for the two-cut boxes with `d <= 6` and `Syl(8), L = 1`, with `b_0 = 1`, i.e. width powers;
  * exactly for the faithful model `M8`, with `b_0 = 3`, and not with `b_0 = 2` (evidence for the last part).
* **It would settle the rational strata.**
  * Under `Loc(b_0)` every rational subspace `M` in `TC` is detected by projection costs of at most `b_0`
    directions.
  * The necessary per-stratum condition `o_M N_2/P(M) < d` (nearstrata Prop 5.1; J5) would then reduce to the
    width / Lemma R1 combinatorics of t1rigour §6. That combinatorics already shows that cheap spans have
    `<= 16` classes.
  * This closes the main open point of t1rigour (Q1) and of nearstrata §7.2(b).
* **It would put the irrational part into a classical framework.**
  * The remaining question becomes: does the base locus of an intersection of few-direction fat-subspace systems
    contain collision-poor real points?
  * Base loci of fat-point systems are governed by linear cycles (here: the cheap spans), Cremona-special loci and
    rational normal curves through the points. There are degeneration tools for them.
* **Either outcome is informative.**
  * A non-local low-degree relation in `J(A)` would be a global constraint coupling many columns, which is
    exactly where nondegenerate tangent-cone points of fake boxes could come from.
  * A proof of locality reduces route 5 to a fat-point base-locus statement plus its quantitative form.

### 7.3 First concrete tasks

1. **Referee J5-J7** (§4). J5 needs an independent check of the two-case argument and of the margins. J6(b) needs
   the quoted Lefschetz theorem checked (strong Lefschetz property of `k[x_1..x_s]/(x_i^3)` with element
   `sum x_i`, characteristic 0).
2. **Locality on exact models, using symmetry to make them feasible.** For each model, find the smallest `b_0`.
   * `M8'` (`L = 3`, `m = 11`): sparse elimination, or the `AGL(3,2)` isotypic decomposition.
   * Paley-12 `L = 2` (`m = 8`): `Z_11 x {+-1}` isotypic components, as in zetc.
   * `Syl(16, 2)` (`rho <= 11`): Singer-cycle components.
   * All two- and three-cut blocks with `d <= 8`.

   Track how `b_0` depends on `m/(2L)`. In the regime, `m'/(2j) ~ 4`.
3. **Projection-cost lower bound (theory).** Show that a family `F` contributes generators in degree `<= m'` only if
   `rho(pi_F(A)) < m'`. Then prove `rho(pi_F(A)) >= 2j(|F| - O(1))` for families of balanced columns, extending
   Lemma R1 (`>= 2j`) by a zonotope or fibration argument inside `pi_F(A)`. Together these would bound `b_0` by
   `m'/(2j) + O(1)`.
4. **Only then: the per-stratum check and the quantitative form.**
   * Verify `o_M N_2/P(M) < (1 - eps) d` for all rational strata.
   * State MIN^pair, or ZE-TC^coll_sigma, for the local system.
   * Test it on the fake pattern boxes (sdc, pal) through fake Thm 2.1. That is the barrier test.

### 7.4 Reading of the outcomes

* **`Loc(b_0)` holds and the strata check passes.** Route 5 reduces to a quantitative base-locus statement for an
  explicit local linear system. Prove that statement first for the fake pattern boxes, which have large symmetry
  groups.
* **`Loc` fails, i.e. `J(A)` has global low-degree generators.** Construct real tangent-cone points from them on fake
  pattern boxes.
  * If such a point is nondegenerate and satisfies (2.2), fake Thm 2.1 caps route 5 as formulated.
  * Record this as an addition to the barrier map.

**Not recommended:**
* uniform-`Q_0` collision hypotheses (J5);
* MIN\* at unbounded levels (false at `n = 4096`);
* deg-type ZE-TC results (no height or order control);
* products or few-cut models as evidence (decomposable or zonotopal);
* floating-point landscape searches for overdetermination (unreliable, per referee_fake);
* further delta-free or braid-only reformulations (R11).

## 8. File index

All paths are relative to `round12/`. Each script writes `<name>_out.txt` next to itself.

| write-up | notes | report | scripts (author / referee) |
| --- | --- | --- | --- |
| `t1rigour.md` | `t1rigour_notes.md` | `referee_t1rigour.md` | `t1rigour_checks/` / `referee_t1rigour_checks/` |
| `zetc.md` | `notes.md` | `referee_zetc.md` | `zetc_checks/` / `referee_zetc_checks/` (notes `referee_zetc_notes.md`) |
| `nearstrata.md` | `nearstrata_notes.md` | `referee_nearstrata.md` | `nearstrata_checks/` / `referee_nearstrata_checks/` |
| `fake.md` | `fake_notes.md` | `referee_fake.md` | `fake_checks/` / `referee_fake_checks/` (notes `referee_fake_notes.md`) |

**Judge's checks** (`judge_checks/`; working notes in `judge_notes.md`). All are unrefereed.

| script | what | status |
| --- | --- | --- |
| `uniformQ0_sylvester.py` | J5: Sylvester `k = 12..20`. Exact R1 test (`8 binom(m'+d-1,d) <= 3^d`, big integers, `k <= 17`), Fraction collision counts, failure levels of uniform-`Q_0` ZE-TC^coll_eps (cases A and B), monotonicity grid of `g*` | exact integers; 60-digit logs |
| `persubspace_sylvester.py` | J5: per-subspace version, coordinate spans `b = 1..4`, `k = 12..14` | exact counts; 50-digit logs |
| `jw_test.py` | J6: `J(A)` vs the width ideal `J_w` (mod-`p` sandwich) on one-cut, two-cut, `Syl(8)` `L = 1, 2` | exact via sandwich where equal |
| `jw_M8_p2.py`, `jw_M8_exact6.py` | J6c: the `M8` degree-6 gap; second prime, then exact over `Z`: `616 < 637` | exact |
| `jproj_test.py`, `jproj3_test.py` | J6c: `J_w` plus 2-, 3- and 4-column projection ideals for `M8`. 3 columns give `J(A)` exactly in degrees `<= 8` | exact (sandwich); the 2-column negative is mod `p` |
| `onecut_slp.py` | J6b: `J(A_S(L)) = (x^3, sigma^(2L+1))`, Hilbert function = Lefschetz truncation, 6 boxes, all degrees | exact via sandwich |

**Caveats.**
* J5-J7 are unrefereed.
* J6(b) quotes the strong Lefschetz property of monomial complete intersections; it is not reproved here.
* "Locality" (J6c, §7) is a hypothesis supported by one faithful model only.
* The repository is read-only. Corrections to apply when round 12 is merged:
  * `round11/README.md` §8.1-8.2: the cheap-span form of ZE-TC is false. "ZE-TC holds in `M8`" holds only in the
    deg form.
  * `round11/README.md` §5.3 J3: the thresholds and `W_1` are now proved (t1rigour §2-3).
  * `t1rigour.md` Prop 5.5 is withdrawn.
  * `t1rigour.md` Thm 6.8 (Q1)+(Q2), and the referee's ZE-TC^coll_eps, are false at `n <= 2^17` with uniform `Q_0`
    (J5). Use per-subspace radii.
  * `nearstrata.md` §4.3: the minimality claim is false; MIN^pair is the weakest form.
