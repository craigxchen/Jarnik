# Round 13 — Judge README (route 5: locality Loc(b0))

## 0. Bottom line

**No growth improvement and no uniform bound were proved this round**, either unconditionally or from a standard
conjecture. `M <= (2/3 + eps) log R/loglog R` (`two_thirds/proof.md`) is still the best bound for actual circles.

**Route 5 is stuck.** It is not capped: no refutation of MIN^pair, and no nondegenerate tangent-cone point of a
regime fake box, is known. But three rounds (11-13) have not moved its residual hypothesis, and this round removed
the last structural handle the record proposed for it:
* **Loc(b0)** (README-12 §7.1) is **neither proved nor refuted for any regime profile**.
  * Projection locality onto columns is **refuted exactly** (Paley-12, degree 7, every `b <= 5`).
  * The README's restriction `rho(pi_F A) < ceil(4w)` rests on a contribution criterion that is **false** (an
    index-3 image with `rho = 10` contributes in degree 7). The correct criterion is linkage: `alpha_F <= e <= s_F`.
  * In toy models, the smallest generating families found get larger as the congruences on the image get more
    complex: 1 direction (cut boxes, M8'), 3 (M8, where 2 provably fail for columns), 6 (complete Fano block
    `Y_3(2)`), 7 (`AG(3,2)`, `j = 1`). For the last two these are upper bounds. In Paley-12, 50 invariant degree-7
    dimensions are explained by no family tried.
  * No model is LP-supercritical, and none has the regime ratio `m'/j ~ 8`.
* **Even a proof of Loc would not decide the barrier.** In the fake model sdc7, `Loc_8(1)` holds exactly, and the
  local (width) system still has certified nondegenerate, well-separated real base points.
* **MIN^pair is open**, with no evidence either way at regime sizes. The one new sufficient condition, Theorem 6.4
  of fatpoint (BLOCK implies MIN^pair for Sylvester-type profiles), has two limits:
  * J9: its headline application is already an unconditional theorem (walsh-extraction Theorem B);
  * J8: its smallest candidate block (`t0 = 5`) is **false** for every `j0`.

What survives refereeing is a body of correct but local algebra:
* the linkage criterion `alpha_F`;
* locality of parity sections, through circuits;
* the string theorem for atom-constant cut boxes;
* a self-contained proof of the strong Lefschetz property, giving `Loc(1)` for disjoint signed cuts;
* the unconditional inclusion `Bs_R(local system) ⊂ TC`;
* the no-go theorem for collision-form certificates;
* the Fano-block cost `3j`, which extends J5.

None of it bears on actual circles.

**Scores** (rubric of `round10/README.md` §1; details in §6):

| write-up | (a) | (b) | (c) | **overall** |
| --- | --- | --- | --- | --- |
| refJ | 7 | 2 | 1 | **4** |
| loctheory | 6 | 1 | 1 | **4** |
| fatpoint | 6 | 2 | 1 | **3** |
| loccompute | 5 | 1 | 1 | **3** |

**Recommendation (§7).** Stop investing in route 5 until someone has a new idea for the fake test. Return to the only
open local class that needs no conjecture, P3-all\*, and spend the next round on its **construction side in the
regime that growth needs** (all `p_j <= M^(h_0)`, moduli up to `M^A`, `A` chosen after `c_0`). No fake of any kind
is known there.
* A member there caps the last local class, completing the barrier map. Only global inputs (Vojta, route 5) would
  then remain.
* A forced lower bound there would be the first local growth mechanism.

## 1. Referee-corrected statements

Conventions are those of `round12/README.md` §1.
* *Proved* means proved with all quantifiers in a refereed write-up.
* *Exact* means integer or rational computation, or mod-`p` computation used only in the valid direction.
  * For equalities this is the mod-`p` sandwich with exact generators.
  * For non-membership or a rank deficit, it additionally needs `rank_p = rank_Q` of the evaluation matrix
    (referee_fatpoint §2.5 gives a counterexample to the converse).
* *Evidence* is anything else.

Unless a model is named, the setting is `T = 1`, `w = sqrt(8 d ln 4r)`, `j = floor(w/2)`, `m' = ceil(4w)`, and the
regime `d >= 316 ln 4r`. **Every exact model this round is LP-subcritical**, and the Hadamard ones are also
ball-bound. **No exact model is a regime box.**

### 1.1 refJ (referee of J5-J7): survives, with one major error removed

1. **J5 (proved, finite `n`).**
   * `ZE-TC^coll_eps(A, m'; Q_0)` is false for every `Q_0 > 0` on Sylvester `n = 2^k`, `12 <= k <= 20`
     (`r = d`, `C = 1`), at every `abar >= max(2 ln 15, g*(Q_5))`. This includes all `W >= 3 n ln n - 14 n`.
   * `ZE-TC^coll_sigma` (per-subspace radii) is false for `k = 12, 13, 14` at `abar >= 3.89, 5.42, 5.42` (the referee's
     sharper thresholds; refJ's 4.89 / 6.20 / 6.42 are also valid). That is all near-extremal `W`.
   * Tools:
     * `rho(P_j) <= 3j` for the parity box `P_j = {|y_i| <= j, y_1 + y_2 + y_3 even}` (slice peeling);
     * a Fano block `{a, b, a+b}` of Sylvester columns has image in `P_j`, so it costs `3j`, not `6j`;
     * the 32-class span `L_{2,1}` lies in `TC(A, m')` for every `k >= 5`.
   * Conventions: `P(pi) >= 1`, and `C = 1` (case B depends on `ln(1/C)`). Nothing asymptotic follows.
2. **J6(a) (proved).** Over both `C` and `R`, `TC(A,m) ⊇ Bs(L_m)`, with multiplicity `m - tau_A(v)` at integer `[v]`
   with `tau_A(v) < m`.
   * Cheap spans are exactly the positive-multiplicity BDP linear cycles.
   * The base locus is strictly larger: in `M8` it contains the 3-dimensional Fano spans, which are not cheap
     spans.
3. **J6(b) (proved, self-contained).** The strong Lefschetz property of `Q[x]/(x_i^3)` is proved by an `sl_2` argument,
   and the Stanley 1980 citation was checked. For pairwise disjoint signed cuts, `J(A) = (x_i^3, l_j^(2L_j + 1))` in all
   degrees. So `Loc(1)` holds there, but these boxes are decomposable toys.
4. **J6(c) (exact, toy `M8`: `Syl(8)`, `L = 2`, `m = 8`; not LP-supercritical, `m/j = 4`).**
   * `J(A) = J_w +` the pullbacks from the 7 Fano-line triples, in all degrees `<= 8`.
   * `dim (J_w)_(6,7,8) = 616, 1485, 2940`, against `dim J(A)_(6,7,8) = 637, 1541, 2961`.
   * Column pairs and non-line triples add nothing over `Q`.
   * `Loc(2)` fails exactly for column pairs and for 51,933 integer planes (`630 < 637`). For all planes it is open;
     the evidence (75,545 planes searched) says it fails.
5. **J7(a) (proved; with `P(pi) >= 1`).**
   `coll_eps(Q_0) => coll_sigma => (MIN) => MIN^pair`, and `MIN* => (MIN)` iff `C <= 1`.
6. **J7(b), corrected.**
   * (1) If `max eps_xy < 2 ln((n-1)/2)`, the gap region `G` is empty.
   * (2') If `min eps_xy >= 2 ln(n/2)`, MIN^pair fails at every order in `(2jk, m_max]` for `k >= 16`. This needs
     `abar >= 2 ln(n/2)`, i.e. `W >~ 2 W_min`.
   * **refJ's Prop 5.2(2) (the max-eps version) is false.** The proof reverses an inequality, and four points with
     mutual `eps = 0` are a counterexample. The claimed "window of width `< 2/(n-1)`" and "no escape on non-empty
     gap regions" are withdrawn.
   * At near-extremal `W`, whether MIN^pair is weaker on non-empty gap regions is **open**.

### 1.2 loccompute (exact computations of Loc): computations survive, headline does not

1. **Locality onto columns, or onto columns plus non-base coordinates (`L_col(b)`, `L_mixed(b)`, with `J_w`).**
   * **`b_0 = 1` holds (exact)** in:
     * `Syl(8,1)`;
     * `M8' = Syl(8,3)`, for all degrees `<= 11`;
     * the 400 completed two-cut blocks with `d <= 7`;
     * 827 of the 953 three-cut types with `d <= 6`.
   * **`b_0 = 3` (exact)** in `M8`; `b = 2` fails.
   * **Paley-12 (`L = 2`, `m = 8`), degree 7, fails for every `b <= 5` (exact, independently re-verified).** In the
     `Z_11`-invariant part, `dim_Q L_col(5)_7 = 82 < 137 = dim_Q J(A)_7`. Mod `p` the failure persists to `b = 6`.
   * **`Syl(16,2)`.**
     * Degree 6 needs `b >= 3` (exact), and the 35 line families suffice exactly through degree 6.
     * In degree 7, among column families, the 7-column plane families are needed and leave a residual of 78 (mod-`p`
       evidence).
2. **README Loc(b0), with arbitrary integer directions: neither proved nor refuted in any model.** Non-column
   directions lower the observed `b_0`:
   * pairs give `b_0 = 2` in every three-cut block tested (exact);
   * in `Syl(16,2)`, pairs inside the line spans give `>= 505` of 540 dimensions in degree 6 (exact);
   * in Paley, `{k_base, x_a, x_b, x_c}` gives exact new generators outside `L_col(5)` (rank 82 to 84), yet **50
     invariant degree-7 dimensions remain unexplained by every family the referee tried**: pairs from 320 directions,
     triples from 89, quadruples from 23;
   * in the `Syl(16,2)` plane ring, families of `<= 3` directions give 147 of the 231 plane-only generators (mod `p`).
3. **Lemmas (proved).**
   * Grid images contribute nothing.
   * Congruences give integer directions.
   * The cube-and-slab vanishing ideal is `(x^3 - x, g_p(y_p))`.
4. **Withdrawn.**
   * "Loc refuted as a model-independent principle": only column locality is refuted.
   * "Regime `b_0 >= 7`, plausibly 15": it counts columns, not directions.
   * "Nothing smaller supplies the line generators": pairs supply 8 of the 9.
   * The count of 440 syzygies is a mod-`p` upper bound. Paley has 1833 point orbits.
   * The README §7.3(3) column projection-cost route fails at its first step. In Paley, the contributing non-column
     families also have `rho = m`.

### 1.3 loctheory (theory of locality): restricted theorem survives, headline claims do not

* **(A) Linkage criterion (proved, any profile, any finite `P` inside a grid `Gamma`).**
  * `H_P(e) + H_(Gamma \ P)(s - e) = H_Gamma(e)`.
  * `J(P)_e = J(Gamma)_e` for `e < alpha(P) = s - rho(Gamma \ P)`, and the inclusion is strict for
    `alpha <= e <= s`.
  * So a family contributes beyond its grid exactly in degrees `[alpha_F, s_F]`, and the family set
    `{alpha_F <= m'}`, together with single directions, suffices for Loc.
  * **The literal criterion "F contributes only if `rho(pi_F A) < m'`" is false.** Witness:
    `P = {y in [-2,2]^5 : sum y = 0 mod 3}` has `rho = 10` and `alpha = 7` (exact).
  * Whether restricting Loc's family set to `rho < m'` actually loses generators is **open**. Referee evidence
    suggests it does.
* **(B) Parity sections (proved, any profile).**
  * Full-Venn sign families have image `[-j,j]^b` and `rho = 2jb`, and contribute only width powers.
  * If `pi_F(A)` is a parity section `P_C`:
    * `H` is an orbit count;
    * `J_F` = width powers + the sum over circuits `s` of `J_(F_s)`;
    * a circuit enters exactly at degree `j|s|`;
    * `rho = alpha = jc` for one `c`-circuit;
    * `rho(P_C) >= (j-1)b`.
  * In degrees `<= m'` only circuits of size `<= 8` matter. **This is an upper bound only.** At toy `j`, 4-circuits
    are generated by 2-3 non-column directions (certified), so circuit size is not a lower bound on `b_0`.
* **(C) String theorem (proved).** For boxes cut by atom-constant columns:
  * `H` and `J` reduce to lattice sets in atom space;
  * `J = (x^3) + tau^* J(T)` iff a box-intersection property (BIP) holds. BIP is exact in 5 two-cut cases and open
    in general.
  * The theorem generalises (and uses) the one-cut theorem.
  * The equivalence with lattice-set locality holds only for atom-space directions.
* **(D) Sylvester images (proved).**
  * If `2^t ln(2n) < j`, `n >= 2^(t+1) j` and `j >= 2^(t-2)`, every column family of `F_2`-rank `t` has image
    `Lambda_F ∩ [-j,j]^F`.
  * **"Binary unless `F` contains a complete `>= 3`-block" is false.** There are 890 block-free non-binary
    families at `t = 4`. The smallest is `AG(3,2)` (the 8 odd-weight vectors of `F_2^4`, with `sum y = 0 mod 4`).
* **(E) Complete XOR blocks (proved).**
  * `rho(Y_t(j))/j` has asymptotic lower bounds `3.00, 4.16, 5.87, 8.30, 11.76, 16.64` for `t = 2..7`.
  * `rho(Y_6(j)) > 9j + 9` for `j >= 91`.
  * So refJ's sufficient `rho_t` failure condition never holds for `t >= 6`, `2^t <= n/(2j)`. This does not verify
    the strata condition for block spans.
* **(F) Toy data (exact, `j = 2` and `j = 1`).**
  * `Y_3(2)`: `J_7` exceeds grid + lines + quadrangles + the 43 cheap width powers by 378 over `Q`.
  * **`J(Y_3(2))` in degrees `<= 8` is generated once the `GL(3,2)`-translates of one 6-direction family are added**
    (certified). The flagship non-locality example therefore dissolves.
  * `AG(3,2)` at `j = 1` is local with 7-direction families (certified).
  * **No computed example refutes `Loc(b_0)` with `b_0 = |F| - 1`.**
* **Withdrawn.**
  * "`b_0 = 8`, not `~4`; README too optimistic by a factor 2": unproved.
  * "Complete blocks are the precise candidate for the failure of Loc": the data point dissolves.
  * "`Y_3(2)` is a regime image": Lemma 5.1 never applies at `j = 2`.

### 1.4 fatpoint (base loci, certificates, barrier test): proofs survive, decisive-test reading does not

1. **Base loci (proved, any finite `A` containing 0, any `m`, any family set; `I_F` = the ideal generated by
   `pi_F^* J(B_F)`).**
   * `Bs_R(∩_F (I_F)_m^perp) ⊂ TC(A,m)` **unconditionally**, with equality under `Loc_m` in this generated-ideal form.
     That form is weaker than the README span form, which implies it.
   * Single-family base loci lie in `span_Q F`, and joins with `sum nu < m` are base points.
   * Span principle: a nondegenerate base point needs a `Q`-irreducible base component spanning `C^d`.
2. **Product certificates (proved).**
   * `prod <z_i,u>` lies in `K_m^perp` iff the polarisations of `K_m` vanish at `(z_1, ..., z_m)`.
   * Width counts are exact for `b = 1`.
   * `N_F < alpha_F` with grid counts suffices, and `N_F <= rho(B_F)` is necessary.
3. **Regime no-go results (proved).**
   * Thm 3.1: no collision-form certificate exists when `W > 8 ln C'`.
   * Prop 3.2: generators are heavy on average.
   * Thm 3.3: the uncertainty budget does not bind, though it is within about 1.03 at the regime boundary.
4. **Parity families impose only counts (proved, given loctheory (B)).**
   * `<= j|s| - 1` per circuit suffices, and `<= j|s|` is necessary.
   * Complete blocks reduce to class-sum products.
   * Complete blocks are **not** the only non-binary families.
5. **`Y_3(2)` (exact, toy set, not a regime image).**
   * Among the 91,272 heavy multisets, a product lies in `J^perp_7` **iff it is a lattice zonotope fitting in
     `Y_3(2)`** (1536 members).
   * The negatives are rigorous, because `rank_Z K_(<=7) = rank_p = 2353`.
   * "Products with zonotope cost `<= 2j` are members" is wrong: half of them are not.
   * 6,492 multisets already fail a `b_0 = 1` width count. Extrapolation to the regime is heuristic.
6. **Conditional theorems (proved implications).**
   * Thm 6.1/Cor 6.2: Loc + local certificate + height + tight characters imply `n < n_1`. The asymptotics need
     `W >= r ln r`.
   * **Thm 6.4/Cor 6.6:**
     * `BLOCK(t_0, j_0, m_0)` (ZE-TC^braid for one box `Syl^box(2^(t_0), j_0)` at `m_0 > 8 j_0`) implies MIN^pair
       for Sylvester-type profiles, hence `n < n_0` for them.
     * It does not apply to fake or generic (full-trace) profiles.
     * The regime threshold is `n >= 2966`.
     * **See J8 and J9 (§8): `BLOCK` is false at `t_0 = 5`, and the headline application is already unconditional.**
7. **Exact models (all LP-subcritical).**
   * Certified `Loc_8(1)` for sdc7, and `Loc_10(1)`, `Loc_11(1)` for `M8'`.
   * **The sdc7 width system has certified nondegenerate, well-separated real base points**. The model is
     underdetermined: `H_A(8) = 1`.
   * Within the search pool (support `<= 4`, entries `±1`) every certificate found is a zonotope. Non-existence for
     `M8` (`m = 8`) and `M8'` (`m = 10`) holds only within that pool.

## 2. Status of Loc(b0)

The statement under test is README-12 §7.1. For every admissible regime profile and every `i <= ceil(4w)`,
`J(A)_i` is the sum of `(pi_F^* J(pi_F A))_i` over families `F` of `<= b_0` integer directions with
`rho(pi_F A) < ceil(4w)`, with `b_0` independent of `n`.

| statement | status | source |
| --- | --- | --- |
| Loc(b0), README form, any regime profile, any `b_0` | **open**: neither proved nor refuted | all four |
| The restriction `rho(pi_F A) < m'` in the family set | its justification (the literal contribution criterion) is **false**. Whether the restriction loses generators is **open**, and evidence (Paley) says yes | loctheory Prop 2.4; referee_loctheory §7.1; loccompute §4 |
| Correct contribution criterion: `F` contributes exactly in degrees `[alpha_F, s_F]` | **proved** | loctheory Thm 2.2, Cor 2.3 |
| Loc(1) for one cut and for disjoint signed cuts | **proved** (SLP) | refJ Thm 3.3 |
| Locality of parity sections through circuits (all degrees) | **proved** | loctheory Thm 3.4 |
| Atom-constant cut boxes: Loc reduces to lattice sets in atom space | **proved** given BIP, and for atom-space directions only. BIP is open beyond 5 exact cases | loctheory Thm 4.1, Prop 6.2 |
| Column (and non-base coordinate) locality | **refuted, exact**: Paley-12, degree 7, every `b <= 5` (mod `p`: `b <= 6`) | loccompute Prop 4.4; referee re-verification |
| `b_0` in toy models, arbitrary directions | 1 (Syl(8,1), M8', cut blocks, sdc7); 2 (three-cut blocks); 3 (M8); `<= 6` (`Y_3(2)`, degrees `<= 8`); `<= 7` (`AG(3,2)`, `j = 1`); Paley-12: `> 5` columns, and 50 invariant dimensions are unexplained by every family tried | exact, except where marked mod `p` |
| Boundedness of `b_0` in `n` | **no evidence**. Every model has `m/j <= 4.5`, against `~8` in the regime, and is LP-subcritical | — |
| Lower bounds on `b_0` from circuit size, "`b_0 >= 8`" or "`>= 7`, plausibly 15" | **unproved / withdrawn** | referee_loctheory §5.3; referee_loccompute |

**Reading.** On the toy evidence, the smallest generating families found grow with the arithmetic complexity of the
congruences on the image. Only M8 and Paley give lower bounds; the other counts are upper bounds. In order:
* none (grids, cut boxes);
* mod 2 lines (M8);
* complete `F_2`-blocks;
* mod 4 (`AG(3,2)`, Sylvester planes);
* mod 3 five-column relations (Paley).

Regime boxes have images with congruences of every small order, at a ratio `m'/j` twice as large as in any model.
Nothing computed or proved suggests that `b_0` stays bounded. Nothing refutes it either. The honest status is
**open, with weak evidence against a small `b_0` and none on boundedness.**

**Loc would not finish the job even if it held.** Under Loc, `TC(A,m')` is the real base locus of an intersection of
local inverse systems (fatpoint §2). The barrier question then becomes whether that base locus has nondegenerate,
collision-poor real points on fake pattern boxes. The only certified data point goes the wrong way: sdc7 satisfies
`Loc_8(1)` exactly and has nondegenerate, well-separated real base points. That model is underdetermined, so it says
nothing about the regime. But it shows that locality is not itself a non-degeneracy mechanism.

## 3. Status of the base-locus / barrier test

The barrier test, from fake Thm 2.1 (round 12): a nondegenerate real `TC_H` point of a (Mar) fake pattern box that
satisfies the Diophantine condition (2.2) yields gap-admissible fakes that refute MIN^pair there.

| statement | status | source |
| --- | --- | --- |
| `Bs_R(local system) ⊂ TC(A,m)` for every family set (no Loc needed); `TC ⊇ Bs(L_m)` with multiplicity `m - tau_A(v)` | **proved** | fatpoint Lemma 2.1; refJ §2 |
| Cheap spans = positive-multiplicity linear base cycles; the base locus is strictly larger (Fano 3-spaces in M8) | **proved** / **exact** | refJ Prop 2.3 |
| Certificates of TC points by products of collision forms, in the regime | **impossible** for `W > 8 ln C'` | fatpoint Thm 3.1 |
| Product-certificate generators must be heavy characters | **proved** (on average) | fatpoint Prop 3.2 |
| Certificates in exact models | zonotopes only (within the search pool). In `Y_3(2)`, heavy products are members **iff** they are fitting lattice zonotopes | fatpoint §5; referee_fatpoint §2.7 |
| Nondegenerate real base points of local systems | **exist** in sdc7 (underdetermined, `Loc_8(1)`, certified); **unknown** in every overdetermined model | fatpoint §5.2 |
| Nondegenerate real TC points of regime-size fake pattern boxes (the barrier question) | **open, no evidence either way** | — |
| Running the test computationally | **impossible** at regime size. The smallest LP-supercritical Sylvester model has `n >= 32`; regime boxes need `n >= 2966` | zetc §5; fatpoint §5.1 |

**Reading.** Round 13 gives the right language for the barrier test. Base loci of local systems are genuine TC points,
and certificates are polarisation conditions. It also proves two no-go results for the cheapest certificates
(collision forms, light characters). It does not run the test. In the toy models every certificate found is a
zonotope, and round 11 (compute Prop 2.2) showed that zonotope certificates are absent in the regime. So **neither
direction of the barrier test has a viable method**:
* to cap route 5, one needs a nondegenerate certificate beyond zonotopes, of which no example exists;
* to beat the barrier, one needs a quantitative non-existence statement on overdetermined boxes, for which no tool
  exists.

## 4. Status of MIN^pair

MIN^pair: `min{kappa_H(u) : |u|_inf = 1, |u_x - u_y| >= e^(-eps_xy/2) for all pairs} > eps_0` at the cluster's own
data. It is the weakest known sufficient form, and with Step 0, separation and rigidity it gives `n = O(loglog R)` on
all norms.

| statement | status | source |
| --- | --- | --- |
| MIN^pair for regime boxes (generic, fake, or actual profiles) | **open**; no evidence either way | — |
| Hierarchy `coll_eps(Q_0) => coll_sigma => (MIN) => MIN^pair`; `MIN* => (MIN)` iff `C <= 1` | **proved** (with `P(pi) >= 1`) | refJ Prop 5.1; referee_refJ |
| `coll_eps(Q_0)` | **false**, every `Q_0`, `n = 2^12..2^20`, near-extremal `W` | refJ §1 |
| `coll_sigma` | **false**, `n = 2^12, 2^13, 2^14`, at `abar >= 3.89 / 5.42 / 5.42` | refJ §1.6; referee_refJ |
| MIN^pair at Sylvester large orders `(2jk, m_max]`, `k >= 16` | **fails** if `min eps_xy >= 2 ln(n/2)` (`W >~ 2 W_min`). The max-eps version (refJ Prop 5.2(2)) is **false**. At near-extremal `W`: open | referee_refJ §5.2 |
| `BLOCK(t_0, j_0, m_0) => MIN^pair` for Sylvester-type profiles | **proved implication** | fatpoint Thm 6.4 |
| `BLOCK(5, j_0, m_0)`, any `j_0`, any admissible `m_0 > 8 j_0` | **false** (J8, unrefereed; exact inputs) | §8 |
| `BLOCK(6, j_0, m_0)` | **false** for `m_0 > 9 j_0`, and for `j_0 = 2` at every `m_0 > 16` (J8). **Open** for `j_0 >= 3`, `8 j_0 < m_0 <= min(9 j_0, 2 rho(Y_3(j_0)))` | §8 |
| Sylvester-type actual clusters (the class Thm 6.4 serves) | already **uniformly bounded, unconditionally** (walsh-extraction Theorem B: `M <= 24` for `C <= 1`) (J9) | `growth/walsh-extraction.md` §4 |

**Reading.** MIN^pair is exactly as open as it was after round 12. The two new handles do not change that.
* The block-product reduction applies only to profiles with low trace complexity on small blocks. Fakes and generic
  profiles have full traces, so it does not touch the barrier.
* The corrected J7(b) dichotomy needs `W >~ 2 W_min`, so it does not reach the near-extremal range.

## 5. What it means for route 5 and the growth rate

**Growth rate: unchanged.** `M <= (2/3 + eps) log R/loglog R` is the best bound for actual circles. Nothing this round
improves the growth rate or touches the uniform bound, unconditionally or from a standard conjecture. Constants are
not counted.

**Conditional theorems on record after round 13** (only the changes against round 12 are listed):

| hypothesis | conclusion | implication | hypothesis status |
| --- | --- | --- | --- |
| MIN^pair, cluster's own `T = 1` box at `ceil(4w)`, all regime profiles | `n = O(loglog R)`, all norms | proved (R12) | **open**, unchanged |
| `ZE-TC^coll_sigma` | same | proved | **false** at `n = 2^12..2^14` (near-extremal `W`); asymptotically open |
| `ZE-TC^coll_eps(Q_0)` | same | proved | **false** for every `Q_0`, `n = 2^12..2^20` |
| Loc(b0) + local certificate + height bound + tight characters | `n < n_1` | proved (fatpoint Thm 6.1) | every hypothesis open; certificates beyond zonotopes unknown |
| `BLOCK(t_0, j_0, m_0)` | `n < n_0` for Sylvester-type profiles | proved (fatpoint Thm 6.4) | false at `t_0 = 5` (J8). The class it serves is already bounded unconditionally (J9) |

**Meaning for route 5.**
1. **The residual hypothesis has not moved in three rounds.**
   * Round 11 reduced ZE to the tangent cone, round 12 made the reduction rigorous on all norms, and round 13 studied
     the algebra of `J(A)`.
   * After all three, the open step is still the same single inequality, and it is still equivalent (fake Thm 2.1) to
     the barrier question on fake pattern boxes.
2. **The structural handle proposed in round 12 did not deliver.**
   * Loc(b0) is undecided, and its only justification (the `rho < m'` contribution criterion) is false.
   * Even a proof of Loc would leave a base-locus question of the same difficulty (§2, sdc7).
   * The per-stratum check that Loc was meant to enable has been done only for binary column spans (loctheory
     Cor 5.3'), and is trivial there.
3. **No method reaches the regime.**
   * Exact computation stops at LP-subcritical toys (`n <= 16`, `m/j <= 4.5`), while the regime starts at
     `n = 2966` with `m'/j ~ 8`.
   * General Lojasiewicz bounds cannot give the size (nearstrata Prop 2.7).
   * Certificates beyond zonotopes are unknown, and zonotopes are absent in the regime.
   * Block products work only for low-trace profiles, which are exactly the profiles fakes avoid.
4. **Fake-violated input (rubric b): unchanged.**
   * The input fakes violate is rigidity at the actual point.
   * The step that must use it, "no nondegenerate fake direction in the near tangent cone", is untouched.
   * Theorem 6.4 isolates a profile property (small blocks with few traces and an LP-supercritical tangent cone) that
     separates Sylvester-type profiles from fakes. But those profiles are already bounded by walsh-extraction
     Theorem B (J9).
5. **Verdict: route 5 is stuck, not capped.**
   * Stuck: no proof strategy for MIN^pair or for ZE-TC on overdetermined regime boxes exists after three rounds.
   * Not capped: nothing shows that fake pattern boxes have nondegenerate real tangent-cone points.
   * Route 5 should be **paused**. It should be resumed only with a genuinely new handle on the barrier question for
     overdetermined fake boxes. More algebra of `J(A)`, more toy models and more reformulations of the final
     hypothesis are not new handles.

## 6. Scores

Rubric of `round10/README.md` §1:
* (a) correctness after refereeing;
* (b) fake-violated input;
* (c) movement of the growth rate or of uniformity.

The overall score is a judgement and weighs (b) and (c) above (a). No write-up moves the growth rate.

| write-up | (a) | (b) | (c) | **overall** | reason |
| --- | --- | --- | --- | --- | --- |
| **refJ** | 7 | 2 | 1 | **4** | Confirms J5-J7 with independent code, and supplies the round's most useful tool: the parity-box bound `rho(P_j) <= 3j` and the Fano-block cost `3j`. These extend J5 to `n <= 2^20`, refute `coll_sigma` at `2^13`, `2^14`, and give J8 below. It also gives a self-contained SLP proof with signed and disjoint cuts, and the exact M8 structure. Deductions: one major false proposition (Prop 5.2(2), a reversed inequality, which turned "open" into a false "dichotomy"); the toy scope of M8 and of the cut boxes is not flagged. |
| **loctheory** | 6 | 1 | 1 | **4** | Real theorems: the linkage criterion `alpha_F` (which corrects README §7), locality of parity sections through circuits with the orbit Hilbert function, the string theorem, and the complete-block `rho` bounds. Deductions, all from the headline: the binarity census is false (`AG(3,2)`); "`b_0 = 8`" is unproved, since circuit size is not a lower bound; the flagship non-locality example (`Y_3(2)`) dissolves under 6-direction families; `Y_3(2)` is mis-scoped as a regime image. |
| **fatpoint** | 6 | 2 | 1 | **3** | Correct proofs throughout: the unconditional base-locus inclusion, polarisation criteria, the no-go Thm 3.1, Prop 3.2, the block-product reduction Thm 6.4 with a correct averaging argument, and certified sdc7/M8' locality. Deductions: the mod-`p` rigour argument is false (the negatives survive only by the referee's rank check); the decisive `Y_3(2)` test is misread (membership = fitting zonotopes; an omitted width count); refuted sibling claims are repeated as facts; Loc is silently redefined. Thm 6.4's smallest block is false (J8), and its headline application is subsumed (J9). |
| **loccompute** | 5 | 1 | 1 | **3** | A large, careful exact effort. It gives the exact refutation of column locality in Paley-12, exact `b_0` in M8, M8', `Syl(16,2)` (degree 6) and the cut blocks, and three proved lemmas. Deductions: every `b_0` is measured on columns while the conclusions are drawn at README level, which non-column pairs overturn; "planes needed" and "`b_0 >= 7`, plausibly 15" are unsupported; there are labelling slips. |

Partial credit in (c) is only for sharpening the map of what route 5 can and cannot do. That covers refJ's extension
of J5 and fatpoint's no-go theorems. No (c) score exceeds 1 because no statement about actual circles moved.

## 7. Single most promising next step

### 7.1 Route 5: pause it

Route 5 is stuck (§5). The record offers no route-5 task that could move the growth rate within a round:
* **Loc** is undecided, and even a proof would not decide the barrier.
* **BLOCK** serves only profiles that are already bounded.
* **MIN^pair** is the barrier question itself.
* **The regime** is beyond exact computation.

If route 5 is kept alive at low cost, one task is worth it: the **first overdetermined data point for the barrier
test**, i.e. the unfinished support-4 certificate search on `pal11` (`L = 2`, `m = 8`) with its `Z_11` symmetry.
* Use fatpoint's polarisation criterion.
* Certify non-membership with `rank_p = rank_Q`.
* Look for any certificate that is **not** a zonotope.

A non-zonotopal nondegenerate certificate on an overdetermined fake box would be the first sign that route 5 is
cappable. A complete negative would be weak evidence the other way. Either way it is toy data.

### 7.2 Recommended alternative: P3-all\* on its construction side, in the regime that growth needs

**Why this line.** The whole record (`docs/claude_review/README.md`; rounds 7-10) gives three reasons.
1. **It is the only local input class that is valid for actual circles, needs no conjecture, and is not known to
   be capped** (round 7 §2.1 row 5; round 8 §3).
2. **Growth needs it only in regime (1.4)**: primes `p_j <= M^(h_0)` with `h_0 = 2 c_0 - 3 + o(1)`, and moduli up to
   `X = 3 M^A`, with `A` chosen after `c_0` (round 9 §0 item 1). **No fake of any kind is known there.**
   * Every construction of rounds 5-7 has `p_j > X`.
   * Every known `O(M log M)` fake needs `P >~ r^2 log r`.
   * P3-G fakes are proved only for `log p >= 2.08 X`.
   * So it is the one place where the fake barrier has not been established. Route 5 had to beat a barrier that is
     known to exist in full strength; this line does not.
3. **Either outcome is decisive, and the first step is finite.**
   * *Construction succeeds* (C1 in regime (1.4)): low-demand, parity-respecting factored residue data good for
     **every** monomial, together with the slab conditions. This caps P3-all\* where growth needs it and completes
     the local barrier map. Then every local input is capped, and only global inputs (Vojta; route 5) remain. That
     is a theorem-level redirection, and it would justify returning to route 5 with a clear conscience.
   * *Construction is obstructed* (all-norm goodness forces genuine labels on a growing range of moduli): this is
     the first local rigidity mechanism. It reduces P3-all\* in the regime to P3-G, where only the label-angle step
     L7 is left. That is the first concrete route to (T) from round 8 §3, i.e. `M = o(log R/loglog R)`
     unconditionally.

**Statement to decide.**

> **Rig(K).** Regime (1.4), `r <= rho M`. Let a parity-respecting label system be good, i.e. meet the per-unit budget
> `m_(v,s) <= nu_s^(-1) e^(w(v)/2)`, for **every** `v` with `||v||_1 <= K`, at every modulus `q <= Q_K`. Does the
> range of moduli on which it must be system-level coherent with the genuine labels (the "rigid core" of
> round 10 §2 item 5) grow with `K`, and polynomially in `M`?

Round 10 (`L4K.md` Thm 1-2) settled the bounded-`K` version negatively: bounded norms force nothing above
`M^(floor(K/2) + o(1))`. It also found exact evidence that the rigid core moves up as `K` grows:
* at the norm-5-forced moduli 167 and 199, there are 96 and 40 non-genuine survivors, against 90,048 and 500,160 at
  norm 4;
* 211 is flexible.

Rig(K) for all `K` is the all-norm question that P3-all\* actually imposes.

**First finite test (exact; evidence).** This is round 10 §8.6 item 2, extended.
* Take `r = 4..8` split primes `<= 400` and `K = 4..7`, with moduli up to the norm-`K` forced/free threshold.
* Compute, by ILP or SAT over signatures with exact budgets, the largest modulus up to which every good system is
  coherent with the genuine one.
* Fit its growth in `K`, and compare it against `M^(floor(K/2))`.
* Use the harness in `round10/L4K_checks/` and `referee_L4K_checks/`.
* Certify negatives (infeasibility) exactly, and treat positives (non-genuine survivors) as explicit witnesses.

**Then.**
* If the rigid core grows, prove Rig(K) at the Prop. H scale (Linnik + Heath-Brown already give L6 in this regime).
* If it stalls, take the flexible tops as the seed of a C1 construction, and add the slab conditions with angles
  chosen after the labels.

**Not recommended.**
* Further algebra of `J(A)`, and Loc in any form.
* Reformulations of route 5's final hypothesis.
* BLOCK (it serves bounded profiles).
* Toy LP-subcritical models as evidence about the regime.
* Robust analytic inputs (capped, round 10).

## 8. Judge's own findings (unrefereed)

Two findings. J8 is backed by exact computation (`judge_checks/block_t0.py`, output `block_t0_out.txt`); J9 is a
comparison with the record. **Neither has been refereed.**

### 8.1 J8: the smallest block of fatpoint Thm 6.4 fails; `BLOCK` needs `t_0 >= 6`

**Setting.** `Syl^box(N_0, j_0) = {k in {-1,0,1}^(F_2^(t_0)) : sum k = 0, |hat k(a)| <= 2 j_0 for all a != 0}`, where
`N_0 = 2^(t_0)`. Since `sum k = 0`, `hat k(a) = 2 k(S_a)` with `S_a = {x : a.x = 0}`, so the constraints read
`|k(S_a)| <= j_0`. `BLOCK(t_0, j_0, m_0)` with `m_0 > 8 j_0` is, by fatpoint Lemma 6.5, ZE-TC^braid for this box at
order `m_0`: every real point of `TC(Syl^box, m_0)` has two equal coordinates.

**Proposition J8.**
* (a) For every `j_0 >= 1` and every `m_0 > 8 j_0`, `TC(Syl^box(32, j_0), m_0)` contains a real point with 32 distinct
  coordinates. Hence `BLOCK(5, j_0, m_0)` is false.
* (b) For `t_0 = 6`, the same holds for every `m_0 > 9 j_0`. At `j_0 = 2` it holds for every `m_0 > 16 = 8 j_0`.
  Hence `BLOCK(6, j_0, m_0)` can hold only if `j_0 >= 3` and `8 j_0 < m_0 <= min(9 j_0, 2 rho(Y_3(j_0)))`.

*Proof.*
1. **Fano blocks.** For `a, b` independent, every `x` lies in all three of `S_a, S_b, S_(a+b)` or in exactly one. So
   `k(S_a) + k(S_b) + k(S_(a+b)) = 2 k(S_a ∩ S_b)` when `sum k = 0` (refJ Lemma 1.4). The image of the box under these
   three columns therefore lies in the parity box `P_(j_0)`, and `rho(P_(j_0)) <= 3 j_0` (refJ Lemma 1.3, refereed).
2. **The family for (a).** Take `F = {e_1, e_2, e_1+e_2, e_3, e_4, e_3+e_4, e_5}`. Then
   `pi_F(Syl^box) ⊂ P_(j_0) x P_(j_0) x [-j_0, j_0]`. Interpolation degree is monotone under inclusion and additive on
   products, so `rho_A(F) <= 3 j_0 + 3 j_0 + 2 j_0 = 8 j_0`.
3. **From cost to tangent cone.** For `u` in the span of the indicators `1_(S_c)`, `c in F`, we have
   `<k,u> = <pi_F k, y>`, a linear form on the image. All forms of degree `> rho` lie in the leading-form ideal of the
   image; pulling back gives `<k,u>^(8 j_0 + 1) in J(A)`. So `nu_A(u) <= 8 j_0 < m_0` (zetc §2.2(1)), i.e.
   `u in TC(A, m_0)`.
4. **A collision-free point in the span.** The indicators of a Fano line, together with the constants, span all
   functions of the two bits involved. The span of `F` plus constants therefore contains `x -> x_1 + 2x_2 + 4x_3 +
   8x_4 + 16x_5`, which takes 32 distinct values (checked exactly by linear algebra). Since `sum k = 0`, constants do
   not change `<k,u>`, and normalising at a base point gives the required point. This proves (a).
5. **(b), general `j_0`.** Use three Fano lines on 6 independent `a`'s, which cost `<= 9 j_0`.
6. **(b), `j_0 = 2`.** Use two complete 3-dimensional blocks `V_1 ⊕ V_2 = F_2^6`. Their images lie in
   `Y_3(2) x Y_3(2)`, and `rho(Y_3(2)) <= 8` (exact, below). Their span plus constants contains an injective
   function on `F_2^6` (checked exactly).
7. **The window.** Any surviving `m_0` must also be `<= 2 rho(Y_3(j_0))` by the same argument. QED

**Exact data** (`block_t0_out.txt`). Here `rho <= e` is certified when the mod-`p` rank of the degree-`<= e` evaluation
matrix equals the number of points. That is the rigorous direction, since `rank_p <= rank_Q <= |Y|`.
* `rho(P_j) <= 3j` for `j = 1..5` (`|P_j| = 13, 63, 171, 365, 665`).
* `Y_3(2)` (the lattice section, lattice index 32): `|Y_3(2)| = 2605`, matching fatpoint and loctheory.
* `rho(Y_3(2)) <= 8`, with `rank_p(deg <= 7) = 2353`, matching referee_fatpoint.
* `rho(Y_3(1)) <= 4`.
* The Fano-block identity was tested on 2000 random exact samples.
* The 5-bit and 6-bit codes lie in the respective spans (exact `sympy` solve).

**Consequences.**
* fatpoint's remark that `BLOCK` needs only `t_0 >= 5` (from the coordinate-column span, cost `2 j_0 t_0`) and the
  referee's remark that mixed spans "do not change `t_0 >= 5`" are both superseded. The Fano-block spans of refJ cost
  about `1.5 j_0` per bit, against `2 j_0`.
* The "smallest candidate block `Syl^box(32, j_0)`" of fatpoint §6.3 item 6 is excluded.
* Complete blocks cost at least about `1.39 j` per bit asymptotically (loctheory (E): `rho(Y_3(j)) >~ 4.16 j`). So
  for large `j_0` a window `(8 j_0, ~8.3 j_0]` at `t_0 = 6` is not excluded by any span known here. That window is
  narrow, and nothing suggests that ZE-TC^braid holds in it.

### 8.2 J9: Theorem 6.4's headline application is already an unconditional theorem

fatpoint Cor 6.6 concludes, conditionally on `BLOCK`, that clusters whose split primes realise exactly the
Sylvester columns (`n = 2^k`, `r = n - 1`, `e_p = 1`) have `n < n_0`. These clusters are pure Hadamard-core clusters
in the sense of `growth/walsh-extraction.md` §4. There, `z_x = g eps_x prod_a G_a^(H(x,a))` with `H` Sylvester and
`G_a` conjugate-primitive, possibly 1, with disjoint prime supports; a single split prime per column is the case
`G_a = pi_a`. walsh-extraction Theorem B, which is refereed, unconditional and effective, gives `M <= 24` for `C <= 1`
and `M < max(28, C^(2+o(1)))` in general. Missing columns (`G_a = 1`) are allowed there too.

So Cor 6.6 adds nothing for actual Sylvester clusters. Its "more generally" clause covers sign profiles admitting
random coset decompositions with Sylvester traces. That clause might reach a few non-Hadamard profiles, for instance
character columns on a union of cosets. It is a restricted class and conditional on `BLOCK`, which is open and,
by J8, false at `t_0 = 5`.

The value of Theorem 6.4 is structural only. It shows how MIN^pair could be inherited from one fixed block by slab-
budget splitting. It also shows that this mechanism needs low-trace blocks, which fakes do not have.

## 9. File index

All paths are relative to `round13/`.

| write-up | notes | report | scripts (author / referee) |
| --- | --- | --- | --- |
| `refJ.md` | `refJ_notes.md` | `referee_refJ.md` | `refJ_checks/` / `referee_refJ_checks/` |
| `loccompute.md` | `loccompute_notes.md` | `referee_loccompute.md` | `loccompute_checks/` / `referee_loccompute_checks/` |
| `loctheory.md` | `loctheory_notes.md` | `referee_loctheory.md` | `loctheory_checks/` / `referee_loctheory_checks/` |
| `fatpoint.md` | `fatpoint_notes.md` | `referee_fatpoint.md` | `fatpoint_checks/` / `referee_fatpoint_checks/` |

**Judge's checks** (`judge_checks/`, unrefereed): `block_t0.py` and `block_t0_out.txt`, for J8 (§8.1).

**Corrections to apply when round 13 is merged.** The repository is read-only.
* `round12/README.md` §7.1-7.3. The family restriction `rho(pi_F A) < ceil(4w)` rests on a false contribution
  criterion; replace it with `alpha_F <= m'` (loctheory Thm 2.2). §7.3(3) "`b_0 <= m'/(2j) + O(1)`" has no proof,
  and its column version fails (Paley).
* `round12/README.md` §4.1 J5. Extend it to `n = 2^12..2^20` (uniform `Q_0`) and refute `coll_sigma` at `2^13`, `2^14`
  (refJ). Use `C = 1` and `P(pi) >= 1`. The near-extremal constant is `13.91 n`.
* `round12/README.md` §4.2 J6(a). "Cheap spans are exactly the linear base cycles" holds only for positive-
  multiplicity BDP cycles; the base locus also contains Fano 3-spaces (M8). J6(b): the SLP is now proved
  self-contained (refJ §3).
* `round12/README.md` §4.3 J7(b). It stays **open**; refJ's Prop 5.2(2) is false.
* `fatpoint.md` §5.1. Mod-`p` non-membership needs `rank_p = rank_Q`. §6.2: `t_0 >= 6` (J8). Cor 6.6: subsumed by
  walsh-extraction Thm B (J9).
* `loctheory.md` §5.2. The binarity census is false (`AG(3,2)`). §5.3: `Y_3(2)` is local with 6-direction families.
* `docs/claude_review/README.md` §3. Add round 13 as follows: Loc(b0) undecided; column locality refuted
  (Paley-12); MIN^pair unchanged; route 5 paused; next round on P3-all\* in regime (1.4).

**Repository note.** Five generated data files over 1 MB are not committed:
`refJ_checks/m8_step7.pkl`, `referee_fatpoint_checks/enum.pkl`,
`referee_loccompute_checks/paley_V.npy`, `referee_loccompute_checks/paley_Einv.npy` and
`referee_loccompute_checks/syl16_A.npy`. The scripts next to them regenerate them.
