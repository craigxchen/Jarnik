# Referee report on round12/zetc.md (ZE-TC: tangent degree, one cut, few cuts, products, small models)

Referee: adversarial. All checks by independent scripts in referee_zetc_checks/.

## 0. Verdict

**The write-up survives with corrections. No growth improvement or uniform bound is claimed, and none follows. ZE-TC in the regime remains open.**

* **Correct.** Lemma 1.1 (all parts) and the statement of Theorem 2.1 are correct. The slice step of the proof has a degree error (C3), which I repaired and checked exactly. Cor 2.2(1)-(2), Props 3.1 and 3.3, Lemmas 3.2, 4.1, 4.3, 4.4 and 5.1, and Theorem 4.2 are correct, the last conditional on J4 (unrefereed but soundly derived from refereed data).
* **Exact claims I re-ran with my own code** (exact or valid-direction mod `p`) **all agree**:
  * one-cut formula, 8237 exact cases;
  * Cor 2.2(3);
  * every two-cut block type with `|S_1 u S_2| <= 7`;
  * 843/843 three-cut types;
  * `Z6`;
  * all Hadamard counts and `rho` values;
  * the Paley-12 `L = 2` spanning certificate;
  * the `Syl(16,3)` fibration bound;
  * the Paley-12 section octics. My check is exact, where the author's sign test was floating point.

  * the `Syl(16,2)` bound, independently and in sharper form (`rho <= 11`, against the author's `<= 15`), so "Sylvester needs `n >= 32`" stands.
* **The headline "the cheap-span form of ZE-TC is false" is right, but the write-up's witness is weak.** The witness is a zonotopal, i.e. unfaithful, model (C1).
* **New finding: the refutation strengthens to faithful models.** By monotonicity of `nu_A` under inclusion together with the one-cut formula:
  * ZE-TC^cheap already fails in the zonotope-free model `M8` at `m = 8` (exact check);
  * it fails in the LP-supercritical, zonotope-free, overdetermined, twin-free product `M8 x M8` at `m = 16`;
  * so README §8.1's "ZE-TC holds in M8 (J4)" is true only in the deg form.

  All failing points are degenerate, and the regime case remains undecided.
* **The deg-type results do not supply the near-strata step** (C2). These are the Sylvester tower, the products and the zonotopal models. They give neither height control nor transversal order. Only the rat form (one cut, disjoint cuts) does.
* There are several quantifier slips ("every order" must be "every `m <= rho`", C4), numerical slips in the section count (C6), and summary-level overstatements about the Paley sections (C7).


## 1. Definitions and Lemma 1.1 (TC = {nu < m})

Re-derived all seven parts. They are correct, with two quantifier remarks.

* (1) Level sets: `J` is an ideal, so `l_u^i in J` gives `l_u^(i+1) in J`. `TC(A,m) = {nu_A < m}` is compute.md R2 restated. Correct.
* (2) Semicontinuity: correct (finite union of complements of hypersurfaces `Phi_c^(i) = 0`).
* (3) Subadditivity: binomial expansion. Correct.
* (4) Projections: for `i > rho(pi(A))`, `J(pi(A))_i` is all forms of degree `i`, so `(sum t_j y_j)^i` agrees on `pi(A)` with a polynomial of lower degree. Correct. The bound `nu(tv) <= h_A(v)` and "cheap span => tangent-cheap" follow with (3).
* (5) Products: Buchberger's first criterion (coprime leading monomials) plus degree-compatibility, then the bigraded tensor argument. Correct over `C`.
* (6) Zonotope lower bound and (7) transversal order: correct. In (7), `l_{u0}^(m-i) = g` on `A` with `deg g < m-i` because `F in J_j` means `F` agrees on `A` with `F - f` for some `f in I(A)` with top form `F`.

**Quantifier remark Q1 (minor).** The chain `cheap => rat => deg` needs `V(A,m) != 0`, i.e. `m <= rho(A)`. For `m > rho(A)`, `TC(A,m) = C^d`, the whole space `M = R^d` is "tangent-cheap" (`nu_A(R^d) = rho < m`), so ZE-TC^rat holds vacuously and ZE-TC^deg is false. Several later statements say "at every order" (Cor 2.2(1), §3.2 "ZE-TC^deg holds ... at every order", Lemma 3.2, Prop 3.3 Remark (1)). These must read "every `m <= rho(A)`". Taken literally, the hypothesis of Lemma 3.2's deg-clause, "ZE-TC^deg for `A_1` at every order", is never satisfied (for `d_1 >= 1`). The proof works verbatim for `m <= rho(A)` once the factor hypotheses are restricted to orders `<= rho(A_i)`: if `a >= rho(A_1)` then `b <= rho(A_2) - 1` and the second factor is degenerate. Nothing downstream depends on the literal form.

**Remark Q2 (major, about what the proved statements deliver).** §0 item 1 calls ZE-TC^deg "what the tangent-cone reduction needs". §1.1 says ZE-TC^rat is what the near-strata step uses, and deg is "the weakest statement compatible with" the reduction. The second reading is the accurate one. ZE-TC^deg puts real zeros on *some* rational hyperplane. It controls neither the height `|v|_1` of the normal vector, which the character/pair inequality needs (`X_v >= W/2 - 2 log(C|v|_1/...)`), nor the transversal order `m - nu_A(M)` (which can be `<= 0` on a degenerate stratum `M` with `nu_A(M) >= m`). So the deg-type results proved in §3.4-§4 and §5.4 (Sylvester tower, products, zonotopal models) do not feed the near-strata step of README §8.2 as stated. Only the one-cut and disjoint-cut results are rat-type.

**Independent exact check.** `referee_zetc_checks/rlib.py` builds `V(A,i)` exactly as integer nullspaces (`fmpz_mat`) and computes `nu_A(u)` from its definition.


## 2. Theorem 2.1 (one cut) and the ZE-TC^cheap counterexample

**Theorem 2.1: statement correct; one step of the proof is wrong as written and easily repaired.**

* *Error in the slice lemma.* The proof asserts that `prod_{x in Z}(1 - k_x^2) prod_{x in N}(k_x^2 - k_x)/2` has degree `|Z| + 2|N| = s - j`. Its degree is `2|Z| + 2|N|`, which exceeds `s - j` whenever `Z` is nonempty. So the bound "every function on `Lambda_j` has degree `<= s - |j|`" is not proved as written.
* *Repair (checked).* Use `prod_{x in Z}(1 - k_x) prod_{x in N}(k_x^2 - k_x)/2` (for `j >= 0`). It has degree `|Z| + 2|N| = s - j`. It forces `k_x = -1` on `N` and `k_x in {-1, 0}` on `Z`. Then `sum_{Z u P} k = j + |N| = |P|` forces `k = +1` on `P` and `k = 0` on `Z`, so it is the indicator of `kappa` on `Lambda_j`. `cor22_slices.py` verifies this indicator property for every `kappa` and every `j`, `s <= 6`. It also verifies exactly that `rho(Lambda_j) = s - |j|` for all `s <= 7`, `0 <= j <= s`.
* *Newton gluing and lower bound.* Re-derived. The gadget table is correct: I listed the 8 lattice points of the triangle and of the hexpair by hand. The case analysis is also correct, including the excess count (`<= L`), the multipartite matching, and the odd case (`>= 3` classes, with no fourth class of maximal size, since four classes of size `(s-L-1)/2` exceed `s - L` once `s - L >= 3`).
* *Independent exact check of the formula* (`t21_check.py`, exact `fmpz` nullspaces, not mod `p`). For every `1 <= s <= 5` and `0 <= L < s` (excluding `(1,0)`), `nu_S(u)` equals `min(2 sigma, 2L + 2t, s + L)` for 8237 integer vectors `u`. These are all of `{-2..3}^s` for `s <= 4`, plus random vectors including additive relations such as `u = (1,2,3,...)`. There were 0 mismatches, and `rho(A_S(L)) = s + L` in every case. The author's check was mod `p`; this one is exact.

**Corollary 2.2: (1) and (2) are correct; the reading of (3) is overstated.**

* (1) The three families of subspaces are right, and ZE-TC^rat holds for one cut (for `m <= rho`; see Q1).
* (2) Re-derived. The smallest rational subspace containing `u` is `R^S x R^Y` for `u_S` generic. Every integer direction has width `>= 2` because `+-e_x in A`. The cost is therefore `>= 2(s+q) >= m`. Correct.
* (3) Exact check (`cor22_slices.py`): for `S = {0,1}`, `L = 1`, `d = 6`, I get `|A| = 567`, `rho(A) = 11`, Hilbert function `[1,7,28,77,162,274,389,480,535,559,566,567]`, `width k(S) = 2` (so `h_A(beta) = 4 = m/2` at `m = 8`), and `nu_A(u) = 7` for `u = (1,2,3,5,0,0), (2,7,-3,11,0,0), (1,2,1,1,0,0), (1,-1,3,5,0,0)`. So `u` lies in `TC(A,8)` and is in no cheap span. **The counterexample is valid as a statement about this `(A, m)`.**

**Issue C1 (major, scope of the headline and of its justification).** The write-up's counterexample model is **zonotopal**: `zeta(A) = rho(A) = 11 >= m = 8` (Theorem 2.1 plus product). By compute.md Props 2.1-2.2 it is "unfaithful": ZE holds there through one character, and no actual regime cluster has such a box. README §8.2 states ZE-TC for the `T = 1` slab-box of an *admissible profile* at `m = ceil(4w)`, which is zonotope-free. So the write-up's own witness does not show that the README's form fails in any model the README cared about.

**A stronger refutation, missed by the author (referee, proved and checked exactly).** The same mechanism works in zonotope-free models, through monotonicity under inclusion: `A subset A'` implies `I(A') subset I(A)`, hence `J(A') subset J(A)` and `nu_A <= nu_(A')`.

* **`M8 = Syl(8,2)`, `m = 8`** (zonotope-free with `zeta <= 7`; 42 tangent forms in 7 variables; twin-free; LP-subcritical).
  * `M8` lies in the one-cut box of any of its columns `S` (`|S| = 4`, `L = 2`).
  * By Theorem 2.1, every `u in R^S` has `nu_(M8)(u) <= s + L = 6 < 8`, so `R^S` lies in `TC(M8, 8)`.
  * Take `u_S = (1, sqrt 2, sqrt 3, sqrt 5)`. Its rational hull is `R^S`, of dimension 4. Since `+-e_x in M8`, every nonzero integer direction has width `>= 2`, so every family whose span contains `u` costs `>= 8 = m`.
  * Hence **ZE-TC^cheap(M8, 8) is false.** Exact check (`m8_cheap.py`): `rho(M8) = 8`, Hilbert function `[1,8,36,113,274,533,820,995,1037]` (matching compute.md), `nu = 6` for three integer `u` supported on a column, and `nu = 8` at a generic point.
  * This also shows that README §8.1's sentence "In the two zonotope-free `T = 1` models it holds (J4)" is true only in the **deg** form. J4 places `TC(M8, 8)` in rational hyperplanes, not in cheap spans.
* **`M8 x M8`, `m = 16`** (LP-supercritical with `h = 8 = m/2`; zonotope-free with `zeta <= 14`; `1764` tangent forms; twin-free; the write-up's own §5.4 model).
  * Take `u = (u_1, u_2)` with `u_i` generic in `R^(S_i)` for a column `S_i` of block `i`.
  * By Lemma 1.1(5), `nu(u) = nu_1(u_1) + nu_2(u_2) <= 12 < 16`.
  * The rational hull has dimension 8, so every spanning family costs `>= 16 = m`.
  * So **ZE-TC^cheap fails in a model that meets every criterion the round-11 README and the task state** (LP-supercritical, zonotope-free, overdetermined, twin-free). The model is decomposable and its columns are unbalanced.
* In all these examples the failing points are **degenerate**: they vanish outside `S` (`S_1 u S_2`). ZE-TC^rat and ZE-TC^deg are untouched.

So the headline "ZE-TC^cheap is false" is **true in a stronger sense than the write-up proves**: it fails in faithful models, both LP-subcritical (`M8`) and LP-supercritical (`M8 x M8`). It is still **not shown for a regime box at `m = ceil(4w)`**. There, every column has `|S| + L >> m`, so this one-cut mechanism produces nothing outside cheap spans. The *conceptual* conclusion stands and is the main contribution of §1-2: the right cost of a rational subspace `M` is `nu_A(M)`, not a width sum, and Lemma 1.1(7) gives transversal order `m - nu_A(M)` along every rational `M`.

## 3. Zonotopes, Prop 3.1, the few-cut claims, products and sections

* **Prop 3.1** is correct. It is compute.md Prop 2.1 read at the lowest order: `Phi_c = m! prod <v_i,u>^(alpha_i)`. The proviso of Q1 applies: "ZE-TC^deg at every order" means every `m <= zeta(A)`.
* **Two crossing or nested cuts, `|S_1 u S_2| <= 7`: confirmed independently.** I did not reuse the author's zonotopes, which are not in the output file. `zonocheck.py`/`twocut_ref.py` compute `rho_p` from the mod-`(2^61-1)` rank profile, which is a valid upper bound for `rho_Q`. Separately, an integer program over edge multiplicities (edges in `{-1,0,1}^d`) finds a zonotope, and every one of its lattice points is checked to lie in `A`. I ran all `(a,b,c)`, `b >= 1`, `a+b+c <= 7`, `1 <= L_1 < a+b`, `1 <= L_2 < b+c`, **without** the `P <-> R` reduction. Results: `twocut_ref_out.txt` (754 blocks, run stopped by my time limit after every type except `(3,1,3)`) and `threecut_ref2_out.txt` (the 9 blocks of type `(3,1,3)`). There is no gap: the zonotope length equals `rho_p` in every block, so `zeta = rho` exactly. Every block type up to mirror symmetry is covered.
* **Three cuts, block size `<= 6`, `L_j in {1,2}`.** I enumerated block types myself: three distinct columns of size `>= 2` with union the block, up to permutations of coordinates and columns. That gives **843** types (3, 40, 197, 603 at block sizes 3-6). The author reports 953 and probably uses a weaker identification, which does not matter for a "for all" claim. Result (`threecut_ref2_out.txt`): **843/843 types have a verified zonotope of length `rho_p`, with no gap.** The three-cut claim is confirmed independently on a complete list of types.
* *Logic of the certificate:* `zeta <= rho_Q <= rho_p = (verified zonotope length) <= zeta`. This is correct.
* **Reading of §3.5.** "An induction on cuts that propagates zonotopes cannot reach the regime" is a correct meta-observation (compute.md Prop 2.2). The only evidence for "zonotope-free needs many cuts" is the survey (Z6 has 6 distinct columns on `d = 6`), and it is labelled as evidence.
* **Z6 (independent, exact; `z6_ref.py`).** `|A| = 247`, cumulative Hilbert function `[1,7,28,75,147,214,244,247]` (the author's increments), `rho = 7`, `dim V(A,7) = 3`, twin-free with 6 distinct columns. A verified zonotope of length 6 (found by my own integer program) and an exact LP certificate `zeta <= 13/2` give `zeta = 6`. Confirmed. Two remarks: the column `{0..5}` (the base-point box `|k_n| <= L`) is listed twice, with `L = 2` and `L = 1`, so the `L = 2` copy is redundant; and the model is underdetermined, as the author says.
* **Plane test (§3.2).** The logic is right. A real irrational root `t` of an irreducible gcd factor gives `u` with rational hull `M` and `nu(u) < nu(M)`, a counterexample to ZE-TC^rat at order `nu(M)`. A gcd with only linear factors over `Q` excludes this on that plane. It is exact per plane and evidence for the statement. I did not rerun it.

### Products (Lemma 3.2) and sections (Prop 3.3)

* Lemma 3.2 is correct, with Q1. The LP-additivity step is correct: if `(v_1, v_2)` is primitive and `v_1 != 0`, then `v_1` is a nonzero integer multiple of a primitive difference direction of `A_1`.
* Prop 3.3 is correct: pull back `l_u^i = g` along `t -> k_0 + Vt` and expand. The degeneracy transfer `w -> Vw` uses `rank V = j`. Correct.
* **Remark (3) of §3.4 has two numerical slips. The qualitative conclusion is unaffected, and it is labelled "a count".**
  1. A section supported on a point set `Y` is slab-slack only when `|Y| <= w/2`, not `|Y| <= w`, because `|<k,beta>| = 2|k(S)| <= 2|Y|`. The coordinate certificates therefore stop near `nu ~ w ~ m/4`, not `2w ~ m/2`.
  2. The per-block usage bound `4||c||_1 <= 4 sqrt 8 ~ 11.4` is a crude upper bound, and using it to say "this allows about `w/11.4` blocks" reads an upper bound on usage as if it were the usage itself. Exact computation (`m8usage.py`) on the 8-point Hadamard block with all eight coordinates boxed (`|B| = 869`, compute.md's `sym8`) gives the maximum of `|<k,beta>|` over the 254 nonconstant sign vectors `beta`. The distribution is `{2: 16, 4: 70, 6: 112, 8: 56}`: maximum 8, mean 5.64. So the guaranteed number of blocks is `w/8` (total top degree `w`, short by a factor 4, not 5.7), and with typical usage about `w/5.6` (short by about 2.8). The claim "M8 is a factor ~5.7 short" should read "a factor 2.8-4 short". This still supports the conclusion that M8-sections cannot reach `m >= 4w`. (ZE-TC for the boxed variant `sym8` is itself not established: compute.md §5.2 left its saturation unfinished. The author notes this.)


## 4. Sylvester tower (Theorem 4.2), symmetry, ball effect

* **Lemma 4.1** is correct. For `k` supported on a subgroup `U`, `hat k(a)` depends only on `chi_a|_U`. Characters trivial on `U` give `hat k = sum k = 0`, and every nontrivial character of `U` occurs.
* **Theorem 4.2** is a correct consequence of Prop 3.3 and Lemma 4.1, **given J4**.
  * *J4's status.* J4 is an **unrefereed** judge finding (README §5 heading). It is not refereed, as the write-up's wording ("rests on the refereed compute.md Thm 5.2") might suggest. I re-derived J4 from the refereed data of compute.md Thm 5.2(b).
    * For the 24 form-type basis polynomials, the orders of the two-term factors and the forms add up to `m` exactly: `2 + 3 + 3` for 11 of them and `4 + 4` for 13. So the lowest-order part of `P_c(e^z)` is `Phi_c/m!` with no cancellation.
    * The saturation of the ideal of the 24 non-binomial parts by the collision forms is the unit ideal.
    * Hence `TC(M8, 8)` lies in the collision hyperplanes together with the 27 two-term tangent hyperplanes, all rational.
    * The same argument works for `M8'` (7 forms, orders adding up to 11).
    * The deduction is sound. I did not recompute the saturation, which the judge's own attempt also abandoned. The theorem should therefore be cited as "conditional on J4, which follows from the refereed Thm 5.2 data".
  * The zonotope-freeness for `L = 2` (`zeta <= 8(n-1)/n < 8`) and `V(A,m) != 0` (from `nu_A >= nu_(A')`) are correct.
  * As the author says, the statement is **LP-subcritical** (`m = 8 < 16`, `11 < 24`) and says nothing at regime orders. It is a deg-statement only (Q2): no height bound and no transversal order.
* **Lemma 4.3** is correct and trivial. A fixed real point has two equal coordinates, or a coordinate equal to `u_0 = 0`.
* **Lemma 4.4** is correct. `H^T H = n I` and `<k, 1> = 0` give `n|k|^2 <= (n-1) w^2`. Here `|k|^2` includes the base coordinate, so "fewer than `4L^2` nonzero entries" counts the base point as well.
* *Normalisation of Hadamard models (minor).* A `T = 1` Hadamard model depends on two choices: the row normalised to all-ones, and the base point. The base point is unboxed and is the reference for the patterns. The stated reason, "automorphism groups are transitive on points", is not quite enough. What is needed is that `Aut(H)` is transitive on rows and that the stabiliser of the normalising row is transitive on points. Both hold for `H_8` (`AGL(3,2)`) and for `H_12` (the row stabiliser acts as `M_11` on 12 points), so the conclusion stands.



## 5. Small-model computations (§5 of the write-up, item (c) of the summary)

**Direction of the mod-`p` logic.** A mod-`p` rank is at most the rational rank. So "`H_p(D) = |A|`" implies "`H_Q(D) = |A|`", i.e. `rho_Q <= D`. In the other direction, mod-`p` data give only evidence. The write-up uses mod `p` only in the valid direction. One caveat: the spanning certificates computed through isotypic components use roots of unity in `F_p`, so they are genuine `F_p`-linear algebra, not reductions of an integer computation. The implication still holds, because the spanning statement is about the integer evaluation matrix, whose `F_p`-rank is basis-independent.

| claim | referee's independent check | result |
| --- | --- | --- |
| `H_8`: `|A| = 57/1037/1877` (`L = 1/2/3`), `rho(L=1) = 3` | `small_had_ref.py` (exact `fmpz` ranks; own Sylvester construction) | confirmed (`H = [1,8,36,57]`) |
| Paley-12 `L = 1`: `|A| = 133`, `rho = 3` | own Paley-I construction, exact | confirmed (`[1,12,78,133]`) |
| Paley-12 `L = 2`: `|A| = 20153`, `rho <= 8` | `cyc_rank.py` + `paley_ref.py` (`p = 1048433`), and `inc_cert.py` + `paley_inc.py` (`p = 699953`): `Z_11 x {+-1}` isotypic spanning test, all 22 components of full rank, ranks summing to 20153 | confirmed. The per-degree gains give `H(8) - H(7) = 11 * 395 = 4345` (mod `p`), matching the author's 4345 |
| Syl(16,1): `|A| = 241`, `rho = 3` | exact | confirmed (`[1,16,136,241]`) |
| Syl(16,2): `|A| = 213937`, `rho <= 15` | `syl16_inc.py`: own `F_16` model, Singer cycle `Z_15 x {+-1}`, incremental mod-`p` echelon | **confirmed, and sharpened to `rho <= 11`.** The 10 Frobenius-class representatives `j in {0,1,3,5,7}` x `eps in {+,-}` all reach full rank using monomials of degree `<= 11` (`p = 699961`). The Frobenius `x -> x^2` maps `A` to itself (checked) and conjugates `sigma` to `sigma^2`, so it equates the ranks within the classes `{0}, {1,2,4,8}, {3,6,12,9}, {5,10}, {7,14,13,11}`. The row counts, weighted by class size, sum to `106969 + 106968 = 213937 = |A|`. Degree gains of the trivial even component: `[1, 8, 189, 1877, 4962, 100]` at degrees 0-10 |
| Syl(16,3): `|A| = 3608377`, `rho <= 22` (fibration) | `fib_ref.py`: own exact Newton degrees on the 3-column fibration (ordering by `|beta|_1`) | confirmed: `max(e_i + delta_i) = 16`, hence `rho <= 22`. `L = 2` gives 21 |
| Syl(16, `L >= 4`): one-cut bound `22 + L < 8L` | arithmetic | confirmed |
| Lemma 5.1 (fibration bound) | re-derived (Newton form over the fibres) | correct |

Theorem 2.1's slice lemma enters the fibration bound. My exact slice check (`rho(Lambda_j) = s - |j|`, `s <= 7`) covers the case `s = 4` used there, so the fibration bound is independent of the proof slip found in §2.

**LP-supercriticality thresholds.** For a Hadamard `T = 1` box the slab `|k(S)| <= L` is attained whenever `L <= |S|` (take `+1` on `L` points of `S`). So LP-supercritical means `m >= 8L`, as the write-up uses. Each table row checks `rho < 8L` against the right threshold.

**Paley-12 sections (§5.5).** The author's sign test (`section_sign.py`) is **floating point**: sampled 8th powers in `float64`. So "indefinite" was evidence, although the text treats it as fact ("hence have nondegenerate real tangent-cone points"). My checks are exact.
* `paley_section_sign.py`, on the author's three classes: `|A_Y| = 909, 901, 915`, `dim V(A_Y, 8) = 1`, top form `= -8!` times an octic that is **irreducible over `Q`** (sympy `factor_list`). Exact integer evaluation finds points of both signs.
* `paley_sections_all.py`, on **all 30 translation classes** (330 sections): exact `dim V(A_Y, 8) = 1` and `rho(A_Y) = 8` in every class. The top form is **irreducible over `Q`** in 30/30 and **indefinite by exact integer evaluation** in 30/30. So every 8-point section has nondegenerate real tangent-cone points. The summary claim is now proved; the write-up had shown it for only 3 classes, with a floating-point sign test.
* Irreducible over `Q` plus indefinite does imply nondegenerate real zeros. The real zero set of an indefinite form has dimension `d - 1`. If countably many rational hyperplanes covered it, one of them would contain a Zariski-dense piece of a `(d-1)`-dimensional component, and its rational linear form would divide the octic.
* The numerical search for common zeros of the 330 octics is correctly labelled evidence. I did not rerun it.

**Summary-level overstatements.** The summary says "each of its 8-point sections has nondegenerate real zeros". The body checked only 3 of the 30 classes, with a floating-point sign test. See the line above for what the exact all-class check gives.

**Products (§5.4).** `M8 x M8`: `rho = 16` and `H_A(16) = 42^2 = 1764` follow from `S/J = S_1/J_1 (x) S_2/J_2` (Lemma 1.1(5)). LP-supercriticality (`h = 8 = m/2`), zonotope-freeness (`zeta <= 7 + 7`) and twin-freeness are correct. The point that "zonotope-free + LP-supercritical + overdetermined" must be supplemented by indecomposability is correct and useful. A product of TC-good pieces passes every stated criterion while carrying no new information.


## 7. Issues list

| # | severity | location | issue | fix |
| --- | --- | --- | --- | --- |
| C1 | major | §0 item 1, §2 Cor 2.2(3), summary | The witness refuting ZE-TC^cheap is a **zonotopal** (unfaithful) model, so it does not show that the README's form fails in any faithful model. | Use the referee's witnesses: `M8` at `m = 8` (zonotope-free) and `M8 x M8` at `m = 16` (zonotope-free, LP-supercritical, overdetermined, twin-free). Both use inclusion monotonicity plus the one-cut formula. State the failure as "in faithful but non-regime models"; the regime remains undecided. |
| C2 | major | §0, §1.1, §4 | ZE-TC^deg is called "what the tangent-cone reduction needs". It controls neither the height of the normal vectors nor the transversal order, so the deg-results (Sylvester tower, products, zonotopal models) do not feed the near-strata step. Only one cut and disjoint cuts are proved in the rat form. | Say "ZE-TC^rat (with height control) is what the near-strata step needs; deg is a necessary condition". |
| C3 | minor (proof gap, statement true) | Thm 2.1 proof, slice lemma | `prod_Z (1 - k_x^2)` has degree `2|Z|`, so the stated degree count `|Z| + 2|N|` is false. | Replace it by `prod_Z (1 - k_x)`. Referee-verified: indicator property for `s <= 6`, exact slice degrees for `s <= 7`. |
| C4 | minor | Cor 2.2(1), §3.2, Lemma 3.2, Prop 3.3 Rmk (1) | "at every order" is wrong for `m > rho` (then `TC = C^d`, and ZE-TC^deg is false while ZE-TC^rat holds vacuously). As literally stated, the deg-hypothesis of Lemma 3.2 is unsatisfiable. | Read "every `m <= rho(A)`" (factors: every order `<= rho(A_i)`). The proofs go through. |
| C5 | minor | §4.1 Thm 4.2 | It rests on J4, which is unrefereed. Its deduction from the refereed Thm 5.2 data checks out (orders add up to `m`; the saturation is by collision forms). | Cite it as conditional on J4. |
| C6 | minor | §3.4 Rmk (3) | The slack-section bound is `|Y| <= w/2`, not `w`. The per-block usage `11.4` is a crude upper bound; the exact maximum is 8 and the mean 5.64. "M8 is a factor ~5.7 short" should read "~2.8-4". | Correct the numbers; the conclusion is unchanged. |
| C7 | minor | §5.5, summary | "Indefinite" rested on a floating-point sign test (`section_sign.py`). The summary's "each of its 8-point sections has nondegenerate real zeros" was claimed for all sections but checked only for 3 of 30 classes. | Now exact for all classes (§5 above). |
| C8 | minor | §5.3 | Uniqueness of the Hadamard models needs transitivity of the row-stabiliser on points, not only transitivity on points. | It holds for `H_8` and `H_12`; add the remark. |
| C9 | cosmetic | §5.2 Z6 | The column `{0..5}` is listed twice (`L = 2` and `L = 1`), so the `L = 2` copy is redundant. | Drop it. |

No error was found in Lemma 1.1, Theorem 2.1 (statement), Cor 2.2(1)-(2), Prop 3.1, Lemma 3.2 (with C4), Prop 3.3, Lemmas 4.1, 4.3, 4.4, Lemma 5.1 or the fibration bound. All the computer-assisted exact claims I re-ran agree.


## 8. Corrected statement

All of the following concern `T = 1` slab-boxes `A` in the ternary cube, cut by `|k(S_p)| <= L_p` (`L_p >= 1` unless stated), at orders `1 <= m <= rho(A)`. "Degenerate" means lying on a rational hyperplane.

1. **(Proved.)** `TC(A,m) = {u : nu_A(u) < m}`. The function `nu_A` is subadditive, additive on products, monotone under inclusion (`A subset A'` gives `nu_A <= nu_(A')`), bounded above by `rho` of projections and below by zonotopes. Along every rational subspace `M` with `nu_A(M) < m`, every `Phi_c` vanishes transversally to order `>= m - nu_A(M)`. For `m <= rho(A)`: ZE-TC^cheap implies ZE-TC^rat, which implies ZE-TC^deg.
2. **(Proved.)** One cut, `0 <= L <= s - 1`: `nu_S(u) = min(2 sigma(u), 2L + 2t(u), s + L)`. The slice step needs the corrected indicator `prod_Z (1 - k_x)`. For one cut and for disjoint cuts, ZE-TC^rat holds at every `m <= rho`.
3. **(Proved, strengthened by the referee.)** ZE-TC^cheap (README §8.2) is false for general `T = 1` models:
   * in zonotopal models (one cut, `|S| + L < m <= rho`);
   * in the zonotope-free model `M8` at `m = 8`;
   * in the zonotope-free, LP-supercritical, overdetermined, twin-free product `M8 x M8` at `m = 16`.

   In every case the failing points are degenerate. Whether ZE-TC^cheap fails for a regime box at `m = ceil(4w)` is **open**; the one-cut mechanism does not apply there.
4. **(Proved / computer-assisted, exact in the direction used; independently re-verified.)** `zeta = rho` holds, and hence ZE-TC^deg holds for every `m <= rho`, in the following models:
   * two-cut models with `|S_1 u S_2| <= 7`, `L_j >= 1`;
   * three-cut models with block size `<= 6`, `L_j in {1,2}` (843 isomorphism types).
5. **(Proved; deg form only, conditional on the unrefereed but sound J4.)** ZE-TC^deg holds:
   * for `Syl(2^k, 2)` at `m = 8` and `Syl(2^k, 3)` at `m = 11`, every `k >= 3` (LP-subcritical);
   * for products of pieces that are TC-good at their top degree, at every `m <= rho` (e.g. `M8 x M8` at `m = 16`).

   These statements give no height bound and no transversal order. They therefore do not supply the near-strata step.
6. **(Rigorous upper bounds on `rho`.)**
   * No Hadamard `T = 1` model with `n <= 12` and one `L` is LP-supercritical with `V(A,m) != 0`.
   * `Syl(16, L)`: `rho <= 22 < 24` at `L = 3`, and the one-cut bound covers `L >= 4`. `L = 2`: `rho <= 11 < 16` (referee's spanning certificate; the author's `<= 15` is also correct). Hence the smallest LP-supercritical Sylvester `T = 1` model has `n >= 32`. Confirmed.
7. **(Exact.)** Paley-12, `L = 2`, `m = 8`: `rho = 8`, `H(8) - H(7) = 4345` (mod `p`), zonotope-free, LP-subcritical. All 330 eight-point sections have `dim V = 1` and an irreducible, indefinite top octic (exact). Each therefore has nondegenerate real tangent-cone points, so Prop 3.3 cannot certify ZE-TC from any single section. The claim that the full system has no nondegenerate real zero is **evidence only**.
8. **Not proved:** ZE-TC in any form for any LP-supercritical *indecomposable* zonotope-free model; ZE-TC in the regime; any growth improvement or uniform bound. `M <= (2/3 + eps) log R/loglog R` remains the best bound.


## 9. Referee scripts (`referee_zetc_checks/`; no author code imported)

| script | content | status |
| --- | --- | --- |
| `rlib.py` | exact `V(A,i)` (integer nullspaces) and `nu_A(u)` | library |
| `t21_check.py` | Theorem 2.1, `s <= 5`, all `L` (8237 cases) | exact |
| `cor22_slices.py` | slice degrees `s <= 7`, corrected indicator, Cor 2.2(3) | exact |
| `m8_cheap.py` | `nu = 6` on column spans in `M8`; ZE-TC^cheap fails in `M8` | exact |
| `zonocheck.py`, `twocut_ref.py`, `threecut_ref2.py`, `tc_enum.py` | `rho_p` and verified zonotopes, two cuts (all types) and three cuts (843 types) | exact in the direction used |
| `z6_ref.py` | Z6 | exact |
| `small_had_ref.py` | `H_8`, Paley-12 and `Syl(16)` counts; `rho` at `L = 1` | exact |
| `fib_ref.py` | fibration bound, `Syl(16,2/3)` | exact |
| `cyc_rank.py`, `paley_ref.py`, `inc_rank.py`, `inc_cert.py`, `paley_inc.py`, `syl16_inc.py`, `syl16_inc2.py` | isotypic spanning certificates (Paley-12 `L = 2`: `rho <= 8`; `Syl(16,2)`: `rho <= 11`) | mod `p`, rigorous upper bounds |
| `paley_section_sign.py`, `paley_sections_all.py` | all 30 Paley-12 section classes: `dim V = 1`, irreducible, indefinite | exact |
| `m8usage.py` | slab usage of an 8-point Hadamard block | exact |
