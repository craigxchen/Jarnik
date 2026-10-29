# Round 27: adversarial review (killer) of lane interp (`round27/interp.md`)

Status: complete. Killer: independent re-derivation and recomputation; nothing imported from `interp_checks/`.
Scripts: `round27/kill_interp_checks/` (absolute path in the Appendix).

## 0. Verdict

**The mathematics survives. The novelty framing and the model evidence do not.** No growth or uniformity is claimed,
and none follows. `M <= (2/3+eps) log R/loglog R` stands.

1. **Correct (re-derived):**
   * Lemma 3.0;
   * Theorem 3.1 (1)–(4);
   * the equivalences (3.2);
   * Corollary 4.1 as an implication;
   * Lemma 2.1, Theorem 2.2, Props 3.2, 3.4, 3.5.

   I found no mathematical error in any proof.
2. **Model claims that are false** (§5.2):
   * *Zonotope order.* "In every computable model it reaches exactly the zonotope order (6 for sdc7 and syl8)" is
     wrong. Both sdc7 and syl8 (`L = 2`) contain **order-7 lattice zonotopes with all generators in `L_pair`**,
     verified exactly. So one product form certifies order 7 under the criterion's own hypothesis, while the
     criterion stops at 6.
   * *Shift sets.* "The failure is genuine" (independent of the shift set) is false for syl8: with all
     `{-1,0,1}^7` shifts, `(SHIFT_6)` holds and the reach is 7.
3. **A favourable fact the write-up missed** (§5.2(d)). On M8' (syl8, `L = 3`), `(SHIFT_e)` w.r.t. `L_pair` holds
   for all `e <= 10` (mod p, two primes). So the criterion proves ZE-braid at `rho = 11` by linear algebra alone,
   reproducing round 11's Gröbner theorem. The exactly verified zonotope order there is 10. This is the only slab box
   I found where the criterion beats zonotopes. In 27 random multi-cut slab boxes it never does: below `zeta` in 20,
   equal in 7, with either shift set.
4. **Novelty corrections** (§7):
   * With `L_pair`, the conclusion of Theorem 3.1 **is round 11's "ZE-braid" target** (compute §5.3), which the
     write-up never cites.
   * Prop 3.4 is the round-11 referee's ZE-braid refutation (`referee_compute.md` §2.8) in its `T = 1` form.
   * Theorem 2.2 is the record's zonotope cap (rounds 11, 13, 16) plus a standard classification. Its moral is in
     round 11 compute §5.2 Remarks.
   * "The kind round 13/16 named as the missing tool" misreads round 13, where the phrase concerns cap-side
     certificates.
   * Non-zonotopal exact certification on slab models already existed (round 11 Thm 5.2, Gröbner saturation).
5. **Tests.** T1 passes (rigidity is fake-violated). T2 and T3 are not applicable, and I agree. T4 passes formally.
   The models are mixed (§6): mostly at or below the Liouville-capped zonotope class, with one exception (M8'). T5
   passes.
6. **Novelty verdict: genuinely new, modestly.** Theorem 3.1 is an elementary sufficient criterion, new to the record
   and not found in the literature after three searches; it may be folklore. It aims at a mapped target (ZE-braid_H /
   ZE^all_H), and Corollary 4.1 reduces route 5 to a hypothesis at least as hard as the record's open lemma.
   Whether `(SHIFT)_A` holds on regime boxes is open. The models are LP-subcritical and point both ways.

## 1. Re-derivation of Lemma 3.0 and Theorem 3.1

**Lemma 3.0: correct.** I checked three points. First, `supp T_l c` lies in `S_l + l`, which is contained in `S`.
Second, for `|alpha| <= e` the difference `(k+l)^alpha - k^alpha` is a polynomial of degree `<= |alpha| - 1 <= e - 1` in
`k`, so it is annihilated by `c in V(S_l, e)`. Third, both identities in (3.1) are direct expansions; the terms
`i >= 2` of `<k+l,u>^{e+1} - <k,u>^{e+1}` are moments of `c` of order `<= e - 1`.

**Equivalences in (3.2): correct.** Since `V(A_l, e)` is contained in `V(A, e)`, the kernel of
`c -> Phi^{(e)}_c` on `V(A,e)` is `V(A,e+1)` (over a field of characteristic 0, the powers `<k,u>^e` span the
degree-`e` forms), and `M_e = J_e^perp`. Hence

```text
V(A,e) = sum_l V(A_l,e) + V(A,e+1)   <=>   M_e(A) = sum_l M_e(A_l)   <=>   J(A)_e = cap_l J(A_l)_e .
```

**Theorem 3.1 (1)-(4): correct as stated.**
* (1) The induction is sound. At each step one term of the finite decomposition is nonzero at `Y`, and either it
  already lies in `V(A,e+1)`, or `(Y^l - 1) != 0` carries it there.
* (2) The same induction works with `Phi^{(e)}`.
* (3), (4) If `(SHIFT_e)` holds over `Q`, every decomposition in the induction can be chosen over `Q`: the base case
  `delta_k` is rational, and the shift preserves `Q`. So part (1) for `m - 1` yields a **rational**
  `c in V(A, m-1)` with `P_c(Y) != 0`, and (3.3) applies to it.
* In (3), `T_l g - g` is an integer vector of `V(A,m)` with `|.|_inf <= 2(H/2) = H`, so it lies in the certified
  set. In (4), the `h_i` terms are killed by `Phi^{(m-1)}`.

I found no gap. The theorem is elementary: it is the induction "if every element of `V(A,e)` splits along the
shift-intersections, then multiplication by `(Y^l - 1)` raises the order without killing the value at `Y`".

**Corollary 4.1 (the implication): correct.**
* `Z^{e_x} = z_x/z_n` and `Z^{e_x - e_y} = z_x/z_y`, so for distinct points `Z^l != 1` for every `l in L_pair`.
* Theorem 3.1(3) produces an integer `c in V(A, m')` with `|c|_inf <= U^2` and `P_c(Z) != 0`.
* `P_c(Z) = N^{-D/2} F_c(z)`. Rigidity (t1rigour Thm 3.1, which I re-read: Liouville plus forced divisibility, any
  norm, cost `<= wW <= m'W/4`) forces `F_c(z) = 0`. This is a contradiction.
* The size statement uses `r <= W/ln 5` (split primes are `>= 5`) and `W >= 3n ln n - O_C(n)` (two_thirds Thm 7.5).
  Both are right.

**What the corollary is.** It is the composition of three arrows:

```text
(SHIFT)_A  ==(Thm 3.1(3))==>  no L_pair-nondegenerate Y in Bs(V_H(A,m'))  ==>  ZE at Z  ==(rigidity; rounds 12, 26)==>  n - 1 < 316 ln 4r.
```

* Only the first arrow is new.
* The middle statement is round 11's **ZE-braid** target in height-cut form (round 11 compute §5.3; see §7). It is
  **strictly stronger** than the record's open lemma ZE^all_H (round 26 §1.1 item 4).
  ZE^all_H excludes only points with `Y^v != 1` for all `v != 0`. The middle statement excludes every point off the
  `|L_pair|` hypersurfaces `{Y^l = 1}`, including points on other proper subtori.
* So Corollary 4.1 reduces route 5 to a hypothesis that implies the record's open lemma in a stronger form. It is
  not a reformulation of that lemma.
* The write-up says this itself ("at least as hard", §6 item 1). That remark is accurate, and it limits the status
  of Corollary 4.1 to a new sufficient condition.

## 2. Lemma 2.1 and Theorem 2.2 (single-form barrier)

**Lemma 2.1: correct.**
* With `alpha = prod_{v_x>0} z_x^{v_x}`, `beta = prod_{v_x<0} z_x^{-v_x}` and `s = sum v^+ = sum v^-`, the exact
  valuations give `v_pi(g) + v_{bar pi}(g) = s e_p - |<v,a_p>|` at split `p`, with full powers at inert primes and
  at `1+i`. So `N(g) = N^s e^{-X_v}`.
* The upper bound `|Z^v - 1| <= |<v',u>| <= |v'|_1 C e^{-W/4}` is right.
* `Lambda(v) = 2X_v/W` holds because `<v,beta_p> = 2<v,a_p>/e_p` on `Z^n_0`.
* No Gaussian gcd cost is discarded. This is the record's character inequality ([D] Prop 5.3(iii)) in Gaussian-gcd
  form.

**(2.2) and "fitting in A": correct.** By Ostrowski, `Newt(prod) = sum Newt`. Vertices carry nonzero
coefficients, so `h_{supp c}(beta) = sum_i |<v_i,beta>|`.

**Theorem 2.2: correct.** I re-derived each step.
* *Baire.* An irreducible complex hypersurface is not a countable union of proper closed subsets. So if all its
  points are degenerate, it lies in one `{Y^v = 1}`, hence equals one torsion coset `{Y^{v'} = zeta}`.
* *UFD.* `Y^{v'} - zeta` is irreducible for primitive `v'`. This gives (2.3).
* *Ostrowski and the cube.* These give `sum_j m_j |v'_{j,x}| <= 2` and
  `Lambda(c) = sum_j m_j Lambda(v'_j) >= ord_e(P) min_j Lambda(v'_j)`.
* *(CG_m)* then contradicts `ord_e >= m`.
* *Regime.* For full-rank actual clusters, Lemma 2.1 with `|v'|_1 <= 2d` gives `Lambda(v) > 1/2` once
  `W > 8 ln(2dC')`.

**Scope and novelty of Theorem 2.2.** The theorem is right, but it is the record's zonotope cap restated, not a new
barrier.
* *Order bound.* The step that bounds the order is round 11 compute Prop 2.2, round 13 fatpoint Thm 3.1 and round 16
  Prop 5.7. The only change is that the binomials may carry roots of unity `zeta_j != 1`. Such factors have order 0
  at `e` and only add cost, so the extension is immediate.
* *Classification.* Laurent polynomials with only degenerate zeros are monomials times products of torsion-coset
  binomials. This is a standard Baire/UFD fact.
* *Reading.* "Every single form has nondegenerate zeros; only the intersection can be toric" was already observed in
  the record. Round 11 compute §5.2 Remarks: "Each single `P_c` therefore has non-toric components through `e` with
  nondegenerate real points ... It is the **intersection** that is toric."
* Round 11 Thm 5.2 is also the record's first **system** certificate on an overdetermined, zonotope-free slab model.
  Gröbner saturation proves ZE on M8 at `m = 8 > zeta = 7` (§5 below).
* So the claim that "a proof must use the system" is consistent with the record. That the record lacked a system
  argument is not: it had a computational one.

## 3. Propositions 3.2, 3.4, 3.5

* **Prop 3.2 (cube): correct.** The leading-form ideals are monomial ideals, and the intersection argument is
  right. Exact check (`k0_small.py`): with both `{+-e_x}` and `L_pair`, `(SHIFT_e)` holds for `e < 2d` and fails at
  `e = 2d`, for `d = 2, 3, 4`.
* **Prop 3.4 (Sylvester): correct as a contrapositive of Thm 3.1.**
  * The `k` coordinate columns are width-cheap at `m' > 2jk` (t1rigour Prop 6.5(b)). Projection cost is at most
    width cost, so `K_F` lies in `Bs(V(A,m'))`.
  * A generic point of `K_F` is `L_pair`-nondegenerate.
  * The conclusion is only "fails for some `e < m'`". It does not locate the failure.
* **Prop 3.5 (`e <= j-1`): correct.**
  * The faithfulness lemma is a correct Gröbner-reduction plus Combinatorial-Nullstellensatz argument.
  * The subgrids `{-1,0,1}^{supp alpha} x {0}` (and their `+-e_x` shifts) satisfy
    `|<k,beta_p>| <= 2(|supp alpha| + 1) <= 2j <= w`.
  * As the write-up says, this range lies inside the zonotope range and proves nothing new about `Z`.
* **The `F_Q` remark** (degree `2j|Q| >= 2j log_2 n`) is right as far as it goes. I checked the analogous simple
  product obstructions `prod_{x in X} k_x^2 prod_{p in Q} s_p^{2j}`: every such form lying in all `J(A_l)`,
  `l in L_pair`, has degree `>= min over Q of (2j|Q| + 2(n - 2^{|Q|}))`, which is `>= 2j log_2 n`.
  * So no **product** obstruction exists below `8j`.
  * But the failures actually observed in the models are not of product type: syl8 at `e = 6` with `L = 2` gives
    `2L log_2 n = 12 > 6`. "Not obstructed by any known base point" therefore does not mean "not obstructed".

## 4. Corollary 4.1: quantifiers, rigidity input, uniformity

* **Uniformity in `R`: correct.** `(SHIFT)_A` depends only on the column directions and on `r`; no constant depends
  on `R`; `n_0(C)` is explicit.
* **Uniformity over profiles: correct.** The hypothesis is stated for every regime profile, but only the cluster's
  own box is used, so no template threshold is passed off as uniform.
* **Point hypothesis: correct.** With `L_pair`, only distinctness is used (`Z^{e_x} = z_x/z_n`,
  `Z^{e_x - e_y} = z_x/z_y`). A larger `L` would need `Z^l != 1` for every `l in L`.
  * This holds on full-rank clusters (round 26 §1.1 item 5).
  * It is not automatic on others.
  * It matters in practice: in syl8 a larger shift set raises the reach from 6 to 7 (§5.2(b)). A rich-`L` variant of
    the corollary would therefore need full rank or a reduction to full-rank subclusters, which the write-up does not
    supply.
* **Height condition (3.3): unproved, as the write-up says.**
  * It asks that small-height integer bases of the `~d^2` spaces `V(A_l, m'-1)` span `V(A, m'-1)` modulo
    `V_H(A, m')`.
  * With `H = U^2`, the record only has `dim V_H >= U/2` (t1rigour Thm 2.1(d)), against `dim V(A,m') ~ 3U/4`. So (3.3)
    is a genuine extra hypothesis.
  * With `H' = e^{m'W/40}` (t1rigour Remark 3.3) the codimension of `V_{H'}` in `V(A,m')` is `O(1/(w ln n))`, and (3.3)
    becomes much more plausible. It is still unproved.
* **Strength of the hypothesis.** By §7, `(SHIFT)_A` implies round 11's ZE-braid_H at `ceil(4w)`. So the corollary
  is: (new criterion) ⇒ (round-11 target) ⇒ (round-12/26 implication) ⇒ growth.

## 5. Independent recomputation of numerical claims

All code is mine: `kill_interp_checks/klib.py`, `klib2.py`, `klib3.py`. I read the record's profile **definitions**
(round 11 `model8.py`, `bigger.py`, `core.py`) and imported nothing from `interp_checks/`. "Exact" means fmpz/sympy
over `Q`. "mod p" means `p = 2147483629` and `p = 1000000007` (evidence). Heavy runs used
`prlimit --as=3500000000`, with at most two at a time.

### 5.1 Claims that reproduce

| claim (interp) | my result | script / output | label |
| --- | --- | --- | --- |
| Prop 3.2: cube, `(SHIFT_e)` for `e < 2d` | holds for `e < 2d`, fails at `2d` (`d = 2,3,4`; `{+-e_x}` and `L_pair`) | `k0_small.py` → `k0_small_out.txt` | exact |
| hexagon: `H = 1 2 3 1`, reach 3, `V(A,3)` alternating | yes. Note: with `{+-e_x}` only, the reach is 2 | `k0_small.py` | exact |
| Prop 3.3 1-D: gcd `(y-1)^m` (`m <= 4`); `m = 5` cofactor of degree 8 with `|prod roots| = 1008/462` | yes; `1008/462 = 24/11`; a cofactor root has modulus `2.6774` (nondegenerate) | `k1_onedim.py` → `k1_onedim_out.txt` | exact |
| 1-D: `max |S_l| = 3` (`l = +-10`); criterion 3; rational cyclotomic products reach 2 | yes (complete enumeration of products of degree `<= 13`) | `k1_onedim.py` | exact |
| 2-D: `S'` Sidon; `(SHIFT_e)` over `Q` for `e <= 2`, fails at 3; rational products reach 2 | yes. `H_S = 1 2 3 4 4`. **Only the shifts `+-(7,3)` are needed.** My enumeration runs over all centrally symmetric lattice polygons with vertices in `S` (112) and all cyclotomic assignments (2024 forms) | `k2_twodim.py` → `k2_twodim_out.txt` | exact |
| sdc7 `L = 2`: `H_A = 1 7 28 77 161 251 188 43 1`; `(SHIFT)` w.r.t. `L_pair` fails at 6 (188 vs 166) | yes | `k3_models.py` → `k3_sdc7_out.txt` | mod p (2 primes) |
| syl8 `L = 2`: `H_A = 1 7 28 77 161 259 287 175 42`; fails at 6 (287 vs 286) | yes | `k3_syl8_out.txt` | mod p (2 primes) |
| sdc7 at `e = 6` with all difference shifts: still fails | yes: with all 2002 shifts in `{-1,0,1}^7` having `|A_l| > 6`, 188 vs 166 | `k10_rich.py` → `k10_sdc7_e6_out.txt` | mod p |
| syl8 with the 98 shifts `{+-e_x, +-e_x+-e_y}`: reach 6 | yes. Adding all 280 weight-3 shifts still gives 286 vs 287 | `k12_98.py` → `k12_98_out.txt` | mod p |

### 5.2 Claims that do not reproduce

**(a) The zonotope order of sdc7 and syl8 is 7, not 6. In both models a single product certificate beats the
criterion.**
* The write-up says (§0 item 6; §4.1 table and bullets; §6 item 1) that in every computable model the criterion
  "reaches exactly the zonotope order (6 for sdc7 and syl8)", and that the first failure sits "where `zeta` also
  sits". The record had only `zeta in 6..7` for sdc7 and `zeta <= 7` for syl8 (round 12 `fake.md` §3.1 table).
* An exhaustive DFS over lattice zonotopes (`k6_zono.py`) finds **order-7 zonotopes in both**. I verified them
  exactly (`k7_zono_verify.py` → `k7_zono_verify_out.txt`):

```text
syl8 (L=2): k0 = (0,0,-1,-1,1,-1,0), generators e_6, e_5-e_6, e_5, e_3-e_4, e_3, e_2-e_4, e_2
sdc7 (L=2): k0 = (0,-1,0,0,0,0,0),   generators e_6, e_5-e_6, e_4-e_5, e_3-e_4, e_2-e_3, e_1-e_2, e_1
```

* For each: all 128 points lie in `A`; `c = prod (T_{v_i} - 1) delta_{k0}` lies in `V(A,7)` but not in `V(A,8)`;
  `P_c = Y^{k0} prod (Y^{v_i} - 1)` identically; `|c|_inf = 2` (syl8) and `1` (sdc7).
* **All generators lie in `L_pair`.** So this single form certifies `Y notin Bs(V(A,7))` for every `Y` with
  `Y^l != 1` (`l in L_pair`). That is the criterion's own hypothesis, at one order higher than the criterion reaches.
* The LP bound `zeta <= 7` is attained. Comparison per model:

| model | criterion (`L_pair`) | criterion (all `{-1,0,1}^7` shifts) | `zeta` | truth (max `m` with ZE^all) |
| --- | --- | --- | --- | --- |
| sdc7 `L = 2` | 6 | 6 | **7** | 7 (zonotope; nondegenerate base points at 8) |
| syl8 `L = 2` (M8) | 6 | **7** | **7** | **8** (round 11 compute Thm 5.2(b), Gröbner saturation) |

So in sdc7 the criterion is strictly weaker than one product form. In M8 it is weaker than a product form w.r.t.
`L_pair`, and below the order 8 at which the record **already certified ZE^all by a system argument** (round 11).
The M8 base locus at 8 includes points of two-term subtori, so the `L_pair` target ZE-braid may genuinely stop below
8 there. No whole two-term subtorus lies in `Bs(V(A,m))` for `m <= 8` (`k5_lightsubtori.py`; mod p, two primes).

**(b) "The failure is genuine" (§4.1) is false for syl8.** With all 2128 shifts `l in {-1,0,1}^7` having
`|A_l| > 6`, `(SHIFT_6)` **holds** in syl8 (287 = 287; `k9_allshifts.py` → `k9_syl8_cube_out.txt`). `(SHIFT_7)`
then fails (175 vs 141; `k10_syl8_e7_out.txt`). The `L_pair` failure at 6 is therefore an artefact of the shift
set; the rescue needs shifts of weight `>= 4`. The write-up tested the full shift set only on sdc7.

**(c) No slab box where the criterion beats zonotopes was found.** I tested random multi-cut sign profiles at
`n = 6, 7` (`T = 1`, `L in {1,2}`, `r` cuts, round-11 convention; `k8_random.py`, `k11_richrandom.py`). Reach was
computed mod `p` with `L_pair` and with all `{-1,0,1}^d` shifts. `zeta` came from an exhaustive zonotope DFS with
target "reach + 1".

| shift set | boxes | zonotope order > criterion reach | equal | criterion > zonotope order |
| --- | --- | --- | --- | --- |
| `L_pair` | 27 | **20** | 7 | **0** |
| all `{-1,0,1}^d` | 27 | **20** | 7 | **0** |

**(d) One slab box where the criterion beats zonotopes and reaches the truth: M8' (found by me; the write-up did
not compute it).**
* M8' = syl8 with `L = 3`: `|A| = 1877`, `rho = 11`. Round 11 compute Thm 5.2(a) proved by Gröbner saturation that
  every point of `Bs(V(A,11))` has `w_x = 1` or `w_x = w_y`. The record bounds `zeta <= 10`.
* `(SHIFT_e)` w.r.t. `L_pair` holds for **every** `e <= 10` (`H_A(e) = 1 7 28 77 161 266 357 386 329 195 63`, sum
  matched each time) and fails only at the top degree `e = 11` (7 vs 0). This is mod p, two primes (`k13_M8p.py`,
  `k14_M8p_confirm.py`).
* So Theorem 3.1(1) **reproves round 11 Thm 5.2(a) by linear algebra**: `Bs(V(A,11))` lies in the `L_pair`
  hypersurfaces, which is exactly the record's statement.
* An order-10 zonotope exists (`k0 = (-1,0,0,-1,0,-1,0)`, ten generators including `e_0 - e_1 - e_2`). I verified
  it exactly: all 1024 points lie in `A`, `c` is in `V(A,10)` but not `V(A,11)`, and `|c|_inf = 2`. So
  `zeta = 10 < 11` = criterion = truth.
* This is the one slab box I found where the multi-shift criterion is strictly non-zonotopal. It is overdetermined
  (`H_A(11) = 7 >= d`) and zonotope-free at its top order, though LP-subcritical. It supports the write-up's
  "non-zonotopal" claim better than the write-up's own 1-D and 2-D examples, which are not slab boxes.

**Net picture on slab boxes (w.r.t. `L_pair`).**
* *Criterion reaches the top degree `rho`:* cube; hexagon; three single-cut boxes; syl8 `L = 1` (`H = 1 7 28 21`,
  reach 3); M8'; 7 of 27 random boxes.
* *Criterion falls short of `rho` and of `zeta`:* sdc7 and syl8 `L = 2`; 20 of 27 random boxes.
* *Criterion above `zeta`:* M8' only.

The write-up's summary, "in every computable model it reaches exactly the zonotope order (6 for sdc7 and syl8)", is
wrong in both directions: below `zeta` in sdc7 and syl8, above `zeta` in M8'. The models do not show a trend that
would decide `(SHIFT)_A` at regime size.

## 6. Pre-registered tests T1-T5 (rerun by the killer)

**T1 (fake test): passes formally; no growth is claimed.**
* Corollary 4.1 has two inputs. Rigidity is exact vanishing of the certified system at `Z`, and it is violated by
  almost every fake because `Z_P` is Lebesgue-null (t1rigour §3 Remark (ii)). `(SHIFT)_A` is a property of the
  profile alone.
* Additional check: suppose a gap-admissible fake with pairwise distinct positions lies **exactly** on `Z_P` at
  regime size (as fake.md Cor 2.3 achieves in the underdetermined sdc7 model). Then Theorem 3.1(3) forces
  `(SHIFT)_A` to fail for that profile.
* So the corollary is consistent with the fake-test theorem because its hypothesis excludes exact fakes. The
  write-up says this correctly: `(SHIFT)_A` is at least as strong as the negative answer to the exact barrier
  question.

**T2 (`k = 4` negative control): not applicable, and I agree.**
* No Question-B statement is made. Theorem 3.1 is a true statement about finite sets, so it cannot "prove" the false
  `k = 4` analogue; the hypothesis would simply fail there.
* The analogous route-5 negative controls are the known base points, run under T5. The criterion correctly fails
  before each of them.

**T3 (split test): not applicable.** The criterion is Galois-blind. The write-up's flag is right: any transplant to
Question B would treat split and twisted forms alike, so it would also have to prove absolute Erdős #887.

**T4 (criticality): passes formally. On the models, the reach is mostly at or below the capped zonotope class,
with one clear exception.**
* Theorem 3.1 is an algebraic statement about base loci, so Lemma K1, the `tau_e = 10/3` cap, class L, the
  multi-point = one-point principle and fixed-`S` Subspace theorems do not enter.
* The one arithmetic input is the record's rigidity. It is LP-supercritical, and Gaussian gcd costs are kept.
* Zonotope (product) certificates are exactly the Liouville-tight single forms of Lemma 2.1, capped at
  `2w(1+o(1)) < m'` in the regime.
* On the slab boxes I computed (§5.2), the criterion's reach is
  * strictly below the zonotope order in sdc7, in syl8 `L = 2` (w.r.t. `L_pair`) and in 20 of 27 random boxes;
  * equal to it in the remaining random boxes, the cube, the hexagon and the single-cut boxes;
  * strictly above it in M8' (11 vs 10).
* So the mechanism is **not** a disguised single-form/Liouville tool: M8' is a genuine non-zonotopal success. But
  nothing in the data indicates that its reach on multi-cut boxes grows beyond the zonotope scale.
* The models are LP-subcritical with `rho/L <= 4` (M8': `rho/L = 11/3`), so the regime ratio `m'/j ~ 8` is out of
  reach. The question stays open, in both directions.

**T5 (calibration): passes.**
* sdc7, `L = 2`: nondegenerate base points at 8; the criterion stops at 6. Recomputed.
* Sylvester cheap subtori at `m' > 2jk`: Prop 3.4. Correct.
* 1-D example: base points at 5; the criterion stops at 3. Recomputed exactly.
* M8 = syl8: no nondegenerate base point up to 8 (record). The criterion stops at 6 (`L_pair`) or 7 (all cube
  shifts), so nothing true is excluded.
* The round-24 real family, the split Pell curves, the 2026 five-divisor construction, Huang's `Y_4` and `P^2` belong
  to Question B and are not touched. The 8- and 10-point families lie outside the regime.
* T5 cannot fail for this write-up, because Theorem 3.1 is correct. The calibration that matters is the opposite
  one: whether the criterion reaches true statements.
  * It under-reaches in sdc7 (6 vs truth 7) and in M8 (6, or 7 with all cube shifts, vs ZE^all true at 8).
  * Caveat for M8: round 11's base locus at 8 lies on pair subtori **or on 27 two-term subtori** such as
    `{w_0 = w_1 w_6}`. So ZE-braid, the `L_pair` target, may fail there before 8. I checked only that no whole
    two-term subtorus lies in `Bs` (`k5_lightsubtori.py`).
  * It reaches the truth in M8' (11 = truth, reproducing round 11 Thm 5.2(a)), in the single-cut boxes, the cube and
    the hexagon.

**Owner constraints: no violation found.**
* No fixed-`S` Subspace argument is used.
* Lemma 2.1 keeps `gcd(alpha, beta)` in full.
* Only degree-0 characters are used.
* `(SHIFT)_A` is assumed profile by profile, not as a uniform template threshold.

## 7. Novelty and literature

**Literature search.** I ran three WebSearch queries:
* de Boor–Ron least space, interpolation of exponentials, translates, leading forms;
* exponentials on finite lattice sets, difference operators, Hilbert function;
* least interpolation with Newton recursion over `S cap (S - l)`.

These returned the de Boor–Ron papers (["On multivariate polynomial interpolation"](https://ftp1.cs.wisc.edu/Approx/multiint.pdf),
[polintelim](https://ftp1.cs.wisc.edu/Approx/polintelim.pdf), [newideal](https://ftp1.cs.wisc.edu/Approx/newideal.pdf),
[two polynomial spaces](https://ftp1.cs.wisc.edu/Approx/twopolspaces.pdf)). Those papers define the least space
`Pi_Theta` as the span of least terms of exponential combinations, and `I^up` as the span of leading forms. The
searches also returned the box-spline and difference-operator literature
([de Boor–Dyn–Ron](https://apps.dtic.mil/sti/pdfs/ADA210559.pdf); Ben-Artzi–Ron,
[Trans. AMS 1988](https://www.ams.org/tran/1988-309-02/S0002-9947-1988-0961608-7/);
[the algebra of the box spline, arXiv math/0602019](https://arxiv.org/pdf/math/0602019)).
* **No source states the `S cap (S - l)` shift criterion.** The search tool's summaries say so explicitly. I read
  only titles and summaries, not full texts.
* Theorem 3.1 is therefore new to the record and not found in the literature. It is elementary: the identity
  `Delta_l e_Y = (Y^l - 1) e_Y` plus an induction on degree. It could well be folklore in least-interpolation or
  box-spline theory, so literature novelty is **unclear**.

**Component verdicts.**

| item | verdict |
| --- | --- |
| Lemma 3.0 / Theorem 3.1 (multi-shift zero estimate) | correct; **new to the record**; elementary; literature status unclear (not found) |
| Theorem 2.2 (single-form barrier) | correct; **repackaged known barrier**: the zonotope/product cap (round 11 compute Prop 2.2, round 13 fatpoint Thm 3.1, round 16 Prop 5.7) plus the standard classification of Laurent polynomials whose zeros are all degenerate. Its moral ("single forms have nondegenerate zeros; only the intersection is toric") is in round 11 compute §5.2 Remarks |
| Corollary 4.1 | correct implication; **repackaging of the known open lemma**. Its hypothesis implies a strengthened ZE^all_H (round 26) and is at least as hard as the exact barrier question (rounds 12–16), as the write-up concedes. The step from ZE to growth is the record's (rigidity, rounds 12 and 26) |
| "a genuinely non-zonotopal certification mechanism, the kind round 13/16 named as the missing tool" | **attribution wrong** (round 13 means the cap side); non-zonotopal on a slab box only through M8' (§5.2(d), not in the write-up); regime reach unknown |
| conclusion of Thm 3.1 with `L_pair` | **= round 11's ZE-braid target** (not cited); Prop 3.4 = the round-11 referee's ZE-braid refutation at `T = 1` |

**The conclusion of Theorem 3.1 with `L_pair` is round 11's "ZE-braid" target, which the write-up does not cite.**
* Round 11 compute §5.3 proposed **ZE-braid(A, m)** as "a purely algebraic form of ZE that implies it (and is
  stronger: it is a statement on the complex torus)". It says that every point of `Bs(A,m)` lies on
  `{x_x = 0} u {x_x = x_y}`, i.e. on `{Y^l = 1}` for some `l in L_pair`.
* The round 11 referee (`referee_compute.md` §2.8, **major**) proved ZE-braid **false** at Sylvester orders
  `n = 2^k >= 2^16` via the `log_2 n`-column cheap subtori. Round 11 README §§3, 7 lists it among the failed
  "algebraic reformulations that exclude rational subtori".
* t1rigour Prop 6.5 is the `T = 1` version: false for `m' > 2jk`, and not refuted at `m' = ceil(4w)`.
* Hence:
  * Theorem 3.1(1)/(3) with `L = L_pair` is a **sufficient linear-algebra criterion for ZE-braid(_H)**.
  * The write-up's Prop 3.4 ("the criterion is limited where it must be") **is** the round-11 referee's ZE-braid
    refutation, transported through Theorem 3.1.
  * Corollary 4.1's hypothesis implies ZE-braid_H`(A, ceil(4w))`. That statement is open, and it is at least as hard
    as round 12's open question on collision-free rational subspaces of `TC_H`.
* The target is therefore a mapped one. What is new is the criterion for reaching it.
* The criterion is not empty. On M8' it proves ZE-braid at the top order by rank computations alone. That is
  exactly round 11 Thm 5.2(a), which round 11 obtained by Gröbner saturation.

**Why the "missing tool" attribution is wrong as stated.**
* *Round 13 used the phrase on the cap side.* Round 13 README §7.1: "A non-zonotopal nondegenerate certificate on an
  overdetermined fake box would be the first sign that route 5 is **cappable**." There, "certificates" certify
  tangent-cone points (§3 table), the opposite direction.
* *Round 16 asked for something at regime size.* Round 16 §0 and §7.2 ask for "the first non-zonotopal lower-bound
  certificate at an irrational, collision-free direction of an overdetermined **regime** box".
* *The write-up's evidence is weak for that purpose.* Its examples are a 1-D set and a 2-D union of two translates,
  not slab boxes. Its slab-model data were misread (§5.2(a)). The best slab-box evidence for the mechanism is M8',
  found here; M8' is LP-subcritical.
* *Non-zonotopal exact certification on slab models was already in the record.* Round 11 Thm 5.2 used Gröbner
  saturation on M8 (`m = 8 > zeta = 7`) and M8' (`m = 11 > zeta = 10`). The multi-shift criterion matches round 11
  on M8' and falls short of it on M8.

**Overall novelty: genuinely new, modestly.**
* Theorem 3.1 is a new elementary criterion (new to the record; not found in the literature). It turns
  ZE-braid-type statements into rank equalities, and on one overdetermined, zonotope-free slab model it certifies
  beyond every zonotope.
* Everything around it is mapped territory:
  * the target (ZE-braid, round 11; ZE^all_H, round 26);
  * the obstruction at large orders (round 11 referee §2.8 / t1rigour Prop 6.5);
  * the single-form barrier (Theorem 2.2 = the zonotope cap).
* Corollary 4.1 is a reduction to a hypothesis at least as hard as the record's open lemma. It moves neither growth
  nor uniformity.


## 8. Corrected statement

> **Round 27 interp (killer-corrected).** No growth improvement and no uniform bound;
> `M <= (2/3+eps) log R/loglog R` stands.
>
> 1. *(proved; new to the record; elementary; not found in the literature)* **Multi-shift criterion (Thm 3.1).** For
>    finite `A` in `Z^d` and nonzero shifts `L`, suppose that for `e < m` one has
>    `J(A)_e = cap_{l in L} J(A cap (A-l))_e`. Then `Bs(V(A,m))` lies in `Union_{l in L} {Y^l = 1}` and
>    `TC(A,m)` lies in `Union l^perp`. Height version: (SHIFT_e) over `Q` for `e <= m-2`, plus (3.3).
>    * With `L = L_pair`, the conclusion is round 11's **ZE-braid(A,m)** (resp. ZE-braid_H).
>    * Prop 3.4 is the round-11 referee's refutation of ZE-braid at Sylvester orders `m' > 2j log_2 n`, transported.
> 2. *(proved; implication only)* **Cor 4.1.** If the cluster's `T = 1` box satisfies `(SHIFT_e)` w.r.t. `L_pair` for
>    `e <= ceil(4w) - 2`, together with (3.3) at `H = U^2`, then `n - 1 < 316 ln 4r`, hence `M = O_C(log log R)`. The
>    hypothesis implies ZE-braid_H`(A, ceil(4w))`. That is stronger than the record's open lemma ZE^all_H and at least
>    as hard as the exact barrier question. It is **open**: proved only for `e <= j - 1` (Prop 3.5), and false at
>    Sylvester orders `> 2j log_2 n` (Prop 3.4).
> 3. *(proved; a repackaging of the record's zonotope cap)* **Thm 2.2.** Every nonzero form of `V(A, m')` has a
>    Zariski-dense set of nondegenerate zeros when `Lambda(v) > Lambda(A)/m'` for all primitive `v in {-2..2}^d`.
>    This holds for full-rank actual clusters with `W > 8 ln(2dC')`. So no single form proves ZE^all.
> 4. *(evidence, mod p unless marked exact)* **Reach in models.** The criterion with `L_pair` reaches:
>    * sdc7 `L = 2`: 6 (zonotope order 7, exact; truth 7);
>    * syl8 `L = 2`: 6 (7 with all cube shifts; zonotope order 7, exact; ZE^all true at 8);
>    * M8' = syl8 `L = 3`: **11 = `rho` = truth** (zonotope order 10, exact). This reproves round 11 Thm 5.2(a) by
>      rank computations;
>    * 27 random multi-cut boxes: never above the zonotope order, strictly below it in 20.
>
>    Exact: the 1-D example (products 2 < criterion 3 < truth 4) and the 2-D two-translate example (products 2 <
>    criterion 3, using only the translate). The models are LP-subcritical and cannot decide `(SHIFT)_A` at regime
>    size.
>
> **Status:** conditional (new sufficient criterion for a mapped target). Route 5's open lemma is unchanged.

## Appendix: scripts

All under `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/kill_interp_checks/`.
Each script writes `<name>_out.txt` or a `*_log.txt`. Notes are in `notes.md`.

| script | content | label |
| --- | --- | --- |
| `klib.py`, `klib2.py`, `klib3.py` | own library: least spaces `M_e(S)`, Hilbert functions, `(SHIFT_e)` (exact over `Q`; mod p; standard monomials on the ternary cube) | — |
| `k0_small.py` | cube `d = 2,3,4`; hexagon; `V(hex,3)` | exact |
| `k1_onedim.py` | 1-D example: gcds, cofactor, `max |S_l|`, cyclotomic products | exact |
| `k2_twodim.py` | 2-D example: Sidon check, `(SHIFT_e)` (all shifts and `+-t`), complete enumeration of rational binomial-cyclotomic products | exact |
| `k3_models.py` | sdc7, syl8 `L = 2`: `H_A`, `(SHIFT_e)` w.r.t. `L_pair` | mod p, 2 primes |
| `k4_singlecut.py` | single-cut boxes: reach = `rho` | mod p |
| `k5_lightsubtori.py` | M8: no two-term subtorus `{Y^v = 1}` lies in `Bs(V(A,m))`, `m <= 8` | mod p, 2 primes |
| `k6_zono.py`, `k6b_zono.py` | exhaustive lattice-zonotope DFS (sdc7, syl8, M8', single cut) | exact search |
| `k7_zono_verify.py` | exact verification of the order-7 zonotopes in sdc7, syl8 | exact |
| `k8_random.py`, `k11_richrandom.py` | 27 random slab boxes: reach (`L_pair` / all cube shifts) vs `zeta` | mod p / exact search |
| `k9_allshifts.py`, `k10_rich.py`, `k12_98.py` | syl8 and sdc7 with 98, 378 and all cube shifts at `e = 6, 7` | mod p |
| `k13_M8p.py`, `k14_M8p_confirm.py` | M8': reach 11 (two primes), exact order-10 zonotope; syl8 `L = 1` | mod p / exact |

Literature (WebSearch summaries only, not full texts): de Boor–Ron
[multiint](https://ftp1.cs.wisc.edu/Approx/multiint.pdf), [polintelim](https://ftp1.cs.wisc.edu/Approx/polintelim.pdf),
[newideal](https://ftp1.cs.wisc.edu/Approx/newideal.pdf), [twopolspaces](https://ftp1.cs.wisc.edu/Approx/twopolspaces.pdf),
[limorg](https://ftp1.cs.wisc.edu/Approx/limorg.pdf); de Boor–Dyn–Ron [ADA210559](https://apps.dtic.mil/sti/pdfs/ADA210559.pdf);
Ben-Artzi–Ron [Trans. AMS 309 (1988)](https://www.ams.org/tran/1988-309-02/S0002-9947-1988-0961608-7/);
[arXiv math/0602019](https://arxiv.org/pdf/math/0602019); [arXiv math/0611940](https://arxiv.org/pdf/math/0611940).
None states the `S cap (S - l)` criterion.
