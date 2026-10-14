# Referee report on round13/fatpoint.md (adversarial)

Scripts (all mine, written from scratch; no author code imported): `referee_fatpoint_checks/`
(`r_y32.py`, `r_enum.py`, `r_member.py`, `r_rankQ.py`, `r_sdc7.py`, `r_loc1.py`, `r_budget.py`, outputs `*_out.txt`).
Working notes: `referee_fatpoint_notes.md`.

## 0. Verdict

**The headline survives:** no growth improvement and no uniform bound; `M <= (2/3 + eps) log R/loglog R` stands. Nothing in
the write-up claims more, and I found no false theorem among the proved items.

**The proofs check out.** I re-derived every proof in §2, §3, §4.1-4.2 and §6. They are correct with the quantifiers stated:
Lemma 2.1, Props 2.2-2.6, Cor 2.7, Thm 3.1, Prop 3.2, Thm 3.3, Prop 4.1, Cor 4.2 (as an upper bound), Lemma 4.3, Thm 6.1,
Cor 6.2, Lemmas 6.3 and 6.5, Thm 6.4 and Cor 6.6. Theorem 6.4 in particular is a correct conditional reduction:
Sylvester-profile `MIN^pair` follows from `BLOCK` for one fixed box.

**Independent exact checks reproduce every number I tested:**
* `|Y_3(2)| = 2605` and `rho = 8`;
* the 91,272 multisets, split 3072 / 88,200;
* 1536 members, and 0 members above axis count 2;
* sdc7: `|A| = 757`, `Delta H = [1,7,28,77,161,251,188,43,1]`, and `Loc_8(1)` certified;
* M8': `Loc_10(1)` and `Loc_11(1)` certified;
* M8: deficit 21 at `b_0 = 1`;
* the Theorem 3.3 thresholds and the `n = 2^12` budget numbers.

**What must change:**
1. **(major, method)** The §5.1 justification of "mod-`p` non-membership is rigorous" is false as stated. There is a
   two-point counterexample (§2.5). The correct criterion is a rank condition. I checked it exactly for the decisive
   `Y_3(2)` run (`rank_Z = rank_p = 2353`), so those negatives are in fact rigorous. The justification must be replaced,
   and the other mod-`p` negatives (the `t = 2` control rows) are not certified.
2. **(major, interpretation of the decisive test)** §5.3 and summary item 5 misdescribe the exact outcome:
   * Exactly: a heavy product is in `J(Y_3(2))^perp_7` **iff it is a lattice zonotope that fits in `Y_3(2)`** (1536/1536).
   * The other 1536 multisets of "zonotope cost `<= 2j`" are **not** members. So "products with zonotope cost `<= 2j` do"
     is wrong.
   * 6,492 of the 91,272 multisets (6,300 of the "0/88,200") already violate a **`b_0 = 1` width count**, in the eight
     order-4 directions `v_w = (chi_w(a))_a` (`tau = 6`). The author never imposed this count.
   * So the claim that the failures happen "although every grid and line condition holds" hides that part of the new
     cost at `3.5j` is itself a plain count.
3. **(major, scope)** `Y_3(2)` is **not** the image of `A` for any `n`: Lemma 5.1 needs `8 ln 2n < 2`. It is a toy set at
   `j = 2`, where `3.5j = rho - 1`. The claims "regime-faithful at any `j`", "decisive" and the §5.3 dichotomy consequence
   are therefore evidence about one toy set, not about regime blocks.
4. **(moderate)** Unrefereed or refuted sibling claims are repeated as facts:
   * "`b_0 >= 8` forced for Sylvester" (loctheory Cor 5.3);
   * "`Loc(1)` is false for Sylvester profiles";
   * "the README cost restriction is wrong";
   * "what remains are the non-binary families, i.e. those containing complete XOR blocks".

   referee_loctheory refutes the last one with block-free non-binary `AG(3,2)` families. The others are unproved.
5. **(moderate)** "No certificate" for `M8` (`m = 8`) and `M8'` (`m = 10`) is limited to the search pool (support `<= 4`,
   entries `±1`). Theorem 3.3 does not exclude these orders. The summary's "certificates exist exactly where the box
   contains a zonotope" is a pool statement.
6. **(moderate)** The write-up silently changes the README's `Loc(b_0)`:
   * from the span of pullbacks over families with `rho < ceil(4w)`;
   * to the sum of generated ideals over all families.

   Its version is **weaker**. Implications "`Loc => ...`" stay valid for the README version, but "`Loc` holds" checks
   (sdc7, pal7, M8') are checks of the weaker version only.

Minor items are listed in §3.

## 1. Setup re-derived

I used the write-up's notation (`T = 1`, `w = sqrt(8 d ln 4r)`, `j = floor(w/2)`, `m' = ceil(4w)`, regime
`d >= 316 ln 4r`). I checked it against README §1:
* `kappa_c(u) = |Phi_c(u)|/(|c|_1 d^{m'})` (README l. 97);
* Step 0 is `kappa_c(hat u) <= eps_0 = C d e^{-W/4}/(m'+1)`;
* nearstrata Lemma 3.1 gives `d_xy >= W/2 - 2 ln C'`. Every use in Thm 3.1, Thm 6.1 and Thm 6.4 matches.

The tangent-space dictionary `TC(A,m) = {u : l_u^m in J(A)_m} = Bs(J(A)_m^perp)` is standard apolarity, and Lemma 2.1 is
exactly that statement for a subideal. The direction is right: a larger inverse system has a smaller base locus, so
`Bs(V_F) ⊂ TC` holds unconditionally.

## 2. Section-by-section checks

### 2.1 §1.4: the hypothesis

**The definition was changed.** The task (README §7.1) defines `Loc(b_0)` as
`J(A)_i = sum_F (pi_F^* J(pi_F A))_i`: the span of **pullbacks** (polynomials in `<= b_0` linear forms), over families with
`rho(pi_F A) < ceil(4w)`. The write-up instead uses `(I_F)_i`, the **ideal generated** by the pullbacks (§1.3), with no
cost restriction. Since `pi_F^* J ⊂ I_F` and the family set only grows, README-`Loc` implies fatpoint-`Loc`. Consequences:
* every implication "`Loc(b_0) => X`" (Lemma 2.1, Cor 2.7, Cor 6.2) remains valid for the README version;
* every verification "`Loc_m(b_0)` holds" (§5.2(a)) verifies only the weaker version.

Example: in sdc7 at `m = 8` there are only 421 primitive directions with `tau <= 7` in `{-2..2}^7`, while
`dim J(A)_8 = 3002`. The pure powers `l_v^8` cannot span `J(A)_8`, so the README's literal `b_0 = 1` statement fails there
for those directions, while the ideal version holds.

The write-up should state the change explicitly. Its justification, "loctheory shows the restriction is wrong", is the
overclaim that referee_loctheory flagged (moderate 3). Prop 2.4 there refutes one sentence of README §7.3(3). It does not
show that generators are missing from the restricted sum.

The remark "nothing is lost by restricting to `alpha_F <= i`" (via loctheory Cor 2.3) is correct: for `i < alpha_F`,
`J(B_F)_{<= i} = J(Gamma_F)_{<= i}` is generated by grid powers.

**"loctheory Cor 5.3 forces `b_0 >= 8` for Sylvester profiles"**: not proved. referee_loctheory (major 2) shows that only
the upper bound is proved. It also shows that 4-circuits are generated by 3 non-column directions at `j = 2`.

### 2.2 §2: local base loci and product certificates. All correct.

* **Prop 2.2(b).** In the bigraded splitting `S = Q[y_F] ⊗ Q[y_rest]`, `I_F = J(B) ⊗ Q[y_rest]` is bihomogeneous. The
  `(0,m)` component of `l_u^m` is `(beta.y_rest)^m`, which is not in `I_F` because `J(B)_0 = 0`. Correct over `C`.
* **Prop 2.3 (joins).** Correct. `nu_B <= rho(B) < infinity` because `J(B)_e = S_e` for `e > rho(B)`.
* **Prop 2.4 (span principle).** Correct. `ell^m in V_F(m)` because `g(D_F) ell^m = 0` for every `g` of positive degree.
  The Galois argument for (b) is fine, since `Bs` is cut out by rational forms.
* **Prop 2.5.** The identity `<f, P> = m! Pol(f)(z_1..z_m)` holds on powers `l_w^m`, hence for all `f`. Correct.
* **Prop 2.6.**
  * (a) The `t`-homogeneous components have `t_q`-degree `<= N_{F,q}`.
  * `J(Gamma_F) = (y_q^{|I_q|})` holds for any product of finite sets.
  * `J(B)_e = J(Gamma)_e` for `e < alpha_F`: I re-derived this independently from the Davis-Geramita-Orecchia linkage
    formula for complete intersections, `Delta H_Gamma(e) - Delta H_B(e) = Delta H_{Gamma\B}(s - e)`. Correct.
  * (b) The top component has exact `t_q`-degree `N_{F,q}`, since `C[t]` is a domain. Correct.
* **Cor 2.7.** For families with `alpha_F > m'`, the condition `N_F < alpha_F` is automatic, and the grid part is the
  count at the member directions. Correct.

### 2.3 §3: the width system

**Theorem 3.1** is correct, for actual clusters and for gap-admissible fakes:
* every column direction `1_{P_p}` has `tau <= 2j < 8j <= m`, so the count applies;
* summing `log p` over **primes**, not distinct patterns, makes **duplicate columns** harmless;
* twins are distinct directions and are counted separately;
* pairs through the base point are covered by `z = e_x` in `k'`-coordinates;
* the contradiction needs `j >= 1` and `W/2 - 2 ln C' >= 0`, both implied.

**Proposition 3.2** is correct. It even holds with `sum g_i` instead of `2 sum g_i`.

**Theorem 3.3** is correct (Hadamard uncertainty plus Cauchy-Schwarz). My own boxes give:

| model | `tau_0` | `T_core` | certificates possible only for `m <=` |
| --- | --- | --- | --- |
| sdc7 | 4 | 28 | 7 |
| `Syl(8,2)` | 6 | 28 | 8 |
| `Syl(8,3)` | 8 | 42 | 10 |

All three agree with the author (`r_sdc7_out.txt`, `r_budget.py`).

Regime numbers reproduced: at `n = 2^12`, `m' = 2256`, `j = 281`, ratio `1.153`, bounds `9058.3 <= sum sigma <= 10446`.

Sharper than the write-up says: for general Hadamard-core profiles, the ratio RHS/`m'^2` is about
`sqrt(d/(512 ln 4r)) + 1/4`. That is about `1.03` at the regime boundary `d = 316 ln 4r` (e.g. `d = 2965`, `r = d`: `1.032`).
So the budget is non-binding but nearly tight at the boundary for every `r`, not only "near `n = 2^12`". The claim "not
binding in the regime" is correct.

**§3.1, sdc7.** I verified `Loc_8(1)` independently and rigorously:
* `rank_p` of the width-power generators in `S_8` is 3002;
* this equals `dim S_8 - Delta H_A(8) = 3003 - 1`, with `Delta H` exact (`fmpz`);
* since mod-`p` rank `<=` `Q`-rank `<=` `dim J(A)_8`, equality holds.

The Baire argument for nondegenerate points near `u*` is correct: a rational hyperplane containing a 6-dimensional
analytic piece of `{Phi = 0}` would divide `Phi` over `Q`. It does depend on R11 Prop 6.2 (`Q`-irreducibility of the
octic) and fake Prop 3.1 (regularity of `u*`), both quoted, not re-checked here.

**§3.4.** "`Loc(1)` is false for Sylvester profiles" is not proved. Lines enter at `3j`, but whether width powers of
signed line sums generate them depends on the parity of `j`:
* `tau(y_a +- y_b +- y_c) = 3j - 1` for odd `j`;
* at `j = 1` they do generate (referee_loctheory);
* at `j = 2` they do not (mod-`p` evidence only).

### 2.4 §4: parity families and complete blocks

* **Prop 4.1** is correct, given loctheory Thms 3.3 and 3.4, which referee_loctheory confirmed.
* **Cor 4.2** is correct as stated, for binary families.
* **The closing sentence of §4.1 is false:** "What remains are the non-binary families, i.e. those containing complete XOR
  blocks". Block-free non-binary column families exist: `AG(3,2)` and 889 others at `t = 4` (referee_loctheory major 1).
  They carry a mod-`2^{t-2}` congruence and are not generated by circuits plus width powers. Non-column families also
  remain (the author's §6.3 item 2 says so).
* **Lemma 4.3** is correct. The independence argument in (b) works because distinct `I` give distinct monomials in the
  algebraically independent `c_i`. The heavy/light class-sum formulas and the primitivity of `c(z)/|H_0|` are correct.

**Missed necessary condition.** Prop 2.6(b) applied to `F_U` gives `N_U <= rho(Y_t(j))` for every complete block. By
loctheory, `rho(Y_t(j))` is known only from below (`>= 4.16j`, `5.87j` for `t = 3, 4`). If the true values were `< m' ~ 8j`,
this would be a rigorous count constraint that every certificate must meet: each Fano block avoided by `>= m' - rho`
generators. It belongs in §4.2.

### 2.5 §5.1: certification direction (major)

**The claim.** The write-up says: "since that kernel contains the reduction of the rational kernel,
`T_m(mod p) ⊃ T_m(Q) mod p`, so a mod-`p` non-membership is a rigorous non-membership over `Q` (all products have
`p`-integral coefficients)". **This is false in general.**

**Why.** A rational `c` realising a `p`-integral moment vector can need `p` in its denominator.

**Counterexample.**
* `A = {0, p} ⊂ Z`, `m = 1`, `P(u) = u`.
* Over `Q`: `J(A) = (x^2)`, so `J(A)_1 = 0` and `P in J(A)_1^perp`.
* Mod `p`: `V(A,1) = span(1,-1)` and `mu(1,-1) = -p ≡ 0`, so `T_1(mod p) = 0`.
* So mod `p` declares `P` a non-member, wrongly.

**The correct criterion.** Membership is solvability of `[K_{<m}; E_m] c = [0; v]`. With
`K_{<= m} := [K_{<m}; E_m]` the evaluation matrix in degrees `<= m`: if `rank_p K_{<= m} = rank_Q K_{<= m}`, then
rational solvability implies mod-`p` solvability. Only then is a mod-`p` non-membership rigorous. The author's own library
docstring (`fplib.tangent_space_mod`) says "used for search, never as a certificate", which contradicts §5.1.

**Status of the reported negatives:**
* **`Y_3(2)`, `e = 7`:** `rank_Z K_{<=7} = 2353` (exact `fmpz`, 9 s), equal to the rank mod `2147475331` and mod
  `10^9+7` (`r_rank7_*_out.txt`). So the 0/88,200 and 0/4705 negatives are rigorous.
* **`Y_3(2)`, `e = 8`:** `rank = |Y|`, so trivially rigorous.
* **`t = 2` control rows** ("0/1748", "0/2000"): not certified, since no rank comparison is given. They are evidence.

My own membership computation (`r_member.py`, prime `10^9+7`, my own `T_7` and annihilator) reproduces the author's
counts exactly.

### 2.6 §5.2: exact models

I reproduced these independently (`r_sdc7.py`, `r_loc1.py`):
* sdc7: `|A| = 757`, and the `Delta H` vector above;
* M8: `|A| = 1037`, `Delta H = [1,7,28,77,161,259,287,175,42]`;
* M8': `|A| = 1877`, `Delta H = [...,63,7]`;
* `Loc_8(1)` for sdc7, and `Loc_10(1)`, `Loc_11(1)` for M8': certified;
* M8 `b_0 = 1`: `rank_p = 294` against target 315, a deficit of 21 (evidence, correctly labelled).

I did not recheck the `b_0 = 3` claims (pal7, M8). They agree with judge J6(c).

**Pool limits.** The certificate table is limited to the pool (balanced, support `<= 4`, entries `±1`) and to width
counts over `{-2..2}^d`.
* For `M8` at `m = 8` and `M8'` at `m = 10`, Theorem 3.3 leaves room. For `M8` at `m = 8`: `sum sigma <= 20`,
  `sum 1/sigma <= 3.5`, and `64 <= 70`.
* So "M8 has no certificate yet no nondegenerate TC point" and "certificates exist exactly where the box contains a
  zonotope" are pool statements.
* The positive half is rigorous: every certificate found is a fitting zonotope, checked exactly.

### 2.7 §5.3: the complete Fano block (major reinterpretation)

**Scope.** `Y_3(2)` is not `pi_U(A)` for any `n`. Lemma 5.1 needs `2^t ln 2n < j`, i.e. `8 ln 2n < 2` (referee_loctheory
minor 5). In the regime the block image is `Y_3(j)` with `j ~ sqrt(2 n ln n)`. Calling the test "regime-faithful at any
`j`" because "the relevant ratio is internal to the block" is an unproved assertion. At `j = 2` the "`3.5j`" degree is
`rho - 1`, so one cannot separate "enters at `3.5j`" from "enters just below `rho`".

**What the exhaustive data actually say** (`r_y32.py`, `r_enum.py`, `r_member.py`, all exact or rank-certified):
1. **Cheap directions.** Exactly 43 primitive directions in `{-2..2}^7` have `tau_Y <= 6`:
   * 7 axes (`tau = 4`);
   * 28 signed line sums (`tau = 6`);
   * 8 **order-4 directions** `v_w = (chi_w(a))_{a != 0}`, with `tau = 6` and values `{-12,...,12}` in steps of 4.

   This matches loctheory's "43 of 3025".
2. **The pool is simple.** The heavy generators with `|H'| = 2` are exactly the 42 forms `2(t_a ± t_b)` (pairs of axes,
   both signs). So 91,272 multisets with axis counts `<= 4` and line counts `<= 5` is reproduced:
   * 3072 with max axis count 2;
   * 88,200 with max axis count 3 or 4.
3. **The omitted width counts bite.** 6,492 multisets violate the `b_0 = 1` count `#{i : <v_w, g_i> != 0} <= 6` in some
   order-4 direction:
   * 192 with max axis count 2;
   * 6,048 with max axis count 3;
   * 252 with max axis count 4.

   The author's filter did not include these.
4. **Membership is exactly the fitting-zonotope condition.** Among all 91,272, a product is in `J^perp_7` iff the lattice
   zonotope `k_0 + sum_i [0,1] lambda_i g_i` (`lambda_i in {2,4}`) fits in `Y_3(2)` for some `k_0`. That is 1536 members,
   all of them fitting zonotopes; no non-zonotope member exists. Of the 3072 multisets with "zonotope cost `<= 2j`", 1536
   are **not** members.
5. **Not all failures are counts.** 83,244 width-count-feasible multisets are non-members. They are not explained by any
   `b_0 = 1` count, consistent with referee_loctheory's certified 6-direction-family generators in degree 7.

**Correct reading.** In this toy block at degree 7 = `rho - 1`, the only product certificates in the heavy pool are
genuine zonotopes. Part of the new obstruction at `3.5j` is a plain width count, the order-4 directions. Whether
`tau(v_w) = 3.5j + O(1)` explains why non-parity generators "start at `3.5j`" in general would need `Y_3(j)` for larger
`j`.

**Consequences for the summary:**
* "products with zonotope cost `<= 2j` do (1536/3072)": only half do.
* "heavy generators pay the zonotope cost": better stated as "only zonotopes survive".
* "although every grid and every line condition holds": true, but the omitted order-4 width conditions are `b_0 = 1`
  counts that already kill about 7% of the multisets.
* The qualitative message (no certificates beyond zonotopes) is, if anything, **strengthened**. The dichotomy (`W_t`) and
  its consequence remain heuristics resting on one toy set.

### 2.8 §6: conditional theorems. Correct, with minor slips.

* **Thm 6.1.** Correct. The constants check: `|<z_i, hat u>| >= C^{-1} e^{-(X_{z_i} - W/2)/2}` from `|u|_inf <= Theta`.
* **Cor 6.2.** Correct as an implication. The final step "`m' ln(Cd) = O(sqrt(n log n) log n)`" assumes
  `ln r = O(log n)`. In the regime `ln 4r` can be as large as `d/316`, which makes `m' ~ 0.64 d`; then `W >= 3 n ln n`
  alone does not dominate `m' ln d` for small `eta`. Fix: use `W >= sum_{p <= p_r} log p ~ r ln r`, so `m' ln d = o(W)`
  for every `r`. Minor.
* **Lemma 6.3.** Correct:
  * the convolution annihilates degree `< sum m_b` because the blocks use disjoint non-base coordinates;
  * `Phi_c` factorises with the multinomial;
  * `nu` is additive on products and monotone.
* **Lemma 6.5.** Correct. Lojasiewicz on the compact normalised set with `v_{x_0} = 0`, and
  `dist_inf(v, Br) >= sigma/2`.
* **Thm 6.4.** Correct:
  * room: `B j_0 + r <= j` for `j delta_0 >= 9 + m_0(8 + delta_0)`;
  * `P(x+y in U) = (N_0-1)/(n-1)`;
  * coset probability `<= B N_0/(n - N_0)`;
  * the `Q_1` bookkeeping with `N_2 = binom(n,2)`;
  * `H_0^B <= |A|^2`, since `ln|A| >= 2j ln(n/2j)`.

  The decomposition depends only on the pair data, as needed for `MIN^pair`. The note that `BLOCK` is false unless
  `m_0 <= 2 j_0 t_0` (binary-expansion span) is right, so `t_0 >= 5`. Sharper: costs `2 j_0 s + 2(N_0 - 2^s)` for mixed
  spans give further necessary conditions, but these do not change `t_0 >= 5`.
* **Cor 6.6.** Correct as an implication. "For `n >= 2900` the box is in the regime" is a slip: `d >= 316 ln 4d` needs
  `n >= 2966`. 2900 is t1rigour's exact threshold for `m_max >= 4w`, which is what is actually used. Harmless.
* **§6.2 remark on generic profiles.** Correct and conservative: for a full-trace block,
  `A_b = {balanced, |k|_1 <= 2 j_b}` has `rho <= 4 j_b`, not only `<= 8 j_b`.

## 3. Issues

| # | severity | where | issue | fix |
| --- | --- | --- | --- | --- |
| 1 | major | §5.1, appendix | The certification argument for mod-`p` negatives is false. Counterexample: `A = {0,p}`, `m = 1`. | Replace by `rank_p K_{<=e} = rank_Q K_{<=e}`. Verified for `Y_3(2)`, `e = 7` (2353) and `e = 8`. Relabel the `t = 2` control negatives as evidence. |
| 2 | major | §5.3, summary item 5, §0 row 6 | "Zonotope cost `<= 2j` suffices" is wrong: 1536 of 3072 fail. Exact: membership iff a lattice zonotope fits in `Y_3(2)`. 6,492 multisets fail a `b_0 = 1` order-4 width count that the filter omitted. | State the exact equivalence. Separate width-count failures (6,492) from the rest (83,244). |
| 3 | major | §5.1 item 2, §5.3, summary | `Y_3(2)` is not `pi_U(A)` for any `n` (Lemma 5.1 needs `8 ln 2n < 2`). "Regime-faithful at any `j`" is unproved, and at `j = 2`, `3.5j = rho - 1`. | Call it a toy lattice set. Keep the dichotomy (`W_t`) explicitly heuristic. |
| 4 | moderate | §4.1 end | "What remains are the non-binary families, i.e. those containing complete XOR blocks" is false: block-free non-binary `AG(t-1,2)` families exist. | Add affine (odd-coset) blocks and non-column families to what remains. |
| 5 | moderate | §1.4, §3.4 | "Cor 5.3 forces `b_0 >= 8`" and "`Loc(1)` is false for Sylvester" are unproved (referee_loctheory major 2; parity of `j`). | Drop or label as open. |
| 6 | moderate | §1.4, summary | The `Loc` definition is changed: ideal generated instead of span of pullbacks, and no `rho < ceil(4w)` restriction. "README restriction is wrong" is an overclaim. | State that fatpoint-`Loc` is implied by README-`Loc`, and that §5.2(a) verifies only fatpoint-`Loc`. |
| 7 | moderate | §5.2(b), summary | "No certificate" for `M8` (`m = 8`) and `M8'` (`m = 10`), and "certificates exist exactly where zonotopes exist", hold only within the pool. | Scope to the pool. |
| 8 | minor | Cor 6.2, Thm 6.4 | Asymptotic `O(sqrt(n log n) log n)` assumes `ln r = O(log n)`. | Use `W >= (1+o(1)) r ln r`. |
| 9 | minor | Cor 6.6 | "`n >= 2900` in the regime": the regime inequality needs `n >= 2966`. | Cite the `m_max >= 4w` threshold instead. |
| 10 | minor | Thm 3.3 reading | Tightness is about 1.03 at the regime boundary for all `r`, not only "near `2^12`". | Rephrase. |
| 11 | minor | §4.2 | The necessary count `N_U <= rho(Y_t(j))` (Prop 2.6(b)) is not used, and could be a rigorous constraint once upper bounds for `rho(Y_3(j))`, `rho(Y_4(j))` are known. | Add it as an open lever. |

**Not found:**
* No non-local generator, cheap-subtorus effect or twin/duplicate-column effect invalidates a proved statement. Duplicate
  columns are handled by summing over primes (Thm 3.1, Prop 3.2), and Lemma 6.3 allows any multiplicities.
* No model is claimed to be LP-supercritical. The write-up correctly says that all exact models are LP-subcritical and
  that `BLOCK` needs `N_0 >= 32`.
* The only "regime-like" claim is the `Y_3(2)` one (issue 3).

## 4. Corrected statement

At `T = 1`, `w = sqrt(8 d ln 4r)`, `j = floor(w/2)`, `m' = ceil(4w)`, regime `d >= 316 ln 4r`, unless a model is named.
**No growth improvement and no uniform bound; `M <= (2/3+eps) log R/loglog R` stands.**

(1) **Proved, any finite `A ∋ 0`, any `m`, any family set.** Write `I_F` for the ideal generated by `pi_F^* J(B_F)`.
* `Bs_R(∩_F (I_F)_m^perp) ⊂ TC(A,m)`.
* Equality holds under `Loc_m` in the generated-ideal form. That form is implied by the README's pullback-span form with
  the `rho < ceil(4w)` restriction, so implications from it hold for the README version.
* Single-family base loci lie in `span_Q F`, and joins with `sum nu < m` are base points.
* Span principle: a nondegenerate base point needs a `Q`-irreducible base component spanning `C^d`; an algebraic one has
  degree `>= d`.

(2) **Proved: product certificates.** `prod <z_i,u> in K_m^perp` iff `Pol(f)(z_1..z_m) = 0` for all `f in K_m`. The local
conditions are:
* the width count `#{i : <v,z_i> != 0} <= tau_A(v)` (exact for `b = 1`);
* `N_F < alpha_F` with grid counts suffices;
* `N_F <= rho(B_F)` and the grid counts are necessary.

(3) **Proved, regime.**
* Thm 3.1: for sign profiles with the pair inequalities and `W > 8 ln C'`, no product of `m >= 4w` collision forms
  satisfies the `b_0 = 1` counts. This covers actual clusters and gap-admissible fakes, and is insensitive to duplicate
  or twin columns.
* Prop 3.2: certificate generators are heavy on average.
* Thm 3.3: for Hadamard-core profiles, `m^2 <= (2d + min(tau_0, m)) T_core/n`. It is not binding in the regime, but
  within a factor of about 1.03 at `d = 316 ln 4r`.

(4) **Proved, given loctheory Thms 3.3-3.4.**
* Parity (binary) families impose only counts: `<= j|s| - 1` per circuit suffices, `<= j|s|` is necessary.
* Complete XOR blocks reduce to class-sum products on `Y_t(j)`.
* Not proved: that complete blocks are the only non-binary column families. Block-free `AG(t-1,2)` families exist.

(5) **Exact, toy set `Y_3(2)`.** This is `Lambda_U ∩ [-2,2]^7`, which is **not** the image of any regime or Sylvester box.
Take the 42 heavy forms `2(t_a ± t_b)` and the 91,272 multisets of 7 with axis counts `<= 4` and line counts `<= 5`.
* A product lies in `J(Y_3(2))^perp_7` **iff** it is a lattice zonotope fitting in `Y_3(2)`: 1536 members, all zonotopes.
  The negatives are rigorous via `rank_Z K_{<=7} = rank_p K_{<=7} = 2353`.
* 6,492 multisets already fail a `b_0 = 1` width count in an order-4 direction (`tau = 6`).
* The remaining 83,244 non-members are not explained by any width count.
* Mixed light/heavy samples agree (0/4705 above zonotope cost; negatives rigorous).
* Heuristic only: if fitting-zonotope behaviour persisted in regime blocks `Y_t(j)`, coset-difference certificates would be
  impossible for gap-admissible data.

(6) **Proved implications.**
* Thm 6.1 / Cor 6.2: `Loc(b_0)` + local certificate + height `psi <= (1-eta)W/8` + tight characters (Ar) imply `n < n_1`,
  hence `n = O(loglog R)`. The asymptotic step needs `W >= (1+o(1)) r ln r`.
* Thm 6.4 / Cor 6.6, no locality needed: `BLOCK`, i.e. `ZE-TC^braid` for one `Syl^box(2^{t_0}, j_0)` at `m_0 > 8 j_0`
  (needs `t_0 >= 5`), implies `MIN^pair` for Sylvester-type profiles at every large `n`. So Sylvester-profile clusters
  have `n < n_0`. This does not reach fake (Paley-type) or generic profiles.

(7) **Exact models, all LP-subcritical; certificate claims limited to the search pool.**
* No exact fake model reaches `rho/(2L) > 2`.
* Certified, generated-ideal form: `Loc_8(1)` for sdc7 (`L = 2`), and `Loc_10(1)`, `Loc_11(1)` for M8'.
* Reported, not rechecked: `Loc_7(3)` for pal7 and `Loc_8(3)` for M8.
* So the sdc7 width system's real base locus equals `TC` and contains certified nondegenerate, well-separated points.
  This is an underdetermined system, `H_A(8) = 1`.
* In the pool (support `<= 4`, entries `±1`), every certificate found is a fitting zonotope (24, 57, 696). None exists
  for `M8` at `m = 8` or `M8'` at `m = 10` within the pool; neither order is excluded by Thm 3.3.

**Method.** A mod-`p` non-membership is rigorous only when the mod-`p` and rational ranks of the degree-`<= e`
evaluation matrix agree.
