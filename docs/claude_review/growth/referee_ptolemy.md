# Referee report on `growth/ptolemy.md` (Ptolemy matching relations, Theorem D)

## Verdict

**Theorem D survives as a correct conditional theorem.** Assume Vojta's Main Conjecture with
`D = 0` for `Mbar_{0,8}` over `Q`. Then the uniform `C sqrt(R)` bound holds. I checked every
step and every constant. I found no hidden uniformity assumption, no discarded gcd or denominator
cost, and no drifting phase.

The unconditional identities are all correct. I re-verified each of them with independent code,
which shares nothing with the author's `gi.py`:

* Theorem A, the general Gaussian Ptolemy relation with mixed units;
* Theorem B (i) and (iv), with the tree-edge/block dictionary `(3.2)`;
* the algebraic core of Lemma 6.2;
* the `k = 5` family of Theorem C1;
* the anchoring Lemma 5.2;
* the canonical class and `K.beta`;
* every constant in Steps 3-5.

The write-up needs correcting in three places. None of them affects Theorem D.

1. **(Major, scope of the no-go.)** Section 7 says "`(L_k)` is false for `k = 5, 6` (Theorem C)"
   and that `(L_7)` "would be false if a balanced rational curve over `Q` exists". Both claims
   are wrong as justified. `(L_k)` lets a fixed proper Zariski-closed set `Z` be discarded.
   Theorem C1 lives on one rational curve, and item 109 on one elliptic curve, and either can be
   placed inside `Z`.
   * For `k = 5` the claim can be repaired. Infinitely many distinct balanced bounded-residue
     curves exist (Section 3.3, verified below), so `(L_5)` is indeed false.
   * For `k = 6`, `(L_6)` remains **unproved**.

   Accordingly, the slogan "no Ptolemy-only argument gives `H -> infinity` for `k <= 6`" holds
   only for arguments that apply to *every* configuration. That is the content of Corollary C3,
   which is correct as stated. An argument shaped like Theorem D (counting lemma plus an
   exceptional set) is not excluded at `k = 6` by anything in the note. The "exact threshold 8"
   is likewise established only in the `Z`-free sense.
2. **(Minor, novelty.)** The headline conclusion "Vojta's Main Conjecture implies the uniform
   theorem" is not new to the project. `growth/conditional.md` (Theorem 5.1, confirmed in
   `referee_conditional.md`) derives a *stronger* conclusion from another `D = 0` instance: an
   absolute `M_0`, and arcs `C R^alpha` for every `alpha < 1`. Theorem D is a genuinely different
   instance, on a different variety and using only cross-ratio data. Its value is structural: the
   form-blind Ptolemy content already suffices under Vojta. The phrase "logically independent
   instances" is informal; neither implication is known, but independence is not proved.
3. **(Minor, novelty of B and C.)** Several pieces were already in the notes:
   * Ptolemy is equivalent to Plücker on half-angle slopes (`ptolemy_squareclass_plucker_audit.md`).
   * The `k = 5` boundary-height dictionary is item 88 (`five_row_del_pezzo_arithmetic.md`).
   * A balanced, pairwise-coprime realization on a fixed elliptic curve in `Mbar_{0,5}` is item 89
     (`five_row_elliptic_height_compatibility.md`).

   The new parts are the following:
   * the general-`k` nested-block Lipschitz dictionary `(3.2)`;
   * the mixed-unit Gaussian form of Ptolemy;
   * a genus-0 `k = 5` family with `O(1)` balance and `|Y| = 1`;
   * the anchoring lemma;
   * `K.beta = 2^(k-3)(k-8)+k+2`.

   These are modest but correct.

All checks are in `growth/referee_ptolemy_checks/` and use exact integer or `Fraction`
arithmetic. Floats are used only for logarithms in balance and height comparisons.

---------------------------------------------------------------------------------------------

## 1. Theorem D, step by step

**Hypothesis `V_8`.** This is Vojta's Main Conjecture (Vojta LNM 1239, Conj. 3.4.3) for
`X = Mbar_{0,8}` with `D = 0` and `A = K + B`. The divisor `A` is ample (Keel–McKernan: the
pair `(Mbar_{0,n}, B)` is its own log canonical model).

* With `D = 0`, no set of places enters. `Z_eps` and `C_eps` depend only on `(X, A, eps)`.
  This is exactly the uniformity that failed in the earlier Subspace/Ru–Vojta attempts, and it
  is legitimately available here.
* `X` is smooth and projective over `Z` (Knudsen), and `Pic = N^1` is spanned by boundary
  divisors. Hence `h_K = sum_S c_S h_{D_S} + O(1)` for any fixed Weil functions.
* The hypothesis is a standard named conjecture. It is open: no nontrivial case of Vojta with
  `D = 0` is known on `Mbar_{0,n}` for `n >= 8`.
* I also checked that `V_8` is not trivially false at "naive" dense families. For example, take
  four points `p`-adically clustered mod a large `n` against four generic ones. The product
  formula forces pair contacts (`s = 2`, coefficient `-2/7`) of total size about `12 log n`, so
  `h_K < 0`.
* Oganesyan's disproof of the `R^alpha` (`alpha > 1/2`) conjecture would otherwise have
  contradicted `conditional.md`. It was withdrawn (arXiv 2107.09991, 22 Aug 2021), as
  `endpoint_recent_literature_check.md` records. It has no bearing on `V_8`.

**Step 0-1 (item 532 plus Markov).** I re-read `exact_nested_profile_residual_extraction.md`.

* The constants agree with item 532 at `k = 8`:
  - `E A <= 2(2^k - 1)/(M - k + 1) = 510/(M - 7)`;
  - `E B <= k(k-1)W/(4(M-1)) = 14W/(M-1)`;
  - `lambda = log(2/C_0)`. This is item 532's normalization; the handoff's `log(sqrt2/C_0)` is
    the older one, and the difference is harmless for `C_0 = 1`.
  - `W >= 4(M-1) lambda`, from (1) and `delta >= 2 lambda`.
* Markov gives `P(A > 4a) <= 1/4` and `P(B > 4b_0) <= 1/4`, since `A, B >= 0`. So at least half
  of the ordered 8-tuples of distinct rows are good.
* Fourier inversion with Cauchy–Schwarz gives
  `eta = sqrt((2^k - 1) * 4a) = 2 sqrt2 * 255/sqrt(M - 7)`.
* Item 532 (6) gives `|log n_T - w| <= eta w` for every nonempty `T`.
* Item 532 (9) with `C_0 < 2` gives `H_Sigma <= B/2 <= 28W/(M-1) = 3584 w/(M-1)`, using
  `W = 128 w`.

All correct.

**Lemma 6.1 (exceptional set counting).** The lemma is correct. The proof is terse in one place.
"`Z_1` proper" needs the set `B = {b : pi^{-1}(b) subset Z°}` to be constructible. It is: its
complement is `pi(M_{0,k} \ Z°)`, which is constructible by Chevalley. If `B` were dense, it
would contain a dense open set whose preimage lies in `Z°`, which is impossible. Over
`b not in B` the fibre `Z° cap pi^{-1}(b)` is a proper closed subset of an irreducible curve, so
it is finite. Its size is uniformly bounded, because the map is quasi-finite and of finite type.

Given the first `k-1` distinct points, the configuration determines the `k`-th point. So the
count is `(c(Z_1) + d(Z)) M^(k-1)`. The configuration of a tuple of circle points is intrinsic:
the Cayley map circle(`Q`) to `P^1(Q)` changes under re-anchoring by a `Q`-rational Möbius map.
So Lemma 6.1 applies with `Sigma` the image of the class.

**Lemma 6.2 (heights of a balanced Boolean point).** This is the heart of the proof. I checked
each ingredient.

* **Finite places.** The model Weil function at `p` is `(Phi.D_T)_p log p`, the thickness of
  the node of type `T` in the stable reduction, which equals the tree edge `l_T(p)`. I verified
  `l_T(p)` vs `v_p(n_T)` with a *different algorithm* from the author's. The author uses the
  quartet min-formula on Gromov products. I used `p`-adic cluster discs with the anchor sent to
  `infinity`:
  `l_T = max(0, min_{i,j in T} v(x_i - x_j) - max_{i in T, j notin T, j != 0} v(x_i - x_j))`.
  - The check covered 470 actual common-unit tuples (`k = 5..8`, half of them consecutive points
    in angular order), with item-532 blocks from threshold layers (nested chains sharing primes)
    and random anchors, at 6747 primes.
  - Every case gave `l_T(p) = v_p(n_T)` when `p` divides no residual, and
    `|l_T(p) - v_p(n_T)| <= 2 max v_p(t)` in general. The maximum observed deviation was 8.
  - Every case gave `|sum_p l_T(p) log p - log n_T| <= 2 H_Sigma`.
  - The checks also confirmed item 532 (7), `N gcd(P_i : i in S) = prod_{T >= S} n_T` for
    `|S| <= 3`, and Theorem B(i).
* **Archimedean place.** `lambda_{D_T,inf}` is bounded below because `D_T` is effective. Since
  `chi_Q^*[0] - D_T` is effective, `lambda_{D_T,inf} <= log^+ (1/|chi|) + O(1)`. I verified
  exactly, on 405,570 (tuple, `T`, crossing quartet) instances, the factorization
  `chi = prod_T n_T^{eps_T} * t_ab t_cd/(t_ac t_bd)`. It has exactly `2^(m-3)` blocks with
  `eps = +1` and `2^(m-3)` with `eps = -1`, all of them boundary blocks. Invisible blocks always
  cancel, anchored quartets included. Hence
  `lambda_{D_T,inf} <= 2^(m-2) eta w + H_Sigma + O(1)`.
* The bound is therefore
  `log n_T - 2H_Sigma - c_0 <= h_{D_T} <= log n_T + 3H_Sigma + 2^(m-2) eta w + c_0`,
  with `c_0` depending only on the fixed Weil functions. It is correct.

**Steps 3-5.** All constants check in exact arithmetic (`ref_check_anchor_K.py`, part D):

* the coefficients `s(8-s)/7 - 2 = -2/7, 1/7, 2/7` on 28, 56 and 35 divisors;
* positive mass `18`, negative mass `8`;
* `h_K >= 10w - 282 eta w - 60 H_Sigma - c`;
* `A`-mass `129`;
* `eps = 1/43` gives `3w + 99 eta w + 9 H_Sigma`;
* `(7 - 381 eta) w <= 69 H_Sigma + C'`;
* `M - 1 > 247296/5`;
* `eta <= 1/381` once `M >= 7.551e10`.

The final bound on `M` is independent of `R` and of the arc position. The transfer back to the
original arc costs a factor `4 ceil(C/C_0)`.

**Hidden-assumption audit (the user's pitfalls).**

* *Template-dependent thresholds:* none. `M_Z`, `C'` and `lambda` are fixed.
* *Discarded gcd or denominator costs:* none. The finite part is exact, and the residual
  valuations enter through the Lipschitz bound.
* *Drifting phase set to zero:* none. The archimedean term is bounded through the cross-ratio
  using only balance and the residues.
* *Non-uniform exceptional set:* none. `Z_{1/43}` is fixed, and Lemma 6.1 avoids it at a cost
  `O(M^7)`, uniformly in `Sigma`.
* *Extraction assumed:* no. It is item 532 plus Markov.

Remarks 2 and 3 are also correct.

* For `k <= 7`, `K = -(1/3) B_2` at `k = 7` and `K = -(2/5) B_2 - (1/5) B_3` at `k = 6`. So
  `-K` is effective and supported on the boundary, and `V_k` with `D = 0` is trivial on
  `M_{0,k}`.
* `k = 8` is the first `k` with `K.beta > 0`.

## 2. Unconditional parts

* **Theorem A.** `ref_check_A.py` uses actual circles with up to four split primes and exponents
  up to 3, all four units, true angular order, and every consecutive 4-window. It checks (1)-(5)
  exactly on 32,300 quadruples, 32,017 of them with mixed units:
  - `z_i - z_j = g_ij v_ij`;
  - `g g = Gamma X_k` in `Z[i]`;
  - the Gaussian identity `X2 v13 v24 = X1 v12 v34 + X3 v14 v23`;
  - the three products lie on one positive primitive ray;
  - the `|v_ij|` formulas for each unit ratio.

  Correct. The common-unit case is `ordered_residue_growth.md` (5), with the same `X_k` built
  from threshold layers. The mixed-unit Gaussian version is new but routine.
* **Theorem B.**
  - (i) and (iv) were verified as above. (ii) is standard Grassmannian algebra; the author's
    rank check passes. (iii) is trivially correct.
  - Novelty: the identification "Ptolemy = Plücker on rational half-angle slopes, and nothing
    more without size or divisor input" is `ptolemy_squareclass_plucker_audit.md`.
  - The notes already work on `Mbar_{0,n}` with boundary pullbacks (items 88, 104-112). Item 88
    already gives the `k = 5` boundary-height dictionary with explicit finite losses.
  - New: the general-`k` statement for the nested (shared-prime) item-532 blocks with the explicit
    `2 max v_p(t)` Lipschitz bound, and the `2 H_Sigma` aggregate.
* **Theorem C1.** `ref_check_C1.py` is an independent check. It confirms:
  - the factorization and all five relations as polynomial identities;
  - 11 distinct roots;
  - cyclic order `infinity, X_1 < X_2 < X_3 < X_4`;
  - balance within `log 20`;
  - `|Y| = 1` and `max |t'| = 144`, where the top block's `{2,3}`-part is `3`.

  It also *proves* pairwise coprimality of the reduced visible blocks for every `s` in the class.
  The primes dividing pairwise resultants are `{2,3,5,7,11,13,19,41}`. On the class, no visible
  block is divisible by `5, 7, 11, 13, 19` or `41`. The `{2,3}`-valuations are constant, because
  `v_2 < 4` and `v_3 < 4` on the class.
* **Lemma 5.2 (anchoring).** `ref_check_anchor_K.py` reimplements it with cluster trees on 300
  configurations and 8170 primes. It confirms `det(P_0, P_j) = +-1`, unchanged cross-ratios, and
  `v_p(det(P_i, P_j)) = delta_p + sum_{T >= ij} l_T(p)` with one `delta_p >= 0`. Correct.
  For C2, I re-derived the cancellation `det(P_i, P_j) = prod_{T not ni inf, T >= ij} d_T` times
  a `Sigma`-bounded factor. This conditionally confirms C2. It rests on item 109, a notes claim
  with its own audit that I did not redo.
* **Canonical class.** Independently, for `n = 5..12`, the boundary formula
  `sum (s(n-s)/(n-1) - 2) D_S` equals Kapranov's
  `K = -(n-2)H + sum (n-3-|I|) E_I` after the substitutions
  `D_{I u {n}} = E_I` (`|I| <= n-4`) and `D_{I u {n}} = H - sum_{J < I} E_J` (`|I| = n-3`).
  The same computation confirms `K.beta = 2^(n-3)(n-8)+n+2` and
  `(K+B).beta = 2^(n-3)(n-4)+1`. I also verified the identity analytically. The relation
  "heuristic exponent = `1 - K.beta`" is an identity about the author's heuristic count. It is
  heuristic and is not used in any proof.

## 3. The overclaim about `(L_k)`, and its repair

### 3.1 The problem
`(L_k)` reads: there exist `c > 0`, `eta_0 > 0` and a proper Zariski-closed `Z` such that every
balanced Boolean configuration off `Z` has `H_Sigma >= c w - O(1)`.

* The C1 family is the image of `s -> Phi_s`, a single rational curve in `M_{0,5}`.
* The item-109 points lie on a single elliptic curve in `Mbar_{0,6}`.

Taking `Z` to be that curve defeats each counterexample. So Theorem C does **not** show that
`(L_5)` or `(L_6)` is false. The claim "`(L_7)` would be false if a balanced rational curve over
`Q` exists" fails for the same reason. Corollary C3 is unaffected, because it has no `Z`-clause.

### 3.2 Why it matters
Theorem D itself is of the "discard a fixed `Z`, then count" type. The note uses Theorem C to say
that the Ptolemy route has "exactly one consistent target, the inequality at `k >= 8`". That
statement needs `(L_6)` and `(L_7)` to fail. Neither is established.

### 3.3 Repair for `k = 5`
Run the C1 construction with general `(a, tau)` whose 11 roots are distinct.

* The polynomial identities are formal.
* A CRT class fixes the valuation of every block at every prime dividing a pairwise resultant or
  a block content, so the residues are bounded on the class.
* Every visible block has degree 1 on the curve, so the map is birational and the labelled
  boundary parameters `tau_1..tau_4` give a curve invariant
  `cr(tau_1, tau_2, tau_3, tau_4)`.
* Varying `tau_4` therefore gives infinitely many distinct curves. Their union is Zariski dense
  in the surface `M_{0,5}`. Any proper `Z` misses some curve, and along that curve `H = O(1)`
  while `w -> infinity`. Hence `(L_5)` is false.

`ref_check_L5_dense.py` verifies this exactly for `tau_4 = 2, 11, 13`:

| `tau_4` | invariant | `max |t'|` | `|Y|` |
|---|---|---|---|
| 2 | `15/14` | 144 | 1 |
| 11 | `24/35` | 3000 | 1 |
| 13 | `36/49` | 140 | 1 |

### 3.4 `k = 6` and `k = 7`
Refuting `(L_6)` needs a Zariski-dense supply of such configurations in `M_{0,6}`, for example
infinitely many item-109-type curves. Items 104-112 show that such curves are delicate. This is
open.

## 4. Novelty against the literature

A bounded web search found no derivation of the `C sqrt(R)` bound from Vojta's conjecture on
`Mbar_{0,n}`. Literature on Vojta for `(Mbar_{0,n}, boundary)` concerns integral points. The
earlier referee found that Cilleruelo–Granville use Bombieri–Lang and `abc` only for sumsets and
squares in arithmetic progressions. The method is standard in spirit: Vojta's conjecture on
blow-ups yielding GCD-type bounds (Silverman 2005; McKinnon; McKinnon–Roth), here on Kapranov's
iterated blow-up of `P^5`, and the note should cite these.

## 5. Corrected statement

**Theorem D (as stated, correct).** Assume Vojta's Main Conjecture for rational points with
`D = 0`, `A = K + B`, `eps = 1/43` on `Mbar_{0,8}/Q`. Then for every `C > 0` there is an
ineffective `M(C) <= 4 ceil(C) M_*` such that every arc of length `<= C sqrt(R)` on
`x^2 + y^2 = R^2` (`R^2` in `Z`, any position) contains at most `M(C)` lattice points. Here `M_*`
is absolute.

**Unconditional (correct).**
* Theorem A in all parts.
* Theorem B (i)-(iv), with `|sum_p l_T(p) log p - log n_T| <= 2 H_Sigma`.
* Theorem C1 and Lemma 5.2.
* Corollary C3 for `k = 5` and, given item 109, `k = 6`: there is no `f -> infinity` with
  `H >= f(w)` for *all* balanced coprime-block solutions.
* `K.beta = 2^(k-3)(k-8)+k+2`, and `V_k` with `D = 0` is vacuous on `M_{0,k}` for `k <= 7`.

**Corrections.**
* `(L_5)` is false, by the multi-curve argument of Section 3.3, not by C1 alone.
* `(L_6)` and `(L_7)` are open.
* A single balanced curve never refutes `(L_k)`.
* The no-go "no Ptolemy-only argument for `k <= 6`" applies only to arguments that do not discard
  an exceptional Zariski-closed set. At `k = 6` an exceptional-set argument is not excluded.
* Vojta implying the uniform theorem is already in `conditional.md`, with a stronger conclusion.
  Theorem D adds a second, form-blind instance on `Mbar_{0,8}`.

## 6. Files

All files are in `growth/referee_ptolemy_checks/`:

| File | Purpose |
|---|---|
| `ref_common.py` | Independent Gaussian helpers |
| `ref_check_A.py` | Theorem A on actual circles, mixed units |
| `ref_check_B.py` | Theorem B(i),(iv), `(3.2)`, item-532 (7), and the Lemma 6.2 cross-ratio core, by the cluster algorithm |
| `ref_check_C1.py` | `k = 5` family; coprimality proved for every `s` in the class |
| `ref_check_anchor_K.py` | Lemma 5.2; canonical class vs Kapranov; Step 3-5 constants |
| `ref_check_L5_dense.py` | Repair of `(L_5)`: several distinct balanced bounded-residue curves |
| `orig/` | Copies of the author's scripts; all re-run and passing |
