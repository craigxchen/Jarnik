# Referee report: `round4/flipped.md` (fully flipped Hadamard phase lemma)

**Verdict.** The result survives, with corrections. It does not prove the fully flipped phase
lemma for any non-Walsh core, and it does not claim to.

What I re-derived and checked independently, and found correct:

* the two-logarithm reduction (Lemma 1.1, Props 2.1-2.2);
* the conditional Theorem L;
* the integer prime-field uncertainty Lemma 6.2;
* Prop U, the sector lemma and Prop 6.4;
* Prop 7.1.

Four claims are overstated or need repair.

1. **Prop B is mis-read.** "Known Baker-type bounds miss by a factor of order `c log M`" is not
   what Prop B proves.
   * The `log M` factor is proved only when `W_F <= (1/2 - delta)(W+D)`, with constant `delta/2`.
     The constant `1/8` holds only for `W_F <= (W+D)/4`.
   * When `W_F >= (W+D)/2` the proof gives a constant factor only, and weight data show that the
     product of heights really can be `Theta(W+D)`.
   * The formal statement of Prop B, with its `min(., 0.87)`, is correct. The qualitative conclusion
     "no proved two-logarithm bound can work" also survives.
2. **Liouville "misses by exactly `F/2 - D/4`" should read "by at most".** Equality holds only for
   balanced labels. "Otherwise it gives nothing" (Cor 3.3) is false: (3.1) is a lower bound on
   *every* label weight.
3. **The "heavy-flip regime is open even conditionally" conclusion is misleading. This is the main
   new finding of this report.**
   * For the canonical Prop-P assignment, the **pair inequalities alone** (from row infinity) give
     `W_F <= W/2 + (M-1)(4 log C - D)/2`.
   * So the heavy regime `W_F >= (1/2 + delta)W` contains only orders `M <= (e/4) C^(2/delta)`.
     Prop 7.1 is therefore vacuous for large `M` on that class.
   * Exact LP certificates show the same cap `W_F/W <= 1/2` (asymptotically) for every Paley-12
     assignment tried, including adversarial ones with capacity 2 to 5, and for random Paley-20
     assignments.
   * Together with Theorem L, what remains of the canonical Paley lemma under Lang–Waldschmidt is
     the **critical band `W_F/W -> 1/2`**, not a half-line.
   * The cap is not universal across Hadamard types. A Walsh-8 assignment allows `W_F/W -> 1` under
     the pair inequalities.
4. **Smaller points.**
   * The Lang–Waldschmidt citation needs care.
   * Theorem L's nonvanishing hypothesis is omitted from the summary.
   * Prop U concerns *linear* consequences only.
   * There is a factor-4 slip in one data remark.
   * The novelty claims should be calibrated.

Details follow. Checks are in `round4/referee_flipped_checks/` (run from that directory). Outputs
are `*_out.txt`.

---------------------------------------------------------------------------------------------

## 1. The exact reduction (Lemma 1.1, Props 2.1-2.2): correct

**Lemma 1.1, re-derived.** `sum_x H(x,a) s_xj = sigma_j [M delta_(a,a(j)) - 2 c_j(a)]`, because each
modified row reverses one term of `sum_x H(x,a)H(x,a(j)) = M delta`. Multiplying (1.1) by `H(x,a)`
and summing then gives (1.2). `theta'` cancels by column balance.

* The error is `|sum_x H(x,a) delta_x| <= M Delta/2`.
* For (1.3): `H(x,a) - H(x,b)` is `+-2` on exactly `M/2` rows, so the error is again `<= M Delta/2`.
* `d_j(a,b)` is an integer because each summand of `c_j(a) - c_j(b)` lies in `{0, +-2}`.

**Prop 2.1, re-derived.**
* With `b_3 = -2(m_a + 2n)`, the form is `Lambda_a = 2 i eta`, so `|Lambda_a| <= M Delta`.
* `|b_3| (pi/2) <= M pi + 2 pi + M Delta`, which gives `|b_3| <= 2M + 5` once `M Delta <= pi/2`.
  The pair case has `+ 4 pi`, giving `|b_3| <= 2M + 9`.
* Nonvanishing:
  * `Lambda_a = 0` forces `M arg G_a - 2 arg P_a` into `(pi/4)Z`.
  * `Gamma_a = G_a^M conj(P_a)^2` is a *positive* rational integer times the conjugate-primitive
    `Gamma'_a`, with exponent `M - 2c_j(a) != 0` at a column of label `a` with `|F_j| < M/2`.
  * The axis gap then gives a contradiction.
  * The pair case works the same way with `+-M - 4 d_j`, `|d_j| <= |F_j| < M/4`.
* Heights: `N X^2 - 2 Re(gamma^2) X + N` is primitive (a common odd prime would divide `x` and `y`).
  So `H <= 2N` and `h = (1/2) log N`.

**Independent exact check** (`ref_two_log.py`, my own code, not the author's `fcommon.py`). The
profiles have:

* arbitrary modification sets, with up to `M/2` rows per column and modified label-0 columns;
* a non-trivial common content `g = (2+i)(4+i)(1+4i)*3`, including an inert factor;
* random orientations and row units;
* Walsh-8/16 and Paley-12/20 cores.

Results:

* identity (1.4) and its pair version hold exactly in all 1,512 cases;
* the **angle bookkeeping** holds in all 1,512 cases: with real lifts and the actual `delta_x`,
  `(M Log u - 2 Log v - 2 i eta)/(i pi/2)` is an integer. The negative controls (coefficient 3
  instead of 2, or `eta` instead of `2 eta`) fail, as they should;
* `Gamma/Gamma'` is a rational integer times a unit, and `Gamma'` is off the axes and diagonals;
* 3,000 minimal-polynomial height checks pass.

The author's own checkers (`check_two_log_identity.py`, `check_sector_lemma.py`, `check_margins.py`,
`check_paley_uncertainty.py`) were rerun. Their outputs are byte-identical to the stored ones.

**Novelty calibration.**
* (1.2) is Lemma 7.2 of walsh-extraction.md multiplied by `M`. The observation that
  `sum_x H(x,a) sigma_x psi_x` is the argument of *one* Gaussian integer `P_a` is a regrouping.
* What is new is the reading of the label condition as a form in **two varying** numbers. It is
  different from item 276's `M log alpha - m log(-1)` (one varying number). It also leads to the
  sum-versus-product-of-heights diagnosis.
* I confirmed by grep that no note mentions Lang–Waldschmidt, and none reduces flipped pinnings to
  two-number forms.

## 2. Liouville / column character (Prop 3.1, Rem 3.2, Cor 3.3): correct inequality, overstated reading

**(a) (3.1) is valid but weaker than the column character.**
* `gamma_a = G_a^(M/2) conj(P_a)` is not conjugate-primitive. `G_a` and `P_a` share the flipped
  primes of label `a`.
* The half-height inequality of the column character `c = H(., a)` has
  `V_c = (M/2) W_a + W_F - 2 f_a`, where `f_a = sum_{a(x)=a} log q_x`. I re-derived
  `|T_j| = M, M-2, 2` on unflipped columns of `a`, flipped columns of `a`, and other flipped columns.
* So the true one-number inequality is
  `(M/2) W_a + W_F - 2 f_a >= (W+D)/2 + log 8 - 2 log(MC)`, which is stronger than (3.1) by `2 f_a`.
  Checked exactly on 260 labels (`ref_column_char.py`).
* "Prop 3.1 is exactly the column character" should read "is implied by the column character".

**(b) "Liouville misses by exactly `F/2 - D/4 + O(W/M)`" (Rem 3.2, summary item 3) is an upper bound.**
* The deficit at the best label is `(M/4) min_a W_a + F/2 - (W+D)/4` (minus `f_a/2`, by (a)).
* This is `<= F/2 - D/4 + W/(4(M-1))`, with equality only when `min_a W_a ~ W/(M-1)`.
* For unbalanced labels Liouville does better.

**(c) Cor 3.3: "Otherwise it gives nothing" is false.**
* Summing (3.1) over `a` throws information away. Each instance separately gives
  `W_a >= (W + D - 2 W_F)/M - (4/M) log(MC/(2 sqrt 2))` for **every** label.
* That is a balancing constraint: no label may be lighter than `(W + D - 2W_F)/M`.
* It does not exclude the balanced Prop-P profiles, so the barrier statement survives. Prop B in
  fact uses exactly this per-label bound.
* The summed bound `D <= 2 W_F + W/(M-1) + 4 log(MC/(2 sqrt 2))` is correct as stated.

## 3. Theorem L (conditional): correct as a conditional theorem

**Re-derived.**
* (L.1): `log` of `MC e^{-(W+D)/4} >= |Lambda_ab| >= c_* B^{-kappa'} (4 e^{W_a+W_b+F_ab})^{-kappa}`
  with `B <= 2M+9`.
* (L.2):
  * `sum_{a<b}(W_a+W_b) = (M-2) sum_a W_a`;
  * a row lies in `D_ab` for `P_x N_x <= (M-1)^2/4` pairs, so the average of `F_ab` is
    `<= F (M-1)/(2(M-2))`;
  * the minimum is at most the average.
* Conclusion:
  * the arithmetic `(kappa delta/2) W <= W(1/4 + 2 kappa)/(M-2) + E_M` is correct;
  * so is `W >= (M-1) log 5`;
  * the trivial case `M Delta > pi/2` is handled.
* The averaging inequality was also tested on 210 random weighted profiles (`ref_heights.py` (1)).

**Hypotheses that must be kept in any summary.**
* **Theorem L needs every nonconstant label to have a column with `|F_j| < M/4`.** This gives
  nonvanishing of all pair forms. It also guarantees `W >= (M-1) log 5`. The task summary's
  "any H, any flips" omits it. Profiles with empty labels, as allowed in Theorem B, are not covered.
* The flip condition is `F <= (1/(2 kappa) - delta) W` for a **fixed** `delta > 0`, with `M_0`
  depending on `delta`. "Covers `W_F < W/2`" means `W_F <= (1/2 - delta) W` for each fixed `delta`.
  Families with `W_F/W -> 1/2` are not covered. Section 5 below shows this is exactly the case that
  matters.

**The Lang–Waldschmidt citation needs care.**
* The commonly quoted form of the Lang–Waldschmidt conjecture is for nonzero rational numbers or
  positive integers. One example is Conjecture 2.5 of Waldschmidt's survey "Open Diophantine
  Problems" (Moscow Math. J. 2004), as reported by a web search.
* That form is `|a_1^{b_1}...a_n^{b_n} - 1| >= C(eps)^n B/(|b_1...b_n| a_1...a_n)^{1+eps}`.
* `(H_kappa)` concerns unimodular elements of `Q(i)` and `log i`. Those cannot be reduced to the
  rational case: they are angles, not real logarithms.
* I could not consult the algebraic-number formulation in the 2000 book (network access to the
  sources failed).
* So "`(H_kappa)` is a special case of the Lang–Waldschmidt conjecture" should be weakened to "the
  natural degree-two (naive-height) analogue of the Lang–Waldschmidt conjecture". Theorem L itself
  is stated conditional on `(H_kappa)`, so its correctness is unaffected.

**Height normalisation.**
* "With the absolute height the statement is false for these forms" is supported only by numerics
  (`check_lw_scale.py`) and the random-model count (`~ XY` pairs, so a minimum of order `1/(XY)`).
* It is not proved. A fixed-`G` Dirichlet argument does *not* give it: the target direction
  `arg G^(M/2) + k pi/4` always has rational slope, so approximations by `P` are only Liouville-good.
* Pigeonhole over pairs only reaches the borderline `|Lambda| << (e^{h_1} e^{h_2})^{-1}`.
* The claim should be labelled heuristic or numerical.

**Context.**
* growth/conditional.md (refereed) already derives the *full* uniform theorem from one instance of
  Vojta's conjecture.
* Theorem L reaches a much smaller conclusion, a subclass of Hadamard-core profiles, from a
  different, classical two-logarithm hypothesis.
* Its value is diagnostic. It locates the Prop-P barrier as a "sum versus product of heights"
  problem for two logarithms. It is not a stronger conditional result. The write-up should say this.

**Consistency with the concurrent `round4/hardness.md`.**
* hardness.md (Thm 3 and a closing remark) says single linear forms carry no Baker or
  Lang–Waldschmidt content. Its remark says a single small form "produces at most a product or
  progression pattern".
* Theorem L is a counterexample to the remark's generality. In Hadamard-core profiles a single
  label-pair form, whose large coefficient `M` multiplies a low-height number, excludes the class
  under Lang–Waldschmidt.
* There is no contradiction: hardness.md Prop 3.2 concerns pair angles with exponents `<= 2`. But the
  two write-ups should cross-reference each other.

## 4. Proposition B (the "sharp obstruction"): formal statement correct, interpretation overstated

**Re-derived.**
* `h(u_a) = W_a/2` and `h(v_a) = W_F/2`.
* `W_a >= log(5*13*17*29*37) = 13.986`.
* `W_F >= M log(4M/e)`, because the `k`-th split prime is `>= 4k+1` and `n! >= (n/e)^n`.
* (3.1) gives `W_a >= (V - 2W_F)/M - (4/M) log(MC/(2 sqrt 2))`.
* So `h_1 h_2 >= (W_F/2) max{(V - 2W_F)/(2M) - ..., 6.99}`.
* The pair version (`F_ab >= (M/2) log(2M/e)`, threshold `V/8`, `W_a + W_b > 27.9`) is also correct.

**The overstatement.** Prop B's display is `h_1 h_2 >= V min((1/8 - o(1)) log M, 0.87)`. For
`log M > 7` that minimum is the constant `0.87`. So the summary's "every known Baker-type bound
misses by a factor of order `c log M`" does not follow. Nor does the status paragraph's
"`h(u)h(v) >= (1/8 - o(1))(W+D) log M` in any configuration".

What the proof actually gives:

* **Corrected Prop B.**
  * If `W_F <= (1/2 - delta)(W+D)`, then
    `h(u_a) h(v_a) >= (delta/2 - o(1)) (W+D) log(4M/e)`. Proof: `W_a/2 >= delta V/M - o`, times
    `W_F/2 >= (M/2) log(4M/e)`. This was also checked on random admissible data,
    `ref_heights.py` (3).
  * The constant `1/8` corresponds to `delta = 1/4`.
  * In all cases `h_1 h_2 >= 3.49 W_F`.
  * For pairs, replace `delta/2` by `delta/4`.
* **In the regime `W_F >= (W+D)/2` there is no `log M` gain.** `ref_heights.py` (2) builds weight
  data with:
  * one label on the five smallest split primes and no flip;
  * the other labels on nearby primes;
  * flipped primes `~ M^32`.

  These data satisfy (3.1) for every label, with `W_F/V ~ 0.62`. The ratio
  `min_a h(u_a) h(v_a)/V` stays at `2.19, 2.17, 2.16, 2.16, 2.15` for `M = 44, ..., 800`, while
  `(1/8) log(4M/e)` grows. These are weight data, not a short-arc configuration. They show the
  `log M` reading is not implied by the constraints Prop B uses.
* **Qualitative conclusion: survives.**
  * A bound `log(1/|Lambda|) <= c h_1 h_2 L` excludes nothing unless `cL < 0.29`, or
    `cL < (2 + o(1))/log M` when `W_F <= V/4`.
  * All published two-logarithm constants have `cL >> 1`. Laurent's two-log bound gives a bounded
    `L` here, because `b' ~ 1/W_a + M/W_F` is small. Matveev gives `L ~ log M`.
  * The "factor `~ c log M`" holds for Matveev-type bounds through `L`. Through the heights it holds
    only for `W_F <= (1/2 - delta) V`.
* **A gap in "every Baker-type bound".** Prop B bounds the heights of the particular pair
  `(u_a, v_a)`. The same `Lambda_a` can be written in other bases of the rank-two group generated by
  `u_a, v_a`.
  * A short argument closes this approximately. Any genuine two-logarithm representation needs one
    number involving the unflipped primes of `a` and one involving the flipped primes of other labels.
  * So the product of heights is still `>= (W_a - f_a)(W_F - f_a)/4`.
  * This should be added.
* **Data remark (Section 5).** "A Baker-type bound needs `min_a W_a W_F/W` below `1/(4c)`" should
  read `1/(cL)`. From `c (W_a/2)(W_F/2) L < V/4` one gets `W_a W_F/V < 1/(cL)`.

## 5. The heavy-flip regime is (largely) empty: pair inequalities cap the flip mass

This is the main new finding of the report. It changes the description of what remains open.

**Star identity.** In any fully flipped Hadamard-core profile (one flip per row), for every row `x`,

```text
sum_{y != x} G_xy = -W + 2 l_x + 2 sum_{z != x} eps_zx l_z,      eps_zx = -H(z,a(z)) H(x,a(z)),
```

where:

* `G_xy = sum_j w_j s_xj s_yj`;
* `l_z` is the log of row `z`'s flipped prime;
* all columns carry nonconstant labels (label-0 unmodified columns are content).

*Proof.* `sum_{y != x} G_xy = sum_j w_j s_xj (S_j - s_xj)`, with column sum `S_j`.
* Unflipped columns are balanced: `S_j = 0`.
* The flipped column of row `z` has `S_j = -2 sigma H(z,a(z))`. Hence `s_xj S_j = 2 eps_zx` for
  `z != x`, and `= 2` for `z = x`.

QED. Checked exactly against literal sign matrices: 1,520 (row, profile) instances, Walsh-8/16 and
Paley-12/20, random capacity 2-4 assignments (`ref_star_bound.py`).

**Proposition S (pair inequalities cap the flips; rigorous).**
* *Hypotheses.* Some row `x0` is all ones, and every other row `z` is flipped at a label with
  `H(z, a(z)) = -1`.
* *Conclusion.* The pair inequalities `G_xy <= 4 log C - D` imply
  `W_F <= W/2 + (M-1)(4 log C - D)/2`.
* *Proof.* `eps_{z x0} = -H(z,a(z)) = +1` for all `z != x0`, so the star sum at `x0` is
  `-W + 2 W_F <= (M-1)(4 log C - D)`. The pair inequality is walsh-extraction Sec. 1, from
  `Norm gcd = e^{D + (W + G)/2}` (checked on 66 literal Paley-12 pairs with exact Gaussian gcds).
* *Where it applies.* The **canonical Prop-P / item-392 assignment** of every prime Paley core
  (`q = 3 mod 4`): finite row `x` is flipped at the diagonal label `x`, with `H(x,x) = -chi(0) - 1 = -1`,
  and row infinity is all ones.

*Consequence.* For that class:

* `W_F >= (1/2 + delta) W` together with `W >= W_F >= M log(4M/e)` forces
  `log(4M/e) <= 2 log^+ C/delta`, i.e. `M <= (e/4) C^(2/delta)`.
* So Prop 7.1, whose hypothesis is `W_F >= (1/2 + delta) W`, is vacuous for large `M` on the
  canonical class.
* Combined with Theorem L under `(H_{1+eps})` for all `eps`: **any infinite family of canonical fully
  flipped Paley clusters at fixed `C` must have `W_F/W -> 1/2`.**
* So the residual case, even conditionally, is the critical band, not "Paley with `W_F >= W/2`".

**LP evidence for general assignments (`ref_pair_lp*.py`).**
* *Method.* Maximise `W_F/W` subject to all pair inequalities in the scale-free limit `G_xy <= 0`
  (`gamma/W -> 0` at fixed `C`). Variables are the unflipped label weights and the flipped weights.
* *Exactness.* Exact rational simplex (Bland's rule). Every optimum comes with an exact dual
  certificate, re-verified, which gives `W_F <= rho W + K (4 log C - D)` with explicit `K`.
* *Results.*

  | case | assignments | `rho = max W_F/W` |
  |---|---|---|
  | Paley-12/20/24 canonical | 1 each | `1/2` exactly (`K = (M-1)/2`, the star certificate) |
  | Paley-12, random capacity 2, 3, 5 | 100 | `0.29` to `0.50` |
  | Paley-12, adversarial hill-climb, capacity 2, 3, 4, 5 | at most ~1,100 LPs | best found exactly `1/2` |
  | Paley-20, random capacity 2 and 4 | 10 | `0.26` to `0.28` |
  | Walsh-16, random | 14 | `<= 0.35` |
  | **Walsh-8, capacity 2**, assignment `[4,7,2,1,1,2,7,4]` | 1 | **`1`** |
  | Walsh-8, capacity 3 | 1 | `7/10` |

* *Coverage gaps.* The Paley-20 adversarial run (`ref_pair_lp_adv20.py`) was killed at its time limit
  before producing a result, so there is no adversarial data at order 20. The random sweep was cut by
  its time limit during Paley-24 random assignments.
* *Interpretation.* For every Paley case tested, the pair inequalities cap the flip mass at
  `(1/2 + o(1)) W`, exactly the threshold of Theorem L. For Walsh-8 they do not, but flipped Walsh
  is settled unconditionally by item 333.
* *Status.* Evidence only, apart from Proposition S. A general proof for all Paley assignments
  is open. The non-canonical certificates are not stars (inspected).

**Corrected Section 7 / summary item 5.**
* "Open even under Lang–Waldschmidt" applies to the band `|W_F/W - 1/2| <= delta + O(log^+ C/log M)`
  (rigorously for the canonical class). It also covers the heavy regime only insofar as that regime
  is non-empty, which it is not on the canonical class and is not in any tested Paley case.
* Prop 7.1 remains a correct implication.

## 6. Prop U, Lemma 6.2, sector lemma, Prop 6.4: correct

**Prop U, re-derived.**
* `lambda = H_A n/M` and the column-image formula.
* (1) Every label has an unflipped column.
* (2) `|supp H_A n| <= t + Bs`, combined with the `l2/l_inf` uncertainty `|supp n||supp H_A n| >= M`.
* (3) Uses Lemma 6.2. If `s > (q+1)/2`, the bound holds trivially.

**Lemma 6.2, re-derived.**
* Euler's criterion gives `chi(y) = y^d mod q`, including `y = 0`.
* The coefficient of `X^k` is `(-1)^k binom(d,k) sum_a n_a a^(d-k)`, with `binom(d,k) != 0 mod q`.
* An invertible Vandermonde on `supp n` follows, needing `d + 1 >= s` and `0^0 = 1`.
* Dividing by a power of `q` handles entries divisible by `q`.
* Root count: `#{x notin supp n : (Hn)_x = 0} <= d`.

**Independent check** (`ref_paley_uncertainty.py`, numpy):

* all of `{-3..3}^7` at `q = 7`;
* all 2,461,404 vectors with support `<= 6` and entries `+-1, +-2` at `q = 11`;
* 4,000 random vectors at each of `q = 19, 23, 31, 43, 47, 59, 67, 71, 79, 83`, with supports up to
  `(q+1)/2` and entries including multiples of `q` and `q^2` and values up to `10^6`.

No failure. The minimum of `#nz + s` is `(q+1)/2 + 1` throughout, so the lemma is off by one from
the truth, as item 393 is in its ternary form.

**Sector lemma (6.3) and Prop 6.4, re-derived.**
* Consecutive triangles in a convex sector have area `>= 1/2` each.
* `Norm(u) Norm(v) >= sin(2 eps)^(-2)`.
* At most `G` points of norm `< 1/(2 eps)`.
* The numerical consequence `W <= 20 log(8MC)` is right.

**Scope and novelty.**
* Prop U is about **rational linear** consequences: no rational combination isolates flipped angles.
  The summary's "the phase system provides no such pinning" is stronger than what is shown.
  Non-linear or arithmetic consequences, such as the pair inequalities used in Section 5, are not
  excluded.
* The Hadamard uncertainty `|supp n||supp Hn| >= M` is classical (Donoho–Stark type).
* Lemma 6.2 is a modest extension of item 393 (`paley_primefield_polynomial_character_gap.md`): the
  same polynomial and root-counting method, with a Vandermonde and `q`-adic normalisation so that
  arbitrary integer vectors are allowed. The write-up says this.

## 7. Section 8 and the remaining claims

* **Section 8 item 6 (general item-532 profiles have support weight `>= (1/2 - o(1)) W`):**
  re-derived. The column value is `-2 lambda(T)`. `lambda(T)` and `lambda(T xor {i})` differ by
  `lambda_i != 0`.
* **(L.4) arithmetic** (`W_F/W <= 0.2496` for `b = 5`, `M >= 6`, `1/r <= 1/25`) and the margins table:
  the rerun is identical.
* **(L.6)** is a correct reformulation.
* **The Section 7 Walsh coset remark** (five-log Lang–Waldschmidt excludes flipped Walsh with
  comparable flipped primes) is plausible. I re-derived the column values `{0, 1, sqrt M/2 - 1,
  sqrt M/2, sqrt M/2 + 1}` after scaling by `sqrt M/2`, and the support weight `~ (B+2) W/sqrt M`.
  It is only sketched, so it should be labelled a remark.
* **Novelty grep** ("Waldschmidt", "Lang", "two logarithm", "moving target", "uncertainty",
  "abc"): confirms the author's account. Only fixed-target `N log alpha - m log(-1)` forms appear in
  the notes (items 276-277, `minimal_prime_rank_width_bound.md`).

## 8. Corrected statement

1. **Exact reduction (new, correct).** For any normalised Hadamard core and any modification sets:
   * each label pinning is the nonzero two-logarithm form `M Log(G_a/conj G_a) - 2 Log(P_a/conj P_a) + b Log i`;
   * each label-pair pinning is `M Log u_ab - 4 Log w_ab + b Log i`;
   * `|b| <= 2M+9`, error `<= M Delta`;
   * nonzero when a column of the label(s) has `|F_j| < M/2` (pairs: `< M/4`);
   * naive heights `<= 2e^{W_a}` and `<= 2e^{F_a}` (pairs: `2e^{W_a+W_b}`, `2e^{F_ab}`).

   This is Lemma 7.2 of walsh-extraction regrouped. The new content is the two-varying-number
   reading.
2. **Theorem L (conditional, correct).** *Hypotheses:*
   * `(H_kappa)`, the natural `Q(i)` naive-height analogue of Lang–Waldschmidt;
   * every nonconstant label has a column with `|F_j| < M/4`;
   * `F <= (1/(2 kappa) - delta) W` for a fixed `delta > 0`.

   *Conclusion:* `M < M_0(kappa, kappa', c_*, delta, C)`. No bound on `W`, the capacity or the prime
   sizes is needed. Under `(H_{1+eps})` for all `eps` this covers `W_F <= (1/2 - delta) W` for every
   fixed `delta`, including the Prop-P class.
3. **Liouville.** (3.1) is implied by the column character, which is stronger by `2 f_a`. It forces
   `W_a >= (W + D - 2W_F)/M - O(log(MC)/M)` for every label and `D <= 2W_F + O(W/M)`. On
   balanced-label profiles it misses a contradiction by `F/2 - D/4 + O(W/M)`; in general by at most
   that.
4. **Prop B (corrected reading).**
   * `h(u_a) h(v_a) >= (W_F/2) max{(W + D - 2W_F)/(2M) - O(log(MC)/M), 6.99}`.
   * Hence `>= (delta/2 - o(1))(W+D) log(4M/e)` when `W_F <= (1/2 - delta)(W+D)`, and
     `>= 3.49 W_F` always (pairs: `delta/4`).
   * Product-of-heights bounds `log(1/|Lambda|) <= c h_1 h_2 L` with `cL > 0.29` give no
     contradiction from these forms.
   * The deficit grows like `log M` only for `W_F <= (1/2 - delta)(W+D)`, or through `L` for
     Matveev.
5. **Prop U, Lemma 6.2, Prop 6.4: correct.** They apply to rational linear consequences only.
6. **New (this report), Proposition S.**
   * For the canonical Prop-P Paley assignment, the pair inequalities give
     `W_F <= W/2 + (M-1)(4 log C - D)/2`.
   * Hence the heavy regime `W_F >= (1/2 + delta)W` is empty for `M > (e/4) C^(2/delta)`, and
     Prop 7.1 is vacuous there.
   * Under `(H_{1+eps})` for all `eps`, every infinite family of canonical fully flipped Paley
     clusters has `W_F/W -> 1/2`.
   * Exact LP certificates indicate, without proof, the same `1/2` cap for all tested Paley
     assignments (orders 12 and 20, capacities 2 to 5, adversarial at order 12). They do not hold for
     Walsh-8, where the cap is `1`.
   * The precise residual case, even conditionally, is the **critical band `W_F/W ~ 1/2`**.
7. **Unchanged.** The fully flipped phase lemma remains unproved unconditionally for every non-Walsh
   core.

## Files (`round4/referee_flipped_checks/`)

| file | content | result |
|---|---|---|
| `ref_two_log.py` | independent exact check of (1.4) and its pair version, angle bookkeeping (with negative controls), nonvanishing, minimal-polynomial heights; general modification sets, content, units | PASS, 1,512 + 1,512 + 1,512 + 3,000 |
| `ref_paley_uncertainty.py` | Lemma 6.2 | exhaustive at `q = 7, 11`, random to `q = 83` with `q`-divisible entries; PASS |
| `ref_star_bound.py` | star identity on literal sign matrices; `H(x,x) = -1` for Paley-I; gcd/pair-slack identity on literal Gaussian points | PASS |
| `ref_pair_lp.py`, `ref_pair_lp_sweep.py`, `ref_pair_lp_adversarial.py`, `ref_pair_lp_adv20.py` | exact rational LP (pair inequalities to max `W_F/W`), verified dual certificates | Section 5 table |
| `ref_heights.py` | Theorem L averaging inequality; heavy-flip weight data for Prop B; corrected Prop B bound | PASS / data |
| `ref_column_char.py` | Sec. 2(a): conjugate-primitive part of `G_a^(M/2) conj(P_a)` has log-norm `(M/2)W_a + W_F - 2f_a` | PASS, 260 labels |

Outputs are in the `*_out.txt` files. The author's checkers were rerun from `flipped_checks/` with
outputs written to the scratchpad. All seven are byte-identical to the stored outputs. The
repository `/home/user/Jarnik` was not modified.
