# Merged referee report on `round6/sharp.md`

This report merges two independent referee reports and adds the merge referee's own checks:
* `referee_sharp_arith.md` (arithmetic and class definitions; scripts in `referee_arith_checks/`);
* `referee_sharp_analytic.md` (analytic estimates and union bounds; scripts in `referee_analytic_checks/`);
* the merge referee's own checks (scripts in `referee_merge_checks/`).

Each issue raised by either referee was re-derived. Where the two referees disagreed, or where a
claim could be computed, it was recomputed with new code. That covers severity, the counterexamples
and the numbers. Each script runs from its directory with `python3 <script>`, and its output is
stored next to it as `*_out.txt`.

* `m_height.py`: height condition, per unit against aggregated. Exact check on actual Gaussian
  primes, plus a Monte Carlo of the random fibre.
* `m_misc.py`:
  * a Lemma 6.2 counterexample;
  * the §7.3 count and the Mertens factor;
  * the `mu_M` chain and its max-`omega` repair;
  * a full scan of the pair demand for every `M` in `[20, 2000]`, `A = 1, 2`;
  * the coherent example `zeta_x = n_x + i`.
* `rerun/`: the author's five scripts, rerun unchanged (§5).

`/home/user/Jarnik` and `round5/` were not modified.

---------------------------------------------------------------------------------------------

## 0. Verdict

**The main results survive.** These are Theorem 5.1 (P3-char), Corollary 5.2 and Theorem 6.4
(P3-all). Both referees re-derived every step they rest on, including every step of `lower.md` that
`sharp.md` uses, and found them correct. Every point of disagreement was re-checked here, and none
changes this conclusion.

**One major error, which can be repaired.** It is in the sharpness part (§6.4) and its Verdict
item 5.
* The "necessary height condition" is stated in aggregated form, `m_v <= sqrt2 e^(w(v)/2)`. That
  form is **false for actual clusters**: 8976 of the 8977 Gaussian primes below `2*10^5` violate it.
* The proof of Prop. 6.5 concludes from a violation of the aggregated form, so as written it proves
  nothing.
* A per-unit argument restores Prop. 6.5 with the same threshold. It then holds for every `A`, but
  only for the uniform random fibre.

The arithmetic referee rated this minor and the analytic referee major. I rate it **major**: it is a
false mathematical statement, advertised in the Verdict, inside the statement of a numbered
proposition, and built into the evidence script. The repair is short and the conclusion survives.

Everything else is minor: quantifier and hypothesis slips, overstated scope, wording, and numeric
slips in heuristic paragraphs. After merging duplicates, 14 minor items remain (issues 2-15 in §2;
issue 15 bundles four small wording points).

| item | status after merging |
|---|---|
| Def. 1.1; Prop. 1.2 (actual clusters lie in `F^all_X` for every `X`) | correct (both referees; exact check on 15 genuine clusters) |
| Prop. 1.3 (`two_thirds` inside `F^char_X`, `X >= 2M`) | correct. `two_thirds` uses the configuration only through (2.5) and (3.6), and (3.6) only at `m_(q^a) < M`, i.e. `q^a < 1.25M` |
| `F^char_(3M) ⊂ III-small` | correct, relative to the informal definition in `lead6/notes.md` §1, §6. Verdict item 1 states the inclusion backwards (issue 15) |
| Prop. 2.1, Cor. 2.2 | correct. Exhaustive over all 184 odd `Q <= 1000`; distinct twin columns are not needed for (2.2) |
| Lemma 3.1 | correct (multiplicity `<= 4` exactly to `2^27`, `<= 9` by Rosser-Schoenfeld beyond) |
| Def. 3.2, Lemma 3.3 | correct, except the written demand chain (issue 4); the claim itself is true |
| `lower.md` Lemma 2.1, 2.2, S; (3.1); Lemma 3.1; Prop. 3.3; Lemma 5.3; Prop. 1.2, 1.3 | correct (re-derived by both referees and here) |
| `lower.md` Lemma 5.2 = `sharp.md` Lemma 5.1b | needs an extra hypothesis (issue 3); harmless in every use |
| Theorem 5.1 | correct. `M_0` is astronomically large with the proved constants (issue 5) |
| Corollary 5.2 | correct |
| Lemmas 6.2, 6.3; Theorem 6.4 | correct (Lemma 6.2's parenthetical needs `e = 1`, issue 11; step 2 wording, issue 14) |
| Prop. 6.1 | correct as stated, per unit |
| Prop. 6.5 | **wrong as written** (major issue 1); true after the per-unit repair |
| Prop. 6.6 (sufficient for "yes") | correct at `X = 3M`: the trivial bound `m <= e^(psi(3M))` covers non-pairs and off-unit pairs |
| Prop. 7.1 | correct, given `v(c) != 0` (issue 10) |
| Prop. 7.4 | needs its pair hypothesis for all units and the same-level clause (issue 6) |
| §7.2/§7.4 "coherent designs fail" | overstated (issue 7) |
| §6.6, §7.3 heuristics | numeric slips (issues 8, 9); qualitative conclusions stand |
| Appendix A | a sketch, presented in the Verdict and Appendix B as established (issue 12) |

---------------------------------------------------------------------------------------------

## 1. The major issue

### Issue 1 (major). The height condition is per unit; the aggregated form is false

*Location:* Verdict item 5; Prop. 6.5 (statement and last bullet of the proof); §8 row for
`check_p3all.py`; the `check_p3all.py` docstring ("`m_(e_j) <= sqrt(2 p_j)` is necessary") and its
Parts A/B.

* Prop. 6.1 is correct per unit: `m_(v,s) <= nu_s^(-1) e^(w(v)/2)` for each `s`. For actual
  clusters this is `m_(v,s) | D_s` with `D_s in {Y, X-Y, X, X+Y}` of `A_v = X + iY`, and
  `|D_s| = nu_s^(-1) |A_v| |sin(arg A_v - s pi/4)|`.
* The aggregated modulus `m_v = prod_q q^(max_s a_(q,s)(v))` collects the odd parts of all four
  `D_s`. These are pairwise almost coprime, so for actual data `m_v` is about `|XY(X^2 - Y^2)|`,
  of order `e^(2w)`.
* **Merge check H1** (`m_height.py`, independent of both referees):
  * Among the 8977 primes `p = 1 mod 4` below `2*10^5`, 8976 violate `m_v <= sqrt2 e^(w/2)` (only
    `p = 5` does not), and none violates Prop. 6.1.
  * With moduli truncated to `q <= 3000` the counts are unchanged.
  * For 4000 two-prime monomials there are 4000 aggregated violations (3829 at `q <= 3000`), and
    no per-unit violation.
  * Example: `pi = 2 + 3i` has `D_s = 3, -1, 2, 5`, so `m_v = 15 > sqrt26 = 5.10`.
* So "some `e_j` has `m_(e_j) > sqrt2 e^(w_j/2)`, so the height condition fails" contradicts
  nothing. As written, Prop. 6.5 is not proved.

**Repair (verified).** Fix a non-twin column `j`.
* For `q in Q_0` (`q = +-3 mod 8`, so `m_q = 4 mod 8`, `i = g^(+-m/4)` and `m/4` odd), the label
  `ell_j(q)` is uniform on the `m_q/2` classes of parity `nu_j(q)` (Lemma 6.3(2)).
* So `e_j` collides at level 1 with the unit `i^s` with probability `2/m_q` when
  `s = nu_j(q) mod 2`, and with probability 0 otherwise.
* Put `lambda_(j,s) = sum 2/m_q` over the `q in Q_0` with `nu_j(q) = s mod 2`. Then
  `lambda_(j,0) + lambda_(j,1) = lambda/2`, so some `s_j in {0, 1}` has
  `lambda_(j,s_j) >= lambda/4`, a positive constant.
* The Poisson lower bound for `kk = floor((1 - eta_0) log M/loglog M)` hits is still
  `exp(-(1 - eta_0 + o(1)) log M)`. The constant mean enters only `O(kk) = o(log M)`.
* Independence over columns is unchanged, so w.h.p. some `(j, s_j)` has
  `log m_(e_j, s_j) >= (1 - eta_0)^2 (1 - o(1)) log^2 M/loglog M > tau/2 + 1`. This is a genuine
  violation of Prop. 6.1.
* With `Q_0 = {q = +-3 mod 8 : X^(1-eta_0) <= q <= X}` each hit contributes
  `>= (1 - eta_0) A log M`. The threshold then becomes `(2A - eps) log^2 M/loglog M` for
  `X = 3M^A`, so Theorem 6.4 is sharp for this method at **every** `A`, not only at `X = 3M`. The
  analytic referee observed this; I re-derived it.

**Size of the error** (merge check H2, Monte Carlo of the fibre model, `b = 33`, `X = 3M`, 3 seeds).
The table gives the least `tau` permitted by the samples, in units of `log^2 M/loglog M`.

| `M` | per-unit (correct) | aggregated (as written) |
|---|---|---|
| `10^3` | 3.47-3.81 | 6.38-6.96 |
| `3*10^3` | 3.43-3.53 | 6.39-6.74 |
| `10^4` | 3.34-3.66 | 6.14-6.49 |

This agrees with the analytic referee's exact DP (3.48 against 6.42 at `M = 10^3`; 3.17 against
5.70 at `10^6`). The aggregated quantity overstates the obstruction by about 1.8x at every
accessible `M`. Both columns tend slowly to 2.

**Severity.** The two referees agree on the substance, the repair and the threshold. They differ
only in the label. I keep **major**, for three reasons:
* the Verdict calls a false inequality "the necessary height condition (Prop. 6.1)";
* the statement of Prop. 6.5 contains the invalid inference;
* the §8 evidence for this section measures the wrong quantity.

The result survives after the repair.

---------------------------------------------------------------------------------------------

## 2. Minor issues (merged; duplicates combined)

2. **Prop. 6.5 scope** (analytic). Prop. 6.5 is proved for the uniform random fibre of Def. 6.3
   only.
   * Its witnesses are non-twin columns. Their labels are free inside their parity class, so
     choosing them away from the 4-torsion removes every single-prime witness at no cost.
   * The `2M` twin columns would still carry Poisson tails, and nobody has analysed a coordinated
     choice.
   * So Verdict item 5's "with per-modulus independent residue choices the height condition fails
     w.h.p." is broader than what is proved. "For the uniform random fibre" is accurate. §6.4's own
     closing sentence ("a property of the method, not of the class") is accurate.
3. **Lemma 5.1b and `lower.md` Lemma 5.2: restriction gap** (both referees; re-derived).
   * Suppose a deleted constant column has `p_l <= X'`. Then `p_l` does not divide `N'`, so the
     restricted model needs point residues at `p_l` and the full condition (1.1) there, for every
     character and every unit.
   * The original model carried only level residues at `p_l`, with the same-level **pair** clause.
     "The collision modulus can only drop" fails there.
   * In `lower.md` the step "multiply by `beta` of norm `N_new/N`" is not even defined mod `p_l`.
   * It is harmless in every use. In Theorems 5.1 and 6.4 every `p_j > X = 3M'^A >= X'`.
   * **Fix:** add the hypothesis "every prime of a deleted column exceeds `X'`" (or simply "all
     `p_j > X`").
4. **Lemma 3.3(1): the demand line is derived through a non-monotone bound** (both referees;
   re-checked).
   * The per-`d` form `(10 + 26/loglog d) log d` is 315, 124 and 104 at `d = 3, 4, 5`.
   * Merge check MU: its maximum over `d < M` exceeds `mu_M log M` at every tested `M`
     (113.5 against 35.5 at `M = 16`; 28.6 against 20.8 at `M = 60000`, in units of `log M`).
   * **Correct route:** `Lambda'(d) + log d <= 10 log d + 18.72 omega(d)`. Then use
     `max_(d<M) omega(d) <= 1.3841 log M/loglog M` for `M >= 16`. This follows from Robin's bound,
     `omega(d) <= 2` for `d < 30`, and the monotonicity of `log x/loglog x` on `[e^e, inf)`; MU
     also checks it for every `M <= 60000`, with minimum slack 0.116. This gives
     `log m_(xy,1) <= (10 + 26/loglog M) log M`, so the claim is true.
   * The true maximum of `Lambda'(d) + log d` is `6.2-6.7 log M`.
   * Nit (arithmetic referee): in Def. 3.2 the structured levels run over
     `2 <= a <= min(a*(q), A_q)`.
5. **Theorem 5.1 / 6.4: `M_0` and parameter ranges** (analytic; re-derived). This is not an error,
   since the statements are asymptotic. With the proved constants (`b = 33`, `C = 1`):
   * the pair family needs `mu_M < 12.5`, i.e. `loglog M > 10.4`. Computing its exponent
     `(2 + mu_M) log M + log(2(2 + pi/C)) + 1/4 - (b-4) tau/4` with the smallest admissible
     `tau`: it is negative for `log M < 16.1`, positive for `16.1 <= log M <= 3.2*10^4`, and
     negative from then on. The analytic referee gave the band as `50 <= log M <= 3*10^4`; the
     computed band is wider;
   * the ternary family needs `log M >~ 450`;
   * in Theorem 6.4, `beta > 0` only for `log M >= 4612` (`A = 1`).

   The "explicit" bound at Paley orders holds only for `M >= M_0`, and this should be said next to
   it.
6. **Prop. 7.4: the pair hypothesis must cover every unit** (arithmetic referee; extended here).
   * The conclusion is about `F^char`, where (1.1) holds for all four units.
   * With moduli up to `N^A`, the trivial bound used in Theorem 5.1 for off-unit pairs is no longer
     available. A unit-1 bound does not control `eta = -1, +-i`.
   * Merge check K: the coherent design `zeta_x = n_x + i` with `n_x^2 + 1` prime `> X`. At
     `X = 10^8` the unit-1 pair moduli are `<= e^(4.84)`, while the other units reach `e^(18.47)`.
   * **Fix:** require `log m^gen_xy <= C_1 log M`, the maximum over all units, or
     `log m^off_xy <= W/4 - 1`. Since all `p_j <= X` in this regime, the same-level pair clause
     at the `p_l` must be bounded too.
7. **"Coherent ones fail by §7.2" (§7.4; §0 item 6) is overstated** (arithmetic; re-checked).
   * §7.2 shows only that the a-priori bound `m^gen <= Y^(2n)/4` of `lower.md` Lemma 4.1
     certifies nothing once `4 log X` exceeds the average pair slack `W/(2(M-1))`. For general
     characters this happens once `A >= 1/4`.
   * It says nothing about the actual moduli of a coherent design. In check K, unit-1 pair moduli
     are `<= e^(4.8)` against the bound `4 log Y = 73.9`.
   * Failure through units other than 1, and through characters of support 4, is plausible but
     heuristic.
   * **Corrected sentence:** "the Lemma 4.1 bound certifies nothing for large `X`; whether some
     coherent design works is not decided here."
8. **§7.3 numbers** (arithmetic; re-checked in C).
   * `#{c : ||c||_1 <= 4M}` is `e^((3.085+o(1))M)`: computed 3.054, 3.080, 3.084 at
     `M = 10^2, 10^3, 10^4`, with limit 3.0846. The value 5.27 is the crude
     `2^n binom(M+n-1, n)`.
   * The per-collision factor is `(log 2)^kk`, not `(log 2/y)^kk`: Mertens gives
     `sum_(e^y<q<=e^2y) 1/q = 0.685, 0.690` at `y = 6, 8`.
   * Corrected: `kk ~ 3.1M/log M`, with demand about `3.1 M y/log M`.
   * Against a heuristic slack `M^(3/2) log M`, the heuristic failure point moves to moduli
     `e^(M^(1/2+o(1)))`. That is still inside `[M^A, N^A]`, so the qualitative claim stands.
9. **§6.6:** `pi(X) log 2` is `~ 3M log 2/log M` at `X = 3M` (and `~ 3M^A log 2/(A log M)` at
   `X = 3M^A`), not `M log 2/log M` (arithmetic).
10. **Prop. 7.1 needs `v(c) != 0`** for (1.1) to apply (both referees). This is automatic on the
    §4 profile (twins, `rank = M`); state it for general profiles.
11. **Lemma 6.2's parenthetical "or any profile with `2M` distinct twin columns" needs `e = 1`**
    (arithmetic; corrected here).
    * The last step uses `|(c'^T a)_j| <= ||c'||_1 e_j/2 < G/2`, which needs `e_j = 1`.
    * The arithmetic referee's `M = 2`, `G = 10`, `||v||_1 = 2` instance is not valid as stated.
      With `M = 2`, normalised twin columns force exponents `>= 2` and `||v||_1 >= 6`.
    * A valid counterexample (merge check G):
      * `M = 2`, with columns `(2,0), (1,0), (0,2), (0,1), (100,0)`;
      * `v = a_1 - a_2 - 100 e_5 = (2, 1, -2, -1, 0)`;
      * then `v in (L + 100 Z^5) \ L` and `||v||_1 = 6`, while `100 > 12`.
    * Every use in `sharp.md` has `e = 1`, so nothing downstream changes.
12. **Appendix A is a sketch, presented as established** (analytic). The appendix itself says "a
    referee pass is still owed". Verdict item 7 ("also works for III-small and for P3-char") and the
    Appendix B rows should say *sketch*.
    * The Green-Sanders use is correct: Boolean theorem via level sets, norm `sum |fhat|`, maximum
      level-set norm 225, 12 copies, cosets versus subgroups a factor 2.
    * The Assembly union bound must run over **multiset classes**. Taken over characters with
      `Mult >= binom(M, ceil(s/4))` it diverges.
    * `m_c <= e^(2.6M)` holds for III-small moduli (`q^a < 2.5M`, `psi(2.5M) < 2.6M`). It fails
      for P3-char at `X = 3M^A`, where the accessible product is `e^(psi(3M^A))`. So "works for
      P3-char" would need Theorem 5.1's theta-moment design.
13. **Evidence wording** (arithmetic; re-checked in D). Verdict item 3 says "measured `<= 6.8 log M`
    up to `M = 16384`". That holds at the five sampled `M` only. The full scan over every
    `M in [20, 2000]` gives a maximum of `6.999 log M` (`M = 111`, `A = 1`) and `7.157 log M`
    (`M = 23`, `A = 2`). This is evidence only; the proof uses `mu_M`.
14. **Theorem 6.4 step 2 wording** (analytic). (1.2) for `v in L` is justified via `m^gen`, which
    pairs do not have. For pairs, `s = 0` follows from `Fail_xy` (the `nu_0 = 1` form) and
    `s != 0` from `Off_xy`, using `sin(pi/8) > e^(-1)`. The argument is complete; only the
    citation is wrong.
15. **Small wording points** (both referees).
    * Verdict item 1, "This contains the task's class III-small", reverses the inclusion. It should
      read `F^char_(3M) ⊂ III-small`, as the table and §1.1 say.
    * Prop. 6.5 has `(b-2)M - b` non-twin columns, not "at least `(b-2)M`".
    * Scope of "exactly one bit per prime and modulus" (§0 item 2, Remark 2.3(2), Cor. 5.2). This
      is exact for the norm-only class. Genuine residues also satisfy quartic reciprocity, so
      `ell_j mod 4` is data at the modulus `p_j`; the arithmetic referee checked this for
      `q = 7 mod 8`. For `X < p_j` (Theorem 5.1) this is invisible. In the §7 regime it is not
      encoded, and one sentence should say so.
    * References: Rosser-Schoenfeld (1962), Robin (1983, `omega(n) <= 1.38402 log n/loglog n`),
      Siegel's lemma (for example Bombieri-Gubler, Lemma 2.9.1) and Burgess are named without
      precise citations. `sharp.md` makes no novelty claims, and Prop. 2.1 is rightly not claimed
      as new.

**Disagreements between the referees and how they were resolved.**
* *Severity of issue 1:* major (see §1).
* *Lemma 6.2 counterexample:* the arithmetic referee's instance was wrong in its details but right
  in substance; it is replaced by check G.
* There were no other disagreements. Each referee raised some issues the other did not, and all of
  them were confirmed:
  * arithmetic referee only: issues 6-11 and 13;
  * analytic referee only: issues 2, 5, 12 and 14.

---------------------------------------------------------------------------------------------

## 3. Dependence on `lower.md` and other inputs

**Steps of `lower.md` used by the proofs in `sharp.md`.**
* Definition 1.1(a) and (e) (normalisation, same-level pair clause).
* Prop. 1.2 and Prop. 1.3.
* Lemma 2.1 (items 1-6), Lemma 2.2 and Lemma S (= `sharp.md` Lemma 4.2).
* The profile of §3.1 with identity (3.1), and Lemma 3.1 (= Lemma 4.1).
* Prop. 3.3 (= Prop. 4.3).
* Lemma 5.3 (= Lemma 5.1a) and Lemma 5.2 (= Lemma 5.1b).

Both referees re-derived each of these line by line, with exact checks. I re-read Lemma 2.1, 2.2, S,
Prop. 3.3, Lemma 5.2, 5.3 and Prop. 1.3 against `two_thirds/proof.md` (2.5), (3.6) and (5.1).
* All are correct except Lemma 5.2, which needs the hypothesis of issue 3.
* In Lemma S, the `kappa = 2` window gives `|T| > 5M/16` and `phi > M/8` each, so `F > 5M/4`.
  The small range is fine for `M >= 44`.
* `lower.md` Lemma 4.1, Theorem 5.1, Cor. 5.4 and Prop. 8.1 are cited only in discussion (§6.6,
  §7.2, Remark 2.3) and carry no proof.

**So `sharp.md` does not depend on any unrefereed part of `lower.md`.** Its dependence on
`lower.md` is exactly the list above, re-verified here, with Lemma 5.2 restricted to `p_j > X`.

**Other inputs.**
* `outside.md`: Lemma 0.1 and Theorem A.1(ii). The latter needs character admissibility on `L` and
  `Pi < 1.75`, which `p_j >= 21 r^2` gives.
* `referee_outside.md` §3.1 (Haar-positivity in the full torus).
* `two_thirds`: Lemma 1.1, Lemmas 3.1-3.4, (2.5), (3.6) and Theorem 7.5.
* Analytic number theory:
  * PNT for primes `= 1 mod 4` in `[Z, 2Z]`, and for primes `= 3 mod 4`, which gives Paley orders
    `M' = (1+o(1))M`;
  * Mertens in progressions mod 8;
  * Rosser-Schoenfeld (`pi`, `psi`);
  * Robin's bound on `omega`.
* The inclusion in III-small is relative to the informal definition in `lead6/notes.md` §1, §6.

---------------------------------------------------------------------------------------------

## 4. Corrected statement of what `sharp.md` proves

**Classes.** Fix `M`, `C > 0`, `X >= 3`, and for odd primes `q` put
`A_q = max{a >= 0 : q^a <= X}`. A member of `F^char_X(M, C)` is
`(p, e, (a_x), phi, theta_0, delta, k, r)` with:
* (a) a normalised profile: distinct `p_j = 1 mod 4`, distinct rows, `min_x a_xj = 0`,
  `max_x a_xj = e_j`;
* (b) `<2a_x - e, phi> = theta_0 + delta_x + (pi/2) k_x` with `delta in [0, C e^(-W/4)]^M`;
* (c) both Lemma 0.1 gaps for every `v != 0`;
* (d) for every odd prime `q <= X` and every `j` with `p_j != q`, a unit `r_j(q) mod q^(A_q)` with
  `Norm r_j(q) = p_j`. Point residues
  `rho_x(q) = i^(-k_x) prod_j r_j^(a_xj) conj(r_j)^(e_j - a_xj)` for `q` not dividing `N`, and level
  residues (the product over `j != l`) at `q = p_l`;
* (e-char) (1.1) for every zero-sum `c` with `v(c) != 0` and every unit `eta`:
  `|sin((c.delta - arg eta)/2)| >= 2^(-1/2) m_(c,eta) e^(-w(v(c))/2)`, with the actual collision
  modulus of `prod rho^c = eta` at `q <= X`, `q` not dividing `N`. Plus the same-level pair clause at
  the `p_l <= X`.

`F^all_X(M, C)` adds (e-all): (1.2) `|sin(<v,phi> - s pi/4)| >= nu_s m_(v,s) e^(-w(v)/2)` for every
`v != 0` and every `s in Z/4`. Here `nu_s = 1` for even `s` and `2^(-1/2)` for odd `s`, and
`m_(v,s)` is the per-unit collision modulus of `prod_j (r_j/conj r_j)^(v_j) = i^s`.

**What is proved** (all asymptotic in `M`, for fixed parameters):

1. **Validity.**
   * Every normalised actual cluster on an arc of length `<= C sqrt R` lies in `F^all_X(M, C)` for
     every `X`.
   * `F^all_X ⊂ F^char_X ⊂ A_X` (`lower.md`), and `F^char_(3M) ⊂ III-small`.
   * For `X >= 2M`, every member of `F^char_X(M, C)` satisfies
     `3 M log M <= W + (14 + 4 log^+ C) M`.
2. **Structure.**
   * `T_(q^a)` is cyclic of order `(q - chi(q)) q^(a-1)`.
   * Factored data force exactly `ell_j = nu_j(q) mod 2`, with `nu_j(q) = [(p_j/q) = -1]`, and
     nothing else among norm-only data.
   * `-1` is always a square, and `i` is a square iff `q = +-1 mod 8`.
   * With twins, a point class map is realisable iff `kappa - sigma_q` is constant mod 2, with the
     explicit preimage (2.2). Distinct twin columns are not needed for this.
3. **Theorem 5.1 (P3-char).** Fix `A >= 1`, `b >= 33`, `C > 0`. There is `M_0(A, b, C)` such that
   for all `M >= M_0`, `F^char_(3M^A)(M, C)` has a member with `e = 1`, all `p_j > 3M^A`, and
   `W <= (b max(A,2) + o(1)) M log M`.
   * At prime Paley orders `M = q_0 + 1 >= M_0`:
     `W <= b(M-1)(max(A,2) log M + 2 loglog M + log(42 b^2)) + 1`.
   * The pair step uses `log m_(xy,1) <= (10 + 26/loglog M) log M`, via the max-`omega` route.
   * With the proved constants `M_0` is astronomically large (pairs need `loglog M > 10.4`).
4. **Corollary 5.2.** An argument deriving `M <= f(R)` only from (I) the exponent profile and
   identities, (II) the Lemma 0.1 gaps, (III_F) (1.1) for some factored residue assignment of the
   genuine norms at odd prime powers `<= 3M^A`, and (IV) Haar-generic angle properties, has
   `f(R) >= (2/K - o(1)) log R/loglog R` along a sequence, with `K = 33 max(A, 2)`. For `A <= 2`
   the best constant from these inputs lies in `[1/33, 2/3]`.
5. **Theorem 6.4 (P3-all).** Fix `A >= 1`, `b >= 33`, `C > 0`. For all large `M`,
   `F^all_(3M^A)(M, C)` has a member with `e = 1` and
   `W <= (2Ab + o(1)) M log^2 M/loglog M`, where `o(1) = O(logloglog M/loglog M)`. The parameters
   are defined for `log M >~ 4600`. Hence such arguments have
   `f(R) >= (1/(Ab) - o(1)) log R logloglog R/(loglog R)^2` along a sequence.
6. **Prop. 6.5 (corrected).** Fix `eps > 0` and `A >= 1`. Take the construction of Theorem 6.4 with
   the **uniform random fibre**, `X = 3M^A`, and `tau <= (2A - eps) log^2 M/loglog M`. Then with
   probability `-> 1` some non-twin column `j` and some unit `s` have
   `m_(e_j, s) > nu_s^(-1) e^(w_j/2)`. This violates the necessary per-unit height condition
   (Prop. 6.1), so no `delta`, `phi` completes the data. Theorem 6.4 is therefore sharp for this
   method at every `A`. The aggregated inequality `m_v <= sqrt2 e^(w/2)` is **not** necessary and
   must not be used.
7. **Lemma 5.1b** holds under the added hypothesis that every deleted prime exceeds `X'`. This holds
   in all uses.
8. **Moduli up to `N^A`: nothing is proved.**
   * Prop. 7.1 (fixed labellings collide everywhere) holds for `c` with `v(c) != 0`.
   * Prop. 7.4 holds with the pair hypothesis taken over all units and the same-level clause.
   * §7.2 shows only that the `lower.md` Lemma 4.1 bound is useless for large `X`.
   * §7.3 is a heuristic, corrected to `e^(3.085M)` characters and Poisson mean `log 2`, which puts
     its failure point at moduli `e^(M^(1/2+o(1)))`.
9. **Appendix A** (Sylvester route) is an unrefereed sketch. It plausibly covers III-small, with a
   union over multiset classes. It does not cover P3-char as written.

None of the corrections changes Theorem 5.1, Corollary 5.2, Theorem 6.4, the constants
`33 max(A,2)` and `66A`, or the bracket `[1/33, 2/3]`.

---------------------------------------------------------------------------------------------

## 5. Merge referee's checks (`round6/referee_merge_checks/`)

| script | what it checks | result |
|---|---|---|
| `m_height.py` H1 | per-unit and aggregated height on all 8977 primes `p = 1 mod 4` below `2*10^5` and 4000 two-prime monomials, moduli unrestricted and `q <= 3000` | per-unit: 0 violations; aggregated: 8976/8977 and 4000/4000 (3829/4000 at `q <= 3000`) |
| `m_height.py` H2 | Monte Carlo of the random fibre (`b = 33`, `X = 3M`, all levels): least admissible `tau`, per-unit against aggregated | 3.3-3.8 against 6.1-7.0 (units `log^2 M/loglog M`), `M = 10^3 ... 10^4` |
| `m_misc.py` G | Lemma 6.2 with general exponents | counterexample at `M = 2`, `G = 100 > 2||v||_1 = 12` |
| `m_misc.py` C | `l^1`-ball count; Mertens factor | exponent 3.0846; `sum 1/q = 0.685, 0.690` |
| `m_misc.py` MU | Lemma 3.1(3) to `6*10^4`; per-`d` chain against max-`omega` chain against `mu_M`; Robin max-`omega` form | per-`d` chain fails at every tested `M`; max-`omega` chain holds; minimum slack 0.116 |
| `m_misc.py` D | worst-case pair demand of `D_A`, every `M in [20, 2000]`, `A = 1, 2`, plus `M = 4096, 16384` | max 6.999 (`M = 111`) and 7.157 (`M = 23`, `A = 2`); reproduces 5.734 ... 6.535 at the author's `M` |
| `m_misc.py` K | coherent `zeta_x = n_x + i`, `X = 10^6, 10^8`, `M = 20` | unit-1 pair `log m <= 4.84`; other units up to 18.47; bound `4 log Y = 73.9` |
| `rerun/` | the author's five scripts, unchanged (`summary.txt`, `timings.txt`) | `check_parity`, `check_lemmaS6`, `check_p3all`, `check_realise`: output identical to the stored one; all PASSED. `check_design`: identical except that Part 1 now runs to `q < 2^25` (default `KMAX`) instead of the stored `2^23`, with the same maximum multiplicity 4; ALL DESIGN CHECKS PASSED |
