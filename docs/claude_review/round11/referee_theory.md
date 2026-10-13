# Referee report on round11/theory.md

Referee: adversarial. Scripts in round11/referee_theory_checks/ (written independently; author code not imported).

## 0. Verdict

**The headline does not survive as stated. Hypothesis SL(n) is false, so Theorem 3.3 is vacuous.** Every
proof in the write-up is correct as an implication, and the parts that are not about SL survive with wording
fixes: Lemmas 1.1-1.2, Props 2.1-2.5, Theorem 3.1, Cor 3.2, Props 3.4, 3.5 and 4.1. The new mechanism (exact
`pi`-adic valuation of `F_c(z)` on a box widened in one direction) is also correct. But the hypothesis that is
supposed to replace ZE fails for two independent, elementary reasons. A repair (grouped widening, §7) gives a
hypothesis SL\* that these arguments do not refute. SL\* is open.

1. **(Fatal for Thm 3.3.) The barrier profiles refute SL(n)** (Proposition R1, §4.1). SL(n) quantifies over all
   admissible profiles, and nothing in the definition rules out duplicated or twin columns.
   * If column `p` has a duplicate `q`, or a twin `q` with `a_p - a_q` supported on at most two non-base points,
     then the slab of `q`, which stays narrow in `A_p`, caps `l_p`: `h_{A_p}(a_p) <= w + 4T < m`. The level-sum
     lemma (the author's own Prop 3.5(i)) then kills `UM_p`.
   * The round-5 Paley fake profile `P(H_q, pi, 8, a*)` has eight identical copies of every label plus one-row
     twins. The skew double core has two-row twins `S'_a = S_a xor {0,a}`.
   * Both are admissible in the regime: rank `n - 1`, pair separation `> 1/2`, and `n - 1 >= 300 ln 4r` from
     `M = 3528` and `M = 3068` respectively.
   * They fail SL at all but at most 1 (resp. 2) columns, against the 8 that SL(n) requires. This holds at every
     Paley order `q + 1` (`q = 3 mod 4` prime) above those thresholds, which is an infinite set.
   * Hence "SL(`n'`) for all `n' >= n_1`" fails for every `n_1`, and Theorem 3.3 is true only because its
     hypothesis is never satisfied.
   * The heuristic after Prop 3.5 ("(ii) is not expected to obstruct admissible profiles") is wrong. Pair
     separation does not prevent duplicate or twin columns, and the barrier constructions are built from them.
2. **(Major.) The parity of `m`** (Lemma R2, §3.2). On `Z^n_0` the widened slab `|<k,beta_p>| <= m/e_p` is the
   same as `|l_p(k)| <= m/2` with `l_p` an integer.
   * If `m` is odd, `l_p` takes at most `m` values on `A_p`. The level-sum lemma then excludes `UM_p` for
     **every** profile and every `p`, so SL(n) is false at every such `n` where an admissible profile exists.
   * *Exact:* `m(n)` is odd for 4076 of the 8137 values `n in [64, 8200]`. This includes `n = 2137`, `3068`,
     `3528`, `4096` and `8192`, and 130 of the 269 Paley orders in the regime of the round-5 profile up to 8200.
   * If `m` is even, `UM_p` is still possible only at the bottom level `l_p(k*) = -m/2` (new Lemma R2, proved).
   * So the claim "on `A_p` the level obstruction of round-10 Prop 5.5(ii) disappears" is false for odd `m`, and
     true only at one level for even `m`.
   * Theorem 3.1 is correct, but for odd `m` it carries no information (every profile satisfies its conclusion).
3. **(Major, evidence.) No evidence model uses the box the theorem uses.**
   * `had_sl.py` widens to `|<k,beta_p>| <= 2m`, which is `2m + 1` levels at LP cost `2m log p`. Theorem 3.1 as
     stated does not license this. It would need `log p <= W/16`.
   * The Boolean models use `l(t) = |t|` with `M + 1 >= 9 > m + 1` levels.
   * My exact rerun on the faithful box (Hadamard `n = 8`, `T = 1`, `|l_p| <= m/2`):
     * `m = 3` and `m = 5`: no `UM_p` point at all;
     * `m = 4`: `UM_p` holds exactly at the 24 bottom-level points.
   * The author's "failed at `l* = -1, 0`" came from a sufficient-only restricted-support test. In the author's own
     `2m` box an exact full scan finds `UM_p` points at levels `-4, -3, -2` and one at level `0`.
4. **(Minor.)**
   * Prop 3.4(ii)-(iii) are equivalent forms of SL-full (`V(S,m)`), not of SL (`V_0`).
   * Lemma 1.1's "answers the referee's worry" applies to `H = e^{mW/40}`, not to the `H_0` used in SL. At `H_0`,
     Lemma 1.1 only gives codimension `<~ E/2`.
   * Prop 2.5's "decided by higher-order terms" is imprecise: near a cheap stratum *every* order is `O(eta)`.
   * The `had_sl` reading "held exactly at level `-2`" refers to 10 sampled points only.

**What survives (corrected statement, §8).**
* Lemmas 1.1-1.2 and Props 2.1-2.5, as stated (wording fixes).
* Theorem 3.1 as stated, with the reading that it says nothing when `m` is odd.
* The referee's grouped version, Theorem R3, proved by the author's argument. Widen a set `Q` of columns by
  factors `lambda_q`, with `sum_{q in Q} lambda_q log q <= W/8`. Then actual clusters in the regime admit no
  `UM_p` vector in `V_0(A_{Q,lambda})` for any `p in Q`.
* A conditional growth theorem from the corresponding Hypothesis SL\*(n), which neither refutation touches. SL\* is
  open, and it must first be tested on the barrier profiles.

No growth or uniformity statement is proved, and none is claimed by the write-up.

## 1. Lemmas 1.1-1.2 (large-height Siegel; dichotomy on general supports)

**Verdict: correct. One interpretive caveat.**

*Lemma 1.1.* I re-derived the bound.
* Pigeonhole: `M_0 {0..H}^U` lies in `[-UHB, UHB]^{E_S}`, so some fibre has `>= (H+1)^U/(2UHB+1)^{E_S}` points.
* A `d`-dimensional subspace meets `[-H,H]^U` in at most `(2H+1)^d` points.
* With `L = ln(2H+1)`: `ln(H+1) >= L - ln 2`, and `ln(2UHB+1) <= L + ln(UB)` (as `UB >= 1`). This gives
  `d >= U - E_S - (U ln 2 + E_S ln(UB))/L`, and `dim V(S,m) = U - E_S`.
* Constants:
  * `40 ln 2 = 27.73` (stated `27.8`);
  * `ln(UB) <= (1 + 0.0261) m ln n` at `n = 64`, using `ln(2n+1)/(0.7 n ln n) = 0.0261` (stated `1.03`);
  * `41.2/3 = 13.73` (stated `13.8`).

*Lemma 1.2.* Re-derived from (1.1).
* The right side of (1.1) is at least
  `mW/2 - mW/20 - 2(n-1) ln(2n+1) - 2m ln 3.884 - 2m log^+ C`.
* The margin over `3mW/8` is `3mW/40`. It suffices that
  `W > (80/3)(0.1085 + ln 3.884 + log^+ C) = 39.09 + 26.67 log^+ C`, which `W >= 27 log^+ C + 40` implies.
* The step to the whole `Q`-span is by linearity of `c -> sum c_k e^{i<k,delta>}`. Correct.

*Caveat on "answers the referee's worry".*
* With `H = e^{mW/40}`, `V_H` has relative codimension `O(1/n)` in `V(A,m)`. This is right, and it shrinks the base
  locus relevant to ZE.
* The new hypothesis SL, however, uses `V_0 = V_{H_0}` with `H_0 = (|S| n^{m-1})^2`. There
  `ln(2H_0 + 1) ~ 2m ln n`, so Lemma 1.1 only gives `codim <= (U ln 2 + 1.03 E m ln n)/(2m ln n) ~ E/2`. The worry
  is therefore answered for ZE but not for SL.
* The write-up notes this itself (open point 2(b)). If one used `H = e^{mW/40}` in SL, the hypothesis would depend
  on `W`.

## 2. Props 2.1-2.5 (angle a)

**Verdict: correct. Two wording points.**

* **Prop 2.1.**
  * (i) is the annihilator duality.
  * (ii): the remainder bound `|c|_1 (|t| D_u)^{m+1}/(m+1)!` is uniform for bounded `u`. So secant directions
    `u_i -> u_hat` of `Z_P` satisfy `Phi_c(u_i) = O(s_i) -> 0`. Correct.
  * (iii) The least-space identification is right: the lowest-order term of `sum c_k e^{<k',x>}` has degree `m`
    iff `c in V(S,m)`, and it is then `Phi_c/m!`. It is not used.
* **Prop 2.2.**
  * (i): tensor Lagrange interpolation on the value grid of `y = (<k',v_j>)`.
  * (ii): `n - 1` independent cheap directions make `k' -> y` injective, so `Poly_{<m}|_S = Q^S`.
  * Correct.
* **Prop 2.3.**
  * (ii): the multilinear representative of `(c_0 + <alpha,t>)^m` has coefficient `m! alpha^{J'}` at `t^{J'}` for
    `|J'| = m`. The reduction `t_i^2 -> t_i` does not raise degree, and the representative is unique.
  * The remark "each slab spent by at most `w` of the `v_i`" is right: the `<v_i, b_q>` are even and their sum is
    `<= 2w`.
* **Lemma 2.4.** Correct (`D_{u_hat} <= |k'|_1 <= (n-1)T = D/2`).
* **Prop 2.5.** Correct.
  * Mean value theorem with `|<k',y>| <= T(n-1)` and `|grad Phi_c| <= m|c|_1 (D/2)^{m-1} T sqrt(n-1)`.
  * `m < (n-1)(2n+1)` follows from `binom(m+n-2,n-1) >= (m/(n-1))^{n-1}`. The arithmetic gives
    `0.977/9.1 = 0.107 >= 0.1` at `n = 64`.
  * The arithmetic bound via the corrected Prop 5.3(iii) is right, and so is the threshold
    `X_v < W + 2 ln(14 n^{3/2}/(C^2|v|_2))`.
  * *Wording 1.* "Membership in `Z_P` is decided by higher-order terms" is imprecise. If `L` is a cheap span, then
    `sum_k c_k e^{it<k',u>} = 0` identically for `u in L`, so `sum_k c_k <k',u_L>^j = 0` for **every** `j`. Near `L`
    all Taylor orders are `O(eta)`, so what decides membership is the transverse behaviour, not the next order.
  * *Wording 2.* "Provably inconclusive" means that the sufficient condition (2.1) provably fails there. It is
    not shown that actual clusters do sit there ("may sit" in the verdict paragraph is the accurate phrasing).

These results are all negative or structural. None of them is used in §3.

## 3. Theorem 3.1 (unconditional necessary condition on the profile)

### 3.1 The proof is correct

* (i) Cost:
  * for `q != p`, `h_{supp c}(a_q) = (e_q/2) h(beta_q) <= e_q w`;
  * for `p`, `|<k,beta_p>| <= m/e_p` gives `h(a_p) <= m`;
  * so `cost <= wW + m log p <= mW/4 + mW/8`, using `m >= 4w` (round-10 Thm 2.1(d)).
* Height: `ln H_0(S) <= 2(n-1) ln(2n+1) + 2(m-1) ln n <= mW/40` as soon as `W >= 80 ln n + 8.68`. The
  stated `+ 9` suffices, since `80 ln(129)/(0.7*64) = 8.677`.
* Lemma 1.2 applies. `V_0(S)` is contained in `V_0(A_p)` for `S` in `A_p`, because `H_0` is monotone in `|S|`.
* (ii) I re-derived the factorisation behind round-10 Prop 5.5(i):
  * `G_k(z) = pi^{eD/2 + <k,a>} conj(pi)^{eD/2 - <k,a>} (N/p^e)^{D/2} prod_x y_x^{k_x}`, using
    `conj(y_x) = (N/p^e)/y_x`.
  * Term `k` of `F_c(z)` therefore has `pi`-valuation exactly `eD/2 + v_p(c_k) + l_p(k)`, times a `pi`-unit.
  * A unique minimiser gives `v_pi(F_c(z)) < infinity`. Rescaling by an integer keeps `UM_p`.
* Correct. Cor 3.2 (at most seven heavy primes) is correct.

### 3.2 But the widened box gives the level obstruction back

On `Z^n_0`, `<k, beta_p> = (2/e_p) l_p(k)` with `l_p(k) = <k, a_p>` an integer. So
`|<k,beta_p>| <= m/e_p` iff `|l_p(k)| <= m/2`, and `l_p` takes at most `2 floor(m/2) + 1` values on `A_p`.

**Lemma R2 (levels on the widened box; proved).** Let `S` be finite in `Z^n_0`, `c in V(S,m)` rational, and let
`l` be an integer-valued linear form whose values on `supp c` lie in `{j_0, ..., j_0 + L}`.
* (a) If `L <= m - 1`, no `c in V(S,m)` has `UM_p` with respect to `l`.
* (b) If `L = m` and `c` has `UM_p` at `k*`, then `l(k*) = j_0`: the unique minimiser is on the bottom level. By
  symmetry, at `conj(pi)` the unique maximiser of `l - v_p(c)` is on the top level.

*Proof.* (a) is round-10 Lemma 5.4, applied as in Prop 5.5(ii).

(b) Let `mu_j = sum_{l(k) = j} c_k`. Since `l^t` has degree `t < m`, the vector `(mu_j)_{j_0 <= j <= j_0+m}` is
orthogonal to all polynomials of degree `< m` in `j` on `m + 1` consecutive integers. That annihilator is
one-dimensional and spanned by the `m`-th difference, so `mu_{j_0+i} = kappa (-1)^i binom(m, i)`.
* If `kappa = 0`, the argument of Prop 5.5(ii) at the level `l(k*)` gives a contradiction.
* So `kappa != 0`. Put `i* = l(k*) - j_0` and `V = v_p(c_{k*})`.
* Every `k != k*` has `v_p(c_k) + l(k) > V + l(k*)`.
* At level `l(k*)` this gives `v_p(mu_{l(k*)}) = V`, so `v_p(kappa) = V - v_p(binom(m,i*))`.
* At level `j_0` it gives `v_p(mu_{j_0}) = v_p(kappa) > V + i*`.
* Hence `-v_p(binom(m,i*)) > i* >= 0`. This forces `i* = 0`, because otherwise the left side is `<= 0 < i*`.

QED

**Consequences for the write-up.**
1. **Odd `m`.** `l_p` takes at most `m` values on `A_p`. By (a), SL(P,p) fails for every profile and every `p`.
   This is exactly the author's own Prop 3.5(i), since `h_{A_p}(a_p) <= m - 1 < m`.
   * So "on `A_p` the level obstruction of round-10 Prop 5.5(ii) disappears" (§0 item 3, §3.1) is false whenever
     `m` is odd.
   * Theorem 3.1(ii) and Cor 3.2 are then true but empty: every profile satisfies their conclusion.
2. **Even `m`.** Only `k*` on the bottom level `l_p = -m/2` can carry `UM_p`, and the top level must also meet
   `supp c`. So one level, not all of them, is freed.
3. **How often `m` is odd** (exact, `ref_m_parity.py`, big-integer `binom`). `m(n)` is odd for 4076 of the 8137
   integers `n in [64, 8200]`. Examples:
   * `m(2137) = 3362393`, `m(3068) = 6929205`, `m(3528) = 9162363`, `m(4096) = 12349471` and
     `m(8192) = 49388999` are all odd;
   * `m(16384) = 197534058` is even.

   Among the Paley orders `M = q+1` (`q = 3 mod 4` prime) in the regime:
   * round-5 profile, `r = 9M - 8`: 130 of the 269 up to 8200 have `m` odd, the first being `M = 3528`;
   * skew double core, `r = 2(M-1)`: 144 of the 297 up to 8200 have `m` odd, the first being `M = 3068`.

   At each such `n`, the Hadamard profile is admissible, so SL(n) is false there.
4. **Exact confirmation** (`ref_bottom_level.py`, `ref_bottom_control.py`, `sl_exact.py`; my own exact
   rational-nullspace plus `p`-saturation code).
   * Random models in dimension 2-3 with `|l| <= m/2`: no `UM_p` point at odd `m`, and at even `m` all 323
     `UM_p` points lie on the bottom level (60 instances, `p = 2, 3, 5`).
   * Control: with `|l| <= m` (wider), 456 non-bottom `UM_p` points appear, so the tester is not biased.
   * 118 of 118 sampled `UM_p` points were confirmed by an explicit certificate `c`, verified from scratch
     (`Mc = 0` exactly, unique `p`-adic minimiser).

**Repair (easy).** The theorem has slack. Either replace `m` by the even number `m_e = 2 floor(m/2)`, or widen to
`|l_p| <= lambda m/2` with `lambda log p <= W/8` (for example `lambda = 2` and `log p <= W/16`, which gives `2m + 1`
levels).
* For the first option, `m_e >= m - 1 >= max(0.7 n(n-1), 4w)` in the regime (checked exactly on `[64, 8200]` for
  `r = n, 9n-8, n^2`).
* Both options leave the proof of Theorem 3.1 unchanged.
* Neither option removes the obstruction of §4, which concerns duplicate and twin columns, not parity.

### 3.3 Interpretation

Theorem 3.1 is a correct unconditional statement. As a necessary condition on profiles it is weak:
* it is empty for odd `m`;
* it is automatically met at all but at most 1-2 columns by the barrier profiles (§4.1);
* it concerns only bottom-level `k*` for even `m`.

It is not a fake-violated input: fake profiles satisfy it at every column that §4.1 covers.

## 4. Theorem 3.3 (SL(n) implies n = O(loglog R))

The deduction is correct as an implication. I checked:
* Lemma 3.1 and Lemma 3.2 give `J` with `n' >= n/4 - 1` and exponent rank `n' - 1`;
* (I4) and `W' >= 4n log^+ C'` give pair separation `>= 1/2 - 1/(2n) >= 1/2 - 1/(2n')`;
* (3.1) holds for `J`;
* Cor 3.2 leaves at most 7 SL columns, against the 8 that SL(`n'`) requires.

The problem is the hypothesis.

### 4.1 SL(n) is false: duplicate and twin columns

**Proposition R1 (twins kill SL; proved).** Let `n` satisfy `n - 1 >= 300 ln(4r)`, so that `m >= 4w`. Let `P` be a
profile with `e_q = 1` for all `q`, and let `p` be a column. Suppose some column `q != p` satisfies
`a_p - a_q = sum_{x in X} s_x e_x` with `X` a set of at most two points that does not contain the base point `n`.
(`X` empty means a duplicate.) Then SL(P,p) fails.

*Proof.* On `A_p` the slab of `q` is narrow, so `|l_q| <= w/2`. Also `l_p = l_q + sum_{x in X} s_x k_x` with
`|k_x| <= T`. Hence `l_p` takes values in an interval of length `<= w + 4T`. Now `w + 4T < 4w <= m`, because
`3w = 3n sqrt(8(n-1) ln 4r) > 4n`. Lemma R2(a) then excludes `UM_p` on `V(A_p, m)`, hence on `V_0(A_p)`. QED

**Corollary R1' (the barrier profiles violate SL(n)).** Each of the following is admissible in the sense of
theory.md §3.3 and satisfies SL(P,p) at no more than the stated number of columns. SL(n) requires 8.
* **(a) Round-5 Paley fake profile** `P(H_q, pi, 8, a*)` (`round5/construct.md` §2, Thm 2), at every
  `M = q + 1 >= 3528` with `q = 3 mod 4` prime (`r = 9M - 8`): at most **1** column.
  * Every unflipped column has 7 identical copies, so `X` is empty.
  * The flipped copy of row `x` differs from an unflipped copy of label `a(x)` exactly in row `x`.
  * With the base point at row `x_inf`, only the flipped copy of that row has no twin.
* **(b) Skew double core** (`round5/exact.md` Thm 4.1), at every skew-Hadamard order `M >= 3068`
  (`r = 2(M-1)`): at most **2** columns.
  * `S'_a = S_a xor {0, a}` is a two-row twin of `S_a`.
  * With the base point at row 1, only the pair `a = 1` can escape.
* **(c) Any Hadamard profile with every column doubled**: **0** columns, at every Hadamard order `n` with
  `n - 1 >= 300 ln(8(n-1))`.

*Admissibility.*
* *Rank `n - 1`.* For (a) and (b) this is the twin argument of round 5 ("rank `S = M`"). For (c) the Hadamard
  rows are independent.
* *Pair separation.* Round-5 Thm 1 gives `G >= 4` at pairs, i.e. at least `(r + 4)/2` separating columns for (a).
  For (b), exactly `M` of the `2(M-1)` cuts separate each pair. For (c) the fraction is `n/(2(n-1))`.
  * With weights `log q` in a window `[tau, tau + eta]`, `eta <= 1/(4r)`, as in the round-5 constructions, every
    fraction is `> 1/2 > 1/2 - 1/(2n)`.
* *Regime.* `n - 1 >= 300 ln(4r)` holds at `M = 3528` (`300 ln 126976 = 3525.5`) and at `M = 3068`
  (`300 ln 24536 = 3032.4`), and at every larger `M`.

*Exact check* (`ref_fake_sl.py`; profiles rebuilt from the round-5 definitions, Hadamard and skew identities
asserted). Instances: Paley fake at `q = 7, 11, 19, 23` (`r = 64, 100, 172, 208`) and skew double core at the same
`q` (`r = 14, ..., 46`).
* Exponent rank is `n - 1` in every case.
* The minimum pair separation is `0.514-0.571`, against the required `0.44-0.48`.
* Minimising over the base point, the number of columns without a twin is **1** (Paley fake) and **2** (skew
  double core), as stated.
* At regime sizes the twin structure is the same by construction. The inequality `w + 4T < m` was checked
  exactly for every `n in [64, 8200]` in the regime (`ref_m_parity.py`, empty failure list).

**Consequence.** Paley primes `q = 3 mod 4` are infinite in number (Dirichlet), so for every `n_1` some
`n' >= n_1` has an admissible profile failing SL(`n'`). **The hypothesis of Theorem 3.3 is never satisfied.**
Theorem 3.3 is a vacuous implication. Independently, §3.2 shows that SL(`n'`) fails at every `n'` with `m(n')`
odd at which an admissible profile exists, which is about half of all `n'` in the computed range.

**Where the write-up goes wrong.** After Prop 3.5 it argues that pair separation "keeps the weight on columns
[the obstructions] do not obviously hit", and that (ii) "is not expected to obstruct admissible profiles". That
reasoning addresses sparse columns and the span of coordinate vectors only. It does not address **duplicates and
twins**, which pair separation allows and which the barrier constructions use systematically. They need twins
for rank (round-5 "twins give rank `M`") and copies for weight. A delta-free hypothesis quantified over all
admissible profiles must hold on the barrier profiles. Checking it there should have been the first test.

**Small-scale illustration** (`ref_group_widen.py`, exact; `q = 7`, `T = 1`, `w = 2`, `m = 4`; evidence).
* Widening one unflipped column of the Paley fake leaves `l_p` with only 3 levels (its 7 copies stay narrow), and
  `V(A_p, 4) = 0`.
* Widening the whole label class (8 copies plus the flipped twins) gives 5 levels and 12 bottom-level `UM_p`
  points, each confirmed by an explicit verified certificate.
* With `lambda = 2` there are 9 levels and `UM_p` points at levels `-4, -3, -2`.

### 4.2 Is SL(n) at least the right kind of statement?

Even after both defects are repaired (§7), the remaining hypothesis SL\* is a statement about **all** admissible
profiles, including the barrier ones. It is plausible only if grouped widening beats every cheap-span obstruction
of Prop 3.5(ii) on such profiles.
* I found no general obstruction to SL\*.
* I also found no reason to believe it beyond the small illustration above. That illustration is at `n = 8`,
  where the model is not supercritical (`m = 4 < 4w = 8`).

Status: open.

## 5. Props 3.4-3.5, Remark 3.6, Prop 4.1

* **Prop 3.4: correct. Its scope is mislabelled in the summary.**
  * (i) `Lambda` is a `Q_p`-subspace intersected with `Z_p^S`, hence saturated. `UM_p` at `k*` corresponds to
    `e_{k*} in Lambda mod p`. Rational vectors are dense in `V (x) Q_p`.
  * (ii) For saturated `Lambda`, `(Lambda mod p)^perp = Lambda^perp mod p`, and
    `Lambda^perp = {(g(k) p^{-l(k)})} cap Z_p^S` with `g in Poly_{<m}`.
  * (iii) The cone condition gives `l(k) >= l(k*) + 1` off `k*`. The generating-function translation is right.
  * But (ii) and (iii) are stated for `V = V(S, m)`, i.e. for **SL-full**, not for SL, which uses `V_0`. Since
    `V_0 subset V(S,m)`, SL implies SL-full, and the converse is open point 2(b). "Equivalent forms of SL" (task
    summary, §0 item 3) should read "equivalent forms of SL-full; (i) holds for any `V`".
  * The corner form also never applies to `S = A_p` itself. For even `m`, `k*` must be on the bottom level
    (Lemma R2), which contains other points of `A_p`, and those lie outside any cone with `l(g_i) >= 1`. So the
    corner form is a sufficient criterion on subsets only.
  * My exact runs computed `R = Lambda mod p` directly (`sl_exact.py`). In the wide controls they found `UM_p`
    points that a restricted-support test misses (see §6).
* **Prop 3.5: (i)-(iii) correct.**
  * (iii): the minimal-valuation term dominates. (ii): `prod_j z_j^{<k',v_j>} = p^{l_p(k)}` with
    `z_j = varpi^{mu_j}`, and tensor Lagrange interpolation works over `K`.
  * The paragraph after it lists "duplicated columns" among the obstructions, then concludes that the
    obstructions do not affect admissible profiles generically. §4.1 shows they hit the barrier profiles at
    essentially every column.
  * The balance computation `sum_q (log q/W) s_q(n - s_q) >= (n(n-1)/4)(1 - 4 log^+ C/W)` is right, but it is
    irrelevant to twins.
* **Remark 3.6** is a heuristic reformulation. Nothing to check.
* **Prop 4.1: correct.**
  * `v_pi(w_x) = a_{p,x} + (e_p - a_{p,n}) - e_p`.
  * Inert primes and `1+i` divide `z_x`, `z_n` equally and `N` twice as often (the ramified prime `1+i` is not
    mentioned, but behaves the same way).
  * SL-res and the summary of arithmetic tools make no claims. One sentence needs care: "SL(P, p) is the statement
    that no point with valuation vector `a'_p` lies in the base locus ... in the strong sense that one element of
    the system has a monomial initial form" is a definition of SL. Prop 3.5(iii) gives only one direction (such a
    base point implies not SL).

## 6. Evidence sections

All of §5 is labelled evidence, and no proof uses it. Good. But it tests a different problem from SL as defined.

1. **The level structure is wrong in every model.** In SL, `l_p` has at most `m + 1` levels on `A_p` (exactly
   `m + 1` only for even `m`).
   * *Boolean models* (`bool_sl*.py`): `l(t) = |t|` on subsets of `{0,1}^M` with `M = 8..12` and `m = 3..5`, so up
     to `M + 1 >= 9` levels. The slab vectors are random `{0, +-1}` vectors and do not cap `|t|`.
   * *Fourier-uniform models*: the level is `|Y|`, again with far more than `m + 1` levels.
   * *Hadamard `n = 8`* (`had_sl.py`): the code and §5 both say `|<k,beta_p>| <= 2m`, i.e. `|l_p| <= m`, with
     `2m + 1 = 9` levels. That is twice the widening of §3.2. Its LP cost `2m log p` is not licensed by
     Theorem 3.1 as stated: it needs `log p <= W/16`, which is a correct variant but not the stated theorem.

   So "SL held iff a product certificate existed" is evidence about wide boxes, where Lemma R2 does not bite.
2. **Faithful rerun** (`ref_had_faithful.py`, exact rational nullspace plus `p`-saturation; every `k*`, no
   sampling; `p = 3, 5`). Sylvester `n = 8`, `T = 1`, narrow slabs `<= 2`.

   | box | `m` | `|S|` | levels | `dim V` | `UM_p` points by level |
   | --- | --- | --- | --- | --- | --- |
   | theory `A_p` (`|l_p| <= m/2`) | 3 | 57 | -1..1 | 21 | none |
   | theory `A_p` | 4 | 105 | -2..2 | 25 | `{-2: 24}` (all bottom points) |
   | theory `A_p` | 5 | 105 | -2..2 | 1 | none |
   | `had_sl` box (`|l_p| <= m`) | 3 | 137 | -3..3 | 101 | `{-3: 16, -2: 24, -1: 16, 0: 1}` |
   | `had_sl` box | 4 | 145 | -4..4 | 54 | `{-4: 4, -3: 16, -2: 24, 0: 1}` |
   | `had_sl` box | 5 | 145 | -4..4 | 10 | `{-4: 1}` |

   One explicit certificate per row with `UM_p` points was rebuilt and verified from scratch. The theorem's box
   behaves exactly as Lemma R2 predicts.
3. **A misreading in the `had_sl` paragraph.** The write-up says SL "held exactly at the `k*` on level `-2` ...
   and failed at `l* = -1, 0`". Two things are wrong with this.
   * Only 10 `k*` were sampled.
   * "Failed" there means that the *restricted-support* test failed. That test uses
     `S = {k*} cup {higher levels}`, which is only a sufficient criterion.

   In the same box (`m = 4`) the exact criterion gives `UM_p` points at levels `-4` and `-3` (not sampled) and at
   one level-0 point. The conclusion "SL is decided by product structure" is therefore drawn from a sample and a
   one-sided test.
4. **The char-0 Boolean claims** rest on ranks modulo `2^31 - 1`. A "yes" can be a modular artefact. The write-up
   says so, which is correct.
5. **Scale.** None of the models is in the supercritical regime of Theorem 2.1 (`m >= 4w` at `T = n`). For example,
   the Hadamard model has `m = 4 < 4w = 8`. The write-up says so too.

The reading I can support is: in small faithful boxes without twins, SL-full holds at the bottom level for even
`m`, and fails for odd `m`.

## 7. Repair: grouped widening (Theorem R3, Hypothesis SL\*)

The author's mechanism survives once the widening is applied to a whole group of columns and made at least two
levels wider. Notation:
* `P` is a profile with base point `n`;
* `W_P = sum_q e_q log q`, so `W_P <= W = log N` for a cluster;
* for a set `Q` of columns and integers `lambda_q >= 1` (`q in Q`), put
  `A_{Q,lambda} = {k in B_T : |<k,beta_q>| <= w (q notin Q),  |l_q(k)| <= lambda_q m/2 (q in Q)}`, with `T`, `w`,
  `m` from round-10 Thm 2.1;
* `(Q, lambda)` is *light* if `sum_{q in Q} lambda_q log q <= W_P/8`.

**Theorem R3 (proved, by the author's argument).** Let a cluster satisfy (3.1) of theory.md, and let `(Q, lambda)`
be light. Then:
* `F_c(z) = 0` for every `c in V_0(S)` and every `S` contained in `A_{Q,lambda}`;
* for no `p in Q` does `V_0(A_{Q,lambda})` contain a vector with `UM_p`.

*Proof.* For `q in Q`, `h_{supp c}(a_q) <= lambda_q m`. Hence
`cost(c) <= w sum_{q notin Q} e_q log q + m sum_{q in Q} lambda_q log q <= wW + mW/8 <= 3mW/8`. The height bound
`ln H_0 <= mW/40` and Lemma 1.2 are unchanged. The `UM_p` step is round-10 Prop 5.5(i), exactly as in
Theorem 3.1(ii). QED

`Q = {p}` with `lambda_p = 1` is theory.md's Theorem 3.1.
* `lambda_p >= 2` gives at least `2m + 1` levels, so Lemma R2 no longer constrains the position of `k*` and
  parity no longer matters.
* Taking `Q` to contain all duplicates and twins of `p` removes Proposition R1.

**Hypothesis SL\*(n).** Every profile `P` on `n` points is covered if it has
* exponent rank `n - 1`,
* `n - 1 >= 300 ln(4r)`, and
* `sum_q (log q/W_P) |a_{q,x} - a_{q,y}| >= 1/2 - 1/(2n)` for all `x != y`.

For every such profile there are a light `(Q, lambda)` and a `p in Q` such that `V_0(A_{Q,lambda})` contains a
vector with `UM_p`.

**Theorem R4 (conditional, proved).** Assume SL\*(`n'`) for all `n' >= n_1 >= 64`. Then every normalised cluster on
a squarefree norm with `W' >= 4n log^+ C'`, `W' >= 27 log^+ C' + 40` and `W' >= 80 ln n + 9` satisfies
`n < 4 max(n_1, 300 ln(4r) + 1) + 4`.

*Proof.* As for theory.md Thm 3.3. The subcluster `J` has an admissible profile `P_J`:
* pair separation holds with `W_{P_J} <= W'` in the denominator, which only helps;
* `W_{P_J}/8 <= W'/8`, so light for `P_J` implies light for the cluster.

SL\* gives a light `(Q, lambda, p)` with `UM_p`, and Theorem R3 forbids it. One witness suffices, so no count of
columns is needed. QED

**Status of SL\*.**
* Open.
* Unlike SL(n), it is not refuted by Proposition R1 or Lemma R2. Grouping a twin class costs weight
  `O(1/n)` on the barrier profiles, far below `1/8`.
* The exact `q = 7` illustration (§4.1) shows grouped widening producing verified `UM_p` vectors on the round-5
  fake profile, where single widening provably cannot. This is evidence only, and outside the supercritical
  regime.
* The obstructions that remain are those of Prop 3.5(ii)-(iii) for `A_{Q,lambda}`: `a'_p` lying in the rational
  span of directions that stay cheap after widening `Q`, or a `p`-adic base point with valuation vector `a'_p`. A
  profile whose columns are short rational combinations of other columns in many *disjoint* ways would defeat
  every light `Q`. I do not know whether such admissible profiles exist.
* **The first test for SL\* is the barrier profiles themselves** (Paley fake and skew double core, `Q` = a label or
  twin class). Theorem R4 can say nothing unless SL\* holds there.
* Note also that SL\* concerns profiles only. Like SL, it does not distinguish actual clusters from fakes. A fake's
  profile satisfying SL\* is not a contradiction, because fakes have no `F_c`. So SL\* is consistent with the
  barrier, but nothing about the barrier makes it true.

## 8. Corrected statement

**No growth improvement and no uniform bound are proved. ZE is not proved. The proposed replacement SL(n) is false
as stated.** What is proved:

1. **Lemma 1.1** (`codim(V_H in V(S,m)) <= (U ln 2 + E_S ln(UB))/ln(2H+1)`) and **Lemma 1.2** (the Thm 4.1
   dichotomy on any `S` in `B_T` with `cost <= 3mW/8`, heights `<= e^{mW/40}`, `W >= 27 log^+ C + 40`). As stated.
   The `O(1/n)`-codimension consequence applies to `H = e^{mW/40}`. At the height `H_0` used in SL it gives only
   about `E/2`.
2. **Props 2.1-2.5.** As stated, with two wording fixes: near cheap strata every Taylor order is `O(eta)`, and
   "inconclusive" refers to the sufficient condition (2.1).
3. **Theorem 3.1 and Cor 3.2.** As stated (unconditional, delta-free), with two qualifications.
   * By Lemma R2 its conclusion holds trivially for every profile when `m` is odd, which is about half of all `n`,
     including `n = 2137, 3068, 3528, 4096, 8192`.
   * For even `m` it concerns only bottom-level `k*`.
   * Repair: Theorem R3 (grouped widening, `lambda_q >= 2` allowed).
4. **Hypothesis SL(n) is false** whenever `m(n)` is odd and an admissible profile exists (Lemma R2). It is also
   false at every Paley order `M >= 3528` (round-5 fake profile) and every skew-Hadamard order `M >= 3068` (skew
   double core), regardless of parity, because duplicate and twin columns cap `l_p` (Prop R1).
   * **Theorem 3.3 is therefore vacuous.** Its proof is a valid deduction, but its hypothesis never holds.
5. **Conditional replacement (Theorem R4).** Hypothesis SL\*(n) (grouped widening, one light witness) implies
   `n = O(log W) = O(loglog R)` on squarefree norms. SL\* is open, delta-free and decidable per profile. Prop R1 and Lemma R2 do not
   refute it, but it has been tested on the barrier profiles only at `q = 7`. It must be tested on them first.
6. **Props 3.4-3.5 and 4.1.** As stated. Prop 3.4(ii)-(iii) characterise SL-full (`V(S,m)`), not SL (`V_0`). The
   heuristic after Prop 3.5 is withdrawn.
7. **Evidence (§5 of the write-up)** concerns boxes with `>= 2m + 1` levels, not `A_p`. On the faithful box the
   exact behaviour is as Lemma R2 predicts.

## 9. Issue list and scripts

| # | severity | issue |
| --- | --- | --- |
| 1 | fatal (for Thm 3.3) | SL(n) is false: Prop R1. Duplicates and twins (round-5 Paley fake, skew double core, doubled Hadamard) are admissible and fail SL at all but `<= 2` columns at every Paley order in the regime. Thm 3.3 is vacuous. |
| 2 | major | Parity: for odd `m`, `|l_p| <= m/2` gives `<= m` levels, so SL fails for every profile (Lemma R2(a), the author's own Prop 3.5(i)). For even `m`, only the bottom level is possible (Lemma R2(b)). "The level obstruction disappears" is false. `m(n)` is odd for about half of all `n`. |
| 3 | major | Evidence boxes have `>= 2m + 1` (`had_sl`) or `M + 1 >= 9` (Boolean, Fourier) levels and do not model `A_p`. `had_sl`'s widening is not licensed by Thm 3.1 as stated. The "failed at `l* = -1, 0`" comes from a sufficient-only sampled test. |
| 4 | minor | Prop 3.4(ii)-(iii) are forms of SL-full, not SL. The corner form never applies to `A_p` itself. |
| 5 | minor | Lemma 1.1's codimension gain is for `H = e^{mW/40}`, not `H_0`. The "referee's worry" is answered for ZE, not for SL. |
| 6 | minor | Prop 2.5 wording ("higher-order terms"; "provably inconclusive"). |
| 7 | minor | The heuristic after Prop 3.5 ignores twins. The balance argument is irrelevant to them. |

Scripts (`round11/referee_theory_checks/`; all mine, nothing imported from `theory_checks/`):

| script | what | status |
| --- | --- | --- |
| `ref_m_parity.py` | exact `m(n)` (big-integer `binom`) on `[64, 8200]` plus `n = 8192, 16384`; parity; Paley orders of both barrier families in the regime; repair inequalities `m_e >= max(0.7n(n-1), 4w)` and `w + 4T < m` | exact |
| `sl_exact.py` | exact SL-full tester: rational nullspace, `p`-saturation, `R = Lambda mod p`; explicit certificate rebuilt and verified (`Mc = 0`, unique `p`-adic minimiser) | exact library |
| `ref_bottom_level.py` | Lemma R2 on 60 random models with `|l| <= m/2`: 0 `UM_p` points at odd `m`, 323 at even `m`, all on the bottom level | exact |
| `ref_bottom_control.py` | control with `|l| <= m`: 456 non-bottom `UM_p` points; 118/118 certificates verified | exact |
| `ref_had_faithful.py` | Hadamard `n = 8` on the theory box versus the `had_sl` box, `m = 3, 4, 5`, `p = 3, 5`, all `k*` | exact |
| `ref_fake_sl.py` | rebuilds the round-5 Paley fake and skew double core (`q = 7..23`): rank, pair separation, twin count per base point | exact |
| `ref_group_widen.py` | Paley fake `q = 7`: single versus grouped widening | exact; evidence |

Each script writes the `*_out.txt` of the same name (`python3 <script> > <script>_out.txt`; `ref_m_parity.py` takes
the arguments `64 8200`).
