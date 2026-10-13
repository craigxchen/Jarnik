# Referee report on round10/L4K.md

Adversarial referee, round 10. Scripts are in `referee_L4K_checks/`; working notes are in
`referee_L4K_checks/notes.md`. All code is independent of the author's `L4K_checks/` (no imports).
Labels: *proved*, *exact* (finite computation, script named), *evidence*, *heuristic*.

## 0. Verdict

**The refutation survives. The corrected question does not survive as displayed.**

* **Theorem 1 (forced/free boundary): correct.**
  * I re-derived every inequality.
  * Independent exact check: `ALL OK`.
  * The forced/free thresholds are reproduced exactly for every first and last prime set (`K = 2, 3`,
    `r = 4..8`; `K = 4`, `r = 4`).
* **Theorem 2 (L4(K) false for every `K >= 2`): correct, after minor repairs.**
  * Repairs: (2.1) is undefined at `q | 2N`; the non-unit convention; an unsupported but unnecessary
    "`r = Theta(M)`"; the range of `Y`.
  * The union bound and the asymptotics re-derive. The hypothesis can be weakened to
    `r >= C(K, h, delta) log M`.
  * The systems are not fakes. Their norm-`(2h+2)` monomials collide at **every level** at a positive proportion
    of the free moduli (verified).
* **Lemma P: correct.** It is a level-1 statement.
* **All exact data that I re-ran reproduce exactly:** the `5,066,549,580,791,808` / `2^37` count, the
  signature/class counts, the `K = 4` one-modulus deformations at six moduli, the `K = 4` joint core
  (1,966,080, all genuine) and the `K = 4` flexible top. The instance coherence shares are consistent with mine.
* **Major M1.** The threshold "`K >= 2L(h_0 + h)` (twice the README)", labelled proved, is wrong as a necessity
  claim, and it is aimed at the wrong theorem.
  * The Hamming conclusion of L4 feeds round-9 **Theorem B** (Heath-Brown), which needs only
    `Y >= p_j Norm zeta`. So `K ~ 2(h_0 + h)` suffices.
  * The rigid-core form concludes the genuine class itself and needs neither Theorem A nor Theorem B.
* **Major M2.** The displayed corrected question L4^forced(K) is not right.
  * It uses `F_K = {n_q <= A_K}`, which is mostly free moduli. My instances give coherent share
    `<= 0.13-0.25` there.
  * It keeps only norm-`<= K` budgets. Adding the **norm-5** budgets (which actual circles satisfy) shrinks
    the `K = 4` "flexible top" drastically (§6.3):
    * one modulus at a time, 4 of 6 tested top moduli become rigid, and the other two keep 40-96 non-genuine
      survivors against `10^5`;
    * jointly, no randomized search finds a non-genuine top.

    So the reported failure of the "most of the mass" form at `K = 4` concerns the norm-truncated statement
    only; with norm 5 the test is undecided.
  * My corrected question is L4\*(K) (§9): the exact forced set, budgets of norms `> K` included.
* **Fake test.** The write-up proposes no mechanism. I prove (§7) that **every** fake circle, given the genuine
  labels of its own primes, satisfies all hypotheses and conclusions of L4(K), L4^forced and L4\*. So the
  surviving lemma uses **no fake-violated input**: its proof would not by itself beat the barrier, and route 3's
  fake-separating power must come from the label-angle step (L7 / slabs).

No growth improvement and no new cap; the barrier map is unchanged.

## 1. What is claimed

The write-up (`round10/L4K.md`) makes four kinds of claims.
* **Proved.** Theorem 1 (forced/free boundary), Theorem 2 (refutation of L4(K) as posed in
  `round9/README.md` §6) and Lemma P (a class is determined by a short-generated collision lattice).
* **Proved by assertion.** §5 item 2: an inverse step feeding Thm A "needs `K >= 2L(h_0+h)`".
* **Exact.** Forced/free thresholds, a count of all good level-1 systems, finite Theorem-2-style instances,
  and the `K = 4` rigidity data.
* **Corrected question.** L4^forced(K) (§4.1), plus a rigid-core form (§6.1).

It claims no growth improvement and no cap. Its "fake-violated input" is "(none stated)": it proposes no
growth mechanism.

## 2. The definitions re-derived

* **beta-good.** `round9/L4.md` §1 calls a system *beta-good up to norm K* if
  `log m_(v,s) <= beta k log M` for every `v` outside the character lattice with `1 <= ||v||_1 = k <= K` and
  every `s`. The write-up's "good up to norm `K`" means per-unit demand `<= B(v,s) = w(v)/2 + [s odd](log 2)/2`
  for **every** `v in V_K`. In regime (1.4), `w(v)/2 <= (k h_0/2) log M`. So "good" implies beta-good for every
  `beta >= h_0/2 + (log 2)/(2 log M)`, which is the genuine `beta`. A system that is good and not coherent is
  therefore a counterexample to L4(K) as posed. **Direction of the refutation: correct.**
* **Level-1 collisions.** "`v` collides at `q` with some unit iff `<v, c> = 0 mod n_q`, `c = ell mod n_q`": correct,
  because `iota` is `m/4` or `3m/4`.
* **Higher levels.** The map `T_(q^a) -> T_(q^(a-1))` sends `log i` to `log i` or `3 log i`. So "a collision at
  `q^a` implies one at `q^(a-1)`" holds, up to the harmless swap of `s = 1, 3`, whose budgets are equal.
* **Genuine collisions** (exact, independent code `referee_L4K_checks/rlib.py`, `ref_genuine.py`). The genuine
  labels collide at `q^a` with unit `s` iff `q^a | D_s(Z_v)`. There are 0 mismatches in 74,080 checks: three
  prime sets, `q <= 109`, levels 1-2, norms `<= 3`, all units. Genuine parity equals `nu_j(q)` at every tested
  `q`, with 0 mismatches.
* **Coherence.** Level-1 coherence is `c(q) = u_q L^zeta(q) mod n_q`. This is exactly coherence mod `m_q` up to
  `(m_q/4) w`, so the write-up's reduction is right. Lemma P and the "genuine class" of §4 are level-1 notions.
  "Coherence with `pi` at `q` is the same as having the genuine all-unit collision set" (§4.1) should say
  *level-1* coherence (minor m8).

## 3. Theorem 1 (boundary): re-derivation and checks

**(a) Forced side.** The pigeonhole is correct. The `|S_h| > n_q` values `<v, c>`, `v in S_h`, collide, and
their difference has norm `<= 2h <= K`; the collision is level-1 with some unit. For odd `K = 2h + 1` the set
`S_h ∪ (e_1 + S_h)` has all pairwise differences of norm `<= K`, so the forced bound improves by a factor of
about 2 (constant only; minor m6).

**(b) Free side.** I re-derived every inequality.
* `2 Lambda = 2K(K+1)B + 2K > 2(K+1)KB + K`, which bounds `|2(K+1)<v, beta^sigma> + <v, nu>|`.
* If `sigma_v = 0` and `<v, beta^sigma> != 0`, then `|2(K+1)<v, beta^sigma>| >= 2(K+1) > K >= |<v, nu>|`.
* If `sigma_v = 0` and `<v, beta^sigma> = 0`, then `v^+` and `v^-` have disjoint supports. Padding both to
  size `h` contradicts `B_h`.
* `max_j a_j = 2(K+1)^2 B + 2K + 1`, so `K max_j a_j = A_K`.

Hence `<v, a>` is a non-zero integer with `|<v, a>| <= A_K < n_q`, and there is no collision at any level.
* **Correct.**
* **Exact check (`ref_thm1.py`).** Independent code with greedy `B_h` sets: `K = 1..6`, `r = 2..6`
  (`r <= 4` for `K >= 5`), every parity pattern, 6 random permutations each. Zero vanishing values; the maximum
  of `|<v, a>|` equals `A_K` in every case. `ALL OK`.
* **Bose-Chowla, as cited from memory:** `P` integers mod `P^h - 1` with distinct `h`-fold sums. This is the
  standard statement. A shift by 1 makes the elements positive, with `B <= P^h <= (2r)^h`. It affects
  constants only.

**(c) Exact thresholds (`ref_free.py`, `ref_free_K3_out.txt`, `ref_free_K4.py`).** The search is an independent
DFS over parity-respecting `c in (Z/n)^r`, with `c_1` reduced to divisors of `n` (unit orbits). Every returned
vector was re-verified. I checked every prime `q <= 4A_K + 5`, above which Theorem 1(b) makes `q` free.
* **`K = 2`** (`r = 4..8`), first and last sets: `q* = 41/31, 41/47, 47/71, 73/79, 79/89`.
* **`K = 3`** (`r = 4..8`), first and last sets: `q* = 103/137, 151/137, 151/137, 167/167, 173/199`.
* **`K = 4`** (`r = 4`, `q <= 700`): the forced sets agree exactly with `free_K4_out.txt`, 40 moduli each, with
  `q* = 199/271`.

**All match the write-up's table.** The non-initial-segment example also checks: 37 is free and 41 is forced for
`{5, 13, 17, 29}`. The "random" prime sets could not be reproduced without the author's seed and were not
checked.

**Verdict on Theorem 1: correct as stated.** It is a genuine sharpening of README Prop. H, from `M^(K+o(1))`
to `M^(floor(K/2)+o(1))`. The criticism of the README §6 sentence ("pigeonhole forces about
`#V_K/q` collisions") is also correct: by (b), zero collisions are possible once `n_q > A_K`.

## 4. Theorem 2 (refutation): re-derivation and checks

**(i) Goodness.** Correct.
* At `q > q_T = 4A_K + 5` we have `n_q >= (q-1)/4 > A_K`, so Theorem 1(b) applies for every `sigma_q`.
* At `q <= q_T` the labels are genuine, so `m_(v,s)` divides the `q_T`-smooth part of `D_s(Z_v)`.
* `D_s(Z_v) != 0` for `v != 0`, by unique factorisation over distinct split primes. So
  `m_(v,s) <= |D_s(Z_v)| <= nu_s^(-1) e^(w(v)/2)`.

**Gap (minor m1).**
* The construction is stated "at every prime `q > q_T`". At `q = p_j > q_T` the vector (2.1) is undefined,
  since `nu_j(p_j)` is a Legendre symbol of 0. Moduli dividing `N` carry level data, not labels. The same goes
  for `q = 2`.
* This happens whenever some `p_j > q_T`. That is possible in regime (1.4) when `K/2 < h_0`, e.g. for the small
  `K` of the finite tests.
* Fix: "every prime `q > q_T` not dividing `2N`", with genuine data at `q | 2N`. Then (i) still gives
  `m_(v,s) | D_s(Z_v)`. `Q_7` never meets `2N`, so nothing else changes.

**(ii) The union bound.** Correct; I re-derived every step.
* Coherence at the top level implies level-1 coherence `c^sigma = u L^zeta mod n_q`.
* The `r` integers `2(K+1) beta_i` are distinct residues mod `n_q`, since `2(K+1)B < A_K < n_q`. So each unit
  `u` admits at most one `sigma`, and `P <= phi(n_q)/r! <= (Y+1)/(4 r!)`.
* The events are independent over `q` for fixed `zeta`.
* There are `<= 2^(pi_7(Y))` sets `S`, and `<= pi(sqrt H + 1)^2 <= 13H` Gaussian integers per column.

**Convention (minor m3).** "If some `zeta_j` is not a unit mod `q`, coherence at `q` fails by definition."
Theorems A and B of round 9 instead *exclude* `q | Norm zeta` from the requirement. Under that convention a
coherent set may also contain the moduli dividing `prod_j Norm zeta_j`. Their mass is
`<= sum_j log Norm zeta_j <= r log H = O(r log M) = o(theta_7([1, Y]))`, so (iii) survives with `delta`
replaced by `delta + o(1)`. The write-up should say so; its finite instances already count non-unit moduli as
coherent.

**(iii) Asymptotics.** Correct. I re-derived:
* `k_0 >= (delta/4)(1 - o(1)) Y/log Y`;
* `log(4r!/(Y+1)) >= (1/2) r log r` once `log Y = o(r log r)`;
* the second term is `<= -c delta r^2 log r/log(16 r)`, against `r(h log M + log 13)`;
* `theta_7([1, q_T])/theta_7([1, Y]) = O_K(r^(-K/2))`.

(My constant is `delta` where the write-up has `delta/4`; theirs is weaker but valid.) Two remarks.
* **Hypothesis on `r` (minor m5).** The proof only uses `r -> infinity` with `r >= C(K, h, delta) log M` and
  `log Y = O(log M)`. "`r >= cM`" is more than needed. The parenthetical "near-extremal clusters have
  `r = Theta(M)` columns (`round8/stability.md`)" is not supported by that file: `stability.md` mentions
  `r = Theta(M)` only in its heuristic §8.4. `p3lower.md` Lemma D.0 gives only `r = O(M)`. The refutation does
  not need it, because L4(K) quantifies over all `r <= rho M`; but the sentence should be removed or sourced.
* **Range of `Y` (minor).** The abstract says "below `Q_K` (or below any `Y <= M^(O(1))`)". The proof needs
  `Y/q_T -> infinity`, e.g. `Y >= Q_K` as in (iii). For `Y <= q_T` the system is genuine and hence coherent. The
  abstract's "delta-fraction" should read "2 delta-fraction" (harmless).

**Verdict on Theorem 2: correct, after the minor repairs m1, m3, m5.** L4(K) as posed in `round9/README.md` §6
is false for every `K >= 2`. The counterexample is parity-respecting, good for every norm `<= K` at every level
and every `X`, and system-level non-coherent with every `zeta` of norm `<= M^h` on all but a `2 delta + o(1)`
fraction of `theta_7([1, Y])`. This holds whenever `r -> infinity` with `r >= C log M`, `Y >= Q_K` and
`log Y = O(log M)`. For `K = 1`, every `q` with `n_q >= 3` is free: take `c_j` odd, or `c_j = 2`. Theorem 2
does not apply there, since `q_T > Q_1`. The write-up's "trivially" needs a different free family, e.g.
uniformly random collision-free vectors, with the same union bound. I did not write that out.

**Remark 2 checked (`ref_remark2.py`).** "The system is not in P3-all\*" is right, and the wording understates
it.
* For `K = 2` the vector (2.1) has a norm-4 monomial with `<v, a> = 0` **as an integer** at 63% (`r = 4`) and
  100% (`r = 8`) of the primes in `(q_T, 20000]`.
* Such a `v` collides at **every level** there, so its demand grows like `log X`. "About as often as for random
  labels" is wrong: it is far more often.
* The general proof is short. A `B_h` set is not `B_(h+1)` once `r^(h+1)/(h+1)! > (h+1)B`. So some `v` of norm
  `2h + 2` has `sigma_v = 0` and `<v, beta^sigma> = 0`, and `<v, nu(q)> = 0` holds at a positive proportion
  of `q`.

**Finite instances (`ref_instance.py`, `ref_instance_out.txt`).** These are an independent re-implementation of
the §3 instances: genuine labels at `q <= q*`, and a random collision-free vector at every `Q_7` prime in
`(q*, Y]`. The coherent mass is the exact maximum over all pairs `(zeta_a, zeta_b)` of norm `<= H` (1573 or 6284
per column, non-units counted coherent), for column pairs `(1,2)` and `(3,4)`.

| `K` | `r` | set | `Y` | `H` | genuine share | max pair-coherent share of the free part (mine) | write-up |
|---|---|---|---|---|---|---|---|
| 2 | 4 | first | 1297 | 2000 | 0.0274 | 0.083 / 0.082 | 0.076-0.083 |
| 2 | 4 | last | 1297 | 2000 | 0.0274 | 0.078 / 0.078 | 0.076-0.083 |
| 2 | 4 | first | 1297 | 8000 | 0.0274 | 0.098 / 0.100 | 0.086-0.102 |
| 2 | 8 | first | 4625 | 2000 | 0.0187 | 0.034 / 0.026 | 0.025-0.027 |
| 3 | 4 | first | 1625 | 2000 | 0.0601 | 0.065 / 0.058 | (Y = 4700 there) |

These are consistent: the instances are random draws, so only ranges can be compared.

**Labelling (minor m2).** The summary calls these "the Theorem-2 instances", but they are a different
construction. They use the exact boundary `q*`, not `q_T`, and random collision-free vectors, not (2.1). With
Theorem 2's own construction at `r = 4`, `q_T = 621` and the genuine part would carry 48% of
`theta_7([1, 1297])`. With `r! = 24` the union bound (3.1) is empty. The body of §3 describes the construction
correctly; the summary should too.

## 5. Numerical claims re-computed

All of the following use independent code (`referee_L4K_checks/`). The only shared inputs are the definitions.

1. **Exhaustive count (§4.2), `ref_forced_sig.py`, `ref_forced_sig_out.txt`.** Reproduced exactly, in 3 s:
   * per-modulus counts of parity vectors, signatures and genuine-signature vectors: 16/16/1, 256/128/2,
     1296/336/4, 10^4/568/32, 20736/640/128, 65536/616/1024, 160000/990/4096;
   * **5,066,549,580,791,808** good level-1 systems in **36,864** signature tuples;
   * **2^37** with the genuine signature everywhere;
   * non-genuine shares 0.938/0.984/0.984/0.938/0.875/0.750/0.500.

   The genuine signature computed from discrete logs equals the `D_s` signature at every `q`.
2. **Genuine signature vs genuine class (`ref_sigclass.py`).** Reproduced: 16/32, 64/128, 256/1024 and 512/4096
   at `q = 19, 23, 31, 41`, and all of them at `q <= 11`.
   **Observation (minor m4).** At `q = 3` (`n = 1`), `q = 5` (`n = 1`) and `q = 7` (`n = 2`), the genuine class
   is *every* parity-respecting vector: class sizes 16 and 256 equal the parity counts. So "rigidity" at those
   moduli is vacuous. Consequences:
   * "`K = 2` is rigid only at `q_0 in {3, 7}`" (and the loose set's `{3, 5, 7}`) means that **`K = 2` has no
     non-vacuous rigidity at all**.
   * The `K = 3` joint-search statement "outside the genuine class at every forced modulus except `{3, 7}`"
     means *at every forced modulus where a non-genuine label exists*.
   * The write-up should say this; it strengthens its own conclusion.
3. **`K = 4` one-modulus deformations (`ref_deform.py`, `ref_deform_out.txt`).** Level-1 relaxation at `q_0`,
   genuine labels (all levels) at the other 39 forced moduli, zero demand above 199. Survivors / survivors
   outside the genuine class, mine and the write-up's:

   | `q_0` | 47 | 83 | 97 | 101 | 103 | 107 |
   |---|---|---|---|---|---|---|
   | mine | 8 / 0 | 36 / 24 | 16 / 0 | 20 / 0 | 48 / 0 | 216 / 144 |
   | write-up | 8 / 0 | 36 / 24 | 16 / 0 | 20 / 0 | 48 / 0 | 216 / 144 |

   Exact agreement.
4. **`K = 4` joint core (`ref_core.py`, `ref_core_out.txt`).** Core `{3, ..., 47}`, the other forced moduli
   genuine at all levels, core levels `>= 2` not charged. The search finds **1** good signature tuple and
   **1,966,080** good level-1 systems, **none** outside the genuine class at any core modulus. This matches
   `core_4_4_first_rest.txt`. (The non-genuine-top variant was not re-run.)
5. **`K = 4` flexible top (`ref_top.py`, `ref_top_out.txt`).** The system is genuine at all levels at forced
   `q <= 103` and collision-free above 199. A random-sampling search, independent of `forced_prefix.py`, found
   on its first seed a vector outside the genuine class at **every one of the 18** forced moduli in
   `[107, 199]`, with:
   * level-1 budgets met (0 violations);
   * level-2 lifts avoiding every level-2 collision at every top modulus, so the system is good at all levels;
   * **no** pair `(zeta_a, zeta_b)` of norm `<= 2000` coherent on all five top `Q_7` moduli
     `{127, 151, 167, 191, 199}`, for any of four column pairs.

   The top `Q_7` moduli carry `0.499` of the forced `Q_7` mass. Claim reproduced. Note, however, that the
   maximal pair-coherent share of the top `Q_7` mass is 0.41-0.60, depending on the column pair. So the
   system-level coherent share of the *whole* forced `Q_7` mass is somewhere between 0.5 (the core alone) and
   about `0.5 + 0.5 * 0.41 = 0.71` (pair upper bound). "Fails in its most-of-the-mass form" is right only if
   "all but a small part" means a fraction below about 0.29. The Theorem-B threshold `(1 + delta)/2` (§6) is
   not violated by this instance (minor m10).
6. **Thresholds (§2 table).** Reproduced for the first and last sets at every `(K, r)`; see §3(c).

## 6. Lemma P, forced range, L4^forced(K)

**6.1 Lemma P: correct.**
* `Lambda_c ⊂ Lambda_(c')`, so `v -> <v, c'>` factors through `Z^r/Lambda_c ≅ Z/n`, which is cyclic because
  `gcd(c, n) = 1`. Hence `c' = u c`, and `u` is a unit iff `gcd(c', n) = 1`.
* Two remarks:
  * Equality of the *full* collision lattices already gives `c' = u c` without the short-generation
    hypothesis. Short generation is needed only to pass from the `V_K`-collision sets to the lattices.
  * Lemma P is a level-1 statement (m8).
* §5 item 3's identification `Lambda_(L^zeta(q)) = {v : q | Im((Z^zeta_v)^4)}` is correct, since
  `Z^4 - conj Z^4 = 2i Im Z^4`.

**6.2 The displayed L4^forced(K) is not the right corrected statement (major M2).**
1. *Definition of `F_K`.* §4.1 defines `F_K = {q : n_q <= A_K}`, then adds "in practice the exact forced set". By
   Theorem 1, `F_K` contains every free modulus with `n* < n_q <= A_K`. Here `n*` is the exact boundary,
   `~ 2.5r` at `K = 2`, against `A_2 = 36r + 10`. Those free moduli carry almost all of `theta_7(F_K)`.
   * The mechanism of Theorem 2 applies to them (as evidence; Theorem 2's proof needs `n_q > A_K`). My exact
     instances (`ref_instance_out.txt`, last line of each block) bound the coherent share of `theta_7(F_K)` by:
     * forced share plus maximal pair-coherent free share `= 0.058 + 0.15-0.17` (`K = 2`, `r = 4`, `H = 2000`);
     * `0.058 + 0.19` (`H = 8000`);
     * `0.073 + 0.08-0.09` (`K = 2`, `r = 8`);
     * `0.060 + 0.06` (`K = 3`, `r = 4`).
   * So the displayed statement fails at the first test.
   * The statement must use the **exact** forced set `F_K^ex` (every parity-respecting vector collides). That set
     depends on the parity patterns, and Theorem 1 only brackets it: `|S_h| <= n* <= A_K`.
2. *The hypothesis is weaker than what is available.* L4^forced(K) keeps only the norm-`<= K` budgets. Actual
   circles, and every member of P3-all\*, satisfy the budgets of **every** norm. Every modulus of `F_K^ex` is a fortiori forced for
   every `K' > K`, and there many more monomials collide. The write-up's own pattern "the one-modulus rigid range at norm `K` is
   about the forced range at norm `K - 1`" says the higher norms are exactly what fixes the top of the forced
   range. My test of that pattern is below (6.3). The "flexible top" of §4.4 is therefore evidence against the
   truncated statement only, not against an inverse step that uses norms `<= K + 1` (or all norms) on `F_K^ex`.
3. *The threshold `K >= 2L(h_0 + h)` (major M1).* §0 item 3, §4.1, §5 item 2 (labelled *proved*) and §6.1
   assert that an inverse step feeding Thm A "must certify coherence on `Q_7` up to `C_L(8 p_j Norm zeta)^L`" and
   hence "needs `floor(K/2) >= L(h_0 + h)(1 + o(1))`, twice the README estimate". Three corrections.
   * **"Must" holds only as a statement about Thm A's hypothesis.** `C_L k^L` is Linnik's *upper bound* for the
     least prime in a class. A proof that goes through Thm A *as stated* needs that range. Nothing shows that
     route 3 needs it.
   * **The conclusion of L4(K) and L4^forced(K) is a Hamming statement** ("all but a small part of the mass").
     Such statements feed **Theorem B** of round 9 (Heath-Brown's quadratic large sieve), not Theorem A.
     * Theorem B needs only `D = p_j Norm zeta <= X`, with `delta X >= 30 log D`.
     * For a column agreeing with a sub-genuine `zeta` on `(1 + delta)/2` of `theta_7([1, X])`, the number of
       exceptional columns is `<= C_eps X^eps log X log D/delta^2`, which is `o(M)`.
     * So coherence on most of the `Q_7` mass up to `Y = M^((h_0 + h)(1 + eps))` already gives genuine height
       for all but `o(M)` columns.
     * Here "genuine height" is in the round-9 sense, `v_(p_j)(Norm zeta_j)` odd (`round9/README.md` §1 item 3).
       The stronger Thm-A conclusion ("`p_j Norm zeta_j` a square or twice a square") does need the Linnik range.
     * The analytic requirement is therefore **`floor(K/2) >= (h_0 + h)(1 + o(1))`, i.e. `K ~ 2(h_0 + h)`**. This
       is *smaller* than the README's `L(h_0 + h)` for every admissible `L >= 2`, not twice it.
   * **The rigid-core form (§6.1 of the write-up) needs neither Thm A nor Thm B.** It concludes the *genuine*
     class, i.e. coherence with `zeta = pi` itself. Its required length "`q_core(K) >= C_L(8 M^(h_0+h))^L`" is
     therefore not the relevant target. The length that route 3 needs from a genuine core is whatever L7 (the
     label-angle step) needs, and nobody has quantified it.

**6.3 Does the `K = 4` "flexible top" survive the norm-5 budgets? (exact per test; evidence overall).**
Setting: `ps = {5, 13, 17, 29}` and the `K = 5` forced set, independently computed (`ref_free_K5_out.txt`: 65
moduli, `q* = 349`, identical to the author's list).

1. *One-modulus deformations with all norms `<= 5`* (`ref_deform.py 5`, `ref_deform_K5_out.txt` for 107-151; `ref_deform2.py 5`, `ref_deform2_K5_*.txt`
   for 167, 199, 211; same code, lower memory). The labels
   are genuine at every other `K = 5`-forced modulus, collision-free above 349; the relaxation is level 1 at
   `q_0`. Survivors / outside the genuine class:

   | `q_0` | 107 | 113 | 127 | 151 | 167 | 199 |
   |---|---|---|---|---|---|---|
   | norm `<= 4`, reference genuine on `F_4` (write-up, reproduced at 107) | 216 / 144 | 336 / 240 | 6592 / 6080 | 8208 / 8064 | 161088 / 160320 | 500800 / 500160 |
   | norm `<= 4`, reference genuine on `F_5` (`ref_deform_mix.py`) | 36 / **0** | -- | 1152 / 640 | 1728 / 1584 | 90816 / 90048 | -- |
   | norm `<= 5`, reference genuine on `F_5` | 18 / **0** | 24 / **0** | 32 / **0** | 36 / **0** | 120 / 96 | 120 / 40 |

   Effects of the extra norm:
   * At 127 and 151 the extra genuine moduli alone leave the modulus flexible (640 and 1584 non-genuine
     survivors), and adding the norm-5 monomials makes it rigid.
   * At 107 the extra genuine demand on `(199, 349]` already suffices.
   * At 167 and 199 the norm-5 monomials cut the non-genuine survivors by factors of about `10^3`
     (90048 to 96) and `10^4` (500160 to 40), but **not to 0**. These two `Q_7` moduli stay one-modulus flexible
     at norm 5.
2. *Joint searches with norm `<= 5`* (`ref_top5.py`, `ref_top5b.py`, `ref_topQ7.py`; random sampling of `10^7`
   vectors per modulus, greedy, level-2 lifts checked).
   * Non-genuine at all 19 forced moduli of `[107, 199]`, with genuine labels at `q <= 103` and on `(199, 349]`:
     every run (3 seeds, plus 2 seeds of the min-collision variant) finds one non-genuine vector at 107
     (13 collisions) and **none** at 109.
   * Non-genuine only at the five top `Q_7` moduli, genuine elsewhere: a non-genuine vector at 199
     (10 collisions), then **none** at 191 (2 seeds).
   * The same norm-`<= 4` search succeeds on its first seed at all 18 moduli (§5 item 5). It still succeeds with
     genuine labels at every prime of `(199, 349]` or `(199, 1000]` (`ref_top2.py`, 2 seeds each).
3. *Reading (evidence).* The norm-4 flexible top is robust to genuine labels above the forced range. The
   norm-5 budgets, which every P3-all\* member and every actual circle also satisfy, shrink it drastically.
   * One modulus at a time: four of the six tested moduli become rigid, and the other two keep 40-96
     non-genuine survivors, against `10^5`.
   * Jointly: no randomized search found a non-genuine system even on two top moduli.

   So the write-up's §4.4 conclusion ("L4^forced(4) fails at this test in its most-of-the-mass form") is a
   statement about the norm-truncated hypothesis. With norm 5 added, the test is undecided: the one-modulus
   flexible moduli 167 and 199 (two of the five top `Q_7` moduli) might still combine, but my searches did not
   find it. At `q_0 = 211`, the first norm-5-forced modulus above the norm-4 range, there are 884 survivors, 832 of them
   non-genuine (`ref_deform2_K5_211.txt`). This is consistent with a flexible top that moves up with the norm, as
   the write-up's pattern `q_core(K) ~ q*(K - 1)` predicts.

## 7. Fake-circle test

The write-up proposes no growth mechanism ("fake-violated input: none stated"), so the mandatory test applies in
two ways. Are the counterexamples fakes? And would the surviving lemma, if proved, use an input that the fakes
violate?

1. **The Theorem-2 systems are not fakes (verified).** `ref_remark2.py` shows that they violate norm-`(2h+2)`
   budgets at **every level** at a positive proportion of the free moduli: for `K = 2`, 63% (`r = 4`) and 100%
   (`r = 8`) of `q in (q_T, 20000]`. The demand of the colliding monomial grows like `log X`, so for large `X`
   they fail the per-unit height condition that every member of P3-all\* satisfies. Actual circles satisfy it
   (genuine collisions are `q^a | D_s(Z_v)`, `ref_genuine.py`). **The input they falsify** is goodness of
   monomials of norm `> K`, all levels, at moduli above `M^(floor(K/2) + o(1))`. This is the write-up's
   statement, confirmed and sharpened.
2. **L4^forced(K) is not a fake-violated input, and cannot be (proved here).** Every fake circle has genuine
   primes (round 4/5 definition: angle models have genuine primes, a genuine exponent profile and real angles).
   This includes the skew double core of `round5/exact.md` and the angle models of `round5/construct.md`. Equip
   such a fake with the **genuine labels** `L_j = log(pi_j/conj pi_j)` of its own primes. By the identity
   "collision at `q^a` with unit `s` iff `q^a | D_s(Z_v)`" and `D_s(Z_v) != 0`, the system is:
   * parity-respecting;
   * good for **every** norm at every level and every `X` (`m_(v,s) <= |D_s(Z_v)| <= nu_s^(-1) e^(w(v)/2)`);
   * system-level coherent with `zeta = pi` (norm `p_j`) at every `q`.

   So the label hypotheses *and* the conclusion of L4(K), L4^forced(K), L4\*(K) and the rigid-core form all hold
   for every fake circle whose primes meet the regime hypotheses (`p_j <= M^(h_0)`, `r <= rho M`). This needs
   `h >= h_0`, which L4 needs anyway: by Thm A every coherent `zeta_j` has `Norm zeta_j >= p_j`, so for
   `h < h_0` the conclusion would fail even for actual circles. Prime sizes play no other role. Consequences:
   * no argument whose only inputs are the label conditions of L4^forced (parity plus per-unit budgets on the
     forced range) can separate a fake from an actual circle;
   * in route 3 (L4 + Thm A/B + L7), the fake-violating input must come from the step that links the labels to
     the **angles** (the slab conditions of P3-all\*, i.e. P3-G/L7), with residue moduli up to `M^A`, `A`
     chosen after `c_0` (round 9 R1).

   **The write-up's §5 item 5 ("nothing here changes the barrier map") is correct but should say this
   explicitly.** L4^forced is a structural lemma about labels, not a barrier-beating input.
3. **Round-5 fakes with residue data.**
   * The `construct.md` angle models carry **point-level** residue data. Their `M = 12` instance violates the
     prime-level parity law (`round5/README.md`, referee of `construct.md`), so they are not parity-respecting
     label systems and lie outside L4's hypotheses altogether.
   * The P3 fakes of `sharp.md` / `lower.md` have every `p_j > X` (round 9 R1).
   * At `p_j > X` the per-monomial budgets `w(v)/2 >= (||v||_1/2) log X` exceed the average forced load
     `~ (K/2) log r` by a factor `A`. There, non-coherent good systems on the forced range should be easy
     (heuristic, not checked). This is consistent with the write-up's load accounting (§5 item 4): rigidity can
     only come from tight budgets, i.e. small primes, which is regime (1.4).
4. **Skew double core.** It is class `A_0` (no residue data). By item 2, giving it genuine labels satisfies
   everything L4^forced sees. Whether it can also meet the P3-all\* slab conditions is C2 of round 9, still
   undone. Nothing in this write-up bears on it.

**Conclusion of the fake test.** The write-up's own claims pass it: the Theorem-2 systems are not fakes, and no
growth is claimed. The surviving question L4^forced uses no fake-violated input (item 2). Proving it would not by
itself beat the barrier. Any growth proof through route 3 has to get its fake-violating power from the
label-angle link (L7 / slab conditions) in regime (1.4), with `A` chosen after `c_0`.

## 8. Issues list

**Fatal:** none. Theorems 1 and 2 and Lemma P are correct. The refutation of L4(K) as posed stands.

**Major.**
* **M1 (threshold).** "`K >= 2L(h_0 + h)` (twice the README estimate)", labelled proved in §5 item 2, is wrong
  as a necessity claim. Linnik's `C_L k^L` is an upper bound for the least prime, not a lower bound. It is
  also the wrong comparison: the Hamming conclusion of L4/L4^forced feeds Theorem B, which needs only
  `Y >= p_j Norm zeta`. So `K ~ 2(h_0 + h)(1 + o(1))` suffices for the analytic step. The rigid-core form
  concludes the genuine class directly and needs neither theorem. Its required length (§6.1 of the write-up:
  `q_core >= C_L(8 M^(h_0+h))^L`) is set by L7, not by Linnik. Affects §0 item 3, §4.1, §5 item 2, §6.1.
* **M2 (the corrected statement).** The displayed L4^forced(K) uses `F_K = {n_q <= A_K}`, which is mostly free
  moduli. Independent instances give coherent share `<= 0.13-0.25` of `theta_7(F_K)`, so it fails the first
  test. It also keeps only norm-`<= K` budgets, discarding the norms `K' > K` that actual circles and P3-all\*
  provide and that are also forced on `F_K^ex`. In §6.3, norm 5 makes 4 of the 6 tested `K = 4` top moduli
  one-modulus rigid, cuts the non-genuine survivors at the other two by factors `10^3`-`10^4`, and defeats
  every randomized joint search. So the "most of the mass fails at `K = 4`" reading (§0 item 4, §4.4, §6.1)
  concerns the norm-truncated statement only. The corrected question should be stated on the exact forced set,
  with all norms (or norms up to `K + O(1)`) as hypothesis (§9).

**Minor.**
* **m1.** (2.1) is undefined at `q = p_j > q_T` and at `q = 2`. State "`q > q_T`, `q` not dividing `2N`", with
  genuine data at `q | 2N`.
* **m2.** The "Theorem-2 instances" of the summary are a different, random construction. Theorem 2's own
  construction is useless at `r <= 8`.
* **m3.** Non-unit convention: Thms A/B exclude `q | Norm zeta`. The theorem survives with `delta + o(1)`; say
  so.
* **m4.** Rigidity at `n_q <= 2` (`q = 3, 5, 7`) is vacuous, since the genuine class is every parity vector
  there. "`K = 2` rigid only at `{3, 7}`" means no non-vacuous rigidity at all.
* **m5.** "Near-extremal clusters have `r = Theta(M)` (`stability.md`)" is unsupported. Only `r = O(M)` is
  proved (`p3lower.md` D.0). It is also unnecessary: Theorem 2 needs only `r >= C(K, h, delta) log M`.
* **m6.** Theorem 1(a) for odd `K` improves by a factor of about 2 via `S_h ∪ (e_1 + S_h)` (constants only).
* **m7.** Remark 2: the norm-`(2h+2)` collisions are not "about as often as for random labels". They occur at
  every level at a positive proportion of free moduli (verified). The conclusion is stronger than stated and has
  a two-line proof.
* **m8.** Lemma P and "genuine class" are level-1 notions; say so where coherence is equated with the collision
  set (§4.1, §5 item 3).
* **m9.** The abstract's "(or below any `Y <= M^O(1)`)" needs `Y/q_T -> infinity` (e.g. `Y >= Q_K`), and
  "delta-fraction" should be "2 delta".
* **m10.** In the `K = 4` flexible-top instance, the system-level coherent share of the forced `Q_7` mass lies
  between 0.5 and about 0.71 (pair bound, `H = 2000`). "Fails in most-of-the-mass form" needs "most" to mean
  more than about 0.71. At the Theorem-B threshold `(1 + delta)/2` this instance is not a counterexample.
* **m11.** §5 item 5 should state the fake test explicitly (§7 item 2). L4^forced's inputs hold for every fake
  circle with genuine labels, so its fake-violating power, if any, must come from the label-angle step.

## 9. Corrected statement

**Theorem 1 (proved; as in the write-up).** Let `h = floor(K/2)`.
* If `n_q < |S_h|`, every label vector at `q` has a level-1 collision of norm `<= 2h`. For odd `K` the bound
  improves by a factor of about 2.
* If `n_q > A_K = K(2(K+1)^2 B + 2K + 1)`, the explicit parity-respecting vectors (2.1) have no collision of norm
  `<= K` at any level, for every permutation `sigma`.
* `B = r` for `K <= 3`. For larger `K`, `B < (2r)^h` by Bose-Chowla, or `B <= h(r + h)^(2h-1) + 1` by greedy.

Bounded norms `<= K` therefore force nothing above `M^(floor(K/2) + o(1))`, which sharpens README Prop. H.

**Theorem 2 (proved, with repairs m1, m3, m5, m9).** Fix `K >= 2`, `h >= 0` and `delta in (0, 1)`. Let `r` and
`Y` satisfy `r <= rho M`, `r >= C(K, h, delta) log M`, `Y >= Q_K` and `log Y = O(log M)`. Let `p_1, ..., p_r` be
distinct primes `= 1 mod 4`.

Build the system as follows: genuine data at every `q <= q_T = 4A_K + 5` and at every `q | 2N`; at every other
`q > q_T`, the vector (2.1) with independent uniform `sigma_q`. For `M >= M_0(K, h, delta, rho)` some
realisation has the following properties.
* It is parity-respecting.
* It is good, i.e. per-unit demand at most the genuine demand, for every monomial of norm `<= K`, at every level
  and every `X`.
* It is not system-level coherent with any `zeta` of norm `<= M^h` on any set of `Q_7` moduli carrying more than
  `2 delta + o(1)` of `theta_7([1, Y])`. Moduli dividing `Norm zeta` may be counted either way.

Hence **L4(K) as posed in `round9/README.md` §6 is false for every `K >= 2`**, with any primes in regime (1.4).
These systems violate the norm-`(2h+2)` budgets at every level at a positive proportion of the free moduli, so
they are not in P3-all\* and are not fakes.

**Lemma P (proved; level 1).** As in the write-up. Equality of the full collision lattices suffices without
short generation.

**What survives (corrected question; open).** Not the displayed L4^forced(K), which fails at the first test
(M2). The question consistent with every computation is an inverse step for the **exact** forced set that uses
the budgets the program actually has:

> **L4\*(K).** Regime (1.4), `r <= rho M`, `h >= h_0`. Suppose a parity-respecting system is good for
> **every** monomial of norm `<= K'` at all moduli of `F_(K')^ex`, the exact forced set of the largest norm
> used, where `K' >= K + 1` (or all norms at all moduli `<= X`, as in P3-all\*). Must it be in the genuine class
> (one shared unit, conjugation per column allowed) on an initial segment `q <= q_core` of `F_K^ex`? Or, in the
> Hamming form, coherent with some `zeta` of norm `<= M^h` on `(1 + delta)/2` of `theta_7(F_K^ex)`?
>
> The budgets must be charged on `F_(K')^ex`, not only on `F_K^ex`. The demand of the moduli in
> `F_(K')^ex \ F_K^ex` matters: `ref_deform_mix.py` at `q_0 = 107` is rigid when the genuine demand on
> `(199, 349]` is charged, flexible when it is not.

What each form would feed:
* **Hamming form.** It feeds **Theorem B**, which needs only `max F_K^ex >= M^((h_0 + h)(1 + eps))`, i.e.
  `K ~ 2(h_0 + h)`. It does not need `2L(h_0 + h)`.
* **Genuine-class form.** It needs neither Theorem A nor Theorem B. Its required length is set by L7 (the
  label-angle step), not by Linnik.

Evidence for L4\*(K):
* `K = 4` to `5`, `r = 4`, smallest primes (see §6.3): norm-5 budgets make 4 of the 6 tested moduli of the
  `K = 4` "flexible top" one-modulus rigid, and leave 40-96 non-genuine survivors at the other two (167, 199).
* Randomized joint searches for a non-genuine top fail there (*evidence*).

Evidence against the norm-truncated version:
* the write-up's `K = 4` flexible top, which I reproduced and which survives genuine labels above the forced
  range.

Not known:
* whether the rigid segment grows polynomially in `M`;
* whether rigidity survives primes with `log p_j/log r` bounded away from 1 (the loose `K = 4` set is far less
  rigid);
* whether rigidity holds jointly.

**Files** (`referee_L4K_checks/`; each runs as `python3 <script> [args]` from that directory):

| script | content |
|---|---|
| `rlib.py` | tori, labels, parity, `D_s`, budgets |
| `ref_genuine.py` | genuine collision iff `q^a \| D_s` |
| `ref_thm1.py` | Theorem 1(b) |
| `ref_free.py` | exact forced/free status (`K = 2, 3`; `_K4`, `_K5` variants) |
| `ref_forced_sig.py` | exhaustive `K = 2` count |
| `ref_sigclass.py` | genuine signature vs genuine class |
| `ref_instance.py` | Theorem-2-style instances, including `F_K = {n <= A_K}` shares |
| `ref_remark2.py` | the Theorem-2 systems violate norm-4 budgets at all levels |
| `ref_deform.py`, `ref_deform2.py`, `ref_deform_mix.py`, `ref_deform_ext.py` | one-modulus deformations |
| `ref_core.py` | `K = 4` joint core |
| `ref_top.py`, `ref_top2.py`, `ref_top5.py`, `ref_top5b.py`, `ref_topQ7.py` | flexible-top searches |

Outputs are in `*_out.txt`, `ref_deform2_K5_*.txt` and `ref_top*_out.txt`.

**Barrier status.** Unchanged.
* No growth improvement and no new cap. `M <= (2/3 + eps) log R/loglog R` stands. The P3-all\* bracket
  `3 M log M - O(M) <= W_min^all <= (66A + o(1)) M log^2 M/loglog M` is unchanged.
* L4^forced and L4\* use no fake-violated input. Every fake circle with its genuine labels satisfies their
  hypotheses and conclusions (§7 item 2). Route 3's fake-violating power must come from the label-angle step
  (L7 / slab conditions) with `A` chosen after `c_0`.
