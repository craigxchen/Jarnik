# Referee report on `round7/p3G.md` (route G: genuine residues, every monomial strengthened)

Referee: adversarial, independent. Scripts are in `round7/referee_p3G_checks/`. None of them imports the
author's code. The one exception is `author_code_on_adaptive.txt`, which runs the author's closure on
the referee's prime sets as a cross-check. The repository and the round-5/6/8 folders were not
modified.

References: `§n` is a section of `p3G.md`, and `R§n` is a section of this report.

Labels: **proved** means re-derived here with all quantifiers. **Exact** means an exact finite
computation, every certificate re-checked. **Evidence** means a computation no proof uses.

---------------------------------------------------------------------------------------------

## 0. Verdict

**The mathematics is correct. The headline is not established.**

* *Correct.* Every lemma, theorem and proposition of §§1-3, §5.1-5.2 and §6.1 was re-derived. No
  mathematical error was found.
* *Computations reproduce.* Theorem 4.1 reproduces exactly with independent code: same prime sets,
  same tight counts, complete closure, every certificate re-checked by an independent lattice
  checker. The completion sums 15190 and 4052 also reproduce exactly.

**What does not survive is the headline.** The headline reads: "Route G does not give fake circles
at `W = O(M log M)`", or "on prime sets of polynomial size the class contains no fakes at all".
* *What is proved.* The collapse is shown for six explicit random prime sets, used with all of
  their columns. Beyond that it holds with probability `1 - o(1)` in a model with independent
  tightness and **non-adaptive** prime choice.
* *What else the lower bound needs.* It also uses the heuristic near-relation count (Remark 3.5),
  and it holds for squarefree members only.
* *Counter-evidence (new; exact for each stated sub-closure, evidence for the question).* A builder
  who **selects** primes, instead of drawing them, defeats the closures actually computed in the
  paper, at the very parameters where random sets collapse.
  * At `M = 644`, `r = 643`, `p in [M^4, 2M^4]`, `X = 2M`, a greedily selected set has an empty
    weight-`<= 2` closure: all 643 columns are free. The random set with the same parameters is a
    verified collapse case of Theorem 4.1. The same holds at `M = 1284`, `K = 4` (random set:
    full collapse).
  * At `M = 164` the same holds down to `p in [M^(2.75), 2M^(2.75)]`, i.e. `lambda/H = 0.88`, with
    pair degree 82 (random: 88, full collapse). That set needs 2 columns deleted to clear a few
    double-pair `k = 3` certificates. The paper's measured threshold is `lambda*/H = 1.29`.
  * "Empty" here means three things:
    * no certificate of the author's `k <= 2` closure;
    * no closed-walk circuit of any length or multiplicity among single and pair relations;
    * after deleting at most 3 columns where needed, no certificate with `<= 3` terms among all
      weight-`<= 2` relations (R§4.2).

    The first was confirmed with the author's own closure code on the five selected sets tested
    (0 certificates each).
  * Weight-3 relations restore much of the pinning. The weight-`<= 2` greedy sets at `K <= 3.5`
    are fully pinned once all tight triples are added. But deleting 4-15 more columns again gives
    an empty weight-`<= 3` (`k <= 3`) closure on 148-159 of 163 columns at `K = 3.0`-`3.5`. Random
    sets are fully pinned there, and still at `K = 3.75` (R§4.3).

So the thresholds of §4.4 (`lambda* ~ 1.23-1.34 H`) are properties of **random** prime sets, not a
size a builder must reach. The premise of the paper's adaptive heuristic (§5.3(3)) is the wrong
combinatorial condition. Each added weight restores some pinning, and the full closure on selected
sets is not computed. Whether `P3-G` has fakes at `W = O(M log M)` therefore remains **open**. No
fake is constructed here.

Severity summary:

| # | issue | severity |
|---|---|---|
| A | Headline and §0/§8 claims ("Route G does not give fakes at `O(M log M)`", "prime sets of polynomial size contain no fakes", "a fake needs `~M/2` columns of size `c H`") stated as established. They rest on a non-adaptive model, a heuristic near-relation count and the squarefree restriction. | major |
| B | The paper's evidence is for random sets only. Selected prime sets have an empty weight-`<= 2` closure at the verified parameters (`M = 644`, `K = 4`) and down to `lambda/H = 0.88` (`M = 164`). After pruning, they have an empty weight-`<= 3` (`k <= 3`) closure where random sets collapse. The premise of Heuristic 5.3(3) (an independent set in the tight graph) is the wrong condition: builders keep pair degree 26-82 with zero certificates. | major |
| C | §4.2 gloss: Theorem 4.1 "covers every `b`-copy Paley profile with `b(M-1) <= r` built on these primes, and every round-5/6 construction transplanted onto them". False for `b(M-1) < r`: 163-subsets of the verified `M = 164, b = 8, K = 4` set leave 35-54 columns free. | major (scope of a verified claim) |
| D | §7: "`P3-G-char` does **not** collapse" is asserted but not proved. Only "the collapse mechanism is unavailable" is shown. | minor |
| E | "`lambda* ~ H_X`, in rough agreement with the measured `1.06-1.22 H_X`": this compares a pair-only measurement (own heuristic limit `(2/3) H_X`, measured ratio rising) with the full-closure heuristic. | minor |
| F | Citation: `p3data.md` gives `D_max/H` = 2.2-3.2 at `X = 2M` (2.9-3.7 at `X = 4M`), not "3.2-3.7 at `X = 2M`". "Within 4%": one row is off by 4.05%. Prop. 6.1's constant is `66 + o(1)` (since `psi(X) ~ X`); 68.7 is valid but not sharp. | minor |

A positive addition from the referee is the **arc quantisation lemma** (R§5). Relative to
`Pi_0 = {0}`, the arc sum of any zero combination is a multiple of `pi/2`. So level-0 certificates
have arc sum exactly `pi/2`, and the level-0 certificate test is purely combinatorial.

---------------------------------------------------------------------------------------------

## 1. Proofs re-derived

Each item was re-derived from the definitions (`sharp.md` Def. 1.1, `outside.md` (1.1), Lemma 0.1
and Thm A.1, `p3lower.md` Def. 1.1). Verdict: **correct** unless noted.

* **Definition 1.1** agrees with `p3lower.md` Def. 1.1 (P3-gen). The inclusion
  `P3-G_(X') ⊂ P3-G_X` for `X' >= X` holds: a member at `X'` reduces to one at `X`.
* **Lemma 1.2.** This is the identity `A_v / conj(A_v) = prod_j (pi_j / conj pi_j)^(v_j)`, together
  with:
  * `A - i^s conj A` = `2iY`, `(1 - i)(X - Y)`, `2X`, `(1 + i)(X + Y)`;
  * `q` odd, so `2` and `1 +- i` are units mod `q`;
  * a rational integer divisible by `q^a` in `Z[i]` is divisible by it in `Z`.

  Hence `m_(v,s) = gcd(D_s, L'_X)`. `D_s != 0` by conjugate-primitivity, and `X + Y` is odd.
* **Remark 1.5.** With every `b_j` even, `Y_v` is even for all `v`. So `D_0` is never odd and unit 0
  is never tight. Correct.
* **Lemma 2.1** (planted point). Immediate from `m | D_s != 0`.
* **Lemma 2.3.**
  * (1) checked for all four units, including `|sin(t - pi/4)| = |X - Y|/(sqrt2 |A|)` and
    `|sin(t - 3pi/4)| = |X + Y|/(sqrt2 |A|)`.
  * (2) `sin 2 alpha = Im(A^2 (-i)^s)/|A|^2`. Expanding gives the four sign formulas.
  * (3) With `sigma = +1`, `alpha = d mod pi`, and `dist(d + x, pi Z) >= d` iff
    `x mod pi in [0, pi - 2d]`. The case `sigma = -1` is symmetric.
  * The referee's checker uses `Im(A^2 (-i)^s)` and `Im(A (1-i)^s)^2 / (2^s |A|^2)`. These are
    formulas independent of the paper's, and they agree on every certificate.
* **Lemma 2.4.**
  * `sigma_i x_i = y_i + pi k_i` with `y_i in [0, pi - 2 d_i]`, and `sum n_i sigma_i x_i in pi Z`.
    So `0 <= sum n_i y_i < pi` forces every `y_i = 0`. Correct; it needs `n_i >= 1`, as stated.
  * The special cases `sin^2 d_1 + sin^2 d_2 > 1` (`k = 2`) and `sin^2 d > 1/2` (self) are
    equivalent to the arc condition, since `d_i in (0, pi/2)`.
* **Corollary 2.6.** Tightness at `X` implies tightness at `X' >= X`, with the same
  `(sigma, d)`, because `L'_X | L'_(X')` for the same `N`.
  * *Invariance (checked).* Replacing `pi_j` by `i pi_j` shifts `s` by `2 v_j` and leaves `alpha`,
    `d`, `sigma` unchanged. Replacing `pi_j` by `conj(pi_j)` maps `v -> R v` with `R` a sign flip,
    and `A_(Rv)(pi') = A_v(pi)`.
  * So `Pi*` is invariant under unit changes and covariant under conjugation, and the conditions
    `2 e_j in Pi*` and `e_j in Pi*` depend only on the rational prime set. The paper uses this
    implicitly; it is true.
* **Theorem 3.1** (collapse).
  * `phi_j = theta_j + (pi/2) t_j` lifts `arg(i^(t_j) pi_j)`.
  * The product has modulus `N^(1/2)` and argument `<2 a_x - e, phi>`. Unique factorisation in
    `Z[i]` gives distinctness.
  * The arc has length `R Delta = C R^(1/2)`.

  Correct, for **every** exponent vector `e`.
* **Theorem 3.2.**
  * The summed identity is `<2 v(c), phi> = c.delta + (pi/2) c.k`, and `|c.delta| <= (n/2) Delta`.
  * `dist(arg A', s pi/4 + pi Z) = arcsin(nu_s |D_s|/|A'|)`, and `sin x <= x`.

  Correct, for every `e`.
* **Proposition 3.3.**
  * PSD rank bound `(tr G)^2 / tr G^2`.
  * `|mu| <= W/2` and `mu >= -c_C` (pair inequality).
  * Cut identity `sum_(x<y) mu_xy <= M W/4`.
  * `rank G <= rank A + 1 <= rank V + 2`.
  * `W >= 3 M log M - O(M)` (`sharp.md` Prop. 1.3, which needs `X >= 2M`).

  Correct, for `e = 1`.
* **Corollary 3.4.** `C_F -> V cap {v|_F = 0}`, with rank `>= rank V - |F|`. Correct, for **squarefree**
  members only, as stated in the corollary.
* **Lemma 5.1.** Both halves were checked: `d(beta,s) + d(beta',s') = pi/2 + beta' - beta`, with
  opposite `sigma`.
* **Lemma 5.2.** The closed form `(1-rho)^(m-1)(1+(m-1)rho)` was re-derived by first-success
  decomposition. Exact enumeration with Fractions (`m <= 7`, three values of `rho`) gives 0
  mismatches (`lemma52_quant_out.txt`). The monotonicity ratio bound is correct.
* **Theorem 5.3** (in the model).
  * *Seeds.* Units (1,2) on `(3pi/4, pi)` and (2,3) on `(0, pi/4)`: opposite `sigma`, and
    `d + d' = 3pi/4`.
  * *Probability.* Seeds have probability `rho_1^2/2`, independently over `j in A`. Chernoff with
    `mu >= 24 log f` gives the `f^(-3)` term.
  * *Propagation.* `sigma_a v_a + sigma_b v_b = e_a - e_b in <E_0>`. The bound
    `(1 + x) e^(rho_2 - x)` with `x >= 3 log f` is valid.
  * *Union bound.* `e^(1/2)(1 + 3 log f)/f^2 + f^(-3) <= (5 log f + 2)/f^2` for `f >= 3`.

  Correct. The "bounded below" extension needs a monotone coupling and **independence** across
  relations; the model assumes both.
* **Corollary 5.4** (in the model).
  * Hildebrand's range covers `u <= log X`. The powersmooth correction is `o(1)`, and
    `rho_1^2 rho_2 >= M^(-1 + 2 eps + o(1))`.
  * Theorem 5.3 needs `rho_2 <= 1/2`. For small primes `tau_X = 1`, which is handled by the
    monotonicity remark with `rho_2 := 1/2`.
  * The final `W` bound needs Corollary 3.4 with `k = O(M/log M)` near-relations. That is a
    **hypothesis**, supported only by the heuristic count of Remark 3.5, and it restricts to
    squarefree members.

  Correct as a conditional statement in the model.
* **Proposition 6.1** (fakes at `W = O(M X)`).
  * Step 1: pigeonhole on `[Z, 2Z]`.
  * Step 3: the `Fail_c` sufficiency covers all `eta` and all `s` via
    `dist(c.delta/2, (pi/4) Z)`. The union-bound arithmetic, with `tau >= 2 log L`, `L >= e^((2 - o(1))M)`,
    and the pairs, ternary and amplitude classes, was checked.
  * Step 4 uses `outside.md` Thm A.1 on the connected coset `phi_0 + U`, where `<v, .>` is
    Haar-uniform for `v` outside `L`. Summing over all `v != 0` only over-counts.
    `4 L (e^(2r/(sqrt Z - 1)) - 1) ~ 1/8 < 1/4`.
  * Steps 5 and 6 are correct.

  The constant: `psi(X) = (1 + o(1)) X`, so the bound is `(66 + o(1)) M X`. The stated `68.7`
  (from `psi < 1.04 X`) is valid but not sharp.
* **§7, `P3-G-char`.**
  * Proved: integral characters have `w(v(c)) >= W/2 - 2 log(nC/2^(3/2))`. So weight-`<= 2`
    monomials are not characters and carry no residue strengthening.
  * Not proved: that `P3-G-char` "does not collapse". Showing that would need a non-actual member
    (or at least an incomplete closure in the `P3-G-char` sense). The text should say "the
    collapse mechanism of §2 is unavailable".

---------------------------------------------------------------------------------------------

## 2. Computations re-run (exact)

**Code.** `rcore.py` and `rclosure.py` were written from scratch:
* primality by deterministic Miller-Rabin, cross-checked against `sympy.isprime` on `2*10^4`
  random integers below `10^15` with 0 mismatches;
* `a^2 + b^2 = p` by Hermite-Serret;
* tightness by the exact test `L'_X mod |D_s| == 0`;
* `sigma` and `sin^2 d` exact;
* a closure using a quotient-group representation, with signed components of status
  `Z`, `Z/2` or `0`, and buckets keyed by the class modulo `Pi`.

The checker re-verifies every certificate from scratch:
* exact Gaussian products, and tightness;
* `sigma` from `Im(A^2 (-i)^s)`;
* `sin^2 d` from `Im(A (1-i)^s)^2 / (2^s |A|^2)`;
* the arc condition as exact rationals;
* membership of the circuit vector in the lattice generated by all earlier added vectors, using a
  sparse echelon (HNF-type) basis over `Z`;
* finally `e_k` and `2 e_k` membership for every column.

**Theorem 4.1 reproduced** (`run_cases.py` -> `cases_A.jsonl`, `cases_B.jsonl`). The prime sets
come from the same generator description and coincide with the author's: `pmin`, `pmax` and the
single/double/pair tight counts are identical in every row.

| M | b | K | r | tight S / T / P | referee certificates | checker: `e_k` in lattice | collapse |
|---|---|---|---|---|---|---|---|
| 44 | 2 | 3 | 86 | 149 / 27 / 1235 | 89 | 86 / 86 | yes |
| 84 | 8 | 3 | 664 | 1130 / 225 / 76601 | 1484 | 664 / 664 | yes |
| 164 | 8 | 3 | 1304 | 2288 / 463 / 294621 | 1710 | 1304 / 1304 | yes |
| 164 | 8 | 4 | 1304 | 1240 / 139 / 38391 | 1307 | 1304 / 1304 | yes |
| 324 | 2 | 4 | 646 | 638 / 67 / 9911 | 646 | 646 / 646 | yes |
| 644 | 1 | 4 | 643 | 589 / 58 / 9131 | 643 | 643 / 643 | yes |
| 644 | 1 | 4.25 (seed 1) | 643 | 557 / 53 / 5159 | 636 | 636 / 643 | no (7 free) |
| 1284 | 1 | 4.5 | 1283 | 816 / 56 / 11623 | 1279 | 1279 / 1283 | no (4 free) |

The certificate counts differ from the author's (89 vs 86, and so on) because the firing order
differs. The final lattices agree. **Theorem 4.1, as stated (members whose prime set is exactly
the listed set), is confirmed.**

**Completion sums** (`completion_check.py`). The tight weight-`<= 2` slab measures `sum 2d/pi` are
15190.1 (`M = 164, b = 8, K = 4`) and 4051.9 (`M = 644, b = 1, K = 4`). Both reproduce the paper.

**Density table (§4.3).** This was not re-run; it is evidence only. Arithmetic on the quoted
numbers gives a largest relative deviation of 4.05% (355 vs 370), not "within 4%". The difference
is immaterial.

---------------------------------------------------------------------------------------------

## 3. Issue C: Theorem 4.1 does not cover profiles on a subset of a verified set

§4.2 says: "Since the closure does not depend on the profile, Theorem 4.1 covers, for instance,
every `b`-copy Paley profile with `b(M-1) <= r` built on these primes, and every round-5/6
construction transplanted onto them."

The closure does not depend on the profile. It does depend on the **prime set**.
* A member whose profile uses `b(M-1) < r` columns has a smaller prime set `S' ⊂ S`.
* Its closure `Pi*(S', X)` is built only from relations supported on `S'`. This set is in general
  much smaller than the restriction of `Pi*(S, X)`, because the relations through the dropped
  columns are gone. (`L'_X` is unchanged, since all dropped primes exceed `X`.)
* Theorem 3.1 needs completeness of `Pi*(S', X)`.

*Exact counter-check* (`subset_test.py` -> `subset_out.jsonl`). Start from the verified set
`M = 164, b = 8, K = 4` (`r = 1304`, complete closure).

| subset | size | tight relations | certificates (all re-checked) | free columns |
|---|---|---|---|---|
| first 163 primes | 163 (`b = 1` Paley at `M = 164`) | 863 | 111 | **52** |
| random 163 (3 draws) | 163 | 765-795 | 109-128 | **35-54** |
| first 326 | 326 (`b = 2`) | 3033 | 325 | 1 |
| random 326 | 326 | 2789 | 317 | 9 |

So the computed closure does not certify collapse for `b = 1` or `b = 2` Paley members built on
these primes. Free in the `k <= 2` sub-closure does not prove free in the full closure. But the
gloss claims these members are *covered by Theorem 4.1*, and they are not.

"Every round-5/6 construction transplanted onto them" fails for a further reason. The round-5/6
constructions at `M = 44` use `r = 33 * 43 = 1419` primes, more than any verified set at that `M`
has (86). Transplanting them means using a different prime set.

**Fix.** Delete the sentence, or restrict it to `b(M - 1) = r`, i.e. members whose prime set is
the whole verified set.

---------------------------------------------------------------------------------------------

## 4. Issue B: selected prime sets (evidence; exact for each stated sub-closure)

The paper's evidence for "prime sets of polynomial size are pinned" consists of three parts:
* Theorem 4.1, six random sets;
* §4.4, random sets, pair-only closure;
* Theorem 5.3, a random model.

A builder of fakes does not draw primes at random. It selects them. The paper treats this case in
one heuristic paragraph, §5.3(3):

> "A builder may pick primes from the pool ... to avoid certificates. With C1-type certificates
> only, this is at best an independent set in a random graph of density `~ 6 rho_2` on the pool ...
> So selection changes `lambda*` only by a factor `1 + o(1)`."

**The premise is the wrong combinatorial condition.**
* A builder need not avoid tight pairs. It must avoid **certificates**.
* By Lemma R.1, a certificate at level 0 is a zero combination with arc sum exactly `pi/2`.
* Only a small fraction of the cycles in the tight graph are certificates.

The experiments below show that a builder keeps the tight-pair degree at 26-82, essentially the
random value, with **zero** weight-`<= 2` certificates.

### 4.1 Method

`adaptive.py`, `adaptive_ties.py`.
* Draw random primes `p = 1 mod 4` in `[P, 2P]`, `P = M^K`, with `X = 2M`.
* Accept a prime iff the accepted set plus it has **no level-0 certificate** of these kinds:
  * same-`v` pairs, for singles, doubles and pairs (the author's `k <= 2` certificates at level 0);
  * double/single with `n = (1, 2)`;
  * any closed walk, of any length and multiplicity, of single and pair relations in the signed
    double cover `{columns} x {+-1} ∪ {ground}`, with arc sum `< pi`.
* A closed walk has `sum n_i sigma_i v_i = 0`, so it is a Lemma-2.4 certificate relative to `{0}`
  iff its arc sum is `< pi`, i.e. exactly `pi/2` by Lemma R.1. Floating error is therefore
  harmless.
* `adaptive.py` also rejects "ties", whose arc sum is exactly `pi`. This is unnecessary;
  `adaptive_ties.py` accepts them.

If none of these exists, then `Pi_1 = {0}` for the closure built from these certificate families,
and every column is free.

**Independent confirmations.**
* `eval_set.py`: the referee's `k <= 2` closure, with the exact checker.
* `author_code_on_adaptive.txt`: the author's own `gtight.closure` gives 0 certificates and all
  columns free on every selected set tested: `M = 164` at `K = 3.0`, `3.5`, `2.75`, `M = 324` and
  `M = 644`. Its tight counts equal the referee's.
* `brute_k3.py`: brute force over all level-0 certificates with `<= 3` terms, `n_i in {1, 2}`,
  among all weight-`<= 2` relations, doubles included. It finds 43 certificates on a random set
  (`K = 3.5`), **0** on the `adaptive.py` sets, and 4 resp. 3 on the two tie-accepting sets.
  These are treated below the table in R§4.2.
* `gclosure.py`: a closure modulo an arbitrary sublattice (canonical coset representatives), with
  `k <= 2` and `k = 3` certificates, every certificate re-checked.

### 4.2 Weight `<= 2`: selection defeats the computed closure

| M | K | `lambda/H` | candidates drawn | accepted | pair degree (random) | free, selected | free, random (same `M, K`) |
|---|---|---|---|---|---|---|---|
| 164 | 3.75 | 1.20 | 252 | 163 | 11.6 (12.0-12.8) | 163/163 | 3-8/163 (author) |
| 164 | 3.5 | 1.12 | 404 | 163 | 18.7 (21.2) | 163/163 | 0 |
| 164 | 3.25 | 1.04 | 2360 | 163 | 30.7 (35) | 163/163 | 0 (author) |
| 164 | 3.0 | 0.96 | 79625 | 163 | 46.6 (56.9) | 163/163 | 0 |
| 164 | 2.75 (ties accepted) | 0.88 | 13456 | 163 | 82.2 (87.7) | 163/163 (\*) | 0 |
| 164 | 2.5 (ties accepted) | 0.80 | 12516 (pool exhausted) | 96 | -- | 96/96 | -- |
| 324 | 3.5 | 1.06 | 9961 | 323 | 37.7 (43.4) | 323/323 | 0 |
| 644 | 4.0 | 1.15 | 2093 | 643 | 26.4 (28.4) | 643/643 | **0** (Theorem 4.1 row) |
| 1284 | 4.0 (ties accepted) | 1.10 | 3558 | 1283 | 54.8 (57.4) | 1283/1283 (\*) | **0** (collapse) |

"Free, selected" is the author's `k <= 2` closure. It is also free for every closed walk.
* For the sets of `adaptive.py` (ties rejected: `M = 164` at `K >= 3.0`, `M = 324`, `M = 644`),
  `brute_k3.py` finds **no** level-0 certificate with `<= 3` terms either, doubles included.
  So `gclosure.py` (weight `<= 2`, `k <= 3`) is empty on them too.
* (\*) The two tie-accepting sets contain a few `k = 3` certificates that combine a double `2 e_j`
  with pairs: 4 at `M = 164`, 3 at `M = 1284`. These are outside the author's `k <= 2` closure,
  and the closed-walk check ignores doubles.
  * On the `M = 164` set they seed a full collapse of the general closure.
  * Deleting 2 (resp. 3) columns removes them (`prune2.py`). That leaves **161/163** (resp.
    **1280/1283**) columns with an empty weight-`<= 2`, `k <= 3` closure.

In the weight-`<= 2` closure the paper computes, the measured threshold is `lambda*/H = 1.29` at
`M = 164` (`K* = 4.03`) and `1.33` at `M = 644`. Selection pushes it below `0.88` and below `1.15`
respectively, and below `1.10` at `M = 1284` (random: `1.34`). The `M = 644` row has the same parameters as a verified collapse case of Theorem 4.1:
the random set collapses, while the selected set is entirely free.

### 4.3 Weight 3: higher weight restores much of the pinning

The greedy controls only weight `<= 2`. Adding all tight weight-3 relations `e_i +- e_j +- e_k`,
the general closure (`gclosure.py`, `k <= 3`, all certificates re-checked) gives:

| M = 164 set | tight relations (weight `<= 3`) | `k = 3` certificates used | free columns |
|---|---|---|---|
| selected, `K = 3.75` | 2254 | 1 | **162**/163 |
| selected, `K = 3.5` | 4629 | 4 (`k = 3`) | 0 |
| selected, `K = 3.25` | 10215 | 10 | 0 |
| selected, `K = 3.0` | 22137 | 0 (seeds: same-`v` triple pairs, `k = 2`) | 0 |
| random, `K = 3.75` | 2376 | -- | 0 |
| random, `K = 4.0` | 1216 | -- | 4 (weight `<= 2`: 34) |
| random, `K = 4.5` | 371 | -- | 160 |
| selected `K = 3.5`, then pruned (`prune.py`) | 4345 | 0 | **159/159** (4 columns deleted) |
| selected `K = 3.25`, then pruned | 9225 | 0 | **157/157** (6 columns deleted) |
| selected `K = 3.0`, then pruned | 16994 | 0 | **148/148** (15 columns deleted) |

So weight-3 relations supply the seeds that the weight-`<= 2` greedy removed, and the dense tight
pair graph then propagates them to every column. The paper predicts exactly this: "higher weight
helps", §5.3(2). A builder can answer by deleting a few more columns. At `M = 164` the
weight-`<= 3` (`k <= 3`) sub-closure is then empty on 159 columns at `K = 3.5`, on 157 at
`K = 3.25`, and on 148 at `K = 3.0`. Random sets pin fully at `K = 3.75`, and they are 98% free only at `K = 4.5`.

### 4.4 What this shows, and what it does not

* **Shown** (evidence, exact for each stated sub-closure):
  * the thresholds of the paper's §4.4 are properties of **random** sets;
  * a builder lowers the critical size `lambda` of the weight-`<= 2` closure by a factor
    `<= 2.75/4.03 = 0.68` at `M = 164` (on 161 of 163 columns), `<= 4.0/4.61 = 0.87` at `M = 644`,
    and `<= 4.0/4.86 = 0.82` at `M = 1284` (on 1280 of 1283 columns). These are upper bounds on the factor, since lower `K` was not pushed at
    `M >= 324`.
  * For the weight-`<= 3`, `k <= 3` sub-closure at `M = 164`, the factor is `<= 3.0/4.0 = 0.75`,
    at a cost of 15 of 163 columns. Random sets at `K = 4.0` still have only 4 free columns there.
  * the independent-set premise of Heuristic 5.3(3) is wrong. Builders keep near-random tight
    degrees.
* **Not shown:**
  * that the **full** closure (all weights, all `k`) is evaded. Each added weight restores
    pinning, and weight-4 relations were not computed;
  * that any selected set carries a fake. Completion is a separate problem, §8.4(c) of the paper.
* **Asymptotics.** The paper's conclusion "selection changes `lambda*` by `1 + o(1)`" may well be
  right for another reason. Heuristically, a builder choosing each column from a pool of
  `~ P/log P` primes can absorb only `O(log P)` expected certificate events per column. That gains
  at most a factor polylogarithmic in `M` in the admissible tight density. Since the density decays like
  `exp(-u log u)`, that gain moves `u log u` by `O(loglog M)`, i.e. `lambda` by a factor
  `1 + O(loglog M/log M)`. But this is heuristic too, and the data show the effect is large at
  all computable `M`.

The headline's "prime sets of polynomial size are pinned" is therefore supported only for random
sets. For selected sets it is open, both by the paper's own open point 8.4(b) and by these data.

---------------------------------------------------------------------------------------------

## 5. Referee lemma: level-0 arcs are quantised

**Lemma R.1.** Let `(v_i, s_i)` be tight relations with data `(sigma_i, d_i)` (Lemma 2.3), and let
`n_i >= 0` be integers with `sum_i n_i sigma_i v_i = 0`. Then `sum_i n_i (pi - 2 d_i) in (pi/2) Z`.

*Proof.*
* By Lemma 2.3, `alpha_i := <v_i, theta> - s_i pi/4` satisfies `sigma_i alpha_i = d_i + pi k_i` with
  `k_i in Z`.
* The linear dependence gives the real identity
  `sum_i n_i sigma_i alpha_i = <sum_i n_i sigma_i v_i, theta> - (pi/4) sum_i n_i sigma_i s_i = -(pi/4) S`,
  with `S in Z`.
* Hence `sum_i n_i d_i in (pi/4) Z`, and `sum_i n_i (pi - 2 d_i) in (pi/2) Z`. QED

**Consequences.**
1. *Level-0 certificates have arc sum exactly `pi/2`.* Relative to `Pi_0 = {0}`, a family is a
   certificate (2.2) iff it is a zero combination with arc sum exactly `pi/2`, since the sum is
   positive and `< pi`. Every non-certificate zero combination has arc sum `>= pi`.
   * So the level-0 test is combinatorial and immune to rounding: any floating error below `pi/4` is
     harmless.
   * Examples: the seeds of Theorem 5.3 have `d + d' = 3pi/4`, i.e. arc sum `pi/2`. Units
     `s, s + 2` on the same `v` always give `d + d' = pi/2`, i.e. arc sum exactly `pi` (a "tie").
   * Checked on all same-`v` opposite-orientation combinations of three random sets: arc`/(pi/2)`
     takes only the values 1, 2, 3, within `1.4 * 10^(-13)` (`lemma52_quant_out.txt`).
2. *Farkas-type certificates of non-pinning cannot work at level 0.* Ties are zero combinations
   with arc sum exactly `pi`, and halving them formally gives `pi/2`. So no LP relaxation separates
   certificates from ties. The referee tried two such LPs (`farkas.py`, `farkas2.py`); both are
   infeasible even on random sets at `K = 5`, where nothing is pinned. Non-pinning at level 0 is
   an **integer** property.
3. At levels `t >= 1` (with `Pi_t != {0}`) the identity picks up `<u, theta>` for `u in Pi_t`, and
   there is no quantisation.

The lemma is elementary. It is not stated in `p3G.md`, whose verifier instead evaluates `sigma`
numerically (60-digit `mpmath`). That is safe, but unnecessary at level 0.

---------------------------------------------------------------------------------------------

## 6. Issue A: what the headline actually rests on

The paper states as established:
* (§0) "Route G does not produce fake circles at `W = O(M log M)`";
* (§0) "on prime sets of polynomial size the class contains no fakes at all";
* (§0) "A fake needs about `M/2` columns of size `log p >~ c log^2 M/loglog M` (proved with
  `c = 1/2 - eps` in the model)";
* (§3, end) "Sections 4 and 5 show that prime sets of polynomial size are pinned";
* (§8.1) "Prime sets of polynomial size are complete";
* (summary) "Route G therefore cannot answer P3-all\*".

What the paper proves toward these:

| ingredient | status | what it covers |
|---|---|---|
| Thm 3.1 (collapse given a complete closure) | proved | any prime set with a complete closure, any `e` |
| Thm 4.1 | exact | six explicit random sets, members using **all** `r` columns |
| Thm 5.3 | proved **in a model** | independent tightness given the angles; non-adaptive primes |
| model `=>` actual primes | **unproved** (open point 8.4(a)) | -- |
| adaptive prime selection | **unproved**: §5.3(3) is a heuristic count, and its premise is wrong (R§4) | -- |
| Cor. 3.4 / 5.4, "about `M/2` free columns" | proved for **squarefree** members, given `<= k` near-relations | `e = 1` only |
| "`O(M/log M)` near-relations" | **heuristic** (Remark 3.5) | -- |

So the negative answer for route G is established only on the six sets of Theorem 4.1. It is a
theorem in the random model, conditional on the near-relation count and on `e = 1`. It is not
established for the prime sets a builder would actually use, which are selected, not drawn.
R§4 shows that selection matters a great deal at every `M` tested (164 to 1284).

The statements "`P3-G` fake rate heuristically `Theta(M log^2 M/loglog M)`" and "`(1/2 - eps) H_X`
in the model" are fine as **labelled** heuristic and model statements. The §0 headline and §8.1
must be weakened to match the table above. The paper's own open point 8.4(b) (adaptive selection)
already concedes the decisive case. It cannot then be cited as answered in §0.

Non-squarefree members are not addressed by Cor. 3.4 at all, because Prop. 3.3 uses `e = 1`.
Theorem 3.1 does cover them, but only on complete prime sets.

---------------------------------------------------------------------------------------------

## 7. Classes, validity for actual circles, attribution

**Classes.**
* Every proved statement of §§1-3 is about members of `P3-G_X ⊂ P3-all ⊂ P3-char ⊂ A_X`.
* Actual clusters are members for every `X`, with `phi = theta` (Prop. 1.3, re-checked against
  `sharp.md` Prop. 1.2 and `p3lower.md` Prop. 1.2). So every inequality used is valid for actual
  circles.
* Theorems 3.1 and 3.2 are trivially true for actual clusters. They obstruct only non-actual
  members, as the paper says.
* The model statements (Thm 5.3, Cor. 5.4) are about the model `T`, not about actual primes. The
  paper labels them correctly in the body. They are not correctly scoped in §0 and §8.1 (Issue A).

**Validity of the inputs.**
* `X = 2M` in Theorem 4.1 covers every member at level `X' >= 2M`, in particular `sharp.md`'s
  `X = 3M^A`. It does not cover members at levels `X' < 2M`, and the paper claims nothing there.
* Theorem 4.1 is stated for members whose prime set is exactly the listed one. Correct as stated;
  the gloss is wrong (Issue C).

**Attribution.**
* Hildebrand's range `y >= exp((loglog x)^(5/3 + eps))`, the Janson-Luczak-Turova-Vallier
  bootstrap threshold `~ 1/(n p^2)` for 2-neighbour bootstrap percolation on `G(n, p)`,
  Rosser-Schoenfeld `psi(x) < 1.03883 x`, and Hecke equidistribution of Gaussian prime arguments
  are all cited correctly.
* The local estimate `Psi(Y/d, X)/Psi(Y, X) <= d^(-1) u^(O(log d/log X))` is the standard
  `d^(-alpha(Y,X))` behaviour (de Bruijn / Hildebrand-Tenenbaum), with `1 - alpha ~ log(u log u)/log X`.
  It is fine in the stated range.
* `p3data.md` is misquoted (Issue F).
* The novelty claim ("the mechanism is new") is plausible within this programme. The planted-point
  and saturation idea does not appear in rounds 4-6. No outside-literature claim is made.

---------------------------------------------------------------------------------------------

## 8. Minor issues (D, E, F)

* **D (§7).** "`P3-G-char` does **not** collapse" should read "the collapse mechanism of §2 is
  unavailable in `P3-G-char`, since weight-`<= 2` monomials are not characters". No non-actual
  member of `P3-G-char` at `W = O(M log M)` is exhibited, and no incomplete-closure argument for
  that class is given.
* **E (§5.3, §6).** "`lambda* ~ H_X`, in rough agreement with the measured `1.06-1.22 H_X`" is not a
  fair comparison.
  * The measurement is the `k <= 2`, weight-`<= 2` closure. By the paper's own item 1, its
    asymptotic heuristic is `(2/3) H_X`.
  * The measured ratio rises with `M` (1.06 to 1.22), away from `2/3`.
  * So the data neither confirm `lambda* ~ H_X` for the full closure nor fit the pair-only
    heuristic. They are pre-asymptotic, as the paper partly concedes ("we do not extrapolate").
* **F.**
  * (i) `p3data.md` §2.1 gives `D_max/H` = 2.2-3.2 at `X = 2M` and 2.9-3.7 at `X = 4M`. Heuristic
    6.2 quotes "3.2-3.7 at `X = 2M`".
  * (ii) "Every row within 4%": 355 vs 370 is -4.05%.
  * (iii) Prop. 6.1: `(66 + o(1)) M X` holds, since `psi(X)/X -> 1` as `X >= 2M -> infinity`. The
    stated `68.7` is valid but loose.
  * (iv) Paper §4.4: at `M = 164`, `K = 4.0` the three seeds give free fractions 0.25, 0.23 and 0.82. The
    interpolated `K* = 4.03` therefore carries a large seed-to-seed spread, which the table does
    not show. This is evidence only, and harmless.
* **Not an issue.** `verify_certs.py` takes `sigma` from 60-digit `mpmath` rather than from the
  exact integer formula of Lemma 2.3(2). It is safe: `XY != 0` and `X^2 != Y^2` always hold.
  The referee's checker uses the exact formula and agrees on every certificate.

---------------------------------------------------------------------------------------------

## 9. Corrected statement

**Class `P3-G_X(M, C)`.** Genuine residues `r_j(q) = pi_j mod q^(A_q)`, every monomial strengthened.
It contains every actual cluster, for every `X`.

1. **(Proved) Mechanism.**
   * `m_(v,s) = gcd(D_s(A_v), L'_X)`.
   * A tight `(v, s)`, i.e. one with `|D_s(A_v)|` dividing `L'_X`, is equivalent to the one-sided
     constraint `<v, phi - theta> in sigma [0, pi - 2d] + pi Z`, with `sigma` and `sin^2 d` exact
     integer data.
   * Circuits `sum n_i sigma_i v_i in Pi`, `sum n_i (pi - 2 d_i) < pi`, pin `<v_i, phi - theta>` to `pi Z`.
   * The certified closure `Pi*(pi, X)` depends only on the rational prime set and `X`, up to
     conjugation covariance.
   * Relative to `Pi_0 = {0}`, arc sums of zero combinations lie in `(pi/2) Z` (Lemma R.1). So
     level-0 certificates have arc sum exactly `pi/2`.
2. **(Proved) Collapse.** Suppose `2 e_j in Pi*(pi, X)` for every `j`. Then every member of
   `P3-G_(X')`, `X' >= X`, whose prime set is **exactly** `{p_j}` is an actual cluster after
   `pi_j -> i^(t_j) pi_j`. This holds for any `M`, `C`, exponents and rows.
3. **(Proved) Dichotomy**, for squarefree members with `X' >= max(X, 2M)`. The members have
   `>= rank V - k >= (1/2 - o(1)) M - k` columns outside the closure, or `k + 1` independent
   genuine near-relations.
4. **(Exact, independently re-verified)** The six random prime sets of Theorem 4.1 have complete
   closures. Hence `P3-G` restricted to members whose prime set is one of these six sets consists
   of actual clusters. Profiles on proper subsets are **not** covered: 163-subsets of the
   `M = 164, b = 8, K = 4` set leave 35-54 columns free.
5. **(Proved in the model `T`, non-adaptive primes)** Under the assumptions of Theorem 5.3, all
   columns of size `log p <= (1/2 - eps) H_X` are pinned with probability `1 - o(1)`. Hence a
   squarefree non-actual member with `O(M/log M)` near-relations has `W >= (1/4 - eps - o(1)) M H_X`.
   This depends on the model, on non-adaptivity, and on the near-relation hypothesis.
6. **(Proved) Fakes exist** in `P3-G_X` at `W <= (66 + o(1)) M X` (Prop. 6.1), for every `X >= 2M`.
7. **(Heuristic)** Fakes at `W ~ (66 + eps) M H_X` (Heuristic 6.2).
8. **(Open)** Whether `P3-G` has non-actual members at `W = O(M log M)`.
   * Builder-selected prime sets of polynomial size are not covered by 4 or 5.
   * Selected sets have empty weight-`<= 2` closures at sizes where random sets collapse:
     `M = 644` and `M = 1284` at `K = 4`, and `M = 164` down to `K = 2.75`.
   * After pruning, selected sets have empty weight-`<= 3`, `k <= 3` closures on 148-159 of 163
     columns at `K = 3.0`-`3.5`, where random sets collapse (R§4).
   * The full closure (all weights, all `k`) of selected sets is not computed, and each added
     weight restores some pinning.

   So route G has neither been shown to give fakes at `O(M log M)` nor shown not to. The
   heuristic critical size `Theta(H_X)` is plausible. Its measured constants are random-set
   constants. P3-all\* is untouched, and no growth statement for actual circles follows.

---------------------------------------------------------------------------------------------

## 10. Files (`round7/referee_p3G_checks/`)

All scripts run from that directory with `python3 <script> ...`. None of them imports the author's
code.

| file | content | type |
|---|---|---|
| `rcore.py` | Miller-Rabin, Hermite-Serret, `L'_X`, exact tight relations of weight `<= 2`, exact `sigma` and `sin^2 d` | library |
| `rclosure.py` | `k <= 2` closure (signed components with `Z`/`Z2`/`0` status) and an independent checker (sparse echelon lattice over `Z`, orientation from `Im(A^2(-i)^s)`) | library, exact |
| `run_cases.py` -> `cases_A.jsonl`, `cases_B.jsonl` | Theorem 4.1 re-run (all 6 collapse rows and 2 partial rows) | exact |
| `subset_test.py` -> `subset_out.jsonl` | Issue C: subsets of the verified `M = 164, b = 8, K = 4` set | exact |
| `completion_check.py` -> `completion_check_out.txt` | completion sums 15190.1 and 4051.9 | float sum of exact data |
| `lemma52_quant.py` -> `lemma52_quant_out.txt` | Lemma 5.2 by exact enumeration; Lemma R.1 quantisation check | exact / evidence |
| `adaptive.py`, `adaptive_ties.py` -> `adaptive_*.txt`, `adaptive*_primes_*.json` | greedy selection avoiding level-0 certificates (same-`v`, double/single, closed walks) | evidence (sets saved) |
| `eval_set.py` -> `adaptive_eval.jsonl` | `k <= 2` closure on the selected sets. Its `k2w3` fields are **not** meaningful: classes of weight 3 are skipped by the `rclosure` quotient. See `gclosure.py`. | exact for weight `<= 2` |
| `author_code_on_adaptive.txt` | the author's `gtight.closure` on the selected sets (cross-check only) | -- |
| `brute_k3.py` -> `brute_w3.jsonl` | all level-0 certificates with `<= 3` terms, `n_i in {1,2}`, weight `<= 2` or `<= 3` | exact (arc quantised) |
| `gclosure.py` -> `gclosure_w3*.jsonl`, `gclosure_rand375.txt` | closure modulo an arbitrary lattice (canonical cosets), `k <= 2` and `k = 3`, weight `<= 3`; certificates re-checked (50-digit arcs, lattice membership) | exact sub-closure |
| `prune.py` -> `prune_K*.txt`, `*_pruned.json` | delete hitting columns until no level-0 weight-`<= 3` certificate is left, then `gclosure` | evidence (sets saved) |
| `prune2.py` -> `prune_w2_*.txt`, `*_pruned_w2.json` | the same with weight-`<= 2` relations only (tie-accepting sets) | evidence (sets saved) |
| `brute_w2_more.jsonl`, `gclosure_w2_ties.jsonl` | brute force on all selected sets; the general closure on the `M = 164`, `K = 2.75` tie-accepting set (full collapse before pruning) | exact |
| `farkas.py`, `farkas2.py`, `farkas2_debug.py` | LP attempts (infeasible because of exact ties, Lemma R.1) | negative result |
| `notes.md` | scratch notes | -- |
