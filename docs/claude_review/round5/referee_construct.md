# Referee report: `round5/construct.md` (cheap angle models, `W_min = Theta(M log M)`)

Referee checks are in `round5/referee_construct_checks/`. Each script runs from that directory with
`python3 <script>`, and its output is stored next to it as `*_out.txt`. Every check is independent of the
construct's `lib.py`, `gauss5.py` and `instance.py`, and uses integer, `Fraction` or 50-digit `Decimal`
arithmetic. The repository `/home/user/Jarnik` was not modified.

## 0. Verdict

**The result survives as precisely stated.** The following are correct:

* Theorem 1 (character margins);
* Lemma 2.2 (decomposition);
* Lemma 4.1 (residue design);
* Theorem 2 (fake circles in `A_1(Q_0)` and `A_0`);
* Prop 1.1 (the literal A.3 class);
* Remark 2.3 (no skew pairing for Walsh cores of order `>= 16`);
* Corollaries 6.1 and 6.2.

I re-derived every step and found no error that affects a conclusion. Every constant in the
union bound checks exactly. An independent instance at the smallest admissible order (`M = 8`, at
`P = P_*` exactly) gives a failure probability `<= 0.447`.

So, over the class `A_1` (point-level residue data with axiom (R) at odd prime powers `<= 2M`):

```text
(3 - o(1)) M log M  <=  W_min(M, C)  <=  (27 + o(1)) M log M + 9 M log(1/C).
```

`W_min/(M log M)` does **not** tend to infinity. No argument whose inputs are exponent identities, all
Lemma 0.1 gaps, (R) at moduli `<= 2M` and Haar-generic angle properties can prove
`M <= (2/27 - eps) log R/loglog R`.

Defects found:

1. **Lemma 1.2 is false as stated (minor).** The claim is that (R) holds on actual clusters "for every
   nonzero zero-sum `c`". It fails when `v(c) = sum c_x a_x = 0`, i.e. on multiplicative rectangles.
   (R) must be restricted to `v(c) != 0`. The fakes are unaffected, since they have rank `M`, and for
   `C <= sqrt 2` actual clusters never have `v(c) = 0`. (§2.1)
2. **Scope of "residue" (major, framing).** Every actual cluster obeys a *prime-level parity law* for
   its point residues. I proved it for twin profiles and checked it exactly on 211,136 actual
   instances. The construct's own `M = 12` instance violates it in all 11 rows, for every choice of
   units.
   * So the barrier covers point-level residue-collision arguments (axiom (R)), not "local residue
     information" in general.
   * The construct discloses this in §7.1. But the headline ("the local + character + residue class
     ... cannot give a growth improvement") needs that qualifier.
   * For the prime-level class, `W_min/(M log M) -> infinity` is not excluded. That is now the precise
     remaining question. (§3)
3. Minor points (§2.2-2.7):
   * the `3 M log M` lower bound on `A_1` needs the level-wise cofactor data (S|), which the class
     definition mentions only in words; without it only `(2 - o(1)) M log M` follows;
   * the `Q_0 = M^K` constant `9(2+K)` is valid only for `K >= 1`;
   * 2-adic moduli, which are part of A.3's `Q_M`, are omitted, though the design extends verbatim;
   * item 394 imposes only the axis gaps, not the diagonal gap of Lemma 0.1;
   * the summary omits `q >= 7`;
   * `check_axiomR.py` silently skips `v = 0`, exactly where Lemma 1.2 fails.

The corrected statement is in §6.

---------------------------------------------------------------------------------------------

## 1. Re-derivation of the main line

### 1.1 Lemma 2.2 and Theorem 1: correct

**Lemma 2.2.** Let `y = S^T c` and `T_a = (H^T c)_a`.
* Unflipped columns give `b sum_{a>=1}(|T_a| - 1) = b(F - M + 1)`.
* The flipped column of row `x` gives `T_{a(x)} - 2 c_x h_x`.
* The map `a` hits every label once and `a*` twice, so `sum_x |T_{a(x)}| = F + |T_{a*}|`.

Summing gives `G = (b+1)(F - M) + b + |T_{a*}| + sum_x Delta_x`. With equal weights,
`V_c - W/2 = (1/2) sum_j w_j (|y_j| - 1) = (tau/2) G`. Correct.

**Theorem 1.**

* `F >= M` follows from `sum_a T_a^2 = c^T H H^T c = M ||c||_2^2`, `T_0 = 0`, `|T_a| <= n` and
  `||c||_2^2 >= n`.
* **(i)** For a pair, rows `x, y` disagree on exactly `M/2` labels, so `F = M` and
  `G = b + |T_{a*}| + Delta_x + Delta_y >= b - 4`.
* **(iii)** `c = H T/M` gives `|c_x| <= F/M`, so `F >= KM`. Then `(K-1)M >= (K-1)n/K >= n/2` and
  `sum Delta >= -2n`. Correct.
* **(ii)** Each step re-derives:
  * the defect identity `sum_a |T_a|(n - |T_a|) = n(F - M)`;
  * `T_a` is even;
  * `|L'| <= n(F-M)/(2(n-2))`;
  * the big labels give `H[X, A] = c|_X sigma^T`, which has rank one;
  * two paired rows `x != y` of `X` with `pi(x), pi(y) in A` would give a rank-one principal `2x2`
    minor, contradicting Definition 2.1;
  * `Delta_x = 2` whenever `T_{a(x)} = 0`;
  * `N_X <= 2 + |L'|`. I checked the case where `x_inf` maps into `L'` rather than `A`: then the
    rows in `A` number at most 1 (paired), and the rows in `L'` at most `|L'| + 1`;
  * `(b + 1 - 2n/(n-2)) >= 0` for `n >= 4`, `b >= 3`.

  Correct.

**Pair sharpness (Remark 3.3).** I proved that the pair minimum is exactly `b - 2`:
* two beneficial paired flips on a pair would violate the skew minor;
* a non-beneficial pair flip always has `Delta = +2`;
* `Delta_{x_inf} = -2` forces `|T_{a*}| = 2`.

**Paley-I skew pairing.** `H_n(i,i) H_n(j,j) - H_n(i,j) H_n(j,i) = 1 - chi(-1) = 2`. Correct.

More generally, a normalised skew-Hadamard matrix `K = I + J` has the pairing `x -> x`, since the
minor is `K(x,0)K(y,0)(1 + J_xy^2) != 0`.

**Exact check (R1, `r1_theorem1.py`).** The independent build checks Theorem 1 on:
* Paley-I for `q = 7, 11, 19, 23`;
* a **non-Paley, non-Sylvester** skew-Hadamard matrix of order 16, from the doubling
  `[[A, A], [A - 2I, -A + 2I]]` of Paley-8. This tests the general skew-pairing statement, not just
  Paley.

| core | `b` | exhaustive range | min slack over the bound |
|---|---|---|---|
| Paley 8, 12 | 3, 5, 8 | all `+-1` characters, all `n` | 2 (pairs) |
| Paley 20, 24 | 3, 5, 8 | all `+-1`, `n <= 6` (1,345,960 at `M = 24`) | 2 (pairs) |
| doubled skew 16 | 3, 5, 8 | all `+-1`, `n <= 8` (450,450 at `n = 8`) | 2 (pairs) |
| Paley 8, 12, 20; doubled 16 | 5 | amplitude `>= 2`, entries in `[-2,2]`, `n <= 8` or `6` | 38-110 |
| Paley 44, 68 | 3, 8 | 120 hill-climbs from random, column, half-sum and amplitude-2 starts | 2 (pairs) |

The pair minimum `G = b - 2` holds in every case. No violation was found.

### 1.2 Lemma 4.1: correct, and it extends to `q = 2`

* `prod rho_x^{c_x} = Lambda_c/conj(Lambda_c)` since `sum c = 0`.
* `conj(Lambda_c)` is a unit modulo `q^a`, because `l_x > Q_0 >= q`.
* `Lambda_c != u conj(Lambda_c)`, by the `omega_x`-valuation, using distinct split `l_x`.
* Both terms have odd norm, so `(1+i)` divides the difference.

Hence `sqrt2 D_u <= |Lambda_c - u conj Lambda_c| <= 2 L^(n/2)`.

**R4 (`r4_residue_design.py`)** checks this in the strong form, as Gaussian **divisibility**
`(1+i) D_u(c) | Lambda_c - u conj Lambda_c`, plus the two norm inequalities, on:
* `M = 8, 12` with `n <= 6`;
* `M = 20` with `n <= 4`;
* 434,448 (character, unit) instances in all.

It also verifies the 2-adic extension. Put `rho_x^(2^a) = s_2 omega_x/conj omega_x`, with `s_2` of
norm `N mod 2^a`, which exists because `N = 1 mod 4`. Then `2^(a_2) | Lambda - u conj Lambda`, so the
same (5.1) covers collision strengthening by powers of 2 (max `a_2 = 5` seen).

### 1.3 Theorem 2: correct

1. **Step 1.** Zero coordinates lose at most `eta/2` each, so `V_c - W/2 >= (tau/2)G - 1/8`.
2. **Step 2.** (R) depends only on `c.delta` and `arg u`, because `psi_u` absorbs `c.k`. Since
   `psi_u in (pi/4)Z`, the bound `|sin(X - psi)| >= sin dist(X, (pi/4)Z)` holds, so (5.1) implies (R)
   for all four units. Since `L >= 1`, it also implies both Lemma 0.1 inequalities on `L`. Saturation
   holds by the twins: a flipped column and its unflipped sibling differ by `-h_x e_x` in `a`.
3. **Step 3.** The density is `<= 2/Delta`, from conditioning on one coordinate with `c_x != 0`.
   The number of windows is `<= 2 Delta n/pi + 2 <= n + 2`.
4. **Step 4.** The count `(2M)^n` is correct. The constant chain in R6 (`r6_union_bound.py`, 50-digit):

   ```text
   arcsin t / t <= 1.0367476            (t <= 5^(-1/2))
   4 * 1.0368 * e^(1/16) = 4.41467      <= 4.415
   sum_{k>=2} (2k+2) x^k <= 6.2065 x^2  (x <= 1/40)
   4.415 (4 + 6.3/40)/40 = 0.45888      < 0.46
   ```

5. **Step 5.** This step completes the realisation. `outside.md` Theorem A.1(ii) (refereed) applies
   because `p_j >= 21 r^2`. The referee's Fubini argument gives positivity in the full torus, since the
   residue data do not depend on `delta`.
6. **Size.** The constant is `r log P_* / (M log M) -> 27` for `Q_0 = 2M` (L is about `2M log M`), and
   18 in `A_0`. Measured values are 50.2, 41.8, 37.5 and 35.1 at `M = 12, 104, 1004, 10004` (R6),
   which matches the construct's `rates.py`.

**Independent instance (R6).** `M = 8` (`q = 7`), `C = 1`, `Q_0 = 16`, `L = 89`.
* `P = P_* = 911,360` exactly, with 64 primes `= 1 mod 4` in `[P, P e^(1/256)]`.
* `1.3201(Pi - 1) = 0.189`.
* Union bound `<= 0.4467`.
* `W = 878.3`, which is `52.8 M log M`.

The construct's five instances (`instance_*_out.txt`) use the same computation. I read
`instance.py`, and its union bound applies the proved margins correctly.

**Consistency with Theorem C.** This is a real test of the result, not a formality.
* Theorem C (`walsh-extraction.md`) is proved from angle-model inputs: characters, aligned flats and
  `k`-residue pigeonhole.
* So if Theorem 1 held for every flip assignment, Theorem 2 would build Walsh-core fakes that
  contradict Theorem C.
* Theorem 1 uses the skew pairing essentially, and Walsh cores of order `>= 16` have none (R3,
  below). So there is no conflict.

### 1.4 Prop 1.1 and Remark 2.3: correct

**Prop 1.1.**
* (i) Property 3 on pairs forces `d_xy - W/2 >= 2 log Q_M + 2 log(2/C)`. The cut identity gives
  `sum_{x<y}(d_xy - W/2) <= WM/4`. Together, `W >= 4(M-1)(log Q_M + log(2/C)) = (8 - o(1)) M^2`.
* (ii) For `2+i` and `1+2i` on `N = 5`, property 3 needs `arcsin(12/sqrt5)`.

So the literal uniform-`Q_M` reading is not a valid modelling class. `A_1` (with residue data) is the
right formalisation of "residue strengthening as in A.3". I agree with this correction of the task's
framing.

**Remark 2.3 (R3, `r3_walsh_skew.py`).** Independent backtracking confirms:
* pairings exist for `t = 2, 3`;
* there is **no** injective `pibar : F_2^4 \ 0 -> F_2^4` (value 0 allowed) with
  `(x+y).(pibar(x) + pibar(y)) = 1` (exhaustive, 17,041 nodes).

The restriction argument is correct: `x + y in U`, and a collision `pibar(x) = pibar(y)` would give 0.
All normalisations of Sylvester are again Walsh cores, so the statement is complete. As a by-product,
the doubled order-16 matrix of R1 is not equivalent to Sylvester.

---------------------------------------------------------------------------------------------

## 2. Defects

### 2.1 Lemma 1.2 fails when `v(c) = 0` (minor; the fix costs nothing)

The proof's step 6 uses "`A` is a conjugate-primitive nonunit", which needs `v(c) != 0`. If
`v(c) = 0`, then `A = 1`, `Z_+ = E Z_-`, and for `u = E`:
* `D_u(c)` is the full product of the odd prime powers `<= Q_0`;
* but `EA - u conj A = 0`.

So (R) demands `0 >= D/sqrt2`.

**Counterexample (R2, `r2_axiomR_v0.py`, exact).**
* On `x^2 + y^2 = 65`, the points `1+8i, 8+i, 4+7i, 7+4i` have gcd 1 and arc constant 3.754.
* `(1+8i)(8+i) = 65 i = (4+7i)(7+4i)`, the exponent rows satisfy `a_1 + a_2 = a_3 + a_4`, and
  `c = (1,1,-1,-1)` has `v(c) = 0`.
* With `Q_0 = 8`, `D_1(c) = 21`; (R) would need `0 >= 14.85`.
* On `N = 32045`, 12 such violations occur among the 4-subsets of an 8-point arc.
* Every instance with `v(c) != 0` (23,792) satisfies (R).

The construct's `check_axiomR.py` skips `v = 0` (`if not any(v): continue`), so its 6,336 instances
never test the failing case.

**Fix.** Impose (R) only for `v(c) != 0`. With that change:
* Lemma 1.2 holds;
* the fakes still satisfy it, since `rank S = M` gives `v(c) != 0` for every `c != 0`;
* `A_U` is contained in `A_1`.

For `C <= sqrt 2`, actual clusters have affinely independent rows (`linear_allocation_affine_rigidity.md`),
so `v(c) != 0` automatically. Rectangles need `C >= 2 sqrt2`.

### 2.2 Lower bound on `A_1` and the (S|) data (minor)

* Corollary 6.1's `(3 - o(1)) M log M` is the `two_thirds` bound. Its split-prime part for small
  `p | N` needs the same-level cofactor collisions (S|), which feed the chain lemma.
* Definition 1.2 excludes `q | N` and says only that "one adds the level-wise cofactor data".
  Without (S|), a model can put every small split prime into `N` with balanced cuts. Then only the
  inert part survives: `W >= (2 - o(1)) M log M`.
* The class should state an (S|)-type axiom formally: an (R)-type inequality for pairs on a common
  level, using cofactor residues modulo `p^a`.
* The fakes satisfy it vacuously, because all their primes exceed `Q_0`. The upper bound is
  unaffected.

### 2.3 `Q_0 = M^K` (minor)

`L` is at least the `M`-th prime `= 1 mod 4`, which is about `2M log M`, whatever `Q_0` is. So the
constant is `9(2 + max(K, 1))`, not `9(2 + K)`, for `K < 1`. The intended range (the referee's
`M^k/k!`, `k >= 1`) is fine.

### 2.4 2-adic moduli (minor)

* A.3's `Q_M = lcm(1..2M)` includes powers of 2. Class `A_1` uses odd moduli only.
* For actual circles the 2-adic version of (R) is valid: `K` has odd norm, so `2^a | EA - u conj A`.
* As R4 shows, Lemma 4.1 covers it with no change to (5.1).
* The corrected class should include `q = 2`.

### 2.5 Attribution of item 394 (minor)

I confirm that item 392 (`paley_canonical_all_character_half_height.md`) is the skew-tournament
argument of Theorem 1: at most one beneficial full-magnitude flip, plus the Parseval defect (its
(2)-(5)).

Item 394 (`paley_compatible_phase_coordinate_gap.md`) is not quite as described:
* It proves simultaneous **axis** gaps `min(|Re G_c|, |Im G_c|) >= 1`, not the diagonal gap
  `|X +- Y| >= 1` of Lemma 0.1.
* So "item 394 + Theorem A.1 = fake circles" needs the diagonal gaps added. Its union bound absorbs
  them with a constant factor, so the `A_0` fakes "were already in the notes" up to this routine
  addition.

The note `lower.md` does state `(48 + o(1)) M log M` (Lemma S route), as cited.

### 2.6 Statement hygiene (minor)

* The summary's Theorem 2 omits `q >= 7`, which the window lemma needs for `r >= 34`.
* `instance.py` asserts `D_u <= 2 L^(n/2)`, which is weaker than the `sqrt2 L^(n/2)` the proof uses.
  The proof's bound is correct (R4 checks `2D^2 <= |diff|^2`).
* The construct correctly says that the `M = 44` `nmax = 4` run and the `M = 24` amplitude-2 run are
  unused.

---------------------------------------------------------------------------------------------

## 3. The questions the task asked

### 3.1 Is every input valid for actual circles?

**Yes, after §2.1.** The inputs are:
* the exponent profile and exponent identities;
* Lemma 0.1 gaps (`outside.md`, refereed);
* (R) for `v(c) != 0` (Lemma 1.2: proof re-derived, and R2 re-verifies the integer form on 23,792
  instances including all four units);
* the arc condition.

Theorem A.1 and the Fubini step are statements about the model, not inputs.

### 3.2 Does the class contain all actual clusters?

**Yes**, with (R) restricted to `v(c) != 0` and (S|) data added (§2.2).
* Normalisation (`two_thirds` Lemma 1.1) does not increase `W`.
* Residues of actual points lie in the conic fibres.
* Units `eps_x = i^(-k_x)` match `psi_u`.

### 3.3 But the class is strictly larger than "local residue information": the prime-level parity law

The issue is in the other direction: the class contains fakes that violate a valid *local* constraint.

**Lemma (parity law; proved here).** Let `z_x = eps_x prod_j pi_j^(a_xj) conj(pi_j)^(e_j - a_xj)` be
actual points, and `q` an odd prime not dividing `N`. Let `K` be the norm-one subgroup of
`(Z[i]/q^a)^*`. It is cyclic of order `(q - chi(q)) q^(a-1)`, so `K/K^2 = Z/2`. Then

```text
class(z_x/z_0)  =  class(eps_x/eps_0) * prod_j chi_q(p_j)^(a_xj - a_0j)      in K/K^2,
class(i) = (2/q),   class(-1) = +1.
```

*Proof.* Write `r_j = pi_j mod q^a`. Then `z_x/z_0 = (eps_x/eps_0) prod_j t_j^(a_xj - a_0j)` with
`t_j = r_j/conj(r_j) in K`.

* For inert `q`, `t_j = r_j^2/p_j`. Since `r_j^(q+1) = p_j`, we get
  `t_j^((q+1)/2) = p_j^((1-q)/2) = chi_q(p_j)`.
* For split `q`, `t_j = (alpha/beta, beta/alpha)` with `alpha beta = p_j`. Its class is
  `(alpha^2/p_j)^((q-1)/2) = chi_q(p_j)`.
* `i^(|K|/2) = (2/q)` by `|K|/2 mod 4`.

QED

**Converse for twin profiles (proved here).** Point residue data `rho_x` (norm `N`) are induced by
*some* prime residues `r_j` (norm `p_j`) and units iff the parity law holds.

* Write `r_j = r_j^0 kappa_j` with `kappa_j in K`. Then
  `rho_x / rho_x^0 = (unit) prod_j kappa_j^(S_xj)`.
* The twins (`kappa_f = kappa`, `kappa_s = kappa^(-1)`) put all of `(K^2)^M` in the image.
* Modulo `K^2`, every `kappa_j` contributes the same class to every row, because `S_xj` is odd. So
  the image is `{y : all y_x in one K^2-coset}`.

**Exact checks.**
* *R5 (`r5_parity_actual.py`).* The parity law holds on all 211,136 (pair, `q <= 59`) instances on
  three actual circles, with `eps_x` computed from the factorisation.
* *R4(c).* For the construct's own `M = 12` instance (`p_j = 10000121..10003529`,
  `l_x = 29..137`, design `rho_x = s omega_x/conj(omega_x)`), take the parity vectors over
  `q in {3,5,7,11,13,17,19,23}`.
  * **All 11 rows violate the law for every choice of units.**
  * The design's ratio has class `chi_q(l_x l_0)`. That class bears no relation to
    `prod chi_q(p_j)^(a_xj - a_0j)`.

**Consequence.**
* The barrier excludes arguments that use point residues only through collisions, i.e. axiom (R).
* It does not exclude arguments that use the factorisation of point residues through prime residues.
  That input is local (modulo `q <= 2M`) and valid for every actual circle, and the construct's fakes
  fail it.
* The construct lists this in §7.1 ("first arithmetic input beyond the class"). The summary's
  "local + character + residue method class cannot give a growth improvement" must therefore read
  "... point-level residue-collision ...".

**What is known for the prime-level class.**
* *Sketch, not written out or claimed.* Choose `l_x = 1 mod 4` with prescribed Legendre symbols
  modulo every odd `q <= 2M`, by CRT and Linnik-Xylouris over distinct residue classes. By the
  converse above, the Lemma 4.1 design is then prime-level consistent, with `log L = O(M)`. This
  gives `W = O(M^2)`.
* *Heuristic.* Taking `log L` about `pi(2M) log 2` gives about `12.5 M^2/log M`, matching
  construct §7.1's pseudosquare route. The construct's random-assignment heuristic is
  `O(M log^2 M/loglog M)`.
* So for the prime-level class `W_min` lies in `[3 M log M - O(M), O(M^2)]` (upper end sketched), and
  **whether `W_min/(M log M) -> infinity` there is open**.
* Because the parity law is valid for actual circles, a proof that it forces
  `W_min/(M log M) -> infinity` would be a general growth improvement. Nothing here excludes that.
* Forced pigeonhole collisions do not increase under the parity split: two cosets of `K^2`, each of
  size `|K|/2`, give the same `M^2/(2|K|)` leading term. So any gain would have to come from
  constraining the *designs*, not from `two_thirds`-type counting.

---------------------------------------------------------------------------------------------

## 4. Novelty and attribution

* Theorem 1 is item 392 in rank-one-rectangle form for a general skew pairing. The construct says so
  correctly. The generalisation beyond Paley is real: R1 checks a non-Paley skew core.
* `walsh-extraction.md` Prop P(2) already quotes the item-392 margin, *individually*. The new step,
  in item 394 and here, is simultaneity plus completion.
* New relative to the notes and `outside.md`:
  * the residue version (Lemma 4.1 and Theorem 2 in `A_1`);
  * the constants 27 and 18;
  * Remark 2.3;
  * Prop 1.1.
* `lower.md` reaches `48` by an independent route, as cited.

---------------------------------------------------------------------------------------------

## 5. Files (`round5/referee_construct_checks/`)

| script | checks | type |
|---|---|---|
| `r1_theorem1.py` | Theorem 1 (i)-(iii), pair minimum `b-2`; Paley 8/12/20/24 and a doubled (non-Paley, non-Sylvester) skew core of order 16, exhaustive; hill-climbs at `M = 44, 68` | exact integers |
| `r2_axiomR_v0.py` | Lemma 1.2 counterexample at `v(c) = 0` (`N = 65`, `N = 32045`); (R) on all `v != 0` instances | exact integers |
| `r3_walsh_skew.py` | Remark 2.3: pairings for `t = 2, 3`; no relaxed injection for `t = 4` | exhaustive |
| `r4_residue_design.py` | Lemma 4.1 as Gaussian divisibility (434,448 instances); 2-adic extension; parity test of the `M = 12` instance | exact integers |
| `r5_parity_actual.py` | prime-level parity law on actual points (211,136 instances) | exact integers |
| `r6_union_bound.py` | Theorem 2 constant chain; independent `M = 8` instance at `P = P_*`; growth of `r log P_*/(M log M)` | 50-digit Decimal; last item floating |

---------------------------------------------------------------------------------------------

## 6. Corrected statement

1. **Class.** Let `A_1(Q_0)` be the angle models (Lemma 0.1 for all `v != 0`) with point residue data
   at all prime powers `<= Q_0` coprime to `N`, `q = 2` included, satisfying:
   * (R) for every zero-sum `c` with `v(c) != 0` and every unit;
   * the analogous same-level cofactor axiom (S|) at small `p | N`.

   Every normalised actual cluster lies in `A_1(Q_0)` (Lemma 1.2 with `v(c) != 0`). The literal A.3
   class (uniform `Q_M`) is not a modelling class (Prop 1.1).
2. **Theorem 1** (any skew-paired fully flipped profile, `b >= 3`) holds as stated.
3. **Theorem 2** holds as stated, for `q >= 7`. The constant is `9(2 + max(K,1))` for `Q_0 = M^K`.
   The design also covers 2-adic moduli.
4. **Corollaries 6.1 and 6.2** hold for `A_1(2M)`:
   `(3 - o(1)) M log M <= W_min <= (27 + o(1)) M log M + 9M log^+(1/C)` along `M = q + 1`, and
   `W_min = O_C(M log M)` for all `M`. No argument using only exponent identities, all monomial gaps,
   (R) at moduli `<= 2M` and Haar-generic angle properties proves
   `M <= (2/27 - eps) log R/loglog R`.
5. **Scope.** "Residue" means point-level collision data. Actual circles also satisfy the prime-level
   parity law (§3.3), which the fakes violate.
   * For the class with prime-level residue data, the construct's method gives `W = O(M^2)` (sketch)
     and heuristically `o(M^2)`.
   * Whether `W_min/(M log M) -> infinity` there is **open**. A proof would be a general growth
     improvement.
   * The statement "the local + character + residue method class cannot give a growth improvement" is
     established only for point-level residue collisions.
