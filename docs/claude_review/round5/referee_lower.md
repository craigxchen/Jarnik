# Referee report: `round5/lower.md` (the angle-model threshold `W_min(M, C) = Theta(M log M)`)

**Checks.**
* My checks are in `round5/referee_lower_checks/`. Each script runs from that directory with
  `python3 <script>`, and its output is stored next to it as `*_out.txt`.
* They do not import the author's `paley.py` or any other file in `round5/checks/`. The Paley matrix, the
  profiles, the residue design and the Gaussian arithmetic are rebuilt from the text of `lower.md`.
* All mathematical comparisons use exact integers or Python big integers. The one exception is
  `r5_threshold.py`, which evaluates the written union bound in floating point and is labelled as an
  illustration.
* I re-ran the author's `check_fake.py` and `check_canonical.py`. Both outputs are byte-identical to the stored
  ones (`rerun_check_fake.txt`, `rerun_check_canonical.txt`).
* The repository `/home/user/Jarnik` was not modified.

Notation follows `lower.md`. `F(c) = ||H_A^T c||_1`, `2B(c) = ||c^T S||_1`, `r = b(M-1)`, `n = ||c||_1`,
`h = |supp c|`.

---------------------------------------------------------------------------------------------

## 0. Verdict

**The result survives.** On the class `A_X(M, C)` of Definition 1.1,

```text
3 M log M - O(M)  <=  W_min^(2M)(M, C)  <=  (48 + o(1)) M log M        (fixed C, M -> infinity),
```

so `W_min/(M log M)` does not tend to infinity.

I re-derived every step, including the two questions the task singles out:
* **Every input is valid for actual circles.** Prop. 1.2 holds, and its restriction to `v(c) != 0` is
  necessary. I checked (1.1) exactly on 5,523 instances with unit 1 and 22,092 with all four units.
* **The class contains every actual cluster.** After the `two_thirds` normalisation, which does not
  increase `W` or `C`, every cluster is a member.

I found no error that changes a conclusion. The following were re-derived line by line:
* Lemma 2.1, Lemma 2.2 and **Lemma S**;
* Prop. 3.2 and Prop. 3.3;
* Lemma 4.1, Lemma 5.3, Theorem 5.1(a)-(d) and Cor. 5.4;
* Prop. 1.3 and Remark 1.4.

Lemma 5.2 needs one extra hypothesis, which the witnesses satisfy.

**Major (framing): "residue" must mean point-level character residues.** The verdict (§0) and the summary
say that "the local + character + residue + GAP toolkit is sharply capped" and "cannot improve the growth
rate". That is established only for the residue data of Def. 1.1. Those data are point residues, used
only through collisions of zero-sum characters, at odd moduli up to `M^A`.
* Actual circles also satisfy the **prime-level parity law** (Prop. 8.1 of the write-up). I checked it
  exactly on 2,619 instances.
* The write-up's own residue design violates that law. On the `M = 12` instance, 9 of the 11 non-base
  rows violate it.
* For the class with prime-level (factored) residue data, the best rigorous bound known is
  `W_min <= (36 + o(1)) M^2` (`referee_construct.md` §3.4). Whether `W_min/(M log M) -> infinity` there is
  **open**.
* `lower.md` §8 (P3) says this correctly. The verdict and the summary need the qualifier.
* The task's own question defines the class "as in A.3", which is point-level. So the posed question is
  answered correctly. The overreach is in the gloss.

**Minor points.**
1. **Lemma S is correct, but Theorem 5.1(b) does not need it, and it gives a weaker result than material
   already in the notes.** Research-notes item 393 proves two inequalities for *every* capacity-two
   assignment. The write-up cites item 393 only for its classification of exact flats. Together the two
   inequalities give `2B - r >= 0.6 M - 22` for every non-pair ternary character when `b >= 5` and
   `M >= 70`, with no exceptional set (§3.3; verified exactly). Two consequences follow:
   * Theorem 5.1(b) holds for every fixed `b >= 17`, not only `b >= 35`.
   * With `b = 16`, every capacity-two assignment, the canonical one included, gives
     `W_min^(2M) <= (32 + o(1)) M log M`. That beats both 48 and `2b >= 70`.

   The novelty claim "an independent route through Lemma S, which covers every capacity-two assignment"
   should therefore credit item 393. Lemma S remains new as a statement about the pure Paley transform.
2. **Def. 1.1 omits some valid local inputs** (§1.3):
   * collisions against the units `i, -1, -i`;
   * twisted point collisions with `sum c != 0`;
   * 2-adic moduli.

   The witnesses satisfy all three at `X = 2M` (proofs and checks in §1.3). The twisted extension is
   **not** established for Theorem 5.1(c) at `X = M^A` with `A > 1`.
3. **Lemma 5.2 is stated too generally.** If a deleted column's prime satisfies `p_j <= X`, the restricted
   tuple needs (R) at `p_j` for every character, and the original member does not supply it. Add the
   hypothesis "every `p_j > X`", which all witnesses satisfy.
4. **The theorems are purely asymptotic, with very large thresholds.** The union bound as written in the
   proof of 5.1(a) drops below 1 only at `log M ~ 85` for `C = 1` (`log M ~ 209` for `C = 0.1`). Part (b)
   needs `tau > 68`. This is legitimate for an `o(1)` statement, but it should be said. No finite instance
   of Theorem 5.1 is certified by any check.
5. **Cosmetic.**
   * In the proof of 5.1(a), the ternary prefactor `8pi/C` should be `4pi/C`, and `12pi/C` overcounts the
     pairs. Both are harmless upper bounds.
   * The author's `unionbound.py` is correctly labelled as an illustration.

The corrected statement is in §8.

---------------------------------------------------------------------------------------------

## 1. The class: validity, containment, and what it leaves out

### 1.1 Prop. 1.2 (validity): correct

I re-derived the proof.

* `U = prod_(c>0) z^c` and `V = prod_(c<0) z^(-c)` both have norm `N^(n/2)`.
* `G = gcd(U, V)` has `|G| = N^(n/4) e^(-w(v)/2)`. In `U` the `pi_j`-exponent exceeds that of `V` by
  `v_j`, and the `conj(pi_j)`-exponent falls short by `v_j`. So `U/G ~ A_v` and `V/G ~ conj(A_v)`.
* A collision `prod rho^c = 1 mod q^a` gives `q^a | U - V = G D`. Since `q` is prime to `Norm G`, this gives
  `q^a | D`. Distinct `q` are coprime, and `(1+i) | D` because all odd-norm Gaussian integers are
  `1 mod (1+i)`. So `|D| >= sqrt2 m_c`, and `D != 0` because `A_v` is not associate to `conj(A_v)`. This step
  needs `v != 0`.
* `|U - V| = 2 N^(n/4) |sin(c.delta/2)|`, because `U/V = e^(i c.delta)` for `sum c = 0`.
* Same-level pairs at `p_j`: `z_x - z_y = D_j (r_x - r_y)` and `g_xy = D_j h_xy` with `p_j` not dividing
  `Norm h_xy`. These give `p_j^a | c_xy`, which is `two_thirds` Lemma 3.3 (S|).

**Exact check** (`r3_residues.py` (A)). I used four circles:
* `N = 5^2*13*17*29`, `5*13*17*29*37`, `13^2*17*41` and `5^3*13^2`;
* random subsets of 2 to 8 lattice points;
* random zero-sum `c`;
* `X = 40..60`;
* level residues at `p_j <= X` for same-level pairs.

I tested the integer form `|U - uV|^2 prod p_j^|v_j| >= 2 m_u^2 N^(n/2)`. Results:
* it holds in all 5,523 instances with `v(c) != 0` and `u = 1`, and in all 22,092 instances with any unit;
* 19,004 level-residue factors entered;
* among 199 instances with `v(c) = 0`, the inequality fails 12 times.

So the restriction "`v(c) != 0`" in Def. 1.1(e) is necessary, and the author has it. (`referee_construct.md`
§2.1 found the same defect in the sibling write-up, whose definition lacked it.)

### 1.2 Containment and Prop. 1.3: correct

* Normalisation (`two_thirds` Lemma 1.1) does not increase `W` or `C`, and keeps `min_x a_xj = 0` and
  `max_x a_xj = e_j`, so every prime of `N'` varies. Residues `z_x mod q^(A_q)` have norm `N`. Level residues
  have norm `N/p_j^(e_j)`.
* **Prop. 1.3.** `two_thirds` (I read §§1-5) uses the configuration in exactly two places:
  * the summed pair inequality (2.5). In the class it follows from (1.1) for `c = e_x - e_y` and
    `|delta_x - delta_y| <= Delta`, which give `log(sqrt2 m_xy) <= log C - W/4 + d_xy/2`. Summing with the cut
    identity reproduces (2.5) exactly;
  * the pigeonhole (3.6). The residue classes there are `z mod q^a` with unit 1 only. Their counts depend
    only on the norm, and the collisions it uses have `q^a < (5/4) M < 2M <= X`. Same-level data at
    `p <= X` are covered by Def. 1.1(d).

  The rest of `two_thirds` is combinatorics and prime sums. So `3 M log M <= W + (14 + 4 log^+ C) M` holds
  on `A_X` for `X >= 2M`.
* **Remark 1.4.** Correct. For pairs, `arcsin(Q_M e^(-d/2)) <= Delta/2` gives
  `d >= W/2 + 2 log Q_M - 2 log(C/2)`. Summing against `sum d <= M^2 W/4` gives
  `W >= 4(M-1)(psi(2M) - log(C/2))`. For the 10-point cluster (`C = 8.24`) the right side is
  `36 * (19.27 - 1.42) = 643`, while `W = 20.9`. The per-assignment reading of A.3 (III) is the right one,
  and it agrees with `referee_construct.md` §1.4.

### 1.3 Valid local inputs that Def. 1.1 omits (minor; the witnesses survive at `X = 2M`)

All three of the following are valid for actual circles. Each enlarges the class, and so strengthens the
barrier when the witnesses still pass.

**(a) Collisions with the units `i, -1, -i`.** If `prod z^c = u mod q^a`, then
`|sin(c.delta/2 - arg(u)/2)| >= 2^(-1/2) m_u e^(-w/2)`. I checked this exactly in the 22,092 instances
above.

The witnesses satisfy it, because (5.1) bounds `dist(c.delta/2, (pi/4)Z)` from below by `m^gen`, where
`m^gen >= m_u` for every `u`:
* `|sin(t - k pi/4)| >= sin dist(t, (pi/4)Z)`;
* and `m^gen e^(-w/2) >= 2^(-1/2) m_u e^(-w/2)`.

Def. 1.1(e) should simply quantify over all four units.

**(b) Twisted point collisions** (`sum c = 2k != 0`). The input is
`prod rho^c = u N^k mod q^a`, which implies `|sin(<v,phi> - psi_u)| >= D' e^(-w(v)/2)/sqrt2`, where
`v = sum c_x a_x - k e`. This is `referee_construct.md` §2.3. For the canonical profile:
* `v = S^T c/2` lies outside `L`, by rank `M`. So `<v, phi>` is Haar-uniform on the realisation torus.
* `||S^T c||_1 >= max(2n, (b-2)(M-n))`, so `w(v) >= tau max(n, 3(M-n))` for `b = 8`. The bound `2n` comes
  from flipped columns with their *distinct* siblings. The bound `(b-2)(M-n)` comes from at least `b - 2`
  unflipped copies per label, together with `F >= M - n`, which follows from
  `sum T^2 = M ||c||_2^2 - (sum c)^2`. Checked in `r2_profiles.py` (twisted), with margins 20, 36 and 84 at
  `q = 11, 19, 43`.
* At `X = 2M`, every `D' <= e^(psi(2M)) <= e^(2.08 M)`. The added bad proportion is at most
  `sum_n (4M)^n e^(2.08M) e^(-(tau/2) max(n, 3(M-n)))`. With `tau = 6 log M` this is
  `exp(-(1.5 - o(1)) M log M)`.

So Theorem 5.1(a) survives (b).

At `X = M^A` with `A > 1`, however, `D'` can only be bounded by `e^((1+o(1)) M^A)`, which beats
`e^(-w/2) = e^(-O(M log M))`. I see no design that controls twisted collisions there. The extension of 5.1(c)
to twisted inputs is therefore **not established**.

**(c) 2-adic moduli.** A.3's `Q_M` contains powers of 2, and actual circles satisfy `2^a | D` when
`prod z^c = 1 mod 2^a`. Extend the design by `rho_x(2^a) = beta_2 varpi_x/conj(varpi_x)`. Write
`Omega = P + iQ`; exactly one of `P, Q` is even, and `P +- Q` are odd. Then:
* a collision at `2^a` against `+-1` means `2^(a-1)` divides that even coordinate;
* `+-i` give no 2-adic collision.

Hence `2^(a2-1) m^gen | Im(Omega^4)/4`, and the bound `Y^(2n)/4` is unchanged. Checked on all 1,499
characters of `r3_residues.py` (B).

---------------------------------------------------------------------------------------------

## 2. Lemma 2.1, Lemma 2.2 and Lemma S: correct

**Lemma 2.1.**
* (1)-(3) are Parseval and inversion `c = (1/M) H T`.
* (4): with `u_a = |T_a|/n`, `sum u(1-u) = (F - M m_2/n)/n`, and both one-sided sums are `<= 2 sum u(1-u)`.
* (5) and (6): `T(c - sH_a) = T - sM e_a` and `T(g) = (M/2)(s e_a + s' e_b)`, the residual is a nonzero
  integral zero-sum vector, and `F >= M` applies to it.

All correct.

**Lemma 2.2.** Correct.
* After dividing by `q^t`, `T_a = -G(a) mod q` off `supp c'`, where `G(Y) = sum c_x (Y-x)^d - c_inf`
  (Euler's criterion, including `y = 0`).
* The coefficients of `Y^1..Y^d` vanish only if the moments `sum c_x x^m` vanish for `m < d`. With `s' <= d`,
  the Vandermonde matrix then forces `c' = 0 mod q`.
* So at most `s' + d` labels have `T_a = 0`.

**Lemma S.** Correct. I re-derived each step with `theta = 1/16`.
* *Step 1.* `|F - kappa h| < 2 theta M` gives `kappa h` in `((7/8)M, (19/16)M)`.
* *Step 2* (`kappa >= 3`). Here `h < 19M/48 < (q-1)/2`. Lemma 2.2, Parseval (`#{|T_a| = h} <= M/h`) and the
  defect identity give `F - M >= (2(h-2)/h)(M/2 - h - M/h)`. Its values are:
  * `M/4 - 4` at `h = 4`;
  * `(4/3)(M/3 - 6)` at `h = 6`;
  * by concavity, `>= (3/2) min(3M/8 - 8, 5M/48 - 48/19)` on `[8, 19M/48]`.

  All three are `>= M/16` exactly when `M >= 40.4`, so `q >= 43` suffices.
* *Step 3* (`kappa = 1`). `1 - u_a < 2 theta M/h` gives `|T_a| > h - 2 theta M > (1 - 4theta) M`. Column
  peeling then gives `F > (3/2) M`.
* *Step 4* (`kappa = 2`). `|T_a|, |T_b| > (1 - 6theta) M/2`, and
  `phi(t) = t - |t - M/2| >= min(2t - M/2, M/2)`. So `F > (5/4) M`.

**Exact checks** (`r1_lemmaS.py`):
* (A) `K = DH = I + J` with `J` skew, at `q = 7, 11, 19, 43, 47`.
* (B) Lemma 2.2 on 22,604 vectors. At `q = 7` I tested all `c` with finite support `<= 3`, entries in
  `[-2,2]` and any `c_inf`. For `q = 11..59` I used random vectors and adversarial ones: quadratic-residue
  supports, arithmetic-progression supports, and `q`-divisible entries. The minimum slack is 1 at every `q`,
  as the author reports.
* (C) Lemma 2.1 (1)-(6) on 12,000 vectors.
* (D) Lemma S at `q = 43, 47, 59`. I used:
  * all zero-sum ternary vectors of support 4 (407,253, 583,740 and 1,462,905 vectors);
  * 6,000 random ternary vectors of every even support;
  * 300 greedy descents of `F` over non-exceptional ternary vectors, started from perturbed columns,
    perturbed half-sums, and sign patterns of 3- and 4-column combinations.

  The minimum of `(F - M)/M` is 0.364, 0.333 and 0.400, always at support 4. This matches the author's
  census, and the proved `1/16` is far from tight.

---------------------------------------------------------------------------------------------

## 3. Profiles: Props. 3.2 and 3.3 correct; item 393 already suffices for every assignment

### 3.1 Identity (3.1), twins, Prop. 3.2: correct

* (3.1): the `b - #flips` unflipped copies of label `a` contribute `|T_a|` each, and the flipped copy of row
  `x` contributes `|T_alpha(x) - 2 c_x H(x, alpha(x))|`.
* Twins: `S_f(x) - S_s(x) = -2 H(x, alpha(x)) e_x`. In the 0/1 normalisation this forces `lambda` to be
  integral. (In the `+-1` normalisation of walsh-extraction Lemma F.1 the same fact reads "`2 lambda`
  integral".) So Theorem F is vacuous for every `m`, as stated.
* Prop. 3.2, skew step. With `T = d - Jd`, a lowering finite flip with `|T_x| = h` needs `d_x J_xy d_y = -1`
  for all `y` in the support. Two such rows contradict `J^T = -J`. The infinity flip at `a*` adds at most
  one more. The remaining lowering flips sit at distinct labels with `2 <= |T| <= h - 2`, each contributing
  `>= 2(h-2)` to the defect.
* For `b >= 5` and `h >= 4`, `b(1 - 2/h) - 2 > 0`, which gives `B - r/2 >= b/2 + h - 4`.
* Amplitude: `F >= MA` and the flips cost at most `2n <= 2MA`.

**Exact check** (`r2_profiles.py`), canonical assignment at `q = 11, 19, 43, 47` and `b = 5, 8`:
* (3.1) agrees with the direct `c^T S` on 400 random vectors per case;
* twins, rank `M`, nonconstant columns and distinct rows all hold;
* tested characters: all pairs, all support-4 ternary vectors (`q <= 47`), 3,000 random ternary vectors, all
  columns and half-sums, and 1,500 amplitude-`>= 2` vectors;
* every bound of Prop. 3.2 holds;
* the pair and ternary bounds are attained (margin 0). The ternary minimum is attained at columns:
  `B - r/2 = b/2 + M - 4`.

### 3.2 Prop. 3.3: correct

I re-derived each row from (3.1):
* pairs, `F = M` and two supported flips: `b - 4`;
* columns, at most two rows with label `a`: `b + 2M - 8`;
* half-sums, at most four rows with labels in `{a, b}`: `b + M - 16`;
* other ternary vectors, by Lemma S;
* amplitude: `bM(A/2 - 1) >= 0` gives `>= AM(b-4)/2 + b`.

**Exact check**: random capacity-two assignments, three each at `(q, b) = (43, 8), (47, 8), (59, 8), (43, 35)`.
* The pair, column and half-sum bounds are attained exactly.
* The intermediate inequality `2B - r >= b(F - M + 1) - 2n` is attained.
* Every other row holds.

### 3.3 Item 393 already gives a uniform `Omega(M)` bound (minor; affects (b), constants and novelty)

`paley_primefield_polynomial_character_gap.md` (item 393) proves, for **every** one-flip capacity-two
assignment with `b >= 5`, two bounds for a ternary zero-sum `c` of support `h`:

```text
(i)  h < q:   F >= h(2M - 2h)/(h + 2)          [#{T_a != 0} >= M/2 - h, and (t-2)(h-t) >= 0]
              => 2B - r >= b(F - M + 1) - 2h;
(ii) h >= 8:  2B - r >= b + 2h - 8M/h          [k <= 2 #{T_a != 0} (capacity two),
                                                 #{T_a != 0} <= M/h + M delta/2 + M delta/h].
```

**Consequence (referee; proof).** Let `M >= 70` and `b >= 5`. Then every non-pair zero-sum ternary `c`
satisfies `2B - r >= 0.6 M - 22`.

*Proof.*
* *Case `4 <= h <= 0.3M`* (so `h < (q-1)/2`). By (i), `F - M >= g(h) := (M(h-2) - 2h^2)/(h+2)`. In
  `u = h + 2`, `g = M + 8 - 2u - (4M+8)/u` is concave. So `g >= min(g(4), g(0.3M)) = M/3 - 16/3` for
  `M >= 70`. Hence `2B - r >= b(M/3 - 13/3) - 0.6M >= 1.06 M - 22`.
* *Case `h >= 0.3M`* (so `h >= 8`). By (ii), `2B - r >= 5 + 0.6 M - 80/3`. QED

Columns and half-sums are included: no exceptional set is needed, because the capacity bound (ii) handles
large supports directly. Compare Lemma S's row of Prop. 3.3, `(b/16 - 2) M + b`, which is positive only for
`b > 32`.

**Exact check** (`r6_item393.py`): random capacity-two assignments at `q = 43, 47, 59, 83` and `b = 5, 8`,
including all support-4 vectors at `q = 43`, random ternary vectors of every support, all columns and all
half-sums.
* (i) and (ii) hold. (ii) is attained with slack 0.
* `min_h max((i), (ii))` equals `0.34 M` to `0.66 M` at these orders, and `0.714 M` at `M = 10^5` for `b = 5`.
* The observed `min (2B - r)/M` over non-pair ternary characters is 0.75 to 0.91.

**Effect on Theorem 5.1** (referee sketch; same residue design, `X = 2M`, `Y = 3 M log M`).
* Pairs need `(b-4) tau/4 > log(M^2 Y^4)`.
* Every other ternary character has slack `>= tau (0.3 M - 11) - 1/2`, which beats `3^M e^(2.08 M)` once
  `tau > 22`.
* Amplitude `>= 2` is unchanged.

With `tau = 2 log M + 2 loglog M + O(1)`, this gives:
* **Theorem 5.1(b) for every fixed `b >= 17`**, with `W <= (2b + o(1)) M log M`;
* at `b = 16`, the pair sum is `O((log M)^(-2))`, so `W_min^(2M)(M, C) <= (32 + o(1)) M log M` for every
  capacity-two assignment, the canonical one included.

So Cor. 5.4's `1/24` improves to `1/16`. The sibling `construct.md` (refereed) has 27, giving `2/27`.

**Novelty.** Lemma S is correct, and new as an `l^1`-stability statement about the pure Paley transform.
Item 393's bound (i) degenerates near `h = M/2`, where Lemma S still gives `F >= (17/16) M`. But Lemma S is
not "what (P1) needed", and the claim "an independent route through Lemma S, which covers every capacity-two
assignment" should credit item 393, which already covers every capacity-two assignment. `outside.md` §2.2
("Nothing handles `||c||_1` in `[M^(1/2), 2M]`") and its referee overlooked item 393's capacity bound (ii) as
well as item 394.

---------------------------------------------------------------------------------------------

## 4. Lemma 4.1 (residue design): correct

* `prod rho^c = Omega/conj(Omega)`.
* Collision against `eta` means `q^a` divides `Q`, `P`, `P - Q` or `P + Q`.
* At most one unit per prime: `q | eta - eta'` is impossible for odd `q`.
* The four groups are coprime, so `m^gen | PQ(P^2 - Q^2) = Im(Omega^4)/4`.
* Nonvanishing: `v_varpi_x(Omega) = c_x^+ != c_x^- = v_varpi_x(conj Omega)`.

All correct.

**Exact check** (`r3_residues.py` (B)): `M = 44`, `X = 88`, the 44 primes `= 1 mod 4` in `(88, 613]`, and
`N` the product of 30 actual primes `= 1 mod 4` above `10^6`.
* Norms hold at all 23 odd primes.
* On 1,499 characters (pairs, random sums, amplitudes up to 3): `m_c <= m^gen | Im(Omega^4)/4 <= Y^(2n)/4`,
  with at most one unit per prime.
* `max log m^gen / log(Y^(2n)/4) = 0.917`.

---------------------------------------------------------------------------------------------

## 5. Theorem 5.1, Lemmas 5.2-5.3, Cor. 5.4

* **Lemma 5.3.** Correct. Conditioning on one coordinate gives density `<= 2/Delta`. At most
  `2 n Delta/pi + 2` windows of length `<= pi eps` are hit. When `eps > sin(pi/8)` the stated bound already
  exceeds 1, so (5.1) is never claimed in an infeasible case.
* **(5.1) implies (G) on `L` and (R).** Correct, since `m^gen >= 1` and `(pi/4) Z` contains both target
  sets.
* **Theorem 5.1(a).** Correct. I checked the following:
  * the weights: there are more than `r^2` primes `= 1 mod 4` in `[Z, 2Z]`, a window of log-width `1/r`
    exists, and `W <= r tau + 1`;
  * `p_j > Y > X`, so `N` is prime to all `q <= X` and no level data arise;
  * the three class sums: pairs `81 (12 pi e^(1/4)/C)/(log M)^2`, ternary with `x = 18/log M`, and
    amplitude `e^(-(3 - o(1)) M A log M)`;
  * completion by `outside.md` A.1(ii), with `Pi - 1 = O(r/sqrt Z)`;
  * the restriction to all `M`.

  *Threshold* (`r5_threshold.py`, illustration). The written bound first drops below 1 at `log M = 85` for
  `C = 1` (209 for `C = 0.1`, 630 for `C = 0.01`), where `W/(M log M) ~ 50.5`. This is fine for an
  asymptotic theorem, but it means that no computed instance supports 5.1(a).
* **Theorem 5.1(b)** is correct as stated, given Lemma S. (`tau > 68` is needed.) See §3.3 for the stronger
  form.
* **Theorem 5.1(c)** is correct. Pairs need `tau > (2 + 4A) log M + O(1)`, ternary vectors need
  `tau/2 > (1 + 2A) log M + O(1)`, and amplitudes have margin `log M`. The twisted caveat of §1.3(b)
  applies.
* **Theorem 5.1(d)** is correct as an attribution. Item 394's (4)-(6) impose only the axis gaps; the bad
  probability is `< 0.224`. Adding the two diagonal events at `eps/sqrt2` at most doubles it, to
  `0.448 < 1`.
* **Lemma 5.2** (minor gap). Suppose a deleted column has prime `p_j <= X`. Then `p_j` does not divide
  `N_new`, so Def. 1.1(d)-(e) require residues of norm `N_new` at `p_j` and (R) at `p_j` for *all*
  characters. The original member's level data only constrain same-level pairs. Fix: add the hypothesis
  "all `p_j > X`", which the witnesses satisfy.

  The rest of the proof is right:
  * deleted columns are constant on the kept rows, so they add a common constant to (1.1) and `v(c)`
    vanishes on them;
  * residues are rescaled by `beta` with `Norm beta = N_new/N`, which leaves `prod rho^c` unchanged because
    `sum c = 0`.
* **Cor. 5.4.** Correct for inputs (I), (II), (III) as defined (zero-sum characters, unit 1, odd moduli
  `<= M^A`) and (IV). The witnesses satisfy (IV) in both senses, through the Fubini argument of
  `referee_outside.md` §3.1. Also `loglog R = (1 + o(1)) log M`, because `log R >= (3/2 - o(1)) M log M` by
  Prop. 1.3.

**Answers to (P1)/(P2) and the phase lemma.**
* (P1) is answered affirmatively. Items 392, 394 and Theorem A.1 already did so for the canonical
  assignment, and item 393 suffices for every capacity-two assignment (§3.3).
* (P2) is false as a universal statement over fully flipped capacity-two Paley profiles at
  `W = O(M log M)`: admissible `delta` exist with `k = 0`.
* The angle-model analogue of the walsh-extraction §10 phase lemma fails for every `K >= 48` (and for every
  `K >= 32` by §3.3). The lemma is quantified over all `K`.

---------------------------------------------------------------------------------------------

## 6. Scope: prime-level residue data (major, framing)

**The law** (`r4_parity.py` (A)). For odd `q` not dividing `N`, and actual points
`z_x = eps_x prod pi_j^(a_xj) conj(pi_j)^(e_j - a_xj)`, the class of `z_x/z_0` in `U_q/U_q^2` equals

```text
class(eps_x/eps_0) + sum_j (a_xj - a_0j) nu_j(q),        nu_j(q) = [p_j is a non-residue mod q].
```

I checked this exactly on all primitive points of three circles (`5*13*17*29*37`, `5^2*13*41`,
`13*17^2*53`) at every odd `q < 60` not dividing `N`: 2,619 instances, no exception. It is Prop. 8.1 of the
write-up, and §3.1 of `referee_construct.md`.

**The design violates it** (`r4_parity.py` (B)). On the canonical `q = 11`, `b = 8` profile with actual
primes, at `q = 7, 17, 23` (`+-1 mod 8`, where the units are squares), only 2 of the 11 non-base rows are
consistent. Six different Legendre patterns would be needed for the auxiliary primes `ell_x`.

At `M = 12` one could pick the `ell_x` by pattern. Asymptotically, however, the patterns live modulo
`prod q ~ e^((1+o(1)) M)`. The write-up's §8 states this obstruction correctly.

**Conclusion.**
* The class with factored residue data (`A_2` of `referee_construct.md`) contains every actual cluster.
* Its best rigorous upper bound is `(36 + o(1)) M^2`, from A.3' on the `b = 8` profile.
* So the statement "the local + character + residue + GAP toolkit cannot improve the growth rate" is
  established only for the point-level class `A_X`, with `X <= M^A`.
* Whether `W_min^(A_2)/(M log M) -> infinity` is open (the write-up's (P3)). A positive answer would still be
  a general growth improvement.

**Possible extension (not checked; not claimed).** The bound `m^gen | Im(Omega^4)/4` does not depend on `X`.
Using two disjoint families of auxiliary primes, one for `q` outside `{ell}` and one for `q` in `{ell}`,
might extend the point-level barrier to residue data at **all** odd moduli. I have not written this out.

---------------------------------------------------------------------------------------------

## 7. Attribution

* §7.1 is right that items 392 and 394, together with `outside.md` A.1, already gave canonical fakes at
  `O(M log M)` without residues. Item 394 imposes only the axis gaps, as the write-up says. `outside.md` does
  not list items 392-394 among the notes it searched, and its referee missed them.
* The write-up in turn cites item 393 only for the exact-flat classification. It misses item 393's theorem
  "every capacity-two assignment ... all saturated characters above primitive half-height" and the two
  inequalities of §3.3, which make the every-assignment route independent of Lemma S.
* New relative to the notes:
  * Lemma S as a transform statement;
  * the residue design (Lemma 4.1);
  * the class definition with validity (Prop. 1.2) and its 2/3-consistency (Prop. 1.3);
  * Remark 1.4;
  * the combination with residues;
  * Prop. 8.1, also found independently in `referee_construct.md`.

---------------------------------------------------------------------------------------------

## 8. Corrected statement

1. **Class.** Use `A_X(M, C)` as in Def. 1.1. Preferably quantify (e) over all four units, and optionally add
   2-adic moduli and twisted collisions (`sum c != 0`), all of which are valid for actual circles.
   * Every normalised actual cluster is a member.
   * The restriction `v(c) != 0` is necessary.
   * For `X >= 2M` every member satisfies `3 M log M <= W + (14 + 4 log^+ C) M`.
   * The literal A.3 class forces `W >= (8 - o(1)) M^2` and contains no small actual cluster.
2. **Lemma S** holds as stated (`q >= 43`, `theta = 1/16`). The census suggests `theta ~ 1/3`.
3. **Theorem 5.1.** For fixed `C`, as `M -> infinity`:
   * `W_min^(2M)(M, C) <= (48 + o(1)) M log M` (as stated);
   * using item 393's capacity bound: `<= (32 + o(1)) M log M`, for every capacity-two assignment with
     `b = 16`, and (b) holds for every `b >= 17`;
   * `W_min^(M^A) <= (16 + 32A + o(1)) M log M`.

   The witnesses also satisfy the all-unit, 2-adic and twisted inputs at `X = 2M`. The twisted part is open
   at `X = M^A`, `A > 1`. The thresholds are astronomically large (`log M ~ 85` at `C = 1` for the written
   proof of (a)).
4. **Lemma 5.2** with the hypothesis "every prime of `N` exceeds `X`".
5. **Barrier.** No argument using only
   * exponent identities,
   * all Lemma 0.1 gaps,
   * point-level residue collisions of zero-sum characters at odd moduli `<= M^A`, and
   * Haar-generic angle properties

   proves `M = o(log R/loglog R)`. The best constant of this class lies in `[2/(16 + 32A), 2/3]`. At `X = 2M`
   that is `[1/16, 2/3]` by §3.3, and `[2/27, 2/3]` with `construct.md`.
6. **Scope.** "Residue" means point-level data. Prime-level (factored) residue data are valid for actual
   circles and violated by these fakes. For that class only `3 M log M - O(M) <= W_min <= (36 + o(1)) M^2` is
   known, and whether `W_min/(M log M) -> infinity` there is open (P3). The summary sentence "the local +
   character + residue + GAP toolkit cannot improve the growth rate" must carry this qualifier.

---------------------------------------------------------------------------------------------

## 9. Files (`round5/referee_lower_checks/`)

| script | what it checks | type |
|---|---|---|
| `rl_common.py` | Paley-I matrix from the text of §2 (Hadamard asserted), transform, exceptional set by definition | helper |
| `r1_lemmaS.py` | skew structure; Lemma 2.2 (22,604 vectors, adversarial supports); Lemma 2.1 (12,000); Lemma S search (exhaustive support 4 at `q = 43, 47, 59`, random, 300 greedy descents) | exact; (D) is evidence |
| `r2_profiles.py` | (3.1), twins, rank; Prop. 3.2 (canonical, `b = 5, 8`, exhaustive pairs/support 4); Prop. 3.3 (random capacity-two, `b = 8, 35`); twisted bound `max(2n, (b-2)(M-n))` | exact |
| `r3_residues.py` | Prop. 1.2 on actual points (all units, level residues, `v = 0` failures); Lemma 4.1 design at `M = 44`, `X = 88`; 2-adic extension | exact (big integers) |
| `r4_parity.py` | parity law on actual points (2,619); the `M = 12` design violates it (2/11 rows consistent) | exact |
| `r5_threshold.py` | where the written union bound of 5.1(a) drops below 1 | floating point, illustration |
| `r6_item393.py` | item 393's bounds (i), (ii) on random capacity-two assignments; uniform `0.6M`-type bound; comparison with Lemma S's row | exact |
| `rerun_check_fake.txt`, `rerun_check_canonical.txt` | re-runs of the author's scripts: identical outputs | — |
