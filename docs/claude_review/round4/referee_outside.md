# Referee report: `round4/outside.md` (GAP/Fourier, energy, incidences, equidistribution; fake-circle barrier)

Referee checks are in `round4/referee_outside_checks/`. Each script runs from that directory with
`python3 <script>`, and its output is stored next to it as `*_out.txt`. The repository
`/home/user/Jarnik` was not modified.

## 0. Verdict

**The theorems survive. The headline framing does not.**

* Theorem A.1 (completion), Lemma A.6 (balancing), Lemma A.7 (character margins), Theorem A.3
  (fake circles), Theorem A.3' and Corollary A.4 as formally stated with classes (I)-(IV) are
  correct. I re-derived every step and found no error that affects a conclusion.
* The side propositions A.4, A.5, B.1(i), B.1(iii), B.2, C.1, C.2 (near-direction case), D.1 and
  D.2 are correct under their stated hypotheses.
* Defects found:
  * one constant is too small by a factor of about 5;
  * one rounding is false in the fourth decimal;
  * the remark "all chord lengths are distinct" (C.3) and the unqualified "multiplicatively Sidon" in
    the summary are false for `C >= 2 sqrt 2`, with the record 10-point cluster as a counterexample;
  * the general-arc bound of C.2 omits an `O(C^2)` term;
  * the parenthetical description of class (III) is inaccurate.
* **Major (framing):** the summary says "All four directions reduce to barred classes". For (a) and
  (b) that is not established.
  * Theorem A.1 reduces (a) and (b) to *simultaneous* character admissibility.
  * The refereed character barrier (Prop P of `walsh-extraction.md`) checks characters only
    *individually*. Theorem A.1 itself shows that individual checks give only necessary conditions.
  * The proved ceiling for the reduced class is `c (log R/loglog R)^(1/2)` (Theorem A.3). That is
    far below `log R/loglog R`.
  * So the write-up does **not** show that (a)/(b) cannot yield a growth improvement
    `M = o(log R/loglog R)`. Its own problem (P2) is exactly such a route. The body of the
    write-up says this correctly in §2.2: "proved ceiling between `c (log R/loglog R)^(1/2)` and
    `(2/3) log R/loglog R`". The verdict and the summary overstate it.

The corrected statement is in §7.

---------------------------------------------------------------------------------------------

## 1. Theorem A.1 (completion): correct

**(i) Determinacy.** From (1.1), two realisations differ by an element of `U = A~^(-1)(R 1)`.
Then

```text
U^perp = A~^T(1^perp) = {2 sum lambda_x a_x : sum lambda = 0} = span_R{a_x - a_y},
```

using `(T^(-1)L')^perp = T^*(L'^perp)` in finite dimension. So `<v,.>` is constant on the coset iff
`v` lies in `U^perp cap Z^r = L`. The value (1.2) follows because `theta_0` and `<e,phi>` cancel
when `sum lambda = 0`. Correct.

**(ii) Completion.**

* `U` is a rational subspace, so `U cap 2 pi Z^r` is a lattice in `U` and `T_U` is a compact
  torus.
* For `v` not in `L`, `u -> <v,u> mod 2 pi` is a well-defined nontrivial continuous character of
  `T_U`. It is surjective, so it pushes Haar measure to Haar measure.
* The union of the two bad events of Lemma 0.1 lies in `{dist(<v,phi>, (pi/4)Z) < arcsin e^(-w/2)}`.
  That event has probability `<= (8/pi) arcsin(e^(-w/2))`. Overlap only helps the union bound.
* `sum_(v in Z^r) e^(-w(v)/2) = prod_j (1 + 2/(sqrt p_j - 1)) = Pi`.
* Pairing `v` with `-v` gives a bad proportion `<= 1.32 (Pi - 1)`.
* The sufficient condition `p_j >= 21 r^2` gives `Pi - 1 <= 0.558`, with the maximum at `r = 1`. So
  `1.33 (Pi - 1) <= 0.742`.

Correct.

**Cosmetic rounding.** `arcsin y <= 1.0367 y` is false at `y = 5^(-1/2)`. The sharp ratio is
`1.036745` (`r2_constants_out.txt` (e)), and `(8/pi) * 1.036745 = 2.64003 > 2.64`. Use `1.0368`,
`2.641` and `1.3201`. The theorem's constant `1.33` absorbs this.

**Scope.** (ii) needs `Pi < 1.75`, so every varying prime must be large. All known extremal
clusters use small primes (5, 13, 17, ...). There, monomials outside `L` have small weight, and
their gaps may cut the realisation torus nontrivially. "Monomial gaps carry exactly the information
of the characters" is therefore proved only in the large-prime regime. The text says so. The
summary's "(a) ... the reduction is now proved" omits the restriction. The barrier (Theorem A.3)
does not need it, since the fakes are built with large primes.

## 2. Theorem A.3 (fake circles): correct, one constant to repair

### 2.1 Construction and Lemma A.6

* *Flips.* `2(M-1)` flipped columns over `M` rows, each row receiving 1 or 2 flips. This is possible
  for `M >= 2`. Designated flips are injective, and siblings are a bijection `F_a -> S_a`.
* *Normalisation.* Every column is nonconstant, since a balanced `+-1` column with one entry
  flipped is still nonconstant for `M >= 4`. So Lemma 1.1(d) of `two_thirds` holds with `e = 1`.
* *Lemma A.6 (greedy pigeonhole).* Correct: `binom(n,k) = (n/k) binom(n-1,k-1)`, and the condition
  `n > k^2 (m-1) J` is exactly what the greedy step needs.

**Constant error (minor).** The text says "about `150 M^5` primes ... we may take
`P ~ 400 M^5 log M`". This is too small.

* Primes `= 1 mod 4` in `[P, 2P]` number about `P/(2 log P)`.
* At `P = 400 M^5 log M` that is about `0.21` of the requirement for every `M` tested.
* An exact sieve confirms the deficit: at `M = 8` there are `7.8e5` available against `3.32e6`
  needed, and at `M = 4` there are `2.1e4` against `6.8e4`.
* The PNT estimate needs `P ~ 2000 M^5 log M` (`r2_constants_out.txt` (a), (b)).

This changes only the `O(1)` in `tau = 5 log M + loglog M + O(1)`. It does not change
`(30 + o(1)) M^2 log M`, `(15 + o(1))`, `p_j >= 21 r^2`, or the step-3 bound. A larger `P` only
helps, because `P^0.43 >= 13 M^2.15` remains true.

### 2.2 Lemma A.7

* **(i)** By Parseval, `sum_(a>=0) T_a^2 = M ||c||_2^2`. Also `T_0 = 0`, `|T_a| <= n` and
  `||c||_2^2 >= n`. Correct.
* **(ii)** `S_(f(x)) - S_(s(x)) = -2 H(x,a) e_x`. So `v = sum lambda_x a_x` integral forces
  `lambda` integral, and `v = (1/2) lambda^T S` because `sum lambda = 0`. `c^T S_j` is even because
  `sum c = 0`.
  * This is the same mechanism as Lemma F.1 / item 392 of `walsh-extraction.md` (there in the `+-1`
    normalisation, where it gives `2 lambda` integral). The attribution should be added.
* **(iii)** Re-derived line by line:
  * Twins give `w_f |T - 2c_x h| + w_s |T| >= 2 tau |c_x|`.
  * `O`-columns give `(1/2) Omega'_a |T_a|`.
  * The parity of `T_a` and the window give
    `sum_a Omega'_a (|T_a| - 1) >= Omega'_0 (F - M + 1) - (M-1)/M`.
  * Then `(1/2)(b-4) tau = (3M-2) tau`, and the `F_a cup S_a` loss is `<= 2(M-1)(tau + log 2)`.
  * Together, `V_c - W/2 >= tau(n + M) - 1/2 - 2(M-1) log 2 >= tau(n + 0.86M)` for `tau >= 10`.

  Correct.

**Independent checks** (`r1_lemmaA7.py`, `r5_aligned_flats.py`). My own build uses *random* flip
rows, designated flips and siblings, and *actual* primes `= 1 mod 4` in `[P, 2P]`. The `O_a`
log-sums are balanced into windows `< 1/M`, and logs are re-verified at 50 digits whenever the float
slack is `< 1e-6`.

* (i), (ii) (twins and rank) and the **weighted** (iii) hold on:
  * Sylvester 8 and 16;
  * Paley 12 and 20;
  * all characters with `||c||_1 <= 6` (`M = 8`) or `<= 4`;
  * columns and two-column half-sums;
  * 1,000 to 3,000 random characters;
  * adversarial hill-climbing on the slack.
* The minimum weighted slacks are 11.7, 20.6, 26.0 and 38.0, always at pairs.
* Every aligned Walsh flat was tested, on every subspace, coset and nonzero functional, for
  Sylvester 8, 16 and 32 with 5 random flip assignments each: 385, 3,825 and 47,585 flats. All are
  flat (`F = M`), and all satisfy the equal-weight (iii), with minimum slack 10, 24 and 54.
* The author's `check_fake_barrier_out.txt` reports 389,603 characters per `b`, which matches the
  claimed "about 389,600".

### 2.3 Probability, union bound, completion

* *Step 1.* `X_c = (1/2) c^T delta` has density `<= 2/(Delta ||c||_inf)` and range length
  `Delta n/2`. Correct.
* *Step 2.* The count of points of `(pi/4)Z` within `eta` of the range is `<= 2 Delta n/pi + 2` for
  `eta <= pi/8`. The algebra to `1.04 Q_M P^(-(n+0.86M)/2)(2.55 n + 8/C)` is correct.
* *Step 3.* The count `#{c : ||c||_1 = n} <= 3^M 2^(n-1)` is correct. The other inputs are also
  correct:
  * `Q_M <= e^(2.08 M)`: checked exactly for `M <= 2000`, with maximum `log Q_M / M = 2.0615`;
    in general it follows from `psi(x) < 1.03883 x`;
  * the final sum is `< 1` for `M >= M_0(C)`; it is `<= 5e-8` already at `M = 4`, `C = 1e-3`, with
    the corrected `P` (`r2_constants_out.txt` (c)).
* *Step 4.* Rank `M` (from the twins) gives realisability for every `delta`. Theorem A.1(ii) gives
  property 2 for `v` outside `L`. Pair gaps give distinct points. Correct.
* *Step 5.* `log R = W/2 <= (15+o(1)) M^2 log M`. Since `x/log x` is increasing,
  `log R/loglog R <= (7.5 + o(1)) M^2`. Correct.

**Theorem A.3'.** With `P >= 35 r^2 Q_M^2`:

* `Q_M e^(-w/2) <= 1/(sqrt 35 r)`;
* `Pi - 1 <= 2.2 r/sqrt P`;
* `1.32 Q_M (Pi - 1) <= 0.49`;
* `tau <= (4.16 + o(1)) M` and `log R <= (12.5 + o(1)) M^3`.

Correct. As an upper bound, `log Q_M ~ 2M` gives `4 + o(1)`.

## 3. Corollary A.4: correct as formally stated; two remarks on the classes

### 3.1 Class (IV) can be strengthened to "Haar-null in the full angle torus"

The summary says "any property of the angles that holds off a Haar-null set" without naming the
torus. The text proves the per-realisation-torus version. The full-torus version is stronger,
because a null set of `(R/2 pi Z)^r` may contain entire realisation tori. It also holds, because
Theorem A.3 randomises `delta`.

*Proof.*

* Let `D' = {delta : delta_1 = 0}` (a copy of `R^(M-1)`).
* Let `Psi : D' x T_U -> (R/2 pi Z)^r` be `(delta', u) -> S^+ delta' + u`, where `S S^+ = I`.
* The differential of `Psi` is injective, since `S^+ d delta' in U` forces `d delta' in R 1 cap D'`,
  which is `0`. Dimensions match, because `dim U = r - M + 1`. So `Psi` is a local diffeomorphism,
  and preimages of null sets are null.
* The good `delta` have positive measure (union bound `< 1`), and characters depend only on
  `delta mod R 1`. So their projection to `D'` has positive measure.
* Each fibre has a good proportion `>= 1 - 1.33(Pi - 1)`.
* By Fubini the good set in `D' x T_U` has positive measure, and it is not contained in
  `Psi^(-1)(Z)` for any null `Z`. QED

So Corollary A.4 holds with (IV) read in either sense.

### 3.2 Class (III): the gloss "the range in which pigeonhole forces collisions" is accurate only for single-point residues

* For points and pairs it is right. `two_thirds` uses classes `(q - chi(q)) q^(a-1) < M`, hence
  `q^a < 2M`. The collision product of a pair divides `Q_M`. The fake gives `|c_ij| >= 2 Q_M`,
  which dominates `sqrt 2 prod q^(a_q)`. I re-checked this against `two_thirds/proof.md`, Lemmas
  3.1-3.4 and Cor. 3.6.
* Pigeonhole on **products** of `k` points (`binom(M+k-1,k)` products in
  `(q - chi(q)) q^(a-1)` classes) forces collisions of `2k`-term characters modulo prime powers up
  to about `M^k/k!`. For example, 4-term collisions are forced modulo primes up to about `M^2/2`.
  Those moduli are outside (III).
* The fake still satisfies any *single* such strengthening numerically, since `Q_M ~ e^(2M)`
  exceeds every forced modulus with `k = O(M/log M)`. It does not satisfy "any residue assignment"
  at those moduli.
* Optional strengthening. The margin `P^(-n/2)` easily absorbs a per-character factor
  `M^(2||c||_1)`: the step-3 sum with `Q_M M^(2n)` in place of `Q_M` is still `< 1e-5` for
  `M >= 4`, `C >= 1e-3` (`r2_constants_out.txt` (f)). So (III) can be enlarged to "per-character
  collision moduli `<= Q_M M^(2||c||_1)`".
* I grepped the notes ("pigeonhole", "product collision") and found no certificate that uses
  product-residue pigeonhole. The notes' product collisions are exact equalities, which lie in
  class (I). So the claim "covers every residue-pigeonhole certificate in the notes" stands. The
  parenthetical should be corrected.

### 3.3 Consistency with refereed results (no contradiction found)

* **`two_thirds`.** `3M log M <= W + K_C M` holds, with `W ~ 30 M^2 log M`.
* **Theorems B and B''.** They need pure cores or `<= M/(40 log M)` modified rows. Every fake row
  is flipped.
* **Theorem C and item 335.** These give a `log*` growth only.
* **`linear_support_character_height_obstruction.md`.** Repeated-Walsh profiles need `b >> sqrt M`
  and `M = O((log N/loglog N)^(2/3))`. The fakes have `b = 6M` and exponent `1/2`.
* **Sub-clusters.** The twins survive restriction to any row subset, so saturation, and every
  character gap, persist for sub-clusters.
* **`round4/flipped.md`, Theorem L.** Conditionally on Lang-Waldschmidt, it excludes every
  Hadamard-core profile with flip mass `F < (1/2 - delta) W` when the angles are those of actual
  Gaussian primes. The fakes have `F/W ~ 1/(3M)`.
  * So a two-logarithm transcendence input would exclude genuine versions of the fake profile.
  * This is consistent with §6 item 2, and it sharpens it: the fake profile class is exactly one
    transcendence statement away from exclusion.

### 3.4 What the class covers is an interpretation, not a theorem

The formal content is "(I)-(IV) cannot prove `M = o((log R/loglog R)^(1/2))`". The list
"Fourier/Riesz-product, Bohr-set, inverse LO, Delsarte, multiplicative energy" is a classification
of methods by their arithmetic input.

* An argument in these families that uses integrality of additive expressions lies outside the
  class (§6 item 1). An example is the power sums `sum_x z_x^k in Z[i]` or `z_a + z_b`.
* "every Fourier/Riesz-product ... argument" should read "every such argument whose arithmetic input
  is (I)-(IV)".
* Likewise, "This is the best inequality obtained from the GAP structure alone" (after Prop A.4) is
  unproved. The proved statements are A.4 (upper) and A.3 (ceiling).

## 4. The headline "all four directions reduce to barred classes" (major, framing)

For (a) and (b), Theorem A.1 lands on **simultaneous** rational-character admissibility. This
class is not barred against a growth improvement.

* The refereed barrier for character arguments is Prop P (walsh-extraction §7.1, and its referee
  §4, item 4). It establishes that fully flipped Paley profiles pass each character inequality
  individually, at `W = Theta(M log M)`.
* By Theorem A.1, individual passes are only necessary for realisability.
* The only proved ceiling for the simultaneous class is Theorem A.3:
  `M >= c (log R/loglog R)^(1/2)`.
* Between that and `(2/3) log R/loglog R`, nothing excludes a growth improvement from (a)/(b).
* The write-up's own (P2), non-simultaneous admissibility of fully flipped Paley profiles at
  `W = O(M log M)`, is precisely such a route. It would be a purely real argument in direction (a).

The correct statement is the one in §2.2 of the write-up: "a proved ceiling between
`c (log R/loglog R)^(1/2)` and `(2/3) log R/loglog R`". It is not "(a), (b) reduce to barred
classes". For (c) and (d) the "orthogonal or local" assessment is reasonable. It is an
interpretation, not a theorem.

## 5. Side propositions

* **Prop A.4.** Correct. Its content is (9) of `linear_allocation_affine_rigidity.md`
  (`M <= ceil(C/sqrt2)(r+1)`, proved there from the pair bound (5) alone) plus `r log(4r/e) <= W`.
  It is not new, and it is weaker than `two_thirds`. The write-up does not claim novelty, which is
  right.
* **Prop A.5.** Correct and trivial.
* **Prop B.1.** (i) and (iii) are correct; I re-derived `z_a z_b/(z_c z_d) = u A_v/conj A_v`.
  (ii) is correctly restricted to `C < 2 sqrt 2` in the text, but the summary says "The cluster is
  Sidon both additively and multiplicatively" without the restriction.
  * **Counterexample:** the record 10-point cluster (`N = 1176852625`, `C = 8.2393`). It is symmetric
    under `z -> -i conj z` (reflection in the diagonal, which it contains). So
    `z_0 z_9 = z_1 z_8 = ... = -iN`, verified exactly (`r4_misc_out.txt`).
  * The barrier conclusion for (b) survives, because exact product relations are exponent identities
    in class (I). "Sum-product inequalities are vacuous" is then unproved for `C >= 2 sqrt 2`.
* **Remark C.3.** "Within a cluster all chord lengths are distinct, by the Sidon property" is
  **false** for `C >= 2 sqrt 2`. The 10-point cluster has only 25 distinct squared chord lengths
  among its 45 chords. Equal chords are equivalent to `z_a z_d = z_b z_c`, which is the
  multiplicative Sidon property, valid only for `C < 2 sqrt 2`.
* **Prop B.2.** The identity is correct; I re-derived it and checked it exactly for `M <= 60`. The
  "average slack `W/4`" holds for **balanced** binary profiles. For unbalanced columns the pair
  average is `< W/2`, and "the pair average is at most `W/2`" gives no lower bound. The summary
  omits "balanced".
* **Prop C.1.** Correct (exponent count `3e - 2(max - min)` and the circumradius formula). I
  re-checked both identities exactly on all 120 triples of the 10-point cluster with independent
  Gaussian gcd code. The 10-point cluster does contain a diagonal direction, so its bound of 98 is
  right.
* **Prop C.2.** The proof is correct, including the at-most-one-point-per-line step for `eta > 0`
  and the range bound `D(eta + D/2)`. The general-arc bound `2 C^(3/2) R^(1/4) + 2` ignores that
  Dirichlet's parameter `Q` is an integer. Taking `Q = ceil(sqrt(2 sqrt R / C))` gives
  `2 C^(3/2) R^(1/4) + sqrt2 C^2 + 2`. This is immaterial, since the bound is worse than the
  divisor bound anyway.
* **Prop D.1.** The proof is correct (constant `8 pi^2/6 = 13.16 <= 14`). An independent recount at
  `X = 10^6` (`r3_density_out.txt`) reproduces the author's numbers exactly:
  `E2 = 28594, E3 = 25` for `C = 1`, and `57046, 1061, 6` for `C = 2`.
* **Prop D.2.** Correct, for `C <= sqrt 2`. The summary drops this hypothesis. For general `C` the
  same proof gives `M 2^(-M/ceil(C/sqrt2) - 1)`.
* **§2.2 "Open sharpening".** The claims that `b = O(1)` handles pairs (equal-weight margin `b - 8`)
  and characters with `||c||_1 >= 2.4M`, and the identity `G(H_a) = b + 4M - 12`, re-derive
  correctly. The Weil sketch is labelled as unclaimed.

## 6. Novelty

* `grep` over `research/docs`, `round4/*.md` and `docs/claude_review` found no completion theorem,
  no simultaneous realisation with all monomial gaps, and no fake circle.
* Closest prior material, which the write-up does not cite:
  * `affine_cube_walsh_rigidity.md` §8: "abstract equal-angle and arbitrary-real-angle models";
    these have no gap inequalities;
  * `strict_obtuse_prime_box_countermodels.md`: genuine-prime countermodels to the *pair/norm*
    package at `W = O(M log M)`; no angles, no characters.
* Theorems A.1 and A.3 strictly strengthen both of these and are new relative to the notes.
  Lemma A.7(ii) is walsh-extraction Lemma F.1 / item 392. C.1 is the classical circumradius formula
  with a gcd count. D.1 and the near-axis case of C.2 are elementary, as the author says.

## 7. Corrected statement

1. **Theorem A.1** holds as stated. The constant `1.0367` should be `1.0368`. The "same
   information" claim is restricted to the regime `prod_j (1 + 2/(sqrt p_j - 1)) < 1.75`.
2. **Theorem A.3 and A.3'** hold as stated, with `P ~ 2000 M^5 log M` (any `P >= c M^5 log M` with
   `c` large enough for Lemma A.6) in place of `400 M^5 log M`. The conclusions are unchanged.
3. **Corollary A.4** holds for the classes (I)-(IV) as defined.
   * (IV) may be read as "off a Haar-null subset of the full torus `(R/2 pi Z)^r`" (§3.1).
   * (III) covers residue strengthenings modulo prime powers `<= 2M` (all pair-type pigeonhole
     information, hence every input of `two_thirds`). With the margin available it can be enlarged
     to per-character strengthenings `<= Q_M M^(2||c||_1)`.
   * It does **not** cover arbitrary residue assignments at the larger moduli where pigeonhole on
     products forces multi-point collisions.
4. **Directions (a) and (b)** reduce, in the large-prime regime, to simultaneous rational-character
   admissibility. Arguments using only (I)-(IV) cannot prove
   `M = o((log R/loglog R)^(1/2))`. Whether they can prove `M = o(log R/loglog R)` is **open**: it
   is decided by (P1)/(P2). The claim that these directions cannot give a growth improvement is not
   established.
5. Remark C.3 should be restricted to `C < 2 sqrt 2`, and B.1's summary likewise. The general bound
   in C.2 should read `2 C^(3/2) R^(1/4) + sqrt2 C^2 + 2`. D.2 needs `C <= sqrt 2`, and B.2's slack
   is for balanced profiles.

## 8. Files (`round4/referee_outside_checks/`)

| script | checks | type |
|---|---|---|
| `r1_lemmaA7.py` | Lemma A.7 (i)-(iii) on an independent random-flip build with actual primes in `[P,2P]` and balanced windows; exhaustive small, structured, random and hill-climbed characters | exact integers; weighted slacks via floats, re-verified in 50-digit decimals when small |
| `r2_constants.py` | Lemma A.6 prime requirement versus an exact sieve (`M = 4, 8`) and PNT; the corrected `P`; the step-3 union bound; `Q_M <= e^(2.08M)` (exact, `M <= 2000`); the `arcsin` constant; the per-character `M^(2n)` extension | exact / numeric |
| `r3_density.py` | independent recount of the Prop D.1 data at `X = 10^6` | exact counts |
| `r4_misc.py` | 10-point cluster: arc constant, diagonal, equal chords and `z_a z_d = z_b z_c` (exact); C.1 on all 120 triples; B.2 identity for `M <= 60` | exact |
| `r5_aligned_flats.py` | every aligned Walsh flat on the Sylvester 8/16/32 fakes (`b = 6M`, 5 random flip assignments each): flatness and Lemma A.7(iii), equal-weight form | exact |
