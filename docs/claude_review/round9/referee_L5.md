# Referee report on round9/L5.md (twin channel lemma L5)

Referee role: adversarial. Repository is read-only; my scripts are in
`round9/referee_L5_checks/` and were written independently (the author's code was not imported).

## 0. Verdict (summary)

**The write-up survives as a partial reduction, with six major corrections and no fatal error.** It
claims no growth and proves none, and L5 stays open. `M <= (2/3 + eps) log R/loglog R` is unchanged.

**What I re-derived and confirmed by independent exact computation:**
* the formulation (Def. 1.2 with `log q` weights) and its validity for near-extremal actual twin-profile
  clusters with `O(log M)` margins;
* the identities of §2 (translate form, product structure, class-mate transfer);
* Lemma 3.1 as a per-modulus statement;
* Theorem 3.2 and Cor. 3.3 (switching among ideal candidates kills pair demand on `Q_large`);
* the volume inequality (Thm. 4.1) and its Paley corollary;
* Theorem 5.1.

**What does not survive as written:**
1. The claim that the unweighted wording "is false for every system because every monomial collides at
   `q = 3, 5`" is wrong. Two moduli are absorbed by `O(||w||_1 + 1)`. The real reason for the weights is
   different: genuine labels can collide at `~ log p/(2 loglog p)` small moduli (exact examples), so the
   unweighted wording is not known to be valid for actual circles.
2. Hypothesis IFL is misdescribed. Its part (a) is **not** the hypothetical object of L4. Genuine labels
   of any primes `<= M^(O(1))` satisfy it, by the write-up's own Remark 1.3. The localisation is
   therefore conditional only on the Hall condition (b), which held at all 194 genuine moduli I tested
   (evidence).
3. "No per-modulus or counting argument proves L5" is not what Lemma 3.1 proves. The lemma lives on
   `n_q > M + 3`, i.e. essentially `Q_large`, and says nothing on `Q_med`, where L5 is localised. On
   `Q_med` the statement is a heuristic.
4. Prop. 5.3(2)-(3) is false as written. It makes `kappa'` injective instead of `kappa` (0/1
   counterexample), and it uses the top level `n_Q` instead of level 1 (`q = 3` counterexample). The
   corrected level-1 statement holds, with a greedy proof.
5. The Hadamard-core examples of §5 are not profiles (constant column). The valid ones have rank `M - 1`
   and are outside §5. `d_M = M/2` is a general theorem about a non-profile.
6. The evidence "generic interfaces are dense" is misread. Exact or constructive integral sections of
   random 0/1 profiles have row norm 2-3 at `M = 8` (proven optimal), `<= 11` at `M = 16`, `r = 32`, and
   `<= 6` at `M = 16`, `r = 64`. The square minors gave medians of 407 and 532. At these sizes the twin localisation transfers to generic
   profiles.

Minor points: L5+ is *weaker* than L5, not stronger; Cor. 3.6 needs a hypothesis on sparse elements of
`L`; Lemma 3.4 does not give simultaneous Hall as stated (a sharper count does); and seven more (§2).

**Scope.** The validity labels in `L5.md` §6 are correct, except as amended by Majors 2-6 and Minors 1
and 3. Nothing in the write-up, nor in this report, bounds actual circles beyond `two_thirds`.

## 1. What was checked and how

* Every numbered statement of `L5.md` (§1 Remarks 1.1 and 1-4, Props. 2.1-2.3, Cor. 2.4, Lemma 2.5,
  Lemma 3.1, Hypothesis IFL, Thm. 3.2, Cor. 3.3, Lemma 3.4, Remark 3.5, Cor. 3.6, Thm. 4.1, Cor. 4.2,
  §4.2, Remark 4.3, Thm. 5.1, Props. 5.2-5.3, the evidence of §5 and the table of §6) was re-derived
  by hand. The derivations are summarised in §3.
* The inputs were checked against the sources: `round7/README.md` §3 (wording of L5),
  `round7/p3Dlower.md` Thm. 2.1, (2.2), (2.3), §5.2 (M1') and open points 2 and 4, `sharp.md`
  Def. 1.1 and §4 (Paley twin profile), `round8/stability.md` Cor. 3.6 and Prop. 4.1, and
  `round9/L4.md` §1 (parity-respecting systems).
* Eight independent scripts are in `referee_L5_checks/`, with outputs in `*_out.txt`. They use exact
  integer arithmetic except where floating point is stated. None imports the author's code. The author's
  `check_nontwin.py` was read only to identify which matrix it calls a "Hadamard core".

Labels used below: **proved** (re-derived with all quantifiers), **exact** (finite exact computation),
**evidence** (finite computation that no proof uses), **model**, **heuristic**.

## 2. Issues, ranked

No issue is fatal. The write-up claims no growth and proves none. Its main structural results survive:
the identities of §2, Theorem 3.2 with Cor. 3.3, the volume inequality, and Theorem 5.1. Six claims are
false, unproved or misdescribed as written (Major 1-6). Ten are minor.

### Major

**Major 1. The "q = 3, 5" argument against the unweighted wording is wrong.** (§1 Remark 1; summary
item 1; "the two corrections".)
* *The claim.* The README wording "avoid at all but `O(||w||_1 + 1)` moduli" is "false for every
  system", because every monomial collides at `q = 3` and `q = 5`, and because of pigeonhole at the
  `q` with `m_q/4 < |V|`.
* *Why it fails.* `q = 3, 5` are two moduli. An exception count `O(||w||_1 + 1)` with constant `>= 2`
  absorbs them (`ref_unweighted.py` (a) confirms the trivial fact that every norm-one ratio is a unit
  there).
* The pigeonhole gives, at each small `q`, *some* difference of two members of `V` that collides. That
  difference depends on `q`. It is moreover a free monomial or a `c' = 2e_y` monomial, not a member of
  `V`. Nothing forces a *fixed* monomial to collide at many moduli. No proof is given that the unweighted
  reading fails for every system, and I see none.
* *What is true, and is the real reason the `log q` weights are needed* (exact, `ref_unweighted.py` (b)).
  Genuine labels can collide at unboundedly many **small** moduli. If `pi = a + b i` is a Gaussian
  prime with `5*7*...*q_k | b`, then `pi = conj(pi) mod q` for each of these `q`. The single-label monomial
  then collides with unit 1 at `k` moduli. Examples:
  * `k = 12` with `log p = 70.6`;
  * in general `k ~ log p/(2 loglog p)`, which is unbounded for `p <= M^(h_0)`.

  The weighted demand stays below `(1/2) log p`, as the height bound requires. Placed on a flip column,
  this is the twin monomial `(0, eps e_y)`. So the unweighted reading is **not known to be valid for
  actual circles**, while the weighted (T) is.
* *Correction.* Replace Remark 1 accordingly. Definition 1.2 (weighted) is the right formulation and is
  unaffected.

**Major 2. Hypothesis IFL is misdescribed. Its part (a) is not hypothetical, and only the Hall
condition (b) is open.** (§3.2 Status; §6 table; summary item 4.)
* *The claim.* The write-up says IFL(a) "is exactly the hypothetical object of L4 (an ideal,
  uncertified system)".
* *Why it fails.* Take genuine labels `ell_j(q) = log(pi_j/conj pi_j)` for **any** distinct split primes
  `p_j <= M^(h_0)`, avoiding the moduli. They are parity-respecting: the forced Legendre parity was
  verified on all `(j, q)` with `q <= 1500` in `ref_hall.py` B. By the height bound the write-up itself
  uses in Remark 1.3, `delta_s(v) <= w(v)/2 + (log 2)/2`. So `f` satisfies (F), and `lambda_y = eps_y
  ell_(f(y))` satisfies (T), for all `w` and every `k_max`, with `beta_0 = O(h_0)`.
* This holds on every twin profile in the regime `p_j <= M^(O(1))`. By `stability.md` Cor. 3.6 that is
  the only regime growth needs. Coherent designs with multiplicatively independent `zeta_j` of height
  `O(log M)` (the designs of §4.1) satisfy IFL(a) as well.
* *Consequence.* Theorem 3.2 and Cor. 3.3 hold **conditional only on IFL(b)** (Hall) for genuine
  candidate sets. This makes the result stronger than stated. It also means §6 misstates what is open.
* *Evidence that (b) holds for genuine candidates* (`ref_hall.py` B). On a random twin profile with
  `M = 12`, `r = 36`, genuine primes in `[13, 432]`, `t = 5`, `c = 2` and `S_y = {0, e_(j_1), ..., e_(j_4)}` with random free coordinates (`c_0 = 1`), the
  Hall condition holds at all 194 moduli `q <= 1500` with `n_q >= 2M`. This is evidence only. It concerns
  one profile, and (b) for genuine labels is a statement about coincidences of dense Gaussian integers
  mod `q`.

**Major 3. The barrier "no per-modulus or counting argument proves L5" is overstated.** (§3.1 after
Lemma 3.1; summary item 3.)
* Lemma 3.1 is correct, but it applies only at moduli with `n_q > M - 1 + 4|V|`. Since `|V| >= 1`,
  these are moduli with `n_q > M + 3`. That is essentially `Q_large`, where Cor. 3.3 already says (under
  IFL) that nothing happens.
* On `Q_med`, where the write-up itself localises L5 and where pigeonhole forces pair collisions, Lemma
  3.1 says nothing.
* The "single modulus" half is trivial at **every** modulus: one modulus adds at most
  `A_q log q <= log X = O(A log M)` to any demand, which is within `beta log M` once `beta > A + o(1)`.
* "Counting per modulus cannot prove L5" on `Q_med` is the heuristic "counting forces only averages"
  (`sharp.md` §6.7). It must be labelled heuristic, not attributed to Lemma 3.1.

**Major 4. Proposition 5.3(2)-(3) is false as written: wrong object, and wrong level.**
* *(i) `kappa'` instead of `kappa`.* The pair data are `kappa = A nut + 2 kappa' (mod 2n)`. The recipe
  "`kappa' = d_M u` with `u` injective mod `n/gcd(d_M, n)`" makes `kappa'` injective, not `kappa`.
  Exact counterexample (`ref_nontwin_01.py`):
  * the 0/1 profile `A = [[1,1,0,1,0,0,0],[1,1,1,1,1,1,1],[0,0,0,1,1,0,0],[1,1,0,0,1,0,0]]`, with
    invariant factors `(1,1,1,2)` and every column nonconstant;
  * `n = 8`, `nut = (1,1,1,0,0,1,0)`;
  * `kappa' = (0,2,4,6)`, which is injective, gives `kappa = (2,8,8,14) mod 16`: rows 2 and 3 collide.

  The author's brute-force "part 1" tests only that `H` contains an injective `kappa'`, so it does not
  test the claim.
* *(ii) Top level instead of level 1.* The proposition uses `n = n_Q` (top level). A unit-1 collision
  at `q` is already counted at level 1, i.e. mod `m_q`. Take `q = 3`, `X = 81`, `n_Q = 54`, `M = 5`, and a
  0/1 profile with invariant factors `(1,1,1,1,2)`. The hypothesis `n_Q/gcd(d_M, n_Q) = 27 >= M` holds.
  Yet every `kappa` collides mod `m_3 = 4` (pigeonhole: 2000/2000 random draws, `ref_nontwin.py` (2)).
* *Corrected statement (proved).* If `n_q/gcd(d_M, n_q) >= M`, a parity-compatible `kappa` with no unit-1
  pair collision at `q` (any level) exists. *Proof.*
  1. At level 1, fix a parity class.
  2. Write `(A nut)_x = pi + 2 e_x`.
  3. Choose `kappa'_x = d_M u_x` greedily so that `e_x + d_M u_x` avoids the at most `M - 1` earlier values
     in the coset `e_x + d_M Z/n_q`, which has `n_q/gcd(d_M, n_q)` elements.
  4. Lift to the top level, using that `A(Z/n_Q)^r -> A(Z/n_q)^r` is onto.

  Brute force confirms it on 8172 instances `(A, n)` (`ref_nontwin.py` (1)). The summary's sentence
  (with `n_q`) is the corrected one; the write-up's item 3 (with `n_Q`) is not.

**Major 5. The Hadamard-core examples of §5 are not profiles, and the answer to open point 4 covers
rank-`M` profiles only.**
* *The examples are not profiles.* The matrices `(J + H)/2` with normalised `H` used in §5 (the
  author's `check_nontwin.py` part 4) have a constant all-ones column. That violates `min_x a_xj = 0`
  (`sharp.md` Def. 1.1(a)).
* *The valid Hadamard-core profiles fall outside §5.* These are the nonconstant columns, with any number
  of copies and no flips. They have **rank `M - 1`**, so `Delta_M = 0`, and neither Theorem 5.1 nor
  Prop. 5.3 applies (exact SNF, `ref_nontwin.py` (3)). Flips were introduced in `sharp.md`/`lower.md`
  precisely to reach rank `M`.
* *`d_M = M/2` is a theorem, not five data points.* It holds for every normalised Hadamard matrix of
  order `M >= 4`. By Sherman-Morrison, `((J + H)/2)^(-1) = (2/M) H^T - e_1 e_1^T` (verified exactly on
  all five matrices), so `d_M = M/2`. But it concerns a non-profile.
* *Rank-deficient profiles are not treated.* Rows with a relation `sum c_x a_x = 0` are possible for
  actual clusters, and §5 says nothing about such profiles.
* *Correction.* Restrict the answer to open point 4 explicitly to rank-`M` profiles. Remove "Hadamard
  cores have `d_M = M/2`" as an example of a profile without a unimodular minor.

**Major 6. The evidence on interfaces is misread.** (Prop. 5.2 evidence; summary item 6.) The
write-up reports the row norms of `A_S^(-1)` for unimodular square minors `S`, with medians up to 532 at
`M = 16`, and concludes "generic profiles therefore seem to have dense interfaces".
* *Exact MILP over all integral sections* (`ref_sections.py`, scipy/HiGHS):
  * at `M = 8`, `r = 16`, `sigma(A) = 2` or `3` on all six random profiles, each proven optimal, against
    5-8 from minors;
  * at `M = 12`, `r = 24`, time-limited runs found sections of norm 6 and 8 (median from minors: 25);
  * a cheaper construction assembles column by column `l_1`-minimal integral solutions of `A s_x = e_x`
    (`ref_sections2.py`, every section checked exactly). It gives `sigma(A) <= 5`-`7` on 8 profiles at
    `(M, r) = (12, 24)`, and `sigma(A) <= 7`-`11` on 8 profiles at `(16, 32)`, where the minors gave a
    median of 407. At `(16, 64)` it gives `sigma(A) <= 3, 3, 4, 6` on 4 profiles, where the minors gave a median of 532. At `(24, 48)` the solver found nothing within its time limit, so that case is inconclusive.
* So the reported numbers measure a poor choice of section, not `sigma(A)`. Up to `M = 16` random
  profiles have sections of small norm (2-11). Whether `sigma(A)` stays bounded as `M -> infinity` is
  open (evidence only). That bounded case is exactly the one in which Prop. 5.2(1) transfers the twin
  localisation, so the qualitative answer to open point 4(ii) is reversed at these sizes.
* *Correction.* The gloss must be withdrawn. What is known is an upper bound on `sigma(A)` from minors.
  It has no bearing on whether generic interfaces are dense, and at small `M` the true `sigma(A)` is
  small.

### Minor

**Minor 1. L5+ is a weaker lemma than L5, not a "strengthened" or "stronger" one.** (Cor. 3.6;
summary.)
* (T+) implies (T): take `c' = eps e_y`. So "no data satisfy (F), (T), (P)" implies "no data satisfy
  (F), (T+), (P)", i.e. **L5 implies L5+**.
* The recommendation is still right, for a different reason. L5+ suffices for the route, since actual
  clusters satisfy (T+) by the height bound. It is easier to prove, and its content cannot be localised
  below `q ~ M^(2+eps)` by the switching device.

**Minor 2. Cor. 3.6 has a gap: the switched monomial may lie in `L`.**
* For `c'` with `sum c' = 0`, the shifted monomial `(w - 2 sum_y c'_y w_i^(y), c')` can be a character
  (twin part `c'`, free part `Phi^T c'`) even though `(w, c')` is not.
* The ideal bound of the `lambda`-system is then not available unless "ideal for all sparse monomials"
  includes the sparse elements of `L`. Otherwise one must assume that `L` has no nonzero element of norm
  `<= (1 + 2c_0) k_max`.
* Theorem 3.2 is unaffected, because there `c' = +-e_y`, which is never in `L`.

**Minor 3. Lemma 3.4 is correct but too weak for its use.**
* The algebra and the three examples were re-derived (`2^4 = 16`, `1.5^7 = 17.1`, `1.1^29 = 15.9`, all
  `>= 2e^2 = 14.78`).
* Theorem 3.2 needs Hall at **every** `q in Q_large(c)` simultaneously. The per-modulus bound `1.5/n_q`
  sums to 1.52, 2.61, 1.61, 2.72 over `Q_large` for `(M, A) = (50, 2), (50, 3), (100, 2), (100, 3)`, which
  is more than 1 (`ref_hall.py` A).
* The exact union bound is dominated by `s = 2`, i.e. `O_t(M^(4-2t))/n`. Summed over all moduli it gives
  `< 10^-14` at `t = 5`, and `3 * 10^-7` at `(M, t) = (100, 3)`. Simultaneous Hall therefore holds with
  probability `1 - o(1)` for `t >= 3`, in the model.
* Monte Carlo (4000 trials) agrees: 0 failures at `(M, c, t) = (20, 2, 5)`. Where the lemma's hypothesis
  fails, e.g. `(20, 2, 2)`, there were 9% failures.
* Lemma 3.4 is a model statement and should keep that label.

**Minor 4. "`Q_med` is exactly where pigeonhole forces pair collisions" is inaccurate.**
* Collisions are forced at `q` only when some `sigma_q`-parity class has more than `n_q` rows. For
  balanced classes that means `n_q < M/2`.
* For `c > 1` the range `M/2 <= n_q < cM` lies in `Q_med`, but nothing is forced there.

**Minor 5. Cor. 2.4 has an unnecessary hypothesis.** "With column weights within a bounded factor of
each other" is not needed in the near-extremal regime.
* `w(a_x - a_y) = d_xy >= W/2 >= (3/2) M log M - O(M)`, and `log p_max <= (4/eps - 3 + o(1)) log M`
  (`stability.md` Cor. 3.6). Together they give `||a_x - a_y||_1 >= (3/(2(4/eps - 3)) - o(1)) M`, which is
  dense since `r = Theta(M)`.
* The inequality `||phi_x - phi_y||_1 >= (||a_x - a_y||_1 - 4)/2` was re-derived and checked on all 30
  ordered pairs of a random twin profile (`ref_identities.py`).

**Minor 6. A claim in §4.2 is unsupported but moot.** The sentence "if the free labels are additive
too, the relation `<w, mu> = -eps lambda_y` with sparse `w` violates (T) itself" asserts that such a
sparse `w` exists, without proof. Pigeonhole yields relations among the free labels, which violate
(F), not a relation hitting a prescribed `lambda_y`. The point is moot, because additive free labels
already fail (F) (`p3Dlower.md` Prop. 3.1).

**Minor 7. "Only unit 1 matters for pairs" is a clarification, not a correction.** The README's
"maximal pair agreement" already means `kappa_x = kappa_y`, i.e. unit 1. Remark 1.1 itself is correct:
the budget for units `eta != 1` is `d_xy/2 - O(Delta) >= W/4 - O(Delta)`.

**Minor 8. The threshold in Thm. 4.1 / Cor. 4.2 is astronomical.**
* The geometric-mean bound `M^(1/2)/(C log M)^(b-1)` exceeds `K log M` only once
  `(1/2) log M > (b - 1) log(C log M) + log(K log M)`. For `b = 33` and `C ~ 2.5` that is about
  `M > e^450`. This is harmless for a growth statement but should be said.
* Cor. 4.2(1) is very lossy. Exact computation gives `lambda_min(Gram_A) = 362.97, 395.97, 494.97`
  against the bound `99, 108, 135`, and `max |E_xy| = 2` against 6 (`ref_volume.py`).

**Minor 9. Quantifiers in Cor. 3.3(2).** "On such profiles L5 holds iff `L5_med(c)` holds" should
read as follows. *For every sequence of twin profiles (`M -> infinity`) carrying IFL with constants
independent of `M`, good data exist for all moduli along the sequence iff good data exist on `Q_med(c)`
along it.* L5 is a statement about all profiles, and the equivalence holds only along profiles that
carry IFL.

**Minor 10. Interpretive phrases are presented as conclusions.**
* "Profiles without a unimodular minor ... cannot help a construction" (Prop. 5.3 discussion).
* "An obstruction must therefore control dense combinations ..." (§4.2).

  Neither is proved. The constructor *may* choose twin profiles, but that does not show other profiles
  are useless for the monomial height conditions. Prop. 2.3 is an identity, so the same event is equally
  a statement about the sparse twin monomial of `y`. Both phrases should be labelled heuristic.

## 3. Section-by-section re-derivation

**§1.**
* **Remark 1.1 (correct).** For `eta != 1`, `|(arg eta)/2 - (delta_x - delta_y)/2|` is within `Delta/2`
  of `pi/4`, `pi/2` or `3pi/4`. So `|sin| >= sin(pi/4) - Delta/2`, and (1.1) of `sharp.md` holds once
  `log m_(c,eta) <= d_xy/2 + log(1 - Delta/sqrt 2)`. The budget is `>= W/4 - O(Delta)`.
* **The (P) bound for actual clusters (correct).** `m_xy | beta_xy` gives `m_xy^2 | |beta_xy|^2`. With
  `mu_xy = log(|beta_xy|^2/t_xy^2)` and `t_xy <= C` this gives `2 delta(x, y) <= mu_xy + 2 log C`.
* **Remark 1.3 (validity for actual circles, correct).**
  * Under `H(eps)`, `P_+(N) <= M^(4/eps - 3 + o(1))` (`stability.md` Cor. 3.6(1)).
  * For genuine labels, `m_(v,s) | D_s(Z_v)` with `|D_s| <= sqrt 2 |Z_v|`. Here `D_0 = Im`,
    `D_2 = Re`, `D_1 = x - y`, `D_3 = x + y`, up to units and the factor 2 invisible at odd `q`.
  * Hence `delta_s(v) <= w(v)/2 + (log 2)/2`. With `||v||_1 <= 2||vt||_1` this gives (F) and (T) with
    `beta = 4/eps + o(1)`, for every `k_max`.
  * (P) follows from Remark 1.1 when every margin is `<= K' log M`.
  * So L5 does imply growth for near-extremal actual clusters with a twin profile (`e = 1`) and
    `K'`-bounded margins, and `r <= rho M` holds by Cor. 3.6(3).
* **Remark 1.2 (correct).** The twin part of `sum c_x a_x` is `c` (the twin coordinates of `a_x` are
  `(phi_x, e_x)`, checked in `ref_identities.py`). A zero-sum `c` cannot be `+-e_y`.
* **Remark 1.1 of the list (wrong):** Major 1.

**§2 (all correct; class: proved, valid for actual clusters).**
* *Prop. 2.1.* `<w, f> + eps Lambda_y = s iota` iff `kappa_y = <phi_y - eps w, f> + eps s iota`, since
  `eps^2 = 1`.
* *Prop. 2.2.* This is `p3Dlower.md` Thm. 2.1(1),(4). The inverse map is `ell_(f(x)) = eps_x Lambda_x`,
  `ell_(s(x)) = U_x - eps_x Lambda_x`. (F) depends on `f` only, (P) on `kappa` only, and (T) on both.
* *Prop. 2.3.* If `kappa_x = kappa_y`, then `Lambda_y = kappa_x - F_y`, so
  `<w, f> + eps Lambda_y = <w + eps(phi_x - phi_y), f> + eps Lambda_x`.
* *Cor. 2.4.* For `z != x, y` the twin relation gives `v_(f(z)) = v_(s(z))`, hence
  `||v||_1 <= 2||(v_J, v_s)||_1 + 4`. At most `n_q` rows per parity class are alone.
* *Lemma 2.5.* Demands are sums of nonnegative terms over `q`.
* *Exact.* `ref_identities.py` checks Prop. 2.1 (17,600 instances), Prop. 2.3 (8,440), Prop. 2.2 round
  trips (440), (2.2) and `kappa = F + Lambda`. It uses all odd `q < 60` at levels 1-2, genuine Legendre
  parities, and 0 failures.

**§3.**
* **Lemma 3.1 (correct as stated).**
  * A collision at level `a` implies one at level 1, with some unit. The reduction maps `T` onto `T`, though it may
    send `i` to `-i` depending on generators, which is harmless here.
  * At most `M - 1 + 4|V_y|` values are excluded in a parity class of size `n_q`.
  * Its reading: Major 3.
* **Theorem 3.2 (correct).**
  * Hall per parity class gives distinct level-1 values, hence no unit-1 collision at any level.
  * `Lambda_y - lambda_y = -2<w_i, f>` is even, so the parities are preserved.
  * At each `q` the twin monomial `(w, eps e_y)` of the switched system has exactly the collision pattern
    of `(w - 2 eps w_(i(y,q)), eps e_y)` in `(f, lambda)`. Summing over `q`, with the sum over `i < t`
    dominating the term at `i(y,q)`, gives `delta_s <= t beta_0 (||w||_1 + 2c_0 + 1) log M`.
  * This needs (T) for `lambda` up to `k_max + 2c_0`, which IFL(a) provides.
* **Cor. 3.3 (correct)** up to the quantifier remark (Minor 9).
* **Status of IFL:** Major 2.
* **Lemma 3.4:** Minor 3.
* **Remark 3.5.**
  * `Q_large(c) ⊂ (2cM - 1, 3M]` for `A = 1`; its reciprocal mass is `O(1/log M)`.
  * The reciprocal mass of `(2cM, 3M^A]` is `log A + o(1)`.
  * The union bound `M^2/(2 n_q)` for injectivity is correct.
  * All are checked.
* **Cor. 3.6:** Minors 1-2.

**§4.**
* **Theorem 4.1 (correct).** Every step re-derived:
  * `prod s_i = det(G^T G)^(1/2) >= 1` for an injective integer `G`;
  * `prod_(i <= k) s_i = ||Lambda^k G||_op <= ||Lambda^k G||_F`, and Cauchy-Binet with Hadamard gives
    `||Lambda^k G||_F^2 <= binom(r, k) B^(2k)`;
  * Cauchy interlacing for the compression of `G^T G` to `V` gives `t_i >= s_(i+k)`;
  * Hadamard's inequality for `vol(G b_x)`.

  Tested on 542 injective instances, each with a 60-step search minimising the slack. The minimum slack
  was 0.347 (`ref_volume.py` (2); floating point, evidence).
* **Cor. 4.2 (correct).**
  * `Gram_A = (bM/4)(I + J) + E`, with each flip perturbing `4 Gram_A` by at most 24.
  * Weyl's inequality, `binom(b(M-1), M-1) <= (eb)^(M-1)` and `r - M + 1 = (b-1)(M-1)` complete the
    proof.
  * The twin coordinates of `e_(f(y))` and `e_(s(y))` are as stated.
  * *Exact* (`ref_volume.py` (1), independent Paley-I construction):
    `(1/2) log det Gram_A = 128.6203, 142.4975, 185.0806` for `q_0 = 43, 47, 59` (`b = 33`), against
    `((M-1)/2) log(bM/4) + (1/2) log M = 128.6218, 142.4988, 185.0816`.
  * Scope as stated by the author (construction method only); see Minor 8.
* **§4.2 and Remark 4.3** are correctly labelled heuristic, except for the phrases in Minors 6 and 10.

**§5.**
* **Theorem 5.1 (correct, rank-`M` profiles).**
  * Onto mod `p` lifts to onto mod `p^e` (Nakayama), then use CRT.
  * A surjection onto a free module splits.
  * `<v, ell> = <v, nut> + 2<s^T v, kappa'> + 2<v, k>`.
  * `s(e_x) = eps_x(e_(f(x)) - e_(s(x)))` is a section for twin profiles.
* **Prop. 5.2 (correct, definitional).** The evidence gloss is wrong: Major 6.
* **Prop. 5.3.**
  * (1) is correct (Smith form). The non-split example `A = (2)`, `n = 4` is correct.
  * (2)-(3): Major 4.
  * The Hadamard examples: Major 5.

**§6 table.** Its rows must be amended by Majors 2-6 and Minors 1, 3. Every other row is accurate.

## 4. Independent computations

All scripts are in `referee_L5_checks/` and run with `python3 <script>` from that folder. Outputs are in
`<script>_out.txt`.

| script | what | result |
| --- | --- | --- |
| `ref_identities.py` | random twin profile (`M = 6`, `r = 20`), genuine Legendre parities, all odd `q < 60`, levels 1-2: `kappa = F + Lambda`, (2.2), Prop. 2.1, Prop. 2.2 inverse map, Prop. 2.3 (with `Lambda_y` re-chosen to force `kappa_x = kappa_y`), Cor. 2.4 norm inequality, twin coordinates of `a_x` | 0 failures (exact) |
| `ref_unweighted.py` | (a) at `q = 3, 5` every norm-one ratio is a unit; (b) Gaussian primes `a + b i` with `5*7*...*q_k \| b` | (a) true; (b) `k` unit-1 collisions at small moduli for `k = 1..12` (`log p` up to 70.6), weighted demand `<= (1/2) log p` (exact) |
| `ref_hall.py` A | Hall failure for `t` uniform candidates; exact union bound; sums over `Q_large` | Monte Carlo agrees with the lemma where it applies; `sum 1.5/n_q` = 1.5-2.7 over `Q_large` vs exact union bound `< 10^-14` (`t = 5`) (model) |
| `ref_hall.py` B | genuine labels on a random twin profile (`M = 12`, `r = 36`, primes `<= 432`), Hall at every `q <= 1500` with `n_q >= 2M`, `t = 5`; forced Legendre parity of genuine labels | Hall holds at all 194 moduli; parity holds at every `(j, q)` (evidence; parity exact) |
| `ref_volume.py` | independent Paley-I construction; exact `det Gram_A` (Bareiss) for `(q_0, b) = (43,33), (47,33), (59,33), (43,25), (47,26), (11,30)`; `lambda_min`, `max\|E\|`; inequality (4.1) under adversarial local search | matches `(bM/4)^(M-1) M` to `< 0.006` in `log`; `lambda_min ~ bM/4`; `max\|E\| = 2`; (4.1) minimum slack 0.347 on 542 instances |
| `ref_nontwin.py` | Prop. 5.3(2) recipe vs corrected claim (8172 instances `(A, n)`, brute force over `H`); level example at `q = 3`; own Smith normal form of `(J+H)/2` for Sylvester 8, 16 and Paley 12, 20, 24; constant columns; rank of the valid core profile; Sherman-Morrison inverse; `Delta_M` of random profiles | recipe fails, corrected claim holds; 2000/2000 level-1 collisions; invariant factors agree with the author; constant column 0 in all five; valid core profiles have rank `M - 1`; inverse identity exact; `Delta_M = 1` in 80/80 random profiles |
| `ref_nontwin_01.py` | 0/1 counterexample to the Prop. 5.3(2) recipe | found: `M = 4`, `r = 7`, `n = 8` |
| `ref_sections.py` | exact MILP for `sigma(A)` (scipy/HiGHS, time-limited) | `M = 8`: `sigma = 2, 3, 2, 3, 2, 3` (optimal); `M = 12`: 6, 8, 50 found in 60 s; `M = 16`: none found in 90 s |
| `ref_sections2.py` | column-wise `l_1`-minimal integral sections (upper bounds on `sigma(A)`) | `M = 12, r = 24`: 5-7; `M = 16, r = 32`: 7-11; `M = 16, r = 64`: 3-6 (every section verified exactly); `M = 24, r = 48`: no section within the 20 s per-column limit (inconclusive) |

## 5. Corrected statement

**L5 (the twin channel) is open.** No construction and no obstruction is proved, and no growth
statement is proved or claimed. `M <= (2/3 + eps) log R/loglog R` is unchanged. What is established:

1. **Formulation (proved; valid for actual circles).**
   * L5 is stated with `log q`-weighted demands, (F), (T), (P) as in Def. 1.2, with (P) for unit-1
     collisions.
   * (F), (T), (P) hold for near-extremal actual clusters with a twin profile (`e = 1`) and margins
     `<= K' log M`, with `beta = 4/eps + o(1)` and `K = K'/2 + O(1)`. So L5 would give growth for that
     subclass.
   * The unweighted wording is **not known to be valid for actual circles**: genuine labels can collide
     at `~ log p/(2 loglog p)` small moduli (exact examples up to 12 moduli). It is **not** shown to be
     false for every system. The `q = 3, 5` argument is void.
2. **Identities (proved; valid for actual clusters).** Props. 2.1-2.3, Cor. 2.4 and Lemma 2.5 as
   written. In Cor. 2.4 the bounded-weight-ratio hypothesis is unnecessary in the near-extremal regime.
3. **Per-modulus freedom (proved).** Lemma 3.1 as written, valid at moduli with
   `n_q > M - 1 + 4 max_y |V_y|`. It says nothing about `Q_med`. That a single modulus cannot prove L5
   is trivial (one modulus adds `<= log X` to any demand). That counting cannot prove L5 on `Q_med` is a
   **heuristic**.
4. **Localisation (proved, conditional on Hall only).**
   * *Setting.* Let a twin profile have primes `p_j <= M^(h_0)`. Then genuine labels (or any
     multiplicatively independent coherent design of height `O(log M)`) satisfy IFL(a) with
     `beta_0 = O(h_0)`.
   * *Hypothesis.* For some row sets `S_y` with `||w_i^(y)||_1 <= c_0`, the Hall condition IFL(b) holds at
     every `q in Q_large(c)`.
   * *Conclusion.* There are data on `Q_large(c)` with zero pair demand, (F) with `beta_0`, and (T) with
     `t beta_0 (1 + 2c_0)`. Along any sequence of such profiles with uniform constants, good data exist
     iff they exist on `Q_med(c)` (L5 iff `L5_med(c)`).
   * *Status of the hypothesis.* Hall for genuine candidates is open. Evidence: it holds at all 194 moduli
     tested. In the random-candidate model, Hall holds simultaneously at all `q in Q_large` with
     probability `1 - O_t(M^(4-2t) log(A+1))`, using the corrected Lemma 3.4 with the `s = 2` term.
5. **L5+ (proved implication, conditional localisation).** L5 implies L5+, and L5+ is what the route
   needs. Under a common-index injectivity hypothesis on `n_q >= M^(2+eps)`, its content lies on
   `n_q < M^(2+eps)`. This also needs the extra hypothesis that the ideal bound covers the sparse elements
   of `L` (or that `L` has none of norm `<= (1 + 2c_0) k_max`).
6. **Certification barrier (proved; construction method only).** Thm. 4.1 for every profile with
   affinely independent rows and multiplicatively independent coherent columns. Cor. 4.2 on the Paley twin
   profile. The threshold beyond which the pair-height bound exceeds `K log M` is astronomical
   (`M >~ e^450` at `b = 33`).
7. **Non-twin profiles (proved for rank `M`).**
   * Thm. 5.1 holds for every rank-`M` profile.
   * Prop. 5.3 holds **with level-1 quantities**. If `n_q/gcd(d_M, n_q) >= M`, a collision-free
     parity-compatible `kappa` exists. The proof is a greedy choice in the cosets `e_x + d_M Z/n_q`, not
     "`kappa' = d_M u` injective". Pair collisions can therefore be forced, by the restriction to `H`, only
     where `n_q < M gcd(d_M, n_q)`.
   * Rank-deficient profiles, which include the Hadamard-core profiles without flips, are **not**
     covered.
   * `d_M((J + H)/2) = M/2` holds for every normalised Hadamard `H`, but `(J + H)/2` is not a profile.
   * The claim that generic interfaces are dense is **withdrawn**. Exact or constructive sections of
     random 0/1 profiles have row norm 2-3 (`M = 8`, proven optimal), `<= 7` (`M = 12`), `<= 11`
     (`M = 16`, `r = 32`) and `<= 6` (`M = 16`, `r = 64`) (evidence). At these sizes Prop. 5.2(1) transfers the twin localisation to
     generic profiles. Whether `sigma(A)` stays bounded as `M -> infinity` is open.
