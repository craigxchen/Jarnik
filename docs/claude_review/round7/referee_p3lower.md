# Referee report on `round7/p3lower.md` (lower bounds from factored residue data, P3-all / P3-gen)

Referee: adversarial, independent. Every proof in §§1-5 of the write-up was re-derived by hand. Every
numerical claim that a conclusion leans on was recomputed with new code. The author's scripts were not
imported or rerun as a substitute. New scripts and their outputs are in
`round7/referee_p3lower_checks/`; each runs from that directory with `python3 <script>`. The repository,
`round5/`, `round6/` and `round8/` were not modified.

Conventions follow the write-up: `W = log N`, `mu_xy = d_xy - W/2`, `D_xy = 2 log m_xy`,
`c_C = 2 log C - log 2`, `L = log M`, `l = log L`.

---------------------------------------------------------------------------------------------

## 0. Verdict

**The result survives. No theorem is wrong, no growth is claimed that is not proved, and the honest
headline ("target not reached") is correct.** Every proved statement checked out under re-derivation:

* Lemmas 1.3 and 1.4;
* Theorem A (LP duality, Perron and hypermetric lower bounds, Goemans-Williamson upper bound);
* Corollaries A.1-A.3, Lemma A.4 and Corollary A.5;
* Propositions B.1-B.3 and B.5;
* Lemma C.1, Theorem C, Corollary C.4, and Lemmas C.2 and C.3;
* Lemma D.0 and Proposition D.1.

Independent exact or LP computation confirms the sandwich, the `D_A` point-total bound, the parity
mechanism on actual Gaussian integers, and the clique-block Gram matrices.

There are **no fatal and no major errors**. There are **eleven minor issues**. Most of them are
overstatements in the prose or tables around correct theorems:

* one stated claim is false (the "PSD = Perron" remark when `C > sqrt 2`);
* one is false as phrased (parity "at most a factor 2 per modulus");
* one numerical slip of a factor 2 (`D_max` of `D_A` in §6);
* the rest are scope and labelling problems.

The single most important caveat for users of the write-up is a matter of scope, not of correctness.

* **Theorem C carries no arithmetic content beyond its hypothesis.** It is Proposition B.1 applied to
  one pair with `>= k ~ 2L/l` window collisions.
* Hypothesis `H` as stated (for *every* `k`) is far stronger than what the proof uses (*one* `k`).
* The "growth" is assumed, not derived. The conditional result is correct, but it does not reduce the
  growth problem to a more accessible statement. It restates the Poisson tail as an assumption.

| # | where | severity | issue (one line) |
|---|---|---|---|
| 1 | §2 after Cor. A.1; §0 item 4; §0 table A.1; mechanism table | minor (false claim, `O(M)` effect) | "the strongest consequence of PSD + `mu >= D - c_C` is exactly the Perron bound" is false when `c_C > 0`: explicit `M = 5` counterexample (Horn matrix), PSD gives `W >= 12/5`, Perron-type sup gives 2, Cor. A.1 gives 0 |
| 2 | §0 table B.3; §2 Cor. A.3 bullet; §7 item 1 | minor (false as phrased) | "parity gains at most a factor 2 per modulus" is false per modulus (ratio unbounded when `M <= m_Q < 2M`; 4 and 6.5 in finite cases). Aggregated, the ratio is 2.15-2.56 at `M = 10^2..10^5`. The proved statement is the additive bound B.3(5) |
| 3 | §6, row "Prop. B.1" | minor (numerical slip) | `D_max <= (10 + o(1)) log M` should be `(20 + o(1)) log M` (Lemma 3.3(1) bounds `log m`, and `D = 2 log m`). The author's own §2 table shows `D_max/log M = 12.4`. The consistency conclusion is unaffected |
| 4 | §4 (Theorem C, `H`), §0 table C, §7 (H-actual) | minor (scope / labelling) | `H` is assumed for all `k`, but one `k` is used. Theorem C is B.1 plus "one pair with `>= k` window collisions", and the growth is in the hypothesis |
| 5 | §4, description of `H` for actual clusters | minor (imprecise) | "Equivalently ... `Im(z_x conj(z_y) u)/|g_xy|^2`-type integers" is not an equivalence. `q | Im(z_x conj z_y)/|g|^2` characterises collisions with unit `+-1`; unit 1 needs in addition `Re(.) = |B|^2 mod q`. Exact check: 10,892 unit-`(-1)`-only instances |
| 6 | Lemma C.3 and its use (§0 table, Remark 1 on Thm C) | minor (scope) | the independent model ignores the forced Legendre parity, so it is an `A_X` model, not a P3 model. Remark 1 ("without uniform margins `H` gives nothing beyond Cor. A.1") needs a demand matrix satisfying `H` with `lambda_max = O(M log M)`. That is evidence only (no proof of `lambda_max` for the model) |
| 7 | §0 item 3; §6 row "sharp.md §6.6" | minor (scope) | "at the real-weight level the obstruction disappears for profiles adapted to the demand pattern": in P3 (and P3-gen) the demand matrix is a function of the profile (`kappa = A ell`, `sigma_q(x) = sum_j a_xj nu_j(q)`), so "adapting the profile to `D`" is not well defined. The supporting `lambda_max(D) ~ 2 M log X` is evidence |
| 8 | §3 after Prop. B.5; mechanism table ("height condition"); §7 item 4 | minor (overstated) | "Pigeonhole consequences of the heights are subsumed by the trivial count of primes" is proved only for the two-term (`v = +-e_j +- e_l`), single-modulus pigeonhole. The "norm `<= k`" sentence is unproved |
| 9 | Remark B.4 | minor (non sequitur) | "Since `alpha >= 2^(-|S|)`" should read "since `alpha <= 1`". "`(2 + o(1)) log M` at best" bounds this CRT argument, not what can be forced |
| 10 | §2 evidence and §6 "measured" values | minor (evidence scope) | `W_cut/(2 lambda) in [1.000, 1.224]` is specific to the author's 48 matrices: independent LP gives 1.443 on a residue matrix (`M = 4`) and 1.500 adversarially. `2 lambda_max(D_A) ~ 6.4 M log M` depends on the colouring: constant parity gives `~ 13 M log M` (still below the proved `72 M log M + 222 M`) |
| 11 | Remark 2 on Thm C; §7 item 2; attribution | minor (update / attribution) | "Theorem L already penalises actual Hadamard-core clusters" needs the round-8 qualification (`round8/lw.md` Thm 5.1: the Paley skew double core is immune to every single LW form). Theorem A(3) is Nesterov's `pi/2` theorem / Goemans-Williamson rounding in LP-dual form, and the cut cone and hypermetric inequalities are Deza-Laurent; these should be cited |

Sections 1-4 below give the derivations, the computations, the class check and the corrected statement.

---------------------------------------------------------------------------------------------

## 1. Statement-by-statement re-derivation

### 1.1 Classes and inputs (§1)

* **Prop. 1.2 (actual clusters lie in P3-gen, hence in P3-all, P3-char and `A_X`).** Correct. It is
  `sharp.md` Prop. 1.2 with `r_j(q) = pi_j mod q^(A_q)`, which is genuine by Definition 1.1. "`F^char ⊂ A_X`"
  is an induced map on residue data (`rho_x(q)` has norm `N`), not a literal inclusion of tuples. That is
  harmless.
* **P3-gen facts.** Checked.
  * The reduced chord `u B - u' conj(B)` is nonzero, because `B` is conjugate-primitive and not a unit.
  * It has modulus `<= 2 e^(d/2)`.
  * The heights are automatic, and the slabs are not.
* **Lemma 1.3.** Re-derived. Take (1.1) of `lower.md` with `c = e_x - e_y` and `|delta_x - delta_y| <= Delta`.
  This gives `m e^(-d/2) <= 2^(-1/2) C e^(-W/4)`, i.e. `mu_xy >= D_xy - c_C`. Correct, including the
  sign of `c_C`.
* **Lemma 1.4.** Correct.
  * Items 1-2: normalisation, and `|a - b| = #{t : exactly one of a, b >= t}`.
  * Item 3: `sum_(y != x) (-1/2) s_x s_y = 1/2 - (1/2) s_x sigma(s)`.
  * Item 4: sum item 3 over `x` and halve.
  * The threshold form `sigma^2/4 = (S_jt - M/2)^2` is right.

### 1.2 Theorem A and its corollaries (§2)

* **A(1), LP duality.** Correct.
  * There are `(2^M - 2)/2 = 2^(M-1) - 1` cuts.
  * The uniform cut measure separates each pair with probability `2^(M-2)/(2^(M-1) - 1) > 1/2`, so the
    primal is feasible. It is bounded below by 0.
  * The dual constraint is `cut_nu(s) - |nu|/2 <= 1`, i.e. `surplus(nu) <= 1`.
  * `surplus(nu) > 0` holds because the average over *nonconstant* cuts exceeds `|nu|/2`.
  * The ratio form follows by scaling whenever the optimum is positive.
* **A(2), Perron.** Correct. For `nu = 4 u u^T` the constraint value is `|u|^2 - (u.s)^2 <= 1` and the
  objective is `2 u^T B u` (zero diagonal). The hypermetric version is also correct: `b.s = sum b mod 2`
  is odd, so `(b.s)^2 >= 1`.
* **A(3), GW.** Correct.
  * `B >= 0` with zero trace and `lambda_max = 0` forces `B = 0`.
  * `2 B_xy = (e_x + e_y)^T B (e_x + e_y) <= 2 Lambda`, so `rho_xy in [-1, 0]`.
  * Sheppard gives `P(s_x != s_y) = arccos(rho)/pi`. Dropping constant cuts keeps
    `mu_xy = pi Lambda (P(s_x != s_y) - (1 - pi_0)/2) >= Lambda arcsin(-rho) >= B_xy`.
  * `W(lambda) <= pi Lambda`.
* **Cor. A.1.** Correct. `u` is the Perron vector of `D`, and `1 <= (sum u)^2 <= M` because `u >= 0`
  and `|u| = 1`. The remark that follows it is **not** correct in general (issue 1, §2.1 below).
* **Cor. A.2.** Correct, and tautological by design. A pair-cut certificate is *defined* to hold on all
  cut-realisable `(W, mu)` with `mu >= D - c_C`, so `Phi(D) <= W_cut(D - c_C(J-I))`.
  * The entrywise bound `(D - c_C(J-I))^+ <= D + c_C^-(J-I)` holds.
  * The Perron root is monotone, and Weyl subadditivity gives the stated cap.
  * **Scope:** `Phi` sees only `D`. Arguments that use the *coupling* between `D` and the profile are
    not covered. Examples are the same-level clause at `p_j <= X` (which pairs share a level is a
    property of the cut decomposition), the chain lemma of `two_thirds`, and in P3 the dependence of
    the point classes on the profile. The list "This includes ..." should say so (issue 7).
* **Cor. A.3.** Correct. It uses the row-sum bound `lambda_max(D) <= max_x D_x`.
* **Lemma A.4.** Re-derived line by line.
  * Unit-1 pair collisions occur only at small `q` and levels `a <= a*(q)` (`sharp.md` Lemma 3.3(1)).
  * This gives `log m_xy <= Lambda'(d) + log d`, and each `d` occurs at most twice in a row.
  * `p'(q) | d < M` needs `q < 8 p'(q) <= 8M`, so the first sum is `<= 8M sum_(q <= 8M) log q/q`.
  * The elementary bound gives `<= 8M(log M + log 32)`.
  * `4(9 M log M + 8M log 32) = 36 M log M + 110.9 M`.
  * Correct for every parity colouring and every random choice. Confirmed numerically (§3.2).
* **Cor. A.5.** Correct. The Theorem 6.4 members have `p_j > X`, so there is no same-level clause.
  `sharp.md` Lemma 6.3(1) says the random fibre shifts all point classes by a common `2t`, so the pair
  collisions are those of `D_A`. Combined with Cor. A.2 and Lemma A.4, every pair-cut certificate gives
  `<= pi(36 M log M + 111 M + (M-1) c_C^-)` on those members, for all large `M`, at `X = 3M^A`
  (`A = 1` included). So no pair-cut certificate proves `W_min^all >= f(M) M log M` with `f -> infinity`.
  This is the one genuinely member-level barrier in the write-up, and it is sound.

### 1.3 Restricted bounds (§3)

* **`E(n, Q)` facts.** Correct: the closed form, `n^2/(2Q) - n/2 <= E <= n^2/(2Q)`, and superadditivity
  (exhaustively rechecked, §3.4).
* **Prop. B.1.** Correct. It needs only the cut identity and the pair inequality.
* **Remark 2 after B.1** (Paley margins). Re-derived. For a pair, `|tau_j| in {0, 2}`, and
  `2B - r = b + (flip terms) in [b - 4, b + 4]`. With `w_j in [tau, tau + 1/r]` this gives
  `mu in [(b-4) tau/2 - 1/2, (b+4) tau/2 + 1/2]`, so `K <= 37/29 + o(1)` for `b = 33`. Correct.
* **Prop. B.2.** Correct. It uses Lemma 1.4(3) and `|s_x sigma| <= |sigma|`.
* **Prop. B.3.** Items 1-5 are correct.
  * Item 1 is `sharp.md` Cor. 2.2(2),(4) and Prop. 2.1(6). Since `ell_j = nu_j mod 2` and
    `iota = lambda_q mod 2`, `kappa_x = sigma_q(x) mod 2`, and `m_(q^a)` is even.
  * Item 3: the lower bound uses `n_0^2 + n_1^2 = M^2(1 + beta^2)/2`.
  * Item 4: collisions are nested in `a`, and the cut identity is applied.
  * Item 5: for `beta = 0` and `M` even the parity count *equals* `E(M, m)`. This was checked:
    `2E(M/2, m/2) = E(M, m)`.
  * What is **wrong is the gloss** "at most a factor 2 per modulus" (issue 2).
  * Verified on actual lattice points (§3.3).
* **Remark B.4.** The argument is right up to a non sequitur (issue 9).
* **Prop. B.5.** Re-derived. Two points need checking.
  * `v_P != 0`, including `j = l`, where `v_P = 2 eps e_j`.
  * The weight sum over unordered pairs of the `2r` signed generators: 4 sign patterns per `{j, l}`
    with `w(v_P) = w_j + w_l`, plus one pair per `j` with `w = 2 w_j`. So
    `sum_P w(v_P)/2 = 2(r-1) sum w + sum w = (2r - 1) sum w`.

  The final bound uses Mertens' theorem and `theta(Y) <= Y log 4`. Correct. The claim drawn from it
  is overstated (issue 8).

### 1.4 Theorem C and the `H` lemmas (§4)

* **Lemma C.1.** Correct. It uses ordered selection and `lambda - i p_max >= lambda - k p_max`. The
  exact window `e_k` satisfy both bounds (§3.5).
* **Theorem C.** Re-derived.
  * *Step 1.* `k <= L` and `p_max < 1/(M^(1-eta) - 1)`.
  * *Step 2.* `log(delta^k E_k) >= 2L - log 4 - k log k + k log(delta lambda_0/2)`, and
    `k log k <= (2 - eps/4)(L + O(L/l))`. Correct.
  * *Step 3.* `F_k >= 1` gives some `c_xy >= k`, and `log q > (1 - eta) L` on the window.
  * *Step 4.* `2(4 - eps/2)(1 - eta) = 8(1 - eta) - eps(1 - eta) > 8(1 - eta) - eps`.

  All correct. `M_0` depends on `(eta, delta, lambda_0, C, eps)`.
* **Theorem C is not vacuous at moderate `M`.** For `eta = delta = eps = 1/2` and the exact window
  primes, `delta^k E_k >= 1` already at `M = 10^4` (`k = 7`, `log(delta^k E_k) >= 0.89`; §3.5).
* **Cor. C.4.** Correct. It uses `M <= log R` (from `two_thirds`) to get `loglog M <= logloglog R`,
  and the case split at `(log R)^theta`.
* **Lemma C.2.** Correct.
  * `q` in the window is a small prime of the design.
  * A collision needs `p'(q) | d`, and `p'(q) > M^(1-eta)/8`.
  * At most `2/(1-eta)` such primes divide `d < M`, with at most 9 preimages each, so
    `c_xy <= floor(18/(1-eta))`.
  * Meanwhile `delta^k E_k -> infinity` at the fixed `k = c_0 + 1`.
  * The random fibre is a common shift.

  Measured `max c_xy` for `eta = 1/2` is 7-8 for `M <= 1024` (§3.2), well inside `c_0 = 36`.
* **Lemma C.3.** Correct.
  * *Pairwise independence.* A shared point gives `P(kappa_x = kappa_y = kappa_z) = 1/m^2`.
    Coordinatewise independence across `q` then makes the vectors `(1[P at q])_q` and `(1[P' at q])_q`
    independent.
  * *Second moment.* `E binom(c,k)^2 <= sum_i e_i e_(k-i)^2` follows by dropping disjointness.
  * *Chebyshev* with `1 - delta^k >= 1 - delta`.
  * *Large `k`.* `binom(k,i)^2 i!/k! = binom(k,i)/(k-i)!`, so the bound is
    `(lambda/lambda')^(2k)(delta(1 + lambda))^k`.
  * *Beyond `k > log M`.* `E_k < 1` there.

  The model's scope is issue 6.

### 1.5 Discreteness (§5)

* **Lemma D.0.** Correct. The sparse-cut inequality `W >= min_r max(f_D(r), r log X)` is valid when
  all `p_j > X` and `r` counts threshold cuts (`sum_j e_j`), not columns. The text says "columns". With
  `e_j > 1` a column carries `e_j` cuts, so "number of cuts" is the right reading. Minor wording.
* **Prop. D.1.** Correct.
  * Item 1 needs `k` with zero diagonal (`[same group]` must mean `x != y` in the same group).
  * Item 2: the Gram entries are `(g-1) - m = -L`, `-1`, and `(g-1) + (t-1)m = r`.
  * Item 3: the eigenvalues are `M - 1 + (L-1)(t-1) = r`, `-1 + (L-1)(t-1)` and `-L`.

  Exactly rechecked in 10 Sylvester cases up to `M = 128` (§3.6).
* **Question Q.** It is a question, so nothing to refute. Note that it uses `k = max(1, ceil(2D/tau))`,
  which drops `c_C` and forces correlation `<= -1` even where `D_xy = 0`. That is a stronger
  requirement than the pair inequality. It should be stated as a modelling choice.

---------------------------------------------------------------------------------------------

## 2. The issues in detail

### 2.1 Issue 1: "PSD-ness carries exactly the Perron bound" is false when `c_C > 0`

**Claim** (after Cor. A.1): "the strongest consequence of '`W I - 2 mu` is PSD and `mu >= D - c_C`' is
exactly the Perron bound, **PSD-ness of the margin Gram matrix carries no information beyond
Corollary A.1**." It is repeated in §0 item 4, in the §0 table row A.1 ("contains ... every PSD
consequence"), and in the mechanism table.

**What is true.**
* The dual of `min{W : W I - 2 mu ⪰ 0, mu_xy >= B_xy}` is `2 max{<Z, B> : Z ⪰ 0, Z >= 0, tr Z = 1}`,
  a maximum over *doubly nonnegative* (DNN) matrices.
* For `B >= 0` this equals `2 lambda_max(B)`, so the claim holds when `c_C <= 0` (`C <= sqrt 2`).
* For `B` with negative entries, DNN is strictly larger than the cone generated by `u u^T` with
  `u >= 0` once `M >= 5`. So PSD can beat `2 sup_(u>=0) u^T B u`.

**Counterexample** (`psd_horn.py`, exact in `Q(sqrt 5)`). Data:
* `M = 5`, `D = 2 x` (adjacency of the 5-cycle), `c_C = 1`;
* then `B = D - (J - I) = I - Horn`;
* `Z = circ(1, 3/5, 0, 0, 3/5)/5`. It is DNN: its eigenvalues are `11/25` and `(7 +- 3 sqrt 5)/50 > 0`.

Every `(W, mu)` with `W I - 2 mu ⪰ 0` and `mu >= B` satisfies `W >= 2 <Z, mu> >= 2 <Z, B> = 12/5`.

| bound | value |
|---|---|
| PSD bound | `>= 12/5` |
| Perron-type bound `2 sup_(u>=0,|u|=1) u^T B u` | `2` (Horn is copositive; attained at `(e_1 + e_2)/sqrt 2`) |
| Corollary A.1, `2 lambda_max(D) - 2(M-1) c_C` | `0` |
| `W_cut(B)` | `10/3` |

**Effect.** Always `|PSD bound - Perron bound| <= 2(M-1)|c_C|`, which is `O(M)`. Nothing about growth
changes.

**Correction.** "For `C <= sqrt 2` PSD-ness gives exactly `2 lambda_max(D + |c_C|(J-I))`. In general it
gives `2 max_(Z DNN) <Z, D - c_C(J-I)>`, which lies within `2(M-1)|c_C|` of `2 lambda_max(D)`. It is
never more than `W_cut`."

### 2.2 Issue 2: parity is not "at most a factor 2 per modulus"

The proved statement B.3(3)-(5) is right:
* `E(M, m) <= E(n_0, m/2) + E(n_1, m/2) <= E(M, m/2)`;
* the additive gain per modulus is `<= M^2/(2m) + M/2`, and it vanishes for `m >= 2M`;
* the total gain is `O(M log M)`.

The phrases in the §0 table ("This gains at most a factor 2 per modulus"), in the Cor. A.3 bullet and in
§7 item 1 ("Parity raises averages by at most a factor 2") are false as stated.

* **Per modulus** (`e_counts.py`):

  | `(M, m)` | parity-free `E(M, m)` | worst parity count | ratio |
  |---|---|---|---|
  | `(10, 12)` | 0 | 4 | unbounded |
  | `(10, 8)` | 2 | 8 | 4 |
  | `(40, 36)` | 4 | 26 | 6.5 |

  The window `M <= m_Q < 2M` always gives an infinite ratio under full bias.
* **Aggregated B.3(4) main terms** at `X = 3M` (full bias / parity-free):

  | `M` | `10^2` | `10^3` | `10^4` | `10^5` |
  |---|---|---|---|---|
  | ratio | 2.556 | 2.297 | 2.200 | 2.150 |

  It tends to 2 only asymptotically.
* **Correction.** "Parity raises the pigeonhole lower bound on `sum_(x<y) D_xy` additively by at most
  `sum_(m_Q < 2M) log q (M^2/(2 m_Q) + M/2) = (1/2 + o(1)) M^2 log M`. That is, at most a factor
  `2 + O(1/log M)` in aggregate, and the gain is positive only through parity bias."

The §6 constants `(4 - o(1))` and `(8 - o(1))` agree with the computation (3.21 and 6.90 at `M = 10^5`).
They converge slowly, at rate `O(1/log M)`.

### 2.3 Issue 3: factor-2 slip for `D_max(D_A)` (§6, row "Prop. B.1")

`sharp.md` Lemma 3.3(1) bounds `log m_(xy,1) <= mu_M log M` with `mu_M = 10 + 26/loglog M`. Hence
`D_max = 2 log m <= (20 + o(1)) log M`, not `(10 + o(1)) log M`. Prop. B.1 then gives at most
`(40 + o(1)) M log M/K`, not `(20 + o(1))`.

The author's own §2 table (`D_max/log M = 12.4` at `M = 1024`) and mine (10.0-13.6, §3.2) show the
quantity is above 10 at accessible `M`. The consistency conclusion (below the fakes'
`W >= 66(M-1) log M`) survives.

### 2.4 Issue 4: Theorem C assumes its growth; `H` is stronger than used

The proof of Theorem C uses `H` at **one** `k = floor((2 - eps/4) L/l)`, and only in the form
`F_k >= 1`, i.e.

```text
(H_1)  some pair (x, y) has unit-1 collisions at >= (2 - eps/4) L/l primes of (M^(1-eta), M].
```

Then `D_max >= (4 - eps/2)(1 - eta) L^2/l - O(L)`, and B.1 finishes.

So Theorem C is exactly "B.1 + (`H_1`)". It carries no arithmetic content beyond the hypothesis, and the
`log M/loglog M` gain is assumed in `H_1`.

* The stated `H(eta, delta)` asks for `F_k >= delta^k E_k` for **every** admissible `k`. That is much
  stronger than `H_1`. The open problem `(H-actual)` in §7 should be posed for `H_1`, the weakest
  hypothesis that suffices.
* The wording "the hypothesis the task asked for" and the §0 table status "proved, conditional on `H`"
  should add "and on `K`-uniform margins and `lambda_Q >= lambda_0`". The table row lists them, but the
  status cell does not.
* This is not an error. It is a scope warning: Theorem C does not *reduce* the growth question to an
  easier one. It only says that uniform margins convert one heavy pair into `W`.

### 2.5 Issue 5: the `Im(.)` description of `H` for actual clusters

For an actual cluster with `q` prime to `N`, put `g = gcd(z_x, z_y)`, `z_x/g = u B`,
`z_y/g = u' conj(B)`, and `w = z_x conj(z_y)/|g|^2 = u conj(u') B^2`. Then:

```text
q^a | (z_x - z_y)/g   <=>   w = |B|^2  mod q^a   <=>   q^a | Im w   and   Re w = +|B|^2 mod q^a.
```

`q^a | Im w` alone gives `Re w = +-|B|^2` (since `Re^2 + Im^2 = |B|^4`). So it characterises
collisions with unit `+-1`, not unit 1. This was checked in exact arithmetic (`parity_actual.py`):
* the three-way equivalence holds in all 215,670 pair-modulus instances;
* `Im`-divisibility without a unit-1 collision occurred 10,892 times.

So "Equivalently it concerns the smooth parts of `Im(z_x conj(z_y) u)/|g_xy|^2`-type integers" should
read "it concerns the joint divisibility of the reduced chords `(z_x - z_y)/g_xy`, i.e. the congruences
`z_x conj(z_y)/|g|^2 = |B|^2 mod q`".

### 2.6 Issue 6: the independent model of Lemma C.3 is not a P3 model

* Lemma C.3 is correct for classes uniform on `Z/m_q`.
* In P3-char and smaller classes, the class of row `x` is forced to have parity `sigma_q(x)` (Prop. B.3).
  So that model is a member-level model only for `A_X`.
* Consider the parity-respecting version: classes uniform among the `m_q/2` classes of the forced
  parity.
  * Pairwise independence survives. With all three parities equal,
    `P(kappa_x = kappa_y = kappa_z) = (2/m)^2`; otherwise both sides are 0.
  * With `sigma` *fixed*, however, `E F_k = sum_P e_k((2/m_q) [sigma_q(x) = sigma_q(y)])_q`. This is not
    `E_k` in general.
  * It equals `E_k` only on average over independent uniform `sigma`, since each `k`-set of primes
    survives with probability `2^(-k)` at doubled weight.
  * So for P3-gen data, `H` also needs a non-degeneracy statement about the Legendre-parity pattern.
* Remark 1 on Theorem C ("Without uniform margins, `H` gives nothing beyond Corollary A.1 at the level
  of pair-cut certificates") also needs a demand matrix that satisfies `H` *and* has
  `lambda_max(D) = O(M log M)`. The write-up proves `H` for the independent model but not the
  `lambda_max` bound. That bound is only measured (1.91-1.93 `M log M`; rechecked, §3.7).
* **Correction:** label Remark 1 as evidence-backed, or add the easy proof (`E D_x ~ 2 M log X` plus
  concentration).

### 2.7 Issue 7: "adapted profiles" in P3 (scope of Cor. A.2)

Theorem A fixes `D` and varies the cut measure. That is legitimate for `A_X`, where point residues are
free up to the norm. In P3-char the point classes are `kappa = A ell`, with forced parity
`sigma_q(x) = sum_j a_xj nu_j(q) + lambda_q k_x`. In P3-gen and in the coherent designs of
`sharp.md` §6.6 they are classes of genuine integers built *from the profile*. Changing the profile
changes `D`.

Hence:
* §0 item 3 and §6 row "sharp.md §6.6" ("At the real-weight level it disappears for profiles adapted
  to the demand pattern") describe the relaxation, not a member. The figure `~ 2 pi M log X` there is
  evidence.
* Cor. A.2 is still a correct statement about certificates that see only `D`. Cor. A.5 is correct
  because it exhibits actual members.
* **Correction:** "In the real-weight relaxation with `D` held fixed, the heavy pairs cost at most
  `pi lambda_max(D)`. Whether a P3 or P3-gen member can realise an adapted profile *together with* its
  own demand matrix is open (it is not the same as Question Q)."

### 2.8 Issue 8: heights "subsumed by the count of primes" is overstated

Prop. B.5 proves it for one specific pigeonhole: two-term monomials `+-e_j +- e_l` at single moduli.
Not covered:
* simultaneous (CRT) pigeonhole over several moduli;
* longer monomials;
* the link between the forced parity `ell_j = nu_j(q) mod 2` and the size of `p_j` (both are functions
  of the actual prime).

The sentences "The analogous count with monomials of norm `<= k` forces only
`max_j w_j >= log(2r/k) - O(1)`" and "Pigeonhole consequences of the heights are subsumed by the
trivial count of primes" (also in the mechanism table and §7 item 4) are unproved generalisations.

**Correction:** "the single-modulus two-term height pigeonhole (Prop. B.5) is subsumed by the count of
primes."

### 2.9 Issue 9: Remark B.4

* "Since `alpha >= 2^(-|S|)`, the largest forced CRT modulus grows by a factor `alpha 2^(|S|) <= 2^(|S|)`".
  The inequality needs `alpha <= 1`, not `alpha >= 2^(-|S|)`.
* "So pigeonhole over many moduli, with the forced parities, gives `D_max = (2 + o(1)) log M` at best"
  is a statement about this CRT argument. It does not bound what any argument can force.
  `D_A` shows only `D_max <= (20 + o(1)) log M` is attainable in P3-char.

### 2.10 Issue 10: evidence ranges are family-specific

* **LP ratio.** `cutlp.py`, my own LP, rechecked the sandwich on 36 matrices at `M = 4..9`. It also
  searched adversarially at `M = 5..7`.
  * `W_cut/(2 lambda_max)` reached 1.443, on a random residue demand at `M = 4`, `X = 12`.
  * Hill-climbing reached 1.500.
  * The author's "`[1.000, 1.224]`" is correct for his 48 matrices, but it is not a property of
    residue demands. Nothing depends on it, and the proved range is `[1, pi/2]`.
* **Colouring.** The §6 value "`2 lambda_max(D) ~ 6.4 M log M` for `D_A`" holds for random parity
  colourings.
  * For the constant colouring (`sigma_q = 0` for all `q`, which Lemma A.4 also covers) I measure
    `lambda_max/(M log M) = 5.56, 5.91, 6.17, 6.36, 6.52` at `M = 64..1024`. That is
    `2 lambda_max ~ 13 M log M`.
  * This is well inside the proved `72 M log M + 222 M`.

### 2.11 Issue 11: update and attribution

* **Theorem L.** Remark 2 on Theorem C and §7 item 2 say `flipped.md` Theorem L "already penalises
  actual Hadamard-core clusters". The later `round8/lw.md` (Thm 5.1, Cor. 5.2; refereed in
  `round8/referee_lw.md`) shows that the Paley skew double core of `round5/exact.md` is immune to every
  single LW form. That core is a Hadamard-core, Plotkin-optimal, near-uniform-margin profile. The
  sentence should say "Theorem L penalises Hadamard-core profiles with flip mass `F < W/2`; the skew
  double core is immune".
* **Theorem A(3).** It is Goemans-Williamson random-hyperplane rounding applied to `rho = I - B/Lambda`.
  The resulting `pi/2` gap between the PSD (Perron) bound and the cut value for nonnegative objectives
  is Nesterov's `2/pi` theorem in LP-dual form (Yu. Nesterov, *Semidefinite relaxation and nonconvex
  quadratic optimization*, 1998).
* **Cut cone.** The cut cone, the hypermetric inequalities and the LP over cut measures are standard
  (M. Deza and M. Laurent, *Geometry of Cuts and Metrics*, 1997).
* **Novelty.** The application to margin profiles is the new part. The sandwich itself is classical and
  should be cited as such.

---------------------------------------------------------------------------------------------

## 3. Independent computations (`round7/referee_p3lower_checks/`)

All scripts are new. LP results are floating-point (HiGHS) and are evidence. The claims marked
"exact" use integer, rational or `Q(sqrt 5)` arithmetic. No conclusion of this report rests on a
floating-point value alone, except where it is labelled evidence.

### 3.1 Theorem A: cut LP, sandwich, duality (`cutlp.py`, output `cutlp_out.txt`)

* **Families.** All cuts were enumerated for `M = 4..9`, and `W_cut(B)` and its dual were solved
  separately. The 36 matrices are:
  * uniform;
  * heavy pair;
  * random exponential;
  * sparse;
  * independent residue demands at `X = 3M` and at `X = M^2`, with random lifts at prime powers.
* **Results.**
  * Primal = dual in every case.
  * `2 lambda_max <= W_cut <= pi lambda_max` in every case.
  * The cut identity holds at the optimum.
  * Optimal supports range from 3 to `binom(M, 2)`.
* **Uniform demand.** For `B = J - I`, `W_cut = 2M` for odd `M` (the hypermetric bound with `b = 1` is
  tight) and `2(M-1)` for even `M` (Hadamard).
* **Dual certificates.** The Perron and hypermetric certificates were checked feasible over all cuts
  (600 random instances).
* **Adversarial search.** Hill-climbing on `W_cut/(2 lambda_max)` at `M = 5, 6, 7` reached 1.500, 1.500
  and 1.499. That is below `pi/2 = 1.571`, consistent with Theorem A(3).
* ALL THEOREM-A CHECKS PASSED.

### 3.2 Lemma A.4 on the design `D_A` (`design_DA.py`, `design_DA_out.txt`)

**Construction, from `sharp.md` Def. 3.2 literally.**
* The map `p'(q)` of Lemma 3.1. Asserted for all `q < 8M`: `2p'(q) <= m_q` and `q < 8p'(q)`.
* Level-1 classes `sigma_q + 2(x mod p'(q))` in `Z/m_q`.
* CRT levels up to `a*(q)`. Injectivity at `a*` was asserted.
* Unit-1 pair moduli, then `D_xy = 2 log m_xy`.

| `M` | colouring | `D_max/log M` | `max_x D_x/(M log M)` | `lambda_max/(M log M)` | `max_x D_x` | `4 sum(Lambda' + log d)` | `36 M log M + 111 M` |
|---|---|---|---|---|---|---|---|
| 64 | constant | 11.47 | 6.37 | 5.56 | 1,696 | 4,037 | 16,686 |
| 64 | random | 10.02 | 3.47 | 2.76 | 924 | 4,037 | 16,686 |
| 256 | constant | 11.89 | 6.85 | 6.17 | 9,728 | 23,248 | 79,520 |
| 256 | random | 10.98 | 3.70 | 3.10 | 5,255 | 23,248 | 79,520 |
| 1024 | constant | 13.21 | 7.08 | 6.52 | 50,267 | 121,084 | 369,186 |
| 1024 | random | 11.88 | 3.63 | 3.26 | 25,734 | 121,084 | 369,186 |

Lemma A.4's chain `lambda_max <= max_x D_x <= 4 sum(Lambda'(d) + log d) <= 36 M log M + 111 M` holds in
every case (also at `M = 128, 512`). The per-pair inequality `log m_xy <= Lambda'(d) + log d` was
checked for every `d`. The window counts satisfy `max c_xy <= 8` for `eta = 1/2` (Lemma C.2 allows 36).
ALL D_A CHECKS PASSED.

### 3.3 Prop. B.3 and the reduced-chord criterion on actual Gaussian integers (`parity_actual.py`)

**Points.** Exact Gaussian arithmetic, with
`z_x = i^(u_x) prod_j pi_j^(a_xj) conj(pi_j)^(e_j - a_xj)`:
* 5-7 split primes from `{5, ..., 113}`;
* `e_j in {1, 2}`;
* 30-60 distinct random rows and random windings.

**Moduli.** 143 odd prime powers `Q = q^a <= 400`, with `q < 60` prime to `N`. In total 215,670
pair-modulus instances.

**Checked.**
* `t = z_x/z_y mod Q` has norm 1, and `t^(m_Q/2) = (-1)^(sum_j (a_xj - a_yj) nu_j(q) + lambda_q (u_x - u_y))`.
  This is the forced-parity formula, at every level.
* Every one of the 10,840 unit-1 collisions joins equal parities.
* The collision counts are `>= E(n_0, m/2) + E(n_1, m/2)`.
* Unit-1 collision `<=>` `Q` divides the reduced chord `<=>` `z_x conj(z_y)/|g|^2 = |B|^2 mod Q`.
* `Q | Im(z_x conj z_y)/|g|^2` is equivalent to a collision with unit `+-1`. There were 10,892 instances
  with unit `-1` only (issue 5).

ALL ACTUAL-POINT PARITY CHECKS PASSED.

### 3.4 `E(n, Q)` facts and the factor-2 claims (`e_counts.py`)

* **Exact checks.**
  * Closed form against brute force (`Q <= 5`, `n <= 10`).
  * Bounds `n^2/(2Q) - n/2 <= E <= n^2/(2Q)` (`Q <= 12`, `n < 60`).
  * Superadditivity (`n_0, n_1 < 40`).
  * The B.3(3) chain and the `beta` lower bound (`m <= 24`).

  All OK.
* **Per-modulus and aggregated ratios.** As in §2.2.

### 3.5 Hypothesis `H`, Lemma C.1 and the Theorem C threshold (`h_stats.py`)

* **Independent model** (window `(sqrt M, M]`, 20 trials), mean `F_k/E_k` (evidence):
  * at `M = 128`: `1.001, 0.993, 0.975, 0.966` for `k = 1..4`;
  * at `M = 512`: `1.000, 1.001, 1.001, 0.996, 0.953, 0.897` for `k = 1..6`.

  This is consistent with `E F_k = E_k`.
* **Lemma C.1** holds for the exact window `e_k`.
* **Threshold.** With exact window primes, `eta = delta = eps = 1/2` and
  `k = floor((2 - eps/4) L/l)`, the Lemma C.1 lower bound for `log(delta^k E_k)` is:

  | `M` | `10^4` | `10^5` | `10^6` | `10^7` |
  |---|---|---|---|---|
  | `k` | 7 | 8 | 9 | 10 |
  | lower bound for `log(delta^k E_k)` | 0.89 | 2.83 | 4.43 | 5.78 |

  So Theorem C's `H` is used at accessible `M`. The guaranteed `D_max` is only about
  `0.6 L^2` (`2(1-eta) k L`), and this coefficient decreases like `1/l`.

### 3.6 Prop. D.1 (`clique_D1.py`, exact integers)

* **Cases.** Ten Sylvester cases `(g, t, L)` from `(2, 2, 1)` to `(32, 4, 97)`, `M` up to 128.
* **Checked.** All Gram entries `-L` / `-1` / `r`, nonconstant columns, and `r = lambda_max(k)`. For
  example, `M = 128`, `t = 4`, `L = 97` gives `r = 415`, against 12,319 for uniform margins.
* ALL D.1 CHECKS PASSED.

### 3.7 Independent residue spectra (`indep_lmax.py`, evidence)

At `X = 3M`, with random lifts at prime powers:

| `M` | `lambda_max/(M log M)` | `max_x D_x/(M log M)` | `D_max/log M` |
|---|---|---|---|
| 64 | 1.909 | 2.43 | 8.7 |
| 256 | 1.916 | 2.28 | 13.2 |
| 1024 | 1.927 | 2.08 | 15.1 |

This agrees with the author's §2 table row (1.905-1.926).

### 3.8 Horn counterexample (`psd_horn.py`)

As in §2.1. The eigenvalues of `Z` are exact, and `<Z, B> = 6/5` is exact. `W_cut(B) = 10/3` (LP).
COUNTEREXAMPLE CONFIRMED.

---------------------------------------------------------------------------------------------

## 4. Class check: are the inputs valid for actual circles?

| statement | class claimed | inputs used | valid for every actual normalised cluster? |
|---|---|---|---|
| Lemma 1.3, 1.4, Cor. A.1 | `A_X` and below | pair inequality (from (1.1) of `lower.md`, proved for actual clusters in `lower.md` Prop. 1.2); profile normalisation | **yes** |
| Thm A, Cor. A.2, A.3 | relaxation | none about circles (statements about LPs or certificates) | n/a. The barrier concerns certificates that see only `D` (issue 7) |
| Cor. A.5 | P3-all at `X = 3M^A` | `sharp.md` Thm 6.4 members (refereed) | n/a (barrier). The members exist, and their `D` is `D_A`'s |
| B.1, B.2 | `A_X` | cut identity, pair inequality, plus a hypothesis on the profile | yes, conditional on the stated profile hypotheses |
| B.3 | P3-char and below | forced Legendre parity | **yes**. Actual residues factor through `pi_j`, and the parity formula is verified exactly on actual points (§3.3) |
| B.5 | P3-all | heights (`sharp.md` Prop. 6.1) | yes (heights hold for actual clusters) |
| Thm C, Cor. C.4 | `A_X` with `K`-uniformity, `lambda_Q >= lambda_0`, `H` | B.1, plus a hypothesis on window collisions | yes, conditional on three hypotheses, none known for any actual cluster |
| Lemma C.2 | the `D_A` design | `sharp.md` Def. 3.2, Lemma 3.1 | n/a (concerns fakes) |
| Lemma C.3 | independent uniform classes | a model | **not** a P3 model (no forced parity; issue 6) |
| D.0, D.1 | profiles with distinct primes; `+-1` matrices | PNT in progressions; Hadamard matrices | D.0 yes. D.1 concerns equal-weight sign matrices, and no circle is claimed |

No statement that is labelled "proved" uses an input that fails for actual circles. No class used for
a valid inequality omits actual clusters: Prop. 1.2 puts them in P3-gen at every `X`.

---------------------------------------------------------------------------------------------

## 5. Corrected statement

**Round 7, p3lower (corrected).** All statements are about the growth question, with constants
irrelevant.

1. **(Cut LP; classical tools, new application.)** For every member of `A_X(M, C)`, hence every
   normalised actual cluster, every member of P3-char, P3-all and P3-gen, and every `X`:

   ```text
   W  >=  W_cut(D - c_C(J - I))  >=  2 lambda_max(D) - 2(M-1) c_C^+ .
   ```

   For `B >= 0`, `2 lambda_max(B) <= W_cut(B) <= pi lambda_max(B)`. The upper bound is
   Goemans-Williamson rounding / Nesterov's `pi/2` theorem, and the lower bounds are Perron and
   hypermetric dual certificates.
   * A certificate that sees only the pair demand matrix `D` (through the pair inequalities and the
     cut structure with real weights) cannot exceed `pi(lambda_max(D) + (M-1) c_C^-)`.
   * PSD arguments give exactly `2 lambda_max(D + |c_C|(J-I))` when `C <= sqrt 2`. In general they
     give a DNN value within `2(M-1)|c_C|` of `2 lambda_max(D)`. They do *not* always reduce to the
     `u >= 0` Perron bound (Horn counterexample at `M = 5`).
2. **(Barrier for P3-all; proved.)** The members of `sharp.md` Thm 6.4 (P3-all, `X = 3M^A`) have
   `lambda_max(D) <= max_x D_x <= 36 M log M + 111 M`. So no certificate of the type in item 1 proves
   `W_min^all >= f(M) M log M` with `f -> infinity`.
3. **(Elementary restricted bounds; proved.)**
   * `K`-uniform margins give `W >= (2(M-1)/K)(D_max - c_C)`.
   * Column-balanced profiles give `W >= 2 max_x D_x - 2(M-1) c_C`.
   * In P3-char and below, unit-1 collisions respect the forced Legendre parity. The parity-refined
     pigeonhole raises the lower bound on `sum D_xy` **additively** by at most
     `(1/2 + o(1)) M^2 log M`: a factor `<= 2 + O(1/log M)` in aggregate, and unbounded for single
     moduli with `M <= m_Q < 2M`. With all `p_j > X >= M` it gives `W >= (4 - o(1)) M log M`, and
     `(8 - o(1)) M log M` under full parity bias.
   * The single-modulus two-term height pigeonhole gives only `W >= r log min(r, X) - O(r)`, which the
     count of distinct primes already gives. Other uses of the heights are not addressed.
4. **(Conditional; proved.)** Suppose a member of `A_X(M, C)` (in particular an actual cluster), with
   `X >= M`, has:
   * `K`-uniform margins;
   * `lambda_Q >= lambda_0`;
   * one pair with unit-1 collisions at `>= (2 - eps/4) log M/loglog M` primes of `(M^(1-eta), M]`.
     This is `H_1`, which is implied by the author's `H(eta, delta)`.

   Then `W >= ((8(1-eta) - eps)/K) M log^2 M/loglog M`, and so
   `M <= (K/(4(1-eta)) + o(1)) log R logloglog R/(loglog R)^2`. This is B.1 applied to the assumed
   heavy pair, and the growth is in the hypothesis.
   * The round-6 design `D_A` violates `H` (and `H_1` at fixed `k`) for large `M`.
   * Independent uniform classes (an `A_X` model without forced parity) satisfy `H` with probability
     `-> 1`.
   * No actual cluster is known to satisfy `K`-uniformity or `H_1`.
5. **(Discreteness; proved parts.)** Distinct primes force `r = O(M)` threshold cuts at
   `W = O(M log M)`. For equal weights, a sign matrix with correlations `<= -k_xy` needs
   `r >= lambda_max(k)`. Clique-partition patterns attain this exactly. Whether prime weights and
   `O(M)` cuts can absorb the *sum over moduli* of such patterns (Question Q) is open. The small-`M`
   data are evidence only.
6. **Not proved** (interpretive or evidence-backed):
   * that `lambda_max(D) = O(M log X)` for random-like or actual demands;
   * that heavy pairs are "absorbed" by profiles of an actual P3 or P3-gen member (the demand matrix
     depends on the profile there);
   * that `H` gives nothing beyond Cor. A.1 without uniform margins;
   * that parity gains at most a factor 2 per modulus (false);
   * that all height pigeonholes are subsumed;
   * that Theorem L penalises all Hadamard-core clusters (false for the skew double core;
     `round8/lw.md`).

**Status.** The target `M = o(log R/loglog R)` is not reached, as the write-up says.

---------------------------------------------------------------------------------------------

## 6. Recommendation

**Accept after minor revision.**
* Make the eleven textual corrections above. Issues 1-3 change sentences, not theorems.
* Restate Theorem C with the weakest sufficient hypothesis `H_1`, and pose `(H-actual)` for `H_1`.
* Add the citations for Theorem A.
* Add a sentence that Cor. A.2 concerns certificates that see `D` only, and that in P3 classes `D`
  depends on the profile.

The core contributions survive refereeing:
* the cut-LP sandwich and its corollaries as valid inequalities for actual clusters;
* the P3-all barrier for pair-only certificates (Cor. A.5);
* the parity and height pigeonhole bounds;
* the conditional Theorem C.
