# Referee report on `round4/hardness.md` (barrier theorems: why the uniform `C sqrt R` bound is hard)

Referee: adversarial pass, round 4. Inputs read: `round4_context.md`; `round4/hardness.md`; its three
check scripts and logs; `docs/claude_review/growth/{conditional,referee_conditional,ptolemy,families,
referee_families,walsh-extraction}.md`; and the research notes `endpoint_recent_literature_check.md`,
`affine_cube_walsh_rigidity.md` and `codex_uniform_bound_research.md` (Coppersmith mention).
The referee's scripts and logs are in `round4/referee_hardness_checks/`.

## 0. Verdict

**The result survives in a weaker form.** Every numbered Theorem, Proposition, Lemma and Corollary
in the body is correct as stated, or becomes correct after a one-line repair (listed in §2). I found
no fatal error. Part of the explanatory layer, however, is presented as proved, or as a consequence
of proved results, when it is not:

* the summary,
* the "What this shows" paragraphs,
* Section 7.

Two of these sentences are false as worded (§1, items M1 and M4). Specifically:

* "Single linear forms carry nothing" and "every single linear form … is controlled exactly, up to a
  constant, by the elementary bound" are **false** (check R12). What is proved is optimality *as a
  function of `n` alone*, and only along the sequence `n = A^2 + (A-1)^2`.
* "Any content of U lies in entangled systems of Liouville-extremal forms" and "such a construction
  must produce clusters whose pair forms are within bounded factors of the elementary bound" are
  **unproved**. The record cluster itself has 20 of its 45 pairs at ratio above 1.5, and 3 above
  100 (check R11).
* "For every modulus **at or above** the Gaussian Lenstra exponent 1/4, each residue class holds at
  most one cluster point" is **false** in two places:
  - at `alpha = 1/4` exactly (check R6);
  - at finite `N` for `1/4 < alpha < 1/4 + log C/log N`, where 22 of the 45 pairs of the record
    cluster share a residue class (check R1).

  The body's Prop. 6.1 has the correct quantifier.
* "This closes the Cilleruelo–Granville link to Rudin's conjecture" claims something about a link
  the author states they never read (§9 of the write-up). Theorem 2 covers only consequences at the
  level of energy (second moment) and sumsets.
* Two items are heuristic dictionary statements, not theorems, yet appear under "What is proved":
  - "On the circle the margin is identically zero" (the Gaussian Coppersmith lemma is "not proved
    here, and not used");
  - "the whole circle problem corresponds to the residual near-balanced range of the split
    problem". Prop. 4.1 proves that *no* homomorphism dictionary exists, so it cannot be the
    source of a proved correspondence.
* `ptolemy.md` Thm D gives `U`, not `U_abs` (its bound is `4 ceil(C)` times an absolute constant).
  So "Vojta on `Mbar_{0,8}` implies `U_abs`" (overview §0.1 and §7 item 3) is a mis-citation.

The headline negative finding stands: no implication from `U`, or from a growth improvement, to a
recognised open problem was found. That is a search report, and it is stated as one.

---------------------------------------------------------------------------------------------

## 1. Major issues (claims stated as proved, or as following from proofs, that are not)

**M1. "Single linear forms carry nothing" (overview §0.3, §4 "What this shows", summary item 3).**

What is proved (Prop. 3.2(1)–(3)):

* every pair satisfies `|lambda| >= sqrt2 n^(-1/2)`;
* equality holds up to `1 + O(1/n)` along the pairs `(g u, g i ubar)` with `u = A + (A-1)i`.

So no lower bound **depending on `n` alone** can beat the elementary one. That is all.

The text goes further. It says that "every single linear form attached to two points of a circle is
controlled exactly, up to a constant, by the elementary bound. No transcendence statement, proved or
conjectured, can say more about it", and that any such consequence "is either already elementary or
false". This is wrong.

* **The ratio can be large.** The ratio `|lambda|/(sqrt2 n^(-1/2))` reaches 763 in the write-up's
  own cluster. Its a-priori cap is `C sqrt(n/(2R))`, which is about `1079` there.
* **Prime powers.** For `n = 5^a` (with `u = (2+i)^a` and `gcd(u, ubar) = 1`), Baker's theorem forces
  the ratio to infinity. Exactly: the best pair chord over `51 <= a <= 100` gives ratio `>= 10^16.9`,
  and over `201 <= a <= 400` it gives ratio `>= 10^69.8` (R12).
* **A third option.** There are therefore single-form statements that are true, non-elementary and
  theorems (Baker on prime-power pair norms). The write-up's own `auxiliary.md` citation ("Baker
  gains on non-saturated perfect powers") says the same.
* **What survives.** Matveev is weaker than Liouville when all exponents are `<= 10^7` (verified:
  `c(k) >= 7.48*10^5`, threshold `a_max <= 1.19*10^7`). Lang–Waldschmidt-shaped bounds do not
  improve the exponent when all exponents are `<= 2`. These are correct and should be the stated
  content.

**M2. Section 7, second bullet, and summary: "any content of U lies in entangled systems of
Liouville-extremal forms" / "must produce clusters whose pair forms are within bounded factors of
the elementary bound … `Theta(M^2)` simultaneously near-extremal approximations".**

This is introduced by "By the results above, such a construction must have all of the following
properties", but nothing in §§2–6 proves it.

* Prop. 3.2(2) bounds each pair's ratio only by `C sqrt(n/(2R))`, which is unbounded once `n >> R`.
  Coprime pairs have `n = N = R^2`.
* Lemma 3.3 excludes only exact **products** with at least 3 independent blocks.
* **Data (R11):** in the record 10-point cluster, 25 of 45 pairs have ratio `<= 1.5`, 20 have ratio
  `> 1.5`, and 3 have ratio `> 100`.

"Entangled (non-product)" is supported by Lemma 3.3. "Liouville-extremal" and "`Theta(M^2)`
simultaneously near-extremal" are not supported. They should be marked as heuristics.

**M3. The Rudin barrier is stated more broadly than Theorem 2 supports (overview §0.2, Prop. 2.3
"What this shows" item 4, summary).**

Theorem 2 and Corollary 2.2 are correct (re-derived, and checked exhaustively in R2–R4). They show
that every consequence of `U` that passes through **additive energy or sumset size** is
unconditional. Two stronger statements are made that this does not support.

* **"No statement … that follows from U through representation counts is open."** The `L^infinity`
  bound also gives higher moments, `sum_n r(n)^k <= g(c)^(k-1) |A||B|` for `k >= 3`. Energy does not
  give these. Unconditionally one has only `max r * E`, which carries the `log R/loglog R` loss.
  Read literally, the sentence is also false, because `U` restated as a "squares in two
  progressions" statement is itself such a statement.
* **"This closes the Cilleruelo–Granville link to Rudin's conjecture as a hardness route at exponent
  1/2."** The author says (§9) that the paper could not be read, "so I make no claim about which of
  their implications coincide with Corollary 2.2". The overview sentence contradicts §9.
  - The research notes (`endpoint_recent_literature_check.md`) already record the endpoint
    difference-multiplicity bound.
  - They also record that `Lambda(4)` and energy control at the endpoint are elementary, and the
    abstract `L^2`/`L^infinity` gap.

  Theorem 2 is an immediate corollary of those notes. The write-up acknowledges this for Lemma 2.1.

A correct statement is: "every consequence of `U` for squares from square-root intervals that
factors through the second moment (energy, `Lambda(4)`, sumset size, squares in far-regime
progressions) is unconditional; `U` is equivalent to the two-interval `L^infinity` statement, whose
one-interval case is trivial."

**M4. Heuristics presented under "What is proved".**

* **The circle margin.** "On the circle the margin is identically zero, because `|z| = |zbar|`" is
  derived in Remark 4.6(2) from a Gaussian analogue of Lemma 4.4 that is "not proved here, and not
  used". The supporting refereed statement is `threshold/construct.md` (exact criticality of
  forced-divisibility Liouville arguments). That is a different, though related, statement.
* **The split dictionary.** "The whole circle problem corresponds to the residual near-balanced
  range of the split problem" (Remark 4.8) rests on a dictionary that Prop. 4.1 proves cannot be a
  torus homomorphism. It is an analogy.
* **"Every standard input below Vojta's conjecture is provably elementary, vacuous or exactly
  critical there."** This is a claim about an open-ended class. What is proved is that the eight
  listed inputs fail, and the Coppersmith item among them is unproved on the circle.

Both of the first two statements are reasonable remarks. They are not theorems, and the summary
lists them as "proved".

## 2. Minor issues (each has a one-line repair)

**m1. Residue-class quantifier (overview §0.5, summary item 5).** The phrases "at or above the
critical exponent 1/4" and "for every modulus in or above the Gaussian Lenstra range" should read:
for fixed `alpha > 1/4` and `N > C^(1/(alpha - 1/4))`. The body's Prop. 6.1 is correct.

* **At `alpha = 1/4` exactly (R6).** Take `z = u^2` and `w = i u ubar`, with `u = A + (A-1)i` and
  `N = Norm(u)^2`. Both points lie in the class `0 mod u`. Here `Norm(u) = Norm(N)^(1/4)` and
  `|z - w| = sqrt2 N^(1/4)`. Verified exactly for `A = 2, 3, 10, 101, 1000`.
* **At finite `N` (R1).** In the record cluster, 22 of 45 pairs share the class `0` modulo their
  gcd `g`, where `Norm(g) = Norm(N)^alpha` with `alpha` up to `0.3322 > 1/4`. This is consistent
  with Prop. 6.1, whose threshold here is `1/4 + log C/log N = 0.3510`.

**m2. Data E misidentifies "the eight closest pairs."**

* `check_forms.py` sorts the pairs by **pair norm**. The pairs with norms 1105, 1625 and 4901 are
  among the *farthest*: their angles of `0.020–0.043` are near the full span `0.0445`.
* The true eight closest pairs (R1) have pair norms 47074105, 1065025 (×2), 240125 (×2),
  95485 (×2) and 47125. Their ratios are `1.000` or `1.414`.
* The qualitative conclusion (Liouville attained at the closest pairs) still holds.

**m3. Mis-citation of `ptolemy.md`.** `ptolemy.md` Thm D (Vojta on `Mbar_{0,8}`) proves
`M_orig <= 4 ceil(C/C_0)` times an absolute bound. That is `U`, not `U_abs`. Only `conditional.md`
Thm 5.1 (on `Bl_Y X_7`) is known to give `U_abs` and hence (v). The places to fix are overview
§0.1, §2 Remark 2 and §7 item 3.

**m4. Cor. 4.5 proof.** It takes `P = D`, but Lemma 4.4 needs `P` to be an integer, and `D` is
real. Put `P = ceil(D)`; every divisor in the window is then `P + x` with `0 <= x <= X`. The rest of
the chain checks.

* R8 checks the chain at the hypothesis threshold for 11970 triples `(eta, C, beta)` and finds no
  failure. The smallest relative margin is `6*10^-4`, at `eta = 10^-3`.
* By hand: `w <= 3m <= 6/eta + 3`, `t = w - m >= 0`, and the lower bound
  `2 beta m/w - m(m+1)/(w(w+1)) > beta^2 (m-1)/(m+beta)` is correct.

**m5. Remark 4.6(1).** "It then gives `S <= 6 log n/log lambda + 4`" should read
`max(C + 1, 6 log n/log lambda + 4)`.

**m6. Remark 4.8.** "Lemma 4.2(a), Cor. 4.5 and Cor. 4.7 leave open `U^split` only in
`1 + n^(-eps) < lambda < n^eta`" treats Cor. 4.7 as unconditional, but it assumes #886.

* Unconditionally, the vertex trick covers only `lambda - 1 = O_C(n^(-1/4))`.
* The near-vertex range `n^(-1/4) << lambda - 1 <= n^(-eps)` is open unless #886 is assumed.
* Also, Cor. 4.5's bound depends on `eta`, so a uniform `K(C)` is not available as `eta -> 0`.

**m7. Lemma 4.2(b).** The "Consequently" clause applies the first part with
`psi = (kappa + C/2) R^(-1/2)`, which needs `psi <= pi/4`, that is,
`R >= (4(kappa + C/2)/pi)^2`. The hypothesis should be added.

* R5 checked all `N <= 3000`, five directions `w`, and grids of `kappa` and `C`: 134850 cases,
  36990 of them with `psi > pi/4`. It found no violation.
* So the clause is very likely true in general, but it is unproved outside the stated range.

**m8. Prop. 4.1, step 1.** "The `psi_j` span `2L`" is inaccurate: `psi_j = (2e_j, -1)` is not in
`2L`. The `psi_j` span a lattice that contains `2L` with index 2 (R10). The conclusion
`sigma = -1` on `L` follows directly from `sigma(u^a t^c) = u^(-a) t^(-c)`, using `sum a = -2c`.
Harmless.

**m9. Prop. 3.2(3), parenthetical.** The family `u = A + i` needs `A` even. For odd `A`, `1 + i`
divides `gcd(u, ubar)` (R7).

**m10. Overview §0.3, "Any Lang–Waldschmidt-shaped bound is weaker on squarefree pair norms."**
This needs `n >= (kappa_k/sqrt2)^(2/(1 + 2 eps))`. The body's item 5 states the condition correctly.

**m11. Prop. 2.3(2).** The proof gives `2(c^2/sqrt2 + 1)`. The stated `2(c^2/sqrt2 + 2)` is weaker
but correct. R3 checks the sharper constant for all `n <= 2*10^5` and `c` in `{1/2, 1, 2}`.

**m12. Dependencies.** "`A_5(kappa)` is infinite for every `kappa`" rests on the notes' asymptotic
proof for the cyclotomic family. `referee_families.md` explicitly did not re-audit that proof. The
conclusion "(v) needs `M >= 9`" does not need it: the golden family alone makes `A_8(5.9)` infinite
(confirmed, decreasing to 5.8977 from above), and subtuple projection then handles every `M <= 8`.

**m13. Novelty.**

* "Empty" in Thm 1(a) is already in the converse of `conditional.md` Prop. 3.2(ii). What is new is
  the finite form, the chord-to-arc Lemma 1.1, and part (b).
* Theorem 2 is an immediate consequence of facts already in `endpoint_recent_literature_check.md`.
* Lemma 3.3 for `m >= 5` blocks is a special case of `walsh-extraction.md` Thm B: a pure
  Hadamard core of order `2^m >= 28` with blocks on the singleton columns. For `m = 3, 4` it is new,
  though elementary.
* Prop. 4.1 is the textbook anisotropy of the norm-one torus.
* Lemma 4.4 and Cor. 4.5 are the standard Coppersmith/Howgrave-Graham bound for divisors near a
  known approximation. They are probably folklore, and no novelty is claimed in the literature
  sense.
* The split-twin framing (`U^split`, Cor. 4.7, Remark 4.8) appears in neither the notes nor
  `docs/claude_review`. I confirmed this by grep.

## 3. What was re-derived and confirmed

| item | status | how |
|---|---|---|
| Lemma 1.1 (chords to arcs) | correct | by hand; `arcsin x <= (pi/2) x` on `[0,1]` by convexity |
| Thm 1(a), all four implications | correct | by hand, including the small-`R` cases `R_P <= kappa^2/3` and the content reduction `C sqrt(R'/d) <= C sqrt R'` |
| Thm 1(b) | correct | by hand: `d <= C^2 R_P` from `max|z'_j - z'_k| >= 1` |
| Remark 1: (v) needs `M >= 9` | correct | golden family plus subtuple projection (m12) |
| Lemma 2.1, Thm 2(1)–(3), Cor. 2.2 | correct | by hand; R2 (all `X <= 3000`, five values of `c^2`), R4 (`q <= 30`, `K <= 40`, 3000 values of `a` from the threshold) |
| Prop. 2.3(1) (`U` iff the box statement) | correct | by hand (`R > 2c^2/3` for Lemma 1.1; at most two quadrants for `R > (8C/pi)^2`) |
| Lemma 3.1, Prop. 3.2(1)–(5) | correct as stated in the body | by hand; R9 (Matveev constant and the `0.8`-product inequality) |
| Lemma 3.3 | correct | by hand: the lifted-angle argument needs `phi < pi`, which holds; `m >= 3` gives `N <= (C^2/2)^6`; the hypotheses are vacuous when `C^2 < 2` |
| Prop. 4.1 | correct (m8 wording) | by hand: `X*(T)^Gamma = 0` and `2 phi* = 0` |
| Lemma 4.2(a) | correct | by hand: `s(d) = d + n/d` |
| Lemma 4.2(b) | correct for `psi <= pi/4` (m7) | by hand; R5 |
| Lemma 4.4 | correct | by hand: triangular determinant; Minkowski for the cube; `d^m \| h(x)` and `\|h(x)\| < d^m` |
| Cor. 4.5 | correct after m4 | by hand; R8 |
| Cor. 4.7 (conditional on #886) | correct | by hand |
| Prop. 6.1 (body) | correct | trivial; the summary's quantifier is wrong (m1) |
| Author's checks E, E', F, G, H, I, J and A–D | reproduced byte-identically | re-run from copies in `referee_hardness_checks/rerun/` |
| Problem statements #886 and #887, and the bound `1 + C^2` | consistent with erdosproblems.com | web-search snippets only (direct fetch blocked) |

## 4. Corrected statement

> **No implication from `U` (or from `M = o(log R/loglog R)`) to a recognised open problem was found**
> (search result). Proved barrier lemmas:
>
> 1. **Fixed variety.** `U` iff for every `kappa` some `A_M(kappa)` is empty, equivalently finite,
>    equivalently non-dense. `U_abs` iff, for one fixed `M`, `A_M(kappa)` is finite for every `kappa`.
>    The golden family forces `M >= 9`. The only known implicant of `U_abs` is Vojta (`D = 0`) on
>    `Bl_Y X_7`. Vojta on `Mbar_{0,8}` gives `U` only.
> 2. **Energy.** For squares from two intervals of length `<= c sqrt X`,
>    `E <= (2 + c^2)|A||B|` unconditionally. So every consequence of `U` that factors through
>    energy or sumset size (`Lambda(4)`, sumset lower bounds, squares in far-regime progressions,
>    e.g. `<= sqrt(6K)` when `a >= (qK/2)^(4/3)`) is already a theorem. `U` is equivalent to the
>    two-interval `L^infinity` statement, whose one-interval case is trivial. Consequences of higher
>    moments are not covered. The Cilleruelo–Granville link itself was not examined.
> 3. **Pair forms.** `|lambda| >= sqrt2 n^(-1/2)` for every pair, with equality up to `1 + O(1/n)`
>    along `n = A^2 + (A-1)^2`. So no bound in terms of `n` alone improves it. Matveev is weaker on
>    pair forms whose exponents are `<= 10^7`. Bounds of shape `prod M(alpha_p)^(-1-eps)` do not
>    improve the exponent when all exponents are `<= 2`. On individual forms, for example
>    prime-power pair norms, the true ratio can be arbitrarily large, and Baker's theorem can beat
>    Liouville there. Product clusters with at least 3 independent blocks force `N <= (C^2/2)^6`.
> 4. **Split twin.**
>    * `Hom_Q(T_M, split) = Hom_Q(split, T_M) = 0`.
>    * At the vertex, `S <= 1 + C^2`.
>    * For `D <= n^(1/2 - eta)` and `n >= n_0(C, eta)`, `S(n, D, C) <= max(C + 1, 6/eta + 4)`
>      (Coppersmith, with `P = ceil D`).
>    * Near the vertex, Ruzsa's conjecture #886 implies `S <= 2K(eps) + 1`.
>
>    That the circle has "zero margin" everywhere, and that it corresponds to the residual
>    near-balanced split range, are heuristics.
> 5. **Residue classes.** For fixed `alpha > 1/4` and `N > C^(1/(alpha - 1/4))`, every residue class
>    modulo `sigma` with `Norm(sigma) >= N^(2 alpha)` contains at most one point of a cluster on an
>    arc of length `C sqrt R`. This fails at `alpha = 1/4`, and for smaller `N`.

## 5. Referee checks (`round4/referee_hardness_checks/`)

All integer arithmetic is exact. Floating point is used only for angle comparisons, logarithms
and the parameter grid in R8.

| script | content | result |
|---|---|---|
| `ref_cluster.py` (R1) | 45 pairs of the record cluster, sorted by angle; residue classes modulo pair gcds | the true eight closest pairs differ from Data E (m2); 22 of 45 pairs share a class with `alpha > 1/4`, max `0.3322` (m1) |
| `ref_elementary.py` (R2) | Lemma 2.1 exhaustively, `X <= 3000` | no violation |
| (R3) | Prop. 2.3(2) with constant `2(c^2/sqrt2 + 1)`, `n <= 2*10^5` | no violation |
| (R4) | Cor. 2.2, `q <= 30`, `K <= 40` | no violation; worst ratio `0.577` |
| (R5) | Lemma 4.2(b) clause, including `psi > pi/4` | 134850 cases, no violation (the proof gap stays, m7) |
| (R6) | class sharing at `alpha = 1/4` | verified (m1) |
| (R7) | `u = A + i` needs `A` even | verified (m9) |
| `ref_coppersmith.py` (R8) | Cor. 4.5 inequality chain | 11970 cases, no failure |
| (R9) | Matveev `c(k) >= 7*10^5`; product inequality | `min c(k) = 7.48*10^5` at `k = 1`; minimum of `prod/sum` is `0.8` |
| (R10) | the lattice `span(psi)` | `psi_1` is not in `2L`; `2L` lies in the span (m8) |
| `ref_ratios.py` (R11) | Liouville ratios in the record cluster | 25/45 `<= 1.5`; 20/45 `> 1.5`; 3 `> 100` (M2) |
| `ref_primepower.py` (R12) | `u = (2+i)^a` | ratio `>= 10^16.9` for all `51 <= a <= 100`, and `>= 10^69.8` for `201 <= a <= 400` (M1) |
| `rerun/` | author's `check_energy.py`, `check_forms.py`, `check_split.py` | outputs byte-identical |

Literature caveat: as for the author, direct fetches of arXiv, UCL and erdosproblems.com were
refused by the proxy, so the statements of #886 and #887 were confirmed from search snippets only.
Search results describe the Cilleruelo–Granville paper (math/0608109) as posing more than twenty
conjectures related to Rudin's problem. I could not determine which of their links go through
energy.
