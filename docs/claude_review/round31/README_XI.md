# Round 31, Part XI: judge report on the cocycle on one fixed dyadic parent, and closing assessment

Source: Part XI of `astra_free_XI/answer.md` (heading "# Part XI. One fixed dyadic parent", line 2500, to the end,
line 2599), her summary `final11.txt`, her script `xi_circle_check.py` (read, not run or imported), and the question
`reply11.txt`. Record: `round31/README_X.md` (format), `expert_note/expert_note.md` §1 and §5 item 12,
`round31/ref8_extend.md` Prop 3.1, `growth/auxiliary.md` l.55-77, `round31/README_IXL.md`, `round31/README_XL.md`,
and `astra_free_X/reply10.txt` for (T1), (T3), (T5).

Two referees each took one block. **Both reports are present; neither failed.**

| key | report | scope |
|---|---|---|
| rough | `ref11_rough.md` | XI.1 (the dilated collision bound (XI.3)-(XI.5)), the application (XI.1), the exact rows of XI.3, the twin and function-field tests |
| cocycle | `ref11_cocycle.md` | XI.2 (canonical identity, rigidity (XI.6), the s/4 bound for configurations, residue pinning), plus its own computations on `q_200` |

The judge's scripts are in `judge_XI/` (own exact code; nothing of Astra's or the referees' imported):
* `j1_triple.py` -> `j1_triple.out`: rebuilds all `4 * 2^16` lattice points of `q_200` from scratch, locates the
  cocycle referee's triple, recomputes rows, units, `H_ij`, `Q_ij`, `T_ij`, `t_ij`, the minimal `C`, and the
  canonical-residue test; brute-force class counts modulo `l^{a+1}` for `l <= 13`, `a <= 2`.
* `j2_dmin.py` -> `j2_dmin.out`: exact least chord on `q_200` (both families).
* `j3_small.py`, `j4_gen.py` -> `j3_small.out`, `j4_gen.out`: exhaustive integral-cofactor configurations in one
  unit class on `q_30`, `q_50`, `q_100` at `C = 8`, against genuine cluster sizes.

**Every claim labelled [J] is the judge's own and has not been refereed.**

Definitions. `q` is an odd squarefree product of `s` primes `p = 1 mod 4`, `W = log q`. For each `p | q` a Gaussian
prime `pi_p` of norm `p` is fixed. A *canonical mark* is `A = eta prod_p f_p` with `f_p in {pi_p, conj pi_p}` (its
*row*) and `eta` a unit; marks have norm exactly `q`. For two marks `i, j`: `H_ij` is the product of the common
factors `f_ip = f_jp`, `Q_ij` the product of the primes where the rows differ (so `|H_ij|^2 = q/Q_ij`), `X_ij` the
product of the factors of row `i` where they differ, and `T_ij = -Im X_ij`. `q_P` is the product of all primes
`p = 1 mod 4` (family 1) or all primes `p = 5 mod 8` (family 2) in `[P, 2P]`. `M_max(q, C)` is the largest number of
lattice points of `|z|^2 = q` on an arc of length `C q^{1/4}`. An **integral-cofactor configuration** at constant
`C` (the object of (A2) in `reply11.txt`) is a list of distinct rows `r_1..r_M` in one unit class `eta` and Gaussian
integers `z_1..z_M` with `z_j - z_i = 2 i eta H_ij t_ij`, `t_ij` a nonzero rational integer,
`|t_ij| <= (C/2) sqrt(Q_ij) q^{-1/4}` (equivalently `|z_j - z_i| <= C q^{1/4}`); the norms `|z_i|^2` are free. A
*fake* is such a configuration that is not a translate of a genuine cluster.

## 0. Verdict in five lines

1. **No mathematical error in Part XI.** Both referees passed every proof step. I rechecked the dilated collision
   bound (XI.5) and the rigidity (XI.6) line by line and found no gap [J] (§1.6); I reproduced the `q_200` row and
   the cocycle referee's triple from scratch [J].
2. **(A1): a constant, not a rate, and not new.** `M_max(q_P, C) <= (1/4 + o(1)) s` for `P >= P_0(C)`. For genuine
   circles this is the record's rough-parent bound (`growth/auxiliary.md` l.74, `ref8_extend.md` Prop 3.1) plus
   `M < P`. So README_X's "proved bound of order `s/3`" was already out of date. Truth (computed): `M = 1`.
3. **(A2): no fake with `M >= (1/4 + eps) s` exists (Proved); small fakes do exist** (`M = 3` on `q_200` at `C >= 3.574`,
   where genuine clusters have `M = 1`; `M = 4` on `q_30` and `q_50` at `C = 8` [J]). Unbounded fakes are undecided
   (heuristically none).
4. **Rigidity (XI.6) settles what (T1) is.** With the genuine `H_ij`, every configuration of `M >= 3` points is
   `w + (1/v) S`, where `S` is a set of genuine norm-`q` points lying on an arc `v` times the critical length,
   pairwise congruent modulo `v`. The genuine circle is encoded in the `H_ij`, so dropping `|z_i|^2 = q` dropped only
   a scale. (T1) is the genuine problem over the dilations `v >= 1`. It is not a weaker input a proof could use.
5. **For the integer problem nothing moves.** The rate `(2/3 + eps) log R/loglog R`, the uniform bound and Questions A
   and B are unchanged. The thread closes here. Score on the goal: **1.5/10**.

## 1. Referee-corrected statements

Labels: **Proved**: a referee checked a complete proof line by line, or an exact computation was reproduced with
independent code. **Computed**, **Heuristic**, **Gap**, **False**, **Withdrawn** as usual. No GAP or FALSE verdict
was returned on any claim of Astra's. The one CORRECT_WITH_CORRECTIONS from each referee concerns quantifiers, scope
and wording.

### 1.1 The dilated collision bound (XI.3)-(XI.5) (rough: CORRECT; cocycle: CORRECT)

**Proved.** *Let `q` be as above, `m >= 2`, and `A_1..A_m` distinct canonical marks (units `eta_i` may vary). Let
`lambda in Q^x` and suppose that for all `i < j`, `B_ij = lambda (A_j - A_i)/H_ij in Z[i]` and
`|lambda (A_j - A_i)| <= C q^{1/4}`. If no prime `<= m` divides `q`, then*

    W >= 4(m-1) log m - O(m) - 4 m log^+ C,

*with an absolute `O`, uniformly in `lambda`.* At `lambda = 1` this is verbatim the record (`growth/auxiliary.md`
l.74; `ref8_extend.md` Prop 3.1, last paragraph). The new content is the denominator `v` of `lambda`. It is handled
by Bezout (`v | D_ij`), by the class count modulo `l^{a+1}`, `a = ord_l v` (exactly `l` norm-`q` lifts of the
common class modulo `l^a`), and by `lv | D_ij => l | B_ij`.

### 1.2 The bound on the explicit parents (XI.1) (rough: CORRECT WITH CORRECTIONS, quantifier)

**Proved.** *For every `C > 0` there is `P_0(C)` such that for all `P >= P_0(C)`, on both families,
`M_max(q_P, C) <= (s/4)(1 + O_C(loglog P/log P))`. Hence for every fixed `c < 1/12` (in particular `c = 1/24`),
`M_max(q_P, C) <= (1/3 - c) s` for `P >= P_0(C, c)`.* The threshold is needed: for tiny `P`, `s` can be 0 while `M`
can be 4. **Novelty for genuine circles: none.** Theorem 1.1 gives `M <= (1/3 + o(1)) s < P`, so no prime `<= M`
divides `q_P` and the record's coefficient 4 applies. Astra's "already implicit in the record" understates this: the
bound is explicit there and needs only `M < P`. Astra's preliminary step (no arc has `floor(P/2)` marks) is not
needed for genuine circles (the rough referee is right). It **is** needed for fakes, where Theorem 1.1 is not
available, and it is valid because subsets of a configuration are configurations with the same `lambda` (the cocycle
referee is right). [J]

### 1.3 Canonical identity (cocycle: CORRECT)

**Proved.** *In one unit class, `A_j - A_i = 2 i eta H_ij T_ij` with `T_ij` a nonzero rational integer,
`0 < T_ij^2 < Q_ij`, `|T_ij| = gcd(Re(A_j - A_i), Im(A_j - A_i))/2`, and `p` does not divide `T_ij` for every `p`
separating the pair. If `q > C^4/16`, every pair of an integral-cofactor configuration has distinct rows.*

### 1.4 Rigidity (XI.6) and the characterisation (cocycle: CORRECT WITH CORRECTIONS, quantifiers and a sharpening)

**Proved.** *Let `M >= 3` and `q > C^4/16`. If `z_1..z_M` is an integral-cofactor configuration, there are
`lambda = u/v in Q^x` (lowest terms, `v > 0`) and `w in v^{-1} Z[i]` with `z_i = w + lambda A_i` for every `i`.* The
proof uses only that each `z_j - z_i` is a real multiple of `A_j - A_i` (exact chord direction). Rationality of
`lambda` comes from integrality. For `M = 2` the statement is empty: every pair with `Q_ij >= 4 sqrt(q)/C^2` is a
2-configuration (`lambda = 1/|T_ij|`).

**Proved (referee's sharpening).** *Let `S` be a set of `M >= 3` canonical marks in one unit class and
`G(S) = gcd_{i<j in S} T_ij`. `S` carries an integral-cofactor configuration at constant `C` if and only if
`diam(S) <= C q^{1/4} G(S)`.* The admissible `lambda` are the `u/v` with `v | G(S)` and `|u| diam(S) <= v C q^{1/4}`,
and `lambda = 1/G(S)` is optimal. Plainly: a fake is a translate of `1/v` times a set of genuine lattice points of
`|z|^2 = q` in one unit class, pairwise `T`-congruent modulo `v` (for odd `v` this is ordinary congruence
`A_i = A_j mod v Z[i]`), lying on an arc `v` times longer than critical. It is genuine exactly when `v = 1`. Exact
exponent rows can always be imposed afterwards by a CRT translation, so the row condition in (A2) adds nothing.

*Scope.* (T1) is defined inside one unit class. Under the extension "`z_j - z_i in Q (A_j - A_i)` and
`H_ij | z_j - z_i`" across classes, the same proof and the same bound hold.

### 1.5 The s/4 bound for configurations, and residue pinning (cocycle: CORRECT)

**Proved.** *For `P >= P_0(C)`, on both families, every integral-cofactor configuration on `q_P` has at most
`(s/4)(1 + O_C(loglog P/log P))` points, uniformly in `lambda` and `w`. In particular no fake with
`M >= (1/4 + eps) s` exists for `P >= P_0(C, eps)`.*

**Proved (pinning).** *Let `M >= 2`, `C^2 < 4 sqrt q`, and let `S` be a set of primes of `q` with
`prod_{p in S} p >= min_{i<j} Q_ij`. If for every `p in S` some pair separated at `p` keeps its canonical squared
chord, `|z_j - z_i|^2 = |A_j - A_i|^2 mod p`, then `lambda = +-1` and `z_i = w +- A_i` with `w in Z[i]`.* Under the
weaker "local vertex value" reading of the signed residues (free values on the two isotropic lines whose products
give the squared chords), dilated fakes keep the law, because every squared chord is multiplied by the same `v^{-2}`.

### 1.6 Judge's line-by-line checks [J]

*(XI.5) with the denominator.* (1) `D_ij = eta_j X_j - eta_i X_i in Z[i]`. (2) `alpha u + beta v = 1` gives
`D_ij/v = alpha B_ij + beta D_ij in Z[i]`. No coprimality is used. (3) For an odd `l <= m`, `l` does not divide
`q`, `a = ord_l v`: since `v | A_j - A_i`, all marks agree modulo `l^a`; any representative `(x, y)` of that class
satisfies `x^2 + y^2 = q mod l^a`, because the marks have norm exactly `q`. The lift condition is linear modulo `l`
because `2a >= a + 1`. It is nondegenerate because `(x, y) != (0, 0) mod l`. So there are exactly `l` classes for
`a >= 1`, and `l - (-1/l)` classes for `a = 0` (brute force for `l <= 13`, `a <= 2`, `j1_triple.out`). (4) A
collision modulo `l^{a+1}` gives `l^{a+1} | H_ij D_ij`. No Gaussian prime above `l` divides `H_ij`, so
`l^{a+1} | D_ij`. With `v = l^a v'` and `v' | D_ij` this gives `lv | D_ij`, hence `B_ij = l u (D_ij/(lv))`. (5)
Distinct collision primes are coprime rational integers dividing the nonzero `B_ij`, so `log|B_ij| >= sum log l`.
(6) The colliding fraction at `l` is at least `(m/(m-1))/(l+1) - 1/(m-1)`. Summing these lower bounds, including
negative ones, is legitimate. Mertens gives `log m - O(1)`. (7) The upper bound is
`log|B_ij| <= log C + W/4 - W/2 + (1/2) log Q_ij`, and the cut count gives a mean of at most
`log C + W/(4(m-1))`. **No gap.** The only place the exact norm enters is step (3): it is the norm **of the canonical
marks**. This is why dropping `|z_i|^2 = q` costs nothing (§3(a)).

*(XI.6).* The fake differences are `lambda_ij (A_j - A_i)` with `lambda_ij = t_ij/T_ij in Q^x`. `T_ij != 0` because
all rows are distinct (1.3). In a triangle,
`(lambda_ij - lambda_ik)(A_j - A_i) + (lambda_jk - lambda_ik)(A_k - A_j) = 0` with real coefficients. Three distinct
points of a circle are not collinear, so the two vectors are independent over `R` and the three ratios are equal.
Triangles through one fixed mark give a single `lambda`. **No gap.** The pinning argument was also checked: `p` does
not divide `v`, because `v | T_ij`, and `p` does not divide `T_ij`; `u | t_ij` holds for every pair.

*Computations.* `j1_triple.out`: on `q_200 = 4810942114471671677256758751310693840381` exactly one lattice point `A`
has `A + G z_2` and `A + G z_3` on the circle (`G = 4391525450`, `z` as in 1.7). All three points have unit 1. Their
rows differ at 9, 13 and 8 primes. `T = -105396610800, -7724047712308850, -4391525450 = G * (-24, -1758853, -1)`.
The minimal `C` is `3.5732`. The configuration passes at `C = 4` and fails at `C = 3`. `gcd(G^2 - 1, q) = 1`, so no
canonical residue is kept. `j2_dmin.out`: on `q_200`, `d_min^2 = 15219819238844501369583650` and the floor is 468
(family 1); on family 2 the floor is 25. Both agree with the referees. `j3_small.out` and `j4_gen.out` (exhaustive,
one unit class, `C = 8`):

| `P` | `s` | marks per class | 2-configurations | 3-sets | 4-sets | max fake `M` | genuine max (all units) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 30 | 3 | 8 | 28 of 28 | 10 | 2 | 4 | 4 |
| 50 | 5 | 32 | 416 of 496 | 28 | 1 | 4 | 3 |
| 100 | 10 | 1024 | 326,656 of 523,776 | 424 | 0 | 3 | 2 |

So the integral-cofactor cocycle admits slightly more than the circle already on small parents, and its yield falls
with `P`, as the heuristic 1.8 predicts.

### 1.7 Computed (reproduced)

* **Exact rows (Proved certificate; rough reproduced both rows, judge the first).** `P = 200`: `s = 16`, 65,536 unit
  orbits, `d_min^2 = 15219819238844501369583650`, `floor(d_min/q^{1/4}) = 468`. `P = 250`: `s = 20`, 1,048,576
  orbits, `d_min^2 = 675416293412257838261903282453746`, floor 4556. Family 2 (rough): floors 25, 20, 174 at
  `P = 200, 250, 300`. Each row certifies `M_max(q_P, 8) = 1`.
* **A fake triple on `q_200` (cocycle; reproduced [J]).** Rows 436, 12399, 50763 (`eta = 1`), `G = 4391525450 =
  2 * 5^2 * 13 * 101 * 151 * 443`, `lambda = 1/G`, `z = 0, (3439204080, 29559507552), (11203893610, 21049952704)`,
  `t = -24, -1758853, -1`. It is valid for every `C >= 3.574` and invalid at `C = 3`. Its exponent rows become exact
  after a translation.
* **Search coverage (cocycle; partial).** The smooth-modulus search over all 65,536 bases at `C = 8` found 4 copies of
  this shape and no 4-set. `10^5` random pairs over all moduli gave 0 hits. 59.9% of pairs are 2-configurations at
  `C = 4` or `8`. This does not prove that no 4-set exists on `q_200`.

### 1.8 Heuristic

The expected number of `M`-point configurations is about `E_M = C^{M-2} 2^{sM}/M! (2 pi q^{1/4})^{-(M-2)}`. It is
dominated by `v ~ q^{1/4}/C`. The referee corrected the congruence cost: it is about `1/v` per extra point, not
`v^{-2}`. The heuristic predicts `M = 3` up to `P ~ 2500` and `M = 2` beyond. It overestimates the observed triples
by one to two orders of magnitude. **Not a proof that fakes are bounded.**

### 1.9 Tests (rough: CORRECT)

**Proved.** In the integer twin the same count gives coefficient 2, i.e. `m <= (1/2 + o(1)) log q/loglog q`,
consistent with (T5). For VIII.3's angle fakes, `23/3 > 4` is consistent. On squarefree rough split parents in a
constant quadratic extension of `F_b[t]`, `deg q >= 4(m-1) log_b m - O_b(m)` (undilated). Not re-checked: the degree
orders quoted for the Part V-VI families.

### 1.10 Wording corrected

* "already implicit in the record" (XI, opening; `final11.txt`): **explicit**, up to the one-line specialisation
  `M < P`.
* "(XI.1)" without a threshold: add `P >= P_0(C)`.
* The summary "the integer-cofactor cocycle is not weaker than the norm equation up to dilation and translation"
  (put to the cocycle referee) is **correct structurally for `M >= 3` and false as a count**. The configurations are
  exactly the `w + (1/v) S`, a union over `v` of problems on longer arcs with congruence modulo `v`, and the union is
  strictly larger (3 against 1 on `q_200`).
* README_X §0 line 5 and §3(e) ("a known bound of order `s/3`", "proved bound of order `s/3`"), repeated in
  `reply11.txt` and in the expert note: **out of date**; the record gave `(1/4 + o(1)) s`.
* `final11.txt` links to `/scratchpad/astra/free/answer.md`; the refereed copy is `round31/astra_free_XI/answer.md`.

Nothing is withdrawn.

## 2. What is new relative to the record, and correct

* **Not new: the coefficient 4 on rough parents, hence `s/4` on `q_P` for genuine circles.** The record states
  `W >= 4M log M - O_C(M)` whenever no prime `<= M` divides `N` (`growth/auxiliary.md` l.74, [Sieve]), and again as
  `ref8_extend.md` Prop 3.1 ("If no split prime `<= M` divides `N'`, this gives `W' >= 4M log M - O_C(M)`"). With
  `M < P` (Theorem 1.1) this gives (XI.1) for circles at once. The record had simply not specialised it to `q_P`.
* **New, small (Proved): the collision bound is uniform under a rational dilation** (1.1). This is the class count
  modulo `l^{a+1}`, about five lines beyond [Sieve].
* **New (Proved): rigidity (XI.6)** and the referee's **iff criterion** `diam(S) <= C q^{1/4} G(S)` (1.4). Elementary,
  and it changes the map (§3(a)).
* **New (Proved): pinning by canonical squared-chord residues** (1.5). Elementary.
* **New (Computed): the first integral-cofactor fake beyond pairs** (`q_200`, `M = 3`), and the small-parent table
  [J]. **New (Heuristic):** the size count 1.8.
* **New data (Computed):** family-2 rows at `P = 200, 250, 300` (rough).

## 3. Rulings

**(a) What of the final question was answered.**

* **(A1): a constant, not a rate.** `M_max(q_P, C) <= (1/4 + o(1)) s`, so every `c < 1/12` works. The answer comes
  from the record, not from Part XI. `M = o(s)` on `q_P` was not obtained. No Diophantine inequality on the
  `theta_p` was isolated, and no fake on the same prime set was built to show that (T1) fails to imply one.
* **(A2): answered in its "`M >> s`" form, negatively; open in its "unbounded" form.** No integral-cofactor
  configuration has `M >= (1/4 + eps) s` (Proved). Configurations with `M = 3` exist on `q_200`, and with `M = 4` on
  `q_30` and `q_50` (Computed). Whether `M` is bounded on `q_P` is undecided (heuristically yes: 3, then 2). The
  question asked what the proof uses, and whether it is more than "such `z_i` would be an actual cluster". It uses
  exact chord directions (rigidity), integrality (rational `lambda`, plus the collision count with denominator), and
  the exact norm **of the canonical marks** that the `H_ij` determine. It is slightly more than "an actual cluster":
  the conclusion is "an affine image of genuine points on a `v`-fold longer arc, congruent modulo `v`". If the
  signed residues are read as canonical numerical values (as the question asked: "kept"), it is **exactly** that
  statement. Then `lambda = +-1`, every configuration is a translated genuine cluster, and on `q_200` and `q_250`
  none with `M >= 2` exists at `C <= 8`.
* **The map of inputs: which input carries the separation.** Bluntly: **the genuine Gaussian primes do.** (T1) is
  posed with the actual `H_ij`, the products of actual Gaussian prime factors. Once those are fixed, the exact
  direction of every chord pins the configuration to `w + lambda * (genuine marks)` (XI.6). Integrality makes
  `lambda` rational and keeps the collision count. The only thing left to the norm equation is the scale
  `lambda = +-1`, worth at most a union over `v` and at most two extra points in every case computed. Each
  ingredient fails alone:
  * real cofactors with exact directions: rigidity holds with `lambda` real, and `lambda = C q^{1/4}/(2 sqrt q)` puts
    all `2^s` marks into one configuration;
  * integral cofactors without exact directions are valuation/twin-type data, which stop at the larger-sieve
    constant `1/2` (IX of the other thread, (T5));
  * angle models have no exact directions at all.

  So the reading of `reply10.txt` (T3), "the input that separates is (T1)", is confirmed but reduced. (T1) separates
  because it contains the Gaussian factorisation of the actual points, i.e. the circle itself. It is not an
  intermediate hypothesis weaker than the norm equation. **A rate from (T1) on `q_P` is at least as hard as a rate for
  genuine circles on `q_P`** (the case `v = 1` is included). The expert note's entry "the Gaussian-integer cocycle of
  the reduced chords" in the uncapped list moves to "the problem itself, up to a rational dilation".

**(b) Consequences for the integer problem.** **None for the rate or the uniform bound.** `M <= (2/3 + eps) log R/loglog R`
for general `R`, `M(C) < infinity` and Questions A and B are unchanged. On the explicit dense parents `q_P` the
record's bound is `s/4` (not `s/3`), and the truth is `M = 1` for `200 <= P <= 450`. The gap between `s/4` and 1 is
a question about signed subset sums of about `s/2` arguments `2 theta_p` near multiples of `pi/2` at the Liouville
threshold. That is linear forms in logarithms with unboundedly many terms, beyond Baker-type methods. Part XI
shows that no relaxation through the integer cocycle sidesteps it. The cocycle thread, the last internal lead named
by the coordinator and the judges, is closed: it is the problem itself.

## 4. Recommended edits to the expert note (recommendations only; nothing was written to the repository)

1. **§1, after "Actual clusters satisfy $W\ge3M\log M-O(M)$."** insert: "If no prime $\le M$ divides $N'$, then
   $W'\ge4M\log M-O_C(M)$ (`growth/auxiliary.md` l.74; `round31/ref8_extend.md` Prop 3.1). In particular, on
   $q_P=\prod_{p\equiv1\,(4),\,P\le p\le2P}p$ (or the primes $\equiv5\ (8)$ there), $M_{\max}(q_P,C)\le(\frac14+o(1))s$
   for $P\ge P_0(C)$, $s=\omega(q_P)$, and the bound holds uniformly for the integral-cofactor configurations below
   (round 31, Part XI.1; **proved**)."
2. **§1, "Not everything is capped" paragraph: replace** ", Lang–Waldschmidt, and the Gaussian-integer cocycle of the
   reduced chords." **by:** ", and Lang–Waldschmidt. The Gaussian-integer cocycle of the reduced chords is not a
   separate input: with the genuine common factors $H_{ij}$, any $M\ge3$ Gaussian integers with
   $z_j-z_i=2i\eta H_{ij}t_{ij}$, $t_{ij}\in\mathbb Z\setminus0$, satisfy $z_i=w+\lambda A_i$ with $\lambda\in\mathbb Q^\times$
   and $A_i$ the genuine points, and a set $S$ of genuine points of one unit class carries such a configuration iff
   $\operatorname{diam}S\le Cq^{1/4}\gcd_{i<j}T_{ij}$ (round 31, Part XI.2; **proved**, `round31/README_XI.md` §1.4).
   So the cocycle is the genuine problem on arcs $v$ times longer, with a congruence modulo $v$. It admits slightly
   more than the circle ($M=3$ on $q_{200}$ at $C\ge3.574$, where genuine clusters have $M=1$; **exact**), and it is
   exactly the genuine problem when canonical squared-chord residues are kept."
3. **§5 item 12, header:** "ten exchanges" -> "eleven exchanges"; "108 referee blocks and seventeen judge reports" ->
   "110 referee blocks and eighteen judge reports"; add `README_XI.md` to the list and `1.5/10` to the scores.
4. **§5 item 12, Part X bullet: replace** "against a proved bound of order $s/3$ points" **by** "against the record's
   bound $(\frac14+o(1))s$ (§1; Part XI.1)".
5. **§5 item 12, add after the Part X bullet:**
   "* ***Part XI: the cocycle on one fixed parent (final)*** (`round31/README_XI.md`; two referee blocks; nothing
   false; judge's score on the goal 1.5/10). (A1): $M_{\max}(q_P,C)\le(\frac14+o(1))s$ for $P\ge P_0(C)$ on both
   dyadic families. For circles this is the record's rough-parent coefficient 4 together with $M<P$. New is only its
   uniformity under a rational dilation. (A2): rigidity, i.e. integral-cofactor configurations with the genuine
   $H_{ij}$ are $w+\frac1vS$ with $S$ genuine, $v$-congruent, on a $v$-fold arc. So none has
   $M\ge(\frac14+\varepsilon)s$, but small ones exist ($M=3$ on $q_{200}$, $M=4$ on $q_{30},q_{50}$; **exact**), and
   whether they are bounded is **open** (heuristically $\le3$). Canonical squared-chord residues at primes of product
   $\ge\min Q_{ij}$ force $v=1$. *Judge's reading*: (T1) separates fakes from circles only because it contains the
   genuine Gaussian factors. It is the problem itself up to a dilation, and a rate from it on $q_P$ is no easier than
   a rate for circles. Do not write "Part XI improves the rate", "the $s/4$ bound on $q_P$ is new", "integer
   cofactors force a genuine cluster" (without "with canonical residues" or "up to dilation"), "no integral-cofactor
   fake exists", "integral-cofactor fakes are bounded", or "$M_{\max}(q_P,C)\le s/4$ for all $P$"."
6. **Pointer table:** add "round 31 Part XI: cocycle rigidity on $q_P$ | `round31/README_XI.md`, `round31/ref11_*.md`";
   add `README_XI.md` to the round-31 row.

## 5. Closing assessment of the whole exchange (eleven parts on the main thread, ten on the Question B thread)

(i) **Theorems about integer circles gained.** *Theorem 1.2* (Part VIII.1): the affine rank bound
`M <= 1 + omega(N') + K_C`, hence `P_+(N') >= (2 - o(1)) M log M` and the coefficient `2/kappa` (`1/2` for squarefull
`R`). It is the one result that changes the stated rate for a natural class of radii. *Theorem 1.3* (Part IX.2): a
projection/larger-sieve count, position-dependent, with nothing for generic arcs. *The random-parent rate* (Part III,
extended in IX of the other thread): `O(log s)` at the smooth scale and `O(1)` at `P = s^{1+eps}`, outside an
exceptional fraction, over a density-zero family. *The `s/4` bound on rough parents*: **not new** (record), now known
to hold uniformly for dilated integral-cofactor configurations. None of these moves the worst-case `2/3`.

(ii) **Barrier results.** Window fakes (VIII.3, extended to `delta < 1/6` in IX.4); P3-char capped at every modulus
(IX.3); cross-gap fakes (IX.5); a rational-direction window fake (IX.6); a partial P3-all\* fake, with the height
route closed (X.1-X.3); the same-norm comparison (fakes `~ s/2` against genuine `O(1)` on the same `N`); the twin
relaxation at the sieve constant `1/2`; the sign-free and integral-cofactor relaxations collapse to the genuine
problem (X of the other thread, XI). In the function field: a uniform bound in characteristic 0 (`8 <= sup M(D,0) <= 24`),
unbounded clusters over `F_2(t)` (also squarefree), over `F_p(t)` for `p = 3 mod 4` and at cap two, and filter
F3+++. On Question A: `(L_8)` implies Rudin's conjecture, and the critical twin is Cilleruelo-Granville's
Conjecture 9. Net effect: every relaxation studied either admits fakes at `W = O(M log M)` or collapses to the
problem itself. Nothing lies strictly between.

(iii) **Corrections to the record.** The printed `(L_8^{pc})` is vacuous; round 8's QS at an absolute threshold is
false; the cocycle was mislabelled as capped (narrowed, and now moved to "the problem itself"); the Rudin
attribution was corrected; the "17 points" headline was superseded; the bound on `q_P` was `s/4`, not `s/3`; the
coordinator's "twist must be `+-1`", "weight `D/2` is the lattice minimum" and F3+ exemptions were corrected; the
fallback of Question A was corrected.

(iv) **Open, and where to point an outside expert.** Open: the growth rate below `(2/3) log R/loglog R` for general
`R`, the uniform bound, and Questions A and B. The judges' estimates for the remaining internal items (point
residues to `3 + eta`, about 4%; interpolation for `2 < kappa < 3`, about 2%) are not worth more rounds. If one
direction is offered, it is **the explicit-parent problem**: prove `M_max(q_P, C) = o(omega(q_P))` for the
product of all split primes in `[P, 2P]`. It is concrete, the truth is computed (`M = 1`), and it is extremal for
Theorem 1.1. Every relaxation is now known either to admit fakes or to be the problem itself. It is a question about
small signed subset sums of arguments of Gaussian primes at the Liouville threshold, for a specialist in linear forms
in logarithms or in additive combinatorics of Gaussian prime angles. Honest probability that a serious expert effort
yields `o(s)` on `q_P`: about 3%. That it moves the worst-case rate for general `R`: at most 1%. The uniform bound:
below 0.5%. The alternative (Question B, i.e. `(L_5^o)`, for function-field and Vojta specialists) is stated in §3 of
the note and needs no further internal work.

(v) **Further internal rounds of this kind: not worth running.** Twenty-one parts produced two theorems on circles,
one of them restricted, a strong map of barriers and many corrections. Since Part VIII the score has stayed between
1.5 and 3. The last three parts on each thread ended in "this relaxation is the problem itself". The marginal
round now returns a structural statement or a constant that the record already had. Close both threads, apply §4,
and hand the note to outside experts with the single pointer in (iv).

## 6. Scores (0-10, on the integer goal)

| block | correctness | score on goal |
|---|---|---|
| XI.1 dilated collision bound (XI.3)-(XI.5) | correct; new only through `lambda` (five lines) | 1 |
| (XI.1) on `q_P` | correct with a threshold; for circles already in the record | 0.5 |
| XI.2 rigidity (XI.6), characterisation, pinning | correct; elementary; settles what (T1) is | 2 |
| XI.3 exact rows, certificate; twin and function-field tests | reproduced; consistent | 0.5 |
| referee additions (`q_200` fake triple, iff criterion, heuristic) | computed/proved; the triple reproduced [J] | 1 |
| **Part XI overall** | **no error; honest scope; novelty overstated only by understating the record** | **1.5** |

As mathematics the part is short, clean and correct (about 7/10). It answers the final question as well as it can be
answered internally: the cocycle is not a relaxation, so the thread had nowhere further to go.
