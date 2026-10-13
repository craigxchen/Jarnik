# Referee report on round9/L1L2.md (adversarial)

Referee: adversarial pass, 2026-10-13. Scripts and outputs are in `referee_L1L2_checks/`. All of them use my
own Gaussian-integer code (`refgauss.py`); none of the author's code is imported. The author's `s9.pkl` was
read only as data, after its pickle opcodes were checked to contain no `GLOBAL`/`REDUCE`, so that the two
cluster sets could be compared.

## 0. Verdict

**The mathematics survives. The headline does not.** I re-derived every proof in §§1-3 and §5 by hand and
found no error in any displayed inequality or identity. I also tested the exact statements (Thm 1.1, the
integer form of (1.1), the gcd formula, the bisector step of Lemma 1.3, the congruence statement, the
discrete-log code of Prop 3.1, the area identity, and Lemma 5.3) with my own exact code: on every pair of
all 22,907 scan clusters, and on four record clusters, including the 10-point record that the write-up
omits. Nothing failed. I re-ran the scan with an independent method and got the author's cluster set
exactly.

These things are wrong or overstated:

1. **"On near-extremal actual clusters L1 and L2 are each other's negation"** (summary, §0, §5.1 item 4) is
   **false as logic.** Thm 3.2 proves only that the two properties are *mutually exclusive* on any single
   near-extremal cluster with `M >= M_1`. Non-uniformity does not imply a heavy pair. If growth holds at
   constant `c_0`, both lemmas are **true** (vacuously), so they are not each other's negation. The same
   holds for (R1-arch)/(R1-arith): (R1-arith) is far stronger than the negation of (R1-arch). The
   conclusion drawn in §5.1 item 4, that L1 and L2 "cannot be pursued as independent lemmas", is a non
   sequitur. In a proof by contradiction the two lemmas are always incompatible on the hypothetical object.
   (Major.)
2. **Prop 2.2 and open point 6.3 are superseded by a refereed theorem the write-up itself cites.**
   `p3lower.md` Thm A (valid on `A_X`, hence for actual clusters) gives
   `W >= 2 lambda_max(D) - 2(M-1) c_C^+`. So every near-extremal cluster has
   `lambda_max(D) <= (c_0/2) M log M + (M-1) c_C^+`. The "window
   `omega(M log M) <= lambda_max(D) <= O(M^(3/2) log M)`" in which the alternative form of L1 is said to
   "live" is **empty**, and open point 6.3 ("is `lambda_max(D) <= O(M log M)`?") is already answered "yes".
   Prop 2.2 survives only as a bound on the point totals `max_x D_x`. (Major.)
3. **Two evidence readings are false** (§4 and the summary's evidence bullet). (Major for the evidence
   section; no proof depends on it.)
   * "Without complementary pairs the median `K` drops to 1.17 / 3.28" is an artefact. The statistic
     removes complementary pairs from the maximum but keeps their margin `W/2` in the mean. Natural
     versions give median `K` = 2.67 (`M = 4`) and 3.54 (`M = 5`) on clusters that have **no**
     complementary pair. Removing those pairs from both the maximum and the mean gives 5.00 and 4.17,
     with 8,372 four-point clusters where the remaining mean is `<= 0`.
   * "Every margin grows linearly in `W`" in the golden family is false. In the `n = 7` member, 24 of the
     28 margins are **identical** to those of the `n = 4` member, and only 4 grow. `K -> 7 = 28/4` is a
     counting artefact of extreme non-uniformity.

   So the data do not support the reading that small-cluster non-uniformity is "driven by symmetry, a
   mechanism Lemma 1.3 confines to bounded `M`".
4. Smaller overstatements:
   * **Thm 1.1:** novelty is modest. It is the pair case, with equality, of `sharp.md` Prop 1.2 and the
     `p3lower.md` Def 1.1 bullets. Its `B_xy` also clashes with `p3lower.md`'s `B_xy`, which is the
     Gaussian integer `zhat_x/g_xy`.
   * **Thm 3.3:** "without `H`" is misleading. It assumes a heavy pair outright, which is what `H` was
     used to produce.
   * **"Lemma 1.3 is the only known upper bound on an individual margin":** the triple identity already
     gives another one.
   * **"Integrality gives only lower bounds":** Lemma 1.3 is itself an integrality argument that gives
     an upper bound.
   * **Prop 3.5(d):** the bold heuristic conclusion contradicts its own caveat.
   * **Cor 3.4.1:** it drops the `M >= M_0(C)` of Thm 2.1(b). At `M = 4` there are near-extremal
     clusters with mean margin 0, on which L1 fails for every `K`.

**What is proved, with all quantifiers, for actual clusters** (all re-derived): Thm 1.1, Cor 1.2, Lemma 1.3,
Thm 2.1, Prop 2.2 (as a point-total bound), Prop 2.3(a)(b), Prop 3.1, Thm 3.2, Thm 3.3, Cor 3.4 items 1
and 3 (item 2 with "together" read as conjunction, not negation), Prop 3.5(a)-(c), Props 5.1-5.3. **No
growth, no new barrier.** Thm 3.2 and Cor 3.4 are the contrapositive of Thm C's two-line proof, written in
the variables `B_xy`, `t_xy`. They locate the difficulty but do not obstruct any argument. The status
"partial" is correct.

## 1. Theorem 1.1 (reduced chord) and Lemma 1.3 (bisector)

**Normalisation.** Correct. Inert primes `q^(2f) || N` force `q^f | z` for every point. `(1+i)^(v_2(N))`
divides every point. Both are absorbed by `G`, so `N'` is a product of split primes, and after division
`min_x a_xj = 0`, `max_x a_xj = e_j`. Also `C' = C|G|^(-1/2)`, because the arc length scales by
`1/|G|` and `sqrt R` by `|G|^(-1/2)`.

**Thm 1.1: correct, proved, actual clusters only.** I re-derived the four cases:
* `A - conj A = 2iY`;
* `A - i conj A = (X-Y)(1-i)`;
* `A + conj A = 2X`;
* `A + i conj A = (X+Y)(1+i)`.

I also checked the edge cases.
* `v = 0`: then `z_y = i^s z_x` and `B = +-1`.
* `v != 0`: `A_v` is neither real nor imaginary nor on a diagonal, because of unique factorisation with
  `pi_j`, `conj pi_j` non-associate. So `B != 0` follows independently of `z_x != z_y`.

The congruence part needs `q` not dividing `N`. Then `q` does not divide `G` or `g_xy`, since their norms
divide `N`.

*Exact tests* (`records.py`, `stats_ref.py`). On all 137,862 pairs of the 22,907 scan clusters, and on the
9-, 10-, 7- and 8-point records, the following hold exactly:
* `beta = unit * gamma_s * D_s(A_v)`, with the `k_x` and `a_xj` computed by my own Gaussian valuations;
* the gcd formula `Norm g = prod p^(e - |a_x - a_y|)`, with `g` from the Euclidean algorithm;
* `|u_x - u_y|^2 = Norm g * nu_s * B^2`.

On the records I also checked the following, with no failure:
* `u_x = u_y mod q^a <=> q^a | B <=> z_x = z_y mod q^a`, for all odd `q <= 200` not dividing `N` and
  `a <= 3`. The last form uses the un-normalised points.
* Prop 3.1 (`w_x(q) = w_y(q) <=> q | B`) for `q <= 60`, with explicit discrete logarithms in the norm-one
  torus. The torus orders `n_q = q - (-1/q)` were confirmed.

*Novelty.* The monomial statement "`q^a | D_s(A_v)`, a nonzero rational integer with
`|D_s(A_v)| = nu_s^(-1)|A_v||sin(<v,phi> - s pi/4)|`" is `sharp.md` Prop 1.2, re-read in `p3lower.md` §1.1.
The reduced chord and its collision meaning are in `p3lower.md` Def 1.1 and `round7/README.md` §3 (L2
bullet). Theorem 1.1 adds the **equality** `beta_xy = i^(k_x) gamma_s D_s(A_v)` for `v = a_x - a_y`. This is
a useful normalisation, not a new theorem.

*Notation clash.* In `p3lower.md` Def 1.1, `B_xy` denotes the Gaussian integer `zhat_x/g_xy`, with
`|B_xy| = e^(d_xy/2)`. Here `B_xy` is the rational integer `D_s(A_v)`. Rename one of them before merging.

**Cor 1.2: correct.** (1.2) is `p3lower.md` Lemma 1.3, with the cofactor `|B|/m >= 1` kept.

**Lemma 1.3 (bisector): correct and, as far as I can find, new.** No `bisector` argument occurs anywhere in
`docs/`. My re-derivation:
* `u_x + u_y = g i^(k_x)(A + i^s conj A)`, which is `g i^(k_x)` times `2X`, `(X+Y)(1+i)`, `2iY` or
  `(X-Y)(1-i)`. So the bisector direction is that of `g`, `ig` or `g(1+-i)`. This dual of Theorem 1.1 is
  confirmed exactly on all 137,862 pairs.
* The arc contains the direction `arg v` (`eta = 0`). Every point therefore has
  `Re(z conj v)/|v| in [R cos(L/R), R]`. That is at most `|v| L^2/(2R) + 1` integer values, with at most
  2 points per line.

So `M <= |v| C^2 + 2`, whatever the position of the arc. This is `outside.md` Prop C.2, case `eta = 0`;
primitivity of `v` is not needed. The lemma had no violation; its smallest slack
`sqrt2 C^2 |g| + 2 - M` over the scan is `4.16`, so it is never close to tight there.

Two qualifications to the surrounding prose:
1. *"The only pair-level upper bound on margins that I know"* (Remark 1; the summary drops "that I know").
   The triple identity (Prop 5.1(b)) gives `mu_xy + mu_yw + mu_xw = W/2 - 2 log Norm g_xyw <= W/2`. With
   `mu >= log 2 - 2 log C` for the other two pairs, `mu_xy <= W/2 - 2 log 2 + 4 log C`. This is another
   individual upper bound, weaker than Lemma 1.3 for `M >= 5`. Lemma 1.3 is the only one with an
   `M`-dependence. Restate it that way.
2. *"Integrality and divisibility give only lower bounds on `|B_xy|`"* (§5.1 item 1, and the summary's
   "Integrality gives only lower bounds").
   * Restricted to divisibility of the `B_xy` themselves (cocycle, residues), this is right.
   * Lemma 1.3 is itself an integrality argument: `Re(z conj v) in Z`. Through (1.1) it gives the
     **upper** bound `|B_xy|/t_xy <= sqrt2 C^2 N^(1/4)/((M-2) sqrt(nu_s))`.
   * The summary sentence should be restricted to "integrality of the reduced chords".

## 2. Theorem 2.1 (L1 <=> polynomial size), Props 2.2, 2.3

**Thm 2.1: correct, proved, actual clusters.**
* (a) `mu <= K mubar <= K W/(2(M-1)) <= 2b log M`. The `|B|` and `t` bounds then follow from (1.1),
  `nu_s in {2, 4}`, `t <= C` and `|B| >= 1`. All re-derived.
* (b) is right as stated, with two qualifications:
  * the `O_C(1/log M)` also depends on `a + b`;
  * it holds only for `M >= M_0(C)`, where the lower bound (F) on `mubar` is positive.

*Content.* On near-extremal clusters (F) gives `mubar = Theta(log M)`, so L1 is literally "`max mu = O(log M)`".
By (1.1) that is "`max |B|/t = M^(O(1))`". Theorem 2.1 is this one-line rewrite. It is correct and
clarifying; it is not a structural result.

*Quantifier slip in Cor 3.4.1.* Cor 3.4.1 states the equivalence L1 <=> (P_b) on all near-extremal clusters.
That drops the `M >= M_0(C)` of 2.1(b). `neg_mubar.py` finds 26 four-point clusters in the scan with
`sum_(x<y) mu_xy = 0` exactly (to rounding), some with `A = W/(M log M) = 2.27-3.5`. They are near-extremal
for every `c_0 >= 3.5`, and L1 fails for them with every `K`, while (P_b) holds for large `b`. Harmless
for the asymptotic statement; add "for `M >= M_0`".

**Prop 2.2: the inequalities are correct, but the conclusion drawn from them is false.**
* The first line is `p3lower.md` Lemma 1.4(3) plus Lemma 1.3. That is Prop B.2, which bounds the same
  quantity.
* Cauchy-Schwarz gives `sum lambda |sigma| <= sqrt(W * 4Q)`.
* `stability.md` (4.1) gives `Q <= ((c_0 - 2)/4 + o(1)) M^2 log M`.

I re-derived all three steps. So `max_x D_x <= (1/2) sqrt(c_0(c_0-2)) M^(3/2) log M (1 + o(1))` is proved
for near-extremal actual clusters.

The write-up then transfers this to `lambda_max(D)` and says the alternative form of L1
(`lambda_max(D) >> M log M`) "is not excluded, but lives in the window
`omega(M log M) <= lambda_max(D) <= O(M^(3/2) log M)`". Open point 6.3 asks whether
`lambda_max(D) <= O(M log M)` on near-extremal actual clusters. **Both are already settled by `p3lower.md`
Theorem A**, which is refereed, holds on all of `A_X` (hence for actual clusters), and is cited in the
write-up:

```text
W  >=  W_cut(D - c_C(J - I))  >=  2 lambda_max(D) - 2(M-1) c_C^+ .
```

Re-derivation: take the Perron vector `u` of `D`. Then `u^T (J - I) u <= M - 1`, and part 2 of Theorem A
applies. So on every near-extremal actual cluster

```text
lambda_max(D)  <=  (c_0/2) M log M + (M-1) c_C^+ ,
```

which is better than Prop 2.2 by a factor of about `sqrt(M)`. The "window" is empty for `M >= M_0`.
The alternative form of L1 is excluded on near-extremal clusters in exactly the way Theorem 3.2 excludes
L2 under L1: it is again equivalent to growth at constant `c_0`.

On the records (`records2.py`, exact `D` at all odd `q <= 10^4` not dividing `N'`), Theorem A's ceiling is
already the smaller of the two:

| cluster | `lambda_max(D)` | Thm A ceiling | Prop 2.2 first line |
|---|---:|---:|---:|
| 9-point | 14.2 | 37.9 | 57.4 |
| 10-point | 14.7 | 42.2 | 61.6 |
| 7-point | 10.8 | 25.4 | 37.9 |

**Correction.**
* Prop 2.2 should be stated as a bound on the point totals `max_x D_x` only.
* The sentence "is not excluded, but it lives in the window ..." and summary item 6, second bullet, should
  be replaced by: "by Theorem A, `lambda_max(D) <= (c_0/2) M log M + O(M)` on near-extremal clusters, so the
  alternative form of L1 holds only vacuously, i.e. it is the growth statement."
* Open point 6.3 should be deleted.
* `max_x D_x` can exceed `lambda_max(D)`. Whether `max_x D_x >> M log M` can occur on near-extremal
  clusters is a legitimate, different question. Per `p3lower.md` Cor A.3 it is necessary, not sufficient,
  for the alternative form.

**Prop 2.3: correct.**
* (a) Markov, with strict inequality when `mubar + m_0 > 0`; trivial otherwise.
* (b) Re-derived: `K(mubar + m_0) >= 2 log(|B|/t) + log 2 + m_0`, then `W >= 2(M-1) mubar`.
* (c) A model statement, correctly labelled. I recomputed the "about 0.69":
  `log(binom(M,2) lambda^k/k!)/log M = 0.692` at `log M = 100`, `eps = 0.4`, `eta = 1/2`,
  `lambda = log 2`, `k = 41.26` (0.702 and 0.661 at `k = 41, 42`).

The reading "uniformity is needed at one rare pair" is a fair interpretation, not a theorem.

## 3. Theorem 3.2 (L1 excludes L2), Theorem 3.3, Cor 3.4, Prop 3.5

**Prop 3.1: correct and checked.** `z_x/z_y = i^(k_x - k_y) prod rho_j^(v_j)` in `Q(i)`. Both points are
units mod `q`. `T_q` is cyclic of order `q - (-1/q)`: for `q = 3 mod 4` it is the norm-one subgroup of
`F_(q^2)^*`, and for `q = 1 mod 4` it is `{(a, a^(-1))}` in `F_q^* x F_q^*`.

**Thm 3.2: correct, proved, actual clusters.** A product of `c` distinct primes `> M^(1-eta)` is
`> M^((1-eta)c)` and divides `B`. Then `|B| <= (C/sqrt2) M^b` by Thm 2.1(a).

*What it is.* It is the contrapositive of the proof of `p3lower.md` Thm C (Prop B.1 on a heavy pair),
restricted to near-extremal clusters. Any two hypotheses whose conjunction contradicts `W <= c_0 M log M`
are mutually exclusive on every near-extremal cluster. That is how every proof by contradiction is built.
It is not a barrier and not evidence against either lemma.

**The headline "L1 and L2 are each other's negation" is false.** The relevant places are the summary,
§0 ("It shows that, on near-extremal actual clusters, L1 and L2 are each other's negation"), §5.1 item 4
("each is the other's negation there") and Cor 3.4.2 ("say that near-extremal clusters satisfy `(P_b)` and
also violate it").
* *Per cluster.* Thm 3.2 gives `U => not H` (`U` = `K`-uniform, `H` = has an L2 pair) on near-extremal
  clusters with `M >= M_1`. The converse `not U => H` is neither proved nor claimed. No argument even
  suggests it: a non-uniform cluster may simply have no heavy pair. At `M <= 10` this is what the data show
  (records: `K = 6.6-13.7`, at most 2 odd primes `q <= M^2` in any `B_xy`), though L2 is not meaningful at
  that size. Mutually exclusive is not complementary.
* *As lemmas.* L1 and L2 are universal statements over near-extremal actual clusters. If growth holds at
  constant `c_0`, i.e. no such clusters exist for `M >= M_1`, then **both L1 and L2 are true**. A statement
  and its negation cannot both be true. The write-up itself says this in Cor 3.4.2 ("both statements hold
  vacuously").
* *(R1-arch) vs (R1-arith).* (R1-arch) is `exists b forall cluster: P_b`. Its negation is
  `forall b exists cluster: not P_b`. (R1-arith) is `forall b forall cluster (M >= M_1(b)): not P_b`, which
  is much stronger than the negation. They are contraries, and their conjunction is equivalent to the
  nonexistence of near-extremal clusters for large `M`.

*Consequence for §5.1 item 4.* "L1 and L2 cannot be pursued as independent lemmas about the same
near-extremal clusters in the hope of verifying either one" is a non sequitur. In route 1 one would prove
(R1-arch) from archimedean/lattice structure and (R1-arith) from arithmetic, each about a hypothetical
near-extremal cluster. The final sentence of the same item ("Any progress must take the form (R1-arch) or
(R1-arith)") says exactly this. What is true, and worth keeping, is narrower:
* neither lemma can be *tested* on actual data at large `M`;
* a proof of L1 alone (without growth) would imply that near-extremal clusters, if any exist, have no
  heavy pair, so L2 would then have to be proved vacuously.

**Thm 3.3: correct.** I re-derived `W >= 2(M-1) mubar`: column by column,
`sum_(x<y)|a_x - a_y| <= e M^2/4` by the threshold-cut decomposition. Then `mubar >= mu_xy/K` and
`mu_xy >= 2 log(|B|/t) + log 2`. With L2, the leading constant is
`4(1-eta)(2 - eps/4) = 8(1-eta) - eps(1-eta) >= 8(1-eta) - eps`, and the slack `eps eta` absorbs the
`O(M)` losses.

*Scope of the claim.* "Recovers Theorem C without its hypotheses `X >= M`, `lambda_Q` and `H`" is
technically true but misleading.
* For actual clusters `X` is free.
* In Thm C, `lambda_Q` and `H` are used only to **produce** a heavy pair (`p3lower.md` Thm C, proof steps
  1-3). Thm 3.3 assumes a heavy pair (or the weaker one large `|B|/t`) directly.
* The `round7/README.md` already states Theorem C with `H_1` in place of `H` and notes "the growth sits in
  the hypothesis".

Thm 3.3 is Prop B.1 with `D_xy` replaced by `2 log(|B|/t)`. That is a mild and correct strengthening (one
super-polynomial `|B|/t` suffices; window structure is irrelevant), not the removal of a hypothesis.

**Cor 3.4: items 1 and 3 correct; item 2 correct only once "each other's negation" is removed; item 4
correct.**
* The derivation "(R1-arch) and (R1-arith) imply no near-extremal clusters for `M >= M_1`, hence
  `M <= (2/c_0 + o(1)) log R/loglog R`" is right. It uses `W = log N' <= 2 log R`.
* "Route 1 is not weaker than the growth statement" holds for **every** route whose lemmas are universally
  quantified over near-extremal actual clusters, by vacuous truth. It is not a defect specific to route 1.
* The round-7 formulation of L1 ("a valid input forcing ...") was meant to be proved inside a class
  (e.g. P3-all\*), where it need not be vacuous. The actual-cluster reformulation loses that, and the
  write-up should say so.

**Prop 3.5.**
* (a) is correct. The count uses `stability.md` Prop 4.1(e) with `|beta|^2 >= 2B^2 > M^(2(1-eta)k)`.
* (b) is correct. The points lie in one coset of `T_q`, `E(M, n_q) >= M^2/(2n_q) - M/2`, and
  `|Q_eta|/(M-1) -> 0`.
* (c) is a correct counting statement. "Cannot be prescribed" means "not every prescription is
  realisable"; at most `P` of the `>= e^((1+o(1))M)` vectors are.
* (d) is labelled heuristic, but its bold conclusion ("L2 cannot follow from genuineness of the labels plus
  profile structure alone ... it must use the short-arc condition") contradicts the paragraph's own caveat:
  "whether such adapted rows can also be near-Plotkin-tight ... is open". Near-Plotkin-tightness **is**
  profile structure (`stability.md` §4). If adapted rows cannot be near-Plotkin-tight, profile structure
  plus genuineness might force L2.
  * Correct reading: "genuineness of the labels alone does not force L2 in the independent-label model;
    whether genuineness plus near-Plotkin profile structure does is open; the short-arc condition is the
    only other link."
  * The summary's "Heuristic: L2 must use the short-arc condition" should be weakened accordingly.

**Status paragraph of §3.**
* "`B_xy` of size `M^(Theta(log M/loglog M))`" is only the **minimum** size of an L2 pair's integer. The
  a-priori range is up to `|B| <= C e^(W/4) = M^(O(M))`.
* "As `p3lower.md` Remark 3 to Theorem C notes, this is beyond equidistribution of Gaussian primes even
  under GRH." That remark is explicitly labelled "(Heuristic remark; nothing depends on it.)" and should
  be cited as such.

## 4. Sections 5-7 (identities, relations, barriers)

**Prop 5.1: correct.**
* (a) is transitivity of congruence. As a statement about the cocycle *at a single `q`* it is right. The
  cocycle is also an exact identity in `Z[i]` with archimedean content, which (a) does not address. The
  sentence "its only consequence ... is transitivity" should keep the qualifier "at an odd prime `q`".
* (b) is `outside.md` Prop C.1. I re-derived the circumradius step and the triple-range identity
  `e - (max - min)`. `records2.py` checks `4(2 Area)^2 = Norm(g_xyw)^2 |b_xy|^2 |b_yw|^2 |b_xw|^2` exactly
  on all 295 triples of the four records.
* (c) is elementary.

**Prop 5.2: correct.**
* The rank `M - 3` is right (three points of a circle are not collinear).
* The covolume formula `covol(ker) = sqrt(M) J/covol(phi(Lambda_0))` is right.
* Hadamard is right, and in fact gives `J <= (1/2) sum |z_x - zbar|^2`, a factor 2 better than stated.
* Minkowski, and the Taylor expansion for the `s_x`-moments, are right (with `|s_x| <= C/2`).

`records2.py` computes the integer relation lattice of the records by LLL on a weighted embedding:
* the ranks are `M - 3`;
* the covolumes (`2.2e4` for the 9-point cluster) are far below the bound (`5.8e7`), because
  `covol(phi(Lambda_0))` is large;
* the shortest relations have `|c|_2 = sqrt 8` (9-, 10-, 8-point) and `6.93` (7-point), against Minkowski
  bounds 35-106.

The remark that the content is lattice-ness, not circularity, is right.

**Lemma 5.3 and the 4-point remark: correct.**
* `v_(pi_j)(c_(x0) z_(x0)) = a_(x0 j)` because `p_j` does not divide `c_(x0)` implies `pi_j` does not
  divide `c_(x0)` (`c` rational).
* The high side uses `conj(pi_j)^(e_j - t + 1)`.
* For four sorted values, the triple ranges are `a_4 - a_2`, `a_4 - a_1`, `a_4 - a_1`, `a_3 - a_1`; they are
  all equal iff `a_1 = a_2` and `a_3 = a_4`.

No violation on 155 (relation, column, threshold) checks over the LLL relations of the records.

**§5.1 Barrier summary.**
* Item 1 is right once restricted to integrality of the reduced chords (§1 above).
* Item 2 is a fair reading. The cited `growth/ptolemy.md` Thm C1 and Thm C (bounded-residue families at
  `k = 5, 6`) match.
* Item 3: the heuristic must be weakened (Prop 3.5(d) above).
* **Item 4 is wrong** (§3 above). It should read: "On a near-extremal cluster with `M >= M_1`, `K`-uniformity
  and an L2 pair cannot coexist (Thm 3.2). So L1 and L2, as universal lemmas, are jointly equivalent to the
  nonexistence of near-extremal clusters for large `M`. Neither can be tested on actual data. Each may be
  pursued separately as a statement about hypothetical near-extremal clusters."

**§6 Open points.**
* Point 3 is answered by Theorem A (§2 above) and should be removed or replaced by the point-total question.
* Point 5's reading of the data is not supported (§5 below).
* Points 1, 2 and 4 are reasonable.

## 5. Independent computations

All Gaussian arithmetic is exact (Python integers). Logarithms and angles are floating point and are used
only for `C`, margins and statistics. **Everything in this section is evidence (exact computation), not
proof.**

| script | what | result |
|---|---|---|
| `refgauss.py` | own toolkit: Euclid gcd, Gaussian primes via `sqrt(-1) mod p`, valuations, windows | -- |
| `records.py` | 9-pt (`C = 7.8804`), **10-pt (`C = 8.2393`, omitted by the author)**, 7-pt, golden-8 records: Thm 1.1 formula, exact (1.1), gcd formula, bisector parallelism, congruences mod `q^a` (`q <= 200`, `a <= 3`, normalised and original points), Prop 3.1 code (`q <= 60`) | all pass; table values reproduced (9-pt `A = 1.056`, `mubar = 0.803`, `max = W/2 = 10.443`, `K = 13.01`; 7-pt `K = 6.56`; golden `K = 7.61`); 10-pt `K = 13.74`, 2 complementary pairs |
| `scan_ref.py` | independent re-scan, different method (numpy generation, octant-midpoint canonicalisation) | 2,621,071 admissible `N`; 22,907 clusters (22,802 + 105); **identical cluster set** to `s9.pkl`, max `|Delta C| = 1.6e-13` (`compare_1e9_out.txt`) |
| `stats_ref.py` | 137,862 pairs | Thm 1.1: 0 failures; bisector: 0; Lemma 1.3: 0 violations, min slack 4.16. `min C`, `min A`, median `K` (3.00, 4.38), percentiles, complementary fractions (0.684, 0.314) reproduced |
| `stats_ref.py` | "K without complementary pairs" | author's statistic (max over non-complementary pairs / mean over **all** pairs) reproduces 1.17 and 3.28. Natural versions: clusters with **no** complementary pair have median `K` 2.67 (`M = 4`, `n = 7194`) and 3.54 (`M = 5`, `n = 72`). With complementary pairs removed from both max and mean: 5.00 and 4.17, with 8,372 (`M = 4`) and 12 (`M = 5`) clusters whose remaining mean is `<= 0` |
| `neg_mubar.py` | clusters with `mubar <= 0` | 26 four-point clusters with `sum mu = 0` exactly (to rounding), `A` from 2.27 |
| `golden7.py` | golden member `n = 7`, `N' = 55728887765525`, own point generation | `C = 5.8977`, `K = 7.375`; **24 of 28 margins identical to the `n = 4` member** (`-2.73 ... 2.73`), only 4 grow (12.61 and three of `W/2 = 15.83`) |
| `records2.py` | Thm A vs Prop 2.2; area identity; LLL relations; Lemma 5.3 | §2 and §4 above |

**Consequences for §4 of the write-up.**
1. *Reading 2 is false.* "In the fixed-`M` families `K` tends to a constant because every margin grows
   linearly in `W`". In the golden family 24 of the 28 margins are constant along the family and only 4
   grow. `K -> 28/4 = 7` because the mean is dominated by the 4 growing pairs. This is the extreme
   *non-uniform* ("adapted") pattern: large margins on few pairs, bounded margins elsewhere. The author's
   own `known_out.txt` shows `min_mu = -2.73` constant along the family, and `-1.76` along Pell.
2. *"Without complementary pairs the median `K` drops to 1.17 (`M = 4`) and 3.28 (`M = 5`)"* misdescribes
   the statistic. Measured in either natural way, removing complementary pairs leaves the margins clearly
   non-uniform: median 2.67/3.54 among clusters with no complementary pair, and 5.00/4.17 if those pairs
   are removed from the mean too.
3. Consequently the summary's evidence bullet and open point 5 are unsupported. They imply that
   small-cluster non-uniformity is a symmetry effect confined to bounded `M` by Lemma 1.3. Correct reading
   (evidence only, `M <= 10`): actual small clusters are non-uniform with or without complementary pairs.
   In the families, all growth of `W` sits on a bounded number of pairs. This is the "adapted margins"
   picture of `p3lower.md` Theorem A, the opposite of what L1/(R1-arch) asserts. Nothing here bears on
   `M -> infinity`, as the write-up rightly says.
4. Minor items.
   * `n = 22784` in `margins_detail_out.txt` vs 22,783 here: one cluster with `mubar` at rounding level.
   * `K(max) = inf` in `stats_s9_out.txt` reflects the 19-26 clusters with `mubar = 0`; the write-up does
     not mention them.
   * The 10-point record (`C = 8.2393`) should be added; it behaves like the 9-point one.

## 6. Issues list

No fatal issue: no displayed theorem is false, and no growth or decisive barrier is claimed.

**Major.**
* **M1 (logic, headline).** "On near-extremal actual clusters L1 and L2 are each other's negation"
  (summary; §0; §5.1 item 4; Cor 3.4.2 wording) is false.
  * Thm 3.2 gives mutual exclusivity on each near-extremal cluster with `M >= M_1`, not complementarity.
  * If growth holds at `c_0`, both lemmas are true (vacuously).
  * (R1-arith) is strictly stronger than not-(R1-arch).
  * The derived claim that the two cannot be pursued as independent lemmas is a non sequitur, and
    contradicts the next sentence of §5.1 item 4.
* **M2 (superseded result, wrong status claim).** Prop 2.2's transfer to `lambda_max(D)`, the "window
  `omega(M log M) <= lambda_max(D) <= O(M^(3/2) log M)`", summary item 6's "ceiling", and open point 6.3 are
  all superseded by `p3lower.md` Theorem A. That theorem is refereed, valid for actual clusters, and cited
  in the write-up. It gives `lambda_max(D) <= (c_0/2) M log M + (M-1) c_C^+` on near-extremal clusters, so
  the window is empty and open point 6.3 is answered. Prop 2.2 survives only as a bound on `max_x D_x`.
* **M3 (evidence misreported).** Two statements in §4 are false.
  * "Without complementary pairs median `K` drops to 1.17/3.28": the statistic keeps the complementary
    margins in the mean. Natural versions give 2.67/3.54 and 5.00/4.17.
  * "Every margin grows linearly in `W`" in the families: 24 of 28 golden margins are constant.

  The interpretation built on them (non-uniformity is a symmetry effect confined by Lemma 1.3) is
  unsupported, and is repeated in the summary and in open point 5.

**Minor.**
* **m1.** Thm 1.1 is the pair case, with equality, of `sharp.md` Prop 1.2 / `p3lower.md` Def 1.1. Credit
  them. `B_xy` clashes with `p3lower.md`'s `B_xy` (the Gaussian `zhat_x/g_xy`).
* **m2.** Thm 2.1 is a one-line rewrite of (1.1) + (F). In (b), the error term depends on `a + b`.
  Cor 3.4.1 drops `M >= M_0(C)`. At `M = 4`, 26 scan clusters have `mubar = 0` (some with `A < 3.5`), so L1
  fails for every `K` there.
* **m3.** Thm 3.3's "without `X >= M`, `lambda_Q` and `H`" is misleading: it assumes the heavy pair (or a
  large `|B|/t`) that `H` and `lambda_Q` were used to produce. It is Prop B.1 with `D_xy -> 2 log(|B|/t)`.
* **m4.** Thm 3.2 and Cor 3.4 are the contrapositive of Thm C's proof on near-extremal clusters.
  * "Route 1 is not weaker than growth" holds for every route whose lemmas quantify over near-extremal
    actual clusters.
  * The write-up should say that the actual-cluster reformulation discards the round-7 intent: L1 was to
    be derived from valid inputs in a larger class.
* **m5.** "Lemma 1.3 is the only known upper bound on an individual margin" (summary, without "that I
  know"). The triple identity gives `mu_xy <= W/2 - 2 log 2 + 4 log C`. Lemma 1.3 is the only one with
  `M`-dependence.
  * "Integrality gives only lower bounds on `|B_xy|`" is contradicted by Lemma 1.3 itself, an integrality
    argument giving an upper bound on `|B|/t`. Restrict it to divisibility of the reduced chords.
* **m6.** Prop 3.5(d): the bold heuristic conclusion contradicts its own caveat (near-Plotkin-tightness of
  adapted rows is open, and is profile structure). The summary's "L2 must use the short-arc condition"
  overstates even the heuristic.
* **m7.** The §3 status paragraph:
  * the size `M^(Theta(log M/loglog M))` is only the minimum for an L2 pair;
  * the "beyond GRH" claim cites a remark that `p3lower.md` labels heuristic.
* **m8.** Prop 5.1(a)'s "only consequence" should keep "at a single odd `q`". The cocycle is an exact
  `Z[i]` identity with archimedean content.
* **m9.** The 10-point record (`C = 8.2393`) is omitted from §4. The 26 clusters with `mubar = 0`
  (reported as `K(max) = inf`) are not mentioned.

## 7. Corrected statement

*Status: partial. No growth improvement and no new barrier. All items are for normalised actual clusters
unless labelled.*

1. **Thm 1.1 (proved, actual).** `beta_xy = i^(k_x) gamma_s B_xy`, with `B_xy = D_s(A_(a_x - a_y))` a
   nonzero rational integer. This is the pair case, with equality, of `sharp.md` Prop 1.2. For odd `q` not
   dividing `N`, `z_x = z_y mod q^a` iff `q^a | B_xy`, and
   `mu_xy = 2 log|B_xy| + log nu_s - 2 log t_xy` exactly. Window collisions are agreements of the
   discrete-log code `w_x(q) = <a_x, ell(q)> + k_x u(q)` (Prop 3.1, proved).
2. **Thm 2.1 (proved, actual).** For `M >= M_0(C)`, L1 on near-extremal clusters is equivalent to
   `max_(x,y) |B_xy|/t_xy <= M^(O(1))`, i.e. to every pair gcd having norm `>= N^(1/2) M^(-O(1))`. This is a
   one-line consequence of (1.1) and `mubar = Theta(log M)`.
3. **Thm 3.2 (proved, actual).** On a near-extremal cluster with `K`-uniform margins, every `B_xy` has fewer
   than `K c_0/(4(1-eta)) + o(1)` prime factors `> M^(1-eta)`. So `K`-uniformity and an L2 pair cannot
   coexist on any near-extremal cluster with `M >= M_1`.
   * This is the contrapositive of the proof of `p3lower.md` Thm C. It makes L1 and L2 **mutually
     exclusive on existing clusters**, not each other's negation.
   * As universal lemmas, L1 and L2 are both true if growth holds at `c_0`, and their conjunction is
     equivalent to the nonexistence of near-extremal clusters for large `M`.
4. **Thm 3.3 (proved, actual).** `K`-uniform margins and one pair with `|B_xy|/t_xy >= M^f` give
   `W >= (2(M-1)/K)(2f log M + log 2)`. This is Prop B.1 with the pair demand replaced by `2 log(|B|/t)`.
   The heavy pair is assumed, not derived.
5. **Cor 3.4 (proved, actual).** Route 1 for actual clusters splits into two contrary statements:
   * (R1-arch): all `|B_xy|`, `1/t_xy` are `M^(O(1))`;
   * (R1-arith): some is super-polynomial.

   Their conjunction is the growth statement at constant `c_0` (true of any route phrased about actual
   near-extremal clusters). Each half may be pursued separately; neither can be tested on data. L2 is one
   mechanism for (R1-arith).
6. **Lemma 1.3 (proved, actual; new).** Every pair has `M <= sqrt2 C^2 |g_xy| + 2`, so
   `mu_xy <= W/2 - 2 log M + O(log C)`. Complementary pairs need `M <= sqrt2 C^2 + 2`. This is the only
   known **`M`-dependent** upper bound on an individual margin, far from the `O(log M)` that L1 needs.
7. **Point totals (Prop 2.2, proved).** On near-extremal actual clusters,
   `max_x D_x <= (1/2) sqrt(c_0(c_0-2)) M^(3/2) log M (1 + o(1))`. For the alternative form of L1,
   `p3lower.md` Thm A already gives `lambda_max(D) <= (c_0/2) M log M + O(M)` on near-extremal clusters. So
   that form, too, holds only vacuously: it is equivalent to growth.
8. **Other results (proved).**
   * Averaged (Markov) uniformity is automatic (Prop 2.3(a,b)). The Poisson rarity of L2 pairs is a model
     statement (Prop 2.3(c)).
   * Exact relations `sum c_x z_x = 0`, `sum c_x = 0`, with `|c| <= M^(c_0/2 + 1/2 + o(1))` exist (Prop 5.2),
     and their content is lattice-ness.
   * Cocycle/area/4-point identities reduce to profile + `B_xy` + residue partitions (Props 5.1, 5.3).
9. **Heuristic.** In the independent-label model, genuineness of the labels alone does not force L2 for
   rows adapted to the labels. Whether genuineness plus near-Plotkin profile structure forces it is open.
   The short-arc condition is the only other link.
10. **Evidence (exact computation, `M <= 10`).** On the 22,907 scan clusters (reproduced independently) and
    on the 7-, 8-, 9- and 10-point records, Thm 1.1, the bisector step and Lemma 1.3 hold with no
    exception.
    * Small actual clusters have non-uniform margins, with or without complementary pairs. The median `K`
      among clusters with no complementary pair is 2.67 (`M = 4`) and 3.54 (`M = 5`); records 6.6-13.7.
    * In the golden family, 24 of 28 margins stay bounded while 4 grow linearly in `W`: the adapted-margin
      pattern, not uniformity.
    * None of this bears on `M -> infinity`.

## 8. Files

In `referee_L1L2_checks/`:

| file | contents |
|---|---|
| `refgauss.py` | own toolkit |
| `records.py` → `records_out.txt` | record clusters |
| `records2.py` → `records2_out.txt` | Thm A vs Prop 2.2, area identity, LLL relations, Lemma 5.3 |
| `scan_ref.py` → `scan_ref_1e9.log`, `ref_1e9.pkl` | independent re-scan to `N = 10^9` |
| `compare.py` → `compare_1e9_out.txt` | identity with the author's cluster set |
| `stats_ref.py` → `stats_ref_out.txt` | pair checks and uniformity statistics |
| `neg_mubar.py` → `neg_mubar_out.txt` | clusters with `mubar = 0` |
| `golden7.py` → `golden7_out.txt` | golden member `n = 7` |
| `notes.md` | referee working notes |
