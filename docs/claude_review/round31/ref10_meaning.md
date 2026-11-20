# Round 31, Part X (Question A/B thread, "astra_free_X"): analyst block "meaning"

Scope. I read all of Part X of `astra_free_X/answer.md` (lines 2317-2498), `final10.txt`, the coordinator's question
`reply10.txt`, `README_IX.md` §§0, 3, 5, and the "Not everything is capped" paragraph of `expert_note.md` §1. I did
not re-referee the proofs of Part X; where I state something as **Proved** below it is my own short argument, given
in full. My scripts are in `ref10_meaning/` (own code, exact or 80-digit arithmetic; nothing of Astra's imported):
`brute.py` (full enumeration of lattice points, s <= 22), `mitm.py` (meet-in-the-middle minimal angle gap, s <= 36),
`paircheck.py` (genuine pair collision moduli from residues). Outputs: `brute_out1.txt`, `mitm_out_pi4.txt`,
`mitm_out_pi4_b.txt`, `mitm_out_5mod8.txt`, `paircheck_out.txt`.

Notation. `q_P` is the product of all primes `p = 1 mod 4` in `[P, 2P]`; `s = omega(q_P)`; the circle is
`|z|^2 = q_P`, radius `R = q_P^{1/2}`, so a critical arc has length `C R^{1/2} = C q_P^{1/4}`. `M_max(q_P, C)` is the
largest number of lattice points on one arc of that length. `theta_p = arg pi_p` for a fixed factor `pi_p` of `p`.
All other notation is that of `README_IX.md` (Definitions paragraph) and `sharp.md` Def. 1.1.

## 0. Results in eight lines

1. **Part X does not answer P3-all\*.** It delivers a partial fake (X.2), a finite-cutoff rigidity for bounded-height
   lifts (X.1), a negative result on one fake-construction strategy (X.3) and a weight bound under a hypothesis that
   genuine witnesses violate (X.4). None of these is progress toward a P3-all\* **lower bound**; X.2 is progress
   toward a **fake**.
2. **Honest balance:** P3-all\* is now more likely capped (my estimate 70%) than carrying a lower bound (<= 2%);
   the rest is "open for ever".
3. **The coordinator's point (2) is correct after one correction:** on genuine witnesses every e-all test (every
   monomial, not only label differences) is passed automatically, with equality before truncation; this is a
   property of the prime set and not of the cluster. The correction: it is equally true of every residue test
   (P3-char included), so it is not a special defect of the difference tests. What is special is that label
   differences are **not cluster words**: the cluster words (mark differences) are already in P3-char, which is
   capped, and on genuine witnesses their test is exactly the X-smooth part of the integrality statement (T1).
4. **Consequence (Heuristic, but I think sound):** any P3-all\* lower bound must go "non-cluster tests => rigidity of
   witnesses => integral (T1) data on cluster words => rate on a fixed parent". The last step is question (A).
   So the residue route closes in the sense the previous judge intended: its only payoff runs through (A) or the
   problem itself.
5. **Computed (the most useful item):** on `q_P` the truth is **M_max = 1 for every C <= 8 once P >= 200**, with
   enormous margin: the closest two lattice points are at arc distance `468, 4.6e3, 3.3e4, 6.6e5, 4.8e6, 8.1e7`
   times `q_P^{1/4}` for `P = 200, 250, 300, 350, 400, 450` (s = 16, ..., 32). The minimal gap matches the random
   heuristic `(pi/4) 3^{-s}` in angle to within a factor 7 in every case. The known bound on these circles is of
   order `s/3`.
6. **(A1) on `q_P` is, as far as I can see, not approachable by known methods:** it is a lower bound at the
   Liouville threshold for linear forms in `s -> infinity` logarithms of Gaussian primes with coefficients in
   `{0, +-1}`; Baker-type bounds miss by a factor `exp(c^s)` in the exponent.
7. **Recommendation:** stop the residue route for internal work now (the judge's condition is met in substance).
   Ask this session **one final question**, (A) in the graded form of §3.3, and then stop the thread.
8. **The growth rate `(2/3 + eps) log R/loglog R` and the uniform bound are unchanged by Part X.**

## 1. What Part X settles (question (1))

### 1.1 Classification of each result

| Item | Status as reported by the refereeing blocks (not re-checked here) | Toward a P3-all\* fake | Toward a P3-all\* lower bound | Meaning |
|---|---|---|---|---|
| X.1 finite rigidity | Proved by Astra | no | **no** (see below) | Moves the "problem itself" end of the residue column from the limit `X -> infinity` (IX.7) down to polynomial cutoffs, for lifts of polynomial height, on parents with bounded `T = max v_2(l - 1)` |
| X.2 partial fake | Proved by Astra | **yes, partial** | no | Every single-label ray, every pair sum, all words with a nonzero class sum (a `1 - o(1)` fraction of each bounded length), all powers, all units, all cutoffs, plus P3-char, on one parent at `W = O_{C,H}(m log m)` |
| X.2, 5 mod 8 version | Proved by Astra | yes (it shows what a fake needs) | it blocks one route | Lift heights must diverge with the cutoff on a passing parent: the height hypothesis of X.1 is **not implied** by e-all |
| X.2, deliberate failure | Proved by Astra | neutral | neutral | A choice of witnesses, not a forced failure: it does not show the same-class differences cannot be repaired |
| X.3 finite-lift dichotomy | Proved by Astra | no (closes one fake strategy) | no | IX.3's owner design admits no short exact witnesses; irrelevant to fakes that use no lift (X.2) |
| X.4 weight bound | Proved by Astra | marginal: a design constraint on fakes | no | Its hypothesis (bounded joint rank) is **false for genuine witnesses** whenever `P <= m^{O(1)}` (§1.2), so it is a theorem about a class of fakes |
| X.5 least cutoffs 439, 3881, 8681, 5783 | Computed by Astra | no | no | Calibration of X.1 |

Why X.1 is not progress toward a lower bound. The alternative posed in `reply10.txt` was: a height-quantified IX.7
**followed by** a lower bound `W >= f(M) M log M` for witnesses that pass. X.1 is the first half. The second half
needs the passing witnesses to have bounded-height lifts, and the 5 mod 8 version of X.2 proves the opposite: on one
fixed parent passing all single-label rays at every cutoff, each finite-cutoff rational lift at label `p_j` has twist
norm `> (X/(sqrt2 p_j))^{1/2}`. So the height cannot be moved to the left-hand side, and X.1 + X.2 together show that
the height-quantified route of `README_IX.md` §5 item 1 cannot by itself produce a lower bound. Quantifiers: in X.1
the parent, the label `p`, and the height bound `H` are fixed before the cutoff `X`; in X.2 the arc constant `C` and
the word length `H` are fixed before the weight constant, the parent sequence is then built, and one witness family
works at every cutoff `X` simultaneously.

### 1.2 X.4 on genuine witnesses (Proved, my argument)

Let the labels be split primes in `[P, 2P]` with genuine witnesses `r_j = pi_j`. If labels `j != k` have equal norm-one
ratio tuples at all inert `l <= x`, then every such `l` divides `pi_j conj(pi_k) - conj(pi_j) pi_k = 2i Im(pi_j
conj(pi_k))` (rational primes `l = 3 mod 4` dividing `w - conj w` divide `Im w`), and `0 < |Im(pi_j conj pi_k)| <=
sqrt(p_j p_k) <= 2P`. Hence `theta_3(x) <= log(2P)`. So if `theta_3(x) > log(2P)` all `s` tuples are distinct,
`D(x) = s`, and (X.4) reads `W >= 0`. Moreover distinct tuples generate a subgroup of size `>= s` inside a group of
exponent `L(x)`, so `d >= log s / log L(x)`; with `x = K log m`, `log L(x) = o(x)` (Astra's lemma) and
`P <= m^{K/2 - 1}` this gives `d_K -> infinity`. So the hypothesis of the "restricted growth rate" of X.4 fails on
every genuine parent with polynomially bounded labels: X.4 never applies to an actual circle in the regime of interest.
This is what Astra's sentence "principally forces the support itself to have large weight" means, stated plainly.

### 1.3 Capped or lower bound? (Heuristic)

After X.2 the only tests not shown passable by local, non-liftable witnesses are the balanced words inside a class,
the extreme case being the same-class differences `e_j - e_k`. What a fake needs there is a "residue code": label
tuples in `prod_{l^a <= X} (norm-one group mod l^a)` such that for every pair the product of the prime powers at which
the two tuples collide (at a unit) is at most about `sqrt(p_j p_k)` times the available sine budget, compatibly with
the block sums of X.2. **Such codes exist**: the genuine residues are one (§2.1; collision moduli are divisors of
integers of size `<= sqrt(2 p_j p_k)`). The difficulty is only to combine a code of genuine-like pair structure with
the non-genuine singleton data that X.2 needs and that X.1 forces (genuine residues at all primes below `sqrt(2p)`
pin the singleton angles, X.1 "direction pinning"). No obstruction to this combination is known; X.2's "limit of the
order argument" excludes only a *common* index bound, not per-pair bounds. My estimates: a full P3-all\* fake at
`3M^A`, `W = O(M log M)`, exists: 70%; a lower bound `W >= f(M) M log M` from P3-all\*: <= 2% (the previous judge
gave 2% for a rate; after X.1-X.2 I would put it at 1%).

## 2. The structure revealed (question (2))

### 2.1 Genuine pair collisions (Proved; Computed check)

**Proposition.** Let `pi_j, pi_k` be Gaussian primes of norms `p_j != p_k`, `w = pi_j conj(pi_k)`, and for `n >= 1`
write `w^n = x_n + i y_n`. For a rational prime power `l^a` with `l` odd, `l` not dividing `p_j p_k`, the ratio
`rho = (pi_j/conj pi_j)(pi_k/conj pi_k)^{-1} = w/conj(w)` satisfies `rho^n = i^u mod l^a` if and only if `l^a`
divides `D_u(w^n)`, where `D_0 = y_n`, `D_2 = x_n`, and `D_1, D_3` are `x_n - y_n`, `x_n + y_n` (up to order, depending
on the orientation convention). These four integers are nonzero, of absolute value `<= sqrt2 (p_j p_k)^{n/2}`, and the
strengthened sine inequality (1.2) of `sharp.md` holds for the genuine angle with equality before truncation:
`|sin(n(theta_j - theta_k) - u pi/4 ...)|` equals `|D_u(w^n)|/(sqrt2^{[u odd]} (p_j p_k)^{n/2})`.

*Proof.* `rho^n = i^u` mod `l^a` iff `l^a | w^n - i^u conj(w)^n` (the denominator `conj w^n` is a unit at `l`). With
`w^n = x + iy`: `w^n - conj w^n = 2iy`, `w^n + conj w^n = 2x`, `w^n - i conj w^n = (x - y)(1 - i)`,
`w^n + i conj w^n = (x + y)(1 + i)`. A rational odd `l^a` divides `c(1 +- i)` with `c in Z` iff it divides `c`.
Nonvanishing: `w/conj w` is not a root of unity, since `pi_j` divides `w` and does not divide `conj w`. The sine
identity is `Im(e^{-i u pi/4} w^n) = |w|^n sin(n arg w - u pi/4)`. `square`

So the collision modulus of the ray `n(e_j - e_k)` at genuine witnesses is the X-truncated odd exterior part of a fixed
nonzero integer `D_u(w^n)` attached to the pair of primes, and the test is automatically passed. The coordinator's
"odd part of `Im(pi_j conj pi_k)`" is the unit-`+-1` case at `n = 1` (Im for `u = 0`, Re for `u = 2`); for the other
two units it is `x -+ y`. `paircheck.py` (P = 100, 10 labels, all 90 ordered pairs, 4 units, exterior prime powers
up to 3000) recomputed the collision modulus from the residues `pi_j mod l^a` and compared it with the integer side:
**360 cases, 0 mismatches**; the largest ratio `log(modulus)/log sqrt(p_j p_k)` was 1.068 (< `log(sqrt2 sqrt(p_jp_k))`)
(Computed).

The same argument applies to every monomial `v`: at genuine witnesses the e-all collision modulus of `v` is the
truncated odd exterior part of `D_u(prod pi^{v_+} conj(pi)^{v_-})`, with equality in the sine inequality before
truncation. **So genuine witnesses pass all of e-all at every cutoff, on every parent, whether or not the parent
carries a cluster** (Proved). This is the correction to the coordinator's formulation: "e-all for differences of labels
carries no information about clusters" is right, but it holds for the whole of e-all and for P3-char as well; the
residue tests are always constraints that genuine data meet for free, and their use for a growth proof is always
indirect (they must exclude the non-genuine angle models that carry large clusters).

### 2.2 Where the cluster enters, and why the residual question reduces to (A) (Heuristic, with one Proved step)

*Proved.* For a genuine cluster in one unit class on a squarefree parent, the mark-difference word `a_i - a_j`
(coefficients in `{0, +-1}` on the set `D_ij` of primes where the two marks differ) is a zero-sum mark character, hence
tested in P3-char. Its genuine collision modulus is the truncated odd exterior part of the integer `t_ij = Im X_ij` of
(T1) (`X_ij` the oriented product over `D_ij`), and the sine inequality reads
`C q^{-1/4} >~ |sin(Phi_i - Phi_j)| >= (X-smooth part of t_ij)/sqrt(Q_ij)`. So the cluster-word tests are exactly the
X-smooth shadow of (T1). P3-char is capped (IX.3): the X-smooth shadow of (T1) does not carry a rate.

*Heuristic.* P3-all\* adds to P3-char only the tests on non-cluster words (label rays, label differences, mixed words).
For genuine witnesses these are pair or word properties of the prime set (§2.1). Their only possible role in a growth
proof is rigidity: to force the witnesses near enough to genuine that the cluster words carry not only the X-smooth
part of `t_ij` but its integrality. X.1 and IX.7 show that when full rigidity is reached the remaining task is the
problem itself; and when only integrality of the `t_ij` on one parent is reached, the remaining task is precisely
question (A): turn `binom(M,2)` near-square subproducts `Q_ij = x^2 + t_ij^2` with tiny `t_ij` and the triangle
cocycle into `M = o(s)` on a fixed parent, which (T2)/(T4) of the other thread do only on random parents.

So the residue tests split as follows, with quantifiers.
* **(a) Capped**, at `W = O(M log M)`, one parent sequence, witnesses fixed before all cutoffs: P3-char at every
  cutoff (IX.3); all single-label rays, pair sums, all bounded words with a nonzero class sum, with all powers and
  units (X.2, `H` fixed before the weight constant).
* **(b) The problem itself**: e-all with one rational lift per label at all cutoffs (IX.7); and at polynomial cutoffs
  for lifts of polynomial height on parents with `T` bounded (X.1, `H` fixed before `X`).
* **(c) Residual**: whether non-liftable local witnesses can also pass the same-class balanced words, in particular
  the differences `e_j - e_k`, at cutoffs `3M^A` with `A` chosen after the weight constant.

**Can resolving (c) change the growth rate for genuine circles?** Logically yes, in one direction only: a proof that
(c) fails at `W = O(M log M)`, i.e. `W >= f(M) M log M` with `f -> infinity` for every witness system passing e-all at
`3M^A`, would apply to genuine circles (which pass e-all at every cutoff, §2.1) and give `M = o(log R/loglog R)`. A fake
for (c) changes nothing for genuine circles. By §2.2 such a lower-bound proof would have to contain a fixed-parent
rate argument at least as strong as (A1) on every parent, including `q_P`, where (§3) the truth is `M = 1` and no
method is known. I put the probability of that scenario at <= 1%. In this sense the residue route is closed as the
previous judge intended, although no complete P3-all\* fake was built.

## 3. Question (A) on a fixed parent (question (3))

### 3.1 Exact computation on `q_P` (Computed)

Method. Every lattice point on `|z|^2 = q_P` is `i^k prod_p pi_p^{eps_p}` (`eps_p` meaning `pi_p` or `conj pi_p`),
`4 * 2^s` points. Two points differ in angle by `2 sum_p c_p theta_p mod pi/2` with `c in {0, +-1}^s`, `c != 0`
(the set `D` of primes with `c_p != 0` is the set where the two points differ), so the minimal arc distance is
`2 R min_c dist(sum c_p theta_p, (pi/4) Z)` (chord and arc agree to relative `O(gap^2)`). For `s <= 22` I enumerated all
points and slid windows of angular width `C q^{-1/4}` (`brute.py`, float64, closest pair re-verified at 60 digits). For
`s <= 36` I computed `min_c` by meet-in-the-middle on the `3^{s/2}` half-sums, with angles stored as `uint64` fixed
point (period `pi/4` mapped to `2^64`, so reduction mod `pi/4` is exact wrap-around; accumulated error
`<= s * 2^{-64}` of the period, at least 300 times smaller than every minimum found), keeping the ten best candidates
and recomputing each at 80 digits from its exact sign vector (`mitm.py`). The minimiser is printed as an exact Gaussian
integer `X_D`.

Split primes in `[P, 2P]` (all units):

| P | s | log10 q_P | min arc / q^{1/4} | heuristic `2 q^{1/4} (pi/4) 3^{-s}` | `|D|` of the closest pair | M_max C=1 | C=4 | C=8 | known bound `(2/3) log R/loglog R` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 50 | 5 | 9.3 | 0.759 | 1.37 | | 2 | 2 | 3 | 3.0 |
| 100 | 10 | 21.7 | 3.13 | 6.98 | | 1 | 2 | 2 | 5.2 |
| 150 | 13 | 30.6 | 6.66 | 43.5 | | 1 | 1 | 2 | 6.6 |
| 200 | 16 | 39.7 | 468 | 304 | | 1 | 1 | 1 | 8.0 |
| 250 | 20 | 51.0 | 4,557 | 2,570 | | 1 | 1 | 1 | 9.6 |
| 300 | 22 | 58.0 | 33,327 | 16,000 | | 1 | 1 | 1 | 10.6 |
| 350 | 26 | 70.4 | 662,168 | 247,000 | 17 | 1 | 1 | 1 | 12.3 |
| 400 | 30 | 83.1 | 4,780,072 | 4,540,000 | 19 | 1 | 1 | 1 | 14.0 |
| 450 | 32 | 90.5 | 81,224,068 | 35,000,000 | 23 | 1 | 1 | 1 | 14.9 |
| 500 | 36 | 103.0 | <= 51,430,627 (not certified, see below) | 595,000,000 | 24 | 1 (heur.) | 1 (heur.) | 1 (heur.) | 16.6 |

For `200 <= P <= 450` the entries `M_max = 1` follow from the minimal gap exceeding `8 q^{1/4}` (certified up to the following caveat: the fixed-point error bound `s 2^{-64} pi/4` is below `1.4e-18`, at least 700 times smaller than each minimum, and the ten best fixed-point candidates were all recomputed at 80 digits). At `P = 500`, `s = 36`, the 64-bit precision is insufficient (error bound `1.5e-18` exceeds the minimum `4.5e-19` found; exact/fixed-point ratio 1.18): the value is the exact distance of a real pair, hence only an upper bound for the minimum, and `M_max = 1` there is heuristic. The "known bound" column is
the record's growth bound without the `eps` and the `O(1)`, for comparison only (it is not a rigorous numerical bound).

Primes `5 mod 8` in `[P, 2P]` (the family X.1 singles out as rigid):

| P | s | min arc / q^{1/4} | `|D|` |
|---:|---:|---:|---:|
| 200 | 9 | 25.6 | 8 |
| 400 | 15 | 738.9 | 10 |
| 500 | 19 | 13,346 | 14 |
| 600 | 24 | 860,483 | 15 |

Observations (Computed). (i) The minimal angular gap is within a factor 7 (ratios 0.15 to 2.7) of the random prediction `(pi/4) 3^{-s}`
in every case; there is no visible arithmetic conspiracy. (ii) Already for `P >= 200` no two lattice points of `q_P`
lie within `8 q^{1/4}`; the margin grows like `((1.4P)^{1/4}/3)^s`. (iii) The closest pair always differs on more than
`s/2` primes, as Liouville forces (below). So the truth on these explicit parents is `M_max(q_P, C) = 1` for every
fixed `C` and large `P`, heuristically, against a known bound of order `s/3`.

### 3.2 Feasibility of (A1) (Proved where stated; otherwise Heuristic)

*Reformulation (Proved).* Two points at arc distance `<= C q^{1/4}` exist iff for some `D` and orientation the
Gaussian integer `X_D = prod_{p in D} pi_p^{+-}` (norm `Q_D`) satisfies `|t| <= (C/2 + o(1)) sqrt(Q_D) q^{-1/4}` for
one of its four "axis coordinates" `t in {Re X_D, Im X_D, (Re -+ Im) X_D}` (the last two scaled by `1/sqrt2`). The
integer `t` is nonzero, so `|t| >= 1` (**Liouville**), and every close pair has `Q_D >= (2/C)^2 q^{1/2}`, i.e. `D`
carries at least half the log weight of `q` minus `O(log C)`. With Plotkin's bound for binary codes of minimum distance
about `s/2` this gives the order-`s` bound and nothing better: the Liouville step is exactly critical.

*Heuristic count.* Treating `t/sqrt(Q_D)` as uniform, the expected number of close pairs is about
`sum_D 2^{|D|} (C/2) q^{-1/4} ~ (C/2) 3^s q^{-1/4}`, which is `< 1` once the geometric mean of the primes exceeds
81, i.e. `P > ~58`. The table confirms the prediction quantitatively.

*Why it is not approachable by known tools.* Proving `M_max(q_P, C) <= 2` (or even `M = o(s)`) requires, for all
`3^s` vectors `c`, a lower bound `|sum c_p theta_p - k pi/4| >> q^{-1/4} = exp(-(1/4) sum log p)` for a linear form in
`s + 1` logarithms (`theta_p = (2i)^{-1} log(pi_p/conj pi_p)`) with coefficients in `{0, +-1}`. Matveev's bound is
`exp(-c^s prod_p log p * log s)` = `exp(-c^s (log P)^s)`: useless by a tower in the exponent. Liouville gives
`|t| >= 1`, exactly the critical exponent `1/4` only for `|D| <= s/2`. Nothing between these two is known for
`s -> infinity`. A pair-free statement for one explicit sequence of radii and fixed `C` would be a special case of the
Cilleruelo-Granville conjecture for an explicit sequence with `omega -> infinity`; I know of no such case in the
literature (not searched systematically in this block). The weaker target `M = o(s)` on `q_P` asks to use many pairs at
once (`binom(M,2)` near-squares and `binom(M,3)` triangle identities on one prime set), which is the (T1)-(T2) route
without the entropy of a random parent; (T4) records that counting exponent patterns and cofactor tables cannot do it
when `P <= s^2`, and `q_P` has `P ~ 2 s log s`. My estimate for a proof of `M = o(s)` on `q_P`: 3%; for an explicit
improvement of the constant `1/3` on `q_P` (a bound `M <= (1/3 - c) s`): 5%.

### 3.3 Proposed final question to this session

**(A1, graded).** On the explicit parents `q_P` (all split primes in `[P, 2P]`; or all primes `5 mod 8` there), with
`C` fixed before `P`: (i) prove any bound `M_max(q_P, C) <= (1/3 - c) s` with `c > 0` absolute, or (ii) `M = o(s)`,
or (iii) state precisely which Diophantine inequality on the arguments `theta_p` would suffice and show it is not
implied by the triangle cocycle and the near-square identities alone (a fake on the same prime set). The computation
above shows the truth is `M = 1` with margin `((1.4P)^{1/4}/3)^s`; any proof must use the dense block structure of
the triangles and cannot use randomness. Deliverables: a theorem, or a fake that keeps all of (T1), the cocycle with
exact divisibility, and the signed residues, and drops only the exact norm.

**(A2).** An integral-cofactor fake on `q_P` itself: Gaussian integers `z_1, ..., z_M`, `M >> s` (or `M` larger than
any fixed constant if that is all that can be done), with `z_j - z_i = 2 i eta H_ij t_ij` as in (T1), `t_ij` nonzero
rational integers with `|t_ij| <= (C/2) sqrt(Q_ij) q^{-1/4}`, the triangle identities exact in `Z[i]`, all
divisibilities exact, the signed residues kept; and with `|z_i|^2 = q_P` dropped (state the norm defect). If it exists,
the cocycle is not the missing input on fixed parents; if not, say what the proof used.

I recommend asking (A1, graded) with (A2) as the explicitly acceptable fallback, and labelling it the final question.

## 4. Recommendation (question (4))

(a) **Stop the residue route for internal work: yes.** The previous judge's condition ("if that pass returns another
fake") is met in substance: X.2 is a fake for every test except the same-class balanced words, and §2.2 shows the
remaining tests can only pay off through (A). Do not commission a further P3-all\* pass. Record (c) of §2.2 as open.

(b) **Ask this session one final question, (A1 graded) with (A2) as fallback (§3.3), then stop the thread.** Honest
probabilities of moving the **worst-case** growth rate `(2/3 + eps) log R/loglog R` for general `R`:
* further P3-all\* work: <= 1%;
* (A1)/(A2): <= 1% for the worst-case rate (a rate on `q_P` is a statement about one explicit family, not about all
  `R`; it moves the worst case only if the method extends to all dense profiles); 3-5% for a restricted theorem on
  `q_P`; about 30% for a useful barrier (an integral-cofactor fake), which would cap the cocycle input on fixed
  parents.
* The uniform bound `M(C) < infinity`: < 0.5% from either.

## 5. Proposed wording for the expert note (question (5); recommendations only, nothing written to the repository)

**§1, in the "Not everything is capped" paragraph, replace** "Not capped: P3-all$^*$ (prime-level residues tested on
every monomial at moduli $3M^A$), ..." **by:**
"Not capped: the same-class balanced words of P3-all$^*$ (see below), twisted point collisions ($\sum c\ne0$) at moduli
$M^A$, the 2-adic clause, ... [rest unchanged]. Within P3-all$^*$ (prime-level residues tested on every monomial at
moduli $3M^A$), round 31 Part X.2 caps, on one parent sequence at $W=O_{C,H}(M\log M)$, every power of every single
label, every pair sum and every signed word of length $\le H$ with nonzero coefficient sum in one of $b\to\infty$
label classes (a fraction $1-o(1)$ of the words of each length), at every unit and every modulus at once, together
with P3-char. Open: the balanced words inside a class, in particular label differences $e_j-e_k$. With rational lifts
of height $\le H$, two ray tests at one dyadic exponent force genuine witnesses at modulus $\sqrt2(pH^2)^{n/2}$,
$n=2^{T-1}$, $T=\max v_2(\ell-1)$ (Part X.1; polynomial on the product of the primes $5\bmod 8$ in $[P,2P]$); on a
passing fake the lift heights must diverge with the modulus (Part X.2), so this does not yield a lower bound. Genuine
witnesses pass every e-all test on every parent (the collision moduli divide fixed integers of size
$\le\sqrt2\prod p^{|v_j|/2}$), so the residual tests can only help through rigidity; on genuine clusters the
cluster-word tests are the truncation of the integrality of the (T1) cofactors, which is the fixed-parent question."

**§5 item 12, add a bullet:** "*Part X* (`round31/README_X.md`): finite-cutoff rigidity for bounded-height rational
lifts (polynomial cutoff on the $5\bmod8$ dyadic parent); a partial P3-all$^*$ fake passing all single-label rays and
almost all bounded words at every modulus; no short exact witnesses on the IX.3 design; a weight bound under a
joint-rank hypothesis that genuine witnesses violate. Computed: on $q_P=\prod_{p\equiv1(4),\,P\le p\le2P}p$ no two
lattice points lie within $8q_P^{1/4}$ for $200\le P\le450$ (closest pair $468$ to $8\cdot10^7$ times $q_P^{1/4}$),
against a known bound of order $\omega(q_P)/3$. The residue route is closed for internal work. No change to the growth
rate."

**Sentences to avoid.**
* "P3-all$^*$ is capped" or "P3-all$^*$ is answered" (the same-class balanced words are open).
* "Part X proves a finite-cutoff lower bound" or "X.1 gives a lower bound for P3-all$^*$" (it applies only to
  bounded-height lifts, which passing witnesses need not have).
* "X.4 proves a growth rate" (its hypothesis fails for genuine witnesses with polynomially bounded labels).
* "the omitted difference must fail" (X.2 chose the failure; it does not exclude repairs).
* "e-all on label differences is the only test genuine data pass automatically" (genuine data pass all residue tests).
* "the twist must be $\pm1$" (false even on genuine circles: $\gamma=\bar\pi$, $\gamma=1+i$).
* "$M_{\max}(q_P,C)=1$ is proved" (computed for $P\le450$, heuristic beyond).

## Appendix: quantifiers and cautions

* Computations: each `P` fixes the parent; `C in {1, 4, 8}` (and any `C` below the tabulated ratio) is then tested;
  nothing is claimed uniformly in `P`.
* No Subspace-theorem exceptional sets, no discarded gcd or denominator costs, no common phase set to zero are used:
  the angle-gap computation is over all units and all sign vectors; the common phase is irrelevant to pair gaps.
* §1.2 and §2.1 are complete elementary proofs. §1.3, §2.2 (except the marked Proved step), §3.2 (except the
  reformulation and the Liouville step) and all probabilities are Heuristic.
