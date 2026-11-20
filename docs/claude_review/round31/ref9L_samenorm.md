# Referee report: the same-norm P2 comparison and the Legendre-profile tests

**Block refereed.** Two passages of Part IX of `astra_lemma_IX/answer.md`:
* §IX.0, the three paragraphs that start at "PROVED (a common norm for the P2 test)" (lines 1990-1996);
* §IX.6, the first paragraph, "PROVED (Legendre profiles P2)" (line 2116).

**Also read.**
* `README_VIII.md` §§1.3-1.6, 2 and 3, and the referee reports `ref8_fakecore.md` and `ref8_fakedefs.md` (the passages on the prime window, J1 and the area test). I read `ref8_faketail.md` by grep only, for every use of the prime interval.
* The parallel thread's VIII.3 in full: `astra_free_VIII/answer.md`, lines 1741-1903.
* `README_III.md` §§1 and 3.2, and expert note lines 885-903, for the record's random-parent theorem and its caveat about narrow-bin parents.
* The companion report `ref9L_powers.md` §1, for the corrected statement of the cap-E Gaussian theorem ("Theorem G").

**Scripts.** All scripts are my own and import nothing from Astra. They are in `ref9L_samenorm/`:
* `paley_checks.py`, output `paley_checks.out`: exact integer checks of the normalisation, parity, Sidon and Reed-Muller statements. It also runs a floating-point simulation of the permutation margin.
* `bincount.py`, output `bincount.out`: a segmented sieve that counts split primes in one bin, at l of about 10^4, 10^5, 10^6 and 10^7.

**Literature.** I did no web search. Every claim in this block is internal to the record.

Labels: **Proved** means I checked every step myself. **Computed** means I reproduced the result with my own code. **Heuristic**, **Gap**, **False** and **Not checked** have their usual meanings.

## 0. Verdict

The same-norm comparison is **correct** (Proved). It depends on two results that were refereed elsewhere: VIII.3 of the parallel thread (refereed in Part VIII) and Theorem G (refereed in `ref9L_powers.md`, verdict Proved).

The widened bin changes exactly one step of VIII.3, the pair margin, and that step survives with room to spare. No other step of VIII.3 uses the width of the prime interval. The normalisation claim is correct. Both claims of the IX.6 paragraph are correct.

The corrections below concern wording, quantifiers and constants. None is mathematical.

1. **The fake half holds for every s-subset of the bin, not only for the subset chosen.** Astra says this in passing. It should be part of the statement, because it determines what the result means (§5).
2. **The genuine half is purely asymptotic.** The bins hold fewer than s primes for every l up to 10^7 (Computed, §2.3). Theorem G as refereed needs J >= 1600 log 54, about 6383, at E = 2. Here J ~ 0.24 log l, so the constants as written need log l >= about 26,640 (Computed). That is legitimate for an "arbitrarily large" statement, but it should be said.
3. **Constant in the bin count.** Some bin, namely the average one, holds at least (1+o(1)) l^{1.24} / (2 · 1.99 · log 2 · log l) ≈ l^{1.24}/(2.76 log l) split primes. The lowest bin holds only about log 2 times that (Computed). The bound K >> l^{1.24}/log l is right.
4. **"The older quartet exposure can likewise be used for the genuine-count half": Not checked.** The record's quartet theorem (Part III, III.9) is stated for dyadic pools [r^D, 2r^D]. Its transfer to a bin pool would have to be re-derived. The result does not need it.

**Novelty.** This is new relative to the record. The record states explicitly that its random-parent theorem "says nothing ... about the narrow-bin parents on which fakes are built (they lie inside the exceptional set)" (expert note, line 901). The Part III judge also wrote: "Whether angle fakes with growing m exist on typical dyadic parents is not established in the record; without that, no separation is shown even on typical parents" (`README_III.md` §3.2(iii)).

Astra's comparison gives the first separation of fake from genuine on a common norm. The norms are typical parents of a bin pool. They are not typical dyadic parents: a random dyadic s-subset lands in a single bin with probability about exp(-0.75 s log s). So the record's question about typical dyadic parents is still open in the form the Part III judge asked it.

**Consequence for the integer problem: none.** The growth rate (2/3 + eps) log R/loglog R and the uniform bound are unchanged. The result sharpens the barrier map. The A_{2M} inputs of P2 cannot certify M = o(s) even on norms where the truth is M = O(1), and they cannot tell such norms apart from any other s-subset of the bin.

## 1. Corrected statement (quoted in full)

**Theorem (same-norm comparison). Proved, conditional on VIII.3 (refereed in Part VIII) and Theorem G (refereed in `ref9L_powers.md`).**

*Setting.* Fix a real rho in (0.99, 1). There is l_0(rho) with the following property. Let l >= l_0 be a prime with l ≡ 3 mod 4. Put s = l - 1 and P = l^{1.99}. Then there is an interval I_l = [P_b, P_b e^{sigma}] contained in [P, 2P], with sigma <= l^{-3/4}, that contains K >= (1 + o(1)) l^{1.24}/(2 · 1.99 · log 2 · log l) primes p ≡ 1 mod 4.

*(i) Fake half, for every subset.* Let S be any s-element set of split primes in I_l, and put N = prod_{p in S} p^2. Then there is an angle model with normalised norm exactly N and with m = (rho/2 + o(1)) l distinct marks on an arc of constant one. The model is the Legendre/Paley model of VIII.3: exponent rows 1 + K_i, cap two, the primes of S assigned to the columns I by a suitable permutation. It satisfies:
* both sharp monomial gaps for every nonzero integral exponent vector, including all rational characters of the profile;
* axiom (R) of round 5 at every odd prime power <= 2m and for all four units;
* positive Haar measure of good angles on a realisation component;
* P_+(N) <= 2P and W = log N = (3.98 + o(1)) l log l.

*(ii) Genuine half, for almost every subset.* There is an absolute constant M_1. All s-subsets S of the split primes of I_l satisfy the following, except a fraction at most exp(-(c + o(1)) l log l) of them, where c = 0.24/1000 with the constants of `ref9L_powers.md`. For every C > 0, every arc of length C sqrt(R), with R = sqrt(N), of the circle x^2 + y^2 = N contains at most ceil(C) · M_1 lattice points.

*Conclusion.* Let l run through the primes ≡ 3 mod 4 and choose S as in (ii). This gives arbitrarily large norms N = R^2, with R squarefree, on which both (i) and (ii) hold for the same N.

*Quantifiers.* rho is fixed first. l_0 depends on rho and on the threshold of Theorem G at E = 2, C = 1. The norms do not depend on C. M_1 is the bound of Theorem G at C = 1 and E = 2; it is absolute and of order 10^7 with the referee's constants. Along a diagonal sequence rho_k -> 1 the fake has ~ s/2 marks.

## 2. The prime input of VIII.3 (item 1 of the task)

### 2.1 Every use of the prime set (Proved, by reading VIII.3 in full)

I went through VIII.3 (lines 1741-1879) and VIII.10 (1881-1903) line by line. The construction uses the following properties of its q - 1 primes (q = l), and no others.

| property | where it is used |
|---|---|
| exactly q - 1 = s **distinct** primes, one per column of I | definition of the profile and of w |
| all primes **split** | the model is a model of a Gaussian circle; a Gaussian integer of norm N exists (residue s_0) |
| p >= P | radii 2 exp(-w(v)/2) <= 2 P^{-h/2} via w(v) >= norm1(v) log P; b_0 = log L_0/log P < 0.51 in (VIII.9); the tail inequalities P/B >> H^2, B < P/q^{0.28}, 1/2 + 1.708/4 < 1.99/2; L_0 = o(P) in the pair instance of (R) |
| p <= 2P | 1/Delta <= 10 (2P)^{q/2}; P_+ <= 2P <= m^{1.995}; W <= 2(q-1) log 2P |
| W = (2d + o(1)) q log q, d = 1.99 | Delta = e^{-W/4}/10, the arc width, the weight (VIII.8) |
| p > 2m and p ≠ auxiliary primes | the point residues are units modulo every odd b <= 2m; the auxiliary primes ℓ_i <= q^{1.01} < P are disjoint from S |
| primes in a window of log-width sigma, randomly permuted | **only** the pair margin w(K_i - K_i') - W/2 >= 2 log P - o(1) |

No step uses any of the following:
* that the primes are consecutive;
* the specific interval of the original pigeonhole;
* the count of primes in the window beyond q - 1;
* any Gaussian argument of the parent's primes (the model replaces those angles);
* any arithmetic of the primes beyond the existence of one Gaussian integer of norm N.

The matrix steps never mention the primes: the uncertainty estimate, the gcd and norm tests on S, the Vandermonde step, the forced characters, the moment lifts and the digit tail. The same holds for VIII.10.

### 2.2 The pair margin in a widened bin (Proved; simulation Computed)

The coefficients. Take two rows i ≠ i' of the class A. For entries in {-1, 0, 1}, at most one of which is zero, |a - b| = 1 - ab. Hence c_j = |K_ij - K_i'j| - 1 = -K_ij K_i'j on I. These lie in {-1, 0, 1}, and their sum is

  -(KK^T)_{ii'} + K_i0 K_i'0 = 1 + 1 = 2.

So w(K_i - K_i') - W/2 = sum_j c_j x_{pi(j)}, where x is the vector of log-primes and pi is the assignment of primes to columns.

**Mean.** For pi uniform, E f = 2 x̄ >= 2 log P_b >= 2 log P.

**Bounded differences for permutations.** Let f(pi) = sum_j c_j x_{pi(j)}. A transposition of the values at two positions j, k changes f by (c_j - c_k)(x_{pi(k)} - x_{pi(j)}), which is at most 2 sigma in absolute value.

Reveal pi(1), pi(2), ... in order. Two continuations of the same history that differ in the value at step k are coupled by the transposition of those two values. Therefore the Doob martingale differences are at most 2 sigma. Azuma's inequality gives

  Pr(|f - E f| >= t) <= 2 exp(-t^2/(8 n sigma^2)), n = l - 1.

Take t = sigma sqrt(48 n log n). The failure probability is then 2 n^{-6}. A union bound over fewer than l^2/8 pairs of A leaves a permutation for which every pair satisfies

  |f - E f| <= sigma sqrt(48 l log l) = O(l^{-1/4} sqrt(log l)) = o(1).

Astra's coupling sentence is exactly this argument. The bound holds for every s-subset of the bin, because it uses only the values in the subset and their spread sigma.

**How much margin VIII.3 needs.** The pair instance of (R) and the two ordinary pair gaps need e^{(margin - 2 log P)/2} · P / L_0 -> infinity (fakecore §2.6). This holds for any deviation up to (1.96 - eps) log l. An o(1) deviation is therefore far more than enough.

**Remark (Proved, same estimate; not needed).** Bins of log-width sigma = sqrt(log l / l)/10 also work. Their deviation is at most 0.7 log l. They give K of order l^{1.49}/log l and J ~ 0.49 log l, hence log(2P)/J ~ 4.06 instead of 8.3.

**Simulation (Computed, floating point, `paley_checks.py`).** Values uniform in [0, sigma], and also a two-point profile, with a random permutation and the maximum over all pairs of A:

| q | max over pairs of abs(dev)/(sigma sqrt(n log n)) | max abs(dev) at sigma = q^{-3/4} |
|---|---|---|
| 227 | 0.46-0.62 (two-point: 0.80) | 0.27-0.37 (0.48) |
| 1019 | 0.48-0.50 (0.82) | 0.22-0.23 (0.38) |
| 2003 | 0.48-0.53 (0.87) | 0.20-0.22 (0.36) |
| 4003 | 0.50-0.55 (0.90) | 0.18-0.20 (0.33) |

The deviation scales as sigma sqrt(n log n) and decreases slowly at sigma = q^{-3/4}, as the proof predicts.

### 2.3 The later estimates of VIII.3 with sigma = l^{-3/4} (Proved; bin counts Computed)

* **"Weighted-bin error O(q^{1.5-d} log^{3/2} q)".** This is not a separate estimate. It is the same permutation deviation sigma sqrt(q log q), evaluated at the original sigma ~ q^{1-d} log q. With the widened bin it becomes O(l^{-1/4} sqrt(log l)). That is o(1) for every d, including the VIII.10 range d = 2 - xi.
* **Offsets.** They depend on W only, through Delta = e^{-W/4}/10, and on the auxiliary primes. W stays in [2s log P, 2s log 2P]. No change.
* **Digit expansion and tail.** They use P only through p >= P (P/B >> H^2, B < P/q^{0.28}, the exponent 1/2 + 1.708/4 < 0.995) and p <= 2P (1/Delta). No change.
* **Uncertainty estimate, gcd and norm tests, Vandermonde step, forced characters, (VIII.9), Haar measure.** These do not involve the primes. No change.
* **Side remark on VIII.11 (Proved, arithmetic).** The planar area test detects the first fake through (1 + sqrt(P) sigma)(1 + sqrt(P) eta). With sigma = l^{-3/4} and eta = O(l^{-0.243}) this product is about l^{0.997}. The test still excludes the first fake asymptotically, but only barely. This is about genuine angles and does not bear on the claim.

**PNT count (Proved).**
* By the prime number theorem modulo 4, [P, 2P] contains (1 + o(1)) P/(2 log P) split primes.
* There are ceil(log 2 · l^{3/4}) bins.
* Pigeonhole gives a bin with at least the average, (1 + o(1)) l^{1.99}/(2 · 1.99 · log l · log 2 · l^{3/4}) ≈ l^{1.24}/(2.76 log l).
* No short-interval prime theorem is needed.

**Bin counts (Computed, `bincount.py`, exact segmented sieve).** First bin [P, P e^{l^{-3/4}}]:

| l | P | split primes in first bin | PNT average bin | K/s |
|---|---|---|---|---|
| 10007 | 9.13e7 | 2,430 | 3,516 | 0.24 |
| 100003 | 8.91e9 | 34,640 | 49,078 | 0.35 |
| 1000003 | 8.71e11 | 500,825 | 712,646 | 0.50 |
| 10000019 | 8.51e13 | 7,461,769 | 10,636,348 | 0.75 |

The first bin holds about log 2 ≈ 0.69 of the average, as expected, because bins higher in [P, 2P] are longer. For l <= 10^7, even the average bin holds fewer than s primes, so no s-subset of a bin exists yet. The statement is asymptotic.

## 3. The genuine half (item 3)

**Theorem G as refereed (`ref9L_powers.md` §1).** Take a pool of K split primes in [P, 2P] and J = log(K/s), and assume J >= 400 E^2 log(E(E+1)^3). Then for s >= s_0(C, E), outside a fraction exp(-sJ/(250 E^2)) of the s-subsets, every arc of length C N^{1/4} on x^2 + y^2 = N contains at most 4(b-1) lattice points, where

  b - 1 = max(216 E^2, ceil(90000 E^5 log(2P)/J)),

for all exponents 1 <= e_p <= E chosen after the sample.

**Application. Proved.**
* The pool is the K split primes of I_l. They lie in [P, 2P], and Theorem G allows an arbitrary pool there.
* E = 2, and e_p = 2 for every p.
* J = log(K/s) = 0.24 log l - log log l + O(1) -> infinity, so J >= J_0(2) for l large. As computed, J >= 6383 needs log l >= about 26,640.
* log(2P)/J -> 1.99/0.24 ≈ 8.3, so b - 1 is bounded by an absolute constant.
* The exceptional fraction is exp(-sJ/1000) = exp(-(0.24/1000 + o(1)) l log l).
* Arc length: C N^{1/4} = C sqrt(R).
* Every lattice point of norm N = prod p^2 (all p split) has the form eta prod pi_p^{a_p} conj(pi_p)^{2-a_p}, which is exactly the form Theorem G covers. No Gaussian gcd is discarded.

**"One may take C = 1 once for all". Proved.** An arc of length C sqrt(R) is a union of ceil(C) arcs of length at most sqrt(R). Each of these carries at most M_1 = 4(b-1) points (at C = 1, E = 2), so the arc carries at most ceil(C) M_1. The good set of S is therefore chosen once, independently of C. The only condition is s >= s_0(1, 2), which is absorbed into l_0.

**Hypotheses.** Every hypothesis holds. The IX.0 squarefree condition K/s >= e^{40} is irrelevant here: the cap-two theorem is used, not the squarefree one.

## 4. Normalisation (item 4) (Proved; Computed)

* Column sums of K vanish, (K^T K)_{0v} = -1 for v ≠ 0, K_0v = chi(v), and K_i0 = -chi(i). Together these give sum_{i in A} K_iv = (-1 - chi(v))/2, which lies in {0, -1}.
* Column v has exactly one zero on all of F_q, at i = v. That zero lies in A exactly when v is in A.
* Hence each column v ≠ 0 has at least (q-5)/4 entries of each sign on A.
* The selection deletes (1 - rho)(q-1)/2 + O(q^{2/3}) < 0.005 q + o(q) rows of A. So both signs, that is both exponents 0 and 2, survive in every column of I.
* This argument is deterministic, so it holds for any selection of that size.
* Column 0 carries no prime.

So no prime of S is lost, and the normalised norm of the model is exactly N. **Computed** (`paley_checks.py`): q = 43, 67, 131, 227, 383. The minimum count of each sign on A per column is floor((q-1)/4), and the identity for sum_A K_iv holds for every column.

## 5. Meaning (item 5)

* **New relative to the record: yes.**
  * Before this result, fakes lived on narrow-bin parents. For those the record knows only VIII.1 (M <= 1 + omega + K_C ≈ s) and VIII.4 (weaker there). The record says they lie in the exceptional set of its random-parent theorem.
  * The record's dyadic theorem cannot reach bin parents. Their measure among dyadic s-subsets is about exp(-0.75 s log s). The dyadic exceptional fraction is exp(-c s log s) with c of order 10^{-2} to 10^{-3}, which is far larger.
  * Astra's device is to apply the random-support theorem with the bin itself as the pool. This is legitimate because Theorem G allows any pool in [P, 2P].
  * The novelty is moderate: it combines two refereed results with a one-line change of the bin width.
* **What it shows about typical parents.** On a common norm, the A_{2M} inputs of P2 hold with ~ s/2 marks while the genuine count is O(1). These inputs are monomial gaps, rational characters and point residues at odd prime powers <= 2M. So they carry no information about the genuine count, even on norms where the truth is bounded.
  * More precisely, the fake exists on **every** s-subset of the bin. The genuine count is O(1) on all but a fraction exp(-Omega(s log s)) of them.
  * A proof of M = o(s) that uses only these inputs fails on typical bin parents, not only on hypothetical bad ones.
  * It remains open whether fakes exist on typical **dyadic** parents. The construction needs a window of width o(sqrt(log l / l)) for its pair margins. A full dyadic interval gives deviations of order sqrt(l log l), far larger than 2 log P.
* **What it shows about worst-case parents: nothing.**
  * The fake is blind to the choice of subset. It exists equally on the exceptional subsets, where the genuine count is not controlled and might be large.
  * The result neither produces a large genuine cluster nor bounds one.
  * The exceptional parents, the squarefree comparable-prime worst case, and the growth rate for general R are all untouched.

## 6. IX.6, first paragraph (item 6) (Proved; Computed)

* **Parity is the singleton indicator of the owner.** The exponent of row i at column v is 1 + chi(v - i). It is odd exactly when chi(v - i) = 0, that is, when i = v.
  * The deleted column 0 has owner 0, which is not in A.
  * For a selected owner v, the column restricted to the selected marks is the indicator of {v}. It is nonconstant once there are at least two marks.
  * Subtracting a common minimum at a column adds a constant mod 2. With §4, no subtraction occurs anyway.
  * Parity equals [a >= 1] + [a >= 2] mod 2 (check 0+0, 1+0, 1+1).
  * Computed, q = 43 to 383: parity equals the owner indicator, and equals the sum of the threshold bits, at every (row, column).
* **Dimension count. Proved.**
  * A Boolean function g(L_1(x), ..., L_d(x)) of d affine forms on F_2^r has algebraic normal form of degree at most d. So every threshold bit, and hence every parity function, lies in RM(d, r).
  * RM(d, r) has dimension sum_{j <= d} C(r, j) <= (r+1)^d.
  * Restricted to the m marks (injectively labelled), this space contains the m singleton indicators. Hence m <= (r+1)^d.
  * With 2^r = O(s), r = log_2 s + O(1). This contradicts m ≍ s once s is large in terms of d.
  * Stronger form: any such labelling needs n = 2^r >= exp(c_d m^{1/d}). For P2, IX.5's bound n^{1 - 2^{1-2d}} is then trivial.
  * Computed: for d = 2 and s = 10^3, 10^4, 10^6, 10^9, dim RM = 67, 121, 232, 497, against m ~ s/2. For d = 4 and s = 10^3 the inequality is not yet violated (562 > 500). It is violated from s = 10^4 on.
* **Sidon property. Proved.**
  * The F_2-sum of the unary signatures of i and j has parity coordinates e_i + e_j. These are linear functions of the sum, and they determine {i, j}. They are nonzero for i ≠ j.
  * So D = |A+A| = 1 + C(m, 2) exactly.
  * Computed: distinct pair sums for every pair, for q = 43, 67, 131, 227, 383 (210 to 17,955 pairs).
  * Consequently IX.3's bound sqrt(D) ≈ m/sqrt(2) and IX.5's bound are both at least of order m on P2, so P2 evades both tests. Here "evades" means the tests give no contradiction. It does not mean P2 violates them.
* **The remaining sentences are correct.** The pool theorem's probability is over the prime subset, not over the assignment of a prescribed subset to columns. The primes have size polynomial in m (P_+ <= 2 l^{1.99}).

## 7. Standing cautions

* There is no exceptional set from a fixed-S Subspace theorem. The exceptional set comes from an explicit union bound over subsets of a fixed pool.
* No Gaussian gcd is discarded. Theorem G covers all points of norm N, and the model's normalised norm is exactly N.
* No drifting phase is replaced by zero. The model's common phase -(2/q) sum epsilon is carried exactly in VIII.3.
* No template-dependent threshold is treated as uniform. l_0 depends on rho and on the fixed E = 2, C = 1 threshold. It does not depend on C, which is handled by subdivision.
