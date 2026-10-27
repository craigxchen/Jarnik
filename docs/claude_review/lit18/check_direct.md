# lit18 / check_direct: adversarial check of `direct.md`

Report checked: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/lit18/direct.md`
(and its log `direct_log.txt`, its fetched cards in `gh/`, its arXiv corpus `tools/db.json`).

Method:
* 21 WebSearch queries of my own (listed in §4). After that the shared web-search budget for this turn ran out
  (two further queries were refused), so some secondary citations could not be re-verified by me. They are marked
  **UNVERIFIED (budget)** below and should not be cited on the strength of `direct.md` alone.
* My own regex sweep of the locally downloaded arXiv corpus (92,926 records, primary categories math.NT, CA, CO,
  AG, SP, MG, PR, AP, DS, GM, HO; first-seen dates 2024-03-05 to 2026-10-26; I confirmed these numbers).
* Every repo file and theorem that `direct.md` cites was opened and checked.
* Arithmetic checks in `chk/arith_check.py` (all pass; output quoted in §2).

---

## 0. Bottom line

**I agree with the verdict: nothing in `direct.md`, and nothing I found in addition, beats the barrier.** No
paper gives `o(log R/loglog R)`, or even `o(log R)`, for every `R`. Every lead fails. I found no surviving lead.

The bottom line is right, but four points in the reasoning are wrong or misleading, and the report missed some
items.

1. **Gabdullin is not "past `sqrt R`" in the sense that matters.**
   * His squares in `[N_0, N_0+N_0^gamma]` are lattice points near the **diagonal**, which is a direction of height
     `sqrt 2`. That is the same restricted position as Chan's near-axis result.
   * In that position, a general-position arc of length `R^lambda` corresponds to `gamma = (1+lambda)/2`.
   * His range `gamma < 0.618` therefore corresponds to general-position arcs of length `R^0.236`. That is already
     inside Cilleruelo-Córdoba's bounded range.
   * The program's `C sqrt R` would need `gamma = 3/4`.
   * So even a pointwise version of his divisor lemma would be useless. Its full-text request should drop from
     priority 1 to low.
2. **The Letendre/#886 "`O(1/delta)` in windows `n^{1/4-delta}`" statement is trivial, not "the CC analogue".** At
   the vertex, an elementary argument gives at most 2 divisors in any window of length `o(n^{1/4})`.
3. **Missed papers in the report's own area.**
   * T. H. Chan, arXiv:2507.07094 (2025), and Chan-Holmes-Liu-Villarreal, arXiv:2508.02818 (2025), on close
     factorizations and lattice points near the centre of `xy = N`.
   * Wang-Xiao, arXiv:2604.26497 (2026), and Huang-Levin-Xiao, arXiv:2406.18879, on families of `S`-unit
     equations. Wang-Xiao was found by the searcher but not analysed.

   All of them fail (§2.9-2.10).
4. **Hegyvári, arXiv:2504.13230, says in its abstract that the result is essentially Mei-Chu Chang's (JFA 2004)
   and that the author will withdraw the file.** `direct.md` does not mention this.

Smaller citation issues are in §1.

---

## 1. Citation verification

Legend:
* **OK**: authors, title, year, venue and main statement were confirmed by my own search snippet, or by an abstract
  in the local corpus.
* **OK(corpus)**: the abstract was read in `tools/db.json`.
* **PARTIAL**: some fields were confirmed.
* **UNVERIFIED (budget)**: I could not check it before the search budget ran out.
* **ISSUE**: an error or overstatement.

| # | Citation in `direct.md` | Status | What I confirmed / problem |
|---|---|---|---|
| R1 | K. Oganesyan, "Lattice points on small arcs", arXiv:2107.09991, withdrawn | PARTIAL | Title and number confirmed (W2). My Temur query (W12) returned the claim "for any α∈(1/2,1) the number of lattice points on an arc of length R^α is not uniformly bounded". The author line and the **withdrawal** come only from the searcher's snippets (Q1, Q4); I could not confirm them. See §2.8 for why the withdrawal is plausible. |
| R2 | M. R. Gabdullin, "Trigonometric polynomials with frequencies in the set of squares [and divisors in a short interval]", arXiv:2205.13611, v3 Dec 2023 | OK, **ISSUE in framing** | W1 confirms `‖f‖_4 ≪ ε^{-1/4}‖f‖_2` for frequencies `{n^2 : N ≤ n ≤ N+N^{γ_0-ε}}`, `γ_0=(√5-1)/2`, the reduction to divisors, and the Ruzsa link. The arXiv listing title is the short form; the long title is the PDF header. "CANT 2022 talk, 27 May 2022" is UNVERIFIED. The log's Q2 lists a different CANT 2022 talk ("Lattice points close to ovals, arcs, and helixes"), so the two may be conflated. Follow-up Gabdullin-Konyagin (cubes) is published in Math. Notes 115:3 (2024) 336-340 (W1). Framing problem: §2.1 below. |
| R3 | P. Letendre, "Divisors of an integer in a short interval", arXiv:2503.12146 (15 Mar 2025) | OK, **ISSUE** | Confirmed: title, author, date, `D_n(X,Y)`, Thm 1 of the form `τ(n)^{1-ξ(θ,η)}…`, Conjecture 1 (Ruzsa-type, window `n^{1/2-ε}`), `D_n(X,X) ≪ τ(n)/sqrt(Ω_2(n))` (W3, W5, corpus). The forum gloss (W6) is "for δ in (0,1/4) the divisor count in (n^{1/2}, n^{1/2}+n^{1/4−δ}) is ≪ 1/δ". The snippet calls this "above the n^{1/4} scale", but the window is below `n^{1/4}`. In that window the bound is trivial (≤ 2; §2.2). So `direct.md`'s "exactly the CC analogue" is wrong. |
| R4 | Erdős #886 (Ruzsa); Erdős-Rosenfeld (Acta Arith. 1997): ≤ 1+C^2 divisors in `[n^{1/2}, n^{1/2}+C n^{1/4}]`; infinitely many n with 4 divisors in `(n^{1/2}, n^{1/2}+16n^{1/4})` | OK | Confirmed via erdosproblems.com/latex/886 snippets (W4, W6). The Acta Arith. 1997 venue is UNVERIFIED (budget). The formal-conjectures PR #7218 ("rosenfeld_4 false, C=31, five divisors") is UNVERIFIED. Also, one numerical example cannot refute an "all sufficiently large `n`" statement unless it is an infinite family; treat it as noise. |
| R5 | T. H. Chan, "Factors of almost squares and lattice points on circles", IJNT 11 (2015) 1701-1708, arXiv:1406.2230; Thm 3: `n=N^2` ⇒ ≤ 10 points with `|b| < n^{1/4}(log n)^{1/7}` | OK | W7 confirms the arXiv paper and topic. The `gh/` card quotes Thm 3 verbatim from p. 2 and gives IJNT 11 (2015) **no. 5**, doi 10.1142/S1793042115400205. One search summary called 10 "a conjectured bound", but the card's verbatim quotation settles it as a theorem. "Five divisors for perfect squares, arXiv:1303.2069": W6 confirms the five-divisor statement; the arXiv number is UNVERIFIED. **Note:** Chan's 2025 paper (2507.07094) says it "correct[s] a mistake of a previous result" on three close factorizations. Which result is not stated, so the full text should be checked (§3). |
| R6 | F. Temur, Turkish J. Math. 46(7) (2022) 2867-2889, arXiv:2012.10784 | OK | W9: venue, pages and DOI 10.55730/1300-0098.3306. W12: Conjecture 3 (`N^η`, `η<1/2`, `≤ C_η`; known for `η<1/4`), and Theorem 5 via simultaneous Pell plus Schmidt. **W12 also confirms the kill:** "Schmidt's theorem is ineffective, which forces us to fix h_i beforehand, so the result gives no new results on Conjecture 1 for any values of N." |
| R7 | P. Kurlberg, S. Lester, "Poisson spacing statistics for lattice points on circles", Invent. Math. 240 (2025) 537-585, arXiv:2112.08522 | PARTIAL | W8: Invent. Math., 2025, pp. 537-585, DOI 10.1007/s00222-025-01324-1; density-one Poisson spacing. Volume 240 is not confirmed. The quote "at scales N^(−γ), ½<γ<log(π/2)/log 2, not i.i.d." is **not confirmed**. W11 places `log(π/2)/log 2 ≈ 0.6515` in an equidistribution estimate, and its summary reads the range the other way (i.i.d.-like), so the direction is unclear. The sentence "the project's paper adds that some subsequences have spacing measure degenerating to a point mass at 0" has no source. The kill (density one) does not depend on any of this. |
| R8 | A. Granville, I. Wigman, CMP (2017), arXiv:1612.07819 | OK | W10: CMP 355 (2017) 767-802, DOI 10.1007/s00220-017-2953-3; density-one subsequence. |
| R9 | Sartori arXiv:2005.04698; Klurman-Sartori arXiv:2104.13038; Benatar-Marinucci-Wigman JAM 141 (2020) 707-749; Cammarota-Klurman-Wigman CMP 2020 | OK | W18: Sartori's title and "generic integers … no quasi-correlations". W19: the Klurman-Sartori title is "On the behavior of nodal lines near the boundary for Laplace eigenfunctions on the square" (MPIM 2021-32), whose main new ingredient is lattice points in small arcs. "Lemma 1.5, density-one separation" is UNVERIFIED (budget). W21: BMW JAM 141(2) 707-749, DOI 10.1007/s11854-020-0114-7. W20: CKW CMP 376(2) (2020) 1261-1310, DOI as stated, arXiv 1903.10602. |
| R10 | Kurlberg-Wigman, Math. Ann., arXiv:1501.01995 | UNVERIFIED (budget) | Irrelevant to the verdict. |
| R11 | Bombieri-Bourgain, IMRN 2015 no. 11, 3343-3407, doi 10.1093/imrn/rnu005 | OK | W17. The content (sum-product, Diophantine approximation, elliptic curves) is not confirmed by my snippet, only the toy-model framing. |
| R12 | Cushman-Demeter-Wu arXiv:2608.12316; Li-Yung arXiv:2609.33949 | OK(corpus) | The abstracts match `direct.md`, and `local_smoothing/chains.md` row E already records CDW. |
| R13 | Bourgain-Rudnick arXiv:1012.3843 (AHP 2011), 1105.0018, 1308.4247; "Lemma 2.1: O(log m) points on arcs ≤ (√m)^{1/2}"; "B_λ ≪ log λ individually" | PARTIAL | W13 confirms the titles: 1012.3843 "On the Geometry of the Nodal Lines of Eigenfunctions of the Two-Dimensional Torus", 1105.0018 "Restriction of toral eigenfunctions to hypersurfaces and nodal sets", 1308.4247 "Nodal intersections and L^p restriction theorems on the torus". It also confirms the definition of `B_λ` and "known for almost all λ". The O(log m) arc lemma is confirmed (W2, W14). The AHP venue and the separation count `O(X^{1-ε/3})` are UNVERIFIED (budget). |
| R14 | Squares in AP: González-Jiménez 2602.03251; González-Jiménez-Tho 2602.01380; Hegyvári 2504.13230; Elsholtz-Wurzinger QJM 75(4) (2024) 1243-1254; González-Jiménez-Xarles 1301.5122 | OK / **ISSUE** | The two González-Jiménez items are OK(corpus). Elsholtz-Wurzinger is OK (W16, DOI 10.1093/qmath/haae044). **Hegyvári:** the corpus abstract ends "This result essentially is proved -- with the same tools -- by Mei-Chu Chang. See in J. Funct. Anal. 207 (2004) … So the author will withdraw this ArXiv file". `direct.md` omits this, and its "uses BZ as input" is not in the abstract. González-Jiménez-Xarles is UNVERIFIED (budget). |
| R15 | V. Ziegler, "On the number of solutions to S-unit equations", arXiv:2609.15487 | OK(corpus), **ISSUE in arithmetic** | Abstract as stated. `direct.md` §2.7 writes `(s-1)/2`, with `s` the number of split primes. Ziegler's count is `(|S|-1)/2`, where `S` is a set of places of `Q(i)`; each split prime contributes 2 places, so `(|S|-1)/2 ≥ s`. The kill only gets stronger (§2.6). |
| R16 | Ferreira-Ribas arXiv:2602.23597; Calegari-Dimitrov-Tang arXiv:2510.04156 | OK(corpus) | Both abstracts as described. CDT also gives irrationality measures for `L(2,χ_{-3})`, the 2-adic `ζ(5)`, and a proof of the transcendence of π. Nothing about many logarithms. |
| R17 | Bloom-Lichtman 2312.12890; "Lattice points on a curve via ℓ^2 decoupling" 2205.08206; Bhattacharjee 2607.24048; BCN | PARTIAL | Bloom-Lichtman and Bhattacharjee are OK(corpus). 2205.08206 is UNVERIFIED (budget; no authors given). |
| R18 | Guo-Ilyin 2309.10971; Waxman-Yesha 2305.02060; Lester-Waxman-Yesha 2609.23214; Burq-Germain-Sorella-Zhu 2507.16599; Zhang-Zhu 2606.08650; Korsky 2606.01780; DGH Ann. of Math. 194 (2021) 237-298; Huang-Liu-Rudnick 1903.04005; Greenfeld-Iliopoulou-Peluse 2401.10821 | PARTIAL | The first six are OK(corpus); Zhang-Zhu are Cheng Zhang and Zhifei Zhu. DGH matches my memory (not searched). Huang-Liu-Rudnick and Greenfeld-Iliopoulou-Peluse are UNVERIFIED (budget). All are irrelevant to the verdict. |
| — | Garaev-Garayev-Konyagin arXiv:2602.08807 | OK | W1 and the corpus. The abstract's results are for cubes and quartics; none is about squares. |
| — | Repo citations: `round10/analytic.md` Thm B; `local_smoothing/chains.md` rows E, F, G, H; `round8/README.md` LW row (line 365) and `lw.md` Thms A-B; `round8/stability.md` Prop 6.1; `round16/README.md` §7.4; `threshold/README.md` | OK | All exist and say what `direct.md` attributes to them. One detail: the round-16 "algebraic independence of ~n logarithms" remark is labelled "heuristic, unverified, from memory" by the judge. `direct.md` cites it as a remark, which is fine. |

**Missed items in the report's area (found in the corpus):**
* **T. H. Chan, "Numbers with three close factorizations and central lattice points on hyperbolas",
  arXiv:2507.07094 (9 Jul 2025).** Abstract: it "correct[s] a mistake of a previous result"; relates three close
  factorizations to lattice points close to the centre `(sqrt N, sqrt N)` of `xy = N`; proves optimal lower bounds
  for the `L^1` distance of these points to the centre; and gives examples from polynomials and Pell equations.
  `direct_log.txt` Q29 reported "nothing after 1406.2230" for Chan.
* **T. H. Chan, L. Holmes, M. Liu, J. Villarreal, "Numbers with Four Close Factorizations", arXiv:2508.02818
  (5 Aug 2025).** Abstract: `n = AB = (A+a_i)(B-b_i)` for `i = 1,2,3`, with `B ≤ A` and `1 ≤ a_i, b_i ≤ C`,
  forces `A ≤ 0.04742… C^3 + O(C)`. The proof uses generalized Pell equations.
* **J. T.-Y. Wang, Z. Xiao, "Families of Unit Equations and Exponential Diophantine Problems via Integral Points",
  arXiv:2604.26497 (2026).** It gives degeneracy results for one-parameter families of unit equations, via Ru-Vojta
  and a generalized Huang-Levin-Xiao inequality. It is in the searcher's log only and was not analysed.
* **K. Huang, A. Levin, Z. Xiao, "A New Diophantine Approximation Inequality on Surfaces and Its Applications",
  arXiv:2406.18879.** GCD inequalities and families of `S`-unit equations.

---

## 2. Kill attempts, lead by lead

Notation as in `direct.md`:
* `N = R^2` is the norm, and `z_1..z_M` are in `Z[i]` with `|z_j|^2 = N`.
* The cluster lies on an arc of length `L = C sqrt R`.
* `s` is the number of split primes of `N`, with `s ≤ ω(N) ≤ (1+o(1)) log N/loglog N = (2+o(1)) log R/loglog R`.

### 2.1 Gabdullin (R2): fails twice over

**(a) Energy vacuity (`direct.md`'s argument).** I checked it.
* Take `A = {a_j} ∪ {b_j}` and unit weights. Then `‖f‖_4^4 ≥ (2M)^2`, from the ordered pairs representing the one
  integer `N`, and `‖f‖_2^4 = |A|^2 ≤ 4M^2`.
* So `‖f‖_4 ≤ c‖f‖_2` is vacuous for every `M`, with or without weights, whenever `c^4 ≥ 1`.
* Capped: `README.md` §3, energy class; `local_smoothing/chains.md` row F.

**(b) Wrong position, so the exponent is too small even pointwise.** This is new, and it removes the report's "still
worth checking" item.
* Squares in `[N_0, N_0+H]` that solve `a^2+b^2 = c^2+d^2` are pairs of lattice points in an `H x H` box at
  `(N_0, N_0)`. That is an arc near the diagonal direction `(1,1)`, which has height `sqrt 2`.
* In that position the sagitta argument applies (`direct.md` §2.3 with `H_dir = sqrt 2`). With `u = a+b` and
  `|a-b| ≤ H`, at most `2(H^2/(sqrt2 N_0) + 2)` points lie on any one circle.
  * The check (`chk/arith_check.py`, brute force over the boxes) is consistent: the maximum is 2-4 points per
    circle, against a bound of 7-27, for `N_0 = 10^4..2·10^5` and `γ = 0.54..0.65`.
  * So in this position the `sqrt R` scale (`γ = 1/2`) is trivial pointwise, and "past `sqrt R`" is the wrong
    yardstick.
* **The correct comparison uses the Cilleruelo-Granville rotation** (`gh/` CG 2007 card, p. 11).
  * A general-position cluster on an arc of length `R^lambda` of `|z|^2 = N` becomes, after multiplying by
    `\bar z_1` and by `1+i`, a near-diagonal cluster.
  * The new circle has radius `~ R^2`, so `N_0 ~ R^2/sqrt2`, and the box side is `H ~ R · R^lambda`.
  * Hence `gamma = (1+lambda)/2`, i.e. `lambda = 2gamma - 1`:

  | γ (Gabdullin normalisation) | general-position arc exponent λ = 2γ−1 | meaning |
  |---|---|---|
  | 0.5 | 0.0 | Sidon threshold |
  | 0.618 (Gabdullin, averaged) | **0.236** | inside Cilleruelo-Córdoba's proved bounded range (any λ < 1/2) |
  | 0.667 | 0.333 | Jarník scale |
  | **0.75** | **0.5** | the program's `C sqrt R` |

* So even a **pointwise** form of Gabdullin's divisor lemma at his full range gives nothing in general position. The
  program's arcs need `γ = 3/4`, and his averaged method stops at 0.618.
* Equivalently, in CG's near-axis normalisation `(5.1)`, with window `n^alpha` and `alpha = gamma/2`:
  * pointwise results reach `alpha = 1/4 + o(1)` (sagitta; Chan's `(log n)^{1/7}`);
  * averaged results reach `alpha ≈ 0.309` (Gabdullin);
  * the target is `alpha = 3/8`.
* The sentence in `direct.md` §0, "averaged (fake-compatible) statements extend past `sqrt R`, while the pointwise
  statement is stuck below it", should be restated in this normalisation.
* **Verdict: fails.** It is capped by the energy class, and it falls short in exponent even if made pointwise.

### 2.2 Letendre and Ruzsa #886 (R3, R4): fails; the report's reading is wrong

* **Elementary vertex bound.** If `d | n` and `d ∈ [sqrt n, sqrt n + Y]`, then `s = d + n/d` is an integer in
  `[2 sqrt n, 2 sqrt n + Y^2/sqrt n]`, and each `s` gives one such `d`. So the count is at most
  `Y^2/sqrt n + 2`.
  * Checked on 3,987 smooth `n < 10^30`: no violation. The maximum count in a window `n^{0.20}` was 1.
* For `Y = n^{1/4-δ}` this gives at most 2. The forum's `O(1/δ)` is therefore weaker than trivial, and it is not
  "exactly the CC analogue (`M(δ) ~ 1/δ`)".
* The chopping arithmetic in `direct.md` §2.2 (giving `O(log n)`) is moot: Erdős-Rosenfeld's `1 + C^2` is the same
  elementary argument.
* Above `n^{1/4}`, Letendre gives `τ(n)^{1-ξ}`, which is `exp(c log n/loglog n)`. That is useless, as `direct.md`
  says.
* The vertex `θ = 1/2` is the height-1 position (the analogue of Chan's case), not general position.
* **Verdict: fails.** No input. The rational problem in general position is open (#886 open; Letendre's Conjecture 1
  is unproved).

### 2.3 Chan 2015, Erdős-Rosenfeld, and Chan 2025 (R4, R5 and the missed papers): restricted position only

* `direct.md` §2.3's transfer computation is correct, and I re-derived it.
  * `w_j = z_j \bar z_1` has norm `n' = N^2 = R^4` and satisfies `|Im w_j| ≤ C R^{3/2} = C n'^{3/8}`, while
    Chan's window is `n'^{1/4}(log n')^{1/7}`.
  * The shortfall is the factor `R^{1/2}/(4 log R)^{1/7}`. In CG's normalisation that is `alpha = 1/4 + o(1)`
    against the target `3/8`.
* Dividing the `w_j` by their common factor undoes the rotation exactly: the common factor is
  `\bar z_1 · gcd(z_j)`, and dividing by it returns the `z_j`. So there is no gcd-reduction loophole.
* **Chan 2507.07094** (central lattice points on `xy = N`): the vertex of the hyperbola, which is the height-1
  position. The output is lower bounds for the distance of three points to the centre: a fixed number of points,
  no growth. Same class as R5.
* **Chan-Holmes-Liu-Villarreal 2508.02818:** four points of `xy = n` in a `C x C` box force `A ≲ 0.0474 C^3`, i.e.
  `C ≳ 2.76 A^{1/3}`.
  * The radius of curvature at the vertex is `~ A`, so this is the **Jarník/`R^{1/3}` scale**, the analogue of
    Cilleruelo-Granville's "Close lattice points on circles" (CJM 2009). It is not the `sqrt R` scale.
  * It is a fixed count (four), with no growth content.
* **Verdict: fails.** These are low-height-position or Jarník-scale results. The program already treats near-axis
  clusters as easy (`direct.md` §2.3; CG 2007 p. 11: `(5.1)` "simple for alpha ≤ 1/4").

### 2.4 Density-one results (R7, R8, R9, R13 separation): fail

* These are statistical or almost-all statements. They are capped by `round10/analytic.md` Thm B and by
  `local_smoothing/chains.md` row G.
* Arithmetic is fine: `r(N) ≤ exp((log 2+o(1)) log N/loglog N) = R^{o(1)}`, so the natural angular scale
  `1/r(N)` is `R^{1/2-o(1)}` times coarser than the cluster's `C R^{-1/2}`.
* The uncertain Kurlberg-Lester quotation (§1) does not affect this.

### 2.5 Cushman-Demeter-Wu (R12): fails

* `M^{3-ε} ≤ c_ε((9/8)C^3 + O(C)) sqrt R`, so `M ≪ R^{1/6+ε}`. It holds for every finite set on a strictly convex
  curve, so for fakes too.
* Matches `chains.md` row E exactly.

### 2.6 Ziegler (R15): fails, with corrected arithmetic

* Using `u_jk = z_j/z_k`, an `S`-unit of `Q(i)` with `|u| = 1`, the cluster gives `u + u^{-1} = 2 - h_jk/N`, where
  `h_jk = |z_j - z_k|^2 ≤ C^2 R`.
  * The right-hand side changes with `h_jk`, so Ziegler bounds only the multiplicity of each chord length, not `M`.
  * The height threshold depends on `S`; that is an `S`-dependent threshold, excluded by the constraints.
* Even if it applied, `|S| ≥ 2s+1` gives `(|S|-1)/2 ≥ s`, which can be as large as `(2+o(1)) log R/loglog R`,
  three times the proved `(2/3+ε)`.
* `direct.md`'s "`(s-1)/2 ≤ (1+o(1)) log R/loglog R`" understates this, but the kill only gets stronger.

### 2.7 Ferreira-Ribas, Calegari-Dimitrov-Tang (R16): fail

* Both are single-form measures. A Baker-type `|Arg(z_j/z_k)| ≥ N^{-c'}` is weaker than the trivial Liouville bound
  `|sin(θ_j - θ_k)| ≥ 1/N`, since `Im(z_j \bar z_k)` is a nonzero integer.
* Both bounds are far below the cluster scale `C N^{-1/4}`.
* Capped: `round8/README.md` line 365 and `lw.md` Thms A-B, `stability.md` Prop 6.1.
* The CDT abstract claims irrationality and linear-independence measures and reproves the classical two-variable
  `S`-unit and Thue-Mahler results. Those are fixed-`S`, with counts `c^{|S|}`. It claims nothing about many
  logarithms.

### 2.8 Oganesyan (R1): irrelevant to the upper bound; the withdrawal is consistent with the repo

* Arcs `R^alpha` with `alpha > 1/2` give no upper-bound content at `sqrt R`.
* New observation: `growth/README.md` (conditional theorem) shows that Vojta's Main Conjecture for one fixed blow-up
  gives bounded counts on arcs `C R^alpha` for every `alpha < 1`.
  * So a correct Oganesyan theorem would have **refuted that instance of Vojta's conjecture**.
  * That is a strong reason to expect the claim to have been wrong. It is consistent with the withdrawal, which I
    could not independently confirm.
* No action needed.

### 2.9 Temur (R6): fails

* My snippet (W12) confirms that Schmidt's theorem forces the parameters `h_i` to be fixed beforehand and "gives no
  new results … for any values of N".
* This is a template-dependent, ineffective threshold, which the constraints exclude.

### 2.10 Families of unit equations and GCD inequalities (Wang-Xiao, Huang-Levin-Xiao; missed or not analysed): fail

* The natural formulation is the family `x + x^{-1} = 2 - h/N`, with `h` varying. These papers give degeneracy, or
  finiteness outside exceptional sets, for **fixed `S`**, via Ru-Vojta or Second-Main-Theorem inequalities.
* Here `S` = the primes of `N`, which vary with `R`. Their exceptional sets depend on `S`.
* This is exactly the case `threshold/README.md` records as failing ("Ru-Vojta and Subspace routes fail … with `S` =
  primes of `N` (exceptional sets depend on `S`)"), and it is excluded by the constraint on fixed-`S` Subspace
  exceptional sets.
* Even where they are effective, the solution counts are `c^{|S|}`-type, which is useless at
  `s ~ log R/loglog R`.

### 2.11 Squares in AP, DGH (R14): fail

* Only the direction arcs ⇒ squares is known (CG 2007; `chains.md` row F).
* DGH counts are `c(g,d)^{1+rank}`, and in any transfer the rank grows with `s`: the forbidden `c^s` regime.
* Hegyvári is self-withdrawn and is in any case a sumset statement.

---

## 3. Papers to request in full text (ranked; none is expected to beat the barrier)

1. **Chan, arXiv:2507.07094 (2025).** This is a record check, not a lead. It says it corrects "a mistake of a
   previous result" on three close factorizations. The full text is needed to see whether the corrected result is
   one the program quotes: the IJNT 2015 Theorem 2 (18 divisors), its Theorem 3 (≤ 10 near-axis points), or
   arXiv:1303.2069 (five divisors). Theorem 3 is quoted in `direct.md` §0 and §2.3.
2. **Letendre, arXiv:2503.12146.** Only to confirm the exact form of Thm 1 and Prop 1 and the forum's
   specialisation. §2.2 shows that the window quoted for the specialisation is in the trivial range.
3. **Gabdullin, arXiv:2205.13611.** Downgraded from rank 1 in `direct.md`. Per §2.1(b), a pointwise version of his
   lemma at `γ < 0.618` corresponds to general-position arcs `R^{0.236}`, so it cannot matter for growth. It is
   worth reading only for the record.
4. **Ziegler, arXiv:2609.15487; Wang-Xiao, arXiv:2604.26497.** Only to document how the thresholds and exceptional
   sets depend on `S`, for the barrier map.
5. **Kurlberg-Lester, Invent. Math. 2025.** Only to correct the quotation about i.i.d. behaviour at scales
   `N^{-γ}` in `direct.md` R7.

`direct.md`'s entries 6-7 (Temur, Bombieri-Bourgain) are not needed; Temur's limitation is confirmed by snippet (W12).

---

## 4. Search log

**WebSearch (mine).** [std] means standard mode, [ext] means extended.
* W1 [std] Gabdullin "Trigonometric polynomials with frequencies in the set of squares and divisors in a short
  interval"
* W2 [std] Oganesyan "Lattice points on small arcs" arXiv 2107.09991
* W3 [std] Letendre "Divisors of an integer in a short interval" arXiv 2503.12146
* W4 [ext] erdosproblems 886 Letendre theta = 1/2 divisors n^{1/2} n^{1/4} uniform bound forum
* W5 [ext] Letendre "Divisors of an integer in a short interval" abstract D_n(X,Y) conjecture … Proposition 1
* W6 [std] "886" erdos problem Ruzsa divisors "n^{1/2-\varepsilon}" Letendre Chojecki comment
* W7 [std] Tsz Ho Chan "Factors of almost squares and lattice points on circles" IJNT 2015
* W8 [std] Kurlberg Lester "Poisson spacing statistics for lattice points on circles" Inventiones
* W9 [std] Temur "Discrete fractional integrals, lattice points on short arcs, and diophantine approximation"
  Turkish J. Math.
* W10 [std] Granville Wigman "Planck-scale mass equidistribution of toral Laplace eigenfunctions" CMP
* W11 [std] Kurlberg Lester "do not behave like" i.i.d. scales log(π/2)/log 2
* W12 [std] Temur 2012.10784 short arcs conics simultaneous Pell Schmidt subspace conjecture N^eta C_eta
* W13 [std] Bourgain Rudnick B_λ log λ lattice points arcs λ^{1/2} "conjecture is known for almost all"
* W14 [ext] lattice points arc circle sqrt R "log R/log log R" Cilleruelo. Nothing at `log R/loglog R`; confirms
  `O(log m)` (Bourgain-Rudnick) and `≪ c log λ` for arcs `c sqrt λ`.
* W15 [std] Cilleruelo Granville "Lattice points on circles, squares in arithmetic progressions and sumsets of
  squares" conjecture R^{1-ε}
* W16 [std] Elsholtz Wurzinger "Sumsets in the set of squares" QJM 2024
* W17 [std] Bombieri Bourgain "A problem on sums of two squares" IMRN 2015 3343
* W18 [std] Sartori "spectral quasi-correlations" 2005.04698
* W19 [std] Klurman Sartori 2104.13038 "lattice points in small arcs"
* W20 [std] Cammarota Klurman Wigman "spectral semi-correlations" CMP doi
* W21 [std] Benatar Marinucci Wigman JAM 141
* Refused because the budget was exhausted: Klurman-Sartori lemma [ext]; Kurlberg-Wigman Math. Ann. [std].
* WebFetch of erdosproblems.com and arxiv.org failed (DNS), as expected.

**Local corpus regexes** (`chk/q.py` over `tools/db.json`; hit counts in brackets):
* `lattice point` × `circle|arc` [19]: only CDW is relevant (energy).
* `concyclic|cocircular` [9]: none relevant.
* `gaussian integer` × `angle|argument|arc|sector` [8]: none.
* `S-unit|unit equation` × `number of solutions|uniform|at most` [6]: Ziegler; Miyazaki-Pink; Le-Miyazaki
  (`a^x+b^y=c^z`; fixed bases, irrelevant).
* `subspace theorem` [18]: Adiceam-Shirandami, "Probabilistic effectivity" (2411.01247), is a density statement and
  so capped. The rest are irrelevant.
* `linear forms in logarithms` × `independen|many|measure` [8]: CDT; nothing on many logarithms.
* `Vojta` [27]: Wang-Xiao 2604.26497, Huang-Levin-Xiao 2406.18879, Yu-Yuan-Zhou 2602.01820 (explicit uniform
  Mordell). All are fixed-`S` or `c^{rank}` (§2.10, §2.11).
* `sums of two squares` × gap/short/close/arc/angle [2]: Guo-Ilyin, Iyer 2505.23428 (irrelevant).
* `simultaneous Pell|Pell equations` [17]: **Chan 2507.07094 and Chan-Holmes-Liu-Villarreal 2508.02818 (missed by
  `direct.md`)**.
* `Rosenfeld|close factori|almost square` [13]: the same two Chan papers.
* `Ruzsa` × `divisor` [3]: none relevant.
* `squares` × `arithmetic progression` [17]: nothing beyond `direct.md`.
* `Bombieri-Pila|determinant method|dimension growth` [41]: nothing for circles at fixed degree 2.
* `arithmetic random wave|toral eigenfunction|nodal…torus` [8]: Zhang-Zhu, Wu-Xiao 2610.23836 (restriction lower
  bounds); none pointwise.
* `algebraic independence` × `logarithm` [7]: function-field or elliptic variants only. Nothing gives an
  algebraic-independence measure for ~n classical logarithms, which the round-16 remark would need.
* `Gaussian primes` [11]; `gcd` × `S-unit|unit` [24]: none relevant beyond §2.10.

**Arithmetic checks:** `chk/arith_check.py`. The output: no violations of the vertex divisor bound; the maximum
number of divisors in `[sqrt n, sqrt n + n^{0.2}]` is 1; the near-diagonal maximum per circle is 2-4, below the
sagitta bound; and the `γ → λ = 2γ−1` table.
