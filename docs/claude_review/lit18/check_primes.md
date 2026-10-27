# check_primes: adversarial check of `lit18/primes.md` (genuine-prime statistics: H^gen, L1, L7)

**Inputs.**
* The report under review: `.../scratchpad/lit18/primes.md`, with its arithmetic in `primes_checks/arith.py` and
  `arith_out.txt`.
* The barrier map: `README.md` §3, `growth/README.md`, `threshold/README.md`, `round10/README.md` §7.1-7.2,
  `round16/README.md` §7.2.
* The repo files the report cites: `round9/README.md` §3-4, `round7/p3lower.md` §4 (H, Thm C, Lemmas C.1-C.3),
  `round15/rigSLtheory.md` §4, `round10/analytic.md` §3-5 and `round10/referee_analytic.md`,
  `local_smoothing/chains.md` row E, `two_thirds/proof.md`.

**Method and a limitation you need to know.**
* **WebSearch was unavailable.** The first query returned "this turn's web search budget is used up (limit: 200
  WebSearch calls per turn, shared by every agent in it)". I did not work around the limit. **No WebSearch query was
  run by me.** If you want the pre-2024 citations re-verified by search, send a follow-up message after this one.
* **Citation verification therefore uses:**
  * the local arXiv RSS mirror built by the sibling agent (`.../lit18/tools/db.json`: 92,926 records, titles and
    abstracts, first announcements and replacements from 2024-03-05 to 2026-10-27), read-only, through my own scripts
    in `.../lit18/check_primes_work/` (run with `python3 -I`);
  * for older papers, the independent snippets recorded by the sibling checker in `.../lit18/check_direct.md`
    (W-numbered entries), where they overlap.
  * Everything else is marked **UNVERIFIED (budget)**.
* **My arithmetic:** `.../lit18/check_primes_work/arith_check.py`, with output in `arith_check_out.txt` (same
  folder).

---------------------------------------------------------------------------------------------------------------

## 0. Bottom line

**The report's verdict survives: nothing in the genuine-prime literature plausibly beats the barrier.** No result
found, by the report or by my extra mirror searches, gives H^gen, L1, L7, or any `o(log R/loglog R)` statement.

The report singles out three things as closest. All three fail.

| lead | report's kill | my verdict |
|---|---|---|
| Croot-Mao-Pohoata-Sheffer-Yip, combinatorial large sieve (2606.17487) | same mechanism as `two_thirds`; fixed-`A` residue data, so capped | **Fails.** The report's prime-range argument rests on an unconfirmed (and probably conflated) "(log N)^3" claim. My kill does not depend on the range: counting, pigeonhole or entropy can force collisions only at moduli below `M`, i.e. at most `j < 1/(1-eta)` window primes jointly. Theorem C needs `k ~ 2 log M/loglog M -> infinity`. `D_A` satisfies every counting-forced statement automatically (§2.1). |
| Ford (2408.03803) and Gorodetsky (2608.14190), Kubilius/Poisson models | right shape, impossible strength | **Fails** (confirmed). The `k`-fold joint modulus exceeds the size of the pair family by a power `(1-eta)(1-eps/8) log M/loglog M -> infinity` for every fixed `eta`. The family is finite and adversarial, so it has no level of distribution beyond its own size (§2.2). The report's asymptotic formula for this ratio is off by a factor of 2 (harmless). |
| Merikoski, Gaussian primes in sparse sets (2302.11331), "uniform over arbitrary `B`" | `\|B\| >= X^0.4` vs selection `<= X^(1/h_0)`, "always too sparse" | **Fails, but the report's reason is incomplete.** In the low-height regime (cluster primes `<= M^2.5`, which includes the uncapped `P_+(N) < M^(2-delta)` window) the size condition *can* be met. The kill must instead be: the theorem is a count over all primes with `b in B` (main term `>= X^0.9/log X >> M`), it says nothing about angle relations, and for arbitrary `B` it gives only a lower bound (§2.3). |

**Citation and attribution problems** (details in §1):
1. "Its entropic version only uses residues modulo products of primes up to `(log N)^3`" (summary and §0 table) is
   not in the 2606.17487 abstract. `(log N)^3` appears in the abstract of a *different* paper, 2607.15311 (the CMPY
   inverse theorem, Selberg-sieve variant). This is probably a conflation. §2.6 of the report hedges ("`<= log N` or
   `<= (log N)^3`"), but the summary states it as fact.
2. **Mis-attributed cap.** The report caps Kuś (2605.19105, Halász over *all* ideals), Lester-Waxman-Yesha
   (2609.23214, *all* lattice points in sectors), Järviniemi-Teräväinen (almost primes, Matomäki-Radziwiłł method)
   and Green-Sawhney (Type I/II sums) under `round10/analytic.md` Thm B. But Thm B covers **prime-supported** inputs
   only. `round10/README.md` §7.2 and `referee_analytic.md` "Slip 3" list composite-supported inputs (the all-ideal
   large sieve, Type I/II sums) as **unclassified**, because they need a multiplicative graft, which is unproved.
   * The items still fail, but for other reasons: averages over `N` and over sector position (`analytic.md` §4.2),
     "almost all" exceptional sets, the existence direction, or irrelevance.
   * The barrier-map entry for this class stays "unclassified", not "capped".
3. **Arithmetic slip.** The log-ratio `log Q_k / log binom(M,2)` is stated as `~ (1-eta) log M/(2 loglog M)`. It is
   `(1-eta)k/2 ~ (1-eta)(1-eps/8) log M/loglog M`, twice as large. The report's own table (2.56 at `M = 10^6`)
   matches the corrected formula, not the stated one (1.32). The conclusion only gets stronger.
4. **Overstated or unconfirmed content:**
   * Kurlberg-Lester's "not i.i.d. at scales `N^(-gamma)`, `1/2 < gamma < log(pi/2)/log 2`": the sibling checker's
     snippet reads the range the other way.
   * Sartori's "generic = i.i.d.-angle model": not in any independently seen snippet.
   * Merikoski "robust to adversarial selection": for arbitrary `B` the abstract gives only a lower bound, with a
     weaker `delta` range.
   * Thm B's hypothesis is quoted without `A >= 4`.
   * Gou's GCD theorem is quoted without `k >= 3`.
   * Brumley-Thorner-Zaman's content is described as "Chebotarev least-prime bounds" (unverified).
   * The search log says 46 queries but lists 48.
5. **Unverified (budget):** about 25 pre-2024 citations (§1.3). None of them carries the verdict.

**Extra searches** (mirror only): 9 further items in the area (§1.4). None survives. The relevant ones:
* Fan 2508.06005 (Hardy-Ramanujan for sifted sets: upper tails, the fake-satisfied direction);
* Kai 2306.16983 (linear patterns of prime elements in number fields: existence);
* Gao-Li-Qi 2608.05934 (hybrid Weyl subconvexity for Hecke characters: robust, so Thm B);
* Bosio-Riblet-Tarr 2607.02728 (`B_2[g]` subsets of squares of size `n^(2g/(2g+1))`). It shows that the
  combinatorial-large-sieve saving is far from the truth even in its home problem. It is a generic mechanism, not a
  rate-changer.

---------------------------------------------------------------------------------------------------------------

## 1. Citation verification

Legend:
* **OK(mirror)**: authors, title, date and main statement confirmed by me from the RSS-mirror abstract.
* **OK(sib)**: confirmed by an independent snippet recorded in `check_direct.md`.
* **PARTIAL**: some fields confirmed.
* **UNVERIFIED (budget)**: no independent check was possible.
* **ISSUE**: a discrepancy with what the source says.

### 1.1 Items the verdict actually leans on

| # | citation as given in `primes.md` | status | notes |
|---|---|---|---|
| 33 | Croot, Mao, Pohoata, Sheffer, Yip, "A combinatorial large sieve for Sidon sets, distances, and norm forms", arXiv 2606.17487 (v2 25 Jun 2026) | **OK(mirror) / ISSUE** | Mirror: first 2026-06-17, v2 2026-06-25, math.NT, authors exact. Abstract confirms: the quoted mechanism sentence; Sidon subsets of `{1^2..N^2}` and repeated-distance-free subsets of `[N]^2` have size `<= N exp(-c log N/loglog N)`; the entropic version for `B_2[g]` in squares, norm forms over arbitrary number fields, `B_3[g]` in cubes and `B_4[g]` in fourth powers. **Not in the abstract:** which primes are used. The phrase "primes up to `(log N)^3`" occurs in 2607.15311, not here. The report's §0 claim "its entropic version only uses residues modulo products of primes up to `(log N)^3`" is therefore unsupported. The "Shearer's inequality" detail comes from the report's web snippet (search 42), which I cannot re-check. |
| 34 | Croot, Mao, Pohoata, Yip, "A sharp inverse theorem for the quadratic large sieve", arXiv 2607.15311 | OK(mirror) | As stated. Also: restrictions for primes up to `(log N)^3` give `exp(c sqrt(log N/loglog N))` via a Selberg-sieve-type argument. This is the source of the `(log N)^3`. |
| 34b | Croot, Mao, Yip, "An inverse theorem on sets with rich additive structure modulo primes", arXiv 2510.08862 | OK(mirror) | As stated. |
| 31 | K. Ford, "Poisson approximation of prime divisors of shifted primes", arXiv 2408.03803 (v4 Jan 2026), IMRN 2025(7), DOI 10.1093/imrn/rnaf079 | PARTIAL | Mirror: title, author, v4 2026-01-21; abstract confirms the Kubilius analogue, the TV estimate, independent Poisson behaviour on disjoint sets, and the transference principle. **The IMRN venue and DOI are UNVERIFIED (budget).** "Hypothesis Z(gamma)" is not in the abstract (it comes from the report's search 48). |
| 32 | O. Gorodetsky, "A Kubilius model for sieve-theoretic sequences", arXiv 2608.14190 (Aug 2026) | OK(mirror) | Abstract as stated: positive level of distribution, expected qualitatively optimal, recovers Ford, and Tenenbaum's optimal bound up to `x^(o(1))` and `u^(o(u))`. **The TV shape `exp(-c u log u)` with `u = log D/log y` used in §2.5 is the report's reconstruction, not the abstract.** I only use the theorem-independent part of the §2.5 argument (§2.2 below). |
| 15 | J. Merikoski, "On Gaussian primes in sparse sets", Compositio Math. (2025), arXiv 2302.11331 (v3 Sep 2024) | PARTIAL / **ISSUE (scope)** | Mirror: v3 2024-09-24, abstract confirms: `delta > 0` with `\|B cap [1,Y]\| >> Y^(1-delta)`; a quasi-explicit formula for `\|B\| = X^(1/2-delta)`, `delta < 1/10`, `B ⊂ [eta X^(1/2), (1-eta) X^(1/2)]`; an asymptotic unless `B` has a large subset of multiples of a fixed large integer. The abstract also says: **"For an arbitrary B we obtain a lower bound ... with a weaker range for delta."** The report's "uniformity over arbitrary `B` ... robust to adversarial selection" overstates this: arbitrary `B` gives only a lower bound. **The Compositio venue is UNVERIFIED** (the mirror carries no journal reference). |
| 17 | B. Green, M. Sawhney, "Primes of the form p^2 + nq^2", arXiv 2410.04189 (v3 Jun 2026) | OK(mirror) | Exact (v3 2026-06-08; `n = 0, 4 mod 6`; FI Gaussian-primes conjecture; Kuca, Kuca-Kravitz-Leng, Leng-Sah-Sawhney). |
| 26 | K. Oganesyan, "Lattice points on small arcs", arXiv 2107.09991, withdrawn | PARTIAL | Title and number: OK(sib, W2). The withdrawal is from the searchers' snippets only, and is already recorded in `round4/hardness.md:872` and `growth/referee_conditional.md:143` (both checked). It is irrelevant to upper bounds. |

### 1.2 Other 2024-2026 items (all confirmed in the mirror)

| # | citation | status | notes |
|---|---|---|---|
| 6 | Lester, Waxman, Yesha, "The variance of the number of lattice points in narrow sectors", 2609.23214 | OK(mirror) | As stated. **Composite-supported** (all lattice points), so Thm B does not apply (see §0 item 2). |
| 9 | J. Kuś, "Halász theorems for Gaussian ideals in sectors and short intervals", 2605.19105 | OK(mirror) | As stated. **Composite-supported** (all ideals). The report's class-1 cap via Thm B is mis-attributed. |
| 11 | Grimmelt, Merikoski, 2505.00493 | OK(mirror) | As stated. |
| 13 | Lin, Tang, Guo, 2412.14792; Devin 2105.02492 | OK(mirror) | Both in the mirror (Devin as a 2025-10 replacement). |
| 19 | Sabuncu, "Missing digits and sums of two prime squares", 2601.20897 | OK(mirror) | v3 2026-08. |
| 25 | Parkkonen, Paulin, 2602.13058 | OK(mirror) | As stated. |
| 27 | Ferreira, Ribas, 2602.23597 | OK(mirror) | As stated (also: `alpha/pi` is rational or transcendental). |
| 30 | H. Gou, 2604.21122 | OK(mirror) / minor | The GCD result is stated for **`k >= 3`**, which the report omits. The abstract also contains a large-sieve-type complementary estimate and an LCM theorem. The wrong-direction kill stands. |
| — | the RSS-checked "irrelevant" list (2609.13894, 2610.18281, 2309.10971, 2606.29185, 2504.09650, 2508.18044, 2610.34556, 2512.04081, 2604.20014, 2605.24909, 2610.20554) | spot-checked | 2606.29185, 2610.34556 and 2610.20554 read in full: irrelevant as claimed. |

### 1.3 Pre-2024 items (outside the mirror; WebSearch unavailable)

| # | citation | status |
|---|---|---|
| 22 | Kurlberg, Lester, Invent. Math. (2025) 537-585, DOI 10.1007/s00222-025-01324-1, arXiv 2112.08522 | PARTIAL: venue, pages and DOI OK(sib, W8). **The quoted "not i.i.d.-like at scales `N^(-gamma)`, `1/2 < gamma < log(pi/2)/log 2`" is flagged by the sibling checker (W11) as possibly reversed.** Irrelevant to the verdict (density one). |
| 20 | Bombieri, Bourgain, "A problem on sums of two squares", IMRN 2015(11) 3343-3407 | PARTIAL: bibliographic data OK(sib, W17). Content (the 6-correlation `O(N^(7/2))` attribution, and Li's Mathematika 2020 unconditional case) UNVERIFIED. Already capped in the repo (`local_smoothing/chains.md` row E, checked; that row also records Cushman-Demeter-Wu 2608.12316, which the report does not list). |
| 23 | Sartori, arXiv 2005.04698 | PARTIAL: title and "generic integers ... no quasi-correlations" OK(sib, W18). **"Generic is defined through a model in which the angles of the Gaussian primes dividing `n` are i.i.d. uniform" is UNVERIFIED.** The report uses it rhetorically ("generic = fake"); the kill (density one) does not need it. |
| 1 | Huang, Liu, Rudnick, Acta Arith. 193 (2020) 183-192, arXiv 1903.04005 | UNVERIFIED (budget) (the sibling checker could not verify it either) |
| 2 | Järviniemi, Teräväinen, Rev. Mat. Iberoam. 40(4) (2024), arXiv 2303.05822 (widths `(log X)(loglog X)^19.2` and `(log X)^15.1`) | UNVERIFIED (budget) |
| 3-5, 7-8 | Rudnick-Waxman 1705.07498; Chen et al., Exp. Math. 2023; Katz IMRN 2016; Stucky; Harman-Lewis; Harman-Kumchev-Lewis; Matsui | UNVERIFIED (budget) |
| 10, 12, 14 | Musicantov-Zehavi QJM 75(2); Halupczok Funct. Approx. 57(1); Khale et al. 2008.09677; Thorner-Zaman 1510.08086; Brumley-Thorner-Zaman JEMS 24 (2022) | UNVERIFIED (budget). BTZ's described content ("Chebotarev least-prime bounds") is doubtful but irrelevant. |
| 16, 18, 19 | Merikoski PLMS 127(4); Friedlander-Iwaniec JEMS 24 (2022); Pratt 1806.02699; Maynard Invent. 2019; Heath-Brown-Li Invent. 208 (2017) | UNVERIFIED (budget) |
| 21, 24, 28, 29 | Kurlberg-Wigman Math. Ann.; Granville-Wigman CMP 355; Krishnapur-Kurlberg-Wigman Annals 177; Erdős-Hall Discrete Math. 200; Pappalardi-Sha-Shparlinski-Stewart TAMS 370; Corvaja-Zannier background | UNVERIFIED (budget) |

None of the unverified items is used for anything but "irrelevant" or "capped". If any were misattributed, the
verdict would not change.

### 1.4 Items in the report's area that it missed (found by my mirror regex searches)

| arXiv | authors, title (mirror) | what it is | kill |
|---|---|---|---|
| 2508.06005 | S. Fan, "The Hardy-Ramanujan inequality for sifted sets and its applications" | a weighted Hardy-Ramanujan **upper** bound `#{n : omega(n) = k}` on sifted sets; large deviations of `omega` on subsets | Upper tails are the **fake-satisfied** direction (`D_A` has bounded `c_xy`, so it satisfies any upper bound on collision counts). It is also an average over a sifted set with a level of distribution; our family is finite and adversarial (§2.2). |
| 2512.01739 | T. Tao, J. Teräväinen, "Quantitative correlations and some problems on prime factors of consecutive integers" | `omega(n+k) << k` infinitely often; Erdős-Pomerance-Sárközy asymptotic; probabilistic method | Averages over `n`; no worst-case content. |
| 2506.03215 | S. Das, W. Kuo, Y.-R. Liu, "A subset generalization of the Erdős-Kac theorem over number fields" | Erdős-Kac (central range) for subsets of ideals | Central limit range, not the tail `k ~ log/loglog`; an average. |
| 2306.16983 | W. Kai, "Linear patterns of prime elements in number fields" | the number-field Green-Tao-Ziegler theorem (pairwise independent linear parts) | Existence (the wrong direction). Bounded complexity. As an external witness it becomes PCG/GPA (Thm D). |
| 2408.17441 | Frei, König, Tinková, "Linear constellations in primes with arithmetic restrictions" | existence of constellations | Existence (the wrong direction). |
| 2608.05934 | Gao, Li, Qi, "Hybrid Weyl subconvexity over imaginary quadratic fields" | `GL_2 x GL_1` subconvexity in the archimedean aspect of the Hecke character | Feeds prime sums with power-saving error terms, so robust (R2), capped by Thm B in prime-supported form. |
| 2505.00653 | A. Pascadi, "On the exponents of distribution of primes and smooth numbers" | level `x^(5/8 - o(1))` | An ensemble level of distribution. It cannot apply to the finite pair family (§2.2). |
| 2607.02728 | Bosio, Riblet, Tarr, "Large `B_2[g]` subsets of the first squares" | lower bound `>> n^(2g/(2g+1)) (log n)^((2-2^g)/(2g+1))` | Context for CMPSY. Sidon-in-squares is `n^(2/3)` from below vs `n exp(-c log n/loglog n)` from above, so the sieve rate is a generic artefact of the mechanism, not the truth. Consistent with "the mechanism caps at `log/loglog`". |
| 2604.23918 | P. Gao, "On the Gauss circle problem over smooth numbers" | an `r(n)`-weighted count of smooth `n <= x` | An average. Irrelevant. |

Outside the primes area, and already handled by the sibling sweeps (`direct.md` R15/R16, `check_direct.md`): Ziegler
2609.15487 (S-unit equation counts `(\|S\|-1)/2` beyond an `S`-dependent height threshold, i.e. a template-dependent
threshold); Calegari-Dimitrov-Tang 2510.04156; Wang-Xiao 2604.26497.

---------------------------------------------------------------------------------------------------------------

## 2. Kill attempts on each lead (with arithmetic)

The targets are verified against the repo:
* H and Theorem C: `round7/p3lower.md` lines 530-600;
* Theorem T: `round15/rigSLtheory.md` §4.1;
* L1 = (R1-arch), and the split (R1-arch) + (R1-arith) = growth: `round9/README.md` §3 items 2-5.

The report quotes all of them correctly, including:
* "Theorem C uses only one heavy pair at `k = floor((2 - eps/4) log M/loglog M)`";
* "`D_A` has `c_xy <= floor(18/(1-eta))`";
* "RigS is vacuous for actual clusters";
* "H^gen + L1 is a split".

### 2.1 Croot-Mao-Pohoata-Sheffer-Yip (combinatorial large sieve): fails, by a range-independent argument

**What the report claims.** It is the `two_thirds` mechanism, so it is capped. The entropic version uses primes up to
`(log N)^3 = M^(3+o(1))`, which is fixed-`A` residue data.

**Problem.** The prime range is not confirmed (§1.1). If the entropic argument used primes up to `N^c`, the "fixed
`A`" cap would not apply as stated, since `exp(c c_0 M log M) >> M^A`.

**Kill that does not depend on the prime range** (`arith_check.py` parts (c)-(e)).
1. *What counting can force.* A sieve of this type gets its information from collisions forced by local branching
   against the size of `A`. In our setting `A` is the cluster, of size `M`.
   * A collision mod `Q` is forced only if the conic has fewer than `M` classes mod `Q`, i.e. `Q < M` up to constants.
   * A product of `j` window primes is `>= M^(j(1-eta))`. So forced joint collisions need `j < 1/(1-eta)`. This gives
     `j_max = 1, 1, 3, 9` at `eta = 0.25, 0.5, 0.75, 0.9`, independent of `M`.
   * H needs `k = 10, 16` at `M = 10^6, 10^12`, and `k -> infinity` in general.
2. *Entropy adds nothing at moduli `> M`.* For *any* `M`-point set, `log #classes(Q) - H(X mod Q) >= log(Q/M)`, which
   is 13.8 nats at `Q = M^2` and 55 nats at `Q = M^5` (`M = 10^6`). A Shearer-type "deficit" at such moduli is
   automatic for fakes and for actual clusters alike. A collision mod `Q` requires `H(X mod Q) < log M`
   (non-injectivity), and per-prime marginals do not force it.
3. *Fakes satisfy every counting-forced statement automatically.* Pigeonhole and entropy inequalities are theorems
   about *every* residue assignment. So the round-6 design `D_A` satisfies them. `D_A` has `c_xy <= 18/(1-eta)`, so
   `F_k = 0` for `k > 18/(1-eta)`. The mechanism therefore cannot produce H's heavy pair.
4. *Without H, the mechanism gives only the current rate.* The forced log-mass is
   `sum_(p^a < M) E(M, p^a) log p ~ (M^2/2) log M`, against the global bound `sum log\|B_xy\| <= M W/8`. This gives
   `W >= c M log M`: the `two_thirds` inequality `3 M log M <= W + O(M)`, constant only.
   * Joint collisions do not add log-mass beyond the per-prime count, because `v_q(B_xy)` is counted prime by prime.
   * Equivalently, the average pair integer has size `M^(c_0/4)` (`M^0.75` at `c_0 = 3`), with at most
     `c_0/(4(1-eta))` window primes per average pair (part (e)). CMPSY's saving `exp(-c log N/loglog N)` is the
     `omega(n) <= log n/loglog n` count: the same rate.

**Verdict:** **lead fails** (it can only change the constant). The report's conclusion is right. Its stated reason
(fixed `A` via `(log N)^3`) is unconfirmed, and it should be replaced by points 1-4.

### 2.2 Ford / Gorodetsky (Kubilius-Poisson): fails (confirmed, with a corrected formula)

**Fake test.** H is fake-violated (`p3lower` Lemma C.2), so the shape is right. The question is strength and scope.

**Arithmetic** (`arith_check.py` parts (a)-(b); `eta = 1/2`, `eps = 0.1`):

| `M` | `k` | `log E_k` | `log Q_k / log binom(M,2)` | report's formula `(1-eta)L/(2l)` | corrected `(1-eta)(1-eps/8)L/l` | Kubilius `u` | `log10 D` |
|---|---|---|---|---|---|---|---|
| 1e3 | 7 | 2.03 | 1.84 | 0.89 | 1.76 | 6.13 | 18.4 |
| 1e6 | 10 | 8.17 | 2.56 | 1.32 | 2.60 | 8.69 | 52.1 |
| 1e9 | 13 | 13.44 | 3.31 | 1.71 | 3.38 | 11.28 | 101.5 |
| 1e12 | 16 | 18.03 | 4.05 | 2.08 | 4.11 | 13.89 | 166.7 |

* The report's table values (columns 1-4, 7-8) reproduce exactly. Its asymptotic formula does not (factor 2).
* **Other `eta`.** Theorem C holds for every fixed `eta < 1`. The ratio `(1-eta)k/2` at `M = 10^6, 10^9, 10^12, 10^30`
  is:
  * 4.5, 5.9, 7.2, 14.4 for `eta = 0.1`;
  * 1.25, 1.6, 2.0, 4.0 for `eta = 0.75`;
  * 0.50, 0.65, 0.80, 1.6 for `eta = 0.9`.

  It is `<= 1` only while `eta >= 1 - 2/k`, i.e. only if `eta -> 1`. That is not allowed (and Theorem C's constant
  `8(1-eta)` would then vanish). So no choice of window rescues the reduction.
* **The decisive point does not depend on the TV-bound shape** (which is unconfirmed). Some pair with `c_xy >= k`
  means one pair integer divisible by a product of `k` window primes, of size `>= M^((1-eta)k)`. That is `10^30` at
  `M = 10^6`, while the family has `binom(M,2) ~ 5 * 10^11` members.
  * Any distribution statement about a finite family at moduli above its size has counts 0 or 1 against expectations
    `< 1`, so it cannot be "Poisson".
  * The family is adversarial (selected by the cluster), and Ford/Gorodetsky are statements about an ensemble with a
    natural level of distribution.

**Verdict:** **lead fails.** It is the correct shape (H), but it would need to *be* L7/H for adversarial families,
which is the open problem.

### 2.3 Merikoski (sparse `B`): fails, but the report's size argument has a hole

* **The report says the cluster selection is "always too sparse"**: `M = X^(1/h_0)` with `h_0 = 3, 5, 9` against
  `\|B\| > X^0.4`. Since `c_0 >= 3` (`two_thirds`), `h_0 = 2c_0 - 3 >= 3`, so this holds at the **top** height.
* **The hole.** The cluster's primes need not sit at the top height. If they are all `<= M^h` with `h <= 2.5`, then
  `M >= X^0.4` and `B` can be padded to the required size (`arith_check.py` part (f): in range for `h = 1.2, 2.0, 2.5`).
  This includes the uncapped low-height window `P_+(N) < M^(2-delta)` (`analytic.md` §6, open point 1), which is the
  one place where a size mismatch would matter.
* **Correct kill:**
  1. The theorem counts primes `a^2 + b^2 <= X` with `b in B`. Its main term is `~ \|B\| X^(1/2)/log X >= X^0.9/log X`:
     `2.3 * 10^9` vs `M = 10^6` at `h = 2`, and `> M` once `h > 1.12`. A padded `B` hides the cluster's `b_j`.
  2. It constrains the imaginary coordinates of primes, not linear relations among their angles mod `pi/2` at
     precision `e^(-W/4)`, which is what clusters carry.
  3. For **arbitrary** `B` the abstract gives only a **lower bound** (with a weaker `delta`). That is the existence
     direction: used as external witnesses, it becomes PCG/GPA and collapses to the problem (`analytic.md` Thm D).

**Verdict:** **lead fails.**

### 2.4 Class 1 (ensemble statistics): the conclusion holds, but the cap is mis-attributed for composite-supported items

* **Prime-supported items: Thm B applies.** Huang-Liu-Rudnick sector counts, Harman-Lewis, Bombieri-Vinogradov in
  sectors (Halupczok, Khale et al.), Musicantov-Zehavi, Hecke prime sums, and Gao-Li-Qi (§1.4).
  * The `A` condition is `A > max(2, 1/gamma)` **and `A >= 4`** (`analytic.md` Thm B; `round10/README.md` §7.1 writes
    `A > max(4, 1/gamma)`). The report drops `A >= 4`. This is harmless.
* **Composite-supported items: Thm B does not apply.** These are Kuś (all ideals), Lester-Waxman-Yesha (all lattice
  points), Järviniemi-Teräväinen (almost primes, Matomäki-Radziwiłł), and Green-Sawhney (Type I/II).
  * `referee_analytic.md` Slip 3: "only prime-supported inputs are covered". A multiplicative graft changes
    `~ X s/P_0` ideals, which is not `o(X^gamma)` for large `X`. `round10/README.md` §7.2 lists the class as
    unclassified.
  * Each item still fails:
    * LWY is a variance averaged over sector position and over all norms `<= X`;
    * Kuś is a mean value over all ideals in a sector;
    * JT is an "almost all sectors" statement about `P_2`/`P_3`, objects the cluster does not contain (cluster points
      have `~ M` prime factors, cluster primes are primes);
    * Green-Sawhney is existence.
  * The averages exclude a single `N` only if `Err < 1` (`analytic.md` §4.2). The cluster's circle carries
    `r(N) = N^(o(1))` points, since `log r(N)/log N ~ (M/4) log 2/(c_0 M log M) -> 0`.
* **Consequence for the record.** "Composite-supported analytic inputs" must stay **unclassified** in the barrier map.
  No found item is a lead, but the class is not capped.

### 2.5 Classes 3-4 (generic `n`; Diophantine and gcd): conclusions hold

* **Density one.** Kurlberg-Lester, Sartori, Granville-Wigman and Erdős-Hall are silent on the density-zero
  near-extremal set (`analytic.md` §4.2, verified).
* **Energy upper bounds** (Bombieri-Bourgain, and Cushman-Demeter-Wu 2608.12316 for every finite `X` on every strictly
  convex curve) are fake-satisfied: `local_smoothing/chains.md` row E (verified).
* **Gou / Green-Walker.** At `k = 2`, `D = N^(1/2) M^(-b)`: `\|A\|^2 << N^2/D^2 = N M^(2b)`, trivial for `\|A\| = M`.
  (Gou's own theorem is `k >= 3`.) Wrong direction for L1.
* **Corvaja-Zannier/Levin:** fixed `S`, so forbidden by the constraints. Ferreira-Ribas: single form, capped.
* **Constraint compliance of the report:** checked against the four project constraints (no fixed-`S` uniform
  exceptional sets, Gaussian gcds kept, no zeroed drifting phase, no template-dependent threshold). No violation
  found.

### 2.6 The structural claim (H^gen + L1 is a split)

Verified.
* `round9/README.md` Thm 3.2: on near-extremal clusters with `K`-uniform margins, each `B_xy` has fewer than
  `K c_0/(4(1-eta)) + o(1)` window prime factors.
* Cor. 3.4: (R1-arch) + (R1-arith) equals growth.

Directly: L1 gives `\|B_xy\| <= M^b`, so `c_xy <= b/(1-eta)`, so `F_k = 0` at `k ~ 2 log M/loglog M`, so H fails. So
the literature would have to prove one half against adversarial clusters:
* the fake-violated half is H (heavy-pair existence: §2.1-2.2 show that counting and Kubilius cannot reach it);
* the other half, L1 (no pair gcd smaller than `N^(1/2) M^(-O(1))`), is fake-satisfied and is not supplied by any
  gcd-extremal result (§2.5).

---------------------------------------------------------------------------------------------------------------

## 3. Papers to request in full text (re-ranked)

None is expected to change the verdict. The order is by how much the record would gain.

1. **Croot-Mao-Pohoata-Sheffer-Yip, arXiv 2606.17487.** Only to settle which primes and moduli the entropic (Shearer)
   argument uses, and to confirm that its forcing comes from moduli below `\|A\|` (the dense regime). That would fix
   the report's unsupported "`(log N)^3`" sentence. §2.1 kills the lead either way.
2. **Gorodetsky, arXiv 2608.14190.** Only to confirm the TV-bound shape behind the `u`/`D` columns. The kill in §2.2
   does not depend on it.
3. **Sartori, arXiv 2005.04698, and Kurlberg-Lester, arXiv 2112.08522.** Record-keeping only: the genericity model,
   and the direction of the `log(pi/2)/log 2` range.
4. **Merikoski, arXiv 2302.11331.** Only the arbitrary-`B` lower-bound statement and its `delta` range, to fix the
   "robust to adversarial selection" wording.

Dropped from the report's list:
* Bombieri-Bourgain: already capped, row E;
* Rudnick-Waxman: the lattice fact of Thm D;
* the remainder of item 4: the kill is independent of method details.

---------------------------------------------------------------------------------------------------------------

## 4. Search log

**WebSearch** (three calls attempted, all refused for budget, no results):
* [E] `Croot Mao Pohoata Sheffer Yip "combinatorial large sieve" Sidon sets distances norm forms`
* [S] `Huang Liu Rudnick "Gaussian primes in almost all narrow sectors" Acta Arithmetica`
* [S] `Järviniemi Teräväinen "Gaussian almost primes in almost all narrow sectors"`

**Local RSS mirror** (`tools/db.json`; my scripts `check_primes_work/show.py`, `grep.py`, `grepnt.py`, `titles.py`;
NT/CO unless stated):
1. ID look-ups for every arXiv number in `primes.md` (47 IDs): 2024+ items found; pre-2024 IDs absent (outside the
   window).
2. `(Erdős-Kac|Poisson|Kubilius)` and `(binary form|polynomial values|differences|norm form|quadratic form)`, all
   categories: noise only.
3. `(Kubilius|Erdős-Kac|number of prime factors|prime divisors of)` and
   `(polynomial|binary form|shifted|sieve|level of distribution|Gaussian|norm)`: 41 hits; found 2508.06005,
   2512.01739, 2506.03215.
4. `lattice points` and `(arc|circle|x^2+y^2)`: 6 hits; nothing new beyond Waxman-Yesha, Cushman-Demeter-Wu (both in
   the sibling report) and Gauss-circle numerics.
5. `(S-unit|unit equation)` and `(uniform|number of solutions|bound)`: 11 hits; Ziegler, Wang-Xiao and CDT, all
   handled by the siblings; Ostafe-Pomerance-Shparlinski 2503.03985 (counts over equations, average).
6. `gaussian primes?` (all math): 11 hits; nothing missed except irrelevant items (Goldbach averages, moat, Chen
   analogue).
7. Hecke characters/L-functions with `(Q(i)|Gaussian|imaginary quadratic|angle|sector)`: 18 hits; only Gao-Li-Qi
   2608.05934 is pertinent (robust).
8. `(sums of two squares|x^2+y^2|a^2+b^2)`: 111 hits; found Bosio-Riblet-Tarr 2607.02728 and P. Gao 2604.23918;
   others irrelevant.
9. `(large sieve|larger sieve)`: 72 hits; nothing pertinent beyond the report's items and Pascadi 2505.00653.
10. `(linear equations|constellation|Gowers|patterns)` with `(Gaussian|number field|imaginary quadratic|Z[i])` and
    `prime`: Kai 2306.16983 and Frei-König-Tinková 2408.17441.
11. `(smooth|prime factors|greatest prime factor)` with `(binary form|values of ... form|norm form|resultant|Lucas|Lehmer)`:
    nothing pertinent.
12. `(abc|a+b=c|S-unit equation)` with `(smooth|uniform|count)`: Croot-Mao-Yip 2603.14654 (Hilbert cubes; ABC-conditional
    bounds for `k`-th powers in AP, `k >= 4`): outside the primes area, not about squares.

**Sibling records consulted** (independent snippets): `check_direct.md` W2, W8, W11, W17, W18 (Oganesyan title,
Kurlberg-Lester venue, the Kurlberg-Lester direction caveat, Bombieri-Bourgain bibliographic data, Sartori title).
