# lit18: judge's summary of the 2015-2026 literature sweep for growth-rate inputs

**Inputs.** I read three sweeps and three adversarial checks:
* `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/lit18/direct.md` (the direct problem;
  query log in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/lit18/direct_log.txt`)
* `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/lit18/tools.md` (height-based tools
  uniform in `S`)
* `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/lit18/primes.md` (genuine-prime
  statistics: H^gen, L1, L7)
* `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/lit18/check_direct.md`
* `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/lit18/check_tools.md`
* `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/lit18/check_primes.md`

**My own checks.**
* I re-read the abstracts of Ziegler, Chan (two papers), Hegyvári, Croot-Mao-Pohoata-Sheffer-Yip, Letendre and
  Cushman-Demeter-Wu in the offline arXiv mirror (script: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/lit18/judge/lk.py`).
* I re-derived Prop B of `check_tools.md` and the Gabdullin and Letendre corrections of `check_direct.md`.
* I ran a crude random-model count for the one surviving hypothesis, in
  `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/lit18/judge/qsprod_heur.py` and
  `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/lit18/judge/qsprod_heur2.py`.
* Repo files were read only, never written.

**Verification limits.**
* About 200 web searches were run in total: 60 + 69 + 48 by the three searchers and 21 by `check_direct`. The shared
  per-turn budget then ran out, so `check_tools` and `check_primes` ran no web searches.
* Post-2024 arXiv items were verified against an offline mirror of the daily arXiv RSS feed. It holds 92,926
  abstracts from 2024-03-05 to 2026-10-27 in 11 math categories, and was downloaded from GitHub
  (`ehijano/rss_fetch@30c967b`).
* Most pre-2024 journal details (volumes, pages) rest on one searcher's snippet. They are marked **S** below.
* No item that carries the verdict depends on an unverified detail.

---

## 0. Bottom line

**Nothing found plausibly beats the barrier, and no literature lead survives the checks.**

1. No paper proves `o(log R/loglog R)`, or even `o(log R)`, for every `R`. The best literature bound valid for all
   `R` is Bourgain-Rudnick's `O(log R)`, which is weaker than the program's own `(2/3+eps) log R/loglog R`.
2. Every item found falls into one of five groups:
   * bounded only near a direction of bounded height: Chan, Erdős-Rosenfeld, and Gabdullin once the position is
     corrected;
   * holds only for almost all `R`, or is an energy, correlation or ensemble statement that fake circles satisfy;
   * fixed-`S` or template-dependent, which the project constraints exclude;
   * the same counting mechanism as the proved bound, at the same rate (the 2026 "combinatorial large sieve");
   * the right shape but impossibly weak (Kubilius/Poisson models for H).
3. The sweep's real output is three barrier-map facts:
   * the round-8 Subspace route `QS` is dead for every threshold that depends on the prime set only through its size;
   * composite-supported analytic inputs remain unclassified, not capped;
   * one open, non-standard hypothesis, `QS^prod`, survives, and it is at least as strong as the target.

---

## 1. Results found, with status after the checks

**Status codes.**
* **W**: confirmed by an independent checker's web snippet (`check_direct.md` W1-W21).
* **A**: abstract read in the offline arXiv mirror, by a checker or by me.
* **F**: verbatim from a fetched file (the GitHub paper card).
* **S**: searcher's snippet only, not re-verified.
* **R**: already in the repo record.
* **!**: a check corrected the searcher's reading.

### 1.1 The direct problem: arcs, close lattice points, divisors near `sqrt n`

| item | status | content as confirmed | verdict and what caps it |
|---|---|---|---|
| Bourgain-Rudnick, arXiv:1012.3843, 1105.0018, 1308.4247 | W (titles; `O(log m)` points on arcs `<= m^(1/4)`; "known for almost all λ, individually only `B_λ << log λ`"); S (separation count `O(X^(1-eps/3))`, AHP venue) | best all-`R` bound in the literature is `O(log R)` | weaker than ours. The separation result is almost-all: `round10/analytic.md` Thm B, `local_smoothing/chains.md` row G |
| K. Oganesyan, "Lattice points on small arcs", arXiv:2107.09991 | W (title, number, claim of unboundedness on arcs `R^alpha`, `alpha in (1/2,1)`); S, three sweeps independently (withdrawn, "mistake in Lemma 2"); R | lower-bound claim above `sqrt R`, withdrawn | irrelevant to upper bounds. A correct version would have contradicted the Vojta instance in `growth/conditional.md` (arcs `C R^alpha`, `alpha < 1`) (`check_direct` §2.8) |
| M. R. Gabdullin, arXiv:2205.13611 (JFAA 30 (2024) per `round4/hardness.md`) | W; R (`chains.md` row F, `round4/hardness.md`); **!**; CANT 2022 talk unverified | `||f||_4 << ||f||_2` for frequencies `{n^2 : N <= n <= N + N^gamma}`, `gamma < (sqrt5-1)/2`, via divisors in short intervals | **dead twice.** (a) It is an energy statement, vacuous on one cluster (`M^2 <= 4c^4 M^2`). (b) **!** It concerns arcs near the diagonal, a direction of height `sqrt 2`. Under the Cilleruelo-Granville rotation, `gamma` corresponds to general-position arcs `R^(2gamma-1)`, so `gamma < 0.618` gives `R^0.236`, already inside Cilleruelo-Córdoba. The target `C sqrt R` needs `gamma = 3/4`. A pointwise version of the lemma would also be useless |
| P. Letendre, "Divisors of an integer in a short interval", arXiv:2503.12146 | W, A; Thm 1 shape S; **!** | `D_n(n^theta, n^eta) << tau(n)^(1-xi) V(n) log tau(n)`. The forum's specialisation gives `<< 1/delta` divisors in `(sqrt n, sqrt n + n^(1/4-delta))` | **!** That window is trivial: `d + n/d` is an integer in an interval of length `Y^2/sqrt n`, so there are at most 2 divisors. It is not the Cilleruelo-Córdoba analogue. Above `n^(1/4)` the bound is only `tau(n)^(1-xi)`. The rational problem at its critical scale is open too (#886) |
| Erdős #886 (Ruzsa) and #887; Erdős-Rosenfeld `1 + C^2` | W (erdosproblems snippets); S (Acta Arith. 1997 venue; formal-conjectures PR #7218, which is noise); R | open in general position | height-1 position only |
| T. H. Chan, IJNT 11(5) (2015) 1701-1708, doi 10.1142/S1793042115400205, Thm 3 | F (verbatim), W; R (`round4/outside.md`, `round4/hardness.md`, `growth/conditional.md` §5); **!** | at most 10 points with `\|b\| < n^(1/4)(log n)^(1/7)` for `n = N^2` | near-axis only. After rotation it misses general position by a factor `R^(1/2)/(log R)^(1/7)`. **!** Chan's arXiv:2507.07094 "correct[s] a mistake of a previous result", without saying which one (see §3) |
| T. H. Chan, arXiv:2507.07094; Chan-Holmes-Liu-Villarreal, arXiv:2508.02818 | A (missed by `direct.md`) | close factorizations; lattice points near the vertex of `xy = N`; `A <= 0.04742 C^3 + O(C)` for four factorizations | vertex (height 1), or the Jarník `R^(1/3)` scale; fixed counts. Dead |
| F. Temur, Turkish J. Math. 46(7) (2022) 2867-2889, DOI 10.55730/1300-0098.3306 | W, including the author's "gives no new results on Conjecture 1" | Schmidt's theorem on fixed gap templates | template-dependent, ineffective threshold: excluded by the constraints. Dead |
| Cushman-Demeter-Wu, arXiv:2608.12316; Li-Yung, arXiv:2609.33949 | A; R (`chains.md` row E) | `J_3(gamma(X)) << \|X\|^(3+eps)` on strictly convex curves | for a `C sqrt R` arc this gives `M << R^(1/6+eps)`, and fakes satisfy it. Capped (energy class) |
| Kurlberg-Lester, Invent. Math. (2025) 537-585, DOI 10.1007/s00222-025-01324-1 | W; **!** volume 240 unconfirmed; the quote "not i.i.d. at scales `N^(-gamma)`" is unconfirmed and possibly reversed | Poisson spacing along a density-one subsequence | density one. Natural scale `R^(-o(1))`, against the cluster's `R^(-1/2)`. Thm B / row G |
| Granville-Wigman, CMP 355 (2017) 767-802; Sartori, arXiv:2005.04698; Klurman-Sartori, arXiv:2104.13038; Benatar-Marinucci-Wigman, JAM 141 (2020) 707-749; Cammarota-Klurman-Wigman, CMP 376 (2020) 1261-1310 | W (bibliographic data, "density one" / "generic"); S (Sartori's i.i.d.-angle model, Klurman-Sartori Lemma 1.5) | statements for generic or density-one `n` | silent on the density-zero set of `N` that carries large clusters. Dead |
| Bombieri-Bourgain, IMRN 2015(11) 3343-3407 | W (bibliographic data); S (content) | 6-correlations of lattice points on circles | energy class (row E) |
| Squares in arithmetic progressions | W (Elsholtz-Wurzinger, QJM 75(4) 2024); A (González-Jiménez 2602.03251, González-Jiménez-Tho 2602.01380, Hegyvári 2504.13230) | no improvement on Bombieri-Zannier `k^(3/5)` found | only the direction arcs ⇒ squares is known. **!** Hegyvári's abstract says the result is essentially Mei-Chu Chang's (JFA 207 (2004)) and that the file will be withdrawn |

### 1.2 Height-based Diophantine tools

| item | status | content | verdict and what caps it |
|---|---|---|---|
| V. Ziegler, "On the number of solutions to S-unit equations", arXiv:2609.15487 | A (I read the abstract); **!** | for fixed `a_1, a_2`, `a_1x + a_2y = a_3` has `<= (\|S\|-1)/2` solutions in `S`-units once `h(a_3)` is large; effective | dead. The cluster's right-hand sides vary: `u + u^(-1) = 2 - h_jk/N`, and the triangle identities have varying small-height coefficients. The threshold depends on `S`. **!** The count is `>= s` (two places per split prime), up to `(2+o(1)) log R/loglog R`, so even a valid application would be 3 times worse than the proved bound |
| Calegari-Dimitrov-Tang, arXiv:2510.04156; Ferreira-Ribas, arXiv:2602.23597 | A | new effective method; arguments of unimodular algebraic numbers are Diophantine | single-form, capped (`round8/lw.md` Thms A-B; `round8/stability.md` Prop 6.1). The trivial bound `\|Im(z_j conj z_k)\| >= 1` is already stronger than a Baker-type `N^(-c')`. **No measure of algebraic independence for about `n` logarithms appears anywhere in the mirror** (the round-16 remark) |
| gcd / Vojta-for-blowups: Levin (Invent. Math. 2019); Corvaja-Zannier; Bugeaud-Corvaja-Zannier; Grieve-Wang (TAMS 2020); Levin-Wang; Wang-Yasufuku (ANT 2021); Huang-Levin 2207.14432; Huang-Levin-Xiao 2406.18879; Wang-Xiao 2506.03534 and 2604.26497; Xiao 2110.01751; Rousseau-Turchet-Wang (Math. Ann. 2024); Bugeaud-Evertse-Győry (2018); Chen-Garcia-Fritz-Mathur-Pasten 2602.06956 | S for the journal items; A for 2406.18879, 2506.03534, 2604.26497, 2602.06956; **!** | gcd and degeneracy bounds for fixed `S` | fixed `S`: the exceptional sets depend on `S` (constraint), and the quantitative forms are exponential in `s`. **!** These theorems bound gcds at places *outside* `S`. A cluster's large quantity is proximity at the archimedean place, which is in `S`, so it is governed by the Subspace system (5.1), not by them |
| single-place theorems: Heier-Levin (AJM 2021; Crelle 2024); Ru-Vojta (AJM 2020); Ru-Wang (IJNT 2022); Vojta (JTNB 2023); Yasufuku (Monatsh. 2025; BLMS 2026); Barrios 2411.06050 and Correction; Zhao 2502.08049 | S; A for Barrios and Zhao (Zhao missed by the searcher) | Seshadri- or beta-weighted approximation to subschemes with `S = {inf}` | capped by `tau* < 2` (`threshold/upper.md`). Heier-Levin would need `eps_Y > M+1`, but `eps_Y <= 1`. Barrios needs an ample normal bundle, but `deg N_(Y/X_M) = 0` |
| quantitative Subspace: Evertse-Ferretti (Ann. of Math. 177 (2013) 513-590; survey 1008.2268); Evertse-Schlickewei; Grieve 2306.16583 | S / C; **!** | snippet: "in many cases R can be taken independently of s" | as worded, the snippet contradicts Props A and B (§2.2) unless `s` enters the threshold or the count. Full-text request #1 |
| unit equations: Evertse-Schlickewei-Schmidt (Ann. 2002); Amoroso-Viada (Duke 2009); Beukers-Schlickewei (Acta Arith. 1996). Lower bounds: Erdős-Stewart-Tijdeman, Konyagin-Soundararajan, Granville, Harper | S / C | counts exponential in the rank; lower bounds super-polynomial | polynomial-in-rank counts are false. Capped |
| dimension growth: Salberger; Castryck-Cluckers-Dittmann-Nguyen; Binyamini-Cluckers-Novikov; Binyamini-Cluckers-Kato 2401.03982. Effective Pila-Wilkie and Wilkie: Binyamini-Novikov(-Zak), Binyamini-Jones-Schmidt-Thomas | S; A for BCK; **!** arXiv:2302.04209 is labelled BCN in `tools.md` but CCDN in `direct_log.txt` | polynomial in the degree | trivial at degree 2; reduces to pair data (`growth/auxiliary.md`). Pila-Wilkie does not apply: the transcendental part is empty |
| Z. Huang (ANT 2021); Lehmann-McKinnon-Satriano 2403.02480 | S; A | approximation constants at a fixed point | fixed target point, but the cluster's phase drifts (constraint). Dead |

### 1.3 Genuine-prime statistics

| item | status | content | verdict and what caps it |
|---|---|---|---|
| prime-supported ensemble statistics: Huang-Liu-Rudnick (Acta Arith. 193 (2020)); Rudnick-Waxman; Chen et al. (Exp. Math. 2023); Katz (IMRN 2016); Stucky; Musicantov-Zehavi (QJM 2024); Halupczok; Khale et al.; Thorner-Zaman; Brumley-Thorner-Zaman; Grimmelt-Merikoski 2505.00493; Lin-Tang-Guo; Devin; Gao-Li-Qi 2608.05934 (missed) | S for pre-2024 items; A for 2024+ items | sector counts, variances, Bombieri-Vinogradov in sectors, Hecke zero density | robust, so capped by `round10/analytic.md` Thm B (`A > max(4, 1/gamma)`). An "almost all sectors" exceptional set contains any `O(M)`-point selection (Prop C) |
| composite-supported inputs: Kuś 2605.19105; Lester-Waxman-Yesha 2609.23214; Järviniemi-Teräväinen (almost primes); Green-Sawhney 2410.04189 | A (2024+); **!** | Halász over all ideals; variance over all lattice points; Type I/II sums | **!** Thm B does **not** cover these (`round10/README.md` §7.2; `referee_analytic.md` Slip 3), so the class stays **unclassified**. Each item fails on its own terms: it averages over position or norm, holds for almost all sectors, or is an existence statement |
| J. Merikoski, arXiv:2302.11331 (Compositio venue S) | A (partial); **!** | quasi-explicit formula for `\|B\| = X^(1/2-delta)`, `delta < 1/10`; for arbitrary `B`, **only a lower bound** | dead. It counts primes with `b in B` (the main term is far larger than `M`) and says nothing about angle relations. It is an existence statement, which as an external witness collapses to the problem (`analytic.md` Thm D, Prop D.1). **!** In the low-height window `P_+(N) < M^(2-delta)` its size condition *can* be met, so "too sparse" is not the reason it fails |
| thin-set existence: Merikoski (PLMS 2023); Friedlander-Iwaniec (JEMS 2022); Pratt; Sabuncu 2601.20897; Maynard; Heath-Brown-Li; Kai 2306.16983; Frei-König-Tinková 2408.17441 | S / A | primes in thin sets; linear patterns | wrong direction (existence). Dead |
| Kubilius/Poisson models: Ford, arXiv:2408.03803 (IMRN venue S); Gorodetsky, arXiv:2608.14190; Fan 2508.06005; Tao-Teräväinen 2512.01739; Das-Kuo-Liu 2506.03215 | A; **!** | Poisson laws for prime factors at a positive level of distribution | **right shape for H, impossible strength.** One heavy pair needs a pair integer divisible by `>= M^((1-eta)k)`, with `k ~ 2 log M/loglog M`. That is `10^30` at `M = 10^6`, against a family of `binom(M,2) ~ 5*10^11` adversarially chosen integers. **!** The log-ratio is `(1-eta)(1-eps/8) log M/loglog M`, twice `primes.md`'s formula. Upper-tail results such as Fan's are in the direction fakes satisfy |
| combinatorial large sieve: Croot-Mao-Pohoata-Sheffer-Yip, arXiv:2606.17487; Croot-Mao-Pohoata-Yip 2607.15311; Croot-Mao-Yip 2510.08862; Bosio-Riblet-Tarr 2607.02728 | A (I read 2606.17487); **!** | "local congruence branching produces many modular collisions, while global bounded-multiplicity hypotheses force these collisions to be rare"; Sidon subsets of squares `<= N exp(-c log N/loglog N)` | the same mechanism as `two_thirds/proof.md`, at the same rate. **!** "Primes up to `(log N)^3`" belongs to 2607.15311, not to this paper. The kill does not depend on the prime range: counting and entropy force collisions only at moduli `< M`, i.e. for at most `1/(1-eta)` window primes jointly, while H needs `k -> inf`. The round-6 fake design `D_A` satisfies every counting-forced statement |
| Pappalardi-Sha-Shparlinski-Stewart (TAMS 2018); Gou 2604.21122 (stated for `k >= 3`); Parkkonen-Paulin 2602.13058; Kurlberg-Wigman; Krishnapur-Kurlberg-Wigman; Erdős-Hall | S / A | averages; wrong direction for L1 | dead |

---

## 2. Surviving leads

### 2.1 From the literature: none

All three sweeps and all three checks agree, and I found no reason to disagree. Each sweep's closest candidate dies at
a named statement:
* **Gabdullin:** wrong position, and averaged.
* **Croot et al.:** the same mechanism as the proved bound.
* **Ford/Gorodetsky:** the needed modulus exceeds the size of the family by an unbounded power.
* **Ziegler:** varying right-hand side, an `S`-dependent threshold, and a count `>= s`.

### 2.2 Barrier-map facts produced by the sweep (unrefereed; not leads)

1. **Prop A** (`tools.md`, recounted independently by `check_tools.md`). Fix a threshold `H_0` that does not depend
   on `S`. Then round 8's `QS(m, delta, f, H_0)` forces `f(s) >= exp(s^(1/2-o(1)))`. Witness: smooth `S`-integers
   in common small squares. This corrects `round8/README.md` §2.3/§2.5, which says "for `m >= 4` neither refuted nor
   supported".
2. **Prop B** (`check_tools.md` §2.2; I re-derived it).
   * *Statement.* Fix any threshold `H_0(s)` that depends on `P` only through `s = |P|`. Then for every `m` and every
     `delta < (m-1)/2` there is a `P` for which the solutions of (5.1) need at least `s - m + 1` subspaces.
   * *Witness.* `s` Gaussian primes in an octant square of side `x^(theta')`, with `x >= H_0(s)` and
     `(m-1)(1-2theta') > delta`.
   * *Check of one case.* With `m = 2`, side 15 and `log H = 22`, `delta_eff ~ 1 - 2 ln15/22 = 0.75`. This matches
     the computed 0.749.
   * *Consequence.* Through `weakest.md` Thm 5.2, the whole `QS` route can give at best `M <= omega(N')`. That can be
     `(2+o(1)) log R/loglog R`, worse than what is proved.
   * *What would escape.* An `S`-integer input must use equality of norms, or a threshold tied to `prod_(p in P) p`.
3. **Composite-supported analytic inputs** (all-ideal large sieve, Type I/II sums) remain **unclassified** (§1.3 row 2).
   The sweep found nothing in this class that bears on a single `N`.

### 2.3 The one open hypothesis: `QS^prod` (not a literature result; at least as strong as the target)

**Statement** (`check_tools.md` §2.5). For some `m >= 2`, `0 < delta < (m-1)/2`, `0 < c < 1/2` and `f(s) = o(s)`:

> for every finite set `P` of split primes with `s = |P|`, the solutions of round 8's system (5.1) with
> `H(x) >= max(H_0, (prod_(p in P) p)^c)` lie in at most `f(s)` proper subspaces of `Q(i)^m`.

**What it would give for our clusters.**
* Cluster tuples satisfy the threshold: `H = N/|g|^2 >= N^(1/2)/C^2 >= (prod_(p | N) p)^c` once `N^(1/2-c) >= C^2`.
* `weakest.md` Thm 5.2 then gives `M <= f(omega(N')) + m - 1`. So `f = o(s)` gives `M = o(log R/loglog R)`, and a
  bounded `f` gives the uniform bound.

**Which barriers it escapes, and why.**
* *Fakes.* Fake circles are not `S`-integers, so it is fake-violated in principle: the finite part of (5.1) is
  exactly `H^(-m)` only through the product formula.
* *Prop A.* The dense tuples have height `<= exp(O(sqrt(s log s)))`, far below `(prod p)^c = exp(c s log s)`.
* *Prop B.* Gaussian-prime clusters have height `~ (prod p)^(1/s)`.
* *Conjugation patterns.* These have relative distance `~ H^(-1/(2J))`, too large to count as close.
* *Constraints.* It is a hypothesis with an explicit count, not an exceptional set imported from a fixed-`S` theorem.

**Why it is not progress.**
* It covers unequal-norm tuples, so it says strictly more than Theorem 5.2 needs.
* No literature on it was found.
* It belongs to the genuine-integer class (H^gen; `round16/README.md` §7.2).

**Judge's heuristic check** (`judge/qsprod_heur2.py`; crude, unrefereed).
* *Model.* Take `m = 2`, `P` = the first `s` split primes, and a Rankin/saddle-point count of `P`-smooth Gaussian
  integers. Then take the expected number of pairs within `|x|^(1-delta)` of each other, maximised over heights
  `>= (prod p)^c`.
* *Results* for `c = delta = 0.49`:

  | `s` | expected log-count |
  |---|---|
  | 50 | `+1.41 s` |
  | 1000 | `+0.43 s` |
  | 3000 | `+0.13 s` |
  | 10^4 | `-0.20 s` |

  For `c = 0.25` it is still `+0.21 s` at `s = 10^4`.
* *Reading.* The model supports `QS^prod` only asymptotically. The margin is `c delta s log s` against `O(s)`, a
  single `log s` factor. In this model any admissible `f` must be `exp(Omega(s))` up to about `s ~ 5*10^3`, i.e.
  `log N >~ 5*10^4`.
* *Small structured families.* Prouhet-Tarry-Escott pairs give one solution per `P` at `H ~ (prod p)^(1/2)`.
  `J`-chains reach only `(prod p)^(1/J)`. Neither refutes the hypothesis.

**First concrete task (cheap, decisive in one direction).** Try to refute `QS^prod(2, delta -> 1/2, c -> 1/2)`.
* *Goal.* Find a family of `P` with `>> s` distinct projective solutions at height `>= (prod p)^c`.
* *Candidates:*
  * Erdős-Stewart-Tijdeman / Konyagin-Soundararajan many-solution `S`-unit configurations lifted to `Z[i]`, after
    checking their heights against `(prod p)^c`;
  * multiplicative combinations of PTE chains;
  * the random model above, with `P` chosen adversarially rather than as the first `s` primes.
* *Outcome.* A refutation closes the `QS` route completely. Failure leaves `QS^prod` as a recorded open,
  non-standard hypothesis. **Do not spend proof effort on it.**

---

## 3. Papers to fetch in full text (ranked; none is expected to change the verdict)

1. **Evertse-Ferretti, Ann. of Math. 177 (2013) 513-590, with the survey arXiv:1008.2268.** Find where `s = |S|`
   enters the quantitative Subspace corollary: the threshold, the count, or the number of forms. Props A and B
   predict an `s`-dependent threshold or count. Any other reading means our new barrier propositions are misread.
2. **Croot-Mao-Pohoata-Sheffer-Yip, arXiv:2606.17487.** Find which primes and moduli the entropic (Shearer) argument
   uses, and whether it uses joint factorization data of norm-form values beyond per-prime residues plus one global
   multiplicity bound. Only the latter is fake-satisfied. Expected result: constant-level only.
3. **T. H. Chan, arXiv:2507.07094.** Find which earlier result is corrected. The record quotes Chan's IJNT 2015
   Thm 3 in `round4/outside.md`, `round4/hardness.md` and `growth/conditional.md` §5. This is record hygiene.
4. **V. Ziegler, arXiv:2609.15487.** Find the exact dependence of the height threshold on `S` and `(a_1, a_2)`, and
   whether any version lets `a_3` range over a family (`2N - h`). This documents the barrier map. Expected result:
   dead.
5. **Wang-Xiao, arXiv:2604.26497 (with Huang-Levin-Xiao, arXiv:2406.18879).** Check whether the degeneracy for
   one-parameter families of unit equations controls its exceptional set uniformly in `S`, or by heights tied to
   `prod p`. A `QS^prod`-type statement would be the only interesting outcome. Expected result: fixed `S`.
6. **J. Merikoski, arXiv:2302.11331.** Read the arbitrary-`B` theorem (a lower bound, and its `delta` range) and
   check whether anything controls angles jointly. This matters only for the uncapped low-height window
   `P_+(N) < M^(2-delta)`. Expected result: no.

**Dropped.**
* Gabdullin: dead even pointwise.
* Letendre: the specialisation sits in the trivial window.
* Temur: the author's own statement suffices.
* Bombieri-Bourgain: row E.
* Kurlberg-Lester and Sartori: record hygiene only.
* Ford and Gorodetsky: the kill does not depend on the shape of the total-variation bound.

---

## 4. Is any lead stronger than the program's current options?

**No.**
* **Route 3 (H^gen + L1).** Nothing found supplies either half.
  * Ford and Gorodetsky have the shape of H but fall short by an unbounded power.
  * Croot et al. reproduce the proved bound's mechanism.
  * Merikoski and the thin-set results point in the existence direction.
  * Every gcd theorem is either fixed-`S` or points the wrong way for L1.
  * `QS^prod` is in the same genuine-integer class and is stronger than what route 3 needs. It is an alternative
    packaging, not an easier target.
* **Vojta-conditional.** Unchanged. It is still the only route from a standard conjecture. Nothing proves the
  needed instance (`Bl_Y X_7` or `Mbar_(0,8)`):
  * Chen-Garcia-Fritz-Mathur-Pasten is fixed-`S` and not shown to cover `Mbar_(0,n)`;
  * Huang and Lehmann-McKinnon-Satriano treat a fixed target point.

  `QS^prod` is non-standard and ranks below Vojta as a conditional route.
* **Route 5 (being stopped).** Nothing reopens it. No algebraic-independence measure for many logarithms exists in the
  sweep, and Calegari-Dimitrov-Tang and Ferreira-Ribas are single-form. Stopping route 5 loses nothing that this
  literature could have supplied.

**Recommended record changes** (the repository is read-only for me; to apply on merge):
1. In `round8/README.md` §2.3/§2.5, change `QS(m, delta, o(s))` from "open for `m >= 4`" to "false for every `m`".
   Prop A covers fixed `H_0`; Prop B covers any `H_0(s)`, and with Prop B the whole `QS` route caps at
   `M <= omega(N')`. Add `QS^prod` as open and non-standard.
2. Keep composite-supported analytic inputs **unclassified** in the barrier map.
3. Add the Gabdullin position correction to `local_smoothing/chains.md` row F: `gamma` corresponds to
   general-position arcs `R^(2gamma-1)`.
4. Check the Chan 2015 citations once the full text of arXiv:2507.07094 is read.

Further literature sweeps on this barrier are not recommended unless the full texts in §3, items 1-2, change a
classification.
