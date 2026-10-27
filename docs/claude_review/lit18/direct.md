# lit18 / direct: literature sweep on the direct problem (lattice points on arcs of length C sqrt(R))

Scope: work from 2015 to 2026 (mainly 2022-2026) on lattice points on short arcs of circles and conics; the
Cilleruelo-Granville (CG) and Cilleruelo-Cordoba (CC) conjectures; close lattice points on circles; squares in
arithmetic progressions; integer points on convex arcs; and points on `x^2+y^2=n` in small caps.

Method:
* WebSearch: about 60 queries, standard and extended.
* GitHub (reachable):
  * GitHub code search.
  * Raw downloads of `ehijano/rss_fetch@30c967b`, which mirrors the daily arXiv RSS listings. I took 964 daily
    listings from 2024-03-05 to 2026-10-27 for each of math.NT, CA, CO, AG, SP, MG, PR, AP, DS, GM and HO. That gives
    92,926 distinct arXiv records with titles and abstracts, which I searched with regular expressions.
  * The paper cards in `plasma-ai/erdos@05a6a79`, which quote CG 2007 and Chan 2015 page by page.
* arXiv, zbMATH and the other databases were unreachable, as stated.
* Nothing is cited that I did not see in a search result or a fetched file.
* Tools and raw data are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/lit18/`:
  `tools/`, `rss/`, `gh/`, and the query log `direct_log.txt`.

---

## 0. Bottom line

**Nothing found beats the barrier.**

* **Nobody has improved on `log R/loglog R`.** In everything I found (2015-2026), the best unconditional bound for
  arcs of length `C sqrt R` that holds for every `R` is Bourgain-Rudnick's `B_lambda << log lambda`. The project's own
  `(2/3+eps) log R/loglog R` is better than anything in the literature. No paper gives an `o(log R)` bound for every
  `R`, let alone `o(log R/loglog R)`.
* **No new uniform result at the `sqrt R` scale.** The only bounded results at or slightly beyond `sqrt R` are
  restricted by **position**, not by a class of `R`:
  * Chan (IJNT 2015): for `n = N^2` a perfect square, at most 10 points with `|b| < n^{1/4}(log n)^{1/7}`. These are
    arcs centred at the lattice point `(R,0)`.
  * Erdős-Rosenfeld, hyperbola analogue: divisors of `n` near `sqrt n` within `C n^{1/4}` number at most `1 + C^2`.

  Both are the "sagitta" argument at a direction of height `O(1)`. After a Gaussian rotation they do not reach
  general position: they miss by a factor `R^{1/2}`, computed in §2.3.
* **Results that hold for almost all `R` are not uniform.** These are:
  * Bourgain-Rudnick separation;
  * Granville-Wigman;
  * Kurlberg-Lester (Invent. 2025);
  * Sartori's spectral quasi-correlations, and Klurman-Sartori.

  They are capped by `round10/analytic.md` Thm B (statistical inputs), and they say nothing about the
  density-zero set of `N` that carries large clusters.
* **The one claimed advance on the CG conjecture was withdrawn.** Oganesyan, arXiv:2107.09991, claimed that arcs
  of length `R^alpha`, `1/2 < alpha < 1`, hold unboundedly many points, refuting CG Conjecture 14. Search snippets
  show the arXiv page marked "withdrawn".
* **The one proved consequence of CG beyond the `sqrt R` scale is averaged.** Gabdullin (arXiv:2205.13611, v3 Dec
  2023) proved the Cilleruelo-Cordoba `L^4` conjecture for squares in `[N, N+N^gamma]` for all
  `gamma < (sqrt5-1)/2 ≈ 0.618`. This is an energy statement. I show in §2.1 that it gives no information about a
  single cluster, so it is in the capped energy class. It is still useful evidence for the barrier diagnosis:
  averaged (fake-compatible) statements extend past `sqrt R`, while the pointwise statement is stuck below it.
* **No improvement on squares in arithmetic progressions.** Bombieri-Zannier's `k^{3/5+o(1)}` (2002) still stands.
  The 2024-2026 work is over number fields (González-Jiménez) or uses BZ as an input (Hegyvári 2025). In any case
  only arcs ⇒ squares is known (CG 2007), not the reverse.
* **New transcendence and unit-equation tools (2025-2026) are capped.** Calegari-Dimitrov-Tang (holonomy bounds),
  Ferreira-Ribas (arguments of unimodular algebraic numbers) and Ziegler (`<= (|S|-1)/2` solutions of
  `a_1x + a_2y = a_3`) are single-form or fixed-`S` statements. Ziegler's bound, even if it applied, is
  `(s-1)/2 <= (1+o(1)) log R/loglog R`, which is no better than what is already proved.
* **No restricted class of `R` with a new result.** The only `R`-restricted results are density-one, plus Chan's and
  Temur's special families near low-height directions.

An honest summary: the literature has not moved on the pointwise `sqrt R` problem since Bourgain-Rudnick and
Cilleruelo-Granville. The two papers worth reading in full are Gabdullin (to see whether its divisor lemma is a
pointwise statement with a Gaussian analogue) and Letendre (the rational analogue of the problem). Neither is
expected to break the barrier.

---

## 1. Results found

Confidence: **H** means the title, authors and venue were seen in several snippets or fetched files. **M** means
one snippet or a garbled snippet. Results are ordered roughly by relevance.

**R1. K. Oganesyan, "Lattice points on small arcs", arXiv:2107.09991 (2021), WITHDRAWN.**
* Snippet (Q1, Q4): "for any α∈(1/2,1) the number of lattice points belonging to an arc of length R^α of the
  circle of radius R centered at the origin is not uniformly bounded in R, which disproves the corresponding
  conjecture of Cilleruelo and Granville"; "estimates for L4-norms of Gauss sums"; "This paper has been withdrawn
  by Kristina Oganesyan".
* Another snippet describes arcs "of length φ(R) around the horizontal axis".
* Confidence: H for the existence of the paper and its withdrawal; the reason for the withdrawal is unknown.
* Relevance: the lower-bound side above `sqrt R` only. It has no bearing on upper bounds at `sqrt R`.

**R2. M. R. Gabdullin, "Trigonometric polynomials with frequencies in the set of squares and divisors in a short
interval", arXiv:2205.13611 (v1 May 2022, v3 Dec 2023); CANT 2022 talk, 27 May 2022.**
* Snippets (Q19, Q20): proves the Cilleruelo-Córdoba conjecture `||f||_4 <= c(gamma)||f||_2` for `f` with
  frequencies in `{n^2 : N <= n <= N + N^gamma}`, for every `gamma < (sqrt5-1)/2 = 0.618...`. The case
  `gamma <= 1/2` is easy (Sidon).
* The proof reduces to upper bounds for the number of divisors of an integer in a short interval. It also links
  the problem to Ruzsa's conjecture.
* Follow-ups: Gabdullin-Konyagin for cubes (arXiv:2311.14937); Garaev-Garayev-Konyagin, "On Sidon sets with
  squares, cubes and quartics in short intervals", arXiv:2602.08807 (2026).
* Confidence: H.
* Relevance: medium. It is the only theorem I found that proves a consequence of CG past the `sqrt R` scale. By
  CG 2007 Thm 16, Conjecture 13 implies the `L^4` Conjecture 19. But it is averaged.

**R3. P. Letendre, "Divisors of an integer in a short interval", arXiv:2503.12146 (15 Mar 2025).**
* Abstract (RSS corpus): studies `D_n(X,Y) = #{d | n : X <= d <= X+Y}`, `Y <= X`.
* Snippets (Q24, Q27): Thm 1 is `D_n(n^theta, n^eta) << tau(n)^{1-xi(theta,eta)} V(n) log tau(n)/(theta(1-theta))`
  for `0 < eta < theta < 1`; also `D_n(X,X) << tau(n)/sqrt(Omega_2(n))`; Conjecture 1 concerns intervals at
  `sqrt n`.
* Erdős-problems forum, #886, as quoted in a snippet (Q28): specialising to `theta = 1/2` gives `O(1/delta)` for
  the window `(n^{1/2}, n^{1/2} + n^{1/4-delta})`, and "for eps <= 1/4 the existence of a uniform bound remains
  open".
* Confidence: M (theorem shape from snippets).
* Relevance: the rational (hyperbola) analogue of our problem. See §2.2.

**R4. Erdős Problem #886 (Ruzsa) and Erdős-Rosenfeld (Acta Arith. 1997).**
* Snippets (Q21, Q26, Q28); GitHub issue `grandchallenge/MATH-PROGRAMME#2157`, fetched.
* Ruzsa asks whether the number of divisors of `n` in `(n^{1/2}, n^{1/2}+n^{1/2-eps})` is `O_eps(1)`. It is
  listed as open.
* Erdős-Rosenfeld: at most `1+C^2` divisors in `[n^{1/2}, n^{1/2}+C n^{1/4}]` for `n > n_0(C)`. Infinitely many
  `n` have 4 divisors in `(n^{1/2}, n^{1/2}+16 n^{1/4})`.
* A pull request on `google-deepmind/formal-conjectures` (#7218, Oct 2026) claims that the absolute-`K` variant
  (#887) is false with five divisors at `C = 31`. This is unrefereed and I only saw the snippet.
* Confidence: H.
* Relevance: the analogue at a direction of height 1 (§2.3).

**R5. T. H. Chan, "Factors of almost squares and lattice points on circles", Int. J. Number Theory 11 (2015)
1701-1708, arXiv:1406.2230.**
* Fetched card `plasma-ai/erdos .../chan_2015.../theorem_3.md`, quoting the preprint p. 2: "For sufficiently large
  perfect squares n = N^2, #{(a,b): a^2+b^2=n, |b| < n^{1/4}(log n)^{1/7}} <= 10."
* The proof uses simultaneous Pell equations and Turk's bound (1984).
* Also: perfect squares have at most five divisors near `sqrt n` (arXiv:1303.2069).
* Confidence: H.
* Relevance: a restricted position (§2.3).

**R6. F. Temur, "Discrete fractional integrals, lattice points on short arcs, and diophantine approximation",
arXiv:2012.10784; Turkish J. Math. 46(7) (2022) 2867-2889.**
* Snippets (Q34, Q50): bounds for discrete fractional integral operators along bivariate quadratic polynomials,
  with "connections to problems on concentration of lattice points on short arcs of conics".
* Its Conjecture 3: an arc of length `N^eta`, `eta < 1/2`, on the circle of radius `sqrt N` has `<= C_eta` points.
  The paper notes that this is known for `eta < 1/4` (CC).
* Its Theorem 5 uses simultaneous Pell equations viewed as hyperbolas with algebraic asymptote slopes, together
  with Schmidt's theorem. Chan uses Turk's effective bounds instead.
* Confidence: M on the theorem content.
* Relevance: special patterns only (§2.4).

**R7. P. Kurlberg, S. Lester, "Poisson spacing statistics for lattice points on circles", Invent. Math. 240
(2025) 537-585, arXiv:2112.08522.**
* Snippets (Q12):
  * Poisson nearest-neighbour spacing along a density-one subsequence of radii;
  * "at scales N^(-γ) with ½ < γ < log(π/2)/log 2, lattice points on circles do not behave like i.i.d. points";
  * the project's paper adds that some subsequences have spacing measure degenerating to a point mass at 0.
* Confidence: H. Relevance: low (§2.5).

**R8. A. Granville, I. Wigman, "Planck-scale mass equidistribution of toral Laplace eigenfunctions", CMP (2017),
arXiv:1612.07819.**
* Snippet (Q8): density-one subsequence; "rests on a number of results about the proximity of lattice points on
  circles, much of it based on foundational work of Javier Cilleruelo".
* Confidence: H. Relevance: low (density one).

**R9. Arithmetic random waves: number-theoretic inputs.**
* A. Sartori, spectral quasi-correlations, arXiv:2005.04698 ("for generic integers representable as a sum of two
  squares, there are no spectral quasi-correlations").
* O. Klurman, A. Sartori, arXiv:2104.13038 ("relies on a new study of lattice points in small arcs"; Lemma 1.5
  gives density-one separation).
* Benatar-Marinucci-Wigman, J. Anal. Math. 141 (2020) 707-749, arXiv:1710.06153 (uses bounds for spectral
  quasi-correlations).
* Cammarota-Klurman-Wigman, CMP 2020, doi 10.1007/s00220-020-03737-w ("spectral semi-correlations").
* Snippets: Q9, Q10, Q11, Q48. Confidence: H. Relevance: low (all generic or density one).

**R10. P. Kurlberg, I. Wigman, "On probability measures arising from lattice points on circles", Math. Ann.
(2016/17), arXiv:1501.01995.**
* Snippet (Q38): attainable measures; the prime-power ones are the extreme points of the symmetric measures.
* Follow-up: Sartori, Bull. Pol. Acad. Sci. Math. 66(2) (2018).
* Confidence: H. Relevance: none for arcs of length `sqrt R` (the scale is `O(1)` in angle).

**R11. E. Bombieri, J. Bourgain, "A problem on sums of two squares", IMRN 2015, no. 11, 3343-3407,
doi 10.1093/imrn/rnu005.**
* Snippet (Q40): the toy model comes from nodal curves on `T^2`; the tools are "sum-product theorems,
  diophantine approximation, and elliptic curves".
* Confidence: H on the citation, M on the content.
* Relevance: energy class; superseded by R12.

**R12. A. Cushman, C. Demeter, S. Wu, "Near optimal three-fold additive energy bound for points on convex
curves", arXiv:2608.12316 (Aug 2026, v3).**
* Abstract (RSS): `J_3(gamma(X)) <<_eps |X|^{3+eps}` for every finite `X` and strictly convex `gamma`. It "settle[s]
  a conjecture from Krishnapur-Kurlberg-Wigman and Bombieri-Bourgain concerning lattice points on dilates of the
  unit circle".
* Reproved by Z. K. Li and P.-L. Yung, "An improved high-low method for the parabola", arXiv:2609.33949.
* Confidence: H. Already in `local_smoothing/chains.md` row E. Capped (§2.6).

**R13. Bourgain-Rudnick (pre-2015, included for status).** arXiv:1012.3843 (AHP 2011), arXiv:1105.0018, 1308.4247.
* Snippets (Q3, Q9, Q31): "Lemma 2.1 gives O(log m) lattice points on arcs of length at most (√m)^{1/2}"; "the
  conjecture is known for almost all λ, but individually we only know B_λ << log λ".
* Separation: "for all but O(X^{1-ε/3}) integers m <= X, the minimal gap ... exceeds (√m)^{1-ε}".
* Confidence: H.

**R14. Squares in arithmetic progressions, 2015-2026.**
* No improvement over Bombieri-Zannier `N^{3/5+o(1)}` was found (Q14, Q15, Q55, Q56, and the RSS corpus).
* Related recent work:
  * E. González-Jiménez, arXiv:2602.03251 (quadratic extensions of number fields; genus-5 curve);
  * González-Jiménez, Nguyen Xuan Tho, arXiv:2602.01380;
  * N. Hegyvári, "Note on the sumset of squares", arXiv:2504.13230 (uses BZ as input:
    `|Q+Q| >= C|Q|(log|Q|)^{1/3+o(1)}`);
  * C. Elsholtz, L. Wurzinger, "Sumsets in the set of squares", Q. J. Math. 75(4) (2024) 1243-1254;
  * González-Jiménez-Xarles (Rudin's conjecture checked for `N <= 52`, arXiv:1301.5122).
* Confidence: H for "no improvement found". This is an absence of evidence, not proof.

**R15. V. Ziegler, "On the number of solutions to S-unit equations", arXiv:2609.15487 (Sep 2026).**
* Abstract (RSS): for fixed `a_1, a_2`, the equation `a_1x + a_2y = a_3` with `x, y` in `U_S` has at most
  `(|S|-1)/2` solutions "provided that the height of `a_3` is sufficiently large"; the result is effective.
* Confidence: H. Relevance: tested in §2.7.

**R16. Transcendence tools, 2025-2026.**
* G. C. G. Ferreira, S. Ribas, arXiv:2602.23597: if `beta` is algebraic with `|beta| = 1` and not a root of unity,
  then `Arg(beta)/2pi` is Diophantine (via linear forms in logarithms).
* F. Calegari, V. Dimitrov, Y. Tang, "Arithmetic holonomy bounds and effective Diophantine approximation",
  arXiv:2510.04156: a new effective method that reproves Baker-type theorems (two-variable S-unit equation,
  Thue-Mahler).
* Abstracts from the RSS corpus. Confidence: H. Capped (§2.8).

**R17. Determinant method and Bombieri-Pila, recent.**
* T. Bloom, J. Lichtman, "The Bombieri-Pila determinant method", arXiv:2312.12890 (exposition).
* "Lattice points on a curve via ℓ^2 decoupling", arXiv:2205.08206 (`|Γ ∩ (N^{-1}Z)^2| <= c(Γ,ε)N^{1/2+ε}` for a
  fixed curve `Γ`).
* S. Bhattacharjee, arXiv:2607.24048 (larger sieve with heights; a Bombieri-Pila-type bound over `O_K`).
* Binyamini-Cluckers-Novikov degree-uniform bounds, quoted in Bloom-Lichtman.
* Confidence: H. Capped (§2.9).

**R18. Not relevant, recorded so they are not chased again.**
* Guo-Ilyin, arXiv:2309.10971: annuli `lambda <= |x|^2 <= lambda + kappa` with `kappa >= C lambda^s`, `s < 1/4`, in
  which lattice points are far apart. This is sharp at `s = 1/4`. It is an existence result for many circles.
* Waxman-Yesha, arXiv:2305.02060, and Lester-Waxman-Yesha, arXiv:2609.23214: lattice points in thin sectors of the
  disc.
* Burq-Germain-Sorella-Zhu, arXiv:2507.16599, and Zhang-Zhu, arXiv:2606.08650: toral restriction, which uses the
  cluster structure as an input.
* Dimitrov-Gao-Habegger, Ann. of Math. 194 (2021) 237-298: uniform Mordell, with `#C(Q)` bounded in terms of
  `g`, `d` and the rank. The bound is exponential in the rank.
* Huang-Liu-Rudnick, arXiv:1903.04005: Gaussian primes in almost all narrow sectors.
* Letendre's other papers; Korsky, arXiv:2606.01780; Greenfeld-Iliopoulou-Peluse, "On integer distance sets",
  arXiv:2401.10821.
* The user's own `craigxchen.com/assets/Jarnik_CraigChen.pdf` also surfaced in searches.

---

## 2. Reductions attempted

Notation: `N = R^2`, the points are `z_1, ..., z_M` in `Z[i]` with `|z_j|^2 = N`, all lying in an arc of length
`L = C sqrt R`, so the angular width is `C R^{-1/2}`. `s` is the number of distinct split primes of `N`. Since
`s <= omega(N) <= (1+o(1)) log N/loglog N = (2+o(1)) log R/loglog R`, any bound linear in `s` is at best
`~ log R/loglog R`.

### 2.1 Gabdullin's L^4 theorem (R2): averaged, gives nothing per cluster

* **What it would use.** Quadruples `a^2+b^2 = c^2+d^2` with `a, b, c, d` in `[N_0, N_0+H]`, `H = N_0^gamma`. These
  are pairs of lattice points on a common circle of radius `R ≈ sqrt2 N_0` inside a box of side `H`. So the theorem
  controls arcs of length `~ R^gamma` with `gamma` up to `0.618 > 1/2`.
* **Reduction to a single cluster.** Take `M` cluster points `(a_j, b_j)` on one circle with all coordinates in the
  window. This is a cluster near the diagonal; CG's rotation trick moves any cluster there, at the cost of changing
  the scale, and the argument below is vacuous anyway. Let `A = {a_j} ∪ {b_j}`, so `|A| <= 2M`, and
  `f = sum_{n in A} e(n^2 x)`. Then
  * `||f||_4^4 >= #{(j,k)} = M^2`, counting the pairs `(a_j,b_j)`, `(a_k,b_k)`;
  * `||f||_2^4 = |A|^2 <= 4M^2`.

  So `||f||_4 <= c||f||_2` reads `M^2 <= 4c^4 M^2`, which is vacuous for every `M`. The `L^4` (Λ(4)) property
  cannot see a single circle.
* **Barrier.** This is the energy class: `README.md` §3 ("character, GAP-Fourier and energy arguments") and
  `local_smoothing/chains.md` row F. Fakes satisfy it.
* **What is still worth checking.** The proof reduces to "divisors of `m ~ X^{1+gamma}` in `[X, X+X^gamma]`". If
  that divisor lemma is pointwise, uniform in `m`, and has a Gaussian-integer analogue, then it would be a genuine
  input that fakes might violate. Writing `a-c = p`, `a+c = q`, `d-b = r`, `d+b = s` gives `pq = rs` with
  `q, s` in `[2N_0, 2N_0+2H]`. These are rational divisors of `m = pq <= 3N_0 H` in a window of length `2H`.
  * Two divisors `d_1 < d_2` of `m` force `gcd(d_1,d_2) >= d_1 d_2/m >= X/(1.5H)` and `gcd | d_2 - d_1 <= H`, so
    `H >= (X/1.5)^{1/2}`. This is why `gamma <= 1/2` is Sidon.
  * For `k` divisors, the lcm bound gives `X^{k-1} <= H^{k(k-1)/2+1}`, which is the CC-type exponent.

  Gabdullin must beat this on average over `m`, which is what reaches `0.618`. Pointwise, nothing past `1/2` is
  claimed in the snippets. **Request the full text** to see the exact lemma.

### 2.2 Letendre and Ruzsa (R3, R4): the rational analogue is just as stuck

* **Dictionary.** Our points are Gaussian divisors `z` of the rational integer `N` with `|z| = sqrt N`, whose
  arguments lie in a window of width `C N^{-1/4}`. The rational analogue is divisors of `n` in `[X, X+Y]`.
* **Where the critical scale sits.** At `X = n^theta`, the hyperbola `xy = n` has radius of curvature
  `rho ≈ (n/X)^3/(2n)`, so an arc of length `sqrt rho` corresponds to `Y ≈ X^{1/2}`. **The critical scale is
  `Y ~ X^{1/2}`.**
  * At `theta = 1/2` (the vertex) this is the window `n^{1/4}` of Erdős-Rosenfeld.
  * Letendre's specialisation gives `O(1/delta)` divisors in windows `n^{1/4-delta}`. That is exactly the CC
    analogue (`M(delta) ~ 1/delta`).
  * Above the critical scale it gives only `tau(n)^{1-xi}`.
* **Arithmetic.** Chopping a `C n^{1/4}` window into `C n^delta` pieces gives `<= C n^delta O(1/delta)`. With
  `delta = 1/log n` this is `O(log n)`, the classical log bound. A power `tau(n)^{1-xi}` is `exp(c log n/loglog n)`,
  far worse.
* **Conclusion.** No gain. The rational analogue in general position (`X = n^theta`, `theta < 1/2`,
  `Y = cX^{1/2}`) is open too: #886 is open and Letendre's Conjecture 1 is unproved. It may be a useful testbed,
  because fakes of the rational problem are easier to write down, but it supplies no input.

### 2.3 Chan, Erdős-Rosenfeld, and the near-axis results (R4, R5): restricted position only

* **What it uses.** The arc is centred at the lattice point `(R,0)` when `R` is an integer, or more generally near a
  direction `(p,q)` of small height `H = sqrt(p^2+q^2)`.
  * Rotate: `u = px + qy`, `v = -qx + py`, so `u^2 + v^2 = H^2 R^2`.
  * The arc condition gives `|v| <= C H sqrt R`, so `u` lies in an interval of length `<= v^2/(2u) <= C^2 H/2`.
  * Hence there are at most `2(1 + C^2 H/2)` points. This is exact integrality of `u`: a genuine lattice property,
    which fakes need not have.

  Chan pushes `H = 1` to a window `(log n)^{1/7}` beyond `sqrt R` by simultaneous Pell equations, but the count is 10
  and does not grow.
* **Transfer to general position.** Multiply by `\bar z_1`. The points `w_j = z_j \bar z_1` lie on the circle of
  radius `N = R^2`, of norm `n' = R^4`, at angular distance `<= C R^{-1/2}` from the axis. So
  `|Im w_j| <= C R^{3/2} = C n'^{3/8}`. Chan's window is `n'^{1/4}(log n')^{1/7} = R (4 log R)^{1/7}`.
  **The gap is a factor `R^{1/2}/(log R)^{1/7}`.** Using a height-`H` direction instead covers only arcs within
  angle `~ C R^{-1/2}` of a direction of height `H`, at a cost `O(C^2 H)`. By Dirichlet, `H` can be as large as
  `R^{1/4}`, which is Jarník-quality.
* **Conclusion.** No growth statement. These are the "near a rational direction of bounded height" clusters, which
  the record already treats as easy.

### 2.4 Temur (R6): special patterns via Schmidt

* **What it uses.** Simultaneous Pell equations, i.e. fixed hyperbolas with algebraic slopes, and Schmidt's
  Subspace/Roth-type theorem for a fixed template.
* **Constraint violated.** The threshold beyond which the pattern is excluded depends on the template and is
  ineffective. That is exactly "a template-dependent threshold treated as uniform in the radius" and "uniform
  exceptional sets from fixed-S Subspace Theorems", both excluded by the owner's constraints. Even for one template
  it says only "finitely often", with no rate in `M`. No growth bound.

### 2.5 Kurlberg-Lester, Granville-Wigman, Sartori, Bourgain-Rudnick separation (R7-R9, R13): density one

* **Scale mismatch (Kurlberg-Lester).** The number of points is `r(N) <= exp((log2+o(1)) log N/loglog N) = R^{o(1)}`,
  so the natural spacing scale is `2pi/r(N) >= R^{-o(1)}`. Our window `C R^{-1/2}` is smaller by a factor
  `R^{1/2-o(1)}`. Their degenerate (clustering) subsequences cluster at scale `o(1/r(N))`, which is still far
  coarser than `R^{-1/2}`. No content for our windows.
* **Bourgain-Rudnick separation.** For `N <= X` outside `O(X^{1-eps/3})` exceptions, the minimal gap is
  `>= R^{1-eps}`, which exceeds `C sqrt R`. So arcs of length `C sqrt R` hold `<= 1` point for those `N`. The
  extremal clusters (`log N ~ c M log M`) lie in the exceptional set.
* **Barrier.** These are statistical inputs, capped by `round10/analytic.md` Thm B (README §7.1 of round 10). The
  `L^4`/energy parts fall in the energy class.

### 2.6 Cushman-Demeter-Wu (R12): energy, cap computed

* Cluster sums `z_i + z_j + z_k` lie in a box of size `3C sqrt R` by `3 x sagitta`, where the sagitta is
  `L^2/(8R) = C^2/8`. The box holds `<~ ((9/8) C^3 + O(C)) sqrt R` lattice points.
* Cauchy-Schwarz gives `J_3 >= M^6/(#sums)`. Combined with `J_3 <= c_eps M^{3+eps}`:

  ```text
  M^{3-eps} <= c_eps ((9/8) C^3 + O(C)) sqrt R,   so   M << R^{1/6+eps}.
  ```
* The theorem holds for every finite set on a strictly convex curve, so fakes satisfy it. Capped
  (`local_smoothing/chains.md` row E).

### 2.7 Ziegler's S-unit count (R15)

* **Cluster formulation.** Let `v_jk = z_j \bar z_k`. It is supported on the primes of `N`, it is an `S`-integer
  with `S` the primes above `N` together with infinity, and `v_jk \bar v_jk = N^2`. It satisfies

  ```text
  v_jk + \bar v_jk = 2N - h_jk,   h_jk = |z_j - z_k|^2 <= C^2 R.
  ```

  This is `x + y = a_3` with `a_1 = a_2 = 1`, but `a_3 = 2N - h_jk` changes with the chord length. Ziegler counts
  solutions for a **fixed** `a_3`, so each value of `a_3` (each chord length) gets at most `(|S|-1)/2` pairs. That
  bounds chord-length multiplicities, not `M`. Normalising by `2N - h` destroys the S-unit property.
* **The threshold.** "Height of `a_3` sufficiently large" depends on `S`, hence on `N`. That is a template- or
  `S`-dependent threshold, excluded by the constraints.
* **Even if it applied.** `(s-1)/2 <= (1+o(1)) log R/loglog R`, which exceeds the proved `(2/3+eps) log R/loglog R`.
  No growth gain.

### 2.8 Ferreira-Ribas and Calegari-Dimitrov-Tang (R16)

Both are single-form transcendence measures, of linear forms in logarithms or effective irrationality.

* For `beta = z_j/z_k`, a measure `|Arg beta| >= exp(-c h(beta) ...)` with `h(beta) ~ (1/2) log N` gives
  `|theta_j - theta_k| >= N^{-c'}`. Here `c'` is a large Baker-type constant, while the cluster needs only
  `C R^{-1/2} = C N^{-1/4}`. This is useless unless `c' < 1/4`.
* In any case, single-form Lang-Waldschmidt is capped even with the optimal constant: `round8/README.md` table row
  "Lang-Waldschmidt, single forms (product or max form, any eps >= 0)", `round8/lw.md` Thms A-B.
* Baker/Matveev: `round8/stability.md` Prop 6.1.
* The round-16 judge's remark (`round16/README.md` §7.4) asks for an algebraic-independence measure for about `n`
  logarithms. Neither paper claims one: CDT give irrationality and linear-independence measures and reprove
  two-variable S-unit results, which are exponential in `|S|`.

### 2.9 Bombieri-Pila and the determinant method (R17)

* On an arc of length `C sqrt R` the determinant method gives only powers of `R`, and Jarník's convexity bound gives
  `L^{2/3}`. Degree-uniform versions (BCN, CCDN) and larger-sieve versions only improve the dependence on the degree,
  which is 2 here. The 2205.08206 decoupling bound is for a fixed curve with shrinking lattice spacing, a different
  regime.
* Interpolation determinants on one circle reduce to Vandermonde times pair data (`growth/auxiliary.md`). The
  single-place ratio is capped at `tau* < 2` (`threshold/README.md`). Capped.

### 2.10 Oganesyan (R1) and squares in AP (R14)

* **Oganesyan.** Even if it were correct, it concerns `alpha > 1/2` and gives lower bounds, so it changes no upper
  bound at `sqrt R`. Through the CG equivalence (Conj 13 ⇔ 14, CG 2007 p. 11) it would refute CG's `R^{1-eps}`
  conjecture (Conjecture 14), but not the `C sqrt R` statement this program targets. It is withdrawn.
* **Squares in AP.** CG prove only the direction arcs ⇒ squares-in-AP and Rudin-type statements. An improvement
  of BZ would not transfer. The BGP/BZ method counts rational points on genus-5 curves. Uniform Mordell (DGH)
  bounds are `c(g)^{1+rank}`. In any transfer, the rank grows with the number of primes, which is the forbidden
  `c^s` regime.

---

## 3. Papers to request in full text (ranked)

1. **Gabdullin, arXiv:2205.13611 (v3, Dec 2023).** I need the exact "divisors in a short interval" lemma: is it
   pointwise in `m`, what are its parameters, and does it have a Gaussian analogue? It is the only proved
   statement past `sqrt R` connected to CG. Read §§2-3 to separate the averaged step from the pointwise step. If
   any step is a pointwise bound uniform in `m` past the Sidon range, it is a candidate genuine input; I expect
   it to be averaged.
2. **Letendre, arXiv:2503.12146.** I need Thm 1 (`xi(theta, eta)`, `V(n)`), Prop 1 (`Y ~ n^{theta^2}`) and
   Conjecture 1. This decides whether the rational analogue at its critical scale `Y ~ X^{1/2}` has anything
   better than `tau(n)^{1-xi}`, and is the sanity check for transferring the program to hyperbolas.
3. **Oganesyan, arXiv:2107.09991 (v1 and the withdrawal note).** What construction was claimed near the axis, and
   what failed? This matters only for the lower-bound side and for whether the `L^4` Gauss-sum estimates survive.
4. **Ziegler, arXiv:2609.15487.** How does the height threshold depend on `S`, and is there a version with the
   right-hand side varying in a bounded-height family? This would only matter if a fixed-`a_3` formulation of
   clusters existed, and I found none.
5. **Calegari-Dimitrov-Tang, arXiv:2510.04156.** Do the holonomy bounds give measures of linear independence for
   `k` logarithms with dependence on `k` better than `c^k`? This is relevant only to the round-16 remark about
   algebraic independence of many logarithms. Low expectation.
6. **Temur, Turkish J. Math. 46 (2022) 2867-2889.** The exact Theorem 5 and its pattern class, to confirm that it
   is template-dependent. Low.
7. **Bombieri-Bourgain, IMRN 2015.** To check whether any pointwise arc statement is hidden in the elliptic-curve
   part. Low.

---

## 4. Search log

The full log with snippets is in
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/lit18/direct_log.txt` (Q1-Q60
plus the GitHub and RSS-corpus searches). Summary of the queries:

* **Direct problem.**
  * "lattice points on short arcs of circles Cilleruelo Granville conjecture" [ext]
  * "lattice points" "short arcs" circle sqrt R bound 2024
  * "Cilleruelo-Córdoba" conjecture arcs lattice points proof
  * Oganesyan "Lattice points on small arcs" [ext]
  * "arcs containing" lattice points circle "k lattice points" Cilleruelo Zygmund
  * lattice points on circles close together sqrt R arc uniform bound 2025 [ext]
  * "Cilleruelo" "Granville" conjecture "short arcs" 2024-2026 [ext]
  * "log R / log log R" lattice points arc [ext]
  * "four/five lattice points" short arc R^{2/5}
  * "close lattice points" follow-up Pell
  * "lattice points" arc "uniformly bounded" 2022-2023
  * lattice points arc R^{1/2} uniform recent [ext]
  * Cilleruelo "Arcs containing no three lattice points" follow-up
  * erdosproblems circle lattice points arc (twice)
  * abc / Vojta / Bombieri-Lang implies arcs bounded (two queries)
* **Close points and toral eigenfunctions.**
  * Granville-Wigman Planck-scale
  * Sartori spectral correlations Bourgain-Rudnick separation
  * Benatar-Marinucci-Wigman
  * Klurman-Sartori [ext]
  * Kurlberg-Lester Poisson
  * small-scale angular distribution 2023-2025 [ext]
  * Kurlberg-Wigman attainable measures
  * Bombieri-Bourgain sums of two squares (two queries)
  * Bourgain-Rudnick restriction caps
  * Zhang-Zhu restriction
  * spherical caps 2023-2025
  * Cammarota-Klurman-Wigman
  * Kátai-Környei / Erdős-Hall
  * minimal gap lattice points circle
* **Divisors and the hyperbola.**
  * Gabdullin L4 squares (two queries)
  * Gabdullin divisors short interval
  * Letendre (three queries)
  * Erdős-Rosenfeld [ext]
  * erdosproblems 886 [ext] (twice)
  * Tsz Ho Chan 2020-2023
* **Squares in arithmetic progressions.**
  * BZ k^{3/5} improvement [ext]
  * Rudin squares in AP 2023-2025
  * uniform Mordell DGH squares in AP
  * "Note on the sumset of squares"
  * "Sumsets in the set of squares"
  * squares in AP N^{3/5} improved [ext]
  * BGP genus 5 uniform Mordell
* **Convex arcs and the determinant method.**
  * integer points on strictly convex arcs Jarník 2022-2024 [ext]
  * Bombieri-Pila uniform BCN
  * Huxley / Swinnerton-Dyer convex arcs new
  * integer points short arcs of conics (Ramana, Shparlinski)
* **Transcendence and analytic inputs.**
  * linear forms in logarithms Gaussian lattice points
  * Gaussian primes narrow sectors 2023-2025
  * Gaussian divisors argument short interval
  * Temur (three queries)
* **GitHub.**
  * code search "Cilleruelo Granville lattice points arc" (50 hits): plasma-ai/erdos, formal-conjectures,
    rss_fetch.
  * code search in `ehijano/rss_fetch` for "short arcs", "Cilleruelo", "Jarnik", "lattice points arc circle".
  * raw downloads: the RSS corpus (11 categories × 964 days) and the plasma-ai/erdos cards for CG 2007 and
    Chan 2015.
* **RSS-corpus regexes** (92,926 abstracts, 2024-03..2026-10):
  * `lattice points` × arc/circle;
  * Cilleruelo; Jarník; squares in arithmetic progression; Rudin; Bombieri-Zannier; Erdős-Rosenfeld;
  * Ruzsa / divisor short interval; sums of two squares; Gaussian integers/primes;
  * toral / arithmetic random waves / nodal; Bombieri-Pila / determinant method; integer points on curves;
  * concyclic; linear forms in logarithms; Vojta; unit equations; subspace theorem;
  * small/short arcs; close lattice points; angles of Gaussian/lattice points;
  * additive energy circle/convex; spherical caps; sectors; Ramana/Vandermonde;
  * squares in short intervals; Sidon squares; Granville; Kurlberg/Lester.

  No 2024-2026 arXiv abstract in these categories concerns lattice points **on** short arcs of a circle except CDW,
  which is an energy statement.
