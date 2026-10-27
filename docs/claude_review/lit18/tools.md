# lit18 / tools: height-based Diophantine tools that are uniform in S (2018-2026 sweep)

Scope: tools that control many algebraic points of bounded height at once, uniformly in the prime set `S`, with
bounds polynomial in `|S|`. These are gcd theorems for `S`-units, Vojta for blowups, Heier-Levin, Ru-Vojta/Vojta beta
constants, Seshadri constants, uniform dimension growth, effective Pila-Wilkie, near-1 counts in finitely generated
groups, and unit equations. For each tool the report checks, with the arithmetic done, whether it can give
`M = o(log R/loglog R)`.

Network: arXiv, ar5iv, par.nsf.gov, Leiden, Annals and unistra all failed with DNS errors through WebFetch, so every
citation below rests on WebSearch snippets. The query log is in §4. Exact scripts and their outputs are in
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/lit18/checks/`.

Conventions. `N = R^2`, `W = log N`, `s = omega(N') <= (1+o(1)) W/log W` with `N'` the normalised norm. The current
bound is `M <= (2/3+eps) log R/loglog R ~ (1/3) W/log W ~ s_max/3`. A bound of the form `M <= f(s) + O(1)`
**improves the growth rate iff `f(s) = o(s)`**. `f(s) <= s/3` only matches the constant, and any `f(s) >= s` is useless.

---

## 0. Bottom line

**Nothing found beats the barrier.** No paper from 2018-2026 that I could find gives any height-based statement that is
uniform in `S`, with a polynomial or sublinear count, in a form that applies to our clusters. Each tool family ends
in one of three places:

1. **Fixed-`S` theorems** (Corvaja-Zannier, Levin 2019, Grieve-Wang 2020, Wang-Yasufuku 2021, Huang-Levin, Xiao,
   Rousseau-Turchet-Wang 2024, Heier-Levin and Ru-Vojta with `S` containing the primes of `N`). Their exceptional sets
   depend on `S`, which the constraints forbid us to ignore. Their quantitative forms are exponential in `s`.
   **New (§2.2, Prop A, unrefereed but essentially rigorous):** exponential-type `S`-dependence is *forced* for any
   statement of the standard shape, meaning one valid for all heights above a threshold fixed independently of `S`.
   Dense smooth `S`-integer tuples solve the round-8 system (5.1) with every `delta < (m-1)/2` and need
   `>= exp(s^(1/2-o(1)))` subspaces. So:
   * a "polynomial in `|S|`" version of any of these theorems is false;
   * round 8's hypothesis `QS(m, delta, f, H_0)` is **false for every `m`**, not only for `m <= 3`;
   * the only usable form is a **sparse-regime** statement, `QS^sparse`: heights `H >= exp(A s log s)`. The
     literature does not treat this regime, and at `m = 2` it is close to a restatement of the problem.
2. **Single-place (`S = {inf}`) theorems** (Heier-Levin with Seshadri constants, Ru-Vojta/Ru-Wang/Vojta beta
   constants, Yasufuku 2025, Barrios 2025). These are capped by `threshold/upper.md` (`tau*(M) = 2 - 1/ceil(M/2) < 2`).
   Heier-Levin is even weaker: `eps_Y(H) <= 1` on `X_M`, because lines through every point of `Y` lie in `X_M`. That
   makes it short by a factor `>= 2(M+1)`. Barrios's corrected theorem needs an ample normal bundle, but
   `deg N_{Y/X_M} = 0`, so it does not apply. My searches found no computation of a Ru-Vojta/Nevanlinna constant
   above 2 for such configurations.
3. **Counting tools** (uniform dimension growth, Bombieri-Pila, effective Pila-Wilkie and Wilkie, unit-equation and
   uniform Mordell-Lang bounds).
   * The determinant-method family reduces to pair data (`growth/auxiliary.md` Prop A) and its counts are trivial at
     `d = 2`.
   * Pila-Wilkie needs a transcendental part, and the circle is algebraic, so that part is empty. The transcendence
     sits in linear forms in logarithms, which are capped (`round8/README.md` §2.4).
   * Unit-equation and Mordell-Lang counts are provably super-polynomial in the rank: the
     Erdős-Stewart-Tijdeman / Konyagin-Soundararajan / Granville / Harper lower bounds.

Honest summary: the property that actual clusters have and fakes lack is that the points are `S`-integers whose
heights are tied to `S` (product formula, sparse regime). Subspace-type tools do use it, but every proved theorem of
this kind is uniform over heights above a fixed threshold. In that formulation, the dense regime (smooth `S`-integers)
forces super-polynomial `S`-dependence. A useful statement would have to be specific to the sparse regime, and no
such tool exists.

---

## 1. Results found

Confidence key:
* **A**: bibliographic data and main statement confirmed by snippets.
* **B**: bibliographic data confirmed, main statement only partly confirmed.
* **C**: statement second-hand, from a citing paper's snippet or from the repo.

| # | Citation (as confirmed) | Statement as confirmed by snippet | Conf. | Relevance |
|---|---|---|---|---|
| 1 | A. Levin, *Greatest common divisors and Vojta's conjecture for blowups of algebraic tori*, Invent. Math. (2019), pp. 493-533, DOI 10.1007/s00222-018-0831-z | Extends gcd bounds for polynomial values at `S`-unit points to any number of variables. Proves special cases of Vojta's conjecture for blowups of toric varieties. Main tool: Schmidt's Subspace Theorem. The form quoted by Xiao: `log gcd(f(u),g(u)) < eps max h(u_i)` outside a Zariski-closed set | A/B | Fixed `S`. §2.2, §2.3 |
| 2 | Corvaja-Zannier (Monatsh. Math. 144, 2004/5) and Hernández-Luca: `gcd(u-1,v-1) < max(|u|,|v|)^eps` for multiplicatively independent `S`-units, except on finitely many proper subtori (`u^p = v^q`, `|p|,|q| <= 1/eps`) plus a finite set. Bugeaud-Corvaja-Zannier, Math. Z. 243 (2003) | as stated (snippets from Levin's notes, Pasten-Wang and Yasufuku) | A | Base case. No `S`-uniform count found |
| 3 | N. Grieve, J. T.-Y. Wang, *Greatest common divisors with moving targets and consequences for linear recurrence sequences*, Trans. AMS 373 (2020) 8095-8126 (arXiv 1902.09109) | Bounds the log-gcd of moving multivariable polynomials at moving `S`-unit arguments | A | Moving targets, still fixed `S` |
| 4 | A. Levin, J. T.-Y. Wang, gcd and Nevanlinna theory on algebraic tori, J. reine angew. Math. (online 2019) (arXiv 1903.03876) | Upper bounds for the counting function of common zeros of meromorphic functions (complex analogue) | A | Analogue only |
| 5 | J. T.-Y. Wang, Y. Yasufuku, *Greatest common divisors of integral points of numerically equivalent divisors*, Algebra & Number Theory 15 (2021) 287-305 (arXiv 1907.09324) | Height of a codim-`>= 2` subscheme is small at integral points when the divisor is a sum of `n+1` effective divisors numerically equivalent to multiples of an ample divisor. Uses Ru-Vojta | A | Fixed `S` |
| 6 | K. Huang, A. Levin, *Greatest common divisors on the complement of numerically parallel divisors* (arXiv 2207.14432). Follow-ups: Huang-Levin-Xiao (arXiv 2406.18879), Wang-Xiao (arXiv 2506.03534) | gcd inequalities at integral points for numerically parallel divisors. Reduces to `G_m^n` via Vojta's theorem on semiabelian varieties. Also studies the exceptional set | A | Fixed `S` |
| 7 | Z. Xiao, *Greatest common divisors for polynomials in almost units and applications to linear recurrence sequences* (arXiv 2110.01751); MSU thesis | Generalises Levin to "almost `S`-units". The exceptional set is "within the control of a logarithmic region", which is necessary | B | Fixed `S` |
| 8 | Y. Yasufuku, *GCD inequalities inspired by Vojta's conjecture*, Monatsh. Math. 206 (2025) 259-279 | gcd inequalities for *all* algebraic numbers. Bounds = height contributions + local heights outside `S`. Main ingredient: Ru-Vojta | A | `S = {inf}` version, capped (§2.4) |
| 8' | Y. Yasufuku, *GCD inequalities arising from codimension-2 blowups*, Bull. LMS (2026), Wiley DOI 10.1112/blms.70349 | Title, journal and year from a link listing only. Content not confirmed | C | Request |
| 9 | B. Barrios, *An upper bound for the generalized greatest common divisor of rational points*, Bull. Aust. Math. Soc. (2025) (arXiv 2411.06050), plus Correction | Upper bound for `h_Y(x)` via a Riemann-Roch-type inequality. The correction says the results hold only if `Y` is positively embedded (normal bundle ample); Matsuzawa gave a blow-up counterexample | A | `S`-free but capped (§2.4) |
| 10 | G. Heier, A. Levin, *A generalized Schmidt subspace theorem for closed subschemes*, Amer. J. Math. 143 (2021) 213-226 (arXiv 1712.02456) | Subspace theorem for closed subschemes in general position, weighted by Seshadri constants with respect to a fixed ample divisor. Recovers McKinnon-Roth | A | §2.4 |
| 10' | Heier-Levin, *A Schmidt-Nochka theorem for closed subschemes in subgeneral position* (arXiv 2308.11460; DOI 10.1515/crelle-2024-0085). He-Ru (arXiv 1910.07966) | Subgeneral-position versions. Uses lower bounds for Seshadri constants of intersections | A | Same cap |
| 11 | M. Ru, P. Vojta, *A birational Nevanlinna constant and its consequences*, Amer. J. Math. 142 (2020) 957-991. Also Ru-Vojta, *An Evertse-Ferretti Nevanlinna constant* (arXiv 2004.12257); Ru-Wang, *The Ru-Vojta result for subvarieties*, IJNT 18 (2022) 61-74; P. Vojta, *Birational Nevanlinna constants, beta constants, and Diophantine approximation to closed subschemes*, JTNB 35 (2023) 17-61 | `Nev_bir(D)` computed via Autissier filtrations. The closed-subscheme versions. The several notions of `beta(L,D)` coincide and are limits | A | Capped: `beta_M <= tau* < 2` |
| 12 | E. Rousseau, A. Turchet, J. T.-Y. Wang, *Divisibility of polynomials and degeneracy of integral points*, Math. Ann. 388 (2024) 1969-1999 (arXiv 2106.11337) | Divisibility framework (Corvaja-Zannier). Ru-Vojta with sharper control of the exceptional set. Simply connected varieties with non-dense integral points. Special cases of Vojta for toric blowups (attribution in the snippet is ambiguous between this paper and Levin's) | B | Fixed `S` |
| 13 | J.-H. Evertse, R. Ferretti, *A further improvement of the Quantitative Subspace Theorem*, Ann. of Math. 177 (2013) 513-590 (arXiv 1008.2340), and the survey *On the quantitative Subspace Theorem* (arXiv 1008.2268) | `n`-dependence improved from `c^(n^2)` to `c^n` and `delta^(-3)(log delta^(-1))^2`. Snippet: "R may be replaced by ns, where s is the cardinality of S ... in many cases R can be taken independently of s". Evertse-Schlickewei's bound is "about `(9n^2 eps^(-1))^(ns)` times" the EF quantity. Repo quote (`round8/weakest.md` §5.2): Evertse 1996, `(2^(60 m^2) delta^(-7m))^s` | B/C | Exact `s`-dependence unresolved; **request** |
| 13' | N. Grieve, *On qualitative aspects of the quantitative subspace theorem* (arXiv 2306.16583) | Parametric Subspace Theorem in the Evertse-Ferretti form; partitions the exceptional set | B | Background |
| 14 | Unit equations. Evertse-Schlickewei-Schmidt, Ann. of Math. 155 (2002) 807-836: bound depending only on `n` and the rank `r`. Amoroso-Viada (Duke Math. J. 150 (2009)): `(8n)^(4n^4(n+r+1))` non-degenerate solutions (quoted by citing papers). Beukers-Schlickewei, Acta Arith. 78 (1996): `2^(8r+8)` for `x+y=1`. Evertse: `3 x 7^(2s+1)` | as stated | A/C | Exponential in rank |
| 14' | Lower bounds. Erdős-Stewart-Tijdeman: `exp((4+o(1))(s/log s)^(1/2))`. Konyagin-Soundararajan (arXiv math/0604453): `exp(s^(2-sqrt2-eps))` (exponent from citing snippets). Granville (via arXiv math/0107219): `exp(c_0 s^(1-1/n)(log s)^(-1/n))` for `n`-term equations. Harper (arXiv 1108.3819): `e^(s^(1-1/(n-1)-o(1)))` | as stated | A/C | Polynomial-in-rank counts are **false** |
| 15 | V. Dimitrov, Z. Gao, P. Habegger, *Uniformity in Mordell-Lang for curves*, Ann. of Math. 194 (2021) 237-298 | Number of points bounded in terms of `g`, `d` and the Mordell-Weil rank. Exact rank dependence not shown in snippets | B | Not applicable (genus 0) |
| 16 | P. Salberger, *Counting rational points on projective varieties*, Proc. LMS 126 (2023) 1092-1133. Castryck-Cluckers-Dittmann-Nguyen (arXiv 1904.13109). Binyamini-Cluckers-Novikov, IMRN (2024) (arXiv 2302.04209): `<< d^2 H^(2/d) (log H)^kappa` on plane curves, optimal in `d`. Binyamini-Cluckers-Kato (arXiv 2401.03982). M. Walsh, IMRN 2015, 5644-5658: `O(N^(2/d))`. Paredes-Sasyk, ANT 16 (2022) 1941 | Uniform dimension growth, polynomial in the degree | A | Capped (§2.6) |
| 17 | G. Binyamini, D. Novikov, *Wilkie's conjecture for restricted elementary functions*, Ann. of Math. 186 (2017) 237-275. Binyamini-Novikov-Zak, *Wilkie's conjecture for Pfaffian structures*, Ann. of Math. 199 (2024) 795-821. Binyamini-Jones-Schmidt-Thomas, effective Pila-Wilkie for Pfaffian sets, JEMS (2026) DOI 10.4171/JEMS/1761 (arXiv 2301.09883). G. Binyamini, *Point counting for foliations over number fields*, Forum Math. Pi 10 (2022) e6 | Polylog-in-`H` counts on transcendental parts. Effective in Pfaffian complexity | A | Inapplicable (§2.7) |
| 18 | Z. Huang, *Rational approximations on toric varieties*, Algebra & Number Theory 15 (2021) 461-512 (arXiv 1809.02001) | Under a geometric hypothesis, best McKinnon-Roth approximations to a *fixed point* of a split toric variety lie on minimal rational curves. Confirms a McKinnon conjecture | A | Fixed target point only (§2.8) |
| 19 | F. Temur, *Discrete fractional integrals, lattice points on short arcs, and diophantine approximation*, Turkish J. Math. 46 (2022) 2867-2889 (arXiv 2012.10784) | Uses Schmidt's theorem for fixed gap patterns `m_0 - m_i = h_i` on circles. "Schmidt's theorem is ineffective, which forces us to fix h_i beforehand ... does not lead to any new results on Conjecture 1"; "if we knew the constant ... of appropriate size, we would obtain Conjecture 1 for some eta > 1/4" | A | Template-dependent (§2.9) |
| 20 | K. Oganesyan, *Lattice points on small arcs* (arXiv 2107.09991) | Claimed unboundedness on arcs `R^alpha`, `alpha in (1/2,1)`. **Withdrawn** (mistake in Lemma 2) | A | None |
| 21 | Bugeaud-Evertse-Győry, *S-parts of values of univariate polynomials, binary forms and decomposable forms at integral points*, Acta Arith. 184 (2018) 151-185 | `S`-part bounds. Effective exponent explicit in `S`. Tools: Baker, `p`-adic Subspace, Barroero-Widmer | B | Fixed `S` |
| 22 | R. C. Chen, N. Garcia-Fritz, S. Mathur, H. Pasten, *Towards Lang-Vojta via degeneration* (arXiv 2602.06956, Feb 2026) | Finiteness and Zariski degeneracy for `S`-integral points | B | Not shown to cover `Mbar_{0,n}`; fixed `S` |

Searched for and **not found**:
* Vojta proved for `Mbar_{0,n}` (only the expectation, arXiv math/0301272) or for blowups of products of conics;
* any `S`-uniform gcd or Subspace count polynomial in `s`;
* any Ru-Vojta/Nevanlinna constant above 2 for point configurations on products of conics;
* any application of Pila-Wilkie to lattice points on circles.

---

## 2. Analysis: which property each tool uses, or which barrier caps it, with the arithmetic

### 2.1 The common reduction

Take a cluster `z_1..z_M` on `|z|^2 = N` within an arc of length `C sqrt R`. Put `u_i = z_i/z_1 in Gamma_S`, where `S`
is the set of Gaussian primes above `N`. Then:
* `|u_i| = 1` and `h(u_i) = (1/2) log Norm(z_1/g) <= (1/2) log N = log R`;
* `|u_i - 1| <= C R^(-1/2)`, so the archimedean proximity to the point `(1,...,1)` is `lambda_inf >= (1/2) log R - log C`;
* the ratio `lambda_inf / h >= 1/2 - o(1)`.

With the projective height `h(1:u_2:...:u_m) = log R - (1/2) log Norm g`, the ratio is still `>= 1/2` **for every
`m`**. This is the only "gcd" information a cluster carries at `inf`. The finite gcd of all differences is factored
out by `g`, so `sum_w lambda_w ~ 0` for full-profile clusters.

In the language of round 8 (`round8/weakest.md` §5.1, Lemma 5.1), the `m`-tuples solve (5.1) for every
`delta < (m-1)/2`. Theorem 5.2 then gives `M <= f(s) + m - 1` if (5.1) has `f(s)` exceptional subspaces with
`P = primes(N')`.

### 2.2 Fixed-`S` gcd/Subspace/Vojta-blowup theorems (#1-7, 10-13, 21 with `S` containing the primes of `N`)

*Property used.* These theorems are fake-separating in principle: they apply only to genuine `S`-integral algebraic
points, through the product formula (finite part exactly `H^(-m)`). *Barrier.* The exceptional set depends on `S`
(constraint). The quantitative versions give `f(s) = exp(O(s))` (`round8/weakest.md` §5.2), which is useless. The
question for this sweep was whether a 2018-2026 result improves the `s`-dependence to polynomial. None does, and the
following shows that none can, in the standard formulation.

**Proposition A (new; unrefereed; proof sketch complete).** Fix `m >= 2`, `0 < delta < (m-1)/2`, `H_0`, and `A > 2`.
For `R -> inf` let `P = (log R)^A`, let `P_s` be the split primes `<= P` (so `s ~ P/(2 log P)`), and let `S` be the
places above `P_s`. Then the solutions of (5.1) with `H(x) >= H_0` are not contained in fewer than `c_A R^(1-2/A)`
hyperplanes. Equivalently, `QS(m, delta, f, H_0)` forces `f(s) >= exp(c_A (s log s)^(1/A))`, and so
`f(s) >= exp(s^(1/2-o(1)))` on taking `A -> 2+`.

*Proof.*
1. **Many `S`-integers.** Count `S`-integers (products of the `2s` Gaussian primes, up to units) of norm `<= X = R^2`.
   They are at least `#{b in Z_{>=0}^(2s) : sum b_k <= K}` with `K = floor(log X/log P)`. That count is
   `C(K+2s, K) >= (2s/K)^K = exp(K (A-1) loglog R (1+o(1))) = R^(2-2/A-o(1))`.
2. **A dense annulus.** By dyadic pigeonhole, a shell `R'/2 < |w| <= R'` with `R' >= R^(1-o(1))` contains
   `U >= R'^(2-2/A-o(1))` of them.
3. **Many tuples.** Cut the shell into about `4R'` squares of side `sqrt(R')/2`. The mean occupancy is
   `mu = U/(4R') >= R'^(1-2/A-o(1))`.
4. **Each tuple solves (5.1).** Take `m` distinct `S`-integers `w_i` in one square, with `g = gcd`, `x = w/g`, and
   `|g| <= sqrt(R')`. The finite part is `prod_i Norm(x_i)^(-1) <= 4^m H^(-m)`. The archimedean part is
   scale-invariant and `<= 4^m (R'^2)^(-(m-1)/2) <= 4^m H^(-(m-1)/2)`. So
   `LHS <= 16^m H^(-m-(m-1)/2) <= H^(-m-delta)`, since `H >= R'/4` is large.
5. **Bad tuples are few.** Tuples with `|g| > sqrt(R')` have primitive parts of size `< sqrt(R')` that lie within `C`
   of each other. There are at most `U R^(o(1)) C^(2m)` of them, which is negligible against `R' mu^m`.
6. **Distinct points.** Two good tuples in the same square give the same projective point only if `w' = lambda w`
   with `|g' - g| <= 1`, so the multiplicity is `<= 5`.
7. **Count hyperplanes.** A hyperplane with `c_j != 0` contains at most as many points as there are
   `(m-1)`-subtuples. By Chebyshev's sum inequality,
   `sum_Q (mu_Q)_m / sum_Q (mu_Q)_(m-1) >= mu_bar - m + 1`. Hence `#hyperplanes >= (mu_bar - m)/5`.
8. **In terms of `s`.** Finally `log R ~ (2A s log s)^(1/A)`. QED.

*Exact and numerical evidence* (`checks/check_smooth2.py`, `checks/count_mu.py`; outputs in `checks/out_C.txt` and
`checks/out_mu.txt`):
* For `m = 4`, 3000 sampled tuples of smooth Gaussian `S`-integers in common squares solve (5.1) with
  `delta_eff in [1.61, 2.62]`, so they cover the whole window `(1, 3/2)` that round 8 left open.
* The projective multiplicity is `1.00`.
* The mean occupancy per square already exceeds `s`: `mu_bar = 68` at `s = 51` (`P = 600`, `R = 10^5`) and
  `mu_bar = 137.6` at `s = 80` (`P = 1000`, `R = 10^5`). So at least about `s` hyperplanes are needed at these finite
  sizes, and `exp(s^(1/2-o(1)))` asymptotically.

**Consequences.**
1. *No polynomial-in-`|S|` quantitative Subspace, gcd or Vojta-blowup theorem exists in the standard
   (fixed-threshold) form.* Specialised to (5.1), or to the gcd of `u_2-1, ..., u_m-1` at `(1,...,1)` in
   `G_m^(m-1)`, any such theorem must allow `>= exp(s^(1/2-o(1)))` exceptional components. Together with the known
   upper bounds `exp(O(s))`, the truth for fixed-threshold statements is between `exp(s^(1/2-o(1)))` and `exp(O(s))`.
   This is the Subspace analogue of the Erdős-Stewart-Tijdeman / Konyagin-Soundararajan lower bounds for unit
   equations (#14').
2. *Barrier map correction.* `round8/README.md` §2.3 lists `QS(m >= 4, delta in (1,(m-1)/2), o(s))` as "neither
   refuted nor supported". As literally stated (fixed `H_0`, all `P`), **it is false for every `m >= 2`**.
3. *The repaired hypothesis.* Theorem 5.2 only ever uses `P = primes(N')` and heights
   `H(x) = N/|g|^2 >= R/C^2 = N^(1/2)/C^2` (`round8/weakest.md` Lemma 5.1). Here
   `N' >= prod_(p|N') p >= exp((1+o(1)) s log s)`, so `H(x) >= exp((1/2+o(1)) s log s)/C^2`. So the usable hypothesis
   is

   > `QS^sparse(m, delta, f, A)`: solutions of (5.1) with `H(x) >= max(H_0, exp(A s log s))` lie in `<= f(s)`
   > subspaces (any fixed `A < 1/2`).

   In that regime the `S`-integers of norm `<= H` number `H^(O(loglog s/log s)) = H^(o(1))`, so both Prop A and
   round 8's Prop 5.3 (pigeonhole pairs at `H <= P^4`) fail to apply. **Round 8's "`QS(m <= 3)` is false" concerns
   only the dense regime. `QS^sparse(m)` is open for every `m >= 2`.**
4. *Why `QS^sparse` is not progress.* At `m = 2`, `delta < 1/2`, `QS^sparse(2, f)` already gives `M <= f(s) + 1`, and
   it is a near-restatement of the problem: two `S`-integers of height `exp(>= A s log s)` within relative distance
   `H^(-delta/2)`. A random model predicts no such coincidences (`H^(o(1)-delta)` expected), which is the
   Cilleruelo-Granville heuristic itself. So in the sparse regime the needed input is "genuine `S`-integers are
   sparse and do not cluster". That belongs to the genuine-prime class (P3-all, H^gen; `round16/README.md` §7.2),
   not to any existing Diophantine theorem.
5. *Arithmetic of the target.* With `f(s) = s^theta`, Theorem 5.2 gives `M <= s^theta + m`. This beats
   `(1/3) s_max` iff `theta < 1`, giving `M <= (W/log W)^theta`. Nothing in the literature reaches even `theta = 1` in
   the sparse regime. In the dense regime `theta < 1` is impossible (Prop A).

### 2.3 Two-function gcd (Corvaja-Zannier, Levin, Grieve-Wang, Xiao, Huang-Levin, Wang-Yasufuku), max-height normalisation

These theorems bound `log gcd(f(u), g(u)) <= eps max_i h(u_i)`. Our cluster gives
`(u, v) = (z_2/z_1, z_3/z_1)` with ratio `>= 1/2`.

But products of pigeonhole pairs `(rho_ij, rho_kl)` on disjoint prime pairs reach ratio `delta` for every `delta < 1`.
Each pair has `rho = eta conj(u)/u`, `u = pi_i conj(pi_j)`, and `h(rho) = (1/2) log(p_i p_j)`. Check (B),
`checks/out_AB.txt`, counts the pairs among the first `s` split primes with `-log|rho-1|/h(rho) > 1/2`:

| `s` | 50 | 100 | 200 | 400 | 800 |
|---|---|---|---|---|---|
| pairs with ratio `> 1/2` | 112 | 320 | 894 | 2449 | 6646 |

The count grows like `s^(1.45)`, and the maximum ratio is 0.93. At any `eps < 1/2` the exceptional set must contain
all products of such pairs. A translate `u^a v^b = c` with `b != 0` holds `<= 4|b|` points per value of `rho_ij`. So
the number of components is `>= #pairs/(4|b|) ~ s^(3/2-o(1))`.

The cluster count only gives `M <= 4|b| f(s)`. With `f(s) >= s^(3/2)` this is useless. *Moreover, the cluster is not
special for two-function gcds:* its ratio `1/2` is beaten by generic products. This is the same input fakes satisfy
(two archimedean proximities). **Capped** for growth even with optimal `S`-dependence, in the dense regime. In the
sparse regime this case is subsumed by `QS^sparse(3)` (§2.2), which is open.

### 2.4 Single-place theorems (`S = {inf}`): Heier-Levin (#10), Ru-Vojta / Ru-Wang / Vojta beta (#11), Yasufuku 2025 (#8), Barrios (#9)

These are `S`-free, so they are uniform, but they are single-place auxiliary-section arguments. **Capped by
`threshold/upper.md`**: `beta_M <= tau*(M) = 2 - 1/ceil(M/2) < 2`. Concretely:

* **Heier-Levin.** For one subscheme `Y`: `eps_Y(H) lambda_Y(P) <= (dim X_M + 1 + eps) h(P) = (M+1+eps) h(P)`.
  * Clusters have `lambda_Y/h -> 1/2`, so the theorem bites only if `eps_Y(H) > 2(M+1)`.
  * But `eps_Y(H) <= 1`. Check (A) verified this exactly for `M = 6`. For `p = (u,...,u; v,...,v) in Y` and a split
    `J1 ⊔ J2`, set `a_j = 0, b_j = c/u` on `J1` and `a_j = c/v, b_j = 0` on `J2`. The line `p + t(a; b)` lies in
    `X_M = {u_1v_1 = ... = u_Mv_M}`, leaves `Y`, and meets it transversally. So `(H.L)/mult = 1`.
  * The theorem is therefore weaker than needed by a factor `>= 2(M+1)`. Several `Y_i` in general position do not
    help: the finite places give `lambda_w ~ 0`.
* **Products of conics / Seshadri constants of diagonals.**
  * In `(P^1)^M`, the small diagonal has `eps <= 1` (fibres), and on `P^1 x P^1`, `eps_Delta(O(1,1)) = 1` exactly.
  * More importantly, `(P^1)^M` forgets the common radius. The cluster ratio there is
    `lambda_Delta/h ~ (1/2 log R)/(M log 2R) = 1/(2M)`.
  * No Seshadri or approximation-constant statement on products of conics can see the clusters.
* **Barrios.**
  * The corrected theorem needs `N_{Y/X}` ample. For `Y` a line in the smooth locus of the complete intersection
    `X_M ⊂ P^(2M-1)` (`growth/conditional.md` Lemma 2.1 item 4), the sequence
    `0 -> N_{Y/X} -> O(1)^(2M-2) -> O(2)^(M-1) -> 0` gives `deg N_{Y/X} = 0`, so the bundle is not ample and the
    theorem does not apply.
  * Even if it did, its constant comes from sections vanishing along `Y`, so it is `>= 1/tau* > 1/2`, against the
    cluster's `1/2`.
* **Ru-Vojta / Ru-Wang / Vojta beta.** `beta_M <= tau* < 2` (`threshold/upper.md`). The searches found no
  computation of a birational Nevanlinna or beta constant above 2 for these configurations (queries 25, 54).
* **Yasufuku 2025.** The height term is a Ru-Vojta beta term, which is capped. With `S = {inf}`, the "local heights
  outside `S`" term is the whole finite gcd, which carries nothing new.

### 2.5 Unit equations, uniform Mordell-Lang (#14, #14', #15)

* Upper bounds are exponential in the rank: ESS, Amoroso-Viada `(8n)^(4n^4(n+r+1))`, and Beukers-Schlickewei
  `2^(8r+8)`.
* Lower bounds are super-polynomial: `exp((4+o(1))(s/log s)^(1/2))` (ESTT), `exp(s^(2-sqrt2-eps))`
  (Konyagin-Soundararajan), and `exp(c s^(1-1/n)(log s)^(-1/n))` for `n` terms (Granville). So **polynomial-in-rank
  bounds are false**.
* The repo's use (`round8/stability.md` Prop 8.2: unit-equation counts against `M^3` triangle identities) would need
  a count `<< s^3`, which is impossible.
* These lower bounds also live in the dense regime (small heights relative to `s`), as Prop A does. No sparse-regime
  unit-equation count was found.
* DGH (#15) concerns genus `>= 2` curves in abelian varieties and does not apply to tori or the circle.

### 2.6 Uniform dimension growth / Bombieri-Pila (#16)

* All of these are determinant-method or real-analytic parametrisation counts. They improve the dependence on the
  degree.
* For the circle (`d = 2`, `H ~ R`), BCN gives `<< 4 R (log R)^kappa` points, which is trivial. On a window of length
  `C sqrt R`, determinant arguments give only `O_C(log R)` (`growth/auxiliary.md` §2.2 table).
* Structurally, every interpolation determinant on one circle is Vandermonde times a symmetric factor with no extra
  arithmetic (`growth/auxiliary.md` Prop A). These methods are therefore capped by pair data, and the single-place
  versions also by `tau* < 2`.
* Counting points of `X_M` of bounded height does not help either: `X_M` is toric, with many rational points, and the
  constraint is archimedean proximity, which these counts do not see.

### 2.7 Effective Pila-Wilkie / Wilkie (#17)

* Pila-Wilkie and its effective versions count rational or algebraic points on the **transcendental part** of a
  definable set.
* Lattice points on circles, and clusters as points of `X_M`, are semialgebraic, so the transcendental part is empty.
* The only transcendence is the angle map `z -> arg z`. There the relevant statements are linear forms in logarithms
  (Baker/Matveev, single-form LW), which fakes satisfy (`round8/README.md` §2.4).
* Quantitatively, the BNZ bounds are `poly_n(g, log H)` with dimension-dependent degree. At dimension `n ~ M` or `~ s`
  this is `>= W^(kappa(M))`, far above `M`.
* Binyamini's foliation counts concern algebraic points on leaves of a fixed foliation. The cluster points do not lie
  on a fixed leaf, and the target direction drifts with `R`.
* **Inapplicable.**

### 2.8 McKinnon-Roth / Huang approximation constants on toric varieties (#18)

* Over `Q(i)`, `X_M` is split toric (`u_j v_j` all equal), so Huang's setting is close.
* But Huang, like McKinnon-Roth, approximates a **fixed point** `x`. The cluster approximates the point
  `(e^(i phi): ...: e^(-i phi))` of `Y`, and `phi` drifts with the cluster.
* Applying a fixed-point approximation constant at the cluster's centre would **replace a drifting common phase by
  a constant**, which the constraints forbid.
* The uniform-in-`phi` analogue for the curve `Y` is the Zariski non-density statement that `growth/conditional.md`
  Prop 3.2 (§9) shows is *equivalent* to the uniform theorem.
* Not usable.

### 2.9 Temur (#19)

* Temur applies Schmidt's theorem to fixed gap templates `h = (h_1, ..., h_l)` and gets finiteness of `N` for each
  template, ineffectively.
* The threshold depends on the template, and the constraints forbid treating it as uniform in the radius.
* Arithmetic: on an arc of length `C sqrt R` near a coordinate axis the x-gaps range up to `C sqrt R`, so there are
  `>= (C sqrt R)^(l-1)` templates. Even an effective constant would have to be summable over that many templates.
* Temur himself says it gives no new result.
* **Capped by the template constraint.**

### 2.10 Counting near-1 elements of a finitely generated group

The cluster's `M(M-1)` ratios lie in `{gamma in Gamma_S : h(gamma) <= log R, |gamma - 1| <= C R^(-1/2)}`. A uniform
bound `F(s) = o(s^2)` on this set would give growth.

* **Dense regime** (`h` small relative to `s`): Prop A, specialised to `m = 2`, shows the set has `>= exp(s^(1/2-o(1)))`
  elements. False.
* **Sparse regime:** a heuristic count for `s` primes near `P` is about `(5.83/P^(1/4))^s << 1`. Clusters are an
  excess over this.
  * A proof would need either Subspace in the sparse regime (§2.2, open) or angle equidistribution of genuine
    Gaussian primes.
  * The robust analytic versions of the latter are capped (`round10/README.md` §7.1, analytic Thm B), and the
    non-robust ones are the P3-all / H^gen class.
  * No literature result found.

### 2.11 Summary table

| tool | uses an actual-circle property? | blocked by | can give `o(log R/loglog R)`? |
|---|---|---|---|
| CZ / Levin / Grieve-Wang / Xiao / Huang-Levin / Wang-Yasufuku / RTW (fixed `S`) | yes (`S`-integrality via the product formula) | `S`-dependent exceptional set (constraint); quantitative forms `>= exp(s^(1/2-o(1)))` (Prop A) | no; would need `QS^sparse`, which is open and a near-restatement |
| two-function gcd, any uniformity | no (cluster ratio 1/2 beaten by generic pairs) | product configurations: `f >= s^(3/2-o(1))` | no |
| Heier-Levin with Seshadri constants | no (`S = {inf}`) | `eps_Y(H) <= 1`, short by `>= 2(M+1)`; `tau* < 2` | no |
| Ru-Vojta / Ru-Wang / Vojta beta / Yasufuku 2025 | no | `beta_M <= tau* < 2` (`threshold/upper.md`) | no |
| Barrios generalised gcd | no | `deg N_{Y/X_M} = 0` (not ample); constant `>= 1/tau*` | no |
| unit equations / Mordell-Lang uniformity | yes | super-polynomial lower bounds in rank | no |
| uniform dimension growth / Bombieri-Pila | no | pair data (`growth/auxiliary.md`), `tau*` | no |
| effective Pila-Wilkie / Wilkie / foliations | no | empty transcendental part; LW capped | no |
| Huang / McKinnon-Roth toric approximation | no | fixed target point (drifting-phase constraint) | no |
| Temur (Subspace on gap templates) | yes | template-dependent threshold (constraint) | no |

---

## 3. Papers to request in full text (ranked)

1. **Evertse-Ferretti, Ann. of Math. 177 (2013) 513-590, and the survey arXiv 1008.2268.**
   * Read the exact statement of the quantitative Absolute/Parametric Subspace Theorem and its classical corollary:
     how the number of subspaces depends on `s`, on `R` (number of distinct forms), and on the "large solution"
     threshold.
   * The snippet "R may be replaced by ns ... in many cases R can be taken independently of s" must be reconciled
     with Prop A. Prop A shows that with a threshold independent of `S` the count must be `>= exp(s^(1/2-o(1)))`, even
     with only `R = 2m-1` distinct forms. So either the threshold depends on `S` or the count does.
   * **Decisive check:** is there any statement whose threshold is of size `prod_(v in S) N(v)^c` with an
     `S`-independent count? That would be a `QS^sparse`. I expect no, since at `m = 2` it would prove the
     Cilleruelo-Granville conjecture on squarefree norms.
2. **Levin, Invent. Math. (2019), DOI 10.1007/s00222-018-0831-z.**
   * The exact codimension-`n` statement (gcd of `n` functions, blowup of the point `(1,...,1)` in `G_m^n`).
   * Whether the exceptional set is quantified in terms of the rank of `Gamma`.
   * Which height normalisation is used: max-height or projective. This decides whether §2.3 or §2.2 applies.
3. **Evertse-Schlickewei, *A quantitative version of the absolute subspace theorem*, J. reine angew. Math. 548 (2002).**
   Same question as 1. It is the source of the `(9n^2/eps)^(ns)` dependence.
4. **Yasufuku, *GCD inequalities arising from codimension-2 blowups*, Bull. LMS (2026), DOI 10.1112/blms.70349.**
   Check whether any inequality uses the codimension of the blown-up subscheme in a way not covered by the `tau*`
   barrier (codim-2 centre vs our codim-`(M-1)` line).
5. **Temur, Turkish J. Math. 46 (2022) 2867-2889.** Exactly how the Subspace argument handles gap templates, and
   whether the template dependence could be made polynomial (it would have to beat `(C sqrt R)^(l-1)` templates).
6. **Huang, ANT 15 (2021) 461-512.** Whether the universal-torsor method extends to approximating a curve `Y` with a
   drifting target point.
7. **Heier-Levin, Amer. J. Math. 143 (2021)**, and **Barrios + Correction, Bull. Aust. Math. Soc. (2025)**. Low
   priority: only to confirm the exact inequalities used in §2.4. Both are capped whatever the constants.

---

## 4. Search log (WebSearch; mode in brackets)

1. greatest common divisor u-1 v-1 S-units uniform bound Levin Wang [std]
2. Wang Yasufuku greatest common divisors integral points blowups algebraic tori [std]
3. Heier Levin generalization of the Subspace Theorem for closed subschemes [std]
4. Levin "Greatest common divisors and Vojta's conjecture for blowups of algebraic tori" Inventiones [std]
5. Grieve Wang "greatest common divisors with moving targets" linear recurrence sequences [std]
6. Levin Wang "Greatest common divisors of analytic functions and Nevanlinna theory on algebraic tori" [std]
7. Levin gcd S-unit bound log gcd(f(u),g(u)) epsilon max height exceptional Zariski closed set depends on S uniform [ext]
8. Bugeaud Corvaja Zannier gcd(a^n-1,b^n-1) uniform in S Corvaja Zannier gcd(u-1,v-1) S-units bound [ext]
9. Ru Vojta "A birational Nevanlinna constant and its consequences" American Journal of Mathematics 2020 [std]
10. Huang Levin "Greatest common divisors on the complement of numerically parallel divisors" [std]
11. Xiao "Greatest common divisors for polynomials in almost units" explicit exceptional set [std]
12. "Counting Polynomial-type Exceptional Units on Algebraic Varieties over Number Fields" [ext] (irrelevant: units mod n)
13. gcd S-units quantitative number of exceptional subvarieties bound in terms of rank of group explicit Subspace Theorem gcd(u-1,v-1) [ext]
14. "Families of Unit Equations" exponential Diophantine 2026 arXiv 2604.26497 [std]
15. Amoroso Viada small points subvarieties torus number of non-degenerate solutions linear equation multiplicative group rank r bound [ext]
16. Konyagin Soundararajan "Two S-unit equations with many solutions" lower bound exp s [std]
17. Castryck Cluckers Dittmann Nguyen dimension growth conjecture polynomial in the degree without logarithmic factors [std]
18. Binyamini Cluckers Novikov bounds for rational points on algebraic curves optimal in the degree [std]
19. Salberger counting rational points on projective varieties uniform dimension growth Proc London Math Soc 2023 [std]
20. Binyamini Novikov Wilkie conjecture restricted elementary functions polylogarithmic bound effective Pila-Wilkie [std]
21. Binyamini Jones Schmidt Thomas effective Pila-Wilkie theorem Pfaffian functions diophantine applications [std]
22. Binyamini "point counting for foliations over number fields" Forum of Mathematics Pi [std]
23. Seshadri constant of the diagonal in a product of curves P^1 x P^1 C x C [std]
24. Vojta's conjecture moduli space M_{0,n} bar integral points proof special case [std]
25. Nevanlinna constant computation blowup exceeds threshold Ru Wang "Nevanlinna constant" examples projective variety divisor [std]
26. gcd S-unit bound "uniform" in S number of primes polynomial dependence subspace theorem exceptional set count rank 2023 2024 2025 [ext]
27. Evertse Ferretti "further improvement of the quantitative subspace theorem" number of subspaces bound dependence on places [ext]
28. Evertse Schlickewei Schmidt "linear equations in variables which lie in a multiplicative group" exp((6n)^{3n}(r+1)) Annals 2002 [ext]
29. Amoroso Viada "On the zeros of linear recurrence sequences" OR "number of solutions" unit equation (8n)^{4n^4 (n+r+1)} improvement [ext]
30. "Amoroso" "Viada" non-degenerate solutions "rank r" bound exponential in r unit equation improvement Evertse Schlickewei Schmidt [ext]
31. Beukers Schlickewei equation x+y=1 in finitely generated multiplicative groups 2^{8(r+1)} solutions [std]
32. "gcd(u-1,v-1)" S-units exceptions bounded number depending only on epsilon and s Corvaja Zannier quantitative [ext]
33. Pasten Wang "GCD bounds for analytic functions" uniform gcd abc conjecture S-units Silverman [std]
34. McKinnon Roth Seshadri constants Diophantine approximation approximation constant product of curves rational surfaces 2020 2022 [std]
35. Vojta "Birational Nevanlinna constants, beta constants, and Diophantine approximation to closed subschemes" main theorem [ext]
36. "Rational approximations on toric varieties" arXiv 1809.02001 author theorem best approximations rational curves unconditional [ext]
37. Heier Levin "A generalized Schmidt subspace theorem for closed subschemes" American Journal Mathematics Seshadri constant epsilon_Y(A) statement [ext]
38. Evertse Schlickewei "absolute subspace theorem" number of subspaces independent of number field and of S bound 2^{3(n+9)^2} [ext]
39. quantitative subspace theorem "independent of S" OR "does not depend on S" number of subspaces large solutions twisted heights [ext]
40. Evertse Ferretti "On the quantitative Subspace Theorem" survey Theorem B number of subspaces "R" distinct linear forms cardinality s of S [ext]
41. "Evertse" "Ferretti" quantitative Subspace Theorem "number of subspaces" "at most" bound "s" cardinality of S exponential dependence [ext]
42. "An upper bound for the generalized greatest common divisor of rational points" arXiv 2411.06050 [ext]
43. Yasufuku "GCD inequalities arising from codimension-2 blowups" Bulletin London Mathematical Society [std]
44. Dimitrov Gao Habegger uniformity in Mordell-Lang for curves bound c(g,d)^{1+rank} Annals [std]
45. lattice points on short arcs of circles Cilleruelo Granville conjecture new results 2020-2025 [ext]
46. Walsh "Bounded rational points on curves" IMRN uniform Bombieri-Pila determinant method [std]
47. Oganesyan "Lattice points on small arcs" abstract alpha in (1/2,1) ... unbounded [ext]
48. "Discrete fractional integrals, lattice points on short arcs, and diophantine approximation" abstract [std]
49. number of elements of a finitely generated multiplicative group close to 1 bounded height counting approximation S-units |1-x| small Ridout quantitative bound rank [ext]
50. bound polynomial in the number of primes s number of S-unit solutions OR exceptional subspaces "polynomial in s" Diophantine [ext]
51. Levin Inventiones 2019 Vojta conjecture blowup of P^n along closed subscheme codimension two S-integral points coordinate hyperplanes theorem statement [ext]
52. S-unit equation non-degenerate solutions upper bound exp(c s^{1-1/n} (log s)^{-1/n}) large s [ext]
53. Seshadri constant of a curve or subvariety in a toric variety computation blowup along torus-invariant line nef cone [std]
54. Ru Vojta beta constant computation toric variety blowup "beta(L,D)" ... Diophantine approximation 2021 2022 2023 [ext]
55. "GCD inequalities inspired by Vojta's conjecture" Monatshefte 2024 authors theorem all algebraic numbers Ru-Vojta [ext]
56. "lattice points" circle short arc "Subspace Theorem" Gaussian integers S-units close points [ext]
57. "On qualitative aspects of the quantitative subspace theorem" 2306.16583 abstract [std]
58. Temur lattice points on short arcs conics Subspace Theorem "fix h_i beforehand" Conjecture 1 eta > 1/4 [ext]
59. Sert Temur discrete fractional integral operators binary quadratic forms lattice points arcs Subspace theorem Turkish Journal [std]
60. Rousseau Turchet Wang "Divisibility of polynomials and degeneracy of integral points" Vojta blowups toric varieties S-integral [std]
61. Binyamini Novikov Zack "Wilkie's conjecture for Pfaffian structures" Annals bound (log H)^ exponent depends on dimension complexity [std]
62. "Vojta's conjecture" implies lattice points short arcs circle bounded Cilleruelo Granville conditional [ext]
63. Corvaja Zannier book "Applications of Diophantine Approximation to Integral Points and Transcendence" 2018 gcd chapter Subspace [std]
64. o-minimal Pila-Wilkie counting lattice points on circles or arcs Gaussian integers application [std]
65. Ru Wang "The Ru-Vojta result for subvarieties" International Journal of Number Theory 2022 closed subschemes statement [std]
66. Bugeaud Evertse Gyory "S-parts of values of univariate polynomials, binary forms and decomposable forms" bound dependence on s [std]
67. S-unit equation "large solutions" number bounded independently of S height larger than threshold depending on S Evertse [ext]
68. "Towards Lang-Vojta via degeneration" arXiv 2602.06956 abstract [std]
69. Heier Levin "Schmidt-Nochka theorem for closed subschemes in subgeneral position" Crelle 2024 main inequality Seshadri [std]

WebFetch attempts, all failing with DNS errors: ar5iv.labs.arxiv.org/html/1902.09109, par.nsf.gov/servlets/purl/10225835,
pub.math.leidenuniv.nl/~evertsejh/dio19-7.pdf, annals.math.princeton.edu (Evertse-Ferretti PDF),
irma.math.unistra.fr (Bugeaud survey).

Checks (all in `checks/`, run with `python3 -I`):
* `check_tools.py`. (A) lines in `X_M` through `Y` (exact, sympy). (B) pigeonhole pair counts. The script's own
  `__main__` loop is too slow; outputs were produced by running its parts, see `out_AB.txt`.
* `check_smooth2.py`. Prop A tuples: counts, sampled exact (5.1) verification, projective multiplicity
  (`out_C.txt`).
* `count_mu.py`. Mean occupancy `mu_bar` against `s` (`out_mu.txt`).
