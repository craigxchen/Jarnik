# Adversarial check of `lit18/tools.md` (height-based tools uniform in S)

Checker's working files: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/lit18/check_tools_work/`
(`gp_cluster.py` and `gp_cluster_out.txt`, `propA_recount.py`, `lookup.py`). All scripts were run with `python3 -I`.

**Limitation, stated first.** The shared WebSearch budget (200 calls per turn, across all agents) was already used up
when this check started. All three WebSearch calls I made were refused (§6). So **I could not re-verify any citation
with my own web search.** In its place I used:
* (i) an offline arXiv metadata mirror that the searcher had downloaded: `lit18/tools/db.json`, built from
  github.com/ehijano/rss_fetch at commit 30c967b. It holds 92,926 arXiv RSS records, first-announced 2024-03-05 to
  2026-10-26, in math.AG/AP/CA/CO/DS/GM/HO/MG/NT/PR/SP. It contains titles, authors and abstracts. I searched it with
  my own regexes and read it as data only;
* (ii) the independent search logs of the sibling reports `lit18/direct.md`, `lit18/direct_log.txt` and
  `lit18/primes.md`;
* (iii) the repository files;
* (iv) mathematical consistency checks.

Anything verified only by (iv), or only from my background knowledge, is marked **UNVERIFIED** below.

---

## 0. Verdict

Everything below is the checker's own work and is unrefereed. Prop B is elementary, proved in full in §2.2 and checked
numerically.

1. **The bottom line of `tools.md` stands: no literature tool found beats the barrier.** Every lead is fixed-`S`,
   capped by a repo barrier, or inapplicable, and I found nothing in the offline corpus that changes this.
2. **Prop A is correct.** I recounted it independently, and it does refute round 8's `QS(m, delta, f, H_0)` for every
   `m` (correction to `round8/README.md` §2.3/§2.5).
3. **The report's main open "repaired" hypothesis is false.** That hypothesis is `QS^sparse`. The report's §2.2
   item 3 calls it "open for every `m >= 2`", and §2.3 and §2.10 lean on its sparse regime.
   * It is refuted by an elementary construction the report missed: **clusters of Gaussian primes in a small square,
     at arbitrarily large height** (§2.2 below, Prop B).
   * Any hypothesis about tuples of `S`-integers whose threshold depends on `P` only through `s = |P|` must allow
     `f(s) >= s - m + 1` exceptional subspaces. That gives at best `M <= omega(N') + O(1)`, which is worse than the
     proved `(2/3+eps) log R/loglog R` in the worst case.
   * The same construction kills the sparse-regime versions of the near-1 count (§2.10: `F(s) >= s(s-1)`) and of the
     two-function gcd (§2.3).
4. **What an `S`-integer input must use is equality of norms (full support), not "sparse heights".** The report's
   "honest summary" names the wrong separating property. The only form of `QS` that survives both known constructions
   needs a height threshold tied to the *product* of the primes, `H >= (prod_{p in P} p)^c` with `c < 1/2` (§2.5).
   * That form is not in any literature found.
   * It strengthens the target, since it covers unequal-norm tuples as well.
   * It belongs to the genuine-integer class (H^gen; `round16/README.md` §7.2).
   * Verdict: **open, non-standard, not a literature lead.**
5. **New items the searcher missed** (offline corpus, §3). The most interesting is **Ziegler, arXiv 2609.15487
   (Sept 2026)**: an effective, linear-in-`|S|` count for three-term `S`-unit equations once `h(a_3)` is large. It
   contradicts the report's "no sparse-regime unit-equation count was found" and "no `S`-uniform count polynomial in
   `s`". It is **inapplicable**: the cluster's triangle equations (8.1) have coefficients that vary from triangle to
   triangle, have small height and are template-dependent. **Lead fails.**
6. **Citations.**
   * No outright fabrication was detected.
   * 5 arXiv items were confirmed offline (authors, title, abstract).
   * 1 attribution conflicts with a sibling log: arXiv 2302.04209 is BCN (Binyamini-Cluckers-Novikov) in `tools.md`
     but "CCDN" in `direct_log.txt`.
   * About 30 journal-level details (volume, pages, year) are **UNVERIFIED** in this check.
   * The Evertse-Ferretti "`R` independent of `s`" snippet is **logically inconsistent with Prop A** unless the
     threshold or the count carries extra `s`-dependence (§2.8). It is the top full-text request.

---

## 1. Citation re-verification

Status key:
* **OFF-A**: confirmed in the offline arXiv RSS mirror (authors, title, date, abstract).
* **SIB**: independently confirmed by a sibling report's own search log.
* **UNV**: not re-verified here (WebSearch exhausted).
* **FLAG**: inconsistency or overstatement.

| # in tools.md | item | status | evidence and remarks |
|---|---|---|---|
| 1 | Levin, Invent. Math. (2019) 493-533, DOI 10.1007/s00222-018-0831-z | UNV | No contradiction with my background knowledge (Invent. Math. 215). The quoted form ("`log gcd(f(u),g(u)) < eps max h(u_i)` outside a Zariski-closed set") is a **counting function at places outside `S`**; see FLAG at §2.4 |
| 2 | Corvaja-Zannier, Monatsh. Math. 144 (2005); Hernández-Luca; Bugeaud-Corvaja-Zannier, Math. Z. 243 (2003) | UNV / SIB (partial) | `primes.md` #29 lists the same CZ/BCZ statements as "background snippets only". The detail "`|p|,|q| <= 1/eps`" is UNV |
| 3 | Grieve-Wang, Trans. AMS 373 (2020) 8095-8126 | UNV | — |
| 4 | Levin-Wang, Crelle | UNV | — |
| 5 | Wang-Yasufuku, ANT 15 (2021) 287-305 | UNV | — |
| 6 | Huang-Levin arXiv 2207.14432; HLX arXiv 2406.18879; Wang-Xiao arXiv 2506.03534 | 2207.14432 UNV; **2406.18879 OFF-A**; **2506.03534 OFF-A** | 2406.18879 is Keping Huang, Aaron Levin, Zheng Xiao, *A New Diophantine Approximation Inequality on Surfaces and Its Applications* (2024-06-27); the abstract mentions gcd inequalities on complements of three numerically parallel ample curves and families of `S`-unit equations. 2506.03534 is Julie Tzu-Yueh Wang, Zheng Xiao, *Hyperbolicity and GCD for n+1 divisors with non-empty intersection* (v2 2026-03-16). tools.md gives no titles for these two; titles added here |
| 7 | Xiao arXiv 2110.01751 | UNV | — |
| 8, 8' | Yasufuku, Monatsh. Math. 206 (2025) 259-279; Bull. LMS 2026 | UNV | The report already rates 8' as C |
| 9 | Barrios arXiv 2411.06050 + Correction | **OFF-A** (arXiv only) | Benjamín Barrios, *An upper bound for the generalized greatest common divisor of rational points*, v1 2024-11-12. Abstract: "upper bound for the generalized gcd of a point `x` with respect to an irreducible subvariety `Y` ... rather uniform Riemann-Roch type inequality". The Correction ("positively embedded", Matsuzawa counterexample) is **UNV**. The mirror has Matsuzawa arXiv 2507.05027 (gen. gcd along orbits), but its abstract does not mention Barrios |
| 10, 10' | Heier-Levin, AJM 143 (2021) 213-226; Heier-Levin Crelle 2024; He-Ru arXiv 1910.07966 | UNV | — |
| 11 | Ru-Vojta, AJM 142 (2020) 957-991; arXiv 2004.12257; Ru-Wang IJNT 18 (2022); Vojta JTNB 35 (2023) | UNV | — |
| 12 | Rousseau-Turchet-Wang, Math. Ann. 388 (2024) | UNV | The report itself says the attribution is ambiguous |
| 13 | Evertse-Ferretti, Ann. of Math. 177 (2013) 513-590 | UNV + **FLAG** | The snippet "in many cases `R` can be taken independently of `s`" cannot mean "count depends only on `(n, delta, R)` with an `S`-free threshold": Prop A (verified, §2.1) refutes that for `R = 2m-1`. §2.8 |
| 13' | Grieve arXiv 2306.16583 | UNV | — |
| 14 | ESS, Ann. 155 (2002); Amoroso-Viada, Duke 150 (2009); Beukers-Schlickewei, Acta Arith. 78 (1996); Evertse `3 x 7^(2s+1)` | UNV; Beukers-Schlickewei also cited in repo `round8/stability.md` Prop 8.2 | — |
| 14' | ESTT; Konyagin-Soundararajan; Granville (arXiv math/0107219); Harper (arXiv 1108.3819) | ESTT: repo cites Compositio 66 (1988) with the same bound. Others UNV | Granville and Harper are second-hand (C), as the report says. I could not pin what math/0107219 and 1108.3819 are |
| 15 | DGH, Ann. 194 (2021) 237-298 | UNV | Relevance "none (genus 0)" is correct |
| 16 | Salberger PLMS 2023; CCDN 1904.13109; BCN IMRN 2024 (arXiv 2302.04209); BCK 2401.03982; Walsh IMRN 2015; Paredes-Sasyk | **BCK OFF-A**; BCN existence indirectly OFF-A; **FLAG on 2302.04209** | BCK: Binyamini, Cluckers, Kato, *Sharp bounds for the number of rational points on algebraic curves and dimension growth, over all global fields* (v2 2024-11-19). Abstract: "`c d^2 H^(2 d_K/d) (log H)^kappa`", and it generalizes "the results by the first two authors and Novikov from the case `K=Q`", which confirms the BCN result and its shape. `direct_log.txt` Q37 labels 2302.04209 "CCDN ... optimal in degree", while tools.md says BCN. The number-to-author mapping is unresolved; the mathematical content agrees |
| 17 | BN Ann. 186 (2017); BNZ Ann. 199 (2024); BJST JEMS 2026; Binyamini Forum Pi 10 (2022) | UNV | The offline corpus shows related 2026 Binyamini work (arXiv 2609.20668, 2609.24135, 2610.28477), none about lattice points |
| 18 | Z. Huang, ANT 15 (2021) 461-512 | UNV | — |
| 19 | Temur, Turkish J. Math. 46 (2022) 2867-2889 | **SIB** | `direct_log.txt` Q34 has "Turkish J. Math 46(7) (2022)". Pages are UNV |
| 20 | Oganesyan arXiv 2107.09991, withdrawn | **SIB** | `direct_log.txt` Q1/Q4 and `primes.md` #26 confirm the withdrawal ("mistake in Lemma 2") |
| 21 | Bugeaud-Evertse-Győry, Acta Arith. 184 (2018) | UNV | — |
| 22 | Chen-Garcia-Fritz-Mathur-Pasten arXiv 2602.06956 | **OFF-A** | Title, authors and date (2026-02-09) match. The abstract confirms "finiteness and Zariski degeneracy of `S`-integral points ... moduli spaces of curves with extra structure". The relevance verdict (fixed `S`) stands |

**Overstatements flagged in the text of tools.md:**
* §1, "Searched for and not found: any `S`-uniform gcd or Subspace count polynomial in `s`", and §2.5, "No
  sparse-regime unit-equation count was found". Ziegler, arXiv 2609.15487 (2026), is an effective count linear in
  `|S|` in a large-height regime (§3). The *conclusion* (no usable tool) survives; the *search claim* does not.
* §2.2 Prop A, step 8: `log R ~ (2A s log s)^(1/A)` should read `(2 s log s)^(1/A)`, since `P ~ 2 s log P ~ 2 s log s`.
  This does not affect the result.
* §2.2 Prop A, step 5 is vacuous. Distinct points in one square of side `sqrt(R')/2` have `|g| <= |w_i - w_1| <=
  sqrt(R'/2)`, so tuples with `|g| > sqrt(R')` do not exist. Harmless.
* §2.2 Consequence 1, "the truth ... is between `exp(s^(1/2-o(1)))` and `exp(O(s))`". The upper end is Evertse 1996,
  which the repo quotes from memory (`round8/weakest.md` §5.2). UNVERIFIED.

---

## 2. Lead-by-lead kill attempts

### 2.1 Prop A (dense-regime refutation of round 8's `QS`): **verified, correct**

* **Independent recount** (`check_tools_work/propA_recount.py`, my own enumeration over all 4 units and both
  conjugates):

  | `P` | `s` | `R` | `U` | `T_4/T_3` |
  |---|---|---|---|---|
  | 100 | 11 | `10^4` | 207,800 | 2.15 |
  | 200 | 21 | `10^4` | 888,768 | 9.50 |

  This exactly reproduces `checks/out_C.txt`. The rigorous worst case over 200 random cells gives `delta_eff >= 1.542`
  (`P=200`) and `>= 1.591` (`P=100`), above the `m = 4` window `(1, 1.5)`. A common divisor `g != 1` only raises
  `delta_eff`: with `x = w/g`, `delta_x = delta_w log H_w / log H_x >= delta_w`.
* **Proof re-derived:**
  * the finite part is `prod N(x_i)^(-1) <= 4^m H^(-m)` for primitive `S`-integers;
  * the archimedean part is scale-free and `<= (2/R')^(m-1)`;
  * the hyperplane count comes from "a hyperplane with `c_(i0) != 0` holds at most `T_(m-1)` tuples" plus Chebyshev.
    No projective-multiplicity argument is needed, because tuples are counted as vectors.
* **Consequence confirmed:** round 8's `QS(m, delta, o(s), H_0)` is false for every `m >= 2`. The entry in
  `round8/README.md` §2.3 and §2.5 item 2 ("for `m >= 4`, `delta in (1,(m-1)/2)` neither refuted nor supported") must
  be corrected.

### 2.2 `QS^sparse` (§2.2 item 3; heights `>= exp(A s log s)`): **REFUTED, for every `m`, `delta` and `A`**

**Proposition B (Gaussian-prime clusters; elementary).** Fix `m >= 2`, `0 < delta < (m-1)/2`, and any threshold
function `H_0(s)`. For every `s >= m` there is a set `P` of `s` split primes for which the solutions of the round-8
system (5.1) with `H(x) >= H_0(s)` cannot be covered by fewer than `s - m + 1` proper subspaces.

*Proof.*
1. **Choose the square.** Pick `theta' > 0` with `(m-1)(1 - 2 theta') > delta`; this is possible since
   `delta/(m-1) < 1/2`. Gaussian primes of norm in `[x/2, x]` number about `x/log x`. A square of side `x^(theta')`
   inside the open first octant therefore contains, on average, about `x^(2 theta')/log x` of them, so by pigeonhole
   some such square contains `>= s` of them once `x` is large. Take `x >= H_0(s)` as well.
2. **The primes.** Let `pi_1, ..., pi_s` be Gaussian primes in that square. In the open first octant no two are
   conjugate or associate, so their norms `p_j` are `s` distinct split primes. Let `P = {p_j}`.
3. **Each ordered `m`-tuple solves (5.1).** For distinct `pi`'s the tuple is primitive and `S`-integral. The finite
   part is `prod N(pi_i)^(-1) = (1+o(1)) H^(-m)`, and the archimedean part is `<= (2 x^(2 theta')/x)^(m-1)`. So the
   left side is `<= H^(-m-(m-1)(1-2theta')+o(1)) <= H^(-m-delta)`.
4. **Count subspaces.** A hyperplane `sum c_i X_i = 0` with `c_(i0) != 0` contains at most `(s)_(m-1)` of the
   `(s)_m` tuples, because the other coordinates determine `x_(i0)`. Hence at least `s - m + 1` proper subspaces are
   needed. QED

**Numbers** (`check_tools_work/gp_cluster_out.txt`, rigorous worst case over **all** ordered tuples):

| `s` | `A` | `log H` | threshold `A s log s` | square side | `m=2` | `m=3` | `m=4` | `m=5` | subspaces needed |
|---|---|---|---|---|---|---|---|---|---|
| 12 | 0.5 | 22.0 | 14.9 | 15 | 0.749 | 1.498 | 2.247 | 2.996 | `s-m+1` |
| 40 | 0.5 | 78.8 | 73.8 | 49 | 0.897 | 1.794 | 2.690 | 3.587 | 39 / 38 / 37 / 36 |
| 40 | 1.0 | 152.6 | 147.6 | 64 | 0.944 | 1.889 | 2.833 | 3.778 | 39 / 38 / 37 / 36 |
| 100 | 0.5 | 235.3 | 230.3 | 134 | 0.956 | 1.911 | 2.867 | 3.822 | 99 / 98 / 97 / 96 |

The `m=2..5` columns are `delta_eff`. Every `delta_eff` exceeds `(m-1)/2`, so the window for every `delta < (m-1)/2`
is covered with room to spare. Primality is sympy's `isprime` (BPSW above `2^64`).

**Arithmetic of the consequence.**
* Theorem 5.2 gives `M <= f(omega(N')) + m - 1`. With `f(s) >= s - m + 1` forced, the best the `QS` route can ever give
  is `M <= omega(N')`.
* `omega(N')` reaches `(2+o(1)) log R/loglog R`, which is **three times the proved bound** in the worst case. So no
  growth improvement and not even the constant.
* §2.2 item 5's target `f(s) = s^theta`, `theta < 1`, is impossible in every regime, not only the dense one.

**Why the report's heuristic misled.** The report's estimate "`H^(o(1)-delta)` expected coincidences" fixes `S` first
and then asks for clustering. But `QS` quantifies over all `P`, so the adversary chooses the primes after the cluster.
The construction works at every height, so no `s`-dependent threshold (sparse or not) can repair it.

### 2.3 Near-1 elements (§2.10): sparse-regime target `F(s) = o(s^2)` **REFUTED**

Take the same `s` Gaussian primes. All `s(s-1)` ratios `gamma = pi_j/pi_k` satisfy:
* `h(gamma) = log R_gamma`;
* `|gamma - 1| <= C R_gamma^(-1/2)` with `C = 1`.

Verified for 132/132, 1560/1560, 1560/1560 and 9900/9900 ratios in the runs above. So `F(s) >= s(s-1)`, which gives
at best `M(M-1) <= F(s)`, i.e. `M <~ s`. No growth.

The cluster ratios additionally satisfy `|gamma| = 1`. With that restriction added, the set is that of
`S`-integers `alpha` with `N(alpha) <= R^2` in a sector of width `~R^(-1/2)`. That is a linear form in the logarithms
of `pi_j/conj(pi_j)`, so it falls under the Baker/LW caps (`round8/README.md` §2.2 rows 5-6), or else it is the
problem itself.

### 2.4 Two-function gcd (§2.3): capped in **both** regimes, with one framing FLAG

* **Dense regime.** The pigeonhole-pair counts are correct as reported (`out_AB.txt`: 112 ... 6646, slope about
  1.47).
* **Sparse regime.** The report deferred this to `QS^sparse(3)`. Gaussian-prime triples `(pi_j/pi_1, pi_k/pi_1)` have
  archimedean proximity ratio `lambda_inf/h >= 0.897` (`s = 40`) and `>= 0.957` (`s = 100`), against the cluster's
  `1/2`. Lines needed `>= s - 2`. So it is capped in the sparse regime too.
* **FLAG: two notions are mixed.** CZ, Levin, Grieve-Wang, Xiao and Huang-Levin bound the gcd **at places outside
  `S`** (the counting function `N_{Y,S}`). The cluster's large quantity is archimedean proximity, at a place **in**
  `S`. For `S`-units, proximity at `v in S` is controlled by Ridout/Subspace, i.e. system (5.1), not by the gcd
  theorems. The finite gcd of cluster differences outside `S` is a local input, already capped (`README.md`
  growth-rate status, "gcds"). The conclusion is unchanged.

### 2.5 What survives of the `QS` idea (`QS^prod`): open, non-standard, **not a literature lead**

A threshold that both constructions miss is

> `H(x) >= (prod_{p in P} p)^c` with `c < 1/2`.

Clusters satisfy it, because `H = N/|g|^2 >= N^(1/2)/C^2 >= N'^(1/2)/C^2`.

| construction | its height | why `QS^prod` excludes it |
|---|---|---|
| Gaussian-prime clusters | `H ~ (prod p)^(1/s)` | below the threshold once `s > 1/c` |
| dense smooth tuples (Prop A) | `H <= exp(O(sqrt(s log s)))` | far below `exp(c s log s)` |
| full-support equal-norm tuples (conjugation patterns, e.g. `x_1 = conj(pi_1) pi_2 ...`, `x_2 = pi_1 conj(pi_2) ...`) | — | relative distance `~ |phi_i - phi_j| ~ H^(-1/(2J))`, far from `H^(-1/4)` |

`QS^prod(m, delta, o(s))` would give growth through Theorem 5.2. But:
* it is a strict **strengthening** of the target, since it covers unequal-norm tuples;
* at `m = 2` it says that `S`-integers whose height is a positive power of `prod p` do not come in close pairs;
* it is in the H^gen class (`round16/README.md` §7.2);
* no literature was found for it, and the WebSearch budget prevented a further search.

**Lead fails as a literature lead.** Record it in the barrier map as "open, non-standard, equivalent in spirit to the
problem".

### 2.6 Single-place items (§2.4): kills confirmed

* **Heier-Levin.**
  * `dim X_M = 2M-1-(M-1) = M`, so the right side is `(M+1+eps) h`.
  * The cluster has `lambda_inf/h -> 1/2`. Even if every place of `S` were counted, so that the finite gcd of the
    differences is included (`<= log(C sqrt R)`), the total is `<= 1 + o(1)`.
  * `eps_Y(H) <= 1`, because transversal lines in `X_M` meet `Y`: `out_AB.txt` (A) is True/True.
  * So the inequality needs `eps_Y > M+1`. **Capped** (also by `threshold/upper.md`).
* **Barrios.** I recomputed the normal bundle. `X_M = {u_1v_1 = ... = u_Mv_M}` is smooth along `Y`: the Jacobian rows
  `(v_1, u_1, -v_j, -u_j)` have disjoint supports outside the first pair. Then
  `0 -> N_{Y/X} -> O(1)^(2M-2) -> O(2)^(M-1) -> 0` gives `deg N_{Y/X} = 0`, so it is not ample. **Inapplicable,
  conditional on the Correction's hypothesis being "ample normal bundle" (UNV).**
* **Ru-Vojta / Vojta beta / Yasufuku 2025.** `beta_M <= tau*(M) < 2` (`threshold/README.md`). **Capped.**
* **Zhao, arXiv 2502.08049** (OFF-A; a Seshadri-weighted Schmidt theorem for points of bounded degree, missed by the
  searcher). Same cap as Heier-Levin.

### 2.7 Other families: kills confirmed

| family | verdict |
|---|---|
| unit equations / uniform Mordell-Lang | Fixed-threshold counts are super-polynomial (ESTT). The repo's needed count (`stability.md` Prop 8.2: `binom(M,3)` triangle equations) is unreachable. **Capped**; Ziegler in §3 |
| dimension growth (BCN/BCK/Salberger/Walsh/CCDN) | The circle has `d = 2` and `H ~ R`, so the count `4R(log R)^kappa` is trivial. Determinant methods reduce to pair data (`growth/auxiliary.md` §2.2: `O_C(log R)`; checked). **Capped.** Bhattacharjee arXiv 2607.24048 (OFF-A, a larger sieve over `O_K` with Bombieri-Pila-type bounds) is in the larger-sieve class, which `growth/auxiliary.md` also caps |
| Pila-Wilkie / Wilkie / foliations | The transcendental part is empty. **Inapplicable** |
| Huang / McKinnon-Roth | Fixed target point, against a drifting phase (constraint). `growth/conditional.md` Prop 3.2(ii) equivalence checked. Lehmann-McKinnon-Satriano arXiv 2403.02480 (OFF-A; McKinnon's conjecture for split surfaces via Vojta) has the same fixed-point structure. **Inapplicable** |
| Temur | Template-dependent threshold (constraint). **Capped** |

### 2.8 The Evertse-Ferretti snippet: a definite consistency test

* Our system `Sigma_m(P)` uses only `R = 2m - 1` distinct forms (`X_1..X_m`, `X_i - X_1`), with all heights 1.
* Prop A shows that solutions with `H >= H_0` (any fixed `H_0`) need `>= exp(s^(1/2-o(1)))` subspaces. Prop B shows
  that solutions above *any* threshold `H_0(s)` need `>= s - m + 1`.
* Therefore any EF corollary of the form "count `<= F(n, delta, R)`, threshold `max(H(L), n^(4n/delta))`" is false
  for this system. The full text must show either an `s`-dependent factor or an `S`-dependent threshold.
* The ESTT lower bound does **not** give this test by itself: the 3-term unit equation needs a per-place choice of
  forms, which is where `3^s` enters. Prop A and Prop B do give it, because their system is fixed.
* This should be checked first when the PDF arrives. Even the most favourable reading cannot beat Prop B's
  `s - m + 1`.

---

## 3. Items missed by the searcher (offline corpus), each with a kill attempt

All of these are OFF-A: seen in `lit18/tools/db.json` through my regex queries (§6), with the abstract text quoted or
paraphrased.

1. **V. Ziegler, *On the number of solutions to S-unit equations*, arXiv 2609.15487 (2026-09-15).**
   * *Abstract:* "for fixed `a_1, a_2` the equation [`a_1 x + a_2 y = a_3`, `x, y in U_S`] has at most `(|S|-1)/2`
     solutions provided that the height of `a_3` is sufficiently large ... our results are effective."
   * **Kill.** The only `S`-unit equations a cluster supplies are the triangle identities (8.1) of
     `round8/stability.md`: `A_w beta_xy + A_x beta_yw = A_y beta_xw`.
     * The coefficients are the `beta`'s. They change with each triangle and have height `M^(O(1))`
       (`stability.md` §8.1).
     * After normalisation `a_3 = 1` (small), or `a_3` is an `S`-unit multiple, which is not a large-height
       coefficient in any meaningful sense.
     * So there is no fixed `(a_1, a_2)` with many solutions, and the large-`h(a_3)` hypothesis fails.
     * The threshold presumably depends on `(S, a_1, a_2)`: a template-dependent threshold (constraint).
     * Fakes satisfy (8.1) verbatim with complex `beta` (Prop 8.2(2)), so the only fake-separating content is the
       integrality of the `beta`'s, which a per-equation count cannot use.
   * *Hypothetical arithmetic:* if all `binom(M,3)` triangles were instances of one equation, then
     `M^3/6 <~ |S|/2` would give `M <~ (3|S|)^(1/3)`. But they are not.
   * **Lead fails.** Low-priority request (§5), only to see how the threshold depends on `S`.
2. **F. Calegari, V. Dimitrov, Y. Tang, *Arithmetic holonomy bounds and effective Diophantine approximation*, arXiv
   2510.04156 (2025).** A new effective method on `P^1` and `G_m`. It recovers the effective two-variable `S`-unit
   equation, Thue-Mahler and Baker-Feldman results, and gives effective measures of linear independence.
   * **Kill:** any lower bound for a *single* linear form in logarithms, up to Lang-Waldschmidt strength, is satisfied
     by the skew-double-core fakes (`round8/README.md` §2.2 rows 5-6, `lw.md` Thms A-B). Its `S`-unit applications are
     fixed-`S`. **Capped.**
3. **J. T.-Y. Wang, Z. Xiao, *Families of Unit Equations and Exponential Diophantine Problems via Integral Points*,
   arXiv 2604.26497 (2026).** Ru-Vojta and higher-dimensional HLX inequalities, giving degeneracy for one-parameter
   families of unit equations and new gcd estimates. tools.md's query 14 hit it, but it is absent from the results
   table.
   * **Kill:** fixed `S`; capped by the exceptional-set constraint, and its Ru-Vojta part by `tau* < 2` when
     `S = {inf}`.
4. **Grieve-Noytaptim, arXiv 2407.08614; Grieve, arXiv 2609.36300; Salami-Shaska, arXiv 2309.10300; Mello-Yasufuku,
   arXiv 2604.03745; Bilu-Luca, arXiv 2607.28857; Nair-Kumar-Rout, arXiv 2506.02898.** Orbits, Vojta's harder
   implication, weighted projective Vojta, and fixed finitely generated groups. Each is fixed-`S` or about a fixed
   variety. **Irrelevant or capped.**

---

## 4. Corrections to apply

1. **tools.md §2.2 items 3-5 and §0.** Replace "`QS^sparse` is open for every `m`" with "`QS^sparse` is false for
   every `m`, `delta` and `A`, and for every `s`-dependent threshold (Prop B: Gaussian-prime clusters,
   `f(s) >= s - m + 1`)". The only surviving form is `QS^prod` (§2.5): open, non-standard, not in the literature.
2. **tools.md §2.10 and §2.3 (sparse-regime remarks).** Replace "open" with "refuted": `F(s) >= s(s-1)`; two-function
   lines `>= s - 2`.
3. **tools.md honest summary.** The separating property of actual clusters is **equal norms (full support)**, not
   "heights tied to `S`". Gaussian-prime clusters are genuine `S`-integers of enormous height relative to `s`, and they
   satisfy every equal-norm-blind input.
4. **New barrier-map entry** (proposed for `round8/README.md` §2.4, unrefereed). Hypotheses on tuples of `S`-integers
   whose quantifier depends on `P` only through `|P|` and a threshold `H_0(|P|)` are capped at `M <= omega(N') +
   O(1)`. Witness: Gaussian primes in a square of side `x^(theta')`, at any height. To escape, an input must use
   equality of norms or a threshold tied to `prod_{p in P} p`.
5. **round8/README.md §2.3 row "QS" and §2.5 item 2.** Replace "neither refuted nor supported for `m >= 4`" with
   "false for every `m` (tools.md Prop A at fixed `H_0`; Prop B at any `H_0(s)`)".
6. **Minor (tools.md).**
   * Prop A step 8 should read `log R ~ (2 s log s)^(1/A)`.
   * Step 5 is vacuous.
   * Record titles for arXiv 2406.18879 and 2506.03534.
   * Resolve the BCN/CCDN label for arXiv 2302.04209.
   * Soften the "not found" claims (Ziegler).

---

## 5. Papers to request in full text (re-ranked)

1. **Evertse-Ferretti, Ann. of Math. 177 (2013), and the survey arXiv 1008.2268.** Apply the §2.8 test.
   * Find where `s` enters: the threshold, the count, or `R`.
   * If the corollary claims an `S`-free count for a fixed system, it contradicts Props A and B; then re-read the
     hypotheses (absolute values, the twist, "large solutions").
2. **Ziegler, arXiv 2609.15487.**
   * Exact meaning of "height of `a_3` sufficiently large": relative to `S`-units? effective threshold in terms of `S`?
   * Whether any multi-equation or averaged version could take the varying `beta` coefficients of (8.1).
   * I expect no, given §3.1.
3. **Barrios + Correction, Bull. Aust. Math. Soc. (2025).** Only to confirm "ample normal bundle" as the corrected
   hypothesis. The kill otherwise rests on the computation `deg N_{Y/X_M} = 0`.
4. **Levin, Invent. Math. (2019).** Only to confirm that the main theorem bounds the counting function outside `S`
   (the §2.4 FLAG), not proximity at `S`. Low priority: either way it is capped.
5. **Calegari-Dimitrov-Tang, arXiv 2510.04156.** Low priority. Only if someone wants to check whether its
   linear-independence measures are multi-form (simultaneous) in a way that escapes the single-form LW cap. The
   abstract gives no indication of that.

The report's items 3-6 (Evertse-Schlickewei, Yasufuku Bull. LMS, Temur, Huang) are demoted. None can change the
verdict given Props A and B and the template and fixed-point constraints.

---

## 6. Log

**WebSearch.** Three queries, all refused because the shared per-turn budget was exhausted:
1. `Evertse Ferretti "A further improvement of the Quantitative Subspace Theorem" Annals of Mathematics 2013 number of
   subspaces` [ext]
2. `Levin "Greatest common divisors and Vojta's conjecture for blowups of algebraic tori" Inventiones Mathematicae 215
   2019` [std]
3. `Corvaja Zannier "A lower bound for the height of a rational function at S-unit points" Monatshefte gcd(u-1,v-1)`
   [std]

**Offline arXiv-mirror lookups** (`check_tools_work/lookup.py` on `lit18/tools/db.json`):
* ID lookups: 2406.18879, 2506.03534, 2602.06956, 2401.03982, 2411.06050 (found). 2207.14432, 2110.01751,
  2302.04209, 2301.09883, 2308.11460, 2306.16583, 2004.12257, 1910.07966, 2106.11337, 2107.09991, 2012.10784 (not in
  the mirror window).
* Regexes:
  * `Barrios|Matsuzawa` together with `common divisor|gcd|Riemann`
  * `Yasufuku`
  * `subspace theorem`
  * `greatest common divisor|\bgcd\b`
  * `lattice points` together with `arc`
  * `S-unit|S-integer` together with `uniform|independent of S|number of primes|rank`
  * `Vojta` together with `blow`
  * `Nevanlinna constant|beta constant|Ru-Vojta`
  * `\bS-unit`
  * `Vojta`
  * `Seshadri`
  * `Ridout|quantitative Roth|quantitative subspace`
  * `unit equation`
  * `smooth numbers`
  * `Cilleruelo`
  * `lattice points on (the )?circle|points on circles|on a circle`
  * `Gaussian (integers|primes)`
  * `linear forms in logarithms|Lang.Waldschmidt`
  * `Pila|Wilkie|o-minimal`
  * `multiplicatively (in)?dependent`
* Abstracts read: 2507.05027, 2604.26497, 2309.10300, 2411.01247, 2502.08049, 2407.08614, 2506.02898, 2403.20037,
  2609.15487, 2503.03985, 2510.04156, 2403.02480, 2506.04087, 2604.23918, 2603.17732, 2601.11376, 2609.36300,
  2607.28857, 2607.24048, 2602.23597, 2606.07959.

**Computations.**
* `gp_cluster.py` (Prop B, near-1 set, two-function ratio) → `gp_cluster_out.txt`.
* `propA_recount.py` (independent Prop A recount; output in this file, §2.1).
* Repo cross-references read:
  * `README.md`;
  * `growth/README.md`;
  * `threshold/README.md`;
  * `round10/README.md` §7;
  * `round16/README.md` §7.2;
  * `round8/README.md` §1.2 and §2.2-2.5;
  * `round8/weakest.md` §5;
  * `round8/stability.md` §8.1-8.2;
  * `growth/auxiliary.md` §2.2;
  * `growth/conditional.md` Prop 3.2.
