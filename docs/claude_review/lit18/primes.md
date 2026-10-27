# lit18 / primes: genuine-prime statistics (H^gen, L1, L7). Literature sweep, 2018-2026

**Scope.** This sweep looks for results that could supply the route-3 growth inputs H^gen + L1 (and L7 = (I3)).
These are statements about the residues and angles of the **genuine** Gaussian primes that divide `R`.
Topics covered:
* angles of Gaussian primes in very short sectors;
* Gaussian primes in thin or sparse sets;
* joint statistics of the angles of several Gaussian primes;
* linear relations among angles and lattice points;
* multiplicative-independence measures;
* added because they match H in shape: Poisson/Kubilius models of prime factors, and the new "combinatorial
  large sieve".

**Method.**
* WebSearch: 46 queries, standard and extended (§4).
* Regular-expression search of a local mirror of the daily arXiv RSS listings, 2024-03-05 to 2026-10-27, for
  math.NT/CO/CA/PR/DS/MG/SP/GM/AG/AP/HO: 92,926 records with titles and abstracts. The sibling agent built it at
  `.../scratchpad/lit18/tools/db.json` from `ehijano/rss_fetch@30c967b` on GitHub, and I only read it.
* arXiv, ar5iv, Semantic Scholar and the other databases returned 403 (re-checked). No full texts were read.
* Every citation below was seen in a search snippet or an RSS abstract. The confirming snippet is quoted or
  paraphrased in §1.
* Arithmetic checks: `.../scratchpad/lit18/primes_checks/arith.py`, with output in `arith_out.txt` (same folder).

---------------------------------------------------------------------------------------------------------------

## 0. Bottom line

**Nothing found plausibly beats the barrier.** No 2018-2026 result supplies H^gen, L1 or L7 at the strength needed:
uniform over all `R`, worst-case `R`, cluster primes up to `M^(h_0)`, window moduli in `(M^(1-eta), M]`. No result
supplies a weaker form that would give `o(log R/loglog R)` either. Every result found falls into one of six
classes, and each class is already blocked by a theorem in the record or by explicit arithmetic (§2):

| class | examples found | why it cannot give growth |
|---|---|---|
| ensemble statistics of Gaussian-prime angles and residues (sectors, variances, Bombieri-Vinogradov in sectors, Hecke zero density) | Huang-Liu-Rudnick 2020; Järviniemi-Teräväinen 2024; Rudnick-Waxman; Chen et al. 2023; Katz 2016; Stucky 2021; Lester-Waxman-Yesha 2026; Kuś 2026; Halupczok 2017; Khale et al.; Musicantov-Zehavi 2024 | robust and prime-supported, so capped by `round10/analytic.md` Thm B (graft barrier). "Almost all sectors" exceptional sets contain any `O(M)`-point selection (Prop. C, §4.2) |
| existence and asymptotics of Gaussian primes in thin sets | Merikoski (Compositio 2025; PLMS 2023); Green-Sawhney 2024/26; Friedlander-Iwaniec JEMS 2022; Pratt; Maynard; Heath-Brown-Li; Sabuncu 2026 | they give lower bounds for primes in a thin set. The need is non-existence, or worst-case statistics of `O(M)` adversarially chosen primes. Used as external witnesses they become PCG/GPA, which collapses to the problem itself (Thm D, Prop. D.1). The one result robust to adversarial selection (Merikoski, arbitrary `B`) needs `\|B\| >= X^0.4`, while the cluster selection is `<= X^(1/h_0)` |
| angular statistics of lattice points on circles, for generic `n` | Kurlberg-Lester (Invent. 2025); Sartori; Granville-Wigman; Kurlberg-Wigman; Bombieri-Bourgain 6-correlations; Parkkonen-Paulin 2026 | density-one in `n` (near-extremal `N` are density zero), or energy/correlation upper bounds that fakes also satisfy (`local_smoothing/chains.md` row E). Sartori's "generic" model *is* the i.i.d.-angle fake model |
| Diophantine and multiplicative independence | Ferreira-Ribas 2026; Pappalardi-Sha-Shparlinski-Stewart 2018; Corvaja-Zannier/Levin gcd; Gou 2026 | single-form Baker (capped), counting on average, fixed-`S` (forbidden by the project constraints), or the wrong direction for L1 |
| Poisson/Kubilius models for prime factors | Ford (IMRN 2025); Gorodetsky (2026) | **the right shape for H**, and it is a fake-violated kind of statement (round-6 `D_A` violates H). But H needs level of distribution `D ~ M^(8.7)` at `M = 10^6` (log-ratio `-> infinity`) for a family of `M^2/2` integers, and that family is chosen by an adversary (§2.5) |
| combinatorial large sieve (collisions mod many small primes vs a global multiplicity bound) | Croot-Mao-Pohoata-Sheffer-Yip (Jun 2026); inverse-sieve papers by Croot-Mao-(Pohoata-)Yip | **the same mechanism as `two_thirds/proof.md`.** Its saving `exp(-c log N/loglog N)` is the current rate. Its inputs are per-prime residues plus a global size bound, the local class capped by fakes (§2.6) |

**Structural reasons, independent of the literature.**
1. *H^gen + L1 is a split of the growth statement, not a reduction of it* (round 9 Thm 3.2, Cor. 3.4). Under L1
   every `|B_xy| <= M^b`. A heavy pair, which is what H forces, needs `|B_xy| >= M^((1-eta)k)`, and `(1-eta)k` is
   3.5, 5.0, 6.5 at `M = 10^3, 10^6, 10^9` (`arith_out.txt`). The two clash as soon as `(1-eta)k > b`. So any
   literature input would have to prove one half against an adversarial cluster.
2. *H is a statement about the far tail for one adversarial family.* It needs some pair of one cluster with
   `k = floor((2 - eps/4) log M/loglog M)` window primes dividing `B_xy`. For a random pair that event has
   probability `~ M^(-2+eps/4)`. The joint moduli are `>= M^((1-eta)k)`: `10^30` at `M = 10^6`, against
   `M^2/2 ~ 10^11.7` pairs. No equidistribution, large-sieve or Kubilius statement controls a family at moduli
   beyond its own size. The family is also selected worst-case, so averaged theorems do not apply.
3. *Fake test.* A literature input beats the barrier only if actual clusters satisfy it and fakes violate it.
   * The results about the Gaussian-prime ensemble, thin sets, generic `n`, and correlation or energy upper bounds
     are all satisfied by grafted or angle-generic fakes.
   * The fake-violated inputs found are the pointwise lattice facts: Rudnick-Waxman's angle-gap lemma, and the
     arithmetic behind the "lattice regularity" of Gaussian prime angles. These are `round10/analytic.md` Thm D,
     i.e. the original problem.

**Side findings.**
* Oganesyan's claimed disproof of the `R^alpha` (`alpha > 1/2`) Cilleruelo-Granville conjecture (arXiv 2107.09991)
  is withdrawn ("mistake in Lemma 2"). This is already recorded in `round4/hardness.md` and
  `growth/referee_conditional.md`, and it is consistent with the Vojta-conditional theorem.
* Huang-Liu-Rudnick, Järviniemi-Teräväinen and Merikoski are new relative to the record (the repository does not
  cite them), but they change nothing.

---------------------------------------------------------------------------------------------------------------

## 1. Results found

Confidence is about whether the stated content is right, judged from snippets or RSS abstracts:
* **H**: abstract or theorem text seen;
* **M**: partial snippet or secondary description;
* **L**: bibliographic data only.

Relevance to H^gen/L1/L7: **none**, **low** (same objects, wrong kind of statement), **shape** (the right kind of
statement, wrong strength or scope), **mechanism** (the same proof mechanism as something the record already
caps).

### 1.1 Angles of Gaussian primes: short sectors, variance, polar boxes

1. **Huang, Liu, Rudnick, "Gaussian primes in almost all narrow sectors", Acta Arith. 193 (2020) 183-192,
   arXiv 1903.04005, DOI 10.4064/aa190331-23-7.** [H]
   * Statement (Theorem 1 as quoted in the snippet): "Let `0 <= rho < 3/5`. Then (1.1) holds for almost all short
     sectors `I ⊂ [0, pi/2)` of length `\|I\| > X^(-rho)`."
   * It follows from an upper bound for the number variance, via zero-density theorems. Järviniemi-Teräväinen
     restate it as "almost all sectors of area `X^(2/5+eps)`".
   * Relevance: low.
2. **Järviniemi, Teräväinen, "Gaussian almost primes in almost all narrow sectors", Rev. Mat. Iberoam. 40(4)
   (2024), arXiv 2303.05822.** [H]
   * Thm 1.1: almost all sectors `{N(n) <= X, theta <= arg n < theta + h/X}` contain a product of exactly three
     Gaussian primes for `h = (log X)(loglog X)^19.2`.
   * Thm 1.2: products of exactly two Gaussian primes for `h = (log X)^15.1`.
   * The method is Matomäki-Radziwiłł adapted to `Z[i]`.
   * Relevance: low.
3. **Rudnick, Waxman, "Angles of Gaussian primes", arXiv 1705.07498** (journal version not confirmed in the
   snippets). [M]
   * The paper studies the variance of the number of prime angles in a short arc. It gives an RMT-motivated
     conjecture, proves the function-field analogue, and gives a GRH "almost all sectors contain a prime angle"
     statement (the snippet couples this with Parzanchevski-Sarnak).
   * It also notes that the angles repel: the minimal gap is bounded below by a power of `x`, and
     `theta_a >> 1/sqrt(Norm a)`.
   * Relevance: low. The repulsion is the lattice fact of Thm D.
4. **Chen, Kim, Lichtman, Miller, Shubina, Sweitzer, Waxman, Winsor, Yang, "A refined conjecture for the variance
   of Gaussian primes across sectors", Exp. Math. (2023), arXiv 1901.07386.** [H]
   * A ratios-conjecture derivation with a power-saving error, showing two bifurcation points.
   * It proves an unconditional result for very small sectors.
   * Relevance: none.
5. **N. Katz, "Witt vectors and a question of Rudnick and Waxman", IMRN 2016, DOI 10.1093/imrn/rnw130.** [M]
   Equidistribution of L-functions of "super-even" characters of truncated Witt vectors. This is the function-field
   input behind the Rudnick-Waxman variance and the Katz-Sarnak picture for angles. Relevance: none.
6. **Lester, Waxman, Yesha, "The variance of the number of lattice points in narrow sectors", arXiv 2609.23214
   (21 Sep 2026).** [H, RSS abstract]
   * The variance of the number of lattice points in smoothed narrow sectors, averaged over the sector position.
   * Three regimes of opening angle, with an asymptotic in each.
   * Relevance: none (an average over positions).
7. **J. Stucky, "Gaussian primes in narrow sectors"** (arXiv 2020, published 2021, venue not confirmed). [M]
   "We generalize a Theorem of Ricci and count Gaussian primes `p` with short interval restrictions on both the
   norm and the argument of `p`" (Heath-Brown's method). Relevance: low.
8. **Background, not 2018+.** [M]
   * Harman-Lewis, "Gaussian primes in narrow sectors", Mathematika 48 (2001): lower bounds in sectors of area
     `X^0.619`.
   * Harman-Kumchev-Lewis: the nearest Gaussian prime to any `z` is `<< \|z\|^0.53` away.
   * Matsui: least Gaussian prime with `alpha < arg w < beta` (cited in Harman-Lewis's references).
9. **J. Kuś, "Halász theorems for Gaussian ideals in sectors and short intervals", arXiv 2605.19105
   (May 2026).** [H, RSS]
   * A quantitative Halász theorem on ideals of `Z[i]`, with a sectorial analogue under angular non-pretentiousness.
   * A sectorial short-interval version for annular sectors whose radial thickness tends to infinity.
   * Relevance: none.

### 1.2 Equidistribution in residue classes jointly with angles

10. **Musicantov, Zehavi, "Sectorial equidistribution of the roots of x^2+1 modulo primes", Q. J. Math. 75(2)
    (2024), arXiv 2112.07494.** [M]
    * Joint equidistribution of the normalised roots `nu/p` of `nu^2 = -1 (mod p)` and the angles of
      `p = a^2 + b^2`.
    * Power-saving linear and bilinear sums, combined with a sieve adapted from Duke-Friedlander-Iwaniec.
    * Relevance: low.
11. **Grimmelt, Merikoski, "On the greatest prime factor and uniform equidistribution of quadratic polynomials",
    arXiv 2505.00493 (2025).** [H, RSS]
    * `P^+(n^2 + h) >= n^1.312` infinitely often.
    * The root equidistribution of DFI is uniform in `h <= n^(1+o(1))` under a hypothesis on real characters.
    * Relevance: none.
12. **Halupczok, "A Bombieri-Vinogradov theorem with products of Gaussian primes as moduli", Funct. Approx. 57(1)
    (2017), arXiv 1607.07265** [M]; **Khale, O'Kuhn, Panidapu, Sun, Zhang, "A Bombieri-Vinogradov theorem for primes
    in short intervals and small sectors", arXiv 2008.09677** [M].
    * Ensemble equidistribution of Gaussian primes in residue classes, on average over the moduli, also in sectors.
    * This is the distribution from which the cluster labels `ell_j(q)` are drawn.
    * Relevance: low.
13. **Lin, Tang, Guo, "The asymptotic estimation of prime ideals in imaginary quadratic fields and Chebyshev's
    bias", arXiv 2412.14792** [H, RSS]; **Devin, "Discrepancies in the distribution of Gaussian primes", arXiv
    2105.02492** [M]. Prime ideals with congruence and argument conditions, and biases. Relevance: none.
14. **Thorner, Zaman: explicit log-free zero density and zero repulsion for Hecke L-functions (arXiv
    1510.08086)** [M]; **Brumley, Thorner, Zaman, JEMS 24 (2022) 1471-1541** [M]. Chebotarev least-prime bounds.
    Relevance: none (robust, R2-type).

### 1.3 Gaussian primes in thin or sparse sets (existence and asymptotics)

15. **J. Merikoski, "On Gaussian primes in sparse sets", Compositio Math. (2025), arXiv 2302.11331 (v3
    Sep 2024).** [H, RSS]
    * There is `delta > 0` such that every `B` with `\|B ∩ [1,Y]\| >> Y^(1-delta)` gives infinitely many primes
      `a^2 + b^2` with `b in B`.
    * There is a quasi-explicit formula, in terms of zeros of Hecke L-functions, for `\|B\| = X^(1/2-delta)`,
      `delta < 1/10`, `B ⊂ [eta X^(1/2), (1-eta) X^(1/2)]`. It gives an asymptotic unless `B` has a large subset of
      multiples of a fixed large integer.
    * Remark 1.4: restricting `a` to a narrow interval puts the primes in a very narrow sector.
    * Relevance: low. Uniformity over **arbitrary** `B` is the only adversarial-selection feature found anywhere.
16. **J. Merikoski, "The polynomials X^2+(Y^2+1)^2 and X^2+(Y^3+Z^3)^2 also capture their primes", Proc. LMS
    127(4) (2023) 1057-, arXiv 2112.03617.** [H] Infinitely many primes of both forms; a Deligne-type bound for
    correlations of hyper-Kloosterman sums. Relevance: none.
17. **B. Green, M. Sawhney, "Primes of the form p^2 + nq^2", arXiv 2410.04189 (v3 Jun 2026).** [H, RSS]
    * For `n = 0, 4 (mod 6)` there are infinitely many primes `p^2 + nq^2` with `p, q` prime, with an asymptotic.
      `n = 4` is the Friedlander-Iwaniec "Gaussian primes conjecture".
    * Type I/II sums in `Q(sqrt(-n))`. The Type II sums use quantitative concatenation (Kuca; Kuca-Kravitz-Leng) and
      the quasipolynomial inverse theorem (Leng-Sah-Sawhney).
    * Relevance: low.
18. **Friedlander, Iwaniec, "Coordinate distribution of Gaussian primes", JEMS 24 (2022) 737-772, arXiv
    1811.05507.** [H] Gaussian primes with one coordinate a prime square and the other an almost-prime (at most 7
    prime factors, all `> k^(1/49)`). Relevance: none.
19. **K. Pratt, "Primes from sums of two squares and missing digits", Proc. LMS, arXiv 1806.02699** [H]: an
    asymptotic for primes `m^2 + l^2` with `l` missing three fixed decimal digits.
    **C. Sabuncu, "Missing digits and sums of two prime squares", arXiv 2601.20897 (Jan 2026)** [M].
    **J. Maynard, "Primes with restricted digits", Invent. Math. (2019) 127-218** [H].
    **Heath-Brown, Li, "Prime values of a^2+p^4", Invent. Math. 208(2) (2017) 441-499** [H].
    Relevance: none.

### 1.4 Lattice points on circles: angular statistics and correlations for generic `n`

20. **E. Bombieri, J. Bourgain, "A problem on sums of two squares", IMRN 2015(11) 3343-3407.** [M]
    * Later papers attribute `\|C_(6;m)\| = O(N_m^(7/2))` (6-correlations) to it.
    * Li (Mathematika 2020) proves unconditionally, in an extremal case, a result that Bombieri-Bourgain proved
      under BSD + RH for elliptic L-functions.
    * The abstract itself was not retrieved. Relevance: low. It is already in the record as capped
      (`local_smoothing/chains.md` row E).
21. **Kurlberg, Wigman, "On probability measures arising from lattice points on circles", Math. Ann. (2016/17),
    DOI 10.1007/s00208-016-1411-4, arXiv 1501.01995.** [H] The attainable angular measures include all extreme
    points of the symmetric invariant measures and are closed under convolution; some symmetric measures are not
    attainable. Relevance: none.
22. **Kurlberg, Lester, "Poisson spacing statistics for lattice points on circles", Invent. Math. (2025)
    537-585, DOI 10.1007/s00222-025-01324-1, arXiv 2112.08522.** [H] Nearest-neighbour spacings are Poissonian
    along a density-one subsequence of admissible radii. The points are not i.i.d.-like at scales `N^(-gamma)`,
    `1/2 < gamma < log(pi/2)/log 2`. Relevance: none.
23. **A. Sartori, "Spectral quasi correlations and phase-transitions for the nodal length of Arithmetic Random
    Waves", arXiv 2005.04698.** [H] Generic sums of two squares have no spectral quasi-correlations, i.e. no small
    sums of lattice points on the circle. "Generic" is defined through a model in which the angles of the Gaussian
    primes dividing `n` are i.i.d. uniform. Relevance: low. That model *is* the fake model.
24. **Granville, Wigman, "Planck-scale mass equidistribution of toral Laplace eigenfunctions", CMP 355 (2017)
    767-802** [H] (density-one subsequence; proximity of lattice points after Cilleruelo).
    **Krishnapur, Kurlberg, Wigman, Annals 177 (2013) 699-737** [H] (background).
    **Erdős, Hall, "On the angular distribution of Gaussian integers with fixed norm", Discrete Math. 200 (1999)
    87-94** [M] (equidistribution for almost all `n`).
    Relevance: none.
25. **Parkkonen, Paulin, "On the multiplicative pair correlations of sums of two squares", arXiv 2602.13058
    (Feb 2026).** [H, RSS] Pair correlations of the logarithms of norm-form values at various scalings, passing from
    a point mass to Poissonian behaviour. Relevance: none.
26. **K. Oganesyan, "Lattice points on small arcs", arXiv 2107.09991, WITHDRAWN** ("mistake in Lemma 2") [H].
    Already in the record.

### 1.5 Diophantine and multiplicative independence; gcds

27. **Ferreira, Ribas, "Arithmetic properties of arguments of algebraic numbers on the unit circle", arXiv
    2602.23597 (Mar 2026).** [H, RSS] If `beta` is algebraic with `\|beta\| = 1` and not a root of unity, then
    `Arg(beta)/2pi` is Diophantine, via linear forms in logarithms. Relevance: none (single form).
28. **Pappalardi, Sha, Shparlinski, Stewart, "On multiplicatively dependent vectors of algebraic numbers", Trans.
    AMS 370 (2018) 6221-6244, arXiv 1606.02874.** [H] Asymptotic counts of multiplicatively dependent vectors of
    bounded height. Relevance: none.
29. **Corvaja-Zannier and Bugeaud-Corvaja-Zannier** (`gcd(a^n - 1, b^n - 1)`; `gcd(u-1, v-1)` for S-units), with
    the Silverman/Levin interpretation as Vojta on blow-ups of `G_m^2`. Background snippets only [M]. Relevance:
    none (fixed `S`).
30. **H. Gou, "Extremal problems for GCDs and LCMs in higher dimensions", arXiv 2604.21122 (Apr 2026).** [H, RSS]
    Extends Green-Walker: if `gcd >= D` on a proportion `delta` of tuples, then
    `prod \|A_i\| <<_eps delta^(-k/(k-1)-eps) prod X_i/D^k`. Relevance: none (the wrong direction for L1).

### 1.6 Poisson/Kubilius models for prime factors

31. **K. Ford, "Poisson approximation of prime divisors of shifted primes", IMRN 2025(7), DOI
    10.1093/imrn/rnaf079, arXiv 2408.03803 (v4 Jan 2026).** [H, RSS]
    * A Kubilius model for shifted primes, with a total-variation bound.
    * The prime factors in disjoint sets behave like independent Poisson variables, under a
      primes-in-progressions hypothesis Z(gamma).
    * A transference principle between random integers and random shifted primes.
    * Relevance: **shape**.
32. **O. Gorodetsky, "A Kubilius model for sieve-theoretic sequences", arXiv 2608.14190 (Aug 2026).** [H, RSS]
    * A total-variation bound in the Kubilius model for any sequence with positive level of distribution, which the
      author expects to be qualitatively optimal.
    * It recovers Ford, and Tenenbaum's optimal bound up to `x^(o(1))` and `u^(o(u))`.
    * Relevance: **shape**.

### 1.7 Combinatorial and inverse large sieves

33. **Croot, Mao, Pohoata, Sheffer, Yip, "A combinatorial large sieve for Sidon sets, distances, and norm forms",
    arXiv 2606.17487 (v2 25 Jun 2026).** [H, RSS and web]
    * "algebraic splitting modulo many small primes: local congruence branching produces many modular collisions,
      while global bounded-multiplicity hypotheses force these collisions to be rare".
    * Sidon subsets of `{1^2, ..., N^2}` have size `<= N exp(-c log N/loglog N)`, and so do repeated-distance-free
      subsets of `[N]^2`.
    * An entropic version (Shearer's inequality: a poorly distributed set modulo products of primes has an entropy
      deficit) handles `B_2[g]` sets for norm forms over number fields, `B_3[g]` in cubes and `B_4[g]` in fourth
      powers.
    * Relevance: **mechanism**.
34. **Croot, Mao, Pohoata, Yip, "A sharp inverse theorem for the quadratic large sieve", arXiv 2607.15311**
    [H, RSS]: `A ⊂ [N]`, `\|A\| >> sqrt N`, `\|A_p\| <= p/2 + O(1)` for `p << log N` force
    `exp(c sqrt(log N)/loglog N)` elements on one quadratic `m +- x^2`.
    **Croot, Mao, Yip, "An inverse theorem on sets with rich additive structure modulo primes", arXiv 2510.08862**
    [H, RSS].
    Relevance: low (L4-type inverse steps only).

---------------------------------------------------------------------------------------------------------------

## 2. Reductions: what each would use, what caps it, and the arithmetic

### 2.0 The target statements, exactly

These are taken from `round15/rigSLtheory.md` §4.1, `round14/rigtheory.md` §4.2, `round7/p3lower.md` §4 and
`round9/README.md` §3.

* **Pair integers.** `B_xy = D_(k_x - k_y)(Z_(a_x - a_y))` is a nonzero rational integer, and it equals the reduced
  chord (round 9 Thm 1.1). For odd `q` not dividing `N`, `z_x = z_y (mod q^a)` iff `q^a \| B_xy`. Exactly,
  `mu_xy = 2 log\|B_xy\| + log nu_s - 2 log t_xy`.
* **H^gen(eta, delta) and H(eta, delta).**
  * Window: `Q_eta = {q prime : M^(1-eta) < q <= M, q odd, q not dividing N}`.
  * `c_xy = #{q in Q_eta : q \| B_xy}`, `F_k = sum_(x<y) binom(c_xy, k)` and
    `E_k = binom(M,2) e_k((1/m_q)_(q in Q_eta))`.
  * H requires `F_k >= delta^k E_k` for every `k` with `delta^k E_k >= 1`.
  * Theorem C uses only `F_k >= 1` at `k = floor((2 - eps/4) log M/loglog M)`, i.e. *one heavy pair*.
* **L1(K_u, lambda_0).** `K_u`-uniform margins and `lambda_Q >= lambda_0`. For actual clusters this is equivalent to
  every `\|B_xy\| <= M^b` and `t_xy >= M^(-b)` (round 9 Thm 2.1). Equivalently, every pair gcd has norm
  `>= N^(1/2) M^(-O(1))`.
* **L7 = (I3).** Collision counts of dense coherent integers (`B_xy`, twin values) in `(M^(1-eta), M]` have Poisson
  tails. It is open even for genuine primes.
* **Payoff.** RigS (vacuous for actual clusters) + H^gen + L1 give
  `W >= ((8(1-eta) - eps)/K_u) M log^2 M/loglog M`, hence
  `M <= (K_u/(4(1-eta)) + o(1)) log R logloglog R/(loglog R)^2`. That would be growth.

**Requirements any literature input must meet.**
* **(Q1) Fake-violated.** `D_A` violates H (p3lower Lemma C.2). Robust ensemble inputs are satisfied by grafts
  (Thm B).
* **(Q2) Worst-case `N`.** Near-extremal `N` form a density-zero set: `#{N <= x : M_N(C) >= 2} << C x^(3/4) log x`
  (`outside.md` Prop. D.1, quoted in `analytic.md` §4.2). So density-one or "almost all" statements are silent.
* **(Q3) Scales.**
  * The cluster primes are `r >= (M/4)(1 - o(1))` (round 15 Thm S) of the `~ X/log X` Gaussian primes below
    `X = M^(h_0)`.
  * The window moduli are `<= M`.
  * The joint moduli needed by H are `Q_k >= M^((1-eta)k) = exp((2 - eps/4)(1-eta) log^2 M/loglog M)`.
* **(Q4) The project's constraints.**
  * No uniform exceptional sets from fixed-`S` Subspace.
  * Keep the Gaussian gcds. `B_xy` is defined after removing `g_xy`, so nothing may be discarded there.
  * Keep the drifting phase. Only pair *differences* `<2(a_x - a_y), theta>` are phase-free, which is why H is
    stated through `B_xy`.
  * No template-dependent `M_0`.

**Arithmetic of H** (`primes_checks/arith.py`; `eta = 1/2`, `eps = 0.1`, `lambda_Q ~ log 2`):

| `M` | `k` | `log E_k` | `log10 Q_k` (joint modulus) | `log10 binom(M,2)` | ratio | Kubilius `u` needed | `log10` of level `D = M^u` |
|---|---|---|---|---|---|---|---|
| 1e3 | 7 | 2.0 | 10.5 | 5.7 | 1.84 | 6.1 | 18.4 |
| 1e6 | 10 | 8.2 | 30.0 | 11.7 | 2.56 | 8.7 | 52.1 |
| 1e9 | 13 | 13.4 | 58.5 | 17.7 | 3.31 | 11.3 | 101.5 |
| 1e12 | 16 | 18.0 | 96.0 | 23.7 | 4.05 | 13.9 | 166.7 |

The ratio `log Q_k/log binom(M,2) ~ (1-eta) log M/(2 loglog M) -> infinity`. Every statement that controls the
family of `binom(M,2)` pair integers through moduli up to its own size, by pigeonhole, the large sieve or a level of
distribution, falls short by an unbounded power.

### 2.1 Ensemble angle and residue statistics (items 1-14)

**What they use.** Statistics of the genuine angle function `varpi -> arg varpi` (and residues) over *all* Gaussian
primes of norm `<= X`.

**Cap.** `round10/analytic.md` Thm B and Cor. B.1.
* Each result is robust in the sense of Def. 3.3:
  * HLR sector counts are of size `>= X^(2/5 - o(1))` (R3, `gamma = 2/5`);
  * Harman-Lewis lower bounds hold for sectors of area `X^0.619` (`gamma ~ 0.62`);
  * GRH-quality Hecke sums have error `X^(1/2)` (R2);
  * Bombieri-Vinogradov in sectors (Halupczok, Khale et al.) has error `>= X^(1/2)/(log X)^A` per modulus on
    average (the congruence version of R3).
* Hence the graft of a `lw.md` 6.1 fake with primes at height `M^A`, `A > max(2, 1/gamma)` (`A > 5/2` for HLR,
  `A > 2` for GRH and Harman-Lewis), satisfies them up to `1 + o(1)`, with `W <= (2A + o(1)) M log M`.

**Concrete attempt (HLR).** What would "almost all sectors of width `X^(-rho)`, `rho < 3/5`, have the expected
count" say about a cluster?
1. The cluster constrains no individual `theta_j`. It constrains `<2(a_x - a_y), theta> in (pi/2)Z + O(Delta)`
   with `Delta = C e^(-W/4)`, which is `M^(-c_0 M/4)` at `W = c_0 M log M`. That is far below any `X^(-rho)` with
   `X <= M^(h_0)`.
2. "Almost all" means an exceptional set of measure `o(1)`. It contains any `O(M)`-point selection, which has
   measure zero (Prop. C, §4.2).
3. The variance results (Rudnick-Waxman, Chen et al., Katz, Lester-Waxman-Yesha) are second moments over sector
   position, so they are averages (Prop. C).

**Low-height regime `P_+(N) < M^(2-delta)`.** This is the one uncapped window (`analytic.md` §6, open point 1).
There `X^(2/5) = M^(0.8 - 0.4 delta) < M`, so a robust count could see `O(M)` points in size. But (1) and (2)
still apply: no individual angle is constrained, and exceptional sets swallow the selection. No result found
constrains a selection of `~ M` primes among `~ M^(2-delta)/log M` through linear relations among their angles.

**Musicantov-Zehavi.** For `pi = a + bi`, the root `nu (mod p)` determines `pi` up to units, so `theta_p` is a
function of `nu_p`. Joint equidistribution of `(nu_p/p, theta_p)` is an ensemble statement about one map
`p -> (nu_p, theta_p)`, and it is robust. Note also that the cluster's labels are residues of `pi_j` modulo
*small* `q <= M`, not modulo `p_j`.

**BV in sectors.** This gives the distribution from which labels are drawn. But H^gen is open *given exact
knowledge of the genuine labels* (round 9: L7 open even for genuine primes), because H depends on how the profile
and windings combine the labels, not on the labels' marginal distribution.

**Verdict.** Capped (Thm B, Prop. C).

### 2.2 Gaussian primes in thin or sparse sets (items 15-19)

**What they would use.** That the cluster primes, or external witness primes, are Gaussian primes with extra
restrictions.

**Why this fails.**
* These theorems are lower bounds or asymptotics for the number of primes in a thin set `T`. A growth proof needs
  non-existence of a configuration, or worst-case statistics of `O(M)` adversarially chosen primes.
* A thin-set theorem enters only by producing *external witnesses*: Gaussian primes near prescribed lines,
  i.e. PCG/GPA (`analytic.md` Def. 5.1, Prop. D.1). Two problems:
  1. Even a proof of GPA collapses PCG to LCG, which is the problem itself (Thm D). "Fake-killing, but no
     reduction."
  2. The witnesses must exist in strips of width `h^(-1/2)` at norm `~ h`. The expected number of lattice points
     there is `O(1)`, a Bunyakovsky-type regime. Thin-set theorems need expected counts `>= X^c`. Merikoski's sets
     have density `X^(-delta)` with `delta < 1/10`.

**The one adversarial-selection result.** Merikoski's formula holds for **arbitrary** `B` with
`\|B\| = X^(1/2-delta)`, `delta < 1/10`, i.e. `\|B\| >= X^0.4` inside a range of size `X^(1/2)`. The cluster
selection has `\|S\| <= M = X^(1/h_0)`:
* `X^(0.333)` for `h_0 = 3`;
* `X^(0.2)` for `h_0 = 5`;
* `X^(0.111)` for `h_0 = 9`.

So it is always too sparse (`arith_out.txt`). Even in range, the theorem counts primes with `b in B`. It says
nothing about linear relations among their angles.

**Green-Sawhney.** The Gowers-norm Type II machinery controls linear patterns of *bounded complexity* in dense
subsets of `Z[i]`. The cluster's relations are `M - 1` linear forms in the *exponent* lattice with `r ~ M` free
prime angles. Their complexity grows like `log R/loglog R`, and they are multiplicative, not additive, in `Z[i]`.
Absolute-value counting of such tuples cannot certify emptiness anyway (Prop. C).

**Verdict.** Capped or irrelevant.

### 2.3 Angular statistics and correlations of lattice points on circles (items 20-26)

**What they would use.** Fine-scale angular structure of `{z : \|z\|^2 = n}`.

**Why this fails.**
* Kurlberg-Lester, Sartori, Granville-Wigman, Erdős-Hall and Klurman-Sartori hold for density-one (generic) `n`.
  By (Q2) they are silent on near-extremal `N`.
* Sartori's genericity is defined by i.i.d. uniform prime angles, which is the angle model of the fakes. So
  "generic" means "fake-like" there.
* Bombieri-Bourgain and the later 6-correlation bounds are upper bounds on additive energies. They hold for every
  angle vector with only the formal relations, so they hold for angle-generic fakes (`local_smoothing/chains.md`
  row E: capped; energies `<= Theta(log R)`).
* Kurlberg-Wigman's attainable measures are weak limits of the angular measures. An `M`-point cluster among
  `r(N) >= 4 * 2^r` points is invisible to them.
* "No small additive relations among lattice points" for generic `n` would rest on unit-equation counts, which are
  `c^s` with `s ~ r ~ M`. These are useless at this scale, as the task notes. Round 9 Prop. 5.2 shows that actual
  clusters *do* carry exact relations `sum c_x z_x = 0` of height `M^(c_0/2 + 1/2 + o(1))`, and that their content
  is lattice-ness.

**Verdict.** Capped (density, energy class, measure-level).

### 2.4 Diophantine and multiplicative independence (items 27-30)

* *Ferreira-Ribas.* An irrationality measure for one `arg(pi/conj pi)/2pi`: single-form Baker, capped (README §3:
  Baker/Matveev and single-form Lang-Waldschmidt).
* *Pappalardi-Sha-Shparlinski-Stewart.* Counts in height boxes. Actual clusters are multiplicatively dependent by
  construction (all are `S`-units for `S` the primes of `N`), so an average count says nothing about one `N`.
* *Corvaja-Zannier / Levin gcd bounds.* The relevant quantity is `gcd(u - 1, v - 1)` with
  `u = z_x/z_y = unit * prod (pi_j/conj pi_j)^(a_xj - a_yj)`, an `S`-unit. Two obstacles:
  * The bounds hold for fixed `S` with exceptional sets depending on `S`, which the project constraints forbid.
  * Their uniform form is Vojta on a blow-up, already the conditional route (`growth/conditional.md`).
* *Gou / Green-Walker.* These give upper bounds on sets with many large gcds. L1 needs a *lower* bound on all pair
  gcds of a cluster. The Gaussian analogue would give `\|A\| << N/D`, trivial at `D = N^(1/2) M^(-b)`.

**Verdict.** Capped, forbidden, or the wrong direction.

### 2.5 Poisson/Kubilius models (items 31-32): the right shape, impossible strength

**Property used.** That the pair integers `B_xy` of an actual cluster have Poisson prime-factor statistics in the
window. This is exactly H. It is fake-violated: `D_A` has `c_xy <= floor(18/(1-eta))`, so `F_k = 0` for large `k`.

**Attempted reduction.**
1. Treat `A = {B_xy}` (`binom(M,2)` integers) as the "sequence".
2. Theorem C needs `F_k >= 1` at `k = floor((2 - eps/4) L/l)`. Its model probability is
   `pi_k = P(Poisson(lambda) >= k) ~ lambda^k/k! = M^(-2 + eps/4 + o(1))`.
3. A Kubilius/Poisson transfer gives `#{pairs with c >= k} >= binom(M,2)(pi_k - TV)`. So it needs `TV < pi_k`.
4. The known bounds (Tenenbaum; Ford; Gorodetsky for level of distribution `D`) have the shape
   `TV << exp(-c u log u) + (level error)` with `u = log D/log y`, where `y = M` is the largest prime used.
5. `TV < pi_k` needs `u log u >= log(1/pi_k)`, i.e. `u = 6.1, 8.7, 11.3` at `M = 10^3, 10^6, 10^9`. That is a
   level `D = M^u` of `10^18`, `10^52`, `10^101`, for a family of `10^5.7`, `10^11.7`, `10^17.7` elements.
6. A finite family cannot have a meaningful level of distribution above its size: for `d > \|A\|` the counts
   `#{a : d \| a}` are 0 or 1 against expectations `< 1`.
7. Worse, `A` is adversarial. "Level of distribution for the pair integers of every near-extremal cluster" is
   exactly the open L7/H.
8. Ford's hypothesis Z(gamma) concerns *shifted primes*, an ensemble with a natural level of distribution. The
   cluster family has none.

**Verdict.** The statement has the right form, but no route reaches the needed level. Its failure is a power that
grows without bound (§2.0 table).

### 2.6 The combinatorial large sieve (items 33-34): the same mechanism as the proved bound

**The match.** In `two_thirds/proof.md`:
* *Local branching.* The `M` points lie on the conic `x^2 + y^2 = N (mod p)`, which has `p - 1` points for split
  `p` not dividing `N`. So at least `E(M, p - 1)` pairs collide mod `p`.
* *Collisions to divisibility.* Each collision forces `p \| B_xy`.
* *Global bound.* From `mu_xy = 2 log\|B_xy\| + log nu_s - 2 log t_xy` and the cut identity
  `sum_(x<y) mu_xy <= M W/4`, we get `sum_(x<y) log\|B_xy\| <= M W/8 + O(M^2 log^+ C)`.
* Comparing this with `sum_(p <= y) E(M, p-1) log p ~ (M^2/2) log y`, with `y ~ M`, gives `W >> M log M`. The
  sharp form, with inert primes, levels and the chain lemma, is `3 M log M <= W + K_C M`. This is the rate
  `log R/loglog R`.

This is precisely "local congruence branching produces many modular collisions, while global bounded-multiplicity
hypotheses force these collisions to be rare". Croot et al.'s saving `exp(-c log N/loglog N)` comes from the same
count: primes up to `~ log N`, and at most `log N/loglog N` of them dividing a nonzero integer `<= N`. That is the
rate the project already has.

**The entropic version (Shearer).** It works with residues modulo *products* of primes, i.e. the `k >= 2` regime
where pigeonhole fails. But by CRT those joint residues are functions of the per-prime residue data, which members
of P3-char and P3-all carry at moduli `<= 3M^A`.
* Croot et al. use primes `<= log N` or `<= (log N)^3`. For our clusters that is `(c_0 M log M)^3 <= M^(3+o(1))`,
  i.e. fixed `A`.
* So any entropic inequality they prove is a residue-class input at fixed `A`. It is capped by the round 5-7
  residue fakes with `p_j > X`. The round-9 R1 caveat (caps only at fixed `A`) does not rescue it, because their
  `A` is fixed.

**Inverse sieves (item 34).**
* They would serve L4-type inverse steps: labels with low demand are structured. Round 10 §7.1 shows that label
  lemmas of L4 type have no fake-violated input.
* Quantitatively they need `\|A\| >> sqrt(height)`. The cluster's `r <= M` primes have height `M^(h_0)`, and
  `M << M^(h_0/2)` for `h_0 > 2`.

**Verdict.** The same mechanism as the proved bound, capped at the same rate. At most a constant could change (out
of scope).

---------------------------------------------------------------------------------------------------------------

## 3. Papers to request in full text (ranked)

None is expected to beat the barrier. The ranking is by how much a full reading could change the classification
above.

1. **Croot-Mao-Pohoata-Sheffer-Yip, arXiv 2606.17487.** The closest mechanism. Check:
   * whether the entropic (Shearer) argument for norm forms uses anything beyond per-prime residues plus one global
     multiplicity bound, for instance the full factorisation of the norm form jointly across primes in a way that
     P3-all members cannot imitate;
   * which primes it uses: the range `<= log N` or `<= (log N)^3` decides whether it sits inside the
     fixed-`A` caps;
   * a by-product could be a better constant than 2/3, not a better rate.
2. **Gorodetsky, arXiv 2608.14190, and Ford, arXiv 2408.03803 (IMRN 2025).** Confirm the exact TV bound shape. Ask
   whether any version gives lower bounds for the *upper tail* `P(omega_window >= k)` at `k ~ log x/loglog x` with
   a weaker input than a level of distribution. That would be the only way the §2.5 estimate could fail.
3. **Bombieri-Bourgain, IMRN 2015(11) 3343-3407.** Check which statements hold for all `n` (worst case) and which
   only for generic `n`, and whether any result on linear relations among lattice points is uniform in `n`.
   Compare with round 9 Prop. 5.2 (relations of height `M^(c_0/2 + 1/2 + o(1))` exist).
4. **Merikoski, Compositio 2025 (arXiv 2302.11331).** Remark 1.4 and the quasi-explicit formula uniform over
   arbitrary `B`. Check whether the method could reach `\|B\| = X^(o(1))`, and whether it controls angles jointly.
5. **Sartori, arXiv 2005.04698, and Kurlberg-Lester, arXiv 2112.08522.** Only to confirm that "generic" means
   density one, or a full-measure set for the i.i.d.-angle model. Then they are silent on near-extremal `N`.
6. **Rudnick-Waxman, arXiv 1705.07498.** Lemma 2.1, the angle-gap lower bound. Expected to be the lattice fact
   `\|Im(varpi conj alpha)\| >= 1`, i.e. `analytic.md` Thm D.

---------------------------------------------------------------------------------------------------------------

## 4. Search log

**WebSearch** (S = standard, E = extended). Each entry gives the query, then what it confirmed.

1. S `Gaussian primes in narrow sectors almost all sectors 2023 arXiv` → Järviniemi-Teräväinen 2303.05822
   (Thms 1.1/1.2, widths); HLR 1903.04005 (area `X^(2/5+eps)`).
2. S `Rudnick Waxman angles of Gaussian primes variance small sectors` → Rudnick-Waxman 1705.07498; Chen et al.
   1901.07386 (Exp. Math. 2023).
3. S `Green Sawhney primes of the form p^2 + nq^2` → 2410.04189, abstract, v3 June 2026; Xiao extension (CMS
   abstract).
4. S `Bombieri Bourgain "A problem on sums of two squares" linear relations lattice points circle` → IMRN
   2015(11) 3343-3407; IAS seminar abstract (toy model for nodal curves).
5. E `Bombieri Bourgain sums of two squares "6-correlations" ...` → `\|C_(6;m)\| = O(N_m^(7/2))` attributed to BB;
   Benatar-Marinucci-Wigman.
6. S `Merikoski Gaussian primes thin sets sectors 2022 2023` → Compositio 2302.11331 (abstract).
7. S `Kurlberg Wigman lattice points on circles linear relations "additive" correlations Nazarov-Sodin` →
   Math. Ann. paper; Chatzakos-Kurlberg-Lester-Wigman (hyperbolic); Bordignon-Kurlberg 2607.13659.
8. S `Duke Friedlander Iwaniec Hecke Grossencharacter large sieve angles Gaussian primes joint distribution` →
   Devin 2105.02492; FI 1811.05507 listing.
9. S `"Discrepancies in the distribution of Gaussian primes" abstract` → Devin (Chebyshev bias, Hecke characters).
10. S `Friedlander Iwaniec "Coordinate distribution of Gaussian primes"` → JEMS 24 (2022) 737-772; Almost Primes
    Theorem.
11. S `Gaussian primes in sectors with restricted digits Maynard ...` → Merikoski Remark 1.4 (narrow sector);
    Maynard digits.
12. S `Pratt primes sums of two squares missing digits` → Pratt 1806.02699; Sabuncu 2601.20897.
13. E `Cilleruelo Granville conjecture lattice points short arcs ...` → Oganesyan 2107.09991; Temur 2012.10784;
    CG CJM 2009.
14. S `Katz "Witt vectors and a question of Rudnick and Waxman" ...` → IMRN 2016, DOI rnw130; super-even
    characters.
15. S `Gaussian primes with small imaginary part Ankeny problem Harman Lewis 2020s improvement` → nothing new;
    Ankeny/Kubilius under GRH; Harman-Lewis.
16. S `lattice points on circles equidistribution residue classes modulo q ...` → Cilleruelo JAustMS 2002;
    Musicantov-Zehavi 2112.07494; EMV (spheres).
17. E `Oganesyan "Lattice points on small arcs" arXiv 2107.09991 abstract` → **withdrawn** (Lemma 2 error).
18. S `"Gaussian primes" angles "linear relations" OR "linearly independent" ...` → qualitative independence of
    `pi` and the prime angles (math-ph/9812019). No quantitative independence measure.
19. S `Sartori "spectral quasi-correlations" ...` → 2005.04698 (generic `n`, i.i.d.-angle model).
20. S `Kurlberg Lester Poisson spacing statistics ... Inventiones 2025` → Invent. Math. 2025, 537-585.
21. S `"Sectorial equidistribution of the roots of" x^2+1 modulo primes abstract` → Musicantov-Zehavi, QJM
    75(2) 2024.
22. S `least Gaussian prime in a narrow sector Thorner Zaman effective bound` → nothing sector-specific; Matsui
    reference.
23. S `pair correlation angles of Gaussian primes Montgomery Hecke L-functions zeros statistics` → nothing
    specific to angles.
24. S `Gaussian primes short intervals and narrow sectors simultaneously polar boxes Harman 2019 2021` → Stucky
    (Ricci generalisation); Harman-Kumchev-Lewis `\|z\|^0.53`.
25. S `large sieve inequality Hecke Grössencharacters imaginary quadratic field angles 2020 2022 improvement` →
    only explicit constants (Grześkowiak 2020) and power moduli (Gao-Zhao 2022).
26. S `"Discrete fractional integrals, lattice points on short arcs, and diophantine approximation" abstract` →
    Temur, Turkish J. Math. 46(7) (2022).
27. S `Klurman Sartori lattice points small arcs density one ...` → not retrieved; CG `R^(1/2 - 1/(4[k/2]+2))`
    (math/0608109).
28. E `Bombieri Bourgain "sums of two squares" IMRN 2015 abstract ...` → BB results conditional on BSD + RH, made
    unconditional in an extremal case by Li (Mathematika 2020).
29. S `Bombieri-Vinogradov theorem Gaussian primes residue classes sectors ...` → Halupczok 1607.07265;
    Khale et al. 2008.09677; Lam 1812.08823.
30. S `number of distinct prime factors of pairwise differences of a set of integers ...` → nothing relevant.
31. S `Erdős–Kac theorem Lehmer sequences Lucas numbers ...` → Stewart, Acta Math. 2013 (largest prime factor);
    Wagner (Erdős-Kac for the primitive-part count). Nothing on window counts.
32. S `greatest common divisor u-1 v-1 S-units Corvaja Zannier uniform bound ... Levin Vojta blowup` → background
    only.
33. S `"multiplicatively dependent vectors" ...` → Pappalardi-Sha-Shparlinski-Stewart, TAMS 2018.
34. S `Kurlberg Wigman "On probability measures arising from lattice points on circles" abstract ...` → abstract.
35. S `Duke Friedlander Iwaniec "Equidistribution of roots of a quadratic congruence to prime moduli" Annals` →
    DFI 1995; Toth 2012; Grimmelt-Merikoski 2025 (uniformity in `h`).
36. S `Heath-Brown Li "Prime values of a^2+p^4" Inventiones` → Invent. Math. 208(2) 441-499.
37. S `Merikoski "The polynomials X^2+(Y^2+1)^2 and X^2+(Y^3+Z^3)^2 also capture their primes"` → Proc. LMS
    127(4) (2023).
38. S `Thorner Zaman log-free zero density Hecke L-functions least prime Chebotarev ...` → 1510.08086; BTZ JEMS
    2022.
39. S `Maynard "Primes with restricted digits" Inventiones 2019 abstract` → Invent. Math. 2019, 127-218.
40. S `lower bounds linear forms in many logarithms small coefficients improvement 2022 2024 Lang-Waldschmidt ...`
    → nothing new; the LW conjecture as stated by Stewart.
41. S `Huang Liu Rudnick "Gaussian primes in almost all narrow sectors" abstract Acta Arithmetica X^(2/5)` →
    Theorem 1 (`rho < 3/5`), Acta Arith. 193 (2020) 183-192.
42. E `Croot Mao Pohoata Sheffer Yip combinatorial large sieve Sidon sets squares distances norm forms` → 2606.17487
    (v2), mechanism description (entropy, Shearer).
43. S `Krishnapur Kurlberg Wigman nodal length fluctuations ...` → Annals 177 (2013) 699-737.
44. S `"lattice points" circle "no non-trivial" linear relations ... Bourgain Rudnick` → BR `O(log m)` on arcs
    `<= m^(1/4)`; separation for all but `O(X^(1-eps/3))` `m`.
45. S `Granville Wigman "Planck-scale mass equidistribution ..."` → CMP 355 (2017).
46. S `Erdős Hall "angular distribution of Gaussian integers with fixed norm" result` → Discrete Math. 200 (1999);
    Dias 1404.6271.
47. S `bounded gaps between Gaussian primes narrow sectors Maynard Tao ...` → nothing beyond items 1, 2, 15.
48. S `Ford "Poisson approximation of prime divisors of shifted primes" ...` → IMRN 2025(7); Hypothesis Z(gamma).

**Fetch attempts.**
* `emergentmind.com`, `pith.science`: DNS failure.
* `arxiv.org`, `export.arxiv.org`, `ar5iv`, `semanticscholar`: 403.

**Local RSS mirror (regular expressions over titles and abstracts, 2024-03 to 2026-10).**
* Patterns: `gaussian prime`; `narrow sector|short sector|angles? of (gaussian )?prime|Hecke (character|Gr|L)|
  ossencharak|sums of two squares|lattice points on (the )?circle|x^2+y^2|Gaussian integer|Z[i]|Q(i)`;
  `Kubilius`; `roots of (quadratic|x^2)`; `large sieve`; `linear forms in logarithms`;
  `Lucas|Lehmer (sequence|number)`; `primitive divisor`; `Erdős-Kac`; `inverse (large )?sieve`;
  `number of prime factors`; `pairwise differences`.
* New hits used: 2605.19105, 2609.23214, 2602.13058, 2602.23597, 2408.03803, 2608.14190, 2606.17487, 2607.15311,
  2510.08862, 2505.00493, 2604.21122, 2412.14792, 2410.04189, 2302.11331.
* Checked and irrelevant: 2609.13894 (Linnik microsquare), 2610.18281 (elliptic sieve), 2309.10971 (annuli),
  2606.29185 (van der Corput for sums of two squares), 2504.09650/2508.18044 (Diophantine approximation with sums of
  two squares), 2610.34556 (admissible sets in `Z[i]`), 2512.04081 (additive relations in irrational powers),
  2604.20014/2605.24909 (Lucas rank densities), 2610.20554 (large sieve and Ostmann for sums of two squares).
