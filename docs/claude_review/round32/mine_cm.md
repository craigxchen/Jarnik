# Round 32, block cm: CM theta values, their jets, and divisibility at inert primes

Sources: S1 = the 4-page remark (timeline items 2, 9, 18); S2 = the downloaded Lean tree (not compiled, untrusted; only
statements were read, no proofs). Labels: **[Read]** = a statement read in S2 (path given); **[Record]** = in the
review record (file given); **[Inferred]**, **[Heuristic]**, **[Speculative]** as usual. About 30 file reads were used.

## 0. Verdict

1. The timeline's exponent "a/p + a/p^3 + ... - O(log a)" is, in S2, the constant
   `gam l = l * rho l = l/(l^2-1) = 1/l + 1/l^3 + ...` (with `rho l = 1/(l^2-1)`) in the jet bounds of
   `Formalism/GaussJetBounds` [Read]. These are Dwork-style bounds for exponentials of theta-type power series
   written in the logarithm coordinate of a **height-2 (supersingular) Lubin-Tate/Honda formal group**,
   `log(x) = sum_e x^{l^{2e}}/l^e`. The CM setting is a **general imaginary quadratic field** (discriminant `-d`,
   primes `p >= 5` with `(-d/p) = -1`), not specifically Q(i) or Q(omega) [Read].
2. On our side the same number `l/(l^2-1)` already appears, exactly, as the per-pair collision density at an inert
   prime summed over all levels: `sum_{a>=1} 1/((l+1) l^{a-1}) = l/(l^2-1)`. The record's proof of Theorem 1.1
   already counts all levels (`two_thirds/proof.md` Lemma 3.3(I), lines 223-225, and the bound at line 268), and
   states that levels `a >= 2` add only `O(M^2)` (line 656) [Record]. Higher inert levels therefore change
   `W >= 3M log M - O(M)` only in the `O(M)` term. The extra constant is
   `sum_{l = 3 (4)} log l/(l^2-1) = 0.2287...` (computed) times `M^2/2`.
3. "Multiplicity at an inert prime" for points on the circle is a confluent Vandermonde on a one-dimensional torus,
   and that determinant is an exact power of the ordinary Vandermonde times factorials. So in dimension one,
   multiplicity yields nothing beyond the valuations of `z_i - z_j` [Inferred; a short proof is given in §3.2]. The
   S2 gain comes from (i) a height-2 formal group, which the norm-one torus does not have: at inert `l` the torus is
   anisotropic but its formal group has height 1; and (ii) determinants on `E^n` with theta functions of growing level.
   Their counterpart on our side is route 5, the multi-prime torus, whose zero estimate is open (expert note §5 item 3).
4. Statements about `a_{4k}(N) = sum_{|z|^2 = N} z^{4k}` (congruences, Hurwitz/Kummer-type, Chowla-Selberg
   normalisation, Siegel mass) concern the whole circle or averages over `N`. They do not constrain a single short
   arc (§4(c)). The Siegel-Weil/mass side is trivial for `Z[i]` (one class in the genus), and arc counts are cuspidal,
   so the Eisenstein part does not see them.
5. **Nothing in this block transfers to a rate improvement.** Honest probability that any idea here leads to
   `M = o(log R/loglog R)`, even on a restricted class: below 1%. The block's value is a clean dictionary entry: S2's
   "inert multiplicity" is our Lemma 3.3(I) with all levels, already in the record.

## 1. What S2 actually states

### 1.1 The jet bounds with exponent l/(l^2-1): `Formalism/GaussJetBounds/Statement.lean` (namespace `Formalism.DworkPowerSeriesWeights`)

Definitions [Read]:
* `rho l = 1/(l^2-1)`, `gam l = l * rho l` (file `Formalism/FrobeniusBaseCase/Statement.lean`, lines 15-17).
* `letterLog l x` / `oneLog l x`: the power series `sum_{e>=0} X^{l^{2e}} / l^e` (in the `d`-letter version the
  `e`-th term moves the letter by `-2e` in `ZMod d`). This is the logarithm of a height-2 formal group of Honda type
  over an unramified extension, i.e. the formal group of a supersingular elliptic curve, as at a prime inert in the
  CM field.
* `BHyp l d K0 zeta be`: `be` is a power series in `d` letters with coefficients in an unramified subfield `K0`. It is
  invariant under the torus action `X_s -> zeta^{+-l^s} X_s` with `zeta` a primitive `(l^d - 1)`-th root of unity, it
  has no monomials of low bidegree, and the exponentials of its two "biextension differences"
  `de1 be = be(x+y) - be(x) - be(y)` (two variants), transported to the log coordinate, are integral with constant
  term 1. Inferred meaning: `be` is the logarithm of a normalised theta function (the cubical/biextension structure
  of a symmetric line bundle), written in formal-group coordinates.

Main statements [Read]:
* `LemmaR` (`Statement.lean` lines 138-156). Hypotheses: `K` ultrametric of characteristic 0, `l` an odd prime with
  `|l| < 1`, `d >= 2` even, `BHyp`, an integral matrix `M` mapping points to variables, `g` the formal inverse of
  `entLog` (the formal exponential), `U` integral, `epsilon = +-1`, and weights `sig_i >= 0` such that every "pure
  bi-power" monomial `X_i^{l^alpha} X_j^{l^beta}` whose coefficient fails the weight-`1/(l-1)` test satisfies
  `l^alpha sig_i + l^beta sig_j >= 1`. Conclusion:
  `|kappa! * coeff_kappa( U(g(M X)) * exp(epsilon * sum be(M X)) )| <= |l|^( sum_i kappa_i (gam l - sig_i) - |iota| (log_l |kappa| + 1) )`.
  In words: the order-`kappa` jets of the theta product are divisible by `l^{gam(l)|kappa| - (loss) - O(log|kappa|)}`.
* `OneLetterRadius` (lines 171-173): the optimal loss in one letter is `sig = 1/(l+1)`, and
  `gam l - 1/(l+1) = rho l`.
* `OneLetterR` (lines 179-199). This is the CM case `d = 2`, `zeta` of order `l^2 - 1`, i.e. `Z_{l^2}` = completion
  of `O_K` at an inert `l`. The jets of `U(g(M X)) * exp((sgn/2) sum_x be(M_x X, ep M_x X))` satisfy
  `|kappa! coeff_kappa| <= |l|^( rho(l) |kappa| - |iota| (log_l|kappa| + 1) )`.
  So the net divisibility of the jets is `|kappa|/(l^2-1) - O(log)`. The timeline's `a/p + a/p^3 + ...` is the gross
  exponent `gam` before the quadratic-term loss `1/(l+1)` is subtracted.
* `ThmAGauss` (`StatementP.lean` lines 33-51): the determinant form. The coefficients of `frameJet`, a product over
  letters of minors of the matrix of jet pairings of `thetaZ = U(g X) * exp(epsilon * bSum be)`, have norm at most
  `|l|^( sum_{s in TT} (topTerm - gam * Wcap) + Fnum - d N0 N1 (log_l a + 1) )`, with
  `Fnum = gam * sum D_s - maxFlow` (`Formalism/IncomingBounds/IncomingBound.lean` line 28). This is the
  "Stepanov determinant" divisibility at an inert place.
* `ThmBGauss` (`StatementY.lean`): the companion statement for a Gaussian kernel `exp(-S sum u_t X X')` with
  `|S| <= |p|^{-1/2}` (half-integral valuations, Gauss-sum normalisation). The lower bound for the minors is
  `(1/2) sum (size - ...)`.

### 1.2 Which theta function, which curve, which derivation

* `th0 W x0 y0 s2 = exp(-lam(...))` (`Formalism/LocalToGlobal/TorsionIntegrality/RamifiedCurveKernels/StatementT.lean`
  line 18) [Read]: the theta function of a Weierstrass curve `W` at a 2-torsion point (`2y0 + a1 x0 + a3 = 0`),
  written in the formal parameter. `ThetaZeroSquare` gives `X^2 th0^2 = sigma*(s2)^2 * (wp shifted at x0)`: up to
  the quasi-period normalisation `s2`, `th0^2` is `sigma^2 (wp - e_i)`, so `th0` is the classical theta function with
  a half-period characteristic. `ThetaZeroEven`: it is even with constant term 1.
* `ThetaZeroJetsCurve` (same file, lines 53-63) [Read]. Hypotheses: `W` integral; an endomorphism `al` with
  `al * albar = #(ker al \ 0) + 1` (its degree), `al(al - albar)` a unit, and the normalising constant fixed by the
  kernel sum `al(al - albar) s2 = sum_{P in ker al, P != 0} (x(P) + b2/12)`. Conclusion: the jets
  `n! coeff_n` of `th0` in the log coordinate (`JetsIn1`) are integral, and the invariant-derivative jets of the
  associated `phT` are integral. This is a Hurwitz-type integrality theorem for CM theta functions, normalised by
  the quasi-period, which is in turn fixed through an isogeny sum.
* The derivation is the invariant derivation `D = (log')^{-1} d/dX` of the formal group. In S2, `jet l f i` is the
  `i`-th `D`-derivative at 0, and `JetsOfDJets` / `IntDJets` (`StatementH.lean`) say that `D`-jets are integral iff
  `n!`-jets in the log coordinate are integral [Read].
* The global CM setting, `Formalism/LocalToGlobal/InertKernelTheory/StatementInt.lean`, `StatementJets.lean`,
  `StatementCurve.lean` [Read]. For `D : D7` with host imaginary quadratic field `K` of discriminant `-d`, an ideal
  `cc` and `delta` a normalised twelfth root of the discriminant (`IsDiscriminantTwelfthRootOfIdealLattice`, i.e. the
  period), and a prime `p >= 5` with `jacobiSym(-d, p) = -1`: `ModelExists` gives an integral Weierstrass model over
  `C_p`; `PhiJetInt` says the iterated invariant derivatives `dYu^n` of the theta-type series `ePhi1` at 0 are
  `p`-integral; `JetsInert` says the series in the "comb" class built from admitted Hermitian lattices `L, L'` over
  a CM field `E` have integral jets. `CurveTie` identifies the period lattice `delta * cc` with the CM lattice
  `(a0, b0, d)` up to an integer `g`.
* So the field is a **general imaginary quadratic** `K`. The curve is `C/(Omega a)` for ideals `a`, and the inert
  primes are `p >= 5` [Read].

### 1.3 Chowla-Selberg, other named modules

* `GlobalTheory.T.classMean_log_norm_disc_chowlaSelberg`
  (`Development/AbelianVarietyTheory/Proofs/DivisorProofs/ClassMeanLogNormDiscChowlaSelbergStatement.lean`) [Read]:
  for imaginary quadratic `K` of discriminant `-d < -4`, the class-group average of
  `log( N(a_c)^6 |Delta(L_c)| )` equals `12 log 2pi - 6 log d + 6 gamma_E - 6 (L'/L)(1, chi_{-d})` (the last term
  is `dlogL D 1`). This is the averaged Chowla-Selberg formula, i.e. the Faltings height of CM elliptic curves.
  `GlobalTheory.T.chowla_selberg_period` (`Development/AutomorphicSpectrum/Proofs/GeometryProofs/ChowlaSelbergPeriodStatement.lean`)
  states a period of a curve `Curve0` equal to `gamma_E + log(2pi/sqrt d) - L'/L(1, chi)`. Use (Inferred): the
  archimedean normalisation of `Omega` in `D^a theta(0)/Omega^a`, entering the product formula for the
  determinants. Timeline item 9's "proven divisibility via Chowla-Selberg" fits this reading.
* `Formalism/StatementCM.lean` [Read] is **not** about CM elliptic curves. `IsCMType Phi` means `Phi` is a half-system
  of `(Z/l)^*` (`g in Phi <-> -g notin Phi`). `TwoTypes`/`HalfTypes` are combinatorics of such types, and
  `ForcedSumCM` bounds `b Sz(a/b)` minus the von Mangoldt sums `sum Lambda(n)(a/n - b)` over the residue classes
  counted by `Phi` and its dual: `0 <= R <= 3a/l log(a/b)(1 + log(a/b)) + a log l/(l-1)`. The name does not match the
  block's theme.
* `Formalism/GaussValuationBounds/Statement.lean` (namespace `Formalism.NOK`) [Read] concerns Iwasawa / Wach-module
  valuations, not CM. Example, `DelVal`: for `chi in Z_l^*` with `v(chi - 1) = 1`,
  `v( prod_{y<=r} (1 - chi^{(l-1)y}) ) = r + v_l(r!)`. `GapPure1/2`: if `T` is annihilated on `x` by
  `prod_{j in R}(X - c_j)` over two separated ranges, then `l^E x = 0` with
  `E ~ (3a + 2S + 2eta)/(l-1) + 2 log_l(...)`. These are the "Iwasawa tower as pigeonhole on F_{q^n}" ingredients
  (timeline item 11).
* `Formalism/LatticeTheta/Statement.lean` [Read]: Poisson summation, theta inversion and positivity for Euclidean
  lattices (standard). `ReciprocityTheory/ThetaQuasiPeriods/ThetaQuasiPeriodStatement.lean` [Read]: Riemann theta
  with characteristics, quasi-periodicity, holomorphy (standard).
* `TowerValuationBounds` [Read, partially]: valuation bookkeeping along product formulas, nothing CM-specific.
  `PadicSubsetSums` [listed only]: its files concern Hecke-algebra central idempotents and `GL_n` counts, and the module
  name does not describe them. `EllipticTorsion`, `FormalGroupLog`, `FormalGroupWitt`, `ThetaLogCoordinates`
  (292 files), `ThetaSupremumBounds`, `CMTypeGalois`, `GroundTheory/HeightsAndTheta` (6410 files): listed, not read
  beyond file names. **Not found**: any statement on Hurwitz numbers, Kummer congruences for Hecke L-values, or power
  sums `sum z^{4k}` over lattice points of a fixed norm (grep for "Hurwitz number", "Kummer congruence" returned
  nothing).

### 1.4 How it is used (Inferred from the import structure, not verified)

The pattern is an interpolation (Stepanov) determinant whose entries are jet pairings of theta functions on `E^n`
at CM/torsion data:
* an inert-place lower bound for its divisibility (`ThmAGauss`, `ThmBGauss`, fed into `LocalToGlobal/...` and
  `SharpRestriction/Divisibility/JetPolynomialLogs`);
* integrality of jets (`ThetaZeroJetsCurve`, `JetsInert`);
* an archimedean upper bound (`ThetaSupremumBounds`) with periods normalised by Chowla-Selberg;
* the product formula.

The `p`-adic jet divisibility plays the role that the factorial `a!` (multiplicity `a`) plays in function-field
Stepanov. At an inert prime the supersingular formal group gives `gam(l) = l/(l^2-1)` per order instead of
`1/(l-1)`.

## 2. Our side: the closest analogue, stated precisely

Setting: `N = R^2`, a normalised cluster `z_1, ..., z_M` (Gaussian gcd removed), all on `|z|^2 = N'`, on an arc of
angular width `delta ~ C N'^{-1/4}`. Put `u_i = z_i/zbar_i = z_i^2/N'`, which lies on the norm-one torus
`T = R^1_{Q(i)/Q} G_m`.

**Proposition A (inert levels; [Record] + [Inferred]).** Let `l = 3 mod 4`, `l` not dividing `N'`. Modulo `l^a`
the points lie in the norm-`N'` fibre, which has `(l+1) l^{a-1}` classes, and `l^a | (z_i - z_j)` iff
`z_i = z_j mod l^a`. Hence
`sum_{i<j} v_l(z_i - z_j) >= sum_{a>=1} E(M, (l+1) l^{a-1})`, where `E(M, m) >= M^2/(2m) - M/2` (positive part).
The leading term is `(M^2/2) * l/(l^2-1) = (M^2/2)(1/l + 1/l^3 + ...)`. This is the Bhargava-factorial type count
of the circle at `l`: `v_l(k!_{S_N}) = sum_a floor(k/((l+1) l^{a-1})) = k * l/(l^2-1) - O(log k)`, which is
literally the timeline's "a/p + a/p^3 + ... - O(log a)".

Explicit count. Summing over inert `l <= M`:
`sum_{i<j} log|c_ij|` (inert part) `>= (M^2/2) [ sum_{l<=M, l=3(4)} log l/(l+1) + sum_{l=3(4)} log l/(l^2-1) ] - O(M^2/log M * ...)`.
The first sum is `(1/2) log M + O(1)` (the level-1 term used in Theorems 1.1/1.2). The second is the whole gain of
levels `a >= 2`; it converges: `sum_{l=3(4), l<10^6} log l/(l^2-1) = 0.2287...` (computed). Since Theorem 1.1's
inequality is quadratic in `M` on both sides, an `O(M^2)` gain becomes an `O(M)` change in
`W >= 3M log M - O(M)`. **No change of rate.** The record says exactly this: `two_thirds/proof.md` line 656,
"deeper prime powers add only O(M^2) in total".

**Proposition B (multiplicity in dimension one collapses; [Inferred], standard).** Let `x_1, ..., x_M` be distinct
in a field and `a >= 1`. The confluent Vandermonde `det[ (d/dx)^k x^m |_{x = x_i} ]`
(`0 <= m < aM`; `0 <= k < a`; `1 <= i <= M`) equals, up to sign,
`prod_{i<j} (x_j - x_i)^{a^2} * prod_i prod_{k<a} k!`. Consequently any interpolation determinant on the
one-dimensional torus `T` (or on the line) with multiplicity `a` at the `M` cluster points is nonzero for free, and
its product-formula inequality is the `a^2`-th power of the ordinary Vandermonde inequality: the cut identity of
expert note §1. The factorials cancel between places. So `p`-adic or archimedean multiplicity on `T` gives nothing
beyond the valuations of `z_i - z_j`. (With `D = u d/du` in place of `d/du`, the extra factors are powers of the
`u_i`, which are units at every `l` not dividing `N`.)

**Why S2's mechanism does not have a counterpart here [Inferred].** S2's gain needs (i) the height-2 formal group,
where `D^p = 0 mod p` because the Hasse invariant vanishes at a supersingular prime, and (ii) determinants of theta
functions on `E^n` with many independent variables. Our points live on `T`. At inert `l`, `T` is the non-split
torus: compact `T(Z_l)`, but its formal group is a twisted form of the multiplicative group, of height 1, so
`D^l = D mod l`, not 0. Dimension one gives Proposition B. The only higher-dimensional object on our side is the
multi-prime torus `T_M` (monomials `prod pi_p^{a_p}`), which is route 5. The ingredient missing there is the zero
estimate at the trivial subgroup (ZE^all_H, expert note §5 item 3), not divisibility: forced divisibility is
already "exactly critical" (Lemma K1).

**Programme C (CM elliptic curve instead of the torus; [Speculative]).** The points of `|z|^2 = N` are the
endomorphisms of degree `N` of `E = C/Z[i]`. A cluster is a set of `M` such endomorphisms whose pairwise differences
have degree `|z_i - z_j|^2 <= C^2 sqrt N`. S2's `ThetaZeroJetsCurve` uses an endomorphism `al` through its kernel
sum. For `Z[i]` the `i`-symmetry forces `s2 = 0` (the non-holomorphic `G_2^*` of `Z[i]` vanishes), so the identity
reads `sum_{P in ker al, P != 0} (x(P) + b2/12) = 0` for every `al`. One would want theta-jet determinants with
"points" built from `ker z_i` (cyclic subgroups of `E[N]`). The arc condition is archimedean on the generators.
Through the kernels it says only that `z_i - z_j` has small degree, which is the chord again. I found no way to make
the theta jets see the arc. Probability that this leads anywhere: below 0.5%.

## 3. Barrier tests

1. **Fakes.** Proposition A uses only point residues modulo `l^a` with `(l+1) l^{a-1} <= M^2` (beyond that, pigeonhole
   forces nothing). These are point-level residues at moduli `<= M^2`, which residue designs satisfy up to `M^A`
   (expert note §1, lines 131-133, with the scope limits there). Fakes satisfy it, so it is capped. Proposition B is a
   reformulation of the cut identity, which fakes satisfy.
2. **Characteristic p.** Over `F_p[t]` with `i` adjoined, inert primes are the irreducibles that stay inert. The norm
   fibre at an inert prime of degree `f` has `(q^f + 1) q^{f(a-1)}` classes, and the same pigeonhole holds verbatim.
   So it cannot give the uniform bound, which fails there. It is a local input, so it does not touch the rate either.
3. **Single-place Liouville.** Proposition B shows that multiplicity at one place in dimension one is the
   Vandermonde, i.e. the pair/cut inequality: exactly critical.
4. **Already in the record.** Proposition A: `two_thirds/proof.md` Lemma 3.3(I) and lines 268, 380, 392, 656.
   The level-1 version: Theorem 1.2 (expert note §1).

## 4. The three questions of the block

**(a) Congruences for `a_{4k}(N)` modulo powers of an inert `l`, uniformly in `k`.** Trivial congruences hold: for
`l` not dividing `N`, `z^{(l^2-1) l^{a-1}} = 1 mod l^a`, so `k -> a_{4k}(N) mod l^a` depends only on `4k` modulo
`(l^2-1) l^{a-1}`. Also `z^l = zbar mod l` (Frobenius is complex conjugation), so on the circle `z^{l+1} = N mod l`
[Inferred, elementary]. Deeper Kummer-type congruences (supersingular Hurwitz numbers, Katz/Chellali-type) concern
`sum_{z != 0} z^{-4k}`, i.e. `L(lambda^{4k}, s)` at a critical point, summed over all `N` [Heuristic; literature not
checked; nothing in S2 states them]. "Multiplicity at an inert prime" for a cluster could only mean either collisions
modulo `l^a`, which are the valuations of `z_i - z_j` (Proposition A), or vanishing of an auxiliary function to
`l`-adic order `a`, which in dimension one is again the Vandermonde (Proposition B). Nothing lies beyond the
valuation of the Vandermonde. The supersingular structure gives no gain beyond `sum_{l<=M} log l/(l+1)` except the
convergent `0.2287 * M^2/2`, i.e. `O(M)` in `W`.

**(b) Chowla-Selberg, Siegel mass, theta correspondence "for norms".** For `Z[i]` the genus of `x^2 + y^2` has one
class, and Siegel-Weil is exact: `r(N) = 4 sum_{d|N} chi_{-4}(d)`. Arc counts are sums over `k` of coefficients of
the CM cusp forms `theta_{lambda^{4k}}`, which are orthogonal to the Eisenstein (mass-formula) part. So the mass
formula is blind to them. Chowla-Selberg for `Q(i)` gives the lemniscate constant
`varpi = Gamma(1/4)^2/(2 sqrt(2 pi))`, a constant normalisation. On average over `N <= X`, the number of pairs `z != w`
with `|z| = |w|`, `|z - w| <= C|z|^{1/2}` is `<< C X^{3/4} log X` [Heuristic count: for `d = z - w`, `w` lies on a
line with about `sqrt X/|d_0|` lattice points]. So almost all `N` have no close pair. This is an "almost all"
statement, capped by grafts (expert note §5 item 2). It is weaker in kind than the record's random-parent theorem
(`round31/README_IXL.md` §1.1), which covers every arc of almost every parent in a pool. Nothing new.

**(c) Whole circle versus one arc.** `a_{4k}(N)` sums over all `r(N) = 4 prod (e_p + 1)` points, and a congruence or
size statement for it constrains the cluster only after subtracting the other `r(N) - M` points, which are free.
For squarefree `N'`, `a_{4k}(N') = 4 prod_p 2 p^{2k} cos(4k theta_p)` is a function of the angles alone. A fake (real
angles) has the same product structure, and only the integrality of `pi_p^{4k}` separates it from a genuine circle,
but that integrality says nothing about which `M` of the points are close. The cluster statement
`sum_{i in cluster} u_i^{2k} ~ M u_1^{2k}` for `|k| << 1/delta` concerns a sub-sum, which is not the coefficient of
any automorphic object. Verdict: none of it constrains a single short arc.

## 5. Borrowable ideas (each with a first step and a probability)

1. **Inert all-level collision count (= S2's `gam(l)` exponent).** Already in the record; `O(M)` in `W`. First step:
   none needed (it is `two_thirds/proof.md` Lemma 3.3(I)). Probability for the rate: about 0%.
2. **Multiplicity determinants on `T` (Proposition B).** A negative result. First step: write out Proposition B with
   `D = u d/du` and record it as a barrier entry ("one-dimensional multiplicity = cut identity"). Probability for
   the rate: 0% (its purpose is to close the lane).
3. **Theta-jet determinants on the CM curve via kernels of the `z_i` (Programme C).** First step (one session):
   compute, for the explicit parent `q_P` with `P ~ 200`, the matrix of `D`-jets of `th0` for `E = C/Z[i]` at the
   points `zbar_i/q_P` (`q_P`-torsion), for pairs `i, j`, and test whether any entry's inert-`l` valuation
   correlates with `|z_i - z_j|` beyond `v_l(z_i - z_j)`. Expected: no correlation. Probability: below 0.5%.

Bottom line: the block's "D^p = 0 at inert p" is a statement about a height-2 formal group. Its trace on our problem
is the record's inert collision count, already used, with a convergent and rate-irrelevant tail.
