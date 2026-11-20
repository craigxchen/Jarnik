# Block "lifts": iterated theta lifts, cost, "division by a gcd" — what the Lean does, and what transfers to short arcs

Sources: S1 = remark.pdf (timeline items 7-10, 12, 14-16, 18). S2 = Lean tree at
`.../scratchpad/dl_lean/x/lean/src` (abbreviated `src/`); NOT compiled by anyone here, correctness unverified.
Record: `/home/user/Jarnik/docs/claude_review/expert_note/expert_note.md` §§1, 5; `round31/README_XI.md` §5.
Labels: **[Read]** = statement read in S2 (path given); **[Inferred]**; **[Heuristic]**; **[Speculative]**.

## 0. Results first (verdict)

1. In S2 "multiplicity m" is **not** an order of vanishing and not a single p-adic divisibility. It is the
   **number of theta lifts in a chain** U(N(0)) -> U(N(1)) -> ... -> U(N(m)), N(j) = k + j(t+1), between definite
   unitary groups over a CM field. Each lift j contributes its own copy of an inert-prime divisibility ("gain j w")
   to one global Arakelov degree. The m copies add linearly; the error terms that do not scale with m are divided by
   m on normalisation, and m >= x^{1/2-eps} is imposed. **[Read]** for the data, **[Inferred]** for the bookkeeping.
2. The quantity bounded is not "a nonzero integer, divisible and small". It is the **slope (normalised Arakelov
   degree) of a hermitian lattice** that receives near-orthogonal integral maps from a family of lattices indexed by
   "words" (twists by characters ramified at an auxiliary split prime T). An averaged Bost-type slope inequality
   bounds a weighted von Mangoldt sum in residue classes by an explicit main term plus error. **[Read]/[Inferred]**
3. The theta-lift machinery itself in S2 is elementary linear algebra on finite class sets: lifts are finite
   kernel operators with mass weights 1/|Aut|; the "seesaw" is finite Fubini; "Rallis" is "T T* = scalar on a
   subspace". The arithmetic content sits in hypotheses/fields (gains, heights, slopes) packaged in structures.
   **[Read]**
4. "Stop at 2/3": **not found** as a statement in S2. "Division": read as the factors `divDe` (division by powers of
   archimedean normalisers) in each lift kernel and the Gauss operators on the ramified part of the twist family.
   **[Read]** for the objects, **[Inferred]** for the identification.
5. Transfer to arcs: **essentially nothing transferable at the level of the rate.** The S2 mechanism improves an
   error term of an averaged statement (a sum over all primes up to x) whose main term is exact; our problem needs a
   change in the main-term order (W/(M log M) -> infinity) for an individual configuration. The three closest
   analogues (rank-1 Siegel-Weil; the tower z -> z^k / chain amplification; gluing across Hecke characters or sign
   flips) are respectively trivial, exactly critical (homogeneous), and capped by the fake/residue barrier.
   Probabilities of reaching M = o(log R/loglog R) by any of them: at most 0.5% each. Details below.

## 1. What S2 actually contains (statements in ordinary language)

### 1.1 Finite theta lifts and their formal properties — `src/Formalism/ThetaLiftIndependence/Statement.lean` [Read]

* `kerOp w K f (j) = sum_i w_i f_i K(i,j)`; `chainOp w K m` = m-fold composite of such kernel operators between
  finite index sets iota(0), ..., iota(m).
* `ChainInt` (Prop): if `D j * w j i` and `K j a b` lie in a subring S for all j, and f is S-valued, then
  `(prod_{j<m} D j) * chainOp w K m f y` lies in S. In words: **denominators multiply along the chain**; with
  weights w = 1/|Aut L| the denominator per step is a common multiple of automorphism-group orders. `ChainIntM` is
  the module version (lattices Lambda j stable under S, kernels mapping Lambda j into Lambda (j+1)).
* `ClassRepsDefinite`: for a totally definite hermitian form over a CM field and a lattice L, a finite set of class
  representatives of the genus exists, every weight `latWt h M = 1/|Aut M|` is the inverse of a positive integer, and
  there is D > 0 with D * latWt integral for all M in the genus.
* `thetaLift h s f tau = sum_{L in s} |Aut L|^{-1} f(L) theta_L(tau)`; `LatLiftIndep`/`ThetaLiftIndep`: independence
  of the choice of class representatives.
* `LatSeeSaw`: `<sum_L wt(L) f(L) k(L,.), g>' = <f, sum_{L'} wt(L') g(L') conj k(.,L')>` — finite Fubini, i.e. the
  seesaw identity in its trivial finite form.
* `RallisAdj`: if `<theta x, theta y> = c L (prod_v Z_v) <x, y>` on a subspace S stable under theta* theta, then
  `theta* theta = c L prod Z_v` on S. `ChainScalar`: if `T_j T_j* = c_j` on S_j and T_j* S_j in S_{j+1}, then the
  composite satisfies `T T* = (prod_j c_j)` on S_0. **The archimedean norm of the iterated lift is known exactly**
  (as a product of Rallis constants: L-value times local zeta integrals), not just bounded.
* `DetChain`: if A acts on eigenspaces E_s by `kappa * prod_j lam_j(s)` and B_j acts on matching spaces by
  `lam_j(s)`, then `det A = kappa^{dim E} prod_j det B_j`. `DetSwap`, `EigenSwap`: det(T*T) = det(TT*) on the
  orthocomplements of kernels; nonzero eigenspaces of BA and AB have equal dimensions.
* `GaussOp`: for psi a nontrivial additive character of a finite field F and rho an action of F^x, on vectors with
  `sum_c rho(c) v = 0`, `G G' v = |F| v` where `G = sum_c psi(c) rho(c)`, `G' = sum_c psi(-c) rho(c)^{-1}`.
  `FFPlancherel`, `FFRamCancel`: finite Fourier identities.

### 1.2 The chain actually used — `src/Formalism/LocalRingCardinality/` [Read]

* `Params.lean`, `structure Ranks`: integers k, m, t with k, m even, t odd, m + 3 <= k, 3m < t, k + m + 4 <= t,
  `m^2 <= k + m(t+1)`; ranks `N j = k + j(t+1)`. So lift j goes from U(N(j-1)) to U(N(j)), adding t+1
  variables: this is the timeline's `pi -> pi boxplus mu_j[t]` (an Arthur parameter gaining a block of size t, plus
  one) **[Inferred]**. `HCParam`: Harish-Chandra parameters u (half-integers), v, b, e (integers) with interlacing
  u > v > u and gaps `e_j - b_j >= (t+3)/2`, `b_{j+1} - e_j >= (t+3)/2` (the "interlacing weights" of item 7).
* `DataKernels.lean`, `structure ColumnChain R iota hc`: finite class sets `cls i` (genus classes on U(N(i))),
  stabilisers with `|stab i c|` dividing one integer `DGam`; vector-valued forms `cls i -> Y i -> C`; kernel
  `Kmap i L L' = rhoDet(L) rhoDet'(L') . divDe_i o unf(L,L') o sharp_i`, where `sharp_i` multiplies by
  `prod_tau de(i,y,tau)^{sSharp}` and `divDe_i` divides by `prod_tau de(i+1,y',tau)^{sDiv}`, with exponents
  `sSharp = e_i - N(i+1)/2`, `sDiv = e_i - N(i)/2`; `rhoDet` unimodular. `col = chainOpM wt Kmap m` is the m-fold lift.
* `structure ColumnsAtT`: an auxiliary prime T, f places v with residue fields of size T, an action `rhoT` of
  F_T^x on forms (twisting by characters of conductor T), `Mram = Mu0 ∩ ker(sum_c rhoT c)` (the purely ramified
  part), `gauss v = gaussOp`, numbers `Dstar` (dimension, a Siegel mass times archimedean dimensions),
  `Xu = (Dstar/duStar)^{1/(fk)}`, `Tstar = Xu * GamP^{-2/(fk)}`.
* `DataParams₂.lean`, `structure SizeData₂ D`: a scale `a = x b`, `b = exp(x^eps)`; for each lift j in
  `lifts D = {1..m}` and place w a real `gain j w`, zero off a support, zero at non-inert places
  (`not_inert`), and (in `SharpRestriction/GlobalTheoremInputs.lean`, `LinkG_at`) bounded by
  `0 <= gain j w <= 4 |TT| a / p` for inert w above p.
* `DataArith.lean`, `structure DivData`: a target hermitian lattice `RR`, lattices `Nu u` for each word u, maps
  `Bk u : Nu u -> dual(RR)` with height `<= h u`; `sFrame` = rank-weighted mean of `mu(Nu u) - h u`.

### 1.3 The final inequality — `src/Formalism/RankGapBounds/MeanBoundSecondForm.lean`, `src/Formalism/SizeWindowBrackets/SChain.lean` [Read]

* `Data.DIV`: `muCC + h - Fo <= sum_t ( m (ph(a_t) - a_t logOm) - a_t log T ) + errors`, with
  `ph(y) = y log y - y`. The arithmetic degree of the target lattice is bounded by **m times** `a log a - a`
  (= log a! up to lower order), minus `a log T`. **[Read]**
* `Data.SlopeIneq` (averaged slope inequality): `s - log|W| - (1-del)^{-1}(4 uF log(1 + HH/uF) + (eta + 4 sqrt e) HH)
  <= muCC`, i.e. the target's slope is at least the mean over words of `slope(Nu u) - h(u)`, minus `log #words`,
  minus a frame-tightness loss. **[Read]** Its frame input is `CorankShare/StatementFR.lean`, `FrameLA`: if the
  operators `c_i B_i* B_i` are close to identity in Hilbert-Schmidt norm on average (<= e dim), the frame operator
  satisfies `tr(Q^2) <= (1+eta) tr(Q)^2/dim H`, and subspaces of codimension <= eps dim are chosen, then the
  orthocomplement of the joint span has dimension `<= (2 eta + 4 eps + 4 sqrt e) dim H`. `LinkS_at`
  (`SharpRestriction/GlobalTheoremInputs.lean`): off-diagonal HS norms `<= (1+eta0) d_u d_u'/Dstar` — "approximate
  orthogonality of the columns" (item 16).
* `MeanBoundSecondForm`: under DIV, slopes, losses, SlopeIneq etc., `mean_u0 ASK <= mean_u0 Main + ErrX` with
  `Main = sum_{j<m} sum_t a_t (log a_t - 1 - gamma + L'/L(1) - mean_rho log w)`. ASK is a sum over j and t of a
  weighted von Mangoldt class sum. So: **von Mangoldt sum in classes <= explicit main term + error**, then Landau
  (sign-only) gives zero-free half-planes (`FixedEpsilonGrowthBound`: `Gz x <= x^{1/2 + (cTop7+3) eps}`).
* `SChain.lean`, `RestSize`: `x^{1/2 - eps} <= m` and `k <= n^2 x^{1/2-eps}`. This is "m like sqrt X". **[Read]**
  Why sqrt: errors such as `cDiv nT a La` (not proportional to m) become `~ x b La / m` per unit of multiplicity,
  while per-lift costs carry factors `m^2 L / k`, `m B0`; balancing x/m against m gives m ~ sqrt x. **[Inferred]**
* Siegel mass: `HermitianLattices/StatementMass.lean`, `SiegelMassRatio`, `MassTelescope`: mass(U(r)) /
  mass(U(k)) = prod_{M=k+1}^{r} massRatioRightSide K M (Euler factors and L-values). `TranslateWordBoundHypotheses.lean`
  l.125-133: `massr = massk * prod_M ((M-1)! (vE/pi)^M)^f / sqrt((2vE)^f) * L(M)` and `Dstar = massr * Isig * prod dsig`.
  This is item 9 ("norms via Chowla-Selberg, Weyl, Siegel mass"): dimensions and norms of the lifted spaces. **[Read]**
* `EisensteinCokernelLengths/StatementLIFTV.lean`, `LagrangeIdealSg`: if `c_a - c_b` is associated to
  `l^{sg + v_l(a-b)}` for a in S1, b in S1 ∪ S2, then `l^{max_a(sum_b (sg+v_l(a-b)) + sum_{a'} (sg + v_l(a-a')))}` lies
  in the ideal `(prod_{S1}(X - c_a), prod_{S2}(X - c_b))`. `WindowDiff`: for S in a window of length y,
  `sum_{b in S, b != a} v_l(a-b) <= (y+1)/(l-1) + log_l D`. `ExponentChainBound.lean`, `GapG`, `RRVal`: l-adic
  Legendre-type valuation bounds for products `prod_j c_j` where `v(c_j - c_j') = v_l(j - j')`. These are the
  "divisibility of the m-fold lift" bookkeeping (items 2, 10). **[Read]**
* `NormUnitaryBounds/Statement.lean`: Bombieri/Fischer norm on polynomials in N x n matrix variables is U(n)-
  invariant (`NormUnitary`); highest-weight subspaces `H R nu` (invariants under upper shears); a weighted count
  identity (`NormUnitaryCount`). This is the archimedean Fock model of the dual pair. **[Read]**

### 1.4 What I could not find

* No statement producing "2/3". Only `(2/3)^f` in a dimension lower bound (`TranslateWordBoundHypotheses.lean` l.129)
  and `2/3` weights in `SharpRestriction/GrowthBounds/SegmentDegree/Coarse.lean` (`thc`). The a-priori 2/3 of the
  timeline is presumably an abandoned intermediate. **NOT_FOUND**.
* No Siegel-Weil formula with explicit local densities for degenerate (rank-deficient) targets; S2 has abstract theta
  series `thetaN` at matrix arguments (`HermitianLattices/StatementTheta.lean`), the mass formula, and a theta-
  injectivity lemma (`Development/WeilBound/ThetaLifts/LatticeBasics/SiegelWeilAnisoOf.lean`, `theta_inj`). **NOT_FOUND**.
* The final theorem goes through `every_character_of_literature_of_vWin (O : OpenLiteratureStatements)`
  (`Development/AutomorphicArithmetic/ThetaCorrespondence/ThetaPairLifts/ThirdPlaceSeries/EveryCharacterOfLiteratureOfVWin.lean`
  l.810), instantiated in `CharacterBoxes/MainConclusion.lean`. Whether that structure's fields are honestly proved
  was not checked (outside this block).

## 2. Translation to lattice points on short arcs

Notation: cluster z_1..z_M on |z|^2 = N in an arc of length C sqrt R, normalised by the Gaussian gcd; W = log N';
d_ij = z_j - z_i = 2 i eta H_ij t_ij. Target: W/(M log M) -> infinity.

### 2.1 Question (a): representation numbers of near-degenerate hermitian matrices

Statement. Let L = Z[i] (rank-1 definite hermitian lattice, one class in its genus). For an M x M hermitian matrix
T, r(L,T) = #{x in L^M : x_i conj(x_j) = T_ij}. For T = (z_i conj z_j) from a cluster, r(L,T) = 4 (the tuple is
determined up to a unit). The Siegel-Weil average over the genus equals r(L,T) because the genus has one class.
**Nothing survives in rank 1.** [Inferred, elementary]

Real version: G = (Re z_i conj z_j) is an integral PSD matrix of rank 2 with diagonal N and
`N - G_ij = |z_i - z_j|^2/2 in [0, C^2 R/2]`. All 3 x 3 minors vanish. These are the Cayley-Menger/Ptolemy
relations among squared chords; abc and triangle consequences are capped (record item 6). Counting T in the box is
counting clusters. S2 has nothing on near-degenerate local densities (§1.4). **Verdict: no input.**

### 2.2 Question (b): "multiplicity over the integers"

Cleanest S2-faithful formulation (programme P_b). Find a chain of integral linear maps
T_1, ..., T_m between lattices attached to the cluster such that
 (i) `T_j T_j* = c_j` exactly on a subspace (Rallis analogue), so `|det(T_m ... T_1)|` is known exactly;
 (ii) the composite has denominators `prod D_j` (ChainInt analogue);
 (iii) each T_j carries a forced divisibility g_j coming from the cluster;
and show `sum_j (g_j - log D_j - (1/2) dim log c_j) > 0` grows faster than the main term.

Smallest-case test (tower z -> z^k). Points z_i^k lie on |w|^2 = N^k. Chords: `|z_i^k - z_j^k| ~ k R^{k-1} |z_i - z_j|
<= k C R^{k-1/2}`, against the critical scale R^{k/2}: for k >= 2 the image is not a short arc. Algebraically,
`z_i^k - z_j^k = prod_{zeta^k=1} (z_i - zeta z_j)`; for k = 2, `(z_i - z_j)(z_i + z_j)`, and the cut identity for the
k = 2 tower is the cut identity for the 2M points {z_i} ∪ {-z_i}. The new factors |z_i + z_j| ~ 2R satisfy their
lower bound trivially. Gain per step equals cost per step: the identity is homogeneous. This is the record's
"Liouville is exactly critical" (expert_note §5 item 3). [Inferred; hand-checked at M = 2, k = 2.]

Why S2 can win and we cannot [Heuristic]: S2's m lifts carry **different** integral structures (different groups
U(N(j))) whose per-lift archimedean cost is fixed by the exact Rallis constant, and the quantity improved is an
**error term** relative to an exact main term (explicit formula). Our cut identity is also exact, but the
deficit is in the **main-term coefficient** (3 in 3 M log M), not in an error. Amplifying by m multiplies main term
and gain alike.

Speculative candidate (degree-g Siegel/hermitian theta at a CM point attached to a g-subset). The tuple
(z_i)_{i in S} defines `E_i -> E_i^g, x -> (z_i x)`; one could look for a hermitian modular form of degree g with
small coefficients vanishing to high order at the associated CM point, and use the inert-prime divisibility of
Taylor coefficients at supersingular CM points (timeline item 2). Problems: (1) the Gram matrix is rank one, so no
polarisation is built from it directly; (2) the Faltings height of a CM abelian variety depends only on the CM type
(Colmez), so the archimedean side cannot see the arc; (3) inert-prime divisibility holds for every N, so it does not
distinguish clustered from unclustered configurations. **[Speculative; I expect it to be empty.]**

### 2.3 Question (c): "2/3, then division by a gcd across a family"

The numerical coincidence: our 2/3 = 2/c with c = 3 in W >= 3 M log M (expert_note Theorem 1.1); the remark's
2/3 is a gain/cost ratio per lift. They are unrelated beyond the number. [Heuristic]

S2's "division" [Inferred from §1.2-1.3]: (i) each lift divides by powers of archimedean normalisers (`divDe`), and
(ii) on the ramified part of a twist family modulo an auxiliary split prime T, Gauss operators satisfy
`G G* = T` exactly (`GaussOp`), so a common factor sqrt(T) per place can be divided out; near-orthogonality of the
twisted columns (`LinkS`, `FrameLA`) lets the averaged slope inequality pay only `log #words`.

Our analogues and their fate:
* Family = twists of the theta series of Z[i] by characters chi mod a split prime tau. The coefficients are
  `a_chi(N) = sum_{|z|^2 = N} chi(z)` and `sum_chi |a_chi(N)|^2 = (T-1) #{(z,z'): z ≡ z' mod tau (up to units)}`.
  This is a collision count modulo tau: a character-tested point residue. **Capped at every modulus at once**
  (expert_note §1, P3-char, round 31 Part IX.3). Fails the fake test.
* Family = the 2^s sign flips pi_p <-> conj(pi_p) on one norm q_P. A flip is not induced by a field automorphism; it
  preserves the exponent profile and changes the angles. An average over flips sees only the profile, i.e. exactly
  what angle models see. **Fails the fake test** by construction. [Inferred]
* Family = Hecke characters lambda^{4k}, k <= K ~ q_P^{1/4}. Gluing across k is a mean square
  `sum_{k<=K} |sum_{|z|^2=q_P} (z/|z|)^{4k}|^2`; its off-diagonal part is the close-pair count. Non-sharp large
  sieves are capped (record item 2); the sharp-constant forms are listed as uncapped but S2 offers nothing toward
  them. [Inferred]
* Family = circles N |gamma|^2 (multiply by gamma): `|gamma| C sqrt R` against `sqrt(R |gamma|)`: not critical for
  |gamma| > 1. No family. [Inferred]

Averaged Bost slope inequality on P^1 at the cluster points is already in the record and reduces to pair inputs
(`round20/arakelov.md` §3.4, overview §5.4 row 18).

Also: even granting S2's claim (GRH for Dirichlet L-functions), nothing follows for q_P: the residual problem is a
signed subset sum of s fixed angles at precision e^{-s log P/4}, a Baker / Lang-Waldschmidt question, while GRH
(even for Hecke L-functions of Q(i)) controls sums over primes p <= x, not exponentially fine additive relations
among a fixed set of angles. [Inferred]

## 3. Barrier tests (summary table)

| idea | fake test | char p | single-place Liouville? | in record? |
|---|---|---|---|---|
| rank-1 Siegel-Weil / Gram representation numbers | vacuous (trivial count) | trivial there too | n/a | Ptolemy/abc capped (item 6) |
| chain amplification, tower z -> z^k | homogeneous: same as cut identity | cut identity holds over F_q(t) | yes, exactly critical | item 3 |
| l-adic Legendre window bounds at inert l = 3 mod 4 (`LagrangeIdealSg`, `WindowDiff`) | point residues: capped | holds verbatim | local | Thm 1.2 uses inert l; collisions at q <= M |
| twist family mod tau + Gauss division | P3-char: capped at every modulus | holds | no | IX.3 |
| sign-flip family average | sees only profile: fakes pass | holds | no | fakes construction |
| averaged Bost slope + frame lemma | reduces to pairs | holds | no | round20/arakelov.md §3.4 |
| CM-point Taylor divisibility at degree g | unclear; probably blind to the arc | supersingular analogue exists | unclear | not in record |

## 4. Borrowable items, first steps, probabilities

1. **Frame-to-span lemma** (`CorankShare/StatementFR.lean`, `FrameLA`, `LayerCake`, `FrameLast`). A clean
   quantitative statement: near-tight frames of maps have joint span of small corank. Use: any future averaged slope
   argument over a family of evaluation maps. First step: restate `FrameLA` with constants and check it on the
   family of evaluation maps of degree-D polynomials at the M cluster points twisted by the 4 units. Probability of
   rate: < 0.5% (the family it would average is fake-visible).
2. **Exact-norm chain (ChainScalar/RallisAdj/DetChain)**. Principle: compute the archimedean side exactly through
   T T* = scalar instead of bounding it. First step: look for an integral operator built from the exact cofactors
   t_ij (the input fakes lack) with T T* scalar. The obvious candidates (Fourier on angles, Gauss sums on
   residues) satisfy Parseval for fakes as well. Probability: < 0.5%.
3. **l-adic Legendre bookkeeping for inert l with prime powers** (`WindowDiff`, `RRVal`). On the norm-one torus
   over Z_l, l = 3 mod 4: `T^1(Z_l) = mu_{l+1} x (1 + l O)^{N=1} ≅ mu_{l+1} x Z_l`, so for w = z_i/z_1 in one
   mu_{l+1}-class, `v_l(w_i - w_j) = 1 + v_l(a_i - a_j)` with a_i in Z_l. First step: prove
   `sum_{i<j} v_l(z_i - z_j) >= sum_{k>=1} (M^2/(2 l^{k-1}(l+1)) - M/2)` and insert it into the proof of Theorem 1.2
   (expert_note §1) to sharpen K_C. Gains a constant at most; residue-type, capped. Probability of rate: ~0.

**Bottom line.** The lift block of S2, read at the level of statements, is an averaged Arakelov slope argument
whose "multiplicity" is the number m of lifts, each adding a copy of an inert-prime factorial divisibility to a
global degree, with m ~ sqrt x chosen to divide away errors not proportional to m. That is a device for error terms
of averaged prime sums. Our obstruction is a main-term coefficient for an individual configuration, and every
family over which S2-style gluing could run (twists mod tau, sign flips, Hecke characters lambda^{4k}) is either
capped by the record's fake/residue barrier or is the record's own large-sieve problem. The block is too far from the
arc problem to transfer; I recommend no work on it beyond, at most, item 3 as a constant-sharpening exercise.
