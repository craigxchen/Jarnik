# Round 32, block "slopes": the Bost slope inequality and the geometry-of-numbers layer of S2

Labels: **Read** (from S2, path given), **Inferred**, **Heuristic**, **Speculative**. "Record" means
`/home/user/Jarnik/docs/claude_review`. S2 means
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/dl_lean/x/lean/src`. Nothing in S2 was
compiled. About 40 reads and greps were used.

## 0. Verdict

1. S2's slope layer is standard Arakelov geometry of numbers for hermitian lattices over a CM field: Harder–Narasimhan
   filtration, the slope inequality for one injective map, its extension form, Minkowski II, and product-formula
   integrality for Gram determinants. S2 never mentions Bost by name. I found no general *filtered* slope inequality
   (deg-hat E <= sum over the filtration of rank times (mu_max + h)). In the one application I traced to the end
   (`Formalism/ProductFormula/StatementOCCM.lean`), the inequality is used as one injective map `mu_Q <= mu_max(R) +
   h_inf + h_fin`. The filtration is collapsed into a sum of m per-lift archimedean costs and m per-lift p-adic gains.
2. That engine is the Chudnovsky–Bost "arithmetic factorial" scheme. Each step pays an archimedean cost of about
   `a log a - a - a log Omega` (Taylor coefficients of order a). It gains p-adic divisibility of about
   `sum_i (a/p^i - w)_+` at the primes in the non-CM-type classes. This is "D^p = 0 on a CM curve for p inert". Our
   Theorem 1.1 is the same scheme on P^1: archimedean Fekete diameter against pigeonhole (Bhargava p-orderings) at
   q <= M, with the Mertens mass `sum log q / q` playing the role of `sum a/p`.
3. On our side, every slope inequality for evaluation at the cluster points on the conic reduces exactly to the cut
   identity. This holds with or without the vanishing-order filtration: multiplicity K scales both sides by K^2.
   Section 3 below computes it. Record: `round20/arakelov.md` Lemma 3.4 and Prop B.2.
4. No slope-type statement in S2 uses the full successive-minima profile or semistability in a way that is not
   Minkowski/Hadamard plus pair data on our objects. The integral determinants that fakes cannot imitate are the
   evaluation determinants. On P^1 they are the Vandermonde (pair data, gain zero). In rank > 1 they are generalised
   Vandermondes `V * Phi` whose extra factor needs route 5's zero estimate.
5. S2 does not get past the zero-estimate problem by a zero estimate. It gets past it by an exact period formula
   (seesaw/Rallis) for a determinant in rank > 1. The transferable lesson is only this: a rate needs an auxiliary
   determinant at a non-degenerate point of rank > 1 with a closed formula. The record already states that need
   (ZE^all_H, expert note §5 item 3). **Nothing in this block transfers to a rate.** Honest probability that any idea
   here gives M = o(log R / loglog R), even on a restricted class: below 1%.

## 1. What S2 actually states (results first)

### 1.1 Basic slope calculus

* **Read** `Formalism/LatticeSlopes/Statement.lean`, `Formalism.F.SlopeInequality`. Let L be a discrete Z-submodule of
  a Euclidean space E, M a discrete Z-submodule of F, L != 0, and phi: L -> M an injective Z-linear map with
  `||phi x|| <= C ||x||`. Then `slope L <= muMax M + log C`. Definitions: `deg L = -log covol L`,
  `slope = deg/rk`, `muMax` = sup of slopes of nonzero submodules. The module also has
  `MuMaxLeMin : (forall x != 0 in L, r <= ||x||) -> muMax L <= -log r + (1/2) log rk L` and `GramDetLe`
  (`det Gram(w) <= C^{2n} det Gram(v)` when `||sum c_i w_i|| <= C ||sum c_i v_i||`).
* **Read** `Formalism/HarderNarasimhan/Statement.lean`: `DegSupermod` (deg is supermodular), `MuChar`, `HNStep`,
  `HNSlopes` (the HN polygon: antitone slope list, `sum slopes = deg L`, `deg N <= hnPoly L (rk N)` for every
  submodule N), `Duality` (`muMax(L^*) = -muMin(L)`), and `SlopeIneq`. Its conclusion, verbatim:
  `(φ ≠ 0 → muMin ρ L - Real.log C ≤ muMin ρ' (range φ) ∧ muMin ρ L ≤ muMax ρ' M + Real.log C) ∧ (Function.Injective φ → L ≠ ⊥ → muMax ρ L ≤ muMax ρ' M + Real.log C)`.
* **Read** `Formalism/HarderNarasimhan/StatementPl.lean`, `SlopeIneqMean`. This is the O_K version: K is a number
  field (totally complex, with a star), and V and W carry hermitian metrics `X v` at each infinite place. The height of
  phi is `htMap = (1/n_inf) sum_v log ||phi||_v^2` and `muDef = 2 slope - log 2 + log|d_K|/[K:Q]`. Then
  `muMinDef A - htMap φ ≤ muMinDef (φ A)`, and for every B containing φ A, `muMinDef A ≤ muMaxDef B + htMap φ`,
  and, if φ is injective on A, `muMaxDef A ≤ muMaxDef B + htMap φ`. This is the textbook unfiltered Bost inequality.
  There is no finite-place height term in it. Finite places enter elsewhere (§1.3).

### 1.2 The "mean height descent" layer

* **Read** `Formalism/MeanHeightDescent/MuMaxMinkowski.lean`, statement `DescentSlope` (proved in
  `AverageDegreeBounds/DescentSlope.lean`). For φ injective with φ A ⊆ B:
  `rk A * (muDef A - htMap φ) + (rk B - rk A) * muMinDef B <= rk B * muDef B`. This is the extension form: the
  degree of B is at least the degree of the image plus the minimal slope times the corank.
* **Read**, same file, `GoodBadCut`. Let e be a self-adjoint idempotent (an orthogonal projector at every place) and N
  a lattice with `g * e(N) ⊆ N`, g >= 1. Then the ranks of `N ∩ im e` and `N ∩ ker e` add up to rk N, and
  `rk N mu(N) <= rk(N∩im e) mu(N∩im e) + rk(N∩ker e) mu(N∩ker e) + 2 rk(N∩ker e) log g`.
  Splitting a lattice along a projector of denominator g costs `2 log g` per rank.
* **Read**, same file, `SpanMuMin`. If B is K-spanned by vectors of norm `<= e^Λ` at every place, then
  `muMinDef B >= -2Λ`.
* **Read** `MeanHeightDescent/StatementA.lean`, `StatementF.lean`, `StatementJ.lean`. These are long real-number
  assemblies. `CovolumeEstimateOfSlopeChain` and `CovolumeEstimateOfAverageDegree` conclude
  `lcv(XW, M)/D <= sum_ρ pq * w_ρ * logOm + (64(...)+8Cc+2) * B0 * llx`, where `lcv = (2/[K:Q]) log ht`. So they
  give an *upper bound on the covolume* of a rank-D O_K-lattice M (a lower bound on its mean slope). The chain is:
  a twist lower bound `muU`; a period shift `cfl = sum_τ w_τ c1_τ - ...`; an extension step
  `dA (muA - P) + (dB - dA)(-Λ) <= dB muN` with `dA >= (1-ε) dB`; and a good/bad row split
  `dN muN <= dg muG + db (mmax + G)` with `db <= 3 δ dN`. "Descent" here means passing from a sublattice of almost
  full rank to the whole lattice for a cost of about `ε Q`. It is not a height descent in the Mordell–Weil sense.
* **Read** `Formalism/LatticeRigidity/StatementD.lean`.
  - `LatticeRigid`. Let A be a linear endomorphism preserving a full lattice M with `||A|| <= δ` and
    `δ n^2 < λ_1(M) λ_1(M^*)`. Then A = 0. (Integrality: a nonzero lattice vector has length at least λ_1.)
  - `Unipotence`. Let g be squarefree with g > 4n and unramified in E. Let γ be invertible, preserve an O_E-lattice L,
    and satisfy `(γ-1)^n L ⊆ P L` for every P | g. Suppose that at every complex embedding there is a positive H with
    `e^2 H - (γ^j)^* H γ^j >= 0` for 1 <= j <= n. Then `(γ-1)^n = 0`. So congruence plus archimedean
    near-unitarity forces exact unipotence.
  - `StatementA.lean` is linking-form and lattice-classification algebra. `StatementC.lean` holds a recursion lemma, a
    Chebyshev leading-coefficient lemma and a box sum. `PackCount` (in `PackCount.lean`) is a packing count with
    geometric-mean separation.

### 1.3 Where integrality ("a nonzero integer has absolute value >= 1") enters

* **Read** `ReciprocityTheory/GeometryOfNumbers/LatticeCovolumes/HadZSt.lean`, theorem `neg_log_det_le`. Let G' and G
  be positive definite of sizes r' <= r. Let M be an r × r' integer matrix, injective, with `A^2 G' - M^* G M >= 0`,
  and let N be an r × r integer matrix with `det N != 0` whose columns have G-length <= λ. Then
  `-(1/2) log det G' <= -(1/2) log det G + r' log A + (r - r') log λ`. The proof appends columns of N to M to form
  an integer matrix U with `det U != 0`, and uses `|det U| >= 1`. This is the slope inequality with the cokernel
  bounded by its last successive minimum.
* **Read** `Formalism/LatticeRigidity/HeightMinimum/SqrtAbsNormLeProd.lean`, `sqrt_absNorm_le_prod`. For
  `0 != l ∈ b` (a fractional ideal, E totally complex), `sqrt(N b) <= prod_w |l|_w`.
* **Read** `Formalism/ProductFormula/StatementOCCM.lean`: `ProductFormula`
  (`log|Δ| = sum_p v_p(Δ) log p` for Δ ∈ Z \ {0}) and `DetGain` (`v_p(Δ) >= g_p` on S implies
  `sum_{p∈S} g_p log p <= log|Δ|`).
* **Read** `Formalism/QuotientMetricDivisors/Statement.lean`. `GramDeterminantQuotientMetric`: the Gram determinant
  `gramR` of the quotient metric lies in the CM-real subfield and equals the archimedean Gram determinant. `LemmaRb`:
  `lcv(XQ, P M) <= sum_w lw(w) * v_w(gramR(quotGram P H)(e_w))` (finite sum over places, with local bases e_w).
  `LemmaRa`: if `l^g` kills `M^*/M`, then `v_w(gramR) <= g D v_w(l)`. Together these say that the covolume of the
  image lattice is controlled by the finite valuations of one algebraic Gram determinant. That is the product formula
  used as an identity, which is stronger than ">= 1".

### 1.4 The application: what the evaluation map is and how the inequality is used

* **Read** `Formalism/ProductFormula/StatementOCCM.lean`, structure `TypeParams` (fields quoted):
  `slope : muQ ≤ muMaxR + hInf + hFin`, `muQ_eq : muQ = (degN - degKr)/rho`, `rank : r * d ≤ rho`,
  `hInf_le : hInf ≤ (∑ j ∈ range m, sK j) + sP`,
  `sK_le : ∀ j < m, sK j ≤ a * Real.log a - a - a * logOm + C*k*f*b*Lcm a l + C*f*k^2*Lcm a l`,
  `hFin_eq : hFin = -∑ p ∈ S, 2 * Real.log p * vt p`,
  `comp : ∀ p, (∑ j ∈ range m, vK j p) + vP p ≤ vt p`, and
  `CardSubmodulesAndLinearEquivBounds : … gainCM l Ψ x a a1 w k m L0 p ≤ vK j p`, where
  `gainCM = ∑_{i ≤ x, p^i ∉ Ψ (mod l)} max(a1/p^i - w, 0) - a·Lcm/x^2 - (k+2m) log a / log p - 1 - L0/(p-1)`.
  Conclusion `OCCMEstimate`: `Gz x ≤ x ℓ + 2x ∑_{p∈Bad} log p/p + K·MainCM`. Corollary `ZetaQuarterPowerCorollary`:
  with `m ≍ x^{1/4}`, `k ≍ x^{1/2}`, `l ≍ x^{1/4}`, ζ(s) != 0 for Re s > 3/4.
* **Read** `Formalism/SizeWindowBrackets/GrowthBoundGivesRiemannHypothesis.lean`. The RH version takes the chain
  length `m ∈ [x^{1/2-ε}, x^{1/2+ε/2}]` (`SizeWindows`), which is the "m like sqrt X" of the remark.
* **Read** `Formalism/LocalToGlobal/Columns/CoefficientMapHeightBound.lean`, `coefficient_map_height_bound`. The map
  being evaluated is a coefficient map `cfMap λ g` from a tensor product of polynomial (Schur) representations of
  lattices M_τ to monomial coordinates with the Fock metric (`fockN f = sum |c_m|^2 m!`). Its height satisfies
  `htMap(cfMap λ g) <= sum_τ |λ_τ| * mean_v log ||g_τ||_v`, i.e. degree times the log-norm of the group element.
  `Formalism/FockGram/StatementFock.lean` supplies the functoriality: `FockGramSym` (Fock Gram of Sym^a is symGram),
  `FockPerm` and `DetSymPow`.
* **Inferred.** The evaluation map is a map from a lattice of automorphic objects (theta lifts with level and
  weight) to coefficient or jet coordinates. The archimedean height per lift is a Stirling term
  `a log a - a` (ph(u) = u log u - u in `Formalism/LogFactorialIndices/StatementSize.lean`). The finite height is
  negative and adds over the m lifts (`comp`). Each p contributes an arithmetic-factorial gain of about
  `sum_i a/p^i`. The rank of the source grows linearly in the chain (`rank : r*d <= rho`).
* **NOT FOUND.** I could not find the "cost" computation of timeline item 14 ("stop at 2/3 a priori"). The only hit
  for "two thirds" in the slope area is `Real.eulerMascheroniConstant_lt_two_thirds`
  (`SizeWindowBrackets/ChainSlopeWindows/PChain.lean`, line 1098), which is unrelated. The 3/4 of
  `ZetaQuarterPowerCorollary` is a different number.

### 1.5 The rest of the assigned list (statements only)

* **Read** `Formalism/SuccessiveMinima/StatementMink.lean`: Minkowski II, both directions
  (`MinkowskiIIUpper`/`Lower`/`Ball`). `StatementNP.lean`: ultrametric Newton slope, Strassmann
  (`#zeros in the r-disc <= wdeg a r`), Weierstrass preparation, zero formula. All standard.
* **Read** `Formalism/LatticeLocalGlobal/Statement.lean`. Definition names only: `UniqueForm`,
  `OneRationalStructure`, `IntegralOfLocalUnits`, `LatticeLocalGlobal`, `GaussSumUnit`, `GramChange`, `DualGram`.
  This is the local-global theory of lattices; I did not read the bodies.
* `HermitianLattices`, `HermitianCurvature` (1922 files), `MahlerMeasure`, `SelfDualUnimodular`, `ReducedBases`
  (LLL), `L2DbarEstimates` (Hörmander) and `DualLattices` were listed or name-scanned only. Nothing in their names or
  in the statements sampled suggests a slope inequality beyond §1.1–1.3.

**Answer to (a), first half.** S2 has no slope-type inequality that is not Minkowski/Hadamard plus the slope
inequality for one map. The HN polygon (`HNSlopes`: `deg N <= hnPoly L (rk N)`) is the semistability input. It is
used through `muMin`/`muMax` and `GoodBadCut`, never as a profile inequality that would bound counts of short vectors.

## 2. The record's side, as needed here

* `round20/arakelov.md` Lemma 3.4 (**record**). On P^1 over Z[i], the sublattice of sections vanishing to order m_j
  at z_j is `V·H^0(O(D - sum m_j))` with V monic. So the vanishing filtration is explicit and every slope inequality
  for it is pair data. Prop B.2: every genus-0 Arakelov inequality lies in the adelic-capacity class, which fakes
  satisfy. Theorem C(3) *Scope*: only exact global integrality of all pair numbers at once escapes, and that is "the
  problem itself".
* Expert note §5 item 3 (**record**). Single-place auxiliary forms are exactly critical. Free-nonvanishing grid
  determinants are subcritical (`round31/README_II.md` §3.1). Supercritical designs need ZE^all_H. Item 8: genus-0
  Arakelov reduces to pairs, τ^* and rank.
* `growth/auxiliary.md` Prop A (**record**, unrefereed). The Bombieri–Pila determinant on the circle is
  `V(z)·Phi(z)`, where Phi is a Z[i]-combination of Schur polynomials. At near-coincident points, Phi equals the
  confluent Wronskian `Wr(w)` up to `((1+η)^Λ - 1) Phi^+`. The p-adic content of det is `v_p(V) + v_p(Phi)`.
* Gram structure (**record**, `round31/README_VIIIL.md` §1.1, `README_IX.md` §1.6). `G_xy = 1 - 2D_xy/W` is PSD; the
  barycentre and moment step uses `E G >= 0` and the integrality `|t| >= 1`. The IX.6 Gram remark was scored
  "near-tautological".

## 3. Our evaluation map, with the slope inequality computed

### 3.1 Setting

Let z_1,...,z_M be distinct points of Z[i] with |z_j|^2 = N, on an arc with chords at most ℓ = C sqrt R, R = sqrt N,
normalised by the Gaussian gcd. Let `E_D = {f ∈ Z[i][z] : deg f < D}`. At ∞, give it the norm
`||f||_ρ = (sum_k |c_k|^2 ρ^{2k})^{1/2}`, where `f = sum c_k (z - z_1)^k`. Translation by z_1 ∈ Z[i] preserves the
lattice. For a multiplicity K >= 1 and D = KM, filter E_D by
`E^m = {f : ord_{z_j} f >= m for all j}`, m = 0..K. The graded piece `gr^m` maps injectively to
`G_m = Z[i]^M` by `f -> (f^{(m)}(z_j)/m!)_j`, with the standard metric at ∞ and the standard integral structure.

### 3.2 The inequality and its value

**Proposition S (Inferred; standard, consistent with record Lemma 3.4).** For D = KM the Bost inequality
`deg-hat(E_D) <= sum_m rk(gr^m) (mu_max(G_m) + h(φ_m))`, with exact finite heights, is equivalent to
`K^2 * [ sum_{i<j} log|z_i - z_j|_∞^2 + sum_{q finite} sum_{i<j} log|z_i - z_j|_q ] >= O(K^2 M^2) + (terms in ρ that cancel)`.
The bracket is 0 by the product formula. So the multiplicity filtration multiplies the D = M inequality by K^2 and
adds nothing.

*Reason.* On P^1 the evaluation-with-jets map E_{KM} → ⊕_m G_m is an isomorphism over K. Its determinant is the
confluent Vandermonde with Taylor normalisation, `prod_{i<j}(z_j - z_i)^{K^2}`. The filtered slope inequality is the
inequality "deg-hat(source) <= deg-hat(target) + log|det|_∞ + sum_q log|det|_q", and it is an equality. Bost's form
only weakens the finite terms, by replacing `v_q(det)` with lower bounds for it. Those lower bounds are the q-adic
M-point diameters of `round20` Theorem B, i.e. pigeonhole at q <= M and the chain at p | N.

**The same computation on the conic (Bombieri–Pila; Inferred, checked by hand).** Take the 2d+1 = M monomials
`z^a` (0 <= a <= d) and `z̄^b` (1 <= b <= d). Since `z̄_j = N/z_j`,
`det[z_j^a, z̄_j^b] = ± N^{d(d+1)/2} V(z) / prod_j z_j^d`, a nonzero Gaussian integer. Hence
`|det|^2 = |V|^2 N^{-d^2}`, and integrality gives `|V|^2 >= N^{(M-1)^2/4}`. This is exactly the parity part of the cut
identity. For M odd every `(2S_{pτ} - M)^2 >= 1`, so `prod_{i<j}|z_i - z_j|^8 >= N^{M(M-2)+1} = N^{(M-1)^2}`. The
archimedean bound `|V| <= (C sqrt R)^{M(M-1)/2}` leaves a margin of `N^{(M-1)/8}`, i.e. `(M-1)W/8` in logarithms.
That margin is all there is to spend. The small-prime collision mass, of order `M^2 log M`, spends it and gives
`W >> M log M`. Theorem 1.1's constant 3 is the precise version (`round20` Theorem B). Below the critical scale the
same determinants give Cilleruelo–Córdoba (literature, recalled).

**Which term would have to improve (Inferred).** There are three terms:

* The archimedean term, `prod |z_i - z_j| <= (δ_M ℓ)^{binom(M,2)}`, is sharp up to `e^{O(M^2)}`: the record clusters
  are nearly Fekete (`round20` §4.2).
* The finite terms at q <= M are local minima. Fakes attain them (Prop B.2).
* The terms at q > M are zero on both sides (`round20` §6.1).

On P^1 the determinant has no factor that is not a product of pair factors. So the term that must improve is a
**multi-point finite divisibility of an evaluation determinant**. That is impossible on P^1 and needs a rank > 1
target (X_M, the torus of monomials `prod z_i^{n_i}`). There the extra factor Phi must be shown nonzero with forced
divisibility, which is ZE^all_H.

### 3.3 The higher-rank lattices named in the task

* *Congruence lattices* `Λ = {z ∈ Z[i] : z ≡ z_1 mod H}` (H the common Gaussian part) are rank 2. Slopes of rank-2
  lattices are λ_1 and covolume, which is pair data (the cocycle, Part XI.2: "the problem itself").
* *Exponent-row lattice with the form G.* Write `x^T G x = sum_p log p sum_τ (sum_i x_i σ_{ipτ})^2`. For x ∈ Z^M each
  inner sum is an integer, so `q(x) >= (min log p) · #{(p,τ) : inner sum != 0}`. That is integrality on the rows. The
  rows are genuine in every fake, and the off-diagonal bounds `G_ij <= 4 log C' - 4 log|b_ij|` hold for fakes with the
  fake's chords. Any slope, Minkowski II, `GoodBadCut` or `DescentSlope` statement about `(Z^M, G)` therefore holds
  for fakes. **Capped.** Theorem 1.2 is the sharp use of this object.
* *Character lattice of monomials* `prod z_i^{n_i}`: this is route 5. Its auxiliary determinants are `V·Phi` (record
  Prop A). Free nonvanishing of Phi, e.g. by the near-diagonal Schur/Weyl-dimension expansion
  `s_λ(z) ≈ w^{|λ|} s_λ(1^M)`, is subcritical (`round31/README_II.md` §3.1). The hook formula is in S2
  (`Formalism/LogFactorialIndices/StatementSize.lean`, `HookFact`, `HookLog`; **Read**), but that adds nothing.

## 4. Answers to (b) and (c)

**(b) MeanHeightDescent and LatticeRigidity.** MeanHeightDescent bounds the covolume of a lattice from above. It
combines an extension step from an almost-full-rank sublattice (`DescentSlope`), a cut by an almost-integral projector
(`GoodBadCut`, cost `2 log g` per rank) and a min-slope bound from short spanning vectors (`SpanMuMin`). LatticeRigidity
says that small lattice endomorphisms vanish (`LatticeRigid`) and that congruence-unipotent, nearly unitary
endomorphisms are unipotent (`Unipotence`).

Exponent rows are genuine in fakes, so any descent or cut on them is fake-satisfied. The natural projector (onto
`sum x = 0`, g = M) gives `2 log M` per rank, the barycentre cost already used in Theorem 1.2. `Unipotence` would need
an endomorphism of a fixed lattice. The natural candidates `u_ij = z_i/z_j` are norm-one but not integral (they have
denominators at p | N'), and multiplication by them preserves no lattice independent of the cluster. **No transferable
descent step** (Inferred).

**(c) The fake test and the integral determinant.**

* Fakes have no lattice points. So `Norm det ∈ Z_{>=1}` for an evaluation determinant is a genuine input.
* The P^1/conic determinant (`N^{d(d+1)/2} V/prod z^d`, or the confluent version) is a nonzero Gaussian integer only
  for genuine clusters. But each of its consequences that a slope argument can use is a per-pair product formula
  with local lower bounds, and fakes satisfy those via pseudo-adelic completions (`round20` Theorem C(3)). **Gain:
  zero.**
* Triangle-area integrality (`det[1, z, z̄] ∈ 2i Z`) gives Jarník's scale R^{1/3} and is slack at √R by `R^{1/2}`;
  record item 6 lists it as capped.
* The only determinants whose integrality fakes cannot copy and which are not pair data are rank > 1:
  - the generalised Vandermonde factor Phi (route 5, needs ZE);
  - determinants of the integer cocycle `t_ij` with the exact additive relations (the problem itself, Part XI.2).
* Comparison with route 5: the slope language is route 5 restated. "Multiplicity" is the filtration. "Height of the
  map" is the Liouville/τ^* term. "Integrality" is forced divisibility. It contains no zero estimate.

## 5. Borrowable ideas, tested

| # | Idea from S2 | Closest analogue for arcs | Fake test | Char p / Liouville | In record? | First step | P(rate) |
|---|---|---|---|---|---|---|---|
| 1 | Arithmetic-factorial gain at inert primes, `gainCM ~ sum_i a/p^i` (`StatementOCCM.lean`) | Bhargava p-orderings of the conic at ℓ ≡ 3 (4): class counts `(ℓ+1)ℓ^{a-1}`, gain `sum_a E(M,(ℓ+1)ℓ^{a-1}) log ℓ` | Fakes carry point residues at q <= M (Prop B.2) | holds over F_2(t), F_3(t) (Thm 1.2 transfers) | yes: Theorem B, Theorem 1.2 | none needed | ~0 |
| 2 | Chaining m lifts: finite valuations add (`comp`), archimedean costs add, μ_max grows slowly; nonvanishing via period/seesaw identity | Products/tensor powers of the cluster on X_M^m or the torus, with a determinant having a closed formula at the cluster | Pure product versions are local, hence capped; a closed formula in rank > 1 would be global | τ^* < 2: products are exactly critical | ZE^all_H, expert note §5 item 3 | Search for one two-variable determinant at the cluster (monomials in (z_i, z_j) over ordered pairs) whose valuation at some q ∤ N, q > M exceeds the pair prediction on the 10-point cluster | <1% |
| 3 | `GoodBadCut` / `DescentSlope` on a lattice | Row lattice `(Z^M, G)` cut by the barycentre projector, g = M | Rows are genuine in fakes: capped | not tested; irrelevant | Theorem 1.2 barycentre | none | ~0 |
| 4 | Product formula on an algebraic Gram determinant (`LemmaRb`) | `det` of the Hermitian Gram of `(z_j^k)_{k<D}`, which by Cauchy–Binet is `sum_K |V s_{λ(K)}|^2` | Factorwise V·Phi: capped (grid determinants) | — | README_II §3.1, auxiliary Prop A | none | ~0 |
| 5 | Fock metric with exact Sym-functoriality (`FockGramSym`, `coefficient_map_height_bound`) | Natural metric on polynomials in (z, z̄) making `h(Sym^a φ) = a h(φ)` exact; removes `e^{O(M^2)}` slack | Changes only the archimedean constant; record clusters are already nearly Fekete | — | Remark B.1 (constants only) | none | 0 |
| 6 | Transcendental auxiliary function with CM divisibility along the weight direction: jets in k of `k -> sum_{cluster} e^{4ikθ_j}` (Hecke characters `λ^{4k}`) | **Speculative.** Divided differences of order a in k at the cluster are tiny, about `(4 a C R^{-1/2})^a`. Integrality would need the full theta coefficient `a_k(N) = sum_{|z|^2=N} z^{4k}`, which mixes in all `4·2^s` points of the circle | The other points enter as an averaged exponential sum: the graft caps apply (expert note §5 item 2) | unknown | large-sieve and character lanes | Write `Δ_k^a a_k(N)` for q_P and check whether the cluster part can be isolated without an averaged estimate | <0.5% |

**Bottom line.** S2's slope machinery is the same arithmetic-factorial / Fekete scheme that already gives Theorem 1.1.
Its new ingredient, nonvanishing of a rank > 1 determinant through an automorphic period identity, has as its only
analogue for us the missing zero estimate ZE^all_H. That estimate is not easier in slope language. The block yields
no transferable idea that passes the fake filter. What is new for the record is a clean dictionary:

* "Bost slope inequality with multiplicity" on the conic = K^2 × the cut identity.
* The Bombieri–Pila determinant on the conic = `N^{d(d+1)/2} V / prod z^d`, whose integrality is the odd-M parity part
  of the cut identity.
* The S2 "D^p = 0" gain = the inert-prime p-ordering gain already in Theorem B.
