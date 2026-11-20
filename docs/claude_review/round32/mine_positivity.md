# Round 32, block "positivity": sign-only arguments, Weil positivity, Behrend sets and subset sums in the downloaded Lean development, and what transfers to lattice points on short arcs

Labels: **[Read]** means taken from a statement in S2, with its path. **[Inferred]** means a conclusion I drew from statements and imports, not from proofs. **[Heuristic]** and **[Speculative]** mean what they say. **[Proved here]** marks short arguments written out in full below. Paths in S2 are relative to `.../scratchpad/dl_lean/x/lean/src/`. Nothing in S2 was compiled or executed. The Lean is untrusted, and every statement quoted from it is a *claimed* lemma.

## 0. Results in six lines

1. In the Lean, the "sign" is **Landau's theorem**, not a Weil-positivity functional and not a path-quadric count. For a real character χ the coefficients Λ(n)(1+χ(n)) are nonnegative. So a **one-sided** bound `G(χ,x) ≥ −C x^{σ₀}`, together with no real zeros in (σ₀,1), gives no zeros of ζ or L(·,χ) with Re s > σ₀ (§1.1). The "qN = M^k + R" path-quadric proof, Weil positivity `W(1_S) ≥ 0` and its determinant do not appear under any name I could find. I found no lemma that matches them (§1.2).
2. `FockGram` is a determinant identity for Fock (Bargmann–Fischer) Gram matrices of symmetric powers: `det Gram(Sym^a e) = det(Gram e)^{C(n+a−1,n)} · ∏ t!`, together with positivity preservation. It is lattice and covolume bookkeeping, not a positivity argument about W (§1.3).
3. `BehrendSets` constructs Behrend-type sets T ⊂ Z/n with no nontrivial relation αx₁+βx₂+γx₃=0 for small coefficients. Each translate of T meets every arithmetic progression of length 2L in at most two points, and meets every short orbit in at most one point. The group is the 2-power cyclotomic Galois group (Z/2^{s+2})^×/±1, so T selects levels of the Iwasawa tower that Frobenius orbits of small primes rarely revisit. A random-alteration lemma comes with it (§1.4).
4. The subset-sum modules count something else entirely. `JacquetSubsetSums` says a Jacquet-module parabolic type must be a subset sum of block sizes. `PadicSubsetSums` counts subsets of irreducible factors with a given total degree, i.e. degrees of divisors over finite fields. Neither concerns signed sums of real angles (§1.5).
5. Our side, (a): the sign-definite moment `Σ_k w(k)|a_{4k}(N)|^{2j}` is an exact restatement of a near-coincidence count. Its exact diagonal is the additive energy of the cube {±1}^s, namely 4^{2j}·C(2j,j)^s. So the best it can give is `M ≲ C(2j,j)^{s/(2j)} ≥ 2^{s/2}`, which is exponential in s, for every j. The off-diagonal part is exactly the residual problem, and angle fakes satisfy the identity term by term. It is capped (§2). (b): every principal minor and Schur complement of the cluster Gram matrix G holds automatically for any exponent profile, fakes included. Combined with the pair bounds, positivity gives nothing beyond the Perron root of the collision matrix, and the record already has that (§3). (c): Behrend-type constructions can at best change constants in fakes or in K_C, never the rate (§4).
6. Nothing from this block transfers. Every candidate I tested is either fake-transparent or already in the record (§5). The probability that anything here leads to M = o(log R/loglog R) is at most 1%.

## 1. What S2 actually contains for this block

### 1.1 The sign-only step: Landau's theorem for nonnegative coefficients [Read]

`Formalism/StatementOC.lean` defines
* `Sz x = Σ_{n≤x} Λ(n)(x/n − 1)`;
* `Gz x = Sz x − x log x + (1+γ)x`;
* `Gchi χ x = Σ_{n≤x} Λ(n) Re χ(n) (x/n−1) + Re(L'/L(1,χ))·x`.

This is timeline item 5's functional.

`Formalism/Statement.lean` defines `LamF χ n = Λ(n)(1+Re χ(n))` (≥ 0 for real χ), `S χ x = Σ_{n≤x} LamF χ n (x/n−1)` and `G χ x = S χ x − x log x + (1+γ_F(χ))x`. Its key statement is:

```
def LandauStep : Prop := ∀ (N) [NeZero N] (χ : DirichletCharacter ℂ N), χ ≠ 1 → (∀ n, (χ n).im = 0) →
  ∀ σ₀, 0 < σ₀ → σ₀ < 1 → (∀ σ : ℝ, σ₀ < σ → σ < 1 → LFunction χ σ ≠ 0) →
  (∃ C, ∀ x ≥ 1, -C * x ^ σ₀ ≤ G χ x) → ∀ s : ℂ, σ₀ < s.re → riemannZeta s ≠ 0 ∧ LFunction χ s ≠ 0
```

`RHofG χ` is the same statement with σ₀ = 1/2 + ε. `Formalism/StatementU.lean` has `ZetaLandau`: for a sign s₀ = ±1, `s₀·Gz x ≤ C x^θ` for large x implies ζ(s) ≠ 0 for Re s > θ. It also has smoothed variants (`Gzbar`, an average of Gz(x/u) over u ∈ [1,2]) and sampling lemmas (`SampleAverageApprox`, `RiemannSumApprox`). `Formalism/ZetaNonvanishingChebyshevWeights.lean` averages Gz against Chebyshev nodes.

`Formalism/SizeWindowBrackets/PowerBoundGivesRiemannHypothesis.lean` has
`PowerBoundGivesRiemannHypothesis : ∀ C ε₁ > 0, (∀ ε ∈ (0,ε₁), eventually Gz x ≤ x^{1/2 + Cε}) → RiemannHypothesis`.
`GrowthBoundGivesRiemannHypothesis.lean` (`RiemannHypothesisOfReducedGrowthData`) reduces this one-sided bound to the existence, for every large x, of "reduced growth data". These consist of a finite abelian group H, characters ψ of E0 × H, a set TT, a thick set T2 ⊂ H, weights w and a size-window inequality `SizeWindows`. The data include `topExcess ≤ Err`, a mean over translates u₀ ∈ H of a sum over TT of an excess `ASK`.

**What is sign-definite and why [Inferred].** The Dirichlet series with coefficients Λ(n)(1+χ(n)) ≥ 0 is −ζ'/ζ − L'/L(·,χ). The Mellin transform of x^{σ₀}C + G(χ,x) has nonnegative "mass". By Landau's theorem its abscissa of convergence is a real singularity. Real zeros are excluded separately, so the bound continues to Re s > σ₀. Only a lower bound on G is therefore needed. This is the classical observation behind Ingham's and Landau's treatment of one-sided Chebyshev bounds. Timeline item 5's "Guinand–Weil had the sign" refers to it.

The remark's footnote 3 states the same principle over F_q: if `M^k = qN − R` with qN ≥ 0 or of known sign, and |R| ≤ C^k q^{k/2+O(1)} by k-fold Parseval, then one side of M is bounded by q^{1/2+O(1/k)}. Letting k → ∞ (a tensor-power trick) and using rationality of zeta turns this one-sided bound into RH for curves, as in Bombieri–Stepanov. **I could not find the finite-field toy proof in S2.** Greps for Parseval, path quadric, chain sums `x i * x (i+1)`, Guinand and Weil positivity found nothing matching. The two hits for `x i.castSucc * x i.succ` (`Formalism/HermitianCurvature/HermitianCurvatureProofs/OrbitalIntegralReduction.lean`, `def g1`) are orbital-integral phases. `Development/EtaleCohomology/WeilConjectures/*` (`CurveCut`, `PurityCut`, `FermatCount`, ...) is library material for étale cohomology. Its file names do not suggest a sign-only proof, and I did not read its proofs.

### 1.2 Weil positivity W(1_S) ≥ 0, items 4, 6, 7 [Not found]

I found no functional named W, no test function 1_{(0,1]} + X^η 1_{(0,1/X]}, and no "Stepanov with multiplicity one" determinant. The "interlacing" in item 7 appears in S2 only as Gelfand–Tsetlin/branching interlacing (`ReciprocityTheory/ReductiveGroups/DiscreteSeries/LaplaceOperators/InterlacingWalks/*`, with defs `Interlace`, `InterlaceBlocks`, `HalfIntegralInterlacing`). That is representation theory for the GGP and theta-lift block, not mine. **The directory `PositivityTheory/` is misnamed for this purpose.** It holds "published inputs": Jacquet–Langlands (114 files), cuspidal continuity, Hilbert Eisenstein series, cyclic base change, real Galois representations, a unit-killing lemma (`PositivityTheory/PublishedTheorems/UnitKilling/.../ExistsUnitKillsOfIdempotent.lean`, `exists_unit_kills_of_idempotent`), and an open-cone invariance lemma (`PositivityTheory/Representations/BoundaryCones/.../IsOpenConeInvariant.lean`). It contains no Weil-positivity statement. **[Read, by listing and name grep only]**

### 1.3 FockGram [Read: `Formalism/FockGram/StatementFock.lean`, `StatementSymPos.lean`, `StatementLat.lean`, `StatementTensor.lean`]

**Definitions.**
* `fock f g = Σ_m f_m · conj(g_m) · m!` is the Fischer (Bargmann–Fock) inner product on polynomials, where `m! = ∏ m_i!`.
* `fockGram v` is the Gram matrix of a family v under this inner product.
* `symMono e s = ∏_{i∈s} e_i` for a multiset s of size a.
* `symGram a G = symPow a G · diag(t!)`.

**Statements.**
* `FockPerm`: the Fock product of two products of linear forms is a **permanent**: `fock(∏ x_s, ∏ y_t) = perm(fock(x_s, y_t))`.
* `FockGramSym`: for linear forms e_i, `fockGram(symMono_a e) = symGram a (fockGram e)`.
* `DetSymPow`: `det(Sym^a A) = (det A)^{C(n+a−1, n)}`.
* `DetFockSym`: `det fockGram(Sym^a e) = (det fockGram e)^{C(n+a−1,n)} · ∏_{t∈Sym^a} t!`.
* `SymGramPosSemidef` / `SymGramPosDef`: if G ⪰ 0 (resp. ≻ 0) then `symGram a G` ⪰ 0 (resp. ≻ 0).
* Lattice bookkeeping: `GramChangeDet` (det Gram(v) = |det_b v|² det Gram(b)); `SubIndexZ` (det Gram of a sublattice = index² · det Gram); `DualIndexZ` (an integral Gram matrix has index of the lattice in its dual equal to |det Gram|); `CovolGram` (covolume = √det Gram); orthogonal-sum and tensor-product Gram determinants (`OrthSumGram`, `TensorGram`, `HTmulGram`).

It is used in `Formalism/SharpRestriction/GrowthBounds/FockGramWindow*` and `Formalism/HermitianLattices/*`. **[Inferred]** It computes covolumes and indices of lattices of symmetric-power, theta-type sections, which is Bost-slope and Siegel-mass bookkeeping (timeline items 9, 13, 18). It is a determinant identity plus "PSD is preserved by symmetric powers", not a positivity *inequality* about primes.

### 1.4 BehrendSets [Read: `Formalism/BehrendSets/Statement.lean`, `StatementE.lean`]

**Definitions.**
* `OrdGrp s = (Z/2^{s+2})^× / ⟨−1⟩`, cyclic of order 2^s (`OrdCard`), with `frob s l` the class of a prime l.
* `OrdIff`: ord(frob l) ≤ 2^j iff l ≡ ±1 mod 2^{s+2−j}.
* `ShortOrdSq`, `ShortOrdLin`, `OrdSqJ`: primes whose Frobenius has small order are rare, since `Σ log l / l²` over such l is ≤ 80 Lw² log(2^s)/4^s.
* `BehrendFree K T`: no x₁,x₂,x₃ ∈ T, not all equal, satisfy αx₁+βx₂+γx₃ = 0 with α+β+γ = 0, nonzero |α|,|β|,|γ| ≤ K.
* `BehrendBound`: such T ⊂ [0,N) exists with |T| ≥ N exp(−C√(log N · log(K+1))), which is Behrend's bound.

**Statements.**
* `ThickSet`: a Behrend-free T inside [0, n/8L) with |T| ≥ (n/8L)·exp(−C√(log n log(2L+1))).
* `WindowTwo`: every translate `thick n T u = T + u ⊂ Z/n` meets every progression `t₀ + i·h`, i ∈ [j, j+2L), with addOrderOf h > 2L, in **at most 2 points**.
* `OrbitOne`: if the order of h is ≤ 8L, the orbit meets T+u in at most one point.
* `SpaceFree` / `Spacing` (StatementE): a variant with "j(s₂−s₁) ≠ i(s₃−s₁)" for 0 < i < j ≤ Lw.
* `TranslateMean`: averaging over u, each position lies in T+u for exactly |T| translates.
* Random-subset alteration: `ProbSub`, `ExpectCard`, `TopsBad` (P[t ∈ T and some window W ∋ … contains K+1 points of T] ≤ δ|𝒲|(Lδ)^K), `BadOrbit`, `BadOrbitE`, `Select`, `SelectGen` (there is T with |T| ≈ δn and few bad points).
* `StatementE` also contains `CountIdentity` and `CountQuad`. These are mass-ratio identities for unitary groups of rank k → k+m−1 with Weyl dimension factors `vand/∏ i!`, i.e. the "cost of each lift" (item 14), stored in the same file.

**What Behrend sets are for [Inferred].** They choose a "thick" set of levels or classes in the cyclic 2-power Galois group of the Iwasawa tower (items 11 and 19). Frobenius orbits of primes l (progressions in Z/2^s) then hit the chosen set in at most two points per window. This keeps an averaged sum over the chosen classes nearly diagonal (`ThickPartTranslateMean`, `topExcess` in `SizeWindowBrackets/GrowthBoundGivesRiemannHypothesis.lean`). It is the classical use of Behrend: a large set with no additive coincidences of bounded height, used as a test set. It does not construct extremal objects.

### 1.5 Subset sums and the other combinatorial modules [Read]

* `Formalism/JacquetSubsetSums/JacquetSubsetSum.lean`, `JacquetSubsetSum`. Suppose φ on GL_r(k) transforms on the left under a parabolic with block function b by a representation ρ with no fixed vectors under the relevant block unipotents (a cuspidality hypothesis), and is right-invariant under the unipotent radical of type (a, r−a). If φ ≠ 0, then a is a subset sum of the block sizes.
  * `GoodImpliesH` adds complement symmetry.
  * Use [Inferred]: non-existence of parabolic constituents (cuspidality of a transfer) from local degree data.
* `Formalism/PadicSubsetSums/CentralIdempotents/GaloisIrreducibility.lean`:
  * `DivisorDegreeSubSum`: f has a divisor of degree a iff a is a subset sum of the degrees of its irreducible factors.
  * `SubsetSumCount`: #{I : Σ_{c∈I} deg c = i} ≤ (Σ deg)^i; the reformulated version gives ≤ C(Σ deg, i).
  * `CycleTypeProb`, `TypeCountUpper/Lower`: counts of irreducibles and permutations by cycle type over F_g.
  * `TorusCountBounds`: point counts on tori.
  * This is the finite-field "pigeonhole on F_{q^n}/F_q" toolkit (item 11).
* `Formalism/WindowSums/Statement.lean`: divisibility in the Iwasawa algebra O[[X]] for the distinguished element `ℓε + uX^ℓ` (`R2`, `PowMemIdealOfLe`, `LeOfPowMemIdeal`, `C1`).
* `FrobeniusContraction`: Lubin–Tate/Coleman data `φX = uX(X+ℓ)^{ℓ−1}`.
* `FrobeniusOrbits`: host fields Q(a,b,ζ+ζ⁻¹) with q ≡ 3 (8), d ≡ 7 (16), together with the group G0 × HGrp.
* `SubspaceCounting/Statement.lean`: lengths of quotients of DVR-lattices via "walk" filtrations (`WalkSharp`, `LengthLayers`).
* `GrowthSums/StatementA.lean` and `StatementG.lean`: elementary inequalities, among them `MultinomialVal` (v_p(multinomial) ≤ (k−1)·log_p(Σ)), `ShiftedRows`, `BackwardProduct`.
* `WeightedDeltaMeans/LowerBoundLinearPower.lean`: functionals killing power sums, with derivative identities.
* `ThickAverageLoss/Statement.lean`: `DLVP`, the de la Vallée Poussin bound |Gz y| ≤ K y exp(−c√log y) used as input, plus parameter structures.

None of these counts signed subset sums of real numbers near a lattice.

Remark [Read, style only]: many proofs end in long `first | exact Bl.<name>.goal1 .. | …` cascades into `Lemmata/`, i.e. machine-generated glue. I did not assess correctness.

## 2. Question (a): sign-definite moments over Hecke characters

**Setting.** N is squarefree, N = ∏_{p∈P} p with p ≡ 1 (4) and s = |P|. Write π_p = √p e^{iθ_p}. The lattice points are z = ε√N exp(i Σ_p σ_p θ_p), with ε ∈ {±1,±i} and σ ∈ {±1}^s, so r = 4·2^s points in total. Then

a_{4k}(N) := Σ_{|z|²=N} (z/|z|)^{4k} = 4 ∏_p 2cos(4kθ_p).

This is the N-th coefficient (normalised) of the CM form θ_{λ^{4k}}.

**Proposition 2.1 (moment identity and cluster lower bound) [Proved here].** Let F_K(x) = Σ_{|k|<K}(1−|k|/K)e^{ikx} = K^{-1}(sin(Kx/2)/sin(x/2))² ≥ 0 be the Fejér kernel. Then

Σ_{|k|<K}(1−|k|/K)|a_{4k}(N)|^{2j} = Σ_{z_1,…,z_{2j}} F_K(4(φ_1+…+φ_j−φ_{j+1}−…−φ_{2j})).

Here φ_i = arg z_i. If M points lie in an arc of angular width δ and K ≤ π/(8jδ), every 2j-tuple from the cluster and its unit rotations has |4·(signed sum)| ≤ π/K, where F_K ≥ 4K/π². Hence LHS ≥ (4K/π²)(4M)^{2j}.

**Proposition 2.2 (the exact diagonal) [Proved here].** The signed sum of the angles is ≡ 0 mod π/2 iff Σ_{i≤j}σ^{(i)} = Σ_{i>j}σ^{(i)} in Z^s. This holds because ∏_p (π_p/π̄_p)^{c_p} is a unit only for c = 0, by unique factorisation in Z[i] and π_p ≁ π̄_p. Coordinatewise, the number of ±1 solutions of σ₁+…+σ_j = σ_{j+1}+…+σ_{2j} is C(2j,j). So the exact solutions number 4^{2j}·C(2j,j)^s, and each contributes K. Therefore

(4/π²)·M^{2j} ≤ C(2j,j)^s + O_j(N)/(K·4^{2j}),

where O_j(N) is the F_K-weighted count of off-diagonal near-solutions (signed sums within ~1/K of π/2·Z but not exactly 0).

**Consequence.** Even with O_j = 0, the bound is M ≤ (π²/4)^{1/(2j)} C(2j,j)^{s/(2j)}. For every j ≥ 1, C(2j,j)^{1/(2j)} ≥ √2 (equality at j = 1) and → 2. So the moment method yields at best M ≤ 2^{s/2}, which is exponentially worse than the record's M ≤ 1+ω(N')+K_C (Theorem 1.2). **No j beats Theorem 1.1**, not even the trivial bound in s. The cause is the additive energy of the cube {±1}^s, i.e. the multiplicative structure of the full set of lattice points, which a single arc never sees.

**What would be needed.** To reach M = o(s) one must abandon the global moment and localise it. Then one needs a bound for the number of 2j-tuples *inside one arc*, which is the problem itself. Alternatively, one weights points by a Riesz product c(σ) = ∏ c_p(σ_p). The cluster then needs min weight ≥ (something), and the diagonal becomes ∏_p E_j(c_p). Normalised by (min weight)^{2j}, for c_p(+) = 1 and c_p(−) = t ≤ 1 the per-coordinate ratio is Σ_m (C(j,m) t^m)² / t^{2j} ≥ C(2j,j) at t = 1, and it grows as t decreases. So the bound stays ≥ 2^{s/2} for any weights that keep all signatures. **[Heuristic: checked for the one-parameter family above, not for general weights]**

**Spectral form [Proved here].** Let A_{zw} = F_K(4(φ_z−φ_w)). This is a PSD matrix, since the Fejér kernel is a positive-definite function. Then λ_max(A) ≍ K·M_max with absolute constants: the lower bound comes from the cluster submatrix, the upper from the maximal row sum. A = V diag(ŵ) V*, with V the Hecke-character vectors, so λ_max(A) = ‖T_K(N)‖ for the weighted Toeplitz matrix T_K(N) = (ŵ(k)^{1/2} a_{4(k'−k)}(N) ŵ(k')^{1/2})_{|k|,|k'|<K}. Hence **M_max(N, C) ≍ ‖T_K(N)‖/K with K ≍ N^{1/4}/C.** This is the "iterated Parseval on a path" form: (A^{m})-sums count m-chains of consecutive close points, and 1^T A^{m} 1 ≥ (cK)^m M^{m+1}. It is an exact equivalence, so it gains nothing by itself. Trace bounds tr(T^{2m}) bring back the cube energy, with r^{1/(2m)} → 1 only when m ≫ s.

**Barrier test.**
* *Fake test.* An angle model has the same a_{4k} up to replacing θ_p by real φ_p, and every identity above holds verbatim. Any upper bound on O_j must use the arithmetic of θ_p, i.e. signed sums of 2θ_p at the Liouville threshold, which is the residual problem. Multiplicativity a(Nm) = a(N)a(m) is built into fakes. The integrality N^{2k}a_{4k}(N)/4 = ∏(π_p^{4k}+π̄_p^{4k}) ∈ Z gives only |2cos(4kθ_p)| ≥ p^{−2k} when nonzero, which is single-place Liouville and exactly critical.
* *Record.* This is overview §5.4 rows 5–6: "GAP/Fourier, energy, incidences" reduce to monomial gaps and characters (`round4/outside.md` Thms A.1, A.3), and "L^p, energy …, single-circle data are one Laplace eigenspace, so the estimates are identities" (`local_smoothing/README.md` §2). It is also expert note §5 items 1–2.
* *Verdict.* The sign-definite formulation adds nothing. Capped by fakes for the off-diagonal part, and by the cube energy for the diagonal.

## 3. Question (b): other uses of positivity of G

**Setting** (`round31/ref8_extend.md` §2, `round31/ref9_restate.md` §3). G = Σ_{p,τ} log p · σ_{pτ}σ_{pτ}^T, where σ_{pτ}(i) = +1 if a_{ip} ≥ τ and −1 otherwise. Then:
* G_ii = W;
* G_ij = W − 2d_ij with d_ij = Σ_p |a_ip−a_jp| log p;
* the pair inequality gives G_ij ≤ −h_ij := −4 log|b_ij| + 4 log C'.

**2×2 minors [Proved here].** W² − G_ij² ≥ 0 ⟺ 0 ≤ d_ij ≤ W, which is true by definition. Combined with the pair inequality it gives |b_ij| ≤ C'N'^{1/4}, which is trivial because |b_ij| ≤ |z_i−z_j|. No information.

**3×3 minors [Proved here].** With x_ij = 1 − 2d_ij/W ∈ [−1,1], det ≥ 0 reads 1 + 2x₁₂x₁₃x₂₃ − x₁₂² − x₁₃² − x₂₃² ≥ 0. This says that the arccos x_ij satisfy the spherical triangle inequalities. Because d is a weighted cut (ℓ₁, hypercube-embeddable) semimetric, it satisfies the stronger triangle and perimeter inequalities d₁₃ ≤ d₁₂+d₂₃ and d₁₂+d₁₃+d₂₃ ≤ 2W, and all hypermetric inequalities. All of these hold for **any real exponent profile**, and fakes share the genuine profile, so G-positivity in any minor, Schur complement or symmetric power (cf. `SymGramPosSemidef`) is free and carries no information about the angles.

**Positivity combined with the pair bounds.** Two examples:
* For x ≥ 0, x^TGx ≥ 0 gives W‖x‖² ≥ Σ_{i≠j} x_ix_j h_ij, so W ≥ λ_Perron(H).
* A k×k principal minor with constant off-diagonal −h gives (k−1)h ≤ W. This is weaker than the all-ones vector over the whole cluster (k = M), which is Theorem 1.1.

The record already states that "Theorem 1.1 and VIII.1 are two certificates for the same PSD system", and that the optimum is the Perron root, which exceeds the average row sum by O(log M) only for irregular collision patterns (`round31/ref8_extend.md` §2). `round31/ref9_restate.md` §3 shows that positivity on all-ones permits mean log|b_ij| = O(log M) in every fake class, and that "no other eigenvector or minor of G can force a single large real cofactor either".

**"Stepanov with multiplicity one" reading [Speculative].** In the remark, a determinant is nonzero and integral, so it is ≥ 1, while a Hadamard bound caps it. For us the integral Gram matrices available are:
* Re(z_i z̄_j), of rank 2;
* the Hermitian rank-one matrix z_i z̄_j;
* determinants of monomials in the points, i.e. route 5.

G itself has entries in Z-combinations of log p and is not integral, so the integrality-versus-Hadamard comparison does not apply to G. The determinant-of-monomials version is route 5, which is capped at the zero estimate (expert note §5 item 3). **Verdict: nothing new.**

## 4. Question (c): Behrend sets and our problem

1. **Building genuine clusters.** The angles of a genuine cluster are fixed by θ_p. Choosing signatures from a Behrend-type set, or from a set with no affine 3-flats in F_2^s as in IX.1, does not place their angle sums in a common window. Large genuine clusters are not constructible by additive selection. The record's refutation barriers apply (expert note §5 item 9).
2. **Building fakes.** Fakes need only a profile and real angles solving a linear system with gaps. Paley/Legendre profiles already give W ≤ (23/3+o(1))M log M inside the window (expert note §1). A Behrend-type profile could only change the constant. It cannot lower W below O(M log M), and that already caps the rate.
3. **Choosing parents or test sets, as in the Lean.** The analogue of `WindowTwo` and `TranslateMean` is to pick a subset P' of the primes, or a set of inert moduli ℓ ∈ S for Theorem 1.2, such that additive coincidences of bounded height among (2θ_p)_{p∈P'}, or among residue classes, are rare. This helps only for chosen parents, and the record already has the random-parent rate (O(log s) at the smooth scale, `round31/README_XI.md` §5(i)). For the worst case, and for q_P, which uses *all* primes in [P,2P], one cannot choose. For S in Theorem 1.2, a better S only improves K_C, i.e. the constant in M ≤ 1+ω(N')+K_C, not the rate.
4. The Lean's random-alteration lemmas (`TopsBad`, `SelectGen`) are standard deletion-method statements. A use for us would be to select, inside a putative cluster, a large sub-cluster whose pairs avoid some bad event. Every pair-level event is fake-satisfiable (expert note §5 item 1).

**Verdict:** no use for the rate.

## 5. Borrowable ideas, tested

| # | Idea from S2 / remark | Closest analogue for arcs | Fake test | Char p | Record | First step | P(rate) |
|---|---|---|---|---|---|---|---|
| 1 | One-sided bound suffices by positivity (Landau, `LandauStep`) | an upper bound on a nonnegative count; already one-sided | n/a | n/a | the count is already one-sided; no analytic family in which a one-sided bound continues | none | ~0 |
| 2 | Tensor-power amplification M^k = qN − R | product cluster: k-fold products of cluster points on \|z\|² = N^k | fake-transparent (angles add); reduced chords of product pairs equal base reduced chords, so no new arithmetic | works in char p too (F_2(t) clusters power up) | Liouville/route-5 multi-point = one-point toric Liouville (expert note §5 item 3) | check that the k-power cluster sits at scale N^{k/2−1/4} ≫ (N^k)^{1/4}, so it is non-critical (done, §5.1) | <0.5% |
| 3 | Iterated Parseval / PSD kernel over the family of Hecke characters | M_max ≍ ‖T_K(N)‖/K, a Toeplitz matrix of a_{4k}(N) (§2) | identical for fakes | the analogue holds for F_q[t] theta series | overview §5.4 rows 5–6 | try to bound ‖T_K(N)‖ using that a_{4k}(N) are CM coefficients **jointly in N and k** (Rankin–Selberg over weights) | <0.5%; averaged-over-N versions are capped by grafts (expert note §5 item 2) |
| 4 | Behrend thick sets in the tower group | choose S or primes avoiding coincidences | constants only | — | random-parent rate, Theorem 1.2's K_C | none worth doing | ~0 |
| 5 | Fock/Gram determinant of symmetric powers (multiplicity via Sym^a) | interpolation determinant on monomials in the points | — | — | route 5, capped at the zero estimate (round 26 `ze.md`, round 27 `interp.md`) | belongs to the multiplicity block | (not this block) |

### 5.1 The power trick, worked out [Proved here]

Let z_1..z_M lie in an arc of angular width δ = CN^{−1/4}. For multisets I of size k, the products z_I lie on |z|² = N^k in an arc of angular width kδ, i.e. of length kC·N^{k/2−1/4}. The critical length on that circle is (N^k)^{1/4} = N^{k/4}, and k/2 − 1/4 > k/4 for k ≥ 2. So the product cluster is far above the critical scale, where no nontrivial bound is known or expected (Cilleruelo–Córdoba bounds need length R^{1/2−ε} or less).

For I, J differing in one element, z_I − z_J = (z_i − z_j)·z_{I∩J}, and the reduced chord is c_ij. So the product configuration carries exactly the pair data of the base cluster. The amplification that works over F_q works there because the error R is a complete sum with exact orthogonality, which has no analogue at the archimedean place at scale N^{−1/4}.

## 6. Bottom line

Nothing in this block of S2 transfers to the short-arc problem:
* The Lean's "sign" is Landau's theorem for nonnegative Dirichlet coefficients, with an apparatus of averaging over Behrend-thick sets in a 2-adic Iwasawa tower.
* The path-quadric toy proof and the Weil-positivity determinant are not in S2 under any name I could find.
* `FockGram` is Gram-determinant bookkeeping.
* The subset-sum modules concern degrees of factors and parabolic types.

On our side, the sign-definite moment over Hecke characters is an exact restatement. Its exact diagonal (the cube energy C(2j,j)^s) already prevents any bound better than 2^{s/2}, and its off-diagonal is the residual problem, which fakes satisfy. Positivity of G on minors, Schur complements or symmetric powers is automatic for every exponent profile, and combined with the pair bounds it reduces to the Perron-root certificate already in the record. Behrend-type selection affects constants or chosen parents only.

Honest probability that anything from this block leads to M = o(log R/loglog R), even for q_P: below 1%.
