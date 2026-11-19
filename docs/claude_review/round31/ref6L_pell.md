# Referee report, round 31, block "pell": Astra Part VI, section VI.4 (four-point constructions)

Source: `astra_lemma_VI/answer.md` lines 1457-1488, together with `vi_pell_four.py`.
Scripts: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref6L_pell/scripts/`
(`pell4.py` + `pell4.out`, `gcdcheck.py`, `minconst.py`, `closedK.py`, `fib.py`). All of them use exact Gaussian
integer arithmetic. Angles are computed with mpmath at 1500 digits. The CG 2007 text (arXiv math/0608109) is in `ref6L_pell/dl/`.

## Verdict table

| Claim | Verdict |
|---|---|
| (1,1,1) Pell family: four points, equal norms, Δ ≍ H^{-1/3} (arc ≍ R^{1/3}), bounded gcd via the triangle bound | CORRECT |
| Three-point sharpness: Δ ≥ 2^{1/3} H^{-1/3} for any three points (area ≥ 1/2 = abc/4R) | CORRECT (classical, Cilleruelo/Jarník) |
| (6,7,7) family: Δ ≍ H^{-3/10} against primitive norm ≍ λ^{40n} (arc ≍ R^{2/5}); gcd analysis with D | CORRECT |
| Zariski density / escape (rational functions of λ^{-2n}, pair-product map with |det| = 2) | CORRECT |
| Fibonacci one-unit-class variant (β = -1-2i, primary factors) | CORRECT |
| Extension test: endpoint five-tuples force w_ij = W/10, all ten triples sharp, four-mark profile (W/10)^4 (W/5)^3 | CORRECT |
| The four-point families cannot be extended at the endpoint (profile mismatch) | CORRECT (it even follows from the exponent alone) |
| Mixed-rate block-1 overlap is bounded "for generic seeds" | CORRECT WITH CORRECTIONS (true for **all** seeds for the 6-vs-7 pairs; genericity is needed only for the equal-rate pair) |
| "O_Z(R^{1/3}) with a constant depending on Z" vs Cilleruelo-Granville | CORRECT, no conflict. The smallest constant in the family is ≈ 4.6015 > 4.392 > 4.347 |

## (a) The construction, line by line

* A = [[5,2],[2,1]] is symmetric with det 1 and trace 6, so its eigenvalues are λ^{±1} with λ = 3+2√2. The λ-eigenvector has
  slope λ-5)/2 = √2-1 = tan(π/8). The two eigenlines are orthogonal. A = 3I + 2D with D = [[1,1],[1,-1]], and
  D(x,y) = (x+y, x-y) is the column of (1+i)(x-iy). So D commutes with A. Checked.
* The four points are P0 = Ḡ1Ḡ2Ḡ3 and P_jk = -i·G_jG_kḠ_l, so P_jk/P0 = -i g_j g_k. Each g_j tends to e^{iπ/4}, so
  -i g_j g_k tends to 1 **exactly**: the attracting phase is cancelled by the unit -i, with no approximation. All norms equal
  ∏|G_j|^2. Checked exactly.
* Asymptotics. Write G = αλ^n e₊ + βλ^{-n} e₋ (real α, β, with αβ = (√2/4)Q(v) and Q = y²+2xy-x² ≠ 0 on rational seeds).
  Then arg G_j = π/8 + u_j x + O(x²), where u_j = β_j/α_j and x = λ^{-2n}. The ratio angles are 0 and 2(u_j+u_k)x. So for
  generic seeds (u_j distinct and no u_j+u_k = 0) all six chords are ≍ R·λ^{-2n}, against R² ≍ λ^{6n}.
* gcd. For a triple, abc ≍ R³λ^{-6n} ≍ R. After dividing by g, a'b'c' = abc/|g|³ ≥ 2R/|g| gives |g|² ≤ abc/(2R) = O(1). Correct.
* (6,7,7). The block-1 deviation is ≍ λ^{-12n} and the others are O(λ^{-14n}), so the spread is ≍ λ^{-12n}. The raw norm is
  ≍ λ^{2(6+7+7)n} = λ^{40n}, which gives Δ ≍ H^{-12/40} = H^{-3/10} once the gcd is bounded. The triangle bound only gives
  |g|² ≪ λ^{2n} here, so the separate argument is needed. I checked it step by step:
  * Primitive columns have odd prime divisors that split, and the 2-part of their norm is ≤ 2.
  * For the same orientation, p^a | Im(Ḡ2G3) = [G2,G3] = [v2,v3], which is A-invariant because det A = 1.
  * For opposite orientations, apply the same argument to DG3, giving [v2,Dv3].
  * The four words realise all four conjugation patterns on any two blocks. So v_π(gcd) ≤ min(v_p N G2, v_p N G3).
  * The (1+i)-valuation is at most 3.
  * The conclusion holds: the gcd is bounded by a constant depending on the seeds. In 40 random seed triples and
    n ≤ 40 (`gcdcheck.py`) the odd part of gcd(N2,N3) always divided [v2,v3][v2,Dv3], and N(Gaussian gcd) ≤ 26.
* Escape. With fixed integer rates, g_j = e^{iπ/4}(α_j + iβ_j x^{r_j})/(α_j - iβ_j x^{r_j}) is rational in x = λ^{-2n}.
  * Pull back a defining polynomial of Z. If it is nonzero at x = 1 (n = 0), it has finitely many zeros.
  * Rational slopes give an infinite, hence Zariski-dense, set in S¹. The product of three such sets is dense in (S¹)³.
  * The pair-product map has exponent matrix [[1,1,0],[1,0,1],[0,1,1]] with det -2. It is an isogeny, so it preserves density.
  * The nondegeneracy conditions are closed.
  * The argument is correct, and the same argument works for (6,7,7) (powers x^6, x^7).
* Fibonacci variant.
  * tan 2θ = 2 when tan θ = 1/φ, so e^{4iθ} = (-3+4i)/5 and β̄/β = (-3-4i)/5. Their product is 1.
  * F_{6m} ≡ 0 and F_{6m+1} ≡ 1 mod 4, so every factor, conjugate and β^{±} is primary.
  * Binet gives G = φ^{6n}(A_d - φ^{-12n}B_d), which is rational in φ^{-12n}.
  * The block norms are F_{12(n+d)+1}, whose pairwise gcds divide F_{gcd} with gcd | (d_j - d_k). So they are bounded.
  * All points checked.

## (b) Own reproduction (exact; `pell4.out`, `fib.py`)

* `vi_pell_four.py` runs and passes. Its scale assertion is loose (factor up to 2^{40}), and it uses a single seed triple.
* My code, seeds (3,1),(4,1),(5,1), rates (1,1,1):
  * N(gcd) = 4 for all n.
  * log Δ / log H' goes -0.300, -0.312, -0.324 at n = 5, 8, 20. The convergence to -1/3 is slow, but arc/R'^{1/3} is
    exactly stable at 7.0128 from n = 2 on.
  * Cut profile: three 2|2 cuts of weight ≈ 0.336, 0.334, 0.330 (each → 1/3), and singleton weight 0. At some n there are
    transient singleton layers of weight ≤ 0.09 that decay.
* Seeds (2,1),(7,3),(1,3): constant 16.562.
* Rates (6,7,7), three seed triples:
  * log Δ / log H' = -0.29991, -0.30016, -0.29937 at n = 20 (primitive norms about 2040 bits).
  * arc/R'^{2/5} is stable at 1.1408, 0.7942 and 2.4310.
  * N(gcd) ∈ {2, 4, 20}.
  * Profile: 2|2 weights ≈ (0.35, 0.35, 0.30) = (14, 14, 12)/40, singleton weight ≈ 0.
* Fibonacci, five offset triples, n ≤ 12: gcd 1, all points and quotients primary, arc/R'^{1/3} constant (494.0 for d = (0,1,2)).

## (c) Literature, and the constant

* **CG 2007 (arXiv math/0608109, §5, after Thm 13), quoted from the PDF:**
  * Cilleruelo-Córdoba: at most k points on arcs of length R^{1/2-1/(4[k/2]+2)}. This is "sharp for k = 1,2,3".
  * Arcs of length (16R)^{1/3} contain ≤ 2 points.
  * "an arc of length (40+(40/3)√10)^{1/3} R^{1/3}, with R > √65, contains no more than three lattice points, whereas there
    exists an infinite family of circles … containing four lattice points on an arc of length
    (40+(40/3)√10)^{1/3} R_n^{1/3} + o(1)".
  * "Other than in the examples arising from this family, an arc of length (40+20√5)^{1/3} R^{1/3} contains no more than three
    lattice points". There is an explicit Fibonacci four-point family attaining (40+20√5)^{1/3}.
  * Item (4): "Theorem 13 … implies that an arc of length R^{2/5} contains at most 4 lattice points … Are there infinitely many
    circles with four lattice points on an arc of length ≪ R^{2/5}?"

  Read literally, item (4) is answered by their own item (3). It must mean five points, which is the record's reading in
  overview §1.3. Astra's R^{2/5} four-point family therefore answers nothing open.
* **CG 2009 (CJM 61).** The full text could not be fetched (cambridge.org and umontreal.ca are blocked). Search snippets say
  it classifies four-tuples on arcs of length tR^{1/3}, giving finitely many families for each t (Cor. 1.7: six exceptional
  families for R > √325). **Unverified.** If that reading is right, the four-tuples with constant ≤ t lie on finitely many
  curves in T_4. Escape from every Z then **forces** the constant to grow with Z, so Astra's O_Z is necessary, not just
  convenient.
* **No conflict with CG.** Four points on arcs of order R^{1/3} are **known** (CG). Over seeds with coordinates ≤ 10 (all
  336,588 triples screened, exact gcd for every candidate whose g = 1 constant is < 9.3; `minconst.py`), the smallest limiting
  constant of the Pell family is
  `K_min = ((320+256√2)/7)^{1/3} ≈ 4.60154`. It is attained, with gcd 1, for example at seeds (1,0),(0,1),(2,-1) (`closedK.py`).

  This exceeds both CG thresholds, 4.34737 and 4.39202, as it must. Caveat: triples with N(gcd) ≥ 16 and g = 1 constant ≥ 9.3
  were not checked exactly; I saw no such gcds. I did not determine whether the Pell family is one of CG09's families.
* **Novelty.** Four-point families at R^{1/3} are in CG and in the record (the canonical four-point Pell family, `growth/`).
  What is new relative to the record is the Zariski-dense / escape form. The record had A_4(δ) false only for δ < 1/14
  (round 20, Thm 4.1, exponent 3/7). VI.4 gives **A_4(δ) false for every δ < 1/6**. For δ ≥ 1/6, A_4 holds vacuously for large
  R, since √2 R^{1/3}-arcs carry ≤ 3 points (CC92). So the four-point essential exponent equals the Liouville/CC exponent,
  which settles the k = 4 analogue completely. The mathematics is elementary (a product of three Pell orbits with a
  unit-cancelled phase) and is probably folklore-adjacent. I found no published Zariski-density statement.

## (d) Extension test and generic-seed gcd

* The chain |q_i-q_j| ≥ e^{-d_ij/2} → Σ d_ij = 6W - 2W_sing → W_sing = o(W) → d_ij = 3W/5 → w_i = 2W/5 → w_ij = W/10 checks
  line by line.
* Deleting one mark gives four singletons of W/10 and three 2|2 cuts of W/5, with primitive norm H^{1+o(1)}. Deleting two
  marks gives primitive norm H^{9/10}, with chords H^{3/20} = R_3^{1/3}, so all ten triples are sharp. Correct.
* Mismatch. The Pell families have profile (singletons ≈ 0; 2|2 ≈ W/3 each), and the mixed family has (≈ 0; 0.30, 0.35, 0.35).
  Both differ from (W/10 ×4; W/5 ×3). More simply, an endpoint four-subtuple has Δ ≍ H₄^{-3/10}, not H₄^{-1/3}. Ratios
  determine the primitive tuple up to a unit, so "retaining ratios" is meaningful. Correct.
* **Correction to the generic-seed claim.** The norm Laurent polynomials are f_j = a_j t^{2r_j} + b_j t^{-2r_j} with
  a_j = α_j², b_j = β_j² > 0. Suppose f_1 and f_j (r = 6 and 7) had a common root. Then t^{24} = -b_1/a_1 < 0 and
  t^{28} = -b_j/a_j < 0, so t^4 > 0 and t^{24} = (t^4)^6 > 0, a contradiction. So the polynomials are coprime for **every**
  seed, and Bézout bounds gcd(N1, Nj) unconditionally. Genericity is needed only for the equal-rate pair, where it is supplied
  anyway by [v2,v3][v2,Dv3] ≠ 0. Numerically, gcd(N1, Nj) ≤ 7066 for n ≤ 40 in 40 random seed triples, while the norms grow
  like λ^{12n}.

## Standing cautions

* No Subspace theorem is used, so no exceptional sets are involved.
* Gaussian gcds are bounded explicitly, by the triangle bound or by cross-products. They are not discarded.
* The common phase is cancelled exactly by -i (or by β̄/β), not set to zero.
* Constants depend on the template (the seeds), and this is stated honestly. It does not affect exponent statements. By CG09
  it is apparently unavoidable.

## Significance

The section is correct and modest. It closes the k = 4 analogue of Question B **negatively**: for four points, the Liouville
exponent is essentially sharp. That is a caution against expecting an essential improvement at k = 5 by analogy. It says
nothing about (*), and Astra does not claim it does. The R^{2/5} four-point family is redundant: R^{1/3} families already exist
in CG.
