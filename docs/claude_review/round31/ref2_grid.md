# Referee report, round 31, block "grid": astra_free_II/answer.md §II.3 (grid determinants)

Referee scripts and outputs: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref2_grid/`
(`symbolic.py`, `table.py`, `cube_mult.py`, `cyc_width.py`, each with an `*_out.txt` output). They were written independently
of Astra's `ii_grid_check.py`. I also reran Astra's script, and its output agrees.

Setup used throughout (from Part I §2). After one-unit shortening and gcd removal,
z_j = ε ∏_p π_p^{a_jp} π̄_p^{e_p−a_jp}, W = log N, d_ij = ∏ p^{|a_ip−a_jp|}, and log d_ij ≥ W/2 + log 16 on every pair.
That last bound is the "pair bound": it is the cutoff 16√N ≤ d in (1).

## Verdicts, claim by claim

### 1. Normalization of κ_raw — CORRECT
Each permutation term of F_G = det((z_g)^α) is a monomial. Its π_p-valuation is Σ_g c_{g,σ(g)}, and its π̄_p-valuation is
e_pV − Σ c. So the forced norm exponent is at least e_pV + L_p − U_p, which accounts for both places. Collecting
monomials shrinks the support, so the forced exponent can only grow. Translation invariance holds: a common shift
transforms the columns of a lower set unitriangularly, so F_G(z+t) = F_G(z). With homogeneity, the multiplicity at
every diagonal point is exactly V. The archimedean bound is |F_G|² ≤ c_G² (RΔ)^{2V} = c_G² C^{2V} N^{V/2}, so the
critical value is κ = 1/2, and c_G is constant for fixed d. Astra's caveat is right: these numbers are forced
monomial divisibility, not actual valuations.

### 2. Finite-size Vandermonde — CORRECT
With s zeros and binary e_p = 1, I get L = C(n−s,2) and U = C(n,2) − C(s,2), so the raw cost equals the collected cost,
C(s,2) + C(n−s,2). Half-half at even n gives (n−2)/(2(n−1)). Fair bits give exactly 1/2.

### 3. Cyclic six-mark design — CORRECT
- **Circulant identity.** I expanded the 6×6 determinant symbolically in generic z_0..z_5 (198 monomials) and compared it
  with ∏_{k=1}^5 Σ_j ξ^{kj} z_j over Q(√−3). The two agree exactly with sign −1 (`symbolic_out.txt` (1)). The proof is
  sound: the circulant C = (z_{i+h}) has column sum S·1, so det C = ±S·det M, and det C = S·∏_{k≥1} f_k. Cancel S in the
  polynomial ring.
- **Summation by parts.** Verified symbolically mod h^6 − 1.
- **Nonvanishing.**
  - The nonzero h^j − 1 have arguments π/2 + ψ/2 with ψ ∈ {π/3, …, 5π/3}. That is a sector of width 2π/3.
  - The increments z_j − z_{j−1} have arguments (θ_j + θ_{j−1})/2 + π/2. Their spread is (θ_5 + θ_4 − θ_1 − θ_0)/2 < θ_5 − θ_0 < π/3.
  - Adding arguments, the products lie in a closed sector of width < π, hence in an open half-plane. The j = 1 term is nonzero.
  - A random numerical check (20000 ordered sextuples of width < π/3) found no near-zero factor.
- **Exact cost.** The cyclic diagonal argument is right: choose j* and give column h to row j* − h + 1. This yields
  L_p = 5 min a and U_p = 5 max a, so κ_raw = 1 − (1/W) Σ_p (max − min) log p. The coefficient of z_j^5 in the product
  is ±1 (verified), so collected = raw. On all 64 binary patterns I computed the true support exponent from the expanded
  polynomial: it is 5·[pattern constant].
- **Correction (minor).** Astra writes "on exactly three-three cuts κ = 0". In fact κ = 0 on *every* nonconstant cut,
  not only three-three cuts. The fair-bit expectation 1/32 is right; I reproduced it as E(m_S + m_{S^c}) = 5/32.
- **Actual profile.** Σ(max − min) log p ≥ log d_ij ≥ W/2 + log 16, so κ ≤ 1/2 − log16/W. This is correct, and in
  general far from tight: κ = 1 − range/W.
- **Conclusion.** The design really has no pair factors and really has free nonvanishing. It is strongly subcritical:
  0 on the Walsh model, and at most 3/31 on 200 random 6-label sets.

### 4. "Every total degree" design with one non-pair factor — CORRECT
- **Structure.** In the Newton basis, N_α(row b) ≠ 0 forces α ≤ b componentwise. The two removed indices are the only
  rows dominating them. The surviving block is lower triangular, and the two new rows meet the two special columns in
  a 2×2 block.
- **Block formula.** I checked P(u)P(v)(u−v)(u+v−a−b) with a = s_{1,d−1}, b = s_{2,0}. It holds symbolically for
  d = 1..5, with sign +.
- **Whole determinant.** For d = 1 the full symbolic determinant factors as
  (s2_0−s2_1)(s3_0−s3_1)(s4_0−s4_1)(u−v)(u+v−s0_0−s1_0), up to sign.
  For d = 1, 2, 3 I evaluated the full 6×6, 21×21 and 56×56 determinants at random integers in exact arithmetic.
  Each equals −(predicted product), and the predicted product has exactly V linear factors (5, 35, 140).
- **Midpoint argument.** It is correct: a chord's midpoint M ≠ 0 determines the chord, and short arcs contain no antipodes.
- **Fair-cut value.** ((V−1)/2 + 1/8)/V is correct. I reproduced it as E(m_S + m_{S^c}) = 17/8 at d = 1, and the
  uniform half-cut replacements are right.
- **Actual profile.** Each linear factor contributes W − range(its marks) ≤ W/2 − log 16. Correct.

### 5. Model table, d = 1..4 — CORRECT (all four rows reproduced)
- **Method.** I used my own design generator. Assignment optima came from an LP over the assignment polytope (HiGHS),
  cross-checked against LSA. For d = 1 I also brute-forced all 720 permutations.
- **Collected column.** It is computed from the factorization. This is exact because a linear functional's minimum over
  the Newton polytope of a product is additive (Minkowski sum). For d = 1 it also matches the support of the expanded
  determinant.
- **Results.**
  - Raw: 61/155, 464/1085, 127/310, 1259/3255.
  - Collected: 67/155, 513/1085, 522/1085, 524/1085.
  - The collected values approach the pair-monochromatic rate 15/31 ≈ 0.4839 from below.
- **Caveat.** This is a model only. It has 31 equal-weight cuts and is not an actual profile, which Astra states.

### 6. Lower-simplex obstruction, min over nonsingular minors κ ≤ 1/2 − log16/W — CORRECT
The lower simplex design on X^5 (d < m) is nonsingular, and its determinant is a product of V pair differences, each
contributing ≤ 1/2 − log16/W. With prime powers each pair contributes e_p − |a_i − a_j|. So the minimum is at most this
value. The logic is trivial once the design exists, but it is the right observation: Cauchy–Binet existence of *some*
nonzero minor cannot give supercriticality. The same holds for the Gram determinant Σ_G |det_G|², whose forced
exponent is at least the minimum over G and no better. The distinct-marks variant needs m ≥ 5(d+1); correct.

### 7. Fair-cut cap — CORRECT, but not new
- **Multiplicity at the vertex.** The bound mult(F, x_S) ≥ m_S + m_{S^c} is correct:
  - The first ideal comes from the support condition. The second comes from translation invariance, F(x) = F(x − 1).
  - The two are monomial ideals in disjoint variables, so the lowest homogeneous component lies in I^a J^b.
- **The averaging lemma.** It is the Dvir–Kopparty–Saraf–Sudan multiplicity Schwartz–Zippel lemma: the expected
  multiplicity at a random point of S^n is at most d/|S|, here with |S| = 2 (SIAM J. Comput. 42 (2013); arXiv:0901.2529;
  the abstract states it in this form). Astra's induction sketch is exactly the DKSS proof, and I checked it step by step.
- **Counterexample search (none found).**
  - 286 random dense polynomials: maximum average/degree = 1/2, attained only in tight cases such as x^a(x−1)^b.
  - 400 random sums of products of affine forms: maximum 1/2.
  - Hill climbing for m = 3, deg ≤ 4: best 7/32.
- **Novelty.** None in substance. The record already has this as `threshold/construct.md` Theorem B and Cor. C3 (the cube
  profile satisfies the forced-divisibility Liouville inequality of *every* F), and identifies it as the Cayley
  transform of DKSS (`round27/kill_stab.md`). Astra's version is a direct Boolean-cube restatement for translation-invariant F.
  Astra correctly notes that it concerns independent fair cuts, not a specified sparse family. It also covers only e_p = 1.

### 8. Laurent bookkeeping — CORRECT
- **Scaling and forced exponent.** z^{−1} = z̄/N. Every term carries total negative exponent Q, so N^Q F is integral and
  |N^Q F|² ≤ N^{T+2Q} Δ^{2v}·coef. The forced exponent is e_p(T+2Q) − width_p. Supercriticality ⇔ Σ_p width_p log p < vW/2.
- **Relation to the record.** This is the record's mean-width criterion (Theorem B: reg ≤ 2E width). No new content beyond the bookkeeping.
- **One remark.** For Laurent spaces, translation invariance and hence the automatic order v are lost. Astra treats v as
  a parameter, which is fine.

## What remains open, and my assessment

Open: a design G that is **tailored to the specified sparse prime profile** (about W/log W cuts, far fewer than 2^m,
so neither the Thm B/C3 average nor the cap applies), with F_G(z) ≠ 0 *at the actual point* and Σ_p width_p(supp F_G) log p < VW/2.

I see no reason to expect a supercritical design with *free* nonvanishing, for these reasons:
- Every mechanism of automatic nonvanishing exhibited (triangular Newton designs, the 2×2 block, the cyclic half-plane
  argument) makes F_G a product of linear forms. Each linear form with sum of coefficients zero has width at least one
  pair distance at every prime, so it contributes at most W/2 − log 16. **Any product of such linear factors is
  subcritical on every actual profile.**
- So a supercritical design needs an irreducible factor of degree ≥ 2. Its nonvanishing at the actual point is no
  longer free: it is a zero estimate at the point.
- That is precisely route 5's open lemma (ZE^all_H; the record already has supercritical forced vanishing for boxes
  chosen after the profile).
- So the grid setting re-poses route 5's missing step; it does not bypass it.

## Novelty and significance
- **New relative to the record.** The record's interpolation-determinant reduction (`growth/auxiliary.md` Prop. A,
  overview row 16) is one-variable. The X^5 grid versions are new: the cyclic non-pair design with exact cost, the
  one-non-pair-factor design, and the "min over minors" observation.
- **Not new.** The fair-cut cap (Thm B / C3 / DKSS) and the Laurent criterion are already in the record.
- **Significance: peripheral.** Nothing here changes the growth rate (still (2/3+ε) log R/log log R), the uniform bound,
  or Questions A and B. For route 5 and the coordinator's Question 4, it is a clean *negative* answer for the natural
  free-nonvanishing designs: there are no pair factors, but the designs are still subcritical. Route 5's status is unchanged.

Corrections to the text: (i) "on exactly three-three cuts κ = 0" should read "on every nonconstant cut κ = 0";
(ii) the cap should be credited to DKSS and the record's Theorem B / C3 rather than presented as new.
