# Referee report, block "baker" (round 31)

**Scope.** I checked `astra_lemma/answer.md` §2.2, the "conclusion" paragraph of §2, and §2.3, together with the
uniform quantifiers of the §2 theorem. §2.1 (the McKinnon–Roth step on a fixed cover) is outside this block. The
theorem inherits it unchanged; see item 9.

**My scripts** are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref_baker_checks/`:
* `check.py` covers §2.2. It checks 271 random full-rank primitive five-tuples with up to 7 split primes and
  index up to 37500, plus a family y_j = x_j^N with N up to 21 and index up to 194481.
* `check23.py` covers §2.3. It checks the 7e_1 example, 9 tuples in case B and 29 in case C.

Every identity below is checked in exact Gaussian-rational arithmetic, and each inequality at high mpmath precision.
I did not just rerun Astra's scripts.

## Verdicts

### 1. Minkowski step and m ≤ I: CORRECT (with a harmless sharpening)
* L* has covolume 1/I. A cube of half-side t has volume (2t)^4, and Minkowski needs (2t)^4 ≥ 2^4/I. So t = I^{-1/4}
  already suffices; Astra's 2I^{-1/4} is weaker but valid.
* In every test, the LLL-plus-enumeration vector satisfied ||v||_∞ ≤ I^{-1/4}.
* L*/Z^4 is dual to Z^4/L, so it has order I. Hence the order m of v divides the exponent of that group, and m ≤ I.
  I checked I ≡ 0 mod m in every case.
* When I ≤ 16, v may lie in Z^4, giving m = 1. Nothing breaks.

### 2. u integral, nontorsion, and 2h(u) ≤ 4||v||_∞ W: CORRECT
* <v,b_p> is an integer because v ∈ L* and b_p ∈ L.
* If all <v,b_p> were 0, then v would be orthogonal to a rank-4 lattice, so v = 0. Hence some exponent is nonzero.
* The g_p = π_p/π̄_p come from distinct non-associate primes. Unique factorisation therefore makes u non-torsion
  (u has a nonzero π_p-valuation).
* 2h(u) = Σ|<v,b_p>| log p ≤ ||v||_1 Σ e_p log p ≤ 4||v||_∞ W. This uses ||b_p||_∞ ≤ e_p and Σ e_p log p = W, which
  needs primitivity. Checked numerically.

### 3. (y^c)^4 = u^{4m} and |u^{4m}−1| ≤ 4||c||_1 Δ: CORRECT
* y^c = ε^c ∏ g_p^{m<v,b_p>} = ε^c u^m with ε^c ∈ μ_4. So (y^c)^4 = u^{4m}; verified exactly.
* For unit-modulus a, b, w and any integer n: |ab−1| ≤ |a−1| + |b−1| and |w^n−1| ≤ |n||w−1|. With |y_j−1| ≤ Δ (since
  y_5 = 1), this gives |(y^c)^4−1| ≤ 4||c||_1 Δ.
* No phase is discarded: Δ is a diameter of ratios, and the common phase drops out exactly.

### 4. Baker–Wüstholz form (3): CORRECT (cited statement checked only through a secondary source)
**Statement used.** If Λ = b_1 log α_1 + … + b_n log α_n ≠ 0, with α_i in a field of degree d, then
log|Λ| > −c(n,d) h'(α_1)···h'(α_n) log B.
* h'(α) ≍ max(h(α), |log α|/d, 1/d), and B = max|b_i| (at least e).
* c(n,d) = 18(n+1)! n^{n+1} (32d)^{n+2} log(2nd) in the usual citation, which I give from memory.

**Application here.** n = 2, d = 2, α = (u, −1), principal logarithms, so |log u| ≤ π and h'(−1) is an absolute
constant.
* Λ = i(nθ − 2πk), where θ = arg u and k is the nearest integer to nθ/2π. Then |Im Λ| ≤ π and
  |u^n−1| ≥ (2/π)|Λ|.
* Λ ≠ 0 because u is non-torsion.
* B ≤ n + 2, so the bound is −log|u^n−1| ≪ (1+h(u)) log(2n) with an absolute constant.
* BW does not need the α_i to be multiplicatively independent, only Λ ≠ 0. The k = 0 case is covered by the same
  two-logarithm form or by Liouville.

**Source check.** I could not open the Crelle paper. A secondary source states exactly this specialisation:
−log|α^b−1| ≪ D^4 log D · h'(α) log b (arXiv:1611.04664, which cites Baker–Wüstholz 1993). The record's
`round8/lw.md` line 126 quotes the same product-of-heights shape.

**Sanity check.** In my samples, the ratio −log|u^{4m}−1| / ((1+h) log 8m) is at most 0.25. This is not a proof.

### 5. I ≤ 24(W/log 5)^4 and ||c||_1 ≤ 8I: CORRECT
* I is the gcd of the 4×4 minors, so it divides any nonzero minor. A minor is at most 4! ∏ e_{p_k}, and
  e_p ≤ W/log p ≤ W/log 5. I checked that the gcd of minors equals the HNF determinant.
* ||c||_1 = m||v||_1 ≤ 4m · 2I^{-1/4} ≤ 8I^{3/4} ≤ 8I.

### 6. Inequality (2): CORRECT, with absolute constants
The chain is

  −log Δ ≤ log(4||c||_1) + C(1 + 2||v||_∞ W) log(8m) ≤ 4C I^{-1/4} log(8I) W + O(log(2+W)),

using log I ≪ log(2+W) from item 5. Dividing by W gives (2). The factor 2 in ||v||_∞ only changes C.

Consistency check: on the family y_j = x_j^N, the index is N^4 I_0, h(u) stays fixed while W grows linearly, and
(−log Δ)/W → 0, as (2) predicts.

### 7. The conclusion's split [d0, I_*] ∪ (I_*, ∞): CORRECT
* **I_*(g) is finite and independent of the tuple.** The function I^{-1/4} log 2I decreases for I > e^4/2 and tends
  to 0, and C is absolute. Define I_*(g) as the last I with C I^{-1/4} log 2I ≥ 1/(2g). It depends only on g and is
  even effective.
* **Large index needs no exceptional set.** For I > I_* and W ≥ W_0(g), where W_0(g) controls the second term of (2),
  one gets Δ > H^{-1/g} for every tuple.
* **The middle range is finite.** Sublattices of Z^4 with d0 ≤ I ≤ I_* are finite in number (Hermite normal form).
  Each lattice has 4^4 unit translates ε, and §2.1 applies to each pair with g < B_{d0} ≤ B_I.
* **Z is proper closed.** Each Z_{L,ε} is the image of a proper closed set under a finite morphism of
  4-dimensional varieties, so it is closed of dimension ≤ 3. A finite union of these is a proper closed subset of
  the irreducible torus T.
* **H_0 depends only on (d0, g).** It is the maximum of finitely many fixed-cover constants and e^{W_0(g)}.
* **Standing cautions.** The many-lattice problem is handled by Baker, not by a Subspace or MR exceptional set; MR is
  applied with S = {∞} on finitely many fixed covers. So none of the cautions is violated. H_0 and Z are
  ineffective (they come from MR via Faltings–Wüstholz).
* **Numerics, exact.** B_5 > 17/5 because 350 > (17/4)^4 = 326.25. B_7 > 94/25 because 4.7^4 = 487.97 < 490.
  94/25 > 124/33 because 3102 > 3100, and 124/33 = 4/(1+2/31).
* **d = 7 is sharp for this method.** B_6 < 124/33 because 420 < (155/33)^4 ≈ 486.7.

### 8. §2.3, deleting a fixed set S_0 of primes: CORRECT
* **Example.** Columns e_1 at p_0, 7e_1 at p_1, and e_2, e_3, e_4 by singleton allocation at three more distinct
  split primes. The z_j are distinct and primitive. I computed I_full = 1 and I_out = 7.
* **Case A (bounded outside index).** Write b_p = r_p + ℓ_p with ℓ_p ∈ L_out and r_p from a fixed set of coset
  representatives. Then y ∈ t·f_P(T'(Q)) with t = ε ∏_{S_0} g_p^{r_p}, and there are finitely many t. The height
  is still W, and §2.1 applies with degree I_out. Correct, modulo §2.1.
* **Case B (unbounded outside index).** Verified exactly:
  * (y^c)^4 = u^{4m} ∏_{S_0} g_p^{4r_p}, with |r_p| ≤ m/2.
  * The power is not the identity, because its valuation at some outside prime is 4m<v,b_p> ≠ 0.
  * 2h(u) ≤ 4||v||_∞ W + Σ_{S_0} log p.
  * BW applies with |S_0| + 2 logarithms; the coefficients and the winding number are O_{S_0}(m). The constants
    depend on S_0, as stated.
* **Case C (outside rank r < 4).** Verified exactly:
  * The signed r-minors on r+1 rows, extended by zeros, annihilate every outside column. This holds because each
    outside column lies in the Q-span of the r chosen ones.
  * c ≠ 0, since c_{j0} is ± the chosen nonzero minor.
  * ||c||_∞ ≤ r!(W/log 5)^r.
  * (y^c)^4 = ∏_{S_0} g_p^{4<c,b_p>}, which is not the identity by overall full rank. Its coefficients are
    O((1+W)^4).
  * Conclusion: −log Δ ≪_{S_0} log(2+W).

### 9. Dependence outside the block
The theorem's middle range, and case A of §2.3, rest entirely on §2.1. That section needs two things:
* McKinnon–Roth Corollary 4.2 and Theorem 5.1 on the normalised cover X_f, with B_f^4 = 70d.
* The comparison dist(x, q) ≍_{t,f} Δ.

Against `papers_check/read_mr.md`, the cited statements have the right shape:
* Thm 5.1 with S = {∞}: if β R > 1, then there are finitely many solutions of d ≤ H^{-(R+δ)} off a proper Z.
* Cor 4.2: β ≥ (n/(n+1)) (Vol/mult)^{1/n}.

I did not referee §2.1 in detail.

## Novelty
* The record uses linear forms in logarithms only with one logarithm per prime in the support (Matveev,
  `round8/stability.md` Prop 6.1). Its conclusion there is about support size, and it is capped by the fakes.
* Astra's device is new to the record (grep: no "dual lattice" or "valuation lattice" index argument anywhere). A short
  vector of the dual valuation lattice collapses the whole product into one element u, whose height is
  ≪ I^{-1/4} W, so only two logarithms are needed whatever the number of primes.
* It is a standard technique in Diophantine approximation, but it is not a routine consequence of anything in the
  record.
* Inequality (2) is unconditional and effective. A notable by-product: every full-rank five-tuple with
  Δ ≤ H^{-κ} has valuation index I ≤ I(κ), effectively. In particular this holds on √R-arcs, where κ ≈ 1/4.

## Significance
* It does not change the growth rate, the uniform bound M(C), Question A, or route 5.
* For Question B (the A_5 form) it removes all full-rank tuples of index ≥ 7, off a proper Z. Combined with §2.1 it
  also removes those of large index for any exponent.
* The remaining case, index ≤ 6 and especially index 1, is the generic one. Singleton and squarefree profiles
  typically give I = 1, and I expect the fake-circle profiles of F1 to be there too (I did not check this). So the
  hard core of Question B is untouched.
* Verdict: correct, cleanly uniform, and modest. It narrows Question B's search space to small index but is
  peripheral to the main obstruction.

Sources: [Baker–Wüstholz, Crelle 442 (1993)](https://www.degruyterbrill.com/document/doi/10.1515/crll.1993.442.19/html);
[arXiv:1611.04664 (secondary statement of the two-log case)](https://arxiv.org/pdf/1611.04664).
