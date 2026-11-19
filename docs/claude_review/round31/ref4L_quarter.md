# Referee report, round 31, Part IV, block "quarter" (astra_lemma_IV/answer.md, sections IV.1 and IV.2)

Scripts are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref4L_quarter/`. I wrote all of them myself.
* `sym.py`: symbolic checks over Q[x_1..x_5,P,q] with generic marks. It checks the determinant identity, the two shell identities (with their exact s_4 corrections) and the residual formula (IV.1).
* `shapes.py`: the shape-density and dominance claims, i.e. the y = 1/x dictionary, the conic parametrisation, Jacobian ranks and the translation quartic.
* `search2.py`: a sieve plus gmpy2 search for P such that all ten P + q b_ij are primes ≡ 1 mod 4.
* `verify.py`: exact big-integer verification of a concrete instance. Output is in `verify.out`.
* `astra_copy/`: copies of Astra's two scripts, run with `python3 -I`. Both pass. Astra's IV.2 certificate is symbolic only; it contains no prime instance of the first family.

All processes ran under `prlimit --as=3.5e9`.

## 0. Bottom line

| claim | verdict |
|---|---|
| IV.1 lattice / index one | **CORRECT** |
| IV.1 prime specialisation (Tao–Ziegler) | **CORRECT WITH CORRECTIONS**: q is at most (log P)^L, not just P^{o(1)}; P need not be prime; q < 0 is allowed and harmless |
| IV.2 construction, identities, primitivity, shell identities, (IV.1) | **CORRECT** (reproduced symbolically and on two explicit prime instances) |
| IV.2 defect N^{1/10+o(1)} vs n_i ~ N^{2/5}; width N^{-3/10+o(1)}; "e(θ)=10/3 for θ>1/4" | **CORRECT WITH CORRECTIONS**: the exponent is reached only as a limit (log-power cofactors); the constants are astronomically large; the result is tied to the endpoint orientation |
| IV.2 exceptional-set quantifier | **CORRECT** (for each fixed Z, with a shape chosen after Z). What it does not say is important; see §4 |

**Novelty.** Real but modest. For integral five-mark fakes, it lowers the record's absolute norm defect from N^{2/5} (III.11, n^{4/5} at n ~ √N) to N^{1/10+o(1)}. The relative defect drops from N^{-1/10} to N^{-3/10+o(1)}, the same order as the angular width. It also adds the escape-from-every-Z statement for integral fakes.

**Significance.** Peripheral for the goal and for Question B. The models have no isotropy at any block prime and lack the dot-product divisibility. Their actual points of T_5 have pair heights N^{4/5}, against the N^{3/5} of the profile. So they are not near-counterexamples to anything on T_5. They are a barrier only to arguments that use the integer star matrix, the determinant divisibilities and archimedean norms to precision N^{1/10+ε}, and nothing about |V_i|^2 modulo block primes. That class of arguments is narrow: the quarter-power threshold in this orientation is exactly the size of a single block prime, so a local isotropy congruence at one block prime would force δ_i = 0 there.

## 1. IV.1: lattice and index one

With endpoint marking, the valuation column of block {i,j} is e_i + e_j.
* The e_i + e_j generate D_5 = {even coordinate sum}: they give e_j − e_k, and then 2e_i.
* D_5 has index 2 in Z^5, and **1** has odd sum. So D_5 + Z**1** = Z^5, and the image in the anchored lattice Z^5/Z**1** is everything.
* Reversing a side replaces the column by **1** − (e_i + e_j), which is ≡ −(e_i + e_j). So the index is unaffected.

Exact check: the Smith form of the 4×10 anchored matrix is (1,1,1,1) (`verify.py`).

The Liouville/width formula sin(angle_ij) = (1+o(1)) |t_ij| / sqrt(n_i n_j / g_ij^2) with g_ij = J_ij and n_i n_j / g_ij^2 = (product of the six separating blocks) is checked exactly: the ratio equals 1 to 1e-13 and 1e-19 in the instances. **CORRECT.**

## 2. IV.1: prime specialisation

The cited theorem exists: T. Tao and T. Ziegler, *Narrow progressions in the primes*, in *Analytic Number Theory (in honour of H. Maier)*, Springer 2015, pp. 357–379; arXiv 1409.1327. The statement as given in the abstract (confirmed by web search; the arXiv mirrors were not reachable from this session, so I could not read the full text) is:

> For polynomials P_1, …, P_k ∈ Z[m] with P_1(0) = … = P_k(0), and A ⊂ primes in [N] of relative density ≥ δ, there is a + P_i(r) ∈ A for all i with 0 < |r| ≤ log^L N, for N large in terms of k, the P_i and δ. In the linear case L does not depend on δ. Their earlier paper (*Polynomial patterns in the primes*) had |r| ≤ N^{o(1)}.

Application:
* Take P_b(r) = m·b·r for the fixed shifts b (including b = 0) and m equal to the required divisibility modulus.
* Take A = primes ≡ 1 (mod 4) in [N/2, N], which has relative density about 1/4.
* Then P = a, q = m r and |q| ≤ m log^L N, so q ≤ (log P)^{O(1)}.

Astra's variant, which embeds into a long arithmetic progression and spaces the indices by m, is equally valid. Corrections:
1. **Sign of q.** The theorem gives 0 < |r|, so q may be negative. The construction is polynomial in (P, q) and the primitivity argument only uses p | q, so a negative q is harmless.
2. **P need not be prime.** IV.2 only needs P ≡ 1 (mod 4) and P coprime to the primes of q. That is what primitivity uses, since the first entry of V_i is ≡ P^2 mod p.
3. **The bound on q is log-power.** The diagonalisation remark is unnecessary.

**Bounded q is not available unconditionally.** Fixed q would need the ten-tuple prime conjecture (Dickson / Hardy–Littlewood).
* An unconditional fixed-q route with composite squarefree blocks might exist through results on linear forms in sums of two squares (Matthiesen). That route would also need pairwise coprimality and split-only prime factors.
* I have not verified this; it is a suggestion, not a claim.

## 3. IV.2: the construction, checked line by line and on primes

**Symbolic checks (`sym.py`, generic marks).**
* The identity x_i h_j − x_j h_i = (x_i − x_j) b_ij holds identically. Hence [V_i, V_j] = LB(x_i − x_j)(P + q b_ij) exactly, because the V_i are the images of (x_i, U_i) under [[B, P], [0, L]].
* First shell identity: Σ_{j≠i} b_ij = 2h_i, identically.
* Second shell identity: e_2({b_ij}) = h_i^2 + 3s_5 x_i − s_1 s_5 − s_4 x_i^2 + s_1 s_4 x_i, identically. This is exactly Astra's "before setting s_4 = 0" form.
* With l^2 = −s_1 s_5 and b_0 = 3s_5/2, the P^2 q^2 coefficient of |V_i|^2 − n_i is −x_i (s_1 − x_i) s_4. It vanishes iff s_4 = 0.
* The P q^3 coefficient minus s_5 x_i (s_2 x_i − s_1 s_2 − s_3) is a multiple of s_4. So (IV.1) holds on s_4 = 0.
* `verify.py` also checks that the residual is **exactly** (lead)·P q^3 + c_4 q^4 with c_4 = b_0^2 x_i^2 + l^2 h_i^2 − e_4({b_ij}).

**The numerical shape.** For x = (−27, −270, 135, 120, 40):
* s = (−2, −53925, −5042250, 0, 4723920000) and −s_1 s_5 = 97200^2.
* b_0 = 7085880000.
* The ten b_ij are distinct, nonzero and all divisible by 2·3·5.
* The leading residual coefficients are −8.15e17, −2.49e19, −1.50e18, −8.71e17 and +5.25e17 (times q^3). These have both signs, as claimed.

**Explicit prime instances (new; exact big-integer checks).** Take q = 30. Both of the following pass every check:
* P = 1642245131041
* P = 112908008110501

For both values:
* All ten J_ij = P + 30 b_ij are prime (gmpy2, 50 rounds), ≡ 1 (mod 4), and distinct.
* The five columns are primitive.
* Each g_ij = gcd(n_i, n_j) equals J_ij.
* All ten minors equal LB(x_i − x_j) J_ij exactly. The cofactors are l b_0 q^3 (x_i − x_j), with l b_0 q^3 = 1.86e19 and max |t| = 7.53e21.
* The index is one.
* diam{arg V_i} · N^{3/10} = 7.5315e21, matching the prediction l b_0 q^3 · max|x_i − x_j| = 7.5315e21.
* The normalised vector (ξ_i/ξ_5 − 1)_i equals [x_i − x_5]/(x_1 − x_5) = (1, 310/67, −95/67, −80/67) to 4–5 digits.

**Caveat on the constants.** At these sizes log|δ_i| / log n_i is 0.70–0.74 (first instance) and 0.64–0.67 (second), and log diam / log N is −0.12 and −0.14. Both are far from the limits 1/4 and −3/10.
* The defect constant is s_5 x_i (…) q^3, about 10^{24} for q = 30. So |δ_i| ≤ n_i^{1/2} needs P ≳ 10^{24}.
* The width constant is about 7.5e21.
* No prime instance in the asymptotic regime is computable this way: ten simultaneous ~25-digit primes in this pattern are out of reach. The asymptotic statements are proved; the instances certify the exact algebra, not the exponents.

**The resulting statements, corrected.**
* With q ≤ (log P)^{O(1)}:
  * |δ_i| ≪_x P q^3 + q^4 = N^{1/10} (log N)^{O(1)};
  * n_i ~ P^4 = N^{2/5};
  * every δ_i ≠ 0 for large P, since at most one mark can have a vanishing leading coefficient and here none does.
* Hence |δ_i| ≤ n_i^θ eventually, for every fixed θ > 1/4. This fails at θ = 1/4 itself, because of the constant times q^3.
* The width is Δ ≍_x q^3 P^{-3} = N^{-3/10} (log N)^{O(1)}.
* The cofactors l b_0 q^3 (x_i − x_j) grow like (log N)^{O(1)}: they are N^{o(1)} but **not bounded**.
* So the exponent 10/3 is attained **only as a limit**: log(1/Δ)/log N → 3/10. It is not attained in the form Δ ≍ N^{-3/10}. That form holds conditionally on the ten-tuple prime conjecture with q = 30 fixed.
* This does not matter for e(θ), which the coordinator defined with slowly growing cofactors and which is a limiting exponent. It is irrelevant to Question B, where the target 15/4 exceeds 10/3 by a fixed amount.
* It **would** matter if e(θ) were redefined with bounded cofactors. Then row 1 is conditional, unless the composite-block route of §2 can be made to work.
* Liouville gives e ≥ 10/3 for these profiles whenever δ = o(n), so "e(θ) = 10/3 for θ > 1/4" is correct in this limiting sense.

**Orientation.** The defect bound n_i^θ is absolute and n_i ~ N^{2/5} here, so the statement belongs to the endpoint orientation. In terms of N:
* first family: defect N^{1/10+o(1)};
* III.11: defect N^{2/5};
* IV.3: o(N^{3/10}).

Astra states this dependence openly. The coordinator's e(θ) was not orientation-invariant, and the record should state e(θ) with the orientation attached.

## 4. The exceptional-set quantifier

**What was checked (`shapes.py`).**
* x = 1/y turns s_4 = 0 into Σ y_i = 0, and −s_1 s_5 = −e_4(y)/e_5(y)^2.
* With y_5 = −s − d, the identity −e_4(y) = g d(d+s) + s a b c holds.
* The substitution a = bc(b+c)/(k^2 − bc) gives s a b c = (ka)^2.
* The line through (d, v) = (0, ka) with slope m gives d = (g s − 2 k a m)/(m^2 − g), and −e_4 is a square there.
* The resulting rational map (b, c, k, m) → (a, b, c, d) has Jacobian rank 4. Its composite with the shape map to [x_i − x_5] ∈ P^3 has Jacobian rank 3 at a random rational point.
* So the shapes coming from Q-points are Zariski dense in P^3. This is a slightly shorter route to Astra's "density plus dominance", and it confirms it.
* The translation quartic e_4(z + t) has leading coefficient 5, as claimed.
* Scaling x by an even λ keeps the square condition (it scales by λ^6) and makes b_0 integral.
* The bad loci (zeros, coincidences, b_ij collisions or zeros) are proper closed sets, because the displayed example avoids them.

**Leading-term argument.**
* Put u_i = ξ_i/ξ_5 − 1. Then u = ε(w_0 + O(q/P) + O(ε)), with ε = 2i l b_0 q^3 P^{-3} and w_0 = (x_i − x_5)_i real.
* The common scalar is complex, but the lowest Taylor form F_d of any F vanishing on Z is homogeneous. So F(1 + u) = ε^d (F_d(w_0) + o(1)), uniformly, because q/P → 0.
* Choose the shape off {F_d = 0}. Then F does not vanish at the model points for large P.

The argument is **CORRECT**, with the quantifier order Z, then the shape, then P. No single family escapes every Z. That is enough to defeat statements of the form "∃ Z ∀ tuples", which is the form of A_5 and (L_5°).

**What it does not say.**
1. **These are not points of T_5 close to e at exponent 10/3 relative to their own height.**
   * Exactly, in both instances, the rational content of conj(V_5) V_i is 1.
   * So the height of ξ_i/ξ_5 in the norm-one torus is |V_i|^2 |V_5|^2 ≈ N^{4/5}, not the profile's n_i n_5 / g^2 = N^{3/5}. A genuine point with this profile has content g_{i5}.
   * The missing factor is the Gaussian gcd, which in a genuine tuple comes from J | V_i · V_j. In the models, J_ij divides the determinant but never the dot product (0 of 10 pairs). No marked block is isotropic (0 of 20 incidences with J_e | |V_i|^2).
   * Measured against their true toric height, the model points are nowhere near critical. Exponent 10/3 refers to the profile height N, which these points do not have.
2. **They say nothing about genuine clusters**, or about whether genuine tuples off some Z have exponent > 15/4.
3. **The barrier covers a narrow class of arguments.**
   * It defeats arguments that take as input the primitive integer star matrix, all ten determinant divisibilities, index one, and |V_i|^2 = n_i + O(n_i^{1/4+ε}), and nothing else, even after a fixed exceptional set is discarded.
   * In this orientation n_i^{1/4} ≍ P is the size of one block prime. Any argument that uses the norm congruence |V_i|^2 ≡ 0 (mod J_e) at even one marked block excludes these models, because there |δ_i| would have to be a nonzero multiple of a prime of size P.
   * Equivalently, the models fail the coordinator's isotropy/quadratic-form condition (c) at every prime.
   * So the result sharpens "integrality plus divisibilities is not enough" (record item 12) to "…plus archimedean norms to N^{1/10+ε}". It does not reach "…plus local norm information".

## 5. Novelty and significance

* **Novelty.** The record has the III.11 model, with absolute defect N^{2/5} and relative N^{-1/10}. It also has the C1 / form-blind relaxations and the gap-only angle model at 10/3.
  * IV.2 is a genuinely new explicit family. Its Newton-identity matching of the P^3 q and P^2 q^2 coefficients is not in the record (grep finds no "norm defect" outside round 31).
  * It lowers the defect by a factor N^{3/10} and adds the escape statement.
  * Like III.11, it is a fixed rational projectivity of five points [x_i : P + q h_i] of P^1. It is still form-blind.
* **Significance for the goal.** None. The growth rate, the uniform bound and Question B are unchanged.
* **Significance for Question B.** It rules out one weak proof strategy (archimedean norm control plus determinant divisibility) and locates the remaining information in local isotropy, i.e. dot-product divisibility. That is the Gaussian-gcd content that genuine tuples have and these models lack.
  * This is a useful negative calibration.
  * It is not "a new barrier" in the sense the coordinator hoped for in (a). That would need θ → 0, and Astra has no construction below θ = 1/4 in this orientation.
