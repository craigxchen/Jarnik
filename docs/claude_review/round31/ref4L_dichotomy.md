# Referee report, round 31, Part IV, block "dichotomy" (astra_lemma_IV/answer.md, section IV.4)

Scripts and outputs are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref4L_dichotomy/`. I wrote all of them and did not run Astra's scripts.
* `sym.py` / `sym.out`: symbolic checks. They cover the interpolation identity for a generic binary quadratic and five generic vectors, the binary Vandermonde (= product of the three pair minors, sign +1), the Gram identity, the polarization formula, and the profile formulas n_i = D B_i/A_i^2, g_ij = D J_ij/(A_iA_j) and prod_{j≠i} g_ij = D^2B_i/A_i^3 in both orientations. They also count marked vertices per edge.
* `dich.py`: exact (Fraction/sympy Rational) evaluation of the whole IV.4 proof on any model. It computes m_ij, t_ij, the (IV.5) left and right sides and their denominators, Q_* from three columns and its check on all five, den(Q_*) | H_t, E, det E, the Gram identity for all ten pairs, the lower bound 1/(4H_t^2), and the rank. In the rank-one case it extracts λ and ℓ.
* `run_tests.py`/`.out`: the IV.2 family with Astra's quintuple, P from 10^30 to 10^120, several q. Also Astra's fully factored balanced parent (IV.3 certificate).
* `run_tests2.py`/`.out`: genuine Gaussian stars, stretched genuine stars (X, dY), Z[√−2]-stars and perturbed models V_i + n_i u_i, in both orientations.
* `asym.py`/`.out`, `Enum.py`/`.out`: E and det E along the IV.2 family.
* `lattice_c.py`/`.out`: the index of {Q : Q(V_i) ≡ 0 mod n_i}, computed by Smith normal form, and the number of directions mod each block prime.

## Verdicts in one place

| Item | Verdict |
|---|---|
| Interpolation identity Σ Q(V_i)V_i/∏_{j≠i} m_ij = 0 | CORRECT (symbolic, generic) |
| (IV.5), denominator of left side ∣ ∏_{a<b; a,b≠k} t_ab, vanishing for max δ ≤ P^{h−1−6κ−ε} | CORRECT |
| Σ c_iA_iV_i = 0; existence and uniqueness of rational Q_* with Q_*(V_i) = n_i | CORRECT |
| den(Q_*) ∣ H_t (local argument, both orientations) | CORRECT (needs blocks pairwise coprime, which "private" supplies) |
| Bound on E(V_i,V_j), Gram inequality, threshold P^{h−3−12κ−ε} | CORRECT |
| Rank-one alternative δ_i = λℓ(V_i)^2; mixed signs impossible | CORRECT |
| κ = 10η | CORRECT (with J_e = P^{1+o(1)} and relative norm error o(1)) |
| Sharpness against IV.2 (h = 4) | CORRECT, but only in the exponent at κ → 0. No sharpness claim for h = 5 (none is made) |
| Last remark at κ = 0 (one zero defect ⇒ all zero) | CORRECT and understated: it holds for every κ in the theorem's range |
| (c) correction: one marked line per prime, lattice index N | CORRECT (proof and SNF computation) |

I found no gap and no false statement in IV.4. Its scope is narrower than the word "dichotomy" suggests; see the last two sections.

## 1. Line-by-line check

**Setup.** h = 4 is the uncomplemented orientation: S_i is the 4 edges at i, |S_i ∩ S_j| = 1, n_i = P^{4+o(1)} and g_ij = J_ij. h = 5 complements the cycle F = {12,23,34,45,51}: |S_i| = 5 and |S_i ∩ S_j| = 2. Both are confirmed in `sym.out`. With pairwise coprime blocks, g_ij = gcd(n_i,n_j) = ∏_{S_i∩S_j} J_e, and n_in_j/g_ij^2 is the product of the six separating blocks.

**Interpolation identity.** The chart argument is right. With V_i = z_i(1,u_i), m_ij = z_iz_j(u_j−u_i), so ∏_{j≠i} m_ij = z_i^3 Z ∏_{j≠i}(u_j−u_i) with Z = ∏z. The term becomes q(u_i)(1,u_i)/(Z∏(u_j−u_i)), where q(u)(1,u) has degree ≤ 3. The fourth divided difference vanishes, and the sign (−1)^4 is harmless. Verified symbolically for both coordinates (`sym.py`). It needs only pairwise non-parallel V_i, i.e. t_ij ≠ 0. It is the same identity as Part III's (8), applied to quadratics times V instead of cubics; it is not new.

**From I to (IV.5).** n_i/∏_{j≠i}g_ij = A_i/D, which I checked symbolically. So the I-identity reads Σ c_iA_i(1+δ_i/n_i)V_i = 0. Take the determinant against V_k and use A_im_ik = t_ik D J_ik/A_k; cancelling D/A_k gives (IV.5) exactly. `dich.py` checks LHS_k = RHS_k as exact rationals on every test model. The left side is Σ_{i≠k} ±J_ik/∏_{j∉{i,k}} t_ij. Each denominator is a product of three of the six t_ab with a,b ≠ k, so the left side lies in (∏_{a,b≠k}|t_ab|)^{−1}Z. This was checked (`lhs_den_ok`) on all models. On the right, |c_it_ik| = 1/|∏ of three nonzero integers| ≤ 1 and J_ik/n_i ≤ P^{1−h+o(1)}. So vanishing holds for max δ ≤ P^{h−1−6κ−ε}, which is n^{3/4} for h = 4 and n^{4/5} for h = 5. The aside that a defect exponent < 4/5 forces the norm-only (8) on balanced rows is correct. III.11 sits exactly at n^{4/5}, consistent with that.

Perturbed models V_i + n_iu_i keep all ten divisibilities but have defects ≈ n_i^{3/2}. Their (IV.5) left side is nonzero and no Q_* exists (`run_tests2.out`). So the vanishing is a real constraint and not automatic.

**Q_*.** From LHS_k = 0 for all k, det(Σc_iA_iV_i, V_k) = 0 for all k, so Σ c_iA_iV_i = 0 because the V_k span. The map y ↦ Σ y_iV_i/∏_{j≠i}m_ij from Q^5 to Q^2 is onto: one y_i alone gives a multiple of V_i. Its kernel therefore has dimension 3 and contains the 3-dimensional image of quadratic evaluation (≥ 3 distinct directions), so the two are equal. Since y = n satisfies the relations (n_i/∏m_ij = c_iA_i/D), Q_* exists and is unique and rational. Checked: Q_* solved on three columns matches all five n_i on every model with small defect.

**Denominator bound.** The 3×3 evaluation determinant is +m_ab m_ac m_bc (symbolic), and the Cramer numerators are integers, so v_p(den) ≤ v_p(m_ab m_bc m_ca) for every triple. Fix a prime p of a private block J_e. The vertices marking e number 2 (h = 4, and h = 5 off the cycle) or 3 (cycle edges for h = 5); `sym.out` prints "marked per edge [2]" and "[2, 3]". So there are always ≥ 2 unmarked vertices. A triple of two unmarked vertices plus one more has no pair that both mark e. Its g's are therefore units at p, and v_p(den) ≤ v_p(t_ab t_bc t_ca) ≤ v_p(H_t). At other primes all g's are units. Hence den(Q_*) | H_t, confirmed numerically (`den_divides_Ht` True on all valid models). Two hypotheses are used here: the blocks must be pairwise coprime (if p divided two blocks, the triple choice could fail), and the V_i must be integral. Primality is not used.

**E(V_i,V_j) and Gram.** Write V_k = aV_i + bV_j with a = m_kj/m_ij and b = m_ik/m_ij; then 2E(V_i,V_j) = (δ_k − a^2δ_i − b^2δ_j)/(ab) (symbolic check). All g_ab = P^{h−3+o(1)} because J_e = P^{1+o(1)}, and this hypothesis is needed. So |a|^{±1}, |b|^{±1} and |a/b|^{±1} are ≤ P^{κ+o(1)}, and |E(V_i,V_j)| ≤ P^{2κ+o(1)}max|δ|. The Gram identity δ_iδ_j − E(V_i,V_j)^2 = m_ij^2 det E is exact (symbolic, and checked on all ten pairs of every test model). With m_ij^2 ≥ g_ij^2 ≥ P^{2(h−3)−o(1)} this gives |det E| ≤ P^{−2(h−3)+4κ+o(1)}max δ^2. The coefficients of E lie in H_t^{−1}Z and det E = e_1e_3 − e_2^2/4, so 4H_t^2 det E ∈ Z and a nonzero |det E| is ≥ 1/(4H_t^2) ≥ P^{−20κ−o(1)}/4. Combining, det E ≠ 0 forces max δ ≥ P^{h−3−12κ−o(1)}. The exponent bookkeeping is right.

**Rank-one alternative.** det E = 0 and E ≠ 0 gives a rank-one rational symmetric matrix, which is μuu^T with u rational. Rescaling gives λℓ^2 with ℓ primitive integral. Then δ_i = E(V_i) = λℓ(V_i)^2, and all nonzero defects share the sign of λ. Correct.

**Zero defects (Astra's last remark).** If δ_i = 0 then ℓ(V_i) = 0, and primitivity of V_i gives ℓ = ±[V_i,·]. Next, |λ| ≥ 1/H_t: take an integer w with ℓ(w) = 1; then λ = E(w) ∈ H_t^{−1}Z. So for j ≠ i, |δ_j| = |λ|m_ij^2 ≥ g_ij^2|t_ij|^2/H_t ≥ P^{2(h−3)−10κ−o(1)}. Astra states this only at κ = 0. In fact 2(h−3)−10κ > h−3−12κ for every κ ≥ 0, so **in the whole range of the theorem a single zero defect forces all defects to be zero**. This needs V_i primitive; without primitivity, divide by the content, at a cost of P^{o(1)} if the content is small.

**κ = 10η.** |t_ij| = |m_ij|/g_ij ≤ (1+o(1))Δ·sqrt(n_in_j)/g_ij = Δ·P^{3+o(1)} = P^{10η+o(1)}, since N = P^{10+o(1)}. So the thresholds are defect exponents 1/4 − 30η (h = 4) and 2/5 − 24η (h = 5) in terms of n_i. Correct.

## 2. Numerical tests (exact)

Run on every model:
* **IV.2 family (h = 4)**, Astra's quintuple, with P ~ 10^30, 10^60, 10^120 prime, q = 30, 210, 4290, and blocks checked pairwise coprime.
  * The columns are primitive, the defects have mixed signs, and (IV.5) vanishes. Q_* exists, den(Q_*) | H_t, and rank E = 2.
  * **New observation:** det E = −16908544/59049 = −(4112/243)^2, **exactly the same for every P and q tested** (`Enum.out`).
  * So E is a split indefinite form, a product of two rational linear forms, with huge coefficients (the Y^2 coefficient grows like P^2, the X^2 coefficient like P, at fixed q; `Enum.out`) and bounded determinant.
  * This is precisely the sharp case of the Gram inequality: |det E| ≍ P^{−2}·(Pq^3)^2/q^6 = O(1). Mixed signs come from the V_i lying on both sides of a rational null line of E. In the language of the dichotomy, IV.2 is "split rank two at bounded determinant", one step from the rank-one degeneration.
  * At computable sizes the cofactors t_ij = LB(x_i−x_j) ≈ 10^{14}q^3 make κ large (κ ≈ 0.23–0.73), so the theorem is vacuous there. It is purely asymptotic: q fixed, P → ∞.
* **IV.3 certificate (h = 5)**: (IV.5) vanishes and Q_* exists with den(Q_*) = 13, which divides H_t. Rank E = 2, all δ > 0, |det E| ≈ P^{1.27}, so this model is far from the threshold, as expected.
* **Genuine Gaussian stars** (ten split primes ≡ 1 mod 8 near 10^6, both orientations): E = 0, Q_* = I.
* **Stretched genuine stars** (X, dY), d = 7 or 3: rank one with λ = 1 − 1/d^2 and ℓ = Y; den(Q_*) = d^2. **Z[√−2]-stars** (products of primes of x^2+2y^2): rank one with λ = −1 and ℓ = Y. Both satisfy all ten divisibilities, with private blocks and primitive columns. **So the rank-one branch is not vacuous algebraically.**
  * At small defect it requires all five columns in a strip |ℓ(V)| ≤ (H_t max δ)^{1/2} ≤ P^{(h−3)/2−κ} around a rational line.
  * It also requires them to be an *exact* star model for the rational form Q_* = I − λℓ^2. The height of Q_* is not bounded, and its denominator is up to H_t.
  * I found no construction with small defect, nor any obstruction.
* **Perturbed models** V_i + n_iu_i: (IV.5) fails and no Q_* exists. This is consistent, since the defects are far above P^{h−1}.

## 3. The correction of suggestion (c)

Correct.
* At p^a ∥ J_e, the marked columns satisfy p^a | g_ij | m_ij. Primitive vectors with determinant ≡ 0 mod p^a are proportional by a unit, so Q(V_j) ≡ c^2Q(V_i) and the condition is one scalar congruence Q(L_e) ≡ 0 mod p^a. It is surjective because L_e has a unit coordinate. CRT then gives index ∏J_e = N.
* SNF computation on small genuine and perturbed private profiles in both orientations: the index equals N in all 29 tests (`lattice_c.out`).
* The unmarked columns are not parallel mod p. In genuine stars too I count 4 directions mod p (1 marked + 3 unmarked) for h = 4, and 3 or 4 for h = 5.
* So the "two points L_p^±" premise fails even for genuine stars in the star normalization. It is valid only in the original common-norm normalization ζ_i, where both sides of a cut are isotropic.
* For III.11, Astra's claim of distinct directions rests on the blocks being coprime to the cofactors u_i−u_j. That coprimality was checked in `ref3_model.md` and I did not recheck it.
* The coordinator's covolume N^2 (the class λQ + N·forms) is the wrong object. The true lattice has covolume N and contains the Euclidean form only up to the defect, as Astra says.

## 4. Clean statement

**Theorem (IV.4, as verified).** Let h ∈ {4,5} be the orientation (h = 4: each pair cut marks its two endpoints; h = 5: the five cycle edges are complemented). Let J_e (e ∈ E(K_5)) be pairwise coprime positive integers with J_e = P^{1+o(1)}. Put n_i = ∏_{e∈S_i}J_e and g_ij = ∏_{e∈S_i∩S_j}J_e. Let V_1,…,V_5 ∈ Z^2 be primitive with g_ij | m_ij = det(V_i,V_j), t_ij = m_ij/g_ij ≠ 0, max|t_ij| ≤ P^{κ+o(1)}, and δ_i = |V_i|^2 − n_i. Fix ε > 0. As P → ∞:
1. if max|δ_i| ≤ P^{h−1−6κ−ε}, then Σ_i c_iA_iV_i = 0, and there is a unique Q_* ∈ Q[X,Y]_2 with Q_*(V_i) = n_i and H_tQ_* integral (H_t = ∏|t_ab|);
2. if max|δ_i| ≤ P^{h−3−12κ−ε}, then either all δ_i = 0, or E = I − Q_* = λℓ^2 with λ ∈ H_t^{−1}Z∖{0} and ℓ primitive integral. In the latter case δ_i = λℓ(V_i)^2, and every δ_i is nonzero of the sign of λ.
3. If the angular width is Δ ≤ N^{−3/10+η}, one may take κ = 10η.

Sharpness: for h = 4 and κ → 0 the exponent h−3 = 1 cannot be raised, by the mixed-sign IV.2 family. The dependence on κ and the h = 5 exponent 2 are not known to be sharp. The best known balanced model is at P^{3−o(1)}, and its own pencil is blocked at P^{5/2}.

## 5. What it means for Question B (blunt)

* **Genuine tuples satisfy it trivially.** δ = 0 gives E = 0, the "exact" branch. **At defect zero the theorem carries no information**: there is no width bound, no bound on κ, and nothing about Z. It is not an inequality of type (b) that genuine clusters satisfy and fakes violate.
* **No known integral model violates it.** IV.2 (P^{1+o(1)}), IV.3 (o(P^3)) and III.11 (P^4) all lie above the threshold, where the theorem says nothing. So it constrains only hypothetical fakes below n^{1/4} (h = 4) or n^{2/5} (h = 5). Those must be rank one: exact star models of the nearby rational form I − λℓ^2, squeezed into a thin strip around a rational line.
* That reduces "small-defect fake" to an *exact* cluster problem for a model-dependent form of unbounded height and denominator, which is no easier than Question B itself. The standing caution applies: a fixed-ellipse result would not be uniform.
* So it neither proves nor disproves anything about e(0). The useful content is structural: the quarter-power barrier of IV.2 cannot be pushed lower by mixed-sign constructions. Any method tolerating defects < n^{1/4−30η} needs only to handle rank-one defects and exact norms.
* The observation that the IV.2 E has constant split determinant −(4112/243)^2 suggests the next test: can the two null lines of E be made to coincide (rank one) within a pencil construction?
* Novelty against the record: the identity is Part III's (8), and the Gram step reuses the II.9/III det Q = 1 device. The rank-one dichotomy and the index-N correction are new to the record. Significance for the goal is small: it does not touch defect zero.
