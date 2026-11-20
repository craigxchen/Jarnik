# Referee report, block ATTACK on Part VIII: can Astra's consequences (A), (B), (C) refute (L_8)?

Referee scripts are in `scratchpad/round31/ref8L_attack/`. All of them are my own code, and none imports Astra's script. Arithmetic is exact (Python integers and `Fraction`, plus C with 64/128-bit integers for the exhaustive scans; floating-point logarithms are used only to rank results).

## 0. Verdict

**NOT REFUTED.** No family checked here breaks (A), (B) or (C). Every construction I found, and every construction in the record, fits the predictions of (L_8) and of the published Cilleruelo–Granville Conjecture 10.

The main findings:

1. **(A) is a special case of a published conjecture (Proved, trivial).** Cilleruelo–Granville Conjecture 10 with d = 2 implies (A) for every eps < 1/2. So refuting (L_8) through (A) would first require refuting CG Conjecture 10 (d = 2), which CG say they "strongly believe" even at exponent 1/d.
2. **Correct test for (A).** Let R(B) be the supremum of log q / log L over B-point two-class configurations. (L_8) implies (A)(eps), which implies R(B) < 2/(1+2eps) for all B > B_0(eps). So (L_8) predicts limsup R(B) < 2. Under the record's Vojta numbers it predicts R(B) < 2 - 4·10^-16 for large B.
   * The coordinator's phrase "approach 2 from below at a rate like 2 - c/B (consistent with (L_8))" is wrong. Such a rate would **refute** (A), and hence (L_8).
3. **A dictionary that produces all the good constructions (Proved).** Two-class configurations from polynomials correspond to *weighted Prouhet–Tarry–Escott* patterns (§2, Prop. 2.1).
   * Nodes n_1..n_k, weights v_j and r moment conditions give B integer points with log q / log L → k/(k-1-r).
   * The parallelogram, the cube, and the new families below are all instances.
4. **Best exponents found (Proved families, exactly verified).**

   | B | best exponent found |
   |---|---|
   | 3 | 3 |
   | 4–6 | 2 |
   | 7–8 | 7/4 |
   | 9–10 | 5/3 |
   | 11–20 | 3/2 |
   | B ≤ C(k, ⌊k/2⌋) | k/(k-2) |

   So R(B) ≥ 1 + (2+o(1))/log_2 B. The B = 5 and B = 6 families at exponent exactly 2 are new to the record.
5. **The best exponent decreases.** It falls below 2 at B = 7 and tends to 1, roughly like 1 + 2/log_2 B. It does not approach 2. This is what CG Conjecture 10 predicts, and it is far inside what (L_8) allows.
6. **At exponent exactly 2, polynomial families have B ≤ 4d (Proved, Plotkin).** Searches found B ≤ 6 for d ≤ 5 (Computed).
   * In **characteristic 2** the bound is attained and unbounded. Over F_{2^m}[T], the indicators of the affine F_2-hyperplanes of F_{2^m} give 2^{m+1} points at exponent exactly 2. This is the Sylvester–Hadamard code (Proved; checked for m = 3, 4).
   * This explains why the record's positive-characteristic clusters are unbounded, and why Hadamard/Paley "fakes" have no counterpart over Q.
7. **(B): the best near-miss is 16 divisors (Proved).** The integer n = ∏_{i=1}^4 (1+iq), with 4! | q, has 16 divisors ≡ 1 mod q, and log q / log n → 1/4 from below. Nothing better was found, and (B) allows any fixed count.
8. **(C): no family in the record meets the hypotheses of VIII.1 with an unbounded number of points (§4).**
   * The record's families have at most 10 marks.
   * The (L_5)/(L_6) refuters are not two-class, and when enlarged their Gram entries tend to 1.
   * The fakes have no integer realisation.
   * The characteristic-2 family is not over Z.
9. **Side result for B = 4 (Proved).** For four integers, q ≤ √2·L^{5/2}. This improves CG Theorem 9's exponent 3 for B = 4.

## 1. Setting, and the threshold re-derived

**Definitions.** For a finite set y of B ≥ 3 integers let L = max y - min y. Let q_2(y) be the largest q such that, at every exact prime power p^f ∥ q, the elements of y occupy at most two classes mod p^f.

**Lemma 1.1 (Proved; checked against the brute-force definition on 300 random tuples, `tc.py`).**
* q_2(y) = lcm over all bipartitions σ (including the trivial one) of g_σ, where g_σ is the gcd of the differences inside each side.
* If the gcd of all differences is 1, the g_σ are pairwise coprime and q_2 = ∏_σ g_σ.

*Proof.* At a prime p, the largest admissible f is max_σ v_p(g_σ). Suppose p divided g_σ and g_τ with σ ≠ τ. The join of two distinct 2-class partitions is the 1-class partition, so p would divide every difference. □

**The threshold.** Run Astra's VIII.2 matrices on a reduced configuration (W = 2 log Q, D_ij = 2 log(Q/h_ij)). This gives log|t_ij| = log|y_j - y_i| - log h_ij ≤ D_ij/2 - W/4 + λW with
λ = (log L / log Q - 1/2)/2, and the Gram entries satisfy G_ij ≤ 2 log L / log Q - 1 = 4λ.

So the hypothesis of VIII.1 at level λ is "ratio ρ := log Q / log L ≥ 2/(1+4λ)". With Astra's λ ≤ λ_* = eps/2 this is ρ ≥ 2/(1+2eps), that is, an interval of length Q^{1/2+eps}. Under the record's Vojta numbers (eps = 2λ_* = 2·10^-16), the threshold is ρ ≥ 2 - 4·10^-16 + O(10^-31). Astra rounds this to q^{1/2+10^-16}, which is safe.

**What a refutation of (L_8) through (A) needs.** For every eps > 0, configurations with unboundedly many points B and ρ ≥ 2/(1+2eps). Equivalently, limsup_{B→∞} R(B) ≥ 2, where R(B) = sup_y ρ(y) over B-point sets.

R(B) is non-increasing in B (Proved): deleting points can only enlarge q_2 and shrink L.

**Relation to Cilleruelo–Granville (Proved).** I read their Conjecture 10, Theorem 9 and Theorem 8 in `ref6L_pell/dl/cg07.txt`, §4.
* Conjecture 10, with at most d residues per prime power, says |Ω(b) ∩ [0, b^{1-ε})| ≤ N(d, ε).
* Shifting every residue set by -x_0 keeps the hypothesis, so CG Conj 10 (d = 2, ε = 1/2 - eps) implies (A)(eps) for every 0 < eps < 1/2.
* CG also write that they "strongly believe" the exponent 1/d version. For d = 2 that is boundedness at ρ = 2 exactly, which is of the same strength as (A).
* Attribution note: in this extraction, "Conjecture 9 ⇒ Conjecture 1 (Rudin)" is **Theorem 11**. Theorem 12 there is "Conjecture 12 ⇒ Conjecture 1". Astra cites "Conjecture 9 and Theorem 12" for the square-root route, perhaps from a different version of the paper. This does not affect the logic.

## 2. Constructions for (A)

### 2.1 The node dictionary

**Proposition 2.1 (Proved).** Let n_1..n_k be distinct integers, v_j nonzero integers and 0 ≤ r ≤ k-2. Let
𝓔 = {ε ∈ {0,1}^k : Σ_j ε_j v_j n_j^p = 0 for p = 0..r-1}.
For an integer T, put Y_ε(T) = Σ_j ε_j v_j ∏_{l≠j}(T - n_l). Then for all large T:

1. Y_ε(T) ≡ ε_j v_j ∏_{l≠j}(n_j - n_l) mod (T - n_j). So the |𝓔| integers occupy two classes modulo each T - n_j.
2. Y_ε is a polynomial of degree ≤ k-1-r. Indeed Σ_j a_j ∏_{l≠j}(T-n_l) = w(T) Σ_{p≥0} μ_p T^{-p-1}, with w = ∏(T-n_l) and μ_p = Σ_j a_j n_j^p, so the top r coefficients vanish.
3. The Y_ε are distinct, because their values at the nodes differ.
4. q_2 ≥ lcm_j |T-n_j| ≥ ∏_j|T-n_j| / ∏_{i<j}|n_i-n_j|, and L ≤ C·T^{k-1-r}.

Hence log q_2 / log L → k/(k-1-r) as T → ∞, with B = |𝓔| points. The homogeneous version (forms in (α, β), nodes in P^1(Q) including ∞) is identical.

**Converse heuristic.** Conversely, any polynomial family with simple contacts has this form. Its patterns form a binary code inside an evaluation (generalised Reed–Solomon) code, and 𝓔 is the set of 0/1 vectors in its dual description. This direction is not needed for any statement below.

**Corollary 2.2 (Proved, Plotkin).** Two distinct points agree at at most k-1-r nodes. So the patterns form a binary code of length k and minimum distance at least r+1.
* At exponent k/(k-1-r) = 2 (k = 2r+2), this gives B ≤ 2k (Plotkin bound A(n, n/2) ≤ 2n).
* At any exponent > 2, B is bounded by Plotkin. This is the code form of CG Theorem 9.

### 2.2 Families, all verified exactly

The families were checked in `families.py`, `verify_nodes.py` and `polyfam.py`. Each check confirms two classes modulo each node value and computes q_2 or lcm exactly, at T = 10^6 + 7 and T = 10^30 + 7 (families.py: X up to 10^24).

| (k, r) | B | exponent k/(k-1-r) | weights / nodes | ratio at T=10^30 |
|---|---|---|---|---|
| (3,1) | 3 | 3 | the set 0, a, a+b | — (CG Thm 9 bound 3 attained) |
| (2,0) | 4 | 2 | parallelogram 0, a, b, a+b | 1.970 (X=10^24) |
| (4,1) | 6 | 2 | u = (1,1,-1,-1) | 1.975 |
| (4 nodes, forms) | 5 | 2 | 0, α(α+β), β(α+2β), α²+2αβ+2β², (α+β)(α+2β) | 1.983 (X=10^24) |
| (4 nodes, forms) | 6 | 2 | 0, (α-β)(α-2β), α(α-2β), -2β(α-β), -αβ, α²-4αβ+2β² | 1.981 (X=10^24) |
| (7,2) | 8 | 7/4 | vectors (-3,-3),(-3,-2),(0,-1),(1,2),(1,3),(2,0),(2,1) | 1.700 |
| (5,1) | 10 | 5/3 | u = (1,1,1,-1,-1) | 1.649 |
| (8,2) | 10 | 8/5 | see `zs2_R3k8.txt` | 1.556 |
| (6,1) | 20 | 3/2 | u = ±1, three of each | 1.479 |
| (k,1) | C(k,⌊k/2⌋) | k/(k-2) | half +1, half -1 (Erdős–Sperner extremal) | — |
| cube {0,1}^n | 2^n | n/(n-1) | a_j = ∏_{i≠j} m_i | 1.489 (n=3), 1.324 (n=4) |

For the r = 2 rows, a vector (x, y) with x ≠ 0 means node y/x and weight x, after a shear that makes all x ≠ 0. For r = 2 the moment condition is just pairwise non-parallelism.

Best lower bounds on R(B), all as limits along proved families:

| B | 3 | 4–6 | 7–8 | 9–10 | 11–20 | 21–35 | 36–70 |
|---|---|---|---|---|---|---|---|
| R(B) ≥ | 3 | 2 | 7/4 | 5/3 | 3/2 | 7/5 | 4/3 |

The general rule is R(B) ≥ k/(k-2) for B ≤ C(k,⌊k/2⌋), so R(B) ≥ 1 + (2+o(1))/log_2 B (Proved).

Upper bounds:
* CG Theorem 9 (read in cg07.txt, d = 2): R(B) ≤ 2(B-1)/(B-2) for even B and ≤ 2B/(B-1) for odd B. This gives 3, 5/2, 5/2, 7/3, 7/3, 9/4, 9/4, 11/5, 11/5 for B = 4..12.
* For B = 4 there is the sharper bound below.

**Proposition 2.3 (Proved; new to the record).** For any four distinct integers, q_2 ≤ √2·L^{5/2}.

*Proof.* Reduce so that the gcd of the differences is 1. This multiplies q_2 and L by the same factor, and the bound survives because 5/2 ≥ 1.

Write the point indices as i, j, k, l. Let g = g_1 g_2 g_3 be the product of the three 2|2 split moduli, and H = h_i h_j h_k h_l the product of the four singleton-split moduli; by Lemma 1.1 all seven are pairwise coprime.

* Each difference y_i - y_j is divisible by the two 2|2 moduli of splits that put i, j on the same side, and by the singleton moduli h_k, h_l of the other two points.
* For the pair i, j, the integer (y_j - y_i)(y_k + y_l - y_i - y_j) is divisible by g·h_k h_l. Its absolute value is at most 2L².
* It vanishes only for the pairing (if any) with y_i + y_j = y_k + y_l. At most one pairing can vanish.
* The other four pairs give g⁴H² ≤ 16 L⁸.
* The product of all six differences gives g²H³ ≤ L⁶.
* Multiplying, q⁴ = g⁴H⁴ ≤ 4 L^{10}. □

The best exhaustive value is ρ = 2.199, at y = (0, 115, 444, 1410). It is inflated by the small singleton moduli 5, 6, 7. Whether 5/2 is approached is open, and irrelevant here because B = 4 is bounded.

### 2.3 Searches (Computed; exact q_2; finite search spaces only)

* **Exhaustive** over 0 = y_0 < ... < y_{B-1} = L (`exh.c`):
  * B = 4, L ≤ 1500: best 2.199;
  * B = 5, L ≤ 300: best 1.721;
  * B = 6, L ≤ 120: 1.550;
  * B = 7, L ≤ 66: 1.328;
  * B = 8, L ≤ 46: 1.088.

  These sizes are far too small to show the asymptotic exponents. For example, the B = 6 exponent-2 family reaches only 1.87 at L ≈ 10^7. The scans are reported only to exclude small surprises.
* **CRT/LLL lattice search** with random splits and random 24-bit prime moduli (`latsearch.py`, own LLL): best 1.70 (B=5), 1.57 (B=6), 1.37 (B=7), 1.34 (B=8). All are dominated by §2.2.
* **Zero-sum searches at exponent exactly 2** (k = 2r+2):
  * r = 2, k = 6: maximum B = 4, exhaustive over integer vectors in [-4,4]² (`zs2.c`).
  * r = 3, k = 8: B = 6 (nodes in {-4..5, ∞}, two-block seeding, `conic.py`).
  * r = 4, k = 10: B = 4 (same method; I am not sure whether the 1200 s run completed).
  * Together with (k, r) = (4, 1), B = 6, this gives max B = 4, 6, 4, ≥ 6, ≥ 4 for d = 1..5. There is no sign of growth towards the Plotkin bound 4d.
* **Recurrence sequences** (Fibonacci, Pell, 2^n, 3^n, squares; all B-subsets of 13-term windows; `recur.py`): every ratio ≤ 1.42. These sequences are useless here.
* **Arithmetic progressions:** q_2 = 2·step, so ρ → 1.

**Heuristic.** At exponent 2 - δ, Prop. 2.1 needs r ≈ k(1-δ)/(2-δ) moment conditions on 0/1 patterns. A Halász-type count gives at most 2^k k^{-r/2}·C^r patterns, which is < 1 once r ∝ k. This heuristic predicts R(B) → 1, as CG Conjecture 10 does. Nothing found contradicts it.

### 2.4 Characteristic 2 (Proved; computed for q = 8, 16 in `char2.py`)

Over F_q[T] with q = 2^m, take nodes = F_q and weights 1. For a subset A of F_q put Y_A = Σ_{n∈A} ∏_{l≠n}(T-l). Then Y_A(n) = [n ∈ A].

The sets A with deg Y_A ≤ q/2 are exactly ∅, F_q and the 2(q-1) affine F_2-hyperplanes. For q = 8 and q = 16 the computation lists exactly these, and no others.

This gives 2q points, two classes modulo every T - n, and degree ≤ (deg modulus)/2. That is exponent exactly 2 with B = 2q, and it attains Plotkin's 2k: it is the Sylvester–Hadamard code.

The constant field grows with m. Over the fixed field F_2[T], only the Frobenius-stable sets (∅, F_q and the two trace hyperplanes) survive, giving 4 points.

The power-sum identities behind this (Σ_{n∈U} n^s = 0 for s below |U|-1) have no characteristic-0 analogue. This is the mechanism behind the record's positive-characteristic failures (Heuristic interpretation). Hadamard/Paley profiles therefore get genuine data only in characteristic 2.

## 3. (B): divisors in a residue class

**Proposition 3.1 (Proved; exact check in `lenstra16.py`).** Let t ≥ 1, let q be divisible by t!, and let n = ∏_{i=1}^t (1+iq).
* The factors are pairwise coprime, since gcd(1+iq, 1+jq) divides j - i and q is divisible by every prime ≤ t.
* So n has at least 2^t divisors ≡ 1 mod q, and log q / log n → 1/t from below.

Checked values:

| t | q | divisors ≡ 1 mod q | log q / log n |
|---|---|---|---|
| 4 | ≈ 10^202 | 16 | 0.249574 |
| 5 | ≈ 10^202 | 32 | 0.199589 |

So Lenstra's c(α) is at least 2^{⌊1/α⌋} as α → 0 along such n, and at α = 1/4 - o(1) the count is 16.

(B) asserts a bound K for q ≥ n^{1/4-κ}. This family only forces K ≥ 16.

**Why a fifth factor does not help (Heuristic).**
* A fifth factor ≡ 1 forces α ≤ 1/5.
* Factors with other residues smaller than q must be the residues themselves.
* Larger structured counts lead to gcd(a-1, b-1) problems for S-units, where gcds are small (Bugeaud–Corvaja–Zannier type).

I found no family with more than 16 divisors at α ≥ 1/4 - κ. I did not run a systematic search for (B).

**Literature.**
* Search snippets only: Lenstra (1984) proves at most 11 divisors for s > n^{1/3}; H. Cohen has infinitely many n with 6 divisors at s > n^{1/3}; Coppersmith–Howgrave-Graham–Nagaraj (2008) give O((α-1/4)^{-3/2}) for α > 1/4; Hales, arXiv:2410.05030, surveys lower bounds, with an example of 6 divisors at α ≈ 0.358.
* I could not read Hales or Lenstra in full. arxiv.org, pdf.arxiv.org, export.arxiv.org and pub.math.leidenuniv.nl are outside this session's allowed network.
* So I cannot confirm whether boundedness below 1/4 is a stated open problem or conjecture.

## 4. (C) and the record's families

VIII.1 needs a set X with M > M_* points, M_* huge. For every distinct 8-tuple it needs admissible chain data, with Gram entries G ≤ 4λ and cofactors bounded as in (VIII.1).

For each family below, M is the number of points and λ is the level at which the VIII.1 hypotheses hold.

* **(L_5) refuter** (pairwise coprime differences; threshold/L8.md Prop. C): M = 5; λ = 0 is fine.
  * It is not two-class: each prime merges a single pair.
  * Extended to M points with balanced pair-cuts, G_ij = 1 - O(1/M) → 1. So the Gram hypothesis fails and VIII.1 does not apply.
* **(L_6) elliptic family** (Thm B): M = 6. No 8-subsets; same Gram defect if extended.
* **Golden 8-point, axis 6-point, record clusters** (growth/families.md): M ≤ 8; λ = log C / log N → 0. M is bounded, consistent with M(C) < ∞.
* **CG families** (≤ 4 points, arcs ≍ R^{1/3}): M ≤ 4; λ < 0. Bounded.
* **10-point family on arcs R^{3/5}** (referee_conditional.md): M = 10; λ = 1/20 (arc R^{1/2+2λ}). Fixed λ far above 10^-16, and M fixed.
* **Part VII Pell five-mark families**: M = 5. Unbalanced (5 of 15 cuts); not an 8-point family.
* **Fakes/Paley/Hadamard**: unbounded M; λ = 0 combinatorially.
  * There is no integer realisation. In the interval model a Hadamard profile with balanced moduli is the Theorem-9-extreme case, and §2 finds no characteristic-0 realisation even at exponent 2 beyond B = 6.
* **Characteristic-p clusters**: unbounded M; λ = 0. They are over F_{2^m}[T] (§2.4), not over Z.
* **§2.2 families**: M ≤ 6 at λ → 0. M = C(k, k/2) at λ = (k-4)/(4k), for example M = 8 at λ = 1/28 and M = 20 at λ = 1/12.

The best trade-off found over Z is M(λ) ≥ C(k, ⌊k/2⌋) with k ≈ 4/(1-4λ). This gives M = 6 as λ → 0 and growth only as λ → 1/4. That fits VIII.1, for any λ_* < 1/28 with M_* ≥ 6.

**Answer to the escape question (Proved).** Astra's Schwartz–Zippel escape from Z would turn any unbounded family meeting the hypotheses into a refutation. The difficulty is that no unbounded integer family meets the hypotheses at small λ.

## 5. Quantifiers and cautions

* Every family in §2.2 is a fixed polynomial family evaluated at T → ∞. The constants depend on the fixed nodes and weights only, so there is no template-dependent threshold.
* The exponents are limits. At finite T the ratio is below the limit; for example, 1.981 at X = 10^24 for B = 6.
* Gcd costs are not discarded. q_2 is computed exactly, and the lower bound uses the lcm with the resultant loss shown in Prop. 2.1(4).
* No fixed-S Subspace input is used anywhere.

## 6. Labels summary

**Proved**
* (A) ⇐ CG Conj 10 (d = 2).
* (L_8) ⇒ limsup R(B) < 2.
* Lemma 1.1, Prop. 2.1, Cor. 2.2, Prop. 2.3, Prop. 3.1, §2.4.
* All families in §2.2: B = 5 and B = 6 at exponent 2; R(B) ≥ k/(k-2) for B ≤ C(k,⌊k/2⌋).

**Computed** — exhaustive and LLL tables, zero-sum maxima (6, 4, ≥6, ≥4 at exponent 2), recurrence scans; all limited to the stated search spaces.

**Heuristic** — R(B) → 1 (Halász-type count); the characteristic-2 interpretation of the record's positive-characteristic clusters; no (B) family beyond 16.

**Gap** — I could not read the (B) literature beyond snippets. The r = 3 and r = 4 searches are partial (two-block seeding only).

**False** — the coordinator's parenthetical that 2 - c/B is "consistent with (L_8)"; corrected in §0 item 2.
