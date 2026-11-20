# Referee report on Part VIII.3 (Lenstra's divisor problem below the 1/4 threshold)

Block: lenstra. Object: Part VIII.3 of `round31/astra_lemma_VIII/answer.md` (lines 1806-1835), together with the parts of VIII.1 it uses. Scripts: `round31/ref8L_lenstra/` (`check_viii3.py`, `margins.py`, `divsearch.c`, `structured.py`, `family.py`, outputs `out_*.txt`). I did not import or run Astra's `viii_check.py`.

## Results

1. **VIII.3, first assertion: Proved (conditional on (L_8)). Correct.** Assume (L_8). Then there are absolute constants λ > 0 and B < ∞ with the following property. Take any n ≥ 2, any q ≥ 1, any r with gcd(r,q) = 1, and any real interval [A, B_1] with 0 < A < B_1 such that q·sqrt(A B_1)/(B_1 − A) ≥ n^{1/4−λ}. Then at most B positive divisors of n lie in [A, B_1] and are congruent to r mod q. I checked every step (§1–§2). The constants depend only on the (L_8) data c, η_0, Z and the O(1).
2. **VIII.3, second assertion: Proved (conditional). Correct, with a slack constant.** For every κ with 0 < κ ≤ 2λ/3 (Astra's choice κ < min(λ/2, 1/8) is enough; the 1/8 plays no role), at most K = B·⌈1/κ⌉ positive divisors of n lie in one invertible class mod q whenever q ≥ n^{1/4−κ}. Astra's count 1 + ⌈1/κ⌉ of intervals is a valid upper bound.
3. **Numerical margins: Computed (exact rational arithmetic). Correct.** The record's Theorem D (growth/ptolemy.md §6.2, Step 4, lines 634-641) gives c = 6/69 and η_0 = 1/381 from V_8 with ε = 1/43. With a = 2·10^-5 and λ_* = 10^-16 all the inequalities of VIII.1 hold: a < η_0; 28(a/4 + λ_*) = 1.4·10^-4 < c/256 = 3.397·10^-4; and 127·sqrt(8λ_*) < a/4, since 127²·8λ_*/(a/4)² = 16129/31250 < 1. The room for the term 2/M needs M > 2.67·10^15, which is part of M_*. The admissible ranges are wider than Astra claims: λ_* < 1.94·10^-16 works for this a, and λ_* < 1.14·10^-15 works for the best a ≈ 4.85·10^-5. Then ε = 10^-16 satisfies ε/2 ≤ λ_* (VIII.2), and κ = 10^-17 < λ_*/2 (VIII.3).
4. **Literature status: open, as far as I could establish; I found no refutation. Heuristic/literature statement, not a theorem.** For q ≥ n^α with α > 1/4 the count is bounded by c(α) (Lenstra 1984; c(α) = O((α − 1/4)^-2)). Coppersmith, Howgrave-Graham and Nagaraj (2008) make this constructive and improve it to O((α − 1/4)^{-3/2}). I found no source that proves boundedness for any α ≤ 1/4, and none that proves unboundedness for any fixed α > 0. I found no construction with unboundedly many divisors in one invertible class for q ≥ n^{1/4−κ}. **Consequently Astra's implication does not refute (L_8).** The best structured lower bound I know of or found gives 2^k divisors when α < 1/k. Hence c(1/4 − κ) ≥ 16 for every κ > 0 (taking n large), which is bounded. The invertibility hypothesis cannot be dropped. Lenstra's own example (n odd: the divisors of n² that are ≡ n mod 2n correspond to the divisors of n) gives unboundedly many divisors in a non-invertible class at α ≈ 1/2.
5. **Search: Computed.** An exhaustive search over n ≤ 3·10^6, all q ≥ n^{0.2} and all invertible classes, and a structured search with large q, found nothing beyond the 2^k pattern. The records are in §4.
6. **Assessment: the implication is evidence about difficulty, not about truth.** (L_8) implies a statement that lies below the barrier of every known method (the Coppersmith/lattice threshold 1/4). The statement is not known to be false, and heuristics predict that it is true. So it neither supports nor undermines the truth of (L_8). It does show that a proof of (L_8) must contain new input on divisors in residue classes below 1/4.

## 1. The layer measure and the cofactor inequality (VIII.3, first paragraph of the proof)

*Definitions.* Write v_p for the p-adic valuation. The atoms are the layers (p, s) with p | n and 1 ≤ s ≤ v_p(n). Each layer has mass log p, so the total mass is W = log n. The sign of a divisor d on layer (p, s) is + if p^s | d and − otherwise.

* **Separating mass.** Two divisors d_i, d_j differ exactly on the layers s with min(b_i, b_j) < s ≤ max(b_i, b_j), where b = v_p. So D_ij = Σ_p |v_p(d_i) − v_p(d_j)| log p = log(d_i d_j / h_ij²) with h_ij = gcd(d_i, d_j). **Proved.** The script checks the integer identity ∏_p p^{#separating layers} · h² = d_i d_j on every pair.
* **Integrality of t_ij.** Since d_i ≡ r (mod q) and gcd(r, q) = 1, we get gcd(d_i, q) = 1, hence gcd(h_ij, q) = 1. Both q and h_ij divide d_j − d_i, so q·h_ij divides d_j − d_i. So t_ij = (d_j − d_i)/(q h_ij) is a nonzero integer for distinct divisors. **Proved.** The script asserts gcd(h, q) = 1 and the divisibility on every pair.
* **Identity.** log|t_ij| = log|d_j − d_i| − log q − log h_ij = D_ij/2 − log q + log(|d_j − d_i|/sqrt(d_i d_j)). **Proved.** The script checks t² = e^{D}·q^{-2}·(d_j − d_i)²/(d_i d_j) exactly in rationals.
* **Why (VIII.4) gives (VIII.1).** Put ρ = d_max/d_min ≥ 1. Then |d_j − d_i|/sqrt(d_i d_j) = sqrt ρ − 1/sqrt ρ, which is strictly increasing in ρ. Two divisors in [A, B_1] have ρ ≤ B_1/A, so the quantity is at most sqrt(B_1/A) − sqrt(A/B_1) = (B_1 − A)/sqrt(A B_1). This covers divisors of very different sizes inside the interval: the worst case is the pair at the two ends. Hence (VIII.4) gives log|t_ij| ≤ D_ij/2 − (1/4 − λ)W, which is (VIII.1). The lower bound 0 ≤ log|t_ij| holds because t_ij is a nonzero integer. The same bound covers the anchor cofactors t_0j = Y_j, which have the same formula. **Proved.**

*Use of VIII.1.* I re-derived VIII.1 in the form that VIII.3 uses:

* Integrality gives the Gram bound G ≤ 4λ.
* E G ≥ 0 gives E|G| ≤ 2 E G_+ ≤ 2/M + 8λ.
* The moment bound gives E m_A² ≤ E G² for even |A| ≥ 2.
* There are 127 nonempty even subsets of the 8 slots.
* By Markov's inequality, with probability at least 3/4 all 128 anchored pattern masses lie within a·w of w = W/128. This needs Σ|m_A| ≤ a, hence the condition 127·sqrt(·) < a/4.
* Then D_ij ≤ 64(1 + a)w = (1 + a)W/2, so every cofactor satisfies log|t| ≤ (a/4 + λ)W.
* The sum over the 28 pairs is below cW/256, while (L_8) gives at least cW/128 − O(1).
* The Schwartz–Zippel escape from the pullback of Z (probability at most d/M) and the repeated-row bound (at most 28/M) do not depend on n, q or the divisor set.

Balance in the (L_8) sense uses only the 119 visible blocks. The argument balances all 127, which is stronger. **Proved.**

## 2. The columns (VIII.5), the determinant identity and the chain condition

*Statement checked.* Let a = d_0 and g_i = gcd(a, d_i). Set P_i = (a/g_i, (d_i − a)/(q g_i))^T, and let n_T be the product of p over the layers whose cut, meaning the set of slots in [7] whose sign differs from slot 0, equals T. Then:

* each P_i is integral and primitive, and P_0 = (1, 0)^T;
* the slope Y_i/X_i equals (d_i/a − 1)/q, so Φ is the configuration of the d_i themselves (an affine map is in PGL_2);
* det(P_i, P_j) = a(d_j − d_i)/(q g_i g_j) = t_ij · a h_ij/(g_i g_j), and for i = 0 this equals Y_j = t_0j;
* v_p(a h_ij/(g_i g_j)) = b_0 + min(b_i, b_j) − min(b_0, b_i) − min(b_0, b_j), and this equals the number of layers at p on which both i and j are opposite slot 0, so det(P_i, P_j) = t_ij ∏_{T ∋ i,j} n_T;
* at each p, the cuts of the layers s ≤ b_0 are T_s = {i : b_i < s}. They increase with s, so they form a chain inside {b_i < b_0}. The cuts of the layers s > b_0 are T_s = {i : b_i ≥ s}. They decrease with s and form a chain inside {b_i > b_0}. These two sets are disjoint.

**Proved.** By hand, in the three cases (both b ≤ b_0: b_0 − max; both b ≥ b_0: min − b_0; b_0 between them: 0), the valuation formula agrees with the layer count. Primitivity: gcd(a/g_i, (d_i − a)/(q g_i)) divides gcd(a/g_i, d_i/g_i) = 1.

**Computed.** `check_viii3.py` verifies every item above on 1516 eight-tuples. All 28 minors are checked, and so is the chain condition (by brute-force search for a partition into two chains with disjoint unions, and separately for Astra's explicit pair of chains). The test sets are:

* 800 random instances: q ∈ {3, …, 60}; n a product of 3–6 random prime powers with exponents up to 4; with probability 0.3, n also includes a prime dividing q, so non-invertible divisors are present and must be excluded;
* 300 adversarial instances with deep layers (n = 2^12·3^9·5^6·11^3, q = 7, class 1), whose valuations at one prime differ widely;
* 416 reorderings of one fixed eight-set, so that every element serves as anchor (n = 2^6·3^5·5^3·13^2, q = 11, class 4).

All checks passed.

## 3. The second assertion

Cover [1, n] by the intervals [n^{(k−1)κ}, n^{kκ}], for k = 1, …, ⌈1/κ⌉. On each interval (B_1 − A)/sqrt(A B_1) ≤ sqrt(B_1/A) ≤ n^{κ/2}. Therefore q ≥ n^{1/4−κ} implies q·sqrt(A B_1)/(B_1 − A) ≥ n^{1/4−3κ/2}, which satisfies (VIII.4) as soon as κ ≤ 2λ/3. Divisors on a common endpoint are counted twice, which is harmless for an upper bound. **Proved.** Astra's κ < λ/2 meets this condition; the extra condition κ < 1/8 is unnecessary. Small n is absorbed into the constant, since W < W_* allows only finitely many n.

*Notational remark.* VIII.3 uses B both for the count bound and inside B_1 for the interval endpoint. This is harmless.

## 4. Numerical margins

All values are from `margins.py`, in exact rationals:

| inequality | value | holds |
|---|---|---|
| a < η_0 = 1/381 | 2·10^-5 < 2.62·10^-3 | yes |
| 28(a/4 + λ_*) < c/256, c = 6/69 | 1.4·10^-4 < 3.397·10^-4 | yes |
| 127·sqrt(8λ_*) < a/4 | ratio of squares 16129/31250 | yes |
| room for 2/M | M > 2.67·10^15 | part of M_* |
| ε/2 ≤ λ_* (VIII.2), ε = 10^-16 | | yes |
| κ ≤ 2λ_*/3, κ = 10^-17 | | yes |
| 7 − 381η_0 = 6 (ptolemy Step 4) | gives H_Σ ≥ (6w − C')/69 | yes |

The constants c = 6/69 and η_0 = 1/381 appear exactly as Astra quotes them, in `growth/ptolemy.md` lines 634-641. The O(1) term, the exceptional set Z_{1/43} and M_* are ineffective, so the margins are explicit but the counting constants B and K are not. Astra says so herself. The margins are conditional on the record's derivation V_8 ⇒ (L_8) (Theorem D, refereed in the record). I did not re-referee that derivation.

## 5. Literature: what is known for q ≥ n^α, α ≤ 1/4

*What I read.* All full-text fetches failed: arxiv.org, ams.org, pub.math.leidenuniv.nl, dominoweb.draco.res.ibm.com, cr.yp.to and oeis.org do not resolve in this session. **Everything below comes from search-engine snippets only.** The cg07.txt file is not relevant to this block, and I did not use it.

From the snippets:

* **Lenstra, Math. Comp. 42 (1984) 331–340.** If 0 ≤ r < s < n, s > n^{1/3} and gcd(r, s) = 1, then at most 11 positive divisors of n are ≡ r mod s. For α > 1/4 the count is at most c(α) when s > n^α. A snippet attributes to him the order O((α − 1/4)^-2). Cohen found infinitely many n with 6 such divisors for some s > n^{1/3}, and a search up to 3·10^6 found no example with 7. The paper also gives the non-invertible example: for odd n, the divisors of n² that are ≡ n mod 2n correspond to the divisors of n, so their number is not bounded, not even polynomially in log n.
* **Coppersmith, Howgrave-Graham and Nagaraj, Math. Comp. 77 (2008) 531–545.** They give a polynomial-time construction of all such divisors for α > 1/4 and the bound O((α − 1/4)^{-3/2}); one snippet gives the explicit form 2 + 2πα/(α − 1/4)^{3/2} + 4α/(α − 1/4), for n large in terms of α. The sentence "it is not presently known how to place a bound on the number of divisors…" refers, according to the snippet, to bounding only the divisors sx + r with x > 0 by lattice methods, as opposed to both signs. It is not a statement about α < 1/4.
* **arXiv 2410.05030, "Divisors in residue classes revisited" (Hales).** It gives lower-bound examples: 6 divisors ≡ 1 at α ≈ 0.3584, and 7 if negative divisors are counted. It says that it is open whether 6 is optimal for α ≥ 1/3, and that the situation is "much less clear" as α approaches 1/4. A snippet attributes c(1/3) = 12 to it, where Lenstra has 11. Probably the two count different things (for instance, signed divisors); I could not check this.

*Conclusion (literature, from snippets).* No source I saw proves boundedness for any α ≤ 1/4, and none constructs unboundedly many divisors in an invertible class for any fixed α > 0. Whether c(1/4) or c(1/4 − κ) exists is, as far as I can determine, **open**. It is not known to be true and not known to be false. Unboundedness is known only in the non-invertible case, which VIII.3 excludes; its proof uses gcd(r, q) = 1 essentially, through gcd(h_ij, q) = 1. **Astra's implication therefore does not refute (L_8).** A reader with network access should confirm this from the full texts, in particular the closing remarks of Lenstra 1984 and §§1 and 6 of CHN 2008.

*Own constructions (Proved/Computed).*

* **The subset family.** Let p_1, …, p_k be primes ≡ 1 mod q. Every subset product is ≡ 1, so the class of 1 contains 2^k divisors of n = ∏ p_i. Since p_i > q, α = log q/log n < 1/k, and α → 1/k from below when p_i = 1 + m_i q with small m_i. Computed instances for k = 4: q ≈ 10^3, 10^6, 10^12, 10^30 give α = 0.1749, 0.1927, 0.2227, 0.2330, each with 16 divisors (`out_structured.txt`; primality by sympy's BPSW test for p > 2^64). Hence **c(1/4 − κ) ≥ 16 for every κ > 0**. The constant K in VIII.3 must therefore be at least 16. That is consistent with Astra's claim.
* **Heuristic.** Random models predict boundedness for every fixed α > 0. We have τ(n) = n^{O(1/log log n)}, so for structure-free n the expected count in one class mod q ≈ n^α tends to 0. The structured families I examined give at most about 2^{1/α} divisors (subset products), or polynomially many in 1/α (Gaussian-type relations mod t² + 1). A small group exponent λ(q) | m does not help, because log q ≤ τ(m) log(m + 1) = m^{o(1)} while the cost is of order m. This is **heuristic**.
* **The q = 1 case of the interval statement.** The six pair products of x, x+1, x+2, x+3 divide n and satisfy (B_1 − A)/sqrt(A B_1) = 4n^{-1/4} (`out_family.txt`). This is exactly the threshold scale, with 6 divisors, so it is consistent with (VIII.4).

*Exhaustive search (Computed).* `divsearch.c` examines all n ≤ 3·10^6, all q ≥ n^{0.2} and all invertible classes. Taking the smallest divisor of each class as base, it counts the q that divide differences of divisors. Record counts by threshold θ, where α = log q/log n ≥ θ:

| θ | record | example (n, q, r, α) |
|---|---|---|
| 0.20 | 20 | 1244880, 17, 2, 0.2019 |
| 0.21 | 18 | 720720, 17, 1, 0.2101 |
| 0.22, 0.23 | 16 | 2882880, 31, 2, 0.2309 |
| 0.24 | 13 | 1441440, 31, 1, 0.2422 |
| 0.245 | 12 | 1166880, 31, 3, 0.2458 |
| 0.25 | 11 | 120120, 19, 1, 0.2517 |
| 0.26 | 11 | 1441440, 41, 1, 0.2619 |
| 0.28 | 9 | 617760, 43, 1, 0.2821 |
| 0.30 | 8 | 1767304, 75, 2, 0.3001 |
| 1/3 | 6 | 245784, 65, 1, 0.3363 |
| 0.35 | 5 | 2520, 17, 1, 0.3617 |

I re-verified five rows by direct divisor enumeration in sympy. The records stay small and fall with α. At α ≥ 1/3 they respect Lenstra's 11 and match Cohen's 6. The small-n regime (q ≤ 75) is dominated by highly composite n and is not asymptotic evidence. A randomized local search over squarefree n with q ≤ 2003 and α ≥ 0.2 did not do better than 14. Nothing found contradicts boundedness at α = 1/4 − κ.

## 6. Does the implication speak to the truth or to the difficulty of (L_8)?

**It bears on difficulty.**

* The consequence, a bound O(1) for q ≥ n^{1/4−κ}, lies exactly beyond the Coppersmith/Lenstra barrier α = 1/4. All known methods stop at that barrier: Lenstra's combinatorial argument and lattice reduction. The q = 1, √n-scale instance of the interval version shows the same barrier.
* The consequence is not known to be false. The known constructions (2^k at α < 1/k) and the heuristics are consistent with it, so the implication gives no evidence against (L_8).
* It is not known to be true either, so it gives no evidence for (L_8) beyond consistency.
* Since V_8 ⇒ (L_8) (Theorem D), it also yields "Vojta on M̄_{0,8} pushes Lenstra's threshold below 1/4". This is one more instance of the general pattern that Vojta's conjecture implies statements beyond lattice methods.
* A future construction of unboundedly many divisors in one invertible class at some fixed α ≥ 1/4 − κ (for every κ) would refute (L_8), and with it V_8 in the form Theorem D uses. No such construction is known, and none was found here.

## Labels

* VIII.3 first assertion: **Proved** (conditional on (L_8)); correct as stated.
* VIII.3 second assertion: **Proved** (conditional); correct. Minimal condition κ ≤ 2λ/3; the bound 1/8 is superfluous.
* Numerical margins (a = 2·10^-5, λ_* = 10^-16, ε = 10^-16, κ = 10^-17): **Computed**, exact; correct under the record's V_8 ⇒ (L_8).
* The identities (VIII.5), the valuation formula, the chains and the 28 minors: **Proved** and **Computed** (1516 tuples).
* Literature status at α ≤ 1/4: open, as far as the snippets show; no construction refutes (L_8). This is a **literature statement from snippets**, not checked against full texts.
* c(1/4 − κ) ≥ 16 for every κ > 0: **Proved** (subset family) and **Computed**.
* Boundedness below 1/4 is plausible: **Heuristic**.
* The attribution "These are Lenstra's weights" in VIII.3 is unverified: I could not read Lenstra's paper.
