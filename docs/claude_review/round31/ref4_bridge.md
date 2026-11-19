# Referee report, block "bridge": does (IV.2) settle the record's open function-field question?

Scripts and outputs are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref4_bridge/`. All the code is mine.
* `dictionary.py` / `dictionary.out`: the exact dictionary in sympy over Q(i), with every hypothesis of (IV.2) asserted and both sides of (IV.2) printed.
* `anchored.py` / `anchored.out`: the anchored version of the star.
* `corollary.py` / `corollary.out`: the final algebra and the thresholds in D/c.

## 0. Verdict

**Yes, conditional on (IV.2).** Inequality (IV.2) applies verbatim to genuine clusters of the record's function-field model (round 29 §2.1, round 30). The dictionary has these properties:
* E_ij = 0;
* L = deg N = 2D;
* δ = K − D/2 = ⌊D/2+c⌋ − D/2 ≤ c;
* an m-point cluster gives exactly m marks.

Hence:

> **Corollary.** Over C(t), and so over Q(i)(t) and R(t):
> * M(D,0) ≤ 24 for every D ≥ 1, and M(D,0) ≤ 23 for every odd D.
> * More generally M(D,c) ≤ 24 for all D whenever 0 ≤ c < 1/2.
> * For c > 0, M(D,c) ≤ m−1 as soon as 3(m−24)D > 2m(m+3)δ(D,c) for some m ≥ 25. In particular:
>   * M(D,c) ≤ 48 for D > (5096/75)c ≈ 67.95c;
>   * M(D,c) ≤ 47 for D > 68c;
>   * M(D,c) ≤ 24 for D > (1400/3)c ≈ 466.7c.
> * Combined with the record's 2^{2c+2}(D+2): **sup_D M(D,c) ≤ max(47, 2^{2c+2}(⌊68c⌋+2)) < ∞ for every c ≥ 0.**

The record lists this uniform bound as OPEN (round 29 ff.md §0 item 5; round 30 README §0, §4; expert note §5 items 10–11). The corollary settles it over C. It does not settle it over F_{q^2}(t): BMY fails in characteristic p, and Astra's trace construction violates (IV.2) there (ref4_examples §0).

What I checked myself:
* The dictionary and its proof (§1), line by line.
* All dictionary hypotheses, exactly, on every explicit cluster I could rebuild (§3).
* The final algebra and the thresholds (`corollary.py`).

What I did not re-derive: the surface part of (IV.2), i.e. Langer BMY, effectivity, and the blow-up ledger. That part was checked line by line and with independent lattice-level Chern-number engines by the blocks `ref4_bmy` and `ref4_examples`, and both found it correct. Langer's Theorem 0.1 is confirmed only through secondary sources. Pokora (arXiv 1612.05141), as quoted in search results, gives: X normal projective, D with coefficients in [0,1], (X,D) log canonical, some multiple of K_X+D effective ⇒ (K_X+D)^2 ≤ 3e_orb(X,D). arxiv.org is not reachable from this session, so I did not read the original. **The corollary is therefore exactly as reliable as (IV.2).**

Both other Part IV referees noticed the same bridge (ref4_bmy §6, ref4_examples §0). This report supplies the rigorous dictionary, the exact tests and the consequences inside the record.

## 1. The dictionary (task a)

**Setting.**
* Cluster: z_1..z_m ∈ C[t], monic of degree D, with z_i·z̄_i = N (bar = conjugation of coefficients) and deg(z_i−z_j) ≤ K = ⌊D/2+c⌋.
* Normalisation N1 (record): for D > 2c+1 all leading coefficients agree, so the common unimodular constant divides out exactly. Nothing drifts here.

**Real roots and common factors.**
* Let r be a real root of N. The factor t−r is self-conjugate, so v_r(z_i) = v_r(z̄_i) = v_r(N)/2 for every i.
* So (t−r)^f, f = v_r(N)/2, divides every z_i. Writing z_i = (t−r)^f w_i gives a cluster at (D−f, c−f/2), since K−f = ⌊(D−f)/2 + c − f/2⌋. This has a smaller δ: δ' = δ − f/2.
* After removing all real roots, N has none. Non-real common factors do **not** need to be removed: the dictionary below accepts them.

**Orientation.** Choose U ∈ C[t] with U·Ū = N and gcd(U,Ū) = 1. For example take U = ∏_k (t−ρ_k)^{e_k}, choosing one root ρ_k from each conjugate pair, with multiplicity e_k. Such a U exists iff N has no real root.

**Star.** Put P_i := gcd(z_i, U) (monic). Then:
* z_i = P_i Ū / P̄_i. Proof: z_i | UŪ, and comparing U-parts in z_i z̄_i = UŪ gives gcd(z_i,Ū) = Ū/P̄_i.
* Concretely, if z_i = ∏(t−ρ_k)^{a_ik}(t−ρ̄_k)^{e_k−a_ik}, then P_i = ∏(t−ρ_k)^{a_ik}.

The column is (X_i, Y_i) = (Re P_i, Im P_i) := ((P_i+P̄_i)/2, (P_i−P̄_i)/(2i)) ∈ R[t]^2 ⊂ C[t]^2. Further:
* n_i = P_i P̄_i = X_i^2+Y_i^2;
* α_i = deg n_i = 2|a_i|;
* D_ij = X_iY_j − X_jY_i = (P̄_iP_j − P_iP̄_j)/(2i).

1. **Primitivity.** A common root τ of X_i and Y_i gives P_i(τ) = X_i(τ)+iY_i(τ) = 0 and P̄_i(τ) = X_i(τ)−iY_i(τ) = 0. This is impossible, since gcd(P_i, P̄_i) | gcd(U,Ū) = 1.
2. **n_i | N.** P_i | U, so P̄_i | Ū.
3. **D_ij ≠ 0.** D_ij = 0 ⟺ P_i/P̄_i = P_j/P̄_j ⟺ z_i = z_j.
4. **g_ij | D_ij, and E_ij = 0.** Write G = gcd(P_i,P_j), P_i = GA, P_j = GB. Then g_ij := gcd(n_i,n_j) = GḠ, because the P's are supported on ρ and the P̄'s on ρ̄. Also D_ij = GḠ·(ĀB − AB̄)/(2i). So g_ij | D_ij, and E_ij = deg g_ij − deg gcd(g_ij, D_ij) = 0 **for every pair**.
5. **Exact degree identity.** z_i − z_j = Ū·(P_iP̄_j − P_jP̄_i)/(P̄_iP̄_j) = −2i Ū D_ij/(P̄_iP̄_j). Hence
   deg(z_i − z_j) = D + deg D_ij − (α_i+α_j)/2, **exactly**.
   The cluster condition is therefore *equivalent* to
   deg D_ij ≤ (α_i+α_j)/2 − L/4 + δ, with L = deg N = 2D and **δ = K − D/2 = ⌊D/2+c⌋ − D/2**.
   * For integer c: δ = c when D is even, and c − 1/2 when D is odd.
   * At c = 0: δ = 0 (D even) or −1/2 (D odd).
   * Negative δ is legitimate in (IV.2). δ enters only through (IV.3), which is a plain sum of the degree bounds, and the coefficient 2(m+3)/9 of S is positive. ref4_bmy makes the same observation, and I verified the algebra in `corollary.py`.
6. **Layers.**
   * Each root ρ_k of N, and also each ρ̄_k, at each depth l ≤ e_k, is one layer, so there are 2D = L layers in all.
   * The marked set is {i : a_ik ≥ l} at ρ_k, and the *same* set at ρ̄_k, because n_i is real. So the s_h come in equal pairs.
   * Geometrically, marked sections pass through the isotropic point [1:i] over ρ_k and through [1:−i] over ρ̄_k.
   * Over C, t−ρ and t−ρ̄ are distinct primes. This is why each integer layer (p,l) becomes two layers here, and why L = 2D rather than D.
7. **Number of marks.** An m-point cluster gives **exactly m marks**. In the unanchored star the base point Ū need not belong to the cluster.
   * Anchored version (ptolemy.md Thm B(i), Part III III.1): for squarefree N, take U = z̄_0. Then P_0 = 1, i.e. column (1,0), and P_i = gcd(z_i, z̄_0) is exactly the anchored star. It still has m marks, P_0 among them.
   * For non-squarefree N the anchored star of ptolemy.md, P_i = ∏(t−ρ)^{(a_i−a_0)_+}(t−ρ̄)^{(a_0−a_i)_+}, is another valid input, with the same δ (by Thm B(i), with the same degree count).
   * Subsets of clusters are clusters with the same N, so (IV.2) may be applied to any sub-collection while keeping the parent.

Non-squarefree N, prime powers, non-real common factors and either orientation are all covered. Real roots are removed first, which only lowers δ.

## 2. The corollary (task b)

Put E = 0 and L = 2D in (IV.2). The quantity Σ_h(s_h−m/2)^2 ≥ 0 is dropped. Then an m-point cluster with m ≥ 4 forces

(3/2)(m−24)D ≤ m(m+3)δ(D,c).

* **c = 0, D even (δ = 0):** m ≤ 24. **c = 0, D odd (δ = −1/2):** m ≤ 23. Moreover:
  * 23 points need D ≥ 598/3, so M(D,0) ≤ 22 for odd D ≤ 199;
  * 22 points need D ≥ 275/3.
* **0 ≤ c < 1/2:** δ ∈ {0, −1/2}, so M(D,c) ≤ 24 for all D.
* **c > 0.** Apply the inequality to an m-subset. The cluster has at most m−1 points if D > f(m)δ, where f(m) = 2m(m+3)/(3(m−24)).
  * f is minimised at m = 49, with f(49) = 5096/75 ≈ 67.95, and f(48) = 68.
  * f(25) = 1400/3. So M(D,c) ≤ 24 once D > 466.67δ.
  * Sample values from `corollary.out` at c = 1:

    | D | 35 | 68 | 69 | 137 | 467 |
    |---|---|---|---|---|---|
    | bound on M | 42 | 48 | 28 | 25 | 24 |

* **Uniformity.** For D ≤ 68c use the record's 2^{2c+2}(D+2). This gives sup_D M(D,c) ≤ max(47, 2^{2c+2}(⌊68c⌋+2)), which is 1120 at c = 1. (IV.2) gives nothing for D < 67.95δ: the RHS then grows like m^2 and beats the LHS.

**Comparison with the record.**
* The record's bound is 2^{2c+2}(D+2), linear in D. At c = 0 it is 4(D+2), so the new bound is better for all D ≥ 5.
* Known values: sup_D M(D,0) ≥ 8 (cube), and M(8,0) = 8 (rank ≤ 4 at D = 8: round 31 ref_ff). Status now: **8 ≤ sup_D M(D,0) ≤ 24.**
* ref4_examples finds C[t] integral *fakes* (no norm condition) with 12 marks at δ = E = 0, L = 8. So 24 cannot be pushed below 12 by any argument that forgets the norm. For genuine clusters nothing above 8 is known.

## 3. Exact tests on the record's clusters (task c)

`dictionary.py` rebuilds each cluster exactly over Q(i), with roots in Q(i) or a composition. It asserts, for each cluster:
* the cluster property and z z̄ = N;
* z_i = P_iŪ/P̄_i;
* primitivity and n_i | N;
* D_ij ≠ 0 and E_ij = 0;
* the exact degree identity;
* δ_actual ≤ δ.

The layers come from factor_list of N over Q(i).

| cluster | D | K | m | L | δ (allowed = actual) | E | layer sizes | Σ(s−m/2)^2 | LHS | RHS |
|---|---|---|---|---|---|---|---|---|---|---|
| cube, y = (−9/4,−2/3,1/4,8/3) | 8 | 4 | 8 | 16 | 0 | 0 | all 4 | 0 | −192 | 0 |
| same, flipped orientation | 8 | 4 | 8 | 16 | 0 | 0 | all 4 | 0 | −192 | 0 |
| same, anchored at each of the 8 points | 8 | 4 | 8 | 16 | 0 | 0 | all 4 | 0 | −192 | 0 |
| cube ∘ (t^2+t+3) | 16 | 8 | 8 | 32 | 0 | 0 | all 4 | 0 | −384 | 0 |
| vertical, y = (1,2,3,4,5,−6,−7) (M(7,1) ≥ 8) | 7 | 4 | 8 | 14 | 1/2 | 0 | 4^10, 5^4 | 4 | −164 | 44 |
| vertical, y = (1,2,3,4,−5,−6,−7,8) (M(8,1) ≥ 14) | 8 | 5 | 14 | 16 | 1 | 0 | all 7 | 0 | −120 | 238 |
| non-squarefree e = (2,1,1), all 12 points | 4 | 3 | 12 | 8 | 1 | 0 | 4,4,6,6,6,6,8,8 | 16 | −56 | 180 |
| non-squarefree e = (3,2), all 12 points | 5 | 4 | 12 | 10 | 3/2 | 0 | 3,3,4,4,6,6,8,8,9,9 | 52 | −38 | 270 |

For the vertical lines I re-ran a brute-force search of my own and recovered the sizes 8 and 14.

**Clusters that are only certified numerically.** The D = 12 half-cube(2,1) cluster, the c = 1 rank-5 six-point clusters at D = 10 and 12, and the 5-point certified clusters cannot be rebuilt exactly: their .npy files are not in the record. Only the code and δ are needed, and these settle them:
* **Half-cube(2,1) at D = 12.** The 8 points are 0, 1, u_i, 1−u_i on columns u ∈ {0,1}^3. Each column lies in exactly 4 points, so Σ = 0, LHS = (3/4)(−16)(24) = −288 ≤ 0.
* **Automatic cases.** Since (s−m/2)^2 ≤ m^2/4, (IV.2) holds automatically whenever 3(24−m) ≥ m^2, i.e. m ≤ 7, for any code and any δ ≥ 0. With no layer of size 0 or m, it holds automatically for m = 8, since (s−4)^2 ≤ 9 ≤ 12. This covers every m ≤ 7 example (the six-point rank-5 clusters, the 5-point clusters) whatever their codes.

**Blunt conclusion of (c).**
* No cluster in the record violates (IV.2), and the dictionary passes every exact check.
* But the record's examples are **no test of (IV.2)**. Every record cluster with δ = 0 has m ≤ 8, where (IV.2) is automatic. The only larger ones have δ ≥ 1/2 and D ≤ 8, where the RHS is huge.
* (IV.2) bites only at m > 24 with D > 68δ, exactly where nothing exists or is claimed to exist.
* The real test of (IV.2) is the Chern-number verification done by ref4_bmy and ref4_examples, not these examples.
* No refutation of either (IV.2) or the dictionary arises.

## 4. Consequences inside the record (task d)

1. **Bounded rank.**
   * The code rank of a cluster is ≤ M−1, so at c = 0 the rank is ≤ 23 for every D. This is the finiteness that the record's "sup M < ∞ ⟺ bounded rank" (ff.md item 5) was waiting for.
   * **(R_4) itself (rank ≤ 4, hence M ≤ 8 by the LP) remains OPEN.** (IV.2) gives 24, not 8.
   * The record's ff §1 remark "the function-field Lang–Vojta exceptional-set problem" survives only for the sharp constant.
2. **Round 30 §6, "unbounded Q(i)(t)-clusters would refute A_5(δ) for δ < 1/30 and V_E(5)".**
   * Its premise is now false for every c. The case c > 0 needs the D > 68c clause plus the record's bound below it.
   * The conditional route to refuting A_5 from function fields is closed. Nothing is proved about A_5 in the other direction.
   * The same holds for the remark about Cilleruelo–Granville with bounded constants.
3. **"The uniform bound has exactly the shape of truncated function-field Vojta, which is itself open"** (round 30 §3.1 and §4; expert note §5 item 11). This needs reconciling, and **the inference was wrong.**
   * Round 30's Vojta statement concerns (Bl_e X_r, π^*∂ + E), an r-dimensional toric blow-up. That is correct as the predicted shape of (R_4), and it remains open.
   * Finiteness of sup_D M does not need it. Astra uses a *surface*: the ruled surface whose sections are the columns, i.e. the C(t)-point of M_{0,m} given by the m directions.
   * On that surface the needed inequality follows from a **known** theorem, the Miyaoka–Yau / Langer orbifold BMY. Miyaoka–Yau is the classical tool behind function-field Vojta-type inequalities on fibred surfaces: Vojta's canonical class inequality, and Tan's height inequality, arXiv alg-geom/9502011, via search summary.
   * So BMY on the ruled surface is *not* literally the Vojta inequality of round 30 §3.1. It is a different, known, two-dimensional inequality of Hirzebruch type. Its weights a = 2/3, 1/2 give a non-sharp constant (24 against the predicted 8).
   * Corrected statement: "the sharp bound (R_4)/limsup = 8 has the shape of open higher-dimensional truncated Vojta; a uniform bound (≤ 24) follows from orbifold BMY on a surface."
   * The record's 'open' label for the uniform bound was not an obstruction. It was a missed proof.
4. **Polynomial templates (bridge Thm 1.2).**
   * A one-parameter template z_i(t) ∈ Z[i][t] gives integer clusters on arcs C√R(t_0) with bounded C as t_0 → ∞ only if deg(z_i−z_j) ≤ D/2. That means a C(t)-cluster with c = 0, after N1.
   * So **every such template carries at most 24 points** (≤ 23 for odd D).
   * The template lower bounds for M_∞(C) can never exceed 24, and templates can never give growth.
   * This is the genuine-cluster version of Astra's own "no fixed family specialises to 25 marks". It is simpler, because the norm condition holds and E = 0 automatically.
   * Multi-parameter families are covered only after restricting them to a one-parameter curve on which the degrees behave; I did not check this.
5. **The integer problem.** **Nothing changes.**
   * The growth rate (2/3+ε) log R/loglog R, the uniform bound and Questions A and B are untouched.
   * The proof uses an actual complex surface. The arithmetic analogue would be an arithmetic BMY / Parshin–Miyaoka-type inequality on the arithmetic ruled surface P^1_Z with its sections. That is open, and plausibly of abc/Vojta strength (heuristic).
   * The record's standing caution applies: the function-field D/2 exact equations come from the infinite residue field (round 29 §2.1 item 7). (IV.2) does not even use the norm condition, and integer integral fakes are not excluded by it (Part IV's own heuristic section).

## 5. Literature (task e)

I searched for: the Cilleruelo–Granville function-field analogue; lattice points on circles / short arcs over F_q[t] and C[t]; Cilleruelo–Córdoba polynomial analogues; and sums of two squares in polynomial rings.
* I found only density and short-interval results for sums of two squares in F_q[t]: Bary-Soroker–Smilansky–Wolf (arXiv 1504.06809), Bank–Bary-Soroker–Fehm (arXiv 1509.02013), and correlations (arXiv 1701.04092).
* I also found the classical integer papers: Cilleruelo–Granville, Can. J. Math. 2009; "Lattice points on small arcs", arXiv 2107.09991.
* **I found no statement of the uniform function-field short-arc bound.** As far as I can tell, the corollary is not in the literature. This is search-based and not exhaustive.
* Hirzebruch-type inequalities for curve arrangements derived from Langer's BMY are a known genre: Pokora 1612.05141, and Q-conic papers such as 2203.11503. Their application to sections of a ruled surface with tangency chains is not something I found.

## 6. Novelty and significance

**Novelty against the record.**
* The record's only BMY computation (round4/L7.md §7) is a no-obstruction slack check for seven curves on M̄_{0,7}.
* Neither BMY on the section surface nor any uniform function-field bound appears in the record. grep for Miyaoka/Langer/BMY/Hirzebruch hits only round4, round22/lit.md, round27 dp2, papers_check.
* Astra did not draw the bridge. The dictionary is one page of routine algebra, but its consequence answers the question rounds 29–30 left open and closed "for internal work".

**Significance for the goal.**
* For the integer growth rate: **none.**
* For the record's calibration: substantial.
  * The function-field model now has a proved uniform bound, 8 ≤ sup_D M(D,0) ≤ 24.
  * The open part shrinks to the sharp constant (R_4).
  * The "open function-field Vojta" calibration of round 30 must be withdrawn as a reason for pessimism.
  * The remaining barrier to transfer is stated cleanly: the missing ingredient is an arithmetic BMY on P^1_Z, not a higher-dimensional Vojta. That is still out of reach.

## 7. Verdicts

| claim | verdict |
|---|---|
| Dictionary: n_i, P_i = gcd(z_i,U), primitive columns (Re P, Im P); E_ij = 0; δ = ⌊D/2+c⌋ − D/2; L = 2D; m marks | CORRECT (proved in §1, exact checks in §3) |
| Corollary M(D,0) ≤ 24 for all D (≤ 23 for odd D); sup_D M(D,c) < ∞ for every c | CORRECT_WITH_CORRECTIONS. Conditional on (IV.2), which I did not re-derive: it rests on Langer's Theorem 0.1, checked through secondary sources only. Exact form: M(D,c) ≤ m−1 whenever 3(m−24)D > 2m(m+3)δ, giving ≤ 48 for D > 67.95c, ≤ 47 for D > 68c, ≤ 24 for D > 466.7c, sup_D ≤ max(47, 2^{2c+2}(⌊68c⌋+2)) |
| Record clusters vs (IV.2) | CORRECT: all satisfy it with E = 0 and δ_actual = δ. But none is a real test (automatic for m ≤ 8) |
| Record claim that the uniform FF bound has the shape of open truncated Vojta | FALSE as an obstruction; correct only for (R_4)/sharp constant |
| Literature: uniform FF short-arc bound known? | NOT_CHECKED beyond searches; none found |
