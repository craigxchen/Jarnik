# Round 31 referee report, block "min": Astra free thread, Part II §II.1 (exact degree-eight count)

Referee scripts are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref2_min/`. All of them are my own code. Astra's scripts were rerun separately, in `astra_rerun/`.
- `flats.cpp`: BFS census of relation flats. Every closure is computed explicitly and deduplicated by its bitmask. There is no greedy-basis rule and no threshold pruning, so it is a full census at every rank. Run as `./flats tern 11` (output `c_tern.out`, `c_tern.err`) and `./flats mid 11` (output `c_mid.out`, `c_mid.err`).
- `flats.py`: a slower Python version of the same BFS. It cross-checks ranks 1 to 3 (`tern.log`).
- `classify.py`: classification of every record up to S_8, by the minimum over all 40320 permutations of the sorted set of all relations in the flat. Here "all relations" means all 210 pairings plus all 168 midpoint relations in every orientation, not just the ordered generators. Outputs: `classes_tern.txt`, `classes_mid12.txt`, and the matching `.pkl` files.
- `realize.py`, `saturate.py`, `saturate2.py`: my nonlinear exclusion test (described in §3). Outputs: `sat_tern.txt`, `sat_mid.txt`, `sat2_*.txt`, `sat_long*.txt`.
- `solve2d.py`, `resolve.py`, `probe.py`: a two-variable variant, used to cross-check two classes (`s2d_a.txt`).
- `cube.py`: exact check of the rational cube example and of the cube code.
- `pte.py`: an independent meet-in-the-middle enumeration for the minimum-14 example.

Every run used `prlimit --as=3500000000`.

## 0. Verdict summary

| Claim of II.1 | Verdict |
|---|---|
| Reduction to parallelograms: weight 4 forces mass (2,2); distinct moment points; one pairing per support | CORRECT |
| No four collinear points on the moment surface (Vieta), and no three collinear at fixed height | CORRECT |
| Exhaustive reduction: ordering chamber, formal points, linear-image argument, discard rules, rank ≤ 5 | CORRECT |
| Enumeration certificate (greedy lexicographic basis, rank-4 skip rule) and mod-769 exactness | CORRECT |
| Ternary count: maximum 12 up to sign, only the cube flat; three planar classes, all excluded | CORRECT (independently reproduced) |
| Repeated-factor variant: again 12, same equality case; 59 planar shapes, all excluded; the hexagon rank-4 class excluded | CORRECT. Census and hexagon lemma reproduced independently. Exclusion of the 59 shapes: 52 by my own certificates, the other 7 only by Astra's rerun certificates (§5) |
| "First gap": 11 directions imply 12 | CORRECT (independently reproduced) |
| Minimum-14 example (x ∈ {1,5,8,10,12,16,22,23}, roots x+i) | CORRECT (independently reproduced) |

**Most serious issue.** None in the mathematics. The weakest point is a scope matter: everything is at D = 8 only. One reproducibility caveat for the repeated-factor variant: my own Gröbner saturation, in sympy, did not finish on 7 of the 59 planar shapes within the time I allowed (see §5). Those 7 rest on Astra's certificates, which I reran and whose logic I checked. They do not rest on fully independent code.

**Novelty.** New as a statement. The record has the cluster (clique) theorem at D = 8 (Part I §7, refereed CORRECT in `ref_ff.md`) and the cube constructions (round 29 Thm C, round 30 Thm 3.1). It has no count of all minimal-shell vectors. Round 30 `prove.md` §5 item 5 lists "any bound on the number of minimal-distance neighbours" as not obtained. II.1 bounds all weight-4 kernel vectors, whether or not they form a clique. This is a genuinely different and stronger statement than "at most 8 rows in a cluster". It also corrects two slips in the coordinator's Question 2(a):
- the count is 12 distinct divisors, not "28 edges up to sign";
- weight D/2 is only a lower bound for the lattice minimum, not "exactly the minimum".

**Significance.** Peripheral. It settles Question 2(a) completely at D = 8. It also answers Question 2(b) negatively at D = 8: no configuration of degree 8 has more weight-four divisors than the cube, so "many minimal vectors but small cliques" cannot occur at this degree. It says nothing for growing D. The growth rate (2/3+eps) log R/loglog R, the uniform bound, Questions A and B, and route 5 are all unchanged. Astra does not claim otherwise.

## 1. The reduction (checked line by line)

**Setting.** Fix 8 distinct roots ρ_i = x_i + i y_i with y_i > 0. A primitive constant-offset divisor of degree 4 is B = Π_{k∈S}(t − ρ_k^{±}) with Im B constant. By Newton's identities this means the first three power sums of its roots are real. That is the same as m ∈ {0,±1}^8 with ‖m‖_1 = 4 and Σ m_i Im ρ_i^k = 0 for k = 1, 2, 3. Rescaling the second coordinate to h_i = (y, xy, 3x²y − y³) is harmless. "Up to sign" is B versus B̄. Each such vector determines its support S, and hence the norm divisor Π_S |t − ρ_k|². The setting is as stated.

**Mass lemma** (round 30 Lemma 1.2, reused in Part I §7; I rechecked the argument). With U = Π(t−ρ)^{m+} Π(t−ρ̄)^{m−} we have U/Ū = 1 + O(t^{-4}), so deg Im U ≤ ‖m‖_1 − 4. The argument of U runs from −πΣm to 0. Hence |Σm| − 1 ≤ ‖m‖_1 − 4. For ‖m‖_1 = 4 this forces Σm = 0, using parity. The weight-4 vectors are therefore:
- e_i + e_j − e_k − e_l with distinct indices;
- in the repeated-factor variant, 2e_j − e_i − e_k;
- 2e_i − 2e_j is impossible, because h is injective (y = h_1, x = h_2/h_1).

Two pairings on one support would give h_j = h_k by subtraction. So the vectors up to sign are in bijection with the parallelograms {i,j,k,l} with h_i + h_j = h_k + h_l (and with midpoint triples, in the variant). Equivalently, they are pairs of pairs sharing a midpoint. **CORRECT.**

**No four collinear.** Write Y = h_1, V = h_2, W = h_3. Then YW = 3x²y² − y⁴ = 3V² − Y⁴, so the moment surface is S: 3V² = YW + Y⁴ with Y > 0.
- On a line with Y non-constant, V = αY + β and W = γY + δ. Then Y⁴ + (γ − 3α²)Y² + (δ − 6αβ)Y − 3β² = 0 is monic with no cubic term. Four distinct positive roots would have sum 0, which is impossible.
- At fixed height Y: if V varies the condition is quadratic in the line parameter; otherwise it is linear. So there are at most 2 points, i.e. no three collinear.

**CORRECT.** This is a shorter replacement for the Rolle argument of Part I §7, and it gives the same conclusion.

**Ordering chamber.** Take a linear functional f injective on the 8 points, and label the points by increasing f. For i < j < k < l only h_i + h_l = h_j + h_k is possible, because the other two pairings have strictly ordered f-sums. A midpoint must be the middle index. So, in this labelling, every satisfied relation is one of the 70 (or 70 + 56) ordered rows. **CORRECT.**

**The linear-image argument: is it enough? Yes.** Let E be the Q-span of the satisfied ordered rows. The Q-linear map e_i ↦ h_i (into R³) kills E, so it factors through Q^8/E. Hence the actual points are the image of the formal points under one linear map. The number of satisfied relations equals the number of candidate rows lying in E: every row in E is satisfied by linearity, and every satisfied row is in E by definition. So the configuration is counted correctly by the flat E alone.

The discard rules are sound:
- two equal formal points, or four affinely collinear formal points, stay equal or collinear under any linear map (affine dependences are preserved), contradicting distinctness or the collinearity lemma;
- they are inherited by every larger flat E' ⊇ E, because Q^8/E' is a linear image of Q^8/E.

So pruning in a search over increasing flats loses nothing. The rows have coordinate sum 0, so the formal affine dimension is 7 − dim E. Rank ≥ 6 puts all 8 points on a line and is discarded. **CORRECT.**

One more point is needed and is true: the census must count all relations lying in E, including non-ordered pairings that may happen to lie in a flat generated by ordered rows. Such a flat is infeasible for that labelling. Counting those relations only enlarges the candidate set. Astra says this ("extra chambers only enlarge the census"), and my census counts them explicitly (the `n4` column, versus `n_gen`).

**Enumeration certificate (Astra's ii_min_flat.cpp).** I checked the greedy-basis argument. Let r_{k+1} be the least row of F outside span(r_1..r_k). Every row in its proportionality class modulo the current span lies in F, so the class minimum equals r_{k+1}, and r_{k+1} > r_k. Thus each flat is visited exactly once.

The rank-4 skip rule is also correct. A row lies in span(E, r) iff its reduction modulo E is zero or proportional to the reduction of r. So the rank-5 child contains exactly the zero rows plus that one group.

The mod-769 exactness argument is correct. Every rank computed is a rank of a set of integer vectors (rows, of norm ≤ √6, and label differences, of norm √2). All have sum zero, so the rank is ≤ 7. Hadamard gives |minor| ≤ 6^{7/2} ≈ 282 < 769. I also read the code. The forward reduction with an echelon basis is correct, and the collinearity test via 2×2 minors against a fixed pivot is correct.

## 2. Independent census (different algorithm)

My `flats.cpp` does plain BFS over closures. For every admissible flat F of rank k and every generating row r ∉ F, it computes cl(F + r), deduplicates by bitmask, and tests admissibility. Admissibility here is: e_i − e_j ∉ E, and rank{E, e_j−e_i, e_k−e_i, e_l−e_i} ≤ rank E + 1 tested for all 70 quadruples. This is complete because admissibility passes to subflats. Arithmetic is mod 10007, which is exact by the same Hadamard bound. The census has no threshold, so it gives the full distribution.

**Ternary (70 rows).** Admissible flats by rank 1..5: 70, 1799, 20640, 94190, 70168. The first four equal Astra's processed-node counts exactly. The maximum count of all 210 pairings is 6 at rank 3, 12 at rank 4 and 14 at rank 5. My Python BFS agrees at ranks 1 to 3.

Records with ≥ 11 pairings: 8745, falling into only **6 S_8-classes**:

| rank | pairings | midpoints present | records with n_gen ≥ 12 (Astra's 50) | Astra's representative |
|---|---|---|---|---|
| 4 | 12 | 0 | 2 | cube |
| 5 | 14 | 5 | 28 | (0,5,12,14,16) |
| 5 | 12 | 5 | 18 | (0,5,12,14,18) |
| 5 | 12 | 6 | 4 | (2,7,14,23,27) |
| 5 | 11 | 6 | – | (1,10,17,24,29) |
| 5 | 11 | 4 | – | (0,5,12,17,30) |

This reproduces Astra's "three classes with record counts 28, 18, 4", the two rank-four nodes, and the two eleven-classes with her exact representatives. Rank-4 records with 11 ordered generators (3 of them) carry 12 pairings and are cubes. Together these give Astra's "five rank-four records" at threshold 11.

Of the 8745 records, the ones that do not respect the ordering chamber (n_gen < n4) are relabellings of the same 6 classes. This gives an extra check that the chamber reduction loses nothing.

**Cube flat.** The rank-4 class has 12 relations: 6 through the centre (from 4 concurrent pairs) and 2 per face centre. Up to S_8 it is the cube, consistent with Astra's (a,b,c,1−a−b−c) description.

**Repeated-factor variant (126 rows).** Admissible flats: 126, 5887, 124292, 1086867 at ranks 1 to 4, equal to Astra's counts. At rank 5 there are 2055914 admissible flats (11387574 closures). With threshold n_gen ≥ 12 there are 10 rank-4 and 3375 rank-5 records, equal to Astra's 10 and 3375.

Restricting to chamber-consistent records with total ≥ 12 (1246 records) and classifying up to S_8 gives:
- **2 rank-4 classes**: the cube, and a class with 9 pairings + 3 midpoints;
- **59 rank-5 classes**, with relation-count distribution 12:30, 13:15, 14:9, 15:1, 17:2, 18:1, 19:1. This is exactly Astra's table.

The second rank-4 class has three midpoint pairs about point 3 and parallelogram relations among 0..6, with point 7 free. That is exactly her "C, C±U, C±V, C±(U±V)" family.

## 3. Re-derivation of the nonlinear exclusions (my own formulation)

For a rank-5 flat E, take H, an 8×3 rational basis of ker E, normalised so that three independent rows are e_1, e_2, e_3.
- The real columns Y and V lie in ker E, so Y = Hu and V = Hv.
- u = (Y_{s1}, Y_{s2}, Y_{s3}) has positive entries. Dilate so that u = (a, b, 1).
- Translating all real parts by τ sends V → V + τY, i.e. v → v + τu. Use this to set v = (s, t, 0).
- W is eliminated: W_i = (3V_i² − Y_i⁴)/Y_i must also lie in ker E. That gives 5 equations, polynomial after clearing the Y_i.

Add z·ΠY_i − 1. Then saturate in sequence by D_ij = (Y_jV_i − Y_iV_j)² + (Y_i − Y_j)² Y_i² Y_j² = (Y_iY_j)² |ρ_i − ρ_j|², which is positive on any real configuration with distinct roots. If the resulting ideal is the unit ideal, no admissible real configuration exists. Real points survive every saturation, so the argument is sound. It is a radical-membership version of Astra's "zero remainder for a positive product" certificate.

**Ternary results (`sat_tern.txt`):**
- (0,5,12,14,16), 14 relations: unit after saturating pairs (0,1), (0,2), (0,3).
- (0,5,12,14,18): unit after (0,1), (0,2), (1,2).
- (2,7,14,23,27): after saturation the ideal is ⟨a − b − 2, 2b² + 1, ...⟩. Its zeros have b non-real, so there is no real realization.
- Eleven-classes (1,10,17,24,29) and (0,5,12,17,30): unit after (0,1), (0,2) and after (0,1), (0,2), (0,3) respectively.

The pairs used match Astra's certificates: (0,3) and (1,2) for the first two classes, and (0,2) and (0,3) for the eleven-classes.

**Astra's hand argument for the third class (checked).** I checked it directly from her displayed equations:
- Y_1 = 2a − 1 > 0, so 2a + b − 1 > 0.
- If b ≠ 1, then a = b, and the second factor becomes (a − 1)², which forces a = 1 = b, a contradiction.
- If b = 1, then point 3 = 2·point 6 − point 7 at constant height 1, which contradicts the fixed-height parabola.

This argument is correct, and so is my independent non-reality exclusion.

**Hexagon lemma (rank-4 class of the variant; checked by hand).** Translate and dilate so that C = (1, 0, −1), i.e. ρ = i. For C ± (u,v,w) on S, subtracting the two surface equations gives w = −3u − 4u³, and adding them gives v² = u²(1 − u²). Additivity of w over U, V, U+V forces (u+u')³ = u³ + u'³, i.e. uu'(u+u') = 0. The displacement of height 0 then has v = w = 0, so two points coincide. **CORRECT.**

**Attainment (`cube.py`).** Astra's q = (−9/4, −2/3, 1/4, 8/3) gives 8 distinct upper-half-plane roots. By exhaustive search over integer vectors of weight ≤ 4, the kernel has no vector of weight 1 to 3 and exactly 24 = 2·12 vectors of weight 4, all ternary. The code (15,51,60,85,170,195,204,240) has 24 edges of weight 4 (12 supports, each twice) and 4 edges of weight 8. **CORRECT.**

**Conclusion for the main claim.** The maximum is 12. A configuration attaining it has satisfied span equal to the cube flat. If it had any further relation, its span would be a rank-5 flat with ≥ 12 relations, and all of those are excluded. Since every rank-4 flat with ≥ 11 relations is the cube, 11 implies 12.

## 4. Minimum-14 example (`pte.py`)

The roots are x + i with x ∈ {1,5,8,10,12,16,22,23}, so Im ρ^k = (1, 2x, 3x² − 1). My independent meet-in-the-middle enumeration of every integer kernel vector with ‖m‖_1 ≤ 16 finds:
- none of weight < 14;
- 26 of weight 14 (13 up to sign);
- 12 of weight 16.

Astra's claims (minimum 14, thirteen shortest directions) are **CORRECT**. Her derivation is also sound: kernel vectors are equal-size multiset pairs with equal sum and sum of squares, and minimal collisions have disjoint supports. The example correctly shows that the mass lemma gives only a lower bound 4 on the minimum, and that "12" refers to the weight-4 shell.

## 5. Repeated-factor variant: status of my independent exclusion

Of the 59 planar shapes:
- 49 gave a unit ideal under my saturation (`sat_mid.txt`, `sat2_a.txt`, `sat2_b.txt`, `sat2_c.txt`).
- 3 more are excluded by non-real residual components:
  - out-index 3 is the ternary third class, with 2b² + 1 = 0;
  - out-index 15 has b² + 1 = 0;
  - out-index 14 has (b² + 1)(b⁴ + 2b³ + 4b² − 2b + 1) = 0. The quartic has no real roots (sympy real_roots).
- Out-indices 2 and 3 were also cross-checked with the two-variable method (`s2d_a.txt`).
- 7 shapes (out-indices 18, 22, 26, 31, 37, 48, 51 in `classes_mid12.txt`) timed out in sympy. Each was given 300 to 400 s, then retried for 25 to 40 min, the last attempt in `saturate4.py`. That script also saturates by the height differences of the midpoint triples, which is legitimate by the fixed-height parabola lemma. The bottleneck is the initial lex Gröbner basis. I am not allowed to install Singular or Macaulay2.

For these 7, I rely on Astra's `ii_min_all_verify.py`, which I reran successfully. It reports "Ruled out all 59 planar shapes. Rank-four survivors equal the ternary cube." Its positive-product certificate is logically sound, as above. I checked her three displayed hand exclusions for structure (factor signs, positivity identities, constant-height midpoint triples) against her script's assertions. I did not rederive their polynomials independently.

So the repeated-factor claim is CORRECT on the evidence. Its independent verification is complete for the census, for the rank-4 classes and for 52 of the 59 planar shapes. The remaining 7 are a reproducibility gap in my verification, not a gap in Astra's proof. Her certificates for them run and assert exact identities. A future referee with Singular could close it by running `saturate4.py`'s ideals there.

## 6. Minor remarks
- Astra calls her affine invariant (maximum-area triangle normalisation) "complete" without proof. This is not needed: my S_8 classification gives the same class counts.
- The text "the 28 unordered pairs therefore give twelve distinct degree-four divisors" should read: 24 of the 28 pairs have weight 4 and give 12 supports, each twice; 4 pairs have weight 8.
- Scope: characteristic 0 with an ordered real place, as for Part I §7. The Vieta, positivity and argument steps use the order on R.
