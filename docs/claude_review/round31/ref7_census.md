# Referee report, round 31, Part VII, block "census" (VII.3, constant ternary squarefree table)

Scripts and outputs: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref7_census/`.
All arithmetic (F_9 = F_3[i], polynomials over F_9 and F_3, phases, Cantor-Zassenhaus) is my own code. sympy was used only to factor large F_3 polynomials in the structured search (§5) and to cross-check factorizations. I read Astra's scripts only after my own programs gave their results.

Model: z in F_9[t] has degree D; N = z·conj(z) lies in F_3[t]; a cluster at c = 0 needs deg(z_i - z_j) <= floor(D/2). Points must be distinct. Since D/2 < D, all points of a cluster share one leading coefficient, so we may take them monic. Conjugates z and conj(z) count as different points, as they do over Z. Unit multiples have a different leading coefficient, so they never fall in the same cluster. Infinity has degree 1 and is inert.

Correction to a phrase in the task text: in this constant extension the **even-degree** F_3-primes split and the **odd-degree** ones are inert. So a squarefree norm is a product of distinct even-degree primes. Astra's text says only "every irreducible of a squarefree norm splits", which is correct.

## Verdicts in brief

| Claim | Verdict |
|---|---|
| (a) The D = 12 five-mark witness, plus all 16 witnesses in `vii_squarefree_table.out` | **CORRECT** (own exact arithmetic) |
| (b) M_3^sf(D,0) = 2,4,3,4,3,4,4,4,4,4,4,5,4,5,4,5, D = 1..16 | **CORRECT**. D = 1..12 (and 13, see §2) were reproduced by a completely different exhaustive method (head enumeration over all monic marks). D = 6..16 were reproduced by my own implementation of the triple reduction, which uses a different phase group representation and direct counting. Her internal triple and completer counts agree exactly for D = 7..16. |
| Exhaustiveness argument for >= 5 (triple with common block >= ceil(D/10), parametrisation, degree bounds, odd-order phase group) | **CORRECT**, checked line by line. Her C++ code was read and found to implement it correctly. |
| (c) Recurrence at D = 12L | **CORRECT** (proof checked; L = 2 and 3 checked numerically on all nine parents) |
| New, beyond her claim | M_3^sf(D,0) = 4, 5, 4 for **D = 17, 18, 19** (exhaustive). M_3(12,0) = 5 **for all parents** (exhaustive), with 12 parents attaining it, 9 of them squarefree. M_3(D,0) for all parents equals M_3^sf for every D <= 12. |
| (d) Six squarefree marks | **None for D <= 19** (exhaustive). None in structured families at D = 24..36 (§5). |
| (e) Multiplicity cap e_0 = 2 | It equals the squarefree value for D <= 12. At D = 12 a cap-two non-squarefree parent also attains 5. |

## 1. (a) The D = 12 witness and its structure

Script: `witness.py`, run on her `vii_squarefree_table.out`. For every line D = 1..16 it checks the following, using my own F_9[t] arithmetic and her encoding (base 9, digit a + 3b = a + bi, lowest coefficient first):
- the marks are distinct;
- every mark has degree D and is monic;
- all norms are equal and lie in F_3[t];
- N is squarefree (gcd(N, N') = 1);
- every pair difference has degree <= floor(D/2).

**All 16 lines pass.**

For D = 12, N = t^24+t^20+2t^18+t^16+2t^10+2t^2+2 (squarefree, and it equals N_0(t^2)). My own trial division gives
N = (t^2+t+2)(t^2+2t+2)(t^4+2t^3+t+1)(t^4+t^3+2t+1)(t^4+t^2+t+1)(t^4+t^2+2t+1)(t^4+2t^2+2),
which agrees with sympy. These are seven split primes with upstairs degrees 1,1,2,2,2,2,2. The prime set is closed under t -> -t.

The ten gaps are {5,5,5,6,6,6,6,6,6,6}, as she states up to the ordering of pairs.

**Cut structure.** For each base prime I computed gcd(z_j, p) in F_9[t] for every mark. Every cut is of type **2|3**. The two degree-2 primes carry the same cut, so there are 6 distinct cuts, each of upstairs weight 2 (6 × 2 = 12 = D). The minority pairs (marks 1..5) are 12, 12, 24, 45, 14, 13, 15. As a graph with six weight-2 edges this is {12, 13, 14, 15, 24, 45}: two triangles 1-2-4 and 1-4-5 sharing the edge 14, plus a pendant edge 13. The degree sequence is (4,3,2,2,1).

All nine squarefree D = 12 parents, i.e. all 24 five-clusters, have exactly this weighted profile up to S_5 (`analyse5.py`).

The phase condition reads "all vertex sums of the edge phases β_e are equal" in the group of order 3^5.

**D = 14 and D = 16.** At D = 14 there are two profiles, both all-2|3, each with 8 distinct pair cuts. At D = 16 there are 11 profiles, several with 4|1 cuts. I found no common template, and none is a reduction of the record's cube or 8-point templates. Those templates have 8 points and their cuts are not of this weight-2 book type.

**Non-squarefree parents at D = 12.** My exhaustive head census (§2) gives **M_3(12,0) = 5 for all parents too**, and nothing larger. Up to t -> ±t + a there are three 5-cluster parent orbits, of sizes 3, 6 and 3, so 12 parents in all:
- Two orbits are squarefree (9 parents). Astra's N lies in the orbit of size 3, which is the even one.
- The third orbit is not squarefree: t^24+2t^22+2t^16+2t^14+2t^10+2t^8+2t^2+1 = (t^2+t+2)^2 (t^2+2t+2)^2 × four distinct quartics. Its exponents are at most 2.

The record had M_3(D,0) <= 4 for D <= 11 (README_V, ref5_extend), which I reproduce (§2). So **D = 12 is the first five-point cluster over F_9/F_3 for all parents**, against the record's first constructed 5 at D = 568,386 (VI.4). This is new to the record.

## 2. (b) Independent exhaustive census, method 1: heads (`census.cpp`)

**Method.** Every point of a cluster has the same head (its coefficients at degrees > floor(D/2)). For each head H, the program:
- enumerates all 9^(floor(D/2)+1) low parts L;
- computes the norm as H·conj(H) + Tr(H·conj(L)) + L·conj(L), using precomputed linear tables and packed F_3 arithmetic;
- radix-sorts the norms;
- treats each run of equal norms as one maximal cluster with that head and parent.

It applies squarefreeness (gcd with the derivative) per group. Heads are reduced modulo the 12-element group generated by t -> s·t + a (s = ±1, renormalised to monic) and Galois conjugation. These maps preserve degree, contact, monicity and squarefreeness. Attaining parents are counted by expanding the canonical-head parents over their affine orbits.

This method uses no factorization, no pigeonhole and no phase group: it is plain brute force on all monic marks.

| D | all parents (mine) | squarefree (mine) | Astra sf | record (all) | # sf parents attaining | # all parents attaining |
|---|---|---|---|---|---|---|
| 1 | 2 | 2 | 2 | 2 | 3 | 3 |
| 2 | 4 | 4 | 4 | 4 | 3 | 3 |
| 3 | 3 | 3 | 3 | 3 | 1 | 1 |
| 4 | 4 | 4 | 4 | 4 | 33 | 36 |
| 5 | 3 | 3 | 3 | 3 | 30 | 45 |
| 6 | 4 | 4 | 4 | 4 | 246 | 351 |
| 7 | 4 | 4 | 4 | 4 | 3 | 3 |
| 8 | 4 | 4 | 4 | 4 | 2550 | 3450 |
| 9 | 4 | 4 | 4 | 4 | 30 | 55 |
| 10 | 4 | 4 | 4 | 4 | 23295 | 31170 |
| 11 | 4 | 4 | 4 | 4 | 144 | 288 |
| 12 | **5** | **5** | 5 | (none) | **9** | **12** |
| 13 | 4 | 4 | 4 | — | 1095 (with a 4-cluster) | 2118 |

**Agreement for every D <= 13.** At D = 13 the census runs over all 2.1e11 monic marks modulo symmetry and needs no pigeonhole, so it independently confirms the upper bound 4 there. The all-parent column reproduces the record's M_3(D,0) for D <= 11.

## 3. (b) Her exhaustiveness argument, checked step by step, and method 2 (`tri.cpp`)

1. **Only split primes.** In a squarefree norm only split primes (even degree) occur, each with exponent 1, so each mark picks one of π, conj(π). Correct.
2. **A triple with a large common block.** Every 2-colouring of 5 marks has a majority of size >= 3, so each cut is constant on at least one of the 10 triples. Weighting the cuts by upstairs degree (total D) gives some triple with common block G of degree >= ceil(D/10). Correct.
3. **Triple parametrisation** z_i = G·U_i·∏_{j≠i} conj(U_j). On a cut that is non-constant on the triple there is exactly one exceptional row, and U_i collects its primes. Correct. The blocks have squarefree, pairwise coprime norms.
4. **Contact condition.** deg(z_i - z_j) <= floor(D/2) is equivalent to ord_∞(z_i/z_j - 1) >= K = ceil(D/2). Also z_i/z_j = r_i/r_j with r = U/conj(U). The Cayley step holds: r_i/r_j - 1 = 2i(σ_i - σ_j)/unit, with σ = B/A. Hence the condition is equality of slopes mod t^-K. Correct.
5. **Degree bounds.** The numerator U_i·conj(U_j) - U_j·conj(U_i) is nonzero, so deg U_i + deg U_j >= K. With Σ deg U + deg G = D this gives max deg U_i <= D - ceil(D/10) - K and deg G <= D - ceil(3K/2). Correct.
6. **Phase group.** The norm-one truncated units mod t^-K form a group of order 3^(K-1). It is the kernel of the surjective norm from (1 + sF_9[[s]])/(s^K), of order 9^(K-1), onto (1 + sF_3[[s]])/(s^K), of order 3^(K-1). The order is odd, so doubling is invertible. Marks z, w on one parent are in one cluster iff z/conj(z) ≡ w/conj(w), because u = z/w has u·conj(u) = 1, so u/conj(u) = u^2, and u^2 ≡ 1 implies u ≡ 1. Correct. My closure computation also confirms that the group has order 3^(K-1) for K <= 10.
7. **Her code** (`vii_squarefree_table.cpp`), read after my runs:
   - the slope recursion, the tangent addition law (u+v)/(1-uv) and the factor bookkeeping are right;
   - the `valid` blocks exclude real primes and repeated or conjugate pairs;
   - the complement trick Σ_y hist[α_u + y] over G-subsets is legitimate, because subset sums are symmetric under complement;
   - for D <= 5 she uses a direct parent census;
   - `assert(max_degree <= factor_limit)` holds through D = 16.

   I found no bug.

**Method 2.** `tri.cpp` is my own implementation of the triple reduction. It represents phases as truncated unit series π/conj(π) mod s^K in the norm-one group, closed by BFS, with a full multiplication table. For each triple it builds the distribution of the products of all ± orientation choices, and then counts directly against each G (no halving trick). Results:

| D | best (mine) | triples / completers (mine) | hers |
|---|---|---|---|
| 6 | 3 (her fallback then gives 4) | 60 / 144 | 2584 / 2752 (includes fallback) |
| 7 | 4 | 12 / 24 | 12 / 24 |
| 8 | 4 | 1236 / 3732 | same |
| 9 | 4 | 600 / 2040 | same |
| 10 | 4 | 20002 / 80868 | same |
| 11 | 4 | 130 / 4680 | same |
| 12 | 5 | 11218 / 445260 | same |
| 13 | 4 | 4698 / 174600 | same |
| 14 | 5 | 167308 / 7619724 | same |
| 15 | 4 | 102316 / 4079856 | same |
| 16 | 5 | 2305474 / 119231604 | same |
| **17** | **4** | 1526174 / 69386628 | not claimed |
| **18** | **5** | 29215634 / 1666720800 | not claimed |
| **19** | **4** | 20825118 / 1073559000 | not claimed |

"best" is a valid lower bound whenever it is >= 4. Values < 5 certify that no squarefree cluster of size >= 5 exists. The upper bound 4 at D = 13, 15, 17, 19 therefore rests on the pigeonhole step (checked above) together with two independent programs (hers through 16, mine through 19).

For D = 13, method 1 needs no pigeonhole at all (see §2).

I checked the new D = 18 witnesses by reconstructing the clusters with plain polynomial products (`verify_c5.py`): size 5, gaps <= 9, squarefree N. An example D = 18 parent is t^36+2t^34+2t^26+2t^24+t^20+t^18+t^16+2t^12+t^10+1, with blocks U = t^2+(1+i), G = t^2+(2+2i), and V, T of degree 7 satisfying T(t) = -V(-t).

**Verdict on (b): CORRECT.** The table is right as stated through D = 16, and it extends to D = 17, 18, 19 as 4, 5, 4. Per D, the number of distinct squarefree parents carrying a 5-cluster is:

| D | 12 | 14 | 16 | 18 |
|---|---|---|---|---|
| parents with a 5-cluster | 9 | 12 | 69 | 99 |
| of which even (N = N_0(t^2)) | 1 | 0 | 5 | 3 |

## 4. (c) Recurrence D = 12L: CORRECT

Substituting φ = t^L (3 ∤ L) or φ = t^L + t (3 | L) into the marks and into N keeps the marks distinct, monic and of degree 12L, keeps the gaps at <= 6L, and gives norm N(φ). Squarefreeness is preserved:
- If 3 ∤ L, then (N(φ))' = L·t^(L-1)·N'(φ). Since N(0) = 2 ≠ 0, a common root α would make α^L a common root of N and N'.
- If 3 | L, then φ' = 1, and the same argument applies.

Checked numerically for φ = t^2, t^2+1, t^2+2 and t^3±t on all nine D = 12 parents. All results are squarefree with maximum exactly 5. The same argument gives 5 at all D in 14N, 16N and 18N.

## 5. (d) Trend and the search for six

- **Exhaustive data (squarefree):** no 6 for D <= 19. M^sf = 5 at every even D in [12,18] and 4 at every odd D in [7,19].
- **Structured search** (`structsearch.py` and `parentmax.py`/`conv.cpp`, an exact per-parent maximum: Cantor-Zassenhaus splitting, then phase convolution):
  - All 189 squarefree five-parents at D = 12, 14, 16, 18 under φ = t^2, t^2+1, t^2+2 (D = 24..36), and the D = 12 parents under t^3+t and t^3+2t (D = 36): all remain squarefree, and **every one has maximum exactly 5**. The substitution creates no extra marks.
  - Products of two five-parents are never squarefree in the samples tried, because they share small primes such as t^2+1 and t^2±t-1. They are useless.
  - Random search is hopeless at D >= 20, since five-parents form a fraction about 3^-2D of the parents.
- **Patterns:**
  - Some five-parents are even: Astra's D = 12 parent and five at D = 16. The D = 16 and D = 18 witnesses have T(t) = -V(-t).
  - No profile recurs across D, and no link to Rédei functions, the cube family or the 8-point templates is visible.
- **HEURISTIC count, mine.** Treat the norms of marks with a given head as random in the 3^(D+floor(D/2)+1) possible lower parts. The expected number of ordered m-tuples with equal norm is 3^E with
  E = 2(D-f-1) + 2m(f+1) - (m-1)(D+f+1), where f = floor(D/2).
  - m = 5: E = 4 for even D and E = 2 for odd D, both constant in D.
  - m = 6: E = 5 - D/2 for even D.

  So m = 5 is exactly critical in a random model. This explains why fives appear sporadically and only at even D (two fewer free trits at odd D), and it predicts no sixes. But the identical random count applies to **all** parents, where clusters are **provably unbounded** (VI.4: 5 points at D ≈ 5.7·10^5, 6 points at D ≈ 3.6·10^6). Structured families beat the random model only at astronomically large D.

  Moreover M^sf = M^all for every D <= 12, so at this scale squarefreeness is not a restriction at all. The number of five-parents also grows (9, 12, 69, 99), faster than the constant random prediction, which indicates correlations.

**Honest reading:** the data show a plateau at 5 through D = 19, with a clean parity effect. They are equally consistent with boundedness and with very slow structured growth. Since all-parent unboundedness is invisible at this scale too, **the data cannot discriminate and do not suggest boundedness**.

## 6. (e) Multiplicity cap e_0 = 2

M^(<=2)(D) lies between M^sf and M^all, which are equal for D <= 12 (and D = 13 per §2). So the cap-two table is identical: 2,4,3,4,3,4,4,4,4,4,4,5. At D = 12 a cap-two non-squarefree parent also carries 5 (§1).

## 7. Standing cautions

- **gcd:** the squarefree marks have no common factor issue; in the non-squarefree D = 12 example the degrees are the raw degrees.
- **Distinctness:** checked.
- **Equal degree, one leading coefficient, one norm in F_3[t]:** checked.
- **Conjugates:** counted as distinct points (as over Z); unit multiples are excluded by the leading coefficient.
- **Infinity:** inert (degree-1 place in a constant extension).
- **Thresholds:** none depending on m are involved.

## 8. Novelty and significance

- **Novel to the record:**
  - the exact squarefree table (the record had only all-parent values for D <= 11);
  - the first five-point cluster at D = 12, for squarefree and all parents alike;
  - my extension of the table to D = 17..19 and of the all-parent table to D = 12 (and 13).
- **Significance for the integer goal: essentially nil.**
  - These are small-degree data in a function-field model.
  - They agree with the all-parent data, which are known to grow at large D.
  - They decide neither squarefree uniformity nor growth.

  Astra says so herself (VII, "HEURISTIC (research verdict)"). The table is correct, and her exhaustiveness argument is a clean, reusable certificate method that reaches D ≈ 19 cheaply. The question of a squarefree analogue of VI.4 remains the substantive one.
