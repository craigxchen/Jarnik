# Referee report, block heron: Part X, sections X.1 and X.2

Refereed: answer.md lines 2151-2177 (X.1, X.2), the matching items of final10.txt, and the question in reply10.txt. Scripts (my own; x_checks.py was not imported): `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref10L_heron/triples.py` and `.../ref10L_heron/line.py`.

## Results

1. **Definition.** Stated exactly below. It agrees with reply10 (iii) and with her X.1 and X.4. Separated pairs are nonzero mod p, so gcd(d_ij, q) = h_ij holds exactly. X.1 itself never uses that the separated residues are nonzero.
2. **X.1, q^2 g^2 | H(a,b,c).** Proved, correct as stated, every sign case and p = 2 included.
3. **X.1, unit constant.** Proved, correct as stated. If H(a,b,c) != 0 then max(a,b,c) >= sqrt(q g). The inequality is strict unless g = q, and the bound is attained when g = q (the equilateral triple (q,q,q), T = 3). The case H < 0 is covered by the strict inequality H > -l^4.
4. **Corollaries.** Proved. (a) If C = 1 and q > 1 is squarefree, every triple degenerates. (b) For any C, every triple whose cut is constant at some prime p > C^2 degenerates. (c) Triples with g = 1 can be nondegenerate as soon as C >= sqrt(15/7), about 1.464 (computed example below).
5. **Exhaustive check of X.1.** Computed: no failure. It covered all 914 squarefree q <= 1500 and all locally admissible triples of positive lengths up to 3 sqrt q, 5,914 triples in all, 285 of them nondegenerate.
6. **X.2, line lemma.** Proved. Her proof is correct. It is a special case of a known theorem of Menger (1928): a metric space with more than four points embeds isometrically in R if every three of its points do. The only four-point exceptions are the pseudo-linear quadruples. The integrality and "b >= 5" parts are exactly as she states.
7. **Exhaustive check of X.2.** Computed: no counterexample. Searched: all positive symmetric integer tables with every triple degenerate on 5 points (entries <= 20; 290,700 tables), on 6 points (entries <= 14; 720,720) and on 7 points (entries <= 6; 2,520). Every one is a line metric. On 4 points (entries <= 24), 25,116 tables: 24,288 are line metrics, 828 are pseudo-linear quadruples, and nothing else occurs.

No statement in X.1 or X.2 is false. No correction is needed beyond the citation and the minor precisions below.

## 1. The definition (as refereed)

Fix:
- squarefree q > 1;
- a vertex set V with |V| = B;
- for each p | q, a two-class cut: a map beta_p : V -> {0,1} together with a residue delta_p, nonzero mod p (the difference of the two prescribed classes);
- the constant C > 0.

A *sign-free table* is a family of integers d_ij = d_ji (i != j in V) with:
- (T0) d_ij != 0 and |d_ij| <= C sqrt q;
- (T1) if beta_p(i) = beta_p(j), then d_ij = 0 mod p;
- (T2) if beta_p(i) != beta_p(j), then d_ij = epsilon delta_p mod p for some epsilon = epsilon(i,j,p) in {+1,-1}. The sign may be chosen independently for every pair and every prime.

Because delta_p != 0 mod p and q is squarefree, (T1)-(T2) say exactly that gcd(d_ij, q) = h_ij, the product of the primes at which i and j agree. This matches reply10 (iii) and her own wording in X.1 ("residue zero" / "+-delta_p != 0") and X.4 ("nonzero edge lengths, correct gcds and the critical cap").

At p = 2 one has delta_2 = 1 = -1, so the sign is vacuous. At a prime where beta_p is constant on V, delta_p plays no role.

The magnitudes are a_ij = |d_ij|. Condition (T2) is invariant under d -> -d, so the table is really a table of magnitudes: a_ij = 0 mod p on same-class pairs, and a_ij^2 = delta_p^2 mod p on separated pairs.

A triple {i,j,k} is *degenerate* if H(a,b,c) = 0 for its magnitudes a, b, c. Since a + b + c > 0, this means the largest magnitude equals the sum of the other two.

## 2. X.1: the divisibility (Proved)

Fix a triple. At each p | q its three vertices either share a class, or split 2 + 1.

- **Constant at p (p | g).** All three magnitudes are 0 mod p, so p divides all four factors of H, and p^4 | H.
- **Split at p.** Say i and j agree. Then a = a_ij = 0 mod p, and b = a_ik, c = a_jk are both +-delta_p mod p (the same delta_p, since it depends only on p). Hence b = c or b = -c mod p, and there are three cases:
  - b = c, b != -c: then a - b + c = 0 and a + b - c = 0 mod p;
  - b = -c, b != c: then a + b + c = 0 and -a + b + c = 0 mod p;
  - both hold (only possible at p = 2, where b and c are odd and a is even): all four factors are even.

  In every case p^2 | H. Nonzeroness of delta_p is not needed.

Since q is squarefree, the full power dividing H is at least prod_{p | q/g} p^2 * prod_{p | g} p^4 = q^2 g^2. The supplementary bound |H| <= (a+b+c)^4 <= 81 C^4 q^2 is correct. Hence g > 9C^2 forces H = 0. This is superseded by the unit constant below.

## 3. X.1: the unit constant (Proved)

The magnitudes are 0 mod every prime of g and g is squarefree, so g | a, b, c. Write u = a/g, v = b/g, w = c/g, and Q = q/g. Then H(u,v,w) = H(a,b,c)/g^4. Since q^2 g^2 = Q^2 g^4, this gives Q^2 | H(u,v,w), so T = H(u,v,w)/Q^2 is an integer.

**Mod 4.** H(u,v,w) = 4u^2v^2 - (u^2+v^2-w^2)^2. I verified the identity: both sides equal 2(u^2v^2 + v^2w^2 + w^2u^2) - u^4 - v^4 - w^4.
- Hence H = -X^2 mod 4, which is 0 or 3 mod 4.
- If Q is odd, then Q^2 = 1 mod 8, so T = T Q^2 = H mod 4, and T = 0 or 3 mod 4.
- If Q is even, then 2 is not constant on the triple: one magnitude is even and two are odd. All four factors are then even, so 16 | H(u,v,w). Write Q = 2Q' with Q' odd. Then Q'^2 | H/16, so 4 | H/(4 Q'^2) = T.

**Range.** Let l = w = max(u,v,w), x = u^2 and y = v^2, both in (0, l^2]. Then H = 2l^2(x+y) + 2xy - x^2 - y^2 - l^4.
- dH/dx = 2(l^2 + y - x) >= 2y > 0, and symmetrically in y. So H is strictly increasing in each variable on the box.
- The infimum over the box, -l^4, is at x = y = 0, which positivity excludes. Hence H > -l^4.
- The maximum, 3 l^4, is at x = y = l^2 (equilateral). Hence H <= 3 l^4.

**Conclusion.** If l < sqrt Q, then |H| < Q^2 in the negative direction and H < 3Q^2, so T is in {0, 1, 2}.
- On the negative side: T > -l^4/Q^2 > -1, and T is an integer, so T >= 0. This answers the question of the brief: the negative (non-triangle) case is handled by the strict lower bound.

Both congruences then force T = 0. Hence H != 0 implies g l >= g sqrt Q = sqrt(q g).

**Precisions.**
- If g < q, then Q > 1 is squarefree and not a square, so the inequality is strict.
- If g = q, equality holds for (q, q, q), with T = 3.
- The smallest nondegenerate ratio max^2/(q g) found with Q > 1 is 15/7: q = 105, g = 1, sides (7, 13, 15), H = 3 q^2.

**C = 1.** Assume q > 1 is squarefree.
- If g > 1, then max <= sqrt q < sqrt(q g).
- If g = 1, an integer side <= sqrt q is < sqrt q, because q is not a square.
- Either way H = 0. Her two cases are correct.

**General C.** Suppose the triple is constant at a prime p > C^2. Then g >= p > C^2 and max <= C sqrt q < sqrt(q p) <= sqrt(q g), so H = 0. Correct. Quantifier: this holds for each fixed C and each such triple; nothing depends on eps or A.

The computed q = 105 example shows the hypothesis "constant at some p > C^2" cannot be dropped once C >= sqrt(15/7). It is consistent with her own C = 6 example in X.5 (not refereed here).

## 4. X.1 computation (Computed)

triples.py checks every squarefree q <= 1500 (914 moduli). For each q it takes all a <= b <= c <= 3 sqrt q. A triple is admissible at p iff either all three magnitudes are 0 mod p, or exactly one is 0 mod p and the other two are nonzero and congruent up to sign (delta_p is then forced to be that common value up to sign). This is exactly the existence of a cut and a delta_p at p.

Results:
- 5,914 admissible triples, 285 nondegenerate.
- Failures of q^2 g^2 | H: 0.
- Failures of the mod-4 statement (T = 0, 3 mod 4 for odd Q; 4 | T for even Q): 0.
- Nondegenerate triples with max^2 < q g: 0.
- The ratio 1 is attained only at g = q (q = 2, 3, 5, 6, 7 with sides (q,q,q)). The next ratios are 15/7 (q = 105 and its double) and 2.4 (q = 15, sides 4, 5, 6).
- A smaller run (q <= 400, cap 2 sqrt q) also gave 0 failures.

## 5. X.2: the line lemma (Proved; known)

**Statement.** Let b >= 5 vertices carry positive symmetric integer distances with every triple degenerate. Then there are integers x_i with d_ij = |x_i - x_j|. The x_i are distinct, since the distances are positive.

**Her proof, checked step by step.**
1. Take a diametral pair (u, v) with d_uv = l and set x_i = d_ui. Degeneracy of (u, v, i) with l maximal gives d_vi = l - x_i. Positivity gives 0 < x_i < l.
2. Degeneracy of (u, i, j) gives d_ij in {|x_i - x_j|, x_i + x_j}. Degeneracy of (v, i, j) gives d_ij in {|x_i - x_j|, 2l - x_i - x_j}.
3. A pair with d_ij != |x_i - x_j| (an "exceptional" pair) therefore has x_i + x_j = l and d_ij = l. Note that x_i = x_j forces the pair to be exceptional (otherwise d_ij = 0), and then x_i = x_j = l/2.
4. Let (i, j) be exceptional and k a fifth vertex.
   - If x_k is not in {x_i, x_j}, then (i, k) and (j, k) are not exceptional (that would need x_k = x_j or x_i). Also |x_i - x_k| + |x_j - x_k| < l for x_k in (0, l) and x_i + x_j = l: the sum is |x_i - x_j|, 2x_k - l or l - 2x_k, according as x_k lies between, above or below them. This contradicts degeneracy of (i, j, k), in which l is the largest side.
   - If x_k = x_i, then x_i = x_k = l/2, so x_j = l/2, and all three pairs are exceptional with distance l. The triple (l, l, l) is not degenerate.
5. Hence no pair is exceptional, and x_i = d_ui is an integral isometric embedding.

Every step is correct. Equal coordinates of distinct vertices are excluded by positivity, as noted in step 3.

**Known result.** This is Menger's theorem on R-embeddability: a metric space with more than four points embeds isometrically in the line if each of its triples does (K. Menger, *Untersuchungen über allgemeine Metrik*, Math. Ann. 100 (1928); L. M. Blumenthal, *Theory and Applications of Distance Geometry*, 1953, section on pseudo-linear quadruples). "Every triple degenerate" implies the triangle inequality, so the table is a metric and the theorem applies. Integrality is automatic, since the coordinates are distances to one point. She should cite it as known; her proof is a correct self-contained version.

**Four points.** The proof shows exactly what happens with four points. A degenerate 4-point table is either a line metric or a pseudo-linear quadruple: a 4-cycle with sides x, l - x, x, l - x (1 <= x <= l - 1) and both diagonals l. Her unit 4-cycle is the case l = 2, x = 1. The table in her X.5 (sides 17, 29, 17, 29; diagonals 46, 46) is the case x = 17, l = 46.

## 6. X.2 computation (Computed)

line.py runs a backtracking enumeration of all labelled positive symmetric integer tables with every triple degenerate, with entries <= N. It tests line-embeddability via the coordinates from one end of a diameter. On 4 points it tests the pseudo-linear pattern under all labellings.

| points | entries <= | tables | line | pseudo-linear | other |
|---|---|---|---|---|---|
| 3 | 6 | 45 | 45 | 0 | 0 |
| 4 | 24 | 25,116 | 24,288 | 828 | 0 |
| 5 | 20 | 290,700 | 290,700 | 0 | 0 |
| 6 | 14 | 720,720 | 720,720 | 0 | 0 |
| 7 | 6 | 2,520 | 2,520 | 0 | 0 |

## 7. Remarks for the other blocks

- X.1 uses only (T1) and the "+-delta_p" part of (T2). The exact gcd condition and nonzeroness of the separated residues are used later: when the line coordinates are turned into two distinct residue classes at each prime (X.3, X.4).
- Once all triples on >= 5 vertices degenerate, X.1 + X.2 give integers y_i with |y_i - y_j| = a_ij. The table conditions then make y constant mod p on each class, with the two class residues differing by +-delta_p != 0. So the y_i are genuine two-class integers for the cut pattern, with the prescribed classes up to a reflection and translation at each p. That alignment (her line 2209) belongs to the X.3/X.4 blocks and was not refereed here.
