# Referee report, block port: Astra Part X, section X.5 (portability and exact tests)

Scope: section X.5 of `astra_lemma_X/answer.md` (the paragraphs "why the Gaussian triangle is different", "regression tests", "two limits of exact sign recovery"), and summary items 9 and 10 of Part X. Sections X.1 to X.4 are refereed in other blocks. I use X.1 and the line lemma X.2 only where X.5 cites them.

Scripts are in `ref10L_port/` under `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/`. They are `gauss.py`, `gauss2.py`, `gauss3.py`, `census.py`, `paley43.py`, `pte.py`, `limits.py` and `erasure11.py`. All of them are my own code and use exact Python integers or `Fraction`. None imports her scripts. I read `x_checks.py` and the `pte` function of `astra_lemma_IX/ix_checks.py` so that I could match her conventions.

## 0. Verdict

**Every number in X.5 is reproduced exactly. Every identity in X.5 is correct.** I found nothing false. There are four points of wording or scope, none of which affects X.1 to X.4.

1. **The Gaussian triangle (Proved, and Computed on 180,280 triples).** Fix a squarefree q all of whose prime factors are 1 mod 4. Take three marks of norm q in one unit class. Then H(a,b,c) = 16 area^2 = (abc)^2/q = 64 g^2 (t12 t23 t13)^2, with every symbol as in X.5. Equivalently, area = 2g|t12 t23 t13|. The hypothesis "all prime factors 1 mod 4" is implicit in X.5 through IX.0.
2. **The norm-65 counterexample (Computed).** It is exactly as stated. The squared chords are 80, 52 and 4. Modulo 5 they are 0, 2 and 4, so the two separated squares 2 and 4 differ. They even have different quadratic characters.
3. **The correct circle analogue of the separated-square congruence (Proved here; new to the record as far as I know).** At a split prime p = pi pi-bar, the residue of a separated Gaussian difference is a *vertex label*, not a constant. If pi | z_i and pi-bar | z_j, then z_i - z_j = -z_j (mod pi) and z_i - z_j = z_i (mod pi-bar). Hence the matrix of separated squared chords modulo p has rank one: |z_i - z_j|^2 = -psi(z_i) phi(z_j) (mod p). Here phi and psi are the reductions modulo pi and modulo pi-bar. The twin's congruence is the special case in which all labels are equal. On a single triangle the labels impose no condition at all. The first condition they impose is on 2-by-2 cross-cut quadrilaterals.
4. **There is no circle counterpart of X.1's mechanism beyond the classical area argument (Proved).** On a circle a triple of distinct marks is never degenerate (H > 0 always). The exact identity H = 64 g^2 T^2 is the only Heron-type statement available. Its content is that area = 2g|T| >= 2g, so abc = 8 g sqrt(q) |T|. This is the Cilleruelo-Cordoba area argument sharpened by the Gaussian factorisation of IX.0. It gives max chord >= 2 g^{1/3} q^{1/6}, the Jarnik scale, and nothing near the critical arc q^{1/4}. Her sentence "the fixed-modulus argument gives no circle rate" is correct. Under her label PROVED, however, it should read as "this argument's input fails on circles". It is not a proof that no analogue exists, and item 3 shows that a weaker analogue does exist.
5. **Regression (a): reproduced exactly (Computed).** There are 182 squarefree moduli, 252 locally admissible triples, and 13 triples with nonzero Heron product. A separate set of 13 admissible triples has max^2 <= q, and all of those are degenerate. The two sets of 13 are **disjoint**, and the equality of the two counts is a coincidence. Her sentence ("including 13 with nonzero Heron product ... all 13 admissible triples at L <= sqrt q are degenerate") is correct but reads as if it described one set.
6. **Regression (b): reproduced exactly (Computed).** I rebuilt the m = 43 relaxation from its specification. With her prime recipe (P = 40000) I get 86 primes, and all 12,341 triples have H != 0 and fail q^2 g^2 | H. The number of distinct separated squares per prime ranges from 452 to 462, against 462 separated pairs per prime. With my own primes the range is 422-451 at P = 6000 and 440-458 at P = 12000. Remark: here |H| < 1.2·10^-5 q^2, so "fails q^2 g^2 | H" is the same statement as "H != 0". The range 452-462 is what one expects from 462 values that are random modulo p ≈ 5·10^4.
7. **Regression (c): reproduced (Computed), but these checks are tautological (Proved).** There are 8 families, including the 70-mark family (54,740 triples), with 55,026 triples in all. Every separated-square set is a singleton, every triple has g = gcd(q, edges), and every triple has H = 0. A genuine two-class configuration passes these tests automatically, at any length. These families are also far from the critical length: log q / log L lies between 1.07 and 1.88. So the checks test code consistency, not X.1 to X.4. She does not claim more than that, but the record should say so.
8. **The four-point table at q = 11339 (Computed and Proved).** It is a valid sign-free table with a valid two-class cut pattern. Each prime induces one of the three 2|2 splits. All four triples are degenerate, 46 < sqrt(11339) ≈ 106.5, and no signing embeds it on a line. The signatures {000, 011, 101, 110} are not Sidon. Her Sidon remark is correct. I checked it through my own derivation of the four-point case of the line lemma.
9. **The five-mark table at C = 6 (Computed).** Every number is verified: u = 11786, v = 11859, u+v = 23645, |u-v| = 73, q = 23295029, and the lengths 73, 28542, 4897. The triple (1,4,5) has H = -1152 q^2 and g = 1. The maximum length is 28542, below 6 sqrt(q) ≈ 28959. Two **additional findings**:
   * Marks 1-4 of the altered table already form the four-point exception. Opposite edges are equal: 11786, 11859 and 73.
   * The five signatures are not Sidon (s1 + s4 = s2 + s3 = 1111). Deleting mark 1 or mark 4 restores a line realisation, so one deletion suffices.

   Hence "for general C some deletion is necessary" is shown only for non-Sidon patterns. The example does not show that deletion is needed for Sidon cut patterns, which is the case X.4 is used for (Gap in scope, not an error).
10. **Paley deletion check at order 11 (Computed).** This is reproduced with my own kernel computation. There are 67 = 1 + 11 + 55 row sets and 549 = 10 + 11·9 + 55·8 annihilators. Rank T = k + 1 in every case, and every image a = Sx has the required vanishing and orthogonality properties.
11. **Reproducibility nit.** `round31/astra_lemma_X/x_checks.py` imports `ix_checks`, which is not in that directory. It runs only from `scratchpad/astra/lemma/`, where an identical copy of `ix_checks.py` exists. I checked both copies with `cmp`.

## 1. Definitions

* q is squarefree. For circles all prime factors of q are 1 mod 4. For each p | q a generator pi_p is fixed, and its conjugate is written pi-bar_p.
* A *Gaussian mark* in the unit class eta is z = eta · prod_p g_p, with g_p in {pi_p, pi-bar_p}. Its bit at p records the choice. For two marks, H_ij is the product of the common factors g_p, and h_ij = |H_ij|^2 is the product of the agreement primes. The cofactor t_ij is defined by z_j - z_i = 2 i eta H_ij t_ij. It is a nonzero rational integer (IX.0, re-verified below).
* For a triple, g is the product of the primes whose bit is constant on the triple, and T = t12 t23 t13.
* H(a,b,c) = (a+b+c)(-a+b+c)(a-b+c)(a+b-c) = 2(A B + B C + C A) - A^2 - B^2 - C^2, where A = a^2, B = b^2, C = c^2.
* phi_p : Z[i] -> F_p is reduction modulo pi_p, and psi_p is reduction modulo pi-bar_p. If phi_p sends i to r, then psi_p sends i to -r. Moreover phi_p(conj z) = psi_p(z) and N(z) = phi_p(z) psi_p(z) (mod p).
* In the twin a triple of positive lengths is *locally admissible* at p in one of two cases. Either all three lengths are 0 mod p (a constant cut). Or exactly one is 0 mod p and the other two are congruent up to sign. This is the definition used in her census, and I use it.

## 2. Part (1): the Gaussian triangle and circle analogues

### 2.1 The identity (Proved; Computed)

**Statement.** Let q > 1 be squarefree with every prime factor 1 mod 4. Let z1, z2, z3 be distinct marks in one unit class, with chord lengths a = |z1-z2|, b = |z2-z3|, c = |z1-z3|. Then:

* a^2 = 4 h12 t12^2, and similarly for the other two chords;
* H(a,b,c) = 16 area^2 = (abc)^2/q = 64 g^2 T^2;
* 2·area = |det(z2-z1, z3-z1)| = 4 g |T|.

*Proof.*

1. Write z_i = eta H beta and z_j = eta H beta-bar, where beta is the product of the separated factors of z_i. This works because at a separated prime z_j carries the conjugate factor. Then z_j - z_i = eta H (beta-bar - beta) = -2 i eta H Im(beta). So t = -Im(beta) is a rational integer, and |z_j - z_i|^2 = 4 h t^2.
2. Heron's formula gives H = 16 area^2. The circumradius is sqrt(q), so abc = 4 sqrt(q) area, and therefore 16 area^2 = (abc)^2/q.
3. Next, a^2 b^2 c^2 = 64 h12 h23 h13 T^2. In the product h12 h23 h13, a nonconstant prime agrees on exactly one pair of the triple, so it occurs once. A constant prime occurs three times. Hence h12 h23 h13 = (q/g) g^3 = q g^2, and (abc)^2/q = 64 g^2 T^2.
4. Finally t_ij != 0 because the marks are distinct. □

**Computed** (`gauss.py`). I checked the circle in the unit classes 1 and i for 10 prime sets of sizes 2 to 6, with primes up to 397. That is all 2^s marks per class and 180,280 triples. On every triple the following hold: t is real and nonzero; |d|^2 = 4ht^2; H = 4 det^2; Hq = ABC; H = 64 g^2 T^2; |det| = 4 g |T|; and H > 0. At separated primes p never divides t, and agreeing pairs have p | chord^2.

**The norm-65 example (Computed).**

* The marks are (2+i)(3+2i) = 4+7i, (2+i)(3-2i) = 8-i and (2-i)(3+2i) = 8+i, all in the unit class 1. At 5 the first two agree. At 13 the first and the third agree.
* The squared chords are 80 (marks 1,2), 52 (1,3) and 4 (2,3). Modulo 5 these are 0, 2 and 4. Modulo 13 they are 2, 0 and 4.
* The separated squares at 5 are 2 and 4. Their quadratic characters are -1 and +1. This is no accident: the Legendre symbol of a separated chord^2 is (h_ij / p), and h_ij varies.
* H = 256, area = 4 = 2·g·|T| with g = 1 and |T| = 2·1·1, and (abc)^2/q = 256. Every claim in her paragraph is confirmed.

### 2.2 What is determined modulo pi and modulo pi-bar (Proved; Computed)

Let p | q be split, and let pi | z_i and pi-bar | z_j (a pair separated at p).

* z_i - z_j = -z_j (mod pi) and z_i - z_j = z_i (mod pi-bar). Both residues are nonzero because q is squarefree. The first depends only on the mark j, the second only on the mark i.
* Write lambda_j = phi_p(z_j) for marks in the pi-bar class and mu_i = psi_p(z_i) for marks in the pi class. Then N(z_i - z_j) = phi(z_i - z_j) psi(z_i - z_j) = -mu_i lambda_j (mod p). **The separated squared chords modulo p form a rank-one (product) matrix on the cut, with nonzero vertex labels.**
* In the twin the two classes at p are points of F_p, so the labels are all equal. That gives d^2 = delta_p^2 (mod p). On the circle each class is a line {0} x F_p in Z[i]/p = F_p x F_p. The labels are free nonzero elements: they depend on the factors of the mark at all other primes and on the unit.
* **Computed** (`gauss2.py`, `gauss3.py`):
  * The rank-one identity holds on all 145,408 separated pairs tested, across four unit classes and four prime sets. The cross-ratio identity N(z_a - z_c) N(z_b - z_d) = N(z_a - z_d) N(z_b - z_c) (mod p), with a, b in the pi class and c, d in the pi-bar class, holds in 17,600 random quadruples.
  * The number of distinct labels is large. For the primes {53, 113, 257, 293, 317, 337}, the 32 marks of each class at p >= 113 all carry distinct labels. Modulo p the separated chord^2 takes up to 208 distinct values, against 1 in the twin.
* **On a single triangle the labels impose nothing (Proved).** Take a triple with one nonconstant prime p. One pair (i, i') lies in one class and the third mark j lies in the other. The two separated squares are -mu_i lambda_j and -mu_{i'} lambda_j. Since mu_i and mu_{i'} are arbitrary nonzero labels, any ordered pair of nonzero residues occurs. The first sign-free condition coming from the labels lives on a cross-cut quadrilateral (the cross-ratio identity above). So the twin's step "two separated squares agree, hence two Heron factors vanish modulo p" has no circle counterpart.

### 2.3 Heron-type divisibility on circles (Proved)

* **Exact value.** H = 64 g^2 T^2 with T a nonzero integer. Hence H > 0: three distinct concyclic points are never collinear. So the twin's conclusion "the triangle degenerates" is unavailable on circles, whatever the local input.
* **Local shadow.** At a nonconstant prime p, let A be the squared chord of the agreeing pair (so p | A), and let B, C be the two separated squared chords. The identity Hq = ABC gives -(B-C)^2 (q/p) = (A/p) B C (mod p). This holds in 811,888 prime-triple instances. B = C (mod p) holds only when p | A/p, that is, when p | t of the agreeing pair (12,064 instances). So the circle replacement for the twin's "p^2 | H at every nonconstant prime" is v_p(H) = 2 v_p(T) at nonconstant primes and v_p(H) >= 2 at constant primes. That is the factorisation H = 64 g^2 T^2 and nothing more.
* **What it gives.** area = 2g|T| >= 2g, so abc = 4 sqrt(q) area >= 8 g sqrt(q), and max chord >= 2 (g sqrt q)^{1/3}. With g = 1 this is the Cilleruelo-Cordoba bound (2·area a positive integer, so abc >= 2R with R = sqrt q), improved by the constant factor 4 that comes from restricting to one unit class. It controls arcs of length about q^{1/6}, not the critical q^{1/4}.
* **Conclusion (Proved for this mechanism).** X.1 uses a pair-level sign-free congruence to force H = 0 (mod q^2 g^2), and hence degeneracy. On circles the pair-level congruence is replaced by the rank-one label structure, which is vacuous on triangles. Degeneracy is impossible, and the only Heron-type statement is the exact area identity. Beyond the classical area argument refined by the Gaussian gcd g, X.1 has no circle counterpart. Whether the quadrilateral cross-ratio condition gives any circle rate is open (Heuristic: it is a consequence of exact concyclic coordinates and seems unlikely to add information beyond them).
* **Scale remark.** In the twin d^2 <= C^2 q, so a congruence modulo q pins d^2 down to O(C^2) values. On circles chord^2 <= C^2 sqrt(q), so even a full congruence modulo q would carry no extra information for a single chord. The twin and circle edges are not on the same scale relative to the modulus.

## 3. Part (2): regressions reproduced

### (a) Small squarefree census (Computed, `census.py`)

Quantifiers: every squarefree q with 2 <= q <= 300, and every multiset of three lengths in [1, floor(2 sqrt q)] (that is, C = 2), which is locally admissible at every p | q.

* There are 182 moduli and 252 admissible triples.
* Every admissible triple satisfies q^2 g^2 | H. The quotient H/(qg)^2 is 0 or 3 mod 4, and it is 0 mod 4 when q/g is even.
* 13 triples have H != 0. Each of them has max^2 >= qg. Equality occurs at q = 2 with lengths (2,2,2), where g = 2 and H/(qg)^2 = 3. The other twelve are (3; 1,1,3), (3; 2,2,3), (3; 3,3,3), (15; 3,5,7), (15; 4,5,6), (21; 5,7,9), (21; 6,7,8), (105; 7,10,18), (105; 7,13,15), (165; 11,20,24), (195; 13,25,27) and (231; 22,27,28).
* 13 **other** triples have max^2 <= q, and all of them are degenerate. For example (30; 2,3,5), (210; 5,9,14) and (286; 2,11,13).
* The two sets of 13 are disjoint. Her wording is accurate but invites reading them as one set.

### (b) The m = 43 relaxation (Computed, `paley43.py`; own implementation)

* **Construction.** m = 43, with the Paley sign matrix S = K - I. I verified S S^t = 44 I - J exactly. Each column carries two primes 1 mod 4 in [P, 2P], greedily paired with product in [T, T(1 + 1/(8m))], where T = (3P/2)^2. I took k = floor(P/(4m)), w_i = k i + [i = m-1] and d_ij = h_ij (w_j - w_i). This is the IX.4 specification. For P = 40000 it is her recipe, re-implemented. P = 6000 and P = 12000 are my own choices.
* **Checks.** Every edge has d^2 <= q and gcd(d, q) = h. Every triple has exactly (m-3)/4 = 10 constant columns, and g = gcd(h12, h23, h13).

| P | log q | max |d|/sqrt q | triples | H != 0 | fail q^2 g^2 \| H | distinct d^2 mod p over 462 separated pairs |
|---|---|---|---|---|---|---|
| 40000 | 946.2 | 0.163 | 12341 | 12341 | 12341 | 452-462 |
| 6000 | 783.1 | 0.159 | 12341 | 12341 | 12341 | 422-451 |
| 12000 | 842.7 | 0.161 | 12341 | 12341 | 12341 | 440-458 |

Her figures (452 to 462, and 12,341 of 12,341) are reproduced exactly at her recipe. The range depends on the primes; the earlier referee's 426-461 also falls in this pattern. Since max |H| / q^2 < 1.2·10^-5, divisibility by q^2 g^2 holds only if H = 0. So the two failure counts express one fact.

### (c) Weighted Prouhet-Tarry-Escott families (Computed, `pte.py`)

**Families.** I took the node families from `ref8L_attack.md` §2.2 at T = 1009, using her convention: node b/a, modulus aT - b, coefficient a^2. These are the 4-, 5- and 8-node one-moment families with weights (+1, ..., -1, ...) and the 7-node two-moment family with vectors (a + 7b, b). I also took the parallelogram and the 5-mark, 6-mark and cube examples at the radical of their maximal two-class modulus (Lemma 1.1 of that report).

**Results.**

| family | marks | primes | triples | log q / log L |
|---|---|---|---|---|
| 4-node | 6 | 7 | 20 | 1.563 |
| 5-node | 10 | 9 | 120 | 1.314 |
| 8-node | 70 | 13 | 54,740 | 1.066 |
| 7-node two-moment | 8 | 14 | 56 | 1.428 |
| parallelogram | 4 | 3 | 4 | 1.881 |
| 5 projective | 5 | 6 | 10 | 1.615 |
| 6 projective | 6 | 6 | 20 | 1.461 |
| cube | 8 | 4 | 56 | 1.436 |

At every prime the marks occupy at most two classes, and the separated-square set is the singleton {(r-s)^2}. On every triple g = gcd(q, edges), H = 0, and the local test of X.1 passes.

**Comment (Proved).** For any genuine two-class configuration, at any length, a separated difference is ±(r - s) modulo p and a line triple has H = 0. So these checks cannot fail. Also log q / log L < 2 throughout, so none of them is at critical length. "Squarefree part" in her text means the radical of the lcm of the node moduli.

### Paley erasure at order 11 (Computed, `erasure11.py`)

There are 67 sets of at most two missing rows. For T = [1; missing rows of S], rank T = k + 1 in every case, which gives 549 kernel basis vectors in total. Each vector x has an image a = Sx with the following properties: a is nonzero, has zero sum, and vanishes on the missing rows; and a^t S_v = 0 for every v outside supp x. This matches her count. My basis comes from row reduction, not her cofactor construction, and both bases have support at most r + 1.

## 4. Part (3): the two limiting examples

### (a) The four-point table, q = 17·23·29 = 11339 (Computed, `limits.py`; Proved)

* **The table.** d01 = 17, d23 = 17, d03 = 29, d12 = 29, d02 = 46, d13 = 46.
* **Cuts.** At 17 the zero edges are {01, 23}, giving the cut {0,1}|{2,3}. At 23 they are {02, 13}, giving {0,2}|{1,3}. At 29 they are {03, 12}, giving {0,3}|{1,2}. Each is a genuine bipartition whose internal edges are exactly the zero edges, so the cut pattern is a valid two-class pattern.
* **Separated squares.** They are {8} modulo 17, {13} modulo 23 and {28} modulo 29. Each is a singleton, so the table is sign-free with delta_p^2 equal to that value. Also gcd(d_ij, q) = h_ij on every edge (for example gcd(46, q) = 23 = h02).
* **Triples and length.** All four triples have H = 0 (17 + 29 = 46). The maximum length is 46 < sqrt(11339) ≈ 106.5, so the example works at C = 1. This does not conflict with X.1's statement for C = 1, which concerns tables with at least five marks.
* **No line realisation.** Brute force over all signings finds none. With the diameter 0-2 the coordinates are 0, 17, 46, 29, and |17 - 29| = 12 != 46.
* **Signatures.** They are 000, 011, 101 and 110 (the even-weight code). Pair sums repeat: 01 and 23, 02 and 13, 03 and 12. So they are not Sidon.
* **Her Sidon remark (Proved).** In the line-lemma setting, take a diametral pair u, v at distance l. Degeneracy puts every other vertex at an interior coordinate. A non-line pair i, j must have x_i + x_j = l and d_ij = l. So the only non-line four-point all-degenerate metric has opposite edges equal: (x, l-x, x, l-x) cyclically, with both diagonals l. Then a prime dividing one edge divides its opposite edge, so the number of zero edges at every prime is even. That excludes a 3|1 cut, which has 3 zero edges. Every bit is therefore constant or 2|2, the four signatures sum to 0 in F_2^s, and s1 + s2 = s3 + s4, which contradicts Sidon. Correct.

### (b) The five-mark table, C = 6, q = 59·67·71·83 = 23295029 (Computed)

* **Setup.**
  * a = 83 and b = 59. Then u = a(a+b) = 11786 = 2·71·83 and v = b(a+2b) = 11859 = 3·59·67. So 71 and 83 divide u, and 59 and 67 divide v; each prime of q divides exactly one of u, v.
  * The marks are 0, 11786, 11859, 23645 and 28542.
  * Each prime gives two classes: {0, 45} modulo 59, {0, 61} modulo 67, {0, 2} modulo 71 and {0, 73} modulo 83.
* **Original configuration.** It is genuine, with L = 28542 = 5.914 sqrt(q). The signatures 0000, 1100, 0011, 1111 and 1001 are distinct.
* **The change.** Replace d14 = 23645 by |u - v| = 73. If p | u, then u + v = v and u - v = -v (mod p). If p | v, both equal u (mod p). So the zero pattern and the squared residues are unchanged. Also gcd(73, q) = 1 = h14. The altered table passes my sign-free checker with the same cuts. Its maximum length is 28542 < 6 sqrt(q) ≈ 28959.
* **Triples.** The triple (1,4,5) has lengths 73, 4897 and 28542, with g = 1 and H = -1152 q^2. It is nondegenerate, so no signing realises the table on a line. The other nine triples are degenerate, and every triple satisfies q^2 g^2 | H. All her numbers are confirmed.
* **Additional findings.**
  * Marks {1,2,3,4} of the altered table are themselves the four-point exception. Opposite edges are equal (d12 = d34 = 11786, d13 = d24 = 11859, d14 = d23 = 73), and the set is not line-embeddable although all four of its triples are degenerate.
  * Of the five 4-subsets, exactly {1,2,3,5} and {2,3,4,5} are line-embeddable. So one deletion suffices.
  * The signatures are not Sidon: s1 + s4 = s2 + s3 = 1111.
* **Scope.** The example therefore shows that deletion can be necessary for some non-Sidon cut pattern at C = 6. It does not show that deletion is necessary under the Sidon hypothesis, which is what summary items 7 and 8 use. As an illustration that "X.4 allows deletions", it is fine.

## 5. Part (4): is anything in X.5 wrong or overstated?

| Claim in X.5 | Label |
|---|---|
| H = 16 area^2 = (abc)^2/q = 64 g^2 (t12 t23 t13)^2 for one unit class, with all primes 1 mod 4 (the hypothesis is implicit) | Proved; Computed |
| Chord = 2 sqrt(h)|t|, and h12 h23 h13 = q g^2 | Proved |
| Norm-65 marks, the agreement pattern at 5, and the separated squares 52 = 2 and 4 = 4 (mod 5) | Computed |
| "No general q^2 divisibility after dividing by g^2" | Proved (H/(64 g^2) = T^2) |
| "The fixed-modulus argument gives no circle rate" | Correct as a statement about this argument. Its input fails, and the circle analogue (rank-one labels, §2.2) is vacuous on triangles. Not a proof that no circle argument exists; it should be labelled as an assessment |
| Census: 182 / 252 / 13 nonzero / q^2 g^2 and modulo-4 tests / capped triples degenerate | Computed, reproduced. The two "13"s are disjoint sets |
| PTE checks, including 54,740 triples | Computed, reproduced. Tautological for genuine configurations, and not at critical length |
| Paley erasure at order 11: 67 sets and 549 annihilators | Computed, reproduced |
| m = 43: 12,341 / 12,341 failures, and 452-462 distinct squares | Computed, reproduced at her recipe. Divisibility failure is the same as H != 0; the range depends on the primes |
| Four-point table at 11339, and the Sidon exclusion | Computed; Proved |
| Five-mark table: every number, and H = -1152 q^2 | Computed |
| "This explains why X.4 allows deletions for general C" | Correct for non-Sidon cuts only. The example contains the four-point exception and is not Sidon (Gap in scope) |

Nothing in X.5 is false. The overstatements are the PROVED label on "gives no circle rate" and the implicit suggestion that the five-mark example bears on the Sidon case.
