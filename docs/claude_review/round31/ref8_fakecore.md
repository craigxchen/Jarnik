# Referee report on Part VIII, block "fakecore"

This report covers VIII.3 of Astra's answer.md, from its beginning through the end of the subsection "Forced characters, offsets and point residues" (lines 1741 to 1815 of `astra_free_VIII/answer.md`). Scripts are in `round31/ref8_fakecore/`. I wrote all of them myself and imported nothing from Astra's scripts.

## 0. Verdict

Every mathematical step in the block holds. One quantitative statement is wrong as written and needs a correction (C1). Three further points are scope remarks (C2 to C4). Astra asked for two steps to be refereed first. Both hold:

* **The mod-four support bound (uncertainty estimate): Proved.** I checked it line by line. It is the classical square-root bound for extended quadratic residue codes, transplanted to the Paley conference matrix. The proof as written is complete in every residue characteristic other than q, including characteristics 2 and 3 and inert residue fields. The bound can be sharpened slightly, to sqrt(q) + 1 - h.
* **The ordered residue witnesses: Proved.** The inequality Astra checks is exactly axiom (R) of round5/construct.md Def. 1.2, in the angle normalisation she uses: zero windings, all four units, all odd prime powers up to 2m, residues compatible across moduli. Her bounds D_u <= sqrt2 L_0^{s/2} and D_1 <= L_0|beta_i - beta_j| are correct. I verified them exactly on Gaussian integers.

| step | content | status |
|---|---|---|
| (1) | Paley matrix: K^T=-K, KK^T=qI-J, row distance q+1, class A, a=1+K, e=2, column 0 deleted | **Proved**, Computed (q = 43, 67, 127, 131, 383) |
| (2) | q-1 split primes in [P_0, P_0e^sigma] inside [P,2P]; random permutation; pair margin w - W/2 >= 2 log P - o(1); W = (2d+o(1)) q log q | **Proved** (see C4) |
| (3) | uncertainty estimate mod odd l <= h and mod 4 | **Proved**, Computed |
| (4) | random S: gcd((Kv)_{S^c}) <= 2 and the norm bound, union bounds; Vandermonde support mod q | **Proved** |
| (5) | forced characters integral; h >= 1.2q for s >= 4; (s+1)h >= q‖lambda‖² + sq/2 | **Proved**, Computed |
| (6) | offsets, ordered residue witnesses, axiom (R), (VIII.9), density/union bound | **Proved**, Computed |
| (VIII.8) | W <= (7.96+o(1)) m log m | **False as stated**. Corrected statement in C1 |

## 1. Corrections and scope remarks

**C1 (the weight coefficient; False as stated, corrected).** The construction has m = (rho/2+o(1))q marks and W = (2d+o(1)) q log q with d = 1.99. So W = (4d/rho + o(1)) m log m = (7.96/rho + o(1)) m log m. Astra fixes a real rho in (0.99, 1), so the coefficient lies in (7.96, 8.04). It is below the record's 8 only when rho > 0.995.

*Corrected statement.* For fixed rho in (0.99,1), (VIII.8) holds with W <= (7.96/rho + o(1)) m log m. The o(1) depends on rho. The stated form (7.96 + o(1)) needs rho = rho_q -> 1. The steps in this block tolerate rho_q -> 1 provided 1 - rho_q >= q^{-0.3}:
* the gcd union bound needs (1-rho) sqrt q >> q^{0.1} log q;
* the norm test needs (1-rho) q^{0.9} >> q^{0.1} log q;
* the Vandermonde count on S^c needs (1-rho) q/2 - q^{2/3} -> infinity.

I have not checked whether the later subsections (moment lifts, completion, tail) tolerate rho -> 1. They use c_rho q free coordinates on S^c. The same issue affects VIII.10, where 4d = 4·23/12 = 23/3, so its weight is (23/3)/rho. I flag this for the referee of that block.

**C2 (wording in the uncertainty proof; correct after rewording).** "Move a unit coordinate of the eigenvector to infinity" should read "move a coordinate that is nonzero modulo the chosen prime of residue characteristic l". Likewise, "-sum f = lambda y_inf is a unit" should read "is nonzero mod that prime". With this reading the argument is complete. It actually gives |supp ybar| >= sqrt q + 1, because the coordinate at infinity counts too. In the application y_inf = 0, since v_inf = 0 and (Cv)_inf = -sum v = 0. So the finite support of Kv mod l (mod 4 when l = 2) is at least sqrt q + 1 - h. This is slightly stronger than the stated sqrt q - h - 1.

**C3 (auxiliary primes; correct, asymptotic only).** Astra needs m ~ rho q/2 split primes in (2m, q^{1.01}]. The prime number theorem in progressions supplies them once q^{0.01} >~ 1.01 log q, that is for q >~ e^{650}. The statement is asymptotic, so this is legitimate. Any exponent 1+eps works, and only b_0 = (1+eps)/d < 0.51 enters later.

**C4 (the o(1) in the pair margin; correct, and robust).** The deviation term is O(sigma sqrt(q log q)) = O(q^{1.5-d}(log q)^{3/2}). It tends to zero because d > 3/2. Numerically it is still of order 1 at computable q. With q-1 consecutive split primes above q^{1.99}, randomly permuted, the minimum over all pairs of A of w(K_i-K_i') - W/2 - 2 log P is about -1.5 at q = 383, -1.1 at q = 1019 and -0.97 at q = 2003 (`margin.py`). This does not matter: the pair instance of (R) (§2.6) only needs a margin of 2 log P - O(1), because the ratio of the two sides is >> P/L_0 -> infinity.

## 2. The steps in detail

### 2.1 The matrix (Proved, Computed)

Let K_ij = chi(j-i) with chi(0) = 0 and q = 3 mod 4. Then:
* K^T = -K, because chi(-1) = -1.
* Off the diagonal, (KK^T)_{ii'} = sum_j chi(j-i)chi(j-i') = -1, since (j-i)/(j-i') takes every value except 1 on j != i'. So KK^T = qI - J.
* Row distance for i != i': write a and b for the numbers of agreements and disagreements on the q-2 columns where both rows are nonzero. Then a + b = q - 2 and a - b = -1, so b = (q-1)/2, and the distance is 2b + 2 = q + 1.
* The class A consists of the rows with chi(-i) = 1, that is the nonsquares i. It has (q-1)/2 elements, and 0 is not in A.

Every row of A has +1 on column 0, so deleting that column changes no difference of rows in A. The coefficients |K_ij - K_i'j| - 1 on I lie in {-1,0,1} and sum to (q+1) - (q-1) = 2. Also, for k, i in A, sum over j in I of K_kj K_ij = q delta_ki - 2.

`paley_basic.py` verifies all of these exactly for q = 43, 67, 127, 131, 383, together with C^2 = -qI.

### 2.2 Primes and the pair margin (Proved)

**Choice of primes.** There are about P/(2 log P) split primes in [P, 2P]. By pigeonhole on log-length, some subinterval of log-length sigma holds at least q-1 of them once sigma >= 2 log 2 · q log P / P. This is sigma << q log q / P.

**Margin.** The weighted sum sum_j c_j log p_{pi(j)} has mean 2·(average log p) >= 2 log P. A transposition changes it by at most |c_j - c_k| sigma <= 2 sigma. The permutation form of the bounded-difference inequality gives deviation <= C sigma sqrt(q log q) with probability <= 2q^{-4}. A union over at most q² pairs then works. The margin holds for all pairs of A, so S may be chosen afterwards.

**Weight.** W = 2 sum log p_j = (2d+o(1)) q log q. No column is removed by normalisation, because S almost surely meets both signs in every column.

**Consistency with round21 J1.** The primes occupy a log-window of width sigma ~ q^{-0.99} log q, far wider than the equal-weight width 1/(4r) of J1. Astra therefore does not use an equal-weight window, so J1 does not apply. The Paley coefficient sum of 2 supplies the margin, and the random permutation keeps the weight fluctuation o(1). This is the step that brings the fake inside the window.

### 2.3 The uncertainty estimate (Proved, Computed)

The checks below take v with sum v = 0, primitive, and h = ‖v‖_1 < q.

1. **The conference matrix.** The matrix C, with finite block K, C_{i,inf} = 1, C_{inf,i} = -1 and C_{inf,inf} = 0, satisfies CC^T = qI and C^T = -C, hence C² = -qI. Augmenting v by 0 at infinity gives Cv = (Kv, 0). Also (Kv)_i = sum_{j != i} ±v_j ≡ -v_i ≡ v_i mod 2, so Cv ≡ v mod 2.
2. **Integrality.** Put g = sqrt(-q). Since -q ≡ 1 mod 4, (1±g)/2 is integral. Hence y_± = (Cv - v)/2 + v(1±g)/2 are integral vectors. One checks C y_± = ±g y_± and y_+ - y_- = gv.
3. **Nonzero reduction.** Let L be a prime of Q(g) over l, with l != q. Then g is a unit mod L. Since v is primitive, gv is not ≡ 0 mod L, so at least one of y_± is nonzero mod L.
4. **Support of the reduction.**
   * For odd l, 2 is a unit, so a coordinate of y_± is ≡ 0 mod L whenever (Kv)_i ≡ 0 and v_i ≡ 0 mod l.
   * For l = 2, write y_i = (Kv_i - v_i)/2 + v_i(1±g)/2. If v_i = 0 and Kv_i ≡ 0 mod 4, then y_i = Kv_i/2 ≡ 0 mod 2. So supp(ybar) ⊂ supp(Kv mod 4) ∪ supp(v).
   * Modulo 2 alone this fails: Kv ≡ v mod 2 has support at most h. That is why mod 4 is needed. Computed: the mod-2 minimum is 2 at every tested q.
5. **The SL_2 symmetry.** Put s_i = (1,i) and s_inf = (0,1). Then C_ab = chi(det(s_a, s_b)). For gamma in SL_2(F_q), write gamma s_a = c_a s_{sigma(a)}. The signed permutation (Py)_{sigma a} = chi(c_a) y_a commutes with C, because det is gamma-invariant. The action is transitive on P^1(F_q), so a coordinate that is nonzero mod L can be moved to infinity. Computed exactly in `convid.py` for q = 11, 19, 23, 43.
6. **Fourier support.** Let f be the finite part and lambda = ±g the eigenvalue. Row infinity gives sum f = -lambda y_inf, which is nonzero mod L. The finite rows give Kf + y_inf·1 = lambda f. In characteristic 0, over Q(g, zeta_q), the Fourier transform at xi != 0 gives -chi(xi) G fhat(xi) = lambda fhat(xi), where G is the Gauss sum and G² = -q. So the nonzero frequencies of fhat lie in one quadratic class.
7. **The convolution identity.** Let n be a nonsquare. Then f * (f∘n) has transform fhat(xi) fhat(xi/n), which vanishes off xi = 0. Hence f * (f∘n) = (sum f)²/q identically. This identity has coefficients in O_{Q(g)}[1/q] and reduces mod L, where the right side is a nonzero constant at all q places. A convolution of two vectors with support b is supported on at most b² places, so b² >= q. Nothing here depends on the residue characteristic or on the size of the residue field (F_l or F_{l²}), as long as l != q. `convid.py` checks the identity exactly for q = 11, 19, 23, 43.
8. **Transfer to Kv.** Combining 4 and 7, |supp(Kv mod l)| >= sqrt q + 1 - h for odd l, and the same bound holds for Kv mod 4.

**Literature (snippets only).** This is the square-root bound for quadratic residue codes, d² - d + 1 >= q when q ≡ 3 mod 4. The reductions of the integral eigenvectors form the extended QR codes. A web search returned search-engine snippets, and I opened none of the pages: the Wikipedia page "Quadratic residue code" and a MaRDI item "On the square-root bound for QR-codes". The snippets state the bound as n <= d(d-1)+1. Astra's version is a proof of this bound in the Paley-matrix language, valid for all l != q at once. It is not new in substance, but it is correct.

**Computed: exact minimum weights of the reduced eigenvector codes** (`eigcode.py`; exhaustive). I derived the l = 2 parametrisation myself. Every y in Λ = ker(C-g) ∩ O^{q+1} has the form y = (C-1)b/2 + b omega with omega = (1+g)/2, where b in Z^{q+1} satisfies Cb ≡ b mod 2.

| q, l | type | dim | min weight | sqrt q |
|---|---|---|---|---|
| 7, 2 | split | 4 | 4 (extended Hamming) | 2.65 |
| 11, 2 | inert, F_4 | 6 over F_4 | 6 | 3.32 |
| 19, 2 | inert, F_4 | 10 over F_4 | 8 | 4.36 |
| 23, 2 | split | 12 | 8 (extended Golay) | 4.80 |
| 31, 2 | split | 16 | 8 | 5.57 |
| 11, 3 | split | 6 | 6 (ternary Golay) | 3.32 |
| 23, 3 | split | 12 | 9 | 4.80 |
| 11, 5 / 19, 5 | split | 6 / 10 | 6 / 8 | 3.32 / 4.36 |

**Computed: minimum support of Kv** (`support_small.py`). The search is exhaustive up to the affine symmetry x -> ax+b, which maps Kv to ±Kv permuted.

| q | h=2: mod 3,4,5,7 | h=4: mod 3, 4, 5, 7 | mod 2 (h=2,4) |
|---|---|---|---|
| 43 | 23 | 22, 18, 24, 24 | 2 |
| 67 | 35 | 37, 28, 38, 38 | 2 |
| 127 | 65 | **73, 56, 74, 74** | 2 |
| 131 | 67 | 74, 58, 76, 76 | 2 |

The q = 127 row reproduces Astra's values (73, 56, 74, 74) exactly. In `support_sampled.py` I also sampled 3000 vectors for each h in {6, ..., 20}, plus structured vectors (subgroup cosets, progressions, shifted squares) at q = 43 and 127. There was no violation. At these q the bound sqrt q + 1 - h is far below the observed supports, which are of order q/2 - O(sqrt q). The numerics do not stress the bound. The proof is what carries it.

### 2.4 The random sample S (Proved)

**The gcd test.** Rows outside A are never sampled. So "supp(Kv mod l) ⊂ S" has probability 0 unless that support lies inside A, and then probability rho^{|supp|}. The union over at most (2q+1)^{q^{0.1}+1} vectors and at most q^{0.1}+1 moduli tends to 0. Only odd primes l <= h matter, because nonzero entries satisfy |Kv_i| <= h. The case of an all-zero restriction is excluded by the norm test.

**The norm test.** We have ‖Kv‖² = q‖v‖². Rows outside A are retained deterministically, and each row of A with probability 1-rho. Chernoff, after scaling by h², gives failure probability <= exp(-c_rho q‖v‖²/h²) <= exp(-c_rho q/h), using ‖v‖_2² >= ‖v‖_1 = h. The deterministic part only helps (I checked the case split).

**Vandermonde.** (Kv)_i ≡ sum_j v_j (j-i)^{(q-1)/2} mod q is a polynomial in i of degree <= (q-1)/2 with unit binomial coefficients. If it vanished identically, all moments up to degree (q-1)/2 would vanish, which forces v = 0 when |supp v| < q/2. Hence Kv mod q has at least (q+1)/2 nonzero coordinates, and at least (1-rho)q/2 - o(q) of them on S^c. A random test (`forced.py`) gave minima of 23 and 34 against the claimed 22 and 34.

### 2.5 Forced characters (Proved, Computed)

**Integrality.** v_0 = sum lambda_i = 0, and Kv = -K²lambda = qlambda - J lambda = qlambda. If mu = qlambda were nonzero mod q, it would be supported on S with |S| < q/2, so K^T mu = qv would have at least (q+1)/2 nonzero coordinates mod q, which contradicts qv ≡ 0. Hence lambda is integral. The rational characters of the profile (zero-sum combinations, as in the record) are therefore exactly the integral ones, and v determines lambda uniquely. Computed: the Smith invariants of the row-difference lattice are all 1 for S = A and for a random S at q = 43 and 67.

**The inequality.** Combining ‖v‖_2² = q‖lambda‖², y² <= (s+1)y - s for integers 1 <= y <= s, and |supp v| >= (q+1)/2 gives (s+1)h >= q‖lambda‖² + sq/2 >= 3qs/2. This holds also when lambda ≡ 0 mod q, after dividing out the power of q; only the support is used. Since s is even, s = 2 means pairs, and s >= 4 gives h >= 1.2q.

**Computed** (exhaustive over s = 4, S = A): the minimum of h is 60 = 1.395q at q = 43 and 94 = 1.403q at q = 67. The inequality held in every case, and the minimum support was 24 (>= 22) and 38 (>= 34). The q = 103 run exceeded its time limit and has no result.

### 2.6 Offsets, ordered residue witnesses and axiom (R) (Proved, Computed)

**Realisation.**
* By the frame identity on I, 2K_S t = epsilon - (2/q)(sum epsilon)·1, and |t_j| <= Delta/4.
* The point phases are <2a_x - e, phi> = (pi/2)(epsilon_x + const). Their spread is <= (pi/2)Delta < e^{-W/4}, which is arc constant 1 with N = N'.
* The characters satisfy <v(c), phi> = (pi/4) c·epsilon, so v·theta = lambda·epsilon/2. All windings k_x are 0, and all point units are 1.

**Comparison with record axiom (R).** Round5/construct.md Def. 1.2 requires:
1. rho_x^(b) in (Z[i]/b)^*;
2. rho conj(rho) ≡ N mod b;
3. compatibility across b = l^a;
4. for every nonzero zero-sum c and every unit u, |sin(<v(c),phi> - arg(u)/2 - (pi/4)c·k)| >= D_u(c) e^{-w(v(c))/2}/sqrt2, over the moduli Q = odd prime powers <= Q_0 = 2M not dividing N.

Astra's data meet each item:
* rho_i = s_0 omega_i/conj(omega_i) is a unit, because the primes of N are >= P > 2m and l_i > 2m.
* It has norm s_0 conj(s_0) ≡ N.
* It is compatible across moduli, as the reduction of fixed elements.
* With k = 0, her displayed inequality is exactly (R).
* The cofactor clause for small q dividing N is vacuous.

**The bounds.** The product of the residues raised to the c_i equals A/conj A, so D_u divides A - u conj A. This is nonzero, since A has distinct, unpaired Gaussian prime factors. It is ≡ 0 mod (1+i), so sqrt2 D_u <= 2|A| <= 2L_0^{s/2}. For pairs, A - conj A = 2i Im(omega_i conj omega_j) and D_1 is odd, so D_1 | Im(omega_i conj omega_j), which is at most sqrt(l_i l_j)|sin(beta_i - beta_j)| <= L_0|beta_i - beta_j|.

**Separation of the directions.** sin(beta_i - beta_j) = Im/sqrt(l_i l_j) with a nonzero integer numerator. It is nonzero because l_i l_j is not a square. So |beta_i - beta_j| >= 1/L_0. The radius Delta/(100L_0) <= 0.01 Delta|beta_i - beta_j| gives |epsilon_i - epsilon_j| >= 0.48 Delta|beta_i - beta_j|.

**Pair instance.** At u = 1, the left side of (R) is at least |epsilon_i - epsilon_j|/2 >= 0.024 e^{-W/4}|beta_i - beta_j|. The right side is at most L_0|beta_i - beta_j| e^{-W/4} P^{-1} e^{o(1)}/sqrt2. The ratio is >= 0.034 P/L_0 -> infinity, uniformly in the pair. This is why the ordering is needed: both sides scale with |beta_i - beta_j|. For u != 1, the sine is >= 1/sqrt2 - o(1).

**Other characters (s >= 4).** I checked (VIII.9) by hand with b_0 <= 1.01/1.99 < 0.51:
* s = 4, 6: h - 0.51s >= 1.2q - 3.06.
* 8 <= s <= q/2: h - 0.51s >= (4/3 - 0.255)q = 1.078q.
* s >= q/2: h >= (8/9)(2s + q/2), so h - 0.51s >= 1.078q.
* Second inequality: h >= qs/m >= (2/rho - o(1))s, so b_0 s/h <= 0.255 + o(1) < 0.26.

Avoiding {0, 1/2} mod 1 by 2L_0^{s/2}P^{-h/2} implies (R) for all four units, since the zeros of sin((pi/2)x - arg(u)/2) lie at x ∈ {0, 1/2} + Z. It also implies both ordinary gaps.

**Density and union bound.** The density of lambda·epsilon/2 mod 1 is O(1 + L_0/Delta). The number of v with ‖v‖_1 = h is <= 2^q C(h+q, q) = e^{O(h)} for h >= 1.2q. For h <= 2q, the total failure is e^{O(q)} · L_0 2^{q/2} P^{-0.03q} -> 0. For h > 2q it is e^{O(h)} P^{-0.12h}, which is summable and tends to 0.

**Computed** (`residues.py`, exact Gaussian arithmetic). For m = 20, 40, 60 I used the first m split primes above 2m and omega_i with argument in (0, pi/4). For all pairs and 400 random characters with s in {4,6,8}, I computed D_u for all four units from the record's definition, over all odd prime powers <= 2m. In every case D_u | A - u conj A, (1+i)D_u | A - u conj A, D_u <= sqrt2 L_0^{s/2} (largest ratio 0.43), and D_1 <= |Im| <= L_0|beta_i - beta_j|. The minimum of L_0|sin(beta_i - beta_j)| was 1.000.

## 3. Quantifiers

* q -> infinity through primes ≡ 3 mod 4, with d = 1.99 fixed.
* rho is fixed in (0.99, 1); with C1 it may tend to 1 at the rate given there. The o(1) terms depend on rho and d.
* The thresholds "q sufficiently large" are effective but astronomically large (C3). They do not depend on any radius: the object is an angle model, not a sequence of circles.
* No exceptional sets, no discarded gcds, and no drifting phases occur in this block. The common phase is the explicit constant -(2/q) sum epsilon, carried exactly.

## 4. Scripts

All scripts are in `round31/ref8_fakecore/`: `paley_basic.py`, `support_small.py`, `support_sampled.py`, `eigcode.py`, `convid.py`, `forced.py`, `residues.py`, `margin.py`.
