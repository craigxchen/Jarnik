# Round 31, Part VI -- referee report, block "bounds" (VI.2 complete Hadamard cores, VI.3 nonconstant translation)

Source: `astra_free_VI/answer.md`, Part VI, sections VI.2 (lines 1295-1318), VI.3 (lines 1320-1350) and the
sharpness paragraph of VI.5 (line 1439). Record context: `README_V.md`, `ref5_extend.md` (Theorem E),
`ref5_lower.md`, `overview/overview.md` §3.1.

Scripts are in `ref6_bounds/`. They use my own F_2[t]/F_4[t] arithmetic (`f4lib.py`: integers as bit vectors,
F_4[t] as pairs (A,B) = A + theta B). Nothing was imported from Astra or the other referees. sympy is used only
to factor small real polynomials in the per-prime check. Outputs are in `*.out`.

| script | purpose |
|---|---|
| `f4lib.py`, `clusters.py` | arithmetic; Part V clusters; my own CRT construction of the VI.5 squarefree clusters |
| `test_vi2.py` → `test_vi2.out` | VI.2 on Part V (m=4,8,16) and VI.5 (m=4,8,16): exact products, depths, s_v, B_v, every inequality |
| `test_vi3.py` → `vi3_d3.out`, `vi3_d4.out`, `vi3_d4rand.out` | VI.3 on all side choices: (VI.3) for every instance; for critical ones, the difference identity, the inert bound and the margin identity |
| `test_vi3_prime.py` → `vi3_prime.out` | the per-prime range inequality of VI.3 (29,557 primes); depth amplification (VI.2) |
| `minsearch.py`, `exh4.py`, `exh5.py`, `exh6.py` → `exh4a/b.out`, `exh5_all_*.out`, `exh6a/b/c.out` | exact pruned search (within the stated T ranges) for the least reduced degree of all-form designs with translation B |
| `verify_gen.py` | from-scratch verification of any found design (d=5, 6) |
| `verify68.py`, `verify_d3.py` | independent re-verification of the new small clusters (norms, distinctness, all pair differences) |

## 0. Verdicts in brief

| claim | verdict |
|---|---|
| VI.2 exact products prod z_i^{chi_v} = r_v^{m/2} and prod z_i^{(chi_v-chi_w)/2} = (r_v/r_w)^{m/4} | **CORRECT** (proved line by line; checked as exact polynomial identities on 6 clusters) |
| VI.2 depth calculus (2^a b multiplies depth by 2^a exactly) | **CORRECT** (proof checked; 0 failures in 7,200 random tests) |
| VI.2: m is a power of two; depth(r_v) >= D/m; depth(r_v/r_w) >= 2D/m | **CORRECT** |
| VI.2: deg B_v <= s_v, deg(A_vB_w + A_wB_v) <= s_v + s_w, the counting, (VI.2) | **CORRECT**, with the hypotheses made explicit in §1.4. The disjointness hypothesis can be weakened (§1.4). |
| VI.5 "sharpness for VI.2" (Theta(m^2 log m) within complete Hadamard cores) | **CORRECT** (the VI.5 clusters satisfy every VI.2 hypothesis; checked m=4,8,16) |
| VI.3 (VI.3): D' >= (1/2) sum_var h_v >= (Q/8 - 1/2) deg B | **CORRECT** (proof checked; 0 failures on 618,000 instances and 29,557 single primes) |
| VI.3 second half (Theorem E with P_J(B); margin identity) and the conclusion D' = Omega(s^2) | **CORRECT** (asymptotic, with no explicit threshold given; the constant is about 1/400) |
| Summary line 6 ("full hyperplane designs with nonconstant translation have degree Omega(m^2)") | **CORRECT**, but only for designs that use **all** Q/2 forms. My data show that nonconstant B really does beat B = 1 (§2.5). |

**New computational facts (mine; exact, independently re-verified).**
1. *8 points at point degree 68.* V = F_2[t]_{<4}, B = t^2+t, Y = t^11 + theta(t^2+t), all 8 forms with l(B)=1, and
   sides (0,1,0,1,1,0,1,0). This gives 8 distinct monic points of degree 68 with a common norm of degree 136 and
   maximum pair difference 34 = D/2. The record's best 8-point cluster had point degree 104 (`ref5_lower.md`, B = 1),
   and over all T with deg T <= 16 the exact minimum for B = 1 is 104 (my search reproduces it).
2. *16 points at point degree 472* (B = t^3+t or t^3+t^2, T = t^33), against 672 for B = 1 (T = t^45) and 720 for
   Part V. The searches cover T = t^L + x with 10 <= L <= 60 and deg x < 6. *32 points at point degree 2607*
   (B = t^4+t^2, T = t^86), against 3646 for B = 1 and 3968 for Part V.
3. The best B are exactly those divisible by **both degree-1 (inert) primes** t and t+1. This is the escape route
   that the VI.3 proof itself points to: B kills the small inert primes that carry Theorem E's loss.

Significance for the goal: **low (1/10)**. Both theorems are correct. Each is a statement about one rigid construction
class, a complete orthogonal core or a full hyperplane design, and the quadratic cost comes from that rigidity. Neither
bears on W(M) for unrestricted binary parents, and nothing transfers to Z.

---

## 1. VI.2: complete Hadamard cores

### 1.1 Setting, as I reconstruct it

There are m >= 4 marks (rows) and m-1 labels v (columns) with signs chi_v(i) = ±1. Adjoining the all-ones column
gives an m×m Hadamard matrix, so every column is balanced and distinct columns are orthogonal. For each label there
is a monic nonconstant f_v in F_4[t] with f_v ≠ fbar_v. The reduced marks are z_i = u_i·prod_v g_v^{(chi_v(i))}, where
g^{(+)} = f_v and g^{(-)} = fbar_v. All other prime material is common to the marks and has been divided out.
Write d_v = deg f_v and D = sum d_v. The cluster is critical: deg(z_i - z_j) <= D/2.

* *Units.* Since D > 0, deg(z_i - z_j) <= D/2 < D forces equal leading coefficients, so the u_i agree. Every product
  prod z_i^{e_i} with sum e_i = 0 lies in U_K, K = ceil(D/2), because it is a product of (z_i/z_0)^{e_i} with
  z_i/z_0 in U_K. Correct.
* *Common factors.* If the original points are G·z_i, every exponent vector with sum 0 kills G and the common unit.
  Removing G turns deg(Gz_i - Gz_j) <= (D+g)/2 into deg(z_i - z_j) <= D/2 - g/2. So the reduced cluster is still
  critical, and a lower bound on the reduced D bounds the original. Correct.

### 1.2 The exact products (CORRECT)

Over the full product, mark i contains f_v^{[chi_v(i)=+]} fbar_v^{[chi_v(i)=-]}.
* prod_i z_i^{chi_v(i)}: label v gives f_v^{m/2}·fbar_v^{-m/2} = r_v^{m/2}. A label w ≠ v contributes
  f_w^{n(++)-n(-+)}·fbar_w^{n(+-)-n(--)}, where n(ab) counts rows with (chi_w, chi_v) = (a,b). Balance and
  orthogonality give n(ab) = m/4 for all four sign pairs, which needs two distinct nonconstant columns, i.e. m >= 3.
  So 4 | m and these contributions are 1.
* prod_i z_i^{(chi_v - chi_w)/2}: the exponents are 0, ±1. Label v gives r_v^{m/4} and label w gives r_w^{-m/4}. A
  third label u needs #{chi_u=+, chi_v=+, chi_w=-} = #{chi_u=+, chi_v=-, chi_w=+}. This is the identity
  sum_i chi_u(chi_v - chi_w) = 0 restricted to the two blocks of m/4 rows. Astra calls this "orthogonality and
  balance"; it is exactly that, and it is a statement about pairs of columns, not triples.
* A cleaner way to see both at once: z_i = sqrt(N)·prod_w r_w^{chi_w(i)/2}, a formal identity, so
  prod z_i^{e_i} = N^{(sum e)/2}·prod_w r_w^{<e,chi_w>/2}.

*Checked.* `test_vi2.py` verifies both products as **polynomial identities** (cross-multiplied, no division) on the
Part V clusters (m = 4, 8, 16) and on my own VI.5 squarefree clusters (m = 4, 8, 16). All are True.

### 1.3 The depth calculus and the inequalities (CORRECT)

* *Depth.* For r = 1+x with ord_inf x = delta: r^2 = 1+x^2 has depth 2·delta. For odd b,
  (1+x)^b - 1 = x(b + O(x)) with b ≡ 1, so the depth is unchanged. Hence depth(r^{2^a b}) = 2^a·depth(r) exactly.
  The proof is correct, and the random test (300 pairs × 24 exponents, `test_vi3_prime.py`) gave 0 failures. In the
  clusters, depth(r_v^{m/2}) = (m/2)·depth(r_v) exactly: 8, 64, 384, 32, 256, 1536.
* *depth(r_v) <= d_v.* depth(r_v) = d_v - deg(f_v + fbar_v) = d_v - deg B_v, because f + fbar = B for
  f = A + theta B (theta + thetabar = 1). Correct, and asserted in the script.
* *Odd part u = 1.* Write m = 2^{a+1}u. Then depth(r_v) >= uD/m. Summing gives (m-1)uD/m <= sum d_v = D, so u = 1
  for m >= 3. A corollary that Astra does not state but which follows directly: **complete Hadamard cores of
  non-power-of-two order (Paley, order 12, 20, ...) never occur over F_2[t], at any degree.**
* *Bounds.* depth(r_v) >= D/m gives deg B_v <= s_v := d_v - D/m. Then s_v >= 0 and sum s_v = D/m. Since m/4 is a
  power of two, depth(r_v/r_w) >= 2D/m.
* *The numerator.* f_v fbar_w + fbar_v f_w = A_vB_w + A_wB_v (char 2, theta+thetabar = 1, theta·thetabar = 1). The
  script cross-checks this identity on every pair. The remark "before common cancellation" is correct and also
  harmless: ord_inf of a quotient is deg(denominator) - deg(numerator) whether or not common factors cancel. No
  coprimality of f_v with fbar_v is needed at this step. So deg(A_vB_w + A_wB_v) = d_v + d_w - depth(r_v/r_w)
  <= s_v + s_w. The numerator is nonzero because r_v ≠ r_w, which holds since f_v ≠ f_w.
* *Counting.* By Markov, at least (m-1)/2 labels have s_v <= s = 2D/(m(m-1)). Among them there are at most
  2^{floor(s)+1} - 1 nonzero B_v of degree <= s. Within a class B_v = B_w = B, the bound
  deg(A_v + A_w) + deg B <= 2s leaves at most 2^{2s+1} values of A_v. Distinct labels have distinct f_v, so
  (m-1)/2 <= 2^{3s+2}, which gives s >= (log_2(m-1) - 3)/3 and hence (VI.2). Correct. The bound is vacuous for
  m <= 8 (value -2.8 at m=4, -1.8 at m=8, 36.3 at m=16).

### 1.4 Exact hypotheses, a weakening, and what is not covered

**Hypotheses actually used.** (H1) The sign matrix with the all-ones column adjoined is Hadamard. "Normalized" means
only that the constant column is the one adjoined; column negation is free (swap f_v and fbar_v), but row negation
is not, since it would replace a point by its conjugate. (H2) m >= 4. (H3) The f_v are monic and nonconstant, with
f_v ≠ fbar_v and **f_v pairwise distinct**. (H4) The displayed block products z_i, of degree D = sum d_v, satisfy
the cluster condition at that D.

Astra's hypotheses, namely f_v coprime to fbar_v and pairwise disjoint base-prime supports, are a sufficient
condition for (H3)+(H4) to hold *after* removing the common factor. With them the block products have gcd 1, since
every column is nonconstant and the blocks are coprime. Otherwise disjointness is used only for distinctness.
**Strengthened statement (proved by the same argument):** (VI.2) holds whenever (H1)-(H4) hold. Part V satisfies
(H1)-(H4) undivided, at D = (m-1)L, although its supports are not disjoint (checked for m = 8, 16). It therefore
belongs to the same theorem. Without (H4) for the undivided products the argument needs care: if blocks share
primes, dividing by a common factor that cuts into the blocks breaks the block form.

**Duplicated or negated columns** reduce to (H1): merge the blocks f_v f_w, or swap f and fbar.

**What VI.2 does not cover.**
1. *Incomplete cores.* Suppose only k of the m-1 Hadamard columns carry blocks, with the restricted rows still
   distinct; for Walsh cores this means the k labels span F_2^r. The same proof gives sum s_v = D(1 - k/m) and
   s = 2D(1-k/m)/k, hence only **D >= k(log_2 k - 3)/(6(1 - k/m))**. For k <= m/2 this is order m log m, the growth
   scale. **The quadratic cost comes entirely from completeness**: the free depth D/m is shared by m-1 labels.
2. *Flipped or extended profiles* (rows with extra non-core primes, the record's "fully flipped" Paley profiles),
   and any profile whose columns are not orthogonal and balanced.
3. *Non-common division*, and labels sharing primes in a way that makes (H4) fail undivided.
4. Constant fields larger than F_4, and odd characteristic. The 2^a-rule is the characteristic-2 version. In
   characteristic p the p-power plays this role, and the "u = 1" step changes.

### 1.5 Numerical test (as requested): Part V and VI.5 at m = 4, 8, 16 (`test_vi2.out`)

| cluster | D | squarefree | disjoint supp. | deg B_v | s_v | min depth(r_v) vs D/m | tightest s_v+s_w - deg(A_vB_w+A_wB_v) | (VI.2) bound |
|---|---|---|---|---|---|---|---|---|
| Part V m=4 | 12 | no | yes | 0 | 1 | 4 vs 3 | 0 | -2.8 |
| Part V m=8 | 112 | no | **no** | 0 | 2 | 16 vs 14 | 1 | -1.8 |
| Part V m=16 | 720 | no | **no** | 0 | 3 | 48 vs 45 | 2 | 36.3 |
| VI.5 m=4 | 48 | yes | yes | 0 | 4 | 16 vs 12 | 0 | -2.8 |
| VI.5 m=8 | 448 | yes | yes | 0 | 8 | 64 vs 56 | 4 | -1.8 |
| VI.5 m=16 | 2880 | yes | yes | 0 | 12 | 192 vs 180 | 8 | 36.3 |

All inequalities hold. All six clusters have a normalized Walsh profile and an all-plus row. My CRT construction of
VI.5 reproduces Astra's moduli degrees 6, 26, 60 and parent degrees 96, 896, 5760, with gcd(N, N') = 1 checked. In
both families B_v = 1 for every label: f_v = Y + v, so A_v = T + v. The pair inequality then reads
deg(v+w) <= 2L/m. This is the binding mechanism, namely m-1 distinct polynomials with pairwise differences of small
degree, and it is the same additive cost as in Part V. Using B_v ≡ 1 the counting improves to
D >= m(m-1)(log_2(m-1) - 2)/4. Part V attains m(m-1)(log_2 m - 1), so within (H1)-(H4) the order is pinned to a
factor of about 4. Under Astra's disjoint-support hypothesis the best known example is VI.5, at
4m(m-1)(log_2 m - 1), about 24 times the bound (VI.2).

### 1.6 Novelty and significance of VI.2

The products prod z^{chi} are the function-field form of the record's "integral character" (item 276;
`overview.md` §3.1, which proves pure Hadamard cores **uniformly bounded over Z**). The new ingredients are the exact
depth multiplication by 2^a and the finite count of bounded-degree traces B_v. The consequences in the F_2 model are
new: no non-power-of-two cores, and Theta(m^2 log m) for complete cores together with VI.5. The record's contrast
with Z is instructive: over Z complete cores are bounded, over F_2[t] they are unbounded at quadratic cost. For the
goal the significance is low. The class is far from growth scale, and (§1.4.1) the argument says nothing once the
core is incomplete.

---

## 2. VI.3: full hyperplane designs with translation B

### 2.1 Setting

V ⊂ F_2[t] is a subspace with |V| = Q, B ∈ V \ {0}, T ∈ F_2[t], and Y = T + theta B, so conj(Y + v) = Y + v + B.
**All** s = Q/2 forms with l(B) = 1 are used, each with an arbitrary side S_l = l^{-1}(c_l). The marks are
F_l = prod_{v∈S_l}(Y+v), and F_l·Fbar_l = P_V(Y) for every l. Unit normalisation and any common division are allowed.

### 2.2 First half, (VI.3): line by line (CORRECT)

* *Content.* With c = gcd(T, V), every Y+v = c·(...), so all F_l share c^{Q/2}. The new data are again of this form,
  and passing to the full gcd only lowers D' and preserves criticality. Correct.
* *At most one fixed pair.* Each S_l contains exactly one of v, v+B. If v and w represent two fixed pairs, every
  used form has l(v+w) constant (0 or 1). The forms with l(B) = 1 form an affine hyperplane of V* that spans V*
  (Q >= 4), so v+w ∈ {0, B}. Correct.
* *Height.* Write Y + v = (T+v) + theta B. Its degree is max(deg(T+v), deg B), since theta ∉ F_2 prevents
  cancellation. gcd(Y+v, Y+v+B) = gcd(T+v, B). The reduced numerator a_v is coprime to its conjugate: a common
  prime would divide (T+v)/g and B/g, contradicting g = gcd. So a_v involves only split primes and deg a_v = h_v.
  Correct; Astra states this without proof, and I supply it here.
* *Level sets.* v_f(F/Fbar) = sum_j (2|S∩A_j| - |A_j|).
  * For j <= b = ord_p B, A_j is B-invariant, so the term is 0.
  * For j > b, if v, w ∈ A_j then f^j | v+w. Since v+w is real, p^j | v+w. So A_j is a coset of V ∩ p^j F_2[t],
    contained in K_f = V ∩ p^{b+1}F_2[t], and B ∉ K_f.
  * A form l_0 with l_0|K_f = 0 and l_0(B) = 1 exists and is used, because all forms are present. It is constant
    on A_{b+1}, which contains every deeper level, so it takes all or none of each A_j with one sign.
  * A mark l_1 with the opposite choice at v_* takes at most half (or at least half) of each affine A_j, since an
    affine functional on a coset is constant or balanced.
  * Deeper levels contain only labels from the (at most one) fixed pair, with the same contribution for every mark.
  * |v_f((Y+v)/(Y+v+B))| equals the number of levels j > b containing v or v+B, at most one of which ever lies in
    such a level.

  Range >= sum_{variable} |v_f(ratio_v)|. Correct.
* *Summation.* v_f(z_l/zbar_l) = 2v_f(z_l) - e_p ∈ [-e_p, e_p], and common division and units only shift it. So
  2e_p >= range. Multiply by deg f and sum over conjugate pairs: sum_p |v_f(a_v/abar_v)| deg f = deg a_v = h_v, by
  the coprimality above. Since D' >= sum_split e_p deg f, we get D' >= (1/2) sum_var h_v. **No criticality is used.**
  This is a statement about common gcds, and I tested it on non-critical instances too.
* *Second inequality.* For p | B, the set {v : p | T+v} is empty or a coset of V ∩ pF_2[t]. It is never all of V,
  by content one, so it has at most Q/2 elements and is B-invariant. That leaves at least Q/4 pairs with p ∤ T+v,
  at least Q/4 - 1 of them variable, each contributing ord_p(B)·deg p to deg B - deg gcd(T+v, B) <= h_v. Summing
  gives (Q/4-1) deg B <= sum_var h_v. Correct.

### 2.3 Second half (CORRECT, asymptotic)

* *Common unit.* s' >= s/3 by pigeonhole over F_4^*. (s' = s whenever deg T > max(deg V, deg B).)
* *Choice of e.* e is odd with s'/8 < 2^e <= s'/2. Case deg B >= 2^{e-2}: (VI.3) gives
  D' >= (s/4 - 1/2)·s/96 = Omega(s^2).
* *Otherwise.* Odd-degree primes of degree e not dividing B have total degree >= 2^e - O(2^{e/2}) - deg B >= 2^{e-1}
  for large e. Correct.
* *Theorem E with B in place of 1.*
  * P-special forms (vanishing on K_p) factor through ev_p(V) and send Bbar ≠ 0 to 1, so there are at most 2^{e-1}
    of them.
  * Two non-special forms attain (1,1) on K_p: equal restrictions or independent ones.
  * With k ∈ K_p, l_i(k) = l_j(k) = 1, we have B + k ∈ J and B + k ≡ B (mod p). So {v ∈ J : v ≡ B mod p} is a coset
    of J ∩ K_p of size >= |J|/2^e = Q/(4·2^e), and ord_p P_J(B) >= Q/(4·2^e).
  * P_J(B) ≠ 0 since B ∉ J.

  Correct.
* *Difference identity.* ker l_i = J ∪ (u+J) and ker l_j = J ∪ (u'+J). Since l_i(B) = l_j(B) = 1, u' ≡ u + B mod J,
  and additivity of P_J gives F_j = P_J(Y+a)(X_0 + P_J(B)). So F_i - F_j = P_J(B) P_J(Y+a). An odd-degree real
  prime with p ∤ B divides no A + theta B with B ≢ 0, so it divides neither X_0 nor gcd(X_0, P_J(B)). Hence
  delta_ij >= sum_{p odd, p∤B} ord_p(P_J(B)) deg p. Correct.
* *Margin identity.* For each split p, deg gcd(z_i, z_j) >= deg f·(e_p - |a_i - a_j|), and each layer separates at
  most s'^2/4 pairs. Inert primes are fully common. So sum_{i<j} deg gcd >= D'(C(s',2) - s'^2/4) and
  sum delta <= s'D'/4. Correct; it is the record's `ref5_extend.md` identity.
* *Conclusion.* D' >= (4/s')·C(3s'/4, 2)·Q/8 ~ (9/64)Qs' >= (3/32)s^2 asymptotically. In both cases
  D' = Omega(s^2), with an overall constant of about 1/400 and "s' sufficiently large" unspecified. **No threshold
  is given; the result is a pure Omega statement.**

### 2.4 Computational checks of VI.3 (`vi3_*.out`)

| run | instances | critical | (VI.3) first ineq. failures | (VI.3) second | identity F_i-F_j | inert bound | margin sum delta <= s'D'/4 |
|---|---|---|---|---|---|---|---|
| d=3, V=F_2[t]_{<3}, all 7 B, all T with deg T <= 9, all 16 sides | 114,688 | 54,640 | 0 | 0 | 0 fail | 0 fail | 0 fail (max ratio 0.83) |
| d=4, V=F_2[t]_{<4}, 10 B, 167 T, all 256 sides | 427,520 | 8,324 | 0 | 0 | 0 fail | 0 fail | 0 fail (max ratio 0.87) |
| d=4, 25 random V ⊂ F_2[t]_{<7}, random B, 12 T each, all sides | 76,800 | 0 | 0 | 0 | -- | -- | -- |
| per-prime range ≥ sum_var and 2e_p ≥ range, random d=3,4 instances | 29,557 primes | -- | 0 | -- | -- | -- | -- |

The largest observed ratio sum_var h_v / (2D') is 0.67 for d=3 and 0.55 for d=4, so the first inequality has real
slack. The per-prime range inequality is attained with equality at 488 primes. The "at most one fixed pair" claim
held on every instance.

### 2.5 New data: nonconstant translation beats B = 1 (exact searches, `exh4*.out`, `exh5_all_*.out`, `exh6*.out`)

`minsearch.py` finds the least full-gcd reduced degree D' of a c=0 cluster formed by all s = Q/2 marks, over all
side choices. It uses exact monotone pruning: adding a mark only lowers deg G and raises the maximum difference.
It reproduces the record's 104 for d=4, B=1 (T = t^15).

**d = 4 (8 points), all T with deg T <= 16:**

| B | 1 | t | t+1 | t^2 | t^2+1 | **t^2+t** | t^2+t+1 | t^3 | t^3+1 | t^3+t | t^3+t+1 | t^3+t^2 | t^3+t^2+1 | t^3+t^2+t | t^3+t^2+t+1 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| least D' | 104 | 86 | 86 | 84 | 84 | **68** | 108 | 104 | 104 | 80 | none | 80 | none | 104 | 104 |

The 68-cluster has Y = t^11 + theta(t^2+t), sides (0,1,0,1,1,0,1,0), gcd degree 20, and pair-difference degrees
29..34 <= 34 = D/2. N is not squarefree (deg gcd(N,N') = 4). No point's conjugate is in the set. Re-verified
independently by `verify68.py`. **This 8-point cluster on a parent of degree 136 improves the record's 208.**

**d = 5 (16 points), T = t^L + x, 10 <= L <= 60, deg x < 6, all 31 B:**
B = 1: 672 (T = t^45; Part V has 720). B = t^3+t and B = t^3+t^2: **472** (T = t^33). B = t^2+t: 476.
B = t, t+1: 574. B = t^2, t^2+1: 540. Every B with D' < 540 is divisible by **both** t and t+1. B containing the
split prime t^2+t+1, or an inert cubic, is *worse* than B = 1 (652-848).

**d = 6 (32 points), T = t^L + x, 60 <= L <= 130, x ∈ {0,1}, selected B:**

| B | 1 | t | t^2+t | t^3+t | t^3+t^2 | t^4+t^2 | t^4+t | t^4+t^3+t^2+t | t^5+t |
|---|---|---|---|---|---|---|---|---|---|
| least D' | 3646 | 3198 | 2716 | 2615 | 2647 | **2607** | 2800 | 2736 | 3168 |

Part V at m = 32 has 3968. The B = 1, t^3+t and t^4+t^2 clusters were re-verified from scratch (`verify_gen.py`):
32 distinct points with equal norms and maximum difference <= floor(D/2). The normalised cost D'/(m^2 log_2 m) of the
best design is 0.35, 0.46, 0.51 at m = 8, 16, 32, against 0.54, 0.66, 0.71 for B = 1. The gain is a constant factor
of about 0.7, and the growth looks like m^2 log m inside V = F_2[t]_{<d}, as predicted below.

**Reading.** Theorem E locates the B = 1 loss at small inert primes dividing P_J(1). With B divisible by t(t+1),
the degree-1 inert primes no longer divide P_J(B) on most pairs, and the cost drops by 30-35%. This is exactly the
dichotomy in the VI.3 proof: B either kills the small inert primes, or has small degree and pays at them. It is
also why VI.3 proves only Omega(s^2) and not Theorem E's s^2 log s.
* Inside V = F_2[t]_{<d}, deg B < d = log_2 Q, so B can kill only the inert primes of degree up to about
  log log s. The log factor should survive there.
* In a general V containing a B of degree about s that is divisible by all odd-degree primes up to log s, the
  inert obstruction disappears, and (VI.3) gives only D' >= (Q/8) deg B ~ s^2.

**So the true order of full designs with translation may be Theta(s^2) rather than Theta(s^2 log s).** I conjecture
this but have not proved it. It is still quadratic and gives no route to growth scale.

### 2.6 What VI.3 does not cover

1. **Sparse subfamilies.** The first half needs the form annihilating K_f to be present at every split prime, and
   "at most one fixed pair" needs all forms. The second half alone survives for subfamilies, and it is defeated
   exactly by B divisible by the small odd-degree primes, which §2.5 shows is the profitable choice anyway. A
   subfamily of s' marks with such a B is covered by neither half.
2. Non-hyperplane transversals (S_i \ S_j not an affine coset), with 2^{Q/2} candidate points of norm P_V(Y).
3. B ∉ V, which needs V + <B> and is then a subfamily of a larger design.
4. Independent orientation choices at individual split primes, non-common division, unions of designs on one
   parent, and deliberate prime powers.

---

## 3. What remains open for growth scale after VI.1-VI.3, and my view

**Closed in the F_2 model (all proved and checked):**
* High-rank parity profiles: VI.1 (not my block), for example Paley double cores.
* Complete Hadamard cores: VI.2, Theta(m^2 log m), and m a power of two.
* Full hyperplane designs with translation 1: Theorem E, Theta(s^2 log s).
* Full hyperplane designs with any translation: VI.3, Omega(s^2).

**Open (no argument in the record or Part VI touches them):**
1. *Incomplete orthogonal cores*, with k <= m/2 Walsh columns spanning the rows. The VI.2 method gives only
   D >~ k log k here (§1.4.1).
2. *Sparse subfamilies of hyperplane designs* with B divisible by all small odd-degree primes (§2.6.1).
3. *Non-hyperplane transversals* of V ∪ (V+B).
4. *Independent prime choices*: per-prime orientation patterns that are not columns of one design. VI.1 is the only
   general constraint, and it allows low-rate codes (Astra's x, x^3 double Walsh passes it at W = O(m log m)).
5. *Unions* of several designs, or of designs with different (V, B, T), on one parent.

**My view on W = M^{1+o(1)} over F_2[t].** I do not find it supported by anything in hand. Every concrete family
built so far, including the improved ones above, sits at W ≍ M^2·(log M)^{0 or 1}:
* Part V and VI.5 are additive (subspace-polynomial) designs.
* Each lower bound (Theorem E, VI.2, VI.3) shows that the quadratic cost is forced by the rigidity that makes those
  designs work: all forms present, or a complete orthogonal core. A growth-scale cluster would have to give up that
  rigidity, and with it the exact algebraic identities (F_i - F_j = P_J(B)P_J(Y+a), exact integral characters) that
  make the known clusters provable.

Equally, nothing excludes M^{1+o(1)}:
* The exhaustive range (M_2(D,0) = 4 for D <= 17) is far too small to discriminate.
* The only unconditional lower bound is still the transferred W >= 3M log_2 M - O(M), the same order as over Z.
* VI.1 is the only parent-wide constraint, and it is satisfied by low-rank profiles.

On balance I would put the probability of growth-scale clusters over F_2[t] at well under even. My reason is
heuristic: it is not a theorem, and the classes above are not known to be exhaustive. A random-model count also
gives no clusters beyond bounded size, so any cluster at all needs structure, and the only structure found
(additivity) costs M^2. Proving W(M) = omega(M log M) for unrestricted parents looks out of reach with the present
tools, since every proof in hand is a statement about one construction class. In any case none of VI.1-VI.3 uses a
step available over Z except the margin identity, which is already in the record. So **none of this changes the
integer growth rate**, and the integer significance is nil.

---

## 4. Side remarks

* **The record's d=3 table** (`ref5_lower.md` line 160: "s = 2, 3, 4 need point degree 1, 8, 12") is for untwisted
  T. With T = 0 (T ∈ V, unit-twisted marks) and B = 1, the 4 marks of V = F_2[t]_{<3} reduce, after a common factor
  of degree 1, to 4 distinct monic points of degree **4** with pair differences of degree <= 2 and a common norm
  N = t^8+t^6+t^5+t^3+t^2+t+1 (`verify_d3.py`). This agrees with M_2(4,0) = 4 in the exhaustive table. It is a minor
  inaccuracy in a table, and the d=4 claim ("twisted points add nothing beyond s=2") is untouched.
* Astra's sentence in VI.2, "These are conditions on the prime profile, with no additive hypothesis on the roots",
  is accurate. Nonetheless every known instance (Part V, VI.5) is additive with B_v ≡ 1.
* VI.3 defines "variable" relative to the chosen sides. With all forms present, at most one pair is fixed whatever
  the sides, so the qualifier costs at most one pair. Correct as stated.
