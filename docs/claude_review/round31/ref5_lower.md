# Referee report, round 31, Part V, block "lower": section V.2 (the Omega(s^2 log s) bound inside the construction)

Source: `astra_free_V/answer.md` lines 1192-1241 (section V.2) and summary lines 1251-1253, `v_fixed2_check.py`
(function `check_residual_gcd`). My scripts are in `ref5_lower/`. I wrote them from scratch. F_4[t] is coded as
coefficient lists, not as Astra's pairs of bit masks. Gcds are computed two ways: by Euclid, and by valuations after
factoring every Y+v over F_4.

| script | what it checks | output |
|---|---|---|
| `f4poly.py` | F_4[t] / F_2[t] arithmetic (tables built from theta^2=theta+1) | - |
| `check_formula.py` | norm identity; raw difference formula; deg P_J(1) bound; deficit sum bound; h <= Q/(2s) | `check_formula.out` |
| `enum_gcd.py d4` | Astra's d=4 enumeration (Euclid and valuations), extended to 64 classes of T | `enum_d4.out` |
| `enum_gcd.py d5 16` | **exhaustive** d=5: all 43,046,688 selections x 16 classes of T (deg T <= 8) | `enum_d5.out` |
| `enum_gcd.py d5euclid` | Euclid = valuation gcd on 1,600 random d=5 selections | `d5_euclid.out` |
| `family_min.py` | exact least point degree of a c=0 cluster inside the family, d=3,4 | `family_d4*.out` |
| `twist.py`, `twist_identity.py` | unit-twisted (monic-normalised) points at deg T < d; a referee identity | `twist*.out` |
| `fullparent.py` | every monic point on a small V.1 parent, and the largest c=0 cluster among them | `fullparent.out` |
| `family_fast.py`, `verify_s8.py` | exact family minimum at d=4 for all deg T <= 18; explicit 8-point cluster | `family_d4_fast.out`, `verify_s8.out` |

## 0. Verdict

**V.2 is correct as stated.** I checked (V.3), (V.4) and the final combination line by line and found no error.
Every numerical ingredient was reproduced by independent code. The d=4 enumeration is reproduced exactly: 26,176
cases, maximum residual gcd degree 4, by two independent gcd methods. At d=5 I ran an **exhaustive** enumeration,
not a sample: 688,747,008 cases, maximum residual gcd degree **19**. The unsimplified prime bound is 22 and Qd/4 is
40, so (V.4) holds with slack at least 21 in every case.

Three caveats concern the scope of the statement, not its correctness:

1. **The thresholds d >= 48, s >= 6 make the theorem purely asymptotic.** It says nothing for Q < 2^48, so nothing
   for any family member that could be computed. Since Q >= 2s, the Omega(s^2 log s) conclusion is a statement only
   about s >= 2^47 (§3).
2. **Unit twists at deg T < d are a missing case (GAP, repaired here).** For T in V the raw marks F_S are not monic:
   their leading coefficient is theta or theta^2, depending on the side. The record's model uses monic points (or
   allows multiplication by the norm-one units F_4^*), and the twisted pair differences are not covered by (V.3).
   I prove that the twisted difference is A*O, where O is the odd-index coefficient sum of P_J. Its degree is exactly
   deg P_J(1) - a_1 with a_1 in {1,2}. So (V.3) holds for twisted points with loss at most 4, and the theorem survives
   unchanged asymptotically (§2.6).
3. **It is a statement about one construction only.** It is not a lower bound for the model, for this parent, or
   even for this norm form (§4). It does not limit other additive spaces, non-linear transversals, other translation
   elements, prime-power exponents, or other choices of F_4-primes above the same F_2-primes of N.

A by-product: inside the family, V.1 is not optimal. At d=4 an 8-point cluster with point degree 104 (parent 208)
exists at Y = t^15 + theta with mixed sides. V.1 needs 112 (parent 224). I verified it explicitly (§3).

**Novelty:** new, but by definition. It is a negative result about Astra's own Part V construction. Nothing in the
record overlaps it; grep finds no subspace-polynomial construction before Part IV/V.

**Significance for the goal: low (1/10).** It closes off one route to growth scale for one construction in a
characteristic-2 function-field model. That model is already known to behave differently from characteristic 0
(Part IV: M(D,0) <= 24 over C(t)), and its transfer to Z is negative (filter F3). It has no consequence for the
integer bound (2/3+eps) log R/loglog R.

## 1. The statement checked

Fix V = F_2[t]_{<d}, Q = 2^d. Select s >= 2 distinct forms l with l(1)=1 and one affine side S_i of each. Take any
T in F_2[t], Y = T + theta, F_i = prod_{v in S_i}(Y+v). Let D_0 = deg F_i, g = deg gcd_i F_i, h = |cap_i S_i|,
B = max(d, deg T). Claims: (V.3) g >= Q(d/2-6); (V.4) g <= hB + Qd/4. For d >= 48 and s >= 6, any common division
giving a c=0 cluster leaves point degree >= Qsd/16.

## 2. Line-by-line check

### 2.1 Setup: CORRECT
* Norm. conj(Y+v) = Y+v+1, and S_i, S_i+1 partition V (because l(1)=1). So F_i * conj(F_i) = P_V(Y) is in F_2[t]
  and all F_i have the same degree D_0. Verified exactly (`check_formula.py`, d=3..6, seven values of T each,
  every side).
* Division. With points F_i/G, deg G = a, the condition deg(F_i-F_j) - a <= (D_0-a)/2 is equivalent to
  a >= 2deg(F_i-F_j) - D_0. Any valid a satisfies a <= g, so the point degree is >= D_0 - g. Correct.
* h <= Q/(2s). If the selected forms span rank r', then s <= 2^{r'-1}, and a consistent system has Q/2^{r'}
  solutions. Correct. No violation in 15,000 random selections for d=4..8.

### 2.2 The raw difference formula F_i - F_j = P_J(1) P_J(Y+a): CORRECT
Proof. Take a in S_i cap S_j (it exists because the forms are independent). Then F_i = P_{H_i}(Y+a). The
recurrence P_{J+<u>}(X) = P_J(X)^2 + P_J(u)P_J(X) gives the difference (P_J(u)+P_J(u'))P_J(Y+a) = P_J(u+u')P_J(Y+a).
Since u+u' = 1 mod J and P_J is J-periodic, this is P_J(1)P_J(Y+a). Verified as an exact polynomial identity for
d=3..6, several T including T in V, and every side combination: 24 to 240 side pairs per (d,T).

### 2.3 Degree counts in (V.3): CORRECT (both bounds are lossy, harmlessly)
* deg P_J(1) >= (Q/4)(d-4). Correct. The proof gives (Q/4)(d-4)+2. Exhaustively over all pairs of forms, d=3..9,
  the true minimum is (Q/4)(d-3)+1 (d=9: 769 against the claimed 640).
* sum_v (B - deg(Y+v)) <= 2Q. Correct, with the convention deg(theta+c) = 0. The exact maximum is 2^{d+1}-2 = 2Q-2,
  attained at T in V, verified exhaustively for all T of degree < d+2, d=2..10. For deg T >= d the sum is 0.
* deg P_J(Y+a) >= QB/4 - 2Q and D_0 <= QB/2. Then 2deg(F_i-F_j) - D_0 >= Q(d-4)/2 + QB/2 - 4Q - QB/2 = Q(d/2-6).
  The arithmetic checks.
* Remark (referee). For deg T >= d every factor has degree B, and the requirement is exactly
  g >= 2 max_pairs deg P_J(1). It is independent of B and of the sides, and it is >= (Q/2)(d-3)+2.

### 2.4 The gcd analysis (V.4): CORRECT
* Removing prod_{v in C}(Y+v) (degree <= hB) is legitimate, since that product divides every F_i.
* Even-degree contraction. If f | Y+v then theta = T(alpha)+v(alpha) lies in F_2(alpha), so [F_2(alpha):F_2] = e is
  even, and the minimal polynomial of alpha over F_4 has degree e/2. So P = f*fbar with f != fbar, and for
  w in F_2[t] we have f^j | w iff P^j | w. Correct.
* Affine structure. For j <= M := max_{v not in C} ord_f(Y+v), the set {v in V : f^j | Y+v} is the coset
  v_* + W_j, with W_j = V cap P^j F_2[t]. This has size 2^{d-ej} if ej < d and is {0} otherwise. Correct, and this
  is where V = F_2[t]_{<d} is used essentially.
* "At most half lie in S_i." S_i is a level set of l. On the coset, l is either constant, equal to its value at
  v_* (which is not S_i's value, so the intersection is empty), or balanced. Correct. Above M, no v in S_i \ C has
  f^j | Y+v. Elements of C with higher order were removed with C. Correct.
* Prime count. The number of elements of F_{2^e} of exact degree e is e*pi_2(e) <= 2^e. Each base prime contributes
  e*(Q/2) sum_{j: ej<d} 2^{-ej}, so the total is <= (Q/2) sum_{e even, 2<=e<d} 1/(1-2^{-e})
  <= (Q/2)((d-1)/2 + sum_k 1/(4^k-1)). Here sum_k 1/(4^k-1) = 0.4207 <= 4/9, so the total is <= Qd/4. Correct.

### 2.5 Combination: CORRECT
hB >= Q(d/4-6), together with h <= Q/(2s), gives B >= s(d/2-12) >= sd/4 for d >= 48. If h = 0 the inequality
already fails for d > 24, so no cluster exists. For s >= 6 (s >= 5 would do), B > d, so D_0 = QB/2 and the factors
over C have degree exactly hB. Then D_0 - g >= (Q/2-h)B - Qd/4 >= (Q/2)(1-1/s)(sd/4) - Qd/4 = Qd(s-3)/8 >= Qsd/16
for s >= 6. Since Q >= 2s, the parent degree is >= s^2(1+log_2 s)/4. "Sharp order" against V.1 (point degree about
Q^2 d/4 at s = Q/2): up to a factor of about 8. All correct.

### 2.6 Missing case: unit twists when deg T < d (GAP; repaired by the referee)
If T is in V (WLOG T = 0), the factors Y+T = theta and Y+T+1 = theta^2 are units, and each side contains exactly one
of them. So lc(F_S) is theta or theta^2, half the sides each (`twist.out`). Points in the model are monic. In any
case u*z with u in F_4^* has the same norm, since N(u)=1. The legitimate points are therefore u_i F_i/G, and (V.3)
was proved only for u_i all equal. For deg T >= d all F_i are monic, unequal units give differences of full degree,
and nothing is lost there.

**Referee identity (proved, and checked exactly on 114 pairs for d=3..6).** Let 0 be in S_i and 1 in S_j (distinct
forms). Write A = prod_{v in S_i cap S_j}(theta+v) and P_J(X) = sum_k c_k X^{2^k}. Then
F_i/theta - F_j/theta^2 = A * O, where O = sum_{k odd} c_k.
Proof. Split theta^{2^k} into theta (k even) and theta^2 (k odd). This gives P_J(theta) = theta E + theta^2 O and
P_J(theta^2) = theta^2 E + theta O. Then theta P_J(theta) + P_J(theta^2) = (1+theta) O = theta^2 O.
Newton polygon. Take a basis of J with distinct leading degrees a_1 < ... < a_{d-2}. Removing the 2^k smallest roots
is a strict break, so deg c_k = deg c_0 - sum_{i<=k} 2^{i-1} a_i, without cancellation. Hence
deg O = deg c_1 = deg P_J(1) - a_1, with 1 <= a_1 <= 2 (J meets span(1,t,t^2)).
Exhaustively for d <= 8, max(deg P_J(1) - deg O) = 2. So the twisted version of (V.3) is g >= Q(d/2-6) - 4, and the
conclusion of V.2 is unchanged: the thresholds move by O(1/Q).

Data: at T=0, among monic-normalised points with any sides (both sides of one form allowed), the largest c=0 cluster
has 2 points for d=4,5,6. For d=3 it has 4 points at point degree 4, better than the untwisted family (12). So twists
do matter at tiny d.

### 2.7 Both sides of one form (not in the statement; covered by the proof)
z_H and conj(z_H) differ by P_H(1), a constant in Y. This is the closest possible pair, but C = empty for any
selection containing both. (V.4) with h = 0 gives g <= Qd/4 < Q(d/2-6) for d > 24. So any selection with two
distinct forms and an opposite pair is not a cluster. The statement may be widened to "any set of sides".

## 3. Computations

**d=4 (reproduces Astra exactly).** There are 8 forms and 3^8-1-16 = 6,544 selections with s >= 2, for each of the
4 classes T in {0, t^4, t^5, t^4+t^5} mod V. Total **26,176 cases, max residual 4**, by Euclid and by valuations,
identical case by case. The histogram of residuals is the same for every class: {0:1574, 1:1296, 2:3072, 3:216,
4:386}. The unsimplified prime bound is exactly 4 (only P = t^2+t+1, j=1). Extended to all 64 classes with
deg T < 10: max residual still 4. The minimum slack in g <= hB + Qd/4 is 12. Astra's justification ("translation by
V permutes the sides") is correct.

**d=5 (exhaustive, not sampled).** There are 16 forms, giving 43,046,688 selections per class, for 16 classes of T
mod V (deg T <= 8, including T=0). The method is a meet-in-the-middle valuation computation. Euclid and valuation gcds
agree on 1,600 random selections.

| T class | primes | max residual | min slack in (V.4) |
|---|---|---|---|
| 0, t^5, t^6, t^6+t^5 | 22-42 | 18 | 22 |
| deg 7 (4 classes) | 44-46 | 19 | 21 |
| deg 8 (8 classes) | 40-54 | 18 or 19 | 21 or 22 |

Overall: **max residual 19 <= 22 (unsimplified bound) <= 40 = Qd/4.** The residual depends on T only modulo the small
split primes, so the histograms repeat with period 2^3 in the class index, as the proof predicts.

**Least point degree of an actual c=0 cluster inside the family (exact, `family_min.py`).** These d are far below the
thresholds; the data illustrate the family, not the theorem.

d=3 (all T classes with deg T <= 12): s = 2, 3, 4 need point degree 1, 8, 12 (12 = V.1 for m=4).

d=4. The run covers T = 0 (`family_min.py`, `family_d4.out`) and all 32,767 classes with 4 <= deg T <= 18
(`family_fast.py`, `family_d4_fast.out`). The two programs agree on deg T <= 12. For deg T > 18 the minimum cannot
improve: the point degree is >= (Q/2-h)B - 4, because the residual is <= 4 at d=4 (proved bound). Twisted points add
nothing beyond s=2 (`twist.out`). So these are the exact family minima at d=4:

| s | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| least point degree | 10 (T=0) | 32 | 32 | 80 | 80 | 96 | **104** |

The s=8 entry is **better than V.1** (point degree 112, parent 224). Take Y = t^15 + theta with sides
(f0 ker, f1 ker, f2 co, f3 co, f4 ker, f5 ker, f6 co, f7 co), where h = 1 and the full gcd has degree 16 = 15 + 1
(one residual linear prime). This gives 8 points of degree 104 on a parent of degree 208, with maximum pair difference
exactly 52 = D/2. I verified it explicitly: Euclidean gcd, all norms, all 28 differences (`verify_s8.py`).
This does not contradict V.1, whose "least L" claim is for the kernel-only marks with Y = t^L + theta. It shows that
V.1 is not optimal inside its own family, by a constant factor; the theorem's order s^2 log s is unaffected.

**Full point set of a small parent.** `fullparent.py` lists every monic point on the V.1 parent at d=3, L=4..8. That
is 16 to 144 points, of which 8 are family points. The largest c=0 cluster is 4 in every case, and it consists of
family points. So at this tiny size the extra points on the parent do not help. This is not evidence for large d.

## 4. Exactly what V.2 allows and what it does not

**Covered (with the referee additions of §2.6-2.7):** the fixed V = F_2[t]_{<d}; any T in F_2[t], which is the same
as any Y in F_4[t] with conj(Y) = Y+1; any subfamily of the 2^{d-1} forms with l(1)=1; any side of each (also both
sides); multiplication by units of F_4^*; division of all points by one common polynomial (the full gcd is optimal);
c = 0. Fixed c also works: replace (V.3) by g >= Q(d/2-6) - 2c, harmless for c = o(Q).

**Not covered. None of these is excluded by V.2:**
1. **Other additive spaces V**: any d-dimensional F_2-subspace of F_2[t] containing 1 that is not a degree
   truncation. Both ingredients fail. The count of low-degree elements (deg P_J(1), the deficit sum) and the fibre
   sizes |V cap P^j F_2[t]| = Q/2^{ej} are specific to F_2[t]_{<d}; fibres can be larger. Astra flags this herself.
2. **Other translation elements.** Y = A + B*theta with B a non-constant F_2[t] polynomial has conj(Y) = Y + B. One
   then needs B in V and forms with l(B) = 1. This cannot be reduced to B = 1.
3. **Non-linear transversals.** Every S with S and S+1 partitioning V gives a point prod_{v in S}(Y+v) with the same
   norm P_V(Y). That is 2^{Q/2} points, of which V.2 treats only the Q affine half-spaces.
4. **Other F_4-primes above the same F_2-primes.** Y+v is usually reducible over F_4, and a point on N may take f or
   fbar independently at each split prime and each multiplicity. The family uses one specific choice.
5. **Non-common division**, e.g. F_i/G for some i and F_j/conj(G) for others (conjugating some points after a common
   division), or general G_i with G_i*conj(G_i) independent of i. Here the cluster condition changes, so this is a
   genuinely different configuration.
6. **Prime powers and extra factors**: exponents n_v on Y+v, or extra factors w_i with w_i*conj(w_i) independent of i.
7. **Mixing evaluations or spaces** on one parent, larger constant fields F_{2^k}, odd characteristic, other
   anisotropic forms.
8. **Small d.** Because of the thresholds d >= 48, s >= 6, nothing is proved for Q < 2^48. The exact table in §3
   is the only information at small d.

**So V.2 must not be read as a lower bound for the model M_2(D,0) or for its growth scale.** The question of whether
O(m log m) parent degree is reachable over F_2(t) stays open.
