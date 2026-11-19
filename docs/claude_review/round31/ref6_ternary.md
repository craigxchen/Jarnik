# Referee report, round 31, Part VI, block "ternary" (VI.4: unbounded clusters for A^2+B^2 over F_3(t))

Source: `astra_free_VI/answer.md` lines 1352-1410, `final6.txt` items 7-9. I did not run Astra's
`vi_ternary_check.py`; I only read it. All scripts are mine and are in `ref6_ternary/`. My own arithmetic is
`f9lib.py` (F_9 = F_3[i], F_9[t] as numpy pairs with exact FFT convolution mod 3, Euclid, gcd) and `bigfield.py`
(F_729 = F_{q^2} by log tables). python-flint was used only as a gcd/factoring engine in three cross-checks
(`valgcd.py`, `valgcd32.py`, `factorN.py`, `extragcd.py`). Each of these is marked where it is used.

| script | purpose | output |
|---|---|---|
| `descent_check.py` | q=27, k=2, four marks at X=t^11: products of the 729 factors (t^11 - v)^{e_j(v)} over F_729. Checks F_9 descent, (VI.5), norms, full gcd, D, differences | `descent_check_out.txt` |
| `descent_k3.py` | the same direct-product check for k=3 (19,683 factors), Astra's five forms, L=3 | (console, recorded in §1) |
| `marks.py`, `test_marks.py` | raw marks for any k from subspace polynomials (Frobenius for q-th powers) and (VI.5); checked against the direct products | |
| `formal.py`, `subsets.py`, `q243.py` | exponent tables e_j(v) = 1+chi(eps_j lambda_j(v)), formal g, best form subsets and signs | `formal_out.txt` |
| `topcoef.py` | delta_ij = deg_t(c_i^{13} - c_j^{13}), the top X-coefficient of F_i - F_j | `topcoef_out.txt` |
| `five.py` | exact t-degrees of (F_i - F_j)(t,t^L) at degree ~10^6-10^7, with certification by the formal G only | `five_*.txt` |
| `fullgcd.py`, `fullgcd3.py` | exact degree of the FULL gcd (formula of §2.3), checked against actual gcds | |
| `fullscan.py`, `valgcd32.py` | least L with the full gcd; flint gcd and norm cross-check of the 5-point record | `fullscan_*.txt`, `valgcd32.txt` |
| `smallsearch.py` | all marks (both signs), all subsets, X = t^L, full subset gcd, exact cluster test (q=3, k=2,3; q=27, k=2) | `ss_*.txt` |
| `factorN.py`, `extragcd.py`, `multN.py` | prime structure of N; the extra gcd; multiplicities | (console, §3) |
| `p7.py`, `halfsys.py` | p = 7 instance; half-system lemma for q = 1 mod 4 | (console, §4) |

## 0. Verdict

* **VI.4 (theorem): CORRECT.** I checked every step by hand (§1) and reproduced it with my own arithmetic. For
  every m there are m distinct monic z_j in F_9[t] of one degree D, with one norm N = z_j z̄_j in F_3[t] and
  deg(z_i - z_j) <= D/2. Here deg N = exp(O(m log m)). So **sup_D M_3(D,0) = infinity** for the anisotropic form
  A^2+B^2 of F_9(t)/F_3(t), which is the reflection-type conjugation of Z[i]. The proof works verbatim for every
  p = 3 mod 4.
* **Four-mark computation: CORRECT, reproduced exactly** (gcd 1,038, D = 6,981, parent 13,962, differences
  3471, 3457, 3457, 3471, 3471, 3457). It is also the least D in the q=27, k=2, X=t^L family.
* **Five-mark certificate: CORRECT, and improvable by a factor of 2.2.** My exact expansion at L = 67 gives
  g = 1301, h = 729, D = 1,231,594 and maximum reduced difference exactly 615,589 < 615,797. With L = 66 it fails.
  A better choice of forms and signs, together with the full gcd, gives **5 points at D = 568,386** (checked
  independently, §2.4).
* **New computed instances:** **6 points at D = 3,582,084** (full gcd). With the formal G only, the certificate
  is rigorous at D = 6,880,250 (§2.5). q = 27 cannot give 7 points. q = 243 gives 7 and 8 at the formal level
  (§2.6).
* **Corrections and additions (none affects the theorem).** (1) The full gcd is much larger than the formal G
  when k >= 3 (about 17,500 extra, all at places in F_{q^2}), so Astra's L is far from minimal. (2) N is far
  from squarefree. Most primes occur squared, and for k >= 3 the F_{q^2}-rational places carry multiplicity
  about q^{k-2} (27-31 in the 5-point example), so multiplicities are **not** bounded as m grows. (3) Most of
  deg N sits on primes whose residue fields contain F_{q^2}. These fields grow with m (§3).
* **Exact rate (my statement, §5):** M_3(D,0) >= (log_2 3 - o(1)) log D / log log D for **all** large D. Astra's
  L = 2(k-1)n+1 gives the constant (1/2)log_2 3. Taking L ≈ (k-1)(n-h)/(g-h) instead doubles it. That is
  m ≍ log W / log log W with W = deg N, doubly logarithmic in the analogue of R.
* **Novelty: high relative to the record.** It answers the record's open question 1 of `README_V.md` §5 ("is
  sup_D M_3(D,0) finite?"). It also bypasses Lemma P: the cofactor differences are nonconstant, and the
  exponents 0/2 come from quadratic characters.
* **Significance: function-field model, substantial. Integer goal, low (2/10).** The exemption "arguments
  needing odd characteristic or a reflection-type conjugation" in the expert note's F3+ is dead. What is left
  is "squarefree or bounded-multiplicity N" (and characteristic 0). The mechanism uses q-linearised subspace
  polynomials, Frobenius descent and prime squares, none of which has an integer counterpart. Nothing changes
  for Z.

## 1. Line-by-line check of VI.4 (all CORRECT)

Notation: q = 3^r with r odd, F_{q^2} ⊃ F_9 ⊃ F_3, i^2 = -1. V = {i Σ_{b<k} x_b t^b : x in F_q^k}, n = q^k,
h = n/q. σ is the coefficient map α ↦ α^9 and τ is α ↦ α^q.

1. **V is an F_q-space stable under σ and τ.** σ(i x_b) = i x_b^9 because i^9 = i, so σ acts on V as x ↦ x^9.
   τ(i x_b) = i^q x_b = -i x_b because q = 3 mod 4 (r odd), so τ(v) = -v. ✓
2. **Descent to F_9[t,X].** e_j(σv) = 1 + chi(eps_j lambda_j(x^9)) = 1 + chi((eps_j lambda_j(x))^9) = e_j(v).
   This holds because lambda_j has F_3 coefficients, eps_j^9 = eps_j, and chi(y^9) = chi(y)^9 = chi(y). So
   σ(F_j) = F_j. Fix(σ) on F_{q^2} = F_{9^r} is F_9. ✓ Note that F_q ∩ F_9 = F_3 for r odd, so F_3-coefficients
   are forced: a form with coefficients outside F_3 is either not F_q-linear or not σ-stable.
3. **Conjugation equals τ on F_9.** α^q = α^{3^r} = α^3 on F_9 for r odd. ✓ So z̄ = τ(z) for z in F_9[t,X], and
   τ(F_j) = Π(X+v)^{e_j(v)} = Π(X-v)^{e_j(-v)}. chi(-1) = (-1)^{(q-1)/2} = -1 because (q-1)/2 is odd. So
   e_j(v) + e_j(-v) = 2 for every v, including the kernel (1+1). Hence F_j F̄_j = P_V^2. ✓ P_V is σ- and
   τ-invariant, so P_V lies in F_3[t,X].
4. **Distinctness.** The F_q-kernels of distinct F_3-lines are distinct. At v in H_i \ H_j, e_i = 1 and
   e_j ∈ {0,2}. Unique factorisation in F_{q^2}(t)[X] does the rest. ✓
5. **Linearity and (VI.5).** For P_{J+F_q u} the recursion is standard (Π_a(Z - aw) = Z^q - w^{q-1}Z with
   Z = P_J(X), w = P_J(u)). On V, P_H is F_q-linear with kernel H, so P_H(v) = c·eps lambda(v). Grouping v by
   a = eps lambda(v) gives F_j = Π_a(A - ac)^{1+chi(a)} = A·(Π_{s square}(A - sc))^2 =
   A(A^{(q-1)/2} - c^{(q-1)/2})^2. Expanding gives P_V + 2c^{q-1}A - 2c^{(q-1)/2}A^{(q+1)/2}. ✓
   **COMPUTED:** the closed form equals the direct product over F_729 for the four marks (k=2, L=11) and for
   Astra's five marks (k=3, L=3, products of 19,683 factors). All coefficients lie in F_9 and the norms are
   equal and lie in F_3[t] (`descent_check.py`, `descent_k3.py`). A side fact used below: P_H has F_3[t]
   coefficients and c ∈ i·F_3[t].
6. **Degree of differences.** The top X-term of F_j - P_V is c^{(q-1)/2}A^{(q+1)/2}, with X-degree
   h(q+1)/2 = (n+h)/2 (an integer). So deg_X(F_i - F_j) <= (n+h)/2. ✓ The leading X-coefficient of F_i - F_j is
   c_i^{(q-1)/2} - c_j^{(q-1)/2}. Its t-degree δ_ij is about (k-1)(n-h)/2 (`topcoef_out.txt`: 337-351 for k=2,
   18,562-18,954 for k=3). This is what fixes the least L.
7. **Common factor.** G = Π(X-v)^{min_j e_j(v)} is σ-invariant, so it lies in F_9[t,X]. G Ḡ is τ-invariant, so
   P_V^2/(G Ḡ) lies in F_3[t,X]. ✓
8. **Expectation (VI.6).** If all m values at v are nonzero, then P(min = 2) = 2^{-m}, since chi(eps·y) is a
   uniform sign independent across j. If a < m values are nonzero, the min is 1 with probability 2^{-a} and 0
   otherwise. So E[g] >= 2^{1-m}n, and g > h iff q > 2^{m-1}. The q-condition is exact for this argument,
   which counts as "strict slack (g-h)/2 >= 1/2". It is **sufficient, not necessary**: q = 27 < 2^5 already
   gives g = 833 > 729 for m = 6 (§2.5). ✓ Also k = ⌈log_3(2m+1)⌉ gives (3^k - 1)/2 >= m, and
   q <= 9·2^{m-1}. ✓
9. **Specialisation X = t^L.** The X^a-coefficient of Z_j has t-degree <= (k-1)(n-g-a). For L > k-1 the points
   are monic of degree (n-g)L. For L > (k-1)(n-g) the intervals [aL, aL + (k-1)(n-g-a)] are disjoint, so
   nonzero differences stay nonzero. The difference bound ((n+h)/2-g)L + (k-1)n < (n-g)L/2 holds for
   L = 2(k-1)n+1, since (g-h)L/2 >= (k-1)n + 1/2. ✓ The sharper sufficient condition used in the five-mark
   certificate is also correct: deg <= a_0 L + (k-1)(n-h)/2 with a_0 = (n+h)/2 - g, which suffices for
   L >= (k-1)(n-h)/(g-h). Caveat: below L ≈ (k-1)(n-g), distinctness after substitution is **not** covered by the
   interval argument. In all my instances it was verified directly. Substitution commutes with coefficient
   conjugation, so the norm identity is exact. ✓
10. **Count.** log D <= log(nL) = O(k log q) = O(m log m). ✓

## 2. Reproductions and smallest instances (COMPUTED)

### 2.1 Four marks (q=27, k=2, forms x, y, x+y, x-y, eps = +, X = t^11)
Formal g = 92 and h = 27 (my enumeration). The raw marks have degree 8019. The **full gcd has degree 1,038 =
92·11 + 26**: it is G(t,t^11) (which I computed directly, F_9 coefficients) times **t^26** (flint factorisation of
the quotient). The extra comes from the place t: the fibre {x_0 = 0} = ker x, where t^11 - i x_1 t has
ord_t = 1. Reduced D = 6,981 and parent degree 13,962. The four norms are equal and lie in F_3[t]. The
differences are **3471, 3457, 3457, 3471, 3471, 3457** < 3490.5. All leading coefficients are 1, no point is the
conjugate of another, and the gcd of the reduced points is 1. Identical to Astra.
**Minimality inside the family** (`ss_27_2.txt`: all 8 marks with both signs, all subsets, full subset gcd,
L = 1..12): 4 points first appear at L = 11, D = 6,981. 3 points need D = 2,704 (L = 5).

### 2.2 Small m against the exhaustive table
q = 3, k = 2 (`ss_3_2.txt`, L <= 40) gives at most 3 points, the best at D = 4. q = 3, k = 3 (L <= 17) gives
only 2 points. The formal reason (`subsets.py`) is that for q = 3 one has g <= h as soon as m >= 3. For m <= 4
the construction is irrelevant: the record already has M_3(2,0) = 4. There is no conflict with
M_3(D,0) <= 4 for D <= 11 (`ref5_extend.md`). The first 5-point instance of this mechanism is five orders of
magnitude above the exhaustive range.

### 2.3 The full gcd for k = 3 (new, proved and checked)
Common roots of t^L - v and t^L - v' are roots of v - v' ∈ V\0, an F_q-polynomial of degree <= k-1 = 2 times i.
So they lie in F_{q^2} (only α = 0 among F_q, since α^L ∈ F_q ∩ iF_q = 0). Therefore
gcd = G(t,t^L)·Π_{α∈F_729}(t-α)^{x_α}, where
x_α = min_j Σ_{v∈fib(α)} e_j(v) ord_α - Σ_{v∈fib(α)} min_j e_j(v) ord_α, and fib(α) = {v : v(α) = α^L}.
For k = 3 each α ∈ F_729\F_27 has a fibre of q^{k-2} = 27 elements, and Σ_fib e_j = 27 whenever lambda_j is not
constant on the fibre. That gives about 702·25 ≈ 17,500 of extra gcd, independent of L. For k = 2 the fibres are
single points and only t contributes (the 26 above). `fullgcd3.py` matches actual gcds: 21548, 22442, 21845,
22858 (L = 3, 4, two form sets) and 61,470 (best5, L = 32). The last one came from flint's half-gcd on degree
~630k.

### 2.4 Five points
* Astra's set (x, y, z, x+y, x+z, all eps = +): **g = 1301** (my enumeration; all 32 sign patterns give at most
  1301). Exact expansion at L = 67 gives D = 1,231,594 and reduced differences 615238, 615589, 615224, 615197,
  615589, 615238, 615238, 615589, 615589, 615224. The maximum **615,589 equals Astra's bound exactly** and is
  <= 615,797. At L = 66 the formal certificate fails (606,684 > 606,606). **CORRECT.** With the full gcd, Astra's
  set already works at L = 37 (D = 662,632). L <= 36 fails (`fullscan_astra5b.txt`).
* **Best found: forms x_2, x_1, x_1+x_2, x_0, x_0+x_1-x_2 (coordinates on 1, t, t^2), signs (+,+,-,+,+): g = 1379
  (the maximum over all 5-subsets and signs, `subsets.py`).**
  * With the formal G only, L = 59 is certified: D = 1,079,936, max difference 539,747 <= 539,968. L = 58 fails.
  * With the full gcd, L = 32 works: **D = 568,386, parent degree 1,136,772.** The full gcd is 61,470 = 44,128 +
    17,342. The ten differences are 284076 (×4), 283725 (×5) and 283711, all <= 284,193. L = 26-31 fail
    (`fullscan_best5*.txt`). An independent check (flint gcd, exact division, norms) gives the same D and
    differences, all five norms equal, and the norm in F_3[t] (`valgcd32.txt`).
  * Since δ_max = 18,954 and the extra gcd is at most about 18k for every choice, L ≈ (2δ - extra)/(g-h) >= 31.
    So about 5.5·10^5 is the floor for 5 points in this family (q = 27, k = 3, X = t^L). I did not try other
    T(t) or other V.

### 2.5 Six points (new)
Forms x_2, x_1, x_1+x_2, x_0, x_0+x_1-x_2, x_0-x_1+x_2, signs (+,+,-,+,+,+): g = 833 > h = 729, although
q = 27 < 2^5.
* Rigorous with the formal G only (exact expansion at degree 7.2M): L = 365, **D = 6,880,250**, max difference
  3,440,099 <= 3,440,125. L = 364 fails.
* With the full gcd (formula §2.3): L = 191, **D = 3,582,084**, max difference 1,790,931 <= 1,791,042. L = 186-190
  fail. The gcd formula is validated, but an independent flint gcd at degree 3.8M ran out of memory. So this
  number rests on `fullgcd3.py`, while the 6.88M certificate does not.
* q = 27 cannot give 7 points: the best g is 651 < 729, and the extra gcd (≈18k) cannot cover 2δ ≈ 37.9k.

### 2.6 q = 243 (formal only)
k = 3 gives m = 7 with g = 232,684 and m = 8 with g = 116,403 (h = 59,049, all eps = +, the first 7 or 8 forms).
The sufficient L is 165 and 499, so D is about 2.3·10^9 and 7.1·10^9. Not expanded. Distinctness at these
L < q^2 is not covered by item 9 and was not checked.

## 3. Legitimacy as an instance of the record's model (CORRECT, with remarks)

* **Form and infinity.** F_9(t)/F_3(t) is the constant extension. Its degree-1 place at infinity is inert
  (residue F_3 ⊂ F_9), so the completion F_9((1/t)) / F_3((1/t)) plays the role of C/R. A^2+B^2 is anisotropic
  over F_3((1/t)) because -1 is a nonsquare in F_3 and the leading terms cannot cancel. This is exactly the
  setting of `ff.md` and of the M_3 table. ✓
* **Points.** They are monic of equal degree, so they are not unit multiples (a unit u ∈ F_9^* with
  u·monic = monic forces u = 1). Conjugates would count as distinct lattice points anyway (z and z̄ are (A,B) and
  (A,-B)), but none occurs in the four-mark example. ✓
* **Full gcd.** Dividing all points by a common monic d of degree e turns the condition into
  deg - e <= (D-e)/2, i.e. deg <= D/2 - e/2. This is stronger (`ff.md` N2), so the cluster condition only
  tightens. It never helps to cancel less. For the four-mark example the full gcd is G(t,t^11)·t^26. ✓
* **Inert primes.** At an inert p (odd degree), v_p(z) = v_p(N)/2 for every point, so p lies in the gcd. After
  full division, N has only split primes. In the four-mark example (flint factorisation) **all 353 primes of N
  have even degree.** ✓
* **Structure of N (four-mark example, deg 13,962).**
  * Prime degrees: 2, 4, 6, 10, 12, 18, 30, 36, 42, 48, 60. Multiplicity is 2 for 338 primes and 1 for 15; the
    radical has degree 7,124.
  * 13,860 of 13,962 (99.3%) of deg N lies on primes whose degree is divisible by 6. Their residue fields
    contain F_729 = F_{q^2}.
  * In general the primes of N are those below the irreducible factors of t^L - v(t), v ∈ V. Roots in F_{3^d}
    can account for at most about 3^d q^{k-1} of the 2nL roots. So the residue fields carrying the bulk of N have
    size >> qL > 2^{m-1}L, and they grow with m (this is a count, not a theorem about every prime).
  * The constant field is fixed (F_9/F_3). The auxiliary F_q only labels the roots.
* **Multiplicities.** The generic exponent pattern is (π^2, π̄^0) or (π π̄): exponents 1+chi ∈ {0,1,2}. For
  k >= 3 the F_{q^2}-rational places have larger multiplicity. In the 5-point example (L = 32) a point has
  multiplicity 0, 4, 27 or 31 at the 703 places α ∈ F_729 with nonempty fibre (577, 50, 74 and 2 places
  respectively; `multN.py`). This grows like q^{k-2}. So N is **far from squarefree, and its multiplicities
  grow with m.**
* **Squarefree variant.** It is not provided, and I see no route from this identity. The degree drop in (VI.5)
  is the sparsity of Π_{squares}(Y - s) = Y^{(q-1)/2} - 1, which needs the doubled exponents. The natural {0,1}
  version takes half-systems on V, giving F_j = (A^{(q-1)/2} - c^{(q-1)/2})·B_j with B_j a half-system product on
  H_j. Its pair differences drop only O(h) below the top degree (n-1)/2, against the ~n/4 needed. A
  bounded-multiplicity (<= 2) variant holds for k = 2 (the four-mark N), but k = 2 caps m at 4.
  **Open, plausibility low.**

## 4. Other characteristics

* **p = 3 mod 4: CORRECT.** Every step holds with F_p in place of F_3 and σ = p^2-Frobenius. chi(-1) = -1 for
  q = p^r with r odd, and there are (p^k - 1)/(p - 1) forms. **COMPUTED** (`p7.py`, pure Python over F_49): q = 7,
  k = 2, X = t^5, forms x, x+y, 2x+y with signs (+,+,-) give 3 points at D = 174. All marks have equal norms in
  F_7[t].
* **p = 1 mod 4: the mechanism does not transfer (formal obstruction, proved).** The only constant-coefficient
  anisotropic binary forms are similar to the norm A^2 - dB^2 of F_{p^2}/F_p. For r odd the q-power map still
  induces √d ↦ -√d. But q = 1 mod 4, so chi(-1) = +1 and e(v) + e(-v) ≠ 2. Any replacement needs an odd sign
  function s on F_q^* (a half-system S), and then deg_X(F - P_V) = h(1 + (q-1)/2 + δ_S) with
  δ_S = deg(Π_S(Y-a) - Π_S(Y+a)).
  * **Lemma (proved here).** For q = 1 mod 4, every half-system has δ_S >= (q-1)/4. Proof: write f = E + O
    (even and odd parts). Then E^2 - O^2 = Π_S(Y^2 - a^2) = Y^{q-1} - 1. If deg O < (q-1)/4, then
    E = Y^{(q-1)/2} + e and 2eY^{(q-1)/2} + e^2 = O^2 - 1, which forces e = 0 and O^2 = 1. That is impossible
    for odd O. Brute force agrees: min δ_S = 1, 3, 5, 7 for q = 5, 13, 17, 29, and 0 for q = 3 mod 4
    (`halfsys.py`).
  * Consequently the formal condition becomes g >= h(1 + 2δ_S) >= n/2 + h/2. But two marks with distinct kernels
    already have g <= n/2 + h/(2q). So the formal mechanism fails for p = 1 mod 4. This holds up to cancellation
    of top coefficients and up to L-independent extra gcd, which could matter only at very small L.
  * Forms with nonconstant coefficients (e.g. A^2 - tB^2) ramify infinity and leave the record's model.
  * **Status for p = 1 mod 4: open.**

## 5. Exact statement with the best verified rate

**Theorem (VI.4, as verified).** Let p = 3 mod 4 and m >= 2. Put k = ⌈log_p((p-1)m+1)⌉, and let q = p^r with
r odd and q > 2^{m-1}. Then for suitable signs, every L >= (k-1)(n-h)/(g-h) with L > (k-1)(n-g) gives m distinct
monic points of degree D = (n-g)L in F_{p^2}[t] with one norm in F_p[t] and pairwise deg(z_i - z_j) <= D/2.
Here n = q^k, h = n/q and g is the formal common degree, g >= 2^{1-m}n > h. (L = 2(k-1)n+1 is Astra's choice.)
**Corollary (padding, proved here).** Multiplying all points by a common monic polynomial of degree e keeps a
c = 0 cluster as long as e <= D - 2·max deg(z_i - z_j). That slack grows linearly in L. So **M_p(D,0) >= m for
every D >= D_1(m) = exp(O(m log m))**.
**Rate.** Take q >= 2^m (least odd r), so g - h >= 2^{-m}n. Then L <= max((k-1)2^m, (k-1)n) and
log D <= k log q + O(m) = (log 2 / log p + o(1)) m log m. Hence, for all large D,

  **M_p(D,0) >= (log_2 p - o(1)) · log D / log log D**; for p = 3 the constant is log_2 3 ≈ 1.585.

Astra's L = 2(k-1)n+1 has log D ≈ 2 log n and gives the constant (1/2)log_2 3. Caveat: for the smaller L the
distinctness after substitution needs either L > (k-1)(n-g), which costs the factor 2 again, or a direct check.
**Proved rate without caveat: (1/2)log_2 p.** **Plausible: log_2 p.** In the model's own terms (W = deg N ↔ log R)
this is M ≍ log W / log log W, i.e. M ≍ log log R / log log log R. That is doubly logarithmic and far below the
record's upper-bound shape W/log W.

Small values, this mechanism, p = 3:

| m | least D found | how | parent degree |
|---|---|---|---|
| 3 | 4 | q = 3, k = 2 | 8 |
| 4 | 6,981 | q = 27, k = 2, L = 11 (least in family) | 13,962 |
| 5 | 568,386 | q = 27, k = 3, L = 32, full gcd (independently checked) | 1,136,772 |
| 6 | 3,582,084 (formal-only: 6,880,250) | q = 27, k = 3, L = 191 (365) | 7,164,168 |
| 7, 8 | ~2.3·10^9, ~7.1·10^9 (formal estimate) | q = 243, k = 3 | |

## 6. Novelty and significance (blunt)

* **Novelty against the record: genuine.**
  * `README_V.md` §5 item 1 lists the question as open (P ≈ 0.005 of mattering).
  * `expert_note.md` §5 says "Odd characteristic is open: nothing was constructed and nothing excluded". The
    same place says the F_2 mechanism "needs many pairwise cofactor differences that are constant ... (Lemma P)".
  * VI.4 escapes Lemma P with nonconstant cofactor differences (X-degree (n+h)/2 - g), with exponent
    pattern 0/2 from 1+chi, and with a large common gcd.
  * I found no earlier construction of this kind in the record (grep for quadratic-character multiplicities,
    M_3, F_9(t)).
* **For the function-field programme: substantial.**
  * Unboundedness now holds over a fixed field of odd characteristic with the same conjugation type as Z[i].
  * The expert note's F3+ exemption list ("arguments needing odd characteristic, a reflection-type conjugation,
    or squarefree N in an essential way") loses its first two items. Only squarefreeness or bounded
    multiplicity, and characteristic 0, remain as exemptions. If VI.5 (squarefree over F_2, another block) is
    also confirmed, the surviving exemption is the conjunction of odd characteristic, reflection and squarefree
    (or bounded multiplicity) N, plus characteristic 0.
  * The expert note's sentence "Not evidence against the conjecture over Z: the mechanism needs ... constant
    [differences]" must be rewritten. Lemma P is no longer the reason. The reason is the q-linearised
    polynomials and Frobenius descent, which have no characteristic-0 or integer analogue.
* **For the integer goal: low (2/10).**
  * It gives no integer construction and no growth improvement. The growth rate (2/3+ε)log R/loglog R and
    Questions A and B are unchanged.
  * The mechanism needs: (i) P_H additive (F_q-linearised) for the degree drop; (ii) a characteristic-p
    Frobenius to descend from F_{q^2} to F_9; (iii) prime squares, with multiplicities about q^{k-2} at
    F_{q^2}-places; (iv) residue fields that grow with m.
  * It is a sharper filter, not evidence about Z.

## 7. Items for the record

1. **VI.4: CORRECT** (theorem, descent, norm identity, (VI.5), expectation, specialisation, count). It works for
   all p = 3 mod 4, and it fails formally for p = 1 mod 4 (half-system lemma, δ_S >= (q-1)/4).
2. **Computed (mine):**
   * four-mark reproduced exactly; the full gcd is G·t^26; least D in its family;
   * five-mark at L = 67 reproduced exactly (615,589 <= 615,797);
   * **5 points at D = 568,386** (independently checked);
   * **6 points at D = 3,582,084** (full-gcd formula), or 6,880,250 with the formal G only;
   * 3 points at D = 174 for p = 7.
3. **Corrections:**
   * the full gcd for k >= 3 includes about 17,500 of extra degree at F_{q^2}-places (exact formula §2.3);
   * the threshold q > 2^{m-1} is sufficient, not necessary;
   * N is not of bounded multiplicity for k >= 3, and most of its residue fields contain F_{q^2};
   * the rate is (1/2)log_2 3 proved with Astra's L, and log_2 3 with the optimised L (modulo distinctness).
4. **Open:**
   * squarefree or bounded-multiplicity N in odd characteristic (no route visible);
   * p = 1 mod 4;
   * whether sup_D M_3(D,0) grows faster than log D / log log D;
   * whether 5 points occur below D ≈ 5·10^5 by other means (the exhaustive range stops at D = 11).
