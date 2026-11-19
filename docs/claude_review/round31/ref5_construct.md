# Referee report, round 31, block "construct": Part V.1 (Astra), the F_4(t)/F_2(t) cluster

Scripts are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref5_construct/`:
`f4arith.py` (my own F_2[t] and F_4[t] arithmetic: bit-vector carry-less multiplication with an 8-bit window table, an (A,B) pair model of
A+Bθ, list-based F_4[t] Euclid, and F_2[t] factorisation by squarefree decomposition, DDF and EDF), `check_V1.py`, `gcd_small.py` and
`vary_L.py`. Outputs are in `out4.txt`, `out5.txt` and `out6.txt`. I did not import or reuse Astra's script.

## Verdict

**V.1 is CORRECT.** I checked every step by hand and reproduced everything for m = 4, 8, 16, 32 and 64 with independent arithmetic.
This is a legitimate instance of the record's model M(D,c), transposed to the constant quadratic extension F_4(t)/F_2(t). That
extension is the only constant quadratic extension in characteristic 2, and in it the place at infinity is inert, so the form is
"definite". The instance is not degenerate: N is not a square, the points are monic and distinct, they are not unit multiples or
Frobenius twists of one another, and the common factor is tiny. Dividing it out only *strengthens* the cluster (c' = -deg g/2).
The constant field F_2 (F_4 upstairs) is **fixed**, independent of m.

Corrected or strengthened statement (all parts verified):

> **Theorem (V.1).** Let r >= 2, m = 2^r, Smax = (r-1)2^{r-1} = m(r-1)/2, and let L >= max(r+1, 2 Smax) = m(r-1) be an
> integer. Put Y = t^L + θ. Then the m polynomials z_H = prod_{v in H\{0}}(Y+v) in F_4[t], where H runs over the hyperplanes of
> V = F_2[t]_{<r+1} that do not contain 1, are distinct and monic of degree D = (m-1)L. They have the common norm
> N = prod_{v in V\{0,1}}(Y+v) in F_2[t] of degree 2D, and
> max_{H≠H'} deg(z_H - z_H') = D/2 + Smax - L/2.
> At L = m(r-1) this gives m points of degree D = m(m-1)(log_2 m - 1), on a parent of degree 2m(m-1)(log_2 m - 1), with maximum
> difference degree exactly D/2. Hence M_{F_4/F_2}(D,0) >= m for every D in (m-1)Z with D >= m(m-1)(log_2 m - 1), and
> sup_D M(D,0) = ∞. Along D_m = m(m-1)(log_2 m - 1) one has m ~ sqrt(2 D_m / log_2 D_m).

The extension to all L >= L_0 (not just the minimal L) is my addition. It follows from (V.1), and `vary_L.py` checks it
exactly for m = 4, 8, 16 over the whole range d <= L < 2Smax + 12, including the failure (c > 0) for L < 2Smax. A consequence of
the general L is noted in §5: at fixed m and large L the family reaches the Cilleruelo–Córdoba exponent.

## 1. Hand verification of V.1

1. **Additivity and recurrence.** P_{H+<b>}(X) = P_H(X)P_H(X+b), and, if P_H is additive, this equals P_H(X)(P_H(X)+P_H(b)). The
   new polynomial is then additive because squaring is additive in characteristic 2 and P_H(b) is a constant. The base case is
   P_{0}=X. Correct.
2. **Conjugation.** θ maps to θ+1 = θ^2 (Frobenius), and coefficients in F_2[t] are fixed. So Y maps to Y+1, and
   conj(z_H) = prod_{v in H\0}(Y+1+v) = prod_{w in (1+H)\{1}}(Y+w). Since ℓ(1) = 1, the coset 1+H equals V\H. Correct.
3. **Common norm.** z_H·conj(z_H) = prod_{v in (H\0) ∪ ((V\H)\1)}(Y+v) = prod_{V\{0,1}}(Y+v) for every H. This product is
   conjugation-stable (v maps to v+1 permutes V\{0,1}), so N is in F_2[t]. Pairing v with v+1 gives
   N = prod_{v even, v≠0} ((t^L+v)^2 + (t^L+v) + 1). For L >= d every factor has degree L, so deg N = (Q-2)L = 2D with
   D = (Q/2-1)L = (m-1)L. Correct.
4. **(V.1).** Take distinct H, H' and J = H ∩ H', with |J| = Q/4. V/J has four classes; [1] is neither 0 nor [u] nor [u'], so
   [1] = [u+u']. Subtracting the two recurrences gives z_H - z_H' = (P_J(u)+P_J(u'))P_J(Y)/Y = P_J(u+u')P_J(Y)/Y = P_J(1)P_J(Y)/Y.
   The last step holds because P_J is additive and vanishes on J. The factor P_J(1) = prod_{v in J}(1+v) is nonzero since 1 is not
   in J, and P_J(Y)/Y = prod_{J\0}(Y+v) is nonzero. The degree is (Q/4-1)L + Σ_{v in J} deg(1+v). Correct. Any codimension-2 J with
   1 not in J arises from exactly one pair (the two lines of V/J that avoid [1]).
5. **Maximisation.** Take an echelon basis of J with distinct leading degrees a_1 < ... < a_{d-2}. Then a_1 >= 1 (otherwise 1 would
   be in J), so a_i <= i+1. Exactly 2^{i-1} elements of J have degree a_i, adding 1 does not change their degree, and v = 0
   contributes 0. So Σ <= Σ_i 2^{i-1}(i+1) = n2^n with n = d-2, that is, (d-2)2^{d-2}. Equality holds for J = t^2F_2[t]_{<d-2},
   the joint kernel of c_0 and c_0+c_1. The half-degree condition is (Q/4-1)L + S <= (Q/4-1/2)L, that is, 2S <= L, so
   L_min = (d-2)2^{d-1} = m(r-1). This is >= d for d >= 3, and D = (m-1)L is even. Correct.

## 2. Independent computation (`check_V1.py`)

For each m I asserted the following, exactly: all norms equal N, and N is computed separately as an F_2[t] product of
quadratics; conj(z) = z·(complement) and z·conj(z) = (N,0); every point is monic of degree D; all m points are distinct; for
**every** pair, the polynomial identity z_H - z_H' = P_J(1)·prod_{J\0}(Y+v) and the degree formula (V.1) hold; and 2·maxdiff = D.

| m | L | D | parent | pairs | maxdiff | min diff | S_J distribution (S: #pairs) | gcd deg | c after gcd |
|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|
| 4 | 4 | 12 | 24 | 6 | 6 | 5 | 1:2, 2:4 | 0 | 0 |
| 8 | 16 | 112 | 224 | 28 | 56 | 53 | 5:4, 7:8, 8:16 | 1 (t+θ) | -1/2 |
| 16 | 48 | 720 | 1440 | 120 | 360 | 353 | 17:8, 21:16, 23:32, 24:64 | 0 | 0 |
| 32 | 128 | 3968 | 7936 | 496 | 1984 | 1969 | 49:16 ... 64:256 | 8 ((t+θ+1)^8) | -4 |
| 64 | 320 | 20160 | 40320 | 2016 | 10080 | 10049 | 129:32 ... 160:1024 | 30 | -15 |

This agrees with Astra's table (m <= 32) and adds m = 64. The m = 64 run took 20 s.

## 3. Legitimacy checks

* **Model.** The record's model (round29/ff.md §1; README_IV) asks for distinct z_j in k[t] of degree D with z_j·zbar_j = N
  (one N, multiplicities allowed) and deg(z_i - z_j) <= D/2 + c. The F_{q^2}/F_q version with Frobenius as conjugation is already in
  ff.md. Here k = F_4, the base is F_2, and the conjugation is Frobenius. The Gaussian form is degenerate in characteristic 2
  (X^2+Y^2 = (X+Y)^2), so the norm of the unique constant quadratic extension, A^2+AB+B^2, is the natural analogue. It is an
  Eisenstein-type form rather than literally the circle. That form is anisotropic over F_2(t) and over F_2((1/t)) because
  x^2+x+1 stays irreducible: the residue field at infinity is F_2. So the set of points of norm N is finite, as for a circle.
* **Monic and units.** Every z_H is monic. Units F_4^* have norm 1, but ωz_H - z_H has degree D, so unit multiples do not enter
  the cluster. No two z_H are unit multiples of each other, since all differences have degree < D.
* **Frobenius twists.** conj(z_H) (the complementary side) is never in the cluster, and min deg(conj z - z') = D - L
  (8, 96, 672, 3840, 19840). Note that squaring maps a c-cluster to a 2c-cluster in characteristic 2, so it preserves c = 0. That
  is a separate Frobenius phenomenon and is not what this construction uses.
* **Common factor (the standing gcd caution).** I computed g = gcd of all z_H in F_4[t] by Euclid for m <= 32. For m = 64 I used
  valuations at all irreducibles of degree <= d-1 (`gcd_small.py`), which is legitimate because a common irreducible must divide
  some v+v' in V; that program also cross-checks m <= 32. deg g = 0, 1, 0, 8, 30 for m = 4, ..., 64. After division,
  D' = D - deg g and maxdiff' = maxdiff - deg g, so c' = -deg g/2 <= 0. **Dividing makes the condition stronger, and the reduced
  cluster still satisfies it with room to spare** (the record's normalisation N2). No gcd cost is hidden. Pairwise gcds are large
  (about D/2, as (V.1) forces), but that is the structure, not a common factor.
* **Factorisation of N** (own factoriser, product re-checked), written as (prime degree, multiplicity): #primes.
  * m=4: (2,1):1, (4,1):2, (4,2):1, (6,1):1.
  * m=8: (2,3):1, (4,1):2, (4,2):1, (8,1):4, (12,1):2, (12,2):1, (14,1):2, (30,1):1, (32,1):2.
  * m=16: 52 distinct primes, with t^2+t+1 to the 12th power and the largest prime of degree 90.

  N is never a square and never squarefree. **Every prime of N has even degree**, so it is split in F_4(t) and there is no inert
  part. Proof: at a root α of a factor, α^L + v(α) is a root of x^2+x+1, so F_4 ⊂ F_2(α). So the parity condition (odd-degree
  primes to even powers) holds vacuously, and the points are honest exponent-vector choices over split primes. For m = 4 the
  monic points of norm N number 2·2·2·3·2 = 48, and the cluster takes 4 of them.

## 4. Novelty against the record

The record has the C(t)/R(t)/Q(i)(t) bound M(D,0) <= 24 (Part IV) and the characteristic-2 trace construction. That construction
is over F_q(t) with q = m+1 **growing**, and it is a **split**-form (hyperbola) cluster: "no anisotropic norm fits it"
(expert_note §5 Part IV; README_IV line 127, "constant field grows with m"). README_IV item F3' (line 399) asked exactly for
fixed-p constructions. A grep finds no earlier fixed-field anisotropic unbounded cluster. A web search found no function-field
version of the Cilleruelo–Granville arc problem (nearest: Yiasemides on sums UE^2+VF^2 in short intervals,
https://arxiv.org/pdf/2309.01290). The idea (subspace polynomials, with conjugation acting as translation by 1) grows out of
Astra's own trace seed, but the result is new to the record.

## 5. Significance

* **For the model: high.** The uniform function-field conjecture is a theorem in characteristic 0 (<= 24) and **false over the
  fixed field F_2(t), anisotropic form, c = 0**, with m ~ sqrt(2D/log_2 D). This answers the coordinator's (a) for p = 2 in the
  strongest form asked: fixed constant field, anisotropic norm, c = 0.
* **For the integers: a filter, not progress.** It removes any route to uniform CG whose every step also holds for
  F_4(t)/F_2(t): characteristic-free algebraic geometry, small-residue-field pigeonhole arguments, and geometry of numbers over
  function fields. That tightens F3: fixed small residue fields do **not** save a characteristic-free argument. The growth rate is
  untouched. m ~ sqrt(D/log D) is far below the pigeonhole ceiling ~D/log D (round29/moon Prop 4.3 shape), and the integer bound
  (2/3+eps) log R/loglog R is unchanged. Relevance to Z is limited by two facts: the mechanism needs characteristic 2 (additive
  polynomials), and F_4(t)/F_2(t) is unramified everywhere, while Q(i)/Q ramifies at 2.
* **Additional observation (mine; exact for m <= 16 in `vary_L.py`, and proved from (V.1)).** Fix m and let L grow. Then
  maxdiff = (m/2-1)L + Smax and D = (m-1)L. So m points satisfy deg(z_i - z_j) <= α_m D + Smax with
  α_m = (m/2-1)/(m-1) = 1/2 - 1/(2m-2). This is exactly the Cilleruelo–Córdoba exponent for m points
  (k = m-1, 1/2 - 1/(4⌊k/2⌋+2) = 1/2 - 1/(2m-2); overview item 4 and line 932: (k-1)/(2k-1) with k = m/2). The additive constant
  Smax is fixed for fixed m, so in this model the CC exponents are attained, up to R^{o(1)}, for every m = 2^r. Over Z they are
  known to be sharp only for k <= 3. In this model this is evidence that CC-type exponents cannot be improved by
  characteristic-free arguments. Whether a function-field CC theorem holds over F_2(t) with some constant I have not checked.

## 6. Claims outside this block

I did not referee V.2 ((V.3), (V.4), Ω(s^2 log s)). The odd-characteristic remark is correct as stated (conj^2 = id forces 2w = 0
for a translation w fixed by conjugation). It only excludes this mechanism, not other constructions over F_3(t).
