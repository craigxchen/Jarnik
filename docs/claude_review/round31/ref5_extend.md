# Referee report, round 31, Part V, block "extend" (analyst, constructive)

Source: `astra_free_V/answer.md` lines 1138-1256 (Part V), `final5.txt`, `reply5.txt`. All scripts are mine, in
`ref5_extend/`. I did not import or run Astra's `v_fixed2_check.py`.

| script | what it does |
|---|---|
| `f2lib.py` | F_2[t] arithmetic (ints as bit vectors), F_4[t] as pairs (A,B) = A+Bθ, conjugation, norm, F_4[t] Euclid |
| `verify_V.py`, `delta_V.py` | V.1 reproduced for m = 4..32; minimality of L; δ_ij = deg(z_i−z_j) − deg gcd(z_i,z_j); margin identity |
| `m2.c` (+ `brute2.py`) | **exact M_2(D,0)** for the anisotropic norm A²+AB+B² of F_4(t)/F_2(t), D = 1..17 (Python cross-check D ≤ 7) |
| `m3.c` (+ `brute3.py`) | **exact M_3(D,0)** for A²+B² of F_9(t)/F_3(t), D = 1..10 (Python cross-check D ≤ 5); `m3w` = 128-bit version for D = 11 |
| `check_thm.py`, `irrlib.py` | checks of every ingredient of Theorem E below on 6,300 random (V, T, subfamily, sides) instances, 117,215 pairs |
| `thm_table.py`, `otherV.py` | Theorem E's bound against the Part V family; least D for all-hyperplane families of 400 random V |
| `search_trans.py`, `search5.py` | all products of linear factors Y+w over all sub-transversals (any subsets, non-linear) for dim V = 4, dim 5 |
| `search_trans3.py`, `search_conj3.py`, `abstract3.py`, `parity_check.py` | characteristic 3: the ± design, the conjugation-twisted design, abstract coset-clique search, the parity lemma |

Outputs: `m2_out.txt`, `m2_17.txt`, `m3_out.txt`, `m3_11.txt`, `outputs.txt`.

## 0. Verdict

1. **Exponent.** I could not beat the exponent 2. No construction I tried, and no search I ran, reached parent degree
   o(m²) over F_2. The best I can reach is Astra's 2m(m−1)(log₂m − 1). Random additive spaces V tie it but never
   beat it (§3).
2. **New theorem (Theorem E, proved, all ingredients machine-checked).** Astra's V.2 extends to *every* additive
   space V ⊂ F_2[t] containing 1, every T, every subfamily, either side of each form, unit twists, and division by any
   common factor. The extension has no thresholds and a much shorter proof. Any c = 0 cluster of s points of this
   type has
   reduced point degree D' ≥ (3/8)s''(s''−1)⌊⌊log₂s''⌋/2⌋, with s'' = ⌈s/4⌉ (s'' = s for one side per form and T ∉ V).
   The proof uses only the record's margin identity and the fact that the small **inert** primes t, t+1, t³+t+1, …
   divide the cofactor difference P_J(1) but can never be cancelled by a gcd. It needs no bound on split-prime gcds,
   so (V.4) is not used at all. On the Part V family the inert part alone accounts for 67-75% of D. So the quadratic
   loss sits at the small inert primes, which play the role of p ≡ 3 (mod 4).
3. **Exhaustive data.** M_2(D,0) = 2,3,2,4,3,4,3,4,3,4 for D = 1..10, and **M_2(D,0) = 4 for every D = 11..17**. So
   four points already occur at D = 4 (an instance of the same mechanism with two forms and both sides), and a fifth
   point needs D ≥ 18. M_3(D,0) = 2,4,3,4,3,4,4,4,4,4,4 for D = 1..11. In both fields the maximum is
   4 throughout the feasible range. The record's 8-point C(t) family at D = 8 has no analogue over F_2 or F_3 there.
4. **Odd characteristic: obstruction (Lemma P, proved).** In characteristic ≠ 2, if conjugation acts on the linear
   factors Y+w by a reflection w ↦ −w, then a "constant cofactor difference" forces |A| to be odd. Three transversals
   can therefore never be pairwise constant-difference, by a parity count. So the mechanism of V.1 gives at most 2
   points, in F_3(t) and in Z[i] alike. The conjugation-twisted design (Y real, pairs {w, w̄}) escapes the parity
   lemma for cosets of conj-stable F_9-subspaces only in odd size. Searches over it found nothing beyond 3 points.
5. **Integer analogue.** The identity with no characteristic-0 counterpart is not additivity alone. It is
   "Δ_A(X) = ∏_A(X+a) − ∏_A(X+σa) is constant, with |A| even". For a translation σ the X^{|A|−1} coefficient is
   |A|·w, which is nonzero in characteristic 0. For a reflection, Lemma P forces |A| odd and the triangle count kills
   it. Odd-size constant-difference sets do exist in Z (A = {1,2,−3}: Δ_A = −12), but they never form a triangle.
6. **Significance for the goal: low (1/10).** Theorem E closes the "other V / other T / subfamily" escape routes for
   hyperplane designs, and it says exactly where a growth-scale design would have to avoid losing degree. The
   unrestricted fixed-field growth question (b) stays open. Nothing transfers to Z.

## 1. Reproduction of V.1 and the margin reading

**CORRECT.** These are independent checks (`verify_V.py`, `delta_V.py`). For m = 4, 8, 16, 32 every z_H has norm N,
N is real, D = (m−1)L, the difference identity z_H − z_H' = P_J(1)P_J(Y)/Y holds exactly for all 650 pairs, and the
maximum difference degree is exactly D/2. With L − 1 instead of L the family is not a cluster (m = 4, 8, 16), so L in
(V.2) is the least L with T = t^L.

**Margin identity (record's, not new).** For any m points with one norm and degree D, with g_ij = gcd(z_i,z_j):
Σ_{i<j} deg g_ij ≥ D·m(m−2)/4. For each prime pair π, π̄ with exponent e, the shared degree of a pair is
e − |a_i − a_j|, and Σ|a_i − a_j| ≤ e m²/4. Inert primes are shared by every pair. Hence for a c = 0 cluster
Σ_{i<j} δ_ij ≤ mD/4, with δ_ij := deg(z_i − z_j) − deg g_ij ≥ 0, that is D ≥ 2(m−1)·avg δ_ij. This is the
function-field form of the record's "margin sum ≤ mW/4" (`ref_divisor.md` item 5, Parts II-III line 51). In
characteristic 2 one has z_i = g a, z_j = g ā with gcd(a, ā) = 1, so **δ_ij = deg(a + ā) = deg B(a)**, the θ-coordinate
of the coprime cofactor. Growth scale (D = O(m log m)) therefore needs the cofactor traces to have degree O(log m)
on average. Data for Part V: m = 16 gives avg δ = 17.7 (max 24) against the allowance D/(2(m−1)) = 24. The extremal
pair is exactly tight.

## 2. Theorem E: the quadratic loss holds for every additive space, every evaluation, every subfamily

**Setting.** V ⊂ F_2[t] is a finite F_2-subspace with 1 ∈ V. T ∈ F_2[t] is arbitrary, Y = T + θ. Take distinct linear
forms ℓ_1..ℓ_s on V with ℓ_i(1) = 1, c_i ∈ F_2, S_i = ℓ_i^{-1}(c_i), F_i = ∏_{v∈S_i}(Y+v), and u_i ∈ F_4^* the leading
coefficient of F_i (u_i = 1 unless T ∈ V). Put J_ij = ker ℓ_i ∩ ker ℓ_j and P_J(X) = ∏_{v∈J}(X+v).

**Theorem E (PROVED).** Let g be a common divisor of F_1..F_s, and suppose the monic points w_i = u_i^{-1}F_i/g form a
c = 0 cluster of degree D'. Let s'' be the size of a subfamily with one side per form and a common u_i, so
s'' ≥ ⌈s/4⌉ always, and s'' = s if T ∉ V and the forms are distinct. Then
```
D' ≥ (4/s'') Σ_{i<j in subfamily} ι(P_{J_ij}(1))
   ≥ (|V|(s''−1)/s'') Σ_{e odd, 2^{e−1}<s''/2} (e·π(e)/2^e)(s''/2 − 2^{e−1})
   ≥ (3/8) s''(s''−1) ⌊⌊log₂ s''⌋/2⌋ ,
```
where ι(x) is the degree of the part of x composed of irreducibles of odd degree, and π(e) is the number of
irreducibles of degree e. The parent degree is 2D'.

*Proof.* A subfamily of a cluster is a cluster with the same D', so we pass to the subfamily. There all u_i agree
and the u_i^{-1} factor cancels in differences.
(i) **Difference.** Choose a ∈ S_i ∩ S_j; it exists because distinct nonzero forms over F_2 are independent. Write
ker ℓ_i = J ∪ (u+J) and ker ℓ_j = J ∪ (u'+J). Since ℓ_i(1) = ℓ_j(1) = 1, we have u' ≡ u+1 (mod J). P_J is additive
with kernel J, so F_i = P_J(Y+a)·X and F_j = P_J(Y+a)·(X + P_J(1)), where X = P_J(Y+a+u). Hence
F_i − F_j = P_J(1)P_J(Y+a) and gcd(F_i,F_j) = P_J(Y+a)·gcd(X, P_J(1)). Since δ is unchanged by dividing both points by
g, this gives δ_ij = deg P_J(1) − deg gcd(X, P_J(1)).
(ii) **Inert primes are never cancelled.** An odd-degree irreducible P ∈ F_2[t] stays prime in F_4[t], because F_4
⊄ F_{2^{deg P}}. It never divides a factor Y + w = (T+w) + θ, since that would force P | 1. So v_P(gcd(X, P_J(1))) = 0,
and δ_ij ≥ ι(P_J(1)).
(iii) **Forced inert divisibility.** Let K_P = ker(ev_P : V → F_2[t]/P). The elements v ∈ J with P | 1+v are those
with ev_P(v) = 1. They form either the empty set or a coset of J ∩ K_P, of size ≥ |J|/2^e. The set is nonempty unless
ℓ_i or ℓ_j vanishes on K_P. To see this, the linear map k ↦ (ℓ_i(k), ℓ_j(k)) on K_P misses (1,1) only if one
coordinate is identically 0. Then (1+k) ∈ J for k with ℓ_i(k) = ℓ_j(k) = 1. Call a form P-special if it vanishes on
K_P. P-special forms with ℓ(1) = 1 factor through ev_P(V) ⊂ F_{2^e}, so there are at most 2^{e−1} of them. So at
least (s''−1)(s''/2 − 2^{e−1}) pairs have v_P(P_J(1)) ≥ |J|/2^e.
(iv) **Sum.** Apply the margin identity (D' ≥ (4/s'')Σδ). Use |J| = |V|/4 ≥ s''/2 (there are only |V|/2 forms with
ℓ(1) = 1). Use e·π(e)/2^e ≥ 3/4 for odd e (values 1, 3/4, 15/16, 63/64, …). For odd e with 2^{e+1} ≤ s'', bound
s''/2 − 2^{e−1} ≥ s''/4. There are ⌊⌊log₂s''⌋/2⌋ such e. □

**Checks (COMPUTED, `check_thm.py`, 4 seeds × 1,500 instances).** Random V of dimension 3-5 inside F_2[t]_{<8} were
used, with random T (small, T ∈ V included, and large), random subfamilies and random sides. Steps (i) and (ii) held on
all 117,215 pairs, with no violation of δ_ij ≥ ι(P_J(1)). Step (iii) held for every odd-degree P with deg P ≤ e and
every non-special pair. The gcd inequality Σ deg gcd ≥ D's(s−2)/4 held in every instance.
**Sharpness (`thm_table.py`).** On the Part V family:

| m | D (Part V) | (4/m)Σ ι(P_J(1)) | closed form (3/8)m(m−1)⌊log₂m/2⌋ |
|---:|---:|---:|---:|
| 4 | 12 | 8 | 4.5 |
| 8 | 112 | 84 | 21 |
| 16 | 720 | 490 | 180 |
| 32 | 3968 | 2790 | 744 |

So the first inequality captures 67-75% of the true D. The construction is within a factor of about 1.4 of the
inert-prime obstruction.

**Relation to V.2.** Astra's V.2 is CORRECT but narrower. It covers only V = F_2[t]_{<d}, holds only from d ≥ 48 and
s ≥ 6, and relies on the split-prime gcd bound (V.4). Theorem E has no thresholds. It covers all V ⊂ F_2[t] (with
1 ∈ V), and it covers the unit-twisted case that `ref5_lower.md` flagged as a gap: when T ∈ V, restrict to a common
unit. It never needs (V.4). The constants are comparable. Astra's conclusion D' ≥ Qsd/16 gives about s²(1+log₂s)/8.
Theorem E gives about (3/16)s²log₂s. Both are of order s² log s.

**What Theorem E does not cover.** Non-hyperplane transversals, where S_i \ S_j is not an affine coset. Division by
non-common factors (F_i/G next to F_j/Ḡ). Prime powers or extra factors. Several (V, T) mixed on one parent. Subsets
V ⊂ F_4[t] that are not real. For the last case: a conj-stable hyperplane avoiding 1 needs 1 ∉ (c−1)V, which forces
V = U ⊕ θU' with 1 ∉ U'. That only enlarges J for the same number of forms. I did not prove anything in these cases.

## 3. Searches for other designs over F_2 (COMPUTED)

* **All-hyperplane families, random V (`otherV.py`, 400 random V of dimension 3-5 in F_2[t]_{<10}).** I took
  T = t^L with the least admissible L and verified each cluster directly. The least D found for s = 4, 8, 16 is
  12, 112, 720, exactly Part V's values. L − 1 never works. No V beats the full space. This is consistent with
  Theorem E: every 2^j distinct nonzero polynomials have total degree ≥ (j−2)2^j, against Part V's j·2^j, so other V
  could gain at most a factor 1 − O(1/log m) in deg P_J(1).
* **Every non-linear sub-transversal family, dim V = 4 (`search_trans.py`).** For V = F_2[t]_{<4} with all 32 classes
  of T of degree ≤ 8 and T = t^L for L = 4..20, and for 3,000 random 4-dimensional V with random T, I enumerated all
  3^8 products ∏_{w∈S}(Y+w), with S a transversal of any union of pairs {w, w+1}. These were grouped by norm and top
  half. The best clusters found: 4 points at D = 4 (T = 0, W = t² + F_2[t]_{<2}, i.e. a 2-form family with both
  sides), 3 at D = 18, 8 at D = 112, nothing else. **No 5, 6, 7-point cluster and no 8-point cluster below D = 112
  arises from any product of linear factors in these spaces.**
* **dim V = 5 (`search5.py`).** All transversals of V\{0,1} minus at most one further pair, for T = t^L with
  L = 5..69 and 63 further T. Clusters found: 4 at D = 294, 8 at D = 322, 16 at D = 720 (Part V). Nothing new.

**Best exponent reached: 2** (parent ≈ 2m² log₂m). I found no route to o(m²).

## 4. Exhaustive small data (COMPUTED; `m2.c`, `m3.c`)

Method. Every cluster has a common leading coefficient (D > 2c+1), so WLOG all points are monic. For c = 0,
deg(z_i − z_j) ≤ ⌊D/2⌋ holds iff the coefficients of degree > ⌊D/2⌋ agree, which is an equivalence relation.
Therefore **M_p(D,0) = the largest fibre of z ↦ (N(z), top half of z)** over all monic z of degree D in F_{p²}[t]. The
C programs fix the top half, enumerate all bottoms (the norm is the norm of the top, plus a cross term linear in the
bottom, plus a precomputed norm of the bottom), sort, and take the longest run. Independent pure-Python brute force
agrees: p = 2 for D ≤ 7, p = 3 for D ≤ 5. For p = 2 I also re-verified the maximal cluster directly.

**p = 2, anisotropic norm A²+AB+B² (F_4(t)/F_2(t)):**

| D | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| M_2(D,0) | 2 | 3 | 2 | 4 | 3 | 4 | 3 | 4 | 3 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
| # maximal classes | 1 | 1 | 16 | 1 | 2 | 6 | 28 | 32 | 222 | 140 | 8 | 728 | 40 | 2556 | 424 | 12502 | 956 |

D = 17 means 4^17 ≈ 1.7·10^10 points (813 s). Example at D = 4: z = (t²+θ)(t²+t+θ) = t⁴+t³+1 + (t+1)θ and its three
partners. This is the V.1 mechanism with V' = {0,1,t,t+1}, T = t², two forms and both sides. **Reading.** The true
maximum is 4 from D = 4 through 17. Part V's 4-point member (D = 12) is not extremal in size or in D. A 5-point
cluster needs D ≥ 18, and Part V's first 8-point cluster is at D = 112. The growth of M_2 is invisible in the feasible
range.

**p = 3, anisotropic norm A²+B² (F_9(t)/F_3(t)):**

| D | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| M_3(D,0) | 2 | 4 | 3 | 4 | 3 | 4 | 4 | 4 | 4 | 4 | 4 |
| # maximal classes | 3 | 3 | 2 | 42 | 90 | 411 | 6 | 3768 | 110 | 34077 | 594 |

For even D the extremal N include t^{2D}+1. D = 10 means 9^10 ≈ 3.5·10^9 points (243 s). D = 11 means 9^11 ≈ 3.1·10^10
points (128-bit build, 3563 s, `m3_11.txt`).

So over both fixed small fields the exhaustive maximum in this range is 4. Over C(t) the record has 8 points at D = 8
(round 29, Thm C, also over F_q for q = 19..43). **Small residue fields visibly suppress clusters at small D**, which
fits the pigeonhole/inert-prime mechanism of §2.

## 5. Odd characteristic (p = 3)

**Lemma P (PROVED; checked by `parity_check.py`).** Let k have characteristic ≠ 2, and let A ⊂ k be finite with
A ∩ (−A) = ∅. If Δ_A(X) = ∏_{a∈A}(X+a) − ∏_{a∈A}(X−a) is constant in X, then |A| is odd.
*Proof.* For even |A|, ∏(X−a) = g_A(−X), so Δ_A = 2·(odd part of g_A). Then constant ⇒ zero ⇒ g_A is even ⇒
{−a} = {a} as multisets, a contradiction. □

**Corollary (triangle-free).** Let S_1, S_2, S_3 be transversals of the pairs {w, −w} of a symmetric set W, with
disagreement sets D_ij. Since D_23 = D_12 Δ D_13, we get |D_12| + |D_13| + |D_23| ≡ 0 (mod 2). So the three
differences S_i \ S_j cannot all have constant Δ. This is the only way the V.1 mechanism produces contact
(cofactor difference constant in Y), so in characteristic ≠ 2 with Ȳ = −Y it yields at most 2 mutually admissible
points. In characteristic 3, cosets J+u of F_3-subspaces do give constant Δ (= 2P_J(u)), but always with |A| = 3^j
odd. In characteristic 2, |J| = 2^j is even and there is no obstruction.

Checks. Over Z, all A ⊂ {±1..±7} of size ≤ 6: constant Δ occurs only for |A| = 1 or 3, e.g. {1,2,−3}, {1,3,−4}.
Over F_3[t]_{<3}, all A of size ≤ 4: only sizes 1 (13 sets) and 3 (52 sets, all affine lines). The abstract clique
search `abstract3.py` (large-L condition |A_ij| > n'/2 with coset differences) returns max 2 for F_3², max 2 for
F_3³, and max 2^{n−1} for the F_2 control (n = 3, 4), as it should.

**Other odd-characteristic designs tried.**
* Concrete ± design, Y = iT, V = F_3[t]_{<2} and F_3[t]_{<3}, T = t^L (L ≤ 13) and small T, all sub-transversals
  (≤ 1 dropped pair for dim 3) (`search_trans3.py`). The best found were 4 points at D = 2 (= M_3(2,0)), 3 at D = 3,
  and for dim 3 only 3 points at D = 20. For T = t^L there is **no cluster of size ≥ 3 at all**.
* Twisted design, Y = T real, factors T + w with w ∈ F_9[t] non-real, pairs {w, w̄} (`search_conj3.py`, 300 random
  conj-stable W of 4-12 pairs from F_9[t]_{<2}). Max 3 points (D = 10). In this design the coset mechanism needs
  A = β + J with J a conj-stable F_9-subspace; then Δ_A = P_J(β) − P_J(β̄) is constant. But |A| = 9^j is odd, so the
  same parity count applies. Even-size constant-difference sets exist (e.g. {a+bi, c−bi}), and they escape Lemma P.
  The Plotkin count then requires pairwise disagreement on just over half the pairs, together with constant Δ.
  For F_3-subspace cosets twisted by α ∉ F_3, Δ has degree |J|/3 in T, which forces disagreement ≥ 3/4 of the pairs,
  so at most 3 points.
* **Precise obstruction met.** In odd characteristic conjugation is an involution, and an affine involution of the
  Y-line is a reflection, never a translation (Astra's remark). Lemma P then makes the constant-difference graph
  triangle-free for every ± design. A Plotkin-tight family (distances n'/2 + O(1)) would need even-size
  constant-difference sets that do not come from subspace cosets, and I have no construction of such sets. Neither
  existence nor non-existence of unbounded clusters over F_3(t) is settled. The exhaustive maximum is 4 for D ≤ 11.

## 6. The integer analogue

Write the template in the form of V.1: factors T + v + i, or Y + v with Y = iT, over Z[i]. Conjugation acts on the
factor set by a reflection. Which identity of V.1 then fails?
* **(V.1a) Conjugation as translation**, Ȳ = Y + 1 with v̄ = v. This is impossible outside characteristic 2: applying
  the involution twice gives 2w = 0 (Astra; CORRECT).
* **(V.1b) Constant cofactor difference.** For a translation the X^{|A|−1} coefficient of Δ_A is |A|·w. It vanishes
  iff char k divides |A|, and additivity of P_J kills the lower coefficients. In characteristic 0 this coefficient is
  never 0, so Δ_A has degree |A|−1. In one-parameter templates with factors linear in T, the cluster condition then
  forces n' ≤ 2, i.e. at most 4 points. For the reflection Lemma P applies. Constant differences need odd |A|, as in
  the Prouhet-Tarry-Escott-type set {1,2,−3}, where (X+1)(X+2)(X−3) − (X−1)(X−2)(X+3) = −12. Such sets never form a
  triangle.
* So **the missing identity is "Δ_A constant for many pairwise differences A of even size"**. In characteristic 2 it
  comes from |J| = 2^j together with additivity. Over Z it would require a family of ideal symmetric PTE-type
  solutions closed under symmetric difference. By Part IV (orbifold BMY, ≤ 24 points for one-parameter templates in
  characteristic 0) no such family can be large. The record's polynomial templates (`growth/families.md`,
  `round30/bridge.md` Thm 1.2) are consistent with this: the 8-point records come from unit templates, not from
  linear-factor products.

## 7. Novelty and significance

* **Theorem E: new.** It is a strict generalisation of V.2: all V, all T, no thresholds, unit twists included, no
  (V.4). The margin identity it uses is the record's. Applying it with inert primes is new to the record. It is
  consistent with `ref5_reconcile.md`'s transferred pigeonhole bound, in which inert primes also carry the leading
  term, but it is different: that bound is linear in m log m and holds for all clusters, while E is quadratic and
  holds within hyperplane designs.
* **Exhaustive M_2, M_3 tables: new to the record** (as far as grep shows; `ref5_lower.md` enumerated within the
  family only).
* **Lemma P and its triangle corollary: new.** It is elementary, and it gives the precise reason why Part V does not
  port to p = 3 or to Z.
* **Significance for the integer conjecture: low.** None of this changes the growth rate (2/3+ε) log R/log log R, and
  none of it bears on Questions A or B. For the function-field side, the escape routes "other additive space, other
  evaluation, subfamily, common division" are now closed for hyperplane designs. What is left for growth scale over
  F_2 is a non-hyperplane design in which the cofactor traces B(a_ij) avoid small odd-degree primes. Nothing I found
  suggests one exists.
