# Referee block "obstruct": Part VII, section VII.2 (obstructions at multiplicity one)

Referee report, round 31, Part VII. I wrote all the code myself: finite fields with Zech logarithms or direct
digit arithmetic, and polynomial gcd, division and products. I imported none of Astra's scripts.
Scripts and outputs are in
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref7_obstruct/`:

| file | purpose |
|---|---|
| `gf.h`, `findpoly.py` | own GF(p^N) (log/Zech and direct), polynomial gcd/divmod/derivative; primitive polynomials |
| `testA.cpp`, `A_q*.out` | VII.2(a): specialisations X=T(t) of VI.4 marks on V=iF_q[t]_{<k}; full gcd, reduced norm, squarefree test |
| `testA_nr.cpp` | the same family with non-real T (norm identity check) |
| `fracB.py` | VII.2(b): exponent-two fraction in the q=27 cap-two example |
| `formal.cpp`, `D_q11.out`, `E_q11_part.out` | VII.2(d),(e) exhaustive / sampled kernel fills; mode 2 = class F1 |
| `f1clust.cpp`, `f2cubic.cpp` | (f) cluster searches with ONE common divisor (whole-fibre translates; anisotropic cubic characters) |
| `t3.cpp`, `t3b.cpp` | (f) the characteristic-3 hyperplane-coset mechanism "T3" (identity, ρ-coincidence scans) |

A bug note, given here for transparency. My first gcd routine padded with the log-zero element (which is 1),
and I found this through a self-test (`tgcd.cpp`). It was fixed before every result reported here.
My first cluster search used a different gcd for each pair. That is illegitimate for m >= 3, since one divisor
must be common to all marks (standing caution). Every cluster number below uses one common divisor.

---

## Verdicts in one table

| claim | verdict |
|---|---|
| (a) no specialisation of VI.4 on the packed V (m>=16, q>=27) has squarefree norm after any common division | **CORRECT** (line by line). Threshold improvable: m >= 6 for q in {19,23,27}, m>=8 for q=11, m>=14 for q=7 with the exact Chebyshev constant. Numerically no squarefree case for any m>=2 with T nonconstant. |
| (b) unpadded VII.1 examples: exponent-two primes carry >= 1-2(q+1)/(m(q-1)) of deg N | **CORRECT**. Exponent-two labels are *exactly* the mixed labels. In the q=27 four-mark example the true fraction is 0.980 (bound 0.46). |
| (c) editing bound \|E\| >= d_*/5 | **CORRECT**, but **vacuous for q <= 41**: labels in exactly one kernel always lie in E, so \|E\| >= 2h(1-1/q) >= d_*/5 automatically when q <= 41. |
| (d) arbitrary fills of a single kernel fail for q >= 11 | **CORRECT**. Exhaustive over all fills for q=11, k=2: no critical pair. q=7 also has none (k=2 exhaustive, k=3 sampled). Critical pairs exist at q=3. |
| (e) two conjugate translates with arbitrary kernel fills fail for q >= 11 | **CORRECT** (expansion checked). Exhaustive q=7: none. q=11 (4.8M sampled pairs): none, min excess +38. Critical pairs exist at q=3. |
| (f) scope | Everything excluded is either the exact 1+χ exponent pattern, or a near-character pattern in formal X. One unexcluded squarefree mechanism (T3, char 3) reaches formal critical contact 4n/9 < n/2, *provided* hyperplane constants coincide. I could not produce the coincidences, and every other template I tried was bounded (<= 4). |

---

## (a) Specialisation obstruction: line-by-line check

Setting: V = iF_q[t]_{<k}, k >= 2, labels u ∈ F_q[t]_{<=b}, b=k-1. The raw marks are
F_j(T) = Π_u (T-iu)^{1+χ_j(u)}, with χ_j(u) = χ(ε_j ℓ_j(u)). B = max(b, deg T).

1. **Levels.** E_s(f) = {u : f^s | T - iu} is the zero set of an F_q-affine condition (u ↦ iu mod f^s is
   F_q-linear), and these sets are nested in s. ✓
2. **Character sums on a level.** ℓ_j is F_q-linear (its coefficients are in F_p ⊂ F_q). On an affine level u+W,
   either ℓ_j is constant on u+W, or it takes each value of F_q equally often, and then Σχ = 0. ✓
3. **Range >= 2 at a singleton mixed level.** We have v_f(F_j(T)) = Σ_s Σ_{u∈E_s} (1+χ_j(u)). Take a form i that
   is positive at u and a form j that is negative at u. Every level contributes a difference >= 0, since each
   term is either 0 or |E_s|·(±1) with the right sign. The singleton level contributes 1-(-1) = 2. ✓
   Dividing by a common G does not change the range. ✓
4. **Squarefree forces range <= 1.** With G = full gcd one has v_f(N) = v_f(P_V(T)^2) - min_j v_f(F_j) -
   min_j v_{f̄}(F_j). Using F̄_j(T) = F_j^{(-ε)}(T) (valid because T is real) and F_j^+F_j^- = P_V^2, this equals
   **exactly** the range max_j v_f(F_j) - min_j v_f(F_j). This is slightly stronger than the inequality Astra
   uses, and it gives a clean criterion. N ∈ F_p[t] is squarefree iff it is squarefree over F_{q^2} (perfect
   field). ✓
5. **Pairwise independence.** For nonproportional F_q-linear forms the map u ↦ (ℓ_i(u), ℓ_j(u)) is onto F_q^2,
   so the pair is uniform for uniform u. Hence the indicators "χ_i = +1" and "χ_j = +1" are pairwise independent,
   each with mean a = (q-1)/(2q). The signs ε only permute values. Chebyshev gives
   P(no positive) <= (1-a)/(ma), and the same holds for negative. So unmixed <= 2(q+1)/(m(q-1)) <= 4/m. ✓
6. **Incidence lower bound (VII.3).** If deg u = b, then deg(T - iu) = B, also when deg T = b, because i ∉ F_q. ✓
   Union bound: mixed labels of degree b have density >= 1 - 4/m - 1/q, and each contributes B incidences, all of
   which must lie at shared levels. ✓
7. **Upper bound.**
   - The direction of E_s(f) is {w : P^s | w, deg w <= b}. For f ≠ f̄ use f^s | w and f̄^s | w for real w. The
     level is shared iff es <= b, and then its density is q^{-es}. ✓
   - Even e: Σ_{P: deg e} e·q^{-e}·(…) <= 1 per even e <= b, using e·π_q(e) <= q^e. The s >= 2 tail is
     <= Σ_e q^{-e}/(1-q^{-e}) <= q^{-2}/(1-q^{-2})^2. ✓
   - Odd e: f = P is real, so P^s | T and P^s | u, and the total is <= Σ_P e·v_P(T)·q^{-1} <= deg T/q. ✓
     This needs T ≠ 0.
8. **Contradiction.** (1-4/m-1/q)B <= (1/2+1/q)B + q^{-2}/(1-q^{-2})^2 fails for m >= 16, q >= 27, B >= 1:
   the margin is 0.176B against 0.0014. ✓
   *Improvement:* with the exact constant 2(q+1)/(m(q-1)) in place of 4/m, it already fails for m >= 6
   (q = 19, 23, 27), m >= 8 (q = 11) and m >= 14 (q = 7). This also extends the theorem to q = 7, 11, 19, 23
   (computed by the one-line script in this report's console).
9. **T = 0.** F_j(X)/X evaluated at 0 is (unit) × (real polynomial). Equal norms then force unit multiples. ✓

**Is "real T ∈ F_q[t]" restrictive?** No. The hypothesis is wider than needed, since a legitimate cluster needs
z_j ∈ F_{p^2}[t], hence T ∈ F_p[t]. For non-real T = T_1 + iT_2 the formal norm identity breaks: the norms
differ by Π_u N(T-iu)^{χ_i-χ_j}, where N(T-iu) = T_1^2 + (T_2-u)^2 is no longer symmetric in u ↦ -u.
`testA_nr` (q = 7, 27) found unequal raw norms in **30/30 and 10/10** random tests.

**Numerics.** `testA` builds F_j(T) over F_{q^2}[t], divides by the **full** gcd, checks that all norms are equal
and lie in F_p[t], and tests gcd(N, N') = 1.
- Totals: q=3 (k=2,3,4), q=7 (k=2,3), q=11 (k=2) and q=27 (k=2) over 252 random (subset of forms of size >= 2,
  signs, T) instances. Every T nonconstant: **0 squarefree.**
- The only squarefree cases (q=3, k=2) have T constant, and there the mixed labels have degree < b. This is
  exactly the 1/q term.
- Conclusion: the theorem is true far below its threshold. The m >= 16 condition is an artefact of the
  Chebyshev route.

**Novelty and scope.** ref6_ternary (round 31) already *observed* that VI.4 parents carry multiplicity q^{k-2},
forced at F_{q^2}-rational places for this V. Astra's (a) is the first **proof** that no specialisation plus
common division of the exact VI.4 exponent pattern on the packed V is squarefree. It is narrow: it concerns one
space and one exponent formula. A simpler general remark covers the other natural space. On any V whose
T - v are pairwise coprime (the VII.1 device), every mixed label is automatically doubled (see (b)). So the
exact 1+χ pattern can never be squarefree there either.

## (b) Fraction bound

At a mixed label v, min_j e_j(v) = 0, and at -v as well, since the signs flip. So the exponent of T - v in
N_0 = R^2/(G(T)Ḡ(T)) is 2 - g(v) - g(-v) = 2. Conversely, exponent 2 forces g(v) = g(-v) = 0, i.e. v is mixed.
**Hence exponent-two labels = mixed labels exactly.**
- The disjointness and simplicity of the zeros of the T - v (VII.1) turn this into prime exponents.
- The fraction equals #mixed/(n-g), which is >= #mixed/n >= 1 - 2(q+1)/(m(q-1)). ✓
- `fracB.py` (own F_27): for forms x, y, x+y, x-y with signs + it gives g = 92 (matching the record) and 624
  mixed labels.
- So the true exponent-two fraction of deg N_0 is **0.9796** (bound 0.4615). For x, y (g = 365) it is 0.929;
  for x, y, x+y (g = 183) it is 0.952.

**CORRECT.** It is routine, and it quantifies how far cap two is from squarefree.

## (c) Editing bound

Every step checks:
- Q is supported on E, with |exponents| <= 4 and exponent sum 0, so its reduced denominator degree is <= 2|E|.
- Hence ord_∞(Q-1) <= 2|E| for nonconstant Q. Nonconstancy holds because labels in exactly one kernel carry an
  odd exponent.
- The labels with signs (+,-) outside E force a reduced degree >= d_* - |E|, so the depth is >= (d_*-|E|)/2.
- F_i/F_j has depth >= (n-h)/2 >= d_*/2.

**CORRECT.** **Observation (new, limits its content):** a label in exactly one of the two kernels has
χ_i - χ_j = ±1, which is odd, while 2(a_i-a_j) is even. So *every* such label lies in E, and
|E| >= 2h(1-1/q) holds for **every** selector pair. This is >= d_*/5 = n(q-1)^2/(20q^2) iff q <= 41. So the
editing bound says nothing for q <= 41. It has content only for q >= 43, which includes the q >= 2^m regime of
VII.1. Like (d) and (e), it is a **formal-X** statement: extra t-gcd after specialisation is not covered.

## (d) A single kernel

The verification:
- U_j^2/(F_j/X) = K_j/K̄_j, where the numerator and denominator each have degree (h-1)/2.
- So the denominator of Q is <= h - 1.
- Critical contact would need depth >= min(d_*/2, (n-h)/2) = d_*/2, and n(q-1)^2/(8q^2) > n/q - 1 for q >= 11.

**CORRECT.** Brute force (`formal.cpp` mode 0, K = F_{q^6}, V = i·span(1,θ) and similar, 6 form pairs × 4 sign
pairs, criticality deg(U_i-U_j) <= (D+g)/2 with g the exact common-root count, which is legitimate for pairs):

| q, k | fills | critical pairs | min excess |
|---|---|---|---|
| 3, 2 | all | 72 / 96 | -1 |
| 3, 3 | all | 144 / 6144 | -3 |
| 7, 2 | all | 0 / 1536 | +3 |
| 7, 3 | 2048 sampled per form | 0 / 7.2M | +63 |
| 11, 2 | all (32 per form) | 0 / 24576 | +15 |

So the theorem holds, and q >= 11 is sufficient but not sharp (q = 7 also fails). The record (ref6_ternary §3)
had already stated this heuristically ("drops only O(h)"). Astra's (d) is the proof.

## (e) Two conjugate translates

I checked the following:
- F_j = K_j((A-b)^s - c^s)((A+b)^s - c^s), since the positive fibres on V ± τ have A-values ±b + cλ, λ a square.
- (A-b)^s(A+b)^s = (A^2-b^2)^s, so the b/A terms cancel and the correction is O(X^{-2h}).
- The c^s terms are O(X^{-sh}) with s >= 5. P_V/A^q = 1 + O(X^{-(q-1)h}). So F_j/P_V = K_j/A_j + O(X^{-2h}). ✓
- d_* (here 2·n(q-1)^2/(4q^2)) >= 4h iff (q-1)^2 >= 8q, i.e. q >= 10. ✓
- deg(K_iA_j - K_jA_i) <= 0 and the expression vanishes at 0, hence it is 0. Taking norms: Ā_j = A_j because
  -ker = ker, so A_i = A_j. ✓

**CORRECT.** Brute force (mode 1, all 2^h kernel fills per form):

| q | k | critical pairs | min excess |
|---|---|---|---|
| 3 | 2 | 192/1536 | -5 |
| 7 | 2 | 0/393216 | +10 |
| 11 | 2 | 0/4,800,000 (6 form pairs × 4 sign pairs × 200k pairs; all 2048 fills of one mark against ~98 of the other; `E_q11_part.out`) | +38 |

## (f) Scope: what is excluded, what is not, and what I tried

**Excluded (precisely).**
1. The exact exponent pattern 1+χ(ε_jℓ_j):
   - by (a), on the packed V for every real T and every common division, with m >= 6..16 depending on q;
   - by (b) and the remark under (a), on any V where the T - v are pairwise coprime.
2. In the formal variable X only:
   - squarefree selector pairs on V\{0} deviating from the 2a = χ pattern on fewer than d_*/5 labels (content
     only for q >= 43);
   - character off one kernel with arbitrary kernel halves (q >= 11);
   - character off the kernel on V ± τ with arbitrary kernel choices (q >= 11).

**Not excluded.**
- Odd sign functions far from characters.
- Root sets other than V\{0} and V ± τ.
- Selections that use whole fibres with a non-character pattern.
- Any construction whose contact appears only after specialisation, through t-gcd.
- Anything not of the form "polynomial in X at T(t)".

**What I tried (≈ 1.5 h).** All cluster sizes below use ONE common divisor for the whole cluster.

* **F1: whole fibres on V ± τ with an arbitrary S ⊂ F_q per mark** (the natural generalisation of (e)).
  Max cluster: **3** at q=3 and **4** at q=7 (exhaustive DFS over 1008 distinct marks). The 4-cluster uses a single
  form: the reduced marks are polynomials in one A, so this is a composition, not growth.
* **F2: anisotropic cubic characters** σ = χ(Q(x)), Q an irreducible binary cubic. These are odd and nonzero on
  V\{0}, so genuinely squarefree. Max cluster: **4** at q=7 (224 marks) and **4** at q=11 (880 marks). No growth
  with q.
* **Higher characters with {0,1} exponents.** Odd half-systems that are unions of μ_d-cosets need d odd. They give
  f_S(Z) = g(Z^d) with a degree drop d < (q-1)/2, so they are strictly worse than χ, and the kernel stays
  self-conjugate. Cubic characters on their own are even (ψ(-1) = 1), so they are not admissible.
* **Additive characters / Artin–Schreier** (exponent from Tr ℓ). The self-conjugate set {Tr ℓ = 0} has size n/p,
  larger than the kernel. Even the doubled drop is only to (p+1)n/(2p), which needs common division of order
  n/(2p) and cannot hold for growing m. Worse than χ on every count.
* **Products over cosets of multiplicative subgroups (X^d - c), Dickson/Chebyshev, Rédei functions.**
  - All of these compose with **one** inner function: X^d, D_n, or the torus isogeny R_n = ((X+i)/(X-i))^n,
    whose fibres are the roots of (X+i)^n - c(X-i)^n.
  - Selections of whole fibres are g(inner). The contact problem becomes the same problem in the inner variable,
    with every degree multiplied by the inner degree. That gives no amplification.
  - For p ∤ deg, Ritt/Engstrom uniqueness of same-degree right factors (Lemma 1 of the record) excludes several
    inner functions. The only source of many close right factors is additive (subspace) polynomials. That is VI.4
    and T3 below.
  - The plain Rédei/Pell template z_a = (X+i)^a(X-i)^{n-a} has norm (X^2+1)^n. It is never squarefree, and its
    contact needs a - b ≡ 0 mod p^e with p^e ≳ n/2, giving <= 3 marks.
* **T3: the simplest unexcluded mechanism with genuine formal critical contact** (characteristic 3 only). This is
  new as far as I can grep.
  - Setup: let V be an anti-invariant F_3-space, n = 3^d, and τ anti-invariant ∉ V. Use the root set
    (V+τ) ∪ (V-τ), which is squarefree with no self-conjugate fibre. For an F_3-hyperplane H ⊂ V, let μ ∈ V\H,
    Y = P_H(X), ν = P_H(μ), so that P_V = Y^3 - ν^2Y.
  - Mark: select the coset H+τ+aμ and the two cosets H-τ-jμ (j ≠ a). This is one root per conjugate pair.
  - **Identity** (verified for all 39 and 120 marks, d = 3, 4, in `t3b`):
    U_{H,a} = P_V + ρY^2 - ρ^2Y - ρ^3 + ρν^2, with ρ = P_H(τ+aμ).
  - Since (Y_H - Y_{H'})^3 = ν_H^2Y_H - ν_{H'}^2Y_{H'}, we get deg(Y_H - Y_{H'}) <= n/9. Hence if
    ρ_{H,a} = ρ_{H',a'}, then **deg(U - U') <= 4n/9 < n/2 = D/2 with NO common division.** That is critical, and
    the norm is squarefree.
  - This is (VI.5) for "q = 3", with the self-conjugate kernel replaced by two conjugate translates.
  - **The obstruction is the coincidence of ρ.** If ρ_H ≠ ρ_{H'}, the difference has degree 2n/3.
    - A cluster of size m needs m hyperplanes with one common ρ. Equivalently, the F_3-additive equations
      (P_{H_1} - P_{H_j})(τ) = a_jν_j - a_1ν_1 must hold (j = 2..m) in the N-dimensional anti-invariant space,
      with τ having only N - d free coordinates mod V. Generically this is overdetermined already for m = 2.
    - Exhaustive scans over all τ (N = 7, d = 2..5) and over all hyperplanes and all a found **no coincidence at
      all**. The maximum ρ-multiplicity is 1.
    - With pair-specific cancellation only the common-divisor route remains, and it gives clusters of size 3
      (`t3`, d = 2, 4).
  - In the polynomial-label model needed for descent to F_3, ρ-equality becomes a divisor identity in F_3[t]. By
    the level argument of (a) it would force almost all labels' zeros to be shared. The only symmetries
    (t ↦ ±t + c) give at most 6 equal ρ's.
  - For p >= 7 the analogue needs e_1, …, e_{(p-1)/2} of the selected fibre values to agree: more conditions, so
    char 3 is the favourable case.

**Bottom line for (f).** Three things are now dead:
- the exact quadratic-character exponent pattern;
- its kernel fills;
- its two-translate fills.

The live classes are non-character odd sign functions and contact created by specialisation. I found no growing
squarefree mechanism. T3 shows that the doubled exponents in (VI.5) are *not* logically necessary for critical
formal contact. What it needs is a coincidence condition I could not satisfy.

## Significance for the integer problem

Low. VII.2 consists of negative results about one function-field construction family. Every argument uses
characteristic-p linearised polynomials, quadratic characters of F_q, or the Chebyshev count over F_q-points,
and none has an integer form. Squarefree uniformity over F_3[t] remains open, and so does the integer squarefree
case. This block neither helps nor hurts it beyond confirming that the known char-p escape needs high
multiplicity in its present form. VII.2's results are honest and correctly labelled. Its "HEURISTIC (scope)"
paragraph is accurate but understates how narrow the class is: (c)–(e) are formal-X only, and (c) is vacuous
below q = 43.
