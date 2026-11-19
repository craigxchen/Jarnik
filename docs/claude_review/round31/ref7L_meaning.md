# Referee report, round 31, Part VII (Astra, five marks over C(t)): block "meaning" (analyst)

Source: `astra_lemma_VII/answer.md`, from line 1565 (`# Part VII`) to the end (my copy is `ref7L_meaning/partVII.md`), together with `final7.txt` and `reply7.txt`.
My scripts are in `ref7L_meaning/code/`. They use my own exact Gaussian-rational polynomial arithmetic and Gaussian-integer Euclid, and no sympy:
* `pell_spec.py`, with output `pell_spec.out`;
* `jac.py`.

The literature PDF was downloaded to `ref7L_meaning/dl/`. It was only text-extracted, not copied anywhere else:
* Guo–Nguyen–Sun–Wang (GNSW), arXiv 2106.15881v3;
* Lehmann–McKinnon–Satriano (LMS), arXiv 2403.02480.

My own claims are marked **[R]**.

## 0. Verdict in six lines

1. Nothing Astra claims in VII.7 is false. The two-place Chebyshev version, the "what transfers" paragraph, the reading of GNSW Theorem 4 and the specialisation lemma are all correct. One correction is needed: the sparse gcd step **does** transfer to integers for a rank-one unit group (a single integer base), and it fails only for rank ≥ 2.
2. **No published function-field theorem implies the unrestricted five-point statement.** In our normalisation (W = 2D, |S| up to 2D) every candidate is vacuous or weaker than trivial. GNSW Theorem 4 is non-vacuous only when N has at most about D/c_0 distinct roots. In that range Mason–Stothers already gives h ≤ s−1, which is much stronger [R].
3. **The VII.2 gcd lemma is known.** It is the case S = {0,∞}, χ_S = 0, h̃ = 0 of Guo–Sun–Wang (JNT 231 (2022), Thm 8 = GNSW Thm 17). Its constants depend only on ε, n and the degrees, so it is uniform in constant coefficients. The VII.2 theorem as a whole is new to the record, and I found no published statement of it. Its method is standard: Schinzel/Zannier lacunary reduction, plus the GSW gcd, plus the record's round-29 Theorem D. Only the "one exponent vector tests the frozen family" device is a small new idea.
4. **VII.2 is not stronger than VI.1 in Diophantine content.** "Arbitrary complex coefficients" means height-0 coefficients, which is the analogue of *fixed* coefficients over Q. VII.2 is stronger only in dropping the honest-block (raw gcd) hypothesis. VI.1 is stronger in allowing any fixed K and S and moving coefficients of positive height. The best integer analogue I believe is provable is stated in §3.
5. **New consequence over Q [R], checked by my own code.** Specialise the VII.5 Pell family at integers T and let n vary. This gives **A_5(δ) false for every δ < 0**, i.e. α^ess_e(−K) ≤ 4 (unrestricted profile). So α^ess ∈ [10/3, 4]: Vojta predicts the upper end exactly, and Question B asks for > 15/4. This is not "trivially known". The record says explicitly that it is open whether a dense family reaches exponent 4 (`ref6L_subspace.md` l.92), and pigeonhole on r(N) gives only arcs R^(1−o(1)). It says nothing about the balanced (L_5°) version.
6. **Assessment of Question B: no change.** Characteristic 0 reproduces the integer status exactly: true on structured (fixed-S / sparse) classes, open in general, with sharp families exactly at Vojta's exponent. **Recommendation:** close the Astra Part VII thread. Record items 3–5. Put the characteristic-0 five-point question in the expert note as a clean function-field question for specialists (wording in §5). It is not a stepping stone the project should keep working on.

## 1. (a) VII.7, line by line

### 1.1 Two-place Chebyshev version: CORRECT (structural check; it inherits VII.2's frozen-family step)

* The pullback at + b_0 = (u+u^{-1})/2 turns each T_j into (u^j+u^{-j})/2. The numerators and denominators become Laurent polynomials in u with at most 2b terms, and the image curve does not change. The map u ↦ x has degree 2, so the height doubles to 2W.
* The place x = ∞ has the two preimages u = 0 and u = ∞, each unramified with local parameter ~ x^{-1}. So the contact is h at both.
* **The extra zero count checks [R].**
  * At u = ∞ (w = 1/u): d/du = −w²d/dw. The source log-derivative is O(w^{h+1}) and the parameter vectors are O(w^h), so the determinant vanishes to order ≥ 4h+1.
  * At u = 0: the source vector is O(u^{h−1}) and the parameters are O(u^h), so the order is ≥ 4h−1.
  * Neither place is charged.
  * Total zero count ≥ 8h, against ≤ 2W(1+o(1)) simple poles. A dominating frozen family therefore needs 8h ≤ 2W + o(W). This is the same threshold h > W/4 as one-place VII.2 at height W.
* Plain VII.2 applied to the pullback would need h > (2W)/4 and would fail, so the two-place refinement is necessary. Astra is right.
* **Meaning for integers [R].** At an integer x, T_j(x) = (λ^j+λ^{-j})/2 with λ = x + √(x²−1). For varying x this is an S-unit sum in **varying real quadratic orders**, the class where VI.1's cancellation estimate is false (README_VIL §1.1, the v_n example).
  * Over C(t) the "varying field" is one fixed double cover of P¹ of genus 0, so the function field cannot see this obstruction.
  * The Chebyshev version is therefore a pure function-field phenomenon with no integer counterpart. Fixed x with growing j is already VI.1 (fixed order). Varying x with bounded j is bounded degree (gauss Thm 4.3(ii), per degree only).

### 1.2 "What transfers": CORRECT WITH CORRECTIONS

Astra says the gcd step's integer counterpart is "divisibility by the value of a resultant", and that this does not give "four independent derivatives". The second part is right: the covering inequality (VII.2) is a degree count on a family of maps from P¹ and has no arithmetic analogue. In VI.1 its role was played by the Subspace Theorem with fixed S plus a logarithmic Jacobian.

**Correction (sharper than her wording) [R].** The resultant argument transfers **verbatim to a rank-one unit group**:
* Take u_i = T^{n_i} for one fixed integer (or Gaussian) base T, with coefficients of height o(L log T).
* Every place sees T^k as k times the same log-vector. So the n-weight is still one real number, and Dirichlet's approximation n ≈ L v/b' controls the archimedean size of R(T^n) by T^{O(L/Q)}·H(coefficients).
* The prime-to-denominator part of gcd(F(T^n), G(T^n)) divides R(T^n). This gives log gcd = o(L log T) off persistent relations n·u = 0.

It **fails for rank ≥ 2** (bases 2 and 3). There the weights are vectors in R^{|S|}, and the elementary bound is lost. This is exactly why Bugeaud–Corvaja–Zannier needed the Subspace Theorem for gcd(a^n−1, b^n−1), and Levin (Invent. 215 (2019)) for general S-unit gcds.

### 1.3 GNSW Theorem 4, read in the PDF (arXiv 2106.15881v3, pp. 3–4): description CORRECT

**Exact statement.** k is algebraically closed of characteristic 0, K = k(C) with C of genus g, and S ⊂ C(k) is finite. G ∈ K[x_1..x_n] is non-constant with no monomial and no repeated factors. Then for every ε > 0 there are c_0 and a proper Zariski-closed Z such that every u ∈ (O_S^*)^n \ Z satisfies either
* (i) max h(u_i) ≤ c_0 h̃(G) + max{1, 2g−2+|S|} (as displayed; the extraction loses brackets, and the proof's (4.4) has the form c_1h̃ + c_2 max{1, χ_S}); or
* (ii)(a) N_S(G(u)) − N_S^{(1)}(G(u)) ≤ ε max h(u_i), and (ii)(b) a truncated lower bound.

Uniformity:
* c_0 depends only on ε, n and deg G.
* Z depends only on ε and G, and is over k if G is (Z3), so it is **uniform in S**.

Astra's description, "alternative (i) bounds the unit height by c_0 max(1, 2g−2+|S|)", is right for constant G (h̃ = 0) under either bracketing.

**In our normalisation [R].** The anchored ratios q_j are S-units for S = roots of N. Since q_j(∞) = 1, ∞ need not be in S. So |S| = s ≤ 2D, χ_S = s − 2 and h(q_j) ≤ D.
* For squarefree N, s = 2D. Alternative (i) then reads max h(u_i) ≤ c_0(2D−2), which always holds once c_0 ≥ 1. **Vacuous.**
* (ii)(a) with G = x_1 − 1 would give h − 1 ≤ εD. That would be enormous, which is why it cannot be invoked: Astra is right that 4(ii) is unavailable.
* GNSW becomes non-vacuous only when s ≲ D/c_0. In that range **Mason–Stothers alone** gives a much stronger bound [R]. Take a reduced pair a = z_2/g, c = z_1/g of degree d, with deg(a−c) ≤ d−h. Mason gives d ≤ (s + d − h) − 1, so **h ≤ s − 1**.
* Hence h ≥ (1/2+δ)D forces N to have ≥ (1/2+δ)D + 1 distinct roots. GNSW adds nothing beyond Mason here. (The record has "Mason vacuous for squarefree N", `round29/ff_notes.md` l.23 and `round30/prove.md` l.365; the explicit h ≤ s−1 is a one-liner.)

GNSW **Theorem 1** (toric Vojta) is worse:
* Its proof ends with N_{D,S} ≥ h_D − (1+|S|)ε h_A − O(1) (pp. 23–24, (5.7)ff.). The O(1) and the |S|ε factor are absorbed only for fixed S, so with |S| ~ 2D it says nothing.
* It is also a statement for divisors D in normal crossings with the toric boundary. The exceptional divisor E of Bl_e is not of that form, so it does not even address approximation to e.

### 1.4 Other function-field results: what each gives here [R]

| result | what it gives with W = 2D, contact h, s = #roots(N) ≤ 2D | status |
|---|---|---|
| Mason–Stothers (pairs) | h ≤ s − 1 | vacuous for s ≥ D/2 + 1, e.g. squarefree N |
| Brownawell–Masser (S-unit sums, GNSW Thm 5) | height ≤ Σ N^{(n)} + n(n−1)/2·χ_S | vacuous: χ_S ~ 2D ~ W |
| Corvaja–Zannier 2008, JAG 17 (P² minus conic and two lines; gcd for n=2) | h(x) ≤ c·χ_S | vacuous; dimension 2 |
| Guo–Sun–Wang 2022 / Levin–Wang gcd | N_{S,gcd}(F(u),G(u)) ≤ ε max h(u), up to c·χ_S | with S = roots of N: vacuous. With S = {0,∞}: it **is** the VII.2 gcd lemma |
| GNSW Thms 1–4 | as in §1.3 | vacuous, or weaker than Mason |
| Capuano–Turchet (EJM 8 (2022), non-split surfaces over G_m²) | degree of φ(C) ≤ c·χ_S | vacuous; dimension 2 |
| Noguchi–Winkelmann–Yamanoi (holomorphic C → semiabelian) | applies to entire curves. Our curves are P¹ minus 2D+1 points: hyperbolic, not images of C. Yamanoi's function-field versions carry χ_S errors | inapplicable / vacuous |
| Wang's function-field Subspace (constant hyperplanes, S = {∞}, error ≤ 0) | non-vacuous (χ_{∞} = −1). Exceptional set: z_i linearly dependent over C (a dense union of hyperplanes, not one closed set; VII.4 shows that at 2h > D any relation is unique, real and of full support). Using X_i − X_1 (i = 2..5) and X_1 gives 4h ≤ 5D | weaker than trivial h ≤ D, and far from Liouville's h ≤ 3D/5 |

**Conclusion.** The unrestricted five-point statement over C(t) is **not** a consequence of any published result I could identify. The structural reason is that every truncated-counting result for S-integral points pays χ_S, and here χ_S is of the order of the height. The S = {∞} formulation (approximation at one place, D = E on the blow-up) is function-field Vojta for a non-toric divisor on a 4-fold, which is open.

### 1.5 Specialisation lemma (VII.7): CORRECT

Checked line by line:
* The Bézout identity Σa_iz_i = c ∈ Z[i] bounds the specialised Gaussian gcd by a divisor of c.
* Hence R'' ≍ T^D and arc ≍ T^{D−h} = R''^{1−h/D}.
* A curve escaping Z has all but finitely many integer points escaping Z.
* The strict slack absorbs the constants.
* The caveat "examples over C[t] alone give no specialisation" is correct and necessary: a dense family over C may have no Q(i)-points compatible with the real structure.

## 2. (b) Is VII.2 new?

* **The gcd estimate: known.**
  * GNSW quotes Guo–Sun–Wang JNT 231 (2022) 401–432, Thm 8, as Thm 17. For coprime F, G ∈ K[x], every u ∈ (O_S^*)^n with max h(u_i) ≥ c_1(h̃F+h̃G) + c_2χ_S^+ satisfies either h(u^m) ≤ c_3(…) + c_4χ_S^+ for some 0 ≠ m with |m| ≤ m_0, or N_{S,gcd}(F(u),G(u)) ≤ ε max h(u_i). The constants depend only on ε, n and the degrees.
  * With K = C(t), constant coefficients (h̃ = 0), S = {0,∞} (χ_S = 0) and u = t^n, this says: either m·n = 0 for some short m, or the common vanishing on C^* is ≤ εL. That is VII.2's lemma, uniform in the coefficients.
  * Earlier: Corvaja–Zannier 2008 (n = 2), the Ailon–Rudnick/Ostafe line, and Amoroso–Leroux–Sombra (cited by Astra).
  * Astra's Dirichlet + weight-homogeneous resultant proof is a nice elementary proof. I note it is correct in outline, but the lemma should be attributed.
* **Passage to a fixed torus:** standard. This is Schinzel's reduction for lacunary polynomials F(x^{n_1},…,x^{n_r}), as used by Zannier (Invent. Math. 2008), Fuchs–Zannier (JEMS 2012, "composite rational functions expressible with few terms") and Corvaja–Levin–Zannier (TAMS 375 (2022), intersections with cosets of one-parameter subtori via S-unit points over function fields). Removing the persistent exponent lattice is the usual Bombieri–Masser–Zannier-style structure step.
* **Covering inequality (VII.2):** the record's round-29 Theorem D, with d_+ and d_- kept separate.
* **New:** the device that one admissible exponent vector n_* makes the frozen family a covering family, whose evaluation image at t = 1 is the whole image of the fixed-shape map Φ, independent of n. I found no published statement of VII.2 itself (lacunary curves in G_m^4 with contact > W/4 at a point lie in one fixed closed set). Fuchs–Zannier and Zannier treat composition and reducibility, not approximation to a point.
* **Verdict:** new as a statement to the record and, as far as I can tell, to the literature. Its novelty is moderate in method and low in significance.

**Standing cautions applied to VII.2 (my reading of the quantifiers; the proof proper is another block's).**
* The contradiction scheme is legitimate. The l-th curve escapes the union of the first l Q-zero-sets, with W ≥ l. Every subsequence extraction (fixed combinatorial type, diagonal on the countably many integer vectors) keeps that property. The closure of Φ's image is Q-defined and proper, so it is the zero set of some enumerated polynomial, from which index on the curves escape it: contradiction.
* So Z_{b,ε} is one set for all degrees, coefficients and exponent vectors. It is ineffective, and W_0 is ineffective too, but not template-dependent.
* The gcd cost is **not** discarded: (VII.3) is about the primitive W.
* There is no drifting phase: constants are coefficients of the family.
* The fixed-S caution is respected because S = {0,∞} is intrinsic to monomials.

## 3. (c) VII.2 versus VI.1, and the integer statement worth stating

| | VII.2 (C(t)) | VI.1 (number field) |
|---|---|---|
| places of the additive ingredients | S = {0,∞} only (monomials t^n; one affine shift) | any fixed finite S, any fixed K ⊃ Q(i) |
| coefficients | arbitrary in C = **height 0** | height ≤ ρ log H (positive, small) |
| exponent vectors / degree | arbitrary, one Z | arbitrary, one Z(K,S,b,ε) |
| raw gcd / common factors | **arbitrary** (resultant shows they are o(W) off persistent relations) | honest blocks needed: ∏\|G_e\|² ≤ H primitive |
| conclusion | contact ≤ (1/4+ε)W off Z, i.e. "α ≥ 4" | Δ ≥ H^{−1/4−ε} off Z, α ≥ 4 |
| tool | covering-family degree count | Subspace (fixed S) + logarithmic Jacobian |

Blunt reading:
* "Moving complex coefficients" sounds stronger than it is. Over a function field with an infinite constant field, constants have height 0. Their number-field analogue is fixed (bounded-height) coefficients, of which there are finitely many (Northcott). VI.1 already allows positive-height coefficients.
* The single genuine gain of VII.2 over VI.1 is **no hypothesis on the raw common factors**.

**Strongest integer analogue I believe true and provable, (I).** Fix K ⊃ Q(i), a finite S, b and ε. Take a primitive five-tuple of norm H whose anchored ratios are q_j = A_j/B_j, where A_j and B_j are sums of ≤ b terms a·u with u ∈ O_{K,S}^* and **fixed** a ∈ K. Put no hypothesis on gcd(A_j, B_j) or on cancellation, and take H to be the primitive norm. Then off an ineffective proper closed Z(K,S,b,ε), Δ ≥ H^{−1/4−ε}.

Route:
1. Formal torus of the exponents, as in VII.2(ii).
2. A number-field gcd bound for S-unit arguments, to show that raw cancellation is H^{o(1)} off persistent multiplicative relations. That bound is Levin (Invent. 215 (2019)), with Corvaja–Zannier 2005 for n = 2. For a single base T it is the elementary resultant of §1.2.
3. VI.1's Subspace/Jacobian argument on the reduced blocks.

With moving coefficients of height ≤ ρ log H, one would need a moving-target gcd theorem over number fields. I do not know one, and did not verify the literature on it. Significance of (I): it removes the "honest block" caveat from the expert-note scope sentence. It does not reach generic Gaussian blocks.

**Not recommended:**
* "S = primes of the bases" with **varying** bases, e.g. z_j = F_j(T) with F_j b-sparse in an integer T → ∞ with growing exponents. This is the direct specialisation class of VII.2. It needs S-uniform Subspace or gcd (standing caution), and integer gcd(F(T),G(T)) can be as large as Res(F,G), with no fixed bound.
* With fixed bases, "primes of the bases" *is* a fixed S, and VI.1/(I) already cover it.

## 4. (d) Status in characteristic 0, and the essential exponent over Q

**Characteristic-0 status after Part VII.** The status list in the task is accurate:

| item | status |
|---|---|
| uniform bound | M(D,0) ≤ 24, (IV.2) |
| five-point A_5 analogue | OPEN for every δ ∈ (0, 1/10) |
| four points | false for δ < 1/6 |
| sparse class | true for every δ > 0 |
| dense five marks | at h = (D−1)/2, W = 4h+2 |

The same pattern holds over Q:
* VI.1 is true on fixed-S unit classes.
* A_4 is false for δ < 1/6.
* (*) and A_5 are open.
* Unit families sit at exponent exactly 4 on curves.

**Does it change the assessment of Question B? No.**
* Characteristic 0 is open at exactly the same place. The structured classes are decided in both settings by fixed-S arguments, and the sharp dense families sit at Vojta's exponent.
* This is what one expects if Vojta is true. It is not evidence that a proof exists.
* In characteristic 2 even the sparse class fails (VII.6). This confirms the record's calibration: a proof must use characteristic 0, and over Q, presumably Subspace-type arithmetic.
* The characteristic-0 function-field problem is itself an open instance of function-field Vojta (§1.4). It is a clean question for specialists, not an easier route.

**Essential exponent over Q: α^ess_e(−K) ≤ 4 (unrestricted profiles). New to the record, reproduced [R].**
* The five-mark family of VII.5 has D = 8n+1 and h = 4n, with 4 rational seeds, over Q(i)[t] with the real structure.
* My code (`pell_spec.py`) confirms (D, deg gcd, h) = (9,0,4), (17,0,8), (25,0,12) by its own polynomial gcd over Q(i).
* Integer specialisation t = T gives Gaussian gcd norms in {160, 2080, 4160}, i.e. bounded.
* log Δ / log H tends to −h/(2D) = −2n/(8n+1), as follows:

| n | T = 10 | T = 10^2 | T = 10^3 | T = 10^4 | target |
|---|---|---|---|---|---|
| 1 | −0.1905 | −0.1996 | −0.2074 | −0.2147 | −0.2222 |
| 2 | −0.2207 | −0.2241 | −0.2279 | −0.2317 | −0.2353 |
| 3 | −0.2305 | −0.2326 | −0.2351 | −0.2377 | −0.2400 |

  The approach is like a constant over log H.
* `jac.py`: at fixed T = 3, 7 and n = 1, 2, the Jacobian of the four anchored arguments with respect to the four seeds has four singular values of the same order. So the seed map is dominant at fixed T (a 3-mark control gives rank 3).

Argument:
1. By VII.7 with slack, for every proper closed Z and every n there are infinitely many primitive five-tuples off Z with Δ ≤ C H^{−2n/(8n+1)}, i.e. α ≤ 4 + 1/(2n).
2. Hence **A_5(δ) is false for every δ < 0**: five lattice points lie on arcs of length R^{1/2+ε}, Zariski densely, for every ε > 0.
3. Equivalently α^ess_e(−K) ≤ 4. With Cilleruelo–Córdoba, **α^ess ∈ [10/3, 4]**.
4. If Vojta (V_E(5) or LMS Prop 1.4) holds, then α^ess = 4 exactly.

Is it trivial?
* **No, as far as I can tell.** Pigeonhole on r(N) = R^{o(1)} points gives only arcs R^{1−o(1)}. Products of m Gaussian primes give 2^m points on arcs R^{1−1/m}, i.e. 4 points at R^{1/2} for m = 2 but not five, and these lie on special subvarieties.
* The record says (`ref6L_subspace.md` l.92): "Whether a *dense* unit family reaches exponent 4 is open". That remains open *at* exponent 4. The Pell sequence reaches 4 + 1/(2n) for every n, i.e. it attains 4 only as a limit over a sequence of families.
* LMS has no general upper bound α_ess ≤ dim. I read its §1–3: it only proves lower bounds from Vojta and curve upper bounds.
* The bound is consistent with the killer's dimension count (`round20/kill_gauss.md` §4: one-variable polynomial-block families can be dense only for arc exponent e > 1/2). The Pell family realises that count with e = 1/2 + 1/(2D).

Scope:
* The Pell profile uses 5 cuts, including two singleton cuts. It is **not balanced**, so it gives no upper bound for the balanced constant in (L_5°).
* It does not touch Question B, which concerns the lower bound > 15/4.
* It is a calibration: Question B's window (15/4, 4] is now known to be the whole remaining window, and no method can prove more than 4.

## 5. (e) Wording for the expert note, and recommendation

**Add to §3, after "Unconditional on a thin class":**

> *Function-field analogue (round 31, Part VII; refereed).* Over $\mathbb C(t)$ with the conjugate norm (distinct primitive $z_i$ of degree $D$, one norm, contact $h=D-\max\deg(z_i-z_j)$, $W=2D$), the five-point analogue of $A_5(\delta)$ is **open** for every $\delta\in(0,1/10)$. It is proved on the sparse class: ratios that are quotients of Laurent polynomials with $\le b$ terms each, with arbitrary complex coefficients, exponents and common factors, have ratio curves in one proper closed $Z_{b,\varepsilon}$ once $h\ge(1/4+\varepsilon)W$. The proof combines the Guo–Sun–Wang gcd estimate (JNT 2022, Thm 8; case $S=\{0,\infty\}$), Schinzel's lacunary reduction and the round-29 covering inequality $4h+1\le W$. This is the function-field counterpart of VI.1, with complex coefficients playing the role of fixed coefficients; it fails over $\mathbb F_4(t)/\mathbb F_2(t)$ with a fixed term budget. No published function-field theorem (Mason–Stothers, Brownawell–Masser, Corvaja–Zannier, Guo–Nguyen–Sun–Wang, Capuano–Turchet, Noguchi–Winkelmann–Yamanoi, Wang's Subspace) gives anything here: each pays $\chi_S=|S|-2$ with $|S|$ up to $2D$, or, for $S=\{\infty\}$, has exponent $n+1$. Mason alone gives $h\le\#\{\text{roots of }N\}-1$.

**Add to the approximation reading (§3, bullet list):**

> * *Upper bound* (round 31, Part VII.5 with the VII.7 specialisation; reproduced with exact arithmetic). $A_5(\delta)$ is **false for every $\delta<0$**: the polynomial Pell family $A(t)=\begin{pmatrix}t^2+1&t\\t&1\end{pmatrix}$, $\beta=t+2i$, $G_j=(1,i)A(t)^n(1,s_j)^t$, marks $\beta\prod\bar G_l$ and $\bar\beta G_jG_k\prod_{l\ne j,k}\bar G_l$ for the edges $12,13,23,14$, has $D=8n+1$, $h=4n$, primitive gcd and dense ratios. Specialised at integers $t=T$ it gives Zariski-dense primitive five-tuples with $\Delta\ll H^{-2n/(8n+1)}$. Hence $\alpha^{\rm ess}_e(-K)\le4$ for unrestricted profiles, so $\alpha^{\rm ess}_e(-K)\in[10/3,4]$, and Vojta predicts the upper end. The profile is unbalanced, so the balanced constant is not bounded. Exponent exactly 4 on a dense set remains open.

**Correct the scope sentence of VI.1** ("A large raw Gaussian gcd is covered only if…") by adding:

> (The raw-gcd hypothesis is probably removable for fixed coefficients using Levin's S-unit gcd theorem, or the elementary resultant bound when the unit group is generated by one base; not carried out.)

**Recommendation.** Close the Part VII thread and the characteristic-0 A_5 analogue as internal work.
* Its unrestricted case is open function-field Vojta for the exceptional divisor on a 4-fold, with no tool in reach. Its structured cases mirror VI.1, and its examples now pin the window from above.
* If anything continues, it should be one of the following.
  1. **(Outreach, not internal work)** Pose to function-field Vojta specialists (e.g. the Guo–Wang and Corvaja–Zannier school): "Distinct primitive $z_1,\dots,z_5\in\mathbb C[t]$ of degree $D$ with $z_j\bar z_j=N$ and $\deg(z_i-z_j)\le(1/2-\delta)D$: do their ratio curves lie in one proper closed subset of $\mathbb G_m^4$ for large $D$?" The trivial bounds are $h\le\min(3D/5,\ s-1)$, and the conjecture is $h\le(1/2+o(1))D$ off $Z$.
  2. **(Optional, low value)** Prove the integer statement (I) of §3 for one fixed base. This removes a caveat from VI.1 and touches nothing generic.
* No further rounds on characteristic 2 or on new sparse classes.

## 6. Literature checked

* **GNSW (arXiv 2106.15881v3):** Theorems 1–4, Theorem 5 (Brownawell–Masser form), Theorem 17 (Guo–Sun–Wang gcd), Theorem 18 and the proof of Theorem 1 (the (1+|S|)ε step), all read in the extracted text.
* **LMS (arXiv 2403.02480):** Def. 1.3, Prop. 1.4, Thm. 1.5, Cor. 3.4 and Lemma 3.6, read.
* **Search results only, full texts not read:** Corvaja–Levin–Zannier (TAMS 2022), Fuchs–Zannier (JEMS 2012), Zannier (Invent. 2008), Capuano–Turchet, Noguchi–Winkelmann–Yamanoi, Wang's function-field Subspace Theorem. Their statements are from GNSW's introduction and my memory.
* **Not verified:** whether a moving-target S-unit gcd theorem over number fields exists.
