# Round 31, Part VII: judge report on multiplicity two, ramified models and squarefree tables

Source: Part VII of `astra_free_VII/answer.md` (heading "# Part VII", line 1461 to the end; scripts
`vii_cap2_check.py`, `vii_ramified_check.py`, `vii_squarefree_table.cpp/.out`, `vii_squarefree_witness.py`,
`vii_double_core.cpp/.out`; coordinator's question `reply7.txt`; Astra's summary `final7.txt`; her Part VI starts at
line 1258 of the same file). Six referees/analysts each took one block:

| key | report | scope |
|---|---|---|
| cap2 | `ref7_cap2.md` | VII.1, constant model: root space, separation, marks, width, padding, rate, certificates |
| ramified | `ref7_ramified.md` | VII.1, ramified model eta^2 = Delta; the ramified tables of VII.3 |
| obstruct | `ref7_obstruct.md` | VII.2: no squarefree specialisation of VI.4; exponent-two share; editing bound; kernel and translate replacements |
| census | `ref7_census.md` | VII.3: squarefree tables over F_9/F_3, the D=12 five-mark witness, the triple reduction, recurrence at D=12L |
| hadamard | `ref7_hadamard.md` | VII.4: incomplete Hadamard cores over F_2[t]; the x, x^3 two-core profile at e=10..16 |
| meaning | `ref7_meaning.md` | analyst: filter update, integer meaning, portability, next question, literature, wording |

The judge's scripts are in `judge_VII/`. They use the judge's own polynomial arithmetic over F_p and F_9 = F_3[i]
(plain Python lists, no libraries) and import nothing from Astra or the referees:
* `witness12.py`: the D=12 squarefree five-mark witness (decoding, norms, squarefreeness, distinct-degree
  factorisation, gaps, unit and conjugate coincidences). Output `witness12.out`.
* `separation.py`: the VII.1 separation step, R = prod_v (T - v) squarefree, constant and ramified models, p = 3, 7,
  both shapes of T. Output `separation.out`.
* `control.py`: negative controls that remove one ingredient of the device at a time. Output `control.out`.

**Every claim labelled [J] is the judge's own and has not been refereed.**

Labels follow Parts III-VI. "Proved" means a referee checked a complete proof line by line, or an exact computation
was reproduced with independent code. **No referee found an error in any claim Astra made, in the form she made
it.** The one GAP verdict (obstruct (f)) concerns a scope question the referee posed, not a claim of Astra's. I
rule on it in §1.3. Several statements need corrections of scope or threshold. They are quoted exactly below.

## 0. Verdict in five lines

1. **Theorem (VII.1).** For every p ≡ 3 mod 4, in the constant model (A^2+B^2 over F_{p^2}(t)/F_p(t)) and in the
   ramified model eta^2 = Delta (deg Delta even and positive, leading coefficient a nonsquare, infinity inert),
   critical clusters on parents with **every prime exponent at most two** are unbounded, with
   liminf M^{(<=2)}(D) loglog D/log D >= (1/2)log_2 p. The cap is reached only in a thin sense: about 98% of
   deg N lies on squared primes at m=4. The judge checked the separation step [J].
2. **Exactly squarefree N in odd characteristic is still open.** Tables: M_3^sf(D,0) = 2,4,3,4,3,4,4,4,4,4,4,5,4,5,4,5
   (D=1..16; referee extension 4,5,4 at D=17..19). It equals the all-parent maximum for every D <= 13, so the
   tables carry no information about boundedness. The judge verified the D=12 witness [J].
3. **VII.2 is correct but narrow.** It kills the exact VI.4 exponent pattern on the packed space under every
   specialisation, and some near-character edits in the formal variable. The editing bound is vacuous for q <= 41.
4. **VII.4 is correct, sharp at m=8 (false at m=4), and narrow.** It excludes single-matrix incomplete Hadamard
   cores with all rows. The binary growth scale stays open: 3M log_2 M - O(M) <= W(M) <= 2M^2 log_2 M.
5. **For the integer problem nothing changes.** The filter becomes F3+++ (§3(e)). No known integer argument used the
   escape that is now closed. Score for Part VII on the goal: **1.5/10**. Recommendation: close the function-field
   programme with one write-up exchange.

## 1. Referee-corrected statements per block

### 1.1 cap2: VII.1, constant model (referee: CORRECT on every item)

**Proved**
* *Root space.* Let r be the least odd integer with q = p^r >= 2^m, and k = ⌈log_p((p-1)m+1)⌉, n = q^k, h = n/q.
  V = { i^p C^p Σ_{j<k} x_j t^{pj} : x_j ∈ F_q } has q^k elements, all p-th powers, and is anti-invariant under the
  q-power conjugation (i^p = -i). The exponents 1+χ(ε_j ℓ_j(v)) are invariant under the p^2-power descent. φ and
  σ generate Gal(F_{q^2}/F_p) because r is odd.
* *Collisions.* v - w = i^p C^p u^p with 0 ≠ u ∈ F_q[t]_{<k}. Its roots lie in F_{q^d}, d <= k-1, so they divide
  Π_{j<k}(t^{q^j} - t), which divides C. At every zero of C all labels vanish.
* *Separation.* With s = L - p deg C, T = 1 + C^p t^s (p ∤ s) or T = 1 + C^p(t^s + t) (p | s). Then T = 1 on the
  collision set, T' = s t^{s-1} C^p or C^p, and d(T - v) = T' dt because dv = 0. So the finite zeros of the T - v
  are simple and pairwise disjoint, and R = Π_v (T - v) ∈ F_p[t] is squarefree. Referee additions: "since t | C" is
  needed only when p ∤ s and s >= 2. Astra's certificate uses a third shape, T = 1 + tC^36, also valid since
  T' = C^36. Exact certificates: R' = (Π_{v≠0}(-v))·T' (P_V is additive) and R ≡ 1 mod C.
* *Marks.* F_j = A(A^{(q-1)/2} - c^{(q-1)/2})^2 with A = P_{ker ℓ_j}; P_V = A^q - c^{q-1}A. G is the formal common
  divisor with g >= 2^{1-m}n. z_j = (F_j/G)(T) descend, are distinct, have degree (n-g)L, and have norm N_0 with
  N_0·Norm(G(T)) = R^2 (exact identity). The exponent at T - v is 2 - g(v) - g(-v). Exponent 2 occurs exactly at the
  mixed labels. **Correction (referee):** the degree statement needs L > b, with b = p(deg C + k - 1); it fails for
  L <= b.
* *Width.* 2 deg(z_i - z_j) <= (n+h-2g)L + b(n-h). **Exact threshold (referee):**
  L_min = max(b+1, p deg C + 2, ⌈b(n-h)/(g-h)⌉). The bound is attained by the worst pair in every computed instance,
  and L_min - 1 fails each time.
* *Padding and rate.* Reserved primes of degrees 1, 2, 4, ..., 2^{⌈log_2 n⌉} divide C. R ≡ 1 mod C, so they are
  coprime to N_0, and e^2 N_0 keeps cap two. The construction is valid for L >= 2^m(b+3) + p deg C + 2, that is,
  D >= C_p 2^m n^2. The log of the threshold is (2 log 2/log p) m log m + O(m), which gives (VII.1). The m-dependence
  is explicit and is not treated as uniform.

**Computed** (referee, independent F_9 FFT and own F_27)
* *Four marks* (q=27, k=2, forms x, y, x+y, x-y, all signs +, C = t^27 - t, T = 1 + tC^36, L = 973): D = 619801,
  parent degree 1239602, g = 92, deg R = 709317. Fully expanded: norms equal and in F_3[t], marks distinct and monic,
  pair degrees 307762 / 307720 against D/2 = 309900.5. **Astra's "upper bound" 307762 is the exact maximum.** Her own
  script did not expand the marks. It relied on the proof for norms and widths, so the referee's expansion goes beyond
  her certificate.
* *Smallest instances:* 2 marks at D=144 (p=q=3, L=36); 3 marks at D=206388 (q=27, L=378); 4 marks at D=578396
  (q=27, L=908, smaller than Astra's); p=7, 3 marks at D=14112. q=3 cannot give m >= 3 because g <= h.
* *Exponent-two share of deg N* (the rest has exponent one): 0.500 (q=3, m=2), 0.929 (q=27, m=2), 0.952 (m=3),
  0.980 (m=4). Astra's lower bound 1 - 2(q+1)/(m(q-1)) gives 0.46 at m=4.

**Judge's check of the step Astra asked to have prioritised [J].** (`separation.py`, `control.py`)
* *Argument.* (i) A collision v(α) = w(α), v ≠ w, forces C(α) = 0, by the root-field argument above. (ii) At C(α) = 0,
  T(α) = 1 and v(α) = 0 for every v, so no T - v vanishes there. (iii) Off C = 0, T'(α) ≠ 0: if p ∤ s,
  T' = s t^{s-1} C^p and t | C; if p | s, T' = C^p. (iv) (T - v)' = T' because v = (i C Σ x_j^{1/p} t^j)^p is a p-th
  power in F_{q^2}[t]. In the ramified model, v = (η C Σ x_j^{1/p} t^j)^p, and η^p = Δ^{(p-1)/2}η, so v is a p-th
  power in the function field, and Δ | C. Hence every zero of T - v is simple, and zeros of different factors are
  distinct. R has no repeated linear factor over the algebraic closure. Since R ∈ F_p[t], it is squarefree. I find
  no gap. Every zero of Δ is a zero of C, where R = 1, so ramified points never meet R.
* *Computation.* With q = p (labels x_j ∈ F_p, so each pair {v, -v} gives the real factor T^2 - v^2), exact
  gcd(R, R') = 1 and R ≡ 1 mod C in every run: p=3 constant model at L = 11, 12, 13, 15 (both shapes of T);
  p=3 ramified model, Δ = 2t^2+1, C = Δ(t^3-t), at L = 17, 18, 22, 24; p=7 constant model at L = 51, 57, 63
  (deg R up to 3087).
* *Controls.* Each ingredient is needed. With C = t (the factor t^p - t dropped), R fails to be squarefree for
  16/24 (p=3) and 12/24 (p=7) of s = 1..24. With non-p-th-power labels v = iCw (the full C kept), it fails for 7/24
  and 6/24. (A first draft of my script encoded t^p - t as t^p - 1 and reported false failures. The outputs in
  `separation.out` are from the corrected script.)

### 1.2 ramified: VII.1, ramified model; ramified tables (referee: CORRECT, two items with corrections)

**Proved**
* *Model.* η^2 = Δ, deg Δ even and positive, leading coefficient a nonsquare. Infinity is inert, the norm form is
  anisotropic, and the pole order is max(deg a, deg b + deg Δ/2). **Correction (qualification, referee):**
  "infinity has degree 2, so |Cl(O)| = 2h is even and O is never a PID. For η^2 = 2t^2+1, Cl(O) = Z/2, and the
  degree-1 ramified primes over t-1 and t+1 are non-principal (every norm has even degree)." Unlike (1+i) in Z[i].
  The construction does not depend on this, because N_0 is coprime to Δ. A non-maximal order is harmless for lower
  bounds.
* *Descent, anti-invariance, separation, width.* All re-derived. b = p(deg Δ/2 + deg C + k - 1). conj(v) = -v
  because η^p ↦ -η^p. The eta-part of R is 0, and N_0·Norm(G(T)) = R^2 holds exactly. The threshold
  L >= b(n-h)/(g-h) is exact for these T: q=3, m=2 at L=45; q=27, m=3 at L=392; q=27, m=4 at L=940 (L=930 fails).
  **Ramification plays no role in the mechanism.**

**Computed**
* *Four marks, Δ = 2t^2+1, D = 619801:* reproduced with marks expanded (her script certifies bounds only). The
  maximum pair degree is 308815, exactly her number. Also 3 marks at D=214032, and 2 marks at q=3, D=180.
* *Ramified tables (η^2 = 2t^2+1, D=1..10):* all parents 2,3,2,4,3,4,3,4,4,4; squarefree parents 1,2,2,3,3,4,3,4,3,4.
  Reproduced exactly with the referee's own C++ census, and by max-clique brute force for D <= 5. Definitions:
  D = max(deg A, 1 + deg B); one norm with any leading coefficient; ±z and conjugates count as distinct points.
* *Ramified infinity (referee's extension, computed only):* Δ = t (O = F_3[u], u ↦ -u, a PID) is critical at L=385
  (predicted 384.75), and Δ = t^3 - t at L=399. This is **outside the project's model** (infinity is not inert). It
  is the literature's model A^2 + tB^2. **Computed at m=3 only, not proved.**

**Correction to Astra's scope claim (g):** VII.1 closes README_VI §5 item 2 as posed (ramified finite primes, a
conjugation that is not a Frobenius). It does not cover p ≡ 1 mod 4, squarefree parents in any ramified model,
Artin-Schreier models at p = 2, or a proof for ramified infinity.

### 1.3 obstruct: VII.2 (referee: four CORRECT, one CORRECT WITH CORRECTIONS, one GAP on scope)

**Proved**
* *(a) No squarefree specialisation of VI.4.* Exact statement: for V = i F_q[t]_{<k}, k >= 2, q ≡ 3 mod 4, q >= 27,
  m >= 16 nonproportional forms, any signs, and every nonzero real T ∈ F_q[t], the specialised VI.4 marks do not
  have a squarefree common norm after division by any common polynomial. Referee sharpenings: with the full gcd,
  v_f(N) **equals** the valuation range, so the criterion is exact. With the exact Chebyshev constant the threshold
  drops to m >= 6 for q ∈ {19, 23, 27}, m >= 8 for q = 11, and m >= 14 for q = 7. Computed: 0 squarefree cases in
  252 random instances (q ∈ {3,7,11,27}, any m >= 2), apart from constant T at q=3. Non-real T breaks the norm
  identity (40/40). Scope: one space, one exponent formula.
* *(b) Exponent-two share >= 1 - 2(q+1)/(m(q-1)).* Correct. True values in §1.1 are far above the bound.
* *(c) Editing bound |E| >= d_*/5,* with d_* = n(q-1)^2/(4q^2). **Correction (referee; [J] rechecked):** labels in
  exactly one kernel always lie in E (χ_i - χ_j = ±1 is odd, 2(a_i - a_j) is even), so |E| >= 2h(1 - 1/q) for
  every selector pair. Since d_*/5 = h(q-1)^2/(20q), the bound has content only for q >= 43 (2(q-1)/q >= (q-1)^2/(20q)
  iff q <= 41). It is a statement in the formal variable X only. Extra t-gcd after specialisation is not covered.
* *(d) Single-kernel fills fail for q >= 11.* Correct. Brute force: 0 critical pairs of 24576 at q=11, k=2
  (minimum excess +15). Critical pairs exist at q=3.
* *(e) Two conjugate translates V ± τ fail for q >= 11.* Correct. q=7 exhaustive, 0 of 393k; q=11 sampled, 0 of 4.8M.

**Gap (referee's own scope question (f)): ruling.** The GAP verdict does not concern a claim of Astra's. Astra
wrote that squarefree ternary uniformity is undecided. The referee asked whether VII.2 excludes every mechanism for
growing squarefree clusters, and found that it does not. I rule: **no gap in Part VII; a gap in coverage, correctly
stated by Astra herself.** VII.2 excludes (i) the exact 1+χ exponent pattern on the packed space under any
specialisation, and (ii) near-character selections on V or V ± τ in the formal variable for q >= 11 (q >= 43 for
(c)). It does not exclude non-character odd sign functions, other root sets, or contact created by specialisation.
The referee's mechanism "T3" (characteristic 3, F_3-hyperplane cosets on V ± τ, formal contact 4n/9 < n/2 with a
squarefree norm, *provided* m hyperplanes share one constant ρ_H) is the referee's own claim. It was **not
reproduced** by me or by any other referee, and no instance with m >= 3 was found. Record it as an untested lead,
not as a result.

### 1.4 census: VII.3 (referee: CORRECT on every Astra claim)

**Proved / Computed**
* *Five-mark squarefree witness at D=12* on N = t^24+t^20+2t^18+t^16+2t^10+2t^2+2, and all 16 witnesses in
  `vii_squarefree_table.out`. Checked with the referee's own F_9[t] arithmetic. N = N_0(t^2) has 7 split primes of
  degrees 2,2,4,4,4,4,4. Gaps {5×3, 6×7}. All 24 five-clusters on all 9 squarefree D=12 parents have the same cut
  profile.
* **[J] Witness verified** (`witness12.py`, own F_9 = F_3[i] arithmetic, decoding digit = a + 3b ↦ a + bi): five
  monic marks of degree 12; norms A^2 + B^2 all equal to exactly t^24+t^20+2t^18+t^16+2t^10+2t^2+2; gcd(N, N') = 1;
  distinct-degree factorisation gives degree 4 in degree-2 primes and degree 20 in degree-4 primes (two plus five
  primes, all of even degree, hence split in the constant extension); pair gaps 5,5,5,6,6,6,6,6,6,6 <= 6 = D/2;
  no mark is a unit multiple of another or of a conjugate of another. Infinity is inert (constant extension, even
  degree).
* *Table* M_3^sf(D,0) = 2,4,3,4,3,4,4,4,4,4,4,5,4,5,4,5 for D = 1..16. Method 1 (full independent brute force,
  D <= 13) and method 2 (the referee's own triple reduction, D = 6..16) agree. Astra's triple and completer counts
  match exactly. Her exhaustiveness argument (split primes only; common block >= ⌈D/10⌉; z_i = G U_i Π conj(U_j);
  phase group of order 3^{K-1}) was checked line by line, and her C++ has no bug.
* *Recurrence:* five marks at every D = 12L, via t ↦ t^L (3 ∤ L) or t^L + t (3 | L). The proof uses N(0) = 2 ≠ 0
  or φ' = 1. The same argument gives five at multiples of 14, 16 and 18.
* *Referee extensions (new):* M_3^sf = 4, 5, 4 at D = 17, 18, 19. No six for D <= 19 (exhaustive). **M_3^sf = M_3 (all
  parents) for every D <= 13.** At D = 12 the all-parent maximum is 5, and the first five for any parent over F_9/F_3
  is at D = 12 (record: <= 4 for D <= 11; first constructed five at D = 568,386).
* *Cap two at small D:* identical to the squarefree table for D <= 13.

**Interpretation (census, endorsed).** Squarefreeness imposes no restriction at this scale, and all-parent clusters
are provably unbounded but only become visible at D ~ 10^5 - 10^6. **The plateau at 5 is no evidence of
boundedness.** Do not cite the tables as evidence either way.

### 1.5 hadamard: VII.4 (referee: CORRECT; one CORRECT WITH CORRECTIONS on the integer transfer)

**Proved**
* *(a) Theorem.* Over F_4[t]/F_2[t]: let the profile be a Hadamard matrix of order m >= 8 with an all-ones column,
  with all rows retained and a proper subset of the nonconstant columns. If the blocks are monic and nonreal with
  pairwise disjoint base-prime supports, the m products are never a critical cluster, after any common division.
  Checked line by line: the identity Π z_i^{(χ_v - χ_w)/2} = r_v^{m/4}; depth doubling by 2^a; the Gram matrix is a
  weighted Laplacian, so rank >= m/2 = n, which forces u = 1 and pairs of opposite sign; restriction to an n×n
  Hadamard matrix; contradiction. **Corrections (referee):** (1) the threshold m >= 8 is **sharp**: at m = 4 the
  statement is false (an explicit 4-point critical cluster with a missing column, f1 = t^2+t+θ, f2 = t^2+θ =
  (t+θ^2)^2, D = 4; referee-computed, not reproduced by me). (2) The one-line summary (Astra's item 9) drops the
  block-form and disjointness hypotheses and must state them. (3) A stronger form holds without disjointness if
  the undivided products are critical and the blocks are pairwise distinct.
  Brute force at m=8 (all GL(3,2) orbits, up to 1.68·10^9 assignments) found nothing. The m=4 controls found clusters.
* *(b) Scope.* With VI.2, only the complete set of m-1 columns survives for a single complete set of m >= 8 rows,
  at cost Ω(m^2 log m) (asymptotic; VI.2 is vacuous for m <= 8). Not covered: row subsets of larger Walsh
  matrices, several Hadamard matrices on the same rows (x, x^3), cores plus extra columns, shared supports,
  multi-level exponents, unions of designs.
* *(c) Criterion.* deg H' <= e/2 - 1 on the seven four-prime supports; key ⌊t^{2e-b-1}B/A⌋ with P = A^2 + tB^2.
  Derived independently: (PQ)' = (AE + BC)^2.

**Computed**
* *Two-core x, x^3 profile:* impossible at e = 10, 12, 14, 16. Counts (99,67,10,0), (335,319,34,0),
  (1161,4382,226,0), (4080,40546,1453,0) reproduced exactly (e = 10, 12 by direct exhaustion). Referee
  extension: also impossible at e = 18. A random model fits the counts within 3-10% and predicts the first
  admissible order-8 cores near e = 22-24. **The non-existence is consistent with chance and is not a structural
  obstruction.**

**Corrected statement (d):** "the Gram step is available over Z" is withdrawn. Over Z[i] the off-diagonal entries
are bounded by 4 log C > 0, not by 0, so the Laplacian sign argument fails. Only the divisibility and the identity
transfer. Rank >= M/2 is weaker than the record's existing rank results (round8 referee_stability Prop 4.1(b);
rank >= M-1 under (AI)). Depth doubling has no archimedean analogue because the unit circle has torsion.

### 1.6 meaning (analyst)

* *Filter (a):* endorsed with corrections. The F3+++ wording is in §3(e).
* *Integer meaning (b):* correct. No squarefree-specific step in the record uses odd characteristic. Lemma 3.1
  (rank >= n/4 - 2), ω(N') > M/4, height isolation (Lemma B), the round-7 dichotomy and the 2/3 chain lemma all hold
  verbatim over F_2[t] with squarefree N, where clusters are unbounded (VI.5). Route 5's Thm 5.7 (squarefree ⇒
  O(loglog R)) draws its power from ZE, which is archimedean.
* *Portability (c):* correct. The analyst verified P_V(T) = A_j(T)^7 - c_j^6 A_j(T) for all three forms in a p = 7
  instance (L=469, D=16884). The mechanism still rests on one parent with m right composition factors of the same
  degree, which is impossible in characteristic 0 (Lemma 1; Astra's own Laurent-root argument in VII.1 reproves it).
* *Literature (e):* nothing found on worst-case short-arc clusters over F_q[t]; the statistical literature uses
  A^2 + TB^2 (ramified infinity). Based on search snippets only.

### Withdrawn

Nothing of Astra's mathematics is withdrawn. Withdrawn or corrected wording: "the Gram step is available over Z"
(§1.5); the item-9 summary without block hypotheses (§1.5); the scope of the editing bound (vacuous for q <= 41);
"COMPUTED" for the four-mark certificate, which certified bounds and did not expand the marks (now fully expanded
by two referees, so the claim stands).

## 2. New relative to the record, and correct

Checked by grep in `docs/claude_review`: the record mentions bounded multiplicity only as an open escape
(`ref6_ternary.md` l.269, `ref6_meaning.md` l.21; expert note F3++ "squarefree or bounded-multiplicity N in odd
characteristic (**open**)", "the model A^2+tB^2 is untested").
1. **Cap-two unboundedness (VII.1)** for every p ≡ 3 mod 4, constant and ramified (finite) models, rate
   (1/2)log_2 p. The separation device is new: p-th-power labels with dv = 0, a factor C killing all collisions, and
   T ≡ 1 on the collision set.
2. **The first proof** that the VI.4 pattern cannot be specialised to a squarefree norm (VII.2(a)). The record had
   only observed forced multiplicities.
3. **The squarefree ternary table** D = 1..16 (and the referee's D = 17..19), and the **first five-point cluster over
   F_9/F_3 at D = 12**, for all parents too.
4. **The ramified tables** for η^2 = 2t^2 + 1, and the class-group remark (O is never a PID with inert infinity).
5. **The exclusion of incomplete single-matrix Hadamard cores** over F_2[t] (closes README_VI 3(d)), sharp at m=8.
   The two-core impossibility for e <= 18.
6. Referee additions: exact thresholds L_min; the smallest four-mark instance D = 578396; ramified infinity at m=3
   (computed); the q <= 41 vacuity of the editing bound; the m=4 counterexample.

None of these is a statement about Z.

## 3. Rulings

**(a) Is "unbounded with every prime exponent at most two, p ≡ 3 mod 4, constant and ramified models" a THEOREM?**
**Yes.** Exact statement: *Let p ≡ 3 mod 4. Take O = F_p[t, η] with η^2 = -1, or with η^2 = Δ(t), Δ ∈ F_p[t] of
positive even degree with nonsquare leading coefficient (infinity inert of degree 2; conjugation η ↦ -η fixes
F_p). Let M^{(<=2)}(D) be the largest number of distinct points of O of pole order D with one norm in F_p[t] whose
prime exponents are all at most 2, and with pairwise differences of pole order <= D/2. Then
liminf_{D→∞} M^{(<=2)}(D) loglog D/log D >= (1/2)log_2 p.* The m-point construction is valid for all
D >= C_{p,Δ} 2^m n^2, with log of the threshold (2 log 2/log p + o(1)) m log m. Caveats: (i) it is a lower bound
only; no upper bound of this order is known. (ii) The cap is reached in a thin sense: exponent-one primes carry a
fraction of deg N that tends to 0 as m grows (0.020 at m=4). The parents are almost squares. (iii) N_0 is coprime
to the ramified primes and to infinity, so ramification is filtered only weakly. (iv) Ramified infinity (A^2 + tB^2)
is computed at m=3, not proved. (v) p ≡ 1 mod 4 remains open.

**(b) Exactly squarefree N in odd characteristic: OPEN.** Tables: M_3^sf(D,0) = 2,4,3,4,3,4,4,4,4,4,4,5,4,5,4,5
(D=1..16, exhaustive by Astra's triple reduction, confirmed twice), 4,5,4 at D = 17..19 (referee); five recurs at
D = 12L, 14L, 16L, 18L. The table equals the all-parent table for D <= 13. Ramified η^2 = 2t^2+1: squarefree
1,2,2,3,3,4,3,4,3,4 (D=1..10). **The data decide nothing.** All-parent clusters are unbounded but only beyond
D ~ 10^5. The only proved obstruction is VII.2, which covers one family.

**(c) Obstructions of VII.2: correct, narrow.** (a) excludes the exact 1+χ pattern on V = iF_q[t]_{<k} for any real
specialisation (m >= 16, q >= 27; m >= 6 for q ∈ {19,23,27} with exact constants). (b) quantifies the doubling.
(c) has content only for q >= 43 and in the formal variable. (d)-(e) exclude single-kernel fills and two-translate
replacements for q >= 11, formally. None excludes non-character sign functions, other root sets or
specialisation-created contact. Squarefree uniformity over F_3[t] is not decided in either direction.

**(d) Binary: the exclusion is a theorem in block form with m >= 8 (sharp); the growth scale is open.**
3M log_2 M - O(M) <= W(M) <= 2M^2 log_2 M is unchanged. The two-core profile's absence for e <= 18 is consistent
with a random model predicting cores near e = 22-24.

**(e) Filter F3+++ (replaces F3++).** *Call an argument function-field valid if every step holds verbatim, with
degree for logarithm and the inert place at infinity for the archimedean place, in one of the settings:
(i) F_2[t] with the norm of F_4(t)/F_2(t), for all N or for squarefree N only; (ii) F_p[t], p ≡ 3 mod 4, with
A^2+B^2, for all N; (iii) F_p[t], p ≡ 3 mod 4, either with A^2+B^2 or with the order F_p[t, η], η^2 = Δ (deg Δ even
and positive, nonsquare leading coefficient, so infinity is inert; finite ramification allowed; conjugation fixes the
constants), for N with every prime exponent at most 2. No function-field-valid argument proves the uniform bound.
One valid in (i) cannot prove M <= f(W) with f = o(√(W/log W)), even for squarefree N. One valid in (ii) or (iii)
cannot prove M <= a log W/loglog W for any a < (1/2)log_2 p (it could still prove M ≪ W^ε, or M ≪ log W).
Admissible, hence insufficient on their own, in addition to the F3++ list: bounded prime multiplicity (cap two) in
odd characteristic; finite ramification; a conjugation that is not a coefficient Frobenius. Still available to a
proof over Z (none known to suffice): characteristic 0 used essentially (Newton linearisation, BMY, uniqueness of
composition factors, Lemma 1); archimedean order and curvature; Baker-type and transcendence input (ZE); exactly
squarefree N in odd characteristic (**open**); p ≡ 1 mod 4 (**open**); ramification at infinity (computed at m=3
only, not filtered). Growth-rate arguments are not filtered at the integer growth scale: W >= 3M log M - O(M)
transfers to all settings, and no function-field family reaches it.*

**(f) Integer problem: nothing changes.** The best growth bound stays (2/3+ε) log R/loglog R (W >= 3M log M - O(M)).
The uniform bound M(C), Questions A and B, and route 5 are unchanged. Bluntly: Part VII closes an escape clause that
no integer argument in the record was using. Every squarefree-specific step in the record already holds over F_2[t]
with squarefree N, where clusters are unbounded. The device (p-th-power labels, dv = 0, linearised root
polynomials, descent from F_q) has no integer analogue. The one escape that remains, exactly squarefree N in odd
characteristic, has no client among known strategies. Part VII is not evidence against the conjecture over Z,
for the same reason as Parts V-VI: the mechanism needs m same-degree right composition factors of one polynomial,
which is impossible in characteristic 0.

## 4. Expert-note sentences to update (recommendations only)

1. §5 item 12, F3++ bullet: replace by F3+++ (§3(e)). Remove "squarefree or bounded-multiplicity N in odd
   characteristic (**open**)" and "ramification" from the "still available" list. Insert "exactly squarefree N in odd
   characteristic (open)", "p ≡ 1 mod 4 (open)" and "ramification at infinity (computed at m=3 only)".
2. "the model A^2+tB^2 is untested" becomes: "ramified models with inert infinity are unbounded at cap two (Part
   VII); A^2+tB^2 (ramified infinity) is computed at m=3 only."
3. Add a Part VII bullet: "Part VII (`round31/README_VII.md`; six referee blocks; nothing false): unbounded at
   prime-exponent cap two for p ≡ 3 mod 4, constant and ramified models, rate (1/2)log_2 p, with about 98% of deg N
   on squared primes at m=4; exactly squarefree N in odd characteristic open; squarefree table over F_9/F_3 for
   D <= 16 (first five at D=12, equal to the all-parent table for D <= 13); VI.4 pattern provably not squarefree
   under specialisation; incomplete single-matrix Hadamard cores excluded for m >= 8 (false at m=4); binary growth
   scale open; no change for integers."
4. Round-31 table row (line 1012): append "cap two unbounded (constant and ramified), F3+++".
5. Sentences to avoid: "unbounded for squarefree parents in odd characteristic"; "false for A^2+tB^2"; "the census
   suggests boundedness"; "squarefreeness is the last hope for Z"; "evidence against the conjecture over Z"; "the
   Gram step ports to Z".

## 5. Open sub-questions, ranked by honest probability of affecting the integer problem

All probabilities are the judge's estimates.
1. **Squarefree mechanism at large D in F_3[t] (T3 or similar)**, P(affects Z) ≈ 0.5%. If a squarefree odd-
   characteristic family grows, the last algebraic escape closes and the filter is complete. If it provably cannot,
   it would have to use something no integer argument uses yet. *First step:* reproduce the T3 identity
   U_H = P_V + ρY^2 - ρ^2Y - ρ^3 + ρν^2 independently, then search exhaustively at q = 27 for three F_3-hyperplanes
   with a common ρ_H. If m = 3 is provably impossible, record T3 as excluded.
2. **Ramified infinity (A^2 + tB^2) proof**, P ≈ 0.2%. Only a write-up completion; it aligns the filter with the
   literature's model. *First step:* write the VII.1 proof for deg Δ odd. The width bound at the single
   ramified infinite place is the step to check.
3. **p ≡ 1 mod 4**, P ≈ 0.2%. The referee's lemma kills the mechanism formally (δ_S >= (q-1)/4). *First step:* an
   exhaustive census over F_5[t] with A^2 - 2B^2 type inert models to D ~ 12, compared with F_3.
4. **Binary growth scale**, P ≈ 0.3%. *First step:* meet-in-the-middle search at e = 20-24 for an admissible order-8
   Walsh core, as the hadamard referee proposes.
5. **Fixed exponent-one fraction** (Astra's proposed follow-up), P < 0.1%. It must be posed for normalised clusters,
   since common nonreal factors add exponent-one primes for free. *Do not pursue.*

**Recommendation: close the function-field programme.** Spend one more exchange with Astra on a self-contained
write-up: one theorem covering V.1, VI.4, VI.5 and VII.1 with exact thresholds; the composition-factor mechanism as a
lemma beside Lemma 1; the ramified-infinity case proved or explicitly left open; and an open list (exactly squarefree
odd characteristic with the tables to D=19, p ≡ 1 mod 4, Artin-Schreier models, binary growth scale). Integer effort
should return to characteristic-0 or archimedean input (BMY, Lemma 1, ZE).

## 6. Scores (0-10, on the integer goal)

| block | correctness | score on goal |
|---|---|---|
| VII.1 cap two, constant model | correct; thresholds made exact by referee | 2 |
| VII.1 ramified model | correct; class-group qualification | 1 |
| VII.2 obstructions | correct; narrow; (c) vacuous for q <= 41 | 1 |
| VII.3 squarefree tables | correct; uninformative at this scale | 1 |
| VII.4 Hadamard exclusion | correct; sharp at m=8; integer claim withdrawn | 0.5 |
| **Part VII overall** | no error found; careful and honest about what is undecided | **1.5** |

Quality as function-field mathematics is high (about 7/10). The low scores reflect distance from the integer goal,
not defects.
