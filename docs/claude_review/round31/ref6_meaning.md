# Referee report, round 31, Part VI, block "meaning": Part VI against the record and the literature

Scripts are in `round31/ref6_meaning/`. All the code is mine. I used no galois, sympy or flint.
* `fstar_q.py`, `fstar_big.py`: the exact configuration-free bound of `two_thirds/proof.md` Thm 5.2, transferred to F_p[t] with the norm of F_{p^2}(t)/F_p(t), for p = 2 and p = 3. Outputs: `fstar_q.out`, `fstar_big.out`.
* `formal5.py`, `formal5b.py`: exact X-level structure of the VI.4 clusters, using my own F_27 and F_243 arithmetic (log tables). It computes the exponent tables, the formal gcd g, the pair gcds, the cut weights of every 4- or 5-subtuple, δ_S, α_S, the balance, and H_Σ/w̄. Outputs: `formal5.out`, `q243.out`, `k5.out`.
* `f9gcd.py`: F_9[t] Euclid on lookup tables. I applied it to the ternary referee's saved four-mark points (`ref6_ternary/four_mark_points.npy`). I used that file as data only. Output: `four_actual.out`.

Other referees cover the correctness of VI.4, VI.5 and VI.1 (`ref6_ternary.md`, `ref6_squarefree.md`, `ref6_code.md`). All three report CORRECT on the main statements. Where this report relies on them, it says so.

## 0. Verdict

* **(a) The growth bound transfers to F_3[t] with the same constant:** W ≥ 8F_3*(M)/M − 4(M−1)c = **3M log_3 M − O(M)**, i.e. M ≤ (1/3+o(1)) W/log_3 W.
  * No prime is ramified. The discriminant −4 is a unit, and the constant extension is unramified at every place, including ∞, which is inert. So the term K log 2/2 of the integer proof is absent.
  * Exact minima: W_min = 4.8, 8.0, 42.0, 151.5, 440.2, 1021.1, 2516.1, 6345.1, 14120.5, 31326.2 at M = 5, 6, 16, 32, 64, 128, 256, 512, 1024, 2048.
  * My code reproduces the F_2 values of `ref5_reconcile` exactly (223, 674, 1611, 4210).
  * The ternary clusters lie **exponentially** far above the bound: W = exp(Θ(m log m)), so m ≍ log W/loglog W. The binary clusters (Part V, VI.5) lie only polynomially far above it: W ≍ m^2 log m, so m ≍ (W/log W)^{1/2}.
* **(b) Filters.**
  * VI.4 closes two of the Part V escape clauses: "odd characteristic" and "reflection-type conjugation". This holds for all N.
  * VI.5 closes "squarefree N", but only in characteristic 2.
  * Lemma P stays true but is no longer a barrier.
  * **Still open:** squarefree N, or bounded multiplicity, in odd characteristic. In VI.4 this is structurally excluded: a generic split prime of N carries the doubled pattern π² / π̄².
  * The F_3 filter is **weaker** than the F_2 filter. It blocks uniform bounds and bounds f(W) = o(log W/loglog W), but it does **not** block, say, M ≪ W^{1/2} or M ≪ W^ε.
* **(c) The F_3(t) analogues of A_5(δ) fail for every δ < 1/30 (PROVED modulo the char-free J4 + Schwartz–Zippel step, as in ref5_reconcile). So does (L_5°).**
  * The ternary clusters are *more* balanced than the binary ones. By Weil, the quadratic-character signs of any five pairwise non-proportional forms are asymptotically independent, so **every** 5-subtuple, dependent or not, has 15 nearly equal cuts.
  * δ_S = 1/30 − 8/(15q) + o(1/q). This is the binary formula with q in place of m.
  * Computed at q = 243 on all 56 five-subtuples of 8 rank-3 forms: δ_S ∈ [0.0310, 0.0316], balance ≤ 1.025, H_Σ/w̄ = 0.328 ≈ 80/q.
  * Four-mark actual data: 7 cuts of weights 1001 (×5) and 988 (×2).
  * The (L_8) analogue fails asymptotically, pending the dictionary, as in characteristic 2.
* **(d) Literature: not found.**
  * The function-field short-interval and sector literature asks a different, statistical question. It uses the ramified analogue A²+TB² (F_q[√−T]) and large-q limits (Bary-Soroker–Smilansky–Wolf, Bank–Bary-Soroker–Fehm, Gorodetsky–Rodgers, Rudnick–Waxman, Yiasemides).
  * I found no worst-case short-arc result over F_q[t].
  * The ingredients of VI.4 are classical: Ore's subspace/linearised polynomials, BKR's "subspace polynomials with common top coefficients", and Π_{s square}(Y − sc) = Y^{(q−1)/2} − c^{(q−1)/2}. The combination with quadratic-character multiplicities and coefficient-Frobenius descent was not found.
* **(e) The heuristic "fakes exist over Z because Z has no derivation, and abc would kill them" is wrong as stated.**
  * VI.1 uses the char-2 fact that d/dt kills squares, so only the *odd* valuation support counts. That is a mod-2 Mason inequality, and its integer analogue is false (ref6_code's Pell family).
  * Ordinary Mason (characteristic 0), or abc over Z, gives only a support condition. The barrier profile (the skew double core) **passes** it (record `lw.md` Thm 5.1; ref6_code (e3)).
  * The record has already tested abc against the fakes: `growth/conditional.md` Prop 8.3 (abc-type inequalities are vacuous on squarefree blocks), `round8/weakest.md`, the abc pencils of rounds 18–19, and this round's ref6_code.
* **(f)** Wording is in §6. Do **not** write "the conjecture fails in odd characteristic" without qualifiers. The correct statement is "for the norm form of the constant extension F_{p²}(t)/F_p(t), p ≡ 3 mod 4, all N allowed".

## 1. (a) Transfer of the 2/3 argument to F_3[t] with A²+B²

**Places.**
* A prime P of F_3[t] of degree e has residue field F_{3^e}. Since −1 is a square in F_{3^e} iff e is even:
  * e **odd**: P is inert, and the fibre {z mod P : z z̄ ≡ N} has 3^e+1 elements.
  * e **even**: P splits, and the fibre has 3^e−1 elements.
* Lemma 3.2 of `two_thirds` goes through word for word. The level-1 counts come from F_{9^e}^* → F_{3^e}^* and from uv = n. Hensel lifting works because the gradient (2x, 2y) is nonzero mod P (characteristic ≠ 2). This is *closer* to the integers than the F_2 transfer, which needed the polar form (y, x).
* **Ramification: none.** The discriminant of X²+1 is −4 = 2, a unit in F_3, so no finite prime ramifies. F_9(t)/F_3(t) is a constant extension, so it is unramified everywhere. The place ∞ (degree 1, odd) is inert. So Lemma 3.3(R) and the term K log2/2 disappear.
* This is a structural difference from Z[i], where 2 ramifies. It affects only O(M²) in F*.

**Steps.**
* Lemma 3.1: deg c ≥ Σ_P v_P(c) deg P, by unique factorisation in F_9[t].
* Lemmas 3.3 and 3.4 and Cor 3.6 carry over, with log p replaced by deg P and log R by D.
* The cut identity (§2), the chain lemma (§4, pure real analysis), Lemma 4.3 and Thm 5.1/5.2 carry over without change.
* Result: F_3*(M) ≤ MD/4 + C(M,2)c, where
  F_3*(M) = Σ_{e odd} e·I_3(e)·Σ_a E(M,(3^e+1)3^{e(a−1)}) + Σ_{e even} e·I_3(e)·min{Σ_a E(M,(3^e−1)3^{e(a−1)}), Ψ_{3^e}(M)}.

**Leading constant.**
* e·I_3(e) = 3^e − O(3^{e/2}), and E(M,Q) ≥ M²/2Q − M/2.
* Each odd degree e ≤ log_3 M contributes M²/2 + O(M·3^e). Each even degree contributes ≥ (c_Q M²/4) with c_Q = 1/Q + O(Q^{−2}).
* Summing over the ≈ (log_3 M)/2 degrees of each parity gives F_3*(M) ≥ (3/8)M² log_3 M − O(M²).
* Hence **W ≥ 3M log_3 M − O(M) − 4(M−1)c.**
* In natural units (|N| = 3^W) this reads log|N| ≥ 3M log M, the integer constant exactly. The inert/split degree parity plays the role of the mod-4 split, and Σ_{deg ≤ E, odd} e I(e)/3^e ≈ E/2 replaces Mertens.

**Exact values** (`fstar_q.out`, `fstar_big.out`; W_min = 8F*/M at c = 0):

| M | p=2: W_min | p=3: W_min | 3M log_3 M | p=3 ratio | ternary cluster (best known parent degree) |
|---|---|---|---|---|---|
| 4 | 4.0 | 0.0 (E(4,4)=0; indeed M_3(2,0)=4) | 15.1 | — | 4 (exhaustive), VI.4: 13,962 |
| 5 | 6.8 | 4.8 | 22.0 | 0.22 | 1,136,772 (ref6_ternary) |
| 6 | 8.0 | 8.0 | 29.4 | 0.27 | 7,164,168 (ref6_ternary, full-gcd formula) |
| 32 | **223.0** | 151.5 | 302.8 | 0.50 | — |
| 64 | **674.0** | 440.2 | 726.8 | 0.61 | — |
| 256 | **4210.0** | 2516.1 | 3876.4 | 0.65 | — |
| 2048 | — | 31326.2 | 42640.8 | 0.735 | — |

The bold p = 2 entries agree with ref5_reconcile. For p = 3, (W_min − 3M log_3 M)/M = −4.6, −5.1, −5.5 at M = 512, 1024, 2048. That is consistent with −O(M) plus slow, Mertens-like convergence, as the integer table (6.00 at M = 200) also shows.

**Comparison.**
* Ternary: log W = Θ(m log m), with M_3(D,0) ≥ ((1/2)log_2 3 − o(1)) log D/loglog D proved (ref6_ternary §5). The gap to the transferred upper bound W/log W is **exponential**.
* Binary: W ≍ 2m² log_2 m (Part V) or 8m² log_2 m (VI.5, squarefree), a factor ≍ m above the bound.
* So over F_3 the true growth rate of sup M_3 is known only to lie between log W/loglog W and W/log W. Over F_2 it lies between √(W/log W) and W/log W.
* Neither model says anything about the growth constant 3.

**Verdict (a): CORRECT.** The transfer is routine and the code is independent. It is new only as a fixed-field odd-characteristic statement.

## 2. (b) The filters after VI.4 and VI.5

**F3++ (ternary filter, recommended wording).** Call an argument *F_3-valid* if every step holds verbatim under all of the following:
* F_3[t] with the form A²+B² (norm of F_9(t)/F_3(t), conjugation i ↦ −i, i.e. coefficient Frobenius);
* all N, prime powers allowed;
* deg in place of log, and the inert place ∞ in place of the archimedean place.

Then **no F_3-valid argument proves the uniform bound, nor any bound M ≤ f(W) with f(W) ≤ (½log_2 3 − ε) log W/loglog W for all large W.**
* This follows from VI.4 together with ref6_ternary's padding corollary. The constant log_2 3 is plausible with a better L, but not proved.
* The same holds for every p ≡ 3 mod 4.
* Unlike F3+ over F_2, this does **not** exclude sub-polynomial growth bounds such as M ≪ W^{1/2} or M ≪ W^ε.

**F3+ squarefree (from VI.5).** No F_2-valid argument that is valid *only for squarefree N* proves the uniform bound or any f(W) = o(√(W/log W)). This depends on the squarefree referee's CORRECT verdict.

**The Part V judge's "not affected" list (README_V §3(d)), clause by clause:**

| escape clause | status now |
|---|---|
| odd characteristic | **closed** (VI.4, p = 3 and every p ≡ 3 mod 4; all N) |
| reflection-type conjugation | **closed** (VI.4: i ↦ −i, the Z[i] type) |
| squarefree N "in an essential way" | **closed in char 2** (VI.5); **open in odd char**. In VI.4, N = P_V²/(GḠ) has exponent max_j e_j(v) − min_j e_j(v) ∈ {0,1,2} at each X-level pair. Exponent 2 with marks at π² and π̄² is the generic case, so squarefree N is incompatible with the quadratic-character mechanism. After X = t^L, multiplicities grow like q^{k−2}. The same holds for every specialisation X = T(t) with this V: for α ∈ F_{q²}∖F_q, v ↦ v(α) is F_q-linear and onto, so q^{k−2} factors T − v share each such root. (Proof is mine, two lines; the referee's data 27/31 match.) |
| Lemma P | still true. It now only says that the *constant-difference* mechanism needs characteristic 2; VI.4 has nonconstant cofactor differences of X-degree (n+h)/2 − g. It is no longer evidence for anything about Z[i]. |
| p ≡ 1 mod 4 (A² − dB², d a constant nonsquare) | open. ref6_ternary's half-system lemma kills the formal mechanism. |
| growth-rate arguments | unaffected (§1) |
| number-field Vojta, and J4/(L_8) ⇒ uniform as implications | unaffected (they are implications) |
| archimedean, characteristic-0, Baker inputs | unaffected (below) |

**What a proof over Z could still use**, and whether the characteristic-0 function-field proof (Part IV) uses it:

| feature of Z[i]/Z | used by Part IV (C(t), R(t), Q(i)(t); Langer BMY)? | comment |
|---|---|---|
| archimedean order: the real arc, endpoints, convexity/curvature (Jarník) | **no** (ultrametric degree condition at ∞) | The record's 2/3 proof does not use it either. It uses only the diameter. Every FF counterexample is an ultrametric ball, so a proof that uses linear order or curvature escapes every FF filter. |
| characteristic 0 | **yes, essentially** (Bogomolov–Miyaoka semipositivity inside Langer's orbifold BMY; (IV.2) fails in char 2, ref5_reconcile) | Both FF counterexamples need char p. Over Z, characteristic 0 is present but its arithmetic BMY is at least of Szpiro strength (README_IV). |
| no constant field; conjugation is not a Frobenius; Q(i)/Q is "geometric" and ramified at 2 | **no** (C(t)/R(t) is itself a constant extension, and Part IV works over an algebraically closed constant field) | VI.4 needs coefficient Frobenius twice: conjugation = q-power map, and descent from F_{q²} to F_9. V.1 needs the Artin–Schreier translation. A *geometric* model with inert ∞, such as F_3(t)(√d) with d of even degree and nonsquare leading coefficient, so that finite primes ramify, is **untested by any construction**. That is a genuine open FF escape clause. The BSSW literature uses the ramified model A²+TB². |
| height versus degree (log|a+b| ≤ max + log 2; height counts are polynomial, not q^H) | **yes**: degree as intersection number on a surface, contact orders at the place ∞ | The arithmetic analogue needs Arakelov intersection theory plus arithmetic BMY, which is not available. |
| Baker-type lower bounds for linear forms in logarithms | **no** | Used in the record (ref_baker; Matveev on Hadamard cores), but too weak at the needed scale. Char-p analogues of Baker fail without extra hypotheses (not checked in detail). |
| a derivation (Mason / abc) | char 0 FF has one, and BMY implicitly goes beyond it | Z has none; abc is the conjectural replacement. It has been tested and is vacuous on the fakes (§5). |

Summary for the note: **every algebraic feature of Z[i] that holds verbatim in F_9[t] is now known not to suffice for the uniform bound.** The surviving candidates are:
* characteristic 0;
* archimedean structure (order or curvature, not just an ultrametric ball);
* the absence of a Frobenius on coefficients, or the geometric/ramified nature of Q(i)/Q;
* squarefree or bounded-multiplicity N combined with odd characteristic (a model-level escape, P small).

## 3. (c) Questions A and B in the ternary model

**Formal structure (proved here; checked).**
* At a conjugate pair {v, −v}, a mark has level e_j(v) ∈ {0,1,2}. This gives two layers: A_1 = {e ≥ 1} and A_2 = {e = 2}.
* For any S of five marks with pairwise non-proportional forms, the sign vector (χ(ε_jλ_j(v)))_{j∈S}, with v uniform in V, is uniform on {±1}^5 up to O(q^{−1/2}). This follows from Weil, because a product of distinct non-proportional linear forms is never a square, and it holds **whether or not the forms are independent**.
* Zeros of forms contribute O(1/q). Hence:
  * g_S/n = 1/16 + O(q^{−1/2});
  * every one of the 15 cuts has weight (n/16)(1 + O(q^{−1/2}));
  * deg c_ab = (n+h)/2 − g_ab ≈ h/2 in X-units, so H_Σ/w̄ ≈ 80/q;
  * δ_S = (g_S − h)/(2(n − g_S)) = 1/30 − 8/(15q) + O(q^{−3/2}), so α_S = 4/(1+2δ_S) → 15/4 from above.
* This is ref5_reconcile's binary formula with q in place of m. After X = t^L the corrections are O(k/L) and an L-independent extra gcd, so the limit is the same for L → ∞.
* **Difference from F_2:** there, dependent 5-tuples had only 7 effective cuts (δ → 1/14). Here **all** 5-subtuples are balanced with 15 cuts.

**Computed** (formal, exact enumeration):

| data | δ_S | α_S | max/min cut | H_Σ/w̄ |
|---|---|---|---|---|
| q=27, k=3, Astra's 5 marks | 0.0156 | 3.879 | 1.077 | 2.86 |
| q=27, k=3, ref6_ternary best 5 | 0.0178 | 3.863 | 1.262 | 2.88 |
| q=27, k=5, 56 five-subtuples of 8 forms (5 independent) | [0.0141, 0.0154] | [3.881, 3.891] | ≤ 1.086 | 2.86 |
| **q=243, k=3, all 56 five-subtuples of 8 forms (rank 3)** | **[0.0310, 0.0316]** | **[3.762, 3.767]** | **≤ 1.025** | **0.328** |
| prediction 1/30 − 8/(15q) | q=27: 0.0136; q=243: 0.0311 | | | 80/q: 2.96, 0.329 |

* **Four marks, actual specialised polynomials** (q=27, L=11, D=6981; own Euclid on the saved points):
  * pair gcds 2990, 2990, 2990, 2977, 3003, 3003;
  * triple gcds 988, 1001, 1001, 1001;
  * quadruple gcd 0.
  * Solving for the 7 cuts gives **1001 ×5 and 988 ×2**, balance 1.013. This is 11 × 91 = 1001, the formal X-level weight 91 (all equal), less 13 at two cuts from the extra t^26.
  * Chords deg c_ab = 481, 467, 467, 494, 468, 454, so H_Σ/w̄ = 2.84. At L = 11 the top-coefficient t-degrees (≈ 340) dominate.
  * The four-point δ is formally 0.051 (Liouville value 1/14). Actually it is (3490.5 − 3471)/6981 = 0.0028.
* **Five-mark actual data** (ref6_ternary, L=32, D=568,386): max difference 284,076, so δ = 2·10^{−4}. Small L hides the formal δ. The asymptotic statement needs L → ∞ and q → ∞, which VI.4 provides.

**Conclusion (c).**
* The F_3(t)/F_9(t) analogue of **A_5(δ) is false for every δ < 1/30**, and so are **(L_5°)** and **V_E(5)**: take q > 8/(15(1/30 − δ)), L large, and m → ∞ for the Schwartz–Zippel escape. The J4 steps and Schwartz–Zippel are characteristic-free (ref5_reconcile §(c)). I did not recheck them in characteristic 3, but nothing in them uses the characteristic.
* The balanced Liouville exponent 15/4 is attained by *every* 5-subtuple.
* The **(L_8)** analogue fails asymptotically: balance → 1 and H_Σ/w̄ = O(1/q) for 8-subtuples, by the same Weil argument. This is **CORRECT_WITH_CORRECTIONS**: asymptotic only, star/Plücker dictionary not re-derived in characteristic 3 (it is plausibly easier than in characteristic 2), and no 8-mark instance computed. ref6_ternary's q = 243 formal 8-mark example has D ≈ 7·10^9.
* So over F_3(t), as over F_2(t), the answer to the analogues of Questions A and B is *no*. The number-field questions are untouched.

## 4. (d) Literature

**Searched:**
* function-field short arcs and sectors;
* Bary-Soroker–Smilansky–Wolf (Landau analogue, A²+TB²);
* Bank–Bary-Soroker–Fehm (sums of two squares in short intervals, large q);
* Gorodetsky–Rodgers (variance);
* Rudnick–Waxman (angles of Gaussian primes; the FF model F_q[√−T], super-even characters);
* Yiasemides (UE²+VF², thin annuli);
* Car (norms from F_{q^h}[T]);
* Sedunova (Bombieri–Pila over F_q[T]);
* Entin–Landsberg (Cilleruelo's lcm problem over FF);
* BKR (subspace polynomials) and Ore;
* "Lattice points on small arcs" (arXiv 2107.09991, already in the record: CG disproved for arcs R^α, α > 1/2).

**Findings.**
* The FF literature on sums of two squares is about **averages** (density and variance in short intervals or sectors, mostly q → ∞). It uses the ramified model F_q[√−T], not the constant extension.
* I found **no worst-case short-arc bound and no failure of uniformity** over F_q[t] for fixed q ≡ 3 mod 4. Searches return snippets only, and the search is not exhaustive.
* **Prior art for the mechanism:**
  * Ore's linearised polynomials, and the subspace-polynomial recursion P_{J+F_q u} = P_J^q − P_J(u)^{q−1}P_J;
  * Ben-Sasson–Kopparty–Radhakrishnan (IEEE TIT 2010), "many subspace polynomials with common top coefficients", for the degree drop;
  * the classical identity Π_{s∈(F_q^*)²}(Y − s) = Y^{(q−1)/2} − 1, for the sparsity behind (VI.5).
* Astra's addition of the quadratic-character multiplicities 1+χ, which give exact norms under reflection, and of descent by coefficient Frobenius was not found anywhere.
* "Paley-type" in the record refers to integer fake profiles (round 5), a different object. I found no relevant Dickson or Rédei-function construction.
* Status: **no prior art found for the theorem; ingredients classical.**

## 5. (e) VI.1, Mason, and the integer barrier F1

* VI.1 is Mason–Stothers in characteristic 2. A product R of point ratios lies in U_K, so ord_∞(R′/R) ≥ K+1.
  * In characteristic 2, (π/π̄)′/(π/π̄) = p′/p, and squares have zero derivative. So only primes of **odd** valuation enter.
  * This gives weighted support ≥ K+1 > W/4 for every **binary** word of the parity code. That excludes the twin-based fake profiles (skew double core, fully flipped cores) over F_2[t] at any degree (ref6_code: CORRECT).
* **In characteristic 0** (C(t)), the same computation gives only the **Z-support** condition: the total degree of primes with v ≠ 0 is > W/4. That is Mason proper.
* **The integer analogue under abc** is ref6_code's (abc-Λ): w_0(v) ≥ W/4 − O(δ)w_1(v) − O(log(|c|C)). This is exactly the record's Lang–Waldschmidt support test.
* **The skew double core passes it** for every δ ≥ 0 (record `round8/lw.md` Thm 5.1: w_0 ≥ W/4 + ½log P − 1/16). abc through character relations kills only the Theorem-L (fully flipped, b ≥ 4) families, as LW does.
* The **mod-2** strengthening has no integer analogue. ref6_code (e1) gives an exact Pell family with support weight log 5 against W → ∞.
* **Record check (grep).** abc and Szpiro have been tested against the fakes:
  * `growth/conditional.md` Prop 8.3: every abc-type inequality is vacuous on squarefree near-extremal blocks, and all fakes are squarefree;
  * `round8/weakest.md` (LW, abc, Subspace, Bombieri–Lang barriers);
  * `overview.md` barrier row 15;
  * rounds 18–19: abc on auxiliary pencils of M̄_{0,8}, unit-pencil lane closed by proof;
  * ref4_arith: (IV.2)/BMY is stronger than abc because it uses all m sections at once; arithmetic BMY is at least of Szpiro strength;
  * ref6_code (e3): Szpiro (κ = 3/2) excludes none of the flipped families.
* **What this says about F1.** The FF fakes are excluded by a *characteristic-2* identity, not by "having a derivation". A characteristic-0 derivation, and abc over Z, give only a support condition that the barrier profile satisfies.
  * So the right reading is: the integer fakes survive every *pair/character* consequence of abc. What kills large configurations in characteristic 0 FF is BMY, a global inequality on all m sections at once. Its arithmetic version is unproved and at least of Szpiro strength.
  * "An abc-type input would exclude the fakes" is **false for abc on character relations** (proved at profile level, record + ref6_code). It is **open** for a hypothetical arithmetic orbifold BMY.
  * Even over F_2[t], VI.1 does **not** settle growth: low-rank Walsh two-core profiles pass at W = O(m log m) (ref6_code (c)).

## 6. (f) Recommended wording for the expert note (recommendations only; the record is read-only)

1. New bullet after Part V in §5 item 12:
   > ***Part VI: unbounded over the fixed field F_3 with A²+B²*** (`round31/README_VI.md`). **Theorem** (colleague, VI.4; ternary referee CORRECT). Let p ≡ 3 mod 4 and conjugate F_{p²}[t] by i ↦ −i (the norm form A²+B² of the constant extension F_{p²}(t)/F_p(t), the Z[i] type). For every m there are m distinct monic z_j ∈ F_{p²}[t] of one degree D, with one norm N ∈ F_p[t] and deg(z_i − z_j) ≤ D/2. Here deg N = exp(O(m log m)), and M_p(D,0) ≥ (½log_2 p − o(1)) log D/loglog D for all large D. Mechanism: quadratic-character multiplicities 1+χ(ελ_j(v)) on the roots of an F_q-subspace polynomial (q = p^r, r odd), descent by coefficient Frobenius, and a large common divisor. N carries squared primes and multiplicities that grow with m. The squarefree variant is open in odd characteristic. The transferred growth bound W ≥ 3M log_p M − O(M) (no ramified place) leaves an exponential gap. **Squarefree over F_2** (VI.5; referee CORRECT): squarefree N of degree 8m(m−1)(log_2 m − 1) carry m = 2^r points.
2. Replace the sentence "Not affected: … arguments needing odd characteristic, a reflection-type conjugation, or squarefree N in an essential way" by:
   > "Not affected: growth-rate arguments; number-field Vojta and the implications A_5 ⇒ uniform and (L_8) ⇒ uniform; arguments that use characteristic 0, archimedean order or curvature, the absence of a Frobenius on coefficients, or squarefree (bounded-multiplicity) N in odd characteristic. Odd characteristic and reflection conjugation alone no longer suffice (Part VI), nor does squarefree N alone (Part VI, characteristic 2)."
   Add the F_3 filter of §2. It is weaker than F3+ for growth: it excludes only f(W) = o(log W/loglog W).
3. Replace "*Not evidence against the conjecture over Z.* The mechanism needs many pairwise cofactor differences that are constant … Odd characteristic is open" by:
   > "*Not evidence against the conjecture over Z.* Both function-field mechanisms need coefficient Frobenius (Artin–Schreier translation in characteristic 2; q-linearised subspace polynomials plus Frobenius descent in characteristic p ≡ 3 mod 4). The characteristic-0 model, which shares the reflection conjugation, the anisotropy at ∞ and exact norms with Z[i], is bounded (Part IV). Lemma P now explains only why the characteristic-2 mechanism does not port. Over F_9(t) the analogues of A_5(δ < 1/30), (L_5°) and V_E(5) are false, with every 5-subtuple balanced at exponent → 15/4. (L_8) fails asymptotically, pending the dictionary."
4. Permitted: "**the uniform statement is false over F_p(t), p ≡ 3 mod 4, for the Gaussian norm form of the constant extension**" (a theorem; ternary referee CORRECT).
   Not permitted:
   * "false in positive characteristic" (p ≡ 1 mod 4 and geometric quadratic extensions are untested);
   * "false for sums of two squares over function fields" (BSSW's A²+TB² model is untested);
   * "evidence against the conjecture over Z";
   * "the growth rate over F_q is log W/loglog W" (only a lower bound for sup M is known).
5. The §5 item 10 sentence on the transfer: add "every algebraic feature of Z[i] that holds verbatim in F_9[t] is insufficient for the uniform bound".
6. On abc (barrier 6), refine to: "abc via character relations is the LW support test. It excludes the Theorem-L families, not the skew double core. Its char-2 FF strengthening (VI.1) has no integer analogue."

## 7. Novelty and significance

* **Novelty against the record.**
  * The F_3 transfer constant is new but routine.
  * Weil equidistribution makes *all* ternary 5-subtuples balanced, so the FF failure of A_5 now holds without a "rank-5 tuples" proviso. This is new to the record, though a direct consequence of known facts.
  * The structural squaredness of N in VI.4 (squarefree is impossible for this mechanism, and multiplicities are forced for every polynomial specialisation in this V) is new to the record.
  * The ledger of Z-features against what Part IV uses is a consolidation.
  * VI.4 itself is the main novelty (ternary block).
* **Significance for the goal: 2/10.**
  * It removes the last algebraic escape clauses for uniform proofs: odd characteristic, reflection, and (in characteristic 2) squarefree N.
  * It does not touch the integer growth rate (2/3+ε) log R/loglog R, the constant 3, or Questions A and B over Z.
  * It is not evidence about Z.

Sources:
* [Bary-Soroker–Smilansky–Wolf, arXiv 1504.06809](https://arxiv.org/pdf/1504.06809)
* [Bank–Bary-Soroker–Fehm, arXiv 1509.02013](https://arxiv.org/pdf/1509.02013)
* [Gorodetsky, arXiv 1603.02890](https://arxiv.org/pdf/1603.02890)
* [Rudnick–Waxman, arXiv 1705.07498](https://arxiv.org/pdf/1705.07498)
* [Yiasemides, arXiv 2309.01290](https://arxiv.org/pdf/2309.01290)
* [Variance of sums of two squares over intervals in F_q[T], arXiv 2204.04459](https://arxiv.org/pdf/2204.04459)
* [Car, Portugaliae Math.](https://www.maths.tcd.ie/EMIS/journals/PM/51f4/8.html)
* [Sedunova, arXiv 1506.08757](https://arxiv.org/pdf/1506.08757)
* [Entin–Landsberg, arXiv 2310.04164](https://arxiv.org/pdf/2310.04164)
* [BKR, Subspace polynomials and limits to list decoding](https://www.Math.Toronto.Edu/swastik/rsld.pdf)
* [Lattice points on small arcs, arXiv 2107.09991](https://arxiv.org/pdf/2107.09991)
* [Cilleruelo–Granville, Close lattice points on circles](http://matematicas.uam.es/~franciscojavier.cilleruelo/Papers/closelattice.pdf)
