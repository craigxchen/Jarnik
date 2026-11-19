# Referee report, round 31, Part V, block "lit": literature and meaning

Object: Part V of `astra_free_V/answer.md` (from line 1138). Over the fixed field F_2 with the anisotropic norm A^2+AB+B^2 of F_4(t)/F_2(t), it gives m = 2^r points at c = 0 on a parent of degree 2m(m-1)(log_2 m - 1).

Scripts: `round31/ref5_lit/` (`check.py`, `sqf.py`). All code is mine: F_2[t] as bit-vector integers and F_4[t] as pairs A + B·θ.

The construction itself belongs to other blocks (`ref5_construct.md`, `ref5_reconcile.md`). I re-checked it only enough to talk about it responsibly:
* `check.py` rebuilds z_H for m = 4, 8, 16. Every norm is equal and in F_2[t], the points are distinct, deg N = (Q-2)L = 2D, and the maximum pair-difference degree is exactly D/2 (6, 56, 360). This agrees with Astra's table.
* `sqf.py` shows that N is **not squarefree** (deg gcd(N, N') = 8 at m = 4 and 36 at m = 8). This is consistent with Astra's phrase "the model allowing prime powers". The integer conjecture allows every N, so this is legitimate, but the record should state it.

**Access caveat.** arxiv.org and math.toronto.edu do not resolve from this session. Every literature statement below comes from WebSearch snippets (abstracts, MathSciNet-style listings, secondary citations) or from memory, and is marked accordingly. None was checked against an original text.

## Verdict table

| # | Question | Verdict |
|---|---|---|
| a1 | Is the function-field short-arc problem (worst case, fixed q, uniform in N) in the literature? | **Not found.** What exists is statistical: averages and variances over short intervals or sectors, mostly as q → ∞. No uniform worst-case bound and no counterexample. Part V's statement appears new. |
| a2 | Is a positive-characteristic failure of it known? | **Not as such.** The ingredients are classical: subspace (linearised) polynomials whose top coefficients agree (Ben-Sasson–Kopparty–Radhakrishnan), and Artin–Schreier conjugation. The failure belongs to a well-documented genre: Mahler's counterexample to Roth, Davenport–Stothers failing in char p, Easton's and Lang's BMY violations, Gasbarri's isotrivial Lang–Vojta counterexamples. |
| a3 | Is isotriviality the reason? | **No, not by itself.** The integer problem and the char-0 function-field model are equally "isotrivial": one fixed norm-one torus serves every N. Part IV proves a bound in the isotrivial char-0 case. The operative mechanism is additive (Artin–Schreier plus F_2-linear subspace polynomials). The classifying point is separable (`ref5_reconcile.md` (d)), so this is not the Frobenius-twist failure either. |
| b | Integer problem: warnings, and evidence beyond heuristics | **No published warning found** that the uniform conjecture may be false. Evidence consists of theorems at smaller scales (Jarník; Cilleruelo–Córdoba below R^{1/2}), almost-all results (Bourgain–Rudnick), and conditional implications (Vojta ⇒ CG Conj. 14 ⇒ critical case, in the record). The record's own fakes are the only "warning", and they are fakes. |
| c | What Part V suggests about Z | It is a **filter, not evidence**. Every load-bearing feature lacks an integer counterpart (list in §(c)). The features it shares with Z (fixed small residue fields, an anisotropic norm inert at the infinite place, the pigeonhole bound) are what make it a sharp filter. The reading "no characteristic-free proof" is right only with two qualifications, in §(c). |
| d | Wording for the record | Proposed text in §(d). |

## (a) Literature

### (a1) Function-field sums of two squares, norms and divisors in short intervals

All of the following are statistical: densities, means or variances, mostly in the large-q limit. None bounds the worst-case number of representations of one fixed N in one short "arc".

* **Bary-Soroker–Smilansky–Wolf**, "On the function field analogue of Landau's theorem on sums of squares" (arXiv 1504.06809; FFA 2016). Their "sum of two squares" is A^2 + T B^2 in F_q[T], q odd, which is a ramified extension. They count in the large-degree and large-q limits. (Snippet.)
* **Bank–Bary-Soroker–Fehm**, "Sums of two squares in short intervals in polynomial rings over finite fields" (arXiv 1509.02013; Amer. J. Math. 2018). The density in intervals deg(f - f_0) ≤ εn is 4^{-n} binom(2n,n) as q → ∞. (Snippet.)
* **Gorodetsky**, "A polynomial analogue of Landau's theorem and related problems" (arXiv 1603.02890). This covers the limit q^n → ∞. (Snippet.)
* **Gorodetsky**, arXiv 1810.00483. He defines the F_q[T] Hooley Δ-function and shows that its short-interval means match the global mean for h/n > 1/2. A later paper calls this the only known short-interval result for Δ over F_q[T]. (Snippet.)
  * This is the function-field version of the "divisors in short intervals" problem for the **split** form. It is again an average, not a worst case.
* **Rudnick–Waxman**, "Angles of Gaussian primes" (arXiv 1705.07498), with **Katz**, "Witt vectors and a question of Rudnick and Waxman" (IMRN 2016), and Perret-Gentil (ANT 2020).
  * The function-field "Gaussian integers" are F_q[S] with S^2 = -T, and the norm is f(S) f(-S).
  * Sectors are cosets of higher unit groups of the norm-one torus, detected by super-even characters. These are exactly the record's U_k sectors at infinity in another model.
  * The results are equidistribution and variance as q → ∞. This is the **closest existing object** to the record's model. It appears in the record already as lit18 `primes.md` ("robust, capped by the graft barrier").
* **Yiasemides**, arXiv 2309.01290. It treats the mean and variance of representations by U E^2 + V F^2 in short intervals. It was already cited by `ref5_construct.md`.
* **Chebotarev in short intervals** (arXiv 1810.06201) is for *geometric* extensions. The constant extension F_4(t)/F_2(t) is not geometric, so it does not apply (contrary to a search-engine summary).
* **Searches that returned nothing on target**: "lattice points on circles function field", "short arcs F_q[t]", and "x^2+xy+y^2 characteristic 2 short intervals". The only on-topic hit in char 2 is a paper on representing elements by x^2 + xy + a y^2 over F_q(t) in characteristic 2 (arXiv 2203.04068), which is algorithmic and irrelevant here.
* Round 29 (`kill_ff.md`) and Part IV's `ref4_lit.md` also searched and found no function-field Cilleruelo–Granville problem. **I confirm that.**

**Model mismatch, worth one line in the record.** The literature's "Gaussian integers" over F_q[T] are ramified (S^2 = -T, q odd). The record's model and Part V use a constant extension, inert everywhere.
* In characteristic 2, x^2 + y^2 = (x+y)^2 is degenerate, and A^2 + T B^2 is an inseparable extension.
* So the only separable quadratic "circle" over F_2(t) with an inert infinite place is the Artin–Schreier one.
* Hence in characteristic 2 the circle problem has no Kummer model at all. **Part V lives in the case where the analogue of Z[i] is least like Z[i].**

### (a2) The ingredients are classical

* **Subspace polynomials with agreeing top coefficients.** Ben-Sasson, Kopparty and Radhakrishnan, "Subspace polynomials and limits to list decoding of Reed–Solomon codes" (FOCS 2006; IEEE Trans. IT 56 (2010) 113–120).
  * Snippet: "a received word that agrees with a super-polynomial number of distinct degree K polynomials". The technique is subspace polynomials, which are sparse linearised polynomials.
  * This is the same phenomenon as (V.1). Many polynomials ∏_{v∈H}(X+v), indexed by subspaces, have pairwise differences of much smaller degree, because P_H(X+b) = P_H(X) + P_H(b) and only X^{2^i} coefficients survive.
  * The record (README_IV line 265, ref4_trace) already calls the trace seed "elementary and classical: the top-coefficient sparsity of linearised polynomials". It did not name this paper.
  * **Recommend citing BKR as the prior art for the mechanism.** Only the specialisation to equal norms via Artin–Schreier conjugation is Astra's.
* **Prouhet–Tarry–Escott reading** (my observation, elementary). Two products ∏_{v∈S}(Y+v) and ∏_{v∈S'}(Y+v) differ in Y-degree < |S| - k exactly when the elementary symmetric functions e_1, ..., e_k of S and S' agree. This is the polynomial form of PTE.
  * In characteristic 2, any subspace H and its translate b + H form an *ideal* solution of every size, since e_j agree for all j < |H|.
  * Part V needs the weaker, multigrade statement for the m sides H \ {0}, with agreement up to e_{Q/4}.
  * Over Z, ideal PTE solutions are known only for small sizes (n ≤ 10 and n = 12, from memory; surveys arXiv 2304.11254 and 2506.11429 were found but not read). A 2026 snippet (arXiv 2610.05813) says a general PTE theory over finite fields "remains to be developed".
  * The record's integer polynomial templates are of exactly this product-of-shifts type: the CG six-point family ∏_{j≤4}(a + j + iσ_j) and the cube family. So the integer counterpart of Part V's engine is a hard open combinatorial problem, not a known construction.

### (a3) Known characteristic-p failures of Diophantine statements that hold in characteristic 0

Part V's statement belongs to this genre:
* **Roth.** It holds over char-0 function fields (Uchiyama, Osgood; from memory). Over F_q(t) it fails. Mahler's α = Σ t^{-p^n}, which satisfies the Artin–Schreier-type equation α^p - α = -1/t (a snippet gives β^q - β - t^{-1} = 0), has exponent equal to the Liouville bound (Kim–Thakur–Voloch, de Mathan, Lasjaunias survey; snippets).
  * **This is the closest parallel.** The record reformulates CG as Roth-type finiteness at Dirichlet exponent 1/2 (overview Prop. 3.2), and both failures are Artin–Schreier/additive.
* **Mason–Stothers / abc.** These hold in characteristic p only for triples not all with vanishing derivative (Wikipedia/Lean snippet). Counterexamples come from p-th powers.
* **Szpiro.** Pesenti–Szpiro (Compositio 120 (2000)) prove the function-field discriminant inequality in characteristic p with a factor p^d, the inseparability degree. Schütt–Schweizer (math/0608427) use it to show that **Davenport–Stothers fails in characteristic p > 3** and to correct it. (Snippets.)
  * Davenport's deg(f^3 - g^2) ≥ deg f/2 + 1 is a "differences cannot be too small" statement of the same flavour as the cluster condition.
* **BMY.** Lang (1983) and Easton, PAMS 136 (2008) 2271–2278, "Surfaces violating BMY in positive characteristic", give covers branched over line configurations that exist only in characteristic p. Urzúa (Duke 2017) shows Chern slopes are dense in (3, ∞). (Snippets.)
  * This matters because Part IV's bound is exactly orbifold BMY, and `ref5_reconcile.md` (a) finds that the only step of (IV.2) failing for the Part V data is Langer's theorem.
  * So **Part V is to (IV.2) what Easton's surfaces are to BMY**: a configuration that exists only in characteristic p and violates the char-0 Chern inequality.
* **Lang–Vojta.** Gasbarri ("On some differences between number fields and function fields", arXiv 1606.06517; snippet) builds positive-characteristic counterexamples to Lang/Vojta-like statements from isotrivial varieties, inseparable covers and Frobenius. Voloch's survey says naive transposition of Vojta to characteristic p "runs into pathological behaviour". arXiv 0910.3417 gives an isotrivial surface over F_q(t) with a large set of integral points.

**Is isotriviality the reason?** Not by itself:
1. The norm-one torus T = R^1_{F_4/F_2} G_m is constant, so the conics x^2 + xy + y^2 = N are torsors under one fixed torus. The same is true of x^2 + y^2 = N over Z: "one fixed variety serves every N (isotriviality)" (overview Prop. 3.2). It is also true of the char-0 model over R(t)(i) or Q(i)(t), where Part IV *proves* M(D,0) ≤ 24. So isotriviality is shared by the bounded char-0 case and the integer case, and cannot separate them.
2. What is special in characteristic p, for constant (isotrivial) targets, is that the ring F_2[t] has the additive Frobenius. This gives additive polynomials, finite additive subgroups of bounded height, and Artin–Schreier Galois action.
3. `ref5_reconcile.md` (d) reports that the classifying point of the Part V cluster is separable. So the "Frobenius iterate / inseparability" failure, which abc and Szpiro repair by a p^d factor, is **not** what happens.
   * Consequence: **a separability-corrected uniform statement in the style of Pesenti–Szpiro is not rescued by Part V's data.** The correction would have to exclude additive configurations, as in Easton's case, not inseparable points.
   * This is my reading and it depends on the reconcile block's separability check, which I did not redo.

## (b) The integer problem: warnings and evidence

Sources: record plus snippets; CG 2007 = arXiv math/0608109; CG 2009, CJM 61.

* **The literature conjecture.** Cilleruelo–Granville, Conj. 14 in the 2007 paper (cited as such by Bourgain–Rudnick, snippet): arcs of length R^{1-ε} contain at most B_ε points. The critical C√R case is the record's sharpening. It follows from Conj. 14 with ε = 1/4 (expert_note, absolute form of the A-branch).
* **Theorems.**
  * Jarník: arcs of length < c R^{1/3} contain at most 2 points.
  * Cilleruelo–Córdoba: arcs of length R^{1/2 - 1/(4⌊k/2⌋+2)} contain at most k points.
  * CG 2009 classify 4-tuples on arcs of order R^{1/3} (snippet).
  * O(log R) on arcs of length R^{1/4}·O(1) (Bourgain–Rudnick, Lemma 2.1, snippet) and ≪_C log R at √R.
  * The record's (2/3+ε) log R / log log R.
* **Almost-all.** Bourgain–Rudnick: for all but O(X^{1-ε/3}) of m ≤ X, distinct points are separated by R^{1-ε} (snippet). This is density-zero-blind, which the record notes repeatedly.
* **Conditional.** In the record, Vojta's conjecture implies CG Conj. 14 (overview Thm 4.2) and hence the critical case. Cilleruelo–Granville themselves connect related questions (squares in arithmetic progressions, Rudin's conjecture) to Bombieri–Lang (math/0608109 abstract).
* **Warnings.** I found **no published suggestion that the uniform conjecture is false**. The only counter-evidence anywhere is the record's own fakes: integer-matrix models satisfying every divisibility and the norms only approximately (item 12). The integer constructions have bounded size; the record's best lower bounds come from the golden and cube families.
* **Reading.** The conjecture rests on (i) theorems below the critical exponent, (ii) almost-all results, and (iii) implication from Vojta. Item (iii) is precisely the type of statement known to fail naively in characteristic p. So Part V is *what one should expect* of a consequence of Vojta, and is not evidence against it over Q.

## (c) What a fixed-field F_2(t) counterexample suggests about Z

### Features of Part V with no integer counterpart (load-bearing)

1. **Additive polynomials.** P_H(X+b) = P_H(X) + P_H(b) is what makes every pair difference factor as P_J(1)·P_J(Y)/Y (V.1). It comes from the additive Frobenius of F_2[t]. Z has no ring endomorphism of this kind; the Frobenius lift exists only p-adically, as Buium's δ, which the record already examined and found point-level and faked.
2. **Large finite additive subgroups of bounded height.** V = F_2[t]_{<d} has 2^d elements of height < d and is closed under addition, and hyperplanes H give Q/2 elements of the same kind. Z has no nontrivial finite additive subgroups, and integer sets of size n have sumsets of size at least 2n - 1.
3. **Conjugation acting as translation** (Artin–Schreier): conj(Y) = Y + 1. Astra proves that a translation needs characteristic 2. In Z[i] and Z[ω], conjugation is the reflection Y ↦ -Y or Y ↦ -1-Y. That reflection is exactly why integer templates need sign vectors with vanishing moments.
4. **Ideal PTE solutions of every size for free** (consequence of 1–2). Over Z these are open beyond small sizes.
5. **Infinite supply of evaluation points of low height** (Y = t^L + θ for every L). Over Z, specialising a one-parameter template is controlled by Part IV: templates with bounded C carry at most 24 points. So an integer "Part V" would have to be non-template.

### Features shared with Z (why the filter is sharp)

6. **Fixed, small residue fields.** F_2 here, F_p over Z. The record's pigeonhole lower bound transfers exactly, as W ≥ 3M log_2 M - O(M) (`ref5_reconcile.md` (b)), and Part V satisfies it. This is the feature Part IV's trace seed lacked, since there q grew with m.
7. **An anisotropic norm with the infinite place inert.** This is the analogue of x^2 + y^2 being definite. The trace seed lacked this too, being split.
8. **Exact norm and exact divisibilities.** A genuine cluster, not a fake.
9. **The same isotrivial torus structure** (see (a3)).

### The right reading

"The uniform conjecture has no function-field-over-finite-field analogue in characteristic 2, so no proof can be characteristic-free" is **correct with two qualifications**:

* **(i) Scope.** It refutes the char-2 analogue only for **p = 2**, through an Artin–Schreier circle that has no Kummer-type sibling. Astra shows the translation trick needs characteristic 2, and ref4_trace found the odd-p subspace analogue misses by (p-2)^2/(4p^2).
  * So the statement "no proof can be valid over every F_p(t)" is proved.
  * The statement "no proof can be valid over F_p(t) for p odd, for example the ternary F_9/F_3 circle" is **open**.
  * A proof using properties that hold over Z and over F_p(t) for odd p but fail for p = 2 is not excluded. Example: the fact that 2 is the only ramified prime of Z[i], or the existence of a Kummer description.
* **(ii) Only uniform statements are excluded.** Part V's parent degree is Θ(m^2 log m), far above the growth scale m log m. So growth-rate arguments are untouched (V.2 and reconcile (e)). Any argument proving o(log R / log log R) is not filtered by Part V.

**Comparison with abc/Szpiro/Vojta.**
* Those are theorems over char-0 function fields, and theorems in characteristic p after a **separability correction** (p^d inseparability factors; non-p-th-power hypotheses).
* Uniform CG is now in the same position for the char-0 half: a theorem over C(t) and Q(i)(t) (Part IV).
* For the char-p half the analogy is **weaker**. No separability correction is known to restore it, and Part V's point is separable.
* The better analogues are **BMY and Davenport–Stothers in characteristic p**, which fail for configurations that exist only in characteristic p (Easton, Schütt–Schweizer), and **Roth**, which fails through an Artin–Schreier element (Mahler).
* In none of these cases did the char-p failure cast doubt on the number-field statement. Roth over Q is a theorem. BMY and abc over number fields are believed (abc: Mochizuki's claimed proof is disputed; not relevant here).

**Honest summary for Z.** Part V is a **barrier result for methods, with zero evidential weight on the truth of the integer conjecture.**
* It removes, for the *uniform* bound, every argument whose steps all hold verbatim for F_4(t)/F_2(t): geometry of numbers, pigeonhole over small residue fields, Hodge index, Riemann–Roch, characteristic-free Chern inequalities (filters F3/F3+).
* What remains is arithmetic input that sees characteristic 0. Examples: archimedean convexity and real analysis of the arc, Baker-type transcendence, Faltings/Vojta-type Diophantine geometry, or 2 ∈ Z^× ⊗ Q.

## (d) Recommended wording for the record

> *Part V (char 2, fixed field; proved, refereed).* Over F_2(t) with the anisotropic norm A^2+AB+B^2 of the constant extension F_4(t)/F_2(t), the function-field cluster model has unbounded clusters at c = 0. For m = 2^r there are m points of degree D = m(m-1)(log_2 m - 1) with all pair differences of degree ≤ D/2 (equality attained); N is not squarefree. Hence the char-0 theorem 8 ≤ sup_D M(D,0) ≤ 24 (Part IV) has no analogue over F_2(t). The mechanism is additive: subspace polynomials (top-coefficient sparsity of linearised polynomials, as in Ben-Sasson–Kopparty–Radhakrishnan's limits to list decoding of Reed–Solomon codes) combined with Artin–Schreier conjugation, which acts as translation and is possible only in characteristic 2. It is not an isotriviality effect (the integer and char-0 problems are equally isotrivial) and not a Frobenius-twist effect (the point is separable). It is of the same type as Mahler's char-p counterexample to Roth and Easton's char-p violations of BMY. Parent degree is Θ(m^2 log m), so the growth scale m log m is not reached within the construction (Ω(s^2 log s) for every subfamily, evaluation and common division of the same V).
>
> *Meaning.* (1) No proof of the uniform bound can consist only of steps valid over every F_p(t) (Filter F3+: fixed small residue fields and an anisotropic norm do not suffice). (2) Odd p (e.g. F_9/F_3) is open; a proof may use features of Z[i] that fail for Artin–Schreier circles. (3) Growth-rate arguments are untouched. (4) This is not evidence against the integer conjecture: like Roth and BMY, the statement is a char-0 phenomenon whose char-p analogue fails for additive reasons. **Nothing changes for the integer problem.**

Avoid these wordings:
* "the conjecture is false in function fields" (it is a theorem in characteristic 0);
* "the conjecture fails in positive characteristic" (proved only for p = 2);
* "isotriviality explains it".

## Novelty and significance

* **Novelty.** The statement is new as far as I can find; grep of the record and web search found no fixed-field anisotropic counterexample. The method is classical (BKR subspace polynomials, Artin–Schreier). The record should cite BKR as prior art for the mechanism, which neither Astra nor the earlier referees did.
* **Significance.**
  * For the function-field model: high. It completes the picture: bounded in characteristic 0, unbounded over F_2.
  * For the integer goal: low (filter only). It sharpens what any uniform proof must use. It does not touch the growth rate. It says nothing about truth over Z.
