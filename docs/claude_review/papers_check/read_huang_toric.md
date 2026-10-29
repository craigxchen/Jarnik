# Full-text reading: Z. Huang, "Rational approximations on toric varieties" (arXiv 1809.02001v6)

Read in full: all 64 pages of the uploaded PDF (`84b58274-1809.02001v6.pdf`), pp. 1-64, including the references and
the index of notation. Page numbers below are the numbers printed on the PDF pages (they coincide with the PDF page
index). "Quote" means verbatim from the PDF up to typesetting (LaTeX symbols rendered in plain text). Anything marked
**[inference]** is mine, not Huang's. Scripts: `huang_hyp_check.py` -> `huang_hyp_check.out`, `huang_lines_check.py` ->
`huang_lines_check.out`, both in this directory.

## 1. Bibliographic data as printed in the PDF

* Title (p. 1): **RATIONAL APPROXIMATIONS ON TORIC VARIETIES**.
* Author (p. 1): **Zhizhong Huang**, Institute of Mathematics, Academy of Mathematics and Systems Science, Chinese
  Academy of Sciences, 100190 Beijing, China; e-mail zhizhong.huang@yahoo.com.
* arXiv stamp (p. 1, margin): **arXiv:1809.02001v6 [math.NT] 4 Aug 2025**.
* **No journal reference is printed anywhere in the PDF** (no "Algebra Number Theory", volume, year or page range).
  Two footnotes show post-publication edits: footnote 10 (p. 28) "Edited on 06.12.2021." and footnote 18 (p. 55)
  "Edited after publication on 06.12.2021." So the record's "Algebra Number Theory 15 (2021) 461-512" is **not
  confirmed by this source** (it is neither contradicted); it remains a snippet-level citation.
* Acknowledgements (p. 9): the paper "grew out of part of my Ph.D. thesis realised at Université Grenoble Alpes";
  thanks to E. Peyre, D. McKinnon, M. Brion (idea of primitive collections), an anonymous referee; support from ANR
  GARDIO, a Riemann fellowship, DFG grant DE 1646/4-2.
* Contents (p. 1): 1 Introduction (2); 2 Geometric preliminaries (10); 3 Universal torsors and Cox coordinates (21);
  4 Approximation constants (26); 5 The canonical embedding of a number field (33); 6 Determination of α-constants and
  locally accumulating subvarieties (34); 7 Toric varieties with Picard number 2 (48); References (60); Index (63).
* Abstract (p. 1), quoted: "Using the universal torsor method due to Salberger, we study the approximation of a
  general fixed point by rational points on split toric varieties. We prove that under certain geometric hypothesis the
  best approximations (in the sense of McKinnon-Roth's work) can be achieved on rational curves passing through the
  fixed point of minimal degree, confirming a conjecture of McKinnon. These curves are also minimal in the sense of
  deformation theory, and they correspond, according to Batyrev's terminology, to the centred primitive collections of
  the structural fan."

## 2. Definitions (with normalisations)

**Field, places, absolute values** (§1.5, p. 9). "We fix throughout this paper a number field K." For finite
ν | p: |x|_ν = |N_{K_ν/Q_p}(x)|_p. For ν archimedean: |·|_ν = |·| if ν real, |·|_ν = |·|^2 if ν complex. So heights
are **relative to K** (not absolute) and the distance at a complex place is the **square** of the usual one.
**[inference]** Both scale by the same factor under base change, so the exponent γ in (2) is unchanged; e.g. for
K = Q(i), H_{Q(i)} = H_Q^2 on Q-points and d_ν = dist^2, hence d_ν^γ H_{Q(i)} = (dist^γ H_Q)^2.

**Setting of §4** (p. 26). "Let X be a projective variety over K ... We fix a rational point Q ∈ X(K̄), ν ∈ M_K, and L a
line bundle, to which we associate a height function H_L : X(K) → R_{>0}." Distance (33): with an embedding
j : X → P^N_K, F ⊃ K with Q ∈ X(F), j(Q) = (x_1, ..., x_N) in an affine chart,
d_ν(·, Q) = min(1, max_{1≤i≤N} |X_i(·) − x_i|_ν̄). "Any other distance functions arising from different embeddings are
equivalent" [MR15, Prop. 2.4, Lemma 2.5]; −log d_ν(·, Q) is a local Weil height for the exceptional divisor of Bl_Q X.

**Definition 4.1** (McKinnon-Roth [MR15, Defs. 2.7-2.9]; pp. 26-27), quoted: "For any subvariety Y ⊂ X, we define the
(best) approximation constant α_{L,ν}(Q, Y) (depending on L and ν) to be the infimum of A_{L,ν}(Q, Y), where
A_{L,ν}(Q, Y) is the following set
{γ > 0 : ∃C > 0, ∃(P_i) ∈ (Y(K) \ {Q})^N, d_ν(P_i, Q) → 0, and for every i, d_ν(P_i, Q)^γ H_L(P_i) < C}."
Points are **K-rational** (P_i ∈ Y(K)); the target Q may be algebraic in general, but all main theorems take
Q ∈ T(K).

**Remarks 4.2** (p. 27): (1) α is independent of the choice of d_ν and H_L. (2) "If X is toric, then height functions
are equivalent under torus action, and the same holds for distance functions as well. So α_{L,ν}(·, X) is constant on
T(K)."

**Examples 4.3** (p. 27): α_{O(1),∞}(θ, P^1) = 1/μ(θ) = 1 if θ ∈ P^1(Q), 1/2 otherwise (θ ∈ P^1(Q̄) ∩ P^1(R)); for
abelian varieties α = ∞.

**Proposition 4.4** ([MR15, Lemma 2.13, Prop. 2.14]; p. 27): (1) Q ∈ P^n(K) ⇒ α_{O(1),ν}(Q, P^n) = 1 for every ν;
(2) α_{mL,ν}(Q, Y) = m α_{L,ν}(Q, Y); (3) Y_1 ⊂ Y_2 ⇒ α_{L,ν}(Q, Y_1) ≥ α_{L,ν}(Q, Y_2). And (p. 28): for a smooth
rational curve l through Q, α_{L,ν}(Q, X) ≤ α_{L,ν}(Q, l) "which is equal to deg_L l if Q ∈ X(K)".

**Proposition 4.5** (p. 28), quoted: "For any closed subvariety Z of X such that Y ⊄ Z and inf_{P∈Z(K)} d_ν(P, Q) > 0,
consider the set B^Z_{L,ν}(Q, Y) = {γ ≥ 0 : ∃C > 0, d_ν(P, Q)^γ H_L(P) ≥ C, for every P ∈ (Y \ Z)(K) \ {Q}} and let
b^Z_{L,ν}(Q, Y) = sup B^Z_{L,ν}(Q, Y). Then (1) α_{L,ν}(Q, Y) ≥ b^Z_{L,ν}(Q, Y). (2) Assume moreover that L verifies
the Northcott property ... on Y \ Z ..., then α_{L,ν}(Q, Y) = b^Z_{L,ν}(Q, Y)." Remark 4.6 (p. 29): big L has Northcott on
a dense open U; for X toric, Q ∈ T, Y meeting T, α(Q, Y) = b^{X\T}(Q, Y).

**Definition 4.7** (Pagelot [Pag08]; p. 29), quoted: "we define the essential constant of Q (with respect to Y) to be
α_ess(Q, Y) = sup_{V⊂Y} α(Q, V), where V ranges over all Zariski open dense subvarieties of Y such that
inf_{P∈V(K)} d_ν(P, Q) = 0." α_ess(Q) := α_ess(Q, X). Footnote 12: in [Hua19, Déf. 2.3] and [Hua20, Déf. 2.2] the sup
was over dense constructible subsets; equivalent. **[inference]** So α_ess is exactly the record's α^ess: the best
exponent of Zariski-dense sequences of K-points (sup over removing proper closed sets).

**Definition 4.8** (p. 30), for Z a proper closed subvariety: (1) α(Q, Z) = α(Q, X): best approximations "can be
achieved" on Z; (2) α_ess(Q, Z) < α_ess(Q): Z is "locally accumulating"; (3) α(Q, X) = α_ess(Q, Z) < α(Q, X \ Z): best
approximations "are properly achieved" on Z; (4) α(Q, Z) = α_ess(Q): "the generic best approximations (of Q) can be
achieved on Z". Remarks 4.9 (pp. 30-31) unpack these.

**Toric/geometric definitions.** Fan △ in N ≅ Z^n, T = Spec K[M] ≅ G^n_{m,K} (p. 10: "We suppose throughout this paper
that X is projective and smooth"); relations Σ a_ρ ρ = 0, positive relations (Def. 2.1, p. 11); curve classes meeting T
↔ positive relations (Thm 2.2, pp. 11-12); L-degree of a relation deg_{O(D)} P = −Σ a_ρ φ_D(ρ) (Def. 2.4, p. 13);
**primitive / centred primitive collection** (Batyrev, Def. 2.5, p. 14: centred if Σ_{ρ∈I} ρ = 0); d-free, free, very
free curves (Def. 2.9, p. 16); Theorem 2.11 (p. 17): "A positive relation P represents very free rational curves if and
only if Vect_Q{ρ : ρ ∈ P(1)} = N_Q."

**Hypothesis (\*)** (p. 5), quoted: "The cone of pseudo-effective divisors Eff‾(X) is simplicial." **Lemma 6.3**
(p. 35): (\*) is equivalent to "(\*\*) there exists σ_0 ∈ △_max such that all generators in △(1) \ σ_0(1) are linear
combinations of those in σ_0(1) with negative integer coefficients."

**β** (35) (p. 34): β = min_{P centred primitive} deg_L P.

**Universal torsor and heights** (§3). Thm 3.1 (p. 22): the quasi-affine 𝔗 ⊂ A^{△(1)} is a universal torsor over X
under the split Néron-Severi torus T_NS (footnote 7: "Since T_NS is split, the unicity follows from the Hilbert 90").
Thm 3.2 (p. 23): every K-point lifts to an O_K-point of a twisted torsor 𝔗̃_c, c ∈ C^r (C = ideal class
representatives), with coordinates X_ρ ∈ c^{D_ρ} satisfying the coprimality condition (26). Proposition 3.4 (p. 24):
H_L(P) = Norm(c^D)^{-1} Π_{ν∈M_K^∞} sup_{σ∈△_max} |X(P_0)^{D(σ)}|_ν (Salberger's formula, only archimedean factors,
for L globally generated).

## 3. Theorems, with exact statements and pages

**Conjecture 1.1** (McKinnon [McK07, Conj. 2.7]; p. 3), quoted: "Let X be a variety over a number field K, L be an
ample line bundle and ν be a place of K. Suppose that Q ∈ X(K) and that there exists a rational curve defined over K
passing through Q on X. Then there exists a rational curve C on X passing through Q achieving the best approximation
constant at Q with respect to L and ν, i.e., α_{L,ν}(Q, X) = α_{L,ν}(Q, C)."

**McKinnon-Roth, as quoted by Huang** (p. 3, eq. (3)): "One version [MR15, Theorem 6.3] of their main results is that
for any rational point Q ∈ X(K̄) and any ample line bundle L, α_{L,ν}(Q, X) ≥ (1/2) ε_L(Q), where ε_L(Q) is the Seshadri
constant of L at Q. Moreover, the inequality above is an equality if and only if both α and ε are computed on some
rational curve on X passing through Q."

**Principle** (p. 5), quoted: "Assume that the variety is rationally connected and the point to be approximated is
general. Then the best (resp. generic) approximations should be achieved on subvarieties swept out by free (resp. very
free) rational curves of small degree." (A heuristic, not a theorem.)

**Theorem 1.2** (p. 5; detailed version Thm 6.1), quoted: "Let X be a split smooth projective toric variety defined
over a number field K and equipped with a line bundle L. Let Q ∈ T(K). Suppose that X verifies Hypothesis (\*).
(1) If L is nef, then the best approximations for Q can be achieved on any free rational curve through Q of minimal
L-degree. (2) If L is ample, then the best approximations for Q are properly achieved on the subvariety swept out by
the free rational curves through Q of minimal L-degree." Context (p. 5): "smooth projective toric varieties X of
dimension ≥ 2"; "We assume that the point to be approximated is general, so that it lies on the open orbit".

**Theorem 6.1** (p. 34). Standing assumptions of §6 (p. 34): D = Σ 𝔞_i D_{ρ_i}, L = O_X(D) on X = X(△), "which we
assume to be smooth projective of dimension at least two and split over K. By torus action, we can assume that the
point to be approximated is the unit element (1, ..., 1) by Remark 4.2 (2)." Quoted: "Suppose that X verifies
Hypothesis (\*). Let Q_0 = (1, ..., 1) ∈ T(K). (1) Suppose that L is nef. Then for every place ν ∈ M_K, we have
α_{L,ν}(Q_0, X) = β. (2) Suppose that L is ample and X ≠ P^n. Then the constant α_{L,ν}(Q_0, X) is properly achieved on
a proper closed subvariety Y which is a finite union of Y_i ≅ P^{N_i}, each one being the fibre P^{N_i} × {1} of an open
toric subvariety of X isomorphic to P^{N_i} × G_m^{n−N_i}. Furthermore, if there exist two different such Y_i, Y_j, then
Y_i ∩ Y_j = Q_0." Remark 6.2 (p. 35): each Y_i is swept out by minimal rational curves of a centred primitive collection
of L-degree β and cardinality N_i + 1.
*Scope:* any number field K, any place ν (archimedean or not), Q_0 any point of T(K) (all equivalent), K-rational
approximants, X smooth projective **split over K**, (\*) required. What (2) gives about α_ess: in the proof (p. 43)
α(Q_0, X \ Y) ≥ β + δ for some unspecified δ > 0 (61), hence α_ess(Q_0) ≥ α(Q_0, X \ Y) > β = α(Q_0, X)
(p. 44: "α(Q_0, X \ Y) > β = α_ess(Q_0, Y) = α(Q_0, X)").

**Theorem 1.3** (p. 6; detailed version Thm 7.5), quoted: "Under the assumption of Theorem 1.2, suppose that X has
Picard number ≤ 2 (in particular they verify (\*)), and that L is nef. Then the generic best approximations for Q can be
achieved on every very free rational curve through Q of minimal L-degree."

**Theorem 7.5** (pp. 54-55), quoted: "Let X be a split projective smooth toric variety over K with Picard number 2.
Assume that L = O_X(D) is nef. (1) We have α_ess(Q_0) = deg_L C_3 if for all i ∈ {1, ..., t}, b_i ≠ 0;
deg_L C_1 + deg_L C_3 if there exists i ∈ {1, ..., t}, b_i = 0. In all cases, the generic best approximations can be
achieved (see Definition 4.8 (4)) on general lines (with respect to the cone σ_{t,n+1} (77)) passing through Q_0, which
are very free of minimal L-degree. Consequently, the essential constant equals the minimal L-degree of very free
rational curves [footnote 18: Edited after publication on 06.12.2021]. (2) Assume moreover that L is big. If either
a_t ≥ 2 or b_{t−1} = 0, then minimal free curves through Q_0 form a locally accumulating subvariety." (X =
P(O ⊕ O(a_1) ⊕ ... ⊕ O(a_t)) over P^s, Kleinschmidt classification, §7.1, pp. 48-49; C_1, C_2, C_3 relations (79),
(80), (82), pp. 50-51; Prop. 7.1 AE_1(X) = N C_1 + N C_2.) **This is the only place in the paper where α_ess is
computed exactly.**

**Section 4.3 example** (pp. 31-33): S_7 (toric dP7) over Q, L = ω^{-1}, ν = ∞, Q = [1:1:1]: α(Q, S_7) = 2 (attained on
l_1 ∪ l_2), α_ess(Q) = 3 = α(Q, l) for a general line l through Q.

## 4. The lower-bound mechanism (how dense sequences are kept away)

Everything below is read off the proofs, pp. 39-48 (Thm 6.1) and pp. 54-60 (Thm 7.5).

**Parametrisation** (47), p. 40. On the chart U_{σ_0} ≅ A^n given by the cone σ_0 of (\*\*):
π(X_1, ..., X_{n+r}) = (X_1/**Y**_1, ..., X_n/**Y**_n), **Y**_i := Π_{j=1}^r X_{n+j}^{𝔟_{i,j}}. By (\*\*) all 𝔟_{i,j} ≥ 0 are
integers, so every Cox ideal c^{D_ρ} can be chosen integral and every P ∈ T(K) has a lift with all X_ρ ∈ O_K \ {0}
(p. 40). Distance (48): d_ν(P, Q_0) = min(1, max_i |X_i(P_0)/**Y**_i(P_0) − 1|_ν).
Huang's own sketch (p. 40): "for P to approximate Q_0 ... the ν-adic values of the denominators (**Y**_i) and the
numerators (X_i) both tend to infinity, and they have almost equal sizes ... Salberger's height formula ... furnishes
us some flexibility of selecting maximal cones so as to control the growth at archimedean places of the numerators and
the denominators at the same time."

**Step 1: one non-zero integer** (51), p. 41: since P ≠ Q_0 some y_{i_0} ≠ 1, so **Y**_{i_0}(P_0) − X_{i_0}(P_0) ≠ 0,
an element of O_K. Product formula (52) (archimedean ν) / Lemma 5.1 (58) (finite ν):
|**Y**_{i_0} − X_{i_0}|_ν ≫ Π_{ν'∈S_1} |**Y**_{i_0}|^{-1}_{ν'} × Π_{ν'∈S_2} |X_{i_0}|^{-1}_{ν'}.

**Step 2: two monomial sections of L** (53)-(56), p. 42: X^{D(σ_{i_0})}/X^{D(σ_0)} = (X_{i_0}/**Y**_{i_0})^{deg_L(P_{n+j_{i_0}})}
(55), with deg_L(P_{n+j}) ≥ 𝔟_{i_0,j} deg_L(P_{n+j_{i_0}}) (56, from Prop. 6.7). Combining with the height formula
(Prop. 3.4): (57) Norm(c^D) H_L(P) d_ν(P, Q_0)^{deg_L(P_{n+j_{i_0}})} ≫ 1, and (59) for finite ν, "uniform for
P ∈ T(K) \ {Q_0}". By Prop. 4.5, α(Q_0, X) ≥ min_{i_0} deg_L(P_{n+j_{i_0}}) = β (p. 43). Upper bound α ≤ β by a smooth
rational curve in the class of a centred primitive collection (Prop. 4.4, Thm 2.2).

**Step 3: beating β off Y** (pp. 43-48, Thm 6.1(2)). Goal (61): d_ν(P, Q_0)^{β+δ} H_L(P) ≥ C_{K,δ,△} for
P ∈ (T \ Y)(K). Case (I) (p. 44): some coordinate y_{i_0} ≠ 1 with ρ_{i_0} in no CPC of degree β ⇒ the same one-integer
argument with deg_L(P_{n+j_{i_0}}) ≥ β + δ_1 (Prop. 6.7 strictness). Case (II) (pp. 44-48): two coordinates
y_1, y_2 ≠ 1 belonging to two different CPCs; uses the **product of two non-zero integers** X_i − **Y**_i (i = 1, 2) via
the elementary inequality (64) "max(Z_1, Z_2) ≥ Z_1^{λ_1} Z_2^{λ_2}" with λ_1 = λ_2 = 1/2 in (65), and five monomial
sections X^{D(σ_m)}, 0 ≤ m ≤ 4, in the height (66)-(68), split into Subcases (II.1)-(II.3) according to which of
|X_i|_{ν'}, |**Y**_i|_{ν'} is larger at each archimedean place. The gain δ_2 > 0 is unspecified ("for any δ > 0
sufficiently small such that (69)").

**Step 4: α_ess for Picard number 2** (pp. 56-59). Exceptional set (89): Z = ∪_{i=1}^n closure{y_i = 1} (the n
"coordinate hypersurfaces" through Q_0 of the chart U_{σ_{t,n+1}}). For P ∉ Z: one non-zero integer X_{t+j_0} − X_{n+1}
(92) plus the two monomial sections X^{D(σ_{t,n+1})}, X^{D(σ_{t,t+j_0})} = X_t^{deg C_1} X_{n+1 or t+j_0}^{deg C_3} (91)
give "Norm(c^D) H_L(P) d_ν(P, Q_0)^{deg_L C_3} ≫ Π |X_t^{deg_L C_1}|_{ν'} ≫ 1" (p. 58), hence (93)
α(Q_0, X \ Z) ≥ deg_L C_3. If some b_i = 0: two non-zero integers (y_{t+j_0} − 1 and y_t − 1) and four monomial
sections (94) give (95) α(Q_0, X \ Z) ≥ deg_L C_1 + deg_L C_3 (p. 59). Upper bound (87) from general lines
t ↦ (m_1 t + 1, ..., m_n t + 1) in U_σ (p. 54), which are very free (Lemma 7.2).

**What kind of argument this is** ([inference], but a direct reading of (51)-(59), (65)-(69), (92)-(95)).
* Each lower bound is a **Liouville / product-formula bound for a non-zero algebraic integer of the form
  X_i − **Y**_i**, the value at the Cox lift of the global section χ-difference "X_i − Π X_{n+j}^{𝔟_{i,j}}" of
  O(D_{ρ_i}), equivalently the toric character y_i − 1 = χ^{ρ_i^*} − 1, which vanishes to order 1 at Q_0. The archimedean
  size of this integer is compared with H_L through Salberger's formula, i.e. through a choice of **monomial sections
  of L** (maximal cones). In Case (II) and in (95) the integer is a **product of two such differences** (exponent
  weights from (64)). There is no counting, no equidistribution, no Roth/Schmidt/subspace input, no auxiliary
  polynomial of high degree, and no use of Zariski density beyond deleting the zero loci {y_i = 1} of these sections.
* The exceptional sets are fixed: Y (Thm 6.1(2)) and Z = ∪{y_i = 1} (Thm 7.5), depending only on (X, L, Q_0, △).
* Strength: in Thm 6.1 the exponent obtained is exactly β (and β + small δ off Y); in Thm 7.5 it is exactly the L-degree
  of the minimal very free lines, because there the upper bound from lines and the section bound **coincide**. The paper
  never proves a bound of the form α_ess ≥ dim X, and never proves α_ess > (best section bound) anywhere; all its exact
  values of α_ess are of the form "section bound = curve bound".
* **[inference]** For L = −K and Picard number 2, α_ess(Q_0) = minimal −K-degree of very free curves ≥ n + 1 (a very
  free f has f^*T_X = ⊕ O(a_i) with all a_i ≥ 1 and one a_i ≥ 2), so in Huang's setting α_ess(−K) > dim X: the points
  are *repelled* beyond the Vojta value, consistent with the record's split calibrations (P^1 × P^1: 4; dP6: 3). This
  is a feature of those varieties, not of a mechanism that targets dim X.

## 5. Relevance to Question B

Target (expert note §3): α^ess_e(−K) > 15/4 on X = Lbar_5^tw at the identity e, over Q-points T_5(Q), at least on
balanced sequences; McKinnon-Satriano canonical boundedness of X' would give ≥ 4.

**(a) Hypotheses of the theorems vs X. None of Huang's theorems applies.**
1. *Splitness is not the binding obstruction* **[inference]**. X_{Q(i)} is the split Losev-Manin 4-fold over K = Q(i),
   e ∈ T(Q(i)), and T_5(Q) ⊂ T(Q(i)). A lower bound d^γ H_{Q(i)} ≥ C on Q(i)-points off a closed Z restricts to Q-points
   (with H_{Q(i)} = H_Q^2, d_ν = dist^2, same γ), and α_ess over Q ≥ α_ess over Q(i) (a dense open over Q(i) contains its
   intersection with its conjugate, a dense open over Q). So any theorem for the split form over Q(i) would transfer.
2. *Hypothesis (\*) fails.* Checked with Huang's own criterion (\*\*) (Lemma 6.3, p. 35) on the permutohedral fan
   (rays 1_A, ∅ ≠ A ⊊ [k], in Z^k/Z1): for k = 3, 4, 5, **no** maximal cone satisfies (\*\*) (0 of 6, 24, 120; e.g. at
   k = 5 the flag cone {1_{s1}, 1_{s1s2}, 1_{s1s2s3}, 1_{s1..s4}} gives 1_{s2} = 1_{s1s2} − 1_{s1}, mixed signs; the
   representation is unique since the cone is a basis). So Eff(Lbar_k) is not simplicial for k ≥ 3 (geometric Picard
   rank 4, 11, 26; `huang_hyp_check.out`). Theorems 1.2/6.1 do not apply even to the split form, nor to split dP6 = Lbar_3.
3. *Picard number.* Theorem 7.5 (the only exact α_ess result) needs geometric Picard number 2; X has 26 (15 over Q).
4. *Smoothness / positivity.* Huang assumes X smooth projective (p. 10). X' (the Gorenstein Fano of the A_4 root
   polytope) is singular, so outside scope; on the smooth X, −K is nef and big but not ample, so Thm 6.1(2) (ample)
   would not apply even with (\*); Thm 6.1(1) and Thm 7.5 need only nef.
5. Consistency check **[inference]**: centred primitive collections of the permutohedral fan are exactly {1_A, 1_{A^c}}
   (sum 0), with −K-degree 2, so Huang's formula would give β = 2 = α_e(−K), matching the record's proved
   ε_e = α_e = 2 on the diagonal conics. (Not a theorem here, since (\*) fails; the record proved it independently.)

**(b) Could the mechanism transfer? Yes, but it is capped below the target** **[inference]**.
* Run by hand over Q(i), the mechanism is: a non-zero Gaussian integer value of a character difference χ^m − 1 (or a
  product of two), compared with monomial sections of L via Salberger's formula. For T_5 the characters are the ratios
  ζ_a/ζ_b and their products; the single-difference bound is the record's pair bound |t_ab| ≥ 1 (Liouville), and the
  product version is an instance of a product of section values. All of these lie in the record's class 𝓛
  (round 24 "Anisotropy": toric Liouville for one section, including isotypic sections and products), capped by the
  refereed width theorem (round 23, τ_e(−K) = 10/3 exactly; Prop 1.3(b) of round22 lit.md): **at most 10/3 for arbitrary
  profiles and exactly 15/4 on balanced profiles.** Huang's method therefore reproduces Cilleruelo-Córdoba / the
  balanced Liouville constant and cannot reach the window (15/4, 4).
* What made Huang's lower bounds sharp is that, in his (\*) / Picard-2 setting, the best section bound equals the
  degree of a covering family of very free rational curves through Q_0 (so upper and lower bounds meet). On X this
  coincidence is impossible: τ_e(−K) = 10/3 < 4 = dim X < 5 = n + 1 ≤ minimal very free −K-degree, so in the
  window (10/3, 5) neither sections nor smooth very free curves decide α_ess. Huang himself flags (pp. 6-7) that
  without (\*) "it is not always true that smooth rational curves of minimal degree contribute to the best or the
  generic approximations", and that singular rational curves (with α = degree / multiplicity at Q, [MR15, Thm 2.16])
  can do better; the record's balanced 15/4 curves on the twisted form (B·C ≥ 30, mult_e ≥ 8, round 24) and the k = 4
  dense family are of this singular type.
* **[inference, conditional on the record's gauss Thm 4.1]**: on the twisted Lbar_4 at e, the dense balanced family at
  δ < 1/14 gives α^ess_e(−K) ≤ 4/(1 + 1/7) = 7/2 < 4 ≤ minimal very free −K-degree. So Huang's Picard-2 identity
  "α_ess = minimal L-degree of very free curves" (Thm 7.5(1), last sentence) **fails outside his hypotheses** in this
  very family of varieties; it cannot be extrapolated to X.
* Huang's general lines (p. 54) on the permutohedral chart have −K-degree 2(k − 1) (= 8 at k = 5;
  `huang_lines_check.out`); their class is that of the record's curves C_v with five distinct entries (γ = 8). Over Q
  they are not available on the twisted form (Galois sends σ to −σ, so there is no Q-rational chart U_σ ≅ A^4; only
  the countably many Galois-stable subtorus curves C_v survive). Irrelevant to a lower bound in any case.
* Twisted-specific input: Huang's method uses no global Galois structure, so it cannot supply the "global
  twist-sensitive input" that round 24 says a T_5-only proof needs.

**Verdict for Question B: no help.** The paper neither implies α^ess_e(−K) > 15/4 nor provides a mechanism beyond
section-type Liouville; its hypotheses ((\*), smoothness, Picard number 2 for α_ess) fail for X and X', independently
of splitness.

## 6. Corrections to the record

1. **round22 lit.md §0 item 3 and README §1.1 item 5 / row "McKinnon conjecture ... all proved cases": "they give no
   lower bound on α_ess beyond α_e" is wrong for Huang.** Thm 6.1(2) (p. 34, proof p. 43-44) gives
   α_ess(Q_0) ≥ α(Q_0, X \ Y) ≥ β + δ > α for L ample, X ≠ P^n, under (\*) (δ > 0 unspecified), and Thm 7.5 (pp. 54-55)
   computes α_ess exactly for Picard number 2 (nef L). The verdict ("too weak / inapplicable") is unaffected because
   these lower bounds are section-type (item 3) and the hypotheses fail (item 4).
2. **Bibliographic.** The PDF (v6, 4 Aug 2025) prints no journal reference; "ANT 15 (2021) 461-512" is not confirmed
   by this source (footnotes 10 and 18 only show edits "after publication on 06.12.2021"). Def. 4.7 = essential
   constant credited to Pagelot [Pag08] is **confirmed** (p. 29), with the numbering of the arXiv v6.
3. **"Section-type" (round22 README row "Huang torsor lower bounds", lit.md §5 verdict 2, §8 table) — upgrade from
   assessment to verified by reading (pp. 40-48, 56-59).** Precisely: product-formula lower bounds for one non-zero
   algebraic integer X_i − **Y**_i (value of the character difference y_i − 1 on the Cox lift), or a product of two,
   compared with monomial sections of L via Salberger's height formula; exceptional sets = zero loci {y_i = 1} or the
   union Y of CPC fibres. This is inside the record's class 𝓛 and capped at 10/3 (arbitrary) / 15/4 (balanced) on X.
4. **"Split only" (lit.md §5 verdict 1) is true but not the decisive reason.** Splitness can be bought by base change
   to Q(i) for lower bounds (§5(a)1 above). The decisive hypotheses that fail are (\*) (Eff simplicial: fails for
   Lbar_3, Lbar_4, Lbar_5, checked with Huang's Lemma 6.3 criterion), smoothness (for X'), Picard number 2 (for
   Thm 7.5) and ampleness (for Thm 6.1(2) on X).
5. **lit.md §5 H3 summary, minor.** (i) Thm 6.1(2)'s "finite union of P^{N_i} ... pairwise meeting only at Q" needs
   **L ample and X ≠ P^n**; Thm 6.1(1) (α = β, every place ν) needs L nef only. (ii) "general rational point Q of the
   open torus": the theorems hold for **every** Q ∈ T(K) (Remark 4.2(2), torus action), K any number field, ν any place
   (archimedean or not). (iii) The "simplicial pseudo-effective cone" hypothesis is confirmed verbatim (Hypothesis
   (\*), p. 5).
6. **McKinnon-Roth citation detail (round22 README §1.1 "unconfirmed").** Huang p. 3, eq. (3), quotes
   "[MR15, Theorem 6.3]: for any rational point Q ∈ X(K̄) and any ample line bundle L, α_{L,ν}(Q, X) ≥ (1/2) ε_L(Q)",
   with equality iff both are computed on a rational curve through Q. This is secondary-source corroboration of the
   record's "α ≥ ε/2 (Thm 6.3)" form, not a reading of MR15 itself.
7. **Possible over-generalisation to avoid.** Huang's "α_ess = minimal L-degree of very free curves" (Thm 7.5) and the
   Principle (p. 5) should not be cited as predicting α^ess_e(−K) on X: on the twisted Lbar_4 the record's own
   Thm 4.1 gives α^ess ≤ 7/2 < 4 ≤ that degree ([inference] §5(b)).
