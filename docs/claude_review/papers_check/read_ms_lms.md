# Full-text check: McKinnon-Satriano erratum and Lehmann-McKinnon-Satriano (arXiv 2403.02480) vs Question B

Both PDFs were read in full: the erratum (2 pages), and LMS (12 pages, pp. 1-12). Quotations are verbatim from the PDFs,
up to TeX rendering. Statements marked **[inference]** are mine and appear in neither paper.

## 1. Bibliographic data as printed in the PDFs

**[E] Erratum.** Title: "ERRATUM FOR APPROXIMATING RATIONAL POINTS ON TORIC VARIETIES". Authors: David McKinnon and
Matthew Satriano, both of the University of Waterloo, Department of Pure Mathematics. It is 2 pages: p.1 has the text,
p.2 the references. The PDF prints **no date, journal, volume or DOI**. The upload's file name
(`coba-erratum-03-29-23`) suggests 29 March 2023, but that is not printed. Funding footnote: NSERC Discovery Grants,
and for MS a Math Faculty Research Chair. References: [MS21] McKinnon, D. and Satriano, M., *Approximating rational
points on toric varieties*, Trans. Amer. Math. Soc. 374 (2021), no. 5, 3557-3577; [Vo87] Vojta, *Diophantine
Approximations and Value Distribution Theory*, LNM 1239, Springer 1987.
LMS cite it (p.12) as: "[MS21a] McKinnon, D.; Satriano, M., 'Erratum to Approximating rational points on toric
varieties', Trans. Amer. Math. Soc., to appear."

**[LMS]** "APPROXIMATING RATIONAL POINTS ON SURFACES", Brian Lehmann (Boston College), David McKinnon and Matthew Satriano
(Waterloo). Margin stamp: "arXiv:2403.02480v1 [math.AG] 4 Mar 2024". 12 pages. Sections: 1 Introduction (p.1),
2 Preliminaries (p.4), 3 General Results (p.6), 3.1 Running the MMP (p.8), 3.2 Finding rational curves (p.9),
4 Surfaces (p.10), References (p.11).

## 2. Definitions (verbatim)

**"Canonically bounded" is not defined in either PDF.** The erratum uses the term without a definition; the definition
is in [MS21], which I did not have. The erratum's proof ends (p.1): "This implies α_P(−K_X) ≥ n − ε for all ε > 0. We
conclude that α_P(K_X) ≥ n, as desired." The second α has a sign typo, K_X for −K_X. The inequality proved there holds
for Q ∈ X(k) outside a proper closed Z (which depends on ε), so what is shown is the essential (Zariski-dense) bound
α_ess(P, −K_X) ≥ dim X. **[inference]** "X is canonically bounded at P" therefore means at least
α_ess(P,−K_X) ≥ dim X, which agrees with the record's reading. The exact MS21 wording is still unverified.

**LMS Definition 1.1 (p.2).** "Let (X, A) be a projective variety over a number field k equipped with a Q-Cartier
divisor A. For any algebraic point P, we let
• the *approximation constant* α(P, A) denote the infimum of α({x_i}, P, A) over all sequences of k-rational points
  {x_i} approximating P,
• the *essential approximation constant* α_ess(P, A) denote the infimum of all α({x_i}, P, A) over all sequences of
  k-rational points {x_i} approximating P such that {x_i} is Zariski dense.
Lastly, if P lies on a subvariety i: Z ↪ X, then considering P as a point on Z, we define α(P, A|_Z) := α(P, i^*A)."

**Setting of §2 (p.4).** "Throughout the rest of the paper, k is a fixed number field, and v a fixed place of k. We
further fix an extension of v to k̄, which we abusively also call v. We use this place only to fix a distance function
dist(·,·) on the rational points of X … height functions defined below are only defined up to multiplication by a
positive bounded function. These heights are only ever used to compute the value of α, where this ambiguity is
unimportant, so we will suppress precise definitions of the height functions." The relative-vs-absolute
normalisation of heights is not specified; for k = Q the two coincide.

**Definition 2.1 (p.4).** "Let X be a projective variety, P ∈ X(k̄), D a Q-Cartier divisor on X. For any sequence
{x_i} ⊆ X(k) of distinct points with dist(P, x_i) → 0, we set
A({x_i}, D) = {γ ∈ R | dist(P, x_i)^γ H_D(x_i) is bounded from above}.
If {x_i} does not converge to x we set A({x_i}, D) = ∅." The text adds that a nonempty A is an interval unbounded to
the right.

**Definition 2.2 (p.4).** "For any sequence {x_i} we set α({x_i}, P, L) to be the infimum of A({x_i}, L) (in particular
if A({x_i}, L) = ∅ then α({x_i}, P, L) = ∞). We call α({x_i}, P, L) the approximation constant of {x_i} with respect to
L." (This is the record's "α_{P,{x_i}}". LMS also write it with subscripts in Lemma 2.6: α_{P,{x_i}}(π^*D').)

**Definition 2.3 (p.4).** "The approximation constant α(P, L) of P with respect to L is defined to be the infimum of all
approximation constants of sequences of points in X(k) converging v-adically to P. If no such sequence exists, we set
α(P, L) = ∞." **Remark 2.4:** α(P, aL) = aα(P, L).

**Definition 2.7 (p.5).** "The essential approximation constant of a variety V over a number field k with respect to a
point P and a Q-Cartier divisor L on V is α_ess(P, L) = inf_{{x_i} Zariski dense} α({x_i}, P, L) where {x_i} ranges over
all Zariski dense sequences that converge to P. In other words, it is the infimum of all the approximation constants
of Zariski dense sequences."

Summary of the hypotheses in these definitions:
* the target P may be algebraic (P ∈ X(k̄));
* the approximating points are k-rational and distinct;
* the sequences in α_ess must be Zariski dense (as sets);
* a single place v is used;
* "approximation" is v-adic convergence.

**Definition 1.3 (p.2), essentially bounded.** "Let X be a smooth algebraic variety defined over a number field k, and
let P ∈ X be any k-rational point. Let K_X be the canonical divisor class on X. Then P is *essentially bounded on X* if
for every ample divisor A and every rational ε > 0, we have α_ess(P, −K_X + εA) ≥ dim X."

**Definition 4.1 (p.10).** "A smooth algebraic surface X is said to be *split* over a field k if the Galois action on
Pic(X) is trivial." By p.3 this "is a more general definition than other ones in the literature; for example, it
includes Brauer-Severi varieties."

**Conjecture 2.9 (Vojta's Main Conjecture, p.5), as LMS state it.** "Let X be a smooth, projective variety defined over
a number field k, with canonical divisor K_X. Let S be a finite set of places of k. Let A be a big divisor on X, and let
D be a normal crossings divisor on X. Choose height functions h_{K_X} and h_A … and define a proximity function
m_S(D, P) = Σ_{v∈S} h_{D,v}(P) … Choose any ε > 0. Then there exists a nonempty Zariski open set U = U(ε) ⊆ X such
that for every k-rational point Q ∈ U(k), we have
(2.10) m_S(D, Q) + h_{K_X}(Q) ≤ ε h_A(Q) + B
where B is a constant that does not depend on Q."

**Conjecture 1.2 (= [McK07, Conj. 2.7], p.2).** "Let X be a smooth projective variety defined over a number field k,
and P ∈ X a k-rational point. Let A be an ample Q-Cartier divisor on X. If α(P, A) < ∞, then there is a curve of best
A-approximation to P. In other words, there is a curve C such that α(P, A|_C) = α(P, A)."

**Conjecture 1.7 (p.3).** "… Let A be an ample line bundle on X. If α(P, A) < ∞, then there is a sequence of best
A-approximation to P that is not Zariski dense. In other words, there is a proper subvariety V such that
α(P, A|_V) ≤ α_ess(P, A)." **Remark 1.8:** Conjectures 1.2 and 1.7 are equivalent for surfaces.

## 3. Theorems, with pages

### 3.1 Erratum

Opening sentence (p.1): "We correct an error in [MS21, Proposition 1.3]: Vojta's Conjecture does not imply canonical
boundedness in general. However, the implication does hold for toric varieties. A corrected statement and proof of
Proposition 1.3 follows."

**Proposition 0.1 (p.1).** "Let X be a smooth projective variety over a number field k, and let P ∈ X(k) be a k-rational
point. If the anticanonical divisor class −K_X is big, then Vojta's Main Conjecture implies that X is canonically
bounded at P."

* **Scope.** The hypotheses are:
  * X is smooth and projective over an arbitrary number field k;
  * P is k-rational;
  * −K_X is big.

  The proposition is **not** restricted to toric varieties. "Toric" appears only in the motivating sentence (a
  complete toric variety has big −K). It has **no splitness hypothesis**.
* **The Vojta instance used.** S = {v}, a single place, and D is "the union of any n normal crossings divisors that
  intersect properly and transversely at P". The proof shows:
  * (0.2) m_S(D,Q) ≥ −n log dist_v(P,Q) + C for Q outside a proper closed Z. This uses φ: Y → X, the blow-up at P, the
    fact that φ^*D − nE is effective, and [Vo87, Lemma 1.3.3(b),(d)].
  * Vojta then gives (0.3) dist_v(P,Q)^n H_{−K_X}(Q) ≥ C_1 H_A(Q)^{−ε}, where A is ample with −K_X − A effective, so
    that H_{−K_X−A} > C_2 > 0.
  * Hence dist^{n−ε'} H_{−K_X} ≥ C_4.
* **What it gives for Question B.** It gives α_ess(P, −K_X) ≥ n, conditional on Vojta. Nothing is proved
  unconditionally.

### 3.2 Lehmann-McKinnon-Satriano: every numbered result

| # | page | statement (verbatim core) | conditional on |
|---|---|---|---|
| Prop 1.4 | 2 (proof pp.5-6) | "Let X be a smooth, projective variety defined over a number field k, and let P ∈ X(k). Then Vojta's Conjecture for X implies that P is essentially bounded." ("a variant of [MS21a]") | **Vojta** |
| Thm 1.5 | 3 (proof p.11) | "Let X be a split, smooth, projective surface defined over a number field k. Let A be an ample Q-Cartier divisor on X and let P ∈ X(k). If P is essentially bounded on X, then Conjecture 1.2 is true for P on X with respect to A. In particular, for an arbitrary smooth projective surface X, there is a finite extension ℓ/k such that Conjecture 1.2 is true for all A and all essentially bounded points of X ×_k ℓ." | **essential boundedness** (a hypothesis) |
| Cor 1.6 | 3 | "Let X be a split, smooth, projective surface … If Vojta's Conjecture is true for X, then Conjecture 1.2 is true for all P and A on X." | **Vojta** |
| Prop 2.5 | 4 | ([MR15, Prop 2.14(b)]) for {x_i} k-rational converging to P and positive rationals a, b: α({x_i},P,aL_1+bL_2) ≥ aα({x_i},P,L_1) + bα({x_i},P,L_2), except when {α(…,L_1), α(…,L_2)} = {−∞,∞} | unconditional (formal) |
| Lemma 2.6 | 4 | ([MS21, Lemma 4.7]) "Let π: X → Y be a surjective birational morphism of projective Q-factorial varieties over a number field k. Let P ∈ X(k) be a point which is not in the exceptional locus Exc(π) and let D' be a Q-Cartier Q-divisor on Y. Suppose either that {x_i} is a Zariski dense sequence on X converging to P and let x'_i := π(x_i), or suppose {x'_i} is a Zariski dense sequence on X' converging to π(P) and let x_i := π^{-1}(x'_i) whenever x'_i ∉ π(Exc(π)). Then α_{P,{x_i}}(π^*D') = α_{π(P),{x'_i}}(D')." | unconditional |
| Cor 2.8 | 5 | α(P,aL_1+bL_2) ≥ aα(P,L_1)+bα(P,L_2), and α_ess(P,aL_1+bL_2) ≥ aα_ess(P,L_1)+bα_ess(P,L_2), unless the pair is {−∞,∞} | unconditional |
| Prop 3.1 | 6 | Conjecture 1.2 for A implies it for yA (y ∈ Q_{>0}) | trivial |
| Prop 3.2 | 6 | "Let X be a smooth, projective variety …, P a k-rational point of X, A an ample Q-Cartier divisor on X such that K_X + A is in the closure of the effective cone. Then if P is essentially bounded, we have α_ess(P, A) ≥ dim X." (proof: [BCHM10, Thm 1.1], [HS00, Thm B.3.2(e)]) | essential boundedness |
| Thm 3.3 | 6-7 | ([MR15, Thm 2.16]) "Let C be any rational curve defined over k and φ: P^1 → C the normalization map. Then for any nef line bundle L on C, and any x ∈ C(k̄) we have the equality: α_{x,C}(L) = min_{q∈φ^{-1}(x)} d/(r_q m_q), where d = deg(L), m_q is the multiplicity of the branch of C through x corresponding to q, and r_q = 0 if κ(q) ⊄ k_v; 1 if κ(q) = k; 2 otherwise." | unconditional |
| Cor 3.4 | 7 | X projective over k, C a rational curve over k, A Q-Cartier with A·C ≥ 0; for P ∈ C(k) a v-adic limit of other k-points of C: α(P, A|_C) = M A·C with M rational, 0 ≤ M ≤ 1, independent of A. Remark 3.5: M = 1 if C is smooth at P | unconditional |
| Lemma 3.6 | 7 | (variant of [MS21, Lemma 4.3]) X projective Q-factorial, essentially bounded at P ∈ X(k), A ample with K_X + A nef, C irreducible through P with −K_X·C ≤ dim X and α(P,(K_X+A)|_C) ≤ α_ess(P,K_X+A) < ∞ ⇒ α(P,A|_C) ≤ α_ess(P,A) | essential boundedness |
| Def 3.7, Prop 3.8 | 8 | (K_X+A)-negative birational contractions; if φ: X ⇢ X' is an isomorphism near P and C' ∋ P' is a k-rational curve, its strict transform C satisfies α(P,A|_C) ≤ α(P',(φ_*A)|_{C'}) | unconditional |
| Thm 3.9 | 8-9 | X smooth, geometrically integral, projective over k; P ∈ X(k) essentially bounded; A ample with K_X+A pseudo-effective; φ a (K_X+A)-negative contraction that is an isomorphism near P; C' rational through φ(P) with φ_*A·C' ≤ dim X and φ(P) a v-adic limit of k-points ⇒ α(P,A|_C) ≤ α_ess(P,A) | essential boundedness |
| Lemma 3.10 | 9-10 | ψ: X → P^n birational morphism, P essentially bounded not in the contracted locus ⇒ there is a smooth irreducible C through P with −K_X·C ≤ dim X, hence α(P,A|_C) ≤ α_ess(P,A) | essential boundedness |
| Question 3.11 | 10 | lengths of extremal rays for terminal Q-factorial X over k̄ (−K·C ≤ n+1, resp. ≤ n if X ≇ P^n); answered in dimension ≤ 3 [CT09], for "many" threefold contractions [Be85],[Ka05], and "for split toric varieties in [MS21]" | open question |
| Question 3.12 | 10 | "Let X be a smooth del Pezzo surface over a number field k. Suppose that P ∈ X(k) is a k-rational point. Is there a curve C birational to P^1_k through P satisfying −K_X·C ≤ 6?" ("To the best of our knowledge the following question is open") | open |

LMS's account of the literature (p.2, verbatim in substance):
* [McK07] proved Conjecture 1.2 "for projective space, for split smooth rational surfaces of Picard rank at most four,
  and for a number of other split rational surfaces";
* "In [MR15], the second author and Roth proved the conjecture for split smooth cubic surfaces and gave the first cases
  of the conjecture for non-split surfaces";
* "Castañeda in [Ca19] proved many more instances … including all split smooth rational surfaces of Picard rank at most
  five" ([Ca19]: Castañeda, D., "Rational approximations on smooth rational surfaces", PhD thesis, University of
  Waterloo, 2019);
* "Huang proved Conjecture 1.2 for some split toric varieties in [Hu21]";
* "In [MS21], the second and third authors proved that Vojta's Conjecture implies Conjecture 1.2 for all Q-factorial
  split toric varieties."

### 3.3 Is any canonical/essential boundedness, or any α_ess > τ, proved unconditionally?

**No.** In both papers, every lower bound on α_ess for an actual Zariski-dense set of rational points either:
* comes from Vojta (erratum Prop 0.1, LMS Prop 1.4, Cor 1.6), or
* is the hypothesis "essentially bounded" (Thm 1.5, Prop 3.2, Lemma 3.6, Thm 3.9, Lemma 3.10, and the
  Kodaira-dimension-≥0 case of Thm 1.5 via [McK07, Thm 4.2]: "all that is used in the proof is essential boundedness
  at P", p.11).

The only unconditional statement about α_ess is degenerate (p.11): "If X is geometrically ruled over a curve B of
positive genus, then α_ess(P, A) = ∞ … because the rational points of X are not Zariski dense." Every other
unconditional result is one of the following:
* formal (Prop 2.5, Cor 2.8 concavity; Lemma 2.6 birational invariance);
* an **upper** bound on α coming from curves (Thm 3.3, Cor 3.4, Prop 3.8).

Neither paper says anything about non-split tori, Liouville or Seshadri lower bounds, or τ_x.

## 4. Mechanism

1. **Lower bounds on α_ess (conditional only).** Vojta's Main Conjecture is applied on X itself, with S = {v} and D a
   normal-crossings divisor of multiplicity dim X at P. LMS (p.5): "D could be the union of dim X transversely
   intersecting smooth divisors containing P in the class of mA for some large m". The proximity satisfies
   m_S(D,Q) ≥ −dim X · log dist(P,Q) − O(1), and Vojta then gives dist^{dim X} H_{−K_X+εA} ≥ N on U(ε) (LMS), or, after
   absorbing εA into −K_X using bigness, dist^{n−ε'} H_{−K_X} ≥ C_4 off Z (erratum).
   * **[inference]** This is a *different instance* of Vojta from the record's V_E(5) = Vojta(Bl_e X, E, {∞}).
   * Write Y = Bl_P X, φ^*D = D̃ + nE and K_Y = φ^*K_X + (n−1)E. Then
     h_{K_Y} + m(D̃+E) = h_{K_X} + m(D) + (n−1)N_S(E). So Vojta(Y, D̃+E) implies both Vojta(X,D) and Vojta(Y,E).
   * Each of the two yields h_{K_X} + n·λ_{P} ≤ ε h + O(1) after dropping nonnegative terms. Neither obviously implies
     the other. So the record's V_E(5) is a legitimate, and possibly weaker, sufficient instance.
2. **Upper bounds on α, and the existence of curves of best approximation.**
   * Curve computation: Thm 3.3 / Cor 3.4.
   * MMP for K_X + A [BCHM10], together with lengths of extremal rays: Mori's bend-and-break (−K·C ≤ dim X+1),
     Kawamata's bound (−K·C ≤ 2 dim X). These are used to find a k-rational curve C' with φ_*A·C' ≤ dim X.
   * Combined with Prop 3.2 (α_ess ≥ dim X under essential boundedness), this gives
     α(P,A|_C) ≤ dim X ≤ α_ess(P,A).

   This is the *opposite* direction to Question B: the essential lower bound is an input, and the curve supplies the
   upper bound. Splitness is used only to find rational curves over k. Its only role is in Thm 1.5/§4, for surfaces.

## 5. Relevance to Question B (X' at e: twisted, non-split, 4-dimensional)

1. **Conditional statement: confirmed, and with fewer hypotheses than the record suggests.** Our X = L̄_5^tw is smooth
   and projective over Q, e ∈ X(Q), and −K_X is nef and big. So erratum Prop 0.1 applies **as printed**:
   * Vojta's Main Conjecture for X (S = {∞}) ⇒ X is canonically bounded at e ⇒ α_ess(e,−K_X) ≥ 4 ⇒ A_5(δ) for every
     δ > 0 (record's Step 1).
   * LMS Prop 1.4 gives the same with −K_X + εA and needs no bigness.
   * **[inference]** Since −K_X = F + A' with A' ample and F effective (supported off e), Cor 2.8 gives
     (1+t)α_ess(e,−K_X) ≥ α_ess(e,−K_X+tA') + tα_ess(e,F) ≥ n + 0, so essentially bounded ⇒ α_ess(e,−K_X) ≥ n here.
   * Neither result needs splitness. The twisted, anisotropic structure is **irrelevant** to the conditional
     implication.
2. **X' versus X.** Prop 0.1 assumes smoothness, and X' is singular. **[inference]** X' is also not Q-factorial: a small
   birational morphism between Q-factorial varieties is an isomorphism. So LMS Lemma 2.6 does not apply literally to
   φ: X → X'. The transfer is nonetheless immediate:
   * φ is an isomorphism near e;
   * −K_X = φ^*(−K_{X'}), with −K_{X'} Cartier, so H_{−K_X} = H_{−K_{X'}}∘φ;
   * the distances are comparable near e;
   * Zariski density is preserved.

   Hence α_ess(e,−K_X) on X equals α_ess(e,−K_{X'}) on X'. **The correct literature statement of Question B is about the
   smooth weak-Fano X.** Whether MS21's definition of "canonically bounded" even applies to the singular X' is
   unverified (MS21 not read).
3. **Unconditional content applicable at e.**
   * Thm 3.3 / Cor 3.4 recompute the record's curve values. On C_v ≅ P^1, smooth at e (Remark 3.5, M = 1), they give
     α(e,−K|_{C_v}) = −K·C_v = 2(max v − min v), hence α_e(−K) ≤ 2 on diagonal conics.
   * **[inference, consistency check]** Thm 3.3 also explains round 24's numbers. A rational curve over Q of class
     B·C = 30 with two conjugate **real-quadratic** branches of multiplicity 4 at e has r_q = 2 (κ(q) ⊂ k_v = R,
     κ(q) ≠ Q), so α = 30/(2·4) = 15/4. Rational branches have r_q = 1 and give 30/4 = 15/2. This matches round 24's
     "alpha = 15/4 … (else the branches are rational and alpha = 15/2)", provided that round 24's "contact 4" is the
     branch multiplicity m_q.
   * These are **upper** bounds on α along curves. They feed the exceptional set Z, never a lower bound on α_ess.
4. **No unconditional case and no method.** Neither paper proves any case of canonical or essential boundedness
   unconditionally, for any variety, split or not. Neither proves α_ess > τ_x at any point. They offer no
   unconditional mechanism for lower bounds on α_ess.
   * **[inference]** The only reduction they offer is Question B ⇐ Vojta(X, D = D_1+…+D_4, {∞}) with D_i transverse
     through e. With n = dim X components through a *rational* point, the known unconditional Vojta-type inputs give
     at most the Liouville/τ value; this is the record's β ≤ τ argument, gauss §4.4 and lit §3. Those inputs are
     Schmidt subspace, Evertse-Ferretti, Corvaja-Zannier, Ru-Vojta and Levin.
   * One could choose D_i as closures of codimension-1 Q-subtori of T_5 through e. These exist because Galois acts by
     −1 on N, so every saturated sublattice is defined over Q. That makes λ_{D_i}(y) ≈ −log|χ_i(y) − 1| for characters
     χ_i. But T_5(Q) is not finitely generated, so S-unit/subspace methods (|u − 1| ≥ H(u)^{−ε} for u in a finitely
     generated group) do not apply. This is a remark, not a route.
   * The MMP/curve strategy of LMS §3 concerns where best approximations live (Conjecture 1.7). It is the dual
     problem, not Question B. Question 3.12 shows that even finding low-degree k-rational curves through rational
     points of del Pezzo *surfaces* is open.
5. **Verdict.** The papers give no help towards proving Question B. They confirm the record's conditional chain with
   cleaner hypotheses: smooth, −K big, any number field, no splitness. They also give exact definitions that match the
   record's Step 1 normalisation (dist^γ H_D bounded from above; α scales linearly in D; infimum over Zariski-dense
   sequences of distinct k-points).

## 6. Corrections to the record

1. **Scope of the erratum** (round22/README.md J1, expert_note §3). Proposition 0.1 is stated for **every smooth
   projective X over a number field with −K_X big, P ∈ X(k)**, not for toric varieties only. "Toric" occurs only in the
   introductory sentence "the implication does hold for toric varieties". The proposition requires **smoothness**, so
   the record's "Vojta implies [canonical boundedness of X'] when −K is big" should be routed through the smooth X
   (see §5 item 2).
2. **No splitness in the conditional implication.** The round22/README.md §2 table row ("McKinnon-Satriano,
   Monahan-Satriano, LMS … conditional on Vojta; split") and lit.md §2 item 4 should be corrected. Neither erratum
   Prop 0.1 nor LMS Prop 1.4 assumes split; splitness enters only LMS Thm 1.5 and Cor 1.6 (surfaces, Conjecture 1.2).
   In LMS "split" means trivial Galois action on Pic(X) (Def 4.1, p.10), which includes Brauer-Severi surfaces.
3. **The Vojta instance used is not V_E(5)** (lit.md §2 item 3: "comes from Vojta's conjecture on Bl_P X with D = E;
   that is exactly the record's V_E(5)").
   * The erratum and LMS use Vojta on **X** with D a normal-crossings divisor of multiplicity dim X at P (LMS p.5;
     erratum p.1). The erratum blows up P only inside the local-height estimate (φ^*D − nE effective).
   * Both instances follow from Vojta(Bl_P X, D̃ + E) and both give exponent dim X **[inference]**. The record's
     V_E(5) remains a valid sufficient hypothesis.
4. **Definition of α_ess.** lit.md §1.1 writes α_ess,x(L) = sup_U α_x(L; points of U) ("Pagelot's notion").
   * LMS Def 2.7 (p.5) defines it as the **infimum of α({x_i},P,L) over Zariski-dense sequences** converging to P,
     with α({x_i},P,L) = inf{γ : dist(P,x_i)^γ H_L(x_i) bounded above} (Defs 2.1-2.2, p.4).
   * The record's self-contained hypothesis (CB_e^norm1) and Step 1 match the LMS form exactly.
   * **[inference]** The sup_U form is ≤ the LMS form, so a bound "off a proper closed Z" implies the LMS bound.
5. **"Essentially bounded" ≠ "canonically bounded" literally.** LMS Def 1.3 (p.2) is
   α_ess(P, −K_X + εA) ≥ dim X for all ample A and rational ε > 0, and it needs no bigness. The erratum's
   "canonically bounded" is not defined in the erratum. Its proof concludes α_ess(P,−K_X) ≥ dim X, and this is why −K big
   is needed. **[inference]** For −K_X big with P off the augmented base locus, the two coincide via Cor 2.8.
6. **"San19" in lit.md §2** is **Castañeda, D., "Rational approximations on smooth rational surfaces", PhD thesis,
   University of Waterloo, 2019** (LMS [Ca19], p.11). LMS's description: "many more instances … including all split
   smooth rational surfaces of Picard rank at most five".
7. **What MS21 proved.**
   * LMS (p.2) state that "[MS21] … proved that Vojta's Conjecture implies Conjecture 1.2 for all Q-factorial split
     toric varieties". The record's snippet version (lit.md §2: "smooth split toric surfaces" plus "better
     approximated on a rational curve" for Q-factorial terminal split toric varieties) differs. LMS's one-line summary
     may be loose, so MS21 must be read to settle this.
   * LMS (p.10) also say that both parts of Question 3.11 (extremal-ray lengths) "are answered affirmatively for split
     toric varieties in [MS21]".
8. **Cubic surfaces attribution.** LMS (p.2) attribute "the conjecture for split smooth cubic surfaces and … the first
   cases of the conjecture for non-split surfaces" to **[MR15]** (Invent. Math. 200 (2015)), while the record (lit.md
   §2, §4) attributes cubic surfaces to MR16 (Eur. J. Math. 2016). This is unresolved here. Note that LMS say MR
   gave **non-split** cases, which the record's "all proved cases are split" (lit.md §2 item 4) overlooks. Those are
   surfaces and concern Conjecture 1.2 (best approximations), so they do not bear on α_ess.
9. **Erratum publication data.** The PDF prints no year or journal. "2023" in the record comes from the file name, and
   LMS (March 2024) list it as "Trans. Amer. Math. Soc., to appear".
10. **LMS's main result is conditional on essential boundedness, not directly on Vojta.** Thm 1.5 assumes P is
    essentially bounded. Only Cor 1.6 is the Vojta-conditional form. The record's "full conjecture for split surfaces
    (conditional on Vojta)" is accurate only for Cor 1.6.
