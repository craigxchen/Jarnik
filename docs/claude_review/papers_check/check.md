# Adversarial check of the four paper readers (Question B)

Method. Every page cited by a reader was reopened: text extraction to stdout (`pdftotext -f p -l p`, no copies of the
PDFs written), plus page images (Read tool) where an overline or exponent mattered (MR-I pp. 4, 33; Huang ANT p. 3;
LMS p. 6; Huang II p. 5 via glyph bounding boxes). Page numbers are PDF page indices. In all seven files the printed
page number, where there is one, equals the PDF index. Record files checked:
`docs/claude_review/round22/{lit.md,README.md,kill_lit.md}`, `expert_note/expert_note.md` §3, and
`round23/README.md`, `round24/README.md`. Labels: **OK** means verified word for word. **MISQUOTE**, **OVERREACH** and
**MIS-TARGETED** are flagged problems. **INF-OK** means a reader's inference that I re-derived. **UNCHECKED** means I
did not re-verify it.

Abbreviations: MR-I = 1306.2976v2, MR-L = 1306.2977v1, ERR = McKinnon-Satriano erratum, LMS = 2403.02480v1,
H-ANT = 1809.02001v6, H-I = 1703.01772v2, H-II = 1805.03920v2.

## A. Verdicts

1. All four readers give `no_help`. **Upheld.** No theorem in the seven files gives `alpha^ess_e(-K) > 15/4`, or
   any `A_5(delta)` with `delta < 1/30`, unconditionally. No reader claims "settles" or "narrows". The only route
   above `tau_e = 10/3` that appears anywhere is conditional on Vojta (ERR p. 1 Prop 0.1; LMS p. 2 Prop 1.4).

## B. McKinnon-Roth I (1306.2976v2), reader 1

2. Defs 2.7-2.9, p. 13: **OK**. Text: "We fix a single place v, archimedean or nonarchimedean". "A({x_i},L) = {γ ∈ R |
   d_v(x,x_i)^γ H_L(x_i) is bounded from above}". "Unless otherwise specified all height functions in this paper are
   multiplicative and relative to k." L is any line bundle (not only ample), and x ∈ X(k̄).
3. Def 3.1, p. 20 (Seshadri, x ∈ X(k̄), L nef): **OK**. γ_eff, f(γ) and β_x(L), p. 25: **OK**. β is defined for
   **ample** L. "β < γ_eff strictly" is INF-OK: f ≤ 1, f is continuous and Vol(L_{γ_eff}) = 0 (p. 25).
4. Cor 4.2, p. 27, "For any ample L and x ∈ X we have β_x(L) ⩾ n/(n+1)·(Vol(L)/mult_x X)^{1/n} ⩾ n/(n+1) ε_x(L)": **OK**.
   Cor 4.4, p. 28: **OK**. The bound β_e(-K) ∈ [2.31, 2.61] is **INF-OK**, recomputed as follows. Lower:
   (4/5)·70^{1/4} = 2.314. Upper: f = 1 − γ^4/70 on [0,2], which is equality by the Remark on p. 27. On [2, 10/3],
   f^{1/4} lies below its tangent at γ = 2. Together these give ≤ 2.604. *Caveats:* the upper bound uses γ_eff = 10/3,
   which is the record's round-23 theorem and not something in MR. And β needs L ample, so it must be computed on X'
   (where −K_{X'} is ample), not on X. The reader omits this second point. The value is unchanged, because the
   sections of m(π*L − γE) agree on Bl_e X and Bl_e X'.
5. Thm 5.1, pp. 28-29: **OK**. "X be an irreducible n-dimensional variety over Spec(k). For each v ∈ S choose a point
   x_v ∈ X(k̄)". L is an ample Q-bundle. (5.1): "There is a proper subvariety Z of X such that for all infinite
   sequences {x_i} of distinct points of X(k) ∖ Z(k), there is at least one v ∈ S so that α_v({x_i},L) ⩾ 1/R_v". It is
   essential-type. "α^ess ⩾ β" for S = {∞} is INF-OK: a Zariski-dense sequence has an infinite subsequence off Z, and by
   p. 13 Remark (b) a subsequence can only lower α.
6. Thms 6.1-6.3, p. 33, and their intro versions on p. 4: **OK**, with one precision point. *Precision:* as stated on
   p. 33, Thms 6.1 and 6.2 are dichotomies about α_x, the infimum ("either α_x(L) ⩾ β_x(L) or … α_x(L) is computed on a
   proper subvariety"). They are not literally essential statements. The essential forms are Thm 5.1, the "Theorem
   (6.2, alternate statement)" on p. 4 ("there is a proper subvariety Z ⊂ X so that for each δ > 0 there are only
   finitely many solutions y ∈ X(k) ∖ Z(k) to d_v(x,y) < H_L(y)^{−((n+1)/(n ε_x(L))+δ)}"), Thm 7.4 (p. 36), and the
   *proof* of 6.1. Reader 1's correction list calls "Thms 6.1, 6.2" essential-type without this distinction. That is a
   **minor overstatement**; the conclusion is unaffected.
7. Thm 6.3 equality clause, p. 33: **OK**. "equality if and only if both α and ǫ are computed on a k-rational curve C
   such that (1) C is unibranch at x, (2) κ(x) ≠ k, (3) κ(x) ⊂ k_v, and (4) ǫ_{x,C}(L|_C) = ǫ_{x_v,X}(L)". At a rational
   x, (2) fails, so the inequality is strict (INF-OK). Cor 6.4 (pp. 4, 33), exponent 2/ε_x(L) + δ: **OK**.
8. p. 24 Liouville: **OK**. "This proves the inequality α_x(L) ⩾ ǫ_x(L), a stronger inequality than α_x(L) ⩾ ½ǫ_x(L),
   one of the main goals of this paper." It is stated for x ∈ X(k). The "(See also [15, Theorem 3.3] …)" pointer is OK,
   and [15] is MR-L, "submitted" (p. 54 bibliography).
9. Thm 8.7, p. 47; Cor 8.9 and the table, p. 48: **OK**. "No gain from étale covers at e" is INF-OK in substance. Both X
   and X' are simply connected: the rays 1_A generate N, so N/⟨rays⟩ = 0. Constant-field extensions give nothing new,
   by Prop 2.12 (p. 15). That last step is my own addition to the reader's argument.
10. Thm 10.1, pp. 53-54: **OK**. "X an irreducible n-dimensional variety … such that −K_X is ample, D a finite subset
    …, If ǫ_x(−K_X) > (n+1)/n for every x ∈ D, then … Σ − log d_v(x,P) + h_{K_X}(P) < δh_A(P) + O(1)". The proof uses
    Thm 7.4(b) with R_v = 1, so at e the exponent is 1 (INF-OK). The reader's caution is correct. *It is pre-emptive:*
    grep finds no citation of MR Thm 10.1 anywhere in round22/23/24 or the expert note. So it is not a correction of
    an existing error.

## C. McKinnon-Roth Liouville (1306.2977v1), reader 1

11. Header, p. 1: "arXiv:1306.2977v1 [math.AG] 12 Jun 2013" and "Date: Revised on January 18, 2013." **OK.** The
    intro sentence: "On P1, except for the case that x ∈ P1 is a rational point of the base number field, Liouville's
    theorem is weaker than Roth's." **OK** (the reader dropped "∈ P1", which is harmless).
12. Lemma 3.2 (p. 6), Thm 3.3 (p. 7), Remark, Cor 3.4 (p. 7): **OK**. Thm 3.3 hypotheses: "x ∈ X(k̄) any point, d =
    [K∶k] … L an ample line bundle on X, and γ > 0 a rational number such that L_γ ∶= π*L − γE is in the effective cone
    … B′ the asymptotic base locus … B = π(B′)". (a): "if infinitely many points of {x_i} are outside B then α({x_i},L) ⩾
    γ/d". The Remark: "still holds if we replace B by the Zariski closure of B(k)". *Caveat on reader 1's relevance
    claim ("alpha^ess ⩾ 10/3 via the Vandermonde"):* L must be ample, so the theorem must be applied on X' with
    −K_{X'}, not on X. The value is unchanged (same reason as item 4). Apart from that caveat, INF-OK.
13. Non-split quadric example, p. 8: **OK**. "We assume that no lines on X are defined over k … π*L−2E is effective
    with base locus the proper transform of the two lines … x is the only k-point of B … αx(L) = 2 for all k-points of X."
14. Def 4.1 (p. 8) and Conj 4.2 (p. 9): **OK**. Conj 4.2 reads: "L any ample divisor … x any k-rational point … assume
    that there is a rational curve defined over k passing through x".
15. Thm 4.5, pp. 9-10: **OK**. "X be a smooth cubic surface in P3 defined over k, and isomorphic over k to the blowup of
    P2 at six k-rational points in general position. Let x ∈ X(k) be any k-rational point". ε ∈ {1, 3/2}; α ∈ {1, 3/2,
    2} exactly as the reader states. The proof uses "the asymptotic base locus of π*L − 2E is C̃_x … Theorem 3.3(b)"
    and Theorem 2.8. Thms 4.6 and 4.7 (p. 11) require "x … that does not lie on a (−1)-curve" and D nef. p. 12: "Conjecture
    4.2 holds for every x ∈ X not on a (−1)-curve". **OK.**

## D. Erratum (2 pp.), reader 2

16. Opening and Prop 0.1, p. 1: **OK, verbatim**. "Vojta's Conjecture does not imply canonical boundedness in general.
    However, the implication does hold for toric varieties." "Proposition 0.1. Let X be a smooth projective variety over
    a number field k, and let P ∈ X(k) be a k-rational point. If the anticanonical divisor class −K_X is big, then
    Vojta's Main Conjecture implies that X is canonically bounded at P." The proof takes "S = {v} … D be the union of any
    n normal crossings divisors that intersect properly and transversely at P", with the statement "for all Q ∈ X(k)
    outside a proper closed subset Z". The blow-up appears only to compare proximity functions ("φ*D − nE is
    effective … [Vo87, Lemma 1.3.3.(b)]"). The final line prints "αP(KX) ≥ n": the sign typo is confirmed.
    "Canonically bounded" is **not defined** in ERR (confirmed).
17. "Applies as printed to smooth X with −K big; no splitness": **OK** for the statement. The step to
    "α^ess_e(−K) ⩾ 4" is INF-OK *from the proof*. The proof's inequality holds off a proper closed Z, which gives an
    essential bound. The *statement*, by contrast, is in terms of the undefined notion "canonically bounded".
18. Bibliography: ERR p. 2 prints "[MS21] … Trans. Amer. Math. Soc. 374 (2021), no. 5, 3557–3577". So the erratum PDF
    does confirm MS21's journal data. It prints no date and no venue for the erratum itself (confirmed). "2023" comes
    from the file name only.

## E. Lehmann-McKinnon-Satriano (2403.02480v1), reader 2

19. Def 1.1 (p. 2), Defs 2.1-2.3 (p. 4), Def 2.7 (p. 5): **OK**. Def 2.1 has P ∈ X(k). Def 1.1 says "For any
    algebraic point P". Def 2.7: "inf over {x_i} Zariski dense of α({x_i},P,L) where {x_i} ranges over all Zariski dense
    sequences that converge to P". The phrase "Heights … only defined up to multiplication by a positive bounded
    function" (p. 4): **OK**.
20. Def 1.3 and Prop 1.4 (p. 2); proof pp. 5-6: **OK**. "Then Vojta's Conjecture for X implies that P is essentially
    bounded." D is "a normal crossings divisor defined over k with multiplicity dim X at P (For example, D could be the
    union of dim X transversely intersecting smooth divisors containing P in the class of mA …)". Conj 2.9 (p. 5) needs
    A big and D normal crossings, with U(ε) open. Prop 1.4 needs no bigness of −K: **OK**.
21. "Essentially bounded ⇒ α^ess(e,−K) ⩾ 4 for our X" is INF-OK, but the reader's justification is loose. A clean
    derivation: take A ample with −K − cA effective. Then (1+ε/c)(−K) = (−K+εA) + (ε/c)(−K−cA). Cor 2.8 (p. 5) and
    α^ess(F) ⩾ 0 for effective F (MR-I p. 14 Remark (c)) give α^ess(−K) ⩾ 4/(1+ε/c) for every ε > 0. The reader's
    hypothesis "e off the base locus" is not needed.
22. Thm 1.5, Cor 1.6 (p. 3) and Def 4.1 (p. 10): **OK**. Thm 1.5 has "If P is essentially bounded on X" as a
    hypothesis. Split means "the Galois action on Pic(X) is trivial", which "includes Brauer-Severi varieties" (p. 3,
    p. 11).
23. Thm 3.3 (p. 6, quoting [MR15, Thm 2.16]): **OK as a quote of LMS**: "Let C be any rational curve defined over k
    … for any nef line bundle L on C, and any x ∈ C(k̄)" (overline confirmed from the image). *Primary source differs:*
    MR-I Thm 2.16 (p. 19) reads "Let C be any **singular k-rational** curve … for any **ample** line bundle L". Both
    readers 2 and 4 quote secondary forms. The formula d/(r_q m_q) and the r_q cases are identical in all three
    sources. INF-OK: 15/(2·2) = 15/4 for the split Thm 2 curve (round23 README: "two conjugate real-quadratic branches
    of order 2 at e, so mult_e = 4"). Also 30/(2·4) = 15/4 for the round-24 class, and 30/4 = 15/2 when the branches
    are rational (matching round24 README §3).
24. Lemma 2.6 (p. 4, "projective Q-factorial varieties … P ∈ X(k) … not in the exceptional locus"): **OK**. "X' is not
    Q-factorial" is INF-OK. A birational morphism onto a Q-factorial variety has divisorial exceptional locus, and φ is
    small and not an isomorphism, since −K_X is not ample.
25. Prop 3.2 (p. 6), Lemma 3.6 (p. 7), Thm 3.9 (pp. 8-9), Lemma 3.10 (p. 9): **OK**. Essential boundedness is a
    hypothesis in all of them. p. 11: "α^ess(P,A) = ∞ for all ample divisors A simply because the rational points of X
    are not Zariski dense": **OK**.

## F. Huang, Rational approximations on toric varieties (1809.02001v6), reader 3

26. p. 1: "arXiv:1809.02001v6 [math.NT] 4 Aug 2025". The abstract says "general fixed point … split toric varieties".
    **OK.**
27. p. 3, eq. (3): **OK as a quote**. "for any rational point Q ∈ X(K̄) and any ample line bundle L, α ⩾ ½ε … equality
    if and only if both α and ε are computed on some rational curve on X passing through Q." *Flag:* Huang's equality
    clause omits MR's conditions (1)-(4) (item 7). It is a loose secondary paraphrase, so it should not be used for the
    equality case. The primary confirmation of "α ⩾ ε/2 (Thm 6.3)" is MR-I pp. 4 and 33 (item 6), and that supersedes
    this.
28. Hypothesis (∗), "The cone of pseudo-effective divisors Eff(X) is simplicial" (p. 5); Thm 1.2 (p. 5); Thm 1.3
    (p. 6): **OK**. Thm 1.3 says "Picard number ⩽ 2 (in particular they verify (∗))". p. 6-7: "Without Hypothesis (∗),
    the situation is noticeably more complicated …": **OK**.
29. Thm 6.1, p. 34: **OK**. "Suppose that X verifies Hypothesis (∗). Let Q0 = (1,…,1) ∈ T(K). (1) Suppose that L is nef.
    Then for every place ν ∈ M_K, we have α_{L,ν}(Q0,X) = β. (2) Suppose that L is ample and X ≠ Pn …" X is "smooth
    projective of dimension at least two and split over K" (p. 34). The proof (pp. 43-44) gives "α(Q0, X∖Y) = b … ⩾
    β + δ > β", so α^ess(Q0) ⩾ β + δ (INF-OK, since X∖Y is dense open; the intro on p. 8 states this too).
30. Lemma 6.3 / (∗∗), p. 35: **OK**. The claim "(∗) fails for L̄_3, L̄_4, L̄_5" is INF-OK. I re-derived it by hand,
    independently of the reader's script. In the basis of a maximal cone f_i = 1_{A_i} (flag A_1 ⊂ … ⊂ A_{k−1}), the ray
    1_{A_2∖A_1} = f_2 − f_1 has a positive coefficient. So no maximal cone satisfies (∗∗), for every k ⩾ 3.
31. Defs 4.1 and 4.7, Prop 4.5, Remarks 4.2 (pp. 26-29); Def 4.8 and Remarks 4.9 (p. 30): **OK**. Def 4.7: "α_ess(Q,Y)
    = sup_{V ⊂ Y} α(Q,V), where V ranges over all Zariski open dense subvarieties of Y such that inf_{P∈V(K)}
    d_ν(P,Q) = 0", credited to Pagelot [Pag08]. Footnote 12 says this is equivalent to the "constructable" version.
32. Thm 7.5, pp. 54-55: **OK**. "split projective smooth toric variety over K with Picard number 2 … L nef". Footnote
    18: "Edited after publication on 06.12.2021." **MISQUOTE (minor):** footnote 10 (p. 28) reads only "Edited on
    06.12.2021.", without "after publication". The reader attributed that phrase to both footnotes.
33. Mechanism, pp. 40-44: spot-checked, **OK**. (51) "Y_{i0}(P0) − X_{i0}(P0) ≠ 0", the product-formula bound (52),
    the monomial ratio (55), the conclusions (57) and (59), and the Case I / Case II split on p. 44. Equations (64)-(69)
    and (92)-(95): **UNCHECKED**.
34. "γ is unchanged under Q → Q(i)" is INF-OK, now with a primary source. MR-I Prop 2.12 (p. 15) gives
    "α_x({x_i},L)_K = (d/m_v) α_x({x_i},L)_k", and d/m_v = 2/2 = 1 for Q(i)/Q at ∞.

## G. Huang, Acta / BSMF papers (1703.01772v2, 1805.03920v2), reader 4

35. H-I §2.1, p. 7: **OK**. "On fixe un point rationnel Q ∈ X(k̄), ν une place réelle … k(Q) … contenu dans k_ν".
    Prop-Déf 2.2 cites [21, Prop 2.11]. Déf 2.3 uses sup over "parties constructibles denses". "on suppose que le fibré
    en droites L est gros". **OK.**
36. H-I Thm 1.2 (p. 3) and the η < 1/35 restriction (p. 5: "uniquement pour tous les η < 1/35"): **OK**. Y_4 is weak
    dP4: p. 36 "gros et engendré par ses sections globales mais il n'est pas ample". The 4A_1 contraction is on p. 49,
    and ω^{−1}_{Y4} = f*ω^{−1}_V = h*O(1,1) is on p. 50. **OK.**
37. H-I Prop 5.1 and Remark (p. 37), Lemmas 5.2-5.3 (p. 37), (5.10)-(5.11) and Thm 5.4 (p. 38): **OK**. Thm 5.4 reads
    "α_ess(Q) = α(Q,Y4) = 2. Par conséquent, il n'y a pas de sous-variétés localement accumulatrices". The remark
    printing zoom "1/2" (p. 37, "δ_{U,Q,B,1/2}(χ(ε)) = 0") is confirmed.
38. H-I Remarque 2.5 (p. 8): **OK** (C over k, L ample, Q ∈ C(k̄)). Thm 4.1 (p. 13): **OK**. Thm 5.5 (p. 41): **OK**.
    Thm 4.8 (p. 18) has 3/π² while Thm 1.3 (p. 4) has 6/π²: the discrepancy is confirmed. Cor 3.4 (p. 11): the free
    action and the orbit bijection are confirmed; the count "at most 3τ(m)" is **UNCHECKED**.
39. H-I Prop 7.1, pp. 50-51: **OK**. The lower bound is a product of the two Liouville bounds, (5.8) times the P^n bound
    (p. 50). The upper bound comes from the nodal curves D_{a,b,l} of degree 2n+6. *New, minor:* p. 51 calls the
    upper bound "la borne inférieure α_ess ⩽ n+3", a misprint in the source.
40. H-II Formules 1.2/1.3 (p. 3), Thm 1.4 (p. 4), Def 2.1 / Déf 2.2 / Thm 2.3 (pp. 10-11), Props 4.1-4.2 and Cor 4.3
    (pp. 15-17): **OK**. Def 2.1: "L un fibré en droites gros tel qu'il induise une application birationnelle … bien
    définie sur un ouvert U0 contenant Q0".
41. H-II Thm 1.5(1), p. 5: **MISQUOTE (minor).** The normalisation is B^{−(1−1/r)}, not B^{−(1−2/r)}. The glyph bounding
    boxes show numerator "1" over "r", and the proof on p. 17 prints "{B^{−(1−1/r)} δ_{Z,Q,B,r}}". Part (2), B^{1/5} at
    r = 5/2: **OK**. "Lemma 7.17 extra log B" and the Hooley/Pell description in Part II: **UNCHECKED**.
42. Reader 4's hand check "F = (x−y)²(t−s)²(xs−yt)yt ∈ H⁰(2(−K_{Y3})), ord_Q F = 5" is **INF-OK**. F has bidegree
    (4,4) = 2·(2,2). The blown-up points are read off from the gcds in Prop 4.2: x=t=0, y=s=0 and y=t=0. F vanishes
    to order exactly 2 at each of them, and ord_Q = 2+2+1 = 5. With MR-L Thm 3.3(a) this gives τ_Q = α_ess = 5/2. The
    same reasoning gives τ = α_ess = 2 on Y_4 (section xy(t−s)²) and n+3 on Y_4 × P^n.

## H. Proposed corrections to the record: verdicts

43. R1-1 (lit.md §3: "no exceptional set" holds only for Thm 6.3 / Cor 6.4; MR's Thm 5.1 is essential-type).
    **Upheld, with the precision of item 6.** Record text checked: lit.md §3 Statements ("No exceptional set appears in
    the inequality itself") and Hypotheses ("Exceptional sets: none in the Roth-type bound"). lit.md §3 verdict item 3
    already caps β ⩽ τ, so the conclusion stands.
44. R1-2 ("α ⩾ ε/2 (Thm 6.3)" confirmed). **Upheld.** MR-I p. 4 and p. 33. Close the "unconfirmed" item in round22
    README §1.1 (Killed/corrected, last bullet).
45. R1-3 (β < γ_eff strictly; MR's γ_eff = the record's τ). **Upheld** (MR-I p. 25). The identification γ_eff = τ is
    INF-OK: L_γ is effective iff some s ∈ H⁰(mL) has ord ⩾ mγ.
46. R1-4, R1-5. **Upheld** (MR-L p. 1; MR-I p. 24; MR-L Def 4.1 p. 8, Conj 4.2 p. 9).
47. R1-6 (the cubic-surface entry is split; Thm 4.5 covers all k-points for −K). **Upheld.** MR-L pp. 9-12. The record
    texts checked are lit.md §2 "cubic surfaces, points off lines (McKinnon-Roth 2016)" and the §4 table row "x off
    lines".
48. R1-7 (any single place). **Upheld** (MR-I p. 13). Immaterial.
49. R1-8 (asymptotic base locus; the Remark on B(k)). **Upheld** (MR-L p. 7).
50. R1-9 (MR Thm 10.1 is not V_E(5)). **Upheld as a caution.** The record never cites Thm 10.1 (item 10).
51. R1-10 (journal references not printed). **Upheld for the MR PDFs themselves, superseded by other uploaded files.**
    LMS p. 12 prints "[MR16] … Eur. J. Math. 2 (2016), no. 4, 929–959. DOI:10.1007/s40879-016-0113-5" and "[MR15] …
    Invent. Math. (200), 513–583 (2015)". H-ANT's bibliography (p. 62) prints the same data with "200(2)". The record's
    citations (lit.md §3) are therefore confirmed by secondary sources.
52. R2-1 (Prop 0.1 holds for every smooth projective X with −K big; route through smooth X). **Upheld** (ERR p. 1).
    Record texts checked: expert_note §3 ("by their 2023 erratum Vojta implies it when −K is big", stated for X') and
    round22 README J1. Add: whether MS21 defines "canonically bounded" for singular X' is unknown, because MS21 is not
    in the upload.
53. R2-2 (no splitness in the conditional implication). **Upheld for round22 README row "McKinnon-Satriano (TAMS
    2021), Monahan-Satriano, Lehmann-McKinnon-Satriano … conditional on Vojta; split"**, which lumps the implication
    with the split main theorems. **MIS-TARGETED for lit.md §2 item 4.** That item is about proved cases of McKinnon's
    conjecture, not about the boundedness implication.
54. R2-3 (the Vojta instance is on X with D through P, not on Bl_P X with D = E). **Upheld** (ERR p. 1 proof; LMS p. 5
    proof). The same correction applies to round22 README §1.2 item 7 ("This is gauss Thm 4.2 and … McKinnon-Satriano's
    corrected Prop 1.3"). The conclusion is the same, but the Vojta instance is not. Reader 2's further claim that both
    follow from Vojta on Bl_P X with D̃+E is **UNCHECKED** and not needed.
55. R2-4 (lit.md defines α_ess as sup_U, LMS as an infimum over dense sequences). **Rejected as a correction.** lit.md
    §1.1 attributes the sup_U form to Pagelot/Huang, and that matches H-ANT Def 4.7 (p. 29) verbatim. The two forms are
    **equal** (my inference; the reader only noted ⩽). There are countably many proper closed subsets over Q, and a
    diagonal sequence through U_n = X ∖ (Z_1 ∪ … ∪ Z_n) with d^{a+2ε}H ⩽ 1 is Zariski dense. Hence inf over dense
    sequences ⩽ sup_U + 2ε. No change is needed.
56. R2-5 ("canonically bounded" is not defined in ERR; LMS's notion is "essentially bounded"). **Upheld** (ERR pp. 1-2;
    LMS p. 2). "Agree for −K big" is INF-OK only as "essentially bounded ⇔ α^ess(P,−K) ⩾ dim" (item 21 plus Cor 2.8). Its
    relation to MS21's definition is unverifiable here.
57. R2-6 ("San19" = Castañeda thesis). **Plausible, not proved.** LMS p. 11 prints "[Ca19] Castañeda, D., 'Rational
    approximations on smooth rational surfaces', PhD thesis, University of Waterloo, 2019". LMS p. 2 credits it with
    "all split smooth rational surfaces of Picard rank at most five". The identification with lit.md's "San19" (which is
    "blow-ups of Hirzebruch surfaces at special configurations", cited via 2310.01640) is an inference.
58. R2-7 (LMS's summary of MS21 differs from the record's). **Upheld.** LMS p. 2: "[MS21] … proved that Vojta's
    Conjecture implies Conjecture 1.2 for all Q-factorial split toric varieties". LMS p. 10: "Both questions are
    answered affirmatively for split toric varieties in [MS21]".
59. **R2-8 (MAJOR, must not be adopted as phrased).** Reader 2 reports that LMS attribute "split smooth cubic
    surfaces and the first cases of the conjecture for non-split surfaces" to [MR15] (Invent.), "whereas the record
    attributes cubic surfaces to MR16". The quote is accurate (LMS p. 2), but **LMS mis-cite and the record is right**.
    MR-I (Invent., 1306.2976) contains no cubic-surface theorem; it has only the Seshadri example on p. 3. The cubic
    surface theorem is MR-L Thm 4.5 (pp. 9-10), and MR-L is LMS's own [MR16] (EJM). The "non-split" case is the
    non-split quadric example, MR-L p. 8. So the substantive half of R2-8 survives: lit.md §2 item 4 ("All proved cases
    are split …") overlooks MR16's non-split quadric. It does so only weakly, because the quadric is homogeneous, and
    round22 README §1.1 item 5 ("split or homogeneous") already covers it. The attribution half would introduce an error
    if adopted.
60. R2-9, R2-10. **Upheld** (ERR p. 2; LMS p. 12 "[MS21a] … Trans. Amer. Math. Soc., to appear"; LMS p. 3).
61. R3-1 ("no lower bound on α_ess beyond α_e" is wrong for Huang). **Upheld.** Record text checked: lit.md §0 item 3
    ("they give no lower bound on alpha_ess beyond alpha_e") and README §2 row "McKinnon conjecture, Coba … statements
    about alpha_e". H-ANT Thm 6.1(2) proof (p. 43) gives α^ess ⩾ β + δ, and Thm 7.5 (p. 54) gives exact values. H-II
    Thm 1.4 gives α_ess = 5/2 > α = 2. The verdict is unchanged.
62. R3-2 (bibliography). **Partly superseded.** "ANT 15 (2021) 461-512 not confirmed by this PDF" is true of H-ANT
    itself. But LMS p. 12 prints "[Hu21] Huang, Z., 'Rational approximations on toric varieties', Algebra and Number
    Theory, vol. 15 (2021), no. 2, 461–512". **The record's unconfirmed ANT page range is thereby confirmed (secondary
    source).** Footnote wording: item 32.
63. R3-3, R3-5. **Upheld** (items 29, 31, 33).
64. R3-4 (splitness is not the decisive obstruction). **Upheld** (items 30, 34). Record text checked: lit.md §5
    verdict 1.
65. R3-6 (Huang p. 3 is secondary confirmation of Thm 6.3). **Upheld, but superseded** by the primary source (item 27).
66. R3-7 (do not extrapolate "α_ess = minimal very-free degree"). **Upheld.** INF-OK: the record's k = 4 family has
    α = 14/4 = 7/2 (round22 README §1.2 Thm 4.4), while very free curves on a 3-fold have −K·C ⩾ 4. Huang himself warns
    on pp. 6-7.
67. R4-1 (Acta Arith. 189 (2019) 1-94 is not printed in either PDF; keep it snippet-level). **Rejected in part.** It is
    true that neither H-I nor H-II prints it. But H-ANT's bibliography (p. 61) prints "[Hua19] … Acta. Arith., 189:1–94,
    2019". It also prints **"[Hua20] … surface torique II. Bull. Soc. Math. France, 148(2):189–235, 2020"**. So the
    record's citation of Part I is confirmed by the author's own later bibliography. The record (lit.md §5 H2; §0 item 5;
    README row "Huang torsor lower bounds (ANT 2021, Acta Arith 2019)") should name Part II as BSMF 148 (2020) 189-235,
    not fold it under Acta Arith.
68. R4-2 (Part II is Y_3, not Y_4). **Upheld** (H-II p. 4 Thm 1.4; Figure 1).
69. R4-3, R4-4. **Upheld** (H-I Déf 2.9 p. 9; items 37, 39, 40, 42). Record text checked: README row "assessment,
    proofs unread".
70. R4-5 (three new α_ess = τ examples). **Upheld** (item 42).
71. R4-6 (Pell is a branch phenomenon on split surfaces, not anisotropy). **Upheld as inference.** It is supported by
    H-I p. 3 ("les pentes … sont réelles et irrationnelles … définis sur une extension quadratique réelle"), (5.11),
    and §6 p. 50.
72. R4-7 (refine "natural scale"). **Upheld** (H-I Prop 2.7 p. 9; Thm 5.5 p. 41; Thm 1.2(2)/η < 1/35).
73. R4-8 (kill_lit "stated for ample L" is inaccurate). **Upheld, but narrowed.** kill_lit.md §2 "Consequence" is
    correct for McKinnon's conjecture (ample; MR-L p. 9, H-ANT p. 3) and for MR's theorems (ample; MR-I pp. 25, 29,
    33). It is wrong for MR's *definition* (any line bundle, MR-I p. 13) and for Huang's constants (big: H-I p. 7, H-II
    p. 10; nef: H-ANT Thm 6.1(1), 7.5).
74. R4-9 (internal discrepancies in H-I). **Upheld** (items 37, 38).

## I. Summary of flags

* **Major (1):** item 59, R2-8. LMS's [MR15] label for the cubic-surface and non-split cases is a mis-citation. The
  record's MR16 attribution is correct and must not be "corrected".
* **Minor:**
  * Precision on Thms 6.1/6.2 (item 6).
  * L must be ample, so β and Thm 3.3 must be applied on X', not X (items 4, 12).
  * The R2-2 target (item 53).
  * R2-4 is not a real discrepancy (item 55).
  * Footnote 10 misquote (item 32).
  * H-II Thm 1.5(1) exponent misquote (item 41).
  * Huang eq. (3) equality clause (item 27).
  * Journal data the readers marked unconfirmed are confirmed by other uploaded PDFs (items 18, 51, 62, 67).
  * Part II's venue is BSMF, not Acta (item 67).
  * Primary vs secondary hypotheses of MR Thm 2.16 (item 23).
* **Unchecked:**
  * H-ANT eqs. (64)-(69) and (92)-(95);
  * H-I Cor 3.4's "3τ(m)" bound;
  * H-II Lemma 7.17 and its Hooley/Pell description;
  * reader 2's Bl_P X remark.
