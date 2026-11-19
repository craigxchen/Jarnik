# Referee report, round 31, block "kodaira": Astra Part VI, section VI.3

Source: `.../scratchpad/round31/astra_lemma_VI/answer.md`, lines 1414-1455 (VI.3 and the subsection "Which genuine curves does this exclude?"). I also read II.3 (lines 292-316) and III.9 (lines 897-941, the last two paragraphs), and the record (`round27/split_notes.md` N2, `round31/README_III.md` §1.8, expert note §3).
Scripts: `.../scratchpad/round31/ref6L_kodaira/` (`h0M.py`, `h0M3.py`, `checkF.py`, `ordJ.py`).

## Verdict summary

| claim | verdict |
|---|---|
| (a) the cover is w^3 = F, F = prod_{0<=a<b<=5}(y_a-y_b), degree 15 on P^4 | CORRECT |
| (b) resolution B, discrepancies, -K_B, branch divisor, cyclic-cover formulas, canonical singularities | CORRECT (B is Kapranov's model of Mbar_{0,7}) |
| (c) h^0(B,jM) = 1 if 3 \| j (j >= 0), else 0; h^0(mK_Y) = 1 for m >= 1; kappa = 0 | CORRECT (proof checked line by line; j = 1..5 reproduced) |
| (d) not uniruled; no dominant family of rational curves with a constant cube class on each curve, even with varying twist | CORRECT |
| (e) for each D_0, the fixed rational carriers of degree <= D_0 lie in a proper closed Z_{D_0}; genus >= 1 carriers do not exist | CORRECT WITH CORRECTIONS (scope only; see §5) |

## 1. (a) The cover

II.3: J(u) = (2i)^{-10} prod_{i<j}(u_i-u_j)/(prod u_i)^2, and J = prod t_ij / H^3 at a genuine tuple (II.4). III.9: if J/lambda is a cube in Q(i), then the point lifts to z^3 = (c/lambda) prod q_i prod_{i<j}(q_i-q_j). The step is (prod q_i)^{-2} = prod q_i * (prod q_i)^{-3}. This is exact: J = c F/(prod q_i)^3 with F = prod q_i * prod_{i<j}(q_i-q_j). If we adjoin y_0 = 0, F becomes the product of all 15 differences of six labels. It is homogeneous of degree 5 + 10 = 15 on P(C^6/C1) = P^4. Over C all twists are isomorphic. Correct.

## 2. (b) The model B and the cover Y

* **Identification.** The coincidence flats of 5, 4 and 3 of the six finite labels are the following. The 6 points e_a lie in general position in P^4, because the e_a are 6 vectors of C^6/C1 summing to 0. There are 15 lines through two of them and 20 planes through three. The 15 diagonals y_a=y_b are the hyperplanes through four of them. So B is Kapranov's iterated blow-up, B = Mbar_{0,7}, with the seventh point at infinity. This is the minimal De Concini-Procesi wonderful model, so the boundary is SNC. Astra's local inductive argument is the standard one and is fine.
* **Discrepancies.** A size-s flat has codimension s-1, so its discrepancy is s-2. The order of F along it is C(s,2), the number of diagonals containing it. Hence K_B = -5h + E_3 + 2E_4 + 3E_5 and M = -K_B = 5h - E_3 - 2E_4 - 3E_5. For iterated blow-ups in increasing dimension, the order along the final exceptional divisor E_G is the multiplicity at the generic point of G. Independent check: the record (round 27, N2) has K_{Mbar_{0,7}} = -(1/3) sum_{|T|=2} D_T. That agrees with the Keel-McKernan coefficients k(n-k)/(n-1) - 2, which give -1/3 for k = 2 and 0 for k = 3 at n = 7. The 21 divisors D_T are the 15 strict diagonals plus the 6 divisors E_5 (label i collides with infinity, i.e. the other five finite labels coincide).
* **Branch divisor.** pi^*F = sum H~ + 3E_3 + 6E_4 + 10E_5. The reduced branch divisor is therefore Br = sum H~ + E_5 = B_2 (all 21 two-point boundary divisors), and Br = pi^*F - 3(E_3 + 2E_4 + 3E_5) lies in |3M|. Correct. So Y is the S_7-symmetric cyclic triple cover of Mbar_{0,7} branched on the whole two-point boundary, and B_2 lies in |-3K|.
* **Cover formulas.** The local equation is w^3 = x_1...x_k times a unit. This hypersurface is normal (S2, and its singular locus has codimension 2), so the simple cyclic cover is already normal. It gives pi_*O_Y = O + M^{-1} + M^{-2} and K_Y = pi^*K_B + 2R with 3R = pi^*Br, so K_Y = pi^*(K_B + (2/3)Br) = pi^*(K_B + 2M) = pi^*M. Correct.
* **Singularities.** C[x,w]/(w^3 - x_1...x_k) is the invariant ring of G = {(zeta^{c_i}) : sum c_i = 0 mod 3} acting on C^k. The invariants are t_i^3 and t_1...t_k, and both rings are domains of dimension k. G lies in SL_k, so the quotient is Gorenstein canonical (Reid-Tai/Watanabe). Astra's toric check is equivalent. Pluricanonical sections pass to a resolution. Correct.

## 3. (c) h^0(B, jM) and kappa(Y) = 0

* **Reduction to forms.** A section is a degree-5j form with multiplicity at least j, 2j, 3j along the planes, lines and points. This holds because the order along the final E_G equals the multiplicity at the generic point of G. Correct.
* **Restriction argument.** I checked it line by line. On H = {y_a=y_b} = P^3, the four planes {a,b,c} become four distinct planes of H. So f|_H is divisible by prod_c (y_c - y_a)^j, and the residual r has degree at most j. The 4-flat of the complementary labels {c,d,e,f} meets H in one point, which lies off those four planes. There f|_H, hence r, has multiplicity at least 2j > deg r, so r = 0. The 15 diagonals are distinct primes, so F divides f.
* **Induction.** The induction uses only the plane and line conditions (j-3 and 2j-6 after dividing by F; the point condition 3j-10 is not of the right form, and Astra correctly drops it). Call this space V_j. Then V_j = F V_{j-3}, and V_j = C F^{j/3} if 3 | j, j >= 0, and 0 otherwise. F^{j/3} has point multiplicity 10j/3 >= 3j. Therefore h^0(mK_Y) = h^0(mM) + h^0((m-1)M) + h^0((m-2)M) = 1 for every m >= 1. Correct; the proof is complete and elementary.
* **Independent computation** (`h0M.py`, `h0M3.py`, exact integer coefficients, rank mod 1000003). Method: the multiplicity conditions are imposed exactly, by changing to coordinates centred at a base label of each cluster. A mod-p rank is at most the rational rank, so the kernel dimension is a rigorous upper bound for h^0 over Q (and over C, by flat base change).

| j | degree | monomials | conditions | kernel mod p | h^0 |
|---|---|---|---|---|---|
| 1 | 5 | 126 | 258 | 0 | 0 |
| 2 | 10 | 1001 | 2149 | 0 | 0 |
| 3 | 15 | 3876 | 8423 | 1 | 1 (F satisfies all conditions exactly over Z: `checkF.py`) |
| 4 | 20 | 10626 (1251 after exact column elimination) | 4655 | 0 | 0 |
| 5 | 25 | 23751 (2751) | 10420 | 0 | 0 |

The variant without the point conditions (Astra's V_j) gives the same kernels for j = 1, 2, 3. j = 6 exceeded the 3 GB memory limit and was not run.

* **Literature.** The statement amounts to kappa(Mbar_{0,7}, -K) = 0: -K is Q-effective, and |{-3mK}| = {m B_2}. A web search did not find it stated explicitly; the search pointed to alg-geom/9607009 (Keel-McKernan), which says -K is not effective for n >= 8. It may well be folklore. Astra's proof is self-contained and correct.

## 4. (d) Consequences

* kappa = 0 implies the plurigenera are nonzero, so Y is not uniruled. Correct.
* Rational curves C_t with J|C_t = b(t) w^3: adjoin b^{1/3} on a finite cover of the base. Each curve then splits into three components birational to C_t on Y. A dominant family on T_5 would give a dominant family on Y, because pi is finite. So the varying twist is legitimately covered. Correct.
* **What it does not give.** kappa = 0 is the borderline case. Varieties with kappa = 0 (K3, abelian, Calabi-Yau) can have Zariski-dense rational points, so Lang-type conjectures predict nothing here. Astra says this, and it is the most useful message of VI.3: the cubic-cover route to non-density via general type is dead.

## 5. (e) "Which genuine curves does this exclude?"

I checked each step.

* **Heights.** The compactification has no real boundary points. E_A and E_{A^c} are conjugate and disjoint, since A and A^c are not nested. So the archimedean local heights are bounded and log n_A = h_{E_A+E_{A^c}} + O(1). Functoriality on P^1 gives d_A h(t) + O(1). The ten-pair profile along an infinite sequence forces d_A = 0 on singleton cuts and d_A = D/10 = l on pair cuts.
* **Degree of g_ij.** 2 deg g_ij = sum_P |a_i(P) - a_j(P)| = sum over cuts separating i and j of d_A. Exactly 6 pair cuts separate i and j, so deg g_ij = 3l. Correct.
* **Limit point.** Any algebraic irrational p gives Roth exponent 2, and a rational p gives exponent 1, so r is in {1, 2}. Note that r = 2 is not restricted to quadratic p; Astra's wording "r = 2 otherwise" is right. Then mu r >= 3l. Conjugates over Q(i) have equal multiplicity, a real p has the same degree over Q(i) as over Q, and the total zero degree is 3l. So p is rational (mu = 3l) or real quadratic (mu = 3l/2, l even). In both cases 3 | m, and the zeros of every g_ij - 1 are exhausted. Correct.
* **div J.** At a boundary point, ord J = sum min(a_i,a_j) - 2 sum a_i = sum over layers of c_A (C(\|A\|,2) - 2\|A\|). A layer of size 2 or 3 contributes -3; a layer of size 1 or 4 would contribute -2, but none occur since d_A = 0 there (`ordJ.py`, 2000 random chains). Equal orders with equal leading terms would create an extra zero of g_ij - 1, which is excluded. Above p, ord J = 10m. So div J = 3D', and J = b w^3 on P^1. Correct.
* **Bounded-degree family.** Maps P^1 -> P^4 of degree <= 3D_0/10 with J∘f = b w^3 form a constructible finite-type family. A dominant component would contradict (d). Z_{D_0} is the union of the closures of the finitely many non-dominant component images. Correct.
* **Two-sided endpoint.** In the rational case, |t - p| ~ e^{-h}; in the quadratic case, continued-fraction convergents (bounded partial quotients) give |t - p| ~ e^{-2h}. Either way Δ ~ H^{-3/10}, arc length ~ R^{2/5}, and the cut-weight errors are O_f(1). Correct.
* **Genus >= 2.** Faltings. Correct.
* **Genus 1.** The argument (Mordell-Weil cosets B + nQ, ĥ(P) = n^2 ĥ(Q) + O(sqrt(ĥ(Q))), Roth on x) is the standard proof of the elliptic-curve approximation theorem (Silverman, AEC IX.2.2: log d_v(P,Q)/h(P) -> 0). In fact it excludes far more than (*): it excludes any power proximity Δ <= H^{-epsilon}. Correct.

**Corrections and scope (what is and is not excluded).**

1. The closed set depends on D_0. Z_{D_0} is uniform over all curves of degree <= D_0 but not over D_0, and the union over D_0 is a countable union of proper closed sets, which may be Zariski dense. This is exactly the template-dependent-threshold situation, and Astra states it honestly ("a dense countable union of isolated curves" remains allowed). The corrected statement should say: *for each D_0 there is a proper closed Z_{D_0} (depending only on D_0) containing every irreducible curve of ratio degree <= D_0 that contains infinitely many (*)-points.*
2. Excluded: fixed carrier curves of bounded degree, and every fixed curve of genus >= 1. Not excluded:
   * sporadic points;
   * infinitely many curves of growing degree, each carrying finitely many points or an infinite sequence (the degree may grow with the index);
   * any sequence that is not contained in finitely many curves.

   A Zariski-dense set of (*)-points is in no way excluded. (*) itself remains undecided, as Astra says.
3. On effectivity and the standing cautions:
   * Z_{D_0} is defined geometrically, without Subspace exceptional sets or S-dependence, so the "uniform exceptional set" caution does not arise.
   * Which curves are carriers is ineffective (Roth), but the containment statement does not need that.
   * There are no gcd or denominator costs and no common phase in this section; the computation is purely geometric on the ratio torus.
4. Comparison with README_III §1.8 (III.6). The record shows that strict ten-pair templates of fixed block degree have union closure of dimension <= 2. VI.3 has a weaker conclusion (dimension <= 3, i.e. merely proper). Its hypothesis is broader:
   * it applies to every rational carrier, not only strict, simple-root, polynomially presented templates;
   * it allows quadratic contact points;
   * the presentation is derived, not assumed: the rigidity shows that every carrier is in effect a ten-pair template with all collisions at p and its conjugate, plus J = b w^3.

   It is a modest complement to III.6, not a strengthening of it.

## 6. Novelty and significance

* **Novelty.** A grep of the record finds no Kodaira-dimension computation for this cover. The record has the divisor relation (N2), the genus-4 fibre (III.9, verified in ref3_linear) and the cubic cover with no volume gain (README_II). The identification Y = triple cover of Mbar_{0,7} branched on B_2 in |-3K|, and the result kappa = 0, are new to the record, though likely known or easy for specialists. The derivation that carriers are rigid (only pair layers, all collisions at a rational or real quadratic point, and 3 | div J) is new and clean.
* **Significance for the goal: low.**
  * The positive part excludes only fixed bounded-degree carriers, which are essentially the template class III.6 already handles better.
  * The main value is negative and clarifying. Fixing the cofactor cube class gives a kappa = 0 fourfold (p_g = 1, not of general type), so no Bombieri-Lang or Vojta-type argument on this cover can yield A_5(delta) or decide (*).
  * The route should be recorded as closed.
