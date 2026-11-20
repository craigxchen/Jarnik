# lit8L_cg: Cilleruelo-Granville, arXiv math/0608109, checked against the full text

## Answer first

* arXiv has **one version only**. The abs page lists just "[v1] Thu, 3 Aug 2006 20:31:28 UTC (23 KB)", with the comment
  "21 pages, preliminary version. Comments welcome". It gives no journal reference. Both https://arxiv.org/abs/math/0608109v2
  and https://arxiv.org/pdf/math/0608109v2 return 404. So there is nothing to compare between versions. All numbering
  below is v1 numbering.
* I could not open the published CRM Proc. Lecture Notes 43 text. ams.org and the author's homepage at dms.umontreal.ca
  both failed with "CONNECT tunnel failed, response 403". The published numbering therefore remains **NOT FOUND**.
* I downloaded the v1 PDF into `lit8L_cg/dl/cg_v1/v1.pdf`. It is byte-identical (same md5) to `ref6L_pell/dl/cg07.pdf`.
  I read all of `v1.txt` (pdftotext -layout, 1232 lines). I also rendered pp. 7, 10 and 15 as images to check the
  displayed formulas and the two flowcharts.
* **(a)** In v1, no arc conjecture (Conj. 13, 14, 15) is linked to Rudin's conjecture or to squares in arithmetic
  progressions, in the text or in either flowchart. Rudin's conjecture follows from Conj. 9 (Theorem 11) and from
  Conj. 12 (Theorem 12). Both are congruence statements. The arc conjectures are linked only to each other and, via
  Theorem 16, to the L^4 Conjecture 19. **The record's "no" is confirmed (READ).**
* **(b)** Under Bombieri-Lang, CG get sigma(k) << k^{4/5} (through Theorems 2, 3, 4) and sigma(k) << k^{5/7} (a direct
  count). Both rest on the Caporaso-Harris-Mazur uniformity [8]. CG add that both bounds are "improved unconditionally
  in [3] and [4]". (READ)
* **(c)** Theorem 10 is about **k** roots, with ell the largest odd integer <= k. It is not about "k+1 roots". The
  statement as printed fails for k = 2, but checks numerically for k = 3..6. CG state **no** uniform arc conjecture at
  length R^{1/2}. Their arc conjectures are at length R^{1-eps}. (READ)
* **(d)** Conjecture 10 is stated at exponent **1 - eps for every eps > 0**, which covers all exponents < 1. It is not
  stated at q^{1/d - eps} or at q^{1/d}. "TR(1/2+eps) is implied by Conjecture 10 (d = 2)" is **right**: it is the
  instance eps_CG = 1/2 - eps, up to an evident translation from CG's initial segment [0, b^{1-eps}) to arbitrary
  intervals. The record's gloss "at a fixed exponent above their critical 1/2" is accurate but can mislead. "Critical"
  is the record's word, not CG's. CG only say they "strongly believe" the version with exponent 1/d. Suggested wording
  is in (d) below. (READ)

## (a) Statements and numbering (v1; the only arXiv version). All READ.

**Conjecture 9** (§4, p. 8): "There exists a constant N such that there are no more than N solutions
0 < x_1 < x_2 < ... < x_N < x_1 + b^{1/2} to the equation x_i^2 ≡ a (mod b), for any given a and b."
* It is introduced with "A slightly improvement on the theorem above would have interesting consequences."
* It holds for all a, nonunits included. (Read literally, it bounds the count by N - 1.)
* Flowchart box 9 (p. 10): "∃ m such that x_i^2 ≡ r (mod q), 1 ≤ i ≤ m ⟹ max |x_i − x_j| > q^{1/2}".

**Theorem 11** (p. 8): "Conjecture 9 implies Conjecture 1." Conjecture 1 (§1) is Rudin's "σ(k) = O(k^{1/2})".

**Conjecture 10** (§4, p. 8): "Let integer d ≥ 1 be given, and suppose that for each prime power q we are given a set
of residues Ω(q) ⊂ (Z/qZ) which contains no more than d elements. Ω(b) is determined for all integers b using the
Chinese Remainder Theorem, as described at the start of this section. Then, for any ǫ > 0 there exists a constant
N(d, ǫ) such that for any integer b there are no more than N(d, ǫ) integers n, 0 ≤ n < b^{1−ǫ} with n ∈ Ω(b)."
* It is introduced with "Conjecture 9 can be strengthened and generalized as follows".
* Next comes "In theorem 9 we proved such a result with the exponent '1 − ǫ' replaced by '1/d − ǫ'. We strongly
  believe Conjecture 10 with '1−ǫ' replaced by '1/d', analogous to Conjecture 9."
* Poonen (1995 email) "asked Conjecture 10 with '1 − ǫ' replaced by '1/2' for d = 4". This would imply uniform
  boundedness for rational preperiodic points of quadratic polynomials.
* **Exponent:** the conjecture itself is at q^{1-eps} for every eps. q^{1/d - eps} is Theorem 9, which is proved, with
  eps = eps_d(k). q^{1/d} is the separate "strongly believe" version.
* Flowchart (p. 10): box 10 "|Ω(b) ∩ [0, b^{1−ε}]| < C(d, ε)", with arrows 10 ⟹ 9 ⟹ 1 and 10 ⟹ 11.

**Theorem 9** (p. 7), the proved version. For d ≥ 2 and |Ω(q)| ≤ d, "for any k ≥ d, there are no more than k integers
x ∈ Ω(n) in any interval of length n^{α_d(k)}, where α_d(k) = (1−ε_d(k))/d > 0 with 0 < ε_d(k) = (d−1)/k + O(d^2/k^2)."
I checked the fraction against the rendered page.

**Conjecture 11** (p. 9) is the monic-polynomial version: degree d, b^{1-eps}, with p^2 | b ⟹ p ∤ disc f.

**Conjecture 12** (§5, p. 9): "There exists δ > 0 and integer m > 0 such that if a_i^2 + b_i^2 = n with a_i, b_i > 0
and a_i^2 ≡ a_1^2 (mod q) for i = 1, ..., m then q = O(n^{1−δ})."
* Distinctness of the m pairs is not written but is clearly intended. Without it the statement is false.

**Theorem 12** (p. 9): "Conjecture 12 implies Conjecture 1". The proof works with n = r_j^2 + s_j^2 and the congruences
r_j, s_j ≡ ±w (mod p^e) for p^e ∥ b. It is a congruence argument, not an arc argument.

**Conjecture 13** (§5, p. 11): "For any α < 1/2, there exists a constant C_α such that for any N we have
#{(a, b), a^2 + b^2 = n, N ≤ |b| < N + n^α} ≤ C_α." The special case N = 0 is (5.1).
* CG: "It is simple to prove (5.1) for any α ≤ 1/4 (and Conjecture 13 for α ≤ 1/4 with N ≪ n^{1/2−α}), but we cannot
  prove (5.1) for any α > 1/4."

**Conjecture 14**: "The number of lattice points {(x, y) ∈ Z^2 : x^2 + y^2 = R^2} in an arc of length R^{1−ǫ} is
bounded uniformly in R."

**Conjecture 15**: the same "in an arc of length R^{1−ǫ} around the diagonal".

**What the arc conjectures are tied to:**
* "Conjecture 13 and the special case 5.1 are equivalent to the following conjectures respectively". Also: "Conjectures
  13 and 14 are simply a rephrasing of one another, and obviously imply (5.1) and Conjecture 15."
* A rotation argument ("In the other direction ...") derives Conj. 14 from (5.1) and from Conj. 15. In effect, all four
  are presented as equivalent.
* **Theorem 16** (§6, p. 14): "Conjecture 13 implies Conjecture 19." This is an implication, not an equivalence.
  Conjecture 19: "For any α < 1, for any trigonometric polynomial f with frequencies in the set {N^2, ..., (N + N^α)^2},
  we have ‖f‖_4 ≪_α ‖f‖_2." It is "trivial for α ≤ 1/2, yet is completely open for any α > 1/2".
* Second flowchart (p. 15): arrows 15 ⟹ 14 ⟹ 19 ⟹ 21.
  * The boxes there are labelled with 14 and 15 **swapped** relative to the text. Box "14." says "around the diagonal".
    Box "15." is the plain arc.
  * No arrow leaves this cluster toward Conjectures 1, 2, 9, 10 or 12.

**Arc ⟹ Rudin / squares in AP?** No, in v1, which is the only arXiv version. I searched the whole text and both
flowcharts.
* Section 5's lattice-point material on Rudin is Conjecture 12 / Theorem 12, which is congruence based.
* Theorem 13, from Cilleruelo-Córdoba [13], and items (1)-(4) after it are about arcs only, with no consequence for
  Rudin.
* The published CRM version could not be opened (NOT FOUND).

## (b) Bombieri-Lang. READ, §3, pp. 3-4.

* **Theorem 2**: "Assume the Bombieri-Lang conjecture. Then Σ_n r^2_{E+E}(n) ≪ |E|^{11/4}". The proof starts "One
  consequence of [8] is that there exists an integer B, such that if the Bombieri-Lang conjecture is true then for any
  polynomial f(x) ∈ Z[x] of degree five or six which does not have repeated roots, there are no more than B rational
  numbers m for which f(m) is a square." Here [8] is Caporaso-Harris-Mazur, "Uniformity of rational points", JAMS 1997.
* **Theorem 3**: "Conjecture 4 implies Conjecture 5 (with the same ε)." **Theorem 4**: "Conjecture 5 implies Conjecture 2
  (ε → ε/(4−2ε))".
* Then: "In particular, theorems 2, 3 and 4 show that the Bombieri-Lang conjecture implies σ(k) ≪ k^{4/5}, which is easy
  to obtain by directly applying the Bombieri-Lang conjecture to our arithmetic progression. To do better than this ..."
  * The argument counts rational points on y^2 = x ∏_{j=1}^5 (x + n_j), summing over r mod s and s > σ/10.
  * It ends "and we obtain σ(k) ≪ k^{5/7}. Anyway this upper bound was improved unconditionally in [3] and [4]."
  * [3] is Bombieri-Granville-Pintz, k^{2/3+o(1)}. [4] is Bombieri-Zannier, k^{3/5+o(1)} (§1).
* Flowchart: "Bombieri-Lang Conjecture ⟹ Conjecture 5 for ε = 3/4". Check: ε = 3/4 gives 3/10 in Theorem 4, so
  1/2 + 3/10 = 4/5, consistent. Also Bombieri-Lang ⟹ Conj. 6 (Solymosi).
* Attribution: k^{4/5} and k^{5/7} are CG's own remarks in §3 of v1, not numbered theorems. They are conditional on
  Bombieri-Lang through the uniformity of [8].
* The record's "Bombieri-Lang gives σ(k) ≪ k^{5/7} (Cilleruelo-Granville)" is accurate. Better: "Bombieri-Lang, via
  Caporaso-Harris-Mazur uniformity". A citation should give "§3", not a theorem number.

## (c) Theorem 10 and arcs of length R^{1/2}. READ.

**Theorem 10** (p. 7): "If x_1 < ··· < x_k are solutions to the equation x_i^2 ≡ a (mod b), then
x_k − x_1 > b^{1/2 − 1/(2ℓ)}, where ℓ is the largest odd integer ≤ k."
* The first proof applies Theorem 9 with d = 2 and k−1 in place of k, noting "α_2(k − 1) = 1/2 − 1/(2l)". Theorem 9
  needs k − 1 ≥ d, so k ≥ 3.
* The second proof uses a Vandermonde determinant divisible by b^{[(k−1)^2/4]}.
* The task's phrase "k+1 roots ... b^{1/2−1/(2l)}" is the same statement with k shifted: ℓ = largest odd ≤ k+1 when
  there are k+1 roots. ref8L_meaning R6 states it in that form correctly.
* **Slip (Computed, my script `lit8L_cg/scripts/t10.py`).** For k = 2, ℓ = 1, the claim is x_2 − x_1 > 1. This fails:
  for b = 3, a = 1 the roots are 1 and 2. There are 749 such cases for b < 1500.
* For k = 3, 4, 5, 6 there are no violations for b < 1500, counting runs of consecutive roots and allowing wrap-around.
  So the statement should read k ≥ 3.
* The Remark after Theorem 10 derives "a weak version of Conjecture 9": at most log(4b)/log 2 roots in an interval of
  length b^{1/2}. It credits Konyagin-Steger [20] with the degree-d analogue at 1/d.

**Arcs of length R^{1/2}. NOT FOUND in CG** (READ: absent).
* The only arc conjectures are 13-15, at length R^{1−ǫ} on x^2 + y^2 = R^2. Equivalently, Conj. 13 bounds b-windows of
  length n^α, α < 1/2, with n = R^2.
* At the R^{1/2} scale CG give only Theorem 13 (from [13]): "There no more than k lattice points ... in an arc of length
  R^{1/2 − 1/(4[k/2]+2)}". Then come the sharpness items (1)-(3) for k = 1, 2, 3.
* Item (4) asks: "Are there infinitely many circles x^2 + y^2 = R_n^2 with four lattice points on an arc of length
  ≪ R_n^{2/5}?"
* So the record's "U ... the uniform Cilleruelo-Granville bound M(C)" (critical arcs C R^{1/2}) does **not** come from
  this paper. If it is attributed to Cilleruelo-Granville, the source must be another paper, probably [15] "Close
  lattice points", listed here as "In preparation". I did not read it.

## (d) Does "TR(1/2+eps) is the case d=2 of CG Conjecture 10 at a fixed exponent above their critical 1/2" match? READ.

The literal Conjecture 10 is quoted in (a). Its exponent is 1 − ǫ for **every** ǫ > 0.

* **"TR(1/2+eps) is implied by Conjecture 10 (d = 2)": right.** Take d = 2 and ǫ_CG = 1/2 − eps.
  * Then N(2, 1/2 − eps) bounds the count in [0, b^{1/2+eps}). Any family of at most 2 residues per exact prime power is
    allowed, including residues that coincide modulo a proper divisor, since CG only require |Ω(q)| ≤ d.
  * CG count only in the initial segment 0 ≤ n < b^{1−ǫ}, while TR uses arbitrary real intervals.
  * Translation bridges this, since Ω(q) − t is again an admissible family and CRT commutes with shifts. That is an
    evident step, not written by CG.
  * CG's flowchart arrow 10 ⟹ 9 likewise uses the instance ǫ = 1/2.
* **"The case d = 2 of Conjecture 10": slightly loose.** TR(1/2+eps) for one fixed eps is a single instance of
  Conjecture 10 (d = 2), which asserts all exponents < 1. It is not the whole d = 2 case.
* **"Their critical 1/2".** CG never use the word "critical". The exponent 1/2 = 1/d appears only in "We strongly believe
  Conjecture 10 with '1−ǫ' replaced by '1/d', analogous to Conjecture 9." Literally, that 1/d version is itself the
  instance ǫ = 1 − 1/d of Conjecture 10.
* The chain is: Conj. 10 (d = 2, all ǫ) ⟹ TR(1/2 + eps) ⟹ TR(1/2), which is the "strongly believe" form. No converse is
  known. This matches the judge's correction in README_VIIIL §1.3 and ref8L_meaning R5.
* **Suggested wording:** "(a special case, d = 2 and exponent 1/2 + ε, of Cilleruelo–Granville's Conjecture 10, which
  they state at every exponent 1 − ε < 1; they 'strongly believe' the version at exponent 1/d = 1/2, which TR(1/2+ε)
  implies)".

## Other checks against the record

* README §1.2: "Theorem 11 is 'Conjecture 9 implies Conjecture 1', and Theorem 12 is 'Conjecture 12 implies Conjecture
  1'". **READ, correct.**
  * The sentence "The published numbering could not be checked" still stands, since the CRM text is unreachable.
  * However, arXiv has only v1, so "arXiv math/0608109v1 numbering" is the only arXiv numbering.
* README §1.7, "no arc conjecture is linked to Rudin": **READ, correct.**
* ref8L_meaning line 95, "13, 14 and 15 are shown equivalent to one another": acceptable. CG assert 13 ⟺ 14 and
  (5.1) ⟺ 15, and they sketch 15 ⟹ 14.
* CG nowhere derive Rudin, Conj. 9 or Conj. 10 from Vojta, Bombieri-Lang or abc. §9 uses abc and Baker's (9.1) only
  for the upper bound a + b ≪ exp(k^{9/2+o(1)}). **READ.**

Files: `lit8L_cg/dl/abs/abs.html` (arXiv abs page), `lit8L_cg/dl/cg_v1/v1.pdf` and `v1.txt`, `lit8L_cg/img/p-07.png`,
`p-10.png`, `p-15.png`, `lit8L_cg/scripts/t10.py`.
