# Referee report, round 31, block "gcd": Part VII, section VII.2, paragraph "PROVED (the gcd estimate)"

Source: `astra_lemma_VII/answer.md`, lines 1597-1614.
Scripts: `round31/ref7L_gcd/` (gcdlib.py, exp.py, exp2.py, exp3.py, resultant_check.py).

## Verdict

**CORRECT.** I checked every step line by line. The constants can be made explicit. The bound is uniform in the
coefficients. I tested the resultant step symbolically and the size of the gcd numerically. There are only
cosmetic corrections, listed in section 3. One extra point: the argument never uses the characteristic, so the
lemma holds over any field. That matters for how VII.6 is read: whatever fails in characteristic 2, it is not
this lemma.

Precise form of the statement that is proved. Fix r, a box B ⊂ Z^r of width w_B, and a sequence n_l with
L_l = ||n_l||_∞ → ∞ and, for every fixed u ≠ 0, n_l·u ≠ 0 for all large l. Then for every η > 0 there is
l_0(η, B, (n_l)) such that for all l ≥ l_0, all k, and all lists F_1..F_k ⊂ K[Z^{±1}] (K any field) with support in
B and no common nonunit factor:

  deg gcd_{K[t^{±1}]}(F_1(t^{n_l}),…,F_k(t^{n_l})) ≤ η L_l.

Concretely, for each Q, once l ≥ l_1(Q, B, (n_l)) the left side is ≤ 2 r² w_B² L/Q. The threshold depends on the
coefficients in no way.

## 1. Line-by-line check

**(a) r = 1.** In one variable, "no common nonunit factor in K[Z^{±1}]" means no common root in K̄^*. A common
root t of the F_i(t^n) would make t^n a common root. So the common vanishing is 0, with multiplicity. Correct.
(The case k = 1 with r ≥ 2 is also trivial: F_1 is then a unit, i.e. a monomial.)

**(b) Dirichlet.** Some entry of n/L equals ±1. Take j = 0..Q^{r-1}, which is Q^{r-1}+1 values, and put the
fractional parts of the other r−1 entries of j n/L into Q^{r-1} cubes of side 1/Q. Two values j₁ < j₂ land in the
same cube. Then b' = j₂ − j₁ ∈ [1, Q^{r-1}] and some v ∈ Z^r satisfy ||v − b'n/L||_∞ < 1/Q. The ±1 coordinate is
exactly the integer ±b'. In fact ||v||_∞ ≤ b', which is better than the stated b'+1. Correct. v ≠ 0, and v ranges
over a finite set that depends only on Q.

**(c) Generic A, C.** For each fixed l the coefficients are fixed. The α for which A = Σα_iF_i has support smaller
than U = ∪ supp F_i lie in finitely many hyperplanes. A has finitely many irreducible factors P_j, and P_j | Σγ_iF_i
is a proper linear condition on γ because P_j does not divide every F_i. Avoiding finitely many proper subspaces
gives C. Over an infinite field such A, C exist for every l. "Uniformly generic" is not needed: everything later
uses only supp A = supp C = U ⊂ B and (A, C) ⊂ (F_1..F_k). So the choice may depend on l without any cost. Full
support is used only to make the X-degrees of the two stripped polynomials equal, with nonzero leading
coefficients (the max-v-weight parts). It is convenient, not essential.

**(d) Substitution.** φ: Z_j ↦ Z_j X^{v_j} is an automorphism of K[Z^±, X^±]. A and C are coprime in K[Z^±] and are
free of X, so they are coprime in K[Z^±, X^±], and hence so are φ(A) and φ(C). Put Ã = X^{-a}φ(A) and
C̃ = X^{-a}φ(C), with a = min_{u∈U} v·u. Both lie in K[Z^±][X], neither is divisible by X, and both have X-degree
d = max v·u − min v·u ≤ ||v||_1 w_B ≤ r b' w_B. Suppose a common nonunit factor exists in K[Z^±][X]. It is either a
nonunit upstairs, which is impossible, or c Z^u X^m with m ≥ 1, in which case X would divide both, also
impossible. Gauss's lemma over the UFD K[Z^±] then gives coprimality in K(Z)[X]. Correct.

**(e) Resultant.** Both polynomials have exact degree d, so the Sylvester resultant R ≠ 0. Sylvester gives
R = P·Ã + S·C̃ with P, S ∈ K[Z^±][X]. Setting X = 1 gives Ã(Z,1) = A(Z) exactly: the substitution is undone and
X^{-a} becomes 1. So R ∈ (A, C) ⊂ (F_1..F_k) in the Laurent ring, and the stripped powers of X cause no trouble.
Every Sylvester entry has support in U, so supp R ⊂ 2d·B, which has width ≤ 2 d w_B. That is O(b') with a constant
depending only on (r, B). v-homogeneity: Ã(ZS^v, X) = S^a Ã(Z, SX), and Res(P(SX), Q(SX)) = S^{d²} Res(P,Q) for
degree-d polynomials. Hence R(ZS^v) = S^{d²+2ad} R(Z), so every monomial of R has v-weight d² + 2ad. **Checked
symbolically** (`resultant_check.py`, sympy, r = 2, 3, several v including negative entries): R ≠ 0; the
observed v-weight set is the single value d²+2ad in all 8 trials; and R reduces to 0 modulo a Gröbner basis of
(A, C), so membership holds already in the polynomial ring. For d = 0, R = A has all of the needed properties.

**(f) Nonvanishing of R(t^n), uniformly.** For fixed Q we have b' ≤ Q^{r-1}, so d ≤ r Q^{r-1} w_B. Every possible
exponent of every possible R lies in a fixed finite set E_Q = ∪_{d ≤ d_max} 2d·B. The difference set D_Q = E_Q − E_Q
minus {0} is finite. The hypothesis gives l_1(Q) beyond which n·w ≠ 0 for all w ∈ D_Q. After that, distinct
monomials of any admissible R have distinct n-weights, so R(t^n) ≠ 0. **Yes, the tail is uniform**: it is the
maximum of finitely many thresholds, and it depends on (n_l), B and Q only. It does not depend on the
coefficients, on v, or on b'.

**(g) Zero count with multiplicity.** For exponents u, u' of R we have v·(u−u') = 0. So

  |n·(u−u')| = |(n − Lv/b')·(u−u')| ≤ (L/(b'Q))·||u−u'||_1 ≤ (L/(b'Q))·2 r² b' w_B² = 2 r² w_B² L/Q.

R(t^n) is therefore a nonzero Laurent polynomial in t whose exponent spread is ≤ 2r²w_B² L/Q. Its zeros on K̄^*,
counted with multiplicity, number at most that spread. R(t^n) = P(t^n,1)A(t^n) + S(t^n,1)C(t^n) is an identity in
K[t^±]. So g = gcd(F_i(t^n)) divides gcd(A(t^n), C(t^n)), which divides R(t^n), and deg g ≤ 2r²w_B²L/Q. This bound
is **with multiplicities**. Correct.

**(h) Independence of Q.** The constant 2r²w_B² does not involve Q, and Q enters only through the threshold
l_1(Q). So limsup_l deg g / L ≤ 2r²w_B²/Q for every Q, which gives o(L). This is uniform over all coefficient
choices and all k, in the quantifier order stated above. Correct.

**The hypothesis is sharp.** For n = (l, l+1, l+2), u = (1,−2,1) has n·u = 0. Take F = Z_1Z_3 − Z_2², G = Z_1 − 1.
These are coprime, but F(t^n) ≡ 0, so the gcd is t^l − 1, of degree L. A fixed direction n = m(p,q) gives
gcd(t^{mp}−1, t^{mq}−1) of degree m ≍ L.

## 2. Standing cautions

* *Uniform exceptional sets from fixed-S theorems:* not used. The lemma is elementary and uses no Subspace or
  Bombieri-Zannier input.
* *Gcd / denominator costs:* this lemma **is** the gcd-cost estimate, and it counts multiplicities. It controls
  only the loss for lists with no common factor on the formal torus. Turning this into "primitive height
  W = d_F + d_G − o(·)" happens in the next paragraph, (VII.3), which is outside this block. There the lists are
  made coprime after the lattice quotient, and n is re-expressed on the subgroup generated by support
  differences. In that paragraph the estimate is applied with L ≍ d_F + d_G, which the record must check
  separately.
* *Drifting phase:* not relevant here.
* *Template-dependent thresholds:* the threshold l_1(Q) depends on the sequence (n_l). In VII.2 the sequence is a
  diagonal subsequence of a hypothetical counterexample sequence, chosen before the lemma is applied. Using
  "eventually" along that fixed subsequence inside a proof by contradiction is legitimate. The final
  exceptional set in VII.2 does not come from this threshold. It comes from one frozen family at one n_*, which
  is outside this block. The lemma gives no rate: see the numerics, where the decay can be as slow as
  O(1/log L). Any downstream use that needs a rate, rather than "eventually", would be a GAP. VII.2 does not
  need one.

## 3. Cosmetic corrections

1. ||v||_∞ ≤ b' (not b'+1). Harmless.
2. Explicit constant: the common vanishing is ≤ 2r²w_B² L/Q for l ≥ l_1(Q, B, (n_l)).
3. "Generic" A, C means: chosen separately for each l, from a nonempty Zariski-open set of α, γ. No uniformity
   is claimed or needed.
4. The proof works over every field, including characteristic p: Sylvester, Dirichlet and counting zeros by
   exponent spread are all characteristic-free. So VII.6's failure in characteristic 2 has to come from the
   covering inequality (VII.2), as Astra herself says. It cannot come from this lemma.
5. The inline citation of Amoroso-Leroux-Sombra is decorative. The argument does not use it.

## 4. Numerics (worst observed common degree against L)

All computations are exact: gcds over F_P with P = 2147483629, plus sympy over Q for the resultant step. Random
coefficients mod P stand in for generic complex coefficients. The mod-P gcd degree is an upper bound for the
characteristic-0 gcd degree of integer data.

* r = 2, F = Z_1−1, G = Z_2−1: gcd degree = gcd(a,b) on (12,18), (100,75), (97,1000), namely 6, 25, 1. As
  expected.
* r = 3, box [−2,2]³. For each n I computed the worst binomial pair, max over independent w, w' of
  gcd(n·w, n·w'), which a planted pair realizes. I then planted it in random moving-coefficient pairs
  F = a_1(Z^w−1) + a_2(Z^{w'}−1), G likewise, with the a_i random 2-term polynomials re-drawn for each n. The
  measured gcd degree equalled the binomial prediction every time (`exp.py`). The coefficients did not create
  anything beyond the torsion-coset contribution.
  - Generic sequence n = (l, ⌊l√2⌉, ⌊l√3⌉): the worst ratio deg/L is about 0.78-1.21 at L ≤ 2·10³. The ratio can
    exceed 1 because n·w ranges up to 5L. The maximum over windows l ∈ [l_0, 2l_0] is 1.53, 1.21, 0.21, 0.061,
    0.0185 for l_0 = 10², 10³, 10⁴, 10⁵, 10⁶ (`exp2.py`). This decays roughly like L^{-1/2}, as expected from
    shortest orthogonal vectors of size ≍ √L.
  - Adversarial sequence satisfying the hypothesis, n = (m, mp, mp−m/2+1) with p ≈ log₂L: the worst ratio is
    0.44, 0.38, 0.38, 0.21, 0.22, 0.23, 0.14 for L = 10³…10⁹. It comes from pairs supported on (n_1, n_2), with
    gcd ≍ c·m = cL/p and c ≤ 6. So it decays only like 1/log L. **o(L) is all one gets, and it can be
    arbitrarily slow.**
* Moving coefficients tuned to n (`exp3.py`): supports of size s in [0,3]³, with coefficients solved so that
  F(t^n) and G(t^n) vanish at t = 1 to the maximal order s−1. The gcd degree was s−1 = 3, 5, 7 at L = 200, 400.
  At L = 800 it was 2(s−1), because gcd(n) = 2 doubles the root. The coefficients buy only O(|B|·gcd(n)).
  This is consistent with the lemma: the tuned zeros contribute O(1) per point, and the large contributions
  are torsion-coset (binomial) resonances.
* Hypothesis violated, n = (l, l+1, l+2): the gcd degree is L, as computed in section 1.

No counterexample was found. The worst observed values are consistent with the bound 2r²w_B²L/Q = 72L/Q at
r = 3, w_B = 4 (box [−2,2]), which is weak for small Q.

## 5. Literature and novelty

I could not fetch arXiv, ar5iv, iris.unito.it or eaca2012.web.uah.es: DNS was blocked in this session. So the
comparison below rests on search snippets and my own memory, and I have **not verified** any exact statement.

* Amoroso-Leroux-Sombra, "Overdetermined systems of sparse polynomial equations", FoCM 15 (2015) 53-87,
  arXiv:1307.5788. Per the snippets, they solve univariate lacunary systems in time quasi-linear in log-degree,
  relying on Bombieri-Zannier on multiplicatively dependent points. As far as I know, their structural results
  are for fixed coefficients (algebraic numbers, or a fixed system) and separate out the torsion/cyclotomic
  part. Astra cites this paper but does not use it.
* Filaseta-Granville-Schinzel (Math. Comp. 2008), "Irreducibility and gcd algorithms for sparse polynomials",
  and Bombieri-Masser-Zannier, "Anomalous subvarieties" (IMRN 2007). For fixed coprime F, G over Q̄ these give,
  outside finitely many exceptional subgroups of exponent vectors, a bounded **non-cyclotomic** part of
  gcd(F(x^n), G(x^n)). The cyclotomic part is left to separate elementary arguments. The constants and
  exceptional sets depend on F, G.
* Schinzel, *Polynomials with special regard to reducibility* (CUP 2000), Ch. 6, with the appendix by Zannier:
  the same fixed-polynomial setting.

None of these, as I know them, states a bound uniform over **moving complex coefficients** with fixed support.
That uniformity is exactly what VII.2 needs. On the other hand the target here is only the total o(L), not O(1)
for the non-cyclotomic part. For that weaker target, Dirichlet plus a v-homogeneous resultant is a short,
correct, elementary argument. It belongs to a familiar genre: the "small orthogonal vector" quantity ||n||/δ(n)
appears throughout this literature. I would call it **new in the record** (grep finds no such lemma in
docs/claude_review) and **plausibly folklore-level in the literature**, but I cannot cite a published statement
with moving coefficients. Neither the fixed-coefficient theorems nor BMZ imply it directly, because their
constants depend on the coefficients.

## 6. Significance

The lemma is correct but low-level. Its only role is to control primitive-height loss inside the sparse-class
theorem of VII.2. VII.2 does not touch A_5(δ) for general tuples, as Astra herself says ("could not decide the
unrestricted question for any δ"). It has no integer analogue: the integer counterpart of "zeros of R(t^n)" is
divisibility by a resultant value, which gives no o(L) saving. So for the goal (δ < 1/30 in A_5, or the uniform
Cilleruelo-Granville bound) its significance is essentially nil. It is a sound, reusable piece of the
function-field toolkit.
