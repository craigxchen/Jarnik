# Round 31 referee report, block "window": astra_lemma_II Part II, sections II.4 and II.5

Scope: the six-norm short-interval lemma and its use through (II.5) to get exponent 4; the congruence (II.6)-(II.7); and in II.5 the normalization and specialization, the literal wording, the moments (II.8), the no-pencil theorem, the reconstruction from four blocks, the generalized-circle exclusion, and the mod-7 computation. I did not referee the two-moment slice, the quadric-determinant paragraph or the fixed-plane construction at the end of II.5. They are outside this block, and I mention them only where they bear on it.

Referee scripts are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref2_window/`. All runs used `prlimit --as=3500000000`.
- `window_ref.py`: symbolic check of Q, plus an exhaustive search for small-|t| representations (output in `window_runs.txt`).
- `window_sharp.py`, `r2zero.py`, `r2sign.py`: an attempt to build sharpness examples for the exponent 1/6.
- `gauss_ref.py`: exact checks of (II.5), (II.6), (II.7), (II.9) and (8) on 3000 random block tuples; the line-graph matrix facts.
- `modp_ref.py p KMAX [--loose]`: my own finite-field enumeration. It is a different algorithm from Astra's `tenpair_modular.py`, and it includes positive controls. Output is in `mod7_kmax_runs.txt` and `mod7_runs.txt`.

Notation as in Part I §6. There are ten edge blocks β_ij with norms n_ij, and P_i = ∏_{j≠i} β_ij. Write [A,B] = Im(conj(A)B), t_ij = [P_i,P_j]/n_ij, u_i = P_i/conj(P_i), and H = ∏ n_e.

## Summary table

| claim | verdict |
|---|---|
| II.4 six-norm lemma | CORRECT (checked by hand; no counterexample found even with much wider windows) |
| II.4 application, exponent 4 when the other four norms are ≍ X | CORRECT |
| II.4 "pointwise gain over 10/3" with no hypothesis on the other four norms | CORRECT. Astra gives no exponent; I supply one: g ≥ 17/5 |
| (II.5), (II.6), (II.7) | CORRECT (checked by hand and with exact random tests) |
| II.5 normalization (monic) and specialization (exponent 10/3 template) | CORRECT |
| II.5 literal wording (what t_ij ≠ 0 excludes) | CORRECT |
| (II.8) moments | CORRECT. At d = 1 it is the record's statement (cotangent.md §3, line ~264) |
| no pencil in any degree | CORRECT |
| reconstruction from four base blocks | CORRECT |
| no generalized circle (d = 1) | CORRECT |
| mod-7 exclusion | CORRECT as a computation, reproduced independently. What it proves is very weak (see §7) |

No claim in the block is false and none has a gap. All of it is peripheral to the goal (§8).

## 1. The six-norm lemma (II.4)

**Proof, line by line.**
- Reduce X to the least integer among the x_ν. The window only shrinks relative to X.
- Write F(Y) = ∏(Y + k_ν), with one k_ν = 0 and 0 ≤ k_ν ≤ K.
- Q is the truncated polynomial square root. Matching the Y^5, Y^4, Y^3 coefficients of Q^2 against e_1, e_2, e_3 gives exactly Astra's a, b, c. I confirmed with SymPy that 16Q ∈ Z[k][Y], that deg_Y(F - Q^2) = 2, and that the Y^2, Y^1, Y^0 coefficients have total degree 4, 5, 6 in k.
- If |t| ≤ c_1 √X, then |s - Q(X)| ≤ (|R(X)| + t^2)/Q(X) ≪ K^4/X + K^5/X^2 + (K^6 + t^2)/X^3. With K ≤ c_0 X^{1/6} this tends to 0.
- 16s and 16Q(X) are integers less than 1 apart, so s = Q(X).
- Because one k_ν = 0, X divides F(X) = Q(X)^2 + t^2. Also 16Q(X) ≡ 16Q(0) mod X. Hence X divides (16Q(0))^2 + (16t)^2.
- This integer is positive (t ≠ 0) and is ≪ K^6 + 256 t^2 ≤ (C c_0^6 + 256 c_1^2) X < X. Contradiction.

The quantifiers are absolute, and nothing depends on a template. Correct.

**Search for counterexamples (`window_ref.py`).** The search is exhaustive over every multiset {k_ν} ⊂ [0,K] containing 0, for every X in the given range. A t with |t| < X^{3/2} can only come from s ∈ {⌊√F⌋, ⌊√F⌋-1, ⌊√F⌋-2}, so I tested those s exactly.

| X range | K | smallest \|t\|/√X found |
|---|---|---|
| [10^3, 3·10^3) | X^{1/6} | none below ~X |
| [10^5, 10^5+300) | X^{1/6} | none below ~X |
| [10^6, 10^6+40) | X^{1/6} | 2·10^6 |
| [10^3, 1300) | X^{1/4} | 357 |
| [10^4, 10^4+60) | X^{1/4} | 2·10^4 |
| [10^3, 1100) | X^{1/3} | 136 |
| X ≤ 500 | X^{1/2} | 5.97 |

So there is no counterexample even with constant 1 and windows far wider than X^{1/6}.

**Sharpness of the exponent 1/6.** I tried to build sharpness examples from the case r_2 = 0, where R = F - Q^2 = r_2Y^2 + r_1Y + r_0, so that R would be linear. Every integral zero of r_2 with K ≤ 30 that I found (641 of them) has F a perfect square and R ≡ 0, so it gives nothing. The exponent 1/6 is an artefact of the proof, and I could not show it is sharp. This does not affect correctness.

## 2. The application via (II.5)

- (II.5): s_ij^2 + t_ij^2 = |P_i|^2|P_j|^2/n_ij^2 = ∏_{k≠i,j} n_ik n_jk, and |u_i - u_j| = 2|t_ij|/√D_ij. Correct, and checked exactly on 3000 random tuples.
- The six factors of D_ij are the six edges adjacent to ij. When they lie in the window, the lemma gives |t_ij| ≫ √X, and √D_ij ≍ X^3, so Δ ≥ |u_i - u_j| ≫ X^{-5/2}.
- With the other four norms (the edge ij and the three edges of the complementary triangle) ≍ X, H ≍ X^{10} and Δ ≫ H^{-1/4}. Exponent 4, with no exceptional set. Correct.

**Supplement: the exponent with no hypothesis on the other four norms.** Astra asserts a gain over 10/3 but gives no exponent. Here is one.
- Write log H = λ log X. The (1/2)log D_e over all ten pairs sum to 3 log H.
- Pair ij contributes 3 log X, so some other pair has (1/2)log D ≤ (λ - 1)/3 · log X. This gives Δ ≫ X^{-(λ-1)/3}.
- The lemma gives Δ ≫ X^{-5/2}.
- So 1/g ≤ min(5/2, (λ-1)/3)/λ. The maximum of the right side is 5/17, at λ = 17/2.
- Hence **g ≥ 17/5 in all cases.**

## 3. The congruence (II.6)-(II.7)

- For any three plane vectors, [P_j,P_k]P_i - [P_i,P_k]P_j + [P_i,P_j]P_k = 0.
- Take squared norms and use |P_i|^2 = ∏ n_i·. Dividing by n_ij n_ik n_jk gives exactly (II.6). I re-derived it by hand and checked it exactly, with signed t, for all triangles and orientations.
- Multiplying (II.5) by 4t_ik^2 t_jk^2 gives 4t_ik^2 t_jk^2 D_ij = A^2 + 4(t_ij t_ik t_jk)^2, which is (II.9).
- If b | D_ij, then A(n) ≡ A(r) mod b, so b divides A(r)^2 + 4(t_ij t_ik t_jk)^2. This number is positive, and |A(r)| ≤ 3T^2K^3, so b ≤ 9T^4K^6 + 4T^6. This is (II.7).
- I checked the divisibility exactly in 360,000 instances, with b a single norm or a product of norms.
- The exponent min(1/6, (1-6δ)/4) is correct.
- Note that only the nine norms appearing in A need small residues.

## 4. II.5: normalization, specialization, literal wording

**Normalization.** A constant t_ij forces the degree-8d leading term of conj(P_i)P_j to be real. So the leading coefficients of the P_i have a common argument mod π, and rescaling multiplies each t_ij by a nonzero real constant. Correct.

**Specialization.**
- The z_i = P_i ∏_{e∌i} conj(β_e) are polynomials.
- Pairwise Bézout bounds the gcds of the twenty factor values by fixed constants c_fg.
- Every factor is absent from some z_i, so the common Gaussian gcd is bounded. At a resultant prime, at most one factor exceeds the fixed valuation bound, and some z_i misses that factor.
- Δ ≍ x^{-6d} and H ≍ x^{20d}, so the exponent is exactly 10/3.

The caveats are correctly stated: a solution over C gives no rational points, and a single curve can always be absorbed into a fixed exceptional set.

**Literal wording.** The argument is that B_ij = ∏_{k≠i,j} conj(β_ik) β_jk exactly (n_ij divides factorwise, so nothing cancels), and B_ij - conj(B_ij) = 2i t_ij ≠ 0. I re-derived each case:
- a real root in a block: excluded;
- a block overlapping its own conjugate: excluded;
- two different oriented blocks sharing a root: excluded (choose i ∈ e∖f, j ∈ f∖e);
- conjugate overlap of incident blocks: excluded (shared vertex v and an outside vertex w).

I also checked the converse claim. For disjoint e, f, no pair ij puts β_e in one of B_ij, conj(B_ij) and conj(β_f) in the other, so this argument cannot exclude that case, exactly as Astra says. Correct.

## 5. (II.8), and comparison with the record

- B_ij is monic of degree 6d with Im B_ij = t_ij. That is equivalent to e_1, ..., e_{6d-1} of its roots being real. Newton's identities are triangular with rational coefficients, so this is equivalent to p_1, ..., p_{6d-1} being real.
- Im p_k(B_ij) = Σ_{h≠i,j}[Im p_k(β_jh) - Im p_k(β_ih)] = S_j - S_i, because the shared edge cancels.
- The count of 20d real unknowns against 4(6d-1) equations is stated honestly as not deciding anything.

Correct.

**Novelty.** At d = 1 this is the record's statement. `growth/cotangent.md` (B4)-(B5) and "Polynomial families" (lines ~259-268) take C_T = α_T(n + γ_T) and require "Im Σ_{T∋i} γ_T^k independent of i for k = 1..5". The record counts this as 20 equations in 18 effective unknowns and reports 450 Levenberg-Marquardt runs, all degenerate. Astra cites this correctly. The general-d form and the monic reduction are routine extensions.

## 6. No pencil; reconstruction; generalized circle

**No pencil.** I checked the proof line by line.
- (II.9) is a polynomial identity.
- A 1-dimensional span forces all n_e = G, and then the left side of (II.9) is (const)·G^6. It cannot equal a positive constant.
- In the 2-dimensional case the binary sextic is nonzero. Each complex linear factor evaluated at (p, q) is a nonzero constant, so 1 lies in the real span, and monicity gives n_e = G + c_e.
- (II.6) then gives s_ij = S(G). Since G is transcendental, (II.5) becomes ∏(Y + c)(Y + c') = S(Y)^2 + t^2. At Y = -c the left side is 0 and the right side is positive.
- The F_e(G) generalization follows the same way: every edge is adjacent to some pair.
- The argument uses only (II.5) and (II.6), so it is valid under the literal wording.

CORRECT. Novelty: "pencil" in the record never refers to this. New, but elementary. Context: the integer fixed-plane family n_e = k^2 + c_e^2 at the end of II.5 is an integer pencil, and this theorem shows it has no bounded-cofactor polynomial version. That is consistent with II.6.

**Reconstruction from four base blocks.**
- D ↦ [U_i, D] is a real-linear map between spaces of real dimension 6d. It is injective, because U_i and conj(U_i) are coprime and deg D < deg U_i.
- Necessity: C_i - U_i has degree < 3d, and [U_i, C_i - U_i] = t_5i, so C_i - U_i = t_5i D_i.
- Sufficiency: the three disjoint degree-d gcds exhaust the monic cubic C_i.
- I checked the asymptotic: arg(P_j/P_i) = (h_j - h_i)x^{-6d} + O(x^{-6d-1}) and |P_i||P_j|/n_ij ~ x^{6d}. So the polynomial [P_i,P_j]/n_ij tends to h_j - h_i and is therefore constant.

CORRECT. My own count for d = 1: 12 real parameters, 12 real conditions, and a 2-dimensional real affine symmetry, so the system is overdetermined by 2. This matches the record's "20 equations in 18 effective unknowns". As Astra says, it decides nothing.

**Generalized circle.**
- Nonhorizontal line: G_ij(a) is real at the real intersection point a.
- Horizontal line: from Im p_2, p_3, p_4 (degree ≤ 5 = 6d - 1), Newton gives equal abscissa triples on disjoint edge sets, contradicting distinctness.
- Circle through a, b on the real axis: by inscribed angles, the ratio G(a)/G(b) is real. Both values have the same nonzero imaginary part, so the ratio is 1. Then M·log|w| = 0, so all roots lie on the perpendicular bisector.
- Tangent circle: G'(a) is real and G(a) is not, so G'(a) = 0. Then M·Re(1/(a - r)) = 0, so all roots lie on a vertical line.
- Apollonian case: |G(a)| = |G(ā)| forces Im R(a) = 0 (since t ≠ 0), so G(a) = G(ā). Then M·arg w ≡ 0 mod 2π, so every w/k is a sixth root of unity, and six values cannot hold ten distinct Möbius images.

I verified the matrix facts with `gauss_ref.py`: det M = -192, eigenvalues 6, 1, -2 with multiplicities 1, 4, 5, and 6M^{-1} = 3I + 3M - 2J, which is integral with entries 1 / 1 / -2.

CORRECT. Novelty: new. Significance: tiny. Ten points on a circle is a codimension-7 condition.

## 7. The mod-7 computation: what it does and does not prove

**Reproduction.** My `modp_ref.py` uses a different algorithm from Astra's. It buckets stars by their moment vector, finds 5-cliques under exactly-one-intersection and conjugate-compatibility, then re-verifies the edge labelling and the moment equations from scratch.
- At p = 7, KMAX = 5: 95,760 stars, 16,660 classes, largest class 14, as Astra states. **Zero solutions**, both strict and in the loose variant (conjugate collisions allowed on disjoint edges).
- Astra's script, rerun, also gives 0 / 0.

**Positive controls (mine).** I counted solutions with fewer moments. The random-model expectation is (number of unordered configurations ≈ 1.09·10^13) / 7^{4·KMAX}.

| KMAX | expected | found |
|---|---|---|
| 3 | 789 | 1638 |
| 4 | 0.33 | 0 |
| 5 | 1.4·10^{-4} | 0 |

So the code does find solutions when they should exist. It also shows that **emptiness at KMAX = 5 is exactly what a random system predicts**. The computation carries essentially no structural information. This is unsurprising: the system is overdetermined by 2 modulo the affine group.

**What it proves.** Take a degree-one solution with roots in Q(i). The prime 7 is inert, and conjugation reduces to Frobenius. Suppose a real rational affine change of variable makes the roots 7-integral. Then the twenty reduced roots and conjugates cannot all be distinct in F_49. This is Astra's statement, and it is correct.

**What it does not prove.**
- It is not 7-adic insolubility.
- It says nothing about solutions over R or C, or over number fields where 7 does not stay inert of degree 2.
- It says nothing about degree d ≥ 2.
- It does not exclude rational solutions: bad reduction at 7 is a positive-density local condition, which large-height configurations meet freely.
- t_ij ≠ 0 is not used mod 7.

Astra's own wording ("a restriction on reduction, not a 7-adic nonsolvability claim") is accurate.

## 8. Significance (blunt)

- **Six-norm lemma.** Correct, but it rests on an additive window of relative width X^{-5/6}. Neither the cut balance of the growth problem nor any multiplicative information forces such a window, and Astra concedes this. It is peripheral: it changes neither the growth rate (2/3+ε) log R/loglog R nor the uniform bound, and it does not change the status of Question A, Question B or route 5. The same applies to (II.7). It is a correct and slightly more flexible congruence version, but the hypothesis "the norms have small residues mod a divisor b of D_ij" is again not supplied by the problem.
- **Polynomial (function-field) side.** The central question is still open: does the ten-pair system have a nondegenerate solution in some degree? A solution would be a genuine 10/3 template outside the known components. Neither Astra nor the record decides it, even at d = 1. What II.5 adds:
  - an exact, correctly argued normalization and reconstruction;
  - one clean all-degree exclusion (no pencil; more generally, no common G with real-rooted F_e);
  - two d = 1 exclusions (generalized circles; good reduction at 7). Both are codimension or local conditions and are far from a decision.
- **Comparison with the record.** The record's numerical evidence (450 LM runs, all degenerate) for d = 1 is at least as informative as everything here.
- **Recommended next step.** A rigorous decision of the d = 1 system: 4 complex base roots and 4 real h_i, with six cubic-resultant conditions, after quotienting the real affine group. This is a small elimination problem; for example, use Gröbner bases or certified homotopy continuation on the 12-parameter system with the gcd/noncollision saturation. That would settle the record's open polynomial test at d = 1.
