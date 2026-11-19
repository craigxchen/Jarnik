# Referee report, round 31, block "descent": Astra Part VI.5 (moving affine pencil G_e = a_e B + c_e)

Source: `astra_lemma_VI/answer.md`, lines 1489-1540 (VI.5), together with VI.6 and the summary item 9.
Scripts: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref6L_descent/`.
I did not modify the record.

## 0. Verdict in one paragraph

**VI.5 is correct.** I checked every step of the descent line by line and found no gap. The exact identities are
reproduced with my own code. The one finite input, the "sixteen-point" two-moment lemma of II.5, is reproduced
independently by a Gröbner basis and Hilbert-function computation, not only by rerunning Astra's Jacobian
certificate. The exceptional set really is independent of the moving coefficients: it is the image closure of one
fixed finite-type algebraic family. No Subspace Theorem is used, and the threshold is effective. The weak points are
scope, novelty and significance. The class covers about P^{2+o(1)} of the P^{20} ten-block configurations. The
uniformity half is the record's III.6 for block degree d=1, extended to coincidence strata. The two-moment input is
II.5. The new part is the descent: a polynomial division in a reduced pair word forces the moving Gaussian parameter
B onto a rational line of subpower height. This says nothing about (*) for independent blocks, about A_5(delta), or
about the growth rate.

## 1. Setting, and the exact statement verified

Notation. Edges e of K_5 are the ten pair cuts, oriented towards their two-mark side. Then
z_k = (unit) * (singleton part) * prod_e (conj G_e if k in e else G_e). P = H^{1/20} and Delta = max|z_i-z_j|/sqrt H.

**Theorem (VI.5, as verified).** Fix nu in [0,1). Consider primitive genuine five-tuples of norm H in the
ten-pair-cut profile: every canonical pair block has |G_e| = P^{1+o(1)}, and the singleton blocks and units have
total size P^{o(1)}. Suppose the canonical blocks satisfy G_e = a_e B + c_e with B in Z[i] arbitrary (unrestricted
primes and argument) and a_e, c_e in Q(i) of absolute height o(log P). All twenty coefficients may move.

Let Z be the Zariski closure, in the complexified ratio torus, of the ratio curves of all degree-one constant-cofactor
templates. A template is ten complex linear functions of a real variable T, with conjugates treated as independent,
and five row scalars, such that all ten reduced cofactors [U_iP_i, U_jP_j]/|G_ij|^2 are constant in T. Strict
templates and templates on coincidence strata are both included. Then dim Z <= 2, and Z depends on nothing.

Every such tuple with Delta <= P^{-6+nu} and P >= P_0 has its ratio point in Z. Here P_0 depends only on nu and on
explicit bounds for the o(1) terms (profile errors and coefficient heights), and it is effective. Equivalently: on
this class, Delta <= H^{-1/4-eps} forces the ratio point into Z. The essential exponent on the class is at least 4.

Remarks on the statement.
* The class is defined *inside* the ten-pair-cut profile. Exponent 4 is not claimed for affine-pencil tuples in other
  profiles, and none is needed for (*). In this profile, Delta >= P^{-6-o(1)} automatically, because |t_ij| >= 1. So
  the theorem concerns exactly the window P^{-6} <= Delta <= P^{-5-eps}.
* Z may be empty. Existence of strict d=1 templates is the record's open "strict linear system" (README_III §1.7).
  Section 4(d) below adds numerical evidence that even the weaker four-moment system has no strict real solution.
  This is evidence only.

## 2. Line-by-line check

| Step | Verdict | Notes |
|---|---|---|
| a_e != 0 and \|B\| = P^{1+o(1)} | CORRECT | If a_e = 0 then \|G_e\| = \|c_e\| = P^{o(1)}. Both bounds follow from \|a_e\|, \|c_e\| in [P^{-o(1)}, P^{o(1)}]. |
| Distinct roots r_e = -c_e/a_e | CORRECT | A shared root gives G_e/G_f = a_e/a_f, of height o(log P). |
| No common oriented Gaussian prime between distinct pair blocks; h(G_e/G_f) >> log P | CORRECT | For an oriented split prime, the cuts whose block contains it are superlevel sets of one nested chain, and two distinct 2-sets cannot both lie in a chain. Opposite orientations can share a rational prime, which Astra notes. The ramified prime 1+i divides a primitive z at most once and equally in all five, so it goes into U. No rational prime divides a primitive z. So gcd(G_e,G_f) = 1 in Z[i], and h(G_e/G_f) = log P + o(log P). |
| U_i Gaussian integers of size P^{o(1)}, exact | CORRECT | Explicitly U_k = (singleton block of k) * u_k with u_k in {1, i, 1+i, 1-i}, chosen so that conj(u_k)/u_k is the required unit. No division or denominator occurs. |
| t_ij = [U_iP_i, U_jP_j]/\|G_ij\|^2 nonzero integers, \|t_ij\| <= P^{nu+o(1)} | CORRECT | [Gx, Gy] = \|G\|^2 [x, y], so the division is exact and t_ij = [U_iP_i', U_jP_j'] in Z, where P' has three factors. Nonzero because the points are distinct. \|q_i - q_j\| = 2\|[V_i,V_j]\|/(\|V_i\|\|V_j\|) gives \|t\| = \|V_i'\|\|V_j'\|\|q_i-q_j\|/2 <= P^{6+o(1)} P^{-6+nu}. Checked exactly in `gauss.py` and `descent_identities.py`. |
| t_ij = Im(kappa conj(A(B)) D(B)) with monic root cubics A, D and kappa in Q(i)^x of subpower height | CORRECT | kappa = conj(U_i prod a) * U_j prod a. Verified exactly on random pencils (`descent_identities.py`, 4 trials x 10 pairs). |
| kappa real | CORRECT | The term Im(kappa) Re(conj(A)D) is at least P^{-o(1)} P^{6-o(1)}. Re(kappa) Im(conj(A)D) is O(P^{5+o(1)}), since the arguments of A(B) and D(B) differ by O(P^{-1+o(1)}). This contradicts \|t\| <= P^{nu+o(1)} with nu < 6. A nonzero rational of height h has absolute value >= e^{-h}. |
| Some pair has m = deg(A-D) in {1,2} | CORRECT | A != D because the two edge triples are disjoint and the roots distinct. If m = 0 for all ten pairs, adding back r_ij shows that the first two star power sums agree at all five vertices. The two-moment lemma then allows at most 4 distinct values. Section 3 has the reproduction. |
| Division D = (A-D)Q + R; Im Q(B) = 0 (VI.6) | CORRECT | With kappa real, Im(D/(A-D)) = Im(conj(A-D) D)/\|A-D\|^2 = t/(kappa \|A-D\|^2), because Im(\|D\|^2) = 0. So \|Im(D/(A-D))\| <= P^{nu-2m+o(1)}. Also \|R/(A-D)\| <= P^{-1+o(1)}, and Q has Q(i) coefficients of subpower height, hence common denominator P^{o(1)}. **Bookkeeping: (VI.6) needs only nu < 2m, so nu < 2 suffices for both m.** nu < 1 is used later. Im Q is a nonzero real polynomial in (Re B, Im B), so (VI.6) is a genuine curve. Checked on an exact nu ~ 3 family: every m = 2 pair has real Q, so Im Q(T) = 0 on the real line (`k2_family_test.py`). |
| Linear Q (m = 2): B on a rational line, all ten cofactors constant | CORRECT | Take T = the larger integer coordinate of B, so \|T\| = P^{1+o(1)}. A coordinate that is constant on the line has subpower height. Each f_ij(T) is a real polynomial of degree <= 6 with Q-coefficients of subpower height. If it were nonconstant, \|f(T)\| >= P^{1-o(1)} > P^{nu+o(1)}. **Uses nu < 1, and this is the sharp point of the method.** |
| Quadratic Q (m = 1): completing the square | CORRECT | Im((alpha+i beta)(X+iY)^2) = beta X^2 + 2 alpha XY - beta Y^2, with determinant -\|a\|^2. The factorization beta L_+ L_- has L_pm = X + (alpha pm \|a\|)Y/beta, which I checked by expansion; for beta = 0 the form is 2 alpha XY. |
| k = 0 | CORRECT | Split: two rational lines through z of subpower height, reducing to the linear case. Nonsplit: L_+ = 0 implies L_- = 0 by conjugation, so B = z, which contradicts \|B\| = P^{1+o(1)}. |
| k != 0, E = Q | CORRECT | L_pm are rationals with denominator <= N = P^{o(1)} and product of subpower height. So \|L_pm\| <= N\|c\| = P^{o(1)}, and the inverse change has determinant -2\|a\|/beta of subpower height. This contradicts \|B\| = P^{1+o(1)}. |
| k != 0, nonsplit: (L_+, L_-) = (T, c/T), sigma B(T) = B(c/T), Laurent cofactors constant | CORRECT | The sample T_0 is in E (real) with \|T_0\| = P^{1+o(1)}, and f_ij(T_0) = t_ij exactly, since T_0 is real and E is in R. Nonzero coefficients lie in a real field of degree <= 2 and have subpower height, so they are >= P^{-o(1)} uniformly in E. Astra's "degree at most four" (for E(i)) is harmless. Positive powers are excluded by size (nu < 1). From sigma f(T) = f(c/T) we get sigma(f_n) = c^{-n} f_{-n}, which kills the negative powers. |
| Exclusion of the conic via S(u,v) | CORRECT | The three-vector identity [V_j,V_k]V_i + [V_k,V_i]V_j + [V_i,V_j]V_k = 0, divided by G_ijG_jkG_ki, gives exactly S. Three affine functions of v with distinct roots have a one-dimensional relation space with no zero coefficient. So S = 0 would make the three quadratics in u (root pairs {kl,km}, {il,im}, {jl,jm}) proportional, which is impossible. deg_v S <= 1, while the smooth conic a(u-z)^2 - conj(a)(v - conj z)^2 = 2ik (k != 0) is irreducible of v-degree 2, and S vanishes on all of it. |
| Uniformity: rank of the ratio map <= 2 on every component of the template variety | CORRECT | (i) The template variety W is finite type: 20 + 20 coefficients and 5 + 5 row scalars, with the vanishing of the T^1..T^6 coefficients of the ten cofactors, plus open conditions. (ii) The source derivative is nonzero: q -> 1 at T = infinity, so a constant q would be identically 1 and force t = 0. (iii) Choose an indicator b_{e0} outside the 3-span and project along it. This kills both roots of block e0, leaving at most 18 poles. (iv) The polar parts at a root are proportional to one total valuation vector in all three rows, also on a coincidence stratum, where coincident roots move together. Hence the poles are simple. (v) At infinity, q_j/q_i - 1 = -2i t G conj(G)/(V_j conj(V_i)) = O(T^{-6}), uniformly in the parameters. So d log q = O(s^5)ds and D log q = O(s^6), and the order of vanishing is >= 5+6+6 = 17. (vi) deg div = -2 on P^1, while 17 - 18 = -1 > -2, so the determinant vanishes identically. This is the record's III.6 argument at d=1 (the record's margin "+1 from the source derivative" is the same), now with coincidence strata. |

**Is the exceptional set independent of the moving coefficients?** Yes. Each individual tuple, with its own moving
coefficients, produces a point of W x A^1 whose image is its ratio point. Z is the union of the image closures over
the finitely many components and strata of W, which is a fixed variety. The only sequence- or template-dependent
quantity is the threshold P_0, and it depends only on the explicit height and profile rates, not on any template.

## 3. The two-moment ("sixteen-point") input

* Rerunning Astra's `tenpair_two_moments.py` gives cycle minor 32 and tetrahedral minor 16w+8, as stated.
* **Independent check (`two_moment_independent.py`).** I shift so that the star sums are 0, which puts the r_e in a
  5-dimensional kernel. There are four homogeneous quadrics. Their grevlex Gröbner basis has Hilbert function
  1, 5, 11, 15, 16, 16, ..., constant 16 through degree 12. So the projective scheme is zero-dimensional of degree 16,
  with no positive-dimensional component. In a random affine chart there are 16 solutions, and their numbers of
  distinct values are 2 (6 times) and 4 (10 times). The maximum is 4.
* This confirms the claim: any solution is c + lambda*(one of 16 points), with at most 4 distinct values.

## 4. Numerical tests

(a) **Astra's scope example** (`scope_test.py`). I took G_e = B + c_e with fixed distinct nonreal c_e (c_e != conj c_f)
and real integer B from 10^3 to 10^24, computed exact primitive tuples by Gaussian gcd, and set
nu := 6 + log Delta / log P. The values are nu = 5.43, 5.20, 5.10, 5.05, tending to 5. So Delta ~ P^{-1}, which
matches the first-order claim. Random pencils with real a_e and Gaussian B also tend to nu = 5. With unaligned a_e,
nu is about 6 (Delta ~ 1). Everything is consistent with exponent >= 4.

(b) **Search for low exponent** (`moment_search2.py`). Templates on a real line, G_e = T + c_e, have cofactor
polynomials of degree <= 5-k exactly when the Im star power sums agree for orders 1..k (II.8). Then
Delta ~ P^{-(k+1)}, i.e. nu = 5 - k. I enumerated all integer c_e with |Im c_e| <= 3, |Re c_e| <= 6 (1453 vectors y up
to sign, x normalised). There are **976,152 strict k = 2 templates and none with k = 3.** So in this box the best is
nu = 3, i.e. exponent 20/3 > 4.

(c) **Off the line** (`k2_offline.py`). For one k = 2 template, nu is about 3.17 on the rational lines B = T and
B = T + i. It is 4.43-4.46 for B = T + i*isqrt(T), and 4.93-4.97 for B = T(7+i), where the template's Q are real but
Im B != 0 along the line. This is the mechanism of (VI.6) in action. **No affine-pencil tuple with exponent below 4
was found.**

(d) **Real moment systems** (`four_moment_real.py`, `four_moment_barrier.py`). The gauge is a real affine map of T and
there are 18 real unknowns. k = 3 has strict real solutions (Jacobian rank 12, strictness about 0.1). For k = 4
(the case nu = 1, the boundary of VI.5), 492 of 500 starts converged, but **every converged point is a collision**
(strictness <= 3e-9). With a collision barrier, the residual stalls at 2e-3, against 2e-6 for k = 3. This is
evidence, not proof, that the strict k = 4 system has no real points. If so, the strict part of Z would be empty and
the method might reach nu < 2 (exponent 5) on this class; coincidence strata were not tested. This agrees with the
record's "at p = 7 even four moments have no strict twisted point".

## 5. Standing cautions

* **Fixed-S Subspace exceptional sets:** not used. Z is an algebraic image closure, uniform by construction. This is a
  genuine advantage over VI.1.
* **Gaussian gcd and denominator costs:** handled. The theorem is stated for the canonical (post-gcd) blocks, the U_i
  are exact, t_ij is an exact integer, and the denominators of Q(B) are P^{o(1)} from the coefficient heights. If one
  starts from raw blocks a_eB + c_e, the canonical blocks differ by Gaussian factors of subpower size, so they are
  again of the form a_e'B + c_e' with Q(i) coefficients (as in Astra's resultant remark). This is legitimate.
* **Drifting common phase:** arg B drifts freely and cancels exactly, because every reduced word has degree 3 in B.
  "kappa real" is proved, not assumed. Nothing is replaced by zero.
* **Template-dependent thresholds:** none. The template is produced after the threshold is fixed. The threshold
  depends only on nu and the o(1) rates, and it is effective, which Astra does not claim. Z is not computed and may be
  empty.
* **What the hypotheses cover:** here K = Q(i) and there are no S-units. The hypothesis is one Gaussian parameter B
  common to all ten blocks, with coefficients of height o(log P). The conclusion is effective in its threshold and
  inexplicit in Z. This is honest. It does not cover two or more independent parameters, coefficients of height
  >= c log P, or any unrestricted ten-block configuration. Astra says so in the scope paragraph.

## 6. How much it covers

The ten blocks range over a 20-real-dimensional space (P^{20} lattice assignments of size P). The class is a moving
2-real-parameter slice, B in Z[i] with |B| ~ P, with P^{o(1)} choices of coefficients. That is P^{2+o(1)}
configurations, codimension 18. Generic members already have Delta ~ P^{-1}. Heuristically, members with
Delta ~ P^{-s} should not exist for s > 3/2 except on structured sublocus (count P^{2-4(s-1)}). So the proved
exponent 4 is far from the class's expected truth, though rigorous. Astra's scope paragraph ("Zariski-dense genuine
approximants at slower widths", rank-four first-order map) is correct: (a) reproduces the P^{-1} width, and the
first-order map is c -> (star sums of Im c_e) - (star sum at 1), which has rank 4.

## 7. Novelty and significance

* **Novelty.** The grep finds no affine-pencil or moving-B result in the record. The fixed-plane theorem of II.6-II.7
  (norm vectors on a fixed rational plane, exponent >= 38/11, constants ineffective and U-dependent) concerns a
  different class. The judge's open sub-question README_II §4.1 ("fixed-plane rigidity without a fixed plane ...
  planes of height X^{o(1)}") is close in spirit. VI.5 is a *uniform, effective, moving-coefficient* rigidity result,
  but for the affine one-Gaussian-parameter class, not for norm planes. It does not answer §4.1 as posed.
  * New: the reduction of a numerical tuple to a d = 1 template (division, rational line, conic exclusion).
  * From the record: the uniformity step (III.6 at d = 1; Astra adds coincidence strata) and the two-moment input
    (II.5).
* **Significance for the goal: low.** VI.5 is a correct toy-class theorem. It does not decide (*). It gives no
  A_5(delta) for unrestricted genuine tuples and leaves the growth bound and the uniform arc question unchanged.
  Astra says all of this plainly (VI.6, "remaining step"), and I agree. Its main value is methodological: a moving,
  arbitrary-prime Gaussian parameter can be descended to a polynomial template when every block is affine in it. For
  ten independent blocks there is no common argument at which to divide, as Astra notes.

## 8. Corrections (minor)

1. State that the class lies inside the ten-pair profile, so the window is P^{-6} <= Delta <= P^{-5-eps}.
2. (VI.6) needs only nu < 2m; nu < 1 is needed only in the line case and the conic positive-power step.
3. The threshold is effective; Z is inexplicit and possibly empty.
4. The uniformity paragraph should cite III.6 (record) as the d = 1 case.
5. "Liouville in degree at most four": the coefficients lie in E, of degree <= 2. This is harmless.
