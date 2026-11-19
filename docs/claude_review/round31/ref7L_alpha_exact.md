# Referee report (exact computation lens): the five-mark Pell family and alpha^ess_e(-K) <= 4

Verdict: **CORRECT WITH CORRECTIONS.** The claim survives every exact check below. The corrections are about
quantifiers and scope (sections 6 and 7). None of them changes the conclusion.

Scripts and outputs are in
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref7L_alpha_exact/`.
* `build.py` checks the polynomials.
* `spec.py` does the integer specialisation.
* `scan.py` scans T in [-300, 3000].
* `jac.py` computes exact Jacobian ranks.
* The outputs are `spec_out.txt`, `scan_out.txt`, `jac_out.txt` and `assoc_poly_out.txt`, plus one `poly_n*.json` per
  family.

All of this code is my own. All arithmetic is exact: sympy over QQ_I, and Python integers for Gaussian integers.
Floating point is used only to print logarithms. Each run was made under `prlimit --as=3500000000`.

## 1. Definitions used

* **The marks.** A(t) = [[t^2+1, t], [t, 1]] and beta = t + 2i. Write each seed as s_j = p_j/c_j in lowest terms and
  put v_j = (c_j, p_j)^T in Z^2. Then G_j = X_j + iY_j, where (X_j, Y_j)^T = A^n v_j lies in Z[t]^2. The bar denotes
  coefficientwise complex conjugation; t is real.
  * z_0 = beta * Gbar_1 Gbar_2 Gbar_3 Gbar_4.
  * z_jk = betabar * G_j G_k * prod_{l != j,k} Gbar_l, for jk in {12, 13, 23, 14}.
  * Using v_j instead of (1, s_j) multiplies every mark by c_1c_2c_3c_4. This leaves the class [zeta] unchanged, and
    the marks now lie in Z[i][t]. This is the honest clearing of denominators. No further factor is discarded before
    the Gaussian gcd is computed.
* **The record's conventions.** These are taken from `expert_note.md` §3, lines 210–240.
  * A tuple is *primitive* if it consists of distinct lattice points with gcd_{Z[i]} = 1.
  * A_k(delta) asks for a proper closed Z in T_k and an R_0 with the following property: every primitive k-tuple on a
    circle of radius R'' >= R_0 with all mutual distances <= R''^{1/2-delta} has [zeta] in Z.
  * **A_k(delta) has no unit-class condition.** The unit-class restriction appears only in (L_k°).
  * The record's height is H = H_{-K} = R''^2, and its exponent is gamma = 4/(1+2delta). If the largest distance is
    R''^e, then delta = 1/2 - e and gamma = 2/(1-e).
* **Unit-class conventions that I test.**
  * (U1) The primary-prime convention: each pi_p is congruent to 1 mod (2+2i). With it, the unit of an odd zeta is
    zeta mod (2+2i).
  * (U2) A convention-free pairwise test: zeta_a/zeta_b = ±P/Pbar with gcd(P, Pbar) = 1. This holds exactly when
    Re(zeta_a zetabar_b) is odd.
  * (U3) The integrality test for t_ab = |c_ab|/2, where c_ab = (zeta_a - zeta_b)/gcd(zeta_a, zeta_b).
  * Under the record's definition, ε_a = ε_b depends on how the pi_p are chosen: replacing pi_p by i·pi_p multiplies
    ε_a/ε_b by (-1)^{a_a(p) - a_b(p)}. Only the relation ε_a/ε_b ∈ {±1} is intrinsic, and it is exactly (U2).

## 2. Seeds and degrees

* The seed vectors are:
  * S1 = (1, 2, 3, 5), integer;
  * S2 = (1/2, 1/3, 2/7, -5/3);
  * S3 = (3/4, -1/5, 7/2, 0);
  * S4 = (2, -3, 5/6, 11/9).
* Each seed vector is run with n = 1 (D = 9) and n = 2 (D = 17).

## 3. Polynomial level over Q(i)[t] (`build.py`): all claims hold

Each of the eight families has the following properties:
* All five marks have degree exactly D = 8n+1.
* z_a·zbar_a is the same polynomial for all five marks.
* The gcd of the five marks in QQ_I[t] is 1 (degree 0).
* **All ten** differences z_a - z_b have degree exactly D - 4n = 4n+1, that is, 5 for n = 1 and 9 for n = 2.

Pairwise gcds are large. The pair (z_0, z_jk) has a gcd of degree 4n, the other pairs a gcd of degree 4n+1, and the
pair (z_23, z_14) a gcd of degree 1. This does not matter, because only the gcd of all five marks is relevant.

**Bezout constant.** I computed an exact identity sum_i a_i z_i = 1 with a_i in Q(i)[t], and I asserted the identity
after computing it. Let N be the lcm of the denominators of the real and imaginary parts of all coefficients of the
a_i. Then sum_i (N a_i) z_i = N with N a_i in Z[i][t]. **Hence for every integer T, g(T) := gcd_{Z[i]}(z_0(T), ...,
z_14(T)) divides N in Z[i], so |g(T)|^2 <= N^2.** This bound is uniform in T for fixed seeds. It is a proof, not an
observation. The values of N are:

| n | seeds | N |
|---|---|---|
| 1 | S1 | 2576870400 = 2^10·3·5^2·13·29·89 |
| 1 | S2 | 4098593409217114204800 |
| 1 | S3 | 424897165448046125744000 |
| 1 | S4 | 7225640452961037507792000 |
| 2 | S1 | 1135431323520000 = 2^10·3·5^4·13·2777·16381 |
| 2 | S2 | 2.97·10^34 |
| 2 | S3 | 4.06·10^39 |
| 2 | S4 | 2.55·10^40 |

N grows with the height of the seeds. **The gcd bound is not uniform in the seeds.** The argument in section 6 does not
need uniformity in the seeds.

## 4. Integer specialisation (`spec.py`, `scan.py`)

The test points are as follows.
* Every integer T with 2 <= |T|, -300 <= T <= 3000, for all eight families.
* The values T = 10^k for k up to 60, and T in {11, 13, 17, 57, 101, 997, 4321, 99991, -7, -1000}.

At each T the script does the following.
1. It evaluates the five Gaussian integers exactly.
2. It computes their Gaussian gcd g by the Euclidean algorithm, normalised to Re > 0, Im >= 0.
3. It checks that g | N.
4. It divides by g to get the five quotients zeta_a.
5. It checks that the five zeta_a have one norm R''^2, that their gcd is 1, that they are distinct, and that no two
   are associates or conjugate-associates.
6. It computes e = log(max distance)/log R'', e_min = log(min distance)/log R'', and gamma = 2/(1-e).

**Results.**
* **g(T) | N held at every T tested in every family.** The observed gcd norms include the following:
  * n=1, S1: {1, 5, 13, 32, 64, 65, 89, 160, 320, 416, 445, 832, 1157, 2080, 2848, 4160, 5696, 14240, 28480, 74048}.
  * n=1, S2: values from 2 up to 3755600.
  * n=1, S4: values up to 13712800.
  * n=2, S1: {2, 10, 26, 32, 64, 130, 160, 320, 416, 832, 2080, 4160}.
  * The full lists are in `scan_out.txt`.
* The gcd is nontrivial at most T, and it is never discarded. Every quotient tuple has gcd 1, and its five points all
  have the same norm. R''^2 is odd at every T tested.
* **Distinctness fails at isolated T.** Two marks coincide (z_0 = z_12) at T = -2 for S1 (both n), at T = 2 for S2
  (both n), and at T = -14 for S4 (both n). These are integer roots of a difference polynomial. As polynomials, z_a ≠
  u·z_b for every a ≠ b and every unit u. Also z_a/zbar_b is non-constant in Q(i)(t) for all a and b, so
  z_a·gbar = u·g·zbar_b cannot hold identically for any of the finitely many divisors g of N (`assoc_poly_out.txt`).
  Hence equalities, associates and conjugate-associates occur only at finitely many T on each curve.
* **Observed conjugate-associates.** These occur only at |T| <= 222. The full list is in `scan_out.txt`. The largest
  values are T = 222 (S4) and T = 60 (S2).
* **Convergence.** The table gives max-distance exponents e and record exponents gamma. The targets are 1/2 + 1/(2D)
  and 4 + 1/(2n).

| n | seeds | T=10 | T=10^3 | T=10^6 | T=10^12 | T=10^30 | T=10^60 | target e / gamma |
|---|---|---|---|---|---|---|---|---|
| 1 | S1 | .6365 / 5.502 | .5839 / 4.807 | .5695 / 4.646 | .5625 / 4.571 | .5583 / 4.528 | .5569 / 4.514 | .5556 / 4.5 |
| 1 | S2 | .6767 / 6.186 | .6001 / 5.001 | .5737 / 4.691 | .5671 / 4.620 | .5602 / 4.548 | .5579 / 4.524 | .5556 / 4.5 |
| 1 | S3 | .6838 / 6.326 | .5943 / 4.930 | .5817 / 4.781 | .5688 / 4.638 | .5609 / 4.554 | .5582 / 4.527 | .5556 / 4.5 |
| 1 | S4 | .7091 / 6.876 | .6101 / 5.130 | .5834 / 4.801 | .5696 / 4.647 | .5612 / 4.558 | .5584 / 4.529 | .5556 / 4.5 |
| 2 | S1 | .5703 / 4.655 | .5438 / 4.384 | .5366 / 4.316 | .5330 / 4.282 | .5308 / 4.263 | .5301 / 4.256 | .5294 / 4.25 |
| 2 | S2 | .5995 / 4.993 | .5543 / 4.488 | .5393 / 4.341 | .5358 / 4.308 | .5320 / 4.273 | .5307 / 4.262 | .5294 / 4.25 |
| 2 | S3 | .6028 / 5.035 | .5502 / 4.446 | .5437 / 4.383 | .5366 / 4.316 | .5323 / 4.276 | .5309 / 4.263 | .5294 / 4.25 |
| 2 | S4 | .6161 / 5.210 | .5594 / 4.540 | .5446 / 4.392 | .5370 / 4.320 | .5325 / 4.278 | .5309 / 4.264 | .5294 / 4.25 |

* The exponent of the minimum distance tends to the same target, for example 0.5555 (n=1) and 0.5294 (n=2) at
  T = 10^60. So all ten distances are of the same order.
* The ratio max distance / (R''^{1/2}|T|^{1/2}) becomes constant along each residue pattern of g(T). Examples: 5.045
  for S1, 22.450 for S2, 31.965 for S3 and 39.600 for S4, at large T.
* The convergence is from above, at rate O_s(1/log T). The reason is that log(max distance) = (4n+1)log T + O_s(1) and
  log R'' = D log T + O_s(1), where the O_s(1) terms include log|g(T)| in [0, log N].

**Unit classes.**
* Under the primary convention (U1), **no T in the scan gives five points in one unit class**, in any family.
* Under the intrinsic pairwise test (U2), typically 4 or 6 of the 10 pairs pass. All 10 pass only for n=1, S1, at
  1648 of the 3297 values of T.
* The t_ab-integrality test (U3) gives the same counts as (U2) at every T printed.
* So the specialised tuples are in general **not** in one unit class. A_5(delta) does not require one unit class
  (section 1), so this does not affect the claim. It does mean the family says nothing about (L_5°).

## 5. Zariski density (`jac.py`)

* Let Phi: (s_1, ..., s_4, t) -> (arg(z_j/z_0))_{j=1..4} be the map into the real 4-torus T_5(R).
* I computed its exact Jacobian at (S1, t=3), (S2, t=7) and (S3, t=10), for n = 1 and n = 2.
* At each of these six points the 4x5 Jacobian has rank 4. The 4x4 seed block also has rank 4, with a nonzero rational
  determinant, for example -2.25·10^-9 at n=1, (S1, 3).
* Hence the Q-morphism Phi: A^4 x A^1 -> T_5 (off the zero loci) has an image that contains a nonempty open subset of
  T_5(R). It is therefore dominant, because the real points of a proper closed subvariety have dimension <= 3.

## 6. The argument, with quantifiers

Fix n >= 1 and a proper Zariski-closed Z ⊂ T_5.
1. The seeds s in A^4 with Phi({s} x A^1) ⊂ Z form a closed set B_Z. B_Z is proper, since otherwise Phi^{-1}(Z) = A^5,
   contradicting dominance. Rational points are Zariski dense, so there is s in Q^4 \ B_Z.
2. For this s, Phi(s, T) lies in Z for only finitely many T. The tuple fails to be distinct (or contains an associate
   or conjugate-associate pair) for only finitely many T (section 4).
3. For the remaining integers T, divide by g(T). This gives a primitive tuple of distinct points, with
   1 <= |g(T)| <= N_s.
4. Then R'' ≍_s |T|^D and every distance is ≍_s |T|^{4n+1}, so dist([zeta], e) ≍_s |T|^{-4n} and H = R''^2 ≍_s |T|^{2D}.
   The constants depend on s, which is allowed: the McKinnon–Roth constant asks only that dist^gamma·H stay bounded on
   one infinite sequence. With gamma = 2D/(4n) = 4 + 1/(2n), dist^gamma·H = O_s(1) on an infinite sequence outside Z.

**Consequences.**
* alpha^ess_e(-K) <= 4 + 1/(2n) for every n, hence alpha^ess_e(-K) <= 4.
* For delta < 0, choose n with 1/(2(8n+1)) < |delta|. Then max distance <= R''^{1/2-delta} for all large T. So for
  every Z and every R_0 there is a counterexample, and **A_5(delta) is false for every delta < 0**.
* At the boundary value delta = -1/(2D) for a fixed n, the inequality is not reached: e decreases to the target from
  above. This is why each n covers only delta < -1/(2(8n+1)).
* Together with the Cilleruelo–Córdoba lower bound (record: arbitrary points, delta > 1/10, alpha >= 10/3),
  10/3 <= alpha^ess_e(-K) <= 4.

## 7. Corrections to the claim as stated

1. "Primitive" and "distinct" hold for all but finitely many T on each curve, not for every T. Equal marks were found
   at T = -2 (S1), T = 2 (S2) and T = -14 (S4), and conjugate-associates at several small |T|.
2. "Gaussian gcd bounded for fixed seeds" is true, with an explicit Bezout constant: g(T) | N_s. The bound N_s grows
   with the seeds. Gcd norms up to 1.4·10^7 occur for small seeds.
3. The tuples are generally not in one unit class, under either the primary convention or the intrinsic ± test. This
   is harmless for A_5(delta), which has no unit-class clause, but the family gives nothing for (L_5°). The family
   is also unbalanced, since it uses 5 of the 15 cuts.
4. The convergence of the exponent to 1/2 + 1/(2D) is from above and only O(1/log T). The statement "exponent
   4 + 1/(2n)" is correct as a statement about bounded dist^gamma·H. It is not correct as an attained inequality
   max distance <= R''^{1/2+1/(2D)}.

## 8. Relation to the record

* The function-field family is already recorded in `expert_note.md` lines 356–358 ("*Dense families:* polynomial Pell
  families (A = [[t^2+1, t], [t, 1]]) give five marks with D = 8n+1, h = 4n, i.e. W = 4h+2").
* `round20/kill_gauss.md` line 177 records the necessary condition e > 1/2.
* The record does not state the consequence over Q, alpha^ess_e(-K) <= 4, or the falsity of A_5(delta < 0). Its
  sentence "the open window alpha^ess in (15/4, 4) is untouched" (line 334) concerns the lower bound and is consistent
  with this.
* The consequence follows from the recorded dense family by two routine steps:
  * specialisation, with Bezout control of the gcd;
  * the standard rational-curve upper bound alpha <= (degree)/(multiplicity) = 2D/h.
* So the result is new as a statement in the record, but it is a short corollary of a family the record already
  contains. The round-31 referee files `docs/claude_review/round31/ref7L_pell.md` (around line 153) and
  `ref7L_meaning.md` state it.
