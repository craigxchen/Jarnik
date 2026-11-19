# Referee report, block "examples": Part IV, inequality (IV.2) and the 24-mark bound

Scripts and outputs are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref4_examples/`:
`engine.py`, `front.py`, `part_a.py`, `part_b_ff.py`, `part_c.py`, `part_d_trace.py`, `part_e_elliptic.py`,
`part_e_kummer.py`, `count4.py`, `findE.py`, and the `*.out` files. Everything is exact (Python `Fraction`, and
sympy over Q, Q(i) and Q(sqrt(-3))), except the GF(2^d) code, which uses its own exact field arithmetic. I did not
use Astra's scripts.

## 0. Verdict

* **I could not break (IV.2), and I now believe it.** I checked the proof line by line: the normalisation at
  infinity, (IV.3), effectivity, the blow-up ledger, the function phi, the layer step and the final algebra. Every
  step is correct.
  * The one external input is Langer's orbifold BMY: (X,D) log canonical, coefficients in [0,1], K_X+D
    Q-effective, implies (K_X+D)^2 <= 3 e_orb(X,D).
  * I confirmed that statement, and the Q-effective (not nef) hypothesis, only through secondary restatements:
    Pokora (arXiv 1612.05141), the abstract, and a ResearchGate copy. The arXiv text could not be fetched, so I did
    not read Langer's text itself.
  * Every crossing in Astra's boundary has coefficients 2/3 = 1-1/3 and 1/2 = 1-1/2. So (Y,B) is a genuine
    snc orbifold with local groups Z/3, Z/2, Z/3xZ/3 and Z/3xZ/2, and the local orbifold Euler numbers are the
    classical ones.
* **Direct Chern numbers.** On 28 configurations over Q, Q(i) and Q(sqrt(-3)) I computed (K_Y+B)^2 and e_orb
  independently. I used an explicit intersection form on pullback and total-transform classes, wrote out the
  strict transforms, and checked the snc conditions.
  * The direct value of 3e_orb-(K+B)^2 equals Astra's ledger exactly every time.
  * The ledger is at most Astra's bound every time.
  * 3e_orb-(K+B)^2 >= 0 every time in characteristic 0.
* **Characteristic 2 (d).** The trace construction violates (IV.2) for q >= 32: m = 31, 63 with delta = E = 0.
  The quantity 3e_orb-(K+B)^2 on Astra's surface equals (q-1)(4-2q/9) exactly: -868/9 at q=32 and -644 at q=64.
  The failing step is Langer/Miyaoka BMY. Every other step is characteristic-free.
* **Largest m found with delta = 0 and all E = 0 over C[t]: m = 12**, exact over Q, at L = 8.
  * It is the cube on a quadric, equivalently a coset of Z/2 x Z/4 on an elliptic curve in P^1 x P^1.
  * Earlier values were 8 (the record's cube cluster) and 9 (the dual Hesse arrangement, also new here).
* **Can 24 be improved?** Not with this method's weights: a = 2/3 is optimal, and exceptional weights cannot help
  at balanced layers.
* **The most important finding is not in Astra's text.** (IV.2) applies verbatim to the record's function-field
  model (rounds 29-30). It gives **M(D,0) <= 24 for every D**, and **sup_D M(D,c) < infinity for every fixed c**.
  The record lists exactly this uniform bound as OPEN.

## 1. The proof, checked step by step

The notation is Part IV's. I checked each step by hand and, where possible, numerically (sections 2-5).

1. **Sections and intersections.**
   * Finite intersection multiplicity equals the vanishing order of D_ij, since the columns are primitive.
   * On a Hirzebruch surface, C_i ≡ C_0 + c_iF. Hence C_i·C_j = -n + c_i + c_j and 2C_iC_j = C_i^2 + C_j^2.
   * The engine asserts this identity for every pair of every configuration, which validates the contact data.
2. **Normalisation at infinity.** In lattice terms an elementary transformation at an infinity point carrying a
   class c of k sections changes the contact matrix at infinity as follows:
   * pairs inside c lose 1;
   * pairs outside c gain 1;
   * mixed pairs stay 0.
   So S changes by m-2k. Since S = 2I/(m-1) >= 0, the process terminates. The cyclic pairing at distance
   floor(m/2) never puts a pair at a common infinity point when every infinity class has size <= m/2. **(IV.3) is
   correct.**
3. **Effectivity.** K_{Y0}+2F+C_i+C_j is linearly equivalent to the intersection divisor of C_i and C_j times F.
   * After blow-ups, the discrepancy is (1-l_k)E_k^*. The negative terms occur only at shared centres, and there
     are at most mu_t(i,j) of them above t.
   * Each E_k^* <= pi^*F_t, because F's total transform passes through every centre over t.
   * Averaging over pairs gives K_Y + 2F + (2/m) sum C_i' >= 0, and a = 2/3 >= 2/m for m >= 3.
   **Correct.**
4. **The ledger.** I re-derived the change in (K+B)^2, which is -(1+b-b0-ak)^2, and the change in e_orb. This gives
   * F(k,0,0) = 4 - 8k/9;
   * F(2,0,a) = 8/9;
   * F(2,a,a) = -4/9.
   A chain therefore totals (12-4e)/9 <= 4e/9. With general a, a k-centre contributes 4 - a(2-a)k in total. **The
   ledger is correct, and the engine reproduces it from scratch in every example.**
5. **phi and the layer step.**
   * phi(u+v) - phi(u) - phi(v) <= 2uv holds in every case: uv+9 when u,v >= 3; v+2 and 2v+5 otherwise; and 0, 0, 5
     in the small cases (1,1), (1,2), (2,2).
   * The sum over h of E_h equals the sum over pairs of E_ij, since E_ij counts the layers where both rows are
     marked but share no centre.
   * Distinct groups and depths give distinct centres, and phi is nondecreasing.
   * The final algebra reproduces 9 - 3m/8 - (s-m/2)^2/2 per layer, which is (IV.2) after multiplying by 2.
   **Correct.**
6. **Langer.** Langer, Proc. LMS 86 (2003), Theorem 0.1: X normal projective, D = sum a_iD_i with 0 <= a_i <= 1,
   (X,D) lc, K_X+D Q-effective, implies (K_X+D)^2 <= 3e_orb(X,D), with equality forcing K_X+D nef. Here K_Y+B is
   often **not nef**: (K+B)·C = -2/3 on some curves in A2, A3, A4 and cube D=16 below. So the Q-effective form is
   needed, and that is the form quoted.

**Gap: none.** One remark: Astra's `iv_section_check.py` "direct" check for y = u_i - t^2/u_i only evaluates
closed-form formulas, and on an unnormalised surface (two weight-0 blow-ups at infinity). That is a valid but
different surface. It is not an independent Chern computation; mine (section 2) is.

## 2. (a) Direct Chern numbers on explicit configurations (`part_a.py`)

Q = 3e_orb - (K+B)^2 is computed directly; "ledger" is Astra's telescoped formula; "bound" is
(8/9)S + (4/9)I - (4/9) sum phi.

| configuration | contact pattern | S | blow-ups | Q direct = ledger | bound | (IV.2): L, delta, sum E, lhs, rhs |
|---|---|---|---|---|---|---|
| A1: y=u_i-t^2/u_i, m=5 | 2 shared centres at infinity (2 elementary transformations), all finite pairs transverse twice | 10 | 0 | 160/9 | 160/9 | 20, 0, 0, -280, 0 |
| A1, m=9 | same | 18 | 0 | 48 | 48 | 72, 4, 0, -360, 432 |
| A2: 7-fold point with nested groups 7>5>3>3 and a pair chain (e=2) | k >= 5 contact | 27 | 8 | 368/9 | 124/3 | 32, 3, 0, -330, 210 |
| A2': 9-fold point plus 3-groups | k = 9 | 27 | 17 | 52 | 52 | 17, 17/4, 0, -103, 459 |
| A3: (t-1)^2 (t+2)^3 prime powers, a tangency of order 3, 4- and 5-groups | prime powers, chains e=2,3 | 30 | 10 | 376/9 | 392/9 | 40, 3, 0, -494, 162 |
| A4: polynomial sections with a 3-group and a pair chain left at infinity after 3 transformations | contact at infinity | 18 | 8 | 296/9 | 304/9 | 27, 5/4, 0, -679/2, 135/2 |
| A4': rational sections, mixed infinity points | chain e=3 | 24 | 4 | 140/3 | 48 | 55, 0, 0, -1375/2, 0 |

In every case the direct value equals the ledger, the ledger is at most the bound, and Q > 0. The strict
transforms satisfy the snc checks: C_i'·C_j' equals the number of unblown transverse points, and every crossing
with a positive-weight exceptional curve is 0 or 1.

## 3. (b) Near-tight families

* **Columns (c + u_iN, 1), every row marked everywhere (n_i = N).** Here D_ij = (u_i-u_j)N, delta = L/4, E = 0.
  * LHS = (3/4)(m-24)L + m^2L/4 and RHS = m(m+3)L/4. The leading terms m^2L/4 match, and **slack = 18L exactly**
    for every m.
  * Checked for N = (t-1)(t-2)^2(t^2+1) and m = 6, 12, 30: slack 90 = 18·5.
  * The Chern side is Q = 4L = 20, and (K+B)^2 = -5 < 0, so K+B is effective but not nef.
* **Lagrange interpolation with prescribed cuts** (L = 12, random cuts of size 6 and 8, m = 10, 16): delta = 2
  and E = 0, with slacks 360, 640 and 578.
* **The record's function-field clusters, passed to exact star matrices** (`part_b_ff.py`; the columns are
  (Re P_F, Im P_F) with P_F the product of the base roots kept). In every case I checked the cluster condition
  deg(z_i - z_j) <= K exactly.

| cluster | m | L | delta | sum E | (s_h) | (IV.2) lhs / rhs / slack | Q (Chern) |
|---|---|---|---|---|---|---|---|
| cube D=8, c=0 (exact, y = (-9/4,-2/3,1/4,8/3)) | 8 | 16 | 0 | 0 | all 4 | -192 / 0 / 192 | 128/3 |
| cube D=16 (composed with t^2+t) | 8 | 32 | 0 | 0 | all 4 | -384 / 0 / 384 | 764/9 |
| D=12 half-cube(2,1) (code only; roots are certified numerically in the record) | 8 | 24 | 0 | 0 | all 4 (computed from the code) | -288 / 0 / 288 | — |
| vertical-line c=1, D=8, Y=(1,2,3,4,-5,-6,-7,8) (exact) | 14 | 16 | 1 | 0 | all 7 | -120 / 238 / 358 | 56 |

## 4. (c) Search for counterexamples: structured families (`part_c.py`, `part_e_*.py`)

Layers are chosen greedily: at each contact point, the nested chain of largest classes, which keeps E = 0. In
all rows, Q direct equals the ledger.

| family | m | L | delta | sum E | (IV.2) slack |
|---|---|---|---|---|---|
| **dual Hesse (Ceva(3)) 9 lines, 12 triple points, generic pencil vertex** | **9** | 12 | **0** | 0 | 108 |
| Hesse minus the 4 lines through one point (Möbius–Kantor 8_3) | 8 | 12 | 0 | 0 | 120 |
| Hesse 12 lines | 12 | 21 | 5/4 | 0 | 186 |
| B3 (9 lines) | 9 | 13 | 1/4 | 0 | 126 |
| A3 (6 lines) | 6 | 7 | 0 | 0 | 183/2 |
| zero-sum parabolas y = t^4 - prod_S(t-r), R = ±1..±4 | 8 | 8 | 0 | 0 | 96 |
| torsion translates y = zeta_6^i t^2 + t | 6 | 2 | 1/2 | 0 | 36 |
| Fermat type y = (t ∓ w^i)^3 | 6 | 19 | 7/4 | 0 | 333 |
| graphs of the dihedral group D6 in PGL_2 | 12 | 13 | 15/4 | 0 | 600 |
| graphs of the octahedral group S4 in PGL_2 (24 sections) | 24 | 25 | 27/4 | 0 | 2206 |
| **elliptic Z/2xZ/4 coset (cube on a quadric), (1,1)-sections** | **12** | 8 | **0** | 0 | 72 |
| Kummer 16_6 via a Heisenberg orbit on a quadric | — | — | — | — | degenerates: the 16 points lie 4 per ruling (only 4 distinct t) |

Counting argument for line arrangements. A pair with deg D = 1 needs alpha_i + alpha_j >= 2 + L/2. Combining the
pair count over layers with the point-pair count gives m <= min(4L^2/(L-4)^2, L(L-1)/(alpha(alpha-1))). This is 9
at L = 12, so **the dual Hesse arrangement is optimal among line arrangements**. For (1,d)-sections on a (2,2)
elliptic curve, the same argument forces L <= 8.

**No counterexample.** The octahedral group gives m = 24 but needs delta = 27/4. Every finite-group and Fermat
family I tried is far from delta = 0.

## 5. (d) The characteristic-2 trace construction (`part_d_trace.py`, own GF(2^d) arithmetic)

**Incidence check, q = 8, 16, 32, 64.** For every pair, D_ab = f_a - f_b has exactly the common marked x as its
roots. All roots are simple, lie in F_q^*, and number q/4 - 1 = deg D_ab. Hence E = 0 and delta = 0 (with slack
1/4 in the degree bound). Also deg f_a = q/4 - 1.

**Geometry.**
* All m sections share q/4 - 1 centres at infinity, removed by as many elementary transformations.
* There is one (q/2-1)-fold centre over each x in F_q^*, and no other contact.
* Hence Q = (q-1)(4 - 2q/9), and the computed values are:

| q | m | Q = 3e_orb - (K+B)^2 | (IV.2) lhs - rhs |
|---|---|---|---|
| 8 | 7 | 140/9 | -175/2 |
| 16 | 15 | 20/3 | -195/2 |
| 32 | 31 | **-868/9** | **+341/2 (violated)** |
| 64 | 63 | **-644** | **+3717/2 (violated)** |

The ledger equals the direct value even in characteristic 2. So the failing step is exactly Langer's inequality.
The pair (Y,B) is still snc with standard coefficients, and the effectivity argument is still valid; Miyaoka–Yau
fails in characteristic p, and q = 2^5 is an explicit failure on Astra's own surface. The proof cannot be
characteristic-free. For the record's function-field model this means: (IV.2) gives nothing over F_{q^2}(t).

## 6. (e) How large can m be with delta = 0 and E = 0?

**Best found: m = 12, exact over Q** (`part_e_elliptic.py 2 5`).
* E: y^2 = x(x+4)(x+25), with G = Z/2xZ/4 and P0 = (-20, 40) of infinite order.
* T = P0 + G, mapped into P^1 x P^1 by two pencils of lines through -sigma_a and -sigma_b, where sigma_a = -P0 and
  sigma_b = 4P0 + (0,0) - sigma_a. The 8 points lie in distinct t-fibres.
* The marks are the 12 four-subsets of T in the best sum class, each cut out by a Möbius graph.
* Exactly 12 of the 70 four-subsets lie on (1,1)-curves.
* Result: alpha = 4, s_h = 6 for all h, all 66 minors of degree 2, L = 8, delta = 0, E = 0, lhs = -72, slack 72.
  On the Chern side, Q = 16, and the ledger and the bound are both 16.

Equivalently, this is the 8 vertices of a cube on a generic diagonal quadric a x^2 + b y^2 + c z^2 = (a+b+c)w^2,
with the 12 planes containing 4 vertices (6 faces and 6 diagonal rectangles). In characteristic 2 the affine
space AG(3,2) has 14 such planes. Among 8-subsets of groups that embed in an elliptic curve, the maximal number of
4-subsets with a fixed sum is 12 (`count4.py`: Z/8 gives 10, Z/2xZ/4 gives 12, (Z/2)^3 gives 14 but only in
characteristic 2).

**Where the slack sits.** For the 12-mark and dual-Hesse examples, every step after BMY is tight: Q equals
(2/9)×(slack of (IV.2)). So **all of the slack is in BMY itself**.

**Is 24 improvable by weights? Not within this scheme.**
* With section weight a and weight-0 exceptionals at k >= 3 centres, the per-layer constant is
  4 - m(2a - 3a^2/2)/4. This is maximised at a = 2/3 and gives exactly m <= 24.
* An exceptional weight b on an isolated k-centre changes its contribution by b(ak-4) + b^2. This helps only for
  k <= 5 (savings 1, 4/9, 1/9 at k = 3, 4, 5).
* The extremal profile has layers of size m/2 = 12, where b = 0 is optimal.
* At m = 24 every inequality in the chain would have to be an equality, including Langer's, which forces K+B nef
  and a ball-quotient orbifold. Excluding that would give 23. A real improvement needs an inequality stronger than
  BMY.

**The honest range for the true maximum is 12 <= max <= 24.** A 16-mark configuration (a biplane of order 4 cut
by (1,1)-curves on a quadric) would be compatible with (IV.2), with slack 32. The Heisenberg-orbit attempt
degenerates, and I neither found nor excluded such a configuration.

## 7. Novelty and significance

* **Novelty against the record.** The record's only BMY computation (`round4/L7.md`) is on M_{0,7}, a different
  surface and a different question. A ruled-surface orbifold-BMY bound for sections with prescribed contact layers
  is new to the record.
* **The 12- and 9-mark delta = 0 configurations are new data.** The record's best was the 8-point cube; integral
  models with full divisibility were known only at 5 marks (III.11).
  * As a polynomial family over Q[t], the 12-mark example meets full polynomial divisibility at the critical
    degree.
  * Its pair separation is borderline: pairs sharing two points have |Δ| = L/2.
  * Integer specialisations would also need control of extra common integer factors. I did not check this.
* **Unstated consequence (significant for the model).** In the record's function-field model, z_F = Q·P_F/conj(P_F),
  so deg(z_i - z_j) = L/2 + deg D_ij - (alpha_i+alpha_j)/2. The cluster condition is therefore exactly the degree
  bound with delta = floor(D/2) + c - D/2 <= c. Divisibility holds automatically (E = 0), and the columns are
  primitive. I checked this on the exact cube and c=1 clusters: computed delta = 0 and 1, E = 0. Applying (IV.2):
  * **M(D,0) <= 24 for every D.**
  * For c >= 1, apply (IV.2) to any 25 marks: 1.5D <= 700c. So **M(D,c) <= 24 for D > 1400c/3**. With the
    record's M(D,c) <= 2^{2c+2}(D+2) this gives **sup_D M(D,c) < infinity for every c.**
  * This settles the record's open "uniform in D" question (the finiteness half of (R_4)), and does it with a
    BMY-type argument, which is what "truncated function-field Vojta" was expected to need.
  * This holds over C and Q(i), not over F_{q^2}.
  * Since (IV.2) never uses norms, the function-field analogue of Astra's (C) is false in characteristic 0: the
    divisibility-plus-degree formulation alone bounds function-field "integral fakes".
* **Significance for the integer goal: low, as Astra says.** Nothing transfers. The arithmetic analogue of
  Langer's inequality is not available, and an integer parent need not lift to C[t] (Astra's item 11). The growth
  rate is unchanged.

## 8. Verdicts

* (IV.2) and the "at most 24 marks" corollary over C[t]: **CORRECT**, conditional only on Langer's Theorem 0.1 as
  quoted (secondary-source verification). Numerically confirmed on 28 configurations; no counterexample.
* The ledger (F-values, phi, the per-layer constant): **CORRECT**, reproduced by direct intersection theory.
* (IV.3) and the normalisation at infinity: **CORRECT.**
* The characteristic-2 trace data (delta = 0, E = 0, m = q-1): **CORRECT** for q <= 64. It violates (IV.2) for
  q >= 32, and the failing step is BMY.
* Astra's "direct" check for u_i - t^2/u_i: **CORRECT WITH CORRECTIONS.** It is a formula check on an unnormalised
  surface, not a direct Chern computation.
