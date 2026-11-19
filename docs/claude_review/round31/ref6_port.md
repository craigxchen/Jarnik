# Round 31, Part VI, block "port": can the ternary construction VI.4 be ported to Z?

Referee/analyst: block `port`. Scripts and outputs are in `ref6_port/` next to this file. All the code is mine:
F_27 uses the modulus a^3+2a^2+1, not Astra's; Gaussian integers, F_9[t] and F_729[t] arithmetic are
hand-written; nothing is imported from Astra's scripts. Labels: **PROVED** means a proof given here or checked
line by line. **EXACT** means a finite exact computation with my own code. **CITED** means a refereed record
result used as a black box. **HEURISTIC** means no proof.

## 0. Verdict in six lines

1. **VI.4 is correct.** I checked the proof line by line. The four-mark certificate reproduces exactly with
   independent code: g=92, full gcd 1038 = 92*11 + 26 (the extra factor is t^26), D=6981, W=13962, and pair
   differences 3471, 3457, 3457, 3471, 3471, 3457. All norms are equal and lie in F_3[t]. The record's cut
   identity holds exactly on it.
2. **It cannot be ported as a polynomial template, and its mechanism cannot be ported at all.** The one step with
   no characteristic-0 counterpart is the existence of one parent P_V with m *different right composition
   factors of the same degree* (P_V = f_H o P_H for every hyperplane H). In characteristic 0, or whenever the
   characteristic does not divide deg f, f o A = f' o A' with deg A = deg A' forces A' = alpha A + beta
   (Lemma 1, elementary, proved below). The characteristic-0 version of the construction therefore yields
   exactly one reflection pair {z, -zbar}. The whole class it lands in (polynomial templates in Z[i][u], with any
   number of parameters under a weighted degree) has at most 24 points by Part IV.
3. **Descent, specialisation and common division all port to Z.** Descent becomes Galois descent from O_K[i],
   specialisation is `round30/bridge.md`, and common division is the record's normalisation. None of them is
   the essential feature. The ultrametric is not essential either. The essential feature is characteristic p:
   additive (Frobenius-linear) polynomials, through which Ritt/Engstrom uniqueness and Newton's identities fail.
4. **Integer computations (EXACT):** Dirichlet-character exponent patterns give only {Z, Zbar} for one character,
   because chi(lambda x) = chi(lambda)chi(x). Several characters, and the literal digit-vector port (F_ell-linear
   forms on base-ell digits), give **width exponents 0.94-1.00** after exact Gaussian gcd removal, over every X
   scanned exactly (X <= 3000) and with the formal gcd up to 2*10^5. A cluster needs 1/2. As polynomial templates
   the points already differ at the first or second odd moment (X-degree n-1 or n-3).
5. **In the record's language** the ternary cluster is of **template type**: a specialised polynomial template
   whose coefficient field F_9(t) has characteristic 3. It is not of covering type. Its integer counterpart class
   is exactly the one Part IV bounds. Everything char-free in the record holds on it verbatim: unique
   factorisation, pair bound, cut identity, chain lemma, Liouville/toric divisibility, small-residue pigeonhole,
   Lemma P. F3+ therefore extends to **F_3 with the reflection form A^2+B^2**, for non-squarefree N.
6. **This is not evidence against the uniform conjecture over Z.** Score for VI.4 on the goal: **2/10**. It closes
   the record's top-ranked odd-characteristic question (README_V §5 item 1) in the negative and kills
   "odd characteristic / Kummer-type reflection" as a route to uniformity. Nothing is proved or suggested for
   integer parents.

## 1. VI.4 checked (prerequisite for porting)

**PROVED (line-by-line check of Astra's proof).**
* Descent: the map alpha -> alpha^9 sends v = i x to i x^(9) (i^9 = i) and preserves 1+chi(eps lambda(x)), because
  lambda has F_3 coefficients and chi(y^9) = chi(y). The fixed field of the map in F_{q^2} is F_9, so F_j is in
  F_9[t][X].
* Conjugation: on F_9 the map alpha -> alpha^q is alpha^3, because r is odd. It sends v to -v (i^q = -i,
  x^q = x), and chi(-1) = -1 because q ≡ 3 (4). So e(-v) = 2 - e(v) and F_j Fbar_j = P_V^2.
* Identity (VI.5): P_H is F_q-linear with kernel H, so P_H(v) = c·eps·lambda(v). Grouping the roots by fibres of
  A = P_H gives F_j = A·prod_{a square}(A - ac)^2 = A(A^s - c^s)^2. Expanding gives P_V + 2c^{q-1}A -
  2c^sA^{s+1}, with P_V = A^q - c^{q-1}A. The pair differences have X-degree <= h(s+1) = (n+h)/2.
* Averaging: the expected minimum exponent is 2^{1-m} if all m forms are nonzero, and 2^{-a} >= 2^{1-m} if
  exactly a < m are nonzero. Hence g >= 2^{1-m}n > n/q = h for q > 2^{m-1}.
* Specialisation: coefficient t-degrees are <= (k-1)(n-g-a) < L, so X -> t^L is injective and keeps the degree
  (n-g)L. The difference bound needs (k-1)n < L(g-h)/2, which holds for L = 2(k-1)n+1 and g-h >= 1.
* Remark (not an error): D = (n-g)L is the *formally* reduced degree. The actual gcd can be larger (t^26 in the
  example). Dividing it out only tightens the cluster (normalisation N2), so no gcd cost is discarded.

**EXACT (`tern4.py`, `tern4.out`; my own F_9[t]).** For q=27, k=2 and the forms x, y, x+y, x-y, the orientation
(+,+,+,+) maximises g, and g = 92 (h=27). There are 249 Frobenius-9 orbits of labels (240 of size 3, 9 fixed).
After formal division every point has t-valuation 26, and stripping t^26 gives the full gcd (gcd of the reduced
points: degree 0). The four points are monic of degree 6981, with one norm of degree 13962 in F_3[t]. Pair
data, with c_ij = (z_i - z_j)/gcd(z_i,z_j):

| pair | deg(z_i - z_j) | deg gcd(z_i,z_j) | deg c_ij |
|---|---|---|---|
| (1,2) | 3471 | 2990 | 481 |
| (1,3) | 3457 | 2990 | 467 |
| (1,4) | 3457 | 2990 | 467 |
| (2,3) | 3471 | 2977 | 494 |
| (2,4) | 3471 | 3003 | 468 |
| (3,4) | 3457 | 3003 | 454 |

These are identical to Astra's six degrees. The pair gcd degrees computed by Euclid agree with the
combinatorial formula sum_p deg(pi)(e_p - |a_ip - a_jp|).

`decomp.py` checks the multiple-decomposition identity directly. The four polynomials A_H^q - c_H^{q-1}A_H,
one per kernel line H, are *the same* element of F_729[t][X], and it vanishes on V.

## 2. Port (1): the characteristic-0 polynomial analogue, and the exact failing step

**Characteristic-0 version (EXACT, sympy, `decomp.py`).** Take A = i·a with a real and c real, and s odd. Then
z = A(A^s - c^s)^2 has norm a^2(a^{2s}+c^{2s})^2, as the coordinator wrote; for s = 3 sympy returns
a^2(a^6+g^2)^2. Also z = P + 2c^{2s}A - 2c^sA^{s+1} with P = A^{2s+1} - c^{2s}A. Allow a complex parameter
g0·zeta in place of c^s. The norm then differs from the g0 norm by a product of two factors. The first vanishes
identically only for zeta = ±1. The second never vanishes identically, because of its 2a^6 term. So for a fixed
A, exactly two points share the norm: z and A(A^s + c^s)^2 = -zbar, a reflection pair.

**Lemma 1 (PROVED; the failing step).** Let k be a field, and let f, f' in k[Y] have degree q >= 2 with
char k not dividing q. Let A, A' in k[X] have degree h >= 1. If f(A) = f'(A'), then A' = alpha·A + beta.

*Proof.* Choose alpha with lc(A') = alpha·lc(A) and beta so that u := A' - alpha A - beta has no constant
term. Then u = 0 or 1 <= e := deg u < h. Suppose u ≠ 0 and put g(Y) = f'(alpha Y + beta). Taylor expansion gives

  f(A) - g(A) = f'(alpha A + beta + u) - f'(alpha A + beta) = f'_Y(alpha A + beta)·u + (terms of degree <= (q-2)h + 2e).

The X-degree of the right side is exactly (q-1)h + e, because lc(f'_Y) = q·lc(f') ≠ 0 and (q-2)h + 2e <
(q-1)h + e. So the right side is nonzero, and its degree is not a multiple of h. The left side is a polynomial
in A, so if it is nonzero its degree lies in hZ. Contradiction; hence u = 0. ∎

The lemma is where VI.4 breaks. Over F_{3^r}, f(Y) = Y^q - c^{q-1}Y has f'_Y = -c^{q-1}, a constant: the
leading coefficient q·lc(f) is 0. The four distinct right factors P_H of `decomp.py` are exactly what Lemma 1
forbids in characteristic 0. This is the polynomial-composition form of Engstrom's lattice theorem, whose
wild-ramification failure in characteristic p is classical. Equivalently, Newton's identities cannot be inverted
in characteristic 3, so the record's linearisation (round 29 `ff.md` Lemma 1.1) does not hold there.

**Consequence (PROVED).** Run the construction in characteristic 0 with common parent P = f_j(A_j),
f_j(Y) = Y^{2s+1} - C_jY:
* Lemma 1 gives A_j = alpha_j A + beta_j.
* Comparing coefficients of f_j(alpha Y + beta) = f_1(Y) gives beta = 0, alpha^{2s+1} = 1 and C_j alpha = C_1.
  Sympy lists the 2s+1 solutions.
* Conjugation compatibility (conj A_j = -A_j for A_j = i a_j) forces alpha real, so alpha = 1.

So the characteristic-0 construction produces **only the pair {z, -zbar}**. Without the special shape, any
C(t)/Q(i)(t) cluster has at most 24 points (CITED: Part IV). The step that fails over C(t) is precisely
"many decompositions of one parent with same-degree right factors". The identity (VI.5) for a single H, the
character exponents, descent, common division and specialisation all exist in characteristic 0.

## 3. Port (2): integer analogues

**3.1 Sign patterns and the odd PTE system (PROVED + EXACT).** Take z_eps = prod_{x in T}(X - i·eps_x·x).
As polynomials, deg(z_eps - z_eps') <= D/2 holds iff sum (eps_x - eps'_x)x^r = 0 for all odd r <= L =
D-1-floor(D/2): even moments cancel by the norm, and Newton applies in characteristic 0. This is the
vertical-line case of `growth/families.md` Lemma P. It is a Chebyshev system in x^2, so a nonzero difference
needs >= ceil(L/2) sign changes, and it is the mirror image of `round30` §3.2 (horizontal line, M <= 2 for
D >= 7).

EXACT exhaustive maxima (`scan.py`, squarefree, S in [1,12..16]):

| D | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|
| max points | 5 | 2 | 2 | 3 | 3 | 2 | 2 |

`pte3.py` covers the ternary-style version: exponents (1+s, 1-s), s in {-1,0,1}, with common division
(exhaustive, d <= 5, S in [1,9]). The maxima are 4, 3, 4, 5 points at reduced degrees 2, 6, 6, 6, the record's
small-D values (cf. 5 at D=6 in families.md). Multiplicities and common division buy nothing visible in
characteristic 0.

**3.2 Dirichlet-character patterns (EXACT, `intport.py`, `scan.py`).** Write Z = prod_{|x|<=B}(X - ix)^{e(x)}.
In every family e(x) + e(-x) = 2, so the common norm X·prod_{x>0}(X^2+x^2)^2 was checked exactly.
* *One odd quadratic character, multipliers lambda.* chi(lambda x) = chi(lambda)chi(x), so there are only two
  patterns (6 multipliers mod 7 give 2), and Z_- = conj(Z_+) exactly. A conjugate pair is a trivial
  (axis-touching) configuration. This is the integer shadow of "one linear form".
* *Several odd real characters* (conductors 3, 7, 11, 19, 15; B = 15, n = 31), and the *literal digit-vector port*
  (balanced base-ell digits, F_ell-linear forms, Legendre symbol; ell = 3, k = 3; ell = 7, k = 2; ell = 11, k = 2,
  m = 4 or 5, orientations chosen to maximise g). Results after the **exact Gaussian gcd**:

| family | m | template X-degree of differences (n) | template width | exact width at X=10^6 | best exact X<=3000 | best X<=2·10^5 (formal gcd) |
|---|---|---|---|---|---|---|
| chars 3,7 | 2 | 30 (31) | 0.941 | 0.948 | 0.943 | 0.953 |
| chars 3,7,11 | 3 | 30 (31) | 0.957 | 0.962-0.966 | 0.970 | 0.970 |
| chars 3,7,11,19 | 4 | 30 (31) | 0.960 | 0.965-0.971 | 0.978 | 0.975 |
| digits ell=3,k=3 | 4 | 26 (27) | 0.952 | 0.960-0.968 | — | 0.972 |
| digits ell=7,k=2 | 4 | 46 (49) | 0.929 | 0.937-0.942 | 0.946 | 0.949 |
| digits ell=11,k=2 | 5 | 118-120 (121) | 0.990 | 0.977-0.996 | — | 0.997 |

  The width exponent is log|Z'_i - Z'_j| / log|Z'|; a cluster needs 1/2 + O(1/log R). The best X are at the top
  of the range: widths fall monotonically towards the template exponent (deg diff - g)/(n - g), and stay near
  0.94 or above. Small X give about 1, as for random points. The first odd moment
  sum (chi_i - chi_j)(x)x is already nonzero (X-degree n-1), whereas in characteristic 3 the identity makes the
  moments agree up to (n-h)/2.

**3.3 Multi-parameter and number-field variants (PROVED, modulo Part IV).**
* *Roots in a totally real field.* Take x in O_K, with T and the exponent pattern Gal(K/Q)-invariant. Then the
  points lie in Z[i][X]: Galois descent replaces the F_q -> F_9 descent verbatim. They are one-parameter
  templates, so at most 24 points.
* *Several parameters.* If z_j in Z[i][X_1..X_s] satisfy deg_w(z_i - z_j) <= deg_w(z)/2 + c for weights w,
  substitute X_k = a_k u^{w_k} with integers a_k off the finitely many hypersurfaces where a weighted-top part
  vanishes. This keeps distinctness, the degrees, the norm (a real) and c, and gives a one-parameter Q(i)(u)
  cluster. So at most M(D,c) points, which is <= 24 at c = 0.
* What remains are families with no polynomial identity behind them: X depending on T, or arithmetic choices.
  These are the general problem. The arc condition becomes m-1 simultaneous conditions
  |theta_ij| <~ |Z'|^{-1/2} on the genuine angles theta_ij = 2 sum d_ij(x) arctan(x/X) (the fake-circle angle
  system with genuine angles). Nothing in VI.4 supplies these. One integer X gives one degree of freedom against
  (m-1)(n-g)/2 orders of precision (HEURISTIC count, consistent with 3.2).

## 4. Port (3): the mechanism in the record's language

**Profile (EXACT).** All primes dividing N are split: pi | Q_O and pibar | Q_{-O}, and -O ≠ O. For each label
orbit O,
* the normalised exponent is a_j = e_j(O) - min_j e_j(O),
* the norm exponent is e_p = max_j - min_j (the range formula of `astra_lemma_V` V.1),
* there are 22 row types; e_p = 2 carries 98% of W (only 286 of 13962 has e_p = 1), so **N is far from
  squarefree**.

`primes4.py` gives the F_9-prime degree histogram (degree: count) {1:6, 2:8, 3:232, 5:8, 6:240, 9:52, 15:32,
18:52, 21:52, 24:52, 30:52}, total 7982 = 728·11 - 26. In particular all three degree-2 primes of F_3[t]
(residue field F_9, of size < 2m) divide N.

**Cut identity (EXACT).** 8·sum deg(z_i - z_j) = 166272 = M(M-2)W + sum_p deg p · sum_tau (2S_ptau - M)^2 +
8·sum deg c_ij = 111696 + 31928 + 22648. Arc budget: cut + 8·sum c = 54576 <= MW = 55848, and the slack
1272 = 8·sum(D/2 - deg diff).

**Cofactors (EXACT).** For each pair, z_i = g_ij U_i and z_j = g_ij U_j with deg U = D - deg g_ij ≈ 3990.
U_i and U_j are coprime with equal norms, so U_j = (unit)·Ubar_i, and deg(U_i - U_j) = deg c_ij ≈ 470 ≈ 0.12
deg U. The integer counterpart is a cofactor at angle about |U|^{-0.88} from an axis: an extreme
"nearly real cofactor" in the sense of the record's star dictionary. The closeness is carried almost entirely by
huge pair gcds (≈ 0.43 D each), not by small chords.

**Classification.** Template type. The X-level points form a cluster over the coefficient field F_9(t) (X-degree
637, differences <= 286 < 318.5), specialised by X = t^L. This is exactly the characteristic-3 counterpart of the
record's Z[i][u] templates (`families.md` §4, `ff.md` §2), with roots on the "imaginary axis" i·F_q[t]. It is
not a covering/dominating family (Theorem D concerns characteristic 0 and families). It is not a fake: the
angles are genuine. Its integer counterpart class is bounded by Part IV.

**What holds verbatim over F_9[t]/F_3[t], and is therefore insufficient for uniformity** (they hold on an
unbounded family):
* unique factorisation and exact norms;
* the pair bound (`ff.md` Lemma 1.2, char-free): w ≈ 3990 >= L+1 = 3491 here;
* the cut identity (verified above);
* the chain lemma (combinatorial on exponent levels);
* Liouville/monomial gaps and toric forced divisibility (valuation statements);
* pigeonhole at small split primes (residue fields F_9, F_81, ... are present);
* Lemma P of Part V (char ≠ 2; VI.4 escapes it through non-constant cofactor differences and exponents {0,2});
* Schwartz-Zippel and intersection theory short of BMY.

The (2/3+eps) argument should transfer as it did to F_2[t] (Part V reconcile). I did **not** redo that
transfer for F_3; either way VI.4 has W = exp(O(m log m)), astronomically above any c·m log m.
Cilleruelo-Cordoba: the example is consistent with its conclusion (relative deficit (g-h)/(2(n-g)) ≈ 2^{-m}).
I did not check whether its proof transfers.

**What does not hold over F_3[t]**, and so survives as possible input:
* Newton linearisation (`ff.md` Lemma 1.1: equal e_1..e_L <=> equal power sums, Q-linear in exponents);
* the balance half of the saturated balance lemma (real place, argument principle);
* Theorem D (generic smoothness);
* Part IV / Langer's BMY;
* Lemma 1 above (Engstrom/Ritt).

**Filter (recommendation).** **F3+ extends to p = 3 with the reflection form A^2+B^2**: no argument valid
verbatim for F_9(t)/F_3(t) (deg for log, the inert degree-one infinite place in place of the real place)
proves the uniform bound. This removes the README_V §3(d) exemption "arguments valid only in odd
characteristic or only for Kummer-type circles".

The exemption for arguments that *essentially* use squarefree N survives over F_3. Over F_2 it is challenged
by Astra's VI.5', which is not refereed in this block. Growth-rate arguments are untouched.

## 5. Port (4): an honest answer

**Is VI.4 evidence against the uniform conjecture over Z? No.**

* The feature F_3[t] has and Z lacks, used essentially: **nonlinear additive polynomials, hence one parent with
  many same-degree decompositions**. In the record's linear language: Newton's identities cannot be inverted,
  so closeness is not a set of Q-linear moment conditions on exponents. In characteristic 0 Lemma 1 forbids
  it, and over Z there is no additive structure at all.
* **Not essential:**
  * the constant-field extension and descent (they port to number fields);
  * the ultrametric (integer specialisation of templates works archimedeanly with bounded C, `bridge.md`);
  * common division by G (archimedean division by a common Gaussian factor improves the width exactly as
    in degree terms).
* **Archimedean/char-0 obstruction visible:** yes, at the level of mechanisms. Every characteristic-0 tool that
  bounds function-field clusters (linearisation, the balance lemma's real-place half, BMY) fails in
  characteristic 3, and VI.4 lives in that gap. Over Z the arc condition is a single real-linear condition on
  exponents weighted by genuine, non-torsion angles. That condition is char-0-like, but it is also what the
  record's fakes exploit, so it is no proof route by itself.
* VI.4 also says nothing at growth scale.

**Significance.** It is model-level, but it is the most Z[i]-like failure yet: odd characteristic, the form
A^2+B^2, conjugation a reflection, inert infinity. It closes README_V §5 item 1: sup_D M_3(D,0) = infinity
for non-squarefree N, against the exhaustive M_3(D,0) <= 4 for D <= 11, which shows again that small-D data
are no guide. It sharpens the filter: a uniform proof must use characteristic 0 or archimedean input in an
essential way.

**Novelty.** The record has no odd-characteristic anisotropic construction; README_IV called the
p-linearised analogue open. The ingredients are classical: subspace/linearised polynomials (Ben-Sasson,
Kopparty, Radhakrishnan, already cited for Part V), quadratic-character multiplicities, Engstrom/Ritt. The
assembly appears new.

**Correction to Astra's "transfer" paragraph.** She lists coefficient-Frobenius descent among the added
inputs. Descent is not specific to characteristic p (it is Galois descent). The single non-portable input is
the multiple decomposition of P_V, i.e. Lemma 1's failure. Her conclusion, "neither an integer growth theorem
nor an integer construction", is correct.

## 6. Files

All in `ref6_port/`:
* `tern4.py`: four-mark reproduction, pair data, cut identity, profile.
* `primes4.py`: prime degrees.
* `decomp.py`: the characteristic-3 multiple-decomposition identity; characteristic-0 norm and symmetry checks.
* `intport.py`: integer ports E1-E3.
* `scan.py`: exact X-scan; vertical-line odd-PTE exhaustive search.
* `pte3.py`: ternary-exponent vertical templates.
* `*.out`: outputs.

All runs are under prlimit 3.5 GB and single-process.
