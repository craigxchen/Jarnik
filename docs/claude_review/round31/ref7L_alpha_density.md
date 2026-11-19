# Referee report (lens: density and quantifiers): does the polynomial Pell five-mark family give alpha^ess_e(-K) <= 4?

Verdict: **CORRECT_WITH_CORRECTIONS.** The bound alpha^ess_e(-K) <= 4 over Q holds, and so does "A_5(delta) is false for
every delta < 0". The argument the judge wrote down rests on the wrong density statement ("the tuples are Zariski
dense as seeds and T vary", "already their values at t = 1 are dense"). That statement does not imply the
conclusion. The correct argument goes through the seed-by-seed curves C_s and is given in Section 4. Five smaller
corrections are listed in Section 7.

Scripts (all mine, exact arithmetic, nothing imported from the referees):
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref7L_alpha_density/`
(`fam.py`, `check1.py`, `jac.py`, `spec.py`, `dbg.py`, `gcdbound.py`, `bezout.py`).

## 1. Definitions used (from the record)

* T_5 = (R^1_{Q(i)/Q} G_m)^5 / diag. It is a 4-dimensional Q-torus and the open torus of the smooth weak Fano 4-fold
  X = Lbar^tw_5 (expert_note.md §3, "Setting"; l. 336).
* A primitive k-tuple is a tuple of distinct lattice points of |zeta| = R'' with Gaussian gcd 1. Its class is
  [zeta] = (zeta_j/zeta_1)_j in T_k(Q). The class is invariant under zeta -> gamma*zeta for gamma in Q(i)^*.
* H = H_{-K}([zeta]) = R''^2, where R'' is the radius of the **reduced** (primitive) tuple (gauss.md Lemma 1.1;
  round22/lit.md Lemma 1.1). dist(y, e) is any fixed archimedean distance on X(R). Near e it is comparable to
  max_j |zeta_j/zeta_1 - 1|, which is comparable to diam(zeta)/R''.
* A_k(delta) (expert_note.md l. 238): there exist a proper Zariski-closed Z in T_k and an R_0 such that every primitive
  k-tuple of distinct lattice points of a circle of radius R'' >= R_0 with all mutual distances <= R''^{1/2-delta}
  has [zeta] in Z. The boxed definition places no restriction on delta. The reduction theorem (gauss Thm 3.4')
  uses 0 < delta < 1/(2^k-2).
* alpha^ess. The record has two definitions.
  * (D1) round22/toric.md §1.2, following McKinnon–Roth: for a sequence of distinct y_i in X(Q) with
    dist(y_i, e) -> 0, alpha({y_i}) = inf{gamma >= 0 : dist(y_i,e)^gamma H(y_i) is bounded};
    alpha_e(Omega) = inf over such sequences in Omega; alpha^ess_e = sup over proper closed Z of
    alpha_e(Omega \ Z). Pagelot/Huang Def 4.7 is the same with sup over dense open U.
  * (D2) papers_check/README.md l. 65, Lehmann–McKinnon–Satriano Def 1.1: alpha_ess = inf of alpha({x_i}) over
    Zariski-dense sequences approximating e.
  * The record notes (D1) = (D2) (papers_check item 55). Both are treated below.
* The dictionary (round22/toric.md Lemma 1.3): an A-point at delta has exponent gamma = 4/(1+2 delta).

## 2. The family (reproduced)

A(t) = [[t^2+1, t],[t, 1]], beta = t + 2i, G_j = X_j + iY_j with (X_j, Y_j)^T = A(t)^n (1, s_j)^T. The marks are
z_0 = beta * prod_l conj(G_l) and z_jk = conj(beta) G_j G_k prod_{l != j,k} conj(G_l) for jk in {12, 13, 23, 14}.
"conj" is coefficientwise conjugation, with t and s real.

`check1.py` uses three seed vectors, (1/2,-2,7/3,4), (-1/3,5/2,6,-7) and (0,1,3/5,-11/4), with n = 1, 2. In all six
cases the five marks are monic of degree D = 8n+1 and have one norm z*conj(z) in Q[t]. All ten pairwise differences
have degree exactly 4n+1. The gcd over Q(i)[t] is 1.

## 3. Dominance (attack point 1): exact Jacobian

Write G_j = c + s_j d, where (c, d) = (1, i)A^n. Put u_j = G_j/conj(G_j) and b = conj(beta)/beta. Then
q_jk := z_jk/z_0 = b u_j u_k. Let a_j = d log u_j / d s_j. Then:

* d conj(c) - c conj(d) = 2i det(A^n) = 2i. So a_j = 2i/(G_j conj G_j), which is nonzero at every real t.
* The 4x4 block of the log-Jacobian of (q_12, q_13, q_23, q_14) with respect to (s_1, ..., s_4) is
  [[a1,a2,0,0],[a1,0,a3,0],[0,a2,a3,0],[a1,0,0,a4]]. Its determinant is -2 a1 a2 a3 a4 = -32 / prod_j |G_j|^2, which
  is nonzero.

`jac.py` checks this symbolically. For n = 1, 2 and four rational points (s, t), it computes the full 4x5 Jacobian
exactly from the rational functions q_jk(s, t). In all eight cases the rank is 4, and the seed-block determinant
equals -2 prod (2i/|G_j|^2) exactly. For example, at s = (1/2,-2,7/3,4), t = 1, n = 1 the determinant is
-576/278953.

Consequence. For every fixed n and every fixed real t, the map s -> [z(s,t)] is dominant (indeed étale) from the
seed space onto T_5. In particular the map (s, t) -> [z] is dominant. This settles point 1 positively, without
floating point.

## 4. The quantifier argument (attack points 3 and 4)

**Logical form.** By (D1), alpha^ess_e(-K) <= gamma_0 means: for every proper Zariski-closed Z in X, there is a
sequence of distinct points y_i in T_5(Q) \ Z with dist(y_i,e) -> 0 and dist(y_i,e)^{gamma_0} H(y_i) bounded. The
bound may depend on Z. This is a sup over Z of an inf over sequences, so an upper bound needs one good sequence for
**each** Z. That is the direction proved here. It suffices to do this for gamma_0 = 4 + 1/(2n) and every n; then
alpha^ess <= inf_n (4 + 1/(2n)) = 4.

**Why "the tuples are Zariski dense" is not enough.** A dense set of points with bounded T does not converge to e.
For example, the values at t = 1 are dense (Section 3), but they give no information about alpha^ess. What is
needed is density of points with T -> infinity **along one curve at a time**.

**Correct argument.** Fix n and a proper closed Z (over Q, or even over C). Let phi(s, t) = [z(s, t)]. This is a
rational map A^4 x A^1 -> T_5, defined on a dense open set. For s with pairwise distinct entries, let C_s be the
closure of t -> phi(s, t).

1. The set S_Z = {s : C_s ⊂ Z} equals the intersection over t of {s : F(phi(s,t)) = 0 for all F in I(Z)}. It is
   therefore Zariski-closed in A^4. It is proper, because otherwise the image of phi would lie in Z, contradicting
   dominance (Section 3).
2. Q^4 is Zariski dense in A^4, as is the set of rational seed vectors with distinct entries. Hence some rational
   seed vector s has distinct entries and satisfies C_s ⊄ Z.
3. C_s is an irreducible curve (the image of A^1). It is not a point: the marks are distinct polynomials, so
   phi(s, t) is not constantly e, while phi(s, t) -> e as t -> infinity. Therefore C_s ∩ Z is finite.
4. For all but finitely many integers T, the point y_T = phi(s, T) is outside Z, the five integer marks are distinct,
   and the y_T are pairwise distinct (the map t -> phi(s,t) is finite-to-one). So there are infinitely many integer
   T with y_T in T_5(Q) \ (Z ∪ Z_diag).
5. Exponent for this fixed s. Let L be the least common multiple of the seed denominators. The cleared tuple z'(T)
   lies in Z[i]^5, and its radius R' satisfies R' ≍_s T^D. Every pair difference satisfies
   |z'_a - z'_b| ≍_s T^{4n+1} (degree exactly 4n+1). Hence dist(y_T, e) ≍_s T^{-4n}. Let g(T) be the Gaussian gcd.
   Then H(y_T) = R'^2/|g(T)|^2 <= R'^2 ≪_s T^{2D}. Therefore
   dist^{gamma_0} H ≪_s T^{-4n gamma_0 + 16n + 2}, which is O_s(1) at gamma_0 = (16n+2)/(4n) = 4 + 1/(2n).
6. Every constant here depends on s, and s depends on Z and n. (D1) allows this, since the bound is per sequence.
   No uniformity in the seeds is used or needed.

**Form (D2).** Choose a sequence s^(1), s^(2), ... of rational seed vectors whose set is Zariski dense in A^4
(for example, all rational seed vectors with distinct entries). At stage k, for each m <= k, choose a new integer
T_{m,k} so that dist(y, e) < 1/k and dist^{gamma_0} H <= C_{s^(m)}. The bound C_{s^(m)} is not uniform in m, so it
is normalised away by taking gamma slightly larger: use gamma_0 + epsilon and require
dist^{gamma_0+epsilon} H <= 1. This is possible because the left side tends to 0 as T -> infinity for fixed m.
* The resulting sequence converges to e.
* It has dist^{gamma_0+epsilon} H <= 1.
* It is Zariski dense. Each curve C_{s^(m)} contains infinitely many terms, so the closure contains every
  C_{s^(m)}. Their union is dense, because the s^(m) are dense and phi is dominant.

So alpha_ess(D2) <= 4 + 1/(2n) + epsilon for every n and every epsilon > 0, hence <= 4.

**Numerical confirmation (`spec.py`, exact Z[i] arithmetic, reduced tuples, H = R''^2 after dividing by the gcd).**
For n = 1, 2, 3, the three seed vectors, and T up to 10^15, the quantity log(dist^{gamma_0} H) stays bounded:

| n | gamma_0 | seed (1/2,-2,7/3,4) | seed (-1/3,5/2,6,-7) | seed (0,1,3/5,-11/4) |
|---|---|---|---|---|
| 1 | 4.5 | 9.8 to 11.4 | 13.3 to 16.4 | 8.1 to 11.5 |
| 2 | 4.25 | 8.0 to 10.8 | 10.3 to 15.4 | 7.6 to 11.0 |
| 3 | 4.1667 | 7.8 to 10.7 | 12.4 to 15.1 | 7.5 to 10.9 |

The effective exponent -log H / log dist decreases towards gamma_0. At T = 10^15 it is 4.58, 4.29 and 4.19 for
n = 1, 2, 3. The approach is slow because the constants depend on s.

## 5. Denominators and gcd (attack point 2)

* Clearing denominators means replacing (1, s_j) by (b_j, a_j), where s_j = a_j/b_j. This multiplies G_j by b_j, and
  so multiplies every mark by the same factor prod_l b_l: each mark contains each G_l or conj G_l exactly once. The
  cleared tuple lies in Z[i]^5, its five entries have one norm, and its class [z] is unchanged.
* Dividing by the Gaussian gcd g(T) does not change the class. It lowers H from R'^2 to R'^2/|g|^2 and leaves dist
  unchanged. For an **upper** bound on alpha this works in our favour, so even an unbounded gcd would be harmless.
* The gcd is in fact bounded for fixed seeds, not uniformly. `bezout.py` computes, for n = 1, a Bezout identity
  sum a_i z'_i = 1 over Q(i)[t]. Its common denominator c is a rational integer: for seeds (1/2,-2,7/3,4),
  c = 2^6·3·5^3·13·17^2·29·193·5009·70709. Hence g(T) divides c for every integer T, and N(g(T)) divides c^2. This
  holds for all |T| <= 300, where the observed maximum N(g) is 24650 (38480 for the second seed vector). The bound
  grows with the seeds. The judge's "observed values {160, 2080, 4160}" refers to a particular sample, not a bound.
* Distinctness of the marks holds only for all but finitely many T. For example, z_0 = z_14 at T = 34 for seeds
  (1/2,-2,7/3,4), for every n = 1, 2, 3, and z_0 = z_12 at T = 2 for seeds (0,1,3/5,-11/4) (`dbg.py`). The claim
  "five Gaussian integers of one norm" therefore needs "for all but finitely many T".

## 6. The meaning of "A_5(delta) false for every delta < 0" (attack point 5)

* The boxed A_k(delta) puts no restriction on delta. For delta < 0 the condition "mutual distances
  <= R''^{1/2-delta}" is weaker, so A_5(delta) is a **stronger** statement than A_5(delta') for delta' > delta. The
  negative range is meaningful. It is trivially false for delta <= -1/2, since T_5(Q) is dense (T_5 is Q-rational)
  and Northcott applies.
* Correct wording. *For every delta < 0, every proper Zariski-closed Z ⊂ T_5 and every R_0, there is a primitive
  5-tuple of distinct lattice points on a circle of radius R'' >= R_0, with all mutual distances <= R''^{1/2-delta}
  and [zeta] not in Z.*
* Proof. Choose n with 1/(2(8n+1)) < -delta, and take the sequence of Section 4. After reduction, the distances are
  ≪_s R''^{(4n+1)/(8n+1)} = R''^{1/2+1/(2(8n+1))}, because g is bounded. The constant is absorbed by the strict
  exponent gap once R'' is large.
* The reduction theorem gauss Thm 3.4' (0 < delta < 1/30) is not affected. A_5(0), that is, arcs of length R''^{1/2},
  is not decided, since every n gives an exponent strictly above 1/2.
* Equivalence. By round22/toric.md Lemma 1.3(a)/(b), "A_5(delta) false for all delta < 0" is equivalent to
  alpha^ess_e(-K) <= 4. The judge's "i.e." is therefore correct.

## 7. Corrections to the claim as stated

1. The density sentence must be replaced by the curve-by-curve argument of Section 4. "Zariski dense in T_5 as
   seeds and T vary", and the judge's "already their values at t = 1 are", are true but irrelevant: points with
   bounded T do not approach e.
2. All constants (the gcd bound, the implied constants in R ≍ T^D and |z_a - z_b| ≍ T^{4n+1}, and the finite
   exceptional T) depend on the seeds. Uniformity is neither available nor needed. The final quantifier is
   "for every n, for every Z, there exist s, and then T -> infinity".
3. Distinctness and "outside Z" hold for all but finitely many T, not for all T.
4. The exponent for fixed seeds is <= 4 + 1/(2n) with no o(1), because g is bounded.
5. The bound concerns the unrestricted essential constant. The profile has 1|4 cuts, so it is unbalanced, and nothing
   follows for the balanced constant or for (L_5°).
   * The known window is 10/3 <= alpha^ess_e(-K) <= 4. The lower end is McKinnon–Roth base-locus Liouville, as in
     the record.
   * Question B, in its unrestricted form, now asks whether alpha^ess_e(-K) lies in (15/4, 4].
   * The literature form of Question B (Lehmann–McKinnon–Satriano essential boundedness, alpha^ess >= 4) becomes
     equivalent to alpha^ess_e(-K) = 4 exactly.

**Geometric reading.** For fixed s, C_s is a rational Q-curve with -K·C_s = 2D = 16n+2 (from H ≍ T^{2D} with bounded
gcd). It has a branch at t = infinity of multiplicity 4n at e, since all ten pair contacts equal 4n. Its integer
points realise McKinnon's curve constant d/m = (16n+2)/(4n). So the upper bound is the standard "covering family of
rational curves with -K·C/mult_e C -> 4" bound. What is new is the existence of the family. Very free smooth curves
give only -K·C >= 5. I found no statement of alpha^ess_e(-K) <= 4 in the record before round 31: expert_note.md
l. 338 calls the window (15/4, 4) "untouched", and README_VIIL.md l. 213 says the literature was not checked. I did
not check the external literature.

## 8. Statement I sign

For every integer n >= 1 and every proper Zariski-closed Z ⊂ X = Lbar^tw_5 (or Z ⊂ T_5, over Q or over C), there are a
rational seed vector s, with pairwise distinct entries and depending on Z and n, and a constant C(s) such that:
* for all but finitely many integers T, the reduced five-tuple obtained from z(s, T) is a primitive five-tuple of
  distinct Gaussian integers of one norm, with class y_T in T_5(Q) \ Z;
* y_T -> e as T -> infinity;
* dist(y_T, e)^{4+1/(2n)} H_{-K}(y_T) <= C(s).

Consequently alpha^ess_e(-K) <= 4 over Q, in both forms (D1) and (D2), and A_5(delta) is false for every delta < 0.
Together with the record's lower bound, 10/3 <= alpha^ess_e(-K) <= 4. Nothing follows for balanced sequences,
for (L_5°), for A_5(delta) with delta >= 0, or for whether the exponent 4 is attained on a Zariski-dense set.
