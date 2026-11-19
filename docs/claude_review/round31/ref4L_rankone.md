# Referee/analyst report, round 31, Part IV, block "rankone" (the residual branch of IV.4)

Scripts are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref4L_rankone/`. I wrote all of them and ran none of Astra's.

| Script | What it does |
|---|---|
| `iv2_rankone.py` | Exact reduction of the IV.2 defect δ_i modulo ∏(X−x_j), with symbolic s1, s2, s3, s5, P, q |
| `iv2_exactform.py`, `iv2_disc.py` | IV.2 models re-read as exact clusters for their interpolating form Q_*. Computes disc(Q_*) and checks splitting of the block primes |
| `split_cycle.py` | Exact Gröbner basis for constant ℓ, split Q_* (polynomial families) |
| `sigma_sys.py`, `sigma_root*.py`, `sigma_complex*.py`, `tau_solve.py` | The σ-system of §2.3: real and complex Newton / least-squares searches |
| `reduction_check.py` | Validates the reduction (G1) → σ-system at 40 digits |
| `search_s3y.py`, `count_s3y.py`, `real_s3.py` | Integer and real search on the IV.2 rank-one surface |
| `z_sqrt2_example.py` | An exact indefinite rank-one model: genuine Z[√2] stars |
| `kernel.py` | Integer kernel of the K5 incidence matrix |

Notation is as in IV.4. h ∈ {4,5} is the orientation, D_δ = max|δ_i|, E = I − Q_*, and A = |ℓ|². Choose an integral basis e1, e2 with ℓ(e1) = 0 and ℓ(e2) = 1. Put B = e1·e2 and C = |e2|², so AC − B² = 1. Write V_i = x_i e1 + y_i e2, so that y_i = ℓ(V_i) ≠ 0.

## Verdicts

| Item | Verdict |
|---|---|
| (a) Arithmetic of Q_* = I − λℓ²: discriminant, definite / indefinite / split / degenerate, size of disc, ideal-theoretic dictionary | CORRECT WITH CORRECTIONS. The precise statements are Theorem R1 below. Astra did not state them; this part is analysis, not refereeing. |
| (b1) Liouville 10/3 for an arbitrary binary form | CORRECT, and uniform in the reduced (discriminant-normalised) metric. Theorem R2. |
| (b2) Exact five-point clusters at reduced exponent 10/3 for forms of growing discriminant | EXIST. Every Part IV.2 model is one: Q_* is positive definite with \|disc\| ≍ P³ ≍ N^{3/10}, and the shapes escape any fixed tangent cone. Theorem R3, exact. A uniform-in-the-form Question B is therefore false. Question B for a fixed form is untouched. |
| (b3) The task's explicit question: can IV.2's leading residual be made rank one and the lower terms killed? | Leading term rank one ⇔ s3 = −s1s2. Then the whole defect is exactly rank one, δ_i = (q³s5/4)(4Ps2 − 4qs1s2² + 5qs5)x_i². Going below P needs s2 = 0, which is impossible (Rolle). With s2 ≠ 0 the defect stays ≍ P. Exact. |
| (b4) Any one-parameter polynomial pencil (h = 4, linear blocks, constant cofactors) with defect o(P) | Constant-ℓ cases: PROVED impossible. Linear-ℓ case: reduced exactly to the "σ-system" (§2.3). Every real and complex solution found is degenerate. So the quarter-power barrier of IV.2 is optimal for this class, modulo one numerical step (EVIDENCE, not proof). |
| (c) Literature | Uniformity in the form for k ≥ 3 points: nothing found. Two-point bounds with explicit d-dependence exist (snippet level only). See §3. |
| (d) Anything new for genuine circle tuples | Nothing, for a structural reason. There is one new barrier (F-form): hypotheses that use only "exact values of *some* rational quadratic form of bounded denominator" cannot beat 10/3. Any proof must use that the form is the Euclidean one, or at least of bounded discriminant. |

## 1. (a) The residual form, made precise

**Theorem R1 (proved here; elementary).** Assume the rank-one branch of IV.4: δ_i = λℓ(V_i)², with λ ∈ H_t^{−1}Z∖{0} and ℓ primitive. Put u_i = Ax_i + By_i and k = Aλ − 1. Then:

1. **Normal form.** A·n_i = u_i² − k·y_i² and A·|V_i|² = u_i² + y_i². So Q_*(x,y) = ((Ax+By)² − ky²)/A, det Q_* = 1 − λA = −k, and disc Q_* = 4k. (Checks: |V|² = Ax² + 2Bxy + Cy² = ((Ax+By)² + y²)/A; and Q_* = |V|² − λy².)
2. **Signature.**
   * λ < 0 ⇒ k < −1 ⇒ Q_* is positive definite (imaginary quadratic).
   * 0 < λ < 1/A ⇒ definite as well. This is possible only when A < H_t.
   * λ = 1/A ⇒ k = 0 ⇒ Q_* = (Ax+By)²/A is degenerate. Then A·n_i is a square for all i, so n_in_j is a square. On private profiles the n_i are distinct squarefree numbers, so **the degenerate case is impossible**.
   * λ > 1/A ⇒ indefinite: real quadratic if k is not a square, split if k = w² (then A·n_i = (u_i − wy_i)(u_i + wy_i)).

   The cone condition and n_i > 0 force nothing here. Q_*(e1) = A > 0, and the V_i cluster around e1 because |λ|y_i² ≤ D_δ ≪ n_i. All three signatures occur algebraically:
   * Z[√−2] stars (λ = −1), as in ref4L_dichotomy;
   * stretched stars (λ = 1 − 1/d², definite);
   * genuine Z[√2] stars (Q_* = X² − 2Y² = I − 3Y²; exact, all ten divisibilities, primitive; `z_sqrt2_example.py`).

   None of these has small defect.
3. **Size of the discriminant.** For every pair, |m_ij| = |u_iy_j − u_jy_i|/A ≤ 2·sqrt((n_max + D_δ)·D_δ/(A|λ|)). Hence

   |k| ≤ 1 + A|λ| ≤ 1 + 4(n_max + D_δ)·D_δ / max_{ij} m_ij².

   Below the IV.4 threshold (D_δ ≤ P^{h−3−12κ−ε}, m_ij ≥ g_ij = P^{h−3}) this gives |disc Q_*| ≤ P^{3−ε+O(κ)} = N^{3/10−ε'} in both orientations.
4. **Ideal-theoretic dictionary.**
   * Let λ ∈ Z and put 𝔞 = [A, B+√k] ⊂ O = Z[√k]. Then N(B+√k) = B² − k = A(C − λ), so 𝔞 is an O-ideal of norm A and Q_*(x,y) = N(x·A + y(B+√k))/A. That is, Q_* is the norm form of 𝔞.
   * For λ ∉ Z, scale by den(λ) | H_t; the same holds with a bounded fixed factor.
   * The five columns are ξ_i = u_i + y_i√k ∈ 𝔞, with 𝔟_i := (ξ_i)𝔞^{−1} integral of norm n_i. All five 𝔟_i lie in the single class [𝔞]^{−1} of Pic(O). Nontrivial class group is allowed and costs nothing.
   * A block prime p ∤ 2k·den cannot be inert: p | n_i gives u_i² ≡ ky_i² with p ∤ y_i (|y_i| < p).
   * If p | k, then v_p(k) = 1 and p ramifies. Such a block contributes nothing to the toric height of the K-point.
   * At split blocks the Gaussian orientation argument of ref3_model §1 applies verbatim, so the same prime 𝔭_e divides ξ_i for every i marking e.
   * Consequence. With K = Q(√k), the class y = (ξ_i/ξ̄_i)_i is a point of the twisted torus T_5^K = (R^1_{K/Q}G_m)^5/diag, with the same 15-cut profile (minus ramified blocks) and height N·P^{O(κ)}. **Below the threshold, an integral model is exactly a genuine K-cluster for O = Z[√k] with |k| ≤ N^{3/10−ε'}.** The order may be non-maximal, may have nontrivial class group, and may be real quadratic. In the real case the units of O act by automorphs: they preserve the n_i and all minors, but not the Euclidean defect.

   Correction to the framing of the task: the toric fourfold depends only on K = Q(√k), not on the order or the ideal class. τ_e = 10/3 for every quadratic twist (round 23). Vojta's prediction is the same for each K, but with constants Z_K, R_0(K) that depend on K.

## 2. (b) Liouville, uniformity, and the constructions

**Theorem R2 (Liouville for any form; proved).** Let F be primitive integral of discriminant D ≠ □, with F(V_i) = n_i (up to a common fixed factor) and g_ij | det(V_i,V_j) ≠ 0. The F-angle (sinh of the hyperbolic angle if D > 0) satisfies

sin θ_ij = √|D|·|det(V_i,V_j)| / (2√(F(V_i)F(V_j))) ≥ (√|D|/2)·g_ij/√(n_in_j).

So in T_5^K, dist(y,e) ≥ c√|D|·N^{−3/10} on balanced ten-pair profiles. Liouville 10/3 holds uniformly in the **reduced** metric dist/√|D|, and the absolute bound improves with |D|.

**Theorem R3 (exact clusters at reduced exponent 10/3 for growing discriminant; exact).** Take Astra's IV.2 models. The two interpolation relations hold identically, because the fourth divided difference of the cubic h vanishes. So there is a rational Q_* with Q_*(V_i) = n_i exactly. The checks below are in `iv2_exactform.py` and `iv2_disc.py`, symbolic in P and q:

* Q_* = aX² + bXY + cY² with a = (719P − 3562392150q)/(141717600q), b and c of degrees 2 and 3 in P/q.
* disc(Q_*) = (−719P³ + 3704109750P²q − …)/(3.347·10^17·q³), a **cubic in P/q**. It is negative for large P, so Q_* is positive definite with |disc| ≍ P³ ≍ N^{3/10}.
* den(Q_*) is fixed for fixed q, and divides H_t.
* All block primes are split or ramified, with no inert ones (35 prime factors checked at P ≈ 10^15).
* det E = −(4112/243)² (consistent with ref4L_dichotomy).

Hence these are genuine exact clusters for the imaginary quadratic orders of the forms F_P = den·Q_*, with reduced width ≍ N^{−3/10} (Liouville-sharp). By IV.2's own leading-vector argument, [x_i − x_5] can be chosen off any prescribed homogeneous G = 0. **Consequence: no version of A_5(δ), δ ≤ 1/10, with exceptional data uniform over all imaginary quadratic orders in the reduced metric can hold.** In absolute terms these points sit at dist ≍ N^{−3/20} from e, so they contradict nothing for a fixed K, and they contradict nothing in the record. In function-field terms, the IV.2 exact family lives on the genus-1 curve w² = disc(P). The rank-one families below would live on a conic.

**(b3) Killing IV.2's lower terms (exact, `iv2_rankone.py`).** In IV.2 the only primitive ℓ with ℓ(V_i) small is ℓ ∝ (L, −P), with ℓ(V_i) = LB·x_i. So rank one means δ_i ∝ x_i². Reducing δ_i modulo ∏(X − x_j), with s4 = 0, b0 = 3s5/2, L² = −s1s5q² and B = b0q², gives

δ_i = −Pq³s5(s1s2+s3)·x_i + (q³s5/4)(4Ps2 − 4qs1s2² + 5qs5)·x_i² + q⁴s1s5(s1s2+s3)·x_i³ − q⁴s5(s1s2+s3)·x_i⁴.

* So δ is rank one **iff s3 = −s1s2**. The defect is then exactly (q³s5/4)(4Ps2 − 4qs1s2² + 5qs5)·x_i², at every order.
* Defect o(P) would need s2 = 0. Then the marks are roots of X⁵ − s1X⁴ − s5, whose derivative X³(5X − 4s1) has only two distinct real roots. So there are at most 3 real roots, a contradiction.
* With s2 ≠ 0 the defect stays ≍ Pq³, so there is no gain. Such a family would at least be exact at the threshold with |k| ≍ P³ (the top of R1.3).
* Even this needs integral points on the degree-7 surface e1(y) = 0, e2(y)e5(y) + e3(y)e4(y) = 0 (y_i = 1/x_i), with −e4(y) a square.
* The surface has plenty of real points: 19733 of 35601 sampled points have −e4 > 0. There are **no** integer points with distinct nonzero |y_i| ≤ 30 (`count_s3y.py`). The surface is of general type, so this is an arithmetic obstruction, not a proof of emptiness.

**(b4) All one-parameter polynomial pencils with defect o(P) (h = 4, linear blocks J_e(P), constant cofactors).** The defect is then constant. By IV.4 (κ = 0, applied for all large P) E(P) has rank one over Q(P), so E = λℓ² with ℓ ∈ Q[P]² primitive. The δ_i are constant, and degenerate J's are excluded, so ℓ(V_i) is constant. Classification by deg ℓ:

1. **deg ℓ = 0, Q_* non-split.** The linear J_e are inert in the constant-field extension, so J_e | N(ξ_i) forces J_e | y_i, a constant. Impossible (proved).
2. **deg ℓ = 0, Q_* split.**
   * Orientation consistency makes the "+" edges a 5-cycle.
   * The conditions become ∏_{C+ at i}(P−ρ) − ∏_{C− at i}(P−ρ) = const, with the same constant for all i.
   * Exact Gröbner basis (`split_cycle.py`): only solutions with at most two distinct ρ. Impossible (proved).
3. **deg ℓ = 1 (the general case, which contains IV.2).** After a real rotation, V_i = (x_iP + y_i, L·x_i) with x_i linear and y_i constant. N(ξ_i) = A·n_i holds over the conic C: w² = μP² + (μ−1)L² (verified symbolically). Parametrising C and counting zeros (deg 6 = 2 over A + 4 over the J_e) turns the problem into the **σ-system**: ten edge parameters σ_e with, at every vertex i,

   ∏_{e∋i} σ_e = K  and  Σ_{e∋i} (σ_e + c/σ_e) = 0,

   the latter being Σ_{e∋i} w(σ_e) = 0.
   * The reduction was validated at 40 digits: at solutions y_i is constant and x_i is linear, and at random non-solutions y_i is not constant (`reduction_check.py`).
   * Expected dimension: 1. Real search: 64 Eulerian sign patterns × 60 starts gave 720 solutions, all degenerate (|σ_e| constant).
   * Complex search in the kernel-torus coordinates (τ = ω·z^B, 5 unknowns): about 200 distinct solutions for each of three β values. **Every one has coincident ρ's**.
   * Independent cross-check without the reduction: a direct least-squares solve of all 25 coefficient equations of (G1) (`g1_quick.py`, output in `g1_quick.out`). It found 28 zero-residual solutions. All 28 are degenerate: leading minors ≈ 0, coincident J-roots, μ = 1.
   * Verdict: EVIDENCE that no such family exists. This is not a proof (no certified solve; Gröbner was out of reach for sympy).
4. **deg ℓ ≥ 2.** The x_i are constant, so all J are equal. Impossible.

So, modulo the numerical step in item 3, **IV.2's quarter power is optimal among one-parameter polynomial pencils at h = 4**. For h = 5 there are no polynomial families with linear blocks, by parity of degrees: IV.3 uses |α|² = p_* linear. That class is not analysed here.

No non-polynomial (Pell/unit-driven) construction was attempted beyond the observation in R1.4 that automorphs preserve exactness but not the Euclidean defect. A random-model count (≈ P^{−12}D_δ⁵ per dyadic block) predicts no small-defect models beyond structured families.

## 3. (c) Literature

Fetches of arXiv, core.ac.uk and geodesic.mathdoc.fr failed (hosts not allowed), so everything here is from search snippets. Within that limit I found:

* **Cilleruelo–Córdoba.** ≤ l points on arcs of length √2·R^{1/2 − 1/(4⌊l/2⌋+2)}. Circles only.
* **Cilleruelo–Granville, *Close lattice points on circles*.** Classification of four close points; circles only.
* **Ramana, Acta Arith. 130 (2007), arXiv math/0703286.** At most 2 integer points on arcs of length R^{1/6}/|d|^{1/3} of X² + dY² = R. The d-dependence makes shorter arcs for larger d, consistent with R2.
* **Cilleruelo–Garaev–Ostafe–Shparlinski, arXiv 1012.3569 (attribution from memory).** Cites a lemma: on x² − Dy² = n, D squarefree, an arc of length n^{1/6} has at most 2 points (a case of an earlier theorem, likely Cilleruelo–Jiménez-Urroz; NOT verified).

I found no uniform-in-the-form statement and no counterexample for k ≥ 3. Theorem R3 appears to be the first explicit non-uniformity, in the reduced metric. The record has no uniformity statement either (grep for "uniform in the form", "growing discriminant", "Q_*": no hits outside unrelated files). The split analogue (round 23/27) is the D = □ case; its known curves sit at 15/4 and 7/2, not 10/3.

## 4. (d) Consequences

* **For genuine circle tuples (δ = 0): nothing new.** IV.4 at E = 0 is vacuous. R1 to R3 are statements about forms other than x² + y², and nothing transfers back, because Euclidean exactness is the hypothesis rather than a conclusion.
* **New barrier (F-form).** Consider an argument whose inputs are:
  * a primitive 2×5 integer matrix;
  * ten minor divisibilities with cofactors P^{o(1)};
  * index one;
  * exact values Q(V_i) = n_i for **some** positive-definite rational binary form Q with den(Q) = P^{o(1)}.

  Such an argument cannot beat 10/3 at five points off any fixed tangent cone, since the IV.2 models satisfy all of these (R3). So the exactness of norms is worthless unless the argument uses the **height/discriminant of the form**: either Q = x² + y², or at least |disc Q| = N^{o(1)}.
  * Combined with IV.4: below defect n^{1/4} (resp. n^{2/5}), an integral fake is a Liouville-sharp genuine cluster for a form with N^{o(1)} < |disc| ≤ N^{3/10−ε'} (if |disc| is bounded, Question B for those finitely many K would kill it; note Question B for Q(√−2) and the like is *not* implied by Question B for Q(i)), and one within rank one of the Euclidean form.
  * Thus the residual case of the dichotomy is exactly the "middle-discriminant" regime that Question B, being a fixed-K statement, does not address. That is the precise reason it cannot help genuine tuples.
* **Significance.** For the goal: nil. For Question B: modest but real. It identifies which feature of the circle any proof must use (the discriminant of the form). It shows that the IV.2 fakes are honest clusters of growing-discriminant forms. It shows that the quarter-power barrier is optimal for polynomial pencils (numerically). It closes the "kill the lower-order terms" route exactly.
* **Novelty.**
  * R1 and R3: new to the record. R3 is the main new item.
  * The (b3) formula: new.
  * The σ-system reduction and its emptiness evidence: new.
  * R2 is routine.
