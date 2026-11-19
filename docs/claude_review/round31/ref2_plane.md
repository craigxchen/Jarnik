# Referee report ref2_plane: astra_lemma_II, sections II.6 and II.7 (fixed rational planes of norm vectors)

Scope: `astra_lemma_II/answer.md` II.6 (finitely many bounded-cofactor tuples on a fixed rational plane, the conic
reduction, and the real five-star lemma) and II.7 (the power lower bound `T >>_U (max n)^{1/6}`, or `1/4` for
pairwise nonproportional forms, the exponents 38/11 and 40/11, the arithmetic-progression assignments, the prime
progression `34913+2100r`, the generic quartet, and the moving progression).
Scripts: `round31/ref2_plane/`. All of them are my own code. I did not rerun Astra's scripts.

| script | what it checks | output |
|---|---|---|
| `ap_and_identities.py` | (a) the 10,080 anchor triples; (b) a second, independent search over all 9! = 362,880 slope assignments with slope 0 on edge 12; (c) primality and residues of `34913+2100r`; Gaussian blocks and nonzero `t_ij` in 5 random edge assignments; (d) (II.10), (II.6) for all ordered triples, (8) and all five Plücker relations, exactly, on 200 random genuine ten-pair tuples | `ap_out.txt` |
| `fivestar.py` | symbolic: `H'/H`, the equal-critical-value determinant, `b=-16/9`, `H'` at `b=-16/9` (the numerical part of this script was superseded by the next two) | first block of `fivestar_out` (printed in the log below) |
| `fivestar_search.py` | vectorised Gauss–Newton search for real and complex solutions of "constant star sums and constant reciprocal star sums" | `fivestar_search_out.txt` |
| `fivestar_square.py` | the same with two coordinates pinned (`y_e0=1`, `y_e1=a`, 56 values of `a`, 4 pin choices), which removes the attraction to the constant solution | `fivestar_square_out.txt` |
| `genus.py` | the genus-five curve of II.7: 3×3 minors, a Jacobian criterion via a Gröbner basis, F_p point counts, Riemann–Hurwitz | `genus_out.txt` |
| `p3_algebra.py` | symbolic check of `prod_j t_ij = λ^4 u_i^3 (prod u) D_i`, the perfect-matching identity for all three matchings at every vertex, and the minors of `W_i` | stdout (all True) |

## 0. Verdict in brief

* **II.6 and II.7 are correct.** I checked every step by hand and found no gap. The two steps Astra flags hold:
  the transfer of (8) to (P3) with the perfect matching, and the Pell exclusion by the two cusps. Both are correct
  as written.
* The genus is 5. I recomputed it two ways: as a smooth complete intersection of three diagonal quadrics in P^4, and
  by Riemann–Hurwitz.
* Exponents. Pairwise nonproportional forms give cofactor exponent 1/4. Proportional forms give 1/6, but there the
  approximation bound is `H^{-2/7}`. The approximation exponents are 38/11 in general and 40/11 for comparable norms.
* Arithmetic progressions: 10,080 triples and 0 solutions, confirmed by a fuller independent search. The ten numbers
  `34913+2100r` are primes ≡ 1 mod 4.
* Corrections, all cosmetic:
  1. The determinant in the five-star lemma has the opposite sign: it is `-b^2(9b+16)(24b+1)/16`. The zeros are
     the same.
  2. Roth is needed only for a squared irreducible cubic factor and for simple roots of degree ≥ 3. Liouville covers
     the quadratic factors.
  3. The exponent statements need the blocks to be private, so that `H=∏ n_e` is the primitive norm. This is
     automatic in Part I's model but is not repeated in the II.7 headline.
* **Significance: peripheral.**
  * "Norm vector on a fixed rational plane" is a codimension-8 linear condition on the ten block norms. The
    constants are not uniform in the plane, and every single tuple lies on some rational plane. So the theorem
    constrains only sequences that stay on one plane. Nothing in the growth regime produces such sequences.
  * The ten-pair profile itself excludes singleton cuts, so it lies outside the balanced domain of `(L_5°)`.
  * The nondegenerate case gives only 40/11 < 15/4. Question B (`A_5(δ)` for some `δ<1/30`) needs exponent above
    15/4. Only the elementary Plücker-nonvanishing cases reach 4.
  * No change to the growth rate, the uniform bound, Questions A and B, or route 5.

## 1. Setting (checked)

For the one-unit-class ten-pair model of Part I §6 the definitions are:
* `P_i=∏_{j≠i}β_ij`, with pairwise coprime block norms `n_e`.
* `t_ij=det(P_i,P_j)/n_ij` and `s_ij=Re(P̄_iP_j)/n_ij`, both integers.
* (II.10): `s_ij^2+t_ij^2=∏_{a∉{i,j}}n_ia n_ja`. The six edges in this product are exactly those with
  `|e∩{i,j}|=1`.
* (II.6): `2t_ik t_jk s_ij=A_ijk`.
* (8): `Σ_i P_i/∏_{j≠i}t_ij=0`.

I verified all of these, with signed `t`, on 200 random genuine tuples. I also verified the five Plücker relations
for `M_ij=t_ij n_ij=det(P_i,P_j)`. (8) is the degree-four Lagrange coefficient identity
`Σ F(P_i)/∏_{j≠i}det(P_i,P_j)=0` for the cubic `F(P)=|P|^2P`, combined with `|P_i|^2=∏_{j≠i}n_ij`. It is correct
(see also ref2_window).

## 2. II.6, from bounded cofactors to a conic: CORRECT

Fix the signed `t` along an infinite subsequence. There are finitely many choices.

* **Dimension one.** At most one positive pairwise coprime vector lies on a rational ray.
* **Dimension two, setup.** `F(p,q)=∏L_iaL_ja-s_ij(L)^2` is a rational binary sextic. It is nonzero because it
  takes the value `t_ij^2≠0`.
  * Its values are infinitely many points on infinitely many rays: a degree-6 homogeneous form takes a fixed
    nonzero value at most twice on a line through 0.
  * Each primitive irreducible integral factor divides a fixed integer. So on a further subsequence each factor
    takes a fixed value.
* **Factor of degree ≥ 3.** Thue: finitely many solutions (Thue 1909, for irreducible binary forms of degree ≥ 3).
* **Two independent linear factors.** These pin `(p,q)`.
* **One linear factor.** The points lie on an affine line. All identities (II.10) and (II.6) become polynomial
  identities on that line.
  * Some `L_ab` is nonconstant there, because the parameter map is injective.
  * At its real root, (II.10) for the pair `(a,j)` with `j∉{a,b}` reads `s^2+t^2=0` with `t≠0`. This is
    impossible. I checked that `n_ab` does occur in that product.
* **Two nonproportional quadratics with fixed values `d_1,d_2`.** `d_2Q_1-d_1Q_2` vanishes on infinitely many
  rays, so `Q_1` and `Q_2` are proportional.
* **Definite quadratic.** Finitely many representations.
* **What remains** is exactly `F=κQ^3` with `Q` indefinite and irreducible over Q, so `Q(p,q)=d` (the Pell case).

**Conic step.**
* With `Q=uv` over R, `(p,q)=αu+β/u`, so `n_e=a_e u+b_e/u`.
* `a_e b_e≠0`, because the isotropic directions are irrational and `L_e` is rational and nonzero.
* Positivity on the chosen branch: if `a_e b_e<0`, then `n_e` has a positive real zero, and (II.10) for a pair
  separated by `e` gives a contradiction. Since `n_e>0` at the samples, `a_e,b_e>0`.
* Proportional forms are excluded by coprimality. Hence the ten ratios `b_e/a_e` are distinct and positive.
* The identities extend to the whole conic because they hold at infinitely many points of the irreducible curve.

**Divisor step.**
* `V_ij=(s+it)/(s-it)` equals `(P_j/P̄_j)/(P_i/P̄_i)` at the samples, so it is a cocycle.
* `s` has poles of order exactly 3 at 0 and ∞, because `s^2=∏n-t^2` has order 6. So `V_ij=1+O(u^{∓3})`.
* The zeros of `(s+it)(s-it)=u^{-6}∏(a_eu^2+b_e)` are the twelve simple points `±i√(b_e/a_e)`. Each is a zero of
  exactly one factor, and real conjugation swaps the zero and the pole.
* A cocycle with value 0 on same-side pairs and ±1 across the cut of `e` is `ε_e(1_{j∈e}-1_{i∈e})`. This gives
  the displayed product with `V(∞)=1`.
* Expanding `log((u-iy)/(u+iy))`:
  * At ∞ it is `-2iy/u+O(u^{-3})`, an odd series, so the `u^{-2}` term vanishes automatically.
  * At 0 it is `iπ+2iu/y+…`. The constant terms cancel because `deg_j-deg_i=0`.
  * So exactly the star sums and the reciprocal star sums of the real `y_e=±√(b_e/a_e)` are constant.
* The `|y_e|` are distinct, as the lemma requires.

## 3. The real five-star lemma (flagged step 2): CORRECT

I went through every step.
* **Quartics.** `P_i(z)=z^4-sz^3+c_iz^2-rd_iz+d_i`, because `e_3/e_4=Σ1/y`. So `P_i(z)=0` exactly when
  `C(z)=[1-rz:z^4-sz^3:z^2]` lies on the line `w_1+c_iw_2+d_iw_0=0`.
* **Birationality.**
  * The forms `[w^4-rzw^3 : z^4-szw^3 : z^2w^2]` have no base point.
  * If the map is not birational, its degree must be 2 and `X=z^{-2}-rz^{-1}` must be invariant under
    `z↦s-z`. This forces `s=0`, and then `r=0`.
  * `r=s=0` makes every `P_i` even, which contradicts "distinct up to sign".
  * So `C` is an irreducible rational quartic, and its total δ is 3.
* **The two unibranch singular points.**
  * At `z=0` the local orders are (≥3, 2) in the chart `w_0=1`. At `z=∞` they are (≥3, 2) in the chart `w_1=1`.
    So the differential vanishes at both.
  * `w_2=0` only at these two parameters, which map to the distinct points `[1:0:0]` and `[0:1:0]`.
  * Both points are real.
* **The five lines.**
  * They are distinct, since the stars share exactly one value.
  * `ℓ_i∩ℓ_j=C(y_ij)`.
  * No three are concurrent, because `Y=z^2-sz` has degree 2 and the `y` are finite and nonzero.
  * So there is a unique nondegenerate real conic tangent to all five lines, and its dual has real points.
* **Sym² step.**
  * `D(h_i,y)=k_if(y)/(y-h_i)`. Interpolation gives `D=(f(x)g(y)-f(y)g(x))/(x-y)`, and I checked the identity at
    `x=h_i`.
  * `gcd(f,g)=1`, since `g(h_i)=k_i≠0`. So `h=f/g` is real of degree 5, with five distinct real points over 0.
* **Local equations.**
  * At a pair with local degrees `(m,n)` the local equation is `X^m=Y^n`. With `m+n≤5`, the only unibranch
    singular case is `(2,3)`, an ordinary cusp. `(2,2)` gives a node.
  * At a diagonal point: `m=2,3` are smooth. `m=4` gives `w(w^2-2v)` (two branches). `m=5` gives
    `w^4-3w^2v+v^2` (two branches, since `x^2-3x+1` has distinct roots).
  * So each of the two unibranch points is a `(3,2)` fibre, and the two lie in different fibres. Both critical
    points in each such fibre are real, because conjugation preserves local degree.
* **Ramification budget.**
  * The two `(3,2)` fibres use 3+3 of the 8 units. δ still owes 1.
  * Of the ways to spend the last 2 units, only a `(2,2,1)` fibre (a node) produces a singularity. The others are
    `(3,1,1)` (smooth points only, since diagonal `m=3` and `(3,1)` are smooth) or two simple critical points
    over distinct values (smooth).
* **Normal form.**
  * Real Möbius maps on the source and target give `H=x^3/((x-1)^3(x-b)^2)`.
  * `H'/H` has numerator `3b-2x^2-x`. My sympy run confirms this.
  * The equal-value condition is `a_1b_0-a_0b_1=-b^2(9b+16)(24b+1)/16`. **Sign correction:** Astra writes `+`.
  * `b=0` is excluded, and `b=-1/24` gives a double critical point, which means no node. So `b=-16/9`.
  * Then the critical points are `-1/4±5√15 i/12`, which are nonreal, and `H'∝x^2(6x^2+3x+16)`.
* **Conclusion.** The only real critical values come from the cusp fibres. Near 0 there are 1 + (0 or 2) real
  preimages on each side. So every regular real fibre has at most 3 real points, which contradicts the five real
  points over 0.

**Numerical cross-checks.** These are evidence, not proof.
* Unpinned search: 19,881 converged real runs.
* Pinned search: 13,244 converged real runs, over 56 values of the ratio `y_e1/y_e0` and 4 pin choices.
* Every converged real solution was degenerate: some `|y_e|` repeated, with largest normalized minimal gap
  `≈2e-9`. Examples: the constant vector, or ± patterns on two values.
* The complex search finds nondegenerate solutions. So the reality hypothesis is essential and the lemma is not
  vacuous over C. This matches the proof, which over C forces `h≅H_{-16/9}`, a one-parameter family up to scaling.

Novelty: I found no five-star or reciprocal-star-sum statement in the record. The nearest is the linear-block
moment problem (growth/cotangent.md). **New.**

## 4. II.7, the power lower bound: CORRECT

**Proportional forms.**
* Coprimality fixes the two norms, so the parameters lie on a fixed affine integral progression, determined by `U`.
* (P1) is an integral polynomial in `v` with coefficients `O_U(T^6)`.
* It is nonzero: its value at the real root of `n_ab` is `-A^2-4t^6<0`.
* Cauchy's bound then gives `B<<_U T^6`. This is effective.

**Plücker step.**
* Each quartet relation is an integral binary quadratic with coefficients `O_U(T^2)` that vanishes at the primitive
  `(p,q)`. `gcd(p,q)=1` follows from pairwise coprimality.
* If a relation is nonzero, the rational-root theorem gives `B<<T^2`.
* Otherwise `M=pM_1+qM_2` is decomposable identically. Then:
  * `M_1` and `M_2` are independent and nonzero, since otherwise the forms would be proportional.
  * `M_1∧M_2=0`, so `M=u∧(pv+qw)` with `u` rational.
  * `u_i≠0`, since otherwise the four forms at vertex `i` would be proportional.
  * Hence `t_ijL_ij=u_iu_j⟨A_j-A_i,(p,q)⟩`.
* Ten distinct directions determine the configuration up to translation and homothety. This gives (P2), with
  `c_ij` and `A^0` depending only on `U`.

**(8) to (P3) (flagged step 1).** This step is correct.
* `det(W_i,W_j)=M_ij/λ`. So the two rank-2 matrices `W` and `P` differ by a matrix in `GL_2(R)`, and the linear
  relation (8) transfers.
* `∏_{j≠i}t_ij=λ^4u_i^3(∏u)D_i`, which gives both coordinates of (P3).
* For each `i` and any perfect matching `{ab,cd}` of the other four vertices:
  `(t_ab/c_ab)(t_cd/c_cd)=λ^2(∏u)/u_i`. I checked this symbolically for all three matchings at every vertex.
* So `1/u_i^2∝m_i^2`, where `m_i` has bounded denominators and size `O_U(T^2)`.
* The second (P3) equation becomes an integral linear form in `(p,q)` with coefficients `O_U(T^4)`. If it is
  nonzero, primitivity gives `B<<T^4`. This case gives the exponent 1/4.

**Genus five (recomputed).**
* The points `[1/u_i]` are rational points of `{[y]: Σ y_i^2 col_i=0}` with `col_i=(1,A_i^0)/D_i`. This is a
  complete intersection of three diagonal quadrics in P^4.
* It is smooth iff every 3×3 minor of the coefficient matrix is nonzero. If at most two coordinates are nonzero the
  equations force `y=0`; with three or more, the Jacobian has rank 3.
* Here all 3×3 minors are nonzero because no three of the `A_i^0` are collinear.
* Adjunction: `K=O(2+2+2-5)=O(1)` and degree 8, so `2g-2=8` and `g=5`.
* Riemann–Hurwitz on the `(Z/2)^4` cover of `Λ≅P^1`, branched at the five coordinate zeros, also gives
  `-32+40=8`. The cover is connected because the four squareclasses are independent.
* On a sample configuration I confirmed smoothness with a Gröbner basis (`y_i^12` lies in the ideal of the quadrics
  plus the Jacobian minors). F_p point counts are consistent with the Weil bound.
* Faltings (Mordell conjecture over number fields, genus ≥ 2) gives finitely many points. These depend only on `U`.

**Size for a fixed pattern.**
* `F=λ^2` with `F` fixed. Near a real root of multiplicity `m`:
  * a rational root gives distance `>>1/B`, which is enough for `m≤5`;
  * an irrational root gives distance `>>B^{-9/4}` by Roth with ε=1/4, which is enough for `m≤2`.
* Liouville already suffices for quadratic factors. Roth is needed only for `f^2` with `f` an irreducible cubic,
  and for simple roots of degree ≥ 3. The statement of Roth used, `|α-p/q|>c(α,ε)q^{-2-ε}`, is the standard one.
* So `|λ|>>B^{1/2}`, except in the shapes `cℓ^6` and `cQ^3`.
* For those two shapes, rescale by `k=|λ|^{-1/3}`. This multiplies `n` by `k`, `P` by `k^2` and `t` by `k^3`, and
  (II.6) and (II.10) are homogeneous. The rescaled points lie on a fixed real line (the root argument applies) or
  on a fixed hyperbola (the II.6 conic argument and the five-star lemma apply, since integrality is not used there).
  So only finitely many samples remain.

**Uniformity and ineffectivity.**
* The constants depend on `U` (its lattice basis and edge labelling) and on nothing in the tuple.
* Ineffectivity enters only in the stratum where all Plücker forms vanish and `V=0`. There it comes from Faltings
  (the finite pattern list), Roth's constant, and the finite exceptional sample sets.
* The proportional, Plücker-nonzero and `V≠0` cases are effective.

**Exponents (recomputed).** I write `X`, `B`, `Y` for the largest norm, the parameter size and the smallest norm.
* Comparable norms: `Δ>>X^{1/4}/X^3=X^{-11/4}=H^{-11/40}`, so exponent 40/11 ≈ 3.636.
* General nonproportional case:
  * nine norms are `≍B` and `H≍B^9Y`;
  * the maximum pair gives `Δ>>B^{-11/4}`, and a pair separated by the cut of norm `Y` gives
    `Δ>>B^{-5/2}Y^{-1/2}`;
  * the weights 8/19 and 11/19 give `(B^9Y)^{-11/38}`;
  * the worst case is `Y=B^{1/2}`, where the two bounds coincide.
  * So the exponent is 38/11 ≈ 3.455.
* Proportional case: some pair is separated by at most `⌊3h/5⌋≤4h/7` nonconstant cuts. For `h=5`, `M1_S=3·1` would
  force `1_S=1/2`. This gives `Δ>>H^{-2/7}`, which is stronger than `H^{-11/38}`. CORRECT.

**Arithmetic progressions.**
* My anchor enumeration gives **10,080 distinct triples and 0 solutions**. My second search does not use the anchor
  reduction: it puts slope 0 on edge 12 and tries all 9! assignments of the other slopes. It also finds **0**
  configurations.
* So some Plücker form is nonzero in every edge assignment, and `T>>√X` with an absolute constant gives exponent 4
  under comparability. CORRECT.
* The generic-quartet determinant is 32 (checked by hand), and the `X_4^{1/2}K^{-3/2}` bound is correct.
* The moving progression bound `T>>√X/K` is correct, with the stated caveat that the three complementary products
  must be independent.

**The prime progression.** `34913+2100r`, for `0≤r≤9`, are all prime and all ≡ 1 mod 4, and `gcd(34913,2100)=1`.
In 5 random edge assignments with random Gaussian generators, all ten `t_ij` are nonzero. CORRECT.

**Green–Tao.** "Linear equations in primes" (Annals 2010), with the Möbius–nilsequence and inverse theorems of
Green–Tao and Green–Tao–Ziegler (2012), covers the ten forms `4x+4ry+1`. Their linear parts are pairwise
independent and the local factors are positive. The statement used is the one that exists. It is needed only for
nonemptiness at all scales, not for any bound.

## 5. Novelty and significance

**Novelty.** The fixed-plane theorem is new to the record, and so are the five-star lemma and the genus-five
reduction. grep finds nothing on fixed planes of norm vectors, star sums, or this curve. The record's "Faltings"
hits concern other curves. Note also that `tenpair_norm_pencil_power.md` in Astra's directory has an older weaker
version (cofactor exponent 1/10, approximation exponent 100/29). `answer.md` supersedes it through the
perfect-matching trick.

**How restrictive is "norm vector on a fixed rational plane"?**
1. **Codimension 8.** Inside the ten-dimensional space of block-norm vectors this is a linear condition of
   codimension 8. A fixed `U` contains `≍X^2` of the `≍X^10` vectors of size `X`.
2. **No uniformity in `U`.** Every single tuple lies on rational planes, for example `span(n,1)`, but these have
   height about `X`. All constants are ineffective functions of `U` (Faltings, Roth). So the theorem says nothing
   about a sequence of tuples unless their norm vectors stay on **one** plane.
   * In the growth regime the ten block norms of a five-subtuple are products of many primes of uncontrolled
     size. No mechanism imposes eight fixed linear relations among them.
   * Astra itself notes that (II.2) controls only logarithmic sizes. I agree with Astra's EVIDENCE paragraph:
     there is no implication from a small cofactor product to a plane of controlled height.
3. **Narrow profile.** The ten-pair profile is a narrow subclass of index-one tuples. II.2 shows that almost all
   growth-scale five-tuples have index exactly one, but with all fifteen cuts and with prime powers.
   * The ten-pair model needs every prime's five-row allocation to be a two-mark cut. So there are no singleton
     cuts, and the five singleton divisors are absent.
   * Hence the model is outside the balanced 15-cut domain of `(L_5°)`.
   * With singleton blocks, `P_i` acquires a factor `γ_i` and (8) changes. II.6–II.7 do not cover this.
4. **Below the Question B threshold.** Question B needs `A_5(δ)` for some `δ<1/30`, that is, exponent above 15/4.
   The theorem's hard case gives 40/11 < 15/4 with comparable norms, and 38/11 in general. Only the Plücker-nonzero
   cases (arithmetic progressions, the generic quartet), which are elementary, give 4.

**Significance for the goal.** None of the following changes:
* the growth rate `(2/3+ε)log R/log log R`;
* the uniform bound;
* Question A or Question B;
* route 5.

The theorem is a correct, carefully proved, peripheral rigidity statement. Its real interest is the observation in
the polynomial and rigidity direction: on a fixed plane, the bounded-cofactor configurations collapse to a
star-sum problem with no real solution. That also explains why the polynomial no-pencil result (II.5) has an
arithmetic shadow. **Score on the goal: 1/10. Mathematical quality: high.**

## 6. Claim-by-claim

| claim | verdict |
|---|---|
| II.6 Thue / linear / two-quadric / definite reductions to Pell | CORRECT |
| II.6 conic step and `V_ij` divisor to constant star and reciprocal star sums | CORRECT |
| Real five-star lemma | CORRECT, with the determinant sign corrected to `-b^2(9b+16)(24b+1)/16` |
| II.7 proportional forms `B<<T^6` | CORRECT (effective) |
| Plücker line to (P2) | CORRECT |
| (8) to (P3), perfect matching, `B<<T^4` | CORRECT |
| Genus 5 and Faltings | CORRECT (genus recomputed two ways) |
| Roth and exceptional shapes | CORRECT (Liouville covers the quadratic factors) |
| `T>>_U max n^{1/6}`, or `^{1/4}` for nonproportional forms; constants depend only on `U` | CORRECT; ineffective only in the degenerate stratum |
| Exponents 38/11 and 40/11 (and `H^{-2/7}` for proportional forms) | CORRECT, given private blocks so that `H` is the primitive norm |
| Arithmetic progressions in every assignment: exponent 4 | CORRECT (10,080 and 0 reproduced; full 9! search also 0) |
| `34913+2100r` | CORRECT |
| Generic quartet (determinant 32) and moving progression | CORRECT |
| Significance | peripheral |
