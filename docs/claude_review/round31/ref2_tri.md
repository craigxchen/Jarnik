# Round 31 referee report, block "tri": astra_free_II/answer.md, section II.5

Scope. Section II.5 of `astra_free_II/answer.md` (lines 607-637):
* the fake analogue (distinct labels in the angle model);
* the triangle completion formula (II.1), with c | N/(ab);
* the reverse parametrisation;
* the profile condition abc | N, gcd(a,b,c) = 1;
* the "six ordered records" enumeration;
* the Vieta descent and where it loses the parent;
* the computed example at N = 131490079981.

I also comment briefly on the PTE remark in the closing HEURISTIC paragraph.

Scripts are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref2_tri/`. All of them are my own code:
* `symbolic.py`: sympy proofs of every identity in the block.
* `genuine.py`: Gaussian arithmetic, lattice points of a circle, and a per-triple checker for (II.1), the reverse
  formulas, the profile condition and the determinant dictionary. It is run on the record's clusters and on random circles.
* `short.py`: Astra's triangle as an actual short triple, the 417-digit four-point cluster and the Vieta example.
  The construction recipe for the 417-digit cluster is copied from Astra's `near_square_check.py`. All the analysis is mine.
* `vieta.py`: the descent and the seed bound, checked by brute force.
* `reflect.py`: the reflected triple is a genuine triangle at a different parent.

I read Astra's `ii_triangle.py`. Its assertions agree with my independent output.

## 0. Summary verdict

* Every PROVED and COMPUTED statement in the block is correct. I found no false statement. Two points of wording
  are noted below; neither changes a conclusion.
* **(II.1) is not new.** It is the identity proved in `round27/stab.md` Prop 1.1(iii), the "binary quadratic form
  identity" `M_{b\a} t_a^2 - 2 y_ab t_a t_b + M_{a\b} t_b^2 = G_ab t_ab^2` of discriminant `-4 t_ab^2`. That
  identity was confirmed by the killer in `round27/kill_stab.md` (lines 69, 83, 455) and listed in `round27/README.md`
  §3 item 5. Astra's version is the same identity under the relabelling
  (G_ab, M_{a\b}, M_{b\a}, t_a, t_b, t_ab, y_ab) -> (c, a/b, b/a, u, v, w, s).
  The reverse formulas are the real-part relations `x_a x_b + t_a t_b = G_ab y_ab` and
  `x_b t_a - x_a t_b = G_ab t_ab` of the same proposition, solved for the coordinates.
* **Yes: (II.1) is the three-point case of the record's anchored determinant data.** It is not the three-point case
  of Question A's data (§2 below):
  * (II.1) is the condition det(P_1,P_2) = t_12 n_12 *plus* the Gaussian norm condition |P_i|^2 = prod_{T∋i} n_T.
  * The norm condition is exactly what Question A forgets (`growth/ptolemy.md` Thm B(iii)).
  * At k = 3, Question A's data are vacuous: there is no visible block with 2 <= |T| <= k-2.
  * gcd(a,b,c) = 1 is the chain condition at k = 3.
* The fake-analogue statement is correct. It concedes the coordinator's point 1 and extends the squarefree
  observation of `ref_divisor.md` to prime powers.
* The Vieta descent is correct and elementary. The reflected triple is in fact a *genuine* triangle at a different
  parent (new check below), which makes the "loss of the parent" concrete.
* **Significance: peripheral.** Nothing here changes:
  * the growth rate (2/3+eps) log R/loglog R;
  * the uniform bound;
  * the status of Question A, Question B or route 5.

  (II.1) plus its reverse form is an *equivalent* restatement of "these three labels form a genuine triangle at
  parent N". It is an exact-arithmetic condition that fakes fail (they have no integer s). It is not an inequality,
  and Astra gives no bound for the proposed count (HEURISTIC, correctly labelled).

## 1. Claim-by-claim

### 1.1 Fake analogue: distinct labels in the angle model (PROVED). Verdict: CORRECT

Setting. In the angle model, consider the edges v = a_i - a_j of the one-unit, shortened cluster that have
w(v) <= W/2 + B log W. There are >> m^2 of them, by Part I's Markov step.

* Since <v,phi> is within Delta/2 of pi Z and Delta/2 is about N^{-1/4}/4 = e^{-W/4}/4, the ratio to the Liouville
  bound e^{-w(v)/2} >= e^{-W/4} W^{-B/2} is at most W^{B/2}/4. So "within a factor O(W^{B/2})" holds.
* Suppose |v_p| = |v'_p| at every p and v != ±v'. Then each coordinate of (v±v')/2 is either v_p or 0.
  * The two halves are integral, have disjoint supports and are both nonzero.
  * Their weights add up to w(v).
  * So one half has weight <= W/4 + (B/2) log W. Astra's "<= W/2" is weaker than this but suffices.
  * Its angle is within Delta/2 of (pi/2)Z, while the axial gap gives distance >= e^{-w/2} >= e^{-W/8} W^{-B/4},
    which is far larger than Delta/2. Contradiction.
* v = ±v' with different edges is an affine relation among three or four rows: a_i - a_j - a_k + a_l = 0, or
  2a_i - a_j - a_k = 0. Rows are affinely independent in the fake model by II.4 (line 554). I checked that proof:
  * |x-y| <= x+y-2xy on [0,1]^2, because the difference is 2 min(x,y)(1-max(x,y)) >= 0;
  * the mean is 2 mu(1-mu) <= 1/2;
  * so some cross pair has weight <= W/2, against w(a_i - a_j) > W/2 from the gap.

Novelty: small. `ref_divisor.md` §1 already says that in the squarefree case D_B is the angle-model count, and
that the flipped-Paley fakes make it of order m^2. The prime-power extension (labels prod p^{|v_p|} distinct) is new
but routine. Significance: negative, and it confirms F1. Distinct labels and square-product triangle counts of
angle-defined labels are fake-satisfied. Only the literal integer identity d - floor(sqrt d)^2 = t^2 (with primitivity)
goes beyond fakes, as Astra says.

### 1.2 Labels ab, bc, ca with abc | N (PROVED). Verdict: CORRECT

* (ab)(ca)/(bc) = a^2, and similarly for the others.
* The local exponents are A+B, B+C, C+A. The range of three points equals A+B+C when min = 0.

My checks:
* Every genuine triple checked (§3) satisfies a = sqrt(d_01 d_02 / d_12), gcd(a,b,c) = 1 and abc | N.
* There a = N(gcd(P_1,P_2)) = n_12.

### 1.3 Formula (II.1) and its proof (PROVED). Verdict: CORRECT; not new

Statement: c = (a u^2 + b v^2 ± 2 s u v)/w^2, with c | N/(ab).

Symbolic checks (`symbolic.py`):
* det(Y,Z)X + det(Z,X)Y + det(X,Y)Z = 0.
* The determinant of the first two vectors is the product of the three dets.
* (A+B-C)^2 = 4AB - 4U.
* AB - U = s^2 u^2 v^2.
* Solving det(P_1,P_2) = ±a u with P_1 = (s,w), Y_2 = ±v gives exactly (II.1).

The sine-formula step is correct: |det(X,Y)| = sqrt(ab) sin psi_ab = w. It uses the closure psi_j = psi_k + psi_l
from Part I's triangle theorem, refereed CORRECT in `ref_divisor.md`. The closure is needed only for labels that
are not given as points. For actual cluster triples, the determinant derivation in §2 below is direct.

Both signs occur, so "±" is needed. Example: at N = 131490079981 the three anchorings give signs -, -, +
(`short.py`). In the 7- and 8-point record clusters, sign -1 occurs at the anchors sampled.

On offsets: the bound "at most W^{B/2}/4" follows from 16 t^2 sqrt N <= d <= sqrt N W^B. Correct.

### 1.4 Reverse parametrisation (PROVED). Verdict: CORRECT

* w^2(ac - v^2) = (au + eps s v)^2 and w^2(bc - u^2) = (su + eps b v)^2 are verified symbolically.
* The integrality argument is correct: a rational number with integral square is an integer.
* On every genuine triple, the recovered X_2 and Re Q match the actual |Re P_2| and |Re A_12| (`genuine.py`).

Wording point: the recovered representations need not be primitive or eligible. Astra imposes (1) afterwards, so
this is consistent.

### 1.5 Exact local profile condition abc | N and gcd(a,b,c) = 1 (PROVED). Verdict: CORRECT

* The largest pair sum equals the sum of the other two iff 2 min = 0.
* Realisation on [0, e_p]: put the extremes at 0 and A+B+C and the middle mark accordingly.
* For squarefree N, divisibility gives pairwise coprimality, so gcd(a,b,c) = 1 is automatic.
* This is the k = 3 chain condition. At p, the blocks {1}, {2}, {1,2} form two chains with disjoint unions unless
  p divides all three.

### 1.6 Enumeration and "six ordered records" (PROVED). Verdict: CORRECT

* There are 3 choices of which label is ab and 2 orders of (a,b).
* a, b, c come from square roots of quotients, the offsets from floors (Part I).
* The sign is unique because suv > 0, so the two values of c w^2 differ.
* Exhaustiveness follows from the triangle theorem.

"Little-oh of order (W/log W)^3 suffices" is Part I's triangle refinement. That was refereed CORRECT in `ref_divisor.md`.

Comment. The enumeration gives only the trivial bound
#triangles <= sum_{d eligible} tau(d) · 2 (W^{B/2}/4)^2,
which is at least D_B(N) · W^B/8. That is *worse* than the target. So the parametrisation is a search recipe, not an
estimate. Astra does not claim otherwise.

### 1.7 Vieta descent (PROVED). Verdict: CORRECT

The equation (A+B-C)^2 = 4AB - 4U is quadratic in C. Its roots have sum 2(A+B) and product (A-B)^2 + 4U > 0
(sympy). So C' is a positive integer.

If C > A+B (C largest), then C' < A+B < C, and the entry sum decreases strictly.

At the stopping point:
* 0 <= A+B-C <= A, so 4AB - 4U <= A^2 <= AB, giving 3AB <= 4U.
* With A >= 1 this gives A+B <= 8U/3.

The step A >= 1 is implicit in Astra's text. It holds because all entries stay positive integers.

Brute force (`vieta.py`):
* For U <= 200 there are no seed-bound violations.
* The seed counts for U = 1..12 are 1,1,2,2,2,2,2,3,3,2,3,4.
* Astra's triple (5668, 20381, 4553) descends to the seed (1,4,5) in 13 steps.

The "loss of the parent" is correct, and it can be made sharper:
* The reflected entry need not even keep the factored form c' w^2. Mod w^2 one has C' ≡ ∓4suv.
* When it does keep the form, the reflected triple is a genuine triangle at a *different* parent.
  * In Astra's example, (a, b', c) = (1417, 61, 4553) has labels 86437 = 294^2 + 1, 277733 = 527^2 + 2^2 and
    6451601 = 2540^2 + 1, all primitive.
  * `reflect.py` finds 8 genuine realisations among the lattice points of x^2 + y^2 = 393547661 = 13·29·61·109·157.
  * The angular cutoff 16 t^2 sqrt N' <= d fails there for 86437.

So the descent moves along a tree of genuine triangles whose parent changes at each step. A count organised by seeds
therefore loses the common-parent constraint, which is the whole content.

### 1.8 COMPUTED example at N = 131490079981. Verdict: REPRODUCED

* N = 1417·20381·4553 = 13·29·89·109·157·229.
* (a,b,c) = (1417, 20381, 4553) and (u,v,w) = (2,1,1).
* (A,B,C) = (5668, 20381, 4553), U = 4, and Heron's relation holds.
* The reflection of B is 61, with N mod 61 = 40. Since v = 1, b' = 61 would have to divide N/(ca) = 20381 = 89·229,
  which it does not.

I also located this triangle as an actual short triple of lattice points on x^2 + y^2 = N. All three anchorings pass
the checker, and the records are (1417,20381,4553; u,v,w = 2,1,1; s = 5374; -), (20381,4553,1417; 1,1,2; 9633; -) and
(4553,1417,20381; 1,2,1; 2540; +).

### 1.9 PTE remark (HEURISTIC paragraph). Verdict: CORRECT (trivial)

Im B constant means B - conj(B) is a nonzero constant. So B and conj(B) share no root and have equal e_1, ..., e_{k-1},
hence equal power sums p_1, ..., p_{k-1}: an ideal PTE pair over Q(i). This is peripheral and I did not grade it further.

## 2. Is (II.1) the three-point case of the determinant formulation?

**Yes, but of the full Gaussian anchored data, not of Question A's data.**

Take three points of one unit class. Anchor z_0 at the vertex where the edges labelled ab (offset w) and ca
(offset v) meet. In the notation of `growth/ptolemy.md` §1 with m = 2:
* P_1 = A_ab = (s, ±w) and P_2 = A_ca = (X_2, ±v);
* the blocks are n_12 = N(gcd(P_1,P_2)) = a, n_1 = b, n_2 = c;
* the residues are t_01 = Y_1 = ±w, t_02 = Y_2 = ±v, and t_12 = det(P_1,P_2)/n_12 = ±u.

So the offset u of the third edge is the residue t_12. This is verified on every triple below, including
u = |Im A_12| computed independently from z_2/z_1.

Then (II.1) is precisely
|P_2|^2 = n_12 n_2 after eliminating X_2 = (s Y_2 ∓ n_12 t_12)/w from det(P_1,P_2) = t_12 n_12,
together with |P_1|^2 = n_12 n_1. Equivalently it is Lagrange's identity
|P_1|^2 |P_2|^2 = (P_1·P_2)^2 + det(P_1,P_2)^2.
gcd(a,b,c) = 1 is the chain condition, and abc | N is prod n_T | N.

Two caveats.
* Question A (L_8), the determinant statement in expert_note §2, keeps only the minors and the chain condition.
  Ptolemy Thm B(iii) states that this system forgets the Gaussian norm condition |P_i|^2 = prod_{T∋i} n_T, the
  isotropy and the singleton blocks. (II.1) *is* the norm condition, so it lies entirely in the forgotten part.
* At k = 3 there are no visible blocks (2 <= |T| <= m-1 = 1 is empty), and M_{0,3} is a point. Question A's data
  have no content there.

The exact identity, with the same discriminant -4t^2 and the same roles, is already in the record as
`round27/stab.md` Prop 1.1(iii). That statement was proved for one-class tuples on (Ch) and confirmed in `kill_stab.md`.

What the divisor/quadratic-form phrasing adds:
1. It anchors the third block to the fixed parent: c must divide N/(ab), not just be some positive integer. In the
   record, N is implicit in the blocks.
2. It gives a reading as a *completion count*. Fix one eligible d = ab and its factorisation. Each pair of small
   offsets (u,v) <= W^{B/2}/4 gives at most two candidate c. The question becomes how many of the (a,u,v) values of
   [a, ±2s, b] are w^2 times a divisor of N/d.
3. It is exhaustive for the Part I triangle target (via the triangle theorem), with multiplicity <= 6.

That is a clean search parametrisation, and it is genuinely arithmetic: fakes have no integer s. But it adds no
inequality, and no counting bound is offered. It is an equivalent restatement of genuineness, in the same sense as
Part I's five-vertex reconstruction.

## 3. Computations on genuine triangles (own code)

`genuine.py` takes the anchored primitive P from the half-angle vector (N + Re G, Im G) / content, with
G = z_i conj(z_0). Triples whose relative unit is ±i are skipped (even norm). For each remaining triple, it asserts all
of the following:
* the label identities, a = n_12, gcd(a,b,c) = 1 and abc | N;
* u = |det|/a = |Im A_12|;
* (II.1) with a determined sign;
* the reverse formulas, with the recovered coordinates equal to the actual ones;
* Heron's relation.

Results:

| source | triples passing (same unit class) | failures |
|---|---:|---:|
| record 8-point cluster, R^2 = 32988780417125 (`growth/families.md`, C = 5.674) | 24 anchored triples | 0 |
| record 7-point cluster, R^2 = 235370525 (C = 5.270) | 30 | 0 |
| 417-digit genuine four-point cluster (Part I) | 12 of 12 anchored triples | 0 |
| N = 131490079981: Astra's short triangle (3 anchorings) + 4834 random triples | all | 0 |
| random circles N = 1185665, 359125, 20480525, 359125 | 4585 | 0 |

In total 9473 triples were checked and none failed. The offsets in the 417-digit cluster (4096, 9216, 35840, 1024, ...)
square to the values in Astra's Part I list.

## 4. Novelty and significance

* (II.1): **not new**. It is `round27/stab.md` Prop 1.1(iii), and Astra does not cite it. The reverse formulas are
  that proposition's linear relations. The profile condition is the k = 3 chain condition.
* Fake analogue for prime powers: new but routine, and it confirms F1.
* Vieta descent: new to the record (grep finds no Vieta or Markov-type descent on this equation), correct, elementary,
  and of no use for counting. As shown above, it walks to genuine triangles at other parents.
* Significance for the goal: **peripheral**.
  * There is no change to the growth rate, the uniform bound, or Questions A, B or route 5.
  * The one research suggestion is to bound the number of small quadratic-form values dividing one fixed complement,
    summed over d. It is unproved, labelled HEURISTIC, and needs an estimate of a kind that neither Astra nor the
    record has.
