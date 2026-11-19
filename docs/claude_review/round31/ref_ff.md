# Round 31 referee report, block "ff": Astra free audit §7 (squarefree D=8, c=0 function-field test)

Referee scripts (all independent of Astra's and the judge's code): `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref_ff/`
- `census_ref.py`: orbit BFS of pair-admissible row sets under S_8 x global flip. Canonical form: min over all 40320 column permutations and the flip of the sorted row tuple.
- `filter_ref.py`: a complete, exact Lemma-1.1 test (not a box test), plus comparison with `round30/judge_checks/classes_D8_c0.pkl`.
- `count_w.py`: statistics of the rejection witnesses.
- `nonlin_ref.py <g>`: my own nonlinear exclusion. The gauge row g is chosen by hand (not RREF). Compatibility polynomials come from the left kernel of the quadratic design matrix.
- `control.py`: positive control, the same algebra run on H built from actual rational roots.
- `size_ref.py`, `size3_ref.py`: BFS over sets of rank <= 4 (resp. <= 3).
- `eight.py`: the 8-row orbits.
- `cube.py`: an exact check of the attaining cube cluster.
- `mapidx.py`: maps my class indices to the judge's.

All runs used `prlimit --as=3500000000`. Each run takes seconds.

## 0. Verdict

**CORRECT.** Astra's theorem holds as stated: every cluster in the squarefree D=8, c=0 model (over C(t), real t, ordinary conjugation) has affine rank <= 4 and at most 8 points. 8 is attained, and the 8-point sign code is unique up to S_8, row order and the global flip. The proof has two parts:
- Elementary lemmas, each of which I checked line by line.
- Finite exact computations, each of which I reproduced with independent code. The census, the 21 survivors, the exclusion of all 21, the rank-<=4 maximum and the uniqueness all agree.

**So the judge's M(8,0)=8 (squarefree N) now stands as a theorem.** It now has a better basis than the judge's: size <= 8 for rank <= 4 needs only pairwise balance, and rank >= 5 is excluded. The judge's saturated-filter maximum is no longer needed.

**Scope (correction of emphasis, not of content).** The result covers squarefree N in characteristic 0 with an ordered real place: C(t), and hence Q(i)(t). It does not cover the F_q(t) version, because the Rolle and positivity steps use the order on R. It does not cover non-squarefree N. Clusters with a nontrivial common factor are covered whenever N is squarefree: they appear as codes with constant columns, and the enumeration includes those.

**Novelty.** New. Round 30 §7 left exactly this as the open follow-up. Expert note line 634 says "what excludes them, if anything, is non-linear". The two geometric third-moment lemmas (proportional rows; four collinear rows of H) are new tools in the record.

**Significance.** Peripheral. It is the first full-D rank theorem in the function-field model, but only at one D. It is consistent with (R_4) and does not prove it for any other D. It changes neither the growth rate (2/3+eps) log R/loglog R nor the uniform bound. It does not change the status of Questions A or B or of route 5. Astra says so herself (EVIDENCE, relevance), and I agree.

## 1. Re-derivation of the 21 surviving classes (task item 1)

**Setting checked against `round29/ff.md` §1.** K=4, L=3. By round-29 Lemma 1.1, the cluster condition for monic degree-8 polynomials is that the first three power sums agree, i.e. sum_j (s_j - s'_j) Im(rho_j^k) = 0 for k=1,2,3. At L=3, round-30 Lemma 1.1 says every nonzero integer a in the span V of the differences has min(a+,a-) >= 2 (from >= 3/2), and therefore ||a||_1 >= 4.

**Census (my code).** Pair-admissible orbit counts for 1 to 6 rows: **5, 14, 26, 50, 65, 72**. Of the 72 six-row orbits, **64** have affine rank 5. Astra's counts are identical. Completeness: removing a row from an admissible k-set gives an admissible (k-1)-set, and the pair condition is invariant under S_8 and the global flip. Hence extending each orbit representative by every compatible row reaches every k-orbit.

**Filter.** I did not use a box test. Instead I decided exactly whether V ∩ Z^8 contains a nonzero a with min(a+,a-) <= 1:
- One-signed a: LP (Gordan).
- a = b - e_k with b >= 0 integer and b_k = 0: an exhaustive DFS bounded by sum_j b_j y_j = y_k, where y = Hu > 0 is an exactly verified positive vector.

**Results.**
- 43 classes are rejected. For each, the witness was checked exactly (H^T a = 0 and min mass <= 1). Every witness has the form b - e_k, with entries of absolute value <= 2 (41 with max entry 1, 2 with max entry 2).
- **21 classes survive.** As a set of canonical keys they equal the 21 classes of `classes_D8_c0.pkl`. They also agree with Astra's regenerated census (`exact_plane.py` prints "Independent census matches the 21 archived classes").

**Strengthening.** Because my filter is complete, the 21 survivors have **no** Lemma-1.1 violation of any size. Both the judge's box [-4,4] and Astra's box [-2,2] were therefore already sufficient here. The linear tool is exhausted at D=8.

I did not rerun the judge's `sym8.py` itself. It is not needed: my census plus the exclusions below proves the judge's conclusion by a different route.

## 2. The exclusion argument (task item 2)

### 2.1 Reduction to three moments — CORRECT

This is Newton's identities. The real parts cancel coordinate by coordinate because Re(rho^k) = Re(rho-bar^k).

### 2.2 Mass lemma — CORRECT

This is round-30 Lemma 1.1 at L=3, re-proved with the same argument-principle proof. I checked each step:
- deg Im P <= w - 4.
- Im P is not identically 0, because no root has its conjugate in the root multiset.
- The total argument change is -pi·sum(m).
- This gives |sum(m)| - 1 distinct real zeros.
- Hence 2·min(m+, m-) >= 3.

### 2.3 Quartic identity — CORRECT

y = Hu, xy = Hv, 3x²y - y³ = Hw. With Q = 3vv^T - (uw^T + wu^T)/2 we get h^T Q h = 3x²y² - y(3x²y - y³) = y^4. For z ⊥ u, z^T Q z = 3(v·z)², so Q restricted to u^⊥ has rank <= 1. This is a necessary condition, and only necessity is used.

**Gauge.** Positive dilation rho -> lambda·rho preserves clusters, so u_3 = 1 can be imposed when the gauge row h_g = e_3 has y_g > 0. Translation x -> x + c leaves u unchanged, so it does not interfere: the (a,b) system is 3 equations in 2 unknowns, consistent with expected dimension -1.

### 2.4 Proportional rows — CORRECT

h_j = lambda·h_k gives lambda > 0, x_j = x_k and lambda³ = lambda, hence rho_j = rho_k. (Remark: lambda = 1 or 1/2 or 2 would already be caught by the mass lemma, so the surviving proportional cases have lambda like 3/2.)

### 2.5 Four collinear rows — CORRECT

Let f(t) = Y^4 - h^T Q h, a quartic. Then f''' = 24(Y')³Y, and Y is positive on the hull of the four parameters. Four roots would contradict Rolle unless Y' = 0. Once Y is constant, the quadratic 3V² - YW - Y^4 vanishes at 4 points, so V' = 0 and all four roots coincide. Rows were taken distinct; coincident rows are the proportional case.

### 2.6 Gröbner certificates — CORRECT

- The interpolation is legitimate: the 8x6 design matrix has rank 6 in all remaining cases. I asserted this in my code for every gauge used.
- A unit ideal excludes every complex solution with y_g ≠ 0, and y_g > 0 holds on a realization.
- Positive-product certificate: B_jk = y_j² y_k² |rho_j - rho_k|² because (y_k h_j - y_j h_k) ⊥ u. So a product of y's and B's that lies in the ideal vanishes on any real solution, which is impossible on a realization.
- The positive control (`control.py`, three random rational root configurations) confirms the algebra: the compatibility polynomials and det(FQF^T) vanish at the true (a,b), and B_jk equals lambda^6 y_j² y_k² |rho_j - rho_k|² exactly.

### 2.7 Reproduction

The table below gives my certificates in the judge/Astra index order, for gauge rows g = 7, 0, 3, 5.

| certificate | Astra | mine (all four gauges) |
|---|---|---|
| proportional rows | 3, 5 | 3, 5 |
| four collinear rows | 8, 9, 11, 12, 13, 15, 16, 17, 19, 20 | the same 10 classes (sometimes different quadruples) |
| unit ideal or zero-remainder positive product | unit 0, 1, 4; positive 2, 6, 7, 10, 14, 18 | all of 0, 1, 2, 4, 6, 7, 10, 14, 18 excluded in every gauge |

- The split between unit-ideal and positive-product certificates depends on the gauge. That is expected, because a gauge discards the complex solutions with y_g = 0.
- With g = 5 my split is exactly Astra's: unit ideal for 0, 1, 4.
- **Every gauge excludes all 21 classes.**
- Astra's own `exact_plane.py` also runs (6 s) and prints the same certificates.

**Individual check, class 0.** With Astra's H, I verified that D·H = 0 and rank D = 5. My left-kernel compatibility polynomials are -2d(3d² - 5) and (3d^4 + 12d³ - 30d² - 20d - 45)/32 with d = a - b, exactly Astra's pair. Their values at d = 0 and at d² = 5/3 are -45 and -260/3. Class 0 is therefore excluded even without the det equation.

### 2.8 Size and equality — CORRECT

- **Rank <= 3 argument.** Each output coordinate on the 3-cube is an affine 0/1 function, hence constant, an input bit, or its complement. Each of the 3 input directions needs >= 4 changing coordinates, and these sets are disjoint, so >= 12 > 8 coordinates would be required. The argument is correct. Independently, my rank-<=3 BFS gives maximum **6**, which agrees with the round-30 search.md §2.2 table.
- **Rank 4.** Astra intersects affine spans of 5-row bases with the 256 binary rows. I used a different method: an orbit BFS restricted to sets of rank <= 4 (complete, since subsets of rank-<=4 sets have rank <= 4). Orbit counts:
  - 6 rows: 8 orbits;
  - 7 rows: 3 orbits;
  - 8 rows: 2 orbits;
  - 9 rows: 0 orbits.
- **Maximum 8.** Both 8-row orbits have rank exactly 4. In canonical form they are (15,51,60,85,170,195,204,240) and (15,51,85,105,150,170,204,240). Astra's three representatives fall into exactly these two orbits: her first two are equivalent.
- **The second orbit is not realisable.** Its H has four pairs of identical rows, so e_j - e_k ∈ V, which has mass (1,1) and gives a repeated root.
- **Uniqueness of the 8-point sign code is proved.** This upgrades round-30 search.md §2.2, where "every 8-point rank-4 code is half-cube" was "verified for the solutions found, not proved".
- **Attainment.** `cube.py` checks exactly that, with q = (-9/4, -2/3, 1/4, 8/3) and eps = (1,-1,-1,1):
  - the 16 numbers ±r_j and their conjugates are distinct;
  - all eight polynomials have the same norm N;
  - all pairwise differences have degree <= 4 (they are 3 or 4);
  - the sign code is the first canonical orbit.

## 3. What is proved, what is computed, and the status of the conclusion

**Proved by hand (checked):**
- the moment reduction;
- the mass lemma (L=3);
- the quartic identity and the rank-one restriction;
- the proportional-row lemma and the four-collinear-row lemma;
- positivity of B_jk;
- the gauge reduction;
- the rank-<=3 size bound;
- the cube attainment.

**Exact finite computations (each reproduced independently):**
- the census (72 orbits, 64 of rank 5);
- 43 rank-5 classes rejected by exact Lemma-1.1 witnesses;
- 21 classes excluded by the proportional, collinear, unit-ideal or positive-product certificates;
- the rank-<=4 maximum of 8, with two orbits, one of them killed.

**Conclusion.** "Every cluster at D=8, c=0 (squarefree N, over C(t) or Q(i)(t)) has rank <= 4 and size <= 8" **stands**. So does M(8,0) = 8 with a unique extremal sign code.

**Not covered:**
- non-squarefree N;
- F_q(t);
- the judge's other claims M(7,0) <= 5 and M(9,0) <= 5 (not in this block; still unrefereed);
- (R_4) at any other D.

## 4. Novelty and significance

**Novelty: new.** It resolves the round-30 §7 follow-up. The judge estimated the chance of finishing this in a round at 0.4. The finish works without homotopy software because Astra:
- eliminates (x, w) through the quartic identity, reducing each class to 3 equations in 2 unknowns;
- disposes of 12 of the 21 classes with two elementary geometric lemmas.

**Significance: peripheral for the goal.** It is a single-D rank theorem in the model. The record already treated this line as closed for internal work. It does not touch the integer problem: as Astra notes, no counterpart forces three exact moment equations at one real radius. The growth rate, the uniform bound, Questions A and B, and route 5 are all unchanged.

The only forward value is methodological. The proportional and collinear lemmas hold for any kernel dimension whenever three moments agree, so the pattern "saturated filter plus third-moment geometry plus small Gröbner bases" could be tried at D=12. That is optional, and it is low-value for the integer problem.
