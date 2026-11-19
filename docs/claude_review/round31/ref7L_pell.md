# Referee report, round 31, block "pell": Part VII.5 and the specialisation lemma of VII.7

Source: `astra_lemma_VII/answer.md`, lines 1662-1693 (VII.5) and 1723 (specialisation with slack).
My scripts (my own exact Q(i)[t] arithmetic, Fractions only; sympy used only for two resultants and a rank) are in
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref7L_pell/`:

| script | what it does |
|---|---|
| `pell_core.py` | Gaussian-rational polynomials: mult, divmod, gcd, conjugation, 2x2 matrix powers; `analyse` (equal norms, gcd, D, contact, distinctness) |
| `run_a.py` | (a) five marks n=1..5 and four marks n=1..5, three seed sets each; three marks d=1,2,3,4,6,9 |
| `run_b.py` | exact identities, raw contact without normalising leading coefficients, special seeds s_j+s_k=0 |
| `invariants.py` | the two n-independent bilinear invariants (below), n=1..6 |
| `merge.py`, `shared.py`, `shared2.py` | which root coincidences could lower the primitive degree; resultants for real seeds (n=1,2,3) |
| `search_b.py` | exhaustive exact search over multiplicative Pell-block templates (all 5-subsets of each charge class) |
| `other_A.py` | A = [[a^2+1,a],[a,1]] for a = t, t^2, t+1, t^2+t-3, t^3, and mixed exponents n_l |
| `spec_e.py` | (e)/(f) integer specialisation, Gaussian gcd, Bezout identity, arc exponent |

I read Astra's `vii_pell_check.py` only afterwards. It uses integer seeds (1,3,5,8) and n = 1,2,3,5, and gives the same
(D,h) values as mine.

## Verdicts

| claim | verdict |
|---|---|
| (a) five marks: equal norms, gcd 1, distinct, D = 8n+1, h = 4n | **CORRECT** (reproduced n=1..5, 3 seed sets; strengthened: true for **all** distinct real seeds) |
| (a) four marks D = 6n+1, h = 4n; three marks D = 3d, h = 2d | **CORRECT** (reproduced) |
| (b) asymptotics v, g^2 = g_0, g_j/g = 1 + 2 i s_j t^(-4n) + ..., role of beta | **CORRECT**. In fact everything is exact (two invariants below) |
| (b) "h = 4n when at least one edge has s_j+s_k != 0" | **CORRECT WITH CORRECTIONS**: h = 4n for **every** choice of distinct seeds. The condition on s_j+s_k is irrelevant, because the pair (z_12, z_13) always has contact exactly 4n |
| (b) attains the covering maximum h = floor((2D-1)/4) | **CORRECT** |
| (b) can the contact be pushed to 2h > D at five marks? | **No, in every family I tried.** Best 2h-D is -1 (dense) and 0 (non-dense) |
| (c) Zariski density for each fixed n | **CORRECT** |
| (d) h/D <= 2/3, 2/3, 3/5; h = 2D/3 impossible at four marks | **CORRECT** (checked line by line) |
| (e) integer specialisation, exponent 1/2 + 1/(2D) | **CORRECT** (computed). Consequence over Q: alpha_ess(e) <= 4 (stated precisely in §5) |
| (f) specialisation lemma (Bezout bounds the gcd) | **CORRECT** |

Nothing in this block is false. All claims are elementary, and I verified them exactly.

## 1. (a) Exact reproduction

Setup: A = [[t^2+1, t],[t,1]], beta = t+2i, G_j = X_j + iY_j with (X_j,Y_j)^T = A^n(1,s_j)^T, z_0 = beta prod conj(G_l),
z_jk = conj(beta) G_j G_k prod_{l != j,k} conj(G_l). The seed sets were (1/2, 1/3, 2/7, -5/3), (3, -1/5, 7/4, 2/9) and
(1,2,3,5).

* **Five marks** (r = 4; edges 12, 13, 23, 14), n = 1..5, all three seed sets: the norms are equal and the gcd over
  Q(i)[t] has degree 0. The marks are distinct and D = 8n+1 = 9, 17, 25, 33, 41. The contact is h = 4n, and **all ten
  pair contacts equal 4n**. All five marks are monic, so the raw contact (no leading-coefficient normalisation) is also
  4n.
* **Four marks** (r = 3; triangle), n = 1..5: D = 6n+1 and h = 4n, with gcd 1, equal norms and distinct marks.
* **Three marks** G_j = t^d - a_j + i, d = 1,2,3,4,6,9, a = (0,1,2) and (1/2,-3,7/5): D = 3d and h = 2d, with gcd 1.

**Strengthening (mine): primitivity for every choice of distinct real seeds.** Two exact identities hold for all n
(verified n = 1..6 in `invariants.py`; the proof is det A = 1 plus A-invariance of a quadratic form):

    G_s conj(G_u) - conj(G_s) G_u = 2i(s-u)
    G_s G_u conj(beta) - conj(G_s) conj(G_u) beta = 4i((s+u)t/2 - 1 + su)

`merge.py` lists the only root coincidences that could lower D. They are G_j with conj(G_k), G_j with G_k, a block with
beta, and a block with its own conjugate.
* G_j with G_k, and G_j with conj(G_j), would make A^n(alpha) singular.
* G_s(-2i) = ±(1+is) != 0, by the second identity with u = s at t = -2i.
* A common root alpha of G_s and conj(G_u) satisfies (s+u)alpha/2 - 1 + su = 0 by the second identity, so alpha is real
  or s = -u with su = 1, which is impossible. But G_s has no real root, since X^2 + Y^2 > 0 there.

Hence D = 8n+1 exactly for every pairwise distinct real seed vector, not only for generic seeds. The resultants in
`shared2.py` (n = 1,2,3) confirm that no real (s,u) works. The first identity also gives z_12 - z_13 =
2i(s_2-s_3)·conj(beta) G_1 conj(G_4), of degree 4n+1. **So h = 4n exactly whenever s_2 != s_3**, with no condition on
s_j + s_k.

## 2. (b) Derivation, role of beta, and attempts to push the contact

**Derivation: CORRECT.** v = 1/(t+v), so v^2 + tv - 1 = 0 and 1 - v^2 = tv. The eigenvalue is
lambda = t^2+1+tv = t^2+2+O(t^-2). Write (1,s) = a(1,v) + b(-v,1). Since -v+i = i(1+iv),

    G_j = (1+iv)(a lambda^n + i b lambda^-n),
    g_j/g = (1+i rho_j)/(1-i rho_j),   with rho_j = lambda^(-2n)(s_j - v)/(1 + s_j v).

This gives 1 + 2i s_j t^(-4n) + O(t^(-4n-1)), as claimed. g^2 = (1-v^2+2iv)/(1-v^2-2iv) = (t+2i)/(t-2i) = beta/conj(beta)
exactly. So beta cancels the moving eigenphase **exactly** (Astra's "including all its lower Laurent terms" is right).
For symmetric A = [[p,q],[q,r]] in SL_2(Q[t]), the same computation gives
beta = ((p-r) + 2iq)/gcd_R(p-r,q). deg beta >= 1 for nonconstant A: if beta were constant, v and then lambda would be
constant.

**Pushing the contact.** In every case, 2h - D = -deg beta at the dense optimum.
* *Other A* (`other_A.py`): A = [[a^2+1,a],[a,1]] gives 2h - D = -deg a, namely -1, -2, -1, -2, -3 for a = t, t^2, t+1,
  t^2+t-3, t^3. So a = t (or any linear a) is optimal.
* *Mixed exponents n_l* only lower h to 4 min n_l while D grows (e.g. n = (2,1,1,1): D = 11, h = 4).
* *Special seeds* s_j + s_k = 0 raise only the anchor pair: n=1 gives pair (0,1) contact 5, others 4. Seeds that make
  two marks share a block (s_3 = s_4, or s_4 = s_2) are not distinct, and the gcd then drops (D = 7 at n=1, marks
  coincide).
* *Exhaustive exact template search* (`search_b.py`): marks = beta^e0 conj(beta)^(w0-e0) prod G_l^e_l conj(G_l)^(w_l-e_l),
  all 5-subsets of each phase-charge class. Templates were r = 4 (w = 1111, w0 = 1, 2), r = 5 (w0 = 1), r = 3 with
  w = 111, 211, 221, and r = 4 with w = 2111.
  * Best dense (rank 4): **2h - D = -1** (Astra's family, D = 9, h = 4).
  * Best non-dense: **2h - D = 0**: five edges of K_4, no beta, D = 8n, h = 4n, rank 3. So **h = D/2 exactly is reached**
    at five marks, but only on a 3-dimensional family.
  * No 5-subset anywhere has 2h > D.

**Why (proof sketch, mine).** For generic seeds, a pair's contact is at most 4 min n_l over the blocks that separate it.
The primitive degree is D = sum_l 2n_l w_l + w_0. 2h > D therefore needs W = sum w_l <= 3. The phase-charge condition
sum e_l + 2e_0 = const then admits at most half of prod(w_l+1) <= 8 vectors, so at most 4 marks. Special seeds cannot
raise all ten contacts: (12,13)-type pairs are pinned by the first identity.

**Answer to (b).** The best five-mark value found is h - D/2 = -1/2 for dense families and 0 for non-dense ones. No
five-mark example with 2h > D exists in these families. This agrees with VII.2 (4h+1 <= 2D for dominating families) and
shows the slack of 1 there is attained.

## 3. (c) Zariski density: CORRECT

At t = 1, A(1)^n is a real invertible matrix, so s -> (X:Y) is a bijection of RP^1 onto RP^1 minus a point. The phase
map (X:Y) -> (X+iY)/(X-iY) is a bijection RP^1 -> S^1, so rational seeds give a dense set of g_j(1) in S^1. The edge
map (g_1..g_4) -> (g_1g_2, g_1g_3, g_2g_3, g_1g_4) has integer matrix of determinant -2 (triangle: -2; appended row
(1,0,0,1) is block-triangular). So it is an isogeny, the image of a dense subset of (S^1)^4 is dense, and the division by
the fixed g_0(1) is a translation. (S^1)^4 is Zariski dense in the complex torus.

The precise sense: **for each fixed n, the union over rational seeds of the ratio curves is Zariski dense in T, and
already the set of their values at t=1 is.** No single curve is dense (each is a curve). This is exactly what the
specialisation lemma needs: for any proper closed Z, some curve of the family is not contained in Z. The same holds for
the four-mark family, and for the three-mark family in T_3.

## 4. (d) Pair-degree bounds and the four-mark endpoint: CORRECT

* For a root pair (rho, conj rho) of N with multiplicity m, the mark multiplicities are a_i at rho and m - a_i at
  conj rho. Primitivity gives min a_i = 0 and max a_i = m, so D = sum of thresholds.
* deg(q_i/q_j) = sum over root pairs of |a_i - a_j| = the number of thresholds separating i from j. No cancellation is
  possible, since rho is not real.
* Contact h <= deg(q_i/q_j) for any nonconstant ratio.
* A threshold of type a|(k-a) lies in a(k-a) pairs. Hence C(k,2) h <= max_a a(k-a) · D, i.e. 3h <= 2D, 6h <= 4D,
  10h <= 6D, which gives 2/3, 2/3, 3/5. 3/5 is delta = 1/10, Cilleruelo-Cordoba.

**Four-mark equality.**
* Equality forces only 2|2 thresholds and all six pair degrees equal h. Two distinct 2|2 cuts cannot be nested, so each
  root pair carries one type.
* Then d_1 = d_2 = d_3 = D/3, and q_1 = g_2g_3, q_2 = g_1g_3, q_3 = g_1g_2. The leading constants are ±1 with a common
  sign, which can be absorbed.
* So g_1^2 = q_2q_3/q_1 = 1 + O(t^-h). In characteristic 0 the square root with g_1(∞) = 1 is 1 + O(t^-h).
* g_1 is a product of (t-rho)/(t-conj rho) powers with no cancellation, so deg g_1 = D/3 < h = 2D/3, a contradiction.

The argument needs characteristic != 2 (the square root), as Astra says. Nonsquarefree parents are covered because
multiplicities are counted by thresholds.

## 5. (e) Integer specialisation (`spec_e.py`): CORRECT, and what it means

Run on the five-mark family: n = 1,2,3; seeds (1,2,3,5) and (1/2,1/3,2/7,-5/3), with the common denominator cleared;
T = 10, 10^2, 10^3, 10^4, 10^6, 10^9.
* The Gaussian gcd of the five values has norm in {16, 32, 64, 272, 416, 656, 928, 102992}. Each divides the integer
  Bezout constant computed from an exact identity sum a_i z_i = 1 over Q(i)[t] (asserted for every T).
* All five values have the same primitive norm H = R^2.
* log(arc)/log R decreases to 1/2 + 1/(2D). Examples:
  * n=3: 0.5216 at T = 10^9 (target 0.52).
  * n=1: 0.5648 at T = 10^9 (target 0.5556); convergence is O(1/log T).
* arc/sqrt(R·T) stays bounded (2.2 to 22.7), so arc ≍ R^(1/2) T^(1/2) = R^(1/2+1/(2D)).

**Precise consequence over Q (mine, via the specialisation lemma at delta < 0).**
* Fix n. The integer points of the curves lie off any prescribed proper closed Z ⊂ T_5 for all large T, on some curve.
* They satisfy dist(y,e) ≍ H^(-(1/4 - 1/(4D))). In McKinnon-Roth form, with H = H_{-K} = R''^2, this is approximation
  exponent alpha = W/h = 2D/h = 4 + 1/(2n). This equals d/m for a rational curve of -K-degree 2D whose branch at e has
  multiplicity h.
* Letting n -> ∞: **alpha_ess,e(-K) <= 4 unconditionally on X = twisted LM_5**. Equivalently, for every eps > 0 the
  primitive five-point tuples on arcs of length <= R^(1/2+eps) are Zariski dense in T_5.
* So Vojta's prediction alpha_ess >= 4 (= A_5(delta) for all delta > 0), if true, is sharp, and the record's "open
  window (15/4,4)" is really (15/4, 4] with 4 the maximum possible value.
* It does **not** give exponent exactly 4 densely: each fixed-n family is at 4 + 1/(2n). It says nothing about A_5(delta)
  for delta > 0 or about (L_5°). The family has 1|4 cuts (beta, block 4) and only 5 of the 15 cuts, so it is far from
  balanced.

**Novelty against the record.**
* grep finds no D = 8n+1 or beta = t+2i family.
* The golden 8-point family reaches exponent exactly 4 (arcs C√R, Delta H^(1/4) = 3.70246) but **on a curve**
  (`ref6L_subspace.md`).
* That report says "Whether a dense unit family reaches exponent 4 is open". The Pell result does not settle that
  (exactly 4), but it gives the first record statement that **4 is the infimum of exponents attained densely**, i.e.
  alpha_ess <= 4.
* This is the expected outcome: the heuristic count says delta = 0 is counting-critical, and covering families of
  rational curves are the standard upper-bound mechanism for alpha_ess (McKinnon's conjecture). I know of no
  literature statement for this variety. I did not check the Cilleruelo-Granville paper for a dense five-point
  construction at R^(1/2+eps); it is plausibly known to experts in some form.
* The four-mark function-field family is the polynomial analogue of the record's integer Cilleruelo-Granville Pell
  families (expert note §3, A_4 false for delta < 1/6). It adds the impossibility of equality at four marks (new, easy).
* The three-mark equality family is the classical Jarnik-sharpness configuration (arcs ≍ R^(1/3)). Its density in T_3 is
  new to the record but trivial.

## 6. (f) Specialisation lemma (VII.7): CORRECT

* Suppose z_i ∈ Q(i)[t] are primitive with no common factor. Then sum a_i z_i = 1 with a_i ∈ Q(i)[t]. Multiplying by
  the coefficient denominators M_a and M_z gives sum (M_a a_i)(M_z z_i) = M_a M_z in Z[i][t]. So the Gaussian gcd of the
  M_z z_i(T) divides the fixed M_a M_z (checked numerically above).
* Hence R'' ≍ T^D, width ≪ T^(D-h), and arc ≤ R''^(1 - h/D + o(1)) ≤ R''^(1/2 - delta) for h/D ≥ 1/2 + delta_0 > 1/2 + delta.
* The ratio points stay off Z for large T, because some equation of Z pulls back to a nonzero rational function of t.
  Distinctness holds for large T.
* The curve-dependent constants are absorbed by strict slack, as she says.
* Quantifier check:
  * The lemma needs one curve escaping each given Z, with a uniform bound h/D ≥ 1/2 + delta_0. It does not need one
    curve escaping all Z.
  * Its hypothesis needs curves over Q(i)[t] with the conjugate structure, so that z_i(T) conj(z_i(T)) = N(T) at real T.
  * No standing caution is violated: the gcd cost is bounded by Bezout, and the height is the primitive one.
* As she says, VII.5 approaches 1/2 from below and gives no delta > 0.

## 7. Standing cautions

* Primitive height: used throughout (gcd computed and removed; integer gcd bounded by Bezout).
* No fixed-S Subspace input appears in this block.
* No drifting phase is replaced by zero: beta cancels the phase g^2 exactly, not approximately.
* The exceptional-set quantifier is respected: density is per fixed n, and the specialisation lemma needs only "for
  each Z some curve".

## 8. Significance (blunt)

VII.5 is a correct and clean **sharpness/calibration** result.
* It shows the covering-family inequality 4h+1 <= W of VII.2 is attained with slack exactly 1.
* It shows that dominating families can never reach Question B's side: they stop at h = (D-1)/2.
* Through integer specialisation it gives the new, modest, unconditional fact alpha_ess(e) <= 4 over Q.

It provides no evidence either way on A_5(delta) for delta > 0 (the actual Question B), on (L_5°), or on the growth
rate or M(C). The only open direction it suggests is a non-dominating (isolated or low-dimensional) mechanism. Here the
best I found is h = D/2 exactly on a 3-dimensional family (five edges of K_4, no beta), still not 2h > D.
