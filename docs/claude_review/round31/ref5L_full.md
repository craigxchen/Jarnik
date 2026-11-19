# Referee report, round 31, Part V (Question B thread), block "full"

Scope: Part V of `astra_lemma_V/answer.md`, the §V.1 claims (full isotropy, multipliers, formula (V.1), the common
multiplier, two-sided width) and §V.2 (full isotropy with no size bound). It also covers the full-isotropy part of
the computed 191-digit certificate (§V.4 "COMPUTED", `v_isotropic_star.py`). Partial isotropy (§V.3 four- and
five-block families, §V.4 (V.2) and the P^(2-eps) rigidity) belongs to other blocks. I touch it only where §V.2
depends on it.

Scripts are in `ref5L_full/`. All are my own code except `astra_copy_v_isotropic_star.py`, a verbatim copy that I
reran.
* `gauss.py`: a Gaussian-integer toolkit (Euclid in Z[i]). Its `genuine(V)` computes the least-norm primitive tuple
  zeta_i = lambda V_i/conj(V_i) through a Gaussian gcd, with no factorisation.
* `test_V1.py`: randomised exact tests of every V.1 formula.
* `base_search.py`: my own integer specialisation of the V.3 four-star at a different parameter.
* `v2_full.py`: my own CRT/SL2 lift, shear and cycle flip on that parent.
* `v2_full_astra_parent.py`: the same code on Astra's 191-digit parent.

## 0. Verdict in five lines

1. Every full-isotropy claim in V.1 and V.2 is correct. I checked each proof line by line, and exact random tests
   confirmed every formula. The size-free construction I rebuilt independently on a new 193-digit parent, and again
   on Astra's 191-digit parent.
2. One overstatement needs correcting. The coordinator's summary says "the two-sided question is Question B itself".
   The two-sided fully isotropic configurations are *genuine* tuples, at height N^(1+o(1)) and Liouville-critical
   width H^(-3/10+o(1)), with a ten-cut (pair-cut) profile. A_5(delta) forbids them off Z for every delta < 1/10, so
   they are covered by A_5 with the same ratio point and the same Z. Excluding them is strictly weaker than Question
   B (delta < 1/30). (L_5°) does not apply, since the profile is not balanced over the 15 cuts.
3. The factorisation V_i = G_i U_i is already in the record, in the free thread's Part III (`astra_free_IV/answer.md`
   l.866-868, refereed in README_III). So is the range formula "primitive parent exponent = max - min of signed
   valuations" (`astra_free_III/iii_graph_notes.md` l.75, `iii_tri_section.md` l.22). The genuinely new items are:
   * the two-sided sandwich N/C^2 <= H <= NC^2 with an explicit circle;
   * c = 2 being equivalent to c = 1, and c^4 | prod t_ij;
   * the two-sided-width recovery c_i = N^(o(1));
   * the SL2/shear size-free construction.
   All are elementary.
4. Plainly, the size-free construction contradicts nothing because N is not the height. On my instance the genuine
   circle has H ~ 10^2104 against N ~ 10^193, and along the shear H grows like m^10 while the width falls like m^-2.
   So the width is about H^(-1/5), far above H^(-1/4). The genuine tuple is not even an A-point at delta = 0, and its
   points are about H^(0.3) apart, not on a sqrt-length arc. The violated hypothesis is |V_i|^2 ~ n_i, i.e.
   c_i = N^(o(1)). Here c_i ~ 10^382.
5. Significance for Question B: none beyond calibration. Full local isotropy, with or without size, is genuineness
   in disguise: with subpower multipliers it is the genuine problem at height ~N, and without size it is the genuine
   problem with the height decoupled from N. Kill criterion (i) applies to (a) and to (c). Local isotropy is not a
   new lever. My score for this block on the goal: 0.5/10.

## 1. Claim-by-claim

### 1.1 Propagation of isotropy along the marked line: CORRECT

Let p | J_e and let S_p be the marked side. For i, j in S_p, p | g_ij | det(V_i,V_j). Both columns are nonzero mod p
by primitivity, so they span one line mod p. If one marked column has p | |V_i|^2, the line is isotropic. Every
other marked column is then a nonzero multiple lambda V_i, so |V_j|^2 = lambda^2 |V_i|^2 = 0 and V_i·V_j = 0 mod p.
J_e is squarefree, so everything holds prime by prime for composite blocks. The same works in the balanced
orientation, since g_ij contains every block on which both i and j are marked.

### 1.2 Factorisation V_i = G_i U_i, |U_i|^2 = c_i, local multiplier rule: CORRECT (not new)

The isotropic lines mod a split p are exactly ker(Z[i] -> Z[i]/pi) and ker(Z[i] -> Z[i]/pibar). A primitive column
cannot be divisible by both pi and pibar, or p would divide it. Neither can it be divisible by an inert prime, or by
(1+i)^2 = 2i. This gives c_i = 2^b · (split primes), with b in {0,1}, single orientation per prime, and orientation
pi_p forced when i is in S_p.

The converse also holds. If pi_p divides V_i and V_j, then conj(V_i)V_j is divisible by pibar_p·pi_p = p, so p
divides both the minor and the dot product. My random test checks this.

The remark "|V_i|^2/n_i -> 1 forces c_i = 1 eventually" is trivial integrality. Precedent: free thread Part III,
"full-isotropy endpoint" and "with exact star norms the local isotropy is automatic"
(`round31/astra_free_IV/answer.md` l.866-868).

Test: `test_V1.py`, 300 random fully isotropic models. Each trial uses ten blocks of one or two random split primes
below 3000, random endpoint or flipped sides, and random U_i, including parent primes, extra primes, squares and
(1+i). The test checks primitivity, n_i | |V_i|^2, g_ij | det and g_ij | dot.

### 1.3 Formula (V.1) and N/C^2 <= H <= NC^2: CORRECT

H is the least norm of lambda in Q(i)* with lambda V_i/conj(V_i) in Z[i] for all i. At pi_p this requires
v_pi(lambda) >= -min_i s_ip, and at pibar_p it requires v_pibar(lambda) >= max_i s_ip. The prime (1+i) imposes
nothing, because (1+i)/(1-i) = i, and inert primes impose nothing either. Equality is attained, and the resulting
tuple has Gaussian gcd one.

Write s_ip = s(G_i) + s(U_i). Here s(G) has range 1 at p | N and range 0 elsewhere, while |s(U_i,p)| = v_p(c_i)
gives range at most 2 v_p(C). Range subadditivity in both directions then gives
p^(e_p - 2v_p(C)) <= p-part of H <= p^(e_p + 2v_p(C)), hence N/C^2 <= H <= NC^2.

For the explicit circle, B(G_i/conj G_i)(C/c_i)U_i^2 = (BC)·V_i/conj(V_i). It is integral, has norm N C^2, and has
the right doubled ratios.

Tests (`test_V1.py`, same 300 models): the formula (factorised) equals H computed independently by Gaussian gcd in
every case, the sandwich holds, and the explicit points are integral of norm NC^2 and proportional to
V_i/conj(V_i).

Gcd costs are fully accounted for; nothing is discarded. The ratio point is anchored, so no common phase is set to
zero.

### 1.4 Common multiplier: H = N prod_{p in D} p^(2v_p(c)); removal; c = 2 equivalent to c = 1: CORRECT

* At p dividing c but not N: s = ±v, so the range is 0 or 2v.
* At p dividing both c and N: marked rows have s = 1+v and unmarked rows have s = ±v, so the range is 1 or 1+2v.
* Single orientation everywhere gives U_i = eps_i·gamma, and gamma cancels from the doubled ratios.
* For c = 2: U_i = eps_i(1+i). Multiplying exact stars by (1+i) keeps them primitive (odd norms) and doubles every
  minor (det((1+i)a,(1+i)b) = 2 det(a,b)).

Tests: 400 random common-c models. The predicted H equals the gcd-computed H in all of them. In 200 exact-star
models, multiplying by (1+i) gives primitive c = 2 models with doubled minors and the same H = N.

Added remark: mixed c_i in {1,2} is exactly the unit-class phenomenon of the record. V_i = (1+i)G_i contributes the
unit i to the doubled direction. Bounded unequal multipliers of this kind are genuine tuples with arbitrary units,
which is why (V.1) needs no one-unit-class restriction.

### 1.5 c^4 | prod_{i<j} |t_ij|: CORRECT

Take an odd p with p^b || c. The five rows split into two orientation classes, which leaves at least four
same-orientation pairs (min over 5 = a + b of C(a,2) + C(b,2) is 4). For such a pair, p^b | conj(V_i)V_j.

* If p | g_ij, both rows are marked, so pi^(b+1) divides each and p^(b+1) divides the minor, while p || g_ij.
* If p does not divide g_ij, the factor p^b passes straight to t_ij.

If 2 | c, all ten minors are even and g_ij is odd, so 2^10 | prod t. Tests: all 400 common-c models satisfy the
divisibility. Consequence: |t| = N^(o(1)) forces c = N^(o(1)).

### 1.6 Two-sided width forces c_i = N^(o(1)): CORRECT

Hypotheses actually used:
* equal-weight blocks (each J_e = N^(1/10+o(1)));
* 1 <= |t_ij| <= N^(o(1));
* full isotropy, so c_i is an integer >= 1;
* the lower bound Delta >= N^(-3/10-o(1)) (two-sided).

Let W_i = V_i/sqrt(n_i). Then |[W_i,W_j]| = |t_ij| / sqrt(prod of the six separating blocks) = N^(-3/10+o(1)) for
every pair. The doubled chord is 2|sin angle| = 2|[W_i,W_j]|/(|W_i||W_j|), correct mod pi. For the pair attaining
Delta, the lower bound gives |W_i||W_j| = N^(o(1)), and each factor is >= 1. Write W_k = aW_i + bW_j with Cramer
coefficients that are ratios of minors, all N^(-3/10+o(1)), so |a|, |b| = N^(o(1)) and |W_k| = N^(o(1)).

Remark: if c_i = N^(o(1)) is given, the lower bound on Delta is automatic by Liouville (the sum argument of 1.7).
So "two-sided" matters only in the size-free direction.

### 1.7 Converse: c_i = N^(o(1)) and Delta <= N^(-3/10+o(1)) force equal weights: CORRECT

The chord bound gives S_ij >= (3/5 - o(1))W. Each pair cut separates 6 pairs, so sum S_ij = 6W, and every
S_ij = (3/5 + o(1))W. The vertex sum is sum_{j != i} S_ij = 4w_i + (2W - w_i) - 2w_i = w_i + 2W, which gives
w_i = 2W/5 and then w_ij = W/10. The separating-block identity n_i n_j / g_ij^2 holds in any orientation, so the
argument survives reversals.

### 1.8 "Covered by the Gaussian A_5 of Question B with the same ratio point and exceptional set": CORRECT WITH CORRECTIONS

What is right: with C = N^(o(1)), (V.1) gives a primitive genuine tuple zeta_i = lambda V_i/conj(V_i) of norm
H = N^(1+o(1)). Its points are distinct, because nonzero minors mean the doubled directions are distinct. Its
anchored ratio point y is literally the model's, so Z is the same. Its mutual distances are
sqrt(H)·Delta = H^(1/5+o(1)) = R''^(2/5+o(1)). So an essentially dense two-sided family would refute A_5(delta) for
every delta < 1/10, and with it V_E(5).

Corrected statement. *A fully isotropic model with c_i = N^(o(1)), |t_ij| = N^(o(1)) and width N^(-3/10+o(1)) is a
genuine primitive 5-tuple at height N^(1+o(1)). Its profile carries weight (1/10+o(1))W on each of the ten pair cuts
and o(W) on the five singleton cuts, and its width is H^(-3/10+o(1)), which is Liouville/Cilleruelo-Córdoba
critical. Such tuples are A-points for every delta <= 1/10 - o(1), so A_5(delta) for any delta < 1/10 (in particular
Question B's delta < 1/30) puts them in Z. The converse fails: excluding them off Z does not give Question B.*

Further notes:
* "The two-sided question is Question B itself" overstates this. It is a sub-question of the genuine problem,
  namely essential attainment of the Liouville exponent 10/3 on pair-cut profiles. Question B concerns 15/4 on
  balanced profiles.
* (L_5°) does not apply, because the profile is unbalanced: the five 1|4 cuts are empty.
* The "each block retains a prime outside C, index one" step is true but not needed for A_5, which has no index
  hypothesis.

### 1.9 V.2: full isotropy, fixed cofactors, arbitrarily small width off any Z: CORRECT (reproduced twice)

The proof checks out.
* SL2(Z) -> prod_p SL2(F_p) is onto via elementary matrices and CRT.
* The matrices gamma·A with gamma in Gamma(N) send every marked line to an isotropic line, which uses p = 1 mod 4.
  They preserve minors and primitivity.
* The coset is Zariski dense, and the preimage of the complement of Z is a nonempty Zariski-open set: it contains
  the identity because the V.3 start avoids Z.
* The shear (x,y) -> (x + mNy, y) lies in Gamma(N). f(y(m)) is a nonzero rational function of m, so finitely many m
  are bad, even when e lies in Z. The chords are of order m^(-2).

The cycle flip also checks out: the minors and the shared factors both scale by J_e/(a_i a_j), the doubled phases
change by the common factor gammabar/gamma, |V|^2/n is preserved, and the degree count is 4 - 2 + 3 = 5.

Dependencies:
* The escape from Z rests on V.3's four-star escape argument (another block).
* Primitivity after the flip needs, for p | J_e, a in e and k not in e, that p does not divide t_ak.
* **Minor discrepancy.** In Astra's 191-digit instance, 13 divides both t_25 and J_25 (marks 1-indexed; edge
  index 6, block factors 13·113·6898679189888741). The V.3 specialisation lemma claims blocks coprime to the
  cofactors and excludes primes <= 23. The discrepancy is harmless: only the cross condition above is used, and it
  holds. My own parent meets the full coprimality.

Reproduction 1 (`v2_full.py`, own parent). The base is the V.3 four-star with satellites (-9,-7,2,4), D = 6350400,
P = 4445280047 (Astra used 3175200023). The parent has 193 digits and 27 distinct split primes, factored by sympy
and certified prime. The base has isotropy on exactly 8 of 20 incidences (blocks 12, 13, 14, 15). My own CRT
lower/upper unipotent lift and shears m in {0, 1, 10^6, N, N^3} give:

| m | all 20 endpoint incidences | cofactors | log10 width | log10 H | log10 c_i |
|---|---|---|---|---|---|
| 0 | yes | unchanged | -396 | 2104 | 382 |
| 1 | yes | unchanged | -397 | 2107 | 383 |
| 10^6 | yes | unchanged | -408 | 2164 | 394 |
| N | yes | unchanged | -782 | 4034 | 768 |
| N^3 | yes | unchanged | -1554 | - | - |

Here N^(-3/10) = 10^-57.9, so the width is far below N^(-3/10). Dot products are also divisible by the shared
blocks. The cycle flip gives all 25 balanced incidences, the same t_ij, the same |V|^2/n, and the same genuine tuple
(equal H, proportional zeta).

Reproduction 2 (`v2_full_astra_parent.py`): my code on Astra's parent gives the same results (191 digits, width
10^-389 at m = 0, H = 10^2069.5). Astra's own script also reruns cleanly.

### 1.10 Why V.2 contradicts nothing

The hypothesis it violates is the archimedean one, |V_i|^2 ~ n_i, i.e. c_i = N^(o(1)). Here c_i ~ 10^382 to 10^768.
The quantity that matters for Question B is the genuine height H, and (V.1) says N/C^2 <= H <= NC^2 gives no
control once C is huge. Measured against H:
* Delta ~ H^(-0.19) at m = 0, tending to H^(-1/5) along the shear (log H grows by ~10 log m, log 1/Delta by
  ~2 log m).
* That is well inside Liouville (H^(-3/10)) and above H^(-1/4), so the genuine tuple is not even an A-point at
  delta = 0.
* The points are ~H^(0.3) apart, not on an arc of length C sqrt(R'').

Every fully isotropic size-free model is a genuine tuple (V.1 dictionary), and every genuine tuple is such a model
(c = 1). So "isotropy without size" is the genuine problem with height decoupled from N, and its parent-normalised
width carries no information. Astra's last remark in V.2 says the same, more softly.

### 1.11 Computed 191-digit certificate: CORRECT, but uninformative about the critical regime

The certificate reruns, and my code reproduces it. Two cautions.
* At P = 500D the "fixed" cofactors are ~10^43, larger than the blocks (~10^19). The base width is 10^-13.2 against
  N^(-3/10) = 10^-57. So the instance certifies identities and congruences, but it is nowhere near the asymptotic
  regime it is meant to illustrate.
* "Width < N^-3 after shearing" is the size-free statement, which is trivial by 1.10.

## 2. Novelty and significance

Novelty:
* Already in the record: the factorisation and the bounded-slack genuineness (free III, refereed), and the range
  formula for the primitive parent (free III notes).
* New and correct: the sandwich N/C^2 <= H <= NC^2 with an explicit circle; the common-c formula; c = 2 equivalent
  to c = 1; c^4 | prod t; the two-sided recovery c_i = N^(o(1)); the converse equal-weight forcing; the Gamma(N)-coset
  size-free construction. All are short, elementary arguments.

Significance for Question B: none. The block shows that full local isotropy adds nothing independent of
genuineness. With subpower multipliers the models are genuine tuples at height N^(1+o(1)), on the critical
pair-cut profile, and are covered by A_5. Without size, they are genuine tuples whose height is not N, and the
parent-normalised width is meaningless. Kill criterion (i) applies: stop the full-isotropy line. Only partial
isotropy (other blocks) could still separate models from genuine tuples, and V.4's (V.2) transfer already shows how
quickly partial isotropy turns into genuine data.

## 3. Standing cautions

* No uniform exceptional set from fixed-S Subspace was assumed anywhere.
* Gaussian gcd and denominator costs are exact in (V.1) and checked by an independent gcd computation.
* The ratio point is anchored, and no drifting phase was set to zero.
* The shape (hence the cofactors) is chosen after Z in V.2/V.3. That template dependence is acknowledged by Astra
  and is harmless for an existence statement.
