# Round 31, Part V -- judge report: an unbounded circle cluster over the fixed field F_2

Source: Part V of `astra_free_V/answer.md` (heading "# Part V", lines 1138-1256; script `v_fixed2_check.py`;
coordinator's question `reply5.txt`, Astra's summary `final5.txt`). Five referees/analysts each took one block:

| key | report | scope |
|---|---|---|
| construct | `ref5_construct.md` | V.1: subspace polynomials, common norm, (V.1), (V.2), legitimacy of the instance |
| lower | `ref5_lower.md` | V.2: (V.3), (V.4), the Omega(s^2 log s) combination, the d=4 enumeration |
| reconcile | `ref5_reconcile.md` | Part V against the record: (IV.2), the transferred 3M log M bound, A_5/(L_5°)/(L_8), F3+ |
| extend | `ref5_extend.md` | analyst: other V/T/designs, exact M_2(D,0) and M_3(D,0), odd characteristic |
| lit | `ref5_lit.md` | literature, prior art of the mechanism, what Part V means |

Judge scripts are in `judge_V/` (own F_2[t]/F_4[t] code, integers as bit vectors; nothing imported from Astra or the
referees): `f4.py` (arithmetic), `checkV.py` (V.1 for m=4,8, conjugation, (V.1) identity, degree formula, minimality
of L, separability of ratios), `twist.py` (lower's unit-twist identity, d=3..6), `s8.py` (lower's 8-point cluster of
point degree 104), `m2.py` (exact M_2(D,0), D<=9), `lemmaP.py` (Lemma P over Z). Outputs `*.out`.
**Every claim labelled [J] is the judge's own and is unrefereed.**

Labels as in Parts III-IV. "Proved" means a complete proof checked line by line by a referee, or an exact
computation reproduced with independent code. **No referee found an error in any claim Astra made in the form she made
it.** The GAP and FALSE verdicts in the summaries are of three kinds: scope remarks on V.2 (lower: thresholds; the
unit-twisted case, repaired), questions the referees posed themselves (reconcile (d): "Frobenius/inseparable failure",
FALSE; lit (a3): "isotriviality is the reason", FALSE), and an analyst's own extension (extend: odd characteristic,
GAP = nothing constructed). I rule on each in §1.

## 0. Verdict in five lines

1. **Theorem (now in the record's function-field model): for the anisotropic norm A^2+AB+B^2 of F_4(t)/F_2(t),
   sup_D M_2(D,0) = infinity.** For every r>=2, m=2^r, and every L >= m(r-1), there are m monic points of degree
   D=(m-1)L on one parent of degree 2D with max pair difference D/2 + m(r-1)/2 - L/2. Rate m ~ sqrt(2D/log_2 D), along
   D in (m-1)Z only. Five independent reproductions plus the judge's (m=4, 8) [J]. N is not squarefree.
2. **V.2 (Omega(s^2 log s) inside the construction) is correct** but is a statement about Astra's own family with
   thresholds d>=48, s>=6. It is superseded by the analyst's **Theorem E** (all additive V containing 1, all T, unit
   twists, any common division, no thresholds), which I checked line by line [J]. Neither bounds M_2(D,0).
3. **Nothing better was found:** no o(m^2) parent, no growth scale, no odd-characteristic construction. Exhaustive data:
   M_2(D,0)=4 for D=11..17, M_3(D,0)<=4 for D<=11 (judge reproduced M_2 for D<=9 exactly, counts included [J]).
4. **Filter F3+ (precise in §3(d)):** no argument valid verbatim over F_2[t] with this norm (fixed small residue
   fields, exact norm, anisotropy, all N) can prove the uniform bound. It does not touch growth-rate arguments, and
   number-field Vojta (the record's conditional routes) is not refuted.
5. **For the integer problem: nothing changes, and this is not evidence against the uniform conjecture over Z** (a
   characteristic-2 additive mechanism with no integer counterpart). Score for Part V on the goal: **2/10**.

## 1. Referee-corrected statements per block

### 1.1 construct: V.1 (referee: correct; statement strengthened)

**Proved**
* *Subspace polynomials.* For a binary subspace H of F_2[t], P_H(X)=prod_{v in H}(X+v) is additive, and
  P_{H+<b>}(X)=P_H(X)^2+P_H(b)P_H(X). Conjugation theta -> theta+1 fixes v in F_2[t] and sends Y=t^L+theta to Y+1, so
  conj(z_H) = prod_{v in (1+H)\{1}}(Y+v) = prod_{v in (V\H)\{1}}(Y+v). [J] Checked exactly in `checkV.py` for every H,
  m=4, 8.
* *Common norm.* N = prod_{v in V\{0,1}}(Y+v) lies in F_2[t], deg N=(Q-2)L=2D, D=(m-1)L. Equivalently
  N=prod over even v != 0 of (t^L+v)^2+(t^L+v)+1. [J] reproduced (N computed separately as a product, all z_H
  have norm N, all distinct, all monic of degree D).
* *(V.1).* z_H - z_H' = P_J(1)·P_J(Y)/Y with J=H cap H', and deg(z_H - z_H') = (Q/4-1)L + sum_{v in J} deg(1+v).
  [J] the identity and the degree formula hold exactly on all 6 (m=4) and 28 (m=8) pairs.
* *(V.2).* max_J sum_{v in J} deg(1+v) = S_max = (d-2)2^{d-2}, attained at J=t^2 F_2[t]_{<d-2}; hence the least L
  for T=t^L is L_m=(d-2)2^{d-1}=m(log_2 m-1), and some pair attains D/2.
* **Corrected/strengthened statement (construct + judge):** *For r>=2, m=2^r, d=r+1 and every integer L>=d, the m
  points z_H have max_{H != H'} deg(z_H - z_H') = D/2 + S_max - L/2, with D=(m-1)L. Hence they form a c=0 cluster iff
  L >= 2S_max = m(r-1), and M_2(D,0) >= m for every D in (m-1)Z with D >= m(m-1)(r-1) (odd L included, where the
  condition reads deg <= floor(D/2)).* [J] checked at L = L_m-1, L_m, L_m+1, L_m+5 for m=4, 8: the excess equals
  S_max - L/2 exactly, and L_m - 1 fails.
* *Legitimacy.* The form is anisotropic over F_2(t) and at infinity (a field norm; the infinite place is inert). All
  points monic (units of F_4^* have norm 1, so this is no restriction). Common gcd degrees 0, 1, 0, 8, 30 for m=4..64;
  dividing by it tightens the cluster (c' = -deg g/2), so **no gcd cost is discarded**. Every prime of N has even
  degree and splits; none is inert. No two points are unit multiples; conjugates are not in the cluster.

**Computed.** (L, D, maxdiff) = (4,12,6), (16,112,56), (48,720,360), (128,3968,1984), (320,20160,10080) for
m=4..64; every norm equal, every point distinct, maxdiff = D/2 exactly. Reproduced independently by construct,
reconcile, extend, lit (m<=16) and the judge (m=4, 8) [J].

**Corrections to Astra's wording.**
1. "N ... This uses the model allowing prime powers." Correct, but the record should say explicitly that **N is not
   squarefree** (lit: deg gcd(N,N') = 8 at m=4, 36 at m=8; [J] confirmed N not squarefree at m=4, 8). The record's
   model (round 29 `ff.md` §1, exponents e_k) allows this, as does the integer conjecture. The squarefree variant over
   F_2 is **open**.
2. The asymptotic m ~ sqrt(2D/log_2 D) holds along D in (m-1)Z, D >= m(m-1)(log_2 m-1); no monotonicity of M_2(·,0)
   in D is proved, so for general D only "limsup" statements follow.
3. The form is the Eisenstein-type norm of F_4/F_2, not X^2+Y^2 (degenerate in characteristic 2), and the extension
   is unramified everywhere, unlike Q(i)/Q.

### 1.2 lower: V.2 (referee: correct; scope corrected; one missing case repaired)

**Proved**
* Setup: all raw marks F_i=prod_{v in S_i}(Y+v) (Y=T+theta, S_i an affine side of a form with l(1)=1) have norm
  P_V(Y); division by a common factor of degree a gives c=0 iff a >= 2deg(F_i-F_j)-D_0 for all pairs, so the full gcd
  is optimal; h=|cap S_i| <= Q/(2s).
* Raw difference F_i - F_j = P_J(1)·P_J(Y+a), a in S_i cap S_j.
* (V.3) g >= Q(d/2-6) (lossy but valid: exact minima deg P_J(1) = (Q/4)(d-3)+1 for d=3..9; deficit sum max 2Q-2).
* (V.4) g <= hB + Qd/4 (even-degree contraction, affine fibres v_*+(V cap P^j F_2[t]), at most half in one side,
  e·pi_2(e) <= 2^e, sum 1/(4^k-1) = 0.4207 <= 4/9). Every step checked.
* Combination, **corrected statement:** *Fix V=F_2[t]_{<d}. For d >= 48 and s >= 6 (s >= 5 suffices), any c=0
  cluster obtained from s raw marks of this family (any T, any sides, unit twists allowed) by dividing by a common
  polynomial has point degree >= Qsd/16, hence parent degree >= s^2(1+log_2 s)/4.* Nothing is asserted for
  Q < 2^48, i.e. the Omega statement concerns only s >= 2^47 in practice.
* **Missing case, repaired (lower).** For T in V the raw marks have leading coefficient theta or theta^2; monic
  normalisation gives twisted differences not covered by (V.3). Identity: if 0 in S_i, 1 in S_j, then
  F_i/theta - F_j/theta^2 = A·O, A = prod_{v in S_i cap S_j}(theta+v), O = sum of the odd-index coefficients of the
  linearised polynomial P_J(X)=sum c_k X^{2^k}, with deg O = deg P_J(1) - a_1, a_1 in {1,2}. So twisted (V.3) reads
  g >= Q(d/2-6)-4 and the conclusion stands. **[J] Ruling: GAP in Astra's V.2, repaired; the repair is correct.** I
  redid the decisive computation (`twist.py`): the identity holds exactly on all 12/56/240/992 ordered pairs for
  d=3,4,5,6, and deg P_J(1) - deg O takes only the values 1, 2.
* Opposite sides of one form are also covered (h=0 gives no cluster for d > 24).

**Computed**
* d=4: 26,176 cases, max residual gcd degree 4 (Astra's number), histogram {0:1574, 1:1296, 2:3072, 3:216, 4:386},
  two gcd methods; extended to 64 classes of T: still 4.
* d=5, **exhaustive**: 688,747,008 cases, max residual 19 against the unsimplified bound 22 and Qd/4 = 40.
* Exact family minima at d=4, s=2..8: point degree 10, 32, 32, 80, 80, 96, **104**. The s=8 entry beats V.1 (112):
  Y=t^15+theta, mixed sides, h=1, gcd degree 16, max reduced difference 52 = D/2. **[J] reproduced** (`s8.py`: over
  all 256 side choices at T=t^15 the minimum reduced degree is 104, sides (0,0,1,1,0,0,1,1), gcd degree 16, max
  reduced difference 52). So V.1 is not optimal inside its own family; the order is unaffected.
* Small parents (d=3, L=4..8, 16-144 monic points): largest c=0 cluster 4, all family points.

**Gap (scope; not an error).** V.2 is not a lower bound for M_2(D,0), for this parent, or for this norm. It does not
cover: other additive V (but see Theorem E, §1.4); conj(Y)=Y+B with B non-constant; non-linear transversals
(2^{Q/2} points with norm P_V(Y)); independent choices of f or fbar at each split prime; non-common division;
prime powers and extra factors; larger constant fields; d < 48 (Theorem E removes this last one).

### 1.3 reconcile: Part V against the record (referee: correct; two corrections)

**Proved**
* *(IV.2) fails in characteristic 2 at a fixed field.* Through the star dictionary (U = product of big parts and one
  prime from each small conjugate pair; P_i = gcd(z_i,U); E_ij = 0 since g_ij = G·Gbar divides D_ij), the Part V
  cluster satisfies every hypothesis of (IV.2) with delta=0 (all pairs at m=8, 16, 32; 300 sampled pairs at m=64), and
  (IV.2) fails: LHS = 62016 vs RHS 0 at m=32, 1375616 vs 0 at m=64; the violation persists after removing the global
  gcd and at odd D (L=129). Every step of Part IV's proof except Langer's Theorem 0.1 holds in any characteristic, so
  the blown-up rational surface over F̄_2 has 3e_orb - (K+B)^2 <= -(2/9)(LHS-RHS) < 0. Not redone by the judge;
  the m=32/64 numbers are the referee's.
* *Transferred growth bound.* The record's `two_thirds` argument transfers to F_2[t] with this norm (odd-degree
  primes inert, fibre 2^e+1; even-degree split, fibre 2^e-1; no ramified prime; Hensel by the nondegenerate polar
  form): **W = deg N >= 3M log_2 M - O(M) - 4(M-1)c, i.e. M <= (1/3+o(1)) W/log_2 W**. Exact W_min = 223, 674, 1611,
  4210 at M = 32, 64, 128, 256. Part V sits above it by a factor of about (2/3)m (35.6x at m=32). The exact cut
  identity holds on the cluster (30960 = 30960 at m=32; 318432 = 318432 at m=64).
* *Question B analogues fail in this model.* Independent five-subtuples have 1/30 - 8/(15m) <= delta_S <= 1/30,
  near-exact balance, H_Sigma/w <= 40/m; J4's averaging + Schwartz-Zippel holds over F_4(t). Hence the F_2(t)
  analogues of **A_5(delta) for every delta < 1/30, (L_5°), and V_E(5) are false**; the Liouville exponent 15/4 is
  attained Zariski-densely.
* *Separability (reconcile (d)).* The claim "Part V is an instance of the inseparable/Frobenius failure of Vojta" is
  **FALSE**: all tested ratios and cross-ratios have nonzero derivative. **[J] ruling: concur.** `checkV.py`
  confirms that every ratio z_a/z_b has nonzero derivative (all pairs, m=4, 8). The mechanism is additive
  (Artin-Schreier: Galois acts on Y as translation by 1).

**Corrected**
1. *(L_8) analogue fails* only asymptotically: eight independent marks need m >= 128, and failure needs
   m >> 896/c; the identification of H_Sigma with reduced-chord degrees (ptolemy Thm B) was not re-derived in
   characteristic 2, and m=128 was not computed. Status: **plausible, not a theorem of the record**.
2. *F3+* is correct only for uniform (or sub-sqrt) statements; see §3(d). [J] Minor slip in the report: "m ≈
   sqrt(W/(2 log_2 W))" should read m ≈ sqrt(W/log_2 W) (W = 2m(m-1)(log_2 m-1) and log_2 W ≈ 2 log_2 m); the
   construct referee's m ~ sqrt(2D/log_2 D) is the same statement.

### 1.4 extend: analyst (other designs, exhaustive data, odd characteristic)

**Proved**
* **Theorem E** (analyst). *Let V be any finite F_2-subspace of F_2[t] with 1 in V, T in F_2[t] arbitrary,
  Y=T+theta; take distinct forms l_i (l_i(1)=1), sides S_i, F_i=prod_{v in S_i}(Y+v), leading units u_i. If g is a
  common divisor of the F_i and the monic points u_i^{-1}F_i/g form a c=0 cluster of degree D', and s'' is the size of
  a subfamily with one side per form and a common u_i (s'' >= ceil(s/4); s''=s if T not in V), then*
  D' >= (4/s'') sum_{i<j} iota(P_{J_ij}(1)) >= (3/8)s''(s''-1) floor(floor(log_2 s'')/2),
  *where iota(x) is the degree of the odd-degree (inert) part of x.* **[J] I checked the proof line by line**:
  (i) F_i = P_J(Y+a)X, F_j = P_J(Y+a)(X+P_J(1)), so delta_ij = deg P_J(1) - deg gcd(X,P_J(1)), invariant under common
  division; (ii) an odd-degree prime of F_2[t] stays prime in F_4[t] and never divides (T+w)+theta; (iii) the
  elements v in J with P | 1+v form a coset of J cap K_P (size >= |J|/2^e) unless one of the two forms vanishes on
  K_P = ker(ev_P), and at most 2^{e-1} forms with l(1)=1 do; the count C(s''-k,2) >= (s''-1)(s''/2-k) holds since the
  difference is k^2-k >= 0; (iv) the record's margin identity sum_{i<j} delta_ij <= s''D'/4 (from
  sum deg g_ij >= D' s''(s''-2)/4 and the c=0 condition), |V| >= 2s'', e·pi(e)/2^e >= 3/4 for odd e. Correct.
  Scope: hyperplane designs with translation element 1 only (not conj(Y)=Y+B, non-hyperplane transversals,
  non-common division, mixed (V,T), prime powers). It strictly generalises V.2 and needs no (V.4).
  On Part V: (4/m) sum iota = 8, 84, 490, 2790 against D = 12, 112, 720, 3968.
* **Lemma P** (analyst). *In characteristic != 2, if A is finite, nonempty, A cap (-A) empty, and
  prod_{a in A}(X+a) - prod_{a in A}(X-a) is constant, then |A| is odd.* [J] proof checked (for even |A| the
  difference is 2·(odd part of g_A); constant forces it to vanish, so g_A is even and A = -A as multisets);
  `lemmaP.py`: over Z with A in {±1..±6}, |A|<=6, constant differences occur only for |A| = 1 (12 sets) and
  |A| = 3 (12 sets, e.g. {1,2,-3}). *Corollary:* three transversals of a symmetric set cannot be pairwise
  constant-difference (|D_12|+|D_13|+|D_23| is even), so the V.1 mechanism under any reflection-type conjugation
  gives at most 2 points, in F_3(t) and in Z[i] alike.

**Computed**
* Exact M_2(D,0) (anisotropic, F_4(t)/F_2(t)): 2,3,2,4,3,4,3,4,3,4 for D=1..10, and 4 for every D=11..17 (C,
  cross-checked in Python for D<=7). **[J] reproduced for D<=9 with my own code, including the number of maximal
  classes (1,1,16,1,2,6,28,32,222)** (`m2.py`, `m2.out`). Reduction used (and checked): for D > 2c+1 all points of a
  cluster share their leading coefficient, units have norm 1, so monic WLOG, and at c=0 the relation "same norm and
  same coefficients above floor(D/2)" is an equivalence.
* Exact M_3(D,0) (A^2+B^2, F_9(t)/F_3(t)): 2,4,3,4,3,4,4,4,4,4,4 for D=1..11 (not redone by the judge).
* Searches: 400 random V of dim 3-5 tie Part V (D = 12, 112, 720), never beat it; all sub-transversal (non-linear)
  products for dim V = 4 give only 4 points at D=4, 3 at D=18, 8 at D=112; dim 5 nothing new. Best exponent: 2.

**Gap (the analyst's own extension, not Astra's claim).** Odd characteristic: nothing constructed and nothing
disproved; even-size constant-difference sets that are not cosets are the open escape from Lemma P. **[J] ruling:**
Astra claimed only that *this implementation* needs characteristic 2 (2w=0 for a conjugation-invariant translation
w), which is correct; she did not claim anything about F_3(t). The GAP is correctly labelled as an open question,
not as a defect of Part V.

**Corrected.** "Which identity has no characteristic-0 counterpart": not merely conjugation-as-translation, but
"Delta_A constant for many pairwise differences A of even size". With a translation the X^{|A|-1} coefficient is
|A|·w (nonzero in characteristic 0); with a reflection Lemma P forces |A| odd and kills triangles.

### 1.5 lit: literature and meaning (referee: correct; wording risks)

**Checked (search snippets only; arxiv.org unreachable for that referee)**
* The worst-case, fixed-q, uniform-in-N function-field short-arc problem is **not in the literature**. Existing work
  (Bary-Soroker-Smilansky-Wolf; Bank-Bary-Soroker-Fehm; Gorodetsky; Rudnick-Waxman and Katz on sectors, the closest
  object; Yiasemides) is statistical, mostly as q -> infinity, and uses the ramified model S^2 = -T. Part V's
  statement appears new.
* **The mechanism is prior art and must be cited**: subspace polynomials with agreeing top coefficients,
  Ben-Sasson, Kopparty, Radhakrishnan, "Subspace polynomials and limits to list decoding of Reed-Solomon codes"
  (FOCS 2006; IEEE Trans. Inf. Theory 2010). It is also the elementary-symmetric (ideal PTE) form of additivity.
  Only the equal-norm specialisation via Artin-Schreier conjugation is Astra's.
* In characteristic 2 the only separable quadratic "circle" over F_2(t) with an inert infinite place is the
  Artin-Schreier one; there is no Kummer model. Part V lives where the analogue of Z[i] is least like Z[i].
* No published warning that the integer uniform conjecture is false was found.

**FALSE (a question the coordinator/referee posed, not Astra's claim):** "isotriviality is the reason for the
characteristic-2 failure". **[J] ruling: concur.** The characteristic-0 model is equally isotrivial (one norm-one
torus for every N) and is bounded (Part IV, <= 24); the operative feature is additivity plus Galois-as-translation,
and the point is separable (§1.3). The right parallels are Mahler's Artin-Schreier counterexample to Roth, the
characteristic-p failures of BMY (Lang; Easton, PAMS 136 (2008)), and Davenport-Stothers in characteristic p.

**Corrected reading** (lit): *no proof of the UNIFORM bound can consist only of steps valid verbatim for
F_4(t)/F_2(t); proved for p=2 only; growth-rate arguments are untouched.*

### 1.6 Withdrawn or corrected in Part V (collected)

1. V.2 does not treat unit-twisted (monic-normalised) points at deg T < d. Repaired by lower; loss <= 4 in (V.3).
2. V.2 thresholds d >= 48, s >= 6: the theorem is vacuous below Q = 2^48. Superseded by Theorem E (no thresholds).
3. "In particular the quadratic order in V.1 is sharp for this family": order right; V.1 itself is not optimal in its
   family (8 points at point degree 104 < 112 at d=4).
4. "m points on a parent of degree 2m(m-1)(log_2 m-1)" is the least L for T = t^L; the full statement holds for
   every L >= m(log_2 m - 1).
5. The mechanism (subspace polynomials) is not new: cite Ben-Sasson-Kopparty-Radhakrishnan.
6. N is not squarefree; say so.
No claim of Astra's is withdrawn.

## 2. What is genuinely new relative to the record, and correct

* **An unbounded c=0 cluster family over a fixed constant field with an anisotropic norm** (V.1). New: the record's
  characteristic-2 example (Part IV trace seed) had a growing field F_{2^d} and the split form; expert note §5
  item 12 records "no anisotropic norm fits it". grep finds no subspace-polynomial construction before round 31.
  Method classical (BKR), specialisation new.
* **A separable characteristic-2 rational-surface failure of the orbifold BMY bound (IV.2)** coming from genuine
  clusters (reconcile). New to the record; char-p BMY failure as such is known.
* **Refutation of the F_2(t) analogues of A_5(delta<1/30), (L_5°), V_E(5)** (reconcile). New calibration.
* **Transferred growth bound** W >= 3M log_2 M - O(M) over F_2[t] (reconcile). Routine; fixed-field form new.
* **Theorem E** (extend): quadratic loss for every hyperplane design with translation 1, located at small inert
  primes. New.
* **Exact tables** M_2(D,0), D <= 17, and M_3(D,0), D <= 11 (extend). New.
* **Lemma P** and its triangle corollary (extend): the precise reason the mechanism needs characteristic 2. New,
  elementary.
* **V.2** (Astra; lower): correct, new by definition, superseded by Theorem E.
* By-products: CC exponent 1/2 - 1/(2m-2) attained up to an additive degree constant at fixed m, large L
  (construct); 8-point family member at point degree 104 (lower).

Not new: the margin identity (`ref_divisor.md` item 5), the star dictionary (`growth/ptolemy.md`), J4's averaging
argument (round 23), subspace polynomials themselves.

## 3. Rulings

### (a) Is "sup_D M(D,0) is infinite for the anisotropic form over F_4(t)/F_2(t)" a THEOREM?

**Yes.** Statement for the record:

> **Theorem (round 31, Part V; referees construct, reconcile, extend, lit; judge [J] for m=4, 8).** Let M_2(D,c) be
> the largest number of distinct z_j in F_4[t] of degree D with z_j·zbar_j = N for one N in F_2[t] (any
> multiplicities), zbar the coefficient conjugation theta -> theta+1, and deg(z_i - z_j) <= D/2 + c. For r >= 2,
> m=2^r, d=r+1, V=F_2[t]_{<d}, every integer L >= m(r-1), Y=t^L+theta, the m points
> z_H = prod_{v in H\{0}}(Y+v) (H ranging over the hyperplanes of V not containing 1) have common norm
> N = prod_{v in V\{0,1}}(Y+v), degree D=(m-1)L, and max deg(z_H - z_H') = D/2 + m(r-1)/2 - L/2. Hence
> M_2(D,0) >= 2^r for every D in (2^r-1)Z with D >= 2^r(2^r-1)(r-1); in particular sup_D M_2(D,0) = infinity, and
> limsup_{D->infinity} M_2(D,0)/sqrt(2D/log_2 D) >= 1.

Exact rate: along D_m = m(m-1)(log_2 m - 1), M_2(D_m,0) >= m ~ sqrt(2D_m/log_2 D_m); parent degree
W = 2D_m = Theta(m^2 log m). Upper side (reconcile, transfer of the record's 2/3 argument; not redone by the judge):
M_2(D,0) <= (1/3+o(1))·W/log_2 W with W = 2D, i.e. <= (2/3+o(1)) D/log_2 D. So
**sqrt(2D/log_2 D)·(1+o(1)) <= limsup M_2(D,0) <= (2/3+o(1)) D/log_2 D**; the truth between is open. For general D
(not a multiple of 2^r - 1) only the exhaustive values are known; monotonicity in D is not proved.
Caveats: N is not squarefree (the squarefree analogue is open); characteristic 2 only; the form is
A^2+AB+B^2, not a "sum of two squares".

### (b) Is the Omega(s^2 log s) statement proved, and what is its scope?

**Yes, with corrections.** Astra's V.2 is correct (lower, line by line; d=4 reproduced; d=5 exhaustive), with the
unit-twist case repaired (judge reproduced the identity, d=3..6). Its scope is narrow: V = F_2[t]_{<d}, Y = T+theta,
affine half-space sides, a **common** divisor, and d >= 48, s >= 6, so it is empty in every computable range.
**Theorem E (§1.4) supersedes it**: for every additive V in F_2[t] with 1 in V, every T, unit twists, any common
divisor, no thresholds, reduced point degree D' >= (3/8)s''(s''-1)floor(floor(log_2 s'')/2), s'' >= ceil(s/4).
Within *hyperplane designs with translation element 1* the parent degree is therefore Theta(s^2 log s) (V.1 gives the
upper side). **Neither statement is a lower bound for M_2(D,0), for the Part V parent, or for the norm form**:
non-hyperplane transversals, other translation elements B, independent prime choices, non-common division, prime
powers and mixed designs are untouched. The fixed-field growth question stays **open**.

### (c) Did the analyst reach a better rate, growth scale, or odd characteristic?

**No on all three.** Best exponent 2 (parent ~ 2m^2 log_2 m); the only improvement is a constant (104 vs 112 at
m=8, lower). No growth-scale design; Theorem E says any growth-scale design must be non-hyperplane, with cofactor
traces of average degree O(log m) avoiding small inert primes. Odd characteristic: nothing constructed; Lemma P shows
that the V.1 mechanism with any reflection-type conjugation gives at most 2 points (in F_3(t) and Z[i]).
**Small exhaustive data:** M_2(D,0) = 4 for all D in [11,17] (5 points need D >= 18), M_3(D,0) <= 4 for D <= 11,
against M(8,0) = 8 over C(t). Small residue fields visibly suppress clusters at small D, yet sup_D M_2(D,0) is
infinite. The data are therefore no guide to the asymptotics, in either direction.

### (d) The strengthened filter for integer arguments

> **F3+ (judge's wording, combining reconcile and lit).** Call an argument *F_2-valid* if every step holds verbatim
> for F_2[t] with the anisotropic norm of F_4(t)/F_2(t), all N (prime powers allowed), deg in place of log, and the
> inert degree-one place at infinity in place of the archimedean place. Admissible steps include unique
> factorisation, exact norms, residue counts and pigeonhole at small primes (all residue fields F_{2^e}, e small and
> fixed), anisotropy and inertness, cut/chord/Plücker identities, the chain lemma, Liouville/toric forced
> divisibility, Schwartz-Zippel, and intersection theory on surfaces short of BMY. **No F_2-valid argument proves
> the uniform bound, nor any bound M <= f(W) with f(W) = o(sqrt(W/log W)).**

What it **does** affect:
* F3's escape clause "or collisions in residue fields smaller than m" is dead for uniform statements: here the
  residue fields are F_{2^e}, e < log_2 m, and clusters are unbounded.
* Any char-free proof of A_5(delta) (delta < 1/30), of (L_5°), or of (L_8) (the last asymptotically, pending the
  char-2 dictionary). The implications J4 (A_5 => uniform) and (L_8) => uniform remain valid; only *char-free proofs
  of their hypotheses* are excluded.
* Liouville/K1-class gap arguments and anisotropy/parity lemmas as stand-alone proofs of uniformity.

What it **does not** affect:
* Growth-rate arguments: the (2/3+eps) proof transfers and Part V is about (2/3)m above it. Any improvement of
  the constant 3 or of the log R/loglog R order is untouched (fixed-field growth scale is open).
* **The record's conditional routes via Vojta are over number fields and are not refuted.** Part V refutes only
  their F_2(t) analogues, which fail in characteristic p for known reasons (as Roth and BMY do); the point is
  separable, so "assume separability" is not a sufficient repair.
* Part IV's (IV.2): characteristic 0, passes the filter, but has no integer analogue.
* Arguments valid only in odd characteristic or only for Kummer-type circles (reflection conjugation): p odd is
  open, and Lemma P shows the mechanism itself does not port. Arguments that need squarefree N in an essential way
  (not merely as a reduction) are also not excluded.
* Archimedean, characteristic-0 and Baker-type inputs.

### (e) What changes for the integer problem

Bluntly: **nothing**.
* **Growth rate:** unchanged, (2/3+eps) log R/loglog R. Part V is consistent with the transferred bound and far from
  it.
* **Uniform bound M(C):** unchanged. Nothing is proved or disproved for any integer parent.
* **Questions A and B:** unchanged over Z. Their function-field (char 2) analogues fail, which only says that a
  proof must use number-field input.
* **Is it evidence against the uniform conjecture over Z? No.** (i) The mechanism (additive subspace polynomials,
  Galois acting as translation, ideal PTE solutions of every even size 2^j) has no integer counterpart. Over Z,
  constant cofactor differences require odd |A| (Lemma P) and cannot form triangles, and one-parameter polynomial
  templates carry <= 24 points (Part IV). (ii) The characteristic-0 function-field model, which shares
  isotriviality, anisotropy and the exact norm with Z, is bounded (Part IV). (iii) Characteristic-p counterexamples
  of this genre (Mahler for Roth, Lang/Easton for BMY) carry no evidential weight on their characteristic-0
  statements. (iv) Part V is far from growth scale, so it does not even show that small primes fail to constrain
  growth. Its value is purely methodological: it identifies a class of proofs that cannot work.

## 4. Sentences of the expert note that should be updated (recommendations only; the record is read-only)

1. §5 item 12, header: "four exchanges ... 34 referee blocks and five judge reports" -> add Part V (five more
   blocks, `round31/README_V.md`).
2. §5 item 12, Part IV bullet *Characteristic 2*: after "The configuration is a genuine split-form (hyperbola)
   cluster; no anisotropic norm fits it." add: "Part V: over the **fixed** field F_2 with the anisotropic norm
   A^2+AB+B^2 of F_4(t)/F_2(t), subspace-polynomial clusters give M_2(D,0) >= 2^r for every D in (2^r-1)Z,
   D >= 2^r(2^r-1)(r-1); so sup_D M_2(D,0) = infinity (N not squarefree; mechanism Ben-Sasson-Kopparty-Radhakrishnan;
   parent degree Theta(m^2 log m), far above the transferred bound W >= 3m log_2 m - O(m))."
3. §5 item 12, *Theorem (function-field uniform bound)*: add "Characteristic 0 is essential: false over F_2(t) with
   an anisotropic norm (Part V); (IV.2) itself fails for a separable rational-surface configuration over F̄_2."
4. §5 item 12, *Filter F3*: replace "must use characteristic 0 or collisions in residue fields smaller than m" by
   F3+ as in §3(d) above: for the **uniform** bound, small-residue collisions do not suffice; a proof must use
   characteristic 0, archimedean structure, or a feature failing for Artin-Schreier circles. Add: "growth-rate
   arguments are not affected; number-field Vojta is not refuted."
5. §5 item 10, *Function-field model*: "The transfer to Z is negative: the D/2 independent equations come from the
   infinite residue field" -> add "; over the fixed finite field F_2 the model is nevertheless unbounded (Part V), by
   an additive mechanism with no integer counterpart."
6. §3 (Question B) and §2 (Question A): a one-line calibration: "the F_2(t) analogues of A_5(delta<1/30), (L_5°),
   V_E(5) are false (round 31 Part V, reconcile); (L_8) asymptotically, pending the char-2 dictionary; number-field
   versions untouched."
7. Do **not** write "the conjecture fails in positive characteristic" (proved for p=2 only; p odd open) or
   "isotriviality explains it" (false).

## 5. Open sub-questions, ranked by honest probability of affecting the integer problem

All probabilities are the judge's estimates [J]; all are small.

1. **Odd characteristic: is sup_D M_3(D,0) finite for A^2+B^2 over F_9(t)/F_3(t)?** A bounded answer proved with
   small residue fields and no characteristic-0 input would be the first candidate mechanism surviving F3+ that
   could plausibly port to Z[i] (both have reflection conjugation). P ≈ 0.005. *First step:* exhaustive search for
   even-size constant-difference sets A in F_9[t]_{<3} that are not cosets (the only escape from Lemma P), then
   extend the M_3 table to D = 12, 13 restricted to parents of the form t^{2D}+1 and products of symmetric factors.
2. **Fixed-field growth scale over F_2:** clusters with deg N = O(m log m), or a proof of deg N >= omega(m)·m log m.
   A proof would be an F_2-valid growth improvement and would need some step beyond pigeonhole; a construction would
   extend F3+ to the growth rate. P ≈ 0.003. *First step:* exact max-clique among all monic points (all prime
   choices, non-linear transversals) on the d=4 parents of lower's family_min table, starting at Y=t^15+theta.
3. **Squarefree N over F_2:** is sup_D M_2(D,0) still infinite for squarefree N? If not, prime powers carry the
   phenomenon and squarefree-N integer arguments escape F3+. P ≈ 0.002. *First step:* compute whether a subfamily of
   V.1 with V chosen so that no two v, v' differ by a multiple of a small even-degree prime yields squarefree N
   (gcd(N, N') = 1), and test the cluster condition at d = 4, 5.
4. **Which step of the record's growth proof (W >= 3M log M) could be improved over F_2[t]?** The transferred bound
   and Part V differ by a factor of about m. P ≈ 0.002. *First step:* compute the exact optimum of the F_2[t] cut
   problem (reconcile's F2*) for M = 16, 32 and compare with M_2 data at D ~ 100-300 from randomised search.
5. **FF-(L_8) at a stated m:** re-derive the star/Plücker dictionary (ptolemy Thm B) in characteristic 2 and compute
   m = 128. Calibration only. P ≈ 0.001. *First step:* run reconcile's `dict_check.py` logic on eight-tuples of the
   m = 128 cluster (D ≈ 9.75·10^4) with sampled cuts.
6. **A function-field Cilleruelo-Córdoba theorem over F_2(t)**, to turn construct's observation (exponent
   1/2 - 1/(2m-2) attained up to an additive constant) into a sharpness statement. P ≈ 0.0005. *First step:* port
   the CC argument (k points on an arc force a divisibility among Gaussian factors) to F_4[t] and identify the
   additive constant.

## 6. Scores (0-10)

| block | score | reason |
|---|---|---|
| construct (V.1 reproduction, legitimacy) | 3/10 | Thorough; strengthened statement for all L >= m(r-1); confirms no gcd cost and genuine anisotropy. Model-level only. |
| lower (V.2) | 1/10 | Careful, exhaustive at d=5, repaired the unit-twist gap, found the 104 < 112 improvement; but certifies a negative result about one construction. |
| reconcile (against the record) | 3/10 | The most informative reading: (IV.2) fails separably, FF A_5/(L_5°) refuted at 15/4, transferred bound shows Part V is far from growth scale, F3+ correctly limited to uniform statements. (L_8) part asymptotic and dictionary unchecked. |
| extend (analyst) | 2/10 | Theorem E (correct, general, threshold-free), exact M_2/M_3 tables, Lemma P. No better rate and no odd-p construction. |
| lit | 2/10 | Found the uncited prior art (BKR), fixed the wording risks (isotriviality, "positive characteristic"); originals unread. |
| **Part V, on the goal** | **2/10** | A clean, correct theorem in the function-field model (fixed field, anisotropic norm, unbounded) and a sharp method filter (F3+) that kills small-residue-field routes to the uniform bound. It changes nothing for integer parents: growth rate, uniform bound, Questions A and B keep their status, and it is not evidence against the conjecture over Z. |
