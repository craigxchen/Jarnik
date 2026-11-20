# Referee report on Part IX, block "crossgap": section IX.5 (the max-scale cross-gap on a growth-scale fake)

Object refereed: `astra_free_IX/answer.md` lines 2138-2229 (section IX.5), with the work notes
`ix_sharp_work.md` (an earlier copy of the same text, identical in content) and `ix_crossgap_work.md`
(an earlier two-copy, Frostman-only construction that IX.5 supersedes). Record statements used:
`round21/fakes.md` section 4 (Theorem LS and the paragraph after it), `round31/ref8_spread.md` section 4,
`round31/README_VIII.md` section 3(e) and section 5 item 3. Referee code is in
`round31/ref9_crossgap/` (`schottky_check.py`, `twelve_check.py`). Exact integer arithmetic was used unless
a line says otherwise. Astra's scripts were not imported.

## 0. Verdict

1. **(IX.8) is stronger than every record form of the cross-gap. Genuine angles do not satisfy it as written.**
   For the genuine angle, w = pi_j (or any element of Z[i]^* Z pi_j) gives sine 0. Genuine angles satisfy all
   three lines once w ranges over w not in Z[i]^* Z pi_j. The fake satisfies the lines with the owner's
   direction included, so it satisfies every record form a fortiori. The proof part of Astra's claim is
   CORRECT. Her "genuine angles satisfy it" remark needs the exclusion.
2. **The Schottky measure (IX.5, "sharp max-scale cross-gap") is CORRECT.** I checked every step. The
   congruence and growth claims were reproduced exactly on random words. The Fourier-decay input is the
   Bourgain-Dyatlov theorem (GAFA 27 (2017) 744-771). Its hypotheses hold. I could not read the full text:
   arxiv.org and dspace.mit.edu are outside this session's network. The citation is therefore checked from
   memory, plus the abstract-level search snippet. One measure works for all labels p >= p_0.
3. **The conditional measures on prescribed sums are CORRECT.** I rederived both cases and the weak-limit step.
4. **The application and the replicated profile are CORRECT WITH CORRECTIONS.** g is not an arbitrary fixed
   integer. It is a fixed integer at least g_*(mu). The height exponent d > max(2, 2/beta) is unspecified,
   because beta < epsilon(Gamma) and Bourgain-Dyatlov give no explicit epsilon. The IX.3 machinery (integral
   right inverse, deletion at parent primes) carries over to the column-replicated matrices, as she says.
   Its correctness is inherited from the IX.3 referee and was not re-refereed here.
5. **The finite sharp test (line 1 on the same model) and the elementary twelve-copy variant are CORRECT.**
   The constants (.38, 0.779, r^10/.38 = 0.216 < 1, b = 0.048) were recomputed. The exact union measure of the
   good set is 0.4326 for p ~ 10^6.
6. **Consequence (Proved here, new):** the fake also satisfies round 21's "sharp union large sieve" with
   constant c = 1. This is the separation form of Theorem LS item 2 with D(K,X) in place of 2D(K,X).
7. **What the cap covers.** It covers A_{2M} plus P3-char at every cutoff (as given by IX.3) plus the
   constant-one max-scale cross-gap in all three forms. The fakes have M = m ~ (rho/(2g)) omega(N'),
   W/(m log m) -> 4gd/rho, and prime heights P_+ ~ m^d with d a large unspecified constant. That is far
   outside the low-height window P_+ <= M^{2-delta}. In the window the cross-gap item of README_VIII stays
   open. The twelve-copy variant needs d > 2 and gives line 1 only, so it is also outside the window.

## 1. (IX.8) stated precisely

Angles phi_j are taken modulo pi/2 (quarter-turns). "For all w" always includes all four associates, so each
line requires both |sin| and |cos| of phi_j - arg w to be bounded below. Labels are rational primes
p_j = 1 mod 4.

* **(L1) finite sharp lattice test (the spread referee's test, `ref8_spread.md` section 4):**
  |sin(phi_j - arg w)| >= 1/p_j for all 0 != w in Z[i] with Norm w <= p_j.
  The ref8 test excludes w associate to the model's own factor. Astra's version does not exclude it.
* **(L2) prime-only, all heights (round 21 section 4, last paragraph):**
  |sin(phi_j - arg varpi)| >= c'/max(p_j, Norm varpi) for Gaussian primes varpi. Round 21 says "against external
  primes", with c' unspecified. Astra takes c' = 1, and every Gaussian prime, own primes included.
* **(L3) lattice, all heights (README_VIII section 5 item 3, "a fixed c"):**
  |sin(phi_j - arg w)| >= c_*/max(p_j, Norm w) for all 0 != w in Z[i].

**Proposition 1.1 (Proved).** Let phi_j = arg pi_j with Norm pi_j = p_j. Then
|sin(phi_j - arg w)| >= 1/sqrt(p_j Norm w) >= 1/max(p_j, Norm w) for every 0 != w not in Z[i]^* Z pi_j.
Hence (L1), (L2) and (L3) hold with constant 1 for every w outside Z[i]^* Z pi_j. For w in that set the sine or
the cosine is 0, so (L1)-(L3) as written in (IX.8) FAIL for genuine angles.

*Proof.* Im(pi_j conj w) is a rational integer, and so is Re(pi_j conj w). Im vanishes iff w/pi_j is real. Re
vanishes iff w/(i pi_j) is real. Since pi_j is primitive (p_j is a split prime), either case gives
w in Z[i]^* Z pi_j. Otherwise |sin(arg pi_j - arg w)| = |Im(pi_j conj w)|/(|pi_j||w|) >= (p_j Norm w)^{-1/2}.
The same bound holds for |cos| with Re in place of Im. The bound (ab)^{-1/2} >= 1/max(a,b) is elementary. QED

**Correction to the text.** Astra's sentence "The first line is the exact finite-scale test of the spread
referee, including the owner's genuine direction" is accurate about her fakes. It is not a test genuine
circles pass. For the purpose of capping, a fake that passes the stronger statement passes the record's
statement as well, so nothing in the cap is affected.

**Measure remark (Proved, standard).** (L1) is a positive-Haar-measure condition: its good set has measure
> .38 (section 5). (L2) and (L3) are Lebesgue-null conditions. (L3) implies badly-approximable-type behaviour
of the direction with a fixed constant. By Khintchine's theorem almost every direction violates it. (L2)
implies |sin| >= 1/Norm varpi for all Gaussian primes of large norm, and the sum over split primes of 1/p
diverges. So a fake satisfying (L2) or (L3) cannot be Haar-generic. A measure with a Lebesgue-null support is
necessary, and that is the role of the Schottky measure.

## 2. The Schottky measure (Proved, all steps checked)

**Statement (corrected wording only).** There are an absolute p_0, an absolute c_* > 0, and a Borel
probability mu on R/Z (y in quarter-turns, angle (pi/2)y) with |hat mu(k)| <= C_0 (1+|k|)^{-epsilon} for some
absolute C_0, epsilon > 0. For every y in supp mu:

* |sin((pi/2)y - arg varpi)| >= 1/Norm varpi for every Gaussian prime varpi of sufficiently large norm. Hence
  |sin((pi/2)y - arg varpi)| >= 1/max(p, Norm varpi) for all Gaussian primes varpi and all p >= p_0.
* |sin((pi/2)y - arg w)| >= c_*/Norm w for every 0 != w in Z[i].

The support is independent of p. One measure serves all labels p >= p_0.

**What the support is.** Take the reflections g_i(x) = a_i + 1/(x - a_i) in the unit circles centred at
a_i = iD (i = 1, 2, 3). Here D is a multiple of Q = prod l_{a,b}, and the l_{a,b} are 120 distinct split primes
> 50, one for each (a,b) in [-5,5]^2 \ {0}. The support is the image of the limit set Lambda of the index-two
(free, rank-two, Schottky, convex co-compact, non-elementary) subgroup Gamma under x -> (2/pi) arg B(1,x),
with B = [[1+7L, L],[7, 1]] unimodular. Every point of Lambda is a signed continued fraction
a_{i1} + 1/((a_{i2}-a_{i1}) + 1/(...)) with all partial quotients in {D, 2D, 3D} or {+-D, +-2D}, hence = 0 mod Q.
So the support is a set of **badly approximable directions with respect to Gaussian rationals (all lattice
directions)**, with constant c_*. In addition, every approximation of quality t^2|sin| < 1 is a combination
A P_k + B' P_{k-1} with |A| <= 5 and |B'| <= 4. By the CRT choice of L, such a combination has norm divisible by
some l_{a,b}. So no Gaussian prime of norm > max l is a good approximant. This gives the constant 1 against
**Gaussian primes**. (Astra's "|A| <= 5" covers both cases. I get |A| <= 4 in the first case and |A| <= 1 in the
second.)

**Steps verified (each by hand).**
* l | ||B(a,b)||^2: here u + Lv = sqrt(-1) v mod l, and v = 7a + b != 0 with |v| <= 40 < l. det B = 1.
* The recurrence for P_n is the standard convergent recurrence of the signed continued fraction. B carries
  convergent vectors (q_n, p_n) to P_n and preserves determinants. The condition D >= 102 + ||Be_1||/||Be_2||
  gives h_1 >= 100 h_0. Then h_{n+1} >= D h_n - h_{n-1} >= 100 h_n and h_{n+1} <= (2D+1) h_n.
* |sin(phi - arg P_n)| <= sum_{j>=n} 1/(h_j h_{j+1}) <= 2/(h_n h_{n+1}).
* The three ranges of t in the composite-norm lemma. The determinant bound h_n/t + 2t/h_{n+1} is < 1 in the
  middle range, at most 1.08 in the low range and at most 4.02 in the high range.
* P_n = B e_2 or B e_1 mod Q, alternating, because all partial quotients are 0 mod Q.
* Irrationality of limit directions, compactness, and the choice of p_0 (larger than the reciprocal of the
  minimal sine distance to finitely many prime directions and the axes).
* The lattice constant: n minimal with h_n >= 4t gives |sin| >= (199/200)/(4(2D+1) t^2) for t > h_1/4. The
  finitely many shorter primitive directions are handled by compactness.

**Computed (`schottky_check.py`, seeds 1 and 7, 20 random reduced words of length 14 each; L has 333 digits,
D has 336 digits).** All assertions passed: det(P_n, P_{n+1}) = +-1; h_{n+1}^2 >= 10^4 h_n^2;
h_{n+1} <= (2D+1) h_n; the residues of P_n mod Q alternate; every A P_k + B' P_{k-1} with |A|, |B'| <= 5,
(A,B') != 0, has norm divisible by some l_{a,b}. I also ran an independent test. For the direction of P_14, I
enumerated the standard continued-fraction convergents and the extreme intermediate fractions (j <= 8 or
j >= a_n - 8) of its slope, with both the vector and its i-rotation. Every w with h_1^2 <= Norm w <= 10^{-6} h_14^2
and t^2 |sin(arg w - arg P_14)| < 1 had norm divisible by a chosen l (about 26-28 such w per word, none
failing). This test is evidence: it covers only these candidates and finite words.

**Fourier decay.** The theorem used is Bourgain-Dyatlov, "Fourier dimension and spectral gaps for hyperbolic
surfaces", GAFA 27 (2017) 744-771 (arXiv:1704.02909). From memory: Theorem 2 there says the following. Let
Gamma be in SL(2,R), convex co-compact and non-elementary, with infinity not in Lambda, and let mu be the
Patterson-Sullivan measure on Lambda. Then there is epsilon = epsilon(Gamma) > 0 such that
|int e^{i xi Phi(x)} g(x) dmu(x)| <= C |xi|^{-epsilon} for |xi| > 1. This holds for every C^2 phase Phi with
Phi' != 0 on a neighbourhood of Lambda and every C^1 amplitude g, with C depending on Phi and g. Theorem 1 is
the case Phi(x) = x.

What I actually saw: only a search-result snippet confirming the title, venue and that "the key new ingredient
is a Fourier decay bound for the Patterson-Sullivan measure" using nonlinearity of Gamma and Bourgain's
sum-product exponential sums. The theorem number and exact hypotheses are from memory and were not read in full
text (arxiv blocked here).

Hypotheses in IX.5:
* Gamma has integral matrices of determinant +1, so Gamma is in PSL(2,Z). It is free of rank 2 and Schottky,
  hence convex co-compact and non-elementary.
* Lambda is a subset of the union of [a_i - 1, a_i + 1], so infinity is not in Lambda.
* Phi(x) = (2/pi) arg B(1,x) has Phi'(x) = (2/pi)/||B(1,x)||^2 > 0 (because det B = 1). It is real-analytic on a
  bounded interval containing Lambda, where arg has a continuous lift.
* The amplitude is a bump equal to 1 near Lambda.
* xi = 2 pi k.

So the hypotheses apply. Alternatives that would also serve: Sahlsten-Stevens (Fourier decay for non-linear
IFS / Gibbs measures, cited from memory) and Li-Naud-Pan (Kleinian Schottky groups; covers Fuchsian ones as a
special case, from memory). Kaufman's measures on numbers with bounded partial quotients and the
Queffelec-Ramare digit-restricted measures do not directly fit: the digits here are restricted to a sparse
signed set {+-D, +-2D} that is 0 mod Q. The constants epsilon(Gamma) and C_0 are not explicit.

## 3. Conditional measures on prescribed sums (Proved)

**Statement.** Let mu on T satisfy |hat mu(k)| <= C_0 (1+|k|)^{-epsilon}, and fix beta < min(1, epsilon).
1. There is g_0 = g_0(mu) with sum_{k != 0} |hat mu(k)|^r <= 1/2 for all r >= g_0.
2. For every fixed g > g_0 + 2/epsilon (Astra's g > g_0 + 3/epsilon + 2 is sufficient), there is A = A(mu, g)
   such that for every Theta in T there is a probability nu_Theta on (supp mu)^g intersected with
   {sum x_i = Theta}. For every nonconstant v in Z^g with h = ||v||_1 and every arc J,
   nu_Theta(sum v_i x_i in J) <= A (1+h)^2 |J|^beta.

**Checks.**
* The Fejer kernel gives mu(J) << |J|^{min(epsilon,1)} << |J|^beta.
* |hat mu(k)| < 1 for k != 0, because mu is atomless. Dominated convergence in r gives item 1, since
  |hat mu|^r is monotone in r. So mu^{*r} has a continuous density in [1/2, 3/2].
* **Case 1** (some value a has multiplicity >= g_0). Condition on the coordinates outside S_a. Then
  P(Sigma in I | rest) <= 3|I|/2. The event {Sigma in I, V in J} forces sum_{i not in S_a} (v_i - a) x_i to lie
  in J - aI, of length <= |J| + h|I|. Choose one i with b = v_i - a != 0, |b| <= 2h. Its Frostman bound under
  x -> bx gives K |b|^{1-beta} (|J| + h|I|)^beta.
* **Case 2** (all multiplicities < g_0). Distinct values v != v' have |(k + vl) - (k + v'l)| >= |l|. So at most
  one value class has |k + v_i l| < |l|/2, and at least g - g_0 + 1 factors have argument >= M/(8(h+1)). The
  frequency sum is sum_M 8M min(1, (C(h+1)/M)^{(g-g_0)epsilon}) = O((1+h)^2) once (g - g_0) epsilon > 2. So
  the joint density is O((1+h)^2).
* Divide by P(Sigma in I) >= |I|/2 and let |I| -> 0 along a subsequence with a weak limit. Portmanteau on open
  arcs (enlarge closed ones) gives the stated bound for every Theta, uniformly in v and J. The support claim
  holds because (supp mu)^g and {Sigma = Theta} are closed.

## 4. Application and the replicated profile (Proved with corrections)

**Corrected statement.** Let mu, epsilon, beta, g_0 be as above, fix an integer g > g_0 + 2/epsilon, fix
d > max(2, 2/beta), and fix .99 < rho < 1. For each C > 0, along Paley primes q = 3 mod 4 (q > q_0(C, rho, g, d)),
there are angle models with the following data:
* labels (j,a), with j in F_q^* and a = 1..g;
* exponent rows a_{i,(j,a)} = 1 + K_{ij} and e = 2;
* m = (rho/2 + o(1)) q marks;
* s = omega(N') = g(q-1), so m/s -> rho/(2g);
* distinct split primes in a single bin [P_0, P_0 e^sigma] inside [P, 2P], P = q^d, with
  sigma = O_{C,g}(q log q/P);
* W = 2 g (q-1) log P (1 + o(1)), so W/(m log m) -> 4gd/rho, independent of C.

These models satisfy everything IX.3 gives (P3-char at every finite cutoff, with the same labels, marks, angles,
windings and residues) together with (L2) at constant 1 and (L3) at the absolute c_*. Every angle phi_{(j,a)}
lies in (pi/2) supp mu.

Corrections:
* "for fixed g" in the opening statement must read "for some fixed g >= g_*(mu)". Any larger g also works.
* The height exponent d is a fixed but unspecified constant, because epsilon(Gamma) is not explicit.
* The O_C in "W = O_C(m log m)" can be dropped: the coefficient is 4gd/rho.

**Checks.**
* Replication: D_g = D (x) 1_g^T and H_g likewise. Taking U_g = U on the first copy and 0 elsewhere gives
  D_g U_g = I. V and F are formed as in (IX.4): H_g - D_g = 0 mod 2, and D_g F = H_g.
* The ratio of marks i and i' is Gamma_i/Gamma_{i'} with Gamma_i = omega_i/conj(omega_i) and
  omega_i = prod_a pi_{ia}. Distinct owners have disjoint supports, so Im(omega_i conj(omega_{i'})) != 0.
* Physical sectors of width b_*/g put arg omega_i in an interval of width b_*.
* The pair surplus is 2g log P_0 - o(1), since sigma sqrt(q log q) = o(1) for d > 3/2.
* The only odd exponent of row i sits at its own columns, so "own row as unique middle mark" holds. Deletion at
  a parent prime is column-local and therefore unchanged. Its correctness is the subject of the IX.3 block and
  was **not re-refereed here**.
* All row characters and mark phases are block-constant (2 sum_j K_{ij} sum_a phi_{ja}). They depend only on
  the coarse block sums, which the fibre measures preserve exactly. Integer shifts of one coordinate match real
  lifts, and the cross-gap is invariant under quarter-turns.
* Coarse noncharacters: Haar measure on the closed connected subtorus (the image of the rational subspace
  ker_R D_g) gives failure probability << sum_h (2q P^{-g/2})^h.
* Mixed words: conditional on the other blocks, a nonconstant block gives risk
  A(1+h)^2 (C_1 P^{-h/2})^beta, where w(v) >= h log P. Summed over at most (2gq)^h words this is o(1) iff
  q P^{-beta/2} -> 0, i.e. iff d > 2/beta.

## 5. The finite sharp test (L1) and the twelve-copy variant (Proved; constants computed)

**(L1) on the same model.**
* Good range: P in [2 h_n h_{n+1}, h_{n+1}^2/32]. The parallel direction has sine >= (1 - 2e-4)/(h_n h_{n+1}).
  Other w with ||w|| <= sqrt(2P) <= h_{n+1}/4 satisfy |sin| >= 1/(2 h_n ||w||) >= 1/P, using P >= 8 h_n^2.
* Log-fraction: b = (1 - log 64/log 100)/2 = 0.0485 > 0 (computed).
* Fubini over log P in [T, 2T] gives some P with mu(G_P) >= b/2.
* Uniformity for nu = mu restricted to G_P and normalised, with nu <= mu/a:
  - Plancherel through nu^{*g_0} <= a^{-g_0} mu^{*g_0} gives sum_k |hat nu(k)|^{2 g_0} <= (9/4) a^{-2 g_0}.
  - The uniform bound sup_{k != 0} |hat nu(k)| <= (1 + cos 2 pi t)/2 comes from uniformly small arc masses of
    k_* mu. For large |k| this follows from Fejer plus decay, giving <= c't + C'|k|^{-epsilon}. For the finitely
    many small k it follows from atomlessness.
  - The joint bounds carry over with a factor a^{-g}.
  - The support lies in the closed set G_P^g.
* A Paley prime q in [P^{1/d}, 2P^{1/d}] exists by PNT in progressions.
* Labels in [P, 2P] pass (L1), since Norm w <= p <= 2P and 1/P >= 1/p.

**Twelve-copy variant.**
* The forbidden set for label p, in quarter-turns, is the union over primitive directions modulo units with
  Norm <= p of intervals of length 4 arcsin(1/p)/pi.
* Count: (6/pi^2)(pi p/4) + O(sqrt p log p). Computed: 477,448 directions at p = 1,000,003 against a predicted
  477,466.
* The union bound gives Haar measure a_p >= 1 - 6/pi^2 - o(1) = 0.392 - o(1). The computed union-bound value is
  0.3921 at p = 1,000,003 and the exact union measure is 0.4326. The union bound is below .38 at p = 101 (0.370),
  but labels are >= P, which is large.
* Bathtub: |hat f_p(k)| <= sin(pi a)/(pi a) <= sin(.38 pi)/(.38 pi) = 0.7788.
* Sum over k != 0 of prod_{a=1}^{12} |hat f_a(k)| <= r^10 ||f_1||_2 ||f_2||_2 <= r^10/.38 = 0.216 < 1. So the
  12-fold convolution has a continuous density >= 0.78.
* The fibre density is <= .38^{-12}/0.78, uniformly. Each nonconstant word is uniform on the Haar fibre.
* Union bound: sum_h (24q)^h P^{-h/2} = o(1) iff d > 2.

This gives (L1), with the owner's direction included, jointly with the IX.3 conclusions, using elementary
estimates only, at any fixed d > 2. It does not give (L2) or (L3), which are Lebesgue-null (section 1).

## 6. Consequence for round 21's sharp union large sieve (Proved, new)

**Proposition 6.1.** Let phi be any angle model whose cluster angles satisfy the pair and doubling gaps (as in
Lemma LS1) and (L2) with constant 1, with all p_j <= X. Then the multiset of the cluster angles +-4 phi_j
together with the genuine angles +-4 arg varpi of all other split primes of norm <= X is (4/X)-separated in
R/2piZ. Therefore the union large sieve holds with constant D(K,X) = 2K + pi X/2, the same constant as for
genuine angles. Theorem LS item 2 had 2D(K,X).

*Proof.* Within the cluster: Lemma LS1. Cluster against external: (L2) for varpi and for i varpi gives
dist(phi_j - arg varpi, (pi/2)Z) >= arcsin(1/max(p_j, Norm varpi)) >= 1/X. Conjugates are also Gaussian primes.
Within the external family: Lemma LS1 (genuine). Then apply Montgomery-Vaughan with delta = 2/(pi X). QED

So the IX.5 fakes also satisfy the "1 <= c < 2" union sieve of round 21 section 4 at c = 1. This holds at
their heights, not in the window.

## 7. Verdict on the record items

* **Capped (conditional on IX.3 being correct, refereed elsewhere):** A_{2M}, plus P3-char at every finite
  cutoff, plus the max-scale cross-gap in each of the forms (L1), (L2) with c' = 1, and (L3) with an absolute
  c_*, plus the sharp union large sieve at constant D(K,X). The cap holds at growth scale
  M ~ (rho/(2g)) omega(N') with W/(m log m) -> 4gd/rho, and at prime heights P_+ ~ m^d.
  - Unspecified constants: epsilon(Gamma), beta, g_0, g, d (requires d > max(2, 2/beta)), c_*, p_0. All are
    absolute (independent of C, of the radius and of the cutoff), but none is explicit.
  - With (L1) only, the twelve-copy variant works at any d > 2 with g = 12.
* **README_VIII section 3(e)**, "union large sieve at the genuine constant (max-scale cross-gap)": this is no
  longer uncapped at large fixed prime height. It **remains uncapped in the low-height window**
  P_+ <= M^{2-delta}, where both constructions are inapplicable (d > 2 in both).
* **README_VIII section 5 item 3** asked for a positive-measure modification of the spread fake satisfying
  |sin(phi_p - arg w)| >= c/max(p, Norm w) at fixed c. What IX.5 achieves:
  - The form (L3) is achieved, but outside the window and not as a modification of the spread fake.
  - It cannot be "positive Haar measure": (L3) is Lebesgue-null (Khintchine). The law used is a product of
    fibre measures on a Lebesgue-null fractal.
  - What remains of item 3 is: (L1), or (L3), in the window P_+ <= M^{2-delta}, on window fakes (spread or
    concentrated), or a scarcity theorem there.
* Integer problem: unchanged. The best bound for general R is M <= (2/3 + eps) log R/loglog R, and M(C) is
  untouched.

## 8. Ledger

| Item | Label |
|---|---|
| (IX.8) as an assertion about the fakes | Proved |
| "genuine angles satisfy (IX.8)" | False as stated; Proved with w outside Z[i]^* Z pi_j (Prop. 1.1) |
| Schottky composite-norm lemma, lattice constant, irrationality | Proved (checked) and Computed (exact, 40 random words) |
| Power Fourier decay of the pushed-forward PS measure | Proved modulo Bourgain-Dyatlov, cited from memory, hypotheses checked; full text not read |
| Conditional fibre measures | Proved |
| Application and replicated profile | Proved with corrections (g >= g_*(mu); d unspecified; W coefficient 4gd/rho); IX.3 deletion inherited |
| (L1) on the same model | Proved |
| Twelve-copy variant | Proved; constants Computed |
| Union large sieve at constant D(K,X) for these fakes | Proved (new, Prop. 6.1) |
| Window status of the cross-gap | Open |
