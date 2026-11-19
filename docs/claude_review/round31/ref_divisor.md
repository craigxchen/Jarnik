# Round 31 referee report: block "divisor" (Astra, astra_free/answer.md sections 1, 2 and 4)

Referee scripts and outputs: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref_divisor/`
(`census_int.py`, `census_sqf.py`, `census_pp.py`, `four_vertex.py`, `genuine.py`, `recon4.py`, `fan.py`,
`nonprim.py`, `fake_diffs.py`; outputs `census_sqf16_out.txt`, `census_sqf22_out.txt`, `census_pp13_out.txt`,
`census_pp_astra_out.txt`). All of the code is mine. The single exception is the 417-digit four-point cluster in
`recon4.py`: I copied its construction from Astra's `near_square_check.py` and analysed it with my own code. I
also ran Astra's `scalar_chord_check.py` and `near_square_check.py`, and they print what the report says.

## Summary verdict

Every PROVED statement in the block survives line-by-line checking. I found no mathematical error, and the
constants are right. I reproduced every census exactly with independent code: the 10^7 census, the 22-prime
squarefree census (and the 16-prime one), and the 13-prime prime-power census on Astra's prime list. I also
extended the prime-power census to the first 13 split primes. There I found 218 incidences and 113 triples, with
no phase failures.

On significance the block is **peripheral**:
- It changes neither the growth rate, the uniform bound, nor the status of Question A, Question B or route 5.
- The divisor reduction (§2) is a clean repackaging of facts already in the record into a new counting target.
- The target **does not escape F1**. In the squarefree case it is literally the angle-model count of
  near-Liouville-extremal signed divisors, evaluated at the genuine angles. The record's flipped-Paley fakes
  make that count of order m^2 (answer to question (1) below).
- The triangle theorem is a single-monomial Liouville gap. Fakes satisfy it too.
- The reconstruction theorems (§2 five-vertex, §4 scalar) are correct, and they are new in exact form. The
  scalar theorem gives the first proof of the converse that round 8 asserted ("integrality of reduced chords is
  equivalent to the genuineness of the prime angles"), and the round-8 referee had flagged that converse as
  unproved (`round8/referee_stability.md` §5.4). These theorems are equivalences with the original problem.
  They are not inputs that would bound it.

## 1. Section 1

**Logical scope (PROVED): CORRECT, and trivial.** Suppose the models satisfy W <= b m log m. If m <= W, then
log m <= log W, so W/log W <= b m. If m > W, then W/log W < m. Either way W/log W <= max(b,1) m, so no
o(W/log W) bound can follow. The two quantifier remarks are correct logic:
- Separate models for two axiom lists say nothing about the conjunction of the lists.
- Constructions with b = b(A) for each fixed A do not rule out an argument that chooses A = A(epsilon).

**"Literal slip" in the expert note (EVIDENCE): CORRECT, cosmetic.** `expert_note.md` §1 says that no o(log
R/loglog R) bound is known "from any standard conjecture". The next sentence says that Vojta gives the uniform
bound, and the uniform bound implies an o-bound. The intended meaning is "no intermediate growth argument from a
weaker standard conjecture".

**"Missing pair input" (PROVED): CORRECT, but not new in substance.** That integrality of reduced chords is what
separates actual circles from fakes is already stated in `round8/stability.md` §7 ("What separates actual circles
from these models"). The precise converse is new, and it is the §4 theorem treated below.

**F2 remark (EVIDENCE).** The remark is that the unqualified claim "every determinant is critical" is
contradicted by route 5. That is the record's own statement: `expert_note.md` around line 387 already says
route 5 is supercritical. The overstatement is in the brief's paraphrase of F2, not in the record. Nothing new.

The "which openings I would fund" paragraph is opinion. I did not grade it.

## 2. Section 2: the reduction to primitive near-square divisors

**Reduction (PROVED): CORRECT, with all constants checked.** The steps are as follows.
1. Subdivide the arc into O(C) pieces of length <= (1/2) sqrt(R_orig) and keep one unit class (loss factor 4).
2. Divide by the Gaussian gcd. The arc becomes <= sqrt(R)/(2 sqrt|g|) <= N^{1/4}/2, and the residual N is odd and
   split-supported.
3. For an edge, z_i/z_j = A/conj(A) with A the conjugate-primitive Gaussian integer of norm
   d_ij = prod p^{|a_ip - a_jp|}. With the sign chosen so that arg A is half the positive angle,
   2t/sqrt d = |z_i/z_j - 1| <= 1/(2 N^{1/4}). Hence 16 t^2 sqrt N <= d, and t >= 1 gives d >= 16 sqrt N.
4. Floor uniqueness: 0 < sqrt d - s <= t^2/sqrt d <= 1/16, using d <= N.
5. Margin sum: sum_{i<j} L_ij <= m^2 W/4 - binom(m,2) W/2 = mW/4, because S(m-S) <= m^2/4 at every threshold
   level. This is valid for every profile. The record already has it in `round8/stability.md` (4.2): mean
   margin <= W/(2(M-1)).
6. Markov: at least half of the edges have L <= W/(m-1).
7. Failure of the rate: m >= M_orig/(4(2C+1)), W <= 2 log R_orig, and u/log u is monotone. Hence
   m >> W/log W and W/(m-1) <= B log W.
8. Injectivity: given d, the pair (s,t) is determined. An oriented edge is (g'A, g'conj A) with N(g') = N/d, and
   g' lies on an arc of length <= (1/2) N^{1/4}/sqrt d <= 1/8, which holds at most one Gaussian integer. So
   distinct edges carry distinct d. (Astra writes the edge as (g' conj A, g' A); that is only a change of
   orientation.)

D_B(N) >= (1/2) binom(m,2) follows.

The genuine data show that both normalisations are necessary:
- **Unit class** (`genuine.py`). On the record's 3-point cluster n = 276108509 (C = 0.4996), two of the three
  edges join different unit classes. Their half-angle representations have s^2+t^2 = 2 d_ij (for example
  965^2 + 1 = 2·465613), so they fail the test. Only the same-class edge passes.
- **Arc length <= sqrt(R)/2.** On the record 10-point cluster (C = 8.24), the edges (0,1) and (8,9) carry the
  same label d = 1032^2 + 1, and only 18 of the 45 labels are distinct. This agrees with the round-4 referee's
  finding that multiplicative Sidon holds only for C < 2 sqrt 2.

**Square-product triples (PROVED): CORRECT.** Fewer than 1/6 of the edges are discarded, each discarded edge
spoils m-2 triangles, and so at least half of the triangles survive. The exponent sums satisfy
|a-b| + |b-c| + |c-a| = 2(max - min), and each quotient (product of two labels over the third) is a square by
the triangle inequality and parity. Injectivity of the labels gives distinct triples. The count is
>> (W/log W)^3.

**Average estimate (PROVED): CORRECT.** sum_{X<=N<=2X} D_B(N) << B X^{3/4} log log X. The ranges are
s in [~4X^{1/4}, X^{1/4}(log X)^{B/2}], #t <= (s+1)/(4X^{1/4}), and #N <= 2X/s^2, giving a harmonic sum of
order (B/2) log log X. The determinant relation N gcd(alpha, beta) | Im(conj(alpha) beta) is correct. As Astra
says, neither statement gives anything pointwise.

**Nonprimitive counterexample (PROVED): CORRECT.** I checked it exactly at y = 200 and y = 400 (`nonprim.py`):
- 45 and 120 distinct d respectively.
- Every lower, angular and B = 5 upper cutoff holds.
- s = kx, t = k and gcd(s, t) = k.

Asymptotically there are ~ (y/4 log y)^2/2 values of k and W ~ y, so the count is >> (W/log W)^2 with an
implied constant near 1/32.

**Census to 10^7 (COMPUTED): REPRODUCED exactly** (`census_int.py`, independent smallest-prime-factor sieve for
the split-support test, explicit upper-cutoff test, with 0 rejections). The parents split 18,180 / 113 / 2 by
number of candidates (1 / 2 / 3). The two 3-candidate parents are 2,859,545, with
(77285, 278, 1), (219965, 469, 2), (2859545, 1691, 8), and 8,424,065. No triple has a square product. At 10^6 the
split is 2867 / 14 / 0, again with no triples. The completeness bound m <= sqrt(X/256) and the automatic upper
cutoff are both correct.

**Pointwise core restriction (PROVED): CORRECT, peripheral.** |k_j s_i - k_i s_j| = |k_i e_j - k_j e_i| <=
sqrt(H) T^2/sqrt D < 1. The two vectors are then proportional, so they are equal. Astra correctly says this
gives no saving at squarefree parents.

**Squarefree census (COMPUTED): REPRODUCED exactly** (`census_sqf.py`). Over the first 22 split primes (5..229),
with no upper cutoff:
- the maximum number of candidates at a parent is 6;
- there are 404 parent/triangle incidences, from 339 distinct triples;
- the exact Gaussian phase relation holds on all 339;
- there are 0 four-cliques with six distinct labels.

With 16 primes (5..149): 26 incidences and 26 triples, all closing.

My census also shows something Astra does not state. All 404 incidences give the four-vertex parallelogram
profile 0, v1, v2, v1+v2. This is forced: at a squarefree parent, a square-product triple is exactly such a
profile, with each label used twice. So **four-vertex exceptions occur at every squarefree triangle**.

**Triangle theorem (PROVED): CORRECT.**
- The signed valuations satisfy |k_jp| = v_p(d_j), by primitivity and odd norm.
- With local distances a, b, a+b, the absolute signed sums over sign patterns are 2(a+b), 0, 2a, 2b, whose
  average is a+b <= e_p.
- P_eps = r eta B^2 with N(B) = prod p^{|h|/2}, and some sign pattern has Q_eps <= sqrt N.
- sqrt(Q)|sin(u/2)| is an absolute coordinate of B rotated by a multiple of pi/4, so it is either 0 or
  >= 1/sqrt 2.
- |sin(u/2)| <= 3 pi/(16 N^{1/4}) < 1/(sqrt 2 N^{1/4}), because 3 pi/16 = 0.589 < 0.707. Hence u = 0.

Floor conditions are not needed, as Astra states.

Novelty and significance: this is a **single Liouville gap**, for the Gaussian integer B of norm <= sqrt N in
the axial and diagonal directions. The angle model imposes the same gaps (axial and diagonal monomial gaps), so
the theorem holds verbatim for fake angles and is inside F1. It explains the census, and nothing more.

**Prime-power census (COMPUTED): REPRODUCED** on Astra's 13 primes (`census_pp.py`): 53 incidences, 52 of them
at nonsquarefree parents, 7 with an imprimitive factor, 42 distinct triples, and no phase failures. On the first
13 split primes (5..109) the numbers are 218, 214, 76 and 113, again with no failures.

**Five-vertex reconstruction (PROVED): CORRECT.** The dichotomy psi_ij = |x_i - x_j| or (psi_ij = D and
x_i + x_j = D) is right. A third interior vertex rules out the second case, because there
|x_k - x_i| + |x_k - x_j| < D. Then:
- w_i is in Q(i);
- |h_i - h_j| = |a_ip - a_jp| on the complete graph, so h = ±a + c;
- the norm-N products and units give the points, and the arc length is 2 sqrt N · D.

The unique-diameter variant is correct as well. I verified it on a genuine cluster (`recon4.py`): the 417-digit
four-point cluster. All 6 edges pass the primitive floor test, the 6 labels are distinct, all 4 triangles close
exactly, the diameter is unique, and the line embedding is consistent to 10^-500. From the edge data alone, the
points z_0 · w_i reproduce the actual z_i with unit 1.

The record's 8- and 10-point clusters (`growth/families.md`: C = 5.67, 5.27, 8.24) **cannot test the theorem**.
Their largest sub-clique on which every edge passes the test has 2 vertices. The record's scan says the same
thing: at most 3 points on arcs of length sqrt(R)/2 for n <= 10^12. A genuine five-vertex test would need
radii like those of the cyclotomic family (~10^3428). I did not construct one. The reconstruction theorem itself
rests on the proof.

**Imprimitive variant (PROVED): CORRECT.**
- The magnitudes satisfy |k| <= v_p(d), with the same parity.
- (|x+u| + |x-u|)/2 = max(|x|, u) is monotone in u, so the sign-average bound a+b <= e_p and the evenness
  survive.
- The range of h at each prime is <= max |a_ip - a_jp| <= e_p.

The test "on every edge, d - floor(sqrt d)^2 = t^2 > 0 and 16 t^2 sqrt N <= d" is therefore, for m >= 5,
equivalent to a genuine one-unit cluster at norm N. That is a restatement of the original problem as an integer
clique problem.

**Four-vertex exceptional profile (COMPUTED+PROVED): CORRECT** (`four_vertex.py`):
- N = 13·29·89·109·157·229;
- d1 = 2540^2 + 1, supported on {13, 29, 109, 157};
- d2 = 5374^2 + 1, supported on {13, 89, 109, 229};
- d3 = 9633^2 + 2^2, supported on the symmetric difference {29, 89, 157, 229};
- d1 and d2 together cover the support of N;
- all three are primitive, pass the angular cutoff and divide N.

The exact relation is psi1 = psi2 + psi3, with psi1 = arg(2540 + i). (A1 A2 conj A3 is not real: the diameter is
d1, not d3.) The diameter label appears on two opposite edges, so no line embedding exists. By the census remark
above, such exceptions are as common as triangles.

**Search consequence (EVIDENCE).** This is a reasonable research suggestion. I did not grade it.

## 3. Section 4

**Scalar reconstruction (PROVED): CORRECT.**
- t_i s_j - t_j s_i = t_ij sqrt(d_i0 d_j0/d_ij), and the radicand is a perfect square (triangle inequality and
  parity).
- Squaring gives 2k t_j s_i in Z. With k != 0 and t_j != 0 (distinct phases), s_i is rational, hence integral.
- The second centre supplies the triangles (0, 0', i). These are needed because edges between two non-centres
  are untested, and the proof uses them correctly.
- The bound |h_i - h_j| <= |a_ip - a_jp| on tested edges, together with the four-case path bound through the
  straddling centres, gives a range <= e_p at every prime.
- The norm-N products, with units, realise the phases.

The hypothesis m >= 3 is essential. The round-8 referee's two-point example on |z|^2 = 13 (sqrt 12 ± i, with
t = 1 integral) is not genuinely realisable, which is consistent with this.

**Novelty:** this is the first proof in the record of the converse of round 8's claim "integrality of reduced
chords is equivalent to the genuineness of the prime angles", which the round-8 referee had flagged as unproved
and as false for two points. It is new and correct. Significance: it characterises genuineness. It does not
bound anything.

**Few edges suffice (PROVED): CORRECT.**
- The identity is sum (2S-m)^2 log p = m^2 W - 4 sum log d_ij <= mW.
- Each bad prime contributes >= m^2 log 5/4, so there are at most 4W/(m log 5) bad primes.
- The union bound holds: 2r (3/4)^{4 log 4r} = (1/2)(4r)^{-0.15} < 1.
- The number of incident edges is O(W + m log r + m) = O(m log m) when W = O(m log m).

Peripheral.

**Counterexample with arbitrary centres (COMPUTED): REPRODUCED** independently (`fan.py`). The scaled sines are
2, 1, 12, 1, 2 on the five tested edges and 7/5 on edge 23. The centres {000, 010} do not straddle at 5 or 17,
so the straddle hypothesis is needed.

## 4. Question (1): does D_B have a fake analogue? Yes, and the fakes violate the target

Take N squarefree and let phi be the genuine prime angles. A primitive Gaussian divisor of norm d corresponds,
up to units, to a pair ±v with v in {0,±1}^omega supported on d, and its angle mod pi/2 is <v, phi>. The floor
condition "d - floor(sqrt d)^2 = t^2 with 16 t^2 sqrt N <= d" is exactly sin ||<v, phi>||_{pi/2} <= N^{-1/4}/4.
Floor uniqueness shows that at most one class ±v per d passes. Therefore, **exactly**,

D_B(N) = #{±v != 0 : w(v) <= W/2 + B log W, sin ||<v, phi>||_{pi/2} <= e^{-W/4}/4}.

Since Liouville gives sin ||.|| >= e^{-w(v)/2}, this is the count of signed divisors whose angle lies within
W^{B/2}/4 of its Liouville limit. The prime-power case is the same, with |v_p| <= e_p.

In an angle-model fake the edge vectors v = a_i - a_j satisfy the same two conditions:
- the arc condition gives |<v, phi>| <= Delta/2;
- the Markov selection of §2 gives L_ij <= B log W on half the pairs, once W = O(m log m).

So the fake count is at least (1/2) binom(m,2) divided by the largest multiplicity of a difference vector. For
the refereed flipped-Paley fake (round5 construct, M = 12, r = 100; `fake_diffs.py`) all 66 differences are
distinct up to sign. In general the flipped-Paley family has distinct differences: H_x + H_v = H_y + H_w with
rows orthogonal off the all-ones column would give M - 2 = -2, and the degenerate coincidences force equal rows.

Hence fakes with W = O(m log m) have D^fake_B >> (W/log W)^2. **The target D_B(N) = o((W/log W)^2) fails in the
angle model, so it does not by itself escape F1.** The integer formulation hides the angles, but the
near-square condition *is* the genuine angle condition. Any proof of the target must use properties of the
actual arg pi that fakes lack. That is exactly F1's conclusion for the original problem.

What the reformulation does buy is the removal of all compatibility between edges: each counted object is a
single form. The F1 obstruction already lives at the level of single-form counts, though, because fakes satisfy
every single-form Liouville gap and still have about m^2 near-extremal forms. Compare `round8/stability.md` §8.4:
equidistribution counting reaches the Liouville scale and no further, even for one form.

## 5. Question (2): novelty of the reduction

The ingredients are all in the record:
- The per-pair statement: each close pair gives a divisor n = |u|^2 of N with s^2 + t^2 = n or 2n, gcd(s,t) = 1
  and one coordinate tiny. This is `growth/conditional.md` Prop 8.1, used there for abc. It includes the unit
  cases that produce 2n.
- The margin mean: `round8/stability.md` (4.2).
- Distinct ratios on short arcs: `round4/outside.md` Prop B.1(ii) and its referee, valid for C < 2 sqrt 2.
- Smooth Gaussian integers near squares at the fixed norm N'^2: `round8/stability.md` §8.3, which reports that
  no reduction was found.

New in Astra:
- normalising to one unit class so that s = floor(sqrt d) is explicit and unique;
- using injectivity to count *ordinary divisors of N*;
- the polylogarithmic window;
- stating the parent-only counting target D_B(N) = o((W/log W)^2), its triangle refinement, the average bound and
  the nonprimitive counterexample.

I found no such count in the record. Greps for "near-square", "divisor count", "Sidon" and "floor(sqrt" turned
up only the items listed above.

Verdict on novelty: **new as a formulation, routine as mathematics.** On significance it is peripheral: the
target is stronger than growth (it implies growth, not conversely), it may be false, and by §4 above it carries
the F1 obstruction in its single-form version.

## 6. Standing cautions

- Gaussian gcd and unit-class costs are handled explicitly; the genuine 3-point check shows the unit class
  matters.
- No fixed-S Subspace theorem is used.
- No common phase is set to zero. The reconstructions are "up to common rotation and point units", which is
  correct.
- No template-dependent threshold appears.
