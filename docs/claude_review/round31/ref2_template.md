# Referee report, round 31, block "template"

Scope: `astra_free_II/answer.md` section II.2. This covers the two rational degree-eight configurations C1 (seven
directions) and C2 (nine directions), the specialization cost, the CRT certificates, the progression argument, the
four-point exclusions, the minimum-shell cliques, and the ellipse and `t -> t^k + a` obstructions. It also covers the
lines of the Part II preamble that restate these results.

My scripts are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref2_template/`:
* `comb.py`: the direction census, done directly from the definition, without using the mass/balance lemma. For every
  vector in {-1,0,1}^8 (and, for midpoints, in {-2..2}^8) of l1-weight 4, I expand
  B_v(t) = prod (t - rho_j)^{v_j>0} (t - conj rho_j)^{v_j<0} over Z[i] or Q(i) and test whether the coefficients of
  t, t^2, t^3 are real. The script also computes the pair-sum classes, all 256 subset-sum fibres, and maximum cliques
  in {0,1}^8 whose pair differences are +-directions (all directions or only the selected ones, with or without
  distinct labels). It also lists the signed zero-sum triples.
* `sel.py`: maps Astra's selected labels onto my enumeration by support and runs the selected-label tests.
* `spec.py`: my own local search and CRT, then a check of every condition of Part I eq. (1). The candidate primes are
  those dividing some b_v or some resultant Res(f_i,f_j). For each prime p I search n mod p^k for the largest local
  contraction v_p(G), subject to every d_v being odd at p = 2. I then glue with modulus p^{1+v_p(N0(n_p))} and take n
  = residue + period*K, with random K of 7 and 31 digits; these are not Astra's representatives. I check (1) exactly
  with integer arithmetic, and W = ln N with 60-digit mpmath. The script also tracks G along 39 progression points.
* `alldirs.py`: the same pipeline using all 7 (resp. 9) directions.
* `tri_n.py`: an exact Gaussian-rational test of which signed triples of selected ratios B_v(n)/conj B_v(n) close at
  an explicit large n.

I also ran Astra's `ii_template.py` and `ii_cliques.py`. Their asserts pass and their outputs agree with mine.

## 0. Verdict

Every claim in the block is correct as stated, with the provisos Astra states itself. I found no error. I make two
small additions:
1. "For each fixed B > 0" holds only beyond a threshold in n that depends on B, and the threshold is astronomical for
   small B. Since d/sqrt N stays bounded (by about sqrt G = sqrt L, roughly 1.7e10 for C1), the condition
   d <= sqrt N W^B needs roughly W >= (1.7e10)^{1/B}. At the parents I computed (773 to 1436 digits), the least
   admissible B is between 1.02 and 2.42. This does not affect the claim. It does show the example is far from the
   asymptotic regime.
2. The C2 certificate is not locally optimal. My search finds further local contractions at 3 (exponent 2), 13 and 37,
   which multiply L by 4329 = 3^2 * 13 * 37. Astra's L already exceeds 256 max b^4, so nothing changes.

**Significance: peripheral.** The block exhibits a bounded number of divisors (5, resp. 7) at one parent. Section II.2
itself calls this "a fixed lower bound and gives no violation of the growth conjecture". It does not move the growth
rate, the uniform bound, or the status of Question A, Question B or route 5. It answers the coordinator's Question
2(b) only in a finite, non-decisive sense. Such configurations exist: 7 or 9 minimum-shell directions, with a largest
minimum-shell clique of 3. But the test the coordinator called decisive would need a family whose number of eligible
directions grows with log N at controlled specialization cost. The two obstructions block only two particular growth
mechanisms (axis-parallel ellipses, and Eisenstein compositions over Q(i)[t]).

## 1. Configurations, pair sums, fibres (COMPUTED). CORRECT

I worked from the polynomial definition, so the weight-four mass argument is not assumed.

| | C1 | C2 | cube (Part I) |
|---|---|---|---|
| ternary weight-4 directions, up to sign | 7 | 9 | 12 |
| all weight-4 integer kernel vectors (adds the midpoints) | 9 | 11 | 12 |
| fibre-size histogram | 1:120, 2:32, 3:8, 4:12 | 1:104, 2:32, 3:16, 5:8 | 1:104, 2:24, 4:24, 8:1 |

The pair-sum classes match Astra's exactly:
* C1: 04=13, 05=27, 06=17, 15=26=34, 36=47.
* C2: 03=14=27, 04=17, 05=16, 13=24, 15=26, 36=45, 46=57.

The midpoints are 2h_3 = h_2+h_7 and 2h_3 = h_0+h_5 in C1, and 2h_4 = h_3+h_7 and 2h_1 = h_0+h_2 in C2.

The constants b_v = |Im B_v| for the C1 selection 04=13, 06=17, 15=34, 26=34, 36=47 are 4800, 9600, 13200, 13200,
4800, as claimed. For the C2 selection (Astra's indices 0,1,2,3,5,7,8, that is 03=14, 03=27, 14=27, 04=17, 13=24,
36=45, 46=57) they are 49200, 98400, 49200, 67200, 31200, 13200, 4800. The maximum, 98400, is the number in Astra's
inequality. I also checked that b_v is the same at n = 0, 1, 7, -13 and 1000.

## 2. Specialization cost. CORRECT

d_v = N B_v(n)/g_v^2 and t_v = b_v/g_v, so d_v/(16 t_v^2 sqrt N) = sqrt(G) prod_{supp} f_j / (16 b_v^2 sqrt N0). The
factor prod_{supp} f_j / sqrt N0 tends to 1. With N = N0 (G = 1) the ratio tends to 1/(16 b_v^2) < 1, so the natural
parent fails the angular cutoff, exactly as stated.

## 3. CRT certificate and passage to the progression. CORRECT

**Persistence.** Suppose n ≡ n_p mod p^k with k > v_p(N0(n_p)). Then v_p(N0(n)) and v_p(N B_v(n)) are unchanged,
because N B_v = prod_{supp} f_j divides N0. The gcd g_v = gcd(Re B_v, b_v) is unchanged too. Indeed
v_p(g_v) <= v_p(N B_v)/2 < k. If v_p(Re B_v(n_p)) < k it persists; otherwise the minimum is v_p(b_v) at both points.
I checked this case split myself. So v_p(G) persists, and G is divisible by L along the progression.

**Bound on G.** Each f_j lies in some selected support, which is true for both selections. Then
f_j | g_v^2 d_v | b_v^2 N. Also prod f_j / lcm f_j divides prod_{i<j} gcd(f_i, f_j). The prime-wise inequality
sum e_i - max e_i <= sum_{i<j} min(e_i, e_j) holds: sort the e_i and keep the terms min(e_1, e_j). Each
gcd(f_i(n), f_j(n)) divides Res(f_i, f_j) ≠ 0. Hence G <= (lcm b)^2 prod Res. This bound has 142 digits for C1 and
162 for C2. Numerically, G = L exactly at all 39 progression points I tested, for both configurations.

**Distinctness, floor and parity.**
* Distinct d_v for large n: there are finitely many g_v | b_v, and distinct products of distinct monic irreducible
  quadratics over Q are never proportional.
* Floor: t_v stays bounded while s_v grows like n^4.
* Oddness: v_2 is fixed by the congruence modulo 2^16.
* Split support: each d_v is odd and a primitive sum of two squares.

**Explicit specializations (my own CRT).**
* C1: my local maxima coincide with Astra's table: (2,2,15), (3,0,2), (5,10,9), (13,5,1), (37,12,1), (41,17,1),
  (53,0,1), (509,238,1), giving L = 306439728192000000000 > 7.77e18 = 256 * 13200^4. At n of 50 and 74 digits, N has
  773 and 1164 digits.
* C2: my L is 541642686091490880000000000 > 2.40e22 = 256 * 98400^4. At n of 67 and 92 digits, N has 1041 and 1436
  digits.

At each of these four parents, every selected d satisfies every condition of eq. (1):
* d | N, with N odd;
* s = floor(sqrt d), d - s^2 = t^2, t > 0 and gcd(s, t) = 1;
* 16 sqrt N <= d, and 16 t^2 sqrt N <= d, tested exactly as d^2 >= 256 t^4 N;
* d <= sqrt N W^B for B = 2.5, and for B = 1.5 in the C1 case.

The divisors are pairwise distinct. All small prime factors of N below 2*10^5 are ≡ 1 mod 4, as the proof requires.
The offsets are t = (2, 12, 11, 11, 2) for C1 and (123, 24, 3, 84, 12, 11, 2) for C2, the same as Astra's. The ratio
d/(16 t^2 sqrt N) ranges over 6.3 to 47.5 for C1 and 150 to 6.3e4 for C2.

**Why only 5 (resp. 7) directions.** With all 7 (resp. 9) directions, my local search reaches only L = 7.5e12
(resp. 1.1e16), well below 256 max b^4. Every specialized ratio then fails the angular cutoff, with values 1e-4 to
0.29. Including every direction puts every pair (f_i, f_j) into a common support, which kills the resultant primes.
So the choice of a subset is forced by the method, not by convenience. This is consistent with Astra's choice, though
it is not a proof that no simultaneous specialization of all directions exists.

## 4. "No four lattice points realize the five selected labels" (C1). CORRECT

The logic is as follows.
1. If all five labels occur among the six edges of four points, five edges carry distinct labels. They form K4 minus
   an edge, which contains exactly two triangles. These share one edge and together carry all five labels, so they
   are two different triples of selected labels.
2. A triangle forces exact closure of r_a^{±1} r_b^{±1} r_c^{±1} at the integer n, with r_v = B_v(n)/conj B_v(n).
   - This identification is legitimate in the normalized short-arc setting of Part I. There the floor condition makes
     (s, t) a function of d, and the near-real condition excludes the unit -1.
   - Astra restricts the claim to this setting explicitly.
3. With phi_j = (t - rho_j)/(t - conj rho_j), r_v = prod phi_j^{v_j}. The phi_j are multiplicatively independent
   rational functions (distinct zeros and poles). So a signed triple with nonzero exponent sum is a nonconstant
   rational function and equals 1 at only finitely many n.

I verified that among the five selected C1 directions there is exactly one signed zero-sum triple: 04=13, 06=17,
36=47, that is, Astra's v0 - v2 + v6 = 0. At an explicit n of about 1.1e43 (`tri_n.py`), the only exactly closing
signed triples are that one and its global inverse. So two distinct closing triples cannot coexist, and the claim
follows. It concerns only the selected labels at that parent, as Astra says.

## 5. C2 has no four-clique among the seven selected ratios. CORRECT

The argument runs in four steps.
1. There are finitely many signed labellings of the 6 edges. So infinitely many bad n would give one labelling for
   which all triangle identities hold identically.
2. Normalizing one vertex to 1 gives exponent vectors in the free group on the phi_j.
3. Freeness turns every pair ratio into a pair difference ±v.
4. All differences ternary implies coordinate range at most 1, so a translate lies in {0,1}^8. This yields a clique in
   my graph.

The clique number with selected labels is 3 (computed). The C2 selection has two signed zero-sum triples:
{03=14, 03=27, 14=27} and {03=27, 04=17, 13=24}. This is why Astra uses the clique argument here instead of the
single-triangle argument. Astra's own script also shows that one fibre covers 6 of the 7 selected labels. That is
consistent with "no polynomial fibre accounts for all seven".

## 6. Minimum-shell cliques 4, 3, 3; fibre maxima 8, 4, 5. CORRECT

I reproduced these by brute force over {0,1}^8. The values are the same when distinct labels are also required, and
when only the selected directions are allowed in C1 and C2. In the cube, the four weight-8 edges form the antipodal
matching. So 4 is the obvious maximum and is attained, for example, by an explicit clique with six distinct labels.
The translation argument that extends this from {0,1}^8 to all of Z^8 is correct.

## 7. Ellipse obstruction. CORRECT (proof checked line by line)

The steps:
1. Im (x0 + Aw + B/w)^j is triangular in the Im w^i, with i ≡ j mod 2, i <= j and diagonal A^j - B^j. Here
   A > |B| because the semiaxes are positive. Real translation is triangular as well, so the moment equations
   transfer to the unit circle.
2. F has real interior coefficients by Newton's identities. Its constant term a_0 is nonreal: the w_k are distinct and
   lie in the open upper half-plane, so the root set is not closed under conjugation.
3. Self-inversiveness gives a_{m-j} = a_0 conj(a_j), which is equivalent to Astra's form. For real interior a_j and
   nonreal a_0 this forces a_j = 0. So F = z^m + a_0 and the selection is a whole fibre.
4. The upper parts of fibre(eta) ∪ fibre(conj eta) for distinct pairs {eta, conj eta} are disjoint sets of size m
   inside 2m roots. Hence there are at most two directions.

Novelty: the record's only "ellipse" (round10/additive.md) is unrelated. The result is new and peripheral. It rules
out only axis-parallel elliptical symmetry. The cube family is evidently not of this type (it has 12 directions).

## 8. Substitution t -> t^k + a. CORRECT, nearly trivial

The CRT choice of a with v_{pi_j}(a - rho_j) = 1 is valid. Z/p^2 → Z[i]/pi^2 is an isomorphism for split p, and the
conjugate valuation follows because a is rational. The choice makes every t^k + a - rho_j (and its conjugate)
Eisenstein, hence irreducible over Q(i). The divisors upstairs are then exactly F(t^k + a), and
Im F(t^k + a) = (Im F)(t^k + a) is constant iff Im F is. The statement is limited to Q(i)[t]-divisors and to this
choice of a, as Astra says. It is new to the record but carries little weight: irreducibility of generic compositions
is standard.

## 9. Novelty against the record

The template mechanism is round30/bridge.md Theorem 1.2. That theorem specializes a *cluster* and gains local
contraction from the Gaussian gcd g(t_0). Astra's variant specializes *unrelated* divisors. It replaces the parent N0
by lcm_v d_v and gains through G = N0/N. This adaptation is new and correctly executed. No earlier round counts
primitive near-square divisors at a common integer parent. Part I's census to 10^7 found at most 3; Astra now
exhibits 5 and 7 at roughly 700-digit parents.

## 10. Significance

The results are peripheral. Bounded examples cannot bear on D_B(N) = o((W/log W)^2). The construction cost (CRT
modulus, a contraction L of size 10^20 to 10^27, and B-dependent thresholds) grows with the number of directions,
and no family is offered. The block changes neither the growth rate nor the uniform bound, and it does not change the
status of Questions A, B or route 5.
