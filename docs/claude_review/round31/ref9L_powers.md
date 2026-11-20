# Referee report: IX.0, bounded powers (twin and Gaussian circles)

Block refereed: the two paragraphs "PROVED (bounded powers in the twin)" and "PROVED (bounded powers on Gaussian circles)" of Part IX, §IX.0, of `astra_lemma_IX/answer.md` (lines 1973-1988). The squarefree main theorem (IX.0) is taken as given; every step used here was re-derived. Scripts (own code, exact arithmetic): `ref9L_powers/gauss_check.py`, `ref9L_powers/twin_check.py` in this directory.

## 1. Results

**Verdict.** Both paragraphs are correct. The exclusive-prime extraction at bounded exponents is sound in every step I checked. The corrections are minor: one step needs b >= 4 (it is supplied by b >= 216E^2), a label count is an overcount, and the constants can be made explicit. Nothing in the argument is uniform in E: the threshold J_0(E), the exceptional exponent and the bound all degrade polynomially in E.

**Theorem T (bounded powers in the twin). Proved.** Setting: a pool of K distinct primes in [P, 2P]. Draw s of them uniformly without replacement, write J = log(K/s), and fix C > 0 and an integer E >= 1. Assume J >= 40E log(1+3E) (in particular J >= 40). Then there is s_0(C, E) such that for s >= s_0 the following holds outside a set of samples S of probability at most exp(-sJ/(8E)). Let the exponents 1 <= e_p <= E (p in S), the two residue classes modulo each p^{e_p}, and an interval of length L <= C sqrt q, q = prod p^{e_p}, all be chosen after S is seen. Then the interval contains at most
  B <= b - 1, with b - 1 = max(2, ceil(100 E^2 log(2P)/J)),
admissible integers. So B <<_E 1 + log(2P)/J, with an implied constant that is an absolute multiple of E^2. C enters only through s_0.

**Theorem G (bounded powers on Gaussian circles). Proved.** Same pool, now consisting of split primes, with the same s, J, C and E. Assume J >= 400E^2 log(E(E+1)^3). Then there is s_0(C, E) such that for s >= s_0 the following holds outside a set of samples of probability at most exp(-sJ/(250E^2)). Let N = prod_{p in S} p^{e_p} with 1 <= e_p <= E chosen after the sample. Then every arc of length C N^{1/4} on x^2 + y^2 = N contains at most 4(b - 1) lattice points, where
  b - 1 = max(216E^2, ceil(90000 E^5 log(2P)/J)).
So M <<_E 1 + log(2P)/J, with implied constant O(E^5). The constants are mine; Astra's "sufficiently large E-dependent multiple" is consistent with them. The exponents, the generators pi_p (a fixed function of p, chosen before sampling), the unit class and the arc are all chosen after the sample.

**Non-uniformity in E.** In the twin, J_0 = 40E log(1+3E), the exceptional exponent is sJ/(8E) and the bound is O(E^2(1 + log 2P/J)). On circles, c_E = 1/(100E^2), J_0 is of order E^2 log E, the exceptional exponent is of order sJ/E^2 and the bound is O(E^5(1 + log 2P/J)). Both results are valid for each fixed E. Neither extends to exponents growing with s without a loss. Example: at E of order log s with P = s^{1+eps}, J of order log s fails J >= J_0(E).

**Consequence for the same-norm comparison (not refereed here).** That comparison uses Theorem G with E = 2, J of order 0.24 log l (so J -> infinity) and log(2P)/J = O(1). Theorem G as corrected supplies O(1) marks with exceptional fraction exp(-Omega(s log s)). The widened bin itself is outside this block.

## 2. Item-by-item verdicts

### (1) Twin with exponents 1 <= e_p <= E

- **Triangle division with prime powers. Proved.** Write "agree at p" for "same class modulo p^{e_p}", and h_ij = prod_{p agree} p^{e_p}. Then h_ij | y_j - y_i and t_ij = (y_j - y_i)/h_ij is a nonzero integer. With only two classes at p, three marks either all agree at p or exactly one pair agrees. Let G be the product of p^{e_p} over the primes where all three agree, and A_ij = h_ij/G. Then each remaining prime divides exactly one A_ij, and it does so with the full exponent e_p. The case where the two classes agree modulo a lower power p^v (v = v_p(r - r') < e_p) changes nothing in the blocks: p^v then divides t for each of the two disagreeing pairs, and it stays inside the cofactors. The identity t_12 A_12 + t_23 A_23 = t_13 A_13 holds. Computed: all ordered triples for q = 9·25·49·11 (48 marks, 103,776 triples; the classes mod 9 and mod 25 agree mod 3 and mod 5) and for q = 8·27·5. In every case: exactly one block, with exponent e_p; p^v divides each disagreeing-pair cofactor.
- **Pair averaging. Proved.** A prime separates at most b^2/4 of the b(b-1)/2 pairs. Hence mean log h >= W(1/2 - 1/(2(b-1))) and mean log|t| <= log C + W/(2(b-1)). The bound is independent of the exponents. Some triple has log|t_12 t_23 t_13| <= 3 log C + 3W/(2(b-1)), and each log|t| >= 0.
- **Heavy block. Proved.** The primes separating marks 1 and 2 carry log mass >= W/2 - log C, so one of A_13, A_23 has log >= W/4 - (log C)/2. Each prime contributes e_p log p <= E log(2P), and W >= s log P. So the block has >= s log P/(4E log 2P) - O_C(1) >= s/(5E) slots for s large (P >= e^40 s is automatic).
- **W <= E s log(2P). Proved.**
- **Labels (1 + 3E). Proved (overcount).** A slot is unused, or lies in one of the three blocks with an exponent. The hidden slots need no exponent, because the forced integer prod_{hidden} p^{e_p} with all e_p >= 1 determines its set of prime factors whatever the exponents are. So "the forced positive product still determines the subset" is correct, by unique factorisation.
- **Union bound, recomputed. Proved.** The log of the union bound is at most
  log 3 + s log(1+3E) + [3 log 2 + 9 log C + 9Es log(2P)/(2(b-1))] - sJ/(5E).
  With b - 1 >= 100E^2 log(2P)/J the coefficient term is <= 0.045 sJ/E. With J >= 40E log(1+3E) the label term is <= 0.025 sJ/E. The total is <= -0.13 sJ/E + O_C(1) <= -sJ/(8E) for s large. This confirms "a strictly negative log bound of order sJ/E". The conditional probability <= binom(K-s+k, k)^{-1} <= (s/K)^k with k >= s/(5E) is unchanged: the exponents enter only the known exposed blocks.
- **Simultaneity over exponents. Proved.** The exceptional event (label assignment, coefficient triple) is defined without reference to the e_p, so post-sample exponents are covered.
- **"In the Gaussian triangle, three distinct exponent levels can put a prime and its conjugate in two different blocks." Proved (computed).** If the values are u < v < w, the pair at (u, v) gets conj(pi)^{w-v} and the pair at (v, w) gets pi^{v-u}. At 5^2 13^2 17^2 there are 13,122 such slot events over ordered triples. Every three-level slot has non-trivial factors in exactly two blocks.

### (2) Gaussian factorisation with exponents

- **z_j - z_i = 2i eta H_ij t_ij with t_ij a nonzero rational integer. Proved.** Fix pi_p. Every point of norm N = prod p^{e_p} (split p only) is uniquely z = eta prod pi_p^{a_p} conj(pi_p)^{e_p - a_p}, with eta a unit and 0 <= a_p <= e_p. Fix eta, and put m = min(a_ip, a_jp) and M = max(a_ip, a_jp). Then z_i = eta H_ij u_i, z_j = eta H_ij u_j, where u_i has factor conj(pi)^{M-m} or pi^{M-m} at p according as a_ip = m or a_ip = M, and u_j has the other factor. Thus u_j = conj(u_i) exactly, with no unit. Then z_j - z_i = eta H_ij (conj u_i - u_i) = 2i eta H_ij (-Im u_i), and t_ij = -Im u_i is in Z. It is nonzero because z_i != z_j. Computed: every pair at 5^2 13^2 17^2 (27 allocations, 351 pairs, units 1 and i), at 5^3·13·17^2, at 5^3·13^3 and at 5^4·29^2. In each case: u_j = conj(u_i), t rational and nonzero.
- **log|H_ij| = (W - D_ij)/2. Proved and computed.**
- **D_ij >= W/2 - O_C(1). Proved.** Exactly D_ij >= W/2 - 2 log(C/2), from 2|H||t| <= C N^{1/4} (chord <= arc) and |t| >= 1.
- **Mean of log|t_ij| <= log C_0 + W/(4(b-1)). Proved.** Write |a_i - a_j| = sum_{theta=1}^{e_p} 1[exactly one of a_i, a_j >= theta]. Each threshold cut separates at most b^2/4 pairs. So sum over pairs of D_ij <= W b^2/4, and the mean D <= W/2 + W/(2(b-1)). Hence mean log|t| <= log(C/2) + W/(4(b-1)), which is sharper than Astra's log C_0. Computed: sum |a_i - a_j| <= e b^2/4 on 20,000 random configurations.

### (3) The extraction

- **sum_p w_p delta_p >= E D_12 >= W/3. Proved, with b >= 4.** The direction is right: |a_1 - a_2| <= e_p 1[a_1 != a_2], so E D_12 <= sum w_p delta_p. With replacement, E D_12 = (1 - 1/b)(mean over distinct pairs) >= (1 - 1/b)(W/2 - 2 log(C/2)). This is >= W/3 once b >= 4 and W is large. At b = 3 it fails by a constant. Astra's "b larger than an absolute constant" covers this, and the later b >= 216E^2 supplies it.
- **P(a_1 = a_2 != a_3) = sum alpha_a^2 (1 - alpha_a) >= delta_p^2/E. Proved.** By Cauchy-Schwarz with sum_{a=0}^{e_p}(1 - alpha_a) = e_p. Computed exactly with fractions on 20,000 random configurations, including the identity with the direct triple count.
- **Expected weight >= W/(9E). Proved.** Weighted Cauchy-Schwarz gives sum w delta^2 >= (sum w delta)^2/W >= W/9.
- **Expected count >= s/(18E^2). Proved.** Use w_p <= E log 2P, W >= s log P, P >= 2.
- **Repetition cost <= 3s/b, and the c_E statement. Proved.** For b >= 216E^2, E[count; distinct] >= s/(18E^2) - s/(72E^2) = s/(24E^2). So the conditional mean on uniform distinct ordered triples is >= s/(24E^2) >= 2c_E s, with c_E = 1/(100E^2). Since the count is <= s, the fraction f with count >= c_E s satisfies f >= c_E.
- **The Markov step. Proved.** The log|t_ij| are >= 0 because the t_ij are nonzero rational integers. The average over all distinct ordered triples of the triple sum is 3 × (mean over distinct pairs). The good set has relative density >= c_E. So some good triple has sum log|t| <= 3(log C_0 + W/(4(b-1)))/c_E. The averaging is over the right set.

### (4) Division and enumeration

- **Case check at a_1 = a_2 = x != a_3 = y. Proved and computed.** Let G_3 have pi-exponent min_i a_i and conj-exponent e_p - max_i a_i. After division, H_ij/G_3 = pi^{min(a_i,a_j) - m} conj(pi)^{M - max(a_i,a_j)}, where m and M are the triple's min and max (e_p cancels).
  - Case y > x: H_12/G_3 = conj(pi)^{y-x}, and H_13/G_3 = H_23/G_3 = 1.
  - Case y < x: H_12/G_3 = pi^{x-y}, and H_13/G_3 = H_23/G_3 = 1.

  Computed by p-adic valuation at both pi and conj(pi) for every ordered triple in all five test norms (17,550 ordered triples at 5^2 13^2 17^2 = 2,925 unordered, 6 × 2,925; 12,144 at 5^3·13·17^2; 3,360 at 5^3·13^3; 2,730 at 5^4·29^2).
- **The other slots are known after exposure. Proved.** All-equal slots are trivial. Three-level slots and the pairs (1,3) or (2,3) agreeing give known factors once the prime and (a_1, a_2, a_3) are fixed, even when a prime lies in two blocks.
- **Enumeration. Proved (overcount).** [E(E+1)^3]^s is an upper bound. Since e_p cancels, (E+1)^3 per slot suffices.
- **Recovering the hidden set. Proved and computed.** Write X for the hidden part of A_12 and U for the exposed part. Then X = (A_13 t_13 - A_23 t_23)/(U t_12) is determined. Its norm prod_{hidden} p^{|x_p - y_p|} has every exponent >= 1, so its prime support is the hidden set. Computed: in every triple with a hidden slot (9,720 ordered triples at 5^2 13^2 17^2, and 6,720, 1,152 and 960 in the other norms), X was recovered from the equation exactly, and factorint(N(X)) has support equal to the hidden primes.
- **Final union bound and parameter order. Proved, with explicit constants.**
  - Each |t| <= Y, where log Y = 300E^2(log C_0 + W/(4(b-1))). So 3 log(2Y) <= O_{C,E}(1) + 225E^3 s log(2P)/(b-1).
  - Conditional probability <= exp(-c_E sJ).
  - With J >= 400E^2 log(E(E+1)^3) and b - 1 >= 90000E^5 log(2P)/J, the total is <= O_{C,E}(1) - c_E sJ/2.
  - Order of choice: E and C, then J_0(E), then b as a function of (E, P, J), then s_0(C, E). The O_{C,E}(1) term does not depend on b.

## 3. Method

Own scripts, Python 3 with exact integer and Fraction arithmetic; sympy for CRT and factorisation. Astra's script was not imported. All assertions passed (runtime about 2 s each). I read lines 1931-1997 of answer.md, and grepped final9.txt and reply9.txt for the exponent claims. No literature was consulted for this block.
