# Referee report, round 31, block "subspace": Part VI.1 of Astra's answer

Source: `astra_lemma_VI/answer.md`, lines 1329-1400 (VI.1). Scripts:
`round31/ref6L_subspace/checks.py`, `round31/ref6L_subspace/golden.py` (run under `prlimit --as=3.5e9`).

## Verdict summary

| Item | Verdict |
|---|---|
| (a) Fixed-place estimate (VI.3): no recurring fixed relation; lower bound | CORRECT (it is an Evertse-type nondegenerate S-unit sum bound for "almost S-units") |
| (b) Size reduction (units of height O(L)) | CORRECT, with one wording fix |
| (c) Main proof (escaping sequence, Lambda, k, deg_v, derivations, log-Jacobian, dependence over K) | CORRECT |
| (d) Honesty and scope against the project cautions | CORRECT WITH CORRECTIONS (gcd caveat on the scope paragraph) |
| (e) One-chord test; varying-orders example v_n | CORRECT (both checked) |
| (f) Record's golden/Pell families | consistent: they sit exactly at exponent 4 (on a curve) |

I found no gap in VI.1. Within its stated hypotheses it is a theorem.

## (a) The fixed-place estimate (VI.3)

Hypotheses: a fixed number m of sequences x_j in K^x; h(x_j) = O(L); sum over v not in S of |log|x_j|_v| = o(L); liminf h(x_j/x_k)/L > 0 for j != k.

**Subspace Theorem used.** This is the projective number-field form (Schmidt; Schlickewei; see Bombieri-Gubler, *Heights*, Thm 7.2.2; Evertse-Schlickewei 2002): fix K, a finite S, and for each v in S linearly independent forms L_{0v},...,L_{nv} with algebraic coefficients. Then the x in P^n(K) with
prod_{v in S} prod_i |L_{iv}(x)|_v / |x|_v <= H(x)^{-n-1-delta}
lie in finitely many proper subspaces. Here the absolute values are normalised (product formula) and H is the absolute height. **No integrality is needed**: x ranges over all of P^n(K). So the "almost S-unit" hypothesis is legitimate. It enters only through two identities: sum_{v in S} log|x_j|_v = -sum_{v not in S} log|x_j|_v = o(L), and sum_{v in S} log|x|_v = h(x) + o(L). Nothing is enlarged and nothing is hidden in epsilon except o(L) against an L-sized margin. I could not fetch the original statements (arxiv.org and math.leidenuniv.nl are not reachable from this session), so the statement above is from memory of the standard references. It is the standard form.

**Minimal relation.** I checked this line by line. Take a fixed relation sum c_j x_j = 0 of minimal support s on some subsequence. s = 1 is impossible. s = 2 contradicts ratio separation. For s >= 3 put y_j = c_j x_j and y' = (y_1:...:y_{s-1}). At v in S take the s-1 of the s forms {Y_1,...,Y_{s-1}, -sum Y} other than the one with the largest value. Any s-1 of these s forms are independent, and there are finitely many choices, so a subsequence fixes them. Then
log prod = sum_{v in S} sum_{j<=s} log|y_j|_v - sum_{v in S} log max_j|y_j|_v = o(L) - h(y) + o(L).
This uses y_s = c_s x_s being an almost S-unit and |log max_j|y_j|_v| <= sum_j |log|y_j|_v| outside S. Dividing by prod_{v in S}|y'|_v^{s-1} gives log LHS = -s h(y') + O(1) + o(L). Since h(y') >= h(y_1/y_2) >> L, this is <= -(s-1+eps')h(y') for any fixed eps' < 1 once L is large. The Subspace Theorem then gives a fixed nonzero relation among y_1..y_{s-1}, of support < s, on a further subsequence. That contradicts minimality. Correct.

**Lower bound.** Suppose log|sum x|_{v0} <= max_j log|x_j|_{v0} - delta L. At v0 replace a largest coordinate form by sum X_j (still independent), and keep coordinate forms elsewhere in S. Then log LHS <= -delta L - m h(x) + o(L). Because h(x) <= C L, this is <= -(m + delta/(2C)) h(x). So a fixed exponent increment exists, and here h = O(L) is genuinely needed. The resulting fixed relation was excluded above. Correct. The estimate is the classical Evertse (1984) lower bound for nondegenerate S-unit sums at one place. The pairwise-separation hypothesis replaces "no vanishing subsum" and in fact implies it.

## (b) Size reduction

Correct. One wording fix: "remove zero subsums" should read "choose on a subsequence a representation of minimal length". Then no proper subsum vanishes, including subsums containing -G_e. This is what the nondegeneracy step needs, because the group containing the unit 1 is then also nonzero.

The rest checks:
- h(G_e) = log|G_e| <= L/2 = o(M).
- Grouping by "ratio height o(M)" is transitive.
- The group coefficients have height o(M) and are eventually nonzero.
- Inter-group ratio heights divided by M have positive limits, so (VI.3) at scale M applies.
- With a single group, 1 lies in it, which forces every unit to have height o(M) and contradicts the definition of M.

## (c) Main proof

1. **Escaping sequence.** Negate "there exist Z, rho, H_0". Enumerate the countably many nonzero Laurent polynomials f_1, f_2, ... over the fixed normal closure. Take Z_n = V(f_1...f_n), rho_n = 1/n and H_0 = n, which gives t_n outside Z_n. Every fixed nonzero f is then nonzero on t_n for n >= index(f). This is correct, including rho -> 0, which gives coefficient heights o(L). The final relation P has coefficients in K, so it lies in the enumeration.
2. **Lambda.** For an S-unit, h(u) = (1/2) sum_{v in S} |log|u|_v|. So h(U^nu)/L -> (1/2) sum_v |<w_v,nu>|, and Lambda is the integer kernel, hence saturated. The weights may be irrational, so the real kernel can be larger than Lambda tensor R. The text correctly claims separation only "for each fixed exponent difference", and that is all (VI.4) uses. Correct.
3. **Coefficient field k.** The diagonal subsequence gives an eventual zero/nonzero dichotomy for all countably many polynomials. The ideal of eventual zeros is prime: if pq is eventually 0 and p is not, then p is eventually nonzero, so q is eventually 0. k is finitely generated over K, and nonzero elements evaluate eventually with height o(L). Correct.
4. **(VI.4).** The terms c_nu X^nu satisfy the hypotheses of (VI.3):
   - heights O(L);
   - outside S they contribute <= 2h(c_nu) = o(L);
   - ratio heights >= h(X^{nu-mu}) - o(L) >> L.

   So deg_v is the limit, deg_v is multiplicative, deg_v of a sum is at most the max, and sum_{v in S} w_v = 0 by the product formula. Correct.
5. **Derivations.** D(Q_i) = D(Q_i - 1) because D(1) = 0. For Laurent f, deg_v Df <= deg_v f, since no new monomials appear and D maps k into k. Quotient rule: D(f/g) = (g Df - f Dg)/g^2 has deg <= deg f - deg g. This holds whatever the sign of the degree. Random symbolic test in checks.py, 0 violations. Hence deg_v(DQ_i/Q_i) = deg_v D(Q_i - 1) - deg_v Q_i <= deg_v(Q_i - 1). This uses deg_v Q_i = 0, which holds **at infinite places only**: |sigma(q_i)| = 1 there, but q_i is not a unit at the split primes in S. The text uses the identity only there, but should say so. Transfer to all infinite places: q_i is in Q(i), so sigma_v(q_i) is q_i or its conjugate, and |q_i bar - 1| = |q_i - 1| <= Delta, with weight e_v. Correct.
6. **Determinant.** J is multilinear in rows that are integer combinations of dlog B_e, dlog B~_e and dlog eta_i. A repeated dlog B gives zero, and dlog eta has entries in k. So F = J prod B_e B~_e is a Laurent polynomial. This holds for any integer exponents b_ie; the cut structure is not used.
   - Infinite v: deg_v F <= -4(1/4+eps)e_v + e_v (log prod|G_e|^2)/L <= -4 eps e_v.
   - Finite v in S: deg_v F <= 0, using only integrality at the finite places of S.

   For any monomial nu_0 in the support, sum_{v in S} deg_v F >= sum_v <w_v,nu_0> = 0, while the bounds give sum <= -4eps. So F = 0. Correct.
7. **The delicate point: does J = 0 give a fixed numerical relation?** Yes. The derivations {X_j d/dX_j} together with {d/dt_l} (dual to a separating transcendence basis of k/K, extended by DX_j = 0) form a basis of Der_K(k(X)). All 4x4 minors vanishing means dQ_2,...,dQ_5 are dependent in Omega_{k(X)/K}. In characteristic 0 this gives algebraic dependence **over K**, not merely over k. So there is P in K[Y] nonzero with P(Q) = 0 in k(X). Clearing denominators gives a formal Laurent polynomial whose coefficients lie in the prime ideal of eventual zeros. Its evaluation, which equals P(q(n)) times eventually nonzero denominators, is therefore 0 for large n. So P is a fixed nonzero polynomial vanishing on the tail of the sequence, which is the contradiction. Including the K-derivations of k is what makes the relation fixed. Without them one would get only a relation with moving coefficients. The argument is correct.

## (d) Project cautions and scope

- **Fixed-S exceptional sets.** K, S and b are fixed hypotheses. Z, rho and H_0 depend on (K, S, b, eps) and on the finitely many cut and orientation lists, and they are ineffective. The text says so ("neither an effective height threshold nor an equation"). This is honest. Nothing is claimed for all parents with S depending on the parent. The parent primes are unrestricted, but the blocks are not.
- **Drifting phase.** The phases are kept exactly in eta_i, as elements of k with height o(L). Correct.
- **Template-dependent threshold.** Uniformity holds inside the class (K, S, b) only.
- **Gaussian gcd. Correction to the scope paragraph.** The hypothesis is prod|G_e|^2 <= H with H the *primitive* norm, and the G_e Gaussian integral. A raw construction with a large common Gaussian factor, for example traces or Pell powers whose rows share primes, is covered only if some block decomposition of the *primitive* tuple meets (VI.2) with prod|G_e|^2 <= H_prim. If one keeps the raw blocks, the theorem only gives the weaker threshold Delta <= H_raw^{-1/4-eps}. Corrected scope sentence: "covers traces/products of units of finitely many fixed orders, fixed-base sparse digits, etc., **provided the gcd-reduced blocks are themselves bounded S-unit sums with subpower coefficients, or the Gaussian gcd is H^{o(1)}**." In the (*) profile itself the blocks are the primitive layers, so (*)-tuples are covered exactly when those layers are bounded unit sums.
- **What it covers.** Pell, unit and trace families in fixed orders with independent exponents per block; small-height seeds and linear maps; bounded-length polynomial expressions; coordinates with boundedly many nonzero digits in a fixed base.
- **What it does not cover.**
  - Varying fields or orders. Example v_n below.
  - An unbounded number of summands.
  - S growing with H.
  - Generic Gaussian integers. These are the actual content of Question B and of (*).

## (e) One-chord test and varying orders

**One-chord test.** Correct.
- t = Im(A prod(x_j +- i y_j)) is a nonzero integer with |t| <= P^{o(1)}, since 2|Im C| = |C||q-1| <= P^6 P^{-6+o(1)}.
- Splitting into minimal zero subsums and grouping by (VI.3) over Q with fixed S forces equal archimedean weight inside each subsum.
- A weight-0 term must exist to cancel -t, so every block has coordinate weights {0, 1}.
- The weight-6 and weight-0 terms are unique.
- Parity: the all-large and all-small products use k and 6-k imaginary factors, so they carry the same component (Re A or Im A). checks.py verifies this symbolically for all 64 sign patterns and all 64 choices of the large coordinate. Hence both coefficients vanish, a contradiction.
- A block with a zero coordinate has no weight-0 coordinate, which contradicts at once.

The negation and quantifiers ("some c > 0, eventually") are handled correctly by a diagonal subsequence.

**Varying orders.** Correct. v_{n+1} + v_{n+1}^{-1} - v_n - v_n^{-1} - 1 = 0 holds exactly. Numerically (all four conjugate sign choices, n = 10, 100, 1000) there is no proper zero subsum. h(v_{n+1}/v_n) is about (1/2) log n, while |v_{n+1} - v_n|/v_n is about 1/n. So the cancellation is of size about 2L, and (VI.3) fails in moving quartic fields. This is a correct and well-known kind of illustration.

## (f) Sanity test on the record's unit families

The golden family of `growth/families.md` (Prop. G) has J_r = f_{r+1} + i f_r = a phi^r - (i a/phi) q^r phi^r, a sum of two S-units of Q(i, sqrt5) with fixed coefficients. Its 4 blocks have product of norms 4H/5 <= H_prim. `golden.py` takes n = 100, 160, 220, all 56 five-point subsets, divides by the exact Gaussian gcd, and finds min Delta H^{1/4} = 3.70245917364 at all three n (log Delta/log H -> -1/4). So the record's unit five-point clusters sit **exactly at exponent 4** (arcs C sqrt R), on a curve. This is consistent with VI.1, which excludes only Delta <= H^{-1/4-eps} off Z. It also shows that 4 cannot be raised by a statement of the form "no points at all". Whether a *dense* unit family reaches exponent 4 is open; the volume heuristic predicts about one tuple there.

## Novelty and significance (blunt)

**Novelty.** The record has no comparable result. Earlier rounds avoided fixed-S Subspace because S would be the primes of N; here S is a fixed hypothesis and the use is legitimate. The method is standard Corvaja-Zannier-style Subspace machinery: an Evertse sum bound, then a formal Laurent model with degrees taken from limiting weight vectors, then a logarithmic-Jacobian degree count. The specific statement is new as far as I can tell.

**Significance.** Low to moderate. It proves essential exponent >= 4 > 15/4, so A_5(delta) holds for every delta > 0, but **only for blocks that are bounded fixed-field S-unit sums**. It closes the most natural *constructive* route to a dense counterexample to the A_5 route (fixed-order Pell, trace and sparse-digit families). It gives nothing for arbitrary genuine tuples, which are what Question B and (*) require, and Astra says so. It is not a step toward a proof of Question B.

Note: arxiv.org and math.leidenuniv.nl could not be reached from this session (network allow-list); the Subspace statement is quoted from standard references from memory.
