# Referee report, block "faketail": Part VIII.3, second half, and VIII.10

Scope: in `astra_free_VIII/answer.md`, lines 1817-1903. That is "Exact integer moment lifts", "Completion and the tail estimate", the final normalisation paragraph, the positive-Haar-measure paragraph, and the range statement (VIII.10) with its parameter changes. Scripts I wrote (exact integer arithmetic unless marked) are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref8_faketail/`: `lifts.py`, `newton.py`, `haar2.py`, `exponents.py`, `sprop.py`, `atoms.py` (the last is floating point). I read Astra's two scripts but did not import them. I did no literature search; nothing in this block needs one.

## Results

**Verdict.** Every step in this block is correct as written. I found no gap and no false statement. The step Astra asked to have checked first is the count of intermediate-length words whose initial moments vanish. That count is correct, including the move from multisets to vectors, unique recovery, and independence across digits. The block depends on two inputs from the first half of VIII.3, which other referees cover and which I took as given:
- (U) the uncertainty estimate: support of Kv mod l and mod 4 is at least sqrt(q) - h - 1;
- the character/offset estimates (VIII.9) and the residue witnesses.

The only open dependency is the gcd property "gcd((Kv)_{S^c}) <= 2". It comes from (U). VIII.10 uses (U) up to h <= q^tau, where tau can be close to 1/4, instead of h <= q^{0.10}. The proof of (U) as written has no restriction on h beyond sqrt(q) - h - 1 > 0. The union bounds close for every tau < 1/2.

| # | Statement | Label |
|---|---|---|
| 1 | Lifts d^(s), 1 <= s <= J = floor(q/4): d^(s) = j^s mod q, K_S d^(s) constant, sup norm O(q^{3/2} sqrt(log q)) | Proved; Computed at q = 103, 127 |
| 2a | D_s = d^(s) + q K_{S^c}^T xi_s keeps the congruence and constant K_S-image; sup norm <= E with probability 1 - o(1) | Proved |
| 2b | Atom bounds (A sqrt(n_0))^{-1} and gcd(w)/(A norm(w)_2) + O(exp(-c(A/T_0)^2)) | Proved (sanity: Computed, floating point) |
| 2c | Pr(v.D_s = 0) << q^{-0.62} h^{-1/2} (primitive balanced, h <= q^{0.10}) and << q^{-0.62} (balanced, 0 < h < q/2) | Proved, constants depend on rho; first bound conditional on (U) |
| 3a | Every term of theta except t has constant K_S-image on I | Proved |
| 3b | Newton: a balanced v != 0 with h < q/2 has v.D_s != 0 for some s <= h/(2 g_0) <= J | Proved; the bound s = h/2 is attained |
| 3c | All words with h <= H are safe on the bounded-digit event | Proved |
| 4a | Tail: unbalanced words; balanced noncharacters with h >= q/2; nonprimitive words with H < h < q/2 | Proved |
| 4b | Tail: primitive balanced words, the eligible count, the factor a^{k-1}, the risk from X_k, the exponent F | Proved; all exponents recomputed exactly |
| 5a | Final normalisation | Proved |
| 5b | Positive Haar measure on the realisation component, also after normalisation | Proved; tangent-space identity Computed at q = 43, 59, 103, 127 |
| 6 | (VIII.10): the parameter system closes exactly at xi = 1/12; limsup log P_+/log m <= 23/12; W/(m log m) <= 4d/rho -> 23/3 | Proved |
| 7 | Count over F_11: 40370, 3520, 440, 0 | Computed (reproduced) |

Quantifiers, stated once:
- rho in (0.99, 1) is fixed. In (VIII.10) only 0 < rho < 1 is needed.
- c_rho and the implied constants in the atom bounds depend only on rho. They scale like (1 - rho)^{-1/2}.
- All other implied constants are absolute.
- The starting threshold is q_0 = q_0(rho) for (VIII.8), and q_0 = q_0(rho, xi, kappa, tau, a, epsilon, lambda) for (VIII.10). It is not uniform in rho or xi. Astra says so, and the diagonal sequence is taken afterwards.
- Every "o(1)" holds as q tends to infinity with these parameters fixed.

The margins are asymptotic only:
- B/(HE) = q^{0.003};
- E/(q A sqrt(q log q)) = q^{0.005}/sqrt(log q);
- P_+ <= m^{1.995} needs q^{0.005} > 8.1, that is q > 10^{180} roughly.

So (VIII.8) is a statement for astronomically large q. This is not an error.

## 1. Exact integer moment lifts (Proved; Computed)

**Integrality of tau = Kd/q.** Modulo q, (Kd)_i = sum_j (j - i)^{(q-1)/2} j^s. This is the sum over F_q of a polynomial in j of degree (q-1)/2 + s. Such a sum vanishes mod q when the degree is below q - 1, that is s <= (q-3)/2. J = floor(q/4) is well inside this range. Negative control: at s = (q-1)/2, integrality fails at q = 103 and q = 127, as it must.

**Completion.** The congruence needed is K_{S^c}^T tau'_{S^c} = -K_S^T tau_S mod q. Row i of K mod q is j -> sum_k binom((q-1)/2, k) j^k (-i)^{(q-1)/2-k}. Hence the row space is the span of 1, j, ..., j^{(q-1)/2}. Any (q+1)/2 rows form a unit-scaled Vandermonde in the row labels, so they span that row space. |S^c| >= (q+1)/2 always holds, because S^c contains the (q+1)/2 rows outside the class A. So the completion exists.

**Height and constancy.**
- The deterministic part satisfies |sum_{i in S} K_ij tau_i| <= q max|tau| = O(q^{3/2} sqrt(log q)).
- Rounding with mean zero and step q makes the S^c part a sum of at most q independent terms of range at most 2q. Hoeffding plus a union bound over q coordinates gives O(q^{3/2} sqrt(log q)).
- d^(s) = d - K^T tau' satisfies d^(s) = d mod q.
- K_S d^(s) = q tau_S - (q tau'_S - (1^T tau') 1) is constant, because tau'_S = tau_S.
- After column 0 is deleted, the image changes by -d^(s)_0 times the constant column K_{S,0} = 1. So it stays constant.

**Computation (`lifts.py`, rho = 0.99, all s <= J):**

| q | |S| | max |tau| | sqrt(q log q) | max |d^(s)| | q^{1.5} sqrt(log q) |
|---|---|---|---|---|---|
| 103 | 50 | 10 | 21.8 | 713 | 2250 |
| 127 | 62 | 11 | 24.8 | 1027 | 3150 |

The check also confirms, in exact integers:
- K K^T = qI - 1 1^T;
- every tau is integral;
- every completion exists;
- every congruence d^(s) = j^s mod q holds;
- K_S d^(s) is constant on all q columns and on I.

## 2. Randomised digits and atom estimates (Proved)

**Digits.** K_S K_{S^c}^T = -1 1^T because S and S^c are disjoint. So q K_{S^c}^T xi_s has constant K_S-image, also on I, and it is 0 mod q.

**Sup norm.** The lattice Gaussian satisfies E e^{u xi} <= e^{A^2 u^2/2}. Indeed E e^{u xi} = e^{A^2u^2/2} theta(A^2 u)/theta(0), where theta(y) = sum_n e^{-(n-y)^2/(2A^2)}. The function theta is maximal at integers, because its Fourier coefficients are positive. A coordinate of q K_{S^c}^T xi_s is q times a sum of at most q such terms. With probability 1 - O(q^2 exp(-c q^{0.01})) all J(q - 1) coordinates are at most q^{1.625}. This is because q^{0.625}/(A sqrt(q)) = q^{0.005}.

**Atom 1.** Pr(w.xi = b) <= integral over R/Z of prod_i psi(w_i x), where psi = theta-normalised characteristic function, which is positive by Poisson summation. Apply Hölder with n_0 factors and use that x -> w_i x preserves Haar measure on R/Z for w_i != 0. This gives the bound integral of psi^{n_0}, which is at most (1 + O(e^{-cA^2}))^{n_0} / (A sqrt(2 pi n_0)) + O(e^{-c A^2}). Since n_0 <= q and A^2 = q^{0.24}, the first factor is 1 + o(1).

**Atom 2.** I re-derived it.
- If ||w_i x|| < 1/(3T_0) for all i, write w_i x = k_i + delta_i. Then |w_j k_i - w_i k_j| = |w_i delta_j - w_j delta_i| < 2/3, so w_j k_i = w_i k_j.
- Then r = k_i/w_i is common to all i. Since r w_i is an integer for every i, Bezout gives r in (1/g)Z with g = gcd(w).
- Near each of the g points the integrand is at most e^{-2 pi^2 A^2 ||w||^2 y^2}(1 + o(1)). Each contributes 1/(sqrt(2 pi) A ||w||).
- The complement contributes at most sup_{||y|| >= 1/(3T_0)} psi(y) = O(e^{-c(A/T_0)^2}).

Floating-point sanity check (`atoms.py`, A = 3, 6, several w): both bounds hold, essentially with equality in the extremal cases.

**Consequences.**
- v.D_s = v.d^(s) + q (w.xi_s) with w = (Kv)_{S^c}, where v is extended by zero on D. Hence Pr(v.D_s = 0) <= sup_b Pr(w.xi_s = b).
- (a) Primitive balanced v with h <= q^{0.10}. We have T_0 <= h, gcd(w) <= 2, and ||w||_2^2 >= (1-rho) q ||v||_2^2 / 2 >= (1-rho) q h / 2, using ||v||_2^2 >= ||v||_1 for integer vectors. So the atom is O((1-rho)^{-1/2} q^{-0.62} h^{-1/2}) + O(exp(-c q^{0.04})).
- (b) Balanced v with 0 < h < q/2. All |v_j| < q, so v is nonzero mod q, with support below q/2. Then (Kv)_i mod q is a polynomial of degree at most (q-1)/2 in i. It is not identically zero, by the moment Vandermonde. So it has at least (q+1)/2 nonzero values. At least (q+1)/2 - |S| >= (1-rho) q/2 - q^{2/3} of them lie in S^c. So n_0 >= c_rho q and the atom is O(q^{-0.62}).

**Where primitivity and the gcd enter.** Primitivity is used only in (a), to get gcd((Kv)_{S^c}) <= 2. The factor 2 cannot be removed: (Kv)_i = sum_j v_j - v_i = v_i mod 2 for balanced v. So (Kv)_{S^c} is all even unless S^c meets the odd support of v. That is why the test is made mod 4. Any prime dividing a nonzero w is at most T_0 <= h. So testing the odd l <= h and 4 suffices, once w != 0, which follows from the norm bound. For nonprimitive v = g_0 v', the gcd is g_0 times larger, and Astra treats these words separately (§4).

Sanity check (`sprop.py`, q = 103, 127, rho = 0.9, all h = 2 words and 2*10^4 sampled h = 4 words): the largest gcd is 2, and the norm inequality holds with ratio at least 6. This is only a sanity check: at these q, q^{0.10} < 2.

## 3. Angle formula, Newton step, short words (Proved)

**Constant images.** Each term other than t has constant K_S-image on I: c_0 1_I has image -c_0 1, since row sums vanish and K_{i0} = 1 on A; the D_s have constant image by §§1-2; and K_S K_{S^c}^T zeta = -(1^T zeta) 1.

**Newton.** Let v be balanced with h < q/2 and t = h/2. Its positive and negative label multisets are distinct multisets of size t over F_q. The labels are distinct elements of F_q, and the two supports are disjoint.
- Newton's identities up to degree t < q recover e_1, ..., e_t, because only 1, ..., t are inverted.
- So equal p_1, ..., p_t would make the multisets equal.
- Hence sum_j v_j j^s != 0 mod q for some s <= t, and then v.D_s != 0.
- Since t < q/4 and t is an integer, t <= floor(q/4) = J.

For v = g_0 v' with g_0 >= 2: g_0 <= h/2 < q, so q does not divide g_0, and the index is at most h/(2 g_0) <= h/4.

The bound s = t is attained: for a coset of the order-3 subgroup of F_{103}^* against another coset, the first nonzero moment is s = 3 = t (`sprop.py`). Over 3000 random balanced v at each of q = 103, 127, the first nonzero index never exceeded t.

**Short words, h <= H, on the event max|D_{s,j}| <= E.**
- Unbalanced: |(sum v) c_0| is in [1/(20H), 1/10]. The digit series is at most 2hE/(100H^2(B-1)) = o(1/H), because E/B -> 0.
- Balanced: the leading term is at least B^{-k}/(100H^2). The later terms are at most 2hE/(B-1) = o(1) relative to it, since HE/B = q^{-0.003}. The other perturbations are tiny:
  - h Delta, with Delta = q^{-(0.995+o(1)) q};
  - h q B^{-J-2} <= h q B^{-k-2}.
- Hence the gap is at least B^{-h/2}/(200H^2) > 2 P^{-h/2}, because P/B = q^{0.282} > 400 H^2 = 400 q^{0.16} and h >= 2. The total phase is o(1), so the half-turn is also avoided.

## 4. The tail estimate without conditioning on digit size (Proved)

**Counts.** The number of integer vectors with l^1 norm h, put h tokens into 2(q-1) signed bins, is at most binom(h + 2q, h). This is at most (Cq/h)^h for h <= q, and exp(O(h)) for h >= q/2. Correct.

**Unbalanced words.** Condition on everything except c_0. Then sum(v) c_0 mod 1 has density at most 1 + 20H. The risk is at most 8(1+20H) P^{-h/2}. Summed over h > H, the per-token exponent is 1 - gamma - 0.995 < 0 for gamma >= 0.08, and exp(O(h)) P^{-h/2} for h >= q/2. Correct.

**Balanced noncharacters with h >= q/2.** If (Kv)_{S^c} = 0, then qv = K^T K v = K_S^T (Kv)_S, using K^T K = qI - 1 1^T and sum v = 0. Also 1^T K v = 0, so v is a zero-sum rational combination of the S rows, that is a character. Otherwise zeta_i has a nonzero integer multiplier times B^{-J-2}, and the density is at most 1 + B^{J+2}. In base-q logarithms, B^{J+2} P^{-h/2} exp(O(h)) has exponent at most 0.427 q - 0.995 h + o(h), which is at most -0.141 h + o(h). Correct.

**Nonprimitive balanced words, H < h < q/2.** There are at most h (Cq)^{h/2} of them. Conditional on all digits, X_k has coefficient (v.D_k) B^{-k}, which is a nonzero integer multiple of B^{-k}, with k <= h/4. So the density is at most 1 + B^{h/4} and the risk is at most 8(1 + B^{h/4}) P^{-h/2}. The exponent per token is 1/2 + 1.708/4 - 0.995 = -17/250 < 0. Correct.

**Primitive balanced words: the count Astra asked about first.** Fix h, t = h/2, and k <= t. A word is *eligible* at (h, k) if p_1, ..., p_{k-1} agree mod q.
- *Vectors versus multisets.* v corresponds bijectively to its ordered pair (P, N) of disjoint t-multisets: P is the positive part with multiplicities, N the negative part. Counting pairs counts vectors exactly.
- *Recovery.* Let N_1 be the t - k + 1 smallest labels of N, with multiplicity, and N_2 the rest, of size k - 1. Then p_s(N_2) = p_s(P) - p_s(N_1) mod q for s < k. Newton up to degree k - 1 < q determines e_1, ..., e_{k-1} of N_2, hence N_2 as a multiset. So (P, N_1) determines v. The count is at most binom(q+t-1, t) binom(q+t-k, t-k+1).
- *Computation (`newton.py`).* For q = 7 (t <= 4), q = 11 (t = 3, 4) and q = 13 (t = 3, 4), over the actual label set F_q \ {0}, the map (P, N_1) -> v is injective, every count is below the bound, and the count is 0 when t moments agree.
- *Binomial bound.* binom(q+t, t) <= (e(q+t)/t)^t, and (h/u)^u <= e^{h/e}. This gives the stated exp(O(h)) (q/h)^{h-k+1}.
- *Probability that the earlier digits vanish.* Given v, the events {v.D_s = 0} for different s depend on the independent xi_s, the d^(s) being fixed. So they are independent, and the probability that digits 1, ..., k-1 vanish is at most a^{k-1}. If a congruence fails at some s < k, that probability is 0.
- *Risk.* The X_k are independent of all xi. So Pr(digits s < k vanish, digit k nonzero, bad) <= a^{k-1} * 8(1 + 100 H^2 B^k) P^{-h/2}. Only integrality and nonvanishing of v.D_k are used, as she says.

**Exponent.** Recomputed exactly with rationals in `exponents.py`, with gamma = log_q h and nu = k/h in [0, 1/2]. The base-q exponent per token is

F = (1 - gamma)(1 - nu) + (b - alpha(gamma)) nu - d/2 + o(1),

with alpha = 0.62 + gamma/2 in the strong range and 0.62 in the weak range. The coefficient of nu is 0.088 + gamma/2 (strong) or 0.088 + gamma (weak). Both are positive, so the maximum is at nu = 1/2.
- Strong range: F = 49/1000 - 3 gamma/4 + o(1), which is at most -11/1000 at gamma = 0.08.
- Weak range: F = 49/1000 - gamma/2, which is at most -1/1000 at gamma = 0.10.

Both of Astra's displayed inequalities (< -0.009 and < -0.0005) are therefore correct. The o(1) collects exp(O(h)), the q^{0.62} h^{1/2} from a^{k-1} against a^k, and H^2. All of these are q^{O(1)}, and h >= q^{0.08}, so their per-token contribution is o(1). The sum over the at most q^2 pairs (h, k) is o(1).

**Completeness of the case split.** The cases cover every nonzero integer v: unbalanced (all h); balanced h <= H; balanced H < h < q/2 nonprimitive; balanced H < h < q/2 primitive; balanced h >= q/2 noncharacters; characters (offsets, other block). No character is balanced with h < q/2: pairs have l^1 norm q + 1 on I, and the others have norm at least 1.2q.

## 5. Normalisation and Haar measure (Proved; Computed)

**Normalisation.**
- Subtracting the column minimum and replacing e_j by the range changes 2a - e by the row-independent column constant 2 - max - min. That is a common phase.
- W does not increase, so the arc Delta = e^{-W/4}/10 is still at most e^{-W'/4}/10.
- w(v) is unchanged for words on the retained columns. Row differences are unchanged there and vanish on deleted columns, so the character set is unchanged under zero extension.
- Pair l^1 distance q + 1 on I is preserved, because deleted columns contribute 0. So at least (q+1)/2 active primes remain.
- In the computed instances no column needed deletion.

**Haar measure.** Let T = {x in R^I : K_{S,I} x in R 1}, the tangent space of the realisation component t + T, and let V = span(1_I, columns of K_{S^c,I}^T).
- *V contained in T.* Shown in §3.
- *V^perp contained in T^perp.* Let v be in V^perp, extended by zero. Then sum v = 0 and (Kv)_{S^c} = 0, so qv = K_S^T (Kv)_S, with (Kv)_S of sum 0. So v lies in the zero-sum row span, which is T^perp. Hence V = T.
- *Computation (`haar2.py`).* At q = 43, 59, 103, 127, containment holds over Z, and rank_p V = |I| - rank_p(K_S row differences) at p = 2^61 - 1. That forces V = T over Q.
- *Absolute continuity.* (c_0, zeta) -> c_0 1 + B^{-J-2} K_{S^c}^T zeta is a linear surjection onto T, applied to a box law. So given everything else, theta has an absolutely continuous law on t + T. A positive-probability good event therefore has positive Haar measure.
- *After normalisation.* Projection onto the retained columns maps T onto the new tangent space T': extend by zero on the deleted columns, which are constant on S. This was also computed.
- *Non-arithmetic angles.* Removing countably many arithmetic angle vectors keeps positive measure. Correct.

**Weight and height.** W <= 2(q-1) log(2P) = (2d + o(1)) q log q, and m = (rho/2 + o(1)) q, so W/(m log m) <= 4d/rho + o_rho(1). Correct.

## 6. The range statement (VIII.10) (Proved)

**Re-derivation.** Each of the following conditions is forced by a specific estimate:
1. b > 1.5 + a + kappa (+epsilon) is needed for HE/B -> 0. Here E = q^{1.5+a+epsilon} bounds q A sqrt(q log q) and the lifts.
2. b < d is needed for (P/B)^{h/2} >= 4 on short balanced words, for the zeta tail (bq/4 - dh/2 < 0 when h >= q/2), and for the nonprimitive count, whose per-token exponent (1-kappa)/2 + b/4 - d/2 < xi/4 - kappa/2.
3. a > tau is needed so that exp(-c(A/T_0)^2) with T_0 <= q^tau is negligible.
4. tau < 1/2 is needed for the support and norm union bounds: (2q+1)^{q^tau} against rho^{sqrt q} and exp(-c q^{1-tau}).
5. In the strong range, at nu = 1/2 the exponent is (xi + kappa + 2 epsilon)/2 - 3 gamma/4. At gamma = kappa this is negative if and only if kappa > 2 xi + 4 epsilon.
6. In the weak range it is (xi + kappa + 2 epsilon)/2 - gamma/2. At gamma = tau this is negative if and only if tau > xi + kappa + 2 epsilon.
7. The unbalanced words need xi/2 < lambda < kappa/2. The per-token exponent is xi/2 - lambda for h > H_u, and the digit series is O(H_u E/B) = O(q^{lambda - kappa - epsilon}) = o(1/H_u) for h <= H_u.

The sympy and fraction computation confirms the coefficients of nu: 2 epsilon + kappa + gamma/2 and 2 epsilon + kappa + gamma.

**Closure.** Combining the conditions gives 1.5 + a + kappa < 2 - xi with a > kappa + xi and kappa > 2 xi. At the lower endpoints this reads 1.5 + 6 xi < 2, so xi < 1/12. This threshold is exact for this system. At xi = 1/12 - 1/1000 an explicit feasible point exists, with d - b = 0.0056.

**Changed estimates.**
- The permutation error is q^{1.5-d} (log q)^{3/2} = o(1).
- 1.01/d <= 0.527 < 0.53.
- The middle case of (VIII.9) gives 4/3 - 0.53/2 = 1.0683 > 1.06.
- The case s >= q/2 gives (16/9 - 0.53)/2 + 4/9 = 1.0683.
- The second bound becomes 1 - 0.53/2 = 147/200 >= 0.73.

**Limits.** d -> 23/12, and 4d/rho -> 23/3 as xi -> 1/12 and rho -> 1 along a diagonal sequence. For 0 < delta < 1/12, choose xi in (delta, 1/12). This gives fakes in the window P_+ <= m^{2-delta}.

**Dependency.** VIII.10 uses (U) and the gcd property up to h <= q^tau, with tau as large as about 1/4. As written, (U) holds for every h with sqrt(q) - h - 1 > 0. Its referee should confirm this range.

## 7. The count over F_11 (Computed)

`newton.py` enumerates ordered pairs (P, N) of disjoint 3-multisets over all of F_11, by the number r of initial moments p_1, ..., p_r that agree. It gives 40370, 3520, 440, 0 for r = 0, 1, 2, 3, which reproduces Astra's numbers. Her bounds binom(13, 3) binom(13 - r, 3 - r) are 81796, 18876, 3146, 286.

Her enumeration includes the label 0. The construction uses only F_q \ {0}, where the counts are 22530, 1980, 240, 0. This difference does not matter.

## Assessment

The faketail block is correct. Nothing in it needs repair, and no exponent needs to worsen. The fake construction of (VIII.8) and (VIII.10) stands or falls with the first half of VIII.3: the uncertainty estimate (U), including its mod-4 form and its range in h, the offset and character estimates, and the ordered residue witnesses.
