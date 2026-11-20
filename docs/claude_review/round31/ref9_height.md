# Referee report, block "height": Part IX.4 (deeper window 11/6, coefficient 22/3) and its "PROVED (quantifiers)" paragraph

Scope: `astra_free_IX/answer.md`, section IX.4 (lines 2087-2136). I also read the corresponding draft in `ix_atom_work.md`, which uses a different residue argument (see item 6). Scripts I wrote are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref9_height/`:
- `exponents.py`: sympy and exact rationals;
- `lp.py`: linear programme for where the parameter system closes;
- `atoms.py`: floating-point exact convolution of discrete Gaussians at q = 103 and 199.

I did not import Astra's scripts. I did no literature search.

This block depends on four inputs that I took as given:
- VIII.3 as refereed, with the corrections of `README_VIII.md` §§1.3-1.4. In particular, the uncertainty estimate is: for primitive balanced v and l != q, at least sqrt q + 1 - h coordinates of Kv are nonzero mod l, for each odd prime l <= h and for l = 4.
- VIII.4, the spread variant, as refereed in `ref8_spread.md`.
- The VIII.10 proof, refereed in `ref8_faketail.md` §6.
- IX.3: the offsets, the factored witnesses and inequality (IX.7). Another referee checks these. IX.4 inherits whatever that referee rules about the P3-char property.

## Results

**Verdict.** The new atom bound is correct, and so is the parameter system with its closure at xi = 1/6. Every changed estimate is correct. The conclusions limsup log P_+/log m <= 11/6 and W/(m log m) -> 22/3 (on a diagonal) are **Proved**, conditional on IX.3. I found no gap. The quantifier paragraph needs three small corrections (item 7):
- the threshold q_0 also depends on the arc constant C;
- q runs over primes q = 3 mod 4;
- "11/6" and "22/3" are values along a diagonal sequence, while at fixed (xi, rho) the values are 2 - xi and 4(2 - xi)/rho.

| # | Statement | Label |
|---|---|---|
| 1 | Support property: for every primitive balanced v with h = \|\|v\|\|_1 <= q^tau (tau < 1/2) and n = 4 or an odd prime n <= h, at least 2 c_rho sqrt q of the w_i = (Kv)_i, i in S^c, are nonzero mod n; also \|\|w\|\|_2^2 >= c_rho q h | Proved (union bound closes exactly when q^tau log q = o(sqrt q)) |
| 2 | Atom bound sup_b P(sum w_i X_i = b) <= q^{o(1)} (A sqrt(qh))^{-1}, for A = q^a with any fixed a > 0 and h <= q^tau | Proved; sanity check Computed (floating point) |
| 3 | The support property covers every modulus d_0 > 2 that occurs | Proved (reduction to an odd prime factor of d_0, or to 4) |
| 4 | Parameter system: feasible if and only if xi < 1/6 (closed form: sup xi = 1/6); the old system with a > tau closes at 1/12 | Proved; Computed (LP and exact rationals) |
| 5 | Exponents per token (xi+kappa+2eps)/2 - 3gamma/4 and (xi+kappa+2eps-gamma)/2; coefficient of nu positive; all unchanged ranges | Proved; recomputed with sympy |
| 6 | Residues: (IX.7) does not depend on d. The draft variant with b_0 <= 1.01/d < .56 also closes | Proved (given IX.3 and VIII.9) |
| 7 | Conclusion and quantifier paragraph | CORRECT WITH CORRECTIONS (corrected statement below) |

### Corrected statement (quantifiers)

**Theorem (IX.4, corrected wording).** Fix C > 0, 0 < xi < 1/6 and rho in (0.99, 1). Fix auxiliary parameters chosen as explicit functions of xi:

    kappa = 2xi + 5eps,  tau = kappa + xi + 3eps,  a = eps,  lambda = (xi + kappa)/4,

with eps = eps(xi) > 0 so small that tau < 1/2 and 1.5 + a + kappa + 2eps < 2 - xi. (Such an eps exists exactly when xi < 1/6.)

Then there is q_0 = q_0(C, xi, rho) such that for every prime q = 3 mod 4 with q >= q_0, the construction of IX.3 + VIII.3 has a realisation with positive probability, and with positive Haar measure on its realisation component. This realisation is a P3-char fake in the sense of IX.3, with:
- m = (rho/2 + o(1)) q marks;
- P_+(N') <= 2 q^{2-xi};
- W/(m log m) <= 4(2 - xi)/rho + o_{C,xi,rho}(1).

Consequently, for each fixed 0 < delta < 1/6 and each C, these fakes eventually lie in the window P_+(N') <= m^{2-delta}: choose xi in (delta, 1/6), then use m^{2-delta} >> q^{2-xi}. Along a diagonal sequence xi_k -> 1/6, rho_k -> 1, q_k >= q_0(C, xi_k, rho_k), one gets limsup log P_+/log m <= 11/6 and limsup W/(m log m) <= 22/3. No threshold is uniform in xi, rho or C. The construction does not reach delta = 1/6, and the limsup statement does not give the endpoint.

The paragraph as written says "take q sufficiently large depending on both [xi, rho]". That omits C, which enters through Delta = min(1, C) e^{-W/4}/10, the sector width b_* = min(1, C)/1000 and the radial bins sigma = O_C(q log q/P) of IX.3. It also omits that q is a Paley prime. Otherwise it is correct.

## 1. The support property (Proved)

**Before sampling.** By the corrected uncertainty estimate, at least sqrt q + 1 - h >= sqrt q (1 - o(1)) coordinates of Kv are nonzero mod n. This holds for n = 4 and for each odd prime n <= h; note n <= h < q, so n != q. Astra's sqrt q - h - 1 is weaker and suffices.

**Sampling.** Rows outside the class A are always in S^c. Each row of A is independently in S^c with probability 1 - rho. So the number of detecting indices in S^c stochastically dominates Bin(n_det, 1 - rho), with n_det >= sqrt q/2. Chernoff's lower tail gives

    P(count < (1 - rho) sqrt q/5) <= exp(-(1 - rho) sqrt q/16)

for large q. This is the required 2 c_rho sqrt q with c_rho = (1 - rho)/10.

**Union bound.** The vectors with \|\|v\|\|_1 <= q^tau number at most (2q+1)^{q^tau+1}, and there are at most q^tau moduli. The logarithm of the count is O(q^tau log q). This is o(sqrt q) if and only if tau < 1/2. So exp(-Omega_rho(sqrt q)) beats the count exactly in the stated range.

**Norm.** The VIII.3 norm test gives \|\|w\|\|_2^2 >= (1 - rho) q \|\|v\|\|_2^2/2 >= (1 - rho) q h/2 >= c_rho q h. Its failure probability is exp(-c_rho q/h). This also beats (2q+1)^{h+1} for h <= q^tau, tau < 1/2, since q/h >> h log q. The events are combined with the VIII.3 and IX.3 events, all of probability 1 - o(1).

## 2. The atom bound (Proved)

Notation:
- X_i are independent, with P(X_i = n) proportional to exp(-n^2/(2A^2));
- A = q^a;
- Z = (log q)/A;
- F(t) = sum_i \|\|w_i t\|\|^2.

I checked each step.

1. **Characteristic function.** Poisson summation gives psi_A(u) = sum_k e^{-2 pi^2 A^2 (u-k)^2} / sum_k e^{-2 pi^2 A^2 k^2}.
   - On \|\|u\|\| <= 1/4 the nearest term dominates up to a factor 1 + O(e^{-cA^2}).
   - On 1/4 < \|\|u\|\| <= 1/2 one has psi <= 2 e^{-2 pi^2 A^2 \|\|u\|\|^2} <= e^{-c A^2 \|\|u\|\|^2} once A is large, after reducing c.
   - The product over at most q factors of (1 + O(e^{-cA^2})) is exp(O(q e^{-c q^{2a}})) = 1 + o(1), for any a > 0.
   - Off {F <= Z^2} the integrand is at most exp(-c log^2 q). This is superpolynomially small against the target q^{-a-1/2-tau/2}.
2. **Excluded indices.** Each index with \|\|w_i t\|\| >= 1/(3h) contributes at least 1/(9h^2) to F. So there are at most 9 h^2 Z^2 = O(q^{2tau - 2a} log^2 q) such indices, which is o(q). The nonzero w_i number >> c_rho q, by the degree-(q-1)/2 polynomial support mod q of balanced v with h < q/2. So >>_rho q indices are good.
3. **Common ratio.** Write w_i t = k_i + delta_i with \|delta_i\| < 1/(3h). Then

       k_i w_j - k_j w_i = w_i delta_j - w_j delta_i,

   which has absolute value < 2/3 because \|w_i\| <= \|\|v\|\|_1 = h. Being an integer, it vanishes. So r = k_i/w_i is common to all good indices. Write r = b_0/d_0 in lowest terms. Then w_i b_0/d_0 is an integer, so d_0 divides every good w_i.
4. **|t - r|.** F(t) >= (t - r)^2 sum_good w_i^2 >= (t - r)^2 d_0^2 c q. Hence |t - r| <= Z/(d_0 sqrt(cq)).
5. **F(r) <= 2Z^2.** Take an index with d_0 not dividing w_i. Then \|\|w_i r\|\| >= 1/d_0, while \|w_i (t - r)\| <= hZ/(d_0 sqrt(cq)). This is o(1/d_0) because h Z/sqrt q = q^{tau - a - 1/2} log q -> 0. So \|\|w_i r\|\| <= (1 + o(1)) \|\|w_i t\|\|. Indices with d_0 dividing w_i contribute 0. Hence F(r) <= (1 + o(1)) F(t).
6. **d_0 = 1, 2.** If d_0 = 2 and some w_i is odd, then F(r) >= 1/4 > 2Z^2, which is impossible. So all w_i are divisible by d_0. Then \|\|w_i t\|\| = \|w_i (t - r)\| because h\|t - r\| -> 0. Hence F(t) = (t - r)^2 \|\|w\|\|_2^2, and each of the at most three centres 0, 1/2, 1 contributes length <= 2Z/\|\|w\|\|_2 = O_rho(Z/sqrt(qh)). (All w_i even is possible: (Kv)_i = v_i mod 2.)
7. **d_0 > 2.** Here d_0 divides a nonzero w_j, so d_0 <= h.
   - *Covering of the modulus (the question raised).* The support property is stated for n = 4 and odd primes n <= h. It covers d_0. Either d_0 has an odd prime factor p, with p <= d_0 <= h, and w_i != 0 mod p implies w_i != 0 mod d_0. Or d_0 = 2^k with k >= 2, and w_i != 0 mod 4 implies w_i != 0 mod d_0. So N_0 >= 2 c_rho sqrt q. Correct.
   - *Uniformity.* Let g_i = gcd(w_i, d_0) and d_i = d_0/g_i > 1. Then w_i b_0/d_0 = (w_i/g_i) b_0/d_i with w_i/g_i a unit mod d_i. Reduction (Z/d_0)^* -> (Z/d_i)^* is onto with equal fibres. So the image of a uniform unit b_0 is uniform on the reduced fractions of denominator d_i. Correct.
   - *Hölder.* For admissible b_0, prod_{i <= N_0} exp(-\|\|w_i b_0/d_0\|\|^2/Z^2) >= exp(-F(b_0/d_0)/Z^2) >= e^{-2}. Hence #admissible <= e^2 phi(d_0) avg_b prod_i f_i(b). Generalised Hölder with N_0 exponents equal to N_0 bounds the average by the geometric mean of avg f_i^{N_0}. Correct.
   - *Average.* avg f_i^{N_0} is (1/phi(d_i)) sum over k coprime to d_i of exp(-N_0 \|\|k/d_i\|\|^2/Z^2). Since d_i > 1, k = 0 is not a reduced residue. So this is at most (2/phi(d_i)) sum_{k=1}^{floor(d_i/2)} exp(-N_0 k^2/(d_i Z)^2). That is at most (2/phi(d_i)) times the integral from 0 to infinity, which equals sqrt(pi) (d_i/phi(d_i)) Z/sqrt(N_0) <= q^{o(1)} Z q^{-1/4}. The absence of k = 0 is exactly what makes the integral comparison valid. Correct. (When tau < 1/4 + a the sum is even exponentially small. The stated bound holds in all cases.)
   - *Number of centres.* At most e^2 phi(d_0) q^{o(1)} Z q^{-1/4} <= d_0 q^{o(1)} Z q^{-1/4}. Their total length at radius O(Z/(d_0 sqrt q)) is q^{o(1)} Z^2 q^{-3/4}, independently of d_0.
   - *Anchors.* Let E = {t : F(t) <= Z^2, d_0(t) > 2} and G(t) be the set of good indices. Then

         |E| <= (cq)^{-1} sum_j |{t in E : j in G(t)}|.

     For t in the j-th set, d_0(t) divides w_j, so only tau(\|w_j\|) = q^{o(1)} denominators occur. This gives |E| <= q^{o(1)} Z^2 q^{-3/4}. The anchor step is needed. Summing over all d_0 <= h instead would give h Z^2 q^{-3/4}, which is absorbed only for tau < (1/4 + a)(2/3).
8. **Absorption.** The ratio (q^{o(1)} Z^2 q^{-3/4}) / (Z/sqrt(qh)) equals q^{o(1)} Z h^{1/2} q^{-1/4} <= q^{-a + tau/2 - 1/4 + o(1)}, which tends to 0 since tau < 1/2. The total is O(log q/(A sqrt(qh))) + o(1/(A sqrt(qh))) = q^{o(1)} (A sqrt(qh))^{-1}. The q^{o(1)} is uniform in v: it comes from log q, the divisor bound for integers <= q^tau, and d/phi(d) << log log d.
9. **No hypothesis a > tau.** The only conditions used are 2tau - 2a < 1 and tau - a - 1/2 < 0, both implied by tau < 1/2. Correct. This is the improvement: the old atom-2 error exp(-c(A/T_0)^2) required A >> T_0.

**Computed (`atoms.py`, floating point, exact convolution; a sanity check only).** Setup:
- q = 103 and 199, rho = 0.9, random S;
- 6 random primitive balanced v for each h in {2, 4, 6, 8, 12};
- A in {0.6, 1, 2, 4}.

Results:
- The maximum of atom * A * \|\|w\|\|_2 * sqrt(2 pi) was at most 1.36 (attained at A = 0.6), and 1.02-1.13 for h >= 6.
- atom * A * sqrt(qh) was at most 0.74.
- Counts of nonzero coordinates of w mod n (n = 4, 3, 5, 7, 11) ranged from 22 to 95. This is consistent with the support property, although the asymptotic hypotheses (sqrt q > h, q^tau log q << sqrt q) are vacuous at these q.
- An adversarial non-Paley control illustrates why N_0 is needed: w = (3, ..., 3, 1^{n_1}) with n_1 = 0 gives ratio 3, the gcd, and the ratio decreases to 1 as n_1 grows.

## 3. The parameter system (Proved; Computed)

The new conditions are:
- kappa > 2xi + 4eps (strong range at gamma = kappa);
- tau > kappa + xi + 2eps (weak range at gamma = tau);
- tau < 1/2 (support and norm union bounds);
- a > 0;
- b = 1.5 + a + kappa + 2eps < d = 2 - xi (HE/B -> 0 with E = q^{1.5+a+eps}; short words; zeta tail; nonprimitive count);
- xi/2 < lambda < kappa/2 (unbalanced words).

At the lower endpoints, tau < 1/2 reads 3xi < 1/2, and b < d reads 1.5 + 2xi < 2 - xi, that is 3xi < 1/2. Both constraints bind at the same point.

- `lp.py`: the closed system has sup xi = 1/6, attained at (xi, kappa, tau, a) = (1/6, 1/3, 1/2, 0). Adding the old constraint a >= tau gives 1/12, which reproduces ref8_faketail §6.
- `exponents.py`, exact rationals: an explicit feasible point exists at xi = 1/6 - 10^{-5}; at xi = 1/6 the conditions tau < 1/2 and b < d both fail.

So 1/6 is exactly where this system closes. Any further gain would need both of the following, since each alone is binding:
- lifts with sup norm below q^{3/2};
- a support (uncertainty) property for h >= sqrt q.

**Exponents per token (sympy).** For primitive balanced words with count (1 - gamma)(1 - nu), digit density b nu and atom cost alpha per erased digit:
- strong range (alpha = a + 1/2 + gamma/2): the coefficient of nu is 2eps + kappa + gamma/2 > 0, and F(1/2) = (xi + kappa + 2eps)/2 - 3gamma/4;
- weak range (alpha = a + 1/2): the coefficient of nu is 2eps + kappa + gamma > 0, and F(1/2) = (xi + kappa + 2eps - gamma)/2.

These are Astra's expressions. Their maxima over the ranges gamma in (kappa, tau] and (tau, 1) are xi/2 - kappa/4 + eps < 0 and (xi + kappa + 2eps - tau)/2 < 0. These are fixed negative numbers, so q^{o(1)} per atom, exp(O(h)) and H^2 are absorbed, since h > q^kappa. The sum over the O(q^2) pairs (h, k) is o(1).

**Unchanged ranges.**
- Nonprimitive words: (1 - kappa)/2 + b/4 - d/2 = a/4 + eps/2 - kappa/4 + xi/2 - 1/8 < -1/8 + a/4 + eps < 0.
- Unbalanced words, h > H_u: the exponent is xi/2 - lambda < 0. For h <= H_u, H_u E/B = q^{lambda - kappa - eps} = o(1/H_u), since 2 lambda < kappa.
- Short balanced words: hE/B = q^{-eps}.
- zeta tail: (b - d) q/4 < 0.
- Lifts O(q^{1.5} sqrt(log q)) and q A sqrt(q log q) are both <= E for every a, eps > 0. The digit sup-norm event has failure probability q^2 exp(-c q^{2eps}).
- Atom 1, used for h > q^tau: (1 + O(e^{-cA^2}))^q = 1 + o(1) for any a > 0.

**Weighted-bin error.** sigma = O_C(q log q/P) from the IX.3 bins. The deviation is O(sigma sqrt(q log q)) = O_C(q^{1.5-d} log^{3/2} q) = o(1), since d > 11/6 > 3/2. Correct.

**Spread variant.** |T| = O(q^{xi/2} + 1) = O(q^{1/12}) = o(sqrt q). So 3^{|T|} rho^{sqrt q - |T| - 1} -> 0, and the spread tails (per-token exponent xi/2 - kappa < 0, and factors q^{O(n_0)}, B^{n_0}) are absorbed. Correct.

## 4. Residues (Proved, given IX.3)

answer.md uses (IX.7), h - b_0 n >= (1 + a_rho) q, with b_0 = log L_0/log P_0 = 1 + o(1). This needs only d > 3/2 and is independent of the digits. The residue part therefore does not constrain xi.

The draft in `ix_atom_work.md` instead keeps the VIII auxiliary primes l_i <= q^{1.01}. There b_0 <= 1.01/d <= 0.5509 < .56 at d >= 11/6. The cases give:
- 4/3 - .28 = 1.0533;
- (16/9 - .56)/2 + 4/9 = 1.0533;
- h - .56 s >= .72h.

This is enough, because the requirement is h - b_0 s >= (1 + c) q against (L_0/Delta) L_0^{s/2} P^{-h/2} exp(O(h)) with 1/Delta << (2P)^{q/2}. The threshold 1.06 written in VIII.9 is not needed. Both versions close.

## 5. Conclusion (Proved; quantifiers corrected above)

- P_+ <= 2P = 2q^{2-xi} and m = (rho/2 + o(1)) q give log P_+/log m -> 2 - xi.
- W <= 2(q - 1) log(2P) gives W/(m log m) <= 4(2 - xi)/rho + o(1).
- On the diagonal the limits are 11/6 and 22/3.
- The window P_+ <= m^{2-delta} is capped for every fixed 0 < delta < 1/6, for each fixed C. It is not capped at delta = 1/6.

**Novelty.** The atom bound is a Halász-type small-ball estimate via rational approximation of the frequency. Its structure (common nearest-integer ratio, denominators, Hölder over reduced residues) is standard in style. I did no search. Its use here is correct and improves the record's window from 1/12 to 1/6. This is a quantitative improvement of a fake, not an integer result.
