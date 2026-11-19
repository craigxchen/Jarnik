# Referee report, round 31, Part VI, block "code" (section VI.1, the derivative code condition)

Scope: VI.1 (the theorem, the packing bound, the application to the skew double core, the two-Walsh-core
low-rank test), a check on the genuine Part V clusters, and the integer and abc analogue. I did not referee
VI.2 to VI.5.

Scripts: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref6_code/`.
All are my own code. `gf2.py` holds the F_2[t] arithmetic, factorisation and F_4[t] arithmetic. The
factorisation was cross-checked against sympy on 60 random polynomials (`test_gf2.py`). The other scripts are
`partV_code.py`, `profiles.py`, `equalweight.py`, `walsh_two.py`, `pell.py` and `abc_chars.py`, with `*.out`
logs. Each run took under a minute and used less than 3.5 GB.

## Verdicts

| claim | verdict |
|---|---|
| VI.1 weighted-distance theorem (every nonzero word of the anchored parity code has weight >= K+1 > W/4) | **CORRECT** (line by line). Minor sharpening: weight >= 2*ceil((K+1)/2) + deg H'. |
| VI.1 packing bound r <= (1-H_2(1/8)) n + O(log n), and the m = 167 numerics | **CORRECT** (reproduced, plus 191 and 251) |
| Application: the skew double core is excluded over F_2[t] | **CORRECT**. With arbitrary weights the certificate works for q >= 151 (first prime q = 3 mod 4). With equal weights the profile is excluded already at q = 19, 23. |
| Low-rank test: two Walsh cores on x, x^3; direct sum; e >= 4d-2; W = O(m log m) still allowed | **CORRECT** |
| (d) Part V genuine clusters satisfy the condition | **CORRECT** (consistency check, m = 4, 8, 16) |
| (e) Integer analogue of the parity condition | **FALSE over Z**, already for 2 points (exact Pell family). The abc analogue is a Z-support (radical) condition. That is the Lang–Waldschmidt-shaped support test, which the double core passes (record round8/lw.md Thm 5.1). No growth improvement follows. |

---------------------------------------------------------------------------------------------------------

## (a) The proof of VI.1, line by line

**Setting.** Points z_0..z_{m-1} lie in F_4[t], have degree D, and satisfy z_i zbar_i = N in F_2[t] and
deg(z_i - z_j) <= D/2. Then W = deg N = 2D.

1. **Monic points.** Every u in F_4^* has norm u^3 = 1, so the leading coefficients are not fixed by the norm.
   Since deg(z_i - z_j) <= D/2 < D (D > 0), all leading coefficients are equal, and after one common unit
   scaling all points are monic. With monic pi_p, the conjugate pibar_p is also monic. Hence
   z_i = prod_p pi_p^{a_ip} pibar_p^{e_p - a_ip} * prod_{q inert} q^{f_q/2} with no unit factor. The extension is a
   constant-field extension, so no prime ramifies. A prime p of F_2[t] splits iff deg p is even. OK.
2. **K.** deg(z_i - z_0) <= floor(D/2), so ord_inf((z_i - z_0)/z_0) >= D - floor(D/2) = ceil(D/2) = K, and
   R_i in U_K = 1 + t^{-K} F_4[[1/t]]. U_K is a group, since (1+x)^{-1} = 1 - x + x^2 - ... OK.
3. **Logarithmic derivative.** Write R = 1+u with ord u >= K. Then R'/R = u'/(1+u), and d/dt(t^{-j}) = -j t^{-j-1}, so
   ord(R'/R) >= K+1, or R'/R = 0. *Characteristic-2 subtlety, a sharpening:* only odd j survive. So u' lies in
   t^{-2}F_4[[t^{-2}]] and ord(R'/R) is **even**. For K even this gives ord >= K+2, and in general
   ord >= 2*ceil((K+1)/2). In the Part V data (K = 6, 56, 360, all even) the minimum of deg H - deg H' over all
   words is exactly K+2 (below). No asymptotic effect.
4. **Identity.** (pi/pibar)'/(pi/pibar) = pi'/pi - pibar'/pibar = pi'/pi + pibar'/pibar = p'/p in characteristic 2.
   Integer coefficients act through their parity. OK.
5. **Reduction to H.** For c in {0,1}^{m-1}, R = prod R_i^{c_i} and R'/R = sum_p v_p p'/p, where
   v_p = sum_i c_i (a_ip - a_0p). Inert primes cancel (ord_q z_i = f_q/2 for all i). Only v_p mod 2 matters, so
   R'/R = H'/H with H = prod{p : v_p odd}. The coordinate of the codeword at p is the pattern of column p
   evaluated at c. **Bookkeeping:** two primes with the same anchored pattern take the same value on every word,
   so they always enter H together. Zero patterns never enter. So the weight of a pattern is the sum of the
   degrees of its primes, and deg H is the weighted weight of the word. Correct.
6. H is squarefree over the perfect field F_2. If H is nonconstant it is not a square, so H' != 0. Then
   ord_inf(H'/H) = deg H - deg H' <= deg H. With step 3 this gives deg H >= K+1 + deg H' >= K+1 > D/2 = W/4. OK.
   The statement also covers prime powers (parities of a_ip - a_0p) and arbitrary degrees.
7. **Version with c > 0.** With deg(z_i - z_j) <= D/2 + c, K drops to ceil(D/2) - c. The test still bites
   asymptotically if c = o(W).

**Packing.** The total weight of all patterns is at most sum_split e_p deg p <= W. Balls of weighted radius W/8
around codewords are disjoint, because two codewords in one ball would be at distance <= W/4 < K+1. Take a
uniform j-subset. Its expected weight is <= jW/n, so by Markov a fraction >= 1 - 8j/n of the j-subsets has
weight <= W/8. Hence 2^r (n - 8j) C(n,j) <= n 2^n for every j < n/8. With j = n/8 - 1 and Stirling this gives
r <= (1 - H_2(1/8)) n + O(log n). Here 1 - H_2(1/8) = 0.4564 < 0.457. **CORRECT.**

**Numerics, reproduced** (`profiles.py`, independent Paley construction; the skew-Hadamard property K K^T = M I
was asserted for q < 120). The table gives the smallest j that works, found by searching j upward.

| profile | rows | distinct nonzero anchored patterns n | anchored binary rank | certificate j |
|---|---|---|---|---|
| DCo(167) (ordinary rows) | 167 | 334 | 166 | 40 (Astra's) |
| DCo(191) | 191 | 382 | 190 | 46 |
| DCo(251) | 251 | 502 | 250 | 59 |
| DC(q) with border row, same q | q+1 | 2q | q | same j |

All 2q anchored columns are distinct, so Astra's script's use of n = 2q without deduplication is harmless. Its
script `vi_profile_check.py` imports `ii_paley` from `astra_free_III`, so it does not run from its own
directory. My code is independent of it.

## (b) The application, and which record fakes are excluded over F_2[t]

**The profile identified.** The double core is round8/lw.md §5.1 (exact.md Thm 4.1): K = I + S Paley
skew-Hadamard of order M = q+1, columns s_a^(1) = K(.,a) and s_a^(2) = K(.,a) with rows 0 and a negated, one
prime per column, e = 1. Astra's Part III.2 uses the q ordinary rows (border deleted). Deleting a row of a
cluster leaves a cluster, so excluding the q-row profile also excludes the full one. The full one is excluded
anyway: rank q on 2q columns.

**Rank.** The two columns of a label differ only at row a ("twins"). So e_a lies in the binary span, the
unanchored rank is m, and the anchored rank is m-1. Reproduced exactly for every q tested.

**Exclusion thresholds** (any weights, Markov certificate):
* double core: first certificate at q = 151;
* fully flipped core FF(q,3) of construct.md §2: first certificate at q = 139 (M = 140);
* single Paley core, q = 3 mod 8: rank q-1 on q columns, certificate already at q = 11;
* single Paley core, q = 7 mod 8: rank (q-1)/2 on q columns, first certificate at q = 311.

**Equal weights** (one prime of a common even degree per column, the natural F_2[t] analogue of primes in
[P, P e^{1/(10r)}]). The test becomes: minimum Hamming weight > (#columns)/4. Exhaustive enumeration in
`equalweight.py`:

| q | 7 | 11 | 19 | 23 |
|---|---|---|---|---|
| min Hamming weight | 4 | 8 | 8 | 8 |
| needed | > 3.5 | > 5.5 | > 9.5 | > 11.5 |
| result | passes | passes | **excluded** | **excluded** |

**Which of the record's families are excluded** (over F_2[t] only; profile level, at any degree):
* **Excluded for large M:** every profile with the twin mechanism. This covers the Paley skew double core
  (lw.md §5, exact.md §4) and the fully flipped Hadamard cores of construct.md §2 and lower.md §3, for any number
  b of copies. Their unflipped copies merge into one pattern, so n <= (M-1) + M and the rank is M-1. The ratio
  is about 1/2 > 0.4564. The round-5 residue designs are built on these profiles, so their profiles are excluded
  too. These are all the record's growth-scale integer fake profiles with integral characters.
* **Pass:**
  * Sylvester/Walsh cores: rank t on 2^t - 1 columns, Hamming weight 2^{t-1} > n/4. This is exactly the
    genuine Part V shape.
  * t >= 3 row-permuted Paley bases (lw.md §7): rank <= M-1 on about tq columns, ratio <= 1/3 < 0.457. The rank
    bound does not exclude them, and I did not compute their weighted distance.
  * Low-rank two-core profiles (part (c)).
* **Not checked:** grafts (round10/analytic.md, round20/composite.md).

Astra's remark that "more distinct patterns or a smaller parity rank can pass" is accurate.

## (c) The low-rank test (two Walsh cores on x and x^3)

* **Direct sum.** If a functional (A,B) vanishes on all (x, x^3), then Tr(Ax + Bx^3) = 0. Polarising gives
  Tr(x[(By)^{2^{d-1}} + By^2]) = 0. I checked this identity numerically on 200 random triples for d = 3, 5, 7.
  The nondegenerate trace then gives (By)^{2^{d-1}} = By^2. Squaring gives By = B^2 y^4, a polynomial of degree
  4 vanishing on more than 4 points, so B = 0 and then A = 0. A shorter route: the exponents 2^i and 3*2^i are
  distinct mod 2^d - 1 for d >= 3. Cubing is a bijection iff gcd(3, 2^d - 1) = 1, i.e. d odd.
* **Computed** (`walsh_two.py`, d = 3, 5, 7): anchored rank 2d, and the nonzero words have
  (core-1 weight, core-2 weight) in {(m/2, 0), (0, m/2), (m/2, m/2)}. With W = 2(m-1)e and K = (m-1)e/2, the
  weighted weight m e/2 >= K+1 holds for e >= 2. So the profile passes.
* **e >= 4d-2.** Each word u on core 1 alone has deg H_u = me/2, so deg H_u' <= e/2 - 1. Distinctness: if
  H_u' = H_v', then H_u + H_v is a nonzero square of degree < me/2 (the H's are monic). It is divisible by the
  squares of the m/4 common primes, which have total degree me/2, a contradiction. In characteristic 2 every
  derivative is a square. The number of nonzero squares of degree <= e/2 - 1 is 2^{ceil(e/4)} - 1, checked for
  even e <= 22. So m - 1 <= 2^{ceil(e/4)} - 1, which gives e >= 4d-2 exactly (the least even e for d = 3, 5, 7, 9
  is 10, 18, 26, 34).
* The sharpening of step 3 does not change this: e = 4d-2 is 2 mod 4, so K is odd. W >= (2m-2)(4d-2) = O(m log m).
  **CORRECT**, and it is only a necessary condition.

## (d) The genuine Part V clusters (`partV_code.py`)

I built the clusters myself: z_H = prod_{v in H\0}(t^L + theta + v), L = (d-2)2^{d-1}. For each I checked
z zbar = N for every point and max deg(z_i - z_j) = D/2. I factored N over F_2 (N has no inert prime factors),
took pi = gcd_{F_4}(p, Y+v) at each split p, checked a + abar = e_p, and enumerated all 2^{m-1} words.

| m | D | K | split base primes | distinct patterns n | parity rank | min weighted weight | K+1 | min over words of deg H - deg H' |
|---|---|---|---|---|---|---|---|---|
| 4 | 12 | 6 | 5 | 2 | 2 | 8 | 7 | 8 = K+2 |
| 8 | 112 | 56 | 16 | 6 | 3 | 92 | 57 | 58 = K+2 |
| 16 | 720 | 360 | 52 | 13 | 4 | 556 | 361 | 362 = K+2 |

Every word satisfies deg H - deg H' >= K+1 (asserted for each word). The parity code has rank log_2 m: low rank,
like the Walsh cores, and the condition is far from tight in raw weight. The bound on deg H - deg H' is attained
at K+2, which confirms the parity sharpening in (a)3.

## (e) The integer analogue, and what abc gives

**(e1) The parity statement is false over Z.** A 2-point exact family (`pell.py`): pi = 2+i,
gamma = x+iy with (x+2y)^2 - 5y^2 = +-1, z_1 = gamma^2 pi, z_0 = conj(z_1). Then z_1 - z_0 = 2i(x^2 + 4xy - y^2) =
+-2i. The points are primitive and coprime, and N is odd. The anchored parity word is supported only at p = 5
(odd exponent difference at 5, even differences at the primes of gamma). So the odd-support weight is
log 5 = 1.609 for every member, while W = log N grows without bound. The first member is gamma = 4+17i,
N = 5^3 * 61^2, W = 13.05, W/4 = 3.26. The eleventh member has W = 128.5, and its arc constant is
|z_1 - z_0|/sqrt(R) = 2.2e-14.

The same happens mod 3: z_1 = gamma^3 pi with |z_1 - z_0| <= N^{1/4}, for example gamma = 29+36i and
gamma = 460+571i. Mod-3 support weight log 5 against W/4 = 10.3. By Dirichlet, |Im(gamma^l pi)| << |gamma|^{l-2} is
attainable, which is <= N^{1/4} for l <= 4.

So the mechanism behind VI.1 (d/dt kills squares, giving a mod-2 support condition) has no integer counterpart.
This is the char-p Mason phenomenon: p-th powers are "constants".

**(e2) The correct integer statement under abc.** Assume abc over Q(i): for coprime a + b + c = 0 in Z[i],
max|.|^2 <= K_delta N(rad(abc))^{1+delta}. Take a normalised cluster on an arc of length <= C sqrt(R) (so
|z_i/z_0 - 1| <= C N^{-1/4}; the task's "N^{-1/2}" should read N^{-1/4} when N = R^2). Take a zero-sum c in Z^m and
v = sum_i c_i a_i in Z^J, nonzero. Then R = prod z_i^{c_i} = beta/conj(beta), with
beta = eta * prod_{v_p>0} pi_p^{v_p} prod_{v_p<0} pibar_p^{-v_p} and eta in {1, i, 1+-i} absorbing the unit.
Apply abc to beta - conj(beta) = 2i Im(beta), bounding N rad(Im beta) <= |beta - conj(beta)|^2. This gives:

> **(abc-Lambda)** w_0(v) := sum_{v_p != 0} log p >= W/4 - (delta/(2(1+delta))) w_1(v) - log(|c|_1 C) - O_delta(1),
> where w_1(v) = sum_p |v_p| log p.

Liouville alone gives w_1(v) >= W/2 - 2 log(|c|_1 C) - O(1), i.e. |beta|^2 >= N^{1/2}/(|c|_1 C)^2. The function-field
counterpart of (abc-Lambda) is Mason: Z-support >= K+1. VI.1 is the mod-2 strengthening, available only in
characteristic 2. A weak abc with exponent kappa gives w_0 >= W/4 - ((kappa-1)/(2kappa)) w_1 - ...; Szpiro implies
only kappa = 3/2 + eps (via the Frey curve and |abc| >= c^2/2). Vojta's height inequality on P^1 minus 3 points
is equivalent to abc, so it gives the same statement.

**(e3) Consequences, with numbers** (`abc_chars.py`, equal weights, q = 43):
* **Double core: not excluded.** Record lw.md Thm 5.1 gives w_0(v) >= (M/2) min w_j >= W/4 + (1/2) log P - 1/16 for
  every nonzero character. So (abc-Lambda) holds for every delta >= 0, every K_delta >= 1 and every c. abc
  through these relations cannot exclude the double-core profile, at any P, hence at W = Theta(M log M).
* **Fully flipped families FF(q,b), b >= 4: excluded.** The light character c = H(.,a) has |c|_1 = M and:

  | b | w_0/W | w_1/W |
  |---|---|---|
  | 3 | 0.272 | 0.751 |
  | 4 | 0.222 | 0.704 |
  | 6 | 0.166 | 0.649 |
  | 10 | 0.114 | 0.599 |

  So abc with delta below about 0.08 (b = 4) or 0.25 (b = 6) excludes them. This is exactly what LW Theorem L
  does (lw.md §0 item 2). Szpiro (kappa = 3/2) excludes none of them: 1/4 - 0.704/6 = 0.133 < 0.222.
* **Growth rate.** No improvement conditional on abc from these relations. The barrier profile, the double core,
  satisfies every character consequence at profile level (proved, from Thm 5.1). For monomials outside the
  character lattice, abc gives gaps e^{-w_0 - delta*w_1/2}. These are weaker than Liouville's e^{-w_1/2} unless
  most |v_p| >= 2. A union bound sum_v e^{-w_0(v) - delta*w_1(v)/2} ~ sum_j 4/(delta p_j log p_j) suggests that
  angle models satisfying them exist once P >> M/delta. That is cheaper than the Liouville requirement in
  lw.md Thm 6.1, so the record's fakes at about 8 M log M plausibly survive abc. This last point is
  **heuristic**: a simultaneous abc-gap version of Thm 6.1 is not proved.
* **The constant 3 in W >= 3M log M.** It comes from pigeonhole at primes q <~ M. (abc-Lambda) is insensitive to P,
  so nothing here moves the constant. Whether the constant can rise somewhere between 3 and 8 stays undecided:
  the record's fake constant is 8.
* **Novelty against the record.** The record has abc only for pairs (growth/conditional.md Prop 8.1, the case
  c = e_x - e_y) and its vacuity on chord, triangle, Ptolemy and n-term relations (round8/weakest.md Prop 4.2).
  The relation beta - conj(beta) = 2i Im(beta) for general characters is a product of points, not of chords, so
  Prop 4.2 does not cover it. Its consequence (abc-Lambda), identified with the LW support test, is new as a
  formulation but routine. The record's statement "abc gives nothing" (barrier 6) should be refined: abc via
  character relations equals the LW support test. It kills the round-5/6 Theorem-L families, as LW does, and not
  the double core.

## Novelty and significance

* **Novelty.** VI.1 is new to the record. A grep for parity/binary code, weighted distance, H_2(1/8) and
  logarithmic derivative finds nothing like it. round29/ff.md's binary codes are a different object (unary
  encodings over C(t)). The proof is a short Mason-type argument; the binary-code packing on top is standard.
* **Significance for the goal: low (about 1/10).** VI.1 is a correct, clean theorem for the fixed field F_2. It
  removes the record's twin-based fake profiles from the binary polynomial model and explains why Part V's
  genuine clusters are low-rank. But its gain is the characteristic-2 identity (d/dt kills squares), which is
  provably false over Z (e1). Its honest integer analogue, under abc, is the LW-type support test. That test
  excludes the Theorem-L families but not the skew double core (Thm 5.1). So it changes neither the growth rate
  nor the constant 2/3, conditionally or unconditionally. Astra's own sentence "None of these growth arguments
  improves the integer bound" is accurate.
