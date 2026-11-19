# Referee block "cap2": Part VII.1, constant model (eta = i, eta^2 = -1), p = 3 and general p = 3 mod 4

Referee scripts and outputs: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref7_cap2/`.
All arithmetic is my own. F_{p^2} = Z[i]/(p) polynomials are complex numpy arrays multiplied by FFT, with the rounding error asserted below 0.2 (observed at most 7e-10). F_27 = F_3[a]/(a^3+2a+1) is a different modulus from Astra's a^3-a+1. F_3[t] Euclid is in my own `f3.py`. I did not import Astra's scripts. I read `vii_cap2_check.py` only afterwards, to compare.

## Verdict in brief

**VII.1 in the constant model is CORRECT for every p = 3 mod 4.** I checked every step line by line. I also expanded every polynomial of Astra's own four-mark certificate (degree 619801) with independent code, and all checks pass. Her certificate never expanded the marks. It relied on the proof for the norms and widths, so this is strictly more than she computed. Her "certified upper bound" 307762 is the exact maximal pair degree. Smaller instances:

| p | q | k | m | g | L (min) | D | max pair degree | D/2 |
|---|---|---|---|---|---|---|---|---|
| 3 | 3 | 2 | 2 | 5 | 36 | 144 | 72 | 72 |
| 3 | 27 | 2 | 2 | 365 | 175 | 63700 | 31759 | 31850 |
| 3 | 27 | 2 | 3 | 183 | 378 | 206388 | 103194 | 103194 |
| 3 | 27 | 2 | 4 | 92 | 908 | 578396 | 289172 | 289198 |
| 3 | 27 | 2 | 4 (Astra, T=1+tC^36) | 92 | 973 | 619801 | 307762 | 309900.5 |
| 7 | 7 | 2 | 2 | 25 | 131 | 3144 | 1569 | 1572 |
| 7 | 7 | 2 | 3 | 13 | 392 | 14112 | 7056 | 7056 |

At L = L_min - 1 every one of these fails criticality. The proof's width bound is **attained exactly** by the worst pair in every instance, so the proof's threshold is exact for this family. It is not merely sufficient.

Novelty: new. The record had cap two only at k = 2 with m <= 4, in the VI.4 family (`ref6_ternary.md`, line 198: "k = 2 caps m at 4"). It listed "squarefree or bounded-multiplicity N in odd characteristic" as **open** and out of VI.4's reach (README_VI items 1 and 4; expert note lines 913 and 958). VII.1 closes the bounded-multiplicity half negatively: any cap e_0 >= 2 permits unbounded clusters, at the same rate (1/2) log_2 p as VI.4. The squarefree half (cap 1) stays open.

Significance for Z: none directly. The device uses p-th-power labels with dv = 0, so T - v is separable with T' independent of v, together with F_q-linearised root polynomials. Neither has an integer analogue. The result extends filter F3++: an argument that proves the uniform bound for bounded-exponent integer n must fail over F_p[t] at cap two. Concretely, it has to use squarefreeness (cap one) or characteristic 0. The integer growth rate and the uniform bound are unchanged.

## (a) The root space V of (VII.2): CORRECT

V = { i^p C^p sum_{j<k} x_j t^{pj} : x_j in F_q }, q = p^r, r odd, so q = 3 mod 4 and chi(-1) = -1.

* **Size n = q^k.** The map x -> v is injective because C != 0 and the monomials t^{pj} are distinct.
* **p-th powers.** Each element equals (i C sum y_j t^j)^p with y_j = x_j^{1/p} in F_q, and i^p = -i for p = 3 mod 4. The label set is symmetric, so the sign is harmless. Fixing the coordinate convention v <-> x through i^p (not i) only flips all epsilon_j, which replaces the whole family by its conjugate.
* **Anti-invariance.** sigma is the q-power map on constants, and r odd makes it the nontrivial automorphism of F_{p^2}. It fixes F_q and t and sends i -> -i, so sigma(v) = -v.
* **Descent.** phi is the p^2-power map. It fixes i and C and sends x -> x^{p^2}. Since the forms have F_p coefficients, l(x^{p^2}) = l(x)^{p^2}, and chi is Frobenius-invariant. So every exponent function and every minimum is phi-invariant. phi and sigma generate Gal(F_{q^2}/F_p), because p^r is not in <p^2> when r is odd. Hence F_j and G lie in F_{p^2}[t][X], and P_V(T) lies in F_p[t].
* **Computed.** For q = 27 I grouped the labels into the 249 phi-orbits (sizes 1 and 3), multiplied each orbit in F_9[t][a]/(a^3+2a+1), and asserted that the a and a^2 components vanish. They do, in every instance. I also asserted that the exponents are constant on orbits.

## (b) Collision claim: CORRECT

Take v != w. Then v - w = i^p C^p (sum y_j t^j)^p with 0 != sum y_j t^j in F_q[t] of degree <= k-1. A root alpha of the inner polynomial generates an extension of F_q of degree d <= k-1, so it is a root of t^{q^d} - t, which divides C. So every common zero of T - v and T - w is a zero of C. At a zero of C every label vanishes, because C^p | v, while T = 1. So different T - v have disjoint zero sets.

## (c) T, T', simple zeros, R squarefree in F_p[t]: CORRECT (one clarification)

There are three cases.

* **T = 1 + C^p t^s, p ∤ s.** Then T' = s C^p t^{s-1}. Its zeros are the zeros of C together with t = 0 (when s >= 2). "Since t | C" is exactly what puts t = 0 inside the collision set. That holds because k >= 2, which is forced for m >= 2 nonproportional forms, and t | t^q - t.
* **T = 1 + C^p(t^s + t), p | s.** Then T' = C^p.
* **Astra's certificate, T = 1 + t C^36.** This is a third form, not covered by the displayed rule. It is equally valid, since T' = C^36 and T = 1 mod C.

In all three cases T' vanishes only where T = 1 and every label vanishes. And (T - v)' = T' because v' = 0. So all finite zeros are simple and pairwise disjoint, R = prod_v (T - v) is squarefree, and R is in F_p[t] by (a).

There is also an equivalent clean certificate. Since P_V is additive, R' = P_V'(T) T' = (prod_{v != 0} (-v)) T'. Every root of the right-hand side is a root of C. So gcd(R, R') = 1 if and only if gcd(R, C) = 1, and R = 1 mod C.

**Computed:**
* R' = ±a_0 T' as an exact polynomial identity, where a_0 = prod_{v != 0} v is computed by my own products.
* R mod (t^q - t) = 1, computed by folding exponents.
* In the D = 148 instance, gcd(R, R') = 1 by direct Euclid.

These hold in every instance in the table.

## (d) The marks: CORRECT

* **Norm identity.** sigma F_j = prod_v (X - v)^{1 - chi(eps_j l_j(v))}, so F_j sigma(F_j) = P_V^2. Kernel labels have exponent 1 + 1.
* **Root-polynomial identity.** A = P_{ker l_j} is F_q-linear on V with kernel ker l_j, so A(v) = c eps_j l_j(v). The selected fibres are those where A(v) = c s with s a nonzero square. Since prod_{s square} (Y - c s) = Y^{(q-1)/2} - c^{(q-1)/2}, we get F_j = A(A^{(q-1)/2} - c^{(q-1)/2})^2 = P_V + 2c^{q-1}A - 2c^{(q-1)/2}A^{(q+1)/2}. I checked this line by line. It gives deg_X(F_i - F_j) <= (n+h)/2.
* **Common divisor G.** The expected minimum is 2^{1-m} at all-nonzero labels and 2^{-a} with a < m nonzero, so there is a choice with g >= 2^{1-m} n. With q >= 2^m this gives g - h >= n/2^m. My independent counts are g = 365, 183, 92 for m = 2, 3, 4 (q = 27, forms x, y, x+y, x-y, all signs +). These are maxima over sign choices, and 92 agrees with Astra.
* **Distinctness.** Two exponent rows differ at some label. Because T - v has simple zeros disjoint from all other factors, the corresponding z_j differ even after common division. They are monic, so they are not unit multiples of each other.
* **Degree.** D_0 = (n-g)L needs L > b. **Computed:** for L <= b the degree drops, and the q = 3 run at L = 11 fails the degree assertion.
* **Norm divides R^2.** N = R^2 / (G sigma G)(T). **Computed** as the exact identity N * Q = R^2 with Q = prod (T - v)^{g(v)+g(-v)}. The norms of all marks are equal and real in every instance.
* **Exponents.** A prime over a zero of T - v has exponent 2 - g(v) - g(-v) in {0, 1, 2}. Exponent 2 occurs exactly at "mixed" labels, where some form has character -1 and some form has character +1. Two independent forms give such labels.

## (e) Width: CORRECT and sharp

* The coefficient of X^a in Z_j is an elementary symmetric function of n-g-a labels, so its degree is <= b(n-g-a), with b = p(deg C + k - 1).
* With L > b the maximum over a <= (n+h)/2 - g occurs at the top, which gives 2 deg(z_i - z_j) <= (n+h-2g)L + b(n-h).
* The slack is (g-h)L - b(n-h), so the cluster condition holds **exactly from**
  L_min = max(b+1, p deg C + 2, ceil(b(n-h)/(g-h))).
* **Computed:** the worst pair attains the bound with equality in every instance, and L_min - 1 fails each time (table above). One exception is the p = 7 triple x, x+y, 2x+y with signs (+,+,-), whose pairs lie below the bound and which is already critical at L = 391.

## (f) Padding and rate: CORRECT

* **Padding.** The reserved primes divide C, and R = 1 mod C, so the reserved primes do not divide N_0. Hence e^2 N_0 keeps cap two. The marks e z_j need slack >= u < n. That holds for L >= 2^m(b+3) + p deg C + 2, because (n/2^m)L - bn >= 3n. Binary expansion of u works because the reserved degrees sum to at least n. The padding is a monotonicity device that changes nothing structural.
* **Threshold.** D >= C_p 2^m n^2, with deg C = O(n) and b = O(pn).
* **Rate.** log q <= m log 2 + 2 log p and k = log m / log p + O(1), so the log of the threshold is (2 log 2 / log p) m log m + O(m). Taking m = a log D / loglog D with a < (1/2) log_2 p gives (VII.1). The threshold depends on m, but m grows only through this explicit choice. There is no hidden uniformity.

## (g) Reproduction: done in full

Every instance in the table was built from scratch:

* T was formed.
* All n linear factors T - v were built, with the F_27 descent done by orbits.
* Every z_j was expanded to its full degree.

**Verified exactly** in each instance:

* equal norms in F_3[t]
* distinctness, monic marks, degree (n-g)L
* all pair-difference degrees against D/2
* the identity N * Q = R^2
* R' = ±a_0 T' and R = 1 mod (t^q - t), which together prove R squarefree; for D = 148 also gcd(R, R') = 1 and gcd(N, R) = rad N (degree 222) by Euclid
* deg N against the label census

**Not verified:** a full Euclid gcd for the large instances. The certificate above replaces it, and that certificate is a proof.

**Smallest instances:**
* 2 marks: D = 144 (p = q = 3, L = 36, T = 1 + C^3(t^27 + t)).
* 3 marks at p = 3: D = 206388 (q = 27, L = 378). q = 3 cannot give m >= 3 here, since g <= h for every sign choice.
* 4 marks at p = 3: D = 578396 (L = 908), smaller than Astra's 619801.
* General p: p = 7 gives 3 marks already at D = 14112.

Running times: under 10 s each; Astra's instance took 9.9 s.

## (h) "The remaining double factors come from the character formula itself": CORRECT

R is squarefree and every T - v is coprime to every other. So the exponents of N are exactly the formal ones, 2 - g(v) - g(-v). There are no accidental ones, and the doubled primes are precisely those at mixed labels, where the character exponents 2 and 0 sit on opposite sides.

Share of deg N carried by primes of exponent exactly 2 (counted with multiplicity), and by exponent exactly 1:

| instance | labels with exponent 2 / 1 / 0 | share of exponent 2 | share of exponent 1 |
|---|---|---|---|
| q = 3, m = 2 | 2 / 4 / 3 | 0.500 | 0.500 |
| q = 27, m = 2 | 338 / 52 / 339 | 0.929 | 0.071 |
| q = 27, m = 3 | 520 / 52 / 157 | 0.952 | 0.048 |
| q = 27, m = 4 | 624 / 26 / 79 | 0.980 | 0.020 |

For the q = 3, m = 2 instance this was also checked by Euclid: deg N = 296, deg rad N = 222, and the exponent-2 part has degree 74.

These shares are consistent with, and much stronger than, Astra's lower bound 1 - 2(q+1)/(m(q-1)), which is 0.46 at m = 4. So the constructed parents are cap-two only in a thin sense: about 98% of deg N sits in squares at m = 4. This confirms that the squarefree question is untouched.

## Standing cautions

* **Gcd costs.** The degrees used are those of the actual marks after division by G(T). Padding multiplies by a real e, which is allowed and is counted in D.
* **Points and norm.** The marks are distinct, of equal degree, monic, and have one norm in F_3[t].
* **Conjugates and unit multiples.** Monic marks with distinct zero sets cannot be unit multiples of one another. Conjugates need no separate treatment here.
* **Infinity.** It is inert, since this is the constant extension.
* **Thresholds.** L_min depends on m, but it is never used as if it were uniform.

**Ramified model.** It was not refereed in this block. Note only that the descent there uses the p-power map, and that b contains (1/2) deg Delta.
