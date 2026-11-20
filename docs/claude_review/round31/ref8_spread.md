# Referee report on Part VIII, block "spread": §VIII.4 and the COMPUTED paragraphs

Refereed text: `astra_free_VIII/answer.md`, lines 1905–1938 (§VIII.4, the two COMPUTED paragraphs and the verdict), together with the parts of §VIII.3 and (VIII.10) that the spread proof uses. My scripts and their outputs are in `round31/ref8_spread/`. Every script is my own code; none imports Astra's scripts. Labels: **Proved** (I checked every step), **Computed** (exact arithmetic, reproduced by my code), **Illustration** (floating point, a single random sample, never used as a proof), **Heuristic**, **Gap**, **False**.

## 0. Verdict

No error in the spread block. Every formal statement in §VIII.4 is correct as stated. Every number in the two COMPUTED paragraphs is reproduced exactly by independent code. There are three clarifications. None of them changes a theorem.

1. **The area test (VIII.11) is correct, and its novelty is minimal.** It is round 21's sector count (`round21/prove.md` Lemma 3.1), refined by a radial bin, composed with (VIII.1). To detect the first fake one needs only the sector count applied to the fake's prime angles; (VIII.1) is not needed for that.
2. **The claim that the test has "no verbatim counterpart at ultrametric infinity" is correct only after a clarification.** An ultrametric sector count does exist and was computed in F_9[t]. What has no counterpart is a radial bin narrower than one degree. The detection of the first fake depends on exactly that, since it needs sqrt(P) sigma = q^(0.005+o(1)).
3. **The spread construction is correct as stated, conditional on the §VIII.3/(VIII.10) estimates it inherits.** These are the mod-four support bound, the erasure count and the ordered residue witnesses, which are refereed in other blocks. The spread fakes are *L*-spread in the sense of `round21/prove.md` Definition 5.1. This **discharges the fake-side hypothesis of round-21 Corollary 5.4** (the conditional cap of the window for A_0 plus robust analytic inputs). Astra does not note this consequence. The remaining hypothesis of that corollary is that the *genuine* angles of the chosen primes are spread (round 21 Caveat 3), and it stays open.

The verdict paragraph says that "a genuinely different input is the association of each individual prime norm with its actual Gaussian argument or its prime-level residues". The record already says this, in the same terms: round 21 §6.3 and Theorem C, Lemma T, P3-all*, and round 31 §1 ("forces genuine angles mod pi/4"). See §5. It is a correct restatement, not a new diagnosis.

## 1. The area test (VIII.11)

**Statement (Proved, as stated).** Let a normalised genuine cluster on an arc of length <= C sqrt(R) have all active primes in [P, P e^sigma], with 0 <= sigma, eta <= 1. Suppose the chosen Gaussian factors pi_p have arguments mod pi/2 in an interval of length eta (any choice of pi_p versus its conjugate). Then

  M <= 1 + K_C + c_abs (1 + sqrt(P) sigma)(1 + sqrt(P) eta),

with c_abs absolute and K_C from (VIII.2), depending only on C.

*Check.* The pi_p have distinct norms, so they are pairwise non-associate. Take the representative of each class with argument in a half-open lift of the interval. These representatives lie in the annular sector {sqrt(P) <= |z| <= sqrt(P) e^(sigma/2), arg z in J}. The unit squares around them are disjoint. They lie in the 2^(-1/2)-neighbourhood of the sector, which is contained in the annular sector with radii sqrt(P) e^(sigma/2) + 2^(-1/2) and sqrt(P) - 2^(-1/2), and opening |J| + O(P^(-1/2)). The area of that sector is (1/2)(eta + O(P^(-1/2)))(r_2 - r_1 + O(1))(r_1 + r_2) = O((1 + sqrt(P) sigma)(1 + sqrt(P) eta)). Hence omega(N') = O(...), and (VIII.1) gives the bound on M. Bounded P is absorbed into c_abs. This is round 21's Lemma 3.1 (which counts the full disc of radius Y^(1/2)), with a radial bin added, and the Steiner step replaced by the enclosing annular sector, since the region is not convex. **Novelty: minimal.**

**The first fake of VIII.3 (Proved; exponents recomputed).** Parameters: d = 1.99, H = q^0.08, E = q^1.625, B = q^1.708. On the bounded-digit event, which the construction must realise because the short words use it, the digit series has sup norm at most 2E/(100 H^2 (B-1)) = O(q^(1.625 - 0.16 - 1.708)) = O(q^(-0.243)) in quarter turns. The term c_0 is common to all coordinates. The offsets t_j (at most Delta) and the term B^(-J-2) zeta are smaller. So all angles lie within eta = O(q^(-0.243)) of the direction c_0 + const.

The pigeonhole bin has sigma << q log q / P (2q log P/P suffices). Then sqrt(P) sigma << q^(1 - 0.995) log q = q^0.005 log q, and sqrt(P) eta << q^(0.995 - 0.243) = q^0.752. The product is O(q^0.757 log q), as stated. At least (q+1)/2 active primes lie in this bin with these angles. So the sector count alone (no (VIII.1) needed) is violated by a factor >> q^0.243/log q.

For the (VIII.10) parameters (d = 2 - xi; E/B = q^(-kappa-eps) once the factor 1/(100H^2) is dropped): eta = O(q^(-kappa-eps)) and sqrt(P) sigma = q^(xi/2+o(1)). The product is q^(1-kappa-eps+o(1)) = o(q), as stated.

**Function-field status (Correct with clarification).** Astra's sentence: "a count with no verbatim counterpart at ultrametric infinity". Corrected statement: *the ultrametric sector count exists, but radial bins narrower than one degree do not, and the detection of the first fake needs sqrt(P) sigma << q^(0.243)/log q, i.e. a radial resolution sigma = o(1).*

The computation is in `ultrametric_sector.py`. Take O = F_9[t] (p = 3, infinity inert). The "argument" of monic z of degree n is u(z) = z/zbar in F_9((1/t)), and a sector of depth k is a ball of radius 3^(-k-1). Each sector of depth k <= n contains exactly 9^n/3^k classes. For every k >= n the maximum is 3^n = sqrt(9^n): these are the "real" polynomials, the analogue of the lattice points on an axis. So the ultrametric analogue reads count <= O(max(P eta, sqrt(P))), which is (1 + sqrt(P) sigma)(1 + sqrt(P) eta) with sigma ≍ 1. At sigma ≍ 1 the archimedean bound is q^0.995 · q^0.752 >> q, so it would not detect the first fake either. The archimedean content of the test is therefore the ordering of the prime norms by size inside one dyadic range. "Function-field-free" is correct in this precise sense, and only in this sense.

## 2. The spread construction (§VIII.4, "spread at the same weight")

All items below are **Proved** (each step checked by hand). They are conditional on §VIII.3/(VIII.10) as refereed elsewhere. Quantifiers: the parameters xi, kappa, tau, a, eps, lambda, n_0 and rho are fixed first, and q -> infinity afterwards. Every threshold q_0 depends on all of them, and no uniformity in rho is claimed, as Astra says.

(a) **Background and windings.** Let b_j = (q + 4j + 4q z_j)/(12q). Row sums of K vanish, and (K J_1)_i = sum_j chi(j-i) j ≡ sum_j (j-i)^((q-1)/2) j ≡ 0 (mod q), because the degree (q+1)/2 is less than q - 1. This gives 2(Kb)_i = (2/3)(tau_1 + Kz)_i exactly. Since 2 is a unit mod 3, the differences on S are integers if and only if tau_1 + Kz is constant mod 3 on S. Column 0 is constant (+1) on S, so the condition is the same with or without it.

(b) **The set T and its size.** For fixed z the grid has spacing 1/(3q) and meets each 10P^(-1/2)-neighbourhood of (1/2)Z in at most 60 q P^(-1/2) + 1 points. The z = 0 grid (numerators in [q, 5q)) meets none. So |T| <= 2(60 q^(xi/2) + 1) = O(q^(xi/2) + 1). For j in T, z_j = 0 gives b_j in [1/12, 5/12).

(c) **Full column rank of K_{S^c,T} over F_3.** C is the Paley conference matrix (C^2 = -qI, g = sqrt(-q)). Take a lift v of a nonzero vector mod 3 supported on T. Then Cv ± gv are eigenvectors with eigenvalues ±g. Their difference 2gv is nonzero modulo every prime above 3 (q ≠ 3), so one of them reduces to a nonzero vector. Its support lies in supp(Kv mod 3) ∪ T ∪ {infinity}.

I checked the support lemma itself. SL_2(F_q) acts by signed permutations that commute with C, since C_ij = chi(det(s_i, s_j)), and this moves a 𝔭-unit coordinate to infinity. On the finite part, f-hat is supported on one quadratic class together with 0. Hence f * f(n ·) = (sum f)^2/q, where n is a nonsquare. The constant is a 𝔭-unit, because sum f = -lambda y_infinity. A convolution of two vectors of support b is supported on at most b^2 points. So b >= sqrt(q).

This gives |supp(Kv mod 3)| >= sqrt(q) - |T| - 1. The union bound 3^|T| rho^(sqrt(q) - |T| - 1) -> 0 needs |T| = o(sqrt(q)), which holds since xi/2 < 1/24. Rows outside A are always in S^c, which only helps.

The lemma is a form of the classical square-root bound for quadratic-residue codes. This is from memory (MacWilliams–Sloane, ch. 16); I did not re-read it. Computed: the F_3-eigenspaces of C have minimum weight 6 at q = 11 (the ternary Golay code) and 9 at q = 23, both >= sqrt(q). At q = 127 and 383, every v with |supp v| <= 3 satisfies the bound by a wide margin (`out_support.txt`).

(d) **The choice of u and z.** Put u_S = tau_1 mod 3, solve (K^T u)_T = 0 by (c), and set z = -q^(-1) K^T u. Then Kz = -u + q^(-1)(1^T u) 1, since K K^T = qI - 11^T. So tau_1 + Kz is constant on S, and z_T = 0.

(e) **Separation.** b_j - b_k = ((j-k) + q(z_j - z_k))/(3q), and j ≢ k (mod q). Hence the distance from Z is at least 1/(3q). b_j + b_k = (q + 2(j+k) + 2q(z_j+z_k))/(6q) has an odd numerator, so its distance from Z is at least 1/(6q). The case j = k is included. Directions modulo units are theta mod 1, and conjugation is theta -> -theta. On the bounded-digit event the final angles move by o(1/q). So the separation is (1 - o(1))/(3q) >= c_rho/m, and (1 - o(1))/(6q) up to conjugation.

(f) **Modified angle formula and short words.** The perturbation is divided by 24q. For h <= H balanced, or h <= H_u unbalanced, it contributes at most (0.1 + o(1))/(24q). This uses H E/B = q^(-eps) and H_u/(10 H_u) = 1/10. The margin off the centres is 1/(12q), because v·b has denominator 12q. The required gap 2P^(-h/2) <= 2P^(-1) = o(1/q) for h >= 2. For h = 1, b_j is never in (1/2)Z (12q b_j is odd), and T gives margin 10P^(-1/2), with P^(-1/2) = q^(-1+xi/2) >> 1/q.

A short unbalanced word at a centre (this forces sum v even) gets a gap >> 1/(q H_u) from c_0. The digits contribute o(1/(q H_u)) since H_u^2 E/B = q^(2 lambda - kappa - eps) -> 0. This suffices since d > 1 + lambda.

A balanced word at a centre satisfies 12q v·b ≡ 4 sum v_j j (mod q), so its first power sum vanishes. With n = h/2 >= 2 and Newton's identities (n < q), the first differing moment s lies in [2, n]. Under the reordering 2, 3, ..., n_0, 1, n_0+1, ..., moment s sits at position s - 1 if s <= n_0, and at position s if s > n_0. So the position is at most n - 1 when n < n_0, and at most n otherwise. The leading-digit bound gives >> B^(1-n)/q and >> B^(-n)/q respectively. Both exceed 2P^(-n), since B >> q and (P/B)^(n_0) >= q^2 by the choice (d - b_*) n_0 > 2. Non-primitive words reduce to their primitive part (g_0 < q).

(g) **Winding shift and axiom (R).** For a character v = K_S^T lambda (zero sum), v·theta = lambda·k/2 + lambda·eps/2 exactly. The background contributes lambda·k/2 because 2K_S b = const + k. Every perturbation has constant K_S-image. In (R), psi_u = arg(u)/2 + (pi/4) lambda·k (`round5/construct.md` §1.2, with eps_x = i^(-k_x)). So <v, phi> - psi_u = (pi/4) lambda·eps - arg(u)/2, which is identical to VIII.3. For pairs whose shift is 1/2, the phase sits at the centre 1/2, whose gap radius is the smaller one. It clears that centre by >> Delta/L_0 >> e^(-w/2). The units u ≠ 1 keep the sine bounded below. So characters and (R) are unchanged.

(h) **Tails.** The densities from X_k, c_0 and zeta acquire a factor 24q. For k > n_0 the vanished positions carry exactly the moments {1, ..., k-1}, so the erasure count and the atoms are unchanged. For k <= n_0 the unrestricted count gives a per-token exponent of at most xi/2 - kappa < 0, times q^(O(n_0)). The non-primitive shift costs B^(n_0). Since h > q^kappa, all these fixed polynomial factors are absorbed. Weight, prime height and positive conditional Haar measure are unchanged.

**Consequence not stated by Astra (Proved).** The spread fakes are L-spread in the sense of `round21/prove.md` Definition 5.1, with L = 2 log_2 M + O(1).

- (i) An arc of length eta contains at most (6/pi) q eta (1 + o(1)) + 1 of the angles. This is at most (16/pi) s eta + 4 (Ls)^(1/2), because s >= (q+1)/2.
- (ii) The 2s angles psi_j and pi/2 - psi_j are pairwise >= (pi/12q)(1 - o(1)) apart, which is >= (pi/2) 2^(-L).

The primes lie in [M^(1+eta'), M^(2-delta)] for every fixed eta' < d - 1 and delta < xi. We also have W = O(M log M) and s = O(M). So the window fakes of Corollary 5.4 now exist, with *their* angles spread. The window is capped for A_0 plus (c, gamma)-robust prime-supported inputs (gamma > 1/(1+eta')) once one also has the hypothesis of Caveat 3: the **genuine** angles of the chosen primes are L-spread.

That last hypothesis is **open (Heuristic)**. The construction may choose its primes from any norm interval with sigma sqrt(q log q) = o(1), not only from the pigeonholed bin. So the question becomes an upper-bound sieve for Gaussian primes in polar boxes of angle >= P^(-1/2) in short norm intervals. I did not check it.

## 3. Computed checks: reproduced

All exact. Scripts: `paley_spread.py`, `interp_lift.py`, `integer_example.py`, `erasure_repro.py`.

| claim (Astra) | my result | status |
| --- | --- | --- |
| q = 383: 191 positive rows of column 0 | 191 = (q-1)/2 | Computed |
| 17,766 pairs (her S = first 189 rows), all l^1-length 384 on 382 labels | 17,766 pairs, distances {384}; also all 18,145 pairs of A | Computed |
| 95 moment lifts, congruences, constant row evaluations, max coefficient 7,558 | (i) the proof's own construction (random mean-zero rounding, completion mod q, lattice-Gaussian term added): all 95 OK, with and without column 0; max 9,215 (her S) and 8,661 (random S), below q^(3/2) sqrt(log q) = 18,280; max\|tau\| = 36 (sqrt(q log q) = 47.7). (ii) independent reimplementation of her interpolation completion: max exactly **7,558** | Computed |
| spread: \|T\| = 124, rank 124 mod 3, clearance, separation, windings | \|T\| = 124, rank_F3 K_{S^c,T} = 124, for her S and for a Bernoulli(0.99) S; windings integral (range [-18, 42]); min \|b_j - b_k\| = 1/1149 = 1/(3q) and min \|b_j + b_k\| = 1/2298 = 1/(6q) (both attained); clearance 0.02763 > 10P^(-1/2) = 0.02690; 11,102 length-4 centre words, every one with vanishing first power sum | Computed |
| further instances | q = 1019 (d = 1.99): \|T\| = 124, rank 124, all checks pass. q = 383, d = 1.95: \|T\| = 140, rank 140. q = 47, 127: T is all or most labels (\|T\| = 47 and 122 > \|S^c\|), so the rank fails, as expected outside the asymptotic range | Computed |
| N = 1176852625, ten marks, C < 33/4 | brute-force enumeration of all representations; best ten-point arc C = 8.23932686587 (the associate i·z of her marks); Gaussian gcd 1 | Computed |
| affine rank 5; 6 with indicators mod 3; 9 with mod 7 | 5; 6 (4 of 4 classes used); 9 (8 of 8); mod 11 and mod 19 stay 9 | Computed |
| Smith 2,2,2,2,24; orders 24,8,6,4,24; t* = (7/24, 3/8, 1/6, 1/4, 11/24) | identical; windings k = (0,3,0,1,0,1,-2,1,1) | Computed |
| erasure: q = 127 minima 73 (mod 3), 56 (mod 4), 74 (mod 5), 74 (mod 7); mod 2 at h = 2 gives 2 | identical (my affine-normalised enumeration) | Computed |
| F_11 counts 40,370, 3,520, 440, 0 | identical | Computed |

Remarks on the computed paragraphs.

- **q = 383 is outside the asymptotic regime of the rank argument.** There |T| = 124 > sqrt(q) = 19.6, so the union bound proves nothing and the full rank is an observation. Likewise the support bound (c) is vacuous at that size. The union bound needs roughly q^((1-xi)/2) log(1/rho) >> |T| log 3, i.e. q of order 10^8 at rho = 0.99. This is no error, since no uniform threshold is claimed.
- **The erasure script's thresholds are stricter than the proof's.** It uses isqrt(q-1) + 1 - h ("needs 10" at q = 127), while the proof's bound is sqrt(q) - h - 1 = 8.27. This is harmless: the mod-2 failure (2) is below both.
- **Reference only:** (VIII.1) needs S = {primes = 3 mod 4 up to 283}, K_C = 4168, at C = 8.2394. The example reaches full rank already with {3, 7}. This is consistent: the theorem's condition is sufficient, not necessary.
- **Minor:** `final8.txt` links to `scratchpad/astra/free/…`, but the files are in `round31/astra_free_VIII/`.

## 4. Illustrations (floating point, one sample; not proofs)

These are from `box_cross_test.py` and `box_cross_test2.py`, at q = 383 and P = q^1.99. The 382 consecutive split primes from 138,241 give sigma = 0.062 and sqrt(P) sigma = 23. The tests compare four angle assignments: genuine; concentrated (VIII.3-like); spread (the actual background plus perturbation); Haar.

- **Exact sector counts (constant 1, no O).** "Lattice" below means the number of Gaussian-integer classes in the same box.
  - At the finest scale (boxes of 1–8 primes and 1–3 angles, lattice count 0), spread and Haar both fail by a factor of 6.
  - At this scale the test *is* genuineness: an empty sector cannot contain a genuine pi_p. So every non-arithmetic model fails it, and it is not an input of the record's kind.
  - Restricted to boxes with lattice count >= 5 or >= 20, the maximum ratios are:

    | | lattice >= 5 | lattice >= 20 |
    | --- | --- | --- |
    | genuine | 0.80 | 0.30 |
    | spread | 0.80 | 0.36 |
    | Haar | 0.80 | 0.40 |
    | concentrated | 25.6 | 13.6 |

  - So coarse box counts separate only the concentrated fake.
- **Max-scale cross-gap |sin(phi_p - arg w)| >= c/p, for all w with N w <= p (not associate to the model's own factor).**
  - Genuine angles: 0 of 382 fail at c = 1. This is proved: Im(pi_p wbar) is a nonzero integer.
  - Spread: 226 fail at c = 1, and 25 at c = 0.1.
  - Haar: 198 and 14.
  - Concentrated: 213 and 25.

## 5. What the spread result means (question 5)

**Proved.** Angular concentration is not a necessary feature of window fakes. Spread window fakes pass every coarse sector count (also illustrated in §4). They satisfy round 21's spreading Definition 5.1, which supplies the fake half of Corollary 5.4.

**Natural further tests (Heuristic or open, ranked).**

1. **The max-scale cross-gap / sharp (X+K) union large sieve.** Genuine circles satisfy it with c = 1, and it involves the angles of *all* Gaussian integers, not only the cluster's. Each assigned angle fails it at a fixed c with probability ≍ c. This is a first-moment count: Σ_{N w <= p} 2c/p ≍ c. So constructions of the present type (independent-looking perturbations at scale >> 1/p) satisfy it only on a set of relative measure e^(-Θ(q)) at best. That is a round-21-Theorem-A-type scarcity, not a non-existence proof. The §4 sample fails it heavily.
   - Record status: uncapped and open (`overview/overview.md` §5.5 item 2; `round10/README.md` §3 item 6; `round21/README.md` §1.2 item 7, Theorem LS).
   - It does not force genuineness: composite lattice directions of norm in (p, 4p] satisfy it (round 10).
   - Astra lists "the sharp union sieve" as left open.
2. **P3-all\*: factored prime-level residues.** These are residues r_j(q) with Norm r_j ≡ p_j and the forced Legendre parity, imposed on every monomial (`overview/overview.md` §§5.5, 6.3, 7.4; target (T′)/C1). The spread fake's residues are point-level (s_0 omega_i/conj(omega_i)), with no factorisation. Astra says so explicitly. Still uncapped.
3. **A sparse form of round 21's Theorem C** (rigidity: a genuine background forces phi_p = ± arg pi_p), which is round 21 §6.4 item 1.
4. **Exact fine-scale sector counts** are genuineness itself (§4), and Lemma T then says arithmetic models are actual clusters. They are the conclusion, not an input.

**Is "the association of each prime norm with its actual Gaussian argument" already identified by the record as the missing input? Yes.**

- `round21/fakes.md` §6.3, last paragraph: "an argument proving growth in the window must decide (Q) itself, i.e. use the supercritical gap structure of the cluster's Gaussian primes (a sparse-prime version of Theorem C …)".
- `round21/README.md` §1, Theorem C and Lemma T.
- `overview/overview.md` §5.5: P3-all\*, "prime-level residues with the genuine norm".
- `round31/README.md` §1: the multiplicative sharp sieve "forces genuine angles mod pi/4, so it is the problem itself".
- `expert_note/expert_note.md` §1, lines 68–76: "a new mechanism must be … tied to the actual arithmetic of N".

Astra's verdict agrees with the record. It adds nothing new except the observation that the low-height restriction alone no longer distinguishes fakes from actual clusters.

## 6. Files (`round31/ref8_spread/`)

- `paley_spread.py q dn dd seed mode`: Paley identities, pair distances, the proof's moment-lift construction, T, F_3 rank, u and z, exact windings, clearance, separation, centre words. Outputs: `out_383_astra.txt`, `out_383_random.txt`, `out_383_d195.txt`, `out_1019.txt`.
- `interp_lift.py`: interpolation completion; reproduces the maximum 7,558 (`out_interp.txt`).
- `support_check.py`: F_3 eigencode minimum weights (q = 11, 23); small-support check at q = 127 and 383 (`out_support.txt`).
- `integer_example.py`: the N = 1176852625 example (`out_integer.txt`).
- `erasure_repro.py`: erasure minima, Newton counts, K_C (`out_erasure.txt`).
- `ultrametric_sector.py`: F_9[t] sector counts (`out_ultrametric.txt`).
- `box_cross_test.py`, `box_cross_test2.py`: the illustrations of §4 (`out_box_cross.txt`, `out_box_levels.txt`).
