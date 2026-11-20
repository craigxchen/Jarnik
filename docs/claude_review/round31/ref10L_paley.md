# Referee report, block "paley": Part X.3 (Paley cuts on a fixed modulus)

Refereed: Part X of `astra_lemma_X/answer.md` (section X.3, with the parts of X.1 and X.2 it uses), `final10.txt` items 1, 4, 5, 11.
Scripts (own code, exact integers): `ref10L_paley/paley_checks.py`, `thresholds.py`, `line_lemma.py`, `circuits.py`; outputs `*.out` in the same directory.

## 0. Verdicts

| Claim | Verdict |
|---|---|
| Pair separation (m+1)/2, exactly (m-3)/4 constant columns in every triple | **Proved** and **Computed** (m = 7, 11, 19, 23, 43, all triples) |
| m = 3 mod 8 needed? | **No.** m prime, m = 3 mod 4 suffices (verified at m = 7, 23, which are 7 mod 8). 3 mod 8 was only for the Sidon property (e) of IX.4. |
| X.1 inputs used here (q^2 g^2 \| H; nondegenerate needs max side >= sqrt(qg)) | **Proved** (rechecked by hand) |
| X.2 line lemma (b >= 5) | **Proved** (rechecked) and **Computed** (exhaustive: b = 5, D <= 10; b = 6, D <= 8; b = 7, D <= 6; 0 failures; b = 4 has 45 failures with D <= 6, e.g. the unit 4-cycle) |
| Sign-free Paley table => genuine two-class configuration with the Paley cut | **Proved**, on **any** row subset of >= 5 rows, whenever L < sqrt(q g_min) |
| Full core: no sign-free table at L = C sqrt q | **Proved**, with corrected quantifiers: m >= m_0(C, t); no dependence on A (no upper bound on P is used) |
| Quantitative margin q^{5/8 - eps} | **Proved** (equal number t of primes per column, m >= m_1(eps, t)); for unequal numbers t_v in [1, t] the proof gives q^{1/2 + 1/(8t) - eps} |
| Missing rows: (3/5 + eps) m retained (t = 2), (2/3 + eps) m (t = 1) | **Proved** (inequality recomputed; general t: retained >= ((t+1)/(2t+1) + eps) m) |
| m = 11 erasure check: 67 row sets, 549 annihilators | **Computed** (reproduced; I also verified the sup-norm bound, the support bound and that the images span the zero-sum space on retained rows, which her script does not check) |
| Comparison with IX.4 | Stronger in hypothesis, weaker in length; the genuine-configuration step is IX.4's argument verbatim (section 3) |
| Circle fakes (about rho m/2 rows) | **Not covered**: rho m/2 < m/2 < (3/5) m < (2/3) m. **Heuristic**: no version of this divisor method reaches retained fraction <= 1/2 |

Overall: X.3 is correct. The corrections are of quantifiers and constants only; none affects the conclusions.

## 1. Setting and definitions

* m is a prime with m = 3 (mod 4); chi is the Legendre symbol mod m; K_iv = chi(v - i); S = K - I (so S_ii = -1, S_iv = chi(v - i) otherwise). Since chi(-1) = -1, K^t = -K and K K^t = mI - J, whence S S^t = S^t S = (m+1) I - J, S 1 = S^t 1 = -1, |det S| = (m+1)^{(m-1)/2} (Computed for the five orders; Proved in general).
* Column v carries Q_v, a product of t_v distinct primes, 1 <= t_v <= t; all s = sum t_v primes are distinct and lie in [P, 2P]; q = prod Q_v (squarefree). Any assignment of primes to columns is allowed. At each p | Q_v the two-class cut of the rows is the sign of S_iv.
* For rows i, j: h_ij = product of the primes at whose column S_iv = S_jv. For a prime p, delta_p is a fixed nonzero residue.
* A **sign-free table** on a row set R is a family of integers d_ij (i, j in R, i != j) with |d_ij| <= L, gcd(d_ij, q) = h_ij, and d_ij = +-delta_p (mod p) whenever p separates i and j (sign free per pair and per prime). Since distinct Paley rows are never equal, h_ij != q and d_ij != 0 automatically.
* For a triple, g = product of the primes of the columns constant on the triple; g_min = min over triples of R.

## 2. The reduction to genuine configurations (Proved)

**Lemma A (counting).** Two distinct rows have inner product -1, so they disagree in exactly (m+1)/2 columns. In a triple a column is either constant or splits 2-1, and then separates exactly two of the three pairs. Counting separations: 2 (nonconstant) = 3(m+1)/2, so the number of constant columns is m - 3(m+1)/4 = (m-3)/4. This needs m = 3 (mod 4) for integrality and the Gram identity; it does not need m = 3 (mod 8). Consequently g_min >= P^{t_min (m-3)/4} >= P^{(m-3)/4} on **every** row subset.

Computed (`paley_checks.out`): pair separations {4}, {6}, {10}, {12}, {22} and constant-column counts {1}, {2}, {4}, {5}, {10} for m = 7, 11, 19, 23, 43 over all 35, 165, 969, 1771, 12341 triples.

**Lemma B (X.1, rechecked).** For three rows with edge magnitudes a, b, c: at p not dividing g one edge is 0 mod p and the other two are congruent up to sign, so two of the four Heron factors vanish mod p; at p | g all three edges and all four factors vanish. Hence q^2 g^2 | H(a, b, c). Dividing by g (g divides all three edges since g | h_ij), Q = q/g, T = H(u, v, w)/Q^2 is an integer with T = 0 or 3 mod 4 (Q odd, as all primes are odd once P >= 3), and -l^4 < H(u, v, w) <= 3 l^4 (monotonicity check of 2xy - x^2 - y^2 + 2l^2(x + y) - l^4 on [0, l^2]^2 redone). So max(a, b, c) < sqrt(q g) forces H = 0, i.e. the triangle degenerates. This uses only |d|, the gcd and the squared residues.

**Lemma C (X.2, rechecked).** Five or more points with positive integer distances, all triples degenerate, embed isometrically in Z. Her proof is correct (diametral pair; the only exceptional pairs are complementary with d = l; a fifth point rules them out). Computed exhaustively (`line_lemma.out`).

**Proposition D.** Let R be any set of at least 5 rows and L < sqrt(q g_min(R)). Every sign-free table on R equals |y_i - y_j| for distinct integers y_i of diameter <= L, and y is a genuine two-class configuration with the Paley cut on R: at each p | Q_h, y_i = y_j (mod p) when i, j are in the same class of column h (p | h_ij | d_ij), and y_i != y_j (mod p) across classes (p does not divide d_ij since gcd(d_ij, q) = h_ij). Within a class all residues agree, so there are at most two residues r_0, r_1, distinct when both classes meet R, and r_1 - r_0 = +-delta_p by any cross edge. Precision: the sign of r_1 - r_0 may be -delta_p at some primes, i.e. the configuration realises the cut with the two prescribed classes possibly interchanged at those primes; this is irrelevant for nonexistence.

## 3. Full core (Proved, with corrected quantifiers)

For columns v != w put a = S_v - S_w. Then a^t 1 = 0 and a^t S_h = -1 - (-1) = 0 for h != v, w, so a sums to zero on both classes of column h; hence p | a^t y for every p | Q_h, and q/(Q_v Q_w) | a^t y. The vector a has exactly (m+1)/2 nonzero entries, each +-2, half of each sign (Computed: supports {4, 6, 10, 12, 22}); after translating y into [0, L], |a^t y| <= ((m+1)/2) L. (Her bound 2mL is valid but crude; the round-9 referee note "m + 1 nonzero entries" should read (m+1)/2.) If ((m+1)/2) L Q_v Q_w < q for all v != w, then S^t y = c 1, hence y = -c 1 (S invertible, S^t 1 = -1): all marks coincide, contradicting d_ij != 0.

**Theorem 1 (corrected).** Let C >= 1 and t >= 1. There is m_0(C, t) such that for every prime m >= m_0 with m = 3 (mod 4), every P, every choice of distinct primes in [P, 2P] with 1 <= t_v <= t primes per column and every assignment, there is no sign-free table on all m rows with L = C sqrt q. A sufficient explicit condition is P^{(m-3)/4} > C^2 and C ((m+1)/2) (2P)^{2t} < P^{m/2} (using sqrt q >= P^{m/2}, Q_v Q_w <= (2P)^{2t}). Since tm <= s <= P, this holds for m >= m_0(C, t), uniformly in P. **No dependence on A**: the hypothesis P <= s^A is not used.

Computed (`thresholds.out`, worst case P minimal): the condition holds for all m from 7 when (t, C) = (1, 1), (2, 1), (2, 10), (3, 1000); last failures m = 7 (t = 1, C = 10), 11 (t = 1, C = 1000), 23 (t = 1, C = 10^6), 11 (t = 2, C = 10^6), 7 (t = 3, C = 10^6), checked for m < 3000.

**Comparison with IX.4.** IX.4 (with the round-9 strengthening) excluded **genuine** configurations with the Paley cut up to length q^{1 - O(1/m)}, for every residue choice and every prime assignment (its proof never used the bins). Theorem 1 excludes the larger class of **sign-free tables** (magnitudes, gcds and squared residues only; signs free per pair and per prime) but only at lengths below sqrt(q g_min) = q^{5/8 + o(1)}. Its proof is Proposition D followed by IX.4's column argument verbatim. The genuinely new content is Proposition D (sign-free => genuine on Paley cuts), which also holds on every row subset. In the range q^{5/8 + o(1)} <= L < q^{1 - O(1/m)} sign-free Paley tables are not excluded by anything in Part X (Open).

## 4. The quantitative margin (Proved)

With exactly t primes per column: log g_min >= t ((m-3)/4) log P and log q <= tm log(2P), so log g_min / log q >= (1/4)(1 - 3/m) log P / log(2P) -> 1/4 (P >= tm/2 -> infinity), and an upper bound 1/4 + o(1) is attained by the triple whose constant columns carry the largest primes. Hence for every eps > 0 there is m_1(eps, t) such that for m >= m_1, L = q^{5/8 - eps} < sqrt(q g_min) and ((m+1)/2) L (2P)^{2t} < q, so no sign-free table exists. The exponent is right: the binding step is the triangle step (sqrt(q g) = q^{5/8 + o(1)}), not the column step (which tolerates q^{1 - 2/m + o(1)}). With unequal t_v in [1, t] the same proof gives log g_min / log q >= (1 + o(1))/(4t) and exclusion at L = q^{1/2 + 1/(8t) - eps}. Computed (minimal P): min log g / log q = 0.221 (m = 43), 0.241 (m = 503), 0.245 (m = 19963) for t = 1, converging slowly from below; at finite m the excluded length is q^{(1 + ratio)/2}.

## 5. Missing rows (Proved; m = 11 reproduced)

Let M be the k deleted rows, R the rest, T the (k+1) x m matrix [1^t; S_M], r = rank T <= k+1. For x in ker T, a = Sx vanishes on M and has 1^t a = -1^t x = 0; a^t S_v = x^t S^t S_v = (m+1) x_v - 1^t x = (m+1) x_v, so a^t S_v = 0 for v outside supp x and q / prod_{v in supp x} Q_v divides a^t y (y from Proposition D, extended by anything on M, since a vanishes there). Cramer vectors with pivot set J (r columns, nonsingular minor) have support <= r + 1 and sup norm <= r^{r/2} (Hadamard, entries +-1). They span ker T (free coordinates independent). |a_i| <= (r+1) r^{r/2}, so |a^t y| <= (m (r+1) r^{r/2}/2) L; her bound without the 1/2 is valid. With exactly t primes per column the sufficient inequality is

t (r+1) log(2P) + (r/2) log r + log(C m (r+1)) < (1/2) log q,  with (1/2) log q >= (tm/2) log P.

Put theta = limsup log m / log P; theta <= 1 always (tm <= P). The left side is (t + theta/2) k log P (1 + o(1)) + O(log C), so the inequality holds for k <= (t/(2t + theta) - eps) m, m >= m_2(eps, t, C). Worst case theta = 1 (P of order m log m, as in IX.4): k <= (t/(2t+1) - eps) m, i.e. retained rows >= ((t+1)/(2t+1) + eps) m: **3/5 for t = 2, 2/3 for t = 1**, as she states. Finally S is invertible and maps ker T onto {a : a_M = 0, 1^t a = 0}; a^t y = 0 for a spanning set forces y constant on R, a contradiction for |R| >= 2.

Computed (`paley_checks.out`): m = 11, k <= 2: 67 row sets, 549 annihilators (10 + 11*9 + 55*8); m = 19, k <= 2: 191 sets, 3077 annihilators; m = 11, k <= 4: 562 sets, 3684 annihilators. In every case: T x = 0, support <= r + 1, sup <= r^{r/2} (observed maxima 1, 2, 4, 16, 48 for r = 1..5), a_M = 0, sum a = 0, a^t S_v = (m+1) x_v, and the images span a space of dimension m - 1 - k. Computed thresholds (`thresholds.out`, C = 1, minimal P): t = 2, m = 19963: k_max/m = 0.423, matching 2/(4 + theta) with theta = log m/log P = 0.75; t = 1: 0.367.

**Remark (Computed / Heuristic).** The Cramer vectors are not optimal: for small k, two columns of T may coincide, giving x = e_v - e_w of support 2 (Computed: every row set at m = 11, k <= 5 and m = 19, k <= 3 has a kernel vector of support <= r; `circuits.out`). A Siegel-lemma choice of windows would replace (r/2) log r by o(r log P) if one could control the **largest** successive minimum of the window lattices, giving retained >= (1/2 + eps) m for every t; I could not prove this (Gap). Any version of the method needs |supp x| < m/2 (the divisor q/prod Q_v must exceed sqrt q), and generically a kernel vector of T has support >= r + 1 = k + 2, so the method should not reach a retained fraction <= 1/2 (Heuristic).

## 6. What remains

* **The circle fakes are not covered.** Their marks are a subset of {i : K_i0 = 1} of relative density rho < 1, i.e. about rho m/2 of the m rows. Since rho m/2 < m/2 < (3/5) m (t = 2) and < (2/3) m (t = 1), they lie **outside** the excluded range, and by the heuristic above outside the reach of any refinement of this column argument. On such row sets Proposition D still applies: any sign-free table is genuine (Proved), so the question there is exactly the genuine-configuration question for restricted Paley cuts (Open).
* Sign-free Paley tables at lengths q^{5/8 + o(1)} <= L < q^{1 - O(1/m)} (Open).

## 7. Corrected theorem (all quantifiers)

Let m be prime, m = 3 (mod 4), S = K - I the Paley core; let q = prod_v Q_v with Q_v a product of t_v in [1, t] distinct primes, all distinct and in [P, 2P] (P arbitrary; any assignment; any nonzero delta_p).

1. (Reduction, any R, |R| >= 5.) If L < sqrt(q g_min(R)) (in particular if L <= C sqrt q with C^2 < P^{(m-3)/4}), every sign-free table on R is |y_i - y_j| for a genuine two-class configuration y on R with the Paley cut (classes possibly interchanged at some primes), of diameter <= L.
2. (Full core.) For C >= 1 there is m_0(C, t) such that for m >= m_0 no sign-free table on all m rows exists with L = C sqrt q. Sufficient: P^{(m-3)/4} > C^2 and C (m+1)(2P)^{2t} < 2 P^{m/2}.
3. (Margin.) If t_v = t for all v, for every eps > 0 there is m_1(eps, t) such that for m >= m_1 none exists with L = q^{5/8 - eps}; for t_v in [1, t], none with L = q^{1/2 + 1/(8t) - eps}.
4. (Missing rows.) If t_v = t for all v, for every eps > 0 and C >= 1 there is m_2(eps, t, C) such that for m >= m_2 no sign-free table with L = C sqrt q exists on any row set R with |R| >= ((t+1)/(2t+1) + eps) m (3/5 for t = 2, 2/3 for t = 1); if log m <= theta log P, the bound improves to |R| >= (1 - t/(2t + theta) + eps) m.

None of m_0, m_1, m_2 depends on A or on P (beyond the existence of the primes).
