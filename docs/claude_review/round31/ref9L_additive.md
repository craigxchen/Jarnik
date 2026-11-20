# Referee report, round 31, Part IX, block ADDITIVE (IX.1, IX.3, and two paragraphs of IX.6)

Source: Astra, `astra_lemma_IX/answer.md`, Part IX, lines 1998-2017 (IX.1), 2029-2037 (IX.3), 2122-2130 (IX.6,
"a criterion on the modulus alone" and "ordinary density of moduli"). Record consulted: expert note §§1, 2, 5;
`round31/README_VIII.md`, `README_VIIIL.md`, `README_II.md`, `README_VI.md`, `ref2_tri.md` §1.1,
`ref7_hadamard.md`, `ref8L_attack.md` §§0-2; `round8/stability.md` §3 (Cor 3.2, Thm 3.3); `round21/prove.md`.
I did not read Astra's script. My scripts (own code, exact integer arithmetic for all arithmetic; floating point
only for angles of exact Gaussian integers) are in
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref9L_additive/`:
`twin.py`, `twin2.py` (twin), `circle.py`, `para.py` (circles), `sharp.py` (sharpness constructions); outputs
`run_*.txt`. No web search was made. Out of block and not checked: IX.0, IX.2, IX.4, IX.5 and the other paragraphs
of IX.6.

## 0. Verdicts

| claim | verdict | new to the record? | bears on worst-case growth? |
|---|---|---|---|
| IX.1 no affine 3-flat in the twin signature set, `q > (4C)^6`, `q > C^2` | **Proved** (checked line by line) | new as a theorem; exponent 2/3 is sharp, matched by the record's cube family | no |
| IX.1 sumset estimate `B(B-1) <= (D-1)(1+sqrt(1+8(D-1)))`, `B << D^{3/4}` | **Proved**; superseded by IX.3 | new, low | no |
| IX.3 twin: at most `(1+C^2)(1+2C^2)` unordered pairs per signature difference; `B <<_C sqrt D` | **Proved** | new: twin analogue of the record's distinct-label lemma | no |
| IX.3 twin corollaries (`D = o(s^2)` gives `B = o(s)`; affine dimension `r <= (2-eps)log_2 s` gives `B << s^{1-eps/2}`) | **Proved** (trivial from the above) | new, trivial | no (a priori `D` bound needed) |
| IX.3 twin alternative ("no parallelogram below `sqrt q/2`") | **Proved**; sharp up to factor 4 | new, low | no |
| IX.3 circle: `M <<_C sqrt D` for unary threshold signatures | **Proved**, with the constant improved (`c < 2` suffices; Astra's `c < 1/10` is safe) | **not new in substance**: a direct corollary of the refereed distinct-label lemma (`README_II.md`, `ref2_tri.md` §1.1) | no; fakes are Sidon |
| IX.6 twin criterion `B <<_C 1 + W/log P_*` | **Proved** | new as a twin statement; elementary transcription | no |
| IX.6 circle criterion `M <<_C 1 + W/max e_p log p` | **Proved** | **not new**: `round8/stability.md` Cor 3.2 and Thm 3.3(2), with better constants | no |
| IX.6 density: `B_max(q)/omega(q) -> 0` in natural density | **Proved** (weak form of a stronger fact) | new, trivial | no (typical `q`, not smooth `q`) |
| Computations (twin small `q`; circles `N <= 10^6` and structured) | **Computed**: 0 violations of every proved inequality | | |

Overall: every statement in this block is correct. None of it touches the growth rate for worst-case parents. The
bound `B << sqrt D` (twin) and `M << sqrt D` (circle) is a structural statement, namely that the signature set is a
bounded-multiplicity Sidon set. It turns into a cardinality bound only when `D` is bounded a priori. The circle
fakes and the Part VIII Legendre fakes P2 have Sidon signatures (`D = 1 + binom(M,2)`), so the bound is saturated
there and gives nothing better than `M << s`.

## 1. IX.1: no affine three-flat of signatures (Proved)

### 1.1 Setting
`q >= 2`. At each exact prime power `m = p^f || q` two residues `r_0(m), r_1(m)` mod `m` are prescribed (a single
residue counts as `r_0 = r_1`, and then the bit is constant). A mark is an integer `y` in an interval of length
`L = C sqrt q` with `y = r_0(m)` or `y = r_1(m)` mod `m` for every `m`. Its signature is `sigma(y) in F_2^{omega(q)}`,
whose bit at `m` is 0 if `y = r_0(m)` and 1 otherwise. `A` is the set of signatures and `D = |A + A|` (this includes 0).

### 1.2 Statement (corrected only cosmetically)
**Theorem (IX.1).** Let `C > 0` and `q > max(C^2, (4C)^6)`. Then `A` contains no affine 3-dimensional subspace of
`F_2^{omega(q)}`. More precisely, if eight marks have signatures forming an affine 3-space, then
`q <= (4L)^3`, so that `L >= q^{2/3}/4`.

### 1.3 Check of the proof, step by step
1. *Injectivity.* Equal signatures give `y = y'` mod every `m`, so mod `q`; since `|y - y'| <= L < q`, `y = y'`.
   This uses `q > C^2`. Correct.
2. *Local bits.* If `sigma(y_x) = sigma_0 + Mx` for `x in F_2^3` (an affine 3-space, `M` injective), the bit at `m`
   is `sigma_0(m) + u_m . x` with `u_m` the row of `M` at `m`. It is constant or affine linear, and the nonzero
   linear part `u_m` is unique. Assign `m` to `u_m`. Let `Q_u` be the product of the prime powers assigned to `u`.
   Distinct `u` receive disjoint sets of coprime prime powers, so `prod_u Q_u | q`. Correct.
3. *Divisibility.* `F_u = sum_x (-1)^{u.x} y_x`. Take `m` with `u_m = v != u`.
   * If `v = 0`, all `y_x` are congruent mod `m`. Since `sum_x (-1)^{u.x} = 0` for `u != 0`, `m | F_u`.
   * If `v != 0, u`, then `y_x = r_{c + v.x}(m)` mod `m`. Each fibre `{v.x = e}` is a coset of `v^perp`, of size 4,
     and `sum_{x in coset of v^perp} (-1)^{u.x} = 0` exactly when `u` is not in `(v^perp)^perp = {0, v}`. Here
     `u` is neither, so `m | F_u`.
   * *Classes agreeing modulo a lower power.* The argument uses only the congruence `y_x = r_bit mod p^f` for the
     full exact power. Whether `r_0 = r_1` mod `p^g` for some `g < f` plays no role. Correct as stated.
   Hence `q/Q_u | F_u`.
4. *Size.* The four `+` terms and the four `-` terms can be paired, and each paired difference has absolute value
   `<= L`. So `|F_u| <= 4L`, and `F_u != 0` forces `q/Q_u <= 4L`.
5. *Three independent frequencies.* Walsh inversion over `Q` gives `8 y_x = sum_u F_u (-1)^{u.x}` (with
   `F_0 = sum y_x`). If the `u != 0` with `F_u != 0` span a subspace of dimension `<= 2`, choose `h != 0`
   orthogonal to all of them. Then `y_{x+h} = y_x`, contradicting distinctness. So there are three linearly
   independent, hence distinct, `u_1, u_2, u_3` with `F_{u_i} != 0`, and `Q_{u_1} Q_{u_2} Q_{u_3} | q`.
6. *Conclusion.* `q >= prod Q_{u_i} >= (q/(4L))^3`, so `L >= q^{2/3}/4`. With `L <= C sqrt q` this needs
   `q^{1/6} <= 4C`, which is false for `q > (4C)^6`. Correct.

The proof is characteristic-zero in an essential way: step 5 uses inversion over `Q`. In characteristic 2 the
signs `(-1)^{u.x}` are all 1, and the record's Sylvester–Hadamard families over `F_{2^m}[T]` (`ref8L_attack.md`
§0 item 6) are unbounded at exponent exactly 2. They contain many affine 3-flats. IX.1 is therefore the
integer-side reason why those families have no counterpart over `Z`. This agrees with Astra's own remark in IX.6.

### 1.4 Sharpness (Computed)
Take pairwise coprime `Q_1, Q_2, Q_3 ~ q^{1/3}`, `q = Q_1 Q_2 Q_3`, and `y_x = x_1 Q_2Q_3 + x_2 Q_1Q_3 + x_3 Q_1Q_2`.
Modulo each prime power of `Q_i` the residue depends only on `x_i`, so it takes two classes. The signatures form an
affine 3-space, and the span is `Q_2Q_3 + Q_1Q_3 + Q_1Q_2 ~ 3 q^{2/3}`. `sharp.py` checks
`(Q_i) = (101,103,107), (2^7,3^4,5^3), (1009,1013,7133)`: span/`q^{2/3}` = 3.00, 3.07, 4.11. The three
coordinate `F_u` are nonzero and divisible by `q/Q_u`. So the exponent 2/3 in (IX.1) is sharp, and the constant
lies between 1/4 and 3. This is the `n = 3` cube family of `ref8L_attack.md` §2.2 (exponent `n/(n-1) = 3/2`).

### 1.5 Sumset estimate (Proved; superseded)
Fix `0 != h in A + A`. The set `A cap (A+h)` is a union of `m_h` pairs `{a, a+h}`, and `m_h` is the number of
unordered pairs of `A` with difference `h`. Let `E_h` be its image in `V/<h>`, so `|E_h| = m_h`. Two different
unordered pairs from `E_h` with equal sums are disjoint. In characteristic 2, four distinct points with zero sum
form an affine plane: the three differences from one of them are nonzero and pairwise distinct, so any two are
independent. The preimage of that plane is an affine 3-space inside `A`, contradicting IX.1. So the
`binom(m_h, 2)` pair sums are distinct. They are nonzero images of elements of `A + A`, and both `0` and `h` map to
`0`, so `binom(m_h,2) <= D - 2 <= D - 1`. Summing `m_h` over the `D - 1` nonzero `h`, and using
`sum_h m_h = binom(B,2)`, gives `B(B-1) <= (D-1)(1 + sqrt(1 + 8(D-1)))`, hence `B << D^{3/4}`. Correct.

The sentence "Such bounds transfer to the log-height rate since `s! <= q`" has content only together with an a
priori bound on `D`. With only `D <= 2^s` it gives nothing.

## 2. IX.3, twin: bounded multiplicity of signature differences (Proved)

### 2.1 Statement
**Theorem.** In the twin with `L = C sqrt q < q`, for every `h != 0` at most `(1 + 2C^2)(1 + C^2)` unordered pairs
of marks have `sigma(y) + sigma(z) = h`. More precisely, with `H = H(h)` the product of the exact prime powers at
the coordinates where `h = 0` and `K = q/H`, there are at most `(1 + floor(2L/K))(1 + floor((L/H)^2/K))` such
pairs. Consequently `B(B-1)/2 <= c_C (D - 1)` with `c_C = (1+C^2)(1+2C^2)`, and `B <= 1 + sqrt(2 c_C (D-1))`.

### 2.2 Check
* At the coordinates where `h = 0` the two marks lie in the same class, so `H | y - z`. Since `y != z`, `H <= L`.
* At each `m | K` the bit differs, so `{y, z} = {r_0, r_1}` mod `m` with `r_0 != r_1`. Hence `y + z = r_0 + r_1`
  and `(y - z)^2 = (r_0 - r_1)^2` mod `m`. Since `gcd(H, K) = 1`, `t = (y-z)/H` has `t^2` in one class mod `K`.
* `y + z` lies in an interval of length `2L`. Its class mod `K` is fixed and `K = q/H >= q/L`, so there are at most
  `1 + 2L/K <= 1 + 2L^2/q = 1 + 2C^2` values. `t^2` is an integer in `[0, (L/H)^2]` in a fixed class mod `K`, so there
  are at most `1 + L^2/(H^2 K) = 1 + C^2/H <= 1 + C^2` values. Both ranges are right.
* `(y + z, H|t|) = (y + z, |y - z|)` determines `{y, z}`. Summing over `h` gives the rate. Correct.

In squarefree language, `h` is the same as the label `H_{yz} = gcd(y - z, q)`. The theorem says that each value of
`gcd(y - z, q)` occurs for at most `c_C` pairs. This is the integer twin of the record's distinct-label lemma on
circles (§3.4 below).

### 2.3 What it says and when it is nontrivial
Always `D <= binom(B,2) + 1`. The bound `B <<_C sqrt D` is the reverse inequality `D >>_C B^2`, so it holds with
content exactly as the structural statement **"`A` is a `c_C`-fold Sidon set"**. It becomes a bound on `B` in terms
of `s = omega(q)` only when `D` is known a priori to be small:
* `D = o(s^2)` gives `B = o(s)`;
* if `A` lies in an affine subspace of dimension `r`, then `A + A` lies in a linear subspace of dimension `r`, so
  `D <= 2^r`. Hence `r <= (2 - eps) log_2 s` gives `B <<_C s^{1 - eps/2}`, and `r <= log_2 s + O(1)` gives
  `B <<_C sqrt s`.

**Class killed.** It kills signature sets of small doubling, `D = o(B^2)`. Examples are large subsets of
low-dimensional affine subspaces, that is, Hadamard/Walsh profiles in which each local bit is an affine linear
function of a label `x in F_2^r` (at most `O_C(2^{r/2})` marks out of `2^r` labels).

**Sylvester–Hadamard profile, checked.** Let the bits be `u.x` for `u in F_2^r \ {0}` (each column may be carried by
a block of primes) and use all `n = 2^r` labels. Then `A` is a linear space of dimension `r`, `D = n = B`, and the
theorem gives `n <= 1 + sqrt(2 c_C (n-1))`, so `n <= 2 c_C + 1`. IX.1 gives more, uniformly in `C`: for
`q > (4C)^6` the profile is impossible as soon as `r >= 3`. With only a subset `X` of the labels, `B <<_C 2^{r/2}`.
These are exactly the profiles that realise unbounded clusters in characteristic 2 (`README_VI.md`: "Walsh/Sylvester
... the shape of all genuine `F_2` clusters"; `ref8L_attack.md` §0.6). So over `Z` they are excluded, with no
characteristic-free argument possible (§1.3).

**What it does not kill.** Sidon signature sets. Astra's IX.4 relaxation is stated to have `D = 1 + binom(m,2)` (I did
not check IX.4), and nothing forces genuine smooth-parent configurations to have small doubling. So IX.3 does not
bear on the worst-case rate.

### 2.4 The alternative mechanism (Proved)
On a 2-flat (`x in F_2^2`), `|F_u| <= 2L`. At least two independent `u` have `F_u != 0`, so `q <= (2L)^2`: there
is no signature parallelogram when `L < sqrt q/2`. Sharpness (Computed, `sharp.py`): `{0, Q_2, Q_1, Q_1 + Q_2}` with
`q = Q_1Q_2` is a parallelogram of span `Q_1 + Q_2 ~ 2 sqrt q`. The constant is therefore sharp up to a factor 4.

## 3. IX.3, circle: `M <<_C sqrt D` (Proved; constant recomputed; not new in substance)

### 3.1 Setting and statement
Let `z_1, ..., z_M` be distinct lattice points of `|z|^2 = N` on an arc of length `C sqrt R` (`R = sqrt N`), so of
angular diameter `<= C N^{-1/4}`. Divide by the Gaussian gcd. Each point becomes
`w_i = eps_i prod_{p | N'} pi_p^{a_ip} conj(pi_p)^{e_p - a_ip}` with norm `N' <= N`, for a fixed choice of `pi_p`.
The unary signature of `i` has the bits `1{a_ip >= j}`, `1 <= j <= e_p`. `A` is the set of these signatures in
`F_2^{Omega(N')}` and `D = |A + A|`.

**Theorem.** `M <= c'_C sqrt D`. More precisely: in one unit class (same `eps_i`), on any subarc of angular diameter
`< 2 N''^{-1/4}` (with `N''` the normalised norm of the points of that subarc), distinct unordered pairs have
distinct signature differences. Since `C N^{-1/4} <= C N'^{-1/4}`, at most `4 ceil(C/c)` such pieces cover the
cluster. No gcd cost is discarded: normalisation only enlarges the angular room relative to `N'^{-1/4}`.

### 3.2 The gap, derived
For a nonzero `v in Z^{omega(N')}` put `beta = prod_{v_p > 0} pi_p^{v_p} prod_{v_p < 0} conj(pi_p)^{-v_p}`. Then
`alpha(v) = prod_p (pi_p/conj pi_p)^{v_p} = beta/conj(beta)` and `|beta|^2 = e^{w(v)}`, where
`w(v) = sum_p |v_p| log p`. Since `beta` contains at most one of `pi_p` and `conj pi_p` at each `p`, `beta` and
`conj beta` are coprime. Hence `beta != +- conj(beta)`, so `Im beta` and `Re beta` are nonzero integers, and

`|alpha(v) - 1| = 2|Im beta|/|beta| >= 2 e^{-w(v)/2}`, `|alpha(v) + 1| = 2|Re beta|/|beta| >= 2 e^{-w(v)/2}`.

Astra's `e^{-w/2}` is correct but loses a factor 2. (Computed: 0 violations of `2 sin(angle/2) >= 2e^{-w/2}` over
all same-unit pairs in every window of §5.)

### 3.3 The two cases
*Case `v != +- v'`.* Equal unary differences mean equal intervals `(min, max]` at every prime, so `|v_p| = |v'_p|`.
Hence `b_+- = (v +- v')/2` are integral, and at each `p` exactly one of them is nonzero (when `v_p != 0`), with
absolute value `|v_p|`. They are nonzero because `v != -+ v'`, and `w(b_+) + w(b_-) = w(v) <= log N''`.
`alpha(b_+-)^2 = alpha(v) alpha(v')^{+-1}` has angle `<= 2 delta`, where `delta` is the subarc diameter, so
`alpha(b_+-)` is within angle `delta` of `+1` or `-1`. The cheaper one has `w <= (1/2) log N''`, and the gap gives
distance `>= 2 N''^{-1/4}`. Contradiction once `delta < 2 N''^{-1/4}`. Astra's `c < 1/10` is more than enough.
(The cross combination actually has angle `<= delta`, which would allow `c < 4`. I prove only `c < 2`.)

*Case `v = +- v'`, made explicit.* First, the rows are affinely independent. A pair in the subarc has
`2 e^{-w/2} <= delta < 2 N''^{-1/4}`, so `w(a_i - a_j) > (1/2) log N''`. An affine dependence
`sum lambda_i a_i = 0`, `sum lambda_i = 0`, `lambda != 0` gives two probability vectors `mu, nu` on disjoint sets of
rows with equal barycentre. For independent `i ~ mu`, `j ~ nu`, normalise each coordinate to `X = a_ip/e''_p` and
`Y` in `[0,1]` with equal means `m`. Then `E|X - Y| <= E(X + Y - 2XY) = 2m(1-m) <= 1/2`, because
`|x - y| <= x + y - 2xy` on `[0,1]^2` (for `x >= y` the difference is `2y(1-x) >= 0`). So `E w(a_i - a_j) <= (1/2) log N''`,
and some cross pair violates the strict lower bound. Now suppose `a_i - a_j = a_k - a_l` (after possibly swapping
`k, l`). Then `a_i - a_j - a_k + a_l = 0`. If `{i, l}` and `{j, k}` are disjoint as index sets this is a nontrivial
affine dependence with coefficients summing to 0. Otherwise `i = k` (since `i != j` and `l != k`) or `l = j`, and the
relation forces the other equality, because distinct points in one unit class have distinct rows. So `(k, l) = (i, j)`.
This is the "algebra of signed coordinate vectors". Correct.

Distinct points in one unit class have distinct rows, hence distinct signatures. So the retained `m >= M/(4 ceil(C/c))`
points have `binom(m,2) <= D - 1`. Correct.

### 3.4 Status relative to the record, and the fakes
* The argument uses only the monomial (Liouville) gap and the archimedean angle bound, so it holds verbatim in every
  angle model that has the monomial gap. That is Astra's own remark in IX.5, and it is correct.
* **It is not new in substance.** The record already has a refereed stronger statement: near-Liouville edges have
  distinct labels `prod_p p^{|v_p|}` in the angle model, including prime powers (`README_II.md` lines 187-188,
  `ref2_tri.md` §1.1, Astra free II.5). Equal unary differences imply equal labels, so IX.3-circle is an immediate
  corollary. The `b_+-` device is literally the same.
* *Reconciliation with the fakes.* The fakes satisfy the monomial gap, so they satisfy `M << sqrt D`. They do so by
  being Sidon: P2 has all pair sums of unary signatures distinct (Astra IX.6, consistent with the distinct-label
  lemma), so `D >= 1 + binom(M,2) ~ s^2/8` at `M ~ s/2`. There is no contradiction, and the inequality cannot improve
  the rate unless some genuinely global input forces `D = o(s^2)`. Local inputs cannot force it, by the fakes.
* The Sylvester–Hadamard profile on circles: all `2^r` rows give pairs `{x, x+h}` with equal signature difference, so
  `2^r <= c'_C^2`. This is the integer counterpart of the characteristic-2 Walsh clusters. The char-2 failure has a
  clean explanation: over `F_2[t]`, `|alpha^2 - 1| = |alpha - 1|^2`, so `alpha(b)^2` near 1 only gives `alpha(b)` within
  the square root of the distance. In characteristic 0, `|alpha^2 - 1| = |alpha - 1||alpha + 1|` and one factor is
  about 2.

## 4. IX.6: the modulus-alone criterion and the density statement

### 4.1 Twin criterion (Proved)
**Statement.** For `q > 1`, `W = log q`, `P_* = max_{p^f || q} p^f`, and `C >= 1` (enlarge `C` if needed): every
admissible configuration has `B <= 2b`, where either `b <= 1`, or `b (log P_* - 2 log C) <= W - 2 log C`. Hence
`B <<_C 1 + W/log P_*`; when `P_* <= C^4`, `q <= lcm(1..C^4)` and `B <= C sqrt q + 1` are bounded.

**Check.** At `P_*` one class holds `b >= B/2` marks. For a pair among them, the product of the agreeing prime powers
divides `y_i - y_j != 0`, so it is `<= C sqrt q`, and the separating log-weight is `D_ij >= W/2 - log C`. Each other
prime power separates at most `b^2/4` pairs, and `P_*` separates none. So
`binom(b,2)(W/2 - log C) <= (b^2/4)(W - log P_*)`, which rearranges to the displayed inequality. Correct. The
corollaries are correct: `log P_*/log log q -> infinity` gives `B = o(log q/log log q)`, and
`s log P_*/W -> infinity` gives `B = o(s)`.

**Novelty and significance.** This is the twin transcription of the record's circle results. It is the usual
collision averaging (Cilleruelo–Córdoba type, larger sieve). On the smooth parents that matter (comparable primes,
`log P_* ~ log s`) it gives `B << s`, exactly the larger-sieve order. No gain.

### 4.2 Circle criterion (Proved; in the record)
`M <<_C 1 + W/max_p e_p log p` with `W = log N'`. Astra's proof is correct: one unit class, a subarc with `c < 1`
(gap `w > W/2`; `c < 2` suffices with §3.2), split `[0, e_p]` at `floor(e_p/2)`, unary cut capacities `b^2/4`, which
gives `b < 2W/(e_p log p)`. **It is in the record**, with better constants: `round8/stability.md` Thm 3.3(2)
(`e log p < 4W'/M + 8 log^+ C` for every `p^e || N'`) and Cor 3.2 (`M < 2 log N/(log P_+(N) - 4 log^+ C)`). It is
also dominated by the Theorem 1.2 consequence `P_+(N') >= (2 - o_C(1)) M log M` (expert note §1). Not new.

### 4.3 Density (Proved, but weak and irrelevant to the hard case)
**Statement.** For each `C`, `B_max(q)/omega(q) -> 0` in natural density, where `B_max(q)` is the maximum over
intervals of length `C sqrt q` and over residue choices.

**Check.** For `q <= X` and `y = X^{1/k}`, the average of the log of the `y`-smooth part is
`(1/X) sum_{p^a <= X, p <= y} log p floor(X/p^a) <= sum_{p <= y} log p/(p-1) = (log X)/k + O(1)` (Mertens). By
Markov, the smooth part exceeds `sqrt X` on a fraction at most `2/k + o(1)`. Off that set, and for `q > sqrt X`, `q`
has a prime factor `> y >= q^{1/k}`, so `log P_* >= W/k` and `B_max <<_C k`. `omega(q) -> infinity` in density
(Hardy–Ramanujan; Astra's CRT plus Chebyshev argument is also correct). Letting `k -> infinity` proves the claim.
Correct.

**Plainly:** this concerns typical integers. They have a prime factor of size `q^{1/k}` except on a set of density
`rho(k)` (Dickman), so in fact `B_max(q)` is bounded in probability, which is much stronger than `o(omega(q))`. The
moduli that matter for the conjecture and for the larger-sieve endpoint are smooth with comparable prime factors (and
on circles, the `N` attaining `log R/log log R` growth). These have density zero, and the statement says nothing
about them. Astra says as much ("does not exclude the balanced smooth case").

## 5. Computations (Computed)

**Twin.** `twin.py`: exhaustive over all class pairs for `C = 1`, every `q` in `[2, 400]` (20,717 distinct maximal
configurations). `twin2.py`: exhaustive DFS over maximal configurations with `B >= thr` for
`q in {30030, 9240, 60060, 4620, 46189, 65536, 59049}` (`C = 1`, thr 3), `{30030, 60060, 510510, 360360, 720720}`
(`C = 4`, thr 6) and `{30030, 510510}` (`C = 8`, thr 8). I stopped the `q = 9699690`, `C = 8` run unfinished. Results:
signature map injective whenever `L < q`; **0 violations** of the multiplicity bound `(1+C^2)(1+2C^2)`, of its sharp
form `(1 + floor(2L/K))(1 + floor((L/H)^2/K))`, and of the sumset bound; no affine 3-flat in any configuration found.
These `q` are below `(4C)^6` for `C >= 4`, so IX.1 is not tested there. The largest multiplicities observed are
3 (`C = 4`) and 4 (`C = 8`), against bounds 561 and 8385. Max `B/sqrt D`: 1.51 (`C = 1`), 2.12 (`C = 4`), 2.14
(`C = 8`). Example `C = 4`, `q = 30030`: `S = {0, 105, 182, 495, 572, 677}`, `B = 6`, `D = 8` (six points of a 3-space).

**Circles.** `circle.py`: all `N <= 10^6` composed only of primes `1 mod 4` with `omega >= 3`, at `C = 2, 8`
(104,280 and 323,192 windows). Structured `N`: products of the first `k` split primes (`k = 3..10`), squares
(`k = 2..6`), cubes (`k = 2..4`), `5^6 13^3 17^2 29 37`, `5^4 13^4 17^2 29^2`, and nine primes in `[97, 181]`, at
`C = 1, 4, 16, 64` (skipping `N` with window `>= 1.5` rad). Results:
* gap inequality: 0 violations;
* distinct signature differences in one unit class on subarcs of diameter `c N''^{-1/4}`: 0 violations for
  `c = 1.99, 2.5, 4`; violations occur at `c = 8` (400 at `C = 8`, range; 32 at `C = 16`), so the constant matters;
* `para.py`: over all same-unit signature parallelograms on 13 structured `N`, the minimum of
  `diam * N''^{1/4}` is 4.0038, consistent with the proved `>= 2`;
* `M` against `sqrt D`. Largest `M/sqrt D`: 1.41 (`C = 1`, range) / 1.50 (`C = 1`, structured), 2.00 (`C = 2, 4`),
  2.83 (`C = 8`), 3.25 (`C = 16`), 4.25 (`C = 64`). At `C = 64`, `N = 5*13*17*29*37*41`: `M = 34`, `D = 64 = 2^s`.
  Here `D` is saturated, the bound is trivial, and `M <= 1 + s + K_C` is the binding bound. No cluster with
  `D = o(M^2)` and large `M` appears. Small-doubling configurations are what the theorem excludes, and they do not
  occur.

## 6. Bearing on the project
* IX.1 and IX.3 (twin) are correct new structural lemmas: no affine 3-flats at the critical length, and bounded
  Sidon multiplicity. They exclude the Walsh/Sylvester-Hadamard shape over `Z`, which is the shape of the unbounded
  characteristic-2 clusters. This is useful as an explanation of the char-2 / char-0 split.
* None of them bounds `B` or `M` for worst-case parents. Doing so would need `D = o(s^2)`. The fakes (circle) and,
  per Astra, the IX.4 relaxation (twin) are Sidon, so no local input supplies that.
* The IX.6 criteria are correct. The circle one is in the record; the twin one is its routine transcription. The
  density statement concerns typical, non-smooth moduli.
