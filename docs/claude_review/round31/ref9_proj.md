# Referee report on Part IX, block "proj": IX.2 and the first paragraph of IX.6

Refereed text: Part IX of `astra_free_IX/answer.md`, section IX.2 (lines 1980-2031) and the first paragraph of IX.6
(lines 2231-2236). Scripts (all my own; Astra's scripts were read, not imported):
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref9_proj/`
(`gz.py`, `cert.py`, `fibre.py`, `record10.py`, `scan.py`, `structured.py`, `structured2.py`, `pairmid.py`,
`field.py`, `paley.py`; outputs `*.log`, `s2.json`, `s4.json`, `s9.json`).

## 0. Verdict

All mathematical statements of the block are **correct**. I rederived each step and tested every inequality with
exact arithmetic, with 0 failures. The corrections are minor: one quantifier, one overstatement about field
examples, and one point of scope. The result is a genuine unconditional theorem. Its content, however, is
Lemma 4.2(b) of `round4/hardness.md` (the trivial projection bound near rational directions) improved by
Gallagher's larger sieve, from exponent 1 to 2/3 in the height of the direction, together with a new and useful
translation of "the arc is near a rational direction of small height" into a check on the exponent rows. The
hypothesis is **exceptional**. Exactly symmetric clusters satisfy it (these include every record cluster with
`M >= 8`), and there the bound is classical. Generic pairs have `T >= N'^{c}`, and the data of Theorem 1.1's
extremal profile admit no small certificate at all. The theorem does not move the growth rate, the constant 2/3
or the deep window as a whole. What it does prove is a restricted rate, in a class where the classical
observation gives only `P_+ >> M`.

| Claim | Label |
|---|---|
| (IX.1) fibre counts, lifts, larger-sieve inequality, Mertens step, `2+X_b >> M^2 exp(-sum) >> M^{3/2}`, absolute constants | **Proved** (checked line by line) |
| (IX.2) factorisation, units, `eta`, angle of `b`, derivative bound, `(1+||c||_1) sqrt T >>_C M^{3/2}` | **Proved** |
| `min (1+||c||_1)^2 T = o(s^3)` implies `M = o(s)` | **Proved** (along any sequence with `s -> infinity`) |
| Pair corollary `M <<_C 1 + P_+^{2/3}`; bounded `M` for `P_+ <= M^{3/2-eta}`; squarefree reading | **Proved** |
| Orientation criterion and single-mark criterion (IX.6, first paragraph) | **Proved** |
| "Tests and scope": Paley fakes have `log T >> q log q` for every `c` with sum 2 | **Proved**; sharp at `(q-1)/2` nonzero coordinates (Computed) |
| Function-field check `p^{deg T} >>_p M^3` | **Proved** for constant quadratic extensions. "Every field example, including ... VII.1" is an **overstatement**, because VII.1's second model `eta^2 = Delta(t)` is not a constant extension |
| COMPUTED paragraph (6,079 fibres; 10-point example; five `T=1` pairs; projections 48503...48515) | **Computed**, reproduced independently |
| Field COMPUTED (D=12 ternary witness: degrees 2, 4, 6 occurring 1, 4, 5 times; four traces per pair) | **Computed**, reproduced |

## 1. The integer projection theorem (IX.1)

**Setting.** Let `z_1..z_M` be the normalised marks with `|z_i|^2 = N'`, where `N'` is odd and all of its primes
are split. Fix `0 != b in Z[i]`. Put `w_i = conj(b) z_i`, `x_i = Re w_i` and `D_b = N'|b|^2`. Let `n` be the number
of distinct `x_i`, and `X_b = max x - min x`.

* **Lifts.** `|w_i|^2 = D_b`, so `x` determines `Im w = ±sqrt(D_b - x^2)`, and `z = w/conj(b)`. Hence each value
  of `x` lifts to at most two marks, and `M/2 <= n <= M`. Correct.
* **Fibre counts.** For odd `p` not dividing `D`, the conic `x^2+y^2 = D` has `p - chi_4(p)` points mod `p`.
  Each `x` with `D - x^2` a nonzero square gives two points, and each `x` with `D = x^2` gives one. So the number of
  first coordinates is `nu_p = (p - chi_4(p) + 1 + (D/p))/2`, which is at most `(p+3)/2`. For an inert
  `p | D`, only `x = 0` occurs, so `nu_p = 1`. For a split `p | D`, all `x` occur, so `nu_p = p`. For `p = 2`,
  `nu_2 <= 2`, which is harmless. **Computed:** all 6,079 pairs `(p, D mod p)` with `3 <= p <= 251` agree with the
  formula, including `D = 0`. The point count `p - chi_4(p)` was checked for `p < 60` (`fibre.py`).
* **Larger sieve (derivation).** Fix any `Q`. For each prime `p <= Q`, Cauchy-Schwarz over the `nu_p` classes
  shows that the number of ordered pairs `i != j` of distinct projections with `x_i = x_j mod p` is at least
  `n^2/nu_p - n`. This holds even when the bound is negative. The distinct primes dividing a difference divide it,
  so
  `sum_{p<=Q} log p (n^2/nu_p - n) <= sum_{i!=j} log|x_i - x_j| <= n(n-1) log X_b`.
  Dividing by `n^2` gives `(1-1/n) log X_b >= sum_{p<=Q} log p/nu_p - theta(Q)/n`.
  - The inequality is valid for **every** `Q`. The choice `Q = M` is correct and natural: `theta(M)/n <= 2 theta(M)/M = O(1)` because `n >= M/2`. Terms with `nu_p > n` (possible for `p` near `M` when `n ~ M/2`) are negative but harmless.
  - `log X_b` is what the sieve gives. Replacing it by `log(2+X_b)` is legitimate because `1 - 1/n >= 0`. The inequality needs `n >= 2`, which holds whenever `M >= 3`; for `M <= 2` the statement is trivial.
* **Mertens step.** At good odd primes, `log p/nu_p >= 2 log p/(p+3) = 2 log p/p - O(log p/p^2)`. Split bad primes
  contribute `log p/p`, and inert bad primes contribute `log p >= 0`. So the right side is at least
  `2 log M - sum_{p<=M, p|D_b, p=1(4)} log p/p - O(1)`. When it is positive, `log(2+X_b) >= (1-1/n)log(2+X_b) >=` the
  right side. The worst case is `sum_{p<=M, p=1(4)} log p/p = (1/2)log M + O(1)`. Hence
  `2 + X_b >> M^2 exp(-sum) >> M^{3/2}`, with **absolute** constants. Correct as stated.
* **Optimality within the method.** If every split `p <= M` divides `N'`, then `sum_{p<=Q} log p/nu_p = (3/2)log Q + O(1)`.
  Optimising `Q` then gains at most logarithms, so the exponent 3/2 is the method's limit. If no small split prime
  divides `D_b`, the same proof gives `2 + X_b >> M^2`. Theorem 1.1's extremal data have `e_p = 1` at every split
  `p <= M`, so the worst case is the relevant one.

## 2. The exponent-row certificate (IX.2)

**Factorisation (Proved).** Write `z_i = eps_i prod_p pi_p^{a_ip} conj(pi_p)^{e_p - a_ip}` and `N' = prod (pi_p conj(pi_p))^{e_p}`.
Let `c in Z^M` with `sum c_i = 2`. Then `prod z_i^{c_i} / N' = u prod_p (pi_p/conj(pi_p))^{v_p}`, where
`v_p = sum c_i a_ip - e_p` and `u = prod eps_i^{c_i}`. This uses `sum c_i(e_p - a_ip) - e_p = -v_p`, and it holds
for negative `c_i` because `z^{-1} = conj(z)/N'`. Put `alpha = prod_{v>0} pi^{v} prod_{v<0} conj(pi)^{-v}`. Then
`alpha/conj(alpha) = prod (pi/conj pi)^{v_p}`, `|alpha|^2 = T(c)`, and `alpha` is primitive.
- **Units.** `eta/conj(eta)` equals `1, -1, i, -i` for `eta = 1, i, 1+i, 1-i` respectively. So `b = eta alpha` gives `b/conj(b) = prod z^c/N'`, and `|b|^2 in {T, 2T}`.
- **Bad primes.** The odd split bad primes of `D_b` are exactly those of `N'`.
- **Computed:** this identity was checked exactly, with the unit found uniquely, for every certificate in the scans below (pairs, `2e_i`, and every `c in {-1..2}^M` with `||c||_1 <= 4` for `M <= 9`).

**Angle (Proved).** `b/conj(b) = e^{2i arg b}`, so `2 arg b = sum c_i arg z_i (mod 2pi)`. This holds for any real lifts,
because a change of lift of `phi_i` by `2pi` moves the half-sum by `c_i pi`. Write `phi_i = theta + delta_i` with
`|delta_i| <= Delta/2`, where `theta` is the midpoint of the marks' own angular span. Then `sum c_i = 2` gives
`(1/2) sum c_i phi_i = theta + (1/2) sum c_i delta_i`. The distance from `theta` mod `pi` is therefore at most
`||c||_1 Delta/4`. No phase is set to zero: `theta` is the actual midpoint.

**Derivative bound (Proved).** `x_i = |b| R' cos(phi_i - arg b)`. Let `psi = theta - arg b`. The mean value
theorem on an interval of length `Delta` about `psi` gives `X_b <= |b| R' Delta (|sin psi| + Delta/2)`. The angular
width satisfies `Delta <= C sqrt R / R = C/sqrt R <= C/sqrt R'`, since `R = |g| R' >= R'`. So `R' Delta^2 <= C^2`,
and `X_b <= C^2 |b| (||c||_1/4 + 1/2) <<_C |b|(1 + ||c||_1)`. Combined with (IX.1) and `|b| <= sqrt(2T)`:
`(1+||c||_1) sqrt T >>_C M^2 exp(-sum_{p|N', p<=M} log p/p) >> M^{3/2}`. The implied constant is of the shape
`c_0/C^2` with `c_0` absolute. **Computed:** the angle claim and the derivative bound hold for every tested
certificate, with 0 failures.

## 3. Corollaries

* **`o(s^3)` implies `M = o(s)`. Proved.** `M^3 <<_C (1+||c||_1)^2 T = o(s^3)`. The rate statement then follows from
  `s << log N/log log N`. Quantifier: this holds along any sequence of clusters with `s -> infinity` (if `s` stays
  bounded, `M` is bounded by Theorem 1.2 anyway). The constants depend on `C` only.
* **Pair corollary. Proved.** With `c = e_i + e_j` we have `||c||_1 = 2` and `T = prod p^{|a_ip+a_jp-e_p|} <= P_+^{sum|v_p|} <= P_+^2`,
  so `M <<_C 1 + P_+^{2/3}`. In the squarefree case `|a_ip + a_jp - 1| = 1` exactly where the rows agree, so the
  hypothesis says the rows agree in at most two columns. Correct.
* **Bounded `M` for `P_+ <= M^{3/2-eta}`. Proved.** We get `M <<_C M^{1-2eta/3}`, so `M <= M_0(C, eta)`.
  - Comparison with Theorem 1.2: `P_+ >= (2-o(1)) M log M` holds for every cluster, so the window `2M log M <= P_+ <= M^{3/2-eta}` is nonempty. In the pair class IX.2 gives the stronger `P_+ >>_C M^{3/2}`. Inside that class this is a genuine improvement on Theorem 1.2.
  - The trivial projection bound gives only `P_+ >>_C M` there. That is weaker than Theorem 1.2, so the larger-sieve gain is exactly what makes the pair corollary nonvacuous.
  - Plainly: for large `M`, pair-class clusters with `P_+ <= M^{3/2-eta}` do not exist. This includes `P_+ <= M(log M)^A`.
* **Orientation and single-mark criteria (IX.6). Proved.** `2+X_b <<_C |b|(1 + sqrt(R')|sin(theta - arg b)|)`
  follows from the derivative bound, using `R' Delta <= C sqrt R'` and `R' Delta^2 <= C^2`. For `c = 2e_i`:
  `T = prod p^{|2a_ip - e_p|} = N'/gcd(Re z_i, Im z_i)^2`, because `e - 2 min(a, e-a) = |2a - e|` and the
  rational content of `z_i` is `prod p^{min(a_ip, e_p - a_ip)}`.

## 4. Meaning, novelty and how typical the hypothesis is

**What the certificate says (Proved, and Computed for pairs).** `b/conj(b) = prod z^c/N'`, and `arg b` lies within
`||c||_1 Delta/4` of the arc. A small certificate is therefore a rational direction of small height (`|b| ~ sqrt T`)
that lies essentially on the arc, and the theorem bounds clusters lying at such directions.
- For a pair: `z_j = (b/conj b) conj(z_i)`. So `z_j` is the reflection of `z_i` in the line through `b`, and `z_i + z_j` is a real multiple of `b`. **Computed:** `|b|^2 = |primitive part of (z_i + z_j)|^2` for all 281,400 pairs tested (`pairmid.py`).
- So the pair hypothesis says that the chord midpoint `z_i + z_j` has rational content at least about `R'/s^{3/2}`.
- `T = 1` means exact symmetry about an axis or a diagonal.

**A general upper bound for pair certificates (Proved here, one line).** For `a, a' in [0, e]` one has
`|a + a' - e| + |a - a'| <= e`. So `log T(e_i+e_j) <= W - d_ij`. Combined with the pair inequality
`d_ij >= W/2 + 2 log|b_ij| - 2 log C'`, this gives `T(e_i+e_j) <= C'^2 sqrt(N')/|b_ij|^2` for every pair. So pair
certificates always have `T <= C^2 sqrt(N')`. The hypothesis of the corollary asks for `T <= P_+^2`, which is
polylogarithmic in `N'` in the window.

**Novelty.**
- Lemma 4.2(b) of `round4/hardness.md` already gives `M <= 2|w|(kappa + C)^2 + 2` for arcs within `kappa R^{-1/2}` of `arg w`. That is the trivial bound `M <= 2(X_b + 1)`, i.e. `M <<_C |b|`. Round 4 Prop 2.3(2) is the case `b = 1 - i`. Round 4 Remark 4.8 already notes that such directions cover a negligible set of arcs.
- IX.2 adds two things. The first is the larger sieve on the projected norm fibre: `M <<_C (|b|(1+||c||_1))^{2/3}`, and `<< X_b^{1/2}` when no small split prime divides `D_b`. The second is the exponent-row form of the certificate.
- For bounded `|b|` (axes and diagonals) both bounds are `O_C(1)`, and nothing is new except the power of `C`.
- The gain matters only when `|b|` grows. There, `M << |b|^{2/3}` against the trivial `M << |b|` is exactly what yields the pair corollary.
- Literature: a single web search found no source for the larger sieve applied to `{x in I : D - x^2 is a square}` on short arcs. I saw only snippets; the nearest hits were [Pandey, squarefree numbers in short intervals](https://arxiv.org/pdf/2401.13981), [the Wikipedia entry on the larger sieve](https://en.wikipedia.org/wiki/Larger_sieve) and [a paper on lattice points on short arcs](https://arxiv.org/pdf/2012.10784). Novelty is unverified beyond this. I expect the sieve step itself to be folklore.

**How typical (Computed and Heuristic).**
- *Small circles (Computed; uninformative about asymptotic typicality).* All maximal windows on all circles `N <= 10^6` were scanned, deduplicated under the 8 lattice symmetries:
  - `C = 2`, `M >= 3`: 1,175 windows;
  - `C = 4`, `M >= 3`: 12,067 windows;
  - `C = 9`, `M >= 5`: 6,413 windows.
  - An exactly symmetric pair (`T = 1`) occurs in 27-46% of `M = 3` windows, 74% of `M = 4` windows at `C = 4`, and 60-100% of `M >= 5` windows at `C = 9` (94% at `M = 8`, 9 of 9 at `M = 9`).
  - The hypothesis `sum|v| <= 2` held in every window. This is forced at these sizes by the bound `T <= C^2 sqrt N'` with `s <= 7`, so it says nothing about large `N'`.
- *Records (Computed).* The 10-point cluster (`N = 1176852625`), the golden 8-point cluster (`537606725`) and the 8-point record (`32988780417125`) have 5, 4 and 4 pairs with `T = 1`. They are symmetric about a diagonal or an axis, where the classical bound already gives `O_C(1)`. The sporadic 7-point cluster (`235370525`, not symmetric) has a pair with `T = 5` (direction of height `sqrt 5` or `sqrt 10`).
- *Generic large circles (Computed).* I took `N'` = product of the first `k = 6..10` primes `= 1 (mod 4)` (`N'` up to `10^15`) and all windows at `C = 2, 4, 8` (`structured2.py`):
  - the minimum over pairs of `log T/log N'` has median 0.18-0.29 and minimum 0.07-0.26 for `k >= 7`;
  - the minimum over all `c` with `||c||_1 <= 4` of `log((1+||c||_1)^2 T)/log N'` has median 0.17-0.39;
  - no symmetric pair occurs for `k >= 7`.
  - So generic certificates have `T ~ N'^{0.2-0.3}`, which is exponential in `s`, against the required `o(s^3)`.
- *Heuristic.* Directions of height `<= H` with `H = s^{3/2}` occupy a fraction `~ s^3/sqrt R'` of the angle, which tends to 0. Small certificates are exceptional outside reflection-symmetric configurations.
- *Theorem 1.1's extremal profile (Proved).* Take the data of `round31/ref8_extend.md` Prop 3.2: `e_p = 1` at the split `p <= M`, plus `M` primes near `M^3` whose balanced 0/1 columns form a basis `B`. Vanishing of `v` at all large primes forces `B^T c = 1`. Its unique solution is `c = (2/M) 1`, which is not integral for `M > 2`. So every `c` has `T(c) >=` the least large prime `>> M^3 ~ s^3`, and the hypothesis fails. (IX.2) requires `T >> M^3/(1+||c||_1)^2` anyway, so IX.2 is consistent with that profile and does not touch the constant 3 or 2/3.

## 5. Numerical tests of (IX.1) and (IX.2) on genuine clusters

For every window of the three scans and every `b` in a test set I checked the exact chain
`sum_{p<=M} log p(n^2/nu_p - n) <= sum_{i!=j} log|x_i - x_j| <= n(n-1) log X_b`. The test set was: all primitive `b`
with `|b|^2 <= 50`, plus the `b` of every tested certificate. I also checked `n >= M/2`, the angle claim and the
derivative bound. **0 failures** in all scans.

| scan | min `(2+X_b)/M^{3/2}` | min `(2+X_b)/(M^2 e^{-sum})` | min `(1+‖c‖_1) sqrt(T)/M^{3/2}` |
|---|---|---|---|
| `N <= 10^6`, `C = 2` | 0.500 | 0.250 | 0.577 |
| `N <= 10^6`, `C = 4` | 0.358 | 0.188 | 0.268 (N = 325, M = 5, T = 1) |
| `N <= 10^6`, `C = 9` | 0.265 (N = 5, the whole circle) | 0.111 | 0.111 (N = 65, M = 9, T = 1) |
| 10-point record, `C' = 8.24` | 0.443 | 0.19 | 0.095 |

The ratios decrease with `C`, as `>>_C` allows (the constant is about `C^{-2}`). On the record, `(2+X_b)` multiplied
by `C^2/M^{3/2}` is about 6. **10-point example reproduced:** the Gaussian gcd is a unit (the record's points
divided by `-i`). Five pairs `(0,9), (1,8), (2,7), (3,6), (4,5)` have complementary rows and `T = 1`. With
`b = 1 - i` in my rotated frame (`1 + i` in Astra's frame) the projections are `-48515, -48511, -48509, -48505, -48503`
(her values, with sign changed by the rotation). These are the only `x` in that interval with `2N' - x^2` a
square. The sieve chain there reads `13.44 <= 32.80 <= 49.70`.

## 6. The function-field check

**Identity (verified).** `x_i - x_j = conj(b)(z_i - z_j)(1 - N b/(conj(b) z_i z_j))`, expanded using `N/z = conj z`.
`N b/(conj(b) z_i z_j) = prod_l z_l^{c_l}/(z_i z_j)` has exponent sum 0. It is therefore a product of ratios in the
group `U_K`, so `1 - (...)` has degree `<= -K`. Hence `deg(x_i - x_j) <= deg b + (D - K) - K <= deg b` for
`K >= D/2`. Also `deg T = 2 deg b`, because the extension is constant (Hilbert 90 for the unit).

**Sieve.**
- Inert fibres have at most `(Q+3)/2` traces, and split fibres at most `(Q+1)/2`.
- `sum_{deg lambda <= k} deg lambda/|lambda| = k + O_p(1)`, of which `k/2 + O_p(1)` comes from even degrees (the split primes).
- The worst case gives `(3/2)k - O_p(1) <= deg b`, so `p^{deg T} >> (M/p)^3 p^{-O(1)} >>_p M^3`. Correct.
- Ultrametricity removes the factor `(1 + ||c||_1)`.

**Relation to Astra's field families.** There is no contradiction with VI.5 or VII.1. The statement only
lower-bounds `T`, and those families have large `T` for every `c` (for VI.5 Astra's orthogonality argument shows
this; I did not recheck it). **Correction:** the proof covers the constant extension `F_{p^2}(t)/F_p(t)`. VII.1's
second model `eta^2 = Delta(t)`, with a ramified finite prime, is not covered as written. "A direct test on every
field example, including ... VII.1" should read "on every field example over a constant quadratic extension
(VI.4, VI.5 and the first model of VII.1)". An adaptation via Chebotarev density 1/2 for the split primes is
plausible but not written (**Gap**, minor).

**Computed (`field.py`, F_9[t] arithmetic via python-flint as a library).** I took the D = 12 five-mark ternary
witness (monic, `deg N = 24`, pair difference degrees `<= 6 = D/2`).
- For each of the 10 pairs: `b` with `z_i z_j conj(b) = N b`, `deg b` in `{2, 4, 6}` occurring `{1, 4, 5}` times, four distinct traces, and maximum trace-difference degree equal to `deg b`. This reproduces her paragraph.
- In addition, all 155 sum-two certificates `c in {-1..2}^5` satisfy `deg(x_i - x_j) <= deg b`.

## 7. "Tests and scope": Paley fakes

The claim is correct. With `a = 1 + K` and `e = 2`, `v = K_S^T c` on `I`. Modulo `q`,
`(K^T c)_j = sum c_i (j - i)^{(q-1)/2}`. This is a polynomial of degree `(q-1)/2` in `j` with leading coefficient
`sum c_i = 2`, so it is nonzero at `>= (q+1)/2` points, and at `>= (q-1)/2` after deleting column 0. Hence
`log T >= ((q-1)/2) log P_0 >> q log q`. **Computed (`paley.py`):** `q = 7, ..., 43`, with all pairs, all `2e_i`
and 20,000 random integral `c` with sum 2. The minimum number of nonzero `v_j` is exactly `(q-1)/2`, so the bound is
sharp. The remark that round 5's point residues provide no integer interval is heuristic and fair: IX.6's fake in
the `2+i` window shows that the orientation criterion uses genuine integrality.

## 8. Corrected statements in full

* **(IX.1).** Let `M >= 3` genuine normalised marks and `0 != b in Z[i]`. Then
  `log(2+X_b) >= sum_{p<=M} log p/nu_p - 2 theta(M)/M`, with `nu_p` as above. Consequently
  `2 + X_b >= c_0 M^2 exp(-sum_{p<=M, p|D_b, p=1(4)} log p/p) >= c_1 M^{3/2}`, with absolute `c_0, c_1 > 0`. Unchanged in substance.
* **(IX.2).** With `C` the arc constant, `(1+||c||_1) sqrt(T(c)) >= c_2 C^{-2} M^{3/2}` for every integral `c` with
  `sum c_i = 2` and `M >= 3`, where `c_2` is absolute. Unchanged.
* **Function-field check.** The check holds in `F_{p^2}[t]/F_p[t]` (constant extension). It applies to the field
  examples built over that extension. It does not apply as written to the ramified model of VII.1.
* **Scope sentence (to add).** The hypothesis of IX.2 holds only for clusters lying within `O(||c||_1 Delta)` of a
  rational direction of height `O(sqrt T)`, with `T = o(s^3)`. Pair certificates always satisfy `T <= C'^2 sqrt N'`.
  Generic clusters have `T = N'^{Omega(1)}`, and the extremal data of Theorem 1.1 have `T >> s^3` for every `c`.
  For bounded height the result is the classical axis and diagonal observation (round 4 Lemma 4.2(b)).
