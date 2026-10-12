# Stability of near-extremal clusters: the norm is `(log R)^(O(1))`-smooth, and the structure is still fake-consistent

Round 8. Target: a rigorous structure theorem for clusters with `M >= eps log R/loglog R` lattice points on an
arc of length `C sqrt R` of `x^2 + y^2 = R^2`, and an attempt to turn the structure into a growth improvement.
Scripts: `round8/stability_checks/` (each runs from that directory with `python3 <script>`; outputs `*_out.txt`).

Conventions. *Theorem/Proposition/Lemma*: proved here with all quantifiers, or quoted with an exact source.
*Exact check*: integer arithmetic (logarithms at 60 digits only where an inequality between logarithms is
checked). *Evidence*: a computation that no proof uses.

## 0. Verdict

**No growth improvement for general circles. The stability theorem is proved in a strong form: every cluster
lives on a norm whose prime factors are all `< N^(2/M) C_+^4`, so a near-extremal cluster lives on a
`(log R)^(O(1))`-smooth norm with `Theta(M)` prime factors. This gives an unconditional growth theorem for every
circle whose norm is not `(log R)^(O(1))`-smooth. On smooth norms the whole structure (and more) is satisfied by
the refereed fake circles of `round5/construct.md` at `W = O(M log M)`, so growth there needs the actual angles of
Gaussian primes of norm `(log R)^(O(1))`, beyond what proved linear-forms bounds give.**

| # | statement | type |
| --- | --- | --- |
| Prop 2.2 | **gcd lemma:** any `s` points of a cluster have a common Gaussian divisor of norm `< N^(1/s) C_+^4`; precisely `log Norm gcd <= W/s + (4(s-1)/s) log C - 8F*(s)/s^2`. Extends the notes' obtuse Bessel bound (`C <= sqrt 2`, cut classes only) to every `C` and every subset, with the arithmetic gain `~ 3 log s` | proved, new in this form |
| Thm 3.1, Cor 3.2 | **smoothness:** every prime `q | N` has `q < N^(2/M) C_+^4` (and `log q <= (2/M) log N - 3 log(M/2) + K_0 + 4 log^+ C`). Hence `M < 2 log N/log(P_+(N)/C_+^4)` | proved, new |
| Cor 3.2 (consequence) | **restricted growth:** along any sequence of circles with `log P_+(N)/loglog R -> infinity`, `M = o(log R/loglog R)`; a prime factor `>= N^lambda C_+^4` gives `M < max(4, 2/lambda)` | proved |
| Thm 3.3, Prop 3.4-3.5 | prime powers `p^e || N'` are `< N'^(4/M) C_+^8`; `omega(N') > M/4` (`C <= 1`, all `M >= 4`) and `>= (M/2)(1 - M^(-1/3))^2`; `> M/2` distinct cut partitions; any two rows differ in `> M/4` cuts; primes `<= M^(1-delta)` carry `O(M^(-delta)/log M)` of `W'` | proved |
| Cor 3.6 | near-extremal (`M >= eps log R/loglog R`): `P_+(N) < C_+^4 (log R)^(4/eps)` (indeed `(log R)^(4/eps - 3 + o(1))`); `M/4 < omega(N') <= Omega(N') <= (2/eps + o(1)) M` | proved |
| Prop 4.1 | Plotkin: average margin in `(log M - 7.3 - 2 log^+ C, W'/(2(M-1))]` (`M >= 120`), i.e. `Theta(log M)`; cut slack `Q < MW'/4 - (M^2/2) log M + O(M^2)`; geometric mean of `|beta_xy|^2` in `[M^(1-o(1)), C^2 M^(1/eps + o(1))]` | proved |
| Prop 5.1, Thm 5.2 | collision, cofactor and level excesses each `<= (A-3)/3 + O(1/log M)` of the main term; near the 2/3 bound (`A <= 3 + delta`) all but a `delta`-fraction of split primes `<= M` divide `N'` | proved |
| Prop 6.1 | **genuine input (Matveev + smoothness):** no nonzero character of a near-extremal cluster is supported on `<= (1 - eta) log M/loglog M` primes | proved (Matveev quoted) |
| Prop 7.1, Cor 7.2 | everything in Sections 2-5 holds for every fake circle in `A_1^>(5M/4)`; the refereed fakes (`W <= 27 M log M`) satisfy all of it; Prop 6.1's conclusion holds for the `construct.md` profiles (Remark 6.2(3)) and the skew double core | proved / barrier |
| Section 8 | three-term identities (Lemma 8.1), S-unit counting (Prop 8.2: bounds are `2^(Theta(M))` vs `M^3` identities), small-`S` relations = parallel chords (Prop 8.3), abc, smooth numbers near squares, GRH counting: none gives growth; exact reasons given | proved / heuristic (8.4) |

**What is new relative to the notes and earlier rounds.** The class bound `V_S <= 2W/M` for `C <= sqrt 2` is the
notes' obtuse Bessel inequality; the general sub-cluster inequality, its arithmetic gain, the treatment of the
content `G` and of inert primes, the smoothness theorem for all `C`, the restricted growth corollary, the prime-count
bounds, the stability bounds of Sections 4-5 and Propositions 6.1, 7.1, 8.2, 8.3 are, to my knowledge, not stated
elsewhere in this review or in the notes. The literature could not be searched in this session.

---------------------------------------------------------------------------------------------

## 1. Setting, notation, and the inputs used

**Configurations.** `N = R^2 >= 2` is an integer and `W = log N = 2 log R`. A *configuration* is a set of
`M >= 2` distinct Gaussian integers `z_1, ..., z_M` with `|z_x|^2 = N` and

```text
|z_x - z_y| <= C sqrt R      for all x, y.                                            (1.1)
```

This chord bound is the only geometric input; it holds for points on an arc of length `<= C sqrt R`.
`log^+ C = max(0, log C)`, `C_+ = max(1, C)`, `K_s = binom(s, 2)`.

**Normalisation** (`two_thirds/proof.md` Lemma 1.1). `G = gcd(z_1, ..., z_M)`, `u_x = z_x/G`,
`N' = N/|G|^2`, `W' = log N' <= W`, `R' = sqrt N'`. The `u_x` are distinct, `gcd(u) = 1`, `N'` is odd and has
no prime `= 3 mod 4`, and `u_x = eps_x prod_{p | N'} pi_p^(a_xp) conj(pi_p)^(e_p - a_xp)` with
`min_x a_xp = 0`, `max_x a_xp = e_p`. The `u_x` satisfy (1.1) on the circle of radius `R'` with the constant
`C' = C |G|^(-1/2) <= C`, because `|u_x - u_y| = |z_x - z_y|/|G| <= C sqrt R/|G| = C' sqrt R'`.

**Cuts.** A *cut* is a pair `(p, tau)` with `p | N'`, `1 <= tau <= e_p`; its weight is `w = log p`, so the cut
weights sum to `W'`. Its sides are `{x : a_xp < tau}` (size `S_(p tau)`) and `{x : a_xp >= tau}`; by
normalisation `1 <= S_(p tau) <= M - 1`. The *partition* of a cut is the unordered pair of its sides.
`r = Omega(N') = sum e_p` is the number of cuts and `omega(N')` the number of distinct primes.
`Q = sum_(cuts) w (S - M/2)^2` is the cut slack.

**Pair data** (two_thirds Section 2). `g_xy = gcd(z_x, z_y)`; the *reduced chord* is
`beta_xy = (z_x - z_y)/g_xy in Z[i] \ {0}` (it is `c_xy` in `two_thirds`; it does not change under
`z -> u`). `d_xy = sum_p |a_xp - a_yp| log p`, so `log |g_xy/G|^2 = W' - d_xy`. The *margin* is
`m_xy = d_xy - W'/2`, and `C_xy = |u_x - u_y|/sqrt R' <= C'`.

**Sub-configurations.** For `B` a set of `s >= 2` indices: `G_B = gcd_(x in B) z_x`, `gamma_B = log |G_B|^2`,
`C_B = max_(x,y in B) |z_x - z_y|/sqrt R <= C`.

**Inputs quoted from `two_thirds/proof.md`** (refereed twice, `README.md` Section 2).

* (T1) Theorem 5.2 together with (2.3): for every configuration of `s >= 2` points with chord constant `c` whose
  gcd is `1` and whose norm is `n`,
  `F*(s) <= sum_(pairs) log|beta| + Q/2 = sum_(pairs) log|z_x - z_y| - (s(s-2)/8) log n <= (s/8) log n + K_s log c`.
  Here `F*(s) > 0` is the explicit configuration-free quantity of Theorem 5.2 (`stab_lib.Fstar` evaluates it).
* (T2) From the proof of Theorem 7.5 (step 2): `8 F*(s)/s >= 8 D(s)/s > 3 s log s - K_0 s + 3 + 2 log 2` for
  every `s >= 2`, with `K_0 < 13.9092`. Also `F*(s) >= K_s (log 2)/2`, from the ramified prime.
* (T3) Theorem 7.5: `3 M log M <= W' + (K_0 + 4 log^+ C) M`.

Write `K' = K_0 + 4 log^+ C` and call `A = W'/(M log M)` the *efficiency*; (T3) says `A >= 3 - K'/log M`.

**The near-extremal hypothesis.** `H(eps)`: `R >= e^e` and `M >= eps log R / loglog R` (`0 < eps <= 1`).
Under `H(eps)`,

```text
W' <= W = 2 log R <= (2/eps) M loglog R,     so     2W/M <= (4/eps) loglog R.               (1.2)
```

Moreover `log M >= loglog R - logloglog R + log eps` (take logarithms in `H(eps)`), and by (T3)
`log M <= loglog R` once `R >= R_0(C)`. Hence `loglog R = (1 + o(1)) log M` as `R -> infinity` with `eps, C`
fixed, and `A <= (2/eps)(1 + o(1))`. All `o(1)` below refer to this limit; every inequality is also stated
in a non-asymptotic form.

---------------------------------------------------------------------------------------------

## 2. The sub-cluster inequality

**Lemma 2.1 (exact sub-cluster identity).** Let `z_1, ..., z_M` be any configuration (normalised or not) and
`B` a set of `s >= 2` indices. Put `u^B_x = z_x/G_B` (`x in B`), a configuration of norm `N_B = N/|G_B|^2` with
gcd `1`, and let `Q^B` be its cut slack. Then

```text
sum_(x<y in B) log|beta_xy| + Q^B/2 + (s^2/8) gamma_B = sum_(x<y in B) log|z_x - z_y| - (s(s-2)/8) W.      (2.1)
```

Equivalently, in integers,
`prod_B |z_x - z_y|^8 = N^(s(s-2)) Norm(G_B)^(s^2) prod_(cuts of u^B) p^((2S - s)^2) prod_B |beta_xy|^8`.

*Proof.* `gcd(u^B_x, u^B_y) = g_xy/G_B`, so the reduced chord of `u^B` at `(x, y)` is `beta_xy`. The identity
(2.3) of `two_thirds`, applied to the normalised configuration `u^B` (radius `R_B = sqrt N_B`), reads
`sum log|beta| + Q^B/2 = sum log|u^B_x - u^B_y| - (s(s-2)/8) log N_B`. Substitute
`log|u^B_x - u^B_y| = log|z_x - z_y| - gamma_B/2` and `log N_B = W - gamma_B`; the coefficient of `gamma_B`
becomes `K_s/2 - s(s-2)/8 = s^2/8`. QED

For a normalised configuration (`G = 1`) the cuts of `u^B` are the cuts of `N'` that are non-constant on `B`,
and `gamma_B = sum_(cuts constant on B) w`: at `p`, the constant cuts are `tau <= min_B a_xp` and
`tau > max_B a_xp`, and `G_B` has `p`-part `pi^(min a) conj(pi)^(e - max a)`.

*Exact check* (`check_subcluster.py`): the integer form of (2.1) holds on all 148,370 subsets tested (all
subsets of windows of `M <= 6` consecutive lattice points with chord constant `<= 6` on every circle
`N <= 3*10^4`; see the script for the window selection).

**Proposition 2.2 (gcd lemma).** For every configuration satisfying (1.1) and every `B` with `|B| = s >= 2`,

```text
gamma_B <= b(s) := W/s + (4(s-1)/s) log C_B - 8 F*(s)/s^2.                                  (2.2)
```

Consequently `gamma_B < W/s + 4 log^+ C` and `gamma_B < W/s - 3 log s + K_0 + 4 log^+ C`. In words: **any
`s` points of a cluster have a common Gaussian divisor of norm `< N^(1/s) C_+^4`.**

*Proof.* The configuration `u^B` has `s` points, gcd `1`, norm `N_B = e^(W - gamma_B)` and chord constant
`c_B = C_B e^(-gamma_B/4)`, since `|u^B_x - u^B_y| = |z_x - z_y|/|G_B| <= C_B sqrt R/|G_B| = c_B N_B^(1/4)`.
By (T1), `F*(s) <= (s/8)(W - gamma_B) + K_s (log C_B - gamma_B/4)`. The coefficient of `gamma_B` is
`s/8 + K_s/4 = s^2/8`, which gives (2.2). For the consequences use `C_B <= C`, `(s-1)/s < 1`, `F*(s) > 0`, and
(T2): `8F*(s)/s^2 > 3 log s - K_0`. QED

**Remarks 2.3.**

1. *Sharpness.* Equality in (2.2) occurs for `s = 2` exactly when `|beta_xy|^2 = 2` (no odd collision): then
   `b(2) = 2 log|z_x - z_y| - log 2 = gamma_B`. The exact checks meet this equality case (slack `0` to 60 digits)
   and no violation. For `B = [M]`, (2.2) says `|G|^2 < N^(1/M) C_+^4`: the content of a cluster is tiny, a much
   stronger statement than `W' <= W`.
2. *Relation to the notes.* For `C <= sqrt 2` with arbitrary units (or `C <= 2` with one common unit) the notes
   prove the *obtuse Bessel* bound `V_S <= 2W_var/M` for the weight `V_S` of every whole cut class
   (`isolable_pattern_bessel_phase_gap.md`, display (6); `GaussianChain/ObtuseBessel.lean`). Proposition 2.2 with
   `B` the larger side of the class gives `V_S <= 2W'/M - 3 log(M/2) + K'` for **every** `C`, and for every
   subset `B`, not only cut sides. It uses the arithmetic of `F*` instead of obtuseness (pair margins can be
   negative when `C > sqrt 2`).
3. *What it is.* Dividing `s` points by their gcd puts them on a circle of radius `R/|G_B|` with the chord
   constant shrunk by `|G_B|^(-1/2)`; Proposition 2.2 is the 2/3 inequality for that smaller cluster. Applied to
   all subsets at once, it is the complete content of "recursion on sub-clusters" (Section 8.5).

---------------------------------------------------------------------------------------------

## 3. Arithmetic structure: smoothness, prime powers, prime counts

Throughout, `b_C(s) = W/s + (4(s-1)/s) log C - 8F*(s)/s^2` and `b'_C(s)` is the same with `W'` in place of `W`.
By (T2) and `F*(s) > 0`, for `s >= M/2`:

```text
b_C(s) < 2W/M + 4 log^+ C,        b_C(s) < 2W/M - 3 log(M/2) + K'.                        (3.0)
```

**Theorem 3.1 (the norm is smooth).** Let `z_1, ..., z_M` (`M >= 4`) satisfy (1.1). Then every rational prime
`q | N` satisfies

```text
log q <= max_(ceil(M/2) <= s <= M) b_C(s)  <  min{ 2W/M + 4 log^+ C ,  2W/M - 3 log(M/2) + K' }.      (3.1)
```

That is, **`P_+(N) < N^(2/M) C_+^4`** for the largest prime factor `P_+(N)` of `N` (content included), and
more precisely `log P_+(N) <= (2/M) log N - 3 log(M/2) + K_0 + 4 log^+ C`.

*Proof.* Normalise (Section 1). Two cases.
1. `q | N/N' = |G|^2`. This covers `q = 2`, every `q = 3 mod 4`, and every split `p` at which all points have
   the same level. Then `log q <= gamma_([M]) <= b_C(M)` by Proposition 2.2 with `B = [M]` (and `C_([M]) <= C`).
2. `p | N'`. Let `U = {x : a_xp >= 1}` and `L = {x : a_xp <= e_p - 1}`. Every point lies in `U` or `L`, so one of
   them, call it `B`, has `s = |B| >= M/2 >= 2`. Each `u_x` with `x in U` is divisible by `pi_p`, each with
   `x in L` by `conj(pi_p)`, so the gcd of `{u_x : x in B}` has norm divisible by `p`. Proposition 2.2 applied to
   the normalised configuration (norm `N'`, chord constants `C'_B = C_B |G|^(-1/2) <= C`) gives
   `log p <= b'_C(s) <= b_C(s)`, since `W' <= W`.
Finally apply (3.0). QED

In case 1 with `q = 3 mod 4` one even gets `2 log q <= b_C(M)`, since `q^(2k) || N` forces `q^k | G`.

**Corollary 3.2.** If `M >= 4` and `P_+(N) > C_+^4`, then

```text
M < 2 log N / (log P_+(N) - 4 log^+ C),    and    M (log P_+(N) + 3 log(M/2) - K') <= 2 log N.      (3.2)
```

In particular: if `N` has a prime factor `>= N^lambda C_+^4`, then `M < max(4, 2/lambda)`; a circle whose norm has a
prime factor `> N^(1/2)` (e.g. `N` prime, `C <= 1`) carries at most 3 points on such an arc.

**Theorem 3.3 (prime powers and prime counts).** Let `M >= 4`, normalise, and write `beta' = max_(ceil(M/2) <=
s <= M) b'_C(s) < 2W'/M + 4 log^+ C`.

1. Every cut partition class (all cuts with the same partition) has total weight `<= beta'`. Hence the number of
   distinct partitions and the number of cuts `r` are both `>= W'/beta' > W'/(2W'/M + 4 log^+ C)`.
2. Every prime power `p^e || N'` satisfies `e log p <= 2 beta' < 4W'/M + 8 log^+ C`.
3. `omega(N') >= W'/(2 beta') > W'/(4W'/M + 8 log^+ C)`. For `C <= 1`: `omega(N') > M/4`, and there are more than
   `M/2` distinct partitions and more than `M/2` cuts.
4. Any two rows are separated by more than `(W'/2 + log 2 - 2 log^+ C)/beta'` cuts (more than `M/4` if `C <= 1`)
   and differ at more than `(W'/2 + log 2 - 2 log^+ C)/(2 beta')` primes (more than `M/8` if `C <= 1`).

*Proof.*
1. A class with sides `L, U` (`|L| + |U| = M`) has every one of its cuts constant on both sides. Let `B` be the
   larger side (`|B| >= M/2`); the class weight is at most the sum of the weights of the cuts constant on `B`,
   which is `gamma'_B` (Section 2, after Lemma 2.1), and `gamma'_B <= beta'` by Proposition 2.2.
2. Let `t* = min{t : #{x : a_xp <= t} >= M/2}` (it exists: `t = e` gives `M`). The set
   `L* = {a_xp <= t*}` has `>= M/2` points and `conj(pi)^(e - t*)` divides each `u_x`, `x in L*`. If `t* >= 1`, the
   set `U* = {a_xp >= t*}` is the complement of `{a_xp <= t* - 1}`, which has `< M/2` points, so `|U*| > M/2`, and
   `pi^(t*)` divides each `u_x`, `x in U*`. Hence `(e - t*) log p <= gamma'_(L*) <= beta'` and
   `t* log p <= gamma'_(U*) <= beta'` (the second is trivial if `t* = 0`).
3. Sum item 2 over `p | N'`: `W' = sum e_p log p <= omega(N') 2 beta'`. Strictness comes from (3.0). For `C <= 1`,
   `beta' < 2W'/M`.
4. By the exact margin identity (Proposition 4.1(a) below), `d_xy = W'/2 + m_xy >= W'/2 + log 2 - 2 log C'`, and
   `C' <= C`. Each separating cut contributes at most `beta'` (item 1) and each prime at most `2 beta'` (item 2).
   QED

**Proposition 3.4 (balanced primes: `omega >= (M/2)(1 - o(1))` for `C <= 1`).** If `C <= 1` and `M >= 8`, then
`omega(N') >= (M/2)(1 - M^(-1/3))^2`.

*Proof.* By (T1) for the normalised configuration and `|beta| >= sqrt 2`, `Q < M W'/4` (for `C <= 1`). Put
`h = M^(2/3)/2` and call a cut *balanced* if `|S_(p tau) - M/2| <= h`. The unbalanced cuts weigh
`<= Q/h^2 < M W'/(4h^2) = W' M^(-1/3)`. For one prime `p`, the balanced cuts form an interval
`tau_1 <= tau <= tau_2` (the `S_(p tau)` increase with `tau`). Every `x` in `L = {a_xp < tau_1}` has
`conj(pi)^(e - a_xp)` with `e - a_xp >= tau_2 - tau_1 + 1`, so `gamma'_L >= (tau_2 - tau_1 + 1) log p`, while
`|L| = S_(p tau_1) >= M/2 - h >= 2` and so `gamma'_L < W'/(M/2 - h)`. Summing over primes,
`W' < W' M^(-1/3) + omega(N') W'/(M/2 - h)`, i.e. `omega(N') > (M/2)(1 - M^(-1/3))^2`. QED

(For general `C` the same argument gives `omega(N') > (M/2 - h)(1 - (MW'/4 + 2K log^+ C)/(h^2 W'))/(1 + 4 log^+ C (M/2 - h)/W')`.)

**Proposition 3.5 (where the weight sits).** Let `W'_(<=Y) = sum_(p <= Y) e_p log p`. For every `Y >= 2`,

```text
W'_(<=Y) < pi(Y) (4W'/M + 8 log^+ C),       r < W'/log Y + pi(Y) (4W'/M + 8 log^+ C)/log 5,       (3.3)
```

with `pi(Y) < 1.25506 Y/log Y` (Rosser-Schoenfeld, Illinois J. Math. 6 (1962), Cor. 1 to Thm 2). Taking
`Y = M^(1-delta)`: the primes below `M^(1-delta)` carry less than a fraction
`5.03 M^(-delta) (1 + 2M log^+ C/W')/((1-delta) log M)` of `W'`. Taking `Y = M/(log M)^2`:
`r <= (1 + O(loglog M/log M)) W'/log M + O(M log^+ C/log M)`.

*Proof.* Item 2 of Theorem 3.3 bounds each `e_p log p`, and `e_p <= (4W'/M + 8 log^+ C)/log 5`. Primes `> Y`
contribute at most `W'/log Y` cuts. QED

**Corollary 3.6 (near-extremal clusters).** Fix `eps in (0,1]` and `C > 0`. Under `H(eps)`, with `M >= 4`:

1. **Every prime factor of `N`** is `< C_+^4 (log R)^(4/eps)`; more precisely
   `log P_+(N) <= (4/eps) loglog R - 3 log(M/2) + K' <= (4/eps - 3) loglog R + 3 logloglog R + O_(eps,C)(1)`,
   i.e. `P_+(N) <= (log R)^(4/eps - 3 + o(1)) = M^(4/eps - 3 + o(1))`.
2. Every prime power `p^e || N'` is `< C_+^8 (log R)^(8/eps)`.
3. `M/4 (1 + o(1))^(-1) < omega(N') <= r <= (2/eps + o(1)) M`; for `C <= 1`, `omega(N') > (M/2)(1 - M^(-1/3))^2`.
   So `omega(N')`, `r` and the number of cut classes are all `Theta(M)`.
4. For each fixed `delta in (0,1)`, all but a fraction `O_(eps, C)(M^(-delta)/log M)` of `W'` comes from primes
   in `[M^(1-delta), C_+^4 (log R)^(4/eps)]`.

*Proof.* (1.2) gives `2W/M <= (4/eps) loglog R`; insert in Theorem 3.1. `H(eps)` gives
`3 log(M/2) >= 3 loglog R - 3 logloglog R - O(1)`. Item 2: Theorem 3.3(2) with (1.2). Item 3: Theorem 3.3(3),
`W' >= 3M log M - K'M` (T3) for the lower bound, Proposition 3.5 and `loglog R = (1+o(1)) log M` for the upper.
Item 4: Proposition 3.5. QED

*Exact checks.* `check_subcluster.py`: Theorem 3.1 (both forms) and Theorem 3.3(1)-(3) on every window of
`M >= 4` consecutive lattice points with chord constant `<= 6` on every circle `N <= 3*10^4` (12,148 windows).
`check_known.py`: Theorem 3.1 and Theorem 3.3(2),(3) on all windows with `C <= 9` of the record clusters of
`growth/families.md`. `check_scan_fast.py` (double precision with a `10^-9` guard; no case comes near it): the
refined form of Theorem 3.1 and Theorem 3.3(3) on every window with chord constant `<= 4` on every circle
`N <= 10^6` (9,436 windows). No violation. The tightest smoothness ratio
`log P_+(N)/max_s b_C(s)` is `0.96` (`N = 654385`, `M = 4`, `C = 3.94`); for the first 4-point cluster with
`C <= 1` (`N = 412333805 = 5 * 13^2 * 197 * 2477`, `C = 0.881`) it is `0.871`. So (3.1) is close to sharp for
small `M`.

---------------------------------------------------------------------------------------------

## 4. Profile structure: Plotkin balance, margins, reduced chords

**Proposition 4.1.** Let `z_1, ..., z_M` (`M >= 2`) satisfy (1.1); normalise.

(a) *Exact margin identity.* `m_xy = 2 log|beta_xy| - 2 log C_xy`. Hence
`m_xy >= log 2 + 2 sum_(q odd) v_q(beta_xy) log q - 2 log C'` (the collision-strengthened pair bound).

(b) *Gram form.* Put `f_x = (sqrt(w) s_(cut)(x))_(cuts) in R^r`, `s(x) = +1` if `a_xp >= tau`, else `-1`. Then
`|f_x|^2 = W'`, `<f_x, f_y> = -2 m_xy` for `x != y`, `|sum_x f_x|^2 = 4Q`, and `sum_(x<y) m_xy = M W'/4 - Q`.

(c) *Plotkin balance.* `Q <= M W'/4 + 2K log C' - K log 2 - 2 sum_(l = 3 (4)) E(M, l+1) log l`, hence

```text
Q < M W'/4 - (M^2/2) log M + 3.9626 M^2 + 2K log^+ C .                                   (4.1)
```

(d) *Average margin.* `mbar = (1/K) sum_(x<y) m_xy` satisfies

```text
(M/(M-1)) (log M - 7.9252) + log 2 - 2 log^+ C  <  mbar  <=  W'/(2(M-1)).                  (4.2)
```

(e) *Reduced chords.* `sum_(x<y) log|beta_xy|^2 <= M W'/4 + 2K log C' - Q`, and
`sum_(x<y) log|beta_xy|^2 > (M^2/2) log M - 3.9626 M^2 + K log 2`. So the geometric mean `|beta|^2_geo` of the
`|beta_xy|^2` satisfies `e^(-7.93 M/(M-1)) 2 M^(M/(M-1)) < |beta|^2_geo <= C^2 e^(W'/(2(M-1)))`, and for every
`lambda > 0` at most `K (W'/(2(M-1)) + 2 log^+ C)/(lambda log M)` pairs have `|beta_xy|^2 > M^lambda`.

*Proof.*
(a) `|z_x - z_y|^2 = |g_xy|^2 |beta_xy|^2` and `log|g_xy/G|^2 = W' - d_xy`, so
`2 log|u_x - u_y| = W' - d_xy + 2 log|beta_xy|`; subtract `W'/2 = 2 log sqrt R'`. The lower bound is
`two_thirds` Lemma 3.1 with Lemma 3.3(R), and `C_xy <= C'`.
(b) `<f_x, f_y> = sum_cuts w s(x)s(y) = W' - 2 d_xy` since a cut contributes `-w` exactly when it separates `x, y`.
`|sum f_x|^2 = sum_cuts w (2S - M)^2 = 4Q`. The last identity is `two_thirds` Lemma 2.2 summed with weights.
(c) By (T1), `Q/2 <= M W'/8 + K log C' - sum log|beta|`, and `sum log|beta| >= (K/2) log 2 + sum_l E(M, l+1) log l`
(depth-one inert collisions, `two_thirds` (3.6)). This is the first claim. Next, `2E(M, l+1) >= M^2/(l+1) - M`
(Lemma 3.4 there; the right side is `<= 0` when `l + 1 >= M`), so with `A(x) = sum_(l <= x) log l/(l+1)` and
`theta_3(x) = sum_(l <= x) log l` (both over `l = 3 mod 4`), `2 sum_l E(M, l+1) log l >= M^2 A(M) - M theta_3(M)`.
By Lemmas 7.1-7.4 of `two_thirds`: `A(M) >= S_3(M) - a_3 >= (log M - beta_0 - xi)/2 - a_3` with
`beta_0 = 1 + sigma + (log 2)/2 <= 2.10195`, `xi < 2.70189`, `a_3 <= 0.17438`, and `theta_3(M) < M log 4`. Since
`(beta_0 + xi)/2 + a_3 + log 4 < 3.9626`,

```text
2 sum_l E(M, l+1) log l  >  (M^2/2) log M - 3.9626 M^2 .                                  (4.3)
```

Insert in the first claim and drop `-K log 2`; this is (4.1).
(d) Upper: `Q >= 0` in (b). Lower: by (a), `sum m_xy >= 2 sum log|beta| - 2K log C'`
`>= K log 2 + 2 sum_l E log l - 2K log^+ C`; use (4.3) and divide by `K`.
(e) The upper bound is (T1) again; the lower bound is `2 sum log|beta| >= K log 2 + 2 sum_l E log l` with (4.3).
The count is Markov's inequality for the nonnegative numbers `log|beta_xy|^2`. QED

**Reading (4.1)-(4.2) for near-extremal clusters.** Under `H(eps)`, `mbar = Theta(log M)`: precisely
`log M - O(1) < mbar <= (1/eps + o(1)) log M`, and the Plotkin bound `sum d_xy <= M^2 W'/4` is attained up to a
relative error `O(1/M)`. The cut slack is `Q <= (A/4 - 1/2 + o(1)) M^2 log M`, i.e. the weighted mean square
imbalance `sum w (2S - M)^2/W' = 4Q/W'` is at most `(1 - 2/A + o(1)) M`. The reduced chords have
`|beta|^2_geo in [M^(1 - o(1)), C^2 M^(1/eps + o(1))]`.

*What (4.2) does not say.* Individual margins are not controlled: `m_xy` ranges over `[log 2 - 2 log C', W'/2]`
(the upper end is attained by complementary rows, `u_y = unit * conj(u_x)`), and only `O(K/lambda)` pairs can
have `m_xy > lambda log M`. In all known fakes the margins are uniform (`check_fake_structure.py`:
`48.4 <= m_xy <= 80.6` around `W/(2(M-1)) = 73.3` for the `M = 12` instance); nothing forces this for actual circles.

---------------------------------------------------------------------------------------------

## 5. Collision stability, and the sharp regime near the 2/3 bound

Normalise. For an odd prime power `q^a` let `kappa_(q^a)` be the number of pairs with `q^a | beta_xy`, counting
for `p | N'` only pairs on a common `p`-level (pairs on different levels have `v_p(beta) = 0`, `two_thirds`
Lemma 3.3(S|)). Let `E_(q^a)` be its pigeonhole minimum: `E(M, (q - chi(q)) q^(a-1))` for `q` not dividing `N'`,
and `sum_t E(n_(pt), (p-1) p^(a-1))` for `p | N'`. Define three nonnegative *excesses*:

```text
X_chord = sum_(x<y) [ log|beta_xy| - (1/2) log 2 - sum_(q odd) v_q(beta_xy) log q ]      (cofactors of beta)
X_coll  = sum_(q odd) log q sum_(a>=1) (kappa_(q^a) - E_(q^a))                           (collisions above pigeonhole)
X_prof  = [left side of two_thirds (5.1)] - F*(M)                                   (levels and cuts above optimum)
```

`X_chord >= 0` is `two_thirds` Lemma 3.1, `X_coll >= 0` is pigeonhole (Lemmas 3.2-3.4 there), `X_prof >= 0` is
step 1 of the proof of Theorem 5.2 there.

**Proposition 5.1 (stability of the collision structure).** `sum log|beta| + Q/2 = F*(M) + X_chord + X_coll + X_prof`,
and each excess is at most

```text
Xi := M W'/8 + K log C' - F*(M)  <=  (M^2/8) (W'/M - 3 log M + K')  =  (M^2 log M/8) (A - 3 + K'/log M).     (5.3)
```

*Proof.* `sum_(x<y) v_q(beta_xy) = sum_a kappa_(q^a)`, so
`sum log|beta| = (K/2) log 2 + sum_q log q sum_a E_(q^a) + X_coll + X_chord`, and
`(K/2) log 2 + sum_q log q sum_a E_(q^a) + Q/2` is exactly the left side of `two_thirds` (5.1). (T1) bounds the
total by `M W'/8 + K log C'`, and (T2) bounds `F*(M)` below. QED

So in aggregate the residue collisions exceed the pigeonhole minimum by at most `(A - 3)/3 + O(1/log M)` times the
main term `F*(M) ~ (3/8) M^2 log M`, the reduced chords carry at most that much beyond their forced collision
part, and the level profiles are that close to the chain-lemma optimum. For bounded `A` this is a
constant-factor statement; it is a genuine stability statement only near `A = 3`:

**Theorem 5.2 (sharp regime: almost every small split prime divides `N'`).**

```text
sum_(p <= M, p = 1 (4), p does not divide N') (log p)/(p - 1)  <=  (W' - 3 M log M + K' M)/(2M)  =  ((A - 3) log M + K')/2.
```

*Proof.* Compare the left side of `two_thirds` (5.1) with the explicit lower bound `D(M)` of its Section 6 term
by term: every term of `two_thirds` (5.1) is at least the corresponding term of `D(M)`, and for split `p <= M` with `p` not
dividing `N'` the difference is at least
`[E(M, p-1) - (c_p M^2/4 - M/2)] log p >= (M^2/(2(p-1)) - c_p M^2/4) log p >= M^2 log p/(4(p-1))`, because
`E(M, p-1) >= M^2/(2(p-1)) - M/2` and `c_p <= 1/(p-1)`. Hence
`(M^2/4) sum_(p not dividing N') log p/(p-1) <= M W'/8 + K log C' - D(M)`, and `8D(M)/M > 3M log M - K_0 M` (T2). QED

By Mertens' theorem in the progression `1 mod 4` (e.g. K. S. Williams, J. Number Theory 6 (1974) 353-359),
`sum_(p <= M, p = 1 (4)) log p/(p-1) = (1/2) log M + O(1)`. So if `A <= 3 + delta`, i.e. `M` is within a factor
`1 + O(delta)` of the two-thirds bound, then all but a fraction `delta + O(1/log M)` (in the measure `log p/p`) of
the split primes below `M` divide `N'`; and by (5.3) the inert collisions, the reduced-chord cofactors and the
level profiles are within a relative `delta/3 + O(1/log M)` of their minima. For `A >= 4 + o(1)` Theorem 5.2 is
empty. Every known fake has `A >= 8 - o(1)` (`round5/README.md`; the draft `exact.md` gives 8, the refereed
`construct.md` 27 in `A_1`), so this regime concerns the constant `2/3`, not the
growth rate, and is recorded only as the stability statement the task asked for.

---------------------------------------------------------------------------------------------

## 6. A genuine (non-local) consequence: no character lives on few primes

Smoothness makes Baker's method usable on characters with a few primes. This is the only place where the actual
angles `arg pi_p` enter.

**Quoted input (Matveev, Izv. Math. 64 (2000) 1217-1269, Cor. 2.3).** Let `alpha_1, ..., alpha_n` be nonzero
elements of a number field of degree `D`, `b_j` integers, `B >= max(1, |b_j|)`, `A_j >= max(D h(alpha_j),
|Log alpha_j|, 0.16)`, and `Lambda = alpha_1^(b_1) ... alpha_n^(b_n) - 1 != 0`. Then
`log|Lambda| > -c(n) D^2 (1 + log D)(1 + log B) A_1 ... A_n`, with an effective `c(n) <= c_0^n` for an absolute
`c_0`. (Only the exponential dependence on `n` is used; the explicit `c(n)` is quoted from memory.)

**Proposition 6.1.** Fix `eps in (0,1]`, `eta in (0,1)` and `C > 0`. There is an effective `R_1(eps, eta, C)`
such that for `R >= R_1` the following holds. Let `z_1, ..., z_M` be lattice points on an arc of length
`<= C sqrt R` of `x^2 + y^2 = R^2` with `M >= eps log R/loglog R`; normalise. Let `c in Z^M`, `sum c_x = 0`,
`||c||_1 <= exp((log R)^(eta/2))`, and `v = sum_x c_x a_x in Z^(primes of N')`. If `v != 0`, then

```text
|supp v|  >  (1 - eta) loglog R / logloglog R  =  (1 - eta - o(1)) log M / loglog M.
```

*Proof.* Write `u_x = R' e^(i phi_x)`. Since `u_x/conj(u_x) = eps_x^2 prod_p (pi_p/conj pi_p)^(2 a_xp - e_p)` and
`sum c_x = 0`,

```text
Pi := prod_x (u_x/conj u_x)^(c_x) = s prod_(p in supp v) (pi_p/conj pi_p)^(2 v_p),   s = prod_x eps_x^(2 c_x) in {+1, -1}.
```

All `phi_x` lie (mod `2 pi`) in an interval of length `Delta = C R^(-1/2)` (arc length over radius; normalisation
rotates all points equally); centring it, `|sum c_x phi_x| <= ||c||_1 Delta/2`, so
`|Pi - 1| = |e^(2i sum c_x phi_x) - 1| <= ||c||_1 Delta`. If `v != 0`, then `Pi != 1`: otherwise
`prod (pi_p/conj pi_p)^(2v_p) = +-1`, which unique factorisation in `Z[i]` forbids (`pi_p`, `conj pi_p`
non-associate). Apply Matveev in `Q(i)` (`D = 2`) to `alpha_0 = s`, `alpha_p = pi_p/conj pi_p`: `n = |supp v| + 1`,
`h(pi_p/conj pi_p) = (1/2) log p` (minimal polynomial `p X^2 - 2(a^2 - b^2) X + p` for `pi_p = a + bi`), so
`A_p = max(log p, pi)`, `A_0 = pi`. By Corollary 3.6(1), `log p <= (4/eps) loglog R + 4 log^+ C`, hence
`A_p <= (5/eps) loglog R` for `R >= R_1`. By Corollary 3.6(2), `|v_p| <= ||c||_1 e_p <= ||c||_1 (log R)` for large
`R`, so `1 + log B <= 2 (log R)^(eta/2)`. Therefore

```text
(1/2) log R - (log R)^(eta/2) - log C  <=  -log|Pi - 1|  <  8 pi (1 + log 2) c_0^n (log R)^(eta/2) ((5/eps) loglog R)^(n-1).
```

Taking logarithms, with `n - 1 = |supp v|` and `L = loglog R`:
`L - O(1) < (eta/2) L + |supp v| (log L + log(5 c_0/eps)) + O(1)`. If `|supp v| <= (1 - eta) L/log L`, the right
side is `(1 - eta/2) L + O(L/log L)`, a contradiction for `R >= R_1`. The `M`-form uses
`loglog R = (1 + o(1)) log M` (Section 1). QED

**Remarks 6.2.**

1. *Smoothness is what makes this work.* Without Theorem 3.1 the product `A_1 ... A_n` may be as large as `W^n`, and
   Matveev's bound is useless already at `n = 2`. With all primes `<= (log R)^(4/eps)` it reaches supports
   `~ log M/loglog M` and arbitrarily large coefficients `|v_p| <= exp((log R)^(eta/2))`.
2. *What it excludes.* Profiles in which an integer row combination cancels on all but `o(log M/loglog M)`
   primes; e.g. a pure Hadamard core `H` (where `c = H e_a` gives `v = (M/2) e_a`, support 1). Pure cores are
   already uniformly bounded (`growth/walsh-extraction.md` Theorem B), so this is not new information about them.
3. *It does not separate the known fakes.* For the profile of `round5/construct.md` (`b` unflipped copies of
   each column of a normalised Hadamard matrix `H`, plus one flipped column per row, with the identity pairing):
   if `c != 0` is zero-sum and `T = supp(H^T c)`, then `v` is nonzero on the `b|T|` unflipped copies of the
   labels in `T` and on the flipped column of every `x in supp c` whose label is not in `T` (there
   `v = -H_(x a_x) c_x`); the pairing is injective except for one repeated label, so
   `|supp v| >= b|T| + |supp c| - |T| - 1 >= 2 sqrt((b-1) M) - 1` by the uncertainty principle
   `|supp c| |T| >= M` (from `H^T H = M I`: `M ||c||_2^2 = ||H^T c||_2^2 <= |T| ||c||_1^2 <= |T| |supp c| ||c||_2^2`).
   `check_fake_structure.py` finds minimum supports 53 (pair characters) and 20 (`c = ` a Hadamard column) at
   `M = 12`, `b = 8`, against the bound `2 sqrt 84 - 1 = 17.3`. For the skew double core of `round5/exact.md` (draft) every nonzero
   character is nonzero on at least `M/2` types (round-8 drafts `lw.md` Theorem 5.1, `weakest.md` Lemma 3.2).
   Both are far above `log M/loglog M`. More generally Proposition 6.1 is implied by any single-form
   Lang-Waldschmidt inequality, and `lw.md` Theorem 6.1 (draft) builds fakes satisfying those.

---------------------------------------------------------------------------------------------

## 7. Consistency: the whole local structure holds for fake circles

**The class.** `A_1(Q_0)` is the class of angle models with residue data of `round5/construct.md` Section 1
(refereed): genuine primes `p_j = 1 mod 4` with exponents `e_j` and a normalised profile (`min_x a_xj = 0`,
`max_x a_xj = e_j`), real angles, points `R i^(-k_x) exp(i <2a_x - e, phi>)` with angular positions in an
interval of length `Delta = C e^(-W/4)`, and residue data `rho_x^(m) in (Z[i]/m)^*` on the conic
`rho conj(rho) = N` for every odd prime power `m = q^a <= Q_0` with `q` not dividing `N`, satisfying the pair case
of axiom (R). Write `A_1^>(Q_0)` for the models in which every `p_j > Q_0`. The refereed models of
`construct.md` Theorem 2 are of this kind (their primes exceed `160 M^2 L/C > Q_0`).

**Proposition 7.1 (transfer).** Let `Q_0 >= 5M/4`. For every model in `A_1^>(Q_0)`, with `G = 1`, `W' = W`,
`C' = C`, and `beta_xy` *defined* by `|beta_xy| = |z_x - z_y| e^(-(W - d_xy)/2)`, the following hold verbatim:
Lemma 2.1, Proposition 2.2 (even with `F*` replaced by the larger `F_free(s)`, which uses the `p`-free branch at
every split `p`), Theorem 3.1 for the primes `p_j`, Theorem 3.3, Propositions 3.4 and 3.5, Corollary 3.6 (1)-(4)
under `H(eps)`, Proposition 4.1 and Proposition 5.1. Theorem 5.2 holds and reads `A >= 4 - o(1)`, since no
`p <= M` divides `N`.

*Proof.* The proofs above use exactly four inputs, each available in the class.
1. *Algebra.* (2.1) and Proposition 4.1(a),(b) are identities in the profile, the weights and the numbers
   `|z_x - z_y|`, `|beta_xy|` (cut identity, `two_thirds` Lemma 2.2); they hold by the definition of `|beta|`.
2. *Chord bound.* `|z_x - z_y| = 2R |sin((delta_x - delta_y)/2)| <= R Delta = C sqrt R`, which is (1.1).
3. *Pair lower bound.* `|beta_xy| = 2 e^(d_xy/2) |sin((phi_x - phi_y)/2)|`, so the pair case of (R) is
   `|beta_xy| >= sqrt 2 D_xy`, with `D_xy` the product of the collision moduli `q^a <= Q_0`; this replaces
   `two_thirds` Lemmas 3.1 and 3.3.
4. *Pigeonhole.* The residues of a fixed modulus lie in the conic fibre of size `(q - chi(q)) q^(a-1)`
   (`two_thirds` Lemma 3.2, valid since `q` does not divide `N`), so any `s` points have at least
   `E(s, (q - chi(q)) q^(a-1))` colliding pairs. The moduli used by `F*(s)`, `s <= M`, satisfy
   `(q - chi(q)) q^(a-1) < s`, hence `q^a < 5s/4 <= Q_0`.
Inputs 3-4 give `sum_B log|beta| + Q^B/2 >= F_free(s) >= F*(s)` for every sub-model (the cuts constant on `B`
removed, `W_B = W - gamma_B`, chord constant `C_B e^(-gamma_B/4)`), which is (T1) for sub-models. Everything else
in Sections 2-5 is derived from (T1), (T2) and inputs 1-2. For Theorem 5.2 the sum runs over all split
`p <= M`, which is `(1/2) log M + O(1)`. QED

**Corollary 7.2 (the structure is consistent at `W = O(M log M)`).** For every `C in (0, 1]` and every
`M >= M_0`, there are models in `A_1^>(2M)` with `W <= (27 + o(1)) M log M + 9 M log(1/C)`
(`round5/construct.md` Theorem 2 and Corollary 6.1, refereed in `round5/referee_construct.md`; the primes exceed
`Q_0` by construction). They satisfy `H(eps)` for every `eps < 2/27`, and by Proposition 7.1 every conclusion of
Sections 2-5 (in particular `N` is `M^(O(1))`-smooth, `omega = Theta(M)`, the average margin is `Theta(log M)`, the
reduced chords are `M^(O(1))` on average, the collisions are within a bounded factor of pigeonhole). So **no
argument whose inputs are the axioms of `A_1(Q_0)` plus the conclusions of Sections 2-5 can prove
`M = o(log R/loglog R)`**: Sections 2-5 add no information to the local class.

**Exact illustration** (`check_fake_structure.py`, the `M = 12` instance of `construct.md`, `C = 1`, primes
`10000121..10003529`):

| quantity | fake | theorem |
| --- | --- | --- |
| `W`, `A = W/(M log M)` | 1611.83, 54.05 | `A >= 3 - K'/log M` |
| `max log p` | 16.12 | `<= 2W/M - 3 log(M/2) + K_0 = 277.2` (Thm 3.1) |
| sub-cluster bound, all 4083 subsets | max `gamma_B - b(s) = -47.7` | `<= 0` (Prop 2.2) |
| `omega = r` | 100 | `> M/4 = 3` (Thm 3.3) |
| distinct cut partitions | 23 | `> M/2 = 6` |
| cut slack `Q` | 193.4 | `< MW/4 = 4835.5` (4.1) |
| margins | min 48.4, mean 70.3, max 80.6 | mean `<= W/(2(M-1)) = 73.3` (4.2) |
| min Hamming distance of rows | 53 | `> M/4` (Thm 3.3(4)) |
| min character support (pairs; Hadamard columns) | 53; 20 | `>= 2 sqrt(7M) - 1 = 17.3` (Remark 6.2(3)) |

**What separates actual circles from these models.** In a model, the "Gaussian primes" `sqrt(p) e^(i phi_p)` and
hence the `beta_xy` are complex numbers. In an actual cluster `beta_xy = u - eps conj(u) in Z[i]` with
`u = z_x/g_xy in Z[i]`; conversely a model whose angles are the arguments of actual Gaussian primes of norm `p_j`
*is* an actual configuration (its points are the Gaussian integers `eps_x prod pi^(a) conj(pi)^(e - a)`), and its
`beta_xy` are then automatically integral. So, given the local structure, the integrality of the reduced chords is
equivalent to the genuineness of the prime angles; the only inputs left for a growth proof are properties of the
specific numbers `arg pi`, `N(pi) <= C_+^4 (log R)^(4/eps)`, beyond the Liouville gaps, Proposition 6.1, and the
residue data at moduli `<= Q_0` (and, per `referee_lw.md`, data the models do not represent at all, such as
residues at all moduli with the prime-level parity law, or additive relations among the points).

---------------------------------------------------------------------------------------------

## 8. Exploiting the structure: what each suggested route gives

Each route below was tried with the structure of Sections 3-6 in hand. None gives growth; for each, the exact
reason is stated.

### 8.1 Three-term identities and S-unit equations

**Lemma 8.1 (triangle identity).** For distinct `x, y, w` of a normalised configuration let
`G_3 = gcd(u_x, u_y, u_w)` and `A_w = g_xy/G_3`, `A_x = g_yw/G_3`, `A_y = g_xw/G_3` (gcds of the `u`'s). Then
`A_w, A_x, A_y` are pairwise coprime products of Gaussian primes above primes of `N'`, and

```text
A_w beta_xy + A_x beta_yw = A_y beta_xw,                                                    (8.1)
W'/6 - m_xy - (8/3) log^+ C + (4/3) log 2  <=  log|A_w|^2  =  W'/2 - m_xy - gamma_3  <=  W'/2 - log 2 + 2 log^+ C.
```

*Proof.* At a Gaussian prime `varpi` let the three valuations be `t_1 <= t_2 <= t_3` (in some order of the
points). The pairwise gcds have valuations `t_1, t_1, t_2`, and `G_3` has `t_1`, so only the gcd of the pair
avoiding the minimiser can exceed `G_3` at `varpi`. (8.1) is `(u_x - u_y) + (u_y - u_w) = u_x - u_w` divided by
`G_3`. The norm identity is `log|g_xy|^2 = W' - d_xy`; the lower bound is Proposition 2.2 for `s = 3`
(`F*(3) = (3/2) log 2`), the upper bound is `m_xy >= log 2 - 2 log^+ C` (Proposition 4.1(a)). QED

So every triangle with typical margins gives a three-term equation whose S-unit parts have norm
`>= N'^(1/6 - o(1))` while (Proposition 4.1(e)) the coefficients `beta` are `M^(O(1))` on average.

**Proposition 8.2 (counting S-unit solutions cannot bound `M`).**
1. Dividing (8.1) by `A_y beta_xw` gives `X + Y = 1` with `X, Y` in the group generated by `i`, the Gaussian
   primes over `N'` and the `beta`'s, of rank `>= omega(N') > M/4` (Theorem 3.3, `C <= 1`). The cluster supplies
   at most `binom(M, 3)` such solutions. The general upper bounds for the number of solutions in a group of rank
   `r` are exponential in `r` (Beukers-Schlickewei, Acta Arith. 78 (1996) 189-199: at most `2^(8(r+1))`), so they
   are compatible with `binom(M, 3)` solutions; and growth faster than any polynomial in the number of primes is
   unavoidable for any bound depending only on `|S|`: Erdos-Stewart-Tijdeman (Compositio Math. 66 (1988)) give sets of `s` primes with at least
   `exp((4 - o(1)) (s/log s)^(1/2))` coprime S-unit solutions of `a + b = c`, which is `> M^3` for `s >= M/4`
   large.
2. (8.1) holds verbatim, with complex `beta`, in every model of Section 7 (it is the triangle identity of the
   points divided by fake gcds). It therefore carries no information beyond the integrality of the `beta`'s.

So finiteness and counting theorems for S-unit equations, whose bounds depend on `|S| = Theta(M)`, cannot see
the cluster.

**Proposition 8.3 (S-unit relations with small `S` are parallel chords).** Under the hypotheses of Proposition
6.1, let `x, y, w, v` be distinct. Either `a_x + a_y = a_w + a_v`, and then `u_x u_y = unit * u_w u_v` and
`phi_x + phi_y = phi_w + phi_v` exactly (the chords `xy` and `wv` are parallel); or the vector
`a_x + a_y - a_w - a_v` is nonzero at more than `(1 - eta) loglog R/logloglog R` primes.

*Proof.* In the first case the two products have the same factorisation and norm, so their ratio is a unit; its
argument `phi_x + phi_y - phi_w - phi_v` lies in `(pi/2) Z` and has absolute value `<= 2 Delta < pi/2`, so it is
`0`. Otherwise apply Proposition 6.1 to `c = e_x + e_y - e_w - e_v`. QED

The exact relations are satisfied by the models too (there the same angle sum is a multiple of `pi/2` and small).
Hence multiplicative relations among cluster points involving few primes are either exact and harmless, or do not
exist.

### 8.2 abc over `Q(i)`

Take a triangle with `N'` squarefree and suppose, for simplicity, that the three terms of (8.1) are coprime (common
factors only lower the height). The radical then contains every prime of `A_w A_x A_y` once, so
`Norm(rad) >= |A_w A_x A_y|^2 = e^(3W'/2 - (m_xy + m_yw + m_xw) - 3 gamma_3)`, while the largest term has
`|A_w beta_xy|^2 = |u_x - u_y|^2 e^(-gamma_3) <= C^2 e^(W'/2 - gamma_3)`. With `gamma_3 <= W'/3 + O(1)`
(Proposition 2.2), the radical exceeds the largest term by a factor `>= e^(W'/3 - (m_xy + m_yw + m_xw) - O(1))`,
so abc holds with room for every triangle with `o(W')` margins. This is the triangle version of `weakest.md`
Propositions 4.1-4.2 (round-8 draft: abc is vacuous on squarefree profiles; for general profiles it gives only a
restricted theorem when `rad(N) <= N^(1/4 - eta)`). In the near-extremal regime the radical cannot be that small
unless `A` is large: by Theorem 3.3(2) every prime `p > M^(1-delta)` has
`e_p < (4W'/M + 8 log^+ C)/((1 - delta) log M) = (4A + o(1))/(1 - delta)`, and by Proposition 3.5 the smaller
primes carry `o(W')`, so `rad(N') >= N'^((1 - delta)/(4A) - o(1))`.

### 8.3 Smooth numbers near squares

For each pair, `w_xy = u_x conj(u_y)` is a `P_+(N')`-smooth Gaussian integer of norm `N'^2` with
`|arg w_xy| <= Delta`; equivalently `X = Re w_xy` satisfies `0 < N' - X = |u_x - u_y|^2/2 <= C^2 N'^(1/2)/2` and
`(N' - X)(N' + X) = (Im w_xy)^2`. Writing `N' - X = d a^2`, `N' + X = d b^2` gives `2N' = d(a^2 + b^2)` with
`d a^2 <= C^2 N'^(1/2)/2`. So the cluster is the same thing as `M - 1` smooth Gaussian integers of one fixed norm
`N'^2` inside one sector of width `Delta`; a smooth-numbers-in-short-intervals theorem would have to count
elements of a single norm fibre in a sector, which is the original problem rotated by `conj(u_1)`. Theorems about
`Psi(x + y, P) - Psi(x, P)` count integers in intervals and give lower bounds; no upper bound of the needed kind
follows. No reduction was found.

### 8.4 Counting and equidistribution (heuristic; evidence only)

A counting proof would show that the expected number of clusters among all choices of `r = Theta(M)` primes below
`P = (log R)^(4/eps)` is `< 1` *and* control the error. The main term is tiny (the angle vector must lie within
`Delta` of a subtorus of codimension equal to the rank of the profile). For the error, consider even one pair
form `psi = sum_(p in D) +-theta_p` (`|D| = n`, weight `d = n log P ~ W/2`): Fourier expansion at precision `Delta`
needs the Hecke sums `S_k = sum_(N pi <= P) e^(4 i k arg pi)` for `|k| <= 1/Delta`. Under GRH for Hecke
L-functions of `Q(i)`, `|S_k| << P^(1/2) log(|k| P)`. The smoothed indicator of `[-Delta, Delta]` has
`~ 1/Delta` Fourier coefficients of size `~ Delta`, so the error is `<< (P^(1/2) W)^n` against the main term
`Delta pi(P)^n`; the ratio is `~ (W log P)^n e^(-d/2)/Delta`, which is `>= 1` whenever `Delta <= e^(-d/2)`, in
particular for an actual pair (`Delta = C e^(-W/4)`, `d ~ W/2`). So square-root cancellation reaches the
Liouville scale and no further, even for a single form; simultaneous forms are out of reach. (This is the counting counterpart of `round4/hardness.md` Theorem 3: Liouville is optimal for
pair forms.)

### 8.5 Recursion on sub-clusters and archimedean (additive) data

Dividing `s` points by their gcd shrinks the chord constant by `|G_B|^(-1/2)` (Remark 2.3(3)); applying the 2/3
inequality to every such sub-cluster is exactly Proposition 2.2, and iterating along a chain of halvings only
re-proves Theorem 3.1. Additive data (integer combinations `sum c_x z_x` small by Taylor expansion at the
archimedean place) is the single-place auxiliary-polynomial method of `threshold/`, barred by
`tau*(M) < 2` there; smoothness of `N` does not enter it.

### 8.6 What a growth proof must use

Combining Corollary 3.6, Proposition 6.1 and Corollary 7.2: under `H(eps)`, an actual cluster is a configuration
on a `C_+^4 (log R)^(4/eps)`-smooth norm with `Theta(M)` primes, near-Plotkin profile and near-pigeonhole
collisions, all of which fake circles share. A proof of `M = o(log R/loglog R)` must therefore use the actual
values `arg pi`, `N(pi) <= C_+^4 (log R)^(4/eps)`, through simultaneous approximation properties stronger than:
every single-form Liouville gap, Matveev on supports `<= (1 - eta) log M/loglog M` (Proposition 6.1), and residues
at moduli `<= Q_0` used through pairs (Section 7). By `round8/lw.md` (refereed in `referee_lw.md`: survives as a
barrier relative to the angle-model inputs), single forms do not suffice even under the Lang-Waldschmidt
conjecture in its max form; by `round8/weakest.md` (draft), five simultaneous forms do. **The structure theorem
lets any such hypothesis be restricted to logarithms of Gaussian primes of norm `<= C_+^4 (log R)^(4/eps)`, with
coefficients `|v_p| <= ||c||_1 e_p`, `e_p log p < (8/eps) loglog R + 8 log^+ C`, and supports of more than
`(1 - eta) log M/loglog M` primes**, which is a strictly smaller (and more natural) range than "all algebraic numbers of the relevant height".

---------------------------------------------------------------------------------------------

## 9. Open points

1. **A local class beyond `A_1`.** Proposition 7.1 shows that Sections 2-5 are consequences of the `A_1` axioms.
   The open local question remains the one of `round5/README.md`: does adding genuine residues at all moduli
   (P3-all, with the prime-level parity law) force `W/(M log M) -> infinity`? Theorem 3.1 adds one fact to that
   question: the primes of `N` lie below `C_+^4 (log R)^(4/eps)`, so every residue datum is that of a product of
   Gaussian primes from a set of size `(log R)^(O(1))`.
2. **Matveev beyond `log M/loglog M`.** Proposition 6.1 is the reach of proved linear-forms bounds. Supports
   `>= M/2` (the skew double core) need simultaneous forms (`weakest.md`) or Vojta (`growth/conditional.md`).
3. **Prime-weight concentration.** All fakes put their primes in one window `[P, P e^(1/(4r))]`. Nothing in
   Sections 2-5 forces this for actual clusters (the weights only need `log p <= 2W'/M + O(1)` and most weight above
   `M^(1-delta)`). Whether near-extremal actual clusters must concentrate their prime weights is open; if they had
   to, the problem would concern angles of Gaussian primes in very short intervals.

## 10. Files (`round8/stability_checks/`)

| script | what it checks | output |
| --- | --- | --- |
| `stab_lib.py` | Gaussian arithmetic, lattice points by factorisation, exact `F*(s)` (from `two_thirds/finite_inequality.py`), `b(s)` | (library) |
| `check_subcluster.py 30000 6 6 12` | exact identity (2.1) on 148,370 subsets; Prop 2.2; Thm 3.1 (both forms); Thm 3.3(1)-(3) on 12,148 windows, `N <= 3*10^4`, `C <= 6` | `check_subcluster_out.txt` |
| `check_known.py` | Prop 2.2 (60,448 subsets), Thm 3.1, Thm 3.3(2),(3) on the record clusters of `growth/families.md` and the 10-point cluster | `check_known_out.txt` |
| `check_scan_fast.py 1000000 4` | refined Thm 3.1 and Thm 3.3(3), every window with `C <= 4`, `N <= 10^6` (double precision) | `check_scan_fast_1e6_C4_out.txt` |
| `check_fake_structure.py` | every local conclusion on the `M = 12` fake of `construct.md` (table of Section 7) | `check_fake_structure_out.txt` |
| `check_constants.py 3000` | (4.3) for `M <= 3000`, the constant `3.9626`, (T2) for `s <= 60`, the threshold `M >= 40` of Prop 3.4 | `check_constants_out.txt` |

All scripts are pure Python 3 with `mpmath` (and `numpy` for the fake profile); each prints
`ALL CHECKS PASSED: True`. The checks are sanity checks of proved statements; no proof depends on them.
