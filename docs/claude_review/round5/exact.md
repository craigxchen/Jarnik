# Round 5 (exact): the angle-model problem as an exact finite optimisation

Target: turn the angle-model problem of `round4/outside.md` into an exact finite optimisation,
compute `W_min`-type data for small `M`, report the trend of `W/(M log M)`, the optimal profiles, and
whether Hadamard-like profiles are optimal.

Context read first: `round4_context.md`, `round4/outside.md` (§§1-2) and its referee report,
`growth/walsh-extraction.md` (Theorem F, Lemma F.1, Prop P), `two_thirds/proof.md`, and the two
concurrent round-5 write-ups `round5/lower.md` (`(48 + o(1)) M log M`) and `round5/construct.md`
(`(27 + o(1)) M log M`, refereed in `round5/referee_construct.md`).

Scripts and outputs are in `round5/exact/`. Every script runs from that directory with
`python3 <script> [args]`; logs are in `exact/logs/`, solver results in `exact/results/`.

Conventions, as in `outside.md` and `construct.md`:
* *Theorem*, *Proposition*, *Lemma* mean proved here with all quantifiers, or cited with the exact
  source.
* *Exact* means integer / `Fraction` arithmetic, or 50-digit `Decimal` logarithms with explicit
  margins.
* *Solver-optimal* means: a feasible point verified exactly, whose value equals the dual bound reported
  by the MILP solver HiGHS (double precision, relative gap `1e-9`) on a valid relaxation. This is
  strong evidence, **not** a proof. It is never used in a proof.
* `W = log N = 2 log R`; an arc of length `<= C sqrt R` has angular width `Delta <= C e^(-W/4)`.

---------------------------------------------------------------------------------------------

## 0. Verdict

**`W_min/(M log M)` stays bounded, and the bound improves.** The exact finite optimisation identifies
the optimal profiles: for the margin structure that fake-circle constructions need, the optimum is a
**skew double Hadamard core**, and it meets the Plotkin bound exactly. This gives:

| quantity | before this note | here |
|---|---|---|
| `limsup W_min/(M log M)`, class `A_1(2M)` (point residues, axiom (R), moduli `<= 2M`) | `27` (construct), `48` (lower) | **`12`** (Theorem 5.1) |
| same, residue moduli `<= M^A` (`A >= 1`) | `9(2 + A)` | **`8 + 4A`** |
| same, no residue data (`A_0`) | `18` (construct) | **`8`** |
| `liminf W_min/(M log M)`, `A_1(2M)` | `3` (`two_thirds`) | `3` (unchanged) |
| best constant `kappa` in "`M <= kappa log R/loglog R`" provable from the class | `[2/27, 2/3]` | **`[1/6, 2/3]`** |

The main results:

1. **Exact reformulation** (§1). In the large-prime regime of `outside.md` Theorem A.1, an angle
   model exists iff, for every rational character `lambda`, one explicit inequality holds
   (Prop 1.1). Which inequality depends on the parity of `lambda.k`, where `k` are the windings.
   For fixed combinatorial data (exponent rows, windings `k`, normalised positions `u = delta/Delta`)
   the admissible log-weight vectors are sandwiched between two polyhedra that differ by
   `O(e^(-W/2))` (Prop 1.2). So "minimal `W` for fixed data" is a linear programme in the column
   log-weights.
2. **The normalised problem** (§2). In the limit where all primes grow (`log p_S = tau omega_S`,
   `tau -> infinity`, `M` fixed), angle models reduce exactly to a finite mixed-integer problem.
   * The variables are the binary column patterns `S` (complement classes) with real weights
     `omega_S in {0} cup [1, infinity)`, and a homomorphism `h : G_T -> Q/Z` that encodes the windings.
   * The constraints are `sigma_lambda(omega) >= r(lambda)` for every rational character in `ker h`.
   * Its value `Omega_r(M)` is finite and exactly verifiable (Lemmas 2.2-2.4).
   * **Theorem 2.5**: for fixed `M` and `C`, the least `W/tau` of angle models whose primes all exceed
     `e^tau` tends, as `tau -> infinity`, to a number in `[Omega_0(M), Omega_+(M)]` (zero margin,
     closed and strict). It equals `Omega_0(M)` for `C >= C_0(M, Q_0)`.
3. **Exact data, `M = 4..10`** (§3, Table 3.1). Zero margin:
   `Omega_0 = 3, 4, 7, 8` and `Omega_+ = 3, 7, 7` for `M = 4..7`. These are exact for `M <= 6`
   (exhaustive exact certificates) and solver-optimal for `M = 7`. So the character constraints
   alone cost only about one weight unit per row.

   Norm-proportional margins `r = rho ||lambda||_1/2` (the margin structure of every union-bound
   construction): **the Plotkin bound `2(M-1) rho` (`M` even) or `2 M rho` (`M` odd) is attained for
   every `M` in `5..9` and every `rho >= 1`**. These values are exact: the Plotkin bound is a theorem
   (Lemma 2.2(e)), the optimal points are verified exactly at `rho = 1`, and homogeneity
   (Lemma 2.2(d)) carries attainment to every `rho >= 1`. `M = 4` is the exception (`10 > 6` at
   `rho = 1`), caused by the aligned 4-flat.
4. **Theorem 4.1 (skew double core).** Let `K` be a skew-Hadamard matrix of order `M >= 8`. The
   `2(M-1)` balanced cuts `S_a = {x : K(x,a) = -1}` and `S'_a = S_a xor {0, a}`, all with weight 1, give:
   * a cyclic group `G_T = Z/(M/2)`, so all rational characters can be excused by windings;
   * `sigma_c >= ||c||_1/2` for **every** integral character.

   Hence `Omega^prop_rho(M) = 2(M-1) rho` exactly for every `rho >= 1` and every skew-Hadamard order
   `M >= 8`.

   **Hadamard-like profiles are optimal, but in double form.** The single-core profiles of
   `construct.md` and `lower.md` (Hadamard core, copies, one flip per row) are not optimal. At
   `M = 8` the best weights on those supports give `30` (skew twins) and `44` (random twins), against
   the optimum `14`, at `rho = 1` (exact LP values on fixed supports, §4.3).
5. **Theorem 5.1 (fake circles at `12 M log M`).** The skew double core plus the residue design of
   `construct.md` Lemma 4.1 gives angle models in `A_1(Q_0)` with
   `W <= 2(M-1)(4 log M + 2 log L + 2 log(1/C) + 12.7) + 1/4`, where `L ~ 2M log M` for `Q_0 = 2M`.
   * Exact instances: at `M = 8`, `W = 409.9` (`24.6 M log M`; refined union bound `356.6`), against
     `878.3` for the `b = 8` Paley profile (`referee_construct.md` R6). At `M = 12`, `W = 698.8`
     (`23.4 M log M`), against `1611.8`.
   * **Corollary 5.3:** no argument whose inputs lie in `A_1(2M)` proves
     `M <= (1/6 - eps) log R/loglog R`.
6. **Where the remaining gap `[3, 12]` lives** (§6). Every normalised profile with norm-proportional
   margins obeys the Plotkin bound (Lemma 2.2(e), Prop 6.1). Translated through the counting of the
   union bound (a heuristic, made rigorous here only for Theorem 5.1), every construction that
   handles characters by a union bound over random positions has constant `>= 8` (`A_0`) or `>= 12`
   (`A_1`). So Theorem 5.1 is optimal for that route. A smaller constant needs structured positions or
   a different residue design. A larger lower bound needs inputs beyond the pair/residue counting of
   `two_thirds`. Zero-margin admissibility alone is cheap (`Omega_0(M) <= M + 1` for `M <= 7`).

Not covered, as before: prime-level residue data and the parity law (`referee_construct.md` §3.3),
transcendence of the actual angles, exact algebraic identities, and multi-place approximation.

---------------------------------------------------------------------------------------------

## 1. Exact form of the admissibility conditions

### 1.1 Setting

We use the notation of `outside.md` §1. A *profile* is `P = (p, e, (a_x), delta, k)`:
* distinct primes `p_j = 1 mod 4`;
* exponents `e_j >= 1`;
* rows `a_x in prod_j [0, e_j]`, normalised so that `min_x a_xj = 0` and `max_x a_xj = e_j`;
* positions `delta_x in [0, Delta]` and windings `k_x in Z`.

Put `W = sum e_j log p_j`, `Delta = C e^(-W/4)` and `u = delta/Delta in [0,1]^M`. The character
lattice is

```text
Lambda_P = { lambda in Q^M : sum lambda = 0,  v(lambda) := sum_x lambda_x a_x in Z^r }.
```

Write `w(v) = sum_j |v_j| log p_j` and `s_lambda := w(v(lambda)) - W/2` (the *slack*). For every
realisation `phi`, `outside.md` (1.2) gives

```text
chi(lambda) := <v(lambda), phi> = lambda.delta/2 + (pi/4) lambda.k .
```

**(AI)** We assume that the rows are affinely independent, i.e. `v(lambda) != 0` for `lambda != 0`.
For `C <= sqrt 2` this holds for every angle model. That is the pair-gap argument used in
`outside.md` Prop A.4, step 2. Under (AI) the realisation equations are solvable for every
`(delta, k)` (`outside.md` Theorem A.1(i)).

### 1.2 The exact inequalities

**Proposition 1.1.** Let `P` satisfy (AI) and `Pi(P) := prod_j (1 + 2/(sqrt p_j - 1)) < 1.75`. Then
`P` has a gap-admissible realisation (an angle model) iff every `lambda in Lambda_P \ {0}` satisfies

```text
d_A(chi(lambda)) >= arcsin(e^(-w/2))   and   d_D(chi(lambda)) >= arcsin(e^(-w/2)/sqrt 2),     w = w(v(lambda)),
```

where `d_A = dist(., (pi/2) Z)` and `d_D = dist(., pi/4 + (pi/2) Z)`. Put `y_lambda = |lambda.delta|/2`
and assume `y_lambda <= pi/16`. This holds whenever `||lambda||_1 Delta <= pi/4`, since
`y_lambda <= ||lambda||_1 Delta/4`. Then the pair of inequalities is **equivalent** to:

```text
(E) lambda.k even:        s_lambda >= 2 log( 2       / (C |lambda.u|) ) + eps_lambda,
(O) lambda.k odd:         s_lambda >= 2 log( sqrt 2  / (C |lambda.u|) ) + eps_lambda,
(N) lambda.k not in Z:    w(v(lambda)) >= -2 log sin( (pi/4) ||lambda.k|| - y_lambda )    (sufficient), and
                          w(v(lambda)) >= -2 log( sqrt 2 sin( (pi/4) ||lambda.k|| + y_lambda ) )   (necessary),
```

Here `eps_lambda = 2 log(y_lambda / sin y_lambda)` lies in `[0, y_lambda^2/2]`, and `||t||` is the
distance from `t` to the nearest integer.

*Proof.*
1. *Reduction to characters.* Necessity is `outside.md` Lemma 0.1 for `v = v(lambda)`.
   Sufficiency: the displayed inequalities for all `lambda` mean that `P` is character-admissible.
   Theorem A.1(ii) (refereed) then gives a gap-admissible realisation, because `Pi < 1.75`.
2. *Monotonicity.* `sin` is increasing on `[0, pi/2]` and every distance is at most `pi/4`. So
   `d >= arcsin(t)` iff `sin d >= t`.
3. *Case (E).* `chi(lambda) = +- y + (pi/2) m`. So `d_A = y` and `d_D = pi/4 - y`. The axis inequality
   `sin y >= e^(-w/2)` rearranges to (E), using `-2 log sin y = W/2 - 2 log((C/2)|lambda.u|) + eps`.
   The diagonal inequality follows from it: `e^(-w/2)/sqrt 2 <= sin(y)/sqrt 2 <= sin(pi/4 - y)` for
   `y <= pi/16`.
4. *Case (O).* `chi(lambda) = pi/4 +- y + (pi/2) m`. So `d_D = y`, and the diagonal inequality is
   (O). The axis inequality follows from it: `e^(-w/2) <= sqrt 2 sin y <= sin(pi/4 - y)` for
   `y <= pi/16`.
5. *Case (N).* `(pi/4) lambda.k` has distance `(pi/4) ||lambda.k||` from `(pi/4) Z`. So both
   distances are at least `(pi/4)||lambda.k|| - y`. Both inequalities hold if
   `sin((pi/4)||lambda.k|| - y) >= e^(-w/2)`. QED

So the windings enter only through the parities `lambda.k mod 2` (a factor `sqrt 2` in (E)/(O)) and
through the *excused* set `{lambda : lambda.k not in Z}`. Case (N) costs only a constant, about
`2 log(4D/pi)` for denominator `D`: it is automatic once `log p_min > 2 log(4D/pi) + 1`.

### 1.3 Fixed combinatorial data: a linear programme in the log-weights

**Proposition 1.2.** Fix the exponent rows `(a_x)`, the windings `k` and the normalised positions
`u`, with `lambda.u != 0` for every `lambda in Lambda \ {0}` with `lambda.k in Z`. Treat the
log-weights `w_j = log p_j` as real variables. Admissibility in the regime of Prop 1.1 means:

```text
sum_j w_j ( |v_j(lambda)| - e_j/2 )  >=  b_lambda + eps_lambda(W)      for every lambda != 0 with lambda.k in Z,
b_lambda = 2 log( c_lambda / (C |lambda.u|) ),   c_lambda = 2 (lambda.k even) or sqrt 2 (odd),
```

together with the (N)-conditions and `Pi < 1.75`. Since `0 <= eps_lambda <= y^2/2` and `y` decreases
in `W`, the admissible set lies between the polyhedron `P^-` (drop `eps`) and the polyhedron
`P^+(W_0) = {... >= b_lambda + (C ||lambda||_1 e^(-W_0/4))^2/8}` for weight vectors with `W >= W_0`.
So the minimal `W` for fixed combinatorial data is the value of a linear programme in `(w_j)`, up to an
additive error `O(e^(-W/2))`.

*Proof.* Combine Prop 1.1 with `W = sum e_j w_j` and `w(v(lambda)) = sum_j w_j |v_j(lambda)|`.
`eps_lambda` is increasing in `y`, and `y <= C ||lambda||_1 e^(-W/4)/4`. QED

The left side grows linearly in `||lambda||`. The right side grows like `log(1/|lambda.u|)`, which is
`O(log ||lambda||)` for almost every `u` (Khintchine). So for almost every `u` only finitely many
constraints bind. Positions, windings and patterns are the integer and nonlinear part of the problem.
§2 removes positions and parities exactly in the large-prime limit.

---------------------------------------------------------------------------------------------

## 2. The normalised problem

### 2.1 Definitions

Rows are `x = 0..M-1`. A *pattern* is a nonempty `S subset {1..M-1}`. It stands for the binary column
`1_S` and for its complement; these give identical constraints on zero-sum characters. There are
`2^(M-1) - 1` patterns.

* A *support* `T` is a set of patterns, and *weights* are `omega in [1, infinity)^T`.
  Write `Omega = sum omega_S`.
* For `lambda in Q^M` with `sum lambda = 0`, put `lambda(S) = sum_(x in S) lambda_x`,
  `N_omega(lambda) = sum_S omega_S |lambda(S)|` and

  ```text
  sigma_lambda(omega) = N_omega(lambda) - Omega/2 = sum_S omega_S (|lambda(S)| - 1/2).
  ```

* `L_T = {lambda : lambda(S) in Z for all S in T}`, and `G_T = L_T / Z^M_0`. (AI) on `T` means that the
  `|T| x (M-1)` matrix `[x in S]` has rank `M-1`; then `G_T` is finite.
* A *winding homomorphism* is any `h : G_T -> Q/Z`. A character `lambda in L_T \ {0}` is *relevant* for
  `h` if `h(lambda) = 0`.
* A requirement is `r(lambda) = rho0 + rho1 ||lambda||_1 / 2` with `rho0, rho1 >= 0`.

**Definition 2.1.** `(T, omega, h)` is *admissible* (for `r`) if (AI) holds on `T` and
`sigma_lambda(omega) >= r(lambda)` for every relevant `lambda`. It is *strictly admissible* if `>`
holds. Put

```text
Omega_r(M) = inf { Omega : (T, omega, h) admissible },      Omega_0 = Omega_(0,0),
Omega_+(M) = inf { Omega : (T, omega, h) strictly admissible for r = 0 },
Omega^prop_rho = Omega_(0, rho)     (norm-proportional margin),     Omega^unif_rho = Omega_(rho, 0).
```

The weight `omega_S` is the total log-weight of the primes carrying pattern `S`, in units of
`tau = log p_min`. Repeated columns merge, so `omega_S` may be any real `>= 1`.

### 2.2 Reductions and elementary bounds

**Lemma 2.2.**
(a) *Binary reduction.* A multi-level profile with weights `omega_j` maps to binary threshold patterns
`S_(j,t) = {x : a_xj >= t}`, `t = 1..e_j`, each of weight `omega_j`. Every character of the binary
profile is a character of the multi-level one, and its slack does not decrease. So multi-level
profiles never beat binary ones.

(b) *Windings are homomorphisms.* Every `h : G_T -> Q/Z` is `lambda -> lambda.k mod 1` for some
`k in Z^M` (computed exactly from the diagonal form, `core.Group.kappa_of`). Only `ker h` matters. A
smaller kernel imposes fewer constraints, so it suffices to range over homomorphisms with
inclusion-minimal kernels. If `G_T` is cyclic, an injective `h` exists, and only integral characters
remain.

(c) *Monotonicity.* `Omega_r(M') <= Omega_r(M)` for `M' <= M`.

(d) *Homogeneity.* `Omega_(t rho0, t rho1)(M) <= t Omega_(rho0, rho1)(M)` for `t >= 1`.

(e) *Plotkin and rank bounds.* Every admissible point has `Omega >= M - 1` and

```text
Omega >= 2(M-1)(rho0 + rho1)   (M even),        Omega >= 2M(rho0 + rho1)   (M odd).
```

Equality in the Plotkin bound forces every pattern to be balanced (`M` even) or near-balanced
(`M` odd), and every pair to have `sigma = rho0 + rho1` exactly.

*Proof.*

(a) Write `a_xj = sum_t [x in S_(j,t)]`. Then `lambda(a_j) = sum_t lambda(S_(j,t))`, so
`|lambda(a_j)| <= sum_t |lambda(S_(j,t))|` and `W` is unchanged. Integrality of all `lambda(S_(j,t))`
implies integrality of `lambda(a_j)`. (AI) passes to the binary profile because its `v` determines the
multi-level `v`.

(b) `Q/Z` is injective, so homomorphisms of the finite group `G_T` lift to `Q^(M-1) -> Q/Z`.
Explicitly, write `L_T = V . prod (1/d_i) Z` (exact diagonal form). Then `h(V n/d) = sum n_i eta_i/d_i`
is realised by `kappa = V^(-T) eta`. This is checked exactly in `test_core.py`.

(c) Delete rows, and delete the patterns that become constant (their `-omega/2` contributions
disappear). A zero-sum character on the kept rows, extended by zero, is a character of the original
profile with the same `|lambda(S)|` on the surviving patterns and the same `||lambda||_1`. Its lattice
condition and its relevance (restricted windings) are inherited. (AI) restricts.

(d) Replace `omega` by `t omega`: every `sigma` scales by `t`, and `omega >= 1` is kept.

(e) `Omega >= |T| >= M - 1` by (AI). For the Plotkin bound, sum over pairs `c = e_x - e_y`; these are
integral, hence always relevant, with `||c||_1/2 = 1`:

```text
sum_(x<y) sigma_(e_x - e_y) = sum_S omega_S ( |S|(M - |S|) - M(M-1)/4 ) <= Omega * max_s ( s(M-s) - M(M-1)/4 ),
```

which is `Omega M/4` (`M` even) or `Omega (M-1)/4` (`M` odd). It is at least `binom(M,2)(rho0 + rho1)`.
QED

### 2.3 Finite exact verification

**Lemma 2.3 (finite reduction).** Let `A` be any exact left inverse of the pattern matrix `Q`
(`A Q = I_(M-1)`). Put

```text
k_1 = max_S ( sum_x |A_xS| + |sum_x A_xS| ) / omega_S,     k_inf = max_S max( max_x |A_xS|, |sum_x A_xS| ) / omega_S.
```

Then `||lambda||_1 <= k_1 N_omega(lambda)` and `||lambda||_inf <= k_inf N_omega(lambda)` for all real
zero-sum `lambda`. Hence, if `rho1 k_1 < 2`, every `lambda` with `sigma_lambda < r(lambda)` satisfies
`N_omega(lambda) <= N_max := (Omega/2 + rho0)/(1 - rho1 k_1/2)`. It therefore lies in the box
`||lambda||_1 <= k_1 N_max`, `||lambda||_inf <= min(k_inf, k_1/2) N_max`.

A sharper `k_1` is
`max_s min { max_S |alpha_S|/omega_S : sum_S alpha_S 1_S = s + c 1 }` over sign vectors `s`. This is LP
duality for `max{ s.lambda : N <= 1 }`. Each `alpha` is rationalised and repaired to satisfy the
equality exactly (`nu1.py`).

*Proof.* `mu = (lambda_1..lambda_(M-1)) = A (lambda(S))_S` and `lambda_0 = -sum mu`. Take absolute
values. `||lambda||_inf <= ||lambda||_1/2` for zero-sum vectors. The dual bound is weak LP duality.
QED

So admissibility of a given `(T, omega, h)` is decided exactly by enumerating, for every class of
`G_T` in `ker h` and for the integral class, all lattice points in that box (`verify.py`,
`inner.py`). This uses integer arithmetic throughout. `test_core.py` checks the group computation,
the homomorphism realisation and the coercivity bounds against brute force on 386 supports, and the
verifier against naive enumeration on 53 cases. All pass.

**Lemma 2.4 (covering families are intrinsic).** Let `F` be a finite set of nonzero classes in
`Q^M_0/Z^M_0`. Call `F` *covering* if every homomorphism `<F> -> Q/Z` vanishes on some member of `F`.
If `F subset G_T` is covering, then every homomorphism `h` of `G_T` kills a member of `F`. The converse
also holds.

*Proof.* Restriction maps homomorphisms of `G_T` onto homomorphisms of `<F>`, because `Q/Z` is
injective. QED

So when the bad classes of a support cover, the cut "if `T subset Good(F)` then some member of `F`
satisfies its constraints" is valid for **every** support and every weight vector (milp2.py).

### 2.4 The large-prime limit

For `tau > 0` let `W^(tau)(M, C)` be the infimum of `W` over angle models satisfying (AI), with an
`M`-point cluster on an arc of length `<= C sqrt R`, and with every prime `p_j >= e^tau`. Let
`W_1^(tau)(M, C; Q_0)` be the same over the class `A_1(Q_0)` of `construct.md`, as corrected in
`referee_construct.md` §6:
* point residue data at prime powers `<= Q_0`, with `q = 2` included;
* axiom (R) for integral `c` with `v(c) != 0`;
* the level axiom (S|).

**Theorem 2.5.** Fix `M >= 3`, `C > 0` and `Q_0 >= 3`.

(a) `liminf_(tau -> infinity) W^(tau)(M, C) / tau >= Omega_0(M)`.

(b) `limsup_(tau -> infinity) W_1^(tau)(M, C; Q_0) / tau <= Omega_+(M)`.

(c) There is `C_0 = C_0(M, Q_0)` such that for `C >= C_0`, `lim_(tau) W_1^(tau)/tau = lim_(tau) W^(tau)/tau = Omega_0(M)`.

*Proof.*

(a) Take an angle model with all `p_j >= e^tau`.
* Split it into binary threshold patterns (Lemma 2.2(a)) and put
  `omega_S = (sum of log p_j over thresholds with pattern S)/tau >= 1`. Then `W = tau Omega` exactly.
* Its windings induce `h` on `G_T`. For relevant `lambda` (so `(pi/4) lambda.k in (pi/4) Z`),
  Lemma 0.1 gives `e^(-w/2)/sqrt 2 <= sin dist(chi(lambda), (pi/4)Z) <= dist(chi(lambda), (pi/4)Z)
  <= |lambda.delta|/2 <= ||lambda||_1 Delta/4`, whatever the size of `y`.
* So the multi-level slack satisfies `s_lambda >= -2 log^+(C ||lambda||_1/(2 sqrt 2))`. By
  Lemma 2.2(a), `tau sigma_lambda(omega) >= s_lambda`.
* Hence `(T, omega, h)` is admissible for `r_tau(lambda) = -(2/tau) log^+(C ||lambda||_1 / (2 sqrt 2))`.
* *Compactness.* Suppose, for contradiction, that points admissible for `r_(tau_i)` have
  `Omega <= Omega_0 - eps` along `tau_i -> infinity`.
  * There are finitely many supports and homomorphisms, and the weights lie in `[1, Omega_0]`.
  * Pass to a subsequence with fixed `(T, h)` and `omega_i -> omega*`.
  * For each fixed relevant `lambda`, `sigma_lambda(omega*) = lim sigma_lambda(omega_i) >= lim r_(tau_i)(lambda) = 0`.
  * So `(T, omega*, h)` is closed-admissible with `Omega <= Omega_0 - eps`, a contradiction.

(b) Take `(T, omega, h)` strictly admissible with `Omega <= Omega_+ + eps`.
* By Lemma 2.3 only finitely many relevant `lambda` have `sigma_lambda <= ||lambda||_1/(2k_1)`. So
  there is `eta > 0` with `sigma_lambda >= eta (1 + ||lambda||_1)` for all relevant `lambda`.
* For large `tau`, choose distinct primes `p_S = 1 mod 4` with `log p_S in [tau omega_S, tau omega_S + 1]`
  (prime number theorem in progressions), and windings `k` realising `h`.
* Use the residue design of `construct.md` Lemma 4.1 (auxiliary primes `l_x > Q_0`, independent of
  `tau`, with `L = max l_x`).
* Draw `u` uniformly. For relevant integral `c`, the requirement
  `dist(c.delta/2, (pi/4) Z) >= arcsin(L^(n/2) e^(-w/2))` (`construct.md` (5.1)) implies (R) and
  Lemma 0.1. It fails with probability `<= (4/C) L^(n/2) e^(-s_c/2)(n Delta + pi)`
  (`lower.md` Lemma 5.3; re-derived in §5). Here `s_c >= tau sigma_c - |T|/2`.
* For relevant rational `lambda` only Lemma 0.1 is needed, with the same bound and `L = 1`, and density
  `<= 2M/(Delta ||lambda||_1)`.
* The number of `lambda in L_T` with `||lambda||_1 <= n` is polynomial in `n`. The bound decays like
  `e^(-tau eta n/2)`. So the total failure probability tends to 0 as `tau -> infinity`.
* *Excused (non-relevant) classes.* They have `||lambda.k|| >= 1/D` and `w(v(lambda)) >= tau`, so (N)
  holds deterministically while `||lambda||_1 Delta <= pi/(2D)`. Beyond that, `||lambda||_1` is at
  least `e^(W/4)/(2DC)`, and `sigma_lambda >= ||lambda||_1/k_1 - Omega/2` makes their total
  contribution negligible.
* Completion: `outside.md` Theorem A.1(ii), with `Pi - 1 <= 2.2 |T| e^(-tau/2) -> 0`.
* The points are distinct by the pair cases.
* Hence `W <= tau(Omega_+ + eps) + |T|`.

(c) Let `(T, omega, h)` attain `Omega_0`. The minimum exists: there are finitely many `(T, h)`, and
each feasible set is closed with compact sublevel sets.
* The set `Z` of relevant `lambda` with `sigma_lambda = 0` is finite (Lemma 2.3), and
  `s_lambda >= -|T|/2` on `Z`.
* For generic `u`, `min_(lambda in Z) |lambda.u| > 0`. The failure probability of the finitely many
  `Z`-constraints is `<= K(M, omega, L)/C`.
* The others are handled as in (b).
* So for `C >= C_0` and large `tau` the union bound is `< 1`, and `W <= tau Omega_0 + |T|`.
* Combine with (a); `A_1 subset A_0`. QED

So at fixed `M` the normalised optimisation is exactly the large-prime asymptotics. `Omega_0` and
`Omega_+` can differ: at `M = 5` they are 4 and 7. Then the limit depends on `C`.

---------------------------------------------------------------------------------------------

## 3. Exact data for `M = 4..10`

### 3.1 Methods and certification levels

* **Search** (`milp3.py`). The MILP has semi-continuous weights `omega_S in {0} cup [1, U]` and
  support indicators `y_S`, with the rank cut `sum y_S >= M-1`. Its constraints are:
  * integral characters (lazy, with their `S_M`-orbits);
  * covering-family cuts (Lemma 2.4) found by the exact inner check (`inner.py`), together with their
    random `S_M`-images.

  The solution is rationalised and checked exactly. Either some homomorphism avoids all bad classes
  (feasible: done), or the bad classes cover (new cut). An earlier formulation with explicit integer
  winding variables (`milp.py`) is correct but much slower. It agrees on `M <= 6`.
* **E (exact).** The value is attained by an exactly verified point, and it equals a proven lower
  bound. The proven lower bounds are:
  * the rank bound or the Plotkin bound (Lemma 2.2(e));
  * the exhaustive certificate `exhaust.py`. It enumerates all supports with `|T| < v` up to `S_M`
    (orderly generation), and all minimal-kernel homomorphisms. For each it gives an exact Farkas,
    Gordan or LP-dual certificate (Fractions) over a finite set of relevant characters (all integral
    `c` with entries in `[-2,2]`, and all rational lifts with `||lambda||_inf <= 3/2`).
* **S (solver-optimal).** The value is attained by an exactly verified point, and HiGHS reports a
  dual bound equal to it on a valid relaxation.
* Every reported optimum was re-verified independently with `verify.py` and its windings
  (`reverify.py`; all pass).

### 3.2 Table

**Table 3.1.**

| `M` | `Omega_0` | `Omega_+` | `Omega^unif_(1/2)` | `Omega^prop_(1/2)` | `Omega^prop_1` | `Omega^prop_2` | `Omega^prop_4` | Plotkin(`rho=1`) |
|---|---|---|---|---|---|---|---|---|
| 4 | 3 E | 3 E | 3 E | 6 S | 10 S | 16 S | 32 S | 6 |
| 5 | 4 E | 7 E | 7 E | 9 S | 10 E | 20 E | 40 E | 10 |
| 6 | 7 E | 7 E | 7 E | 9 S | 10 E | 20 E | 40 E | 10 |
| 7 | 8 S | TBD | | | 14 E | 28 E | | 14 |
| 8 | TBD | `<= 14` E (Thm 4.1) | | | 14 E | 28 E | | 14 |
| 9 | | | | | 18 E | | | 18 |
| 10 | | | | | TBD | | | 18 |

Uniform margins `Omega^unif_rho` (`rho = 1, 2`) equal Plotkin for `M = 4, 5, 6` (6/12, 10/20, 10/20; E).

**Optimal profiles** (rows `0..M-1`; patterns listed by their rows, all weights 1 unless stated).
* `M = 4`, `Omega_+ = 3`: the Hadamard core `{12, 13, 23}`. `G = Z/2`, killed by the winding
  `k = (0,0,0,1)`.
* `M = 5`, `Omega_0 = 4`: `{14, 24, 34, 1234}`. Some pairs have `sigma = 0` exactly; it is not strictly
  admissible. `Omega_+ = 7`: `{12, 13, 23, 24, 124, 34, 134}`, with `G` trivial.
* `M = 6`, `Omega_0 = Omega_+ = 7`: `{12, 134, 234, 135, 235, 145, 245}`. `G = Z/4` (cyclic), with an
  injective winding.
* `M = 7`, `Omega_0 = 8`: `{1234, 135, 245, 126, 36, 346, 256, 1456}`. `G = Z/2`, minimum slack 0.
* Margins, `M = 5, 6`: all balanced (`M = 6`) or 2|3 (`M = 5`) cuts, weight `rho`.
* Margins, `M = 7`: fourteen 3|4 cuts, every pair separated 8 times, `G` trivial.
* Margins, `M = 8`: two disjoint Hadamard cores (14 balanced cuts, every pair separated 8 times,
  `G = Z/4`). This is the skew double core of §4.
* Margins, `M = 9`: eighteen 4|5 cuts, every pair separated 10 times, `G` trivial.
* `M = 4`, `rho1 = 1`: `{2 (w 1), 12, 13, 23 (w 3)}`, value 10 against Plotkin 6. The aligned 4-flat
  `(1,-1,-1,1)` of the order-4 Hadamard core has `sigma = omega/2` but needs `2 rho`.

---------------------------------------------------------------------------------------------

## 4. The skew double core

### 4.1 Theorem

**Theorem 4.1.** Let `K` be a skew-Hadamard matrix of order `M >= 8`, i.e. `K + K^T = 2I` and
`K K^T = M I`. Normalise `K` so that `K(x,0) = 1` for all `x`: replace `K` by `DKD`, which keeps
skewness. Write `K = I + S` with `S^T = -S`, so `S(a,0) = 1` and `S(0,a) = -1` for `a >= 1`. Put

```text
S_a = { x : K(x,a) = -1 },     S'_a = S_a xor {0, a}        (a = 1..M-1),     T = {S_a} cup {S'_a},
```

with all weights 1. Then:

(i) the `2(M-1)` cuts are distinct and balanced, and every pair of rows is separated by exactly `M` of
them;

(ii) `G_T` is cyclic of order `M/2`, generated by `lambda* = (2/M) 1 - 2 e_0`;

(iii) `sigma_c >= ||c||_1 / 2` for every nonzero integral zero-sum `c`, with equality at pairs.

Consequently, with an injective winding (e.g. `k = e_1`), `(T, 1, h)` is admissible for
`r = ||lambda||_1/2`. So for every skew-Hadamard order `M >= 8` and every `rho >= 1`:

```text
Omega^prop_rho(M) = Omega^unif_rho(M) = 2(M-1) rho    (the Plotkin bound),     Omega_+(M) <= 2(M-1).
```

For odd `M` with `M + 1` a skew-Hadamard order `>= 8`, restricting to `M` rows (Lemma 2.2(c)) keeps
all `2M` cuts. Then `Omega^prop_rho(M) = 2M rho`, again the Plotkin bound.

*Proof.*

(i) *Cuts.* `K(a,a) = 1` and `K(0,a) = -1`, so `0 in S_a` and `a notin S_a`. Thus
`S'_a = (S_a \ {0}) cup {a}` has the same size `M/2`. (Column `a >= 1` is orthogonal to the all-ones
column `0`, so it is balanced.)
* `S_a = {x : K(x,a) != K(x,0)}` is the cut of columns `a` and `0` of `K`. Likewise
  `S'_a = {x : K^T(x,a) != K^T(x,0)}`, the cut of columns `a` and `0` of the Hadamard matrix `K^T`.
  To check this, compare row by row: `x = 0` (equal entries), `x = a` (`1` against `-1`), and
  `x != 0, a`, where `K^T(x,a) = -K(x,a)` and `K^T(x,0) = -1`.
* Two rows `x != y` are separated by the cut of columns `a, 0` iff
  `H(x,a)H(y,a) != H(x,0)H(y,0)`. Orthogonality of rows `x, y` of a Hadamard matrix `H` makes this
  happen for exactly `M/2` columns `a`. So each family separates every pair `M/2` times.
* *Distinctness.* `|S_a xor S_b| = M/2` for `a != b`. Moreover:
  * `S'_a xor S'_b = S_a xor S_b xor {a, b}`;
  * `S'_a xor S_b = S_a xor S_b xor {0, a}`;
  * `S'_a xor S_a = {0, a}`.

  All these sizes lie in `{2, M/2 - 2, M/2, M/2 + 2}`, never 0 or `M` for `M >= 8`. So the `2(M-1)`
  cuts are distinct as partitions.

(ii) *The group.* Let `lambda in L_T`.
* Then `lambda(S'_a) - lambda(S_a) = lambda_a - lambda_0` is an integer for every `a`, so
  `lambda = lambda_0 1 + c` with `c in Z^M`.
* `sum lambda = 0` gives `M lambda_0 in Z`.
* `lambda(S_a) = lambda_0 M/2 + c(S_a) in Z` gives `lambda_0 in (2/M) Z`.
* Conversely, every such `lambda` lies in `L_T`.
* So `lambda -> lambda_0 mod 1` is an isomorphism `G_T -> (2/M) Z/Z`.
* With `k = e_1`, `lambda*.k = 2/M` has order `M/2`, so `h` is injective.

(iii) *The margin.* For zero-sum `c`, `|c(S_a)| = |T_a|/2` with `T_a = (K^T c)_a`, because
`sum_x c_x K(x,a) = -2 c(S_a)`. Also `c(S'_a) = c(S_a) - c_0 + c_a`. One checks:

```text
2 sigma_c = Phi(c) - 2(M-1),     Phi(c) := ||u||_1 + ||v||_1,     u = K^T c = (I - S) c,     v = K c~,
```

where `c~` is `c` with `c_0` negated.
* *Identity.* `T_a + 2c_0 - 2c_a = -(Kc)_a + 2c_0 = -(K c~)_a`, because `K(a,0) = 1`. Also
  `u_0 = sum c = 0` and `v_0 = (Kc)_0 - 2c_0 = 0`.
* *Norms.* `||u||_2^2 = ||v||_2^2 = M ||c||_2^2`. Also `|u_a|, |v_a| <= ||c||_1 = n`.
* *Defect identity.* `||u||_1 = ||u||_2^2/n + (1/n) sum_a |u_a|(n - |u_a|)`, and the same for `v`.

We need `Phi >= 2(M-1) + n`.

*Amplitude `A = ||c||_inf >= 2`.*
* `c = K u/M` and `c~ = K^T v/M`, so `||u||_1, ||v||_1 >= MA`.
* Then `Phi >= 2MA >= 2M - 2 + AM >= 2M - 2 + n`, since `n <= AM`.

*Ternary `c`, `n = ||c||_1` (even).* `Phi = 2M + D/n` with
`D = sum_a [|u_a|(n - |u_a|) + |v_a|(n - |v_a|)]`. We need `D >= n(n-2)`.
* All `u_a, v_a` are even, as sums of `n` terms `+-1`.
* So every nonzero term `|t|(n - |t|)` with `|t| in [2, n-2]` is at least `2(n-2)`.
* `n = 2` is trivial.

*Case `c_0 = 0`.* For `a` in the support, `u_a = c_a - (Sc)_a` and `v_a = c_a + (Sc)_a`. So
`u_a + v_a = 2c_a = +-2`.
* Both `|u_a|, |v_a|` lying in `{0, n}` would force `n = 2`.
* So each of the `n` support rows contributes `>= 2(n-2)`, and `D >= 2n(n-2)`.

*Case `c_0 = eps != 0`.* Let `X = supp(c) \ {0}` and `r_a = sum_(x in X \ {a}) S(a,x) c_x`. For
`a in X`:

```text
u_a = (c_a - eps) - r_a,       v_a = (c_a - eps) + r_a.
```

The rows `a in X` with `c_a = -eps` number exactly `n/2`. There `u_a + v_a = -4 eps`.
* *`n >= 6`:* both `|u_a|, |v_a|` in `{0, n}` is impossible. So `D >= (n/2) 2(n-2) = n(n-2)`.
* *`n = 4`:* every row `y notin supp(c) cup {0}` has `u_y = -(eps + m_y)` and `v_y = m_y - eps`, with
  `m_y = sum_(x in X) S(y,x) c_x` odd. Exactly one of `|u_y|, |v_y|` equals 2. So each such row
  contributes 4, and `D >= 4(M - 4) >= 8 = n(n-2)` for `M >= 6`.

Equality holds at pairs, where `Phi = 2M`. QED

### 4.2 Checks

`skew_double.py` builds the normalised Paley skew matrices for `q = 7, 11, 19, 23, 43` and the
doubling skew matrices of orders 8, 16, 32 (exact). It then checks:
* (i) and (ii);
* the identity `2 sigma_c = Phi(c) - 2(M-1)` on random integer `c`;
* `min (Phi - (2M - 2 + n))` over exhaustively many ternary `c`, and over hill-climbed vectors with
  entries in `[-2,2]`.

| `M` | `G_T` | exhaustive ternary range | min `Phi - (2M-2+n)` | exact verifier (all classes) |
|---|---|---|---|---|
| 8 (Paley, doubling) | `Z/4` | all supports (553) | 0 (pairs) | admissible, injective winding |
| 12 | `Z/6` | all supports (36,894) | 0 (pairs) | admissible, injective winding |
| 16 (doubling) | `Z/8` | support `<= 8` (536,110) | 0 | not run |
| 20 | `Z/10` | `<= 6` (402,325) | 0 | not run |
| 24 | `Z/12` | `<= 6` (1,378,114) | 0 | not run |
| 32 (doubling) | `Z/16` | `<= 4` (108,376) | 0 | not run |
| 44 | `Z/22` | `<= 4` (408,199) | 0 | not run |

Hill climbing (300 starts each) never goes below 0. These checks are consistent with the proof;
the proof does not depend on them.

### 4.3 Are Hadamard-like profiles optimal?

Fixed-support optima at `M = 8` (`hadamard_eval.py`, exact LP value and exact dual certificate; the
twin supports have trivial `G_T`):

| support | `|T|` | zero margin (`rho0 = 1/1000`) | `rho1 = 1` | `rho1 = 2` |
|---|---|---|---|---|
| Paley core alone | 7 | infeasible (bad classes cover; Theorem B type) | infeasible | infeasible |
| core + skew twins (`construct.md` profile, optimal weights) | 15 | `16.014` | 30 | 44 |
| core + random capacity-two twins (`lower.md` type) | 15 | `28.014` | 44 | 64 |
| **skew double core** | 14 | `<= 14` | **14** | **28** |
| optimum | | `[8 (S; 7 E), 14]` | 14 | 28 |

* For the margins that constructions need, **the double Hadamard core is optimal, and single-core
  profiles are not**. They lose a factor of about 2 to 3.
* For zero margin, Hadamard-type profiles are far from optimal at small `M`. The optima use about
  `M + O(1)` cuts.
* At `M = 4` the single Hadamard core is optimal for zero margin, but not for margins.
* Random double cores (a second core from a random row permutation) can fail. At `M = 12` the first
  two random choices had integral violations; a screened search found an admissible one
  (`double_search.py`). The skew construction succeeds for every skew order.

---------------------------------------------------------------------------------------------

## 5. Fake circles at `(12 + o(1)) M log M`

**Lemma 5.0 (one character; `lower.md` Lemma 5.3, re-derived).** Let `c != 0` be zero-sum and
integral, `n = ||c||_1`, `0 < Delta <= 1`, and `delta` uniform on `[0, Delta]^M`. Then
`Pr[dist(c.delta/2, (pi/4) Z) < arcsin(min(1, eps))] <= 4 eps (n + pi/Delta)`.

*Proof.*
* Conditioning on all coordinates but one with `c_x != 0`, `c.delta/2` has density `<= 2/Delta`.
* Its range is an interval of length `n Delta/2` containing 0.
* At most `2n Delta/pi + 2` bad windows meet it, each of length `2 arcsin(min(1, eps)) <= pi eps`.
QED

**Theorem 5.1.** Let `K` be a skew-Hadamard matrix of order `M >= 8`, `C in (0,1]` and `Q_0 >= 3`.
Let `L` be the `M`-th prime `= 1 mod 4` exceeding `Q_0` (for the class `A_0` put `L = 1` and use no
residue data). Let `r = 2(M-1)` and

```text
P >= P_* := max{ 160000 M^4 L^2 / C^2,  21 r^2 },     p_1 < ... < p_r primes = 1 mod 4 in [P, P e^(1/(4r))].
```

Then there is an angle model in `A_1(Q_0)` (in `A_0` if `L = 1`) with the following properties:
* the profile of Theorem 4.1, with the primes assigned to the `2(M-1)` cuts, `e = 1`, and windings
  `k = e_1`;
* `M` distinct points on an arc of length `<= C sqrt R`;
* Lemma 0.1 for every monomial;
* the residue design of `construct.md` Lemma 4.1, satisfying (R) for every integral character;
* `W <= r (log P + 1/(4r))`.

With the window lemma of `construct.md` §5 (`P <= 2 P_**`) one gets, as `M -> infinity` along skew
orders:

```text
A_1(2M):    W <= (12 + o(1)) M log M + 4 M log(1/C)      (L ~ 2 M log M),
A_1(M^A):   W <= (8 + 4A + o(1)) M log M + 4 M log(1/C)   (A >= 1),
A_0:        W <= (8 + o(1)) M log M + 4 M log(1/C).
```

For all `M`: Paley orders `q + 1` with `q = 3 mod 4` prime are skew, and `q + 1 = (1 + o(1)) M` by the
prime number theorem in progressions. Monotonicity (`lower.md` Lemma 5.2) transfers the bound.

*Proof.* We follow `construct.md` Theorem 2 step by step, with Theorem 4.1 in place of its Theorem 1.

1. *Weights.* `w_S = log p_S in [tau, tau + 1/(4r)]` with `tau = log P`. For integral `c`,
   `s_c = sum_S w_S(|c(S)| - 1/2) >= tau sigma_c - r/(8r) >= tau n/2 - 1/8`, by Theorem 4.1(iii).
2. *Requirement.* For integral `c` we demand `dist(c.delta/2, (pi/4) Z) >= arcsin(L^(n/2) e^(-w(v(c))/2))`.
   This implies (R) for every unit (`construct.md` step 2, with the referee's
   `sqrt 2 D_u <= 2 L^(n/2)`), and Lemma 0.1 on these characters.
3. *Union bound.* `e^(-w/2) = (Delta/C) e^(-s_c/2)`. Lemma 5.0, `(n Delta + pi) <= (pi + 1) n` and
   step 1 give

   ```text
   Pr[c fails] <= (4(pi+1) e^(1/16)/C) n (L^(1/2) e^(-tau/4))^n.
   ```

   There are at most `(2M)^n` integral `c` with `||c||_1 = n`. With `x = 2M L^(1/2) e^(-tau/4) <= sqrt(C)/10`
   (this is exactly `P >= 160000 M^4 L^2/C^2`):

   ```text
   sum_c Pr <= (17.64/C) sum_(n>=2) n x^n <= (17.64/C) * 2.346 x^2 <= 0.414.
   ```

4. *Excused rational classes.* `G_T = Z/(M/2)` and the winding is injective. So every non-integral
   `lambda in L_T` has `||lambda.k|| >= 2/M`.
   * Moreover `w(v(lambda)) >= tau`, so `arcsin(e^(-w/2)) <= 1.04/sqrt P <= 1.04/(400 M^2)`.
   * If `||lambda||_1 Delta <= pi/M`, then `dist(chi(lambda), (pi/4)Z) >= pi/(2M) - pi/(4M)`, which
     exceeds this. So Lemma 0.1 holds deterministically.
   * Otherwise `||lambda||_1 > pi/(M Delta)`, which is astronomically large. Then
     `N(lambda) >= ||lambda||_1/2`: for every real `lambda` by Theorem 4.1(iii), scaling and
     density.
   * The density of `lambda.delta/2` is `<= 2M/(Delta ||lambda||_1)`. The number of
     `lambda in (2/M) Z^M` with `||lambda||_1 <= nu` is `<= (M nu + 1)^M`.
   * So these classes contribute at most `sum_(nu >= pi/(M Delta)) (M(nu+1)+1)^M 4.2 M^2 e^(-tau nu/4) < 10^(-100)`.
5. *Completion.* `P >= 21 r^2` gives `1.33(Pi - 1) < 1`. `outside.md` Theorem A.1(ii) supplies
   Lemma 0.1 for monomials outside `L`, and the referee's Fubini argument puts this off a Haar-null
   set.
   * The residue data are those of `construct.md` Lemma 4.1. Every `p_S > Q_0`, so (S|) is vacuous.
   * 2-adic moduli are covered as in `referee_construct.md` R4.
   * The rank is `M`, so (R) is only ever required for `v(c) != 0`.
6. *Points.* The pair cases give distinct points on an arc of angular width `Delta`, i.e. of length
   `C sqrt R`.
7. *Size.* `W <= r(tau + 1/(4r))`. With `tau = log P_* + O(1)`,
   `log P_* = 4 log M + 2 log L + 2 log(1/C) + log 160000`. Also `L = (2 + o(1)) M log M` for
   `Q_0 = 2M` (a sieve gives `L = 89, 137, 257, 613` for `M = 8, 12, 20, 44`), or `L ~ M^A`. QED

**Corollary 5.2.** For every `C > 0` and `A >= 1`:
* `limsup_M W_min^(A_1(M^A))(M, C)/(M log M) <= 8 + 4A`;
* `limsup W_min^(A_1(2M))/(M log M) <= 12`;
* `limsup W_min^(A_0)/(M log M) <= 8`.

Combined with `two_thirds` (`>= 3` on `A_1(2M)`, with level data):

```text
3  <=  liminf W_min^(A_1(2M))/(M log M)  <=  limsup  <=  12.
```

**Corollary 5.3 (class barrier).** Consider an argument whose inputs about a cluster are:
* exponent identities;
* all Lemma 0.1 gaps;
* (R) at moduli `<= 2M`;
* Haar-generic angle properties.

No such argument proves `M <= (1/6 - eps) log R/loglog R`. (`log R = W/2 <= (6 + o(1)) M log M`.) The
best constant for this class lies in `[1/6, 2/3]`. With moduli `<= M^A` the lower end is `1/(4 + 2A)`.

**Exact instances** (`instance_skew.py`). All values are exact. Primes are the `r` smallest
`= 1 mod 4` above `P_*`, logarithms are 50-digit. "Refined" re-runs the union bound with the exact
slack of every character with `||c||_1 <= 6` (`M = 8`) or `<= 4` (`M = 12`), plus the theorem's tail,
and bisects on `P`.

| `M` | class | `C` | `L` | `W` (theorem) | `W/(M log M)` | refined `W` | refined ratio | `b = 8` Paley profile (construct / referee) |
|---|---|---|---|---|---|---|---|---|
| 8 | `A_1(16)` | 1 | 89 | 409.89 | 24.64 | 356.59 | 21.44 | 878.3 (52.8) |
| 8 | `A_1(16)` | 1/2 | 89 | 429.30 | 25.81 | 375.35 | 22.56 | |
| 8 | `A_0` | 1 | 1 | 284.21 | 17.08 | 230.91 | 13.88 | |
| 12 | `A_1(24)` | 1 | 137 | 698.78 | 23.43 | 624.06 | 20.93 | 1611.8 (54.1) |
| 12 | `A_0` | 1 | 1 | 482.30 | 16.17 | 407.58 | 13.67 | |

In every row the theorem's union bound is `0.414 < 1/2`, the window has log-width `<= 1/(4r)`, and
`1.33(Pi - 1) <= 1.5 * 10^(-3)`.

**Growth of the constant** (`rates_skew.py`; arithmetic on the proved bound, `C = 1`, exact `L` by
sieve):

| `M` | 8 | 12 | 20 | 44 | 104 | 1004 | 10004 | 100004 | 1000004 |
|---|---|---|---|---|---|---|---|---|---|
| `A_0` | 17.10 | 16.18 | 15.20 | 14.01 | 13.03 | 11.46 | 10.60 | 10.08 | 9.73 |
| `A_1(2M)` | 24.65 | 23.44 | 22.24 | 20.64 | 19.36 | 17.20 | 16.00 | 15.27 | 14.76 |

`construct.md` (`rates.py`) has `50.0, 41.7, 35.1, 32.8` (`A_1`) and `49.2, 37.3, 28.3, 25.1` (`A_0`)
at `M = 12, 104, 10004, 1000004`. The limits are `12` and `8`, approached at the rate
`(loglog M + O(1))/log M`.

---------------------------------------------------------------------------------------------

## 6. The trend of `W/(M log M)`, and where the gap `[3, 12]` lives

### 6.1 Normalised data as `W/(M log M)` proxies

Write `tau = t log M` (all primes about `M^t`). Then `W ~ t Omega log M`, so
`W/(M log M) ~ t Omega/M`. Two finite-`M` proxies follow.

* **Necessary side (any construction with all primes `>= M^(1-o(1))`):** `t >= 1` and margins
  `>= -o(1)`. This gives `K_nec(M) := Omega_0(M)/M`:
  `0.75, 0.80, 1.17, 1.14` for `M = 4..7`. This is only a heuristic for the `M log M` scale: o(M) small
  primes can saturate lattices cheaply (Lemma 2.2(b)), and the compactness argument of Theorem 2.5(a)
  is not uniform in `M`.
* **Union-bound side (random positions).** The characters of 1-norm `n` number about `M^n`, and with
  the residue design each carries a factor `L^(n/2) = M^(n/2 + o(n))`. So `Pr < 1` needs, up to
  `1 + o(1)` factors,
  `sigma_lambda >= rho ||lambda||_1/2` with `rho = 4/t` (`A_0`) or `6/t` (`A_1(2M)`). This gives

  ```text
  K_UB(M) = min_rho (4/rho) Omega^prop_rho(M)/M   (A_0),        (6/rho) Omega^prop_rho(M)/M   (A_1).
  ```

  Table 3.1 gives `K_UB^(A_0) = 8, 8, 6.67, 8, 7, 8` for `M = 4..9`, i.e. exactly `8(M-1)/M` (even
  `M >= 6`) or `8` (odd `M`). For `A_1`: `12(M-1)/M` or `12`.

**Proposition 6.1 (Plotkin floor of the union-bound route).** Every normalised profile with
`sigma_lambda >= rho ||lambda||_1/2` for all pairs has `(4/rho) Omega >= 8(M-1)` (`M` even), resp.
`>= 8M` (`M` odd). So no profile can push the union-bound route below `8` (`A_0`) or `12` (`A_1`).
Theorem 5.1 attains both.

*Proof.* Lemma 2.2(e) with `rho0 + rho1 = rho` on pairs. QED

The translation from "union bound over random positions" to "margin `rho = 4/t`" is a counting
heuristic. It is made rigorous only for the skew double core, in Theorem 5.1. Prop 6.1 itself is
exact.

### 6.2 What would move the constant

* **Upper side (below 12).** Random positions force every pair to have slack `>= 2 log(M^2)`, because
  the closest of `M` random points are `M^(-2)` apart. Structured positions (gaps about `1/M`) would
  need only about `2 log M` for most pairs. If all characters could be handled that way, the Plotkin
  floor would drop to about 4 (`A_0`). No such construction is known. Higher characters then need
  pigeonhole-aware position design, of Sidon/`B_h` type.
* **Lower side (above 3).** `two_thirds` uses pairs and residue pigeonhole only. The normalised data
  show that character admissibility *alone* costs about one weight unit per row
  (`Omega_0(M) <= M + 1` for `M <= 7`). So the gap between 3 and 12 lies in the positional
  simultaneity (union bound versus packing) and in the residue design. Exact identities, prime-level
  residues and transcendence remain outside the class.

### 6.3 Small-`M` `W_min` brackets (the class contains actual circles)

For small `M` the large-prime fakes are far from optimal: actual circles are members of the class
(`construct.md` Lemma 1.2 / `lower.md` Prop 1.2), with small primes.

| `M` | `C` | upper bound `W_min <= W(actual)` | source | rigorous lower bound |
|---|---|---|---|---|
| 4 | 1 | 19.84 (`n = 412333805`) | `growth/families.md` data | 7.01 (rank: `log 5 + log 13 + log 17`, `C <= sqrt 2`) |
| 5 | 3 | 16.16 (`n = 10414625`) | same | `two_thirds` `F*(5)`: `W >= 7.63 - 16 log 3 < 0`, vacuous |
| 7 | 5.27 | 19.28 (`n = 235370525`) | same | vacuous |
| 8 | 5.674 | 31.13 (`5^3 13^2 17^2 89 109 557`) | same | vacuous |
| 10 | 8.24 | 20.89 (`n = 1176852625`) | `referee_outside.md` | vacuous |

* At `C = 1` the exact `two_thirds` function `F*(M)` gives `W >= 8F*(M)/M = 4.16, 7.63, 9.86, 12.31, 14.10, 19.04, 23.91` for
  `M = 4..10` (`tt/finite_inequality.py`, copied from the repository).
* For `C <= sqrt 2` the rank bound `sum_(j < M) log q_j` (`q_j` the primes `= 1 mod 4`) is larger at
  these `M`.
* So `W/(M log M)` of the best known actual clusters is between about `0.9` (`M = 10`) and `3.6`
  (`M = 4`). At small `M` the class optimum is dominated by small primes, and says nothing about the
  asymptotic constant.

---------------------------------------------------------------------------------------------

## 7. Status, attribution, and what is not covered

* *New here:*
  * the exact reformulation in the form of Props 1.1-1.2;
  * the normalised MILP with covering cuts (Lemma 2.4) and its large-prime theorem (2.5);
  * the exact data of Table 3.1;
  * Theorem 4.1 (skew double core, Plotkin-tight) and Theorem 5.1 (`12`, `8`, `8 + 4A`);
  * Prop 6.1.

  The `two_thirds`, `outside.md` (A.1), `construct.md` (Lemma 4.1, the instance template) and
  `lower.md` (Lemma 5.3, Lemma 5.2) inputs are cited.
* *Relation to `construct.md` Theorem 1.* That theorem uses one Hadamard core with `b` unflipped copies
  and one flipped copy per row. The double core replaces the `b` copies by a second, transposed core
  whose cuts are the first core's cuts with the two rows `0` and `a` flipped. So every row gets `M-1`
  "flips". The skew structure enters as in item 392 / `construct.md`, through
  `u_a + v_a = 2(c_a - c_0)`.
* *Prime-level residues* (`referee_construct.md` §3.3). The parity law is not addressed. The fakes of
  Theorem 5.1 use the same point-level design as `construct.md`, so they are not claimed to satisfy
  it. Whether `W_min/(M log M)` is unbounded for the prime-level class remains open.
* *Solver dependence.* Entries marked S, and the M = 7, 8 zero-margin values, rely on HiGHS. No
  theorem in this note uses them.

---------------------------------------------------------------------------------------------

## 8. Files (`round5/exact/`)

| file | content | type |
|---|---|---|
| `core.py` | patterns, exact diagonal form and `G_T`, homomorphisms and windings, exact left inverse, coercivity, box enumeration | exact |
| `verify.py`, `inner.py`, `nu1.py` | exact admissibility verifier (Lemma 2.3); all-class analysis; covering test; exact `k_1` via dual LPs | exact |
| `test_core.py` | brute-force tests of `core`/`verify` (386 + 53 cases) | exact |
| `milp.py`, `milp2.py`, `milp3.py` | MILP with windings; Benders/covering cuts; with support indicators | search (HiGHS) + exact checks |
| `exhaust.py` | exhaustive exact lower-bound certificates (orderly generation up to `S_M`, Farkas/Gordan/dual) | exact |
| `reverify.py` | independent exact re-verification of every reported optimum | exact |
| `fixed.py`, `hadamard_eval.py` | fixed-support LP over weights and windings; Hadamard-type supports | exact certificates |
| `structure.py` | pattern sizes, pair separations, Hadamard cores inside a support | exact |
| `skew_double.py` | Theorem 4.1 checks, `M = 8..44` | exact |
| `double_test.py`, `double_search.py` | random double cores (`M = 8, 12`) | exact on survivors |
| `instance_skew.py`, `rates_skew.py` | Theorem 5.1 instances (`M = 8, 12`) and growth table | exact; rates are arithmetic |
| `tt/finite_inequality.py` | copy of `two_thirds/finite_inequality.py` (small-`M` `F*`) | exact integers |
| `logs/`, `results/` | run logs and JSON optima | |

The repository `/home/user/Jarnik` was not modified.
