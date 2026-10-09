# Round 5: the angle-model threshold W_min(M, C) is Theta(M log M)

Question. Let `W_min(M, C)` be the least `W = log N` over angle models carrying an `M`-point
cluster on an arc of length `C sqrt(R)`, with genuine primes `= 1 mod 4`, a genuine exponent
profile, real angles, all Lemma 0.1 gaps, and residue-collision strengthening at small moduli.
Known: `3 M log M - O(M) <= W_min <= (30 + o(1)) M^2 log M`. Does `W_min/(M log M)` tend to infinity?

Checks are in `round5/checks/`. Each script runs from that directory with `python3 <script>`, and its
output is stored next to it as `*_out.txt`. Conventions are those of `round4/outside.md`:
* `W = log N`, `R = N^(1/2)`;
* the arc has angular width `Delta = C e^(-W/4)`;
* `w(v) = sum_j |v_j| log p_j`;
* `L` is the character lattice, and `chi_P` is (1.2) of `outside.md`.

*Theorem*, *Lemma* and *Proposition* mean proved here with all quantifiers, or cited with the exact
source and re-checked. *Exact check* means integer arithmetic. *Illustration* or *evidence* means
numerics that no proof uses.

---------------------------------------------------------------------------------------------

## 0. Verdict

**The target lower bound is false. `W_min(M, C) = Theta(M log M)`.** The class of angle models with
residue data (§1) carries fake circles at the rate `W = O(M log M)`. Since the class contains every
actual circle, the local + character + residue + GAP toolkit is **sharply capped at
`log R/loglog R`**: no argument whose inputs lie in the class can prove `M = o(log R/loglog R)`.
This settles problems (P1) and (P2) of `round4/outside.md` §6: (P1) affirmatively, (P2) negatively.

| result | statement | status |
|---|---|---|
| Def. 1.1, Prop. 1.2 | the class `A_X(M, C)`: angle models with point-residue data at odd moduli `<= X`; it contains every actual normalised cluster | proved |
| Prop. 1.3 | for `X >= 2M`, every member satisfies `3 M log M <= W + (14 + 4 log^+ C) M` (`two_thirds` inside the class) | proved (re-reading of `two_thirds`) |
| Remark 1.4 | the literal A.3 strengthening (factor `Q_M` for every character) forces `W >= (8 - o(1)) M^2`, but it is **not** a valid input for actual circles | proved; explains the definition |
| **Lemma S** | **`l^1`-stability of the Paley transform.** For prime `q = 3 mod 4`, `q >= 43`, `M = q + 1`, every nonzero zero-sum integer `c` other than a pair, a signed column or a signed two-column half-sum has `||H_A^T c||_1 >= (17/16) M`. If some `|c_x| >= 2`, then `>= 2M` | proved; new. Exact census: the true constant is about `M/3` |
| Prop. 3.2 | canonical five-or-more-copy Paley profile: `B - r/2 >= b/2 - 2` (pairs), `>= b/2 + h - 4` (ternary support `h >= 4`), `B >= M A (b/2 - 1)` (amplitude `A >= 2`) | research-notes item 392, re-proved and re-checked exactly |
| Prop. 3.3 | every one-flip, label-capacity-two Paley profile: slack lower bounds for all characters, from Lemma S | proved; new |
| Lemma 4.1 | residue design: point residues from auxiliary Gaussian primes; every generous collision modulus divides `Im(Omega_c^4)/4` | proved; new; exact check |
| **Theorem 5.1** | **`W_min^(X=2M)(M, C) <= (48 + o(1)) M log M`** (canonical route). Every capacity-two assignment gives `(2b + o(1)) M log M` for every fixed `b >= 35`. Residue data up to modulus `M^A` give `(16 + 32 A + o(1)) M log M` | proved (asymptotic in `M`, for fixed `C`) |
| Cor. 5.4 | no argument with inputs in `A_X` proves `M <= (1/24 - eps) log R/loglog R` for all large `R`. So the best constant of the class lies in `[1/24, 2/3]` | proved |
| §7.1 | **attribution.** Research-notes item 394 (simultaneous character admissibility of canonical Paley at `W = 20 M log M`) combined with `outside.md` Theorem A.1 already gave fake circles at `O(M log M)` for the inputs without residues. Neither `outside.md` nor its referee noticed this | observation |
| §8, (P3) | not covered: residue data for the **individual primes `pi_j`** (point residues that factor through prime residues). A quadratic-residue parity obstruction blocks the residue design used here | open; precise obstruction stated |

What remains for a growth improvement is listed in §8. It is: exact algebraic identities among the
points; transcendence of the prime angles; and residue data of the individual Gaussian primes (P3).

---------------------------------------------------------------------------------------------

## 1. The class, its validity, and why the strengthening must be per-assignment

Notation follows `outside.md` §1. `pi_j` is a fixed Gaussian prime over `p_j`. For a row vector
`c in Z^M` with `sum c = 0`, put `v(c) = sum_x c_x a_x in L`, and `n = ||c||_1`.

**Definition 1.1 (the class `A_X(M, C)`).** Fix `X >= 1`, and for an odd prime `q` put
`A_q = floor(log X / log q)`. A member is a tuple `(p, e, (a_x), phi, theta_0, delta, k, rho)` with
the following properties.

(a) *Profile.* `p_1, ..., p_r` are distinct primes `= 1 mod 4`, and `e in Z_(>=1)^r`. Put
    `N = prod p_j^(e_j)` and `W = log N`. The rows `a_x in prod_j [0, e_j]`, `x = 1..M`, are
    pairwise distinct, and `min_x a_xj = 0`, `max_x a_xj = e_j` for every `j`.

(b) *Angles.* `phi in R^r`, `theta_0 in R`, `delta in [0, Delta]^M` with `Delta = C e^(-W/4)`, and
    `k in Z^M`, satisfying (1.1) of `outside.md`:
    `<2a_x - e, phi> = theta_0 + delta_x + (pi/2) k_x`.

(c) *(G) Gaps.* Every `v in Z^r \ {0}` satisfies both inequalities of Lemma 0.1.

(d) *Residue data.* For every odd prime `q <= X` with `q` not dividing `N`, and every `x`, a unit
    `rho_x(q)` of `Z[i]/q^(A_q)` with `Norm rho_x(q) = N mod q^(A_q)`. For every `p_j <= X`, and
    every `x`, a unit `rho_x(p_j)` of `Z[i]/p_j^(A_(p_j))` with norm `N/p_j^(e_j)` (level residues,
    as in `two_thirds` Lemma 3.3 (S|)).

(e) *(R) Strengthened character gaps.* For every `c in Z^M` with `sum c = 0` and `v(c) != 0`:

    ```text
    |sin(c . delta / 2)|  >=  2^(-1/2) m_c(rho) e^(-w(v(c))/2).                          (1.1)
    ```

    Here `m_c(rho) = prod_q q^(a_q(c))`. For odd `q <= X` not dividing `N`,
    `a_q(c) = max{a <= A_q : prod_x rho_x(q)^(c_x) = 1 mod q^a}`. For `q = p_j <= X`, `a_q(c) = 0`
    unless `c = e_x - e_y` with `a_xj = a_yj`, in which case
    `a_q(c) = max{a : rho_x(p_j) = rho_y(p_j) mod p_j^a}`.

`W_min^X(M, C)` denotes the infimum of `W` over `A_X(M, C)`.

The exponent identities (gcd norms, cut identity, levels and chains) hold automatically in every
profile. Theorem F of `walsh-extraction.md` is (G) applied to `v in L`. The Haar-generic properties
(class (IV) of `outside.md` Cor. A.4) are a property of an argument, not of a model. They are dealt
with in Cor. 5.4.

**Proposition 1.2 (validity).** Take a normalised actual cluster `z_1, ..., z_M` on an arc of length
`<= C sqrt R` (`two_thirds` Lemma 1.1), with arc parametrisation `arg z_x = theta* + delta_x`,
`delta_x in [0, Delta]`, and `k_x` from (1.1). Put `rho_x(q) := z_x mod q^(A_q)` and
`rho_x(p_j) := r_x mod p_j^A`, where `r_x = z_x / (pi_j^(a_xj) conj(pi_j)^(e_j - a_xj))`. Then the
tuple lies in `A_X(M, C)` for every `X`.

*Proof.*
* (a) and (b) are `outside.md` §1, and (G) is Lemma 0.1.
* *(R) for odd `q` not dividing `N`.* Let `c != 0` with `sum c = 0` and `v = v(c) != 0`. Put
  `U = prod_(c_x>0) z_x^(c_x)` and `V = prod_(c_x<0) z_x^(-c_x)`. Both have norm `N^(n/2)`.
* *`U` and `V` against `A_v`.* Compare exponents of `pi_j` and `conj(pi_j)` in `U` and `V`. They
  differ exactly by `v_j`. So with `G = gcd(U, V)` one gets `U/G = u_1 A_v` and `V/G = u_2 conj(A_v)`
  for units `u_1, u_2`. Then `|G| = N^(n/4) e^(-w(v)/2)`.
* *Collisions divide `D = U/G - V/G`.* If `prod_x rho_x(q)^(c_x) = 1 mod q^a`, then
  `U = V mod q^a`, because `V` is a unit mod `q`. Since `q` is prime to `Norm G` (which divides
  `N^(n/2)`), `q^a` divides `D`. Distinct `q` are coprime, so `m_c` divides `D`.
* *`D != 0` and `(1+i) | D`.* `D = 0` would make `A_v` associate to `conj(A_v)`. That is
  impossible, since `A_v` is a conjugate-primitive nonunit. Every odd-norm Gaussian integer and every
  unit is `1 mod (1+i)`, so `(1+i)` divides `D`. Hence `|D| >= sqrt(2) m_c`.
* *Conclusion.* `U/V = prod z_x^(c_x) = e^(i c.delta)`, because `sum c = 0` and all points lie at
  `R e^(i(theta* + delta_x))`. So `|U - V| = 2 N^(n/4) |sin(c.delta/2)|`. Also
  `|U - V| = |G| |D| >= N^(n/4) e^(-w/2) sqrt(2) m_c`. This is (1.1).
* *Same-level pairs at `p_j`.* `two_thirds` Lemma 3.3 (S|) gives `p_j^a | c_xy`. Lemma 3.1 there
  multiplies the collisions, and the factor `(1+i)` divides `c_xy` (Lemma 3.3 (R)). Then
  `|z_x - z_y| = |g_xy| |c_xy|` and `|g_xy| = R e^(-d_xy/2)` give (1.1) for `c = e_x - e_y`.
QED

**Proposition 1.3 (the 2/3 bound holds on the class).** If `X >= 2M`, every member of `A_X(M, C)`
satisfies `3 M log M <= W + (14 + 4 log^+ C) M`.

*Proof.* The proof of `two_thirds` Theorem 7.5 uses the configuration in exactly two places. It needs
no other input.
1. *Inequality (2.5).* In the class, (1.1) for `c = e_x - e_y` and `|delta_x - delta_y| <= Delta`
   give `log(sqrt(2) m_xy) <= log C - W/4 + d_xy/2`. Summing over pairs and using the cut identity
   `sum_(x<y) d_xy = (M^2/4) W - Q` gives
   `sum_(x<y) log(sqrt2 m_xy) + Q/2 <= (M/8) W + binom(M,2) log C`.
   This is (2.5) with `log|c_ij|` replaced by `log(sqrt2 m_ij)`.
2. *Inequality (3.6).* This is a pigeonhole on residues with the class counts of Lemma 3.2. Those
   counts depend only on the norm conditions in (d). The collisions it uses are mod `q^a` with
   `(q - chi(q)) q^(a-1) < M`, hence `q^a < 2M <= X`. The ramified factor is the `sqrt 2` in step 1.
The chain lemma and Sections 4-8 there are pure combinatorics and prime sums. QED

So the lower bound `W_min^X >= 3 M log M - O(M)` holds for `X >= 2M`, as the question states.

**Remark 1.4 (the universal strengthening of A.3 is not a valid input).** Theorem A.3's fakes satisfy
`dist(chi(v), (pi/4) Z) >= arcsin(Q_M e^(-w(v)/2))` for every character, with
`Q_M = lcm(1, ..., 2M)`. Suppose this condition is imposed on a class.
* For pairs it gives `d_xy >= W/2 + 2 log Q_M - 2 log(C/2)`.
* Summing over pairs against `sum d_xy <= (M^2/4) W` gives
  `W >= 4(M-1)(psi(2M) - log(C/2)) = (8 - o(1)) M^2`.
* So under the universal strengthening, `W_min/(M log M) -> infinity` holds trivially.
* But the condition fails for actual clusters. The 10-point record cluster has `W = 20.9`, while
  `8 M^2 = 800`. The valid input is (1.1), with the collision modulus of the actual residues.
* The question's implication "`W_min/(M log M) -> infinity` gives a growth improvement" needs a valid
  class. That is why Definition 1.1 quantifies over residue data **existentially**: an actual circle
  supplies one residue assignment.

---------------------------------------------------------------------------------------------

## 2. Lemma S: `l^1`-stability of the Paley transform

Let `q = 3 mod 4` be prime, `M = q + 1`, and `H` the normalised Paley-I matrix of
`round4/flipped.md` Lemma 6.2. Its rows are `x in F_q` and `infinity`. Its columns are the constant
column and the labels `a in F_q`, with `H(x,a) = -chi(a - x) - [a = x]` and `H(infinity, a) = 1`.
For `c in Z^M` with `sum c = 0`, put

```text
T = T(c) = (sum_x c_x H(x,a))_(a in F_q),    F(c) = ||T||_1,    n = ||c||_1,    m_2 = ||c||_2^2.
```

The **exceptional set** `E` consists of:
* the pairs `+-(e_x - e_y)`;
* the signed columns `+-H_a`;
* the signed half-sums `(s H_a + s' H_b)/2`, with `a != b` and `s, s' = +-1`.

Each of these has `F = M`.

**Lemma 2.1 (elementary facts; any normalised Hadamard matrix).** Let `c != 0` be an integer vector
with `sum c = 0`.
1. *(Parseval)* `sum_a T_a^2 = M m_2`. Every `T_a` is even, `|T_a| <= n`, and `n` is even.
2. *(defect identity)* `F = M m_2 / n + (1/n) sum_a |T_a| (n - |T_a|)`. In particular `F >= M`.
3. *(inversion)* `F >= M max_x |c_x|`, and `F >= n` (and `F >= M m_2/n` by 2).
4. *(integrality)* Let `kappa = #{a : |T_a| > n/2}`. Then `|F - kappa n| <= 2 (F - M)`.
5. *(column peeling)* If `c != +-H_a`, then `F >= min(2|T_a|, 2M)`.
6. *(half-sum peeling)* Let `a != b`, `s = sgn T_a`, `s' = sgn T_b` (with `sgn 0 := 1`), and
   `g = (s H_a + s' H_b)/2`.
   If `c != g`, then `F >= M + phi(|T_a|) + phi(|T_b|)`, where `phi(t) = t - |t - M/2|`.

*Proof.*
1. Parseval is `H^T H = M I`, with `T` at the constant column equal to `sum c = 0`. Parity holds
   because `T_a = sum c_x = 0 mod 2`.
2. The identity is `F - sum T^2/n = (1/n) sum |T|(n - |T|)`, and `m_2 >= n`.
3. `c = (1/M) H T` (including the zero constant coordinate) gives `|c_x| <= F/M`. Cauchy-Schwarz
   gives `m_2 >= n^2/|supp c| >= n^2/M`, so `F >= M m_2/n >= n`.
4. Put `u_a = |T_a|/n` in `[0,1]`. Then `sum_a u_a (1 - u_a) = (F - M m_2/n)/n <= (F - M)/n`, and
   `sum_a u_a - kappa = sum_(u<=1/2) u - sum_(u>1/2) (1 - u)`. Both sums are at most
   `2 sum u(1-u)`. Multiply by `n`.
5. Let `s = sgn T_a`. For `b != a`, `T_b(c) = T_b(c - s H_a)`, and `T_a(c - s H_a) = T_a - s M`.
   Hence `F(c) = |T_a| + F(c - sH_a) - ||T_a| - M|`. Since `c - sH_a != 0` is integral and zero-sum,
   `F(c - sH_a) >= M` by 2. So `F(c) >= M + |T_a| - ||T_a| - M| = min(2|T_a|, 2M)`.
6. The same computation with `T(g) = (M/2)(s e_a + s' e_b)` and `g` integral. QED

**Lemma 2.2 (transposed prime-field uncertainty).** Let `c in Z^M`, `sum c = 0`, `c != 0`, and
`s' = |supp c cap F_q| <= (q-1)/2`. Then `#{a in F_q : T_a != 0} >= (q+1)/2 - s'`.

*Proof.* This is the argument of `flipped.md` Lemma 6.2 / research-notes item 393, for `H^T` instead
of `H`.
1. *Reduction.* Write `c'` for the finite part. If `q^t` exactly divides `c'`, then it divides
   `c_infinity = -sum c'` and every `T_a`. So divide by `q^t`. We may assume `c' != 0 mod q`
   (`c' = 0` forces `c = 0`).
2. *The polynomial.* For `a` outside `supp c'`, `T_a = c_infinity - sum_x c_x chi(a - x)`. Euler's
   criterion gives `chi(y) = y^d mod q` (`d = (q-1)/2`, also at `y = 0`). So
   `T_a = -G(a) mod q`, where `G(Y) = sum_x c_x (Y - x)^d - c_infinity` has degree `<= d`.
3. *`G` is nonconstant.* For `1 <= k <= d`, the coefficient of `Y^k` is
   `binom(d,k) sum_x c_x (-x)^(d-k)`, and `binom(d,k) != 0 mod q`. If all of them vanished, then
   `sum_x c_x x^m = 0` for `m = 0..d-1`. Since `s' <= d`, the Vandermonde matrix on `supp c'` is
   invertible, which would force `c' = 0 mod q`.
4. *Count.* So `G` has at most `d` roots. Every `a` outside `supp c'` with `T_a = 0` is a root.
   Hence at most `s' + d` labels have `T_a = 0`. QED

*Exact check* (`check_lemmaS.py`, Part 1): exhaustive at `q = 7, 11` (supports up to 8 resp. 6), and
43,735 vectors in all, including random ones up to `q = 83` with entries divisible by `q` and `q^2`.
The minimum of `#nonzero + s' - (q+1)/2` is 1, so the bound is off by one from the truth, as in
item 393.

**Lemma S.** Let `q = 3 mod 4` be prime with `q >= 43`, put `M = q + 1`, and let `theta = 1/16`.
Every nonzero zero-sum `c in Z^M` outside `E` satisfies `F(c) >= (1 + theta) M`. If moreover
`max |c_x| >= 2`, then `F(c) >= 2M`.

*Proof.* The second claim is Lemma 2.1(3). So let `c` be ternary, with support `h = n` (even). If
`h = 2`, `c` is a pair. Assume `h >= 4`, `c` not in `E`, and `F < (1 + theta) M`. We derive a
contradiction.

*Step 1 (windows).* By 2.1(4), `|F - kappa h| < 2 theta M`. With `M <= F`, this gives
`kappa h in ((1 - 2theta) M, (1 + 3theta) M)`, so `kappa >= 1`.

*Step 2 (`kappa >= 3`: the small range).* Then `h < (1 + 3theta) M/3 = (19/48) M`, which is below
`(q-1)/2`.
* By Lemma 2.2, at least `M/2 - h` labels have `T_a != 0`.
* By Parseval, at most `M/h` labels have `|T_a| = h`.
* Every other nonzero `T_a` is even with `2 <= |T_a| <= h - 2`, so `|T_a|(h - |T_a|) >= 2(h - 2)`.
* By 2.1(2),

  ```text
  F - M >= f(h) := (2(h-2)/h) (M/2 - h - M/h).                                          (2.1)
  ```

* `f(4) = M/4 - 4` and `f(6) = (4/3)(M/3 - 6)`.
* For `8 <= h <= 19M/48`, the factor `2(h-2)/h` is at least `3/2`, and `g(h) = M/2 - h - M/h` is
  concave. So `f(h) >= (3/2) min(3M/8 - 8, 5M/48 - 48/19)`.
* For `M >= 44` all three bounds are `>= M/16`. This contradicts `F < (17/16) M`.

*Step 3 (`kappa = 1`).* Then `h > (1 - 2theta) M` and one label has `u_a > 1/2`.
* From the proof of 2.1(4), `1 - u_a <= 2 sum u(1-u) < 2theta M/h < 2theta/(1 - 2theta)`.
* Hence `|T_a| > (1 - 4theta) M`.
* Since `c != +-H_a`, peeling (2.1(5)) gives `F >= min(2|T_a|, 2M) > (2 - 8theta) M = (3/2) M`.
  Contradiction.

*Step 4 (`kappa = 2`).* Then `h > (1 - 2theta) M/2` and two labels have `u > 1/2`.
* The same estimate gives `1 - u < 4theta/(1 - 2theta)` for both.
* So `|T_a|, |T_b| > (1 - 6theta) M/2`.
* Since `c` is not the half-sum `g` of 2.1(6),
  `F >= M + phi(|T_a|) + phi(|T_b|) > M + 2 (1/2 - 6theta) M = (2 - 12theta) M = (5/4) M`.
  Contradiction. QED

**Remarks 2.3.**
1. *What is new.* Research-notes item 393 classified the exact flats (`F = M`, ternary `c`, prime
   `M > 12`). Lemma S is the **quantitative** version that `outside.md` §2.2 asked for, for all
   integer `c`: "if `||H^T c||_1 <= (1 + gamma) M`, then `c` is a pair, a column or a half-sum".
   It holds with `gamma = 1/16` at every prime order `q >= 43`.
2. *Which steps use Paley.* Steps 1, 3 and 4 hold for **every** normalised Hadamard matrix. Near-flat
   characters always have `||c||_1 = M/kappa + O(theta M)`, and in the windows `kappa = 1, 2` they
   are exactly columns and half-sums. Only the small range `kappa >= 3` uses the prime field
   (Lemma 2.2). This is where Sylvester matrices differ: their aligned flats have `||c||_1 = M/2^j`.
3. *Sharpness (evidence; `check_lemmaS.py` Part 3, `check_medium.py`).* The minimum of
   `(F - M)/M` over non-exceptional `c` was searched by:
   * exhaustive search over support 4, and up to support 6 for `q <= 23`;
   * structured families;
   * hill climbing, including band-restricted climbing over `||c||_1 in [M/4, 2.5M]`.

   The minimum found is `0.36, 0.33, 0.40, 0.35, 0.39, 0.40, 0.38, 0.38, 0.41` for
   `q = 43, 47, 59, 67, 71, 79, 83, 103, 107`. It is always attained at support 4. In the medium
   range it is `>= 0.82`. So the proved `1/16` is conservative by a factor of about 6.
4. *Exact checks of every inequality used* (`check_lemmaS.py` Part 2, 31,990 vectors over 8 orders):
   the defect identity, integrality, `F >= n`, both peelings, and (2.1).

---------------------------------------------------------------------------------------------

## 3. Fully flipped Paley profiles and their character slacks

### 3.1 The profile

Fix `q` (Paley, `M = q + 1`), an integer `b >= 4`, and a map `alpha : rows -> F_q` (labels) with
label capacity `|alpha^(-1)(a)| <= 2`.

* *Columns.* Each label `a` has `b` physical columns, all equal to `H_a` before flipping.
* *Flips.* Row `x` is flipped in one physical column `f(x)` of label `alpha(x)`, and distinct rows
  use distinct physical columns.
* *Siblings.* Each flipped column `f(x)` has an unflipped *sibling* `s(x)` in the same label. This is
  possible since `b - 2 >= 2`.
* *Matrices.* Let `S in {+-1}^(M x r)` be the resulting sign matrix, `r = b(M-1)`. The rows are
  `a_x = (1 + S_x)/2 in {0,1}^r`, and `e = 1`.

For `c in Z^M` with `sum c = 0` put `v(c) = sum_x c_x a_x = (1/2) c^T S` (integral, since every
`c^T S_j` is even) and

```text
B(c) = ||v(c)||_1 = (1/2) sum_j |c^T S_j|,
2B(c) = b F(c) + sum_x ( |T_alpha(x) - 2 c_x H(x, alpha(x))| - |T_alpha(x)| ).                (3.1)
```

(3.1) holds because the unflipped columns of label `a` contribute `|T_a|` each, and the flipped one
of row `x` contributes `|T_a - 2c_x H(x,a)|`. It is checked exactly in `check_fake.py` (a).

**Lemma 3.1 (twins, saturation, normalisation).**
1. `S_(f(x)) - S_(s(x)) = -2 H(x, alpha(x)) e_x` for every row `x`. Hence `rank S = M`.
2. `L = {v(c) : c in Z^M, sum c = 0}`, and `v(c) != 0` for `c != 0`.
3. Every column is nonconstant (`M >= 4`), so (a) of Definition 1.1 holds.

*Proof.* Item 1 is immediate. For item 2, if `v = sum lambda_x a_x` with rational `lambda`,
`sum lambda = 0`, is integral, then `v_(f(x)) - v_(s(x)) = -lambda_x H(x, alpha(x))` forces
`lambda in Z^M`. This is walsh-extraction Lemma F.1 / `outside.md` Lemma A.7(ii). Item 3: a balanced
`+-1` column with one entry flipped is not constant. QED

So every rational character is integral, and Theorem F with any denominator is vacuous on these
profiles. With weights `w_j in [tau, tau + 1/r]` and `W = sum w_j <= r tau + 1`, the **slack** of a
character satisfies

```text
s_c := w(v(c)) - W/2 >= tau (B(c) - r/2) - 1/2.                                             (3.2)
```

### 3.2 The canonical assignment (research-notes item 392)

The *canonical* assignment is `alpha(x) = x` for finite `x`, together with `alpha(infinity) = a*`
in a second copy of a fixed label `a*`. Here `H(x, x) = -1`.

**Proposition 3.2.** For `b >= 5` and every nonzero zero-sum `c in Z^M`:
1. pairs: `B - r/2 >= b/2 - 2`;
2. ternary `c` of support `h >= 4`: `B - r/2 >= b/2 + h - 4`;
3. `A = max|c_x| >= 2`: `B >= M A (b/2 - 1)`.

*Proof (re-derived from item 392).* Let `c` be ternary with support `h`, and write
`F = M(1 + delta)`. Each supported flip changes `|c^T S_j|` by `+-2`. Let `k` count the supported
flips that lower it. By (3.1), `B - r/2 = b/2 + h - 2k + (bM/2) delta`.

*The skew structure.* Put `D = diag(+1 at infinity, -1 at finite rows)` and `K = DH`. Then
`K = I + J` with `J` skew: `J(x,a) = chi(a - x)` on the finite block, `J(infinity, .) = 1`,
`J(., infinity) = -1`. With `d = Dc`, `T = K^T d = d - Jd`.

*Lowering flips.* The finite flip at row `x` (label `x`) lowers `|T_x|` iff `d_x T_x > 0`, and
`d_x T_x = 1 - d_x (Jd)_x`. A lowering flip with `|T_x| = h` needs `d_x J_xy d_y = -1` for every
other supported `y`. Two such rows `x, y` would give `d_x J_xy d_y = d_y J_yx d_x = -1`,
contradicting `J_yx = -J_xy`. So at most one finite lowering flip has `|T_x| = h`, and the infinity
flip adds at most one more.

*Ternary `h >= 4`.* The other `k - 2` lowering flips sit at distinct labels with
`2 <= |T_x| <= h - 2`. Each contributes `>= 2(h-2)` to `sum_a |T_a|(h - |T_a|) = M h delta`. So
`delta >= 2(h-2)(k-2)_+/(M h)`, and

```text
B - r/2 >= b/2 + h - 4 + (k - 2)_+ (b(1 - 2/h) - 2) >= b/2 + h - 4    (b >= 5, h >= 4).
```

*Pairs.* `delta = 0` and `k <= 2`.

*Amplitude `A >= 2`.* Lemma 2.1(3) gives `F >= MA`, and the `M` flips change `2B` by at most
`2 sum|c_x| <= 2MA`. QED

*Exact check* (`check_canonical.py`): `q = 11, 19, 23, 43, 47, 59`, `b = 5, 8`. Pairs are checked
exhaustively, ternary support 4 exhaustively (and support 6 at `q = 11`), together with columns,
half-sums, 20,000 random ternary vectors per order and 5,000 vectors of amplitude `>= 2`. All three
bounds hold, and the first two are attained (minimum margin 0).

### 3.3 Arbitrary capacity-two assignments (from Lemma S)

**Proposition 3.3.** Let `q >= 43` and let `alpha` be any capacity-two assignment. For nonzero
zero-sum `c`, with `h = |supp c|` and `n = ||c||_1`:

| character | lower bound for `2B - r` |
|---|---|
| pair | `b - 4` |
| `+-H_a` | `b + 2M - 8` |
| half-sum `(s H_a + s' H_b)/2` | `b + M - 16` |
| any other ternary `c` | `b(M/16 + 1) - 2h >= (b/16 - 2) M + b` |
| amplitude `A >= 2` | `b(AM - M + 1) - 2AM >= AM(b - 4)/2 + b` |

*Proof.* Use (3.1). Each flip term is at least `-2|c_x|`, so the flip sum is `>= -2n`.
* *Pairs.* `F = M`, and only the two supported flips move, each by `+-2`.
* *Columns `c = sH_a`.* `T = sM e_a`. A flip at a row with `alpha(z) != a` adds `+2`, and there are
  at least `M - 2` of these. The at most two rows with `alpha(z) = a` lose `2` each.
* *Half-sums.* The support `E` has `M/2` rows, and `T = (M/2)(s e_a + s' e_b)`. A flip at a row of
  `E` with label outside `{a, b}` adds `+2`. At most four rows have `alpha in {a, b}`.
* *Other ternary `c`.* Lemma S gives `F - M + 1 >= M/16 + 1`.
* *Amplitude `A >= 2`.* `F >= AM` and `n <= AM`. QED

*Exact check* (`check_fake.py`, random capacity-two assignments, `q = 43, 47, 59, 67, 83`):
* identity (3.1) and the twin structure are verified;
* the pair, column and half-sum bounds are attained exactly;
* every one of about 5,700 other vectors per order satisfies its row of the table.

The smallest observed value of `(2B - r)/(bM)` among the other characters is `0.33-0.39`, against
`1/16 - 2/b` guaranteed.

---------------------------------------------------------------------------------------------

## 4. Residue design

**Lemma 4.1.** Let `X >= 3`. Let `ell_1, ..., ell_M` be distinct primes `= 1 mod 4` in `(X, Y]`,
and `varpi_x` a Gaussian prime over `ell_x`. Let `N` be prime to every odd `q <= X`. For every odd
prime `q <= X` pick `beta_q` with `Norm beta_q = N mod q^(A_q)`; this is possible since the norm is
onto `(Z/q^a)^*`. Put

```text
rho_x(q) := beta_q varpi_x conj(varpi_x)^(-1)   mod q^(A_q).
```

1. `Norm rho_x(q) = N mod q^(A_q)`, so this is residue data as in Definition 1.1(d).
2. *Generous collision modulus.* For nonzero zero-sum `c in Z^M`, let
   `m_c^gen := prod_q q^(a_q)` with `a_q = max{a <= A_q : prod_x rho_x(q)^(c_x) = some unit mod q^a}`.
   Put `Omega_c = prod_(c_x>0) varpi_x^(c_x) prod_(c_x<0) conj(varpi_x)^(-c_x)`. Then
   `Im(Omega_c^4) != 0` and `m_c^gen` divides `Im(Omega_c^4)/4`. In particular
   `m_c^gen <= Y^(2||c||_1)/4`, and also `m_c^gen <= e^(psi(X)) <= e^(1.04 X)`.
3. `m_c(rho) <= m_c^gen` (the exact-unit collision of Def. 1.1(e) is a special case).

*Proof.*
1. Item 1 holds since `Norm(varpi/conj(varpi)) = 1`.
2. *Reduction to `Omega`.* `prod_x rho_x^(c_x) = Omega/conj(Omega)` (`sum c = 0`), and `Omega` is a
   unit mod `q`. For a unit `eta`, `Omega/conj(Omega) = eta mod q^a` iff `q^a` divides
   `Omega - eta conj(Omega)`. Writing `Omega = P + iQ`, this means `q^a` divides `Q`, `P`, `P - Q`
   or `P + Q` for `eta = 1, -1, i, -i` respectively.
3. *One unit per prime.* For a fixed odd `q`, at most one `eta` can occur at level `>= 1`. Two would
   give `q | (eta - eta')`, and `|eta - eta'|^2 <= 4`.
4. *Divisibility.* So the prime powers split into four pairwise coprime groups, each dividing one of
   `Q, P, P - Q, P + Q`. Hence `m_c^gen` divides `PQ(P^2 - Q^2) = Im(Omega^4)/4`.
5. *Nonvanishing.* All four factors are nonzero: otherwise `Omega` would be associate to
   `conj(Omega)`. But `v_(varpi_x)(Omega) = c_x^+` and `v_(varpi_x)(conj Omega) = c_x^-`, since the
   `ell_x` are distinct split primes, and these differ for `c_x != 0`.
6. *Bounds.* `|Omega|^4 = prod ell_x^(2|c_x|)`. The second bound is `m_c^gen | lcm(1..X)`, together
   with Rosser-Schoenfeld `psi(x) < 1.03883 x`.

Item 3 is clear. QED

*Exact check* (`check_fake.py` (e), `X = 60`, `M = 24`): the pair collision criterion
`q^a | Im(varpi_x conj(varpi_y))` is verified in 5,520 cases. The divisibility
`m_c^gen | Im(Omega_c^4)/4` is verified on 366 random characters.

---------------------------------------------------------------------------------------------

## 5. The theorem

Say a member satisfies **(R_gen)** if, for every nonzero zero-sum `c`,

```text
dist(c . delta / 2, (pi/4) Z)  >=  arcsin( min(1, m_c^gen e^(-w(v(c))/2)) ).                 (5.1)
```

With `k = 0`, `c.delta/2 = chi_P(v(c))`. Since `m^gen >= 1`, (5.1) implies (G) for every `v in L`.
Since `|sin t| >= sin dist(t, (pi/4)Z)` and `m_c <= m_c^gen`, it also implies (R).

**Theorem 5.1.** Fix `C > 0`. As `M -> infinity`:

(a) `W_min^(2M)(M, C) <= (48 + o(1)) M log M`. The witnesses are canonical Paley profiles with
    `b = 8`, and they satisfy (R_gen).

(b) For every fixed `b >= 35` and every capacity-two one-flip assignment on prime Paley orders, there
    are members of `A_(2M)(M, C)` with `W <= (2b + o(1)) M log M`.

(c) For fixed `A >= 1`: `W_min^(M^A)(M, C) <= (16 + 32 A + o(1)) M log M`.

(d) Without residue strengthening (only (G), i.e. `X < 3`): `<= (20 + o(1)) M log M`. This is
    research-notes item 394 plus `outside.md` Theorem A.1, with the diagonal gaps added as noted below.

Together with Prop. 1.3: `3 M log M - O(M) <= W_min^X(M, C) <= K M log M` for `2M <= X <= M^A`.

**Lemma 5.2 (restriction).** If `A_X(M', C)` has a member with `e = 1` and `W'`, then for every
`M <= M'`, `A_X(M, C)` has a member with `W <= W'`.

*Proof.* Keep any `M` rows.
* *Profile.* Delete the columns that are constant on them. `N` and `W` decrease, and the rows stay
  distinct and normalised.
* *Angles.* Restrict `phi`. The deleted columns add a common constant to (1.1), which goes into
  `theta_0`. `delta` and `k` are unchanged, and `Delta` grows.
* *(G).* It restricts.
* *Residue data.* Multiply by a fixed `beta` of norm `N_new/N`. This leaves every
  `prod rho_x^(c_x)` unchanged, since `sum c = 0`.
* *(R).* For `c` supported on the kept rows, `v(c)` vanishes on the deleted columns, so
  `w(v(c))`, `sin(c.delta/2)` and `m_c` are unchanged. QED

**Lemma 5.3 (one character).** Let `c != 0` be zero-sum with `||c||_1 = n`, `0 < Delta <= pi`, and
`delta` uniform on `[0, Delta]^M`. For `0 < eps`,

```text
Pr[ dist(c.delta/2, (pi/4)Z) < arcsin(min(1, eps)) ] <= 4 eps (n + pi/Delta).
```

*Proof.*
* `c.delta/2` has density `<= 2/Delta` (condition on all coordinates but one with `c_x != 0`).
* It takes values in an interval of length `n Delta/2`.
* The bad set meets that interval in at most `2n Delta/pi + 2` intervals, each of length
  `2 arcsin(min(1, eps)) <= pi eps`. QED

*Proof of Theorem 5.1(a).* Let `q = 3 mod 4` be prime, `M = q + 1`, `b = 8`, `r = 8(M-1)`, and use
the canonical profile.

1. *Weights (dyadic prime cluster, as in `strict_obtuse_prime_box_countermodels.md`).*
   * Put `Z = M^6 (log M)^6`. By the prime number theorem for primes `= 1 mod 4`, `[Z, 2Z]` contains
     more than `r^2` of them for large `M`.
   * Cut `[Z, 2Z]` into `r` intervals of length `Z/r`. One of them, `[A_0, A_0 + Z/r]`, contains `r`
     such primes `p_j`.
   * Put `tau = log A_0`. Then `w_j = log p_j in [tau, tau + 1/r]` and
     `tau = 6 log M + 6 loglog M + O(1)`.
   * `W <= r tau + 1 = (48 + o(1)) M log M`. Since `p_j >= Z > 21 r^2`, the condition of
     `outside.md` Theorem A.1(ii) holds, with `Pi - 1 <= 2.2 r/sqrt(Z) -> 0`.
2. *Residue data.* Take Lemma 4.1 with `X = 2M` and `Y = 3 M log M`. For large `M` there are at least
   `M` primes `= 1 mod 4` in `(2M, Y]` (PNT in progressions). All `p_j > Y > X`, so `N` is prime to
   every `q <= X`, and no `p_j <= X` needs level data.
3. *Positions.* Take `k = 0` and `delta` uniform on `[0, Delta]^M`. By Lemma 5.3 with
   `eps = m_c^gen e^(-w/2) = (m_c^gen Delta/C) e^(-s_c/2)`, each character fails (5.1) with
   probability at most `(4pi/C)(1 + n) m_c^gen e^(-s_c/2)`. Sum over `c` modulo `+-`, using (3.2),
   Proposition 3.2 (`b = 8`), and Lemma 4.1:
   * **pairs** (`s >= 2 tau - 1/2`, at most `M^2/2` of them, `m^gen <= Y^4/4`): total
     `<= (12 pi/C) e^(1/4) M^2 Y^4 e^(-tau) = O(C^(-1) (log M)^(-2))`;
   * **ternary, support `h >= 4`** (`s >= tau h - 1/2`, at most `(2M)^h`, `m^gen <= Y^(2h)`): total
     `<= (8pi/C) e^(1/4) sum_(h>=4) (1 + h) (2M Y^2 e^(-tau/2))^h`. Here
     `2M Y^2 e^(-tau/2) <= 18/log M`, so the total is `O(C^(-1) (log M)^(-4))`;
   * **amplitude `A >= 2`** (`B - r/2 >= 3MA - 4M >= MA`, so `s >= tau M A - 1/2`; at most
     `(2A+1)^M` vectors; `n <= MA`; `m^gen <= e^(2.08M)`): total
     `<= (4pi/C) e^(1/4) sum_(A>=2) (1 + MA) exp(M(log(2A+1) + 2.08) - tau M A/2) -> 0`.

   For `M >= M_0(C)` the sum is `< 1`. So some `delta` satisfies (5.1) for every character
   simultaneously. The points are distinct, by the pair cases.
4. *Completion.* `S` has rank `M` (Lemma 3.1), so the profile is realisable for this `delta`.
   `outside.md` Theorem A.1(ii) gives a Haar-positive set of realisations `phi` satisfying (G) for
   every `v` outside `L`. For `v in L`, (G) follows from (5.1).
5. *Membership.* (a)-(d) of Definition 1.1 hold by Lemma 3.1 and steps 1-2, (G) by step 4, and (R)
   by (5.1). This proves (a) at prime Paley orders.
6. *All `M`.* For general `M`, let `M'` be the least prime Paley order `>= M`. Then
   `M' = (1 + o(1)) M` (PNT for primes `= 3 mod 4`). Lemma 5.2 restricts the witness for `M'` to `M`
   rows. Membership is monotone in `X`: truncating the residue data to moduli `<= 2M` only lowers
   every `m_c`. So `W_min^(2M)(M) <= W_min^(2M')(M') <= (48 + o(1)) M log M`. QED

*Proof of (b).* Same scheme, with Proposition 3.3 in place of 3.2, and `Z = M^2 (log M)^4`, so
`tau = (2 + o(1)) log M` and `W <= (2b + o(1)) M log M`. The class sums are:
* *pairs:* `s/2 >= (b - 4) tau/4 - 1/4 = ((b-4)/2 + o(1)) log M`, which exceeds
  `log(M^2 Y^4) = (6 + o(1)) log M` since `b > 16`;
* *columns:* `s >= tau(M - 4 + b/2) - 1/2` against `M e^(2.08M)`;
* *half-sums:* `s >= tau(M/2 - 8 + b/2) - 1/2` against `M^2 e^(2.08M)`;
* *other ternary:* `s >= (tau/2)((b/16 - 2)M + b) - 1/2`, which for `b >= 35` is `>= (3/32) tau M`,
  against at most `3^M` vectors and `m <= e^(2.08M)`. This needs `tau > 68`, true for large `M`;
* *amplitude `A >= 2`:* `s >= (tau/4)(b - 4) AM - 1/2` against `(2A+1)^M e^(2.08M)`.
All tend to 0. QED

*Proof of (c).* Repeat (a) with `X = M^A` and `Y = 2 M^A` (for `A > 1`; for `A = 1` use
`Y = 3 M log M`), so `log Y = (A + o(1)) log M`, and with `tau = (2 + 4A) log M + 6 loglog M`.
* Pairs need `tau > 2 log M + 4 log Y`.
* Ternary `h` needs `tau h/2 > h log(2M) + 2h log Y`.
* Amplitude `A' >= 2` needs `tau M A'/2 > M log(2A'+1) + 2 n log Y` with `n <= MA'`.
All hold. `W <= 8 M tau`. QED

*On (d).* Item 394 proves, for canonical `b = 5` and `tau >= 4 log M + 4 log(1000/C)`, that some
`delta` gives `min(|Re G_c|, |Im G_c|) >= 1` for all characters, with `W = 20 M log M + O_C(M)`. That
is the axis part of (G) on `L`. The diagonal part (`|Re G_c +- Im G_c| >= 1`) doubles the number of
bad components in its estimate (6), giving `0.448 < 1`. `outside.md` Theorem A.1(ii) then supplies
(G) outside `L`. Neither source states this combination.

**Corollary 5.4 (sharp barrier).** Fix `C > 0`. Let `I` be any set of inputs about a cluster that
consists of:
* (I) the exponent profile and exact exponent identities;
* (II) the Lemma 0.1 gaps for all monomials;
* (III) residue-strengthened character gaps (1.1) for **some** point-residue assignment at odd
  moduli `<= M^A`;
* (IV) any property of the angles valid off a Haar-null set of each realisation torus, or of the
  full torus.

Every actual cluster provides such inputs (Prop. 1.2). Every argument that derives `M <= f(R)` from
`I` alone has `f(R_M) >= M >= (2/K - o(1)) log R_M / loglog R_M` along a sequence `R_M -> infinity`.
Here `K = 16 + 32A`. For the A.3 range `X = 2M`, Theorem 5.1(a) gives `K = 48` and `2/K = 1/24`. In particular **no such argument
proves `M = o(log R/loglog R)`**. The best constant attainable from `I` lies in `[2/K, 2/3]`.

*Proof.* The witnesses of Theorem 5.1 satisfy (I)-(III). They satisfy (IV) because the good `delta`
form a set of positive measure, and Theorem A.1(ii) gives positive Haar proportion in every fibre (the
full-torus version is `referee_outside.md` §3.1). Finally
`log R = W/2 <= (K/2 + o(1)) M log M` gives `log R/loglog R <= (K/2 + o(1)) M`. QED

---------------------------------------------------------------------------------------------

## 6. Answers to the specific questions

**Which characters become cheap at `W = K M log M`?** Only the pairs are forced to be cheap.
* In every member, the cut identity gives `sum_(x<y) (d_xy - W/2) = MW/4 - Q <= MW/4`. So the average
  pair slack is at most `W/(2(M-1)) = (K/2 + o(1)) log M`.
* Prop. 1.3 says that small-prime pair collisions consume a share `3/K` of this average.
* Nothing else is forced cheap. In the canonical profile with `b = 8`:
  * pairs have slack `>= 2 tau - 1/2 ~ 12 log M`;
  * ternary characters of support `h` have slack `>= tau h - 1/2 >= 24 log M`;
  * characters of amplitude `A >= 2` have slack `>= tau M A - 1/2`;
  * columns and half-sums have `Theta(M log M)`.
* For an arbitrary capacity-two assignment, Lemma S gives the same picture with explicit constants
  (Prop. 3.3).

**Can windings and positions rescue all characters simultaneously?** Yes, and windings are not even
needed: take `k = 0`. A uniformly random `delta in [0, Delta]^M` avoids:
* about `M^2/2` pair slabs, each of relative width at most `M^(-2+o(1))` (`M^(-6+o(1))` without
  residue strengthening), of total measure `O((log M)^(-2))`;
* slabs of total measure `O((log M)^(-4))` for all other characters together.

This is Lemma 5.3 and the union bound of Theorem 5.1.

**Theorem F with many denominators.** Every rational character of these profiles is integral
(Lemma 3.1, twins), so the residue pigeonhole on windings has nothing to act on, for any denominator.

**Cut/Plotkin structure.** The unflipped Paley core is Plotkin-tight: every pair distance is exactly
`W/2 + Omega_0/2`, and every column is balanced. Flips perturb this by `O(tau)`. The structure forced
near `W ~ M log M` is therefore realised, not contradicted.

**Small-prime residue data.** Lemma 4.1 gives a Sidon-type design with every pair collision modulus
at most `Y^4/4 = M^(4+o(1))`. This is absorbed by the pair slack `2 tau`. It is the only reason the
constant rises from 20 (no residues) to 48.

---------------------------------------------------------------------------------------------

## 7. Relation to earlier rounds

**7.1 The barrier was already available.** Research-notes item 394
(`paley_compatible_phase_coordinate_gap.md`) constructs, for the canonical five-copy Paley profile at
`W = 20 M log M + O_C(M)`, positions such that every saturated character satisfies the two axis
coordinate gaps. It frames the remaining gap as "joint Gaussian integrality". `outside.md`
Theorem A.1 shows that in the large-prime regime this gap is empty for every Lemma-0.1-type input:
monomials outside `L` impose nothing on a Haar-positive set of realisations. The combination (with the
diagonal gaps added, Theorem 5.1(d)) is a fake circle at `O(M log M)`. Three documents missed it:
* `outside.md` §2.2 ("Nothing handles `||c||_1` in `[M^(1/2), 2M]`") and its (P1);
* `referee_outside.md` §4 ("the decisive open question is whether fake circles with
  `log R = O(M log M)` exist");
* walsh-extraction Prop. P ("individually").

This round adds:
* the residue data (Lemma 4.1, Prop. 1.2), which the question requires;
* the precise class and its validity;
* monotonicity in `M`;
* an independent route through Lemma S, which covers **every** capacity-two assignment and not only
  the canonical one. Robustness matters here, because the canonical assignment is special: its skew
  tournament argument controls the lowering flips.

**7.2 Other consistency checks.**
* *`two_thirds`.* `48 >= 3`.
* *Theorems B, B'' and C of walsh-extraction.* Their hypotheses fail: every row is flipped, and the
  core is Paley, not Walsh.
* *`linear_allocation_affine_rigidity.md`.* `M <= r + 1 = 8(M-1) + 1`.
* *`round4/flipped.md` Theorem L.* It is conditional on a Lang-Waldschmidt two-logarithm bound, and
  it excludes, for **actual** Gaussian prime angles, Hadamard-core profiles with flip mass
  `<= (1/2 - delta) W`. Our fakes have flip mass `W_F <= W/b + O(1) = W/8 + O(1)`. So the genuine
  versions of these profiles are excluded conditionally. This is consistent: transcendence is
  outside the class.
* *Walsh-extraction §10.* The "fully flipped Hadamard phase lemma" is **false for angle models**.
  Theorem 5.1 exhibits fully flipped Paley angle models at `W = O(M log M)` satisfying every
  character gap: canonical with `b = 8` in (a), and every capacity-two assignment with `b >= 35`
  in (b). So that lemma, if true for actual circles, needs an input
  outside `A_X`.

---------------------------------------------------------------------------------------------

## 8. What is not covered, and open problems

The class `A_X` has point-residue data. Actual circles have more structure: their point residues
**factor through residues of the individual primes `pi_j`**. That factorisation is also what makes
the strengthened gaps hold for every monomial, not only for characters (`outside.md` §2.4, §6 item 3).
Theorem 5.1 does not cover this refinement. The obstruction is the following.

**Proposition 8.1 (quadratic parity of factored residues).** Let `q = +-1 mod 8` be an odd prime not
dividing `N`. Then `8 | q - chi(q)`, and the norm-one group `U_q` of `Z[i]/q` has `U_q/U_q^2 = Z/2`,
with all four Gaussian units in `U_q^2`. Suppose
`rho_x = eps_x prod_j r_j^(a_xj) conj(r_j)^(e_j - a_xj)` with `Norm r_j = p_j mod q`. Then the class
of `rho_x/rho_y` in `U_q/U_q^2` is `sum_j (a_xj - a_yj) nu_j(q) mod 2`, where `nu_j(q) = 1` iff
`p_j` is a quadratic non-residue mod `q`.

*Proof.*
* `U_q` is cyclic of order `q - chi(q)`.
* `i` has order 4, so it is a square iff `8 | q - chi(q)`.
* `r/conj(r)` is a square in `U_q` iff `Norm r` is a square mod `q`.
  * *Split case.* `Z[i]/q = F_q x F_q`, `r = (alpha, beta)` with `alpha beta = Norm r`, and
    `U_q = {(u, u^(-1))} = F_q^*`. Then `r/conj(r)` corresponds to
    `alpha/beta = alpha^2/Norm r`, which is a square iff `Norm r` is.
  * *Inert case.* `r/conj(r) = r^(1-q)`, and `u -> u^(1-q)` maps `F_(q^2)^*` onto `U_q` with kernel
    `F_q^*`. This kernel consists of squares, because `q + 1` is even. So `r^(1-q)` is in `U_q^2` iff
    `r` is a square in `F_(q^2)`, iff `Norm r = r^(q+1)` is a square in `F_q`.
* Multiply over `j`. QED

The design of Lemma 4.1 extends to factored residue data if the Legendre vectors `nu(q)` are
orthogonal mod 2 to all row differences `a_x - a_y`, for every `q = +-1 mod 8` up to `X`.
* Then all `rho_x` lie in one coset of `U_q^2`.
* The twins let the square parts be chosen freely per row, for example as
  `(varpi_x / conj varpi_x)^2`.
* For `q = +-3 mod 8`, the unit `i` is a non-square and absorbs the parity in the generous collision.

On actual primes `p_j`, however, that orthogonality is a system of about
`(M - 1) #{q <= X : q = +-1 mod 8}` linear conditions on their quadratic characters. We know no way
to meet it with `p_j <= M^(O(1))`.

One sufficient way is to make every `p_j` a quadratic residue modulo every such `q`. That puts `p_j`
in prescribed classes modulo a number of size `e^((1+o(1)) M)`. Quantitative forms of Linnik's
theorem should then give enough such primes with `log p_j = O(M)`, i.e. `W = O(M^2)`, which is
already the A.3 regime. This is not written out here. Hence:

* **(P3)** Is `W_min = O(M log M)` when the residue data are required to factor through residues of
  the `pi_j`? Such data also give residue-strengthened gaps for every monomial at moduli `<= 2M`.
  Prop. 8.1 shows that factored data carry genuinely new parity information. It is tied to the
  Legendre symbols of the varying primes by quadratic reciprocity. This is the one "local" input
  whose power is now undecided. A growth improvement from it would need to turn the parity
  constraint into forced collisions. Balanced parity classes force no extra pair collisions, so
  pigeonhole alone does not do this.
* **Inputs outside every local class** (unchanged from `outside.md` §6):
  * exact algebraic identities among the points (auxiliary polynomials, capped by
    `threshold/upper.md`; Ptolemy/(L_8));
  * transcendence of the prime angles (`flipped.md` Theorem L would exclude the genuine versions of
    exactly these fake profiles, conditionally);
  * multi-place Diophantine approximation (barred).
* **The constant.** `3 <= liminf W_min^(2M)/(M log M) <= limsup <= 48`.
  * The upper constant can be lowered by sharpening Lemma S (the census suggests `theta ~ 1/3`), by
    using the exact-unit strengthening instead of the generous one, or by a finer pair count.
  * The lower constant is `two_thirds`.
  * Determining the true constant would locate the exact ceiling of the method class.

---------------------------------------------------------------------------------------------

## 9. Files (`round5/checks/`)

| script | what it checks | type |
|---|---|---|
| `paley.py` | normalised Paley-I matrix (verified Hadamard), transform, exceptional set | helper |
| `check_lemmaS.py` | Part 1: transposed uncertainty (exhaustive `q = 7, 11`; 43,735 vectors). Part 2: defect identity, integrality, `F >= n`, both peelings, small-range bound (31,990 vectors, 8 orders). Part 3: minimum of `(F-M)/M` over non-exceptional `c` (exhaustive small support, structured, hill climbing), `q = 7..107` | exact; Part 3 is evidence |
| `check_medium.py` | band-restricted hill climbing for near-flat characters with `||c||_1 in [M/4, 2.5M]` | evidence |
| `check_canonical.py` | Prop. 3.2 bounds for the canonical assignment, `q = 11..59`, `b = 5, 8` | exact |
| `check_fake.py` | identity (3.1), twins, rank, Prop. 3.3 bounds for random capacity-two assignments; Lemma 4.1 collision criterion and divisibility | exact |
| `unionbound.py` | numerical value of the union bound of Theorem 5.1(b) for explicit `b, M`, assuming the needed primes exist | illustration only |

Outputs: `check_lemmaS_out.txt` (ALL CHECKS PASSED), `check_medium_out.txt`,
`check_canonical_out.txt` (ALL CHECKS PASSED), `check_fake_out.txt` (ALL CHECKS PASSED) and
`unionbound_out.txt`. In the last, the bound is `< 1` from `M = 44` for `b = 48, 65`. For `b = 35, 40`
the asymptotic regime starts much later, as the proof predicts (`tau > 68` is needed).

Notes searched before claiming novelty (`research/docs/`):
* `paley_canonical_all_character_half_height.md` (item 392; used as Prop. 3.2);
* `paley_primefield_polynomial_character_gap.md` (item 393; exact-flat classification and the
  uncertainty argument, extended here to Lemma S);
* `paley_compatible_phase_coordinate_gap.md` (item 394; §7.1);
* `hadamard_ternary_support_spectral_gap.md`, `paley_sqrt_support_fourth_moment_gap.md`;
* `endpoint_continuation_status.md` (items 391-395);
* `round4/outside.md` and its referee report.

Not found anywhere:
* the quantitative stability (Lemma S) for all integer `c`, with its integrality/peeling argument;
* the residue design;
* the validity proposition for residue-strengthened characters (Prop. 1.2);
* the combination yielding fake circles at `O(M log M)`.

The repository `/home/user/Jarnik` was not modified.
