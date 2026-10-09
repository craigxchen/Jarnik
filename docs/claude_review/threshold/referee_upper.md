# Referee report on `tau/upper.md`: is the threshold `tau(M) = tau*(M) = 2 - 1/ceil(M/2)`?

Referee: adversarial check, 2026-10-09. All my scripts and logs are in
`scratchpad/tau/referee_upper_checks/`.

* `rlib2.py` is written from scratch. It imports nothing from `upper_checks/` or from other
  agents' code.
* It uses different primes: `2^21 - 9` and `2^21 - 19`.
* It uses a different coordinate projection: it drops `k_1`, whereas the author drops `k_M`.
* It enumerates the point sets directly from the definition. The enumeration is validated
  against an independent counting formula.

Conventions:

* **Re-derived** means I checked the argument line by line, with all quantifiers.
* **Certified** means an exact computation:
  - an upper bound on `reg` comes from full rank modulo a prime, which implies full rank over `Q`;
  - a lower bound comes from an explicit integer measure whose moments are checked in Python
    integers.
* **Evidence** means anything else.

## 0. Verdict

**The result survives in full.** I re-derived every Theorem, Lemma, Proposition and Corollary
in `upper.md`, with all quantifiers, and found no gap. That covers:

* the duality `order = regularity`;
* the corrected reduction to full slices;
* the peeling lemma and the induction over the family `l_j`;
* the alternant lower bound, the cone tiling and the asymptotic exactness;
* the closed forms for `M = 3, 4`;
* the generalised alternant;
* the geometric corollaries: bigness, `beta`, and the movable-class ratio.

I also attacked the upper bound with a rigorous lower bound that does not depend on any
construction: the generic dimension count. I recomputed regularities with independent code up
to `M = 8`. I tested the reduction (Lemma 1.2) directly on raw forms in `(u, v)`, without the
slice formalism. Nothing contradicts the claims.

All issues found are **minor** (§4):

* imprecise wording about which points attain the maximum;
* the exponent statement is a supremum, not an attained value;
* the barrier's scope should say "along `Y`";
* one check is labelled "exact" but uses floating point (harmlessly);
* a citation can replace "as I recall".

The corrected statement (§5) is the author's statement with these clarifications.

## 1. Re-derivation

### 1.1 Lemma 1.1 (`max ord = reg`)

Correct.

* (a) Substituting `w = e^y` shows that `ord_{w=1} f_c` is the least `d` with
  `sum_k c_k (k.y)^d != 0`. This is equivalent to the moment condition.
* (b) `c` annihilates `Poly_{<m}|_A` iff `Poly_{<m}|_A` is a proper subspace of `Q^A`,
  i.e. `m <= reg A`.
* Field independence holds because ranks of rational matrices are unchanged by extension.

### 1.2 Lemma 1.2 (correct reduction to full slices `A_{s,D} = Q^M((D+s)/2, (D-s)/2)`)

Correct. Some details I checked:

* **Basis and uniqueness.**
  - The torus chart `u_j = u w_j`, `v_j = v/w_j` sends `N^{(D-|k|_1)/2}[z^k]` to
    `u^{(D+s)/2} v^{(D-s)/2} w^k`. Distinct `k` give distinct monomials, so the expansion is
    unique.
  - Linear independence therefore comes from the chart, not from projective normality.
  - Projective normality (complete intersection, hence ACM, hence `H^0(O(D))` = forms) is what
    item 4 needs. See issue m1.
* **Parity.** `|k|_1 = sum k (mod 2)`, so the parity condition in `A_{s,D}` is automatic once
  `s = D (mod 2)`.
* **Chart.** `{u_1 v_1 != 0}` in `X_M` is isomorphic to a torus, with coordinates
  `(u, v, w_2, .., w_M)`; there `u_j v_j = u_1 v_1 != 0` forces `u_j != 0`. The Taylor
  coefficients along `Y = {w = 1}` are polynomials in `rho = v/u` with distinct exponents
  `(D-s)/2`. So `ord_Y F = min_s ord G_s` exactly.
* **Generic order versus real points.**
  - The Liouville step needs `F in I_Y^m` near every real point of `Y`. Real points have
    `uv = |z|^2 > 0`, and they are Zariski dense in `Y`.
  - Since `Y` is smooth and lies in the smooth locus, `I_Y^m` is unmixed. So vanishing near one
    point of `Y`, or generically, is the same as vanishing along all of `Y`.
  - Hence the generic `ord_Y` is the right quantity, and no special point of `Y` can do better
    for a form that is independent of the arc position (Remark 1.4 is right).
* **Item 4 (pseudo-effective threshold).**
  - `Y` lies in the smooth locus (`conditional.md` Lemma 2.1), so `pi_* O(-mE) = I_Y^m`.
  - `X_M` is a complete intersection, so `H^0(O(D))` consists of forms.
  - Pseudo-effective plus big is big, so the pseudo-effective threshold equals
    `sup m/D = tau*`.
  - Superadditivity holds because `X_M` is integral and `ord_Y` is a valuation.
* **Direct test, independent of the slice formalism (`t5_direct.py`, log `t5_direct.txt`).**
  - Method: take all raw monomials of degree `D` in `u_1..u_M, v_1..v_M` and Taylor-expand
    them in `x = w - 1` at `D + 1` random values of `rho` (mod `2^21 - 9`). Then compute
    `max{m : rank T_m < dim H^0(O(D))}`.
  - Cases: `M = 2` (`D <= 6`), `M = 3` (`D <= 8`), `M = 4` (`D <= 5`), `M = 5` (`D <= 3`);
    22 cases in all.
  - The direct maximum equals `max_s reg(A_{s,D})` in all 22 cases.
  - The lead's shell-only prediction (`max` over shells with `t <= D`) is strictly smaller in
    7 cases, for example `M = 3`, `D = 8`: 9 against 12. This confirms the correction in
    Remark 1.3 at the level of the maxima, not only for the single `X_2` example.
* **Remark 1.3.**
  - `z_1 zbar_2 + zbar_1 z_2 - 2|z_1|^2 = -|z_1 - z_2|^2 + (|z_2|^2 - |z_1|^2)`, which is
    `-|z_1 - z_2|^2` on `X_2`.
  - In the chart, `(u_1 - u_2)(v_1 - v_2) = -uv (1 - w_2)^2 / w_2`, which has order 2.
  - The shell parts are `2N` and `-2N` on `Y`.

  Correct. (Also found independently by `verify.md` §3.2.)

### 1.3 Lemmas 2.1–2.3 (peeling)

Correct.

* Lemma 2.1: `P = P_1 + L_1 P'` interpolates `g` on `A`, with
  `deg P <= max(reg A_1, reg A' + 1)`. Induction on `r` works with empty `A_i` too.
* Lemma 2.2:
  - the `k_1 = c` sections are `Q^{M-1}(p-c, n)` for `0 <= c <= p` and `Q^{M-1}(p, n-|c|)`
    for `-n <= c < 0`;
  - they are empty outside `[-n, p]`;
  - `reg` is invariant under the affine identification of a section with `Q^{M-1}`.
* Lemma 2.3: `N(r_(i)) >= i` gives `r_(i) + i - 1 <= K` for real `K`.

### 1.4 Theorem 3.1: `reg Q^M(p,n) <= l^M_j(p,n) = (M-1)(p/(M-j) + n/j)` for all `M >= 2`, `p, n >= 0`, `1 <= j <= M-1`

Correct, with all quantifiers.

* **Structure.** The induction is on `M`, carrying all `(p, n)` and all `j` at once.
* **Base case `M = 2`.** `p + n + 1` collinear points, so `reg = p + n = l^2_1`.
* **Case `j = 1`.** `K = p + mn` with `m = M - 1`.
  - `r^+_c <= V - c` and `r^-_c <= V`, where `V = p + (m-1)n`.
  - Sub-case `v >= (m-1)n`: `N(v) <= (V - v + 1) + n = K - v + 1`.
  - Sub-case `v < (m-1)n`: `N(v) <= p + n + 1 < K - v + 1`.
  - Only `j' = 1` is used, which is admissible for every `m >= 2`. So `M = 3` is covered.
* **Case `j = m`.** Apply case `j = 1` under `k -> -k`, using `l^M_m(p,n) = l^M_1(n,p)`. The
  case `j = 1` is proved for all `(p,n)` in the same step.
* **Case `2 <= j <= m-1`** (so `M >= 4`).
  - The indices `j' = j` and `j' = j-1` both lie in `[1, m-1]`.
  - With `alpha = (m-1)/(m-j)` and `beta = (m-1)/(j-1)`: `1/alpha + 1/beta = 1`.
  - By hand, `V+/alpha + V-/beta = p m/(m+1-j) + n m/j = l^M_j(p,n)`.
  - Each slice value is `<= min(V+, V-)`, because both `l^{M-1}_{j}` and `l^{M-1}_{j-1}` bound
    both kinds of slice and are monotone. So `x+, x- >= 0`.
  - Counting gives `N(v) <= floor(x+) + 1 + floor(x-) <= K - v + 1`.
* **Edge cases.** `p = 0`, `n = 0`, `(p,n) = (0,0)` and non-integral `K` all go through.

**Corollary 3.2** is correct.

* `M = 2k`: `l_k = ((2k-1)/k)(p+n)`.
* `M = 2k+1`: `(l_k + l_{k+1})/2 = ((2k+1)/(k+1))(p+n)`, and `k + 1 <= M - 1`.

**The crude recursion is also correct.** I re-derived
`reg <= tau'(p+n) + (2 - tau') min(p,n)` in its three sub-cases, `x <= n`, `n < x <= p` and
`x > p`. This gives `tau_M <= 1 + tau_{M-1}/2`, hence `2 - 2^{2-M}`.

### 1.5 Lower bound

* **Proposition 4.1** is correct. It holds with *exact* order `C(r,2)`, for a shorter reason
  than divisibility. For consecutive `v = (-a..b)`:
  ```text
  sum_sigma sgn(sigma) w^{sigma v}  =  ± (w_1 ... w_r)^{-a} prod_{i<j} (w_j - w_i).
  ```
  Each factor `w_j - w_i` has order exactly 1 at `w = 1`, and lowest-order parts multiply. I
  checked this identity as an exact Laurent-polynomial identity for all `(a,b)` with
  `r <= 7` (`t4_exact.py` (a)).
* **Lemma 4.2** is correct: supports add, `|k + k'|_1 <= |k|_1 + |k'|_1`, the monomial shift
  lands in `Q^M(p,n)`, and orders add.
* **Lemma 4.3** is correct: `l_j(u_j) = l_j(u_{j+1}) = C(M,2)`, and the slopes strictly
  increase. I also checked `max_j C(M,2)/|u_j| = 2 - 1/ceil(M/2)` for `M <= 200`.
* **Theorem 4.4** is correct.
  - (a) `C(M,2) <= reg S(u_j) <= reg Q(u_j) <= phi_M(u_j) <= C(M,2)`.
  - (b) The cone coefficients `lambda, lambda'` are rational. The floor loss is at most
    `2 C(M,2)`.
* **Theorem 5.1** follows: `tau = tau* = 2 - 1/ceil(M/2)`, attained. See issue m2 about the
  list of maximizers.
* **Theorem 5.3** (closed forms) is correct.
  - `M = 3`: `n` hexagons plus `p - n` linear factors.
  - `M = 4`, `p >= 3n`: blocks of cost `(3,1)` plus linear factors.
  - `M = 4`, `n <= p < 3n`: `x = floor((p-n)/2)` blocks of cost `(3,1)`, `n - x >= 0`
    hexagons and `eps` linear factors, for total order `3n + 3x + eps = floor(3(p+n)/2)`.
* **Remark 5.4** is correct.
  - The 4 orbit representatives are
    `(-3,0,1,2,3), (-2,-1,0,1,5), (-2,-1,0,2,4), (-2,-1,1,2,3)`.
  - On `e_1 = 3`, `Sym_{<=3}` is spanned by `1, e_2, e_3`.
  - Rebuilt from scratch: the condition matrix has rank 3, `gamma = (1, 1, -2, 3)`, the
    support has 480 points, and the **exact order is 14**.
  - The semigroup bound is `L_5(6,3) = 13`, as claimed.

### 1.6 Geometry and consequences

* **Corollary 5.2.**
  - Items 1–3 are correct. Big iff `t < tau*` because big classes are the interior of the
    pseudo-effective cone.
  - Item 4 is correct. Each summand is `<= 1` and the sum stops at `m = N tau*`, so
    `beta <= tau*` (in fact `beta < tau*`). The twisted version is the same argument.
  - Item 5: `1 - 1/(2 - 1/k) = (k-1)/(2k-1)`. See m4.
* **The open question in `conditional.md` §6.** There `beta_M` is literally `beta(H,E)` on
  `X^#_M`, with the same Ru–Vojta formula, and `h^0` on `X^#` equals `h^0` on `Bl`
  (normality). So "`beta(H,E) > 2` for some `M`?" is correctly answered **no**.
  - Consistency check: the rigorous lower bound `2n/(n+1)^{1+1/n}` of Prop. 6.2 stays below
    `2 - 1/ceil(M/2)`.
* **Remark 3.3.**
  - Line bound: multiplying by `prod (1 + lambda v_i)^n` gives a polynomial of degree
    `p + (M-1)n`, since every `k_i >= -n`. A general line through `e` meets the torus. So
    `ord_e <= l_1`. Correct.
  - Movable classes: if `H.Gamma = 0` then `E.Gamma <= 0`, because `H - tau* E` is effective.
    So the supremum of `E.Gamma/H.Gamma` is taken over `H.Gamma > 0` and equals `1/tau*`, by
    BDPP on the resolution `X^#`. Correct; see m6.
* **§6.1 (dual polytope).**
  - The vertices satisfy all generator inequalities, by Theorem 3.1 together with
    Proposition 4.1.
  - Conversely, the inequalities at `u_1, .., u_M` give `y.x >= phi_M(x)` on each cone `C_j`,
    which is the support-function description of `conv(vertices) + R^2_{>=0}`.

  Correct. `vertices.py` reproduces as documented.
* **§6.2 (`M = 3` del Pezzo picture).** Conics through the 4 points give `2p + n`, and lines
  through `e` give `p + 2n`. Correct.

## 2. Exact computations (independent code, `referee_upper_checks/`)

| check | what | result |
|---|---|---|
| `t1_enum.py` | `ballQ`/`shellS` enumeration vs the counting formula `sum_{a,b} C(M,a) C(M-a,b) C(p-1,a-1) C(n-1,b-1)`, `M <= 6`; sizes vs `growth/referee_checks/log_tau.txt` | PASS |
| `t2_dimcount.py` | **Adversarial.** `reg Q >= min{d : C(d+M-1, M-1) >= |Q|}` is a rigorous lower bound that uses no construction. Compared with `phi_M` on 66,285 triples (`M <= 40`, `t <= 120/60/30`). | **0 violations**; equality with `floor(phi)` in 13,288 cases |
| `t3_modreg.py` | mod-`q` `reg` of `Q` and `S` with two new primes, dropping `k_1`; `M = 3` (`t <= 12`), 4 (`t <= 9`), 5 (`t <= 7`), 6 (`t <= 6`), 7 (`t <= 5`), 8 (`t <= 4`); 260 rows | **0 rows above `floor(phi_M)`**; agrees with the author's logs in 222/222 overlapping rows; `M = 7, 8` new |
| `t4_exact.py` | exact integer moments: alternant identity (`r <= 7`); `f_j` has exact order `C(M,2)` and support in `Q^M(u_j)` (`M <= 5`); explicit products for `reg Q^5(3,2) >= 7`, `reg Q^6(4,2) >= 9`, `reg Q^3(5,5) >= 15`, `reg Q^4(6,4) >= 15`, `reg Q^4(5,4) >= 13`; generalised alternant `M = 5, (6,3)`, order 14; `L_5(6,3) = 13` | all PASS |
| `t5_direct.py` | Lemma 1.2 on raw forms (§1.2) | 22/22 agree |
| `t6_table_peel.py` | §7 table rows in my range (11/11); closed forms of Thm 5.3 vs my data (48 + 29 slices); re-implemented sorted peeling `B_M <= phi_M` (10,172 cases, `M <= 24`); crude recursion; `max_j` ratio for `M <= 200` | all PASS |
| `author_copy/` reruns | `certify.py` (A)–(E), `gen_alternant.py`, `vertices.py` | reproduce exactly (`rerun_author_*.txt`) |

**Certified exact values.** Combining `t3` (upper bounds) with `t4` (lower bounds):

* `reg Q^5(3,2) = 7 < 8 = phi_5(3,2)`;
* `reg Q^6(4,2) = 9 < 10`;
* `reg Q^5(6,3) = 14 = phi_5(6,3)`, where the upper bound is Theorem 3.1.

The author's claims are confirmed. Slices strictly below `floor(phi_M)` at `M >= 5` also occur at
`M = 6, 7, 8`, for example `(2,1)` and `(2,2)`. That is consistent with Theorem 4.4(b), which is
only asymptotic.

**The author's data claim.** "262 rows, 0 failures; every slice value certified" is correct as
stated. I checked the logs:

* 96 + 70 + 58 + 38 rows, with no `FAIL`;
* `L = reg_q` on every ball row except `M = 5, (6,3)`, which is covered by Remark 5.4;
* shell values are mod-`q` upper bounds only, as `upper.md` says.

## 3. Hidden assumptions and quantifier checks: all fine

* The order is uniform in the arc position because `I_Y^m` is unmixed (§1.2). The barrier is
  genuinely about forms vanishing **along** `Y`; see m5.
* The coefficient field does not matter: `Z[i]`, `Q(i)` and `C` give the same ranks.
* Forms that are inhomogeneous or `N`-weighted (several shells, polynomial in `N`) are already
  forms of one degree `D` on `X_M`. They are covered by `tau*`, as `upper.md` §8.1 says.
* `tau` (shells) `<= tau*` (slices) by monotonicity of `reg`, and both are attained on the
  shell `S(u_j)`.
* The induction in Theorem 3.1 never uses an index outside `[1, m-1]`.
* The counting in the middle case uses only `x+, x- >= 0`, which is guaranteed.

## 4. Issues (all minor; none affects any stated value)

* **m1 (justification).** Lemma 1.2.1 attributes uniqueness of the expansion to projective
  normality.
  - Uniqueness follows from the torus chart (§1.2).
  - Projective normality (complete intersection, hence ACM) is what item 4 needs:
    `H^0(O(D))` consists of forms.
  - Cosmetic.
* **m2 (maximizers).** Theorem 5.1 says "the maximum is attained at the following points" and
  lists `(T_k, T_k)` (`M = 2k+1`) and `(T_k, T_{k-1})`, `(T_{k-1}, T_k)` (`M = 2k`). The list
  is not exhaustive.
  - For `M = 2k`, the ratio `2 - 1/k` is attained on the whole flat cone `C_k` wherever
    `reg = phi`. An example is the alternant on only `2k - 1` coordinates with
    `(a,b) = (k-1, k-1)`: cost `(T_{k-1}, T_{k-1})`, order `C(2k-1, 2)`.
  - Concretely, `M = 4` at `(1,1)` (`t = 2`), and `M = 6` at `(3,3)` (`t = 6`); the latter is
    in the author's own table.
  - This also explains `tau*(2k) = tau*(2k-1)`: the extremal measure for even `M` needs only
    `M - 1` coordinates.
  - Suggest "attained, for example, at".
* **m3 (labelling of a check).** `certify.py` claims "all integer / Fraction arithmetic; no
  floating point decisions", but check (D) compares complex floats with tolerance `1e-6`.
  - It is harmless: the inputs are integers `<= 99`, exact in binary64, and the identity is
    trivial algebra.
  - Check (D) does not test the order statements of Remark 1.3. Those are one-line hand
    computations, confirmed in §1.2 above and, at the level of maxima, by `t5_direct.py`.
* **m4 (exponent wording and the citation).**
  - "The best exponent ... is `(k-1)/(2k-1)`" is a supremum. The Liouville condition
    `m(1 - alpha) > D` holds exactly for `alpha < (k-1)/(2k-1)`, which is never attained.
  - The Cilleruelo–Córdoba match can be confirmed without "as I recall".
    `research/docs/integer_radius_axis_lift_audit.md` (line 88, checked against the primary
    PDF in that note) gives: at most `m` points on arcs shorter than
    `sqrt(2) R^{1/2 - 1/(4 floor(m/2) + 2)}`.
  - With `m = M - 1` and `k = ceil(M/2) = floor((M-1)/2) + 1`, this exponent is
    `1/2 - 1/(4k - 2) = (k-1)/(2k-1)`. It matches, and CC attain the endpoint itself.
* **m5 (scope).** The barrier bounds `ord_Y`, i.e. forms that vanish along all of `Y`. This is
  exactly what a form independent of the arc position needs (§1.2). It does not bound
  point-vanishing orders at a cluster-adapted point of `Y`.
  - Example: on `X_2 = P^1 x P^1`, the tangent hyperplane at a point gives point order `2D`,
    while `tau*(2) = 1`.
  - `upper.md` §8.2 acknowledges this in its last bullet. But the sentence "requirement
    becomes `m > (2 + 2 eta) D`, which is still impossible by Theorem 5.1" applies only to
    forms that still vanish along `Y`.
  - The headline ("the reduction of `tau_problem.md` can never be fed") is correct as stated.
* **m6 (Remark 3.3).** The supremum of `E.Gamma/H.Gamma` should be taken over movable `Gamma`
  with `H.Gamma > 0`. This is automatic: `H.Gamma = 0` forces `E.Gamma <= 0`. Cosmetic.
* **m7 (novelty).** Relative to the project notes, I searched `research/docs` (1170 files)
  and `growth/` for: Cremona, regularity index, peeling, interpolation, Hilbert function,
  pseudo-effective, Seshadri, `2 - 1/ceil`, alternant, root polytope, Córdoba.
  - The peeling bound `phi_M`, the exact value of `tau`, and the negative answer for `beta`
    are new.
  - The full-slice formulation and the pseudo-effective threshold predate this note, in
    `growth/conditional.md` Prop. 6.2 and `growth/referee_checks/tau_threshold.py`.
    `upper.md` says so.
  - `hermitian_polynomial_barrier.md` is a different barrier, about divisibility orders, and
    does not anticipate this one.
  - No outside literature search was done, by the author or by me. The interpolation degree of
    the lattice points of `p Delta - n Delta` may be known in the literature on Hilbert
    functions of lattice-point sets or on vanishing at a general point of toric varieties.
  - The concurrent `referee_verify.md` §6 also re-derived Theorem 3.1 and found no gap. My
    check is independent and agrees.

## 5. Corrected statement (the original survives; wording clarified)

For every `M >= 2` and all integers `p, n >= 0`:

```text
reg Q^M(p,n)  <=  phi_M(p,n) = (M-1) min_{1<=j<=M-1} ( p/(M-j) + n/j ),
lim_k reg Q^M(kp,kn)/k = phi_M(p,n).
```

Equality `reg = C(M,2) = phi_M` holds at the `M` points `u_j = (T_{M-j}, T_{j-1})`. Hence:

* every form `F` of degree `D` on `P^{2M-1}` with `F|X_M != 0` satisfies
  ```text
  ord_Y (F|X_M)  <=  (2 - 1/ceil(M/2)) D  <  2D;
  ```
* `tau_shell(M) = tau*(M) = 2 - 1/ceil(M/2)`. This is the pseudo-effective threshold of
  `pi^*H - tE` on `Bl_Y X_M`, and it is attained by alternating Vandermonde orbits, for example
  at `u_{ceil(M/2)}`. For even `M` it is also attained by the alternant on `M - 1` coordinates.
* The hypothesis `m > 2D` of the `tau_problem.md` reduction is unsatisfiable for every `M`.
* `beta_M = beta(H,E) <= tau*(M) < 2`, which answers the question in `conditional.md` §6
  negatively.
* The single-place Liouville method on `X_M`, with forms vanishing along `Y`, works exactly for
  `alpha < (k-1)/(2k-1)`, `k = ceil(M/2)`: a supremum, not attained. This is the
  Cilleruelo–Córdoba exponent `1/2 - 1/(4 floor((M-1)/2) + 2)`.
* Forms adapted to the arc position that vanish only at a point of `Y` are outside the scope
  of this barrier.

## 6. Reproducibility

Run each script from inside `scratchpad/tau/referee_upper_checks/` (they import `rlib2.py`
from there). Logs are the `*.txt` files beside the scripts.

* `t1_enum.py`: a few seconds.
* `t2_dimcount.py`: about 14 s.
* `t3_modreg.py M tmax [maxsize]`: about 15 min for all runs.
* `t4_exact.py`: about 4 s.
* `t5_direct.py`: about 26 s.
* `t6_table_peel.py`: under 1 s.
* `author_copy/`: unmodified copies of the author's scripts; `rerun_author_*.txt` are the
  reruns.
