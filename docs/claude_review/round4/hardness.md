# Why the uniform `C sqrt(R)` bound is hard: what it implies, what it does not, and its split twin

Round 4, hardness task. Target: an implication from the uniform theorem (or from a growth
improvement `M = o(log R/loglog R)`) to a recognised open problem, or an equivalence with a named
conjecture for a specific variety.

Conventions. *Theorem*, *Proposition*, *Lemma* and *Corollary* mean proved here, with all
quantifiers, or proved in a refereed file that is cited by name. *Data* means an exact finite
computation, never used as proof. *Remark* marks interpretation. Scripts and logs are in
`round4/hardness_checks/` (Section 8).

---------------------------------------------------------------------------------------------

## 0. Verdict

**No implication from the uniform theorem, or from a growth improvement, to an independently
recognised open problem was found.** For each of the five suggested routes, the reason the standard
mechanism fails is proved below. Section 7 gives the resulting explanation of the difficulty. The
status is therefore a set of barrier theorems, not a hardness reduction.

Notation used in this summary: `U` is the uniform theorem (`M(C) < infinity` for every `C`).
`U_abs` is its absolute form (one `M_0` for all `C`, once `R >= R_0(C)`). This is the form that
Vojta's conjecture gives (`growth/conditional.md`, Thm 5.1). `A_M(kappa)` is the set of rational
points of `X_M`, off the big diagonals, whose `M` points lie within `kappa R^(1/2)` of each other.

1. **Fixed-variety form (Theorem 1).**
   * `U` holds iff for every `kappa` there is an `M` with `A_M(kappa)` **empty**. Equivalently,
     `A_M(kappa)` is finite, or is not Zariski dense.
   * `U_abs` holds iff for **one fixed** `M` the set `A_M(kappa)` is **finite for every**
     `kappa`.

   So `U_abs` is a Roth-type finiteness theorem at the Dirichlet exponent `1/2` for the rational
   line `Y` in the fixed toric variety `X_M`. The known families force `M >= 9`. The only known
   implicant is Vojta's Main Conjecture (`D = 0`) for a model of `Bl_Y X_7` or for
   `Mbar_{0,8}`. This is the requested "equivalence with a Vojta-type statement on a fixed
   variety". The statement it is equivalent to is not itself a named conjecture.
2. **The additive-combinatorial shadow is elementary (Theorem 2).**
   * For squares of integers from two intervals of length `<= c sqrt X`, the additive energy
     satisfies `E <= (2 + c^2)|A||B|` unconditionally.
   * Hence the sumset, uniform `Lambda(4)`, and Rudin-type consequences are already theorems. One
     example is `<= sqrt(6K)` squares in `{a + qn : 0 <= n < K}` whenever
     `a >= (qK/2)^(4/3)`.
   * `U` is equivalent to the `L^infinity` (bounded sum-multiplicity) statement for **two**
     square-root intervals. For **one** interval that statement is trivial.

   So no Rudin-type or sumset statement can witness the hardness of `U`. This closes the
   Cilleruelo–Granville link to Rudin's conjecture as a hardness route at exponent `1/2`.
3. **Single linear forms carry no information (Theorem 3).**
   * For two points `z, w` of a circle, the angle `lambda_zw` is a linear form in logarithms of
     Gaussian primes. It satisfies the elementary bound `|lambda_zw| >= sqrt2 n^(-1/2)`, where `n`
     is the pair norm.
   * This bound is **optimal**: it is attained up to the factor `1 + O(1/n)` by the pairs
     `(g u, g i ubar)` with `u = A + (A-1)i`. In the record 10-point cluster it is attained to
     three digits by the closest pairs.
   * Matveev's bound is weaker than the elementary bound on every such form with exponents
     `<= 10^7`. Any Lang–Waldschmidt-shaped bound is weaker on squarefree pair norms.
   * Product (cube) clusters with `>= 3` independent blocks force `N <= (C^2/2)^6`.

   So `U` has no single-form (Baker, Lang–Waldschmidt) content. Its content is the joint
   behaviour of entangled systems of Liouville-extremal forms.
4. **The split twin (Theorem 4).**
   * The divisor problems of Erdős–Rosenfeld (Erdős #887) and Ruzsa (#886) live on the **split**
     `Q`-form `X_M^split = {x_1 y_1 = ... = x_M y_M}` of the same toric variety.
   * The tori admit **no** nonzero `Q`-homomorphism in either direction (Prop. 4.1), so no
     equivariant transfer exists.
   * The `C`-dependent split analogue of `U` is trivial at the vertex `d = sqrt n` (`<= 1 + C^2`).
   * It is a **theorem away from the vertex**, by Coppersmith–Howgrave-Graham: at most
     `max(C + 1, 6/eta + 4)` divisors in the window when `D <= n^(1/2-eta)` and
     `n >= n_0(C, eta)` (Cor. 4.5). The margin is exactly `(1 - beta)(1/2 - beta)` at divisor
     size `n^beta`.
   * It follows from Ruzsa's conjecture in windows within `n^(1/2 - eps)` of the vertex (Cor. 4.7).
   * The circle problem is the case where this margin is **identically zero**: complex conjugation
     forces `|z| = |zbar|` at every point. `threshold/construct.md` proves the matching statement
     for all forced-divisibility arguments on the circle (exactly critical).
5. **Residue classes (Prop. 6.1).** For every modulus in or above the Gaussian Lenstra range, each
   residue class contains at most one point of a cluster. So divisor-in-residue-class bounds
   (Lenstra, Coppersmith–Howgrave-Graham–Nagaraj, Hales) at or above the critical exponent `1/4`
   cannot see clusters.

**Explanation in one paragraph.** `U` is a finiteness statement at the exact Dirichlet threshold:
rational points approximating a rational line to exponent `1/2`. It is posed on a toric variety
whose torus is anisotropic, so every point is "balanced", and every standard input below Vojta's
conjecture provably loses there:

* additive energy (elementary at this scale);
* single linear forms in logarithms (the trivial bound is already optimal);
* residue classes (one point per class);
* forced divisibility and Coppersmith (zero margin at balanced points);
* single-place auxiliary polynomials (`tau* < 2`, `threshold/upper.md`);
* Subspace and Ru–Vojta (`beta_M < 2`, or exceptional sets that depend on `S`);
* `abc` (vacuous on squarefree blocks, `growth/conditional.md` §8).

The recognised neighbours #887 and #886 are vertex statements on the split twin. The split twin's
own analogue of `U` is a theorem except in a thin near-balanced range. That range is exactly the
part that corresponds to the whole circle problem.

| candidate route | verdict | where |
|---|---|---|
| Rudin / squares in APs / sumsets / `Lambda(4)` | every consequence via additive energy is already unconditional at exponent `1/2`; `U` is an `L^infinity` statement | Thm 2 |
| divisors in residue classes at exponent `1/4` (Lenstra, Coppersmith) | non-archimedean bounds see `<= 1` cluster point per class; the archimedean (Coppersmith) mechanism has margin `(1-beta)(1/2-beta)`, which is `0` on the circle | Prop 6.1, Thm 4 |
| quantitative linear independence of Gaussian prime arguments | single forms: the elementary bound is optimal, so nothing beyond Liouville can be implied; product systems are Liouville-bounded | Thm 3 |
| Erdős–Rosenfeld (#887), Ruzsa (#886) | Galois twists (split form, at the vertex); no `Q`-homomorphism of tori; both follow from Vojta, neither follows from `U` by any equivariant map | Thm 4 |
| Vojta-type statement on a fixed variety | `U_abs` iff `A_M(kappa)` is finite for all `kappa` (fixed `M`); implied by Vojta on `X^#_7`; not a named conjecture | Thm 1 |

---------------------------------------------------------------------------------------------

## 1. Setting

For an integer `N >= 1` put `R = sqrt N`, `Gamma_N = {z in Z[i] : |z|^2 = N}`.
`M(C) = sup #(Gamma_N cap arc)`, the supremum over all `N` and all arcs of length `<= C sqrt R`.

* `U`: `M(C) < infinity` for every `C > 0`.
* `U_abs`: there is `M_0` such that for every `C` there is `R_0(C)` with the following property:
  every arc of length `<= C sqrt R` on a circle with `R >= R_0(C)` contains `<= M_0` lattice
  points.

`X_M = {x_1^2+y_1^2 = ... = x_M^2+y_M^2} in P^(2M-1)`, `Y` is the diagonal line, and
`Delta = union_{j<k} {z_j = z_k}` is the union of the big diagonals.

For `P in X_M(Q)` take the primitive integral representative `(x_1, y_1, ..., x_M, y_M)`, with
gcd 1 and unique up to sign. Put `z_j = x_j + i y_j`.

* Then `|z_j|^2 = N_P` is common to all `j`, and `N_P >= 1`, because over `Q`, `x^2 + y^2 = 0`
  forces `x = y = 0`.
* Put `R_P = sqrt(N_P)`.
* For `kappa > 0` define

  ```text
  A_M(kappa) = { P in X_M(Q) \ Delta : |z_j - z_k| <= kappa R_P^(1/2) for all j,k }.
  ```

Since `R_P/sqrt2 <= max|coord| <= R_P`, membership in `A_M(kappa)` is the same, up to changing
`kappa` by bounded factors, as `lambda_Y(P) >= h(P)/2 - log kappa - O(1)`. That is the
approximation condition of `growth/conditional.md` §4. `r_2(N) = |Gamma_N| <= 4 d(N) <= 4N`.

---------------------------------------------------------------------------------------------

## 2. Theorem 1: the uniform theorem as a finiteness statement on a fixed variety

**Lemma 1.1 (chords to arcs).** Let `S` be a set of points on a circle of radius `R` with pairwise
distances `<= l`, where `l < sqrt3 R`. Then `S` lies on an arc of length `<= (pi/2) l`.

*Proof.*

1. Pairwise angular distances are `<= phi := 2 arcsin(l/(2R)) < 2 pi/3`.
2. Fix `s in S`. Measure angles `t in [-phi, phi]` from `s`, and let `a, b` be the points of `S`
   with extreme `t`.
3. Suppose `t_b - t_a > phi`. Since `t_b - t_a <= 2 phi`, we also have
   `2pi - (t_b - t_a) >= 2pi - 2phi > phi`. Then the angular distance of `a` and `b` exceeds
   `phi`, which is impossible.
4. So `S` lies on an arc of angle `<= phi`. Its length is
   `R phi = 2R arcsin(l/2R) <= (pi/2) l`. ∎

**Lemma 1.2 (Bezout grid lemma; `growth/conditional.md` Lemma 3.1 and Prop. 3.2(i), refereed).**
Let `S` be a subset of `Gamma_N`. Suppose every ordered `M`-tuple of distinct points of `S`, read as
a point of `P^(2M-1)`, lies in a proper Zariski-closed `W` of `X_M`. Let `F` be a form of degree `D`
vanishing on `W` but not on `X_M`. Then `|S| <= 2D + M - 1`.

**Theorem 1.**

(a) The following are equivalent:

* (i) `U`;
* (ii) for every `kappa > 0` there is an `M` with `A_M(kappa)` empty;
* (iii) for every `kappa > 0` there is an `M` with `A_M(kappa)` finite;
* (iv) for every `kappa > 0` there is an `M` with `A_M(kappa)` not Zariski dense in `X_M`.

(b) `U_abs` holds if and only if (v): there is an `M` such that `A_M(kappa)` is finite for every
`kappa > 0`.

*Proof of (a).*

* **(i) ⇒ (ii).** Fix `kappa`. Put `C = pi kappa/2` and `M = 1 + max(M(C), floor(4 kappa^4/9))`.
  Let `P` be in `A_M(kappa)`.
  - If `R_P > kappa^2/3`, then `kappa sqrt(R_P) < sqrt3 R_P`. By Lemma 1.1 the `M` distinct points
    `z_j` lie on an arc of length `<= C sqrt(R_P)`, so `M <= M(C)`. That is a contradiction.
  - If `R_P <= kappa^2/3`, the `z_j` are `M` distinct elements of `Gamma_{N_P}`, and
    `|Gamma_{N_P}| <= 4N_P <= 4 kappa^4/9 < M`. Again a contradiction.
* **(ii) ⇒ (iii) ⇒ (iv).** Immediate, since `X_M` is irreducible of positive dimension.
* **(iv) ⇒ (i).** Fix `C`. Take `kappa = C` and `M` as in (iv), and let `W` be the Zariski
  closure of `A_M(C)`.
  - Let `z_1, ..., z_M` be distinct lattice points on an arc of length `<= C sqrt R`. Their
    chords are `<= C sqrt R`.
  - Let `d` be the gcd of their coordinates and write `z_j = d z'_j`, `R = d R'`. Then
    `|z'_j - z'_k| <= C sqrt(R)/d = C sqrt(R'/d) <= C sqrt(R')`.
  - So `[z'] in A_M(C) subset W`.
  - Lemma 1.2 bounds the arc by `2 deg F_W + M - 1`, independently of `N`. ∎

*Proof of (b).*

* **`U_abs` ⇒ (v).** Put `M = M_0 + 1`. Fix `kappa`, and put `C = pi kappa/2`.
  - Suppose `P` is in `A_M(kappa)` with `R_P > max(R_0(C), kappa^2/3)`.
  - By Lemma 1.1 its `M` points lie on an arc of length `<= C sqrt(R_P)` on a circle with
    `R_P >= R_0(C)`. This contradicts `U_abs`.
  - So every `P` in `A_M(kappa)` has bounded `R_P`, hence bounded coordinates. There are finitely
    many such `P`.
* **(v) ⇒ `U_abs`.** Put `M_0 = M - 1`. Fix `C`, and let `F = A_M(C)`, which is finite. Put
  `R_0(C) = 1 + C^2 max_{P in F} N_P` (or `1` if `F` is empty).
  - Suppose an arc of length `<= C sqrt R`, with `R >= R_0(C)`, carries `M` distinct points.
  - As in (iv) ⇒ (i), `P = [z']` lies in `F`, with `R = d R_P`.
  - Since the `z'_j` are distinct, `max|z'_j - z'_k| >= 1`. Hence
    `d <= max|z_j - z_k| <= C sqrt(d R_P)`, so `d <= C^2 R_P`.
  - Then `R = d R_P <= C^2 N_P < R_0(C)`. That is a contradiction. ∎

**Remarks.**

1. *(Which `M` can work in (v).)* The known infinite families give lower limits.
   - By `growth/families.md`, `M(C) >= 5` for every `C > 0`. This comes from the cyclotomic
     family, whose constants tend to `0` as `R` tends to infinity. Hence `A_5(kappa)` is infinite
     for every `kappa`.
   - The golden eight-point family has arc constant tending to `5.8977`. So `A_8(kappa)` is
     infinite for `kappa >= 5.9`.
   - Hence (v) can only hold with `M >= 9`.
   - Under Vojta's conjecture, `growth/conditional.md` Thm 5.1 gives `U_abs`, and with it (v) for
     `M = M_0 + 1`.
2. *(Interpretation.)* Statement (v) says that the rational line `Y` is "badly approximable at the
   Dirichlet exponent `1/2`, for every constant", in the fixed variety `X_M` off `Delta`.
   - For a real number, finiteness of `{p/q : |alpha - p/q| < c/q^2}` can hold only for small `c`,
     by Dirichlet.
   - Here it is required for every `c`. This is possible only because `dim X_M = M` is large. For
     `M >= 7` the expected number of exponent-`1/2` approximations converges
     (`growth/conditional.md` §5, heuristic).
   - It is a Roth-type statement at the exact Dirichlet threshold.
   - The only known implicants are Vojta's Main Conjecture on a model of `Bl_Y X_7`
     (`conditional.md`) or on `Mbar_{0,8}` (`ptolemy.md`). Subspace and Ru–Vojta with
     `S = {infinity}` reach only `beta_M < 2` (`threshold/upper.md`), that is, only exponents
     above `1/2`.
3. *(Not a named conjecture.)* The finiteness statement (v) is not, to my knowledge, a special case
   of any named conjecture weaker than Vojta's. This is a statement about the absence of a
   reduction, not a theorem.

---------------------------------------------------------------------------------------------

## 3. Theorem 2: the additive (`L^2`) shadow of `U` is elementary

Cilleruelo–Granville tie arcs on circles to Rudin's conjecture through sums of squares
(`Lambda(p)` sets and sumsets of squares). This section shows that every link of that type is
already a theorem at exponent `1/2`. So such links cannot carry hardness.

**Lemma 2.1 (difference multiplicity).** Let `X >= 1`, `c > 0` and `I = [X, X + L] cap Z` with
`L <= c sqrt X`. For every integer `d != 0`:

```text
r_I(d) := #{(x, x') in I^2 : x^2 - x'^2 = d} <= 1 + floor(c^2).
```

*Proof.*

1. By symmetry take `d > 0`. Then `x > x'`; put `q = x - x' in [1, L]` and
   `s = x + x' in [2X, 2X + 2L]`. Then `d = qs`, so `d <= L(2X + 2L)`.
2. The value `q` lies in `[d/(2X + 2L), d/(2X)]`. This interval has length
   `2dL/(2X(2X + 2L)) <= L^2/X <= c^2`.
3. So `q` takes at most `1 + floor(c^2)` values, and `q` determines `s = d/q` and hence `(x, x')`. ∎

(The literature check `endpoint_recent_literature_check.md` already records this bound with
`1 + floor(k^2/N)`.)

**Theorem 2.** Let `I = [X, X + L_I]` and `J = [Y, Y + L_J]` with `L_I <= c sqrt X` and
`L_J <= c sqrt Y`. Let `A` be a subset of `{x^2 : x in I}` and `B` a subset of `{y^2 : y in J}`.
Then:

1. `E(A,B) := #{a + b = a' + b'} <= |A||B| + (1 + floor(c^2)) min(|A|,|B|)^2 <= (2 + c^2)|A||B|`;
2. `|A + B| >= |A||B|/(2 + c^2)`;
3. (uniform `Lambda(4)`) for every finitely supported `(a_x)_{x in I}`:

   ```text
   || sum_{x in I} a_x e(x^2 theta) ||_{L^4(T)}^4 <= (2 + c^2) ||a||_2^4.
   ```

*Proof.*

1. Write `E(A,B) = sum_delta r_{A-A}(delta) r_{B-B}(-delta)`.
   - The term `delta = 0` contributes `|A||B|`.
   - For `delta != 0`, Lemma 2.1 gives `r_{B-B}(-delta) <= 1 + floor(c^2)`. So the remaining sum
     is `<= (1 + floor(c^2))(|A|^2 - |A|)`.
   - Exchange the roles of `A` and `B` to get the `min`.
2. Apply Cauchy–Schwarz: `|A+B| >= |A|^2|B|^2/E(A,B)`.
3. The fourth power of the `L^4` norm equals
   `sum_m |sum_{x^2 - x'^2 = m} a_x conj(a_x')|^2`.
   - The term `m = 0` is `||a||_2^4`.
   - For `m != 0`, Cauchy–Schwarz with Lemma 2.1 bounds each term by
     `(1 + floor(c^2)) sum_{x^2 - x'^2 = m} |a_x|^2 |a_x'|^2`.
   - The sum over `m` is then `<= (1 + c^2)||a||_2^4`. ∎

**Corollary 2.2 (Rudin's bound in the far regime, unconditional).** Let `a, q, K` be positive
integers with `a >= (qK/2)^(4/3)`. Then `{a + qn : 0 <= n < K}` contains at most `sqrt(6K)`
squares.

*Proof.*

1. The squares are `y^2` with `y` in `I = [ceil(sqrt a), floor(sqrt(a + q(K-1)))]`. This interval
   has length `<= qK/(2 sqrt a) <= a^(1/4) <= sqrt(min I)`, so `c = 1`.
2. Theorem 2 gives `|A + A| >= |A|^2/3`.
3. Also `A + A` is a subset of `{2a + qm : 0 <= m <= 2K - 2}`.
4. Hence `|A|^2 <= 3(2K - 1)`. ∎

(Gabdullin, JFAA 30 (2024), Thm 1.4 (cited in `endpoint_recent_literature_check.md`), gives
uniform `Lambda(4)` on the much longer intervals `[N, N + N^(0.618 - eps)]`. So Corollary 2.2 is
neither new nor sharp. It is included because it is exactly the implication that
"lattice points on arcs ⇒ Rudin" would deliver at exponent `1/2`.)

**Proposition 2.3 (`U` is an `L^infinity` statement, and only its two-interval case is hard).**

1. `U` is equivalent to the following. For every `c > 0` there is `g(c)` such that for all
   `n >= 1` and all intervals `I, J`, each a subset of `[0, infinity)` of length `<= c n^(1/4)`:

   ```text
   #{(x, y) in (I x J) cap Z^2 : x^2 + y^2 = n} <= g(c).
   ```
2. If `I` and `J` have the same left end, or more generally if `|x - y| <= c n^(1/4)` on
   `I x J`, the count is `<= 2(c^2/sqrt2 + 2)` unconditionally.

*Proof.*

1. **`U` ⇒ the box statement.** The points lie in a box of side `c n^(1/4) = c sqrt R`, so their
   pairwise distances are `<= sqrt2 c sqrt R`.
   - If `R > 2c^2/3`, Lemma 1.1 puts them on an arc of length `<= (pi/sqrt2) c sqrt R`, so there
     are at most `M(pi c/sqrt2)` of them.
   - Otherwise there are at most `4n <= 16 c^4/9` of them.
2. **The box statement ⇒ `U`.** Take an arc of length `<= C sqrt R` with `R > (8C/pi)^2`.
   - Its angular length is `< pi/8`, so it meets at most two closed quadrants.
   - Map each piece into the first quadrant by a symmetry of `Z^2`.
   - In each piece both coordinates range over intervals of length `<= C sqrt R = C n^(1/4)`.
   - So the arc holds at most `2g(C)` points. Arcs with `R <= (8C/pi)^2` hold at most
     `4N <= 4(8C/pi)^4` points.
3. **Item 2.** Put `zeta = z(1 - i) = (x + y) + i(y - x)`. It lies on the circle of radius
   `sqrt(2n)`, with `|Im zeta| <= c n^(1/4)`.
   - Then `sqrt(2n) - (x + y) = (y - x)^2/(sqrt(2n) + x + y) <= c^2/sqrt2`.
   - So `x + y` takes at most `c^2/sqrt2 + 1` integer values.
   - Each value of `x + y` gives at most two points `(x, y)`. ∎

**What this shows.**

1. Rudin's conjecture concerns squares from one progression, so through `A + A` it sees only the
   one-interval case. That case is trivial in `L^infinity` (Prop. 2.3(2)) and elementary in `L^2`
   (Theorem 2).
2. Two progressions with the same difference, both in the far regime, give
   `|A||B| <= (2 + c^2)(K_1 + K_2)` from Theorem 2 already.
3. `U` would improve sumset lower bounds only by replacing the constant `2 + c^2` with `g(c)`.
4. So no statement about squares in progressions or sumsets of squares that follows from `U` through
   representation counts is open at exponent `1/2`.

The content of `U` is the gap between `L^2` and `L^infinity` for **two separated** intervals.
`endpoint_recent_literature_check.md` gives an abstract (non-square) example showing that this gap
is real: bounded difference multiplicity and bounded energy, yet unbounded sum multiplicity.

---------------------------------------------------------------------------------------------

## 4. Theorem 3: single linear forms carry no information; product systems are elementary

**Lemma 3.1 (pair identity; `growth/conditional.md` Prop. 8.1, refereed).**

*Statement.* Let `z != w` be in `Gamma_N`. Put `g = gcd(z, w)`, `u = z/g` and
`n = N/|g|^2 = |u|^2`. Then:

* `w = g eta ubar` for a unit `eta`, and `gcd(u, ubar) = 1`;
* every prime `p | n` splits, so `p = 1 (mod 4)`;
* `|z - w| = |g||u - eta ubar| >= sqrt2 |g|`.

*Proof.*

1. **The form of `w`.**
   - Write `w = g v` with `gcd(u, v) = 1`.
   - Since `u ubar = v vbar` and `u` is coprime to `v`, `u` divides `vbar`.
   - The norms are equal, so `vbar = eta' u` for a unit `eta'`, that is, `v = eta ubar`.
   - Then `gcd(u, ubar) = gcd(u, v) = 1`.
2. **Every `p | n` splits.**
   - `(1+i)` cannot divide `u`, because `1 + i` and its conjugate are associates.
   - An inert prime `q` cannot divide `u`, since `q = qbar` would divide `gcd(u, ubar)`.
3. **The bound `|u - eta ubar| >= sqrt2`.** Write `u = A + Bi`. The value of `u - eta ubar` is:
   - `2Bi` for `eta = 1`;
   - `2A` for `eta = -1`;
   - `(A - B)(1 - i)` for `eta = i`;
   - `(A + B)(1 + i)` for `eta = -i`.

   It is nonzero because `z != w`, and each nonzero value has modulus `>= sqrt2`. ∎

Let `lambda_zw in (-pi, pi]` be the angle from `w` to `z`. Then `z/w = eta^(-1) u/ubar`. Write
`u = eps prod_{p | n} pi_p^(a_p)`, with `pi_p` the prime factor of `u` above `p`. Then

```text
i lambda_zw = sum_{p | n} a_p log(pi_p/pibar_p) + b log i      (mod 2 pi i),
```

which is a linear form in `omega(n) + 1` logarithms. The numbers `alpha_p = pi_p/pibar_p` are
roots of `p X^2 - 2 Re(pi_p^2) X + p`, whose roots have modulus `1`. So their Mahler measure is
`M(alpha_p) = p`, and their absolute logarithmic height is `h(alpha_p) = (1/2) log p`.

**Proposition 3.2 (the elementary bound and its optimality).**

1. `|lambda_zw| >= |z - w|/R >= sqrt2/sqrt(n)` for every pair.
2. If `z, w` lie on an arc of length `<= C sqrt R`, then `|lambda_zw| <= C R^(-1/2)`. Hence
   `n >= 2R/C^2`, and the elementary bound is attained up to the factor `C sqrt(n/(2R))`.
3. **Optimality.** Let `A >= 2`, `u = A + (A-1)i` (so `gcd(u, ubar) = 1`) and
   `n = A^2 + (A-1)^2`. For any nonzero `g`, the pair `z = g u`, `w = g i ubar` lies on the
   circle of norm `N = |g|^2 n`, has pair norm `n`, and satisfies

   ```text
   |lambda_zw| = 2 arcsin(1/sqrt(2n)) = sqrt2 n^(-1/2) (1 + O(1/n)).
   ```

   * So every true lower bound for `|lambda_zw|` in terms of `n` alone is at most
     `(1 + O(1/n))` times item 1.
   * Item 1 holds unconditionally.
   * (The family `u = A + i`, `w = g ubar`, attains it up to the factor `sqrt2`.)
4. **Matveev.** Matveev's theorem (complex case, `D = 2`; Bugeaud's normalisation, as in
   `growth/auxiliary.md` Prop. B2) gives
   `log|lambda| > -c(k) D^2 (1 + log D)(1 + log B) prod_j A_j`.
   * Here `k = omega(n) + 1` and `c(k) = min((1/2)(ek/2)^2 30^(k+3) k^3.5, 2^(6k+20))`, so
     `c(k) >= 7*10^5` for all `k >= 1`.
   * The height parameters are `A_p = max(2 h(alpha_p), |Log alpha_p|, 0.16) >= log p >= log 5`
     and `A_{log i} = pi/2`.
   * For numbers `x_j >= log 5 > 1.6`, `prod x_j >= 0.8 sum x_j`. So the exponent is at least
     `7*10^5 * 4(1 + log 2)(pi/2)(0.8) sum_p log p > 5.9*10^6 sum_p log p`.
   * This exceeds the elementary exponent `(1/2) sum_p a_p log p - log sqrt2` whenever
     `max a_p <= 10^7`.
   * Hence on all such forms Matveev's bound is **weaker** than item 1.
5. **Lang–Waldschmidt shape.** Any bound of the shape `|lambda| >= kappa_k prod_p M(alpha_p)^(-1-eps)`
   is `<= kappa_k rad(n)^(-eps) n^(-1/2)` when all `a_p <= 2`, since
   `prod p^(-1-eps) <= prod p^(-a_p/2 - eps)`.
   * So it does not improve the exponent `1/2` of item 1.
   * Such a bound is weaker than item 1 as soon as `kappa_k rad(n)^(-eps) <= sqrt2`.

*Proof.*

1. **Item 1.** A chord of angle `lambda` has length `2R sin(|lambda|/2) <= R|lambda|`. Combine
   this with Lemma 3.1 and `|g| = R/sqrt n`.
2. **Item 2.** Use `|lambda| <= arc/R`, together with item 1.
3. **Item 3.** Here `u - i ubar = (A - B)(1 - i)` with `A - B = 1`. So `|z - w| = sqrt2 |g|`,
   and `|lambda_zw| = 2 arcsin(|z - w|/(2R)) = 2 arcsin(1/sqrt(2n))`. Also
   `gcd(u, i ubar) = 1`, because `gcd(A, A-1) = 1` and `2A - 1` is odd. So the pair norm is `n`.
4. **Items 4 and 5.** These are the displayed comparisons. ∎

*Data (check E).*

* In the record 10-point cluster (`N = 1176852625 = 5^3 13^2 17 29 113`, `C = 8.2393`,
  `growth/referee_families.md`), the 45 pair norms range from `1105` to `N`.
* The eight closest pairs have pair norms `1105` (two pairs), `1625` (two) and `4901` (four). For
  six of them the ratio `|lambda|/(sqrt2 n^(-1/2))` is `1.00`; for the other two it is `1.41`.
* The largest ratios are `623` and `763`, attained at the coprime pairs, where `n = N`.
* Matveev's lower bound is below `10^(-10^13)` for every pair. This uses exact `A_p` and a lower
  bound for `B`, which is conservative.
* The family of item 3 attains item 1 up to a relative error `<= 1.8*10^(-2)` for
  `2 <= A < 2000`, the worst case being `A = 2` (check E').

**Lemma 3.3 (product clusters).**

*Statement.* Let `N = N_1 ... N_m` with `m >= 3` and pairwise coprime `N_i`. For each `i` let
`a_i != b_i` be in `Gamma_{N_i}`. Suppose the `2^m` points
`z_eps = prod_i (b_i if eps_i = 1 else a_i)` lie on an arc of length `<= C sqrt R` of `Gamma_N`,
and suppose `R > 4C^2`. Then `N <= (C^2/2)^6`.

*Proof.*

1. **Lifting the angles.** Put `phi = C R^(-1/2) < pi/2`. Choose real lifts `phi_eps` of
   `arg z_eps` in one interval of length `phi`.
   - Points differing only in block `i` have ratio `b_i/a_i`.
   - So `phi_(eps + e_i) - phi_eps` lies in `[-phi, phi]` and is congruent to
     `lambda_i = arg(b_i/a_i)` modulo `2 pi`.
   - Hence it equals the representative of `lambda_i` in `(-pi, pi]`, independently of `eps`.
   - Therefore `phi_eps = phi_0 + sum eps_i lambda_i`, and the spread `sum_i |lambda_i| <= phi`.
2. **A lower bound per block.** By Prop. 3.2(1) inside the circle `Gamma_{N_i}`,
   `|lambda_i| >= sqrt2 N_i^(-1/2)`. With step 1 this gives `N_i >= 2 N^(1/2)/C^2` for every `i`.
3. **Conclusion.** Multiplying over the `m` blocks gives `N >= (2/C^2)^m N^(m/2)`, that is,
   `N^((m-2)/2) <= (C^2/2)^m`. For `m >= 3` this gives `N <= (C^2/2)^(2m/(m-2)) <= (C^2/2)^6`. ∎

So "cube" or product configurations obey `U` by the elementary bound alone. The same argument
applies to arithmetic progressions `Pi^j Pibar^(t-j)` on a circle.

*Check F.* On 3000 random block triples the actual arc constant is at least
`11.8 * sqrt2 * N^(1/12)`.

**What this shows.**

* Every single linear form attached to two points of a circle is controlled exactly, up to a
  constant, by the elementary bound. No transcendence statement, proved or conjectured, can say
  more about it.
* Systems that decompose into independent blocks are controlled by the elementary bound alone.
* So `U` (and a growth improvement, see Section 7) cannot imply a "quantitative linear
  independence of Gaussian prime arguments" statement about single forms. Any such consequence is
  either already elementary or false.
* The possible content is about **entangled systems** whose pairwise forms are simultaneously within
  bounded factors of the elementary bound, for example the Hadamard and Paley profiles of
  `growth/walsh-extraction.md` §10.
* `growth/auxiliary.md` (b) made the corresponding observation for Baker certificates: Liouville is
  attained, and Baker gains only on non-saturated perfect powers. Prop. 3.2(3)–(5) and Lemma 3.3 are
  the formal versions for arbitrary single forms and for products.

---------------------------------------------------------------------------------------------

## 5. Theorem 4: the split twin, Erdős–Rosenfeld and Ruzsa, and the balanced locus

### 5.1 Two `Q`-forms, and no homomorphism between their tori

Over `Q(i)`, the substitution `u_j = x_j + i y_j`, `v_j = x_j - i y_j` identifies `X_M` with
`{u_1 v_1 = ... = u_M v_M}`. Complex conjugation then acts by `u_j <-> v_j`. The split form
`X_M^split = {x_1 y_1 = ... = x_M y_M}` is the same variety with the trivial action.

* Rational points of `X_M` are `M`-tuples of Gaussian integers of equal norm, that is, Gaussian
  divisors `z` of `N` with `z zbar = N`.
* Rational points of `X_M^split` are `M`-tuples of factorisations `n = d e`.

Let `T_M` and `T_M^split` be the open tori.

**Proposition 4.1.** Complex conjugation acts by `-1` on the character group of `T_M` (tensored
with `Q`). Hence:

* `T_M` is `Q`-anisotropic, and `Hom_Q(T_M, G_m) = 0`;
* `Hom_Q(T_M, T_M^split) = 0`;
* `Hom_Q(T_M^split, T_M) = 0`.

*Proof.*

1. **The character lattice.** Over `Qbar`, `T_M` is `{(u, v) : u_j v_j = t for all j}/G_m`, where
   `G_m` acts by `lambda(u, v) = (lambda u, lambda v)`, so `t -> lambda^2 t`.
   - With coordinates `(u_1, ..., u_M, t)`, the characters are `u^a t^c` with
     `sum a_j + 2c = 0`. This is a lattice `L` of rank `M`.
   - The characters `psi_j = u_j^2/t = u_j/v_j` lie in `L` and span `2L`, because
     `sum a_j (2 e_j - e_t) = 2(a, c)`.
   - Conjugation maps `psi_j` to `v_j/u_j = psi_j^(-1)`. So `sigma = -1` on `L` (which is
     torsion-free).
2. **No homomorphism.** A `Q`-homomorphism `T -> T'` is a Galois-equivariant map
   `X*(T') -> X*(T)`.
   - For `T = T_M`, `T' = T_M^split` (trivial action on `X*(T')`): the image consists of fixed
     vectors of `-1`, so it is `0`.
   - For `T = T_M^split`, `T' = T_M`: `phi*(-chi) = phi*(sigma chi) = phi*(chi)`, so
     `2 phi* = 0` and `phi* = 0`. ∎

*Remark.* `X_M` and `X_M^split` are both `Q`-rational, hence `Q`-birational. Prop. 4.1 excludes
only **equivariant** transfers. Those are the transfers that would map the multiplicative
structure of Gaussian divisors to that of rational divisors, which is how every reduction we know
between divisor problems works.

### 5.2 The split twin and its vertex

Write `E = n/D` and define

```text
S(n, D, C) = #{ d | n : D <= d <= D + C D^(3/2) n^(-1/2) }      (1 <= D <= sqrt n).
```

For `d` in this window the co-divisors satisfy `|e - e'| = n|d - d'|/(d d') <= C E^(1/2)` and
`|d - d'| <= C D E^(-1/2) <= C E^(1/2)`. So `S` counts `M`-tuples of `X_M^split(Q)` within
`C H^(1/2)` of `Y^split`, with `H = E` the height. This is the same approximation condition as
`A_M` on the split form.

**The split uniform statement `U^split`.** For every `C` there is `K(C)` with
`S(n, D, C) <= K(C)` for all `n` and all `D <= sqrt n`.

Its absolute vertex case, `D = sqrt n` with one `K` for every `C` once `n >= n_0(C)`, is the
Erdős–Rosenfeld question (Erdős problem #887). Ruzsa's question (#886) is the vertex statement
with window `n^(1/2 - eps)`. Vojta's conjecture with `D = 0` on the split model gives both
(`growth/conditional.md` Cor. 5.3, details in `referee_conditional.md` §1).

**Lemma 4.2 (vertex and rational directions: `C`-dependent bounds are trivial there).**

(a) *(split vertex)* `#{d | n : sqrt n <= d <= sqrt n + C n^(1/4)} <= 1 + C^2`.

(b) *(circle, rational directions)* Let `w` be a nonzero Gaussian integer and `0 < psi <= pi/4`.
Then at most `2(R |w| psi^2 + 1)` points `z` of `Gamma_N` satisfy `|arg z - arg w| <= psi`
(mod `2pi`). Consequently, an arc of length `C sqrt R` whose midpoint direction is within
`kappa R^(-1/2)` of `arg w` holds at most `2|w|(kappa + C)^2 + 2` points.

*Proof.*

* **(a)** `s(d) = d + n/d` lies in `[2 sqrt n, 2 sqrt n + C^2]`, because
  `s - 2 sqrt n = (d - sqrt n)^2/d`. Also `s` is an integer and determines `d > sqrt n`. (The
  erdosproblems page for #887 records this bound.)
* **(b)** Put `zeta = z wbar = X + iY`.
  - Then `|zeta| = R|w|` and `|arg zeta| <= psi`.
  - So `X > 0` and `R|w| - X = Y^2/(R|w| + X) <= R|w| psi^2`.
  - `X` is an integer in an interval of length `R|w| psi^2`, `Y` is determined up to sign by `X`,
    and `z = zeta/wbar`.
  - For the arc, all its points are within `(kappa + C/2) R^(-1/2)` of `arg w`. ∎

On both forms, then, the `C`-dependent problem is trivial near rational directions of small height:
at the vertex of the hyperbola, and near Gaussian directions of small norm on the circle.

### 5.3 Away from the vertex the split twin is a theorem (Coppersmith)

**Lemma 4.4 (Coppersmith–Howgrave-Graham counting).** Let `n >= 2` and `P` be integers, and let
`0 < beta <= 1`, `m >= 1`, `t >= 0`, `omega = m + 1 + t` and a real `X >= 1` satisfy

```text
omega^omega n^(m(m+1)/2) X^(omega(omega-1)/2) < n^(beta m omega).          (*)
```

Then at most `omega - 1` integers `x` with `|x| <= X` have a divisor `d >= n^beta` of `n` with
`d | P + x`.

*Proof.*

1. **The lattice.** Put `g_j(y) = n^(m-j)(y + P)^j` for `0 <= j <= m`, and
   `g_j(y) = y^(j-m)(y + P)^m` for `m < j < omega`.
   - The coefficient vectors of `g_j(Xy)` form a lower-triangular basis of a lattice `Lambda` in
     `R^omega`.
   - The diagonal entries are `n^max(m-j,0) X^j`. So `det Lambda = n^(m(m+1)/2) X^(omega(omega-1)/2)`.
2. **A short vector.** Minkowski's theorem for the cube gives a nonzero `v` in `Lambda` with
   `||v||_infinity <= (det Lambda)^(1/omega)`. It is the coefficient vector of `h(Xy)`, where `h`
   is a nonzero integer combination of the `g_j`, of degree `< omega`.
   - For `|x| <= X`: `|h(x)| <= sum_i |h_i| X^i <= omega ||v||_infinity`.
   - By (*) this is `< n^(beta m)`.
3. **Divisibility.** If `d | n`, `d | P + x` and `d >= n^beta`, then every `g_j(x)` is divisible by
   `d^m`, and so is `h(x)`. Since `|h(x)| < n^(beta m) <= d^m`, we get `h(x) = 0`.
4. **Counting.** A nonzero polynomial of degree `<= omega - 1` has at most `omega - 1` roots. ∎

**Corollary 4.5 (`U^split` off the vertex).** Let `0 < eta <= 1/4` and `C >= 1`. Suppose

```text
log n > (16/(5 eta)) ( log C + 2 log(6/eta + 5) ).
```

Then for every `D <= n^(1/2 - eta)`,

```text
S(n, D, C) <= max(C + 1, 6/eta + 4).
```

*Proof.*

1. **Small `D`.** If `D <= n^(1/3)`, the window has length `<= C`.
2. **The parameters.** Otherwise put `beta = log D/log n`, which lies in `(1/3, 1/2 - eta]`. Take
   `P = D`, `X = C D^(3/2) n^(-1/2)` (which is `>= 1`), `m = ceil(2/eta)` and
   `omega - 1 = w = ceil(m/beta) <= 3m + 1 <= 6/eta + 4`.
3. **Reducing (*).** Taking `log_n` and dividing by `omega(omega-1)/2`, condition (*) reads

   ```text
   2 log_n(omega)/w + m(m+1)/(w(w+1)) + 3 beta/2 - 1/2 + log_n C < 2 beta m/w.
   ```
4. **A lower bound for the right side.** From `w in [m/beta, m/beta + 1)`:

   ```text
   2 beta m/w - m(m+1)/(w(w+1)) > beta^2 - beta^2 (1 + beta)/(m + beta) >= beta^2 - 3/(8m).
   ```
5. **Using the margin.** Since `beta^2 - 3 beta/2 + 1/2 = (1 - beta)(1/2 - beta) >= eta/2` and
   `3/(8m) <= 3 eta/16`, it suffices that `log_n C + 2 log_n(omega)/w < 5 eta/16`. That holds
   under the hypothesis on `n`.
6. **Conclusion.** Every divisor in the window is `>= D = n^beta`, so Lemma 4.4 applies. ∎

**Remark 4.6 (the margin is the imbalance, and it vanishes identically on the circle).**

1. *(The split margin.)* At a divisor of size `n^beta` and the natural window
   `n^(3 beta/2 - 1/2)`, the Coppersmith margin is

   ```text
   beta^2 - (3 beta/2 - 1/2) = (1 - beta)(1/2 - beta).
   ```

   It is positive for `beta < 1/2` and zero exactly at the balanced point `beta = 1/2`.
   * Near the vertex, write `D = sqrt(n)/lambda` and apply Corollary 4.5 with
     `eta = log lambda/log n`.
   * Its hypothesis holds once `(C log n)^8 <= lambda <= n^(1/4)` and `n` is large. It then gives
     `S <= 6 log n/log lambda + 4`.
   * At `lambda = (log n)^O(1)` this is `O(log n/loglog n)`, the same shape as the circle's local
     bound.
2. *(The circle.)* Every point of `Gamma_N` is a Gaussian divisor of `N` with
   `|z| = |zbar| = |N|^(1/2)`. That is, `beta = 1/2` **at every point and in every direction**.
   This is the anisotropy of Prop. 4.1.
   * The Gaussian analogue of Lemma 4.4 (same construction over `Z[i]`; not proved here, and not
     used) would therefore have margin `1/4 - 1/4 - log C/log N < 0` at the natural window
     `C|N|^(1/4)`.
   * So it covers only arcs of length `sqrt(R)/K`, with `O(log R/log K)` points.
3. *(The refereed statement.)* `threshold/construct.md` (Thm C2, Cor. C3; refereed) proves the
   corresponding statement for all forced-divisibility Liouville arguments on the circle. They are
   exactly critical (ratio `1`), never supercritical.
4. *(Summary.)* The single-place auxiliary-polynomial barrier `tau*(M) < 2` (`threshold/upper.md`)
   and the Coppersmith barrier are the same phenomenon. On the split form, only the balanced locus
   is critical. On the anisotropic form, every point is balanced.

**Corollary 4.7 (Ruzsa implies `U^split` near the vertex).**

*Statement.* Assume #886: for each `eps > 0` there are `K(eps)` and `n_1(eps)` such that `n` has at
most `K(eps)` divisors in `(sqrt n, sqrt n + n^(1/2-eps)]` for `n >= n_1(eps)`.

Then for `n >= n_1(eps)` with `C n^(1/4) <= n^(1/2-eps)`, and for
`sqrt n/(1 + n^(-eps)) <= D <= sqrt n`, we have `S(n, D, C) <= 2K(eps) + 1`.

*Proof.*

1. **Divisors `d <= sqrt n` in the window.** Their co-divisors lie in
   `[sqrt n, n/D] subset [sqrt n, sqrt n + n^(1/2-eps)]`.
2. **Divisors `d > sqrt n` in the window.** They lie in `(sqrt n, sqrt n + C n^(1/4)]`.
3. Each range holds at most `K(eps)` divisors, apart from `d = sqrt n`, which accounts for the
   `+1`. ∎

**Remark 4.8 (the residual range).** Lemma 4.2(a), Corollary 4.5 and Corollary 4.7 leave open
`U^split` only in the range `1 + n^(-eps) < lambda = sqrt(n)/D < n^eta` (for all fixed
`eps, eta > 0`), at slopes `lambda^2` not close to rationals of small height.

* There it is an "off-centre Erdős–Rosenfeld" problem: divisors in windows of length
  `~ C n^(1/4)` around `sqrt(n)/lambda`.
* I did not find it in the literature.
* Its circle counterpart is `U` in directions badly approximable by Gaussian directions of small
  norm. By Lemma 4.2(b) that is essentially all of `U`.

So, under the dictionary of Prop. 4.1, the **whole** circle problem corresponds to the **residual
near-balanced range** of the split problem. Erdős #887 and #886 are the absolute statements at the
split vertex, where the `C`-dependent statement is trivial.

*Check H.* The proof of Lemma 4.4 was executed with exact LLL in place of Minkowski. The setup was:

* `n ~ 10^37` with three divisors in `[P, P + 24]`, `P ~ n^(1/3)`;
* `m = t = 3`, `omega = 7`.

Each reduced polynomial satisfies `sum |h_i| X^i < P^m`. In all 10 cases it has exactly the three
divisor offsets as its integer roots in `[0, 24]`.

*Check I.* For all `n <= 2*10^5`, the maximum of the vertex count is `1, 4, 4` for `C = 1, 2, 3`.
This is within `1 + C^2`.

---------------------------------------------------------------------------------------------

## 6. Residue classes cannot see clusters

**Proposition 6.1.** Let `z_1, ..., z_M` be distinct points of `Gamma_N` on an arc of length
`<= C sqrt R`, and let `sigma` be in `Z[i]` with `|sigma| > C N^(1/4)`.

* The `z_j` are pairwise incongruent modulo `sigma`.
* In particular this holds for every modulus in the Gaussian Lenstra range
  `Norm(sigma) >= Norm(N)^alpha = N^(2alpha)` with `alpha > 1/4`, once
  `N^(alpha - 1/4) > C`.

*Proof.* If `sigma | z_j - z_k` and `z_j != z_k`, then `|z_j - z_k| >= |sigma| > C sqrt R`. That
exceeds the arc length. ∎

So bounds of Lenstra type (at most `O(1)` divisors of `B` in **one** residue class modulo a modulus
above the critical exponent) give exactly one point per class. They say nothing about how many
classes a cluster occupies. `gaussian_strip_divisor_route_audit.md` §1 made the same observation for
its divisor region.

The natural lift `zeta_j = z_j zbar_1` puts a cluster into the class `0` modulo `zbar_1`. That
modulus has `Norm(zbar_1) = N = Norm(N^2)^(1/4)`, which is exactly critical. But the class is not
coprime to the modulus, which is outside Lenstra's hypothesis. Dividing out returns the original
cluster.

The archimedean analogue of Lenstra's problem at the critical exponent is `U` itself (Remark 4.6).

---------------------------------------------------------------------------------------------

## 7. The explanation, and what an implication would need

**The precise reason `U` is hard.**

1. **Shape.** `U` is a finiteness statement for approximations of the rational line `Y` in
   `X_M` at the Dirichlet exponent `1/2` (Theorem 1). Its absolute form is Roth-type finiteness
   for one fixed `M >= 9`.
2. **What happens below Vojta.** Every standard input below Vojta's conjecture provably loses at
   that exponent on this variety:
   * pair and local data are Plotkin-critical at `log R/loglog R` (`two_thirds/proof.md`);
   * additive energy is elementary at this scale (Theorem 2);
   * single linear forms are elementary and optimally so (Theorem 3);
   * residue classes see one point per class (Prop. 6.1);
   * forced divisibility and Coppersmith are exactly critical, because every point is balanced
     (Remark 4.6, `construct.md`);
   * single-place auxiliary polynomials are capped by `tau* < 2` (`upper.md`);
   * Ru–Vojta and Subspace give `beta_M < 2`, or exceptional sets depending on `S`;
   * `abc` is vacuous on squarefree blocks (`conditional.md` §8);
   * multi-row characters fail on fully flipped Paley profiles (`walsh-extraction.md`).
3. **Vojta suffices.** Vojta (`D = 0`) on `Bl_Y X_7` or on `Mbar_{0,8}` implies `U_abs`.
4. **Why the circle is harder than its twin.** The anisotropy of the torus (Prop. 4.1) is the
   structural reason the circle problem is harder than its split twin.
   * On the split form, the imbalance between the two torus coordinates is a resource worth
     `(1-beta)(1/2-beta)` to Coppersmith's method.
   * The vertex trick handles the balanced rational point.
   * On the circle, complex conjugation sets the imbalance to zero everywhere.
   * The vertex trick works only near Gaussian directions of small norm, which by Lemma 4.2(b)
     cover a negligible set of arcs.

**Why the candidate implications cannot come from the standard mechanism.** An implication
"`U` ⇒ P" needs a construction that turns counterexamples of `P` into short-arc clusters. By the
results above, such a construction must have all of the following properties.

* *It cannot be a homomorphism of tori*, by Prop. 4.1. This excludes transfers from rational divisor
  problems (Erdős–Rosenfeld, Ruzsa, Lenstra) by multiplicative dictionaries.
* *It must produce clusters whose pair forms are within bounded factors of the elementary bound, in
  an entangled (non-product) pattern* (Prop. 3.2, Lemma 3.3). So the counterexamples of `P` must
  carry `Theta(M^2)` simultaneously near-extremal approximations. A single small linear form,
  which is what Baker-type or Lang–Waldschmidt-type problems are about, produces at most a
  product or progression pattern, and those satisfy `U` unconditionally.
* *It must use `L^infinity` information about two separated square-root intervals*
  (Prop. 2.3). Statements about squares in one progression only see one interval.

None of the recognised problems examined has a natural source of such configurations. The
Erdős–Rosenfeld and Ruzsa problems are the closest, but they are the split twin's **vertex**
statements and absolute in `C`.

**Growth improvement.** The same obstructions apply to `M = o(log R/loglog R)`. This is a remark,
not a theorem:

* Theorems 2 and 3 and Prop. 6.1 bound the information those inputs carry, independently of `M`.
* The known local models reaching `log R/loglog R` (cube and Hadamard profiles) are made of
  elementary-extremal pairs.
* Corollary 4.5 shows the split twin has a growth-type bound `O(log n/log lambda)` from
  Coppersmith alone, away from the vertex. The circle has no analogue of this, because
  `lambda = 1` identically.

I did not find a statement that a growth improvement would imply and that is open in the
literature.

**Open points left by this note.**

1. Is the finiteness statement in Theorem 1(b) a special case of a named conjecture weaker than
   Vojta's? For example: a McKinnon-type "best approximations lie on rational curves" statement,
   combined with Bezout on the finitely many curve families. A positive answer would give the
   requested equivalence.
2. The residual near-balanced range of the split twin (Remark 4.8). Is "off-centre
   Erdős–Rosenfeld" open? A proof there by `Q`-form-independent methods would likely transfer to
   the circle.
3. The Lenstra problem at exactly `alpha = 1/4` is the non-archimedean, split, balanced twin. The
   bound of Coppersmith–Howgrave-Graham–Nagaraj, `O((alpha-1/4)^(-3/2))`, diverges there. I could
   not access the sources to confirm its status at `alpha = 1/4`.

---------------------------------------------------------------------------------------------

## 8. Exact checks (`round4/hardness_checks/`)

All integer arithmetic is exact. Floating point is used only for printed angles and logarithms.

| script | claim checked | result (log file) |
|---|---|---|
| `check_energy.py` (A) | Lemma 2.1 on 60 random intervals, `X` up to `4*10^6`, `c` in `{0.5, ..., 3}` | max multiplicity `<= 1 + floor(c^2)` in all cases |
| (B) | Theorem 2(1) on 40 random pairs of subsets | all within bound |
| (C) | Prop. 2.3(2) on 30 one-interval boxes | max `r(n) <= 4`, within bound |
| (D) | Corollary 2.2 on random APs with `a >= (qK/2)^(4/3)` | worst `count/sqrt(6K) = 0.289` |
| `check_forms.py` (E) | Prop. 3.2 on the 10-point cluster (45 pairs) | elementary bound attained (ratio `1.00`) at the closest pairs; Matveev weaker on all pairs |
| (F) | Lemma 3.3 on 3000 random block triples | `C >= 11.8 sqrt2 N^(1/12)` |
| (G) | Lemma 4.2(b), 6072 cases | no violation |
| `check_split.py` (H) | Lemma 4.4 mechanism with exact LLL, `n ~ 10^37` | 30 divisor offsets are exact roots; HG condition holds |
| (I) | Lemma 4.2(a), all `n <= 2*10^5` | max `1, 4, 4 <= 1 + C^2` |
| (J) | data: `max_D S(n, D, 1)` for `n <= 4*10^5` | max `3`, both near and away from the vertex (small range; not used) |

Outputs: `check_energy_output.txt`, `check_forms_output.txt`, `check_split_output.txt`.

---------------------------------------------------------------------------------------------

## 9. Literature and novelty

**Grep of the research notes.** I searched for: Rudin, Rosenfeld, Ruzsa, Coppersmith, Lenstra,
Howgrave, Lang–Waldschmidt, Matveev, anisotropic, twin, sumset, energy.

* The notes contain:
  - Chan's almost-square results (`chan_almost_square_axis_audit.md`);
  - the Erdős–Rosenfeld/Chan comparison (`balanced_clique_residual_factorization.md` §7);
  - Hales' Gaussian Lenstra algorithm and the "one class per divisor" observation
    (`gaussian_strip_divisor_route_audit.md`);
  - Gabdullin's energy theorem and the difference-multiplicity bound
    (`endpoint_recent_literature_check.md`);
  - the Baker-certificate discussion (`growth/auxiliary.md`).
* Not in the notes or in `docs/claude_review`:
  - Theorem 1(b) and the "empty ⇔ finite ⇔ non-dense" form of Theorem 1(a);
  - Prop. 2.3 and the conclusion that Rudin-type links are unconditional at exponent `1/2`;
  - Prop. 3.2(3)–(5) and Lemma 3.3;
  - Prop. 4.1;
  - the split twin `U^split`, with Lemma 4.2(a) and (b) used as twins, Cor. 4.5, Remark 4.6 and
    Cor. 4.7.
* Lemma 4.4 is the standard Coppersmith–Howgrave-Graham argument. It is included so that
  Corollary 4.5 is self-contained.

**Web sources.** Direct fetches of arXiv, UCL, UAM, AMS and Granville's site were refused by the
egress proxy. Literature statements therefore come from search snippets and from the notes, and
should be re-checked against the primary sources.

* Erdős problems #886 (Ruzsa: `O_eps(1)` divisors in `(sqrt n, sqrt n + n^(1/2 - eps))`, open) and
  #887 (Erdős–Rosenfeld: an absolute `K` for `(sqrt n, sqrt n + C n^(1/4))`; the page notes the
  bound `1 + C^2`). Source: erdosproblems.com, via search.
* Chan (arXiv 1406.2230 and a later paper) proves Erdős–Rosenfeld for perfect squares.
* Coppersmith–Howgrave-Graham–Nagaraj, Math. Comp. 77 (2008): at most `O((alpha-1/4)^(-3/2))`
  divisors in a residue class for `alpha > 1/4`. Lenstra, Math. Comp. 42 (1984).
* Cilleruelo–Granville (CRM Proc. 43 (2007)) motivate their conjectures by Rudin's conjecture
  through `Lambda(p)` and sumsets of squares. I could not read the paper, so I make no claim about
  which of their implications coincide with Corollary 2.2.
* Oganesyan's claimed disproof of the `R^alpha` conjecture (arXiv 2107.09991) was withdrawn
  (recorded in the notes).

I found no published implication from the `C sqrt R` arc bound to a recognised open problem, and no
equivalence with a named conjecture. This is a search result, not a proof of absence.
