# Round 7, route G: genuine residues pin the angles. P3-genuine collapses to actual clusters

Status: draft, unrefereed. Companion to `round7/p3lower.md` and `round7/p3data.md`; same conventions.

Conventions are those of `round4/outside.md`, `round5/lower.md` and `round6/sharp.md`:
* `W = log N`, `R = N^(1/2)`, and the arc has angular width `Delta = C e^(-W/4)`;
* `w(v) = sum_j |v_j| log p_j`, `v(c) = sum_x c_x a_x`, `n = ||c||_1`, `s_c = w(v(c)) - W/2`;
* `mu_xy = d_xy - W/2` is the pair margin, `c_C = 2 log C - log 2`;
* `A_q = max{a : q^a <= X}` for odd primes `q`, and `L'_X = prod_(q <= X odd, q not dividing N) q^(A_q)`;
* `H = log^2 M / loglog M`, and `H_X = log X log M / loglog M` (so `H_X = (1 + o(1)) H` when `X = M^(1+o(1))`).

*Theorem*, *Lemma*, *Proposition* and *Corollary* mean proved here with all quantifiers, or cited
exactly and re-derived. *Computer-verified* means an exact finite computation whose every
certificate is re-checked by independent code. *Evidence* means a computation that no proof uses.
*Heuristic* means an argument that uses an unproved independence or equidistribution assumption,
which is named. Scripts are in `round7/p3G_checks/`; each runs from that directory with
`python3 <script>`, and its output is stored next to it. The repository and the round-5/6/8 folders
were not modified.

---------------------------------------------------------------------------------------------

## 0. Verdict

**Answer to the target.** Route G does not produce fake circles at `W = O(M log M)`.

Genuine residues, with every monomial strengthened, **pin the fake angles to the genuine ones**.
* *Saturation.* Genuine data saturate the per-unit height condition, `m_(v,s) = |D_s(A_v)|`, on
  a positive proportion of the low-weight monomials whenever the primes have polynomial size.
* *One-sided constraints.* Each saturated ("tight") monomial turns its slab condition into a
  one-sided constraint through the genuine angle.
* *Pinning.* Short circuits of such constraints force `phi_j - theta_j in (pi/2) Z`.
* *Collapse.* When every column is pinned, **every** member of `P3-G` on that prime set is an actual
  cluster, whatever the profile.

So on prime sets of polynomial size the class contains no fakes at all. This is computer-verified
for the sets of Theorem 4.1 and holds with probability `1 - o(1)` in a random-tightness model
(Theorem 5.3). A fake needs about `M/2` columns of size `log p >~ c log^2 M/loglog M` (proved
with `c = 1/2 - eps` in the model; `c ~ 1` heuristically). That is the same rate `Theta(M log^2 M/loglog M)` as
`sharp.md` Theorem 6.4 (P3-all with fake residues), not `O(M log M)`. Route G therefore does not
answer P3-all\*. It is not a lower-bound shortcut either: where the closure is complete, `P3-G` *is*
the actual problem.

| # | statement | class | status |
|---|---|---|---|
| 1.1-1.3 | definition of `P3-G_X`; `m_(v,s) = gcd(D_s(A_v), L'_X)`; actual clusters are members for every `X` | `P3-G ⊂` P3-all `⊂` P3-char `⊂ A_X` | proved |
| 2.1 | the genuine angles `theta` satisfy every monomial condition (planted point) | `P3-G` | proved |
| 2.3 | a tight `(v, s)` is the one-sided constraint `<v, phi - theta> in sigma [0, pi - 2d] + pi Z`, with exact integer formulas for `sigma` and `sin^2 d` | `P3-G` | proved |
| 2.4-2.6 | circuit pinning; the certified closure `Pi*(pi, X)` depends only on the prime set and binds every member at every `X' >= X` | `P3-G` | proved |
| **3.1** | **collapse:** if `2 e_j in Pi*` for all `j`, every member with these primes (any `M`, `C`, profile, `X' >= X`) is an actual cluster after `pi_j -> i^(t_j) pi_j` | `P3-G` | **proved** |
| 3.2 | dichotomy: every character vanishing on the unpinned columns `F` is a genuine near-relation, `abs(D_s(A'_(v(c)))) <= (sqrt2/4) n C e^(s_c/2)` | `P3-G` | proved |
| 3.3-3.4 | `rank V >= (1/2 - o(1)) M`; a member has `>= (1/2 - o(1)) M - k` unpinned columns or `k + 1` independent near-relations | P3-char / `P3-G`, `e = 1`, `X >= 2M` | proved |
| **4.1** | complete closures for random primes of size `M^3`, `M^4` (`M` = 44 to 644, `r` = 86 to 1304, `X = 2M`); hence `P3-G` on these prime sets = actual clusters | those prime sets | **computer-verified** (exact certificates, independent checker) |
| 4.3-4.4 | tight densities of actual primes = divisor density of `L'_X` (within 4%); pinning transition at pair degree 7-8, `log p* = (1.23-1.34) H` for `M = 84..1284`, growing faster than `log M` | actual primes | evidence |
| **5.3-5.4** | in the random-tightness model: `(f-1) rho_1^2 rho_2 >= 48 log f` pins all `f` columns w.p. `1 - O(log f/f^2)`; hence `>= (1/2 - o(1)) M` free columns of size `log p > (1/2 - eps) H_X`, and `W >= (1/4 - eps) M H_X` for members with `O(M/log M)` near-relations | model `T` | **proved in the model** |
| 5.3 (heur.) | critical size `lambda* ~ H_X` (full closure), adaptive selection changes it by `1 + o(1)` | model | heuristic |
| 6.1 | `P3-G_X` has non-actual members with `W <= (68.7 + o(1)) M X` | `P3-G` | proved |
| 6.2 | non-actual members at `W ~ (66 + eps) M H_X` | `P3-G` | heuristic (smooth-value statistics) |
| 7 | task (ii) is pre-empted in `P3-G`; for the characters-only genuine class `P3-G-char` it is `p3data.md`'s (Q) + H, open | -- | discussion |

**Task (i), completion for monomials outside `L`.** For genuine primes of size `M^(O(1))` the sum
`sum_v m_(v,s) e^(-w(v)/2)` is not small: on two verified sets its tight weight-2 part alone is
4052 and 15190. Worse, by Theorems 3.1/4.1 no completion exists for non-actual members on such
sets. The needed prime size, for about `M/2` columns:
* `log p >= (1/2 - eps) H_X` (model, proved);
* `log p ~ H_X` (heuristic; measured `(1.06-1.22) H_X` for the pair-only closure at `M <= 1284`).

Sufficient: `(2 + eps) H_X` (heuristic) or `2.08 X + O(log M)` (proved).

**Task (ii), pair demands.** Not the binding constraint in `P3-G`. At `W = O(M log M)` the collapse
comes first. At the scale where fakes can exist, uniform Paley margins already exceed the genuine
pair demands. Small-prime repair columns are themselves pinned and act as seeds against the free
columns (§7).

**Validity for actual circles.** Every inequality used is valid for every member of `P3-G`, hence
for actual clusters (Prop. 1.3). For actual clusters the conclusions of Theorems 3.1 and 3.2 hold
trivially, with `phi = theta`. They are obstructions to *fakes*, not to clusters. No growth
improvement for actual circles is claimed.

---------------------------------------------------------------------------------------------

## 1. The class P3-G

### 1.1 Definition

**Definition 1.1 (P3-G, genuine residues; `p3lower.md` Def. 1.1).** Fix `M`, `C > 0` and `X >= 3`.
A member of `P3-G_X(M, C)` is a member `(p, e, (a_x), phi, theta_0, delta, k, r)` of
`F^all_X(M, C)` (`sharp.md` Def. 1.1: conditions (a)-(d), (e-char) and (e-all)) together with
Gaussian primes `pi_1, ..., pi_r`, `pi_j conj(pi_j) = p_j`, such that

```text
r_j(q) = pi_j  mod q^(A_q)      for every odd prime q <= X with q != p_j.
```

Only the column angles `phi in R^r` (hence `delta`, `theta_0`, `k`) are free; they are **not**
required to be the arguments of the `pi_j`. Fix real lifts `theta_j` of `arg pi_j`; these are the
*genuine angles*. For `v in Z^r` put

```text
A_v = prod_j pi_j^(v_j^+) conj(pi_j)^(v_j^-) = X_v + i Y_v,     |A_v| = e^(w(v)/2),     arg A_v = <v, theta> mod 2 pi,
D_0 = Y_v,  D_1 = X_v - Y_v,  D_2 = X_v,  D_3 = X_v + Y_v,      nu_s = 1 (s even), 2^(-1/2) (s odd).
```

Every member of `P3-G` is a member of P3-all, P3-char and `A_X` (`sharp.md` §1.1). If `X' >= X`,
then `P3-G_(X')(M, C) ⊂ P3-G_X(M, C)`: every collision modulus at `X` divides the one at `X'`, so
(1.1) and (1.2) of `sharp.md` at `X'` imply them at `X`.

**Lemma 1.2 (the genuine collision moduli).** In a member of `P3-G_X`, for every `v != 0` and
`s in Z/4`,

```text
m_(v,s) = gcd( D_s(A_v), L'_X ),      D_s(A_v) != 0.
```

*Proof.*
* By definition `m_(v,s) = prod q^(a)` over odd `q <= X` not dividing `N`, where `a <= A_q` is
  maximal with `prod_j (r_j/conj r_j)^(v_j) = i^s mod q^a`.
* With genuine residues, `prod_j (pi_j/conj pi_j)^(v_j) = A_v/conj(A_v)`, and `A_v` is a unit
  mod `q` because its norm is prime to `q`. So the condition is `A_v = i^s conj(A_v) mod q^a`.
* `A - i^s conj(A)` equals `2iY`, `(1 - i)(X - Y)`, `2X`, `(1 + i)(X + Y)` for `s = 0, 1, 2, 3`.
  Since `q` is odd, the condition is `q^a | D_s(A_v)`.
* So the exponent of `q` in `m_(v,s)` is `min(v_q(D_s), A_q)`, which is the exponent of `q` in
  `gcd(D_s, L'_X)`.
* `D_s(A_v) != 0`: `A_v` is conjugate-primitive, of odd norm (`outside.md` Lemma 0.1, proof). QED

**Proposition 1.3 (actual clusters are members).** Let `z_1, ..., z_M` be a normalised actual cluster
on an arc of length `<= C sqrt R`. For every `X >= 3` it defines a member of `P3-G_X(M, C)`, with
`phi = theta`.

*Proof.* `sharp.md` Prop. 1.2 builds a member of `F^all_X(M, C)` with `phi_j = arg pi_j` and
`r_j(q) = pi_j mod q^(A_q)`; these residues are genuine (`p3lower.md` Prop. 1.2, which re-reads
that proof). QED

**Remark 1.4 (heights are automatic, and saturated).** Lemma 1.2 gives
`m_(v,s) <= |D_s(A_v)| = nu_s^(-1) |A_v| |sin(arg A_v - s pi/4)| <= nu_s^(-1) e^(w(v)/2)`. So the per-unit
height condition (`sharp.md` Prop. 6.1, as corrected by the round-6 referees) holds for genuine
data. Section 2 uses the case of **equality** `m_(v,s) = |D_s(A_v)|`. It happens for a positive
proportion of the low-weight monomials when the primes have polynomial size (§4), and it is what
pins the angles.

**Remark 1.5 (parity).** `X_v + Y_v` is odd, since `X_v^2 + Y_v^2` is odd. So `D_1`, `D_3` are odd,
and exactly one of `D_0`, `D_2` is odd. If every `pi_j = a_j + i b_j` is normalised with `b_j` even,
then every `A_v` has `Y_v` even (that set is closed under products and conjugation). Then `D_0` is
even for every `v`. A unit change `pi_j -> i pi_j` permutes the four `s` and changes nothing below.

---------------------------------------------------------------------------------------------

## 2. Tight relations and the pinning closure

**Lemma 2.1 (the planted point).** In a member of `P3-G_X`, the genuine angles `theta` satisfy
(1.2) of `sharp.md` for every `v != 0` and every `s`, and both gap inequalities of `outside.md`
Lemma 0.1.

*Proof.* `|sin(<v, theta> - s pi/4)| = nu_s |D_s(A_v)| e^(-w(v)/2) >= nu_s m_(v,s) e^(-w(v)/2)`, since
`m_(v,s)` divides `D_s(A_v) != 0` (Lemma 1.2). Lemma 0.1 is the case `m = 1`. QED

The fake angles `phi` must lie in the realisation coset of the fake arc (`outside.md` Thm A.1(i)),
and `theta` in general does not. The question is how far `phi` may move away from `theta`.

**Definition 2.2 (tight relations).** A pair `(v, s)`, `v in Z^r \ {0}`, is *tight* (for the primes
`pi` and the level `X`) if `m_(v,s) = |D_s(A_v)|`. By Lemma 1.2 this means that `|D_s(A_v)|` divides
`L'_X`: it is odd, prime to `N`, and every prime power exactly dividing it is a `q^a <= X`. If
`(v, s)` is tight at `X`, it is tight at every `X' >= X`.

**Lemma 2.3 (a tight relation is a one-sided arc).** Let `(v, s)` be tight. Put
`alpha = <v, theta> - s pi/4`, `d = dist(alpha, pi Z)` and `sigma = sign sin(2 alpha)`. Then:
1. `0 < d < pi/2`, and `sin^2 d = nu_s^2 D_s(A_v)^2 / Norm(A_v)`;
2. `sigma = sign(X_v Y_v)`, `sign(Y_v^2 - X_v^2)`, `-sign(X_v Y_v)`, `sign(X_v^2 - Y_v^2)` for
   `s = 0, 1, 2, 3`;
3. for every `phi in R^r`, (1.2) holds for `(v, s)` if and only if

   ```text
   x := <v, phi - theta>  in  sigma [0, pi - 2d] + pi Z.                                       (2.1)
   ```

*Proof.*
1. `sin d = |sin alpha| = nu_s |D_s(A_v)|/|A_v|` by the identity in Lemma 2.1. This is `> 0` since
   `D_s != 0`. It is `< 1`: for `s` even, `|D_s| = |Y|` or `|X|`, both `< |A|` as `X, Y != 0`;
   for `s` odd, `|X -+ Y| = sqrt2 |A|` would force `X = +-Y`, impossible as `X + Y` is odd.
2. `sin 2 alpha = sin(2 arg A_v) cos(s pi/2) - cos(2 arg A_v) sin(s pi/2)`, with
   `sin(2 arg A) = 2XY/|A|^2` and `cos(2 arg A) = (X^2 - Y^2)/|A|^2`.
3. `|sin t|` is an increasing function of `dist(t, pi Z) in [0, pi/2]`. Since `m = |D_s|`, the right
   side of (1.2) is `|sin alpha| = sin d`, so (1.2) says `dist(alpha + x, pi Z) >= d`.
   * If `sigma = +1`, then `alpha = d mod pi`, and `dist(d + x, pi Z) >= d` iff
     `d + x mod pi in [d, pi - d]` iff `x mod pi in [0, pi - 2d]`.
   * If `sigma = -1`, then `alpha = -d mod pi`, and the same computation gives
     `x mod pi in [-(pi - 2d), 0]`. QED

So `phi = theta` lies on the boundary of every tight constraint, and `x` may only move to one side.

**Lemma 2.4 (circuit pinning).** Fix `phi` and a subgroup `Pi ⊂ Z^r` with `<u, phi - theta> in pi Z`
for all `u in Pi`. Let `(v_1, s_1), ..., (v_k, s_k)` be tight, with data `(sigma_i, d_i)` as in
Lemma 2.3, and let `n_1, ..., n_k >= 1` be integers with

```text
sum_i n_i sigma_i v_i  in  Pi      and      sum_i n_i (pi - 2 d_i)  <  pi.                       (2.2)
```

If `phi` satisfies (1.2) for these `k` relations, then `<v_i, phi - theta> in pi Z` for every `i`.

*Proof.* Put `x_i = <v_i, phi - theta>`. By (2.1), `sigma_i x_i = y_i + pi k_i` with
`y_i in [0, pi - 2 d_i]` and `k_i in Z`. By the first part of (2.2),
`sum_i n_i sigma_i x_i = <sum n_i sigma_i v_i, phi - theta> in pi Z`. Hence `sum_i n_i y_i in pi Z`.
But `0 <= sum n_i y_i <= sum n_i (pi - 2 d_i) < pi`. So `sum n_i y_i = 0`, every `y_i = 0`, and
`x_i in pi Z`. QED

All the data in (2.2) are exact: `sigma_i` by Lemma 2.3(2), and the arc condition by
`d_i in (0, pi/2)` and Lemma 2.3(1). For `k = 2`, `n_1 = n_2 = 1` it reads
`sin^2 d_1 + sin^2 d_2 > 1`; for `k = 1`, `n_1 = 2` it reads `sin^2 d_1 > 1/2`.

**Definition 2.5 (certified closure).** For Gaussian primes `pi = (pi_1, ..., pi_r)` and a level `X`,
let `Pi_0 = {0}`. Let `Pi_(t+1)` be generated by `Pi_t` and every `v_i` that occurs in some family of
relations tight at `X` satisfying (2.2) relative to `Pi_t`. Put `Pi*(pi, X) = union_t Pi_t`, and call column `j` *pinned*
if `2 e_j in Pi*` and *exactly pinned* if `e_j in Pi*`. Write `F(pi, X)` for the set of unpinned
columns.

**Corollary 2.6.** Let `pi` and `X` be given. In every member of `P3-G_(X')(M, C)` with these primes,
for any `X' >= X`, any `M`, `C`, any exponents and rows, and any angles:
`<u, phi - theta> in pi Z` for all `u in Pi*(pi, X)`. In particular `phi_j - theta_j in (pi/2) Z` for
every pinned column, and `phi_j - theta_j in pi Z` for every exactly pinned column.

*Proof.* The member satisfies (1.2) for every `(v, s)`, at level `X'`. A relation tight at `X` is
tight at `X'` with the same `(sigma, d)`. Induction on `t` with Lemma 2.4 gives the claim for each
`Pi_t`. QED

The definition does not depend on the profile: it is a property of the **prime set** and of `X`.
The computations below use a sub-closure: relations of weight `||v||_1 <= 2` (singles `e_j`,
doubles `2 e_j`, pairs `e_i +- e_j`; in one experiment also triples), and families with `k <= 2`.
That only makes the computed group smaller, so Corollary 2.6 applies to it as well.

---------------------------------------------------------------------------------------------

## 3. Collapse and the dichotomy

**Theorem 3.1 (collapse).** Let `pi_1, ..., pi_r` be Gaussian primes of distinct norms
`p_j = 1 mod 4`, and `X >= 3`, such that every column is pinned: `2 e_j in Pi*(pi, X)` for all `j`.
Then for every `X' >= X`, `M >= 2`, `C > 0`, every member of `P3-G_(X')(M, C)` with these primes is an
actual cluster. Precisely: there is `t in Z^r` such that, with `pi'_j = i^(t_j) pi_j`, the Gaussian
integers

```text
zeta_x = i^(-k_x) prod_j pi'_j^(a_xj) conj(pi'_j)^(e_j - a_xj)          (x = 1, ..., M)
```

are `M` distinct lattice points of `x^2 + y^2 = N`, with `arg zeta_x = theta_0 + delta_x`
(mod `2 pi`). They lie on an arc of angular width `Delta`, i.e. of length `<= C sqrt R`.

*Proof.*
* By Corollary 2.6, `phi_j - theta_j = (pi/2) t_j` with `t_j in Z`, for every `j`. Then
  `theta'_j := phi_j` is a real lift of `arg pi'_j`.
* `prod_j pi'_j^(a_xj) conj(pi'_j)^(e_j - a_xj)` has modulus `N^(1/2)` and argument
  `<a_x, phi> - <e - a_x, phi> = <2 a_x - e, phi>` (mod `2 pi`). By (1.1) of `outside.md` this is
  `theta_0 + delta_x + (pi/2) k_x`. Multiplying by `i^(-k_x)` gives `arg zeta_x = theta_0 + delta_x`.
* Distinct rows give non-associate Gaussian integers, since `pi'_j` and `conj(pi'_j)` are
  non-associate primes for `p_j = 1 mod 4`. So the `zeta_x` are distinct.
* `delta_x in [0, Delta]` and `R Delta = C R^(1/2)`. QED

So on such a prime set **the class `P3-G` contains no fake circle at all**: its members are actual
clusters, with the member's own `delta` as their actual angles.

**Theorem 3.2 (dichotomy: free columns or genuine near-relations).** Take a member of `P3-G_(X')`
with primes `pi` and `X' >= X`. Let `F = F(pi, X)` be the unpinned columns, and choose `t_j` with
`phi_j - theta_j = (pi/2) t_j` for `j not in F` (Corollary 2.6); put `pi'_j = i^(t_j) pi_j` for
`j not in F`. Let `c in Z^M`, `sum c = 0`, with `v(c)_j = 0` for every `j in F`, and
`A' = prod_(j not in F) pi'_j^(v_j^+) conj(pi'_j)^(v_j^-)` with `v = v(c)`. Then

```text
min_s  nu_s |D_s(A')|  <=  |A'| sin(n Delta / 4),     hence     |D_(s')(A')|  <=  (sqrt2/4) n C e^(s_c/2)   for some s'.
```

If moreover `v(c) != 0`, we call `c` a *genuine near-relation* of the member.

*Proof.*
* Summing (1.1) of `outside.md` with weights `c_x` gives the real identity
  `<2 v(c), phi> = c . delta + (pi/2) c . k`, because `sum c = 0`.
* `v(c)` vanishes on `F`, and `phi_j = theta'_j := theta_j + (pi/2) t_j` off `F`. So
  `<v(c), theta'> = c . delta / 2 + (pi/4) c . k`. Since `theta'_j` lifts `arg pi'_j`,
  `arg A' = <v(c), theta'>` (mod `2 pi`).
* `|c . delta| <= (n/2) Delta`, because `delta in [0, Delta]^M` and the positive part of `c` has
  mass `n/2`. So `dist(arg A', (pi/4) Z) <= n Delta/4`.
* `(pi/4) Z` is the union of the classes `s pi/4 + pi Z`, `s = 0..3`, and
  `dist(arg A', s pi/4 + pi Z) = arcsin(nu_s |D_s(A')|/|A'|)` (Lemma 2.3(1) without tightness).
* Finally `|A'| = e^(w(v)/2)` and `Delta = C e^(-W/4)`, and `nu_s >= 2^(-1/2)`. QED

*Interpretation.* The actual points
`i^(-k_x) prod_(j not in F) pi'_j^(a_xj) conj(pi'_j)^(e_j - a_xj)` (the cluster with the free columns
removed) satisfy the character relation `c` with the same precision `n Delta/4` as the fake. For a
pair `c = e_x - e_y` (rows `x`, `y` equal on `F`) it says: these two actual points are within angle
`Delta` of each other up to a unit. If `F` is empty, every `c` is a near-relation, and this is
Theorem 3.1 again.

**Proposition 3.3 (the character lattice has rank at least about `M/2`).** Let a member of
`F^char_X(M, C)` with `X >= 2M` and `e = 1` be given, and `V = span_R{a_x - a_y}`. Then

```text
rank V  >=  M / (2 + 4 M c_C^+ / W)  -  2  =  (1/2 - o(1)) M.
```

*Proof.*
* Let `S = 2A - J` be the `M x r` sign matrix of the profile (`A` the exponent matrix), and
  `G = S diag(w) S^T`. Then `G` is positive semidefinite, `G_xx = W`, and `G_xy = W - 2 d_xy = -2 mu_xy`.
* Eigenvalue Cauchy-Schwarz: `rank G >= (tr G)^2 / tr(G^2) = M^2 W^2 / (M W^2 + 4 sum_(x != y) mu_xy^2)`.
* `|mu_xy| <= W/2` (as `|G_xy| <= W`). The pair inequality (`p3lower.md` Lemma 1.3) gives
  `mu_xy >= D_xy - c_C >= -c_C`. With the cut identity `sum_(x<y) mu_xy <= M W/4`
  (`p3lower.md` Lemma 1.4(4)):
  `sum_(x != y) |mu_xy| <= M W/2 + 2 M^2 c_C^+`.
* Hence `4 sum mu^2 <= 4 (W/2)(M W/2 + 2 M^2 c_C^+) = M W^2 + 4 M^2 W c_C^+`, and
  `rank G >= M/(2 + 4 M c_C^+/W)`.
* `rank G <= rank S <= rank A + 1`, and `V` is the row space of `(I - J/M) A`, of rank
  `>= rank A - 1`.
* `W >= 3 M log M - (14 + 4 log^+ C) M` (`sharp.md` Prop. 1.3), so `4 M c_C^+ / W -> 0`. QED

For the Paley profiles of rounds 5/6, `rank V = M - 1` (`sharp.md` Lemma 4.1).

**Corollary 3.4 (what a non-actual member needs).** In a squarefree member of `P3-G_(X')`,
`X' >= max(X, 2M)`, the genuine near-relations together with `0` form a group whose rank modulo
`{c : v(c) = 0}` is at least `rank V - |F(pi, X)|`. So for every `k >= 0`, either

* at least `(1/2 - o(1)) M - k` columns are unpinned, or
* there are `k + 1` independent genuine near-relations.

*Proof.* The `c` with `sum c = 0` and `v(c)|_F = 0` form a group `C_F`. The map `c -> v(c)` sends
`C_F` onto a lattice of full rank in the rational subspace `V cap {v : v|_F = 0}`, whose dimension
is `dim V - dim(V|_F) >= rank V - |F|`. Every `c in C_F` with `v(c) != 0` is a genuine
near-relation (Theorem 3.2). Apply Proposition 3.3. QED

**Remark 3.5 (where the near-relations could come from).** A near-relation is an arithmetic
coincidence among actual lattice points: the argument of a genuine monomial of weight
`w(v(c)) >= W/2 - O(log n)` lies within `n Delta/4 = (nC/4) e^(-W/4)` of `(pi/4) Z`. For a fixed
character and "generic" primes this has probability about `(2 n C/pi) e^(-W/4)` (*heuristic*:
equidistribution of the arguments of Gaussian integers; for products of distinct Gaussian primes
the arguments are equidistributed by Hecke, but no uniformity in the needed range is claimed).
A builder who chooses the rows among the `2^r` lattice points of the circle has about `M r log 2`
nats of choice, and each near-relation costs about `W/4`. So a member with `W = O(M log M)` and
`r = O(M)` can carry at most `O(M / log M)` independent near-relations (*heuristic count*). In
that regime Corollary 3.4 forces `|F| >= (1/2 - o(1)) M`.

**What this means for the barrier programme.** Route G was meant to produce fake circles with
genuine residues, hence barriers for P3-all. Theorem 3.1 shows that on any prime set with a
complete closure there are no fakes in `P3-G` whatsoever: every member is an actual cluster.
Theorem 3.2 shows that a non-actual member must keep about `M/2` columns unpinned (or pay with
genuine coincidences). Whether a prime set is pinned is a question about smooth values of the
integers `D_s(A_v)` for low-weight `v`. Sections 4 and 5 show that prime sets of polynomial size
are pinned, by exact computation and in a probability model.

---------------------------------------------------------------------------------------------

## 4. Actual primes: verified collapse, densities and the pinning threshold

### 4.1 The computation

`gtight.py` computes, for Gaussian primes `pi_j = a_j + i b_j` (`b_j` even, random conjugates) and a
level `X`, every tight relation of weight `<= 2`:
* singles `v = e_j` and doubles `v = 2 e_j`, units `s = 0..3`;
* pairs `v = e_i + tau e_j`, `i < j`, `tau = +-1`, units `s = 0..3`.

Tightness is the exact test `|D_s(A_v)|` odd and dividing `L'_X` (vectorised trial division by the
primes `q <= X`, with the caps `A_q`). The closure keeps a signed union-find. It encodes the
generators `e_k - sign(k) e_root(k)`, `2 e_R` (root "half") and `e_R` (root "exact"). In each pass
it looks for the certificates of Lemma 2.4 with `k <= 2`:
* *self*: `2 sigma v in Pi`, `sin^2 d > 1/2`;
* *pair*: `sigma v + sigma' v' in Pi`, `sin^2 d + sin^2 d' > 1`.

It stops when a pass adds nothing.

`verify_certs.py` re-checks every certificate with independent code. It does not use the
union-find:
* it recomputes `A_v` in exact Gaussian arithmetic, `D_s`, and tightness by `L'_X mod |D_s| = 0`;
* it recomputes `sigma` from a 60-digit `mpmath` value of `arg A_v`, and `sin^2 d` exactly;
* it checks the circuit vector's membership in the lattice generated so far by an incremental
  Hermite normal form over `Z`;
* finally it checks `2 e_k` (and `e_k`) for every column.

`lemma_checks.py` brute-forces Lemmas 2.3 and 2.4 and Lemma 5.1 below in floating point (0
mismatches in 20000, 4956 and 100000 trials; evidence).

### 4.2 Verified collapse instances

**Theorem 4.1 (computer-verified).** For each prime set below, the computed sub-closure contains
`2 e_k` (indeed `e_k`) for every column `k`. Hence, by Theorem 3.1, every member of
`P3-G_(X')(M', C)` with that prime set, for any `X' >= X`, any `M'`, `C`, exponents and rows, is an
actual cluster after unit adjustment of the primes.

The prime sets are produced by `collapse_cases.py`: `r` random primes `p = 1 mod 4` in `[P, 2P]`
with `P = floor(M^K)` (generator `random.Random(seed)`), random conjugates, and `X = 2M`.

| M | b | K | r = b(M-1) | X | tight singles / pairs | pair degree | certificates | verified |
|---|---|---|---|---|---|---|---|---|
| 44 | 2 | 3 | 86 | 88 | 149 / 1235 | 28.7 | 86 | all 86 exact |
| 84 | 8 | 3 | 664 | 168 | 1130 / 76601 | 230.7 | 667 | all 664 exact |
| 164 | 8 | 3 | 1304 | 328 | 2288 / 294621 | 451.9 | 1308 | all 1304 exact |
| 164 | 8 | 4 | 1304 | 328 | 1240 / 38391 | 58.9 | 1304 | all 1304 exact |
| 324 | 2 | 4 | 646 | 648 | 638 / 9911 | 30.7 | 646 | all 646 exact |
| 644 | 1 | 4 | 643 | 1288 | 589 / 9131 | 28.4 | 643 | all 643 exact |

("Pair degree" is twice the number of tight pair relations divided by `r`.) Since the closure does
not depend on the profile, Theorem 4.1 covers, for instance, every `b`-copy Paley profile with
`b(M-1) <= r` built on these primes, and every round-5/6 construction transplanted onto them.

Not every prime set collapses (`collapse_cases_out.jsonl`):

| prime set | X | pair degree | pinned / r |
|---|---|---|---|
| `sharp.md` Thm 5.1 primes at `M = 44`, `b = 33` (1419 consecutive primes `= 1 mod 4` from `Z = 605521919`) | 88 | 3.2 | 2 / 1419 |
| same | 132 | 12.7 | 1272 / 1419 |
| random, `M = 644`, `b = 1`, `K = 4.25`, seed 1 | 1288 | 16.0 | 636 / 643 |
| random, `M = 1284`, `b = 1`, `K = 4.5` | 2568 | 18.1 | 1279 / 1283 |
| `p3data.md` window primes (`~8.4 * 10^8`), `M = 44`, `r = 388` | 88 | 0.7 | 0 / 388 |
| same, `M = 84`, `r = 748` | 168 | 9.7 | 18 / 748 |

The partial cases are not collapses: a few free columns remain. Corollary 3.4 then still requires
the profile to see at least `rank V - k` of them on its character lattice (or `k + 1`
near-relations). The window primes of `p3data.md` at `M <= 84` are far above the threshold below
(`K = 5.4` and `4.6`), and their closures are essentially empty.

### 4.3 Tight densities follow the divisor count of `L'_X` (evidence)

*Model prediction.* A relation with value `D` is tight with probability
`tau_X(|D|) := #{d | L'_X : d in [|D|/2, 2|D|]} / #{odd n in [|D|/2, 2|D|]}`. The divisors of `L'_X` are
enumerated exactly. Data: 400 primes in `[M^K, 2M^K]`, all singles, `2 * 10^5` random pairs per
sign, units `s = 1, 2, 3` (`s = 0` is never tight, Remark 1.5). Source: `density.py`,
`density_table_out.txt`.

| M | K | X | singles: observed / predicted | rate | pairs: observed / predicted | rate |
|---|---|---|---|---|---|---|
| 44 | 3 | 88 | 682 / 664 | 0.568 | 66633 / 66296 | 0.0557 |
| 44 | 4 | 88 | 355 / 370 | 0.296 | 8602 / 8721 | 0.0072 |
| 84 | 3 | 168 | 674 / 672 | 0.562 | 69752 / 68666 | 0.0583 |
| 84 | 4 | 168 | 365 / 376 | 0.304 | 8838 / 8991 | 0.0074 |
| 164 | 3 | 328 | 714 / 688 | 0.595 | 69017 / 68391 | 0.0577 |
| 164 | 4 | 328 | 389 / 376 | 0.324 | 9154 / 8918 | 0.0076 |

Two things stand out:
* The values `D_s(A_v)` of actual Gaussian primes are as smooth as random odd integers of the same
  size, to within sampling error (every row within 4%).
* At fixed `K` the densities do not depend on `M`. `tau_X(M^K)` tends to a positive constant (a
  Dickman-type value) as `M -> infinity`, while the number of partners `r` of each column grows.

### 4.4 The pinning threshold (evidence)

`threshold.py`: random primes in `[M^K, 2M^K]`, `b = 1` (`r = M - 1`, the fewest columns any
member with `rank V = M - 1` can have), `X = 2M`, three seeds (two at `M = 1284`). `K*` is the
interpolated exponent at which half the columns stay unpinned (`threshold_table_out.txt`).

| M | X | K* | lambda* = K* log M | H | lambda*/H | pair degree at K* |
|---|---|---|---|---|---|---|
| 84 | 168 | 3.65 | 16.2 | 13.2 | 1.23 | 8.0 |
| 164 | 328 | 4.03 | 20.6 | 16.0 | 1.29 | 6.7 |
| 324 | 648 | 4.29 | 24.8 | 19.0 | 1.30 | 7.9 |
| 644 | 1288 | 4.61 | 29.8 | 22.4 | 1.33 | 7.0 |
| 1284 | 2568 | 4.86 | 34.8 | 26.0 | 1.34 | 7.6 |

* *The transition is sharp.* At `M = 644`, for example, the free fraction is about 0.01 at
  `K = 4.25`, 0.14 at `4.5` and 0.99 at `4.75` (`threshold_fine.txt`).
* *It sits at a fixed pair degree,* about 7-8, for every `M`. That is, the tight pair density at
  the threshold is about `1/M`.
* *`K*` grows with `M`.* No fixed exponent `K` escapes pinning for all `M`. A fixed `K` gives a
  constant tight density (§4.3), hence pair degree `~ const * M`.

*Calibrated degree model* (`extrapolate.py`). Fit the measured pair degree by
`8 (M-1) (3/4) tau_X(c M^K)`, with `tau_X` from an exact-count log-grid DP of the divisors of `L'_X`.
One fitted constant (`log c = -0.30`) gives rms log-error 0.060 over the 110 runs with pair degree
`>= 0.5`. Setting the degree to 7.5 reproduces the measured `K*` to within 0.04
(`extrapolate_out.txt`).

We do **not** extrapolate the constant `lambda*/H` beyond the data. In the measured range,
single-prime seeds are abundant (0.4-1.1 tight singles per column at `K*`, all units together). §5.3 explains why
the transition degree must grow once they thin out. The data support:
* the order `lambda* = Theta(H_X)` (§5);
* the qualitative fact that `lambda*/log M` grows: no fixed exponent `K` escapes.

They cannot separate `H` from other slowly growing laws: `K*/loglog M` is also nearly constant,
2.45-2.47.

*Higher weight helps.* With the weight-3 relations `e_i +- e_j +- e_k` added (`triples_test.py`,
`M = 164`, `r = 163`): at `K = 4` the free columns drop from 41 to 4 and from 37 to 1 (two seeds); at
`K >= 4.25` nothing changes. Asymptotically, weight-`n` relations are heuristically more numerous
than pairs (§5.3), so the true threshold is above the pair-only one.

---------------------------------------------------------------------------------------------

## 5. A model theorem: random prime sets of polynomial size are pinned

### 5.1 The model

**Model `T(f; rho_1, rho_2)`.**
* There are `f` columns, with angles `theta_1, ..., theta_f` independent and uniform on `R/2 pi Z`.
* The relations are singles `(e_j, s)` and pairs `(e_i + tau e_j, s)`, `i < j`, `tau = +-1`, for units
  `s in {1, 2, 3}`. Unit `0` is never tight: this is Remark 1.5 with every `pi_j` normalised to `b_j`
  even.
* Each relation has its genuine angle `<v, theta>`, hence its orientation `sigma` and distance `d`
  (Lemma 2.3).
* Conditionally on `theta`, the tightness indicators are independent. A single is tight with
  probability `rho_1`, a pair with probability `rho_2`.
* `Pi*` is defined from these data exactly as in Definition 2.5.

*Relation to actual primes (heuristic).* For primes of size about `P`, §4.3 suggests
`rho_1 ~ tau_X(P^(1/2))` and `rho_2 ~ tau_X(P)`. The independence of tightness from the angles,
and between relations, is the modelling assumption. It is not used in §§1-4.

**Lemma 5.1 (two angles in one half).** Work modulo `pi`. Let `H_1 = (pi/4, 3pi/4)` and
`H_2 = (3pi/4, 5pi/4)`, and let `beta < beta'` lie in the same `H_i`. Then the unit pair
`(s, s') = (3, 1)` for `H_1`, or `(1, 3)` for `H_2`, gives orientations and distances with

```text
sigma(beta, s) = -sigma(beta', s')     and     d(beta, s) + d(beta', s') = pi/2 + (beta' - beta)  >  pi/2,
```

where `sigma(beta, s)` and `d(beta, s)` are those of Lemma 2.3 for `alpha = beta - s pi/4`.

*Proof.*
* `H_1`, `s = 3` at `beta`: `beta - 3pi/4 in (-pi/2, 0)`, so `sigma = -1` and `d = 3pi/4 - beta`.
  `s' = 1` at `beta'`: `beta' - pi/4 in (0, pi/2)`, so `sigma = +1` and `d = beta' - pi/4`.
* `H_2`, `s = 1` at `beta`: `beta - pi/4 in (pi/2, pi)`, so `sigma = -1` and `d = 5pi/4 - beta`.
  `s' = 3` at `beta'`: `beta' - 3pi/4 in (0, pi/2)`, so `sigma = +1` and `d = beta' - 3pi/4`. QED

**Lemma 5.2 (an ordering bound).** Let `m >= 1` slots be ordered by angle. Each slot `a` carries two
indicators `T_lo(a)`, `T_hi(a)`, and all `2m` indicators are independent with mean `rho`. Then

```text
P( no a < b with T_lo(a) = T_hi(b) = 1 )  =  (1 - rho)^(m-1) (1 + (m-1) rho)  <=  (1 + m rho) e^(-rho (m-1)),
```

and this is nonincreasing in `m`.

*Proof.* Let `a*` be the first slot with `T_lo = 1`. The event says that `T_hi = 0` after `a*`. So
its probability is `(1-rho)^m + sum_(a=1)^m (1-rho)^(a-1) rho (1-rho)^(m-a)`. The ratio of
consecutive values is `(1 - rho)(1 + rho/(1 + (m-1) rho)) <= 1 - rho^2`. QED

### 5.2 The theorem

**Theorem 5.3 (random prime sets are pinned, in the model).** In `T(f; rho_1, rho_2)` with `f >= 3`
and `rho_2 <= 1/2`, if

```text
(f - 1) rho_1^2 rho_2  >=  48 log f,
```

then with probability at least `1 - (5 log f + 2)/f^2` every column is exactly pinned:
`e_j in Pi*` for all `j`.

*Proof.*
1. *Seeds.* Let `A` be the first `floor(f/2)` columns. Call `j in A` a seed if one of these holds:
   * `theta_j mod pi in (3pi/4, pi)` and the singles `(e_j, 1)`, `(e_j, 2)` are tight;
   * `theta_j mod pi in (0, pi/4)` and `(e_j, 2)`, `(e_j, 3)` are tight.

   In the first case unit 2 has `sigma = +1`, `d = theta - pi/2`, and unit 1 has `sigma = -1`,
   `d = 5pi/4 - theta`. In the second, unit 3 has `sigma = +1`, `d = theta + pi/4`, and unit 2 has
   `sigma = -1`, `d = pi/2 - theta`. Either way the two relations have the same `v = e_j`,
   opposite `sigma`, and `d + d' = 3pi/4 > pi/2`. That is a certificate (2.2) relative to `{0}`, so
   `e_j in Pi*`.

   Seeds occur independently over `j in A`, each with probability `rho_1^2/2`. Let `E_0` be the
   set of seeds, and `mu = floor(f/2) rho_1^2/2 >= (f-1) rho_1^2/4`. Chernoff gives
   `P(|E_0| < mu/2) <= e^(-mu/8) <= f^(-3)`, since `(f-1) rho_1^2 >= 96 log f` (the hypothesis
   together with `rho_2 <= 1/2`).
2. *Propagation.* Fix `j not in E_0`. Condition on all angles and on the singles; `E_0` is then
   fixed. The pair indicators remain independent with mean `rho_2`.
   * For `i in E_0`, `e_i in Pi*`, so the pair relation `v_i = e_j - e_i` has class
     `[v_i] = [e_j]` modulo `Pi*`. Its angle is `beta_i = theta_j - theta_i mod pi`; these are distinct
     almost surely.
   * Some half `H_1` or `H_2` contains `m >= |E_0|/2` of the `beta_i`.
   * In that half, take `T_lo(i)` and `T_hi(i)` to be the tightness of `(v_i, s)` and `(v_i, s')`,
     with `(s, s')` from Lemma 5.1. If `T_lo(a) = T_hi(b) = 1` for some `beta_a < beta_b`, then
     `sigma_a [v_a] + sigma_b [v_b] = (sigma_a + sigma_b)[e_j] = 0` and the arcs sum to less than
     `pi` (Lemma 5.1). That is a certificate relative to the group generated by `E_0`, so
     `e_j in Pi*`.
   * By Lemma 5.2, `P(e_j not in Pi* | angles, singles) <= (1 + x) e^(rho_2 - x)`, where
     `x = rho_2 |E_0|/2`.
3. *Union bound.* On `|E_0| >= mu/2` we have `x >= rho_2 (f-1) rho_1^2/16 >= 3 log f`.
   `(1 + x) e^(-x)` decreases for `x >= 0`, and `e^(rho_2) <= e^(1/2)`. So the total failure
   probability is at most `f^(-3) + f e^(1/2) (1 + 3 log f) f^(-3) <= (5 log f + 2)/f^2`. QED

The closure is monotone in the set of tight relations (it is the least fixed point of a monotone
rule), so the theorem also holds when the tightness probabilities are only bounded below by
`rho_1`, `rho_2`.

**Corollary 5.4 (in the model: polynomial primes cannot stay free).** Take `X = M^(1+o(1))` and
assume the model with `rho_1 >= tau_X(P^(1/2))`, `rho_2 >= tau_X(P)` for every set of columns
whose primes are `<= P`. Fix `eps > 0` and put `lambda_eps = (1/2 - eps) H_X`.
* With probability `1 - o(1)`, at most `M^(1 - eps)` columns with `log p_j <= lambda_eps` are
  unpinned, for every non-adaptive choice of the prime sizes.
* In particular, every squarefree member with at most `O(M/log M)` independent genuine
  near-relations (Remark 3.5; actual clusters have `M - 1`) has

  ```text
  W  >=  sum_(j in F) log p_j  >=  (1/4 - eps - o(1)) M H_X.
  ```

*Proof.* Two facts about smooth numbers are needed.
* `tau_X(X^u) >= rho_D(u) e^(-o(u log u))` uniformly for `1 <= u <= log X`.
  * This follows from Hildebrand's `Psi(x, y) = x rho(u)(1 + O(log(u+1)/log y))`, valid for
    `y >= exp((loglog x)^(5/3 + eps))`, which covers this range. Apply it to the window
    `[Y/2, 2Y]`; shifting `u` by `O(1/log X)` changes `rho` by a factor `1 + o(1)`.
  * Restricting to odd divisors of `L'_X` removes the smooth integers divisible by some
    `q^(A_q + 1) > X`. Since `Psi(Y/d, X)/Psi(Y, X) <= d^(-1) u^(O(log d/log X))`, that proportion
    is `<= sum_(q <= sqrt X) u^2/X + sum_(q > sqrt X) u^2 q^(-2) = o(1)` in this range of `u`.
  * The odd restriction changes the density by a bounded factor.
* `rho_D(u) = exp(-u log u (1 + o(1)))`.

Let `S` be the set of columns with `log p_j <= lambda_eps`, and `P = e^(lambda_eps)`. Every pair
relation inside `S` has value `O(P)`, and every single has value `O(P^(1/2))`. Write `u = lambda_eps/log X = (1/2 - eps) log M/loglog M`. Then

```text
rho_1^2 rho_2  >=  rho_D(u/2)^2 rho_D(u) e^(-o(u log u))  =  exp(-2 u log u (1 + o(1)))  =  M^(-1 + 2 eps + o(1)).
```

If `|S| >= M^(1-eps)`, Theorem 5.3 applied to `S` (the closure of `S` is contained in that of all
columns) pins all of `S` with probability `1 - O(log M/M^(2-2eps))`. The free columns then have
`log p_j > lambda_eps`, except at most `M^(1-eps)` of them. Corollary 3.4 with `k = O(M/log M)`
gives `|F| >= (1/2 - o(1)) M`. QED

The proved model threshold `(1/2) H_X` is conservative. It uses only single seeds and one round of
propagation through two specific units.

### 5.3 Heuristics beyond the theorem

All of the following assume the model, i.e. tight densities `~ rho_D(u)` and independence. Write
`log P = 2 kappa H_X`. Then `rho_2 ~ M^(-2 kappa)` (pairs) and `rho_1 ~ M^(-kappa)` (singles).

1. *Pair-only closure.* The pair degree is `~ M^(1 - 2 kappa)`, and single seeds number
   `~ M^(1 - 2kappa)`. Pinning through pairs is a 2-neighbour bootstrap percolation with edge
   density `~ rho_2` started from the seeds. It percolates once the seed set exceeds
   `~ 1/(M rho_2^2) = M^(4kappa - 1)` (the Janson-Luczak-Turova-Vallier threshold for `G(n, p)`),
   i.e. for `kappa < 1/3`.
   * In the measured range (§4.4) seeds are abundant, and the transition sits at a constant
     degree, `kappa ~ 0.55-0.6`.
   * Asymptotically the heuristic pair-only threshold is `lambda* ~ (2/3) H_X`.
2. *Higher weight.* Once an exactly pinned set `E` of size `Theta(M)` exists, consider relations
   `e_j + (terms in E)` of weight `n`. Their class is `[e_j]`, so they act as singles on `j`. Per
   column there are `~ M^(n-1) rho_D(n u/2) ~ M^(n - 1 - n kappa)` of them, which is many for
   `kappa < 1 - 1/n`.
   * So, with seeds, weight-`n` relations pin every column of size `lambda < (2 - 2/n) H_X`.
   * Without seeds, a circuit (2.2) among unpinned columns needs every column of its support in at
     least two relations. Counting gives `~ M^(sum n_i (1/2 - kappa))` such circuits, which tends
     to zero for `kappa > 1/2`.
   * With seeds, there are many relations into `E`; but the seeds themselves (`~ M^(1-2kappa)`)
     only percolate for `kappa < 1/3` (item 1).
   * For `kappa < 1/2`, a cycle of `L` pair relations is a candidate circuit. There are
     `~ M^(L(1 - 2kappa)) c^L/L!` such cycles with consistent orientations and arc sum `< pi`. So for
     `L > 1/(1 - 2kappa)` every column lies on many pinning circuits.
   * Heuristically, then, the full closure (circuits of all lengths and weights) pins prime sets
     with `log p < (1 - eps) H_X`, and pins nothing when all `log p > (1 + eps) H_X`. So
     `lambda* ~ H_X`, in rough agreement with the measured `1.06-1.22 H_X` (§4.4).
   * The computed closure (circuits with `k <= 2`) has the lower heuristic threshold `(2/3) H_X`
     of item 1.
3. *Adaptive selection.* A builder may pick primes from the pool of `~ P/log P` candidates to avoid
   certificates. With `C1`-type certificates only, this is at best an independent set in a random
   graph of density `~ 6 rho_2` on the pool, of size `~ 2 log P/(6 rho_2) = lambda/(3 rho_2)`.
   Holding `(1/2 - o(1)) M` free columns needs `rho_2 <~ lambda/M`, i.e.
   `u log u >= log M - O(loglog M)`. So selection changes `lambda*` only by a factor `1 + o(1)`.
   This is a count, not a proof: actual primes are not independent random samples.

So the heuristic critical size of a free column is `log p ~ c H_X`, with `c` between about 1 (pinning
side) and 2 (union-bound side, §6.2). For `X = M^(1+o(1))` that is `Theta(log^2 M/loglog M)`.

---------------------------------------------------------------------------------------------

## 6. Upper bounds, and the answer to task (i)

**Proposition 6.1 (fakes in P3-G at `W = O(M X)`, rigorous).** Fix `C > 0`. There is `M_0(C)` such
that for all `M >= M_0` and all `X >= 2M`, `P3-G_X(M, C)` has a member with
`e = 1`, all angles different from the genuine ones, and

```text
W  <=  33 (M' - 1) ( 2 psi(X) + 2 log(64 * 33 M') + 2 )  =  (68.7 + o(1)) M X,
```

where `M' = (1 + o(1)) M` is the least prime Paley order `>= max(M, 44)` and `psi` is Chebyshev's
function (`psi(X) < 1.04 X`, Rosser-Schoenfeld).

*Proof.* First let `M = M'` be a prime Paley order. Take the profile of `sharp.md` §4 with `b = 33`,
`r = 33(M-1)`, `e = 1`. Let `L = prod_(q <= X odd) q^(A_q) <= e^(psi(X))`.
1. *Primes.* Put `Z = (64 r L)^2`. Then `[Z, 2Z]` contains more than `r^2` primes `= 1 mod 4`, so
   some interval of length `Z/r` contains `r` of them (`sharp.md` Thm 5.1, step 1). Take any
   Gaussian primes `pi_j` above them, and their genuine residues. Then `p_j > X`, so `L'_X = L`,
   `w_j in [tau, tau + 1/r]`, and `tau <= log(2Z)`.
2. *Moduli.* By Lemma 1.2, every collision modulus (characters and monomials, every unit) divides
   `L`.
3. *Characters, and monomials in `L`.* Use `k = 0` and `delta` uniform on `[0, Delta]^M`. For a
   character `c`, put `eps_c = L e^(-w(v(c))/2)`, and let `Fail_c` be the event
   `dist(c . delta/2, (pi/4) Z) < arcsin(min(1, eps_c))`.
   * *Sufficiency.* Off `Fail_c`, (1.1) of `sharp.md` holds for every unit `eta`, since
     `|sin((c . delta - arg eta)/2)| >= sin dist(c . delta/2, (pi/4) Z)`. (1.2) holds for
     `v = v(c)` and every `s`, since `<v(c), phi> = c . delta/2`. Lemma 0.1 also holds, since the
     moduli are `>= 1`.
   * *All monomials in `L` are characters.* `L = {v(c)}` (`sharp.md` Lemma 4.1).
   * *Union bound.* By `sharp.md` Lemma 5.1a,
     `P(Fail_c) <= 4 eps_c (n + pi/Delta) <= 4 L (n + pi/C) e^(-s_c/2)`.
   * *Slacks.* Use `sharp.md` Prop. 4.3:
     * pairs: `s_c >= 14.5 tau - 1/2`;
     * ternary non-pairs: at most `3^M` of them, `s_c >= tau M/32 - 1/2`;
     * amplitude `h >= 2`: at most `(2h+1)^M` of them, `s_c >= 7.25 tau h M - 1/2`.
   * *Total.* Since `tau >= log Z >= 2 log L`, the sum of the bounds is at most
     `4 e^(1/4) (M^2 (2 + pi/C) L^(-13.5) + 3^M (M + pi/C) L e^(-tau M/64) + sum_(h >= 2) (2h+1)^M (hM + pi/C) L e^(-3.6 tau h M))`.
     This is `< 1/2` for `M >= M_0(C)`. Indeed `L >= prod_(3 <= q <= 2M) q = e^((2 - o(1)) M)`
     (as `X >= 2M`), so `L e^(-tau M/64) <= L^(1 - M/32)` and `L e^(-3.6 tau h M) <= L^(1 - 7.2 h M)`
     beat the counts `3^M` and `(2h+1)^M`. Fix `delta` off every `Fail_c`.
4. *Monomials outside `L`.* For `v` outside `L`, `<v, phi>` is Haar-uniform on the realisation torus
   (`outside.md` Thm A.1(i); `rank S = M`, so every `delta` is realisable). The failure set of (1.2)
   for unit `s` has measure `(2/pi) arcsin(min(1, nu_s m_(v,s) e^(-w/2))) <= L e^(-w(v)/2)`. Summing
   over `s` and `v != 0`:

   ```text
   Haar(bad)  <=  4 L sum_(v != 0) e^(-w(v)/2)  =  4 L ( prod_j (1 + 2/(sqrt(p_j) - 1)) - 1 )  <=  4 L (e^(2r/(sqrt Z - 1)) - 1)  <  1/4.
   ```

   (1.2) with `m >= 1` implies Lemma 0.1. So a positive-measure set of `phi` gives a member. Almost
   every such `phi` differs from `theta` (indeed from every vector of genuine angles).
5. *Weight.* `W <= r tau + 1 <= r (log 2 + 2 log L + 2 log(64 r)) + 1`.
6. *General `M`.* Restrict a member at `M'` to `M` rows and delete the columns constant on them
   (`sharp.md` Lemma 5.1b). Genuine residues stay genuine, and moduli can only drop. QED

The proposition shows that `P3-G` does contain non-actual members (fakes). Its rate `M X >= 2 M^2`
is far above the round-5/6 fakes. It is the price of making completion deterministic (`m <= L`).

**Heuristic 6.2 (fakes at `W = O(M H_X)`).** Take the same profile, with primes of size
`log p = (2 + eps) H_X` (`u = (2 + eps) log M/loglog M` in units of `log X`). Model the moduli of
§4.3 by `m_(v,s) ~ gcd(D, L'_X)` for a random odd `D` of size `e^(w(v)/2)`.
* *Completion.* The expected completion sum over weight `n = ||v||_1` is
  `~ (4r)^n tau_X(P^(n/2)) ~ (4r)^n M^(-n(1 + eps/2) + o(n))`. Summed over `n` this is `o(1)`.
* *Pairs.* The genuine pair demands have maximum `~ 4H` (`p3data.md` §2.1 measures `D_max/H` of
  3.2-3.7 at `X = 2M`). That is below the uniform margin `14.5 (2 + eps) H_X`.
* *Non-pair characters.* They have slack `Theta(M H_X)`.

So `W ~ 33(M-1)(2 + eps) H_X = (66 + 33 eps + o(1)) M H_X`, the same order as `sharp.md` Theorem
6.4 for P3-all with fake residues. Making this rigorous needs smooth-value statistics for the
integers `D_s(A_v)`, uniformly over all low-weight `v`. For Gaussian primes of size `e^(Theta(H))`
these are out of reach of current equidistribution results (Siegel-Walfisz range only).

**Answer to task (i) (completion with genuine data at polynomial prime size).**
* *The completion sum is not small.* For primes of size `M^K`, its tight weight-2 part already has
  `~ 3 r^2 tau_X(M^K)` terms, each `(2/pi) d_e` with `d_e` of order one. That is a positive constant
  times `r^2` (§4.3 measures `tau_X(M^K)`). Exactly, on two verified sets of Theorem 4.1
  (`completion_sum.py`), the slab measures of the tight relations of weight `<= 2` alone sum to
  15190 (`M = 164`, `b = 8`, `K = 4`) and 4052 (`M = 644`, `b = 1`, `K = 4`), against the budget
  `< 1` of a union bound.
* *Worse, completion is impossible.* It is not just the union bound that fails. By Theorems 3.1-3.2
  the angles are pinned to the genuine ones, so no non-actual member exists at all on such prime
  sets: rigorously on the verified sets of Theorem 4.1, and with probability `1 - o(1)` in the
  model (Theorem 5.3).
* *Needed prime size.* A non-actual member needs about `M/2` columns with `log p >= lambda*`, where:
  * `lambda* >= (1/2 - eps) H_X` (model, proved);
  * `lambda* ~ H_X` (heuristic, §5.3);
  * `lambda* = (1.06-1.22) H_X` (measured, `M <= 1284`, pair-only closure).

  Sufficient: `log p >= (2 + eps) H_X` (heuristic), or `log p >= 2.08 X + O(log M)` (proved,
  Prop. 6.1).

---------------------------------------------------------------------------------------------

## 7. Task (ii): the pair channel is not the binding constraint in P3-G

The task asked for an adaptive profile with `r = O(M)` prime columns that absorbs the genuine pair
demands at `W = O(M log M)`, with all character margins. In `P3-G` (every monomial strengthened)
this question is pre-empted by §§3-5. At `W = O(M log M)` most columns have polynomial size, so
typical prime sets collapse, whatever the profile does with the pairs.

* *Repair columns are counterproductive.* A repair column with a small prime (`log p = O(log M)`)
  is itself pinned in typical prime sets. Worse, each pinned column acts as a seed for its free
  neighbours. A free column `j` of size `e^lambda` has `~ |P| tau_X(e^(lambda/2 + O(log M)))`
  tight relations into the pinned set `P`, which is many unless `lambda >~ 2 H_X`. So pinned
  repair columns raise the size the free columns need (heuristic, as in §5.3).
* *At the scale where fakes can exist, pairs are slack.* At `log p ~ c H_X` per free column, the
  uniform Paley margin `(b - 4) c H_X/2` exceeds the heuristic maximal genuine demand `~ 4 H` once
  `b >= 4 + 8/c`. No adaptivity is needed (Heuristic 6.2).
* *Characters-only genuine class.* Let `P3-G-char` be `P3-G` with (e-all) dropped (genuine
  residues, characters strengthened, all other monomials only by Lemma 0.1). It does **not**
  collapse. The tight relations of §2 are low-weight monomials. In a member, every integral
  character has `w(v(c)) >= W/2 - 2 log(n C/2^(3/2))` (by (e-char) with `m >= 1`), so these
  monomials are not integral characters, and in `P3-G-char` they carry no residue strengthening. Completion there is `outside.md` Thm A.1(ii) (`p_j >= 21 r^2`).
  So for `P3-G-char` the task's question (ii) is exactly `p3data.md`'s question (Q) together with
  Hypothesis H, and stays open.
  * `p3data.md` Theorem 5.1 shows that every fixed profile fails with probability `1 - o(1)`
    against residues random relative to it, below `(2 - eps) M H`.
  * Whether prime selection can correlate genuine residues with an adaptive profile is not
    decided here.

---------------------------------------------------------------------------------------------

## 8. Consequences and open points

1. **Route G gives no barrier at `O(M log M)`.** On prime sets whose closure is complete, `P3-G` has
   no fake circles at all (Theorem 3.1). Prime sets of polynomial size are complete:
   * verified on the sets of Theorem 4.1, at `M` up to 644 and `r` up to 1304;
   * with probability `1 - o(1)` in the model (Theorem 5.3, Corollary 5.4).

   A non-actual member needs about `M/2` columns of size `log p >~ c H` (Corollary 3.4, §5.3). So
   route G cannot answer `sharp.md`'s P3-all\* positively. The genuine class sits at the same
   heuristic rate `Theta(M log^2 M/loglog M)` as `sharp.md` Theorem 6.4, by a different mechanism:
   saturation of the height condition at a planted point, instead of Poisson tails of independent
   fibres.
2. **Nor does it give a lower bound for free.** `P3-G` contains every actual cluster
   (Prop. 1.3), and on complete prime sets it contains nothing else. A lower bound for `P3-G` at
   scale `O(M log M)` is therefore, on those prime sets, a lower bound for actual circles. The
   genuine residue data plus every monomial gap determine the angles up to units: they are not a
   relaxation there.
3. **The mechanism is new and exact.** The per-unit height condition is attained with equality,
   `m_(v,s) = |D_s(A_v)|`, on a positive proportion of low-weight monomials (§4.3). Such a monomial
   turns the slab condition into a one-sided constraint through the genuine angle (Lemma 2.3).
   Circuits of one-sided constraints pin (Lemma 2.4). The certificates are exact integer data, and
   every collapse claim here is re-checked by independent code.
4. **Open points.**
   * (a) *Unconditional density.* Lower bounds for the density of tight relations
     (`D_s(pi_i conj(pi_j))` dividing `L'_X`) for random Gaussian primes of size `M^K`, of the
     order `tau_X(M^K)`. With them, Corollary 5.4 holds for random prime sets without the model.
     This is a smooth-values-of-bilinear-forms-at-primes problem.
   * (b) *Adaptive selection.* Can primes be selected (not drawn) so that `>= M/2` columns of
     size `M^(O(1))` stay unpinned? The independent-set count of §5.3 says no, up to a factor
     `1 + o(1)` in `lambda*`, but it is a count in a model.
   * (c) *The window `(lambda*, 2 H_X)`.* Unpinned is not the same as completable. Fakes are only
     known (heuristically) above the union-bound threshold `2 H_X`, and rigorously above `2.08 X`.
   * (d) *`P3-G-char`* is `p3data.md`'s (Q) + H (§7).

---------------------------------------------------------------------------------------------

## 9. Files (`round7/p3G_checks/`)

All scripts run from this directory with `python3 <script>` (`-W ignore` silences numpy overflow
warnings in the screening step, which never decides a certificate).

| file | content | type |
|---|---|---|
| `gtight.py` | Gaussian primes, tight relations (singles, doubles, pairs, triples; vectorised exact screening by trial division with caps), exact orientation and `sin^2 d`, the signed union-find closure with `k <= 2` certificates | library |
| `verify_certs.py` (-> `verify_certs_out.txt`) | independent re-check of every certificate: exact `Z[i]` products, `L'_X mod abs(D)`, 60-digit `mpmath` orientation, exact arc sums, lattice membership by incremental HNF over `Z`; final `2 e_k` / `e_k` membership | exact; ALL CERTIFICATES VERIFIED |
| `collapse_cases.py` -> `collapse_cases_out.jsonl`, `collapse_cases_log.txt` | Theorem 4.1 and the table of partial cases (incl. the `sharp.md` Thm 5.1 primes at `M = 44` and the `p3data.md` window primes) | computer-verified |
| `trend.py` | random prime sets and one closure run | library |
| `threshold.py` -> `threshold_out.jsonl`, `threshold_b1.txt`, `threshold_fine.txt`; `threshold_table.py` -> `threshold_table_out.txt` | §4.4 pinning threshold | evidence |
| `density.py` -> `density_out.txt`; `density_table.py` -> `density_table_out.txt` | §4.3 tight densities vs divisor counts of `L'_X` | evidence |
| `extrapolate.py` -> `extrapolate_out.txt` | §4.4 calibrated degree model (log-grid DP for the divisors of `L'_X`) | evidence |
| `triples_test.py` -> `triples_out.jsonl` | §4.4 effect of weight-3 relations (`M = 164`) | evidence (certificates verified) |
| `completion_sum.py` -> `completion_sum_out.txt` | §6, tight weight-`<= 2` part of the completion sum on two verified sets | exact up to float summation |
| `lemma_checks.py` -> `lemma_checks_out.txt` | brute-force checks of Lemmas 2.3, 2.4, 5.1 | evidence |
| `notes.md` | scratch notes | -- |
