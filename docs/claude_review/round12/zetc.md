# Round 12: the tangent-cone hypothesis ZE-TC(A,m) for T = 1 slab-boxes


## 0. Summary of results

**No growth improvement and no uniform bound.** `M <= (2/3 + eps) log R/loglog R` is still the best bound. The
qualitative tangent-cone hypothesis is **not decided in the regime**. What this round establishes about it, for
`T = 1` slab-boxes (`A` in the ternary cube cut by `|k(S_p)| <= L_p`):

1. **The cheap-span form of ZE-TC (README §8.2) is false; the correct form uses the tangent degree** (§1-2, proved).
   * `TC(A,m) = {u : nu_A(u) < m}`, where `nu_A(u)` is the largest `i` with `<k,u>^i` not in the leading-form ideal
     (Lemma 1.1). `nu_A` is subadditive, additive on products, bounded by `rho` of projections, bounded below by
     zonotopes, and controls the transversal order along *every* rational subspace `M` in `TC`: order `>= m - nu_A(M)`.
   * **One cut, exact (Theorem 2.1):** `nu(u) = min(2|supp u_S|, 2L + 2t(u_S), s + L) + 2|supp u_R|`. Hence
     `TC(A,m)` is a finite union of rational subspaces at every `m` (ZE-TC^rat holds), but for `s + L < m <= rho(A)`
     the pieces `R^S x R^Y` are **not** spans of cheap families: ZE-TC^cheap fails, also in LP-supercritical models
     (`n = 7`, one column `{1,2}` with `L = 1`, `m = 8`; exact check). The relevant cost of a rational subspace is
     `nu_A(M)`, not the sum of widths. In the regime at `m = ceil(4w)` this particular failure does not occur
     (balanced columns have `|S| >> m`).
   * Replacement: **ZE-TC^rat** (real `TC` in the union of rational subspaces `M` with `nu_A(M) < m`), with the weaker
     **ZE-TC^deg** (real `TC` consists of degenerate points), which is what the tangent-cone reduction needs.
2. **Few cuts never test ZE-TC** (§3).
   * Proposition 3.1: if `zeta(A) >= m` then one product certificate puts `TC(A,m)` in finitely many rational
     hyperplanes. One cut: `zeta = rho` (proved). Two disjoint cuts: product, ZE-TC^rat (proved). Two crossing or nested
     cuts with `|S_1 u S_2| <= 7`, and three cuts with `|S_1 u S_2 u S_3| <= 6`, `L <= 2`: `zeta = rho` (computer-assisted,
     exact in the direction used: 421 and 953 blocks). So ZE-TC^deg holds there at every order, for a trivial reason.
   * **An induction on cuts that propagates zonotopes cannot reach the regime**, where `zeta < m` (compute.md Prop 2.2).
3. **Section monotonicity (Proposition 3.3, proved):** ZE-TC^deg passes from any lattice section of `A` to `A`.
   Consequences: products of pieces that are TC-good at their top degree are TC-good; and the **Sylvester tower**
   (Theorem 4.2): ZE-TC^deg holds for `Syl(2^k, 2)` at `m = 8` and for `Syl(2^k, 3)` at `m = 11`, for every `k >= 3`
   (from J4 for `M8`, `M8'`). The `L = 2` statements are zonotope-free and overdetermined, but **LP-subcritical**. A count
   shows sections of this kind cannot certify the regime (they would need top degree `>= ~4x` their slab usage;
   `M8` has `~0.7x`).
4. **Sylvester-Hadamard symmetry (§4):** a nondegenerate point has trivial stabiliser (Lemma 4.3), so orbit-sum and
   fixed-subspace reductions only see degenerate points; symmetry helps the computations (isotypic Hilbert
   functions) but not the real-geometry question. Ball effect (Lemma 4.4): Hadamard boxes lie in `|k|_2 < w`, which
   makes small Hadamard boxes thin.
5. **Computations (§5).**
   * **No Hadamard `T = 1` model with `n <= 12` is LP-supercritical with `V(A,m) != 0`** (rigorous: exact or mod-`p`
     upper bounds on `rho`, which are valid upper bounds). Sylvester `n = 16`: `L = 3` excluded by an exact fibration
     bound (`rho <= 22 < 24`), `L >= 4` by the one-cut bound; `L = 2` excluded by a mod-`p` spanning certificate through the Singer cycle (`rho <= 15 < 16`; all 30
     isotypic components of full rank, row counts summing to `|A| = 213937`). **So the smallest LP-supercritical Sylvester
     `T = 1` model has `n >= 32`.**
   * Random, graph, design and few-cut models with `n <= 9` are zonotopal or LP-subcritical; the smallest zonotope-free
     twin-free model found is `Z6` (`n = 7`, exact `zeta = 6 < 7 = rho`, 3 tangent forms).
   * **Products** (e.g. `M8 x M8`, `n = 15`, `m = 16`, `1764` tangent forms) meet every stated criterion
     (zonotope-free, LP-supercritical, overdetermined, twin-free), and ZE-TC^deg holds there by reduction to the
     factors. So the criteria of the task must be supplemented by *indecomposability*; no indecomposable
     LP-supercritical zonotope-free model was found in computable range.
   * Paley-12 (`L = 2`, `m = 8`, `4345` tangent forms, zonotope-free, LP-subcritical): each of its 330 8-point sections
     has a single top form, and those checked exactly (3 of the 30 classes) are irreducible and indefinite, hence have
     nondegenerate real tangent-cone points; but a numerical search
     (**evidence**) finds common zeros of all 330 section octics only next to degenerate strata.
   * **No nondegenerate real common zero was found in any overdetermined zonotope-free model.**

| statement | status |
| --- | --- |
| ZE-TC^cheap (README §8.2) | **false** in general `T = 1` models (one cut, LP-supercritical, `n = 7`) |
| ZE-TC^rat, one cut / disjoint cuts | **proved** (all orders) |
| ZE-TC^deg, zonotopal models (incl. one cut; two cuts on `<= 7` coordinates; three cuts on `<= 6`, `L <= 2`) | **proved** (trivial: Prop 3.1) |
| ZE-TC^deg, `Syl(2^k,2)` at `m = 8` (zonotope-free), `Syl(2^k,3)` at `m = 11`, all `k >= 3` | **proved** from J4 (LP-subcritical) |
| ZE-TC^deg, products of TC-good pieces (e.g. `M8 x M8`, LP-supercritical) | **proved** |
| indecomposable LP-supercritical zonotope-free `T = 1` models | none with Hadamard profile at `n <= 12`; none with Sylvester profile at `n <= 16`; none found elsewhere in range |
| ZE-TC in the regime | **open** |

## 1. Setup and conventions

Labels as in rounds 10-11. *Theorem/Proposition/Lemma*: proved here with all quantifiers. *Exact*: integer or
exact-rational computation (python-flint `fmpq_mat`/`fmpz_mat`, sympy over `Q`). *Mod-p*: rank computations modulo
a prime `p ~ 2^61` (a mod-p rank is `<=` the rational rank, so mod-p computations give rigorous *upper* bounds on
`rho` and rigorous spanning statements, and only evidence in the other direction; §5.1). *Evidence*: floating point or finite experiment that no proof uses. Scripts: `zetc_checks/`
(`<script>.py` writes `<script>_out.txt`).

**Models.** `d = n - 1`, coordinates `k' = (k_1, ..., k_d)` of `k in Z^n_0` (base point `n`, `k_n = -sum k'`).
A sign column `beta in {+-1}^n` is recorded by its pattern `S = {x : beta_x != beta_n}`, a nonempty subset of
`[d]`; then `<k, beta> = +-2 k(S)`, `k(S) = sum_{x in S} k_x`, and the slab `|<k,beta>| <= w` is `|k(S)| <= L`,
`L = floor(w/2)`. (The base-point pattern is `S = [d]`, `k([d]) = -k_n`.) The `T = 1` slab-box is

```text
A = A(S_1, L_1; ...; S_r, L_r) = { k' in {-1,0,1}^d : |k(S_p)| <= L_p  (p = 1..r) },
```

with possibly different `L_p` (one `L` in the regime). `h_A(beta_p) = 2 width_A(k(S_p))`, so **LP-supercritical**
(`h_A(beta_p) <= m/2` for every column) means `width_A(k(S_p)) <= m/4` for every `p`; when the slab is attained,
`width_A(k(S_p)) = 2 L_p` and the condition reads `m >= 8 L_p` (`= 4w`). `I(A)` is the vanishing ideal,
`J = gr I(A)` its ideal of leading forms in `S = Q[k_1..k_d]`, `rho(A)` the top degree of `S/J`
(`V(A,m) != 0` iff `m <= rho(A)`), and `TC(A,m) = {u in C^d : Phi_c(u) = 0 for all c in V(A,m)}`.
`zeta(A)` is the largest `sum L_i` over lattice zonotopes `k_0 + sum_i [0,L_i] v_i` in `A` (compute.md §2).

Since `A` is contained in the ternary cube, `I(A) = (k_x^3 - k_x : x <= d) + (f_p : p <= r)` with
`f_p = prod_{|t| <= L_p} (k(S_p) - t)`: the cube ideal is the vanishing ideal of the `3^d` points, its coordinate
ring is a product of fields, and in a product of fields every ideal is the vanishing ideal of a subset. So `J`
contains `k_x^3` and `k(S_p)^(2L_p + 1)`, but in general `J` is strictly larger than the ideal they generate.

### 1.1 The tangent degree `nu_A`

**Definition.** For `u in C^d` let `l_u(k) = <k', u>` and

```text
nu_A(u) = max{ i >= 0 : l_u^i is not in J }      ( = max{ i : the function l_u^i on A is not a polynomial of degree < i } ).
```

For a rational subspace `M` of `R^d` put `nu_A(M) = max_{u in M} nu_A(u)`. Call `M` **tangent-cheap for `(A,m)`** if
`nu_A(M) < m`.

**Lemma 1.1 (basic properties).** Let `A` be finite, `u, u' in C^d`.
1. *(Level sets.)* `l_u^i in J` implies `l_u^(i+1) in J`; hence `{i : l_u^i notin J} = [0, nu_A(u)]`,
   `nu_A(u) <= rho(A)`, and `TC(A,m) = {u : nu_A(u) < m}`.
2. *(Semicontinuity.)* `{u : nu_A(u) >= i}` is Zariski open in `C^d`. So on an irreducible variety (e.g. a rational
   subspace `M`, over `C`) `nu_A` takes its maximum `nu_A(M)` on a dense open set.
3. *(Subadditivity.)* `nu_A(u + u') <= nu_A(u) + nu_A(u')`.
4. *(Projections.)* If `u` lies in the complex span of integer vectors `v_1, ..., v_j` and
   `pi(k) = (<k,v_1>, ..., <k,v_j>)`, then `nu_A(u) <= rho(pi(A))`. In particular `nu_A(t v) <= h_A(v)` for an integer
   direction `v`, and for a cheap family `F` (`sum h_A(v_j) < m`) the span `L_F` is tangent-cheap.
5. *(Products.)* If `A = A_1 x A_2` in disjoint sets of coordinates, `u = (u_1, u_2)`, then
   `nu_A(u) = nu_(A_1)(u_1) + nu_(A_2)(u_2)`.
6. *(Zonotope lower bound.)* If `k_0 + sum_i t_i v_i in A` for all integers `0 <= t_i <= L_i`, and `<v_i, u> != 0`
   for every `i`, then `nu_A(u) >= sum_i L_i`.
7. *(Transversal order.)* If `u_0 in C^d`, `c in V(A,m)` and `v in C^d`, then `Phi_c(u_0 + eta v)` vanishes to order
   `>= m - nu_A(u_0)` at `eta = 0`.

*Proof.* (1) `J` is an ideal. The equality for `TC` is compute.md R2 (`u in TC(A,m)` iff `l_u^m in J_m`). If
`V(A,m) = 0` then `J_m = S_m`. (2) For `i >= 1`, `l_u^i notin J_i` iff `Phi_c^{(i)}(u) := sum_k c_k <k',u>^i != 0` for
some `c in V(A,i)` (R2 applied with `m = i`); this is a finite union of complements of hypersurfaces.
(3) `l_{u+u'}^i = sum_j binom(i,j) l_u^j l_{u'}^(i-j)`; if `i > nu(u) + nu(u')` every term has `j > nu(u)` or
`i - j > nu(u')`, so lies in `J`. (4) Write `l_u = sum_j t_j y_j` with `y_j = <k, v_j>`. For `i > rho(pi(A))` the
leading-form ideal of `I(pi(A))` (in `Q[y]`) contains every form of degree `i`, so `(sum t_j y_j)^i = g(y)` on
`pi(A)` with `deg g < i`; then `l_u^i = g(pi(k))` on `A`, of degree `< i` in `k`. For `j = 1`, `pi(A)` is a set of at
most `h_A(v) + 1` integers, `rho <= h_A(v)`; for a family, combine with (3). (5) Fix a degree-compatible monomial
order on `S`; it restricts to degree-compatible orders on the two sets of variables. Reduced Groebner bases
`G_1, G_2` of `I(A_1), I(A_2)` (extended to `S`) have pairwise coprime leading monomials across the two sets, so by
Buchberger's criterion `G_1 u G_2` is a Groebner basis of `I(A_1) S + I(A_2) S = I(A)`, and since the order is
degree-compatible the top forms of a Groebner basis generate `gr I`: `J = J_1 S + J_2 S`. So
`S/J = S_1/J_1 (x) S_2/J_2` as bigraded algebras, and `l_u^i = sum_j binom(i,j) l_{u_1}^j (x) l_{u_2}^(i-j)` has its
`(j, i-j)` components in different bigraded pieces; in characteristic 0 it is nonzero iff some `j` has
`j <= nu_1(u_1)` and `i - j <= nu_2(u_2)`, i.e. iff `i <= nu_1 + nu_2`. (6) With `N = sum L_i`, compute.md
Prop 2.1 gives `c in V(A, N)` with `sum_k c_k e^(<k,z>) = e^(<k_0,z>) prod_i (e^(<v_i,z>) - 1)^(L_i)`. The
lowest-degree part of the right side is `prod_i <v_i,z>^(L_i)`, of degree `N`, and the degree-`N` part of the left
side is `Phi_c^{(N)}(z)/N!`. So `Phi_c^{(N)}(u) = N! prod <v_i,u>^(L_i) != 0`, and `l_u^N notin J`.
(7) The coefficient of `eta^i` is `binom(m,i) sum_k c_k <k,u_0>^(m-i) <k,v>^i`. If `m - i > nu(u_0)` then
`l_{u_0}^(m-i) = g` on `A` with `deg g < m - i`, so `l_{u_0}^(m-i) l_v^i` agrees on `A` with a polynomial of degree
`< m`, which `c` annihilates. QED

**Consequences for the hypothesis.** By (1), ZE-TC is a statement about the sublevel set `{nu_A < m}`. By (4),
spans of cheap families are tangent-cheap, but (4) also shows that the relevant cost of a rational subspace `M` is
`nu_A(M) <= rho(pi_M(A))`, which can be much smaller than the cheapest width-sum spanning `M` (§2). By (7), the
transversal-order statement of README §5.3 holds with `d_F` replaced by `nu_A(M)` for **every** rational
`M` in `TC(A,m)`. This suggests three nested forms of the qualitative hypothesis, for real points:

```text
ZE-TC^cheap(A,m): TC(A,m) n R^d  is contained in the union of the spans L_F of cheap families      (README §8.2)
ZE-TC^rat(A,m):   TC(A,m) n R^d  is contained in the union of the tangent-cheap rational subspaces
ZE-TC^deg(A,m):   every point of TC(A,m) n R^d is degenerate (lies on a rational hyperplane)
```

`cheap => rat => deg`. `ZE-TC^rat` is what the near-strata step of README §8.2 actually uses (transversal order
`m - nu_A(M)` by Lemma 1.1(7)); `ZE-TC^deg` is the weakest statement compatible with the tangent-cone reduction.

## 2. One cut: exact tangent degree; ZE-TC^rat holds, ZE-TC^cheap fails

Fix `S` with `|S| = s` and `L >= 0`, and write `R = [d] \ S`. Then
`A = A_S(L) x {-1,0,1}^R` with `A_S(L) = {k in {-1,0,1}^S : |k(S)| <= L}`. By Lemma 1.1(5) and the cube case
(on the uncut cube `nu = 2|supp|`: the cube is a product of segments `{-1,0,1}`, and on a segment `nu(u_x) = 2` if
`u_x != 0` since `k_x^2` is not affine on three points),

```text
nu_A(u) = nu_S(u_S) + 2 |supp u_R| ,      nu_S := nu_(A_S(L)).                                         (2.1)
```

For `u in C^S` let `sigma(u) = |supp u|` and `t(u) = s - max_a #{x in S : u_x = a}` (the number of coordinates
outside a largest value class).

**Theorem 2.1 (tangent degree of one cut).** Let `0 <= L <= s - 1` and `(s, L) != (1, 0)`. For every `u in C^S`,

```text
nu_S(u) = min( 2 sigma(u),  2L + 2 t(u),  s + L ).
```

In particular `rho(A_S(L)) = zeta(A_S(L)) = s + L` (if `L >= s` the slab does not cut and `nu_S = 2 sigma`).

*Proof. Upper bounds.* `u` lies in the span of the `sigma(u)` coordinate functionals of its support, each of width
`<= 2` on `A_S`; and `u = a 1_S + (u - a 1_S)` with `a` the value of a largest class, where `k(S)` has width `<= 2L`
and `u - a 1_S` is supported on `t(u)` coordinates. Lemma 1.1(3)-(4) give `nu <= 2 sigma` and `nu <= 2L + 2t`.
For `nu <= rho(A_S) <= s + L` (Lemma 1.1(1)) we interpolate along the levels of `sigma = k(S)`.
* *Slices.* For `|j| <= s` let `Lambda_j = {k in {-1,0,1}^S : k(S) = j}`. Every function on `Lambda_j` is a polynomial
  of degree `<= s - |j|`. Indeed, for `j >= 0` and `kappa in Lambda_j` with zero set `Z` and `N = {x : kappa_x = -1}`,
  the polynomial `prod_{x in Z} (1 - k_x^2) prod_{x in N} (k_x^2 - k_x)/2` has degree `|Z| + 2|N| = s - j` (because
  `#{+1} - |N| = j`), equals `1` at `kappa`, and vanishes at every other point of `Lambda_j`: it forces `k = 0` on `Z`
  and `k = -1` on `N`, and then the remaining `s - |Z| - |N| = #{+1}` coordinates sum to `j + |N| = #{+1}`, so they are
  all `+1`. For `j < 0` apply `k -> -k`.
* *Newton gluing.* Order the levels `j_0 = 0, j_1 = 1, j_2 = -1, ..., j_{2L-1} = L, j_{2L} = -L` and write a function
  `f` on `A_S` as `sum_i prod_{i' < i} (sigma - j_{i'}) q_i`, where `q_i` is chosen on `Lambda_{j_i}` to correct `f` there
  (the product is a nonzero constant on `Lambda_{j_i}` and vanishes on the earlier levels). The `i`-th term has degree
  `<= i + s - |j_i|`, and `max_i (i + s - |j_i|) = s + L` (at `i = 2L`). So `rho(A_S) <= s + L`.

*Lower bound* (Lemma 1.1(6)). We use four gadgets on disjoint coordinate sets, placed at the origin of the other
coordinates; each lies in `{-1,0,1}^S` and moves `k(S)` within `[-1,1]` (cost 1) or not at all (cost 0):

| gadget | coordinates | edges (`[0,1]` unless stated) | degree | cost | needs (`<v,u> != 0`) |
| --- | --- | --- | --- | --- | --- |
| single | `x` | `e_x` with length 2, from `-e_x` | 2 | 1 | `u_x != 0` |
| pair | `x, y` | `e_x - e_y` with length 2, from `e_y - e_x` | 2 | 0 | `u_x != u_y` |
| triangle | `x, y, z` | `e_x - e_y, e_y - e_z, e_x - e_z`, from `e_z - e_x` | 3 | 0 | `u_x, u_y, u_z` distinct |
| hexpair | `x, y` | `e_x, e_y, e_x - e_y`, from `-e_x` | 3 | 1 | `u_x, u_y != 0`, `u_x != u_y` |

(The triangle's 8 lattice points `k_0 + sum t_i v_i`, `t in {0,1}^3`, cover the 7 points of `{k in {-1,0,1}^3 : k_x + k_y + k_z = 0}`; the hexpair's are the
7 points of `{|k_x|, |k_y|, |k_x + k_y| <= 1}`.) A disjoint union of gadgets with total cost `<= L` is a lattice
zonotope in `A_S` whose degree is the sum of the degrees. Let `tau = min(2 sigma, 2L + 2t, s + L)`, `z = s - sigma`
the number of zeros, `M = s - t` the size of a largest class, `a` its value.
* *Case `tau = 2 sigma`.* Then `sigma <= z + L` (from `2 sigma <= s + L`). Pair `min(sigma, z)` nonzero coordinates
  with zero coordinates and use the remaining `sigma - min(sigma,z) <= L` nonzero coordinates as singles: degree
  `2 sigma`.
* *Case `tau = 2L + 2t < 2 sigma`.* Then `a != 0` (if `a = 0`, `t = sigma`), and `2t <= s - L` gives `M >= t + L`.
  Pair each of the `t` coordinates outside the class with a coordinate of the class, and use `L` further
  coordinates of the class (value `a != 0`) as singles: degree `2t + 2L`.
* *Case `tau = s + L < min(2 sigma, 2L + 2t)`.* Then `z < (s-L)/2` and every class has size `< (s+L)/2`. Put
  `h = (s - L)/2`. If `s - L = 1`: all coordinates are nonzero and there are two distinct values (`t >= 1`); take one
  hexpair on two coordinates with distinct values and `L - 1` singles (here `L >= 1`): degree `3 + 2(L-1) = s + L`.
  If `s - L >= 2`: choose `L` singles among the nonzero coordinates, first removing the excess `|c| - h` from every
  nonzero class `c` of size `> h`. This is possible because the total excess `sum_c (|c| - h)^+` is `<= L`: if one
  class exceeds, its excess is `< (s+L)/2 - h = L`; if `j >= 2` classes exceed, the excess is
  `<= s - j h <= s - 2h = L`. Then complete the choice arbitrarily among the remaining nonzero coordinates (there are
  `sigma > (s+L)/2 >= L` nonzero coordinates in all). The remaining `s - L` coordinates have all classes of size `<= h` (the zero class
  already has `z < h`). If `s - L` is even, a perfect matching with distinct values exists (complete multipartite
  graph, largest part `<= half`), giving `(s - L)/2` pairs; if `s - L >= 3` is odd, every class has size
  `<= (s - L - 1)/2`, so there are `>= 3` classes; take a triangle with one coordinate from each of three largest
  classes and match the rest (largest part now `<= (s-L-3)/2`). Degree `2L + (s - L) = s + L` in both cases.
In all cases `nu_S(u) >= tau`. With `u` generic this also gives `zeta(A_S) >= s + L >= rho(A_S)`. QED

*Exact check (evidence for the statement, not used in the proof).* `zetc_checks/onecut_test.py`: the formula agrees
with the mod-`p` value of `nu_S(u)` for all `u in {-1,0,1,2}^s` (`s <= 4`, all `L < s`) and 300 random integer `u`
for `s = 5, 6` (all `L < s`): 8157 cases, 0 disagreements.

**Corollary 2.2 (one cut).** Let `A` be a `T = 1` slab-box with one cutting column `S` (`1 <= L < s`, `d >= s`).
1. For every `m`, `TC(A,m)` (over `C`) is the finite union of the rational subspaces
   `{supp u in X}` (`2|X| < m`), `{u - a 1_S supported on X, a in C}` (`2L + 2|X| < m`) and
   `C^S x {u_R supported on Y}` (`s + L + 2|Y| < m`). In particular **ZE-TC^rat(A,m) holds at every `m`.**
2. **ZE-TC^cheap(A,m) fails for every `m` with `s + L < m <= rho(A) = 2d - s + L`**, if `L <= s - 2` and `d > s`.
   Witness: `u = (u_S, u_R)` with `u_S` having `s` distinct nonzero coordinates and `u_R` supported on `q` coordinates,
   `q` the largest integer with `s + L + 2q < m` (so `q <= d - s - 1`). Then `nu_A(u) = s + L + 2q < m`, but every
   family whose span contains `u` must span `C^S x C^(supp u_R)` (if `u_S` is chosen outside every proper rational
   subspace of `C^S`) and so has at least `s + q` members, each of width `>= 2` on `A` (`A` contains `+-e_x` for
   every `x`); its cost is `>= 2s + 2q >= m` because `m <= s + L + 2q + 2 <= 2s + 2q`.
3. The failure occurs inside LP-supercritical models: e.g. `s = 2`, `L = 1` (the hexagon), `d = 6`, `m = 8`:
   `h_A(beta) = 4 = m/2`, `rho(A) = 11`, `V(A,8) != 0`, and `u = (1, sqrt 2, sqrt 3, sqrt 5, 0, 0)` has
   `nu_A(u) = 3 + 4 = 7 < 8` while the cheapest family containing it costs `8`.

*Proof.* (1) By (2.1) and Theorem 2.1, `nu_A(u) < m` iff one of the three expressions plus `2|supp u_R|` is `< m`;
each condition says that `u` lies in one of the listed subspaces (with `X` = support of `u` or of `u - a 1_S` together
with `supp u_R`). Every listed subspace has `nu_A < m` by the upper bounds, so it is tangent-cheap. (2) If
`m > s + L` then `q >= 0` exists and, since `m <= rho(A) = s + L + 2(d - s)`, `q <= d - s - 1`; minimality of `q`
gives `m <= s + L + 2q + 2 <= 2s + 2q` as `L <= s - 2`. A family `{v_j}` whose real span contains `u` contains `u`
in a rational subspace, hence contains the smallest rational subspace containing `u`, which is
`R^S x R^(supp u_R)` for `u_S` outside proper rational subspaces of `R^S`. (3) is (2) with `q = 2`. QED

**Reading.** The single-cut tangent cone is a union of rational subspaces, as ZE-TC needs, but its largest pieces
(`R^S x ...`) are not spans of cheap families: the relevant cost of `R^S` is `nu(R^S) = s + L`, attained by the
non-product interpolation of Theorem 2.1, while every spanning family has width-cost `>= 2s`. These pieces are
degenerate (`u_x = 0` off `S u Y`), and the near-strata step of README §8.2 goes through for them with transversal
order `m - nu_A(M)` (Lemma 1.1(7)) and separating pairs `(x, n)`, `x notin S u Y`. **The formulation ZE-TC^cheap of
README §8.2 is therefore false already for one cut in LP-supercritical models, and should be replaced by
ZE-TC^rat**, with "cost" read as the tangent degree `nu_A(M)` of a rational subspace.

*Scope of the counterexample.* It needs a column with `|S| + L < m`. In the regime at the recommended order
`m = ceil(4w)` every balanced column has `|S| ~ n/2 >> m`, and then Corollary 2.2(1) shows that the one-cut
tangent cone *is* the union of cheap spans (the third family is empty). So the failure of ZE-TC^cheap is a failure
of the formulation in general `T = 1` models (small columns, or orders `m > |S| + L`), not yet an obstruction in the
regime; but there is no reason to expect width-cost to equal tangent degree for the many-cut regime boxes either, and
the near-strata step only needs ZE-TC^rat.

*Exact check.* `zetc_checks/cor22_check.py`: for `S = {0,1}`, `L = 1`, `d = 6`, exact rational rank profiles give
`rho(A) = 11` and `nu_A(u) = 7` for `u = (1,2,3,5,0,0)`, `(2,7,-3,11,0,0)`, `(1,2,1,1,0,0)`.

## 3. Few cuts: zonotopal models, products, and why an induction on cuts cannot reach the regime

### 3.1 Zonotopes settle the tangent cone

**Proposition 3.1 (zonotope certificate at the tangent-cone level).** If `A` contains a lattice zonotope
`k_0 + sum_{i <= s} [0, L_i] v_i` with `sum_i L_i >= m`, then `TC(A,m)` (over `C`) is contained in the union of the
`s` rational hyperplanes `<v_i, u> = 0`. In particular, if `zeta(A) >= m` then ZE-TC^deg(A,m) holds, and no
nondegenerate complex point lies in `TC(A,m)`.

*Proof.* Choose `0 <= alpha_i <= L_i` with `sum alpha_i = m`; by compute.md Prop 2.1 and the computation in the proof
of Lemma 1.1(6), the corresponding `c in V(A,m)` has `Phi_c(u) = m! prod_i <v_i,u>^(alpha_i)`. QED

So **ZE-TC has content only for zonotope-free pairs, `zeta(A) < m <= rho(A)`** (compute.md's "faithful" models),
exactly as for ZE itself; and by compute.md Prop 2.2 (with the referee's restriction to supports in the box) the
`T = 1` slab-box of an actual cluster in the regime is zonotope-free. In a *zonotopal* model (`zeta(A) = rho(A)`)
the tangent cone is settled at every admissible order by a single product certificate.

### 3.2 One and two cuts are zonotopal

* **One cut** (Theorem 2.1): `zeta = rho = 2d - s + L` (block value `s + L` plus `2` per free coordinate). So
  ZE-TC^deg(A,m) holds for every `m`, by one zonotope; ZE-TC^rat holds by Corollary 2.2.
* **Two disjoint cuts** (`S_1 n S_2 = {}`): `A` is a product of two one-cut blocks and a cube, so by Lemma 1.1(5)
  `nu_A` is the sum of the one-cut formulas and `TC(A,m)` is a finite union of products of the rational subspaces
  of Corollary 2.2: ZE-TC^rat holds at every `m` (Lemma 3.2 below), and `zeta = rho`.
* **Two crossing or nested cuts (computer-assisted, exact in the direction used).** Write `P = S_1 \ S_2`,
  `Q = S_1 n S_2`, `R = S_2 \ S_1` with sizes `a, b, c`. For every block with `a + b + c <= 7`, `b >= 1`, and
  `1 <= L_1 < a + b`, `1 <= L_2 < b + c` (`zetc_checks/twocut_zono.py`, all cases up to the symmetry `P <-> R`),
  an explicitly constructed and point-by-point verified lattice zonotope in `A` has total length equal to the mod-`p`
  top degree `rho_p(A)`. Since `zeta(A) <= rho_Q(A) <= rho_p(A)` (a mod-`p` rank is at most the rational rank, so
  the mod-`p` Hilbert function reaches `|A|` no earlier than the rational one), this proves `zeta(A) = rho(A)` for all
  these blocks, hence (free coordinates add `2` to both) for every two-cut `T = 1` model with `|S_1 u S_2| <= 7` and `L_1, L_2 >= 1`
  (slabs that do not cut reduce to one cut; disjoint cuts are products).
  Result: **421 of 421 blocks certified** (`twocut_zono_out.txt`). By Proposition 3.1, ZE-TC^deg holds for every
  such model at every order.
* *Conjectured closed form* (fits all these blocks): `rho = min(2d - |S_1| + L_1, 2d - |S_2| + L_2, d + L_1 + L_2 - [b = 1])`
  with `d = a + b + c` the block size. The first two terms are the one-cut bounds (`A` lies in each one-cut box);
  I have not proved the third.
* **Three cuts** (`zetc_checks/threecut_zono.py`): all Venn types of three distinct columns (each of size `>= 2`) on
  block size `<= 6`, up to permuting the columns, with `L_j in {1,2}` (`L_j < |S_j|`): **953 of 953 blocks have a
  verified zonotope of length `rho_p`**, so they are zonotopal and ZE-TC^deg holds at every order.

* **ZE-TC^rat beyond one cut: an exact plane test** (`zetc_checks/plane_test.py`). For a rational plane
  `M = span(v_1, v_2)`, the points `u = v_1 + t v_2` with `nu_A(u) < nu_A(M)` are the common real roots of the
  univariate polynomials `Phi_c^{(N)}(v_1 + t v_2)`, `c in V(A, N)`, `N = nu_A(M)`. An irreducible factor of degree
  `>= 2` of their gcd with a real root would give a point whose rational hull is `M` but whose tangent degree is
  smaller: an exact counterexample to ZE-TC^rat. Result: 12 random planes (integer generators in `[-3,3]`) in each of five
  crossing/nested two-cut blocks of size 5-6 and in `M8` (72 planes): every gcd factors over `Q` into linear factors
  (degrees 0-8, no irreducible factor of degree `>= 2` at all), so **no violation of ZE-TC^rat on any tested plane**.
  This is exact for the planes tested and evidence for the statement. (For `M8` the gcd was trivial on all 12 planes,
  consistent with J4.)

### 3.3 Products

**Lemma 3.2 (products).** Let `A = A_1 x A_2` (disjoint coordinates). Then for every `m`
`TC(A, m) = U_{a + b = m - 1} TC(A_1, a + 1) x TC(A_2, b + 1)`. Hence if ZE-TC^rat (resp. ZE-TC^deg) holds for `A_1`
and `A_2` at **every** order, it holds for `A` at every order. Moreover `rho(A) = rho(A_1) + rho(A_2)`,
`zeta(A) >= zeta(A_1) + zeta(A_2)`, and `zeta(A) <= zeta_LP(A_1) + zeta_LP(A_2)` for the LP bound of §5.1.

*Proof.* `nu_A(u) = nu_1(u_1) + nu_2(u_2)` (Lemma 1.1(5)), so `nu_A(u) < m` iff `nu_1(u_1) <= a` and
`nu_2(u_2) <= b` for some `a + b = m - 1`. A product of rational subspaces is rational, and a product with a degenerate
factor point is degenerate. `rho` is the case `u` generic. Products of zonotopes are zonotopes. For the LP bound,
a primitive difference direction `(v_1, v_2)` of `A` has `v_1` or `v_2` a nonzero integer multiple of a primitive
difference direction of the factor, so dual certificates add. QED

### 3.4 Sections: ZE-TC^deg is inherited from any section

**Proposition 3.3 (section monotonicity).** Let `V in Z^(d x j)` have rank `j`, `k_0 in Z^d`, and
`A' = {t in Z^j : k_0 + V t in A}`. Then for every `u in C^d`, `nu_A(u) >= nu_(A')(V^T u)`. Consequently, if
ZE-TC^deg(A', m) holds, then ZE-TC^deg(A, m) holds.

*Proof.* If `l_u^i = g` on `A` with `deg g < i`, then on `A'`,
`(<k_0,u> + <t, V^T u>)^i = g(k_0 + V t)`, so `<t, V^T u>^i = g(k_0 + Vt) - sum_{j<i} binom(i,j) <k_0,u>^(i-j) <t,V^T u>^j`
has degree `< i` in `t`. Hence `i > nu_A(u)` implies `i > nu_(A')(V^T u)`. For the consequence: if `u in TC(A,m)`
is real then `nu_(A')(V^T u) < m`, so `V^T u in TC(A',m)`, so (ZE-TC^deg for `A'`) `<w, V^T u> = 0` for some
`0 != w in Z^j`; then `Vw in Z^d \ 0` (rank `j`) and `<Vw, u> = 0`: `u` is degenerate. QED

*Remarks.* (1) When the edges are linearly independent, Proposition 3.1 is the case where `A'` is a box
`prod [0, L_i]` (`V` = the edge vectors), on which
ZE-TC^deg holds at every `m <= sum L_i` (the box is a product of segments; Lemma 3.2). (2) ZE-TC^deg holds for a product
at every order `m <= rho` as soon as each factor satisfies it at its own top degree `rho_i`: if `u` is in
`TC(prod A_i, m)` then `sum nu_i(u_i) < m <= sum rho_i`, so some `nu_i(u_i) < rho_i` and `u_i` is degenerate. So
**any `A` containing, as a section, a product of pieces that are "TC-good at their top degree" (such as segments,
one-cut blocks, `M8` by README §5.4 J4) with total top degree `>= m` satisfies ZE-TC^deg(A, m).** This is the
tangent-cone form of compute.md Prop 2.3, and it is the only inheritance mechanism I found. ZE-TC^rat does not pass
to `A` this way (the pulled-back subspace `(V^T)^(-1)(M')` need not be tangent-cheap for `A`).
(3) *Limits in the regime (a count, not a theorem).* On a set `Y` of points with `|Y| <= w` every slab of the regime
box is slack, so the section `{k supported on Y}` is a full ternary cube: such sections only reproduce the
coordinate certificates `nu_A(u) >= 2|Y|`, which stop near `2w ~ m/2`. A product of `k` copies of an 8-point
Hadamard piece (`M8`, or its variant with the extra base box `|k_(y_0)| <= 1` for the blocks not containing the
actual base point; ZE-TC is known only for `M8` itself) placed on disjoint 8-point blocks fits in the regime box only
while `sum_blocks max |<k, beta_p>| <= w` for every column. Since `|<k, beta|_Y>| = |sum_(a != 0) c_a hat k(a)|
<= 4 ||c||_1 <= 4 sqrt 8 < 11.4` on such a block (`c` the Walsh coefficients of the sign vector `beta|_Y`,
`||c||_2 = 1`), this allows about `w/11.4` blocks, of total top degree `8k < 0.71 w`, far below `m >= 4w`. Certifying
the regime through sections needs TC-good pieces whose top degree is at least about four times their slab usage
in every column; zonotopes cannot do it (compute.md Prop 2.2) and `M8` is a factor `~5.7` short.

### 3.5 Reading: the induction on cuts

Round-10 Prop R1 (no cuts) and §3.2 (one or two cuts) are all *zonotopal* statements: the tangent cone is settled
by one product certificate (Proposition 3.1), i.e. by a character. A zonotope-free model needs many cuts: every
direction `v` of `A - A` must be "paid for" by the slabs (compute.md Prop 2.2 is this LP), and with `r` cuts the
`d - r`-dimensional lattice of directions orthogonal to all columns already produces long zonotopes when `r` is small.
In the computations above the zonotope-free `T = 1` models found all have at least 6 distinct columns on `d <= 7`
coordinates (`Z6` of §5.2, compute.md's `M8`), while every three-cut block on `<= 6` coordinates is zonotopal. So **an induction on the number of cuts that propagates
zonotopality cannot reach the regime**: the step from zonotopal to zonotope-free is exactly where ZE-TC acquires
content, and at that step a single product certificate stops existing (Prop 2.2 of compute.md). Any induction must
carry a non-product invariant (for instance the full system `V(A,m)` restricted to sections), which is the
overdetermination phenomenon of round 11.

## 4. The Sylvester-Hadamard family

Notation: `n = 2^k`, points `x in F_2^k`, base point `0`, columns `chi_a(x) = (-1)^(a.x)` (`a != 0`), patterns
`S_a = {x : a.x = 1}` (`|S_a| = n/2`), and `Syl(2^k, L) = {k in Z^n_0 : |k_x| <= 1 (x != 0), |hat k(a)| <= 2L (a != 0)}`
with `hat k(a) = <k, chi_a> = -2 k(S_a)`. So `width(k(S_a)) <= 2L`, `h_A(beta_a) <= 4L`, and LP-supercritical means
`m >= 8L` (when the slab is attained). `M8 = Syl(8,2)`, `M8' = Syl(8,3)` (compute.md §3.1).

### 4.1 The tower: sections along subgroups

**Lemma 4.1 (subgroup sections).** Let `U` be a subgroup of `F_2^k` of order `2^j`. The section
`{k in Syl(2^k, L) : supp k in U}` is (identifying `U = F_2^j`) exactly `Syl(2^j, L)`.

*Proof.* For `k` supported on `U`, `hat k(a)` depends only on the restriction of `chi_a` to `U`, which is a character of
`U`; it is trivial iff `a in U^perp`, and then `hat k(a) = sum_U k = 0`. Every nontrivial character of `U` arises
(from `2^(k-j)` values of `a`), and the box conditions on `U \ {0}` are those of `Syl(2^j, L)`. QED

**Theorem 4.2 (monotonicity in `k`).** If ZE-TC^deg(Syl(2^j, L), m) holds, then ZE-TC^deg(Syl(2^k, L), m) holds for every
`k >= j`. In particular (from README §5.4 J4, which rests on the refereed compute.md Thm 5.2):

```text
ZE-TC^deg(Syl(2^k, 2), 8)   and   ZE-TC^deg(Syl(2^k, 3), 11)   hold for every k >= 3.
```

In both cases `rho(Syl(2^k, L)) >= rho(Syl(8, L)) = 8` resp. `11` (Prop 3.3 with `u` generic), so `V(A, m) != 0`. For
`L = 2` the pairs are zonotope-free for every `k >= 3` (compute.md Prop 3.1: `zeta <= ((n-1)/n) 8 < 8 = m`). For `L = 3`
Prop 3.1 only gives `zeta <= ((n-1)/n) 12`, which is `< 11` only at `k = 3` (`zeta(M8') <= 10`); for `k >= 4` I have not
decided whether `zeta < 11` (if not, the statement is trivially true by Prop 3.1).

*Proof.* Proposition 3.3 with the section of Lemma 4.1 (`V` = inclusion of `Z^U_0`). QED

*What this does and does not give.* These are the first zonotope-free statements of ZE-TC for an infinite family
with growing `n`; the number of tangent forms is `H_A(8) = 42` at `k = 3` (and `4345`, mod `p`, for the Paley analogue
at `n = 12`, §5.5), so these systems are overdetermined (I did not compute `H_A(m)` for `k >= 4`). But the orders are LP-subcritical
(`m = 8 < 2 h = 16`, resp. `11 < 24`): at these orders the regime's LP comparison is not available, and the
statements say nothing about orders `m >= 8L`.

### 4.2 What the symmetry can and cannot do

**Lemma 4.3 (nondegenerate points have trivial stabilisers).** Let `G` be a group of permutations of the points
fixing the base point and preserving `A` (for `Syl`, `G = GL(k,2)`; for the variant with the box also at `0`,
`G = AGL(k,2)` acting on `Z^n_0`). Then `TC(A,m)` and the set of tangent-cheap subspaces are `G`-invariant, and a real
`u` fixed by some `1 != g in G` is degenerate.

*Proof.* `G` permutes the coordinates of `k` and preserves `A`, hence `I(A)`, `J`, `V(A,m)`, `nu_A` and `TC`. If
`g u = u` with `g x != x` for some point `x`, then `u_x = u_(g x)` (with `u_0 = 0`), a collision. QED

So reductions to fixed-point subspaces or orbit sums (`u` in `Fix(H)`, `H <= G`) only probe **degenerate** points:
they decide which symmetric rational subspaces are tangent-cheap, never whether a nondegenerate point lies on `TC`.
What symmetry does provide is (i) a block decomposition of the *equations* (the isotypic components of
`span{Phi_c} = (S/J)_m^*`), which makes the Hilbert-function computations of §5 feasible (characters of the Singer
cycle or of a cyclic subgroup; `zetc_checks/cyclic_hilbert.py`), and (ii) the tower of Lemma 4.1. Invariant forms
(the trivial isotypic part) cut out a `G`-invariant superset of `TC`; their number is small (one per invariant of
degree `m`), so they cannot by themselves exclude nondegenerate zeros.

### 4.3 The ball effect: small Hadamard boxes are thin

**Lemma 4.4.** For a Hadamard profile of order `n` (rows `h_1 = 1, h_2, ..., h_n` orthogonal) and slabs
`|<k, h_p>| <= w` (`p >= 2`), every `k in A` satisfies `|k|_2^2 <= (n-1) w^2 / n < w^2`.

*Proof.* `sum_{p=1}^n <k,h_p>^2 = n |k|_2^2` and `<k, h_1> = 0`. QED

With `w = 2L` this says that `A` lives on vectors with fewer than `4L^2` nonzero entries. In the regime
(`w^2 = 8(n-1) ln 4r > n`) the ball is slack, which is why `|A| >= 3^(n-1)/2` there; but for small Hadamard models
it is binding (`L = 1`: at most 3 nonzero entries, `|Syl(16,1)| = 241`), and the cuts are so heavy that
`rho(A)` stays far below `8L` (§5.3). Consequently the Sylvester models in which the regime's LP condition
`m >= 8L` can be met are large (§5.3), and Theorem 4.2 does not reach them.

## 5. Computation: where are the LP-supercritical zonotope-free `T = 1` models?

### 5.1 Tools (all in `zetc_checks/`)

* `tclib.py`: point sets, mod-`p` Hilbert functions (`p = 2^61 - 1`), and `nu_A(u)` by one rank-profile computation
  (columns: monomials by degree with `l_u^i` inserted before degree `i`).
* `certify.py`: exact Hilbert functions (`fmpz_mat` ranks) and **exact LP certificates** for `zeta`: rational
  `y_f >= 0` with `sum_f y_f |<v,f>| >= 1` for every primitive direction `v` of `A - A` (checked with `Fraction`s);
  then `zeta(A) <= sum_f y_f width_A(f)` (proof in the module docstring: widths of a zonotope add).
* `zono.py`: lower bounds for `zeta` by MILP over edge multiplicities, **each reported zonotope placed and verified
  point by point**.
* `cyclic_hilbert.py`, `singer16.py`: Hilbert functions / spanning tests through the isotypic decomposition under a
  cyclic group of coordinate shifts (Paley translations `Z_11`; Singer cycle `Z_15` times `k -> -k` for `n = 16`).
  Validated against the direct computation on `M8` (identical Hilbert function `[1,7,28,77,161,259,287,175,42]`).
* **Direction of mod-`p` results.** A mod-`p` rank is `<=` the rational rank, so the mod-`p` Hilbert function reaches
  `|A|` no earlier than the rational one: `rho_Q <= rho_p`. Hence every *upper* bound on `rho` (and every "monomials
  of degree `<= j` span" statement) proved mod `p` is rigorous; lower bounds on `rho` from mod-`p` data are evidence
  unless recomputed exactly.

### 5.2 Small models are zonotopal

| family | range | result | status |
| --- | --- | --- | --- |
| one cut | all | `zeta = rho = 2d - s + L` | Theorem 2.1 |
| two cuts (blocks `P,Q,R`, `b >= 1`, `L_j >= 1`) | `a+b+c <= 7`, all `L` | 421/421 blocks: verified zonotope of length `rho_p` | exact (§3.2) |
| three cuts (all Venn types up to `S_3`) | block size `<= 6`, `L_j in {1,2}` | 953/953 blocks: verified zonotope of length `rho_p` | exact |
| random profiles, `2 <= r <= 2d` columns, `L in {1,2}` | `d = 6`: 400 models (273 twin-free); `d = 7`: 200 models (122 twin-free) | `zeta_LP < rho_p` for 1 model at `d = 6`; none is LP-supercritical with `zeta_LP < m` | LP / mod-`p` |
| twin-free, all `L = 1` (so `m >= 8` is LP-supercritical) | `d = 7` (600 trials), `d = 8` (300 trials), patterns of size 2-4 | every model has `zeta_LP = rho_p` (no certified zonotope-free model); `rho_p <= 9` | LP / mod-`p` |
| graph models `|k_x + k_y| <= 1` | `C_5..C_9`, `K_33`, `K_44`, `Q_3` | `zeta_LP = rho_p` in all; bipartite graphs have the free alternating direction | LP / mod-`p` |

**The smallest zonotope-free model found, `Z6`** (`z6_model.py`, exact): `d = 6` (`n = 7`), columns
`{1,2,3,4,5}, {0,...,5} (L = 2 and L = 1), {0,3,4}, {1,4,5}, {1,2,3,4}, {0,2,4,5}` with
`L = (2,2,1,2,2,1,1)`; twin-free; `|A| = 247`, exact Hilbert function `[1,6,21,47,72,67,30,3]`, so `rho = 7` and
`H_A(7) = 3`; exact LP certificate `zeta <= 13/2`, verified zonotope of length 6, so `zeta = 6 < 7 = rho`. It is smaller
than `M8` (`n = 8`), but underdetermined (`3` tangent forms in `6` variables) and LP-subcritical (`2 h = 16 > 7`).

So below `n ~ 8` zonotope-freeness is rare, needs many cuts relative to `d`, and occurs only at `m = rho` with very few
tangent forms. None of the models found is LP-supercritical.

### 5.3 Hadamard models

Hadamard profiles are automatically zonotope-free at LP-supercritical orders (compute.md Prop 3.1:
`zeta <= ((n-1)/n) 4L < 4L <= m/2`), so for them the only question is whether `rho(A) >= 8L`.

| `n` | profile | `L` | `|A|` | `rho` | needs `rho >=` | verdict |
| --- | --- | --- | --- | --- | --- | --- |
| 8 | `H_8` (unique) | 1 / 2 / 3 | 57 / 1037 / 1877 | 3 / 8 / 11 (exact) | 8 / 16 / 24 | no (`L >= 4` does not cut) |
| 12 | Paley (unique `H_12`) | 1 | 133 | 3 (exact) | 8 | no |
| 12 | Paley | 2 | 20153 | `<= 8` (mod `p`, rigorous upper bound) | 16 | no |
| 12 | Paley | `>= 3` | | `<= 2d - n/2 + L = 16 + L` (one-cut bound) | `8L` | no |
| 16 | Sylvester | 1 | 241 | 3 (exact) | 8 | no |
| 16 | Sylvester | 2 | 213937 | `<= 15` (Singer-cycle spanning test, mod `p`: rigorous) | 16 | no |
| 16 | Sylvester | 3 | 3608377 | `<= 22` (fibration bound, Lemma 5.1) | 24 | no |
| 16 | Sylvester | `>= 4` | | `<= 22 + L < 8L` (one-cut bound) | `8L` | no |

**Lemma 5.1 (fibration bound).** Let `A` be finite and `sigma: Z^d -> Z^j` linear, with fibres
`F_beta = A n sigma^(-1)(beta)` (`beta in B = sigma(A)`), and suppose every function on `F_beta` is a polynomial of degree
`<= delta(beta)`. Order `B = (beta_1, beta_2, ...)` and let `e_i` be the least degree of a polynomial on `Q^j` vanishing at
`beta_1, ..., beta_(i-1)` and not at `beta_i`. Then `rho(A) <= max_i (e_i + delta(beta_i))`.

*Proof.* With such polynomials `N_i`, write `f = sum_i N_i(sigma(k)) q_i(k)`, choosing `q_i` (degree `<= delta(beta_i)`) on
`F_(beta_i)` inductively so that the sum equals `f` there: the terms `i' > i` vanish on `F_(beta_i)` and `N_i(beta_i) != 0`.
QED

Applied to `Syl(16, L)` with the three columns `a, b, a + b`: the 15 nonzero points split into three regions of size 4
(`(a.x, b.x) = (1,1), (1,0), (0,1)`) and three free points, `Syl(16,L) in Block x {-1,0,1}^3`, the block fibres over
the three region sums `beta` in `[-4,4]^3` (with the three slab conditions) into products of slices of degree
`<= 12 - |beta|_1` (Theorem 2.1 proof). With `B` ordered by `|beta|_1` and **exact** Newton degrees
(`zetc_checks/fibration_bound.py`), `rho(Block) <= 16` for `L = 3` and `<= 15` for `L = 2`, so
`rho(Syl(16,3)) <= 22 < 24` and `rho(Syl(16,2)) <= 21`.

(`H_8` and `H_12` are unique up to equivalence, and their automorphism groups are transitive on points, so the base
point does not matter; for `n = 16` only the Sylvester matrix was examined.)
**So no Hadamard `T = 1` model with `n <= 12` is LP-supercritical with `V(A,m) != 0`**, and at `n = 16` (Sylvester)
no `L` works either (`rho <= 15 < 16` at `L = 2`, `rho <= 22 < 24` at `L = 3`, one-cut bound for `L >= 4`): **the
smallest LP-supercritical Sylvester `T = 1` model has `n >= 32`** (`d >= 31`, beyond exact computation of tangent
cones). The reason is the ball effect (Lemma 4.4): at these sizes the slabs
cut the cube down to a set on which every function has degree `~ 4L`, while the LP needs `8L`. In the regime the
ball is slack and `rho` is of order `n` (README §5.3, J3).

### 5.4 Products satisfy every stated criterion and carry no information

`M8 x M8` (two Sylvester blocks on disjoint 7-point sets, `n = 15`, 14 columns, each column constant on the other
block) with `m = 16`:
* LP-supercritical: every column has `h_A = 8 = m/2`;
* zonotope-free: `zeta <= 7 + 7 = 14 < 16` (Lemma 3.2 with compute.md Prop 3.1 for each block);
* `rho = 16`, `dim V(A,16) = H_A(16) = 42^2 = 1764 >> n`: strongly overdetermined;
* twin-free (points of one block are separated from each other and from the base by that block's columns).

And ZE-TC^deg(`M8 x M8`, 16) **holds** (Lemma 3.2(2) with J4). More generally every product of zonotope-free pieces
that are TC-good at their top degree passes all the criteria of the task once enough factors are taken, and ZE-TC
then reduces to the factors. **So "zonotope-free + LP-supercritical + overdetermined" must be supplemented by
indecomposability** (or by a regime-like profile: balanced columns, pair separation close to `r/2`) for a model to test
anything new.

### 5.5 An overdetermined zonotope-free model where the sections are bad but the system is not (evidence)

Paley-12, `L = 2`, `m = 8`: zonotope-free (`zeta <= 7`), `rho_p = 8` with `H_p(8) = 4345` tangent forms in 11
variables (`rho_Q = 8` is proved: `<= 8` mod `p`, and `>= 8` because the exact section measures below push forward to
nonzero elements of `V(A,8)`; the value `4345` is mod `p`), LP-subcritical (`2h = 16`).
* Each of its 330 sections on 8 points containing the base (30 classes under translation) has `rho_p = 8`,
  `H_p(8) = 1` and a verified zonotope of length 7 (`paley_sections.py`). For three classes (one of each size
  `|A_Y| = 909, 901, 915`) the exact top form is `-8!` times an **irreducible, indefinite** octic in 7 variables
  (`section_forms.py`, exact factorisation over `Q`; `section_sign.py`, both signs at sampled points). So, as in compute.md
  Prop 6.2, each of these sections alone has nondegenerate real tangent-cone points: no section certifies ZE-TC (Prop 3.3 gives
  nothing).
* By Prop 3.3, `TC(A,8)` lies in the common zero set of all 330 section octics. A numerical search
  (`paley_numeric.py`, 77 least-squares runs from random starts on the unit sphere of `R^11`; **evidence**) found
  approximate common zeros only near degenerate strata: the 66 runs with residual `< 10^-8` end at two-term distance
  `min(|u_x|, |u_x - u_y|)` with median `1.2 * 10^-4` and maximum `1.6 * 10^-3` (61 of 66 below `5 * 10^-4`;
  coordinates recorded to 4 digits), against a median `5.7 * 10^-3` for random unit vectors. The end points approach
  either a sparse stratum (one large coordinate) or the diagonal (all non-base points colliding). The same pipeline
  applied to one section octic, or to three random sections, finds zeros that are no more degenerate than random
  points (`numeric_validate.py`); with ten random sections the zeros found are already degenerate.
* So, on this evidence, the overdetermined system kills the nondegenerate zeros that every individual section has,
  exactly the "ZE is an overdetermination statement" picture of round 11, now at `n = 12` with `4345` forms. It is
  LP-subcritical and remains evidence.

## 6. Open points

1. **ZE-TC^rat for two crossing cuts in general.** Proved only for disjoint cuts; for crossing/nested cuts only the
   zonotopal (deg) form on `<= 7` coordinates. The conjectured closed form
   `rho = min(2d - |S_1| + L_1, 2d - |S_2| + L_2, d + L_1 + L_2 - [b = 1])` (block size `d`) fits all 421 blocks; a proof
   needs an interpolation argument beyond the tensor Newton scheme (double slices do not have degree
   `d - |j| - |l|`; `zetc_checks/dslice.py` lists 469 violations).
2. **Indecomposable LP-supercritical zonotope-free models.** The smallest size is unknown. Hadamard profiles need
   `n >= 16`, and Sylvester profiles need `n >= 32` (§5.3); the four non-Sylvester Hadamard matrices of order 16 were not examined. Beyond that, deciding ZE-TC exactly (thousands of forms of degree `>= 16` in `>= 15`
   variables) is out of reach of the methods used here; the section method (Prop 3.3) is the only exact tool that
   scales, and it needs TC-good sections of top degree `>= m`.
3. **TC-good pieces with large degree-to-slab ratio.** Prop 3.3 would prove ZE-TC^deg in the regime if the regime box
   contained, as a section, a product of TC-good pieces of total top degree `>= m = ceil(4w)`. This needs pieces whose
   top degree is about four times their slab usage in every column; none is known (`M8`: `~0.7`).
4. **Paley-12 at `m = 8`, exactly.** Decide whether the 330 section octics (or the full 4345-form system) have a
   nondegenerate real common zero. The numerical evidence says no.
5. **Quantitative form.** With ZE-TC^rat the near-strata step uses transversal order `m - nu_A(M)` (Lemma 1.1(7)); the
   Lojasiewicz constants and the arithmetic separation from the tangent-cheap subspaces (which include non-cheap ones
   such as `R^S x R^Y`) remain to be done.

## 7. Files (`zetc_checks/`)

| script | content | status |
| --- | --- | --- |
| `tclib.py` | point sets, mod-`p` and exact `nu_A(u)`, mod-`p` Hilbert functions | library |
| `onecut_test.py` | Theorem 2.1 formula vs `nu` (8157 cases) | mod `p`, evidence for a proved statement |
| `cor22_check.py` | Corollary 2.2(3) example | exact |
| `twocut_rho.py`, `twocut_zono.py` | two-cut blocks: `rho_p` and verified zonotopes (421 blocks) | exact in the direction used |
| `threecut_zono.py` | three-cut blocks (953) | exact in the direction used |
| `gadget_ip.py`, `twocut_formula.py`, `dslice.py` | exploratory (gadget IP, atom bound, double slices) | evidence |
| `zono.py`, `certify.py` | verified zonotopes; exact Hilbert functions; exact LP certificates for `zeta` | library |
| `faithful_search.py`, `fs_*.py` | random / graph / design / `L = 1` surveys | mod `p` + LP, evidence |
| `z6_model.py` | the model `Z6` | exact |
| `small_hadamard.py` | `H_8`, Paley-12, `Syl(16)` at `L = 1` | exact |
| `cyclic_hilbert.py` | Paley-12 `L = 2` Hilbert function via `Z_11` | mod `p` (upper bound on `rho` rigorous) |
| `singer16.py` | `Syl(16, 2)` spanning test via Singer cycle (`singer16_L2_summary.txt`: 30/30 components full rank) | mod `p` (full rank = rigorous upper bound) |
| `fibration_bound.py` | Lemma 5.1 bound for `Syl(16, L)` | exact |
| `paley_sections.py`, `section_forms.py`, `section_sign.py` | Paley-12 sections: `rho`, exact factorisation, sign change | exact / evidence (sign) |
| `paley_numeric.py`, `numeric_validate.py` | numerical common zeros of the section system | evidence |
| `froberg.py`, `had_size.py`, `had_size2.py` | heuristics and counts | evidence / exact counts |
