# Round 24: twisted Pell-type search (two conjugate real-quadratic branches, B.C = 30, balanced, alpha = 15/4)

## 0. Verdict

**No growth improvement and no uniform bound.** `M <= (2/3 + eps) log R/loglog R` stands. Status: **partial**, with a
search-negative outcome. The class searched (round 23 `split.md` §6.2 item 3) is: rational curves over `Q` on the twisted
Losev–Manin 4-fold through `e`, with a pair of conjugate real-quadratic branches of contact `mu_0 = 4` each, of the
balanced minimal class (`E_A.C = 1` for all 30 rays, so `B.C = 30`, every cut of weight 2), hence with `alpha_C = 15/4`.

1. **Normal form (proved, §1).** Galois forces all 30 boundary places to be `Q(i)`-rational points of `P^1_t`. So the
   "Pell conic with Gaussian quadratic blocks" ansatz is the same as `P^1` with 15 Gaussian **linear** blocks `t - r_T`,
   `r_T in Q(i)`, and branch pair `+-sqrt m`. The exact system is `Im(gamma_j P_j) = kappa_j (t^2 - m)^4`, `j = 1..4`.
   Equivalently, after a real Möbius change, the system is (1.3): the 8-point sets `R_j` have `e_k` real for `k != 4`.
2. **Reformulations (proved, §4; symbolic identities checked exactly).** (a) A five-column two-sided Prouhet–Tarry–Escott
   condition. (b) A half-angle "star" form on the 15 vertices and edges of `K_5`, with five degree-5 polynomials `z_a`. (c) A
   Möbius form: five real Möbius maps together with one degree-4 map `gamma` whose real trace touches `RP^1` to order 4 at
   two points. In this form contact `(4,4)` is automatic, and the only condition left is that the graph of `gamma` pass
   through one coincidence point of each of the 10 pairs.
3. **Expected dimension 1 (heuristic count, calibrated, §3).** The count reproduces split Thm 2's 2-dimensional family,
   the record's 3-dimensional `k = 4` family and the 0-dimensional single-branch problem. So **if such curves exist, they
   come in 1-parameter families and sweep at most a surface, not a hypersurface.** The split analogue sweeps a hypersurface
   `W`. The literal transplant of Thm 2's structure (three coordinates pulled back from the conic's `x`-line) is overdetermined
   by 3. This is the first place where the anisotropy enters a **dimension count**, not only a parity.
4. **Search (numerical evidence, §5).** No nondegenerate curve was found. The solvers were LM in five parametrisations, all
   symmetric single-involution ansätze, complexified Newton, staged continuation, lifting from doubled `k = 4` curves, and
   pseudo-arclength homotopies. All are calibrated on problems that do have solutions: twisted `k = 4`, doubled `k = 4`, `k = 5`
   at lower contact, and the split Theorem-2 class itself. Every homotopy path that reaches the target system ends at a
   **collapsed** configuration, with four of the five points identical. This supports **Q24.1**: *the class is empty*. It is
   not a proof.
5. **Bearing on Question B (§6).** Nothing is refuted and nothing is proved. A twisted counterpart of the split hypersurface
   `W` is excluded by the count (heuristic) and absent in every computation. So the Pell mechanism that forces `Z ⊇ W` in the
   split form does not transfer to `T_5`. This is map value for the anisotropy. Curves of higher class
   (`B.C = 60, ...`) and sporadic points are untouched.


## 1. Setting and the curve class

Labels: **proved** (proof here), **exact** (exact computation in `search_checks/`), **numerical** (floating point,
stated precision), **heuristic** (dimension counts that are not theorems), **record** (proved earlier, cited).

**1.1 Objects (record).** `X = Lbar_5^tw` is the permutohedral toric 4-fold over `Q` whose torus is the anisotropic
norm-one torus `T_5 = (R^1_{Q(i)/Q} G_m)^5/diag`; Galois acts on cocharacters by the quadratic character of `Q(i)`.
`E_A` (`A` proper nonempty subset of `[5]`, 30 rays) are the boundary divisors, `B = sum_A E_A = -K`, `e` the identity.
Normalising the fifth coordinate, a point of `T_5` is `(u_1,...,u_4)`, `u_j in Q(i)`, `N(u_j) = u_j ubar_j = 1`. For a
curve `C` through `e` with a branch at a point `t_0` of its normalisation, `mu_0 = min_j ord_{t_0}(u_j - 1)`.

**1.2 The class searched (round 23 `split.md` §6.2 item 3).** Rational curves `phi: P^1_Q -> X` over `Q`, not in the
boundary, with
* (balanced, minimal) `E_A . C = 1` for all 30 rays, so every cut has weight `2` (the parity lemma forces even weight)
  and `B . C = 30`;
* (ray-generic) the boundary places are 30 distinct points, each with valuation vector a ray `1_A`
  (places with chain vectors `1_A + 1_B`, `A` inside `B`, are the limits where two places collide; they keep the class
  balanced but were not searched separately, and none occurred among the computed limits, §5.2);
* (two conjugate real-quadratic branches) `phi^{-1}(e)` contains a pair `t_0, t_0'` of conjugate points of degree 2
  with `t_0` real (equivalently, after `PGL_2(Q)`, `t_0 = sqrt m`, `m in Q_{>0}` not a square), and
  `mu_0 = 4` at each.

Such a curve has `alpha_C = B.C/(2 mu_0) = 30/8 = 15/4` (**proved** as in split Thm 2(iv): `h_{-K}(phi(t)) = 30 h(t) + O(1)`
because `-K` is base-point free, `dist(phi(t), e) ~ |t - t_0|^4`, and Liouville/Dirichlet at the quadratic `t_0` give
`|t - t_0| >> H(t)^{-2}` with equality infinitely often). It is the twisted analogue of split Theorem 2, with the
weights doubled as the parity lemma requires. By the two-branch covering bound (`split.md` 5.2) such curves cannot cover `T_5`;
a family of them can at best force a proper subvariety into every admissible `Z`.

**1.3 Lemma 1 (normal form; proved).** Let `phi` be in the class 1.2. Then there are `r_T in Q(i)` (`T` nonempty subset of `[4]`;
after a `PGL_2(Q)` change of parameter all are finite), pairwise distinct with `r_T != rbar_S` for all `T, S`, and
`gamma_j in Q(i)^*`, such that, with `P_j(t) = prod_{T contains j} (t - r_T)` (degree 8),

    u_j(t) = (gamma_j P_j(t)) / conj(gamma_j P_j)(t),   j = 1..4,                                   (1.1)

where `conj` conjugates coefficients. The place `t = r_T` has vector `1_T` and `t = rbar_T` has vector `1_{[5]-T}`.

*Proof.* Let `sigma` be any element of `Gal(Qbar/Q)`. It maps a boundary place `p` of `P^1` to `sigma(p)`, and the
valuation vector of `sigma(p)` is `chi(sigma) v(p)`, `chi` the character of `Q(i)` (the twisted Galois action on
cocharacters; `split.md` 1.4, referee §4). If `sigma` fixes `i`, `v(sigma p) = v(p)`; since the 30 vectors are
distinct rays, `sigma p = p`. So every place is `Q(i)`-rational, `p = r in P^1(Q(i))`, and complex conjugation maps the
place with vector `1_A` to the place with vector `1_{A^c}`; for `A = T` a subset of `[4]` write `r_T` for it, so the place
`rbar_T` carries `1_{T^c}`. A translation moves the branch pair to `+-sqrt m`; the stabiliser of `{+-sqrt m}` in `PGL_2(Q)`
(the infinite group of norm-one elements of `Q(sqrt m)`) then moves all 30 places off `infinity`. Now
`u_j in Q(i)(t)` has a zero exactly at the places with `v_j - v_5 = 1`, i.e. `1_A` with `j in A`, `5 notin A`:
these are the `r_T` with `T` containing `j`; and poles exactly at `1_A` with `5 in A`, `j notin A`: the `rbar_T`, `T`
containing `j`. All orders are 1. Hence `u_j = c_j P_j/Pbar_j`, and `N(u_j) = 1` gives `N(c_j) = 1`, so
`c_j = gamma_j/gammabar_j` (Hilbert 90 in `Q(i)`). Distinctness is ray-genericity. `□`

**1.4 Lemma 2 (branch condition; proved).** Let `q(t) = t^2 - m`. For (1.1),
`mu_0 >= 4` at both roots of `q` iff

    Im(gamma_j P_j)(t) = kappa_j q(t)^4   for some kappa_j in Q,  j = 1..4,                          (1.2)

where `Im` is taken coefficientwise; `kappa_j != 0` (else `u_j == 1`). Then `mu_0 = 4` exactly. Indeed `u_a/u_b - 1`
has degree 8 (the 8 cuts separating `a, b`), and it vanishes to order `>= mu_0` at both branch points, so `2 mu_0 <= 8`.
Equality means that `u_a/u_b - 1` has **no other zeros** (V-null).

*Proof.* `u_j - 1 = 2i Im(gamma_j P_j)/conj(gamma_j P_j)`, and `conj(gamma_j P_j)` does not vanish at `+-sqrt m`
(those are not boundary places). So `ord_{sqrt m}(u_j - 1) >= 4` iff `(t - sqrt m)^4` divides the real polynomial
`Im(gamma_j P_j)` (coefficients in `Q`); by conjugation over `Q` the same holds at `-sqrt m`, so `q^4` divides it, and
both have degree 8. `□`

**1.5 Lemma 3 (the s-form; proved).** Over `R` put `s = (sqrt m + t)/(sqrt m - t)` (real Möbius, upper half plane to
itself, `sqrt m -> infinity`, `-sqrt m -> 0`), `rho_T = s(r_T)`, `R_j = {rho_T : T contains j}`, `Q_j(s) = prod_{R_j}(s - rho)`.
Then (1.2) holds iff for every `j`

    e_k(R_j) in R  for  k in {1,2,3,5,6,7,8}    (all coefficients of Q_j real except that of s^4).       (1.3)

Equivalently `Im p_k(R_j) = 0` for `k = +-1, +-2, +-3` (`p_k` = power sums) and `prod R_j in R`.
*Proof.* `(t - r)/(t - rbar) = ((s - rho)/(s - rhobar)) (c rhobar + d)/(c rho + d)` for real Möbius coefficients, so
`u_j = c'_j Q_j/Qbar_j`, `|c'_j| = 1`. Contact 4 at `s = infinity` means `deg(c'_j Q_j - Qbar_j) <= 4`, contact 4 at `0`
means `s^4 | (c'_j Q_j - Qbar_j)` (`rho != 0`), so `c'_j Q_j - Qbar_j = kappa' s^4`; the leading coefficient forces
`c'_j = 1` after the phase is absorbed (`gamma` real), which is (1.3). The power-sum form uses Newton's identities at
`infinity` and for `1/rho` at `0`, plus `e_8 = prod R_j`. `□`

The Galois conjugation `tau: sqrt m -> -sqrt m` (fixing `i`) acts by `tau(rho_T) = 1/rho_T`; so `Q`-curves are the
solutions of (1.3) with `rho_T in Q(i, sqrt m)` and `rho_T tau(rho_T) = 1`.

**1.6 Pair structure (proved).** With `g_j = gamma_j P_j/kappa_j = a_j + i q^4` (`a_j` real), `g_j - g_k = a_j - a_k` is
real and vanishes at `r_T`, `T` containing `j,k` (4 points) and their conjugates, so `a_j - a_k = c_{jk} prod_{T
contains j,k} |t - r_T|^2` with `c_jk in Q^*`. Hence on `R` the functions `f_j = a_j/q^4` never cross: the five points
`u_1(t),...,u_4(t), 1` keep their cyclic order on the unit circle for all real `t` and collapse to `1` only at `+-sqrt m`.
In Cayley coordinates (`f = cot` of the half-angle) the condition for the pair `a, b` reads `(A_b - A_a) | (A_a A_b + q^8)`
(tangent addition formula). S_5-symmetrically: for each of the 10 pairs `a,b` of `[5]`, `Im(e^{i phi_ab} S_{a|b}) in R q^4` with
`S_{a|b} = prod_{A contains a, not b}(t - r_A)` (degree 8). These are the pair (bounded-residue) conditions; they are
consequences of (1.2).


## 2. Necessary conditions

All statements in this section are **proved** for every curve of class 1.2 (in the s-form of Lemma 3, `kappa_j > 0`
after replacing `gamma_j` by `-gamma_j` if needed; this does not change `u_j`).

**2.1 Parity and rationality of places.** Lemma 1: each cut carries exactly two places, `r_T` and `rbar_T`, both
`Q(i)`-rational. There is no freedom of orientation (unlike the split case, where each cut is hit once and one must
choose which side carries the place).

**2.2 Half-planes.** `Im Q_j(s) = kappa_j s^4 >= 0` on `R`, so `arg Q_j(s)` stays in `[0, pi]` along `R`; by the argument
principle each `R_j` has exactly **4 points in each half-plane** (the same as toric §7(i)). Consequently the
half-plane sign pattern `epsilon_T = sign Im rho_T` satisfies `sum_{T contains j} epsilon_T = 0` for `j = 1..4`.
There are 278 such patterns, in 15 classes under `S_4 x {global flip}` (`signpat.py`, exact enumeration).

**2.3 The argument sum.** Following `arg Q_j` along `(-infinity, 0)`: each root contributes exactly its principal
argument, so `sum_{rho in R_j} Arg rho = arg Q_j(0) in {0, pi}` **as a real number** (not only mod `pi`), since
`arg Q_j` never leaves `[0, pi]`.

**2.4 Pair structure and ordering (1.6).** For `j != k`, `a_j - a_k = c_jk prod_{T contains j,k} |t - r_T|^2`, with no real
zeros. So the real functions `f_j = a_j/q^4` are strictly ordered on `R - {+-sqrt m}`, the five circle points never
collide for real `t != +-sqrt m`, and each `u_j(t)` avoids `1` there. These are exactly the record's bounded-residue
(V-null) conditions (`toric.md` §7); in the twisted setting they follow from (1.2).

**2.5 What is not a restriction.** The conditions above leave a nonempty open set of configurations (they are open
conditions on the half-plane pattern and on arguments). No further real obstruction was found; §5 records what the
search shows.


## 3. The exact polynomial system and its expected dimension

**3.1 The system (exact).** Unknowns: `rho_T in C` for the 15 nonempty `T` in `[4]`. Equations (28 real, polynomial):

    Im e_k({rho_T : T contains j}) = 0,    j = 1..4,  k in {1,2,3,5,6,7,8}.                         (3.1)

Equivalently (Newton, for nonzero roots) `Im p_{+-1}, Im p_{+-2}, Im p_{+-3}` of `R_j` vanish and `prod R_j` is real.
Gauge: `s -> lambda s`, `lambda in R^*` (the stabiliser of the branch pair). Nondegeneracy: `Im rho_T != 0`,
the 30 points `rho_T, rhobar_T` distinct, `kappa_j = Im e_4(R_j) != 0`.
The `t`-form is (1.2): unknowns `r_T`, `gamma_j`, `m`; complexified, the conjugate roots `r'_T` are independent and
(1.2) reads `e_k(R_j) = e_k(R'_j)` for `k != 4`.

**3.2 Expected dimension (heuristic count; checked against three known cases).**

| problem | unknowns (real) | equations | gauge | expected dim |
|---|---|---|---|---|
| twisted, two conjugate branches, `mu_0 = 4`, `B.C = 30` (this class) | 30 | 28 | 1 | **1** |
| twisted, one rational branch, `mu_0 = 8` (`toric.md` §7) | 30 | 28 | 2 | 0 |
| twisted `k = 4`, two branches, `mu_0 = 2`, `B.C = 14` | 14 | 9 | 1 | 4 |
| twisted `k = 4`, one branch (`gauss` Thm 4.1) | 14 | 9 | 2 | 3 (the record's 3-dim family) |
| split, two branches, `mu_0 = 2`, `B.C = 15` (`split` Thm 2) | 15 | 12 | 1 | 2 (Thm 2's 2-parameter family) |

The split and `k = 4` rows reproduce the known dimensions exactly, so the count method is calibrated. Moreover, at every
solution found for the neighbouring `k = 5` problems with lower contact, `(5,0)`, `(4,1)` and `(4,2)`, the Jacobian has
full rank and the local dimension equals the count (`dimcheck.py`: 14, 14, 10 including gauge). The same holds at the
twisted `k = 4` two-branch solutions (4) and at the doubled `k = 4` solutions (`dbl4.py`: 7). The equations are
transversal wherever they could be tested (**numerical**). **The twisted
analogue of Theorem 2 is expected to come in 1-parameter families of curves**, whose union is a **surface** in `T_5`,
not a hypersurface. In the `(P^1)^4` picture (Cayley coordinates `f_j`, in which `T_5 = (P^1 - {+-i})^4` over `Q`) the
same count reads: curves `s -> (f_j)` with `f_j in <s^{-4},...,s^4>` (35 parameters after the `C^*` gauge) through the
points `(-+i)^4` (3+3 conditions), 8 lines (16) and 12 planes (12): `35 - 34 = 1`.

**3.3 Why the twisted count drops (the anisotropy in the count).** In the split two-branch problem each coordinate
pays 3 conditions (contact 2 on degree-4 blocks) for 15 places; in the twisted one, parity doubles the places (30) and
the degree (8), and contact 4 costs 7 conditions per coordinate. Net: split `15 - 12 - 1 = 2`, twisted
`30 - 28 - 1 = 1`.

**3.4 The `x`-structured ansatz (the literal twisted Theorem 2; heuristic count).** Theorem 2's curve has three coordinates
pulled back from the `x`-line of the conic and one conic block `F_1 = alpha + beta z`. The twisted analogue: the projection
forgetting `u_4` factors through `x` (degree 2), so it is a `k = 4` twisted curve with a rational branch at `x = infinity`,
`mu = 4`, `B.C_x = 14`, i.e. a member of `gauss` Thm 4.1's 3-dimensional family; and `u_4 = H/Hbar` with
`H = alpha(x) + beta(x) z`, `deg alpha <= 4`, `deg beta <= 3`. Contact 4 at both points over `x = infinity` forces
(**proved**: `Im alpha + (Im beta) z` must be a regular function on the affine conic with no pole at infinity, hence constant)

    beta real,   alpha = alpha_r + i kappa  (alpha_r real, kappa constant),

and the norm condition `H H^tau = alpha^2 - beta^2 q(x) = (x - rho_4) prod_S (x - rho_S)/c` with the 7 roots of
the Thm 4.1 curve. The leading coefficient forces `c` real; the `x^7, x^6, x^5` coefficients of `Im` then cost 3 real
conditions on 2 unknowns (`rho_4`): codim 1. The real part must be `-c beta^2 q` up to the free parameter `s = c kappa^2`,
i.e. a degree-8 real polynomial must have the shape (cubic)^2 (quadratic): codim 3 - 1 = 2 (after removing the gauges
`H -> sH` and `(beta, q) -> (beta/lambda, lambda^2 q)`). Total codim 3 on the
3-dimensional Thm 4.1 moduli: **expected dimension 0** (isolated curves). In the split case the same ansatz is not
overdetermined (`beta` may have degree 1 and `alpha` keeps 3 free coefficients: codim 0, giving Theorem 2's 2-dimensional
family). The twisted contact condition is what kills the freedom: in the twisted world `beta` must be real and
`alpha` real up to a constant.


## 4. Reformulations used by the search (all proved) and the search methods

**4.1 Column form (S_5-symmetric).** Index the 30 places by `A` (proper nonempty subsets of `[5]`), with
`rho_{A^c} = rhobar_A`, and put `C_a = prod_{A contains a}(s - rho_A)` (degree 15). Then the class-1.2 conditions are
equivalent to

    C_a - C_b = kappa_ab s^4 G_ab,   G_ab = prod_{A contains a,b}(s - rho_A)  (degree 7),   all a != b,          (4.1)

i.e. the five column polynomials agree in the coefficients of `s^14, s^13, s^12` and `s^3, ..., s^0`; equivalently the
five **column moment vectors** `sum_{A contains a} rho_A^k`, `k = +-1, +-2, +-3`, coincide and the five column products
coincide (a five-column two-sided Prouhet–Tarry–Escott condition). *Proof.* `u_a/u_b = C_a/C_b` (the phases are 1 because
`u_a/u_b -> 1` at `s = infinity` and the `C` are monic); contact 4 at `0` and `infinity` with no further zeros of
`u_a/u_b - 1` (degree count: `deg(u_a/u_b) = 8 = 4 + 4`) is exactly (4.1), the cofactor being forced to be a constant by
degree. `□` Also `C_a Cbar_a = n(s) = prod_{all 30}(s - rho_A)` for every `a`: the `C_a` are five polynomial lattice
points on the polynomial circle `|z|^2 = n(s)`.

**4.2 Half-angle (star) form.** Label each cut by its small side (5 one-sets, 10 two-sets: the vertices and edges of
`K_5`) and put `z_a = (s - rho_a) prod_{x != a}(s - rho_{ax})` (degree 5: the star of `a`). Then
`C_a = z_a Omegabar / zbar_a` with `Omega = prod` of the 15 small-side factors, hence `u_a/u_b = (z_a zbar_b)/conj(z_a zbar_b)`, and

    Im(z_a zbar_b) = kappa_ab s^4 |s - rho_ab|^2   for all pairs a,b                                            (4.2)

(the real factor `|s - rho_ab|^2` divides `z_a zbar_b` automatically). Equivalently (expanding `log z_a` at `0` and `infinity`) the star sums
`sum_{rho in star(a)} Im rho^k` (`k = +-1, +-2, +-3`) are independent of `a` and the star argument sums agree mod `pi`.
All five stars have the same number of upper-half-plane points (argument principle applied to (4.2)).

**4.3 Grassmannian / Möbius form.** Read `z_a(s)` as a vector in `R^2`; (4.2) says the `2 x 5` real polynomial matrix
`Z(s) = [z_1 ... z_5]` has Plücker coordinates `kappa_ab s^4 N_ab(s)`, `N_ab = |s - rho_ab|^2`. Smith normal form gives
`Z = M(s) Zhat(s)` with `det M = c s^4` and `Zhat` of Plücker degree 2; in the balanced splitting type `Zhat = A s + B`, i.e.
five real Möbius maps `zhat_a`. Writing `(g_1, g_2) = adj(M)(1, i)^T` one gets **`z_a = -i det(zhat_a, g)`** and
`Im(g_1 gbar_2) = -c s^4`. Conversely, for **any** real linear `zhat_a` and any `g` with `Im(g_1 gbar_2) = -c s^4`, the
`z_a := -i det(zhat_a, g)` satisfy (4.2) with `kappa_ab s^4 N_ab = c s^4 det(zhat_a, zhat_b)`: contact `(4,4)` is then
automatic, and **the only remaining condition is the class**: `z_a` and `z_b` must share a root for every pair. Geometrically:
`gamma = [g_1 : g_2]` is a degree-4 map with `gamma(R)` touching `RP^1` only at `gamma(0)`, `gamma(infinity)` (order 4), and
its graph must pass through one of the two (conjugate) coincidence points of `zhat_a, zhat_b` for each of the 10 pairs
(10 complex interpolation conditions for a 9-dimensional space of degree-4 maps). All pairwise ratios `zhat_b^{-1} zhat_a`
are elliptic. The other splitting type (row degrees `(0,2)`) is a further closed subcase, not searched separately. Count: `9 - 2 - 6 = 1`, as in 3.2.

**4.4 Projections.** Forgetting one coordinate gives a doubled `k = 4` curve: class 28 (every `k=4` ray hit twice), contact
`(4,4)`, V-null. Conversely a `k = 5` curve is a doubled `k = 4` V-null curve (on `{1,2,3,5}`) plus a choice of one place in
each of the 7 pairs (the one labelled `B u {4}`) plus `rho_4`, subject to (1.3) for `R_4`: 5 net conditions on the 6-dimensional
doubled family (`dbl4.py`: such doubled `k=4` curves exist, are not pull-backs under `s -> s^2`, and have local dimension
exactly `7 = 6 + gauge`, numerically).

**4.5 Search methods (all in `search_checks/`).**
1. Levenberg–Marquardt on (3.1) in several parametrisations: `t`-form with phases (`num_lm.py`), eliminated `t`-form
   (`num_red.py`, 21 unknowns), `s`-form power sums (`sform.py`), half-plane sign-pattern parametrisation
   `rho_T = a_T + i eps_T e^{b_T}` with analytic Jacobian (`fastK.py`, `fastK_run.py`), with deflation (`num_defl.py`).
2. All single-involution symmetric ansätze `(g, pi, eps)`, `g in {-s, 1/s, -1/s}`, `pi in S_5`, `eps = +-1` (`symm.py`,
   `symm_lm.py`).
3. Complexified Newton (conjugate roots independent), with and without Rabinowitsch variables for `kappa_j != 0`
   (`cnewton.py`, `cnewton_rab.py`), and in the complexified Möbius form (`cmob.py`).
4. Staged continuation through intermediate contact orders `(4,0) -> (4,1) -> (4,2) -> (4,3) -> (4,4)`
   (`staged.py`, `staged2.py`, `testK.py`) and a pseudo-arclength homotopy from nondegenerate `(4,2)` solutions,
   `H(x, tau) = [E_{(4,2)}(x); E_{5,6}(x) - (1 - tau) E_{5,6}(x_0)]` (`homot.py`).
5. Lifting doubled `k = 4` V-null curves (4.4) with all `2^7` place choices, then polishing (`lift.py`, `lift2.py`).
6. Möbius form 4.3: Levenberg–Marquardt over five elliptic Möbius maps, the interpolation defect and the shape of
   `Im(g_1 gbar_2)` (`mob.py`, `mob_lm.py`), and a pseudo-arclength homotopy in the `M`-form, from random
   `(M, zhat)` with `det M = s^4` to the root-sharing conditions (`mhomot.py`).
7. A targeted search of the `x`-structured subclass 3.4 as a square system (23 equations, 23 unknowns: the Thm 4.1 curve,
   `rho_4`, `sigma`, `beta`, `q`) (`xstruct.py`).
8. **Calibration** of every solver on problems with known solutions: twisted `k = 4` two-branch (`redk.py`: 57/200
   nondegenerate), doubled `k = 4` V-null (`dbl4.py`), and the `k = 5` problems with lower contact `(3,3)`, `(4,0)`, `(4,1)`,
   `(4,2)`, `(5,0)` (all solved, see 5.1). A nondegenerate solution means: residual `< 1e-10`, all 30 places at distance
   `> 1e-3` from each other and from `R`, `0 < |rho| < 1e4`.


## 5. Results

**5.1 Calibration: the same solvers do find solutions of the neighbouring problems that have them (numerical).**

| problem (contact at `infinity`, `0`) | expected dim | solver | result |
|---|---|---|---|
| twisted `k = 4`, two branches `(2,2)`, class 14 | 4 | LM, eliminated form (`redk.py`) | 57 nondegenerate / 200 starts |
| doubled `k = 4` (class 28), V-null `(4,4)` | 6 | LM (`dbl4.py`) | found; local dim exactly 7 = 6 + gauge (Jacobian rank 21 of 28) |
| doubled `k = 4`, V-null `(4,4)` | 6 | **homotopy** from `(4,2)` (`homot_dbl.py`) | 49 paths: 5 reach `tau = 1`, **1 nondegenerate** |
| `k = 5`, class 30, `(3,3)` | 9 | LM (`kset_run.py`) | 5 / 15 starts |
| `k = 5`, `(4,0)` | 16 | LM, analytic Jacobian (`fastK_run.py`) | 14 / 30 |
| `k = 5`, `(4,1)` | 13 | same | 31 / 400 (also reached by staging) |
| `k = 5`, `(4,2)` | 9 | same | 1 / 400 direct; routinely by staging `(4,0)->(4,1)->(4,2)` |
| `k = 5`, `(5,0)` | 12 | same, unrestricted half-plane patterns | 29 / 400 |
| **split** `k = 5`, class 15, `(2,2)` (split Thm 2 class, known 2-dim family) | 2 | same LM, real places, all 278 orientation patterns (`split_calib.py`) | 2 / 400 (387 converged) |

Two lessons from the calibration: direct LM success falls quickly with the number of conditions (at `(6,1)`, expected
dim 5, it found 0 / 400), so a direct-LM failure alone means little; and the half-plane split `4/4` is forced only when
`Im(gamma Q_j)` has constant sign (contact `(4,4)` or `(8,0)`), so restricting to `4/4` patterns is legitimate exactly in
the target case.

**5.2 The target class (contact `(4,4)`, class 30): no nondegenerate curve found (numerical evidence).**

| method | runs | converged | nondegenerate |
|---|---|---|---|
| LM, `t`-form with phases (`num_lm.py`) | 40 | 7 | 0 |
| LM, eliminated `t`-form (`num_red.py`) | 3 x ~1500 | 11 | 0 |
| LM, `s`-form power sums + slice (`sform.py`) | 40 | 36 | 0 |
| LM, `s`-form, analytic Jacobian, `4/4` patterns (`fastK_run.py 4 4`) | 400 | 251 | **0** |
| LM with deflation (`num_defl.py`, `sform.py` with `eps > 0`) | ~3000 | 0 | 0 |
| symmetric ansätze, 14 nontrivial `(g, pi, eps)` classes x 40 (`symm_lm.py`) | 560 | 249 | 0 |
| complexified Newton (`cnewton.py`) | 30 | 21 | 0 (all `kappa = 0`) |
| staged LM `(4,0) -> ... -> (4,4)` (`staged.py`, `staged2.py`) | ~50 paths | `(4,2)` reached often | 0 at `(4,3)`/`(4,4)` |
| lifting doubled `k = 4` V-null curves + all `2^7` choices + polish (`lift.py`, `lift2.py`) | 96 + 84 polishes | all | 0 |
| pseudo-arclength homotopy from nondegenerate `(4,2)` points (`homot.py`, `homot_c.out`) | 165 paths | 17 reach `tau = 1` (all at exact coincidences, separation `<= 3e-15`); 46 degenerate earlier; 36 stall; 63 turn back; 3 step limit | **0** |

Every converged point is degenerate. The degenerations are of three kinds: (a) `kappa_j = 0` (some `u_j == 1`), (b) places
on the real line (a cut disappears), (c) coincidences of places. In the homotopy, **every path that reaches `tau = 1` ends
at an exact coincidence configuration** (separation `~1e-16`), and these are *collapsed* configurations: the merged
places all carry the same vector `+-1_{a}` (e.g. all 30 places merge in pairs `{2,345}, {245,3}, {12,1345}, ...` whose
vectors are all `+-1_1`), i.e. four of the five points coincide identically (all 17 endpoints, checked by `collapse_check.py`; clusters whose
vectors sum to `0` are not places). Paths that do not reach `tau = 1` either
turn back or stall with two or more places converging as `tau -> 1`. In the calibration problem (doubled `k = 4`, same
homotopy, same tracker) nondegenerate endpoints *are* reached.

**5.3 Reading.** *Evidence, not proof:* the polynomial system (3.1) has, near every nondegenerate configuration the
searches could reach, only degenerate solutions; its solution set appears to consist of the collapsed and boundary
components. This is the two-branch analogue of `toric.md` §7's single-branch finding, now with a calibrated solver,
a homotopy, five different parametrisations and the symmetric subcases. I could not prove non-existence: the reformulations
4.1–4.3 make the problem a clean incidence problem (five real Möbius maps and one degree-4 map), but I found no
obstruction (the count is `1`, the sign and argument constraints of §2 are satisfiable, the Ptolemy/Plücker relations are
consistent).

**5.4 Family dimension (if curves exist).** By 3.2 the expected dimension of the family of such curves is **1**, so they
would sweep at most a **surface** in `T_5` (heuristic count, calibrated in 3.2 against three known families, but not a
theorem: an excess-dimensional component is not excluded). A surface, or finitely many, can be put into `Z`; only
infinitely many hypersurfaces would bear on refuting Question B. For comparison, the split analogue sweeps a hypersurface
(`split` Thm 2 (v), refereed). Arithmetic: a curve over `Q` in the class needs, in the `s`-form, `rho_T in Q(i, sqrt m)` with
`rho_T tau(rho_T) = 1`; so even a real 1-parameter family would contribute `Q`-curves only through the rational points of
its (twisted) parameter curve.

**5.5 The literal twisted Theorem-2 ansatz.** The `x`-structured subclass (3.4) is overdetermined by 3 on the
3-dimensional `gauss` Thm 4.1 family (expected dimension 0); the variant with the double cover ramified at the branches
(`x = s^2`; three coordinates pulled back from a `k = 4` two-branch `(2,2)` curve) is overdetermined by 1 (count in
`search_notes.md` N15). Both subclasses lie inside the class searched above, but being of expected dimension `<= 0`
they were not targeted separately; random-start searches of the full class would not be expected to hit them.


## 6. Bearing on Question B, and the owner-constraint audit

**6.1 What this round adds to the map.**
* *The twisted analogue of split Theorem 2 is a different object.* The parity lemma doubles the weights; Lemma 1 then
  forces all 30 places to be `Q(i)`-rational points of `P^1_t`, so a "Pell conic with Gaussian quadratic blocks" is
  nothing more than `P^1` with Gaussian **linear** blocks and a real-quadratic branch pair `+-sqrt m`. The class is the
  explicit system (1.2)/(3.1).
* *Count:* the split two-branch problem has expected (and actual) dimension 2, which is why split Theorem 2 sweeps a
  hypersurface `W`; the twisted one has expected dimension **1** (3.2–3.3). The literal transplant of Theorem 2's structure
  (three coordinates through the conic's `x`-line) is overdetermined by 3 (3.4). So even before any search, the
  anisotropy predicts **no twisted `W`**: at most surfaces. This is the first place where the anisotropy enters a
  dimension count rather than only a parity.
* *Search:* in the full class no curve was found by any of the methods of 4.5, all calibrated (5.1–5.2). Numerically the
  solution set of (3.1) near everything reachable consists of collapsed configurations. Together with `toric.md` §7
  (one rational branch, `mu_0 = 8`; re-run here as `fastK_run.py 8 0`: 271 converged, 0 nondegenerate) this is consistent
  with: **the balanced twisted Liouville constant `15/4` is not attained by rational curves of the minimal class through `e`,
  with one or with two branches.** That is a conjecture supported by evidence, not a theorem.
* *Refutation side:* nothing here refutes Question B. A refutation needs a Zariski-dense supply of balanced points at
  exponent `>= 7/15`; families of the searched type (if they existed) sweep at most surfaces, and the searched class is
  empty in every computation. If one tries to *refute* `A_5` by Pell-type curves, the twisted form leaves much less room than the
  split one.
* *Proof side:* nothing here proves Question B either. Absence of minimal-class curves at `15/4` says nothing about
  sporadic points or about curves of higher class (`B.C = 60, 90, ...` with `mult_e = 16, 24, ...`), and the count 3.2
  does not exclude excess components. The reformulations 4.1–4.3 (column PTE, star form, five Möbius maps plus a degree-4
  map) are the natural place to look for a proof of emptiness.

**6.2 New precise questions.**
> **Q24.1 (emptiness).** Are there 15 points `rho_c in C - R` (indexed by the vertices and edges of `K_5`), pairwise
> distinct and distinct from all conjugates, such that the five star sums `sum_{rho in star(a)} Im rho^k`,
> `k = +-1, +-2, +-3`, are independent of `a` and the star argument sums agree mod `pi`? (Equivalently: is the class of 1.2
> empty over `R`?) Evidence says no. A **no** means that no twisted rational curve of the minimal balanced class reaches
> `15/4` through `e` with a real-quadratic branch pair.
>
> **Q24.2 (Möbius form).** Five real Möbius maps with pairwise elliptic ratios, and a degree-4 map `gamma` whose real
> trace touches `RP^1` only at two points, to order 4 each. Can the graph of `gamma` pass through one coincidence point
> of each of the 10 pairs?

**6.3 Owner-constraint audit.**
* *No uniform exceptional set from a fixed-`S` Subspace theorem* is used or claimed; the only exceptional sets mentioned
  are the hypothetical `Z` of `A_5` and the (non-)existent families.
* *No Gaussian gcd or denominator cost discarded:* Lemma 1 keeps the constants `gamma_j` (Hilbert 90), the heights are the
  toric heights of the record, and `alpha_C = 15/4` is computed with `h_{-K} = 30 h(t) + O(1)`, i.e. with all
  denominators. No point counts or radii are claimed.
* *No drifting common phase replaced by zero:* the phases `gamma_j` are unknowns throughout; they are forced to be real only
  in the `s`-normalisation, by the contact condition at `s = infinity` (Lemma 3), not by assumption.
* *No template-dependent threshold treated as uniform:* nothing here is a bound in the radius.


## 7. Files and checks

All under `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round24/`; the repository was
not edited. Heavy runs used `prlimit --as=3500000000`, at most two at a time. Working notes: `search_notes.md` (N1–N21).
(`notes.md` and `theory*.md` in the same folder belong to the parallel theory agent.)

| file (`search_checks/`) | content | label |
|---|---|---|
| `sym_identities.py` | symbolic: Möbius converse identity `Im(z_a zbar_b) = Im(g_1 gbar_2) det(zhat_a, zhat_b)`; half-angle identity `C_a zbar_a = z_a Omegabar` and `C_a Cbar_a` independent of `a` (exact Gaussian rationals) | **exact** (all True) |
| `signpat.py`, `signpat_all.py` | the 278 half-plane patterns (15 classes); all satisfy the 20 ordered-pair constraints | exact enumeration |
| `common.py`, `num_lm.py` | `t`-form LM with phases | numerical |
| `red.py`, `num_red.py`, `num_defl.py` | eliminated `t`-form (21 unknowns), with deflation | numerical |
| `redk.py` (`k4_calib.out`) | general-`k` eliminated form; `k = 4` calibration 57/200 | numerical |
| `sform.py`, `kset.py`, `kset_run.py`, `testK.py`, `fastK.py`, `fastK_run.py` (`fk_*.out`) | `s`-form, contact-pair scans `(mu_inf, mu_0)` | numerical |
| `symm.py`, `symm_lm.py` (`symm_lm_*.out`) | symmetric ansätze | numerical |
| `cnewton.py`, `cnewton_rab.py`, `cmob.py` | complexified Newton (place form, Rabinowitsch, Möbius form) | numerical |
| `staged.py`, `staged2.py`, `warm.py`, `diag50.py` | staged continuation and coincidence diagnostics | numerical |
| `dbl4.py`, `lift.py`, `lift2.py` | doubled `k = 4` V-null curves (local dim 7) and lifting to `k = 5` | numerical |
| `homot.py`, `homot_dbl.py`, `mhomot.py` (`homot_c.out`, `homot_dbl.out`, `mhomot_a.out`) | pseudo-arclength homotopies (target, calibration, Möbius form) | numerical |
| `dimcheck.py` (`dimcheck.out`) | Jacobian rank at found solutions equals the count of 3.2 (`(5,0)`, `(4,1)`, `(4,2)`, `k=4`) | numerical |
| `collapse_check.py` | clusters the coincident places of each homotopy endpoint and sums their valuation vectors: all 17 endpoints are 4-point collapses | exact bookkeeping on numerical endpoints |
| `split_calib.py` (`split_calib.out`) | the same LM on the split Thm-2 class: 2 nondegenerate / 400 | numerical |
| `mob.py`, `mob_lm.py` | Möbius form 4.3: interpolation defect for random elliptic data (`sigma_min/sigma_max ~ 1e-6..1e-3`, never 0) and LM | numerical |

