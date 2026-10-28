# Round 24: twisted Pell-type search (two conjugate real-quadratic branches, B.C = 30, balanced, alpha = 15/4)

## 0. Verdict

**No growth improvement and no uniform bound.** `M <= (2/3 + eps) log R/loglog R` stands. Status: **anisotropy gain (map
value)**.
* The searched class **exists over `R`**: an explicit one-parameter family, sweeping a surface.
* Only **finitely many** of its members are defined over `Q` (proved).
* The class searched (round 23 `split.md` §6.2 item 3) is: rational curves through `e` on the twisted Losev–Manin 4-fold, of
  the balanced minimal class (`E_A.C = 1` for all 30 rays, `B.C = 30`, every cut of weight 2), with a pair of conjugate
  real-quadratic branches of contact `mu_0 = 4`, hence `alpha_C = 15/4`.

1. **Normal form (proved, §1).** Galois forces all 30 boundary places to be `Q(i)`-rational points of `P^1_t`. So the
   ansatz "Pell conic with Gaussian blocks" is `P^1` with 15 Gaussian **linear** blocks `t - r_T`, `r_T in Q(i)`, and
   branch pair `+-sqrt m`. The exact system is `Im(gamma_j P_j) = kappa_j (t^2 - m)^4` (`j = 1..4`). After a real Möbius
   change it becomes (1.3): each 8-point set `R_j` has `e_k` real for `k != 4`. Three further exact reformulations are in §4:
   a five-column Prouhet–Tarry–Escott form, a half-angle form on `K_5`, and a Möbius/Grassmannian form.
2. **Expected dimension 1 (§3, calibrated).** The split analogue has 2 (Theorem 2's hypersurface `W`). So twisted families
   sweep at most surfaces.
3. **Theorem A (§5.1; exact algebra plus 50-digit verification of nondegeneracy).** The literal transplant of split Thm 2
   has real solutions. They form an explicit family `M` parametrised by `m in (-3, -1)`.
   * Its structure: three coordinates through the conic's `x`-line, plus a reflection symmetry.
   * The 8 base points are `p_r = t/b_r + i b_r` and their reflections `-pbar_r`.
   * The `b_r` are the roots of `W^4 + v W^2 - v^2 W + u v^2`, where
     `u = -(m+3)c/(16(m-3)^2(m-1)^2)`, `v = -(m+1)c/(4(m-3)^3(m-1))` and `c = m^3 - 11m^2 + 39m - 13`; `t^2 = prod b_r`.
   * Every member hits each of the 30 rays exactly once and has contact 4 at both real branches.
   * `M` is a whole component of the class: the Jacobian has rank 28 there. It sweeps a real **surface** in `T_5`.
   * All blind searches missed it (§5.4). Their earlier reading as "evidence of emptiness" is **withdrawn**.
4. **Theorem B (§5.3, proved).** `M` contains **only finitely many curves defined over `Q`**. The proof combines Lemma 1, the
   rationality of the invariant curve, an exact discriminant with 8 odd-order points, Riemann–Hurwitz (genus `>= 3`) and
   Faltings. An exact search finds none with `|num|, |den| <= 60`.
   * The split family is rational in its places, which forces `Z ⊇ W`.
   * The twisted family forces at most finitely many curves into `Z`. The anisotropy enters through the rationality of the
     30 places: the `b`-quartic must split completely over the rational `m`-line.
5. **Bearing on Question B (§6).** Nothing is refuted and nothing is proved. The Pell mechanism that forces a hypersurface into
   `Z` in the split form does not transfer to `T_5`, at least not via `M`. Open: other components (Q24.1), the rational
   members of `M` (Q24.2), and higher classes (Q24.3).

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
  balanced but were not searched separately, and none occurred among the computed limits, §5.4);
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
conditions on the half-plane pattern and on arguments). No further real obstruction exists in general: §5.1 exhibits real solutions.


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

**3.4 The `x`-structured ansatz (the literal twisted Theorem 2).** Theorem 2's curve has three coordinates
pulled back from the `x`-line of the conic and one conic block `F_1 = alpha + beta z`. The twisted analogue: the projection
forgetting `u_4` factors through `x` (degree 2, unramified over the branches), so it is a `k = 4` twisted curve with a
rational branch at `x = infinity`, `mu = 4`, `B.C_x = 14`, i.e. a member of `gauss` Thm 4.1's family; and `u_4 = H/Hbar`
with `H = alpha(x) + beta(x) z`, `deg alpha <= 4`, `deg beta <= 3`, on the conic `z^2 = q(x)`. Contact 4 at both points over
`x = infinity` forces (**proved**: `Im alpha + (Im beta) z` is a regular function on the affine conic without pole at
infinity, hence constant)

    beta real,   alpha = alpha_r + i kappa  (alpha_r real, kappa constant),

and the zeros of `H` must lie over the 8 `x`-values of the places: `H H^tau = alpha^2 - beta^2 q = (x - rho_4) prod_S (x - rho_S)/c`.
Writing `n = (x - rho_4) prod_S (x - rho_S)` and `sigma = c kappa^2`, this is: `c` real, `Im n` of degree `<= 4`, and

    Re n - (Im n)^2/(4 sigma) + sigma = -c beta^2 q    (a degree-8 real polynomial of the shape (cubic)^2 (quadratic)).   (3.2)

*Naive count:* codim 3 on the 3-dimensional Thm 4.1 moduli, so expected dimension 0. **This count is wrong by one:** the
actual solution set is a curve (§5.1). The reason is a hidden structure, which §5.1 makes explicit: a reflection
symmetry, under which (II) `deg Im n <= 4` holds automatically. In the split case the same ansatz is not
overdetermined (`beta` of degree 1, three free coefficients of `alpha`), giving Theorem 2's 2-dimensional family.

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

**5.1 Theorem A (curves exist over `R`: an explicit one-parameter family).** For a real parameter `m` let

    c(m) = m^3 - 11 m^2 + 39 m - 13,
    u(m) = -(m+3) c(m) / (16 (m-3)^2 (m-1)^2),     v(m) = -(m+1) c(m) / (4 (m-3)^3 (m-1)).

1. Let `b_0, ..., b_3` be the roots of `W^4 + v W^2 - v^2 W + u v^2`; these are the elementary symmetric functions
   `e_1 = 0, e_2 = v, e_3 = v^2, e_4 = u v^2`.
2. Put `t = sqrt(b_0 b_1 b_2 b_3)`, `a_r = t/b_r` and `p_r = a_r + i b_r`.
3. In a real coordinate `X` on the base line, the 8 base points are
   `rho_123 = p_0`, `rho_4 = -pbar_0`, `rho_j = p_j`, `rho_{[3]-j} = -pbar_j` (`j = 1,2,3`).
4. Let `n(X) = prod_r (X - p_r)(X + pbar_r)` and `sigma = -prod_r |p_r|^2`.
5. Then `L(X) = Re n - (Im n)^2/(4 sigma) + sigma` is even, and `L = Y (Y + e)^2 (Y + D)` with `Y = X^2`. Put
   `beta = X(X^2 + e)`, the conic `z^2 = -(X^2 + D)/sigma`, and `H = Im n/(2 sigma) + i + beta z`.
6. Set `u_j = P_j/Pbar_j` with `P_j = prod_{S contains j}(X - rho_S)` (`j = 1, 2, 3`), and `u_4 = H/Hbar`.

For every `m` in the open interval `(-3, -1)` the four `b_r` are real, and this map `conic -> T_5` is a rational curve
over `R` **in the class 1.2**:
* every one of the 30 rays is hit exactly once (`B.C = 30`, all 15 cuts of weight 2);
* its two branches through `e` are the real points at infinity of the conic, and each has contact exactly 4.

*Proof and labels.*
* (a) **Exact** (`xsym.py`, `xsym2.py`). Under the reflection symmetry `rho_{partner} = -conj(rho)`, the Thm 4.1 conditions
  on the `k = 4` part are equivalent to `sum b_r = 0`, `a_r b_r = t` for all `r`, and `t^2 = b_0 b_1 b_2 b_3`.
* (b) **Exact** (`xsym3.py`). With these, `Im n = c_3 X^3 + c_1 X`, so condition (II) of 3.4 holds identically. Also
  `L(0) = 0` exactly for `sigma = -prod |p_r|^2`, and `L` is even.
* (c) **Exact** (`xsym4.py`, `xsym6.py`). `L/Y` has a double root iff a discriminant vanishes. The relevant factor (degree
  24 in `(b_1, b_2)`) is the `S_4`-invariant
  `F = 32e_2^6e_4^3 - 96e_2^5e_3^2e_4^2 + 37e_2^4e_3^4e_4 - 256e_2^4e_4^4 - 4e_2^3e_3^6 + 1024e_2^3e_3^2e_4^3 - 992e_2^2e_3^4e_4^2 + 288e_2e_3^6e_4 - 27e_3^8 - 256e_3^4e_4^3`.
* (d) **Exact** (`xsym7.py`, `xsym8.py`). In `u = e_4/e_2^2`, `v = e_3^2/e_2^3`, the equation `F = 0` is a quintic
  `Phi(u,v) = 0` with a triple point at the origin. The pencil of lines `v = k u` through that point gives the discriminant
  `(k-4)^5 (9k-4)^3`, so `Phi` is birational to the conic `w^2 = (k-4)(9k-4)`, which has the rational point `(4, 0)`. Hence
  `Phi` is **rational over `Q`**, with the parametrisation above (`Phi(u(m), v(m)) = 0` checked symbolically). The other
  branch is `m -> -m`.
* (e) Given (a)–(d), the identity `(alpha_r + i)^2 - beta^2 q = n/sigma` and contact order 4 at both points at infinity
  hold identically. Indeed `u_j - 1 = const/Pbar_j` and `u_4 - 1 = 2i/Hbar`, and both denominators have poles of order 4.
* (f) The class is an open, nondegeneracy condition: 30 distinct places, each with a ray vector. It was verified with
  50-digit arithmetic at `m = -1.1, -1.5, -1.8, -2, -2.2, -2.5, -2.8, -2.95` (`construct_m.py`, `construct_m.out`; every
  output reads "places: 30, rays hit: 30, all once, contact 4.0 at both branches").
  - The member found first, by the targeted search (`xstruct.py`, solution 76), was refined to residual `1e-119`
    (`identify_x.py`) and verified the same way (`verify_x.py`).
  - Converted to the `s`-form, its 15 places satisfy the general system (3.1) to `5e-51` (`xs_to_s.py`).

So existence over `R` is **proved** up to the nondegeneracy check at specific parameters. That check is numerical with
large margins: the minimal separation of the 30 places, in the conic parameter, ranges from `3.2` at `m = -1.1` to
`9e-5` at `m = -2.95`, against a working precision of `1e-50`; places merge only as `m -> -3`. `□`

**5.2 Dimension and sweep.**
* At the member of 5.1, the Jacobian of the 28 equations (3.1) has rank 28 (`xs_to_s.py`). So the solution set is locally a
  curve modulo the gauge, and the irreducible component `M` of the class through it **is** this family. Over `C` it is
  parametrised by the plane curve `C_b = {F = 0}` of projective `b`-quadruples (`sum b = 0`), which is a degree-24 cover of
  the `m`-line (the 24 orderings of the `b_r`).
* The x-structured `23 x 23` system has a 1-dimensional kernel there (a singular value `1.7e-81` at 130 digits). This
  confirms that the naive count 0 of 3.4 is wrong by one.
* The family **sweeps a real surface** in `T_5(R)`: the map (family parameter, curve parameter) `-> (S^1)^4` has rank 2 at
  45/45 sampled points (`xfam_sweep.py`; continuation `xfam.py` along 382 points).
* This is the dimension predicted by the count 3.2. It is **not a hypersurface**, unlike the split `W`.

**5.3 Theorem B (only finitely many members are defined over `Q`; proved).** The component `M` contains only finitely many
curves defined over `Q`.

*Proof.*
1. Let a curve of `M` be defined over `Q`. By Lemma 1 its places are `Q(i)`-rational.
2. The `x`-coordinate, the centring `X = x + q_1/2` (with `q_1 in Q`) and the scaling that normalises the `b_r` are all
   defined over `Q`. Hence `p_r in Q(i)`, so `b_r in Q` (up to a common rational factor) and `t in Q`.
3. Then `(u, v) in Q^2`, so `m in Q` away from the finitely many points where the birational map `Phi -> P^1_m` is not
   defined.
4. Moreover `Y_r := b_r^2/(-e_2)`, which is invariant under scaling and under `b -> -b`, are four rational roots of
   `G_m(Y) = (Y^2 - Y + u)^2 + v Y`.
5. By `genus_Y.py` (exact), `disc_Y G = -2^47 (m-3)^42 (m-1)^34 (m+1)^2 (m+3) c(m)^3 (m^4 + 70 m^2 - 48 m + 9)` up to the
   content. The leading coefficient `256 (m-3)^9 (m-1)^7` is nonzero at the 8 roots of `(m+3) c(m) (m^4+70m^2-48m+9)`,
   and the discriminant has **odd** order at each of them.
6. So the local monodromy of the roots there is an odd permutation, with ramification index at least 2.
7. Let `N != 1` be the group of the geometric Galois closure `X~ -> P^1_m` of `G_m`. Riemann–Hurwitz gives
   `2g(X~) - 2 >= |N| (-2 + 8 * 1/2) = 2|N|`, so `g(X~) >= |N| + 1 >= 3`.
8. A rational `m` at which `G_m` splits completely over `Q` lifts to a rational point of one of the finitely many geometric
   components of the curve of ordered roots, unless `m` is one of the finitely many branch points. Each component is a copy
   of `X~`, of genus `>= 3`.
9. By Faltings there are finitely many such `m`. Each `m` gives at most `24 x 2` curves. `□`

*Exact search* (`qsearch_m.py`): no `m = p/q` with `|p|, |q| <= 60` gives even **one** rational root `b`. So no curve of `M`
over `Q` is known. The real-quadratic condition on the branches (`-1/sigma` not a square) and the condition
`b_0 b_1 b_2 b_3 in Q^2` would cut the finite set down further.

**5.4 What the blind searches showed, and why they missed `M`.**

**5.4a Calibration of the blind solvers (numerical).**

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

**5.4b Blind searches of the target class (numerical).**

| method | runs | converged | nondegenerate |
|---|---|---|---|
| LM, `t`-form with phases (`num_lm.py`) | 40 | 7 | 0 |
| LM, eliminated `t`-form (`num_red.py`) | 3 x ~1500 | 11 | 0 |
| LM, `s`-form power sums + slice (`sform.py`) | 40 | 36 | 0 |
| LM, `s`-form, analytic Jacobian, `4/4` patterns (`fastK_run.py 4 4`) | 400 | 251 | **0** |
| same, second run (`fk_44_big.out`) | 2000 | 1195 | **0** |
| `x`-structured square system (`xstruct.py`), targeted | ~225 | 84 | 5 passed the filter: **1 genuine (= `M`)**; 4 degenerate on full check (16 places, no contact at the second branch; `verify_x.py`) |
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

**Reading (corrected).** All the blind methods of 4.5 (LM in five parametrisations, symmetric ansätze, complexified Newton,
staging, lifting, homotopies) missed the family `M`, which exists. The homotopy endpoints all collapsed (`collapse_check.py`).
Their failure reflects tiny basins of attraction, not absence. **The provisional "evidence of emptiness" of an earlier draft is
withdrawn.** The targeted square formulation of the literal Theorem-2 transplant (`xstruct.py`) found `M`, once in about 225 starts.
Lesson for the record: in these problems blind searches are not evidence of non-existence, even when calibrated, if the
target set has lower dimension than the calibration problems.

**5.5 Split versus twisted: where the anisotropy enters.**

| | split Thm 2 (round 23) | twisted family `M` (this round) |
|---|---|---|
| structure | three coordinates through the `x`-line, one conic block | the same, plus a reflection symmetry |
| class, branches | `B.C = 15`, two branches of contact 2 | `B.C = 30`, two branches of contact 4 |
| `alpha_C` | `15/4` | `15/4` |
| family over `R` | 2-dimensional, sweeps a hypersurface `W` | 1-dimensional, sweeps a surface |
| parameter space in place coordinates | the places `(a, b, c) in Q^3` themselves (rational) | the curve `C_b` of `b`-quadruples, whose completely split points come from a curve of genus `>= 3` |
| curves over `Q` | Zariski dense in `W`, so `Z ⊇ W` is forced | **finitely many** (Theorem B); none of small height |

The invariants `(u, v)` of `M` form a **rational** curve over `Q`, so the "shape" of the configuration varies rationally, exactly
as in the split case. The anisotropy enters through Lemma 1: the twist forces each of the 30 places to be individually
`Q(i)`-rational, which requires the `b`-quartic to split completely. Over the rational `m`-line that quartic has
nontrivial monodromy, and this is what makes the rational members finite.

## 6. Bearing on Question B, and the owner-constraint audit

**6.1 What this round adds to the map.**
* **The twisted Pell mechanism exists over `R` but dies over `Q`, at least on the component found.**
  - Theorem A: rational curves of the minimal balanced class with two branches of contact 4, i.e. Liouville-extremal at
    `alpha = 15/4`, exist on the twisted 4-fold over `R`.
  - They form a 1-parameter family, the literal transplant of split Theorem 2 plus a reflection symmetry, and sweep a real
    surface.
  - Theorem B: only finitely many of them are defined over `Q`. So this component forces at most finitely many curves into
    an admissible `Z`, never a surface.
  - In the split form the analogous family is rational in its places, and forces the hypersurface `W` into `Z`.
  - This is the first point where the anisotropy of `T_5` changes an **arithmetic** outcome rather than a parity or a count:
    Lemma 1 (all 30 places `Q(i)`-rational) plus `S_4`-type monodromy of the places over the rational invariant curve.
* *Count.* The general class has expected (and, on `M`, actual) dimension 1, against 2 for the split class. So even
  geometrically, twisted Pell-type families sweep surfaces, not hypersurfaces (3.2, 5.2).
* *Refutation side.* Nothing is refuted. A refutation of Question B needs a Zariski-dense supply of balanced points at
  exponent `>= 7/15`. `M` gives finitely many curves over `Q`, if any; none was found up to height 60. Other components of
  the class, or curves of higher class (`B.C = 60, 90, ...` with `mult_e = 16, 24, ...`), are not excluded.
* *Proof side.* Nothing is proved about Question B. Theorem B is specific to `M`.
* *Methodological correction.* Blind numerical searches of this system are **not** evidence of emptiness. All of them
  missed `M` (5.4). `toric.md` §7's single-branch "no nondegenerate solution" should be read with the same caution.

**6.2 New precise questions.**
> **Q24.1 (other components).** Does the class 1.2 over `R` (equivalently, the system (3.1)) have irreducible components
> other than `M` and its `S_5`-images? Each component is a curve, with expected dimension 1 (3.2). Does any component carry
> infinitely many curves over `Q`? That needs its place-curve to have genus `<= 1`.
>
> **Q24.2 (rational members of `M`).** Is there `m in Q` such that `G_m(Y) = (Y^2 - Y + u(m))^2 + v(m) Y` splits completely
> over `Q` and the resulting curve has real-quadratic branches? The set is finite; none has height `<= 60`.
>
> **Q24.3 (higher class).** Balanced twisted curves at `15/4` with `B.C = 60`, `mult_e = 16` (two branches of contact 8):
> expected dimension? Higher-degree families with genus-0 place curves would be the next candidates for forced surfaces.

**6.3 Owner-constraint audit.**
* *No uniform exceptional set from a fixed-`S` Subspace theorem* is used or claimed; the only exceptional sets mentioned
  are the hypothetical `Z` of `A_5` and the (non-)existent families.
* *No Gaussian gcd or denominator cost discarded:* Lemma 1 keeps the constants `gamma_j` (Hilbert 90), the heights are the
  toric heights of the record, and `alpha_C = 15/4` is computed with `h_{-K} = 30 h(t) + O(1)`, i.e. with all
  denominators. No point counts or radii are claimed.
* *No drifting common phase replaced by zero:* the phases `gamma_j` are unknowns throughout; they are forced to be real only
  in the `s`-normalisation, by the contact condition at `s = infinity` (Lemma 3), not by assumption.
* *No template-dependent threshold treated as uniform:* nothing here is a bound in the radius. Theorem B is an
  ineffective finiteness statement (Faltings); no height bound for the rational members of `M` is claimed.


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
| `xstruct.py` (`xstruct_run.out`) | targeted square system for the `x`-structured ansatz; found `M` (solution 76) | numerical |
| `verify_x.py`, `debug_x.py`, `xs_to_s.py`, `identify_x.py` | 50–130-digit refinement; places/vectors/contact of the found curve; conversion to (3.1) (residual `5e-51`) and rank 28; non-algebraicity of low degree; `1`-dim kernel of the `x`-structured system | numerical, high precision |
| `xfam.py` (`xfam_points.npy`), `xfam_rel.py`, `xfam_sweep.py` | continuation of `M` (382 points), numerical relations, sweep rank 2 | numerical |
| `xsym.py`, `xsym2.py`, `xsym3.py` | symbolic: Thm 4.1 part `<=>` `sum b = 0, a_r b_r = t, t^2 = prod b`; `Im n` of degree 3 identically | **exact** |
| `xsym4.py` … `xsym8.py` (`xsym4.pkl`, `f24.pkl`, `F24.pkl`, `uv_branch*.pkl`) | discriminant factorisation; degree-24 factor; its `S_4`-invariant form `F(e_2,e_3,e_4)`; rational parametrisation of `Phi(u,v)` | **exact** |
| `xsym9.py`, `xsym10.py`, `genus_Y.py` | `b`-quartic over `Q(m)`; discriminants; the invariant quartic `G_m(Y)` and its discriminant (Theorem B) | **exact** |
| `construct_m.py` (`construct_m.out`) | builds the curve from `m` and checks the class and contact (8 values of `m`) | exact construction, 50-digit check |
| `qsearch_m.py` | rational `m` with `|p|,|q| <= 60`: no rational `b`-root | **exact** |
| `split_calib.py` (`split_calib.out`) | the same LM on the split Thm-2 class: 2 nondegenerate / 400 | numerical |
| `mob.py`, `mob_lm.py` | Möbius form 4.3: interpolation defect for random elliptic data (`sigma_min/sigma_max ~ 1e-6..1e-3`, never 0) and LM | numerical |

