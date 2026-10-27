# Round 16: Exact reduction of MIN^pair to discrete anti-concentration (AC)

**Status: complete. No growth improvement and no uniform bound is proved.** `M <= (2/3 + eps) log R/loglog R` stands.

## 0. Summary of results

Setting throughout: `T = 1`, `w = sqrt(8(n-1) ln 4r)`, `j = floor(w/2)`, `m' = ceil(4w)`, regime `n - 1 >= 316 ln 4r`,
any profile, any norm; models are named where used.

1. **The exact reduction (Theorem 2.1, Corollary 2.2; proved).** MIN^pair with `kappa_V` in place of `kappa_H` is
   *equivalent* to the sup-norm statement `Cheb_A(u) > d^{m'} eps_0` on the gap sphere `G_1`, and, up to a factor
   `3^{d/2} = e^{0.55 d}` (negligible against the room `e^{(3/4 - o(1)) n ln n}`), to the `L^2` statement
   `Psi_A(u) > d^{2m'} eps_0^2`, where `Psi_A(u)` is the squared `L^2(mu_A)` norm of the **top orthogonal component** of
   `l_u^{m'}` on `A`. All probabilistic versions (any probability threshold, "all or most `k`", uniform on `A` or the cube
   with threshold `mu(A)`) are the same statement up to `e^{O(d)}`. The scale is absolute (`d^{m'} eps_0`), not relative
   to `||Q||_2`, and the distribution's *shape* is irrelevant: **"AC" is an `L^2` size statement about one vector, not an
   anti-concentration statement.** `Psi_A` is an explicit sum-of-squares form of degree `2m'` with zero set `TC_V`. The
   exceptional set is exactly the per-pair braid neighbourhood `R^d \ G` (radii `e^{-eps_xy/2}`), which contains every
   collision stratum with radius `(1/2) e^{-min_coll eps/2}`, never smaller than the averaged radii of nearstrata L3.3.
2. **Modulo the slab ideal (Prop 3.1, Cor 3.2, Prop 3.3; proved).** `Psi_A(u) = q_A(l_u^{m'} mod J(A)_{m'})` for a
   positive definite form `q_A` on the `h_A(m')`-dimensional space `R[k']_{m'}/J(A)_{m'}`; slab polynomials enter only
   through `J(A)`. The target splits (lossily) into conditioning of `q_A` and an algebraic distance to `J(A)_{m'}`, and is
   equivalently a **restriction ratio** statement: the top component must keep an `e^{-(1/2 - o(1))W}` fraction of its
   cube `L^2` mass on the slab event, pointwise on `G_1`.
3. **Passage `H = infinity -> V_H` (Theorem 4.1; proved, Bombieri-Vaaler quoted).**
   `codim(V_H in V(A,m)) <= r_0 ln U/(2 ln H)` with `r_0 = dim Poly_{<m}|_A = e^{o(d)}` at `m'`: relative codimension
   `e^{-Omega(d)}`, improving theory.md Lemma 1.1 (which carries a `U ln 2/ln H` term). **No codimension bound gives a
   pointwise passage**; the two exact forms needed, (P1) span+norm and (P2) pointwise, are open in the regime and hold
   automatically wherever product certificates exist. Since fakes can use `TC_H`, the passage is part of the barrier.
4. **Proved cases (§5).** Uncut cube: exact formula for `Psi` (exact check), MIN^pair for `m' <= theta d`, `theta < 1`.
   Products/disjoint cuts: an exact product formula for `Psi` (exact check), giving quantitative additivity of the
   tangent degree. One cut (every `L >= 0`) and few cuts (`q <= (1/2) log_2(d/ln 4r) - 5.6` sign columns, or any profile with
   `<= 2^q` patterns): MIN^pair itself (with `kappa_H`) with explicit room, via pair-zonotopes in twin classes and a
   light-matching lemma. Thin cuts (`mu(A^c) <= 3^{-3m}/4`): `rho_A >= 1/2` by hypercontractivity. **Regime: nothing**, and
   the boundary is sharp for these methods: pair certificates are absent for every actual cluster and every
   gap-admissible fake once `W > 8 ln C'` (Prop 5.7, same weight count as rigidity's `m' > 2w`), regime slabs are thick,
   and profile-independent subsets of `A` have interpolation degree `<= 4j < m'`.
5. **Literature (§7).** Prop 7.1 (exact): at regime parameters a degree-`(2j+1)` slab polynomial (`465 < m' = 1864` at
   `n = 2901`) vanishes on all but `10^{-13.2}` of the cube, with influences `-> 0`. Hence every anti-concentration or
   invariance bound in terms of degree, variance and influences (Carbery-Wright via MOO, Meka-Nguyen-Vu, Fox-Kwan-
   Sauermann, PTF regularity) is necessarily vacuous at the needed degree; inverse Littlewood-Offord theorems
   (Kwan-Sauermann, Jin-Kwan-Sauermann-Wang) have the right genre but bounded degree and irrelevant thresholds;
   Kane's theorem is Gaussian. None applies.
6. **Minimal open statement (§8).** `(R_H)`: `min_{G_1} Psi_{A,H'} > d^{2m'} eps_0^2` — equivalent to MIN^pair up to
   `e^{O(d)}`, split as `(R_inf)` (restriction ratio `> e^{-(1/2-o(1))W}` on `G_1`) plus the passage (P2). Open for every
   regime box; by fake Thm 2.1 it must distinguish actual from fake pattern boxes through `J(A)` and `V_H` alone.

**Fake-violated input (rubric b): unchanged** — rigidity at the actual point, entering through Step 0. Everything here
is a statement about the box. **Growth/uniformity (rubric c): no movement.** What this round adds is a precise
negative answer to "is MIN^pair a discrete anti-concentration problem amenable to the anti-concentration literature?":
it is an `L^2` restriction problem for top components at degree `Theta(sqrt(n ln n))`, where the discrete
anti-concentration literature is provably vacuous.

## 1. Setting and notation (T = 1, w, j, m', regime)

Conventions are those of `round12/README.md` §1 and `round13/README.md` §1. *Theorem / Proposition / Lemma*: proved
here with all quantifiers, or quoted with an exact source (marked **[quoted]**, and **[memory]** when quoted from
memory without having re-read the source). *Exact*: integer or exact-rational computation (python-flint `fmpq_mat`).
*Evidence*: anything else. Scripts: `reduction_checks/` (each `x.py` writes `x_out.txt`).

**Box data** (every statement below names its setting). `d = n - 1`; points `k' in {-1,0,1}^d` (the base point `n` has
`u_n = 0`, `k_n = -sum k'` unconstrained at `T = 1`); `w = sqrt(8 d ln 4r)`, `j = floor(w/2)`, `m' = ceil(4w)`;

```text
A = { k' in {-1,0,1}^d : |<k, beta_p>| <= w for every column p },   U = |A| >= 3^d/2,
regime:  d >= 316 ln 4r  (so m' <= m_max),   W > W_1 .
```

For a finite `A` in the ternary cube and an order `m >= 1`:
* `l_u(k) = <k', u>`, `g_u = l_u^m |_A`; `Poly_{<m}|_A` = restrictions of polynomials of degree `< m`;
  `r_0(A,m) = dim Poly_{<m}|_A` (the rank of the moment matrix; `r_0 <= E(m) = binom(m+d-1, d)`);
* `V = V(A,m) = (Poly_{<m}|_A)^perp` in `R^A`; `V_H` = Q-span of the integer `c in V` with `|c|_inf <= H`;
* `Phi_c(u) = sum_k c_k <k',u>^m`, `kappa_c = |Phi_c|/(|c|_1 d^m)`, `kappa_V = sup_{c in V} kappa_c`,
  `kappa_H = sup_{c in V_H} kappa_c`; `TC_V = {kappa_V = 0}`, `TC_H = {kappa_H = 0}`;
* `J(A)` = the leading-form ideal of the vanishing ideal `I(A)`; `h_A(m) = dim R[k']_m / J(A)_m`.

**Cluster data.** Pair excesses `eps_xy = d_xy - W/2 + 2 ln C' >= 0`, `sum_{x<y} eps_xy <= N_2 abar`,
`abar = W/(2d) + 2 ln C'`. The **gap region** and its sphere:

```text
G = { u in R^d : |u_x - u_y| >= e^{-eps_xy/2} for all x < y in [n] }   (u_n = 0),     G_1 = G n {|u|_inf = 1}.
```

**Step 0 level** `eps_0 = C d e^{-W/4}/(m'+1)`. **MIN^pair** (referee_nearstrata §4): `min_{u in G_1} kappa_H(u) > eps_0`
(empty minimum = `+infinity`). Proved (nearstrata Thm 4.1 + referee, t1rigour Thm 6.8): MIN^pair for the cluster's
own box at `m' = ceil(4w)` (any order in `[ceil(4w), m_max]` also works) contradicts the existence of the cluster; if it
holds for all regime data on `n >= n_1` points, then `n = O(loglog R)` on all norms.

**Scale.** For `|u|_inf = 1` the Chebyshev ceiling is `d^m kappa_V(u) <= 2^{1-m} d^m` (J2), and the target is
`d^{m'} eps_0`; the logarithmic room is `W/4 - (m'-1) ln 2 - ln(Cd/(m'+1)) = (3/4 - o(1)) n ln n` at `W ~ 3 n ln n`
(t1rigour §4). Every loss of size `e^{O(d)}` below is therefore negligible for the growth question (it shifts `W`
by `O(d)`), but it is stated explicitly each time.

## 2. The reduction: which anti-concentration statement suffices for MIN^pair

Throughout §2, `A` is any finite subset of `{-1,0,1}^d`, `m >= 1`, and `nu` ranges over probability measures on `R^d`.
For `u in R^d` and `delta > 0` define

```text
Cheb_A(u)    = min_{P in Poly_{<m}} max_{k in A} |l_u(k)^m - P(k)|                    (directional Chebyshev error)
Psi_nu(u)    = min_{P in Poly_{<m}} E_nu (l_u^m - P)^2,      Psi_A := Psi_{mu_A}  (mu_A uniform on A)
SB_nu(u;delta) = sup_{P in Poly_{<m}} nu( |l_u^m - P| <= delta )                       (small-ball / concentration)
```

The minima exist (finite-dimensional least squares / LP, coercive modulo `I(supp nu)`), and `Cheb`, `Psi` depend on
`P` only modulo the vanishing ideal of the support: **they are intrinsically defined on the quotient ring
`R[k']/I(A)`**, so the slab polynomials (which lie in `I(A)`) never have to be "handled": they are invisible.

**Theorem 2.1 (all anti-concentration forms are one statement up to `e^{O(d)}`).** For every finite `A` in the ternary
cube, every `m >= 1` and every `u in R^d`:
1. (duality, nearstrata Prop 2.1) `d^m kappa_V(u) = Cheb_A(u)`.
2. (sup versus `L^2`) For every probability `nu` with `supp nu in A`: `Psi_nu(u) <= Cheb_A(u)^2`. Conversely
   `Cheb_{supp nu}(u)^2 <= Psi_nu(u) / min_{k in supp nu} nu(k)`. In particular
   `Psi_A(u) <= Cheb_A(u)^2 <= U Psi_A(u) <= 3^d Psi_A(u)`.
3. (probability thresholds are irrelevant) For `nu` supported in `A`, `delta > 0`, `p in (0,1]`:
   `SB_nu(u;delta) <= 1 - p  ==>  Psi_nu(u) >= p delta^2  ==>  Cheb_A(u) >= sqrt(p) delta`; and
   `Cheb_A(u) > delta  <=>  SB_{mu_A}(u; delta) < 1`.
4. (cube form) With `mu` uniform on `{-1,0,1}^d`: `SB_mu(u;delta) < mu(A)  ==>  Cheb_A(u) > delta`. Also
   `Psi_A(u) = min_P E_mu[(l_u^m - P)^2 1_A] / mu(A)` (the uniform measure on `A` *is* the cube measure conditioned on
   the slab event; conditioning on a subfamily of slab statistics does not give a measure supported in `A`, Remark 2.4).
5. (sum of squares) `Psi_A(u) = (1/U) sum_i Phi_{c_i}(u)^2` for any `l^2(A)`-orthonormal basis `(c_i)` of `V(A,m)`. So
   `Psi_A` is a nonnegative form of degree `2m` whose real zero set is `TC_V(A,m)`.

*Proof.* (1) is nearstrata Prop 2.1. (2) For `Q = l_u^m - P`, `E_nu Q^2 <= max_A Q^2`; minimise over `P` (the same `P`
is admissible on both sides). Conversely `E_nu Q^2 >= min nu * sum_{supp nu} Q^2 >= min nu * max_{supp nu} Q^2`, then
minimise; `mu_A` has `min nu = 1/U` and `U <= 3^d`. (3) If `nu(|Q| <= delta) <= 1 - p` for every admissible `Q`, then
`E_nu Q^2 >= p delta^2` for every `Q`; and `Cheb_A >= Cheb_{supp nu} >= sqrt(Psi_nu)` by (2). For the equivalence:
`Cheb_A(u) <= delta` iff some `Q` (the minimiser) has `|Q| <= delta` on all of `A`, i.e. `mu_A(|Q| <= delta) = 1`; the sup
in `SB` is attained by compactness of `{P : max_A |l_u^m - P| <= delta}` modulo `I(A)` (if nonempty). (4) If every
admissible `Q` has `mu(|Q| <= delta) < mu(A)`, then `{|Q| <= delta}` does not contain `A`, so `max_A |Q| > delta`; this holds
for the minimiser. The identity is `E_mu[f 1_A] = mu(A) E_{mu_A} f`. (5) For `g in R^A`,
`min_{f in V^perp} |g - f|_2^2 = |P_V g|_2^2 = sum_i <c_i, g>^2`, and `<c_i, g_u> = Phi_{c_i}(u)`; divide by `U`. The zero set
is `{P_V g_u = 0} = {g_u in Poly_{<m}|_A} = TC_V`. QED

**Corollary 2.2 (the exact reduction at `H = infinity`).** Fix cluster data and `m = m'`. Consider

```text
(MIN^pair_inf)   min_{u in G_1} Cheb_A(u)  >  d^{m'} eps_0                       [ = MIN^pair with kappa_V for kappa_H ]
(AC_2)           min_{u in G_1} Psi_A(u)   >  d^{2m'} eps_0^2                    [ L^2 form, uniform measure on A ]
(AC_cube)        for all u in G_1:  SB_mu(u; d^{m'} eps_0) < mu(A)                [ steelman's AC, threshold mu(A) >= 1/2 ]
(AC_p)           for all u in G_1:  SB_{mu_A}(u; d^{m'} eps_0 / sqrt p) <= 1 - p   [ any fixed p in (0,1] ]
```

Then `AC_2 ==> MIN^pair_inf`, `AC_cube ==> MIN^pair_inf`, `AC_p ==> AC_2 ==> MIN^pair_inf`, and conversely
`MIN^pair_inf ==> AC_2` **at the level `eps_0 3^{-d/2}`**. Since `ln(3^{d/2}) = 0.55 d` and the room is
`(3/4 - o(1)) n ln n`, the four statements are interchangeable for the growth question (each at level `eps_0` implies
the others at level `eps_0 e^{-O(d)}`, i.e. at `W` replaced by `W - 4.4 d`). *Proof.* Theorem 2.1 (2)-(4) with
`delta = d^{m'} eps_0`. QED

**What this settles about the task's questions.**
* *Pointwise in `u`.* Required, for every `u in G_1`: every point of `G_1` is compatible with the pair data that the
  argument uses (referee_nearstrata §4, J7), so a statement for "most `u`" in any measure does not suffice.
* *Which measure.* Any measure supported in `A`; the uniform one is optimal up to `e^{O(d)}`; the cube measure enters
  only through the threshold `mu(A)` (4). A measure not supported in `A` (Gaussian, cube, product tilts) gives nothing
  without the threshold `mu(A)`: on `TC` points the minimiser vanishes on `A`, which has mass `>= 1/2`.
* *Which scale, relative to which norm.* The scale is **absolute**: `delta = d^{m'} eps_0` for `Q` normalised by its
  leading form (`|u|_inf = 1`). It is **not** relative to `||Q||_2`: Carbery-Wright-type bounds
  `nu(|Q| <= eps ||Q||_2) <= C m eps^{1/m}` are about the *shape* of the distribution of `Q`, while for the minimiser
  `Q* = l_u^m - Proj l_u^m` one always has `max_A |Q*| >= ||Q*||_{L^2(mu_A)} = Psi_A(u)^{1/2}` (shape is irrelevant).
  **So the entire content of the "anti-concentration" statement is the size of `Psi_A(u)`, the `L^2` norm of the top
  orthogonal component of `l_u^m` on `A`; no small-ball information beyond that is needed or useful.**
* *Which exceptional set.* Exactly `R^d \ G`, a union of per-pair slabs `{|u_x - u_y| < e^{-eps_xy/2}}`. For every
  subspace `L` with a collision pair, `u in G` gives `dist_inf(u, L) >= rho_L := (1/2) exp(-(1/2) min_{(x,y) in coll(L)} eps_xy)`
  (`|u_x - u_y| <= 2 dist(u,L)` for a collision pair). Cheap spans (in the braid arrangement at `m' = ceil(4w)`,
  t1rigour Cor 6.4), projection-cheap spans `R^{S'}` and every collision stratum are therefore excluded with their
  **per-pair** radii, which are never smaller than the averaged radii `(1/2) e^{-N_2 abar/(2P(L))}` of nearstrata Lemma 3.3.
  Neighbourhoods of collision-free subspaces of `TC` are not excluded (Prop 5.1(i) of nearstrata: they kill MIN).

**Remark 2.3 (the Gaussian proxy is an `L^2` identity, not anti-concentration).** For `nu = N(0, I_d)`,
`Psi_nu(u) = m! |u|_2^{2m}` (Hermite chaos orthogonality: `l_u^m` has top chaos `|u|_2^m He_m(l_u/|u|_2)`). So the
"continuous proxy is easy" of steelman S3 needs no Carbery-Wright at all, and its discrete analogue on the uncut cube is
the exact formula of Theorem 5.1 below. The discrete difficulty is the passage from the cube to `A` (§3, §7), which is a
*restriction* problem for the top component, not an anti-concentration problem.

**Remark 2.4 (conditioning on slab statistics).** Write `S(k) = (<k,beta_p>)_p`. Then `mu_A` is the mixture of the
fibre measures `mu(. | S = t)` over the admissible values `t` (all `|t_p| <= w`), and every probability measure supported
in `A` is a mixture of measures supported on these fibres. Conditioning on a *subfamily* `P'` of the statistics gives
fibres `{S_p = t_p, p in P'}` that are **not** contained in `A` (the other slabs still cut), so local limit theorems for
conditioning on a few linear statistics can only reduce an `r`-cut problem to an `(r - |P'|)`-cut problem on a fibre.
When the columns together with `1` span `R^n` (e.g. Hadamard/Sylvester and Paley profiles; generic profiles with
`r >= d`), the full statistic is injective on sum-zero vectors (`S(k) = S(k'')` forces `k - k''` orthogonal to `R^n`), so
the fibres are single points. So "the cube with the
conditional measure" is literally `mu_A`, and nothing is gained by phrasing (AC) through conditioning.

## 3. Formulation modulo the slab ideal I(A): the graded form and the restriction ratio

### 3.1 The top component lives on `R[k']_m / J(A)_m`

Let `H_m(A) = Poly_{<=m}|_A  (-)  Poly_{<m}|_A` (orthogonal complement in `L^2(mu_A)`), and for a form `phi` of degree `m`
let `gr(phi) in H_m(A)` be the orthogonal projection of `phi|_A` onto `H_m(A)` (its *top component*).

**Proposition 3.1 (graded form).** `gr : R[k']_m -> H_m(A)` is surjective with kernel exactly `J(A)_m`. Hence
`q_A(phi) := ||gr(phi)||^2_{L^2(mu_A)}` is a positive definite quadratic form on the `h_A(m)`-dimensional space
`R[k']_m / J(A)_m`, and

```text
Psi_A(u) = q_A( l_u^m  mod J(A)_m ).
```

*Proof.* `Poly_{<=m} = Poly_{<m} + R[k']_m`, so `gr` is onto. `gr(phi) = 0` iff `phi|_A in Poly_{<m}|_A` iff `phi - P in I(A)`
for some `deg P < m` iff `phi = 0` or `phi` is the leading form of an element of `I(A)` of degree exactly `m`, i.e.
`phi in J(A)_m`. `Psi_A(u)` is the squared distance of `l_u^m|_A` to `Poly_{<m}|_A`, i.e. `||gr(l_u^m)||^2`. QED

*Normal form.* Fix a degree-compatible monomial order; let `B_m` be the standard monomials of degree `m` for `J(A)`
(for graded orders `in(I(A))` and `in(J(A))` agree). Each `phi` has a unique normal form `phi = sum_{b in B_m} nf_b(phi) k^b
mod J(A)_m`, and with `Gamma_m(A)` = the Gram matrix of `{gr(k^b)}_{b in B_m}` in `L^2(mu_A)` (a Schur complement of the
moment Gram matrix of standard monomials of degree `<= m`),

```text
Psi_A(u) = nf(u)^T Gamma_m(A) nf(u),     nf(u) = ( nf_b(l_u^m) )_{b in B_m}  (polynomial in u, degree m).
```

This is the precise sense in which "P is only defined up to `I(A)`" is handled: the slab polynomials (for a sign column
with side `S_p`, `prod_{|t|<=j}(k(S_p) - t)`, leading form `sigma_{S_p}^{2j+1}`; for a general integer column the product over
the `h_A(a_p) + 1` values of `<k, a_p>`), the ternary relations `k_x^3 - k_x`, and every other element of `I(A)` enter
only through `J(A)_m`, which the normal form quotients out. Along cheap spans `l_u^m` lies in the ideal generated by the width powers
`<k',v>^{h_A(v)+1}` (cheap spans are the linear base cycles of the width system, refJ J6(a)); the open question is what
else lies in `J(A)_m` (rounds 12-13: `J(A)` is strictly larger than the width ideal, e.g. in `M8`).

**Corollary 3.2 (split into conditioning and algebra; sufficient only).** With `lambda_min = lambda_min(Gamma_m(A))`,

```text
Psi_A(u) >= lambda_min * |nf(u)|_2^2 .
```

So `AC_2` follows from **(Cond)** `lambda_min(Gamma_{m'}(A)) >= e^{-a}` and **(Alg)** `|nf(u)|_2 >= e^{-b}` on `G_1`, whenever
`a + 2b < W/2 - 2 ln(C d^{m'+1}/(m'+1))`. (Alg) concerns the ideal `J(A)` and the Veronese map only; (Cond) concerns the
measure only. *Caveats.* (i) The split is basis dependent and can be very lossy. (ii) (Alg) is a Lojasiewicz
inequality for the distance from the Veronese image of `G_1` to a linear subspace; nearstrata Prop 2.7 applies to it
verbatim (bounds depending only on degree, height and qualitative zero set cannot reach the needed exponent).
(iii) The only general lower bound for (Cond) I can prove is the integrality bound
`det(Gamma^Sigma_m) = det G_{<=m}/det G_{<m} >= 1/det G_{<m}` (counting measure, integer Gram matrices), which gives
`ln(1/lambda_min) <= r_0 ln(U r_0) + O(h_A(m) ln U)` with `r_0 = e^{Theta(w ln d)}`: useless at the needed scale `n ln n`.
In the uncut cube `Gamma_m` is diagonal with `lambda_min = (2/9)^{floor(m/2)}(2/3)^{m mod 2}` (Theorem 5.1), so (Cond)
is harmless there.

### 3.2 The restriction ratio (the cube with the slab event)

**Definition.** For `u` with `Psi_cube(u) > 0` let

```text
rho_A(u) = min_P E_mu[(l_u^m - P)^2 1_A]  /  min_P E_mu[(l_u^m - P)^2]   in [0, 1]          (mu uniform on {-1,0,1}^d).
```

(`rho_A <= 1`: take for `P` the cube minimiser. `rho_A = mu(A) Psi_A / Psi_cube`.)

**Proposition 3.3 (dimensionless form of the target).** Let `m <= 2d`. For `u in G_1`, `MIN^pair_inf` at `u` follows from

```text
rho_A(u)  >  eps_0^2 * d^{2m} / Psi_cube(u),
```

and by Theorem 5.1(b) below `d^{2m}/Psi_cube(u) <= exp( 2m ln(ed/m) + 1.46 m + m abar_n )` on `G_1`, where
`abar_n = (1/d) sum_{x<n} eps_xn` is the average excess at the base point (`<= abar` for a suitable base point). At `m = m' = O(sqrt(n ln n))` and `abar = O(ln n)` this factor is `e^{O(m' ln n)} = e^{o(n)}`.
So the target is: **the top component of `l_u^m` keeps at least a fraction `e^{-(1/2 - o(1)) W}` of its cube `L^2`
mass on the slab event `A`**, for every `u` in the gap region. *Proof.* `Psi_A = rho_A Psi_cube/mu(A) >= rho_A Psi_cube`
and Corollary 2.2 (`AC_2 ==> MIN^pair_inf`). QED

`rho_A(u) = 0` exactly on `TC_V(A,m)` (where `Psi_cube > 0`). A restriction ratio bounded below by a constant is false
in the regime (cheap directions give `rho = 0`); the requirement `rho >= e^{-(1/2-o(1))W} = e^{-(3/2 - o(1)) n ln n}` is
enormously weaker than any "restriction theorem" one would normally prove, but it must hold pointwise on `G_1`.

## 4. Passage from H = infinity to V_H

MIN^pair is stated for `kappa_H`, and Step 0 holds only on `V_H` (integrality of `c` is essential: t1rigour Thm 3.1).
`kappa_H <= kappa_V`, so every statement of §2-§3 controls `kappa_V` and needs a passage. In `L^2` form:

```text
Psi_{A,H}(u) := min_{f in V_H^perp} E_{mu_A}(l_u^m - f)^2  =  Psi_A(u) - (1/U) |P_F g_u|_2^2,      F = V (-) V_H,
```

and `kappa_H(u) d^m >= Psi_{A,H}(u)^{1/2}` (Theorem 2.1(2) with `V_H` in place of `V`; same proof).

**Theorem 4.1 (codimension of `V_H`, via Bombieri-Vaaler).** For every finite `A` in the ternary cube, every `m`, and
every `H >= 2`:

```text
codim( V_H  in  V(A,m) )  <=  r_0 ln U / (2 ln H),        r_0 = dim Poly_{<m}|_A .
```

The same bound holds for `codim( Phi(V_H) in Phi(V) )`, where `Phi(V) = {Phi_c : c in V} = J(A)_m^perp` (apolar) has
dimension `h_A(m)`; `ker(Phi|_V) = V(A, m+1)`.

*Proof.* Let `M_0` be `r_0` linearly independent rows of the moment matrix (rows `(k^alpha)_{k in A}`, `|alpha| < m`, entries in
`{-1,0,1}`); `V = ker M_0`, of dimension `D = U - r_0`. **[quoted, memory]** Bombieri-Vaaler (*On Siegel's lemma*, Invent.
Math. 73 (1983), Thm 2): an `r_0 x U` integer matrix of rank `r_0` has `D` linearly independent integer kernel vectors
`x_1, ..., x_D` with `prod_i |x_i|_inf <= sqrt(det(M_0 M_0^T)) / g`, `g >= 1` the gcd of the maximal minors. (Its proof is
Minkowski's second theorem plus Vaaler's cube-slicing inequality — every central section of `[-1,1]^U` has volume
`>= 2^D` — both confirmed by the search, see §10.) Hadamard: `det(M_0 M_0^T) <= prod_rows |row|_2^2 <= U^{r_0}`. Each
`|x_i|_inf >= 1` (nonzero integer vectors), so at most `ln(U^{r_0/2})/ln H` of them exceed `H`; the others lie in `V_H`
and are independent. For the second claim, `Phi` is linear and `dim Phi(V_H) >= dim V_H - dim ker`. `c in V` has
`Phi_c = 0` iff `c` annihilates all degree-`m` forms, i.e. `c in V(A, m+1)`; `Phi(V) = J(A)_m^perp` is compute R2 / refJ J6(a). QED

**Comparison.** `theory.md` Lemma 1.1 gives `(U ln 2 + E_S ln(UB))/ln(2H+1)`; Theorem 4.1 removes the term `U ln 2/ln H`,
which dominates (`U >= 3^d/2` against `r_0 <= E(m')`). At `m' = ceil(4w)`, `r_0 <= binom(m'+d-1, d) = e^{O(w ln d)} = e^{o(d)}`, so
the relative codimension is `e^{-Omega(d)}`: at the record's `H = U^2` it is `<= r_0/4`, and at the largest height allowed by
rigidity, `ln H' = m'W/40` (t1rigour Remark 3.3), it is `<= 20 r_0 d ln 3/(m' W)`.

**What the codimension does not give.** `Psi_{A,H}(u) = Psi_A(u) - |P_F g_u|^2/U` with `dim F <= e := r_0 ln U/(2 ln H)`.
The vectors `P_V g_u` range over a space of dimension `h_A(m)` (Theorem 2.1(5), Prop 3.1); `e < h_A(m)` does not prevent a
given `u` from having `P_V g_u` almost inside `F`. **No dimension count yields a pointwise passage.** Two exact
reformulations of what is needed:

* **(P1) span form.** `Phi(V_H) = Phi(V)`, i.e. `V_H + V(A, m+1) = V(A,m)`, together with a norm comparison: a constant
  `K` such that for every `c in V` some `c' in V_H` has `Phi_{c'} = Phi_c` and `|c'|_1 <= K |c|_1`. Then
  `kappa_H >= kappa_V/K` everywhere. (Proof: `kappa_{c'}(u) = |Phi_c(u)|/(|c'|_1 d^m) >= kappa_c(u)/K`.)
* **(P2) pointwise form.** For every `u in G_1`: `Psi_{A,H'}(u) >= e^{-o(n ln n)} Psi_A(u)` with `ln H' = m'W/40`.

(P1) holds in the uncut cube, in one-cut boxes and wherever the certificates of §5 apply: those certificates are
products of differences with `|c|_1 <= 2^m <= H` (nearstrata Lemma 2.3), hence in `V_H`. **Neither (P1) nor (P2) is
known for any regime box.** The reason it matters: fake Thm 2.1 produces gap-admissible fakes from any nondegenerate
real point of `TC_H`, and `TC_H` can a priori be larger than `TC_V`; so the passage is part of the barrier question,
not a technicality (skeptic §2.3 item 3, made precise).

**Remark 4.2 (why small-support integer certificates do not help).** Circuits of the degree-`<m` evaluation matroid of
`A` have Hadamard heights `<= (s-1)^{(s-1)/2}` (`s` = circuit size; entries are `(s-1)`-minors of a `{-1,0,1}` matrix) and would lie in `V_{H'}` for `s <~ 0.15 n`. But every
circuit has `s >= 3^{m/2}`: a set `S` in the ternary cube with interpolation degree `>= m` has a standard monomial
`k^alpha` (`alpha_x <= 2`) of degree `>= m` whose divisors are all standard, so `|S| >= prod(alpha_x + 1) >= 3^{|alpha|/2}`.
At `m' = ceil(4w)`, `3^{m'/2} >> n`, so no circuit has a small height bound from Hadamard alone.

## 5. Proved cases of AC

Setting for all of §5: `T = 1`, any `d`, any order `m` as stated, probability normalisation for `Psi` unless "counting"
is said. Each theorem gives `AC_2` (hence `MIN^pair_inf`) and, where a product certificate is used, `MIN^pair` itself
with `kappa_H` (`H >= 2^m`), i.e. the passage (P1) is automatic there.

### 5.1 The uncut cube

**Theorem 5.1 (uncut ternary cube `B_1 = {-1,0,1}^d`, any `1 <= m <= 2d`).**
* (a) **Exact formula.**
  `Psi_cube(u) = sum_{alpha in {0,1,2}^d, |alpha| = m} (m!/alpha!)^2 (2/3)^{n_1(alpha)} (2/9)^{n_2(alpha)} u^{2 alpha}`,
  `n_i(alpha) = #{x : alpha_x = i}`. In particular `TC = {u : |supp u| < m/2}` (nearstrata Thm 2.4), and
  `Gamma_m(B_1)` is diagonal in the standard monomials `k^alpha` (`alpha_x <= 2`) with
  `lambda_min = (2/9)^{floor(m/2)} (2/3)^{m - 2 floor(m/2)}`.
* (b) **Lower bound on the gap region.** For `m <= 2d` even and `u in G_1`,
  `Psi_cube(u) >= (m!)^2 2^{-m} (2/9)^{m/2} exp(-m abar_n)`, `abar_n = (1/d) sum_{x<n} eps_xn`; hence
  `d^{2m}/Psi_cube(u) <= exp(2m ln(ed/m) + 1.46 m + m abar_n)`. (Odd `m`: same with one `alpha_x = 1`, factor `3/2`.)
* (c) **MIN^pair for the cube box.** For `H >= 2^m`, `u in G_1` and `m <= 2d` even,
  `kappa_H(u) >= (m!/(2d)^m) exp(-(m/2) abar_n)`. So MIN^pair holds whenever
  `(d/2)(abar - 2 ln C') - (m/2) abar_n > m ln(2ed/m) + ln(Cd/(m+1))`; with the base point chosen so that `abar_n <= abar`
  this holds for `m <= theta d`, `theta < 1`, as soon as
  `abar >= (2/(1-theta)) (theta ln(2e/theta) + ln C') + O(ln d / d)`.

*Proof.* (a) Under the uniform measure on `{-1,0,1}` the polynomials `p_0 = 1`, `p_1 = k`, `p_2 = k^2 - 2/3` are orthogonal
with `||p_1||^2 = 2/3`, `||p_2||^2 = E k^4 - (4/3) E k^2 + 4/9 = 2/9`. The products `p_alpha = prod_x p_{alpha_x}(k_x)`,
`alpha in {0,1,2}^d`, are an orthogonal basis of `L^2(B_1)`, and `Poly_{<=D}|_{B_1} = span{p_alpha : |alpha| <= D}`
(triangularity: `k^alpha = p_alpha + ` lower for `alpha <= 2`, and `k_x^a = k_x^{a-2}` on `B_1` for `a >= 3`). Expand
`l_u^m = sum_{|beta| = m} (m!/beta!) u^beta k^beta`: terms with some `beta_x >= 3` restrict to degree `< m`, and the others
have top component `p_beta`. So `gr(l_u^m) = sum_{beta <= 2, |beta| = m} (m!/beta!) u^beta p_beta`, whose squared norm is the
formula. The zero set: all terms are nonnegative and some term is positive iff `u` has `>= m/2` nonzero coordinates.
(b) Keep the single term `alpha = 2 * 1_X`, `X` = the `m/2` coordinates `x < n` with smallest `eps_xn`: its value is
`(m!/2^{m/2})^2 (2/9)^{m/2} prod_{x in X} u_x^4`, and `|u_x| = |u_x - u_n| >= e^{-eps_xn/2}` on `G`, while
`sum_{x in X} eps_xn <= (m/(2d)) sum_x eps_xn = (m/2) abar_n`. The second claim uses `m! >= (m/e)^m` and
`ln 2 + (1/2) ln 4.5 < 1.46`. (c) The zonotope `-sum_{x in X} e_x + sum_{x in X} [0,2] e_x` lies in `B_1`; nearstrata
Lemma 2.3 with `alpha_x = 2` gives `kappa_H(u) >= (m!/(2d)^m) prod_X u_x^2`. Compare with `eps_0 = (Cd/(m+1)) e^{-W/4}` and
`W/4 = (d/2)(abar - 2 ln C')`. Averaging `abar_n` over the `n` choices of base point gives `(2/(n d)) sum_{x<y} eps_xy <= abar`, so some base point has
`abar_n <= abar` (all statements of t1rigour hold for every base point).
QED (Exact check of (a): `check_cube_product_out.txt`, `d = 3, 4, 5`, every `m <= 2d+1`, random rational `u`.)

### 5.2 Products and disjoint cuts

**Theorem 5.2 (product formula).** Let `A = X x Y` (coordinates split as `u = (u_X, u_Y)`) with the product of the
uniform measures. Then for every `m`,

```text
Psi_A^{(m)}(u) = sum_{i=0}^m binom(m,i)^2 Psi_X^{(i)}(u_X) Psi_Y^{(m-i)}(u_Y)        (Psi^{(0)} := 1).
```

Consequently `nu_A(u) = nu_X(u_X) + nu_Y(u_Y)` for the tangent degree `nu = max{i : Psi^{(i)} > 0}` (zetc §1, additivity
on products), now with a quantitative form, and `AC_2` for a product follows from `AC_2`-type lower bounds on the factors
at complementary orders.

*Proof.* With `H_a(X) = Poly_{<=a}|_X (-) Poly_{<a}|_X`, `L^2(X x Y) supseteq Poly_{<=m}|_{XxY} = (+)_{a+b<=m} H_a(X) (x) H_b(Y)`
orthogonally (polynomials in separate variables; tensor products of orthogonal decompositions). Expand
`l_u^m = sum_i binom(m,i) l_{u_X}^i l_{u_Y}^{m-i}`; the `H_a (x) H_b` component of `l_{u_X}^i l_{u_Y}^{m-i}` vanishes unless
`a <= i`, `b <= m - i`, so its component in total degree exactly `m` is `gr_i(l_{u_X}^i) (x) gr_{m-i}(l_{u_Y}^{m-i})`. These
are orthogonal for different `i`; take norms. `Psi_A^{(m)} > 0` iff some `i <= nu_X` has `m - i <= nu_Y` iff
`nu_X + nu_Y >= m` (and `Psi^{(i)} > 0` iff `i <= nu`, since `l^i in J` implies `l^{i+1} in J`). QED
(Exact check: `check_cube_product_out.txt`, `X` = one-cut set `s = 3, L = 1`, `Y` = cube `d = 2`, `m <= 6`.)

Disjoint signed cuts give products of one-cut blocks and a cube block (at `T = 1` the base coordinate is free), so
Theorems 5.1-5.3 combine by Theorem 5.2. These are the boxes for which refJ J6(b) proves `Loc(1)` via the strong
Lefschetz property.

### 5.3 Monotonicity, pair certificates, one cut, few cuts

**Lemma 5.3 (monotonicity).** If `B` is contained in `A` then `Cheb_A >= Cheb_B`, `kappa_{V(A,m)} >= kappa_{V(B,m)}`,
`TC_V(A,m) in TC_V(B,m)`, and (counting measure) `Psi^Sigma_A >= Psi^Sigma_B`. *Proof.* `V(B,m)` extended by zero lies in
`V(A,m)`; a max over `A` is at least the max over `B`; a sum over `A` of squares is at least the sum over `B`. QED
(Exact check of the counting form on 489 one-cut instances, `check_onecut_zero_out.txt`.)

So AC may be proved on any well-understood subset `B` of `A` whose tangent cone avoids `G`. For pairs `(x, y)` of points
of `[n]` write `Delta_p(x,y) = beta_{p,x} - beta_{p,y}` and `v_xy = e_x - e_y` (in `k`-space; when `y = n` this is `e_x` in
`k'`-coordinates).

**Theorem 5.4 (pair-zonotope certificates).** Let `M = {(x_i, y_i)}_{i <= s}` be pairwise disjoint pairs with
`sum_{i<=s} |Delta_p(x_i,y_i)| <= w` for every column `p`. Then `Z_M = {sum_i t_i v_{x_i y_i} : t_i in {-1,0,1}}` lies in `A`, and
for `m <= 2s` even, `H >= 2^m`, every `M' in M` with `|M'| = m/2` and every `u`:

```text
kappa_H(u)  >=  (m!/(2d)^m) prod_{(x,y) in M'} (u_x - u_y)^2   >=   (m!/(2d)^m) exp( - sum_{(x,y) in M'} eps_xy )   on G.
```

*Proof.* Points of `Z_M` have ternary coordinates (disjoint pairs) and `<k, beta_p> = sum_i t_i Delta_p(x_i,y_i)`, of modulus
`<= w`. `Z_M = k_0 + sum_i [0,2] v_i` with `k_0 = -sum_i v_i`; nearstrata Lemma 2.3 with `alpha_i = 2` on `M'` gives the first
bound (`<v_xy, u> = u_x - u_y`), and `G` gives `(u_x - u_y)^2 >= e^{-eps_xy}`. QED

**Lemma 5.5 (light matchings inside classes).** Let `C_1, ..., C_q` be disjoint subsets of `[n]` with `|C_a| >= c_0 >= 2`,
`s_* = sum_a floor(|C_a|/2) >= m/2`. Then there is a set `M'` of `m/2` disjoint pairs, each inside some `C_a`, with
`sum_{M'} eps_xy <= (m/2) N_2 abar / ((c_0 - 1) s_*)`. *Proof.* Take independent uniformly random maximum matchings of
the `C_a`; a pair inside `C_a` is used with probability `<= 1/(|C_a| - 1)`, so the expected total excess of the `s_*` matched
pairs is `<= sum_a sum_{pairs in C_a} eps_xy/(c_0 - 1) <= N_2 abar/(c_0 - 1)`. Fix an outcome at most the mean and keep its
`m/2` lightest pairs. QED

**Corollary 5.6 (one cut, every `L >= 0`; few cuts).** Call `C in [n]` a *twin class* if `Delta_p(x,y) = 0` for all
`x, y in C` and all `p`.
* (a) *One cut* `A_{S,L}` (one sign column, side `S`, `|S| = s`, `L >= 0`): `S` and `[n] \ S` are twin classes. With
  `c_0 = min(s, n - s)`, for `u in G_1` and even `m <= n - 2`,
  `kappa_H(u) >= (m!/(2d)^m) exp(-(m/2) abar n(n-1)/((c_0 - 1)(n-2)))`. So MIN^pair holds for one-cut boxes as soon as
  `(d/2)(abar - 2 ln C') - (m/2) abar n(n-1)/((c_0-1)(n-2)) > m ln(2ed/m) + ln(Cd/(m+1))`; for a balanced cut this needs
  only `m <= (1 - eps) n/2` and `abar >= C(eps)(1 + ln C')`. **This holds for every `L >= 0`**, i.e. it uses only the slice
  `{k(S) = 0}` (a superset of `TC(A_{S,L})` that is still inside the braid arrangement), not the exact one-cut tangent cone.
* (b) *Few cuts.* If the profile has `q` sign columns (or any columns with at most `2^q` distinct patterns), then with
  `c_0 = n/2^{q+1}`: `s_* >= (n/2 - 2^q)/2`, and for `u in G_1`,
  `kappa_H(u) >= (m!/(2d)^m) exp(-(m/2) abar 2^{q+2} (1 + O(2^q/n)))`. So MIN^pair holds when
  `m 2^{q+2} (1 + o(1)) < d` and `abar >~ ln(d/m) + ln C'`: e.g. at `m' = ceil(4w)`, for
  `q <= (1/2) log_2(d/ln 4r) - 5.6`.

*Proof.* (a) Pairs inside a twin class have `Delta = 0`, so every matching inside twin classes satisfies the hypothesis of
Theorem 5.4. Lemma 5.5 with `s_* >= (n-2)/2`, then compare with `eps_0` as in Theorem 5.1(c). (b) Points in classes of size
`< c_0` number `< 2^q c_0 = n/2`; Lemma 5.5 on the other classes, `N_2/((c_0 - 1) s_*) <= 2^{q+2}(1 + O(2^q/n))`. The numerical
form: `m' <= 4w + 1` and `4 * 2^{q+2} w < d` iff `2^q < d/(16 w) = sqrt(d/(8 ln 4r))/16`. QED

**Proposition 5.7 (pair certificates are absent for every cluster and every gap-admissible fake).** Let cluster data
(actual, or a gap-admissible angle model) have `W > 8 ln C'`. Then every family `M` of disjoint pairs as in Theorem 5.4
has `|M| <= w W/(W - 4 ln C') < 2w`, so `Z_M` has total length `< 4w <= m'`: Theorem 5.4 never applies at any admissible
order. *Proof.* `|a_{p,x} - a_{p,y}| = (e_p/2)|Delta_p(x,y)|`, so
`sum_{i} d_{x_i y_i} = sum_p (e_p log p/2) sum_i |Delta_p(x_i,y_i)| <= (w/2) sum_p e_p log p <= wW/2`, while each
`d_xy >= W/2 - 2 ln C'` (`eps_xy >= 0`, Lemma 3.1 of nearstrata; fakes satisfy the same pair inequalities). QED

This is the pair case of nearstrata Remark 2.10 and of fatpoint Thm 3.1 (same threshold `W > 8 ln C'`), with an
elementary proof. It is also the exact boundary of the provable cases: Corollary 5.6 works because twin classes exist,
and the regime (no twins, every pair separated by half the weight) is exactly where Proposition 5.7 removes them. The
slab structure that rigidity uses (`m' > 2w`) and the one that kills pair certificates (`|M| < 2w`) are the same weight
count.

### 5.4 Thin cuts (hypercontractivity)

**Theorem 5.8.** Let `A` be any subset of `{-1,0,1}^d` with `mu(A^c) <= eta`. Then for every `u` and `m`,
`rho_A(u) >= 1 - 3^{3m/2} eta^{1/2}`, i.e. `mu(A) Psi_A(u) >= (1 - 3^{3m/2} eta^{1/2}) Psi_cube(u)`.

*Proof.* For `Q = l_u^m - P` (`deg P < m`): `E_mu[Q^2 1_{A^c}] <= ||Q||_4^2 eta^{1/2}`. `Q` has (Efron-Stein) degree `<= m` on the
product space `({-1,0,1}, uniform)^d`, whose atoms have probability `1/3`. **[quoted, memory]** O'Donnell, *Analysis of
Boolean Functions* (2014), Thm 10.21: on a product of finite spaces with minimal atom probability `lambda`,
`||f||_q <= (sqrt(q-1) lambda^{1/q - 1/2})^k ||f||_2` for `deg f <= k`; with `q = 4`, `lambda = 1/3`: `||Q||_4 <= 3^{3m/4} ||Q||_2`.
Hence `E[Q^2 1_A] >= (1 - 3^{3m/2} eta^{1/2}) E Q^2 >= (1 - 3^{3m/2} eta^{1/2}) Psi_cube(u)`; minimise over `P`. QED

Thin cuts inherit the cube's AC (Theorem 5.1) with factor `1/2` once `eta <= 3^{-3m}/4`. For a balanced sign column the
slab `{|<k,beta>| <= w'}` has complement mass `e^{-Theta(w'^2/d)}` (Hoeffding above, the central limit theorem below, since
`Var <k,beta> = Theta(d)`), so thinness at order `m = 4w'` forces `w'^2 >~ c m d`, i.e. `w' = Omega(d)`; at the regime width
`w' = w = O(sqrt(d ln r))` each slab removes mass `e^{-O(ln r)}`, up to `1/(2r)`, against the `3^{-3m'}` allowed. **The
regime slabs are thick at every order that matters.**

### 5.5 What is proved in the regime

Nothing beyond the record. By Lemma 5.3 a proof may pass to any subset `B` of `A`, but:
* profile-independent subsets lie in `S_w = {|k|_1 <= 2j, sum k = 0}` (sup over columns in `[-1,1]^n` of `|<k,beta>|` is
  `|k|_1`), whose interpolation degree is `<= 4j < m'` (standard monomials have support `<= 2j`), so `Cheb_{S_w} = 0`;
* product/slice subsets need twin classes (Corollary 5.6), and pair and zonotope certificates are excluded
  (Proposition 5.7; nearstrata Remark 2.10; fatpoint Thm 3.1);
* thin-cut arguments need `mu(A^c) <= 3^{-3m}` (Theorem 5.8).
So every regime proof must use a non-product, profile-specific certificate or subset. None is known.

## 6. Comparison with exact tangent-cone results (zetc one-cut formula, refJ Lefschetz)

1. **One cut, exact zero set (exact check).** On `A_{S,L}`, `Psi_A(u) = 0` iff `nu_S(u) = 2a + min(2 sigma, 2L + 2t, s + L) < m`
   (zetc Thm 2.1 / nearstrata Thm 2.5, with `eta`-clusters replaced by exact value classes): verified exactly on 489
   planted instances (`d = 5, 6`; `(s, L) = (3,1), (4,1), (4,0), (5,2)`; every `m` up to `rho + 1`), no mismatch
   (`check_onecut_zero_out.txt`). This is the `L^2` form of the zetc formula; by Theorem 5.2 it extends to disjoint cuts
   (`nu` additive).
2. **Lefschetz.** refJ J6(b): `J(A_S(L)) = (k_x^3, sigma_S^{2L+1})` in all degrees. In the language of Prop 3.1 the form
   `q_A` lives on `R[k']_m/(k^3, sigma_S^{2L+1})_m`, of dimension `h_A(m) = max(h_m - h_{m-2L-1}, 0)` convolved with the
   cube part. The Hilbert-function values computed here agree with the record: e.g. for `M8` the computed
   `h_A(6), h_A(7), h_A(8) = 287, 175, 42` equal `binom(12,6) - 637`, `binom(13,6) - 1541`, `binom(14,6) - 2961`, the exact
   `dim J(A)_{6,7,8}` of refJ J6(c) (`cond_lambda_out.txt`, float rank; consistent with the exact record values).
3. **MIN^pair needs much less than the exact tangent cone.** Corollary 5.6(a) proves MIN^pair for every one-cut box
   using only the slice `{k(S) = 0}` (`L = 0`), whose tangent cone `{u : u - a 1_S` has small support`}` is strictly larger
   than `TC(A_{S,L})` for `L >= 1` but still lies in the braid arrangement with `>= binom(s - m/2, 2)` collision pairs per
   stratum. So the exact tangent cone (and the Lefschetz property) is not what MIN^pair uses: any subset whose tangent
   cone avoids `G` suffices (Lemma 5.3). This is the "superset of TC inside collision strata" principle, and it is what
   would have to be found in the regime.
4. **Conditioning (Cor 3.2 (Cond)), evidence only** (`cond_lambda_out.txt`, float64 SVD): `lambda_min(Gamma_m)` in the
   coefficient metric equals the exact cube value `(2/9)^{floor(m/2)}(2/3)^{m mod 2}` on the cube (`d = 6`, `m <= 8`); on
   one-cut boxes it is smaller by factors `4` to `30`; on `M8` (`Syl(8)`, `L = 2`) it is `4.0e-5`, `5.0e-6`, `6.4e-5`, `1.6e-3`
   at `m = 5..8` (cube: `3.3e-2`, `1.1e-2`, `7.3e-3`, `2.4e-3`). So in the only faithful toy the conditioning loss is a factor
   `<= 2.2e3` per order, i.e. `e^{O(m)}` at these sizes — harmless against the room `e^{(3/4) n ln n}` if it persisted
   (`e^{O(m')} = e^{O(sqrt(n ln n))}`), but nothing is proved about it and the toys are LP-subcritical.

## 7. Literature: what the anti-concentration toolbox gives and what it cannot give

Sources actually located this round are listed in §10; statements marked **[memory]** were not re-read.

**Proposition 7.1 (cube anti-concentration is false at the needed degree; exact numbers).** For `S in [d]`, `|S| = s`, and
an integer `j >= 0`, the polynomial `Q_S(k) = prod_{|t| <= j} (k(S) - t)` has degree `2j+1`, is not identically zero on the
cube, vanishes on every slab box whose column has side `S` (so on `A`), and

```text
Pr_mu( Q_S != 0 ) = Pr( |k(S)| > j ) <= 2 exp(-(j+1)^2/(2s)) .
```

At regime parameters (exact, `slab_concentration_out.txt`): `n = 2901`, `r = n`, `j = 232`, `deg Q_S = 465 < m' = 1864`:
`Pr(Q_S != 0) = 10^{-13.18}` (`s = d/2`), `10^{-6.91}` (`s = d`); `n = 6000`, `r = n^2`: `10^{-25.70}` (`s = d/2`). All
normalised influences are `<= (2j+1)/s -> 0` (symmetric in `S`, total influence `<= degree x variance`).
*Consequence.* Any inequality `sup_c Pr_mu(|Q - c| <= eps sigma(Q)) <= F(D, eps, tau)` valid for all degree-`D` polynomials
with maximal normalised influence `<= tau` on the ternary cube must have `F(2j+1, 0, (2j+1)/s) >= 1 - 2e^{-(j+1)^2/(2s)}`.
So every anti-concentration bound whose hypotheses are degree, variance and influences — Carbery-Wright transferred
by an invariance principle (Mossel-O'Donnell-Oleszkiewicz, error `O(D tau^{1/(8D)})` **[memory, the exponent varies
between sources]**), Meka-Nguyen-Vu (Theory of Computing 12 (2016) art. 11; bound decaying like a power `1/(4D+1)` of the
rank **[memory]**), Fox-Kwan-Sauermann (Math. Proc. Camb. Phil. Soc. 171 (2021); nonnegative coefficients), regularity
lemmas for low-degree PTFs — is **necessarily** vacuous at degree `>= 2j + 1 ~ w`, i.e. already at a quarter of `m'`.
This is not a weakness of the proofs: in the Gaussian proxy `Q_S` is a nonzero polynomial in one Gaussian and has
`Pr(Q_S = 0) = 0`. *Proof.* Hoeffding for `k(S)` (ranges `2`); the values are exact trinomial sums. QED

**Reading of the toolbox against the reduction (§2-§3).**
* *Carbery-Wright / Gaussian proxy.* Not needed even there: the Gaussian statement is the `L^2` identity
  `Psi_Gauss = m! |u|_2^{2m}` (Remark 2.3). Its discrete analogue is Theorem 5.1(a).
* *Invariance principles* (MOO; Filmus-Mossel for slices, PTRF 175 (2019), degree `o(n)` harmonic polynomials on Boolean
  slices). The slice result is the closest in spirit to Lemma 5.3/Corollary 5.6 (slices inside `A`), but every
  invariance statement compares distributions of a *fixed* polynomial under two measures, and Prop 7.1 shows the cube
  and `A` (or the cube and the Gaussian) disagree completely at degree `~ w`. What is needed is the opposite direction:
  a lower bound for the *minimum over the coset* `l_u^m + Poly_{<m}`, which no invariance principle addresses.
* *Kane's structure theorem* (Ann. Probab. 45 (2017) 1612-1679): Gaussian, decomposition into a bounded number of
  "diffuse" polynomials with bounds depending on the degree only. Continuous; for fixed degree; vacuous here.
* *Inverse Littlewood-Offord* (Nguyen-Vu linear; Kwan-Sauermann quadratic, Compositio 161 (2025) and arXiv 1909.02089;
  Jin-Kwan-Sauermann-Wang, arXiv 2505.23335, "the Meka-Nguyen-Vu bound can be substantially improved unless the
  polynomial is close to one with special algebraic structure"). These are the right *genre* — concentration forces
  algebraic structure — and in our setting the "algebraic structure" is exactly membership of `Q` in `I(A)` modulo lower
  degree, i.e. `u in TC`. But they are proved for bounded degree (or with losses exponential in the degree), with
  probability thresholds of order `n^{-1/2}`-to-`1`, while §2 shows that the probability threshold is irrelevant and the
  degree is `Theta(sqrt(n ln n))`.
* *Halasz-type bounds* are for linear forms with arithmetic conditions; degree 1.
* *Local limit theorems for conditioning on linear statistics*: Remark 2.4 (they reduce `r` cuts to `r - q` cuts).
* *Remez-type inequalities on discrete sets* (Klein-Slote-Volberg-Zhang, arXiv 2301.01438: dimension-free Remez on
  hypergrids, constants depending on the degree). The restriction ratio of §3.2 is a Remez inequality *modulo lower
  degree for pure powers*; full Remez on `A` is false at degree `>= 2j+1` by Prop 7.1.
* *Pluripotential theory* (steelman S4): the homogeneous Bernstein-Walsh bound is the continuous `Cheb`; on the finite set
  `A` its role is played by `q_A` on `R[k']_m/J(A)_m` (Prop 3.1), which is finite-dimensional linear algebra.

**Net.** The reduction of §2 removes every probabilistic feature (threshold, shape) and leaves an `L^2`
restriction statement for top components; the toolbox addresses shapes of distributions of fixed polynomials, and
its discrete members are provably vacuous at the needed degree (Prop 7.1). No tool in this list applies.

## 8. The minimal open statement

Setting: `T = 1`, `w = sqrt(8 d ln 4r)`, `m' = ceil(4w)`, regime `d >= 316 ln 4r`, `W > W_1`, `H' = e^{m'W/40}` (admissible in
rigidity once `W > 10 d ln 3/m' + 10 ln(edC'/m')`, t1rigour Remark 3.3), cluster data `(eps_xy)`, gap sphere `G_1`.

> **(R_H) — the open statement.** `min_{u in G_1} Psi_{A,H}(u) > d^{2m'} eps_0^2`, where
> `Psi_{A,H}(u) = min_{f in V_H^perp} E_{mu_A}(l_u^{m'} - f)^2` (Theorem 2.1 with `V_H` for `V`).

* **Equivalence.** `(R_H) ==> MIN^pair ==> (R_H)` at level `eps_0 3^{-d/2}` (Theorem 2.1(2) for `V_H`). Since `e^{O(d)}` is
  negligible, `(R_H)` *is* MIN^pair for the growth question; it is the minimal statement of anti-concentration type, and
  no probabilistic weakening exists (§2).
* **Splitting.** `(R_H) <== (R_inf) + (P2)`, where
  * **(R_inf)** (restriction ratio, `H = infinity`): `rho_A(u) > eps_0^2 d^{2m'}/Psi_cube(u)` for `u in G_1`; by Theorem
    5.1(b) it suffices that `rho_A(u) > exp(-W/2 + O(m'(ln n + abar_n)))`: **the top component of `l_u^{m'}` keeps an
    `e^{-(1/2 - o(1))W}` fraction of its cube `L^2` mass on the slab event.**
  * **(P2)** (passage): `Psi_{A,H'}(u) >= e^{-o(n ln n)} Psi_A(u)` on `G_1` (§4; codimension alone cannot give it).
* **Qualitative cores.** `(R_inf)` needs `TC_V(A, m') n G = empty`; `(R_H)` needs `TC_{H'}(A, m') n G = empty`.
* **Growth version.** For `M = o(log R/loglog R)` only clusters with `n >= c log R/loglog R` must be excluded
  (nearstrata §4.3): `(R_H)` is needed only at `abar <= (1/c + o(1)) ln n`, i.e. `W <= (2/c + o(1)) n ln n`, where the target
  `eps_0^2 = e^{-W/2 + O(ln n)}` is *larger*; the hardest case is the near-extremal `W ~ 3 n ln n`. No weaker AC gives a
  weaker growth improvement: the reduction is all-or-nothing in `n` (Cor 4.2 of nearstrata).
* **Fake test.** By fake Thm 2.1, `(R_H)` fails for the data of any (Mar) fake pattern box with a nondegenerate real
  point of `TC_{H'}` in `G` satisfying (2.2) (at the barrier scale). Since `(R_H)` is a statement about the box alone
  (arithmetic is spent in Step 0), it must distinguish actual regime boxes from fake pattern boxes through their
  algebra (`J(A)`, `V_H`) alone. Nothing in §5 does: all proved cases use twin classes or thin cuts, which neither
  actual nor fake regime data have (Prop 5.7, Thm 5.8).

**Status of each piece.**

| statement | status | where |
| --- | --- | --- |
| MIN^pair `<=>` `(R_H)` up to `e^{O(d)}`; probability thresholds irrelevant; scale absolute | **proved** | Thm 2.1, Cor 2.2 |
| `Psi_A` = SOS form; graded form on `R[k']_m/J(A)_m`; split (Cond)+(Alg) | **proved** | Thm 2.1(5), Prop 3.1, Cor 3.2 |
| `codim(V_H in V) <= r_0 ln U/(2 ln H)` (improves theory.md L1.1) | **proved** (Bombieri-Vaaler quoted) | Thm 4.1 |
| pointwise passage (P1)/(P2) | **open** in the regime; automatic for product certificates | §4 |
| cube: exact `Psi`, `MIN^pair` for `m <= theta d` | **proved**, exact check | Thm 5.1 |
| products / disjoint cuts: `Psi` product formula, quantitative `nu`-additivity | **proved**, exact check | Thm 5.2 |
| one cut (every `L >= 0`), few cuts (`q <= (1/2) log_2(d/ln 4r) - 5.6`): MIN^pair with explicit room | **proved** | Cor 5.6 |
| thin cuts `mu(A^c) <= 3^{-3m}/4` | **proved** (hypercontractivity quoted) | Thm 5.8 |
| pair certificates for actual/fake cluster data | **absent** (`W > 8 ln C'`) | Prop 5.7 |
| cube anti-concentration at degree `2j+1` | **false** (exact `1 - 10^{-13}` concentration at `n = 2901`) | Prop 7.1 |
| `(R_inf)`, `(R_H)`, MIN^pair for any regime box | **open** | — |
| growth improvement / uniform bound | **not proved** | — |

## 9. Computations (reduction_checks/)

| script | what | type |
| --- | --- | --- |
| `lib.py` | exact `Psi^Sigma_A(u)` over `Q` (python-flint): pivot monomials of degree `< m` on `A`, Gram inverse, residual | library |
| `check_cube_product.py` | Thm 5.1(a) formula vs exact least squares (`d = 3,4,5`, all `m <= 2d+1`); Thm 5.2 product formula (`X` one-cut `s=3, L=1`, `Y` cube `d=2`, `m <= 6`) | **exact**; all match |
| `check_onecut_zero.py` | `Psi_A(u) = 0 <=> nu_S(u) < m` on one-cut boxes (489 planted instances); counting monotonicity `Psi_A >= Psi_{slice}` | **exact**; 0 mismatches, 0 violations |
| `cond_lambda.py` | `h_A(m)` and `lambda_min(Gamma_m)` (coefficient metric) on cube `d=6`, one-cut boxes, `Syl(8)` `L=1`, `M8` | **evidence** (float64 SVD); `h_A` consistent with refJ's exact `dim J(A)_{6,7,8}` for `M8` |
| `slab_concentration.py` | `Pr(|k(S)| > j)` at regime parameters (`n = 2901, 5000` with `r = n`; `n = 6000` with `r = n^2`) | **exact** (trinomial big integers), printed in log10 |

All runs finished in under one minute each.

## 10. Sources consulted this round (WebSearch; abstracts/listing pages only)

* D. M. Kane, *A structure theorem for poorly anticoncentrated polynomials of Gaussians and applications to the study of
  polynomial threshold functions*, Ann. Probab. 45 (2017) 1612-1679; arXiv 1204.0543.
* R. Meka, O. Nguyen, V. Vu, *Anti-concentration for polynomials of independent random variables*, Theory of Computing
  12 (2016), art. 11.
* J. Fox, M. Kwan, L. Sauermann, *Combinatorial anti-concentration inequalities, with applications*, Math. Proc. Camb.
  Phil. Soc. 171 (2021) 227-248; arXiv 1905.12142.
* M. Kwan, L. Sauermann, *An algebraic inverse theorem for the quadratic Littlewood-Offord problem*, arXiv 1909.02089;
  *Resolution of the quadratic Littlewood-Offord problem*, Compositio Math. 161 (2025) 3089-3139.
* Z. Jin, M. Kwan, L. Sauermann, Y. Wang, *Algebraic aspects of the polynomial Littlewood-Offord problem*, arXiv 2505.23335.
* Y. Filmus, E. Mossel, *Harmonicity and invariance on slices of the Boolean cube*, PTRF 175 (2019) 721-782; arXiv 1507.02713.
* E. Mossel, R. O'Donnell, K. Oleszkiewicz, *Noise stability of functions with low influences*, Ann. Math. 171 (2010)
  295-341 (error exponent quoted from a secondary source; varies between sources).
* J. D. Vaaler, *A geometric inequality with applications to linear forms*, Pacific J. Math. 83 (1979) (central cube sections
  have volume `>= 2^D` in `[-1,1]^U`); E. Bombieri, J. D. Vaaler, *On Siegel's lemma*, Invent. Math. 73 (1983) (statement of
  Thm 2 used in Thm 4.1 quoted from memory; the search confirms it combines Minkowski's second theorem with Vaaler's bound).
* O. Klein, J. Slote, A. Volberg, H. Zhang, *Quantum and classical low-degree learning via a dimension-free Remez
  inequality*, arXiv 2301.01438.
* From memory, not re-read: Carbery-Wright (2001); O'Donnell, *Analysis of Boolean Functions* (2014) Thm 10.21
  (general hypercontractivity); Halasz (1977).
