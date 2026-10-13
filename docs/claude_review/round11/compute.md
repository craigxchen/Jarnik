# Round 11 — the finite test of ZE at model sizes

(skeleton; sections filled in progressively)

## 0. Summary of outcome

Task: the finite test of `round10/README.md` §8.4 for Hypothesis ZE (`determinant.md` §5.5): at model sizes,
exact bases of `V(A,m)`, real tangent cones, decomposition of the base locus `Z_P` near `e`, intersections, and
an adversarial search for fakes on `Z_P`. **No growth improvement is proved, and ZE is neither proved nor
refuted in the regime.** What is established:

1. **The README's model sizes are vacuous (proved + exact).** If the slab-box `A` contains a lattice zonotope
   `k_0 + sum [0,L_i] v_i` with `sum L_i >= m`, an explicit `c in V_H` (a product of finite differences) has
   `sum c_k w^k = w^{k_0} prod (w^{v_i} - 1)^{alpha_i}`, so every point of `Z_P` is degenerate (Prop. 2.1). For
   actual nondegenerate clusters in the regime this can never happen: it would violate a character inequality
   (Prop. 2.2; Prop. 6.1 for every admissible support). A model is faithful to the regime only if
   `zeta(A) < m <= rho(A)`. At `n - 1 = 2` no admissible model is faithful (Prop. 2.4, proved). At `n - 1 = 3`
   all 355 admissible slab-boxes of full-rank sign profiles (`T <= 4`, `L <= 3`) have `zeta = rho` (exact). At
   `n - 1 = 4`, all 74 distinct slab-boxes examined (12 random full-rank profiles, `T, L <= 3`) have `zeta = rho`.
   So ZE "holds" there only for character reasons.
2. **First faithful models (proved criterion + exact).** For a Hadamard profile,
   `zeta(A) <= ((n-1)/n) max_p h_A(beta_p)` (Prop. 3.1). At `n = 8` (Sylvester `H_8`, `T = 1`) the models
   `M8 = (L = 2, m = 8)` (`dim V = 42`) and `M8' = (L = 3, m = 11)` (`dim V = 7`) are faithful, with
   `V_H = V(A,m)` exactly.
3. **ZE holds in both, exactly and on the whole complex torus (Theorem 5.2).** Every point of `Bs(A, m)` lies on
   an explicit subtorus: for `M8'` a collision `delta_x = delta_y` (mod `2pi`); for `M8` a collision or one of six
   four-term relations such as `delta_0 + delta_7 = delta_1 + delta_6`. Method (Lemma 5.1): exact factorisation of
   the basis polynomials into two-term factors and irreducible factors that turn out to be **forms in `w - 1`**,
   then saturation of the ideal of these forms by the collision hyperplanes (Bayer's criterion, exact Groebner
   bases over `Q`), which ends in the unit ideal. Each single form has non-toric zeros through `e`; only the
   intersection is toric. No fake (nondegenerate) position lies on `Z_P` for these profiles.
4. **ZE-real fails in an underdetermined faithful model (Prop. 6.2, proved from exact data).** For the
   skew-double-core profile at `q = 7` (`T = 1`, `L = 2`, `m = 8`), `dim V = 1` and `Z_P` contains a real
   hypersurface through `e` with nondegenerate points arbitrarily close to `e` (the leading form takes both signs,
   exact witnesses, and is irreducible over `Q`). Faithfulness alone does not give ZE; the size of `V_H`
   (`>= |A|/2` in the regime) is essential.
5. **Reformulation (§5.3).** In `M8'`, `Bs(A,m) = e + V_C(Sigma)` where `Sigma` is a space of degree-`m` forms
   ("realisable symbols" `H` with `H(Delta)` supported in `A`). This suggests **ZE-braid**: the realisable
   symbols have no common zero off the braid arrangement. ZE-braid is purely algebraic, implies ZE on the complex
   torus, and is the recommended next target (§8).

Status: `evidence` for ZE in the overdetermined regime (two exact model theorems), together with an exact
counterexample to ZE-real for an underdetermined model, and proofs that the originally proposed model sizes
cannot discriminate.

## 1. Setup recalled, and three reformulations

Conventions. *Proposition/Lemma*: proved here with all quantifiers. *Exact*: integer or exact-rational
computation (python-flint `fmpz_mat`, `fmpz_mpoly`, sympy over `Q`). *Certified*: interval arithmetic
(flint `arb/acb`) with explicit error bounds. *Evidence*: floating point or finite search that no proof uses.
Scripts and outputs are in `compute_checks/` (`<script>.py` writes `<script>_out.txt`).

**Notation (as in `round10/determinant.md`).** `n` points, `d = n - 1`. A frequency is `k in Z^n_0`, written
in the coordinates `k' = (k_1, ..., k_d)`, `k_n = -sum k'`. Positions enter through `u_x = delta_x - delta_n`
(`x < n`), since `<k, delta> = <k', u>` on `Z^n_0`; write `w = e^{iu}` on the torus `T^d = (C^*)^d`, so that
`e^{i<k,delta>} = w^{k'}`. A point `delta` is *nondegenerate* iff `<v, u> != 0` for every nonzero `v in Z^d`,
i.e. iff `u` lies on no rational hyperplane (this is the definition of `determinant.md` §5.5 rewritten).

A squarefree sign column `beta in {+1,-1}^n` with `S = {x : beta_x = +1}` gives `<k, beta> = 2 k(S)` with
`k(S) = sum_{x in S} k_x` (because `sum k = 0`), and `k(S) = -k(S^c)`. So the slab `|<k,beta>| <= w` is
`|k(S)| <= L` with `L = w/2`, and `h_A(beta) = 2 width_A(k(S))`. A *profile* is a set of patterns `S`
(up to complement). The slab-box with box size `T` and slab half-width `L` is

```text
A = A(T, L, profile) = { k' in Z^d : |k_x| <= T (x < n),  |k(S)| <= L for every S in the profile }.
```

(One `L` for all columns, as in Theorem 2.1, where one `w` serves all columns.) `V(A, m)` is the space of
`c in Q^A` annihilating every polynomial of degree `< m` in `k'`; `Poly_{<m}|_A` is its annihilator.

**Reformulation R1 (base locus = exponential interpolation).** Since `V(A,m)` is exactly the annihilator of
`Poly_{<m}|_A` in `Q^A`,

```text
Bs(A, m) := { w in T^d : sum_k c_k w^k = 0 for all c in V(A,m) }
          = { w : the function k -> w^k on A coincides with a polynomial of degree < m in k }.
```

`Z_P` (for `V_H = V(A,m)`, which is the case in every model below, §2) is `Bs(A,m)` intersected with the real
torus `|w_x| = 1`.

**Reformulation R2 (tangent cone = powers of linear forms in the leading-form ideal).** Let `I(A)` be the
ideal of polynomials vanishing on `A` and `J = gr I(A)` the ideal of their top-degree forms; `S = Q[k']`.
Then `V(A,m) != 0` iff `(S/J)_m != 0` (the Hilbert function `H_A(j) = dim Poly_{<=j}|_A` equals
`sum_{i<=j} dim (S/J)_i`, and `S/J` is generated in degree 0). Put `rho(A) = max{j : (S/J)_j != 0}`; then
`V(A, m) != 0` iff `m <= rho(A)`. The forms `Phi_c(u) = sum_k c_k <k,u>^m`, `c in V(A,m)`, span the degree-`m`
part of the Macaulay inverse system of `J` (the least space of de Boor-Ron), and

```text
TC(A, m) := { u in C^d : Phi_c(u) = 0 for all c in V(A,m) } = { u : <., u>^m in J_m }
          = { u : the polynomial <k,u>^m agrees on A with a polynomial of degree < m }.
```

*Proof.* `Phi_c(u) = sum_k c_k <k,u>^m` vanishes for all `c in V(A,m)` iff `<.,u>^m|_A` lies in the annihilator
of `V(A,m)`, which is `Poly_{<m}|_A`; i.e. iff `<k,u>^m - g(k)` vanishes on `A` for some `g` of degree `< m`,
i.e. iff `<.,u>^m` is the top form of an element of `I(A)` of degree `m`. QED

**Reformulation R3 (parity).** Every `A` here is symmetric (`A = -A`), and `k -> -k` preserves degree, so
`V(A,m) = V^+ + V^-` (even and odd `c`). On the real torus the equations of `Z_P` are
`sum_k c_k cos<k,u> = 0` (`c in V^+`) and `sum_k c_k sin<k,u> = 0` (`c in V^-`).

## 2. The zonotope certificate, and which models are faithful to the regime

The first finding of the computations (§3) is structural, and it changes what a "model" must be. It is
proved here.

**Proposition 2.1 (zonotope certificate).** Let `v_1, ..., v_s in Z^d \ {0}` (any number, dependence and
repetition allowed), integers `L_i >= 1`, and `k_0 in Z^d`, such that every point `k_0 + sum_i t_i v_i` with
integers `0 <= t_i <= L_i` lies in `A`. Suppose `sum_i L_i >= m`, choose integers `0 <= alpha_i <= L_i` with
`sum alpha_i = m`, and let `c` be the coefficient vector of the Laurent polynomial
`w^{k_0} prod_i (w^{v_i} - 1)^{alpha_i}`. Then

* (a) `supp c` is contained in `A`, `c in V(A, m)`, and `|c|_inf <= |c|_1 <= 2^m`, so `c in V_H` whenever
  `H >= 2^m` (true for the `H = (U T^{m-1})^2` of Theorem 2.1 once `T >= 2`);
* (b) `sum_k c_k w^k = w^{k_0} prod_i (w^{v_i} - 1)^{alpha_i}` for every `w in (C^*)^d`.

Hence every point of `Bs(A, m)` (and of `Z_P`) satisfies `w^{v_i} = 1` for some `i`; in particular **no
nondegenerate point of the whole complex torus lies in `Bs(A,m)`**, and ZE holds for `(A, m)`.

*Proof.* Expanding, `c_k = sum (prod_i (-1)^{alpha_i - t_i} binom(alpha_i, t_i))` over `t <= alpha` with
`k_0 + sum t_i v_i = k`; so `supp c` is inside `A` and `|c|_1 <= prod_i 2^{alpha_i} = 2^m`. For a polynomial
`P` of degree `< m` put `Q(t) = P(k_0 + sum t_i v_i)`, a polynomial in `t in Z^s` of total degree `< m`. Then
`sum_k c_k P(k) = (prod_i Delta_i^{alpha_i} Q)(0)`, with `Delta_i` the forward difference in `t_i`. The
operator `prod Delta_i^{alpha_i}` kills every monomial `t^gamma` unless `gamma_i >= alpha_i` for all `i`,
which needs `|gamma| >= m`; so the sum is `0` and `c in V(A,m)`. (b) is the definition of `c`. If `w` is in
`Bs(A,m)` the right side of (b) vanishes, so some `w^{v_i} = 1`; for `w = e^{iu}` near `e` this means
`<v_i, u> = 0`. QED

`H >= 2^m`: `U >= 1` and `T^{2(m-1)} >= 2^m` for `T >= 2`, `m >= 2`. For `T = 1` one checks `H` directly.

Put `zeta(A) = max sum_i L_i` over such lattice zonotopes in `A`, and recall `rho(A)` from R2. Proposition 2.1
gives `V(A, zeta(A)) != 0`, so **`zeta(A) <= rho(A)`**, and ZE is automatic for every `m <= zeta(A)`.

**Proposition 2.2 (zonotopes are invisible for actual nondegenerate clusters in the regime).** Let a cluster
satisfy the hypotheses (4.1) of `determinant.md` Theorem 4.1, with nondegenerate positions, and let `(A, m)`
be its own slab-box (Theorem 2.1, `T = n`). Then `zeta(A) < m`.

*Proof.* Suppose `Z = k_0 + sum [0, L_i] v_i` lies in `A` with `sum L_i >= m` (as in Prop. 2.1). View
`v in Z^d` as `(v, -sum v) in Z^n_0`. Widths of zonotopes add: `h_Z(beta_p) = sum_i L_i |<v_i, beta_p>|`.
Since `Z` lies in `A`, `h_Z(beta_p) <= h_A(beta_p) <= 2w` for every `p` (Theorem 2.1(b)). Multiply by
`omega_p`, sum over `p`, and use `sum omega_p <= 1`, `2w <= m/2` (Theorem 2.1(d)) and
`omega_p |<v, beta_p>| = (2 log p / W) |<v, a_p>|` (since `<v, a_p> = (e_p/2) <v, beta_p>` on `Z^n_0`):

```text
sum_i L_i X_{v_i}  <=  m W / 4  <=  (sum_i L_i) W / 4,        X_v := sum_p log p |<v, a_p>| .
```

So some `X_{v_i} <= W/4`. On the other hand `<v_i, delta> != 0` (nondegeneracy), and the character
inequality (`determinant.md` Prop. 5.3(iii), with the constant corrected in `referee_determinant.md` §6) gives
`X_{v_i} >= W/2 - 2 log(C |v_i|_1 / (2 sqrt 2))`. Both `k_0` and `k_0 + L_i v_i` lie in the box, so
`|v_i|_inf <= 2T = 2n` in the first `d` coordinates and `|v_i|_1 <= 4n(n-1) = 2D`. Hence
`W <= 8 log(C D / sqrt 2) < 8 log^+ C + 16 ln n + 3`, contradicting `W >= 16 ln n + 8 log^+ C + 14`. QED

So Proposition 2.1 is exactly the zonotope extension of `determinant.md` Prop. 5.3 (parallelepipeds): it
proves ZE for a profile, but never for the slab-box of an actual nondegenerate cluster in the regime, where
it would contradict the character inequalities. **A model `(A, m)` is therefore faithful to the regime only if
`zeta(A) < m <= rho(A)`.** Models with `m <= zeta(A)` say nothing about ZE in the regime: they are settled by a
character. In particular every model with `m <= 2T` (the coordinate segments have `2T + 1` points) is unfaithful.

**Proposition 2.3 (section principle).** Let `v_1, ..., v_j in Z^d` be linearly independent, `k_0 in Z^d`, and
`A' = {t in Z^j : k_0 + sum t_i v_i in A}`. Let `psi(w) = (w^{v_1}, ..., w^{v_j})`. Then
`psi(Bs(A, m))` is contained in `Bs(A', m)`. Consequently, if near `e` the set `Bs(A', m)` lies in finitely many
proper subtori `{z : z^b = 1}` (`b in Z^j \ 0`), then every `w in Bs(A,m)` near `e` satisfies
`w^{sum b_i v_i} = 1` and is degenerate.

*Proof.* If `w^k = g(k)` on `A` with `deg g < m`, then `psi(w)^t = w^{-k_0} g(k_0 + sum t_i v_i)` on `A'`, a
polynomial of degree `< m` in `t`. QED

In the regime a `j`-dimensional section lies in the box, so its coordinate functionals give a complete
intersection in `J(A')` and `rho(A') <= 2Tj`; a useful section therefore needs `j >= m/(2T) ~ 0.37 n`. This is
the same obstruction as for zonotopes.

**Proposition 2.4 (n - 1 = 2: no faithful model exists).** Let `n = 3` and let `(A, m)` be a slab-box model
with one `L` for all columns, with (i) slabs cutting, (ii) `h_A(beta_p) < m` for every column and (iii)
`V(A,m) != 0`. Then `m <= zeta(A) = rho(A)`, so every such model is settled by Proposition 2.1.

*Proof.* The non-constant patterns are `{1}`, `{2}`, `{3}` (functionals `k_1`, `k_2`, `k_3 = -(k_1+k_2)`).
*Two or more cutting patterns.* Two of these functionals, say `l_1, l_2`, are independent and each takes at
most `2L + 1` values on `A`, so `prod_{|s| <= L}(l_i - s)` vanishes on `A` and `l_i^{2L+1} in J`. The ideal
`(l_1^{2L+1}, l_2^{2L+1})` contains every form of degree `>= 4L + 1`, so `rho(A) <= 4L`. But (ii) needs
`m > h_A(beta_p) = 2 width(l_i)`; and the width is `2L` once the slab cuts (the box is symmetric and contains
`k` with `l_i(k) = L`, else the slab would not cut), so `m >= 4L + 1 > rho(A)`, contradicting (iii).
*One cutting pattern `{1}` or `{2}`.* `A = [-L', L'] x [-T, T]` (`L' = min(L,T)`) is a lattice parallelepiped,
`zeta >= 2L' + 2T`, and `rho <= 2L' + 2T` by the complete intersection of the two coordinate functionals.
*One cutting pattern `{3}` (the base point).* `A` is the hexagon `{|k_1|, |k_2| <= T, |k_1 + k_2| <= L}`,
`L < 2T`. It contains the zonotope `(-T, T - L) + [0, 2T - L](1,-1) + [0, L](1,0) + [0, L](0,1)` (check the
three functionals on the generators), so `zeta >= 2T + L`. Conversely, fibre `A` by the lines
`k_1 + k_2 = s`, `|s| <= L`, which carry `2T + 1 - |s|` points; order them `s_0 = 0, s_1 = 1, s_2 = -1, ...`
and interpolate line by line (Newton form `sum_j prod_{i<j}(k_1 + k_2 - s_i) q_j`, `q_j` of degree
`2T - |s_j|` in `k_1`). The `j`-th term has degree `j + 2T - ceil(j/2) <= 2T + L`, so every function on `A` is a
polynomial of degree `<= 2T + L` and `rho <= 2T + L`. QED

The proof of the "two cutting patterns" case is the general mechanism behind the next observation: for
full-rank profiles `rho(A)` is at most the sum of the ranges of any `d` independent slab functionals, so
condition (ii) needs `d >= 3`.

## 3. Admissible models at n-1 = 2, 3, 4: exact bases, rho, zeta

**Conditions.** A model is `(n, profile, T, L, m)` with one slab half-width `L` for all columns and

* (i) the slabs cut (`A != B_T`);
* (ii) `h_A(beta_p) = 2 range_A(k(S_p)) < m` for every column (the README's "`2w < m`"), or the weaker
  (ii') `range_A(k(S_p)) < m` (this is what puts the one-parameter subgroup of every column direction into
  `Bs`, which is the geometric content of supercriticality; Prop. 5.2 of `determinant.md`);
* (iii) `V(A,m) != 0`, i.e. `m <= rho(A)`;
* (iv) faithfulness: `zeta(A) < m` (§2);
* (v) full rank: the functionals `k(S_p)` span `(R^d)^*`, as for the nondegenerate subclusters of Lemma 3.2.

`rho(A)` is computed exactly (rank of the integer moment matrices; `models.py:rho`), or first by rank modulo
`p = 2^61 - 1` for screening (`scan.py`; `rank_p <= rank_Q`, so this is an upper bound for `rho`, confirmed
exactly for every instance used below). `zeta(A)` is computed by an exhaustive search over lattice zonotopes
(`zono.py`): since `A` is the set of lattice points of the polytope `Q = {|k_x| <= T, |k(S)| <= L}`, a zonotope
lies in `A` iff each of these constraint functionals stays within its bounds on `Z`, a finite check; non-primitive and parallel
generators are merged; and the search prunes with the valid bound `sum L_i <= sum_x (2T - used_x)`. The search
returns a witness zonotope, so `zeta >= s` is certified by the witness; `zeta < s` is certified by exhaustion.

**n - 1 = 2.** Proposition 2.4: no faithful model. Exact base loci of the hexagon models
(`n3_hexagon.py`, output `n3_hexagon_out.txt`, 9 instances with `T <= 6`, `L <= 3`): the curve part of `Bs`
(exact gcd of the Laurent polynomials, factored over `Q`) is `(w_1 - w_2)^j`, times `(w_1 - 1)(w_2 - 1)` when
`m <= 2T + 1`; the tangent cones are `(u_1 - u_2)^j [u_1 u_2]`; for instance at `T = 4, L = 1, m = 9` the single
kernel vector is exactly `(w_1 - 1)(w_2 - 1)(w_1 - w_2)^7`, the zonotope certificate of Prop. 2.1. All curve
components are subtori through `e`. (Isolated points were only located through resultant factors and not
classified; this is not needed, since Proposition 2.4 already puts all of `Bs` in finitely many subtori.)

**n - 1 = 3 (exhaustive over sign profiles).** All `2^7 - 1` sets of non-constant patterns at `n = 4`
(`{x}` for the four points and the three balanced cuts), restricted to the full-rank ones, with
`1 <= T <= 4`, `1 <= L <= 3` (`scan.py 4`, `faithful.py`; outputs `scan4_out.txt`, `faithful4_out.txt`):

* 980 instances satisfy (i), (iii), (v) and (ii'); 339 of them satisfy (ii); 355 distinct sets `A` occur;
* **for every one of the 355, `zeta(A) = rho(A)`** (exact). So there is no faithful model, and ZE holds in each
  of them by Proposition 2.1 (an explicit zonotope with `sum L_i = rho >= m` is printed for each).
* The Hadamard profile at `n = 4`, which is also the skew double core at `q = 3` (its six cuts are the three
  balanced cuts, each with both orientations), gives `A ~ {(a,b,c) in [-L,L]^3 : a+b+c even}` with
  `rho = zeta = 3L < 4L = h_A(beta_p)`: condition (ii) is infeasible there, and only (ii') holds.

**n - 1 = 4.** Twelve random full-rank sign profiles (`r` between `n` and `3n` random columns), `T <= 3`,
`L <= 3` (`scan.py 5`, `faithful.py`; `scan5_out.txt`, `faithful5_out.txt`): 101 instances with (i), (iii),
(v) and (ii') (34 with (ii)); among the distinct sets `A` processed, every one has `zeta = rho` exactly
(witness zonotopes); for one of them, `((0,2),(0,3),(1),(2),(2,4),(3,4),(4))` with `T = L = 3`, `rho = 13`, the
search timed out and the witness comes from the factorisation of the basis instead (§3.2). **All 74 distinct
sets have `zeta = rho`.**
Upper bounds for `zeta` are also available from an LP: every zonotope in `A` satisfies
`sum_i L_i |lambda(v_i)| <= 2 b_lambda` for each constraint `|lambda| <= b_lambda`; a rational dual solution `y`
with `sum_lambda y_lambda |lambda(v)| >= 1` for all candidate generators `v` gives `zeta <= sum 2 b_lambda y_lambda`
(`lpzeta.py`, verified in exact rationals).

Random symmetric polytopes `{k in [-T,T]^d : |lambda_j(k)| <= b_j}` with `lambda_j in {-1,0,1}^d`
(`gap_search.py`; `gap4_T1_out.txt`, `gap4_T2_out.txt`, 125 instances in `d = 3, 4`): `zeta = rho` in every case.

**Reading.** At `n - 1 <= 4` the slab-box is "zonotopal": the top degree of its interpolation space is
attained by a product of finite differences. Every such model is decided by a character-type certificate
(Prop. 2.1), which Proposition 2.2 shows is unavailable in the regime. **The finite test proposed in the README
(§8.4) is therefore vacuous at `n - 1 = 2, 3, 4`**: ZE holds there for a reason that cannot survive into the
regime. A meaningful test must go to a size where `zeta < rho`. §3.1 gives the first such models.

### 3.1 The first faithful models: the Hadamard profile at n = 8

**Proposition 3.1 (Hadamard profiles have zeta < h).** Let `H` be a Hadamard matrix of order `n` with first row
`1`, rows `h_p`, and let `A` be any slab-box (any `T`, `L`) cut by the `n - 1` columns `beta_p = h_p`
(`p >= 2`). Then `zeta(A) <= ((n-1)/n) max_p h_A(h_p)`.

*Proof.* For `0 != v in Z^n_0`: `sum_{p=1}^n <v,h_p>^2 = n |v|_2^2` (orthogonality of the rows),
`<v, h_1> = 0`, and `|<v,h_p>| <= |v|_1 <= |v|_2^2` (integer entries). Hence
`sum_{p>=2} |<v,h_p>| >= sum_p <v,h_p>^2 / |v|_1 >= n`. For a zonotope `Z = k_0 + sum [0,L_i] v_i` in `A`,
`sum_i L_i |<v_i, h_p>| = h_Z(h_p) <= h_A(h_p)`; summing over `p >= 2`,
`n sum_i L_i <= sum_{p>=2} h_A(h_p) <= (n-1) max_p h_A(h_p)`. QED

So for Hadamard profiles every `m` in `(((n-1)/n) max_p h_A(h_p), rho(A)]` is faithful, and condition (ii)
alone (`m > max_p h_A(h_p)`) implies faithfulness. At `n = 4` the window is empty (`rho = 3L`). At `n = 8`
(Sylvester `H_8`; the seven cuts are the affine hyperplanes of `F_2^3`), `T = 1`:

| `L` | `|A|` | `range k(S_p)` | `max h_A(beta_p)` | `rho(A)` (exact) | `zeta <=` | faithful `m` |
| --- | --- | --- | --- | --- | --- | --- |
| 2 | 1037 | 4 | 8 | 8 | 7 | **8** (ii'; not (ii)) |
| 3 | 1877 | 6 | 12 | 11 | 10 | **11** (ii'; not (ii)) |

**The model `M8 = (H_8, T = 1, L = 2, m = 8)`** (`had8.py`, output `had8_T1_L2_m8_out.txt`, exact):
`rank Poly_{<8}|_A = 995`, `rank Poly_{<9}|_A = 1037 = |A|`, so `dim V(A,8) = 42` and `rho(A) = 8`. The
nullspace basis, LLL-reduced, has sup-norm `6`, against the Siegel height `H = (|A| T^{m-1})^2 ~ 1.1 * 10^6`; so
`V_H = V(A, 8)`. All 42 vectors are verified exactly. The cheap directions (range `< 8`) include every coordinate
direction, every column direction `k(S_p)` (5 values) and many short sums, as in the regime.

**The model `M8' = (H_8, T = 1, L = 3, m = 11)`** (`had8.py 1 3 11`, output `had8_T1_L3_m11_out.txt`, exact):
`|A| = 1877`, `rank Poly_{<11}|_A = 1870`, `rank Poly_{<12}|_A = 1877`, so `rho = 11` and `dim V(A,11) = 7`;
LLL basis sup-norm `4` against `H ~ 3.5 * 10^6`, so `V_H = V`.

**Other profiles at `n = 8`, `T = 1`** (`model8.py`, exact `rho`, exact LP bound for `zeta`):

| profile | `L` | `|A|` | `rho` | `zeta <=` | faithful `m` | `dim V` |
| --- | --- | --- | --- | --- | --- | --- |
| skew double core `q = 7` (14 cuts) | 2 | 757 | 8 | 7 | 8 | **1** |
| skew double core `q = 7` | 3 | 1765 | 10 | `= 10` (witness zonotope, `sdc7_L3_zeta.py`) | none | - |
| round-5 Paley-I `q = 7` (15 cuts) | 2 | 807 | 7 | 7 | none | - |
| random sign profile (seed 1, 12 cuts) | 2 | 1223 | 9 | 19/2 | not certified | - |
| random sign profile (seed 2, 21 cuts) | 2 | 905 | 8 | 17/2 | not certified | - |

### 3.2 The one `n - 1 = 4` instance not settled by the search

For `((0,2),(0,3),(1),(2),(2,4),(3,4),(4))`, `T = L = 3` (`|A| = 811`, `rho = 13`, LP bound `zeta <= 27/2`), the
exhaustive zonotope search did not finish in 2 hours. The exact basis of `V(A, 13)` (`n5_hard_basis.py`,
`cone_reduce.py`; output `n5_hard_basis_out.txt`; `dim V = 5`, sup-norm 15) settles it: **each of the 5 basis
polynomials is a product of a monomial and two-term factors `w^a - w^b` of total order 13**, i.e. a zonotope
certificate (the Newton polytope of a product is the Minkowski sum, so the zonotope lies in `conv A` and its
lattice points in `A`). Hence `zeta = rho = 13`, and this model is unfaithful too.

## 4. Tangent cones: real zeros of Phi_c

**Unfaithful models (all admissible models at n - 1 = 2, 3, and all those found at n - 1 = 4): decided exactly.**
If `A` contains a zonotope `k_0 + sum [0,L_i] v_i` with `sum L_i >= m`, the certificate `c` of Prop. 2.1 has
`sum_k c_k e^{<k,z>} = e^{<k_0,z>} prod_i (e^{<v_i,z>} - 1)^{alpha_i}`, whose lowest-order term is
`prod_i <v_i, z>^{alpha_i}` (degree `m`). Comparing degree-`m` parts of `sum_k c_k e^{<k,z>}`:

```text
Phi_c(u) = sum_k c_k <k,u>^m = m! prod_i <v_i, u>^{alpha_i}        (a product of rational linear forms).
```

Hence **every common real zero of the `Phi_c` lies on one of the rational hyperplanes `<v_i, u> = 0`**, i.e.
is degenerate. This answers test (2) of the README positively, exactly, for every model of §3 with
`n - 1 <= 3` and every `n - 1 = 4` model examined, and equally for test (3): `Bs(A, m)` lies in the finite
union of subtori `{w^{v_i} = 1}` through `e` (Prop. 2.1(b)). For the hexagons (`n - 1 = 2`) the tangent cones
were also computed directly (gcd of the binary forms, factored over `Q`, `n3_hexagon_out.txt`): they are
`(u_1 - u_2)^j`, times `u_1 u_2` when `m <= 2T + 1`. By Proposition 2.2 none of this transfers to the regime.

**The faithful model M8** (`had8_tc.py`, output `had8_T1_L2_m8_tc_out.txt`, exact). The 42 forms `Phi_c`
are octics in 7 variables with integer coefficients (LLL-reduced, max coefficient `120960`). Since
`k_x^3 - k_x` vanishes on `A`, only monomials `u^alpha` with `alpha <= 2` (entrywise) occur, and every
monomial that occurs has **at least 4 distinct variables**. Exactly:

* all `7 + 21 + 35` coordinate subspaces of dimension `<= 3` lie in `TC` (Prop. 5.1 of `determinant.md`:
  `2T|S| = 6 < 8`);
* exactly `7` of the `35` coordinate 4-planes lie in `TC`, namely `span{e_x : x in S_p}` for the seven cuts
  `S_p` not containing the base point (there `|k(S_p)| <= 2` removes the corners of the 4-cube, and the image
  has `rho <= 7`);
* no coordinate subspace of dimension `>= 5` lies in `TC`.

The tangent cones of the faithful models are not computed separately: Theorem 5.2 below decides the whole
base locus, and the tangent cone of `Bs` at `e` is then contained in the union of the tangent spaces of the
subtori listed there (all rational). (A direct numerical search for real zeros of the 42 octics,
`tc_search.py`, converged slowly to singular components and is not used.)

## 5. The base locus near the identity

### 5.1 The cone reduction (exact method)

**Lemma 5.1 (cone reduction).** Let `c_1, ..., c_D` be a basis of `V(A, m)` and `P_j(w) = sum_k c_{j,k} w^k`.
Factor each `P_j` over `Q` as `P_j = (const) w^{a_j} prod_b b(w)^{e_b} prod_f f(w)^{e_f}`, where the `b` are the
two-term factors and the `f` the irreducible factors with at least three terms. Suppose every two-term factor
has the form `w^a - w^{a'}` (`a != a'`). For a subset `J` of indices put `G_j(x) = prod_f f(1 + x)^{e_f}`
(`x = w - 1`). Then every `w in Bs(A, m)` satisfies `w^{a - a'} = 1` for one of the two-term factors, or
`w - 1 in V(G_j : j in J)`. If moreover every `G_j` (`j in J`) is a form in `x` and a chain of saturations of
the homogeneous ideal `(G_j : j in J)` by linear forms `l_1, ..., l_K` ends in the unit ideal, then
`V(G_j) subset V(l_1) cup ... cup V(l_K)`.

*Proof.* `Bs(A,m)` is the common zero set of the `P_j` (they span `V(A,m)`), and the monomial is a unit on the
torus. If no two-term factor vanishes at `w`, then for each `j` some `f` with `e_f > 0` vanishes, i.e.
`G_j(w - 1) = 0`. For the second statement, `V(I : l^inf)` contains `V(I) \ V(l)`, so if
`I_0 = (G_j)`, `I_{k+1} = I_k : l_k^inf` and `I_K = (1)`, then `V(I_0)` lies in `V(l_1) cup ... cup V(l_K)`.
QED

The saturations are computed by Bayer's criterion: for a homogeneous ideal and the grevlex order in which `l`
is the last variable (after a linear change of coordinates), the reduced Groebner basis divided by the highest
powers of `l` is a Groebner basis of `I : l^inf` (Eisenbud, *Commutative Algebra*, Prop. 15.12). Groebner bases
are computed over `Q` with sympy (exact rational arithmetic); factorisations with flint (`fmpz_mpoly.factor`).
Scripts: `cone_reduce.py` (factorisation, form check), `cone_sat.py` (saturations).

### 5.2 Theorem: ZE holds, globally, in the two faithful Hadamard models

**Theorem 5.2 (exact computation).** Let `A_L` be the slab-box of the Sylvester-Hadamard profile at `n = 8`
(the seven affine hyperplanes of `F_2^3` as cuts), box `T = 1`, slab half-width `L`.

* (a) `M8' = (L = 3, m = 11)`: faithful (`zeta <= 10 < 11 = rho`), `dim V = 7`, `V_H = V`. Every point of
  `Bs(A_3, 11)` (in the whole complex torus) satisfies `w_x = 1` for some `x` or `w_x = w_y` for some `x != y`.
  In terms of the positions: **every `delta in Z_P` has two coinciding coordinates `delta_x = delta_y`
  (mod `2 pi`)**, `x, y in {1, ..., 8}`.
* (b) `M8 = (L = 2, m = 8)`: faithful (`zeta <= 7 < 8 = rho`), `dim V = 42`, `V_H = V`. Every point of
  `Bs(A_2, 8)` lies on one of the subtori `{w_x = 1}`, `{w_x = w_y}`, or `{w^a = w^{a'}}` for the 27 two-term
  factors listed in `had8_T1_L2_m8_cone_out.txt` (for example `w_0 = w_1 w_6`, i.e.
  `delta_0 + delta_7 = delta_1 + delta_6`).

In both models **no nondegenerate point of the complex torus lies in the base locus**, so ZE holds (not only near
`e` and not only on the real torus), and no gap-admissible fake position lies on `Z_P`.

*Computation.* (a) (`cone_reduce.py`, `cone_sat.py`; outputs `had8_T1_L3_m11_cone_out.txt`,
`had8_T1_L3_m11_out.txt`.) Each of the 7 basis polynomials factors as a monomial, two-term factors of total order
`8` or `7` at `e` (all of the form `w_x - 1`, `w_0 - w_1`, `w_3 - w_4`), and **one** irreducible factor which is a
form in `x = w - 1`: three cubics with 6 terms and four quartics with 12 terms (12 and 24 terms in `w`); for
instance `x_2x_3x_4 + x_2x_3x_5 - x_2x_3x_6 - x_2x_4x_5 - x_3x_4x_5 + x_4x_5x_6`. The orders add up to `m = 11`.
The ideal of the 7 forms has a 30-element Groebner basis; its Hilbert polynomial has degree 3 (a projective
threefold of degree 5, containing exactly two of the coordinate 4-spaces). Saturating successively by
`x_0, ..., x_6, x_0 - x_1, ..., x_2 - x_3` gives the unit ideal. Lemma 5.1 applies with `J` = all.
(b) (`had8_T1_L2_m8_cone_out.txt`.) Of the 42 basis polynomials, 24 have only form-type non-binomial factors
(11 of shape "binomials of order 2 times two cubic forms", 13 of shape "binomials of order 4 times one quartic
form"); the other 18 have a non-homogeneous factor and are not used. The ideal of the 24 forms saturates to the
unit ideal after `x_0, ..., x_6, x_0 - x_1, ..., x_2 - x_3` (237 s). Lemma 5.1 applies with `J` = these 24. QED

**Remarks.**

* The non-binomial factors are not subtori and their zero sets are not toric (they are cones in `x = w - 1`, not
  in `log w`). Each single `P_c` therefore has non-toric components through `e` with nondegenerate real points:
  the README's failure mode exists *for each form separately*. It is the **intersection** that is toric. This is
  exactly the unlikely-intersection mechanism that ZE asserts, now observed in a model where no zonotope
  certificate exists.
* The degeneracies produced are of the simplest kind: coincidences `delta_x = delta_y` (and, in (b), the
  four-term relations of the two-term factors). For an actual cluster, `delta_x = delta_y` is impossible
  (distinct lattice points), and the relations `w^a = w^{a'}` are excluded by the gap lemma whenever the
  subcluster is nondegenerate (Lemma 3.2 of `determinant.md`). So in these models the rigidity relation (4.2)
  would by itself contradict the existence of a nondegenerate cluster with this profile.
* Numerical cross-check (`explain.py`, evidence, superseded by the theorem): 150 Levenberg-Marquardt zeros on
  the spheres `|tan(u/2)| = 0.3, 0.1, 0.02`; 117 lie (to `1e-8`) on subtori verified exactly in `Bs` by the
  push-forward test; the remaining 30 are poorly converged points visibly on rational subspaces (coordinates
  equal to `1e-6` accuracy), where LLL relation detection at tolerance `1e-7` failed. No zero without relations
  was found, consistent with Theorem 5.2.

**A third support for the same profile (README "freedom").** Adding `|k_8| <= 1` for the base point (the
symmetric box, `sym8.py`, output `sym8_out.txt`): `L = 2` gives `|A| = 869`, `rho = 8`, exact LP bound
`zeta <= 7`, so `m = 8` is faithful, with `dim V = 28` (sup-norm 6, `V_H = V`). Exactly
(`had8sym_T1_L2_m8_cone_out.txt`): **all 28** basis polynomials have the form `w^a H(w - 1)` with `H` a form of
degree 8 and **all** their two-term factors are collisions (`w_x - 1`, `w_x - w_y`). The saturation of the
ideal of the forms (all 28, or the 23 with non-binomial factors of degree `<= 4`) by the collision hyperplanes
was still running when this write-up was closed (`cone_sat.py`; sympy's Groebner bases are slow on these
octics), so for this support the conclusion "`Bs` lies on the braid arrangement" is **not** established.
`L = 3` is not faithful (`zeta <= 8 = rho` by LP; not decided).

### 5.3 The mechanism, and an algebraic reformulation suggested by it

In `M8'` every basis polynomial has the form `P_c(w) = w^a H_c(w - 1)` with `H_c` a **form of degree exactly
`m`** in `x = w - 1`: the two-term factors there are `w_x - 1 = x_x`, `w_0 - w_1 = x_0 - x_1`, `w_3 - w_4`, all
linear forms in `x`. (In `M8` the 24 basis polynomials used are of this form *except* for six two-term factors of
the type `w_0 - w_1 w_6`, which are not linear in `x`; this is why Theorem 5.2(b) lists these extra subtori,
e.g. `delta_0 + delta_7 = delta_1 + delta_6`.) Equivalently `c = tau_a H_c(Delta)`: a homogeneous polynomial of
degree `m` in the coordinate forward differences `Delta_x`, translated by `a`. For such `c`:

* `Phi_c(u) = m! H_c(u)`, so the leading form *is* the symbol `H_c`, and
* `{w : P_c(w) = 0} = {w : w - 1 in V(H_c)}`: the zero set of `P_c` is the **translate by `e` of the complex
  tangent cone** of `P_c`, read in the coordinates `x = w - 1`.

So, writing `Sigma_A = {H form of degree m : tau_a H(Delta) is supported in A for some a}` (the "realisable
symbols"),

```text
Bs(A, m)  is contained in  e + V_C(Sigma_A)      (x = w - 1),
```

and in `M8'` the realisable symbols of the 7 basis vectors already have `V_C` contained in the arrangement
`{x_x = 0} cup {x_x = x_y}` (Theorem 5.2(a)); in `M8` the same holds for the non-binomial parts. With
`x_n := 0` for the base point this is the **braid arrangement** `{z_x = z_y}` in the variables
`z_x = e^{i delta_x}`, i.e. the collision locus. This suggests a purely algebraic form of ZE that implies it
(and is stronger: it is a statement on the complex torus):

> **ZE-braid(A, m).** The degree-`m` forms `H` in `x_1, ..., x_d` whose difference operators `H(Delta)` fit
> (after translation) inside `A` have no common zero off the braid arrangement
> `{x_x = 0} cup {x_x = x_y}`.

ZE-braid implies ZE (Lemma 5.1 with all two-term factors linear) and even rules out every non-collision point of
`Z_P`. A product `prod (x_{x_i} - x_{y_i})` is realisable only through a zonotope (Prop. 2.1) and is therefore
unavailable in the regime (Prop. 2.2); the models show that **non-product symbols** do the job instead. It is a
statement about which forms of degree `m ~ 0.74 n^2` are "supported" on the slab polytope, i.e. about the
interplay between the Newton basis `binom(k + T, alpha)` and the slab cuts, and it no longer involves exponentials.
I have not proved it for any family; I record it as the precise next target.


## 6. Intersections over (T, w) and over sub-clusters

**Proposition 6.1 (only faithful supports may be intersected).** Let `A'` be any finite subset of `Z^d` and `m'`
an integer such that the comparison in the proof of Theorem 4.1 of `determinant.md` goes through, i.e.
`sum_p omega_p h_{A'}(beta_p) <= m'/2` (plus the size conditions on `W`). Then for an actual cluster with
nondegenerate positions, `zeta(A') < m'`.

*Proof.* Word for word as Proposition 2.2, with `2w` replaced by `sum_p omega_p h_{A'}(beta_p)`. QED

So the freedom of `determinant.md` open point 2 (other `T`, other `w`, unions of supports) can never be
exploited through a zonotope, i.e. through a character; every support usable in the regime is faithful.

**At model size.**

* *Over `(T, w)`.* For the Hadamard profile at `n = 8`, the two faithful supports `M8` and `M8'` each already
  force degeneracy (Theorem 5.2), so their intersection adds nothing. All supports at `n - 1 <= 4` are
  unfaithful (§3) and are excluded by Proposition 6.1.
* *Over sub-clusters.* A sub-cluster `J` of an `n = 8` model has `|J| <= 7`. Its own slab-box (profile
  restricted to `J`, same `T`, `L`) is zonotopal whenever `|J| <= 5` (§3: every admissible model at
  `n <= 5` examined has `zeta = rho`), so at model size it would "prove" degeneracy of `delta_J` for free. By
  Proposition 6.1 this is an artefact: in the regime a sub-cluster in the regime is itself faithful. Hence
  sub-cluster intersections carry no information at the sizes that can be computed, and I do not use them.

## 6bis. A codimension-one faithful model: the skew double core at q = 7, where ZE-real fails

The skew double core of `round5/exact.md` Thm 4.1 at `q = 7` (`n = 8`, the 14 balanced cuts `S_a`, `S'_a`; this
is the profile of the skew-double-core fakes) gives, with `T = 1`, `L = 2` (`model8.py sdc7 1 2`, output
`sdc7_T1_L2_m8_out.txt`, exact): `|A| = 757`, `rho(A) = 8`, the exact LP dual certificate gives `zeta <= 7`
(`lpzeta.py`), so `m = 8` is faithful; ranges of the 14 cut functionals are `4` ((ii') holds, (ii) fails as for
all `T = 1` models at `n = 8`); full rank. **But `dim V(A, 8) = 1`**: one vector `c`, sup-norm `42`, support
`687`, even (`c_{-k} = c_k`).

**Proposition 6.2 (ZE-real fails for this model).** For every `eps > 0` there is a nondegenerate `u` with
`|u| < eps` and `sum_k c_k e^{i<k,u>} = 0`, i.e. a nondegenerate point of `Z_P` within `eps` of `e`.

*Exact inputs* (`sdc7_refute.py`, `sdc7_tc.py`; outputs `sdc7_refute_out.txt`, `sdc7_tc_out.txt`):
(1) `c` is even, so `F(u) := sum_k c_k e^{i<k,u>} = sum_k c_k cos<k,u>` is real and `Z_P = {F = 0}`.
(2) `Phi(u) = sum_k c_k <k,u>^8 = 8! * g(u)` where `g` has 196 terms with coefficients in `{-2,-1,1,3}` and is
irreducible over `Q` (flint). (3) `g(u_+) = 19849334540098159992214034316 > 0` at
`u_+ = (1020, -1592, -6221, -2363, 5574, 2741, -3682)` and `g(-1,-1,-1,-1,0,0,0) = -1 < 0` (exact integers).
(`g` is negative at all 65604 nonzero-valued points of `[-2,2]^7`; its positive region is thin.)

*Proof.* Since `c` annihilates polynomials of degree `< 8`, `F(tu) = t^8 Phi(u)/8! + R(t,u)` with
`|R(t,u)| <= t^{10} sum_k |c_k| <k,u>^{10}/10!` (alternating cosine series). Choose balls `B_+`, `B_-` around
`u_+/|u_+|`, `u_-/|u_-|` on which `Phi > eta`, resp. `< -eta`, for some `eta > 0`. For `t` small, `F > 0` on
`t B_+` and `F < 0` on `t B_-`. Put `e = u_- - u_+` (normalised points) and let `B'` be a small 6-dimensional
ball in `e^perp` such that `u_+ + s in B_+` and `u_- + s in B_-` for `s in B'`. Each segment from `t(u_+ + s)` to
`t(u_- + s)` contains a zero of `F` (intermediate value theorem). Suppose every zero of `F` in the union `R_t`
of these segments were degenerate, i.e. lay on one of the countably many rational hyperplanes
`H_v = {<v,u> = 0}`. Let `pi` be the orthogonal projection onto `e^perp`; then `B'` is contained in
`union_v pi(Z(F) cap H_v cap R_t)`. Fix `v`. If `F|_{H_v}` is not identically zero, `Z(F) cap H_v` is the zero set
of a nonzero real-analytic function on the 6-dimensional space `H_v`, hence of 6-dimensional measure zero, and
its image under the linear map `pi|_{H_v}` has measure zero in `e^perp`. If `F|_{H_v} == 0`, then
`Phi|_{H_v} = lim_t 8! t^{-8} F(tu) == 0`, so the rational linear form `<v, u>` divides `Phi` in `Q[u]`,
contradicting the irreducibility of `g` (degree 8). So a set of positive measure would be a countable union of
null sets, which is impossible. Hence `R_t`, which lies within distance `O(t)` of `0`, contains a nondegenerate
zero of `F`. QED

**What this shows, and what it does not.**

* `Z_P` near `e` is here a real **hypersurface** through `e` (one real equation) that is not contained in any
  countable union of subtori: the "non-toric component through nondegenerate real points near `e`" of README
  §8.5 exists **in this model**, and in a profile that is the profile of an actual fake family. For this
  `(A, m)` the rigidity relation is a single real equation, and it cannot separate actual clusters from
  positions chosen on it.
* The model is faithful in the sense of §2 (`zeta < m`), but it is **not overdetermined**: `dim V = 1 < d = 7`.
  In the regime `dim V_H >= |A|/2`, which is exponentially larger than `d`. The two overdetermined Hadamard
  models (`dim V = 42` and `dim V = 7 = d`) are toric (Theorem 5.2). So the finite evidence says: ZE is not a
  consequence of faithfulness; it is an unlikely-intersection phenomenon that needs the codimension, exactly as
  `determinant.md` §5.5 describes it, and it does hold in both overdetermined models that could be computed.
* Whether a gap-admissible fake can be placed on this hypersurface cannot be decided at `n = 8`, `T = 1`: the
  fake constructions need `W` far beyond what an 8-point model encodes. The positive-measure good sets of the
  constructions (`round5/exact.md` Thm 5.1) would have to meet a hypersurface through `e`, which is plausible but
  not shown.

## 7. Adversarial search for gap-admissible fakes on Z_P

Fakes have nondegenerate positions (`referee_determinant.md` §4: full-rank profile plus the gap lemma).

* **Hadamard models `M8`, `M8'`.** By Theorem 5.2 every point of `Z_P` (indeed of `Bs` in the complex torus)
  is degenerate. So **no** nondegenerate position, in particular no fake position, lies on `Z_P`, anywhere on
  the torus. The adversarial search is settled exactly; no search is needed. As a sanity check of the
  generic statement (Cor. 4.3 of `determinant.md`), random positions `delta ~ U[0, 10^-3]^8` were certified off
  `Z_P` by "leading Taylor term > remainder bound" at 50 digits (`fakes_m8.py`): 100/100 for `M8`, 80/100 for
  `M8'` (the other 20 are not certified by this crude test at this radius; by Theorem 5.2 they are off `Z_P`
  anyway, being nondegenerate almost surely).
* **Skew double core `q = 7`, `T = 1`, `L = 2`, `m = 8`.** `Z_P` contains nondegenerate real points
  arbitrarily close to `e` (Proposition 6.2), on a real hypersurface. An adversary can therefore choose
  nondegenerate positions on `Z_P` for this `(A, m)`. Whether such positions can also be made gap-admissible
  (realisable by an angle model with the required gaps) is not decidable at this size (see §6bis).
* **Unfaithful models (`n - 1 <= 4`).** Every point of `Z_P` is degenerate (Prop. 2.1), so no fake lies on `Z_P`;
  by Prop. 2.2 this is uninformative for the regime.


## 8. Verdict for the route

**No growth improvement and no proof of ZE.** The finite test decides less than the README hoped, for a
reason that is itself a result, and then gives the first meaningful evidence.

1. **The proposed model sizes cannot test ZE.** At `n - 1 = 2` no admissible model is faithful (Prop. 2.4,
   proved); at `n - 1 = 3` none is, over all full-rank sign profiles with `T <= 4`, `L <= 3` (exact, 355 sets);
   at `n - 1 = 4` none of the 74 distinct sets examined is (§3, §3.2). In all of them ZE holds
   because of a zonotope (product of finite differences), i.e. a character, and Proposition 2.2 / 6.1 show that
   such certificates never exist for actual nondegenerate clusters in the regime. So a positive outcome at
   these sizes would have been evidence of nothing.
2. **Faithful models exist from `n = 8`.** For Hadamard profiles `zeta <= ((n-1)/n) max h_A(beta_p)`
   (Prop. 3.1), so every model satisfying the README's condition (ii) is automatically faithful; at `n = 8`,
   `T = 1`, the windows `m = 8` (`L = 2`) and `m = 11` (`L = 3`) are faithful (with (ii') in place of (ii)).
3. **In both overdetermined faithful models ZE holds, exactly and globally** (Theorem 5.2): the base locus in
   the whole complex torus lies in explicit subtori — collisions `delta_x = delta_y`, plus six four-term
   relations in `M8`. The mechanism is new: no single auxiliary form certifies anything (each has non-toric
   zero sets through nondegenerate real points), but the non-binomial factors are forms in `w - 1` and their
   common zeros lie on the braid arrangement. This is exactly the unlikely-intersection behaviour that ZE
   postulates, observed where characters are unavailable. **This is evidence for ZE** in the sense of README
   §8.5, first bullet, and it suggests the sharper algebraic target ZE-braid (§5.3).
4. **ZE-real fails in a codimension-one faithful model** (Prop. 6.2): for the skew double core at `q = 7`
   (`T = 1`, `L = 2`, `m = 8`, `dim V = 1`), `Z_P` contains a real hypersurface through `e` with nondegenerate
   points. So ZE is **not** a consequence of faithfulness (absence of characters) alone; it depends on the
   linear system being large. In the regime `dim V_H >= |A|/2`, exponentially larger than `d`, which is the
   overdetermined side; both overdetermined models computed are toric. The fake-transfer scenario of README
   §8.5 (second bullet) is therefore realised only in the underdetermined model, where it is not surprising.
5. **Status of route 5.** Not refuted, not proved. The evidence is consistent with ZE being true in the regime
   and identifies what must be proved: an overdetermination statement for the realisable symbols (ZE-braid),
   not a statement about single forms. The README's suggested next step (Seshadri-type bounds
   `L_A . Y >= m mult_e(Y)` for non-subtorus curves `Y` through `e`) gives only a *necessary* condition for `Y`
   to escape `Bs`; the models show that the gap to sufficiency is real: every single section `P_c` of a faithful
   model vanishes on non-toric hypersurfaces through `e`, so a proof must exploit many sections at once (as the
   saturation computation does), not a degree count for one section.

**Recommended next step.** Decide ZE-braid for the Hadamard slab-boxes in a family (`n = 8, 12, 16`, `T = 1`,
`L` with `|A| ~ |B_T|/2`): compute the realisable symbol space `Sigma_A` directly (it is a linear-algebra
problem on Newton coefficients, no exponentials) and test whether its common zeros lie on the braid
arrangement; then look for the structural reason (in `M8'` the symbols are cubics with 6 terms and quartics
with 12 terms in the `x_x`, attached to the Fano-plane structure of `H_8`).

## Appendix: scripts (`compute_checks/`; each writes `*_out.txt` or the named output)

| script | content | status |
| --- | --- | --- |
| `core.py` | slab-boxes, moment matrices, saturated integer kernel (HNF of `[M^T | I]`, LLL), tangent-cone forms | exact |
| `explore.py`, `n3_hexagon.py`, `bs2d.py` | `n = 3`: admissible models; hexagon base loci (gcd, factorisation, resultants) | exact (+ certified acb root balls) |
| `scan.py` | screening of full-rank sign profiles, `rho` mod `2^61 - 1` | exact mod p (confirmed over Z later) |
| `models.py`, `zono.py`, `faithful.py`, `lpzeta.py` | exact `rho`; exhaustive zonotope search; LP dual bound on `zeta` | exact |
| `gap_search.py`, `bigger.py` | random symmetric polytopes `d = 3, 4`; `n = 6, 8` design profiles | exact |
| `had8.py`, `model8.py` | `n = 8` models: exact `rho`, exact basis of `V(A,m)`, `V_H = V` check | exact |
| `had8_tc.py` | exact tangent-cone forms of `M8` and coordinate subspaces in `TC` | exact |
| `cone_reduce.py`, `cone_sat.py` | Lemma 5.1 / Theorem 5.2: factorisation, form test, Bayer saturations | exact (flint, sympy over Q) |
| `factor_basis.py`, `nonbinom.py`, `m8p_reduce.py`, `cone_variety.py` | factor structure, Groebner basis and Hilbert function of the `M8'` cone (degree 5 threefold) | exact |
| `sdc7_refute.py`, `sdc7_tc.py` | Prop. 6.2 inputs: parity, exact sign witnesses, irreducibility of `Phi` | exact (the float minimisation inside `sdc7_tc.py` only locates the positive witness, which is then evaluated exactly) |
| `fakes_m8.py` | random positions certified off `Z_P` (50 digits) | certified |
| `zetest.py`, `tc_search.py`, `bs_search.py`, `bs_search2.py`, `explain.py` | numerical searches for zeros of `TC` and `Bs` | evidence only; superseded by Theorem 5.2 |
| `sdc7_L3_zeta.py`, `n5_hard.py`, `n5_hard_basis.py` | single hard instances | exact |
