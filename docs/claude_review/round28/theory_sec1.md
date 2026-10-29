## 1. Definitions and conventions (checked against the round-27 files)

Labels as in the record: **proved** (proof given here, unrefereed), **exact** (finite computation over `Q` or by
exhaustive search), **mod p** (ranks modulo `P1 = 2147483629`, and where stated also `P2 = 1000000007`; a mod-`p`
equality of ranks is evidence, not proof), **heuristic**, **open**. Scripts: `theory_checks/` (absolute paths in §7).

**Boxes.** `n` points, base point `n` (1-indexed, as in the record; in the scripts the points are `0..n-1` and the base
is index `n-1`, so syl8's base coordinate `k_8` is `k[7]` there), `d = n - 1`, frequencies `k' in {-1,0,1}^d`,
`k_n = -sum k'` (round-27 killer's `klib.build_box`; equivalent to `interp_checks/models.py` up to complementing
patterns). A *slab box* is

```text
A = { k' in {-1,0,1}^d : |k(S_p)| <= L_p  for every slab (S_p, L_p) },     k(S) = sum_{x in S} k_x  (S a set of points).
```

A slab whose set contains the base reads `|k(S)| = |k(S^c)|` (as `sum_all k_x = 0`), so every slab is a cut on a subset
`P_p` of the `d` cube coordinates ("pattern"). In the regime (round 12 `t1rigour.md`), sign columns give
`<k,beta_p> = 2 k(P_p)` with `P_p` the side omitting the base, and `L_p = j = floor(w/2)`, `w = sqrt(8 d ln 4r)`. `C = {-1,0,1}^d` is the ternary cube, `B = C \ A`
the *tail*. Models: `syl8` (7 Fano-complement patterns of size 4, `n = 8`), `sdc7` (14 patterns), `L in {1,2,3}`;
`M8 = syl8, L = 2`; `M8' = syl8, L = 3` (round 11 `compute.md` §5). Every slab box is centrally symmetric, `A = -A`.

**Least spaces, leading forms.** `V(S,e)` = measures on `S` with all moments of order `< e` zero; `M_e(S)` = their
degree-`e` moment vectors; `J(S)_e` = degree-`e` forms `F` with `F|_S in Poly_{<e}|_S` (the degree-`e` part of the
leading-form ideal `J(S)` of `I(S)`); `M_e(S) = J(S)_e^perp`; `H_S(e) = dim M_e(S)`; `rho(S) = max{e : H_S(e) > 0}`.
On subsets of `C` we compute in standard monomials (exponents `<= 2`), since `k_x^3 - k_x in I(C)`.

**Shift property, reach.** `S_l = S cap (S - l)`. `(SHIFT_e)` w.r.t. `L`: `M_e(S) = sum_{l in L} M_e(S_l)`, equivalently
`J(S)_e = cap_{l in L} J(S_l)_e`. `reach(S, L) = max{m : (SHIFT_e) holds for all e < m}`. `L_pair = {+-e_x, +-(e_x - e_y)}`
(all `x != y` in `[d]`; in full coordinates these are the roots `e_x - e_y`, `x != y in [n]`). Because `S = -S` implies
`S_{-l} = -(S_l)` and `M_e(-T) = M_e(T)` (sign `(-1)^e`), only one of `+-l` is needed for centrally symmetric `S`.

**The orders compared.**
* `truth_L(S) = max{m : Bs V(S,m) subset Union_{l in L} {Y^l = 1}}` (ZE-braid order for `L = L_pair`).
* `zeta_L(S)` = max order `t` of a lattice zonotope `{k_0 + sum_{i in T} v_i : T subset [t]} subset S` with all `v_i in L`
  (repetitions allowed). This is the definition in the task statement. **Correction to be aware of:** the round-27
  killer's DFS (`k6_zono.py`) used *all* primitive generators in `{-2..2}^d`; call that `zeta_all >= zeta_pair`. The M8'
  order-10 zonotope of the record uses `e_0 - e_1 - e_2 notin L_pair`. Exhaustive DFS restricted to `L_pair`
  (`t11_zeta_pair.py`, exact search): `zeta_pair = 3, 7, 10, 3, 7, 10` for syl8 `L = 1, 2, 3` and sdc7
  `L = 1, 2, 3`. So `zeta_pair(M8') = 10 = zeta_all(M8')`; nothing in round 27 changes.
* `rho(S)`: the regularity / LP bound.

**Elementary inequalities (proved).** `reach(S,L) <= truth_L(S) <= rho(S)` and `zeta_L(S) <= truth_L(S)`.
*Proof.* Thm 3.1 of `interp.md` gives `Y notin Bs V(S,m)` for `m <= reach` and `Y^l != 1` (`l in L`); `V(S, rho+1) = 0`
has base locus everything, so some `(SHIFT_e)` with `e <= rho` fails and `reach <= rho`. A zonotope form
`Y^{k_0} prod (Y^{v_i} - 1)` with `v_i in L` lies in `V(S,t)` and does not vanish off `Union {Y^l = 1}`. QED

**Recomputed reference data (mod p unless marked).** `t1_profiles.py`, defect `delta(e) = H_A(e) - dim sum_l M_e(A_l)`
w.r.t. `L_pair`, all `e <= rho`:

| box | `|A|` | `H_A` | `delta(e)`, `e = 0..rho` | reach | `zeta_pair` (exact) | rho |
| --- | --- | --- | --- | --- | --- | --- |
| syl8 `L = 1` | 57 | 1 7 28 21 | 0 0 0 21 | 3 | 3 | 3 |
| syl8 `L = 2` (M8) | 1037 | 1 7 28 77 161 259 287 175 42 | 0 0 0 0 0 0 **1** 34 42 | 6 (**exact** failure at 6, §4) | 7 | 8 |
| syl8 `L = 3` (M8') | 1877 | 1 7 28 77 161 266 357 386 329 195 63 7 | 0 (x11) 7 | 11 | 10 | 11 |
| sdc7 `L = 1` | 57 | 1 7 28 21 | 0 0 0 21 | 3 | 3 | 3 |
| sdc7 `L = 2` | 757 | 1 7 28 77 161 251 188 43 1 | 0 0 0 0 0 0 **22** 43 1 | 6 | 7 | 8 |
| sdc7 `L = 3` | 1765 | 1 7 28 77 161 266 357 379 301 153 35 | 0 (x10) 35 | 10 | 10 | 10 |

These reproduce the round-27 numbers (killer: M8' reach 11 at two primes; syl8/sdc7 `L = 2` reach 6) and add sdc7
`L = 3` (reach = rho = zeta_pair = 10) and the full defect profiles. The linkage identity
`H_A(e) = T_d(e) - H_B(2d - e)` holds in all six boxes (mod p; `T_d(e)` = coefficient of `t^e` in `(1+t+t^2)^d`).
