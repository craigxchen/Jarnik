## 3. Part (b): structure that could let the reach exceed `zeta`

### 3.1 Nested shift chains (proved) and what they dominate

**Theorem 3.1 (chain criterion; proved).** Let `L` be a set of shifts and `S_0 subset S_1 subset ... subset S_m`
finite sets with `S_0 != {}` and, for each `e < m`,

```text
(CH_e)    V(S_e, e)  subset  sum_{l in L} V((S_{e+1})_l, e)  +  V(S_{e+1}, e+1).
```

Then every `Y in G_m^d(C)` with `Y^l != 1` (`l in L`) has `P_c(Y) != 0` for some `c in V(S_m, m)`; if `S_m subset A`,
`Y notin Bs V(A, m)`. With rational `(CH_e)` and the analogue of (3.3) at the top step, the height version holds as in
round 27 Thm 3.1(3).
*Proof.* Induction on `e`: `c_0 = delta_k`, `k in S_0`. Given `c_e in V(S_e, e)` with `P_{c_e}(Y) != 0`, split it by
(CH_e) as `sum_l c'_l + r`; either `P_r(Y) != 0` (take `c_{e+1} = r`) or some `P_{c'_l}(Y) != 0`, and then
`c_{e+1} = T_l c'_l - c'_l in V(S_{e+1}, e+1)` (round 27 Lemma 3.0 for the set `S_{e+1}`) has value
`(Y^l - 1) P_{c'_l}(Y) != 0`. QED

Special cases: the constant chain `S_e = A` is round 27's Thm 3.1; the constant chain `S_e = A'` on a **subset**
`A' subset A` is the subset criterion of §2.4; the **zonotope chain** `S_e = k_0 + {sum_{i in T} v_i : T subset [e]}`
with `v_i in L` satisfies (CH_e) trivially, because `S_e subset (S_{e+1})_{v_{e+1}}`. Hence, with
`reach_chain(A, L)` the largest `m` reachable by a chain inside `A`,

```text
max( zeta_L(A),  max_{A' subset A} reach(A', L) )  <=  reach_chain(A, L)  <=  truth_L(A)  <=  rho(A).
```

So the chain criterion **dominates zonotopes by construction**; "reach versus zeta" is then the question whether the
inequality on the left is strict on regime boxes. Mixed chains interpolate: for a product sub-box
`Z x B subset A` (`Z` a zonotope in some coordinates, `B` a box in the others), Theorem 2.1 gives
`reach(Z x B) = ord(Z) + reach(B)` when `Z` is a cube-type zonotope (generators `e_x`), so zonotope steps can be traded
for free shift steps; each zonotope step costs `Lambda(v) ~ 1` of the Liouville budget, shift steps cost nothing, as
round 27 open point 3 anticipated. On M8 the subset `A cap {|sigma| <= 2}` reaches 7 `= zeta_pair` (§3.5); I found no
chain on M8 beating 7 (truth `<= 8`).

### 3.2 The derivative (first-order) structure of `J(A_l)` (proved; data mod p)

**Lemma 3.2 (proved).** For finite `S, S'`: `I(S cap S') = I(S) + I(S')` (as `K^S (x)_R K^{S'} = K^{S cap S'}`). Hence
`I(A_l) = I(A) + tau_l I(A)` with `(tau_l f)(k) = f(k + l)`, and for `f in I(A)` with leading form `F`,
`tau_l f - f in I(A_l)` has leading form `d_l F` (if nonzero). Therefore

```text
J1(A, l) := J(A) + R . d_l J(A)  subset  J(A_l)      (d_l = directional derivative; k_x^2 = d_{e_x}(k_x^3)/3 is included).
```

For slab boxes `I(A) = I(C) + (omega_p)`, `omega_p = prod_{|t| <= L_p} (s_p - t)`, and
`I(A_l) = I(C_l) + (omega_{p,l})` with `deg omega_{p,l} = 2L_p + 1 - |s_p(l)|` (gcd of two univariate slab polynomials).

*Data* (`t3_firstorder.py`, M8, all `e <= 8`, `l = e_0, e_0 - e_1`): `J1(A,l)_e = J(A_l)_e` in every degree except one
dimension at `(l, e) = (e_0, 6)` and `(e_0 - e_1, 5)` (second-order Taylor terms, which depend on the sub-leading
coefficients of `I(A)`, not on `J(A)` alone). The M8 obstruction (§4) already lies in `cap_l J1(A,l)_6` (`t4`).
**Consequence.** `(SHIFT_e)` implies the first-order property `cap_l J1(A,l)_e = J(A)_e`, a statement about `J(A)`
alone, and by §3.3 about the tail ideal `J(B)` alone. Failures of the first-order property are failures of `(SHIFT)`.

### 3.3 The linkage form (proved; verified mod p)

`C = {-1,0,1}^d` is a complete intersection whose leading forms `k_x^3` are a regular sequence; `R/J(C)` is Gorenstein
with socle degree `2d`. The identity `H_A(e) = T_d(e) - H_B(2d - e)` (round 13; verified here for all six model boxes,
`t1`) and the inclusion `J(A) subset J(C) : J(B)` (`LF(f) LF(g) = LF(fg)`) give, by Gorenstein duality,
`J(A) = J(C) : J(B)`. The same holds for `A_l` inside the grid `C_l = C cap (C - l)` (socle degree `2d - |supp l|`;
verified for `l = e_0, e_0 - e_1` on M8, `t14`). The inclusion `(R_e/J(C_l)_e)^* -> (R_e/J(C)_e)^*` is, under the two
dualities, multiplication by `mu_l = k_x` (`l = +-e_x`) or `k_x k_y` (`l = +-(e_x - e_y)`). Hence:

**Proposition 3.3 (linkage form of the criterion; proved given the linkage identity).** With `B_l = C_l \ A_l`,

```text
(SHIFT_e)   <=>   J(B)_{2d-e}  =  sum_{l in L_pair} mu_l . J(B_l)_{2d-e-deg mu_l}      (mod J(C)),
```

where `B_l = C_l cap (B u (B - l))` is the shifted tail inside the grid `C_l`, and `dim J(B)_{2d-e}/J(C) = H_A(e)`.
Verified on M8 (`t14`, mod p): `e = 5`: 259 = 259; `e = 6`: 286 vs 287; `e = 7`: 141 vs 175, matching the direct
computation exactly. So `(SHIFT_e)` says: *every top-degree leading form of the tail ideal is a combination of `k_x`
(resp. `k_x k_y`) times top-degree leading forms of the shifted tails.* The tail is an explicit union of `r` slab tails
`{|s_p| > L}`, each a union of products of lower sets; in the regime the relevant degrees `2d - e`, `e <= 8j`, are
within `8j = o(d)` of the socle. I did not find an argument that decides this form either way; it is the natural
setting for a proof by induction on the number of slabs, but the intersection over `p` of slab-tail ideals is exactly
where the non-univariate obstructions of §3.4 live.

### 3.4 What the failures are: implied-slab powers (data, mod p)

For each failure `(A, e)` I tested all `phi^e`, `phi in {-1,0,1}^d` (`t6_powers.py`, `t9`):

| box, degree | `delta(e)` | `phi^e` that are obstructions | dim of their span mod `J(A)_e` | typical `phi`, range on `A` -> max on `A_l` |
| --- | --- | --- | --- | --- |
| M8, `e = 6` | 1 | 1 | **1** | `sigma = (1,...,1)`: 7 -> 6 values (rigid Fano level) |
| M8 minus `|sigma| = 3`, `e = 7` | 28 | 70 | **28** | many, e.g. `(1,-1,-1,-1,-1,-1,1)` |
| sdc7 `L = 2`, `e = 6` | 22 | 14 | 7 | `sigma - k_x`: 7 -> 6 values |
| random `n=6` seed 101, `e = 6` | 1 | 1 | **1** | `(0,1,-1,1,-1)`: 9 -> 8 |
| random `n=6` seed 100, `e = 4` | 2 | 7 | **2** | `(1,0,0,0,-1)`: 5 -> 4 |
| random `n=7` seed 4, `e = 4` | 1 | 5 | **1** | `(1,1,1,0,0,1)`: 5 -> 4 |
| random `n=7` seed 1, `e = 7` | 9 | 12 | 6 | ranges 9–11 |
| random `n=6` seed 113, `e = 6` | 1 | 8 | **1** | `(1,1,1,0,1)`: 7 -> 6 |

In 6 of 8 cases the whole defect is spanned by powers of linear forms with few values on `A` ("implied slabs"); several
are explained exactly by Theorem 2.3 (range drops by one on every `A_l`), others lie in `cap_l J(A_l)` without a range
drop. The residue (sdc7: 15 of 22 dimensions; random `n = 7` seed 1: 3 of 9) is not spanned by powers. **The failures
of the criterion on the computable boxes are therefore mostly univariate, of the implied-slab type, and tied to
level sets of small-range forms that are rigid or nearly rigid.** In the regime (`e <= 8j = o(d)`), small-range
integer forms are: the slab forms `s_p` (range `2j + 1`; their extreme levels are far from rigid, since moves outside
`P_p` are free), sums of a few slab forms (same), and the base coordinate and its near relatives `sigma - k(X)`, `|X|`
bounded (range `~4j`, rigid only by design coincidence, §2.3).

### 3.5 Why M8' reaches 11 (and M8 does not)

* **M8 (`L = 2`).** The top level of `sigma` is `{1_l : l a Fano line}` (7 points, `h = 3 = floor(7L/4)`, i.e. the
  averaging bound `sum_p s_p = 4 sigma` is attained). It is rigid: replacing one point of a line gives a 3-set lying in
  some pattern (a Fano-line complement), whose slab value 3 exceeds `L`. So `sigma^6` blocks `(SHIFT_6)` (§4).
* **M8' (`L = 3`).** `h(sigma) = 4 < floor(7L/4) = 5`: the averaging bound is not attained, the 28 top points carry
  slack, and 252 pair moves stay in `A` (`t7`). The BCO is absent (the range argument for
  `sigma^8 in J(A_{e_0 - e_1})` no longer applies, and `(SHIFT_8)` holds). Empirically no implied-slab obstruction
  occurs below `rho`, and the criterion runs to `rho = 11` (mod p). The extremes of
  `sigma` in M8' are fixed by the cube geometry of the vertex cuts (`sigma = 5` would force some pattern to take its
  vertex value `4 = |P|`), not by the slab LP (`sum_p s_p = 4 sigma <= 21` allows 5): they are *cube-bound*, not
  *slab-bound*. sdc7 `L = 3` behaves the same way (`h = 4`, not rigid, reach `= rho = 10 = zeta_pair`).
* **What this means for the regime.** M8' works because its cuts are algebraically weak (`L = |P| - 1`: each slab only
  removes the two vertices of a 4-cube) and no design-type extremal configuration exists at `L = 3`. Regime boxes have
  algebraically strong cuts (`j/|P| -> 0`), which makes them slab-bound like M8; but their `L` is generic, which (by
  §2.3) makes them like M8' as far as rigidity is concerned. The two features that separate M8 from M8' pull in
  opposite directions at regime size, so **M8' is not evidence for the regime, and M8 is not evidence against it.**
* **No deeper "reason" for M8' was found.** The staircase (standard-monomial) sufficient criterion
  (`N(A)_e subset Union_l N(A_l)_e` for some degree-compatible order implies `(SHIFT_e)`; proof: initial ideals of an
  intersection lie in the intersection of initial ideals) certifies `(SHIFT_e)` on M8' only for `e <= 7` (deglex and
  degrevlex; `t13`), and on M8 only for `e <= 4`. The reach 11 on M8' is a genuine linear-algebra phenomenon beyond the
  combinatorial (Gröbner) certificates I could produce.
