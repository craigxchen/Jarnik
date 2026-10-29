# Round 28 — theory: multi-shift reach versus the zonotope order

Status: complete (unrefereed). Lane: theory. Date: 2026-10-29. Companion computation lane: `compute.md`.

## 0. Verdict

**No growth improvement and no uniform bound.** `M <= (2/3+eps) log R/loglog R` stands. The comparison of the
multi-shift reach with the zonotope order is **not decided at regime size**; what is proved is below.

1. **(a) Upper bounds.**
   * *No unconditional bound (proved, Thm 2.1, Cor 2.2).* Reach, `rho`, `zeta_pair`, `zeta_all` are all additive under
     products of boxes, so `(M8')^t` has `reach = 11t`, `zeta = 10t` (M8' input mod p, two primes). No inequality
     `reach <= zeta + c` holds on slab boxes. Product (block) profiles violate the character gap, so this does not help
     the route, but it shows that any such bound must use the gap, and that the favourable branch of the decision rule
     can be met vacuously.
   * *A conditional Liouville-capped bound (proved, Thm 2.3, Cor 2.4, Prop 2.5).* **Implied-slab obstruction:** if an
     integer form `phi` takes `R` values on `A` and at most `R - 1` on every `A_l`, then `phi^{R-1} in cap_l J(A_l)`; if
     also `phi^{R-1} notin J(A)`, the criterion stops at `R - 1`. For the **base coordinate** `sigma = -k_n` of a
     balanced sign-column profile, `|sigma| <= 2j` on the `T = 1` box (the averaged slab has the zonotope budget
     `Lambda(A) <= 2w`), so if its top level is rigid under pair moves (and `sigma^{2h} notin J(A)`), `reach(A, L_pair) <= 4j <= 2w < m' - 1`.
   * *But the hypothesis is a design coincidence.* Rigidity holds in M8 (the top level is the 7 Fano lines) and, in
     Sylvester boxes up to `n = 128` (IP), exactly when `h = 2L - 1`, i.e. for `L` a power of 2 (subspace configurations);
     it fails in every decided non-design case, including all decided random profiles. So (a) is **open** in the regime; the exact obstruction to a proof is
     that no other low-degree obstruction family is known (slab-power products cost `>= 2j log_2 n`; non-univariate
     residues exist in models but have no known structure).
2. **(b) Structure.**
   * *The syl8 `e = 6` obstruction is `sigma^6`, the sixth power of the base-point coordinate `k_8` (exact, Thm 4.1).*
     `|sigma| <= 3` because `sum_p s_p = 4 sigma` (Fano incidences); the top level `{1_l : l a Fano line}` is rigid under
     weight-2 moves and not under weight-4 moves, which is why all cube shifts rescue `(SHIFT_6)`.
   * *Why M8' reaches 11:* its base-coordinate top level is not rigid (`h = 4 < floor(7L/4) = 5`; 252 legal pair moves),
     because `L = 3` admits no design and its cuts are cube-bound vertex cuts; no implied-slab obstruction arises below
     `rho`. M8 and M8' differ in exactly the two features (design value of `L`; slab- versus cube-bound extremes) that
     pull in opposite directions at regime size, so neither is evidence about the regime.
   * *Most model failures are implied-slab powers* (6 of 8 defects fully spanned by `phi^e`, `phi` ternary).
   * *Chains and subsets (proved, Thm 3.1):* any chain or subset `A' subset A` gives the same conclusion, so
     `reach_chain >= max(zeta_pair, reach)`; on M8, `A cap {|sigma| <= 2}` lifts the reach from 6 to 7 `= zeta_pair`.
   * *First-order and linkage forms (proved; verified mod p):* `J(A_l) supset J(A) + R d_l J(A)` (equal in models up to
     one dimension), and `(SHIFT_e) <=> J(B)_{2d-e} = sum_l mu_l J(B_l)_{2d-e-deg mu_l}` for the tail `B = C \ A`.
3. **(c) Regime.** `(SHIFT)_A` (or its chain version) at regime size requires the absence of rigid implied slabs of
   range `<= 8j` **and** the absence of non-univariate obstructions in degrees `(4j, 8j)`. Hilbert-function tools
   (Macaulay, Gotzmann, Green, gin, linkage) control dimensions of single leading-form ideals, not the position of
   `~d^2` subspaces, and representation counting gives no obstruction for actual (symmetry-free) profiles. **No known
   asymptotic tool decides it.**
4. **Decision rule.** The theory does not decide it. It advises reading failures that are design effects (rigid
   implied slabs at Sylvester `L = 2^m`) and successes on product or vertex-cut boxes as uninformative about the
   regime, and using `zeta_pair` as the same-hypothesis competitor (§6).

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

## 2. Part (a): upper bounds on the reach

### 2.1 No unconditional bound `reach <= zeta + c` on slab boxes (proved, modulo the mod-p M8' input)

**Theorem 2.1 (products; proved).** Let `A_1 subset Z^{d_1}`, `A_2 subset Z^{d_2}` be finite, `A = A_1 x A_2` (disjoint
coordinate blocks, one shared base point, so `L_pair(d_1 + d_2)` consists of the in-block roots and the cross roots
`+-(e_x - e_y)`, `x` in block 1, `y` in block 2). Then

```text
H_A = H_{A_1} * H_{A_2} (convolution),  rho(A) = rho(A_1) + rho(A_2),
reach(A, L_pair) = reach(A_1, L_pair) + reach(A_2, L_pair),
zeta_pair(A) = zeta_pair(A_1) + zeta_pair(A_2),   zeta_all(A) = zeta_all(A_1) + zeta_all(A_2).
```

*Proof.* (i) `K[A] = K[A_1] (x) K[A_2]`, and the degree filtration of `K[A]` is the tensor product of the degree
filtrations (a polynomial of degree `<= e` is a sum of products of polynomials of degrees `a`, `b` with `a + b <= e`).
Over a field, `gr` of a tensor product of filtered spaces is the tensor product of the `gr`'s, so
`R/J(A) = R_1/J(A_1) (x) R_2/J(A_2)` as graded algebras and, dually, `M_e(A) = (+)_{a+b=e} M_a(A_1) (x) M_b(A_2)`. This
gives the convolution and the additivity of `rho`.
(ii) For `l = (l_1, l_2)` one has `A_l = (A_1)_{l_1} x (A_2)_{l_2}` (with `S_0 = S`). Put
`N_a = sum_{l_1} M_a((A_1)_{l_1})` (in-block roots of block 1) and `N'_b` likewise. In-block shifts of block 1 contribute
`(+)_{a+b=e} N_a (x) M_b(A_2)`, those of block 2 contribute `(+) M_a(A_1) (x) N'_b`, and a cross root
`e_x - e_y = (e_x, -e_y)` contributes `M_a((A_1)_{e_x}) (x) M_b((A_2)_{-e_y}) subset N_a (x) N'_b`. Hence

```text
sum_l M_e(A_l) = (+)_{a+b=e} ( N_a (x) M_b(A_2) + M_a(A_1) (x) N'_b ),
```

and since `(M (x) M')/(N (x) M' + M (x) N') = (M/N) (x) (M'/N')`, `(SHIFT_e)(A)` holds iff for every `a + b = e` with
`H_{A_1}(a) H_{A_2}(b) > 0`, `(SHIFT_a)(A_1)` or `(SHIFT_b)(A_2)` holds. With `r_i = reach(A_i)` (finite, `<= rho`) the
first failure of `A_i` is at `r_i` with `H_{A_i}(r_i) > 0`; so `(SHIFT_e)(A)` holds for `e < r_1 + r_2` and fails at
`e = r_1 + r_2`.
(iii) `>=` for `zeta`: products of zonotopes (block-1 generators are roots of `L_pair(d)`). `<=`: project a zonotope
`{k_0 + sum_{i in T} v_i} subset A` to block 1; the image is the lattice zonotope with generators `{v_i^{(1)} != 0}`, all
of whose subset sums lie in `A_1`, and `v_i^{(1)} in L_pair(d_1)` when `v_i in L_pair(d)` (a cross root projects to
`+-e_x`). Every `i` has `v_i^{(1)} != 0` or `v_i^{(2)} != 0`, so `t <= zeta(A_1) + zeta(A_2)`. QED

**Corollary 2.2.** For the `t`-fold product `(M8')^t` (`d = 7t`, 7t slabs, `L = 3`):
`reach(L_pair) = 11t = rho`, `zeta_pair = zeta_all = 10t`, `rho/L = 11t/3`.
Inputs: `reach(M8') = 11` (**mod p**, two primes: killer `k13/k14`, reproduced here `t1_syl8_L3_out.txt`),
`zeta_pair(M8') = 10` (**exact** exhaustive search, `t11`), `zeta_all(M8') = 10` (record LP bound + killer's exact
order-10 zonotope). So **`reach - zeta` is unbounded on slab boxes**; no theorem `reach <= zeta + c` (or
`reach <= (1 + o(1)) zeta`) holds for slab boxes in general. (A small exactly-verified base box with
`reach > zeta_pair` would make the corollary unconditional; a search of ~400 small boxes, `n <= 6`, found none,
`t12_out.txt`.)

*Why this does not help the route.* Block profiles violate the character gap: in the regime normalisation a root
inside one block is seen only by that block's slabs, of total weight `~1/t`, so `Lambda(e_x - e_y) = O(1/t)`, whereas
every actual cluster has `Lambda(v) >= 1 - o(1)` (interp Lemma 2.1). A bound `reach <= zeta + c` can only hold on
gap-satisfying ("connected") profiles, and must use that hypothesis. *Consequence for the decision rule:* its
favourable branch ("reach > zeta on a growing share of boxes as they grow, with `rho/L` increasing") is satisfied
vacuously by product boxes; it must be read on gap-satisfying multi-cut profiles (Sylvester, sdc, random
non-block profiles), as the judge intended.

### 2.2 The implied-slab obstruction (proved) and the base coordinate

**Theorem 2.3 (implied-slab obstruction; proved).** Let `S subset Z^d` be finite, `L` a set of shifts,
`phi : Z^d -> Z` linear, `V = phi(S)`, `R = |V|`. If `|phi(S_l)| <= R - 1` for every `l in L`, then
`phi^{R-1} in cap_{l in L} J(S_l)_{R-1}`. Hence if moreover

```text
(N_phi)    phi^{R-1} notin J(S)_{R-1},
```

`(SHIFT_{R-1})` fails and `reach(S, L) <= R - 1`.

*Proof.* `omega_l(t) = prod_{v in phi(S_l)} (t - v)` has degree `R_l <= R - 1` and `omega_l(phi) in I(S_l)`; its leading
form is `phi^{R_l}`, so `phi^{R-1} = phi^{R-1-R_l} phi^{R_l} in J(S_l)`. QED

*The hypothesis, shift by shift.* If `phi(l) != 0` then `phi(S_l) subset V cap (V - phi(l))`, which misses `max V` or
`min V`: automatic. If `phi(l) = 0` it is a **rigidity** condition: some level set `{phi = v}` of `S` has no point
`k` with `k + l in S`.

*Certificates for `(N_phi)` (proved).* `(N_phi)` holds iff some `c in V(S, R-1)` has `sum_k c_k phi(k)^{R-1} != 0`
(equivalently, the push-forward `phi_* : V(S,R-1) -> {measures on V}` is nonzero; its image is spanned by the
alternating binomial measure on the `R` values). A sufficient certificate is a lattice zonotope of order `R - 1` in
`S` all of whose generators have `phi(v_i) != 0`: for `c = prod (T_{v_i} - 1) delta_{k_0}`,
`sum c_k phi(k)^{R-1} = (R-1)! prod_i phi(v_i)`.

**Corollary 2.4 (base-coordinate obstruction, BCO; proved).** Let `S subset {-1,0,1}^d`, `sigma = k_0 + ... + k_{d-1}
= -k_n` (the base coordinate), `h = max_S sigma`, `R = |sigma(S)|`. Suppose

```text
(R)  no k in S with sigma(k) = h has k + e_x - e_y in S  (x != y),           (N)  sigma^{R-1} notin J(S).
```

Then `reach(S, L_pair) <= R - 1 <= 2h` (`R = 2h + 1` when `S = -S` and `sigma(S)` is an interval).
*Proof.* `sigma(+-e_x) = +-1 != 0`; for `l = e_x - e_y` (all ordered pairs), (R) removes the level `h` from
`sigma(S_l)`. Apply Theorem 2.3. QED

**Proposition 2.5 (the BCO is Liouville-capped; proved).** For a sign-column profile with weights `omega_p`
(`sum omega_p <= 1`) and sides `P_p` omitting the base, let `f_x = sum_p omega_p [x in P_p]`. On the `T = 1` box,
`|sum_{x < n} f_x k_x| <= j`. If the profile is **balanced at the base**, `f_x = f` for all `x != n`, then
`h(sigma) <= j/f`, and for `f = 1/2`: `2h <= 4j <= 2w`. Consequently, for a balanced profile satisfying (R) and (N),
`reach(A, L_pair) <= 2w < ceil(4w) - 1`: `(SHIFT)_A` fails, i.e. the (constant-chain) hypothesis of round 27's
Corollary 4.1 is false for that profile.
*Proof.* `sum_p omega_p <k, beta_p> = sum_p omega_p 2k(P_p) = 2 sum_x f_x k_x` and `|k(P_p)| <= j`. QED

So the averaged slab `<k, sum_p omega_p beta_p>`, whose range is bounded by the same budget
`Lambda(A) = sum_p omega_p h_A(beta_p) <= 2w` that caps zonotope orders, is a multiple of the integer form `sigma`
exactly when the profile is balanced, and then it caps the constant-chain criterion at the zonotope scale **whenever
its top level is rigid** (and (N) holds).
This is the precise "Liouville-capped" upper bound that part (a) asks for, with its exact hypothesis.

### 2.3 When is the top level rigid? (exact on models; IP evidence up to `n = 128`)

* Models (exhaustive, `t7_bco_scan.py`): (R) holds for syl8 `L = 2` only (top level = the 7 Fano lines `1_l`,
  `h = 3`; §4), and for 3 of the 27 random boxes (in each, `reach = 2h` exactly as Cor 2.4 predicts). (R) fails for
  syl8 `L = 1, 3`, sdc7 `L = 1, 2, 3` and the other 24 random boxes.
* Sylvester `n = 2^k`, `T = 1`, IP with HiGHS (`t10_bco_ip.py`; one pair suffices by 2-transitivity of the base
  stabiliser):

| `n` | `L` | `h = max sigma` | averaging bound `floor(2L(1-1/n))` | `max sigma` on `A_{e_0 - e_1}` | (R) |
| --- | --- | --- | --- | --- | --- |
| 8 | 1, 2, 3 | 1, **3**, 4 | 1, 3, 5 | 1, 2, 4 | no, **yes**, no |
| 16 | 2, 3, 4 | **3**, 4, **7** | 3, 5, 7 | 2, 4, 6 | **yes**, no, **yes** |
| 32 | 3, 5, 8 | 4, 8, **15** | 5, 9, 15 | 4, 8, 14 | no, no, **yes** |
| 64 | 8, 12, 16, 26 | **15**, 21, **31**, 48 or 49 | 15, 23, 31, 51 | 14, 21, 30, 47 or 48 | **yes**, no, **yes**, undecided |
| 128 | 16, 24 | **31**, 45 or 46 | 31, 47 | 30, — | **yes**, undecided |

  (`L = 26` at `n = 64` is the regime value `j(64)`; `n = 64, L = 12` needed a 1500 s IP. Undecided entries: the IP
  did not close in 400–900 s.) In every decided case (R) holds exactly when `h` attains the averaging bound `2L - 1`,
  which in these data happens exactly for `L` a power of 2. **Explanation (proved for the configurations, not for all
  maximisers):** with the base at `0 in F_2^k` and `V` a linear subspace of dimension `m`, `S = V \ {0}` has `|S cap P_u| = 2^{m-1}` or `0`, so `1_S` attains
  `sigma = 2L - 1` when `L = 2^{m-1}`; at `sigma = 2L - 1` the total slack `sum_u (L - k(P_u))` equals `n/2 - L`, and a
  pair move needs slack `>= 1` on the `n/4` patterns separating the pair, so for `L > n/4` (R) is automatic; for `m >= 2`
  the subspace configurations themselves admit no pair move (the affine set `{u.x = 1, u.y = 0}` is not contained in the
  linear space `V^perp`). The IP confirms that *all* maximisers are rigid in the listed cases.
* Random (non-symmetric, only approximately balanced) sign-column profiles, `r = n - 1` random half-size patterns
  (`t16_random_rigidity.py`, one IP for `h` and one feasibility IP for "some top point has a legal pair move"): here
  `h` can exceed `2L` (e.g. `n = 64`: `L = 6, 10, 16, 26` give `h = 14, 21, 33, 51`). (R) **fails** for `n = 16, L = 2`
  (3 seeds), `n = 32, L = 4` (2 seeds), `n = 64, L = 10` and `L = 26`; the feasibility IP did not close for
  `n = 64, L = 6, 16` and `n = 64, r = 126, L = 10` (300 s), and for `n = 128` (`L = 12, 20, 39`) already the IP for
  `h` did not close (`h >= 25, 41, 78`, bounds 27, 42, 79). No random profile was found rigid.
* **Reading.** A rigid top level is a *design coincidence* (Fano lines; subspaces of `F_2^k`). When `h < 2L - 1` the
  maximisers carry slack and are not rigid in any decided case. Regime boxes have `L = j = floor(w/2)` with no
  arithmetic relation to the profile, so the BCO is expected **not** to apply generically (**heuristic**; evidence only
  up to `n = 128`).

### 2.4 Verdict on (a)

* **Unconditional bound: false** (Cor 2.2). Any valid bound must use the character gap.
* **Conditional Liouville-capped bound: proved** (Cor 2.4 + Prop 2.5): balanced + (R) + (N) ⇒ `reach <= 4j <= 2w`.
  Its hypothesis (R) is a design coincidence; it holds in M8 and explains M8's failure, but there is no reason (and no
  evidence) that it holds on generic regime profiles. So it does **not** settle the decision rule against the route.
* **The exact obstruction to an unconditional regime bound.** Any proof of `reach <= 2w(1+o(1)) + O(1)` must produce,
  for every regime profile, a form in `cap_l J(A_l)_e \ J(A)_e` with `e <= 2w(1+o(1))`. The known families cost:
  slab-power products `prod_{p in Q} s_p^{2j} prod_{x in X} k_x^2`, when their membership in every `J(A_l)` is
  certified by the range mechanism, need `|X^c| + 1 <= 2^{|Q|}` (the patterns in `Q` must separate the points outside
  `X`, base included), hence degree `>= min_q (2jq + 2(n - 2^q)) ~ 2j log_2 n` (cf. killer §3); implied-slab powers
  need a rigid level set (design coincidence); forms `phi^e` with `phi` having all coefficients distinct (so that every root shortens it automatically) have range
  `>= ~ j d` (useless). What is left is non-univariate obstructions, of which the models show some (sdc7 `L = 2`: 15 of
  22 defect dimensions are not spanned by powers of ternary forms; `t6`), with no known structure. Alternatively one
  could bound `reach <= truth_L` by exhibiting exact collision-free base points at order `~2w`, which is the open
  barrier question of rounds 12–16. Neither is available; (a) is **open**.
* **Subsets.** Any `A' subset A` certifies `Bs V(A,m) subset Bs V(A',m)`, and Corollary 4.1 holds verbatim with
  `(SHIFT)_{A'}` for some `A' subset A` (the forms it produces lie in `V_H(A', m') subset V_H(A, m')`, which rigidity
  covers). So even where the BCO applies it is not a barrier for the route: `A' = A cap {|sigma| <= h - 1}` removes it.
  On M8 this raises the reach from 6 to 7 `= zeta_pair` (`t8`, mod p), and no further (§3.5).

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

## 4. Part (b): the syl8 `e = 6` obstruction, identified

**Theorem 4.1 (exact).** Let `A` = M8 (syl8, `L = 2`, `d = 7`, the 7 patterns are the complements of the 7 lines of the
Fano plane on the non-base points `0..6`), and `sigma = k_0 + ... + k_6 = -k_8` the **base-point coordinate**. Then
1. `sigma(A) = [-3, 3]` and `sigma^7 in J(A)_7` (exact);
2. `sigma^6 notin J(A)_6` (exact over `Q`, `t15`);
3. `sigma^6 in J(A_l)_6` for every `l in L_pair` (proved: Theorem 2.3, with the ranges below; also checked exactly for
   four shifts, `t15`);
4. the defect `cap_{l in L_pair} J(A_l)_6 / J(A)_6` is one-dimensional (mod p, `t1`; the invariant part is computed
   exactly over `Q`, `t2`), hence it is spanned by the class of `sigma^6`. The "one-dimensional defect that is not a
   product of column forms" of round 27 (`t6_obstruction`) is the **sixth power of the base coordinate**.

*Proof of 3.* `sigma(e_x) = 1`, so `sigma(A_{+-e_x})` has at most 6 values. For `l = e_x - e_y`, `sigma(l) = 0`; the top
level `{k in A : sigma(k) = 3}` is exactly the set of 7 indicator vectors `1_l` of Fano lines (exhaustive enumeration),
and `1_{l'} - 1_l` has weight `|l Delta l'| = 4` for distinct lines, so no top point survives a weight-2 move: the top
level is rigid and `sigma(A_{e_x - e_y}) subset [-2, 2]` (in fact equal). By central symmetry the bottom level behaves
the same way for `-l`. Theorem 2.3 applies with `R = 7`. QED

**Why `|sigma| <= 3`.** Each non-base point lies on 3 lines, hence in `7 - 3 = 4` patterns, so `sum_p s_p = 4 sigma`
with `|s_p| <= 2`: `|sigma| <= 14/4`. This is Proposition 2.5 for the Sylvester profile (`f = 4/7` with uniform
weights; the averaging bound is attained). The top level is the "design": the Fano lines are the 3-sets meeting
every line complement in at most 2 points.

**Why the weight-4 shifts remove it.** For two lines `l, l'` (meeting in one point) the shift `1_{l'} - 1_l` has weight
4 and `sigma = 0`, and carries the top point `1_l` to the top point `1_{l'}`; on `A_{1_{l'} - 1_l}` the base coordinate
again takes 7 values and the range argument no longer applies. This matches the killer's observation that
`(SHIFT_6)` holds with all `{-1,0,1}^7` shifts and that "the rescue needs shifts of weight `>= 4`" (round 27 §5.2(b)).

**Symmetry.** The coordinate symmetry group of M8 is `GL(3,2)` (order 168; the base point is fixed, `t2`), it fixes
`sigma`, and the defect is its trivial representation. Over `Q`: `dim J(A)_6^G = 2`, `dim J(A_l)_6^G = 3` for both
orbit representatives, and the invariant obstruction found by linear algebra is `F == c . sigma^6 mod J(A)_6^G`
(`t5`). It already lies in the first-order closures `J1(A, l)` (`t4`).

**Predictions from the same mechanism (not computed: `d >= 15`).** Sylvester boxes with `L = 2^{m-1}` have the rigid
subspace level (§2.3; IP-verified for `(n, L) = (16, 2), (16, 4), (32, 8), (64, 8), (64, 16), (128, 16)`), hence
`reach(A, L_pair) <= 2h = 4L - 2` whenever `(N)` holds there. For `(n, L) = (16, 2)` that is `reach <= 6`.

**The next failure of M8 (`e = 7`).** `delta(7) = 34` on `A`; after deleting the 14 points `|sigma| = 3` the subset
`A'` has `reach = 7` and `delta(7) = 28`, and that whole defect is spanned by 70 powers `phi^7` of ternary forms
(`t9`). Removing the extreme levels of these forms one at a time does not reduce it (28, 27, 26, 29 over four steps).

## 5. Part (c): what `(SHIFT)_A` on regime boxes would require

### 5.1 The statement

Regime: `n - 1 >= 316 ln 4r` (`n >= 2966`-type sizes), `w = sqrt(8 d ln 4r)`, `j = floor(w/2)`, `m' = ceil(4w) ~ 8j`.
By §2.4 and Theorem 3.1 the hypothesis of round 27's Corollary 4.1 may be weakened to:

> **(SHIFT\*)\_A.** There is a chain `S_0 subset ... subset S_{m'} subset A` (in particular: a subset `A' subset A`, constant
> chain) satisfying `(CH_e)` w.r.t. `L_pair` over `Q` for `e <= m' - 2`, and the top-degree height condition (3.3) with
> `H = U^2`, `U = |A|`.

`(SHIFT*)_A` for every regime profile still implies `M = O_C(log log R)` (the forms produced lie in `V_H(A, m')`, which
rigidity covers). Known: the constant chain on `A` satisfies `(CH_e)` for `e <= j - 1` (round 27 Prop 3.5). Needed:
`e in [j, 8j]`, with the decisive window `(2w(1+o(1)), 4w) = (4j(1+o(1)), 8j)` beyond every zonotope.

### 5.2 Necessary conditions (proved)

1. **No rigid implied slab below `m'`** (Theorem 2.3): for the sets `S_e` used, no integer form `phi` with
   `|phi(S)| - 1 <= m' - 2`, `(N_phi)`, and a level set rigid under every root `l` with `phi(l) = 0`. For a balanced
   profile and the constant chain on `A`, the base coordinate must have a **non-rigid** top level (or
   `sigma^{2h} in J(A)`), because `2h <= 4j < m' - 1` (Prop 2.5). This is a real constraint on `A`, but a subset
   `A cap {|sigma| <= h-1}` evades it.
2. **First-order property** (Lemma 3.2): `cap_l (J(A) + R . d_l J(A))_e = J(A)_e` for `e <= m' - 2`: a statement about
   `J(A) = J(C) : J(B)` alone.
3. **Linkage form** (Prop 3.3): `J(B)_{2d-e} = sum_l mu_l J(B_l)_{2d-e-deg mu_l} (mod J(C))` for `e <= m' - 2`, i.e.
   in degrees within `8j = o(d)` of the socle degree `2d`, for the explicit tail `B = Union_p {|s_p| > j}`.
4. **No representation-theoretic obstruction**: if a group `G` of coordinate symmetries preserves `A` and `L_pair`,
   `(SHIFT_e)` fails on the `rho`-isotypic part whenever
   `mult_rho M_e(A) > sum_{orbits O of L} dim Hom_{Stab(l_O)}(rho, M_e(A_{l_O}))` (Frobenius reciprocity: the sum over
   an orbit is the image of `Ind_{Stab}^G M_e(A_{l_O}) -> M_e(A)`). For actual-cluster profiles (`G = {+-1}`) this is
   the plain count `H_A(e) <= sum_{l in L_pair/+-} H_{A_l}(e)`, which always holds in the regime (`|L_pair| ~ d^2`). So in the regime
   **no counting obstruction exists**: any failure must be structural (implied slabs, designs, or the unexplained
   non-univariate residue of §3.4).

### 5.3 Does any asymptotic tool decide it? No.

* **Hilbert-function tools** (Macaulay's growth bound, Gotzmann persistence, Green's restriction theorem, generic
  initial ideals, Cayley–Bacharach / linkage) compute or bound `dim J(S)_e` for one set `S`. `(SHIFT_e)` is a statement
  about the **position** of `~d^2` subspaces `J(A_l)_e subset R_e` whose codimensions sum to far more than `H_A(e)`;
  generic subspaces of these codimensions would intersect in `J(A)_e`, so dimension data are consistent with
  `(SHIFT)` at every degree and cannot refute it; nor can they prove it, since the intersection of non-generic
  subspaces is not determined by their dimensions (the M8 failure is caused by the rigidity of 7 extremal points, a
  property of the configuration that no dimension count of `A` or of a single `A_l` sees). Gotzmann persistence needs
  generation in low degree with maximal growth (`J(A)` has generators in all degrees `2j + 1 .. rho`); generic initial
  ideals need a generic change of coordinates, which destroys the fixed
  root directions `L_pair`, and `gin` does not commute with intersections.
* **Gröbner (staircase) certificates**: `N(A)_e subset Union_l N(A_l)_e` for a degree-compatible order is a combinatorial
  sufficient condition (§3.5); it is the regime-scalable kind of statement (standard monomials of slab boxes are
  governed by lower sets), but it is weak already in the models (on M8' it certifies `e <= 7`, while `(SHIFT_e)` holds for `e <= 10`).
* **Linkage**: turns `(SHIFT)` into a statement about top-degree leading forms of an explicit union of slab tails
  (Prop 3.3). The intersection over `p` of slab-tail ideals is again an intersection problem; I found no induction on
  the number of slabs that controls it.
* **Asymptotic heuristic (not a proof).** In degrees `e <= 8j = o(d)` the leading-form ideal of `A` is the cube ideal
  plus the slab ideal `(s_p^{2j+1})` plus syzygy corrections. In the idealised model "independent slab forms, no cube"
  `J(A) = (y_p^{2j+1})`, `J(A_l) = (y_p^{2j+1-|y_p(l)|})` are monomial ideals and `(SHIFT_e)` holds exactly up to the
  product degree `2j . sep(L_pair) >= 2j log_2 n >> 8j` (separating families). If the only real obstructions are the
  implied-slab powers of §3.4, and generic regime profiles have no rigid levels (§2.3), this heuristic is favourable.
  It is contradicted in the small models only by rigid level sets of the symmetric profiles (M8 entirely; sdc7 for 7 of
  its 22 defect dimensions) and by the unexplained non-univariate residues, whose regime behaviour is unknown.

**Answer to (c).** `(SHIFT)_A` on regime boxes requires, at the very least, the absence of rigid implied slabs of
range `<= 8j` (a property of extremal configurations of the slab LP, decidable by integer programming profile by
profile) **and** a mechanism excluding non-univariate obstructions in degrees `(4j, 8j)`, for which no tool exists.
No asymptotic theorem known to me (or found in the record) decides it in either direction.

## 6. Consequences for the pre-registered decision rule

The rule is unchanged (it was fixed before the data). What the theory contributes to reading the `d = 7 -> 8 -> 9` data:

1. **The theory does not settle the rule.** It proves neither `reach <= zeta + c` on regime boxes (false in general,
   Cor 2.2; conditional and design-dependent in the form of Cor 2.4) nor any mechanism forcing `reach > zeta`.
2. **Use `zeta_pair`.** The criterion's hypothesis is `Y^l != 1` for `l in L_pair`; the fair comparison is with
   zonotopes whose generators lie in `L_pair`. On the four `L >= 2` models `zeta_pair = zeta_all` (§1), but at `d = 8, 9`
   the two can differ, and only `zeta_pair` is the same-hypothesis competitor.
3. **Design effects.** A failure at `e = 2h(sigma)` with a rigid top level of the base coordinate (Sylvester boxes with
   `L` a power of 2; M8), or more generally a defect spanned by powers of implied slabs with rigid level sets (§3.4),
   is predicted by Theorem 2.3. Such failures are arithmetic coincidences of small symmetric profiles; they are evidence
   about designs, not about generic regime profiles. Recording, for each failure, `h(sigma)`, the rigidity of its top
   level (`t7`/`t10` style) and the dimension of the span of power obstructions (`t6`) makes this visible.
4. **Vacuous favourable cases.** Product (block) profiles satisfy the favourable pattern trivially (Cor 2.2), and
   vertex cuts (`L = |P| - 1`, like M8') are cube-bound rather than slab-bound. Neither resembles the regime
   (`j/|P| -> 0`, character gap). Favourable evidence should come from gap-satisfying, slab-bound (`L <= |P|/2`) profiles.
5. **Constant chain versus chains.** The computed "reach" is the constant-chain reach. Corollary 4.1 works verbatim
   with any subset or chain (Theorem 3.1, §5.1), and `reach_chain >= zeta_L` by construction. So "reach <= zeta" on
   every box would show that the *constant chain* never beats zonotopes, not that the interpolation route is closed.
   I do **not** propose to change the rule; I flag that the cheap extra statistic `reach(A cap {|sigma| <= h-1})`
   (one more rank computation per box) measures whether the BCO is the only obstruction (on M8: 6 -> 7 = `zeta_pair`).
6. **Probabilities (judgement).** The round-27 estimate (favourable trend ~15%; proof of `(SHIFT)` at regime size within
   a few rounds 5–10% conditional on that) is not moved much. The theory removes one reason for pessimism (the M8 failure is a design
   effect, and part of the sdc7 failure is a rigid-level effect of the same kind) and one reason for optimism (M8' succeeds because its cuts are weak and
   design-free). If the `d = 8, 9` data come out unfavourable on slab-bound, non-design, gap-satisfying boxes, the
   referral to experts is the right call; nothing in this note gives an internal route to a proof of `(SHIFT)` at
   regime size.

## 7. Reproducibility: scripts, outputs, labels

All under `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/theory_checks/`.
Heavy runs used `prlimit --as=3500000000`, at most two of my processes at a time. Boxes are built by
`tlib.build_box(n, slabs)` exactly as the round-27 killer's `klib.build_box` (patterns `syl(3)`, `sdc(7)`; random boxes
by the killer's `k8/k11` seeds and generator, re-implemented in `t6/t7`).

| script -> output | content | label |
| --- | --- | --- |
| `tlib.py` | library: standard-monomial least spaces mod p, shift sums, exact `J`-membership (fmpz), power obstructions | — |
| `t1_profiles.py` -> `t1_{syl8,sdc7}_L{1,2,3}_out.txt` | full defect profiles w.r.t. `L_pair`; graded `H_A`, `H_B`; linkage identity | mod p (`P1`) |
| `t2_syl8_invariant.py` -> `t2_syl8_L2_e6_out.txt` | `GL(3,2)` (168), 8 invariant sextic orbits, invariant parts of `J(A)_6`, `J(A_l)_6`; obstruction | **exact** |
| `t3_firstorder.py` -> `t3_syl8_L2_out.txt` | first-order closure `J1(A,l)` vs `J(A_l)`, `e <= 8` | mod p |
| `t4_firstorder_inv.py` -> `t4_out.txt` | obstruction lies in `J1(A,l)_6` (invariant parts) | mod p |
| `t5_identify.py` -> `t5_out.txt` | obstruction `== c sigma^6 mod J(A)^G` | **exact** (sympy) |
| `t6_powers.py` -> `t6_syl8_2_e6_out.txt`, `t6_misc_out.txt` | power obstructions `phi^e`, `phi` ternary | mod p |
| `t7_bco_scan.py` -> `t7_bco_scan_out.txt` | `h(sigma)`, rigidity of the top level, 6 models + 27 random boxes | **exact** (enumeration) |
| `t8_subset.py` -> `t8_syl8_L2_sig1_out.txt` | `A cap {|sigma| <= 2}` on M8: reach 7 | mod p |
| `t9_subset_greedy.py` -> `t9_syl8_L2_out.txt` | greedy removal of extreme levels; `e = 7` defect spanned by powers | mod p |
| `t10_bco_ip.py`, `t10b_h_only.py` -> `t10_k3_out.txt`, `t10_out.txt`, `t10c_out.txt`, `t10b_6_12_out.txt` | Sylvester `n <= 128`: `h`, rigidity, box-zonotope (IP, HiGHS) | exact when status 0 (IP optimality; floating-point solver) |
| `t11_zeta_pair.py` -> `t11_out.txt` | `zeta_pair` by exhaustive DFS over `L_pair` generators | **exact** search |
| `t12_small_search.py` -> `t12_out.txt` | ~400 small boxes (`n <= 6`): no `reach > zeta_pair` | mod p / exact search |
| `t13_staircase.py` -> `t13_*_out.txt` | staircase (standard-monomial) criterion | mod p |
| `t14_linkage.py` -> `t14_syl8_L2_out.txt` | linkage form of `(SHIFT_e)`, `e = 5,6,7`; linkage identity in `C_l` | mod p |
| `t16_random_rigidity.py` -> `t16_out.txt` | (R) on random (unbalanced) sign-column profiles, `n <= 128` (IP) | exact when status 0/2 |
| `t17_product_check.py` -> `t17_out.txt` | Theorem 2.1 sanity: reach and `H` of small products add/convolve | mod p |
| `t15_sigma_exact.py` -> `t15_out.txt` | `sigma^6 notin J(A)_6`, `sigma^7 in J(A)_7`, `sigma^6 in J(A_l)_6`; top level | **exact** |

Notes file: `theory_notes.md`. Everything here is **unrefereed**.

