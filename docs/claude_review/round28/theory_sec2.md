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
