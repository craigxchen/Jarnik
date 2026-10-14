# Round 13 — Referee report on judge findings J5, J6, J7 (round12/README.md §4)

Status: COMPLETE. Referee: independent re-derivation and own code (no judge script was reused).
Settings used: `T = 1`, Sylvester profiles `n = 2^k`, `r = d = n - 1` sign columns (`e_p = 1`), `w = sqrt(8d ln 4d)`,
`j = floor(w/2)`, `m' = ceil(4w)`, `C = 1` unless stated (§1, §5); exact models `M8` (`Syl(8)`, `L = 2`, `m = 8`) and one-/multi-cut
ternary boxes (§2-§4).

## 0. Summary of verdicts

| item | judge's claim | verdict | main change |
| --- | --- | --- | --- |
| J5 uniform `Q_0` | false on Sylvester `2^12..2^17` for every `Q_0`, all `abar >= 5.42` | **confirmed** (own code, all margins) | **extended to `2^12..2^20`** with projection-cheap Fano blocks (32-class span refutes every `Q_0 < Q_5 ~ 32`); `C = 1` needed literally (threshold moves by `<= 8.7e-5 ln(1/C)`); convention `P(pi) >= 1` made explicit |
| J5 per-subspace | fails only at `n = 4096` from `W ~ 2.7 W_min`; not refuted at `2^13, 2^14` | **confirmed as computed, superseded** | **ZE-TC^coll_sigma is false at `n = 2^12, 2^13, 2^14` for all `W >= ~0.95 W_min`** (Fano spans; Lemma 1.3: `rho(P_j) <= 3j`) |
| J6(a) fat points | `TC ⊇ Bs_R(L_m)`; cheap spans = linear base cycles | inclusion **proved** (real and complex, multiplicity `m - tau_A(v)`, integer `v`) | "exactly" only in the BDP sense of linear cycles `L_I`; the base locus contains further rational subspaces (the Fano spans of `M8`, Prop 2.3) |
| J6(b) one cut / SLP | `J(A_S(L)) = (x^3, sigma_S^(2L+1))` via quoted SLP | **confirmed**; SLP **proved here** (`sl_2`), reference statement checked | extended to arbitrary signs and any number of disjoint cuts (Thm 3.3); `rho` formulas; exact checks on 7 boxes |
| J6(c) `M8` | dims `[0,0,0,7,49,203,637,1541,2961]`; triples give `J(A)`; pairs add nothing (mod `p`) | **confirmed exactly over `Q`** | only the 7 Fano triples matter (parity relation); column pairs/non-line triples add nothing **structurally**; `b_0 = 2` fails exactly for column pairs and for 51,933 integer planes, **open for all planes**; general planes through the half-sum direction give 14 of the 21 missing dimensions |
| J7(a) hierarchy | `coll_eps => coll_sigma => MIN => MIN^pair`; `MIN* => MIN` | **confirmed** | `MIN* => MIN` needs `C <= 1`; first implication needs `P(pi) >= 1`; compactness noted |
| J7(b) empty gap region | witness has `G = ∅`; weaker on non-empty `G` open | **confirmed** (`104.94 > 2`) | **dichotomy** (Prop 5.2): at Sylvester orders `(2jk, m_max]`, MIN^pair fails whenever `max eps_xy >= 2 ln(n/2)`, holds vacuously below `2 ln((n-1)/2)`: no escape on non-empty gap regions except in a window of width `< 2/(n-1)` |

**Bearing on `Loc(b_0)`.** Proved: `Loc(1)` for disjoint cuts (all degrees). Exact: `Loc(3)` in `M8` with Fano triples. The
non-width part of `J` comes from lattice relations among columns, which exist in every Sylvester regime box and are
admissible there (cost `3j < ceil(4w)`). New precise open question (§6.4): growth in `t` of the interpolation degree `rho_t`
of XOR-closed column blocks; if `(m' - rho_t) 2^t > d` for some `t`, the per-subspace route fails asymptotically.

**No growth improvement or uniform bound is claimed.**

## 1. J5 — failure of the uniform-Q0 collision hypothesis on Sylvester n = 2^12..2^17

**Verdict: CONFIRMED (proof re-derived, all margins recomputed with independent code), with two small corrections
(the role of `C`, and the near-extremal constant), and two STRENGTHENINGS (§1.6), both from a new projection-cheap
family the judge did not use (Fano blocks `{a, b, a+b}` of Sylvester columns, projection cost `3j`, not `6j`):**
* **the uniform-`Q_0` hypothesis is false on Sylvester `n = 2^12, ..., 2^20`** (not only `2^17`), for every `Q_0`, at
  every level `abar >= 5.4161` for `k <= 19` and `abar >= 9.682` for `k = 20` (both include all `W >= 3n ln n - 14n`);
* **the per-subspace hypothesis ZE-TC^coll_sigma is false on Sylvester `n = 2^12, 2^13, 2^14`** at every level
  `abar >= 4.892 / 6.204 / 6.417` respectively, i.e. from `W ~ 0.89 / 0.95 / 0.85 W_min`, which includes the whole
  near-extremal range. README §4.1's "no failure from coordinate spans at `n = 8192, 16384`" is true but no longer
  informative; its "failure at `n = 4096` only from `W ~ 2.7 W_min`" is superseded.

Scripts: `refJ_checks/j5_recheck.py` (judge's claims), `refJ_checks/j5_extend.py`, `persub_fano.py`,
`parity_box_rho.py` (strengthening); outputs `*_out.txt`.

### 1.1 The hypothesis being refuted (stated precisely)

Setting: `T = 1`, profile columns `beta_p`, `d = n - 1`, `w = sqrt(8 d ln 4r)`, `j = floor(w/2)`, `m' = ceil(4w)`,
`A = {k in Z^n : sum k = 0, |k_x| <= 1 (x != base), |<k,beta_p>| <= w for all p}`, `u`-coordinates with
`u_base = 0`, `<k',u> = sum_x k_x u_x`. `kappa(u) = max_{c in V_H} kappa_c(u)` (any choice of `H`; only the upper
bound of Lemma 5.0, valid for all `c in V(A,m')`, is used). For a partition `pi` of `[n]`,
`Delta_pi = {u : u_base = 0, u_x = u_y whenever x ~_pi y}` and `P(pi)` = number of within-block pairs; note
`P(Delta_pi) = P(pi)` (two distinct non-base blocks take independent values on `Delta_pi`). Distances are sup-norm.

> **ZE-TC^coll_eps(A, m'; Q_0)** (referee_t1rigour §6.2). For every `hat u` with `|hat u|_inf = 1`, if
> `dist(hat u, Delta_pi) >= sigma := (1/2) e^{-Q_0 abar/2}` for every partition `pi` with `P(pi) >= N_2/Q_0`, then
> `kappa(hat u) > eps_0 = C d e^{-W/4}/(m'+1)`.

Convention: "collision subspace" means `P(pi) >= 1` (t1rigour Lemma 6.6 needs it), so the family is
`{pi : P(pi) >= max(1, N_2/Q_0)}`. (Read literally with `Q_0 > N_2`, the family would contain the discrete partition,
whose `Delta_pi = {u_base = 0}` contains every point, and the hypothesis would be vacuous; that reading is not intended
and is not used by Thm 6.8, which applies only to `P(L) >= 1`.)

Here `abar = W/(2d) + 2 ln C'`. Below `C = C' = 1` unless stated; then `eps_0 = (d/(m'+1)) e^{-d abar/2}`.

### 1.2 Two general lemmas (proved; they contain the judge's cases A and B)

**Lemma 1.1 (case A, general form).** Let `L` be a linear subspace (in `u`-coordinates) with `Phi_c|_L = 0` for all
`c in V(A,m')` (e.g. `L = span F` with `rho(pi_F(A)) < m'`, Lemma 5.0). Let `y in L`, `|y|_inf = 1`, have
collision partition exactly `pi(L)` (i.e. `y_x = y_{x'}` iff `x, x'` collide on all of `L`), and let `delta_y` be its
least nonzero gap. Then for every `Q_0 < N_2/P(L)` and every level with `sigma <= delta_y/2`, the hypothesis fails at
`hat u = y`.

*Proof.* `kappa(y) = 0 <= eps_0`. If `P(pi) >= N_2/Q_0 > P(L) = P(pi(L))`, then `pi` does not refine `pi(L)` (a
refinement has at most as many within-block pairs), so some block of `pi` contains `x, x'` with `y_x != y_{x'}`, hence
`|y_x - y_{x'}| >= delta_y`. Every `z in Delta_pi` has `z_x = z_{x'}`, so `|y - z|_inf >= delta_y/2 >= sigma`. QED

For Sylvester take `L = L_4 = span(1_{S_1},...,1_{S_4})` (after subtracting the base value), `S_i = {x : x_i = 0}`,
`y = sum_{i<=4} 2^{i-1}(1_{S_i} - 1_{S_i}(base))/15`, renormalised to `|y|_inf = 1` (renormalisation only enlarges
gaps, so `delta_y >= 1/15`). `<k',y>` restricted to `A` is `sum 2^{i-1} k(S_i)/15` because `sum k = 0`; on `A`,
`<k,beta_{a}> = 2k(S_a)` for the sign column of `a`, so `|k(S_i)| <= j`, the image of `F_4` lies in `{-j..j}^4`, and
`rho(pi_{F_4}(A)) <= 8j`. The script certifies at 80 digits that `w/2` and `4w` are at distance `> 10^-20` from
the nearest integer (so `j`, `m'` are determined rigorously), hence `8j = 8 floor(w/2) <= 4w < ceil(4w) = m'` with
`4w` not an integer, so `o_4 = m' - 8j >= 1`
(exactly `o_4 = 8, 6, 7, 7, 7, 1` for `k = 12..17`). So case A holds for every `Q_0 < Q_4 = N_2/P(L_4) = 16(n-1)/(n-16)`
at every level `Q_0 abar >= 2 ln 15 = 5.4161`; it needs **no** assumption on `C`.

**Lemma 1.2 (case B, general form).** Let `L` be as in Lemma 1.1 with transversal order `o` in the sense of Lemma 5.0
(`kappa_c(hat(y + eta v)) <= 4^{m'} eta^o` for `y in L`, `|y|_inf = 1`, `|v|_inf <= 1`, `0 <= eta <= 1/2`). Suppose
`pi(L)` has two blocks (for simplicity; this is the case used), each of size `n/2`, and `Q_0 >= N_2/P(L)`. With
`c_n = n/(n-1)` and `a = m' ln 4 - ln(d/(m'+1)) + ln(1/C)`, the hypothesis fails at every level
`abar >= g*(Q_0) := 2[a + o ln(c_n Q_0)]/(o Q_0 - d)` provided `o Q_0 > d`.

*Proof.* `y` = the normalised indicator of the non-base class (values `0` and `+-1`). `v`: inside each class
distinct values `2 rho/n`, `rho = 0..n/2-1` the rank, base point rank 0; so `|v|_inf < 1`, `v_base = 0`. Put
`u = y + eta v`, `eta = 2 sigma Q_0 c_n`. Then `|u|_inf in [1-eta, 1+eta]` and `dist(hat u, .) = dist(u, .)/|u|_inf`.
(i) If a block of `pi` meets both classes, `dist(u, Delta_pi) >= (1 - eta)/2`. (ii) Otherwise `pi` refines `pi(L)`;
`P(pi) = sum binom(b_i,2) <= (b_max - 1) n/2`, so `P(pi) >= N_2/Q_0` gives a block of size `b_max >= 1 + (n-1)/Q_0`
inside one class, on which `v` has spread `>= 2(b_max-1)/n`, so `dist(u, Delta_pi) >= eta (n-1)/(n Q_0) = 2 sigma`.
In both cases `dist(hat u, Delta_pi) >= sigma` as soon as `eta <= 1/2` and `(1-eta)/(2(1+eta)) >= sigma`. Lemma 5.0
gives `kappa(hat u) <= 4^{m'} eta^o`, and `4^{m'} eta^o <= eps_0` rearranges (taking logs, using
`ln eta = ln(c_n Q_0) - Q_0 abar/2`) to `abar (o Q_0 - d)/2 >= a + o ln(c_n Q_0)`. Finally the failure inequality
itself forces `eta^o <= 4^{-m'} eps_0 <= 4^{-m'}` (as `eps_0 <= 1` in the range considered), so `eta <= 1/4`, and then
`(1-eta)/(2(1+eta)) >= 3/10 >= sigma`. QED

*Monotonicity.* `d g*/dQ` has the sign of `o(1 - ln(c_n Q)) - d/Q - a`, which is `< 0` for `c_n Q >= e` because
`a > 0`. So on `[Q_4, infinity)` the worst value is `g*(Q_4)`. (The judge's remark "analytic sign agrees with a grid
check" is correct; the sign is in fact immediate.)

### 1.3 Recomputation (`j5_recheck_out.txt`; exact integers for R1, Fractions for `Q_b`, 80-digit logs)

| `k` | `n` | `j` | `m'` | `8 binom(m'+d-1,d) <= 3^d` | `o_1` | `o_4` | `Q_4` | `o_1 Q_4 - d` | `g*(Q_4)` (C=1) | `eta` at `g*` | `(3n ln n - 14n)/(2d)` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 12 | 4096 | 281 | 2256 | exact True | 1694 | 8 | 16.0588 | 23108.6 | 0.67769 | 0.0696 | 5.478 |
| 13 | 8192 | 412 | 3302 | exact True | 2478 | 6 | 16.0294 | 31529.7 | 0.72642 | 0.0475 | 6.517 |
| 14 | 16384 | 602 | 4823 | exact True | 3619 | 7 | 16.0147 | 41574.1 | 0.80446 | 0.0255 | 7.557 |
| 15 | 32768 | 878 | 7031 | exact True | 5275 | 7 | 16.0073 | 51671.7 | 0.94340 | 0.0084 | 8.596 |
| 16 | 65536 | 1278 | 10231 | exact True | 7675 | 7 | 16.0037 | 57293.1 | 1.23794 | 8.0e-4 | 9.636 |
| 17 | 131072 | 1858 | 14865 | exact True | 11149 | 1 | 16.0018 | 47333.4 | 2.17681 | 4.4e-7 | 10.675 |
| 18 | 262144 | 2695 | 21568 | (lgamma) True | 16178 | 8 | 16.0009 | **-3280.2** | — | — | 11.715 |

All numbers agree with `judge_checks/uniformQ0_sylvester_out.txt`. In addition the script checks the failure
inequality **directly** (not via the closed form) at `g*(Q_4) + 10^-40` for `Q_0` on a geometric grid of ratio 1.01
over `[Q_4, 10^6]`, together with both side conditions (`eta <= 1/2`, mixing distance); all pass. Each Sylvester box
is in the regime (`d >= 316 ln 4r` with `r = d`) and the sign columns are realised with `e_p = 1`.

Lemma 1.2's case (ii) uses only `P(pi) >= N_2/Q_0`, so it is unaffected by the convention `P(pi) >= 1` (for
`Q_0 > N_2` it still gives `b_max - 1 >= (n-1)/Q_0`); the monotonicity is analytic, so the grid's upper end `10^6`
plays no role.

**Conclusion (proved, for `C = 1`).** For `12 <= k <= 17`, every `Q_0 > 0` and every level
`abar >= 2 ln 15 = 5.4161` (in particular every `W >= 3n ln n - 14n`), ZE-TC^coll_eps`(A, m'; Q_0)` is false on the
Sylvester box. (`Q_0 < 1`: the family is empty and `y` itself refutes.)

### 1.4 Corrections to the J5 text

1. **`C`.** README says "any primes, any base point, `C <= 1`". Case A is `C`-free. Case B's threshold is
   `g*(Q_4) + 2 ln(1/C)/(o_1 Q_4 - d)`; the coefficient is `<= 8.7e-5` for all `k = 12..17`, so the statement holds for
   `C >= e^{-(5.4161 - g*(Q_4))(o_1 Q_4 - d)/2}` (e.g. `C >= e^{-37000}` at `k = 17`), not literally for every `C <= 1`.
   For `C > 1`, `eps_0` grows and the failure only gets easier, but `abar` is then `W/(2d) + 2 ln C`; stating J5 at
   `C = 1` is cleanest.
2. **Near-extremal constant.** The judge's script uses `W_min = 3n ln n - 13.91 n` while the README says `-14n`; with
   `-14n` the level at `k = 12` is `5.478` (script) and still `> 5.4161`. No change to the conclusion.
3. **What the judge's mechanism shows for `k >= 18`.** With coordinate spans only, Lemma 1.1 gives: for every
   `k >= 12` the uniform hypothesis is false for every `Q_0 < Q_4 = 16(n-1)/(n-16)`; case B stops at `k = 18` because
   `o_1 Q_4 < d`. By Lemma 5.0 any perturbation argument off a single subspace `L` needs `o_L Q_0 > d`, so the *threshold*
   `Q_0` below which case A works is what limits the method. §1.6 raises that threshold from `Q_4 ~ 16` to `Q_5 ~ 32` with
   projection-cheap Fano blocks, which extends the failure to `k <= 20`.

### 1.5 The corrected hypothesis (per-subspace radii)

> **ZE-TC^coll_sigma(A, m').** For every `hat u` with `|hat u|_inf = 1` and `hat u_base = 0`: if
> `dist_inf(hat u, Delta_pi) >= sigma_pi := (1/2) exp(-N_2 abar/(2 P(pi)))` for **every** partition `pi` of `[n]` with
> `P(pi) >= 1`, then `kappa(hat u) > eps_0 = C d e^{-W/4}/(m'+1)`. (Finite family: `B_n` partitions.)

It suffices for Theorem 6.8 (first part with `𝓛 = {Delta_pi}`, radii from Lemma 6.6, since `P(Delta_pi) = P(pi)`).
It is implied by the uniform form for every `Q_0` (J7(a), §5).

**Per-subspace recomputation (Lemma 1.2 with per-subspace radii; own derivation).** For a refinement `pi` of
`pi_b` (classes of size `s = n/2^b`) with `P(pi) = P`, a block of size `>= 1 + 2P/n` lies in one class, so
`dist(u, Delta_pi) >= eta P/(ns)`; the needed `eta >= (ns/P) e^{-N_2 abar/(2P)}` is increasing in `P` for
`P < N_2 abar/2`, so it binds at `P = P_b`. Mixing partitions: `dist >= (1/(2^b-1) - eta)/(2(1+eta))` and every
`sigma_pi <= e^{-abar/2}/2`. Failure iff `o_b Q_b > d` and `abar >= g*_b = 2[m' ln 4 + o_b ln(ns/P_b) - ln(d/(m'+1))]/(o_b Q_b - d)`.

| `k` | `b` | `o_b Q_b - d` | failure from `abar >=` | `W >=` | `W/W_min` |
| --- | --- | --- | --- | --- | --- |
| 12 | 2 | 436.3 | 17.9347 | 146885 | 3.27 |
| 12 | 3 | 472.8 | 14.9028 | 122054 | 2.72 |
| 12 | 1, 4 | < 0 | — | — | — |
| 13-20 | 1..4 | < 0 | — | — | — |

Confirmed (judge: `14.90`, `17.93`, `W >= 1.22e5 ~ 2.7 W_min`). **Superseded by §1.6:** with Fano blocks the
per-subspace form fails at `n = 2^12, 2^13, 2^14` already from `W ~ 0.85-0.95 W_min`. **Scope:** at `k >= 13` the coordinate spans are
not refuting; this says nothing about other subspaces of `TC`. By the same Lemma 5.0 argument the per-subspace form
fails as soon as `TC_H` contains a linear subspace `L` with `o_L N_2/P(L) > d` (nearstrata Prop 5.1(ii) gives this at
large levels; Lemma 1.2's version gives the explicit level), so the per-subspace hypothesis has the necessary
condition **`o_L N_2/P(L) <= d` for every linear `L subset TC_H`, at the levels where it is applied**, which is
exactly the "strata half" of README §7. Status: open; consistent with all spans computed.

### 1.6 Strengthening: Fano blocks are projection-cheap at cost `3j`

**Lemma 1.3 (parity box).** For `j >= 1` let `P_j = {y in Z^3 : |y_i| <= j, y_1 + y_2 + y_3 even}`. Then
`rho(P_j) <= 3j` (every function on `P_j` is a polynomial of degree `<= 3j`). Exact: `rho(P_j) = 3j` for `j <= 4`
(`parity_box_rho_out.txt`; the lower bound holds in general by the hexagon zonotope in the slice `y_1+y_2+y_3 = 0`, not
needed here).

*Proof.* Peeling lemma: if `Y = Y_1 ⊔ Y_2` and a polynomial `g` of degree `e` vanishes on `Y_1` and nowhere on `Y_2`, then
`rho(Y) <= max(rho(Y_1), e + rho(Y_2))` (interpolate `f` on `Y_1` by `p_1`, then `(f - p_1)/g` on `Y_2` by `p_2`; take
`p_1 + g p_2`). Iterating with linear `g`, `rho(Y) <= max_i (i + rho(Y_i))` for any ordered partition into pieces
`Y_0, Y_1, ...` each cut out from the remaining ones by one linear form. Use the linear coordinates
`(y_1, y_2, z)`, `z = (y_1+y_2+y_3)/2 in Z`: `P_j = {|y_1|, |y_2| <= j, |2z - y_1 - y_2| <= j}`. Outer pieces: the slices
`H_z` (`z` fixed), ordered `z = 0, 1, -1, 2, -2, ...`, so `H_{±m}` has index `<= 2m`. Inner pieces of `H_m` (`m >= 0`; `m < 0`
by the symmetry `y -> -y`): the diagonals `y_1 + y_2 = c`, `c in [2m - j, min(2j, 2m + j)]`, each a set of `2j + 1 - |c|`
consecutive points of a line (`rho = 2j - |c|`), ordered by `|c|`. A diagonal with `|c| = t` has index
`<= #{c' in range : |c'| <= t} - 1`, which is `<= 2t` if `t <= j - 2m` and `<= t - 2m + j` otherwise; in both cases
`index + 2j - t <= 3j - 2m`. So `rho(H_{±m}) <= 3j - 2m` and `rho(P_j) <= max_m (2m + 3j - 2m) = 3j`. QED

**Lemma 1.4 (Fano blocks).** Sylvester profile, `n = 2^k`. For `a, b != 0`, `a != b`, the three sign columns
`a, b, a+b` satisfy `k(S_a) + k(S_b) + k(S_{a+b}) = 2 k(S_a cap S_b)` on `A` (each `x` lies in 3 or 1 of the three
hyperplanes, and `sum k = 0`), so `pi_{a,b,a+b}(A) subset P_j` and `rho <= 3j`. More generally, for `s` such blocks and
`c` further columns, all on linearly independent `a`'s (`2s + c <= k`), the image lies in `P_j^s x [-j,j]^c`, so
`rho <= (3s + 2c) j`, and by Lemma 5.0 the span `L_{s,c}` (functions `x -> f_1(a_1.x, b_1.x) + ... + g(a'.x)`, dimension
`3s + c`) has transversal order `o >= m' - (3s + 2c) j`. Its collision partition is the `2^{2s+c}` cosets, all of size
`n/2^{2s+c}`. In particular `L_{2,1}` (cost `8j < m'`) is a **32-class** subspace of `TC(A, m')` for **every** `k >= 5`.

(These spans are projection-cheap, not width-cheap: a Fano block has width cost `6j` but projection cost `3j`. They
are also exactly the source of the non-width generators of `J(M8)`, §4.)

**Consequences (proved; `j5_extend_out.txt`, exact `R1` for `k <= 20`, Fractions for `Q`, 60-digit logs, failure
inequality re-checked directly at the threshold).**
1. *Uniform `Q_0`.* Case A (Lemma 1.1) with `L_{2,1}` (32 classes, gap `1/31`) refutes every `Q_0 < Q_5 = 32(n-1)/(n-32)` at
   `Q_0 abar >= 2 ln 31` (for `Q_0 >= Q_4` this is `abar >= 0.43`); `L_4` handles `Q_0 < Q_4` at `Q_0 abar >= 2 ln 15`; and case
   B (Lemma 1.2, `L_1`) now only needs `Q_0 >= Q_5`, where `o_1 Q_5 - d > 0` holds for `12 <= k <= 20`:

   | `k` | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 |
   | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
   | `o_1 Q_5 - d` | 50526 | 71406 | 99645 | 136193 | 180181 | 225781 | 255614 | 225933 | 37249 | -527080 |
   | `g*(Q_5)` | 0.357 | 0.369 | 0.386 | 0.412 | 0.453 | 0.525 | 0.673 | 1.103 | 9.681 | — |
   | near-extremal `abar` | 5.478 | 6.517 | 7.557 | 8.596 | 9.636 | 10.675 | 11.715 | 12.755 | 13.794 | 14.834 |

   So ZE-TC^coll_eps`(Q_0)` is false for every `Q_0`, `n = 2^12..2^20`, at all `abar >= max(5.4161, g*(Q_5))`. Also, for
   **every** `k >= 12`, any `Q_0` that could work must exceed `Q_5 ~ 32` (Case A), which already exceeds the `Q_0 = 16(1 + O(w/n))`
   that t1rigour Thm 6.7 proves for width-cheap spans.
2. *Per-subspace radii.* The §1.5 argument with the spans `L_{s,c}` (`B = 2s + c` bits, `2^B` equal classes). Test point:
   `y = (N(x) - N(base))/(2^B - 1)` with `N = N_1 + 4 N_2 + 16 N_3` (`N_i` = the integer coded by the bits of block `i`;
   every function of the two bits of a Fano block lies in its span), so `y in L_{s,c}`, its values are distinct with gaps
   `>= 1/(2^B - 1)`; perturbation `v` = rank within each class as in §1.5:

   | `k` | best span `(s,c)` | classes | `o` | `o Q - d` | fails for `abar >=` | `W/W_min` |
   | --- | --- | --- | --- | --- | --- | --- |
   | 12 | (1,1) | 8 | 851 | 2724.7 | 4.8918 | 0.893 |
   | 12 | (1,0) | 4 | 1413 | 1561.1 | 5.2624 | 0.961 |
   | 13 | (1,1) | 8 | 1242 | 1753.5 | 6.2033 | 0.952 |
   | 13 | (2,0) | 16 | 830 | 5113.4 | 6.4161 | 0.985 |
   | 14 | (2,0) | 16 | 1211 | 3010.8 | 6.4161 | 0.849 |
   | 15..23 | none of `(1,0),(2,0),(1,1),(1,2),(0,2),(0,3)` | | | `< 0` | — | — |

   (The thresholds `4.8918 = 2 ln 7 + 1`, `6.4161 = 2 ln 15 + 1` are the crude mixing-partition bound, not the transversal
   one.) So **ZE-TC^coll_sigma is false on Sylvester `n = 4096, 8192, 16384` at every near-extremal `W`**; `n_1 > 2^14` in any
   use of it. Larger XOR-closed blocks (all `7` columns of a 3-dimensional `V`, 8 classes) have `rho <= 3` at `j = 1` and
   `rho <= 8` at `j = 2` (mod-`p` ranks give rigorous **upper** bounds on `rho`; `xor_family_rho_t3_out.txt`); if
   `rho <= 4j` held for all `j` (unproved), the per-subspace form would also fail at `n = 2^14` with more margin, but not at
   `2^15`. Asymptotically all these constructions have `o N_2/P(L) = O(j) = o(d)`, so nothing asymptotic follows.

## 2. J6(a) — the fat-point inclusion

**Verdict: the inclusion is CORRECT (proof below, with all quantifiers), and holds over `C` as well as over `R`, with a
slightly stronger multiplicity (`tau_A(v)` = number of values minus one, in place of the width). The sentence "the
cheap spans are exactly its linear base cycles" is correct only in the technical sense of linear cycles `L_I`
(spans of subsets of the fat points); it is FALSE if read as "the linear subspaces contained in the base locus are
exactly the cheap spans": in `M8` the 3-dimensional Fano spans `U_l` lie in `Bs(L_8)` but are not spans of any cheap
family (Prop 2.3).**

### 2.1 Setting

`A subset Z^d` finite (here `d = n - 1`, `k'`-coordinates), `R = Q[x_1..x_d]`, `I(A)` its vanishing ideal, `J(A)` the
ideal of leading (top-degree) forms of elements of `I(A)`. Dual space: forms in `u = (u_1..u_d)`; `f in R` acts as
`f(partial_u)`. For `v in R^d` write `l_v(x) = <x, v>` and `D_v = l_v(partial_u)` (derivative along `v` in `u`-space).
`V(A,m) = {c in R^A : sum_k c_k k^alpha = 0, |alpha| < m}`, `Phi_c(u) = sum_k c_k <k,u>^m`,
`TC(A,m) = {u in R^d : Phi_c(u) = 0 for all c in V(A,m)}`; the complex version `TC_C` uses `u in C^d`.
For the cluster's `V_H subset V(A,m)` one has `TC_H supset TC(A,m)`, so every inclusion below also holds for `TC_H`.

**Lemma 2.0 (duality; = compute R2, re-proved).** `J(A)_m = {f in R_m : sum_k c_k f(k) = 0 for all c in V(A,m)}`, and
`span{Phi_c : c in V(A,m)} = J(A)_m^perp := {Phi in Sym^m : f(partial)Phi = 0 for all f in J(A)_m}`.

*Proof.* `f in R_m` lies in `J(A)_m` iff `f - g in I(A)` for some `g` of degree `< m`, iff `ev(f) in E := ev(R_{<m})`
(`ev` = evaluation on `A`). `V(A,m) = E^perp` in `R^A`, so the first claim is `ev(f) in E <=> ev(f) perp E^perp`.
For the second, `f(partial) Phi_c = m! sum_k c_k f(k)` (as `f(partial)<k,u>^m = m! f(k)` for `f in R_m`), so
`Phi_c in J(A)_m^perp` for `c in V(A,m)`, and the map `c -> Phi_c` is `ev_m^T` (transpose of `ev_m : R_m -> R^A`) under
the (nondegenerate, characteristic 0) apolarity pairing. Dimension count: `dim ev_m^T(E^perp) = dim(E + im ev_m) - dim E`
and `dim J(A)_m^perp = dim R_m - dim ev_m^{-1}(E) = rank ev_m - dim(E cap im ev_m)`; these are equal. QED

Consequently `TC(A,m) = {u : l_u^m in J(A)_m}` (zetc Lemma 1.1): `Phi_c(u) = sum_k c_k l_u(k)^m` is the pairing of `c` with
`ev(l_u^m)`, and `ev(l_u^m) perp V(A,m) = E^perp` iff `ev(l_u^m) in E` iff `l_u^m in J(A)_m`.

### 2.2 The inclusion

For `v in R^d` let `tau_A(v) = #{<k,v> : k in A} - 1` (so `tau_A(v) <= h_A(v)` for integer `v`, with equality iff the
values have no gaps). Let `J_tau` be the ideal generated by `l_v^{tau_A(v)+1}`, `v in Z^d`, and `J_w subset J_tau` the one
generated by `l_v^{h_A(v)+1}` (the judge's ideal). Define the linear system of degree-`m` forms

```text
L_m = { Phi in Sym^m(u) : D_v^{tau_A(v)+1} Phi = 0 for every v in Z^d with tau_A(v) < m } .
```

**Proposition 2.1.** (i) `J_w subset J_tau subset J(A)`. (ii) `Phi in L_m` iff for every integer `v` with `tau_A(v) < m`
the form `Phi` has multiplicity `>= m - tau_A(v)` at the point `[v] in P^{d-1}` (the point of `u`-space with
coordinates `v`). (iii) `span{Phi_c : c in V(A,m)} subset L_m` and hence

```text
TC_C(A,m)  ⊇  Bs_C(L_m),       TC(A,m) = TC_C(A,m) ∩ R^d  ⊇  Bs_R(L_m) = Bs_C(L_m) ∩ R^d .
```

`L_m` is defined over `Q` (finitely many of the rational conditions already cut it out), so `Bs(L_m)` is a
`Q`-variety. The same holds with `h_A(v)` in place of `tau_A(v)` (a weaker system, larger `L_m`, smaller base locus).

*Proof.* (i) `prod_{t in <A,v>} (l_v - t)` vanishes on `A` and has leading form `l_v^{tau+1}`. (ii) For any `w`,
`Phi(w + t v) = sum_i (t^i/i!) (D_v^i Phi)(w)`; `D_v^{tau+1}Phi = 0` iff this has `t`-degree `<= tau` for all `w`. Write
`Phi(v + w) = sum_i Phi_i(w)` with `Phi_i` homogeneous of degree `m - i` in `w` (homogeneity: `Phi(t v + w)` has the
term `t^i Phi_i(w)`); `t`-degree `<= tau` iff `Phi_i = 0` for `i > tau` iff the local expansion at `[v]` starts in degree
`>= m - tau`, i.e. multiplicity `>= m - tau`. (iii) `D_v^{tau+1}Phi = 0` iff `(l_v^{tau+1} g)(partial)Phi = 0` for all `g in R_{m-tau-1}`
(nondegeneracy of apolarity in degree `m - tau - 1`), so `L_m = (J_tau)_m^perp superset J(A)_m^perp = span{Phi_c}` by (i) and
Lemma 2.0. A common zero of `L_m` is a common zero of every `Phi_c`. Real points: `Phi_c` and a rational basis of `L_m`
are real forms, and the real loci are the real points of the complex ones. QED

**Remarks.** (a) The README's "every `[v]` in `P^(d-1)`" must be read as "every **integer** `[v]` with
`h_A(v) < m`"; for directions with `h_A(v) >= m` the condition is vacuous. (One may also allow real `v` with
`tau_A(v) < m`; the proof is the same.) (b) The inclusion is in general strict:
in `M8`, `(J_w)_8 = 2940 < 2961 = J(A)_8` (§4), so `span{Phi_c} = J(A)_8^perp` is a proper subspace of `L_8` and
`TC(A,8)` may be strictly larger than `Bs(L_8)`. Equality `J(A)_m = (J_tau)_m` is exactly `Loc(1)` in degree `m`.

### 2.3 Linear subspaces of the base locus

**Proposition 2.2 (cheap spans lie in `Bs(L_m)`).** If `v_1..v_s in Z^d` and `sum_i tau_A(v_i) < m` (in particular if
`sum h_A(v_i) < m`, i.e. the family is cheap), then `Lambda = span(v_1..v_s) subset Bs_C(L_m)`.

*Proof.* For `Phi in L_m` put `psi(t_1..t_s) = Phi(sum t_i v_i)`, a form of degree `m`. `partial_{t_i} psi = (D_{v_i}Phi)(...)`, so
`partial_{t_i}^{tau_i+1} psi = 0`, i.e. `deg_{t_i} psi <= tau_i` for each `i`. Every monomial of degree `m` has some
exponent `alpha_i > tau_i` because `sum tau_i < m`. So `psi = 0`. (No independence of the `v_i` is needed.) QED

In the language of fat-point systems (Brambilla-Dumitrescu-Postinghel, "On a notion of speciality of linear systems
in P^n", Trans. AMS 367 (2015)), the *linear cycle* `L_I` spanned by a subset `I` of the (finitely many relevant)
fat points `[v_i]^{mu_i}` is contained in the base locus with multiplicity `max(0, sum_{i in I} mu_i - (|I|-1)m)`; with
`mu_i = m - h_i` this is positive iff `sum_{i in I} h_i < m`. So **"cheap spans = linear cycles `L_I` of positive
multiplicity" is a correct, but definitional, statement.** It does not say that every linear subspace of
`Bs(L_m)` is a cheap span, and that is false:

**Proposition 2.3 (the Fano spans of `M8`).** In `M8` (`Syl(8)`, `L = 2`, `m = 8`; definitions in §4.1) let `l = {a, b, a+b}`
be a line of the Fano plane and `U_l = span(c_a, c_b, c_{a+b})` (3-dimensional). Then
(i) every nonzero integer `v in U_l` has `h_A(v) >= 4`, so no family of integer directions spanning `U_l` (or any
3-dimensional subspace of it) is cheap at `m = 8`; (ii) `U_l subset Bs_C(L_8)`, hence `U_l subset TC(A, 8)`;
(iii) `rho(pi_{U_l}(A)) = 6 < 8` (exact, §4), so `U_l` is projection-cheap with transversal order `>= 2`.

*Proof.* (i) `U_l` is spanned by the indicator vectors `p_1 = 1_{T_a cap T_b}`, `p_2 = 1_{T_a \ T_b}`, `p_3 = 1_{T_b \ T_a}` of
three disjoint 2-sets (`c_a = p_1 + p_2`, `c_b = p_1 + p_3`, `c_{a+b} = p_2 + p_3`). A nonzero `v in U_l` is constant `= alpha != 0`
on some 2-set `{x, y}`; `±(e_x + e_y) in A` (all `k'` with `|k'|_1 <= 2` lie in `A`), so `h_A(v) >= 4|alpha| >= 4`.
(ii) Let `Phi in L_8` and `psi(t) = Phi(t_1 p_1 + t_2 p_2 + t_3 p_3)`, a ternary octic, i.e. a plane octic on `P(U_l)`. By
Prop 2.1(ii) (restricted to `U_l`) it has multiplicity `>= 8 - h` at each integer point of `P(U_l)` with `h < 8`. Exact
widths: `h(p_i) = 4` and `h(p_i + p_j) = 4` (these are `c_a, c_b, c_{a+b}`), `h(p_1 + p_2 + p_3) = 6` (`p_1+p_2+p_3 = z_l`,
the indicator of six non-base points, in the width set of §4). So `psi` has multiplicity `>= 4` at
`P_i = [p_i]`, `>= 4` at `Q_ij = [p_i + p_j]`, `>= 2` at `Z = [p_1+p_2+p_3]`. The six lines
`P_1P_2Q_12, P_1P_3Q_13, P_2P_3Q_23, P_1Q_23Z, P_2Q_13Z, P_3Q_12Z` are rational lines each containing three of these points.
Peel them off in this order, using: *if a line contains points whose current multiplicities sum to more than the
current degree, the line divides the curve, and dividing lowers the multiplicity by 1 exactly at the points on that
line.* The sums are `12 > 8, 11 > 7, 10 > 6, 7 > 5, 6 > 4, 5 > 3`. What remains is a conic with multiplicity `>= 2` at the
three non-collinear points `Q_12, Q_13, Q_23`, which is `0` (the line `Q_12Q_13` carries `2 + 2 > 2`, divide; the remaining
linear form has multiplicity `2` at `Q_23`). Hence `psi = 0`. (Exact check, `fano_baselocus_out.txt`: the ternary octics
with these multiplicities at the 7 points, and also with all integer points of `U_l` with coefficients in `[-2,2]`, form
the zero space; the widths of the 7 points are as stated.) (iii) is computed in §4. QED

So the base locus of the width fat-point system already contains rational subspaces that are not cheap spans, and
in `M8` it **sees the projection-cheap Fano spans**, which are exactly where the non-width generators of `J(A)`
live (§4). Linear base components of fat-point systems through *three or more collinear points* (as here) are
classical (they are the components that BDP's `L_I` with `|I| <= dim + 1` independent points do not list).

**Correction for README §4.2(a):** replace "The cheap spans are exactly its linear base cycles" by "The cheap
spans are exactly the linear cycles `L_I` of positive multiplicity; the base locus contains further rational
subspaces (e.g. the Fano spans of `M8`, through collinear triples of fat points), and possibly non-linear
components."

## 3. J6(b) — the cut ideal and the strong Lefschetz property

**Verdict: CORRECT. The quoted theorem is right and is now proved here (self-contained `sl_2` proof, characteristic
0); the Hilbert-function/dimension argument is correct; it extends verbatim to arbitrary signs and to any number of
disjoint cuts (Theorem 3.3). Exact checks: `refJ_checks/slp_check.py` (`slp_check_out.txt`).**

### 3.1 The quoted theorem

**Reference statement.** R. P. Stanley, *Weyl groups, the hard Lefschetz theorem, and the Sperner property*, SIAM J.
Algebraic Discrete Methods 1 (1980) 168-184: for a field `k` of characteristic 0 and `B = k[x_1..x_s]/(x_1^{a_1},...,x_s^{a_s})`,
the linear form `l = x_1 + ... + x_s` is a strong Lefschetz element, i.e. for all `i, t >= 0` the map
multiplication by `l^t : B_i -> B_{i+t}` has maximal rank (Stanley derives it from hard Lefschetz on `P^{a_1-1} x ... x P^{a_s-1}`, whose
cohomology ring is `B` with degrees doubled). J. Watanabe (*The Dilworth number of Artinian rings and finite posets with
rank function*, Adv. Stud. Pure Math. 11 (1987)) gives the representation-theoretic proof; see also Harima-Maeno-Morita-
Numata-Wachi-Watanabe, *The Lefschetz Properties*, Springer LNM 2080 (2013), Ch. 3 and 4. The quotation in README §4.2(b)
("Stanley 1980 for monomial complete intersections, char 0") is accurate. Here only `a_i = 3` is needed.

**Theorem 3.1 (SLP for `a_i = 3`, self-contained).** Let `char k = 0`, `B = k[x_1..x_s]/(x_1^3..x_s^3)`, `eps in {±1}^s` and
`l = sum eps_i x_i`. With `h_i = dim B_i` (coefficients of `(1+t+t^2)^s`), multiplication by `l^t : B_i -> B_{i+t}` has rank
`min(h_i, h_{i+t})` for all `i, t >= 0`.

*Proof.* (Signs.) `x_i -> eps_i x_i` is a graded automorphism of `B` mapping `sum x_i` to `l`, so take `eps = 1`.
(`sl_2`-structure.) On `V = k[x]/(x^3)` with basis `1, x, x^2` put `E = (x .)`, `H = diag(-2, 0, 2)`, and `F(1) = 0`,
`F(x) = 2`, `F(x^2) = 2x`. Direct check: `[H,E] = 2E`, `[H,F] = -2F`, `[E,F] = H` (e.g. `[E,F]x = E(2) - F(x^2) = 2x - 2x = 0 = Hx`,
`[E,F]1 = -F(x) = -2`, `[E,F]x^2 = E(2x) = 2x^2`). So `V` is the 3-dimensional irreducible `sl_2`-module. On
`B = V^{⊗ s}` the operators `E_tot = sum_i E^{(i)}`, `H_tot = sum_i H^{(i)}`, `F_tot = sum_i F^{(i)}` again satisfy the
`sl_2` relations, `E_tot` is multiplication by `x_1 + ... + x_s`, and `B_i` is the `H_tot`-eigenspace of weight `2i - 2s`.
(Weights.) In characteristic 0 the finite-dimensional module `B` is a direct sum of irreducibles `V(N)` (weights
`-N, -N+2, ..., N`, each with a 1-dimensional weight space, and `E` maps the weight-`mu` line isomorphically onto the
weight-`mu+2` line for `mu < N`). Put `lambda = 2i - 2s`, `mu = lambda + 2t`. If `|mu| <= |lambda|` (then `lambda <= 0`), every
`V(N)` with a nonzero weight-`lambda` component has `N >= |lambda| >= mu`, so `E^t` is injective on `B_i`. If `|mu| >= |lambda|`,
every `V(N)` with a weight-`mu` component has `N >= |mu| >= |lambda|`, so it has a weight-`lambda` component mapping onto its
weight-`mu` component and `E^t` is surjective onto `B_{i+t}`. One of the two cases always holds. QED

(`slp_check_out.txt`: maximal rank of `l^t : B_i -> B_{i+t}` verified exactly for all `i, t`, `s <= 6`, three sign vectors.)

### 3.2 Hilbert function of the truncation

**Lemma 3.2.** Let `h = (h_0..h_{2s})`, `N = 2L + 1`, `0 <= L <= s - 1`, `s >= 2`. Then
`dim (B/l^N B)_i = max(h_i - h_{i-N}, 0)` (`h_j = 0` for `j < 0`), the sum over `i` is `sum_{|i-s| <= L} h_i
= #{k in {-1,0,1}^s : |sum eps_i k_i| <= L}`, and the top nonzero degree is exactly `s + L`.

*Proof.* `dim (l^N B)_i = rank(l^N : B_{i-N} -> B_i) = min(h_{i-N}, h_i)` by Theorem 3.1. `h` is symmetric about `s` and
weakly increasing on `[0, s]`, so `h_a >= h_b` whenever `|a - s| <= |b - s|`. For `i <= s + L` one has
`|i - s| <= |i - N - s|` (if `i <= s` trivially; if `s < i <= s + L` because `2(i - s) <= 2L < N`), so `h_i >= h_{i-N}`; for
`i >= s + L + 1` the reverse inequality holds (if `i - N >= s` trivially; otherwise `2(i-s) >= N`). Telescoping,
`sum_i max(...) = sum_{i <= s+L} (h_i - h_{i-N}) = sum_{i=s+L-N+1}^{s+L} h_i = sum_{|i-s| <= L} h_i`, and `h_i` counts
`k in {0,1,2}^s` with sum `i`, i.e. `k - 1 in {-1,0,1}^s` with sum `i - s`; signs do not change the count
(`k_i -> eps_i k_i`). Top degree: at `i = s + L`, `h_{s+L} = h_{s-L} > h_{s-L-1} = h_{s+L-N}` by **strict** increase of
trinomial coefficients on `[0, s]` for `s >= 2` (induction: `h^{(s)}_{i+1} - h^{(s)}_i = h^{(s-1)}_{i+1} - h^{(s-1)}_{i-2} > 0` for
`0 <= i < s`, using the strict statement for `s - 1 >= 2`, symmetry, and the base `(1,2,3,2,1)`; also checked for `s <= 200`). QED

### 3.3 One cut, signs, disjoint cuts

**Theorem 3.3.** Let `S_1, ..., S_q subset [d]` be pairwise disjoint, `eps^{(j)} in {±1}^{S_j}`, `0 <= L_j`, and
`A = {k in {-1,0,1}^d : |sum_{i in S_j} eps^{(j)}_i k_i| <= L_j, j = 1..q}`. Then

```text
J(A) = ( x_1^3, ..., x_d^3,  l_1^{2L_1+1}, ..., l_q^{2L_q+1} ),      l_j = sum_{i in S_j} eps^{(j)}_i x_i ,
```

in all degrees, and the Hilbert series of `R/J(A)` is `prod_j T_{s_j, L_j}(t) * (1+t+t^2)^{d - sum s_j}`, where
`T_{s,L}(t) = sum_i max(h^{(s)}_i - h^{(s)}_{i-2L-1}, 0) t^i` (`T_{s,L} = (1+t+t^2)^s` if `L >= s`).

*Proof.* Let `J_0` be the right side. `J_0 subset J(A)`: `x_i^3` is the leading form of `x_i^3 - x_i in I(A)`, and
`l_j^{2L_j+1}` the leading form of `prod_{|c| <= L_j}(l_j - c) in I(A)` (Prop 2.1(i)). The generators of `J_0` live in
disjoint sets of variables, so `R/J_0 = ⊗_j k[x_{S_j}]/(x^3, l_j^{N_j}) ⊗ k[x_rest]/(x^3)` and its Hilbert series is the
product; by Lemma 3.2 (Theorem 3.1 for each block) its total dimension is `prod_j #{k_{S_j} : |l_j(k)| <= L_j} * 3^{rest} = |A|`.
For finite `A`, `dim_k R/J(A) = dim_k R/I(A) = |A|` (the affine Hilbert function `dim R_{<=i}/I(A)_{<=i}` equals
`sum_{t<=i} dim (R/J(A))_t`, and is `|A|` for large `i`). Two ideals `J_0 subset J(A)` with quotients of the same finite
dimension are equal. QED

**Corollaries.** (a) One cut (`q = 1`, `d = s`, `0 <= L <= s-1`): `rho(A) = s + L` — the interpolation degree is the top
degree of `R/J(A)` (Lemma 3.2), which is zetc Thm 2.1's `rho` (a one-line proof, as the judge says). With `d - s` free
coordinates, `rho = s + L + 2(d - s)`. (b) Disjoint cuts: `rho(A) = sum_j (s_j + min(L_j, s_j)) + 2(d - sum s_j)`. (c) For
disjoint cuts `J(A) = J_w`-type: `J(A)` is generated by width powers of `b_0 = 1` directions (`Loc(1)` holds in all degrees).
(d) By Lemma 2.0, `span{Phi_c : c in V(A,m)}` is the space of degree-`m` forms `Phi` with `partial_{u_i}^3 Phi = 0` for all `i`
(i.e. multilinear-of-degree-`<= 2` in each variable) and `(sum_{i in S_j} eps_i partial_i)^{2L_j+1} Phi = 0`; and
`TC(A,m) = {u : l_u^m in J_0}`, which is the "explicit ideal" form of zetc Thm 2.1.

**Exact checks** (`slp_check_out.txt`): for 7 boxes (one cut with signs, `d = 5,6,7`; two and three disjoint cuts with
mixed signs, `d = 6, 7`), the Hilbert function of `R/J(A)` computed from exact evaluation ranks equals the product of
Lefschetz truncations in every degree. This is an independent re-check of `judge_checks/onecut_slp_out.txt` (which
covers `eps = 1`, one cut) and of the extension.

**Limits.** The argument needs disjoint supports (tensor decomposition) and `±1` coefficients (so that the values of
`l_j` on the cube are all consecutive integers in `[-s_j, s_j]` and the count matches). Overlapping cuts are not covered;
`M8` shows that there `J(A)` is **not** generated by the cut powers and `x_i^3` (§4: `J_w` already contains more, and is
still too small).

## 4. J6(c) — the M8 locality computation

**Verdict: all numbers CONFIRMED exactly over `Q` with independent code. The "3-column" statement is sharpened: only
the 7 Fano-line triples matter, for a structural reason (a parity relation among the three columns), and the 21 column
pairs and 28 non-line triples add nothing *for a structural reason* (their images are full grids). "b_0 = 2 fails"
is now exact over `Q` for column pairs (in all degrees `<= 8`), and also exact over `Q` for a large family of general
integer pairs (51,933 planes), but it is NOT proved for all integer pairs; and general integer pairs are not useless:
planes through the Fano "half-sum" direction supply 14 of the 21 missing degree-6 dimensions.** Library
`refJ_checks/m8lib.py`; drivers `m8_step1..9.py`, outputs `m8_step*_out.txt`.

### 4.1 The model

`M8 = Syl(8)`, `L = 2`, `m = 8`: points `F_2^3 = {0..7}`, base `0` (the judge used base `7`; translations of `F_2^3`
permute the sign columns up to sign, so the boxes are isomorphic), `k' in {-1,0,1}^7`, `k_base = -sum k'`. On this set
`<k, beta_a> = 2 k(S_a) = -2 k'(T_a)`, `T_a = {x != 0 : a.x = 1}` (4 points), so `A = {k' : |k'(T_a)| <= 2, a = 1..7}` and
the column directions are `c_a = 1_{T_a}`. All `k'` with `|k'|_1 <= 2` lie in `A` (checked).

### 4.2 Exact recomputation

| quantity | value | method |
| --- | --- | --- |
| `|A|` | 1037 | enumeration |
| `H(R/J(A))_i`, `i = 0..8` | `[1, 7, 28, 77, 161, 259, 287, 175, 42]` | exact `fmpz` ranks of cumulative evaluation matrices (monomials with exponents `<= 2` suffice on the cube) |
| `dim J(A)_i` | `[0, 0, 0, 7, 49, 203, 637, 1541, 2961]` | `binom(i+6,6) - H_i`, exact |
| width directions with `tau + 1 <= 8` | 358 primitive `v` (7 with `h = 2`, 49 with `h = 4`, 302 with `h = 6`); `tau = h` for all; `max |v_i| = 2` | enumeration of `|v_i| <= 3`; complete (below) |
| `dim (J_w)_i` | `[0, 0, 0, 7, 49, 203, 616, 1485, 2940]` | **exact over `Q` in every degree `<= 8`** (`fmpz` rank of all generator multiples; degree 8: 9464 x 3003) |

*Completeness of the width set.* `±e_x, ±e_x ± e_y in A`, so `<A, v>` contains `0, ±v_x, ±v_y, ±v_x ± v_y`. If `v` is primitive
with `M = max|v_x| >= 4`, some entry `b = |v_y|` has `0 < b < M` and `M != 2b` (otherwise all entries lie in `{0, ±M, ±M/2}` and
`gcd >= 2`), and then those values are 9 distinct numbers, `tau >= 8`. So `|v_i| <= 3` captures every generator of degree `<= 8`
(in `J_tau` as well as in `J_w`).

### 4.3 Which projections matter, and why

| family | image `Y = pi_F(A)` | `rho(Y)` | `J(Y)` | pullback adds to `J_w` (degrees `<= 8`) |
| --- | --- | --- | --- | --- |
| 21 column pairs | full grid `[-2,2]^2` (25 points) | 8 (`= m`, not admissible in `Loc`) | `(y_1^5, y_2^5)` | **nothing** (pullbacks are `c_a^5`, width powers) |
| 28 non-line column triples | full grid `[-2,2]^3` (125) | 12 | `(y_i^5)` | **nothing** |
| 7 Fano-line triples `{a, b, a+b}` | 57 points of `P_2` (the 63-point parity box minus the 6 points `(±2,±2,±2)` whose signs are not all equal) | **6** | minimal generators: 6 quintics, 3 sextics | **everything**: `J_w + (7 Fano pullbacks) = J(A)` in all degrees `<= 8` |

The grids are products of 1-dimensional sets, so `J(Y)` is generated by the coordinate powers (Theorem 3.3 with no
cut); this makes "pairs add nothing" a theorem, not a mod-`p` observation. The Fano image is constrained by
`y_a + y_b + y_{a+b} = 2 k(S_a cap S_b)` (Lemma 1.4), a **parity relation**: in `u`-space the half-sum
`z_l = (c_a + c_b + c_{a+b})/2 = 1_{T_a ∪ T_b}` (six non-base points) is an integer direction with `h = 6`.

Exact results (degree-6 dimensions exact over `Q`; degrees 7, 8 by the valid sandwich
`dim_p(J_loc)_i <= dim_Q(J_loc)_i <= dim_Q J(A)_i` with the right side exact):

| ideal | dims, degrees 0..8 | equal to `J(A)`? |
| --- | --- | --- |
| `J_w` | `[0,0,0,7,49,203,616,1485,2940]` | no, from degree 6 (exact) |
| `J_w` + 21 column-pair pullbacks | same as `J_w` | no (exact, structural) |
| `J_w` + 35 column-triple pullbacks | `[0,0,0,7,49,203,637,1541,2961]` | **yes**, all degrees `<= 8` |
| `J_w` + 7 Fano-line triples only | `[0,0,0,7,49,203,637,1541,2961]` | **yes** |
| `J_w` + 28 non-line triples only | `[0,0,0,7,49,203,616,1485,2940]` | no |

The judge's `[0,0,0,7,49,203,637,1541,2961]`, `616 < 637` and "3 columns suffice" are confirmed. The 21-dimensional
degree-6 gap is `7 lines x 3` (each Fano triple contributes 3 sextic generators).

### 4.4 Does `b_0 = 2` fail? (general integer directions, as `Loc` is stated)

`Loc(b_0)` allows families of **arbitrary integer directions**, not only columns. For a rational plane `P` the
pullback ideal `J_P = pi_P^* J(pi_P(A))` depends only on `P` (a change of basis acts on the image by an invertible linear
map, and leading-form ideals are equivariant), and in degree `<= 6` only `(J_P)_6` can matter (because
`J(A)_e = (J_w)_e` for `e <= 5`). Results (degree 6, exact over `Q`):

| planes tested | number | contributing | `dim (J_w + sum J_P)_6` |
| --- | --- | --- | --- |
| planes inside the 7 Fano 3-spaces `U_l`, normals with entries `<= 2`, `3`, `5` (basis `c_a, c_b, z_l`) | 343 / 1015 / 4039 | 56 / 112 / 280 sextic generators | **630** in all three cases |
| all planes spanned by two of the 358 low-width directions | 51,933 | 3,213 | **630** |
| (for comparison) `J(A)_6` | | | 637 |

The contributing planes inside `U_l` are exactly the pencil through `z_l`; each contributes one sextic, and the pencil
contributes 2 dimensions per line in total. The **admissible** ones (`rho(pi_P(A)) < 8`) are `span(z_l, c_a)`,
`span(z_l, c_b)`, `span(z_l, c_{a+b})` (`|Y| = 21`, `rho = 6`), and they alone already give the full gain of 14
(`m8_step9_out.txt`).

**Conclusions.**
* Exact over `Q`: `b_0 = 2` fails for column pairs (in degrees 6, 7, 8), and for the 51,933 planes above (degree 6:
  `630 < 637`).
* Not proved: that no rational plane at all fills the remaining 7 dimensions. A finite reduction would need a
  classification of the planes `P` for which `pi_P(A)` lies on a curve of degree `<= 6` (only those contribute in degree
  `<= 6`); with the admissibility condition `rho < 8` one has `|pi_P(A)| <= 36`, which is a strong restriction, but I have not
  completed that classification. So "b_0 = 2 fails in M8" is **exact for the tested families, evidence in general**.
* The judge's reading "pairs add nothing" is an artefact of using only **columns**: general integer pairs through the
  half-sum direction `z_l` add 2/3 of the gap. The irreducible part (one sextic per Fano line) needs the whole 3-space.
* Degrees `7, 8` for general planes were not tested.

### 4.5 What `M8` says about `Loc` in the regime

The non-width part of `J(M8)` comes from a mechanism that exists in **every** Sylvester box: Fano triples of columns
satisfy a parity relation, their image is a sublattice set (`P_j`), and its projection cost is `3j` instead of `6j`
(Lemma 1.3). At `m' = ceil(4w) ~ 8j` these blocks are admissible (`3j < m'`), as are products of two of them with a
column (`8j < m'`). So in the regime the analogue of `M8`'s locality *suggests* (not proved: no regime `J(A)` is
computable) `b_0 >= 3`, and possibly more: XOR-closed blocks of 7 columns have `rho <= 3, 8` at `j = 1, 2` (rigorous
upper bounds via mod-`p` ranks). For README §7.3 item 3: `rho(pi_F(A)) >= 2j(|F| - O(1))` fails with a deficit of `3j`
per Fano block contained in `F` (`|F| = 3s`, `rho <= 3js`), so any such bound must count the XOR relations inside `F`;
projection-cheap families are not captured by width combinatorics, and they are where the low-degree non-width
generators come from.

## 5. J7 — the hierarchy and the empty-gap-region remark

**Verdict: (a) CORRECT, with one missing qualifier (`MIN* => (MIN)` needs `C <= 1`, as in MIN\*'s own definition) and a
compactness remark. (b) CORRECT as stated; STRENGTHENED: on Sylvester boxes at the large orders of nearstrata Prop 5.2,
MIN^pair fails for every pair-data set with non-empty gap region, except in a window of width `< 2/(n-1)` for
`max eps_xy` (Prop 5.2 below). So the referee's witness is the vacuous half of a dichotomy, and MIN^pair gives no escape
from that obstruction on data carried by actual clusters or gap-admissible fakes.**

### 5.1 The hierarchy

Fix the cluster's data (profile, primes, `W`, `C`), `abar = W/(2d) + 2 ln C'`, `eps_xy = d_xy - W/2 + 2 ln C' >= 0`,
`sum eps_xy <= N_2 abar`, an order `m''` in the window, `eps_0 = C d e^{-W/4}/(m''+1)`, and `S = {u in R^n : u_base = 0,
|u|_inf = 1}`. Sets (all subsets of `S`):
* `G = {u : |u_x - u_y| >= e^{-eps_xy/2} for all x < y}` (gap region);
* `Gam = {u : Gamma(u) <= abar}`;
* `R_per = {u : dist_inf(u, Delta_pi) >= sigma_pi for all partitions pi with P(pi) >= 1}`,
  `sigma_pi = (1/2) e^{-N_2 abar/(2P(pi))}`;
* `R_unif(Q_0) = {u : dist_inf(u, Delta_pi) >= (1/2) e^{-Q_0 abar/2} for all pi with P(pi) >= N_2/Q_0}`.

**Proposition 5.1.** `G ⊂ Gam ⊂ R_per ⊂ R_unif(Q_0)` for every `Q_0 > 0`. Hence, with each hypothesis read as
"`kappa_H > eps_0` on the set" (MIN and MIN^pair as minima over the compact sets `Gam`, `G`; empty minimum `= +infinity`):

```text
ZE-TC^coll_eps(Q_0)  =>  ZE-TC^coll_sigma  =>  (MIN) at level abar  =>  MIN^pair ,
MIN*  =>  (MIN)   whenever abar >= 1.5 ln n - 7 and C <= 1 .
```

*Proof.* `G ⊂ Gam`: on `G`, `ln^+(1/|u_x - u_y|) <= eps_xy/2` (as `eps_xy >= 0`), so `Gamma(u) <= (1/N_2) sum eps_xy <= abar`.
`Gam ⊂ R_per`: nearstrata Lemma 3.3(ii) with `g = abar`, applied to `L = Delta_pi` (`P(Delta_pi) = P(pi)`; the lemma's
proof uses only collision pairs of `L` and `u_base = 0`). `R_per ⊂ R_unif(Q_0)`: with the convention of §1.1 the uniform
family `{P(pi) >= max(1, N_2/Q_0)}` is a subfamily of the per-subspace family, and on it
`sigma_pi = (1/2)e^{-N_2 abar/(2P(pi))} >= (1/2)e^{-Q_0 abar/2}`; so a point meeting all per-subspace radii meets the
uniform ones. Implications: a pointwise statement on a larger set implies it on a smaller one;
for MIN and MIN^pair note that `Gamma` is lower semicontinuous (`ln^+(1/|t|) -> +infinity` as `t -> 0`), so `Gam` and `G` are
compact and the minima of the continuous `kappa_H` are attained; pointwise `> eps_0` on the set then gives `min > eps_0`.
`MIN* => (MIN)`: MIN\* gives `Lambda_A(abar) > (d/(m'+1)) e^{-d abar/2} = (d/(m'+1)) e^{-W/4} C'^{-2d}`; this is
`>= eps_0 = C d e^{-W/4}/(m'+1)` iff `C C'^{2d} <= 1`, i.e. iff `C <= 1`. QED

(Without the convention `P(pi) >= 1`, `Q_0 > N_2` would put the discrete partition in the uniform family, make
`R_unif = ∅`, and break the inclusion; so the convention is needed for the first implication as stated.)

All four hypotheses suffice for t1rigour Thm 6.8 / nearstrata Thm 4.1 (each contains `hat u`: `hat u in G` by the pair-gap
lemma, which is the weakest membership). J5 refutes the first for `n = 2^12..2^17`, extended here to `2^12..2^20`; §1.6
refutes the second for `n = 2^12, 2^13, 2^14` at near-extremal `W`. All four are open asymptotically, and (MIN), MIN^pair
are open at every regime `n`. The README table §3.3 entry "ZE-TC^coll_eps(Q_0) => ZE-TC^coll_sigma => (MIN) => MIN^pair; MIN* => (MIN)" should add
"(`C <= 1`)" to the last implication.

### 5.2 The empty gap region, and a dichotomy at Sylvester large orders

The referee's witness (referee_nearstrata §2.3: Sylvester `k = 16`, genuine primes, `eps_xy in [12.4997, 12.8739]`) has
`(n-1) e^{-max eps/2} = 65535 * e^{-6.43695} = 104.94 > 2` (recomputed). Since every pair of `G` needs a gap
`>= e^{-max eps/2}`, and `n` points of `[-1,1]` (base included) have a pair closer than `2/(n-1)`, `G = ∅`. An actual cluster
has `hat u in G` (pair-gap lemma, nearstrata Lemma 3.1 / t1rigour Lemma 6.6 proof), and a gap-admissible fake has its
direction in `G` by definition; so neither carries these data. The judge's reading is correct: the example separates
MIN^pair from (MIN) only as hypotheses on arbitrary data.

**Proposition 5.2 (dichotomy; proved).** Let `n = 2^k`, the Sylvester profile, and an order `m''` in `(2jk, m_max]`
(non-empty iff `k >= 16`, nearstrata Prop 5.2). For any pair data `eps_xy >= 0`:
1. if `max eps_xy < 2 ln((n-1)/2)`, then `G = ∅` and MIN^pair holds vacuously;
2. if `max eps_xy >= 2 ln(n/2)`, then MIN^pair **fails** (the minimum over `G` is `0`).

*Proof.* (1) as above. (2) Let `b_i(x) = 1[x_i != x_{base,i}]` and `u_x = -b_1(x) + sum_{i=2}^k 2^{-(i-1)} b_i(x)`. Then
`u_base = 0`, `u` lies in the span `L` of the `k` coordinate columns (in `u`-coordinates), its `n` values are
`-1 + t 2^{-(k-1)}` (`b_1 = 1`) and `t 2^{-(k-1)}` (`b_1 = 0`), `t = 0..n/2-1`, all distinct with consecutive gaps exactly
`2^{-(k-1)} = 2/n`, and `|u|_inf = 1`. Every pair has `|u_x - u_y| >= 2/n >= e^{-max eps/2} >= e^{-eps_xy/2}`, so `u in G`. By
nearstrata Prop 5.2, `L subset TC_H(A, m'')`, so `kappa_H(u) = 0 <= eps_0`. QED

Since `2 ln(n/2) - 2 ln((n-1)/2) = 2 ln(n/(n-1)) < 2/(n-1)`, MIN^pair can hold on a non-empty gap region at these orders
only if `max eps_xy` lies in a window of width `< 2/(n-1)` (`= 3.05e-5` at `k = 16`; endpoints `20.794385`, `20.794415`).
In particular for **every actual cluster** with a Sylvester profile, `k >= 16`, MIN^pair fails at every order in
`(2jk, m_max]` unless `max eps_xy` falls in that window. (Within the window `G` needs nearly equally spaced points with
gap `> 2/n`; `L`-points with distinct values have gap `<= 2/n` because their values are subset sums of `k` reals in a range
of length `<= 2` containing `0`, so the argument above cannot be pushed into the window; this is not needed.)

**Consequence for README §4.3(b).** "Whether it is weaker on non-empty gap regions is open" should be sharpened: at the
large Sylvester orders (the only regime obstruction to (MIN) known), MIN^pair is *not* weaker on non-empty gap regions,
up to the width-`2/(n-1)` window. At `m' = ceil(4w)` (the order route 5 uses) the question remains open: every linear
subspace of `TC` known at that order has collision pairs (at most 32 classes, §1.6), so none meets `G` (whose points have
all coordinates distinct), and MIN^pair can only fail through near-strata or irrational points (§6).

## 6. Consequences for `Loc(b_0)` and for route 5

1. **Proved instances of locality.** `Loc(1)` (generation by width powers) holds in all degrees for one cut, for any
   number of disjoint cuts, with arbitrary signs (Theorem 3.3), via the strong Lefschetz property (Theorem 3.1, proved
   here). This is the whole content of J6(b), now without a quoted theorem.
2. **`M8`.** `Loc(3)` holds exactly in all degrees `<= m = 8`, with the 7 Fano-line triples (admissible: `rho = 6 < 8`) as the
   only non-width families needed. `Loc(2)` fails exactly over `Q` for column pairs (all degrees 6-8) and for 51,933
   general planes (degree 6); for all rational planes it is open (evidence: fails). General integer pairs are not
   useless (14 of 21 missing dimensions).
3. **The mechanism is arithmetic, not combinatorial.** The non-width generators of `J(M8)` come from a *lattice relation*
   among columns (`c_a + c_b + c_{a+b} in 2 Z^d`): the image of a Fano triple is a sublattice set whose interpolation
   degree is half that of the full box. The same relation holds for every Sylvester box, where Fano blocks cost `3j`
   (Lemma 1.3, proved for all `j`) instead of `6j`. So in the regime:
   * the tangent cone `TC(A, ceil(4w))` contains projection-cheap rational subspaces with up to 32 classes (`L_{2,1}`)
     — twice the 16 classes allowed for width-cheap spans (t1rigour Cor 6.4);
   * these subspaces refute the uniform-`Q_0` hypothesis up to `n = 2^20` and the per-subspace hypothesis up to `n = 2^14`
     at near-extremal `W` (§1.6);
   * any `Loc(b_0)` for Sylvester-type profiles must allow families with lattice relations (by analogy with `M8` this
     suggests `b_0 >= 3`; not proved in the regime), and the cost bound of README §7.3 item 3 must subtract `3j` per Fano
     block.
4. **The strata half, made precise.** By §1.2/§1.5 the per-subspace hypothesis forces `o_L N_2/P(L) <= d` for every
   linear `L subset TC_H` (at the levels used). For XOR-closed families of dimension `t` (2^t classes) this reads
   `(m' - rho_t) 2^t <~ d` with `rho_t = rho(Y_t(j))`, the interpolation degree of the image box
   `Y_t(j) = {kappa in Z^{2^t}_0 : |kappa(S_a)| <= j}`. Known: `rho_1 = 2j`, `rho_2 <= 3j` (proved), `rho_3 <= 3, 8` at
   `j = 1, 2`. **If `(m' - rho_t) 2^t > (1 + delta) d` for some `t <= k` along a sequence `n -> infinity`, the per-subspace
   hypothesis (and with it this form of route 5) fails at all large levels, asymptotically** (Lemma 1.2 with per-subspace
   radii; the span of a `t`-dimensional XOR block has `2^t` equal classes). Since `j ~ sqrt(2 n ln n)`, this needs
   `rho_t <= m' - Omega(n/2^t)`, i.e. `rho_t/j` must stay below `8` up to `2^t ~ sqrt(n/ln n)`. Deciding the growth of
   `rho_t/j` in `t` is therefore the first concrete
   regime question that `Loc` raises; it is a question about the lattice points of the dilated Sylvester polytope
   `{kappa : |hat kappa(a)| <= 2j}` only, and is computable for small `t, j`.
5. **Apolarity picture (J6(a)).** `TC(A,m) ⊇ Bs_R(L_m)` holds with all quantifiers; the base locus contains the cheap spans
   and more (the Fano spans of `M8`, via collinear triples of fat points). Under `Loc(b_0)` the tangent forms are
   `∩_F (pi_F^* J_F)_m^perp`, an intersection of local inverse systems, each a few-direction (possibly lattice-constrained)
   problem; Theorem 3.3 solves the disjoint-cut ones exactly.

**No growth improvement and no uniform bound follow from anything in this report.** The results are a referee
verdict on J5-J7, two strengthened finite-`n` refutations (uniform `Q_0` to `n <= 2^20`, per-subspace radii to
`n <= 2^14`), proofs of the quoted Lefschetz input and of its multi-cut extension, an exact `M8` analysis, and a sharper
formulation of the open strata question.

## Appendix — scripts and outputs

All under `refJ_checks/` (each writes `<name>_out.txt`; `m8_step*` write `m8_step*_out.txt` via `tee`).

| script | what | status |
| --- | --- | --- |
| `j5_recheck.py` | J5 as stated: exact R1 (`k <= 17`), `Q_b` (Fractions), cases A/B, direct check of the failure inequality on a `Q_0` grid, `C`-dependence, per-subspace coordinate spans | exact integers, 80-digit logs |
| `j5_extend.py` | §1.6: uniform `Q_0` with the 32-class span (`k = 12..23`; exact R1 for `k <= 20`), per-subspace with Fano/column spans | exact integers, 60-digit logs |
| `persub_fano.py` | per-subspace failure levels for the Fano span `U` and comparisons | 60-digit logs |
| `parity_box_rho.py` | `rho(P_j) = 3j`, `j <= 4`, Hilbert functions | exact (`fmpz`) |
| `xor_family_rho.py` | `rho` of XOR-closed family images, `t = 2` (`j <= 4`), `t = 3` (`j <= 2`) | mod-`p` (rigorous upper bounds on `rho`) |
| `slp_check.py` | SLP for `Q[x]/(x_i^3)`, `s <= 6`, signs; Hilbert functions of one-/multi-cut boxes vs Lefschetz products; strict unimodality `s <= 200` | exact |
| `m8lib.py` | M8 box, exact Hilbert function, width set, projection ideals, pullbacks, exact/mod-`p` dimensions | library |
| `m8_step1.py` | `|A|`, `H_A`, `dim J(A)_i` (exact), width set | exact |
| `m8_step2.py` | `dim (J_w)_i`, `i <= 8` | exact |
| `m8_step3.py` | column pairs / triples / Fano vs non-line triples | exact deg `<= 6`, sandwich deg 7-8 |
| `m8_step4.py` | explicit generators of `J(Y)` for a Fano triple | exact |
| `m8_step5.py` | planes inside Fano 3-spaces (`R = 2, 3, 5`) | exact deg 6 |
| `m8_step6.py`, `m8_step9.py` | which planes contribute; admissibility (`rho < 8`); image structure (grids, parity) | exact |
| `m8_step7.py` | all 51,933 planes spanned by pairs of low-width directions | exact deg 6 |
| `fano_baselocus.py` | Prop 2.3: `L_8` restricted to a Fano span is zero; widths of the 7 points | exact (sympy) |

Working notes: `refJ_notes.md`.
