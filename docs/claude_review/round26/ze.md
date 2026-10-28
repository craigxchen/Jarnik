# Round 26 — Route 5's zero estimate via zero/multiplicity estimates on G_m^n

Status: complete (unrefereed). Author: round-26 ZE agent. Date: 2026-10-28.

## 0. Summary and verdict

**No growth improvement and no uniform bound.** `M <= (2/3 + eps) log R/loglog R` stands. ZE-col\*, MIN^pair and
`(R_H)` at the actual point are **not** derived. **Verdict: zero/multiplicity estimates on commutative algebraic
groups cannot supply route 5's missing lemma. The reason is proved, and it is not an obstruction subtorus.**

1. **The exact inequality (proved, given the shape of Philippon's theorem, which is a memory-quote; §2).**
   * *Torus picture* (`G_m^d`). Route 5's arithmetic certifies, for one form `P_c`, vanishing at the powers `Z^j`
     only (Lemma 3.1: Liouville plus forced divisibility at `z^j`), with order `>= ord(c) - |j| Lambda(c) - O(E)`.
     The range is `|j| < ord(c)/Lambda(c) <= 1.07 sqrt(n)/lambda` (new spectral LP cap, Prop 3.3, sharpening [D]
     Thm 2.3). At regime `T = 1` boxes the certified range is `j = +-1` at `m' = ceil(4w)`, and `|j| <= 2..4` at
     `m_max` (`n <= 16384`).
   * Philippon needs vanishing on `Sigma(d)`, i.e. progressions of length about `d`. So `Sigma` is a singleton
     whenever `lambda > 2.14 sqrt(n)/d`, and (2.1) is satisfied by the trivial subgroup, with margin `e^{22784}` at
     `n = 2966` (Prop 3.5).
   * A generous variant without the sumset loss is also vacuous, for every full-rank profile in the growth regime.
2. **The obstruction subgroups (proved).**
   * Torus picture: the trivial subgroup.
   * Dual picture (`G_a^d x G_m`, one exponential polynomial vanishing on `A`; Lemma 4.1): the subgroups are
     `V x {1}`, `V x G_m` (Lemma 4.2). For nondegenerate `Z` every class count is `Z`-independent (Prop 4.3).
   * So a zero estimate decisive at `Z` must prove **ZE^all**: no nondegenerate point at all, fake positions
     included, lies on `Bs(V(A,m))`.
   * With degrees the dual estimate is vacuous (Prop 4.4). With Hilbert functions it is **false**, refuted by
     explicit points of one-column cheap subtori (Prop 4.5, exact at `n = 2966, 4096`).
   * The subtori that can legitimately obstruct are exactly the projection-cheap subtori `K_M ⊂ Bs`. These are the
     same objects that refuted ZE-braid, SL, SL\* and ZE-TC^cheap, whose refuting points (or `p`-adic shadows) lay
     **on** them.
3. **The actual point (proved).**
   * `Z` (on a full-rank subcluster) lies on **no** proper algebraic subgroup (Prop 5.1).
   * Its relative distance to a subtorus `K_U` is `>= e^{(W/2 - X_u)/2}/poly` for every defining character `u`
     (Prop 5.2). This is polynomial for light characters, recovering R12's separation lemma, and can be
     exponentially small only if every defining character is heavy.
   * This exclusion does not help: the inequality to beat is a counting inequality at the trivial subgroup, short
     by `~ d!`.
4. **A barrier for every single-form estimate (proved; exact numbers).** For **every** point `Y` there is a nonzero
   form in `V(A,mu)` carrying all vanishing data route 5 can certify, for every `mu <= m_max + 1` at
   `n = 2966..8192` (Thm 3.6). No single-polynomial zero, multiplicity or small-value estimate, with any constants,
   distinguishes the actual `Z` from fake or generic points. ZE is a statement about the whole system (R11's
   overdetermination lesson).
5. **Map consequence (proved implications, heuristic reading in §6.4).** The zero-estimate route to route 5 is
   closed: either vacuous (degree level) or reduced to ZE^all / ZE-tor^glob for overdetermined regime boxes, which
   is the barrier question of R16 §5.1. The open lemma is unchanged.

Status for the record: **obstruction_found** (a proved barrier for the zero-estimate method), with no new growth.

## 1. Setting recalled from the record (with pointers)

Labels: **proved** (all quantifiers, here or in a refereed record write-up), **exact** (integer/rational
computation), **evidence** (floating point / finite experiment, used by no proof), **memory-quote** (a literature
statement I could not re-read this session; see §2), **heuristic**. Everything new here is **unrefereed**.
Scripts: `ze_checks/` (outputs `*_out.txt` next to each script). The repository was only read.

**Cluster and its point.** `z_1..z_n` Gaussian integers, `|z_x|^2 = N`, on an arc of length `C sqrt R`
(`R = sqrt N`, `W = log N`), `d = n-1`. Relative angles `delta_x`, `|delta_x| <= Theta/2`, `Theta <= C e^{-W/4}`
([D] = `round10/determinant.md` §1). The **actual point** is

```text
Z = (z_x / z_n)_{x<n} = (e^{i(delta_x - delta_n)})_{x<n}  in  G_m^d(Q(i)),
```

which is the cluster-torus point `[z] in T_M(Q)` of round 25 (`T_n ≅ G_m^d` over `Q(i)` after choosing `z_n`
as base; round 25 intra Lemma 2.1). For `k in Z^n_0` (identified with `k' = (k_x)_{x<n} in Z^d`),
`Z^k := prod_{x<n} Z_x^{k_x} = e^{i<k,delta>}`.

**Profile.** Split primes `p`, `a_{p,x} = v_pi(z_x)`, `e_p`, `beta_p = (2a_p - e_p 1)/e_p`, `omega_p = e_p log p/W`,
`sum omega_p <= 1`. For a finite `S` in `Z^n_0` put

```text
Lambda(S) := sum_p omega_p h_S(beta_p)          (h_S(b) = max_S <k,b> - min_S <k,b>),
```

and `Lambda(c) := Lambda(supp c)`. Since `<k,a_p> = (e_p/2)<k,beta_p>` on `Z^n_0`,
`sum_p log p * h_S(a_p) = (W/2) Lambda(S)`: **`Lambda(c) W/2` is exactly the forced-divisibility cost of `c`**
([D] (I1), Thm 4.1 proof).

**Auxiliary forms.** For `c` on a finite `A ⊂ Z^n_0`: `P_c(w) = sum_k c_k w^k` (Laurent polynomial on `G_m^d`),
`ord(c) >= m` iff `c` annihilates `Poly_{<m}|_A` (equivalently `P_c` vanishes to order `>= m` at `e`; [D] §1,
R11 theory Prop 2.1). `V(A,m)`, `V_H` (integer span with `|c|_inf <= H`), `Bs(V) = {w : P_c(w) = 0 for all c in V}`.

**Slab boxes.** `T = n` ([D] Thm 2.1) or `T = 1` (`round12` t1rigour): `A = {k in Z^n_0 : |k_x| <= T (x<n),
|<k,beta_p>| <= w}`, `w = T sqrt(8 d ln 4r)`, `D = 2dT`; orders `m in [ceil(4w), m_max]`; `|A| >= (2T+1)^d/2`,
`dim V_H >= |A|/2`, `Lambda(c) <= Lambda(A) <= 2w` for every `c` on `A`.

**What is proved about `Z` (used below).**
* **Rigidity** ([D] Thm 4.1; R12 t1rigour Thm 3.1 at `T = 1`): in the regime (`d >= 316 ln 4r` at `T = 1`,
  `W > W_1`), `P_c(Z) = 0` for every `c in V_H`, i.e. `Z in Bs(V_H)`.
* **Nondegeneracy** ([D] Lemma 3.2; R11 subcluster Lemma 1.0): on a full-rank (affinely independent) subcluster,
  `Z^v != 1` for every `0 != v in Z^d`; equivalently `Z` lies in no proper algebraic subgroup of `G_m^d`, and the
  cyclic group `<Z>` is Zariski dense. Quantitatively (character inequality, [D] Prop 5.3(iii)):
  `|<v,delta>| >= 2 e^{-X_v/2}` (mod `2pi`) with `X_v = sum_p log p |<v,a_p>|`.
* **Separation** (R12 t1rigour Lemma 6.6, nearstrata §3): `|hat u_x - hat u_y| >= e^{-eps_xy/2}` and
  `dist(hat u, L) >= (1/2) e^{-N_2 abar/(2P(L))}` for collision subspaces.

**The open lemma** (R25 §8, R16 §5.1). ZE at the actual point: some `c in V_H` with `P_c(Z) != 0` (ZE-col\* for
the universal family; `(R_H)` / MIN^pair for its quantitative `T = 1` form). Conditional consequences on record:
ZE-col\* gives `M = O_C((loglog R)^2)` (R11 subcluster Thm 5.4 corrected); MIN^pair / `(R_H)` gives
`M = O(loglog R)` (R12 Thm 4.1/6.8, R16).

**Cheap subtori** (the objects that refuted ZE-braid, SL, SL\*, ZE-TC^cheap). For integer directions
`F = (v_1..v_s)` let `rho_A(F)` be the interpolation degree of `F(A) ⊂ Z^s`. If `rho_A(F) < m`, the subtorus
`K_F = {w : w^u = 1 for all u in Z^d with <u, v_i> = 0 for all i}` (the closure of `exp(i span_R F)`) lies in
`Bs(V(A,m))` ([D] Prop 5.2; R12 zetc). The width-cheap families (`sum h_A(v_i) < m`) are special cases.

## 2. The zero estimates quoted (with true hypotheses; labels: [verified-source] / [memory-quote])

**Source status.** A web search confirmed the bibliographic data: P. Philippon, *Lemmes de zéros dans les groupes
algébriques commutatifs*, Bull. SMF 114 (1986) 355–383, Théorème 2.1, with *Errata et addenda*, Bull. SMF 115
(1987) 397–398. Neither numdam.org nor arxiv.org could be reached from this session (DNS blocked by the session's
network policy), so I could **not** re-read the statement. Everything in this section is a **[memory-quote]**. The
results of §3–§5 are therefore stated so that they hold for **every** normalisation of the theorem in the range
described in Remark 2.2. They use only the *shape* of the conclusion, and they also cover a deliberately
over-generous variant.

**Theorem 2.1 (Philippon's zero estimate; [memory-quote], shape as in Philippon 1986 Thm 2.1, Masser–Wüstholz,
and Waldschmidt, *Diophantine approximation on linear algebraic groups*, 2000, Thm 5.1).**

*Hypotheses.*
* `G` is a connected commutative algebraic group of dimension `g` over `C`, with a fixed equivariant compactification
  in a product of projective spaces. The cases used here are `G = G_m^d ⊂ (P^1)^d` and
  `G = G_a^d x G_m ⊂ P^d x P^1`.
* `W ⊂ Lie G` is a linear subspace, `Sigma ⊂ G(C)` is a finite set, `T >= 0` is an integer, and `D` is a
  multidegree.
* `P` is multihomogeneous of multidegree `<= D`, not identically zero on `G`, and vanishes to order `> gT` along
  `W` at every point of `Sigma(g) = {s_1 + ... + s_g : s_i in Sigma}`. The group law is written additively, and
  WLOG `0 in Sigma`.

*Conclusion.* There is a connected algebraic subgroup `H ≠ G` with

```text
binom(T + l_H, l_H) * card((Sigma + H)/H) * H(H; D)  <=  H(G; D),        l_H = dim W - dim(W ∩ Lie H),      (2.1)
```

where `H(X; D)` is the multidegree of the closure of `X` evaluated at `D`, with the factorials included. For
example:
* `H(G_m^d; (D_1..D_d)) = d! D_1...D_d` in `(P^1)^d`;
* `H(G_a^d x G_m; (D_0, D_1)) = (d+1) D_0^d D_1` in `P^d x P^1`;
* `H(point) = 1`;
* `H(V x {1}) = D_0^{dim V}` and `H(V x G_m) = (dim V + 1) D_0^{dim V} D_1`, for a linear subspace `V`.

These degree values are exact intersection numbers (`∫ (sum D_i h_i)^{dim}` on the product of projective
spaces). They are **proved**, independently of the quotation.

**Remark 2.2 (what is used, and robustness).** The only feature of (2.1) used below is this. The conclusion is
**automatically satisfied** by `H = {e}` whenever

```text
binom(T + dim W, dim W) * |Sigma|  <=  H(G; D).                                                        (2.2)
```

The vacuity results (§3.3, §4.2) prove (2.2) with margins of `e^{1.7*10^4}` or more at the smallest regime sizes.
So they survive any misquotation of (2.1) that changes `H(G;D)` by a factor up to that size: for instance
`d!` against `d^d` normalisations, total-degree against multidegree embeddings, or `gT` against `T` in the order
hypothesis.

**Generous variant (GV).** This is the same conclusion (2.1), but under weaker hypotheses: vanishing only on
`Sigma` (not on `Sigma(g)`), and only to order `> T`. GV is **not** a theorem. It is stronger than every zero
estimate I know. It is used only to show that vacuity is not an artefact of the `Sigma(g)` and `gT` losses.

**Hilbert-sharp variant (HS).** This is (2.1) with every degree `H(X;D)` replaced by the Hilbert function
`h_X(D)`, the dimension of the space of restrictions to `X` of the polynomials of multidegree `<= D`. HS is the
form that would follow from a lossless count. §4.3 proves that **HS is false** on regime boxes.

**Related tools, not used in the proofs.**
* Multiplicity estimates of Masser–Wüstholz and Nesterenko type.
* Small value estimates (Roy, Laurent–Roy, for `G_a x G_m`).

These have the same shape: one polynomial (or a sequence of single polynomials, one per degree), a group-structured
set of points, and an obstruction subgroup. The remark in §6.3 explains why the shape, not the constants, is the
problem. This is a **[memory-quote]** at the level of shape only.

## 3. The exact zero-estimate inequality in route 5's setting

A zero estimate needs one polynomial and a *group-structured* set of points at which it is known to vanish, with
multiplicities. Route 5 supplies one point, `Z`, for a whole linear system. The first task is to find every point
and every multiplicity that route 5's arithmetic can actually certify for a single `P_c`. The answer is the powers
`Z^j`, up to `|j| < ord(c)/Lambda(c)`, with orders decreasing linearly (§3.1). The second task is to bound this
range: it is at most `1.07 sqrt(n)/lambda` (§3.2). The third is to substitute the result into (2.1) (§3.3–3.4).

### 3.1 Extrapolation to powers (what route 5 certifies for one form)

**Lemma 3.1 (forced vanishing at `Z^j`; proved).** Let `z_1..z_n` be a cluster (§1; any `n`, any norm), `A ⊂ Z^n_0`
finite, `D >= max_A |k|_1`, `T_A = max_A |k'|_inf`, `C' = max(C,1)`. Let `c in Z^A \ 0`, `m <= ord(c)`,
`j in Z \ 0`, `t >= 1`. Suppose that for every `s in {0, ..., t-1}`

```text
(m - s - |j| Lambda(c)) W/2  >  2 ln|c|_1 + 2 s ln T_A + 2 (m-s) ln(D |j| C'/2) - 2 ln((m-s)!).            (3.1)
```

Then `P_c` vanishes to order `>= t` at `Z^j`.

*Proof.*

*Step 1: one point, order 0.* Let `c'` be a nonzero integer vector on `A`, with `ord(c') >= mu'`. Suppose

```text
(mu' - |j| Lambda(c')) W/2  >  2 ln|c'|_1 + 2 mu' ln(D|j|C'/2) - 2 ln mu'! .
```

We show `P_{c'}(Z^j) = 0`.
* First take `j > 0`. Put `y_x = z_x^j`. These are Gaussian integers of norm `N^j`, with `v_pi(y_x) = j a_{p,x}`,
  and `p^{j e_p} || N^j`. Put
  `G_k(y) = N^{j(D-|k|_1)/2} prod_x y_x^{k_x^+} conj(y_x)^{k_x^-} in Z[i]` (note `|k|_1` is even on `Z^n_0`).
  Then `G_k(y) = N^{jD/2} Z^{jk}`.
* [D] (I1) applied to the tuple `y` gives `v_pi(G_k(y)) = j(e_p D/2 + <k,a_p>)` and
  `v_{conj pi}(G_k(y)) = j(e_p D/2 - <k,a_p>)`. Inert primes and `1+i` divide every `G_k(y)` to the full power.
  The proof is [D] Thm 4.1's, verbatim, with `z` replaced by `y`. The exact check is `ze_checks/powers_valuation.py`:
  1620 checks for `j = 1, 2, 3`, 0 failures.
* So if `F := sum_k c'_k G_k(y) != 0`, then
  `Norm F >= N^{jD} prod_p p^{-j h_{supp c'}(a_p)} = N^{jD} e^{-jW Lambda(c')/2}`.
* For the upper bound, the moments of `c'` of order `< mu'` vanish, and `|<k, j delta>| <= D|j|Theta/2` with
  `Theta <= C' e^{-W/4}`. Hence `|F| <= N^{jD/2} |c'|_1 (D j Theta/2)^{mu'}/mu'!`.
* Since `Norm F = |F|^2`, the two bounds contradict the displayed hypothesis. So `F = 0`, i.e. `P_{c'}(Z^j) = 0`.
* For `j < 0`, `P_{c'}(Z^j)` is the complex conjugate of `P_{c'}(Z^{-j})`, because `|Z_x| = 1` and `c'` is real.

*Step 2: derivatives.* For a multi-index `alpha` on the `d` coordinates `x < n`, the invariant derivations
`theta_x = w_x d/dw_x` satisfy `theta^alpha P_c = P_{c^(alpha)}`, where `c^(alpha)_k = c_k k'^alpha`. Here:
* `ord(c^(alpha)) >= m - |alpha|`, because `k'^alpha P(k)` has degree `< m` when `deg P < m - |alpha|`;
* `|c^(alpha)|_1 <= |c|_1 T_A^{|alpha|}`;
* `supp c^(alpha) ⊂ supp c`, so `Lambda(c^(alpha)) <= Lambda(c)`.

Step 1 applies to `c^(alpha)` (or `c^(alpha) = 0`) with `mu' = m - |alpha|`. The `theta_x` form a basis of
invariant vector fields, so the vanishing of all `theta^alpha P_c` with `|alpha| < t` at a point of `G_m^d` is
equivalent to vanishing to order `>= t` there. QED

**Corollary 3.2 (route-5 data, regime numbers; exact `m_max`, floating budget).** On a `T = 1` slab box every `c`
satisfies `Lambda(c) <= 2w`, `ln|c|_1 <= 3d ln 3` (for `c in V_H`, `H = U^2`), `T_A = 1` and `D = 2d`. So Lemma 3.1
certifies, for every `c in V_H(A,m)`, vanishing at `Z^j` for `|j| <= J_lemma(m)`, where `J_lemma(m)` is the largest
`j` with `m - 2jw > E(m,j)` and `E` is the right side of (3.1) at `s = 0`, divided by `W/2`.

At `W = W_min = 3n ln n - 14n` and `C = 1` (`ze_checks/vacuity_numbers_out.txt`):

| `n` | `m_max` | `w` | `ceil(4w)` | `E(m_max,1)` | `J_lemma(m_max)` | `J_lemma(ceil 4w)` |
| --- | --- | --- | --- | --- | --- | --- |
| 2966 | 1906 | 471.7 | 1887 | 1.69 | 2 | 1 |
| 4096 | 2631 | 563.8 | 2256 | 1.54 | 2 | 1 |
| 8192 | 5260 | 825.4 | 3302 | 1.30 | 3 | 1 |
| 16384 | 10518 | 1205.6 | 4823 | 1.12 | 4 | 1 |
| `10^6` | `>= 641697` (proved floor) | 11027.9 | 44112 | 0.62 | 29 | 1 |

So at the R11/R12 recommended order `ceil(4w)`, route 5 certifies **only `Z` itself** (`j = +-1`, the latter by
conjugation). At `m_max` it certifies a symmetric progression `Z^{-J}..Z^{J}` with `J = 2..4` in the regime range.
The round-11 remark (`theory.md` l. 465: "`j` up to about `m/(2w) ~ sqrt(n/ln r)`, far too few") is confirmed and
made quantitative. The next proposition shows that **no** choice of support can do much better on well-conditioned
profiles.

### 3.2 The spectral cap on the extrapolation range

Let `P = I - 11^T/n` be the projection onto `1^perp ⊂ R^n`, let `bt_p = P beta_p`, and put
`G := sum_p omega_p bt_p bt_p^T`, a PSD operator on `1^perp`. Let `lambda := lambda_min(G|_{1^perp})`.

* For sign columns (`e_p = 1`) with `sum omega_p = 1`, `G` is the projected Gram matrix of R10 Lemma 3.1, with
  entries `1 - 2 d_xy/W'` before projection.
* For the Sylvester–Hadamard profile, `lambda = n/(n-1)` (exact eigenvalue; `ze_checks/spectral_cap_out.txt`
  confirms it at `n = 4, 8, 16`).

**Proposition 3.3 (spectral LP cap; proved).** Suppose `lambda > 0`. Then for every finite `S` in a slice of
`Z^n` and every nonzero complex `c` on `S`,

```text
ord(c)  <=  Omega(S)  <=  sum_p omega_p h_S(beta_p) |G^+ bt_p|_2  <=  (sqrt(n)/lambda) Lambda(S).               (3.2)
```

*Proof.*
1. The first inequality is (I2) (`threshold/construct.md` Thm B, Cor B1; it holds over any field of
   characteristic 0).
2. On a slice, `<k,1>` is constant, so `h_S(v) = h_S(Pv)`.
3. Because `lambda > 0`, the vectors `bt_p` span `1^perp`. So `P eps = G G^+ P eps = sum_p omega_p bt_p <G^+ bt_p, eps>`.
4. `h_S` is a seminorm (subadditive, `h_S(a v) = |a| h_S(v)`). So
   `h_S(eps) <= sum_p omega_p |<G^+ bt_p, eps>| h_S(beta_p)`.
5. Average over `eps` uniform in `{+-1}^n`: `E|<u,eps>| <= (E<u,eps>^2)^{1/2} = |u|_2`. This gives the middle
   inequality.
6. Finally, `|G^+ bt_p|_2 <= |bt_p|_2/lambda <= |beta_p|_2/lambda <= sqrt(n)/lambda`. QED

For Hadamard profiles (3.2) gives `ord(c) <= sqrt(n) ((n-1)/n) Lambda`, which slightly sharpens [D] Thm 2.3.
`ze_checks/spectral_cap.py` checks (3.2) on 49 random profiles (`n <= 7`, random weights, random supports): `rho(S)`
is computed exactly, `Omega(S)` exactly, and the spectral quantities in floating point. There are 0 failures, and
the maximum of `rho/cap` is 0.19 (evidence).

**Corollary 3.4 (extrapolation range; proved).** Assume `W >= 40` and `D >= ord(c)`; the latter holds on slab boxes,
where `ord(c) <= reg(A) <= 2dT = D`. Then (3.1) with `t = 1` can hold only if

```text
|j|  <  ord(c) (1 + 4 ln 2/W) / Lambda(c)   <=   1.07 sqrt(n)/lambda.                                     (3.3)
```

Moreover, if the exponent map `k -> (<k,a_p>)_p` is injective on `Z^n_0` (full exponent rank, as on the
nondegenerate subclusters of [D] Lemma 3.2), then `|j| < 1.07 dT W/log 5` for **every** profile.

*Proof.* At `s = 0` the right side of (3.1) is bounded below: `ln|c|_1 >= 0`, `|j|, C' >= 1`, and
`2m ln(D/2) - 2 ln m! >= -2m ln(2m/D) >= -2m ln 2`, because `m! <= m^m` and `m <= D`. So (3.1) implies
`m - |j|Lambda(c) > -4m ln 2/W`. Apply this with `m = ord(c)`; then (3.2) gives the second bound in (3.3).

For the last claim:
* `ord(c) >= 1` forces `|supp c| >= 2`.
* Injectivity gives some `p` with `h_{supp c}(a_p) >= 1`. So `Lambda(c) >= (2 log p/W) h_{supp c}(a_p) >= 2 log 5/W`,
  since split odd primes are `>= 5`.
* Also `ord(c) <= reg(supp c) <= 2dT`, by Lagrange interpolation in each coordinate. QED

### 3.3 The inequality (2.1) in route 5's setting (torus picture)

**Data.**
* `G = G_m^d ⊂ (P^1)^d` and `W_Lie = Lie G` (dimension `d`).
* `P = w^{T 1} P_c`, a polynomial of multidegree `<= (2T, ..., 2T)`. So `H(G;D) = d!(2T)^d`.
* The certified vanishing set is `{Z^j : |j| <= J}` (Lemma 3.1), with orders `tau_j <= ord(c)(1.07) - |j|Lambda(c)`.
* Translating by `Z^{j_0}` preserves multidegrees and orders on a torus.

The exact zero-estimate inequality for route 5 is (2.1) with these data:

```text
binom(T_0 + l_H, l_H) * card((Sigma H)/H) * H(H; 2T,...,2T)  <=  d! (2T)^d,         Sigma(d) ⊂ {Z^j : |j| <= J},   d T_0 < min_{Sigma(d)} tau.   (3.4)
```

**Proposition 3.5 (the torus-picture zero estimate is vacuous; proved, given the shape of (2.1)).**
1. *(singleton.)* If `J < d/2`, every admissible `Sigma` is a single point up to translation. By Cor 3.4 this
   holds for every profile with `lambda > 2.14 sqrt(n)/d`, a mild well-conditioning hypothesis (`lambda > 0.04` at
   `n = 2966`). Sylvester profiles have `lambda = n/(n-1)` and satisfy it for `n >= 8`. For profiles violating it,
   item 2 applies; it also covers the real theorem, whose hypotheses are stronger than GV's.
   * Then (2.2) holds: `binom(T_0 + d, d) <= binom(2T + d, d) <= d!(2T)^d` for `d >= 3`, because
     `dT_0 < tau <= ord(c) <= 2dT`.
   * So (2.1) is satisfied by `H = {e}`, and the zero estimate yields **no information**.
2. *(generous variant GV; any profile of full exponent rank.)* Under GV, `|Sigma| <= 2J+1` with
   `J < 1.07 dT W/log 5` and `T_0 < 2dT`. Then (2.2) holds as soon as

   ```text
   ln(2.14 dT W/log 5 + 1)  <=  S_d := ln(d!(2T)^d) - ln binom((2T+1)d, d),      where S_d >= d ln(d/(1.5 e^2))  (d >= 12).
   ```

   At `T = 1`, `n = 2966` the right side is `17142` (`vacuity_numbers_out.txt`), so GV is vacuous for every cluster
   with `log_10 W < 7441`. This includes the entire growth regime `W = O(n log n)`, and in general every `W` with
   `ln W <= d ln(d/(1.5e^2)) - ln(2.2 dT/log 5)`.

*Proof.*
1. `Sigma(d)` contains `d s` for each `s in Sigma` (with `0 in Sigma`). So a non-singleton `Sigma` inside a
   progression of `2J+1` terms requires `d <= 2J`.
2. The binomial inequalities are `binom(2T+d,d) <= (2T+1)^d <= d!(2T)^d` (since `d! >= 1.5^d` for `d >= 3`) and
   `binom((2T+1)d, d) <= (e(2T+1))^d`, together with `d! >= (d/e)^d`. QED

**Reading.** Philippon's theorem needs vanishing along a progression of length about `g = d` inside `<Z>`, with
multiplicity `g T_0`. Route 5 certifies a progression of half-length at most `1.07 sqrt(n)/lambda` (`~ 58` for `lambda = n/(n-1)` at
`n = 2966`; the range guaranteed for every `c in V_H` on regime slab boxes is only `1..4`, Cor 3.2), and
multiplicity at most `ord(c) <= 2dT`. The `Sigma(g)` requirement alone kills the estimate
(item 1). Even without it, the `d!` in the degree swamps every progression route 5 can certify (item 2).

### 3.4 A barrier for every single-polynomial zero estimate (Hilbert level)

Prop 3.5 uses the degree normalisation. The next theorem uses only linear algebra. It shows that the vacuity is
not a weakness of Philippon's constants but a property of the data.

**Theorem 3.6 (single-polynomial barrier; proved; exact numbers).** Let `A ⊂ Z^n_0` be finite, `mu >= 1`,
`J >= 1`, and let `tau_j >= 0` be integers for `1 <= |j| <= J`. If

```text
dim V(A, mu)  >  sum_{1<=|j|<=J} binom(tau_j + d - 1, d),                                                  (3.5)
```

then for **every** `Y in G_m^d(C)` there is a nonzero `c in V(A,mu) ⊗ C` such that `P_c` vanishes to order
`>= mu` at `e` and to order `>= tau_j` at `Y^j` for all `1 <= |j| <= J`.

*Proof.* Vanishing to order `>= tau` at a point imposes `binom(tau + d - 1, d)` linear conditions on `c` (one for
each derivative `theta^alpha`, `|alpha| < tau`). The order at `e` is built into `V(A,mu)`. QED

**Application (proved for the theorem; exact for the numbers).** Take the profile with `lambda = n/(n-1)`
(Sylvester) and the `T = 1` box. Put `J = floor(1.07 sqrt(n)/lambda)` and
`tau_j = floor(mu(1 + 4ln2/W - |j|lambda/sqrt n)) + 1`. By Lemma 3.1, Prop 3.3 and Cor 3.4, these orders dominate
everything route 5 certifies for a form `c` with `ord(c) = mu`. Use
`dim V(A,mu) >= |A| - binom(mu+d-1, d) >= ceil(3^d/2) - binom(mu+d-1,d)`. Then (3.5) holds, in exact big-integer
arithmetic (`ze_checks/single_poly_barrier_out.txt`), at

| `n` | `mu = ceil(4w)`: `ln` count / `ln` avail | `mu = m_max`: `ln` count / `ln` avail | largest `mu` for which (3.5) holds |
| --- | --- | --- | --- |
| 2966 | 3205.4 / 3256.7 | 3223.5 / 3256.5 | 1907 (`m_max + 1`) |
| 3000 | 3233.3 / 3294.0 | 3260.1 / 3293.8 | 1929 |
| 4096 | 4090.8 / 4498.1 | 4458.9 / 4498.0 | 2633 |
| 8192 | 6842.2 / 8998.0 | 8942.7 / 8997.9 | 5262 |

So, for every order up to the Siegel range and at **every** point `Y` (fake positions, generic points, the actual
`Z`), some nonzero form in the system satisfies all the vanishing that route 5 can certify for that form.

**Consequence.** Any statement of the form "a nonzero `P` with support `A` and these vanishing data forces `Y` into
an exceptional set `E`" has `E = G_m^d`. Such statements include zero estimates, multiplicity estimates and small
value estimates at a single degree, with any constants. For Philippon-type statements the conclusion is a
`Y`-independent inequality for nondegenerate `Y` (§4.1). So it holds automatically, as Prop 3.5 found. Orders
`mu > m_max + O(1)` (forms whose existence is structural, not forced by counting) are not covered by Thm 3.6. They
are covered at the degree level by Prop 3.5.

**What this says about ZE.** ZE is a statement about the **whole system** `V_H` (dimension `>= |A|/2`) at one point.
Restricting attention to one form loses exactly the overdetermination, which is the only feature that separates ZE
from the false statements of R11 (Prop 6.2: `dim V = 1` has nondegenerate zeros; thin families). Single-polynomial
zero estimates are blind to it by construction.

## 4. Which subtori obstruct; comparison with the cheap rational subtori that refuted SL/SL*

There are two ways to put route 5's data into a group. In the **torus picture** (§3.3) there is one form `P_c` on
`G_m^d` and points in `<Z>`. In the **dual picture** (§4.2) the whole system is replaced by one exponential
polynomial on `G_a^d x G_m`, which vanishes on the large set `A`. The dual picture is the only one in which a zero
estimate sees the overdetermination. It is also the one in which the obstruction subgroups can be listed
completely.

### 4.1 Torus picture: the obstruction is the trivial subgroup

The connected algebraic subgroups of `G_m^d` are the subtori `K_U = {w : w^u = 1 for all u in U}`, one for each
saturated sublattice `U ⊂ Z^d`. Let `Sigma ⊂ <Z>` and `Z` be nondegenerate. Then `Z^i in K_U` forces `Z^{iu} = 1`
for `u in U`. So `card(Sigma K_U/K_U) = |Sigma|` for every proper `K_U`, and (3.4), with `W_Lie = Lie G` so that `l_H = codim H`, reads

```text
binom(T_0 + codim K_U, codim K_U) * |Sigma| * H(K_U; 2T,...,2T)  <=  d!(2T)^d ,      H(K_U; .) >= 1 (a degree),
```

By Prop 3.5 this already holds for `K_U = {e}` (`codim = d`, `deg = 1`). So the "obstruction subgroup" that the
theorem returns for route 5's data is the **trivial subgroup**. That is a counting deficit, not a subtorus.

**Comparison with the cheap subtori.** For a cheap family `F` (`rho_A(F) < m`), `K_F ⊂ Bs(V(A,m))`. For a point
`Y in K_F` and `Sigma ⊂ <Y>`, `card(Sigma K_F/K_F) = 1`, and `P_c` vanishes on all of `K_F`. So `K_F` is the
natural obstruction subgroup for points **on** it. These are exactly the subtori that refuted:
* ZE-braid (Sylvester `log_2 n` coordinate columns; R11 referee);
* SL and SL\* (twins, duplicates, and the column identity of R11 J1: the `p`-adic point `p^{a'_p}` lies on a cheap
  `K_F`);
* ZE-TC^cheap (R12 t1rigour Thm 6.11, zetc C1).

In each case the refuting object was a point, or a `p`-adic shadow, lying **on** a cheap `K_F`. The archimedean
actual point `Z` lies on **no** proper subtorus (nondegeneracy; §5). So the cheap subtori cannot be what blocks a
zero estimate at `Z`. What blocks it is that the data at `Z` are too thin (Prop 3.5, Thm 3.6).

### 4.2 Dual picture: one exponential polynomial, many points

**Lemma 4.1 (proved).** For `Y in G_m^d(C)`: `Y in Bs(V(A,m))` iff there is a polynomial `p in C[X_1..X_d]` of
degree `< m` with `Y^k = p(k')` for every `k in A`. For `V_H ⊂ V(A,m)` only "`⇐` implies `Y in Bs(V_H)`" holds.

*Proof.* `V(A,m)` is the annihilator of `Poly_{<m}|_A` under `<c, f> = sum_A c_k f(k)`. In finite dimension
`(Ann E)^perp = E`, so `(Y^k)_{k in A} ⊥ V(A,m)` iff `(Y^k)_A in Poly_{<m}|_A`. QED

Put `phi_Y : Z^d -> G' := G_a^d x G_m`, `k -> (k, Y^k)`, a group homomorphism, and `Q = X_{d+1} - p(X)`, of
bidegree `(m-1, 1)`. Then `Y in Bs(V(A,m))` iff `Q` vanishes on `phi_Y(A)`. This is the classical setting of zero
estimates for exponential polynomials, with no multiplicities: only values at lattice points are known.

**Lemma 4.2 (proved; standard).** The connected algebraic subgroups of `G_a^d x G_m` over `C` are exactly
`V x {1}` and `V x G_m`, for `V ⊂ C^d` a linear subspace.

*Proof.*
* A commutative linear algebraic group is the product of its unipotent and semisimple parts, and both are closed.
* The unipotent part of a closed connected `H` lies in `G_a^d x {1}`. In characteristic 0 a connected closed
  subgroup of `G_a^d` is a linear subspace.
* The semisimple part is a connected subgroup of `{0} x G_m`, i.e. trivial or `G_m`.
* `Hom(G_a, G_m) = Hom(G_m, G_a) = 0`, so no graph-type subgroups occur. QED

**Proposition 4.3 (`Z`-blindness; proved).** Let `Y` be nondegenerate (`Y^v != 1` for `0 != v in Z^d`), and let
`Sigma ⊂ Z^d` be finite. Then

```text
card(phi_Y(Sigma) + (V x G_m) / (V x G_m)) = |pi_V(Sigma)|,        card(phi_Y(Sigma) + (V x {1}) / (V x {1})) = |Sigma| ,
```

independently of `Y`. Consequently, for any normalisation of (2.1) whose hypothesis is "`Q` of bidegree
`(m-1,1)` vanishes on `phi_Y(Sigma(g))`", the following holds. If it excludes **one** nondegenerate `Y` from
`Bs(V(A,m))`, it excludes **all** nondegenerate `Y` with the same `(A, m)` and `Sigma`. So a zero estimate in the
dual picture can prove ZE at the actual point only by proving

> **ZE^all(A,m):** `Bs(V(A,m))` contains no nondegenerate point of `G_m^d`.

This is a `delta`-free, height-free, global version of round 10's ZE(n). In particular it uses **nothing** about
`Z` except nondegeneracy: not that `Z in G_m^d(Q(i))`, not its height `~ W/2`, not `|delta| <= C e^{-W/4}`, and not
the pair inequalities.

*Proof.*
* For `V x G_m`: `(k, Y^k) ~ (k', Y^{k'})` iff `k - k' in V`.
* For `V x {1}`: `(k, Y^k) ~ (k', Y^{k'})` iff `k - k' in V` and `Y^{k-k'} = 1`. The second condition forces
  `k = k'`, by nondegeneracy.
* Neither count involves `Y`. So the conclusion of (2.1) is the same inequality for every nondegenerate `Y`.
* The hypothesis "`Q` vanishes on `phi_Y(Sigma(g))`" holds for every `Y in Bs(V(A,m))`, with `Sigma(g) ⊂ A`.
* If (2.1) fails for all `H` (a contradiction for `Z`), it fails for all `H` at every nondegenerate `Y`. QED

**Proposition 4.4 (the dual-picture zero estimate is vacuous; proved given the shape of (2.1)).** Let `A` be a slab
box with `m >= 4w` (any `T`). Then (2.2) holds for every `Sigma ⊂ A`, even under GV. Indeed
`|Sigma| <= |A| <= (2T+1)^d <= (d+1)(m-1)^d = H(G'; m-1, 1)`, because `m - 1 >= 4w - 1 >= 2T + 1`
(`w = T sqrt(8 d ln 4r) >= T`). The margins at `T = 1` are `ln 3^d = 3257` against `ln(d+1)(m-1)^d >= 22370` at
`n = 2966` (`vacuity_numbers_out.txt`).

*Remark (`T = n` and `T = 1` are `Sigma` and `Sigma(g)`).* With the real hypothesis `Sigma(d+1) ⊂ A`, the `T = 1`
slab box `A_1` is an admissible `Sigma` for the `T = n` box `A_n`. Its `(d+1)`-fold sums satisfy
`|k_x| <= d+1 = T` and `|<k,beta_p>| <= (d+1) w_1 = w_n`. So the strongest honest use of Philippon at `T = n` takes
`Sigma = A_1`. It is still vacuous: `3^d << (d+1)(m_n - 1)^d`, with `m_n ~ 0.7 n^2`. At `T = 1` no nontrivial
`Sigma` exists, because `Sigma(d+1) ⊂ {-1,0,1}^d` forces `Sigma = {0}`.

**Proposition 4.5 (the Hilbert-sharp variant is false; proved, with exact big-integer verification of the final
inequalities).** Let `A` be a `T = 1` slab box with `w >= 10` and `m in [ceil(4w), m_max]`. Let `p` be a column with
`L := h_A(a_p) in [3, m-1]`. For a sign column with at least two 0's and two 1's, `L >= 4`, and `L <= w < m` by the
slab. Let `theta in R` with `theta/pi` irrational, and put

```text
Y = (e^{i theta (a_{p,x} - a_{p,n})})_{x<n}   (a point of the cheap subtorus K_{(a_p)}, |Y_x| = 1).
```

Then `Q = X_{d+1} - g(<k, a_p>)` vanishes on `phi_Y(A)`, where `g` is the Lagrange interpolant of `y -> e^{i theta y}`
on the `L + 1` values of `<k,a_p>`, so `deg g <= L < m`. Yet **no** connected `H ≠ G'` satisfies
`card(phi_Y(A)H/H) h_H(m-1,1) <= h_{G'}(m-1,1)`. So HS (= GV with Hilbert functions) is false.

*Proof.* The Hilbert functions are `h_{G'} = 2 binom(m-1+d, d)`, `h_{V x G_m} = 2 binom(m-1+e, e)` and
`h_{V x 1} = binom(m-1+e, e)`, with `e = dim V`.

Classes.
* For `V x G_m` there are `|pi_V(A)|` classes.
* For `V x {1}` there are `|pi_U(A)|` classes, where `U = V ∩ a_p^perp`. Here `Y^u = e^{i theta <u,a_p>} = 1` iff
  `<u, a_p> = 0`. This count is also `>= L + 1`, the number of distinct `Y`-values.

Two lower bounds for `|pi_U(A)|`, where `U` has codimension `s` and dimension `e_U`:
* (i) `|pi_U(A)| >= 3^{min(s, floor(w/2))}`. Choose `s' = min(s, floor(w/2))` coordinates `x < n` whose `e_x` are
  independent modulo `U`. The `3^{s'}` vectors `k' in {-1,0,1}` supported there (with `k_n = -sum k'`) lie in `A`,
  because `|<k,beta_q>| <= 2s' <= w`, and they have distinct images.
* (ii) `|pi_U(A)| >= |A|/3^{e_U} >= 3^d/(2 * 3^{e_U})`. An `e_U`-dimensional affine subspace meets
  `{-1,0,1}^d` in at most `3^{e_U}` points, since it projects injectively to `e_U` coordinates.

`ze_checks/hs_refute.py` takes the worst case `L = 3` and verifies, in exact integers for every `e = 0..d`, that
`max(L+1, LB(e)) binom(m-1+e,e) > h_{G'}` and `LB(e) * 2binom(m-1+e,e) > h_{G'}` (`e < d`). Here
`LB(e) = max(3^{min(d-e, floor(w/2))}, ceil(ceil(3^d/2)/3^e))`, and `e_U <= e`, `s >= d - e`. Results at
`n = 2966, 4096` and both ends of the order range: all hold. The tightest case is the hyperplanes `e = d - 1`, with
`ln(LHS/RHS) >= 0.197`.

The general argument for all `n` is as follows.
* At `e = d - 1` the ratio is `4d/(2(m-1+d)) > 1`, since `m - 1 < d`.
* For `d - e <= w/2`, use `3^{d-e} > 2 prod_{i=e+1}^{d} (m-1+i)/i`, which holds once `d - e >= 2`.
* For `d - e > w/2`, `f(e) = binom(m-1+e,e)/3^e` is unimodal with `f(0) = 1`, and
  `f(d - floor(w/2)) >= f(d) * 1.6^{floor(w/2)}`. Both exceed `4 f(d)`.

I state the theorem only for the verified `(n, m)` pairs, and record the all-`n` argument as a sketch. QED

**Reading of §4.2.**
* With the true (degree) normalisation, the dual zero estimate is vacuous: the trivial subgroup obstructs (Prop
  4.4).
* With the lossless (Hilbert) normalisation it is false, because points of cheap subtori violate it (Prop 4.5).
* Any intermediate statement that is decisive at `Z` is decisive at every nondegenerate point (Prop 4.3), i.e. it is
  ZE^all.

**Which subtori obstruct, exactly.** For points `Y` of the torus, the subtori that can legitimately appear as
obstructions are the `K_M` with `Y in K_M ⊂ Bs(V(A,m))`. By Lemma 4.1, `K_M ⊂ Bs` iff every function of `pi_M(k)`
agrees on `A` with a polynomial in `k` of degree `< m` (exponentials with distinct frequencies span all functions on
`pi_M(A)`). Lagrange interpolation shows that this holds when `rho(pi_M(A)) < m`. These are the projection-cheap
subtori of R12 (`rho_A`, tangent degree `nu_A`). These are precisely the
subtori that refuted ZE-braid, SL, SL\* and ZE-TC^cheap. Through `Y = Z` none of them passes. A **sharp** zero
estimate would be

> **ZE-tor^glob(A,m):** `Bs(V(A,m))` is the union of the subtori `K_M ⊂ Bs(V(A,m))` described above, up to torsion
> translates.

This is the global form of R11's ZE-tor (referee_subcluster §7.1 made it local). It implies ZE^all, hence ZE.
Neither zero-estimate theory nor anything in the record proves it. R12 fake Thm 2.2 (exact fakes on `Z_P` in the
underdetermined `sdc7`, `sdc11` models) shows that the analogue fails for underdetermined boxes. So ZE-tor^glob is
again an overdetermination statement, i.e. the barrier question of R16 §5.1.

## 5. Can the actual point Z lie on an obstruction subtorus?

**Proposition 5.1 (exact exclusion; proved).** Let `J` be a full-rank subcluster (affinely independent exponent
vectors, [D] Lemma 3.2; on squarefree norms one of size `>= n/4 - 1` exists, [D] Lemma 3.1). Then its point `Z_J`
lies in **no** proper algebraic subgroup of `G_m^{|J|-1}`, including torsion translates `zeta K_U` of proper subtori.
The same holds for every power `Z_J^i`, `i != 0`. Hence, for every proper `H` and every finite `Sigma ⊂ <Z_J>`,
`card(Sigma H/H) = |Sigma|`. In the dual picture, every `V x {1}` sees `|Sigma|` classes.

*Proof.* Algebraic subgroups of a torus are finite unions of torsion translates of subtori. Suppose
`Z^i in zeta K_U` with `u in U \ 0`. Then `Z^{iu} = zeta^u`, a root of unity of some order `N_0`, so
`Z^{iN_0 u} = 1`. This contradicts nondegeneracy ([D] Lemma 3.2: `sum v_x delta_x != 0` for `0 != v in Z^J_0`;
`Z^v = e^{i<v,delta>}` with `|<v,delta>| < 2 pi`). QED

**Proposition 5.2 (quantitative position; proved, from the character inequality).** Let `K_U` be a proper subtorus
and `L_U = U^perp ⊂ R^d` its real Lie algebra. Then

```text
dist_2(delta', L_U) / |delta'|_inf  >=  max_{u in U}  (4/(C' |u|_2)) * exp((W/2 - X_u)/2),     X_u = sum_p log p |<u,a_p>|,
```

where `delta' = (delta_x - delta_n)_{x<n}`.

*Proof.* [D] Prop 5.3(iii) gives `|<u,delta'>| >= 2 e^{-X_u/2}` for every `u` with `<u,delta'> != 0`. This
applies to all `u != 0`, by nondegeneracy. Moreover `dist_2(delta', L_U) >= |<u,delta'>|/|u|_2` and
`|delta'|_inf <= Theta <= C' e^{-W/4}`. QED

**Reading.**
* For `u` a light character (`X_u = W/2 + O(log)`), the relative distance is polynomial. The pair characters
  `u = e_x - e_y` give exactly R12's per-pair gaps `e^{-eps_xy/2}` and the separation lemma
  `dist(hat u, L) >= (1/2) e^{-N_2 abar/(2P(L))}` for collision subspaces.
* At `T = 1`, `m' = ceil(4w)`, every width-cheap span is a collision subspace with `<= 16` classes (R12 t1rigour Cor
  6.4). So `Z` is polynomially far from all width-cheap subtori (proved, R12).
* For a projection-cheap subtorus whose defining characters are **all heavy** (`X_u >> W/2`), Prop 5.2 gives only
  an exponentially small lower bound. Whether such subtori occur in regime boxes is R12's open point (Q1) /
  ZE-TC^coll. MIN\* fails at `n = 4096` through a 3-column span with collisions (R12 referee_nearstrata), so closeness
  at finite `n` is real.

**Decision on item (3) of the task.**
1. **`Z` lies on no obstruction subtorus**, exactly (Prop 5.1). This is what the record's nondegeneracy gives. Round
   12's rigidity (t1rigour Thm 3.1) places `Z` **in** `Bs(V_H)`; it does not place it on a subtorus.
2. **This does not let the zero estimate close.** For route 5's data, the zero estimate's conclusion (2.1) is
   already satisfied by the **trivial** subgroup (Prop 3.5 torus picture, Prop 4.4 dual picture). So excluding `Z`
   from every positive-dimensional subtorus changes nothing. The inequality to beat is a counting inequality, and
   the counts are short by factors of order `d!`:
   * `e^{22784}` at `n = 2966` in the torus picture (`ln d! 2^d - ln binom(T_0+d,d)`);
   * `e^{19113}` in the dual picture (`ln(d+1)(m-1)^d - ln 3^d`).
3. **Hence neither ZE-col\* nor MIN^pair is derived.** I do not claim either, nor any growth.
4. **The explicit obstruction**, in the three senses in which one exists:
   * (a) *The trivial subgroup.* It is the obstruction returned by (2.1) for every profile, actual or fake, with the
     margins above.
   * (b) *A single-form witness at every point.* For every `Y in G_m^d` (the actual `Z` of any actual cluster, a fake
     position, a generic point), Thm 3.6 gives a nonzero form in `V(A,mu)` that carries **all** vanishing data
     route 5 certifies for forms of order `mu`, for every `mu <= m_max + 1` (exact, `n = 2966..8192`). This witness
     depends only on the profile, which is the data actual and fake clusters share.
   * (c) *An explicit counterexample to the lossless dual zero estimate.* The points
     `Y = e^{i theta (a_p - a_{p,n})}` on the cheap subtorus of one column (Prop 4.5). The construction uses only
     one column of the profile, so it applies to fake and actual profiles alike.
5. **Track of dependence on the profile.**
   * Every threshold in §3–§4 depends on the profile only through:
     * `lambda` (Prop 3.3, Cor 3.4), which is `n/(n-1)` for Sylvester profiles;
     * `log 5` (Cor 3.4, full rank);
     * `w`, `m`, `d`, which are uniform functions of `(n, r)`.
   * No template-dependent threshold is treated as uniform. Where a statement needs `lambda >= lambda_0` it says so
     (Prop 3.5(1)). The version valid for **all** full-rank profiles (Prop 3.5(2)) uses only `log 5`, and is uniform
     in the profile for `ln W <= d ln(d/(1.5e^2)) - ln(2.2dT/log 5)`.

## 6. What closes and what does not (ZE-col*, MIN^pair)

### 6.1 Logical summary

| step | status | where |
| --- | --- | --- |
| route 5 certifies `P_c(Z^j) = 0` with order `>= ord(c) - |j|Lambda(c) - O(E)` for `|j| < ord(c)/Lambda(c)` | **proved** (new as a lemma; anticipated as a remark in R11 theory l. 465) | Lemma 3.1 |
| extrapolation range `<= 1.07 sqrt(n)/lambda` for **every** support | **proved** | Prop 3.3, Cor 3.4 |
| at regime `T = 1` boxes: `J = 1` at `ceil(4w)`, `J = 2..4` at `m_max` (`n <= 16384`) | **exact** `m`; floating budget | Cor 3.2 |
| Philippon on `G_m^d` with route-5 data: `Sigma` a singleton (if `lambda > 2.14 sqrt(n)/d`); (2.1) holds at `H = {e}`; no information | **proved** given the shape of (2.1); margin `e^{2.2*10^4}` at `n = 2966` | Prop 3.5(1) |
| GV (no sumset, low order) on `G_m^d`: vacuous for every full-rank profile and every `W` with `ln W <~ d ln(d/11)` | **proved** | Prop 3.5(2) |
| no single-form statement with any constants separates `Z` from any other point | **proved** (Thm); **exact** numbers `n = 2966..8192`, all `mu <= m_max + 1` | Thm 3.6 |
| dual picture `G_a^d x G_m`: obstruction subgroups are `V x {1}`, `V x G_m`; counts `Z`-free for nondegenerate `Z` | **proved** | Lemma 4.2, Prop 4.3 |
| a dual zero estimate decisive at `Z` must prove ZE^all (all nondegenerate points; `delta`- and height-free) | **proved** | Prop 4.3 |
| dual zero estimate with degrees: vacuous for every slab box (`m - 1 >= 2T + 1`) | **proved** | Prop 4.4 |
| dual zero estimate with Hilbert functions (HS): false at cheap-subtorus points | **proved** for the verified `(n,m)`; all-`n` argument a sketch | Prop 4.5 |
| `Z` on an obstruction subtorus: impossible (nondegenerate); quantitative distance governed by the heaviness `X_u` of the defining characters | **proved** | Props 5.1, 5.2 |
| ZE-col\*, MIN^pair, `(R_H)` at the actual point | **open** (not derived) | — |
| growth `M = o(log R/loglog R)`, or uniform bound | **not proved** | — |

### 6.2 Why the zero estimate cannot close, in one paragraph

A zero estimate converts *many points* (a group orbit, with multiplicities) and *one polynomial* into an obstruction
subgroup. Route 5 has *one point* and *a huge linear system*.

* In the torus picture, the arithmetic of `Z` (Liouville at `Z^j`) produces at most `1.07 sqrt(n)/lambda` further
  points. Philippon needs about `d`. Even an idealised count needs more than the Siegel slack of route 5's own box,
  which is a factor `~8`. So the conclusion is satisfied by the trivial subgroup (Prop 3.5, Thm 3.6).
* In the dual picture, the system becomes one exponential polynomial vanishing on `|A| >= 3^d/2` points. That is
  the right shape. But for nondegenerate `Z` every obstruction count is independent of `Z` (Prop 4.3). So the zero
  estimate could only prove the `delta`-free statement ZE^all for **all** nondegenerate points of the torus,
  including every fake position. Real (degree-normalised) estimates do not reach it (Prop 4.4), and the lossless
  form is false (Prop 4.5).

The only features of the actual point that zero-estimate theory can see are nondegeneracy (used, Prop 5.1) and the
group `<Z>` (exhausted after `O(sqrt n)` steps). Its height, its field of definition `Q(i)`, its size `e^{-W/4}`
and its pair data are invisible to the obstruction inequality. This is the zero-estimate face of R16 §5.1:
"arithmetic is spent; a box-level proof works only if no admissible regime box, fake ones included, has a near-zero
in `G`".

### 6.3 Variants considered and why they do not change the verdict

* **Multiplicity estimates at `e` alone** (`ord_e P_c >= m`). This is the singleton case of Prop 3.5. The trivial
  subgroup obstructs. Moreover a nonzero form with `ord_e >= mu` exists for every `mu <= rho(A)` by the definition
  of `rho(A)`, so no multiplicity estimate at `e` can say anything new.
* **Other points.** Do other points carry certified vanishing? Lemma 3.1 needs a tuple of Gaussian integers of one
  common norm whose relative angles are tiny. Within one cluster these are the powers `z^j` (and their conjugates).
  Torsion twists `i^{s_x} z_x` have large relative angles. Other arcs of the same circle are not near `e`.
  Subtuple points (round 25) live on `T_5`, not on `G_m^d`. So `Sigma ⊂ <Z>` is forced. (This is a complete list
  for Liouville-type certification at the archimedean place; it is **heuristic** as a statement about all
  conceivable arithmetic inputs.)
* **`p`-adic zero estimates.** The `p`-adic analogues use the valuation vectors `a'_p` as points. These lie **on**
  cheap subtori for Sylvester profiles (R11 J1). That is exactly why SL and SL\* failed, and Prop 5.1 has no
  `p`-adic counterpart.
* **Small value estimates** (Roy; Laurent–Roy on `G_a x G_m`; [memory-quote], shape only). They need a *sequence*
  of polynomials, one for each degree, each with small values at a point, and conclude algebraicity. Route 5 has one
  degree with exact vanishing, and the conclusion would again be `Z`-blind in the sense of Prop 4.3. No gain.
* **Ideal-theoretic estimates** (Nesterenko, Philippon's algebraic-independence criteria). These see the whole
  system, but their output is a degree/height bound for isolated zeros, i.e. arithmetic Bezout. That gives
  `h(Z) <~ D^{d-1} h(V_H)`, which is vacuous against `h(Z) ~ dW/2`, as already noted in [D] §5.5 and R11 theory
  §4(a).
* **ZE-col\* (universal family).** Its relations are integer measures on `Z^n_0` pushed to the column space `L` of
  rank `rho_hat` (R11 subcluster Thm 3.1). The arguments of §3–§5 use only "one point, a linear system of Laurent
  forms, nondegeneracy relative to the ambient torus". They therefore apply verbatim on the `rho_hat`-dimensional
  torus of `L`, with `d` replaced by `rho_hat`. The numbers were not recomputed for that setting.

### 6.4 What a successful argument would have to be (heuristic, for the map)

It must be **multi-polynomial**: it must use `dim V_H ~ |A|/2`, since single forms provably carry no information
(Thm 3.6). It must also be **sensitive to more than nondegeneracy**, since otherwise it proves ZE^all, a statement
about all fake positions too, i.e. the barrier question itself (Prop 4.3; R12 fake Thm 2.1/2.2; R16 §5.1). In the
dual picture the correct target is ZE-tor^glob(A,m) or ZE^all(A,m) for overdetermined regime boxes. Zero-estimate
theory is not the tool for either, because its obstruction inequality is a counting inequality, and the counts here
are a factor `>= d!/C^d` short. I see no way to recover that factor; this is a judgement, not a theorem.

### 6.5 Owner constraints

* *No uniform exceptional sets from fixed-`S` Subspace Theorems.* None used.
* *No Gaussian gcd or denominator costs discarded.* Lemma 3.1 keeps the full forced-divisibility cost `j W Lambda/2`
  at every place, including inert primes and `1+i` (full powers). The derivative step tracks `|c^(alpha)|_1` and
  `T_A^{|alpha|}`. No normalisation by gcds occurs, because `Z = (z_x/z_n)` is used only through the integral forms
  `G_k(z^j)`.
* *No drifting common phase replaced by zero.* All evaluations are of degree-0 characters (`sum k = 0`), so
  `theta_0` cancels identically, as in [D] §1. No phase is fixed.
* *No template-dependent threshold treated as uniform.* See §5 item 5. The `lambda`-dependent bounds are stated with
  `lambda` explicit. The uniform versions use only `log 5`, `n`, `r`, `T`.

## 7. Evidence / scripts

All runs used `prlimit --as=3500000000`, one process at a time. Paths are relative to `round26/ze_checks/`.

| script → output | what | label |
| --- | --- | --- |
| `powers_valuation.py` → `powers_valuation_out.txt` | `v_pi(G_k(z^j)) = j(e D/2 + <k,a_p>)`, the conjugate identity, and `G_k(z^j) = N^{jD/2} Z^{jk}` on a 4-point cluster of `N = 5*13*17*29`, `j = 1,2,3`: 1620 checks, 0 failures | exact (valuations); floating (archimedean identity, rel. err. `< 1e-9`) |
| `spectral_cap.py` → `spectral_cap_out.txt` | (3.2) on 49 random profiles (`n = 4..7`, random weights and supports): `rho(S)` exact by Fraction elimination, `Omega(S)` exact, spectral bound floating; 0 failures; max `rho/cap = 0.19`. Sylvester `lambda = n/(n-1)` at `n = 4, 8, 16` | exact / evidence |
| `vacuity_numbers.py` → `vacuity_numbers_out.txt` | `m_max` (exact), `w`, `ceil(4w)`, `E`, `J_lemma` (Cor 3.2); margins of Prop 3.5(1), (2) and Prop 4.4 for `n = 2966 .. 10^6` | exact `m` (`n <= 16384`; proved floor `0.6417d - 2.2137` above); mpmath logs |
| `single_poly_barrier.py` → `single_poly_barrier_out.txt` | (3.5) in exact big integers for `n = 2966, 3000, 4096, 8192`, `mu = ceil(4w), m_max`, and the largest `mu` (`m_max + 1` or `+ 2`) | exact |
| `hs_refute.py` → `hs_refute_out.txt` | Prop 4.5: no connected `H ≠ G'` satisfies the Hilbert-sharp inequality, for every `dim V = 0..d`, at `n = 2966, 4096`, `m = ceil(4w), m_max`, worst case `L = 3`; tightest `ln(LHS/RHS) = 0.197` (hyperplanes) | exact |

**Not done / limits.**
* Philippon 1986 and the errata could not be re-read (network). §2 is a memory-quote. The conclusions use only the
  shape (2.1) with margins `>= e^{1.7*10^4}`.
* Prop 4.5 is verified for the stated `(n, m)` only. The all-`n` argument is a sketch.
* The ZE-col\* transfer (§6.3, last item) is argued, not recomputed.
* No regime-size base-locus computation was attempted. R11 and R12 show that small models are zonotopal or
  LP-subcritical, and the statements here do not need any.
