# Round 11: many rigidity relations at once (sub-cluster intersection of the forced-divisibility base loci)

Task. `determinant.md` Theorem 4.1 applies to every sub-cluster `J` of an actual `M`-point cluster (same `N`, `W`,
`r`, profile restricted to `J`) once `|J| - 1 >= 300 ln(4r)`, and to every admissible slab parameter and height.
So the position vector `delta in R^M` lies in the intersection of all `pi_J^{-1}(Z_{P_J})`. Is this combined base
locus contained in the degenerate set, or does it at least force `rank <= theta log W` (which by `weakest.md`
Thm 1.3 would give growth)? Analyse rigorously; test at small sizes against fakes and random positions.

Conventions as in `round10/determinant.md` (cited as **[D]**) and its referee report (**[RD]**). *Theorem /
Proposition / Lemma*: proved here with all quantifiers, or quoted with source. *Exact*: integer or mod-`p` linear
algebra. *Evidence*: floating point or finite experiments; no proof uses them. *Heuristic*: labelled. Scripts and
outputs: `subcluster_checks/` (each `*.py` writes the `*_out.txt` / `*.txt` named in §6).

## 0. Summary and verdict

**No growth improvement and no uniform bound.** The intersection of the sub-cluster rigidity relations was analysed
rigorously. The question "is it contained in the degenerate set?" reduces *exactly* to a zero estimate on each
cluster's own column torus (ZE-col). I prove structural facts that locate what the intersection can and cannot do,
one new inequality (Theorem 4.4), and a conditional growth theorem that covers all norms (Theorem 5.4).

1. **Universal admissibility (Lemma 1.1, proved).** Every finitely supported integer `c` on `Z^M_0` satisfying the
   LP inequality (1.2), with the order measured *relative to the column space* `L`, gives an exact relation
   `sum_k c_k e^{i<k,delta>} = 0`. All sub-cluster slab-box systems are special cases. `delta` is nondegenerate
   relative to `L` (Lemma 1.0).
2. **The intersection is small but uninformative (Props. 2.1-2.2, proved).** The combined locus has dimension
   `<= n - 2` for every `M`, and each further point is determined up to `2T` choices by `n - 1` others (unless a
   whole fibre lies in the locus). This says
   nothing about whether one actual (algebraic) point lies on it.
3. **Push-forward reduction (Thm 3.1, Prop 3.2, proved).** On `L`, every relation becomes a trigonometric polynomial
   whose coefficient measure is `hat pi_* c` on the exponent lattice. Nondegenerate sub-clusters give the same system in
   other coordinates. **However the sub-cluster systems are combined, any product-difference certificate obtained is
   a character certificate up to `O(log)` per unit**, so the referee's R1 (full box) mechanism never becomes
   available beyond the strength of the character inequalities (capped at `log R/loglog R`).
4. **Visibility (Thm 4.1, Cor. 4.2, proved).** If `rho_hat >= 300 ln 4r`, the cluster's column torus is not contained in
   the combined locus. If `M >= (rho_hat+1)(rho_hat+2)`, it is not contained in the universal locus (linear forms in the
   points, `D = 2`). **A cluster blind to the universal system has `M < (300 ln 4r + 2)^2 = O((loglog R)^2)`.** This
   complements [D] §5.1 ("low-rank clusters lie in `Z_P` automatically"), which describes one slab-box only.
5. **Signed mean-width bound (Thm 4.4, proved, new).** `reg(S) <= E_eps h_S(<eps,.>)` for every finite
   `S subset Z^rho`. This strengthens Theorem B of `threshold/construct.md` and is sharp on boxes and `l^1`-balls.
6. **The rank target (Prop. 5.1, Cor. 4.5).** Assume ZE-tor for every single sub-cluster system. Then every actual
   cluster with an admissible sub-cluster is *blind* to the combined locus. So "the combined locus forces
   `rank <= theta log W`" is **equivalent to a blindness statement B(theta)** about profiles: blind profiles have small
   rank. It has a sufficient form B'(theta) that is pure linear algebra on exponent data. Theorem 4.4 yields a sufficient criterion for blindness (ratio-2 criticality; Corollary 4.5:
   cross-polytope geometries). The examples I have are too small to have admissible sub-clusters, so **B(theta) is
   neither proved nor refuted**. On squarefree norms B is not needed ([D] Lemma 3.1).
7. **Conditional (Thm 5.4).** ZE-col for the universal system gives `M < (300 ln 4r + 2)^2 = O((loglog R)^2)` on
   **all** norms, for `R >= R_0(C)`. By Lemma 5.3, the part of ZE-col beyond ZE-tor for single slab-box systems is a zero estimate for
   linear forms in the points on low-rank tori. Remark 5.5 makes it concrete: actual clusters satisfy
   `z_x = Q(hat a_x - hat a_b)` with `Q` quadratic. This sub-question is finite and checkable, and I did not settle it.
8. **Computations (§6).**
   * *Exact.* Below the regime the honest combined system is **empty**: all 120 tested pair-separated sub-cluster
     systems have `V(A,m) = 0`. Cross-polytope geometries of rank 2, 3 have zero push-forward in every one of up to
     140 systems, and zero universal slab space. Without pair separation (cube profile), product certificates appear
     in single systems.
   * *Evidence.* On the skew-double-core and round-5 Paley-I profiles at `q = 7` (`M = 8`, `n = 4`, 280 systems,
     non-supercritical slab-boxes), 100/100 fake-type and random positions are certified off every single base
     locus, and the single and combined tangent cones have no zeros besides the diagonal. In a supercritical but non-pair-separated model, flatness of the tangent cones makes the
     float classification inconclusive.

| question | answer |
| --- | --- |
| combined base locus contained in the degenerate set? | equivalent to ZE-col (§5.2); not decided. Its dimension is `<= n-2` (Prop. 2.1), which is irrelevant for one point |
| does combining sub-clusters create new nonvanishing certificates? | not of product type beyond characters (Prop. 3.2); nondegenerate sub-clusters add nothing (Thm 3.1(c)) |
| combined locus forces `rank <= theta log W`? | not decided. Under ZE-tor it is equivalent to the combinatorial blindness statement B(theta) (Prop. 5.1); a blindness criterion is proved (Thm 4.4, Cor. 4.5), but no counterexample in range is known |
| what replaces the rank target? | blind-to-universal implies `M < (300 ln 4r + 2)^2` (Cor. 4.2, unconditional); ZE-col(univ) implies `M = O((loglog R)^2)` on all norms (Thm 5.4) |
| squarefree / near-extremal profiles | squarefree: rank `>= M/4 - 2` ([D] Lemma 3.1), so the rank target is moot and [D] Thm 5.7 stands. Near-extremal non-squarefree: B(theta) is open |
| small-size tests (`n = 4..6`, `M = 6..12`) | honest systems empty (exact); fakes certified off; non-honest models: no evidence against ZE, and inconclusive where supercritical |

## 1. Setting, the combined locus and the universal admissible system

### 1.1 Notation

A normalised cluster `z_1, ..., z_M` (`round8/weakest.md` §1.1): `|z_x|^2 = N`, `N` odd with all prime factors
split, `N = prod_{p} p^{e_p}` (`r` primes), arc constant `C`, `C' = max(C, 1)`, `W = log N`, angular width
`Theta <= C e^{-W/4}`, relative angles `delta_x = theta_x - theta_0` with `|delta_x| <= Theta/2`. Exponents
`a_{p,x} in [0, e_p]`, columns `beta_p = (2a_p - e_p 1)/e_p`, weights `omega_p = e_p log p / W` (`sum omega_p = 1`).
`phi_p = arg pi_p in (0, pi/2)` (choose `pi_p` in the open first quadrant).

**Windings.** By `round4/outside.md` (1.1) (as used in [D] Lemma 3.2), there are integers `kappa_x` with
`<2 a_x - e, phi> = theta_0 + delta_x + (pi/2) kappa_x`. Subtracting the equation for a fixed base point `b`,
`(pi/2)|kappa_x - kappa_b| <= 2 sum_p e_p phi_p + Theta < pi sum_p e_p + 1`, so `|kappa_x - kappa_b| <= 2 sum_p e_p + 1`.

**Extended exponents and the column space.** Put `hat a_x = (a_x, kappa_x) in Z^{r+1}` and
`hat phi = (2 phi, -pi/2)`. Then, exactly,

```text
delta_x = <hat a_x, hat phi> - c_0       (c_0 = theta_0 + <e, phi>),                              (1.1)
```

and for `k in Z^M_0 = {k : sum k = 0}`: `<k, delta> = <hat pi(k), hat phi>` with `hat pi(k) := sum_x k_x hat a_x in Z^{r+1}`.
Let `L subset R^M` be the span of `1` and the `r + 1` columns of the `M x (r+1)` matrix `(hat a_x)`; let
`rho_hat = dim L - 1` (the affine rank of `{hat a_x}`; `rho_hat <= rank + 1`, where `rank` is the exponent rank of
`weakest.md` §1.1). `Lambda := hat pi(Z^M_0)`, a lattice of rank `rho_hat`.

**Lemma 1.0 (delta is nondegenerate in L).** `delta in L`, and for `v in Z^M_0`: `<v, delta> = 0` iff `hat pi(v) = 0`
(i.e. `v perp L`). So `L` is the smallest subspace of `R^M` that is spanned by `1` and integer vectors and contains
`delta`, and `delta` satisfies no integer relation beyond `L^perp`.

*Proof.* `delta in L` by (1.1). If `hat pi(v) = 0` then `<v, delta> = <hat pi(v), hat phi> = 0`. Conversely let
`u = sum_x v_x a_x in Z^r`, so `<v,delta> = 2<u,phi> - (pi/2) <v,kappa>`. If `u != 0`, `outside.md` Lemma 0.1 gives
`dist(<u,phi>, (pi/4)Z) > 0`, so `2<u,phi> notin (pi/2)Z` and `<v,delta> != 0`. If `u = 0`, then
`<v, delta> = -(pi/2)<v, kappa>`, which is `0` iff `<v,kappa> = 0`, i.e. iff `hat pi(v) = 0`. QED

**Relative order, LP value.** For a finitely supported integer vector `c` on `Z^M_0` let `hat pi_* c` be its push-forward
(a measure on `Lambda`), and

```text
ord_L(c) = sup { m : sum_k c_k P(hat pi(k)) = 0 for every polynomial P on R^{r+1} with deg P < m },
LP(c)    = sum_p omega_p h_{supp c}(beta_p) = (2/W) sum_p log p * h_{supp c}(a_p).
```

`ord_L(c) >= ord(c)` (the order of [D] §1), because `P o hat pi` is a polynomial of degree `deg P` in `k`.

### 1.2 The universal admissibility lemma

**Lemma 1.1 (any support, relative order).** Let `c != 0` be an integer vector on a finite subset of `Z^M_0`, let
`D` be even with `D >= |k|_1` on `supp c`, and `m = ord_L(c) >= 1`. If

```text
(W/2) LP(c)  <  (m/2) W - 2 ln|c|_1 - 2m ln(D C'/2) + 2 ln m!                                         (1.2)
```

then `F_c(z) = 0`, equivalently `g_c(delta) := sum_k c_k e^{i <k, delta>} = 0`. Call such `c` *admissible*.

*Proof.* This is the chain of [D] Theorem 4.1 with `ord` replaced by `ord_L`. Put `x_k = <k,delta> = <hat pi(k), hat phi>`.
For `j < m`, `k -> x_k^j = P_j(hat pi(k))` with `P_j(y) = <y, hat phi>^j` of degree `j`, so `sum_k c_k x_k^j = 0`.
Hence `|g_c(delta)| = |sum_k c_k (e^{i x_k} - sum_{j<m} (i x_k)^j/j!)| <= |c|_1 (D Theta/2)^m / m!`, using
`|x_k| <= |k|_1 max_x |delta_x| <= D Theta/2`. Suppose `F_c(z) = N^{D/2} g_c(delta) != 0`. By forced divisibility
([D] (I1), Theorem C1 of `threshold/construct.md`), `Norm F_c(z) >= N^D prod_p p^{-h_{supp c}(a_p)} = N^D e^{-(W/2) LP(c)}`.
With `|F_c|^2 = Norm F_c` and `ln Theta <= ln C' - W/4`:
`-(W/2) LP(c) <= 2 ln|c|_1 + 2m ln(DC'/2) - mW/2 - 2 ln m!`, contradicting (1.2). QED

**The two loci.** Fix an actual cluster in the regime of [D] (4.1) for `n = M`.

* *The combined (sub-cluster) locus* `Z_comb`: the set of `delta' in R^M` with `g_c(delta') = 0` for every `c` in
  `V_H(J)` (extended by zero from `Z^J_0` to `Z^M_0`), over all sub-clusters `J` with `|J| - 1 >= 300 ln(4r)` and
  `|J| >= 64`, every choice of base point in `J`, with the slab-box `(T, w, m, H)` of [D] Theorem 2.1 built from the
  profile restricted to `J`. This is the object of the task. Other slab parameters `(T, w, H)` are covered whenever
  their vectors satisfy (1.2); they then belong to the universal family below.
* *The universal locus* `Z_univ`: the common zeros of `g_c` over `c` admissible (1.2) **and** all `c` above.

**Lemma 1.2.** `delta in Z_univ subset Z_comb`.

*Proof.* [D] Theorem 4.1 applies to each sub-cluster `J` (a cluster on the same arc and the same `N`, with `r`
columns restricted to `J`; its hypotheses (4.1) hold with `n = |J|`). Lemma 1.1 handles the admissible `c`. The
inclusion is by definition. QED

Remarks. (i) Sub-cluster systems are genuinely new relations: for `J subset J'` with a common base point, `A_J`
embeds into `A_{J'}` (when `T_J <= T_{J'}`, `w_J <= w_{J'}`), but `V(A_J, m_J)` is not contained in
`V(A_{J'}, m_{J'})` because `m_J < m_{J'}`. (ii) In the slab-box world everything depends on the columns only
through `LP`; Lemma 1.1 is the general form, and the task's intersection over `(J, T, w, H)` is a subfamily of it.

## 2. What intersecting over sub-clusters does to the base locus

Work modulo `R 1` and `2 pi`: positions live on the torus `T^{M-1} = (R/2piZ)^M / (diagonal)`. Fix one size `n`
with `n - 1 >= 300 ln(4r)`, `n >= 64`, `n <= M`, and let `Z_comb^{(n)}` be the intersection over the `n`-subsets `J`,
base points and admissible parameters (finitely many choices of `(T, w, H)` suffice for what follows).

**Proposition 2.1 (the combined locus is small-dimensional).** `Z_comb^{(n)}` is a real algebraic subset of
`T^{M-1}` (in the coordinates `cos, sin` of `delta_x - delta_b`) of dimension at most `n - 2`, for every `M >= n`.

*Proof.* Each `Z_{P_J}` is the common zero set of the finitely many trigonometric polynomials `g_c`, `c` in a basis of
`V_H(J)`, hence real algebraic in `T^{n-1}`. It is a proper subset ([D] Proposition 4.2: Lebesgue null), so
`dim Z_{P_J} <= n - 2`. Finite intersections of preimages under the coordinate projections `pi_J` are real
algebraic. Suppose `dim Z_comb^{(n)} >= n - 1`. A real algebraic set of dimension `d` contains a point at which it is a
`C^infty` (Nash) submanifold of dimension `d`; let `T_0` be its tangent space there. The `M - 1` functionals
`delta -> delta_x - delta_b` (`x != b`, `b` fixed) form a basis of the dual of the tangent space of `T^{M-1}`, so their
restrictions span `T_0^*`; choose `n - 1` of them, `x_1, ..., x_{n-1}`, linearly independent on `T_0`. For
`J = {b, x_1, ..., x_{n-1}}`, the projection `pi_J` restricted to the submanifold is a submersion at that point, so
`pi_J(Z_comb^{(n)})` contains a nonempty open subset of `T^{n-1}`. But `pi_J(Z_comb^{(n)}) subset Z_{P_J}`, which has
empty interior. Contradiction. QED

**Proposition 2.2 (one coordinate over the others).** Let `J_0` be an `(n-1)`-subset containing the base point `b`,
`x notin J_0`, `J = J_0 cup {x}`. For fixed `delta_{J_0}`, the set of `zeta = e^{i delta_x}` with
`(delta_{J_0}, delta_x) in Z_{P_J}` is either the whole circle or has at most `2T` elements.

*Proof.* For `c in V_H(J)`, `g_c = sum_{j=-T}^{T} b_j(c) zeta^j` with `b_j(c) = sum_{k in A, k_x = j} c_k e^{i<k_{J_0}, delta_{J_0}>}`
(`|k_x| <= T` since `x` is not the base point). A nonzero Laurent polynomial of this shape has at most `2T` zeros on
the circle. QED

**What these give, and what they cannot give.** For `M >> n` the combined locus is tiny: its dimension does not grow
with `M` (Proposition 2.1), and each further point is determined up to `2T` choices by `n - 1` others, unless a whole
fibre lies in the locus (Proposition 2.2). But an actual cluster is **one point**, and one point is always
zero-dimensional. For actual clusters the content of Proposition 2.2 is also automatic: given `N` and `a_x`, the point
`z_x` is already determined up to a unit. Small dimension of `Z_comb` therefore carries no information on whether a
given actual `delta` can lie on it. A contradiction needs, as before, one `c` with `g_c(delta) != 0`.

So the question is not the size of `Z_comb` but its *structure near the identity relative to rational subtori*. The
right coordinates for this are those of `L` (Lemma 1.0), because `delta` is a nondegenerate point of `L`.

## 3. Push-forward to the column torus

Points of `L` are `delta'(s) = (<hat a_x, s>)_x + s_0 1` with `s in R^{r+1}`, `s_0 in R` (the actual cluster is
`s = hat phi`, `s_0 = -c_0`). For a measure `mu` on `Z^{r+1}` put `ghat_mu(s) = sum_y mu_y e^{i<y,s>}` and
`LP(mu) = sum_{p <= r} omega_p (2/e_p) h_{supp mu}(y_p)` (the `kappa`-coordinate carries no weight).

**Theorem 3.1 (push-forward).** For every finitely supported `c` on `Z^M_0`:

* (a) `g_c(delta'(s)) = ghat_{hat pi_* c}(s)` for all `s` and `s_0`; hence
  `L subset {g_c = 0}` iff `hat pi_* c = 0`.
* (b) `ord_L(c) = ord(hat pi_* c)`, the order of vanishing of the measure `hat pi_* c` on `Z^{r+1}`; and
  `LP(hat pi_* c) <= LP(c)`.
* (c) If `J` is a sub-cluster whose vectors `{hat a_x : x in J}` are affinely independent, then `hat pi` is injective on
  `Z^J_0`, and for `c` supported in `Z^J_0`: `hat pi_* c = 0` iff `c = 0`, and `ord_L(c) = ord(c)`.

So the combined system restricted to `L` is the family of trigonometric polynomials on the torus `R^{r+1}/2piZ^{r+1}`
(effectively on `L`, of dimension `rho_hat`) with coefficient measures `hat pi_* c`, all supported in `Lambda`,
each of order `>= ord` and `LP <=` the `LP` of the support it came from.

*Proof.* (a) `<k, delta'(s)> = <hat pi(k), s> + s_0 sum_k = <hat pi(k), s>` on `Z^M_0`; group the terms by
`y = hat pi(k)`. Distinct `y in Z^{r+1}` give linearly independent characters of `s in R^{r+1}`, so the sum vanishes
identically iff all coefficients vanish. (b) The first claim is the definition of `ord_L`. For `p <= r` and
`k in Z^M_0`, `<k, beta_p> = (2/e_p) <k, a_p> = (2/e_p) hat pi(k)_p`, and `supp hat pi_* c subset hat pi(supp c)`, so each
width can only shrink. (c) Affine independence of `{hat a_x}_{x in J}` means `sum_{x in J} k_x hat a_x = 0` with
`sum k_x = 0` forces `k = 0`. A linear injection `Z^{|J|-1} -> Z^{r+1}` has a linear left inverse on its image, so
polynomials of degree `< m` correspond in both directions. QED

**Consequence.** For a *nondegenerate* sub-cluster (affinely independent exponents) the push-forward is the same
system in other coordinates: no gain and no loss. All the new information that intersecting over *degenerate*
sub-clusters could bring lives in the push-forward images `hat pi_* V_H(J)`, which are sums of the coefficients of `c`
over the fibres of `hat pi`.

### 3.1 Product push-forwards are character certificates

The only mechanism known to prove a statement of ZE type is the referee's Proposition R1 ([RD] §5): when the box is
not cut, `V_H` is spanned by product differences `Delta^alpha` with `|alpha| >= m`, and `g` of a product difference is
`prod (zeta_i - 1)^{alpha_i}`. One might hope that intersecting many sub-cluster systems generates, after
push-forward, product differences on `Lambda` even though no single slab-box contains one. The next proposition shows
that such a certificate, however it arises, is no stronger than characters.

**Proposition 3.2.** Let `v_1, ..., v_d in Z^M_0` with `lambda_i = hat pi(v_i)` linearly independent, `alpha in N^d`,
`|alpha| >= 1`, `y_0 in Lambda`, and

```text
mu_alpha = sum_{0 <= s <= alpha} (-1)^{|alpha| - |s|} binom(alpha, s) * delta_{y_0 + s_1 lambda_1 + ... + s_d lambda_d}.
```

Let `c` be any integer vector on `Z^M_0` with `hat pi_* c = mu_alpha` (for instance an element of a sum of
sub-cluster systems). Then

* (i) `g_c(delta) = e^{i<y_0, hat phi>} prod_i (zeta_i - 1)^{alpha_i}`, `zeta_i = e^{i<v_i, delta>}`;
* (ii) `ord_L(c) = |alpha|` and `sum_i alpha_i X_{v_i} <= (W/2) LP(c)`, where `X_v = sum_p log p |<v, a_p>|`;
* (iii) hence if `c` is admissible (1.2), then
  `sum_i alpha_i (X_{v_i} - W/2) < -2 ln|c|_1 - 2|alpha| ln(DC'/2) + 2 ln |alpha|!`.

*Reading.* By (i), Lemma 1.1 and Lemma 1.0, an admissible `c` of this kind would force `<v_i, delta> in 2 pi Z` for some `i`, hence
(when `|v_i|_1 Theta < 4 pi`) `<v_i, delta> = 0` and `lambda_i = 0`. So for an actual cluster such `c` cannot exist, and a profile that admits one carries no actual cluster.
Whether a profile admits one is decided by (iii). By Stirling, the right side of (iii) is at most
`2|alpha| ln(2|alpha|/(e D C')) + O(ln|alpha|)`. So the `alpha`-weighted average of the character sums `X_{v_i}` must lie
below `W/2 + 2 ln(2|alpha|/(eDC')) + o(1)`. Every actual cluster, and every fake, has
`X_v >= W/2 - 2 log(C|v|_1/(2 sqrt 2))` whenever `<v,delta> != 0` ([RD] §6, corrected [D] Prop. 5.3(iii)). Hence
product push-forward certificates, however the sub-cluster systems are combined, operate only in an `O(log)` window
around the character threshold `X_v = W/2`, which is the window in which single characters already decide. This is
[D] Proposition 5.3 extended from parallelepiped supports in `k`-space to push-forwards of arbitrary combinations.
Single characters are capped at `log R/loglog R` by the refereed barriers (`round8/README.md` §2.4).

*Proof.* (i) By Theorem 3.1(a), `g_c(delta) = ghat_{mu_alpha}(hat phi)`, and
`sum_{s <= alpha} (-1)^{|alpha|-|s|} binom(alpha,s) prod_i w_i^{s_i} = prod_i (w_i - 1)^{alpha_i}` with
`w_i = e^{i<lambda_i, hat phi>} = e^{i<v_i,delta>}`. (ii) The points `y_0 + sum s_i lambda_i` are distinct (independence),
so `supp mu_alpha` is the whole lattice parallelepiped `Pi`. The measure annihilates every polynomial of degree
`< |alpha|` (iterated differences) and pairs nontrivially with `prod_i ell_i^{alpha_i}`, where `ell_i` are linear forms on
`span(lambda)` dual to the `lambda_i`; so its order is exactly `|alpha|`, and `ord_L(c) = |alpha|` by Theorem 3.1(b).
For `p <= r`: `h_{supp c}(a_p) >= h_{Pi}(y_p) = sum_i alpha_i |lambda_{i,p}| = sum_i alpha_i |<v_i, a_p>|`, because
`Pi = supp hat pi_* c subset hat pi(supp c)`. Multiply by `log p` and sum over `p`. (iii) Insert (ii) into (1.2) with
`m = |alpha|`. QED

**Illustration (exact; `pushforward.py`, `pushforward_out.txt`).** Take the *cube profile* (points `{0,1}^3`, three
columns, `e_p = 1`, equal weights), which is **not** pair-separated: adjacent vertices differ in one column of weight
`1/3`, so `||e_i|| := sum_p omega_p |<e_i, beta_p>| = 2/3 < 1`. For sub-clusters of size `n = 5, 6` with LP-supercritical
model parameters (`m > LP`), every single push-forward `hat pi_* V(A_J, m)` already equals the full space
`V(y-box, m)` and contains every product difference with `|alpha| = m` (dims 4 and 1 at `m = 5, 6`; 280 and 168
systems). Without pair separation, product certificates (and hence the R1 mechanism) appear immediately. With pair
separation, Proposition 3.2 confines them, in every combination, to the `O(log)` window of the character inequalities.

## 4. Blindness and visibility of a cluster's column torus

Call the cluster **blind** to a family of relations if its column space lies in their common zero set (`L subset Z_comb`,
resp. `L subset Z_univ`). By Theorem 3.1(a) this means `hat pi_* c = 0` for every `c` in the family: every relation is
then satisfied by *every* point of `L`, so it carries no information about `delta`. Otherwise call it **visible**.
Visibility does not exclude `delta`; it says that `delta` lies in a *proper* analytic subset of `L`.

### 4.1 Blind clusters are small for the universal system

**Theorem 4.1 (visibility criteria).** Let a normalised cluster satisfy `W >= 16 ln M + 8 log^+ C + 14`.

* (a) *(nondegenerate sub-cluster.)* If `rho_hat >= 300 ln(4r)`, the cluster is visible to `Z_comb`.
* (b) *(linear forms in the points.)* If `M >= (rho_hat + 1)(rho_hat + 2)` and
  `W > 8 ln M + 8 ln(2W) + 12 ln C'`, the cluster is visible to `Z_univ`. The witness is an admissible `c`
  supported on `{e_x - e_b : x}` (`D = 2`) with `ord_L(c) >= 3`, for which `F_c(z) = conj(z_b) sum_x c_x z_x`.

*Proof.* (a) Choose a maximal affinely independent subset `J_0` of `{hat a_x}`; `|J_0| = rho_hat + 1 >= 64`. [D]
Theorem 4.1 applies to the sub-cluster `J_0` (`|J_0| - 1 >= 300 ln 4r` and `W >= 16 ln |J_0| + 8 log^+ C + 14`), and
`V_H(J_0) != 0` ([D] Thm 2.1(c)). By Theorem 3.1(c), every `0 != c in V_H(J_0)` has `hat pi_* c != 0`.

(b) Fix `b`. Let `V` be the affine span of `{hat a_x - hat a_b}`, of dimension `rho_hat`, and choose `rho_hat`
coordinates of `Z^{r+1}` whose restrictions to `V` are a basis of `V^*`. A polynomial of degree `<= 2` on `R^{r+1}`
restricted to `V` is a polynomial of degree `<= 2` in these coordinates, so `ord_L(c) >= 3` for `c = (c_x)` on the
support `{e_x - e_b}` amounts to `E = binom(rho_hat + 2, 2)` integer linear equations in `M` unknowns, with
coefficients bounded by `B = Y^2`, `Y = max_{x,j} |hat a_{x,j} - hat a_{b,j}| <= 2 sum_p e_p + 1 <= 2W/ln 5 + 1`
(§1.1; `e_p log p <= W`, `p >= 5`). Since `M >= 2E`, the pigeonhole (Siegel) argument of [D] Thm 2.1(c) with
`H = 2MB + 1` (then `(H+1)^{M-E} > (2MB+1)^E`) gives an integer `c != 0` with `|c|_inf <= 2MB + 1`, so
`|c|_1 <= 3 M^2 B <= 3 M^2 (2W)^2` (`2W/ln 5 + 1 <= 2W` for `W >= 2`). Here `D = 2` and
`LP(c) <= sum_p omega_p (max_x beta_{p,x} - min_x beta_{p,x}) <= 2`, so the left side of (1.2) is `<= W`. Its right
side, `f(m) = m(W/2 - 2 ln C') + 2 ln m! - 2 ln|c|_1`, is increasing in `m`, and
`f(3) >= 3W/2 - 6 ln C' - 2 ln(3 M^2 (2W)^2) + 2 ln 6 > W` under the stated bound on `W`. So `c` is admissible, and
`g_c(delta) = 0`. Finally `hat pi_* c = sum_x c_x delta_{hat a_x - hat a_b}` with distinct points (distinct exponent
vectors, `weakest.md` Thm 1.3 step 1), so `hat pi_* c != 0`. QED

**Corollary 4.2 (blind clusters are small).** Under the hypotheses of Theorem 4.1 (both parts):

* if `L subset Z_comb` then `rho_hat < 300 ln(4r)`;
* if `L subset Z_univ` then `rho_hat < 300 ln(4r)` **and** `M < (rho_hat + 1)(rho_hat + 2) < (300 ln(4r) + 2)^2`.

Since `r <= W/ln 5`, a cluster blind to the universal system has `M = O(log^2 W) = O((loglog R)^2)`. The bound on `W`
in Theorem 4.1(b) holds for `R >= R_0(C)`, because `M <= W` there (`two_thirds/proof.md`).

This corrects the reading of [D] §5.1 ("low-rank clusters are invisible"): that is a property of one slab-box. Linear
forms in the points see every cluster with `M >= (rho_hat+1)(rho_hat+2)`. (*Heuristic remark.* Forms of degree `d` in the
points, i.e. supports `{sum_{j<=d} e_{x_j} - d e_b}` with `LP <= 2d`, push forward to the sumset `d S`,
`S = {hat a_x - hat a_b}`. By Khovanskii and the Ehrhart bound `rho_hat! vol(conv S) >= M - rho_hat`, for large `d` one
expects `|dS| ~ (M - rho_hat) d^{rho_hat}/rho_hat!` and hence admissible order `m ~ d (M - rho_hat)^{1/rho_hat}`, i.e. LP ratio
`~ (M - rho_hat)^{1/rho_hat}/2`. So low-rank clusters with `M >> 2^{rho_hat}` should be seen *supercritically*. Not proved:
the Khovanskii threshold is not controlled.)

### 4.2 Blindness criteria for the sub-cluster family

**Proposition 4.3 (cheap column subtori lie in the combined locus).** For every set `S` of split primes with
`sum_{p in S} e_p <= 3`, the subtorus `K_S = {(sum_{p in S} t_p a_{p,x})_x + t_0 1}` lies in `Z_comb`.

*Proof.* For every admissible sub-cluster `J` and slab-box, `h_{A_J}(a_p|_J) = (e_p/2) h_{A_J}(beta_p|_J) <= e_p w_J`, so
`sum_{p in S} h_{A_J}(a_p) <= 3 w_J < 4 w_J <= m_J` ([D] Thm 2.1(d)). [D] Proposition 5.2 (with `kappa = 0`) gives
`pi_J(K_S) subset Z_{P_J}`. QED

The next theorem identifies a whole class of blind column geometries. For finite `S subset Z^rho` let
`reg(S)` be the least `d` such that every function on `S` is the restriction of a polynomial of degree `<= d`, and

```text
Omega'(S) := E_eps h_S(<eps, .>),     eps uniform in {+1,-1}^rho.
```

**Theorem 4.4 (signed mean-width bound).** For every `rho >= 0` and every finite nonempty `S subset Z^rho`,
`reg(S) <= Omega'(S)`. Equality holds for boxes `prod [0, a_i]` (both sides equal `sum a_i`) and for `l^1`-balls
`{|y|_1 <= R}` (`Omega' = 2R`, and an axis segment of `2R + 1` points gives `reg >= 2R`).

This strengthens Theorem B of `threshold/construct.md` (`reg <= MW = 2 E_P w_P`): writing `2 1_P = eps + 1`,
`MW(S) = E_eps (1/2)[h(<eps + 1,.>) + h(<1 - eps,.>)] >= E_eps h(<eps,.>) = Omega'(S)` by subadditivity of widths, and
the inequality is strict on `l^1`-balls (`reg = 2R = Omega'`, `MW = 3R` for `rho = 2`).

*Proof.* The induction of Theorem B, with signed instead of `0/1` averages. `rho = 0`: a point, `reg = 0 = Omega'`.
For `rho >= 1` slice by `y_1 = c`; each nonempty slice `S_c` is a translate of a set `S'_c subset Z^{rho-1}`, and by
induction `reg(S_c) <= g_c := Omega'_{rho-1}(S'_c) = E_{eps'} l_c(eps')`, where `l_c(eps')` is the width of
`<eps', y'>` on `S_c`. Put `K = Omega'_rho(S) = E_{eps'} (1/2)[h_S(<eps',y'> + y_1) + h_S(<eps',y'> - y_1)]`.
*Key inequality:* for slices `c <= c'`, let `[a_1,a_2]`, `[b_1,b_2]` be the ranges of `<eps',y'>` on `S_c`, `S_{c'}`. Then
`h_S(<eps',y'> + y_1) >= (c' + b_2) - (c + a_1)` and `h_S(<eps',y'> - y_1) >= (a_2 - c) - (b_1 - c')`. Adding,
`h(.. + y_1) + h(.. - y_1) >= 2(c' - c) + l_c(eps') + l_{c'}(eps')`, and averaging over `eps'`,

```text
K >= (c' - c) + (g_c + g_{c'})/2.                                                                      (4.1)
```

With `c = c'` this gives `g_c <= K`. *Counting:* if `C_v = {c : g_c >= v}` is nonempty, (4.1) at its extreme
elements gives `max C_v - min C_v <= K - v`, so `|C_v| <= K - v + 1`. *Peeling* (`threshold/construct.md` §1.1,
Lemma 2.1, refereed: if the slices are `S_(1), S_(2), ...` in any order, interpolate on `S_(1)` by `P_1`, then write
`P = P_1 + (y_1 - c_(1)) P'` and recurse on the remaining slices) gives `reg(S) <= max_i (reg S_(i) + i - 1)`. Order
the slices by decreasing `g`; then `|C_{g_(i)}| >= i`, so `i <= K - g_(i) + 1` and
`reg(S) <= max_i (g_(i) + i - 1) <= K`. QED

**Corollary 4.5 (cross-polytope geometries are blind).** Suppose `Lambda = hat pi(Z^M_0)` has a lattice basis in whose
coordinates `y` the columns include, for every `eps in {+-1}^rho_hat` modulo sign, a prime `p_eps` with `e = 2` and
`<k, a_{p_eps}> = <eps, y(k)>`, and with weight `omega_{p_eps} >= gamma 2^{1 - rho_hat}` for some `gamma > 1/2`. (This
requires the winding coordinate to depend on the columns on the configuration, `rho_hat = rank`; [RD] §6 shows this
by Cramer's rule when the rank is small.) Then every slab-box sub-cluster system pushes forward to zero, so
`L subset Z_comb`.

*Example.* Take the points `0, +-u_1, ..., +-u_rho` (`M = 2rho + 1`) and `r = 2^{rho - 1}` primes, all with `e_p = 2`
and `a_{p,x} = <eps_p, x> + 1`. Then `omega_p = log p / sum_q log q`, so `gamma > 1/2` holds as soon as all the primes
lie in `[P, P^{1.9}]`. With equal weights the profile is pair-separated with equality:
`sum_p omega_p |<v, beta_p>| = E_eps |<eps, y(v)>| >= max_i |y_i(v)| >= 1` for `y(v) != 0`. With unequal weights,
pairs `e_i - e_j` lose a fraction bounded by the weight imbalance.

*Proof.* For `c in V_H(J)`, [D] Thm 2.1 gives `ord_L(c) >= ord(c) >= m_J >= 4 w_J` and `LP(c) <= 2 w_J`. Suppose
`mu = hat pi_* c != 0`. Here `(2/e_p) h(y_{p_eps}) = h(<eps, .>)`, so
`LP(mu) >= gamma sum_{eps mod +-} 2^{1-rho} h_{supp mu}(<eps,.>) = gamma Omega'(supp mu)`. By Theorem 3.1(b) and
Theorem 4.4, `4 w_J <= ord(mu) <= reg(supp mu) <= Omega'(supp mu) <= LP(mu)/gamma <= LP(c)/gamma <= 2 w_J/gamma`,
so `gamma <= 1/2`. Contradiction. QED

The exact computations of §6.2 (rank `2, 3`: the universal slab space `V(S_slab, m)` itself is `0`) agree.

**Caveat (scope of §4.2).** Proposition 4.3 and Corollary 4.5 are statements about the *locus* `Z_comb`. They do not
yet concern actual clusters to which the sub-cluster family applies. A cluster in a subtorus `K_S` of Proposition 4.3
has rank `<= 3`, hence `M <= (5.5^3 + 1)/2 < 84` (`weakest.md` Thm 1.3). The cross-polytope example of Corollary 4.5
has `M = 2rho + 1`, while `r = 2^{rho-1}` gives `300 ln(4r) > 200 rho`. In both cases `M < 300 ln(4r) + 1`, so the
cluster has **no** admissible sub-cluster, and `Z_comb` is the whole torus for it anyway. What survives is the
mechanism. A column geometry that is critical at ratio `2` in the sense of the proof of Corollary 4.5
(`ord(mu) <= reg(supp mu) < 2 LP(mu)` for every `mu != 0` on `Lambda`) makes the cluster blind to every slab-box
sub-cluster system, whatever its number of points. Theorem 4.4 is a tool for proving such criticality.

## 5. What the combined locus can and cannot force

### 5.1 The rank target reduces to a combinatorial blindness question

Hypotheses ZE-tor and ZE-col are defined in §5.2; Lemma 5.3(i)-(ii) is used here.

**Proposition 5.1 (reduction of the rank target).** Assume ZE-tor for every single admissible slab-box sub-cluster
system, at the given `(N, C)`. Then every normalised cluster with `M >= 300 ln(4r) + 1` and
`W >= 16 ln M + 8 log^+ C + 14` is **blind** to `Z_comb`: `L subset Z_comb`. In particular `rho_hat < 300 ln 4r`.
Consequently the following are equivalent for a given `theta`:

* (i) membership `delta in Z_comb` (together with Lemma 1.0) implies `rho_hat <= theta ln W` for every such cluster;
* (ii) **B(theta)**: every profile of an actual cluster with `M >= 300 ln(4r) + 1` whose column space is contained in
  `Z_comb` has `rho_hat <= theta ln W`.

*Proof.* By Lemma 5.3(i), `Z_comb` is near `e` a finite union of rational subtori through `e`. By Lemma 1.2 and
Lemma 1.0, `delta in Z_comb` is nondegenerate relative to `L`; by Lemma 5.3(ii) (its proof), `L subset Z_comb`. Then
Theorem 4.1(a) gives `rho_hat < 300 ln 4r`. For the equivalence: (ii) => (i) is immediate from blindness. If (ii)
fails for some actual cluster, its `delta in L subset Z_comb`, so membership in `Z_comb` cannot exclude it, and (i)
fails. QED

**Status of B(theta).** B(theta) mentions actual clusters only through their profiles and windings. A sufficient,
purely combinatorial version is **B'(theta)**: every pair-separated exponent profile with `M >= 300 ln(4r) + 1` points
whose column space (for any winding vector `kappa`) lies in `Z_comb` has rank `<= theta ln W - 1`. B'(theta) is a
statement about integer exponent data, slab-boxes and push-forward linear algebra; it involves no transcendence and no
positions. Blindness requires `hat pi_* V_H(J) = 0` for every
admissible sub-cluster `J`, every base point and every slab-box. The proof of Corollary 4.5 gives a sufficient
condition (ratio-2 criticality of the column geometry), and Theorem 4.4 can verify it. The only examples I have are
too small to have admissible sub-clusters (Caveat in §4.2). So **I can neither prove nor refute B(theta) for any
`theta < 1/ln 5.5`.** For `weakest.md` Cor. 1.4(b) one would need B(theta) with `theta < 0.587`. Three remarks:

* on squarefree norms B is not needed: [D] Lemma 3.1 forces `rho_hat >= M/4 - 2`, and Proposition 5.1 gives
  `M < 1200 ln 4r + 8` directly. This is [D] Theorem 5.7 under ZE-tor;
* for the *universal* system the analogue of B holds unconditionally in the strong form "blind implies
  `M < (rho_hat+1)(rho_hat+2)`" (Corollary 4.2);
* the rank is the less natural currency here. The statement that converts blindness into a bound is "blind implies
  few points".

**Remark 5.2 (what ZE alone gives).** By Lemma 1.0 the nondegenerate sub-cluster `J_0` of Theorem 4.1(a) has a
nondegenerate position. If `rho_hat >= 300 ln 4r`, [D] Theorem 4.1 puts `delta_{J_0}` on `Z_{P_{J_0}}`, which
contradicts ZE(`rho_hat + 1`). So **ZE alone gives `rho_hat < 300 ln(4r)` on all norms**, with no intersection over
sub-clusters. This is [D] Theorem 5.7's argument without Lemma 3.1. On general norms `300 ln 4r` is far above
`0.587 ln W`. Proposition 5.1 says that intersecting over sub-clusters lowers it further exactly to the extent that
blind profiles of intermediate rank do not exist.

### 5.2 Hypotheses relative to the column torus, and a conditional theorem on all norms

For a family `F` of relations (`F = comb` or `univ`) and a profile, consider near the identity:

* **ZE-tor(F):** in the ball `max_x |delta'_x - delta'_b| < 1/(4D)`, `Z_F` coincides with a finite union of rational
  subtori through the identity.
* **ZE-col(F):** if `L not subset Z_F`, then no `delta' in L` with `max|delta'_x| <= Theta` that is nondegenerate
  relative to `L` (no relation `v in Z^M_0 \ L^perp`) lies in `Z_F`.

**Lemma 5.3.** (i) ZE-tor for each member of a finite family implies ZE-tor for the intersection of their zero sets.
(ii) ZE-tor(F) implies ZE-col(F). (iii) For `L = R^M` (`rho_hat = M - 1`) and `M - 1 >= 300 ln 4r`, ZE(`M`) of [D]
implies ZE-col(univ).

*Proof.* (i) A finite intersection of finite unions of rational subtori through `e` is a finite union of
intersections of such subtori. Near `e`, a closed subgroup coincides with its identity component, a rational
subtorus. (ii) If `delta' in Z_F cap L` is nondegenerate relative to `L`, it lies in some rational subtorus
`K subset Z_F` through `e`. The smallest rational subspace containing `delta'` is `L`, so `L subset K subset Z_F`.
(Analytic continuation: `Z_F` is analytic and contains an open piece of the connected `K`, hence all of `K`.)
(iii) `Z_univ subset Z_P` for the slab-box of the whole cluster, and "nondegenerate relative to `R^M`" is
"nondegenerate". QED

**Theorem 5.4 (conditional, all norms).** Assume ZE-col(univ) for every profile in the regime of Theorem 4.1, and
`R >= R_0(C)`. Then every normalised cluster on an arc of length `C sqrt R` satisfies

```text
M  <  (rho_hat + 1)(rho_hat + 2)  <  (300 ln(4r) + 2)^2  <=  (300 ln(4W/ln 5) + 2)^2  =  O((loglog R)^2).
```

*Proof.* By Lemma 1.2, `delta in Z_univ`; by Lemma 1.0 it is nondegenerate relative to `L`; and
`max|delta_x| <= Theta/2`. If `L not subset Z_univ`, ZE-col(univ) is contradicted. So `L subset Z_univ`, and
Corollary 4.2 applies. QED

Compared with [D] Theorem 5.7 (ZE implies `M = O(log W)` on squarefree norms only), Theorem 5.4 needs no rank lemma
and covers all norms. Its hypothesis is different, not weaker. By Lemma 5.3, it splits into two parts:

* For clusters visible to `Z_comb`, ZE-col(univ) follows from ZE-tor for the single slab-box systems. That is the
  natural "components through `e` are subtori" form of ZE-real in `round10/README.md` §8.1.
* For clusters blind to `Z_comb` but visible to `Z_univ` (rank `< 300 ln 4r`, `M >= (rho_hat+1)(rho_hat+2)`), it is a
  zero estimate for the **non-slab** supports on a low-rank torus: linear forms in the points, and (heuristically,
  §4.1) sumset supports. Here the information sits in relations that the sub-cluster family does not contain.

**Remark 5.5 (the linear-forms part is concrete).** Under the bound on `W` of Theorem 4.1(b), every integer `c` in the
kernel `C_3` of the degree-`<= 2` moment map of `S = {hat a_x - hat a_b}` with `4 ln|c|_1 < W - 12 ln C' + 4 ln 6` is
admissible (the computation in the proof of Theorem 4.1(b)). By Bombieri-Vaaler, `C_3 cap Z^M` contains `M - E'`
independent vectors (`E' <= E = binom(rho_hat+2,2)` the rank of the moment map) with product of sup-norms
`<= (sqrt(M) B)^{E}`, hence a lattice basis of sup-norm `<= M (sqrt(M) B)^{E}`. For `W` large compared with `E ln(MB)`,
every basis vector is admissible, and by linearity every `c in C_3` gives `sum_x c_x z_x = 0`. So **`z_x = Q(hat a_x - hat a_b)` for one
complex quadratic polynomial `Q` on the affine span**: the second-order Taylor expansion of `z_x` in the exponent
coordinates is exact. (This is the classical "points on short arcs lie on low-degree curves" mechanism, here in
exponent space.) The corresponding piece of ZE-col asks whether a nondegenerate small `s` can make `y -> e^{i<y,s>}`
quadratic on `S`. At the tangent-cone level this asks whether `S` lies on a cubic hypersurface whose top-degree form is
`<y,u>^3` for an irrational `u`. Counting gives `2(M - E)` real conditions on `rho_hat` real unknowns, so the system is
overdetermined once `M > E + rho_hat/2`, which holds in Theorem 4.1(b)'s range. This sub-question is finite and
checkable, and I did not settle it.

## 6. Computations (exact where stated; otherwise evidence)

All scripts are in `subcluster_checks/`. `common.py` and `engine.py` hold the shared code: slab-boxes, mod-`p` ranks,
nullspaces and RREF with `p = 2147483629`, Chebyshev-basis projections, and tangent-cone residuals.

*What "exact" means here.* A rank over `F_p` is a lower bound for the rank over `Q`. So "`V(A, m) = 0`" and
"`V(S_slab, m) = 0`" are certified over `Q` when the evaluation matrix has full rank `|A|` (resp. `|S_slab|`) mod `p`.
A vanishing universal slab space forces every push-forward of order `>= m` supported in the slab region to vanish,
so in §6.2 the push-forward claims are certified over `Q` as well. Positive dimensions and memberships (§6.1) are
computed over `F_p`. Over `Q` they are evidence: they could change only if `p` divides certain minors.

### 6.1 Without pair separation, product certificates appear at once (over `F_p`; `pushforward.py`)

Cube profile (points `{0,1}^3`, `M = 8`, three columns, `e_p = 1`, equal weights, `||e_i|| = 2/3 < 1`). Sub-clusters
of size `n = 5, 6`, all base points, `T = 2, 3`, `w = 2` (`|y_i| <= 1`), `m = 5, 6 > LP = 4`. The output is in
`pushforward_out.txt`:

* for 280 (`n = 5`) and 168 (`n = 6`) systems, `dim hat pi_* V(A_J, m)` is `3`-`4` (`m = 5`) and `0`-`1` (`m = 6`);
* the combined dimension equals `dim V(y-box, m)` (`4`, resp. `1`);
* every product difference `Delta^alpha` with `|alpha| = m` in the `y`-box lies in a single push-forward.

This is the referee's R1 mechanism, available here because the profile violates pair separation (§3.1).

### 6.2 Pair-separated low-rank geometries are blind (exact; `pushforward2.py`, `crosspoly.py`)

`lowrank_profiles.py` (LP over weights, `lowrank_rho2.txt`) shows that rank-2 pair-separated profiles with points
in a `4 x 4` grid and column directions in `[-3,3]^2` have at most 5 points (`M = 5`: the checkerboard, which is the rank-2 cross-polytope, with min pair norm
exactly `1`). The output is in `crosspoly2_out.txt` and `crosspoly3_out.txt` (`pushforward2_out.txt` repeats the
`rho = 2` case). The setup:

* cross-polytope profiles with `rho = 2` (`M = 5`, `r = 2`) and `rho = 3` (`M = 7`, `r = 4`), `e = 2`;
* sub-clusters of size `n = 3, 4, 5`, all base points, `T in {2, 3}`, `w in {2, 3}`, `m in {5, ..., 8}` with `m > LP`;
* 5 to 140 admissible systems per parameter set.

The results:

* in every case `hat pi_* V(A_J, m) = 0` for **every** system, although `dim V(A_J, m)` (over `F_p`) is as large as 646;
* the combined space is `0`;
* the universal slab space `V(Lambda cap slab, m)` is `0`.

Theorem 4.4 proves this for all `rho`. On the slab region every `<eps, .>` has width `<= 2w`, so
`reg(S_slab) <= Omega'(S_slab) <= 2w < m` and `V(S_slab, m) = 0` (the argument of Corollary 4.5 with `gamma = 1`).

### 6.3 The signed mean-width bound (exact; `cross_reg.py`, `cross_reg_out.txt`)

`reg(S) <= Omega'(S)` was checked on 500 random sets (`rho = 2`: sizes 3-15 in `[0,3]^2`; `rho = 3`: sizes 4-19 in
`[0,2]^3`). There were no violations. The maximum of `reg/Omega'` was `1.000` for `rho = 2` and `0.889` for `rho = 3`.
Equality holds on boxes and `l^1`-balls (for example `rho = 3`, `R = 3`: `reg = 6 = Omega'`, `MW = 10.5`), and the
script asserts Theorem B (`reg <= MW`) on every set. Theorem 4.4 is proved in §4.2; this is a consistency check.

### 6.4 Honest small models have no supercritical sub-cluster systems (exact; `honest_small.py`)

Pair-separated full-rank squarefree profiles: binary codes with minimum distance `>= r/2`, equal weights,
`(M, r) = (6,6), (8,7), (8,8), (10,10), (12,12)`. Sub-clusters of size `n = 4, 5`, slab-boxes with `T in {3, 4}`,
`w in {4, 6, 8}`, and `m = floor(LP) + 1`. In all 120 tested systems `dim V(A_J, m) = 0` (`honest_small_out.txt`). So
**at the sizes `n = 4..6`, `M = 6..12` requested, the combined system of an honest profile is empty**: no relation is
forced, and `Z_comb` is the whole torus. The regime of [D] Theorem 2.1 starts at `n = 2137` (exact, `r = n`). Any
small-size test of the combined locus therefore has to use non-honest models, either not pair-separated (§6.1, §6.6)
or not LP-supercritical (§6.5). Their evidential value for ZE in the regime is limited.

### 6.5 Fake profiles with non-supercritical slab-boxes (evidence; `fakes_combined.py`, `fakes_single.py`)

The profiles are the skew double core at `q = 7` (`M = 8`, `r = 14`) and the round-5 Paley-I profile with `b = 8`,
`q = 7` (`M = 8`, `r = 64`), both from the referee's `ref_fakes.py` (code reused read-only). Sub-clusters have
`n = 4`, all bases (280 systems), `T = 3`, a median slab cut (as in the round-10 evidence scripts), and
`m = 6 = 2T` from `4 binom(m+2,3) <= |A|`. These models are **not** LP-supercritical.

* *Fake-type positions* `delta ~ U[0, 1e-4]^8` and random positions `N(0, 1e-3)` (50 each): all 100 are certified off
  **every** one of the 280 single base loci (`|leading Taylor term| > remainder`; skew double core).
* *Tangent cones.*
  * Single systems: 12 sampled, BFGS from 30 starts each. No zeros; `min tau = 9.6e-4` in every system.
  * Combined system: 30 BFGS starts on `S^6`. `min tau = 1.6e-3`, attained near the single-point-moved directions,
    which are not zeros.

  So in these models the tangent cone is trivial already for one system, and intersecting adds nothing.
* *Round-5 profile* (`fakes_combined_r5.txt`, same protocol; `|A|` 225-239, `dim V` 169-183). All 100 positions are
  certified off every single base locus. The combined tangent-cone residual has minimum `1.7e-3` over 30 starts,
  attained near single-point-moved directions and near two-cluster splits. No zeros.

### 6.6 A supercritical but non-pair-separated model: inconclusive (evidence; `supercrit_combined.py`, `single_tc.py`, `basepoint_check.py`)

Here `M = 6` with four random `0/1` columns (one a singleton column), `n = 4`, `T = 4`, `w = 6`, `m = 13 > LP`. There
are 52 admissible `(J, base)` systems with `dim V` between 12 and 79.

* *Flatness.* The single tangent cones vanish to very high order near the cheap rational strata. Near a sparse
  direction `e_x` (`x` not the base), `(A u)^m` differs from `Poly_{<m}|_A` only at order `>= 5` in the normal
  direction, because `k_x` takes `2T + 1 = 9` values: `k_x^{m-j} <k, eta>^j` has effective degree `<= 8 + j < m` for
  `j < 5`.
* *Effect on the numerics.* Float minimisation therefore returns spurious near-zeros at distance `~1e-2` from these
  strata. The PSLQ classification in `supercrit_combined_out.txt` ("93 of 131 single zeros without relation") is an
  artefact of this. Gauss-Newton refinement shows the minimisers within `1e-5`-`1e-7` of the sparse/column strata,
  not converged zeros.
* *Margin-constrained minimisation* over directions at relation-distance `>= eps` (`supercrit_combined_out2.txt`):
  * at `eps = 0.05`, the single systems give `min tau >= 5e-9` (or no feasible minimiser);
  * at `eps = 0.02`, they give values `~1e-13`, consistent with the flatness tail;
  * the combined system never produced a feasible minimiser.

  The data neither show nor exclude nondegenerate tangent-cone points of single systems, so they say nothing reliable
  about a gain from combination.
* `basepoint_check_out.txt`: base-point freedom removes some sparse directions (`tau = 2.3e-8` versus `~1e-29` when
  `x` is the base point) but not the column directions. This agrees with [RD]'s base-point asymmetry.

## 7. Conclusions and open points

**Answer to the question of the round.** Intersecting the rigidity relations of all sub-clusters (and all slab
parameters) gives a base locus of dimension `<= n - 2` independent of `M` (Proposition 2.1). That says nothing about
one actual point. Restricted to the cluster's own column torus `L`, where an actual position is nondegenerate, the
intersection is again a single system of the same kind, with supports in the push-forward lattice (Theorem 3.1).
Combining sub-clusters cannot manufacture the only known nonvanishing mechanism, product differences, beyond the
strength of the character inequalities (Proposition 3.2). Whether the combined locus is degenerate near the identity
is therefore still a zero estimate (ZE-col), not decided here in either direction. Two structural facts are new:

* under ZE-tor for single systems, every actual cluster in the range of the sub-cluster family is blind to the
  combined locus. So whether the intersection forces `rank <= theta log W` is exactly the combinatorial question
  B(theta) (Proposition 5.1). A sharp tool for blindness is proved (Theorem 4.4, a signed strengthening of Theorem B
  of `threshold/`), but no profile in range has been shown blind or visible;
* blindness is harmless once the relations are taken from all admissible supports. A cluster whose column torus lies
  in the universal locus has `M < (300 ln 4r + 2)^2` (Corollary 4.2). Hence ZE-col(univ) gives
  `M = O((loglog R)^2)` on **all** norms (Theorem 5.4).

**Open points.**

1. **The linear-forms zero estimate (Remark 5.5).** For a pair-separated exponent configuration `S` of affine rank `rho`
   with `M >= (rho+1)(rho+2)` points: can `y -> e^{i<y,s>}` be quadratic on `S` for a small `s` that is nondegenerate
   relative to `S`? First the tangent-cone version: does `S` lie on a cubic hypersurface with top form `<y,u>^3`, `u`
   irrational? The tangent-cone version is finite and checkable. A negative answer to the full (non-tangent-cone)
   question, together with ZE-tor for single slab-box systems, would give Theorem 5.4's conclusion
   `M < (300 ln 4r + 2)^2` (Proposition 5.1 and Theorem 4.1(b)).
2. **Sumset supports (heuristic remark in §4.1).** Make "low-rank clusters with `M >> 2^rho` are seen with LP ratio
   `~ (M - rho)^{1/rho}/2`" rigorous. This needs an effective Khovanskii threshold for `|dS|` and control of the
   heights. It would make the low-rank part of ZE-col a *supercritical* unlikely-intersection problem, like ZE itself.
3. **ZE-tor for one slab-box system** (components through `e` are rational subtori). It implies ZE and, by
   Lemma 5.3, ZE-col(comb). It is the component form of `round10/README.md` §8.1's ZE-real.
4. **The blindness question B(theta) (Proposition 5.1).** Are there pair-separated profiles with
   `M >= 300 ln(4r) + 1` points and rank in `(theta ln W, 300 ln 4r)` whose column torus lies in `Z_comb`? Ratio-2
   criticality (proof of Corollary 4.5, with Theorem 4.4) is a sufficient mechanism, and the cross-polytope examples
   have only `2rho + 1` points. If there are none for some `theta < 1/ln 5.5`, then ZE-tor for single systems gives
   `M = O((log R)^{1-eta})` on all norms via `weakest.md` Cor. 1.4(b). The question is pure linear algebra on profiles.
5. **No honest small-size test exists.** Below `n ~ 2137` the honest combined system is empty (§6.4), and non-honest
   models are either trivial (§6.5) or numerically too flat to classify (§6.6). Further progress on route 5 has to be
   theoretical, or use exact symbolic tangent-cone computations in models designed to be both pair-separated and
   supercritical. By Corollary 2.4 of [D] such models need `r << e^{0.0144 n}` and are out of reach at small `n`.
