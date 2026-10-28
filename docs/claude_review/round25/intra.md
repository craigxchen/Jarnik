# Round 25: intra-cluster multi-point Liouville on $T_5$ (shared primes)

*Round 25, write-up `intra`. Nothing here is refereed. Status: **reduces_to_route5**.*

## 0. Summary and verdict

**No growth improvement and no uniform bound.** `M <= (2/3+eps) log R/loglog R` stands.

**Verdict: the intra-cluster multi-point argument is the determinant method (route 5) in disguise.** It stops at
route 5's zero estimate in the growth regime, and on balanced clusters it is capped at the one-point value `15/4`
whatever is shared.

1. **One point, not many** (Lemma 2.1, proved, exact). The `binom(M,5)` gcd-reduced subtuple points
   `y_S in T_5(Q)` are the images `pr_S(Z)` of the single cluster point `Z = [z] in T_M(Q)`. Their heights, their gcds,
   their valuation vectors and the primes they share are all functions of the one profile `v_pi(Z)`. Explicitly,
   `h_B(y_S) = sum_p (max_S a - min_S a) log p`.
2. **Reduction** (Theorem 3.1, proved, exact). Every inference of the intra-cluster `n`-point Liouville method
   (`L_intra`: any blocks, any `F` on `T_5^n` chosen after the profile, Dyson-type weighted index, toric forced
   divisibility using the shared primes, integrality of the gcd-reduced tuples, `iota`-conjugates, group words) is
   implied by one-point toric Liouville (TL_M) at `Z` for the single Laurent polynomial `G = Phi^*F`. The two
   methods have the same value, and the widths satisfy `RHS_G <= RHS_F` with equality absent cancellation. The
   orders satisfy `ord_e G >= idx F`, and the inequality is strict in practice (47/47 exact cases). The gcd reduction
   changes no inequality (Lemma 2.2, the multi-point Lemma K1).
3. **Converse** (Theorem 4.1, proved, exact). Take a star partition: blocks `{o} ∪ B_j` with the `B_j` partitioning
   `[M] \ {o}` into 4-sets. Then `T_M ≅ T_5^{(M-1)/4}` as tori over `Q`. Route 5's slab-box forms are exactly
   intra-cluster `F` of block degree `4T`, with `idx F = ord_e G` and `RHS_F = RHS_G`. Route 5's forced divisibility
   (I1) is literally (TL_M) (Lemma 4.0). So the two classes coincide (Cor 4.2). The `T_5` geometry restricts nothing,
   and the gcd reduction is invisible. The shared-prime gain is route 5's forced-divisibility LP (Remark 3.4; on
   actual clusters the gain is `3.5–64` in `log` units against `W ≈ 21`).
4. **Balanced cap with maximal sharing** (Theorem 5.1, Cor 5.2–5.3, proved, exact). This settles round 24's open
   point 1 for intra-cluster configurations. Suppose the cluster has a simple profile that is `eta`-balanced on all
   `2^{M-1}-1` cuts. Then no intra-cluster inequality excludes it, nor does (TL_M) with any `G`, for any `F` and any
   number of blocks, unless some subtuple used has `mu(y_S) > 4(1-eta)`. That is the one-point cap `15/4`. Yet the
   subtuples of such clusters sit at `mu >= 4 - 2^{3-M} - O(eta)`. The proof is the mean-width theorem on `T_M`
   (`threshold/construct.md` Thm B) through Theorem 3.1. Hence the class, even with an ideal zero estimate, can at
   best give `M = O(log log R)`, and **never the uniform bound** or `(L_5^circ)`.
5. **The zero step** (Prop 6.1, proved; zero lemmas cited). Thue–Roth–Dyson zero lemmas need lacunary degrees
   `d_{i+1} <= omega d_i`. With lacunary degrees the `n`-point inequality is the one-point inequality for the leading
   subtuple up to `O(kappa omega/rho)`, so sharing gives nothing. With comparable degrees no product zero lemma
   applies, and the zero step is a base-locus statement on `T_M`.
6. **Growth regime** (Cor 7.1, proved from refereed route-5 theorems). Here `log r = o(M)`, and the intra-cluster
   method with the auxiliary polynomial chosen after the profile *is* LP-supercritical. Its output is an exact
   relation `F(y_{S_1}, ..., y_{S_q}) = 0` among star subtuples, which is route 5's rigidity. A growth statement needs
   a nonvanishing `F`. That is route 5's zero estimate at the actual point (ZE-col*; MIN^pair / `(R_H)` suffices,
   giving `O(log log R)`), open. Counting in `T_5` cannot supply it, because the base locus is template-dependent.
7. **Tests** (§8). The reduction identities are exact on three actual clusters. Same-profile angle models share every
   input of the method except the value, so they can be separated only by exact relations. The monomial relations
   (50 on the 10-point cluster) are fake-blind. A supercriticality scan at `M = 10` modulo these relations finds
   nothing beyond the critical pair inequality (`max E = -1.5e-4`; evidence, §8).

**Answer to the task's question.** The argument is not genuinely different. Its "new" ingredients (the `T_5` geometry
and the gcd reduction) are invisible to Liouville-type inequalities. Its one real ingredient, forced divisibility
across subtuples sharing primes, is route 5's LP. It beats the `15/4` cap only where route 5 is supercritical, and
there only by producing exact relations whose use is route 5's zero estimate.

## 1. Setting and notation

**Labels.** *Proved* = proved here with all quantifiers, or cited with the place of a refereed (or killer-checked)
proof. *Exact* = finite exact computation in `intra_checks/` (Python integers, `Fraction`, exact Gaussian rationals).
*Evidence* = a computation no proof uses. *Heuristic* = no proof. *Record* = proved earlier in the record. Nothing in
this file is refereed. Paths are relative to `docs/claude_review/` unless they start with `intra_checks/`.

**1.1 Tori.** `K = Q(i)`. For `k >= 2`, `M_k = {s in Z^k : sum s_j = 0}` (characters) and
`T_k(Q) = {x in (K^*)^k : xbar = x^{-1}}/K^*_{N=1}` the norm-one similarity torus (`gauss.md` §1.1; dimension
`k-1`). A character `s in M_k` gives `x^s = prod x_j^{s_j}`, well defined on `T_k` because `sum s = 0`. The identity is
`e`. For `S = (s_1 < ... < s_5) ⊆ [M]`, `iota_S : M_5 -> M_M` puts `s` on the coordinates of `S` (zero elsewhere), and
`pr_S : T_M -> T_5`, `[x] -> [(x_j)_{j in S}]` is the dual homomorphism (`(pr_S x)^s = x^{iota_S s}`). `pr_S` is
defined over `Q` (it commutes with conjugation and inversion). For blocks `S_1, ..., S_n` put
`Phi = (pr_{S_1}, ..., pr_{S_n}) : T_M -> T_5^n` and `Phi^* : M_5^n -> M_M`, `(s_1..s_n) -> sum_i iota_{S_i} s_i`.

**1.2 Clusters.** A *normalised cluster* is `z_1, ..., z_M` distinct in `Z[i]`, `|z_j|^2 = N` (`N` odd, all
`p | N` split, `gcd_{Z[i]}(z) = 1`), with all mutual distances `<= L` (an arc of length `<= C sqrt R`, `R^2 = N`, gives
`L <= C sqrt R`). Write `z_x = eps_x prod_p pi_p^{a_x(p)} pibar_p^{e_p - a_x(p)}` (`two_thirds/proof.md` Lemma 1.1) with
`a_p = (a_x(p))_x in Z^M`, `W = log N`. The **cluster point** is `Z = [z] = (z_j/z_1)_j in T_M(Q)`. For a 5-subset
`S`: `g_S = gcd_{Z[i]}(z_j : j in S)`, `zeta_S = z_S/g_S` (primitive, on `|zeta| = R''_S`), and `y_S = [zeta_S] in T_5(Q)`
— exactly the points of the task (`gauss.md` Lemma 1.2: `[zeta_S] = [z_S]`).

**1.3 Valuations, heights, distances.** For a split place `pi | p`, `v_pi(Z) = a_p mod Z·1` and
`v_pibar(Z) = -v_pi(Z)`; at inert `q` and at `(1+i)` all coordinates of `Z` are units (`round24/theory.md` Lemma 2.1).
For a finite set `Q ⊂ M_k` and `v in Z^k/Z·1`, `width_Q(v) = max_{s in Q}<s,v> - min_{s in Q}<s,v>` (well defined on
`M_k`). `R_5 = conv(e_a - e_b)`; `width_{R_5}(v) = 2(max v - min v)`. The toric height of `y in T_5(Q)` is
`h_B(y) = (1/2) sum_p width_{R_5}(v_pi(y)) log p = log R''^2` (`gauss.md` Lemma 1.1; `round24/theory.md` Prop 2.2(iii)).
`lambda(y) = -log dist(y,e)`, `dist(y,e) = max_{a,b} |y_a/y_b - 1|` (any representative; for `y = y_S` this is
`max_{a,b in S} |z_a - z_b|/R`). For the cluster point, `lambda(Z) >= log(R/L) >= W/4 - log C`.

**1.4 The comparison inequalities.**
* **(TL_k) Toric Liouville on `T_k`** (`round22/toric.md` Thm 3.2 for `k = 5`; the proof is the product formula and is
  repeated in Lemma 3.3 below for every `k`). For `0 != G in Z[i][M_k]` with `G(x) != 0`:
  `log|G(x)| >= -(1/2) sum_p width_{supp G}(v_pi(x)) log p` (exactly, no `O(1)`), against the archimedean Taylor
  bound `log|G(x)| <= -ord_e(G) lambda(x) + C(G)`.
* **(*) The `n`-point inequality** (`round24/theory.md` Prop 4.3): for `0 != F in Z[i][M_5^n]`, `F(y_1..y_n) != 0`,
  `min_{F_J != 0} sum_i j_i lambda(y_i) <= (1/2) sum_p width_{Q_F}((v_pi(y_i))_i) log p + C(F)`, where `F_J` are the
  multihomogeneous Taylor components at `(e,...,e)` and `C(F) << log H(F) + sum_i d_i log(2+d_i)`.
* **Route 5** (`round10/determinant.md` §1, Thm 2.1, Thm 4.1; `round12/README.md` §2.1; refereed): monomials
  `G_k(z) = N^{(D-|k|_1)/2} prod z_x^{k_x^+} zbar_x^{k_x^-}`, `k in M_n` (`n = M`), forms
  `F_c(z) = sum_{k in A} c_k G_k(z)`, the slab box `A = {k : |k_x| <= T (x < n), |<k, beta_p>| <= w}` built from the
  cluster's own profile, the order condition `c in V(A, m)` (`sum_k c_k P(k) = 0` for `deg P < m`), the forced
  divisibility (I1) `Norm F_c(z) >= N^D prod_p p^{-h_{supp c}(a_p)}`, rigidity `F_c(z) = 0` for all `c in V_H`
  (LP-supercritical regime), and the open zero estimate ZE / MIN^pair / (R_H).

**1.5 What "beats the cap" and "reduces" mean.** The one-point balanced cap (`round24/theory.md` Prop 4.3(iii)):
on `eta`-balanced simple points of `T_5`, (TL_5) and the disjoint-prime `n`-point (*) hold whenever
`mu(y) := lambda(y)/w <= 4(1-eta)`; `alpha = 15/4` is `mu = 4`. An intra-cluster argument *beats the cap* if it
produces a contradiction (or an exact relation) for configurations whose subtuples satisfy `mu(y_S) <= 4 - c`,
`c > 0` fixed. It *reduces to route 5* if every inequality it produces is implied by (TL_M) at `Z` for a Laurent
polynomial on `T_M` chosen after the profile — the inequality route 5 uses — and, conversely, route 5's own
inequalities are intra-cluster ones.

## 2. The intra-cluster data are one point of $T_M$ (Lemma 2.1, Lemma 2.2)

**Lemma 2.1 (the subtuple points are projections of the cluster point; proved, exact check).** Let `z_1..z_M` be a
normalised cluster with cluster point `Z`, and `S ⊆ [M]`, `|S| = 5`. Then:
1. `y_S = pr_S(Z)`; for every `s in M_5`, `x^s(y_S) = prod_{j in S} zeta_j^{s_j} = prod_{j in S} z_j^{s_j} = Z^{iota_S s}`.
2. For every split `p`: `v_pi(zeta_j) = a_j(p) - min_{S} a(p)` (`j in S`), so `v_pi(y_S) = v_pi(Z)|_S` modulo `Z·1`, and
   for every finite place `w` and every `(s_1..s_n) in M_5^n`:
   `<Phi^*(s_1..s_n), v_w(Z)> = sum_i <s_i, v_w(y_{S_i})>`.                                             (2.1)
3. `N(g_S) = prod_p p^{min_S a(p) + e_p - max_S a(p)}`, hence the height of the gcd-reduced subtuple is
   `h_B(y_S) = log(N/N(g_S)) = sum_p (max_S a(p) - min_S a(p)) log p = (1/2) sum_p width_{R_5}(v_pi(Z)|_S) log p`,
   and `W/2 <= h_B(y_S) <= W` when `L <= sqrt R` (more generally `h_B(y_S) >= W - 2 log L`).
4. `lambda(y_S) >= lambda(Z)`.

*Proof.* (1) `zeta_j = z_j/g_S` and `sum_j s_j = 0`, so `prod zeta_j^{s_j} = g_S^{-sum s_j} prod z_j^{s_j}`; and
`Z^{iota_S s} = prod_{j in S}(z_j/z_1)^{s_j} = prod_{j in S} z_j^{s_j}`, again because `sum s_j = 0`. (2) `g_S` has
`pi`-valuation `min_S a(p)` (the `pi`-part of a gcd of elements `pi^{a}pibar^{e-a}·unit` is `pi^{min a}`, and
`pibar^{min(e-a)}` for `pibar`). (2.1) is the pairing of (1): `<iota_S s, v> = <s, v|_S>`. (3) As in (2), and
`|g_S|^2 = N(g_S)`; `R''^2_S = N/N(g_S)` (`gauss.md` Lemma 1.2); `max - min` is half of `width_{R_5}`. For the bounds:
`h_B(y_S) <= W` trivially, and `|g_S| <= |z_a - z_b| <= L` for some pair in `S` (`gauss.md` Lemma 1.2), so
`h_B(y_S) = W - log N(g_S) >= W - 2 log L >= W/2` if `L <= sqrt R = e^{W/4}`. (4) `dist(y_S, e) =
max_{a,b in S}|z_a - z_b|/R <= dist(Z, e)`. `□`

*Exact check* (`intra_checks/check1_reduction.py` → `check1_reduction.out`): on the record's 10-point cluster
(`N' = 1176852625`, primes `5^3 13^2 17 29 113`), its 9-point subcluster, and an 8-point window of
`N' = 3159797225`, all 252 / 126 / 56 five-subsets satisfy (1)–(3) exactly, with `h_B(y_S)/W in [0.7737, 1]`,
`[0.7737, 1]`, `[0.8349, 1]`.

*Reading.* The `binom(M,5)` points `y_S` are **not** independent data: they are the images of the single rational
point `Z` of the `(M-1)`-dimensional torus `T_M` under `binom(M,5)` homomorphisms. Their heights `H_S`, their gcds
`g_S`, their valuation vectors and the primes they share are all functions of the one profile `(a_p)_p` of `Z`; their
distances to `e` are bounded below by `lambda(Z)`. "Shared primes" means precisely that the `v_pi(y_S)` are
restrictions of one vector `v_pi(Z)`.

**Lemma 2.2 (the gcd reduction is invisible to forced divisibility; proved).** Let `F` be a polynomial in the
coordinates of the reduced tuples `zeta_{S_1}, ..., zeta_{S_n}` and their conjugates, with `Z[i]` coefficients,
multihomogeneous of bidegree `(a_i, b_i)` in `(zeta_{S_i}, zetabar_{S_i})`. Then
`F(zeta_{S_1}, ...) = prod_i g_{S_i}^{-a_i} gbar_{S_i}^{-b_i} F(z_{S_1}, ...)`, and at every finite place `w` of `K` the
bound `|F(zeta)|_w <= 1` (integrality of the reduced tuples) is implied by the toric bound
`|F(zeta)|_w <= max_{m in supp F} |m(zeta)|_w`. On the circle, `F(z_{S_1},...)` is a monomial in the `z_j, zbar_j`
times a Laurent polynomial `G` on `T_M` evaluated at `Z` (`zbar_j = N/z_j`), with the same Newton data; so every
inequality obtained from integrality or forced divisibility of the reduced tuples is implied by (TL_M) for `G` at `Z`.

*Proof.* Each monomial `m` of bidegree `(a_i,b_i)` satisfies `m(z) = prod_i g_{S_i}^{a_i} gbar_{S_i}^{b_i} m(zeta)`
(Lemma K1, `kill_gauss.md` §6.1). `m(zeta)` is a Gaussian integer, so `|m(zeta)|_w <= 1`, and the ultrametric
inequality gives `|F(zeta)|_w <= max_m |m(zeta)|_w <= 1`. Writing `zbar_j = N/z_j` turns `F(z)` into
`N^{b} prod z^{(\cdot)}` times a sum of characters `Z^k`; the factor in front is the same for all monomials (it is
`z_1^{sum(a_i - b_i)} N^{sum b_i}` after dividing by `z_1`), so the toric bound for `G` at `w` is the toric bound for
`F(z)` divided by the common factor, which is exactly the bound (I1) of route 5. `□`

So "using the gcd reduction" adds nothing to a Liouville-type argument beyond (TL_M): the reduction changes the
**heights** `h_B(y_S)` (which are toric heights of projections, Lemma 2.1(3)), not the inequalities. This is
Lemma K1 in multi-point form. The only objects in the record that see the reduction are Vojta's `D = E` term and the
heuristic count (`kill_gauss.md` §6.1); neither is an auxiliary-polynomial argument.

## 3. The intra-cluster multi-point inequality and the reduction (Theorem 3.1)

**Definition 3.0 (the intra-cluster multi-point method `L_intra`).** Data: a normalised cluster, blocks
`S_1..S_n` (5-subsets of `[M]`, repetitions and overlaps allowed), and `0 != F in Z[i][M_5^n]` with Newton polytope
`Q_F` — chosen in any way, in particular **after** the profile `(a_p)` and the positions are known — with
`F(y_{S_1}, ..., y_{S_n}) != 0`. Archimedean side: any upper bound `log|F(y)| <= U`; the Dyson-type bound is
`U_F = log sum_{J in Min} e^{-sum_i j_i lambda(y_{S_i})} + C(F)` (`Min` = minimal multidegrees of the Taylor support of
`F` at `(e,..,e)`), weights being implicit in `lambda(y_{S_i})`. Non-archimedean side: toric forced divisibility at
every place of `K`, using the shared structure, i.e. the full vector `(v_pi(y_{S_i}))_i` (Lemma 2.1(2)), plus anything
implied by integrality of the gcd-reduced tuples (Lemma 2.2). Conclusion: if `U < -(1/2) sum_p
width_{Q_F}((v_pi(y_{S_i}))_i) log p` then `F(y) = 0`.

**Lemma 3.2 (orders; proved, exact check).** Let `G = Phi^*F := F∘Phi in Z[i][M_M]`. Then
(a) `G(Z) = F(y_{S_1}, ..., y_{S_n})` and `supp G ⊆ Phi^*(supp F)` (equality unless monomials cancel);
(b) `idx(F) := min{|J| : F_J != 0} <= ord_e(G) <= reg(supp G)` whenever `G != 0`;
(c) if `Phi` is an isomorphism of tori (Theorem 4.1), `idx(F) = ord_e(G)`.
*Proof.* (a) Lemma 2.1(1) monomial by monomial. (b) In logarithmic coordinates `x = exp(t)`, `t in Lie(T_M)`, the map
`Phi` is the linear map `dPhi`. So `G(exp t) = sum_J F_J(dPhi t)` and the degree-`m` Taylor part of `G` is
`sum_{|J| = m} F_J∘dPhi`, which vanishes for `m < idx(F)`. `ord_e G <= reg(supp G)` is the definition of `reg`
(`round22/toric.md` Lemma 3.1). (c) If `dPhi` is invertible, `sum_{|J|=m}F_J∘dPhi = 0` iff `sum_{|J|=m}F_J = 0` iff
`F_J = 0` for all `|J| = m` (the `F_J` with `|J| = m` are the multihomogeneous parts of one homogeneous polynomial). `□`

*Exact check* (`check2_cap.py` part (ii)): 94 auxiliary `F` (blocks in `[6]`, `[7]`, `[9]`), built by exact kernels
either with `T_5^n`-index `>= t` (mode 0) or with `T_M`-order `>= t` (mode 1). In every case
`idx <= ord <= reg`; with star blocks `idx = ord`; and in **all** 47 mode-1 cases `idx(F) < ord_e(G)`
(e.g. `idx = 1`, `ord = 3`). So the Dyson index bookkeeping on `T_5^n` is in general strictly weaker than the
order on `T_M`.

**Lemma 3.3 (exact product formula on `T_k`; proved, exact check).** Let `x in T_k(Q)` and `0 != G in Z[i][M_k]`
with `G(x) != 0`. Then `Norm(G(x)) · prod_p p^{width_{supp G}(v_pi(x))}` is a positive integer; in particular
`log|G(x)| >= -(1/2) sum_p width_{supp G}(v_pi(x)) log p`.
*Proof.* At `pi`: `v_pi(G(x)) >= min_{k in supp G}<k, v_pi(x)>` (coefficients in `Z[i]`, ultrametric). At `pibar`:
`v_pibar(G(x)) >= min_k <k, -v_pi(x)> = -max_k <k, v_pi(x)>`. At inert `q` and at `1+i`, `v(G(x)) >= 0` (Lemma 2.1 of
`round24/theory.md`). Hence `G(x) · prod_p pi^{-min}pibar^{max}` is a Gaussian integer, nonzero; its norm is
`Norm(G(x)) prod_p p^{max - min}`. `□` (Exact check, `check1_reduction.out` (c): every random `F` with `G(Z) != 0`
among 120 drawn, integral `>= 1`.)

**Theorem 3.1 (reduction; proved).** Every inference of `L_intra` is implied by the one-point inequality (TL_M) at the
cluster point `Z` for the single Laurent polynomial `G = Phi^*F`, with the same archimedean bound. Precisely, with
`RHS_F = (1/2) sum_p width_{Q_F}((v_pi(y_{S_i}))_i) log p` and `RHS_G = (1/2) sum_p width_{supp G}(v_pi(Z)) log p`:
1. `log|F(y)| = log|G(Z)|`, so any `U` bounding one bounds the other;
2. `RHS_F = (1/2) sum_p width_{Phi^*Q_F}(v_pi(Z)) log p >= RHS_G` (equality iff no cancellation changes a width);
3. hence if `L_intra` concludes `F(y) = 0` (because `U < -RHS_F`), then (TL_M) for `G` concludes `G(Z) = 0`, which is the
   same statement; and (TL_M) may conclude it in cases where `L_intra` does not (`RHS_G < RHS_F`, or a better order
   `ord_e G > idx F` in the Taylor bound).
4. The same holds with the gcd-reduced tuples and integrality in place of toric divisibility (Lemma 2.2), with
   `iota`-conjugates or group words of the `y_S` among the points (they are again images of `Z` under homomorphisms
   `T_M -> T_5`; `round24/referee_theory.md` Lemma R1 is the case of words), and for split points.
*Proof.* (1) Lemma 3.2(a). (2) By (2.1), `<s, (v_pi(y_{S_i}))_i> = <Phi^*s, v_pi(Z)>` for every `s in supp F`, so the two
widths are max-minus-min of the same numbers; `supp G ⊆ Phi^*(supp F)` gives `>=`. (3) Combine (1), (2) with
Lemma 3.3. (4) Lemma 2.2; an `iota`-conjugate is `pr_S∘iota`, a word `y_{S_1}^{±1}y_{S_2}^{±1}...` is the image of `Z`
under the homomorphism `prod_i pr_{S_i}^{±1}`, and the argument above uses only that `Phi` is a homomorphism of tori
defined over `Q`. `□`

*Exact check* (`check1_reduction.out` (c)): on the three actual clusters, 120 random `F` with `n in {2,3,4}` blocks,
block degrees `1..3`, and forced cancellations: `F(y) = G(Z)` exactly in `Q(i)`; `width_{Q_F} = width_{Phi^*Q_F}`
always, `> width_{supp G}` at 109 prime-instances (cancellation).

**Remark 3.4 (where the shared primes enter).** For disjoint prime sets, `RHS_F` splits as the sum of one-point right
sides (`round24/theory.md` Prop 4.3(ii)). For shared primes it does not, and it can be much smaller: on the three
actual clusters, `sum_i RHS_{pr_i Q_F}(y_{S_i}) - RHS_F` ranged over `[3.5, 64.3]` (medians `21–34`, in units of
`log`; `check1_reduction.out` (d)), against `W ≈ 21`. This **is** the gain the task asks about, and Theorem 3.1 says
where it lives: it is the difference between the one-point widths of the projections and the width of one polytope
`Phi^*Q_F` at one vector `v_pi(Z)` — the forced-divisibility LP of route 5 (`round10/determinant.md` §1, (1.1)).

## 4. Route 5 is an intra-cluster multi-point argument (Theorem 4.1: star partition)

**Lemma 4.0 (route 5's inequality is (TL_M); proved from the record).** In route 5's notation (`n = M` points),
`G_k(z) = N^{D/2} Z^k` for `k in M_M`, `|k|_1 <= D` (`round10/determinant.md` §1: `G_k(z) = N^{D/2}e^{i<k,delta>}` and
`e^{i<k,delta>} = Z^k`). Hence `F_c(z) = N^{D/2} G_c(Z)` with `G_c = sum_{k in A} c_k x^k in Z[M_M]`, and (I1)
`Norm F_c(z) >= N^D prod_p p^{-h_{supp c}(a_p)}` is literally Lemma 3.3 for `G_c` at `Z`
(`h_{supp c}(a_p) = width_{supp c}(v_pi(Z))`). The order condition `c in V(A,m)` is `ord_e(G_c) >= m` (the Taylor
coefficients of `G_c(exp(t))` are `sum_k c_k <k,t>^j/j!`), and the archimedean bound of `round10/determinant.md` (4.1)
is the Taylor bound at `Z` with the isotropic scale `Theta = max|delta|`. `□`

**Theorem 4.1 (star partition: `T_M ≅ T_5^q`, and route 5 is an intra-cluster argument; proved, exact check).**
Let `M = 4q + 1`, fix an anchor `o in [M]`, and partition `[M] \ {o}` into 4-sets `B_1..B_q`. Put `S_j = B_j ∪ {o}`.
1. `Phi = (pr_{S_1}, ..., pr_{S_q}) : T_M -> T_5^q` is an isomorphism of tori over `Q`; `Phi^*` is the identity matrix
   in the bases `(e_x - e_o)_{x != o}` of `M_M` and `(e_x - e_o)_{x in B_j}` of the `j`-th `M_5`.
2. For every route-5 support at box size `T` (`|k_x| <= T` for `x != o`), the unique preimage `(s_1..s_q)` of `k` has
   `s_j = sum_{x in B_j} k_x (e_x - e_o) in 4T·R_5` (`sum|s_j| <= 8T`). So `F_c(z) = N^{D/2} F(y_{S_1}, ..., y_{S_q})` with
   `F = (Phi^*)^{-1}G_c in ⊗_j Z[4T R_5 ∩ M_5]`, no cancellation, `RHS_F = RHS_{G_c}`, and `idx(F) = ord_e(G_c)`
   (Lemma 3.2(c)).
3. Consequently route 5's auxiliary-polynomial inferences — LP supercriticality (`round10/determinant.md` Thm 2.1)
   and rigidity `F_c(z) = 0` for all `c in V_H` (Thm 4.1 there; `round12` at `T = 1`) — are inferences of `L_intra`
   about the `q = (M-1)/4` star subtuples `y_{S_1}, ..., y_{S_q}` of one cluster, which all contain the anchor `o`
   (so every prime that separates `o` from a point of `B_j` is shared by `y_{S_j}` with every other block it touches).
   Route 5's remaining steps (Step 0, the separation lemma, MIN^pair) concern the positions `u` and the base locus,
   not the auxiliary inequality, and transfer unchanged.
4. For general `M`, take `q = ⌈(M-1)/4⌉` blocks covering `[M] \ {o}` (the last block padded with indices that carry
   exponent `0`); then `Phi^*` is still injective on route-5 supports and (2)–(3) hold verbatim (`Phi` is then a
   surjection onto a subtorus, not an isomorphism; `idx(F) <= ord_e(G_c)`).
*Proof.* (1) `x^{e_x - e_o}(Phi(Z)) = Z^{e_x - e_o}` for `x in B_j` (Lemma 2.1(1)), so `Phi^*` maps the basis to the basis;
an isomorphism of character lattices compatible with the Galois action (`Phi` is defined over `Q`) is an isomorphism
of tori. (2) `k = sum_{x != o} k_x(e_x - e_o)` (since `sum k = 0`), grouped by blocks; `sum|s_j| = sum_{x in B_j}|k_x| +
|sum_{x in B_j} k_x| <= 8T`, and `sum|s| <= 2m` is membership in `m R_5`. The remaining claims follow from
Lemma 4.0, Theorem 3.1(2) (no cancellation since `Phi^*` is injective) and Lemma 3.2(c). (4) as (2). `□`

*Exact check* (`check2_cap.out` (i)): `M = 9, 13, 17`: `Phi^*` is the identity; 2000 random slab exponents at
`T = 1, 2` each pull back with `s_j in 4T R_5`.

**Corollary 4.2 (the two methods coincide; proved).** As classes of inferences,
`L_intra ⊆ {(TL_M) at Z for some G in Z[i][M_M] chosen after the profile} = route 5's framework ⊆ L_intra`:
the first inclusion is Theorem 3.1, the last is Theorem 4.1 (any `G` is `Phi^*` of `F = (Phi^*)^{-1}G` for the star
partition). Route 5 is the special choice "`G` from a Siegel lemma on the profile's slab box, order `m' >= 4w`";
`L_intra` with general blocks is the same `G` with a weaker order bookkeeping (`idx(F) <= ord_e G`, strict in all 47
mode-1 checks) and a weaker width (`RHS_F >= RHS_G`). **The intra-cluster multi-point argument is the determinant
method in disguise.**

*What the `T_5` framing adds, item by item.*
* *The `T_5` geometry.* Only through the block structure of `Q_F`: `Q_F ⊆ prod_i d_i R_5`. This restricts nothing
  (Theorem 4.1(2): every slab monomial is in `⊗ 4T R_5`), and toric Liouville does not see the compactification
  (`Lbar_5^tw` versus `Lbar_M^tw`): both sides of (TL) depend on `supp G` and the point only.
* *The gcd reduction.* Invisible (Lemma 2.2): it fixes the heights `h_B(y_S)`, which are functions of `v_pi(Z)`, and
  changes no inequality. Route 5 did use it implicitly — (I1) is the same bound.
* *The shared prime structure.* It is the width of `Phi^*Q_F` at `v_pi(Z)` — route 5's LP (Remark 3.4).
* *Per-subtuple distances `lambda(y_S) >= lambda(Z)`.* A refinement of the archimedean bound, already contained in
  the exact Taylor expansion of `G` at `Z` (route 5's `Phi_c(u) = sum c_k <k',u>^{m'}` is the leading term in the
  direction `u` of the positions, which encodes all subtuple diameters). It does not touch the non-archimedean side,
  where the gap is.

## 5. Balanced cap with shared primes (Theorem 5.1)

Round 24 left open whether the `n`-point inequality (*) has a cap when the points share primes
(`round24/README.md` §1.1 item 5; `theory.md` Prop 4.3, Remark (1)). For the configurations that matter here — all
points are subtuples of one cluster — Theorem 3.1 turns the question into a one-point question on `T_M`, where the
mean-width theorem applies.

**Theorem B on the zero-sum slice (record).** `threshold/construct.md` Theorem B (refereed): for finite `S ⊂ Z^M`,
`reg(S) <= 2^{1-M} sum_{P ⊆ [M]} w_P(S)`, `w_P(S) = max_S x(P) - min_S x(P)`. For `S ⊂ M_M`, `x(P) = -x(P^c)`, so
`w_P = w_{P^c}`, `w_∅ = w_{[M]} = 0`, and `w_A(S) = width_S(1_A)`; hence
`ord_e(G) <= reg(supp G) <= 2^{2-M} sum_{cuts {A,A^c}} width_{supp G}(1_A)`.        (5.1)
(`reg` on `Z^M` at `(1,..,1)` equals `reg` on `T_M` at `e`, because a Laurent polynomial supported on `M_M` is the
pull-back of a function on `T_M` along the smooth quotient `G_m^M -> T_M`.) For `M = 5` this is the `2^{-3}` of
`round24/theory.md` §1.3 (MW). *Exact check* (`check2_cap.out` (ii)): (5.1) on all 94 sampled supports, `M = 6, 7, 9`.

**Theorem 5.1 (balanced intra-cluster cap; proved, exact check).** Let `M >= 5`, `0 <= eta < 1`, `w > 0`. Let
`Z in T_M(Q)` be the point of a normalised cluster whose profile is **simple** (`e_p = 1` for all `p | N`, so
`v_pi(Z) = 1_{A_p}` for a cut `A_p`) and **lower `eta`-balanced on all cuts**: `log n_A >= (1-eta)w` for each of the
`2^{M-1} - 1` cuts, `n_A = prod_{p : A_p ∈ {A, A^c}} p`. Let `S_1..S_n`, `F` (any, chosen after the profile) with
`F(y_{S_1}, ..., y_{S_n}) != 0`, and `G = Phi^*F`. Then
`LHS(*) := min_{F_J != 0} sum_i j_i lambda(y_{S_i}) <= idx(F) · max_i lambda(y_{S_i})`, and
`RHS_F >= RHS_G >= 2^{M-3}(1-eta) w · ord_e(G) >= 2^{M-3}(1-eta) w · idx(F)`.
Consequently:
1. (*) holds — `L_intra` derives no contradiction — whenever
   `max_i lambda(y_{S_i}) <= 2^{M-3}(1-eta)w - |C(F)|/max(1, idx F)`;
2. the one-point Taylor form on `T_M` (route 5's form, Lemma 4.0) derives no contradiction whenever
   `lambda(Z) <= 2^{M-3}(1-eta)w - |C(G)|/max(1, ord_e G)`, for **every** `G in Z[i][M_M]`, chosen after the profile or not.
*Proof.* `LHS(*) <= sum_i j^0_i lambda_i <= |J^0| max_i lambda_i` for `J^0` with `|J^0| = idx(F)`. `RHS_F >= RHS_G` is
Theorem 3.1(2). Simplicity and balance: `RHS_G = (1/2) sum_p width_{supp G}(1_{A_p}) log p = (1/2) sum_{cuts}
width_{supp G}(1_A) log n_A >= (1/2)(1-eta) w sum_{cuts} width_{supp G}(1_A) >= 2^{M-3}(1-eta) w ord_e(G)` by (5.1).
`ord_e G >= idx F` is Lemma 3.2(b) (`G != 0` because `G(Z) = F(y) != 0`). (1): a contradiction needs
`LHS(*) - C(F) > RHS_F`, i.e. `idx(F)(max lambda_i - 2^{M-3}(1-eta)w) > -|C(F)|`, excluded by the hypothesis
(if `idx F = 0` the Taylor bound is `>= log|F(y)|`'s trivial bound and nothing is derived). (2) the same with
`ord_e G`, `lambda(Z)`. `□`

*Exact check* (`check2_cap.out` (iii)): cube profile (every cut one prime, `w = 1`), 94 kernel-built `F` with
shared blocks in `[6]`, `[7]`, `[9]`: `RHS(*) >= 2^{M-3} ord_e(G) >= 2^{M-3} idx(F)` in every case.

**Corollary 5.2 (the window of the subtuples; proved).** In the setting of Theorem 5.1, suppose moreover the profile
is upper balanced (`log n_A <= (1+eta)w`), the cluster has diameter exactly `L = C sqrt R` (so
`lambda(Z) = W/4 - log C`), and every subtuple used has diameter `>= L/kappa`. Then:
1. every 5-subset `S` has scale `w_S`: each of its 15 cuts receives exactly `2^{M-5}` cuts of `[M]`, so
   `log n_{A'}(y_S) in [2^{M-5}(1-eta)w, 2^{M-5}(1+eta)w]`, and
   `mu(y_S) = lambda(y_S)/(2^{M-5}w) >= (4 - 2^{3-M})(1-eta) - (log C)/(2^{M-5}w)`;
2. the hypothesis of Theorem 5.1(1) holds as soon as `eta < 1/(2^M - 1)` and
   `w >= 4(log(kappa/C) + |C(F)|/max(1,idx F))/(1 - (2^M-1)eta)`.
So, measured in the subtuples' own units, `L_intra` excludes a balanced cluster only if some subtuple has
`mu(y_{S_i}) > 4(1-eta)`: **the intra-cluster cap is the one-point cap `15/4`, even with maximal prime sharing**,
while the subtuples sit at `mu >= 4 - 2^{3-M} - O(eta) - o(1)`, inside the window that `(L_5^circ)` must close.
*Proof.* (1) A cut `{A', S\A'}` of `S` is induced by the cuts `{A, A^c}` of `[M]` with `A ∩ S ∈ {A', S\A'}`; there are
`2^{M-5}` of them (choose `A \ S` freely, then identify `A` with `A^c`). The other `2^{M-5} - 1` cuts restrict
trivially and form `g_S` (Lemma 2.1(3)). `lambda(y_S) >= lambda(Z) = W/4 - log C` (Lemma 2.1(4)) and
`W >= (2^{M-1}-1)(1-eta)w`. (2) `lambda(y_{S_i}) <= log(R/(L/kappa)) = lambda(Z) + log kappa` and
`W <= (2^{M-1}-1)(1+eta)w`, so `max_i lambda_i <= (2^{M-1}-1)(1+eta)w/4 - log C + log kappa`; and
`2^{M-3}(1-eta) - (2^{M-1}-1)(1+eta)/4 = (1 - (2^M-1)eta)/4`. `□`

**Corollary 5.3 (no uniform bound from `L_intra`, nor from route 5's framework; proved as a statement about the
method).** Fix `M`. Every inference of `L_intra` (or of (TL_M) with any `G`, any order, any choice after the profile,
any zero estimate) is consistent with clusters whose profile is simple and `eta`-balanced on all cuts with
`eta < 1/(2^M-1)` and `w` as in Cor 5.2(2). Such profiles need `r >= 2^{M-1} - 1` primes, i.e. `W >= 2^{M-1}(1-eta)w`;
with `w ≍ M` (enough primes `≡ 1 (4)` in a window of relative width `eta`) this is `M ≈ log_2 log R - log_2 log log R`.
Hence the class can at best prove `M = O(log log R)` — which is route 5's conditional target (`round12/README.md`
§0 item 3) — and **never the uniform bound**, and never `(L_5^circ)`-strength exclusion of balanced subtuples.
*Proof.* Theorem 5.1 and Cor 5.2: the inequality the method needs to violate holds on these data, for every `G`.
A zero estimate only supplies some `G` with `G(Z) != 0`; it cannot make the inequality fail. `□`

*Scope and reading.*
* Theorem 5.1 is `threshold/construct.md` Cor C3 (the cube profile defeats every auxiliary polynomial) transported to
  the intra-cluster `n`-point setting by Theorem 3.1. The new content is that **prime sharing among subtuples cannot
  undo it**: the sharing gain of Remark 3.4 is bounded by the mean width on `T_M`, and the mean width is exactly what
  a cube profile pays.
* It is a statement about the method on data with the stated parameters. It does not assert that actual clusters with
  cube profiles exist on short arcs (that is what the conjecture forbids for large `M`), nor that cube-profile fake
  circles satisfying every monomial gap exist (not claimed here).
* The constants `C(F)`, `C(G)` are template-dependent; the cap is stated with them explicitly, as in
  `round24/theory.md` Thm 7.1(2). For families with `log H + sum_i d_i log(2+d_i) = o(d·w)` they are lower order.
* Congruence divisibility at primes `q ∤ N` (residues, the pigeonhole of Theorem 1.1) is **not** in `L_intra`; it is a
  separate input class with its own record status (residue designs; P3-all* uncapped). Theorem 5.1 says nothing
  about it.
* In the growth regime the profile is not a cube: `r <= W/log 5 = O(M log M)`, so `log r = O(log M) = o(M)`. There the
  bound (5.1) is far from tight for profile-adapted supports and Theorem 5.1 does not apply — §7.

## 6. The Thue–Roth–Dyson zero step: lacunary degrees kill the sharing gain (Prop 6.1)

A Thue/Siegel/Dyson argument has two halves: the Liouville comparison (handled in §§3–5) and a **zero step**
guaranteeing `F(y) != 0` (or small index at `y`). The classical zero steps — Roth's lemma, Faltings' product theorem,
Dyson's lemma — need **lacunary degrees** `d_{i+1} <= omega d_i` with `omega` small in terms of `n` and the target
index (Bombieri–Gubler, *Heights in Diophantine Geometry*, Ch. 6; Schmidt, LNM 785, Ch. V; *cited from memory: only
the shape "degree ratios `<= omega`, `omega -> 0`" is used, no constant*). With comparable degrees no such lemma is
available, and the zero step is the zero estimate on `Phi(T_M)`, which is route 5's (§7). The next proposition shows
the other half of the dichotomy: with lacunary degrees, the shared primes give (almost) nothing.

**Proposition 6.1 (lacunary degrees: the `n`-point inequality is the one-point inequality for the leading point up to
`O(omega)`; proved, exact check of (a)).** Let `y_1..y_n in T_5(Q) \ {e}` (sharing primes or not), and
`0 != F in Z[i][M_5^n]` with `supp F ⊆ prod_i d_i R_5`, `F(y) != 0`. Put `P_i = pr_i(supp F)`, `R'_i = reg(P_i)`,
`RHS_P(y) = (1/2)sum_p width_P(v_pi(y)) log p`.
(a) `RHS_F >= RHS_{P_1}(y_1) - sum_{i>=2} d_i h_B(y_i)`.
(b) `LHS(*) <= sum_i R'_i lambda(y_i)`, `R'_i <= ⌊10 d_i/3⌋`, and `lambda(y) <= h_B(y)/2` for every `y != e`.
(c) Hence a violation of (*) implies `R'_1 lambda(y_1) > RHS_{P_1}(y_1) - E`, with
    `E = (8/3) sum_{i>=2} d_i h_B(y_i) + |C(F)|`.
(d) If `y_1` is `eta`-balanced simple at scale `w_1`, `d_{i} <= omega^{i-1} d_1` and `h_B(y_i) <= kappa h_B(y_1)`,
    a violation implies `R'_1 (mu(y_1) - 4(1-eta)) w_1 > -E`, `E <= 40(1+eta)kappa omega/(1-omega) · d_1 w_1 + |C(F)|`.
    So if the construction gives the leading block regularity `R'_1 >= rho d_1` (as Thue–Roth constructions do, the
    index being spread with `j_1 ≍ d_1`), `L_intra` with lacunary degrees excludes only points with
    `mu(y_1) > 4(1-eta) - 40(1+eta)kappa omega/(rho(1-omega)) - |C(F)|/(rho d_1 w_1)`: the one-point cap up to `O(omega)`.
    For subtuples of one cluster `kappa <= 2` (Lemma 2.1(3)).
*Proof.* (a) Pick `s, s' in supp F` with `<s_1 - s'_1, v_1> = width_{P_1}(v_1)` (`v_i = v_pi(y_i)`). Then
`width_{Q_F}(v) >= <s - s', v> = width_{P_1}(v_1) + sum_{i>=2}<s_i - s'_i, v_i> >= width_{P_1}(v_1) - sum_{i>=2} d_i
width_{R_5}(v_i)`, since `s_i, s'_i in d_i R_5`. Multiply by `(1/2)log p`, sum over `p`, use
`(1/2)sum_p width_{R_5}(v_pi(y_i)) log p = h_B(y_i)`. (b) The first inequality is `round24/theory.md` Prop 4.2 (some
`J` with `F_J != 0` has `j_i <= R'_i`); `R'_i <= t_max(d_i) = ⌊10d_i/3⌋` is `round23/width.md` Thm 2.5. For
`y != e` pick `a, b` with `y_a != y_b`; Prop 2.2(iv) of `round24/theory.md` gives `|y_a/y_b - 1| >= e^{-d_ab/2}`
(the factor `2^{1-(e_a+e_b)/2} >= 1`), and `d_ab <= h_B(y)`, so `lambda(y) <= h_B(y)/2`. (c) A violation is
`LHS(*) > RHS_F + C(F)`; insert (a), (b) and `R'_i lambda_i <= (10/3)d_i h_B(y_i)/2`, so
`sum_{i>=2}(d_i h_B + R'_i lambda_i) <= (8/3)sum_{i>=2} d_i h_B(y_i)`. (d) For balanced simple `y_1`,
`RHS_{P_1}(y_1) >= (1/2)(1-eta)w_1 sum_{cuts} width_{P_1}(1_A) >= 4(1-eta)w_1 R'_1` by (5.1) at `M = 5`; and
`sum_{i>=2} d_i h_B(y_i) <= kappa h_B(y_1) d_1 omega/(1-omega)` with `h_B(y_1) <= 15(1+eta)w_1`. `□`

*Exact check* (`check2_cap.out` (iv)): inequality (a) at the level of each prime, 300 random `F` with degrees
`(8,2,1)` on blocks of the record's 10-point cluster: holds, minimal slack `1`.

**The dichotomy (proved halves, cited zero lemmas).** To exploit shared primes the auxiliary polynomial must couple
blocks of comparable degree (otherwise Prop 6.1 applies and the gain is `O(omega)`); to have a Thue–Roth–Dyson zero
step the degrees must be lacunary. With comparable degrees the only zero step is a statement about the base locus of a
linear system on `T_M` at the actual point `Z` — route 5's ZE. Note also that for `n > (M-1)/4` the tuple
`(y_{S_i})_i` lies on the proper subvariety `Phi(T_M)` of `T_5^n` (dimension `M-1 < 4n`), so many nonzero `F` vanish
at it identically (e.g. `x^{e_a-e_b}(y_S) x^{e_b-e_c}(y_{S'}) - x^{e_a-e_c}(y_{S''})` for `a,b in S`, `b,c in S'`,
`a,c in S''`); any zero step must see `Phi(T_M)`, i.e. work on `T_M`.

## 7. Growth regime: supercritical, exact relations, zero estimate (Cor 7.1)

In the growth regime (`M ≍ log R/log log R`) the profile is sparse: `r <= W/log 5`, so `log r = O(log M) = o(M)`, far
from the cube profile of §5. There the intra-cluster method **does** beat Liouville, but only in route 5's way.

**Corollary 7.1 (intra-cluster rigidity; proved from the refereed route-5 theorems).** Let a normalised cluster have
`n = M >= 64`, `n - 1 >= 300 ln 4r` and `W >= 16 ln n + 8 log^+ C + 14` (`round10/determinant.md` Thm 4.1 as stated in
`overview/overview.md` §3.4; at `T = 1` the hypotheses of `round12/README.md` §2.1: `d >= 316 ln 4r`, `W > W_1`). Let
`S_1..S_q` be a star partition (Theorem 4.1, `q = ⌈(M-1)/4⌉`). Then for every `F` in the image
`(Phi^*)^{-1}(V_H)` of route 5's integral auxiliary space (block degrees `<= 4T`):
`F(y_{S_1}, ..., y_{S_q}) = 0`.
In words: the intra-cluster multi-point inequality, with the auxiliary polynomial chosen after the profile, is
**LP-supercritical** (by the factor `>= 0.12 sqrt((n-1)/ln 4r)` of `round10/determinant.md` Thm 2.1), and its output is
not a contradiction but an **exact algebraic relation** among the `q` star subtuples. *Proof.* Theorem 4.1(2)–(3):
`F_c(z) = N^{D/2} F(y_{S_1}, ...)` and route 5's rigidity gives `F_c(z) = 0`. `□`

**What is still needed (open; identical to route 5's).** A growth statement needs one `F` in that space with
`F(y_{S_1}, ..., y_{S_q}) != 0` for every regime cluster — by Theorem 4.1 this is literally route 5's zero estimate
at the actual point (ZE-col* in `round11/subcluster.md`; the box form MIN^pair / `(R_H)` suffices,
`round12/README.md` §2.3, `round16/README.md` §5.1). Then `M = O(log log R)` on all norms (`round12/nearstrata.md`
Thm 4.1, refereed implication), and by Cor 5.3 nothing better can come from this class. Neither the multi-point
framing nor the `T_5` geometry supplies the zero estimate:
* *Product zero lemmas* need lacunary degrees, which kill the sharing gain (§6).
* *Counting in `T_5`* (`gauss.md` Lemma 3.3: a fixed proper closed `Z_5 ⊂ T_5` contains `<= c(Z_5)M^4` of the 5-tuples)
  would need the projections `pr_S(B)` of the base locus `B = B(V_H) ⊂ T_M` to lie in a **fixed** `Z_5`. But `B`
  depends on the profile, on `H` and on `m'` (it is cut out by the cluster's own slab box), and its degree grows with
  `m'`; using `c(pr_S B)` as if it were uniform would treat a template-dependent threshold as uniform in the radius,
  which the owner's constraints forbid. With honest constants `c(pr_S B) >= deg(pr_S B)`, which is unbounded.
* *Fakes.* By `round10/determinant.md` Prop 4.2 and Cor 4.3 almost every fake (angle model) violates the exact system, so the
  intra-cluster relations are fake-violated in the same sense as route 5's; whether gap-admissible fakes on the real
  tangent cone exist at barrier scale is the undecided CIRCLE(q) question (`round16/README.md` §0,
  `overview/overview.md` §6.2 item 6). In the low-height window, scenario P (`round21/fakes.md`) would supply fakes for
  the monomial-gap class `A_0`, not for these relations.

So, for growth, the intra-cluster multi-point argument stops exactly where route 5 stops: at the zero estimate /
`(R_H)`.

## 8. Tests on actual clusters and fake circles

All scripts run under `prlimit --as=3500000000`, one process at a time; outputs next to the scripts.

| script → output | what | label | result |
| --- | --- | --- | --- |
| `check1_reduction.py` → `.out` | Lemma 2.1 (1)–(3) on all 5-subsets; Theorem 3.1 (values, widths); Lemma 3.3 integrality; sharing gain | exact | record 10-pt cluster, its 9-pt subcluster, an 8-pt window of `N' = 3159797225`: all identities hold; `h_B(y_S)/W in [0.77, 1]`; `Norm(G(Z)) prod p^{width}` integer `>= 1` in every case; sharing gain `sum_i RHS_i - RHS_F in [3.5, 64.3]` (Remark 3.4) |
| `check2_cap.py` → `.out` | Theorem 4.1 (star partition); Lemma 3.2 `idx <= ord <= reg`; (5.1); Theorem 5.1 on cube profiles; Prop 6.1(a) | exact | all hold; `idx < ord` in 47/47 `T_M`-built `F`, `idx = ord` for star blocks |
| `check3_fakes.py` → `.out` | actual 10-pt cluster versus same-profile angle models; monomial relations; slab scan | exact (A) / evidence (B) | see below |
| `check3b_quotient.py` → `.out` | supercriticality scan modulo monomial relations (190 boxes on `T'`) | evidence | `max E = -1.5e-4 < 0` (pair-critical) |

**Actual clusters (subtuple data).** On the record 10-point cluster (`C ≈ 8.24`), the subtuples have approximation
exponents `h_B(y_S)/lambda(y_S)` between `5.07` and `6.71` (evidence; the full cluster has `W/lambda(Z) = 6.71`): small
record clusters are far from the critical window (`15/4`) — they have large `C` and few primes, and say nothing about
the asymptotic question, as the record has noted for every computation of this size (`round18/census.md` Thm E).

**Fake circles.** The relevant fakes here are angle models with the genuine primes and exponent profile
(`overview/overview.md` §5.1). Because every input of `L_intra` other than the value `F(y)` is a function of the profile
(Lemma 2.1, Theorem 3.1(2)), the method separates an actual cluster from a fake with the same profile and positions
only through **exact relations** — and only in a supercritical regime (Cor 7.1). `check3_fakes.out`, Part A:
* the record 10-point cluster has `r = 5`, `rank B = 5`, free dimension `0` (pinned, as `round21/README.md` J3 found
  for the record clusters); the angle model at `phi = arg pi_p` reproduces the actual relative angles to `1.8e-15`;
* three same-profile angle models with `phi' = phi + t psi` (diameter ratio `0.95–1.004`, `|arg Z_f - arg Z| ≈ 1–2e-3`):
  `Z_f != Z`, all widths, heights and gcd data identical by construction, subtuple depths within `0.09`;
* the multiplicative relation lattice: 50 vectors `k` with `k' in {-1,0,1}^9`, `<k, a_p> = 0` for all `p`; `Z^k = 1`
  exactly for all 50, and every angle model satisfies them identically (`arg Z_f^k = sum_p 2<k,a_p>phi'_p = 0`). These
  monomial relations are forced by `L_intra` with `RHS = 0` (supports inside the relation lattice), and they are
  fake-blind (the rational-character class of `overview/overview.md` §5.1). Fake-violated relations need non-monomial
  `F` in a supercritical regime (Cor 7.1).

**Supercriticality at `M = 10` (evidence).** `check3_fakes.out` Part B scanned 104 slab-type supports
`{k : k_x in {-1,0,1}, <k,a_p> = 0 (p in P0), |<k,a_p>| <= b}`. Every support with `E > 0` lies in the multiplicative
relation lattice: `RHS = 0`, and the forced relations are the 50 monomial relations above, which fakes also satisfy.
`check3b_quotient.py` therefore works on the torus `T'` cut out by those relations. Its characters are
`l = (<k,a_p>)_p in Z^5`; 3875 of them are reached from `k' in {-1,0,1}^9`, and `theta(l) = <k,delta>` is consistent
mod `2pi`. The scan covers 190 boxes `|l_p| <= b_p in {0,1,2}`. It computes `reg` (mod-`p` rank, an upper bound),
the depth `L`, `RHS` and `E = reg·L + log reg! - RHS`. The result is **`max E = -1.5e-4 < 0`**. The extremal
supports are `{-l, 0, l}`, that is, `G = (x^{l/2} - x^{-l/2})^2`, which is the pair inequality at `|t_ab| ≈ 1`:
exactly critical, not supercritical (and `log|c|_1 = log 4` makes it strictly subcritical). So at the record's
largest cluster nothing beyond the monomial and pair relations is forced. This is consistent with route 5's need
for `n >= 64`, `n - 1 >= 300 ln 4r` (`r = 5` would need `n >= 900`). It says nothing about the asymptotic regime, and
fakes with the same profile pass every such test.

**Summary of the fake test.** The method's inputs are fake-blind except for exact relations. The exact relations that
exist at testable sizes are monomial, and those are fake-blind too. The fake-violated relations of Cor 7.1 exist
only in route 5's regime, where whether gap-admissible fakes satisfy them is the undecided CIRCLE(q) question. The
fake status of `L_intra` is therefore exactly route 5's.

## 9. Constraint audit, what remains open

**Owner's constraints.**
* *No uniform exceptional set from a fixed-`S` Subspace theorem.* None used. The only exceptional sets are zero sets
  of auxiliary polynomials and the base locus `B(V_H)`, which §7 explicitly refuses to treat as uniform.
* *No Gaussian gcd or denominator cost discarded.* Lemma 3.3 is the full product formula over `K`; the gcd-reduced
  tuples are handled exactly (Lemma 2.2: the toric bound is at least as strong as integrality and the common factor
  `prod g^a gbar^b` is tracked); route 5's `N^{D/2}` is tracked in Lemma 4.0. Inert primes and `1+i` contribute `0`
  (Lemma 2.1 of `round24/theory.md`), a gain not used, not a cost dropped.
* *No drifting common phase replaced by zero.* All auxiliary objects are characters with `sum s = 0` on `T_5` or
  `T_M`, hence independent of the representative; Lemma 2.1(1) verifies this exactly for the reduced tuples. Nothing
  fixes `P_1 = 1` or `z_1 = 1`: `Z = (z_j/z_1)` is a point of the quotient torus, and every evaluation used is of a
  degree-0 character.
* *No template-dependent threshold treated as uniform in the radius.* Theorem 5.1 states its constants `C(F)`, `C(G)`
  explicitly; Cor 5.2(2) gives the threshold on `w` as a function of `C(F)`; Prop 6.1(d) keeps `|C(F)|/(rho d_1 w_1)`;
  §7 identifies the base locus as template-dependent and does not use it as a fixed `Z`. The route-5 regime constants
  (`n >= 64`, `n-1 >= 300 ln 4r`, `W_1`) are the refereed ones and depend on `n, r` only.

**Open (sharpened, nothing new opened).**
1. The zero estimate for the intra-cluster relations of Cor 7.1 = route 5's ZE at the actual point
   (ZE-col*; sufficient: MIN^pair / `(R_H)`). Open, unchanged.
2. Round 24's open point 1 for shared-prime configurations of A-points **not** coming from one cluster (e.g. points on
   different circles sharing primes). Theorem 5.1 covers every configuration that is the image `Phi(Z)` of a balanced
   simple point of some `T_m` under a homomorphism `Phi : T_m -> T_5^n` defined over `Q` (the proof uses only that); a
   general shared-prime `n`-tuple of `T_5(Q)` points need not be of that form with a balanced `Z`.
3. Non-simple balanced profiles (`e_p >= 2`): (5.1) needs `v_pi = 1_A`; for chains `sum_t 1_{A_t}` the width is not
   additive and Theorem 5.1 is not proved (as in `round24/theory.md`, which also assumes simple profiles).
4. Congruence divisibility at `q ∤ N` combined with intra-cluster polynomials (residue inputs) is outside `L_intra`.
