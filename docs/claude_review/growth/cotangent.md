# The integer-cotangent route: dictionary, rigidity at the Cilleruelo–Córdoba scale, direction pinning, and a split-scale growth bound

Growth-attempt write-up, cotangent route, 2026-10-08. All scripts are in
`growth/cotangent_checks/` (run them from that directory with `python3`; they use only
`numpy` and the standard library). Labels used below: **[proved here]**, **[known]** (with
the note in `research/docs/`), **[data]** (finite computation, not evidence of a theorem),
**[heuristic]**.

## 0. Status

* **The uniform theorem is not proved, and the general growth rate is not improved.** The
  best general bound is still the one recorded in the notes.
* **New rigorous results (all restricted)**:
  * **Theorem A [proved here].** Let `L_0` be the intrinsic cotangent scale of the cluster,
    that is, the lcm of all reduced half-angle-cotangent denominators. Suppose `L_0` has no
    split prime factor `p <= M/8`. Then `M <= (sqrt(192)+o(1)) sqrt(log R / log log R)`.
    If instead `v_p(L_0) <= B` for every split `p <= (M/8)^{1/(B+1)}`, then
    `M = O_B(sqrt(log R))`. The theorem comes from an exact inequality that trades the
    prime-power content of `L_0` against level spreading at split primes.
  * **Theorem B [proved here]: rigidity at the Cilleruelo–Córdoba (CC) scale, for every
    `M`.** Suppose `M` points have all mutual distances at most `K N^{theta_M}`, where
    `theta_M` is the CC exponent. Then three quantities are bounded by explicit functions of
    `K` and `M`: the product of all reduced cotangent denominators, the cut deficit, and
    `L_0`. For five points: `sigma * prod_{i<j} s_ij <= K^10/128` and `6 | L_0`. As a
    consequence, five lattice points always have diameter at least `1.943 R^(2/5)` (CC gives
    `sqrt 2`). Also, Cilleruelo–Granville's open question (can five points lie on arcs of
    length `O(R^(2/5))` infinitely often?) is *equivalent* to an integer problem with a
    **bounded** cotangent scale `L`.
  * **Theorem C [proved here]: direction pinning for even `M` at the CC scale.** The product
    of the points equals `N_bal^(M/2)` times a Gaussian integer of bounded norm. Hence the
    cluster lies within `O_K(N^(theta_M - 1/2))` of one of finitely many directions
    `theta`, each with `e^(2iM theta)` in `Q(i)`, and none of them a rational lattice
    direction once `N` is large. This was checked exactly on the four-point Pell family: for
    every member, `prod z_i` is a rational multiple of `-3+4i`.
* **Structured subclasses that were requested.** These are settled, but by results that
  are already known:
  * Arithmetic progressions: fixed projective shape, item 24.
  * `L` prime, a prime power, or not divisible by 6: at most three finite coordinates, by the
    local bounds.
* **Constructions and data.**
  * The five-row cyclotomic family (`C_* -> 0`) has `L_0 ~ N^0.239`, with one residue of size
    `~ N^0.23`.
  * The best five-point clique with `L_0 = 6` and `A <= 10^6` has `C_* = 1.187`. No other
    clique with `C_* <= 1.6` exists in that range.
  * An exhaustive search over `L ∈ 6Z` with `L <= 240`, `A <= 2*10^4`, and also over
    `L <= 120`, `A <= 10^5`, finds no five-point clique with diameter `< 4.94 R^(2/5)`.
    The largest `N` found with diameter `<= 8 R^(2/5)` is `6.8*10^11`.
  * Numerically, the rigid five-point system has no non-degenerate linear-block polynomial
    family.
* **Missing lemma** (Section 8):
  * For the uniform theorem: a slack-absorption lower bound for `Delta * prod s_e^2` that is
    linear in `W = log N`.
  * For the Cilleruelo–Granville five-point question: finiteness of the bounded-`L` rigid
    five-point system (equations (B4) and (B5)).

## 1. The dictionary, re-derived and checked on actual circles

**Setting.** Let `z_0,...,z_k` be distinct Gaussian integers of common norm `N = R^2`. Divide
out their Gaussian gcd. This does not increase the normalized arc. The norm `N` is then odd
and divisible only by split primes `p = pi conj(pi)`. Write
`z_j = eps_j prod_p pi_p^{a_jp} conj(pi_p)^{e_p - a_jp}`, with `min_j a_jp = 0` and
`max_j a_jp = e_p`. Put `M = k+1` and `W = log N`. Order the points counterclockwise inside a
short arc of angular width less than `pi`, and let `theta_ij` be the angle from `z_i` to
`z_j`.

**(D1) Chart.**
`c_ij = cot(theta_ij/2) = (N + <z_i,z_j>)/det(z_i,z_j)` is rational. Let `s_ij` be its
reduced denominator (the *residue*), and let `L_0 = lcm_{i<j} s_ij` be the intrinsic scale.
With the anchor at the endpoint `z_0`, the coordinates `X_j = L_0 c_0j` are positive
integers, and they decrease as the angle increases. The cotangent addition formula gives

    Q_ij := L_0 c_ij = (X_i X_j + L_0^2)/(X_i - X_j)  in Z,
    q_j  := z_j/z_0  = (X_j + i L_0)/(X_j - i L_0).

Any `L` that clears all the `c_ij` is a multiple of `L_0`. Common scaling `(X, L) -> (uX, uL)`
changes nothing below, so we write `L` for any admissible scale.

**(D2) Radius.**
`N = lcm_e n_e`, where for each of the `binom(M,2)` edges

    g_e = gcd(Q_e, L),
    n_e = ((Q_e/g_e)^2 + (L/g_e)^2)/eps_e,
    Q_0j = X_j.

Proof. At a split prime `p`, put `sigma_j = v_pi(q_j) = a_jp - a_0p`. Each reduced edge
denominator has `p`-adic valuation `|sigma_i - sigma_j|`. Primitivity makes the maximum of
these valuations equal to `e_p`. Inert primes and `1+i` never occur in `n_e`. This is the
argument of `integer_cotangent_lcm_height_target.md` §2, which was re-proved here.

**(D3) Arc.** The span is `2 arctan(L/A)`, where `A = min X_j`. Hence

    C_* := arc/sqrt(R) = 2 N^{1/4} arctan(L/A),
    C_*^4 ~ 16 N/(A/L)^4.

**(D4) Chords.**

    |z_i - z_j|^2 = 4 N L^2/(Q_ij^2 + L^2) = 4 N s_ij^2/(eps_ij n_ij),
    n_ij = N/Norm(gcd(z_i, z_j)).

Here `eps_ij = 2` exactly when `z_i` and `z_j` have x-coordinates of different parity.
Consequently `prod_e eps_e <= 2^{floor(M^2/4)}`.

**(D5) Good-prime dictionary [proved here].** Let `p ∤ 2L`, and put
`alpha_j = v_pi(X_j + iL)` and `gamma_j = v_{conj pi}(X_j + iL)`. Then:

* `min(alpha_j, gamma_j) = 0`;
* `v_p(N) = max_j alpha_j + max_j gamma_j`;
* `a_jp = max gamma + alpha_j - gamma_j`;
* `v_p(X_j^2 + L^2) = alpha_j + gamma_j`;
* under the clique condition, `v_p(X_i - X_j) = min(alpha_i, alpha_j) + min(gamma_i, gamma_j)`.

At a good prime, therefore, the clique condition says exactly this: two coordinates that
agree mod `p` are isotropic for the *same* Gaussian orientation, and they agree to exactly
their common depth. The second clause is the statement that same-level points do not
collide. At good primes the integer condition carries **no information beyond the Gaussian
allocations**. Every genuine restriction sits at the primes dividing `2L`, that is, in the
residues. This is why the cotangent route can only reformulate the residue and cut
bookkeeping of the Gaussian route (Section 8).

**(D6)** The pair gcd formula `Norm(gcd(X_i+iL, X_j+iL)) = |X_i - X_j| d`, with `d | L`
**[known]**.

**(D7) Endpoint swap.** `X -> A + (A^2+L^2)/(X-A)`, `infinity -> A`, is an involution. It
preserves `N`, `L` and `A` **[known]**.

**(D8) Cut identity.**

    (prod_e n_e)^4 * prod_p p^{4 D_p} = N^{M^2},
    D_p = sum_{tau=1}^{e_p} (S_tau - M/2)^2,
    S_tau = #{j : a_jp < tau}.

Verification: `check_dictionary.py` checks (D1)–(D8) exactly, including the parity rule and
the good-prime formulas, on more than 2500 consecutive-point clusters cut from random
primitive circles with up to five split primes. `check_theorems.py` checks (D8) in the
integral fourth-power form, together with the cotangent form of `prod_e n_e`, on 3510
clusters. Both pass.

## 2. The master identity and what (H) asks for

Combining (D4) and (D8) gives the following exact identity **[proved here]**. Here
`P = binom(M,2)`. The finite-pair norm satisfies
`Q_ij^2 + L^2 = (X_i^2+L^2)(X_j^2+L^2)/(X_i-X_j)^2`.

    (M1)  prod_{i<j} |z_i - z_j|^2 = 4^P N^{P - M^2/4} prod_p p^{D_p} prod_e s_e^2 / prod_e eps_e,

    (M2)  prod_e n_e = prod_{i=1}^k (X_i^2+L^2)^k / ( prod_{i<j}(X_i - X_j)^2 prod_e g_e^2 eps_e ).

Let `Q_* = floor(M^2/4)`, and let `Delta_* = prod_p p^{D_p - e_p(M^2/4 - Q_*)}`. This is an
integer. For odd `M` its exponent is `sum_tau [(S_tau - M/2)^2 - 1/4]`. It equals `1` iff
every layer is an extremal cut (`M/2|M/2`, or `(M±1)/2`). In the endpoint chart (M1)–(M2)
become

    (M3)  N^{Q_*} = Delta_* prod_e s_e^2 * prod_i (X_i^2+L^2)^k / (L^{2P} prod_{i<j}(X_i-X_j)^2 prod_e eps_e).

Sort the coordinates as `A = X_1 < ... < X_k`. Then
`prod_{i<j}|X_i - X_j| <= prod_i X_i^{i-1}`, which gives
`prod_i (X_i^2+L^2)^k / prod_{i<j}(X_i-X_j)^2 >= A^{2P}`. Hence

    (M4)  N >= (1/2) (A/L)^{2P/Q_*} (Delta_* prod_e s_e^2)^{1/Q_*},
          2P/Q_* = 4M/(M+1)  (M odd),   4(M-1)/M  (M even).

With `Delta_* = s_e = 1` this is exactly the Cilleruelo–Córdoba exponent: 10/3 for five
points and 7/2 for seven. These are the exponents of `affine_content.md` (15). The new
feature is the factor `(Delta_* prod s_e^2)^{1/Q_*}`. **By (M4), (H) follows from the
slack-absorption inequality**

    (SA)  Delta_* prod_e s_e^2  >=  c_M (A/L)^{2 floor(M/2)}  (~ N^{floor(M/2)/2}).

Conversely, a cluster on an arc of normalized length `C` satisfies

    Delta_* prod_e s_e^2 <= (C^2/4)^P prod_e eps_e N^{floor(M/2)/2},

by (M1) and `|z_i - z_j| <= C N^{1/4}`. So the full slack `N^{floor(M/2)/2}` must be
absorbed by cut deficits and residues. The converse direction, from (H) to (SA), holds
only after dividing (SA) by the shape factor
`prod_i (X_i^2+L^2)^k/(A^{2P} prod_{i<j}(X_i-X_j)^2) >= 1` of (M3). (SA) is therefore a
sufficient target, and it is equivalent to (H) on configurations of bounded shape.

## 3. Theorem B: rigidity at the Cilleruelo–Córdoba scale

Let `theta_M = (P - Q_*)/(2P)`, that is, `(M-1)/(4M)` for odd `M` and `(M-2)/(4(M-1))` for
even `M`. Thus `theta_4 = 1/6`, `theta_5 = theta_6 = 1/5` and `theta_7 = 3/14`; in terms of
`R`, the scales are `R^{1/3}` and `R^{2/5}`.

**Theorem B [proved here].** Suppose `M` points of a primitive circle have all mutual
distances at most `K N^{theta_M}`. Then

    (B1)  Delta_* * prod_{i<j} s_ij^2  <=  K^{2P} prod_e eps_e / 4^P  <=  K^{2P} 2^{Q_*} / 4^P.

Proof. By (M1), `prod |z_i-z_j|^2 = 4^P N^{P-Q_*} Delta_* prod s^2 / prod eps`, and the left
side is at most `K^{2P} N^{2P theta_M} = K^{2P} N^{P - Q_*}`. Then use
`prod eps <= 2^{Q_*}`. ∎

Consequences.

1. `L_0 <= prod s_ij <= K^P 2^{Q_*/2}/2^P` is bounded.
2. The non-extremal cut layers have total norm at most `Delta_*`, which is bounded. So all
   but boundedly many layers are extremal cuts.
3. A lower bound on the constant: `K >= 2^{1 - Q_*/(2P)} (prod s_ij^2 Delta_*)^{1/(2P)}`.

**Five points.**

* For `M = 5`, `Delta_* = sigma^2`, where `sigma` is the product of `p` over the singleton
  layers. (B1) reads

      (B2)  sigma * prod_{i<j} s_ij <= K^10/128.

* The local bounds force `6 | L_0` for every five-point cluster. The inert bound with
  `p = 3` gives `m <= 3` when `3 ∤ L`. The dyadic bound gives `m <= 3` when `L` is odd. Both
  are **[known]**, `integer_cotangent_local_bounds.md`. Checked: `6 | L_0` held in all 1411
  random five-point clusters in `check_theorems.py`.
* Hence `prod s_ij >= 6`, and **every five lattice points on a circle of radius `R` have
  diameter at least `(2^14 * 36)^{1/20} R^{2/5} = 1.943... R^{2/5}`**. CC's statement is
  that an arc of length `sqrt 2 R^{2/5}` holds at most four points. Since diameter is at
  most arc length, this improves the constant from `1.414` to `1.943`.
* For general `M`, inert-prime collision counting (`inert_prime_cofactor_bound.md` §3)
  gives `log prod s_ij >= (1/4 + o(1)) M^2 log M`. So `K >= (c + o(1)) sqrt M` at the CC
  scale. This is a routine transfer of that note's argument.

**Block balance (M = 5).**
Let `beta_T = log b_T` be the log-norm of the doubleton block `T`. (B2) and the
pair-distance bounds `d_ij >= (3/5)W - log(K^2/2)` force `d_ij = (3/5)W + xi_ij` with
`|xi_ij| <= 9 log(K^2/2) + 10 log K`. The doubleton parts satisfy
`(I + Pet) beta = W_dbl 1 - d^dbl`, where `Pet` is the Petersen (disjointness) adjacency on
the ten doubletons. This matrix inverts as `(I+Pet)^{-1} = Pet/2 - J/8`, with
`infinity`-norm 2. Therefore

    (B3)  |log b_T - W/10| <= 60 log K    (K >= 2).

**Reduction of the Cilleruelo–Granville five-point question [proved here].** The question
is whether five lattice points lie on arcs of length `O(R^{2/5})` for infinitely many `R`.
Take the endpoint chart of a five-point cluster with diameter at most `K R^{2/5}`.

* `L = L_0 <= K^10/128` and `6 | L`.
* All four finite coordinates satisfy `X_i ∈ [A, c_K A]`. This follows from (B3) and (D4):
  the shortest chord is `>= 2^{1/2}(2/K^2)^{9/2} N^{1/5}`.
* `N <= (K sqrt((A/L)^2+1)/2)^{10/3}`.

Conversely, every integer 4-clique with `N <= c (A/L)^{10/3}` gives five points on an arc of
length `<= 2c^{3/10} R^{2/5}`. Hence

> **CG's five-point question ⟺ for some fixed `L ∈ 6Z`, there are infinitely many integer
> 4-cliques `A = X_1 < X_2 < X_3 < X_4` (with `X_i - X_j | X_i^2 + L^2`) whose all-edge lcm
> satisfies `N <= c (A/L)^{10/3}`.**

In Gaussian form (the rigid model), the question asks for ten blocks `C_T`, indexed by the
edges of `K_5`, together with bounded integers `t_ij`, such that

    (B4)  Im( V_i^{(j)} conj(V_j^{(i)}) ) = t_ij,   V_i^{(j)} = prod_{l != i,j} C_il,

up to a bounded fudge factor coming from the singleton and bad-prime layers. Equivalently,
for every triangle `{i,j,k}` with complement `{l,m}`:

    (B5)  t_jk conj(C_jk) C_il C_im + t_ki conj(C_ki) C_jl C_jm + t_ij conj(C_ij) C_kl C_km = 0.

(B5) is the Plücker relation for the vectors `U_i = prod_{l != i} C_il`, divided by
`C_ij C_jk C_ki`. The points are recovered as `z_i = (fudge) conj(G) U_i/conj(U_i)`, where
`G = prod_T C_T`.

**Heuristic.** There are `N` choices of blocks. There are four independent angle
conditions, each of precision `N^{-3/10}`. The expected number of solutions is therefore
`N^{1 - 6/5} = N^{-1/5}`, which predicts finitely many solutions for each `K`. The same
count gives `N^0` for four points at `R^{1/3}`; there the Pell families exist, as
Cilleruelo–Granville classified.

**Data.**

* *Exhaustive search, `search_cg.py`.*
  * Every 4-clique with `L ∈ {6,12,...,240}`, least coordinate `A <= 20000` and
    `diam <= 6 R^{2/5}`. 528 cliques were found.
  * The smallest diameter ratio is `K5 = 4.9475`, at `N = 65`, `L = 12`,
    `X = (12,18,24,44)`.
  * The largest `N` with `K5 <= 6` is `10,414,625` (`L = 6`, `X = (282,333,594,1608)`,
    `K5 = 5.42`).
  * A second exhaustive run covered `diam <= 8 R^{2/5}`, `L ∈ {6,...,120}` and
    `A <= 10^5` (`cg_K8_L120_A1e5.txt`). It found 1476 cliques. The minimum is still
    `K5 = 4.9475`, and the largest `N` is `683,617,125,425` (`L = 36`). Nothing approaches
    the theoretical floor 1.943, and `K5` does not decrease as `N` grows. This is consistent
    with the heuristic, but it is data, not proof.
* *Polynomial families.* Take blocks linear in a parameter `n`, `C_T = alpha_T (n + gamma_T)`.
  (B4) with constant `t_ij` is equivalent to: `Im sum_{T ∋ i} gamma_T^k` is independent of
  `i` for `k = 1,...,5`. That is 20 real equations in 18 effective unknowns.
  * In 450 Levenberg–Marquardt runs (`poly_family_lm.py`), every exact solution was
    degenerate. In each, some `gamma_T` is real, which gives a rational content that cancels
    from `U/conj U`, or some pair of `gamma`'s is equal or conjugate. Such solutions collapse
    points or freeze the configuration.
  * No generic solution was found. This is numerical evidence, not a proof.

## 4. Theorem C: direction pinning for even M at the CC scale

**Lemma (layer bookkeeping) [proved here].** For any primitive cluster of `M` points,

    prod_{i=1}^M z_i = eps * prod_{layers (p,tau)} p^{min(h, M-h)} * Phi,
    Norm(Phi) = prod_{layers} p^{|M - 2h|},

where `h` is the number of points on the high side of the layer and `eps` is a unit.

Proof. A layer contributes `pi` to its `h` high points and `conj(pi)` to the others. Their
product is `p^{min(h,M-h)}` times `pi^{h-(M-h)}` or `conj(pi)^{(M-h)-h}`. ∎
Checked on 470 actual clusters (`check_pinning.py`).

**Theorem C [proved here].** Let `M >= 4` be even. Suppose the primitive cluster has all
distances at most `K N^{theta_M}`, so that by Theorem B, `Delta = prod p^{D_p}` is bounded.
Balanced layers have `|M - 2h| = 0`. Each unbalanced layer has
`|M-2h| = 2|h - M/2| <= 2(h-M/2)^2`. Hence

    prod z_i = eps N_bal^{M/2} Phi,   Norm(Phi) <= Delta^2 <= (K^{2P} 2^{M^2/4} 4^{-P})^2.

All points lie within angle `delta = 2 arcsin(K N^{theta_M - 1/2}/2)` of `z_1`, so

    | M arg z_1 - arg(eps Phi) |  <=  (M-1) delta    (mod 2 pi).

Consequences.

* For each `K`, CC-tight even clusters lie within `O_{K,M}(N^{theta_M - 1/2})` of a
  **finite set of directions** `D_K = {theta : e^{iM theta} = eps Phi/|Phi|}`. The same
  holds for non-primitive clusters, since the Gaussian gcd is bounded in terms of `K`.
* **Rational directions are excluded [proved here].** Suppose `theta` is the direction of a
  Gaussian integer `omega`. A point with `|z| = R` within angle `eta` of `theta` has
  `Re(z conj omega) ∈ [R|omega| cos eta, R|omega|]`. When `R|omega| eta^2 < 2` this
  interval contains at most one integer, and so the arc contains at most 2 lattice points.
  Here `eta ~ N^{theta_M - 1/2}` and `R eta^2 ~ N^{2theta_M - 1/2} -> 0`. Hence for
  `N >= N_0(K)`, every pinned direction is **irrational**, with `e^{2iM theta}` in `Q(i)`
  but `e^{i theta}` not proportional to any `omega` in `Z[i]`.

**Exact verification on the four-point Pell family.** This is the family of
`integer_cotangent_lcm_height_target.md` §5: `L = 24`, `U + V sqrt 5 = (9+4 sqrt 5)^n`.

* For `n = 1, 11, ..., 61`, the primitive part of `z_0 z_1 z_2 z_3` is exactly `-3+4i`.
  Here `N_u = 5` (one singleton layer) and `Phi = (1+2i)^2`.
* `4 arg z_0 - arg(-3+4i)` vanishes to machine precision for `n >= 11`.
* So the family is pinned to `theta = arg(1+2i)/2`, which has `tan theta = (sqrt 5 - 1)/2`.
  The notes had observed this limiting direction for this one family (item ~344). Theorem C
  explains it and applies to all even `M`.

**Six points.** By Theorem B, a six-point cluster of diameter `<= K R^{2/5}` is rigid: all
but boundedly many layers are `3|3` cuts, and the residues are bounded. By Theorem C it is
pinned to finitely many irrational directions with `e^{12 i theta}` in `Q(i)`. In block form:

* Use the ten `3|3` cuts, with sides containing point 1, labelled by edges `{x,y}` of `K_5`
  on `{2,...,6}`.
* Then `z_1 = (fudge) P` with `P = prod beta_xy`, and the other five points form a rigid
  five-point configuration with blocks `beta`.
* The sixth point adds exactly one condition: the common direction `phi` of the five
  `U_j = prod_{y} beta_jy` satisfies `3 phi ≡ (bounded phase) + O(N^{-3/10})`.

So a CC-tight six-point cluster is a CC-tight five-point cluster plus one more condition of
precision `N^{-3/10}` (heuristically `N^{-1/2}` many). This did not produce a contradiction.

**Why exactness does not follow.** One could hope that some product of blocks is forced to
be *exactly* parallel to a bounded direction. That would happen if a product
`Gamma = prod beta_e^{n_e}` with `sum|n_e| < 6` were known to within angle `O(N^{-3/10})`,
because `|Gamma| < N^{3/10}`. The relations available all come from `arg z_i ≈ arg z_j`.
They span the lattice of edge vectors `(sum a - a_x - a_y)_{xy}` with `a ∈ Z^5`. A complete
search over `a ∈ [-4,4]^5` shows that its shortest nonzero vectors have `l^1`-norm exactly 6.
These are the relations "`Y_j = prod_{e not containing j} beta_e` is nearly real", whose
imaginary parts are the bounded nonzero residues. So the rigid six-point system sits
**exactly at the exactness threshold**. This is the structural reason the method stalls.

## 5. Theorem A: the split part of the cotangent scale

**Level capacity [known; re-proved].** Let `p` be a split prime. Points at one allocation
level are pairwise incongruent modulo `p^{v_p(L_0)+1}` after removing the level factor. So
each level holds at most

    c_p = (p-1) p^{v_p(L_0)}

points (`split_prime_cotangent_capacity.md`). The proof: embed `Q(i)` in `Q_p` with
`i -> r`. On one level `q_j = p^h u_j` with `u_j` a unit, and
`c_ij = r(u_i+u_j)/(u_i-u_j)`. If `u_i ≡ u_j (mod p^{e+1})`, then
`v_p(c_ij) <= -e-1`, while `L_0 c_ij` is an integer.

**Level-spreading deficit [proved here].** Put `n_a <= c` points on consecutive levels, with
`M` points in total. Then

    D = sum_tau (S_tau - M/2)^2  >=  (M - 4c)_+^3 / (12 c).

Proof. The partial sums `S_tau` increase in steps of at most `c`, starting `<= c` and ending
`>= M-c`. So every interval `[M/2 + jc, M/2 + (j+1)c)` with `(j+2)c <= M/2` contains some
`S_tau`, and so does its mirror image. Hence `D >= 2 c^2 sum_{j=1}^{J-1} j^2` with
`J = floor(M/(2c))`, and `D >= (2/3)c^2(J-1)^3 >= (M-4c)^3/(12c)` for `M >= 4c`. ∎
An exact dynamic programme for `M < 90` and all `c` finds no violation; the worst ratio of
`D` to `M^3/(12c)` over `M >= 8c` is `0.656` (`check_level_deficit.py`).

**Theorem A [proved here].** Every cluster of `M` points on a primitive circle, with all
distances `<= C N^{1/4}`, satisfies

    (A1)  sum_{p ≡ 1 (4)} log p * (M - 4 c_p)_+^3 / (12 c_p)  <=  (M/4) W + binom(M,2) log(C^2/2),
          c_p = (p-1) p^{v_p(L_0)}.

Proof.

1. Plotkin (D8): `sum_{i<j} log n_ij = W M^2/4 - sum_p D_p log p`.
2. By (D4) with `s >= 1` and `eps <= 2`, `log n_ij >= W/2 + log(2/C^2)`.
3. Subtracting gives `sum_p D_p log p <= WM/4 + binom(M,2) log(C^2/2)`.
4. Insert the deficit bound. A split `p` with `c_p < M` is forced to divide `N`, because one
   level can hold at most `c_p` points. ∎

Checked on 1883 clusters (`check_theoremA.py`).

**Corollaries.**

* **(A2)** Suppose `v_p(L_0) = 0` for every split `p <= M/8`. Then
  `W >= (M^2/24) sum_{p<=M/8, p≡1(4)} log p/(p-1) - 2(M-1)log(C^2/2) = (1/48 - o(1)) M^2 log M`.
  This uses Mertens' theorem in the progression `1 mod 4`. Since `W = 2 log R`, this gives
  `M <= (sqrt 192 + o(1)) sqrt(log R/log log R)`.
* **(A3)** Suppose `v_p(L_0) <= B` for every split `p <= (M/8)^{1/(B+1)}`, with `B >= 1`.
  The sum `sum log p/p^{B+1}` converges, so `W >= c_B M^2` and `M = O_B(sqrt(log R))`.

**Interpretation.** Small split primes must either collide inside a level, which forces
them into `L_0` and costs residue, or spread over levels, which costs a cubic cut deficit.
The general `(2/3)` chain lemma optimizes this trade-off when collisions are allowed. (A2)
is the extreme case with no collisions allowed. Theorem A does **not** improve the general
rate: a cluster can put high powers of the small split primes into `L_0` at negligible cost,
because one colliding pair suffices to raise `v_p(L_0)`.

## 6. The structured subclasses that were requested

* **Arithmetic progressions `X_j = A + jD` [known: item 24, fixed projective shape].**
  * With `G_j = (A+iL) + jD`, the points are `q_j = G_j/conj(G_j)`, a fixed-shape family
    with `u_j = j`, `v_j = 1`, and anchor `alpha = D`.
  * After normalizing `gcd(alpha, beta) = 1`, three facts follow:
    * `gcd(G_i, G_j)` divides `Delta_ij = j - i`;
    * `N >= prod Norm(G_j')/(2^{k+1} prod|Delta|^{2k+2})`, where `G'` is the primitive part
      of `G`;
    * the product of the contents divides `|Im(alpha conj beta)|` times a bounded factor.
  * Multiplying the chord inequalities `2|Delta_ij| D <= C N^{-1/4} c_i c_j |G_i'||G_j'|`
    over all pairs then gives:
    * `C_* >= c_k > 0` for `k >= 3`;
    * finitely many AP cliques, up to scaling, for `k >= 4` and fixed `C`.
  * I re-derived this independently. It is the item-24 theorem specialized to
    `s = {infinity, 0, 1, ..., k-1}`.
* **`L` prime, a prime power, or bounded [known: local bounds].** Four finite coordinates
  need `6 | L`. A prime or prime-power `L` therefore allows `k <= 3`. More generally, a
  bounded `L` gives `k <=` the least inert prime not dividing `L`.
* **Prime-power differences.** I found no structure that interacts with the all-edge lcm.
  The condition `X_i - X_j = p^a` puts all of a pair's overlap into one Gaussian prime. This
  is a two-valued-model statement already covered by the cut accounting. No result is
  claimed.

## 7. Symmetries, and working backwards from an unbounded family

**Symmetries.**

* A rational Möbius map preserves every clique and its lcm only if it is induced by
  relabelling (re-anchoring) or by reflection. The points are relative phases, and a
  non-isometric map changes `N` (item 53, item 54).
* Among endpoint charts, the only map preserving `(N, A/L)` is the endpoint swap (D7), an
  involution.
* Galois conjugation of `Q(i)` is `X -> -X`.
* Rotation by a non-cluster rational angle moves the anchor, so it is not a chart map.
* Translations `X -> X + t` with closure conditions exist, but they change `N` with no
  monotonicity (`integer_cotangent_translation_obstruction.md`).

Conclusion: **symmetry supplies no descent**.

**Necessary structure of a hypothetical unbounded family** (`M -> infinity` at fixed `C`).
Each item below is proved here or in the notes.

1. Slack (Section 2): `sum_p D_p log p + 2 sum log s_e <= WM/4 + O_C(M^2)`.
2. Inert radical: `v_p(L_0) >= log_p((M)/(p+1))` for inert `p < M`, so
   `log L_0 >= (1/2+o(1))M` (`integer_cotangent_inert_radical_obstruction.md`).
3. Split trade-off, Theorem A: the small split primes enter `L_0` to depth about
   `log_p(M/p)`, or else `W >> M^2 log M`.
4. Cross-ratio heights `R^{1/4+o(1)}` (item 25). So the configuration is not of bounded
   projective shape, and AP-type subclasses are excluded.

**The first contradictions that can actually be proved here:**

* (i) `L_0` bounded implies `M` bounded (inert radical; trivial).
* (ii) At the CC scale, Theorems B and C exclude every cluster pinned to a rational
  direction.

**Constructions showing that the natural next steps fail.**

* **"`C_* -> 0` forces small residues" fails.** I recomputed the five-row cyclotomic family
  (`cyclo_scale.py`, `k = 19`, `d = 1,...,11`). It has `log L_0 / log N = 0.2388`,
  `0.2387`, ... (stable). At `d = 1` the residues are `{1,1,1,e^9,...,e^40}` plus one residue of size
  `N^{0.23}`. For every `d`, six residues are `N^{0.0006}` to `N^{0.0026}` and one is
  `N^{0.23}`. Also `log C_*` decreases linearly in `d` (to `-42` at `d = 11`). Small
  normalized arcs can coexist with very large and very unequal residues.
* **"Small `L` forces `C_*` bounded away from 0" is not refuted.** The best five-point
  cliques found with `L_0 = 6` have:

  | range | best `C_*` |
  |---|---|
  | `A <= 2*10^4` | `2.258` |
  | `A <= 2*10^5` | `1.187` (`X = (163488, 221633, 435366, 546558)`, `N = 6.8*10^16`, `log N/log(A/L) = 3.796`) |
  | `A <= 10^6` (`deep_L6.txt`) | still `1.187`; no other clique with `C_* <= 1.6` |

  Every other `L ∈ 6Z` up to 120 gives either scaled copies or `C_* >= 1.996`
  (`many_m4.txt`). For four points (`many_m3.txt`, `A <= 10^5`), `C_*` is as small as
  `0.28` (`L = 2`). That is the Pell-type phenomenon. So `min C_*` at fixed `L` dropped from
  `A = 2*10^4` to `A = 2*10^5`, then did not move up to `10^6`. Neither the conjecture
  "fixed `L` and five points imply `C_* >= c_L`" nor its negation is supported strongly.
* **Pinning holds exactly in the only known CC-tight family** (Section 4).

## 8. The missing lemma, and why the method stalls

**For the uniform theorem.** By (M3)–(SA), the cotangent target (H) follows from the
slack-absorption inequality, and it is equivalent to it for configurations of bounded shape: `Delta_* prod_e s_e^2 >= c_M (A/L)^{2floor(M/2)}`. In
logarithmic form this reads

    sum_p D_p log p + 2 sum_{i<j} log s_ij  >=  (floor(M/2)/2) W - O_C(M^2)    (M >= M_0).

All known lower bounds for the left side are local. They come from:

* inert collision counting;
* split capacity and levels (Theorem A, the chain lemma);
* pair and four-row characters (`general_four_row_quartic_slack_obstruction.md`, the
  Paley obstruction).

These bounds depend on `M`, not on `W`. They give `(3/4 + o(1)) M^2 log M` at best. The
missing lemma must produce a lower bound that is *linear in `W` per point*. (D5) shows why
the cotangent integrality cannot supply it by itself: away from `2L` the clique condition is
equivalent to the Gaussian allocation data, and at `2L` it is exactly the residue
bookkeeping. The cyclotomic family shows the residue mass can concentrate on a single pair,
so any argument that averages residues over pairs has to cope with that concentration.

**For the CG five-point question (a cleaner and much smaller target).** Prove that the
bounded-`L` rigid system (B4)/(B5) has only finitely many solutions for each bound on the
`t_ij`. Equivalently, for each fixed `L ∈ 6Z` there are only finitely many integer 4-cliques
with `N <= c(A/L)^{10/3}`. This would improve CC's five-point bound to
`diam >= f(R) R^{2/5}` with `f -> infinity`. The method stalls on two points:

* (B5) is a system of three-term `S`-unit-type equations whose "units" are products of
  *unknown* blocks. Evertse-type uniform counts and Baker bounds do not apply when `S`
  varies with the solution. `abc` gives nothing, because the radicals are as large as the
  terms when the blocks are squarefree.
* For polynomial families the system is overdetermined (20 equations, 18 unknowns), but I
  have only numerical evidence that its solutions are degenerate. An elimination proof would
  need a computer-algebra Gröbner computation, which was not available.

## 9. Files (`growth/cotangent_checks/`)

| file | content |
|---|---|
| `cot_lib.py` | Exact Gaussian and cotangent library, written independently of the notes' checkers |
| `check_dictionary.py` | (D1)–(D8), parity rule, good-prime dictionary; random actual circles |
| `check_theorems.py` | (D8) fourth-power form, (M2), level capacity, (B2), `6 | L_0` |
| `check_five_diameter.py` | Brute force over `N <= 20000`: min diam/`N^{1/5}` over five points is 3.24 (at `N = 5`), above the bound 1.943 |
| `check_level_deficit.py` | Exact DP for the level-spreading deficit |
| `check_theoremA.py` | Chain of inequalities of Theorem A on actual clusters |
| `check_pinning.py` | Layer-bookkeeping identity; exact pinning of the Pell family to `-3+4i` |
| `search_fixedL.py`, `search_pruned.py`, `run_many.py`, `many_m3.txt`, `many_m4.txt`, `deep_L6.txt` | Fixed-`L` clique searches (`C_*`) |
| `search_cg.py`, `cg_K8_L120_A1e5.txt` | Exhaustive CG-regime search (diameter / `R^{2/5}`) |
| `cyclo_scale.py` | Intrinsic scale and residues of the cyclotomic five-row family |
| `poly_family_lm.py` | Levenberg–Marquardt exploration of linear-block polynomial families |
| `analyse_clique.py`, `pell_direction.py` | Anatomy of individual cliques; Pell direction |
