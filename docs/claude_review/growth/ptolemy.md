# Ptolemy matching relations and S-unit structure (route report)

Route: the quadruple (Ptolemy) relations `n_A a + n_B b = n_C c` of a primitive short-arc
cluster, treated unconditionally as an S-unit system, with the aim of forcing a residual-height
lower bound `H >= c w`.

## 0. Summary

**Verdict.** The route gives no unconditional growth improvement. Its information content can be
identified exactly. With that identification, the route proves a sharp conditional theorem at
`k = 8` points and a no-go theorem for `k <= 6`.

* **What the Ptolemy system is (Theorem B).** Take an anchored common-unit cluster as in item 532.
  The `C(k,4)` quadruple relations, with their exact matching coefficients, are the Plücker
  relations of the integer vectors `P_i`, after the common factors are cancelled. Every factor
  `conj(P_i)` cancels identically. Equivalently they are the S-unit equations
  `chi + (1 - chi) = 1` satisfied by the cross-ratios. The whole system has rank `(k-2)(k-3)/2`.
  Its content is precisely the rational point `[P_0 : ... : P_m]` of `M_{0,k}` together with its
  boundary contacts. At a prime `p`, the contact with the boundary divisor `D_T` is the `p`-adic
  tree edge length. That length equals `v_p(n_T)` exactly when `p` divides no residual, and in
  general differs from it by at most `2 max v_p(t_ij)`.

  Four things never enter the system: the singleton blocks, the top block `n_[m]`, the Gaussian
  norms `|P_i|^2`, and the isotropy of the block residues (`rho_T^2 = -1 mod n_T`). The system is
  therefore blind to the quadratic form: the split form (divisors in short intervals) and every
  other binary form give the same system. This sharpens two notes:
  `ptolemy_squareclass_plucker_audit.md`, where Ptolemy is shown equivalent to Plücker on
  rational slopes, and `plucker_compatibility_route.md`, which gives the collision-order barrier.
* **Answer to the overdetermination question (Section 4).** The system is never overdetermined
  algebraically, for any `k`. Its solution variety is the rational Plücker variety. Elimination
  does force polynomial identities among block norms with small coefficients; they are the
  anchored Plücker identities `(4.2)`, and they are always solvable. The obstruction is
  arithmetic, and it has a sharp threshold.

  Write `beta` for the balanced curve class on `Mbar_{0,k}`, the class with `beta.D_S = 1` for
  every boundary divisor. Then

  ```text
  K_{Mbar_{0,k}} . beta = 2^(k-3)(k-8) + k + 2,
  ```

  which is negative for `k <= 7` and equals `+10` at `k = 8`. The integer heuristic exponent per
  shift class is exactly `1 - K.beta`, namely `6, 9, 8` at `k = 5, 6, 7` and `-9` at `k = 8`.
* **Theorem C (no-go below 7, unconditional).**
  - For `k = 5` there is an explicit polynomial family. In it all residues are bounded
    (`|Y_i| = 1`, `|t_ij| <= 144`), the 10 visible blocks are pairwise coprime, and every block
    has `log n_T = w + O(1)`. All five quadruple relations hold as polynomial identities.
  - For `k = 6`, item 109 of the notes, combined with an anchoring lemma proved here
    (Lemma 5.2), gives the same conclusion.

  So no argument that uses only the following can force `H -> infinity` for `k <= 6`:
  - the quadruple relations;
  - the exact matching-coefficient structure (coprimality and balance);
  - the cyclic order;
  - the integrality of the residues.

  For `k = 7` the expected dimension is `+8`. My numerical search found no nondegenerate
  balanced rational curve in `Mbar_{0,7}` (Section 5.4), so existence remains open.
* **Theorem D (conditional, sharp).** Assume Vojta's Main Conjecture with `D = 0` for the single
  fixed variety `Mbar_{0,8}` over `Q`, that is, `h_K <= eps h_A + O(1)` off a proper Zariski-closed
  set. Then the uniform theorem holds: every arc of length `C sqrt(R)` contains at most `M(C)`
  lattice points. The proof uses only the Ptolemy (`M_{0,8}`) data of 8-point sub-tuples, the
  item-532 extraction, and a counting lemma that avoids the exceptional set. At `k = 8` it gives
  `H_Sigma >= (6w - O(1))/69`.

  For `k <= 7`, `-K_{Mbar_{0,k}}` is effective and supported on the boundary. Hence the
  corresponding instance of Vojta's conjecture is trivially true and carries no information.
  `k = 8` is the first nontrivial case.

  A different conditional theorem in `growth/conditional.md` uses the same form of Vojta's
  conjecture on another 7-fold, together with the Gaussian norm structure. Theorem D shows that
  the form-blind Ptolemy data alone already suffice, and that `k = 8` is exactly where they start
  to suffice.
* **Missing lemma (Section 7).** Let `(L_8)` be the following statement. For some `c > 0` and
  `eta_0 > 0`, every balanced Boolean Plücker configuration of 8 points, with relative error
  `eta <= eta_0` and lying off a fixed proper Zariski-closed subset of `M_{0,8}`, has
  `H_Sigma >= c w - O(1)`. Then `(L_8)` implies the uniform theorem (Theorem D, steps 1-5).
  Vojta's conjecture implies `(L_8)`, since it is exactly the inequality `h_K <= eps h` at these
  points. No unconditional tool reaches it. The 3-term S-unit and `abc` forms of the relations are
  vacuous on squarefree balanced blocks. The Subspace and Ru–Vojta theorems need a fixed `S`,
  while here `S` is the set of primes of `N`. Fixed polynomial eliminations are capped by the
  collision-order barrier.

**New relative to the notes**, as far as grep shows:

* the mixed-unit Gaussian form of Ptolemy (Theorem A(3), (5));
* the exact `M_{0,k}` contact dictionary for the nested extraction, with the Lipschitz bound
  (Theorem B(iv));
* the anchoring lemma (Lemma 5.2);
* the explicit `k = 5` bounded-residue family (Theorem C1);
* the canonical-class computation `K.beta = 2^(k-3)(k-8)+k+2` and the threshold `8`;
* Theorem D.

Known before: the common-unit relation (`ordered_residue_growth.md` (5)); Ptolemy equals Plücker
(`ptolemy_squareclass_plucker_audit.md`); the CRT large-residue realization and the
collision-order barrier (`plucker_compatibility_route.md`); cut-norm heights of cross-ratios
(`coupled_crossratio_height_norm.md`); and the elliptic `k = 6` contact profile (item 109).

All identities below were checked exactly by computer on actual circles (Section 8). Numerical
(floating-point) experiments are labelled as such and are not used in any proof.

---------------------------------------------------------------------------------------------

## 1. Setting and notation

Fix `0 < C_0 < sqrt 2`. Retain one Gaussian-unit class of a short-arc cluster and divide by its
Gaussian gcd, as in item 532 (`exact_nested_profile_residual_extraction.md`). Then `N = R^2` is
odd with only split prime factors, `W = log N`, and every point has the form

```text
z_a = eps * prod_p pi_p^{a_a(p)} conj(pi_p)^{e_p - a_a(p)},   0 <= a_a(p) <= e_p .
```

For an ordered tuple `(z_0, ..., z_m)`, with `k = m+1` and anchor `z_0`, item 532 defines the
following data.

* The anchored numerators
  `P_i = prod_p pi_p^{(a_i - a_0)_+} conj(pi_p)^{(a_0 - a_i)_+} = X_i + i Y_i`, with `P_0 = 1`.
  They are conjugate-primitive, and `z_i/z_0 = P_i/conj(P_i)`.
* The blocks `H_T`, for nonempty `T` contained in `[m]`, built from threshold layers, with
  `n_T = Norm(H_T)`. Then `P_i = prod_{T contains i} H_T`, and
  `G_ij := Norm(gcd(P_i, P_j)) = prod_{T contains i,j} n_T`.
* The residues `t_ij = det(P_i, P_j)/G_ij`, which are nonzero integers, and
  `t_0j := Y_j = det(P_0, P_j)`. Here `det(u, v) = u_1 v_2 - u_2 v_1`.
* The aggregate residual height `H_Sigma = sum_j log|Y_j| + sum_{i<j} log|t_ij|`, together with
  `w = W/2^m`.

At one rational prime, the blocks with `p | n_T` form two nested chains on disjoint row sets.
Blocks may share rational primes.

**Boundary splits.** A *boundary split* of `{0, ..., m}` is a pair `{S, S^c}` with
`|S|, |S^c| >= 2`. It is labelled by its side `T` not containing `0`, so
`T` is contained in `[m]` and `2 <= |T| <= m-1`. There are `2^m - m - 1` of them. The remaining
nonempty `T` are *invisible*: the singletons `T = {i}` and the top block `T = [m]`.

## 2. Theorem A: the general Ptolemy matching relation

**Theorem A.** Let `z_1, z_2, z_3, z_4` be distinct lattice points of `x^2 + y^2 = N`. The units
`eps_a` are arbitrary, and the exponents and allocations `a_a(p)` are arbitrary. For a pair put

```text
g_ij = prod_p pi_p^{min(a_i,a_j)} conj(pi_p)^{e_p - max(a_i,a_j)},
h_ij = prod_{a_i > a_j} pi_p^{a_i - a_j} prod_{a_i < a_j} conj(pi_p)^{a_j - a_i},
v_ij = eps_i h_ij - eps_j conj(h_ij).
```

At each prime sort the four allocations `b_1 <= b_2 <= b_3 <= b_4`. Put `m_p = b_3 - b_2`, and let
`kappa(p)` be the matching that pairs the two lower points with each other and the two upper
points with each other. Set `X_k = prod_{kappa(p) = M_k} p^{m_p}` for the three matchings
`M_1 = 12|34`, `M_2 = 13|24`, `M_3 = 14|23`. Then:

1. `z_i = eps_i g_ij h_ij`, `z_j = eps_j g_ij conj(h_ij)`, and `z_i - z_j = g_ij v_ij`.
2. **Exact common factor.** With `Gamma = prod_p pi_p^{b_1+b_2} conj(pi_p)^{2e_p - b_3 - b_4}`
   in `Z[i]`:

   ```text
   g12 g34 = Gamma X_1,   g13 g24 = Gamma X_2,   g14 g23 = Gamma X_3   (equalities in Z[i]).
   ```

   Equivalently, `X_k` is the product of the norms of the threshold layers whose cut restricts
   to a 2|2 cut of the quadruple of type `M_k`.
3. **Gaussian Ptolemy.** `X_2 v13 v24 = X_1 v12 v34 + X_3 v14 v23` holds exactly in `Z[i]`.
4. Suppose the points are in cyclic order. Then `v12 v34`, `v14 v23` and `v13 v24` are positive
   integer multiples `a, b, c` of one primitive Gaussian integer `omega`, and

   ```text
   c X_2 = a X_1 + b X_3,      with X_1, X_2, X_3 pairwise coprime.
   ```

5. Write `eps = eps_j/eps_i`. Then `|v_ij| = 2|Im h_ij|`, `2|Re h_ij|` or
   `sqrt2 |Re h_ij -+ Im h_ij|`, according as `eps = 1`, `eps = -1` or `eps = +-i`. In the
   common-unit case, `|v_ij| = 2|t_ij|`.

*Proof.*

* (1) is a direct prime-by-prime comparison.
* (2) At one prime, the matching that pairs the low points with each other and the high points
  with each other contributes
  `pi^{b_1+b_3} conj(pi)^{2e-b_2-b_4} = Gamma_p (pi conj(pi))^{b_3-b_2}`. Each of the other two
  matchings contributes `pi^{b_1+b_2} conj(pi)^{2e-b_3-b_4} = Gamma_p`. Ties between `b_2` and
  `b_3` give `m_p = 0`. The threshold layers `l` with `b_2 < l <= b_3` are exactly the layers
  whose cut is 2|2 of type `kappa(p)`.
* (3) Substitute (1) and (2) into the polynomial identity
  `(z_1 - z_3)(z_2 - z_4) = (z_1 - z_2)(z_3 - z_4) + (z_1 - z_4)(z_2 - z_3)`, then cancel the
  nonzero factor `Gamma`.
* (4) The cross-ratio of four concyclic points is real, so the three products lie on one real
  line through `0`. Ptolemy's theorem for cyclic order gives
  `|v13 v24| X_2 = |v12 v34| X_1 + |v14 v23| X_3`. Together with (3), this forces all three
  products onto the same ray. Gaussian integers on a ray are positive multiples of its primitive
  generator. Coprimality of the `X_k` holds because each prime appears in at most one of them.
* (5) Put `h = x + iy`. Then `h - conj h = 2iy`, `h + conj h = 2x`, and
  `h -+ i conj(h) = (x -+ y)(1 -+ i)`. `[]`

*Relation to the notes.* The common-unit case of (4) is equation (5) of
`ordered_residue_growth.md`.

Items 396-416 extract specific integers from ordered 4- and 5-point relations.

* Some of them are functions of the point of `M_{0,5}` or `M_{0,6}` and its `p`-adic boundary
  contacts (Theorem B below): the positive-pentagon content `F` of item 413 and the
  adjacent-exchange quantities of item 414, which are built from cross-ratios, and the gap `K`
  of item 399, which is the contact of the split {outer pair} | {interior triple}.
* The affine index of item 396 also uses the planar embedding, through triangle areas, and so
  involves the norms.

Parts (3) and (5) of Theorem A, which cover the mixed-unit version and the Gaussian (not merely
modular) identity, are new as far as I could grep. The check
`check_general_ptolemy.py` verifies (1)-(5) and the threshold-layer description exactly. It covers
6258 quadruples in mixed-unit clusters with prime powers, plus 40 genuinely shortest windows on
actual circles.

## 3. Theorem B: what the system contains

**Theorem B.** Work in the setting of Section 1.

**(i) Anchored form.** For all `i, j` in `{0, ..., m}`,

```text
z_i - z_j = -2i z_0 det(P_i, P_j) / (conj(P_i) conj(P_j)).
```

Consequently the Ptolemy relation of any four points is the Plücker relation
`p_ac p_bd = p_ab p_cd + p_ad p_bc` for `p_ij = det(P_i, P_j) = t_ij G_ij` (`t_0j = Y_j`,
`G_0j = 1`). The factors `conj(P_i)` cancel identically.

For a quadruple `Q`, a block `T` behaves as follows:

* if `T` meets `Q` in 3 or 4 rows, it contributes to all three terms;
* if `T` meets `Q` in 2 rows, it contributes to exactly one term;
* if `T` meets `Q` in at most 1 row, it contributes to none.

Dividing by the common part gives the matching form

```text
n_C (t_ac t_bd) = n_A (t_ab t_cd) + n_B (t_ad t_bc),     n_A = prod_{T cap Q in {ab, cd}} n_T, etc.   (3.1)
```

Theorem A(4) is this identity read through `|v_ij| = 2|t_ij|`.

**(ii) Rank and generation.** Over `Q`, the `C(k,4)` relations have Jacobian rank `(k-2)(k-3)/2`.
This is the codimension of the affine cone over `G(2,k)`. On `{p_01 != 0}`, the
`(k-2)(k-3)/2` relations through the fixed pair `{0,1}` generate all the others. For fixed
nonzero residues, the relations are linear in `(G_ij)`. Their solution space is exactly
`{G_ij = (X_i Y_j - X_j Y_i)/t_ij}`, of dimension `m - 1`.

**(iii) Invisibility (form-blindness).** The matching coefficients `n_A, n_B, n_C` involve only
boundary blocks `2 <= |T| <= m-1`. A singleton or the top block meets every quadruple in a set
that is never a 2|2 cut. Conversely, let `n_T` be positive integers and `P_i` primitive integer
vectors with `prod_{T contains i,j} n_T | det(P_i, P_j)`, with no further condition. Then every
relation (3.1) holds, with `t_ij := det/G_ij`. In particular the system is satisfied by Plücker
data that violate all of the following:

* the Gaussian norm condition `|P_i|^2 = prod_{T contains i} n_T`;
* the isotropy condition `X_i = rho_T Y_i mod n_T` with `rho_T^2 = -1`;
* the existence of the singleton blocks.

Explicit examples are in Section 5.

**(iv) Moduli dictionary.** Put `Phi = [P_0 : ... : P_m]` in `M_{0,k}(Q)`. For a prime `p` and a
boundary split `T`, let `l_T(p)` be the `p`-adic tree edge length. It is computed from the Gromov
products `g_ij = v_p(det(P_i, P_j))` by

```text
l_T(p) = max(0, min_{a,b in T; c,d not in T} min(g_ab + g_cd - g_ac - g_bd, g_ab + g_cd - g_ad - g_bc)).
```

Then the following hold.

1. `l_T(p)` is the local intersection number `(Phi.D_T)_p` on Knudsen's smooth model of
   `Mbar_{0,k}` over `Z`.
2. If `p` divides no residual, then `l_T(p) = v_p(n_T)`. In general,
   `|l_T(p) - v_p(n_T)| <= 2 max_{i<j} v_p(t_ij)`.
3. The non-archimedean part of the height of `D_T` therefore satisfies

   ```text
   | sum_p l_T(p) log p - log n_T | <= 2 H_Sigma.                       (3.2)
   ```

**(v) S-unit form.** For a quadruple in cyclic order, (3.1) is

```text
chi + chi' = 1,
chi  = (n_A/n_C) (t_ab t_cd)/(t_ac t_bd),
chi' = (n_B/n_C) (t_ad t_bc)/(t_ac t_bd).
```

That is, the cross-ratio unit equation of `M_{0,k}`. The S-unit parts are the matching
coefficients, and the non-unit parts are bounded by the residues.

*Proof.*

* **(i)** `z_i - z_j = z_0 (P_i conj(P_j) - P_j conj(P_i)) / (conj(P_i) conj(P_j))`, and
  `P_i conj(P_j) - conj(P_i) P_j = 2i Im(P_i conj(P_j)) = -2i det(P_i, P_j)`. Insert this into the
  polynomial Ptolemy identity and multiply through by `prod conj(P)/(-2i z_0)^2`. The block count
  in a single `G_ij` gives the three cases listed in (i).
* **(ii)** This is standard for Grassmannians. The exact check `check_system_rank.py` verifies the
  ranks for `k = 4..9` and the linear solution space for `m = 3..7`.
* **(iii)** The first claim is the case analysis in (i). The second holds because (3.1) is an
  identity for arbitrary vectors once the common factor is removed.
* **(iv)** The forgetful map `chi_Q: Mbar_{0,k} -> Mbar_{0,4} = P^1` satisfies
  `chi_Q^*[0] = sum_{T inducing ab|cd} D_T`, and similarly at `1` and at `infinity`. So for the
  `Z_p`-point extending `Phi`,
  `v_p(chi_Q(Phi)) = sum_{T: ab|cd} (Phi.D_T)_p - sum_{T: ac|bd} (Phi.D_T)_p`.
  The `T` with `(Phi.D_T)_p > 0` are pairwise compatible, because the reduction lies in a nonempty
  stratum. A compatible weighted split system is determined by its quartet valuations, through
  the min-formula above (four-point condition). The Gromov products `v_p(det)` of primitive
  vectors are the Bruhat–Tits Gromov products from the standard vertex. They give a compatible
  weighted split system with the same quartet valuations. Hence `(Phi.D_T)_p = l_T(p)`.
  - With zero residual valuations, `g_ij = sum_{T contains i,j} v_p(n_T)`. Here the blocks at
    `p` form two nested chains on disjoint sets, hence are compatible. This is exactly the
    Gromov-product system of the tree with edges `T` of length `v_p(n_T)`, seen from a base
    vertex on the `0`-side of every edge, so `l_T(p) = v_p(n_T)`. The singleton and top blocks
    are pendant edges and are invisible to every quartet.
  - Adding residual valuations `delta_ij = v_p(t_ij) >= 0` changes each quartet form by at most
    `2 max delta`. The min-formula is 1-Lipschitz in sup norm, which gives the bound in (2).
  - Summing `2 max_{i<j} v_p(t_ij) log p <= 2 sum_{i<j} v_p(t_ij) log p` over `p` gives (3.2).
* **(v)** Divide (3.1) by `n_C t_ac t_bd`. `[]`

*Checks.*

* `check_m0k_tree.py` uses 246 actual common-unit clusters with 5 to 8 points, arbitrary
  exponents and nested blocks, and examines 3971 primes, 3352 of which divide some residual. At
  every prime it verifies:
  1. the positive edges are compatible;
  2. the tree reproduces every quartet cross-ratio valuation;
  3. `|l_T(p) - v_p(n_T)| <= 2 max v_p(t)`, with equality whenever no residual is divisible by
     `p`.
* `check_general_ptolemy.py` (G7) checks
  `det(P_i, P_j) = t_ij prod_{T contains i,j} n_T` and the anchored matching identity on 1363
  quadruples.

## 4. The quadruple system for k = 5, 6, 7, and the overdetermination question

**Explicit form.** Take the anchor `0` and rows `1..m`. For `Q = {0, i, j, l}` with `i < j < l`,
let `x` range over the remaining rows. Relation (3.1) reads

```text
Y_i t_jl N_jl  -  Y_j t_il N_il  +  Y_l t_ij N_ij  =  0,
N_ab := prod_{T cap {i,j,l} = {a,b}} n_T,                              (4.1)
```

up to orientation signs. For nonanchor quadruples it reads
`n_C t_ac t_bd = n_A t_ab t_cd + n_B t_ad t_bc`, with each `n_A, n_B, n_C` a product of
`2^(m-3)` boundary blocks.

* For `k = 5` there are 5 relations of rank 3. The 4 anchored ones have
  `N_ab = n_ab n_{ab x}` (two blocks). The fifth is
  `n_13 n_24 t_13 t_24 = n_12 n_34 t_12 t_34 + n_14 n_23 t_14 t_23`.
* For `k = 6` there are 15 relations of rank 6.
* For `k = 7` there are 35 relations of rank 10.

**Counts.**

| k | visible blocks `2^(k-1)-k-1` | invisible blocks | residues `C(k,2)` | relations | rank |
|---|---|---|---|---|---|
| 5 | 10 | 5 | 10 | 5 | 3 |
| 6 | 25 | 6 | 15 | 15 | 6 |
| 7 | 56 | 7 | 21 | 35 | 10 |
| 8 | 119 | 8 | 28 | 70 | 15 |

**Dependencies.** The only dependencies are the Plücker (Grassmann) syzygies. With base pair
`{0,1}`, every relation follows from the `(k-2)(k-3)/2` relations

```text
Y_1 t_ij G_ij = Y_i t_1j G_1j - Y_j t_1i G_1i      (2 <= i < j <= m),   (4.2)
```

on `{Y_1 != 0}`. After the common factors are cancelled, these are exactly the "polynomial
identities among block norms with small coefficients" that elimination forces. They are
satisfiable for every `k`. The solution variety in the `G`'s is rational, parametrized by
`(X_i)` (Theorem B(ii)). So no `k` makes the system algebraically overdetermined, and no
contradiction can come from elimination alone.

This agrees with the barrier in `plucker_compatibility_route.md` §1. For every fixed bracket
polynomial, the average prescribed vanishing order is at most that of a product of minors.

**Where the overdetermination actually is.** It is arithmetic, and it is controlled by one
intersection number. By Theorem B(iv), a balanced configuration is a point of `M_{0,k}(Q)` with
`h_{D_T} = w(1 + o(1))` for every boundary divisor, and with archimedean proximity
`O(eta w + H)`. Its canonical height is therefore `w (K.beta) + o(w)`, where `beta` is the
balanced curve class. Using the canonical class (Keel–McKernan, Pandharipande),
`K = sum_S (s(k-s)/(k-1) - 2) D_S`:

| k | `(K+B).beta` | `K.beta` | expected dim of rational curves of class `beta` | integer heuristic exponent per shift class |
|---|---|---|---|---|
| 5 | 5 | -5 | 4 | 6 |
| 6 | 17 | -8 | 8 | 9 |
| 7 | 49 | -7 | 8 | 8 |
| 8 | 129 | **+10** | -8 | -9 |
| 9 | 321 | 75 | -72 | -74 |

The closed forms are `(K+B).beta = 2^(k-3)(k-4) + 1`, `K.beta = 2^(k-3)(k-8) + k + 2`, and
heuristic exponent `= 1 - K.beta`. The last column is the heuristic of `check_system_rank.py`:
`exp(e w)` expected configurations per shift class. `check_k_class.py` checks the canonical class
in Kapranov's basis for `k = 5..10`, and checks all the closed forms.

So the "angular-only" system becomes overdetermined exactly at `k = 8`. That is the first `k` at
which `K.beta > 0`, and Vojta's conjecture then predicts a height gap (Section 6).

## 5. Theorem C: bounded-residue balanced solutions for k = 5 and 6

### 5.1 An explicit k = 5 family

Put `a = (1, 2, 3, 4)`, `tau = (tau_0, ..., tau_4) = (4, 6, 7, 0, 2)`, and for `i = 1..4`

```text
X_i(s) = a_i (s - tau_0) prod_{l != i} (s - tau_l),     P_0 = (1,0),  P_i = (X_i(s), 1).
```

The blocks are:

* the top block `n_[4] = s - tau_0`;
* the triples `n_{[4]\{l}} = s - tau_l`;
* the pairs `n_ij = (a_i - a_j) s - (a_i tau_j - a_j tau_i)`.

**Theorem C1.** The following hold.

1. `det(P_0, P_i) = 1` and `det(P_i, P_j) = prod_{T contains i,j} n_T(s)`, as polynomial
   identities, so all `t_ij = 1`.
2. All five quadruple relations hold in the matching form `n_C = n_A + n_B`, as polynomial
   identities. The `n`'s are quadratic.
3. The 11 roots are distinct.
4. Take `s = 52346140003 mod 85900454640`. Every block has constant `{2,3}`-part on this class,
   from `{1, 2, 3, 4, 12}`. Move these parts into the residues; then `|t'_ij| <= 144` and
   `|Y_i| = 1`. The reduced blocks are pairwise coprime, since every common prime divides a
   resultant and all such primes `>= 5` are excluded. Also `|log n'_T - log s| < log 20`.
5. At every block prime, all quartet cross-ratio valuations equal the Boolean ones.

So for `w = log s -> infinity`, this is a balanced, coprime-block, bounded-residue solution of the
whole `k = 5` Ptolemy system (`H = O(1)`).

*Proof.* (1) holds because
`X_i - X_j = (s - tau_0) prod_{l != i,j} (s - tau_l) [a_i (s - tau_j) - a_j (s - tau_i)]`. The rest
is checked by `k5_family.py`, as polynomial identities and at 399 values of `s`. `[]`

The family is not a circle configuration: `X_i^2 + 1` has no reason to factor along the blocks.
It is a point of the form-blind relaxation (Theorem B(iii)). Geometrically, after the top factor is dropped, it is a rational
plane cubic through the four triple points of the `A_3` arrangement. That is an
anticanonical rational curve on the del Pezzo surface `Mbar_{0,5}`.

### 5.2 Anchoring lemma, and k = 6 from item 109

**Lemma 5.2 (anchoring).** Let `w_0, ..., w_m` be primitive integer vectors representing distinct
points of `P^1(Q)`, with `w_0 = (1,0)` (move it there by `SL_2(Z)` first), and write
`w_i = (x_i, y_i)`. Put `D = lcm_j |y_j|` and let `P_i` be the primitive integer vector
proportional to `(D x_i, y_i)`. Then:

* the `P_i` represent the same point of `M_{0,k}`, and `P_0 = (1,0)`;
* `det(P_0, P_j) = +-1`;
* for every prime `p`, with `l_T(p)` the edge lengths of Theorem B(iv),

  ```text
  v_p(det(P_i, P_j)) = delta_p + sum_{T contains i,j} l_T(p)      (i, j >= 1),
  ```

  with one constant `delta_p >= 0` independent of `i, j`.

So after anchoring, every configuration is in exact anchored Boolean form. The constant
`prod p^{delta_p}` is an extra top block, which is invisible. The residues `t_ij` come only from
the non-Boolean part of the `p`-adic trees.

*Proof.* Put `r_p = v_p(D) = max_j v_p(y_j)`. The vertices on the ray from the standard vertex
`Z_p^2` towards the end `(1,0)` are `L_r = Z_p + p^r Z_p`, for `r >= 0`. The geodesic from end `0`
to end `j` meets this ray at distance `(0|j) = v_p(det(w_0, w_j)) = v_p(y_j)` from the standard
vertex. So `L_{r_p}` lies on every geodesic `(0, j)`, that is, on the ray from the branch point
`b_p` towards `0`, at some distance `delta_p >= 0` beyond `b_p`.

The map `diag(D, 1)` sends `L_{r_p}` to the standard vertex, up to homothety, at every `p`. For
primitive vectors, `v_p(det)` is the Gromov product from the standard vertex. From a vertex on
the `0`-ray, this Gromov product is `0` for pairs `(0, j)` and `delta_p + sum_{T contains i,j}
l_T(p)` for `i, j >= 1`. Finally, `gcd(D x_j, y_j) = |y_j|` because `y_j | D` and
`gcd(x_j, y_j) = 1`, so `det(P_0, P_j) = y_j/|y_j| = +-1`. `[]`

`check_anchoring.py` verifies the explicit choice exactly, on 300 random configurations and 7599
primes. It checks that the cross-ratios are unchanged, that `det(P_0, P_j) = +-1`, and that the
Gromov products equal `delta_p + sum l_T(p)`.

**Corollary C2 (k = 6).** Item 109 (`disjoint_mixed_elliptic_integer_contact_profile.md`) supplies
infinitely many rational six-point configurations. Their 25 boundary contacts `d_C` are pairwise
coprime, with `log d_C = alpha N^2 + O(N)`. Their determinant remainders `b_ij` are bounded and
supported on a fixed finite set `Sigma`. Their anchors are `v_infinity = (1,0)`, `v_0`, `v_1`.
By Lemma 5.2, applied with anchor `infinity`, they become anchored Plücker data with:

* `|Y_j|` and `|t_ij|` bounded;
* balanced, pairwise coprime boundary blocks;
* `w = alpha N^2 -> infinity`.

So the `k = 6` Ptolemy system also has bounded-residue balanced solutions. Here `D` is
`prod_{C containing infinity} d_C` up to `Sigma`-factors. The residues are bounded, because the
`p`-adic trees are exactly Boolean off `Sigma` and the `Sigma`-parts are bounded by the `b_ij`.

No extra top block appears (`delta_p = 0`). The reason is that the standard vertex is the median
of the three anchors `infinity, 0, 1`, and that median lies on the `infinity`-ray on the near side
of `b_p`. A top block of any prescribed size `n` (for instance `n ~ e^w`) can be added by
`diag(n, 1)`. Singleton blocks never enter.

### 5.3 Consequence

**Corollary C3.** For `k = 5` and `k = 6` there is no function `f(w) -> infinity` such that every
balanced, coprime-block, integral solution of the full quadruple system, in cyclic order,
satisfies `H >= f(w)`. Hence a lower bound `H >= c w`, or any superlogarithmic bound, at `k <= 6`
must use information outside the Ptolemy relations. Examples are the Gaussian norm/isotropy
conditions or the singleton blocks.

### 5.4 k = 7 (numerical exploration only; not used in any proof)

* For `k = 6`, Newton's method in a pole chart with the 10 triple points fixed easily finds
  nondegenerate balanced rational curves (`rchart5b.py`). At one rational parameter point it
  found at least 5 distinct ones, all real (`count5.py`). None was rational (`ratcheck5.py`).
  So genus-0 balanced curves exist over `R` at `k = 6`, but the degree of the parameter map is
  at least 5. Rational ones were not found. Item 109's elliptic curves cover `k = 6`.
* For `k = 7`, the expected dimension is `+8`. Every chart I tried requires solving for at least
  30 of the 56 special points. The X-chart (`balanced_newton.py`), the homotopy tracker
  (`homotopy.py`) and the pole charts (`rchart6.py`, `rchart6b.py`) all converged only to
  degenerate solutions (merged special points or a degree drop). The pole-chart runs comprised
  about 155 Newton trials with 4 seeds, plus smaller X-chart and homotopy batches. Existence of balanced rational or elliptic curves in `Mbar_{0,7}` therefore remains
  open, consistent with item 112 of the notes.

  If such a curve exists over `Q`, Lemma 5.2 extends Corollary C3 to `k = 7`. By Theorem D, `8`
  would then be the exact threshold.

## 6. Theorem D: Vojta's conjecture on Mbar_{0,8} implies the uniform theorem

**Hypothesis `V_8`.** This is Vojta's Main Conjecture (P. Vojta, *Diophantine Approximations and Value
Distribution Theory*, LNM 1239, Main Conjecture 3.4.3; see also Bombieri–Gubler, *Heights in
Diophantine Geometry*, Ch. 14) for the smooth projective variety `X = Mbar_{0,8}` over `Q`, in
the case `D = 0`. Fix an ample divisor `A`; we take `A = K_X + B`, where `B` is the boundary,
which is ample by Keel–McKernan. The statement is: for every `eps > 0` there are a proper
Zariski-closed `Z_eps` in `X` and a constant `C_eps` such that

```text
h_{K_X}(P) <= eps h_A(P) + C_eps        for all P in X(Q) \ Z_eps .
```

With `D = 0` no set of places `S` enters. `Z_eps` depends only on `(X, A, eps)`. It does not
depend on `R`, on the primes of `N`, or on any template.

The ampleness of `K + B` is not essential. `Pic(Mbar_{0,8})` is spanned over `Q` by the
boundary divisors (Keel). Lemma 6.2 below bounds every `h_{D_T}` above and below. Hence every
ample `A` satisfies `h_A <= c_A (w + H_Sigma) + O(1)` at the points used. That is all Step 4
needs, after shrinking `eps`.

**Theorem D.** Assume `V_8`. Then for every `C > 0` there is `M(C) < infinity` such that every arc
of length at most `C sqrt(R)` on `x^2 + y^2 = R^2` (`R^2` in `Z`, any position) contains at most
`M(C)` lattice points.

### 6.1 Two lemmas

**Lemma 6.1 (exceptional set).** Let `Z` be a proper Zariski-closed subset of `Mbar_{0,k}`. There
is `c(Z)` such that, for every finite set `Sigma` of `P^1(Q)` with `|Sigma| = M`, at most
`c(Z) M^{k-1}` ordered `k`-tuples of distinct points of `Sigma` have their configuration
in `Z`.

*Proof.* Induct on `k`. For `k = 3`, `M_{0,3}` is a point, so `Z` meets it in the empty set.
Let `Z° = Z cap M_{0,k}` and let `pi: M_{0,k} -> M_{0,k-1}` forget the last point. Its fibres are
irreducible curves (`P^1` minus `k-1` points). Let `Z_1` be the closure of
`{b : pi^{-1}(b) subset Z°}`. The preimage of that set lies in `Z°`, which is not dense in
`M_{0,k}`. So `Z_1` is a proper closed subset.

Off `Z_1`, the fibres of `pi|Z°` are finite. Their cardinality is bounded by some `d(Z)`, since
the number of geometric points in the finite fibres of a morphism of finite type is bounded. So
the number of tuples in `Z°` is at most

```text
#{(k-1)-tuples in Z_1} * M  +  d(Z) * M^{k-1}  <=  (c(Z_1) + d(Z)) M^{k-1}
```

by induction. `[]`

**Lemma 6.2 (heights of a balanced Boolean point).** Let `Phi` in `M_{0,k}(Q)` come from an
anchored tuple as in Section 1, with residual height `H_Sigma`, and suppose
`|log n_T - w| <= eta w` for every boundary `T`. Put `m = k-1`. Then for every boundary
divisor `D_T`, with fixed Weil functions (model metrics at the finite places),

```text
log n_T - 2 H_Sigma - c_0  <=  h_{D_T}(Phi)  <=  log n_T + 3 H_Sigma + 2^(m-2) eta w + c_0 .
```

*Proof.*

* The finite part is (3.2).
* At the archimedean place, `lambda_{D_T,inf}` is bounded below because `D_T` is effective.
* For the upper bound, choose `a, b` in `T` and `c, d` not in `T`. Then `D_T <= chi^*[0]` for
  the cross-ratio `chi = det(a,b) det(c,d) / (det(a,c) det(b,d))`, whose zero is the boundary
  point `ab|cd` of `Mbar_{0,4}`. Hence `lambda_{D_T,inf} <= log^+ (1/|chi(Phi)|) + O(1)`.
* Now `log|chi| = sum_T eps_T log n_T + log|t_ab t_cd / (t_ac t_bd)|`. Exactly `2^(m-3)` of the
  blocks have `eps_T = +1` and `2^(m-3)` have `eps_T = -1`, with every block within `eta w` of
  `w`. So the first sum is at most `2^(m-2) eta w` in absolute value.
* The residual ratio has `|log| <= H_Sigma`, because all residues are integers with `|t| >= 1`.
  `[]`

### 6.2 Proof of Theorem D

**Step 0 (reductions, as in item 532).**

* Subdivide the arc into `s = ceil(C/C_0)` pieces with `C_0 < sqrt 2`, say `C_0 = 1`.
* Retain a largest unit class in a largest piece.
* Divide by its Gaussian gcd.

The resulting class has `M >= M_orig/(4s)` points, and its radius and arc constant do not
increase. Put `k = 8`, `m = 7`, `w = W/128` and `lambda = log(2/C_0) > 0`. Then
`W >= 4(M-1) lambda`, so `w >= (M-1) lambda / 32`.

**Step 1 (many good tuples).** For a uniformly random ordered 8-tuple, item 532 gives
`E[A] <= a = 2*255/(M-7)` and `E[B] <= b_0 = 56 W/(4(M-1))`. By Markov's inequality, at least
half of all ordered 8-tuples satisfy `A <= 4a` and `B <= 4 b_0`. For such a tuple, Fourier
inversion on the sign cube gives

```text
|log n_T - w| <= eta w   for every nonempty T in [7],    eta = 2 sqrt2 * 255 / sqrt(M - 7),
H_Sigma <= B/2 <= 28 W/(M-1) = 3584 w/(M-1).
```

The second line uses item 532 (9) with `C_0 < 2`, which gives `log|t| <= delta/2`.

**Step 2 (avoid `Z_eps`).** Fix `eps = 1/43` and let `Z = Z_{1/43}`. By Lemma 6.1, if
`M(M-1)...(M-7)/2 > c(Z) M^7`, which holds for `M >= M_Z`, some tuple from Step 1 has
configuration `Phi` outside `Z`. The configuration is the point of `M_{0,8}(Q)` given by the
eight circle points, or equivalently by the vectors `P_i` (Theorem B(i)).

**Step 3 (heights).** For `k = 8`, `7K = sum_S (s(8-s) - 14) D_S`, and the coefficients
`s(8-s)/7 - 2` are `-2/7`, `1/7` and `2/7` for `s = 2, 3, 4` (28, 56 and 35 divisors). Lemma 6.2
gives

```text
h_K(Phi) >= (8+10)((1-eta)w - 2 H_Sigma) - 8((1 + 33 eta) w + 3 H_Sigma) - c
          = 10 w - 282 eta w - 60 H_Sigma - c,
h_A(Phi) <= 129 ((1 + 33 eta) w + 3 H_Sigma) + c = 129 w + 4257 eta w + 387 H_Sigma + c,
```

because `A = K + B` has positive coefficients `5/7, 8/7, 9/7`, which sum to 129 over the
boundary.

**Step 4 (apply `V_8`).** From `h_K(Phi) <= h_A(Phi)/43 + C_{1/43}`:

```text
10w - 282 eta w - 60 H_Sigma <= 3w + 99 eta w + 9 H_Sigma + C'
  =>  (7 - 381 eta) w <= 69 H_Sigma + C'.
```

If `eta <= 1/381` this gives `H_Sigma >= (6w - C')/69`.

**Step 5 (conclusion).** Combine Step 4 with Step 1:

```text
(6w - C')/69 <= 3584 w/(M-1).
```

This is impossible as soon as all of the following hold:

* `M - 1 > 247296/5`;
* `w >= C'`;
* `eta <= 1/381`, which holds once `M >= 7 + (2 sqrt2 * 255 * 381)^2`, about `7.6 * 10^10`;
* `M >= M_Z`.

The condition `w >= C'` holds once `M >= 1 + 32 C'/lambda`. Hence the class has
`M < max(7.6 * 10^10, M_Z, 1 + 32 C'/lambda)`, a bound independent of `R` and of the arc
position. So `M_orig <= 4 s` times this bound. `[]`

*Remarks.*

1. The bound is ineffective, through `Z_eps` and `C_eps`, as Vojta's conjecture is.
2. Only `V_8` is used. The argument needs `K.beta > 0`, so it fails for `Mbar_{0,7}`, where the
   corresponding `h_K` is about `-7w`. Using `V_k` for larger `k` also works, but no smaller `k`
   does. Indeed, for `k <= 7` every coefficient `s(k-s)/(k-1) - 2` is `<= 0`
   (`K = -(1/3) D_2` on `Mbar_{0,7}`, `K = -(2/5) D_2 - (1/5) D_3` on `Mbar_{0,6}`). So
   `-K` is effective and supported on the boundary, and `V_k` (`D = 0`) holds trivially off the
   boundary. `k = 8` is the first case in which Vojta's conjecture says anything about these
   points.
3. The proof never uses the Gaussian norm structure beyond the extraction. Under `V_8` the same
   conclusion holds for any form-blind analogue that has the same Boolean extraction.
4. The other conditional theorem (`growth/conditional.md`, Theorem 5.1) uses `V` with `D = 0` on
   a resolution of the blow-up of the 7-point circle variety `X_7` along its diagonal. It needs 7
   points and uses the quadric (norm) structure and the archimedean proximity directly. Theorem D
   uses only Ptolemy/cross-ratio data and non-archimedean contacts. The two are logically
   independent instances of the same conjecture.

### 6.3 What Theorem D says about the route

At `k = 8`, Vojta's inequality for these points is literally the following statement:
balanced 2|2-coefficient S-unit systems of the form (3.1) cannot have bounded residues. For
`k <= 7` the same inequality holds automatically, because `h_K < 0`. So no height gap is
predicted there, and for `k = 5, 6` bounded-residue solutions exist (Theorem C). The Ptolemy/S-unit
route therefore has exactly one consistent target, the inequality at `k >= 8`.

## 7. The missing lemma and why the method stalls

**(L_8), the missing lemma.** There are `c > 0`, `eta_0 > 0` and a proper Zariski-closed `Z` in
`Mbar_{0,8}` such that the following holds. Let `Phi` in `M_{0,8}(Q) \ Z` be given by anchored
integer vectors with `det(P_i, P_j) = t_ij prod_{T contains i,j} n_T`. Here the `n_T` lie in
two nested chains at each prime, as in item 532, and satisfy `|log n_T - w| <= eta_0 w`. Then

```text
sum log|Y_i| + sum log|t_ij|  >=  c w - O(1).
```

* `(L_8)` implies the uniform theorem, by Steps 0-5 of Section 6.2 with Step 4 replaced by
  `(L_8)`.
* `V_8` implies `(L_8)`.
* `(L_8)` is a statement about rank-2 integer matrices and divisibility of minors. It involves
  no Gaussian integers and no circle.
* `(L_k)` is false for `k = 5, 6` (Theorem C). For `k = 7` it is open; it would be false if a
  balanced rational curve over `Q` exists, which the expected dimension `+8` suggests.

**Why every unconditional tool I know stalls at (L_8).**

1. **Elimination.** The relations generate the Plücker ideal, and its variety is rational. By
   `plucker_compatibility_route.md` §1, every fixed bracket polynomial has average prescribed
   vanishing order at most that of a product of minors. So no polynomial identity gives a
   positive multiple of `w`.
2. **Single S-unit equations and `abc`.** Each relation is `n_A a + n_B b = n_C c`, with
   `n_A n_B n_C` essentially squarefree on the core primes. So `rad` is comparable to the
   product, and `abc`, `n`-term `abc` and `abc` over `Q(i)` are vacuous; `conditional.md` gives a
   proposition to this effect. The gap appears only when all `C(8,4) = 70` equations are taken
   jointly. That joint statement is the inequality `h_K <= eps h` on the 5-dimensional
   `Mbar_{0,8}`, not a curve statement.
3. **Subspace / Ru–Vojta / Evertse.** In S-unit form, `S` is the set of primes of `N`. It is
   unbounded, and the number of primes is of order `log R/loglog R`. Counting bounds such as
   `exp(c |S|)` are useless, and the exceptional sets depend on `S`. This is the documented
   failure in `claimed_complete_proof_salvage_audit.md`. In the `D = 0` form no `S` appears, but
   the inequality involves the global heights `h_{D_T}`. These add up the contacts over all,
   unboundedly many, primes. That is a GCD-type statement on an iterated blow-up of `P^5`
   (Kapranov's model), and Subspace-type methods with a fixed `S` do not control it. `V_8` is
   the first nontrivial instance (Remark 2 of Section 6.2).
4. **Exact arithmetic of the circle.** The isotropy of residues (`rho_T^2 = -1`), the singleton
   blocks and the Gaussian norms are exactly what the Ptolemy system forgets (Theorem B(iii)).
   They are what distinguish the circle from the `k = 5, 6` examples. Any proof at `k <= 7` must
   use them. That is a different route: the full Gaussian model, items 515-533.

## 8. Files and exact checks

All scripts are in `growth/ptolemy/` and are run with `python3` from that directory. Exact
arithmetic uses Python integers and `Fraction`. Floating point is used only in the scripts marked
NUMERICAL.

| Script | Content | Result |
|---|---|---|
| `check_general_ptolemy.py` | Theorem A (1)-(5) and threshold layers, on mixed-unit clusters with prime powers; 40 shortest windows on actual circles; anchored identity (G7) | 6258 quadruples and 1363 anchored quadruples, all exact |
| `check_anchored_identity.py` | Theorem B(i): `z_i - z_j = -2i z_0 det(P_i,P_j)/(conj P_i conj P_j)` | 2500 ordered pairs, exact |
| `check_system_rank.py` | Theorem B(ii): ranks `(k-2)(k-3)/2` for `k = 4..9`; linear solution space in `G`; heuristic table | passed |
| `check_m0k_tree.py` | Theorem B(iv): tree edge `l_T(p) = v_p(n_T)` (equal when `p` divides no residual), Lipschitz bound, cross-ratio reproduction, nested blocks | 246 clusters, 3971 primes, passed |
| `check_k_class.py` | `K = sum (s(n-s)/(n-1)-2) D_S` in Kapranov's basis for `n = 5..10`; closed forms for `K.beta`, `(K+B).beta`; heuristic `= 1 - K.beta` | passed |
| `k5_family.py` | Theorem C1: factorization, all 5 matching relations as polynomial identities, coprimality class, balance, Boolean trees | passed |
| `check_anchoring.py` | Lemma 5.2 with `D = lcm |y_j|` | 300 configurations, 7599 primes |
| `rchart5b.py`, `count5.py`, `ratcheck5.py` | NUMERICAL: balanced genus-0 curves in `Mbar_{0,6}` exist over `R` (at least 5 per parameter point); none rational | evidence only |
| `balanced_newton.py`, `homotopy.py`, `rchart6.py`, `rchart6b.py`, `prox_newton.py` | NUMERICAL: `k = 7` searches; only degenerate limits found | evidence only |
| `poly_family_newton.py`, `k5_family_search.py` | from the interrupted earlier attempt at this route; superseded | - |
