# Round 18: the lifting test for (L_8)

## 0. Summary

**Task** (round-17 judge, §7.3 task 2). Decide whether the `(L_6)`-refuting families `E_t ⊂ Mbar_{0,6}`
(`threshold/L8.md` Thm B), or any bounded-residue family at `k = 6, 7`, lift to balanced anchored
8-configurations with bounded residues. Do the exact dimension count, then test explicit members.

**Verdict.** The `E_t` mechanism does **not** lift, already at the first step, at the curve level. The
counts are overdetermined at every level. Explicit members show no lift. `(L_8)` itself is untouched.

1. **Theorem 4.2 (new, unconditional, exact computation plus an openness argument).** For `t` in a
   nonempty Zariski-open `U_0` of the parameter space (containing `t_0 = (2,3,-5,1)`, item 109's
   `(1/3,3,2,5/3)` and 7 further explicit points), the curves of class `beta_7` in `Mbar_{0,7}` that
   project onto `E_t` are **exactly the pencil lifts**
   `C'_lam = {f_a(a) = f_b(b) = f_c(c) = f_lam(d)}`, of geometric genus 9 (5 for at most two `lam`).
   * In particular no class-`beta_7` curve of genus `<= 4` lies over any such `E_t`, so the positive-rank
     genus-one mechanism of Theorem B does not extend to a seventh point.
   * The proof reduces the question to a rank condition on a `4 x 3` matrix of sections over `Pic^8(E_t)`
     (Lemma 4.1, eq. (4.1)), and decides it exactly over `Q`.
   * The linear count predicts `-1`; the truth is a 1-parameter pencil, forced by the pencil identity
     `P_lam - z^2 Q_lam = lam prod(z - x_j)`, and nothing else.
2. **Corollary 4.3.** A class-`beta_8` curve of genus `<= 4` over `E_t` must have both 7-point projections
   birational, of class `2 beta_7` and degree 4 over `E_t`. That case is not decided; its expected
   dimension is `-18` per projection and `-42` in all. The pencil lifts to `beta_8` exist, with genus 33.
3. **Counts (Section 3), all exact bookkeeping, conclusions heuristic.**

   | lift | expected dimension or exponent |
   |---|---|
   | one point, genus one, class `beta_7` | `-9` per `E_t` |
   | one point, genus one, class `2 beta_7` | `-18` per `E_t` |
   | two points, genus one, class `beta_8` | `-42` per `E_t` (`-39` over the 3-parameter family) |
   | integer 7-lift of one configuration | `e^{-9w}` |
   | integer 8-lift of one configuration | `e^{-42w}` |

   * The `28 - 15 = 13` new minors `det(P_i, P_new)` refine the 25 old blocks into 100 and add 19 new
     ones (`119 = 4 * 25 + 19`). For each new point they impose one congruence per old contact, modulo a
     divisor of size `≈ d_A^{1/2}` (total modulus `e^{25w}` at the 7-lift scale `w`), on a box of
     `e^{16w}` pairs. The 19 new blocks are free cofactors. The system is over-determined, not
     under-determined.
   * The curve and integer counts agree exactly. The bookkeeping reproduces the `k = 7` and `k = 8`
     exponents of `ptolemy.md` §4 and `L7.md` §8.
   * The counts alone decide nothing: on `E_t` itself the same per-point count predicts `e^{-7W}` while
     `Phi_N` exists for every `N` (measured in Section 5).
4. **Explicit members (Section 5, Data).** For `Phi_N(t_0)`, `N = 3..8` (`w` up to 75), every
   configuration has between 6 and 16 contact integers with no divisor in `[d^{1/4}, d^{3/4}]`. So no
   balanced 7-lift exists (`eta <= 1/5`, `Sigma`-supported residues), hence no 8-lift.
   * On the remaining contacts the lift lattices are random: maximal skew `0.05 w` at `N = 8`, against
     the `>= 4.5 w` a lift needs.
   * The movers of `Phi_N` themselves show the structural skew `7 w`.
5. **Meaning.**
   * A positive lift would have refuted `(L_8)` only if its configurations were Zariski dense in
     `Mbar_{0,8}`. Otherwise it would have produced the first explicit component of `Z`.
   * The obstruction proved here is unconditional, but it lives at the curve (function-field) level. It
     is **not** a piece of `(L_8)` for points. The point-level statement is a non-skewness statement for
     lattices built from elliptic divisibility sequences, as hard as `(L_8)`.
   * Nothing changes for the growth rate or the uniform bound.

## 1. The object: which formulation of (L_8)

Conventions. *Theorem/Lemma/Proposition* = proved here with all quantifiers (computer steps are exact
rational or exact finite-field computations, reproducible by the scripts of Section 7). *Count* = an
expected-dimension or expected-number computation; the arithmetic of the count is exact, but the
conclusion "expected dimension `< 0` hence nothing exists" is a heuristic and is never used as a proof.
*Data* = an exact finite computation on explicit configurations, used as evidence only.

**(L_8), as used here** (verbatim normalisation of `growth/ptolemy.md` §7 and `threshold/L8.md` §1).
There are `c > 0`, `eta_0 > 0` and a proper Zariski-closed `Z ⊂ Mbar_{0,8}` such that: if
`P_0 = (1,0), P_1, ..., P_7` are primitive integer vectors (anchor `0`), `Phi = [P_0 : ... : P_7]` lies in
`M_{0,8}(Q) \ Z`, and

```text
det(P_i, P_j) = t_ij * prod_{T ⊆ [7], T ∋ i,j} n_T      (t_ij ∈ Z \ {0}),     Y_j := det(P_0, P_j) = t_0j,
```

where at every prime the blocks `n_T` with `p | n_T` form two nested chains on disjoint row sets (item 532)
and `|log n_T - w| <= eta_0 w` for every boundary block (`2 <= |T| <= 6`), then
`H_Sigma := sum_j log|Y_j| + sum_{1<=i<j} log|t_ij| >= c w - O(1)`.

This is the statement that Steps 0-5 of `ptolemy.md` §6.2 need (Step 4 replaced by `(L_8)`): Step 1
produces exactly such tuples (balance `eta = 2 sqrt2 * 255/sqrt(M-7)` for all nonempty `T`, including the
invisible singletons and top block, which `(L_8)` does not use), Step 2 avoids `Z` by Lemma 6.1, and
Step 5 compares with `H_Sigma <= 3584 w/(M-1)`. Nothing below changes this formulation; the lifting test
is a test of whether the known counterexamples one level down propagate to it.

**Negation used.** `(L_8)` is false iff for every `(c, eta_0, Z, C_0)` there is such a configuration off
`Z` with `H_Sigma < c w - C_0`. It suffices to have, for every proper closed `Z`, configurations off `Z`
with `eta -> 0`, `w -> infinity` and `H_Sigma = o(w)` (`L8.md` §1).

**The lifting question.** Let `Phi_6` be a balanced anchored 6-configuration with bounded residues (from
`L8.md` Thm B, Section 2 below). An *8-lift* of `Phi_6` is an anchored 8-configuration as in `(L_8)` whose
sub-configuration on six of its labels (containing the anchor) is `Phi_6`, i.e. the same point of
`M_{0,6}(Q)`. A *7-lift* is defined in the same way with 7 labels.

**Lemma 1.1 (sub-configurations).** Let `Phi^+` be an anchored `k`-configuration as above with data
`(n_T, t_ij)`, scale `w`, balance `eta`, and let `K ⊆ {0,...,k-1}` contain the anchor `0`, `|K| = k'`.
Then the sub-configuration on `K` is anchored with the same vectors `P_i` (`i ∈ K`), the same residues
`t_ij` (`i,j ∈ K`), and blocks

```text
n'_S = prod_{T : T ∩ K = S} n_T          (S ⊆ K \ {0}, S ≠ ∅),
```

each boundary block `n'_S` (`2 <= |S| <= k'-2`) being a product of exactly `2^{k-k'}` boundary blocks of
`Phi^+`. Hence it is balanced at scale `2^{k-k'} w` with the same `eta`, satisfies the nested-chain
condition, and `H_Sigma(sub) <= H_Sigma(Phi^+)`.

*Proof.* For `i, j ∈ K`, `prod_{T ∋ i,j} n_T = prod_{S ∋ i,j} prod_{T ∩ K = S} n_T`. A boundary `S` has
`2 <= |S| <= k'-2`, so every `T` with `T ∩ K = S` satisfies `2 <= |T| <= k-2`; there are `2^{k-k'}` of
them (any subset of the forgotten labels). `T ⊆ T'` implies `T ∩ K ⊆ T' ∩ K`, so nested chains restrict to
nested chains, on the same disjoint row sets. The top block `n'_{K\{0}}` collects visible blocks of `Phi^+`,
which is harmless because top blocks are invisible. The residues are a subset. `[]`

**Consequences.** An 8-lift of a 6-configuration `Phi_6` at scale `w` forces `Phi_6` to be balanced at
scale `4w` and both 7-sub-configurations through `Phi_6` to be 7-lifts of `Phi_6` at scale `2w`, with
residues bounded by those of the 8-lift. So *non-existence of balanced bounded-residue 7-lifts implies
non-existence of 8-lifts*, and the integer test of Section 5 is a test of 7-lifts.

## 2. The k = 6 families E_t and their data

Notation of `L8.md` §3.1 (family of `cm_disjoint_mixed_pencil_search.md`, eqs. (1)-(6)). For distinct
nonzero rationals `t = (x_1, x_2, x_3, x_4)` with elementary symmetric functions `e_1..e_4`,

```text
P_lam(z) = (1 + lam e_2) z^2 - lam e_3 z + lam e_4,   Q_lam(z) = 1 + lam e_1 z - lam z^2,   f_lam = P_lam/Q_lam,
P_lam - z^2 Q_lam = lam prod_j (z - x_j)            (so f_lam(x_j) = y_j := x_j^2 for every lam),
f_a = f_0 = z^2,  f_b = f_{lam_b},  f_c = f_{lam_c},
E_t = normalisation of {(a,b,c) : f_a(a) = f_b(b) = f_c(c)},   Psi_t(a,b,c) = (x_1, x_2, x_3, a, b, c) ∈ Mbar_{0,6}.
```

Labels: anchors `A1, A2, A3` (positions `x_1, x_2, x_3`), movers `a, b, c`; the common value is `T = a^2`.

**Space-quartic model (used for all computations).** Put `Delta_b(T) = disc_z(P_b - T Q_b) = kb T (T - s)`
and `Delta_c(T) = kc (T - t')` (`L8.md` B.1(3)). Then

```text
E_t  ≅  { u^2 = kb (a^2 - s h^2),  w^2 = kc (a^2 - t' h^2) }  ⊂  P^3[a : u : w : h],
b = (-A1b(T) + a u)/(2 A2b(T)),   c = (-A1c(T) + w)/(2 A2c(T)),   T = a^2,
```

with `A_k(T)` the coefficients of `P - T Q`. This is a smooth complete intersection of two quadrics (genus
one, degree 4) for `t ∈ U` (`s, t'` distinct and nonzero). `etdata.py` checks that the map is well defined
at all points used.

**The 32 special points.** Over `T = y_j` the 8 points of `E_t` are `(±x_j, ±U_j, ±W_j)`, all rational. For
`j <= 3`, seven of them are the contacts `B_A`, `A = {Aj} ∪ S` (`S ⊆ {a,b,c}` the movers equal to `x_j`,
`S ≠ ∅`), and the eighth (`S = ∅`) is a non-contact point `aux_j`. For `j = 4`, the four with `|S| >= 2`
are the contacts `B_A`, `A = S ⊆ {a,b,c}`. The **cluster value** `v_A` (the common position of the cluster
`A` on the side with at most one anchor, in the frame where the anchors sit at `x_1, x_2, x_3`) is
`v_A = x_j`: for `j <= 3` this is the anchor's position, and for `j = 4` it is the common argument `x_4`,
which is not a marked point. So the 25 incidence data are

```text
(B_A, v_A) ∈ E_t x P^1:   7 points with v = x_1,  7 with v = x_2,  7 with v = x_3,  4 with v = x_4.     (2.1)
```

At `t_0 = (2, 3, -5, 1)`: `kb = -43680/14641`, `s = 369721/10920`, `kc = 660`, `t' = -364/165`,
`(U_1,W_1) = (1142/121, 64)`, `(U_2,W_2) = (1042/121, 86)`, `(U_3,W_3) = (622/121, 134)`,
`(U_4,W_4) = (1198/121, 46)` (`etdata.py`).

**The configurations.** For `t ∈ U(Q) \ Z_tor` and a non-torsion `P ∈ E_t(Q)`, `Phi_N := Psi_t(N P)`
anchored by `L8.md` Lemma B.6 has `Y_j = ±1`, bounded `Sigma`-supported residues, pairwise coprime
`Sigma`-free blocks `d_A(N)` with `log d_A(N) = alpha N^2 + O(N)`. These are the `(L_6)`-refuting
configurations. The union of the curves `E_t` is Zariski dense in `Mbar_{0,6}` (`L8.md` B.4), but for each
fixed `t` the configurations lie on one curve, and for bounded residues `t` must have bounded height.

## 3. Dimension / rank count for lifting to k = 7, 8

### 3.1 Curve level over a fixed E_t

Let `C ⊂ Mbar_{0,k'}` (`k' = 7, 8`) be an irreducible curve meeting `M_{0,k'}`, of class `c beta_{k'}`,
whose image under forgetting the `k'-6` new labels is `E_t`. For a boundary divisor `D'_A` of
`Mbar_{0,6}`, `pi^* D'_A = sum_{B ⊆ new} D_{A ∪ B}` has `2^{k'-6}` terms, each meeting `C` with
multiplicity `c`; since `E_t . D'_A = 1`, the degree of `C -> E_t` is

```text
d = c * 2^{k'-6}     (d = 2c for k' = 7,  d = 4c for k' = 8).
```

Each new label `n` has a coordinate `x_n` in the frame of `E_t` (anchors at `x_1, x_2, x_3`); `x_n` is the
pull-back of a cross-ratio, whose zero divisor on `Mbar_{0,k'}` is a sum of `2^{k'-4}` boundary divisors, so
`deg(x_n | C) = delta = c 2^{k'-4} = 4d`. Of the `d` points of `C` over `B_A` (counted with multiplicity),
exactly `d/2` lie on divisors containing `n`, and there `x_n = v_A` (Section 4.1 proves this incidence in
general).

**Count G1 (C of genus one).** By Riemann-Hurwitz `C -> E_t` is unramified, so over `Qbar` `C` is one of
the `sigma_1(d)` connected unramified covers (3 for `d = 2`, 7 for `d = 4`). Write `x_n = N/D` with
`N, D ∈ H^0(C, M)`, `deg M = delta`: `2 delta - 1` projective parameters plus 1 for `M ∈ Pic^delta(C)`.
Conditions: `d/2` linear conditions `N - v_A D = 0` at chosen preimages, for each of the 25 contacts (a
finite choice of which preimages). For `k' = 8` add, for each of the six old labels `i`, the `c` conditions
that `c` of the `2c` zeros of `x_6 - x_i` not lying over contacts are also zeros of `x_7 - x_i` (the
`D_{i67}` contacts); the `D_{67}` contacts are then automatic by degree.

| lift | class | `d` | parameters | conditions | expected dimension per `E_t` | over the 3-dim family |
|---|---|---|---|---|---|---|
| one point | `beta_7` | 2 | 16 | 25 | **-9** | -6 |
| one point | `2 beta_7` | 4 | 32 | 50 | **-18** | -15 |
| two points | `beta_8` | 4 | 64 | 100 + 6 | **-42** | -39 |
| two points | `c beta_8` | `4c` | `64c` | `100c + 6c` | `-42c` | `3 - 42c` |

(General `c`, one point: `8d - 25 d/2 = -4.5 d`.) This refines the count "overdetermined by 10" of
`L8.md` §5 (which omitted the `Pic` parameter) and agrees with `L7.md` §5 Lemma 5 (`-8` for genus 0 with
fixed branch points, `-9` here for genus one).

**Count L (any genus, class `beta_7`, the linear relaxation).** Section 4.1 shows that any such `C` gives a
curve `Gamma ∈ |L ⊠ O(2)|` on `E_t x P^1` through the 25 points (2.1), `L ∈ Pic^8(E_t)`. Parameters
`3 * 8 - 1 + 1 = 24`, conditions 25: **expected `-1`**. Theorem 4.2 shows the actual answer is `+1`
(a pencil), the excess being forced by the pencil identity `P_lam - z^2 Q_lam = lam prod(z - x_j)`.

**Generic comparison (not over `E_t`).** Genus-one curves of class `beta_k` have expected dimension
`-K.beta_k = 8, 7, -10` for `k = 6, 7, 8`; genus-zero ones `-K.beta + (k-6) = 8, 8, -8`
(`ptolemy.md` §4, `L7.md` §5). So at `k = 8` no balanced curve of genus `<= 1` is expected anywhere, which
is the function-field form of `(L_8)` (`L8.md` §7).

### 3.2 Integer level over a fixed configuration

Fix a balanced 6-configuration `Phi_6` with bounded residues at scale `2w` (for a 7-lift at scale `w`). In
the `E_t` frame, a new point is `x = u/v` with `gcd(u,v) = 1` and, by Lemma 6.2 of `ptolemy.md` applied to
the 8 divisors in the zero divisor of its cross-ratio, `log max(|u|,|v|) = 8w + O(eta w + H_Sigma + 1)`:
a box of `e^{16w}` pairs. At a prime `p ∤ Sigma` with `p | d_A`, the lift joins the cluster `A` to depth
`v_p(n_{A ∪ new})`; since the anchors are `p`-adically distinct this is the congruence

```text
u β_A - α_A v ≡ 0  (mod e_A),     v_A = α_A/β_A,   e_A := n_{A ∪ new} | d_A,   e_A ≈ d_A^{1/2}.      (3.1)
```

The new pair contacts (`x` with one old label) impose no congruence (their blocks are free cofactors). So
the admissible `(u,v)` form, for each choice of `(e_A)`, a rank-2 lattice `Lambda(e)` of index
`prod e_A ≈ e^{25 w}`, and a lift needs a primitive vector of sup-norm `<= e^{8w(1+o(1))}`:

```text
expected number of balanced 7-lifts of one Phi_6  =  e^{(16 - 25) w + o(w)}  =  e^{-9w + o(w)}      (L7.md §8),
expected number of balanced 8-lifts of one Phi_6  =  e^{-42 w + o(w)}       (two such points at scale 2w, -18w each, and 6 joint congruences).
```

(The `e^{o(w)}` is the number of choices of the `e_A`, by the divisor bound.) Bookkeeping check: the
balanced 6-configurations whose blocks factor as required number `e^{32 w}` per shift class without top
block (`8 * 4w`), and `32 - 42 = -10`, which is the `k = 8` exponent `1 - K.beta_8 = -9` without the top
block. The one-point check is `8 * 2w - 9w = 7w`, the `k = 7` exponent without top block (`L7.md` §8).
So the curve counts of 3.1 and the integer counts agree exactly (`-9`, `-42`).

### 3.3 What the counts can and cannot say

* **Over- or under-determined?** Every lift of the `E_t` mechanism is overdetermined: by 9 (one point,
  genus one, per curve), by 18 (class `2 beta_7`), by 42 (two points), and by 1 even in the linear any-genus
  relaxation; the integer counts per configuration are the same numbers in units of `w`.
* **Counts are beaten by structure.** The same per-point count applied to `E_t` itself predicts nothing:
  a mover `a = u_a/v_a` of `Phi_N` has `h(a) = 4W` (box `e^{8W}`) and satisfies congruences modulo the 15
  contacts `d_A` with `a ∈ A` (index `e^{15W}`), expected number `e^{-7W}`, yet `Phi_N` exists for every
  `N`. A curve of genus `<= 1` with positive rank defeats any per-point count. Likewise Count L predicts
  `-1` and the truth is `+1`. **So the counts do not decide the lifting question; Section 4 decides it at
  the curve level by exact computation, and Section 5 tests the integer level on explicit members.**
* **Other bounded-residue families.** The real genus-zero balanced 6-curves (`ptolemy.md` §5.4, `L7.md`
  §7; NUMERICAL, none known over `Q`) lift to `beta_7` with expected dimension `8 + 2 - 8 = 2` in total
  (`L7.md` §5, one involution), and from such a 7-curve to `beta_8` with `33 + 2 - 56 = -21` per curve,
  `-19` in total. No bounded-residue family at `k = 7` is known (`(L_7)` open). The dense integer families
  at `k = 5` (`L8.md` Prop. C) are 5-sub-configurations of everything, so lifting them is `(L_8)` itself
  (count `1 - K.beta_8 = -9`).

## 4. Curve level: the exact classification of class-beta_7 lifts of E_t

### 4.1 The incidence lemma

**Lemma 4.1.** Let `t ∈ U` and let `C ⊂ Mbar_{0,7}` be an irreducible curve meeting `M_{0,7}`, with
`C . D_T = 1` for all 56 boundary divisors, whose image under the map `pi` forgetting the new label `6` is
`E_t` (with its labels `A1, A2, A3, a, b, c`). Let `x_6` be the position of label 6 in the frame where the
anchors sit at `x_1, x_2, x_3` (the forgetful morphism to `Mbar_{0,{6,A1,A2,A3}} ≅ P^1`, normalised so
that the boundary point `{6, Aj} | rest` is `x_j`). Then:

1. `C -> E_t` has degree 2, and `x_6|_C` has degree 8;
2. the image `Gamma` of `C` under `(pi, x_6)` in `E_t x P^1` is irreducible, `C -> Gamma` is birational,
   and `Gamma ∈ |L ⊠ O_{P^1}(2)|` for some `L ∈ Pic^8(E_t)`;
3. `Gamma` contains the 25 points `(B_A, v_A)` of (2.1).

*Proof.* (1) `pi^* D'_A = D_A + D_{A ∪ 6}`, so `pi_*[C] . D'_A = 2`, while `E_t . D'_A = 1`. The zero
divisor of `x_6 - x_j` is `sum D_T` over the `2^3 = 8` sets `T ∋ 6, Aj` missing the other two anchors, so
`deg x_6|_C = 8` (`C` is not contracted, since the sum meets it).
(2) On `M_{0,7}` the map `(pi, x_6)` is injective (a configuration is determined by its first six points and
the position of the seventh relative to three of them), and `C` meets `M_{0,7}`. So `C -> Gamma` is
birational, `Gamma . (E_t x {v}) = 8` and `Gamma . ({e} x P^1) = 2`. Since
`Pic(E_t x P^1) = pr_1^* Pic(E_t) ⊕ pr_2^* Pic(P^1)`, `Gamma ∈ |L ⊠ O(2)|` with `deg L = 8`.
(3) `C . D_{A ∪ 6} = 1` and `C ⊄ D_{A ∪ 6}`, so there is `c ∈ C ∩ D_{A ∪ 6}`. Forgetting 6 maps
`D_{A ∪ 6}` into `D'_A`, and `E_t ∩ D'_A = {B_A}` with `B_A` in the open stratum of `D'_A` (`L8.md` Lemma
B.2(2)), so `pi(c) = B_A`. In the stable curve `c`, label 6 lies in the subtree containing the cluster `A`,
separated from `A^c`, which contains at least two anchors. If `A ∋ Aj`, the 4-point stabilisation of
`(6, A1, A2, A3)` is the boundary point `{6, Aj} | rest`, so `x_6(c) = x_j = v_A`. If `A` contains no anchor,
the component carrying the three anchors survives stabilisation and the subtree containing 6 contracts to
the node, whose position is the cluster value `x_4 = v_A`. `[]`

So the existence of such a `C` forces a nonzero section `sigma ∈ H^0(E_t x P^1, L ⊠ O(2))` vanishing at the
25 points (2.1). Writing `sigma(., v) = sum_{j<=3} l_j(v) sigma_j` (Lagrange basis at `x_1, x_2, x_3`), the
conditions are `sigma_j ∈ H^0(L(-sum_{A: v_A = x_j} B_A))` (a line, since `deg = 1`), spanned by `tau_j`,
and `sum_j l_j(x_4) sigma_j` vanishing at the four `B_A` with `v_A = x_4`. Hence

```text
a nonzero sigma exists   <=>   rank [ tau_j(B_{4,k}) ]_{k = 1..4, j = 1..3}  <=  2.            (4.1)
```

**The pencil solutions.** For every `t` and `lam`, the section
`sigma_lam := h^2 P_lam(v) - a^2 Q_lam(v)` of `O(2H) ⊠ O(2) = T^*O(1) ⊠ O(2)` vanishes at every point
`(P, x_j)` with `T(P) = y_j`, because `f_lam(x_j) = y_j`. So for `L = 2H` the solution space contains the
2-dimensional pencil `Pi_t = <sigma_0, sigma_inf>`.

### 4.2 Theorem 4.2 (classification)

**Theorem 4.2.** There is a nonempty Zariski-open `U_0 ⊆ U` such that for `t ∈ U_0` the following hold.
`U_0` contains `t_0 = (2,3,-5,1)`, item 109's `t = (1/3,3,2,5/3)`, and the seven further parameter points
listed in Section 7.

(a) For every `L ∈ Pic^8(E_t ⊗ Qbar)`, the only sections of `L ⊠ O(2)` vanishing at the 25 points (2.1)
are `0`, unless `L ≅ O(2H) = T^*O(1)`, in which case they form exactly the pencil `Pi_t`.

(b) Every irreducible curve `C ⊂ Mbar_{0,7}` of class `beta_7` meeting `M_{0,7}` whose image after
forgetting the new label is `E_t` is a **pencil lift**

```text
C'_lam = closure of { (x_1, x_2, x_3, a, b, c, d) : f_a(a) = f_b(b) = f_c(c) = f_lam(d) },
```

with `lam ∉ {0, lam_b, lam_c}` and `f_lam` non-degenerate. Its geometric genus is 9, or 5 for at most two
values of `lam` (those for which `f_lam` has `s` or `t'` as a branch value). **In particular no class-`beta_7`
curve of genus `<= 4` lies over `E_t`.**

(c) Conversely, for `lam` outside a finite set, `C'_lam` is a smooth curve of genus 9 and class `beta_7`.

*Proof of (a) at an explicit `t` (exact computation, `beta7lift2.py`).*

* *Parametrisation of `Pic^8`.* Fix the three non-contact points `aux_1, aux_2, aux_3` (Section 2). The map
  `q -> L_q := O(3H - aux_1 - aux_2 - aux_3 - q)` is a bijection `E_t -> Pic^8(E_t)`. Since the space
  quartic is projectively normal, `H^0(O_E(3H))` is the 12-dimensional space of cubic forms in the normal
  monomials `a^i h^k u^e w^d` (`e, d <= 1`), and `H^0(L_q)` is the subspace vanishing at the four base
  points (with multiplicity when `q` is one of them).
* *Generic `q`.* For `q` outside the 28 special points (`aux_j` and the 25 `B_A`), eleven distinct points
  impose independent conditions on cubics, so `tau_j` is the cofactor vector of an `11 x 12` matrix whose
  last row is the monomial vector of `q`. Thus `G(q) := [tau_j(B_{4,k})]` has entries that are cubic forms
  in `q`, and its four `3 x 3` minors `m_i` are degree-9 forms on `E_t` (36 zeros each).
* *Forced zeros, valid for every `t`.* `m_i` vanishes to order `>= 3` at each `aux_j` (all columns vanish:
  the `q`-row equals a base row); to order `>= 1` at each `B_A` with `v_A = x_j`, `j <= 3` (column `j`
  vanishes); to order `>= 1` at `B_{4,k}` for the three minors containing row `k` (`F_j(q) = 0` for all `F`
  vanishing at `q`); and to order `>= 2` at `q* := B_{bc}` (the residual point of the plane through
  `aux_1, aux_2, aux_3`, so `L_{q*} = O(2H)`; there the pencil gives `rank <= 1`). The coplanarity
  `det(aux_1, aux_2, aux_3, B_bc) = 0` is a symbolic identity in `x_1..x_4` (`coplanar.py`).
* *Residual zeros.* With `N(m_i) ∈ Q[a]` the norm from `E_t` to the `a`-line (product over `u -> ±u`,
  `w -> ±w`), the forced zeros account for 34, 34, 34 and 35 of the 36 zeros, and the residual polynomials
  `rho_i = N(m_i)/prod (a - a(P))^{ord}` (exact division) have degrees `2, 2, 2, 1` and
  **`gcd(rho_1, ..., rho_4) = 1` over `Q`**, with no common zero at `h = 0`. So no non-special `q ∈ E_t(Qbar)`
  has `rank G(q) <= 2`.
* *Special `q`.* For each of the 28 special points the dimension of the solution space was computed
  directly (36 unknowns, with order-2 vanishing along the curve where a base point coincides with a
  `B_A`): it is `0` for 27 of them and `2` for `q* = B_bc`, where it therefore equals `Pi_t`.

This proves (a) at the nine parameter points (all PASS, `batch_beta7.json`).

*Openness.* The forced orders, the identity `q* = B_bc`, and the lower bound `dim >= 2` at `q*` hold for all
`t ∈ U` (structural). The remaining conditions are open in `t`: that the binary forms `rho_i(t)` have no
common zero on `P^1` (the bad set is the image of a closed subset of `U x P^1` under a proper map), and that
the 28 direct systems have maximal rank (respectively rank `36 - 2` at `q*`). Each holds at `t_0`, so all
hold on a nonempty open `U_0 ∋ t_0`.

*Proof of (b).* By Lemma 4.1 and (a), `Gamma = {sigma_lam = 0}` for some `lam ∈ P^1`, i.e. `Gamma` is the
fibre product `E_t x_{P^1_T} P^1_v` of `T` and `f_lam`. It must be irreducible of degree 2 over `E_t`.
For `lam ∈ {0, lam_b, lam_c}` it splits into the graphs of `m` and `m-bar` (`m = a, b, c`), and for the
finitely many degenerate `lam` (`P_lam`, `Q_lam` with a common root `z_0`) it contains `E_t x {z_0}`. The
function field of `E_t` over `Q(T)` is cut out by the quadratic characters with branch pairs
`{0,inf}, {0,s}, {inf,t'}`; they generate the characters of the six pairs in `{0, inf, s, t'}`. So
`Gamma_lam` is irreducible iff the branch pair `{w_1, w_2}` of `f_lam` is not contained in
`{0, inf, s, t'}`. Branch value `0` occurs only for `lam ∈ {0, lam_b}` and `inf` only for `lam ∈ {0, lam_c}`
(`L8.md` B.1(3)); the pair `{s, t'}` occurs iff `Delta_lam(s)` and `Delta_lam(t')` have a common root in
`lam`, which is a closed condition on `t` not satisfied at `t_0` (shrink `U_0`). For the other `lam`, the
normalisation of `Gamma_lam` is a double cover of `E_t` branched exactly over the 8 preimages of each new
branch value, so `2g - 2 = 8 * #{new branch values}` and `g ∈ {5, 9}`; `g = 5` iff one branch value is
`s` or `t'` (at most one further root of `Delta_lam(s) = 0` besides `lam_b`, and one of `Delta_lam(t') = 0`
besides `lam_c`). Since `C -> Gamma` is birational, `g(C) ∈ {5, 9}`.

*Proof of (c).* The argument of `L8.md` Lemma B.2 applies with four movers: over `y_j` the 16 points
`(±x_j, {x_j, b'}, {x_j, c'}, {x_j, d'})` give the 15 contacts `{Aj} ∪ S` (`j <= 3`) resp. the 11 contacts
`S`, `|S| >= 2` (`j = 4`), i.e. all `3 * 15 + 11 = 56` boundary divisors once each; two pencil members
meet only at the `x_j` (B.1(2)), and transversality needs the four slopes `f_m'(x_j)` to be nonzero and
distinct, which fails for finitely many `lam` (checked to hold at `lam ∈ {1, 2, -1, 1/7}` for `t_0` and
item 109; `pencil_lift.py`, which also confirms genus 9 there). `[]`

**Remarks.**
1. The pencil lifts are item 112's "four maps of a pencil" (`seven_row_quadratic_quotient_obstruction.md`;
   `L8.md` §5). What is new is that they are the **only** class-`beta_7` lifts of `E_t`, of any genus: the
   linear count `-1` of Section 3.1 is beaten exactly by the pencil identity and by nothing else.
2. Every rational seventh point `d` over a rational point of `E_t` lies on the pencil lift with
   `lam(d) = -(d^2 - T)/(P_1(d) - T Q_1(d))`. The pencil lifts therefore contain all 7-lifts as *points*;
   the theorem says that no curve of genus `<= 1` (indeed `<= 4`) organises them.

### 4.3 Corollary 4.3 (8-point lifts)

**Corollary 4.3.** Let `t ∈ U_0` and let `C ⊂ Mbar_{0,8}` be an irreducible curve of class `beta_8` meeting
`M_{0,8}` whose image after forgetting the new labels `6, 7` is `E_t`. For `n ∈ {6, 7}` let `C_n` be the
image after forgetting the other new label. Then either `C -> C_n` has degree 2 and `C_n` is a pencil lift
`C'_lam` (genus 5 or 9), or `C -> C_n` is birational and `C_n` has class `2 beta_7` and degree 4 over `E_t`.
If `g(C) <= 4`, the second alternative holds for both `n`.

*Proof.* Forgetting one label pulls each 7-point boundary divisor back to the sum of two 8-point ones, so
`pi_*[C] = 2 beta_7` and `deg(C -> C_n) [C_n] = 2 beta_7`. Degree 2 gives a class-`beta_7` curve over
`E_t`, a pencil lift by Theorem 4.2(b). A nonconstant map `C -> C_n` forces `g(C) >= g(C_n) >= 5`. `[]`

The pencil mechanism does lift to `beta_8`: `{f_a(a) = f_b(b) = f_c(c) = f_lam(d) = f_mu(e)}` is a
`(Z/2)^5`-cover of `P^1_T` branched over 8 points, of genus `1 + 2^{5-2}(8-4) = 33` for generic
`(lam, mu)`, and of class `beta_8` by the same contact argument (`3 * 31 + 26 = 119` contacts). Genus-one
8-lifts of `E_t` would have to be of the second kind in Corollary 4.3, which Count G1 puts at expected
dimension `-18` for each 7-projection and `-42` for the curve. That case is **not** decided here.

## 5. Integer level: exact tests on explicit members (Data)

Members: `t = t_0`, `Phi_N = Psi_{t_0}(Q_1 + N (Q_2 - Q_1))` (the orbit of `L8.md` §3.3), `N = 2..8`.
For each `N` (`intlift.py`):

* the 25 contact integers `d_A(N)` are computed by gcd/Möbius inversion (`profile.py` of `L8.md`), their
  `Sigma`-free parts are factored completely (`Sigma` = the structural bad set of `L8.md` §3.3), and
  `w = w_6/2`, `w_6` = mean log of the `Sigma`-free parts (scale of a 7-lift);
* **sanity check of the dictionary (3.1):** each mover `m = u_m/v_m` of `Phi_N` satisfies
  `u_m beta_A - alpha_A v_m ≡ 0 (mod d_A)` for every contact `A ∋ m` (full `d_A`): passed for all `N`;
* **structure beats count, measured:** for the mover `a`, the lattice of (3.1) with the full `d_A`
  (`A ∋ a`, 15 contacts) has index `e^{(≈30) w}`, yet its shortest vector is `(u_a, v_a)` itself, of size
  `e^{(≈8) w}` (column "mover a");
* **the 7-lift test:** for every contact `A`, the divisors `e_A | d_A` with
  `|log e_A/log d_A - 1/2| <= 1/4` are listed. This window contains `e_A = n_{A ∪ 6}` for every 7-lift with
  balance `eta <= 1/5` whose residues are supported on `Sigma`, since then
  `log e_A/log d_A ∈ [(1-eta)/(2(1+eta)), (1+eta)/(2(1-eta))] ⊆ [1/4, 3/4]`. Contacts with no such divisor
  are counted ("unsplittable"; any balanced lift is then impossible for that `N`) and, to keep testing
  the lattice part, dropped (`e_A = 1`, which only weakens the necessary condition). Over all splits (or
  4000 random ones), the shortest vector of `Lambda(e)` is computed exactly by Gauss reduction. A
  balanced 7-lift with `H_Sigma = o(w)` needs `log lambda_1 <= (8 + o(1)) w` for some split.

| `N` | `w` | contacts with no divisor in window | contacts used | splits tested | max skew / `w` | mean skew / `w` | mover `a`: `h/w`, `log lambda_1/w`, index`/w`, skew`/w` |
|---|---|---|---|---|---|---|---|
| 3 | 7.4 | 16 | 9 | 4000 | 1.511 | 0.213 | 8.96, 8.96, 31.83, 6.95 |
| 4 | 15.1 | 11 | 14 | 4000 | 0.285 | 0.039 | 8.74, 8.74, 30.52, 6.53 |
| 5 | 26.4 | 8 | 17 | 4000 | 0.165 | 0.023 | 8.30, 8.30, 30.12, 6.76 |
| 6 | 40.3 | 6 | 19 | 4000 | 0.124 | 0.014 | 8.12, 8.12, 30.09, 6.93 |
| 7 | 56.5 | 8 | 17 | 4000 | 0.092 | 0.010 | 8.09, 8.09, 30.00, 6.91 |
| 8 | 75.3 | 9 | 16 | 4000 | 0.052 | 0.008 | 8.09, 8.09, 29.95, 6.89 |

Skew of a lattice of index `I` with shortest vector `lambda_1`: `s = (log I / 2 - log lambda_1)/w`. A random
rank-2 lattice has `s = O(1/w)` (probability of `s >= sigma` is about `e^{-2 sigma w}`). A balanced 7-lift
with all 25 contacts needs `s >= 25/2 - 8 = 4.5` (up to `O(eta) + o(1)`); with only `k` contacts kept, it needs
`s >= k/2 - 8 + o(1)`.

Reading of the data (evidence only; `N = 2` has all 25 contacts prime or 1 and is omitted):

* **No balanced lift exists for any `N <= 8`, already for factorisation reasons.** In every case between 6
  and 16 of the 25 contact integers have no divisor in `[d^{1/4}, d^{3/4}]`, so no 7-lift with balance
  `eta <= 1/5` and `Sigma`-supported residues exists, and by Lemma 1.1 no such 8-lift either (with
  general bounded residues the statement holds up to the `O(H_Sigma)` mass that residual primes can move). Asymptotically this obstruction only costs a
  constant factor `c(eta)^{25}` per configuration (a positive proportion of integers has a divisor near
  their square root), so it is not the real reason. It is recorded because it is exact.
* **On the splittable contacts the lift lattices are random.** Over 4000 splits per `N`, the maximal skew
  falls from `1.5 w` (`N = 3`, `w = 7.4`) to `0.05 w` (`N = 8`, `w = 75.3`), and the mean is `<= 0.01 w` from `N = 6` on. In
  absolute terms the skew stays `O(log 4000)`, as it does for random lattices. A lift would need about
  `0.5 w` to `1.5 w` even with the reduced contact sets (`k = 17..19`; at `k = 16`, `N = 8`, the reduced test is void and only the factorisation obstruction applies).
* **Contrast: the movers of `Phi_N` themselves.** The lattice of the mover `a` (full `d_A`, the 15 contacts
  `A ∋ a`) has index `e^{30 w}`, but its shortest vector is `(u_a, v_a)`, of size `e^{8 w}`: skew `≈ 7 w`
  for every `N`. That is `30/2 - 8 = 7` exactly, the signature of the curve `E_t`. A curve-level lift
  would produce the same kind of systematic skew in the lift lattices. None is seen, which is consistent
  with Theorem 4.2 and with the counts.

## 6. Interpretation

### 6.1 The answer to the lifting test

| level | lift | result | status |
|---|---|---|---|
| curve, one point | class `beta_7` over `E_t`, any genus | **only the pencil lifts `C'_lam`, genus 5 or 9** | Theorem 4.2 (`t ∈ U_0`, incl. 9 explicit `t`) |
| curve, one point | class `beta_7`, genus `<= 4` (incl. the genus-one positive-rank mechanism of Thm B) | **none** | Theorem 4.2 |
| curve, one point | class `c beta_7`, `c >= 2`, genus one | expected dimension `-4.5 d = -9c` per `E_t` | Count G1; open |
| curve, two points | class `beta_8`, genus `<= 4` | must have both 7-projections birational of class `2 beta_7`; expected dim `-42` | Cor. 4.3 + Count G1; open |
| curve, two points | class `beta_8`, pencil | exists, genus 33 | Section 4.3 |
| integer, one point | balanced 7-lift of a given `Phi_6` | expected number `e^{-9w}` | Count 3.2 |
| integer, two points | balanced 8-lift of a given `Phi_6` | expected number `e^{-42w}` | Count 3.2 |
| integer, explicit members | `Phi_N(t_0)`, `N <= 8` | no balanced lift (see Section 5) | Data |

So: **the `(L_6)` mechanism does not lift at the first step.** A balanced curve of class `beta_7` over
`E_t` exists, but only through the pencil identity, and it has genus 5 or 9; it carries finitely many
rational points (Faltings), and the bounded-residue configurations of Theorem B do not propagate along it.
Every count, at the curve level and at the integer level, is overdetermined (by 9 per added point at
`k = 7` and by 42 for two points at `k = 8`), and the explicit members behave like random lattices,
without the structural skew that the movers of `Phi_N` themselves exhibit.

### 6.2 What a positive lift would have meant

A positive lift means a family of balanced 8-configurations with `H_Sigma = O(1)` and `w -> infinity`
over the `E_t`, e.g. positive-rank genus-one curves of class `c beta_8` over each `E_t` in a 3-parameter
family. Their union would lie in a closed `W ⊂ Mbar_{0,8}` with `dim W <= 3 + 1 + (dimension of the family
of lifts of one E_t)`. If `W ≠ Mbar_{0,8}` (in particular whenever each `E_t` has only finitely many such
lifts, `dim W <= 4`), `(L_8)` would survive with `Z ⊇ W`, and `W` would be the first explicit component of
the exceptional set. Only a dense `W` would make `(L_8)` false, and then `V_8` would be false too
(`ptolemy.md` Thm D shows `V_8 => (L_8)`). Theorem 4.2 rules out the cheapest version of this (class
`beta_7`, genus `<= 4`), and the counts make the remaining versions expected-empty.

### 6.3 What the obstruction means, and what it does not

* **It is unconditional, but it is a curve (function-field) statement.** Theorem 4.2 says: over `E_t`
  (`t ∈ U_0`), every balanced class-`beta_7` curve has genus `>= 5`. That is the function-field shadow of
  "the `(L_6)` counterexamples do not lift": it excludes the only mechanism known to produce infinitely
  many bounded-residue configurations (rational points on positive-rank curves of genus `<= 1`).
* **It is not a piece of `(L_8)` (nor of `(L_7)`) for points.** `(L_8)` concerns rational points, which
  need not lie on curves of small genus. The point-level analogue would be: every balanced 7- or 8-lift of
  `Phi_N(t)` has `H_Sigma >= c w`. By Section 3.2 this is a statement that the lattices `Lambda(e)`, built
  from the contact integers `d_A(N)` (values of elliptic divisibility sequences), are never skewed by
  `e^{4.5 w}`. That is a Diophantine statement of the same kind as `(L_8)` itself. The constraints of the
  task exclude the obvious tools: the primes of `d_A(N)` are unbounded, so a fixed-`S` Subspace theorem does
  not apply with a uniform exceptional set, and `t` (hence `Sigma`) must stay fixed for bounded residues, so
  any threshold depending on `t` is template-dependent, not uniform. No unconditional piece of `(L_8)`
  results, and none of the growth statements of the record changes.
* **Why the obstruction is not vacuous.** The linear count for class `beta_7` is `-1`, and structure has
  already beaten counts twice in this family (Thm B itself, and the pencil, `+1` against `-1`). Theorem 4.2
  is the exact computation showing that the pencil is the only such structure at `k = 7`.

### 6.4 Consequences for the next steps

1. The census of small-height configurations at `k = 8` (round-17 task 1) should treat the pencil lifts of
   genus 9 and 33 as known sources of balanced curves over the `E_t`. They contain only finitely many
   rational points each, so they are not candidates for an infinite bounded-`H` family on one curve.
2. The only curve-level loophole left over `E_t` is a genus-one curve of class `2 beta_7` (degree 4 over
   `E_t`, expected dimension `-18`), and then its 8-point combination (`-42`). Via the quotient `C/H`
   (`H ⊂ Gal(C/E_t)` of order 2) this becomes a linear system of 50 conditions on 48 parameters for each of
   `3^{25}` incidence patterns, before the genus condition. That is not feasible exactly, and not
   recommended, given the count.
3. Nothing here moves the growth rate or the uniform bound. The judge's assessment of `(L_8)`
   (`round17/README.md` §7.3) is unchanged. One concrete mechanism for a counterexample has been removed.

## 7. Scripts and reproducibility

All under `round18/lift_checks/`. Run with `python3` from that directory. Arithmetic is exact (`Fraction`,
`flint.fmpq_poly`, `fmpq_mat`, `nmod_mat`, `fq_default`) except the logarithms used to report sizes.

| script | content | type |
|---|---|---|
| `fam.py`, `ell.py`, `lift.py`, `profile.py` | copied from `threshold/L8/k6fam/` (family, quotient curve, orbit lifts, contact profile) | exact |
| `etdata.py` | space-quartic model of `E_t`; the 25 contacts `(B_A, v_A)`; `aux_j` | exact |
| `beta7lift.py` | first version of the (4.1) test over `F_p` (base points `aux_j`); found the pencil at `q = B_bc` | exact (`F_p`) |
| `probe1.py`, `probe2.py` | identification of the nullity-2 point | exact (`F_p`) |
| `beta7lift2.py` | Theorem 4.2(a) over `Q`: minors, norms, forced zeros, `gcd(rho_i) = 1`, infinity, 28 direct systems | exact (`Q`) |
| `batch_beta7.py` → `batch_beta7.json` | Theorem 4.2(a) at 9 parameter points (all PASS; 3 random draws outside `U` skipped) | exact (`Q`) |
| `coplanar.py` | symbolic identity `det(aux_1, aux_2, aux_3, B_bc) = 0` and the pencil identity | exact (symbolic) |
| `pencil_lift.py` | pencil lifts: transversality slopes, branch values, genus 9 | exact |
| `branch_st.py` | at all 9 parameter points: no pencil member has branch pair `{s, t'}`; the genus-5 members | exact |
| `mktable.py`, `show_int.py` | tables of Section 5 from `intlift.json` | reporting |
| `intlift.py` → `intlift.json` | Section 5 data | exact integers; logs for reporting |

Parameter points of Theorem 4.2 (all PASS): `(2,3,-5,1)`, `(1/3,3,2,5/3)`, `(-2,7/5,-3,1)`,
`(2,9/5,6/5,1)`, `(5/2,-9/5,-7,1)`, `(-3/4,-1/3,1/2,1)`, `(-3,1/2,4,1)`, `(3/5,2,-3/2,1)`,
`(7/5,3/2,7/3,1)`.
