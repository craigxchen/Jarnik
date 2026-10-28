# Round 18 theory: the arithmetic Ptolemy lemma (L_8)

Task (round-17 judge §7.3): theory on `(L_8)`; proof lane or refutation lane. Scripts in
`theory_checks/`. Companion round-18 reports: `census.md` (small-height census) and `lift.md`
(lifting test), cited by section.

## 0. Summary and verdict

**Neither lane decides `(L_8)`. No growth improvement and no uniform bound are claimed.** The
refutation lane got further in the sense that it ends at a sharp statement: refuting `(L_8)` is
the same as refuting Vojta's `D = 0` inequality at balanced points of `Mbar_{0,8}`, and every
family-type or density-type mechanism is excluded by proof. The proof lane ends at a concrete
lattice statement with no structure left beyond the balanced class. One gap in the record's
barrier map (abc on auxiliary pencils) is identified and partly closed.

1. **Formulation (§1).** `(L_8)` is used verbatim from `ptolemy.md` §7 and checked against
   Steps 0-5. *Lemma 1.1 (new, exact check on all 39208 trees with 8 leaves):* the item-532
   "two nested chains on disjoint row sets" condition is anchor-independent and equals
   "the `p`-adic tree is a caterpillar". With `census.md` Prop 1.3, any version with `c w`
   replaced by `f(w) -> infinity` is no easier than `(L_8)` (rate dichotomy).
2. **What the joint system contains (Prop 2.1).** The 70 Ptolemy relations jointly are the
   point `Phi`; anchoring is free; on pairwise coprime profiles (which every proof must cover)
   the two nested chains, the per-prime gcd/valuation bookkeeping and the residues are vacuous.
   **No structure of the joint system other than `K.beta = 10` survives in the worst case, and
   none yields a positive multiple of `w`.**
3. **Lattice form (Lemma 2.2, Prop 2.3; new).** In Kapranov's chart,
   `h_K = sum_{|T|>=3} (|T|-2) ell_T - 6 log|x|` and the congruence lattice of `x` has index
   `prod d_T^{|T|-1}`. A balanced point is a vector of its own congruence lattice `e^{18w}` below
   the Minkowski scale (`e^{31w}` against `e^{49w}`; trivial lower bound `e^{30w}`). Hence every
   density mechanism (boxes, sieves, Ekedahl–Poonen, geometry of numbers) yields only
   `h_K <= -sum ell_T < 0`, and the heuristic exponent `-K.beta = 98 - 108 = -10` is
   (profiles) minus (shortness).
4. **Restricted theorems (§2.3).** Lemma 2.4: the matching coefficients of every quadruple are
   pairwise coprime for every `(L_8)`-datum. Thm 2.5: for a **fixed** prime set `S_0`
   (Ridout), `H_Sigma >= 8(1-eta)(1-eps) w - C(S_0)` with no exceptional set; useless for the
   circle problem because `C` depends on `S_0` (lit18 Props A/B), and stated only as such.
5. **abc lane (§2.4; new gap, partly closed).** The refereed abc barrier covers only relations
   among products of differences; "abc + an auxiliary identity" is not covered by it or by the
   Liouville barrier. Prop 2.6 (conditional on abc): a unit pencil with
   `gain = L.beta - F_red.beta > 0` gives `(L_8^pc)`; equivalently
   `gain = 10 - K^log_{X~/P^1}.beta~`. Excluded: all morphisms (gain `<= -32`, via Bruno–Mella,
   cited from the abstract) and all unit pencils of degree `<= 4` (data: maximum gain `-32` over
   277940 pencils, attained e.g. by `((x_1-x_2)/(x_1-x_3))^D`). Open in general.
6. **Refutation lane (§3).** Expected dimensions `-10d + 2` (genus 0) and `-10d` (genus 1),
   versus `7d + 1` at `k = 7`; in Kapranov's coordinates the count is `(186d + 2) - 196d`, i.e.
   literally `-K.(d beta)`. Covering families need `K.gamma <= 2g - 2`; density mechanisms give
   `h_K < 0`; lifts of the `(L_6)` families fail (`lift.md`). What is left is a Zariski-dense
   supply of rigid curves of genus `<= 1` with `K.gamma > 0` (the rigid curves with `K > 0` that I
   computed, the F-curves of types `(2,2,2,2)` and `(3,2,2,1)`, lie in the boundary) or of sporadic
   exceptionally short lattice vectors, either of which refutes `V_8`. **The formulation audit
   (§3.4) finds no loophole**: `O(1)`, `eta_0`, the exceptional set, residues, anchoring,
   invisible blocks and chains are all handled by the refereed `V_8 ⇒ (L_8)`.
7. **Precise residual (§5).** `(L_8^pc)` in lattice form; and the existence of a pencil of
   positive abc-gain (which would make the uniform theorem conditional on abc, after extending the
   gcd bookkeeping to nested chains).

## 1. The object: exact formulation of (L_8) used here

Conventions as in `L8.md`: *Theorem/Proposition/Lemma* = proved here with all quantifiers, or
cited with the place of proof; *Data* = exact finite computation; *Heuristic* = a count never
used as proof. `census.md` and `lift.md` (this round, concurrent) are cited by section.

### 1.1 Statement (verbatim from `ptolemy.md` §7, as in `census.md` §1.1)

> **(L_8).** There are `c > 0`, `eta_0 > 0` and a proper Zariski-closed `Z ⊂ Mbar_{0,8}` such
> that every datum `D = (P_0 = (1,0), P_1, ..., P_7; (n_T)_{∅≠T⊆[7]}; (t_ij))` of primitive
> integer vectors with `det(P_i,P_j) = t_ij prod_{T∋i,j} n_T` (`t_ij ∈ Z\{0}`),
> `Y_j = det(P_0,P_j) ≠ 0`, blocks in two nested chains on disjoint row sets at every prime,
> `|log n_T - w| <= eta_0 w` for the 119 visible `T` (`2 <= |T| <= 6`), and
> `Phi(D) ∈ M_{0,8}(Q) \ Z`, satisfies
> `H_Sigma(D) = sum_j log|Y_j| + sum_{i<j} log|t_ij| >= c w - O(1)`.

**Normalisations, checked against `ptolemy.md` §1 and §6.2 Steps 0-5.**

* `k = 8`, `m = 7`, anchor `0`. Step 0 gives `w = W/128`; in `(L_8)` `w` is any real for which
  the balance holds.
* Step 1 gives balance for **every** nonempty `T` with `eta = 2 sqrt2 * 255/sqrt(M-7)` and
  `H_Sigma <= 3584 w/(M-1)`. `(L_8)` only needs the 119 visible blocks: the 7 singletons occur in
  no minor and the top block can be inflated by `diag(n,1)` without moving `Phi`
  (`L8.md` Thm B proof; `census.md` §1.1). So "balanced" means: on the 119 visible blocks.
* Step 4 replaced by `(L_8)`: `c w - C <= H_Sigma <= 3584 w/(M-1)` bounds `M` once
  `eta(M) <= eta_0`, `M >= M_Z` (Lemma 6.1) and `w >= 2C/c`. Only **one** fixed `eta_0 > 0` is
  needed, since `eta(M) -> 0`; and `Z` must be independent of `R`, `N`, the primes and the
  template (it is, by the quantifier order).
* `H_Sigma` has `7 + 21` terms, each `log|.| >= 0`.

### 1.2 The nested-chain hypothesis is the caterpillar condition (anchor-free)

At a prime `p` the visible blocks divisible by `p` are edges of the `p`-adic tree of `Phi`
(up to the residue correction of `ptolemy.md` Thm B(iv)). Call a set of pairwise compatible
splits of `{0..7}` a **caterpillar** if one side of each split can be chosen so that the chosen
sides form a single chain under inclusion; equivalently, all internal edges of the tree lie on
one path.

**Lemma 1.1 (caterpillar lemma).** For a compatible split system `Sigma` on `{0, ..., k-1}` the
following are equivalent:
(a) for the anchor `0`, the visible sides (sides not containing `0`) form a union of two chains
whose unions are disjoint (the item-532 condition);
(b) (a) holds for every choice of anchor;
(c) `Sigma` is a caterpillar.

*Proof.* (c) ⇒ (a) for any anchor `a`: let `v` be the spine vertex at which the pendant edge of
`a` is attached (or the nearest spine vertex). The spine edges on either side of `v` have sides
away from `a` that are nested increasingly towards `v`; the two families of leaves beyond `v`
on the two sides are disjoint. (a) ⇒ (c): if `T_1 ⊂ ... ⊂ T_r` is a chain of sides not
containing `0`, the edge of `T_{i+1}` lies on the path from the edge of `T_i` to the leaf `0`;
so each chain lies on a path towards `0`. Two such paths can share an edge only if its side
contains `T_r ∪ S_s` and lies in one chain, impossible since the two unions are disjoint and
every side of a chain lies in its union. So the two paths meet only at a vertex, and their union
is a path containing every internal edge. (b) ⇒ (a) is trivial, and (c) ⇒ (b) is the first
step. `[]`

*Data* (`theory_checks/caterpillar.py`): all compatible split systems (all trees) on
`k = 5, 6, 7, 8` leaves — 26, 236, 2752, 39208 of them, the known counts — were tested; (a),
(b), (c) agree in every case; 26, 221, 2227, 25698 are caterpillars.

So `census.md`'s "admissible for some anchor and labelling" is the same as "caterpillar at every
prime", a condition on `Phi` alone. A caterpillar on 8 leaves has at most 5 internal edges,
which is `census.md`'s bound of `k - 3` visible blocks per prime.

### 1.3 The intrinsic form (from `census.md` Prop 1.3; recalled because §§2-3 use it)

By anchoring (`ptolemy.md` Lemma 5.2; `census.md` Lemma 1.1) every `Phi ∈ M_{0,8}(Q)` is
represented by a primitive `x = (x_1, ..., x_7) ∈ Z^7` modulo translation (anchor `0` at
`infinity`) with

```text
|x_i - x_j| = prod_{T visible, T ∋ i,j} d_T(Phi),      Y_j = 1, t_ij = ±1,          (1.1)
```

where `d_T = prod_p p^{l_T(p)}` are the contact integers and `ell_T = log d_T`. With the
Lipschitz bound `|log n_T - ell_T| <= 2 H_Sigma` this gives (`census.md` Prop 1.3):

* `(L_8^c)(eta_1, Z)` — *no `Phi ∉ Z` has `|ell_T - w| <= eta_1 w` for all 119 `T` with
  `w >= w_1`* — implies `(L_8)` with `eta_0 = eta_1/2`, `c = eta_1/4`;
* `(L_8)` implies `(L_8^c)(eta_0, Z)` for caterpillar `Phi`, in particular for every `Phi`
  whose contact integers are pairwise coprime.

**Consequence used below (rate dichotomy).** Replacing `c w` in `(L_8)` by any `f(w) -> infinity`
does not weaken it in substance: the `f`-version applied to `H_Sigma = 0` data already gives
`(L_8^c)(eta_0, Z)` for caterpillar points, which is the whole existence question. So a
"growth" version of `(L_8)` (`c` replaced by a slowly growing function) is no easier than
`(L_8)` itself; the only freedom is in the constants. Everything below is about the set

```text
Bal(eta) = { Phi ∈ M_{0,8}(Q) : |ell_T(Phi) - w| <= eta w for all 119 visible T, some w },
```

restricted when stated to caterpillar points (`Bal^cat`) or to points with pairwise coprime
contact integers (`Bal^pc ⊂ Bal^cat`).

## 2. Proof lane

### 2.1 What the joint system contains, and which of its structures survive in the worst case

**Proposition 2.1 (joint content).** Let `D` be an `(L_8)`-datum with point `Phi`.

1. *(Plücker, jointly.)* The 70 relations `(3.1)` of `ptolemy.md`, taken together, are
   equivalent to the existence of the point `Phi ∈ M_{0,8}(Q)` with the given minors
   (`ptolemy.md` Thm B(i),(ii): they generate the Plücker ideal, whose variety is the rational
   cone over `G(2,8)`). After cancelling blocks they say nothing about `D` beyond `Phi` and the
   divisibilities `prod_{T∋i,j} n_T | det(P_i,P_j)`.
2. *(Divisibility and valuations, prime by prime.)* At every prime, `v_p` of all minors is
   determined by the `p`-adic tree of `Phi` up to `2 max v_p(t)` (Thm B(iv)). Summed over
   primes, the blocks are the contact integers up to `e^{±2 H_Sigma}`.
3. *(Anchoring.)* Every `Phi` admits an exact representation with `H_Sigma = 0`, `n_T = d_T`
   (§1.3). So the choice of anchor and the residues `Y_j, t_ij` carry no information about `Phi`.
4. *(Nested chains.)* By Lemma 1.1 the chain condition is "caterpillar at every prime". On
   `Bal^pc` (pairwise coprime contact integers) every prime divides exactly one block, the tree
   at `p` is a single edge, and the chain condition, the per-prime gcd bookkeeping and the
   residue corrections are all vacuous.

Hence any argument that uses the 70 relations jointly, the divisibility of minors, the
valuation bookkeeping at every prime of the nested chains, the anchoring and the balance is an
argument about the set `Bal^cat(eta)`, and any proof of `(L_8)` must in particular prove

```text
(L_8^pc):  there are eta_1 > 0, Z, w_1 with  Bal^pc(eta_1) ∩ {w >= w_1} ⊂ Z,
```

a statement in which the only structure left is: a rational point of `M_{0,8}` whose 119 finite
boundary heights are pairwise coprime integers of size `e^{w(1 ± eta_1)}`.

*Proof.* (1)-(3) are the cited theorems and §1.3. (4): pairwise coprime contact integers put at
each prime exactly one edge in the tree; a single edge is a chain, its `l_T(p) = v_p(d_T)`, and
the item-532 two-chain bookkeeping reduces to `v_p(n_T) = v_p(d_T)`. `(L_8) ⇒ (L_8^pc)` is
`census.md` Prop 1.3(2). `[]`

**Answer to the judge's question "which structure of the joint system is genuinely new?"** In
the worst case none of them: the two nested chains per prime and the anchoring are vacuous on
`Bal^pc`, which any proof must cover, and the joint Plücker system *is* the point `Phi`. The
only global input left is the balanced class `beta` and the sign `K.beta = 10 > 0`. So the
proof lane has to establish a positive multiple of `w` from "balanced rational point of
`M_{0,8}`" alone. §2.2 makes that statement as concrete as possible; §2.3 lists what can be
proved; §2.4 examines the one unconditional-flavoured mechanism not covered by any barrier in
the record (abc applied to auxiliary pencils, which is conditional on abc but stronger than
Liouville); §2.5 gives the verdict.

### 2.2 Kapranov/lattice form of (L_8): balanced points are exceptionally short lattice vectors

Use Kapranov's model `Mbar_{0,8} -> P^5 = P(Q^7/Q·1)` for the anchor `0 = infinity`: the point
of §1.3 is `[x] ∈ P^5(Q)`; for `3 <= |T| <= 6` the boundary divisor `D_T` is the exceptional
divisor `E_T` over the linear space `L_T = {x_i = x_j for i, j ∈ T}` (codimension `|T| - 1`), and
for `|T| = 2` it is the strict transform of the hyperplane `{x_i = x_j}`. Then
`H.beta = 31`, `E_T.beta = 1` and

```text
K = -6H + sum_{3<=|T|<=6} (|T| - 2) E_T            (Kapranov; = ptolemy.md's sum_S c_S D_S).
```

For primitive `x` put `log|x| = max_{i<j} log|x_i - x_j|` (the height of `[x]` up to `O(1)`).

**Lemma 2.2 (lattice dictionary).** For primitive `x` with contact integers `d_T`:

1. `Lambda(x) := {y ∈ Z^7 : y_i ≡ y_j mod p^t whenever i, j lie in one residue class of
   x mod p^t of size >= 3}` is a lattice containing `x` and `Z·1`, of index
   `[Z^7 : Lambda(x)] = prod_{3<=|T|<=6} d_T^{|T|-1}`.
2. The finite part of the canonical height is

   ```text
   kappa(x) := sum_{3<=|T|<=6} (|T|-2) ell_T - 6 log|x|
             = log[Z^7 : Lambda(x)] - sum_{3<=|T|<=6} ell_T - 6 log|x|,
   ```

   and `h_K(Phi) = kappa(x) + sum_T (|T|-2) lambda_{E_T,inf}(Phi) + O(1)`.
3. If `Phi ∈ Bal(eta)` at scale `w`, then `log|x| = 31w + O(eta w)`,
   `log[Z^7:Lambda(x)] = 294 w + O(eta w)`, `kappa(x) = 10 w + O(eta w)`, so

   ```text
   |x|  =  [Z^7:Lambda(x)]^{1/6} · e^{-18 w + O(eta w)} :
   ```

   `x` is a vector of `Lambda(x)` (rank 6 modulo `Z·1`) shorter than the Minkowski scale of its own
   congruence lattice by the factor `e^{18w}`.

*Proof.* (1) At a prime `p`, passing from level `t-1` to level `t`, the conditions for a class
`C` of size `>= 3` (already a class at level `t-1` or contained in one) fix the next `p`-adic
digit of `y_i - y_{i_0}` for `i ∈ C \ {i_0}`: `|C| - 1` independent conditions mod `p`. A visible
`T` is a class at exactly `l_T(p)` levels (`census.md` Lemma 1.1), and the top class never occurs
for primitive `x`. So the `p`-part of the index is `p^{sum_{|T|>=3} (|T|-1) l_T(p)}`; multiply
over `p`. (2) `h_H = log|x| + O(1)`; the finite part of `h_{E_T}` is `ell_T` (`ptolemy.md`
Thm B(iv)(1)); insert into `K`, and use (1). (3) For `x` in the integral chart, (1.1) gives
`log|x_i - x_j| = sum_{T∋i,j} ell_T`, a sum of 31 visible terms; `sum_{|T|>=3}(|T|-1)` over the
98 sets is `35·2 + 35·3 + 21·4 + 7·5 = 294`, and `sum (|T|-2) = 196`; `196 - 186 = 10`. `[]`

*Data* (`theory_checks/lattice_index.py`): 60 primitive planted configurations; the index
computed by a Smith normal form equals `prod d_T^{|T|-1}` and (1.1) holds, in all cases.

**Proposition 2.3 (density mechanisms produce only `h_K <= 0`).** Let `e = (e_T)_{|T|>=3}` be
pairwise coprime positive integers, `Lambda_e = {y : e_T | y_i - y_j for i, j ∈ T}`, and let `x ∈ Lambda_e` be
primitive with `|x| >= X/C` for some `X` with `X^6 >= [Z^7 : Lambda_e]` (the range in which
the lattice points of the box `|y| <= X` have positive density, so that sieve, Ekedahl–Poonen
and geometry-of-numbers arguments apply). If the contacts of `x` are not larger than imposed,
`ell_T(x) <= log e_T + delta` for `|T| >= 3`, then

```text
kappa(x) <= - sum_{|T|>=3} log e_T + 196 delta + 6 log C.
```

*Proof.* `kappa = sum (|T|-2) ell_T - 6 log|x| <= sum (|T|-2) log e_T + 196 delta - 6 log X
+ 6 log C`, and `6 log X >= log[Z^7 : Lambda_e] = sum (|T|-1) log e_T` (CRT; each condition
`e_T | y_i - y_j`, `i, j ∈ T`, has index `e_T^{|T|-1}`). For non-coprime `e` the same proof gives
`kappa(x) <= sum (|T|-2) log e_T - log[Z^7 : Lambda_e] + 196 delta + 6 log C`. `[]`

So the mechanism that refuted `(L_5)` (`L8.md` Prop C: a positive proportion of a box) can never
reach `h_K > 0`, let alone `Bal`. By Lemma 2.2(3) a point of `Bal` exists only when its
congruence lattice has a primitive vector `e^{18w}` below the Minkowski scale.

*Heuristic (labelled).* For a random lattice of covolume `V` in `R^6` the probability of a
nonzero vector of length `<= V^{1/6} e^{-18w}` is `≍ e^{-108w}`; the number of balanced
pairwise coprime profiles `(d_T)_{|T|>=3}` at scale `w` is `e^{98w + o(w)}` (98 free integers of
size `e^w`; the 21 pair contacts are then determined by (1.1)). The expected number of balanced
points per unit of `w` is `e^{-10w} = e^{-(K.beta) w}`. This matches `-K.beta`, the exponent of the
other counts in the record (`ptolemy.md` §4 writes `1 - K.beta` per shift class, the extra `1`
being its shift parameter), and shows where the exponent comes from: `-K.beta = 98 - 108` = (profiles) minus (shortness).

The lattice form is a reformulation, not progress: `(L_8^pc)` says that the structured lattices
`Lambda_e` (balanced, pairwise coprime `e`) have no primitive vector `e^{18w}` below their
Minkowski scale, outside a fixed `Z`. The trivial (Liouville) lower bound is: if `y ∈ Lambda_e`
is not in `Z·1`, pick `y_i ≠ y_j`; then `prod_{T ⊇ {i,j}, |T|>=3} e_T` divides `y_i - y_j`
(pairwise coprime `e`), so `|y| >= e^{30w(1 - eta)}`. The balanced `x` sits at `e^{31w}`, only
`e^{w}` above this trivial bound and `e^{18w}` below the Minkowski scale `e^{49w}`. So the
trivial bound cannot exclude it (this is `L8.md` Thm A in lattice language), and nothing
between the trivial bound and Minkowski is known for these lattices.

### 2.3 What can be proved: restricted versions with all quantifiers

**Lemma 2.4 (matching coefficients are coprime for every compatible profile).** Let `D` be an
`(L_8)`-datum and `Q = {a, b, c, d}` a quadruple of labels. The three matching coefficients
`n_A, n_B, n_C` of `ptolemy.md` (3.1) are pairwise coprime.

*Proof.* At a prime `p` the blocks divisible by `p` are pairwise compatible splits (two nested
chains, Lemma 1.1). Their restrictions to `Q` are compatible splits of a 4-set, and two distinct
`2|2` splits of a 4-set are incompatible. So all blocks at `p` meeting `Q` in a `2|2` pattern
induce the same matching, and `p` divides at most one of `n_A, n_B, n_C`. `[]`

(`ptolemy.md` Thm A(4) states this for circle data; the lemma shows it is a property of every
`(L_8)`-datum.)

**Theorem 2.5 (fixed prime set; Ridout).** Let `S_0` be a finite set of primes, `0 <= eta < 1`
and `eps > 0`. There is an (ineffective) constant `C = C(S_0, eta, eps)` such that every datum
`D` as in `(L_8)` — no exceptional set — satisfies

```text
H_Sigma(D)  >=  ( (1 - eps) 16 (1 - eta) w  -  Xi_{S_0}(D)  -  C ) / (2 - eps),
```

where `Xi_{S_0}(D) = sum log(S_0-free part of n_T)` over the 48 visible blocks occurring in the
matching coefficients of one anchored quadruple `{0, i, j, l}`. In particular, if all visible
blocks are `S_0`-units, `H_Sigma >= 8(1 - eta)(1 - eps) w - C'`.

*Proof.* Use the anchored relation `(4.1)` of `ptolemy.md`:
`a N_A + b N_B = c N_C` (signs absorbed), with `a, b, c` products of two of the six distinct
residues `Y_i, Y_j, Y_l, t_jl, t_il, t_ij`, so `|abc| <= e^{H_Sigma}`, and `N_A, N_B, N_C`
products of 16 visible blocks each, pairwise coprime (Lemma 2.4), each `>= e^{16(1-eta)w}`.
Let `g` be the gcd of the three terms; a prime power dividing `g` divides two terms, hence (the
`N`'s being coprime) divides `abc`; so `g <= |abc| <= e^{H_Sigma}`. Put `x, y, z` = the three
terms divided by `g`: coprime, `x + y = z`, `max(|x|,|y|,|z|) >= e^{16(1-eta)w - H_Sigma}`. Their
`S_0`-free part divides `abc` times the `S_0`-free parts of the `N`'s, so it is
`<= e^{H_Sigma + Xi}`. Ridout's `p`-adic Roth theorem for the three points `0, -1, infinity` of
`P^1` at the places `S_0 ∪ {infinity}` (equivalently the three-term Subspace Theorem) gives
`log(S_0-free part of xyz) >= (1 - eps) log max(|x|,|y|,|z|) - C(S_0, eps)`. Combine. `[]`

This is a correct restricted `(L_8)` with explicit `c = 8(1-eta)(1-eps)` and **no** exceptional
set, and it is worthless for the circle problem for exactly the reason recorded in lit18
Props A/B: `C(S_0, eps)` depends on `S_0`, not only on `|S_0|`, and in the application `S_0`
is the set of primes of `N` (at least 24 of them by `census.md` Thm E, and unbounded). No
uniformity in `S_0` is claimed or used. The variant "bounded number `s` of primes, varying" is
**not** provable by the same tools: Evertse's uniform bound `3·7^{3+2s}` counts solutions but
bounds no height, and Baker-type height bounds grow with `max S_0`.

**Restrictions that do not help.** "One-chain profiles" (all blocks at each prime in one chain)
contain `Bal^pc`, so `(L_8)` restricted to them still contains `(L_8^pc)`; "squarefree blocks"
contain the squarefree part of `Bal^pc`, where the heuristic exponent is the same `-10` and where
the record's fake circles live (`round8/weakest.md` Cor 4.3); the restriction to a bounded number of primes is discussed above; replacing `c w`
by `f(w) -> infinity` is no easier (§1.3).

**Other restricted versions already proved this round or before** (cited, not redone):
`(L_8)` holds on every subvariety missing a boundary divisor, every symmetric locus, every
forgetful fibre (`census.md` Prop 7.1, Cors 7.2, 7.3), and on every fixed curve of genus `<= 1`
of class not proportional to `beta` or with `K.C <= 0` (`census.md` Props 7.4, 7.5(4)). In the
function-field setting, no covering family of curves of genus `<= 1` has class with imbalance
`< 5/13` (`census.md` Prop 7.5(2)); for curves of genus `g` through a general point the
deformation argument gives `K.gamma <= 2g - 2`, so a class `d·beta` needs `g >= 5d + 1`.

### 2.4 The one mechanism no barrier in the record covers: abc on auxiliary pencils

The record's abc barrier (`conditional.md` Prop 8.3, `round8/weakest.md` Prop 4.2, refereed)
covers relations whose terms are products of differences with pairwise coprime squarefree blocks
(pair, triangle, Ptolemy, `n`-term relations). Its own "precise scope" paragraph excludes
relations with high-degree cancellation and refers them to the auxiliary-polynomial barrier. That
referral is not a proof: `L8.md` Thm A bars *Liouville* arguments (size versus forced
divisibility), while abc also exploits that a prime counts once in the radical however often it
divides. So "abc + an auxiliary identity" is a gap in the barrier map. Since Vojta implies abc,
closing it would at best give the uniform theorem **conditionally on abc**: a conditional theorem
from a different, more widely studied hypothesis than the two Vojta-conditional theorems of the
record, but not an unconditional one. This section makes the gap
precise and partly closes it.

**Setting.** A *pencil* is `f = s A/C` with `s = ±1` and `A, C ∈ Z[x_1..x_7]` coprime,
homogeneous of the same degree `D`, translation invariant; put `B = C - sA`. In Kapranov's model
let `nu_T(F)` be the order of `F` along `L_T` (`|T| >= 3`) or the exponent of `x_i - x_j`
(`T = {i,j}`), and `mu_T = min(nu_T(A), nu_T(C))` (`= 0` for pairs). Define

```text
L.beta     = 31 D - sum_{|T|>=3} mu_T,
F_red.beta = #{T visible : max(nu_T(A), nu_T(B), nu_T(C)) > mu_T} + sum_Y (31 deg Y - sum_{|T|>=3} nu_T(Y)),
gain(f)    = L.beta - F_red.beta,
```

the last sum over the `Q`-irreducible factors `Y` of `ABC` that are not differences
(interior components). `L` is the class of the fibres of `f` and `F_red` the reduced union of the
fibres over `0, 1, infinity` (on the blow-up of the base locus).

**Proposition 2.6 (abc-gain criterion).** Assume the abc conjecture over `Q`. Let `f` be a
*unit pencil* (`A` and `C` products of differences `x_i - x_j`; `B` arbitrary). If
`gain(f) > 0`, then `(L_8^pc)` holds with `Z = {B = 0}`: there are `eta_1 > 0` and `w_1` such
that no `Phi ∈ Bal^pc(eta_1) \ Z` has `w >= w_1`.

*Proof.* Let `x` be the primitive integral representative (§1.3) of `Phi ∈ Bal^pc(eta)`. By
(1.1), `A(x) = ± prod_T d_T^{nu_T(A)}`, likewise `C(x)`, and (pairwise coprime `d_T`)
`gcd(A(x), C(x)) = prod_T d_T^{mu_T}`. Every `|x_i - x_j|` is `e^{31w + O(eta w)}`, so
`h(f(Phi)) = log max(|A(x)|,|C(x)|) - log gcd = (L.beta) w + O(eta w)`. Write `a, b, c` for
`sA(x), B(x), C(x)` divided by `g = gcd`; `a + b = c`, coprime. A prime of `d_T` divides `abc`
only if `max(nu_T(A),nu_T(B),nu_T(C)) > mu_T`, contributing at most `ell_T <= (1 + eta) w` to
`log rad`. Every other prime divides `Y(x)` for an interior factor `Y` of `B` (the only interior factors,
`A` and `C` being units) beyond its forced part `prod_T d_T^{nu_T(Y)}` (forced because `Y` lies in
the `nu_T(Y)`-th power of the ideal of `L_T` over `Z` and `d_T | x_i - x_j` for `i, j ∈ T`); that cofactor has log at most
`31 deg(Y)(1 + eta) w - sum_{|T|>=3} nu_T(Y)(1 - eta) w`. Hence
`log rad(abc) <= (F_red.beta) w + O(eta w)`. abc with exponent `1 + eps` gives
`(L.beta) w <= (1 + eps)(F_red.beta) w + O(eta w) + C_eps`, impossible for large `w` once
`eps` and `eta` are small compared with `gain(f)`. Points with `B(x) = 0` lie in `Z`. `[]`

For a pencil whose base locus meets the interior (`A, C` with common zeros off the boundary),
`log gcd(A(x), C(x))` contains accidental primes, a height of the base locus that balance does
not control (a Vojta-type gcd problem in its own right). So unit pencils, and pencils whose base
locus lies in codimension-2 boundary strata (where the gcd vanishes on `Bal^pc`), are exactly
the ones for which abc can be used without a further gcd hypothesis.

**Identity (why the gain is the right quantity).** On the blow-up `X~` of the base locus
(inside codimension-2 boundary strata) `f` is a morphism `X~ -> P^1`; with the relative
logarithmic canonical class `K^log_{X~/P^1} = K_{X~} + F_red - f^*(K_{P^1} + [0] + [1] + [inf])`
and `beta~ = pi^* beta` (exceptional degrees `0` on `Bal^pc`),

```text
gain(f)  =  K.beta  -  K^log_{X~/P^1}.beta~  =  10 - K^log_{X~/P^1}.beta~ .
```

(Insert `f^*(K_{P^1} + 3 pts) = L` and `K_{X~}.beta~ = K.beta`.) So abc reaches `(L_8^pc)`
through `f` iff the relative log canonical class of the pencil has `beta`-degree `< 10 = K.beta`.
For curves this is Elkies' Belyi argument (relative term `0`). Campana's theorem that rationally
connected varieties admit no fibration onto a hyperbolic orbifold curve does **not** exclude
`gain > 0`: it constrains the minimal multiplicity of each fibre, while the gain is a
`beta`-weighted average, and a fibre may contain reduced boundary components of `beta`-degree 1
next to highly multiple ones. A covering family of curves of genus `g` and class `d·beta` would
give `gain <= (2g - 2)/d` (Mason on the curves), which is `>= 10` by `K.beta <= 2g - 2`; so no
curve argument can exclude `0 < gain <= 10` either.

**Partial closure (proved and computed).**

1. *Morphisms.* If `f` has no base points, `f: Mbar_{0,8} -> P^1` with Stein factorization
   `h ∘ phi`. By Bruno–Mella (*The automorphism group of Mbar_{0,n}*, JEMS 15 (2013); arXiv
   1006.0987: a dominant morphism with connected fibres `Mbar_{0,n} -> Mbar_{0,r}` factors
   through a forgetful map — read from the abstract and a search summary only, the text was not
   accessible here), `phi` is a forgetful map to `Mbar_{0,4}` up to `Aut(P^1)`. For
   `f = h ∘ chi_Q`: `L.beta = 16 deg h`, the reduced fibres over `0, 1, inf` are
   `chi_Q^*(h^{-1}{0,1,inf})_red` of `beta`-degree `16 #h^{-1}{0,1,inf} >= 16(deg h + 2)`
   (Riemann–Hurwitz on `P^1`). So `gain <= -32`.
2. *Unit pencils of low degree (Data).* `theory_checks/abc_gain.py` enumerates all unit pencils
   with `deg A = deg C = D <= 3` (all `A` up to `S_7`, all coprime `C`, both signs), factors `B`
   over `Q` (python-flint) and computes `nu_T` of the interior factors exactly (two independent
   random specializations, exact interpolation; generic-specialization values, so Data, not
   proof). Results: `D = 1`: 40 pencils, gains `-32` (20) and `-48` (20). `D = 2`: 1180 pencils,
   maximum `-32`. `D = 3`: 19800 pencils, maximum `-32` (attained by 20 pencils, e.g. the
   composition `((x_1-x_2)/(x_1-x_3))^3` of a cross-ratio with `z^3`), next values `-40, -44, -48`. Sanity:
   Ptolemy and anchored Ptolemy both give `-32`. `D = 4` (`abc_gain_D4.py`, all 22 `A`-orbits):
   256920 pencils, maximum `-32` (24 pencils), next `-40, -48`. The `D <= 3` run was repeated
   with an independent computation of `nu_T` (composition into `Z[t]` instead of interpolation),
   with identical output (`abc_gain_D3_v2.out`).
3. *Linear pencils with base locus in a boundary stratum,* e.g.
   `f = (x_1 - x_2)/(x_3 - x_4)` (base locus `D_12 ∩ D_34`): `L.beta = 24`,
   `F_red.beta = 48 + 24` (the interior hyperplane `x_1 - x_2 = x_3 - x_4` has `beta`-degree
   24), `gain = -48` (by hand; it is one of the 20 `D = 1` pencils with gain `-48`).

4. *Multiplicity form and a Mason bound.* Since the three fibres are linearly equivalent,
   `3L = F` (total fibre), so `gain(f) = sum_Y ((m_Y - 3)/3) (Y.beta)` over the components `Y`
   of the three fibres with multiplicities `m_Y`. A positive gain needs components of
   multiplicity `>= 4` to dominate in `beta`-degree in **all three** fibres together; the
   cross-ratio has 48 components of multiplicity 1 (`48·(-2/3) = -32`). For a unit pencil, Mason
   on a general line of `P^5` (it meets only the 21 hyperplanes `x_i = x_j`) gives
   `sum_{pairs in B} (b_ij - 1) + sum_{Y interior} (m_Y - 1) deg Y_red <= 19 - #(pair factors of B)`:
   interior components of the fibre over 1 can be multiple only in bounded degree. The
   exceptional divisors `E_T` (`|T| >= 3`) are invisible to lines, which is why no line or
   rational-curve argument settles the sign of the gain on `beta` (`K.beta > 0`, while
   `K < 0` on free rational curves).

So the abc lane is open only through unit pencils of degree `>= 5` (or non-unit pencils with boundary
base locus) whose three fibres are, in the `beta`-weighted sense, more than triply multiple.
No such pencil is known, and the low-degree data show the gain stuck at `-32`, the value of the
cross-ratios.

### 2.5 Verdict of the proof lane

No proof of `(L_8)`, of `(L_8^pc)`, or of any `f(w) -> infinity` version was found; by §1.3 the
last is not easier than `(L_8)`. The joint system yields no positive multiple of `w` from any
structure other than `K.beta = 10`, because in the pairwise coprime subcase (which every proof
must cover) the nested chains, the per-prime bookkeeping and the anchoring are vacuous
(Prop 2.1). What is proved: the lattice form (Lemma 2.2) with the exponent count `98 - 108`;
the exclusion of density mechanisms (Prop 2.3); a fixed-`S_0` restricted theorem with explicit
`c = 8` and no exceptional set (Thm 2.5), useless for the application; and, conditionally on
abc, the exact criterion of Prop 2.6 for the one mechanism the barrier map does not cover, with
morphisms excluded and unit pencils of degree `<= 4` excluded by computation.

## 3. Refutation lane: dimension counts and Vojta

### 3.1 What a refutation needs (quantifiers)

`(L_8)` fails iff for every `eta_0 > 0` and every proper closed `Z` there are data off `Z` with
`eta <= eta_0`, `w` arbitrarily large and `H_Sigma = o(w)`; by §1.3 this is the same as:
**for every `eta > 0` and every `w_1`, the set `Bal(eta) ∩ {w >= w_1}` is Zariski dense in
`Mbar_{0,8}`** (caterpillar points suffice, and give data with `H_Sigma = 0`; conversely data
with `H_Sigma = o(w)` give points of `Bal(eta + o(1))` by the Lipschitz bound).

### 3.2 Expected dimensions (the analogue of the `+8` at `k = 7`)

In Kapranov's model a curve of class `d·beta` is a curve of degree `31d` in `P^5` meeting each
`L_T` (`3 <= |T| <= 6`, codimension `|T| - 1`) in `d` points. Meeting a codimension-`c` linear
space once costs `c - 1` conditions, so the incidences cost `d · sum_{|T|>=3} (|T| - 2) = 196 d`.

| curves of class `d·beta` | parameters | conditions | expected dimension |
|---|---|---|---|
| genus 0, `k = 8` | `6·31d + 2 = 186d + 2` | `196 d` | `-10d + 2` |
| genus 1, `k = 8` | `186 d` | `196 d` | `-10 d` |
| genus 0, `k = 7` (check) | `5·15d + 1` | `68 d` | `7d + 1` (`= 8` at `d = 1`, as in `ptolemy.md`) |

In every case the expected dimension is `-K.(d beta) + (dim X - 3)(1 - g)`, and the parameter
minus condition count `186d - 196d` is literally `-K.(d·beta)` written in Kapranov's
coordinates (`K = -6H + sum (|T|-2) E_T`). For point counts the same bookkeeping is Lemma 2.2(3)
and its heuristic: `98w` free contacts against `108w` of shortness, net `-10w`.

### 3.3 Why no family-type mechanism can work, and what is left

1. **Covering families of curves are excluded** (proved): a covering family of genus-`g`
   curves has `K.gamma <= 2g - 2`; a class `d·beta` needs `g >= 5d + 1`, and curves of genus
   `>= 2` have finitely many rational points (Faltings). Genus `<= 1` covering families have
   imbalance `>= 5/13` (`census.md` Prop 7.5). This kills the mechanisms of `(L_5)` (rational
   covering families) and `(L_6)` (`L8.md` Thm B: a covering genus-one family, `K.beta_6 = -8`).
2. **Density (box/sieve) mechanisms are excluded** (proved, Prop 2.3): they produce only
   `h_K <= -sum_{|T|>=3} ell_T`. This kills the mechanism of `L8.md` Prop C.
3. **Lifts of the `(L_6)` families** are excluded at the curve level (`lift.md` Thm 4.2).
4. **What remains** for a refutation, each alternative being a counterexample to `V_8`
   (Vojta, `D = 0`, on `Mbar_{0,8}`), because `V_8 ⇒ (L_8)` (`ptolemy.md` Thm D Steps 3-4,
   refereed):
   * infinitely many **rigid** curves of genus `<= 1` (positive rank in genus 1), of classes
     `gamma_n` with `eta(gamma_n) -> 0` (so `K.gamma_n > 0`), whose union is Zariski dense; or
   * a Zariski-dense set of **sporadic** points lying on no such curve, each an exceptionally
     short vector (factor `e^{18w}` below Minkowski) of its own congruence lattice.

   Rigid rational curves with `K.gamma > 0` do exist in `Mbar_{0,8}`, but only in the
   boundary as far as I know: the `F`-curve of type `(2,2,2,2)` has `K.F = 2`, type
   `(3,2,2,1)` has `K.F = 1`, types `(3,3,1,1)`, `(4,2,1,1)` have `0` and `(5,1,1,1)` has `-1`
   (`theory_checks/fcurves.py`). They lie in the boundary and carry no points of
   `M_{0,8}`. I did not survey rigid curves meeting the interior (e.g. the hypertree-type
   constructions of Castravet–Tevelev) for their `K`-degree; none of class near the `beta`-ray
   is known in the record, and none was found (`census.md` §6).

### 3.4 Can (L_8) fail while Vojta holds? The formulation audit

No. `V_8 ⇒ (L_8)` is a refereed theorem, with `Z = Z_{1/43}` and `eta_0 = 1/381`; it uses only
Lemma 6.2 (heights of a balanced Boolean point, which absorbs the residues through the Lipschitz
bound and the archimedean place through balance). Item by item:

* **`O(1)`.** It may depend on `(c, eta_0, Z)`; Theorem D uses it only through `w >= C'`, with
  `w >= (M-1) lambda/32`. Harmless.
* **`eta_0`.** Any fixed `eta_0 > 0` suffices (Step 1 gives `eta(M) -> 0`). `(L_8)` for small
  `eta_0` is the weakest form; Vojta gives it for `eta_0 <= 1/381`. Imbalanced points
  (`eta >= 5/13`) can lie on covering families (`census.md` Prop 7.5) and are irrelevant.
* **Exceptional set.** `Z` is fixed before the data; Lemma 6.1 avoids it at cost `O(M^7)`
  uniformly in the cluster. No fixed-`S` Subspace exceptional set enters.
* **Residues, anchoring, invisible blocks, nested chains.** All handled exactly (§1.3,
  Lemma 1.1, Prop 2.1). The chain hypothesis only shrinks the hypothesis class.
* **Balance on visible versus all blocks.** Step 1 gives all; `(L_8)` asks only visible; the
  invisible ones are free. No gap.

So a refutation of `(L_8)` is exactly a refutation of the `D = 0` Vojta inequality at balanced
points of `Mbar_{0,8}`, in its weakest form (`eps = 1/43`, any `eta_0`), and `(L_8)` cannot be
"mis-stated" in a way that a refutation could exploit. The refutation lane therefore ends at
the conjecture itself: it would need one of the two configurations of §3.3(4), and both are
predicted not to exist by the same heuristic (`-10w`) that correctly predicted existence at
`k = 5, 6` (positive exponents `6, 9`) and is consistent with all `k = 8` data (`census.md`).

## 4. Computations

All scripts are in `theory_checks/`; run with `python3` (each well under the 3.5 GB limit).

| script | content | type | result |
|---|---|---|---|
| `caterpillar.py` | Lemma 1.1 on all trees with 5-8 leaves | exact | 26/236/2752/39208 trees, 0 mismatches |
| `lattice_index.py` | Lemma 2.2(1) index via Smith normal form; identity (1.1) | exact | 60/60 configurations pass |
| `fcurves.py` | `K.F`, `B.F`, `(K+B).F` for the 5 F-curve types of `Mbar_{0,8}` | exact | `K.F = -1, 0, 0, 1, 2`; `(K+B).F = 1` for all |
| `abc_gain.py` | gain of all unit pencils, `D <= 3` (A up to `S_7`) | exact arithmetic; `nu_T` of interior factors by generic specialization (Data) | max gain `-32` for `D = 1, 2, 3` (40, 1180, 19800 pencils) |
| `abc_gain_D4.py` | same, `D = 4`, all 22 `A`-orbits | as above | 256920 pencils, max gain `-32` |

Outputs: `abc_gain_D3.out`, `abc_gain_D3_v2.out` (independent `nu_T` method), `abc_gain_D4.out`,
`fcurves.out`. Driver: `run_abc.sh`.

## 5. Precise residual

1. **Which lane got further.** Neither decides `(L_8)`. The refutation lane ends in a clean
   statement: a refutation of `(L_8)` is a refutation of Vojta's `D = 0` inequality at balanced
   points of `Mbar_{0,8}` (with `eps = 1/43`), and must come from infinitely many rigid curves of
   genus `<= 1` with `K.gamma > 0` that are Zariski dense, or from a Zariski-dense set of sporadic
   exceptionally short lattice vectors (§3.3); all family and density mechanisms are excluded by
   proofs (§3.3(1)-(3)). The proof lane ends at `(L_8^pc)` in lattice form (§2.2), with no
   structure beyond the balanced class.
2. **The single sharpest open statement** (equivalent in substance to `(L_8)`, by §1.3 and
   Prop 2.1): *there are `eta_1 > 0` and a proper closed `Z` such that for pairwise coprime
   integers `e_T ∈ [e^{(1-eta_1) w}, e^{(1+eta_1) w}]` (`T ⊆ [7]`, `3 <= |T| <= 6`), the
   lattice `Lambda_e = {y ∈ Z^7 : e_T | y_i - y_j for i, j ∈ T}` has no primitive vector `y ∉ Z·1`,
   whose point `[y] ∈ M_{0,8}` lies outside `Z`, with `|y| <= e^{31 w (1 + eta_1)}` whose pair cofactors
   `(y_i - y_j)/prod_{T ⊇ {i,j}} e_T` are pairwise coprime of size `e^{w(1 ± eta_1)}`, once
   `w >= w_1`* (`(L_8^pc)` in lattice form; it implies the uniform theorem only together
   with the chain cases, i.e. `(L_8)` itself).
3. **Conditional side result to pursue (new gap in the barrier map).** abc + a pencil of
   positive gain (Prop 2.6) would give `(L_8^pc)`; with an extension of the gcd bookkeeping to
   nested chains it would give the uniform theorem under abc. Excluded: morphisms (Bruno–Mella,
   cited from the abstract), unit pencils of degree `<= 4` (data, 277940 pencils, maximum gain `-32`). The
   criterion is `K^log_{X~/P^1}.beta~ < 10`; equivalently the `beta`-weighted average
   multiplicity of the three special fibres must exceed 3. I do not know whether such a pencil
   exists; Campana's theory and curve arguments do not exclude it (§2.4).
4. **Not claimed.** No growth improvement, no uniform bound, no unconditional case of `(L_8)`
   beyond fixed `S_0` (Thm 2.5) and the special loci of `census.md` §7.


