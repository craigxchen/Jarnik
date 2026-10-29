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
