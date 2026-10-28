# Referee report on round18/lift.md (the lifting test for (L_8))

Scripts: `round18/referee_lift_checks/` (`ref_beta7.py` -> `ref_beta7.json`, `ref_symb.py`, `ref_branch.py`,
`ref_phiN.py` -> `ref_phiN.json`, `ref_split.py`, `ref_lattice.py`, `ref_skew.py`, `NOTES.md`). All are written
from scratch: they do not import the author's code. Arithmetic is exact (`Fraction`, sympy, flint integer
factorisation); logarithms are used only to report sizes.

## 0. Verdict

**The main theorem survives. The headline verdict is overstated, and Section 5's data is almost vacuous.**

* **Theorem 4.2 is correct.** I re-derived it by a different and shorter route (Section 3) and confirmed
  it at all nine parameter points, with exact arithmetic over Q. The coplanarity identity
  `q* = B_bc` holds symbolically in `x_1..x_4`.
* **Lemma 1.1, Lemma 4.1, Theorem 4.2(b), Corollary 4.3 and the genus computations are correct.**
* **The counts of Section 3 are arithmetically correct.** They are heuristics, and the write-up labels
  them as such.
* **No constraint of the project owner is violated.** The write-up makes no growth or uniformity
  claim.
* **Overclaim.** "The `E_t` mechanism does not lift, already at the first step" is stronger than what
  is proved. Only curves of class exactly `beta_7` are excluded. Genus-one lifts of class `c beta_7`
  with `c >= 2` are still the same mechanism, and they are open (the write-up's own table says so).
* **Section 5 data is nearly vacuous.**
  * For `N = 3..7`, the author's configurations `Phi_N` are not 1/5-balanced themselves. By Lemma 1.1,
    no 1/5-balanced 7-lift can exist, whatever the factorisations of the contacts.
  * At `eta = 1/5` the lattice ("skew") test requires a negative skew, so it tests nothing.
  * The only informative datum is the factorisation obstruction at `N = 8`.

## 1. Is the object exactly (L_8)?

* **Formulation.** The write-up's statement in Section 1 matches `ptolemy.md` §7 and `L8.md` §1:
  * anchor `P_0 = (1,0)`;
  * `det(P_i,P_j) = t_ij prod_{T ∋ i,j} n_T`;
  * two nested chains at each prime;
  * `H_Sigma = sum log|Y_j| + sum log|t_ij|`;
  * a proper closed exceptional set `Z`.
* **Balance (small deviation).** The write-up imposes balance on boundary blocks only
  (`2 <= |T| <= 6`). `ptolemy.md` §7 writes `|log n_T - w| <= eta_0 w` without naming the range of `T`.
  * The boundary-only reading has a weaker hypothesis, so it is the **stronger** lemma.
  * It still feeds Steps 0-5 of §6.2, because Step 1 gives balance for every nonempty `T`.
  * The word "verbatim" is therefore slightly inaccurate, but nothing depends on it, because the
    write-up neither proves nor refutes `(L_8)`.
* **Refutation standard.** The write-up claims no refutation, so no Zariski-density check is needed.
  Its Section 6.2 correctly says that a positive lift would refute `(L_8)` only if the resulting
  configurations were dense in `Mbar_{0,8}`.

## 2. Lemma 1.1 (sub-configurations)

Re-derived; the lemma is correct.

* Take `S ⊆ K \ {0}` with `2 <= |S| <= k'-2`. Each `T` with `T ∩ K = S` is `S ∪ B`, where `B` is a set
  of forgotten labels. Then `2 <= |T| <= k-2`, and there are `2^{k-k'}` such `T`.
* The residues and the `Y_j` are a sub-family, so `H_Sigma` can only decrease.
* Chains restrict to chains: `T ⊆ T'` implies `T ∩ K ⊆ T' ∩ K`. Disjoint row sets stay disjoint.
* The consequences stated (an 8-lift at scale `w` gives 7-lifts at scale `2w` and forces `Phi_6` to be
  balanced at scale `4w` with the same `eta`) follow.

**Caveat used in Section 5 below.** The lemma also implies that `Phi_6` must itself be `eta`-balanced
for an `eta`-balanced lift to exist. The write-up never applies this to its own data.

## 3. Lemma 4.1 and Theorem 4.2

**Lemma 4.1.** Re-derived; correct.

* (1) `pi^* D'_A = D_A + D_{A∪6}` gives degree 2 over `E_t`. The zero fibre of `x_6 - x_j` consists of
  the 8 divisors `D_{{6,Aj} ∪ S}`.
* (2) `Pic(E x P^1) = Pic(E) ⊕ Z`.
* (3) Stabilisation of `(6, A1, A2, A3)` on `D_{A∪6}` over `B_A`. Two cases:
  * when `A` contains no anchor, the node position `x_4` on the anchor component;
  * otherwise the boundary point `{6, Aj} | rest`.
* Reduction (4.1) is correct. Each `sigma_j` lies in `H^0(L(-D_j))` with `deg = 1`, which is a line.

**Theorem 4.2(a): independent proof by a different route** (`ref_beta7.py`).

1. Let `V` be the space of forms `F` in `H^0(O_E(3H)) ⊗ H^0(O_{P^1}(2))` (36-dimensional) that vanish
   on `aux_j x P^1`. Let `W` be the subspace of `V` vanishing at the 25 points `(B_A, v_A)`.
2. **Computed exactly at all nine parameter points: `dim W = 2`.**
3. Let `Pi` be the plane through `aux_1, aux_2, aux_3`. Then `Pi*sigma_0` and `Pi*sigma_inf` lie in
   `W` and are independent, so `W = Pi * Pencil`. The fourth point of `Pi ∩ E` is `B_bc` at all nine
   points.
4. **Coplanarity is symbolic.** `det(aux_1, aux_2, aux_3, B_bc) = 0` holds identically in
   `x_1..x_4` (`ref_symb.py`).
5. Every `L ∈ Pic^8` is `L_q = O(3H - aux - q)` for a unique `q`. For `q` outside the 28 special points,
   the solution space is `{F ∈ W : F(q, .) ≡ 0}`. An element `Pi*sigma_mu` lies in it iff `Pi(q) = 0`
   or `sigma_mu(q, .) ≡ 0`.
   * `Pi(q) = 0` happens only at the special points `aux_j` and `B_bc`.
   * Expanding `sigma_mu(q, .) ≡ 0` coefficient by coefficient in `v` gives
     `D(t) := e_3^2 - e_1 e_2 e_3 + e_1^2 e_4 = 0`. That is exactly the case of a constant pencil member.
   * `D ≠ 0` at all nine points; for example `D(t_0) = 1440`.
6. **The 28 special `q` were computed directly**, with order-2 vanishing along the tangent of the
   space quartic. The solution dimension is 0 except at `q = B_bc`, where it is 2. This holds at all
   nine points.

This confirms (a) for `L` over `Qbar` at the nine points. It also confirms the identification
`L_{q*} = 2H`. The write-up's openness argument is valid: rank conditions on matrices whose entries are
regular in `t ∈ U`, and the absence of common zeros of a family of binary forms, are open conditions.

**Recommendation.** Replace the minors-and-norms proof by steps 1-6. They are shorter, and they make
`U_0` explicit up to the rank conditions. Note that `U_0 ⊆ {D ≠ 0}`: when `D = 0`, statement (a)
fails for every `q` with `T(q) = -e_3/e_1`.

**Theorem 4.2(b) and (c).** Correct.

* The irreducibility criterion holds: the characters `T`, `T - s` and `T - t'` generate all six branch
  pairs in `{0, inf, s, t'}`.
* `2g - 2 = 8n` gives `g = 9` or `g = 5`.
* At all nine points, `{s, t'}` is not a branch pair: the resultant in `lam` of `Delta_lam(s)` and
  `Delta_lam(t')` is nonzero (`ref_branch.py`).
* The degenerate `lam` (where `P_lam` and `Q_lam` share a root) are finite: three or four values at
  each point, for example `{1/2, 1/6, 1/30}` at `t_0`.
* Small point on (c): "smooth curve" is asserted, not proved. Only the genus of the normalisation and
  the class are used, so this does not matter.

**Model constants.** The constants of Section 2 (`kb`, `s`, `kc`, `t'`, `U_j`, `W_j` at `t_0`) were
reproduced exactly.

## 4. Corollary 4.3

Correct.

* `pi_*[C] = 2 beta_7`, so `deg(C -> C_n) ∈ {1, 2}`.
* In the degree-2 case, `C_n` is a pencil lift, and Riemann-Hurwitz gives `g(C) >= g(C_n) >= 5`.
* The pencil `beta_8` lift has genus 33: `2g - 2 = 32(-2) + 8 * 16`. It has `93 + 26 = 119` contacts.

## 5. Counts (Section 3)

All the arithmetic re-derived; the numbers are correct:

* `16 - 25 = -9`, `32 - 50 = -18` and `64 - 106 = -42`;
* `24 - 25 = -1` for the linear relaxation;
* the genus-one and genus-zero expected dimensions `-K.beta = 8, 7, -10` and `8, 8, -8`;
* the integer exponents `e^{-9w}` and `e^{-42w}`;
* the bookkeeping `32 - 42 = -10 = -K.beta_8`.

They are heuristics and are labelled as such. Section 3.3 itself shows that structure beats them twice.

## 6. Section 5 data

**Reproduced exactly with independent code** (`ref_phiN.py`, `ref_split.py`, `ref_lattice.py`):

* the author's configurations (`w_6` = 14.890, 30.202, 52.729, 80.638, 113.053, 150.652);
* the unsplittable counts 16, 11, 8, 6, 8, 9;
* the mover-`a` column: `lambda_1 = (u_a, v_a)`, skew 6.95 to 6.89 `w`;
* random-split maximal skews of the same order (0.092, 0.063, 0.047 at `N = 6, 7, 8` with 1000 samples).

**Problems.**

1. **Nearly vacuous.** The balance `eta_min = (max - min)/(max + min)` of the log Sigma-free contacts
   of `Phi_N` itself is:

   | `N` | 3 | 4 | 5 | 6 | 7 | 8 |
   |---|---|---|---|---|---|---|
   | `eta_min` (author's configuration) | 1.00 | 0.46 | 0.44 | 0.22 | 0.26 | 0.157 |
   | `eta_min`, all 8 points over `T(R_N)` | | | | | | `>= 0.147` |

   By Lemma 1.1, a 1/5-balanced lift forces `eta_min(Phi_N) <= 1/5`. So for `N <= 7` no such lift exists
   for this trivial reason, and the "factorisation" obstruction adds information only at `N = 8`.
2. **The lattice test is void at the stated `eta`.** With `k` kept contacts, an `eta`-balanced lift needs
   skew `>= k(1-eta)/2 - 8(1+eta)`.
   * At `eta = 1/5` this is `-2.0`, `-2.8` and `-3.2` for `k = 19, 17, 16`.
   * The bound `k/2 - 8` quoted in the write-up is the `eta -> 0` value.
   * The test bites only when `eta` is below about 0.09. At that balance every `Phi_N` with `N <= 8`
     already fails to be balanced.
   * So the statements "a lift would need about `0.5w` to `1.5w`" and "random lattices" do not test
     lifts at these `N`.
3. **The window `[d^{1/4}, d^{3/4}]` assumes Sigma-free lift blocks.** If blocks are allowed bounded
   Sigma content, `log e_A` shifts by `O(1)`. That is not negligible at `w = 7` to `75`.
4. **The configurations are slightly misdescribed.** `intlift.py` takes the first of the 8 points of
   `E_t` over `T(R_N)`, without checking `W`. It is therefore a point over `±R_N`, not exactly
   `Psi(Q_1 + N(Q_2 - Q_1))`. This is harmless.

The data is correctly labelled as evidence, but its reading in Section 5 and in §0 item 4 overstates
what it shows.

## 7. Constraint check

* No fixed-`S` Subspace argument is used.
* No Gaussian or denominator costs are discarded; the setting is form-blind and needs none.
* No drifting phase is set to zero.
* No `t`-dependent threshold is called uniform. Section 6.3 says this explicitly.
* Evidence and proof are kept separate, except for the overstated readings noted in Sections 0 and 6.

## 8. Corrected statement

For `t` in a nonempty Zariski-open `U_0 ⊆ U ∩ {e_3^2 - e_1 e_2 e_3 + e_1^2 e_4 ≠ 0}`, which contains
the nine listed parameter points (verified independently, exact over Q):

1. The sections of `L ⊠ O(2)` on `E_t x P^1` that vanish at the 25 points `(B_A, v_A)` are zero for
   every `L ∈ Pic^8(E_t ⊗ Qbar)` except `L = O(2H)`. For `L = O(2H)` they form exactly the pencil
   `{h^2 P_lam(v) - a^2 Q_lam(v)}`.
2. Hence every curve of class `beta_7` in `Mbar_{0,7}` that maps onto `E_t` is a pencil lift `C'_lam`,
   of geometric genus 9, or 5 for at most two `lam`. **No class-`beta_7` lift of genus `<= 4` exists.**
3. A class-`beta_8` lift of genus `<= 4` must have both 7-point projections birational, of class
   `2 beta_7`.

What remains open:

* genus-one lifts of class `c beta_7` with `c >= 2`;
* class-`beta_8` lifts whose 7-point projections have class `2 beta_7`;
* all point-level (integer) lifts.

The "`E_t` mechanism does not lift" conclusion holds only for class `beta_7`. The integer data on
`Phi_N(t_0)`, `N <= 8`, is consistent with this, but it is non-trivial only at `N = 8`, and only
through the factorisation obstruction (9 unsplittable contacts). Nothing here bears on `(L_8)` for
points, on the growth rate, or on the uniform bound.
