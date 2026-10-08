# Growth-rate and uniform-bound attempts (2026-10-08)

Every write-up here was produced by an independent agent and then attacked by a separate
adversarial referee, whose report sits next to it (`referee_*.md`). The paths cited inside the
write-ups point to the session scratchpad; the corresponding scripts are copied into the
`*_checks/` and `families/` directories here. **None of these results proves the uniform theorem
or improves the general growth rate `log R/loglog R`.** Statements below are the referees'
corrected versions where they differ from the original.

## Conditional theorems

* `conditional.md` (referee: survives, minor issues). **Vojta's Main Conjecture with `D = 0`,
  `A = H`, one fixed `eps < 1/2`, for a smooth model of the blow-up of
  `X_7 = {x_1^2+y_1^2 = ... = x_7^2+y_7^2} in P^13` along its diagonal line, implies the uniform
  theorem**, with an absolute `M_0`: `M(C) <= max(M_0, c_0 C^22 + 1)`. The same argument with
  `M_alpha = 3 + floor(2/(1-alpha))` points gives bounded counts on arcs of length `C R^alpha` for
  every `alpha < 1`. No uniformity in a finite place set `S` is needed: only the archimedean place
  is used, and the arithmetic of `N` enters through the height of the point of `X_M`.
  Unconditional by-products: the uniform theorem holds iff for each `C` some set of short-arc
  `M`-tuples is not Zariski dense in `X_M` (Bezout on the product of circles); Ru-Vojta with
  `S = {infinity}` reproves the qualitative Cilleruelo-Cordoba theorem below exponent `1/2`.
  Whether the Ru-Vojta constant `beta_M` exceeds 2 for some `M` (which would give the theorem
  unconditionally) is open. Deriving approximation exponents from Vojta on blow-ups is standard
  (Silverman, McKinnon, McKinnon-Roth); the new part is the variety and the uniformity mechanism.
* `ptolemy.md` (referee: survives; one major over-scoping corrected). A second instance:
  Vojta (`D = 0`) on `Mbar_{0,8}` implies the uniform theorem using only Ptolemy/cross-ratio data.
  It isolates the purely arithmetic missing lemma **(L_8)**: for anchored primitive integer vectors
  `P_0 = (1,0), P_1, ..., P_7` with `det(P_i,P_j) = t_ij prod_{T contains i,j} n_T` (balanced nested
  blocks as in item 532), `H_Sigma >= c w - O(1)` outside a fixed proper Zariski-closed subset of
  `Mbar_{0,8}`. (L_8) implies the uniform theorem. (L_5) is false (infinitely many balanced
  bounded-residue genus-0 curves); (L_6), (L_7) are open. The Vojta instance is vacuous for
  `k <= 7` because `K.beta_bal = 2^(k-3)(k-8)+k+2`.

## Restricted unconditional theorems

* `walsh-extraction.md` (referee: survives with novelty corrections). A residue-pigeonhole
  argument for rational characters (Theorem F) gives: **pure Hadamard-core clusters of every
  type (Sylvester, Paley, ...) are uniformly bounded**, `M <= 24` for `C <= 1` and
  `M < max(28, C^(2+o(1)))` in general (Theorem B, effective, elementary; item 276 had all types
  only for `C <= sqrt 2` via Matveev, item 278 Walsh only). Up to `M/(40 log M)` deleted or modified
  rows are allowed for `M >= 4096` and `M >= C^O(1)` (Theorem B''; item 277 allowed about `log_2 M`).
  Repeated-Walsh profiles with flips on an **arbitrary** row subset satisfy
  `M = o(log R/loglog R)` at the item-335 rate (Theorem C; the growth part also follows from
  items 278 and 333). The sharpened barrier: fully flipped capacity-two prime-order Paley profiles
  pass every pair, integral and rational character, aggregate-collision and capacity inequality
  at `W = Theta(M log M)`. The missing input is the "fully flipped Hadamard phase lemma"
  (Section 10): `M-1` simultaneous Diophantine approximations whose targets move with `M` free
  prime angles. The referee notes that the claim that this is *exactly* the barrier is not proved.
* `cotangent.md` (referee: survives; one major correction). Theorem A (new): with `L_0` the
  intrinsic cotangent scale and `c_p = (p-1)p^(v_p(L_0))`,
  `(M-4c_p)^3 log p/(12 c_p) <= (M/4) log N + binom(M,2) log(C^2/2)` for each split `p`; so a single
  split prime with bounded `v_p(L_0)` forces `M = O(sqrt(log R))`, and `L_0` free of split primes
  below `M/8` forces `M <= (sqrt 24+o(1)) sqrt(log R/loglog R)`. It does not improve the general
  rate, since one colliding pair raises `v_p(L_0)` at negligible cost. Theorem B: at the
  Cilleruelo-Cordoba scale residues, `L_0` and cut deficit are bounded; any five lattice points have
  diameter `>= 1.9433 R^(2/5)`; Cilleruelo-Granville's five-point question is equivalent to a
  bounded-`L` integer clique problem. **Correction:** Theorem C's product identity omits a rational
  factor; the corrected form is `prod z_i = eps N_bal^(M/2) m Phi`, `m = prod p^min(h, M-h)`.
* `auxiliary.md` (no new theorem). Interpolation determinants on one circle factor as
  Vandermonde times a symmetric factor carrying no extra arithmetic, so Bombieri-Pila,
  Heath-Brown and Salberger arguments reduce to the pair data; gap principles, larger-sieve
  variants and Stepanov give nothing beyond `O(log R/loglog R)`.

## Constructions and data

* `families.md` (referee: survives; one major correction to a data claim). Exhaustive scan for
  `R^2 <= 10^12`: at most 3 points on arcs of length `(1/2) sqrt R` and 4 on arcs of length
  `sqrt R`; best constants `C` for `k = 3..8` points are `0.2524, 0.5060, 1.9962, 4.0001, 5.2699,
  5.8978`. A golden-unit eight-point template gives `M_inf(C) >= 8` for `C > 5.89771`, about
  `7.6 * 10^6` times smaller than the notes' Pell constant; an axis Pell family gives `M(C) >= 6` for
  `C > 4`. The referee found a primitive **10-point** cluster with `C = 8.2394` on
  `n = 1176852625 = 5^3 13^2 17 29 113`, so `M(C) >= 10` for `C >= 8.2394`. Heuristics predict
  `M(1/2) = M(1) = 5`; a uniform bound is plausible.

## The sharpest open unconditional target

If for some `M` there is a polynomial on `X_M` whose order of vanishing along the diagonal
exceeds twice its degree, the uniform theorem follows unconditionally (Liouville at the
archimedean place, then the Bezout grid lemma). The Vandermonde/Ramana determinant reaches ratio
`(2s+1)/(s+1) -> 2`. Exact computation finds nothing better for `M <= 7`. This threshold is being
investigated separately.
