# Auxiliary-polynomial and Liouville barriers (2026-10-08/09)

This round asked whether a single auxiliary polynomial could prove the uniform theorem
unconditionally. **It cannot: every route tested here is closed by a refereed barrier theorem.**
Each write-up has an adversarial referee report next to it; statements below are the corrected
versions. `tau_problem.md` is the original problem statement; its "shell decoupling" claim was
wrong (see `verify.md`), but the conclusions below use the corrected full-slice invariant.

## The reduction (verified, `verify.md`)

If a form `F` of degree `D` with `Z[i]` coefficients, nonzero on
`X_M = {x_1^2+y_1^2 = ... = x_M^2+y_M^2}`, vanished along the diagonal line `Y` to order `m > 2D`,
then every arc of length `C sqrt(R)` would hold at most `2D` lattice points once
`R > (K C^m)^2`. The proof is a trigonometric Taylor (Liouville) bound plus a grid lemma.
The relevant invariant is the pseudo-effective threshold of `pi^*H - tE` on `Bl_Y X_M`,
`tau*(M) = sup_{D,s} reg(A_{s,D})/D`, computed on full slices `A_{s,D}` (not single `l^1` shells).

## The barrier (proved, `upper.md`; refereed in `referee_upper.md`)

For every `M >= 2` and integers `p, n >= 0`,
`reg Q^M(p,n) <= phi_M(p,n) = (M-1) min_j (p/(M-j) + n/j)`, by a hyperplane-peeling induction.
Hence

```text
tau*(M) = 2 - 1/ceil(M/2) < 2,
```

attained by alternating Vandermonde orbits (the Ramana/Cilleruelo-Cordoba determinant). So:

* no form on any `X_M` vanishes along the diagonal to order `>= 2 deg F`; the single-place
  Liouville method reaches exactly the exponents `alpha < (k-1)/(2k-1)`, `k = ceil(M/2)`, which
  is the Cilleruelo-Cordoba sequence, and never `alpha = 1/2`;
* the Ru-Vojta constant satisfies `beta_M <= tau*(M) < 2` for every `M`, which answers
  `growth/conditional.md` §6 negatively: Ru-Vojta with `S = {infinity}` cannot give the theorem.

`compute.md`: exact isotypic computation of 1002 shells (`M = 3..9`) and 616 full slices, two
primes each, found nothing above the Vandermonde ratio; computer-assisted values
`tau_shell(3) = tau_shell(4) = 3/2`, `tau_shell(5) = 5/3`.

## Forced divisibility does not help (`construct.md`; refereed)

* Mean-width bound `reg(S) <= 2^{1-n} sum_P w_P(S)` for every finite `S` in `Z^n` (equality for
  boxes). The referee notes it is the Cayley transform of the multiplicity Schwartz-Zippel
  inequality already in `hermitian_polynomial_barrier.md`; the proof here is new and elementary.
* Counting the Gaussian divisibility that the cluster forces on `F(z)` (the tropical bound at
  each split prime) makes the criterion degree-free, and the uniform cube profile (item 27's block
  model) satisfies every such inequality. Products of pair differences, including the
  Vandermonde, become exactly critical, never supercritical. Coefficients depending on `N` gain
  nothing proportional to `log N`.

## The Ptolemy lemma (`L8.md`; refereed)

* **(L_6) is false** unconditionally. A rational 3-parameter family of genus-one curves in
  `Mbar_{0,6}` (the disjoint mixed pencils of `cm_disjoint_mixed_pencil_search.md`), with generic
  positive rank, sweeps out `Mbar_{0,6}` and yields balanced anchored configurations with bounded
  residues outside every proper closed set.
* **No auxiliary-polynomial argument proves (L_k)** for any `k`, nor any `H >= f(w)` with
  `f -> infinity`: the balanced class lies in the interior of the movable cone of `Mbar_{0,k}`
  (Keel-McKernan).
* (L_7) is open; heuristics suggest it is false. (L_8) is open and is a Vojta-type statement.

## Where this leaves the uniform theorem

Local inputs (pair chords, gcds, small-prime collisions, cut bound) are capped at
`log R/loglog R`. Multi-row character certificates fail on fully flipped quasi-random Hadamard
cores. Single-place auxiliary polynomials are capped at `tau* < 2`. Adding valuation-forced
divisibility is exactly critical. The Ru-Vojta and Subspace routes fail, with `S = {infinity}`
(`beta < 2`) and with `S` = primes of `N` (exceptional sets depend on `S`). What remains is a
genuinely Vojta-type Diophantine input for a fixed variety (`growth/conditional.md`, (L_8)).
