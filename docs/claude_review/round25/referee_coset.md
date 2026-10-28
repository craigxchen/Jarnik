# Referee report on round25/coset.md (Proposition A' on realisation cosets)

Referee scripts and outputs: `referee_coset_checks/` (written from scratch for this report; the author's
`coset_checks/` code was read but only its RNG definition of the 60 random lines and the stored `b`, `u` of the
`eta`-line were reused, because those define *which* lines are claimed). Notes: `referee_coset_checks/notes.md`.
Labels: **proved** / **exact** / **evidence** / **heuristic**, as in the overview.

## 0. Verdict

**The partial result survives.** I re-derived every proof with all quantifiers and found no error that changes any
theorem. The mathematical core holds as stated: Lemma 1.1, Lemma 1.2, Lemma 2.1, Theorem 2.2, Corollary 2.3,
Theorem 3.1, Theorem 4.1, Corollary 4.2, Theorem 5.1, Corollary 5.2, Proposition 6.1, Proposition 7.1, Lemma 7.2 and
Corollary 7.3. The exact claims also reproduce independently: the 61 line covers, the graph-coordinate formula, the
`M_0` thresholds and the transplant density table. No growth and no uniformity is claimed, and none follows.

What needs correcting is in the statements around the theorems: the summaries, the readings and (Q_seed). In
order of weight:

| # | severity | item |
|---|---|---|
| 1 | major | (Q_seed), §0 item 6 and the reading of Cor 7.3 quote the heavy-window mass as `e^{-W/4+O(r)}` and the heavy-character condition as `w(v) >= W/2 - O(log W)`. Both drop `log(C kappa_inv)`, where `kappa_inv` is the norm of a right inverse of `(lambda_x)`. That is a **denominator cost**, and Hadamard's bound only gives `log kappa_inv = O(r log(r e_max))`, the same order as `W`. This breaks the owner's constraint in the headline open question (§7 below). |
| 2 | minor-major | "The `A_0` window question reduces to (Q_seed)" (summary). Only one direction is proved: a yes to (Q_seed) gives P. A no gives nothing (its "scenario R" consequence rests on heuristic A), and slack points that are not near a Gaussian-seeded `(pi/4)`-base coset are not covered. coset.md §7.3 says "for the seeding mechanism", but the summary leaves this out. |
| 3 | minor | Remark 2.3(2) says fakes' open point 3 is false "for the coset `T^r`". But `T^r` is a realisation coset of no profile with `M >= 2` (`L != 0`). What Corollary 2.3 refutes is the **full-torus** Poisson extrapolation of fakes §5.1. Open point 3 as posed (realisation cosets) is untouched, as coset.md itself concedes one sentence later. |
| 4 | minor | Theorem 2.2 settles the main question of fakes' open point 4, but not its parenthetical single-coordinate version (fibre over a **genuine** background, which is a tight point, not a slack point). "Supersedes Proposition A'" holds only qualitatively: Corollary 2.3's bound `exp(-C r^2 log P)` is far weaker than A' wherever A' applies. |
| 5 | minor | §6 says "by Theorem 3.1 such a [slack] point is also necessary". Theorem 3.1 does not prove that, because its second implication goes to existence **without** slack. The claim is nonetheless true: a Borel-Cantelli argument (Proposition R1 below) makes Theorem 3.1 an **equivalence**, with characters unslacked. This strengthens the paper. |
| 6 | minor | Numerical error in §8.1 and §8.3: for the primes `13..73` with predecessor transplants, `tau^2 = min p/q = 41/37` (`tau = 1.0527`), not `73/61` (`tau = 1.094`). Corollary 2.3 then gives `n_0 = 1859` and `log mu >= -3.2e4`, not `939` and `-1.6e4`. This is evidence-level only, and nothing depends on it. |
| 7 | minor | Proposition 6.1 and "seeding with arithmetic points is circular" hold for transplants with `N(beta_j) <= p_j` (the convention of `referee_fakes.md` item 7). §1.3 defines arithmetic vectors for **all** `beta`, and arguments of Gaussian integers of unbounded norm are dense, so "non-arithmetic slack point" must be read as "not a transplant with `N(beta_j) <= p_j`". |
| 8 | minor | The summary says "61 lines ..., one of which passes every character test". In fact 10 of the 61 pass the **truncated** character tests: 9 random lines at `T = 22`, plus the `eta`-line. No line's characters are checked beyond a finite weight. I extended the `eta`-line to `T = 40` (it still passes). |
| 9 | minor | The barrier (Theorem 5.1 / Corollary 5.2) is a statement about certificates of exactly the form of Theorem 4.1: full-circle Fubini steps with a per-step union bound over fibres, in lattice coordinates (also finite-index ones) and any order. It does not cover exact computation of `mu_1(A_k(theta_{<k}))`, non-uniform (not inf-over-`theta_{<k}`) integration, or localized steps of Lemma 2.1 type. The summary's "no step-by-step union-bound certificate" is accurate only with this scope. |

The task asked me to check any claimed reduction to route 5 or to one-point caps. coset.md claims none (it neither
uses nor mentions route 5, Prop 4.3 or Theorem 7.1 of round 24), so there is nothing to check there.

## 1. Lemma 1.1, Lemma 1.2, Lemma 2.1, Theorem 2.2 (full-torus positivity)

**Lemma 1.1** (coordinates on `H_0`). Correct. `K` is saturated, so `AB = I_m` with `A` integral. Hence `theta -> B theta`
is injective on `T^m`, its image is a compact connected `m`-dimensional subgroup of `H`, so it equals `H_0`, and
`B^T v = 0` iff `v in L_R ∩ Z^r = L^sat`.

**Lemma 1.2** (Gaussian primes of distinct split norms). Correct. `X^2 + Y^2 = prod q_j^{|v_j|}` is odd, so `X ± Y` is
odd. `Y = 0` or `X = 0` gives `conj(gamma_v) = ±gamma_v`, and unique factorisation in `Z[i]` excludes this, because
`conj(beta_j)` is not associate to `beta_j` for split odd `q_j`, and only `beta_j` (or only `conj(beta_j)`) occurs for index
`j`. The two distance formulas, `arcsin(min(|X|,|Y|)/|gamma|)` and `arcsin(min|X±Y|/(sqrt2 |gamma|))`, are right.
*Exact check* inside `seed_check.py` (Gaussian integer arithmetic, see §7).

**Lemma 2.1** (slack points are density points). Re-derived step by step.
* (i) `y'_v = y_v prod (p_j/p'_j)^{|v_j|/2} >= tau^{|v|_1} y_v`, together with `arcsin b - arcsin a >= b - a`, gives an excess of at least
  `(tau^{|v|_1} - 1) y_v` in the `0`-family and at least `sigma_v` in the `pi/4`-family. Correct.
* (ii) Both distance functions are 1-Lipschitz on `T` and `|<mu_v, t>| <= |mu_v|_1 eps`. For `v in L^sat` the value is
  constant on the coset and the condition holds at `psi` for `p' <= p`. Correct.
* (iii) The coordinate with `|mu_{v,k}| = |mu_v|_inf` gives density `<= 1/(2 eps |mu_v|_inf)`, and convolution keeps the sup
  bound. An interval of length `2|mu_v|_1 eps` meets at most `2|mu_v|_1 eps/(pi/2) + 2` translates of an interval shorter than
  `pi/2`, for each of the two families. The window lengths are `<= pi y` and `<= pi y/sqrt2`. Inserting
  `y_v <= sqrt2 |mu_v|_1 eps/(tau^n - 1)` (the definition of `U(eps)`) gives
  `pi(1+sqrt2)/2 (|mu|_1/|mu|_inf)(4|mu|_1 eps/pi + 2)/(tau^n - 1)`. With `|mu|_1/|mu|_inf <= m` and `2|mu|_1 eps/pi < b*|v|_1/2`
  (because `eps < pi/4` and `|mu_v|_1 <= b*|v|_1`), this is `<= c_2 m (1 + b*|v|_1)/(tau^n - 1)` with `c_2 = pi(1+sqrt2)`. Correct.
* (iv) The series converges, since `#{|v|_1 = n}` is polynomial in `n` and `tau > 1`. Each non-character `v` leaves `U(eps)`
  as `eps -> 0`, and dominated convergence on counting measure gives `F(eps) -> 0`. Correct.

**Theorem 2.2** (`mu(G(p)) > 0` iff `5 ∉ p`). Correct. Proposition B (refereed) gives the `5 in p` direction. The predecessor
split primes `q_j` are distinct with `5 <= q_j < p_j`, so Lemma 1.2 gives `arg beta in G(q)` and Lemma 2.1 (`B = I`, `b* = 1`)
gives positivity. There is indeed no condition on `Lambda`. *Evidence* (`density_ref.py`, my own enumeration of all
30,228 monomials of weight `<= 24`, 2000 samples per scale). Survival in cubes around the transplant vector for
`13..73` is `0.001, 0.0915, 0.4845, 0.9195, 1.0` at normalised half-sides `1e-3, ..., 1e-7`, and uniform points give `0`.
This matches §8.1. Note that §8.1's `eps` is in normalised units `x = phi/(pi/2)`, while Lemma 2.1's is in angle units.
For this prime set Proposition A' already fails at `k = 2` (`kappa_1 Pi_{13} = 3.89 > sqrt17 - 1 = 3.12`), so Theorem 2.2
is a genuine extension.

*Scope corrections (finding 4).* fakes' open point 4 has two parts. Theorem 2.2 answers the first ("does some set of
split primes `>= 13` have `mu(G(p)) = 0`? No."). It does **not** answer the parenthetical: is the fibre of admissible
`phi_p` over a **genuine** background `S` null? A genuine background is tight (`min(|X|,|Y|) = 1` occurs), so Lemma 2.1
does not apply. "Supersedes Proposition A'" should read "supersedes A' **as a positivity statement**": quantitatively
A' gives `prod(1 - U_k)`, which is far larger than Corollary 2.3 wherever A' applies.

## 2. Corollary 2.3 and Remark 2.3(2) (Poisson extrapolation)

**Corollary 2.3.** Correct. Rankin's trick with `z^2 = 1/tau` gives, for `|v|_1 >= n_0`,
`tau^{-|v|_1} <= tau^{-n_0/2} z^{|v|_1}`. The identity `sum_v (1+|v|_1) z^{|v|_1} = f + z f' = f(1 + 2rz/(1 - z^2))` holds,
and the sum over `v/±` absorbs the `1/2`. Safety for `|v|_1 <= n_0 - 1` follows from `y_v >= P^{-|v|_1/2}` and
`tau^{|v|_1} - 1 >= tau - 1`. The asymptotic constant `log((1+z)/(1-z))/log tau = 2.5431` at `tau = 2` is reproduced
(`constants.py`).

*Numbers* (`constants.py` -> `constants_out.json`): `n_0 = 64, 547, 12772` and `log mu >= -1776, -2.52e5, -2.21e8` for
`(r, P) = (8, 1e3), (100, 1e4), (2500, 1e6)` at `tau = 2`, as claimed. The Poisson value at `r = 2500`, `P = 1e6` is `-162.15`,
as claimed. **Finding 6:** the `r = 8`, `P = 73` row uses `tau^2 = 73/61`, but the correct value is
`tau^2 = min_j p_j/q_j = 41/37` (ratios `2.6, 1.31, 1.71, 1.28, 1.108, 1.29, 1.15, 1.20`). With it, `tau = 1.0527`,
`n_0 = 1859` and `log mu >= -3.2e4`. §8.1's "small slack `tau = 1.094`" is therefore `1.053`. Nothing depends on this number.

**Remark 2.3(2).** The comparison is right: `Pi >= exp(M^{delta/2}/4)` against a proved
`log mu >= -(5.1 + o(1)) M^2 log M`, with `tau = 2` seeds `q_j <= P/8`, which exist because `r <= M` is far below the
`~P/(16 log P)` split primes `<= P/8`. So the **full-torus** extrapolation of fakes §5.1 is false for large `M`. The
crossover at `Lambda >~ log r` is also right: `1.1 e^{2 Lambda} ~ 2.54 r^2 log P` gives
`Lambda ~ log r + (1/2) loglog P`.

**Finding 3.** fakes' open point 3 concerns `nu(G ∩ coset)` on **realisation** cosets. `T^r` is the realisation coset of
a profile only when `L = 0`, that is `M = 1`. The sentence "fails for the coset `T^r` itself" therefore refutes a
reading that open point 3 never intended. The correct statement is: *the full-torus Poisson extrapolation of fakes
§5.1 is false for large `Lambda`; open point 3 on realisation cosets is open.* coset.md's next sentence says the
second half, so this is wording.

## 3. Theorem 3.1 (sandwich), and its sharpening

**Theorem 3.1.** Correct and short. The first implication is Lemma 2.1 at a point of `G(p') ∩ C`: characters are constant
on `C` and hold for `p' <= p`. The second is trivial.

**Finding 5, and a strengthening.** §6 says a slack point is "also necessary" by Theorem 3.1. Theorem 3.1 does not show
this, since its right-hand side is existence in `G(p)`, not in `G(p')`. The claim is true all the same:

> **Proposition R1 (referee; proved).** Let `C = g + H_0` be any coset with Haar measure `nu`, and let all `p_j >= 5`. Then
> `nu(G(p) ∩ C) > 0` **iff** there are `phi in C` and `tau > 1` such that `phi` satisfies every character condition
> (`v in L^sat`) for the weights `p`, and every non-character condition for the weights `p/tau^2`. Moreover
> `nu`-almost every point of `G(p) ∩ C` is such a point.

*Proof.* (⇐) This is Lemma 2.1 verbatim. Its proof uses slack only for non-characters (step (ii) needs characters only
at `psi`, for the weights `p`).

(⇒) Fix `1 < tau_1 < sqrt5`, put `y^{(1)}_v = min(1, tau_1^{|v|_1} y_v)`, and let `W^{(1)}_v` be the enlarged windows. For a
non-character `v`, `theta -> <v,g> + <mu_v, theta>` pushes Haar measure on `T^m` to Haar measure on `T`, because
`mu_v != 0`. So `nu{phi in C : <v,phi> in W^{(1)}_v} <= (2+sqrt2) tau_1^{|v|_1} y_v`. Also
`sum_v tau_1^{|v|_1} y_v = prod_j (1 + 2 tau_1 p_j^{-1/2}/(1 - tau_1 p_j^{-1/2})) < infinity`. By Borel-Cantelli, `nu`-a.e.
`phi in C` lies in only finitely many `W^{(1)}_v`. The boundaries `∂W_v` are finite sets with null preimages, and there
are countably many `v`, so a.e. `phi` lies on none of them. Take `phi in G(p) ∩ C` outside these null sets. For the
finitely many failing non-characters `v`, `phi in G(p)` and `phi ∉ ∂W_v` give strict inequality in both families, so
some `tau(phi) in (1, tau_1]` keeps them outside the `tau(phi)`-enlarged windows. All other non-characters are outside
the `tau_1`-enlarged windows, hence outside the `tau(phi)`-enlarged ones. QED

So the "sandwich" is an equivalence: **Haar positivity on a coset = existence of a point with valid characters and
uniformly slack non-characters.** This is the precise form of coset.md's "Haar positivity is existence with slack",
and it justifies the "necessary" in §6.

## 4. Theorem 4.1 and Corollary 4.2

**Theorem 4.1.** Correct. Conditions with `mu_v = -mu` are the conditions of `-v` (`W_v` is symmetric), so `B_mu` captures
all of them. `t -> mu_k t + c` preserves Haar measure on `T`, so each bad set has measure exactly `U(mu)` uniformly in
`theta_{<k}`, and Fubini closes the argument.

**Corollary 4.2.** Re-derived. `mu_v = v_J + C^T x`. `last(mu_v) = k` pins `v_{j_k'} = -f(x)_{k'}` for `k' > k` and excludes one
value of `v_{j_k}`. The sums give `Pi_{J,<k}` and `Pi_{j_k} - p_{j_k}^{-|f(x)_k|/2}`, and the `1/2` comes from `v -> -v`.
*Independent exact check* (`graph_check.py`): free primes `17, 29, 37`, dependent `13, 41`, `C = [[1,-2,1],[0,3,-1]]`.
Brute-force sums over a box agree with the formula to relative error `2.9e-14`, `6.0e-13` and `1.4e-7` (box truncation).
The truncated true fibre-union masses `U_k = 0.725, 1.378, 5.129` lie below `(kappa_1/2)` times the formula
(`0.737, 1.439, 5.59`), as Corollary 4.2 says.

The "Reading" (`Lambda <= (1/4) log P - (1/2) log(1+h) - O(1)`) is labelled as a sketch and is used nowhere. I checked
that its order of magnitude is right: dependent primes whose last entry is at step `k` contribute `~ p_i^{-1/2}`
without a `p_{j_k}^{-1/2}` factor. It stays a sketch.

## 5. Theorem 5.1 and Corollary 5.2 (barrier)

**Theorem 5.1.** Correct. For `v in Z^J \ 0`, `v ∉ L_R`, so `mu_v != 0`. `mu_v = ±mu_{v'}` forces `v ∓ v' in L_R ∩ span(e_J) = 0`,
so the map `(Z^J \ 0)/± -> ∪_k F_k` is injective. Then `U(mu_v) >= omega(v) >= kappa y_v`, which sums to `(kappa/2)(Pi_J - 1)`.
The argument works in any basis of `K`, and also for a finite-index sublattice basis (an isogeny chart), because
`B^T v = 0` iff `v in L_R` still holds. *Check* (`graph_check.py`): truncated `sum_k U_k = 7.23 >= (kappa/2)(Pi_J - 1) = 2.53`.

The remark that A' is optimal for its method up to lower-order terms is also right. With `Pi ~ e^{2 Lambda}`, the
barrier is `Lambda <~ (1/2) log r = (1/4) log P + (1/2) log Lambda`.

**Corollary 5.2.** Correct. `|J| = m` comes from Theorem A(e), `Pi_J >= exp(c M^{delta/2})`, and `m <= r < M^2`. The function
`f(M) = c M^{delta/2} - log(1 + 2M^2/kappa)` decreases and then increases (it can be positive at `M = 1` when `c > 0.65`),
so "the last zero" is the correct threshold. *Independent recomputation* (`constants.py`, 40-digit bisection on
`log M`): `M_0 = 4.48e3, 4.55e5, 6.99e8, 4.93e17, 4.33e20` for `(c, delta) = (1/4,1), (1,1/2), (1/4,1/2), (1/2,1/4), (1/4,1/4)`.
These agree with the claimed `4.5e3, 4.6e5, 7.0e8, 4.9e17, 4.3e20` to the stated precision. `M_0` depends on `(c, delta)`
only, so the uniformity constraint holds.

**Scope (finding 9).** The theorem is about certificates of the form "`U'_k >= sum_{mu in F_k} U(mu)`, each `< 1`". Lemma 2.1
itself, which localizes to cubes, escapes it. So would an exact or non-uniform evaluation of `mu_1(A_k(theta_{<k}))`. The
statement "the judge's proposed route cannot reach the window" is accurate for the route as the judge sketched it
(`round21/README.md` §7). It is asymptotic: for `(c, delta) = (1/4, 1/4)` it bites only at `M >= 4.3e20`.

## 6. Proposition 6.1 (arithmetic seeds are actual clusters)

Correct as stated (with `N(beta_j) <= p_j`).
* `arg z_x = <u_x, arg beta>` modulo `2 pi`, and realisation modulo `pi/2` is absorbed by units.
* `z_x/z_y = gamma_{2d}/prod N(beta_j)^{|d_j|}` (checked: `beta/conj(beta) = beta^2/N(beta)`). If `eps_x z_x = eps_y z_y`, then
  `<lambda, psi> in (pi/2)Z`, which `psi in G(p)` forbids.
* The arc constant is `C_beta = C (N_beta e^{-W})^{1/4} <= C`, and Theorem 1.1 of the expert note
  (`3M log M <= log R^2 + (14 + 4 log^+ C) M`, with `R^2 = N_beta`) gives the bound. The common phase `theta_0 = arg z_1 - delta_1`
  is kept exactly.

**Finding 7.** The circularity applies to transplants with `N(beta_j) <= p_j`, which is the definition fixed in
`referee_fakes.md` (item 7). An arithmetic point `arg beta` with large norms is an actual cluster with
`C_beta = C (N_beta/e^W)^{1/4} > C`, which is no contradiction. Such points are dense in `T^r`, so they are not a
meaningful class. §1.3 should restrict "arithmetic" to `N(beta_j) <= p_j`.

The remark that "characters forbid `H_0`" is correct: on `H_0` every `L^sat`-character vanishes identically.

## 7. Proposition 7.1, Lemma 7.2, Corollary 7.3, and (Q_seed)

**Proposition 7.1.** Correct. For a non-heavy non-character, the shift `<v, eta>` (a distance in `T`) uses at most
`(tau^n - tau_1^n) y_v/sqrt2` of the excess, which leaves `sigma_v(tau_1)`. Heavy non-characters have `mu_v != 0`, so their
bad sets have Haar measure exactly `omega(v)` on `T^m`.

**Lemma 7.2.** Correct. `L ⊂ 2Z^r` because `u_x = 2a_x - e`, so `(pi/4)`-points lie in `H`. An odd multiple of `pi/4` swaps the
two families, so both distances are `>= arcsin(1/(sqrt2 |gamma_{mu_v}|))`. Also
`log |gamma_{mu_v}|^2 <= w(v) - |v|_1 log(2 tau^2)` and `2^{|v|_1/2}/sqrt2 >= 1`. The converse is right.
*Exact test* (`seed_check.py`): `r = 4`, `m = 2`, `B = [[1,0],[0,1],[1,1],[1,-1]]`, `q = (5, 13)`, `p = (17, 41, 197, 229)`,
`tau^2 = 1.5`. I enumerated all 256 shifts `g_0 in (pi/4){0..3}^4` and all 588 non-characters with `w <= 22` (one per sign),
computing `gamma_mu` in exact Gaussian integer arithmetic. The minimum relative margin over all weight-`p/tau^2`
conditions is `0.026 > 0`.

**Corollary 7.3.** Correct as stated, with `rho` explicit. The majorant `f(z)(1 + 2b* r z/(1 - z^2))`, safety from
`b*|v|_1 eps_* < sigma_v(tau_1)`, and `tau^n - tau_1^n >= tau - tau_1` all check.

**Finding 1 (major): the denominator cost disappears in the reading and in (Q_seed).** §7 opens correctly: `rho` is
`Delta` times the norm of a right inverse of `(lambda_x)`, "which must not be treated as uniform". The reading and
(Q_seed) then drop it:
* `log H_s(rho) = s log(sqrt2 C kappa_inv/(tau - tau_1)) - sW/4 + O_s(sum_j p_j^{-(1-s)/2})`. That is `-sW/4 + O_s(r)` **only if**
  `log kappa_inv = O(r)`. The right inverse has rational entries with denominator a maximal minor of
  `(2(a_x - a_1))`, and Hadamard only gives `log kappa_inv <= O(r log(r e_max))`. For example, with `e_j = 2` and `a_1 = (1,...,1)`, the
  differences `a_x - a_1` can form the unit upper-triangular matrix with `-1` above the diagonal. Its inverse has entries
  `2^{k}`, so the naive right inverse has norm `2^{Theta(r)}`. (I did not check whether such designs meet the other row-3
  inputs, or whether the `mod pi/2` freedom in choosing `eta` reduces the norm. The point is that nothing in the paper
  bounds it.) At barrier scale `W = Theta(r log M)`, so this term can be of the order of `W/4` itself.
* Lemma 7.2's converse gives `w(v) >= W/2 - 2 log(sqrt2 |v|_1 C kappa_inv)`. That is "`W/2 - O(log W)`" only under the same
  hypothesis. As written, (Q_seed) states a **stronger** necessary condition than the lemma proves.

Corrected (Q_seed): in (1), "every `v in L^sat \ 0` has `arcsin(y_v/sqrt2) <= |v|_1 rho`, `rho = C kappa_inv e^{-W/4}`"; in the
analytic part, "the heavy-window mass `H_s(rho)`, which is `e^{-sW/4 + s log(C kappa_inv) + o(r)}`". The reading
`m <= M^{(2-delta)/5}` then additionally needs `log kappa_inv + log b* = o(W)` as a profile hypothesis. It is a sketch
in any case.

**Finding 2.** (Q_seed) gives a sufficient route to P. It is not a reduction of the `A_0` window question. By
Proposition R1, P on a realisation coset is **equivalent** to a slack point with valid characters on that coset. Gaussian-
seeded `(pi/4)`-base cosets are one source of such points. Others (higher torsion, open point 4; slack points that are
not `Delta`-close to any seed) are not excluded. A no to (Q_seed) implies scenario R only through heuristic A.

## 8. Computations (independent re-runs)

| claim | referee script | result |
|---|---|---|
| 60 random lines at `T = 22` covered by non-character windows | `cover_check.py` (own enumeration, own interval construction and sweep; 17,017 monomials up to sign) | **reproduced, exact**: all 60 certified, minimum overlap `6.016e-4` |
| `eta = 1e-3` line covered at `T = 20` | same | **reproduced, exact**: overlap `4.392e-6`, also at `T = 24`. The `eta = 0` line is not covered (survivors) |
| `eta = 1e-3` line "passes every character" | `char_check.py` (characters only, `v ⊥ b`, up to `T = 40`: 134,310 characters) | holds to `T = 40`, min margin `1.12e-6` (evidence, truncated) |
| density near transplant (§8.1) | `density_ref.py` | reproduced (see §1) |
| Cor 4.2 formula | `graph_check.py` | reproduced (rel err `<= 1.4e-7`, truncation) |
| Cor 5.2 `M_0` | `constants.py` | reproduced |
| Cor 2.3 numbers | `constants.py` | reproduced for the stated `tau`. The `r = 8` row uses the wrong `tau` (finding 6) |
| Lemma 7.2 | `seed_check.py` | exact test passes |

*Rigour of the covers.* `g` is an exact binary rational. Every interval endpoint is computed to within about `1e-14`
(float64 dot products of integer vectors of norm `<= 30`, `asin` to about `1 ulp`, division by integers). The certifying
overlaps are `>= 4.4e-6`, which is `8` orders of magnitude larger, so "exact" is justified. The covers use only
non-characters, so the conclusion `G(p) ∩ line = ∅` does not depend on the truncation of the character tests.

*Significance.* The lines are random 1-dimensional cosets in `T^8` with primes `13..73`. They are not realisation cosets of
any profile, and they are not in the window regime. They are evidence for heuristic A, as labelled, and nothing more.

## 9. Owner constraints

* *Uniform exceptional sets from fixed-`S` Subspace theorems:* none used. Complies.
* *Gaussian gcd / denominator costs:* `L ⊂ 2Z^r`, the torsion of `L^sat/L`, Lemma 1.2's unique-factorisation argument and
  Proposition 6.1's `C_beta` are all exact. **Not complied with** in the reading of Corollary 7.3, in §0 item 6 and in
  (Q_seed), where `log kappa_inv` (a determinant denominator) is absorbed into `O(r)` or `O(log W)` (finding 1).
* *Drifting common phase replaced by zero:* not done. `theta_0` is absorbed exactly into `Real(delta)`, the translate `eta`
  is kept (`Heavy(eta)`), and Proposition 6.1 uses `theta_0 = arg z_1 - delta_1`. Complies.
* *Template-dependent threshold treated as uniform:* `M_0(c, delta)` is uniform. The profile constants of Corollary 7.3
  (`b*`, `kappa_inv`, `tau`) are kept as hypotheses in the proved statements, but the asymptotic reading
  `m <= M^{(2-delta)/5}` silently assumes `log kappa_inv`, `log b* = O(r)`. It is labelled a sketch, so this is minor.

## 10. Corrected statement

**Partial; no growth, no uniformity.** Proved: Lemma 2.1 (slack ⟹ Lebesgue density point, on any coset). Theorem 2.2
(`mu(G(p)) > 0` iff `5 ∉ p`, for every `Lambda`), which settles the main question of fakes open point 4 but not its
genuine-background fibre version. Corollary 2.3 (`log mu(G(p)) >= -(2.54+o(1)) r^2 log P` at `tau = 2`, explicit), which
refutes the **full-torus** Poisson extrapolation of fakes §5.1 for large `Lambda` (crossover `Lambda ~ log r`), while
fakes open point 3 on realisation cosets stays open. Theorem 3.1, sharpened by the referee's Proposition R1 to an
**equivalence**: Haar positivity on a coset holds iff the coset contains a point with valid characters and uniformly
slack non-characters, and a.e. point of `G(p) ∩ C` is such a point. Theorem 4.1 / Corollary 4.2 (Proposition A' on
cosets, exact step sums). Theorem 5.1 / Corollary 5.2: certificates of the form of Theorem 4.1 (full-circle Fubini
steps with per-step fibre union bounds, any lattice basis, any order) cannot certify positivity when `m >= cM` and
`M >= M_0(c, delta)` (`M_0 = 4.6e5` at `(1, 1/2)`). Proposition 6.1 (transplant points with `N(beta_j) <= p_j` on `Real(delta)` are
actual clusters with `log N_beta <= W`, `C_beta <= C`). Proposition 7.1, Lemma 7.2 and Corollary 7.3 (positivity on
`Delta`-translates of a Gaussian-seeded `(pi/4)`-base coset when `(1/2)(4 eps_*/pi)^m > H_s(rho)`, with `rho = C kappa_inv e^{-W/4}`
explicit; such seeds force `arcsin(y_v/sqrt2) <= |v|_1 rho` for every character).

Sketch only: the range `m <~ M^{(2-delta)/5}`, which additionally needs `log kappa_inv + log b* = o(W)`.

Exact: `G(p) ∩ line = ∅` for 61 explicit lines in `T^8` (primes `13..73`); 10 of them pass the truncated character
tests.

Open: (Q_seed) with the corrected heavy condition `arcsin(y_v/sqrt2) <= |v|_1 C kappa_inv e^{-W/4}` and heavy mass `H_s(rho)`.
It is a sufficient route to scenario P, not a reduction of the `A_0` window question. By J2, neither answer gives growth.
