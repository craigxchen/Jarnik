# Referee report: "The integer-cotangent route" (growth/cotangent.md)

Referee run of 2026-10-08. This run continues an earlier referee attempt that was interrupted after
writing its check scripts. All referee scripts are in `growth/referee_cotangent_checks/`, and every
claim marked *checked* below was reproduced there with exact integer or rational arithmetic.

## Verdict

**The result survives as a restricted theorem, with one false statement that must be corrected
and several overstated claims.** No fatal error was found. As the author says, the uniform theorem
is not proved and the general growth rate is not improved.

| Claim | Status |
|---|---|
| Dictionary (D1)–(D8), master identity (M1)–(M4) | Correct. Mostly known material (see §2), re-derived. |
| Theorem A (split content of `L_0` versus level spreading) | Correct. New relative to the notes. The constant can be improved and the hypothesis weakened (§3). Narrow in scope. |
| Theorem B (rigidity at the CC scale; five-point diameter `>= 1.943 R^(2/5)`; CG five-point question reformulated with bounded `L`) | Correct. It is a direct corollary of the classical product identity. The 1.943 constant and the bounded-`L` reformulation are new but modest. |
| Theorem C (direction pinning for even `M`) | **The stated identity is false.** The corrected identity is true and every consequence survives. The mechanism is already in the notes for exactly balanced clusters (item 15). |
| Pell family: `prod z_i` is a rational multiple of `-3+4i` "for every member" | True for every member, by a one-line identity from the notes. The write-up only checked seven members. |
| Data: the five-point diameter ratio "never falls below 4.9475" | False as stated (`N=25` gives 4.2024). It holds only in the span-at-most-90° regime (`A>=L`) that the search actually covered. |

## 1. What was re-run

| Script | Result |
|---|---|
| `check_core.py` (earlier referee run) | On 1323 actual primitive clusters (`M` from 3 to 9, some non-consecutive), every case passes: (D4) chord, residue and parity rule (18518 pairs), (D8) cut identity, (M1) in fourth-power form, level capacity, the Theorem A chain, the layer lemma, and `6 | L_0` (788 clusters with `M>=5`). The *literal* Theorem C fails on 453 even clusters. |
| `deficit_dp.py` (new: exact shortest-path DP over all level vectors) | The cubic deficit lemma `D >= (M-4c)^3/(12c)` has **no violation** over `1<=c<=40`, `4c<M<4c+400` (15960 cases). `min D/bound = 1.0227`. `D_min/(M^3/(12c))` tends to 1 (0.999 at `c=1, M=3000`), so the lemma is asymptotically sharp. |
| `thmA_active.py` | Theorem A's chain `lhsA <= sum D_p log p <= MW/4 + P log(C^2/2)` on 80 actual clusters with `v_5(L_0)=0` and `M` from 17 to 44, so the **cubic term is active**. Capacity 4 is saturated at `p=5`. Note: the author's `check_theoremA.py` has `M<=12<4c_5+1`, so it never exercises the cubic term. |
| `five_point_exhaustive.py 1000000` | **Every** circle with `N<=10^6`: 3,138,560 cyclic five-windows. The diameter bound `diam^20 >= 2^14*36*N^4` holds on all of them. (B2) `128 sigma prod s N^2 <= dmax^10` holds on all 2,315,520 primitive half-plane windows. `6 | L_0` holds on all of them and on 51,331 random non-consecutive five-subsets. |
| `thmC_literal.py 20000` | Even windows (`M=4,6`), `N<=20000`: the literal statement fails on 34,164 of 100,256. The corrected statement holds on all 100,256. |
| `check_M3.py` | Exact chart identity (M3) and inequality (M4) on 1813 actual clusters (`M` from 3 to 8). |
| `check_D5.py` | Good-prime dictionary (D5) on 2593 (cluster, prime) pairs. No failures. |
| `cg_regime_crosscheck.py 1000000 6` | Lists every primitive five-window with `diam/N^(1/5)<6`. See §6. |
| Author's `check_dictionary.py`, `check_theorems.py`, `check_pinning.py`, `check_theoremA.py`, `check_level_deficit.py`, `cyclo_scale.py` | All pass and reproduce the stated numbers (cyclotomic: `log L_0/log N = 0.2388`, largest residue `N^0.2306`). |

## 2. Dictionary and master identity

The proofs were checked line by line.

* (D4) follows from `c^2+1 = 2N(N+<z_i,z_j>)/det^2`. For an antipodal pair, `c=0` (use
  `c = i(z_j+z_i)/(z_j-z_i)`), and the formula still holds.
* (M3) follows from (D8) and (M2) by pure algebra.
* (M4) uses `prod_{i<j}|X_i-X_j| <= prod_j X_j^(j-1)`, and the exponent count gives
  `sum_j (2k-2j+2) = k(k+1) = 2P`. Correct.

**Novelty.**

* (D1)–(D8) are, as the author says, largely in `integer_cotangent_*` notes:
  * `split_prime_cotangent_capacity.md` (item 151);
  * `integer_cotangent_offset_height_cut_residue_dictionary.md`;
  * `intrinsic_cotangent_scale_triangle_divisibility.md`.
* (M1) is in the notes for `M=8`: `balanced_bonus_phase_audit.md` (2) reads
  `prod |z_i-z_j|^2 = 4^28 T^2 N^12 e^D`, with the same cut-count argument.
* The CC exponents are `integer_cotangent_affine_content.md` (15).
* The general-`M` form with `Delta_*`, and the target (SA), are tidy restatements, not new theorems.
  The write-up's labels are broadly fair here.

## 3. Theorem A

**Proof check.**

* Level capacity `c_p = (p-1)p^(v_p(L_0))` is item 151, and its proof is correct: `u_i+u_j` is a
  unit because `p` is odd.
* Step 1 (Plotkin), step 2 (`n_ij >= 2N^(1/2)/C^2`) and step 3 are correct.
* Primes `p` that do not divide `N` contribute 0, because `M <= c_p` for them.
* Deficit lemma: the argument is right. One slip: the condition "`(j+2)c <= M/2`" should read
  "`(j+1)c <= M/2`". The latter is what guarantees that `[M/2+jc, M/2+(j+1)c)` lies below
  `S_e >= M-c`, and it is the condition needed for the stated sum up to `J-1`. The final bound is
  unaffected, and the exact DP confirms it.
* (A2) and (A3) follow (Mertens in the progression 1 mod 4).

**Weaknesses.** These are not errors.

1. *The hypothesis of (A3) is far stronger than needed.* All terms of (A1) are nonnegative, so
   **one** split prime suffices. If `v_p(L_0) <= B` for a single split `p`, then
   `M^2 <= (3c_p/log p) W (1+o(1))`, that is,
   `M <= (sqrt(6 c_p/log p) + o(1)) sqrt(log R)`. For example, `5 ∤ L_0` already gives
   `M <= (sqrt(24/log 5)+o(1)) sqrt(log R) ≈ 3.86 sqrt(log R)`.
2. *The constant in (A2) is loose by a factor `sqrt 8`.* Mertens mass concentrates on
   `p <= M^(1-delta)`, where `(M-4(p-1))^3 = M^3(1-o(1))`. So the same hypothesis gives
   `W >= (1-o(1)) M^2 log M/6` and `M <= (sqrt 24 + o(1)) sqrt(log R/log log R)`, not
   `sqrt 192`.
3. *Scope.* The hypothesis (`L_0` free of small split primes, or of bounded split valuation) is a
   non-generic property. Nothing shows that a hypothetical unbounded short-arc family must have
   it; one colliding pair per level raises `v_p(L_0)`, as the author says. So Theorem A does not
   bear on the uniform theorem or on the general rate.

**Novelty.** The cubic level-spreading deficit and its combination with item 151 do not appear in
the notes. Searches for spread/deficit/cubic/capacity, `inert_prime_residual_growth.md` §5 and
item 151 found nothing like it. The lead author's 2/3 sketch has only the quadratic chain lemma.
**New, but a restricted result.**

## 4. Theorem B

**Proof check.**

* (B1) is (M1) together with `prod eps <= 2^(floor(M^2/4))` and every chord at most
  `K N^(theta_M)`. Correct.
* For `M=5`: `Delta_* = sigma^2`, because singleton layers have exponent `(3/2)^2 - 1/4 = 2`
  and the layers with `S=0` or `S=5` cannot occur. So (B2) is correct.
* `6 | L_0` for five points follows from `integer_cotangent_local_bounds.md`. If `L` is odd, or
  `3 ∤ L`, there are at most 3 finite coordinates. Apply this with `L = L_0`.
* So `prod s >= 6`, `K^10 >= 768` and `K >= 1.9433`.
* Non-primitive five-tuples reduce to primitive ones with a smaller `K`. Antipodal pairs are
  trivial.
* All of this is checked exhaustively for `N <= 10^6`. The minimum half-plane ratio is 4.2024, at
  `N=25`.
* The block-balance bound (B3) was re-derived. The `xi` vector satisfies `|xi| <= 9 log(K^2/2)`,
  and `||(I+Pet)^(-1)||_inf = 2` (Petersen srg identity `A^2 = 2I - A + J`). This gives
  `|log b_T - W/10| <= 58.5 log K`, so 60 suffices.

**The CG reformulation.** It is correct in both directions.

* Forward: fixed `K` bounds `L_0` and the Gaussian gcd, so infinitely many radii give infinitely
  many primitive cliques with one fixed `L_0 ∈ 6Z`.
* Converse: the arc is at most `2 c^(3/10) R^(2/5)`, and infinitely many cliques at fixed `L`
  force `N` to be unbounded.

**Calibration.**

* The plain CC product, with `s>=1`, already gives five-point diameter at least
  `2^(7/10) R^(2/5) = 1.6245 R^(2/5)`, which beats `sqrt 2`. The residue input (`prod s >= 6`)
  only accounts for the step from 1.6245 to 1.943.
* Bounded residues, `L_0` and deficit at the CC scale are immediate corollaries of the classical
  identity.
* I count these as **new but modest**: a constant improvement plus an equivalent reformulation.
  Neither resolves the CG five-point question.
* Literature: the primary Cilleruelo–Granville papers could not be fetched; the egress proxy
  returned 403 for arXiv and Cambridge. The project notes (`endpoint_continuation_status.md`
  around line 13998) record that CG's question concerns five points at `O(R^(2/5))`. That reading
  is consistent with four points being realized at scale `R^(1/3)` (the Pell and Fibonacci
  families).
* Side note for the user: `paper/jarnik_sublog.tex`, line 332, in the repository still says the
  question is open "for four lattice points on an arc of length `o(R^(2/5))`". The notes call
  that a mistaken claim that was corrected, and the Pell family at `R^(1/3)` contradicts it.

## 5. Theorem C

**The stated identity is false.** The write-up's Theorem C and the summary both say
"`prod z_i = eps N_bal^(M/2) Phi`, `Norm(Phi) <= Delta^2`", with `N_bal` the product of `p` over
balanced layers. This drops the rational factor `p^(min(h,M-h))` contributed by every
*unbalanced* layer.

Smallest counterexample: `N = 25`, `M = 4`, points `(-4,-3), (-3,-4), (0,-5), (3,-4)`. Then:

* `prod z_i = 375 - 500i`;
* `N_bal = 5` and `Delta = 5`;
* `prod/N_bal^2 = 15 - 20i`, which has norm `625 > Delta^2 = 25`.

The literal statement fails on 34% of all even windows with `N <= 20000`.

**Corrected statement.** It holds on all 100,256 windows and follows from the author's own layer
lemma:

    prod z_i = eps * N_bal^(M/2) * m * Phi,   m = prod_{unbalanced layers} p^(min(h,M-h)) in Z_{>0},
    Norm(Phi) = prod_{layers} p^(|M-2h|) <= Delta^2.

`m` is a positive rational integer, so `arg prod z_i = arg(eps Phi)`. Every consequence therefore
survives:

* `|M arg z_1 - arg(eps Phi)| <= (M-1) delta`;
* finitely many pinned directions, each with `e^(2iM theta)` in `Q(i)`;
* for `N >= N_0(K)`, no rational direction. The proof via `R eta^2 -> 0` is correct, using
  `theta_M < 1/4`.

**Novelty.** The mechanism, "the product of the points pins the mean phase", is in the notes:

* item 15 / `balanced_pair_products.md` §2: for exactly balanced clusters,
  `prod z_i = eps^M R^M`, with a rational-branch exclusion;
* `codex_uniform_bound_research.md` (6.39m).

The Pell family's limiting direction `arctan((sqrt 5-1)/2)` is item 344. Theorem C's new content
is the extension to bounded deficit at the CC scale, via Theorem B. That is a modest extension.
For `M=4` it may also follow from Cilleruelo–Granville's classification of four points on
`t R^(1/3)` arcs; I could not verify this against the source.

**The Pell claim.** The write-up checks only `n = 1, 11, ..., 61`, but the summary says "for every
member". The claim is in fact true for every member. In
`four_point_bonus_counterexample.md` (4) the points are `-PABC`, `-P̄AB̄C̄`, `-P̄ĀBC̄` and
`-P̄ĀB̄C`, so

    prod z_i = P P̄^3 a^2 b^2 c^2 = 5 a^2 b^2 c^2 (-3+4i)

identically. The write-up should cite this identity instead of finite evidence.

## 6. Data claims

* "**In exhaustive searches the five-point diameter ratio never falls below 4.9475**" is false as
  an unqualified statement. The primitive window `(-4,-3), (-3,-4), (0,-5), (3,-4), (4,-3)` at
  `N=25` has ratio 4.2024, and another `N=25` window has 4.6985. Both have span over 90°
  (`A<L`), which `search_cg.py` excludes because it loops `A` from `L` upward.
* In the `A>=L` regime, my `N<=10^6` scan agrees with the author's list: 4.9475 at `N=65`, then
  5.0712, 5.1765, 5.3622, ...
* `N=65` has two inequivalent minimal windows, with `L_0=12` and `L_0=42`. This is consistent:
  the author's listing dedupes by `N`.
* The search is exhaustive only over `L_0 <= 120` (or 240) in `6Z`. Theorem B allows `L_0` up to
  `K^10/128 ≈ 8.4*10^6` at `K=8`, so the search covers a small slice. That is fine as data, but
  "exhaustive" should be qualified.
* The cyclotomic numbers reproduce.
* The polynomial-family Levenberg–Marquardt runs are numerical only, and are labelled so.
* The six-point "shortest vector has l1-norm 6" box search is actually complete. If
  `s = sum a ≠ 0`, then `||v||_1 >= 6|s|`. If `s = 0`, then `||v||_1 >= 3||a||_1/2`. So
  vectors with `||v||_1 < 6` lie inside the searched box.

## 7. Corrected statement of what is proved

1. **(A)** For every split `p` with `c_p = (p-1)p^(v_p(L_0))` and `M > 4c_p`:

       (M-4c_p)^3 log p/(12c_p) <= sum_q D_q log q <= MW/4 + binom(M,2) log(C^2/2).

   Hence:
   * if `v_p(L_0) <= B` for **one** split `p`, then `M <= (sqrt(6c_p/log p)+o(1)) sqrt(log R)`;
     for example `5 ∤ L_0` gives `M <= (3.862+o(1)) sqrt(log R)`;
   * if `v_p(L_0) = 0` for all split `p <= M/8`, then
     `M <= (sqrt 24 + o(1)) sqrt(log R/log log R)`. The stated `sqrt 192` is also valid.

2. **(B)** As stated:
   * `Delta_* prod s^2 <= K^(2P) 2^(floor(M^2/4))/4^P`;
   * for five points, `sigma prod s <= K^10/128` and `6 | L_0`;
   * any five lattice points on a circle of radius `R` have diameter at least
     `768^(1/10) R^(2/5) = 1.9433 R^(2/5)`;
   * CG's five-point `O(R^(2/5))` question is equivalent to the fixed-`L ∈ 6Z` clique problem
     with `N <= c(A/L)^(10/3)`.

3. **(C)** For even `M` and all distances at most `K N^(theta_M)`:
   * `prod z_i = eps N_bal^(M/2) m Phi`, with `m` a positive integer and
     `Norm(Phi) <= Delta^2 <= (K^(2P) 2^(M^2/4) 4^(-P))^2`;
   * the cluster lies within `O_{K,M}(N^(theta_M-1/2))` of one of finitely many directions, each
     with `e^(2iM theta)` in `Q(i)`;
   * no such direction is a lattice direction once `N >= N_0(K,M)`.

   The Pell family satisfies `prod z_i = 5a^2b^2c^2(-3+4i)` identically.

## 8. The missing lemma

I agree with the author's diagnosis. Neither theorem gives a lower bound on
`sum_p D_p log p + 2 sum log s_e` that is linear in `W` per point, and that is what the uniform
theorem needs. By (D5), away from `2L` the cotangent integrality is equivalent to the Gaussian
allocation data. So the route cannot by itself supply arithmetic beyond the residue and cut
bookkeeping. For the CG five-point question, the open step is finiteness of the bounded-`L` rigid
system (B4)/(B5). It has varying `S`-unit-type supports, so neither Evertse-type counts nor `abc`
applies directly, as the author notes.
