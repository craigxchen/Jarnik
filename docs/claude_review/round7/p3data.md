# Round 7 data: genuine residue demands and the discreteness question (Q)

Status: draft, unrefereed; all sections filled. Companion to `round7/p3lower.md`; same conventions.
A finishing session added §0, §4.2, Prop. 4.1, Remark 4.2 and Example 6.2A, and the `r2/` checks.
It also corrected the counting remark of §6.2 and the ILS calibration sentence of §3, and withdrew the
claim that discreteness is excluded on small-`M` evidence (§6.3).

*Labels.* *Proposition/Lemma/Theorem/Remark/Example* means proved here with all quantifiers.
*Evidence* means a computation that no proof uses. *Heuristic* means an argument with a named
unproved assumption. Solver upper bounds are evidence unless marked optimal (`*`).

## 0. Verdict

**Growth: no proof of `M = o(log R/loglog R)`, and no new barrier theorem for actual circles.** The
discreteness question (Q) is **not decided**.

The round locates the `log^2 M/loglog M` effect of the uniform-margin fakes.
* It is not in the pair LP: adapted real weights cost the Perron bound.
* It is not in separation, and not in equal-weight granularity at the Plotkin level.
* It does return in a proved model when the residues are independent of the profile (coupling).
* Whether adapted near-tight `+-1` designs with `O(M)` prime columns exist for genuine demands is the
  open discrete core of (Q). The heuristic upper bounds found here rise relative to the Perron bound
  with `M`. But the solvers are shown to miss structured codes that exist, so this is no evidence
  either way.

An earlier version of this draft said discreteness was excluded "on small-`M` evidence". That claim
is **withdrawn** (§4.2, §6.3).

| # | statement | class | type |
|---|---|---|---|
| 1 | Genuine unit-1 demands of the Paley profiles with actual primes match the independent residue model in every statistic tested. Mean `~ 2 log X`. `lambda_max(D) = (1.5-2.0) M log M` and `max_x D_x = (2.0-2.8) M log M`, both flat in `M`. `D_max` grows faster than `log M` (`D_max/H` from 2.2 to 3.7, `H = log^2 M/loglog M`, `M <= 84`). Uniform margins pay `(M-1) D_max/lambda_max` times the Perron bound, from 3.6 to 5.6 (§2) | data of actual lattice points (P3-gen members) | exact per trial; aggregates are evidence |
| 2 | Identity (6.1): `W - 2 lambda_max = imbalance + u-weighted slack` for every cut measure and every profile. The real-weight optimum equals `2 lambda_max(D)` iff `2(lambda_max I - D)` is a combination of `u`-balanced cut matrices (Lemma 6.3) | any `D >= 0` | proved |
| 3 | For genuine `D`, `W_cut = 2 lambda_max(D) + O(1)`: excess 0-15, ratio `<= 1.014` for `M >= 20`. From `M = 12` on, the optimal margins are exactly tight and the excess is pure imbalance (§6.2) | (Q), genuine `D`, `M <= 24` | certified LP values (evidence about `M -> infinity`) |
| 4 | Equal-weight granularity: every equal-weight profile has `W >= LB_eq` (Prop. 4.1). On the 120 genuine instances `LB_eq/2 lambda_max` has median 1.12-1.22 for every `(M, X)` with `M = 8..32`: a flat factor | (Q), equal weights | exact bound; flatness is evidence |
| 5 | Ties: two `D = 0` pairs whose margins cancel identically make a sign matrix unrealisable by distinct primes (Remark 4.2) | (Q) with actual primes | proved |
| 5a | One collision layer (a partition into equal classes of a Hadamard order or size 2) is met at the prime floor with `W <= 2 lambda_max(D) + 2 M w_0` (Example 6.2A). So at the floor the difficulty of genuine `D` is merging its `~ X/log X` layers into `O(M)` shared columns | any such `D` | proved |
| 6 | Best profiles found at the prime floor, as median multiples of `2 lambda_max` (`X = 2M` / `4M`): 1.32 / 1.22 at `M = 8` (7 of 12 proved optimal), 1.45 / 1.37 at `M = 12`, 1.59 / 1.64 at `M = 20`, 1.82 / 1.72 at `M = 24`. They use `1.5M-4.7M` columns. Uniform margins: 3.2 / 2.9, 3.5 / 3.7, 4.3 / 4.2, 5.3 / 6.3. So uniform costs 2.3-2.7 times the best found for `M <= 20`, and 2.9 / 3.7 at `M = 24`. Verified actual-prime profiles: 1.7, 1.8 / 1.5, 2.3, 2.7 / 2.3, i.e. 1.8-2.4 times cheaper than the genuine Paley fakes of §2.4. All 40 realised profiles at `M >= 12` pass the exact pair and support-`<= 4` checks, and their least character margin (`|c_x| <= 12`) is at a pair (§4.2, Tables 4.2a/b) | (Q), genuine `D` | upper bounds (evidence); prime profiles verified exactly |
| 7 | Calibration: local search finds Plotkin-tight codes known to exist (stacked Paley) at `M = 8`, sometimes at `M = 12`, never at `M = 20, 24` (§4.2). So item 6 is no evidence about the optimum for `M >= 20` | solver behaviour | evidence |
| 8 | Separation costs `O(M)` at real weights and `O(M log M)` at the prime floor (Prop. 6.2); a term of order `M` is unavoidable (Prop. 3.2) | (Q) / all members of `A_X` | proved |
| 9 | Margins and character margins are additive under concatenation (Lemma 6.1). Character margins of a fitted profile cannot be repaired cheaply that way | any profiles | proved |
| 10 | Against residues uniform relative to a fixed profile, every profile with `W <= (2 - eps) M log^2 M/loglog M` satisfies all pair inequalities with probability `<= 4 M^(-eps/16)`. This covers adapted profiles, any weights and any number of columns (Thm 5.1, Cor. 5.2) | model `U_X` | proved |
| 11 | The fitted `M = 8` profiles are infeasible for their own genuine residues (0 of 400 conjugation choices each). At equal `W`, uniform margins beat Paley profiles, which beat adapted ones (§5.1) | P3-gen data | exact per trial; evidence |
| 12 | Prime selection (choosing among `K >= 4` primes per column) is not excluded by the random-residue method (§6.3 remark) | -- | proved limitation of the method |
| 13 | Counting decides nothing. The demand-side count allows `O(M)` columns. The design-side (local-CLT) count forbids random-like designs only for `M > 2^(8K)`, and the same count wrongly forbids Hadamard matrices (§6.2) | -- | heuristic |

**What this means for the conjecture.**
* *(Q)-order* (`D` fixed, profile fitted). At real weights adaptivity absorbs the
  `log^2 M/loglog M` effect (proved, `p3lower.md` Thm A; item 3). With prime columns the question
  becomes the existence of `+-1` designs whose pair distances follow `D` to within `O(1)` columns on
  most pairs, with `O(M)` columns. That is a structured-design problem, which neither the computations
  nor counting decide.
* *P3-gen order* (`D` comes from the profile's own primes). Independent residues force
  `(2 - o(1)) M log^2 M/loglog M` (item 10).
* *A growth proof through pair data* needs an input about **actual** circles in either order:
  * (a) their prime residues cannot follow their exponent profile (an `H`-type hypothesis for adapted
    profiles); or
  * (b) a prime-floor lower bound for (Q) that holds for every demand matrix an actual cluster can
    produce.

  Route (b) also needs an `H`-type input: structured demand matrices, such as a single collision
  layer, are cheap (Example 6.2A), so a lower bound can only hold for "generic" `D`.
* *A construction at `W = O(M log M)`* needs both an adapted structured design and prime selection
  that steers the residues modulo many `q ~ M` at once.
* Neither direction is settled.

## 1. Setting, conventions and what is exact

Notation is that of `p3lower.md` §1:
* `W = log N` and `w_j = log p_j`;
* `mu_xy = d_xy - W/2` is the pair margin;
* `m_xy` is the **unit-1** pair collision modulus at odd prime powers `Q = q^a <= X`, and
  `D_xy = 2 log m_xy` is the pair demand.

Throughout `C = sqrt 2`, so `c_C = 2 log C - log 2 = 0`. Any other fixed `C` shifts every pair
constraint by the same constant `c_C`.

**Genuine data.** The profiles are the fully flipped Paley profiles:
* `construct`: `round5/construct.md` §2, `b = 8`, `r = 9M - 8`;
* `lower`: `round5/lower.md` §3.1 = `sharp.md` §4, `b = 33`, `r = 33(M-1)`.

`M = q_0 + 1` with `q_0 = 3 mod 4` prime. The primes are actual primes `p_j = 1 mod 4`:
* *window*: consecutive primes from a random start in `[10^7, 10^9]`;
* *iid*: independent random primes in `[10^8, 10^10]`.

Each prime gets a random choice between `pi_j` and `conj(pi_j)`, and the windings are `k = 0`. The
points `zhat_x = prod_j pi_j^(a_xj) conj(pi_j)^(1 - a_xj)` are then actual lattice points of
`x^2 + y^2 = N`. Their residues are genuine, so the data are those of a member of P3-gen
(`p3lower.md` Def. 1.1). For every pair and every `Q <= 4M`, collisions are computed exactly by
discrete logarithms in `T_Q` (`demands7.py`).

`verify_direct.py` (round-7 partial work, re-run here) confirms that computation against direct
Gaussian arithmetic:
* modulo every `Q <= 4M`, for all ordered pairs and all four units, at `M = 12, 20, 24, 32`;
* in exact big-integer `Z[i]` arithmetic (`N` up to 5561 digits) at `M = 12, 20`.

**Which demand enters the margin.** In (e-char) with `eta != 1` the left side is
`|sin((delta_x - delta_y - arg eta)/2)| >= sin(pi/4) - Delta/2`. So a collision with another unit
only forces `d_xy >= 2 log m_(xy,eta) - log 2 - o(1)`, which is `W/2 + mu_xy >= ...`. That is
automatic once `W >= 4 log m_(xy,eta)`. Only the unit-1 moduli constrain margins, so `D` always
means the unit-1 demand. Any-unit moduli are reported in §2 because the task asked for them.

**Inputs used, and whether each is valid for actual circles.**

| input | statement | valid for actual circles? |
|---|---|---|
| pair inequality | `mu_xy >= D_xy - c_C` | yes (`p3lower.md` Lemma 1.3, Prop. 1.2) |
| cut representation | margins are cut combinations with total weight `W`; each cut weight is a sum of `log p_j` | yes (`p3lower.md` Lemma 1.4) |
| prime floor | if every `p_j > X`, every cut weight is `>= w_0(X) := log p_min(X)`, the least prime `= 1 mod 4` above `X` | yes, in the subclass `p_j > X`, which is the setting of every round-5/6 fake |
| `delta`-placement | Lemma 1.1 below | yes (necessary direction) |
| character margins `s_c := w(v(c)) - W/2 >= 0` | imposed as the task asks | **stronger** than valid. The valid condition is `s_c >= 2 log m_c - 2 log(n C/2^(3/2))`, `n = ||c||_1`, i.e. `s_c >= -2 log(n/2)` when `m_c = 1`. So it is imposed only on the construction (upper-bound) side |

**Lemma 1.1 (`delta`-placement).** Let `C = sqrt 2`, a squarefree profile and moduli `m_xy` be
given. Write `Pi_sep(x,y) = prod_(j: a_xj != a_yj) p_j` and `Pi_same(x,y) = prod_(j: a_xj = a_yj) p_j`.

1. *(necessary)* If `delta in [0, Delta]^M` satisfies the unit-1 pair inequalities (1.1) of
   `sharp.md`, then `t = delta/Delta in [0,1]^M` satisfies `|t_x - t_y| >= m_xy e^(-mu_xy/2)` for
   all `x != y`.
2. *(sufficient, exact)* If `t in [0,1]^M` satisfies, for all `x != y`,
   `|t_x - t_y|^4 (Pi_sep - m_xy^2) >= m_xy^4 Pi_same`, then `delta = Delta t` satisfies (1.1) for
   every pair, with `eta = 1`.

*Proof.*
1. `|sin(v/2)| <= |v|/2`. So `|delta_x - delta_y| >= sqrt2 m e^(-d/2)`, and dividing by
   `Delta = sqrt2 e^(-W/4)` gives `m e^(-(d - W/2)/2) = m e^(-mu/2)`.
2. Put `u = 2^(-1/2) m e^(-d/2)`, so `u^2 = m^2/(2 Pi_sep)`; the hypothesis forces `m^2 < Pi_sep`, so
   `u < 1`. For `|v| <= Delta <= pi`, `|sin(v/2)| >= u` holds as soon as `|v| >= 2 arcsin u`. Since
   `arcsin u <= u (1 - u^2)^(-1/2)` and `2u/Delta = m e^(-mu/2)`, it suffices that
   `|t_x - t_y|^2 (1 - u^2) >= m^2 e^(-mu)`. Here `e^(-2mu) = Pi_same/Pi_sep`. Square this, and use
   `(1 - a/2)^2 >= 1 - a` with `a = m^2/Pi_sep`. QED

So, with positions fixed, the placement condition is a set of exact integer inequalities. Without
fixed positions it is a min-span ordering problem, which is solved by a MILP over order variables
(`verify_profile.placement_milp`).

*Remark.* The necessary direction is a Hamiltonian-path condition: some ordering has
`sum_(consecutive) m e^(-mu/2) <= 1`, so `sum_path (mu - D) >= 2(M-1) log(M-1)`. At real weights this
costs only `O(M)` (§4). It is not part of `p3lower.md`'s pair LP.

## 2. Genuine pair demands of the fully flipped Paley profiles

*Class: data of actual lattice points (members of P3-gen). Type: exact per-trial integers, aggregated
statistics are evidence.* There are 4000 trials (`genuine_stats.py 300 100 100`):
* for each `M in {12, 20, 24, 32, 44, 60, 72, 84}`: 300 `construct`/window, 100 `lower`/window and
  100 `construct`/iid trials;
* 12 more `construct` trials at `M = 8` for §4.

Notation: `L = log M`, `H = log^2 M/loglog M`. "Indep." is the independent uniform residue model:
classes uniform on `Z/m_Q`, nested in the level, independent over `q` and `x`. Its exact mean is
`E D = 2 sum_(Q <= X) log q/m_Q` (`demands7.expected_mean_independent`).

**2.1 Unit-1 demands** (`part1_compact_out.txt`; all 500 trials per `M` pooled, since the three
families agree, see 2.2):

| M | X | E D (exact, indep.) | mean D | D_max median [5%,95%] | D_max/L | D_max/H | lambda_max/(M L) | max_x D_x/(M L) | (M-1)D_max/lambda_max |
|---|---|---|---|---|---|---|---|---|---|
| 12 | 24 | 3.761 | 3.782 | 14.8 [11.1, 20.2] | 5.95 | 2.18 | 1.509 | 2.044 | 3.64 |
| 12 | 48 | 5.155 | 5.186 | 19.6 [14.7, 26.9] | 7.87 | 2.88 | 2.050 | 2.773 | 3.54 |
| 20 | 40 | 4.638 | 4.637 | 20.9 [16.8, 26.3] | 6.97 | 2.55 | 1.542 | 2.074 | 4.28 |
| 20 | 80 | 6.121 | 6.112 | 26.7 [21.5, 34.0] | 8.91 | 3.26 | 2.020 | 2.707 | 4.18 |
| 24 | 48 | 5.155 | 5.141 | 23.7 [19.6, 30.3] | 7.45 | 2.71 | 1.610 | 2.135 | 4.46 |
| 24 | 96 | 6.349 | 6.339 | 28.9 [23.6, 36.9] | 9.08 | 3.30 | 1.978 | 2.615 | 4.35 |
| 32 | 64 | 5.651 | 5.645 | 27.3 [23.2, 34.3] | 7.88 | 2.83 | 1.624 | 2.136 | 4.73 |
| 32 | 128 | 7.027 | 7.021 | 33.3 [28.4, 41.2] | 9.61 | 3.45 | 2.015 | 2.642 | 4.64 |
| 44 | 88 | 6.247 | 6.249 | 31.9 [27.8, 39.6] | 8.43 | 2.97 | 1.647 | 2.121 | 5.06 |
| 44 | 176 | 7.658 | 7.672 | 38.9 [32.8, 48.2] | 10.29 | 3.62 | 2.019 | 2.591 | 4.97 |
| 60 | 120 | 6.883 | 6.881 | 37.3 [32.3, 45.0] | 9.11 | 3.14 | 1.676 | 2.114 | 5.33 |
| 60 | 240 | 8.277 | 8.281 | 43.8 [37.6, 53.1] | 10.70 | 3.68 | 2.016 | 2.536 | 5.21 |
| 72 | 144 | 7.244 | 7.239 | 39.9 [35.2, 48.7] | 9.33 | 3.17 | 1.689 | 2.115 | 5.46 |
| 72 | 288 | 8.662 | 8.662 | 47.0 [40.9, 55.6] | 10.99 | 3.74 | 2.019 | 2.513 | 5.36 |
| 84 | 168 | 7.566 | 7.571 | 42.7 [37.3, 52.5] | 9.64 | 3.24 | 1.705 | 2.104 | 5.61 |
| 84 | 336 | 8.904 | 8.913 | 49.3 [43.2, 59.2] | 11.13 | 3.74 | 2.006 | 2.464 | 5.49 |

Reading the table:
* **Mean.** The mean demand equals the exact independent expectation to within `0.5%`. It is
  `~ 2 log X - O(1)`, i.e. `Theta(log M)`, as the pigeonhole of `two_thirds` forces.
* **Perron root.** `lambda_max(D) = (1.5-2.0) M log M` and is flat in `M`. Its Perron vector is
  nearly uniform, and `lambda_max ~ (M-1) * mean`.
* **Point totals.** `max_x D_x = (2.0-2.8) M log M`, also flat. There are no irregular point
  totals, as `p3lower.md` Cor. A.3 would need for a pair-cut growth argument.
* **Maximal demand.** `D_max` grows faster than `log M`: `D_max/L` goes from 6.0 to 9.6 at `X = 2M`.
  Its ratio to `H` also still rises slowly, from 2.2 to 3.2 at `X = 2M` and from 2.9 to 3.7 at
  `X = 4M`. That is consistent with the Poisson heuristic `D_max ~ 4(1 - o(1)) H`, approached
  slowly. These `M` cannot separate `log M` from `H` growth by fitting.
* **Uniform-margin penalty.** `(M-1) D_max/lambda_max` grows from 3.6 to 5.6. This is the ratio of
  the uniform-margin cost `2(M-1) D_max` to the Perron bound `2 lambda_max`, i.e. the size of the
  "`log^2 M/loglog M` effect" at these `M`. (`L/loglog L` only grows from 2.7 to 3.0 over the same
  range. The faster growth here comes from `D_max/H` itself still rising.)

**2.2 Genuine = independent (evidence).** The genuine families agree with each other and with
simulations of the independent model (100 per row). This holds for every statistic above, to
within sampling error (`genuine_table_out.txt`):
* `construct` vs `lower` (different profiles, `r` up to 2739);
* window vs iid primes.

Tail counts, i.e. the mean number of pairs per trial with `D >= tL`, genuine/independent, at
`M = 84`, `X = 168`:

| t | 2 | 3 | 4 | 5 | 6 | 8 |
|---|---|---|---|---|---|---|
| genuine / indep. | 1287/1279 | 662/658 | 302/301 | 126/125 | 46.7/46.2 | 5.2/4.8 |

Window collision counts per pair (level 1, primes in `(sqrt M, X]`) match the exact
Poisson-binomial law of independent `Bernoulli(1/m_q)` to within sampling error. At `M = 84`,
`X = 336` the frequencies of `c = 0..6` are
`.420 .371 .157 .0424 .0084 .00118 .00014`, against the exact law
`.420 .371 .157 .0422 .0082 .00122 .00014`. (Plain Poisson overstates the tail: `.0016` at `c = 5`.)
Forced Legendre parity is present in the genuine data (`p3lower.md` Prop. B.3). With random primes
it produces no visible bias.

**2.3 Any-unit demands** (`eta` in all four units). The mean equals the exact independent
expectation `2 sum log q * 4/m_Q`. That includes the forced collisions at `q = 3, 5` (`m = 4`):
every pair collides there with some unit. Values, `X = 2M` / `X = 4M`:

| M | 12 | 24 | 44 | 84 |
|---|---|---|---|---|
| mean `D_any` | 15.0 / 20.7 | 20.6 / 25.4 | 25.0 / 30.6 | 30.3 / 35.6 |
| median max `D_any` | 29.5 / 41.9 | 47.4 / 59.2 | 63.8 / 78.8 | 83.4 / 99.2 |

As explained in §1, these do not constrain margins.

**2.4 The Paley fakes with their own primes** (`paley_selfcheck.py`; 60 `construct` and 20 `lower`
trials per `M`, primes `~ 4*10^8`, `tau = 19.7`):
* every trial at every `M <= 84` satisfies all unit-1 pair inequalities at `X = 2M` and `X = 4M`;
* but the weight scale they need, `tau_req = max_xy D_xy/nu_xy` with `nu = mu/tau` (the normalised
  margin), grows.

| M | 12 | 20 | 24 | 32 | 44 | 60 | 72 | 84 |
|---|---|---|---|---|---|---|---|---|
| `construct`, `X=2M`: median `tau_req` | 3.67 | 4.97 | 5.84 | 6.68 | 8.11 | 8.94 | 9.82 | 10.46 |
| `W_req = r tau_req` | 367 | 854 | 1215 | 1871 | 3147 | 4758 | 6285 | 7823 |
| `W_req / (2 lambda_max)` (mean `lambda_max` over the trials of 2.1) | 4.1 | 4.6 | 4.9 | 5.2 | 5.7 | 5.8 | 6.0 | 6.2 |
| `lower`, `X=2M`: `W_req` | 305 | 786 | 1154 | 1728 | 2724 | 4289 | 5878 | 7614 |

So with genuine residues the uniform-margin Paley fakes pay a factor growing from 4 to 6 over the
Perron bound in this range. At `tau = 19.7` they are still feasible for `M <= 84`, because
`tau_req <= 12`.

## 3. The discrete problem (Q): formulation

Fix a symmetric demand matrix `D >= 0`. In §4 it is the genuine unit-1 demand of a `construct`
profile with actual primes. "Cut" means a nonconstant `s in {+-1}^M` modulo `s ~ -s`. With a weight
`lambda_s >= 0` on every cut, `W = sum lambda_s` and
`mu_xy = sum_s lambda_s([s_x != s_y] - 1/2)` (`p3lower.md` Lemma 1.4). The problems are:

| name | constraints | role |
|---|---|---|
| `2 lambda_max(D)` | -- | Perron lower bound for all of the below (`p3lower.md` Thm A(2)) |
| `W_cut` | `mu >= D` | real weights (LP) |
| `W_sc` | `mu >= D`, `lambda_s in {0} u [w_0, inf)`, `w_0 = log p_min(X)` | **prime floor**: the relaxation for profiles with all `p_j > X` |
| `W_prime` | `mu >= D`; profile of actual distinct primes `p_j = 1 mod 4`, `p_j > X`, `e = 1` | the discrete problem itself; every upper bound below is such a profile, verified exactly |
| `+sep` | in addition, a `delta`-placement (Lemma 1.1) | the full unit-1 pair system |
| `+char` | all `s_c >= 0` | as the task asks (stronger than valid, §1) |
| `W_unif` | `mu_xy >= D_max` for all pairs | uniform margins. With real weights `W_unif = 2(M-1) D_max` exactly for even `M` (uniform measure on bisections; Perron for the lower bound). This is the cost of the Paley-type fakes up to their factor `b/(b-4)` |

`r` is the number of physical columns (primes). In `W_sc`, `r <= W/w_0`, so `W = O(M log M)` forces
`r = O(M)`, the regime of the task.

**Proposition 3.1 (the prime floor is a valid relaxation).** Let a profile have distinct primes
`p_j > X` and any exponents. Then its cut weights `lambda^prof` (`p3lower.md` Lemma 1.4) satisfy
`lambda_s in {0} u [w_0, inf)`. Hence every member of `A_X(M, C)` with all `p_j > X` and pair
demands `>= D` has `W >= W_sc(D - c_C)`.

*Proof.* Each threshold cut `s^(j,t)` carries weight `log p_j >= w_0`. Several `(j,t)` may give the
same cut; their weights add. QED

**Proposition 3.2 (separation costs at least about `4M`; valid for actual circles).** Every member
of `A_X(M, C)` satisfies

```text
W  >=  (4/M) [ sum_(x<y) D_xy  +  2 S_M ]  -  2(M-1) c_C,
S_M := sum_(1 <= k < M/2) (M-k) log((M-k)/k)  =  (1/4 + (log 2)/2 + o(1)) M^2,
```

where `D` are the member's unit-1 pair demands. The term `(8/M) S_M` is 22.8, 40.2, 76.3, 94.6,
131.6 at `M = 8, 12, 20, 24, 32`, and `(2 + 4 log 2 + o(1)) M = (4.77 + o(1)) M` in general (4.52 M
at `M = 100`). Actual clusters are members (`p3lower.md` Prop. 1.2).

*Proof.*
1. *Placement.* By Lemma 1.1(1), for general `C`,
   `mu_xy >= D_xy + 2 log(1/|t_x - t_y|) - c_C` for the member's own `t = delta/Delta in [0,1]^M`.
2. *Sum over pairs.* Sort `t_(1) <= ... <= t_(M)`. For each `k`, the `M-k` gaps
   `t_(i+k) - t_(i)` sum to at most `k`, since every point lies in at most `k` of the windows. So by
   AM-GM, `sum_i log(1/(t_(i+k) - t_(i))) >= (M-k) log((M-k)/k)`. Every term is `>= 0`, since gaps
   are `<= 1`, so we may keep only `k < M/2`.
3. *Cut identity.* Sum over all pairs and use `sum_(x<y) mu_xy <= M W/4` (`p3lower.md` Lemma 1.4(4)).
4. *Asymptotics.* `S_M ~ M^2 int_0^(1/2) (1-u) log((1-u)/u) du`, and the integral equals
   `1/4 + (log 2)/2 = 0.5966`. QED

*Remarks.*
1. With `D` replaced by its pigeonhole lower bound this is the `two_thirds` mechanism plus `~4M`.
   It is a lower-order improvement (`O(M)` against `3 M log M`) and does not change the constant 3.
2. Item 2 is sharp up to the constant. For `sum_(x<y) log(1/|t_x - t_y|)`:
   * the bound gives `0.597 M^2`;
   * the true minimum (Fekete points, logarithmic capacity `1/4` of `[0,1]`) is `~ (log 2) M^2 = 0.693 M^2`;
   * equal spacing gives `~ (3/4) M^2`.
   So at real weights separation costs between `4.77 M` and `~ 6 M`.
3. What matters for (Q): uniform margins have slack `D_max - D_xy`, which is large for most pairs,
   so they pay essentially **nothing** for separation (§4: `W_unif+sep <= W_unif + 0.5` in every
   instance). Adapted margins are tight and pay the full `Theta(M)`.

**How the solvers work** (`qsolve.py`, `q_run.py`).
* `W_cut`: full LP (`M <= 12`, all `2^(M-1) - 1` cuts) or column generation with exact pricing
  (`M = 20, 24`, Gray-code max-cut enumeration in C, `maxcut.c`, checked against brute force). The
  dual `nu` gives the certified bound `W_cut >= <nu, D>/surplus(nu)`, with `surplus` computed exactly
  by the same enumeration (Theorem A(1)).
* `W_sc`:
  * `M = 8`: MILP over all 127 cuts with indicator binaries, HiGHS, proved optimal where it says so;
  * `M >= 12`: the same MILP restricted to candidate cuts (large-neighbourhood search over LP support
    and low-reduced-cost cuts). This gives an upper bound only; the lower bound is `W_cut`.
* `W_prime`: from a `W_sc` solution, each support cut takes one actual prime from a short menu near
  its weight, all primes distinct (MILP). If that fails, the weight floor is raised and `W_sc` is
  re-solved. Every reported `W_prime` profile is re-checked in integer arithmetic (`verify_profile.py`):
  * distinct primes `= 1 mod 4`, `> X`;
  * nonconstant columns and distinct rows;
  * every pair inequality `prod_sep p >= m^4 prod_same p`.
* `+sep`:
  * positions are fixed first: the min-span placement of the `W_prime` solution, rescaled to
    `[0,1]`, or equal spacing in its order, whichever has the smaller LP value;
  * demands become `D_xy + 2 log(1/|t_x - t_y|) + 10^(-6)`, and the pipeline is re-run;
  * the final profile is verified against Lemma 1.1(2) exactly, with rational `t`.
* `+char`:
  * every ternary `c` of support 4 and every `2e_x - e_y - e_z` is checked exactly
    (`prod p^(2|v_j|) >= N`);
  * the minimum of `s_c` over all zero-sum `c` with `|c_x| <= 12` is computed by `M` MILPs (floating
    point). Pairs are included in that minimum.
* `W_unif + sep`: least uniform margin with a feasible placement (placement MILP, bisection).

At `M = 12, 20, 24` the full pipeline is too slow. The first `M = 12` run of `q_run.py` took 42
minutes per instance and found no actual-prime profile (`superseded/q_results_M12_firstrun.jsonl`).
The second crashed in the nearest-prime rounding: its fixed search window could be exhausted, a bug
since fixed in `q_run.round_nearest`. The lean driver `q_lean.py` is used there instead:
* `W_cut` as above, with certified dual bound;
* `W_sc` upper bound: reweighted `l^1` (LPs with costs `1 + beta/(lambda_s + eps)`, which push the
  weight onto few cuts), then the restricted semi-continuous MILP on that support plus the 150 cuts of
  best reduced cost (`rwl1.py`). LNS rounds (`qsolve.lns`) did not improve these solutions in tests
  and were switched off;
* `W_prime`: at the floor `w_0` the optimum puts many cuts on `w_0` itself, but only one prime has
  that weight. So the floor is raised to `log p_k`, `k = n/2` or `n` (`n` = number of `W_sc` cuts,
  `p_k` the `k`-th prime `= 1 mod 4` above `X`). The demands are shifted by
  `delta = 0.3, 1, 2`, the `rwl1` step is repeated, and each cut is rounded to the nearest unused
  prime (or product of two). The cheaper of this and the nearest-prime menu MILP with the exact
  demands is kept. Every profile is re-verified in integers (`check_pairs`). Also reported:
  `stair(r) = sum_(i <= r) log p_i`, the least possible weight of any `r` distinct primes
  `= 1 mod 4` above `X`;
* characters: support `<= 4` in floating point (`char_small_float`);
* no separation on the discrete profiles. Prop. 6.2 bounds its cost rigorously, and the LP value
  with equal spacing in a greedy order (`W_cut+sep`) is reported for scale.

The `rwl1` bounds degrade quickly with `M`: at `M = 20` they are 2.3-3.0 times `W_cut`. That is a
failure of the heuristic: the restricted cut pool from column generation is a poor source of sparse
solutions. So a second heuristic was added, which works with the sign matrix directly
(`fls.py`, `ils_run.py`):
* *flip local search with iterated kicks (ILS):* fix `r` columns `S in {+-1}^(M x r)`. For fixed
  `S` the best weights `>= w_0` solve an LP (`lp_S`, with an `l^1` violation penalty). Single
  entries of `S` are flipped in the order of their first-order gain from the LP duals, and a flip is
  accepted if the re-solved LP value drops. When no flip helps, 3 random flips are applied to the
  best `S` and the search restarts, 60 times. At the end, columns are dropped while that helps.
  Starts: two and three stacked copies of the Sylvester design restricted to the `M` points;
* *calibration:* at `M = 8` ILS (as rerun in `r2/ils2_run.py`, which adds Paley starts) is within
  0-8.5% (median 2.7%) of the six proved MILP optima of §4.1 (`r2/ils2_M8_log.txt`). An earlier
  version of this sentence cited `ils_M8_log.txt` for "1-11%, median 5%". That file holds a single
  non-optimal instance followed by the crash, so the claim was unsupported and is replaced;
* *realisation:* `S` is kept, and the columns, in increasing order of weight, get the smallest unused
  prime `= 1 mod 4` above `X` with `log p >= w_j - 0.02`. The LP is re-solved with these per-column
  floors until the assignment is stable, with demands `D + delta`, `delta = 0.25, 0.75, 1.5, 3`.
  The result is verified exactly (`check_pairs`).

## 4. Exact and certified results for (Q) at small M

*Class: (Q) with a fixed genuine `D` (§3). Type: LP values with certified dual bounds; MILP optima
where marked `*`; every `W_prime` profile verified exactly (pairs, placement, characters of support
`<= 4`). These are evidence about (Q), not statements about actual circles.*

Instances: genuine unit-1 `D` of `construct` with actual window primes, from the trial files of §2.
`M = 8` uses 6 trials at `X = 16` and 6 at `X = 32`. The floors are `w_0 = log 17` and `log 37`.

### 4.1 `M = 8` (`q_table_out.txt`)

Columns:
* `W_sc [lb, ub]`: `*` marks a proved MILP optimum;
* `W_prime (r)`: verified actual-prime profile, with `r` physical columns;
* `sep LB`: Proposition 3.2 (valid for actual circles);
* `W_cut+sep`: LP with the chosen positions;
* `char min`: the minimum of `s_c` over all `c` with `|c_x| <= 12`, by MILP. Every profile also
  passes the exact support-`<= 4` check, with a minimum non-pair margin of `>= 2.5`.

| X | t | D_max | 2 lam | W_cut | W_sc [lb, ub] (cuts) | W_prime (r) | sep LB | W_cut+sep | W_prime+sep (r) | W_unif | char min |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 16 | 0 | 9.2 | 42.3 | 44.82 | [51.0, 56.8] (17) | 74.3 (17) | 62.2 | 72.9 | 101.2 (22) | 128.7 | 0.42 |
| 16 | 1 | 10.5 | 40.3 | 43.74 | 52.6* (15) | 69.5 (16) | 60.3 | 70.3 | 84.9 (19) | 147.6 | 1.04 |
| 16 | 2 | 9.5 | 40.8 | 44.25 | [48.0, 52.9] (16) | 69.5 (16) | 60.1 | 72.0 | 102.1 (22) | 133.3 | 0.34 |
| 16 | 3 | 9.5 | 52.6 | 55.30 | [57.6, 65.1] (18) | 80.0 (18) | 70.8 | 85.3 | 97.3 (21) | 133.3 | 0.34 |
| 16 | 4 | 8.0 | 36.5 | 40.55 | 48.5* (14) | 84.9 (19) | 54.4 | 68.1 | 101.6 (22) | 112.2 | 0.69 |
| 16 | 5 | 9.0 | 27.5 | 33.27 | 41.8* (12) | 69.1 (16) | 48.4 | 57.1 | 102.8 (22) | 126.3 | 1.50 |
| 32 | 0 | 13.7 | 73.4 | 83.16 | 86.3* (18) | 115.3 (23) | 90.9 | 110.9 | 116.9 (23) | 191.5 | 2.51 |
| 32 | 1 | 12.8 | 66.9 | 70.01 | [75.6, 80.8] (17) | 115.4 (23) | 87.0 | 101.8 | 142.2 (24) | 178.9 | 2.27 |
| 32 | 2 | 13.9 | 62.0 | 67.59 | 76.4* (17) | 134.1 (24) | 73.2 | 93.5 | 131.0 (24) | 194.6 | 0.81 |
| 32 | 3 | 15.8 | 88.1 | 94.73 | [94.9, 97.8] (20) | 152.9 (29) | 101.3 | 119.4 | 154.2 (26) | 221.1 | 3.32 |
| 32 | 4 | 13.9 | 57.3 | 59.21 | [67.1, 73.4] (13) | 120.4 (23) | 74.3 | 89.3 | 147.2 (29) | 194.6 | 2.38 |
| 32 | 5 | 15.8 | 47.9 | 62.91 | 62.9* (15) | 70.0 (15) | 59.2 | 73.5 | 110.4 (22) | 220.6 | 0.16 |

Uniform margins pay nothing for separation in any instance: `W_unif+sep <= W_unif + 0.5`.

**Observations at `M = 8`.**
1. *Real weights.* `W_cut/(2 lambda_max) = 1.00-1.31`, median 1.08. Its certified dual bound
   equals `W_cut` to 4 digits.
2. *Prime floor.*
   * `W_sc/W_cut` has median 1.20 (`X = 16`) and 1.08 (`X = 32`); 7 of 12 values are proved optimal.
   * The optima use **12-20 cuts**, i.e. `1.5M-2.5M` columns, against `binom(8,2) = 28` cuts in the
     LP.
   * The minimum weight per column costs only 8-20%.
3. *Actual distinct primes.* These add 30-75% at this tiny `X`: the primes `= 1 mod 4` just above
   16 or 32 are 17, 29, 37, 41, 53, ..., so 20 distinct ones have mean log about 4.5-5, far above
   `w_0`. This is a small-`X` effect. For `r = O(M)` primes above `X = 2M` it costs
   `log(1 + O(log M)) = O(loglog M)` per column, which is lower order.
4. *Separation.*
   * Adapted profiles are tight. Their margins correlate with `D` (correlation 0.72-0.91), with
     max/mean margin 2.4-4 and up to half the pairs within 0.5 of their demand.
   * So they pay the separation term: `W_cut+sep - W_cut = 11-32`, median 27.6. Compare the term
     `(8/M) S_8 = 22.8` that Prop. 3.2 adds to the averaged bound `(4/M) sum D`; it is not added to
     `W_cut` itself, which is why one instance gains only 11.
   * The no-separation prime profiles all have minimal placement span 1.86-2.73 `> 1`, i.e. they
     are **not** realisable by angles.
5. *Uniform margins* cost `W_unif/(2 lambda) = 2.5-4.6` (medians 3.2 and 2.9). That is `1.10-2.0`
   times the verified adapted actual-prime profiles with separation (medians 1.29 at `X = 16`, 1.46
   at `X = 32`).

### 4.2 `M = 12, 20, 24` (finishing session; `p3data_checks/r2/`)

*Class: (Q) with fixed genuine `D`. Type: `W_cut` certified (largest relative gap to the dual bound
`1.9e-10`); `LB_eq` exact; every other `W` is a heuristic upper bound, so evidence. Every actual-prime profile is
verified exactly (`check_pairs`).*

**What was run.**
* `W_cut`: certified, from `perron_check.json`, `q_lean_M12.jsonl` and `r2/wcut_extra.py` (column
  generation with exact max-cut pricing, for `M = 24`, `X = 96`).
* `LB_eq` (exact; `r2/eqlb.py`): a lower bound for **equal-weight** profiles only (Proposition 4.1
  below). It is not a lower bound for `W_sc`.
* `W_ils` (`r2/ils2_run.py`): the ILS of §3 from four starts (two and three stacked copies of the
  Sylvester rows and of the Paley Hadamard design), 60 kicks each, followed by the fixed-`S`
  actual-prime realisation. The first ILS run crashed when heavy columns needed primes beyond a fixed
  list of 3000. The realisation now generates primes on demand: the smallest unused prime `= 1 mod 4`
  above `max(X, e^(w_j - 0.02))`.
* `W_sa` (`r2/best_ub.py`): equal-weight SA from *structured* starts (`k = floor(r/(M-1))` stacked
  copies of a randomly row-permuted Paley design plus random columns). It is run for
  `tau/w_0 in {1, 1.25, 1.5, 1.75, 2, 2.5}`, from `r = r_P(tau)` upwards in steps of 2 and then
  refined downwards, followed by LP polish of the weights. Then ILS (20 kicks) on
  `D + 0.05 [D = 0]` (Remark 4.2), and the actual-prime realisation of the result.
* rwl1 `W_sc`: the earlier lean pipeline (`q_lean.py`), where available.
* `W_unif = 2(M-1) D_max`: uniform margins at real weights.

**Proposition 4.1 (equal-weight lower bound; exact).** Let a profile have `r` columns, all of the same
weight `tau`, with margins `mu >= D`. Put `T_xy(r, tau) := ceil(r/2 + D_xy/tau)`. Then:
1. `sum_(x<y) T_xy(r, tau) <= r floor(M/2) ceil(M/2)`;
2. hence `W = r tau >= tau r_P(tau)`, where `r_P(tau)` is the least `r` satisfying item 1;
3. `r_P` is nonincreasing in `tau`. So for `w_0 = a_0 < a_1 < ... < a_n`,
   `inf_(tau >= w_0) tau r_P(tau) >= min( min_i a_i r_P(a_(i+1)), a_n r_P(inf) )`, where `r_P(inf)` is
   the least `r` that passes item 1 with `T_xy = floor(r/2) + 1` if `D_xy > 0` and `ceil(r/2)`
   otherwise.

*Proof.*
1. `mu_xy = tau (d_xy - r/2) >= D_xy` with `d_xy` an integer, so `d_xy >= T_xy`. A column with `c`
   entries `+1` separates `c(M-c) <= floor(M/2) ceil(M/2)` pairs, so `sum_(x<y) d_xy` is at most
   `r floor(M/2) ceil(M/2)`.
2. Item 2 is immediate from item 1.
3. `T_xy(r, .)` is nonincreasing in `tau`, so the set of `r` passing item 1 grows with `tau`. On
   `[a_i, a_(i+1)]` this gives `tau r_P(tau) >= a_i r_P(a_(i+1))`. For `tau >= a_n`, `T_xy(r, tau)` is
   at least the `tau = inf` threshold. QED

`r2/eqlb.py` uses the grid `a_(i+1) = 1.01 a_i` up to `4 D_max + w_0`. The thresholds are evaluated
with `D/tau - 10^(-9)` before the ceiling, which can only lower the bound. `LB_eq` can fall below
`2 lambda_max` (minimum 0.95 times it, at `M = 8`): the Plotkin count weights all pairs equally,
whereas the Perron bound weights them by `u_x u_y`. The optimising `tau` is always close to the
floor (median `tau/w_0 = 1.00-1.08`, `eqlb_summary.txt`).

**Remark 4.2 (ties are not realisable by distinct primes; exact).** Let `S` be a sign matrix and
suppose some nonzero `c >= 0` on the pairs with `D_xy = 0` makes `sum c_xy mu_xy(w)` vanish
identically in `w`. For example, two such pairs may have complementary separating sets. Then no
assignment of distinct primes to the columns of `S` satisfies the pair inequalities.

*Proof.* The pair inequalities force `mu_xy >= 0` on these pairs, so every `mu_xy` with `c_xy > 0`
is `0`. That is, `prod_sep p = prod_same p` over the columns. The two products are over disjoint,
nonempty sets of distinct primes, which contradicts unique factorisation. QED

Equal-weight solutions often contain such ties: equal weights make `mu = 0` on every pair at distance
exactly `r/2`. The prime-floor relaxation `W_sc` allows them; the realisation step does not. So
`r2/best_ub.py` refines with ILS on `D + 0.05 [D = 0]` before assigning primes.

**Calibration of the heuristics: local search does not find structured codes.**
* *Equal-weight simulated annealing* (`r2/sacode.c`, `r2/eqw.py`): it searches `S` with
  `d_xy >= T_xy(r, tau)` for `r = r_P(tau), r_P(tau) + 1, ...`, then optimises the weights by LP.
  * At `M = 8` it is within 0-15% of the six proved optima of §4.1 (`r2/test8.py`).
  * At `M = 12` (`X = 24`, `t = 0, 1, 2`) it gives 132.7, 117.9 and 113.4, against ILS 126.9, 118.2
    and 123.0. It needs `r = (1.45-1.56) r_P`.
  * At `M = 20` (`X = 40`, `t = 0`) it finds **no** code up to `r_P + 30`, for any
    `tau/w_0 in {1, ..., 3}` (`r2/test20.py`).
* *Control on known-feasible targets* (`r2/test_unif.py` -> `test_unif_out.txt`). The target is
  `d_xy >= kM/2` for all pairs with `r = k(M-1)`. It is Plotkin-tight and is met exactly by `k`
  stacked copies of the Paley Hadamard design. Best violation over three SA runs of `2 * 10^7` moves:

  | M | 8 | 12 | 20 | 24 |
  |---|---|---|---|---|
  | `k = 2` | 0, 0, 0 | 0, 4, 4 | 40, 42, 46 | 94, 82, 78 |
  | `k = 3` | 0, 0, 0 | 4, 6, 4 | 48, 48, 50 | 92, 90, 88 |

  So local search finds tight codes at `M = 8`, sometimes at `M = 12`, and never at `M = 20, 24`,
  even though such codes exist. Near-Perron discrete solutions must be near-tight codes (Lemma
  6.3(3) and the design-side count of §6.2). So the upper bounds below say nothing about the
  optimum of (Q) at `M >= 20`, and little at `M = 12`.
* *Structured pools* (`r2/hpool.py`): semi-continuous MILP over the columns of `K` random row
  permutations of the Paley design (`M = 20`, `X = 40`, `t = 0`). `K = 10` (190 columns, 300 s)
  gives `W = 439`; `K = 30` (570 columns, 400 s) gives `W = 580`. Both are worse than ILS (363), and
  the MILP's own bound (292) is a bound for that pool only.

**Table 4.2a** (`r2/table42.py` -> `table42_out.txt`). `W_ils`, `W_sa`, rwl1: prime-floor upper bounds from ILS, structured SA and the lean pipeline. "best": the least of them, with method and column count. `W_prime`: the cheapest verified actual-prime profile, with `stair(r)` for its `r`. "char4 min": its least support-`<= 4` character margin.

| M | X | t | D_max | 2 lam | W_cut | LB_eq | W_ils | W_sa | rwl1 | best (method, r) | best/2lam | W_prime (r) | stair(r) | W_unif | char4 min |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 12 | 24 | 0 | 15.8 | 94.2 | 95.3 | 106.2 | 126.9 | 132.3 | 142.4 | 126.9 (ils, 30) | 1.35 | 171.8 (30) | 149.7 | 348.2 | 8.7 |
| 12 | 24 | 1 | 13.1 | 74.5 | 76.3 | 87.5 | 118.2 | 112.7 | 131.4 | 112.7 (sa, 28) | 1.51 | 132.5 (22) | 104.0 | 287.8 | 9.8 |
| 12 | 24 | 2 | 10.5 | 69.9 | 70.4 | 84.9 | 123.0 | 111.6 | 124.4 | 111.6 (sa, 26) | 1.60 | 138.5 (26) | 126.5 | 232.0 | 10.3 |
| 12 | 24 | 3 | 12.2 | 94.5 | 98.1 | 107.8 | 130.2 | 132.5 | 156.5 | 130.2 (ils, 33) | 1.38 | 167.1 (32) | 161.5 | 269.3 | 18.9 |
| 12 | 48 | 0 | 19.2 | 134.3 | 136.8 | 150.6 | 184.2 | 183.4 | 182.5 | 182.5 (rwl1, 34) | 1.36 | 202.9 (33) | 174.7 | 422.3 | 16.3 |
| 12 | 48 | 1 | 23.5 | 117.1 | 121.2 | 125.2 | 160.2 | 166.8 | 182.8 | 160.2 (ils, 30) | 1.37 | 179.3 (30) | 156.7 | 517.6 | 8.9 |
| 12 | 48 | 2 | 13.1 | 103.8 | 109.1 | 124.0 | 154.1 | 148.3 | 164.5 | 148.3 (sa, 34) | 1.43 | 187.6 (33) | 174.7 | 289.1 | 19.8 |
| 12 | 48 | 3 | 23.6 | 121.2 | 124.9 | 141.9 | 165.5 | 168.4 | 181.4 | 165.5 (ils, 30) | 1.37 | 184.9 (30) | 156.7 | 519.4 | 21.5 |
| 20 | 40 | 0 | 21.5 | 201.0 | 201.0 | 231.9 | 362.9 | 321.4 | 593.4 | 321.4 (sa, 74) | 1.60 | 466.7 (74) | 439.2 | 816.5 | 61.9 |
| 20 | 40 | 1 | 24.4 | 204.5 | 207.4 | 244.4 | 362.8 | 329.2 | 475.1 | 329.2 (sa, 74) | 1.61 | 465.5 (74) | 439.2 | 926.9 | 61.1 |
| 20 | 40 | 2 | 20.4 | 203.0 | 204.6 | 239.6 | 363.8 | 316.5 | -- | 316.5 (sa, 82) | 1.56 | 510.5 (82) | 494.7 | 775.5 | 55.5 |
| 20 | 40 | 3 | 25.7 | 209.3 | 209.3 | 234.9 | 367.4 | 331.5 | -- | 331.5 (sa, 78) | 1.58 | 482.4 (78) | 466.9 | 976.6 | 62.0 |
| 20 | 80 | 0 | 27.4 | 256.4 | 259.5 | 304.9 | 451.1 | 428.6 | -- | 428.6 (sa, 92) | 1.67 | 599.7 (92) | 577.6 | 1039.9 | 83.0 |
| 20 | 80 | 1 | 29.4 | 253.0 | 255.0 | 311.4 | 454.7 | 409.1 | -- | 409.1 (sa, 88) | 1.62 | 575.8 (88) | 549.1 | 1116.6 | 87.0 |
| 20 | 80 | 2 | 27.8 | 255.6 | 258.9 | 308.3 | 449.5 | 430.8 | -- | 430.8 (sa, 92) | 1.69 | 596.2 (92) | 577.6 | 1057.8 | 83.9 |
| 20 | 80 | 3 | 29.9 | 262.0 | 262.0 | 306.5 | 452.8 | 422.5 | -- | 422.5 (sa, 90) | 1.61 | 582.4 (90) | 563.3 | 1135.2 | 61.7 |
| 24 | 48 | 0 | 25.5 | 256.1 | 256.1 | 301.7 | 595.0 | 466.5 | -- | 466.5 (sa, 106) | 1.82 | 696.4 (106) | 669.2 | 1171.2 | 108.5 |
| 24 | 48 | 1 | 33.1 | 255.4 | 255.4 | 309.7 | 565.1 | 464.1 | -- | 464.1 (sa, 112) | 1.82 | 683.2 (67) | 394.5 | 1521.2 | 88.9 |
| 24 | 96 | 0 | 36.5 | 329.6 | 330.4 | 390.4 | 742.9 | 560.7 | -- | 560.7 (sa, 112) | 1.70 | 745.2 (112) | 725.8 | 1679.3 | 118.7 |
| 24 | 96 | 1 | 50.0 | 303.5 | 303.5 | 366.0 | 674.6 | 524.8 | -- | 524.8 (sa, 104) | 1.73 | 681.9 (104) | 667.0 | 2299.5 | 100.2 |

**Table 4.2b** (ratios to `2 lambda_max`, median [min, max]; the `M = 8` rows are §4.1, where "best" is the MILP value, optimal for 7 of 12):

| M | X | n | W_cut | LB_eq | best W_sc UB | W_prime | W_unif | W_unif / best UB | best r / M |
|---|---|---|---|---|---|---|---|---|---|
| 8 | 16 | 6 | 1.09 [1.05, 1.21] | 1.20 [1.09, 1.25] | 1.32 [1.24, 1.52] | 1.74 [1.52, 2.52] | 3.17 [2.54, 4.60] | 2.42 [2.05, 3.02] | 1.94 [1.50, 2.25] |
| 8 | 32 | 6 | 1.08 [1.03, 1.31] | 1.10 [0.97, 1.18] | 1.22 [1.11, 1.31] | 1.73 [1.46, 2.16] | 2.91 [2.51, 4.60] | 2.40 [2.21, 3.51] | 2.12 [1.62, 2.50] |
| 12 | 24 | 4 | 1.02 [1.01, 1.04] | 1.16 [1.13, 1.21] | 1.45 [1.35, 1.60] | 1.80 [1.77, 1.98] | 3.51 [2.85, 3.86] | 2.32 [2.07, 2.74] | 2.42 [2.17, 2.75] |
| 12 | 48 | 4 | 1.03 [1.02, 1.05] | 1.15 [1.07, 1.19] | 1.37 [1.36, 1.43] | 1.53 [1.51, 1.81] | 3.71 [2.78, 4.42] | 2.73 [1.95, 3.23] | 2.67 [2.50, 2.83] |
| 20 | 40 | 4 | 1.00 [1.00, 1.01] | 1.17 [1.12, 1.19] | 1.59 [1.56, 1.61] | 2.31 [2.28, 2.52] | 4.30 [3.82, 4.67] | 2.68 [2.45, 2.95] | 3.80 [3.70, 4.10] |
| 20 | 80 | 4 | 1.01 [1.00, 1.01] | 1.20 [1.17, 1.23] | 1.64 [1.61, 1.69] | 2.30 [2.22, 2.34] | 4.24 [4.06, 4.41] | 2.57 [2.43, 2.73] | 4.55 [4.40, 4.60] |
| 24 | 48 | 2 | 1.00 [1.00, 1.00] | 1.20 [1.18, 1.21] | 1.82 [1.82, 1.82] | 2.70 [2.68, 2.72] | 5.27 [4.57, 5.96] | 2.89 [2.51, 3.28] | 4.54 [4.42, 4.67] |
| 24 | 96 | 2 | 1.00 [1.00, 1.00] | 1.20 [1.18, 1.21] | 1.72 [1.70, 1.73] | 2.25 [2.25, 2.26] | 6.34 [5.10, 7.58] | 3.69 [3.00, 4.38] | 4.50 [4.33, 4.67] |

**Reading Tables 4.2a and 4.2b** (evidence, except where marked exact or certified).
1. *Real weights.* `W_cut/2 lambda_max = 1.00-1.05` for `M >= 12` (certified).
2. *Equal-weight granularity.* `LB_eq/2 lambda_max` has median 1.10-1.20 per row of Table 4.2b, and
   1.12-1.22 per `(M, X)` over all 120 saved instances with `M = 8..32` (`eqlb_summary.txt`). That is a
   flat factor (exact bound for equal weights).
3. *Best prime-floor profiles found.* Median `best/2 lambda_max`, at `X = 2M` / `X = 4M`:

   | M | 8 | 12 | 20 | 24 |
   |---|---|---|---|---|
   | median `best/2 lambda_max` | 1.32 / 1.22 | 1.45 / 1.37 | 1.59 / 1.64 | 1.82 / 1.72 |

   At `M = 8`, 7 of the 12 values are proved optimal. The profiles use `1.5M-4.7M` columns.
   * The ratio rises with `M`. But the structured SA lowered every ILS value: by 4-13% at `M = 20`
     (`W_ils/2 lambda_max` from 1.73-1.81 to 1.56-1.69), and by 18-25% at `M = 24` (from 2.21-2.32
     to 1.70-1.82). So the earlier bounds were far from optimal, and nothing indicates that the new
     ones are.
   * The calibration above shows that at `M >= 20` the solvers miss even the known Plotkin-tight
     codes. The rise is therefore no evidence of growth.
4. *Uniform against adapted.* Median `W_unif/best`:

   | M | 8 | 12 | 20 | 24 |
   |---|---|---|---|---|
   | median `W_unif/best` | 2.4 / 2.4 | 2.3 / 2.7 | 2.7 / 2.6 | 2.9 / 3.7 |

   Meanwhile `W_unif/2 lambda_max` grows from about 3 to 5.3-6.3. So the best adapted profiles found
   grow more slowly than uniform margins, but in this range not by much. How much depends on the
   solver: the `M = 24` ratio was 2.3 / 2.8 with ILS alone. The calibration says the solvers are not
   reliable at `M >= 20`, so the comparison of growth rates stays open in both directions.
5. *Actual primes.* `W_prime/2 lambda_max` has median 1.7 (`M = 8`), 1.5-1.8 (`M = 12`), 2.3
   (`M = 20`) and 2.3-2.7 (`M = 24`).
   * At `M = 20`, `W_prime` is within 3.2-6.2% of `stair(r)`, the least weight of any `r` distinct primes
     `= 1 mod 4` above `X`. So the cost is the column count times the log of the `r`-th such prime.
     That is a small-`X` effect: for `r = O(M)` it is `r (log X + O(loglog M))`. At `M = 24` it is
     within 2.2-4.1% for the three SA-based profiles. The one ILS-based profile carries a few heavy
     columns and exceeds `stair(r)` by 73%.
   * Compared with the genuine Paley fakes of §2.4 (median `W_req` = 367, 854, 1215 at
     `M = 12, 20, 24`, `X = 2M`), the best adapted actual-prime profiles cost about 2.4, 1.8
     and 1.8 times less.
6. *Characters* (`r2/char_all.py`). All 40 realised profiles at `M = 12, 20, 24`, from both
   pipelines, pass:
   * the exact pair check;
   * the exact support-`<= 4` character check.

   The minimum of `s_c` over all zero-sum `c` with `|c_x| <= 12` (MILP, floating point) lies in
   `[0.02, 0.33]` and is attained at a pair (support 2) in every case. So the character condition "+char"
   never binds beyond the pairs, as at `M = 8`.

## 5. Self-consistency: the demands of the adapted profile itself

In (Q) the demand matrix `D` is fixed first, and the profile is adapted to it. In P3-gen this order is
not available. `D` is a function of the profile and of the primes: `kappa_x = <a_x, ell>` with
genuine `ell`. This section proves that adaptivity is useless against residues that are random
relative to the profile. The result is a theorem about a probability model, not about actual
circles.

**The uniform residue model `U_X`.** Fix:
* a squarefree profile: distinct primes `p_j > X`, `e = 1`, distinct rows `a_x in {0,1}^r`,
  nonconstant columns;
* arbitrary real weights `w_j`, hence margins `mu_xy` and `W = sum w_j`;
* arbitrary windings `k in Z^M`.

For each odd prime `q <= X` put `A_q = max{a : q^a <= X}` and `m_q* = m_(q^(A_q))`. Let `ell^(q)` be
independent and uniform on `(Z/m_q*)^r`. Pair `(x,y)` collides (unit 1) at level `a` iff
`<a_x - a_y, ell^(q)> = iota_q (k_x - k_y) mod m_(q^a)`; this is the class map of `sharp.md` Cor. 2.2
with point classes `kappa_x = <a_x, ell> - iota k_x`. Let `D_xy = 2 sum_q log q * #{colliding levels}`.

Let `Fbar_X(t) := P(D^ind > t)`. Here `D^ind = 2 sum_q log q * #{a <= A_q : U_q = 0 mod m_(q^a)}`,
with `U_q` independent and uniform on `Z/m_q*`; this is the independent model of §2.

**Theorem 5.1 (random residues defeat every fixed profile).** In the model `U_X`:

```text
P( mu_xy >= D_xy - c_C  for all x != y )  <=  min_x  1 / sum_(y != x) Fbar_X(mu_xy + c_C),
```

and in particular

```text
P( all pair inequalities hold )  <=  2 / ( (M-1) Fbar_X( W/(M-1) + 2 c_C ) ).          (5.1)
```

*Proof.*
1. *Each pair has the independent law.* `u = a_x - a_y in {-1,0,1}^r` is nonzero, so it has a
   coordinate `+-1`. Then `<u, ell^(q)>` is uniform on `Z/m_q*`, and so is its shift by
   `iota(k_x - k_y)`. Reduction mod `m_(q^a)` gives level-`a` collision probability `1/m_(q^a)`,
   nested in `a`. Independence over `q` gives `D_xy ~ D^ind`.
2. *Star pairs are pairwise independent.* Let `y != z`, both different from `x`, and put
   `u = a_x - a_y`, `u' = a_x - a_z`.
   * Choose `j` with `a_yj != a_zj`. Then `a_xj` equals one of them, so `(u_j, u'_j)` is `(0, +-1)`
     or `(+-1, 0)`.
   * In the first case pick `k` with `u_k != 0` (it exists since `x != y`). Then
     `u_j u'_k - u_k u'_j = -u_k u'_j = +-1`. The second case is symmetric.
   * So `ell -> (<u,ell>, <u',ell>)` restricted to coordinates `j, k` is an invertible `2 x 2`
     matrix over `Z/m`. Hence the homomorphism `(Z/m)^r -> (Z/m)^2` is onto, and the image of the
     uniform `ell` is uniform on `(Z/m)^2`.
   * So the collision patterns of `(x,y)` and `(x,z)` at `q` are independent. Over different `q` all
     variables are independent. Hence `D_xy` and `D_xz` are independent.
3. *Second moment.* Fix `x`, let `V_y = {D_xy > mu_xy + c_C}` and `N_x = sum_y 1[V_y]`. The
   indicators are pairwise independent, so `Var N_x = sum Var 1[V_y] <= E N_x`. Chebyshev gives
   `P(N_x = 0) <= Var N_x/(E N_x)^2 <= 1/E N_x`, and `E N_x = sum_y Fbar_X(mu_xy + c_C)` by step 1.
   Feasibility implies `N_x = 0`.
4. *(5.1).* By `p3lower.md` Lemma 1.4(4), `sum_x sum_(y != x) mu_xy = 2 sum_(x<y) mu_xy <= M W/2`.
   So some `x_0` has `sum_y mu_(x_0 y) <= W/2`.
   * If some `mu_(x_0 y) + c_C < 0`, that pair fails surely, since `D >= 0`.
   * Otherwise the `M-1` nonnegative numbers `mu_(x_0 y) + c_C` sum to at most
     `W/2 + (M-1)c_C`. By Markov, at least `(M-1)/2` of them are `<= W/(M-1) + 2c_C`.
   * `Fbar` is nonincreasing, so `E N_(x_0) >= ((M-1)/2) Fbar_X(W/(M-1) + 2c_C)`. QED

**Corollary 5.2 (asymptotics).** Fix `eps in (0,1)` and `C > 0`. There is `M_0` such that for all
`M >= M_0` and `X >= M` the following holds. Every squarefree profile with distinct rows, all
`p_j > X`, any weights, any number of columns, and

```text
W <= (2 - eps) M log^2 M / loglog M
```

satisfies, in `U_X`, `P(all unit-1 pair inequalities hold) <= 4 M^(-eps/16)`.

*Proof.* Put `L = log M`, `l = log L`, `H = L^2/l` and `eta = eps/8`. Let `c` count the level-1
collisions at the window primes `q in (M^(1-eta), M]`. All of these are `<= X` and do not divide `N`.
1. *The threshold.* `t := W/(M-1) + 2c_C <= (2 - eps/2) H` for `M >= M_0`.
2. *Reduction to `c`.* Window collisions alone give `D^ind >= 2(1-eta) L c`. So
   `Fbar_X(t) >= P(c >= k)` with `k = floor(t/(2(1-eta)L)) + 1`.
3. *Size of `k`.* `k <= ((1 - eps/4)/(1 - eps/8)) L/l + 1 <= (1 - eps/8) L/l + 1`, using
   `(1 - eps/8)^2 >= 1 - eps/4`.
4. *Lower bound on `P(c >= k)`.* `c` is a sum of independent `Bernoulli(1/m_q)`. So

   ```text
   P(c = k) >= e_k(p) prod_q (1 - p_q) >= ((lambda_Q - k p_max)^k / k!) e^(-lambda_Q - 1)
   ```

   (`p3lower.md` Lemma C.1). Here `lambda_Q = sum 1/m_q = log(1/(1-eta)) + o(1)`, which lies in
   `[eta/2, 1]` (Mertens), and `k p_max -> 0`.
5. *Conclusion.* `log P(c >= k) >= -k log k - O(k) - 2 >= -(1 - eps/8)L - O(L log l/l)`
   `>= -(1 - eps/16) L`. Now (5.1) gives `<= 2 M^(1 - eps/16)/(M-1) <= 4 M^(-eps/16)`. QED

**What this means.**
* *(Adaptivity.)* The bound holds for **every** fixed profile, adapted or not, with any weights and
  any number of columns. Against residues that are random relative to the profile, non-uniform
  margins buy nothing beyond the constant. The heavy collisions land on uniformly random pairs, so
  every star must carry margin at the `1/M`-quantile of the demand law. For comparison, uniform
  margins at the `1/M^2`-quantile succeed with probability `>= 1/2` (union bound), at cost
  `2(M-1) q_(1/M^2)` with `q_(1/M^2) ~ 4H` (Poisson heuristic, not proved here). So in this model
  `W_needed / (M H)` lies between `2 - o(1)` (proved) and about `8` (heuristic); adaptivity can save
  at most a constant factor.
* *(What it does not say.)* It is not a statement about actual circles and not a lower bound for
  P3-gen. A fake builder may **select** primes. With twin columns, changing the twin prime of row `x`
  re-randomises `kappa_x` alone, and a Moser-Tardos resampling over the pair events needs only
  `Fbar(mu_xy) <~ 1/(2eM)`, i.e. uniform margins at the `1/M`-quantile. That is again
  `Theta(M log^2 M/loglog M)` under the Poisson heuristic, with a smaller constant. Getting below
  this needs residues correlated with the margins, i.e. prime selection with control of residues
  modulo many `q ~ M` at once. No such mechanism with `W = O(M log M)` is known, and no argument
  excludes one.
* *(Relation to `p3lower.md`.)* Cor. A.2 there says that pair-cut **certificates** cannot use heavy
  pairs against adapted margins. Theorem 5.1 is the dual statement for **constructions**: adapted
  margins are only useful if the residues follow them. Theorem C there (uniform margins plus `H`)
  and Theorem 5.1 (random residues, any margins) bracket the same phenomenon from two sides.
  Neither applies to actual circles without a hypothesis on how their residues relate to their
  profile.

### 5.1 Numerical checks (evidence)

**`Fbar_X`.** It is computed by dynamic programming on a log grid of step `10^(-3)`, once with
every term rounded down and once rounded up. The rounded-down table gives a rigorous lower bound
`Fbar_lo <= Fbar`, and every bound below uses it, so the bounds are valid. The two tables agree to 4
digits at `t = 2 log X` (`theorem_check_out.txt`).

**Monte Carlo in `U_X`** (2000 draws per row; `theorem_check.py`, `theorem_check2.py`). Profiles:
* *uniform*: the `M-1` nonconstant columns of the normalised Paley Hadamard matrix, equal weights;
* *paley-b8*: `construct` with equal weights;
* *adapted*: the realised prime profiles of §4, fitted to one genuine `D`, with their weights scaled.

`P_hat` is the observed frequency with which all pair inequalities hold, and "bound" is
`min_x 1/sum_y Fbar_lo(mu_xy)` from Theorem 5.1.

| M, X | W | uniform `P_hat` (bound) | paley-b8 `P_hat` (bound) |
|---|---|---|---|
| 12, 24 | 59.6 | 0 (0.16) | 0 (0.15) |
| 12, 24 | 238.6 | 0.031 (1) | 0.010 (0.68) |
| 12, 24 | 477.1 | 0.98 (1) | 0.88 (1) |
| 24, 48 | 152.5 | 0 (0.078) | 0 (0.069) |
| 24, 48 | 610.2 | 0 (0.61) | 0 (0.27) |
| 24, 48 | 1220.4 | 0.78 (1) | 0.51 (1) |

**Adapted against uniform at equal `W`** (`M = 8`, `X = 16`, three adapted profiles; full output in
`theorem_check2_M8_out.txt`):

| `W / W_0` (`W_0` = adapted profile's own `W`, 69-74) | 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| adapted `P_hat` (three profiles) | 0, 0, 0 | 0, 0.002, 0 | 0.001, 0.009, 0.005 | 0.001, 0.018, 0.006 |
| uniform `P_hat` | 0.002-0.005 | 0.51-0.64 | 0.97-0.98 | 1.00 |
| paley-b8 `P_hat` | 0-0.001 | 0.25-0.33 | 0.84-0.90 | 0.98-0.99 |

So against residues that are random relative to the profile:
* every observed frequency lies below the Theorem 5.1 bound (which is far from tight);
* uniform margins beat the Paley profile, which beats the adapted ones;
* an adapted profile stays below 2% feasibility even at four times its weight, because its
  low-margin pairs are hit at random.

**Actual primes, the adapted profiles' own demands.** For each realised adapted prime profile of §4
at `M = 8`:
* take its own Gaussian primes `pi_j` (norms `p_j` from 17 to a few hundred) and random conjugates;
* compute its genuine unit-1 demands exactly at the same `X`;
* check its own pair inequalities (400 conjugation choices per profile).

Result (16 profiles, with and without separation): on average **5.6 to 12.5 violated pairs**, and
**0 of 400** conjugation choices feasible in every case. So the profiles are pair-feasible for the
`D` they were fitted to, and infeasible for their own residues. Uniform margins at the same `W` are
infeasible too, since `W` is below the random-residue threshold. The point of the previous table
is that adaptation makes things worse, not better.

## 6. Does adaptivity with few columns remove the log^2 M/loglog M effect?

The question has two orders of quantifiers, and the answers differ.
* **(Q)-order:** `D` is fixed first, then a profile is fitted to it. §6.1 shows that separation is
  lower order and that characters are the only extra constraint; §6.2 collects the evidence.
* **P3-gen order:** the profile and its primes come first, and `D` is whatever their residues give.
  §5 answers this for residues that are random relative to the profile; §6.3 puts the two together.

### 6.1 Exact reductions

**Lemma 6.1 (concatenation).** Let `P1`, `P2` be squarefree profiles on the same `M` points, with
disjoint sets of primes, rows `a^1_x in {0,1}^(r1)`, `a^2_x in {0,1}^(r2)` and weights `w^1`, `w^2`.
Let `P = P1 (+) P2` have rows `a_x = (a^1_x, a^2_x)` and weights `(w^1, w^2)`. Then for all
`x != y` and every integer `c` with `sum c = 0`:

```text
W(P) = W(P1) + W(P2),    mu_xy(P) = mu_xy(P1) + mu_xy(P2),    s_c(P) = s_c(P1) + s_c(P2).
```

*Proof.* `d_xy = sum_j w_j |a_xj - a_yj|` and `w(v(c)) = sum_j w_j |<c, a_(.j)>|` are sums over
columns, and the columns of `P` are those of `P1` followed by those of `P2`. Subtract `W/2`. QED

So `W_cut`, `W_sc` and `W_prime` (with disjoint prime pools) are subadditive in the demand matrix.
The same holds for the character margins, which are what "+char" constrains.

**Proposition 6.2 (separation is lower order; class (Q), and every member of `A_X` for the lower
bound).** Let `C = sqrt 2`, `M >= 3`, `D >= 0` symmetric with zero diagonal, and fix any order of
the points with positions `t_x = i(x)/(M-1)`. Put `E_xy := 2 log((M-1)/|i(x) - i(y)|)`, the separation
demand of these positions (Lemma 1.1). Then:
1. `lambda_max(E) <= max_x sum_y E_xy <= 2(1 + log 2)(M-1) - 4`;
2. *(real weights)* `W_cut(D) <= W_cut(D + E) <= W_cut(D) + 2 pi (1 + log 2)(M-1)`;
3. *(prime floor)* if `w_0 <= 4 log(M-1) + 2`, then
   `W_sc(D + E + 1) <= W_sc(D) + (2M-2)(4 log(M-1) + 2)`;
4. a profile with margins `mu >= D + E + 1` (off the diagonal) satisfies the full unit-1 pair system
   with the placement `delta = Delta t`, i.e. it is feasible for "`+sep`".

So the separation requirement changes the real-weight optimum by `O(M)` and the prime-floor optimum
by `O(M log M)`. (Prop. 3.2 shows that a term of order `M` is unavoidable in the averaged bound.) It
cannot produce `log^2 M/loglog M` growth.

*Proof.*
1. `E >= 0`, so `lambda_max(E)` is at most the largest row sum. Put `n = M-1` and
   `F(m) = sum_(k=1)^m log(n/k) = m log n - log m!`. From `log m! >= m log m - m + 1` (compare with
   `int_1^m log t dt`) we get `F(m) <= m log(n/m) + m - 1` for `m >= 1`, and `F(0) = 0`. The point at
   position `i` has row sum `2(F(i) + F(n-i))`. If `1 <= i <= n-1`, this is
   `<= 2(n H(i/n) + n - 2) <= 2(n log 2 + n - 2)`, with `H` the natural-log entropy. If `i in {0, n}`,
   it is `<= 2(n - 1)`.
2. The left inequality holds because `E >= 0`. For the right one, measures add (Lemma 6.1 at real
   weights): if `lambda` is optimal for `D` and `lambda'` for `E`, then `lambda + lambda'` is
   feasible for `D + E`. Finally `W_cut(E) <= pi lambda_max(E)` (`p3lower.md` Thm A(3)) and item 1.
3. Let `n = 2^ceil(log2 M) <= 2M - 2`, take the Sylvester Hadamard matrix `H_n`, restrict it to `M`
   of its rows, and drop the columns that are constant on them. Any two distinct rows of `H_n` agree
   in exactly `n/2` columns, and a dropped column separates no pair. So each pair is separated by
   exactly `n/2` of the `r' <= n - 1` remaining columns. Give every column the weight
   `w := 4 log(M-1) + 2 >= w_0`. Every margin is then `w n/2 - w r'/2 >= w/2 = max(E + 1)`, at total
   weight `r' w <= (2M-2) w`. Add this to an optimal prime-floor solution for `D`; the sum is still
   prime-floor feasible and has margins `>= D + E + 1`.
4. By Lemma 1.1(2) it suffices that `|t_x - t_y|^4 (Pi_sep - m^2) >= m^4 Pi_same`, where
   `m = m_xy = e^(D_xy/2)`. Since `log Pi_sep = W/2 + mu` and `log Pi_same = W/2 - mu`, the condition
   reads `2 mu + log(1 - m^2/Pi_sep) >= 2 D + 2 E`. Now `mu >= D + E + 1` gives
   `log Pi_sep = W/2 + mu >= D + 1 > D + log 2`, so `m^2/Pi_sep < 1/2`. Then
   `2 mu + log(1 - m^2/Pi_sep) >= 2D + 2E + 2 - log 2`, which suffices. QED

*Check* (`prop62_check.py` -> `prop62_check_out.txt`; ALL PROP 6.2 CHECKS PASSED):
* the row-sum bound for `M = 3..400` (equality at `M = 3`);
* the Sylvester separation count in exact integers for `M = 3..130`;
* `W_cut(E)` by the full LP for `M = 6, 8, 10, 12`: `20.9, 31.8, 42.9, 54.3`, against
  `2 lambda_max(E) = 18.9, 29.9, 41.2, 52.6`. So at real weights separation alone costs about
  `4.5 (M-1)`;
* Lemma 6.1 on 200 random profiles with 20 random characters each.

*For actual primes* item 3 needs `r' <= 2M-2` distinct primes `= 1 mod 4` whose logarithms lie in a
window of width `O(1/M)` around `w`. That is a short-interval prime count, which is not proved here.
The uniform-margin profiles of `construct.md` §2 (genuine primes, pair margins `>= 2 tau`,
`W = O(M log M)`) can be used instead. So the conclusion "separation is `O(M log M)`" holds for
actual primes as well, with the constant of that construction.

**Example 6.2A (one collision layer is cheap at the prime floor; exact).** Let the `M` points be
partitioned into `m` classes, all of size `k+1`, where `k+1 = 2` or `k+1` is the order of a Hadamard
matrix. Let `D = h (B - I)`, with `B` the class indicator (`B_xy = 1` iff `x, y` are in the same
class), so `lambda_max(D) = hk`. Then for every `w_0 > 0` there is a profile with all weights
`>= w_0`, margins `mu >= D`, and

```text
W  <=  2 lambda_max(D) + 2 M w_0.
```

*Proof.*
1. Let `n_A = 2^ceil(log2 m) < 2m`, with Sylvester rows `h_1, ..., h_m in {+-1}^(n_A)` (pairwise
   orthogonal).
2. Let `sigma_0, ..., sigma_k in {+-1}^k` have pairwise inner products `-1`, with every coordinate
   balanced over the `k+1` vectors. For `k = 1` take `(1), (-1)`. Otherwise take the rows of a
   normalised Hadamard matrix of order `k+1` with the constant column removed.
3. Give point `i` of class `c` the row `h_c (x) sigma_i in {+-1}^(n_A k)`. Its inner product with the
   row of point `i'` of class `c'` is `(h_c . h_c')(sigma_i . sigma_i')`. That is `-n_A` when `c = c'`,
   `i != i'`, and `0` when `c != c'`.
4. Every column `(a, b)` takes both signs inside each class, since `sigma`'s columns are balanced. So
   every column is nonconstant.
5. With equal weights `tau`, `mu_xy = -tau (s_x . s_y)/2`. That is `tau n_A/2` inside classes and `0`
   across them. Take `tau = max(w_0, 2h/n_A)`. Then
   `W = tau n_A k = max(n_A k w_0, 2hk) <= 2 m k w_0 + 2 lambda_max <= 2 M w_0 + 2 lambda_max`. QED

*What this explains.* At level 1 the genuine demand is a sum of layers
`D = sum_(q <= X) 2 log q (B_q - I)`, with `B_q` the partition into classes of `kappa(q)`. Each layer
alone costs the Perron bound plus `O(M w_0)`, at least when its classes are equal-sized. For such
layers Lemma 6.1 then gives `W <= 2 sum_q lambda_max(D^(q)) + 2 M w_0 #{q}`. The second term is `~ 2 M X` for
`X ~ M`, i.e. of order `M^2`. So at the prime floor the whole difficulty of (Q) is to merge the
`~ X/log X` layers into `O(M)` shared columns. At real weights this merging is free (Lemma 6.1 with
no floor). It also shows that a lower bound for (Q) at the prime floor cannot hold for all `D`:
structured demands such as one layer are cheap.

**Characters do not reduce the same way.** By Lemma 6.1 a character margin can be raised by
concatenation. The `b = 8` Paley profile of `construct.md` Theorem 1 has `s_c >= 2 tau` on pairs and
`s_c >= tau n` on the other characters, with `n = ||c||_1`. But a fitted profile `P1` can have
`s_c(P1)` as low as `-W(P1)/2`, when `v(c) = 0`. Repairing small `n` this way would cost a factor
`~ M` in `W`. So "+char" is a genuine constraint on the fitted profile itself, and it is checked
directly in §4: exactly for supports `<= 4`, and by MILP up to `|c_x| <= 12` at `M = 8`.
*Heuristic* (not used anywhere): for a column `a` uniform on `{0,1}^M` and a `+-1` character of
support 4, `E(|<c,a>| - 1/2) = 1/4`. Pairs give `0`, and by Szarek's inequality pairs are the extreme
case for every zero-sum `c`. So pairs are expected to bind before characters, as observed.

### 6.2 The real-weight optimum is the Perron bound, and what discreteness must preserve

**Lemma 6.3 (Perron attainment; class: any symmetric `D >= 0` with zero diagonal).** Let `D` be
irreducible, `Lambda = lambda_max(D)`, with unit Perron vector `u > 0`. For any `lambda >= 0` on cuts
with `mu(lambda) >= D`, put `G := sum_s lambda_s s s^T`. Then

```text
W(lambda) - 2 Lambda  =  sum_s lambda_s <s, u>^2  +  sum_(x<y) 4 u_x u_y (mu_xy - D_xy),          (6.1)
```

and both sums are `>= 0`. Consequently:
1. `W_cut(D) = 2 Lambda` iff `2(Lambda I - D)` lies in the cone generated by the matrices `s s^T`,
   `s in {+-1}^M` (a constant `s` cannot occur in such a representation, since `<1, u> > 0`);
2. in that case every optimal `lambda` is supported on `u`-balanced cuts (`<s, u> = 0`) and has
   `mu = D` exactly;
3. the same identity holds for every profile (Lemma 1.4): a profile with `W <= (1 + eps) 2 Lambda`
   has `sum_j w_j <s_j, u>^2 <= 2 eps Lambda` and u-weighted total slack `<= 2 eps Lambda`.

*Proof.* `G` has diagonal `W` and off-diagonal entries `-2 mu_xy` (`p3lower.md` §1.3). So
`u^T G u = W - 2 sum_(x != y) u_x u_y mu_xy`, and `u^T G u = sum_s lambda_s <s,u>^2`. Write
`mu = D + (mu - D)` and use `sum_(x != y) u_x u_y D_xy = Lambda`. This gives (6.1). Both sums in
(6.1) are nonnegative, since `mu >= D` and `u > 0`.
1. If `W(lambda) = 2 Lambda`, both sums vanish, so `mu = D` and `G = W I - 2D = 2(Lambda I - D)`.
   Conversely, if `2(Lambda I - D) = sum lambda_s s s^T`, the diagonal gives `W = 2 Lambda` and the
   off-diagonal gives `mu = D`.
2. and 3. are read off from (6.1). QED

So, at real weights, the Perron bound is attained exactly when the shifted demand matrix is a
nonnegative combination of `u`-balanced cut matrices; the margins are then tight on every pair. For a
generic `u` no cut is exactly `u`-balanced, so exact attainment is exceptional. The relevant
quantity is the excess `W_cut - 2 Lambda`, which by (6.1) splits into an imbalance part and a slack
part.
Item 3 says what a discrete profile must preserve to stay within `1 + eps` of the bound: nearly
`u`-balanced columns and nearly tight margins on average.

*Counting does not decide (Q)* (heuristic; not used anywhere). There are two first-moment counts, and
they point in opposite directions.
* *Demand side.* With `r = cM` columns there are `2^(cM^2)` sign matrices. Distinguishing the
  `binom(M,2)` demands to relative precision `eps` needs about `(M^2/2) log2(1/eps)` bits. So the
  number of designs is no obstruction to `W = O(M log M)` with `O(M)` columns.
* *Design side.* By `p3lower.md` Lemma 1.4(4), `sum_(x<y) mu_xy <= M W/4`. So `W <= K M log M`, with
  `sum D ~ M^2 log X`, leaves an average slack `mu - D` of `O(K) log M` per pair. That is `O(1)`
  columns of weight `>= w_0`, so at least half of all pairs have their distance in a window of `O(1)`
  columns. A *random-like* `S` has pair distances spread over `~ sqrt r` columns. If these behaved
  independently (the local-CLT heuristic of de Launey-Levin for partial Hadamard matrices), such `S`
  would number `2^(Mr) (c/sqrt r)^(M^2/4)`. That is `< 1` once `r < (M/8) log2 r (1 + o(1))`. Since
  `r <= W/w_0 <= K M` (as `w_0 >= log X >= log M`), this happens roughly when `M > 2^(8K)`: beyond any
  computation, but it is growth if it is real.
* The design-side count is no evidence either. The same count says that Hadamard matrices of order
  `n` should not exist (`2^(n^2) (c/sqrt n)^(n^2/2) -> 0`). Yet they do for every `M` used here
  (Paley). So it only says that near-optimal designs for large `M`, if they exist, are structured
  and not typical. §4.2 shows that local search does not find such structured designs even when they
  are known to exist.

**Evidence that genuine demands nearly attain the Perron bound** (`perron_check.py` -> `perron_check_out.txt`;
LP values with certified dual bounds; per instance the identity (6.1) is checked to `1e-6`):

| M | X | instances | median `2 lambda_max` | `W_cut/2 lambda_max`: median [min, max] | excess `W_cut - 2 lambda_max`: median [min, max] | largest slack term | max identity error | max relative gap LP vs certified bound |
|---|---|---|---|---|---|---|---|---|
| 8 | 2M | 6 | 40.5 | 1.085 [1.052, 1.212] | 3.46 [2.49, 5.81] | 1.02 | 3e-14 | 3e-14 |
| 8 | 4M | 6 | 64.5 | 1.082 [1.034, 1.313] | 6.08 [1.95, 14.98] | 2.13 | 4e-14 | 8e-14 |
| 12 | 2M | 6 | 81.7 | 1.028 [1.007, 1.054] | 2.21 [0.52, 4.40] | 0.00 | 8e-14 | 7e-11 |
| 12 | 4M | 6 | 117.2 | 1.028 [1.016, 1.051] | 3.20 [1.91, 5.26] | 0.00 | 7e-14 | 2e-10 |
| 20 | 2M | 4 | 203.8 | 1.004 [1.0000, 1.014] | 0.81 [0.00, 2.84] | 0.00 | 1e-13 | 1e-10 |
| 20 | 4M | 4 | 256.0 | 1.010 [1.0000, 1.013] | 2.47 [0.00, 3.33] | 0.00 | 3e-13 | 2e-10 |
| 24 | 2M | 2 | 255.7 | 1.0000 [1.0000, 1.0000] | 0.00 [0.00, 0.00] | 0.00 | 2e-13 | 1e-10 |

Reading (evidence; the LP values are certified to the stated relative gap, the excess is a
difference of certified quantities):
* From `M = 12` on, the slack term of (6.1) is `0` to solver precision. The optimal margins are
  **exactly tight**, `mu = D`, on every pair, and the whole excess is cut imbalance
  `sum lambda_s <s,u>^2`.
* The excess stays `O(1)` in absolute terms (at most 15, typically 0-6), while `2 lambda_max` grows
  like `M log M`. At `M = 20, 24`, 4 of the 10 instances have excess `< 0.01`.
* So for genuine demands the real-weight pair channel costs `2 lambda_max(D) + O(1)` in this range,
  far below the Goemans-Williamson cap `pi lambda_max(D)` of `p3lower.md` Thm A(3). Equivalently,
  `2(lambda_max I - D)` is, up to `O(1)`, a nonnegative combination of nearly `u`-balanced cut
  matrices. Nothing here is proved for general `M`.

**Corollary 6.4 (granularity; equal weights).** Let a profile have `r` columns, all of weight `tau`,
and margins `mu >= D`. Put `phi_xy := ceil(D_xy/tau + r/2) - (D_xy/tau + r/2) in [0, 1)`. Then

```text
W  >=  2 lambda_max(D)  +  tau sum_(x<y) 4 u_x u_y phi_xy.
```

*Proof.* `mu_xy = tau (k_xy - r/2)` with `k_xy` the number of separating columns, an integer. So
`mu_xy >= D_xy` forces `mu_xy - D_xy >= tau phi_xy`. Insert this into (6.1). QED

If the fractional parts `phi` are spread over `[0,1)`, the extra term is about
`tau (1/2) sum 4 u_x u_y ~ (M-1) tau >= (M-1) w_0`. So with equal weights, discreteness costs an
additive `~ M log X`, of the same order as `W_cut ~ 4 M log X`, i.e. a constant factor and not
growth. With unequal weights the lattice is finer, and the M = 8 optima show that the loss can be
much smaller.

### 6.3 Putting the two orders together

*Class of each statement:* §6.1, Lemma 6.3 and Prop. 4.1 are proved for every symmetric `D >= 0`.
Prop. 4.1 is proved for equal-weight profiles only. Theorem 5.1 is proved in the model `U_X`. The
(Q) tables are solver evidence on genuine `D`.

1. **(Q)-order: absorbed at real weights; undecided at the prime floor.** For genuine `D` with
   `M <= 24` (`M <= 32` for Prop. 4.1):
   * The real-weight optimum is the Perron bound up to a factor that tends to 1 (Lemma 6.3 table).
     Uniform margins cost `(M-1) D_max/lambda_max` times that bound, which grows from 3.6 to 5.6 over
     the range of §2. So at real weights the `log^2 M/loglog M` effect is absorbed by adaptivity
     (proved, `p3lower.md` Thm A).
   * The granularity of equal weights costs a flat factor (exact, Prop. 4.1: `LB_eq/2 lambda_max`
     has median 1.12-1.22 for every `M` in 8..32).
   * Separation costs `O(M log M)` at most (Prop. 6.2).
   * The remaining question is whether near-tight `+-1` designs adapted to `D` exist with `O(M)`
     columns. The best profiles found cost, as medians, 1.22-1.32 times the Perron bound at `M = 8`
     (7 of 12 optimal), 1.37-1.64 at `M = 12, 20`, and 1.72-1.82 at `M = 24` (§4.2, Table 4.2b).
     Their ratio to the Perron bound rises with `M`, and uniform margins cost 2.3-3.7 times them. But the solvers demonstrably miss structured codes that exist
     (§4.2 calibration), so this trend is **no evidence** of growth. *Caveat:* at these `M`,
     `D_max/w_0` is only 2-6, and no finite computation of this size decides the asymptotic (Q).
2. **P3-gen order: the effect returns against independent residues.** Theorem 5.1 and Cor. 5.2 hold
   for every fixed profile, however adapted, with any weights and any number of columns. They say
   that `W <= (2 - eps) M log^2 M/loglog M` succeeds with probability `<= 4 M^(-eps/16)` when the
   residues are uniform relative to the profile. At `M = 8` every fitted profile is infeasible for its
   own genuine residues (0 of 400 conjugation choices, §5.1).
3. **Where the effect lives.** The extra `log M/loglog M` is:
   * not a pair-LP effect (Thm A);
   * not a separation effect (Prop. 6.2);
   * not an equal-weight granularity effect at the level of the Plotkin count (Prop. 4.1).

   Two places remain.
   * *Coupling* (proved in a model): the effect appears when the residues are independent of the
     profile (Theorem 5.1 here; Theorem C of `p3lower.md` for uniform margins under `H`).
   * *Discreteness* (open): whether adapted near-tight `+-1` designs with `O(M)` prime columns exist
     for genuine `D`. An earlier version of this draft said discreteness was excluded "on small-`M`
     evidence". That is withdrawn: only `M = 8` is solved to optimality, and the `M >= 20`
     heuristics are not reliable (§4.2).
4. **What this means for a growth proof.** A proof of `M = o(log R/loglog R)` through pair data must
   show that **actual** clusters cannot couple their residues to their exponent profile. The precise
   form would be a statement about how the residues `pi_j mod q`, for `q <= X`, of the Gaussian
   primes dividing an actual cluster are distributed relative to the cluster's own exponent matrix.
   * That is a Hypothesis-H-type statement (`p3lower.md` §4), now needed for adapted profiles and not
     only for uniform ones.
   * Conversely, a fake builder who wants `W = O(M log M)` must select primes whose residues modulo
     the many `q ~ M` follow a prescribed adapted pattern.
   * Neither is decided here.

*Remark (selection is not excluded by these methods).* Let the builder choose each `p_j` from `K`
candidates after seeing their residues. A union bound over the `K^r` choices then needs a failure
bound `<= K^(-r)`.
* Theorem 5.1 gives only `M^(-eps/16)`.
* Inside one star the events are fully independent when the differences `a_x - a_y`, `y != x`, are
  independent mod every `m_q`. Then `P(N_x = 0) = prod_y (1 - Fbar(mu_xy))`. If every margin exceeds
  the median of `D^ind`, this is `>= 2^(-(M-1))`.
* For `C = sqrt 2` margins are `>= 0`. The rows of `S diag(w)^(1/2)`, with `S` the `+-1` form of the
  profile, then have pairwise nonpositive inner products, so `r >= M/2`. Hence no star bound beats
  `K^(-r) <= K^(-M/2)` once `K >= 4`.
* A bound of the form `exp(-c sum_(x<y) Fbar(mu_xy))` over all pairs would suffice for polynomial
  `K` near the threshold of Cor. 5.2. But pairs are not independent: collisions of `(x,y)` and
  `(y,z)` at `q` force a collision of `(x,z)`. No such bound is proved.

So the random-residue argument does not extend to P3-gen with prime selection without a new input.

## 7. Files (`round7/p3data_checks/`)

All scripts run from this directory with `python3`. The repository and the round-5/6 folders were
not modified.

| file | content | type |
|---|---|---|
| `gauss7.py`, `demands7.py` | exact Gaussian arithmetic, `T_Q` discrete logs, Paley profiles, collision moduli (round-7 partial work, unchanged) | library |
| `verify_direct.py` (+`_out.txt`) | collision computation vs direct Gaussian arithmetic mod `Q` and in exact `Z[i]` | exact; ALL DIRECT CHECKS PASSED (re-run) |
| `genuine_stats.py` -> `genuine_trials.jsonl`, `data/` | part 1: 4000 trials, exact moduli; `D`/moduli matrices for `M <= 32` | exact per trial |
| `genuine_table.py` -> `genuine_table_out.txt`, `part1_compact.py` -> `part1_compact_out.txt` | aggregation; independent-model simulation; tails; Poisson-binomial comparison | evidence |
| `paley_selfcheck.py` -> `paley_selfcheck_out.txt`, `.json` | Paley fakes with their own primes; `tau_req`, `W_req` | evidence |
| `maxcut.c` (`libmaxcut.so`), `qcore.py` | exact top-`K` max-cut by Gray code (checked against brute force for `M <= 14`); cut utilities | library |
| `qsolve.py` | LP (full, column generation with exact or heuristic pricing, certified dual bound), semi-continuous MILP, LNS, prime realisation | library |
| `verify_profile.py` | exact checks: pairs, Lemma 1.1(2) placement, characters (support `<= 4` exact; MILP minimum) | exact (character MILP: floating point) |
| `q_run.py` -> `q_results_M8.jsonl`, `q_run_M8_log.txt`, `q_run_M8b_log.txt`, `q_progress_M8.txt` | part 2 full pipeline at `M = 8` (patched by `round7/patch_qrun.py`, `patch_qrun2.py`; `round_nearest` window bug fixed afterwards) | LP certified; MILP optima where marked; upper bounds verified exactly |
| `q_run_M12_log.txt`, `superseded/` | the failed `M = 12` runs of `q_run.py` (see §3) | record only |
| `rwl1.py` | reweighted-`l^1` support selection, restricted semi-continuous MILP, nearest-prime menu MILP | library (upper bounds only) |
| `q_lean.py` -> `q_lean_M{12,20,24}.jsonl`, `q_lean_M*_log.txt` | §4.2: lean (Q) pipeline at `M = 12, 20, 24` | LP certified; prime profiles verified exactly; `W_sc` heuristic upper bounds |
| `fls.py`, `ils_run.py` -> `ils_M{8,12,20}.jsonl`, `ils_M*_log.txt` | flip local search / ILS library (`fls.py`, used by `r2/`); `ils_run.py` is the first ILS run, which crashed in the prime realisation (fixed list of 3000 primes) and is superseded by `r2/ils2_run.py` | library / record only |
| `sc_full_M12.py` -> `sc_full_M12_out.txt` | full semi-continuous MILP at `M = 12` (2047 cuts, 1500 s): no lower bound above `W_cut` was obtained | record |
| `q_lean_table.py` | tables of §4.2 | -- |
| `q_run_large.py` | heuristic-pricing variant for `M >= 32` (not run for this write-up) | library |
| `q_table.py` -> `q_table_out.txt` | the tables of §4.1 | -- |
| `perron_check.py` -> `perron_check_out.txt`, `.json` | Lemma 6.3 on genuine `D`: `W_cut/2 lambda_max`, imbalance and slack terms of (6.1) | LP certified; identity (6.1) checked |
| `prop62_check.py` -> `prop62_check_out.txt` | Prop. 6.2 and Lemma 6.1 checks | exact (Sylvester counts) / float (row sums, LP) |
| `theorem_check.py` -> `theorem_check_out.txt`, `.json` | §5.1: `Fbar` bracket, Monte Carlo vs Theorem 5.1, actual-prime self-consistency | evidence (bracket rigorous) |

**Finishing session, `p3data_checks/r2/`** (run from `r2/` with `python3`; `notes.md` is the working log):

| file | content | type |
|---|---|---|
| `eqlb.py` -> `eqlb_out.txt`, `eqlb_out.json`, `eqlb_summary.txt` | Prop. 4.1: equal-weight lower bound `LB_eq` for all 120 saved instances (`M = 8..32`), medians per `(M, X)` | exact (float only in `D/tau`, rounded down) |
| `sacode.c` (`libsacode.so`), `eqw.py` | equal-weight code search by simulated annealing; Plotkin `r_P`; LP polish | library (upper bounds) |
| `test8.py`, `test20.py` | SA calibration on genuine `D` at `M = 8` (vs proved optima), `M = 12, 20` | evidence |
| `test_unif.py` -> `test_unif_out.txt` | SA on Plotkin-tight targets met by stacked Paley designs, `M = 8, 12, 20, 24` | evidence (calibration) |
| `test_struct.py` | SA from stacked-Paley starts, `M = 20` | evidence |
| `hpool.py` | semi-continuous MILP over permuted-Paley column pools, `M = 20` | upper bounds |
| `ils2_run.py` -> `ils2_M{8,12,20,24}.jsonl`, `_log.txt` (`M = 8`: calibration against the proved optima of §4.1) | ILS from Sylvester/Paley starts; on-demand prime realisation (Remark 4.2 handling); exact pair check; support-4 characters | upper bounds; prime profiles verified exactly |
| `best_ub.py` -> `best_ub_M{12,20,24}.jsonl`, `_log.txt` | structured SA over a `tau` grid, LP polish, ILS on `D + 0.05[D=0]`, prime realisation | upper bounds; prime profiles verified exactly |
| `wcut_extra.py` -> `wcut_extra.json`, `_out.txt` | certified `W_cut` for `M = 24`, `X = 96` | LP certified |
| `char_all.py` -> `char_all_out.txt` (= `char_all_ils2_out.txt` + `char_all_best_1220_out.txt` + `char_all_best_24_out.txt`), `char_all_*.json` | exact pair re-check, exact support-`<= 4` characters, MILP minimum of `s_c` over `|c_x| <= 12` for every realised profile | exact / floating point (MILP) |
| `table42.py` -> `table42_out.txt`, `table42_rows.json` | Tables 4.2a, 4.2b | -- |
