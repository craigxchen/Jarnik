# Referee report on Part X.3 (final): "What minimising the right inverse can buy" (block: inverse)

Referee for Astra's Part X, section X.3 (answer.md lines 2406 to 2445). I read section X.3 line by line and checked
it against README_IX.md sections 1.3, 1.4, 1.8 and 5. I wrote my own integer-programming code; Astra's scripts were
not imported. Scripts are in `.../round31/ref10_inverse/` (`minV.py`, `minV2.py`).

## 0. Verdict

* **X.3, the theorem (finite-lift dichotomy plus the mod-four bound): CORRECT WITH CORRECTIONS.** Every step is
  right. The corrections concern scope, not mathematics. (i) The bound applies only to the **owner rows**
  `i in S` (the `m` marks), not to all `q - 1` labels. (ii) The "owner design" (pair targets exactly
  `u gamma_i/gamma_{i'}`) is IX.3's choice and is not forced by the Legendre profile; I give the robust form
  below (Remark R). It covers every target that differs from the design by even parent-prime valuations, and the
  cost of that difference is charged to the target. (iii) In the exact case "almost every" means: all but
  `O(eta^{-2})` owners. That is a bounded number, which is stronger than `o(m)`.
* **Computed (referee, exact ILP over all integral right inverses).** Take the full positive class (`m = (q-1)/2`)
  and unit weights. Then for q = 11, 19, 23 the minimum of `sum_{i in S} sum_k |V_ik|` over **all** integer `V`
  with `D V = (H - D)/2` is exactly the mod-four parity bound, `m(q-3)/2`. Every owner row has weight exactly
  `(q-3)/2` at the optimum. For q = 31 and 43 the exact minimum (exhaustive, with non-owner entries unbounded) is
  strictly above the bound: 240 > 210 and 504 > 420. Astra's Cauchy–Schwarz form of the bound is weaker but has
  the same leading term `(1/2) m T`.
* **What is closed.** The README_IX §5 item 1 suggestion "minimise the witness heights over right inverses `U`
  (closest-vector problem)" cannot produce small or near-genuine witnesses. Every exact global witness family for
  the IX.3 design, for every `U`, has `log Norm(twist_i) >= (1/2 - eta) T` for all but `O(eta^{-2})` owners. Hence
  `log J_i = log(p_i H_i^2) >= (1 - 2 eta) T`, roughly `q log P`. Against the judge's single-`V`
  value `log J ~ 1.2 q log P`, optimising over `V` gains at most a factor of about 1.2. The goal "forced bound `J_j <= p_j^{1+o(1)}`" is
  therefore unreachable on this design.
* **What is not closed.** (a) Whether some large-height exact witness passes e-all at `X = 3M^A`. Height alone is
  not an e-all failure, and X.3 contains no e-all test. (b) Local (finite-cutoff, non-global) witnesses on the
  Legendre profile. X.3's first alternative merely says that they have `sum log H_k >= (1/2)L_X - O(q log P)`,
  which every finite-cutoff lift is allowed to have. (c) Other designs or profiles. **P3-all\* remains open.**
  Astra says this herself in her interpretation paragraph, and that paragraph is correct.

## 1. The statement, exactly (with quantifiers)

**Setting (fixed first).** A prime `q = 3 mod 4` and `chi` the Legendre symbol mod q. The label set is
`I = F_q^*` (column 0 deleted), `K_ik = chi(k - i)` and `a_ik = 1 + K_ik`. The mark set `S` is a subset of the
positive class `{i : K_i0 = 1}` and `m = |S|`. The distinct split primes are `p_k in [P, 2P]` (`k in I`) with chosen
factors `pi_k`, and `gamma_k = pi_k/conj(pi_k)`. Put `T = sum_k log p_k`, `B = sum_k (log p_k)^2` and
`P_max = max p_k`. Fix the units `u_{ii'}` (allowed by the windings). Next fix the cutoff `X`. Finally take
**one** family `r_k in Q(i)^*` (`k in I`), chosen after all of the above.

**Hypothesis (Loc_X).** For every odd prime `l <= X` with `l` not in `{p_k}`, all `r_k` are `l`-units,
`r_k conj(r_k) = p_k`, and for all `i, i' in S`
`prod_k (r_k/conj r_k)^{a_ik - a_i'k} = u_{ii'} gamma_i/gamma_{i'}  (mod l Z[i])`.

**Definitions.** `r_k/pi_k = u_k G_k/conj(G_k)` with `u_k` a unit and `G_k in Z[i]` conjugate-primitive (no
rational integer factor > 1; `(1+i)` absorbed into `u_k`). `H_k = |G_k|^2`. Every other Hilbert-90 representative
of `r_k/pi_k` is `n G_k` times a unit, so it has norm `>= H_k`. `L_X = sum log l` over the same primes `l`.

**Theorem X.3 (Proved; referee-checked).** Under (Loc_X) at least one of the following holds.
1. `2 sum_{k in I} log H_k >= L_X - (q+3)/2 log(2P) - log 2`. This holds as soon as one pair identity is false in
   `Q(i)`. In fact the sharper bound with `sum_k |d_k| log H_k` in place of `2 sum log H_k` holds.
2. All pair identities hold exactly in `Q(i)^*`. Then `V_jk := v_{pi_k}(r_j/pi_j)` satisfies `D(I + 2V) = H`, and
   `sum_{i in S} log H_i >= (1/2)(mT - sqrt(m q B) - m log P_max)`.
   For `0 < eta < 1/2` with `2 eta T > log P_max`, all but at most `qB/(2 eta T - log P_max)^2` owners `i in S`
   satisfy `log H_i >= (1/2 - eta)T`.

**Asymptotic form (Proved).** Let `m ≍ q` (i.e. `q << m <= q`), `P = q^{O(1)}`, `P -> infinity`. Then `T ~ q log P` and
`qB/T^2 = O(1)`, so case 2 gives `log H_i >= (1/2 - eta)T` for all but `O(eta^{-2})` owners. If moreover
`X/(q log P) -> infinity`, then `L_X >= theta(X) - T - log 2 = (1 + o(1))X` and
`sum_k log H_k >> min(X, q^2 log P)`.

**What the cutoff does (the "dichotomy").** Either the minimal lift's total twist is at least about `X/2`, which is
the size every CRT lift may have (finite-cutoff lift proposition: height `e^{O(X)}`), or the identities are
global. In the global case the parent-prime valuations obey the integral relation `D(I + 2V) = H` and the
mod-four bound applies. So small twists force global identities as soon as
`L_X > sum_k |d_k| log H_k + (q+3)/2 log 2P + log 2`, i.e. once `X` exceeds, roughly, the total log-height. That
threshold is linear in the log-heights.

## 2. Proof check

**(a) Length of `d`. Proved.** Take `d = a_i - a_{i'}` with `i, i' in S`. Over `F_q`,
`sum_k chi(k-i) chi(k-i') = -1`. There are `q - 2` nonzero products, so the coordinates where the two rows
disagree number `(q-1)/2`, each with `|d_k| = 2`. The coordinates `k = i` and `k = i'` each have `|d_k| = 1`. At
`k = 0`, `d_0 = 0` because both marks lie in the positive class. Hence `||d||_1 = q + 1`, `||d||_inf <= 2`, and
`d = h mod 2` with `h = e_i - e_{i'}`.

**(b) Local to global. Proved.** We have
`r_k/conj r_k = u_k^2 gamma_k (G_k/conj G_k)^2`. Put
`B_d = prod_{d_k > 0} (pi_k G_k^2)^{d_k} prod_{d_k < 0} (conj pi_k conj G_k^2)^{|d_k|}`. Then
`prod (r_k/conj r_k)^{d_k} = unit * B_d/conj B_d`. Put `A = B_d conj(pi_i) pi_{i'}`. The pair identity is
equivalent to `A = u conj(A)` for one unit `u`. The element `A` is an `l`-unit at every exterior `l`, so (Loc_X)
gives `l | A - u conj A` in `Z[i]` for every exterior odd `l <= X`. This includes both factors of a split `l`,
because the residue ring is `Z[i]/l`. If the identity fails, then `A - u conj A != 0`, so
`exp(2 L_X) <= Norm(A - u conj A) <= 4 Norm(A)` with
`Norm(A) = p_i p_{i'} prod_k p_k^{|d_k|} H_k^{2|d_k|} <= (2P)^{q+3} prod_k H_k^{2|d_k|}`.
This is exactly case 1. The answer to question (2) is yes: this is the dichotomy, and the height bounds are right.
It is the standard "the difference of numerators is a nonzero integral element of bounded norm" argument.
Nothing is discarded: the Gaussian gcd plays no role because divisibility is by the rational prime `l`.

**(c) Exponent relation. Proved, for general `r_k`.** The witnesses are **not** assumed to lie in the parent's
S-unit group. The `G_k` may contain any primes, exterior ones included. Only the valuations at the parent primes
`pi_k` are used. For `j != k`, `v_{pi_k}(r_j) - v_{conj pi_k}(r_j) = 2V_jk`. For `j = k` the difference is
`1 + 2V_kk`. So the `pi_k`-valuation of `r_j/conj r_j` is `(I + 2V)_jk`, and taking the `pi_k`-valuation of the
exact identity gives `d(I + 2V) = h`. Exterior prime factors and units drop out of this equation.

**(d) "Modulo four". Proved.** `d(I + 2V) = h` is equivalent to `dV = (h - d)/2`. Since `d = h mod 2`, we have
`dV = hV = V_i. - V_i'. mod 2`. Hence `V_ik - V_i'k = (h_k - d_k)/2 mod 2`. This is the content of reducing the
identity mod 4. It fixes, column by column, the parities of all owner entries `V_ik` (`i in S`) up to one free
bit `c_k`. For `i != k`, `(-1)^{V_ik} = sigma_k K_ik`. The diagonal entry is free. Non-owner rows `j not in S` are
not constrained mod 2. All parity patterns allowed by this are realised by integral `V`: `D` has all Smith
invariants 1 (I recomputed this for q = 11, 19, 23), so `ker_Z(D) mod 2 = ker_{F_2}(D mod 2)`. The parity bound is
therefore the exact answer to the mod-2 relaxation.

**(e) Paley variance. Proved.** `sum_{k : V_ik odd} log p_k = (1/2)(T - (K f)_i - eps_i)` with
`f_k = sigma_k log p_k`, `f_0 = 0` and `|eps_i| <= log P_max`. `K` is antisymmetric (`chi(-1) = -1`), so
`K^T K = qI - 11^T` and `sum_{i in S} (Kf)_i^2 <= ||Kf||^2 <= q||f||^2 = qB`. Cauchy–Schwarz gives the sum bound,
and Chebyshev gives the count bound. Since `log H_i >= sum_k |V_ik| log p_k >= sum_{V_ik odd} log p_k`, both
claims of case 2 follow. The role of the Paley structure is that no single choice of column signs `sigma` can
make a large fraction of owner rows nearly even. The mod-four relation forces about half of each owner row to be
odd.

**(f) Remark R (referee's strengthening; Proved).** Suppose the pair target is any element whose parent-prime
valuation vector is some `h'`, for example `gamma_i/gamma_{i0} prod_k gamma_k^{2E_ik}`. Then `h' = d mod 2` is
forced, and the same argument gives
`sum_k (|V_ik| + |E_ik|) log p_k >= sum over a set of k forced odd`, with the same Paley bound. So moving weight
from the twist into the target does not help. The target excess `prod p_k^{2|E_ik|}` enters the P3-char collision
cost `L_0` of IX.3, so it is not free either. The extreme case `V = 0`, `E = -B`, i.e. `r_k = pi_k` and target
`prod gamma^d`, is the genuine witness. It pays the full amount in the target. *Computed:* the joint minimum of
twist plus target excess equals `m(q-3)/2` at q = 11 and 19 (§3).

**(g) Constant check.** The displayed case-1 inequality matches what (b) proves. The case-2 sum bound and the
exceptional count are as stated. The claim "`sum_k log H_k >> min(X, q^2 log P)`" is right under the stated
ranges.

## 3. Computations (referee's own code in `ref10_inverse/`: `minV.py`, `minV2.py`, `minV6.py`; exact integer checks)

**Method.** The columns of `V` decouple: `D V_{.k} = ((H - D)/2)_{.k}`. For each column `k` I computed the exact
minimum of `sum_{i in S} |V_ik|` over **all** integer columns, as follows (`minV6.py`). Non-owner entries are left
completely unbounded and eliminated by the Smith normal form `P D_N Q = diag(d)` of the non-owner block `D_N`. An
owner vector `x` is feasible iff `(P(b - D_S x))_i = 0 mod d_i` for every `i`, and `= 0` past the rank. All owner
vectors with the forced parities are enumerated in order of increasing `l^1` weight, so the first feasible weight
is the exact minimum. The **parity bound** is `sum_k min_c #{i in S : c + B_ik odd}`, the exact optimum of the
mod-two relaxation of §2(d). The **CS bound** is Astra's `(1/2)(mT - sqrt(mqB) - m log P_max)` with unit weights
(comparable parent, `log p = 1`). For q <= 23 an independent MILP (`minV.py`, HiGHS, with `D(I + 2V) = H` verified
exactly on the solution) gives the same minima.

Setup: full positive class (`m = (q-1)/2`), `i0 = min S`.

| q | m | Smith inv. of D (all 1?) | inv. of D_N | exact min of total owner weight | parity bound `m(q-3)/2` | CS bound | mean per owner row | `(q-1)/2` |
|---|---|---|---|---|---|---|---|---|
| 11 | 5 | yes | {2} | 20 | 20 | 10.8 | 4 | 5 |
| 19 | 9 | yes | {2,4} | 72 | 72 | 48.8 | 8 | 9 |
| 23 | 11 | yes | {2} | 110 | 110 | 78.2 | 10 | 11 |
| 31 | 15 | (not recomputed) | {2,58} | 240 | 210 | 158.4 | 16 | 15 |
| 43 | 21 | (not recomputed) | {2,508} | 504 | 420 | 333.1 | 24 | 21 |

**Computed.**
1. At every tested q the exact minimum is at least the parity bound, and the parity bound is at least the CS bound.
   This confirms X.3 numerically.
2. For q = 11, 19 and 23 the bound is attained exactly. At the optimum every owner row has weight `(q-3)/2` and
   every owner entry is 0 or +-1.
3. For q = 31 and 43 the exact minimum is strictly larger. At q = 31, 15 of the 30 columns exceed the parity
   bound by 2. At q = 43, 21 of the 42 columns exceed it by 4. The obstruction is the odd part of the
   non-owner Smith invariants (`58 = 2·29`, `508 = 4·127`): it is a congruence that no sign pattern of minimal
   weight satisfies. The mean owner-row weight is `(q+1)/2` at q = 31 and `(q+5)/2` at q = 43. So the bound is
   a lower bound and is not sharp in lower-order terms. Whether the true optimum is `(1/2 + o(1))T` or `c T` with
   `c > 1/2` is not decided by these data; only the lower bound matters for the negative answer.
4. Consequence: the optimal twist is at least `log H_i`, roughly `(1/2)T`, i.e. `(q/2) log P`, for a comparable
   parent, and `log J_i = log(p_i H_i^2)` is roughly `q log P` or more. The judge's single-`V` value
   `log J ~ 1.2 q log P` is therefore within a factor of about 1.2 of optimal.
5. Free targets (Remark R; MILP over `v` and the target excess `e` with `D v - e = b`): the joint minimum of twist
   plus target excess is 20, 72 and 110 at q = 11, 19 and 23. This is again `m(q-3)/2`, so trading twist for
   target height buys nothing.
6. Side remark: X.3 says nothing about non-owner rows. At q = 19, minimising all rows with equal weight gives an
   owner total of 90 and a non-owner total of 72. Non-owner entries can be large; one optimum at q = 23 reached
   the box bound 60. These rows matter for e-all but are outside the theorem.

## 4. Consequences, and what is and is not closed

1. **Closed (Proved).** On the IX.3 Legendre design with exact global witnesses, the minimal witnesses have
   `log J_i >= (1 - 2 eta) T` for all but `O(eta^{-2})` owners, uniformly over all right inverses `U`, and
   indeed over all integral `V` and all auxiliary exterior twists. The README_IX §5 item 1 alternative
   "a forced bound `J_j <= p_j^{1+o(1)}`" is therefore false for this design. The closest-vector problem has
   optimum at least `(1/2 - o(1))T` per owner. The computations show it is attained at q <= 23 and exceeded at
   q = 31 and 43.
2. **What a large twist means for e-all.** X.1 says that a twist that is small relative to the cutoff (Norm `<= H`,
   `X >= sqrt2 (p H^2)^{n/2}`) and passes the ray tests must be genuine. X.3 shows that the twist is necessarily
   large: `log H` is roughly `(q/2) log P`. X.1 then needs `X` about `exp(n q log P)`, which is far beyond
   `3M^A`. So X.1 gives **nothing** at polynomial cutoffs for these witnesses. A large twist does not by itself
   fail e-all. The judge's earlier observation (6/10 and 16/18 labels fail the singleton test at `X = 3m` for one
   `V`) is evidence, not proof. I did not test e-all on the optimal witnesses (NOT CHECKED); that is the remaining
   half of the judge's suggestion.
3. **Not closed.** A P3-all\* fake on the Legendre profile with **local** witnesses (one lift per cutoff,
   alternative 1 of the dichotomy). Alternative 1 is merely the generic size of CRT lifts, and X.2 uses local
   witnesses. Also not closed: designs whose targets are not of the form "small height times even parent
   valuations", and other profiles. X.3 bounds `sum log H`. P3-all\* charges `W = sum e_j log p_j`, and no
   hypothesis of the record relates the two. Astra's interpretation paragraph says exactly this and is correct.

**Corrected statement (full).** *Fix `q = 3 mod 4`, the Legendre profile, `S` in the positive class of column 0,
split primes `p_k in [P, 2P]` with factors `pi_k`, units `u_{ii'}`, and then a cutoff `X`. Let `r_k in Q(i)^*` be
any family satisfying (Loc_X). Then either (1) `sum_k |d_k| log H_k >= L_X - (q+3)/2 log 2P - log 2` for some
pair `d`, or (2) the pair identities hold in `Q(i)^*`, the parent-prime valuations satisfy `D(I + 2V) = H`, and the
owner twists satisfy `sum_{i in S} log H_i >= (1/2)(mT - sqrt(mqB) - m log P_max)`, with all but
`qB/(2 eta T - log P_max)^2` owners having `log H_i >= (1/2 - eta)T`. The same bound holds for `|V| + |E|` if the
targets are altered by even parent valuations `2E`. The bound concerns owner rows only. It is attained exactly at
q = 11, 19, 23. At q = 31 and 43 the true minimum exceeds it by 1 and 3 per row on average (computed). It says nothing about e-all passing or failing, nor about P3-all\* with
witnesses that are local at the cutoff.*

**Is README_IX §5 item 1 (second half) answered negatively?** Yes, as far as heights go. Minimising over `U`
cannot give witnesses of height `p^{1+o(1)}`, or anything below `exp((1/2 - o(1))T)` per owner. No, as far as
P3-all\* goes. Whether the height-optimal exact witnesses pass e-all at `3M^A` is untested. The large heights make
a pass unlikely on the evidence of the judge's computations, but no theorem excludes it.
