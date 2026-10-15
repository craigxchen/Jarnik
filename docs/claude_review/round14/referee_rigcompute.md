# Referee report on `round14/rigcompute.md` (Rig(K): the first finite test)

Adversarial referee, round 14. All scripts are in `referee_rigcompute_checks/`, run from that directory. They are
written from the definitions in `sharp.md` Def. 1.1 and `L4K.md` §1. Nothing is imported from `rigcompute_checks/`
or the repository: the author's files are read only as data, and pickles are opened with a class-blocking unpickler.
The labels below are *proved*, *exact* (finite computation with exact integer or rational arithmetic), *solver*
(floating-point MILP or branch-and-bound, independent of the author's code) and *evidence*.

## 0. Verdict

**The computations survive. The class does not match P3-all\*, and that changes which conclusions survive.**

What survives:
* **The arithmetic.** Every number I re-derived agrees exactly with the write-up: the forced sets, the one-modulus
  survivor counts, all 10 joint witnesses, all first-FLEX witnesses, the `Q`, `n` and ratio tables, and the
  Lemma M side conditions.
* **One certificate, re-checked independently.** For `first4`, `K = 4`, `Q_core >= 61`: my own MILP on the same
  level-1 relaxation is infeasible at `Q = 61` and feasible at `Q = 67` (solver).

What does not hold as stated:
* **The class is weaker than P3-all\*.** The write-up's class ("good" means the per-unit budget only) is exactly
  Rig(K) as round 13 posed it. It is **strictly weaker** than what P3-all\* imposes on labels at norm `<= K`.
  * Condition (1.2) of `sharp.md` Def. 1.1 (e-all) requires one angle vector `phi`. For each `v`, the same
    `x = <v, phi>` must serve all four units `s`. So a label system can lie under a P3-all\* member only if the
    slab system `|sin(<v,phi> - s pi/4)| >= nu_s m_(v,s) e^(-w(v)/2)` (all `v`, `s`) has a solution `phi`.
  * A necessary label-only consequence, for each single `v`, is that the four windows are compatible ("per-v angle
    compatibility"). It implies, for example, `m_(v,0)^2 + m_(v,2)^2 <= Norm Z_v` and
    `m_(v,1)^2 + m_(v,3)^2 <= 2 Norm Z_v`.
  * Genuine labels satisfy all of this. The per-unit budget does not imply it.
* **Consequence for certificates and witnesses.**
  * Every certificate (lower bound) transfers to P3-all\*.
  * **No joint witness does.** All 10 fail: 8 fail per-v compatibility (exact). The two `last4` witnesses pass
    per-v compatibility, but their slab systems at norm `<= K` have no solution `phi` (complete branch-and-bound,
    solver).
  * In the class that P3-all\* actually constrains, the "exact on both sides" cores are therefore only lower
    bounds.
  * The loose-prime "stall" and the "C1 seeds" reading are not established.
  * For `first4`, `K = 4`, the bound `q_flex1 = 83` is not an upper bound in that class: once per-v compatibility is
    imposed, `83` is one-modulus **rigid** (exact).
* **What I can prove in the stronger class** (§4.6), where `Q_core^sl` is the core of the class "per-unit +
  slab at norm `<= K`":
  * `61 <= Q_core^sl(first4, 4) <= 103` (new slab-feasible one-modulus witnesses at 107);
  * `Q_core^sl(first4, 5) >= 109`;
  * so the core **does grow from `K = 4` to `K = 5` in this class too**;
  * the `r`-saturation bound `Q_core <= 173` (`first6..8`, `K = 4`) survives (its witnesses are slab-feasible at
    `phi = arg pi_j`).

No growth improvement, uniform bound or cap follows, and the write-up claims none. The barrier logic is unchanged.

| item | severity | §|
|---|---|---|
| I1: the per-unit class is strictly weaker than P3-all\* at the label/angle level; witnesses, "exact both sides", the loose stall and the C1-seed reading hold only for the relaxation | major | 1, 4.5, 4.6 |
| I2: `q_flex1(first4, K)` is not an upper bound for P3-all\*: 83 (`K = 4`) and 163 (`K = 5`) are rigid once per-v angle compatibility is imposed; at `K = 4` the one-modulus segment in that class reaches 103 | major (consequence of I1) | 4.6 |
| I3: "local `<=` joint `<=` one-modulus" is not a logical nesting; Lemma K's "cannot grow with `K`" is proved only for `q <= p_min^((K+1)/2)` | minor | 3 |
| I4: "local rigidity vacuous for the loose sets" is false for `rand4` (`q = 17` is locally rigid); for `first5`, `q = 23` is locally rigid outside the initial segment | minor | 3 |
| I5: the `first8` FLEX witness at 179 has no stored lift (I supplied and verified one) | minor | 4.3 |
| I6: regime wording ("tight = `log p/log r` near 1" for `first4`, where it is up to 2.43); finite `r <= 8` cannot be placed in regime (1.4) | minor | 2 |
| I7: all-level goodness vs. levels `q^a <= X`; the deep certificates need `X >= 27, 49, 121, 169` | minor | 2 |

## 1. Class definition check

I compared each ingredient with `sharp.md` Def. 1.1 (e-all) and Prop. 1.2.

| ingredient | P3-all\* (`sharp.md` Def. 1.1) | write-up | verdict |
|---|---|---|---|
| moduli | odd primes `q <= X`, `q` not dividing `N`; levels `a <= A_q` (`q^a <= X`) | odd `q` not dividing `N`, all levels | matches for `X >= ` every modulus used (I7) |
| labels | `r_j(q)` units of norm `p_j` mod `q^(A_q)`, lower levels reductions; label `log(r_j/conj r_j)` | `ell in (Z/m_(q^a))^r`, lifts compatible | equal: `r_j/conj r_j = (r_0/conj r_0) t^2` with `t` in the norm-one torus, so the achievable labels are exactly one parity coset |
| parity | forced (Prop. 2.1) | `ell_j = [(p_j/q) = -1] mod 2` | verified (`t_sanity.py`: genuine `L_j mod 2` equals this for 8 primes, all `q < 200`) |
| units / collisions | `prod (r_j/conj r_j)^(v_j) = i^s mod q^a` | `<v, ell> = s iota mod m q^(a-1)` | equivalent. The author's checker takes `s` from `d mod m` and counts deeper levels by `n q^b \| d`, which is the compatible choice of `iota` at every level. Genuine collision iff `q \| D_s(Z_v)`: 0 mismatches in 21,888 checks (`t_sanity.py`) |
| monomials | **every** `v != 0` | `1 <= \|\|v\|\|_1 <= K` | a relaxation, as Rig(K) asks |
| condition | (1.2): `\|sin(<v,phi> - s pi/4)\| >= nu_s m_(v,s) e^(-w(v)/2)` for all `v`, `s`, with **one** `phi in R^r` | `m_(v,s) <= nu_s^(-1) e^(w(v)/2)` per `(v,s)` | **strictly weaker** (I1) |

**I1 (major).** The per-unit budget is (1.2) with `|sin| <= 1` applied separately to each `(v, s)`. It discards two
couplings that P3-all\* keeps.

1. **Coupling across `s` for each `v`.** The same `x = <v, phi>` must lie in all four allowed sets. With `y = 2x`,
   the four conditions say that `(cos y, sin y)` lies in the box
   `[2 m_2^2/N - 1, 1 - 2 m_0^2/N] x [m_3^2/N - 1, 1 - m_1^2/N]`, where `N = Norm Z_v` and `m_s = m_(v,s)`. So
   per-v compatibility is a rational condition: the box is non-empty and meets the unit circle. Non-emptiness alone
   already requires `m_0^2 + m_2^2 <= N` and `m_1^2 + m_3^2 <= 2N`. Genuine labels satisfy these, because
   `m_s <= |D_s|`, `X^2 + Y^2 = N` and `(X-Y)^2 + (X+Y)^2 = 2N`. The per-unit budget allows `m_0^2 + m_2^2` up to
   `2N` and `m_1^2 + m_3^2` up to `4N`. (`reflib.angle_ok_exact`: exact `Fraction` arithmetic, 0 disagreements with
   a direct trigonometric test on 20,000 random instances.)
2. **Coupling across `v`.** `x_v = <v, phi>` is linear in `v`. This is the slab system.

Rig(K) as written in round 13 §7.2 uses the per-unit budget. So the write-up answers the question as posed, and
its **certificates** are valid for P3-all\* (the relaxation direction is fine). But:
* the "growing" branch of round 13 (P3-all\* reduces to P3-G, then L7) needs only certificates, so it can use the
  write-up's lower bounds;
* the "stalling" branch (flexible tops seed C1, then the slab conditions are added with angles chosen after the
  labels) needs label systems whose slab system is solvable.

None of the 10 joint witnesses has a solvable slab system (§4.5). So the write-up's step "stall, hence C1 seed"
does not go through for any of them. "Exact on both sides" holds for the per-unit relaxation only.

Everything else in the class definition is right:
* the column-sign symmetry `ell -> eps * ell` preserves budgets, parity and (as I add) the slab system, with
  `phi -> eps * phi`;
* the free-moduli reduction (Lemma R1) stays valid in the stronger class, because demands only decrease.

## 2. Regime and parameters

* **Regime (1.4)** is asymptotic: `p_j <= M^(h_0)`, `r <= rho M`, moduli up to `X = 3 M^A`, with `A` chosen after
  `c_0`. A finite instance (`r = 4..8`, `p_j <= 400`, `K <= 7`) lies "in" it only vacuously, for suitable constants.
  The write-up says so: §7 calls everything about Rig(K) in the regime evidence. I agree with that labelling. **I6
  (minor).**
  * The "tight" sets are not tight in the sense of `log p_j/log r -> 1`: for `first4`, `max log p_j/log r =
    log 29/log 4 = 2.43`; for `last4` it is 4.3.
  * The dichotomy is "smallest available primes" versus "larger primes". The link to near-extremal clusters
    (`L4K.md` §5 item 4) is a heuristic there and stays one here.
* **Moduli.** All of `F_K^ex` is used, plus collision-free labels above it. This is the right universe for norm
  `<= K` (Theorem 1(b) of L4K). It is "moduli `<= Q_K`" in Rig(K) for every `Q_K >= q*(K)`, and it needs
  `X >= q*`.
* **Levels (I7, minor).** The write-up imposes goodness at all levels, while P3-all\* counts only `q^a <= X`.
  * For witnesses this is stronger, hence harmless.
  * Certificates that model depth (`3^3, 5^3, 7^2, 11^2, 13^2`) are valid relaxations of P3-all\* only when `X` is
    at least those prime powers. That holds in the regime (`X = 3M^A`), but should be stated.
* **Parameters of each claim.** The tables state `P`, `r`, `K` and `Q` for every entry. I found no claim whose
  parameters are missing or mis-stated, apart from the class issue I1.

## 3. Lemmas R, K, M

**Lemma R (reductions): correct.** I re-derived all six items.
1. *Free moduli.* A collision at `q^a` implies one at `q`. So a level-1 collision-free label has no collision at
   any level, and replacing a label by it only removes factors.
2. *Level-1 relaxation.* A good system's level-1 demand divides its full demand.
3. *Minimal signatures.* Demand is monotone in each signature. A replacement `(sigma', mask')` with
   `sigma' ⊆ sigma` and `mask' ⊆ mask` keeps `AND mask = 0`. For depth tables, dominance must be by multiset
   (signature with depths), which is what `deep` tables use. I did not re-audit `cpc.c`'s dominance code; my MILP
   (§4.4) uses no dominance at all.
4. *Individual admissibility.* A single factor `q` already exceeds the budget. With depth `b`, the same argument
   applies to `q^b`.
5. *Dropping moduli* removes non-negative demand.
6. *Upper bounds.* Immediate from the definitions.

All six remain valid in the stronger class (per-unit plus slab), because the slab conditions are also monotone in
the demand.

**Lemma M: correct.** Its side condition `F_(K+1)^ex ∩ [3, Q] = F_K^ex ∩ [3, Q]` (both equal to all primes `<= Q`
not dividing `2N`) is verified exactly from my own forced sets for `first4` at `(K, Q) = (4, 61), (5, 109),
(6, 193)` and for `first5` at `(4, 101), (5, 163)`. The consequence "the all-norm core is `>= sup_K Q_core(K)`" needs
exactly this condition, and it holds in these instances.

**Lemma K: the proof is correct but the conclusion overreaches (I3, minor).**
* What is proved: for a fixed `q`, the set of individually admissible labels at norm `K'` is constant for
  `K' >= k(q)`. In particular it is constant once `p_min^(K+1) >= q^2`.
* What is claimed besides: "it cannot grow with `K`", and "it ends at 31 for `first4`" for every `K`. That needs,
  for every `q > 31`, a non-class admissible label at `K = k(q)`, i.e. at a `K` that grows with `q`.
* What is checked: `loc_first4_K7.json` checks `q <= 113`. Since `5^8 >= 113^2`, the conclusion is proved for
  `q <= 113` and every `K`. Above 113 it is an extrapolation.
* A short proof is available, but is not given: for large `q` the `O((log q)^r)` monomials with
  `Norm Z_v < q^2` cannot exclude all non-class labels among about `q^r` candidates.

**Nesting "local `<=` joint core `<=` one-modulus segment" (I3).**
* The second inequality is a theorem: a FLEX witness is a good system.
* The first is **not** a logical implication. Local rigidity puts each label in *some* column-sign class, while
  coherence needs *one* class for all `q`. It is an observed comparison of initial-segment ends.
* As a set statement it fails:
  * `rand4`, `K = 4`: `q = 17` (`n = 4`) is locally rigid (`loc_rand4_K4.json`: LOCe true at `{3, 5, 7, 17}`) but
    lies outside the joint core `{3, 5, 7}`.
  * `first5`, `K = 4`: `q = 23` is locally rigid outside the initial segment `{3, 7, 11}`.
* **I4 (minor).** The sentence "local rigidity ... vacuous for the other sets" is therefore false for `rand4`. The
  `q = 17` rigidity is real but elementary: norm-1 collisions are inadmissible there, so `c_j != 0 mod 4`, which with
  the parity forces `c_j = +-L_j`.

## 4. Independent computations

### 4.1 Forced sets (*exact*; `rfree.c`, `rforced.py`, outputs `rforced_*.json`, `rforced.log`)

My own complete DFS works as follows:
* column 1 runs over representatives of the unit orbits on `Z/n` (`d | n`);
* the other columns are exhaustive, with bitsets of partial sums `<v', c>` by norm;
* every FREE vector is re-verified in Python.

I ran it over every prime `q <= Q_scan = 4A_K + 5`. It agrees **exactly** with `thr_*.json` for
* `first4`, `K = 4, 5, 6`;
* `rand4`, `bal4`, `last4`, `K = 4, 5`;
* `first5`, `K = 4, 5`.

For the remaining sets (`r = 4`, `K = 6, 7`; `first6`, `K = 4, 5`; `first7`, `first8`, `K = 4`; `rfreevec.py`):
* every stored FREE vector is parity-respecting and collision-free for all of `V_K`;
* every prime `q <= Q_scan` is classified as forced, aborted or free.

`Q_scan` itself agrees with Theorem 1(b) of L4K for greedy `B_h` sets (e.g. `A_4 = 1636`, `A_7 = 12649` for
`r = 4`).

### 4.2 Joint witnesses (*exact*; `rverify.py`, output `rverify_all.txt`)

My verifier uses integer-label semantics:
* the stored top lift is the label at every level;
* the collision depth of `v` at `q` is `1 + v_q(<v,top>/n)`, and it would be infinite if `<v,top> = 0` (this never
  happens);
* the unit is read from `<v,top> mod m`.

It checks coverage of **my** forced set, parity, lift consistency, and exact budgets for every `v` with
`||v||_1 <= K` (both signs) and every `s`. It also recomputes the genuine logs and the column-sign masks
independently. **All 10 `jwmix_*.pkl` pass**, and they reproduce the write-up's table exactly:
* coherent on `F ∩ [3, Q_core]`;
* non-coherent at the stated `Q`;
* the counts outside every class: `34/36, 62/64, 23/26, 37/40, 74/75, 33/33, 37/40, 67/67, 34/36, 64/65`.

The upper bounds `Q_core < Q` are therefore correct **for the per-unit class**.

### 4.3 One-modulus data (*exact*; `onemod_angle.py`, `rflex2.py`, `rflex8.py`; outputs `om_*.txt`)

* **Survivor counts** (exhaustive over all parity-respecting labels at `q0`, genuine labels on `F \ {q0}`, level-1
  relaxation; orbit-weighted). For `first4`, `K = 4` they reproduce the round-10 referee table (`referee_L4K.md`
  §6.3, row "genuine on `F_4`") exactly: `47: 8/0`, `83: 36/24`, `107: 216/144`, `113: 336/240`,
  `127: 6592/6080`, `167: 161088/160320`, `199: 500800/500160` (the full scan of `F ∩ [61, 199]` is in
  `om_first4_K4.txt`).
  * The number "90,048 at 167, norm 4" quoted in the task and in round 13 §7.2 is the *mixed* variant (genuine on
    `F_5`). Like-for-like, the norm-4 count is 160,320.
* **The RIGID/FLEX pattern** for `rand4`, `bal4`, `last4` at `K = 4` (`q <= 79`) also matches the write-up's §3.2
  table:
  * segment ends 29/37/29, first FLEX 31/41/31;
  * `rand4` rigid above it at `{37, 43, 47}`.
* **First-FLEX witnesses** (`rflex2.py`: genuine labels at every forced modulus `!= q0`, the stored lift at `q0`,
  exact budgets for all `v`, mask-AND on `F ∩ [3, q0]`): all 21 checked (`first4..7`, `rand4`, `bal4`, `last4`,
  `K = 4..7`) are good and non-coherent.
  * Using genuine labels at **all** odd `q` instead gives budget violations, so the collision-free labels at free
    moduli are needed. The write-up uses them.
* **I5 (minor).** For `first8` the log stores only the level-1 label at 179, with no lift. I found the lift
  `[29882, 17478, 9588, 4382, 22423, 1433, 17965, 19968]` (no level-2 collision) and verified the system exactly: 0
  violations, non-coherent. So `q_flex1(first8, 4) = 179` stands.

### 4.4 Independent check of a certificate (*solver*; `rmilp.py`, outputs `rmilp_*.txt`, `rmilpd_*.txt`)

`rmilp.py` sets up the same relaxation as the author's certificate but solves it differently:
* **Relaxation** (Lemma R items 1, 2, 4, 5): the moduli are `F ∩ [3, Q_sub]`, and at each modulus every individually
  admissible parity-respecting label is allowed. Labels are grouped by (signature, column-sign mask), with no
  dominance pruning.
* **Constraints:** the level-1 budgets (as `log` inequalities), plus non-coherence on `C = F ∩ [3, Q]` (for each of
  the `2^(r-1)` classes, some `q in C` whose mask misses it).
* **Solver:** HiGHS (`scipy.optimize.milp`). Its tolerances only enlarge the feasible set, so a reported
  infeasibility is not an artefact of rounding in the budgets. A reported solution is re-checked exactly.

The per-modulus type counts at the deep moduli (`3: 128`, `7: 512` with `3^3`, `7^2`) coincide with the author's
deep tables.

| instance | `Q_sub` | `Q` | depth | result |
|---|---|---|---|---|
| `first4`, `K = 4` | 79 | 61 | level 1 | **INFEASIBLE** (283 s): confirms `Q_core >= 61` independently of `cpc.c`/`cp_mix.py` |
| `first4`, `K = 4` | 79 | 67 | level 1 | feasible; solution exactly re-checked (0 violations, mask-AND 0); sanity: the test is not vacuous |
| `first4`, `K = 4` | 79 | 67 | `3^3, 7^2` | feasible (consistent with the witness at 67) |
| `first4`, `K = 5` | 131 | 109 | `3^3, 7^2` | **feasible** (765 s; exact re-check: 0 violations, non-coherent). The moduli `137..199` are needed for the author's certificate, consistent with their `mech_joint_first4_K5.txt` (outer 151 needed). Neither confirms nor refutes |
| `first4`, `K = 5` | 151 | 109 | `3^3, 7^2` | time limit (4000 s), no feasible point found; inconclusive. A feasible point was found within 765 s at `Q_sub = 131`, so this is weak evidence that adding `137..151` constrains the problem much more |

The per-modulus type counts for `K = 5` (`3: 128, 7: 512, 59: 256, 79: 4096, 97: 18816, 131: 27968`, ...) coincide
exactly with the author's `fj_first4_K5b.log`. These signature tables were computed independently.

The `K = 4` certificate survives this independent test. The deeper ones (`K = 5, 6, 7`, `first5`) I re-checked only
through the author's logs and by consistency (§4.2: every claimed witness is above its claimed core). They remain
single-implementation certificates, apart from the author's own Python cross-check at `K = 4`.

### 4.5 The witnesses fail the condition P3-all\* actually imposes (*exact* and *solver*; `rverify.py`, `slabbnb.py`, `slabcore.py`)

**Per-v angle compatibility** (exact, rational box criterion of §1), computed with the full all-level demand of
each joint witness:

| witness | per-unit | per-v angle-incompatible `v` (of `|V_K|`) | slab system at norm `<= K` |
|---|---|---|---|
| `bal4 K4 Q=11` | good | 2 (`±(1,1,2,0)`: `m = (1, 13589, 1, 9555)`, `m_1^2 + m_3^2 = 2.76 * 10^8 > 2N = 2.33 * 10^8`) | infeasible (implied) |
| `bal4 K5 Q=11` | good | 4 | infeasible (implied) |
| `first4 K4 Q=67` | good | 32 (e.g. `(-2,0,-1,-1)`: `m_1 = 157` with `157^2 = 24649 <= 2N = 24650`, so `x` is pinned and `m_3 = 43` fails) | infeasible (implied) |
| `first4 K5 Q=113` | good | 32 | infeasible (implied) |
| `first4 K6 Q=197` | good | 40 | infeasible (implied) |
| `first5 K4 Q=103` | good | 44 | infeasible (implied) |
| `rand4 K4 Q=11` | good | 2 (`±(1,0,0,1)`: `m_0^2 + m_2^2 = 65^2 + 87^2 = 11794 > N = 9169`) | infeasible (B&B, 133 nodes) |
| `rand4 K5 Q=11` | good | 4 | infeasible (implied) |
| `last4 K4 Q=11` | good | **0** | **infeasible** (B&B, 191 nodes) |
| `last4 K5 Q=11` | good | **0** | **infeasible** (B&B, 1489 nodes) |

**How the slab test works** (`slabbnb.py`). It is a complete interval branch-and-bound over `phi in [0, pi)^r`. A box
is pruned only when some `(v, s)` is violated on the **whole** box, with a `10^-9` margin in favour of feasibility.
So "infeasible" is a floating-point certificate with a safety margin, not an exact one.
* **Sanity check:** with the *genuine* demand on the same moduli (relaxed by `10^-7`), the search is never pruned to
  infeasibility (ABORT at `2 * 10^5` nodes). It does not invent obstructions.
* **The obstruction is small.** A deletion filter (`slabcore.py`) finds a **10-monomial** infeasible subsystem for
  `last4`, `K = 4`: `(1,0,0,0)`, `(1,0,0,-1)`, `(1,0,-1,0)`, `(0,2,0,0)`, `(1,-2,0,0)`, `(2,-1,0,0)`, `(1,-1,0,1)`,
  `(1,0,-1,1)`, `(1,1,1,0)`, `(2,0,1,-1)`, with demand vectors listed in the script output.

So **no joint witness of the write-up is the norm-`<= K` restriction of any P3-all\* member**, whatever angles are
chosen. Round 13's plan "take the flexible tops as the seed of C1 and add the slab conditions with angles chosen
after the labels" fails for every one of them at the first step. This is not a flaw in the write-up's arithmetic.
It is a flaw in reading per-unit witnesses as construction seeds.

### 4.6 What changes in the class "per-unit + angles" (*exact*; `onemod_angle.py`, `angflex.py`, `flexslab.py`)

I ran two refinements of the one-modulus test. ANG is the per-unit budgets plus per-v angle compatibility, a
necessary condition. SL is the per-unit budgets plus the slab system at norm `<= K` with one `phi`; it is stronger,
and is checked here only at `phi = arg pi_j`, in 40-digit arithmetic. At `phi = arg pi_j` the left side equals
`nu_s |D_s(Z_v)|/sqrt(N)`, so a genuine violation has size at least `nu_s/sqrt(N) > 10^-7`; the test flags anything
below `-10^-30`.

* **`first4`, `K = 4`.** Columns are all / outside the genuine class:

  | `q0` | per-unit survivors | ANG survivors |
  |---|---|---|
  | 83 | 36 / 24 | **12 / 0** |
  | 89, 97, 101, 103 | rigid | rigid |
  | 107 | 216 / 144 | 180 / 108 |
  | 109 | 198 / 162 | 108 / 72 |
  | 113 | 336 / 240 | 168 / 72 |
  | 127 | 6592 / 6080 | 4352 / 3840 |

  * So **83 is one-modulus RIGID in ANG** (exact certificate: level-1 relaxation, exact integers, exact rational
    angle test). None of the three stored per-unit FLEX witnesses at 83 is angle-compatible (`rflex4.py`).
  * At 107 three explicit systems (`angflex.py`; e.g. `ell = (1, 50, 61, 70)`, lift `[8533, 3506, 10213, 4930]`,
    genuine on `F \ {107}`) are per-unit good, angle-compatible, satisfy the slab system at `phi = arg pi_j`, and
    are non-coherent on `F ∩ [3, 107]`.
  * Hence, **in SL, `61 <= Q_core^sl(first4, 4) <= 103`**: the lower bound is the author's certificate, which
    transfers; the upper bound is my witness.
* **`first4`, `K = 5`, `q0 = 163`** (the per-unit `q_flex1`): per-unit 160 / 80, ANG **80 / 0**. So 163 is also
  one-modulus rigid under per-v angle compatibility (exact). The pattern repeats at `K = 5`: the first per-unit
  FLEX modulus is an artefact of the relaxation.
* **Growth `K = 4 -> 5` survives in SL.** `Q_core^sl(first4, 5) >= Q_core(first4, 5) >= 109 > 103 >=
  Q_core^sl(first4, 4)`. The per-unit `K = 5` and `K = 6` upper bounds do not transfer: the 113 and 197 joint
  witnesses and the 163 FLEX witnesses all violate ANG or SL. So in SL only lower bounds are known for `K >= 5`.
* **`r`-saturation survives in SL.** The first-FLEX witnesses of `first6`, `first7` (all 20 stored) and `first8`
  (my lift) at 179, `K = 4`, satisfy the slab system at `phi = arg pi_j` exactly. So `Q_core^sl <= 173` for
  `r = 6, 7, 8`, `K = 4`, as claimed.
* **Loose sets.**
  * The FLEX witnesses of `last4` at 31 (`K = 4`) and at 43 (`K = 5`, two of three) satisfy SL at
    `phi = arg pi_j`. So `7 <= Q_core^sl(last4, K) <= 29` resp. `<= 41`.
  * `rand4` (31) and `bal4` (41) at `K = 4`: the stored witnesses violate SL at `phi = arg pi_j`.
    * The ANG-compatible non-class survivors form a single `U_1`-orbit in each case (8 vectors; e.g. `(1,1,13,18)`
      at 31 and `(1,33,3,22)` at 41). `U_1` scaling preserves signatures, hence demands.
    * Their slab systems are **infeasible** (`onemod_slab.py`, complete B&B, 21 and 49 nodes).
    * So, in SL, 31 (`rand4`) and 41 (`bal4`) are one-modulus **rigid** (solver). The loose one-modulus segments
      also move up once angles are imposed.
  * So the claim "the loose core is exactly 7 (vacuous)" is **unproved** in SL. What is proved is
    `Q_core^sl(last4, 4) <= 29` and `Q_core^sl(last4, 5) <= 41`.

## 5. Implications for Rig(K) and the barrier logic

1. **Growth branch (Rig(K) growing ⇒ P3-all\* reduces to P3-G in the regime ⇒ only L7 remains).**
   * This branch needs only certificates, and the certificates transfer.
   * The finite evidence "the per-unit core of `first4` grows `61 -> 109 -> 193 (-> >= 239)`" is correct as a
     statement about the relaxation.
   * In the stronger class SL, the strict growth `K = 4 -> 5` is now also exact (§4.6). Beyond that only lower
     bounds exist, because the SL core at `K = 4` may itself be larger than 61.
   * Nothing here approaches the regime: `r = 4, 5`, `K <= 7`, and the fit to `M^(floor(K/2))` rests on 3-4 points
     whose ratios fluctuate (0.37, 0.66, 0.37). The write-up labels this evidence. I agree, and add that it is
     evidence about the per-unit relaxation.
2. **Construction branch (stall ⇒ flexible tops seed C1).**
   * The write-up's loose-prime witnesses do not seed C1. Eight of ten violate per-v angle compatibility, so no
     angle choice can complete them. The other two (`last4`) have an infeasible slab system at norm `<= K`
     (§4.5).
   * Round 13's plan ("add the slab conditions with angles chosen after the labels") therefore fails at its first
     step for every exhibited top.
   * The stall itself (`Q_core = 7`) is proved only for the per-unit class. In SL it is open between 7 and the
     one-modulus bounds (29 and 41 for `last4`).
   * So this round gives **no** C1 seed. The barrier map is unchanged: P3-all\* is open, uncapped, with no growth.
3. **What the class issue suggests (heuristic, labelled).** The slab coupling is itself a rigidity mechanism. At
   `first4`, `K = 4`, it moves the one-modulus segment from 79 to 103. Any continuation of Rig(K) should state the
   class as "parity-respecting labels **and** angles `phi` satisfying (1.2) for `||v||_1 <= K`", i.e. the norm-`K`
   truncation of P3-all\* itself.
   * Witnesses must then come with a verified `phi`. The one-modulus witnesses at `phi = arg pi_j` are a cheap
     source, valid wherever the extra demand fits under `|D_s|`.
   * Certificates may keep the per-unit relaxation, or add per-v box compatibility, which is exact and rational.
   * Rig(K) as worded in round 13 §7.2 ("meet the per-unit budget") should be corrected accordingly. The per-unit
     wording is also what makes the "stall ⇒ C1" inference invalid.
4. **No growth improvement, uniform bound or cap** follows from the write-up or from this report. The bracket
   `3 M log M - O(M) <= W_min^all <= O(M log^2 M/loglog M)` and `M <= (2/3 + eps) log R/loglog R` are unchanged.

## 6. Corrected statement

**Status: evidence** about Rig(K). Exact finite results, valid only in the stated classes.

Let **PU** be the write-up's class: parity-respecting labels, all levels, per-unit budgets
`m_(v,s)^2 <= Norm Z_v (2 if s odd)` for `1 <= ||v||_1 <= K`. Let **SL** be PU plus the slab system (1.2) of
`sharp.md` for those `v` with one `phi in R^r`. SL is necessary for a label system to lie under a P3-all\* member,
and PU is a strict relaxation of SL.

1. **Lower bounds, exact for PU and hence valid for SL and P3-all\*.** `Q_core(first4, K) >= 61, 109, 193, 239` for
   `K = 4, 5, 6, 7`, and `Q_core(first5, K) >= 101, 163` for `K = 4, 5`.
   * They come from complete searches over the Lemma R relaxation.
   * The `first4`, `K = 4` bound is confirmed independently by MILP (solver evidence). The others are
     single-implementation certificates.
   * By Lemma M (side condition verified), the all-norm core is at least these values.
2. **Upper bounds, exact for PU only.** `Q_core^PU = 61, 109, 193` (`first4`, `K = 4, 5, 6`); `101` (`first5`,
   `K = 4`); `7` for `rand4`, `bal4`, `last4` at `K = 4, 5`; `<= 173` for `first6..8`, `K = 4`. All witnesses are
   verified.
   * **None of the 10 joint witnesses lies in SL.** Eight violate per-v angle compatibility (exact). The two
     `last4` witnesses have infeasible slab systems (complete branch-and-bound).
3. **In SL** (exact, from verified witnesses at `phi = arg pi_j`, plus the transferred lower bounds):
   * `61 <= Q_core^SL(first4, 4) <= 103`; hence `Q_core^SL(first4, 5) >= 109 > Q_core^SL(first4, 4)`;
   * `Q_core^SL <= 173` for `first6..8`, `K = 4`;
   * `7 <= Q_core^SL(last4, 4) <= 29` and `7 <= Q_core^SL(last4, 5) <= 41`;
   * at `first4`, `q = 83` (`K = 4`) and `q = 163` (`K = 5`) are one-modulus rigid under per-v angle
    compatibility (exact). So the per-unit `q_flex1` values 83 and 163 are artefacts of the relaxation.
4. **Not established:**
   * the loose-prime "stall" in SL or P3-all\*;
   * any C1 seed, since no exhibited flexible top admits angles;
   * the scale law "core `~ 0.4-0.7 |S_floor(K/2)|`" beyond PU at `r = 4, 5`;
   * any statement in regime (1.4).
5. **Minor corrections.**
   * Lemma K proves `K`-independence of local admissibility only for `q <= p_min^((K+1)/2)`. "Local rigidity ends at
     31 for all `K`" is checked for `q <= 113`.
   * "local `<=` joint core" is an observed comparison of initial segments, not an implication.
   * Local rigidity is not vacuous for `rand4` (`q = 17`).
   * The `first8` FLEX witness needed a lift, which I supplied and verified.
   * Certificates with depth need `X >= 27, 49, 121, 169`.

No growth improvement, uniform bound or cap is claimed or implied. The P3-all\* bracket and `two_thirds` are
unchanged.

## 7. Scripts (`referee_rigcompute_checks/`)

| script | purpose | output |
|---|---|---|
| `reflib.py` | own exact helpers: Gaussian primes, torus logs, `D_s`, budgets, column-sign masks, exact angle-box test | -- |
| `t_sanity.py` | parity law and `collision <-> q \| D_s` (21,888 checks) | stdout |
| `rfree.c`, `rforced.py` | forced/free decision, complete DFS; comparison with `thr_*.json` | `rforced_*.json`, `rforced.log` |
| `rfreevec.py` | re-check of the author's free vectors and scan coverage (`r = 4` at `K = 6, 7`; `first6..8`) | stdout |
| `rverify.py` | independent all-level verifier of joint witnesses, plus per-v angle test | `rverify_all.txt` |
| `rflex2.py`, `rflex3.py`, `rflex4.py`, `rflex8.py` | one-modulus FLEX witnesses (genuine on `F \ {q0}`), angle test, `first8` lift | stdout |
| `onemod_angle.py` | exhaustive one-modulus survivors in PU and ANG | `om_first4_K4.txt`, `om_loose_K4.txt` |
| `angflex.py` | SL-valid FLEX witnesses at `first4`, `K = 4`, `q0 = 107` | stdout |
| `flexslab.py` | slab slack at `phi = arg pi_j` for FLEX witnesses | stdout |
| `rmilp.py` | independent MILP (HiGHS) joint test, optional depth | `rmilp_*.txt`, `rmilpd_*.txt` |
| `onemod_slab.py` | slab B&B for one-modulus systems (loose ANG survivors) | stdout |
| `slab.py`, `slabbnb.py`, `slabcore.py` | slab system search, complete branch-and-bound, minimal infeasible subsystem | stdout |
| `repair.py` | attempted one-modulus repair of loose witnesses into ANG (none found up to `q = 59`; evidence only) | stdout |
| `notes.md` | working notes | -- |
