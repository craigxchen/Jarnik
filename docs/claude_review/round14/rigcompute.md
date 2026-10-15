# Round 14 — first finite test of Rig(K): rigid cores for r = 4..8, K = 4..7

(round 14; unrefereed; scripts and outputs in `rigcompute_checks/`, run from that directory)

## 0. Summary of results

Labels: *proved* (all quantifiers, here), *exact* (finite computation, certificate or verified witness, script named),
*evidence*. **No growth improvement, no uniform bound and no new cap is claimed.** The bracket for P3-all\* and
`M <= (2/3 + eps) log R/loglog R` are unchanged. Everything concerns label systems in the Rig(K) setting of round 13
§7.2 (parity-respecting labels, per-unit budgets for every `||v||_1 <= K`, all levels); actual circles give the genuine
system, which satisfies everything here.

1. **The first finite test of Rig(K) is done exactly on both sides** for `r = 4` (`K = 4, 5, 6`; lower bound at
   `K = 7`) and `r = 5` (`K = 4`; lower bound at `K = 5`), with one-modulus data and exact forced sets for `r = 4..8`
   (`r = 8`: `first8`, `K = 4` only; 26 moduli above 499 of unknown status, one-modulus scan up to `q = 353`). The rigid core `Q_core(P, K)` (largest `Q`
   such that every good system is system-level coherent with one column-sign conjugate of the genuine labels on
   `F_K^ex ∩ [3, Q]`) is certified from below by a complete search over an exact relaxation (Lemma R) and from above
   by an explicit good system, verified at every level by two independent checkers. New exact tools: a C engine for
   single-modulus label problems, exact multi-level signature tables, and a complete joint search (`cpc.c`) whose
   infeasibility answers are certificates.
2. **Three nested rigidity notions, all exact** (§5): local (intrinsic) rigidity `<=` joint core `<=` one-modulus
   segment. *Lemma K (proved):* local rigidity is independent of `K` once `p_min^(K+1) >= q^2`, so it cannot grow
   with `K`; it ends at `q = 31` (`first4`), `11` (`first5-7`), and is vacuous for the other sets.
3. **Tight primes (`first`, the `r` smallest split primes): the joint core grows with `K` and `r`.**
   `Q_core = 61, 109` (`K = 4, 5`), `193` (`K = 6`), `>= 239` (`K = 7`), `r = 4`; `Q_core = 101` (`K = 4`) and `>= 163` (`K = 5`) for `r = 5`.
   In units `n = m_q/4`: `15, 27, 48, >= 60` (`r = 4`) and `25, >= 41` (`r = 5`, `K = 4, 5`), i.e. `0.37-0.66 |S_floor(K/2)|`: the core
   tracks the pigeonhole scale `M^(floor(K/2))` (§6), a stable fraction `0.69-0.78` (in `n`) of the one-modulus segment, while
   the forced range `n*` runs ahead at even `K`. *Caveat (exact):* at `K = 4` the one-modulus segment, an upper bound
   for the core, stops growing at `r = 6`: `q_flex1 = 179` for `first6`, `first7`, `first8` (the added primes 41, 53,
   61 are no longer small), so the core is `<= 173` (`n <= 43`) there and its `r`-growth at fixed `K` saturates.
   Mechanism (exact relaxation analysis): the forced collisions of the
   smallest moduli (`q = 3, 7`: every monomial collides) plus a few mid moduli, acting through the norm-3..K budgets
   jointly; no single modulus or norm carries it.
4. **Loose primes (`rand4`, `bal4`, `last4`: all `p_j >= 53`): the joint core stalls at the vacuous moduli** (`n_q <= 2`)
   for `K = 4, 5` (six exact cases), although one-modulus rigidity reaches `q = 29-61`. Explicit verified witnesses are
   good for all norms `<= K` and lie outside every column-sign class at `34/36` to `66/67` of the forced moduli above `q = 11`: flexible
   tops covering the whole forced range (C1 seeds at norm `<= K`). At `K = 6` the exact relaxations still admit
   non-coherent systems at `q = 11` (no certificate beyond the vacuous moduli), but no witness was completed
   (deep 3-adic collisions at `q = 3`; §7).
5. **Reading for Rig(K) (evidence).** Rig(K) holds in its "growing" form only where the budgets are tight
   (`log p_j/log r` near 1, where near-extremal clusters put their primes); for loose primes the local rigidity
   mechanism is absent at the tested norms and the flexible tops seed a C1 construction. Neither branch is decided in
   regime (1.4): the core growth is a finite-`r` fit, and the loose-prime systems are not P3-all\* members (they are
   good only up to norm `K`; the all-norm core is `>= sup_K Q_core(K)`, Lemma M).

## 1. The exact finite problem (class, regime, parameters)

Conventions are those of `round10/L4K.md` §1, re-implemented independently (`rigcompute_checks/rl.py`; nothing is
imported from the repository).

* **Columns.** `r` distinct split primes `p_j = 1 mod 4`, `p_j <= 400`, with the Gaussian primes `pi_j = a + bi`
  (`a` odd, `b` even, both positive). `Z_v = prod pi_j^(v_j)` (negative exponents meaning `conj pi_j`),
  `D_s(Z_v) in {Y, X - Y, X, X + Y}` for `Z_v = X + iY`.
* **Labels.** At an odd prime `q` not dividing `N`: `ell(q) in (Z/m_q)^r`, `m_q = q - chi_(-4)(q)`, `n_q = m_q/4`,
  **parity-respecting**: `ell_j = [(p_j/q) = -1] mod 2`. A monomial `v` collides at level 1 with unit `s` iff
  `<v, ell> = s * iota_q mod m_q` (`iota_q = log i`), at level `a` with the lifted labels.
* **Good up to norm `K`.** For every `v` with `1 <= ||v||_1 <= K` and every unit `s`, the product `m_(v,s)` of all
  colliding prime powers satisfies `m_(v,s)^2 <= Norm(Z_v) * (2 if s odd)`, i.e. `log m_(v,s) <= w(v)/2 + [s odd] log2/2`
  (exact integer test throughout).
* **Genuine labels** `L_j(q) = log(pi_j/conj pi_j)` satisfy this at every level: `v` collides at `q^a` with unit
  `s` iff `q^a | D_s(Z_v)` (re-checked at every modulus used, `deform14.py::modulus_data`).
* **Genuine class at `q`** (level 1, as in L4K Lemma P): `c = ell mod n_q` equals `u L(q) mod n_q` for a unit `u`.
  **Column-sign classes**: `u (eps * L(q))` with `eps in {+-1}^r` (conjugating `pi_j` in column `j`); the
  global map `v -> eps * v` preserves all budgets, so the `eps`-conjugate of a good system is good. *System-level
  coherence with the genuine labels on a set `C` of moduli* means one `eps` for all `q in C`.
* **Universe of moduli.** The exact forced set `F_K^ex`: `q` is *forced* if every parity-respecting label vector has a
  level-1 collision of some `v in V_K`. At a free modulus a collision-free vector exists, has zero demand, and can
  replace any label there without increasing any `m_(v,s)`. Hence: **a good system with a prescribed behaviour on
  forced moduli exists iff one exists with collision-free labels at all free moduli.** All tests below therefore
  live on `F_K^ex` with zero demand elsewhere, and every explicit witness is completed by the verified
  collision-free vectors of §2 (and by Theorem 1(b) vectors of L4K above the scanned range).

**Definitions used in the tables.** For a prime set `P` and norm `K`:
* `q*(K)` = largest forced prime (exact), `n* = n_(q*)`.
* **One-modulus test at `q0 in F_K^ex`**: genuine labels at every other forced modulus (all levels), arbitrary
  parity-respecting `ell` at `q0`. `q0` is **RIGID** if every good such system has `ell(q0)` in the genuine class,
  **FLEX** if an explicit good system (verified exactly, all levels) has `ell(q0)` outside it, **VACUOUS** if the genuine
  class contains every parity-respecting label (`n_q <= 2`; referee m4).
* `q_flex1(K)` = the smallest FLEX modulus. **Exact consequence:** the rigid core `Q_core(K)` (the largest `Q` such that
  every good system is coherent with the genuine one on `F_K^ex ∩ [3, Q]`) satisfies `Q_core(K) < q_flex1(K)`,
  because the FLEX witness is itself a good system (genuine elsewhere) that is not coherent at `q_flex1`.
  `q_seg(K)` = the largest forced prime below `q_flex1` (end of the one-modulus rigid initial segment).
* **Joint tests** (everything free) give lower bounds on `Q_core`; see §3.3.

**Parameters actually run.** `r = 4..8`; prime sets per `r` (`sets.py`): `first` (the `r` smallest split primes:
tightest budgets), `rand` (random `r`-subset of the split primes `<= 400`, seed `1400 + r`), `bal` (`r` consecutive
split primes from 97: equal weights, the Plotkin-like near-extremal profile with all `log p_j` equal up to 13%),
`last` (the `r` largest split primes `<= 400`: loosest budgets). `K = 4..7`. Moduli: all of `F_K^ex` (scanned up
to `Q_scan = 4 A_K + 5`, Theorem 1(b) of L4K with greedy `B_h` sets, so that every prime above `Q_scan` is free).

## 2. Harness and certification method

All scripts are in `rigcompute_checks/` (run from that directory; `gcc -O2 -o csp csp.c`, same for `freescan.c`,
`deepsig.c`). Nothing imports from the repository; the round-10 harness was read, copied to `orig/` for reference, and
re-implemented.

| component | what it decides | exactness |
|---|---|---|
| `freescan.c` + `thresholds.py` | forced/free status of every prime `q <= Q_scan` for `V_K`; DFS over `(Z/n)^r`, column 0 over gcd-orbits | exact (complete DFS); every FREE vector re-verified in numpy; node-capped runs report `ABORT` explicitly |
| `csp.c` (modes 0/1/2) + `deform14.py` | one-modulus deformations: all parity-respecting `ell in (Z/m)^r` at `q0` (column 0 over orbits of `U_1 = {u in (Z/m)^*: u = 1 mod 4}`, which preserves the parity, the unit index of every collision and the genuine class), forbidden residues = units `s` with `(q0 * G_other(v,s))^2 > cap(v,s)` | complete DFS = exact certificate of RIGID (level-1 relaxation, hence valid at all levels); FLEX witnesses re-verified independently in Python with exact integers and lifted to level 2 with no level-2 collision |
| `csp.c` (mode 3) + `sigtab.py` | per-modulus tables of all achievable level-1 collision signatures of *individually admissible* labels (`q^2 <= cap(v,s)` for each collision), with the column-sign classes realised | exact enumeration with orbit weights |
| `deepsig.c` + `deep.py` | the same with labels mod `m q^(a-1)` (multi-level signatures with multiplicities) at the smallest moduli | exact |
| `loc.py` | LOCAL rigidity: every admissible label at `q` lies in `class(L)` (LOC) / in some column-sign class (LOCe) | exact (complete DFS) |
| `cp_joint.py` | JOINT test, every forced modulus free: complete DFS over per-modulus minimal signatures with capacity propagation | `INFEASIBLE` is an exact certificate (§3.3); `FOUND` is verified exactly or reported as unlifted |
| `joint_lp.py`, `joint_ilp.py`, `ls_joint.py` | Lagrangian/LP certificates; MILP and min-conflicts search | LP certificates verified in interval arithmetic; MILP/LS failures are evidence only |

**Validation against round 10** (exact agreement): forced sets `F_4^ex`, `F_5^ex` of `{5,13,17,29}` and
`{173,181,193,197}` (`q* = 199, 349, 271, 607`); one-modulus survivor counts (all / outside the genuine class)
`47: 8/0, 83: 36/24, 97: 16/0, 101: 20/0, 103: 48/0, 107: 216/144` at norm 4 and
`107: 18/0, 113: 24/0, 127: 32/0, 151: 36/0, 167: 120/96, 199: 120/40, 211: 884/832` at norm 5 — the referee's
table (`referee_L4K.md` §5 item 3, §6.3) to the last digit.

**Lemma R (the reductions behind every certificate; proved).** Fix `P`, `K`, and write `F = F_K^ex`. Call a
*system* an assignment of parity-respecting labels (all levels) to every odd prime `q` not dividing `N`; it is good
if every `(v, s)`, `v in V_K`, meets its budget.
1. *(free moduli)* Replacing the label at a free modulus by a collision-free one (which exists by definition, and
   has no collision at any level since a level-`a` collision implies a level-1 collision) changes no
   `m_(v,s)` except by removing factors. So a good system with any prescribed labels on `F` exists iff one exists
   with collision-free labels off `F`.
2. *(level-1 relaxation)* Let `Good_1` be the set of label choices on `F` whose **level-1** demand
   `prod_q q^([v collides at q with s])` meets every budget. Every good system restricts to an element of `Good_1`.
   Hence "no element of `Good_1` has property X" certifies "no good system has property X".
3. *(minimal signatures)* The level-1 demand of a choice depends only on the collision signatures
   `sigma_q = {(v,s)}` and is monotone in each `sigma_q`. If a choice in `Good_1` uses `sigma_q`, and some label at `q`
   has signature `sigma' ⊆ sigma_q`, the modified choice is still in `Good_1`. So for any property X that only
   constrains the labels on a set `C`, it suffices to let the moduli outside `C` range over inclusion-minimal
   achievable signatures, and the moduli in `C` over signatures minimal among those realising X (with the realised
   column-sign masks, dominance `sigma' ⊆ sigma` and `mask' ⊆ mask` for the coherence test).
4. *(individual admissibility)* A label at `q` with a collision `(v,s)` such that `q^2 > Norm(Z_v) (2 if s odd)`
   is in no good system. Tables may therefore be restricted to labels all of whose collisions are individually
   admissible.
5. *(dropping moduli)* Removing a modulus from the joint problem removes its (non-negative) demand: a relaxation.
6. *(upper bounds)* An explicit good system (all levels, all moduli: labels on `F` with exact all-level demand,
   collision-free vectors elsewhere) that is not coherent on `F ∩ [3, Q]` proves `Q_core < Q`; in particular a FLEX
   one-modulus witness at `q0` proves `Q_core < q0`.

*Proof.* Each item is immediate from the definitions: all demands are products of prime powers attached to
collisions, the budget is an upper bound on a product, and collisions at a deeper level imply collisions at every
shallower level. In 4, a single collision already costs `log q > B(v,s)`. QED

**Lemma M (monotonicity in `K`; proved).** Let `Q = Q_core(K)` and suppose `F_(K+1)^ex ∩ [3, Q] = F_K^ex ∩ [3, Q]`
(in every instance below both equal the set of all primes `<= Q` not dividing `2N`). Then `Q_core(K+1) >= Q_core(K)`.
*Proof.* A system good up to norm `K + 1` is good up to norm `K` (it meets a subset of the budgets), hence coherent on
`F_K^ex ∩ [3, Q] = F_(K+1)^ex ∩ [3, Q]`. QED. So the all-norm core that P3-all\* imposes,
`Q_core(infinity) = sup_K Q_core(K)`, is bounded below by every finite-`K` core.

**Lemma K (local admissibility is K-stable; proved).** Let `q` be a prime and `k(q)` the largest `||v||_1` among
`v` with `Norm(Z_v) < q^2`. If `K >= k(q)` (e.g. if `p_min^(K+1) >= q^2`), the set of individually admissible
parity-respecting labels at `q` for norm `K'` is the same for every `K' >= K`. *Proof.* A monomial of norm `> K`
has `Norm(Z_v) >= q^2` (by the definition of `k(q)`, or because `Norm(Z_v) >= p_min^(K+1) >= q^2`), so none of its collisions is inadmissible; the constraint set does not grow.
QED

So **local (intrinsic) rigidity** at `q` — every individually admissible label lies in a column-sign class — is a
property of `q` and `P` alone once `K >= k(q)`; adding norms cannot extend it. Any growth of the rigid core with
`K` must come from the **joint** budget interaction (the demand that the forced collisions of the other moduli put
on shared `(v, s)`), which is what §3.3 measures.

## 3. Results: rigid cores per (prime set, K)

### 3.1 Exact forced sets `F_K^ex` (*exact*; `thresholds.py`, outputs `thr_<set>_K<K>.json`)

Every prime `q <= Q_scan = 4 A_K + 5` (`q` not dividing `2N`) was decided by complete DFS; above `Q_scan` Theorem 1(b)
of L4K applies. Columns: `|S_h|` (`h = floor(K/2)`) is the pigeonhole threshold (forced if `n_q < |S_h|`), `A_K` the
Theorem 1(b) threshold (free if `n_q > A_K`, greedy `B_h` set). "Status unknown" moduli hit the node cap `3*10^7`; they
are treated as zero-demand in rigidity certificates (a relaxation) and as genuine-labelled in witnesses (a legitimate
label), so every claim below is valid either way.

| set | r | K | P | |F_K^ex| | q* | n* | pigeonhole |S_h| | Thm 1(b) A_K | scanned to | status-unknown (node cap) |
|---|---|---|---|---|---|---|---|---|---|---|
| bal4 | 4 | 4 | [97, 101, 109, 113] | 40 | 271 | 68 | 41 | 1636 | 6549 | - |
| bal4 | 4 | 5 | [97, 101, 109, 113] | 68 | 367 | 92 | 41 | 2935 | 11745 | - |
| bal4 | 4 | 6 | [97, 101, 109, 113] | 113 | 809 | 202 | 129 | 8310 | 33245 | - |
| bal4 | 4 | 7 | [97, 101, 109, 113] | 174 | 1913 | 478 | 129 | 12649 | 50601 | - |
| first4 | 4 | 4 | [5, 13, 17, 29] | 40 | 199 | 50 | 41 | 1636 | 6549 | - |
| first4 | 4 | 5 | [5, 13, 17, 29] | 65 | 349 | 87 | 41 | 2935 | 11745 | - |
| first4 | 4 | 6 | [5, 13, 17, 29] | 115 | 1039 | 260 | 129 | 8310 | 33245 | - |
| first4 | 4 | 7 | [5, 13, 17, 29] | 173 | 2089 | 522 | 129 | 12649 | 50601 | - |
| last4 | 4 | 4 | [353, 373, 389, 397] | 44 | 263 | 66 | 41 | 1636 | 6549 | - |
| last4 | 4 | 5 | [353, 373, 389, 397] | 71 | 383 | 96 | 41 | 2935 | 11745 | - |
| last4 | 4 | 6 | [353, 373, 389, 397] | 113 | 1231 | 308 | 129 | 8310 | 33245 | - |
| last4 | 4 | 7 | [353, 373, 389, 397] | 174 | 1231 | 308 | 129 | 12649 | 50601 | - |
| rand4 | 4 | 4 | [53, 61, 89, 173] | 40 | 199 | 50 | 41 | 1636 | 6549 | - |
| rand4 | 4 | 5 | [53, 61, 89, 173] | 69 | 607 | 152 | 41 | 2935 | 11745 | - |
| rand4 | 4 | 6 | [53, 61, 89, 173] | 113 | 929 | 232 | 129 | 8310 | 33245 | - |
| rand4 | 4 | 7 | [53, 61, 89, 173] | 174 | 1879 | 470 | 129 | 12649 | 50601 | - |
| bal5 | 5 | 4 | [97, 101, 109, 113, 137] | 54 | 293 | 73 | 61 | 2636 | 10549 | - |
| bal5 | 5 | 5 | [97, 101, 109, 113, 137] | 101 | 647 | 162 | 61 | 4735 | 18945 | - |
| first5 | 5 | 4 | [5, 13, 17, 29, 37] | 54 | 293 | 73 | 61 | 2636 | 10549 | - |
| first5 | 5 | 5 | [5, 13, 17, 29, 37] | 104 | 719 | 180 | 61 | 4735 | 18945 | - |
| last5 | 5 | 4 | [349, 353, 373, 389, 397] | 60 | 337 | 84 | 61 | 2636 | 10549 | - |
| last5 | 5 | 5 | [349, 353, 373, 389, 397] | 102 | 617 | 154 | 61 | 4735 | 18945 | - |
| rand5 | 5 | 4 | [89, 109, 149, 281, 353] | 55 | 293 | 73 | 61 | 2636 | 10549 | - |
| rand5 | 5 | 5 | [89, 109, 149, 281, 353] | 102 | 641 | 160 | 61 | 4735 | 18945 | - |
| bal6 | 6 | 4 | [97, 101, 109, 113, 137, 149] | 75 | 599 | 150 | 85 | 4236 | 16949 | - |
| bal6 | 6 | 5 | [97, 101, 109, 113, 137, 149] | 147 | 967 | 242 | 85 | 7615 | 30465 | - |
| first6 | 6 | 4 | [5, 13, 17, 29, 37, 41] | 74 | 433 | 108 | 85 | 4236 | 16949 | - |
| first6 | 6 | 5 | [5, 13, 17, 29, 37, 41] | 150 | 1039 | 260 | 85 | 7615 | 30465 | - |
| last6 | 6 | 4 | [337, 349, 353, 373, 389, 397] | 75 | 601 | 150 | 85 | 4236 | 16949 | - |
| last6 | 6 | 5 | [337, 349, 353, 373, 389, 397] | 148 | 1039 | 260 | 85 | 7615 | 30465 | - |
| rand6 | 6 | 4 | [37, 101, 257, 317, 373, 389] | 76 | 449 | 112 | 85 | 4236 | 16949 | - |
| rand6 | 6 | 5 | [37, 101, 257, 317, 373, 389] | 147 | 911 | 228 | 85 | 7615 | 30465 | [919, 1801] |
| bal7 | 7 | 4 | [97, 101, 109, 113, 137, 149, 157] | 92 | 599 | 150 | 113 | 6236 | 24949 | [499, 541, 547, 557, 563, 571, 673] |
| first7 | 7 | 4 | [5, 13, 17, 29, 37, 41, 53] | 92 | 599 | 150 | 113 | 6236 | 24949 | [499, 541, 547, 557, 563, 571] |
| rand7 | 7 | 4 | [5, 41, 53, 149, 173, 257, 313] | 88 | 509 | 127 | 113 | 6236 | 24949 | [499, 541, 547, 557, 563, 571, 769] |
| first8 | 8 | 4 | [5, 13, 17, 29, 37, 41, 53, 61] | 97 | 631 | 158 | 145 | 9036 | 36149 | [499, 541, 547, 557, 563, 571, 587, 613, 619, 641, 643, 647, 653, 659, 661, 677, 683, 691, 701, 719, 727, 733, 743, 751, 769, 863] |
| rand8 | 8 | 4 | [61, 97, 173, 197, 229, 337, 353, 373] | 99 | 761 | 190 | 145 | 9036 | 36149 | [499, 541, 547, 557, 563, 571, 587, 613, 619, 643, 647, 653, 659, 661, 677, 683, 691, 701, 719, 727, 733, 743, 751] |


(r = 7, 8 at `K >= 5` and r = 8 at `K = 4` were not completed: the exhaustive proofs of forcedness near the boundary
cost about an hour per set at `r = 7`, `K = 4`.)

* The boundary `n*` sits at `1.2-4.0 |S_h|`, far below `A_K` (as in L4K: `n* ~ 2.5r` at `K = 2`). It jumps between
  odd and even `K` (`h = floor(K/2)`): `n*(K=6)/n*(K=5) = 1.5-3.2` for `r = 4`.

### 3.2 One-modulus rigidity (*exact*; `deform14.py`, outputs `def_wit_<set>_K<K>.json`, logs `log_def_*`)

`q_flex1` = first FLEX modulus (an explicit good system, genuine elsewhere, verified at all levels); hence
**`Q_core < q_flex1` exactly**. "Segment" = all forced `q < q_flex1` are RIGID (or vacuous).

| set | r | P | K | |F| | q* (n*) | one-mod rigid segment: last q (n) | first flexible q | rigid above it | #rigid/#flex | complete |
|---|---|---|---|---|---|---|---|---|---|---|
| bal4 | 4 | [97, 101, 109, 113] | 4 | 40 | 271 (68) | 37 (9) | 41 | [] | 8/29 | yes |
| bal4 | 4 | [97, 101, 109, 113] | 5 | 68 | 367 (92) | 61 (15) | 67 | [] | 14/51 | yes |
| bal4 | 4 | [97, 101, 109, 113] | 6 | 113 | 809 (202) | 83 (21) | 89 | [] | 19/91 | yes |
| bal4 | 4 | [97, 101, 109, 113] | 7 | 174 | 1913 (478) | 127 (32) | 131 | [] | 23/148 | yes |
| first4 | 4 | [5, 13, 17, 29] | 4 | 40 | 199 (50) | 79 (20) | 83 | [89, 97, 101, 103] | 19/19 | yes |
| first4 | 4 | [5, 13, 17, 29] | 5 | 65 | 349 (87) | 157 (39) | 163 | [181, 191] | 32/31 | yes |
| first4 | 4 | [5, 13, 17, 29] | 6 | 115 | 1039 (260) | 263 (66) | 269 | [271, 277, 281, 283] | 53/60 | yes |
| first4 | 4 | [5, 13, 17, 29] | 7 | 173 | 2089 (522) | 331 (83) | 337 | [347, 353, 359, 373, 379, 389, 401, 409, 439] | 69/102 | yes |
| last4 | 4 | [353, 373, 389, 397] | 4 | 44 | 263 (66) | 29 (7) | 31 | [] | 6/35 | yes |
| last4 | 4 | [353, 373, 389, 397] | 5 | 71 | 383 (96) | 41 (10) | 43 | [73] | 10/58 | yes |
| last4 | 4 | [353, 373, 389, 397] | 6 | 113 | 1231 (308) | 41 (10) | 43 | [47, 73] | 11/99 | yes |
| last4 | 4 | [353, 373, 389, 397] | 7 | 174 | 1231 (308) | 53 (13) | 59 | [67, 71, 73] | 15/156 | yes |
| rand4 | 4 | [53, 61, 89, 173] | 4 | 40 | 199 (50) | 29 (7) | 31 | [37, 43, 47] | 9/28 | yes |
| rand4 | 4 | [53, 61, 89, 173] | 5 | 69 | 607 (152) | 47 (12) | 59 | [] | 11/55 | yes |
| rand4 | 4 | [53, 61, 89, 173] | 6 | 113 | 929 (232) | 79 (20) | 83 | [101, 103, 113] | 19/91 | yes |
| rand4 | 4 | [53, 61, 89, 173] | 7 | 174 | 1879 (470) | 79 (20) | 83 | [97, 101, 103, 107, 109, 113, 127, 139, 149, 151] | 26/145 | yes |
| bal5 | 5 | [97, 101, 109, 113, 137] | 4 | 54 | 293 (73) | 53 (13) | 59 | [61] | 13/38 | yes |
| bal5 | 5 | [97, 101, 109, 113, 137] | 5 | 101 | 647 (162) | 107 (27) | 127 | [] | 22/76 | yes |
| first5 | 5 | [5, 13, 17, 29, 37] | 4 | 54 | 293 (73) | 127 (32) | 131 | [137, 139, 149, 193] | 27/25 | yes |
| first5 | 5 | [5, 13, 17, 29, 37] | 5 | 104 | 719 (180) | 293 (73) | 307 | [311, 313, 317] | 57/45 | yes |
| last5 | 5 | [349, 353, 373, 389, 397] | 4 | 60 | 337 (84) | 37 (9) | 41 | [73] | 9/48 | yes |
| last5 | 5 | [349, 353, 373, 389, 397] | 5 | 102 | 617 (154) | 53 (13) | 59 | [61, 67, 71, 73] | 16/83 | yes |
| rand5 | 5 | [89, 109, 149, 281, 353] | 4 | 55 | 293 (73) | 47 (12) | 53 | [] | 11/41 | yes |
| rand5 | 5 | [89, 109, 149, 281, 353] | 5 | 102 | 641 (160) | 71 (18) | 73 | [79, 83] | 18/81 | yes |
| bal6 | 6 | [97, 101, 109, 113, 137, 149] | 4 | 75 | 599 (150) | 89 (22) | 103 | [] | 20/52 | yes |
| bal6 | 6 | [97, 101, 109, 113, 137, 149] | 5 | 147 | 967 (242) | 179 (45) | 181 | [191, 193] | 33/111 | yes |
| first6 | 6 | [5, 13, 17, 29, 37, 41] | 4 | 74 | 433 (108) | 173 (43) | 179 | [181, 191, 193, 197, 199, 211, 223] | 38/34 | yes |
| first6 | 6 | [5, 13, 17, 29, 37, 41] | 5 | 150 | 1039 (260) | 409 (102) | 419 | [421, 431, 433, 439, 443, 449, 457, 461, 467, 491, 599] | 82/66 | yes |
| last6 | 6 | [337, 349, 353, 373, 389, 397] | 4 | 75 | 601 (150) | 53 (13) | 59 | [61, 67] | 14/58 | yes |
| last6 | 6 | [337, 349, 353, 373, 389, 397] | 5 | 148 | 1039 (260) | 127 (32) | 131 | [139] | 28/117 | yes |
| rand6 | 6 | [37, 101, 257, 317, 373, 389] | 4 | 76 | 449 (112) | 79 (20) | 83 | [] | 17/56 | yes |
| rand6 | 6 | [37, 101, 257, 317, 373, 389] | 5 | 147 | 911 (228) | 157 (39) | 163 | [173] | 32/112 | yes |
| bal7 | 7 | [97, 101, 109, 113, 137, 149, 157] | 4 | 92 | 599 (150) | 127 (32) | 131 | [] | 23/66 | yes |
| first7 | 7 | [5, 13, 17, 29, 37, 41, 53] | 4 | 92 | 599 (150) | 173 (43) | 179 | [181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 271, 293] | 42/48 | yes |
| rand7 | 7 | [5, 41, 53, 149, 173, 257, 313] | 4 | 88 | 509 (127) | 103 (26) | 107 | [109, 113, 127, 131, 137, 139, 157] | 28/58 | yes |


Reading (*exact data*):
* For the **tightest** primes the one-modulus segment grows steadily with `K`: `n_seg = 20, 39, 66, 83` (`r = 4`,
  `K = 4..7`), `n_seg = 39, 73, 102` (`K = 5`, `r = 4..6`); `n_seg / |S_h| = 0.49, 0.95, 0.51, 0.64` (`r = 4`) and
  `0.95, 1.20, 1.20` (`K = 5`). In `r` at `K = 4` it grows only up to `r = 6`: `n_seg = 20, 32, 43, 43, 43` for
  `r = 4..8` (`first8`, `K = 4`: `log_def_first8_K4.txt`, run up to `q = 353`: FLEX at `179` and `307`, undecided at `347`
  (level-1 survivors outside the class, none verified when the status-unknown moduli carry genuine labels), every other
  forced `q <= 353` RIGID; the forced set has 97 moduli, `q* >= 631`, 26 moduli above 499 of unknown status).
  So `n_seg/|S_2|` falls from `0.49-0.52` (`r = 4..6`) to `0.38, 0.30` (`r = 7, 8`): the columns added at `r >= 7`
  (`p = 53, 61`) are looser, and the same flexible label at `q = 179` survives. **One-modulus rigidity reaches the
  pigeonhole scale `|S_floor(K/2)|` (i.e. `M^(floor(K/2))`) as long as the added primes stay small**; it is a fraction
  `0.16-0.45` of the forced range `n*`.
* For **loose** primes (`last`, primes near 400) it barely moves: `n_seg = 7, 10, 10, 13` (`r = 4`, `K = 4..7`).
  `rand` and `bal` are in between (`rand4`: `7, 12, 20, 20`; `bal4`: `9, 15, 21, 32`).
* Rigid moduli also occur above `q_flex1` (column "rigid above it"); the segment is not the whole rigid set.

### 3.3 Joint rigid cores: every forced modulus free (*exact on both sides*; `fj.py`, `cpc.c`, `check_system.py`)

**The test.** For `Q` in `F`, put `C = F ∩ [3, Q]` and ask: is there a good system (all forced moduli free, all levels)
whose labels on `C` are **not** system-level coherent with any column-sign conjugate `eps * L` of the genuine labels
(i.e. `AND_(q in C) mask(ell(q)) = 0`)? This is exactly the negation of "`Q <= Q_core`".
* *Certificate of coherence (lower bound).* By Lemma R it suffices that the following relaxation is infeasible:
  moduli `F ∩ [3, Q_sub]` only; level-1 signatures, except at the smallest moduli where labels are taken mod
  `m q^(a-1)` with exact multi-level signatures (`(3,3), (5,3), (7,2)`, then `(11,2), (13,2)`; "depth level" 0/1/2);
  individually admissible labels only; minimal signatures off `C`, all realised (signature, mask) pairs on `C`.
  `cpc.c` decides it by complete DFS with sound propagation (capacity tests with tolerance `+1e-9` in favour of
  feasibility; logs are the only floating-point quantities). Both the C search and an independent Python
  implementation (`cp_mix.py`) return INFEASIBLE with identical node counts on `first4`, `K = 4`.
* *Witness of non-coherence (upper bound).* A solution of the deepest relaxation is converted to explicit labels,
  extended to all of `F` by exact one-modulus DFS (`jointcore.extend`), lifted level by level, and verified exactly
  (`joint_ilp.verify`); then **re-checked from scratch** by `check_system.py` (all `v` with `||v||_1 <= K`, both
  signs, every level, exact integer budgets, coverage of `F_K^ex`, class status of every label).
* Sanity checks: the search finds the known feasible cases (target "outside `class(L)`" at `q0 = 37, 47` — realised by
  conjugate-genuine systems — and the `Q = 67` witness), so it is not vacuously infeasible.

## 4. Non-genuine survivors above the core (explicit witnesses)

**One-modulus witnesses.** Every FLEX entry of §3.2 is an explicit label vector at `q0` (field `witnesses` of
`def_wit_<set>_K<K>.json`: the level-1 label `ell`, its collisions `(v, s)`, and a level-2 lift `lift` with no level-2
collision), completed by genuine labels at all other forced moduli and collision-free vectors elsewhere, and
re-verified with exact integer budgets. Examples (first FLEX modulus; label `ell in (Z/m)^r`):
| set | K | q0 | m | ell (level 1) | #collisions | collisions (v : unit s) | level-2 lift |
|---|---|---|---|---|---|---|---|
| bal4 | 4 | 41 | 40 | [1, 33, 3, 22] | 16 | (-3, 0, 1, 0):0; (-3, 1, 0, 0):3; (-2, 0, 0, 1):2; (-1, -1, 0, 2):1 ... | [681, 993, 1443, 22] |
| bal4 | 5 | 67 | 68 | [1, 3, 43, 61] | 19 | (-3, 0, 0, 2):3; (-3, 1, 0, 0):0; (-2, 0, 1, 1):2; (-1, -3, 0, 1):3 ... | [613, 955, 3579, 4073] |
| bal4 | 6 | 89 | 88 | [6, 27, 24, 43] | 29 | (-4, 0, 1, 0):0; (-3, -1, 1, 1):3; (-3, 0, -1, 2):2; (-2, -1, -2, 1):2 ... | [886, 6803, 2928, 1539] |
| bal4 | 7 | 131 | 132 | [1, 64, 24, 98] | 35 | (-3, -2, 0, 1):3; (-2, -2, 0, 2):2; (-2, 0, -4, 1):0; (-2, 2, 3, 0):2 ... | [10561, 12604, 17184, 3002] |
| bal5 | 4 | 59 | 60 | [1, 11, 39, 51, 32] | 21 | (-2, 0, 0, 0, 1):2; (-2, 1, 0, 1, 0):0; (-1, -1, -1, 1, 0):0; (-1, -1, 0, 2, 0):2 ... | [841, 311, 2499, 1791, 92] |
| bal5 | 5 | 127 | 128 | [1, 23, 121, 18, 113] | 21 | (-3, 0, 0, 1, 1):0; (-2, -1, -1, 1, 0):0; (-2, -1, 1, 0, 0):3; (-2, 0, 0, 0, 2):3 ... | [641, 14871, 10233, 9490, 9329] |
| bal6 | 4 | 103 | 104 | [2, 5, 57, 75, 70, 54] | 25 | (-2, 0, 0, 0, 0, 2):0; (-1, -1, 1, 0, 0, 1):0; (-1, 0, 0, 0, 0, 1):2; (-1, 0, 0, 0, 3, 0):0 ... | [9466, 6037, 9625, 9331, 9846, 2862] |
| bal6 | 5 | 181 | 180 | [3, 126, 29, 7, 12, 19] | 37 | (-4, 0, 0, 0, 1, 0):0; (-2, 0, -2, 0, 0, 1):1; (-1, -1, -2, 1, 0, 0):0; (-1, 0, 0, 0, 4, 0):3 ... | [11343, 3186, 10289, 1087, 2532, 14059] |
| bal7 | 4 | 139 | 140 | [1, 73, 55, 42, 108, 93, 87] | 27 | (-3, 0, 0, 0, 1, 0, 0):3; (-3, 1, 0, 0, 0, 0, 0):2; (-2, 0, 1, 0, 0, 0, 1):0; (-1, -2, 0, 1, 0, 0, 0):1 ... | [1, 19253, 10835, 182, 18868, 12273, 10027] |
| first4 | 4 | 83 | 84 | [1, 75, 46, 48] | 7 | (-2, 0, -1, 1):0; (-1, -2, 1, 0):3; (-1, 1, 1, 1):0; (0, -1, 0, 2):1 ... | [5293, 4947, 970, 1392] |
| first4 | 5 | 163 | 164 | [3, 39, 31, 149] | 8 | (-2, 0, 2, 1):3; (-1, -1, 1, 2):1; (-1, 1, -2, 1):1; (0, -2, 0, 3):3 ... | [22635, 26115, 19219, 24585] |
| first4 | 6 | 269 | 268 | [4, 148, 245, 151] | 6 | (-2, 1, 1, 1):0; (-1, 1, -1, 2):1; (0, -2, 1, 3):2; (1, 0, -2, 1):1 ... | [5364, 3632, 9089, 42763] |
| first4 | 7 | 337 | 336 | [9, 264, 97, 329] | 12 | (-4, 3, 0, 0):3; (-3, -1, 3, 0):0; (-2, 1, 1, 1):0; (-1, 3, -1, 2):0 ... | [34617, 85944, 89473, 62153] |
| first5 | 4 | 131 | 132 | [6, 50, 47, 57, 127] | 9 | (-2, 1, 0, 0, 1):1; (-1, -1, 0, 0, 2):2; (-1, 0, 0, 3, 0):1; (0, 0, -2, 0, 1):1 ... | [10566, 12590, 17207, 2961, 13327] |
| first5 | 5 | 307 | 308 | [11, 17, 134, 93, 90] | 11 | (-3, 0, -1, 0, 1):1; (-3, 1, 0, 1, 0):3; (0, -1, -1, -1, 1):2; (0, -1, 3, 0, 0):3 ... | [80399, 47141, 46642, 82021, 75858] |
| first6 | 4 | 179 | 180 | [2, 18, 48, 62, 103, 173] | 15 | (-1, -2, 0, 0, 0, 1):1; (-1, -1, 1, 1, 0, 0):2; (-1, 0, 0, 2, 1, 0):3; (0, -1, -1, -1, 0, 1):3 ... | [9002, 5598, 2928, 2582, 23863, 14033] |
| first6 | 5 | 419 | 420 | [4, 192, 295, 348, 318, 364] | 15 | (-1, -1, 1, 0, 2, 0):1; (-1, 0, -3, 0, 0, 1):1; (0, -4, 0, 1, 0, 0):0; (0, -2, -1, 0, 0, 1):3 ... | [133984, 146772, 155275, 27228, 28038, 60844] |
| first7 | 4 | 179 | 180 | [2, 18, 48, 62, 103, 173, 145] | 23 | (-1, -2, 0, 0, 0, 1, 0):1; (-1, -1, 0, 0, 0, 0, 2):2; (-1, -1, 1, 1, 0, 0, 0):2; (-1, 0, 0, 2, 1, 0, 0):3 ... | [9002, 5598, 2928, 2582, 23863, 14033, 20485] |
| last4 | 4 | 31 | 32 | [1, 18, 19, 6] | 22 | (-3, 0, 1, 0):2; (-2, 1, 0, 0):2; (-1, -1, -1, 1):0; (-1, -1, 1, 0):0 ... | [385, 146, 691, 38] |
| last4 | 5 | 43 | 44 | [2, 7, 25, 40] | 27 | (-3, 0, 2, 0):0; (-2, -2, 0, 1):2; (-2, 0, -1, 1):3; (-1, -1, 3, 0):2 ... | [90, 799, 421, 1316] |
| last4 | 6 | 43 | 44 | [2, 23, 17, 40] | 59 | (-4, 1, 0, 1):3; (-3, -1, 0, 1):3; (-3, 0, -2, 1):0; (-3, 0, 1, 0):3 ... | [1410, 1211, 1601, 1052] |
| last4 | 7 | 59 | 60 | [1, 20, 20, 5] | 101 | (-5, 0, 0, 1):0; (-5, 0, 1, 0):1; (-5, 1, 0, 0):1; (0, -5, 0, 2):2 ... | [2281, 1640, 860, 2405] |
| last5 | 4 | 41 | 40 | [8, 22, 16, 2, 19] | 38 | (-2, 0, -1, 1, 0):1; (-2, 0, 1, 0, 0):0; (-1, -2, 0, 1, 0):3; (-1, 0, -1, 2, 0):2 ... | [968, 862, 576, 402, 979] |
| last5 | 5 | 59 | 60 | [1, 29, 40, 40, 25] | 62 | (-1, -1, -1, 0, 1):1; (-1, -1, -1, 1, 0):2; (-1, -1, 0, -1, 1):1; (-1, -1, 0, 0, 3):3 ... | [841, 329, 2500, 1780, 85] |
| last6 | 4 | 59 | 60 | [1, 9, 21, 40, 40, 25] | 45 | (-1, -1, 0, 0, 0, 1):1; (-1, -1, 0, 0, 1, 0):2; (-1, -1, 0, 1, 0, 0):2; (-1, 0, 1, 0, 0, 1):3 ... | [841, 309, 2481, 1780, 100, 505] |
| last6 | 5 | 131 | 132 | [2, 5, 62, 121, 51, 36] | 53 | (-3, 0, 0, 0, 0, 2):2; (-2, 1, 1, 0, 0, 1):3; (-1, -1, -1, 0, 0, 1):3; (-1, -1, -1, 0, 2, 0):1 ... | [7394, 14921, 16562, 6853, 7575, 1752] |
| rand4 | 4 | 31 | 32 | [1, 1, 13, 12] | 19 | (-2, 0, 2, 0):3; (-2, 2, 0, 0):0; (-1, -1, 2, 0):3; (-1, 0, 1, 1):3 ... | [385, 129, 685, 44] |
| rand4 | 5 | 59 | 60 | [12, 5, 47, 1] | 20 | (-3, 1, 0, 1):2; (-2, -1, -1, 1):3; (-2, 1, 2, 0):1; (-1, -1, 0, 2):3 ... | [852, 305, 2507, 1741] |
| rand4 | 6 | 83 | 84 | [1, 14, 77, 64] | 33 | (-3, 0, 0, 3):1; (-2, -1, 1, 2):1; (-2, 0, 0, 2):2; (-2, 1, -1, 2):3 ... | [2017, 1862, 3353, 904] |
| rand4 | 7 | 83 | 84 | [1, 14, 77, 64] | 47 | (-3, 0, 0, 3):1; (-2, -3, 0, 2):0; (-2, -2, -1, 2):1; (-2, -1, -2, 2):2 ... | [2017, 1862, 3353, 904] |
| rand5 | 4 | 53 | 52 | [2, 19, 50, 38, 15] | 26 | (-3, 1, 0, 0, 0):3; (-2, 0, -1, 0, 1):3; (-2, 0, 0, 0, 2):2; (-2, 1, 1, 0, 0):3 ... | [2030, 695, 1558, 1702, 2355] |
| rand5 | 5 | 73 | 72 | [8, 24, 32, 11, 24] | 50 | (-4, 0, 1, 0, 0):0; (-3, 0, 0, 0, 1):0; (-3, 1, 0, 0, 0):0; (-2, -1, 0, 2, 0):1 ... | [2528, 1104, 4640, 4403, 5208] |
| rand6 | 4 | 83 | 84 | [12, 19, 5, 8, 74, 43] | 31 | (-2, -1, 0, 0, 0, 1):0; (-2, 0, -1, 1, 0, 0):3; (-2, 0, 0, -1, 1, 0):2; (-2, 1, 1, 0, 0, 0):0 ... | [5304, 4891, 929, 1352, 4106, 379] |
| rand6 | 5 | 163 | 164 | [3, 149, 109, 113, 84, 113] | 46 | (-2, 0, -1, 0, 1, 1):2; (-2, 0, -1, 1, 1, 0):2; (-2, 0, 0, 0, 3, 0):2; (-2, 1, 0, 0, 0, 2):3 ... | [22635, 26225, 19297, 24549, 7300, 20777] |
| rand7 | 4 | 107 | 108 | [1, 64, 40, 46, 67, 63, 62] | 36 | (-2, 1, 0, 1, 0, 0, 0):0; (-1, -1, 0, 2, 0, 0, 0):1; (-1, 0, 0, 1, 0, 1, 0):0; (-1, 1, 0, 0, 0, 2, 0):3 ... | [6589, 10540, 4252, 1882, 5899, 7515, 10430] |


**Joint witnesses** (`jwmix_<set>_K<K>_<Q>.pkl`: labels at every forced modulus and all lifts; each re-checked
from scratch by `check_system.py`, which recomputes every collision of every `v`, `||v||_1 <= K`, both signs, at every
level, checks every budget exactly, and classifies every label):

| witness file | set | K | non-coherent on F ∩ [3,Q], Q = | |F| | coherent (one eps) on forced q <= | moduli with label outside every eps-class | of the forced moduli above Q: outside every eps-class |
|---|---|---|---|---|---|---|---|
| `jwmix_bal4_K4_11.pkl` | bal4 | 4 | 11 | 40 | 7 | 35 (smallest 11) | 34/36 |
| `jwmix_bal4_K5_11.pkl` | bal4 | 5 | 11 | 68 | 7 | 63 (smallest 11) | 62/64 |
| `jwmix_first4_K4_67.pkl` | first4 | 4 | 67 | 40 | 61 | 24 (smallest 67) | 23/26 |
| `jwmix_first4_K5_113.pkl` | first4 | 5 | 113 | 65 | 109 | 38 (smallest 113) | 37/40 |
| `jwmix_first4_K6_197.pkl` | first4 | 6 | 197 | 115 | 193 | 75 (smallest 197) | 74/75 |
| `jwmix_first5_K4_103.pkl` | first5 | 4 | 103 | 54 | 101 | 34 (smallest 103) | 33/33 |
| `jwmix_last4_K4_11.pkl` | last4 | 4 | 11 | 44 | 7 | 38 (smallest 11) | 37/40 |
| `jwmix_last4_K5_11.pkl` | last4 | 5 | 11 | 71 | 7 | 68 (smallest 11) | 67/67 |
| `jwmix_rand4_K4_11.pkl` | rand4 | 4 | 11 | 40 | 7 | 35 (smallest 11) | 34/36 |
| `jwmix_rand4_K5_11.pkl` | rand4 | 5 | 11 | 69 | 7 | 65 (smallest 11) | 64/65 |


Each joint witness is a good system (norm `<= K`, every level, every modulus) which is coherent with one column-sign
conjugate of the genuine labels exactly on the rigid core and then breaks coherence at the first modulus above it;
above the core almost every label lies outside **every** column-sign class. These are the "flexible tops" of the
task, now in the strongest sense (all forced moduli free, all levels): they are the seeds for a C1 construction.

## 5. What makes moduli rigid or flexible

Three nested notions, all exact (`LOCe` range `<=` joint core `<=` one-modulus segment):

1. **Local rigidity (`loc.py`).** Every individually admissible label at `q` lies in a column-sign class. It uses only
   the monomials with `Norm(Z_v) < q^2` and is therefore **independent of `K`** once `p_min^(K+1) >= q^2` (Lemma K).
   Ranges (initial segments; the same for every `K >= 4` tested): `first4` `q <= 31`; `first5`, `first6`, `first7`
   `q <= 11`; `rand4`, `bal4`, `last4` only the vacuous `q <= 7` (`n_q <= 2`). More columns *shrink* it (a label has
   more coordinates to place). **Local rigidity cannot grow with `K`.**
2. **One-modulus rigidity ("spent budgets").** With genuine labels elsewhere, a label at `q0` may only collide on
   `(v, s)` whose slack `B(v,s) - log(genuine demand off q0)` is at least `log q0`. `mech.py` splits the forbidden
   collisions into *intrinsic* (`q0^2 > cap`) and *spent* ones: at `first4`, `K = 4`, the rigid moduli `11..43` are
   already rigid with the constraints of norm `<= 3`, the moduli `47..103` need norm 4 (`k_min` column of
   `mech_first4_K4.json`). Intrinsic constraints alone never give the genuine class itself beyond the vacuous moduli (at
   most the column-sign classes, `q <= 31`); one-modulus rigidity always uses spent budgets. This is round 10's mechanism, and
   it reaches the pigeonhole scale `n ~ |S_floor(K/2)|` (§3.2).
3. **Joint rigidity (the rigid core).** With every modulus free, the spent budgets are no longer fixed; what remains
   is the *forced* demand of the other moduli. `mech_joint.py` (exact; relaxations of the certified instance
   `first4`, `K = 4`, `C = F ∩ [3, 61]`, moduli `<= 79`):
   * removing any one of the core moduli `3, 7, 37, 41, 47, 53` (its demand and its coherence requirement), or the
     demand of the outer modulus `73`, makes a non-coherent system appear in the relaxation; removing any other single
     modulus does not;
   * lifting all norm-3 budgets, or all norm-4 budgets, makes it feasible; norms 1 and 2 are not needed; no single
     norm suffices on its own.

   * the same analysis at `first4`, `K = 5` (`C = F ∩ [3, 109]`, moduli `<= 199`, depth level 1;
     `mech_joint_first4_K5.txt`): needed are the outer moduli `113, 151` and the core moduli
     `3, 7, 11, 19, 23, 31, 37, 43, 61, 67, 71, 73, 79, 83, 107`; lifting the norm-3, norm-4 or norm-5 budgets each
     destroys the certificate; norms 1-2 are not needed; no single norm suffices.

   So the core is held by the **forced collisions of the smallest moduli** (`q = 3`: `n = 1`, every monomial
   collides, `log 3` charged to each `v`; `q = 7`: `n = 2`, every monomial with `<v, c> = 0 mod 2` collides) together
   with a few mid moduli, acting through the norm-3 and norm-4 budgets simultaneously. It is a genuinely joint
   (packing) obstruction: no single modulus and no single norm carries it.

## 6. Growth of the rigid core in `K` (and `r`); comparison with `M^(floor(K/2))` and the forced/free threshold

All `Q_core` entries are exact on both sides unless marked `>=` (lower bound certified, no witness yet).
`n(q) = m_q/4`. `|S_h|` = pigeonhole threshold `~ (2r)^h/h!` (the `M^(floor(K/2))` scale of L4K Theorem 1), `n*` = exact
forced/free boundary, `LOCe` = local rigidity range (§5), `seg`/`q_flex1` = one-modulus segment/first flexible modulus.

| set | r | K | Q_core | n(Q_core) | first non-coherent Q (witness) | LOCe | one-mod seg (n) | q_flex1 | n* | pigeonhole |S_h| | n(Q_core)/|S_h| |
|---|---|---|---|---|---|---|---|---|---|---|---|
| bal4 | 4 | 4 | 7 | 2 | 11 | 7 | 37 (9) | 41 | 68 | 41 | 0.05 |
| bal4 | 4 | 5 | 7 | 2 | 11 | 7 | 61 (15) | 67 | 92 | 41 | 0.05 |
| last4 | 4 | 4 | 7 | 2 | 11 | 7 | 29 (7) | 31 | 66 | 41 | 0.05 |
| last4 | 4 | 5 | 7 | 2 | 11 | 7 | 41 (10) | 43 | 96 | 41 | 0.05 |
| rand4 | 4 | 4 | 7 | 2 | 11 | 7 | 29 (7) | 31 | 50 | 41 | 0.05 |
| rand4 | 4 | 5 | 7 | 2 | 11 | 7 | 47 (12) | 59 | 152 | 41 | 0.05 |
| first4 | 4 | 4 | 61 | 15 | 67 | 31 | 79 (20) | 83 | 50 | 41 | 0.37 |
| first4 | 4 | 5 | 109 | 27 | 113 | 31 | 157 (39) | 163 | 87 | 41 | 0.66 |
| first4 | 4 | 6 | 193 | 48 | 197 | 31 | 263 (66) | 269 | 260 | 129 | 0.37 |
| first4 | 4 | 7 | >= 239 | >= 60 | - | 31 | 331 (83) | 337 | 522 | 129 | >= 0.47 |
| first5 | 5 | 4 | 101 | 25 | 103 | 11 | 127 (32) | 131 | 73 | 61 | 0.41 |
| first5 | 5 | 5 | >= 163 | >= 41 | - | 11 | 293 (73) | 307 | 180 | 61 | >= 0.67 |


**Fits (evidence; tiny `r`, few points).**
* *Tight primes, growth in `K` (`r = 4`).* `n(Q_core) = 15, 27, 48, >= 60` for `K = 4, 5, 6, 7`. Between `K = 4` and `5`
  (same `h = 2`) the core already grows by `x1.8`; at `K = 6` it is `x3.2` the `K = 4` value. Normalised by the
  pigeonhole scale: `n(Q_core)/|S_floor(K/2)| = 0.37, 0.66, 0.37, >= 0.47`; by the forced range:
  `n(Q_core)/n* = 0.30, 0.31, 0.18, >= 0.11`. So the core tracks `|S_floor(K/2)|` (i.e. `M^(floor(K/2))`), not the
  forced range, which grows faster at even `K` (`n*/|S_h| = 1.2, 2.1, 2.0, 4.0`).
* *Tight primes, growth in `r`.* `K = 4`: `n(Q_core) = 15, 25` for `r = 4, 5` (`|S_2| = 41, 61`): ratio `1.67`, against
  `1.49` for `|S_2|` and `1.56` for `r^2` (local exponent `2.3`). `K = 5`: `n(Q_core) = 27, >= 41`: ratio `>= 1.52`.
  For `r = 6, 7, 8` (`K = 4`) only the exact upper bound is known: `Q_core < q_flex1 = 179`, i.e. `n(Q_core) <= 43`,
  against `|S_2| = 85, 113, 145`: the `r`-growth at fixed `K` cannot continue at the `|S_2|` rate once the added primes
  (41, 53, 61) are no longer small.
* *One-modulus segment vs joint core.* `n_seg/n(Q_core) = 1.33, 1.44, 1.38, <= 1.38` (`r = 4`), `1.28` (`r = 5`):
  the joint core is a stable fraction (`0.69-0.78` in `n`, exact cases) of the one-modulus segment, which itself scales like `|S_h|`
  (§3.2: `0.3-1.2 |S_h|`, decreasing in `r` once the added primes are no longer the smallest).
* *Loose primes.* `rand4`, `bal4`, `last4`: `Q_core = 7` (only the vacuous moduli `n <= 2`) at `K = 4, 5`; the first
  non-vacuous forced modulus `q = 11` already carries a non-coherent good system, although one-modulus rigidity reaches
  `q = 29-61`. **No growth at all** from `K = 4` to `K = 5` for these sets.

## 7. Verdict on Rig(K) at this test (evidence only) and open points

Labels: everything in §§3-6 is *exact* for the stated finite instances; what follows is *evidence* about Rig(K) in
regime (1.4), which no finite computation can decide. No growth improvement, no uniform bound and no cap is claimed;
the bracket `3 M log M - O(M) <= W_min^all <= O(M log^2 M/loglog M)` and `M <= (2/3 + eps) log R/loglog R` are
unchanged.

1. **The joint rigid core exists, is computable exactly, and is strictly smaller than the one-modulus segment.** In
   every tight instance it is a non-trivial initial segment (`first4`: 61, 109, 193, `>= 239` for `K = 4..7`;
   `first5`: 101, `>= 163` for `K = 4, 5`), and it is exactly bracketed (except the `>=` entries) by a certificate (Lemma R relaxation + complete search) and a
   verified witness. The round-10 "rigid core `{3..47}` at `K = 4`" survives in the strong sense (everything free):
   the true core is `[3, 61]`.
2. **For the tightest primes the core grows with `K` and with `r`, at the pigeonhole scale.** `n(Q_core)` is
   `0.37-0.66 |S_floor(K/2)|` in all four tight `K`-values and both `r`-values tested; it does *not* stall where it
   could (between `K = 4` and `K = 5`, same `h`), and the one-modulus segment that bounds it from above scales the same
   way, *as long as the added primes are small*: at `K = 4` the upper bound `q_flex1 = 179` is the same for `r = 6, 7, 8`
   (`first6..8`), so growth in `r` at fixed `K` saturates once the new columns are loose. This is the "growing" branch of round 13 §7.2, **for tight primes only**: a first local rigidity mechanism,
   driven by the forced collisions of the smallest moduli (`q = 3, 7`, every monomial collides) acting jointly with the
   norm-3..K budgets (§5).
3. **For loose primes it stalls completely.** With all `p_j >= 53` (`rand4`, `bal4`, `last4`, `r = 4`) the core contains
   only the vacuous moduli `n_q <= 2` at `K = 4, 5` (all six cases exact); the witnesses are good systems whose labels
   at `34/36` to `66/67` of the forced moduli above `q = 11` lie outside **every** column-sign class. These are explicit flexible tops covering
   the entire forced range: C1 seeds at norm `<= K`.
4. **What this means for the two branches of round 13 §7.2.**
   * *Growth branch.* A proof of Rig(K) would have to use the arithmetic of small primes (the data show no rigidity
     without them), i.e. live in the part of regime (1.4) where `log p_j/log r` is close to 1 for a positive share of
     the budget. Near-extremal clusters put their primes there (round 10 L4K §5 item 4), so this is not excluded, but
     a statement for all prime sets in regime (1.4) is not supported: at the tested sizes the loose sets have
     vacuous cores at norms `<= 5` (they could still start growing at larger `K`).
   * *Construction branch (C1).* The loose-prime witnesses are the first explicit label systems that are good for all
     norms `<= K` and non-coherent on essentially all of the forced range. They are **not** P3-all\* members: goodness
     for all norms is needed, and by Lemma M the all-norm core is `>= sup_K Q_core(K)`, which for loose primes was only
     tested up to `K = 5`; at `K = 6` (`rand4`, `last4`) the relaxations still admit non-coherent systems at `q = 11`
     at every depth level tried, but all 100 solutions examined failed to lift because of deep 3-adic collisions
     (`q = 3` at levels 3-5: 215 violations, `q = 5` level 3: 53, `q = 11` level 2: 23; `diag_rand4_K6.txt`), and
     deeper tables exhausted memory, so `K = 6` is undecided for loose primes. Whether the loose cores stay vacuous for all `K` is the
     next finite question; if they start to grow at `K >= 6` the loose sets behave like the tight ones with a delay.
5. **Comparison with the targets of round 13 §7.2.** `M^(floor(K/2))`: the tight-prime core tracks it (ratio
   `0.37-0.66`). Forced/free threshold of L4K Thm 1: the core is a decreasing fraction of it at even `K`
   (`0.30, 0.31, 0.18, >= 0.11`), because `n*` itself grows faster than `|S_h|` at even `K`.

**Open points.**
* `first4`, `K = 7`: exact core. Certified `>= 239` (moduli `<= 263`); at `Q = 241` the relaxations (all depth levels)
  contain non-coherent systems, but every one of 30 solutions fails to extend at the first modulus above the window
  (`q = 269`, the one-modulus FLEX modulus): the window is too small to decide. (`K = 6` was decided only after
  making the lift search exhaustive for `q^r <= 2*10^5`: the earlier failures were level-2 collisions at `q = 19`.)
* `first5`, `K = 5`: certified `>= 163` (moduli `<= 163`); no witness attempted above.
* `r = 6..8` joint cores: the signature-table enumeration does not scale to `r = 6` (the `first6`, `K = 4` tables up to
  `q = 107` did not finish in 40 minutes); `r = 8` forced sets (`first8`, `K = 4`) took 4.8 hours with 26 moduli above 499 left
  undecided at the node cap; `rand8`, `bal8` were not reached. One-modulus data exist for `r = 4..7`.
* Loose sets at `K = 6, 7` (the decisive test of item 4): needs `q = 3` modelled to 3-adic depth `>= 5` with compact
  (C-side) dominance reduction.
* All-norm version: a joint test with budgets of *all* norms `<= K_max` imposed on `F_(K_max)` for several `K_max` (this
  is what P3-all\* imposes); by Lemma M it can only enlarge the cores reported here.
* A proof: the mechanism of §5 suggests an exact statement of the form "the forced demand of the moduli `n_q <= 2`
  plus the pigeonhole demand of the moduli up to `Q` exceeds the slack of the norm-`3..K` budgets unless the labels are
  coherent", whose scale is `|S_floor(K/2)|`. Nothing of this is proved here.

## Appendix: scripts and outputs (`rigcompute_checks/`)

| file | role |
|---|---|
| `rl.py` | shared exact helpers (labels, parity, `Z_v`, `D_s`, exact caps `Norm(Z_v)(2 if s odd)`, torus discrete logs, orbit representatives, `./csp` driver) |
| `csp.c` | single-modulus exhaustive DFS (count / witness / existence / signature-table modes) |
| `freescan.c`, `thresholds.py`, `run_thr.py`, `bh.py`, `sets.py` | exact forced sets `F_K^ex` (`thr_<set>_K<K>.json`), Theorem 1(b) scan bound, prime sets |
| `deform14.py`, `run_def.sh`, `summarize.py` | one-modulus tests (`def_wit_<set>_K<K>.json`, logs `log_def_*`), table `onemod_table.md` |
| `mech.py` | intrinsic vs spent forbidden collisions, minimal norm forcing one-modulus rigidity (`mech_<set>_K<K>.json`) |
| `loc.py` | local rigidity LOC/LOCe (`loc_<set>_K<K>.json`) |
| `sigtab.py`, `deepsig.c`, `deep.py` | level-1 and multi-level signature tables |
| `cp_joint.py`, `cp_mix.py` | Python joint DFS (targets: label outside all eps-classes; non-coherence on a segment) |
| `cpc.c`, `fj.py`, `fjs.py` (streaming, memory-light), `ladder.py`, `queue_ladder.sh` | C joint DFS and the joint-core pipeline (`fj_<set>_K<K>_<Qsub>.log`, `ladder_<set>_K<K>.log`, `jointcore_*.json`, caches `l1cache_*`, `deepcache_*`) |
| `jointcore.py` (`extend`), `joint_ilp.py` (`verify`) | witness extension to all forced moduli; exact all-level verification |
| `check_system.py`, `witness_table.py` | independent re-check of every joint witness `jwmix_<set>_K<K>_<Q>.pkl` |
| `mech_joint.py` | relaxation analysis of a certified core (`mech_joint_first4_K4.txt`, `mech_joint_first4_K5.txt`) |
| `joint_lp.py`, `joint_ilp.py`, `ls_joint.py` | LP/Lagrangian certificates, HiGHS MILP, local search (superseded; evidence only) |
| `diag_wit.py`, `diag_stream.py`, `repair.py` | diagnostics of failed witness constructions (`diag_*.txt`); norm-raising repair (failed; evidence only) |
| `joint_results.json` | the joint-core results with their sources (logs, witness files) |
| `joint_results.json`, `make_joint_table.py`, `thr_table.py`, `wit_examples.py` | result tables used in §§3-6 |
| `orig/` (one level up) | read-only copy of the round-10 harness used for reference |

Reproduction: `gcc -O2 -o csp csp.c && gcc -O2 -o freescan freescan.c && gcc -O2 -o deepsig deepsig.c && gcc -O2 -o cpc cpc.c`;
then e.g. `python3 thresholds.py first4 4 5,13,17,29 6549`, `python3 deform14.py thr_first4_K4.json wit all`,
`python3 fj.py first4 4 103`, `python3 check_system.py jwmix_first4_K4_67.pkl thr_first4_K4.json`.
