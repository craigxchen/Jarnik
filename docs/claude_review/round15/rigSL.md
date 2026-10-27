# Round 15 — Rig^SL: the single-angle-vector rigidity core

(skeleton; sections filled progressively)

## 0. Summary of results

PENDING

## 1. Setting, definitions, class conventions

Labels used throughout: *proved* (all quantifiers, proof here); *exact* (finite computation; infeasibility decided by a
complete search whose only floating-point step is a comparison of sums of logarithms with a tolerance `1e-9` **in
favour of feasibility**, so an INFEASIBLE answer is a certificate; every witness re-checked by the independent checker
`slverify.py` in exact integer arithmetic or in 60-digit interval arithmetic); *evidence*. Nothing here is refereed.

**Conventions** (those of `round10/L4K.md` §1, `round14/rigcompute.md` §1). Columns `j = 1..r` carry split primes
`p_j`, `pi_j = a + bi` (`a` odd, `b` even, both positive), `Z_v = prod pi_j^(v_j)` (negative exponents meaning
`conj pi_j`), `N_v = Norm Z_v = e^(w(v))`, `D_s(Z_v) in {Y, X - Y, X, X + Y}` for `Z_v = X + iY`.
`V_K = {v in Z^r : 1 <= ||v||_1 <= K}`. At an odd prime `q` not dividing `N`: `m_q = q - chi_(-4)(q)`, `n_q = m_q/4`,
labels `ell(q^a) in (Z/m_(q^a))^r`, compatible under reduction, **parity-respecting** (`ell_j = [(p_j/q) = -1]
mod 2`, `sharp.md` Cor. 2.2). `v` collides at `q^a` with unit `s` iff `<v, ell(q^a)> = s iota_(q^a) mod m_(q^a)`;
`m_(v,s) = prod_q q^(a_(q,s)(v))`. Prime sets (round 14, `sets.py`): `first4 = {5,13,17,29}`,
`rand4 = {53,61,89,173}`, `bal4 = {97,101,109,113}`, `last4 = {353,373,389,397}`, `firstr` = the `r` smallest split
primes. `F_K^ex` = the exact forced set of round 14 (`thr_<set>_K<K>.json`; referee-checked for `K = 4, 5`).

**The three classes** (round 14 README §1):
* **PU**: parity-respecting labels at all levels with `m_(v,s)^2 <= N_v 2^[s odd]` for all `v in V_K`, all `s`.
* **SL = Rig^SL(K)**: PU plus ONE `phi in R^r` with, for all `v in V_K` and all `s in Z/4`,

  ```text
  |sin(<v,phi> - s pi/4)|  >=  nu_s m_(v,s) e^(-w(v)/2),      nu_0 = nu_2 = 1, nu_1 = nu_3 = 2^(-1/2).       (1.2)
  ```
  (PU is implied by (1.2).) Equivalently `log m_(v,s) <= B_phi(v,s) := (1/2) log(N_v 2^[s odd]) + log|sin(<v,phi> - s pi/4)|`.
* **GB(K)** ("genuine-angle budgets"): labels with `m_(v,s) <= |D_s(Z_v)|` for all `v in V_K`, `s`. This is SL with
  `phi = arg pi := (arg pi_1, ..., arg pi_r)` (since `|D_s(Z_v)| = nu_s^(-1) e^(w(v)/2) |sin(arg Z_v - s pi/4)|`,
  `sharp.md` Prop. 1.2), an exact integer test.

Every class is closed under replacing labels at free moduli by collision-free ones (round 14 Lemma R(1); demands only
decrease, `phi` unchanged), so all systems below live on `F_K^ex` with collision-free labels elsewhere.

**Symmetry group** `G = {+-1}^r x (Z/2)^r` (column conjugations `eps`, associate twists `t`; referee_rigtheory Prop. T)
together with scalings by units `u = 1 mod 4`. On SL systems: conjugation of column `j` maps
`(ell_j, phi_j) -> (-ell_j, -phi_j)`; the twist maps `(ell_j, phi_j) -> (ell_j + m/2, phi_j + pi/2)` at every level.
Both map SL systems to SL systems (Lemma 2.1).

**Coherence and the core.** As in round 14 (`rigcompute.md` §1, `sigtab.class_masks`): a label at `q` lies in the
column-sign class `eps` iff `ell mod n_q = u (eps∘L(q)) mod n_q` for a unit `u` of `Z/n_q`, `L` the genuine labels.
A system is *coherent on a set `C`* if one `eps` serves every `q in C`. Since `m/2 = 0 mod n` and every `u = 1 mod 4`
is a unit mod `n`, this notion is invariant under the whole group (twists, conjugations, scalings), and it is the
weakest of the coherence notions in use (equality mod `m` up to `G` implies it). **`Q_core^X(P, K)`** (class
`X in {PU, SL, GB}`) = the largest forced `Q` such that every `X`-system is coherent on `F_K^ex ∩ [3, Q]`.
Non-coherence in this weak sense is the strongest flexibility statement; a certificate of coherence in this sense is
what round 14 certified in PU.

Obvious relations: `SL ⊂ PU` (so `Q_core^SL >= Q_core^PU`, round-14 certificates transfer); `GB ⊂ SL` (so a GB
witness is an SL witness, with the verified angle vector `phi = arg pi`).

## 2. Angle pinning: when SL collapses to GB

**Lemma 2.1 (fundamental domain; proved).** Every SL system is `G`-equivalent to an SL system with
`phi in [0, pi/4]^r`, and `G` preserves coherence on every set `C` (in the sense of §1).

*Proof.* (i) `phi_j -> phi_j + pi` changes `<v, phi>` by a multiple of `pi`, which leaves `|sin(<v,phi> - s pi/4)|`
unchanged; labels unchanged. (ii) The twist `t = e_j` maps labels to `ell + (m_(q^a)/2) e_j` at every level (compatible
with reduction since `m_(q^a)/2 = m_(q^b)/2 * q^(a-b)` and `q` is odd, so the reduction of `m_(q^a)/2` is `m_(q^b)/2`
mod `m_(q^b)`; parity kept since `4 | m`). Then `<v, ell'> = <v, ell> + v_j m/2`, so `m'_(v,s) = m_(v, s - 2 v_j)`;
and `|sin(<v,phi'> - s pi/4)| = |sin(<v,phi> + v_j pi/2 - s pi/4)| = |sin(<v,phi> - (s - 2v_j) pi/4)|`, with
`nu_s = nu_(s - 2v_j)`. So (1.2) for `(v, s)` in the new system is (1.2) for `(v, s - 2v_j)` in the old one.
(iii) Conjugation `eps`: `m'_(v,s) = m_(eps v, s)` and `<v, eps phi> = <eps v, phi>`, `w(eps v) = w(v)`, `V_K` is
`eps`-stable. By (i)-(iii) `phi_j` may be moved into `[-pi/4, pi/4]` (mod `pi/2`) and then into `[0, pi/4]`.
Coherence: `m/2 = 2n = 0 mod n`, so twists do not change `ell mod n`; conjugation maps the class `eps'` to
`eps∘eps'`; so "one `eps` for all `q in C`" is preserved. QED

**Lemma 2.2 (tight label-free constraints; proved).** Every system has `m_(v,s) >= 1`, so (1.2) implies the
label-free constraints `(*)`: `|sin(<v,phi> - s pi/4)| >= nu_s e^(-w(v)/2)` for all `v in V_K`, `s`. At `phi = arg pi`
the constraint `(*)` for `(v, s)` holds with equality iff `|D_s(Z_v)| = 1`. In that case `(*)` excludes the open interval
of length `2 alpha`, `alpha = arcsin(nu_s e^(-w(v)/2))`, centred at the nearest point of `s pi/4 + pi Z`, and
`<v, arg pi>` is an endpoint of it; so near `arg pi` the constraint is the half-space `sigma <v, phi - arg pi> >= 0`,
`sigma = sign(<v, arg pi> - s pi/4 - k pi)`.

*Proof.* `|D_s(Z_v)| = nu_s^(-1) e^(w/2) |sin(arg Z_v - s pi/4)|` and `D_s in Z \ {0}`, so `|D_s| >= 1` with equality iff
`(*)` is tight. The rest is the shape of `{x : |sin x| < c}` (`c < 1`), a union of open intervals of length
`2 arcsin c` centred at `pi Z`. QED

**Theorem 2.3 (angle pinning; proved, computer-assisted).** Let `P = first4 = {5,13,17,29}` and `K >= 4`. Every `phi`
satisfying `(*)` for all `v in V_K` is `G`-equivalent to `arg pi`. Consequently every SL system for `(first4, K)` is
`G`-equivalent to a GB system, and **`Q_core^SL(first4, K) = Q_core^GB(first4, K)`** for every `K >= 4`.

*Proof* (`pin.py first4 4`, output `log_pin.txt`). By Lemma 2.1 take `phi in [0, pi/4]^r`; let `g` be the folded
genuine point (`g = (0.46365, 0.58800, 0.24498, 0.38051)`, the image of `arg pi` under `G`; the folded primes are the
associates/conjugates `pi'_j` with `arg pi'_j = g_j mod pi`, and the folded system is GB for the `pi'_j`).
1. *Global exclusion.* A box branch-and-bound over `[0, pi/4]^4` (bisection, rigorous upper bounds of `|sin|` on each box,
   margin `1e-9` in favour of feasibility) discards every box on which some `(v, s)`, `v in V_4`, violates `(*)`; it
   terminates at depth 67 with **all** of `[0, pi/4]^4` discarded except one box contained in the cube
   `|phi - g|_inf <= rho = 10^-6`.
2. *Local cone.* At `g` there are 19 tight pairs (`|D_s(Z'_v)| = 1`, Lemma 2.2). Thirteen of them carry an exact
   positive relation `sum_k lambda_k sigma_k v_k = 0` with `lambda_k in {1, 3/2}` (printed by `pin.py`), and these 13
   vectors span `R^4` (exact rank). Hence a direction `d` with `sigma_k <v_k, d> >= 0` for all `k` has
   `sum lambda_k sigma_k <v_k, d> = 0` with non-negative terms, so all 13 terms vanish and `d = 0`.
3. *The cube.* The smallest tight half-width is `alpha_min = 0.00418`, and `K rho = 4*10^-6 < 2 alpha_min`. So inside the
   cube each tight constraint is exactly its half-space (moving to the wrong side by `|<v,d>| <= K rho < 2 alpha`
   stays inside the excluded interval). By 2, `phi = g`.

Monotonicity in `K` is clear (`V_4 ⊂ V_K`). At `phi = arg pi` (after undoing the group element, which maps
`(v, s)` constraints to `(v', s')` constraints bijectively, Lemma 2.1), (1.2) reads `m_(v,s) <= |D_s(Z_v)|` (§1). QED

*Remarks.* (a) The first column is pinned already at norm 1: for `pi_1 = 1 + 2i`, `D = (2, -1, 1, 3)`, the two tight
constraints `s = 1, 2` pin `phi_1` from both sides (J0 box: `|cos y| <= 3/5` and `|sin y| <= 4/5` force
`|sin y| = 4/5`). At `K = 1` the set `{5,13,17,29}` is not pinned (4 tight pairs; the cone is not trivial). (b) The
pinning mechanism is the existence of Gaussian integers `Z_v`, `||v||_1 <= K`, lying as close to an axis or a diagonal
as an integer point can (`|D_s| = 1`). It is a property of the primes, and it is the angle analogue of the norm-1
pinning of `c1construct` Prop. P (there for labels).

## 3. Loose prime sets (last4, rand4, bal4) at K = 4,5,6 in SL

**Answer: the loose primes still stall once the single angle vector is imposed.** For every loose instance decided,
there is an SL system with the verified angle vector `phi = arg pi` (the genuine angles) whose labels are
non-coherent already at `q = 11`, the first forced modulus with `n_q > 2`; hence `Q_core^SL = 7`, the vacuous value
(coherence on `F ∩ [3, 7]` is automatic because `n_3 = n_5 = n_7 = 1, 1, 2`, so every label is in every class).

**Theorem 3.1 (exact; SL witnesses with verified angle vector).** For the instances in the table there is a
parity-respecting label system on the exact forced set `F_K^ex`, at all levels (collision-free labels at every free
modulus), which satisfies (1.2) for every `v in V_K` and every unit `s` **at `phi = arg pi`** (i.e. `m_(v,s) <= |D_s(Z_v)|`,
checked in exact integers), and whose label at `q = 11` lies outside every column-sign class. So these are GB, hence
SL, hence P3-all\*-truncation systems that are not coherent on `F ∩ [3, 11]`, and `Q_core^SL(P, K) = 7`.

| set | `K` | `|F_K^ex|` | witness file | built by | coherent on | labels outside **every** column-sign class | in the genuine class |
|---|---|---|---|---|---|---|---|
| rand4 `{53,61,89,173}` | 4 | 40 | `W_rand4_K4_gen.pkl` | cpc + deep(3:3,5:3,7:2) + extend + lift | `q <= 7` only | **34/40** | 5/40 |
| bal4 `{97,101,109,113}` | 4 | 40 | `W_bal4_K4_gen_Q11.pkl` | same | `q <= 7` only | **35/40** | 3/40 |
| last4 `{353,373,389,397}` | 4 | 44 | `W_last4_K4_gen_Q11.pkl` | same | `q <= 7` only | **41/44** | 3/44 |
| last4 | 5 | 71 | `W_last4_K5_gen_Q11.pkl` | same | `q <= 7` only | **68/71** | 3/71 |

Each witness was re-checked from scratch by `slverify.py` (own Gaussian arithmetic and discrete logarithms; nothing
imported from round 14): generator check, parity, coverage of `F_K^ex`, absence of collisions at the deepest stored
level (so all levels are determined), all `4 |V_K|` pairs `(v, s)` with both signs of `v` (1280 at `K = 4`, 2720 at
`K = 5`), and the class status of every label. `sanity_verify.py`: the checker accepts the stored witnesses and
rejects 40/40 random parity-preserving perturbations at `q < 60` (first4) and genuine level-1 labels whose collisions
persist at the stored depth.

*How this relates to round 14.* The PU witnesses of rigcompute all fail SL (referee: 8 fail ANG, the two `last4` ones
have no `phi`). The witnesses here are new systems found directly in GB (budgets `log|D_s(Z_v)|` fed to the same
complete joint search `cpc.c`). At `phi = arg pi` the level-1 relaxation on `F ∩ [3, 47]` is non-coherent at `Q = 11`
after 14-24 search nodes for all three sets (`t1.py`, `drv_gen.py`). The one-modulus SL rigidity of rand4 at 31 and
bal4 at 41 (referee, genuine labels elsewhere) is not in conflict: the joint witnesses change almost every label.

## 4. first4: the SL core at K = 4,5,6

By Theorem 2.3, for `first4` and `K >= 3` the SL core is the GB core: a finite problem with **integer** budgets
`m_(v,s) <= |D_s(Z_v)|`. It is decided exactly by the round-14 machinery with these budgets: Lemma R (free moduli,
level-1 relaxation, minimal signatures, individual admissibility `q <= |D_s(Z_v)|`, dropping moduli) holds verbatim
for GB because GB is a per-pair budget class (only the budget numbers change). The complete search `cpc.c` is run on
GB tables (individually GB-admissible labels, `drv_ladder.py`); INFEASIBLE is a certificate (the only float step is the
comparison of log sums with tolerance `1e-9` in favour of feasibility; genuine loads equal budgets up to `9e-16`,
`sanity1.py`); witnesses are explicit, all-level, and re-checked by `slverify.py`.

**Theorem 4.1 (exact).** In the class SL (= P3-all\* truncated at norm `K`, labels and one angle vector):

| `K` | `Q_core^PU` (round 14) | **`Q_core^SL(first4, K)`** | certificate (coherent on `F ∩ [3,Q]`) | witness (verified, non-coherent) |
|---|---|---|---|---|
| 4 | 61 | **79** | `Q = 79` on `F_4^ex ∩ [3, 97]` (cpc: root; independent re-check `indep_cert.py`: root) | `W_first4_K4_gen_Q83.pkl`: coherent on `q <= 79`, non-coherent at 83; 20/40 labels outside every class |
| 5 | 109 | **139** | `Q = 139` on all of `F_5^ex` (65 moduli, `<= 349`; root) | `W_first4_K5_gen_Q149.pkl`: coherent on `q <= 139`, non-coherent at 149; 33/65 labels outside every class |
| 6 | 193 | PENDING | | |

So the angle coupling moves the core from 61 to 79 at `K = 4` and from 109 to 139 at `K = 5`, and the SL core grows
strictly `79 -> 139` from `K = 4` to `K = 5` (round 14 had only `[61, 103]` and `>= 109`). Both witnesses carry the
explicit angle vector `phi = arg pi` and satisfy (1.2) exactly (integer test). The referee's one-modulus statement "83
is rigid under ANG" (genuine labels at all other moduli) is consistent with the joint witness at 83, which changes 20
labels.

*Cross-checks.* (i) `qsubscan.py`: with the moduli `F ∩ [3, Qsub]` only, the level-1 GB relaxation at `Q = 79` is
feasible for `Qsub = 79, 83, 89` and infeasible for `Qsub = 97`: the certificate needs the moduli 83, 89, 97 above the
core. `indep_cert.py` (own enumeration of all parity labels in `(Z/m)^4`, own discrete logs, own class masks, own
propagation) reproduces both answers (non-coherent level-1 system at `Qsub = 79`; certificate at `Qsub = 97`).
(ii) `sanity1.py`: without the non-coherence requirement the same search is feasible, and the genuine level-1
signature is present in every table and fits the budgets.

## 5. r-saturation r = 6..8 in SL

PENDING

## 6. Mechanism

**6.1 Where the angle coupling bites: unit-D monomials pin `phi` (proved for first4, first5, first6).** The label-free
part of (1.2) (`m_(v,s) >= 1`) is tight at the genuine angles exactly on the pairs with `|D_s(Z_v)| = 1`
(Lemma 2.2): Gaussian integers `Z_v` of norm `<= K` that are as close to an axis or diagonal as an integer point can be.
When these pairs positively span `R^r` (and no other component of the feasible region exists), the angle vector is
pinned to `arg pi` up to `G`, and SL collapses to GB, whose budgets `|D_s(Z_v)|` are those of the genuine circle. Counts
of unit-D pairs (sign-normalised `v`, `K = 3, 4, 5, 6`): `first4` 16, 19, 21, 21; `rand4` 2; `bal4` 2; `last4` 0
(`pin.py`, `log_pin.txt`). Pinning needs at least `r + 1` such pairs; small primes have many (for `5 = |1 + 2i|^2` alone,
`D = (2, -1, 1, 3)`, two of them at norm 1), large primes essentially none (heuristically the number of `v in V_K`
with `|D_s(Z_v)| = 1` is `~ sum_v N_v^(-1/2)`). So **the angle coupling is decisive exactly for tight primes, and
there it is total**: no residual angle freedom remains.

**6.2 What the coupling buys for tight primes.** With `phi` pinned, the budget of `(v, s)` drops from the PU value
`(N_v 2^[s odd])^(1/2)` to `|D_s(Z_v)|`, which is smaller by the factor `|sin(arg Z_v - s pi/4)|^(-1)`: on average
over the four units by `~2` (ANG alone gives `m_0^2 + m_2^2 <= N`), and by a large factor for the units at which
`Z_v` is near an axis. This is what moves the core from 61 to 79 (`K = 4`) and from 109 to 139 (`K = 5`).
At `K = 4` (`indep_cert.py first4 4 79 97 base v`): GB-individual admissibility leaves exactly **one** (signature,
class) pair at every forced `q <= 41` (3, 7, 11, 19, 23, 31, 37, 41), always inside a column-sign class containing the
genuine one (masks `0b1`, `0b101`, `0b10001`; 3 and 7 vacuous). From `q = 43` on, labels outside every class are
individually admissible (43: 2 pairs, 53: 3, 59: 8, 61: 7, 67: 36, 71: 31, 73: 17, 79: 134), and they are excluded
only jointly: the mandatory loads of the moduli `<= 97` (above the core: 83, 89, 97 are needed, `qsubscan.py`) leave no
budget for them. At `K = 5` GB admissibility leaves a single signature at every forced `q <= 53` and at 61 (GB table
sizes in `log_lad_first4_K5.txt`). In every certified case the complete search closes **at the root** (pure
propagation, no branching): rigidity in SL for tight primes is "local rigidity at small `q` + joint budget exhaustion
by the forced collisions of the next moduli".

**6.3 How loose primes stay flexible: re-routing of the genuine budget (exact data, `diag.py`).** For loose primes
the budgets `|D_s(Z_v)|` are large (up to `~10^5` at `K = 4`) and mostly smooth over the forced set: the genuine
system already spends almost all of them. The witnesses do **not** live on uncharged slack: they *move* collisions.
In `W_rand4_K4_gen.pkl`, 391 of the 522 charged pairs carry a non-genuine collision; for 180 of them the uncharged
cofactor `U(D_s)` is trivial or tiny and the deviation is paid by **dropping** a genuine collision at a larger forced
modulus. Example: `v = (-1,0,0,1)`, `s = 1`, `D_s = 107` (prime, forced): the genuine system pays 107, the witness
pays 11 (its non-coherent label collides there) and does not collide at 107. Non-genuine collisions are spread over
all forced moduli (34, 19, 61, 44 at `q = 3, 5, 7, 11`, decreasing to 1-3 above 150), on norms 2, 3, 4
(28, 119, 244 pairs). So the flexibility is a *permutation of the collision budget among moduli*, possible because
`log|D_s| >> log q` for most `(v, s)` once all `p_j >= 53`. Rigidity for tight primes is the absence of this room:
many `|D_s(Z_v)|` are small primes or 1.

## 7. Witnesses (explicit phi, verification)

PENDING

## 8. Open points

PENDING
