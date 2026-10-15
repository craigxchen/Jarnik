# Referee report on round14/pal11cert.md (adversarial)

## 0. Verdict

**The main result survives, and I re-proved it by an independent and simpler exact computation.** Theorem 4.3
holds: no 8 nonzero real vectors `z'_1..z'_8 in R^11` satisfy the 1519 width counts of pal11 (`T = 1`, `L = 2`,
`m = 8`, `|A| = 12093`). So `J(A)_8^perp` contains no product of 8 real linear forms.

All of the following re-derive exactly with my own code:
* the box;
* the 1519 test directions and their widths;
* the min-`t` table (Theorem 4.1);
* the attainable types (Lemma 4.2);
* the 25 profiles;
* the automorphism group (Lemma 1.1);
* the twin coordinate spaces (Prop 6.2).

My proof of Theorem 4.3 uses a different and simpler reduction than the author's (§5). Each profile is covered
by a 7-element explicit "rest" whose every count-feasible instance has saturated rank 11. That leaves the 8th
generator completely free, so no `(12,2)` representatives, no symmetry and no knapsack weights are needed.

The round-13 correction is **confirmed** on the two counts that change the record:
* M8' (syl8, `L = 3`) at `m = 10` has 1092 count-feasible pool products, all zonotopes, so `zeta(M8') = 10`;
* M8 at `m = 8` has 42 count-feasible products and 0 certificates (rigorous);
* pal7 at `m = 7` has 2134 count-feasible products and exactly 300 certificates, all zonotopes (rigorous).

**No fatal or major error was found.** The issues are:
* **Scope (not an error, but decisive for this round).** The write-up answers the paused route-5 item, not
  Rig(K) / P3-all\*. None of the P3-all\* class checks apply. It has no bearing on growth, which the author says.
* **Minor:**
  * a false sentence on entry sizes (`(7,5)` generators have entries `3`);
  * fatpoint Prop 2.6(c) is cited for real `z_i`, but it is stated for integer `z_i` (the proof extends verbatim);
  * the M8' model is not named, and the example fails in pal7 coordinates;
  * "436 octics" is a mod-`p` value;
  * "maximal certified" in Prop 6.2 needs qualifying.

## 1. Scope and placement against the round task (Rig(K)) and the barrier logic

* **The write-up does not address this round's question.** Round 14 asks Rig(K) for the residue class P3-all\*
  in regime (1.4). pal11cert.md is the cheap route-5 item of README-13 §7.1 (product certificates for the
  `T = 1` box pal11, `L = 2`, `m = 8`). The round task says route 5 is paused. So the checks the task lists for
  P3-all\* claims have **no object here**: parity, the per-unit budget `m_(v,s) <= nu_s^(-1) e^(w(v)/2)`, which
  monomials, moduli `q <= Q_K`, primes `p_j <= M^(h_0)`, `A` chosen after `c_0`. No statement in pal11cert.md is
  about a label system, a modulus or a prime, and none claims to be. The author says this in §0 and §7 and
  recommends going to §7.2 (Rig(K)). I agree. The write-up is **out of scope for this round**. That is a
  placement issue, not an error.
* **Growth / uniformity.** No growth or uniformity claim is made. None follows: every statement is about one
  fixed box (`n = 12`). Box statements are identical for fakes and actual clusters with this pattern set. The
  author's own §7.2 says the pal11 mechanism (core budget) does **not** extrapolate to the regime (fatpoint
  Thm 3.3: factor 1.15 at `n = 2^12`). I confirm: `M <= (2/3 + eps) log R/loglog R` is untouched.
* **Barrier logic.** Theorem 4.3 is a statement about the width (`b = 1`) system. Fakes satisfy every width
  count by construction, since the box is the same. So the theorem neither caps nor uncaps anything. Its only
  role is to show that the *product-certificate device* cannot decide `TC(pal11, 8)`. The recommendation "do not
  reopen route 5" is a judgement, not a theorem. It is reasonable, but evidence-based (§8 below).

## 2. Re-derivation of the pal11 box (A, tau_A, test directions, symmetry)

All scripts are mine, in `referee_pal11cert_checks/`. Nothing was imported from `pal11cert_checks/` or from the
repository; the Paley-I matrix and `pal_patterns` were re-implemented from their definitions.

* `r_box.py` (exact). 23 patterns: 11 core (size 6) and 12 twins (size 5). `|A| = 12093`. Every slab is active.
  `k_0 = -sum k'` ranges over `[-2, 2]`. The 12 Hadamard columns (all-ones plus `beta_c = 2*1_P - 1`) satisfy
  `H H^T = 12 I`. The test set: `{-1,0,1}^11` and primitive `{-2..2}^11` vectors of support `<= 4` with a `±2`, up
  to sign, give **129,713** directions. Exactly **1519** have `tau_A < 8`, split `11 / 134 / 1374` at
  `tau = 2 / 4 / 6`. All as claimed.
* `r_aut.py` (exact). Backtracking over all signed permutations of the 11 `k'` coordinates, pruned by the pair
  statistics `#{k : k_x = s, k_y = t}` (an invariant), and every complete survivor tested by `g(A) = A`. Result:
  **10** automorphisms, the multipliers `x -> a x` by the five squares, times `±1`. **Lemma 1.1 is confirmed.**
  Remark: the lemma speaks only of signed coordinate permutations, which is what the searches use. Point
  permutations that move the base are not automorphisms of a `T = 1` box anyway, since `k_0` is not in `[-1,1]`.
  The task's "`Z_11` symmetry" does belong to Paley-12, not to pal11. The correction is right.
* Minor: §1.1 says "436 tangent octics". `dim T_8 = 436` is proved only mod `p` (Prop 3.1 says so correctly).
  Over `Q` only `<= 436` is proved. Harmless, since nothing uses the value.

## 3. Theorem 4.1 (min core columns by support): independent check

**Method (`r_flats.c`, exact).** This is a brute force, not a flat enumeration. For every support `R` (`|R| >= 2`,
4083 sets) and **every** subset `Zc` of the 11 core columns (2048), I take `M = {1|_R} ∪ {beta_c|_R : c in Zc}` and
compute its RREF mod `p = 2^31 - 1`. A generic element of `null(M) ⊂ R^R` has full support iff no `e_x` lies in
the row space. Its `t` is `#{c : beta_c|_R notin rowspace(M)}`. `Zc` is recorded only if it is closed, i.e.
`|Zc| = 11 - t`.

**Why mod `p` is exact here.** Every matrix involved has entries in `{0, ±1}` and at most 13 rows. Every minor is
at most `13^(13/2) < 1.8e7 < p` in absolute value. So rank mod `p` equals rank over `Q` for every test used (row
space membership is a rank equality).

**Completeness.** A real `z` with `supp z = R` has a closed zero set `W = Z(z)`. Then `z in null(M_W)` with full
support, and a generic point of `null(M_W)` has zero set `W`. So the brute force attains exactly the `(sigma, t)`
of all real balanced `z`.

**Results** (`r_flats_sum.py`):
* `min t` by `sigma = 2..12`: `6, 9, 6, 7, 2, 5, 4, 3, 4, 4, 1`. **Identical to the table of Theorem 4.1(a).**
* Theorem 4.1(b): each core column `1_P` has `tau_A = 4` and each `e_x` has `tau_A = 2` (`r_box.py`). The
  all-ones direction has `tau_A = 4`. Necessity of the counts for real `z_i` is fatpoint Prop 2.6(c). That
  proposition is stated for `z_i in Z^d`, but its proof (restriction to `u_0 + t v`, top component) is
  verbatim valid for real `z_i`. The write-up should say so. The arithmetic `8 * 6 = 48 > 44` and the
  conclusion "some `t_i <= 5`, hence `sigma_i >= 6`" are correct.
* (c) follows. **Theorem 4.1 stands as stated.**

## 4. Lemma 4.2 and the 25 profiles: independent check

* **Attainable types** (same brute force, all `sigma`). For `sigma <= 8`: `{6}, {9}, {6,8,9,10,11}, {7..11},
  {2,4,6..11}, {5..11}, {4..11}`. **Identical to Lemma 4.2.** For `sigma = 9..12` I have the full attainable sets,
  not only the minimum. So my profile enumeration uses only attainable types and needs no relaxation.
* **Complement dimensions and explicit lists.** The types `(2,6), (4,6), (4,9), (6,2), (6,4), (7,5), (12,1)` occur
  only with a 1-dimensional complement. They give **66, 990, 495, 110, 1320, 9240, 11** lines, as claimed. Exact
  integer generators were recomputed over `Q` (`r_pools.py`, `fmpq_mat` nullspace), and each one's
  `(sigma, t)` was re-checked exactly. `(8,4)` and `(12,2)` occur only with a 2-dimensional complement, and `(8,5)`
  with both dimensions. These also agree.
* **Text error (minor).** Lemma 4.2 says the explicit lists have entries `|z_x| <= 2`. That is false for
  `(7,5)`. Its primitive generators have entries up to **3**, also on non-base rows, e.g.
  `z = (0,0,0,0,0,1,-3,-2,-1,2,2,1)`. The author's pools are built from flats (`univpool.py`), not from an entry
  bound, so I see no evidence that the search lost vectors. My search, which has no entry bound, reaches the same
  conclusion (§5). Still, the sentence must be corrected.
* **Profiles** (`r_profiles.py`). Multisets of 8 attainable `(sigma, t)` with `sum sigma <= 26` and `sum t <= 44`:
  exactly **25**. Of these, 13 are all-explicit, one is `7 pairs + (12,2)`, and 11 have one parametric type.
  Every profile has `>= 4` pairs. **This agrees with §4.2.**

## 5. Theorem 4.3: independent check of the width-count infeasibility

I re-proved Theorem 4.3 by a **different and simpler** route. It uses no symmetry, no knapsack weights, no
reserved weights, no column-requirement pruning and no `(12,2)` representatives.

**Lemma R (free 8th generator).** Let `z_1..z_7` be nonzero real vectors satisfying the 1519 counts, and let
`S = {v : N_v(z_1..z_7) = tau_A(v)}`. If `rank_Q S = 11`, then no nonzero real `z_8` extends them to a count-feasible
8-tuple.

*Proof.* `z_8` must satisfy `<z_8, v> = 0` for every `v in S`, so `z_8 = 0`. For integer `S`, `rank_p S <= rank_Q S <= 11`,
so `rank_p S = 11` certifies it. QED

**Reduction.** Every sub-multiset of a count-feasible multiset is count-feasible. For each of the 25 profiles I
pick a 7-element "rest" made only of explicit types (pairs, `(4,6)`, `(6,2)`, `(6,4)`). Then any count-feasible
8-tuple of that profile contains a count-feasible 7-tuple of that rest. The 8th generator is left **completely
free**: any real vector, any type. Eight rests cover all 25 profiles:

| rest (7 generators) | profiles covered (8th generator) |
|---|---|
| R1 = 7 pairs | `(6,2)`, `(12,1)`, `(12,2)` |
| R2 = 6 pairs + `(6,2)` | `(4,6)`, `(6,2)`, `(6,4)`, `(6,6)`, `(7,5)`, `(7,6)`, `(8,4)`, `(8,5)`, `(8,6)` |
| R3 = 6 pairs + `(6,4)` | `(6,4)`, `(8,4)` |
| R4 = 5 pairs + 2 `(6,2)` | `(3,9)`, `(4,6)`, `(4,8)`, `(4,9)`, `(4,10)` |
| R5 = 5 pairs + `(4,6)` + `(6,2)` | `(4,6)`, `(6,4)`, `(6,6)` |
| R6 = 5 pairs + `(6,2)` + `(6,4)` | `(4,8)` |
| R7 = 5 pairs + 2 `(6,4)` | `(4,6)` |
| (profile 24, done directly) 4 pairs + `(6,2)` + 3 `(4,6)` | all explicit: full 8-multiset enumeration |

**Enumeration (exact, integer only).** `r_pairs.c` lists every count-feasible multiset of `k` pairs `e_a - e_b`
(rows `0..11`), pruned only by the counts. There are 689,590 for `k = 4`, 2,638,339 for `k = 5`, 909,817 for `k = 6`
and 7307 for `k = 7`. `r_rest.c` extends each of them by the explicit slots, again pruned only by the counts, with
a bitset filter that drops candidates touching a saturated direction (a consequence of the counts).
`r_satrank2.py` then re-checks every leaf from scratch: balance, all 1519 counts recomputed, and `rank S` mod
`2^31 - 1`.

| rest | count-feasible 7-multisets | `rank_p S = 11` |
|---|---|---|
| R1 | 7307 | 7307 (all) |
| R2 | 72,210 | all |
| R3 | 58,995 | all |
| R4 | 31,725 | all |
| R5 | 627,895 | all |
| R6 | 69,740 | all |
| R7 | 41,290 | all |
| profile 24 (direct, all 8 explicit) | **0** count-feasible 8-multisets (689,590 pair multisets; 210M nodes in 3 parts) | n/a |

Cross-check with the author: R2 has 72,210 = 5 x 14,442 leaves, and the author's 14,442 used `Z_5` orbit
representatives for the first element. So the counts match exactly.

For profile 24 the rest "4 pairs + 2 `(4,6)` + `(6,2)`" has millions of count-feasible instances (the run was
stopped at > 47 MB of output). Since all 8 types of that profile are explicit, I enumerated the full profile
instead. `r_rest.c` with slots `(6,2), (4,6), (4,6), (4,6)` over all 689,590 count-feasible 4-pair multisets,
split into three parts, gives **0** count-feasible 8-multisets.

**Conclusion.** Each of the 25 profiles is excluded:
* 24 of them because every count-feasible instance of their rest has saturated rank 11 (Lemma R);
* profile 24 by direct exhaustion.

With Theorem 4.1 and the type/profile reduction (both re-derived in §3 and §4), **Theorem 4.3 is independently
confirmed**, including its claim for all **real** generators. The author's run data are consistent with mine
wherever they are comparable (R2: 14,442 x 5).

What my proof uses, explicitly:
* the exact type tables (§3, §4);
* the exact integer explicit lists;
* count-only enumerations;
* mod-`p` ranks only in the direction `rank_p = 11 => rank_Q = 11`.

It uses no symmetry and no `(12,2)` representative lemma.

## 6. Certification (rank_Q = rank_p) and the mod-p direction

* **Theorem 4.3 does not use Prop 3.1.** Its proof contains no membership test; the leaves are empty before one
  is reached. The mod-`p` steps it does use are the free-slot ranks (author: `rank_P(S) = 11` with
  `P = 2147483629`; mine: mod `2^31 - 1`). Both run in the rigorous direction, `rank_p <= rank_Q`. My type tables
  also use mod-`p` ranks, which are exact by the Hadamard bound (§3).
* **Prop 3.1** (`rank_Q K_{<=8} = 12093`) is argued correctly: `rank_p <= rank_Q <= #columns = |A|`. The
  isotypic block decomposition is legitimate. `x -> 3x` permutes `A` and the reduced monomials, the
  automorphism has order 5, and `p = 1 mod 5`, so the character projectors exist over `F_p` and conjugation
  preserves rank. I did **not** recompute the 12093 rank (a `31824 x 12093` matrix). It is not load-bearing for
  any claim that survives here. The "scaling bug fixed before use" cannot be audited from outside, and it does not
  matter for the same reason.
* The criterion "mod-`p` unsolvable implies `Q`-unsolvable when `rank_p K = rank_Q K`" is correctly stated and
  proved. The right-hand sides are `p`-integral because `p > 8!`.
* The claim `dim T_8 = 436` is proved only mod `p`, i.e. only `<= 436` over `Q`. The write-up says this, but
  §0/§1 still state "436 tangent octics" without qualification. Minor.

## 7. Round-13 correction (prodsearch cap)

* **The bug is real.** `round13/fatpoint_checks/prodsearch.py`, line 124: `if any(inc[x] >= 2 for x in S): continue`
  runs over all rows, the base row 0 included, with `capacity = 2*n` (line 105). Its docstring (line 9) states
  "`v = -1` (base incidence `<= 2`)". The true width of the all-ones direction is `tau_A(1)`, which is **6** for
  M8 = syl8 at `L = 2`, **8** for M8' = syl8 at `L = 3`, **6** for pal7 at `L = 2`, and **4** for pal11
  (`r_m8prep.py`, `r_pal7prep.py`, `r_box.py`). So the round-13 searches imposed a constraint that is not
  necessary, and were not exhaustive within their pool.
* **M8' at `m = 10` (the important reversal) is confirmed independently.** `r_msearch2.c` uses the same pool
  (support `<= 4`, entries `±1`; 238 generators) and directions `{-2..2}^7` with `tau < 10` (1570). It finds
  **1092** count-feasible 10-multisets, and `r_zono.py` shows exactly that **all 1092** are zonotopes (a translate
  of the full set of subset sums lies in `A`), hence certificates by fake Lemma 5.1(1). The explicit example
  `-e_1,-e_1,-e_2,-e_2,-e_3,-e_4,-e_5,-e_6,e_3-e_4,e_5-e_6` fits with `k_0 = (1,1,0,1,0,1,0)` (`r_m8p_zono2.py`). So
  `zeta(M8') >= 10`, and with R11's `<= 10`, `zeta(M8') = 10`. Round 13's "`M8'` at `m = 10` is zonotope-free
  within the pool" is **false**. Note: the M8' box is syl8 at `L = 3`, as the author's `mem_m8pm10_out.txt`
  header shows. The coordinates of the example are syl8 coordinates. It is **not** a zonotope in pal7 (`L = 3`)
  coordinates (`r_m8p_zono.py`: no translate). fake.md calls pal7 at `L = 3` "M8' up to a point permutation", so
  a reader who tries the example in pal7 coordinates will fail. The write-up should name the model (syl8, `L = 3`).
* **M8 at `m = 8` is confirmed independently.** There are **42** count-feasible products. `r_m8mem.py` computes
  `rank_Q K_{<=8} = rank_p K_{<=8} = 1037 = |A|` exactly (`fmpz_mat`), and **0 of 42** are members. That is
  rigorous by the rank criterion. "Certificates are not necessary" (M8, `m = 8`) survives.
* **pal7 at `m = 7`.** I find **2134** count-feasible products, of which **300** are zonotopes (exact). `r_m8mem.py` gives `rank_Q K_{<=7} = rank_p` (asserted exactly), 300 members mod `p`, and 1834 non-members, which is rigorous. Since every zonotope is a member, the 300 members are exactly the 300 zonotopes. So "every certificate found is a zonotope" is confirmed for pal7 as well.
* Not re-checked: sdc7 (45/24) and M8 at `m = 7` (23,355/1746).

## 8. Tangent cone statements (Prop 6.1, Prop 6.2, evidence)

* **Prop 6.1** is correct and nearly trivial. At fixed `n` the set of relevant `lambda` with `||lambda||_1 <= 101` in
  the discrete lattice `L_T` is finite. A nonzero zero-sum `lambda` has `lambda' != 0`. Nondegeneracy, stated for
  integer `v`, extends to rational `lambda'` by scaling. So a finite `a(u*)` exists. Quantifiers are fine. As the
  author says, it says nothing about `a` bounded in `n`, which is the real content of (2.2).
* **Prop 6.2** is confirmed independently (`r_coordtc.py`, `r_coordtc45.py`). Over all coordinate sets `S` with
  `|S| = 3, 4, 5, 6`, the sets with `rank K_{<=7}(pi_S A) = |pi_S A|` number **165, 115, 12, 0**. The 12 five-sets
  are exactly the twin patterns. 55 of the 115 four-sets are maximal (not inside a twin). Positives are mod `p`,
  which is the rigorous direction. For the 6-sets I also confirmed failure by exact `fmpz` rank. The logic is
  correct: `rank = |B|` gives `J(B)_8 = S_8`, so `TC(B, 8)` is everything, then fatpoint Prop 2.2(b) and Lemma 2.1
  apply. Two small points:
  * "No set with `|S| = 7, 8` passes" was not re-checked by me. It is a statement about a sufficient test, so
    nothing depends on it.
  * The sentence "the maximal certified coordinate subspaces are ..." means maximal among those certified by
    this test, not maximal coordinate subspaces of `TC`. It should say so.
* **§6.3** is correctly labelled EVIDENCE. Not re-run. The reading "on the degenerate side" is fair as evidence,
  but nothing more. Local minimisation of `F` cannot certify the absence of a nondegenerate component.
  Points within `3.4e-4` of collision strata with `F ~ 1e-20` are, as the author says, not distinguishable from
  nearby genuine zeros. Note also that `TC(A,8)` is a cone defined by 436 octics. A nondegenerate real zero, if it
  exists, can sit on a component of positive codimension. Random-start minimisation finds such components with
  probability that can be arbitrarily small. The evidence is therefore weak, and §7.3 correctly keeps the
  question open.

## 9. Issues list and corrected statement

| # | severity | where | issue | fix |
|---|---|---|---|---|
| 1 | scope | whole write-up | Addresses the paused route-5 item (README-13 §7.1), not this round's Rig(K) for P3-all\*. The P3-all\* class checks (parity, per-unit budget, monomials, moduli `q <= Q_K`, primes `p_j <= M^(h_0)`, `A` after `c_0`) have no object here. No growth consequence. | File as a route-5 record item, not as a round-14 answer. |
| 2 | minor | Lemma 4.2 | "entries `|z_x| <= 2`" is false. `(7,5)` primitive generators have entries up to 3 on non-base rows, e.g. `(0,0,0,0,0,1,-3,-2,-1,2,2,1)`. | Delete the bound, or say "`<= 3`". The author's pools are flat-generated, and my entry-free search agrees. |
| 3 | minor | Thm 4.1(b), Thm 4.3(a) | Necessity of the counts is quoted from fatpoint Prop 2.6(c), which is stated for `z_i in Z^d`. The theorems claim real `z_i`. | Add: the proof of Prop 2.6(b),(c) (top `t`-component on `u_0 + t v`) is valid verbatim for real `z_i`. |
| 4 | minor | §4.5 | The M8' example is in syl8 (`L = 3`) coordinates and is not a zonotope in pal7 (`L = 3`) coordinates, which fake.md calls "M8' up to point permutation". | Name the model: syl8, `L = 3`, `k_0 = (1,1,0,1,0,1,0)`. |
| 5 | minor | §0, §1.1 | "436 tangent octics": `dim T_8 = 436` is a mod-`p` value (`<= 436` over `Q`). | Qualify, as Prop 3.1 already does. |
| 6 | minor | Prop 6.2 | "maximal certified coordinate subspaces" means maximal among those passing this sufficient test. | Say so. |
| 7 | minor | §4.4 proof, step 3 | The free-slot reserves of `csearch8/9` are valid only because every weight vector is constant on its classes (checked in `mkweights.py`: core, twin, coordinates plus all-ones, `tau` classes). The text does not say so. | State it, or replace step 3 by the free-8th-generator argument of §5 here, which needs no reserves. |
| 8 | minor | §1.2 / Lemma 1.1 | Correct, and independently confirmed. Its scope is signed coordinate permutations only. | None needed beyond the remark. |

**Corrected statement.**

> **Theorem 4.3 (exact, computer-assisted; confirmed independently).** Let `A` be the pal11 box (`n = 12`, `T = 1`,
> `L = 2`, `|A| = 12093`), and let `V` be the 1519 directions `v` in `{-1,0,1}^11 ∪ {primitive v in {-2..2}^11 :
> |supp v| <= 4}` with `tau_A(v) < 8`. There are no nonzero real vectors `z'_1..z'_8 in R^11` with
> `#{i : <z'_i, v> != 0} <= tau_A(v)` for all `v in V`. Hence, by the necessity of the width counts (fatpoint
> Prop 2.6(c), whose proof holds for real `z_i`), `J(A)_8^perp` contains no product of 8 real linear forms, and
> neither does `L_8` or any `V_F(8)` whose family set contains `V`. This says nothing about whether `TC(A,8)` has
> a nondegenerate real point; that question is open. The statement is about one fixed box. It holds identically
> for fakes and for actual clusters with this pattern set, and it has no consequence for the growth of `M`:
> `M <= (2/3 + eps) log R/loglog R` remains the best bound. It does not bear on Rig(K) or on P3-all\*.
>
> **Round-13 correction (confirmed).** prodsearch.py capped base incidence at 2. Correct counts:
> * M8' (syl8, `L = 3`) at `m = 10`: 1092 count-feasible pool products, all zonotopes, so `zeta(M8') = 10`;
> * M8 (syl8, `L = 2`) at `m = 8`: 42 count-feasible products, 0 certificates (rigorous, `rank_Q K_{<=8} = rank_p = 1037`);
> * pal7 (`L = 2`) at `m = 7`: 2134 count-feasible, 300 zonotopes, exactly the 300 certificates (rigorous).

## 10. Scripts (`referee_pal11cert_checks/`, all written for this report)

| script | what | status |
|---|---|---|
| `r_box.py` | pal11 patterns, `A`, Hadamard check, 129,713 directions, 1519 with `tau < 8` | exact |
| `r_aut.py` | signed-permutation automorphisms of `A` (10) | exact |
| `r_flats.c`, `r_flats_sum.py` | brute-force `(sigma, t)` table over all supports and core subsets; min `t`, attainable types, line counts | exact (Hadamard-bounded minors `< p`) |
| `r_profiles.py` | the 25 profiles | exact |
| `r_pools.py` | exact integer explicit-type lists (66/990/495/110/1320/9240/11) | exact (`fmpq`) |
| `r_pairs.c` | all count-feasible pair multisets, `k = 4..7` | exact |
| `r_rest.c` | extension of pair multisets by explicit slots (rests R1-R7, profile 24) | exact |
| `r_satrank.py`, `r_satrank2.py` | leaf re-check (balance, all counts) and saturated rank mod `2^31 - 1` | exact (rigorous direction) |
| `r_coordtc.py`, `r_coordtc45.py` | Prop 6.2 coordinate subspaces | positives mod `p` (rigorous direction); 6-set negatives exact |
| `r_m8p_zono.py`, `r_m8p_zono2.py` | the M8' example: fails in pal7 `L = 3` coordinates, holds in syl8 `L = 3` | exact |
| `r_m8prep.py`, `r_pal7prep.py`, `r_msearch2.c` | round-13 models with true widths: M8 `m = 8` (42), M8' `m = 10` (1092), pal7 `m = 7` (2134) | exact |
| `r_zono.py`, `r_m8mem.py` | zonotope test; membership with exact `rank_Q = rank_p` | exact |
