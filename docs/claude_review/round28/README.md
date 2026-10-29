# Round 28 — Judge's README (multi-shift criterion vs zonotopes; pre-registered rule)

## 0. Verdict (short)

**Decision rule (pre-registered in round 27, README §11; restated in the round-28 task; not changed here): UNFAVOURABLE.**

* The favourable branch is **not met**. `reach > zeta` occurs on 2 of 105 boxes: M8' (`d = 7`) and syl16|J9a at
  `L = 3` (`d = 8`). The share is 1/6 at `d = 7`, 1/48 at `d = 8` and 0/24 at `d = 9`. Within the only class where it
  ever happens (Sylvester profiles at `L = 3`) it is 1/1 at `d = 7` and 1/5 at `d = 8` (1/6 with the verifier's
  J9a'). `reach - zeta` does not increase with `rho`.
* The unfavourable branch ("`reach <= zeta` on (essentially) every box beyond M8'") **holds with exactly one declared
  exception**, J9a. That is 1 of the 104 boxes beyond M8' (§4 explains why it is counted as "essentially").
* **Consequence under the rule:** Corollary 4.1 of round 27 (constant-chain `(SHIFT)_A` ⇒ `M = O_C(log log R)`) is
  treated as **presumptively vacuous**. Internal route-5 rounds **stop**. Questions A and B, the lemma ZE^all_H, the
  certified dP2 instance and the `Mbar_{0,7}` reformulation go to experts (§8).
* **Growth rate: unchanged.** `M <= (2/3+eps) log R/loglog R` (two referees) stands. **Uniform bound: unchanged
  (open).** No lane proved, and nothing here implies, any `M = o(log R/loglog R)` statement.
* Round score: **2 / 10** on the goal (§7). The round was decisive, honest and well verified, and it moved nothing.

Labels as in the record: **proved** (proof read and checked by the verifier; unrefereed by a third party),
**exact** (finite computation over `Q`/`Z` or a completed exhaustive search), **mod p** (ranks modulo primes `> 2^30`;
evidence), **heuristic**, **open**, **unrefereed** (judge's own items, §6).

## 1. Inputs and what was checked

| lane | document | role | verifier's finding |
| --- | --- | --- | --- |
| compute | `compute.md`, `boxes.csv`, `trend.csv`, `compute_checks/` | reproduce round 27; 105 boxes (`d = 5..9`), reach mod 2 primes, exact `zeta` | data confirmed, no numerical discrepancy |
| theory | `theory.md` (+ `theory_sec*.md`, `theory_notes.md`, `theory_checks/`) | upper bounds on reach, structure, regime | proofs correct up to minor gaps; every re-run claim reproduces |
| verify | `verify.md`, `verify_checks/` | independent code, own primes `2^31-1`, `2^30+3`; exact-over-`Q` certificates | — |

I read all three documents in full and the round-27 README §§0, 1, 5, 10, 11. My own checks (§6, all **unrefereed**)
cover the CSV bookkeeping, every zonotope witness of all 105 boxes, the face structure of J9a and J9a', and an
independent "no order-11 zonotope" search on J9a.

**Definitions used (as checked by all three lanes and by me).** A slab box is
`A = {k' in {-1,0,1}^d : |kfull(S_p)| <= L for every cut S_p}`, with `kfull = (k', -sum k')` and the base point last.
The base coordinate is **not** itself restricted to `{-1,0,1}`. In the `beta` normalisation the half-width is `2L`.
The task's "`|<k, beta_p>| <= L`" is the same box with `L` rescaled by 2, so `rho/L` must be read in the record's
convention. `reach(A, L)` = the first degree at which `(SHIFT_e)` fails; it is always `<= rho`.
`zeta_pair` uses generators in `L_pair` (the task's definition); `zeta_all` uses all primitive `v in {-2..2}^d`
(the class behind the record's value `zeta(M8') = 10`). On both `reach > zeta` boxes the two coincide, and in the
rule's reading they give identical counts of `reach > zeta` (§6 J1).

## 2. Verified data table

**Status of the numbers.**
* `|A|`, `H_A`, `rho` and reach are **mod p**. Each was computed at two primes by compute (`P1 = 2147483629`,
  `P2 = 2147483587`), and the two agree on all 105 boxes.
* The verifier recomputed 13 boxes plus the 27 round-27 random boxes at two further primes (`2^31-1`, `2^30+3`),
  with independent code. Every recomputed value matches.
* **Exact over `Q`** (verifier, counting certificate): reach(M8', `L_pair`) = 11, and reach(M8, `L_pair`) = 6 with
  an exactly one-dimensional defect spanned by `sigma^6`.
* Every `zeta` is **exact**. Lower bounds come from witnesses: I re-verified all 210 of them on boxes rebuilt with my
  own code (§6 J3). Upper bounds come from an exact rational LP dual, plus a completed exhaustive search where
  `floor(LP)` is not attained. Two lanes did this independently on 16 boxes; for J9a there is a third, my J4.
* Truncated cube-shift reaches (marked ≥) are lower bounds.

**2.1 All structured boxes** (from `boxes.csv`; `zeta` column = `zeta_all`; spec syntax as in `compute.md` §7:
`syl:k:L` Sylvester `H_{2^k}`; `sdc:q:L`; `rsyl:4:J:L` / `rsdc:11:J:L` / `rhad:11:J:L` = Sylvester-16 / skew double core
`q = 11` / Paley-I Hadamard 12 restricted to the ordered point list `J`, last entry = base; `syl:3:3` = M8',
`syl:3:2` = M8, `rsyl:4:0,...,8:3` = J9a, `rsyl:4:0,1,2,4,8,15,3,5,6:3` = J9b):

| spec | d | L | cuts | `|A|` | `rho` | `rho/L` | reach `L_pair` (P1/P2) | first failure | reach cube | `zeta_pair` | `zeta_all` (LP) | reach − `zeta` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `sdc:7:2` | 7 | 2 | 14 | 757 | 8 | 4.0 | 6 / 6 | 166/188 | 6 | 7 | 7 (7) | -1 |
| `syl:3:2` | 7 | 2 | 7 | 1037 | 8 | 4.0 | 6 / 6 | 286/287 | 7 | 7 | 7 (7) | -1 |
| `syl:3:3` | 7 | 3 | 7 | 1877 | 11 | 3.667 | 11 / 11 | none below `rho` | 11 | 10 | 10 (21/2) | +1 |
| `sdc:7:3` | 7 | 3 | 14 | 1765 | 10 | 3.333 | 10 / 10 | none below `rho` | 10 | 10 | 10 (21/2) | +0 |
| `syl:3:1` | 7 | 1 | 7 | 57 | 3 | 3.0 | 3 / 3 | none below `rho` | 3 | 3 | 3 (7/2) | +0 |
| `sdc:7:1` | 7 | 1 | 14 | 57 | 3 | 3.0 | 3 / 3 | none below `rho` | 3 | 3 | 3 (7/2) | +0 |
| `rsyl:4:0,1,2,3,4,5,6,7,8:1` | 8 | 1 | 15 | 73 | 3 | 3.0 | 3 / 3 | none below `rho` | 3 | 3 | 3 (15/4) | +0 |
| `rsyl:4:0,1,2,4,8,15,3,5,6:1` | 8 | 1 | 15 | 73 | 3 | 3.0 | 3 / 3 | none below `rho` | 3 | 3 | 3 (15/4) | +0 |
| `rsdc:11:0,1,2,3,4,5,6,7,8:1` | 8 | 1 | 22 | 73 | 3 | 3.0 | 3 / 3 | none below `rho` | 3 | 3 | 3 (11/3) | +0 |
| `rhad:11:0,1,2,3,4,5,6,7,8:1` | 8 | 1 | 11 | 73 | 3 | 3.0 | 3 / 3 | none below `rho` | 3 | 3 | 3 (11/3) | +0 |
| `rsyl:4:0,1,2,3,4,5,6,7,8:2` | 8 | 2 | 15 | 1837 | 8 | 4.0 | 6 / 6 | 572/580 | 7 | 7 | 7 (15/2) | -1 |
| `rsyl:4:0,1,2,4,8,15,3,5,6:2` | 8 | 2 | 15 | 1841 | 8 | 4.0 | 6 / 6 | 576/577 | 7 | 7 | 7 (15/2) | -1 |
| `rsdc:11:0,1,2,3,4,5,6,7,8:2` | 8 | 2 | 22 | 1367 | 7 | 3.5 | 6 / 6 | 336/381 | 6 | 7 | 7 (22/3) | -1 |
| `rhad:11:0,1,2,3,4,5,6,7,8:2` | 8 | 2 | 11 | 2049 | 8 | 4.0 | 7 / 7 | 12/412 | 7 | 7 | 7 (22/3) | +0 |
| `rsyl:4:0,1,2,3,4,5,6,7,8:3` | 8 | 3 | 15 | 4861 | 11 | 3.667 | 11 / 11 | none below `rho` | — | 10 | 10 (11) | +1 |
| `rsdc:11:0,1,2,3,4,5,6,7,8:3` | 8 | 3 | 22 | 4457 | 10 | 3.333 | 9 / 9 | 604/606 | — | 10 | 10 (1316/121) | -1 |
| `rhad:11:0,1,2,3,4,5,6,7,8:3` | 8 | 3 | 11 | 4997 | 11 | 3.667 | 10 / 10 | 342/368 | — | 11 | 11 (11) | -1 |
| `rsyl:4:0,1,2,3,4,5,6,7,8,9:2` | 9 | 2 | 15 | 3847 | 8 | 4.0 | 6 / 6 | 1250/1251 | 7 (≥, truncated) | 7 | 7 (15/2) | -1 |
| `rsdc:11:0,1,2,3,4,5,6,7,8,9:2` | 9 | 2 | 22 | 2709 | 7 | 3.5 | 6 / 6 | 893/895 | 6 (≥, truncated) | 7 | 7 (22/3) | -1 |
| `rhad:11:0,1,2,3,4,5,6,7,8,9:2` | 9 | 2 | 11 | 4525 | 8 | 4.0 | 7 / 7 | 1112/1227 | 7 (≥, truncated) | 7 | 7 (22/3) | +0 |
| `rsyl:4:0,1,2,3,4,5,6,7,8,9:1` | 9 | 1 | 15 | 91 | 3 | 3.0 | 3 / 3 | none below `rho` | 3 | 3 | 3 (15/4) | +0 |
| `rsyl:4:0,1,2,4,8,15,3,5,6:3` | 8 | 3 | 15 | 4787 | 11 | 3.667 | 10 / 10 | 215/281 | 10 (≥, truncated) | 10 | 10 (11) | +0 |
| `rsyl:4:3,11,2,8,9,13,15,10,7:3` | 8 | 3 | 15 | 4809 | 11 | 3.667 | 10 / 10 | 218/284 | — | 10 | 10 (89/8) | +0 |
| `rsyl:4:13,3,14,2,6,12,10,15,8:3` | 8 | 3 | 15 | 4809 | 11 | 3.667 | 10 / 10 | 218/284 | — | 10 | 10 (89/8) | +0 |
| `rsyl:4:0,6,3,1,12,14,5,2,8:3` | 8 | 3 | 15 | 4809 | 11 | 3.667 | 10 / 10 | 218/284 | — | 10 | 10 (89/8) | +0 |
| `rhad:11:8,11,4,3,7,1,6,0,2:3` | 8 | 3 | 11 | 4981 | 11 | 3.667 | 10 / 10 | 334/359 | — | 11 | 11 (11) | -1 |
| `rhad:11:10,11,4,1,9,5,0,7,3:3` | 8 | 3 | 11 | 4997 | 11 | 3.667 | 10 / 10 | 342/368 | — | 11 | 11 (11) | -1 |
| `rsdc:11:2,3,5,1,0,10,6,11,8:3` | 8 | 3 | 22 | 4465 | 10 | 3.333 | 10 / 10 | none below `rho` | — | 10 | 10 (1646/151) | +0 |
| **`rsyl:4:1,2,3,4,5,6,7,8,0:3` = J9a'** (verifier only; base inside the hyperplane) | 8 | 3 | 15 | 4791 | 11 | 3.667 | 9 / 9 (verifier's primes) | 740/741 | 9 (≥, 150 classes) | 10 | 10 (43/4) | −1 |

**2.2 Random boxes** (the round-27 killer's law: `n` points, `r ~ U{n-1..2n}` distinct random canonical cuts,
`L ~ U{1,2}`, `random.Random(seed)`; full rows in `boxes.csv`):

| group | boxes | zonotopal (`zeta_all = rho`) | reach vs `zeta_all`: > / = / < | mean reach − `zeta_all` |
| --- | --- | --- | --- | --- |
| `n = 6` seeds 100–114 (`d = 5`) | 15 | 15 | 0 / 5 / 10 | −0.73 |
| `n = 7` seeds 1–12 (`d = 6`) | 12 | 12 | 0 / 2 / 10 | −1.17 |
| `n = 9` seeds 1–30 (`d = 8`) | 30 | 28 | 0 / 0 / 30 | −1.87 |
| `n = 10`, first 20 seeds in 1–38 drawing `L = 1` (`d = 9`) | 20 | 20 | 0 / 1 / 19 | −1.30 |
| **all random** | **77** | **75** | **0 / 8 / 69** | — |

**2.3 Summary by `d`** (all 105 boxes; recounted from the CSV in §6 J1):

| `d` | boxes | faithful (`zeta < rho`) | reach > / = / < `zeta_all` | Sylvester `L = 3` boxes with reach > `zeta` |
| --- | --- | --- | --- | --- |
| 5–6 | 27 | 0 | 0 / 7 / 20 | — |
| 7 | 6 | 3 | 1 / 3 / 2 | 1 of 1 (M8') |
| 8 | 48 | 10 | 1 / 10 / 37 | 1 of 5 (J9a); 1 of 6 with J9a' |
| 9 | 24 | 2 | 0 / 3 / 21 | no `L = 3` box computed (memory) |

The counts of `reach > zeta` are the same against `zeta_pair` (M8' and J9a only). Against `zeta_pair` there are
more ties; this is recorded in J1.

**Facts that bear on the rule (all verified):**
1. Exactly two boxes have `reach > zeta`, M8' and J9a, both with `reach = rho = 11` and `zeta = 10`. On J9a, `zeta` is
   exactly 10 although the LP bound is exactly 11: three independent exhaustive searches find no order-11 zonotope
   (compute, verify, and my J4).
2. J9a's 9-point set contains M8''s 8-point configuration (the affine hyperplane `{0..7}`), but its **base** is the
   extra point. So the J9a box is not an extension of the M8' box; its face is `M8' cap {|k_7| <= 1}` (1093 points;
   exact, J2). The box that *does* extend M8' (J9a', base inside the hyperplane, face `{k_8 = 0}` = the M8' box;
   exact, J2) has `reach = 9 < zeta = 10`.
3. The other four `d = 8` Sylvester `L = 3` restrictions tie (`reach = zeta = 10 < rho = 11`). On J9b, the 300
   largest cube-shift classes do not raise this.
4. At `L = 2` the Sylvester boxes repeat the M8 pattern at `d = 7, 8, 9`: pair reach is `zeta - 1` and cube reach is
   `zeta`. The Hadamard boxes tie. The sdc boxes are zonotopal and fall one short.
5. Every box with `rho/L > 4` (59 boxes) is random, and 57 of them are zonotopal. None has `reach > zeta`. The
   verifier's count "53, 51 zonotopal" is a minor miscount (J1); the conclusion is the same. So the decline of
   `reach - zeta` with `rho/L` (compute Table 6.4) is a random-versus-structured comparison, not a within-family
   trend, as the verifier says.

## 3. Theory results and verified status

| item (`theory.md`) | statement (short) | status after verification | relevance |
| --- | --- | --- | --- |
| Thm 2.1 (products) | `H` convolves; `rho`, `reach(L_pair)`, `zeta_pair`, `zeta_all` add under products of boxes | **proved**. Minor gap in (iii): a projected generator may be non-primitive; harmless, because `m v'` splits into `m` copies of `v'` in a slab box. Checked numerically on 4 products | structural |
| Cor 2.2 | `(M8')^t`: reach `11t`, `zeta 10t`, so no bound `reach <= zeta + c` (or `<= (1+o(1)) zeta`) holds on slab boxes | **proved, now unconditional**: its input reach(M8') = 11 is **exact over `Q`** (verifier §6) | vacuous for the route: block profiles violate the character gap. Products are not in the rule's test population (§4) |
| Thm 2.3 (implied-slab obstruction) | if `phi` takes `R` values on `A` and `<= R-1` on every `A_l`, then `phi^{R-1} in cap_l J(A_l)`; with `(N_phi)` the reach is `<= R-1` | **proved** | explains most model failures |
| Cor 2.4 (BCO) + Prop 2.5 | balanced profile + rigid top level of `sigma` + `(N)` ⇒ `reach(A, L_pair) <= 4j <= 2w < m' - 1` | **proved**. Minor gap: "`R - 1 <= 2h`" needs `min sigma >= -h`, which holds for symmetric slab boxes | conditional Liouville cap; its rigidity hypothesis holds only at designs (Sylvester `L = 2^m`; IP up to `n = 128`, partly undecided) — **heuristic** that it fails generically |
| §2.3 rigidity data | (R) holds for syl8 `L = 2`, Sylvester `n = 16`, `L = 2, 4`, and 3 of 27 random boxes (each with reach `= 2h`) | **exact** where re-run (`n = 8, 16`, 27 random boxes); IP rows `n >= 32` not re-run | — |
| Thm 3.1 (chain/subset criterion) | any chain `S_0 ⊂ ... ⊂ S_m ⊂ A` with `(CH_e)` certifies ZE-braid at order `m`; so `max(zeta_L, max_{A'⊂A} reach(A')) <= reach_chain <= truth` | **proved** (verifier) | the zonotope chain is a special case: chains dominate zonotopes *by construction* (§5) |
| Lemma 3.2 (first order) | `I(A_l) = I(A) + tau_l I(A)`; `J(A_l) ⊃ J(A) + R d_l J(A)` | **proved** | — |
| Prop 3.3 (linkage form) | `(SHIFT_e)` ⇔ a statement about the tail `B = C \ A` in degree `2d - e` | identity checked mod p on 3 boxes; **duality step not independently re-derived** | used nowhere else |
| Thm 4.1 | the M8 (syl8 `L = 2`) failure at `e = 6` is `sigma^6`; the defect is exactly 1-dimensional | **exact over `Q`** (verifier: 56 signed shifts, `sigma^7 in J(A)`, defect certified) | a design effect (Fano lines) |
| §3.5 | M8' succeeds: top level of `sigma` not rigid (28 points, 252 legal pair moves) | **exact** | M8' is a weak-cut (vertex-cut), design-free box; not regime-like |
| §3.4 | 6 of 8 model defects spanned by powers of ternary forms; sdc7 `L = 2` only 7 of 22 dims | mod p, **not re-run** | — |
| §5.2 item 4 | "no counting obstruction in the regime" | **heuristic** (verifier: asserted, not proved) | — |
| §5.3 | no known asymptotic tool (Macaulay, Gotzmann, gin, linkage, staircase) decides `(SHIFT)` at regime size | judgement; I agree | — |

**Net theory content.** Nothing here bears on growth. The lane settles the *shape* of the question:
* Any bound `reach <= zeta + c` must use the character gap (Cor 2.2).
* The one proved Liouville-capped obstruction (BCO) needs a design (Cor 2.4 and §2.3).
* The open part is a structureless "non-univariate residue", together with the absence of rigid implied slabs at
  regime size.
* The verifier found a new instance of that residue: J9a' fails at `e = 9` with no implied-slab explanation
  (`t6_J9ap_out.txt`).

## 4. Reading of the pre-registered decision rule

**The rule, as fixed before the data.** From round-27 README §11, restated in the round-28 task:
* *Favourable:* `reach > zeta` on a growing share of boxes as the boxes grow (`d = 7 -> 8 -> 9`, `rho/L` increasing),
  with `reach - zeta` increasing in `rho` (round-27 wording).
* *Unfavourable:* `reach <= zeta` on (essentially) every [multi-cut] box beyond M8'.

Here `reach` is the constant-chain reach `reach(A, L_shift)` of the task's definitions, and the boxes are the
round-11/12 slab-box families (Sylvester, sdc and their subcluster restrictions) plus the multi-cut random boxes.

**Favourable branch: not met.**
* The share of `reach > zeta` is 0/27 (`d <= 6`), 1/6 (`d = 7`), 1/48 (`d = 8`), 0/24 (`d = 9`).
* On faithful boxes (`zeta < rho`, the only boxes where `reach > zeta` is possible) it is 1/3, 1/10, 0/2.
* In the only class where it ever happens, Sylvester at `L = 3`, it is 1/1 at `d = 7` and 1/5 at `d = 8` (1/6 with
  J9a').
* `reach - zeta` takes the values +1, +1 and then 0 or −1 on structured boxes. It is never above +1, and it does not
  increase with `rho`.
* Every one of these restrictions points the same way: the share falls.

**Unfavourable branch: met, with one declared exception.**
* `reach <= zeta` on 103 of the 104 boxes beyond M8'.
* The exception is J9a (`d = 8`, Sylvester-16 restricted to `{0..8}`, `L = 3`), with reach 11 (mod p, four primes)
  and `zeta` 10 (exact, three independent searches).

I count it under "essentially" for four reasons, each verified:
1. It is one box in 104.
2. Its point set is M8''s configuration plus one point; it is the only one of the five `d = 8` Sylvester `L = 3`
   restrictions that contains an affine hyperplane.
3. The genuine one-point extension of the M8' *box* (J9a', which contains M8' as a coordinate face) has
   `reach 9 < zeta 10`.
4. The other four restrictions tie.

So the phenomenon is attached to one small symmetric configuration with weak (vertex) cuts. It does not propagate
to its natural extension, and it is not a feature of growing boxes.

**What I did not count, and why (no goalpost moved).**
* *Product boxes `(M8')^t`* (theory Cor 2.2; reach `11t` > `zeta 10t`, now unconditional). They are not in the
  pre-registered test population (round-27 §11 names the Sylvester/sdc slab families and multi-cut random boxes),
  and they violate the character gap that every actual cluster satisfies. Adding them would be a post-hoc family
  built to meet the favourable branch; the theory lane says the same.
* *Chain or subset reach* (theory Thm 3.1). The rule is stated for the constant-chain reach. A subset statistic
  would change some numbers: J9a' contains the M8' box, so its subset reach is `>= 11 > zeta`. But it would not
  change the picture, because every such success is inherited from M8' by containment. It is reported in §5 and not
  used for the reading.
* *Swapping `zeta_all` for `zeta_pair`*. This changes no `reach > zeta` count (J1).

**Strict reading.** If "essentially" were dropped and "every" read literally, J9a would leave the unfavourable branch
unmet as well. That reading would give **undecided**, never favourable. The pre-registered text includes
"(essentially)", and one box out of 104, explained by containment of M8''s configuration and not surviving the true
extension J9a', is what that word was written for. **Reading: UNFAVOURABLE.**

**Residual uncertainty (stated, not used to soften the reading).**
* (i) There is no `d = 9` box at `L = 3` (memory). This is the class where the phenomenon lives. For the favourable
  branch to have been met, the `d = 9`, `L = 3` share would have had to rise after falling from 1/1 to 1/5 at
  `d = 8`. I put this at under 10% (judgement).
* (ii) Every reach other than M8 and M8' is mod-p evidence. The primes agree on all boxes, and four primes agree on
  J9a.
* (iii) All models are LP-subcritical. Structured boxes have `rho/L <= 4` against the regime ratio `m'/j = 8` (record
  convention); the decisive regime window `e/j in (4, 8)` is never reached by a structured box.
* (iv) 75 of the 77 random boxes are zonotopal, where `reach > zeta` is impossible a priori. The faithful-box samples
  are small (3, 10, 2).

These caveats are why the rule says *presumptively* vacuous. They do not turn the data favourable.

## 5. Consequences for route 5 and for the growth rate

**Growth rate: unchanged.** `M <= (2/3+eps) log R/loglog R` (two referees) stands. Nothing in round 28 proves,
conditionally on anything weaker than what the record already had, any `M = o(log R/loglog R)`. **Uniform bound:
unchanged**; Question B is open in every form.

**Route 5.**
1. *Cor 4.1 (round 27) is treated as presumptively vacuous*, as the rule prescribes. In every computable model,
   the constant-chain multi-shift criterion certifies ZE-braid no further than zonotopes do, except on M8''s own
   configuration. The gap `zeta - reach` is 1–3 on faithful random boxes. A criterion that does not outrun the
   Liouville-capped zonotope order cannot deliver the `ceil(4w)` that Cor 4.1 needs, since `ceil(4w)` is about
   twice the zonotope cap `2w`.
2. *What the round does establish about route 5* (proved or exact):
   * `reach - zeta` is unbounded on slab boxes in general (products), so any regime statement must use the
     character gap.
   * The proved Liouville-capped obstruction (BCO) needs a design.
   * The verified failures include a structureless residue (J9a' at `e = 9`; most of sdc7 `L = 2`).
   * No asymptotic tool decides `(SHIFT)` at regime size.

   This sharpens, but does not change, the record's barrier: route 5 needs ZE^all_H beyond the zonotope order, and
   round 28 found no internal mechanism for it.
3. *The chain/subset version (theory Thm 3.1) is not closed by this rule, and does not justify another internal
   round.*
   * Chains dominate zonotopes by construction, so for chains "beating `zeta`" is not the right test.
   * The only evidence of strict domination is M8' and boxes that contain it. My J2 confirms exactly that J9a'
     contains the M8' box as the face `{k_8 = 0}`. With the exact reach(M8') = 11, this gives
     `reach_chain(J9a') >= 11 > zeta(J9a') = 10` (**unrefereed**, §6).
   * Whether chains reach `~8j` on gap-satisfying regime profiles is a reformulation of ZE^all_H. It goes to experts
     as part of that lemma (§8), not as a new internal route.
4. *Probabilities (judgement).*
   * Round 27 put a favourable trend at about 15%. It came out unfavourable.
   * A proved growth improvement along route 5 by internal rounds: about 1% before, now **< 0.5%**.
   * An expert pointer to `alpha^ess > tau` at a rational point with `tau < dim`: about 3–5%, unchanged.

## 6. Judge's own checks (all **unrefereed**; scripts and outputs in `judge_checks/`)

| id | check | method | result |
| --- | --- | --- | --- |
| J1 | recount of `boxes.csv` (`j1_recount.py` → `j1_out.txt`; table `j1_table.md`) | exact integer bookkeeping | 105 boxes, P1 = P2 and `H` agreement everywhere, every `zeta` flagged exact. `reach > zeta` only on `syl:3:3` (M8') and J9a, against **both** `zeta_all` and `zeta_pair`. Per-`d` counts as in §2.3. All 59 boxes with `rho/L > 4` are random (57 zonotopal); the maximum structured `rho/L` is 4. (The verifier's "53 / 51" is a minor miscount; the conclusion is unchanged.) Random: 75/77 zonotopal, 69/77 `reach < zeta`. |
| J2 | M8', J9a, J9a' rebuilt from scratch (own Sylvester/restriction code; `j2_boxes_witnesses.py` → `j2_out.txt`) | exact enumeration of `{-1,0,1}^d` | `|A|` = 1877, 4861, 4791 (base coordinate unconstrained, as in the record). J9a face `{k_base = 0}` = M8' ∩ `{|k_7| <= 1}` (1093 points), **not** the M8' box. J9a' face `{k_8 = 0}` = the M8' box exactly (1877; also checked to be the translate `x -> x XOR 7` of the base-7 M8' box). The order-10 witnesses of J9a (both classes), M8' (both) and J9a' (pair) have all `2^10` vertices in `A` and generators in `L_pair`. |
| J3 | all 105 boxes rebuilt with own box code from the compute lane's profile lists (`j3_all_witnesses.py` → `j3_out.txt`) | exact | all 105 `|A|` match the CSV. **All 210 `zeta` witnesses** (pair and all) have every vertex in `A` and the recorded order; every pair witness is in `L_pair`. So `zeta >= ` the CSV value is exact on every box, and `reach_pair < zeta_pair` is certified (given the mod-p reach) on 78 boxes. |
| J4 | no order-11 zonotope in J9a, any primitive generators in `{-2..2}^8` (`j4_J9a_no11.py` → `j4_out.txt`) | own LP dual (24 functionals: coordinates, base coordinate, 15 slabs), rationalised, checked in **exact integer arithmetic** over all 155 263 fitting candidates; then own exhaustive DFS over the dual-reduced set | exact dual bound `B = 11` with `min c_v = 1`. Only the 36 `L_pair` classes have `c_v = 1`. Exhaustive DFS over them: **no order-11 zonotope** (3 506 942 nodes, 219 s); the same search finds order 10 at once. So `zeta(J9a) = 10` exactly: a third independent confirmation. |

**New claims of my own (unrefereed):**
* (i) The `rho/L > 4` population is entirely random; the verifier's count needs correcting (J1).
* (ii) `reach_chain(J9a') >= 11 > zeta(J9a') = 10`. This combines theory Thm 3.1 (subset case, verified proof), the
  exact face identity J2 and the verifier's exact reach(M8') = 11. Since `A' = M8'` lies in the hyperplane
  `{k_8 = 0}`, the `L_pair(8)` shifts that move `k_8` give empty `A'_l`, and the remaining ones are `L_pair(7)`.
  The leading forms of `A'` in 8 variables are those in 7 variables. Hence `reach(A', L_pair(8)) = 11`.
  Consequence: truth(J9a') >= 11 = `rho` (mod p).
  This is the M8' phenomenon inherited by containment, not new evidence about growing boxes, and it does not enter the
  rule's reading (§4).

## 7. Scores

Rubric of `round10/README.md` §1, as in rounds 22–27: (a) correctness after verification; (b) fake/barrier status of
the input (for a test round: whether the test was decisive and honestly read); (c) movement of growth or uniformity.
No movement of growth or uniformity caps the overall score at 3.

| lane | (a) | (b) | (c) | overall | reasons |
| --- | --- | --- | --- | --- | --- |
| **compute** | 8 | 6 | 1 | **3** | Strengths: the pre-registered test was executed as specified, with 105 boxes, two primes, exact `zeta` everywhere, a correct refinement lemma, every round-27 value reproduced, and zero discrepancies under independent recomputation. Deductions: J9a was called "M8' plus one point" without checking which box contains M8' (the true extension J9a' fails); `zeta_all` was used where the task defines `zeta_pair` (harmless here); the `rho/L` trend was presented without noting that it is confounded with random versus structured boxes; there is no `d = 9`, `L = 3` coverage. |
| **theory** | 7 | 4 | 1 | **3** | Correct new structure: product additivity (unbounded `reach - zeta`), the implied-slab obstruction / BCO with an exact design explanation of the M8 failure (`sigma^6`), and the chain criterion dominating zonotopes. Honest that the rule is not settled. Deductions: minor proof gaps (non-primitive projections; `min sigma >= -h`); one heuristic stated as fact (§5.2 item 4); Prop 3.3's duality not independently checked; "design effect" reasoning applied asymmetrically, as the verifier notes. |
| **verify** | 9 | 7 | 1 | **3** | Independent code and primes. Exact-over-`Q` certificates upgrade reach(M8') = 11 and the M8 defect to theorems. The adversarial box J9a' is the most informative single data point of the round. Coverage and confounding caveats are correct. Deductions: a minor miscount (53 vs 59 boxes with `rho/L > 4`); J9a reach left mod p (justified by cost). |

Round as a whole: **2 / 10** on the goal. The round did what a decisive test round should: it applied a fixed rule
to verified data, and the rule closed route 5's interpolation branch for internal work. It moved neither the growth
rate nor the uniform bound.

## 8. Next step implied by the rule

The rule's unfavourable branch prescribes the next step, and I follow it without amendment.

1. **Stop internal route-5 rounds.** No further internal computation or theory on `(SHIFT)_A`, its chain version,
   ZE-braid or ZE^all_H. This is also the standing recommendation since round 22, now triggered by a pre-registered
   test.
2. **Send to experts** (one package, an update of the record's expert note; text proposals only, since the
   repository is read-only for this round):
   * **Question A and Question B** as in the record.
   * **Lemma ZE^all_H** (height-cut ZE-braid at the cluster point), with its round-27/28 forms:
     * the multi-shift criterion `(SHIFT_e)` and its chain version (theory Thm 3.1);
     * the exact facts reach(M8') = 11 > `zeta` = 10, M8's `sigma^6` obstruction, and J9a' reach 9 < `zeta` 10
       with an unexplained residue;
     * the product theorem, showing that any bound must use the character gap.

     The precise question: *on gap-satisfying slab boxes at regime size, does any linear-algebra certificate
     (subset/chain criterion) reach order `~8j` beyond the zonotope cap `~4j`, or is there an obstruction family of
     degree `<= 2w(1+o(1))`?*
   * **The certified dP2 instance** (round 27 §7: `Pic = ZK`, `tau = 3/2`, real-cubic `tau`-curve).
   * **The `Mbar_{0,7}` reformulation** (Castravet–Tevelev: `Bl_e LM_5` is a contraction of `Mbar_{0,7}`, so
     `tau_e(-K)`, Theorem B, WTB and Conjecture E are statements about `Eff(Mbar_{0,7})`).
3. **Barrier-map update (text only).** Route 5, interpolation branch, gets a new entry, `closed for internal work
   (round 28, pre-registered test)`, with this content:
   * `reach(L_pair) <= zeta` on 103 of 104 slab boxes beyond M8' (`d = 5..9`);
   * the single exception, J9a, contains M8''s configuration, and the true extension of M8' fails;
   * `reach - zeta` is unbounded only on block (gap-violating) profiles;
   * the BCO caps the criterion at the zonotope scale on design profiles.
4. **Optional, low priority, not a new round.** If compute capacity is free, the `d = 9`, `L = 3` Sylvester
   restrictions (`|A|` ≈ 11.5k–13.4k; would need a sparse/blocked `H_A`) would close coverage gap (i) of §4. The rule
   does not require it, and its outcome cannot make the reading favourable unless the share *rises* from 1/5. It
   should not be run as an internal round.

## Files

All under `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/`:
`compute.md`, `boxes.csv`, `trend.csv`, `compute_checks/`; `theory.md`, `theory_checks/`, `theory_notes.md`;
`verify.md`, `verify_checks/`; this `README.md`, `judge_checks/` (`j1_recount.py`, `j1_out.txt`, `j1_table.md`,
`j2_boxes_witnesses.py`, `j2_out.txt`, `j3_all_witnesses.py`, `j3_out.txt`, `j4_J9a_no11.py`, `j4_out.txt`),
`judge_notes.md`.
