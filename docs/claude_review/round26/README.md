# Round 26 — judge report (zero / multiplicity estimates for route 5's ZE at the actual point)

| key | write-up | referee report | question |
| --- | --- | --- | --- |
| ze | `ze.md` (`ze_checks/`, `notes.md`) | `referee_ze.md` (`ref_ze/`) | Can a Philippon-type zero or multiplicity estimate on a commutative algebraic group prove route 5's zero estimate at the actual point `Z` (ZE-col\*, MIN^pair, `(R_H)`)? This is the one internal item on which rounds 11–17 and 25 converge (round-25 README §8). |

All paths are relative to `scratchpad/round26/`, and the repository was only read. Labels follow the overview:
proved / conditional / exact / evidence / heuristic / unrefereed. §1 gives the **referee-corrected** statement.
Remarks marked *judge* are mine and are **unrefereed**. My computations are in `judge_checks/` (§4).

## 0. Bottom line

**No growth improvement and no uniform bound.** `M <= (2/3 + eps) log R/loglog R` (`two_thirds/proof.md`) is still the
best bound for every `R`. ZE-col\*, MIN^pair and `(R_H)` are **not** derived, and route 5's zero estimate at the
actual point is exactly as open as at the end of round 25. The write-up claims no growth or uniformity, and the
referee found no step that does.

**What the round settled (refereed).** The two natural ways of feeding route 5's data to a Philippon-type zero
estimate are **proved useless**, with the trivial subgroup as obstruction:
* *Torus picture.* One form `P_c` on `G_m^d` with vanishing data at the powers `Z^j`. Route 5 certifies `|j| <
  ord(c)/Lambda(c) <= 1.07 sqrt(n)/lambda` (new spectral LP cap), i.e. `|j| <= 1..4` on the regime boxes, while
  Philippon needs a progression of length about `d`. Conclusion (2.1) holds at `H = {e}`, with margins `>= e^{1.7e4}`.
* *Dual picture.* One exponential polynomial on `G_a^d x G_m` vanishing on the slab box. The obstruction subgroups
  are `V x {1}`, `V x G_m`, every class count is point-independent at nondegenerate points, the degree-level
  estimate is vacuous and its Hilbert-sharp form is false.

**What did not survive as stated (four major corrections).** "Route 5 certifies only the powers `Z^j`" is false
(endomorphism images `Z^M`, in particular profile automorphisms, are certified too: M1). The universal "no
single-form zero / multiplicity / small-value estimate with any constants" barrier is downgraded to powers-only data,
`lambda >= lambda*(n,mu)`, geometric (not small-value) estimates, and four values of `n` (M2, M3, m1). The dual
picture concerns `Bs(V(A,m))`, but route 5 has only `Z in Bs(V_H)`, so "ZE^all ⇒ ZE" is unproved and the right
global target is **ZE^all_H** (M4).

**Status: obstruction_found (partial), confirmed.** It is a proved barrier for one *tool* (degree-level zero
estimates in either natural framing), not for the open lemma. Two judge remarks (unrefereed, §2): the automorphism
data of M1 are system data re-encoded on one form (J-a, identity checked exactly), and every non-singleton
Philippon configuration built from route-5-certified points requires a **light endomorphism** of cost
`< 2 ord(c)/d` (J-b). So M1 does not reopen the torus picture except through the referee's open family (b).

**Score** (§5): ze 3. **Decision** (§7): stop internal rounds, as rounds 19, 20 and 22–25 already recommended; fold §6
into the expert note and put Questions A, B and the route-5 lemma ZE^all_H to experts.


## 1. Referee-corrected statement

### 1.1 ze (`ze.md`; `referee_ze.md`): survives in part; status `obstruction_found` confirmed as **partial**

Setting. `Z in T_M(Q)` is the cluster point, `d = n - 1`, `A` a `T = 1` slab box in `Z^n_0`, `V(A,m)` the forms on
`A` vanishing to order `m` at `e`, `V_H` the span of the integer `c in V(A,m')` with `|c|_inf <= H = U^2`. The cost
of a form along an integer endomorphism `M` of the cluster torus is `Lambda_M(c) = sum_p omega_p h_{supp c}(M beta_p)`
(`Lambda = Lambda_I`), and `lambda = lambda_min(G|_{1^perp})` for the projected Gram operator `G` of the profile.
Philippon's theorem (Bull. SMF 114 (1986) 355–383, bibliography verified) is used in the shape quoted from memory in
`ze.md` (2.1); neither author nor referee could reach numdam, eudml or arxiv from the session.

> **Round 26 ze (referee-corrected).** No growth improvement and no uniform bound; `M <= (2/3+eps) log R/loglog R`
> stands. ZE-col\*, MIN^pair and `(R_H)` are not derived.
>
> 1. *(proved)* **What route 5 certifies for one form.** Lemma 3.1 (Liouville plus forced divisibility, all Gaussian
>    gcd costs kept, common phase cancelled by degree-0 characters) certifies `P_c(Z^M) = 0`, to order decreasing
>    linearly, for every integer endomorphism `M` with `Lambda_M(c) < ord(c) - O(E)`. This includes
>    * the powers `Z^j`, `|j| < ord(c)/Lambda(c) <= 1.07 sqrt(n)/lambda`, by the new spectral LP cap
>      `ord(c) <= (sqrt(n)/lambda) Lambda(supp c)` (Prop 3.3, sharpening [D] Thm 2.3; `lambda <= n/(n-1)` always,
>      with equality for Sylvester profiles, so numbers at that value are best cases);
>    * the images `Z^{j sigma}` under profile automorphisms, at cost exactly `Lambda(c)` (Sylvester: `|AGL(k,2)|`,
>      `ln = 106.89` at `n = 4096`);
>    * possibly light endomorphisms `M = 1 m^T + u v^T` (profile-dependent; not analysed).
>
>    On `T = 1` regime boxes the powers guaranteed for every `c in V_H` are `j = +-1` at `ceil(4w)` and `|j| <= 2..4`
>    at `m_max` (`n <= 16384`).
> 2. *(proved, given the shape of Philippon's theorem; margins `>= e^{1.7e4}`, so the exact constants are
>    irrelevant)* **Torus picture is vacuous.** On `G_m^d` with data in `<Z>`, (2.1) is satisfied by `H = {e}` whenever
>    `lambda > 2.14 sqrt(n)/d` (e.g. `lambda > 0.04` at `n = 2966`). The generous variant GV without the sumset loss
>    is vacuous for every full-rank profile with `ln W <~ d ln(d/11)`, which contains the whole growth regime.
>    Automorphism points do not change this, because their `d`-th powers cost `d Lambda(c) > ord(c)`. Light
>    endomorphisms are not treated.
> 3. *(proved only in this restricted form; exact big-integer counts)* **Hilbert-level single-form barrier.** For every
>    `Y in G_m^d(C)` some nonzero `c in V(A,mu) (x) C` carries the vanishing data at the **powers** `Y^j` with
>    `ze.md`'s `c`-independent orders (Thm 3.6). This is checked at `n in {2966, 3000, 4096, 8192}` only:
>    * for all `mu <= m_max + 1` when `lambda = n/(n-1)`;
>    * at `mu = m_max` iff `lambda >= lambda*`, with `lambda* = 0.0168, 0.0198, 0.0094, 0.0061` (template-dependent;
>      for `lambda = 0` profiles the count fails at `m_max`);
>    * with the `AGL(k,2)` orbit of the Sylvester profile added, only for `mu <= 2559` (`n = 4096`) and `mu <= 5186`
>      (`n = 8192`); beyond, the linear-algebra argument fails (not disproved).
>
>    A barrier covering **all** route-5 data needs a per-form count over the certified set `D(c)`; it is open. Thm 3.6
>    says nothing about small-value estimates (the witness is complex, no height control). Margins are thin
>    (33–56 nats).
> 4. *(proved)* **Dual picture** (`G_a^d x G_m`, Lemma 4.1–Prop 4.5). The obstruction subgroups are `V x {1}` and
>    `V x G_m`; for nondegenerate points every class count is point-independent; the degree-level estimate is vacuous;
>    the Hilbert-sharp GV form is false (exact at `n = 2966, 4096`, tightest `ln` ratio `0.197`; the step "`>= L+1`
>    values on `A`" is repaired for sign columns with two 0's and two 1's off position `n`). **But** the dual picture
>    is about `Bs(V(A,m)) ⊂ Bs(V_H)`. Route 5 knows only `Z in Bs(V_H)` (from rigidity) and needs `Z notin Bs(V_H)`;
>    the record has `dim V_H >= |A|/2` only, against `dim V(A,m) ~ (3/4)|A|` at `m_max`. So ZE^all(A,m) and
>    ZE-tor^glob do **not** imply ZE unless `V_H (x) Q = V(A,m) (x) Q`, which is not on record. The relevant global
>    target is **ZE^all_H: `Bs(V_H)` contains no nondegenerate point**, which is stronger and not group-structured.
> 5. *(proved)* **The actual point.** `Z` on a full-rank subcluster lies on no proper algebraic subgroup (unique
>    factorisation plus affine independence; the bound `|<v,delta>| < 2pi` in `ze.md`'s proof is false and not
>    needed). Its relative distance to a subtorus `K_U` is `>= (2/(C'|u|_2)) e^{(W/2 - X_u)/2}` for each defining
>    character `u` (constant corrected from `4/...`). This is uniform in `R` but does not touch the obstruction,
>    which is a count at the trivial subgroup.
>
> **Verdict.** Both natural zero-estimate framings of route 5 are useless at the degree level, with the trivial
> subgroup as obstruction. The general "no single-form zero or multiplicity estimate can help" is proved only for
> powers-only data and `lambda >= lambda*`. The "short by `~ d!`" headline describes Philippon's degree normalisation
> (`d!(2T)^d` against `(2T+1)^d`), not an intrinsic deficit. Status **obstruction_found (partial)**; the open lemma
> is unchanged.

**Owner constraints.** One violation, corrected by the referee: a template-dependent threshold (`lambda >=
lambda*`) presented as uniform (M3). No fixed-`S` Subspace exceptional set is used; no Gaussian gcd or denominator
cost is discarded (Lemma 3.1 takes the minimum over `supp c` of exact valuations, which is a legitimate sharpening);
no drifting common phase is set to zero (degree-0 characters); Props 5.1–5.2 are uniform in `R` with explicit
constants, and nothing is claimed uniform that depends on the template except through stated `lambda` or `(n, r, T,
log 5)` dependence.


## 2. Status of route 5's zero estimate at the actual point after this round

**Status: open, unchanged in substance; one tool family closed; the global target restated.**

The lemma is still: for the cluster's own regime box, some form in the forced space `V_H` does not vanish at `Z`
(ZE-col\*; equivalently `(R_H)` / MIN^pair at the level of round 16). Round 26 changes the map in three ways.

1. **Closed (proved, degree level).** Philippon-type zero and multiplicity estimates, in either natural framing, give
   no information: the inequality to beat is a count at the trivial subgroup, and route 5's data sit below it. The
   actual point's nondegeneracy (§1.1 item 5) is irrelevant to that count. This joins round 16's proof that
   degree/variance/rank-blind anti-concentration tools are vacuous here, and round 11's refutations of ZE-braid, SL,
   SL\* through cheap rational subtori.
2. **Corrected target (proved implication structure, M4).** The global statement that a group-theoretic argument would
   have to deliver is **ZE^all_H(A,m): `Bs(V_H)` contains no nondegenerate point**, not ZE^all(A,m). Since `V_H` is
   defined by a height cut-off, ZE^all_H is not group-structured, and no zero estimate on a group takes it as input.
   ZE^all(A,m) would suffice only if `V_H (x) Q = V(A,m) (x) Q`; whether a height-`U^2` integer basis of `V(A,m)` exists
   is beyond the Liouville budget (referee §1.6) and is not on record.
3. **Open sub-questions isolated by the referee.** (a) A per-form count over the certified set `D(c)` (all
   endomorphisms `M` with `Lambda_M(c) < ord(c) - E`), which would replace Thm 3.6's application. (b) Light
   endomorphisms, which are the only certified points not yet excluded from the torus picture (J-b below).

*Judge remark J-a (unrefereed; the identity is checked exactly in `judge_checks/j1_aut.py`).* For a coordinate
permutation `sigma`, `P_c(Z^{j sigma}) = P_{sigma.c}(Z^j)`, where `sigma.c` permutes exponent vectors. If `A` and the
order condition at `e` are `sigma`-invariant (true for profile automorphisms that preserve the slab box; I have not
checked the slab cut itself), `sigma.c in V(A,mu)` with the same `Lambda` and `ord`. So the automorphism data of M1 say
only that `Z^j in Bs(span_sigma {sigma.c})`: they are **route 5's own system data at `Z`**, restricted to the
`Aut`-orbit span of one form, re-encoded as extra points for that form. This supports `ze.md`'s moral (single-form
estimates lose the overdetermination) and explains M2: once orbit points are added, the "single-form" count becomes a
small-system count and the overdetermination reappears at `mu` near `m_max`. Orbits under `Aut(G_m^d)` are not
translates, so they are not the `Sigma` of a Philippon estimate except through the sumsets treated in J-b.

*Judge remark J-b (unrefereed; proved here by one triangle inequality).* Let `Sigma` be any finite set of
route-5-certified points for `P_c`, possibly after translating `P_c` (a zero estimate may be applied to any translate).
Translate so that `e in Sigma'` and the translate is by `t = Z^{M_0}`, `t` certified, so `Lambda_{M_0}(c) < ord(c)`.
For `s = Z^{M_1} in Sigma'`, Philippon needs vanishing at `t s^d`, certified only if
`Lambda_{M_0 + d M_1}(c) < ord(c) - E`. Since `M -> Lambda_M(c)` is a seminorm,
`Lambda_{M_0+dM_1}(c) >= d Lambda_{M_1}(c) - Lambda_{M_0}(c) > d Lambda_{M_1}(c) - ord(c)`. Hence a non-singleton
admissible `Sigma` needs an endomorphism `M_1` with `Z^{M_1} != e` and
`Lambda_{M_1}(c) < 2 ord(c)/d <= 4T`. Powers and automorphisms cost `|j| Lambda(c) >= lambda ord(c)/sqrt(n)` by Prop 3.3,
which exceeds `2 ord(c)/d` once `lambda > 2 sqrt(n)/d`; differences `sigma_1 - sigma_0` of automorphisms are not
covered by Prop 3.3. So, beyond the referee's `d`-th power argument (which assumes `e in Sigma` without translation),
**the torus picture can reopen only through endomorphisms of cost `< 2 ord(c)/d`**, i.e. directions nearly invisible
to `supp c`. These are the referee's family (b) and, at `Lambda_M(c) = 0`, the projection-cheap subtori that refuted
SL/SL\* (round 11). Nothing here suggests they help: cheap directions are exactly where fakes and the actual point
look alike.


## 3. Consequences for growth and for the uniform bound

Only growth (`M = o(log R/loglog R)` for every `R`) or uniform progress (`M(C) < infinity`) counts. **Neither moved.**

* **Growth.** No route shortened. Route 5 still needs ZE at the actual point; by round 11, ZE-col\* would give
  `M = O((loglog R)^2)`, and by round 16, `(R_H)` would give MIN^pair and `M = O(loglog R)` on all norms (round 10
  Thm 5.7 on squarefree norms under ZE). Round 26 removes one candidate tool for that lemma (degree-level zero /
  multiplicity estimates) and corrects the global target to ZE^all_H. The `(2/3 + eps)` bound stands.
* **Uniform bound.** Untouched. Route 5 is a growth route (its LP supercriticality factor
  `0.12 sqrt((n-1)/ln 4r)` needs `n` large), and Question B (`(L_5^circ)`, or `A_5(delta)` with `delta < 1/30`)
  is not addressed.
* **Map value.** A proved negative answer to the question "can a Philippon-type zero estimate prove ZE?", which an
  expert would ask first. The answer is robust to the unverified constants of Philippon's theorem (margins
  `>= e^{1.7e4}`). The Hilbert-level part is thin and template-dependent and should be cited only in its corrected
  form (§1.1 item 3).


## 4. Judge's own checks (unrefereed)

All runs use `prlimit --as=3500000000`, one process at a time. The scripts are written from scratch and import nothing
from the author's or referee's code.

| id | script → output | claim checked | result |
| --- | --- | --- | --- |
| J1 | `judge_checks/j1_aut.py` → `j1_aut.out` (exact integers / Gaussian rationals) | referee M1(a): every `x -> Ax+b` maps each Sylvester column to `+-` a column; `Lambda` is invariant; identity `P_c(Z^sigma) = P_{sigma.c}(Z)` (J-a); `ln|AGL(k,2)|` | `k = 3`: all 1344 affine maps OK; `k = 4`: 3000 random affine maps OK; `Lambda` invariance 0/400 failures; identity 0/60 failures; `ln|AGL(12,2)| = 106.89`, `ln|AGL(13,2)| = 124.91` |
| J2 | `judge_checks/j2_orbit_count.py` → `j2_orbit_count.out` (exact big integers for the holds/fails verdicts) | referee M2 at `n = 4096`: orbit-inclusive count `|AGL| x sum binom(tau_j+d-1,d)` against `ceil(3^d/2) - binom(mu+d-1,d)` | `J = 68`, `m_max = 2631`; orbit count holds at `mu = 2559` (`4497.88 < 4498.12`), fails at `2560` (`4498.84 > 4498.12`) and at `m_max` (`4565.80 > 4498.00`); powers-only holds at all five `mu` tested. Agrees with the referee |

J-a's consequence ("automorphism data are system data") and J-b are logical remarks of mine, **unrefereed**; J-b's
proof is the one-line seminorm inequality written in §2, parked with its scale comparison in
`judge_checks/judge_notes.md`. I did not re-derive Prop 4.5 or §5; the referee's independent reproduction stands.


## 5. Score

Rubric of `round10/README.md` §1: (a) correctness after refereeing, (b) fake-violated input, (c) movement of growth
or uniformity. The overall score is a judgement that weighs (b) and (c) above (a).

| write-up | (a) | (b) | (c) | overall | reasons |
| --- | --- | --- | --- | --- | --- |
| **ze** | 5 | 3 | 2 | **3** | **Correct core.** Lemma 3.1 (all gcd costs kept), the spectral cap Prop 3.3, Cor 3.2/3.4, the degree-level vacuity in both pictures (Props 3.5, 4.4), Lemma 4.2/Prop 4.3 and Prop 5.1's conclusion were re-derived by the referee, with every number reproduced; J2 reproduces the referee's orbit count. **Deductions.** Four majors: a false "only the powers" claim (M1); the universal single-form barrier not established once certified symmetric data are counted (M2); a template-dependent threshold presented as uniform, an owner-constraint violation (M3); a broken implication "ZE^all ⇒ ZE" in the map statement (M4). Minor overclaims (small-value estimates, "`n = 2966..8192`", "short by `~ d!`") and two slips in §5. **(b) = 3.** The input is route 5's forced-divisibility system, whose exact relations almost every fake violates; the residual (ZE at `Z`) is the open problem, and nothing new is added on the fake side. **(c) = 2.** A proved barrier for a natural tool (degree-level zero estimates, both framings), robust to the unverified constants, and a proved correction of the global target to ZE^all_H. No growth or uniform movement. |

The work answers honestly the question that rounds 11–25 left as the record's one internal decisive item, but only
in the negative and only for one tool family. As in rounds 22–25, no movement of growth or uniformity caps the overall
score at 3.


## 6. Proposed addition to the expert note (text only)

This is a text proposal only; `expert_note/expert_note.md` was not edited. It uses the note's labels.

**Add to §5, item 3 ("Auxiliary polynomials and single-place Diophantine approximation"), a sixth bullet:**

> * *Zero estimates do not supply route 5's missing lemma* (round 26 `ze.md`, `referee_ze.md`; **proved** at the
>   degree level, given the shape of Philippon's 1986 zero lemma, with margins $\ge e^{1.7\cdot10^4}$ that make its
>   exact constants irrelevant). Route 5's Liouville step certifies, for one form $P_c$, vanishing at $Z^M$ only for
>   integer endomorphisms $M$ with $\Lambda_M(c)<\mathrm{ord}(c)-O(E)$; for the powers $Z^j$ this means
>   $|j|<\mathrm{ord}(c)/\Lambda(c)\le1.07\sqrt n/\lambda$, by the spectral cap
>   $\mathrm{ord}(c)\le(\sqrt n/\lambda)\Lambda(\mathrm{supp}\,c)$, and only $|j|\le1$–$4$ on regime boxes. A
>   Philippon estimate on $\mathbb G_m^d$ needs progressions of length about $d$, so its conclusion is met by the
>   trivial subgroup. The dual framing (one exponential polynomial on $\mathbb G_a^d\times\mathbb G_m$ vanishing on
>   the slab box) has obstruction subgroups whose counts do not depend on the point, is vacuous at the degree level,
>   and is false in Hilbert-sharp form. Moreover it concerns $\mathrm{Bs}(V(A,m))$, while route 5 needs
>   $Z\notin\mathrm{Bs}(V_H)$ for the height-cut forced space $V_H$. The global statement any group-theoretic
>   argument would have to prove is therefore **ZE$^{\rm all}_H$: $\mathrm{Bs}(V_H)$ contains no nondegenerate point**,
>   which is not group-structured. The Hilbert-level claim that no single-form estimate can help is proved only for
>   powers-only data and well-conditioned profiles ($\lambda\ge\lambda^*(n,\mu)$, checked at four values of $n$).
>   Profile automorphisms are certified too, and once they are counted the argument fails near $m_{\max}$ (the
>   claim is then unproved, not refuted), so it should not be cited as general.
>   The actual point lies on no proper algebraic subgroup, but that does not affect a count at the trivial subgroup.

**Add to the pointers table of §6:**

> | round 26: zero estimates vs route 5's ZE; spectral LP cap; target ZE$^{\rm all}_H$ | `round26/README.md`, `round26/ze.md` (with `referee_ze.md`) |


## 7. Honest next step

**Decision: stop internal rounds.** This confirms rounds 19, 20 and 22–25.
* Round 25 named route 5's zero estimate at the actual point as the only internal item on which several routes
  converge. Round 26 tested the most natural tool for it and proved that tool vacuous at the degree level.
* What remains is ZE^all_H (or its per-cluster form ZE-col\* / `(R_H)`), a statement about the base locus of a
  height-cut integer linear system. Rounds 11, 16 and 26 have now shown that cheap-subtorus counterexamples, rank-
  and degree-blind anti-concentration, and group zero estimates all fail to reach it.
* No new input that actual circles satisfy and fakes provably violate was produced.

**The single most promising next step is to fold §6 into the expert note and put Questions A and B, together with
ZE^all_H as a stand-alone base-locus question, to experts in Diophantine geometry and in zero estimates.** The
expert question can be stated without the circle problem: for the forced-divisibility system of a regime box, does
the span of its integer forms of height `<= U^2` have a nondegenerate common zero near `e`? An expert might know a
height-sensitive base-locus tool (arithmetic Bézout / Arakelov-type Hilbert–Samuel for small-height sections) that
this record has not tried. I cannot judge from inside the record whether such a tool exists.

*If one more internal round is run anyway*, the only concrete item left by this round is narrow and is more likely
to yield a further barrier than a proof: classify the **light endomorphisms** (cost `Lambda_M(c) < 2 ord(c)/d`,
J-b) for regime profiles, and do the per-form count over `D(c)` that would turn Thm 3.6 into a full single-form
barrier. A positive outcome would only close the torus picture completely. It would not prove ZE.

*Honest probabilities (judge's estimates; not derived).*

| outcome of one more internal round on ZE / ZE^all_H | probability |
| --- | --- |
| a proved map correction or barrier (e.g. light endomorphisms classified, per-form `D(c)` count) | about 20% |
| ZE-col\*, `(R_H)` or ZE^all_H proved at the actual point | about 1% |
| unconditional growth `M = o(log R/loglog R)` for every `R` | < 0.5% |
| the uniform bound, `(L_5^circ)` or `A_5(delta < 1/30)` | < 0.1% |

No probability is offered for expert circulation, for the reason given in round 25 §8.


## 8. Files

* Write-up: `ze.md`; notes `notes.md`; checks `ze_checks/` (`powers_valuation`, `spectral_cap`, `vacuity_numbers`,
  `single_poly_barrier`, `hs_refute`, each `.py` with `_out.txt`).
* Referee report: `referee_ze.md`; checks `ref_ze/` (`recompute.py` → `recompute_out.txt`, `hs_check.py` →
  `hs_check_out.txt`, notes `ref_ze/notes.md`).
* Judge: this `README.md`; `judge_checks/j1_aut.py` → `j1_aut.out`, `judge_checks/j2_orbit_count.py` →
  `j2_orbit_count.out`, derivation notes `judge_checks/judge_notes.md`.

