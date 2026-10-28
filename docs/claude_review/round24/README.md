# Round 24 — judge report (does the anisotropy of T_5 give anything over the 15/4 cap?)

| key | write-up | referee report | lens |
| --- | --- | --- | --- |
| theory | `theory.md` (`theory_checks/`, `notes.md`) | `referee_theory.md` (`referee_theory_checks/`, `referee_notes.md`) | can the anisotropy of `T_5` (isotypic sections, conjugate products/Thue, point parity, Hilbert 90) give a section/Liouville inequality beating `15/4`? |
| search | `search.md` (`search_checks/`, `search_notes.md`) | `referee_search.md` (`referee_search_checks/`, `referee_search_notes.md`) | the round-23 twisted Pell ansatz: balanced rational curves through `e` with `B.C = 30`, two conjugate real-quadratic branches of contact 4, `alpha = 15/4` |

All paths are relative to `scratchpad/round24/`. The repository was read only. Labels follow the overview:
proved / conditional / exact / numerical / evidence / heuristic / unrefereed. §1 gives the **referee-corrected**
statements. Remarks marked *judge* are mine and are **unrefereed**; my computations are in `judge_checks/` (§5).

## 0. Bottom line

**No growth improvement and no uniform bound.** `M <= (2/3 + eps) log R/loglog R` is still the best bound for every
`R`. Neither write-up reaches `A_5(delta)` for any `delta < 1/30`, `alpha^ess_e(-K) > 15/4`, or `(L_5^circ)`. Both
survive refereeing in their mathematics. Only growth or uniform progress counts, and there is none.

**Does the anisotropy give any gain over the 15/4 cap? No — nothing found this round.**
* *Pointwise/section side (theory, refereed, status `reduces_to_split_cap` confirmed).* Every argument of the class
  `L` — Liouville with toric forced divisibility, enriched by `sigma`/`iota`-isotypic sections, conjugate products and
  group words (Thue-type), point parity `v_pibar = -v_pi`, reality of values, and Hilbert-90 binary forms — is implied
  by toric Liouville for a single section (or its `n`-point version) and so is capped at the **split** caps: `15/4`
  (`delta = 1/30`) on balanced simple profiles, `delta = 1/10` for arbitrary profiles, for every family with
  `log H(F) + sum d_i log d_i = o(sum d_i w_i)` (in particular every fixed finite family). The twist fixes `e` and the
  orders at `e`, so Galois averaging loses nothing; the anisotropy changes **which points exist**, not which section
  inequalities hold.
* *Curve side (search, refereed, status downgraded from `anisotropy_gain` to `map value`).* Balanced twisted curves at
  exactly `15/4` **exist over `R`**: an explicit real one-parameter family `M` (Theorem A). **Only finitely many members
  of its complex component are defined over `Q`** (Theorem B, Faltings, ineffective; none of height `<= 300`). But the
  same finiteness holds for the split `Q`-form of the same complex component, so Theorem B is **not** an anisotropy
  effect. The only twist-specific input on record remains the round-23 **parity lemma** (`B.C >= 30`, `mult_e >= 8`).

**What changed for Question B.** Nothing in its truth value. Map value only:
1. a proved, scoped barrier: anisotropy-aware *pointwise* arguments of Liouville type are exactly as weak as split
   ones (§1.1);
2. the first balanced twisted `15/4` curves (over `R`), whose `Q`-members are finite in number, so `M` forces at most
   finitely many curves into `Z` and is useless for refuting Question B (§3);
3. a methodological correction: every blind numerical solver missed `M`, so blind searches in this class are **not**
   evidence of emptiness (this retroactively weakens round 22 `toric.md` §7's "no nondegenerate solution").

**Scores** (§6): theory 4, search 3. **Decision** (§8): stop internal rounds (confirming rounds 19, 20, 22, 23).
Fold §7 into the expert note and put Question B to experts.

## 1. Referee-corrected statements

### 1.1 theory (`theory.md`; `referee_theory.md`): survives with corrections; status `reduces_to_split_cap` confirmed

Setting: `X = Lbar_5^tw`, `sigma_X = iota∘conj` (Galois acts by `-1` on cocharacters, `1_A <-> 1_{A^c}`), `e` the
identity, A-points `y` with toric scales `w`, archimedean depth `lambda(y) = -log dist(y,e)`, `mu(y) = lambda(y)/w`
(balanced Liouville gives `mu <= 4 + O(eta)`; `(L_5^circ)` is `mu <= 4 - c'` off `Z` on balanced points). `(TL)` is toric Liouville: archimedean Taylor upper bound `>=` product-formula lower bound,
for one section `f` with `f(y) != 0`.

> **Theorem (round 24 theory, referee-corrected; proved except where labelled).** Let `L` be the class of arguments
> "archimedean Taylor upper bound `>=` product-formula lower bound", where the lower bound may use integrality, toric
> forced divisibility at every place of `Q(i)`, point parity and reality of values, and the auxiliary object is
> (a) `f in Z[i][M]` of any `sigma_X`/`iota`-isotypic type; (b) `F in Z[i][M^n]` evaluated at group words (including
> `iota`) in A-points `y_1..y_r`; or (c) any `G in Z[i][P, Pbar]` evaluated at the Hilbert-90 vectors
> (`x_j = P_j/Pbar_j`) of **any** representative of `y`.
> 1. *Reduction.* (a) is exactly `(TL)` for the `Q`-rational centrally symmetric section `f f^*`. (c) is **implied by**
>    `(TL)` for one multihomogeneous component `h_{D,k} in K[S_m]` with `h_{D,k}(y) != 0` (Lemma R2). (b) is implied, at
>    configurations with `lambda(w_k) <= min_{letters} lambda(y_i) + O(1)`, by the `n`-point inequality `(*)` for
>    `phi^*F` at the generating points (Lemma R1, Prop 4.1).
> 2. *Balanced cap.* At `eta`-balanced A-points with **simple** profiles and pairwise disjoint primes (lacunary scales
>    in case (b)), every inequality in `L` holds whenever `mu(y_i) <= 4(1-eta) - C(F)/w_min`, with
>    `C(F) << log H(F) + sum d_i log(2+d_i)`. Hence no argument in `L` whose family satisfies
>    `log H(F) + sum d_i log d_i = o(sum d_i w_i)` proves `(L_5^circ)` or `A_5(delta)` for any `delta < 1/30`; in
>    particular no fixed finite family does.
> 3. *Arbitrary profiles.* Under the same proviso nothing in `L` gives `delta < 1/10`; the witnessing profile depends
>    on the family (not uniform in it).
> 4. *Form-independence.* The same numbers hold for the split form. For every `m`, the `iota`-isotypic part of
>    `H^0(m(-K))` with the parity of `t_max(m) = floor(10m/3)` reaches `t_max(m)`, the other part `t_max(m) - 1`
>    (proved for all `m`; exact for `m <= 7`: `(2,3),(6,5),(10,9),(12,13),(16,15),(20,19),(22,23)`).
> 5. *Not covered.* `(*)` on shared-prime configurations of genuine A-points (open; not needed for 2). Outside `L`:
>    global inputs (curve concentration + parity lemma; McKinnon/Vojta-type, conditional) and twist-sensitive
>    congruence divisibility at inert `q ≡ 3 (4)` and at `(1+i)` (twisted points reduce into
>    `T_5(F_q) = mu_{q+1}^5/mu_{q+1}`), whose effect is `O_F(1)` for fixed families: `c_q(f) = 0` for `q >= 2m`
>    (proved by character independence plus weak approximation; the `O(m log m)` bound for growing `m` is a sketch).

*Corrections adopted.*
* **M1 (major, repaired).** Theorem 7.1 claimed a cap for group words via Prop 4.3, which needs pairwise disjoint
  primes — never true for a tuple containing a word and one of its letters (`(y, iota y)`, `(y_1, y_1 y_2)`). Repaired
  by Lemma R1 (pull back along `phi` to the generating points; `width_{Q_F}(phi_* v) = width_{phi^*Q_F}(v)`,
  `Q_G ⊆ phi^*Q_F`). Lacunarity is a choice of *hypothetical* configuration, which suffices for a barrier.
* **M2 (major, repaired; owner constraint).** Prop 6.2's "fix `P_1 = 1`" for a non-homogeneous `G` evaluates at one
  representative only; at any other representative `x = rho x^(1)` the value `sum_k rho^k g_k` depends on a **drifting
  common phase**, so fixing `P_1 = 1` replaces that phase by zero — exactly what the owner forbids, and the author's §8
  audit asserted the opposite. Lemma R2 (split into multihomogeneous components) removes the dependence;
  "identical" becomes "implied by". Conclusion unchanged.
* **m1–m6 (minor).** `O_F(1)` thresholds are template-dependent (now explicit in item 2); the "exhaustive" list of
  inputs outside `L` omitted inert/ramified congruence divisibility (item 5); Prop 5.3's prose "parity can only weaken"
  is false section by section (`f = x_1/x_2 - 1`, `v = 1_{2}`: split 1, twisted 1/2; twisted strictly smaller in 8285
  of 20000 random cases) and is true only against the better of `f` and `iota^*f`; "simple" missing in Cor 3.3; the
  distance formula after Prop 2.2 holds only up to `≍`; "every argument equals its Galois average" is exact only for
  (a) and the pair `(y, iota y)`, otherwise "is implied by".
* *Also proved/exact and useful.* Hilbert-90 model: `Phi_0 = ab` (split) vs `a^2+b^2` (twisted); the bracket identity
  `t_ab^2 = [P_a,P_b]^2/(prod p^{2c_ab} 2^{e_a+e_b})` explains the record's `|t_ab| in sqrt(odd/2)` across unit
  classes; new exact fact (referee): `e_a = e_b` iff the units agree up to sign. Picard ranks 15 (twisted) and 26
  (split). Prop 4.2 (index additivity at the rational target) has a general three-line proof.

### 1.2 search (`search.md`; `referee_search.md`): mathematics survives; headline and status do not

> **Round 24 search (referee-corrected).** Consider the **ray-generic** minimal balanced class of rational curves
> `C` through `e` on `Lbar_5^tw`: `E_A.C = 1` for all 30 rays (so every cut has weight 2 and `B.C = 30`), 30 distinct
> places, two conjugate real-quadratic branches of contact `mu_0 = 4`, so `alpha_C = 15/4`.
> 1. **Normal form (proved).** For a `Q`-curve of this class all 30 places lie in `P^1(Q(i))`, and
>    `u_j = gamma_j P_j/conj(gamma_j P_j)` with `Im(gamma_j P_j) = kappa_j (t^2-m)^4`, `kappa_j in Q^*`; equivalently
>    `e_k(R_j)` is real for `k != 4`; the class is V-null with `mu_0 = 4` exactly. The split analogue (with `Q` in
>    place of `Q(i)`) holds by the same proof.
> 2. **Theorem A (proved; component claim numerical).** Put `c = m^3-11m^2+39m-13`,
>    `u = -(m+3)c/(16(m-3)^2(m-1)^2)`, `v = -(m+1)c/(4(m-3)^3(m-1))`; let `b_r` be the roots of
>    `W^4 + vW^2 - v^2W + uv^2`, `t^2 = prod b_r`, base points `p_r = t/b_r + i b_r` and `-conj(p_r)`. For **every**
>    `m in (-3,-1)` this gives a real curve of the class (the referee supplied the whole-interval argument; 399-point
>    grid). It is the only nondegenerate family in the reflection-symmetric x-structured ansatz. The general system has
>    Jacobian rank 28 along `M` (80 digits), so `M` is a component of dimension 1 modulo gauge, sweeping a real surface.
> 3. **Theorem B (proved, ineffective).** Only finitely many curves of the complex component `M_C` (hence of `M` and its
>    `S_5`-images) are defined over `Q`. Proof: a `Q`-member forces the scale-invariants `Y_r = b_r^2/(-e_2)` — the roots
>    of `G_m(Y) = (Y^2-Y+u)^2 + vY` — to be rational with `m in Q` (step 2 repaired by the referee to cover complex
>    points of `M_C`; `m -> (u,v)` birational over `Q`, checked by the referee); `disc_Y G_m` has exactly 8 odd-order
>    points, the geometric monodromy is `S_4`, the Galois closure has genus `>= 3`, and Faltings applies. No
>    nondegenerate `m = p/q` with `|p|, q <= 300` (only the degenerate `m = -3`, `b = 0`).
>    **The same finiteness holds for the split `Q`-form of `M_C`**, so Theorem B is not an anisotropy effect.
> 4. **Not proved.** Finiteness for the whole class (other components may exist: Q24.1); the expected dimension 1 is
>    heuristic; existence of even one `Q`-member of `M` is open (Q24.2).
>
> **Status: map value** (one component of the twisted Pell-type class eliminated over `Q`, in every form). Not an
> anisotropy gain.

*Corrections adopted.*
* **Major 1 (scope).** "Only finitely many of the class's members are defined over `Q`" is proved for `M_C` and its
  `S_5`-images only.
* **Major 2 (attribution).** Both ingredients of Theorem B are form-independent: ray-generic places are rational over
  the splitting field in every `Q`-form (in the split form they are even `Q`-rational), and the `S_4` monodromy lives
  on `M_C`, shared by both forms (`u = e_4/e_2^2`, `v = e_3^2/e_2^3` are invariant under every common scaling
  `b -> lambda b`, including `b -> i b`). The contrast "split forces `W`, twisted only finitely many" compares
  `B.C = 15` (split Thm 2) with `B.C = 30` (`M`); the step separating them is the round-23 parity lemma alone.
* **Minor.** Step 2 of Theorem B was written for real points only (repaired); birationality of `m -> (u,v)` was
  unchecked (now checked); "for every `m`" rested on 8 samples (now proved); the rational search silently skipped
  `m = ±1, ±3` (`m = -3` gives the degenerate `b = 0`); discriminant exponents depend on normalisation (parities
  confirmed); "`b_r` rational" holds up to a common factor; the class is the ray-generic one; the two-branch covering
  bound is cited from a source that did not re-derive it; "the other branch is `m -> -m`" has no counterpart.
* **Methodological correction (author, endorsed).** All blind solvers (LM in five parametrisations, 14 symmetric
  ansätze, complexified Newton, staging, lifting, 165 homotopy paths) missed `M`; the author withdrew the provisional
  "evidence of emptiness" reading. Blind searches here are not evidence of non-existence.

## 2. Does the anisotropy give any gain over the 15/4 cap?

**Answer: no gain is known, and the two natural places to look are now closed or neutral.**

| where the anisotropy could enter | finding | label |
| --- | --- | --- |
| isotypic / `Q`-rational sections, reality of values | same maximal orders as the split form (`t_max`, `t_max - 1`); `f -> f f^*` is lossless | proved (theory Prop 3.1–3.2, Thm 7.1) |
| conjugate products, `(y, iota y)`, group words (Thue-type) | single sections of `(d_1+d_2)(-K)` or `n`-point `(*)` at generating points; capped on disjoint-prime balanced data | proved (Props 4.1–4.3, Lemma R1) |
| point parity `v_pibar = -v_pi` | is the Gaussian profile itself (fakes satisfy it); the twisted right side is the mean of split right sides at `±v`, so it never beats the better of `f`, `iota^*f` | proved (Prop 5.1, 5.3; Cor 5.2 record-conditional) |
| Hilbert 90 / binary forms `a^2+b^2` vs `ab` | brackets reproduce exactly the record's pair bound; polynomials in `(P, Pbar)` reduce to sections (Lemma R2) | proved |
| congruence divisibility at inert `q ≡ 3 (4)`, `(1+i)` | twist-sensitive (`T_5(F_q) = mu_{q+1}^4`), but `c_q(f) = 0` for `q >= 2m`: `O_F(1)` for fixed families | proved (referee, m2); growing-`m` bound a sketch |
| curves over `Q` through `e` (parity lemma) | `E_A.C = E_{A^c}.C`, so balanced twisted `15/4` curves have `B.C >= 30`, `mult_e >= 8` | proved (round 23) |
| the `B.C = 30` Pell component `M_C` | finitely many `Q`-members — **in both forms** | proved (search Thm B + referee §5.2) |
| global concentration of A-points on curves (McKinnon/Vojta-type) + parity | would let parity sharpen degrees; but gauss Thm 4.3(ii) already finishes once concentration on bounded-degree curves is known | conditional / open |

**Judgement.** Every pointwise (local, section-theoretic) argument is form-blind for a structural reason: `sigma_X`
fixes `e`, preserves `ord_e`, and maps the product formula to itself, so averaging over `Gal(Q(i)/Q)` turns any
twisted argument into a symmetric split one with the same constants. The anisotropy is visible only to **global**
objects: which rational curves exist (parity), and which points exist at all. The one global object tested this round
(`M_C`) turned out to be killed over `Q` by its own monodromy, not by the twist. So the calibration from round 23 is
sharpened, not changed: a proof of Question B for `T_5` alone that is not also a proof for the split form (hence of
Erdős #887, absolute `K`) would have to use a **global** twist-sensitive input — at present only the parity lemma, whose
only known use is through curve concentration, which is McKinnon/Vojta-strength.

*Judge remark (unrefereed, logical).* The referee's form-independence of Theorem B is the algebraic statement that
`(u, v)` and `Y_r` are invariant under every common scaling of the `b_r` (checked exactly, J1). Twisted `Q`-members
need `b_r in lambda·Q` with `lambda` such that the places are conjugate pairs; split `Q`-members need
`b_r in (lambda/i)·Q`. Both reduce to "`G_m` splits completely over `Q` at rational `m`", which is the same condition.

## 3. Balanced twisted curves at 15/4: status, the exceptional set Z, refuting Question B

**Status (refereed unless marked).**

| question | answer |
| --- | --- |
| Do balanced twisted curves at `alpha = 15/4` exist over `R`? | **Yes** (Theorem A): the real one-parameter family `M`, `m in (-3,-1)`, of the ray-generic class `B.C = 30`, two conjugate branches of contact 4. This is the first explicit object of this kind; round 22's searches and the round-24 blind solvers all missed it. |
| Over `Q`? | **Unknown whether any exists.** From `M_C`: at most finitely many (Theorem B, ineffective); none with `m` of height `<= 300`. A `Q`-member also needs the conic to have a `Q`-point and `-1/sigma` to be a non-square (else the branches are rational and `alpha = 15/2`). Other components of the class: open (Q24.1). |
| Is the finiteness an anisotropy effect? | **No.** The same proof gives finiteness for the split `Q`-form of `M_C` (referee §5.2; judge J1). |
| What is twist-specific? | Only the parity lemma: it forbids the split `B.C = 15` Pell hypersurface `W` and forces `B.C >= 30`, `mult_e >= 8`. |

**Consequences for the exceptional set `Z`.**
* A `Q`-curve of the class with irrational (real-quadratic) branches carries infinitely many A-points at `alpha = 15/4`
  (Pell convergents), hence A-points for every `delta < 1/30`. So **every admissible `Z`** in `A_5(delta < 1/30)` or
  `(L_5^circ)` must contain each such curve. From `M_C` this is at most finitely many curves — harmless.
* *Contrast with the split form:* there, every admissible split `Z` must contain the hypersurface `W` (round 23).
  On the twisted form nothing of dimension `>= 2` is known to be forced into `Z`. That is the only (weak) sense in
  which the twisted form currently looks "better" than the split one, and it rests on the parity lemma alone.
* A component of the class whose place-curve has genus `<= 1` with infinitely many `Q`-points would force a
  **surface** into `Z` (Q24.1). Still harmless for Question B as stated (`Z` may be any proper closed set).

**Consequences for refuting Question B.**
* *Judge (unrefereed, logical consequence of refereed statements).* No fixed curve class can refute Question B. By
  gauss Thm 4.3(ii) (killer-checked; `k >= 5`, `delta > 0`), for each `d` the curves with `B.C <= d` carrying infinitely
  many A-points lie in a proper closed subset of `T_5`. Curves at exactly `15/4` give A-points for every
  `delta in (0, 1/30)`, so for each fixed `delta` they are covered by that statement. Hence a refutation by curves
  needs `Q`-curves of **unbounded** `B.C` whose A-points are Zariski-dense; the `B.C = 30` class (and any fixed higher
  class, e.g. `B.C = 60` of Q24.3) can only enlarge `Z`. This is the same conclusion as round 23's "a refutation needs
  infinitely many `W_i`", now for the twisted form.
* The absolute form (round 23 J4, refereed): a refutation could also come from `M_inf(C) -> inf`. Nothing this round
  bears on it.
* **Net:** balanced twisted `15/4` curves are now known to exist over `R` and to be scarce over `Q` on the one known
  component. This neither proves nor refutes Question B; it removes one candidate source of Zariski-dense A-points and
  confirms that the A-points at `15/4`, as far as they are known, sit on finitely many curves (twisted) versus a
  hypersurface (split).

## 4. Growth and the uniform bound

**Nothing this round moves either.**
* *Growth `M = o(log R/loglog R)`.* Neither write-up contains a growth statement. Theory is a cap on a class of
  methods; search is a statement about finitely many curves. `(2/3 + eps)` stands.
* *Uniform bound.* The routes remain `alpha^ess_e(-K) > 15/4` on balanced one-unit-class sequences (`(L_5^circ)`) or
  `A_5(delta < 1/30)` (which gives even the absolute form), then gauss Thm 3.4′. Both open and Vojta-strength. This
  round adds:
  * a **proved barrier**: anisotropy-aware pointwise Liouville arguments (class `L`) cannot reach them;
  * a **calibration**: the anisotropy is invisible to every pointwise argument and so far visible only through curves
    (parity), so a `T_5`-only proof needs a global twist-sensitive input;
  * **no** new conditional route.
* *Fakes.* No input this round is fake-violated. Theory Prop 5.1/Cor 5.2: point parity is the Gaussian profile, so
  fakes satisfy it (record-conditional on the killer-checked fakes). The search concerns curves over `Q`, which do not
  distinguish fakes.

## 5. Judge checks (unrefereed)

Scripts and outputs are in `judge_checks/`, run under `prlimit --as=3500000000`, one process at a time, exact sympy.
**None of this has been refereed.**

**J1. Theorem B ingredients and form-independence, independent code** (`j1_thmB.py` → `j1_thmB.out`,
`j1b_inf.py` → `j1b_inf.out`).
* `Res_W(b-quartic, Y + W^2/v)` is a `Y`-free multiple of `G_m(Y) = (Y^2-Y+u)^2 + vY` (exact).
* For the `b`-quartic, `e_1 = 0`, `e_2 = v`, `e_3 = v^2`, `e_4 = u v^2`, so `u = e_4/e_2^2`, `v = e_3^2/e_2^3`; both are
  invariant under `b -> lambda b` for symbolic `lambda` (in particular `lambda = i`). This is the algebraic core of the
  referee's form-independence claim.
* `disc_Y` of the cleared `G_m` (leading coefficient `256(m-3)^4(m-1)^4`, all coefficients of degree 8 in `m`) is
  `-2^47 (m-3)^12 (m-1)^16 (m+1)^2 (m+3) c^3 (m^4+70m^2-48m+9)`, identical to the referee's normalisation. Odd-order
  finite points: roots of `(m+3)·c·(m^4+70m^2-48m+9)`, **8 distinct points** (squarefree product; coprime to the leading
  coefficient). Order at `m = infinity`: `6·8 - 44 = 4`, even.
* `Gal(G_m/Q) = S_4` at `m = -2, -3/2, 5, 7/3, -5/4` (sympy `galois_group`).
* *Consequence (judge, unrefereed, elementary).* With geometric monodromy `S_4` (normal in the arithmetic group, contains
  odd elements), the Galois closure is an irreducible degree-24 cover of `P^1_m`, and each of the 8 odd points has
  ramification index `>= 2` above it, contributing `>= 12` to Riemann–Hurwitz: `2g - 2 >= -48 + 96`, so `g >= 25`.
  This strengthens the stated `g >= 3`; it does not change Theorem B.

**J2. No other new claims.** The logical remarks of §2 (form-independence reduces to complete splitting of `G_m`) and
§3 (no fixed curve class refutes Question B, from gauss Thm 4.3(ii)) are deductions from refereed or killer-checked
statements and are labelled unrefereed where they appear.

## 6. Scores

Rubric of `round10/README.md` §1: (a) correctness after refereeing, (b) fake-violated input, (c) movement of growth or
uniformity, each out of 10; the overall score weighs (b) and (c) above (a).

| write-up | (a) | (b) | (c) | **overall** | reason |
| --- | --- | --- | --- | --- | --- |
| **theory** | 7 | 0 | 2 | **4** | **Strengths.** Answers the lens with a proved, scoped barrier: all four anisotropy candidates (isotypic sections, conjugate products/Thue, point parity, Hilbert 90) reduce to toric Liouville for one symmetric section, so the caps are the split caps. Clear structural reason (lossless Galois averaging). Isotypic orders proved for all `m` (exact `m <= 7` after the referee). The Hilbert-90 bracket formula explains the record's `sqrt(odd/2)` anomaly. Honest open-point list. **Deductions.** Two major gaps: group words claimed capped via a disjoint-primes proposition that never applies to them (M1), and a non-homogeneous `G` handled by fixing `P_1 = 1`, which sets a drifting phase to zero — an owner constraint — while the author's audit asserted compliance (M2). Template-dependent `O_F(1)` thresholds stated as uniform (m1); an "exhaustive" list that was not (m2); a false prose slogan (m3). All repaired without changing the conclusion. **(b) = 0:** point parity is fake-satisfied (author's own Cor 5.2). **(c) = 2:** a proved barrier-map item (anisotropy is pointwise-invisible), no movement of the target. |
| **search** | 6 | 0 | 1 | **3** | **Strengths.** Found what every blind search missed: the explicit real family `M` of balanced twisted `15/4` curves (Theorem A, now proved on the whole interval), with exact algebra, an `S_4`-invariant reduction to a rational quintic, and a clean Faltings finiteness theorem (Theorem B). Exemplary honesty on evidence: withdrew its own "evidence of emptiness" reading after finding `M`. **Deductions.** The headline did not survive: finiteness was stated for the class but proved for one component (major), and the status `anisotropy_gain` rested on an attribution the referee refuted — the same proof works for the split form (major). Step 2 covered only real points; birationality unchecked; the rational search silently skipped `m = ±1, ±3`. **(b) = 0.** **(c) = 1:** one component eliminated over `Q`; map value only, and by §3 no fixed class could have refuted Question B anyway. |

*Referee reports.* Both were thorough and independent, with new code for every number.
* referee_theory found both major gaps (one of them a constraint violation the author had certified as compliant),
  supplied Lemmas R1/R2, finished the author's unfinished computations (`m = 6, 7`), gave a general proof of Prop 4.2,
  and found the inert-place omission with a proof that it is `O_F(1)`.
* referee_search rebuilt `M` independently, proved the whole-interval claim, repaired step 2 for complex points,
  checked birationality and the `S_4` monodromy, extended the rational search to height 300, and — the decisive
  contribution — showed Theorem B is form-independent, downgrading the status.

For comparison: round 23 split 4 and width 4; round 22 toric 4 and lit 3; round 20 gauss 5.

## 7. Proposed additions to expert_note.md §3 (text only)

Proposed text only; the repository was not edited. Items (1)–(2) are refereed. Item (3) is the judge's and must
carry its label if adopted.

**(1) New paragraph after "Split calibration": what the anisotropy does and does not do.**

> **Anisotropy** (round 24 `theory.md`, `search.md`, refereed). *Pointwise, it does nothing.* Let $\mathcal L$ be the
> arguments "archimedean Taylor upper bound $\ge$ product-formula lower bound" whose lower bound may use integrality,
> toric forced divisibility at every place of $\mathbb Q(i)$, point parity $v_{\bar\pi}=-v_\pi$ and reality of values,
> with auxiliary objects that are isotypic sections for $\sigma_X=\iota\circ\overline{\phantom{x}}$, polynomials in
> group words of A-points (including $\iota y$), or polynomials in the Hilbert-90 vectors ($x_j=P_j/\bar P_j$). Every
> such inequality is implied by toric Liouville for one section (or its $n$-point form at the generating points), so
> on balanced simple profiles with disjoint primes nothing in $\mathcal L$ beats $15/4$ ($\delta=1/30$), and nothing
> beats $\delta=1/10$ for arbitrary profiles, for every family with $\log H(F)+\sum d_i\log d_i=o(\sum d_iw_i)$. The
> reason is that $\sigma_X$ fixes $e$ and orders at $e$, so Galois averaging is lossless; the isotypic parts of
> $H^0(m(-K))$ have maximal orders $\lfloor10m/3\rfloor$ and $\lfloor10m/3\rfloor-1$ in both forms. Twist-sensitive
> congruence divisibility at inert places contributes $O_F(1)$. *Globally, it is visible only through curves.* By
> the parity lemma a balanced twisted curve at $15/4$ has $B\cdot C\ge30$. Such curves exist over $\mathbb R$: an
> explicit one-parameter real family $M$ (two conjugate branches of contact 4, each ray hit once; Jacobian rank 28,
> so a component). Only finitely many curves of its complex component are defined over $\mathbb Q$ (geometric
> monodromy $S_4$ of the place configuration, Faltings; ineffective; none of height $\le300$). The same finiteness
> holds for the split form, so this is not an anisotropy effect. Whether the class has other components is open.
> Blind numerical solvers all missed $M$, so failed searches in this class are not evidence of emptiness. Consequently
> a proof of Question B for $T_5$ that is not also a proof for the split form (hence of Erdős #887) must use a global
> twist-sensitive input; the only one known is the parity lemma.

**(2) In the "Question B" box, after "bounded-degree covering families are excluded (Thm 4.3(ii))", add:**

> The two-branch twisted Pell ansatz ($B\cdot C=30$, contact 4, round 24) yields real curves at exactly $15/4$ but only
> finitely many over $\mathbb Q$ on the known component; it can only enlarge $Z$.

**(3) Optional, judge remark (unrefereed). Append to (2):**

> *(Judge, unrefereed.)* By Thm 4.3(ii), for each fixed $d$ and $\delta\in(0,1/30)$ the curves with $B\cdot C\le d$
> carrying infinitely many A-points lie in a proper closed set; so no fixed curve class at $15/4$ can refute
> Question B, and a refutation by curves needs $\mathbb Q$-curves of unbounded degree.

**Also.** In the "Approximation reading" bullets, after "What is needed is …", add "(round 24: no pointwise use of the
anisotropy helps; see **Anisotropy**)". Round 22 `toric.md` §7's single-branch "no nondegenerate solution" should be
relabelled "blind-search negative; not evidence of emptiness".

## 8. Next step, probabilities, stop or continue

**Decision: stop internal rounds on growth and uniformity** (confirming rounds 19, 20, 22, 23). Round 23 left the
anisotropy as the one unexplored internal lever (a `T_5`-only proof "may need the anisotropy"). This round tested it
from both sides. The pointwise side is now a proved barrier. The curve side produced a real object whose `Q`-scarcity
is form-independent, and §3 shows that no fixed curve class could have refuted Question B anyway. No internal lever
toward growth or uniformity is left on record.

**Single most promising next step (external).** Fold §7 (1)–(2) into `expert_note.md` §3 and put Question B to
approximation-of-rational-points experts (McKinnon–Roth, McKinnon–Satriano, Huang) in this sharpened form:

> On the twisted toric weak-Fano 4-fold of the `A_4` root polytope (rational points `T_5(Q)`, anisotropic norm-one
> torus), is `alpha^ess_e(-K) > 15/4` on Zariski-dense balanced sequences? Every pointwise Liouville-type argument,
> anisotropy-aware or not, stops at `15/4`; the split form is Erdős-#887-hard; the known balanced `15/4` curves lie on
> finitely many `Q`-curves (twisted) versus a hypersurface (split). Is there a **global** mechanism (e.g. McKinnon's
> curve-concentration for toric varieties, or a Brauer/descent-type obstruction for the anisotropic torus) that
> separates the forms?

*Cheap internal items, only if wanted for the record (map value, not progress):* (i) a referee pass on the judge
remarks of §3/§5 (one short deduction each); (ii) Q24.1, a structured search for other components of the `B.C = 30`
class with genus `<= 1` place-curves. Neither can move growth or uniformity.

**Honest probabilities.**

| outcome | probability |
| --- | --- |
| unconditional growth `M = o(log R/loglog R)` from a further internal round | **< 0.1%** |
| an unconditional uniform bound from a further internal round | **essentially 0** |
| a proof of `(L_5^circ)` or `A_5(delta < 1/30)` internally | **< 0.5%** (unchanged; the pointwise use of the anisotropy is now closed) |
| a global twist-sensitive mechanism (beyond parity) found internally | **about 1–2%**, and even then most likely conditional (McKinnon/Vojta-type) |
| the literature already contains a theorem implying `(L_5^circ)` | **about 1%** (unchanged) |
| a refutation of Question B (Zariski-dense twisted A-points, or `M_inf(C) -> inf`) | **about 1–2%** (unchanged; would not count as progress) |
| some `Q`-member of `M` exists (Q24.2) | **about 20–30%** (finite, ineffective; none to height 300; irrelevant to Question B) |
| the `B.C = 30` class has another component with infinitely many `Q`-curves (Q24.1) | **about 10–20%** (would force a surface into `Z`; map value only) |
| the judge remarks of §§2–3, 5 survive a referee | **about 90%** (short deductions from refereed statements plus exact checks) |
| an expert recognises Question B as known or tractable | a few percent; still the only path with non-negligible value |

**Stop.** The record's lasting value is unchanged except for the barrier map: the `(2/3 + eps)` theorem; the
conditional uniform theorem in its smallest form `alpha^ess_e(-K) > 15/4` on a 4-fold (with the absolute form); the
refereed barrier map, now including "anisotropy is pointwise-invisible" (round 24); and the calibrations split ⇒
Erdős #887 and twisted ⇒ absolute Cilleruelo–Granville.

## 9. Constraint audit

The owner's four constraints, for every claim this README relies on:
* **No uniform exceptional set from a fixed-`S` Subspace theorem.** None used by either write-up or by me. The only
  exceptional sets are the hypothetical `Z`, the zero sets of auxiliary sections, and the explicit curves of `M`.
  Theorem B uses Faltings (finiteness of rational points on one curve), not Subspace.
* **No Gaussian gcd or denominator cost discarded.** Theory's `(TL)` uses the full product formula over `Q(i)`,
  primitive tuples and the `2^{e_j}` factors; the inert-place divisibility omitted in the author's list is a *gain*
  left out, not a cost discarded (referee m2). Search keeps `gamma_j` and full toric heights.
* **No drifting common phase replaced by zero.** **Violated once, in theory Prop 6.2** (fixing `P_1 = 1` for a
  non-homogeneous `G`), while the author's audit claimed compliance; repaired by Lemma R2, conclusion unchanged. In the
  search, `c'_j = 1` is forced by `u_j(branch) = 1`, a property of the curve, not a phase.
* **No template-dependent threshold treated as uniform in the radius.** **Violated in theory's first draft** (the
  `O_F(1)` thresholds of Prop 4.3(iii), Thm 7.1(2)); corrected to the explicit proviso
  `log H(F) + sum d_i log d_i = o(sum d_i w_i)`, and the witnessing profiles are stated as family-dependent. Theorem B
  is ineffective and asserts no radius bound. My §3 remark uses gauss Thm 4.3(ii) for each fixed `d` and `delta`, and
  claims nothing uniform in `d`.

## 10. Files

All under `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round24/`.
* Write-ups: `theory.md` (`theory_checks/`, `notes.md`); `search.md` (`search_checks/`, `search_notes.md`;
  `search_v1_backup.md` is the superseded draft).
* Referee reports: `referee_theory.md` (`referee_theory_checks/`, `referee_notes.md`); `referee_search.md`
  (`referee_search_checks/`, `referee_search_notes.md`).
* Judge: this `README.md`, plus `judge_checks/`: `j1_thmB.py` → `j1_thmB.out`, `j1b_inf.py` → `j1b_inf.out` (J1).
* The repository, the overview and the expert note were not edited.
