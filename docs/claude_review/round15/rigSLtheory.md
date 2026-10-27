# Rig^SL: theory and scale (round 15)

Status: complete draft, unrefereed. Scripts and outputs: `rigSLtheory_checks/`. Working notes: `rigSLtheory_notes.md`.

## 0. Summary of results

**No growth improvement and no uniform bound.** `M <= (2/3 + eps) log R/loglog R` stands; P3-all\* stays open and
uncapped in regime (1.4). Rig^SL is not decided. What is new:

**(1) Theory of the angle coupling (§2; proved unless marked).**
* *Four-unit laws* (Cor. 2.2): in SL, opposite units share one quadratic budget, `m_s^2 + m_(s+2)^2 <= D_s^2 + D_(s+2)^2`
  (`= N_v` or `2N_v`), and `prod_s m_s <= (N_v^2/4)|sin 4<v, phi>|` (PU allows `2 N_v^2`). Lemma J0 is put in closed
  form (OPP / ADJ / MAX).
* *Theorem A (one-modulus dichotomy).* Every PU-affordable extra collision is of cofactor type C (`q_0 <= U_(s')`),
  A0 (archimedean, opposite unit saturated) or A1 (archimedean, opposite unit not saturated). C collisions fit at
  `phi = arg pi_j`; **A0 collisions fit at no `phi`**; for A1 the archimedean term is cut by `1 - U_(s'+2)^(-2)`.
  *Exact criterion for the 83/163 phenomenon:* at `first4`, `q_0 = 83` (`K = 4`) and `163` (`K = 5`) every PU survivor
  outside the class carries an A0 collision and none carries an A1 collision; this is exactly why PU is flexible and
  SL rigid there (exact, `classify.py`).
* *Theorem B (pinning: the angle is determined by the labels).* A monomial whose two odd discriminants
  (`Re(Z_v^2) = D_1 D_3`) are fully charged pins `sin 2<v, phi>`; additive triples of such monomials lock `phi` to
  `phi_0 = (arg pi_j)` modulo `pi Z^r` (the arithmetic input is the Q-independence of `pi` and the prime angles,
  Lemma B0). In the one-modulus problem this is **label-independent**, so whenever the triples generate `Z^r` the
  one-modulus SL question is the integer *cofactor criterion*, with no angles left. Exact: the triples generate
  `Z^4` at every forced `q_0 >= 7` for `first4`, `rand4`, `bal4`, `last4` at `K = 4` and `K = 5`.
  Consequences (exact, integer arithmetic only): `rand4`/31 and `bal4`/41 are one-modulus SL-rigid (round 14:
  solver only; here also re-proved by a rigorous interval B&B); one-modulus SL segments `rand4` `[3, 47]` (first
  flexible 59), `bal4` `[3, 41]` (43), `last4` `[3, 29]` (31), `first4` `[3, 103]` (`K = 4`) and `[3, 163]` (`K = 5`,
  first flexible 167), with lifted all-level witnesses; hence `Q_core^SL(rand4, 4) <= 47`, `Q_core^SL(bal4, 4) <= 41`,
  and `109 <= Q_core^SL(first4, 5) <= 163`.
* *Corollary B2 (joint).* An SL system that keeps the genuine odd collisions on a triple family generating `Z^r` has
  `phi = phi_0` and lies in the genuine-angle class GA (`m_(v,s) <= |D_s(Z_v)|`: no archimedean slack at all); to move
  its angle a system must drop genuine odd collisions on enough monomials to break every such family.
* *Theorem C (SL transport).* `(E_s g_s/delta_s)^2 + (E_(s+2) g_(s+2)/delta_(s+2))^2 <= D_s^2 + D_(s+2)^2`: an extra beyond
  the cofactor slack on a monomial with saturated opposite unit forces a dropped collision at a second modulus
  dividing `Re(Z_v^2)` or `Im(Z_v^2)`. Deviations do not come alone; closing this into joint rigidity is open.
* *Theorem P (pinning closure for P3-all\* itself).* A P3-all\* member that is `(eps, t)`-strict on a set `S` of moduli
  with the triple condition TC(`S`) **is an actual cluster**. So **P3-G collapses** in the regime on prime sets with
  TC (round 14's c1construct claim, which failed for norm-1 pinning, holds with all-norm pinning), and a P3-all\*
  fake must break strictness on every TC family. TC holds (exact) for the first `r <= 64` split primes and the loose
  sets with all moduli charged, and with `S = [3, 4r]` (`Y = 4r = Theta(M)`) for `r = 16, 32, 64`.

**(2) Scale (§3; proved).** *Theorem S:* every actual cluster and every member of P3-char / P3-all\* (any `M >= 4`, no
near-extremality, `X >= 2M`) has `r > W/(4W/M + 8 log^+ C)` columns, i.e. `r > M/4` for `C <= 1` and
`r >= (M/4)(1 - O_C(1/log M))` in general. (For actual clusters this is `stability.md` Thm 3.3(3); the extension to
members checks the four inputs of its Prop. 7.1 without the large-prime hypothesis.) So round 14's open
"`r >= M^gamma`?" has the answer `gamma = 1 - o(1)`; every `q <= 2M` is forced already for norm 2 (Cor. S1), and
bounded `K` (`K = 3`) reaches the window `Q_eta` in size. The window lies in the **dense regime** `n_q <= r`, opposite
to all finite tests so far (`n_q >> r`). What bounded `K` still needs is covering (an equidistribution input shared by
fakes), F1 exclusion, and the joint, coherent, strict upgrade.

**(3) Transfer theorem and outcomes (§4).** The corrected Theorem T (RigS^(eps,t) + H^gen + L1 ⇒
`W >= ((8(1-eta) - eps')/K_u) M log^2 M/loglog M` for members) is restated; its size hypothesis is now Theorem S. New:
coherence mod `n` does not give strictness: one-modulus SL class survivors include local associate twists and unit
shifts (`first4`, `K = 4`: at 79, 89, 103; `K = 5`: none up to 151, some at 157, 163; exact). Outcome table (§4.2):
strict-coherent Rig^SL on `Q_eta` gives genuine pair statistics for P3-all\* (growth then rests on H^gen + L1, which
concern actual Gaussian integers); strict coherence on `[3, Theta(M)]` with TC makes P3-all\* exact (Theorem P), which
removes the barrier but is not a growth proof; a stall would give the first genuine C1 seeds but no cap by itself.
**For actual circles no outcome of Rig^SL changes the growth rate by itself.**

## 1. Setting, conventions, what is imported

Labels: *proved* (all quantifiers, proof here or cited exactly), *exact* (finite computation in exact integer /
rational arithmetic, infeasibility by complete search), *interval* (rigorous interval arithmetic, `mpmath.iv` at 40
digits with outward rounding, complete branch-and-bound; a box is discarded only when a violated constraint is
proved on the whole box), *evidence*, *heuristic*, *conditional*. Nothing here is refereed. **Nothing here bounds
actual circles beyond `M <= (2/3 + eps) log R/loglog R`.**

**Conventions** (`round10/L4K.md` §1, `round14/rigtheory.md` §1). Columns `j = 1..r`, distinct `p_j = 1 mod 4`,
fixed Gaussian primes `pi_j`, `phi_0 = (arg pi_j)_j`. For `v in Z^r`: `Z_v = prod pi_j^(v_j^+) conj(pi_j)^(v_j^-) =
X + iY`, `N_v = Norm Z_v = e^(w(v))`, `x_0(v) = arg Z_v = <v, phi_0>`, and `(D_0, D_1, D_2, D_3) = (Y, X - Y, X, X + Y)`.
`nu_s = 1` (`s` even), `2^(-1/2)` (`s` odd). The identity behind everything below (`sharp.md` Prop. 1.2):

```text
nu_s |D_s(Z_v)| = sqrt(N_v) |sin(x_0(v) - s pi/4)|,   so   nu_s^2 (D_s^2 + D_(s+2)^2) = N_v   (s = 0, 1).     (1.1)
```

Labels `ell(q^a) in (Z/m_(q^a))^r`, parity-respecting (`ell_j(q) = [(p_j/q) = -1] mod 2`); `n_q = m_q/4`; `v`
collides at `q^a` with unit `s` iff `<v, ell(q^a)> = s m_(q^a)/4`. Genuine labels `L`: `v` collides at `q^a` with
unit `s` iff `q^a | D_s(Z_v)`. `V_K = {1 <= ||v||_1 <= K}`, `S_h = {||v||_1 <= h}`.

**The three classes** (round-14 README §1).
* **PU**: parity-respecting labels with `m_(v,s)^2 <= nu_s^(-2) N_v` for every `v in V_K`, every `s`.
* **SL = Rig^SL(K)**: PU together with **one** `phi in R^r` such that (1.2) of `sharp.md`,
  `|sin(<v,phi> - s pi/4)| >= nu_s m_(v,s) N_v^(-1/2)`, holds for every `v in V_K` and every `s`. This is the
  norm-`<= K` truncation of P3-all\*.
* **ANG**: PU plus, for each `v in V_K` separately, some `x in R` satisfying the four inequalities (Lemma J0 of the
  round-14 judge). `SL ⊆ ANG ⊆ PU`.

**Symmetries** (referee_rigtheory Prop. T): column conjugations `eps in {+-1}^r`, associate twists `t in {0,1}^r`
(`ell -> ell + (m/2) t`, `phi -> phi + (pi/2) t`), and scalings `u = 1 mod 4`. Coherence is read modulo this group.

**Regime (1.4).** `W <= c_0 M log M`, `p_j <= M^(h_0)`, `r <= rho M`, moduli up to `X = 3 M^A`, `A` chosen after `c_0`;
`Q_eta = {q prime : M^(1-eta) < q <= M, q not dividing 2N}`.

**Imported and used** (all refereed in round 14 or earlier): Lemma J0 (judge, round 14 §5); rigtheory Lemma 3.1,
Theorem 3.4, Cor. 3.5, Theorem 4.1 with referee correction `(eps, t)`; Lemma R (rigcompute); `sharp.md` Def. 1.1,
Prop. 1.2, Cor. 2.2; `lower.md` Def. 1.1, Prop. 1.3; `round8/stability.md` Prop. 2.2, Thm 3.3, Prop. 3.4, Prop. 7.1.

## 2. Rigidity theory for SL (angle coupling)

### 2.1 Lemma J0 in closed form; the four-unit laws

For `v != 0` and demands `m_0..m_3 >= 1` put `a_s = nu_s^2 m_s^2 / N_v` (so PU is `a_s <= 1` for each `s`).

**Lemma 2.1 (J0, closed form; proved).** The four inequalities (1.2) for `v` hold at `x = <v, phi>` iff
`(cos 2x, sin 2x) in B := [2a_2 - 1, 1 - 2a_0] x [2a_3 - 1, 1 - 2a_1]`. Some `x` works iff all three of

```text
(OPP)  a_0 + a_2 <= 1   and   a_1 + a_3 <= 1;
(ADJ)  max(0, 2a_0 - 1, 2a_2 - 1)^2 + max(0, 2a_1 - 1, 2a_3 - 1)^2 <= 1;
(MAX)  max(|1 - 2a_0|, |1 - 2a_2|)^2 + max(|1 - 2a_1|, |1 - 2a_3|)^2 >= 1
```

hold. (OPP) is non-emptiness of `B`; (ADJ) is `min_B |z| <= 1` (it involves one even and one odd unit, i.e. adjacent
axes); (MAX) is `max_B |z| >= 1`.

*Proof.* The first sentence is Lemma J0 (`sin^2 x = (1 - cos 2x)/2`, `cos^2 x = (1 + cos 2x)/2`,
`sin^2(x - pi/4) = (1 - sin 2x)/2`, `sin^2(x - 3pi/4) = (1 + sin 2x)/2`, squared). `B` is a compact convex set, so it
meets the unit circle iff `B != ∅` and `min_B |z| <= 1 <= max_B |z|` (`|z|` is continuous on the connected set `B`).
For `B = I_C x I_S`: `min |z|^2 = d(0, I_C)^2 + d(0, I_S)^2` and `d(0, [l, h]) = max(0, l, -h)`, which gives (ADJ)
with `l = 2a_2 - 1`, `-h = 2a_0 - 1` (and the odd pair likewise); `max |z|^2 = max(l^2, h^2) + max(l'^2, h'^2)`, which
gives (MAX). QED

**Corollary 2.2 (four-unit laws; proved).** Every SL system satisfies, for every `v in V_K`, with `x_v = <v, phi>`:
1. *(quadratic conservation)* `m_0^2 + m_2^2 <= N_v = D_0^2 + D_2^2` and `m_1^2 + m_3^2 <= 2 N_v = D_1^2 + D_3^2`;
2. *(products)* `m_0 m_2 <= (N_v/2) |sin 2x_v|`, `m_1 m_3 <= N_v |cos 2x_v|`, hence
   `(2 m_0 m_2/N_v)^2 + (m_1 m_3/N_v)^2 <= 1`;
3. *(total demand)* `m_0 m_1 m_2 m_3 <= (N_v^2/4) |sin 4 x_v| <= N_v^2/4`.

The genuine system attains these with `|D_0 D_1 D_2 D_3| = |Im(Z_v^4)|/4 = (N_v^2/4) |sin 4 x_0(v)|`; PU allows
`m_0 m_1 m_2 m_3 <= 2 N_v^2`. So the angle coupling lowers the total per-monomial budget by the factor
`8/|sin 4 x_v|`, where `x_v` is **linear** in `v`.

*Proof.* (1) is (OPP) multiplied by `N_v`, together with (1.1). (2): `m_0 <= sqrt(N) |sin x|`, `m_2 <= sqrt(N) |cos x|`,
`m_1 <= sqrt(2N) |sin(x - pi/4)|`, `m_3 <= sqrt(2N) |cos(x - pi/4)|`; multiply in pairs and use
`2 sin(x - pi/4) cos(x - pi/4) = -cos 2x`. (3): multiply the two products, `|sin 2x cos 2x| = |sin 4x|/2`.
Genuine: `Im(Z^4) = 4XY(X^2 - Y^2) = 4 D_0 D_2 D_1 D_3` and `Im(Z^4) = N^2 sin 4x_0`. QED

(1) is the precise sense in which the units of one monomial **share** their budget: opposite units `s, s+2` have
a common quadratic budget, exactly the one the genuine system uses up when both of its discriminants are
fully charged.

### 2.2 The one-modulus dichotomy: cofactor, A0 and A1 collisions

**Setting.** Fix a set of moduli `F` (in Rig^SL(K): `F = F_K^ex`), `q_0 in F`, `F' = F \ {q_0}`. The genuine labels
sit on `F'` at all levels, an arbitrary parity-respecting label `c'` at `q_0` (level 1; this is the level-1 relaxation
of Lemma R(2), valid for certificates). For `v in V_K`: `g_s(v)` = the `F'`-part of `D_s(Z_v)` (all levels), and the
**cofactor** `U_s(v) = |D_s(Z_v)|/g_s(v) >= 1` (it contains the 2-part, the `p_l`-parts, the `q_0`-part and the parts
at primes outside `F`). Demands: `m_s(v) = g_s(v) q_0^([c' collides v at q_0 with unit s])`. At most one unit `s'`
collides for a given `v`; call `(v, s')` a **`q_0`-collision** and classify it:

| type | condition | meaning |
|---|---|---|
| **C** (cofactor) | `q_0 <= U_(s')(v)` | contains the genuine collisions (`q_0 \| D_(s')` gives `q_0 \| U_(s')`) |
| **A0** | `q_0 > U_(s')(v)` and `U_(s'+2)(v) = 1` | archimedean, opposite unit saturated |
| **A1** | `q_0 > U_(s')(v)` and `U_(s'+2)(v) >= 2` | archimedean, opposite unit unsaturated |

With `theta = x_0(v) - s' pi/4` the PU condition for the collision is `q_0 <= U_(s')/|sin theta|` (rigtheory Cor. 3.2):
types A0/A1 are exactly the PU-affordable collisions that use the archimedean factor `1/|sin theta|`.

**Theorem A (one-modulus SL dichotomy; proved).** In this setting:
1. The system satisfies (1.2) at `phi = phi_0` for all `v in V_K` iff every `q_0`-collision is of type C.
2. If some `q_0`-collision `(v, s')` is of type A0, then **no** `phi in R^r` satisfies (1.2) for `v`: (OPP) fails.
3. For a `q_0`-collision of type A1, (OPP) for `v` is equivalent to
   `q_0^2 <= U_(s')^2 + (D_(s'+2)/g_(s'))^2 (1 - U_(s'+2)^(-2))`, whereas PU is
   `q_0^2 <= U_(s')^2 + (D_(s'+2)/g_(s'))^2`. The archimedean term is cut by the factor `1 - U_(s'+2)^(-2)`.
4. Consequently, if no PU survivor outside the class has an A1 collision, then at `q_0`:
   one-modulus SL-flexible ⇔ ANG-flexible ⇔ some label outside the class has only type-C `q_0`-collisions
   (the **cofactor criterion**); and "PU-flexible but SL-rigid" holds iff every PU survivor outside the class has at
   least one A0 collision.

*Proof.* (1) At `phi_0` condition (1.2) for `(v, s)` reads `m_s <= |D_s(Z_v)|` (by (1.1)). Units without a
`q_0`-collision have `m_s = g_s <= |D_s|`. For the colliding unit, `q_0 g_(s') <= |D_(s')| = U_(s') g_(s')` iff
`q_0 <= U_(s')`.
(2) Since only `s'` collides at `q_0`, `m_(s'+2) = g_(s'+2) = |D_(s'+2)|` (as `U_(s'+2) = 1`) and `m_(s') = q_0 g_(s')`.
By (1.1), `a_(s') + a_(s'+2) = nu^2 (q_0^2 g_(s')^2 + D_(s'+2)^2)/N > nu^2 (U_(s')^2 g_(s')^2 + D_(s'+2)^2)/N =
nu^2 (D_(s')^2 + D_(s'+2)^2)/N = 1`, using `q_0 > U_(s')`. So (OPP) fails, and Lemma 2.1 excludes every `x`.
(3) The same computation with `m_(s'+2) = |D_(s'+2)|/U_(s'+2)`, divided by `nu^2 g_(s')^2/N`; for PU drop the
`(s'+2)` term.
(4) By (1) and (2): a label all of whose collisions are C is SL at `phi_0`; a label with an A0 collision is not even
ANG. Without A1 collisions there is no third case. QED

**Remarks.** (a) The theorem is the exact form of the slogan "the archimedean slack of a unit is the budget of the
opposite unit": PU lets an extra collision live on `1/|sin theta|`, but (1.1) says that this factor is
`|D_(s'+2)|/|D_(s')|`-sized room that the opposite unit's genuine demand already occupies when `U_(s'+2) = 1`.
(b) In the **all-norm** one-modulus test of P3-all\* (genuine labels at every `q <= X` not dividing `2N`,
`X >= max |D_s|`), `U_s(v)` is the `2N`-part of `D_s(Z_v)` together with the `q_0`-part. Since `D_1, D_3` are odd,
for odd `s'` the collision is A0 exactly when `D_(s'+2)` is coprime to `N`; for even `s'` when the opposite even
discriminant is odd and coprime to `N`.
(c) Theorem A concerns one modulus with the other labels genuine. The joint version is Theorem C below.

### 2.3 The exact explanation of the 83 / 163 phenomenon (exact)

`classify.py` enumerates every parity-respecting level-1 label at `q_0` (complete column-wise DFS with exact PU
pruning), keeps the PU survivors, and classifies each of their `q_0`-collisions (C / A0 / A1, and for J0 failures
which of OPP, ADJ, MAX fails), with exact integers and rationals. Counts (all / outside the genuine class / outside
every `eps`-class) reproduce the round-14 J2 table and referee_rigcompute §4.6 to the last digit.

| set, `K`, `q_0` (`n`) | PU survivors | outside every `eps`-class: collision profile (number of labels) | ANG out / SL at `phi_0` out |
|---|---|---|---|
| `first4`, 4, **83** (21) | 36 / 24 / 24 | `{C: 8, A0: 4}` (12), `{C: 2, A0: 2}` (12) | 0 / 0 |
| `first4`, 4, 107 (27) | 216 / 144 / 144 | `{C: 2}` (72), `{C: 4}` (36), `{C: 4, A0: 2}` (36) | 108 / 108 |
| `first4`, 5, **163** (41) | 160 / 80 / 80 | `{C: 12, A0: 2}` (40), `{C: 8, A0: 2}` (40) | 0 / 0 |
| `rand4`, 4, 31 (8) | 32 / 24 / 24 | `{C: 22, A0: 2}` (8), `{C: 24, A1: 2, A0: 2}` (8), `{C: 36, A1: 4}` (8) | 8 / 0 |
| `bal4`, 4, 41 (10) | 16 / 8 / 8 | `{C: 14, A1: 2}` (8) | 8 / 0 |
| `last4`, 4, 31 (8) | 32 / 24 / 24 | `{C: 40}` (16), `{C: 42}` (8) | 24 / 24 |

(`first4 = {5,13,17,29}`, `rand4 = {53,61,89,173}`, `bal4 = {97,101,109,113}`, `last4 = {353,373,389,397}`;
forced sets from the round-14 referees' exact files. Every J0 failure in the table is an (OPP) failure with
`U_(s'+2) = 1`, i.e. an A0 collision; there are no ADJ or MAX failures and no A0 collision passes J0, as Theorem A(2)
requires.)

**Reading (exact, for these instances).**
* **83 (`K = 4`) and 163 (`K = 5`) are PU-flexible and SL-rigid for one reason only**: every PU survivor outside
  the class carries 2 or 4 A0 collisions, and none carries an A1 collision. By Theorem A(4) this is exactly the
  phenomenon, and the cofactor criterion decides the modulus.
* At 107 the 108 SL survivors are exactly the labels with only C collisions; the 36 others die by A0.
* `last4` at 31: all collisions are C; that is why the round-14 witnesses there are SL at `phi_0`.
* `rand4` at 31 and `bal4` at 41 are the only instances with A1 collisions. Their ANG survivors (8 each, one
  `U_1`-orbit) need `phi != phi_0`. **The rigorous interval branch-and-bound (`slab_iv.py`) proves that none of the
  8 + 8 labels admits any `phi in [0, pi]^4`** (21 resp. 49 boxes, every box discarded by a proved violation). So
  `rand4` at 31 and `bal4` at 41 are one-modulus SL-rigid at `K = 4` (*interval*, level-1 relaxation at `q_0`, all
  levels elsewhere); round 14 had this as *solver* only. Sanity: with the genuine label, a box of width `2*10^-12`
  around `phi_0` is never discarded (`sanity_iv.py`).

*Why the A1 labels die.* A rigorous deletion filter (`core_iv.py`, same interval B&B) finds tiny infeasible cores:
`rand4`/31: `{e_2, e_4, (0,2,0,1)}`; `bal4`/41: `{e_1, e_3, e_1 + e_4, (1,0,-1,2)}`. In both, the last vector carries
the A1 collision and is an integer combination of the others, and each of the others has **both odd units
saturated** (`U_1 = U_3 = 1`). That is the mechanism of §2.4: such monomials pin their angle, the pins fix `phi`, and
at the pinned `phi` the A1 collision does not fit. §2.4 turns this into an exact, angle-free theorem.

### 2.4 The angle is determined by the labels (pinning)

**Lemma B0 (independence of prime angles; proved).** If `c in Z^r`, `c_0 in Q` and `sum_j c_j arg pi_j = c_0 pi`, then
`c = 0` and `c_0 = 0`.

*Proof.* Write `c_0 = a/b`. Then `e^(2ib sum c_j arg pi_j) = 1`, i.e. `prod_j (pi_j/conj pi_j)^(b c_j) = 1`, i.e.
`prod_j pi_j^(b c_j) = prod_j conj(pi_j)^(b c_j)` up to moving negative powers across. The `pi_j, conj pi_j` are pairwise
non-associate Gaussian primes (the `p_j` are distinct and split), so unique factorisation gives `b c_j = 0`. Then
`c_0 pi = 0`. QED

**Lemma B1 (opposite-pair pinning; proved).** Let `phi` satisfy (1.2) for `v` with demands `m_s`, and suppose an
opposite pair is at full quadratic budget: `nu_s^2 (m_s^2 + m_(s+2)^2) = N_v` (`s in {0, 1}`). Then
`sin^2(x_v - s pi/4) = nu_s^2 m_s^2/N_v` with `x_v = <v, phi>`. In particular, if `(m_1, m_3) = (|D_1(Z_v)|, |D_3(Z_v)|)`
then `sin 2x_v = sin 2x_0(v)`, and if `(m_0, m_2) = (|D_0|, |D_2|)` then `cos 2x_v = cos 2x_0(v)`.

*Proof.* (1.2) gives `nu_s^2 m_s^2 <= N sin^2(x - s pi/4)` and `nu_s^2 m_(s+2)^2 <= N cos^2(x - s pi/4)`; the sum of the
left sides is `N`, so both are equalities. With genuine demands, (1.1) gives `sin^2(x - s pi/4) = sin^2(x_0 - s pi/4)`,
i.e. `cos(2x - s pi/2) = cos(2x_0 - s pi/2)`. QED

In the one-modulus setting of §2.2 put

```text
P0(q_0) = { v in V_K : U_1(v) = U_3(v) = 1 }      (both odd discriminants D_1 D_3 = Re(Z_v^2) fully charged by F').
```

`P0` is closed under `v -> -v` (`D_1(conj Z) = D_3(Z)`). An **additive triple** is `{v, w, v + w} ⊆ P0` with
`v + w != 0` (`v = w` allowed); `T(q_0)` is the union of all additive triples and `Lambda_T` the group it generates.

**Theorem B (pinning; proved).** In the one-modulus setting at `q_0`:
1. *(label-independent)* For **every** label `c'` at `q_0` and every `phi` satisfying (1.2) for all `v in V_K`:
   `sin 2<v, phi> = sin 2x_0(v)` for every `v in P0(q_0)`.
2. *(the angle is determined)* Every such `phi` lies in `phi_0 + pi Lambda_T^*`, where
   `Lambda_T^* = {mu in R^r : <y, mu> in Z for all y in Lambda_T}`. In particular, if `Lambda_T = Z^r` then
   `phi = phi_0 mod pi Z^r`.
3. *(angle-free criterion)* If `Lambda_T = Z^r`, then a label at `q_0` is SL-feasible iff all its `q_0`-collisions are of
   type C; the one-modulus SL question at `q_0` is the integer **cofactor criterion**, with no angles left, and
   SL = SL at `phi_0`.

*Proof.* (1) Let `v in P0`. If `c'` collides `v` at an odd unit `s'`, the collision is A0 (`U_(s'+2) = 1`, and
`q_0 > 1 = U_(s')`), and Theorem A(2) excludes every `phi`. Otherwise the odd demands are `g_1 = |D_1|`,
`g_3 = |D_3|`, and Lemma B1 applies.
(2) Let `phi` be a solution and `eta = 2(phi - phi_0)`. Since `sin a = sin b` iff `a = b` or `a = pi - b` mod `2pi`,
(1) gives `tau_y in {0, 1}` with `<y, eta> = tau_y (pi - 4 x_0(y)) mod 2pi` for `y in P0`. For an additive triple
`(v, w, u = v + w)`, linearity of `<., eta>` and of `x_0(y) = <y, phi_0>` gives
`sum_j c_j arg pi_j + (tau_u - tau_v - tau_w) pi in 2 pi Z` with `c = 4(tau_v v + tau_w w - tau_u u) in Z^r`. By Lemma
B0, `c = 0` and `tau_u - tau_v - tau_w` is even. `c = 0` reads `tau_v v + tau_w w = tau_u (v + w)`: for independent
`v, w` this gives `tau_v = tau_w = tau_u`; for `w = lambda v` (`lambda != 0, -1`) a check of the eight 0/1 patterns
gives the same. Then `tau_u - tau_v - tau_w = -tau` is even, so `tau = 0`. Hence `<y, eta> in 2 pi Z` for all `y in T`,
hence for all `y in Lambda_T`, i.e. `phi - phi_0 in pi Lambda_T^*`.
(3) By (2), `phi = phi_0 + pi k`, `k in Z^r`, which changes every `<v, phi>` by a multiple of `pi` and leaves (1.2)
unchanged. So (1.2) holds at `phi_0`, and Theorem A(1) applies. QED

**Exact computation** (`triple_gen.py`: exact integer test of `Lambda_T = Z^r` by Smith normal form). For `first4`,
`rand4`, `bal4`, `last4` at `K = 4`, at **every** forced `q_0 >= 7`: `[Z^4 : Lambda_T] = 1` (`|P0|` between 18 and 81).
The only exceptions are `q_0 = 3` (index 2, all four sets) and `q_0 = 5` (`rand4`, `bal4`: index 2; `last4`: rank `< 4`),
where the small prime is itself uncharged; at those two moduli PU already has no survivor outside the class, so
nothing is lost. At `K = 5` the same holds: index 1 at every forced `q_0 >= 7` for all four sets (`tg_<set>_K5.txt`). Consequently, **for these four prime sets at `K = 4`
the one-modulus Rig^SL problem at every forced `q_0 >= 7` is exactly the integer cofactor criterion**: the SL counts
of `classify.py` ("SL at `phi_0`") are the exact SL counts, and the interval certificates of §2.3 for `rand4`/31 and
`bal4`/41 are re-proved in pure integer arithmetic. (`phi_pin.py`, an independent interval-arithmetic enumeration of
all `2^r d` candidate angle vectors from a basis of `P0`, confirms `Phi(q_0) = {phi_0}` at `rand4`/31, `bal4`/41,
`last4`/31, `first4`/83, 107 (`K = 4`) and `first4`/163 (`K = 5`).)

**Exact one-modulus SL segments (angle-free, by Theorem B(3)).** One-modulus test (genuine labels on
`F_K^ex \ {q_0}` at all levels, level-1 relaxation at `q_0`), counts outside every `eps`-class, `scan_onemod.sh`:

| set, `K` | one-modulus SL-rigid at every forced `q_0` in | first SL-flexible `q_0` (SL labels outside every `eps`-class) | PU first flexible (round 14) |
|---|---|---|---|
| `first4`, 4 | `[3, 103]` | 107 (108) | 83 |
| `first4`, 5 | `[3, 163]` | 167 (24); also 173 (42) | 163 |
| `rand4`, 4 | `[3, 47]` | 59 (752) | 31 |
| `bal4`, 4 | `[3, 41]` | 43 (40) | 41 |
| `last4`, 4 | `[3, 29]` | 31 (24) | 31 |

At each first SL-flexible modulus 20 labels were lifted to level 2 without level-2 collisions and re-verified with
the all-level demand: PU-good, (1.2) exact at `phi_0`, outside every `eps`-class (hence every `(eps, t)`-class)
(`lift_check.py`, outputs `lift_out.txt`, `lift_first4_K5_167.txt`). So `Q_core^SL(rand4, 4) <= 47`,
`Q_core^SL(bal4, 4) <= 41`, `Q_core^SL(last4, 4) <= 29` and `Q_core^SL(first4, 5) <= 163` (*exact*; all but the third
are new). With the transferred PU certificate, **`109 <= Q_core^SL(first4, 5) <= 163`** (round 14: only `>= 109`). These are one-modulus statements: they bound the joint SL core from above; they give no
lower bound on it (the joint PU core of the loose sets is 7, and the joint SL core of the loose sets is open).

**Corollary B2 (joint form; proved).** Let `ell'` be any SL system on `F` (all moduli free) with angle `phi`. Let
`P0^F = {v in V_K : U^F_1(v) = U^F_3(v) = 1}` (cofactors relative to all of `F`) and let `R(ell') ⊆ P0^F` be the set of
`v` for which `ell'` **retains** every genuine odd collision (`a'_q(v, s) >= a^gen_q(v, s)` for `s` odd, all
`q in F`). If the additive triples inside `R(ell')` generate `Z^r`, then `phi = phi_0 mod pi Z^r`, and `ell'` lies in
the **genuine-angle class** GA:

```text
m'_s(v) <= |D_s(Z_v)|      for all v in V_K and all units s,
```

i.e. Rig(K) with the integer budget `|D_s(Z_v)|` in place of `nu_s^(-1) sqrt(N_v)`. In GA the slack of the genuine
system at `(v, s)` is exactly `log U^F_s(v)` (Lemma 3.1 of rigtheory without its archimedean term).

*Proof.* For `v in R(ell')`, `m'_1 >= |D_1|` and `m'_3 >= |D_3|`; Cor. 2.2(1) forces equality and Lemma B1 pins
`sin 2x_v`. The rest is the proof of Theorem B(2)-(3) with `phi_0` evaluated through (1.1). QED

**Reading.** Theorem B is the precise sense in which "per-`v` compatibility plus parity" is stronger than PU: a
monomial whose two odd discriminants are fully charged is a **rigid bar** in the angle space (it fixes
`sin 2<v, phi>`), additive triples of bars lock the angle completely (Lemma B0 is the arithmetic input: the prime
angles are independent of `pi`), and once the angle is the genuine one there is **no archimedean slack at all**:
only cofactor slack `U` remains. Corollary B2 is the joint statement: an SL system can move its angle away from
`phi_0` only by **dropping** genuine odd collisions on enough monomials to break every generating family of
triples; each such drop is a non-genuine behaviour at a modulus dividing `Re(Z_v^2)`.

### 2.5 The joint transport inequality in SL

**Theorem C (SL transport; proved).** Let `ell'` be an SL system on `F` (norm `<= K`, all levels, one `phi`), `v in V_K`,
`s in {0, 1}`. For `t in {s, s+2}` and `q in F` let `a'_q(t)`, `a^gen_q(t)` be the colliding levels of `(v, t)`, and put
`delta_t = prod_q q^((a^gen_q(t) - a'_q(t))^+)` (dropped genuine part), `E_t = prod_q q^((a'_q(t) - a^gen_q(t))^+)`
(extra part), `g_t = g^F_t(v)`, `U_t = U^F_t(v)`. Then

```text
(E_s g_s / delta_s)^2 + (E_(s+2) g_(s+2) / delta_(s+2))^2  <=  (U_s g_s)^2 + (U_(s+2) g_(s+2))^2 .        (2.1)
```

In particular, if `U_(s+2) = 1` and `E_s > U_s delta_s` (an extra at unit `s` exceeding the cofactor slack and the
drops of unit `s`), then `delta_(s+2) > 1`: `ell'` drops a genuine `(v, s+2)` collision at some `q_1 in F`, and
quantitatively `g_(s+2)^2 (1 - delta_(s+2)^(-2)) >= g_s^2 (E_s^2/delta_s^2 - U_s^2)`.

*Proof.* `m'_t = R_t E_t` with `R_t = g_t/delta_t` (retained part), so (2.1) is Cor. 2.2(1) with
`|D_t| = U_t g_t`. The special case: `E_(s+2) >= 1`, and if `delta_(s+2) = 1` the left side exceeds
`U_s^2 g_s^2 + g_(s+2)^2`. QED

*Comparison with PU* (rigtheory Cor. 3.5(3)): there an extra of size `E_s` is affordable without any drop as long as
`log E_s <= log delta_s + log U_s + log(1/|sin theta|)`; the archimedean term pays for it. In SL, when the opposite
unit is saturated, the archimedean term is not available: the extra must be paid by the cofactor or by a drop of the
opposite unit **elsewhere**, i.e. by a second non-genuine modulus `q_1 | D_(s+2)(Z_v)`. Combined with rigtheory Theorem
3.4 (at a covered `q_0` every label outside the class has a non-genuine collision `v in (D - D) \ Lambda_c`), (2.1)
gives the **propagation statement**: at a covered `q_0`, a non-genuine label whose forced extra `(v, s')` has
`U_(s'+2)(v) = 1` and `q_0 > U_(s')(v) delta_(s')` forces a deviation at some `q_1 in F \ {q_0}` dividing
`D_(s'+2)(Z_v)` (or dividing `D_(s')(Z_v)`, if `delta_(s') > 1`), so dividing `Re(Z_v^2)` or `Im(Z_v^2)`.
Deviations in SL do not come alone. This is the "transport" half of the mechanism; turning it into joint rigidity
needs a counting argument showing that the chain of forced deviations cannot close up within the slack, which is not
proved (§6).

### 2.6 Pinning closure: strict labels on few moduli turn a P3-all\* member into an actual cluster

Theorem B has an all-norm consequence for P3-all\* itself (not only for its truncations). For `(eps, t)` let
`pi'_j = i^(t_j) pi_j^(eps_j)` (`pi^(-1)` meaning `conj pi`), `Z'_v`, `x'_0(v) = arg Z'_v`, `phi'_0 = (arg pi'_j)`.
A member is **`(eps, t)`-strict at `q`** if `ell(q^a) = u_(q^a) (eps∘L(q^a) + (m_(q^a)/2) t)` with `u = 1 mod 4` at every
level `q^a <= X` (round-14 corrected notion). Then `v` collides at `q^a` with unit `s` iff `q^a | D_s(Z'_v)`.

**Definition (triple condition).** For a set `S` of odd primes not dividing `N`, let

```text
P0'(S) = { v != 0 : gcd(Re(Z'_v^2), N) = 1, every prime power q^a || Re(Z'_v^2) has q in S and q^a <= X }
```

(`Re(Z^2) = X^2 - Y^2 = D_1 D_3` is odd). **TC(S)** holds if the additive triples `{v, w, v + w} ⊆ P0'(S)` generate
`Z^r`.

**Theorem P (pinning closure; proved).** Let a member of `F^all_X(M, C)` (P3-all\*) be `(eps, t)`-strict at every
`q in S`, and suppose TC(S). Then its angle vector satisfies `phi = phi'_0 mod pi Z^r`, and its points
`R e^(i(theta_0 + delta_x))` are the Gaussian integers `+- i^(-k_x) prod_j pi'_j^(a_xj) conj(pi'_j)^(e_j - a_xj)`: the
member **is an actual configuration** of `M` lattice points of norm `N` on an arc of length `<= C sqrt R`.

*Proof.* Let `v in P0'(S)`. By strictness, `(v, 1)` collides at every `q^a || D_1(Z'_v)` (all such `q in S`,
`q^a <= X`, `q` not dividing `N`), so `m_(v,1) >= |D_1(Z'_v)|`, and likewise `m_(v,3) >= |D_3(Z'_v)|`. (e-all) holds for
every `v`, so Cor. 2.2(1) gives equality and Lemma B1 gives `sin 2<v, phi> = sin 2x'_0(v)`. The proof of Theorem B(2)
(with Lemma B0 for the primes `pi'_j`, which are again pairwise non-associate) gives `phi - phi'_0 in pi Z^r` under
TC(S). By Def. 1.1(b), `e^(i(theta_0 + delta_x)) = i^(-k_x) e^(i<2a_x - e, phi>)` and
`e^(i(2a_xj - e_j) phi'_j) = (pi'_j/sqrt p_j)^(a_xj) (conj pi'_j/sqrt p_j)^(e_j - a_xj)`; the shift by `pi Z^r` contributes a
sign. With `R = sqrt N` the point is `+- i^(-k_x) w_x`, `w_x = prod pi'_j^(a_xj) conj(pi'_j)^(e_j - a_xj) in Z[i]`,
`|w_x|^2 = N`. Distinct rows give non-associate `w_x` (unique factorisation), and the angular positions
`theta_0 + delta_x` lie in an interval of length `Delta = C e^(-W/4)`, i.e. on an arc of length `C sqrt R`. QED

**What this says.** c1construct (round 14) claimed "a P3-G member is an actual cluster" from **norm-1** pinning, under a
per-column condition (i) that fails for 80-91% of columns of small split primes (referee). Theorem P replaces norm-1
pinning by pinning through **all** odd-saturated monomials, locked by additive triples, and the hypothesis becomes the
lattice condition TC, which is generic:

*Exact* (`triple_allnorm.py`, Hermite reduction, sanity-tested on lattices of index 2, 3, 4): with `S` = all primes
`<= X` not dividing `2N` (so `P0'` is "`gcd(Re Z_v^2, N) = 1`"), TC holds using triples inside `V_3` for the first
`r = 4, 8, 16, 32, 64` split primes and for `rand4`, `bal4`, `last4`; with `S` = the primes `<= Y` only (`Re(Z_v^2)`
required `Y`-smooth), TC holds using triples of vectors of `V_2` (sums in `V_4`) for the first `r = 16, 32, 64` split
primes at `Y = 4r` (and at `Y = 16r` for `r = 16, 32`), and for `r = 64` already at `Y = 2r`; it fails (rank `< r`, or index 2 for `r = 32`,
`Y = 2r`) at `Y <= 2r` otherwise. Since `M/4 <~ r <= rho M`, `Y = 4r` is `Theta(M)`.

**Consequences (proved).**
1. **P3-G collapses.** Every member of P3-all\* whose labels are `(eps, t)`-strict at all `q <= X` (the class P3-G of
   round 7, up to the symmetry) and whose primes satisfy TC is an actual cluster. In regime (1.4) (`X = 3M^A`, `A`
   large) P3-G contains **no fake** on prime sets with TC. (The only P3-G fakes known have `log p >= 2.08 X`
   (rigtheory §4.2), far outside the regime; this is consistent.)
2. **What a P3-all\* fake must do.** A member that is not an actual cluster must, for every generating family of
   additive triples, fail strictness at some prime factor of some `Re(Z'_v^2)` in the family, for **every** symmetry
   `(eps, t)`. (TC is symmetric: `P0'` for `(eps, t)` is `eps∘P0'` for `(1, 0)`, because `Z'_v = +-i^k Z_(eps∘v)` and
   `Re` changes at most its sign; so the computed instances give TC for every `(eps, t)`.) With `S = [3, Y]` and TC(`[3, Y]`): **a fake cannot be `(eps, t)`-strict and coherent on all of
   `[3, Y]`**. So the rigidity statement that would close P3-all\* (make it exact) is strict coherence on an initial
   segment `[3, Y]` with `Y` of the order of `M` (evidence: `Y = 4r`, i.e. `Y ~ M` in the regime by Theorem S), not
   on the window `Q_eta` alone.
3. Theorem P does not bound actual circles: a member that collapses is an actual cluster, about which it says
   nothing. Its role is in the barrier map: it identifies exactly which label freedom a P3-all\* fake must use.

## 3. Scale: number of split primes in regime (1.4)

The round-14 README (§2 item 3, §8) left open whether near-extremal clusters in regime (1.4) must have `r >= M^gamma`,
and stated that "(1.4) only gives `r <= rho M`; near-extremality forces `sum e_j >~ 3M/h_0`, not many columns". The
answer was already essentially in `round8/stability.md` Theorem 3.3(3), whose proof **uses no near-extremality
hypothesis**. What is needed in addition is that the bound holds for **members** of P3-all\* (the transfer theorem
applies rigidity to members, not only to actual clusters); that is Theorem S(b).

**Theorem S (every cluster and every member has more than `M/4` columns; proved).** Let `M >= 4`, `C > 0`.
* (a) Every normalised actual configuration of `M` lattice points on an arc of length `<= C sqrt R` of
  `x^2 + y^2 = R^2`, and
* (b) every member of `lower.md`'s class `A_X(M, C)` with `X >= 2M`, hence every member of P3-char `F^char_X(M, C)` and
  of P3-all\* `F^all_X(M, C)` (`sharp.md` §1.1: `F^all_X ⊆ F^char_X ⊆ A_X`),

has `r` columns (distinct primes `p_j | N'`) with

```text
r  >  W' / (4W'/M + 8 log^+ C)  =  (M/4) (1 + 2 M log^+ C / W')^(-1).                                  (3.1)
```

In particular `r > M/4` when `C <= 1`, and `r >= (M/2)(1 - M^(-1/3))^2` when `C <= 1`, `M >= 8`. Since
`W' >= 3 M log M - (K_0 + 4 log^+ C) M` (`two_thirds` Thm 7.5, valid on the class by `lower.md` Prop. 1.3),
`r >= (M/4)(1 - O_C(1/log M))` for every fixed `C`. No hypothesis on `R`, on near-extremality or on the size of the
`p_j` is used.

*Proof.* (a) is `stability.md` Theorem 3.3(3) and Proposition 3.4. Their proofs use only `M >= 4`, the chord bound
(1.1) there, and the gcd lemma (Prop. 2.2) for subsets `B` with `|B| >= M/2 - M^(2/3)/2`; the gcd lemma is (T1) of
`two_thirds` applied to the sub-configuration `u^B`. Recall the argument of Theorem 3.3(2)-(3): for `p | N'` let
`t* = min{t : #{x : a_xp <= t} >= M/2}`; the sets `L* = {a_xp <= t*}` (size `>= M/2`) and `U* = {a_xp >= t*}`
(size `> M/2`) have common divisors `conj(pi)^(e - t*)` and `pi^(t*)`, so `e_p log p <= gamma_(L*) + gamma_(U*) <=
2 beta'` with `beta' < 2W'/M + 4 log^+ C`; summing over `p | N'` gives `W' <= 2 beta' omega(N')`, i.e. (3.1).
(b) The same proofs apply verbatim to angle models as soon as (T1) holds for every **sub-model**; this is the content
of `stability.md` Prop. 7.1, proved there for models whose primes all exceed `Q_0 >= 5M/4`. We remove that hypothesis.
For a member and `B ⊆ [M]`, `|B| = s >= 2`, let `alpha_j = min_B a_xj`, `beta_j = max_B a_xj`, and define the sub-model by
the rows `a_x - alpha` (`x in B`), exponents `beta - alpha` (columns with `alpha_j = beta_j` dropped), the same `phi`,
`delta_x`, `k_x`, and `theta_0^B = theta_0 + <e - alpha - beta, phi>`; then `N_B = N e^(-gamma_B)`, the chord constant is
`C_B e^(-gamma_B/4)`, and the reduced chords `|beta_xy| := |z_x - z_y| e^(-(W - d_xy)/2)` are unchanged. The four inputs
of Prop. 7.1:
1. *algebra* (cut identity, (2.1)): identities in the profile, unchanged;
2. *chord bound*: as computed;
3. *pair lower bound* `|beta_xy| >= sqrt2 m_xy`: this is (e) of `A_X` for `c = e_x - e_y` (`w(v(c)) = d_xy`, and
   `|z_x - z_y| = 2R |sin((delta_x - delta_y)/2)|`), where `m_xy` contains the collisions at `q` not dividing `N` and the
   same-level collisions at the `p_j`;
4. *pigeonhole* at every prime power `q^a < 5s/4 <= X` not dividing `N_B`: (i) for `q` not dividing `N`, multiply the
   residues `rho_x(q)`, `x in B`, by one fixed unit of norm `N_B/N mod q^a`: pair collisions are unchanged and the
   residues lie on the conic of norm `N_B`; (ii) for `q = p_j` with column `j` constant on `B`, every pair in `B` is a
   same-level pair, and the level residues `rho_x(p_j)` (norm `N/p_j^(e_j)`, `lower.md` Def. 1.1(d)), multiplied by a
   fixed unit, lie on the conic of norm `N_B`; their pair collisions are exactly the same-level collisions charged in
   (e). At the `p_j | N_B` the same-level clause is the input `two_thirds` uses (Lemma 3.3 (S|) and the level-split
   class counts), exactly as in `lower.md` Prop. 1.3.

So the sub-model satisfies every input of `lower.md` Prop. 1.3 (the gap condition (c) is not used), (T1) holds for it,
and Prop. 2.2, Theorem 3.3 and Prop. 3.4 of `stability.md` hold for members, with `G = 1`, `W' = W` (members are
normalised by Def. 1.1(a)). QED

**Corollary S1 (no free moduli on the window; proved).** Let `r` be the number of columns of a member or cluster,
`K >= 2`, and `q` an odd prime not dividing `N` with `n_q < 2r + 1`. Then `q` is forced for norm 2 (every
parity-respecting label has a level-1 collision of some `v` with `1 <= ||v||_1 <= 2`). Under (3.1) with `C <= 1`
(or `M >= M_0(C)`), every `q <= 2M` not dividing `2N` is forced for norm 2; in particular `Q_eta ⊆ F_K^ex` for
every `K >= 2`.

*Proof.* The `2r + 1` points of `S_1 = {0, +-e_j}` have labels in `Z/n_q`; two coincide, and their difference
`v in S_1 - S_1 ⊆ S_2` collides. For `q <= 2M`: `n_q <= (q + 1)/4 < M/2 + 1 < 2r + 1`. QED

**Corollary S2 (the bounded-`K` size obstruction disappears; proved).** In regime (1.4), every member and every
actual cluster has `M/4 (1 - o(1)) <= r <= rho M`. Hence:
1. By Cor. S1, `F_K^ex ⊇ Q_eta` for `K >= 2`, so the round-14 ceiling (`q_core(K) <= max F_K^ex <= 4A_K + 5`) lies above
   the window; and the anticode range `n_q <= |S_h|`, `|S_h| >= (2r)^h/h!` (`h = floor((K-1)/2)`), contains the window
   already at `K = 3` (`|S_1| = 2r + 1 > M/2 + 1 > n_q` for `q <= M`).
   The hypothesis "`r >= M^gamma`" of the corrected transfer theorem (round-14 README §1.2) **holds with `gamma = 1 - o(1)`**,
   unconditionally.
2. The window lives in the **dense regime** `n_q <= r (1 + o(1))`: the cyclic label group at `q in Q_eta` is no
   larger than the number of columns. The finite Rig/Rig^SL tests of rounds 10-15 (`r = 4..8`, `n_q` up to 260) sit
   in the opposite, **sparse** regime `n_q >> r`, so they probe the shape of the core but not the asymptotic regime.
3. What bounded `K` still needs on `Q_eta` is not size but (i) **covering** (`c(S_h) = Z/n_q`; for `h = 1` this is a
   coupon-collector statement for the `2r + 1` values `{0, +-L_j(q) mod n}`, plausible only for `q <~ 8r/log r`; for
   `h = 2`, `|S_2| ~ 2 r^2 >> n_q log n_q`, it is a mild equidistribution statement for products of two Gaussian
   primes mod `q`; an L7/H-type input, shared by fakes), (ii) F1 exclusion, and (iii) the joint, coherent,
   `(eps, t)`-strict upgrade.

*Proof.* (1) is Theorem S with Theorem 2.5 of rigtheory and Prop. 3.6 there. (2): `n_q <= (M + 1)/4`. (3) is a
reading of rigtheory Theorem 3.4 (heuristic counts labelled as such). QED

**Remark (the large sieve does not prove covering).** Covering of `Z/n_q` by the slopes of the `pi_j` mod `q` for
most `q <= M` would follow from a large-sieve inequality only if the number of points were `>> q^2`-ish per modulus;
here there are `~ M/2` points (or `~ M^2/8` for `S_2`) of height `M^(O(h_0))`, and the two-dimensional large sieve
gives only trivial bounds (heuristic estimate in the notes). Covering stays an open equidistribution input.

## 4. Corrected transfer theorem and what each outcome of Rig^SL means

### 4.1 The corrected transfer theorem (restated)

Fix `c_0, eta, delta, lambda_0, K_u`. Regime (1.4) as in §1. For `(eps, t) in {+-1}^r x {0,1}^r` put
`pi'_j = i^(t_j) pi_j^(eps_j)` and `Z'_v` as in §2.6.

> **RigS^(eps,t)(eta, c_0).** There is `A_1(c_0)` such that for `A >= A_1`, `M >= M_0`: every member of
> `F^all_(3M^A)(M, C)` with `p_j <= M^(h_0)`, `r <= rho M`, `W <= c_0 M log M` is **`(eps, t)`-strict and coherent on
> `Q_eta`**: one `(eps, t)` such that `ell(q^a) = u_(q^a) (eps∘L(q^a) + (m_(q^a)/2) t)`, `u_(q^a) = 1 mod 4`, for every
> `q in Q_eta` and every level.
>
> **H^gen(eta, delta).** For every such member, with its profile `(a_x)`, windings `(k_x)` and its `(eps, t)`, the
> integers `c_xy = #{q in Q_eta : q | D_(k_x - k_y)(Z'_(a_x - a_y))}` satisfy the factorial-moment inequalities of
> `H(eta, delta)` (`p3lower.md` §4).
>
> **L1(K_u, lambda_0).** `K_u`-uniform margins and `lambda_Q >= lambda_0`.

**Theorem T (transfer; conditional; rigtheory Thm 4.1 with the referee's `(eps, t)` correction).** Under RigS^(eps,t),
H^gen and L1, every member as in RigS has `W >= ((8(1 - eta) - eps')/K_u) M log^2 M / loglog M` for `M >= M_0`. Hence
no member has `W <= c_0 M log M` for large `M`; in particular near-extremal actual clusters with `K_u`-uniform margins
satisfy `M <= (K_u/(4(1 - eta)) + o(1)) log R logloglog R/(loglog R)^2`.

*Proof (as in rigtheory, with the fix).* A pair `(x, y)` collides at `q` with unit 1 iff
`<a_x - a_y, ell(q)> = iota (k_x - k_y) mod m_q` (`sharp.md` Cor. 2.2), `iota = log i`. With
`ell = u (eps∘L + (m/2) t)`, `u = 1 mod 4`: `u^(-1) iota = iota mod m`, and `eps∘L + (m/2) t` is the genuine label of
the primes `pi'_j` (`log(i pi/conj(i pi)) = log(-pi/conj pi) = L + m/2`; conjugation negates `L`). So the condition is
`q | D_(k_x - k_y)(Z'_(a_x - a_y))`, the member's pair integers are the `c_xy` of H^gen, and `p3lower.md` Thm C gives
the bound. QED

**What changed this round.**
1. *The size hypothesis is now a theorem.* Round 14 required, for RigS to be reachable at bounded `K`, the extra
   hypothesis `r >= M^gamma`. By Theorem S, `r >= (M/4)(1 - o(1))` for every member and every actual cluster, so
   `gamma = 1 - o(1)`; by Cor. S1, `Q_eta ⊆ F_2^ex` (no free moduli on the window); `K = 3` already reaches the window
   in size.
2. *The class is Rig^SL.* RigS^(eps,t) follows from the statement "every Rig^SL(K) system on the forced set of the
   member's primes is `(eps, t)`-strict and coherent on `Q_eta`" for some fixed `K`, since the norm-`<= K` truncation of
   a member (labels at all levels up to `X`, its own `phi`) is a Rig^SL(K) system. PU statements cannot give it
   (PU witnesses exist where SL is rigid, §2.3).
3. *Strictness is a separate requirement, and SL does not supply it at small `K`* (exact, `strict_check.py`, one
   modulus, `first4`, forced sets, all forced `q_0 in [7, 163]`; `scan_strict.sh`): at `K = 4`, inside the one-modulus
   SL-rigid range, the class survivors at `q_0 = 79, 89, 103` include labels that are strict only for a **local** twist
   `t` at `q_0`, or for no `(eps, t)` at all (unit shifts), and they are SL-feasible at `phi_0` (89: 80 class
   survivors, 40 strict for some `(eps, t)`, 20 strict for the global `(1, 0)`); every other `q_0 <= 101` has only
   `(1, 0)`-strict survivors. At `K = 5` every class survivor at every forced `q_0 <= 151` is `(1, 0)`-strict; local
   twists reappear at 157 and 163 (48 resp. 80 class survivors, half `(1, 0)`-strict). So the strict one-modulus
   segment grows from `[7, 73]` (`K = 4`) to `[7, 151]` (`K = 5`), like the rigid one (exact data; growth is
   evidence). A coherence-only Rig^SL core (mod `n`) does **not** yield RigS^(eps,t); strictness must be certified
   separately.
4. *Theorem P sharpens what RigS would buy.* If strict coherence held not just on `Q_eta` but on an initial segment
   `[3, Y]` with TC(`[3, Y]`), every member would be an actual cluster (§2.6), so P3-all\* would be exact in the
   regime. RigS on `Q_eta` alone gives only the pair statistics (Theorem T).

### 4.2 What each outcome of Rig^SL would mean

| outcome | for the barrier map (P3-all\* in regime (1.4)) | for actual circles' growth rate |
|---|---|---|
| **Rig^SL strict-coherent on `Q_eta`** at some bounded `K` (or `K` growing slowly), for all prime sets of the regime | RigS^(eps,t) holds; P3-all\* members have **genuine** pair statistics on `Q_eta` (Theorem T); a cap of P3-all\* must be a fake whose `Q_eta`-statistics are those of genuine twisted primes. With Theorem S no size hypothesis is left | **No change by itself.** RigS is vacuous for actual clusters (their labels are genuine). Growth `M = o(log R/loglog R)` would follow from H^gen + L1, which are statements about actual Gaussian integers (L7/I3-type), open and unchanged |
| **Rig^SL strict-coherent on `[3, Y]`**, `Y ~ M`, with TC(`[3, Y]`) | every member is an actual cluster (Theorem P): **P3-all\* is exact** in the regime, so it can no longer serve as a fake barrier; the fake-circle barrier for local classes would have to be rebuilt from a strictly smaller class | **No change by itself**: P3-all\* = actual clusters means a growth proof for P3-all\* is a growth proof for actual circles, which is the original problem. It removes a barrier, it does not supply an argument |
| **Rig^SL coherent but not strict** (local twists / unit shifts survive on a positive proportion of `Q_eta`) | RigS^(eps,t) fails; Theorem T does not apply; pair statistics of members are those of a per-modulus **mixed** twist, `c_xy = #{q : q | D_(k_xy + 2<a_xy, t_q>)(Z_(a_xy))}`, with `t_q` chosen per modulus. A growth argument would need H for such mixed-twist pair integers (an adversarial version of H^gen) | none |
| **Rig^SL stalls** (non-coherent SL systems at every bounded `K`, on a positive `theta`-mass of `Q_eta`) | refutes RigS at every bounded `K`. Does **not** cap P3-all\*: a cap still needs completion to all norms, all moduli `<= X`, the pair conditions and, by Theorem P, failure of strictness on every TC family. Stalling SL systems would be the first genuine C1 seeds (labels and angles consistent at norm `<= K`) | none; a cap of P3-all\* would only confirm that local inputs cannot prove growth |

**Bottom line for the growth rate.** No outcome of Rig^SL moves the growth rate of actual circles by itself: the
growth input that remains in every branch is a statement about genuine Gaussian integers (H^gen + L1, or the
original problem itself once P3-all\* is exact). What Rig^SL decides is the barrier map: whether the last uncapped
local class P3-all\* is a genuine barrier or collapses onto actual circles.

## 5. Computations (scripts, exact certificates)

All scripts are in `rigSLtheory_checks/` (own code; `sllib.py` is a copy of the round-14 judge's `jglib.py`, residue-field
discrete logs and `D_s`; the forced sets are copies of the round-14 referees' exact files `forced_K4.txt`,
`forced_K5.txt`, `rforced_<set>_K<K>.json`; nothing is imported from the repository). Run from that directory with
`python3`. Conventions: `L_j(q) = log(pi_j/conj pi_j)` with `log i = m/4`, `pi_j = a + bi` with `0 < a < b`
(`sllib.two_sq`); the round-14 sanity check of these conventions (660,480 label/`D_s` comparisons) applies.

| script | what it decides | type |
|---|---|---|
| `classify.py FFILE q0 K ps` (env `DUMP`, `DUMPSL`, `DUMPPU` write label lists) | complete column-wise DFS over parity-respecting level-1 labels at `q0` with exact PU pruning; classification of every `q0`-collision (C / A0 / A1; OPP / ADJ / MAX for J0 failures); counts PU / ANG / SL-at-`phi_0`, inside / outside the genuine and `eps`-classes. Reproduces round-14 J2 and referee §4.6 counts exactly | exact |
| `slab_iv.py`, `run_angA1.py`, `sanity_iv.py` | rigorous interval B&B for one `phi` (mpmath `iv`, 40 digits, outward rounding; a box is discarded only on a proved violation; a surviving tiny box returns UNDECIDED, never INFEASIBLE). `rand4`/31 and `bal4`/41: all 8 + 8 A1 labels INFEASIBLE (21, 49 nodes) | interval |
| `core_iv.py` | deletion filter for minimal infeasible cores (output `core_rand4_K4_31.txt`, `core_bal4_K4_41.txt`) | interval |
| `phi_pin.py` | enumerates all `2^r d` candidate angle vectors of Theorem B from a basis of `P0`, excludes candidates by interval evaluation; survivors = `{phi_0}` at `rand4`/31, `bal4`/41, `last4`/31, `first4`/83, 107 (`K = 4`), 163 (`K = 5`) | interval |
| `triple_gen.py FFILE K QMAX ps` | exact test `Lambda_T = Z^r` (Smith normal form) per forced `q0`; outputs `tg_<set>_K<K>.txt` (`K = 4, 5`, four sets): index 1 at every forced `q0 >= 7` | exact |
| `scan_onemod.sh` | one-modulus scans (classify + phi_pin) per forced `q0`; outputs `scan1_<set>_K4.txt` | exact |
| `lift_check.py` | level-2 lifts and all-level re-verification of SL witnesses (PU, (1.2) at `phi_0`, outside every `eps`-class) | exact |
| `strict_check.py`, `scan_strict.sh` | strictness of one-modulus class survivors for global `(1, 0)`, any `(eps, 0)`, any `(eps, t)`; outputs `strict_first4_K4.txt`, `strict_first4_K5.txt`, also all-moduli-charged runs | exact |
| `affordable.py` | per modulus, numbers of PU / ANG / cofactor-affordable single extras on non-genuine `v` of norm `<= K-1` (`afford_first4_K5.txt`, `afford_allnorm_first4_K5.txt`) | exact |
| `scan_allnorm.sh` | one-modulus scan with **every** odd prime `< 3000` not dividing `N` charged (the all-moduli reference), `first4`, `K = 4` (`scan_allnorm_K4.txt`) | exact |
| `triple_allnorm.py` | TC of Theorem P for all-moduli charging (`tc_allnorm_K3.txt`), or with `Y`-smooth `Re(Z_v^2)` (`tc_Y.txt`, `tc_smallY.txt`) | exact |

Stored outputs of the §2.3 table: `cls_first4_K4_83.txt`, `cls_first4_K4_107.txt`, `cls_first4_K5_163.txt`,
`cls_rand4_K4_31.txt`, `cls_bal4_K4_41.txt`, `cls_last4_K4_31.txt`; interval runs `iv_all_rand4.txt`, `iv_all_bal4.txt`,
`sanity_iv_out.txt`; Theorem B enumerations `phi_pin_out.txt`; first-flexible witnesses `cls_first4_K5_167.txt`,
`cls_first4_K5_173.txt`, `lift_out.txt`, `lift_first4_K5_167.txt`.

**Further exact data (not used in proofs).**
* *Affordable extras* (`first4`, `K = 5`, forced set): per covered modulus the numbers of single extras on
  non-genuine `v` of norm `<= 4` that are PU- / ANG- / cofactor-affordable fall from `125 / 92 / 77` (`q0 = 11`) to
  `37 / 24 / 22` (`q0 = 113`); with every modulus charged, `110 / 73 / 62` to `20 / 6 / 5`, and `<= 3` ANG-affordable
  from 127 on. SL removes 30-70% of the PU-affordable extras (almost all of the archimedean ones: ANG minus C is
  `<= 15`), but the local certificate of rigtheory Cor. 3.5(2) (no affordable extra in `D - D`) still never applies:
  rigidity at a covered modulus remains a counting fact in SL too.
* *All-moduli reference, one modulus, `first4`, `K = 4`* (every odd prime `< 3000` not dividing `N` charged):
  SL-rigid (up to `eps`) at every `q0 <= 113` (109 is PU-flexible but SL-rigid), first modulus with level-1 labels outside every `eps`-class that are SL at `phi_0`: 127 (512;
  not lifted to level 2 in this run), then 137 (128), 139, 149, 151, 157, ... (ANG and PU counts in the output file). Compare the forced-set reference: SL-rigid
  `[3, 103]`. Charging the free moduli moves the one-modulus SL segment from 103 to 113 (rigtheory §5.7 saw the same
  effect in PU).

## 6. Open points

1. **Joint SL rigidity (the core question).** Theorems A, B, C and P are one-modulus or conditional-joint statements.
   No theorem here says that every Rig^SL(K) system is coherent on a growing range. What is missing is a counting
   argument: by Theorem C (transport) an A0-type deviation forces a second deviation at a divisor of `Re(Z_v^2)` or
   `Im(Z_v^2)`, and by Corollary B2 an angle move forces dropped odd collisions on a set that breaks every generating
   triple family; one must show that these chains cannot close up within the cofactor slack. In the dense regime
   (`n_q <= r`) this is a question about index-`n` lattices whose short vectors avoid the (sparse) affordable set.
2. **Loose sets in the joint SL problem.** One-modulus SL segments are now exact (`rand4` `[3,47]`, `bal4` `[3,41]`,
   `last4` `[3,29]` at `K = 4`), but the joint SL core of the loose sets is only known to be `>= 7` (PU).
   Theorem B suggests a joint certificate scheme: fix `phi = phi_0` (GA class) whenever the odd collisions of a
   generating triple family are retained, and treat drops explicitly.
3. **Strictness.** Strictness up to a global `(eps, t)` is not implied by coherence mod `n`; one-modulus data show
   local twists surviving SL at `K = 4` (79, 89, 103) and at `K = 5` only from 157 on. A proof of strictness, or the
   mixed-twist version of H^gen that would make it unnecessary, is open.
4. **Covering and TC in the regime.** Covering of `Z/n_q` by the genuine labels of `S_h` (needed for rigtheory Thm
   3.4) and the triple condition TC (needed for Theorem P) are statements about the specific Gaussian primes. Both
   hold in every computed instance (TC: first `r <= 64` split primes, loose sets) and are heuristically generic, but
   neither is proved for the prime sets of the regime. Both are shared by fakes (they concern genuine primes only).
5. **The pair conditions.** Everything here uses only the monomial conditions (e-all). The pair conditions (e-char)
   with point residues are where local twists change the pair collision units (§4.1 item 3); their interaction with
   SL rigidity is unexplored.
6. **Growth of actual circles.** Untouched: `M <= (2/3 + eps) log R/loglog R` stands; the remaining inputs in every
   branch (H^gen + L1) are statements about genuine Gaussian integers.
