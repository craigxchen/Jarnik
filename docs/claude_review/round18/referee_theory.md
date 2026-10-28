# Referee report on round18/theory.md

Referee scripts and outputs: `round18/referee_theory_checks/`. All were written independently; none imports the author's code.

## 0. Verdict

**The result survives in its stated weak form: partial, neither lane decides `(L_8)`, and no
growth or uniformity is claimed.** I found no fatal error.

* **What holds.**
  * The object is `(L_8)` of `ptolemy.md` §7 as used in Steps 0-5 (§1).
  * Every formal proof checks: Lemma 1.1, Prop 2.1 as a reduction, Lemma 2.2, Prop 2.3,
    Lemma 2.4, Theorem 2.5, Prop 2.6 with its identity, the morphism bound, the Mason bound, the
    expected dimensions and the F-curve degrees (§2).
  * Every computation I re-ran with my own exact code is reproduced exactly (§3).
* **The one major defect is interpretive (item B).** The refutation lane says that refuting
  `(L_8)` is "exactly" refuting Vojta's `D = 0` inequality at balanced points, and §3.1 says
  "iff". Only one direction is proved. Non-caterpillar balanced points would refute `V_8` and not
  `(L_8)`. So the chain hypothesis, which is vacuous for the proof lane (Prop 2.1 is right about
  that), constrains the refutation lane.
* **One minor-to-major computational overclaim (item F).** "All unit pencils of degree `<= 4`" is
  really "all unit pencils `s A/C` with `s = ±1`". Pencils `lambda A/C` with other rational
  `lambda` are not covered by the author's data, and some have more multiple fibres (`4z(1-z)`). I tested every special
  rational `lambda` for `D <= 3`, and the generic value too: none beats `-32` (§3.4). `D = 4` with `lambda ≠ ±1` remains untested.
* **Six minor items (A, C, D, E, G, H)** of overstated wording; table in §2.9.
* **No refutation of (L_8) is claimed, so no Zariski-density check is needed.** No owner
  constraint is violated (§4).

## 1. Object check: is the object exactly (L_8) of ptolemy.md §7 / Steps 0-5?

**Yes, with two small caveats about how far the reductions reach (§2, items A and B).** I checked
the object item by item against `ptolemy.md` §1, Theorem B, §6.2 and §7.

* **Anchoring and residues.** The write-up uses `P_0 = (1,0)`, `Y_j = det(P_0,P_j)`,
  `t_ij = det(P_i,P_j)/prod_{T∋i,j} n_T` and `H_Sigma` with `7 + 21` terms. This is `ptolemy.md`
  §1 (with `P_0 = 1`) and §7.
* **Nested chains.** `ptolemy.md` §1 says: "at one rational prime, the blocks with `p | n_T` form
  two nested chains on disjoint row sets". This is the write-up's condition (a) in Lemma 1.1,
  word for word (sides not containing the anchor).
* **Balance.** `ptolemy.md` §7 writes `|log n_T - w| <= eta_0 w` without saying which `T`. The
  write-up imposes it on the 119 visible blocks only.
  * That makes the hypothesis weaker, so the lemma is stronger.
  * It is equivalent to the all-`T` version: singletons occur in no minor, and the top block
    can be inflated by `P_j -> (n X_j, Y_j)` with a fresh prime `n`, a Möbius map `z -> nz`
    that fixes `Phi`. A fresh prime that divides a single block is a one-element chain, so the
    chain condition survives. This is the argument of `census.md` §1.1, which the census referee
    has already accepted.
  * Step 1 gives balance on all nonempty `T`, so it feeds Steps 0-5 either way.
* **Exceptional set.** `Z` is quantified before the data. Lemma 6.1 removes it at cost
  `c(Z) M^7` in Step 2, uniformly in `R`. No fixed-`S` exceptional set enters. ✓
* **eta_0, O(1), w.** One fixed `eta_0` suffices, because `eta(M) -> 0` in Step 1. The `O(1)`
  enters only through `w >= C'`, and `w >= (M-1) lambda/32`. ✓
* **The step from (L_8) to the census object.** The write-up passes to contact profiles through
  `census.md` Prop 1.3. I re-read that proposition. Its part (2), "`(L_8)` ⇒ `(L_8^c)`", holds
  **only for admissible (caterpillar) `Phi`**. Two global statements in the write-up drop this
  restriction (items A and B of §2).

So the census is of the right object. The defects are in two "iff/equivalent" claims that
overstate what the reduction gives.

## 2. Item-by-item re-derivation

I re-derived every proof step by hand (derivations in `referee_theory_checks/notes.md`). Items
marked ✓ are correct as stated. Lettered items A-H are the defects, collected in §2.9.

### 2.1 Lemma 1.1 (caterpillar lemma) ✓

* (c) ⇒ (a) for any anchor: correct.
* (a) ⇒ (c): correct. The key step is that an internal edge shared by the two paths towards `0`
  would have a side meeting both chain unions. That side is itself in `Sigma`, so it lies in one
  chain, which contradicts the disjointness of the unions.
* My own code (§3) reproduces the counts 26/236/2752/39208 and 26/221/2227/25698 with 0
  mismatches, using a fourth, independent criterion: every vertex meets at most 2 internal edges.

*Scope (minor, item H).* The lemma is about split systems. The (L_8) hypothesis is about the
**blocks** of a datum. Blocks equal the `p`-adic tree edges only at primes that divide no residue
(Thm B(iv)(2)). So "the chain condition is the caterpillar condition on `Phi`, a condition on
`Phi` alone" is exact only for such primes, in particular for `H_Sigma = 0` data. At primes that
divide residues the tree of `Phi` need not be a caterpillar.

### 2.2 Proposition 2.1 (joint content) ✓ as a reduction; the gloss is overstated

* (1)-(3) are the cited Theorems B(i)-(iv) and Lemma 5.2. (4) is immediate: on `Bal^pc` each
  prime divides exactly one contact.
* `(L_8) ⇒ (L_8^pc)` is `census.md` Prop 1.3(2). Pairwise-coprime points are admissible. ✓
* The sentences "no structure of the joint system other than `K.beta = 10` survives" and "none of
  it yields a positive multiple of `w`" are **not theorems**. What is proved is a reduction: on
  `Bal^pc` the chain, anchoring and residue hypotheses are automatically satisfied. What remains is
  not "the number `K.beta = 10`" but the full arithmetic of a rational point with pairwise-coprime
  contacts, that is, the lattice problem of §2.2. (Minor, item G.)

### 2.3 Rate dichotomy (§1.3) — item A

The proof shows that the `f`-version (`c w` replaced by `f(w) -> infinity`), applied to
`H_Sigma = 0` data, gives `(L_8^c)(eta_0, Z)` **for caterpillar `Phi`**. What is actually proved
is the chain

```text
(L_8^c)^all  ⇒  (L_8)  ⇒  (L_8^f) (f <= c w)  ⇒  (L_8^c)^cat ,
```

and nothing closes it. Data with `0 < H_Sigma <= c w` can represent **non-caterpillar** points.
The residues can carry the non-caterpillar part of the tree, at total cost `O(H_Sigma)`. So
`(L_8^c)^cat` does not visibly imply `(L_8)`. "No easier than `(L_8)`" and "the only freedom is in
the constants" are therefore overstated.

The correct statement is: *any proof of the `f`-version must prove `(L_8^c)^cat`, which is also
the whole `H_Sigma = 0` content of `(L_8)`.* This still kills slowly growing `f` as a shortcut,
which was the point. (Minor.)

### 2.4 Lemma 2.2, Prop 2.3, heuristic ✓ (with item D)

* **Index.** Going from level `t-1` to level `t`, a class `C` with `|C| >= 3` adds `|C|-1`
  independent conditions mod `p`. `T` is a class at exactly `l_T(p)` levels. So the index is
  `prod d_T^{|T|-1}`. This is confirmed by my HNF computation (§3).
* **Kapranov's canonical class.** `K = -6H + sum (|T|-2) E_T`, since `L_T` has codimension `|T|-1`
  in `P^5`. The checks `H.beta = 31` (the visible `T ⊇ {i,j}`), `196`, `294` and
  `K.beta = 196 - 186 = 10` are all correct. Minkowski scale `e^{49w}`, `|x| = e^{31w}`, trivial
  bound `e^{30w}` (30 = 31 - 1: drop the pair itself): correct.
* **Prop 2.3** is correct arithmetic, but it bounds `kappa`, the **finite part** of `h_K`, and not
  `h_K` itself. Summary item (4) and §3.3(2) write `h_K <= -sum ell_T`, which omits the
  archimedean terms `(|T|-2) lambda_{E_T,inf}` (minor, item D).
* **Chart.** Prop 2.3 is proved in the anchored Kapranov chart. The `(L_5)` refutation (`L8.md`
  Prop C) uses the integral chart, with `infinity` unmarked. I checked by hand that the conclusion
  survives there for balanced points at `k = 8`. Let `C_S` be the clustering side of split `S`,
  `a` the number of `2|6` splits clustered on the 6-side, and `b` the number of `3|5` splits
  clustered on the 5-side. Positive density needs `7 log X >= sum(|C_S|-1) w`, and balance forces
  `log X = w · sum C(|C_S|,2)/28`. Together these give
  `101.5 + 3.5a + 1.75b >= 245 + 4a + 2b`, which is impossible. So the meta-claim "every density
  mechanism fails" is right here, but Prop 2.3 as written covers one chart only.
* The heuristic `98 - 108 = -10` is labelled as a heuristic. ✓

### 2.5 Lemma 2.4 and Theorem 2.5 ✓

* **Lemma 2.4.** Blocks at `p` are nested or disjoint (two chains with disjoint unions). Two
  distinct `2|2` restrictions to a 4-set are incompatible. ✓ (It is a trivial lemma.)
* **Theorem 2.5.** I re-derived it from (4.1) of `ptolemy.md`, namely
  `Y_i t_jl N_jl - Y_j t_il N_il + Y_l t_ij N_ij = 0`.
  * Each residue occurs once, so `|abc| <= e^{H_Sigma}`.
  * `|{T : T ∩ {i,j,l} = {j,l}}| = 16`, all visible.
  * `g | abc`, since each prime divides at most one `N`.
  * The Ridout step follows from the Subspace Theorem for `n = 2`. At each `p ∈ S_0` only one of
    the coprime `x, y, z` is divisible by `p`, so `[xyz]_{S_0} <= C H^eps |L_1 L_2|_inf`. Hence
    the `S_0`-free part is `>= H^{1-eps}/C`. The finitely many exceptional ratios are absorbed
    into `C`.
  * Combining gives the displayed bound and `8(1-eta)(1-eps)` when `Xi = 0`. ✓
* `C(S_0)` is not claimed uniform, and the theorem is declared useless for the application. This
  respects the owner's first constraint.

### 2.6 abc lane: Prop 2.6, identity, partial closure ✓ (with item F)

* **Prop 2.6** is correct. I checked the three delicate points.
  1. `Y ∈ I_{L_T}^{nu}` over `Z`, because `y_i = x_i - x_{i0}` is a unimodular change of
     variables. So `prod d_T^{nu_T(Y)} | Y(x)`, and the cofactor is an integer.
  2. A prime `p | d_T` with `max nu = mu_T` can still divide `b`. It then divides an interior
     cofactor and is counted there.
  3. `gcd(A,B,C) = gcd(A,C) = prod d_T^{mu_T}` exactly on `Bal^pc`, so
     `log max(|a|,|c|) = h(f(Phi)) = (L.beta) w + O(eta w)`.
* **The identity.** `K^log = K + F_red - L` gives `gain = 10 - K^log.beta~`. The multiplicity
  form `gain = sum ((m_Y - 3)/3) Y.beta` follows from `3L ≡ F_0 + F_1 + F_inf`. ✓
* **Morphisms.** The Stein factor is a forgetful map (Bruno–Mella, cited from the abstract only,
  and labelled so), so `L.beta = 16 deg h`. Riemann–Hurwitz gives `#h^{-1}{0,1,inf} >= deg h + 2`,
  hence `gain <= -32`. ✓
* **Mason on a general line.** Re-derived: `D <= #(distinct roots) - 2`, with the pairs of
  `A, B, C` disjoint. This gives the displayed bound with `19 = 21 - 2`. ✓
* **Item F (the enumeration is not "all unit pencils").** The scripts take `f = s A/C` with
  `s = ±1` only. A unit pencil is `lambda A/C` with any `lambda ∈ Q^*`. abc applies equally to
  `q C - p A = q B`, because the constants change `rad` by `O(1)`. Special values
  `lambda ≠ ±1` give genuinely different, more multiple fibres. Example:
  `f = 4z(1-z)` with `z = (x_1-x_2)/(x_1-x_3)` gives `B = (x_1 - 2x_2 + x_3)^2`. That pencil is
  **not** among the 277940. So "all unit pencils of degree `<= 4` are excluded" is overstated: the
  data cover `lambda = ±1`. I closed this gap independently for `D <= 3` (§3.4). Every rational
  special `lambda` gives gain `<= -32`.
* **Strengthening (in the author's favour).** The `nu_T` values are computed by generic
  specialization, so they can only be **over**-estimated. An overestimate of `nu_T(Y)` lowers the
  interior `beta`-degree by at least as much as it can add a boundary component. So the computed
  gain is always `>=` the true gain, and the data give a **rigorous upper bound** `gain <= -32`
  (given exact factorization over `Q`). The label "Data, not proof" can be upgraded, provided the
  enumeration is complete. My own exact implementation (§3) needs no such argument.

### 2.7 Refutation lane, §3.1-§3.4 — items B, C, E

* **Expected dimensions.** Re-derived:
  * genus 0: `186d + 2 - 196d = -10d + 2`;
  * genus 1: `186d - 196d = -10d`;
  * `k = 7`: `75d + 1 - 68d = 7d + 1`.

  Note that `L8.md` §7 already records `-10d + 2`. ✓
* **F-curves.** Re-derived with `c(s) = s(8-s)/7 - 2`. `K.F` is `2, 1, 0, 0, -1` for types
  `(2,2,2,2)`, `(3,2,2,1)`, `(3,3,1,1)`, `(4,2,1,1)`, `(5,1,1,1)`. ✓
* **Covering families of genus `<= 1`.** By the normal-bundle argument, `K.gamma <= 2g - 2`, and
  `K` is positive near the `beta`-ray. So no genus-`<= 1` curve carrying infinitely many balanced
  points lies in a covering family. ✓
* **Item C.** "Covering families of curves are excluded (proved)" is wrong for genus `>= 2`.
  Faltings gives finitely many points per member, but a covering family of genus-`>= 2` curves
  with finitely many balanced points on each member can still give a Zariski-dense set. (So can
  `Q`-lines of `P^5`, which are a covering family of genus 0 through every point.) The census
  referee raised the same issue. Here it is partly repaired, because the residual category
  "sporadic points" absorbs these cases. The word "excluded" must still be restricted to
  "covering families whose members carry infinitely many balanced points".
* **Item E.** "Rigid curves" should read **non-covering** curves: curves that move only in
  families sweeping out proper subvarieties are also allowed. The two-way split (curves with
  infinitely many points / sporadic points) is a tautological case split, not a theorem. It is
  harmless once it is labelled as a description.
* **Item B (major, interpretive).** §3.1 states: "`(L_8)` fails **iff** for every `eta` and
  `w_1`, `Bal(eta) ∩ {w >= w_1}` is Zariski dense". §3.4 and §5(1) state: "a refutation of
  `(L_8)` is **exactly** a refutation of `D = 0` Vojta at balanced points". Only one direction
  is proved:

  ```text
  Bal^cat(eta) dense for all eta  ⇒  not (L_8)  ⇒  Bal(eta) dense for all eta  ⇒  not V_8 .
  ```

  A Zariski-dense set of **non-caterpillar** balanced points would refute `V_8` but not `(L_8)`.
  An `(L_8)`-datum for such a point needs the residues to absorb the non-caterpillar part of the
  trees, at cost of order `w`. So the chain hypothesis is vacuous for the **proof** lane (on
  `Bal^pc`), as Prop 2.1 correctly says, but **not** for the refutation lane. A refuter must
  produce caterpillar points, up to `o(w)` of residue.

  "`(L_8)` cannot be mis-stated in a way that a refutation could exploit" is true: `V_8 ⇒ (L_8)`
  is refereed. The converse, which is what "exactly" asserts, is false as stated. The correct
  statement: refuting `(L_8)` refutes `V_8` at balanced points; it is **harder** than that, since
  it needs (near-)caterpillar balanced points.
* **Formulation audit, items `O(1)`, `eta_0`, `Z`, invisible blocks.** ✓ (see §1).

### 2.8 Precise residual §5 ✓ with wording fixes

* §5(2) calls `(L_8^pc)` "equivalent in substance to `(L_8)`". It is a strict consequence of
  `(L_8)`; the parenthesis concedes this.
* The lattice set in §5(2) is slightly larger than `Bal^pc`: pair cofactors may share a prime with
  `e_T` for `T ⊇ {i,j}`. I checked that every point of it is in `Bal^cat` with contacts
  `d_T = e_T` for `|T| >= 3`. So it is still a consequence of `(L_8)`. ✓ (Item G.)
* §5(3) is correctly labelled conditional.

### 2.9 List of defects

| item | where | severity | defect | fix |
|---|---|---|---|---|
| A | §1.3, summary (3) | minor | "f-version no easier than (L_8)" | it implies `(L_8^c)^cat` only |
| B | §3.1, §3.4, §5(1), summary (8) | major (interpretive) | "iff" / "exactly a refutation of V_8" | one direction; refutation needs caterpillar points |
| C | §3.3(1) | minor | genus `>= 2` covering families "excluded" | restrict to members with infinitely many points |
| D | summary (4), Prop 2.3, §3.3(2) | minor | `h_K` written for `kappa` | finite part only |
| E | §3.3(4) | minor | "rigid" | "non-covering"; case split is a description |
| F | §2.4(2), summary (7), §5(3) | minor-major | "all unit pencils of degree `<= 4`" | only `lambda = ±1`; referee closes `D <= 3` |
| G | Prop 2.1 gloss, §5(2) | minor | "no structure survives" / "equivalent in substance" | reduction only |
| H | Lemma 1.1 application | minor | chain = caterpillar "on `Phi` alone" | only at primes dividing no residue |

## 3. Independent computations

All scripts are in `referee_theory_checks/`. They were written from scratch and import none of
the author's code. Arithmetic is exact (Python integers, `flint` `fmpz_mpoly` factorization and
composition, HNF). Every job ran under `prlimit --as=3500000000`.

### 3.1 Lemma 1.1 — `ref_caterpillar.py` → `ref_caterpillar.out`

* All compatible systems of nontrivial splits were enumerated by backtracking over the
  compatibility graph.
* Tested: (a) for anchor 0, (a) for all anchors, (c) "one side per split forms a chain", and an
  independent criterion (d) "every vertex of the rooted tree meets at most 2 internal edges".

| leaves | trees | caterpillars | mismatches among (a), (b), (c), (d) |
|---|---|---|---|
| 5 | 26 | 26 | 0 |
| 6 | 236 | 221 | 0 |
| 7 | 2752 | 2227 | 0 |
| 8 | 39208 | 25698 | 0 |

**Reproduced exactly.** (The tree counts are the Schröder numbers A000311, and
`236 - 221 = 15` snowflake trees.)

### 3.2 Lemma 2.2(1) and (1.1) — `ref_lattice.py` → `ref_lattice.out`

* 80 primitive configurations from three generators, including strongly clustered ones.
* Contacts were computed from residue classes mod `p^t`.
* The index of `Lambda(x)` was computed as `|image of Z^7 in ⊕ Z/p^t|` via an HNF of
  `[M | diag(p^t)]`.
* **80/80**: the index equals `prod_{|T|>=3} d_T^{|T|-1}` and
  `|x_i - x_j| = prod_{T∋i,j} d_T`. **Confirmed.**

### 3.3 abc gain, `lambda = ±1` — `ref_gain.py`, `ref_gain_chunk.py`, `run_D4.sh`

* `nu_T` of each interior factor is computed **exactly**. The substitution
  `x_i -> x_{i0} + e x_i` (`i ∈ T \ {i0}`) is applied symbolically, and the lowest power of `e`
  is read off, with no specialization.
* Factorization is over `Q` (`flint`).
* `A`-orbits under `S_7` are computed by brute-force canonical relabelling.

| D | A-orbits | pencils | top of gain histogram | author |
|---|---|---|---|---|
| 1 | 1 | 40 | `-32` (20), `-48` (20) | same |
| 2 | 3 | 1180 | `-32` (24), `-40` (12), `-48` (152), `-52` (144), `-56` (272) | same |
| 3 | 8 | 19800 | `-32` (20), `-40` (4), `-44` (2), `-48` (208), `-52` (188), `-56` (876) | same |
| 4 | 22 | 256920 | `-32` (24), `-40` (14), `-48` (318), `-50` (10), `-52` (238), `-54` (6), `-56` (1636) | same |

All four rows agree with the author's output, histogram entry by histogram entry. The maximum `-32` is attained, for example, by `((x_1-x_2)/(x_1-x_3))^4`. Outputs: `ref_gain_D12.out`, `ref_gain_D3.out`, `ref_gain_D4_chunks.jsonl`, `ref_gain_D4_summary.out`. (The single-process `D = 4` run ran out of memory. The cause was a growing in-process memory
footprint, not the polynomials. It was therefore split into 66 processes by orbit and by
`C`-chunk.)

### 3.4 abc gain, other `lambda` — `ref_gain_lambda.py` → `ref_gain_lambda_D*.out` (new)

For each `(A, C)` with `D <= 3` I computed every rational `lambda ≠ 0, ±1` at which
`B = C - lambda A` can be special:

* (i) the initial forms of `A` and `C` along some `L_T` (`|T| >= 3`, `nu_T(A) = nu_T(C)`) are
  proportional with ratio `lambda`, so that `nu_T(B) > mu_T`;
* (ii) the restrictions of `A` and `C` to some hyperplane `x_i = x_j` are proportional, so that
  `B` acquires a pair factor;
* (iii) `lambda` is a rational root of the gcd of the `t`-discriminants of `c(t) - lambda a(t)` on
  two random lines. This is a necessary condition for `B` to be non-reduced.

The exact gain was then evaluated at every candidate. Outside this candidate set, `B` is reduced,
has no extra pair factor and no extra boundary order. Its gain is then the generic value, which
can only be lower.

**Result (Data; `ref_gain_lambda_D2.out`, `ref_gain_lambda_D3.out`).**

* All `20 + 590 + 9900 = 10510` pairs `(A, C)` with `D <= 3` were processed.
* 17 special rational `lambda ≠ ±1` occur (6 at `D = 2`, 11 at `D = 3`). Example:
  `lambda = -1/4` for `A = (x_1-x_2)^2`, `C = (x_1-x_3)(x_2-x_3)`, which gives
  `B ∝ (x_1 + x_2 - 2x_3)^2`, the pencil `4z(1-z)` in another guise. Every one of them has
  gain exactly `-32`.
* The generic value. For `lambda` outside the candidate set, `B` is reduced, has no new pair
  factor, and has `nu_T(B) = mu_T`. The gain is then exactly
  `-#{T visible : nu_T(A) ≠ nu_T(C)}`. Its maximum over all 10510 pairs is `-32`.
* So **for `D <= 3`, no unit pencil `lambda A/C` with `lambda ∈ Q^*` has gain above `-32`.**
* Caveats. Step (iii) uses two random lines (fixed seed), so the candidate set is complete only
  with probability 1. A missed `lambda` would need both lines to be special, so this is labelled
  Data. `D = 4` with `lambda ≠ ±1` was not run.

## 4. Evidence vs proof; owner constraints

* **Labels.** Proved, Data and Heuristic are used correctly throughout, with two exceptions.
  * The words "excluded (proved)" in §3.3(1)-(2) and "exactly" in §3.4 claim more than the
    proofs give (items B, C, D).
  * The enumeration is called "Data", which is accurate. By §2.6 it can be upgraded to a
    rigorous upper bound for the pencils it covers, but it does not cover all unit pencils
    (item F).
* **Bruno–Mella** is cited from the abstract only, and the write-up says so. I could not read the
  text either (arXiv is blocked), so the morphism exclusion is conditional on the abstract being
  accurate. Note also that `Aut(Mbar_{0,4}) = PGL_2`, so the factorization is "forgetful map
  followed by an isomorphism `Mbar_{0,4} ≅ P^1`". That isomorphism is absorbed into `h`, and the
  `-32` bound is unaffected.
* **Owner constraints.**
  1. *No uniform exceptional sets from fixed-`S` Subspace.* Thm 2.5 has `C = C(S_0)`, is not
     claimed uniform, and is declared useless. ✓
  2. *No discarded Gaussian gcd or denominator costs.* `(L_8)` is form-blind; Step 0 (Gaussian
     gcd) is untouched; `H_Sigma` keeps all 28 residues. ✓
  3. *No drifting common phase replaced by 0.* No phase enters. ✓
  4. *No template-dependent threshold treated as uniform.* Prop 2.6 has `w_1` depending only on
     the fixed pencil `f`, and `Z = {B = 0}` is fixed. The cited `census.md` Prop 7.5(4) has
     curve-dependent constants, and the write-up does not use it uniformly. ✓
* **No growth or uniformity is claimed.** ✓

## 5. Corrected statement

**Status: partial. Neither lane decides `(L_8)`. No growth improvement and no uniform bound are
claimed or proved.**

Proved (all quantifiers checked by the referee):

1. **Lemma 1.1.** For a compatible split system on `{0..k-1}`, the item-532 two-chain condition
   for one anchor, the same condition for every anchor, and the caterpillar condition are
   equivalent. Applied to an `(L_8)`-datum, it says the chain hypothesis is the caterpillar
   condition on the **block system** at each prime. That is the `p`-adic tree of `Phi` at primes
   dividing no residue, and in particular for `H_Sigma = 0` data.
2. **Prop 2.1 (reduction).** Every proof of `(L_8)` must prove `(L_8^pc)`. On pairwise-coprime
   profiles the chain, anchoring and residue hypotheses are automatically satisfied. This is a
   reduction, not a theorem that "no structure survives".
3. **Rate statement.** For every `f -> infinity`,
   `(L_8^c)^all ⇒ (L_8) ⇒ (L_8^f)` (for `f <= c w`) `⇒ (L_8^c)^cat`. So any slowly-growing
   version still has to exclude dense caterpillar balanced points. It is **not** shown to be
   equivalent to `(L_8)`.
4. **Lattice form (Lemma 2.2).**
   * `[Z^7 : Lambda(x)] = prod_{|T|>=3} d_T^{|T|-1}`.
   * `kappa(x) = sum (|T|-2) ell_T - 6 log|x|`, the finite part of `h_K`.
   * A balanced point is a vector of `Lambda(x)` at `e^{31w}`, against Minkowski `e^{49w}` and the
     trivial bound `e^{30w}`.
   * **Prop 2.3.** In the Kapranov chart, at positive lattice-point density, `kappa <= -sum log e_T`.
     This bounds the finite part, not `h_K`. The analogous count in the integral chart also
     forbids balanced points at `k = 8`.
   * The count `98 - 108 = -10` is a heuristic.
5. **Lemma 2.4** (matching coefficients pairwise coprime) and **Theorem 2.5** (fixed `S_0`,
   Ridout: `H_Sigma >= 8(1-eta)(1-eps) w - C(S_0, eta, eps)` when all blocks are `S_0`-units; no
   exceptional set). `C` depends on `S_0`, so the theorem is useless for the circle problem.
6. **abc lane (conditional on abc).**
   * **Prop 2.6.** A unit pencil `lambda A/C` (`lambda ∈ Q^*`) with
     `gain = L.beta - F_red.beta = 10 - K^log_rel.beta > 0` would give `(L_8^pc)`.
   * Morphisms have gain `<= -32` (Bruno–Mella, abstract only).
   * **Data.** Unit pencils with `lambda = ±1` and `D <= 4` (277940 pencils): maximum gain
     `-32`. All rational `lambda` special for some pencil with `D <= 3`: maximum gain `-32`
     (referee). Pencils with other `lambda` and `D = 4` are untested.
   * Since specialization only overestimates `nu_T`, these data are rigorous upper bounds for the
     pencils enumerated.
   * The abc barrier of the record does not cover such pencils. This gap in the barrier map is
     genuine.
7. **Refutation lane.**
   * Expected dimensions `-10d + 2` and `-10d` (the first already in `L8.md` §7).
   * No genus-`<= 1` curve with infinitely many balanced points is a member of a covering family.
   * The rigid curves with `K > 0` that were computed (F-curves) lie in the boundary.
   * A refutation of `(L_8)` refutes `V_8` at balanced points. The converse fails as stated: a
     refutation of `(L_8)` needs a Zariski-dense supply of balanced points that are caterpillar at
     each prime up to `o(w)` of residue. Non-caterpillar balanced points refute only `V_8`.
   * What such a supply would have to look like (non-covering genus-`<= 1` curves carrying
     infinitely many points, or points with finitely many per curve) is a case split, not a
     theorem.
