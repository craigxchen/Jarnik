# Round 23 — referee report on judge check J4 (absolute form of the A_5 reduction)

Object: `round23/README.md` §5 J4 (Proposition J4 and its 7-step proof, "Status and meaning") and §7 item (3)
(proposed expert-note text). Cited inputs read: `round20/gauss.md` §§1.2–1.3, 3.1–3.2 (Lemma 1.2, Lemma 3.2,
Lemma 3.3, Thms 3.4, 3.4′), `round20/kill_gauss.md` §§2, 7, `round23/split.md` §3.4 (Thm 1) and §6.2,
`round23/referee_split.md`, `expert_note/expert_note.md` §§1, 3, `overview/overview.md` §1 (definition of `M`,
`M_∞`), `growth/ptolemy.md` §6.2 (Steps 0–5), `lit18/direct.md` §2.10 and the record's copy of
Cilleruelo–Granville Conjecture 14. Scripts and outputs: `referee_J4_checks/` (absolute path at the end). The
repository was read only.

## 0. Verdict

**Proposition J4 survives as stated. The proof is correct, with all quantifiers.** `K` depends only on `δ`
(through `η`) and on `Z` (through the degree `D` of one form vanishing on the cone over `Z`). `C` enters only
through `R_1(C)`. The position of the arc and the radius do not enter `K`. `g_S` is divided out exactly, the
reduced radius is exact, and units play no role. I re-derived every step (§2) and tested the identities exactly
on actual clusters (§4): 795,172 five-subsets, 0 failures. The tests include ramified `(1+i)`, inert `q`,
non-normalised clusters and cross-unit-class subsets, which the judge's own check did not cover.

**The gloss around it overclaims in one place (major, text only), and three remarks need correcting (minor).**
1. *(major, text)* "Question B is at least as strong as …" and "would refute Question B". Question B is
   "prove `(L_5^∘)` **or** `A_5(δ)`, `δ < 1/30`". J4 proves the absolute form only from `A_5`. The `(L_5^∘)`
   branch goes through `gauss.md` Thm 3.4′ / `ptolemy.md` Step 0. Step 0 subdivides the arc into `⌈C⌉` pieces
   and uses `C_0 < √2` in the residue bound, so its `K` is **not** shown to be `C`-independent. Correct the
   statement to `A_5(δ)` (and `V_E(5)`).
2. *(minor)* "strictly stronger than the owner's goal `M(C) < ∞`". It is **formally** stronger; strictness is
   not known.
3. *(minor, missing context)* The absolute form `sup_C M_∞(C) < ∞` is a **consequence of the published
   Cilleruelo–Granville Conjecture 14**: arcs of length `R^{1-ε}` hold `≤ B_ε` points; take `ε = 1/4` and
   `R ≥ C^4`. So J4 places the A_5 branch of Question B above a statement that sits between the critical
   conjecture and CG's own `R^{1-ε}` conjecture. The "circle analogue of #887 `parts.ii`" reading is fine. The
   sentence on Oganesyan's withdrawn `R^α` claim is irrelevant at the `√R` scale (`lit18/direct.md` §2.10
   already says so) and is best dropped or rephrased.
4. *(minor, understated)* "What is new is **only** the `C`-independence of `K`" undersells the result. The same
   argument gives the **sharp** exponent `δ < 1/(2^k-2)` for every `k`, by convexity averaging alone, with no
   balanced selection and no unit-class restriction. This corrects the expert note's "Reduction" bullet, which
   says that balanced selection within one unit class is needed beyond `δ ≤ 2^{-k-1}` (§3.3).

Growth or uniform progress: **none.** J4 is a reduction from an open hypothesis. `(2/3+ε) log R/loglog R`
stands.

## 1. Statement, as refereed

Definitions (`overview.md` §1, `expert_note.md` §§1, 3):
* `M(C) = sup_R max_I #(I ∩ E_R)` and `M_∞(C) = limsup_{R→∞} max_I #(I ∩ E_R)`, where `I` ranges over arcs of
  length `≤ C√R` of `|z| = R`, `N = R² ∈ N`.
* `A_5(δ)`: there are a proper Zariski-closed `Z ⊂ T_5` (over `Q`) and an `R_0` such that every primitive 5-tuple of
  distinct lattice points of a circle of radius `R'' ≥ R_0`, with all mutual distances `≤ R''^{1/2-δ}`, has
  `[ζ] ∈ Z`. There is no unit-class condition and no balance condition.

> **Proposition J4 (refereed; correct as stated).** Let `0 < δ < 1/30` and assume `A_5(δ)` with `(Z, R_0)`. Put
> `η = (1/(1+2δ) - 15/16)/2 > 0`, `κ = 1 - (15/16 + η)(1+2δ) > 0`, and let `D` be the degree of a nonzero homogeneous
> `F ∈ Q(i)[u_1..u_5]` vanishing on the cone over `Z_C`. Let `M_1` be the least integer `≥ 9` with
> `(η/2) C(M,5) > D M^4` for all `M ≥ M_1` (one has `M_1 ≤ 938 D/η`), and put `K = max(⌈9 + 2/η⌉, M_1)`. Then for every
> `C > 0` and every `R ≥ R_1(C) := max(C^{2/κ}, (C R_0)^2)`, every arc of length `≤ C√R` of `|z| = R` contains at most
> `K` lattice points. Hence `M_∞(C) ≤ K` for every `C`, so `sup_C M_∞(C) ≤ K`. Also `M(C) < ∞` for every `C`, but
> `M(C)` is **not** claimed to be `≤ K`: radii `R < R_1(C)` are only finitely many circles, and they contribute a
> `C`-dependent bound.

> **Generalisation (same proof; refereed here).** For every `k ≥ 3`: `A_k(δ)` for one `δ < 1/(2^k-2)` implies
> `sup_C M_∞(C) < ∞`, with `K = K(k, δ, Z)`.

The explicit `R_1(C)` above is mine. The judge's step 5 gets `R'' ≥ R_0` by a compactness constant `c_0`
instead. That is also correct, but less explicit (§2, step 5).

## 2. Step-by-step re-derivation (all quantifiers)

Fix `δ, Z, R_0, η, κ, D, K` as in §1. These depend on `δ` and `Z` only. Let `C > 0`, `R ≥ R_1(C)`, and let
`z_1, …, z_M` be the distinct lattice points of an arc of length `≤ C√R`. Assume `M ≥ ⌈9 + 2/η⌉` (otherwise
`M < K`). Every chord satisfies `|z_a - z_b| ≤` arc length `≤ C√R`.

**Step 1 (factorisation; correct).** `N = 2^{e_0} Π q^{2f_q} Π p^{e_p}`, and every point is
`z_j = ε_j (1+i)^{e_0} Π q^{f_q} Π π_p^{a_j(p)} π̄_p^{e_p - a_j(p)}`. The ramified and inert parts are the same for all
points:
* `(1+i)` is the only prime above 2 and `(1+i) ~ (1-i)`, so `v_{1+i}(z) = v_2(N)`;
* an inert `q` divides `z` exactly `f_q` times.

`ε_j` is arbitrary, so the cluster may meet up to four unit classes; the record's 10-point cluster meets all four
(`c2`).

**Step 2 (exact reduced radius; correct).** `g_S ~ (1+i)^{e_0} Π q^{f_q} Π π^{min_S a} π̄^{e_p - max_S a}`, because
the gcd is computed prime by prime and units do not matter. Hence
`|g_S|² = 2^{e_0} Π q^{2f_q} Π p^{min_S a + e_p - max_S a}` and `R''_S² = N/|g_S|² = Π_p p^{max_S a - min_S a}`.
The reduced tuple `ζ = z_S/g_S`:
* has distinct points;
* is primitive (its gcd is a unit);
* lies on `|ζ|² = R''²`;
* satisfies `[ζ] = [z_S]`.

Exact on 795,172 subsets (`c2`), including ramified, inert, non-normalised and cross-unit-class subsets.

**Step 3 (averaging; correct).** For a split layer `(p, s)`, `1 ≤ s ≤ e_p`, with `κ_s = #{j : a_j(p) ≥ s}`, the
probability that a uniform 5-subset meets both sides is `1 - (C(κ_s,5) + C(M-κ_s,5))/C(M,5) ≤ q(M)`. Here
`q(M) ≤ 15/16 + 15/(16(M-4))` for `M ≥ 9`: this is split Thm 1(ii), and I re-checked it exactly for
`9 ≤ M ≤ 3000` (`c1`).

So `E log R''_S² ≤ q(M) Σ_p e_p log p ≤ q(M) log N`. Ramified, inert and common split factors only help, since
they enter `log N` and not `R''`. The judge's proof says `≤ (15/16 + 15/(16(M-4))) log N`, which is correct.

Put `Y = 1 - log R''²/log N ∈ [0,1]`. Then `E Y ≥ 1/16 - 1/(M-4)` and `P(Y ≥ 1/16 - η) ≥ E Y - (1/16 - η)`. This is
reverse Markov, valid because `Y ≤ 1`. It gives `P ≥ η - 1/(M-4) ≥ η/2` once `M ≥ 9 + 2/η`, checked exactly in
`c1`. So at least `(η/2) C(M,5)` subsets have `R''_S ≤ R^{15/16+η}`.

Also checked exactly, per layer and in aggregate, on 120 actual clusters with `M = 9, 12, 16` (`c3`).

**Step 4 (A-points; correct).** Reduced chords are `≤ C√R/|g_S| = C R''_S R^{-1/2}`, and
`C R'' R^{-1/2} ≤ R''^{1/2-δ}` is equivalent to `C² R''^{1+2δ} ≤ R`. With `R'' ≤ R^{15/16+η}` it suffices that
`C² ≤ R^κ`, i.e. `R ≥ C^{2/κ}`. This threshold depends on `C`, `δ` and `η` only, not on the cluster, the
position or the template.

**Step 5 (`R'' ≥ R_0`; correct, and simpler than written).** By Lemma 1.2,
`|g_S| ≤ min_{a≠b} |z_a - z_b| ≤ C√R`, so `R''_S ≥ √R/C ≥ R_0` once `R ≥ (C R_0)²`. The judge's version
(finitely many primitive tuples with `R'' < R_0`, so a distance `c_0 > 0` from `e`, while
`dist([ζ_S], e) ≤ C R^{-1/2}`) is also correct. `c_0` depends only on `R_0`, hence on the hypothesis, not on the
data. `A_5(δ)` now applies: `[ζ_S] = [z_S] ∈ Z` for every good `S`, in any ordering, because `A_5` quantifies
over all ordered tuples.

**Step 6 (Schwartz–Zippel; correct).** Over `C`, `R^1_{Q(i)/Q} G_m ≅ G_m` via `u = x + iy`. A rational point
`w = a + bi` (`a² + b² = 1`) goes to the complex number `w`, so `[z] ↦ (z_j/z_1)_j` in `(C^*)^5/C^*`. The
preimage `Z'` of `Z_C` in `(C^*)^5` is a proper closed cone. The ideal of its closure in `A^5` is homogeneous and
nonzero, so it contains `F` of degree `D`, and `[z] ∈ Z` iff `F(z_1, …, z_5) = 0` with `z_j` read as complex
numbers. `F` depends only on `Z` (and on the fixed choice `u = x + iy`), not on the circle.

Exact test (`c4`) on three choices of `Z`: a group relation, `Re w_2 = Re w_3`, and `Σ_{j≥2} x_j = 0`. On every
ordered tuple of distinct points of three circles, `[z] ∈ Z ⟺ F(z) = 0`, and the zero counts in `Σ^5` are
`≤ D M^4`.

Distinct good subsets, listed in index order, are distinct points of `Σ^5`, where `Σ = {z_1..z_M} ⊂ C`. So
`(η/2) C(M,5) ≤ D M^4`, hence `M < M_1`.

**Step 7 (conclusion; correct).** `M ≤ max(⌈9 + 2/η⌉ - 1, M_1 - 1) < K`. `K` depends on `(δ, Z)`. `C` appears only
in `R_1(C)`, and `R_0` only in `R_1(C)`. For illustration (`c5`): `δ = 1/40`, `D = 1` gives `η = 5/672` and
`K = 32266`; `δ = 1/31`, `D = 1` gives `η = 1/1056` and `K = 253450`. The constants are large but absolute.

**General `k`.** Replace 5 by `k`. By gauss Lemma 3.2 (refereed),
`1 - ρ_{M,k} ≤ 1 - 2^{1-k}(1 - k(k-1)/(2(M-k)))`. The exponent condition becomes
`(1 - 2^{1-k} + η)(1 + 2δ) < 1`, which holds for small `η` iff `δ < 1/(2^k-2)`: equality at the threshold is exact
for `k = 3..8` in `c5`. Schwartz–Zippel gives `(η/2) C(M,k) ≤ D M^{k-1}`.

## 3. Constraint audit and the claims around J4

### 3.1 Owner constraints

* *Template-dependent threshold treated as uniform in the radius.* **Not violated.**
  * The only thresholds are `R_1(C) = max(C^{2/κ}, (C R_0)²)` and `⌈9 + 2/η⌉`.
  * `R_1(C)` depends on `C` and on the hypothesis constant `R_0`. It does not depend on the cluster, the
    template, the phase or the prime support.
  * The conclusion about `R < R_1(C)` is only `M(C) < ∞`. Nothing claims `M(C) ≤ K`.
  * `K` does not depend on `R`, `C`, the position or `N`'s factorisation.
* *Gaussian gcd and denominator costs.*
  * **Not discarded.** `g_S` is the full Gaussian gcd, including `(1+i)^{e_0}`, the inert part and the common
    split part.
  * `R''_S` is exact (`c2`: 0 failures, including `N = 103798401525000 = |γ|² · 1176852625` with
    `γ = (1+i)^3 · 3 · 7 · (2+i)^2`).
  * No residue, denominator or `t_ab` is used, so the cross-unit-class non-integrality of `t_ab`
    (kill_gauss minor 1) cannot arise.
  * `A_5` carries no unit-class restriction, and the 5-subsets are drawn from all unit classes (99.99% of tested
    subsets are cross-class).
* *Phase.* None occurs. `[z]` is similarity-invariant.
* *Fixed-`S` exceptional set.* None. The only exceptional set is the hypothesis `Z`.

### 3.2 The `(L_5^∘)` branch (major, for the text)

The README's status bullet ("So Question B is at least as strong as an absolute CG statement") and note item (3)
("any construction with `M_∞(C) → ∞` … would refute Question B") cover more than J4 proves.

`expert_note.md` §3 states Question B as "Prove `(L_5^∘)`, or `A_5(δ)` for some `δ < 1/30`", and `(L_5^∘)` is
weaker than `A_5` (gauss Prop 3.1; the converse is not proved). The only route from `(L_5^∘)` to a bound is
gauss Thm 3.4′ = `ptolemy.md` §6.2 Steps 0–5. Its Step 0 cuts the arc into `s = ⌈C/C_0⌉` pieces, with `C_0 < √2`
needed for `log|t| ≤ δ/2` in Step 1. The output is `M_orig ≤ 4s · (bound)`, so `K` grows linearly in `C`.

Making it `C`-free would mean redoing item 532's balance and residue estimates with an `O(log C)` term, absorbed
by `W ≥ log(R/C²) → ∞` for `R ≥ R_1(C)`. That is plausible but **not checked**: I did not re-derive item 532's
(2)–(4), (9) with `C` tracked.

**Correction.** In both places, replace "Question B" by "`A_5(δ)` for every `δ < 1/30` (hence `V_E(5)`)". If the
`(L_5^∘)` branch is wanted, it has to be proved separately.

### 3.3 Novelty: stated too narrowly (minor)

"What is new is only the `C`-independence of `K`." There is a second new point.
* The expert note's Reduction bullet and kill_gauss §7 item 1 say that `δ ≤ 2^{-k-1}` needs only the convexity
  averaging, and that larger `δ` needs item 532's balanced selection restricted to one unit class.
* J4's Markov step shows that the convexity averaging of the tuple gcd alone reaches the **sharp** range
  `δ < 1/(2^k-2)` for every `k`. The reason is that Markov is applied at level `1 - 2^{1-k} + η` rather than at
  `2^{-k}`.
* This also holds inside the original Thm 3.4 framework: with subdivision to `L ≤ √R_1` and Lemma 1.2,
  `θ = 2^{1-k} - η` gives `δ ≤ θ/(2(1-θ)) → 1/(2^k-2)`. So:
  * for the `A_k` branch, neither balanced selection nor the unit-class restriction is needed at any admissible
    `δ`;
  * they are needed only for `(L_k^∘)`, whose hypothesis is stated on balanced, one-unit-class tuples.

### 3.4 Strength calibration (minor)

* "Strictly stronger than `M(C) < ∞` for each `C`": the absolute form implies the per-`C` statement. No
  separation is known, so it is "formally stronger", not "strictly stronger".
* **The absolute form follows from Cilleruelo–Granville Conjecture 14** (CG 2007 p. 11; the record's copy is in
  `lit18/gh/…conjecture_13.md`): an arc of length `R^{1-ε}` holds a number of points bounded uniformly in `R`. With
  `ε = 1/4`, `C√R ≤ R^{3/4}` for `R ≥ C^4`, so `sup_C M_∞(C) ≤ B_{1/4}`. So `sup_C M_∞(C) < ∞` sits between the
  critical case and CG's own published conjecture. J4 shows that the `A_5` branch of Question B is at least as
  strong as this intermediate statement. It does not show that Question B is stronger than anything CG
  conjectured. This is what the "circle analogue of #887 `parts.ii`" means, and the note should say it.
* Oganesyan (arXiv 2107.09991, withdrawn per `lit18/primes.md` #26) concerns arcs `R^α` with `α > 1/2`.
  * Correct or not, it says nothing about `√R`-arcs (`lit18/direct.md` §2.10), so it neither supports nor
    threatens `sup_C M_∞(C) < ∞`.
  * In the note text it reads as if it were relevant evidence. Drop it, or say "(it concerns `R^α`, `α > 1/2`,
    and is withdrawn)".
  * A WebSearch this session still lists the preprint without showing the withdrawal. I rely on the record's
    check for the withdrawal and did not re-verify it.
* "The record's families give `M_∞(C) ≥ 8` for `C > 5.89771`" agrees with `overview.md` §1 and
  `referee_families.md`.
* "Contrapositive … would refute `A_5(δ)` for every `δ < 1/30`, and with it `V_E(5)`": correct.
  * By gauss Thm 4.2 (killer-checked), `V_E(5, ε)` with `ε < 2δ/c_A` gives `A_5(δ)` for every `δ > 0`.
  * The construction needs `sup_C M_∞(C) = ∞`. Since `M_∞` is nondecreasing in `C`, this is the same as
    `M_∞(C) → ∞`.

## 4. Exact checks (own code)

Directory: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round23/referee_J4_checks/`.
Python integers and `Fraction` throughout. Floats are used only for the reported means and for the evidence-only
A_5-scale count.

| script → output | content | result |
| --- | --- | --- |
| `c1_hypergeom.py` → `.out` | `q(M) = 1 - min_κ (C(κ,5)+C(M-κ,5))/C(M,5)` against `15/16 + 15/(16(M-4))`; the reverse-Markov step at `M = ⌈9+2/η⌉ .. +200` for six `η`; the exponent threshold | no violation for `9 ≤ M ≤ 3000` (smallest slack `1.05e-4` at `M = 3000`); `P ≥ η/2` in all cases; `(15/16)(1+2δ) = 1` exactly at `δ = 1/30`, `495/496` at `δ = 1/31` |
| `c2_identity.py` → `.out` | steps 1–2: for each 5-subset, the Euclid gcd `g_S`; `N/|g_S|² = Π p^{max a - min a}`; `g_S` an associate of the formula; common `v_{1+i}`, `v_q`; reduced tuple primitive, distinct, on `R''`, with the same class; `|g_S| ≤` min chord | **0 failures in 795,172 subsets.** These are the record's 10-point cluster (`n = 1176852625`, chord span `8.2386√R`, all four unit classes), the same cluster times `γ = (1+i)^3·3·7·(2+i)^2`, and all 5123 circles `N ≤ 20000` with `≥ 5` points (up to 300 sampled subsets each). Of the subsets, 393,328 are ramified, 104,000 inert, 474,001 non-normalised and 794,553 cross-unit-class |
| `c3_average.py` → `.out` | step 3 on actual clusters: per-layer probability `≤ q(M)` (exact) and `E log N''/log N ≤ q(M)`; Lemma 1.2 `R'' ≥ R/chord` | the record cluster has mean `0.9765 ≤ q(10) = 0.9921`; 120 clusters (`M = 9, 12, 16`, ten circles) all pass, with max mean/`q(M)` `= 0.9918`. No subset of these small clusters is yet at the `A_5(1/31)` scale (expected: `R ≪ R_1(C)`; evidence only) |
| `c4_schwartz_zippel.py` → `.out` | step 6: `Z` given by `Q`-equations in `(x_j, y_j)` and the homogeneous `F(u)`; equivalence on all distinct ordered tuples; zero count `≤ D M^4` | equivalence holds for all 3 `Z` × 3 circles (`M = 14, 16, 16`); every zero count is below `D M^4` |
| `c5_constants.py` → `.out` | explicit `η, κ, M_1, K` (exact bisection, verified on `[M_1-1, M_1+2000]`); the general-`k` threshold | for example, `δ = 1/40, D = 1`: `K = 32266`; `δ = 1/31, D = 1`: `K = 253450`; `M_1 ≤ 937.3 D/η` in all cases; `(1-2^{1-k})(1+2/(2^k-2)) = 1` exactly for `k = 3..8` |

**What the checks cannot test.** `A_5(δ)` itself, and any asymptotic statement. All tested radii are far below
`R_1(C)`. The checks confirm the exact identities and inequalities the proof uses, not the hypothesis.

## 5. Corrected text for expert_note §3 item (3)

> *Absolute form (round 23 judge J4, **refereed**).* If $A_5(\delta)$ holds for one $\delta<1/30$, with $(Z,R_0)$,
> then there is $K=K(\delta,Z)$ such that for every $C$, every arc of length $C\sqrt R$ on a circle of radius
> $R\ge R_1(C)=\max(C^{2/\kappa},(CR_0)^2)$ contains at most $K$ lattice points. Here
> $\kappa=1-(15/16+\eta)(1+2\delta)$ and $\eta=(1/(1+2\delta)-15/16)/2$. So $\sup_CM_\infty(C)<\infty$, while
> $M(C)<\infty$ for each $C$ with a $C$-dependent bound.
>
> The proof:
> * draws a random 5-subset of the arc's points;
> * divides it by its full Gaussian gcd, which gives the exact reduced radius
>   $R''^2=\prod_pp^{\max_Sa(p)-\min_Sa(p)}$ (ramified and inert parts cancel), with
>   $\mathbb E\log R''^2\le(15/16+O(1/M))\log R^2$;
> * applies $A_5$ to the good subsets and one equation of $Z$ through Schwartz–Zippel.
>
> No unit class, no balanced selection and no residue is used. The same argument shows that $A_k(\delta)$ for one
> $\delta<1/(2^k-2)$ gives the same absolute conclusion for every $k\ge3$. In particular, the convexity averaging
> of the tuple gcd alone already reaches the sharp range in the "Reduction" bullet.
>
> Consequently, a construction with $M_\infty(C)\to\infty$ would refute $A_5(\delta)$ for every $\delta<1/30$, and
> with it $V_E(5)$. Whether it would also refute $(L_5^\circ)$ is **not** shown: the $(L_5^\circ)$ route
> (Thm 3.4′) subdivides the arc and gives a $K$ linear in $C$.
>
> For calibration, $\sup_CM_\infty(C)<\infty$ also follows from Cilleruelo–Granville's Conjecture 14 (arcs of
> length $R^{1-\varepsilon}$ hold boundedly many points; take $\varepsilon=1/4$). So it lies between the critical
> conjecture and CG's $R^{1-\varepsilon}$ conjecture, and is the circle analogue of Erdős #887 `parts.ii`. The record's
> families give $M_\infty(C)\ge8$ for $C>5.89771$, and nothing unbounded.

Also amend:
* the "Reduction" bullet: "For $\delta\le2^{-k-1}$, only a convexity averaging…" becomes "For every
  $\delta<1/(2^k-2)$, the convexity averaging of the tuple gcd suffices for the $A_k$ branch (round 23 J4);
  balanced selection in one unit class is needed only for $(L_k^\circ)$";
* README §4, "strictly stronger" becomes "formally stronger";
* README §5 J4 status, "Question B is at least as strong as…" becomes "the $A_5$ branch of Question B is at least
  as strong as…".

## 6. Files

* Report: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round23/referee_J4.md`
* Checks: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round23/referee_J4_checks/`
  * `gz.py` (Gaussian helpers)
  * `c1_hypergeom.py`, `c2_identity.py`, `c3_average.py`, `c4_schwartz_zippel.py`, `c5_constants.py`, each with
    its `.out`
