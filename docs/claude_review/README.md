# Independent review of the uniform endpoint program (2026-10-08)

This directory records an independent review of the research notes on branch
`agent/uniform-gcd-formalization` (pushed 2026-10-08, HEAD `5b47737`), requested in
`docs/claude_handoff_uniform_endpoint.md` on that branch. The note names below refer to files in
`docs/` on that branch. Nothing here proves the uniform theorem; it remains open.

## 1. Audit verdicts

Each cluster was audited by an independent agent that read the notes in full, reconstructed
every numbered claim, reran the cited checkers, and wrote its own exact tests. A second agent
then tried to overturn every non-trivial finding. No mathematical error was found that affects
a stated theorem.

| Cluster | Notes | Verdict |
| --- | --- | --- |
| General bound | `inert_prime_cofactor_bound.md` | Verified step by step: normalisation, conic count `(p+1)p^(a-1)`, `E(M,q)`, finite inequality (6), asymptotics (10). Exact tests on 267,102 actual configurations (3.3M pairs, 261,119 with mixed unit classes). The "1,035 configurations / 7,859 checks" sentence names no checker and could not be audited; the new sweep replaces it. |
| Squarefree 2/3 | `inert_prime_residual_growth.md` | Verified in its odd squarefree model. **Correction:** the §5 remark that prime powers dilute the split gain to `m^2/(2(e+1)(p-1))` is true for the residual product alone but false as an obstruction to the final coefficient (see §2 below). Upheld on re-check. |
| Item 532, selection | `exact_nested_profile_residual_extraction.md` | Verified: simultaneous tuple selection (`E A <= 2(2^k-1)/(M-k+1)`, `E B <= k(k-1)W/[4(M-1)]`, one tuple with both), pattern error `eta`, aggregate bound `H_Sigma <= k(k-1)W/[4(M-1)]`, threshold `Q`, and the conditional implication (12) including both growth readings. The obtuse Bessel step matches `GaussianChain/ObtuseBessel.lean`. Cosmetic: only `C_0 < 2` is needed; the handoff's one-line summary of (12) omits the fixed neighbourhood width `epsilon_0` and threshold `w_0` that the note states. |
| Item 532, gcd identities | same note | Verified as exact Gaussian identities (not only norms): unit-class normalisation, (5)-(8) with two nested chains per prime, `K_i = 1`, `n_empty`, anchor invariance, (9)-(10), Hankel divisor (13), and the Cauchy-Binet step the checker does not test. Tested on genuine lattice points up to `N ~ 10^127`, exponents up to 14. The checker's coverage is narrower than the note suggests (one prime exhaustively, one multi-prime profile). |
| Reformulations | `integer_cotangent_lcm_height_target.md`, `ordered_pair_triple_gap_denominator.md`, `ordered_conductor_index_uniformity_criterion.md`, `endpoint_full_profile_quantifiers.md` | All 23 claims verified. Sharpening: (H) with `k` coordinates holds iff some arc length `c_1 sqrt(R)` holds at most `k` points, so the least admissible `k` is `lim_{C->0+} M*(C) >= 5`. The all-edge lcm matches the Gaussian primitive radius on 19,064 actual sub-tuples; the anchor-only lcm is wrong in 4,362 of them. The five-point family has `K = 1` exactly. |
| Arithmetic lemmas | `thin_edge_divisor_matching_criterion.md`, `four_row_one_coordinate_divisor_reduction.md` §5-7, `positive_hankel_full_profile_budget.md` | Verified (143,120 exact thin-division pairs, 600 random 15-block profiles, 148 Hankel profiles). **Correction:** in §7 of the four-row note, "Equivalently, gcd(\|Y_i\|,(\|t_ij\|)_j) = gcd_i \|Y_i\|" should read "Consequently"; explicit data show it is strictly weaker than (17). |
| One-class cyclotomic count | `resonant_fibonacci_cyclotomic_reduction.md`, `cyclotomic_affine_rank_reduction.md`, `cyclotomic_uniform_rank.md`, `cyclotomic_nullity_four_counterexample.md`, `cyclotomic_five_row_subendpoint_family.md` | Verified (`rho <= 1125`, `2251` rows, `1126` when `L < 4h`); nullity-four and five-row counterexamples confirmed exactly. **Corrections:** the §5 sentence of the nullity-four note calling non-squarefree compression open is stale (`cyclotomic_uniform_rank.md` §1-2 proves it); the five-row family works for every odd `k >= 19`, not only `k = 0,1 mod 3`. |
| Packet projections | `mixed_affine_cyclotomic_rank_scope.md`, `common_slope_...`, `frequency_separated_...`, `divisibility_chain_...`, `three_class_affine_uniform_count.md`, `multiprime_packet_cost_obstruction.md` | Verified, including preserved low coefficients (with the factors `A_c S_(c,k) gamma_c^k`), the kernel of the complete composition inside the direct sum of classwise truncated Mobius kernels with 1125 charged once, termination, and the two-term phase separation. **Corrections:** the §4 sign remark in the mixed note depends on the convention (conclusion unaffected); the three-class note should restate that its entry index depends on the template. |

The reports themselves were working files and are summarised here; the checkers they rely on
are the notes' own checkers, all of which pass.

## 2. A general (2/3 + eps) log R / loglog R bound (constant only)

`two_thirds/proof.md` proves, for every arc of length `<= C sqrt(R)` anywhere on
`x^2 + y^2 = R^2`,

```text
3 M log M <= log(R^2) + (14 + 4 log^+ C) M,
```

hence `M <= (2/3 + eps) log R / loglog R` for `R >= R_0(C, eps)`. It extends the exact finite
inequality (6) of `inert_prime_cofactor_bound.md` to split primes, both those not dividing `N`
and those dividing it at equal exponent levels. The new ingredient is a **chain lemma**: at a
split prime `p | N` with exponent `e`, the cut slack and the level collisions together satisfy

```text
sum_tau (S_tau - M/2)^2 + (1/(p-1)) sum_a n_a^2 >= c(1/(p-1)) M^2 / 2,
c(lambda) = (sqrt(1+4 lambda) - 1)/2 = lambda - lambda^2 + O(lambda^3),
```

so spreading points over many exponent levels costs cut slack that more than repays the lost
collisions. Two adversarial referees found no error. **This changes only the leading constant;
it is recorded here and not pursued further.**

## 3. Status of the uniform theorem

**For outside experts:** `expert_note/expert_note.md` is a short (about 5 pages) self-contained statement of the three open questions the record reduces to: `(L_8)` (integer `2 x 8` matrices; implies the uniform theorem; Vojta-predicted), `(L_5^circ)` (five lattice points on a circle; implies the uniform theorem; the 4-point version is false), and positive-gain non-unit pencils on `Mbar_(0,8)`, with what any proof must avoid. Accuracy-checked (`expert_note/check.md`).

**Start here:** `overview/overview.md` is a consolidated account of the whole record for outside experts (problem, the
`(2/3+eps)` theorem with a proof sketch, restricted and conditional theorems, the barrier map, routes and where each stops,
open problems), with a status label and file pointer on every statement. It was checked claim by claim against the record
by an adversarial checker (`overview/check.md`, 4 major and 23 minor findings, all applied).

Open. The growth-rate and uniform-bound work of this review is indexed in `growth/README.md`:
two conditional theorems (Vojta's conjecture for one fixed variety implies the uniform bound),
restricted uniform and growth theorems for Hadamard and Walsh profiles, the cotangent trade-off,
new lower-bound constructions and exhaustive data, and the missing lemmas that would complete
each route.

The auxiliary-polynomial round is in `threshold/README.md`: no form on any equal-norm variety
`X_M` vanishes along the diagonal to order `2 deg F` (`tau*(M) = 2 - 1/ceil(M/2)`), so neither a
single-place Liouville argument nor Ru-Vojta with `S = {infinity}` can prove the theorem, and
forced Gaussian divisibility only reaches criticality.

Round 4 (flipped Hadamard cores, (L_7), hardness, fake-circle barrier) is in `round4/README.md`.

Round 5/6 (fake circles at `log R = O(M log M)` with residue data) is in `round5/README.md`. Round 7
(residue data factoring through the primes, every monomial strengthened: the last local class not
known to be capped) is in `round7/README.md`; that class is still open, with
`3 M log M - O(M) <= W_min <= O(M log^2 M/loglog M)`. Round 8 (growth rate: Lang-Waldschmidt,
simultaneous hypotheses, stability of extremal clusters) is in `round8/README.md`. Round 9 (the
combinatorial missing lemmas L4, L5 and L1/L2 of round 7) is in `round9/README.md`: all partial, no
growth improvement and no new cap. Round 10 (`round10/README.md`) refutes the inverse step L4(K) as
formulated, and finds the most promising new route: in the growth regime the determinant method with
forced divisibility becomes supercritical, so every large actual cluster satisfies an exact
algebraic system that almost every fake violates. Under an unproved zero estimate (ZE) this would
give `M = O(loglog R)` on squarefree norms. Round 11 (`round11/README.md`) attacks ZE. ZE is still open:
the only evidence for it is exact theorems on two small LP-subcritical models. Every stronger or angle-free
replacement proposed (ZE-braid, the separation hypotheses SL and SL\*) is false, through the cheap rational
subtori of Hadamard-type profiles. New: ZE restricted to the actual point (ZE-col\*) would give
`M = O((loglog R)^2)` on all norms, and route 5 also runs at box size `T = 1`, where ZE reduces to a real
tangent-cone statement (unrefereed). Round 12 (`round12/README.md`) makes route 5 at `T = 1` rigorous on all norms,
with no rank or squarefree input: rigidity, the reduction to the tangent cone, the structure of cheap families and the
separation of the actual point from them are proved and refereed. Growth then follows from one open inequality,
MIN^pair: the tangent-cone forms cannot all be tiny on the set of directions that respect every pair gap of the
cluster. MIN^pair would give `M = O(loglog R)` on all norms. It is open, with no evidence either way at regime sizes,
which are beyond exact computation (the smallest LP-supercritical Sylvester model already has `n >= 32`). By the
fake-test theorem, deciding it on the fake profiles is exactly the barrier question. Round 13 (`round13/README.md`)
tests whether the slab-box's leading-form ideal is local (generated by projections onto boundedly many directions).
Locality is neither proved nor refuted in the regime: locality onto columns is refuted exactly in a Paley model, and
even a proof of locality would not settle the barrier (a fake model satisfies it and still has nondegenerate base
points). The judge recommends pausing route 5 and turning to the construction side of the residue class P3-all\*.
Round 14 (`round14/README.md`) tests that class in the growth regime. Rigidity of all-norm-good label systems
(Rig(K)) is undecided: its rigid core grows with K for the tightest primes, but the round's flexibility witnesses
all fail the single-angle coupling of P3-all\*, so they are not fakes. No fake (C1) was built; several
construction routes are proved to fail. P3-all\* remains open and uncapped there.

Round 15 (`round15/README.md`) was interrupted by a usage limit; its theory part shows that no outcome of the
P3-all\* rigidity question changes the growth rate of actual circles by itself. `local_smoothing/README.md`
assesses OpenAI's claimed proof of the local smoothing conjecture in all dimensions (Oct 2026): not useful for
this problem (nothing new in the plane; vacuous on a single circle; L^p and energy chains are capped at
Theta(log R); fake circles satisfy every such estimate). Round 16 (`round16/README.md`) recasts route 5's last
inequality as a statement about the slab box alone: MIN^pair is an L^2 restriction statement about the top component of
<k,u>^m', and the polynomial anti-concentration literature is provably vacuous at the degree it needs. The barrier side
reduces to one quantitative test on the Fourier circles of the skew double core boxes. Route 5 is open, stuck, not capped. `lit18/README.md` is a literature sweep (2015-2026: the direct
problem, height-based tools uniform in `S`, genuine-prime statistics), each part checked by an adversarial referee. Nothing found
beats the barrier: no paper proves `o(log R)` for every `R` (the best all-`R` literature bound is Bourgain-Rudnick's `O(log R)`).
Its barrier-map output: the round-8 Subspace route `QS` is dead for every threshold that depends on the prime set only through its
size (Props A, B), composite-supported analytic inputs remain unclassified, and the one surviving hypothesis `QS^prod` is at
least as strong as the target. Round 17 (`round17/README.md`) runs the Fourier-circle test. H_light is proved for
`sdc_q` (`c_0 = 1/6` for all `q >= 7`, `>= 0.2349` in the regime), so the circle cap test needs no hypothesis, and on the circles
the quantitative test reduces to an exact tangent-cone membership against one witness height. The test itself is undecided at
every regime `q`, so **route 5 is stopped**: neither capped nor unblocked. The judge's remaining target is the arithmetic lemma
`(L_8)` (Vojta on `Mbar_(0,8)`), with an estimated chance of proof well below 1%. Round 18 (`round18/README.md`) examines `(L_8)` by census, lifting test and theory. It stays
**open**, with no concrete exceptional set. Proved: a balanced 8-datum needs at least 24 primes and has all minors `>= 10^49`,
so no feasible computation can test `(L_8)`; the `(L_6)` curves do not lift through class `beta_7`; on pairwise-coprime
profiles the chain, anchoring and residue structure carry no information. One uncovered conditional mechanism (abc on
auxiliary pencils of positive gain) is identified; no such pencil exists in degree `<= 4`. The judge recommends **stopping
internal rounds** on growth and uniformity and putting `(L_8)` before outside experts. Round 19 (`round19/README.md`) decides that one item. A
beta-weighted log Riemann-Hurwitz inequality on `Mbar_(0,8)` proves that **no unit pencil has positive gain**, in any degree
(refereed), and an exact search confirms gain `<= -32` for all unit pencils of degree `<= 6` (refereed). The abc lane is closed;
only non-unit pencils remain, with no proved conditional implication. Internal rounds on growth and uniformity stop here; the
consolidated `overview/overview.md` is updated through round 19. Round 20 (`round20/README.md`) is a second fresh-mechanism panel
(circle-specific Gaussian arithmetic at small `k`, genus-0 arithmetic intersection theory, composite-supported analytic inputs,
unlikely intersections in tori), each attacked by an adversarial killer. All four lenses were killed; no live lead survives. New
on record: robust composite-supported analytic inputs are capped; genus-0 Arakelov inequalities reduce to pair data, `tau*` and
rank; unlikely-intersection inputs are vacuous or capped at `omega(N')`; and a fifth conditional form, a Vojta instance `V_E(5)`
on a 4-fold whose 5-point statement `(L_5^circ)` implies the uniform theorem (the 4-point analogue is false), every known method
being exactly Liouville-critical on it. Round 21 (`round21/README.md`) attacks the one regime with no known fakes, the low-height window
`P_+(N') <= M^(2-delta)`. It stays **open**: no restricted growth theorem inside it and no cap. Proved: restricted growth just below
it (`P_+(N') = o(log R)` implies `M = o(log R/loglog R)`); the large sieve there reduces to arc counts and to the cross-gap item;
the record's Haar-completion fake constructions provably cannot run there (success probability `<= exp(-c M^(delta/2))`); angle
models are genuine on every complete initial segment of split primes. The judge again recommends stopping; the overview is
updated through round 21. Round 22 (`round22/README.md`) restates `(L_5^circ)` in the approximation-constant literature's terms: it is
`alpha^ess_e(-K) > 15/4` at the identity of the twisted weak-Fano Losev-Manin 4-fold, implied by McKinnon-Satriano canonical
boundedness at `e`, which is only a hypothesis (or a consequence of Vojta) in the literature. Every Liouville/section/Seshadri
argument is proved to stop at `tau_e <= 15/4`; no known theorem reaches the target. Round 23 (`round23/README.md`) proves `tau_e(-K) = 10/3`
exactly (so toric Liouville gives exactly Cilleruelo-Cordoba), and shows the split twin of Question B implies the open Erdős
problem #887 (absolute form): it is not easier and supplies no mechanism; a 2026 five-divisor construction is a balanced rational
curve at the Liouville constant that cannot exist on the twisted form.

**Growth-rate status (2026-10-28, after round 23).** No growth improvement over `log R/loglog R` is proved, either
unconditionally or from a standard conjecture. The following are capped at `log R/loglog R` by
refereed barrier theorems:
* every local input: pair chords, gcds, cut, residues at the point level, and prime-level residues
  used through characters. These residue caps hold for residue moduli up to `M^A` with `A` fixed in
  advance; round 9 notes that a growth argument may choose `A` after the cluster's size parameter,
  and caps uniform in `A` are not proved;
* character, GAP-Fourier and energy arguments;
* single-place auxiliary polynomials;
* Baker/Matveev bounds;
* single-form Lang-Waldschmidt.

The only conditional route known to work from a standard conjecture is Vojta's conjecture for one fixed
variety, which gives the uniform bound. Route 5 (rounds 10-12) gives `M = O(loglog R)` on all norms from the
non-standard inequality MIN^pair (round 12), and earlier from ZE or ZE-col\*. All of these are open. Still open: the class P3-all (genuine prime-level residues imposed on every monomial),
and the full Lang-Waldschmidt conjecture with `eps -> 0` used simultaneously.
