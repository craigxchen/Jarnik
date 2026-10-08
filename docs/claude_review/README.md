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

Open. Growth-rate work from this review is recorded separately as it is refereed.
