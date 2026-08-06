# Jarnik

Exact integer tools and Lean formalization for lattice points on short arcs of circles.

## Proven core

The main development proves the current sublogarithmic arc bound and contains certificate-search machinery. The umbrella module `GaussianChain.lean` imports only this proved core.

## Uniform endpoint problem

The target is a uniform bound for lattice points on arcs of length

```text
C * R^(1/2).
```

This remains open in this repository.

The earlier Ru--Vojta/generalized-GCD and determinant equality-case programs are no longer active. They are retained as an audit trail, but they are not imported by the umbrella build and should not be read as the current proof plan.

The inverse-concentration reset begins at:

- `docs/research_reset_inverse_concentration.md`

The current ordered/residual frontier is:

- `docs/short_gap_prime_reveal_renormalization.md`
- `docs/positive_subset_transition_obstruction.md`
- `docs/inert_prime_residual_growth.md`
- `docs/adaptive_residual_diagonal_countermodel.md`
- `docs/transition_content_cocycle.md`
- `docs/interval_content_aggregation_audit.md`
- `docs/high_return_pattern_descent.md`
- `docs/five_product_capacity_audit.md`

The inert-prime note proves aggregate `m^2 log m` growth of the forced-reduced
Plucker residuals in the odd squarefree split model and sharpens the resulting
sublogarithmic coefficient.  The transition-content note identifies the exact
phase/conductor cocycle: prime re-entry is the rational content removed from a
product of primitive transition blocks.  The all-interval audit then proves
that scalar content cancels exactly back to the old pairwise cut bound, so
content size alone cannot close the endpoint.

The high-return descent note is the latest positive reduction.  It shows that
a hypothetical `n`-point counterexample must contain linearly many disjoint
complete orientation-pattern blocks, each of norm at most
`N^(2/n+o(1))`, and linearly many short transitions each involving linearly
many such blocks.  Thus the remaining uniformity problem is a dense Gaussian
phase-cancellation problem.  The countermodels show that support, collision,
run-count, and bounded-depth non-reuse estimates alone cannot supply the
missing contradiction; the next input must quantify critical lcm-height
efficiency of the actual signed Gaussian products.

The five-product audit tests the first proposed fixed-dimensional resultant at
the exact finite valuation level.  It proves that pair determinants, Plucker
cancellation, rank four, and pattern-forced monomial content cannot provide the
required fourth-order decay with a one-use conductor budget.  It replaces that
failed universal lemma by an exact weighted occurrence-pattern capacity test.

The broader program starts from the global exponent-angle model. After factoring the common norm, every point is represented by an exponent vector `a` in a Gaussian conductor box, and its angle is a linear form

```text
Phi(a) = sum_j (2 a_j - e_j) theta_j  mod 2 pi.
```

A hypothetical counterexample gives arbitrarily large sets of exponent vectors whose images lie in an interval of width `O(N^(-1/4))`.

The new program is to combine:

1. a weighted inverse-concentration theorem for these exponent boxes; and
2. arithmetic rigidity of the Gaussian prime angles
   `exp(2 i theta_j) = pi_j / conjugate(pi_j)`.

The intended contradiction is that high concentration forces either a heavy prime block, which can be descended, or a bounded-rank approximate relation that becomes an impossible exact multiplicative relation by Gaussian unique factorization.

## Historical research files

Files concerning finite beta calculations, uniform exceptional sets, adapted sextic bases, Segre embeddings, determinant surplus, and balanced prime-layer cuts are historical diagnostics only. They record failed or conditional routes and are intentionally excluded from `GaussianChain.lean`.

## Source layout

- `GaussianChain/` contains the Lean formalization.
- `src/numerics/` contains the Rust numerical code: exact arithmetic, Gaussian-integer factorization, lattice-point generation, cluster diagnostics, certificate search, and extension checking.
- `outputs/` is reserved for generated numerical JSON/JSONL artifacts.

## Commands

Run numerical commands with `cargo run -- <command>`.

```text
search --n-max N [--arc-scale C] [--min-cluster M] [--out path.jsonl]
arc-families --n-max N [--n-start A] [--arc-scale C] [--min-size M]
analyze-n --n N [--arc-scale C] [--out path.json]
cg-family [--n-start A] [--n-end B] [--arc-scale C] [--out path.jsonl]
random-squarefree --num-primes K --trials T [--arc-scale C]
known-cert
generate <max_z> <max_det> [limit]
extend-known [max_w] [comma_primes]
verify-circle <radius_squared>
sign-search <a> <n> <max_r> <max_t> [min_vectors] [max_results]
arc [min_radius_squared] [max_radius_squared]
```

## Lean

Lean is installed through `elan` and the proved core builds with:

```bash
lake build
```
