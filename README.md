# Jarnik

Exact integer tools and Lean formalization for endpoint-scale lattice-point arcs
on circles.

The current `master` development proves the sublogarithmic bound and contains
certificate-search machinery.  The branch
`agent/uniform-gcd-formalization` adds a separate research track for the
uniform-bound problem.

## Uniform generalized-GCD research track

The branch contains:

- `GaussianChain/FiniteBeta.lean`: an explicit finite Ru--Vojta beta-level
  calculation on the blowup of `P^2` at `[1:1:1]`;
- `GaussianChain/UniformExceptionalReduction.lean`: the conditional
  `1/2 - eta` endgame and grid/fiber cardinality bound;
- `GaussianChain/SlopeClosure.lean`: maximal-slope closure under meet and join
  for modular rank and supermodular degree;
- `GaussianChain/FiniteFlatReduction.lean`: finite-signature range bounds for
  canonical exceptional objects;
- `docs/uniform_gcd_route.md`: the full circle-side reduction and current gap
  audit;
- `docs/canonical_hn_uniformity.md`: a self-contained derivation showing why
  the canonical Harder--Narasimhan subspace ranges over a fixed finite list,
  independent of the number of places.

The branch does not add an axiom for the remaining Diophantine input.  The
remaining task is to match the finite Ru--Vojta max-over-adapted-bases system
to a uniform large-parameter form of the absolute parametric Subspace Theorem.

## Source Layout

- `GaussianChain/` contains the Lean formalization.
- `src/numerics/` contains the Rust numerical code: exact arithmetic,
  Gaussian-integer factorization, lattice-point generation, cluster diagnostics,
  certificate search, and extension checking.
- `src/numerics/legacy_arc.rs` contains the original floating-point arc scanner.
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

Lean is installed through `elan` and the project builds with:

```bash
lake build
```

The GitHub Actions workflow on the uniform-GCD branch runs the full Lean build.
