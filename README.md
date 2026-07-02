# Jarnik

Exact integer tools for experimenting with endpoint-scale lattice-point arcs on
circles, following the Gaussian-chain certificate handoff in this repository's
Codex attachment.

The current working target is the rooted 4-point certificate obstruction:

```text
#(I cap Z^2) <= 4 for arcs |I| = R^(1/2)
```

The code keeps the old floating-point arc scanner, but the main new machinery is
integer-only certificate generation and fifth-point extension checking.

## Source Layout

- `GaussianChain/` contains the Lean formalization.
- `src/numerics/` contains the Rust numerical code: exact arithmetic,
  Gaussian-integer factorization, lattice-point generation, cluster diagnostics,
  certificate search, and extension checking.
- `src/numerics/legacy_arc.rs` contains the original floating-point arc scanner
  used by the `arc` command and Criterion benchmarks.
- `outputs/` is reserved for generated numerical JSON/JSONL artifacts.

## Commands

Run commands with `cargo run -- <command>`.

```text
search --n-max N [--arc-scale C] [--min-cluster M] [--out path.jsonl]
```

Scans `n <= N`, generates lattice points from Gaussian-prime factor choices,
and emits JSONL records for circles whose largest endpoint-scale arc cluster has
at least `M` points.

```text
arc-families --n-max N [--n-start A] [--arc-scale C] [--min-size M] [--max-records K] [--out path.jsonl]
```

Fast discovery scan for endpoint-scale arc families, intended for 5- and
6-point searches. It builds a smallest-prime-factor sieve up to `N`, skips
circles that cannot have enough lattice representations, computes only the
maximum cluster size during the scan, and emits compact JSONL family records
instead of the full diagnostic tables used by `search`.

```text
analyze-n --n N [--arc-scale C] [--out path.json]
```

Analyzes one circle, including all maximum clusters, full relation-height
diagnostics, exponent-layer cut diagnostics, transition-shape summaries, and
half-angle divisor data. The primary statistic is
`full_relation_height.budgets[*].min_h_over_log_r`.

```text
cg-family [--n-start A] [--n-end B] [--arc-scale C] [--out path.jsonl]
```

Generates the Cilleruelo-Granville Fibonacci 4-point family. Family index `7`
is the known example `n = 567454025`; analysis is emitted for generated circles
whose radius square fits the current `u64` pipeline.

```text
random-squarefree --num-primes K --trials T [--arc-scale C] [--seed S] [--prime-limit P] [--out path.jsonl]
```

Samples deterministic squarefree products of primes `p == 1 (mod 4)` and
analyzes the maximum endpoint clusters for each trial.

```text
known-cert
```

Prints the first known endpoint-scale rooted 4-point certificate:

```text
(x,y,z) = (233,1098,1597)
(U,V,w,c,J,h,N) = (6,6,10,2,72,5,22698161)
```

```text
generate <max_z> <max_det> [limit]
```

Searches over `(x,y,z,U,V)` and emits certificates as JSONL. The generator uses
the congruence

```text
V*x == -z*U (mod y)
```

before constructing large certificate expressions. Its stderr counters report
`checked_blocks` for `(z,y,U,V)` blocks and `congruence_hits` for candidate
`x` residues.

Example:

```bash
cargo run -- generate 1597 6
```

rediscovers the known certificate.

```text
extend-known [max_w] [comma_primes]
```

Runs the extension quadratic/discriminant checker for the known certificate. If
`max_w` is omitted, the code finds the finite positive discriminant window. The
default CRT sieve uses:

```text
3,5,7,11,13,17
```

Example result:

```text
upper_w=6020700
witnesses=0
```

```text
verify-circle <radius_squared>
```

Generates all integer points on `x^2 + y^2 = radius_squared`, sorts them around
the circle, and reports the largest endpoint-scale cluster. The final integer
sanity check uses the exact chord condition `chord^4 <= radius_squared`; a
floating arc-width count is printed as a diagnostic.

Example:

```bash
cargo run -- verify-circle 567454025
```

recovers the four lattice points from the handoff.

```text
sign-search <a> <n> <max_r> <max_t> [min_vectors] [max_results]
```

Searches weighted sign-product coefficient collisions for

```text
prod_j (a + r_j + i*sigma_j*t_j)
```

using the first three endpoint-scale angle coefficients.

```text
arc [min_radius_squared] [max_radius_squared]
```

Runs the original floating-point Jarnik/Cilleruelo-Granville style arc-ratio
scanner.

## Implemented Handoff Pieces

- Fast rooted certificate generation over `(x,y,z,U,V)`.
- Exact certificate reconstruction checks, including the `n2` and `n3` square
  identities for geometric certificates.
- Exact fifth-point extension quadratic checker.
- Modular quadratic-residue CRT sieve for extension discriminants.
- Direct lattice-circle sanity verification for reported examples.
- Optimized endpoint-scale arc-family scan for 5- and 6-point numerical
  searches.
- Weighted sign-product collision search scaffolding.
- Lean 4 + Mathlib project with the Gaussian-chain algebra and the log/loglog
  proof split across the `GaussianChain/` modules.

## Lean Status

Lean is installed through `elan`:

```text
Lean 4.30.0
Lake 5.0.0-src+d024af0
```

The Lean project builds with:

```bash
lake build
```

Formalized pieces include the rooted coordinate identities, determinant
identity, certificate structures, extension quadratic, discriminant-square
implication, known example certificate, the geometric leading-discriminant
simplification, and the log/loglog arc-count bounds in
`GaussianChain/MainSublogBound.lean` and `GaussianChain/AsymptoticParameters.lean`.

The endpoint no-fifth-point theorem is separate from the log/loglog result and
is not needed for that proof.
