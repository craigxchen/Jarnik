# Occupancy-weighted coefficient floor

This note records the occupancy count in the proved tensor multiplicity floor
of [`higher_rank_tensor_multiplicity_floor.md`](higher_rank_tensor_multiplicity_floor.md),
which strengthens the Weyl-union floor in
[`higher_rank_weyl_union_floor.md`](higher_rank_weyl_union_floor.md). The exact
arithmetic statement is conditional on the core factorization, coprimality,
and conformal-block annihilator hypotheses in that theorem.

Let `q=2^d`. For a subset `T` of the `nd` labels, let `f(T)` be the number
of coordinate blocks wholly contained in `T`, and let `e(T)` be the number
of coordinate blocks disjoint from `T`. The occupancy count is

```text
M(n,d) = sum_T min(f(T), e(T))
       = sum_(f,e>=0, f+e<=n) n!/(f!e!(n-f-e)!)
           * (q-2)^(n-f-e) * min(f,e).
```

The existing Weyl-union count is

```text
U(n,d) = sum_T 1_(f(T)>=1 and e(T)>=1)
       = q^n - 2(q-1)^n + (q-2)^n.
```

Thus `M-U` counts the extra multiplicity when a subset contains at least two
full blocks and at least two empty blocks. In particular, `M=U` for `n<=3`.
For fixed `n` and `q -> infinity`,

```text
U = n(n-1)q^(n-2) - 6*C(n,3)q^(n-3) + O_n(q^(n-4)),
M-U = C(n,2)C(n-2,2)q^(n-4) + O_n(q^(n-5))  (n>=4).
```

For fixed `q` and `n -> infinity`, view each block independently as full,
empty, or partial with probabilities `1/q`, `1/q`, and `1-2/q`. If `F,E`
are the full and empty counts, then

```text
M = q^n E[min(F,E)]
  = q^n (n/q - sqrt(n/(pi*q)) + o(sqrt(n))).
```

For the last asymptotic, write `min(F,E)=(F+E-|F-E|)/2`.
The centered difference is a sum of independent bounded variables of
variance `2/q`. The central limit theorem, together with the uniform
second-moment bound for `(F-E)/sqrt(n)`, gives convergence of the
absolute first moments and hence the displayed coefficient.

The exact checker compares `M` with the existing determinant threshold
`B(n,d)/h(n,d)`. On the grid `2<=n<=8`, `2<=d<=12`, it finds no strict
gain over that threshold: equality holds for `d=2`, and `B/h>M` for every
`d>=3`. Extra samples through `n=12` for `d=3,...,6` have the same direction.
These finite comparisons do not establish an all-parameter inequality; the
`d=2` equality for every rank is proved in the linked tensor-floor note.

Run `python3 docs/check_occupancy_candidate_floor.py` to reproduce the
finite comparisons and verify the occupancy formula against direct subset
enumeration on small instances.
