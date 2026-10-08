# Exact finite optima for simultaneous coset-weight constraints

This is a finite rational linear-program calculation, not an endpoint
configuration or an asymptotic growth theorem. It tests all coset
characters at once, with one unit-weight physical flip per row and
class weights `W_a>=5`. The remaining weight in each class may be
distributed among its unflipped physical columns.

For each nonzero row subspace V, let `h=|V|`, choose a nonzero
restriction ell, and maximize the matching count k over every row
coset. With `W=sum_a W_a`, the condition that each such character
has height at least `W/2` is exactly

```text
h sum_(a|V=ell) W_a-W >= 4k-2h.                           (1)
```

The full-support constraints `h=M` are included. Equation (1) is a
formal half-height condition. It omits the actual arc phase allowance,
prime-log arithmetic, and common content.

The [checker](check_walsh_simultaneous_coset_weight_lp_20260915.py)
generates every subspace once in binary reduced row-echelon form and
every nonzero restriction. The
[rational fixture](walsh_simultaneous_coset_weight_lp_20260915_fixture.json)
supplies a primal and a nonnegative dual with equal objective values.
Default verification requires only the Python standard library; SciPy
was used to locate the certificates, and remains optional with `--solve`.

| Assignment | Variables | Constraints | Exact minimum W | Best uniform W |
| --- | --- | --- | --- | --- |
| Archived nonlinear 32-row map | 31 | 2,077 | 155 | 155 |
| 64-row shuffle with seed 20260915 | 63 | 23,562 | 496 | 504 |

The archived map is the one in the
[full-character classification](walsh_all_character_half_height_classification.md).
The 64-row map shuffles the labels `1,...,63,1` with the stated seed,
so every label is used at most twice. At 32 rows, all class weights are
five and all constraint slacks are at least one; the lower-bound dual
alone proves optimality. At 64 rows, the exact dual uses 58 coset
constraints and no class-lower-bound constraints. The checker verifies
every primal inequality, nonnegativity, dual stationarity, complementarity,
and equality of primal and dual objectives with rational arithmetic.
Root separately replayed both solves and the stored certificates.

The 64-row primal has a simple Walsh expression. Write
`chi_d(a)=(-1)^(a dot d)`. For every nonzero label,

```text
W_a=63/8+(chi_4(a)+chi_29(a))/16.                         (2)
```

Its values are `31/4`, `63/8`, and `8`, with multiplicities 16, 32,
and 15. The auxiliary value of the right side of (2) at the zero label
is eight.
The two exceptional transform directions are 4 and 29. Subtracting
the flip correction `4n_a/M`, where `n_a` is the assignment count,
gives nonzero-direction transform values

```text
qhat(4)=-4,       qhat(29)=-31/8,
qhat(d)=-8       for other even d,
qhat(d)=-63/8    for other odd d.
```

Thus a small, structured redistribution can change the optimum even
after imposing every coset constraint. The finite improvement from
504 to 496 concerns a diagnostic LP, not the requested growth rate.
It supplies no obstruction for unbounded M and no actual lattice
points on short arcs.
