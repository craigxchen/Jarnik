# K5 variable-ratio identity over Q(t)

Fix the adjacent edge values `x_01=1` and `x_02=t`. Let
`(a,b,c,d)=(x_14,x_23,x_24,x_34)`. Solving the four vertex-sum
differences directly from the `K5` incidence matrix gives pivot determinant
one and the complete substitution

```text
x_03=a+c+d-t-1,       x_04=-a+2b+c+d-2,
x_12=b+c+2d-t-2,      x_13=-a+b+c+t-1.
```

The [checker](check_k5_variable_ratio_qq.py) substitutes these expressions
into the four vertex cube-difference equations and computes a completed
28-element grevlex Groebner basis over the exact field `Q(t)`. It then
reduces successive factors of

```text
P_9=product_(j=1..9)(x_01^2-x_j^2),
```

with the ten edges in lexicographic order. Each step multiplies the
preceding remainder by the next factor, so the final remainder is also the
normal form of the whole product. The first eight remainders are nonzero;
the ninth is zero. Thus `P_9` belongs to the ideal after extension to
`Q(t)`.

This is an exact generic-field identity, stronger than the earlier
`t=2` computation. Field membership alone does not settle every
specialization at poles of a rational certificate. In particular, the
coefficient denominators of a final Groebner basis do not by themselves
control denominators in its representation by the original generators.

The companion [denominator diagnostic](check_k5_variable_ratio_denominator_diagnostic.py)
reported the accumulated denominator from the sequential division
quotients as

```text
764411904000*(t - 3)*(t - 1)^3*(3*t - 1)^3.
```

Inspection of the completed basis gave coefficient denominators dividing

```text
2400*(t - 3)*(t - 1)^2*(3*t - 1).
```

These are pole diagnostics, not a proved complete exceptional-ratio list:
the scripts do not track the transformation coefficients from the computed
basis back to the original four cubics. A separate polynomial-ring
computation over `Q[t,a,b,c,d]` completed a 32-element basis, but its
sequential reduction was interrupted after prefix seven. No complete
polynomial-ring certificate from that attempt is claimed.

The generic computation also concerns zeros over an algebraic closure of
`Q(t)`. The [independent real-order proof](k5_real_first_third_moment_exclusion.md)
now excludes distinct-magnitude nonzero real solutions at **every** ratio,
without using these pole diagnostics. Neither result is a proof of the
uniform lattice-circle bound.

The completed computations used SymPy 1.14.0. Their interpreter in this
workspace was
`/Users/cxc/.cache/uv/environments-v2/high-level-tutorial-14eec61d3278d4de/bin/python`.
