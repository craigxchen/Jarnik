# A whole-clique affine translation has no automatic descent

The identity

```text
ab-c(a+b)=L^2  <=>  (a-c)(b-c)=c^2+L^2
```

suggests a Vieta or frieze mutation. Here `a,b` are successive oriented arc
cotangents and `c=(ab-L^2)/(a+b)` is their sum-angle cotangent. The natural
whole-clique candidate is
an affine translation `X_i' = X_i+t`, with `L` fixed (or with a common scaling
of both `X_i,L`). This note gives its exact closure condition and shows why it
does not furnish a uniform height descent.

For every pair put `d_ij=X_i-X_j` and `c_ij=(X_iX_j+L^2)/d_ij`. Then

```text
c'_ij-c_ij = t(X_i+X_j+t)/d_ij.                 (1)
```

Therefore translation by `t` preserves the entire integer cotangent clique
whenever

```text
d_ij | t(X_i+X_j+t)                              (2)
```

for every pair. A simple sufficient condition is
`D=lcm_(i<j)|d_ij|` dividing `t`; it is exact closure, including all nested
prime powers. The condition is pairwise and does not follow from a single
triangle or from the frieze identity alone.

The anchor is the point at infinity, so its edge quotient is simply the
finite coordinate `X_i'`; there is no cotangent-zero anchor in this chart.
If `L'=uL` and `X_i'=uX_i+v`, the same calculation yields

```text
c'_ij = u c_ij + [uv(X_i+X_j)+v^2]/[u d_ij],
```

so common scaling `v=0` is only projective normalization; a genuine affine
shift again needs the explicit divisibilities, now with the displayed
numerator divisible by `u d_ij`.

## Exact lcm cost

Let `n(q,L)` denote the reduced Gaussian edge norm from the all-edge lcm
formula. Under a closed translation the new radius is exactly

```text
N' = lcm( n(X_i+t,L), n(c_ij+t(X_i+X_j+t)/d_ij,L) ).  (3)
```

There is no divisibility relation between `N'` and `N`: the additive term in
(1) can introduce new prime powers, while deleting no edge occurs. Thus this
whole-clique operation has no general monotonicity. The divisor complement is
the separate circle reflection that preserves both `N` and the normalized
angular span.

For a positive endpoint chart, a translation that lowers `A=min X_i` must use
negative `t` while retaining all `X_i+t>0`. Even when (2) holds, (3) may grow
dramatically. Conversely, a rare decrease in `N'` can leave the positive chart
or simply be a finite replacement of the type already covered by the exact
candidate update in [the mutation audit](luna_one_coordinate_mutation_audit.md).

The companion checker evaluates the exact update on small cliques and records
both signs of the canonical shift `t=D`. It verifies closure and the all-edge
lcm formula. Its examples disprove automatic monotonicity and admissibility
of this canonical shift; they do not rule out every translation-based
existence argument with a configuration-dependent choice of `t`.
