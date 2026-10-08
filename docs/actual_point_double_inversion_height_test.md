# An actual-point double inversion is a nonconformal cotangent scaling

This is an exact test of a concrete inversion construction. It supplies a
nonconformal rational projective map, but does **not** show that its image is
another endpoint configuration. The uniform radius-independent point bound
remains unproved.

Let the source circle have squared radius `N>0`, and choose an actual lattice
point `a` on it. For a center `c` and positive scale `h`, write

```text
I_(c,h)(z)=c+h(z-c)/|z-c|^2.
```

The first inversion `I_(a,tN)` with `0<t<4` maps the circle, which passes
through `a`, to a line. The second inversion
`I_(-a,(4-t)N)` maps that line back to the **same** centered circle.
Both centers are lattice points on the source circle, and `-a` is available
whether or not it belongs to the short source arc. Taking rational `t`
makes both scales rational, so the composition maps rational circle points
to rational circle points.

Rotate coordinates only for the calculation, so `a=sqrt(N)` is on the
positive real axis. Put `w=z/a=q(X)=(X+iL)/(X-iL)`, with the anchor at
`X=infinity`. Directly on the circle,

```text
(w-1)/|w-1|^2=-1/2+iX/(2L).
```

The first image is the line `Re(u)=1-t/2`, where
`u=1-t/2+i tX/(2L)`. Set `k=(4-t)/t`. The second inversion gives

```text
w'=-1+(4-t)(u+1)/|u+1|^2
   =q(k L^2/X).
```

The original anchor maps to `-1`. Multiply the entire image by `-1`
(a Gaussian unit) and conjugate it to restore a positive one-sided
cotangent chart. The resulting map is exactly

```text
q(X) -> q(X/k).
```

Thus the physical two-inversion operation is, after the harmless common
reanchoring and reflection, the projective half-angle matrix
`diag(1/k,1)`; equivalently it keeps `X` and replaces `L` by `kL`.
It is nonconformal when `k!=1`, i.e. `t!=2`. In particular `t=3`
uses integer scales `3N,N` and produces `X -> 3X` with `L` fixed,
compressing the small angle near the anchor by a factor of three.

The exact least squared radius of the transformed tuple is obtained from
the existing **all-edge** integer lcm formula, not from the physical
inversion scales. For `t=3`, the anchor edges use `Q'_0i=3X_i`, and every
other edge uses the rational cotangent

```text
Q'_ij=(9X_i X_j+L^2)/(3(X_i-X_j)).
```

For each rational edge cotangent `Q'=u/v`, reduce the integer pair
`(u,Lv)` by its ordinary gcd. Its norm contribution is
`(a^2+b^2)/epsilon`, where `epsilon=2` precisely when the reduced
`a,b` are both odd. The ordinary lcm of all these contributions is
`N'`, exactly as in `integer_cotangent_lcm_height_target.md`.
The transformed normalized span is

```text
C'_*=2 (N')^(1/4) arctan(L/(3 A)),  A=min_i X_i.
```

For the existing four-point Pell instance at `n=1`, one has
`L=24`, `X=(3723,3456,7184)`, `N=358853785`, and the formula gives

```text
N'=2085984659911668125.
```

Its source normalized span is about `1.912`, whereas the compressed
image's normalized span is about `175.944`. Exact checks at Pell indices
`11` and `21` give decimal digit pairs `(84,169)` and `(160,320)` for
`(N,N')`; the numerical logarithmic ratios are `2.014` and `2.008`.
These are finite checks, not a proof of the Pell-family asymptotic.
They show that angle compression can be overwhelmed by the least-radius
cost, even when both inversion scales are integral and centered at actual
lattice points.

This construction is a special case of the rational projective maps
already bounded in `mobius_radius_content_bound.md`: for `t=3`, its
primitive matrix is `diag(3,1)`, with height three and nonconformal
invariant `K=64`. Existing fixed-unit source estimates prohibit a
comparable-radius image for sufficiently large fixed-unit clusters with
at least five points. That estimate does not supply an image-existence
theorem or settle arbitrary-unit configurations. For a global uniform
route, one would still need a choice of inversion parameters and centers
whose **all-edge** lcm stays sufficiently small for every large endpoint
tuple; the example above supplies no such control.

The [exact checker](check_actual_point_double_inversion.py) verifies 75
rational inversion identities, the three all-edge Pell radius calculations,
and the four-point selective-projective-map example in the companion
[reflection note](chord_reflection_projective_no_go.md). Decimal span and
logarithmic-ratio statements are finite numerical diagnostics only.
