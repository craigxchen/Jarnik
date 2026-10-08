# Affine Walsh flip assignments contain large phase certificates

This note supplies aligned affine planes and larger affine cosets for the
repeated Walsh flip profile studied in
[the linear-support character note](linear_support_character_height_obstruction.md).
The original construction below gives dimension `floor(t/3)` and the
`3/4` growth exponent. A sharper universal binary-bilinear theorem,
proved in [the isotropic-subspace note](walsh_arbitrary_binary_bilinear_isotropic_bound.md),
gives dimension `floor((t-2)/2)` and strengthens the same nearby-prime
affine-profile growth exponent to `2/3`; the final section records this
improvement. The note covers affine assignments and a controlled number
of changed rows. It does not assert that an arbitrary nonlinear
assignment has such a plane, or that the sign profiles themselves occur
on endpoint arcs.

Let `G=F_2^t`, `t>=6`, `M=2^t`, and use the ordinary binary dot product.
Start with an affine label assignment

```text
a(x)=Lx+b,                 rank L >= t-1.             (1)
```

The labels may initially include zero. In particular, every full-rank
affine map does include zero; changing that row to a nonzero label is one
of the exceptions allowed below. A plane `x0+V`, `dim V=2`, is **aligned**
if all four labels have the same nonzero restriction `l` to `V`:

```text
a(x0+r)|V=l != 0            for every r in V.        (2)
```

**Plane certificate.** Assignment (1) has at least `2^(t-3)` disjoint
aligned planes, all parallel to one direction `V`. Changing labels at
fewer than `2^(t-3)` rows leaves at least one of these planes aligned.
The changed labels need not themselves be affine. This counts actual
planes, not just four-row vectors in a character lattice.

To prove it, set `Q(u)=u dot Lu`. Its polar form is
`B(u,v)=u dot Lv+v dot Lu`. Every quadratic form on a three-dimensional
binary space has a nonzero isotropic vector. Indeed, if all seven
nonzero vectors had `Q=1`, then for a basis `e1,e2,e3` the three pair
sums force `B(ei,ej)=1`; the quadratic identity gives
`Q(e1+e2+e3)=3+3=0` in `F_2`, a contradiction. Apply this fact in any
three-dimensional subspace of `G` and choose nonzero `u` with `Q(u)=0`.

Now put

```text
H={v: u dot Lv=0 and v dot Lu=0}.                 (3)
```

We have `u in H` and `dim H>=t-2>=4`. Choose a three-dimensional
subspace of `H` avoiding `u`. It contains a nonzero `v` with `Q(v)=0`.
Thus `V=span(u,v)` has

```text
r dot Ls=0                 for every r,s in V.      (4)
```

For `r,s in V`, equation (4) gives
`a(x+r) dot s=(Lx+b) dot s`. The restriction therefore stays constant
through each coset `x+V`. The induced affine map from the `2^(t-2)`
cosets to `V*` has linear rank at least one: if its linear part were zero,
`V` would lie in `ker L^T`, whose dimension is at most one by (1).
Consequently at most half the cosets have zero restriction, proving the
count. The good planes are disjoint, so one edited row can spoil at most
one of them. Their initial labels are automatically nonzero. A patch of
the single zero label of a full-rank map is harmless for `t>=6`.

The same construction gives a larger certificate. Set
`d=floor(t/3)` and `h=2^d`. There is a `d`-dimensional direction `V`
on which `r dot Ls=0` for every `r,s in V`. At least
`2^(t-d-1)` disjoint cosets `x+V` have a constant nonzero label
restriction to `V`, and changing fewer than that many row labels leaves
one. Indeed, if a totally isotropic subspace `U` has dimension `j`,
the analogue of (3), with both conditions imposed for every `u in U`,
contains `U` and has codimension at most `2j`. Its quotient by `U`
has dimension at least `t-3j`. For `j<d` this is at least three.
Choose a three-space in it, lift it to a three-space disjoint from `U`,
and use the quadratic fact above to add an isotropic vector. Finally
`dim ker L^T<=1<d`, so the restriction map on the cosets has positive
rank and at least half of them have nonzero restriction. The plane
statement is the case `d=2` when `6<=t<=8`.

## Consequence for repeated physical Walsh copies

For each nonconstant Walsh class `h_a(x)=(-1)^(a dot x)`, use `b_copy`
physical columns, with `b_copy>=5` for the strict-obtuse construction.
Flip one **distinct physical column** at every row `x`; its class label
is the resulting `a_x`. A sufficient capacity condition is that each
label is assigned to at most `b_copy` rows. This condition, along
with nonzero labels, is required when applying the plane certificate to
the physical prime profile. The stronger fiber bound `<=2` in the
[linear-support note](linear_support_character_height_obstruction.md)
is used for its separate grouped-column lattice-minimum statements;
the following coset calculation uses the distinct physical columns
directly.

Choose one of the aligned `d`-dimensional cosets `P=x0+V` with common
restriction `l!=0` and set `c_(x0+r)=(-1)^l(r)` on `P`, zero
elsewhere. Then `sum c_x=0`, `||c||_1=h`, and its half-integral row
certificate is integral in every physical column. For an unflipped
Walsh column,

```text
u_a=(1/2)sum_x c_x h_a(x)
   =(h/2) h_a(x0) if a|V=l, and 0 otherwise.      (5)
```

Exactly `M/h` classes qualify. At each of the `h` assigned physical
columns on `P`, the flip changes the corresponding absolute coefficient
from `h/2` to `h/2-1`; flips at other rows contribute nothing. Thus the
exact physical coefficient sum is

```text
sum_(physical j) |v_j| = b_copy M/2-h,
max_j |v_j|<=h/2.                                    (6)
```

Assume the physical columns use distinct odd split primes with
`ell=min_j log p_j` and `max_j log p_j<ell+1/r`, where
`r=b_copy(M-1)`, and their literal Gaussian orientations are
conjugate-primitive as in the linked note. Put `W=sum_j log p_j`.
If the whole source tuple has an additional common Gaussian factor
`d_common`, charge `D=log Norm(d_common)>=0`, so its squared radius is
`N=exp(W+D)`. This factor cancels from the signed row product below.
Equation (6) gives

```text
V(Beta)-W/2 < -(h-b_copy/2)ell+1.                  (7)
```

The certificate has denominator `2`, so its phase targets the
elementary `pi/4` grid. Its `h/2` positive and `h/2` negative row
coefficients give `|sum_x c_x theta_x|<=h Delta/2` on an arc of width
`Delta<=C exp(-(W+D)/4)`. The composite `Beta` is nonunit and
conjugate-primitive: its nonzero exponents use disjoint physical split
primes in only one orientation. The angular distance of `arg Beta`
from the grid is at least `1/sqrt(2 Norm(Beta))`, hence any actual
endpoint realization requires

```text
V(Beta) >= (W+D)/2+log 8-2 log(hC),
b_copy > 2h-2(1-log 8+2 log(hC)-D/2)/ell.        (8)
```

For every fixed `C` and sufficiently large `M`, `ell>=log 5` makes
(8) imply `b_copy>h>=M^(1/3)/2`. Since the `r` rational primes are
distinct, `r=b_copy(M-1)>=M^(4/3)/4` and
`W>=log((r+1)!)>=(r/2)log(r/2)`. Thus
`log N=W+D>=c M^(4/3)log M` for a positive absolute `c` at large
`M`. Inverting this increasing bound gives any actual endpoint
realization of this near-affine, nearby-prime profile

```text
M=O_C((log N/loglog N)^(3/4)).                    (9)
```

No assumption that `ell` is approximately `log M` is needed for (9).
For the original five-copy construction, `b_copy=5`, so (8) instead
excludes all sufficiently large `M` at fixed `C`. At `d=2`, (6)
specializes to `5M/2-4`, and the linked note gives the sharper
quantitative constant `C_source>exp(3ell/4)/(sqrt(2)e)`.

The conclusion applies to any modified assignment for which an aligned
coset survives and the nonzero-label/distinct-physical-column conditions
still hold. The nearby-prime logarithmic spread is also a hypothesis
of (7)--(9). This is a growth constraint for this structured prime
profile; it gives no general latticepoint count.

The [checker](check_walsh_affine_assignment_plane_certificate.py) tests
the quadratic construction, rank and coset count for literal binary
matrices, exception survival, and the physical coefficient formula (6).

## Sharper universal isotropic direction and growth exponent

The original greedy isotropic direction above is valid, but its dimension
is not optimal as a universal lower bound. For any binary bilinear form
`B(u,v)=u dot Lv`, the [universal theorem](walsh_arbitrary_binary_bilinear_isotropic_bound.md)
provides a direction `V` with

```text
d=dim V=floor((t-2)/2),
h=|V|=2^d >=sqrt(M)/(2sqrt(2)),       M=2^t.       (10)
```

Take `t>=6`, `rank L>=t-1`, and let `a(x)=Lx+b` before edits. Since
`B|VxV=0`, every coset `x+V` has a constant label restriction to `V`.
The restriction map on the `M/h` cosets has a nonzero linear part:
otherwise `V subset ker L^T`, whereas `dim ker L^T<=1<d`. Hence at
least half the cosets are disjoint aligned cosets. Changing fewer than

```text
M/(2h)=2^(t-d-1)                                       (11)
```

row labels leaves one intact. As before, the final labels must be
nonzero and each row must flip a distinct physical column; enough
copies must be available for their label fibers.

The physical character and rational-angle argument (5)--(8) apply to
this surviving coset without modification. If actual distinct physical
prime norms obey the nearby-prime log spread in that section and an
endpoint arc has a fixed normalized constant `C`, then (8) forces
`b_copy>h` for all sufficiently large `t`. Thus the number of distinct
physical prime norms satisfies
`r=b_copy(M-1)>=c M^(3/2)` for an absolute positive `c`.
Ordering those primes gives
`log N>=log((r+1)!)>=c' M^(3/2)log M`; any common Gaussian content
only increases `log N`. Inverting yields the strengthened bound

```text
M=O_C((log N/loglog N)^(2/3)).                         (12)
```

This conclusion holds for the stated rank-at-least-`t-1` affine map and
for assignments with fewer than `M/(2h)` row-label edits, under the
same copy-capacity and actual nearby-prime hypotheses. It makes no
claim about arbitrary nonlinear assignments or a generic endpoint
count.
