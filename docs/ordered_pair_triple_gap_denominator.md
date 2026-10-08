# The ordered cut as an exact pair–triple radius denominator

The endpoint-versus-interior forced factor has an index-free description
using two actual primitive subcircles. It is the reduced denominator of
their squared-radius product divided by the full squared radius. This
removes accidental triangle valuations from the sufficient five-point
target. It supplies no improved upper bound for that denominator.

Let `z_0,...,z_4` be five distinct equal-norm Gaussian integers, primitive
as a full tuple, with norm `N=R^2`. Let `O={0,4}` and `I={1,2,3}`. Define

```text
G_O=gcd_G(z_0,z_4),             A_O=Norm(G_O),
G_I=gcd_G(z_1,z_2,z_3),         A_I=Norm(G_I),
N_O=N/A_O,                     N_I=N/A_I.
```

Thus `N_O,N_I` are the squared radii after independently dividing each
group by its Gaussian gcd; this also rotates the two groups independently.
Their original angular spans are preserved. There is no assertion that
the two normalized groups remain on a common circle or a common arc.

At an odd split prime `p`, put `e=v_p(N)` and `t_j=v_pi(z_j)`.
Full primitivity gives `min t_j=0,max t_j=e`. Set

```text
U_O=[min(t_0,t_4),max(t_0,t_4)],
U_I=[min(t_1,t_2,t_3),max(t_1,t_2,t_3)],
h_p=distance(U_O,U_I),
l_p=length(U_O intersect U_I),
K=product_p p^h_p,             L=product_p p^l_p.
```

Empty intersection has length zero. These are integer lengths of real
intervals, so a shared endpoint contributes zero. The `h_p` is precisely
the threshold multiplicity of the cut `O|I` in the
[ordered index theorem](ordered_affine_circuit_conductor_index.md).

## Exact gap–overlap decomposition

Write `r_O=length(U_O),r_I=length(U_I)`. Since their convex hull is
`[0,e]`, elementary interval geometry gives

```text
r_O+r_I-e=l_p-h_p,
h_p=max(0,e-r_O-r_I),          l_p=max(0,r_O+r_I-e).
```

The primitive subgroup exponents are exactly `v_p(N_O)=r_O` and
`v_p(N_I)=r_I`. Consequently

```text
N_O N_I/N=L/K,                gcd(K,L)=1,             K L | N,

K=N/gcd(N,N_O N_I)
 =A_O A_I/gcd(N,A_O A_I),

L=N_O N_I/gcd(N,N_O N_I)
 =N/gcd(N,A_O A_I).                                      (1)
```

These formulas retain nested prime powers and arbitrary row units. No
triangle contents occur. Both `h_p` and `l_p` are at most `e` and cannot
be positive together; this proves the coprimality and divisibility in
(1). In terms of the actual primitive radii,

```text
(R_O R_I/R)^2=L/K.                                      (2)
```

Thus a large real value of this radius ratio controls its numerator
relative to its denominator. It does not control the denominator alone.

## What the actual arc geometry implies

Now assume geometric order on a minor arc of angular span `Delta`, and
write `Delta_I` for the angular span of its three interior points.
The primitive norm of every subgroup is odd: any ramified or inert
factor would still divide all of that subgroup's Gaussian coordinates.
Two points of odd integer norm have squared distance at least two.
Also every triangle determinant on such a circle is even: all rows
have `x+y` odd, so all differences have `x+y` even. Noncollinearity
therefore makes its absolute determinant at least two.

Applying these facts to the two independently primitive groups gives

```text
N_O >= 1/(2 sin^2(Delta/2)) >= 2/Delta^2,
N_I >= 16/Delta_I^3.                                    (3)
```

For the second inequality, the determinant of a triangle on squared
radius `N_I` is at most `N_I Delta_I^3/8`: the product of its three
chords is divided by twice its radius, and the product of successive
angular gaps with their sum is at most `Delta_I^3/4`.

Equations (1) and (3) imply the unconditional, opposite-direction
restriction

```text
L/K >= 32/(N Delta^2 Delta_I^3).
```

In particular, at `Delta<=C N^(-1/4)`,

```text
L >= 32 C^(-5) K N^(1/4).                              (4)
```

At `C=1/2` the constant is `1024`. In particular, replacing the actual
product `N_O N_I` by a quantity at most `N` is impossible on this arc
class. The overlap numerator is necessarily large.

The direct elementary upper bound for the gap factor remains exponent
`1/4`. Every `p^h_p` divides `A_I`: the interior allocation interval lies
entirely to one side of the gap, so `e-r_I>=h_p`. If `D_123` is the
actual interior triangle determinant, division by `G_I` gives
`|D_123|=A_I |D_123^primitive|>=2 A_I`. Hence

```text
K <= A_I <= |D_123|/2 <= N Delta_I^3/16
  <= (C^3/16) N^(1/4).                                (5)
```

This is a factor-two constant refinement of the usual nonzero-integer
triangle estimate, not an exponent improvement.

These exponents do not create new slack after sampling five rows from a
large configuration. For a formal fully fair inherited profile of ambient
weight `W`, the five-tuple's own primitive weight is `15W/16`, while
`log N_O=W/2`, `log N_I=3W/4`, `log K=W/16`, and
`log L=3W/8`. Thus `log(L/K)=5W/16`. If the inherited angular scale is
`exp(-W/4)`, the logarithmic leading term in (4), retaining the actual
smaller five-tuple arc constant, is also `5W/16`. Gaussian primitive
division preserves angular width and decreases the normalized arc
constant; replacing its intrinsic `N` by the ambient norm would lose
this bookkeeping. The fair-profile calculation is an allocation model,
not a claimed actual endpoint realization.

## The weaker sufficient arithmetic target

The [ordered uniformity criterion](ordered_conductor_index_uniformity_criterion.md)
works verbatim under the weaker hypothesis

```text
log K <= delta log N+B,             0<=delta<1/15,      (6)
```

for every primitive ordered five-tuple in the `C=1/2` arc class. Indeed
the inherited cut weight called `A_J` in that proof is **exactly**
`log K_J`, while the removed constant-layer weight is unchanged. Its
averaging argument therefore has the same `beta_M,gamma_M` and the
same explicit point-count consequences. The affine index theorem gives
`K|gcd(Q,N)`, so (6) does not require control of accidental valuations
of `gcd(Q,N)`.

Equivalently, (6) asks that the reduced denominator of
`N_O N_I/N` be at most `exp(B)N^delta`, or that
`gcd(N,N_O N_I)>=exp(-B)N^(1-delta)`. Neither (3), (4), nor (5)
proves this denominator estimate. This is an exact arithmetic
reformulation and a sharper statement of the missing input, not a
uniform endpoint theorem.

There is also an [exact oriented descent](ordered_gap_opposite_twist_descent.md):
the gap defines `alpha` with norm `K`, and opposite division by `alpha`
and its conjugate removes it at norm `N/K`. The two groups can then lie
in separate sectors. Restoring their interlacing is one phase interval,
so this descent alone supplies no stronger denominator bound.

The [checker](check_ordered_pair_triple_gap_denominator.py) exhausts small
allocation vectors and independently verifies the formulas on actual
Gaussian tuples, including nested powers, varying units, nontrivial
overlap and gap factors, and literal ordered fixtures.
