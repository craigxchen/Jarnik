# Small matching residuals reduce to a rational projection

This audits one possible use of the exact cut matching. It is not a general
obstruction to the matching method, and gives no new unrestricted growth bound.

Let distinct `z_i=beta_i w_i` be Gaussian integers of norm `N=R^2`, with
`Norm(beta_i)=Q`, `Norm(w_i)=N/Q`, and `|z_i-z_0|<=L`, with `L>0`. Here `Q|N` can be any
available cut product. For nonzero Gaussian integers `a,b`, put

```
r_i=b w_i-a conjugate(beta_i),
H=(|b| L+|b z_0-a Q|)/sqrt(Q),
v=a conjugate(b),
D=Norm(b) N/Q+Norm(a) Q.
```

Then, exactly,

```
|r_i| <= H,
Norm(r_i)=D-2 Re(conjugate(v) z_i).                 (1)
```

The first identity follows by multiplying `r_i` by `beta_i`:
`beta_i r_i=b z_i-aQ`. The second follows by expanding the norm.
Neither primitivity nor separation of source primes is needed.

For each fixed nonzero residual vector there are at most two distinct
products: the equal-norm condition cuts the beta-circle by a line. More
strongly, each fixed residual **norm** gives at most two products directly
from (1), since `v!=0`. Zero residual gives only `z_i=aQ/b`.
Consequently the easy bound is

```
m <= 2 floor(H^2)+1.                              (2)
```

Every residual norm is also congruent to `D` modulo two. If `D` is even,
the sharper bound is `m<=2 floor(H^2/2)+1`; if `D` is odd, it is
`m<=2 floor((H^2+1)/2)`. These count available integer projection levels,
not independent conditions coming from different cut factors.

## Why varying the cut does not automatically give bounded H

The displayed sufficient functional H imposes a restriction independent of
how Q was selected. Since `|b|L<=H sqrt(Q)`,

```
|ab| <= |b|^2 R/Q+H |b|/sqrt(Q)
     <= H^2 R/L^2+H^2/L.                          (3)
```

The first inequality uses `|a|Q<=|b|R+H sqrt(Q)`.
Suppose the anchor z_0 itself has the same factorization, as it does in the
cut construction. Its residual satisfies (1), and AM--GM gives
`D>=2R|v|`. Therefore

```
2R |v| (1-cos(arg(z_0)-arg(v))) <= H^2.            (4)
```

For `L=C sqrt(R)` and bounded H, (3) bounds the Gaussian normal v by
`H^2/C^2+o(1)`, even for unbalanced Q. Equation (4) forces the anchor
direction within `O(R^(-1/2))` of one of finitely many such normals.
Thus a theorem asserting that some cut always makes this functional
bounded would already have to prove this strong rational-direction
concentration for every sufficiently crowded configuration. A generic
Dirichlet argument cannot supply it.

If one instead passes to k points in a subarc of relative length k/m,
then (3) permits normals of size `O(H^2 m^2/(C^2 k^2))`. The required
angular accuracy in (4) remains; shrinking alone does not produce it.

The bound is a special form of the rational-normal projection estimate in
`multipoint_continuation.md`, Section 1. It does not extend the growing
primitive-radial-height Pell slice of `endpoint_independent_attack.md`:
that result exploits an actual lattice anchor in a rational radial
position and yields a fixed count while its allowed height grows.

The residual vectors retain more information than their norms, and several
residual systems might conceivably interact. Equations (1)--(4) rule out
only treating a bounded single residual functional as a direction-free
new global constraint.

## The split quadratic also retains the original half-angle coordinate

A different proposal is to use the four integer coordinates of `(beta,w)`
and the quadratic equation `Q Norm(w)-(N/Q) Norm(beta)=0` near the matching
plane. This quadratic form is rationally split: the matching plane
`w=(z_0/Q) conjugate(beta)` itself is rational and totally isotropic.
Its large rational denominator should not be replaced by irrationality.

Explicitly, define the Gaussian integers

```
u_i=Q w_i-z_0 conjugate(beta_i),
v_i=Q w_i+z_0 conjugate(beta_i).
```

Then

```
Re(u_i conjugate(v_i))=0,
Norm(u_i)+Norm(v_i)=4QN,
u_i (z_i+z_0)=v_i (z_i-z_0).                     (5)
```

Writing `u=x+iy`, `v=s+it`, the integer matrix `[[x,y],[-t,s]]`
has determinant zero. This is an exact rank-one constraint, but its
scalar ratio is

```
u_i/v_i=(z_i-z_0)/(z_i+z_0),                      (6)
```

whenever the denominator is nonzero. In particular it is a purely
imaginary rational number independent of the chosen cut. It is precisely
the original rational tangent-half-angle coordinate at the anchor.
If `z_i+z_0=0`, then (5) gives `v_i=0` and
`u_i=-2 z_0 conjugate(beta_i)`, which is nonzero. This is the single
projective tangent-half-angle value at infinity omitted by the quotient
notation; it also carries no selected-cut dependence.
Passing to this ratio erases all selected-cut information. The factors
in an integral rank-one decomposition of the matrix can still retain
that information; their compatibility across different rows is not
excluded by this audit. It would be necessary to preserve those factors,
or add another simultaneous arithmetic condition, to obtain a new
constraint from this presentation of the quadratic.

`check_cut_matching_residual_rational_projection.py` exhaustively checks
the norm identities and at-most-two norm fibers on small Gaussian circle
products, together with the integer rank-one and cut-independent ratio
identities. These finite checks supplement the proofs above.
