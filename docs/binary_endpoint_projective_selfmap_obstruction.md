# Large binary endpoint tuples can have no projective self-symmetry

Let `z_0,...,z_(m-1)` be a primitive tuple of distinct Gaussian integral
points on squared radius `N`, in an arc of length at most `2 N^(1/4)`. Suppose

```text
m>=128,       gcd(m,6)=1,       N>2^176,
gcd_Z(Re z_j,Im z_j)=1 for every j.                  (1)
```

Then its rational half-angle configuration has **trivial projective
automorphism group**. Thus permuting the given points cannot provide the
nonconformal map needed by the reciprocal-stretch divisor theorem in
this class. This is a conditional structural obstruction, not a
construction of arbitrarily large endpoint tuples.

More generally, under the cardinality hypotheses in (1), every
nonidentity projective self-map fixes an actual source point. If the map
is nonconformal, that fixed point has ordinary content at least
`2^(-22) N^(1/8)`. If it is conformal, its content is `sqrt(N)`.

## 1. The finite projective group forces an actual fixed point

Let `X` be any set of `m>=3` distinct rational projective directions.
Its projective automorphism group embeds in the finite permutation group
of `X`, since fixing three directions makes a projective map the identity.
The orientation-preserving subgroup acts by cyclic shifts of the cyclic
order of `X`, hence is cyclic. Its action is free: a nonidentity finite
order orientation-preserving real projective map cannot have a real fixed
point. Thus its order divides `m`.

An orientation-preserving finite order element of `PGL_2(Q)` has order
in `{1,2,3,4,6}`. For a nonidentity element of order `n`, the eigenvalue
ratio is a primitive `n`th root `zeta`; the rational matrix invariant

```text
trace(M)^2/det(M)=2+zeta+zeta^(-1)
```

is rational. Since `zeta+zeta^(-1)` is a real algebraic integer in
`[-2,2]`, it belongs to `{-2,-1,0,1,2}`, giving precisely those orders.

When `gcd(m,6)=1`, the orientation-preserving subgroup is trivial.
There is therefore at most one nonidentity automorphism, which must be
an orientation-reversing involution. Since `m` is odd, its permutation
has a fixed point. A nonidentity projective involution has at most two
fixed directions, so its odd number of fixed sampled points is exactly
one. Call that actual point `j`.

## 2. A fixed-point involution is an actual-anchor shear

Reanchor both source and its permuted target at `j`. A primitive integral
representative is triangular. A nonidentity projective involution has
trace zero; changing projective sign if necessary therefore gives

```text
F=((a,b),(0,-a)),       a>0, gcd(a,b)=1.
```

Conjugating the entire physical target tuple gives

```text
M=((a,b),(0,a)).                                    (2)
```

Both least squared radii are still exactly `N`, both endpoint constants
are unchanged, and the anchor remains an actual sampled point. If
`b!=0`, this is precisely the nonconformal actual-anchor parabolic case.

The shear (2) has determinant core exactly `a^2`, with no anchor-scale
factor. For `b` odd, use reduced Gaussian coefficients
`u=2a-ib,v=ib`: their norms are coprime and their difference is `4a^2`,
with normalization factor `epsilon=1`. For `b` even, `a` is odd; use
`u=a-ib/2,v=ib/2`, whose norms are coprime and differ by `a^2`, with
`epsilon=4`. In both cases `k_j=1` and `c=epsilon|Norm(u)-Norm(v)|/4=a^2`.
Equivalently, all coefficient intervals contain the anchor coordinate
zero. This controls the diagonal at this particular fixed point and does
not combine estimates from two different anchors.

The arbitrary-unit full-core theorem gives, with
`F_m=floor((m-1)^2/4)`,

```text
a <= A_m N^h_m,
A_m=2^[m(m-1)/(2F_m)],       h_m=m/(8F_m).           (3)
```

This uses the [full source-and-target parity bound](mobius_arbitrary_unit_retention_addendum.md), applied to the original
tuple and its conjugate permutation.

## 3. The fixed point must have large coordinate content

Let `q_m,u_m` be the erased-conductor and finite-sample condition exponents
from [the critical triangular theorem](mobius_finite_sample_critical_form.md).
The [parabolic anchor divisibility](parabolic_shear_anchor_content_obstruction.md)
and the finite-sample condition bound imply

```text
Nclip | |b|(b^2+4a^2) g_j,
Nclip >= 2^(-m q_m) N^(1-q_m),
|b|(b^2+4a^2) <= 40 a^3 (2N)^[(3/2)u_m].
```

Using (3), define

```text
eta_m=1-q_m-(3/2)u_m-3m/(8F_m)
     =1/4-7/m+O(m^(-2)),
K_m=40 * 2^[m q_m+(3/2)u_m+3m(m-1)/(2F_m)].
```

Then the actual fixed point satisfies

```text
g_j >= K_m^(-1) N^eta_m.                            (4)
```

For completeness, all constants in (1) follow from elementary uniform
estimates. For `m>=128`, the established bounds give

```text
q_m<=8/m,       u_m<=1/2+2/m,
3m/(8F_m)<=2/m,       m(m-1)/F_m<=5.
```

Hence `eta_m>=1/4-13/m>=1/8` and `K_m<2^22`. Thus

```text
g_j >= 2^(-22) N^(1/8).                             (5)
```

Under (1) this contradicts `g_j=1`.

## 4. Conformal reflections also force content at their actual fixed point

If `b=0`, the involution at the actual anchor is conjugation of phases.
At each odd split source prime, its permutation sends every signed
phase valuation `t_i` to `-t_i`. The source interval
`[-e_j,e-e_j]` is therefore symmetric, so `e_j=e/2` at every prime.
Both interval endpoints are attained: collective Gaussian primitivity
gives minimum valuation zero at each of the two orientations, and
`v_pi(z_i)+v_bar(pi)(z_i)=e` then gives maximum valuation `e`.
Consequently

```text
g_j=product p^(e/2)=sqrt(N).                         (6)
```

This also contradicts the binary hypothesis when `N>1`. Equations
(5)--(6) and the finite-group classification prove the theorem.

## 5. What this says about minimizing a projective orbit

The [existing orbit-radius rigidity theorem](projective_orbit_radius_rigidity.md)
controls the radii of two endpoint representatives already known to
exist. It does not prove the existence of a second representative.
The argument here is different: it excludes every nonidentity projective
permutation of one large binary configuration with the stated cardinality.

Taking a minimum over endpoint representatives in an orbit guarantees a
minimum because squared radii are positive integers. It supplies no second
minimizer or nonconformal competitor. Common rotations, reflections, and
actual reanchoring preserve the conformal character of the same underlying
map. They do not escape that missing-existence step. The theorem above
shows that an appeal to a nontrivial permutation of the minimizing tuple
would fail in the binary class (if that minimizer itself satisfies (1)).

Distinct configurations at the same radius in the same projective orbit
are not classified here, and minimizing over the orbit need not preserve
the binary property. No conclusion excludes a construction using a new
image configuration with controlled determinant, clipping, and endpoint
constant. That broader bridge remains open.

The [checker](check_binary_endpoint_selfmap.py) verifies the fixed-anchor
matrix normalization, both coefficient parities, exact core values,
reflection valuation symmetry, and the explicit exponent certificates.
