# Constraints on actual and virtual radial normals

These elementary restrictions test ways to remove the anchor hypothesis
in [the direct radial growth theorem](anchor_radial_deficit_squareclass_growth.md).
They do not supply a radius-independent bound for arbitrary arcs.

## 1. Two actual anchors cannot both have very short primitive directions

Write actual circle points as `z_i=g_i v_i`, where `g_i` is the positive
ordinary coordinate gcd, `v_i` is primitive integral, and `|z_i|=R`.
Put `rho_i=|v_i|=R/g_i`. For distinct nonantipodal points on an arc
of length `L`, the determinant is a nonzero integer and

```text
1 <= |det(v_i,v_j)|
  = rho_i rho_j |sin(theta_i-theta_j)|
 <= rho_i rho_j L/R.                                (1)
```

Thus `rho_i rho_j>=R/L`. For `L<=C sqrt(R)<pi R`, the arc contains
no antipodes, so it contains at most one point with

```text
rho_i < R^(1/4)/sqrt(C).                             (2)
```

In particular, two polylogarithmic primitive radial lengths cannot occur
in the same fixed-constant endpoint arc for large `R`. Coupling the direct
small-deficit theorem at two such actual anchors cannot remove its
hypothesis. This says nothing against coupling a short direction to a
longer one.

## 2. A uniformly short virtual normal need not lie inside the angular arc

Set `epsilon=R^(-1/2)<=1/10`. Consider an arc whose midpoint angle is
`2epsilon` and whose angular half-width is `epsilon`. Its length is
`2sqrt(R)`. Every primitive integer direction in its angular interval
`[epsilon,3epsilon]` has coordinates `(q,p)` with positive integers
`p,q`. Since `tan(3epsilon)<=4epsilon`,

```text
|v| >= q = p/tan(arg v) >= 1/(4epsilon) = sqrt(R)/4. (3)
```

The tangent inequality follows from `sin(3epsilon)<=3epsilon` and
`cos(3epsilon)>=1-(3epsilon)^2/2>3/4`. For every proposed constant
`K`, (3) exceeds `K R^(1/4)` for sufficiently large `R`.
Consequently no uniform rule supplies a primitive normal of length
`O(R^(1/4))` within one endpoint angular half-width of every midpoint.
The arc position may depend on `R`, as required by the original
uniform-in-position problem. This is an angular selection counterexample,
not a many-point circle construction.

## 3. Integer projection still gives a useful alternative test

Let `v` be any nonzero integer vector, with length `rho`, making angle
`delta` with an arc midpoint. Parameterize the arc by angles
`delta+u`, where `|u|<=a=L/(2R)`. Taylor's theorem gives

```text
cos(delta+u)=cos(delta)-u sin(delta)+E(u),
|E(u)|<=u^2/2.
```

The linear term has range at most `2a|sin(delta)|`; the remainder has
range at most `a^2`. Therefore the range of integer projections obeys

```text
width(v dot z) <= rho(L|sin(delta)|+L^2/(4R))=:B_v.  (4)
```

Every projection level meets the circle in at most two points. Hence

```text
m <= 2(floor(B_v)+1).                               (5)
```

For the arc in Section 2, the short normal `(1,0)` is outside the
angular interval but is still useful: `|delta|=2epsilon` and
`L=2sqrt(R)` give `B_v<=5`. Thus that selection counterexample itself
has at most twelve lattice points by (5). Failure to find a short normal
inside the angular arc does not imply failure of a projection bound.

What is missing is a theorem that selects an integer normal with useful
projection width for every possible many-point endpoint configuration,
or supplies a stronger arithmetic argument in the complementary case.
The validity of (5) alone is not such a selection theorem. The earlier
[primitive-chord transform](endpoint_primitive_chord_transform.md) is a
related exact projection construction, with its remaining dependence
on the chord's coordinate gcd explicitly retained.

The [checker](check_radial_virtual_normal_constraint.py) verifies actual
integer-circle determinant and short-direction constraints, the explicit
angular obstruction, and sampled projection-width inequalities. The
proofs above, rather than the finite samples, establish the general
statements.
