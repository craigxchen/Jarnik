# The mirror denominator, edge conductor, and point contents

Reflecting a short cluster about the midpoint of its extreme angles
preserves its angular interval. The reflected points can always be
adjoined after clearing Gaussian denominators. Failure of integer
cotangent integrality at the old common denominator is therefore not,
by itself, an obstruction to this construction. The exact radius cost
below is the quantity that has to be controlled.

The least-union formula and radial counting obstruction are already
proved in the [reflection-union audit](reflection_union_radial_count_obstruction.md).
This note connects that existing cost to the ordinary edge conductor
and the two physical point contents through (2), and to primitive chord
parity through (3). The repeated construction and Pell fixture are
normalization checks, not a new symmetrization strategy or growth result.

Let distinct Gaussian integers `z_i` form a collectively primitive
tuple of common squared norm `N`, and select two distinct points
`z_0,z_a`. Write `g_i=gcd_Z(Re z_i,Im z_i)`. Define

```text
c=z_0 z_a/N=z_a/conjugate(z_0),       |c|=1,
w_i=c conjugate(z_i),
c=alpha/beta in reduced Gaussian terms,
D=Norm(beta).
```

The least squared radius of the combined original and reflected tuple is

```text
N_union=N D.                                          (1)
```

Indeed, `beta z_i` and `alpha conjugate(z_i)` are integral, of norm
`ND`. Their two collective Gaussian gcds are respectively `beta` and
`alpha`; these are coprime. The combined tuple is primitive. Equivalently,
any multiplier making the original primitive tuple integral must be a
Gaussian integer, and making all the reflected points integral then
requires divisibility by `beta`.

There are two useful exact descriptions of `D`. Let `n_0a` be the norm
of the reduced Gaussian denominator of `z_a/z_0`, which is one of the
all-edge norms in the [integer cotangent formulation](integer_cotangent_lcm_height_target.md).
Then

```text
D n_0a gcd_Z(g_0,g_a)^2=N.                            (2)
```

At a split prime `p=pi bar(pi)` of exponent `e` in `N`, put
`u=v_pi(z_0)` and `v=v_pi(z_a)`. The exponents in `D`, `n_0a`, and
`gcd(g_0,g_a)` are respectively

```text
|u+v-e|,       |u-v|,       min(u,e-u,v,e-v).
```

Their combination in (2) is `e`, by
`|u+v-e|+|u-v|+2min(u,e-u,v,e-v)=e`. No inert or ramified
factor occurs in the common norm of a collectively primitive tuple.

Alternatively, write the chord

```text
z_a-z_0=g(u+iv),       gcd_Z(u,v)=1,
epsilon=2 if u,v are both odd, and 1 otherwise.
```

Since `(z_a-z_0)/conjugate(z_a-z_0)=-c`, Gaussian reduction gives

```text
D=(u^2+v^2)/epsilon.                                 (3)
```

Thus this mirror denominator is exactly the squared primitive-chord
length, with the ramified parity factor retained. It is not a new
independent small conductor. Equation (3) identifies the connection to
the [primitive-chord strip transform](endpoint_primitive_chord_transform.md).

If `z_0,z_a` are the extreme points of a minor arc of angular width
`Delta`, the reflected set lies in that same interval, and the
extremes are still attained. Its exact normalized arc constant is

```text
C_union=Delta (N D)^(1/4)=C_source D^(1/4).            (4)
```

Consequently doubling most of the points by reflection does not come
with a uniform geometric cost unless `D` is controlled. After the common
Gaussian multiplication the symmetry axis is an axis or diagonal of the
integer lattice: reduced equal-norm coprime `alpha,beta` satisfy
`alpha=unit*conjugate(beta)`. Counting normal integer levels there
recovers the same primitive-chord estimate, rather than removing its
uncontrolled chord-direction factor.

For the [four-point Pell family](four_point_bonus_counterexample.md),
use its cotangents

```text
L=24,
X=(90UV+30V^2+3, 96UV, 160V^2+128UV+16),
U^2-5V^2=1,
a=1+10V^2+2UV, b=1+10V^2-2UV, c_P=13+130V^2+38UV.
```

At the permitted indices, `N=5ab c_P`, the extreme cotangent is
`A=96UV`, its reduced edge norm is `ab=16U^2V^2+1`, and the two
physical extreme contents are one. Therefore

```text
D=5c_P,        N_union=25ab c_P^2.                    (5)
```

Here `D` grows like `V^2`, even though the source is an actual endpoint
family whose normalized arc constant tends to zero. The original four
points and their reflection have six distinct points: an additional
coincidence would give a nontrivial pair-product collision, excluded by
the [multiplicative rectangle theorem](multiplicative_rectangle_separation.md)
once the source constant is below `2sqrt(2)`. The new normalized
constant stays bounded because `N_union` grows like `V^8` and the
angular width like `V^(-2)`. This is a fixed six-point family, not an
unbounded endpoint cluster or a growth-rate improvement.

The [checker](check_mirror_union_least_radius_cost.py) verifies the
reduced Gaussian fraction, combined primitivity, both exact formulas
for `D`, and the displayed Pell specialization.
