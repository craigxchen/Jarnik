# Overlapping ordered windows have coprime oriented gap factors

Opposite-twist descents of overlapping five-point windows obey an exact
joint divisibility relation. For adjacent windows, the product of their
gap norms divides the norm of the Gaussian gcd of the two common interior
points. The two twist phases can nevertheless cancel completely, as an
actual ordered six-point example shows.

This is a simultaneous arithmetic restriction, not a stronger endpoint
exponent. At arc constant `1/2`, its resulting size bound is weaker than
the product of the two existing interior-triangle bounds.

## Coprimality for crossing endpoint pairs

For an ordered five-subset `U`, first divide its points by their own
Gaussian gcd, and construct its oriented factor `alpha_U` as in the
[opposite-twist descent](ordered_gap_opposite_twist_descent.md). Write
`K_U=Norm(alpha_U)`. All factors below use the same choice of the Gaussian
prime `pi` above each split rational prime; arbitrary unit choices do not
affect the conclusions.

Suppose two five-subsets have crossing endpoint pairs `a<b<c<d`, with
`a,c` the endpoints of `U`, `b,d` the endpoints of `V`, and assume
`b` belongs to the interior of `U` and `c` to the interior of `V`.
Then

```text
gcd_G(alpha_U,alpha_V) is a unit.                       (1)
```

Indeed, use allocations `t_i=v_pi(z_i)` of the original common-circle
points. Separate subset normalization subtracts a constant from all
allocations in each subset, so neither the gap nor its orientation changes.
If both factors contained `pi`, the outer-above-interior inequalities
would imply both `t_c>t_b` and `t_b>t_c`. If both contained
`conjugate(pi)`, the inequalities would both reverse, again contradicting
one another. Thus a common rational prime can occur only in opposite
Gaussian orientations. The factors contain no ramified or inert primes,
so this proves (1) globally.

For consecutive five-row windows starting at indices `j` and `j+d`, this
applies for `d=1,2,3`. It need not apply at `d=4`, when the windows have
only a common endpoint and no crossing interior incidences.

## A product divisor on the common interior

If a point `z_i` is interior to both windows, both
`conjugate(alpha_U)` and `conjugate(alpha_V)` divide it. More precisely,
each divides that point after division by the corresponding subset gcd;
it therefore also divides the original point. Equation (1) gives

```text
conjugate(alpha_U alpha_V) | z_i.                       (2)
```

For adjacent windows `U=(0,1,2,3,4)` and `V=(1,2,3,4,5)`, their common
interior points are `z_2,z_3`. Hence

```text
conjugate(alpha_U alpha_V) | gcd_G(z_2,z_3),
K_U K_V | Norm(gcd_G(z_2,z_3)).                         (3)
```

This holds with separate primitive normalization of each window; no
assumption that their Gaussian gcds coincide is needed. It is also
unchanged by a common rotation or common Gaussian division of all six
points, with (3) then applied to the resulting actual coordinates.

After dividing the pair by its Gaussian gcd, its two distinct equal-norm
Gaussian points have squared distance at least two. Consequently

```text
K_U K_V <= |z_3-z_2|^2/2 <= N Delta_23^2/2,            (4)
```

where `N` is the original common norm and `Delta_23` is the positive
angular gap between the shared interior points. For windows shifted by
two positions, only one common interior point remains, so (2) supplies
an integer divisor but no analogous nonzero-chord bound. For a shift of
three positions there is no common interior point.

Three consecutive windows have pairwise coprime oriented factors and one
common interior point. Their product of conjugate factors divides that
point. Four consecutive windows remain pairwise coprime but have no common
interior point. These are exact local consequences; pairwise coprimality
is not asserted for windows at arbitrary larger distances.

## Actual phase cancellation and sharp pair separation

Consider the ordered primitive six-tuple

```text
N=14365,
z=((-107,-54),(-98,-69),(-91,-78),
   (-78,-91),(-69,-98),(-54,-107)).
```

Every pair in displayed order has positive determinant and positive dot
product, so these points lie in the displayed strict order inside an arc
shorter than `pi/2`. The two adjacent windows have

```text
alpha_U=3+2i,          alpha_V=3-2i,
K_U=K_V=13,           alpha_U alpha_V=13,
gcd_G(z_2,z_3)=13 up to a unit,
|z_3-z_2|^2/2=169=K_U K_V.                             (5)
```

Thus their rational twist phases obey
`(alpha_U/conjugate(alpha_U))(alpha_V/conjugate(alpha_V))=1`.
The simultaneous descents can cancel phases exactly, and the product
divisor and primitive-pair distance bounds are sharp on actual ordered
circle points. The final angular estimate in (4) is not an equality.
This tuple's
normalized arc constant is approximately `6.963`; it is not a fixture
in the `C=1/2` class. The checker certifies that exclusion with the exact
endpoint chord, without relying on the decimal approximation.

## Quantitative comparison with the existing triangle bounds

Put `R=sqrt(N)`, and suppose the entire six-tuple is on an arc of physical
length at most `L=C sqrt(R)`. The individual gap bounds give

```text
K_U <= |D_123|/2,          K_V <= |D_234|/2.
```

Writing `c=|z_3-z_2|`, the circumradius formula therefore yields

```text
K_U K_V <= |D_123 D_234|/4
 = c^2 |z_2-z_1| |z_3-z_1| |z_4-z_2| |z_4-z_3|/(16R^2)
 <= (C^4/16)c^2.                                      (6)
```

At `C=1/2`, this is `c^2/256`, stronger than the size bound `c^2/2`
in (4). Thus merely multiplying or averaging the new pair-size inequalities
does not improve what the individual triangle inequalities already imply
on that arc class. Equation (3) still retains an actual Gaussian divisor
and its orientation, which the real size comparison does not exhaust.
No claim is made that every possible use of overlapping descents reduces
to these inequalities.

The [checker](check_ordered_overlapping_gap_factors.py) exhausts bounded
allocation vectors for all three shifts, including the separate-subset
normalization and shared-interior divisibility, and verifies (5) with
exact Gaussian arithmetic.
