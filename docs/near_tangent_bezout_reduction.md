# Bezout complement audit for a primitive tangent vector

This is an exact congruence refinement of the earlier
[projection audit](fresh_geometric_projection_route.md). No new
circle-point count follows here.

Let `u=(p,q)` be primitive, `A=p^2+q^2`, and choose `w=(r,s)` with
`ps-qr=1`. For `z=(X,Y)`, define

```text
x=pX+qY,       y=pY-qX,       B=pr+qs.
X=(px-qy)/A,  Y=(qx+py)/A.                          (1)
```

Writing `z=m u+n w` proves that the image of `Z^2` is exactly

```text
x=A m+B n,  y=n,       or equivalently x=B y mod A. (2)
```

The identity `A|w|^2-B^2=1` gives `B^2=-1 mod A`. On the circle,

```text
x^2+y^2=A R^2.                                     (3)
```

Replacing `w` by `w+k u` changes `B` by `kA`, so one may choose
`|B|<=A/2`. This does not change the image lattice. In complex
coordinates it is simply

```text
x+i y=(p-i q)(X+i Y).                              (4)
```

Thus the congruence lattice is a rotated square lattice of spacing
`sqrt(A)`, not an arbitrary index-`A` lattice.

## The correct rectangle estimate

For an axis-aligned rectangle of side widths `K,H`, its image lattice
point count satisfies

```text
N(K,H) <= KH/A + sqrt(2)(K+H)/sqrt(A)+2.             (5)
```

Center a fundamental square of area `A` at each counted point.
Their interiors are disjoint and every square is within distance
`sqrt(A/2)` of the rectangle. Their union therefore lies inside the
rectangle of side widths `K+sqrt(2A),H+sqrt(2A)`; area comparison
proves (5).

An earlier draft incorrectly replaced the denominator `sqrt(A)` in
the boundary error by `A`. That estimate is false even with
`B^2=-1 mod A`. For an exact counterexample, take

```text
u=(p,1), w=(-1,0), A=p^2+1, B=-p,
H=floor(p/2), K=pH,
rectangle: -pH<=x<=0, 0<=y<=H.
```

Writing `x=-py+A m`, these ranges force `m=0`: `m>=1` makes
`x>0` and `m<=-1` makes `x<-pH`. Thus the count is exactly
`H+1`, while area divided by `A` is `pH^2/(p^2+1)`. Their
difference is `p/4+O(1)`, whereas `K/A+H/A+1=O(1)`. This example
concerns rectangle lattice points; they are not claimed to lie on a
common circle.

## Scope

For a near-tangent vector, `y` is the narrow coordinate and `x`
varies on the scale `O(|u|L)`. The rectangle estimate retains a
boundary term of order `L`, along with its area contribution. It
supplies no uniform endpoint bound. The exact equation (3) could
contain additional information; this coordinate change has not
extracted it.

A finite enumeration could list integer `y`, test whether
`A R^2-y^2` is an integer square, impose (2), and check arc membership
using (1). That is equivalent to the original lattice-point test,
not a proof of a new count.

The [checker](check_near_tangent_bezout_reduction.py) verifies inverse
coordinates, congruences, circle norms and the exact rectangle family.
The upper bound (5) is proved by the tiling argument above.
