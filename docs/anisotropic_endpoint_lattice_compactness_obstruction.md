# The centered endpoint lattice has a collapsing anisotropic covolume

There is an exact integer coordinate model at an endpoint of a lattice
circle.  It keeps the circle center, the content of the endpoint, and the
congruence lattice simultaneously.  The model gives the expected tangent
scale `\sqrt{\rho}` (where `\rho` is the circle radius), but it also
shows why convergence of the normalized real points alone does not retain
their multiplicities. Any compactness proof must keep more arithmetic
information than this normalization by itself provides.

## 1. Exact dot/cross coordinates

Let `z_0=(a,b)\in\mathbb Z^2` be a nonzero point on the centered circle
`|z|=\rho`, so `N=\rho^2=a^2+b^2\in\mathbb Z`.  Write

```text
g=gcd(a,b),       z_0=g u,       u=(p,q),       A=p^2+q^2.
```

Thus `\rho=g\sqrt A`, and `u` is primitive.  For another lattice point
`z=z_0+\delta` on the same circle, put

```text
x=p delta_x+q delta_y,
y=p delta_y-q delta_x,
d=-x.
```

Both `x` and `y` are integers.  Expanding
`|z_0+\delta|^2=|z_0|^2` gives

```text
x^2+y^2=-2gA x,
d>0,
y^2=d(2gA-d).                                  (1)
```

The positivity of `d` is strict for a distinct point, and

```text
|delta|^2=(x^2+y^2)/A=2g d.                     (2)
```

Consequently an endpoint arc of length `L` gives the exact integer bound
`1<=d<=L^2/(2g)`.  This is the radial arithmetic content of the short arc;
no translated circle has been introduced.

The full image lattice is also explicit.  Choose integers `r,s` with
`ps-qr=1` and set `B=pr+qs`.  Every displacement has a unique form
`delta=m u+k(r,s)`, and hence

```text
y=k,       x=A m+B k,
x=B y mod A,
d+B y=0 mod A.                                  (3)
```

Moreover the Bezout identity gives
`A(r^2+s^2)-B^2=(ps-qr)^2=1`, hence `B^2=-1 mod A`.  Formula (3) is
the centered, rotated-square congruence; it is not a generic index-`A`
lattice assertion.

## 2. The correctly normalized local parabola

Use the physical orthonormal normal and tangent coordinates relative to the
actual endpoint radius:

```text
n = x/sqrt(A),
t = y/sqrt(A),
s = t/sqrt(rho).                                (4)
```

Here `n` is the signed radial displacement and `s` is the tangent
displacement divided by `sqrt(rho)`.  Since `rho=g sqrt(A)`, (1) becomes

```text
s^2=-2n-n^2/rho.                                (5)
```

This is an exact identity, including the finite-radius correction.  If
`L<=C sqrt(rho)`, then (2) and (4) give

```text
-C^2/2 <= n < 0,       |s| <= C.                (6)
```

Thus the usual tangent-by-`sqrt(rho)`, normal-by-`1` scaling is the right
one.  The remaining arithmetic is precisely (3), together with the fact
that the endpoint is the actual centered vector `g(p,q)`.

## 3. What the normalization loses

Let `Lambda_{z_0}` be the image of the displacement lattice `Z^2` under
`delta -> (n,s)`.  The determinant of the integer dot/cross map is `A`,
and the two scalings in (4) have determinant
`1/(A sqrt(rho))`.  Therefore

```text
                 covol(Lambda_{z_0})=rho^(-1/2).  (7)
```

The covolume tends to zero along every sequence of radii tending to
infinity. Hence these normalized lattices do not lie in a compact subset
of the space of full-rank lattices: covolume is a positive continuous
function on that space. Rescaling to covolume one multiplies the local
coordinates by `rho^(1/4)` and changes the fixed parabola scale.
Neither observation rules out a compactification that retains additional
arithmetic or multiplicity data.

There is an actual centered primitive obstruction to interpreting this as
a harmless lattice artifact.  In the four-point Pell family of
[four_point_bonus_counterexample.md](four_point_bonus_counterexample.md),
the four points have common norm `N=\rho^2`, with `\rho\asymp V^3`, while
all pairwise displacements are `O(V)`. For any chosen endpoint, (2) gives
`g d=|delta|^2/2=O(V^2)` and therefore, regardless of its content g,

```text
|n|=d/sqrt(A)=g d/rho=O(V^(-1)),
|s|=|delta_t|/sqrt(rho)=O(V^(-1/2)).               (8)
```

All four distinct points therefore converge to the same point `(0,0)` on
the limiting parabola.  Their collapse is compatible with the exact
congruence (3), the actual centered circle, and primitive Gaussian
content. Counting distinct limiting points would therefore underestimate
even this four-point configuration. A compactness argument must control
the number of original points represented by a single limit point;
real convergence of the normalized coordinates alone loses distinctness.

This is only an obstruction to the tempting compactness step.  It gives no
counterexample to a bound for arbitrarily many points: the Pell family has
four points, and this normalization argument gives no growing cluster.

## 4. A four-point rational-direction height floor

The same primitive Pell family also tests extraction of a short rational
direction from a physical endpoint cluster, without any auxiliary map.  Use
the exact coordinates of
[the four-point construction](four_point_bonus_counterexample.md), with
`U^2-5V^2=1`, `U/V -> sqrt(5)`, and common squared radius `N asymp V^6`.
Its displayed large vector `K` satisfies

```text
K/V^3 = (40t+32t^2) + i(60t+4t^2),       t=U/V,
z_j-K=O(V),                               0<=j<=3.
```

The Pell identity gives `t-sqrt(5)=O(V^(-2))`.  Hence, uniformly in the
four labels,

```text
arg(z_j) = theta_* + O(V^(-2)) = theta_* + O(N^(-1/3)),
theta_* = arctan(alpha),                  alpha=(sqrt(5)-1)/2.
```

Indeed, the ratio of the imaginary and real coordinates of the limiting
vector is exactly
`(60 sqrt(5)+20)/(40 sqrt(5)+160)=alpha`.  These four
points lie in arcs whose length divided by `sqrt(R)` tends to zero, so
the direction estimate concerns genuine centered endpoint tuples.

There is an elementary height floor for approximating `theta_*`.  For any
nonzero Gaussian integer `v=a+ib`, put `Q=|v|` and
`beta=-(sqrt(5)+1)/2`.  The integer

```text
f(a,b)=b^2+ab-a^2=(b-alpha*a)(b-beta*a)
```

is nonzero: for `a!=0` a rational slope cannot equal either quadratic
irrational root, while for `a=0` it equals `b^2`.  Since
`|b-beta*a|<=sqrt(1+beta^2) Q<2Q`, integrality gives
`|b-alpha*a|>=1/(2Q)`.  Taking the determinant with `1+i alpha`
therefore proves the exact angular inequality

```text
dist_mod_pi(arg(v),theta_*)
  >= |sin(arg(v)-theta_*)|
  >= 1/(2 sqrt(1+alpha^2) Q^2).
```

Consequently each physical Pell direction obeys

```text
dist_mod_pi(arg(v),arg(z_j))
  >= c/Q^2 - O(N^(-1/3)),
c=1/(2 sqrt(1+alpha^2)),                  0<=j<=3.       (9)
```

For every fixed `delta>0` and fixed angular constant `K_0>0`, (9)
excludes, at all sufficiently large Pell indices, every rational
direction of length `Q<=N^(1/8-delta)` within `K_0 N^(-1/4)` of any
physical point.  More generally, if `0<=eta<1/12` and the angular
precision is `K_0 N^(-1/4+eta)`, (9) requires
`Q>=c_(K_0) N^(1/8-eta/2)` for sufficiently large indices.  The
restriction on `eta` is a convenient conservative range; the displayed
inequality is the precise statement.

This is a four-point obstruction to a direct low-height conclusion at the
natural `N^(-1/4)` angular scale.  It does not disprove a claim of height
`N^(O(1/m))` for large `m`: four points allow the hidden constant in
`O(1/m)` to exceed `1/8`, and no unbounded endpoint clusters are
constructed here.  Any such extraction must use arithmetic information
that becomes stronger with the number of points; real coalescence and
the short-arc geometry alone do not yield it.
