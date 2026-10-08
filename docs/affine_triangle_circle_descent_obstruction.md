# Affine triangle descent on a lattice circle is only content removal

## Status

A triangle can be used to define many rational affine changes of
coordinates.  This note shows that none of them gives a new descent for
five or more points on a centered lattice circle.  If an invertible rational
affine map sends those points to another lattice circle with integral center,
then the map is a similarity and its denominator divides the common Gaussian
content of all the source points.  After primitive normalization the output
radius cannot be smaller.

The argument is exact and has no height hypothesis.  It does not apply to a
nonaffine fractional-linear transformation of the circle, and it does not
prove a radius-independent point bound.

## 1. Five points force an affine map to be a similarity

Let

```text
z_1,...,z_m in Z[i],       |z_j|=R,       m>=5,
```

be distinct.  Regard a rational affine map as

```text
T(x)=Ax+b,       A in GL_2(Q),       b in Q^2.
```

Suppose the `T(z_j)` are lattice points on a circle of radius `R'` centered
at `c in Z^2`.  Then

```text
b=c,       A^t A=lambda I,       (R')^2=lambda R^2       (1)
```

for some positive rational `lambda`.

Indeed, the source circle and the pullback of the target circle have
equations

```text
G(x)=x^t x-R^2=0,
F(x)=(Ax+b-c)^t(Ax+b-c)-(R')^2=0.                       (2)
```

Both are nonsingular irreducible conics.  If they had no common component,
Bezout's theorem would give at most four intersection points over the
algebraic closure, counted with multiplicity.  The five distinct source
points therefore force

```text
F=lambda G.                                             (3)
```

Comparison of quadratic terms gives `A^t A=lambda I`.  Comparison of linear
terms gives `A^t(b-c)=0`, hence `b=c` because `A` is invertible.  The constant
terms then give the radius identity in (1).

Every rational matrix satisfying (1) is multiplication by some
`alpha in Q(i)`, possibly followed by complex conjugation.  Thus the affine
map has introduced no anisotropic triangle normalization: on the entire
source circle it is a Euclidean similarity.

The threshold five is the natural one for this argument.  Two unrelated
conics can meet in four points, so a triangle, or even four selected points,
does not by itself force the pullback conic to equal the source circle.

## 2. The denominator is exactly common Gaussian content

Let

```text
G_0=gcd_G(z_1,...,z_m),                                 (4)
```

defined up to a Gaussian unit.  In the orientation-preserving case write the
similarity in lowest Gaussian terms as

```text
alpha=u/v,       u,v in Z[i],       gcd_G(u,v)=1.       (5)
```

Since `b=c` is integral, target integrality says

```text
alpha z_j=T(z_j)-c in Z[i]                              (6)
```

for every `j`.  Euclid's lemma in `Z[i]` gives `v|z_j` for every `j`, and
hence `v|G_0`.  The orientation-reversing case gives the same conclusion
with `bar(G_0)`.  Therefore

```text
R'/R=|alpha|=|u|/|v| >= 1/|G_0|.                       (7)
```

Equality is attained by dividing every point by `G_0`.  Consequently the
least radius obtainable by a rational affine map of this kind is exactly

```text
R/|G_0|.                                                (8)
```

In particular, once the whole configuration has been divided by its common
Gaussian gcd, every centered integral affine image has radius at least `R`.

## 3. Why a common difference divisor does not suffice

Suppose instead that all selected points satisfy

```text
z_j=z_0+d h_j,       h_j in Z[i],       d in Z_(>0).    (9)
```

The apparent homothety sends them to the integral points `h_j`, but their
circle has radius `R/d` and center

```text
-z_0/d.                                                 (10)
```

This center is a lattice point exactly when `d|z_0`.  In that event (9)
shows that `d` divides every `z_j`, so the operation is precisely common
content removal.  If the center in (10) is not integral, translating it to
the origin changes `h_j` to `z_j/d`, which are not all lattice points.

Thus a large gcd of chord coordinates, even one shared by every difference,
does not produce a smaller instance of the original centered lattice-circle
problem unless it already comes from common Gaussian content.

The same proof works for a Gaussian difference divisor in place of the
rational integer `d`.

## 4. Adjacent triangle areas constrain chord content in the opposite direction

There is also a direct reason not to expect density to force large contents
of consecutive chords.  Order the points along a minor arc.  Let `s_j` be
the arc length of the gap from `z_j` to `z_(j+1)`, let

```text
e_j=z_(j+1)-z_j,       g_j=gcd(|Re(e_j)|,|Im(e_j)|),    (11)
```

and let the total arc length be `L`.  The doubled area of the consecutive
triangle is the positive integer

```text
D_j=det(e_j,e_(j+1)).                                  (12)
```

It is divisible by `g_j g_(j+1)`.  Exact circle geometry and
`sin x<=x` give

```text
D_j
 =4R^2 sin(s_j/(2R)) sin(s_(j+1)/(2R))
        sin((s_j+s_(j+1))/(2R))
 <=s_j s_(j+1)(s_j+s_(j+1))/(2R).                      (13)
```

Using `s_j s_(j+1)<=(s_j+s_(j+1))^2/4` yields

```text
s_j+s_(j+1)
 >=2 R^(1/3)(g_j g_(j+1))^(1/3).                       (14)
```

Summing (14), and observing that the sum of consecutive pair sums is at
most `2L`, proves the weighted inequality

```text
sum_j (g_j g_(j+1))^(1/3) <= L/R^(1/3).                (15)
```

For `L<=C sqrt(R)`, its unweighted consequence is the classical
`M=O_C(R^(1/6))` estimate.  The added weights show that large adjacent
chord contents consume more of the available arc length.  Hence a putative
dense endpoint cluster is pushed toward the primitive-chord regime rather
than toward an affine gcd descent.

This does not control every nonconsecutive chord and supplies no improved
growth exponent.  Its role here is to identify the direction of the exact
triangle-area constraint, not to claim a uniformity result.
