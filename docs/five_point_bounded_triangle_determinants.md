# Five points rule out a bounded adjacent-determinant alphabet

There is an exact rigidity statement behind the adjacent-chord approach:
five lattice points on a quarter-circle arc cannot have all their
consecutive triangle determinants bounded while their radius tends to
infinity. The proof uses lattice determinants, chord order, and the
uniqueness of a conic through five points. It does not give a uniform
endpoint count, because the determinants themselves may grow.

An explicit integer-parabola family below shows why exact concyclicity is
essential. It has unbounded cardinality, bounded chord recurrence and
determinants, endpoint-normalized angular span tending to zero, and radial
error tending to zero. That error nevertheless exceeds the exact
squared-norm quantization scale.

## 1. Quantitative bounded-alphabet theorem

Let `P_0,...,P_4` be five distinct lattice points, in order on a circle arc
of angular width at most `pi/2`. Put

```text
e_i=P_(i+1)-P_i,       i=0,1,2,3,
D_i=|det(e_(i-1),e_i)|, i=1,2,3.
```

Suppose `D_i<=K`, where `K>=1`. Define

```text
B=2^18 K^10,
H=24(16B^2)^4.
```

Then the circle radius satisfies

```text
R^2 <= 8H^3.                                      (1)
```

The constants are intentionally crude. In particular (1) is a bound
depending only on the determinant alphabet, independent of the circle,
orientation, chord contents, or recurrence coefficients.

### Chord lengths become comparable

Write the four chord lengths as `l_0,...,l_3`, with minimum `m` and maximum
`L`. Consecutive chords of lengths `a,b` and angular gaps `alpha,beta`
satisfy

```text
D(a,b)=ab sin((alpha+beta)/2),
ab(a+b)/(2sqrt(2)R) <= D(a,b) <= ab(a+b)/(2R).      (2)
```

Indeed the half-gaps have sum at most `pi/4`, and
`sin(x+y)=sin x cos y+cos x sin y`, with both cosines at least `1/sqrt(2)`.

Each `D_i` is a nonzero integer, hence lies in `[1,K]`. Comparing two
successive formulas (2), whose middle length is shared, gives

```text
1/(sqrt(2)K) <= l_0/l_2 <= sqrt(2)K,
1/(sqrt(2)K) <= l_1/l_3 <= sqrt(2)K.               (3)
```

For example, `a(a+b)/(c(c+b))>=a/c` when `a>=c`.

Put `T=sqrt(2)K`. If the largest and smallest lengths have the same
parity of index, (3) already bounds `L/m` by `T`. Otherwise both chords
of the smaller parity have length at most `Tm`. They are nonparallel:
their chord-midpoint directions are strictly ordered within an angle
less than `pi`. Their integer determinant therefore has modulus at
least one.

The entire angular span `Delta` satisfies

```text
Delta <= sqrt(2)(l_0+l_1+l_2+l_3)/R <= 4sqrt(2)L/R.
```

The adjacent determinant involving the longest chord, together with (2),
gives `1/R<=2sqrt(2)K/(L^2 m)`. Therefore the determinant of the two
short-parity chords is at most

```text
(Tm)^2 Delta <= 32K^3 m/L.
```

Its nonzero integrality proves

```text
L/m <= 32K^3.                                    (4)
```

This is the step that rules out the unbounded alternating long/short
pattern available to four points: five points supply two distinct short
chords whose determinant must also be a nonzero integer.

### All chord determinants and affine coordinates are bounded

Using (2) once more gives `1/R<=sqrt(2)K/m^3`. Thus any two of the four
chords obey

```text
|det(e_i,e_j)| <= L^2 Delta
 <= 8K(L/m)^3 <= 2^18 K^10 = B.                  (5)
```

Reverse orientation if necessary so `d=det(e_0,e_1)>0`. Then `d<=K`.
With `A` the matrix whose columns are `e_0,e_1`, set

```text
U_i=d A^(-1)(P_i-P_0).
```

The `U_i` are integer vectors, because every coordinate of `d A^(-1)e_j`
is a determinant of two integer chords. Formula (5) bounds each such
coordinate by `B`, so

```text
U_0=0,       |(U_i)_x|, |(U_i)_y| <= 4B.          (6)
```

The physical Euclidean metric in these coordinates is

```text
G=A^T A/d^2,       det G=1/d^2.                   (7)
```

### Five-point conic uniqueness fixes the metric scale

The four nonzero nodes `U_i` impose four linear equations on the five
coefficients of a conic through zero,

```text
a x^2+bxy+c y^2+u x+v y=0.                        (8)
```

Their coefficient matrix has rank four. To see this without an
irreducibility assumption, evaluation of all six quadratic monomials on
the five nodes has rank five: for each chosen node, pair the other four
nodes into two lines and multiply their equations. No three circle
points are collinear, so this quadratic vanishes at the other four and
is nonzero at the chosen node. Removing the zero-node condition leaves
rank four in (8).

Every matrix entry has modulus at most `16B^2`. Its signed `4 x 4`
cofactors give a nonzero integer kernel vector. Dividing its integer
content gives primitive coefficients of (8), each of modulus at most

```text
H=4!(16B^2)^4.
```

Choose their sign so the quadratic matrix

```text
G_0=[[a,b/2],[b/2,c]]
```

is positive definite. Uniqueness implies `G=lambda G_0`, `lambda>0`.
Since `4 det G_0=4ac-b^2` is a positive integer, (7) gives

```text
det G_0>=1/4,       lambda<=2/d<=2.
```

Completing the square in (8), the squared physical radius is

```text
R^2=(lambda/4)(u,v) G_0^(-1) (u,v)^T.
```

Every entry of `G_0^(-1)` has modulus at most `4H`, while
`|u|,|v|<=H`. Thus the displayed quadratic form is at most `16H^3`,
and in fact `R^2<=8H^3/d`; dropping `d>=1` gives (1).

## 2. Exact equality of neighboring determinants is already forbidden

For four ordered circle points in an arc of width less than `pi`, let
their angular gaps be `alpha,beta,gamma`. Equality of the two consecutive
triangle determinants gives

```text
sin(alpha/2) sin((alpha+beta)/2)
 = sin(gamma/2) sin((beta+gamma)/2).
```

For fixed `beta`, the left expression is strictly increasing in its gap
argument throughout this range: its derivative is one half of
`sin(gap+beta/2)>0`. Hence `alpha=gamma`. The four lifted angles then obey
`theta_0+theta_3=theta_1+theta_2`, so their Gaussian representatives satisfy
`z_0 z_3=z_1 z_2`.

By `multiplicative_rectangle_separation.md`, such equality is impossible
on an arc of length less than `2sqrt(2)sqrt(R)`. Conversely equality of
the outer gaps gives equal determinants. The four-point example in
`multipoint_continuation.md` reaches this threshold asymptotically.

This local observation alone would only exclude neighboring repeated
determinant values. The bounded-alphabet theorem in Section 1 is stronger:
even an alternating bounded alphabet cannot persist at unbounded radius
through five consecutive lattice points on a quarter-circle arc.

## 3. Arbitrarily close circles do not suffice

Let `Q>=3` be odd and define

```text
C=((Q^2+1)/2, -Q(Q^2+1)/2),
R=(Q^2+1)^(3/2)/2,
P_j=(Qj+j^2,j)-C,       0<=j<=M<=Q.
```

All `P_j` are integer points; `C` is integral and
`R^2=(Q^2+1)^3/4` is an integer. Direct expansion gives

```text
|P_j|^2-R^2=j^3(2Q+j).                            (9)
```

Consequently they lie in the annulus

```text
R <= |P_j| <= R+3M^3/Q^2.                         (10)
```

Their argument is strictly monotone in `j`. For the interpolating curve
`P(t)=(Qt+t^2,t)-C`,

```text
det(P(t),P'(t))
 = -t^2-(Q^2+1)^2/2-Q(Q^2+1)t < 0.
```

Also `|P_M-P_0|<=3QM`. Since all radii are at least `R`, their angular
span obeys, for the stated range,

```text
Delta <= 2|P_M-P_0|/R,
R Delta/sqrt(R) <= 6sqrt(2)M/sqrt(Q).              (11)
```

There is no winding ambiguity: `P_0 dot P_j=R^2-(Q^2+1)j^2/2>0` for
`j<=Q`, so all these arguments lie in a single interval of width less
than `pi/2`. The chord inequality gives
`Delta<=2arcsin(|P_M-P_0|/(2R))<=2|P_M-P_0|/R`.

The consecutive chords are primitive and satisfy

```text
e_j=(Q+2j+1,1),
|det(e_j,e_(j+1))|=2,
e_(j+1)=2e_j-e_(j-1).                             (12)
```

Taking `M=floor(Q^(1/4))` gives unbounded cardinality, endpoint-normalized
angular span tending to zero, and radial error at most `3Q^(-5/4)` tending
to zero. Thus bounded determinants, bounded recurrence coefficients,
primitive chord directions, and even a shrinking annular error cannot
replace exact concyclicity in Section 1.

The scale distinction is material. Since squared norms are integers,
radial uncertainty smaller than order `1/R` can force an exact common
norm. Here `1/R` has order `Q^(-3)`, much smaller than the error in (10).
More exactly, the next possible larger squared norm is separated by
`sqrt(R^2+1)-R=1/(sqrt(R^2+1)+R)` in radius.
Equation (9) retains the actual nonzero integer norm defects. This example
only obstructs arguments stable at its error scale; it does not obstruct
an argument that uses exact norm quantization.

## 4. Remaining scope

The theorem proves that any unbounded-radius five-point sequence on such
arcs must have an unbounded adjacent determinant. Equivalently, it rules
out all fixed bounded determinant alphabets, not just a fixed recurrence
coefficient. The displayed constants even give a weak explicit
radius-dependent lower bound for the largest determinant.

It does not bound these determinants on endpoint arcs: the elementary
chord estimate still permits them to grow with `R`. A uniform point-count
proof would need a many-point mechanism producing a five-point window
with a radius-independent determinant bound, or different information
that handles the growing alphabet. No such extraction is asserted here.

The standard-library checker `check_five_point_bounded_triangle_determinants.py`
verifies 56 exact five-point conic reconstructions, including the metric
determinant `1/d^2`, its primitive scale, and the original radius. It also
checks the integer norm defects and chord recurrence of the annular family.
These finite checks supplement the proofs above.
