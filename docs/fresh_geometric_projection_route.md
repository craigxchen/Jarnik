# Fresh geometric route: integer projections and their obstruction

This note starts from the circle itself, without prime-pattern hypotheses
or a prescribed residue inequality. It does not prove uniformity or a
better growth rate. It records a complete elementary projection test and
an explicit family that prevents a naive bounded-direction argument.

## 1. An integer-projection certificate

Let an arc on the radius-`R` circle have length `L<=C sqrt(R)` and
angular midpoint with outward unit normal `n`. Assume its half-angle
`h=L/(2R)` is at most `pi`; larger arcs can be subdivided.
Choose a nonzero integer vector `u`, write `Q=|u|`, and let `epsilon`
be its angle from `n`. The integer-valued projection of lattice points
onto `u` is the restriction of

```
F(t)=QR cos(t-epsilon),       -h<=t<=h.
```

Using the cosine addition formula, its range has width at most

```
QR (|cos epsilon|(1-cos h)+2|sin epsilon|h)
 <= Q L^2/(8R)+Q L |sin epsilon|
 <= Q C^2/8+Q C sqrt(R)|sin epsilon|.                (1)
```

The use of `2h` rather than `2 sin h` is valid on the whole stated
half-angle range. There are at most `floor(width)+1` integer levels,
and each level line intersects the circle in at most two points.
Consequently

```
M <= 2(floor(Q C^2/8+Q C sqrt(R)|sin epsilon|)+1).    (2)
```

For a perfectly aligned rational normal this is independent of the
radius. For a general orientation, (2) is a uniform certificate only
if both the denominator cost `Q` and the angular-error cost remain
bounded. Rational approximation must account for both costs.

In particular, a bound `Q C^2/8+Q C sqrt(R)|sin epsilon|<=B`
would require `Q<=8B/C^2`. There are only finitely many such integer
directions, and the arc normal would have to approach one of them
within an angle of order `R^(-1/2)`. This is a substantive arithmetic
restriction on a cluster, not a property of arbitrary thin arcs.

## 2. Explicit small clusters with an irrational limiting direction

Put `A_j(t)=t+j+i`, for `j=0,1,2`, and

```
z_j(t)=A_j(t) product_(l!=j) conjugate(A_l(t)).
```

All three are Gaussian integers of common radius

```
R_t=product_(j=0..2) |A_j(t)| asymptotic to t^3.
```

Their arguments differ by twice the differences of `arg A_j(t)`.
Since the derivative of `arctan(1/x)` has magnitude `1/(x^2+1)`,
their angular span is at most `4/t^2` for `t>0`. Thus they occupy
an arc of length `O(t)=O(R_t^(1/3))`, and are distinct.
Their common limiting direction is the positive real axis.

Now take `t=m^2` and multiply all three by the Gaussian integer

```
g_m=m+i floor(sqrt(2)m).
```

The new radius is of order `m^7`, the containing arc length is
`O(m^3)`, and hence its normalized length is

```
L_new/sqrt(R_new)=O(m^(-1/2)) -> 0.                 (3)
```

The limiting arc direction is `arg(1+i sqrt(2))`, which is not
parallel to any nonzero integer vector. For every fixed bound on
`Q`, the angular distance from all the finitely many permitted
integer directions stays bounded away from zero (up to reversing
the direction). Thus (1), with any fixed positive `C`, cannot give
a uniform certificate on this family using bounded `Q`.

This is not a counterexample to the desired theorem: it has only
three points. It does show that increasingly short normalized arcs
need not approach a rational direction of bounded denominator.
A proof along this route would have to exploit genuinely large
cardinality to force additional structure, not shortness alone.

## 3. Any hypothetical unbounded family can approach any chosen direction

There is a stronger obstruction to using limiting orientation alone.
Suppose an unbounded family exists: arcs of constant `C` contain
`M_n` points, with `M_n -> infinity`. Fix any target angle `theta`.
Put `k_n=floor(M_n^(1/3))` and select a consecutive `k_n`-point
window of arc length at most

```
C (k_n-1)/(M_n-k_n+1) sqrt(R_n).
```

This follows by summing consecutive-window lengths and counting each
gap at most `k_n-1` times. Let `alpha_n` be the window's angular
midpoint. Choose a Gaussian integer `g_n` by rounding both coordinates
of `k_n exp(i(theta-alpha_n))` to the nearest integers. Then

```
|g_n| = k_n+O(1),
arg(g_n)+alpha_n = theta+O(1/k_n) modulo 2pi.
```

Multiplication by `g_n` preserves integrality, equal radius, and
distinctness. The new normalized arc length is at most

```
C (k_n-1) sqrt(|g_n|)/(M_n-k_n+1)
 = O_C(M_n^(-1/2)) -> 0.                            (4)
```

Its cardinality `k_n` still tends to infinity and its midpoint tends
to the prescribed angle. Thus, if uniformity fails at all, it fails
along clusters with arbitrarily small normalized lengths and *any*
chosen limiting direction, including an irrational one.

This does not produce a counterexample: it is a conditional transfer
of a hypothetical counterexample. It means that limiting orientation
alone cannot single out the surviving configurations. The speed of
angular approximation, the radius, and the multiplier cost matter.
Also, the construction introduces a common Gaussian factor. It makes
no claim about the direction after primitive normalization.

## 4. Why a modular-box estimate does not fill the gap

A second possibility was to retain large denominators and count in
the resulting congruence boxes. The primary paper of Cilleruelo and
Garaev, [Concentration points on two and three dimensional modular
hyperbolas](https://arxiv.org/abs/1007.1526), gives, in particular, an
`M^(o(1))` bound for a modular hyperbola in a square of side
`M<p^(1/4)`. This still grows with the interval size and does not
provide a constant. In addition, a reduction of the circle problem
would have to preserve that theorem's prime-modulus and scale
hypotheses. Neither a uniform theorem nor such a complete reduction
has been inferred here.

## 5. Remaining issue for this approach

The elementary projection test identifies an orientation denominator
as the obstacle. Simply selecting a good rational approximation,
rotating the picture, or replacing the arc by a bounded-width strip
does not remove it. A genuinely new step would need to use the exact
quadratic equation and simultaneous integrality of both coordinates
to control the number of points while permitting this denominator
to grow. No such step has been proved in this note.
