# Circle symmetry with the integrality domain retained

This route uses the actual reflection symmetries of a circle. It does
not assume a prime-pattern limiting model. The conclusions below are
prose proofs; the uniform endpoint count remains unproved.

## 1. Every pair determines a rational reflection

Let `a,b` be nonantipodal Gaussian integers with `|a|=|b|=R`, and
let `N=R^2`. Reflection about the radius through the midpoint of
their short connecting arc is

```
sigma_ab(z) = (ab/N) conjugate(z).                   (1)
```

It exchanges `a,b`, is an involution, and preserves the entire circle.
Its coefficient is rational Gaussian and has modulus one. These facts
do not mean it preserves the integer lattice.

To find its exact integrality domain, take the primitive integer vector
`u=A+iB` parallel to `a+b`, with `gcd(|A|,|B|)=1`. Set

```
d = u                  if A,B have opposite parity,
d = u/(1+i)            if A,B are both odd.
```

In the first case put `eta=1`, in the second `eta=i`. Then

```
ab/N = u/conjugate(u) = eta d/conjugate(d),
gcd(d,conjugate(d)) = 1  up to units.
```

The coprimality follows from primitivity: an odd Gaussian prime common
to a number and its conjugate would force a rational prime to divide
both coordinates; the only remaining prime, `1+i`, has been removed.
Consequently

```
sigma_ab(z) belongs to Z[i]  iff  d divides z.        (2)
```

Indeed integrality is equivalent to `conjugate(d)` dividing
`d conjugate(z)`, and coprimality permits cancellation. The domain
is the sublattice `d Z[i]`, whose index is

```
D=|d|^2 = (A^2+B^2)/kappa,    kappa=1 or 2.        (3)
```

Both anchors belong to this domain, so `D` divides their norm `N`.
Globally lattice-preserving reflections arise only when `D=1`;
they are the familiar reflections about the axes and diagonals.

## 2. Integral reflection of a third nearby point is impossible

Suppose a circle arc has angular width `Delta<pi/2` and
`Delta sqrt(R)<sqrt(2)`. Choose distinct points `a,b` on it.
Then at most those two arc points belong to the domain in (2).

For if a third point `c` belongs to the domain, put `c'=sigma_ab(c)`.
All four are Gaussian, have radius `R`, and satisfy

```
a b = c c'.
```

Choose argument lifts of `a,b,c` in the original interval of length
`Delta`. The reflected argument is `arg a+arg b-arg c`. The original
interval together with its reflection has width at most `2Delta`.
The two product pairs are different, since `c` is different from both
anchors, even if `c'=c`. The sharp multiplicative-rectangle separation
proved in `multiplicative_rectangle_separation.md` therefore gives

```
2Delta sqrt(R) >= 2sqrt(2),
```

a contradiction. Thus, on arcs of constant `C<sqrt(2)`, every
pair-bisector reflection domain contains exactly its two anchors among
the selected points. At `C=1/2` the angular hypothesis holds whenever
the arc contains a nonzero lattice point.

If `a,b` are the outermost points of the arc, the reflected interval
is the original interval itself. The same argument then works with
the stronger threshold `Delta sqrt(R)<2sqrt(2)` and `Delta<pi`.

This is an integrality-domain restriction, not a cardinality bound.
It does not justify repeatedly reflecting arbitrary members of a cluster.

## 3. Explicit denominator growth on increasingly short arcs

The three-point family in `fresh_geometric_projection_route.md` makes
the denominator issue concrete. For an odd positive integer `t`, put

```
A_j=t+j+i,    j=0,1,2,
z_j=A_j product_(l!=j) conjugate(A_l).
```

The points have radius of order `t^3`, and angular span at most `4/t^2`.
Their normalized arc length therefore tends to zero. The sum of the
outer points is

```
z_0+z_2 = 2(t+1)^2 conjugate(A_1).
```

The primitive bisector vector is `u=(t+1)-i`. Its coordinates have
opposite parity, so the reflection domain is precisely

```
((t+1)-i) Z[i],       D=(t+1)^2+1.                  (4)
```

The third point `z_1=A_1 conjugate(A_0) conjugate(A_2)` is coprime
to `u`. The differences between `u` and the latter two conjugate
factors are units. A common divisor of `u` and `A_1` would divide
`2i`, whereas `|u|^2` is odd. Thus `sigma_(z_0,z_2)(z_1)` is not
integral. The circle symmetry exists, but cannot supply another
lattice point in this example.

This family has only three points and does not contradict uniformity.
It does prevent treating the reflection denominator as bounded merely
because the arc is increasingly short relative to `sqrt(R)`.

## 4. What a symmetry proof would still have to establish

A possible sufficient result would force, in every sufficiently large
short-arc cluster, a third point into the integrality domain of one
pair-bisector reflection. Section 2 would immediately contradict it.
No such forcing result is proved here; it cannot be assumed from
circle symmetry alone.

In the squarefree Gaussian-factor model, domain membership means that
the third row agrees with the two anchors at every factor where the
anchors agree. Abstract binary families can avoid this for arbitrarily
many rows. Thus a proof must exploit the actual Gaussian values or the
actual angles, not only the incidence of reflection domains.

The next mathematical issue is therefore specific: can the simultaneous
avoidance of all these rational-reflection sublattices persist for an
unbounded cluster on an endpoint arc? The geometric identities and
domain computation above do not resolve it.
