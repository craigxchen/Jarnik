# Fixed pentagon projective radius floor

## Status

For a fixed ordered positive rational pentagon, the simultaneous primitive
Ptolemy content `F` is projectively invariant.  It bounds the endpoint-versus-
interior gap of every primitive integral circle realization, including
realizations obtained from one another by arbitrary rational projective
changes of half-angle coordinate.  It also gives a scale-retaining lower
bound for the exact least squared radius in the projective orbit:

```text
K <= F <= 2 N_*.
```

Here `K` belongs to any chosen primitive realization and `N_*` is the minimum
of the primitive squared radius over the whole rational projective orbit.
Thus projective minimization would imply the desired small gap only after an
independent upper bound for `N_*`.  The arithmetic furnished by the fixed
shape points in the other direction.  If a radius-minimizing image is also an
endpoint image with normalized arc constant `C_*`, the stronger conditional
bound is `F <= (C_*^2/4) sqrt(N_*)`.

This is a bounded no-gain audit, not progress on the growth exponent.  The
note gives the exact minimization formula and the short proof.  It does
not use the universal cross-ratio height norm or a CRT obstruction.

## 1. Exact radius functional on the rational projective orbit

Represent the five distinct rational half-angle parameters by primitive
integer columns

```text
H_i=(x_i,y_i),                 0<=i<=4.
```

Every rational projective change of the half-angle coordinate is represented
by a nonsingular primitive
integral matrix `M`; multiplying `M` by a nonzero rational scalar changes
nothing after the following row normalization.  Put

```text
V_i=M H_i,
c_i=gcd(|Re V_i|,|Im V_i|),
B_i=V_i/c_i.
```

Regard `B_i` as a Gaussian integer.  A general `M` moves the anchor, so one
must first take relative phases.  For `i>0` put

```text
C_i=B_i bar(B_0),
G_i=gcd_G(C_i,bar(C_i)),
D_i=C_i/G_i,
D_0=1.                                                   (1)
```

The Gaussian gcd is chosen compatibly with conjugation; changing its unit
changes `D_i` only by a unit.  Thus `D_i` and `bar(D_i)` are coprime and

```text
(B_i/bar(B_i))/(B_0/bar(B_0))=epsilon_i D_i/bar(D_i)
```

for a Gaussian unit `epsilon_i`.  The unit is harmless: it has a unit or
`1+i` half-angle representative and introduces no odd Gaussian denominator.
Formula (1), in particular its gcd, also handles the odd-odd parity case.
Therefore the exact least primitive squared radius realizing these labelled
relative phases is

```text
N(M)=Norm(lcm_G(D_0,D_1,D_2,D_3,D_4)).             (2)
```

Indeed an integral anchor realizing all relative phases must be divisible by
every `bar(D_i)`, while the conjugate of their Gaussian lcm realizes them.
This realization is primitive: a common nonunit Gaussian divisor would
produce another integral anchor divisible by all `bar(D_i)` with smaller
norm, contradicting the lcm's divisibility property. Units in the phases
do not affect this argument.
More explicitly, if `L=lcm_G(D_0,...,D_4)`, every integral realization
`w_i` of these relative phases has `w_0=alpha bar(L)` with
`alpha in Z[i]`. Its entire tuple is `alpha` times the displayed minimal
realization. Thus every primitive realization of these relative phases
has `alpha` a unit and squared radius exactly (2).
Taking ratios against the transformed anchor is a rational unit-circle
rotation. Clearing their Gaussian denominators is then a common similarity
that preserves all relative phases and cross-ratios. Its rotation after
normalizing to the unit circle need not have rational coordinates; no such
claim is needed for the radius of an integral realization.
Consequently the orbit minimum is exactly

```text
N_*=min_(M in Mat_2(Z), gcd(entries M)=1, det M!=0) N(M). (3)
```

The minimum is attained simply because the nonempty set of values in (3) is
a set of positive integers.

Formula (3) is exact, but it is a discrete Gaussian-lcm minimization rather
than a reduction of a product of binary quadratic forms.

### Why relative normalization is necessary

Let `pi=2+i` and take the five primitive Gaussian rows

```text
H_i=pi^(i+1),                 0<=i<=4,
M=I.
```

Their arguments are strictly cyclic.  The incorrect lcm of the five absolute
row denominators is `pi^5`, of norm `5^5=3125`; it has silently inserted an
extra phase-one anchor.  Relative to the first row, however, the phases are
`(pi/bar(pi))^i`.  Formula (1) gives `D_i=pi^i`, so (2) is `5^4=625`.
It is realized primitively by

```text
z_i=pi^i bar(pi)^(4-i),       0<=i<=4,
```

whose Gaussian gcd is one.  This example also shows why merely applying the
odd-odd division to every absolute transformed row does not repair the
anchor error.

## 2. The invariant content survives minimization

Write the first two positive pentagon coordinates in lowest terms as

```text
u_0=a/b,             u_1=c/d,
H=gcd(b,d),           X=gcd(a,d-c),
Y=gcd(c,b-a),         F=(ad+bc-bd)/(HXY).
```

The five `u_i` are the positive absolute values of labelled real
cross-ratios, so `F` depends only on the labelled projective shape and is
unchanged by every `M in PGL_2(Q)`.  A real projective automorphism of the
circle preserves cyclic order or reverses it globally, even if its affine
formula crosses infinity.  With the row labels retained, each chord matching
still has the same positive absolute cross-ratio; reversal therefore does
not alter the displayed `u_i` or `F`.  The exact primitive-content
identity in `positive_pentagon_primitive_content.md` proves, for every
primitive integral realization of this shape,

```text
K | F.                                                   (4)
```

This remains true when `K` and `F` have accidental equal-level chord
content: primewise, at every gap prime,

```text
v_p(F)=h_p+v_p(c_O)+v_p(c_I).
```

Thus applying a projective transformation does not erase the source gap
from the invariant; it places it inside the same fixed integer `F`.

## 3. A universal radius floor

Consider any primitive circle realization of squared radius `N'=R'^2`.
For each of the three matching ratios

```text
q_ab=(l_04 l_ab)/(l_0b l_a4),
```

the parity/Bezout argument for primitive odd norm gives

```text
num(q_ab) <= l_04 l_ab/2.
```

Every chord is at most `2R'`.  Since `F` is the gcd of these three reduced
numerators,

```text
F <= num(q_ab) <= 2R'^2=2N'.                           (5)
```

Apply (5) to a minimizer in (3) and combine with (4):

```text
boxed:                 K <= F <= 2N_*.                 (6)
```

The factor two only uses the diameter bound.  If the minimizing realization
lies on a minor arc of angular span `Delta_*` with
`Delta_* <= C_* N_*^(-1/4)`, some pair among the three labels `1,2,3`
has angular separation at most `Delta_*/2`. The chord `l_04` is at most
`sqrt(N_*) Delta_*`. Apply the numerator bound to that interior pair to get

```text
F <= (C_*^2/4) N_*^(1/2).                              (7)
```

This argument allows the large complementary gap to occur anywhere in the
labelled cyclic order after the projective map.

No endpoint condition on the original realization implies such a condition
for an orbit minimizer: a rational projective map may enlarge angular span.

## 4. Precise obstruction to a gap exponent

For an original endpoint realization of squared radius `N`, (6) would give

```text
log K <= delta log N+O(1)
```

only after proving the separate orbit-compression statement

```text
N_* <= O(N^delta).                                      (8)
```

Nothing in the fixed-shape divisibility proves (8).  Instead (5) gives the
opposite-direction necessary condition `N_*>=F/2`.  Even preservation of the
endpoint condition by the minimizer, with a uniformly bounded `C_*`, would
merely replace the needed exponent in (8) by `2 delta` through (7).
For a fixed target `delta<1/15`, the sufficient compression statements
remain `N_*<=O(N^(2 delta))` for such an endpoint minimizer, or
`N_*<=O(N^delta)` without that condition. A bound whose exponent only
approaches `2/15` or `1/15` from below does not supply a fixed exponent
margin. An uncontrolled `C_*` likewise supplies no such conclusion.

Accordingly, optimizing projective symmetry does retain the scale through
the exact integer `F` and the Gaussian-lcm radius (2), but it does not by
itself improve the existing `K<=C^3N^(1/4)/16`.  The missing theorem is an
upper bound for the orbit minimum, together with endpoint control if one
wants to use the square-root estimate (7).

## Verification

Root corrected an initial absolute-phase lcm formula using the explicit
five-row example above, and removed an unnecessary claim that the rotation
part of denominator clearing has rational coordinates. The argument uses
relative phases and a common similarity.

The [exact checker](check_fixed_pentagon_projective_radius.py) compares
three computations on 448 nonsingular primitive integer projective matrices:
the relative-phase Gaussian lcm, an independently cleared absolute tuple
followed by its full Gaussian gcd, and the ordinary lcm of all pair
denominator norms. It verifies their agreement, both orientations of the
projective maps, invariant positive pentagon coordinates and `F`, and the
explicit `3125` versus `625` anchor counterexample. These finite checks
supplement the formula's proof; they do not solve the minimization or
supply a uniform circle-point bound.
