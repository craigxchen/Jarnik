# Rational projective symmetry and the cost of moving cross-ratios

This note gives a new algebraic restriction on endpoint clusters. It does
**not** prove a uniform point count. In particular, its constants depend on
a fixed rational configuration; they must not be treated as uniform when
the cross-ratios change.

The arguments are prose proofs, not Lean formalizations. The extension
from unimodular matrices to general integer matrices was developed jointly
with the primary agent; the exact minimal-radius formula and the
five-point height consequence are recorded here as part of that argument.

## 1. Exact arithmetic cost of a rational circle configuration

Write a rational unit-circle direction as

```text
U/conj(U),                 U in Z[i], U != 0.
```

First divide the real and imaginary parts of `U` by their positive
integer gcd. If both resulting coordinates are odd, divide once more by
`1+i`. The resulting Gaussian integer `d` satisfies

```text
gcd(d,conj(d)) = 1,
U/conj(U) = epsilon d/conj(d),        epsilon in {1,i}.
```

For finitely many directions let `L=lcm(d_j)` and `G=gcd(d_j)`, with
Gaussian gcds and lcms understood up to units. A rational Gaussian number
`beta` makes every

```text
z_j = beta epsilon_j d_j/conj(d_j)
```

integral if and only if

```text
beta belongs to (conj(L)/G) Z[i].                      (1)
```

Consequently the least possible positive radius of a centered lattice
circle realizing these rational directions, with an arbitrary common
rational rotation and scale permitted, is exactly

```text
R_min = |L|/|G|.                                      (2)
```

To prove (1), at a Gaussian prime `pi` the required lower bound on
`v_pi(beta)` is

```text
max_j (v_pi(conj(d_j))-v_pi(d_j)).
```

The two valuations in each difference cannot both be positive. Therefore
this maximum equals

```text
max_j v_pi(conj(d_j)) - min_j v_pi(d_j)
 = v_pi(conj(L)/G).
```

This also handles the case that every `d_j` has a common prime factor:
`beta` can have a denominator, and the factor `G` in (2) is essential.
The removed units have no effect on the integrality ideals.

This formula exposes the exact lcm-versus-gcd problem. By itself it is
not a stronger general inequality than the earlier Gaussian-factor
framework. The useful extra input below is that projective images of a
*fixed* rational configuration have controlled pair gcds and controlled
products of individual coordinate contents.

## 2. A theorem for a fixed rational projective configuration

Fix `M>=2` primitive integer vectors

```text
v_j=(a_j,b_j),
Delta_ij = a_i b_j-a_j b_i != 0       for i!=j.
```

Thus they specify distinct rational projective directions. Define the
positive constants

```text
S = product_(i<j) |Delta_ij|,
d_* = min_(i<j) |Delta_ij|,
K = max(1,
        (|Delta_lj|+|Delta_il|)/|Delta_ij| : i!=j, all l).
```

Let `A` be any nonsingular rational two-by-two matrix. Clearing its
denominators gives an integer matrix. Regard its two columns as Gaussian
integers, divide both by their common Gaussian gcd, and use the resulting
integer matrix in place of `A`. Dividing by that common Gaussian number
only rotates all the doubled-angle directions by the same angle, so it
is absorbed into the common factor `beta`.

We may therefore assume the two Gaussian columns of `A` are coprime.
Put

```text
D = |det A| >= 1,
U_j = (A v_j)_1 + i (A v_j)_2,
h_j = |U_j|,
H = max_j h_j,
m = min_j h_j.
```

Suppose all `z_j=beta U_j/conj(U_j)` are Gaussian integers and lie
on an arc of angular width `Delta<pi`, on their common radius-`R`
circle. Write `C_arc=Delta sqrt(R)`. Then

```text
C_arc >= c sqrt(D) H^((M-4)/2),                       (3)

c = 2^(1-M/4) sqrt(d_*) / (S K^((M-1)/2)) > 0.
```

In particular, for a fixed rational configuration with `M>=5`, every
realization with `C_arc<=C` has

```text
H <= (C/c)^(2/(M-4)),
R <= C^2/(4 d_*^2) (C/c)^(8/(M-4)).                  (4)
```

Thus no unbounded-radius endpoint family with five or more points can
have its rational cross-ratios all fixed. This conclusion allows *all*
rational projective transformations, not just a fixed one-parameter
compression or matrices of determinant one.

### Pair gcds after normalizing the columns

Let the two Gaussian columns be `F_1,F_2`. Any Gaussian common divisor
of `U_i=a_i F_1+b_i F_2` and `U_j=a_j F_1+b_j F_2` divides

```text
Delta_ij F_1    and    Delta_ij F_2.
```

Since the columns are Gaussian-coprime, it divides the rational integer
`Delta_ij`. Hence

```text
|gcd(U_i,U_j)| <= |Delta_ij|.                         (5)
```

This is where normalization by the common *Gaussian* column gcd matters.
Merely taking a matrix with coprime integer entries would not establish
(5).

### The product of the individual coordinate contents

Let `r_j=gcd(|Re U_j|,|Im U_j|)`. Applying the adjugate matrix and
using the primitivity of `v_j` shows that `r_j` divides `D`.
For a rational prime `p`, at least one entry of `A` is a `p`-adic unit:
otherwise `p` would divide both Gaussian columns. If `p^k` divides both
`r_i` and `r_j`, the same two linear combinations used above show that
`p^k` divides all entries of `Delta_ij A`. Therefore

```text
min(v_p(r_i),v_p(r_j)) <= v_p(Delta_ij).
```

For nonnegative numbers `e_1,...,e_M`,

```text
sum e_j - max e_j <= sum_(i<j) min(e_i,e_j).
```

Applying this to the content valuations, and using
`max_j v_p(r_j)<=v_p(D)`, gives

```text
product_j r_j <= D S.                                (6)
```

Individual contents can grow; (6) controls their product. Omitting
them would give an incorrect extension from `SL_2(Z)` to arbitrary
rational projective transformations.

### Lower bound for the radius

Use the primitive reductions `d_j` from Section 1. They satisfy

```text
|d_j| >= h_j/(sqrt(2) r_j).
```

Their pair Gaussian gcds still satisfy (5), and their common gcd has
modulus at most `d_*`. The elementary lcm inequality

```text
|lcm(d_j)| >= product_j |d_j| /
                       product_(i<j)|gcd(d_i,d_j)|
```

follows prime by prime from the same exponent inequality just used.
Combining it with (2), (5), and (6) yields

```text
R >= [2^(-M/2)/(S^2 d_*)] product_j h_j / D.          (7)
```

This argument includes fractional `beta`. It does not assume that a
common rotation factor is itself integral.

### Geometry of a fixed projective configuration

Every `v_l` is a linear combination of any two distinct `v_i,v_j`:

```text
v_l = (Delta_lj/Delta_ij) v_i
      + (Delta_il/Delta_ij) v_j.
```

Consequently `H<=K max(h_i,h_j)` for every pair. Taking the two
smallest heights shows that all heights except possibly the smallest
are at least `H/K`, so

```text
product_j h_j >= m (H/K)^(M-1).                      (8)
```

Choose signs of the `U_j` so their arguments occupy a half-angle
interval of length `delta=Delta/2`. This is possible because their
doubled arguments occupy an interval of length less than `pi`.
The determinant of two image vectors is `det(A) Delta_ij`.
Apply this to a smallest-height and a largest-height vector. Since
`|sin t|<=|t|`,

```text
delta >= D d_* /(m H).                               (9)
```

Equations (7)--(9), followed by `m<=H`, give (3). They also give
`Delta>=2D d_*/H^2`, which proves the radius bound in (4).

## 3. The four-point exponent cannot be replaced by a positive one

The disappearance of the power of `H` at `M=4` is real. Take

```text
U_j(T)=T+j+i,             j=0,1,2,3,
z_j(T)=U_j(T) product_(l!=j) conj(U_l(T)).
```

These arise from the fixed primitive vectors `(j,1)` under the
unimodular matrices `[[1,T],[0,1]]`. Their common radius is

```text
R_T = product_j |U_j(T)| ~ T^4.
```

Their angular span is

```text
Delta_T = 2(arctan(1/T)-arctan(1/(T+3))) ~ 6/T^2,
```

and hence `Delta_T sqrt(R_T) -> 6`. Thus a positive power of `H`
cannot replace `H^0` in (3) for four points. This example says nothing
about whether every four-point rational template survives, or whether
four points exist below a specified small arc constant.

In contrast, (3) gives a positive lower bound for normalized span on
*each fixed* four-point rational projective orbit. Thus four-point
families whose normalized spans tend to zero must have varying rational
cross-ratios.

Explicitly, normalize the four directions to
`(1,0),(0,1),(1,1),(a,b)`, with `(a,b)` primitive, and put
`B=max(1,|a|,|b|)`. Then `d_*=1`,
`S=|ab(a-b)|<=2B^3`, and `K<=4B`. Equation (3) gives

```text
Delta sqrt(R) >= 1/(16 B^(9/2)).                     (9a)
```

Thus a four-point arc of normalized length at most `C` has cross-ratio
height at least `(1/(16C))^(2/9)`. In particular, if a hypothetical
unbounded cluster is windowed to make its normalized length tend to
zero, every four-point subconfiguration in those windows has unbounded
rational cross-ratio height. Convergence of its real cross-ratio is
still possible.

## 4. An unconditional cross-ratio height consequence for five points

For any five distinct lattice points on a centered circle, their relative
unit-circle directions have rational half-angle parameters. Normalize
three of these parameters by a rational projective transformation to
the directions `(1,0),(0,1),(1,1)`. Write the two remaining primitive
directions as `(a,b),(c,d)`, and put

```text
B = max(1,|a|,|b|,|c|,|d|).
```

Equivalently, `B` bounds the reduced numerator and denominator heights
of two cross-ratio coordinates. The five directions are distinct, so
all determinants used in the theorem are nonzero. For this canonical
configuration,

```text
d_* = 1,
S = |a b c d (b-a)(d-c)(ad-bc)| <= 8 B^8,
K <= 4 B^2.
```

For `M=5`, the constant in (3) therefore satisfies

```text
c >= 2^(-29/4) B^(-12).
```

If the five points lie on an arc with `Delta<pi` and
`Delta sqrt(R)<=C`, (4) gives

```text
R <= 2^56 C^10 B^96.                                 (10)
```

In particular,

```text
B >= (R/(2^56 C^10))^(1/96).                         (11)
```

The constants and exponent are intentionally crude. The point is that
every five-point endpoint family must have projective arithmetic
complexity growing at least as a positive power of the radius. This
conclusion does not require a fixed template hypothesis: it quantifies
how the template must move.

## 5. What remains missing

Equation (11) is an obstruction to bounded projective complexity, not
a point-count bound. Rational cross-ratios of actual lattice clusters
can have heights growing with `R`; none of the arguments above bounds
those heights independently of `R`, nor forces a large cluster to
contain five points with bounded cross-ratio heights.

Thus it would be incorrect to extract a fixed rational template from
an arbitrary unbounded family by compactness alone. Real cross-ratios
may converge while their reduced rational heights tend to infinity.
Any continuation based on projective symmetry must retain this
denominator growth. The present proof isolates it quantitatively but
does not remove it, and the uniform endpoint bound remains unproved.
