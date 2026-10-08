# A quadratic-composition obstruction for the full three-row family

No real rational quadratic polynomial can simultaneously split the two
Gaussian linear factors with roots `6-3i` and `6-4i` into linear factors
over `Q(i)`. Consequently the full three-row polynomial family described
below cannot be extended to four full Boolean rows by substituting one
real rational quadratic and splitting every old block into two Gaussian
linear blocks.

This is an obstruction to that specific extension mechanism. It neither
excludes other four-row polynomial constructions nor supplies a uniform
short-circle bound. The proof uses an explicit two-isogeny descent and
small exact reduction counts; it does not assume the perfect-cuboid
conjecture or infer nonexistence from a bounded search.

## 1. The two factors already obstruct the composition

The [balanced three-row family](balanced_three_row_gaussian_polynomial_family.md)
has seven Gaussian linear blocks

```text
H_123=t-i,
H_12=t-6-5i,       H_13=t-6+3i,       H_23=t-6+4i,
H_1=t-4+3i,        H_2=t-3+2i,        H_3=t-7-6i.
```

Its two relevant roots are

```text
r_13=6-3i,         r_23=6-4i.
```

Suppose `q(s)=A s^2+B s+C`, with rational `A!=0`, splits both
`q(s)-r_13` and `q(s)-r_23` over `Q(i)`. Completing the square gives

```text
q(s)=A(s+B/(2A))^2+d,       d=C-B^2/(4A) in Q.
```

Thus both `(r_13-d)/A` and `(r_23-d)/A` must be squares in `Q(i)`.
The norm of a Gaussian rational square is a rational square. With
`z=6-d`, this implies

```text
z^2+9 = u^2,       z^2+16 = v^2,       z,u,v in Q.       (1)
```

Sections 2--4 prove that (1) forces `z=0`. In that case the quotient of
the two proposed Gaussian squares would be

```text
((-3i)/A)/((-4i)/A)=3/4.
```

But `3/4` is not a square in `Q(i)`: a square `(x+iy)^2` that is
positive rational has `2xy=0`, hence `y=0`, and must already be a
square in `Q`. This contradiction proves the asserted obstruction,
including nonmonic quadratics and negative leading coefficients.

## 2. A complete descent for the necessary elliptic curve

A solution of (1) gives a rational point

```text
(x,y)=(z^2,zuv)
```

on

```text
E: y^2=x(x+9)(x+16)=x^3+25x^2+144x.
```

Its two-isogenous curve is

```text
E': y^2=x^3-50x^2+49x=x(x-1)(x-49).
```

For clarity, the standard descent input is stated explicitly. On
`y^2=x^3+a x^2+b x`, define `alpha(O)=1`, `alpha((0,0))=b`,
and `alpha((x,y))=x` in `Q*/Q*^2` otherwise. Its image is a subgroup;
every image class has a squarefree representative `b_1|b`. A class
`b_1` occurs exactly when the corresponding homogeneous quartic

```text
N^2=b_1 U^4+a U^2 V^2+(b/b_1)V^4                      (2)
```

has the requisite rational point, equivalently an integral solution
with `gcd(U,V)=1`, allowing the usual point-at-infinity exceptions.
For the isogenous pair, the rank is

```text
2^rank(E(Q)) = |im(alpha_E)| |im(alpha_E')| / 4.         (3)
```

These are the ordinary two-isogeny descent statements, as presented in
[Dujella, Section 15.5, equations (15.29)--(15.31), pp. 507--509](https://web.math.pmf.unizg.hr/~duje/pdf/section15-5.pdf).
No analytic-rank assertion or unproved local-to-global implication is
used below: every excluded class has an explicit local obstruction,
and every retained class on `E` is realized by a displayed point.

## 3. The two descent images

For `E`, the possible squarefree classes are

```text
1,-1,2,-2,3,-3,6,-6.
```

The classes `1,-1,3,-3` occur: use `O`, `(-9,0)`, `(12,84)`,
and `(-12,12)`, respectively. The points lie on `E` by direct
substitution. Classes `2` and `-2` would give, respectively,

```text
N^2= 2U^4+25U^2V^2+72V^4,
N^2=-2U^4+25U^2V^2-72V^4.                              (4)
```

Neither has a primitive integer solution. The possible right sides
modulo 16 are as follows; at least one of `U,V` must be odd.

| Parity of `U,V` | Class `2` | Class `-2` |
|---|---|---|
| odd, odd | `3,11` | `7,15` |
| odd, even | `2,6` | `2,14` |
| even, odd | `8,12` | `8,12` |

None belongs to the square residues `0,1,4,9`. Since the image is a
group containing class `3`, a class `6` or `-6` would give class `2`
or `-2` after multiplication by `3`. Therefore

```text
im(alpha_E)={1,-1,3,-3},       |im(alpha_E)|=4.           (5)
```

For `E'`, the possible classes are `1,-1,7,-7`. The polynomial
`x(x-1)(x-49)` is negative for `x<0`, so rational points have
nonnegative `x` and both negative classes are impossible. Class `7`
would require

```text
N^2=7U^4-50U^2V^2+7V^4.                                (6)
```

Modulo 7 this is `N^2=-(UV)^2`. As `-1` is not a square modulo 7,
one of `U,V` and also `N` must be divisible by 7. Primitivity makes
exactly one of `U,V` divisible by 7. The right side of (6) then has
7-adic valuation exactly one, contradicting its being a square.
Equivalently, there is no primitive solution modulo 49. Hence

```text
im(alpha_E')={1},       |im(alpha_E')|=1.                (7)
```

Equations (3), (5), and (7) prove `rank(E(Q))=0`.

One can also read this conclusion directly from the descent quotients.
For the dual isogenies `phi:E->E'`, `psi:E'->E`, the standard kernel
identities give `E'(Q)=phi(E(Q))` from (7), and then
`psi(E'(Q))=2E(Q)`. Equation (5) says `[E(Q):2E(Q)]=4`.
Since `E` has all four rational two-torsion points, finite generation
forces its free rank to be zero.

## 4. All rational points, and the square-coordinate condition

The discriminant of the displayed integral model is

```text
Delta_E=16*144^2*49,
```

so 5 and 11 are primes of good reduction. Direct counts give

```text
|E(F_5)|=8,       |E(F_11)|=8.
```

For either good prime, the prime-to-that-prime torsion injects into the
reduction group. Using the two primes removes any exceptional primary
part and proves that the rational torsion order divides 8. There are
already eight distinct rational points:

```text
O, (0,0), (-9,0), (-16,0), (12,84), (12,-84),
(-12,12), (-12,-12).                                  (8)
```

Because the rank is zero, these are all of `E(Q)`. Its only nonnegative
finite `x` coordinates are `0` and `12`, and `12` is not a rational
square. The point with `x=z^2` therefore has `x=0`, proving that
(1) forces `z=0`.

This also supplies an explicit proof of the familiar obstruction to a
rational Euler brick with two edges in ratio `3:4`: the prospective
third edge would be a nonzero `z` in (1). Only this fixed-ratio statement
is used; the general perfect-cuboid problem is irrelevant.

The [checker](check_quadratic_extension_gaussian_square_obstruction.py)
exhausts the complete local residue certificates modulo 16 and 49,
checks the two good-reduction counts and all eight rational points,
and verifies their group closure and orders with exact fractions.
These finite calculations certify the finite ingredients of the
descent proof, not a bounded search for rational values of `d`.
