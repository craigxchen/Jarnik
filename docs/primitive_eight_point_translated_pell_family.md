# Eight primitive points attaining the fixed shifted-Pell endpoint bound

There is an infinite family of eight distinct, individually primitive
Gaussian points on one circle, with global Gaussian gcd one, lying in
arcs `C sqrt(R)` for one fixed constant `C`. Four pairwise norm-coprime
large factors have the same logarithmic scale. The eight sign words
form the full odd parity class of the four-dimensional cube.

This attains the existing
[eight-point upper bound for fixed shifted-Pell templates](pell_product_endpoint_obstruction.md).
It establishes sharpness within that specified class, not a new general
construction route or an obstruction to the uniform circle-point count.
The same eight-point family contains a translated four-group
Hadamard quartet with an individually and globally primitive extension.

## 1. The factors and eight points

Use the Pell sequence from
[the canonical four-point family](canonical_four_point_counterexample.md):

```text
lambda=9+4 sqrt(5),
x_r+y_r sqrt(5)=lambda^r,
(x_(r+1),y_(r+1))=(9x_r+20y_r,4x_r+9y_r),
F=1+2i,
H_r=x_r+F y_r=(x_r+y_r)+2iy_r.
```

Choose `n=1 mod60`, `n>=61`, and put

```text
D=H_(n-2), A=H_n, B=H_(n+2), C=H_(n+4).
```

For each sign word on `D,A,B,C` with an odd number `p` of unbarred
factors, use the correction

```text
K_p=F       when p=1,
K_p=bar(F)  when p=3.
Z_word=K_p times the sign word.                         (1)
```

Every point has norm `5 n_D n_A n_B n_C`, where `n_G=Norm(G)`,
and the common circle radius is

```text
R=sqrt(5) |DABC|.                                      (2)
```

There are four words with `p=1` and four with `p=3`, giving eight
points. All four shifts and all correcting factors are independent
of `n`, as required by the fixed-template theorem.

## 2. Primitivity and complete norm disjointness

The Pell identity gives `x_r^2-5y_r^2=1` and `gcd(x_r,y_r)=1`.
For positive `r`, `x_r` is odd and `y_r` is divisible by four. Hence
the coordinates of `H_r` are coprime and have opposite parity: any
common divisor divides `2`, and the real coordinate is odd. Thus
`H_r` is odd and conjugate-primitive, and its norm has only split
rational prime factors.

The original note proves the following useful congruences:

```text
H_(r+d)       =4(2-i)y_d y_r mod H_r,
bar(H_(r+d)) =-4i x_d y_r mod H_r.                     (4)
```

Here `y_r` is a unit modulo `H_r`, since `gcd(x_r,y_r)=1`.
Consequently a common rational norm prime at index gap `d` must
lie over `2`, `5`, or a split prime dividing `x_d y_d`.

For the three possible gaps the exact factorizations are

```text
x_2=161=7*23,                 y_2=72=2^3*3^2,
x_4=51841=47*1103,            y_4=23184=2^4*3^2*7*23,
x_6=16692641=7*23*103681,
y_6=7465176=2^3*3^3*17*19*107.                         (5)
```

The displayed primes other than `2,17,103681` are `3 mod4`
(and `103681` need not be assumed prime in the argument below).
Thus gaps two and four have no possible shared split prime except
five. At gap six the only remaining possibilities are `17` and
prime divisors of `103681`.

Five is absent from all four norms. Indeed, modulo five,

```text
x_r=(-1)^r,       y_r=(-1)^r r.
```

The factor `F` never divides `H_r`, whereas
`bar(F)|H_r` holds exactly when `r=2 mod5`. The four indices here
have residues `4,1,3,0 mod5`.

Only `D,C` have gap six. Modulo `17`, equation (5) gives `y_6=0`,
and the Pell identity gives `x_6=+/-1`. Therefore `Norm(H_r)` is
periodic with period six modulo `17`. Modulo the integer `103681`,
one has `x_6=0` and `5y_6^2=-1`; squaring the Pell pair gives
`(x_12,y_12)=(-1,0)`. Thus `H_(r+12)=-H_r` modulo `103681`, and
its norm has period twelve.

The Pell sequence extends to negative indices by `x_-r=x_r` and
`y_-r=-y_r`. Since `n-2=-1 mod12`, both periodicity statements give

```text
Norm(D)=Norm(H_-1)=Norm(5-8i)=89 mod17 and mod103681.
```

The integer `89` is coprime to `17*103681`. This excludes all
gap-six exceptions without needing primality of `103681`.
Consequently the five norms

```text
Norm(F), n_D, n_A, n_B, n_C
```

are pairwise coprime. Every point in (1) is a product of
conjugate-primitive factors from these disjoint supports, and is
therefore individually conjugate-primitive.

The eight points are distinct: different sign words have opposite
orientations at a nonunit core factor, which cannot be absorbed by
the correcting factor supported over five. Their global Gaussian gcd
is a unit. Each of `D,A,B,C` occurs in both orientations across the
eight points, and the corrections include both `F` and `bar(F)`.

## 3. Balanced factor heights and endpoint geometry

For `r>=1`, the Pell formulas imply

```text
lambda^r/2 <= |H_r| <= lambda^r.
```

All four factor norms therefore satisfy

```text
log n_D, log n_A, log n_B, log n_C=2n log(lambda)+O(1),
(sqrt(5)/16)lambda^(4n+4) <= R <= sqrt(5)lambda^(4n+4).
                                                               (6)
```

The correcting factors in (1) have constant modulus `sqrt(5)`.

Let `theta=arg(F)/2`, so `0<theta<pi/4`. There is an exact phase
formula

```text
H_r=cos(theta) lambda^r exp(i theta)
       (1-i tan(theta) lambda^(-2r)),
arg(H_r)=theta+epsilon_r,
epsilon_r=-atan(tan(theta) lambda^(-2r)).              (7)
```

In particular `|epsilon_r|<=lambda^(-2r)`. A sign word with `p`
unbarred factors has leading phase `(2p-4)theta`. For `p=1` its
correction contributes `2theta`; for `p=3` the correction contributes
`-2theta`. Thus all eight corrected words have leading phase zero.
Their deviations from zero are signed sums of the four errors in
(7), with absolute value at most

```text
E=4lambda^(-2n+4).
```

Thus all eight points lie in one arc of angular width at most `2E`,
and its length obeys the explicit uniform bound

```text
arc length/sqrt(R) <= 8*5^(1/4)*lambda^6.              (8)
```

The same points lie in a translated
rectangle with radial width at most `R E^2/2` and tangential width
at most `2R E`. By (6), these are respectively `O(1)` and
`O(sqrt(R))`, with constants independent of `n`.

This is an endpoint statement for one fixed, possibly large constant
in (8). It does not assert existence for every prescribed positive
constant or a normalized arc length tending to zero.

## 4. The translated Hadamard quartet and its primitive extension

The four points of (1) that leave `D` unbarred are

```text
F D bar(A) bar(B) bar(C),
bar(F) D bar(A) BC,
bar(F) DA bar(B) C,
bar(F) DAB bar(C).                                    (9)
```

Set `A'=bar(A)`, `B'=bar(B)`, and `C'=bar(C)`. The same quartet is

```text
F DA'B'C',
bar(F) DA' bar(B') bar(C'),
bar(F) D bar(A') B' bar(C'),
bar(F) D bar(A') bar(B') C'.                           (10)
```

This is the four-group Hadamard pattern, with common factor `D` and
bounded corrections. The other four points of (1) conjugate `D`
without leaving the endpoint arc. Thus the quartet admits an
individually and globally primitive eight-point extension. No
additional construction or change of radius is required.

## 5. Scope

The [exact checker](check_primitive_eight_point_translated_pell_family.py)
verifies ten indices `n=61+60k`, including all five pairwise coprime
norm supports, eight distinct equal-norm primitive points, their global
Gaussian gcd, the quartet's exact common norm, and the imaginary-coordinate
bound by ordered arithmetic in `Q(sqrt(5))`. It also checks the Pell
resultant congruences and phase identities without floating point.
The proof above covers every stated index, rather than just these fixtures.

The origin-based
[four-witness gap](four_signed_witness_matching_norm_gap.md) is not
contradicted. Small radial deviations from a translated circle point
are different from subpower absolute imaginary coordinates. Here the
tangential width is of order `sqrt(R)`, rather than subpower at the
large-block scale.

The selected quartet (10) has a large common divisor, but the full
eight-point cluster is globally and individually primitive. A proposed
extension based only on adjoining a point reversing that common
factor therefore also fails: all four such additional points in (1)
remain on the same endpoint arc. The eight-point family (1) attains
the fixed shifted-Pell template bound, whose equality classification
requires exactly four width-one active factors and a complete parity
class. This construction has only four large orientation groups
and a fixed number of points. It does not realize the full Boolean
profile required for the general endpoint reduction or contradict
a radius-uniform bound depending on the arc constant.
