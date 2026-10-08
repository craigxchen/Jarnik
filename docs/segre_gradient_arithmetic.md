# Exact Segre polar coordinates and their arithmetic heights

The six-row matching invariants and the balanced-kernel triangle products
are connected by an exact cubic gradient map. In a full centrally truncated
seven-row profile, their primitive projective heights are respectively
`18w+o(w)` and `16w+o(w)`. The common divisors responsible for these values
can be bounded from both sides, including the actual integer corrections.
The inverse polar map has an exact Vandermonde factor and returns the
original matching point. These facts give a simultaneous geometric
reformulation; they do not yet give a new circle configuration, an iterated
height descent, or a uniform bound.

For geometric context, the six-point quotient is the Segre cubic, and its
projective dual is the Igusa quartic. The source also explains the outer
`S_6` action and the interpretation through six points on a conic:
[Howard--Millson--Snowden--Vakil, Section 1](https://math.stanford.edu/~vakil/files/0906.2437v1.pdf).
The formulas and arithmetic estimates below are independently checked by
[check_segre_gradient_arithmetic.py](check_segre_gradient_arithmetic.py),
which uses exact polynomial arithmetic and integer linear algebra.

## 1. A small integral matching basis and its cubic

For six binary vectors `v_i=(x_i,y_i)`, write
`Delta_ij=x_i y_j-x_j y_i`. Use the five matching coordinates

```text
A=Delta_12 Delta_34 Delta_56,
B=Delta_12 Delta_36 Delta_45,
C=Delta_14 Delta_23 Delta_56,
D=Delta_16 Delta_23 Delta_45,
E=Delta_16 Delta_25 Delta_34,
L=A+B+C+D+E.
```

They satisfy the exact cubic identity

```text
Phi(A,B,C,D,E)=BCE-ADL=0.                                 (1)
```

For a direct verification, use the dense projective chart
`(infinity,0,1,a,b,c)`, with vectors `(1,0),(0,1),(1,1),
(a,1),(b,1),(c,1)`. The five coordinates become

```text
((1-a)(b-c), (1-c)(a-b), c-b, b-a, b(a-1)).
```

Substitution proves (1); row homogeneity extends it to all binary
vectors. The fifteen matching products are integral linear combinations
of these five, with coefficient sums at most five. Conversely these five
are among the fifteen. In particular their numerical gcds agree exactly.

The cubic gradient is

```text
G_A=-D(2A+B+C+D+E),
G_B=CE-AD,
G_C=BE-AD,
G_D=-A(A+B+C+2D+E),
G_E=BC-AD.                                                (2)
```

## 2. The gradient is exactly the integral triangle lattice

For a sorted triple `abc`, put `T_abc=Delta_ab Delta_bc Delta_ca`.
Let `T_J` here denote the product of the two triangle invariants for
the partition `J | J complement`, with `1 in J`. In the order
`123,124,125,126,134,135,136,145,146,156`, the exact identities are

```text
T_123 = G_E,
T_124 = -G_D+G_E,
T_125 = G_C-G_D,
T_126 = G_C,
T_134 = G_B-G_D,
T_135 = -G_A+G_B+G_C-G_D+G_E,
T_136 = -G_A+G_C,
T_145 = -G_A+G_E,
T_146 = -G_A+G_B,
T_156 = G_B.                                              (3)
```

There are integral inverse formulas, for example

```text
G_A=2T_123-T_124+T_125+T_134-T_135,
G_B=T_123-T_124+T_134,
G_C=T_123-T_124+T_125,
G_D=T_123-T_124,
G_E=T_123.                                                (4)
```

Thus the gcd of the five gradient values equals the gcd of all ten
triangle-product values exactly. Moreover

```text
max|T_J| <=5 max|G_i|,
max|G_i| <=6 max|T_J|.                                   (5)
```

This is also precisely the balanced kernel discussed in
[small_invariant_relations_at_balanced_cuts.md](small_invariant_relations_at_balanced_cuts.md).
The fifteen quadratic monomials in `A,...,E` are linearly independent.
Their evaluations on all twenty balanced two-direction specializations
have rank ten. The five independent quadrics (2) vanish at all of them.
The exact row reduction has integral entries, and the gradient quadrics
have a signed permutation matrix on the five free coordinates. Hence
the balanced kernel in the integral quadratic matching-monomial lattice
is exactly the integral span of (2).

In particular, a small numerical zero relation already known to lie in
this kernel becomes a small-coefficient hyperplane equation on the
gradient point. This identifies its common geometry; it does not exclude
such a hyperplane equation at an arithmetic point.

## 3. Actual central corrections and exact gcd control

Now take six nonanchor rows from the actual seven-row central extraction.
For each of the `64` subsets `S` of the nonanchor labels, let

```text
n_S=N(H_S),       (1-eta)w<=log n_S<=(1+eta)w.
```

The rational integers `n_S` are pairwise coprime. For every pair of
nonanchor rows the actual central gcd factorization is

```text
g_ij=(product_(S containing i,j) H_S) L_ij,
log|L_ij|<=s.
```

This stronger bound is part of the actual all-anchor central system;
see [all_anchor_primitive_transition.md, Section 3](all_anchor_primitive_transition.md).
Consequently, with `0<|t_ij|<=T`, the real bracket has absolute value

```text
|Delta_ij|=c_ij b_ij,
c_ij=product_(S containing i,j)n_S,
b_ij=N(L_ij)|t_ij| in Z_(>0),
log b_ij <= beta=2s+log T.                                (6)
```

Primes of `b_ij` may overlap any core norm. Nothing below assumes
their coprimality. If only the weaker abstract hypotheses
`P_i=K_i A_i`, `|K_i|<=exp(s)` are available, one may instead take
`beta=4s+log T`. The sharper value in (6) uses the stated actual
central pair factorization.

For either all fifteen matching coordinates or all ten triangle
products, let `G_0` be the gcd computed from their core factors alone,
and `G_actual` their actual numerical gcd. Every graph uses any edge
at most once. Therefore

```text
G_0 | G_actual | G_0 product_(all15 pairs) b_ij.            (7)
```

To prove the upper bound at any rational prime, choose a coordinate
attaining the smallest core valuation. Its additional valuation is at
most the sum of all fifteen edge correction valuations. This argument
also applies to primes absent from the core and permits arbitrary
prime powers. In logarithms, the extra common gcd in (7) is at most
`15 beta=30s+15log T` for each coordinate family.

The core exponents depend only on `a=|S|`. A perfect matching has
at least `g(a)=max(0,a-3)` edges entirely inside `S`, and a matching
attaining that value exists. For the two-triangle graph, the minimum is

```text
a       0  1  2  3  4  5  6
g(a)    0  0  0  0  1  2  3
f(a)    0  0  0  1  2  4  6.                             (8)
```

Indeed its order is
`binom(|J intersect S|,2)+binom(|J complement intersect S|,2)`;
balancing the two counts proves the second line. Thus

```text
G_match,0=product_S n_S^g(|S|),
G_triangle,0=product_S n_S^f(|S|),
sum_(a=0..6) binom(6,a)g(a)=30,
sum_(a=0..6) binom(6,a)f(a)=80.                           (9)
```

Equations (3)--(4) transfer the second gcd statement exactly to the
cubic gradient. There is no unaccounted cancellation among its sums.

## 4. The height change is 18w to 16w

For a nonzero integral vector `v`, use projective height
`h([v])=log max_i|v_i|-log gcd_i(v_i)`.
Each core bracket contains sixteen block norms. Every matching
therefore has logarithmic size between `48(1-eta)w` and
`48(1+eta)w+3 beta`. Every triangle product has logarithmic size
between `96(1-eta)w` and `96(1+eta)w+6 beta`.
They are all nonzero, because the original directions are distinct.

Combining these bounds with (7)--(9), the height using all fifteen
matching coordinates obeys

```text
18w-78eta w-15beta <= h_match,15
                      <=18w+78eta w+3beta.                (10)
```

For all ten triangle coordinates,

```text
16w-176eta w-15beta <= h_triangle,10
                       <=16w+176eta w+6beta.              (11)
```

The chosen five matching coordinates and five gradient coordinates obey

```text
h_match,15-log5 <= h([A:B:C:D:E]) <= h_match,15,
h_triangle,10-log5 <= h([G_A:G_B:G_C:G_D:G_E])
                                      <=h_triangle,10+log6. (12)
```

Hence in the increasingly accurate full-profile extraction,

```text
h([A:B:C:D:E])=18w+o(w),
h([grad Phi])=16w+o(w).                                  (13)
```

The normalized forward gradient has a large exact common factor.
Let `g=gcd(A,B,C,D,E)` and `h=gcd(G_A,...,G_E)`, and let
`X=(A,B,C,D,E)/g`, `Y=(G_A,...,G_E)/h`. Then

```text
grad Phi(X)=(h/g^2)Y,
h/g^2 is a positive integer,
log(h/g^2)=20w+o(w).                                     (14)
```

At the core level its excess exponent `f(a)-2g(a)` is one precisely
on the twenty cuts of size three, and zero elsewhere. Formula (7)
controls the remaining correction factors.

## 5. Exact inverse quartic and the return of the original point

In dual coordinates `(p,u,v,q,w_0)` corresponding in order to
`(G_A,G_B,G_C,G_D,G_E)`, put

```text
Z=pq-uv-u w_0-v w_0,
J=Z^2-4u v w_0(u+v+w_0-p-q).                             (15)
```

The change from `w` to `w_0` in this displayed formula distinguishes
the fifth coordinate from the arithmetic block weight. This is an
equation for the dual Igusa quartic in the chosen matching basis.
For an algebraic irreducibility check, its discriminant as a quadratic
in `p` is `16u v w_0(u-q)(v-q)(w_0-q)`, which is not a square
in `Q(u,v,w_0,q)`. Its coefficients are primitive over
`Q[u,v,w_0,q]`, so the quartic is irreducible.
Write

```text
V=product_(1<=i<j<=6) Delta_ij.
```

There is the exact polynomial identity in the original binary vectors

```text
grad J(grad Phi(A,B,C,D,E))=-4V(A,B,C,D,E).                 (16)
```

The certificate verifies each of the five coordinates by expansion in
six independent affine row variables. Row homogeneity proves the binary
identity as well. Euler's identities then also give `J(grad Phi)=0`:
the scalar product of (16) with `grad Phi` is
`-4V times 3Phi=0`, while the left side is `4J(grad Phi)`.

One direct way to derive (15) is to write `r=AD`. At a generic point
of the cubic, formulas (2) and `BCE=rL` give

```text
Z=2u v w_0/r,
u+v+w_0-p-q=u v w_0/r^2.
```

Eliminating `r` yields (15), and the polynomial identity extends across
the chart exceptions. The separately verified (16) avoids presuming
that any denominator in this derivation remains nonzero.

For the primitive vectors `X,Y` of Section 4, (16) gives

```text
grad J(Y)=mu X,
mu=-4Vg/h^3.                                             (17)
```

Here `mu` is an actual nonzero integer: the left side is integral and
the coordinates of `X` have gcd one. An integral Bezout combination of
those coordinates proves integrality of the scalar. Its absolute value
is exactly the gcd of the inverse gradient coordinates. Since
`log|V|=240w+o(w)`, the exact arithmetic estimates give

```text
log|mu|=30w+o(w).                                        (18)
```

At a core cut of size `a`, its exponent is
`binom(a,2)+g(a)-3f(a)`, equal to one for `a=2,4` and zero
otherwise. The thirty cuts of those sizes account for the inverse
content, just as the twenty balanced cuts account for (14).

Thus the inverse cubic gradient has raw height `48w+o(w)` and common
gcd height `30w+o(w)`, returning the original projective height
`18w+o(w)`. The forward image is on a quartic in the dual projective
space; it is not a new verified circle tuple to which the same forward
construction can be reapplied. Composing the two polar maps returns
the original matching point. No iteration or uniformity conclusion
follows merely from the numerical drop in (13).

## Audit and scope

The exact checker verifies the cubic, both integral span conversions,
the quadratic rank and balanced kernel, all sixty-four cut minima, and
the complete inverse identity (16). The uniformity-audit agent separately
checked the primewise upper bound (7), the actual central correction cost,
and all finite height constants in (10)--(12). The root agent independently
derived the height counts and checked the inverse quartic formulas.
The uniformity-audit agent also read the complete final note and checked
the irreducibility discriminant and the primary-source scope.

The positive conclusion is a fully normalized simultaneous map of the
actual arithmetic data, with controlled forward and inverse common
divisors. A further restriction on its image or on small hyperplanes
through that image remains necessary for a uniform circle-point bound.
