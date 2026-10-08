# Gale interpolation, square quartics, and a seven-row counterexample

For seven rank-two Gale rows, the zero-sum relation and the seven signed
star-square conditions do not force the rows to come from rational circle
data.  The missing condition has a simple form: the polynomial interpolated
from the row scalings must be a square quartic.  An explicit example from
seven multiples of one rational elliptic-curve point satisfies the zero-sum
and primitive star-square conditions but has interpolant `t^3-2`.

This is an inverse-data obstruction only.  It neither constructs a circle
configuration nor proves or disproves any radius-height inequality.

## 1. Zero row sums are equivalent to a degree-four interpolant

Let `t_1,...,t_n` be distinct rational numbers, put

```text
P(t)=product_i(t-t_i),       D_i=P'(t_i),             (1)
```

and let `F` be the unique polynomial of degree at most `n-1` satisfying

```text
F(t_i)=a_i D_i.                                      (2)
```

Partial fractions give the exact identity

```text
F(t)/P(t)=sum_i a_i/(t-t_i).                         (3)
```

At infinity, the right side is

```text
t^(-1) sum_i a_i+t^(-2) sum_i a_i t_i
 +t^(-3) sum_i a_i t_i^2+... .                       (4)
```

For seven rows `a_i(1,t_i)`, equations `sum_i a_i=sum_i a_i t_i=0`
are therefore equivalent to `deg F<=4`.  For eight conic Gale rows
`a_i(1,t_i,t_i^2)`, the three zero-sum equations are again equivalent to
`deg F<=4`.

For the actual seven- or eight-point circle Gale formula,

```text
a_i=c q(t_i)^2/D_i                                   (5)
```

for a nonzero common scalar `c` and a rational quadratic `q`.  Equations
(2) and (5), at more than four distinct points, show that

```text
F(t)=c q(t)^2.                                       (6)
```

Conversely, (6) gives (5).  Thus, after fixing projective coordinates, the
square-quartic condition is exactly the algebraic scaling condition.  To
come from the circle problem, `q` must additionally be positive definite
and its Gram determinant must be a positive rational square. In homogeneous
coordinates, (6) says that a binary quartic is a scalar times the square
of a binary quadratic. This property is invariant under rational
projective changes of the parameter, so a different affine chart does
not evade the condition.

## 2. The square-quartic locus has no quadratic coefficient equations

Write the coefficients of

```text
(A t^2+B t+C)^2
```

in descending order.  The coefficient map is

```text
[A,B,C] -> [A^2, 2AB, B^2+2AC, 2BC, C^2].           (7)
```

There is no nonzero quadratic equation in the five target coordinates
which vanishes on (7).  Indeed their 15 pairwise quadratic products span
all 15 degree-four monomials in `A,B,C`.  In the natural lexicographic
orders, the exact coefficient determinant is `-16384`.  Hence equations
cutting out the square-quartic surface begin in degree at least three.

This explains why merely searching for quadratic relations among the five
coefficients of `F` cannot detect the missing condition.

## 3. Seven exact rational rows with square stars

Consider

```text
E: y^2=t^3-2,       Q=(3,5).                         (8)
```

Using the rational chord-and-tangent law, form the seven points
`Q,2Q,...,7Q`.  Direct exact arithmetic verifies that they are finite and
their seven `t`-coordinates are distinct; no assertion about the order or
infinitude of `Q` is needed.  Write them as `(t_i,y_i)` and set

```text
F_i=y_i^2=t_i^3-2,
a_i=F_i/D_i,
b_i=a_i(1,t_i).                                      (9)
```

The interpolant in (2) is exactly

```text
F(t)=t^3-2.                                          (10)
```

It has degree three, so (3)-(4) give

```text
sum_i b_i=0.                                         (11)
```

Let the rational Plucker coordinates be

```text
p_ij=det(b_i,b_j)=a_i a_j(t_j-t_i),
p_ji=-p_ij.                                          (12)
```

For each label define its signed Gale star

```text
S_i=product_(j!=i) p_ij.                             (13)
```

Since there are seven labels,

```text
product_i D_i=-product_(i<j)(t_i-t_j)^2.             (14)
```

Also, directly from (9), (12), and (13),

```text
S_i=(product_j a_j) F_i^5/D_i^4.                    (15)
```

Every `F_i=y_i^2` is a rational square.  Equations (14)-(15) therefore
show that every `S_i` is the negative of a rational square.

Clear the common denominator of the 21 coordinates `p_ij` and divide by
their integer gcd.  Call the resulting primitive decomposable coordinates
`P_ij`.  This is the Plucker normalization of the saturated rank-two
lattice.  The common rescaling from `p_ij` to `P_ij` changes every
six-factor star by a sixth power.  Hence

```text
-product_(j!=i) P_ij
```

is a rational square and an integer, so it is an integer square.  Thus the
primitive lattice has all seven required negative square stars.

Nevertheless, (10) cannot equal a nonzero scalar times the square of a
quadratic: a nonzero polynomial square has even degree.  By Section 1,
these zero-sum rank-two Gale data cannot be the kernel data of rational
circle points.  The example proves that zero row sums, decomposability,
saturation, and all primitive signed star squares are not sufficient for
the inverse circle reconstruction.

The exact checker
[`check_gale_square_quartic_characterization.py`](check_gale_square_quartic_characterization.py)
computes `Q,...,7Q`, verifies (8)-(15), performs the primitive Plucker
normalization and every four-label Plucker relation, checks all seven
integer star square roots, recovers the interpolant (10), and verifies the
rank-15 assertion following (7).
