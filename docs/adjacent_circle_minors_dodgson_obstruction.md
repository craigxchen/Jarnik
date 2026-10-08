# Adjacent circle minors satisfy Dodgson condensation, but add no size bound

The integer real-basis evaluation matrices on a circle have an exact adjacent
Desnanot--Jacobi recurrence.  Its central minor is again the canonical
lower-degree circle determinant.  Its only unconditional integer separation is
the product of two already-known odd Vandermonde determinants.  The recurrence
therefore does not improve the classical determinant budget by itself.

## 1. Circle-adapted integer matrices

For `k>=1`, order the integer polynomial columns as

```text
B_k=(1; x,y; x^2,xy; ...; x^k,x^(k-1)y).             (1)
```

Given `2k+1` ordered points `P_0,...,P_(2k)` on an origin circle, let `A_k`
be their evaluation matrix on `B_k`, and put `D_k=det A_k`.  These are integer
matrices for lattice points.  Distinct circle points make `D_k` nonzero,
because the restrictions of (1) span the real trigonometric polynomials of
degree at most `k`.

The real-to-Laurent column change gives the exact magnitude

```text
|D_k|=product_(i<j)|P_i-P_j|/(2^(k^2)R^(k^2)).       (2)
```

Thus `|D_k|>=1` is exactly the classical odd-cardinality chord-product lower
bound, not a relaxation of it.

## 2. The adjacent Dodgson identity

Write the last two columns of `A_k` as

```text
c=x^k,       d=x^(k-1)y.
```

For `epsilon in {L,R}` and `eta in {c,d}`, let `E_(epsilon,eta)` be the
`2k` by `2k` determinant obtained by deleting the left or right endpoint row
and the indicated column.  Let `D_(k-1)^circ` be the determinant obtained by
deleting both endpoint rows and both columns `c,d`.  Its rows are the
`2k-1` interior points and its columns are exactly `B_(k-1)`.

Desnanot--Jacobi applied to the two endpoint rows and the last two columns is

```text
D_k D_(k-1)^circ
 =E_(L,c)E_(R,d)-E_(L,d)E_(R,c).                     (3)
```

Every term is an integer.  The identity is valid without genericity
assumptions on the even minors; one or more `E` terms may vanish.

No sign regularity of the four `E` minors is asserted.  The selected real
column systems depend on which top column is omitted, and an even minor may
vanish on an arbitrarily short arc when its leading Wronskian coefficient
vanishes at the arc location.  Thus (3) must be retained as a signed
determinantal identity, not replaced without proof by a positive relation.

## 3. Where the apparent extra angular order goes

If the angular gaps have common small scale `Delta`, equation (2) gives

```text
D_k=O(R^(k(k+1)) Delta^(k(2k+1))),
D_(k-1)^circ=O(R^(k(k-1)) Delta^((k-1)(2k-1))).       (4)
```

The radius powers are the sums of the polynomial degrees in the respective
bases.  Thus the left side of (3) has angular order

```text
k(2k+1)+(k-1)(2k-1)=4k^2-2k+1.                       (5)
```

Each product of bordering even minors has generic order one lower,
`4k^2-2k`.  Their leading jets cancel, and (3) identifies the first surviving
term exactly.  This is the discrete-Toda cancellation suggested by the circle
recurrence, but the surviving quantity is `D_k D_(k-1)^circ` itself.

For lattice points, integrality therefore yields only

```text
|E_(L,c)E_(R,d)-E_(L,d)E_(R,c)|
 =|D_k D_(k-1)^circ| >=1.                            (6)
```

Using (2), (6) is precisely the product of the classical inequalities for the
outer `2k+1` points and the interior `2k-1` points.  It supplies no new power
of `R` or `Delta`.

## 4. Why a hypothetical sign control would not repair the loss

Even in a configuration where all four even minors happen to have controlled
signs, (3) compares two large
positive products.  Total positivity orders their ratio but gives no
radius-independent gap between them.  The integer gap in (6) is already fully
accounted for by the two odd determinants on the left.

Normalizing the `E` minors by their leading tangent-jet powers makes them
bounded quantities and turns (3) into a nontrivial limiting Toda relation, but
those normalizing factors contain powers of `R`, chord scale, and generally
the local slope denominator.  The normalized quantities are not integers.
Applying `|integer|>=1` after this normalization would therefore be a false
extra quantization assumption.

Conversely, keeping the unnormalized even minors integral leaves them at the
large scale whose near-cancellation is measured by (6).  An improvement would
require a divisibility statement forcing the right-hand difference to be a
multiple of an additional circle-dependent integer, or a uniform arithmetic
gap for the ratio of the two even-minor products.  Neither follows from
Desnanot--Jacobi or sign regularity.

## 5. Scope and vanishing minors

Equation (3) is exact for every lattice configuration and remains true when an
even minor vanishes.  Arguments that divide by an even minor, as in quotient
forms of Dodgson condensation, must first exclude this case.  The
determinant form (3) needs no such exclusion.

The result does not say that higher Pluecker relations are useless.  It
identifies the missing ingredient: a new gcd or divisibility law for the even
evaluation minors that is not generated by the odd determinants already in
(3).  Without such arithmetic, adjacent condensation reproduces the usual
Vandermonde bounds and does not force a seven- or nine-point endpoint
contradiction.

## 6. The Gaussian combination of an even-minor pair

There is an exact closed form for the two even minors associated with one set
of `2k` rows.  Let

```text
X=det[B_(k-1),x^k],       Y=det[B_(k-1),x^(k-1)y],
E=X+iY.                                                   (7)
```

Column order is the displayed order; changing conventions only multiplies
`E` by a Gaussian unit.  Modulo the circle equation `x^2+y^2=R^2`, lower
degree terms can be eliminated using `B_(k-1)`, and the top harmonic satisfies

```text
z^k=2^(k-1)(x^k+i x^(k-1)y) modulo lower columns.        (8)
```

Indeed, reducing `y^(2r)` to `(-1)^r x^(2r)` makes every even term in
`Re(x+iy)^k` positive, and their binomial coefficients sum to `2^(k-1)`;
the odd terms give the same coefficient for `i x^(k-1)y`.

Write the row coordinates as nonzero complex numbers `z_1,...,z_(2k)`, put
`P=product_i z_i`, and let

```text
V=product_(i<j)(z_j-z_i).
```

Changing the lower real columns to
`bar(z)^(k-1),...,bar(z),1,z,...,z^(k-1)` has determinant magnitude
`2^(-(k-1)^2)`.  Appending `z^k`, using (8), and then replacing
`bar(z)^j` by `N^j z^(-j)`, where `N=R^2`, gives

```text
E=i^(k-1) 2^(-k(k-1)) N^(k(k-1)/2) V/P^(k-1).         (9)
```

For the displayed ordering the unit is exactly `i^(k-1)`. To check its
sign, put `r=k-1`: permuting the Laurent columns into
`(1,bar(z),z,...,bar(z)^r,z^r)` has sign `(-1)^r`. Each harmonic pair
contributes `-i 2^(1-2j)` on returning to the real basis. Their product
therefore has unit `(-1)^r(-i)^r=i^r`; the appended top column only
contributes the positive factor `2^(-r)`.
The exponent `k(k-1)/2` is integral. Formula (9) is valid when
`X=0` or `Y=0`; since the rows are distinct, `V` and hence `E` are nonzero.

Taking the quotient with the conjugate removes every size factor.  Since

```text
V/bar(V)=(-1)^(k(2k-1)) P^(2k-1)/N^(k(2k-1)),
```

equation (9) yields

```text
E/bar(E)=-P/N^k.                                      (10)
```

Thus the phase of the even
minor pair is exactly the monomial product phase of its `2k` rows.

For the first cases, the powers in (9) are

```text
k=1: E=V,
k=2: E=i*(N/4)V/P,
k=3: E=-(N^3/64)V/P^2.                              (11)
```

These formulas also show what ordinary primitive normalization does.  Dividing
`X,Y` by their positive integer gcd leaves (10) unchanged.  This is the
canonical real two-vector normalization.  Dividing `E` by a nonreal Gaussian
factor is a different operation: it changes the phase by that factor divided
by its conjugate.  It cannot be called primitive normalization without a
separate choice and an independent theorem controlling this oriented Gaussian
content.

All visible factors in (9) come from the ordinary Vandermonde differences,
the source product `P`, the circle norm, and powers of two from the real-basis
conversion.  Consequently the primitive normalization recovers the same
product-phase map already present in monomial descent.  No radius decrease or
new bounded-degree map follows from the even-minor pair alone.  Axis cases
`X=0` or `Y=0` simply make (10) a unit phase and do not invalidate the formula.

The exact checker `check_adjacent_circle_even_minors.py` verifies (3), (9), its
norm, and the phase identity (10) for `k=1,2,3,4` on an integer circle.
It uses integer Bareiss determinants and Gaussian integer arithmetic only.

## 7. The four-row case is the existing conic-pencil factor

For cyclically ordered four rows, use the exact notation of
`moving_quadrilateral_fifth_point_product.md`:

```text
(z_1-z_0)(z_3-z_2)=b gamma,
(z_2-z_0)(z_3-z_1)=a gamma,
(z_3-z_0)(z_2-z_1)=(a-b) gamma,
Gamma=Norm(gamma),    4sN=Gamma a b(a-b),
```

where `a>b>0` are coprime integers and `s` is a positive integer.
This fixes the unit of `gamma` geometrically. The six-difference
Vandermonde is `a b(a-b) gamma^3`, and
`gamma/bar(gamma)=P/N^2`. Hence `gamma^2/P=Gamma/N^2`, and (9) gives

```text
det[1,x,y,x^2]+i det[1,x,y,xy] = i s gamma.            (12)
```

Thus the four-row even-minor vector contains exactly the already retained
matching factor and pencil factor. Its phase is meaningful with the stated
cyclic convention, but it supplies no new arithmetic quantity beyond them.
