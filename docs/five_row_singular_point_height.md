# Coefficient-height gaps at singular and ramified five-row configurations

Let an irreducible degree-two-in-each-row integer invariant `Q` vanish
at an actual five-row full-cut configuration. If the corresponding
interior moduli point is singular on the anticanonical curve of `Q`, then

```text
8 log C_Q + 80 sigma + 20 log T
    >= (10-70 eta) w - log 32.                         (1)
```

Here `C_Q` is the sum of the absolute coefficients of `Q`. The cut
norm logarithms lie in `[(1-eta)w,(1+eta)w]`, row corrections satisfy
`log|K_i|<=sigma`, and the nonzero primitive pair residues satisfy
`|t_ij|<=T`, in the notation of
[the five-row arithmetic dictionary](five_row_del_pezzo_arithmetic.md).
The constants are absolute and independent of the curve and configuration.

There is also a gap at any rational ramification point of a forgetful
map, including smooth points. If forgetting any one of the five rows
is ramified at the actual point, then

```text
14 log C_Q + 80 sigma + 20 log T
    >= (10-70 eta) w - log 128.                        (1a)
```

The bound (1) holds if the curve is irreducible over `Q` but becomes
reducible over an algebraic closure. Consequently a sufficiently small
irreducible relation defines a geometrically integral curve and its
actual moduli point is smooth. This is a bound for a single relation.
It neither excludes a smooth genus-one curve nor excludes a singular
rational curve whose actual point lies in its smooth locus. The uniform
lattice-arc bound remains unproved.

## 1. An elementary height bound for a singular rational point

Normalize the five ordered directions to

```text
(0,1), (1,0), (1,1), (1,a), (1,b).
```

Distinctness gives rational finite `a,b` with
`a b (a-1)(b-1)(a-b)!=0`. Substitute these rows in `Q`, obtaining
`F(a,b)`. Literal substitution uses only zero, one, and the two
variables, so

```text
deg_a F, deg_b F <= 2,       ||F||_1 <= C_Q.           (2)
```

The compactified moduli surface is
`Bl_((0,0),(1,1),(infinity,infinity))(P^1 x P^1)`.
The section `F` has no boundary component: such a component would
make the corresponding bracket divide `Q`, contradicting
irreducibility and its multidegree. Its restriction to the interior
is irreducible. Indeed, normalization of three distinct rows identifies
the open configuration space with the product of the moduli chart and
the frame and row-scaling variables, and identifies `Q` with a unit
times `F`. A factorization of `F` would therefore factor `Q` in the
localized polynomial ring. Since `Q` is prime and divides no inverted
bracket, one factor would have to be a boundary unit. The absence of
boundary components excludes that possibility. Hence `F` is
irreducible over `Q`.

Suppose the actual point `(a_0,b_0)` is singular. The surface chart
is smooth there, so `F=F_a=F_b=0` at that point. Write

```text
F(a,b)=A(b)a^2+B(b)a+C(b).
```

First suppose `deg_a F=2`. The nonzero polynomial

```text
Delta_b(b)=B(b)^2-4A(b)C(b)                            (3)
```

vanishes at `b_0`. To see the vanishing without a leading-coefficient
assumption, use the identity
`Delta_b=F_a^2-4A F`. To see that the polynomial is nonzero, use
Gauss's lemma: irreducibility of `F` gives irreducibility over `Q(b)`;
a quadratic with identically zero discriminant is reducible in
characteristic zero. With `H=||F||_1`, submultiplicativity gives

```text
||Delta_b||_1 <= ||B||_1^2+4||A||_1||C||_1
             <= (||A||_1+||B||_1+||C||_1)^2 = H^2.   (4)
```

For completeness, if `deg_a F=1`, write `F=A(b)a+B(b)`.
The condition `F_a(a_0,b_0)=0` says `A(b_0)=0`, where `A` is a
nonzero polynomial of coefficient sum at most `H`. This already
gives the height estimate below. In fact `F(a_0,b_0)=0` would also
give `B(b_0)=0`; a rational common root contradicts irreducibility.
If `deg_a F=0`, irreducibility and a rational zero make `F` linear
in `b`, which is smooth. Thus that case cannot be singular.

The rational root theorem for a nonzero integer polynomial `G`
implies

```text
h(r) <= log ||G||_1             whenever G(r)=0, r in Q. (5)
```

Here `h(p/q)=log max(|p|,|q|)` for coprime integers `p,q`, `q>0`.
For a nonzero root, remove initial powers of the variable and apply
the numerator and denominator divisibilities; for the zero root,
the height is zero. Equations (3)--(5), and the corresponding
argument with `a,b` interchanged, prove

```text
h(a_0), h(b_0) <= 2 log H <= 2 log C_Q.                (6)
```

No constant depending on `F` is hidden in (6).

## 2. Comparison with the actual full-cut height

Use the twenty-two bracket graph coordinates of
[the arithmetic dictionary](five_row_del_pezzo_arithmetic.md): twelve
pentagons and ten triangles times squared complementary edges.
Each graph has five brackets and degree two in every row. After
normalization, every bracket is a constant, a variable, or a difference
of two of these. Thus every graph polynomial has bidegree at most
`(2,2)` and coefficient sum at most `2^5=32`.

Write `a_0=p/q`, `b_0=r/s` in lowest terms. Clearing all graph
coordinates by `q^2 s^2` makes them integers of absolute value at
most `32 exp(2h(a_0)+2h(b_0))`. Removing their common divisor only
decreases the maximum. Therefore their primitive projective height
satisfies

```text
h_L <= 2h(a_0)+2h(b_0)+log 32
    <= 8 log C_Q+log 32.                              (7)
```

This is exactly the same projective height as for the original
rows. Simultaneous frame changes and separate row scalings multiply
all degree-two-each graph coordinates by the same nonzero scalar.
There is no frame-height error in this comparison.

For the actual rows, write

```text
|D_ij| = b_ij product_(U containing i,j) Norm(H_U),
1 <= b_ij <= exp(beta),       beta=4 sigma+log T.
```

The correction bound follows because the extra Gaussian pair gcd
has norm at most `Norm(K_i K_j)<=exp(4sigma)`. The graph-gcd
calculation in the arithmetic dictionary gives

```text
h_L >= (10-70 eta)w-20 beta.                          (8)
```

Combining (7) and (8) proves (1). This step explicitly retains
correction factors, primitive residues, the upper and lower fair
weight errors, and the chosen projective-coordinate normalization.

## 3. Geometric irreducibility below the same threshold

A rationally irreducible curve in characteristic zero is geometrically
reduced, and the Galois group acts transitively on its geometric
irreducible components. If the curve is geometrically reducible, a
rational point on one component is fixed by Galois, so lies on every
conjugate component. More explicitly, over the algebraic closure write
`F=c product_j F_j` with distinct irreducible factors. Every `F_j`
vanishes at the rational point. If there are at least two factors,
the product rule gives `F_a=F_b=0` there. Thus the actual rational
interior point is singular and (1) applies.

Thus a relation below (1) is geometrically integral and smooth at
the actual point. Combined with
[the boundary-node height gap](five_row_boundary_node_height.md), a
sufficiently small relation also avoids all fifteen boundary nodes
and meets each of the ten boundary lines at a different rational
nonnodal point. Smoothness at the actual point does not imply that
the whole curve is smooth. A geometrically integral anticanonical
curve has arithmetic genus one; the remaining possibilities include
a smooth genus-one curve and an irreducible singular curve with
rational normalization. Neither is excluded here.

## 4. A gap when even one forgetful map is ramified

Permute the labels so that the forgotten row has coordinate `a` and
the remaining four-point moduli coordinate is `b`. Permutation preserves
`C_Q` and the twenty-two-coordinate height. An irreducible anticanonical
section without boundary components has degree two in `a`: a lower
degree would make its bidegree-`(2,2)` homogenization contain the
boundary fiber `a=infinity`. The forgetful map on the curve is the
degree-two map `(a,b)->b`.

Suppose `F(a_0,b_0)=F_a(a_0,b_0)=0` at the rational interior point.
The discriminant argument again gives `h(b_0)<=2 log C_Q`.
Moreover `A(b_0)!=0`: otherwise `F_a=0` gives `B(b_0)=0` and
`F=0` gives `C(b_0)=0`. Since `b_0` is rational, these three
vanishings make its linear minimal polynomial divide `F`, contrary
to irreducibility. Therefore

```text
a_0=-B(b_0)/(2A(b_0)),
h(a_0)<=2h(b_0)+log(2C_Q)<=5log C_Q+log2.              (9)
```

The first height bound in (9) follows by clearing the square of the
denominator of `b_0` in both numerator and denominator. Their absolute
values are bounded by `2C_Q exp(2h(b_0))`; removing their gcd reduces
the height. Hence

```text
h_L<=2h(a_0)+2h(b_0)+log32<=14log C_Q+log128.          (10)
```

Combining (10) with (8) proves (1a). At a smooth point,
`F_a=0` is exactly the ramification condition for projection to `b`:
the tangent equation is `F_a da+F_b db=0` with `F_b!=0`.
This proof also covers a singular actual point, though (1) is stronger
there.

This has an integral derivative formulation on the original rows.
Since `gcd_G(P_i,bar(P_i))=1`, one has `gcd(X_i,Y_i)=1`.
Euler's identity and `Q(P)=0` give unique integers `lambda_i` with

```text
(Q_(X_i)(P),Q_(Y_i)(P))=lambda_i(-Y_i,X_i).            (11)
```

Indeed a rational multiple is forced by perpendicularity, and a
Bezout combination of `X_i,Y_i` proves its integrality. In a normalized
chart (11) vanishes exactly when `F_a=0`. The normalization factor
contributes only terms proportional to `F`, which vanish at the
actual point. Thus (1a) applies if even one `lambda_i` is zero.
Below this threshold, all five integers are nonzero and all five
forgetful maps are unramified at the actual point.

For a concrete smooth ramified example, the earlier elliptic section

```text
F=-12a^2+(12+13b-b^2)a-12b
```

has `F=F_a=0`, `F_b=-2` at `(a,b)=(2,4)`. Its discriminant
is `(b-1)(b-4)(b^2-21b+36)`. Thus the ramified-point statement
handles some configurations not covered by the singular-point statement.
This example is algebraic and is not an endpoint full-cut family.

## Verification

The singular-point discriminant argument and constants were independently
derived by the Astra worker, who also audited the ramification extension,
its constants, and the derivative comparison. The exact algebraic checks are in
[check_five_row_singular_point_height.py](check_five_row_singular_point_height.py).
An independent [interior nodal certificate](check_five_row_singular_interior_certificate.py)
supplies the geometrically integral anticanonical section

```text
F=120a^2-310ab+120b^2-33a-37ab^2+72b+68a^2b,
```

which has an ordinary node at `(2,3)` and avoids all fifteen boundary
nodes. Its discriminant in `a` is
`(b-3)^2(1369b^2-1486b+121)`; the second factor has nonzero discriminant
`1545600`, so is not a square over an algebraic closure. The coefficients
as polynomials in `b` have no common factor, and each prescribed blow-up
point has multiplicity exactly one. This verifies geometric
irreducibility without assuming it from numerical data. In the directed
six-graph basis of the boundary-slope table its coefficients are
`(116,-33,72,-4,-4,0)`; all thirty values
`alpha_S,beta_S,alpha_S+beta_S` are nonzero. Thus the singular-point
stratum is not already contained in the older boundary-node case.
