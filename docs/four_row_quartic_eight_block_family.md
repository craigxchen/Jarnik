# Four quartic rows with eight independent root blocks

There is an exact real one-parameter family of four quartic
constant-imaginary numerators with support pattern

```text
{1}, {2}, {3}, {4}, {123}, {124}, {134}, {234}.              (1)
```

Each row has degree four, each pair has gcd degree two, and their
lcm has degree eight. Thus it meets the polynomial endpoint
least-radius degree exactly. It disproves a proposed extension of
the three-row example's maximality to all real polynomial families.
The displayed symmetric family has no simultaneous rational
normalized rows; that arithmetic obstruction is proved in Section 5.
It is therefore not a Gaussian-integer construction or an unbounded
count counterexample.

A subsequent [nonsymmetric rational construction](four_row_rational_quartic_construction.md)
realizes the same eight supports over `Q(i)`. Thus the arithmetic
obstruction proved here is specific to the displayed symmetric
parameterization, as its proof and scope require.

## 1. The eight roots

Fix a real number `u>1/2`, and put

```text
a = sqrt(u(u+1)),
b = sqrt(u^2-1/4),
K = u(2u+1).
```

All square roots here are positive. Define the four roots of the
triple blocks by

```text
z_1 =  a + i u,            z_2 = -a + i u,
z_3 =  b - i(u+1/2),       z_4 = -b - i(u+1/2).             (2)
```

The block with root `z_j` belongs to all rows except `j`.
Define the private roots by

```text
w_1 =  a u/(u+1) + i(u+1),
w_2 = -a u/(u+1) + i(u+1),
w_3 =  b(u+1/2)/(u-1/2) - i(u-1/2),
w_4 = -b(u+1/2)/(u-1/2) - i(u-1/2).                        (3)
```

The monic numerator for row `j` is

```text
h_j(t)=(t-w_j) product_(l!=j)(t-z_l).                      (4)
```

The roots in (2)--(3) are all nonreal and distinct. Their four
imaginary levels are `u,u+1,-u-1/2,-u+1/2`; none of the two
positive levels equals the negative of a negative level. At each
level the two nonzero real parts have opposite signs. Consequently
no root in the entire eight-root set is the conjugate of another.
These checks hold throughout `u>1/2`, not merely generically.

## 2. Exact expansions and nonzero constant imaginary parts

Direct multiplication, using `a^2=u(u+1)` and
`b^2=u^2-1/4`, gives

```text
h_1(t) = t^4 + [a/(u+1)]t^3 + (2u+1)t^2
               + (2u+1)a t + K^2 + i K^2/a,

h_3(t) = t^4 - [b/(u-1/2)]t^3 - 2u t^2
               - 2u b t + K^2 + i K^2/b,

h_2(t) = bar(h_1)(-t),
h_4(t) = bar(h_3)(-t).                                    (5)
```

The conjugation in (5) acts on coefficients. Thus

```text
Im h_1= K^2/a,     Im h_2=-K^2/a,
Im h_3= K^2/b,     Im h_4=-K^2/b.                          (6)
```

In particular all these imaginary parts are nonzero constants.
The inequalities `0<b<a` show that the four constants in (6)
are distinct.

One way to discover (2)--(3) is to demand that
`Im(w_j^r-z_j^r)` be independent of `j` for `r=1,2,3`.
With the root pairing in (2), the first two conditions reduce to
`Im w_j=Im z_j+1` and
`Re w_j=Re z_j Im z_j/(Im z_j+1)`. The third condition then
gives precisely the two square-root identities for `a,b`.
The explicit expansion (5) verifies all the resulting coefficient
conditions without relying on a numerical search.

## 3. Primitive pairs and the least-radius degree

Let `tau_j=Im h_j` and normalize

```text
F_j=Re(h_j)/tau_j,       h_j/tau_j=F_j+i.                  (7)
```

The `F_j` are real quartics. Their leading coefficients are the
four distinct numbers `1/tau_j`, so every difference has degree
four. Since the eight roots are distinct and their incidences are
exactly (1),

```text
gcd(h_i,h_j)=product_(l notin {i,j})(t-z_l),
deg gcd(h_i,h_j)=2,
L=lcm(h_1,h_2,h_3,h_4)=product_l(t-z_l)(t-w_l),
deg L=8,             gcd(L,bar L)=1.                       (8)
```

All polynomial gcds and lcms are over `C[t]`, up to constants.
Put `G_ij=gcd(h_i,h_j)`. The primitive pair expression

```text
h_j bar(h_i)/(G_ij bar(G_ij))                             (9)
```

is a polynomial of degree four. Its imaginary part is a nonzero
constant: before division that imaginary part is
`tau_j Re h_i-tau_i Re h_j`, a degree-four real polynomial
whose leading coefficient is `tau_j-tau_i !=0`. Division by
the real monic polynomial `G_ij bar(G_ij)` of degree four
therefore gives a nonzero constant.

Equivalently, (7)--(9) give the exact division condition

```text
F_i-F_j divides F_i F_j+1.                                (10)
```

Thus these are four simultaneous rows of the faithful polynomial
division model, with every primitive pair division retained.

For real `t` tending to positive infinity, `|L(t)|~t^8` and
`arg h_j(t)=tau_j t^(-4)+O(t^(-5))`. The circle construction
using anchor `bar L` and the four ratios `h_j/bar h_j` has five
points, angular span of order `t^(-4)`, and radius of order
`t^8`. Its normalized endpoint arc length stays bounded. This
is a real polynomial construction, not a lattice claim.

## 4. Relation to the previous quartic example

The common gcd of all four rows is constant. Its eight supports
are the four private blocks and the four missing-singleton blocks.
This differs from the seven nonempty supports on three rows in
the earlier example. Restricting this construction to any three
rows gives exactly those seven nonempty supports, each once:
the omitted row's private block disappears and its missing-row
triple block becomes the full cut. Thus any three rows here
have a common linear factor and lcm degree seven.

Accordingly there is no conflict with the exact maximality of the
specified triple in
[quartic_extension_maximality.md](quartic_extension_maximality.md).
That theorem concerns a different triple and was explicitly scoped
to that triple. A universal claim that noncollinear quartic systems
cannot have a fourth row is false over the real numbers.

## 5. Rationality obstruction for this symmetric family

There is no rational `u>1/2` for which both `a` and `b` are
rational. To see this, suppose such a triple exists and put

```text
d=a/u,       e=2b/u,
X=-d^2,     Y=d e.
```

The defining equations imply

```text
e^2=(3-d^2)(1+d^2),
Y^2=X(X-1)(X+3),       -3<X<-1.                           (11)
```

The two-isogenous curve

```text
E: y^2=x^3-x^2+x                                          (12)
```

receives (11) under the dual isogeny with `x`-coordinate

```text
x=(X^2+2X-3)/(4X)=b^2/a^2,       0<x<1.                   (13)
```

Here is a short descent proving that `E(Q)` consists exactly of
`O,(0,0),(1,1),(1,-1)`, which contradicts (13). Write
`E':Y^2=X^3+2X^2-3X` for the curve in (11). For the standard
dual two-isogenies `phi:E->E'` and `psi:E'->E`, the squareclass
maps send a nonzero finite point to its `x`-coordinate, infinity
to `1`, and `(0,0)` to the linear coefficient of the cubic.
Their kernels are `psi(E'(Q))` and `phi(E(Q))`, respectively.
These standard identifications and the explicit isogenies are
given in the primary paper
[Rank one elliptic curves and rank stability, Section 2.1](https://pi.math.cornell.edu/~zywina/papers/Rank1WithStability.pdf).

For an integral model `y^2=x(x^2+A x+B)`, the squareclass of a
nonzero rational `x` is represented by a signed divisor of `B`.
Indeed negative prime valuations of `x` must be even by comparing
valuations in the equation, and positive valuations at primes
not dividing `B` must also be even. On `E`, every nonzero real
point has `x>0`, since `x^2-x+1>0`. Its squareclass image is
therefore `{1}`. On `E'`, the image is contained in
`{1,-1,3,-3}`; the points `O,(-1,2),(3,6),(0,0)` show that
all four occur.

For clarity the factor of two in the rank calculation can be
checked directly. Put `A=E(Q)` and `B=E'(Q)`. Then
`[A:psi B]=1` and `[B:phi A]=4`. The nonzero kernel point of
`psi` has squareclass `-3`, so is not in `phi A`. As
`psi phi=[2]`, this gives

```text
[A:2A] = [A:psi B] [B:phi A+ker psi] = 1*(4/2)=2.
```

The only nonzero rational two-torsion point of `E` is `(0,0)`.
The Mordell--Weil theorem therefore makes the left side
`2^(rank E+1)`, proving rank zero.

Finally direct counts give `#E(F_5)=#E(F_7)=8`; both primes
are of good reduction because the discriminant is `-48`.
Torsion injectivity at good primes bounds its order by eight.
The point `(1,1)` has order four, and the unique rational
two-torsion point makes any group of order eight cyclic. Such
a group would contain a half of `(1,1)` or `(1,-1)`. The
duplication formula would force its rational `x`-coordinate to
satisfy

```text
(x^2-1)^2=4x(x^2-x+1),
x^4-4x^3+2x^2-4x+1=0.                                   (14)
```

The rational-root theorem allows only `x=1,-1`, and neither
is a root. Thus the torsion order is four, proving the asserted
complete list of `E(Q)` and the impossibility of (13).

As an external crosscheck, `E'` becomes
`y^2=x^3-x^2-4x+4` under `x=X+1`; the
[LMFDB record for 24.a4](https://www.lmfdb.org/EllipticCurve/Q/24/a/4)
lists its rank as zero and group as `Z/2 x Z/4`, consistently
with the descent. The proof above does not require that database
classification.

This also excludes simultaneous rational coefficients for the
four normalized polynomials (7) in the displayed parameterization,
even if `u` was initially allowed to be irrational. Indeed (5)
gives

```text
F_1(0)=a, F_2(0)=-a, F_3(0)=b, F_4(0)=-b,
u=a^2-b^2-1/4.
```

If all four `F_j` had rational coefficients, these identities
would force `a,b,u` all rational, already ruled out.

No assertion is made that this symmetric ansatz exhausts all real
or rational realizations of the eight-block support design, or
that every possible real change of parameter has been excluded.
The construction settles the real existence test; obtaining a
Gaussian-integer endpoint family requires further arithmetic.

## 6. A necessary common circle for this support design

This last statement applies to any realization of the exact eight
simple, globally conjugate-disjoint root blocks in (1), with real
quartics `F_i+i` and all pair differences of degree four. It does
not assume the symmetric parameterization above.

Let `z_j` be the four missing-row triple roots and put
`Q_j(t)=(t-z_j)(t-bar z_j)`. The exact primitive pair identities
give nonzero real constants `c_ij` with

```text
F_i-F_j=c_ij Q_k Q_l,       {i,j,k,l}={1,2,3,4}.            (15)
```

Telescoping the differences for rows `1,2,3` and cancelling `Q_4`
shows

```text
c_12 Q_3+c_23 Q_1=c_13 Q_2.
```

The triple `1,2,4` similarly puts `Q_4` in the same real linear
span of `Q_1,Q_2`. The four distinct monic quadratics consequently
lie on one affine line in coefficient space. Equivalently the
points `(Re z_j,|z_j|^2)` are collinear. This gives either a
circle with real center containing all four roots, or a vertical
line on which their real parts coincide. This span argument was
derived independently by the root orchestrator.

The vertical alternative is impossible here. Translate the real
variable to make the common real part zero, and write `z_j=i b_j`.
The `b_j` are four distinct nonzero real numbers. All `Q_j` are
even, so (15) says the four `F_i` have the same odd part, say
`C t^3+E t`. At every `z_j`, an incident row satisfies
`F_i(i b_j)=-i`. Taking imaginary parts gives

```text
-C b_j^3+E b_j=-1       for all four j.                    (16)
```

A cubic with constant term one cannot have these four distinct
roots. Hence only the circle alternative remains. In the family
(2), the common circle is explicitly `|z|^2=K`.

This is a support-specific geometric restriction with exact
linear factors. It is not a classification of all realizations,
and does not by itself supply rational root coordinates.

## 7. The remaining root constraint is an explicit second conic

Translate the real center from Section 6 to zero, and write the
common circle as `|z_j|^2=K`, with `K>0`. Each product `Q_k Q_l`
in (15) has constant coefficient `K^2` times its leading
coefficient, and linear coefficient `K` times its cubic
coefficient. Thus the four normalized polynomials have common
real defects

```text
U = constant(F_i)-K^2 leading(F_i),
V = linear(F_i)-K cubic(F_i).
```

For suitable real `A_i,B_i,C_i`, write

```text
F_i(t)=A_i(t^4+K^2)+B_i(t^3+Kt)+C_i t^2+U+Vt.             (17)
```

At an incident root `z_j=x_j+i y_j`, divide `F_i(z_j)=-i`
by `z_j^2`. The first three terms of (17) become real, since
`K/z_j=bar z_j`. The imaginary part of the remaining expression
therefore gives

```text
x_j^2+y_j^2=K,
x_j^2-y_j^2-2U x_j y_j-VK y_j=0.                         (18)
```

This is a necessary two-conic description for every realization
of this exact support design, beyond the symmetric subfamily.

There is also a direct reconstruction from (18). Suppose four
nonreal, distinct, conjugate-disjoint points satisfy these two
conics. Their real parts are distinct. For each omitted index
`i`, let the real quadratic `q_i(x)` interpolate at the other
three real parts the values

```text
q_i(x_j) = -Re((U+i)/z_j^2+V/z_j),       j != i.
```

Write `q_i(x)=4A_i x^2+2B_i x+(C_i-2K A_i)` and use (17).
Its real and imaginary parts show exactly that `F_i(z_j)=-i`
for `j!=i`. If the resulting polynomials are distinct quartics
and all eight triple and private roots together are distinct,
nonreal, and globally conjugate-disjoint, they realize the
desired eight-block system. These
nondegeneracy checks must still be made; they do not follow
just from solving the two conics. For rational root coordinates,
`K,U,V` rational, the interpolation is entirely rational.

On the unit circle one may take rational half-angle parameters
`r_j` and put

```text
z_j=(1-r_j^2)/(1+r_j^2) + i*2r_j/(1+r_j^2).
```

Substitution into the second conic gives the quartic

```text
r^4+(4U-2V)r^3-6r^2+(-4U-2V)r+1=0.                     (19)
```

Consequently four finite rational parameters with elementary
symmetric sums `e_4=1` and `e_2=-6` give candidates with
`U=(e_3-e_1)/8` and `V=(e_1+e_3)/4`. The reconstruction above
then provides exact polynomial candidates and explicit
nondegeneracy tests. This reduction was derived by the root
orchestrator. It does not restrict the general design to the
rank-zero symmetric slice of Section 5.
The common-circle and two-conic derivations, including the
interpolation converse with these explicit nondegeneracy
conditions, were independently audited by the algebraic route.
