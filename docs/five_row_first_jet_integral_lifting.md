# Exact integral lifting of five-row gradient stresses

The common-invariant condition adds no new growing core-prime index to
the five-row gradient stresses. Let `B=product_(i<j) b_ij`, where the
actual pair corrections satisfy `1<=b_ij<=exp(beta)` and
`beta=4sigma+log T`. The image of the integral numerical-zero invariants
in the lattice of core-divisible Veronese stresses has exact index

```text
I=epsilon delta rho,
epsilon in {1,3},       delta | B^2,       rho | B.
In particular,         I | 3B^3,
                       I<=3 exp(30beta).              (1)
```

Thus multiplying any such stress by `I` makes it the gradient of
**one** integer numerical-zero invariant. This includes stresses with
all five coordinates nonzero. There is no coefficient-height bound on
that invariant: its coefficients may be large. Consequently (1) is an
integral-lifting statement, not a coefficient-height gap or an endpoint
descent. It preserves the common-`Q` condition that the earlier
[cofactor obstruction](five_row_gradient_cofactor_stress_obstruction.md)
left open.

This is the first-jet map of the five-row degree-two-each invariant
space at an actual configuration. It is distinct from the ambient
six-row cubic polar map in
[the Segre gradient calculation](segre_gradient_arithmetic.md).
The proof below uses finite integral polynomial identities and the
actual cut and residue factors; no general-position assertion at core
primes is needed.

## 1. The evaluation kernel and its gradient image

Use the conjugate-primitive Gaussian rows, disjoint odd split cores,
and actual residue factorization

```text
P_i=X_i+iY_i=K_i product_(U containing i) H_U,
n_U=Norm(H_U),
|Delta_ij|=b_ij product_(U containing i,j) n_U,
Delta_ij=X_iY_j-Y_iX_j !=0.                            (2)
```

Choose the integral six-graph basis `Q_0,...,Q_5` certified in
[the gradient calculation](five_row_gradient_discriminant_core_count.md).
Its unimodular polynomial coefficient minor implies that these six
graphs are a basis for the full integer invariant lattice. In particular,
the gcd of their values equals the gcd over all twenty-two graphs.
Write

```text
v_alpha=Q_alpha(P),       g=gcd_alpha |v_alpha|,
Lambda={c in Z^6: sum_alpha c_alpha v_alpha=0}.         (3)
```

Choose integer vectors `u_i` with `det(P_i,u_i)=1`, and let
`J_i(Q)=d_i Q(P)[u_i]`. For `Q(P)=0`, the exact gradient identity is

```text
grad_i Q=lambda_i(-Y_i,X_i),       J_i(Q)=lambda_i.
```

Replacing `u_i` by `u_i+kP_i` changes `J_i(Q)` by `2kQ(P)` and
therefore changes neither its restriction to `Lambda` nor any minor
having the evaluation row first. Define

```text
G_i=product_U n_U^max(0,2|U\{i}|-4),
A_i=G_i(X_i^2,X_iY_i,Y_i^2)^t,
K={ell in Z^5: A ell=0}.                              (4)
```

The discriminant proof gives `G_i|lambda_i` without a correction
denominator. The three infinitesimal `SL_2` identities give `A ell=0`
for `ell_i=lambda_i/G_i`. Hence there is an integer map

```text
Phi: Lambda -> K,
c -> (J_i(sum_alpha c_alpha Q_alpha)/G_i)_i.           (5)
```

Every three old directions are distinct, so `A` has rank three and
`K` is a saturated rank-two lattice in `Z^5`. The image in (5) has
rank two, as also follows from the explicit first-jet identities below.

## 2. The twenty first-jet minors are quadratic sections

For `i<j` and an ordered triple of basis indices `alpha<beta<gamma`, put

```text
H_(ij;alpha,beta,gamma)
 =det [ Q_alpha(P)    Q_beta(P)    Q_gamma(P)
        J_i(Q_alpha)  J_i(Q_beta)  J_i(Q_gamma)
        J_j(Q_alpha)  J_j(Q_beta)  J_j(Q_gamma) ].      (6)
```

Although integer frames were used to define (6), the result is an
integer polynomial in the five rows. Indeed one can instead use
`u_i=(-1/Y_i,0)` and `u_j=(-1/Y_j,0)` over the rational function
field. The numerator is the determinant with `partial_(X_i)` and
`partial_(X_j)` rows, divided by `Y_iY_j`. Euler's identity shows
that the numerator vanishes on each of `Y_i=0` and `Y_j=0`, so the
division is exact in the integer polynomial ring. The resulting
invariant has row degree four at `i,j` and six at the other rows.

For the complementary triple `C={1,...,5}\{i,j}`, define the signed
triangle, with indices numbered from one, by

```text
T_ij=(-1)^(i+j+1) product_(r<s, r,s in C) Delta_rs.
```

There are twenty integer invariants `S_alpha,beta,gamma`, each of degree four in
every row, for which

```text
H_(ij;alpha,beta,gamma)
       =T_ij S_alpha,beta,gamma(P)       for every i<j. (7)
```

To prove the common factorization, restrict the two derivative rows
to the kernel of evaluation. Their gradient vectors satisfy the three
Veronese equations, so their two-fold exterior product is proportional
to the signed complementary three-by-three Veronese minors `T_ij`.
This proves the common rational quotient in (7) on the open set of
distinct rows. The ten triangle polynomials have polynomial gcd one:
each bracket factor is absent from at least one triangle. Therefore
the common quotient is polynomial. Primitivity of those triangle
polynomials and the integer coefficients of (6) show that it is an
integer polynomial. Its multidegree is `(4,4,4,4,4)`.

For an explicit finite description, normalize the rows to

```text
(0,1), (1,0), (1,1), (1,a), (1,b)
```

and let `F_alpha(a,b)` be the six graph sections. The complementary
triangle for the last two rows is one, and (7) becomes

```text
S_alpha,beta,gamma(a,b)
 =det [ F_alpha            F_beta            F_gamma
        partial_a F_alpha  partial_a F_beta  partial_a F_gamma
        partial_b F_alpha  partial_b F_beta  partial_b F_gamma ]. (8)
```

The two partial derivatives refer to the displayed chart variables.
The checker verifies all two hundred identities (7) as complete
polynomials in these variables, using constant determinant-one frames
for the five normalized rows. Their simultaneous `SL_2` covariance
and row homogeneity, together with the frame-independence of (6),
extend the identities to the full dense open configuration space.
This supplies a second verification of the homogeneous factorization,
in addition to the Veronese exterior-product argument above.

The twenty polynomials in (8) have rank sixteen and lie in the
integer span of the twenty-one products `F_alpha F_beta`, `alpha<=beta`.
Their lattice has **index exactly three** in that product lattice.
Here is a finite integral certificate:

* Sixteen selected products have a polynomial coefficient minor of
  determinant `-1`, and all twenty-one products are integer combinations
  of them.
* All twenty first-jet polynomials are integer combinations of those
  sixteen products.
* A selected sixteen-by-sixteen minor of their coefficient matrix has
  absolute determinant six.
* The full coefficient matrix has rank sixteen modulo two and rank
  fifteen modulo three.

The index divides six by the selected minor. The two modular ranks
remove its factor two and retain its factor three, proving the stated
index. Every step is checked by complete polynomial equality, not
numerical interpolation. It follows at any integer source tuple that

```text
s=gcd_(alpha<beta<gamma) |S_alpha,beta,gamma(P)| = epsilon g^2,
epsilon in {1,3}.                                      (9)
```

Indeed the gcd of all quadratic products of the six evaluations is
`g^2`; the first two lattice inclusions and index three imply both
`g^2|s` and `s|3g^2`. At least one product value is nonzero, so (9)
also proves that the complete first-jet map has rank three.

## 3. The exact image index

Let `a_0` be the gcd of the ten three-by-three minors of the matrix
`A` in (4). Then the index of (5) is

```text
[K:Phi(Lambda)] = s a_0/(g product_i G_i).             (10)
```

Here all the quantities are ordinary positive integers, and the
displayed ratio is an integer. To see the formula directly, perform
unimodular integer column operations on the evaluation row so that
it is `(g,0,0,0,0,0)`. Its last five columns give an integer basis
of `Lambda`. Let `D_i` be the derivative row on those five columns.
For each fixed pair `i,j`, the gcd of the three-by-three minors in
(6) is then `g` times the gcd of the two-by-two minors of the two
rows `D_i,D_j`. By (7) that first gcd is `|T_ij|s`.

Dividing derivative row `i` by `G_i` gives the integer image matrix
of (5). Thus its pair-`i,j` minors have gcd
`|T_ij|s/(gG_iG_j)`. Taking the gcd over the ten row pairs gives
(10), because the complementary minor of `A` is
`T_ij product_(k notin {i,j})G_k`, up to sign.
Unimodular integer reduction identifies this gcd of maximal minors
with the index of the image in its saturated rational span. That span
is exactly `K`.

## 4. All growing core factors cancel

Set

```text
D_0=product_U n_U^max(0,2|U|-5),
L=product_(|U|=3) n_U^3
  product_(|U|=4) n_U^9 n_{12345}^15.
```

The evaluation core divisor is `D_0`, and the cofactor core divisor
is `L`, as in the earlier gradient notes. Write

```text
g=D_0 delta,       a_0=L rho.                          (11)
```

Both `delta` and `rho` are positive integers. The cut exponents give
the exact identity

```text
product_i G_i=D_0 L.                                  (12)
```

For cut sizes one through five, the two sides of (12) have exponents
`0,0,4,12,20`. At equal core weights the logs of `D_0`, `L`, and
`product_i G_i` are respectively `30w`, `90w`, and `120w`; equality
(12) does not require equal or comparable weights.

The extra gcd factors satisfy the divisibilities

```text
delta | product_(i<j) b_ij^2 = B^2,
rho   | product_(i<j) b_ij   = B.                     (13)
```

For `delta`, at any core norm-prime select one of the twenty-two
graphs attaining its minimum internal-edge order. After removing
`D_0`, its remaining valuation comes only from pair factors `b_ij`.
Each graph uses any edge at most twice, which proves the first
divisibility prime by prime. For primes outside all cores choose any
graph. The gcd over six basis graphs is the gcd over all twenty-two,
so restricting to the basis loses nothing.

For `rho`, at a core norm-prime choose a three-column minor attaining
the minimum core order. Its remaining pair factors are the three
edges of a triangle, each used once. The same argument proves the
second divisibility. The factors `b_ij` may overlap core primes or
one another arbitrarily; neither proof assumes additional coprimality.

Combining (9)--(13) proves (1). In particular, at every odd split
core prime that divides none of the actual pair factors `b_ij`, the
map (5) is already saturated. The possible factor three in (9) is
fixed and is not a split core prime.

## 5. What is and is not lifted with a height bound

A finite abelian quotient of order `I` is killed by multiplication
by `I`. Thus for every `ell in K`, there is an integer coefficient
vector `c in Lambda` with

```text
Phi(c)=I ell.                                         (14)
```

The cofactor construction gives an `ell` with all coordinates nonzero
and `max|ell_i|<=5exp(6(1+eta)w+3beta)`. Its scaled gradient in (14)
therefore satisfies

```text
max|I ell_i|<=15exp(6(1+eta)w+33beta).                 (15)
```

All five row multipliers of the resulting single invariant are
nonzero. However, (15) bounds the normalized gradient, **not** its
coefficient vector `c` and not `C_Q`. The preimage of a prescribed
gradient is an affine translate of a rank-three integer kernel of
the full first-jet map. Neither its index nor its covolume bounds the
distance of that affine lattice from the origin in coefficient space.
Excluding a small nonzero coefficient vector still requires a further
argument about this distance or the evaluation lattice itself.

## Verification

[check_five_row_first_jet_lifting.py](check_five_row_first_jet_lifting.py)
verifies the polynomial index-three certificate, all two hundred
complementary-triangle identities as complete chart polynomials,
and 600 additional first-jet minors on three literal full-cut source tuples,
the exact image-index formula, and both divisibilities in (13).
It also explicitly solves the integer lifting equations for three
full-support cofactor stresses and produces a single numerical-zero
coefficient vector in each case. The fixtures include shared row
correction primes. No coefficient-height bound on those lifts and no
unbounded endpoint family are asserted.
