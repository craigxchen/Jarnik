# Independent audit of the seven-row middle radical determinant

This note independently checks
`seven_radical_hankel_content.md` and the squareclass identification in
`gale_radical_squareclass_identification.md`.  The determinant constant,
the complete rational content, the integer base, and the covariance factors
are correct.  No new height saving follows from the identity.

## 1. Affine determinant and its sign

For seven distinct real nodes, let

```text
P(t)=product_i(t-t_i),       y_i^2=1+t_i^2,
M_(a,b)=sum_i y_i t_i^(a+b)/P'(t_i),   0<=a,b<=2.
```

Cauchy--Binet, or equivalently a Laplace expansion of the evaluation
matrix with columns `1,t,t^2,t^3,y,ty,t^2y`, gives

```text
det M=-det E/V(t),       V(t)=product_(i<j)(t_j-t_i).
```

With `u_i=y_i+t_i`, one has

```text
t_i=(u_i-u_i^(-1))/2,       y_i=(u_i+u_i^(-1))/2.
```

Changing the seven displayed functions to the Laurent basis
`u^(-3),...,u^3` has determinant `1/512`.  Alternation then gives

```text
det E=(1/512) V(u)/product_i u_i^3.
```

Finally

```text
t_j-t_i=(u_j-u_i)(1+u_i u_j)/(2u_i u_j),
```

so the sign and power of two are exactly

```text
det M=-2^12 product_i u_i^3/product_(i<j)(1+u_i u_j).  (1)
```

The independent checker evaluates (1) over `Q` on several half-angle
tuples and on all 128 branch choices.  It also verifies directly that their
product is `2^192/V(t)^64`.  The denominator never vanishes: under the
conic parametrization `1+u_i u_j=0` implies `t_i=t_j`.

## 2. Complete content and the integer base

For a triple `I`, direct substitution of

```text
c_i=Q_i^5/(D_i^2 tau)
```

in the squared Cauchy--Binet term gives

```text
w_I^2=(product_(i in I)Q_i) V_I^2 V_(I^c)^2
      /(tau^3 D^2),                                  (2)
```

where `D=product_(i<j)|d_ij|`.  Therefore, with

```text
G=gcd_Q((product_(i in I)Q_i)V_I^2V_(I^c)^2),
```

one has exactly

```text
g_3=G/(tau^3D^2),
B_0=8 det(q)^(3/2)D/G.                               (3)
```

This calculation retains primes in denominators as well as numerator
primes.  Moreover `w_I^2/g_3` is a nonnegative rational integer for every
triple.  Hence each signed `w_I/sqrt(g_3)` satisfies a monic quadratic
over `Z`, and every formal signed value of `Z_3` is an algebraic integer.
This proves the square-root-content step even when the radicals have field
relations.

The product over all formal signs is rational and equals `B_0^64`.
It is a product of algebraic integers, so it is a rational integer.  Since
`B_0` is positive rational, prime valuations then show `B_0` itself is a
positive integer.  This argument does not require all formal signs to be
Galois automorphisms.  Accordingly, the product is an actual minimal-field
norm only when the squareclass relations permit that interpretation; the
source note correctly keeps these notions separate.

The separate normalization of the individual shared-coefficient forms is
also exact.  For fixed Veronese degree `r`, the square of a term in the
form of numerator degree `j` is

```text
c_i x_i^(2(2r-j)) y_i^(2j)/Q_i^(2r),   0<=j<=2r.
```

For a primitive marked vector, at least one of `x_i,y_i` is a unit at each
prime, and both pure powers occur.  Taking the primewise minimum over all
labels and all `j` gives precisely

```text
e_r=gcd_Q(c_i/Q_i^(2r):i=1,...,7).
```

At `r=3`, the seven forms have a Vandermonde coefficient matrix and the
Lagrange inverse

```text
sigma_i sqrt(c_i)=Q(t_i)^3 sum_j ell_(i,j)L_j.
```

Its leading coefficient `ell_(i,6)=1/P'(t_i)` contains all six small gaps.
The resulting `exp(96w+o(w))` inverse size restores the radical scale from
the small raw form scale, so these seven forms are not uniformly
conditioned.  The checker verifies both claims over `Q`.

For

```text
U_ij=(q(v_i,v_j)+sqrt(Q_iQ_j))/(sqrt(det(q)) |d_ij|),
```

the binary Gram identity gives

```text
q(v_i,v_j)^2-Q_iQ_j=-det(q)d_ij^2.
```

Thus the formal quadratic conjugate of `U_ij` is `-U_ij^(-1)`, and the
covariant form of (1), combined with (2)--(3), gives

```text
Z_3^2=B_0/product_(i<j)U_ij.                         (4)
```

The factors need not be algebraic integers, so (4) is a factorization in
the compositum rather than a factorization into global units.

## 3. Metric and basis covariance

The formulas survive all normalization choices used in the Gale
construction.

* Scaling `q` by a positive rational `lambda` scales `Q_i` and its bilinear
  values by `lambda`, and scales `sqrt(det(q))` by `lambda`.  It scales `G`
  by `lambda^3`.  Thus both `B_0` and every `U_ij` are fixed.  The common
  scaling of the raw `rho_i^2` is absorbed by `tau`, so the primitive `c_i`
  are fixed as well.
* Under a rational basis change `v_i -> T v_i` with the quadratic form
  transported contragrediently, `Q_i` and its bilinear values are fixed,
  while `|d_ij|` scales by `|det T|` and `sqrt(det(q))` by
  `|det T|^(-1)`.  Hence every `U_ij` is fixed.  In (3), `D/G` scales by
  `|det T|^3`, which is canceled by `det(q)^(3/2)`.
* In a saturated Gale basis, the square-star relations give
  `D=s^7 product_i c_i` and `tau=s^(-4)`.  Substitution into (3) yields
  exactly the alternate formula
  `8 det(q)^(3/2)s^5/(product_i c_i g_3)`.

The checker tests a determinant-six diagonal basis change, a shear, and a
rational metric scaling, in addition to random integral marked vectors.

## 4. Squareclasses and field relations

For

```text
ell(x,y)=a x+b y+i sqrt(det(q))y,
beta_ij=ell(v_i) conjugate(ell(v_j)),
```

one has `N(beta_ij)=a^2Q_iQ_j`.  If an integral Gaussian `h_ij` represents
the same full phase ratio, equality of `beta/conjugate(beta)` and
`h/conjugate(h)` forces `beta/h` to be rational.  Therefore

```text
[N(h_ij)]=[Q_iQ_j]=[c_ic_j] in Q*/Q*^2.              (5)
```

The second equality follows exactly from

```text
c_ic_j/(Q_iQ_j)=(Q_i^2Q_j^2/(D_iD_j tau))^2.
```

This includes the unit contribution: the phase `i` is represented by
`1+i`, whose norm contributes the class of 2.  The independent checker
clears arbitrary rational contents in `beta_ij` and verifies (5) pairwise.

If the seven `c_i` have distinct squareclasses, their square roots occupy
distinct character spaces and are linearly independent.  Therefore the
stabilizer of `(sum a_i sqrt(c_i))^2`, for all `a_i` nonzero, consists
exactly of sign actions that give all seven radicals the same sign.  Its
fixed field is generated by the six pair classes `[c_ic_0]`.  Invoking the
already proved six-wise separation theorem makes these six classes
independent and gives degree 64.  This last conclusion is correctly scoped
to the extracted common-unit endpoint configurations satisfying that
theorem; it is not a statement about arbitrary seven-point data.

## 5. Result of the audit

No missing constant, sign, finite-place factor, basis factor, or hidden
independence assumption was found.  The identity explains the zero total
exponent of the middle determinant: its universal formal sign product is
the 64th power of the subpower integer `B_0`.  The explicit collision family
where `v_p(B_0)=3e` also shows that this integer is not absolutely bounded.
Thus the factorization is exact structural information, but it does not
give an additional uniformity estimate.

The bounded independent computation is in
`check_seven_radical_hankel_content_independent.py`.
