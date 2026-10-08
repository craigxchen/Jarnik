# The positive Cauchy difference has the derivative kernel

The positive even/odd Bezoutian has a useful exact Cauchy-matrix
factorization. Its rational nullspace, however, simplifies to the
ordinary polynomial derivative denominators. This audits a possible
positivity refinement of the [integer cotangent height
target](integer_cotangent_lcm_height_target.md); it does not prove a new
conductor inequality or improve the general point-count bound.

Let `0<X_1<...<X_k` be distinct, put `r=floor(k/2)`, and write

```text
P(z)=product_i(z+X_i)=E(z^2)+z O(z^2),
F(z)=product_i(z-X_i),
O_i=O(X_i^2)=product_(j!=i)(X_i+X_j)>0,
D_i=P'(-X_i)=product_(j!=i)(X_j-X_i),
w_i=D_i/(2 X_i O_i),
C_ij=1/(X_i+X_j),       W=diag(w_i),
V_ia=X_i^(2a),          0<=a<r.
```

Let `B` be the positive definite coefficient Bezoutian defined by
`(E(s)O(t)-E(t)O(s))/(s-t)`, as in the [positive-matrix
audit](integer_cotangent_positive_bezoutian_obstruction.md). Then

```text
C-W = diag(O_i)^(-1) V B V^t diag(O_i)^(-1).          (1)
```

For unequal indices this follows from `E(X_i^2)=X_i O_i`. On the
diagonal, differentiating `P` at `-X_i` gives
`D_i=O_i-2X_i(E'(X_i^2)-X_i O'(X_i^2))`, which gives the same
identity. Thus `C-W` is positive semidefinite of rank exactly `r`.

An explicit basis for its kernel is

```text
h_i^(q)=X_i^(2q)/F'(X_i),       0<=q<k-r.             (2)
```

Indeed, put `G(t)=product_i(t-X_i^2)`. The usual interpolation identity
`sum_i (X_i^2)^a/G'(X_i^2)=0` for `0<=a<=k-2` proves that the vectors
`O_i X_i^(2q)/G'(X_i^2)` lie in the kernel in (1). But
`G'(X_i^2)=F'(X_i)O_i`, giving exactly (2). Their independence follows
from the Vandermonde matrix at the distinct squared roots.

The Cauchy inverse has the exact form

```text
rho_i=2X_i product_(j!=i)(X_i+X_j)/(X_i-X_j)
     =(-1)^(k-1)/w_i,
C^(-1)=diag(rho_i) C diag(rho_i),
W C^(-1) W=C.                                       (3)
```

For example, this follows by partial fractions, or by checking the
residues at the simple poles `-X_i` in the rational interpolation
formula. Consequently `C^(-1/2) W C^(-1/2)` is a symmetric involution.
The second matrix `C+W` is positive semidefinite of complementary rank
`k-r`. No approximation to the spectrum is used.

For any subset `S` on which all `w_i>0`, taking its Cauchy determinant
and cancelling the factors in `product_(i in S)w_i` gives

```text
det(C_S)/product_(i in S)w_i
  = product_(i in S,j not in S)(X_i+X_j)/|X_i-X_j| >= 1.  (4)
```

Every factor in the last expression already exceeds one. Thus this
particular principal-minor inequality supplies no extra arithmetic
condition on the positive roots. The kernel calculation supplies the
previously present `F'(X_i)` denominators. Neither assertion rules out
a different argument using the full matrices and the endpoint height
assumptions together.

Two literal integer cotangent fixtures make the normalizations concrete.
For `L=24` and `X=(3456,3723,7184)`, the scalar Bezoutian is
`[833126395920]`; its primitive matrix is `[1]`. A primitive kernel
vector with third coordinate zero is `(10640,-10907,0)`.
For `L=6` and `X=(1146,1257,1842,4104)`, the matrix is

```text
[[531601657320495744, 26730985716],
 [26730985716,                    8349]].
```

Its content is `3`, and its determinant is
`3723796639619822934000`. The all-edge conductor is
`N=5^3*13*17*29*37*73`. The determinant contains its full `5^3`
factor and none of the other five conductor primes; the lower diagonal
entry of the primitive matrix is `2783`, a `5`-unit. These are exact
fixtures, not asymptotic counterexamples.

The [checker](check_cauchy_bezoutian_derivative_kernel.py) verifies (1)--(4)
over exact rational arithmetic, both kernel dimensions and degree
parities, and all displayed fixture values.
