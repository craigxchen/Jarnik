# Higher-rank terminal invariants and a finite height comparison

Let `n>=2`, `d>=2`, and `m=n*d`. Take the multilinear `SL_n`
invariant space on `m` vector rows. It is the rectangular Specht module
`S^(d^n)` and has dimension

```text
r_(n,d)=(nd)! product_(i=0)^(n-1) i! / product_(i=0)^(n-1)(d+i)!.
```

The binary Veronese substitution

```text
(P_i,Y_i) |-> (P_i^(n-1),P_i^(n-2)Y_i,...,Y_i^(n-1))
```

is injective even on the **whole** multilinear polynomial space.
Indeed, a coordinate-assignment monomial chooses one index
`j_i in {0,...,n-1}` for every row `i`; its image is the distinct
binary monomial `product_i P_i^(n-1-j_i)Y_i^(j_i)`. No two assignments
have the same image, so coefficients are preserved one to one.
In particular, a product of `d` disjoint `n`-row determinants maps to
a nonzero product of binary Vandermonde bracket products.

## 1. Exact generalized coherent annihilator

Write the ternary identity of
[the `SL_3` note](terminal_coherent_conformal_block_identity.md) in
`n` coordinates. For a multilinear invariant `Q` in rows
`(a_i^0,...,a_i^(n-1))`, set

```text
D_z=sum_i z_i a_i^0 partial_(a_i^(n-1)),
R_z Q=Q((1,z_i,x_i^2,...,x_i^(n-1))_i).
```

Every invariant has total degree `d` in each coordinate family, so

```text
D_z^d Q=d! Q((a_i^0,a_i^1,...,a_i^(n-2),z_i a_i^0)_i).
```

For `n>=3`, swap coordinate families `1` and `n-1` globally. This
matrix has determinant `-1`; the invariant has determinant weight
`d`. Factoring `a_i^0` out of each row gives the polynomial identity

```text
D_z^d Q=(-1)^d d! (product_i a_i^0)
 R_z Q(a_i^2/a_i^0,...,a_i^(n-2)/a_i^0,a_i^1/a_i^0).  (1)
```

For `n=2`, there is no coordinate swap or free `x` family and the
same top derivative is `d!(product_i a_i^0)Q(1,z_i)`. Thus in every
rank the generalized coherent kernel is exactly the annihilator of
the `d`-th power of the weighted highest-root operator. The
[standard genus-zero conformal-block quotient criterion](https://math.umd.edu/~pbrosnan/Papers/BBMArxiv.pdf)
identifies it, with the dual/opposite-root convention in the linked
`SL_3` note, as the level `ell=d-1` conformal-block functional space
with all marks carrying the defining fundamental weight. Equation
(1) is valid for every `z`; distinctness is needed only to regard
the marks as a smooth pointed curve.

The coefficient of any fixed free-coordinate monomial in `R_z Q`
uses `d` labels in each of the `n-2` free families. The remaining
`2d` labels form an integral multilinear binary `SL_2` invariant.
For one choice of the two binary coordinate families, the earlier
projected-block argument gives the floor

```text
log C_Q >= 2^((n-2)d+1)*(1-eta)*w-o(w).              (2)
```

The first `2^((n-2)d)` comes from blocks with a chosen outside
`P`-support; the equal-magnitude complementary `Y`-support supplies
the second coprime conductor. All coordinate pairs can be used on
**the same coefficient**. Fix a nonzero coordinate-assignment
coefficient of `Q`, partitioning the labels into `n` color blocks
`B_0,...,B_(n-1)` of size `d`. For every ordered pair of distinct
colors `(a,b)`, a global Weyl permutation identifies the highest-root
annihilator with the coherent map using colors `a,b` as its binary
coordinates. All other color blocks are fixed by a free-coordinate
monomial. The same coefficient, up to a fixed sign, receives the core
conductor from every `T` satisfying `B_a subset T` and
`T intersection B_b=empty`. Taking the union over ordered pairs,
each core block is counted exactly once if it contains a full color
block and excludes another. Inclusion-exclusion gives the stronger
support count

```text
U(n,d)=2^(nd)-2*(2^d-1)^n+(2^d-2)^n,
log C_Q >= U(n,d)*(1-eta)*w-o(w).                      (2a)
```

The coefficient-level conductor argument, including the ordered-pair
union and row-correction bookkeeping, is proved in
[the Weyl-union floor note](higher_rank_weyl_union_floor.md) under the
same primitive Gaussian core hypotheses as the three-coordinate
actual-cut theorem. The finite checker here compares exponents; it
does not verify those core hypotheses for actual circle configurations.
The finitely many correction losses from the coordinate pairs remain
`o(w)` for fixed `(n,d)`.

## 2. Fusion and determinant exponents

Let `h` be the conformal-block rank. For a level-`ell` dominant
weight in Dynkin coordinates `a=(a_1,...,a_(n-1))`, the minuscule
fundamental fusion step adds one of the `n` vectors

```text
(+1,0,...,0),
(-1,+1,0,...,0), ..., (0,...,-1,+1),
(0,...,0,-1),
```

retaining exactly the states with all coordinates nonnegative and
`sum a_i<=ell`. This is the type-A alcove Pieri rule of
[Andersen–Stroppel, equation (2.2) and Section 2.2](https://people.math.uni-bonn.de/stroppel/Fusion.pdf).
Let `f_s(a)` count paths of length `s` from zero; then `h=f_m(0)`.
The dual weight reverses Dynkin coordinates. With `c(a)` the Casimir
and `c_omega=(n^2-1)/n`, define

```text
tau_s=[sum_a c(a) f_s(a) f_(m-s)(a*)-s*h*c_omega]
      /[2*(d+n-1)].                                  (3)
```

The stable-quotient and KZ argument of
[the `SL_3` fusion-content note](terminal_kernel_fusion_content.md)
uses only the highest-root identification, factorization, and the
Casimir identity, all valid for `SL_n`. Thus the primitive kernel
Pluecker numerator has coalition orders

```text
nu_s=tau_s-binom(s,2)*tau_2       (2<=s<=m/2).       (4)
```

The pair-subtracted determinant connection has no other finite
poles and no regular polynomial one-form by dilation invariance.
For the smallest case `m=4`, a stable pair collision leaves only two
labels outside, so the three-outside-label chart in Fakhruddin
Section 3.1.4 cannot be used verbatim. The unfixed affine KZ
pair-sum calculation supplies the same residue: in the Casimir sign
convention of Fakhruddin's equation (3.4), the quotient residue on
a coalition channel `a` is
`(c(a)-s*c_omega)/[2*(d+n-1)]`; the dual determinant trace is its
negative. This follows from
`2*sum_(i<j in T)Omega_ij=C_T-s*c_omega`, with the KZ sign
calibrated to equation (3.4). It applies also when only two labels
remain outside. Different publications use opposite signs for the
Casimir or for writing horizontal KZ equations, so the calibration is
part of this statement.
The affine Euler vector field and total Casimir give the common
homogeneous binary-row degree

```text
delta= -h*c_omega/(d+n-1) -(m-1)*tau_2.             (5)
```

For `s>m/2`, binary-coordinate interchange gives

```text
nu_s-nu_(m-s)=delta*(s-m/2).                          (6)
```

When `m` is odd, `m/2` in (6) is the literal half-integer; the exact
checker uses rational arithmetic and confirms that every resulting
order is integral on its finite grid.

These orders are termwise support floors for the polynomial
numerator, hence give universal Gaussian core divisibility under
the full-profile hypotheses. They do not determine any additional
specialized gcd.

The resulting covolume and first-minimum exponents are

```text
B=m*delta*2^(m-3)-sum_(s=0)^m binom(m,s)*nu_s,
B/h.                                                   (7)
```

The exact [finite checker](check_higher_rank_terminal_fusion_grid.py)
computes (3)–(7) for `2<=n<=8`, `2<=d<=12`, using integer path
counts and rational Casimirs. Every one of the 77 cases has

```text
B/h >= U(n,d) >= 2^((n-2)d+1),
```

with equality in the first comparison only at `(n,d)=(2,2)` and
`(3,2)`; all others are strict in the opposite direction from the
desired contradiction. Selected ratios `(B/h)/U(n,d)` are:

| `(n,d)` | ratio |
|:---:|---:|
| `(2,2)` | `1` |
| `(3,2)` | `1` |
| `(4,2)` | `1.054545...` |
| `(3,4)` | `5.591397...` |
| `(4,3)` | `2.403697...` |
| `(6,2)` | `1.247964...` |

The checker also verifies, on pairs inside its grid, the critical
rank-duality identities `h_(n,d)+h_(d,n)=r_(n,d)` and the symmetry of
`delta` and `B` under `(n,d)` interchange. These are finite checks,
not proofs of those identities in general. No all-parameter inequality
or radius-uniform arc bound is inferred from the grid.
