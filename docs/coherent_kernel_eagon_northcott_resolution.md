# A global resolution of the coherent kernel

The determinantal description gives an explicit resolution of the
coherent kernel at every distinct configuration. It also gives an
independent formula for its rank and the degree of its primitive
Pluecker numerator. In particular, the twelve-row kernel is generated
by the localized six-row relations everywhere on the distinct locus.
This is algebraic generation over the ground field. No bounded-height
integral lifting is asserted.

Work in characteristic zero. Put `m=nd`, `n,d>=2`, and `e=d-2`.
Let `A` be the space of multilinear `SL_n` invariants of weight
`det^d` on `m` vector rows. At distinct binary directions
`p_i=(P_i:Y_i)`, let `K_p` be the coherent kernel from the
[determinantal generation theorem](distinct_configuration_determinantal_kernel_generation.md).
Its generators are maximal minors of

```text
M_p(V)=[(P_i v_i^a) | (Y_i v_i^a)]_(i=1..m,a=1..n),
```

multiplied by invariant polynomials on the remaining rows.

## The terms of the resolution

Let `f^lambda` denote the number of standard tableaux of the partition
`lambda`, with `f^empty=1`. Define, for `0<=j<=2e`,

```text
a_j = sum_(a+b=j, e>=a>=b>=0)
       (a-b+1) f^(e,...,e,e-b,e-a),                    (1)
```

where the partition has `n-2` initial parts equal to `e`. Set
`a_j=0` outside this range. On `X=(P^1)^m`, write
`O(-1_S)` for the line bundle with degree minus one at labels
in `S` and zero at other labels, and put

```text
F_j = direct-sum_(|S|=2n+j) O(-1_S)^(a_j).           (2)
```

On the distinct locus there is an exact complex of vector bundles

```text
0 -> F_(2e) -> ... -> F_1 -> F_0 -> K -> 0.           (3)
```

Zero terms at either end are omitted. The map `F_0 -> A` is the
localized determinant construction. All preceding differentials come
from the Eagon--Northcott differential and are linear in the binary
row coordinates. Each term selects a row set `S`, so the binary
line-bundle shift in (2) is one at every selected row.

Here is the representation calculation. Use the convention in which
the degree-one vector-coordinate module is the standard `n`-dimensional
module `U`. The maximal minors have weight `det^2`. The `j`th
Eagon--Northcott generator module has the additional factor
`Sym^j(U tensor k^2)`, and its complementary `m-2n-j=ne-j`
row variables contribute `U^(tensor(ne-j))`. To extract weight
`det^d`, it is therefore enough to find the multiplicity of
`det^e` in their tensor product.

The Cauchy decomposition gives

```text
Sym^j(U tensor k^2)
 = direct-sum_(a+b=j,a>=b>=0) S_(a,b)U tensor S_(a,b)k^2.
```

The second factor has dimension `a-b+1`. The dual of the first
factor, twisted by `det^e`, has highest weight
`(e,...,e,e-b,e-a)`. It is polynomial exactly when `a<=e`.
Schur--Weyl duality says that its multiplicity in the remaining
tensor power is `f^(e,...,e,e-b,e-a)`. This proves (1) and the
vanishing of all invariant terms with `j>2e`. Extracting one row
multidegree is exact; taking the `det^d` isotypical part is exact
in characteristic zero. Thus the Eagon--Northcott resolution of the
maximal-minor ideal gives (3).

For the imported algebra, the complex and expected-grade acyclicity
are recorded in [Kustin, Sections 1 and 3](https://arxiv.org/pdf/1509.04367).
The [Cauchy decomposition in Pragacz, Section 2](https://www.impan.pl/~pragacz/ub.pdf)
and [Schur--Weyl decomposition in Snowden, Lecture 6](https://websites.umich.edu/~asnowden/msri19/course.pdf)
fix the representation conventions. The expected height of this
particular specialized matrix is proved in the linked determinantal
generation note; it is not assumed from the generic-matrix case.

## Rank and the primitive determinant degree

Taking Euler characteristics in (3) gives

```text
h(n,d) = sum_j (-1)^j binom(m,2n+j) a_j.             (4)
```

There is also an extension across every codimension-one collision,
which is needed before reading a primitive numerator degree from
the determinant line. Let `U_n subset X` consist of configurations
in which no direction occurs more than `n` times. The complement
has codimension `n>=2`.

The maximal-minor ideal still has expected height at every point of
`U_n`, although it need not be prime there. To see this, use an
affine chart for the finitely many distinct binary directions. For
`[u:v] in P^(2n-1)`, a rank defect imposes the row equations
`v_i dot (u+z_i v)=0`. With independent `u,v`, all row normals
are nonzero, and the incidence dimension is at most

```text
D=m(n-1)+2n-1.
```

For a dependent pair that makes a normal zero on a group of `r`
equal directions, the incidence dimension is
`m(n-1)+n-1+r<=D`, because `r<=n`. The dependent pairs with
no zero normal have smaller dimension as well. Projection therefore
bounds the rank-defect variety by dimension `D`. The general
maximal-minor height bound supplies the reverse inequality. Thus
the height is `m-2n+1`, as required by Eagon--Northcott.

Extracting the same multidegree and isotypical component now proves
fiberwise exactness of (3) on `U_n`, with its last term interpreted
as the image of `F_0 -> A`. The ranks of all differentials are
constant by induction from the left, so their images are vector
subbundles. In particular this extends `K` as a rank-`h` subbundle
of the constant source bundle on `U_n`.

The determinant of (3) has degree `-delta` in every binary row,
where

```text
delta(n,d)=sum_j (-1)^j binom(m-1,2n+j-1) a_j.       (5)
```

This sign agrees with `K` being a subbundle of a constant bundle.
Its determinant inclusion gives a tuple of sections of
`O(delta,...,delta)` on `U_n`. Normality of `X` and the
codimension-two complement extend those sections uniquely to `X`.
At least one coordinate is a unit in each local trivialization on
`U_n`, so the extended tuple has no common divisor: every divisor
of `X` meets `U_n`. Thus it is the primitive polynomial Pluecker
tuple, up to one scalar. Formula (5) computes its actual row degree,
with no hidden pair-diagonal factor. Clearing one constant
denominator and content gives the primitive integral normalization.

## The twelve-row resolution

At `(n,d)=(3,4)`, formula (1) gives

```text
(a_0,...,a_4)=(5,10,9,4,1),
(rank F_0,...,rank F_4)=(4620,7920,4455,880,66).
```

Consequently

```text
h=4620-7920+4455-880+66=341,
delta=2310-4620+2970-660+55=55.
```

This independently recovers the coherent-kernel rank and primitive
row degree previously obtained by joint differentials and fusion.
It also specifies the dimensions of all higher relations among the
localized generators, rather than just their number and image rank.

## Arithmetic scope

The rows in an extracted core coalition can have multiplicity larger
than `n` after reduction at their Gaussian prime. The extension
argument above intentionally does not cover those special fibers.
Moreover, exactness over a characteristic-zero field does not bound
the denominators or heights needed to lift one integer kernel vector
through (3). Those are precisely the arithmetic issues needed for
a new first-minimum estimate. The resolution supplies a concrete
presentation for that investigation, not the missing estimate.

The later [linked determinant theorem](linked_determinantal_coalition_saturation.md)
and [global primitive-generator index bound](global_primitive_local_generator_index.md)
control local saturation and give subpower global index for the
extracted profile. They do not bound the Archimedean sizes of the
coefficients in a lift, so the first-minimum question remains open.

The [exact checker](check_coherent_kernel_eagon_northcott_resolution.py)
checks tableau dimensions independently, compares (4)--(5) against
the fusion recurrence, and verifies explicit linear Laplace syzygies.
The all-configuration result is proved above and in the linked
generation theorem, not by the finite tests.
