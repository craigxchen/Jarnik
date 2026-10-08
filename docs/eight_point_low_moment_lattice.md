# The joint low-moment lattice and a remaining radius-height target

The primitive central vectors of overlapping six-point subsets belong to
one integer lattice of relations among the first two circle harmonics.
Its covolume has an exact all-prime formula. On the full fair eight-point
profile its logarithm is `53w+o(w)`, compared with squared-radius logarithm
`127w+o(w)`. This retains the joint lattice, rather than just the heights
of separate central vectors.

No radius bound in terms of this covolume is proved here. A bound with
exponent less than `127/53` on the specified fair endpoint configurations
would suffice for the uniform goal. The exact normalization below makes
that a concrete further target, not a conclusion from lattice integrality.

## 1. One lattice containing all six-point relations

Let `n>=6`, let `z_i` be distinct rational unit-circle points, and define

```text
K=ker_Q [ 1 ; Re z ; Im z ; Re(z^2) ; Im(z^2) ],
Lambda=K intersect Z^n,       d=n-5,
H_Lambda=covolume(Lambda) in its Euclidean span.
```

The matrix has rank five: changing to the complex rows with exponents
`-2,-1,0,1,2` gives a Laurent Vandermonde matrix. Thus `Lambda` is a
saturated rank-`d` integer lattice. The common rotation or reflection of
all circle points changes its row space by an invertible transformation,
so `K`, `Lambda`, and `H_Lambda` are intrinsic to the relative configuration.

Represent the nodes by nonzero integer columns `v_i=(x_i,y_i)` and a
positive rational Gram form `q(v)=v^T Q v`, where after primitive integral
normalization `det Q=r^2`, `r>0` is an integer. Set

```text
q_i=q(v_i),       Delta_ij=det(v_i,v_j).
```

All `q_i` are positive and all `Delta_ij` are nonzero. The same row space
is given by homogeneous quartic evaluation divided by `q_i^2`:

```text
A_(j,i)=x_i^(4-j) y_i^j/q_i^2,       0<=j<=4.          (1)
```

Indeed the numerators of the five circle harmonics are five independent
homogeneous quartics after multiplying by `q_i^2`. The real linear change
of variables defining the Gram form is invertible, so this is a basis of
the space of quartics. In particular `ker A=K` over the rationals.

For every five-subset `J`, the corresponding maximal minor is exactly,
with the consistent determinant ordering,

```text
M_J=product_(i<j in J) Delta_ij / product_(i in J) q_i^2. (2)
```

Let `g` be the positive rational gcd of these nonzero rational numbers:
`v_p(g)=min_J v_p(M_J)`. Then `P_J=M_J/g` are primitive integers.
Complementary-minor duality identifies them, up to the standard signs,
with the primitive Plucker coordinates of `Lambda`. Consequently

```text
H_Lambda^2=sum_(|J|=5) P_J^2.                         (3)
```

For completeness, a basis of a saturated integer lattice extends to an
integer basis of the ambient free module. Its maximal minors therefore
have gcd one, by applying Cauchy--Binet to that unimodular basis and its
inverse. Any rational basis of the same subspace has proportional maximal
minors. Primitive normalization gives (3) by the usual Gram determinant
identity. No chosen unsaturated relation basis is substituted for `Lambda`.

## 2. The exact signed-edge formula at every prime

For a rational prime `p`, put

```text
n_i=v_p(q_i),       d_ij=v_p(Delta_ij),       rho=v_p(r),
lambda_ij=n_i+n_j-2rho-2d_ij,
E_J=sum_(i<j in J) lambda_ij,       |J|=5.
```

Each vertex of the complete graph on `J` has degree four, so (2) gives

```text
v_p(M_J)=sum_(i<j in J)d_ij-2sum_(i in J)n_i
         =-E_J/2-10rho.
```

Taking the minimum and subtracting it yields the exact formula

```text
v_p(P_J)=(max_(|J'|=5) E_J' - E_J)/2.                (4)
```

The half-integer notation creates no parity assumption: the differences
are even because they equal twice differences of valuations of rational
minors. Formula (4) includes primes dividing `2r`, accidental cancellation
in node differences, and the full primitive normalization.

If `delta_ij=|z_i-z_j|`, then

```text
delta_ij^2=4r^2 Delta_ij^2/(q_i q_j).
```

Thus `lambda_ij=-v_p(delta_ij^2)` at odd primes. At two, the latter
quantity is `lambda_ij-2`. Every five-subset energy is shifted by the same
`20`, so (4) is unchanged. This is another fully projective description
using just the squared chord data. The least squared radius, by the
earlier exact formula, satisfies

```text
v_p(N)=max_(i<j) max(0,lambda_ij)   for p odd,
v_2(N)=0.                                             (5)
```

Equations (4)--(5) compare two distinct extrema of the same signed local
edge data. They do not by themselves bound their global heights.

## 3. Explicit Gale presentation and deletion compatibility

Write `D_i=product_(j!=i)Delta_ij`. A rational basis of `K` is the
`n` by `n-5` matrix whose row `i` is

```text
(q_i^2/D_i)
 (x_i^(n-6), x_i^(n-7)y_i, ..., y_i^(n-6)).           (6)
```

Multiplying an entry of (1) by one of these entries gives a homogeneous
polynomial of degree `n-2` divided by `D_i`. Its sum over all nodes is
zero by homogeneous Lagrange interpolation. The columns in (6) are
independent by the Vandermonde rank calculation, so they span `K`.

For eight points, let `u_i` be their primitive central integer vector,
whose raw weights are `q_i^3/D_i`. Then (6) is projectively the matrix
with rows

```text
(u_i/q_i) (x_i^2,x_i y_i,y_i^2).                      (7)
```

Its triple minors are
`(product_(i in I)u_i/q_i) product_(i<j in I)Delta_ij`.
The six-point central vector obtained by deleting labels `a,b`, extended
by zero in those two positions, is projectively

```text
u_i Delta_ia Delta_ib/q_i.                            (8)
```

This is a quadratic evaluation in (7). Hence all such deleted central
vectors belong to this one rank-three lattice after integer normalization.
Their individual contents cannot be discarded, but their rational span is
exactly `K`. For example the pairs `{a,b}`, `{a,c}`, `{b,c}` for three
distinct markings give three independent quadratic forms.

This is stronger bookkeeping than treating the six-subset vectors as
unrelated short vectors. It does not imply that these primitive circuits
form a saturated integer basis.

## 4. Exact height on the full fair profile

Suppose a hypothetical actual configuration has one independent cut block
of logarithmic norm `w+o(w)` for each nontrivial cut, identifying a cut
with its complement. There are `2^(n-1)-1` such blocks. Assume also that
the accumulated errors from private factors and cotangent denominators
are `o(w)`, as in the earlier fixed-cardinality truncation. These are the
full-profile hypotheses; no such endpoint family is constructed here.

On a clean cut of size `s`, `lambda_ij` equals its height when the edge
crosses the cut and zero otherwise. If a five-subset `J` has `a` nodes on
one side, then `E_J=a(5-a)` times that height. Its maximum is four for a
singleton cut and six for every other cut when `n>=6`. There are `n`
singleton cuts. Therefore the sum of the maximal-energy coefficients is

```text
6(2^(n-1)-1)-2n.
```

For a fixed five-subset `J`, each of its ten edges crosses exactly
`2^(n-2)` of the cuts, giving total coefficient `10*2^(n-2)`. Equation
(4) consequently gives, for every primitive coordinate `P_J`,

```text
log |P_J|=(2^(n-2)-n-3)w+o(w),
log H_Lambda=(2^(n-2)-n-3)w+o(w).                    (9)
```

The transfer of the error bound is legitimate: (4) is a fixed finite
maximum of linear forms. Its total variation is bounded by a constant
depending only on `n` times the total valuation error. No prime is
discarded merely because its entire cut height is large at a bad place.

Some coefficients are

| Number of points | Lattice rank | `log H_Lambda / w` | `log N / w` |
|---:|---:|---:|---:|
| 6 | 1 | 7 | 31 |
| 7 | 2 | 22 | 63 |
| 8 | 3 | 53 | 127 |
| 9 | 4 | 116 | 255 |
| 10 | 5 | 243 | 511 |

In particular, at eight points the maximal-energy sum is `746w`, the
energy of a fixed five-subset is `640w`, and half their difference is
`53w`. This calculation concerns the joint saturated lattice.

## 5. Sufficient height comparisons, still unproved

The backwards extraction applies to arbitrarily slowly growing clusters.
Thus an estimate, on the indicated eight-point endpoint profiles,

```text
N <= B H_Lambda^gamma,       gamma < 127/53,           (10)
```

with `B` independent of radius, would contradict (9) and
`log N=127w+o(w)`. It would therefore suffice for the uniform objective
when combined with the already proved extraction and truncation. The
analogous seven-point threshold is `63/22`. Neither comparison is proved.

A stronger universal quadratic comparison would give the endpoint bound
directly. To see this, take actual integer circle points `Z_i` of squared
radius `N`, and use the integer moment matrix

```text
[1 ; Re Z ; Im Z ; Re(Z^2) ; Im(Z^2)].
```

Its five-column determinant has modulus exactly

```text
(N^3/4) product_(i<j in J)delta_ij.                    (11)
```

The factor `1/4` comes from the two changes between real/imaginary row
pairs and complex conjugate row pairs. The complex determinant is the
Laurent Vandermonde. Each raw minor is an integer, and dividing by their
integer gcd can only decrease their Euclidean norm. If the unit-circle
arc has angular length `theta`, then `delta_ij<=theta`, giving

```text
H_Lambda <= sqrt(binomial(n,5)) N^3 theta^10/4.        (12)
```

If one could prove `N<=B_n H_Lambda^2` for some fixed `n`, (12) would
imply

```text
theta >= [16/(B_n binomial(n,5))]^(1/20) N^(-1/4).
```

That implication is valid; the proposed height comparison is not known
here. It is false at six points because the previously constructed
fixed-central elliptic family has `H_Lambda=sqrt(102)` fixed and
unbounded `N`. In particular one cannot infer the comparison simply from
the existence of an integer relation lattice.

Even divisibility by `H_Lambda^2` is false at eight points. For the
cotangents `0,1,...,7`, exact calculation gives

```text
N=1022125,
H_Lambda^2=1053718749231052896,
H_Lambda^2 mod N=879896.
```

This does not disprove an inequality `N<=B H_Lambda^2`; it rules out
that particular divisibility shortcut. The remaining work is an actual
radius-height theorem using circle realizability, or a counterexample
to delimit its correct scope.

## 6. Seven points: the relation lattice determines the circle metric

For seven points the relation space has dimension two. Choose a rational
basis matrix `B`, with nonzero rows `b_i=(x_i,y_i)` in seven distinct
projective directions, and put
`D_i=product_(j!=i) det(b_i,b_j)`. These directions are distinct by the
Gale presentation above. Express the original circle metric as a positive
rational quadratic form `q` in these projective coordinates.

The Gale formula says that `B` and
`diag(q(b_i)^2/D_i) B` have the same column space. Consequently a
two-dimensional linear operator has the seven row directions as
eigenlines. An operator with three distinct eigenlines is scalar, so

```text
q(b_i)^2/D_i = a common nonzero constant.             (13)
```

In particular all the `D_i` have the same sign, and every ratio
`D_i/D_1` is a positive rational square. Normalize the metric by
`q(b_1)=1`, and set `c_i=sqrt(D_i/D_1)>0`. Its three coefficients
are determined uniquely by `q(b_i)=c_i` on any three directions:
the quadratic evaluation determinant is, up to sign, the product of
the three nonzero pair brackets. The recovered form fits the other
four rows and has positive rational-square determinant because it is
the original metric up to a positive rational factor. Every squared
chord is therefore recovered by

```text
delta_ij^2 = 4 det(q) det(b_i,b_j)^2/(c_i c_j).        (14)
```

Thus the seven-point relation space alone determines the intrinsic
squared radius `N`, as the least common multiple of the reduced
squared-chord denominators.

For an integer basis of the saturated lattice, let
`p_ij=det(b_i,b_j)`, extend these antisymmetrically, and set
`S_i=product_(j!=i) p_ij` and `d_0=gcd_i |S_i|`.
Equation (13) implies that all `v_p(S_i)` have the same parity.
Subtracting their minimum at each prime proves

```text
|S_i|/d_0 = c_i^2,    c_i positive integers, gcd_i c_i=1. (15)
```

The metric can instead be normalized to these integer values.
In fact all seven `S_i` are negative: their product is
`(-1)^21 product_(i<j)p_ij^2`, and they have one common sign.
Moreover `d_0` is itself a square, since
`product_(i<j)p_ij^2=d_0^7 product_i c_i^2`.
Thus every `-S_i` is an integer square. These stronger necessary
conditions still do not imply the quadratic-fitting conditions below.
The elementary bound `c_i<=H_Lambda^3` follows because each of the
six brackets is at most `H_Lambda`. This does not yet yield the
required exponent below `63/22` for the radius.

There is also an exact converse test. A rank-two rational matrix with
seven distinct projective row directions is realized by such a circle
relation space if its star ratios are positive rational squares and
the quadratic form interpolating their positive square roots fits
**all seven** rows, is positive definite, and has rational-square
determinant. The Gale formula then recovers the given column space.
The star-square condition by itself does not imply those additional
quadratic conditions.

## Verification

[check_low_moment_plucker_independent.py](check_low_moment_plucker_independent.py)
checks direct rational nullspaces, quartic evaluation minors, weighted Gale
minors, the eight-point central presentation, every prime including two,
and the fair-cut enumeration independently. The root reran 18 actual
seven/eight-point configurations and 139 prime checks, including nine
seven-point metric reconstructions. The formal fair
enumeration is kept separate from these actual circle examples.
