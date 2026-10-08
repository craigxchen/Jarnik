# Redundant generators do not force a short relation

This exact lattice model tests what follows from the numerical data
`ambient rank 462`, `kernel rank 341`, `4620` generators, generator
height exponent `1152`, saturated image, and determinant exponent
`171600`. It shows that these data, including substantial redundancy,
do not force a vector below the average exponent `171600/341`.

## A saturated rank-341 lattice in the correct ambient dimension

Put `h=341`, and for an integer `q>=2` define

```text
K_q = {x in Z^(h+1) : x_0 + q*x_1 + ... + q^h*x_h = 0}.
```

Embed it in `Z^462` by appending 120 zero coordinates. It is saturated:
it is the kernel of a surjective homomorphism to `Z`, since the first
coefficient of the defining form is 1. A basis is

```text
b_i = q*e_i - e_(i+1),  0<=i<h.
```

The signed maximal minors of this basis matrix are
`(1,q,q^2,...,q^h)`, so

```text
covol(K_q) = sqrt(1+q^2+...+q^(2h)) = q^h * Theta(1).
```

Every nonzero vector in `K_q` has some coordinate of absolute value at
least `q`. Reduce the defining equation modulo `q`: if every coordinate
had absolute value less than `q`, then `x_0=0`; divide by `q` and repeat
to obtain `x_1=...=x_h=0`. Since `||b_i||=sqrt(q^2+1)`, this gives the
sharp first-minimum scale

```text
q <= lambda_1(K_q) <= sqrt(q^2+1).
```

Choose `q=q_w` with
`log q_w = alpha*w+o(w)`, where
`alpha=171600/341=503.2258...`. Then `K_q` has determinant exponent
`171600` and first-minimum exponent `503.2258...`, well above 90.

## 4620 primitive generators, all at height exponent 1152

Let `U_N` be a unimodular `h`-by-`h` integer matrix whose columns all
have norm `Theta(N^2)`. One explicit construction starts with the block
diagonal matrix consisting of copies of

```text
[[N, N-1], [N+1, N]],
```

whose determinant is 1, and a final `[1]` if `h` is odd. Add `N` times
column 1 to the last column, then add that last column to each other
column. These are unimodular column operations. All entries involved
are nonnegative for `N>=2`, so every resulting column has norm
`Theta(N^2)`.

Map the columns `z` of `U_N` to `Bz`, where `B` has columns `b_i`.
The triangle inequality gives

```text
(q-1)||z|| <= ||Bz|| <= (q+1)||z||.
```

Every `Bz` is primitive in the ambient integer lattice when `z` is a
column of a unimodular matrix: if a prime divided all coordinates of
`Bz`, backward substitution from its last coordinate would make it
divide every coordinate of `z`. Since `U_N` is unimodular, these `h`
vectors generate all of `K_q`.

Now add the `h` columns of `B C_R`, where `C_R=R I+P` and `P` is the
cyclic permutation matrix. Each column of `C_R` is primitive, has norm
`Theta(R)`, and
`det C_R=R^h+1` because `h=341` is odd. These extra columns leave the
generated lattice equal to `K_q`, while one `h`-minor in the first `h`
coordinates has absolute value `q^h(R^h+1)`. Finally append
`4620-2h=3938` copies of a column of `U_N`. The resulting 4620 columns
are all primitive, all have norm `Theta(qN^2)`, span `K_q`, and have
integer relation rank `4620-h=4279`.

Set `R=N^2` and choose

```text
log q = (171600/341)w+o(w),
log R = 2 log N = (1152-171600/341)w+o(w).
```

Then every generator has logarithmic height `1152w+o(w)`, but the
saturated generated lattice still has `lambda_1=exp(503.2258...w+o(w))`.
Thus saturation, generator count, and individual heights do not force a
relation of height `90w`.

The full real generator matrix is well conditioned at this scale as
well. The singular values of `C_R` lie between `R-1` and `R+1`,
and those of `B` lie between `q-1` and `q+1`. Its `B C_R` block
therefore gives a least nonzero singular value of order `qR`;
the fixed number of columns and their norm bounds give the matching
upper bound. Hence all `h` nonzero singular values of the complete
generator matrix are `Theta(qR)`. Uniform Archimedean conditioning
of the normalized full family does not remove the obstruction.

## Why the maximal minor is not the lattice covolume

In coordinates of the basis `B` of `K_q`, the generator matrix contains
the blocks `U_N` and `C_R`. Its maximal-minor gcd is 1 because
`det U_N=1`; this gcd is the index of the generated lattice in `K_q`.
At the same time, the `C_R` block has an archimedean maximal minor of
size `R^h+1`. The ambient minor gains the additional factor `q^h`.
The image remains `K_q`, whose Euclidean covolume is only
`q^h Theta(1)`. A large maximal minor records a large parallelepiped
from one selected tuple of columns; it is not the gcd of all maximal
minors and does not replace the lattice covolume.

For a general full-rank column map `A: Z^N -> L` into an `h`-dimensional
Euclidean space, the exact coarea identity is

```text
covol(L) * covol(ker(A) intersect Z^N) = sqrt(det(A A^T)).
```

The Eagon--Northcott ranks `(4620,7920,4455,880,66)` alternate to 341
and describe ranks of modules in a polynomial resolution. They do not
bound either factor on the left after specialization. The integer
syzygy group of a single specialized map is a free abelian lattice of
rank 4279; its covolume and orientation in `Z^4620` are separate metric
data. Using it would require an independent height or covolume estimate
for that embedded lattice, as well as control of the index of the
local-generator image. Exact local saturation sets that index to one
but still leaves the above model possible.

Even exactness with those numerical ranks alone adds no restriction:
the successive image ranks after the first map can be
`4279,3641,814,66`. Since the first integer kernel is free, choose
a basis for it and extend by split free summands to obtain an exact
free resolution with precisely the ranks
`(4620,7920,4455,880,66)`. This abstract resolution need not have the
specific polynomial differentials or row gradings of the actual
Eagon--Northcott complex. It isolates why those additional structures,
rather than only exactness and ranks, would have to supply a height gain.

This is an abstract obstruction, not a construction from the actual
circle-point generators. Actual coalition identities could impose
additional restrictions on the integer matrix that exclude this model;
those restrictions would need to enter a new proof.
