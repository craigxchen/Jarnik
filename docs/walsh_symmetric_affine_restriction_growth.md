# Symmetric affine assignments tolerate exceptional rows and arbitrary prime sizes

The symplectic and ordinary-dot phase bounds extend to symmetric affine
Walsh assignments after a controlled loss of rows. The resulting growth
bound uses the original radius and imposes no comparison between prime
sizes. It is a theorem within this assignment family; it does not supply
an extraction of such assignments from arbitrary circle configurations.

Use `M=2^t` Walsh rows and `b>=5` physical copies of each nonzero label,
each attached to a distinct split prime. Every row flips one distinct
physical column with nonzero assigned label `a_x`. Allow arbitrary prime
orientations, row units and a common Gaussian factor. Suppose that outside
an exceptional set of at most `K` rows,

```text
a_x=Lx+v,       L=L^transpose,       rank L>=t-k.              (1)
```

Here `k>=0` is an integer. The original physical flip capacity is part of
the hypothesis: there must be enough copies of each assigned label, and
all flipped columns must be distinct. Any zero labels of `Lx+v` must be
among the exceptions, since the physical model has no zero Walsh label.
No bound on the actual prime weights is assumed.

If the actual rows lie on an arc of length at most `C sqrt(R)`, put
`W=log R^2`. There is a constant `A_C` depending only on `C` such that

```text
M <= A_C 2^k (K+1)^2 [1+(W/log(2+W))^(2/3)].                 (2)
```

In particular fixed `k,K` give the growth exponent `2/3`. The input phase
theorems are the [symplectic averaging proof](walsh_symplectic_weighted_phase_growth.md)
and [ordinary-dot parity proof](walsh_ordinary_dot_weighted_phase_growth.md).
The argument below is an exact restriction of the physical source, not
a replacement of its primes or a reduction to a smaller uncharged radius.

## Finding a nondegenerate affine restriction

Let `q=ceil(log_2(K+1))`, with `q=0` when `K=0`. If `t<=k+2q+2`,
then `M<=4*2^k*4^q<=16*2^k*(K+1)^2`; this is absorbed by (2).
Otherwise choose any subspace `H` of codimension `q`. Its `2^q>K`
affine cosets cannot all meet the exceptional set, so choose one
`x_0+H` that avoids it.

Consider the symmetric bilinear form `B(u,w)=u dot Lw`. Restricting its
matrix to `H` deletes `q` rows and `q` columns in an adapted basis.
Consequently

```text
rank(B|H) >= rank B-2q >= t-k-2q.                            (3)
```

Write `H=rad(B|H) direct-sum U`. Symmetry makes `B|U` nondegenerate:
a vector in `U` orthogonal to `U` is also orthogonal to the radical,
hence belongs to the radical and is zero. Set `n=dim U`; then
`n>=t-k-2q`. The whole affine subspace `x_0+U` avoids the exceptions.

Choose a matrix `T` whose columns are a basis of `U`. Restricted label
values on rows `x=x_0+Ty` are

```text
T^transpose a_x = S y+d,
S=T^transpose L T,       d=T^transpose(Lx_0+v).              (4)
```

The symmetric matrix `S` is invertible. Translating the coordinate `y`
by `S^(-1)d` makes (4) equal to `Sy`. This translation only changes a
constant sign in each Walsh column, which is absorbed by conjugating
that column's Gaussian-prime orientation. All `m=2^n` selected points
are actual original rows in the original containing arc.

## The restricted physical model and its constant columns

Every nonzero restricted label has exactly `2^(t-n)` original labels
above it. Thus it has exactly `b*2^(t-n)` distinct physical prime
columns. The original zero label is absent only in the fiber of the
zero restricted label. A flip at a row outside the selected affine
subspace changes no selected entry. At each nonzero selected coordinate
`y`, the assigned restricted label is `Sy`, nonzero and different for
different rows. Its physical column remains distinct from all other
assigned columns.

The selected coordinate `y=0` has restricted assigned label zero.
This is harmless. Every phase certificate used in the two linked
theorems has `c_0=0` as well as `sum c_y=0`. For a constant restricted
column the unflipped contribution to `sum c_y s_yj` is zero, and a
possible flip at `y=0` also contributes zero. Hence these columns
cancel exactly from every signed product and every composite height.

If `W_const` denotes their total prime-log weight and `D` the original
common Gaussian log norm, the certificates have the same height
inequalities as the linked theorems with effective content

```text
D_eff=D+W_const>=0,
W=W_nonconstant+D_eff.                                      (5)
```

Equation (5) is accounting in the exact signed products. It does not
assert that the constant-column product divides the selected zero row
in the same orientation as every other selected row. This distinction
causes no loss: all certificates avoid that row, while their angular
upper bound still uses the original `W`.

## Classification of the restricted symmetric form

A nondegenerate symmetric binary form is either alternating or has an
orthonormal basis. Here is the elementary classification needed here.
If a vector has self-pairing one, split off its nondegenerate line and
continue in its orthogonal complement. If every self-pairing is zero,
nondegeneracy supplies a pair with mutual pairing one; split off that
alternating plane and continue. This decomposes the form into lines
with matrix `[1]` and alternating planes.

If at least one line occurs, it can absorb any alternating plane.
For a line vector `z` and a plane basis `u,w` with `B(u,w)=1`, the
three vectors `z+u`, `z+w`, `z+u+w` are an orthonormal basis of their
three-space. Repeating yields the identity matrix. If no line occurs,
the form is nondegenerate alternating. An invertible change of selected
row coordinates and its dual relabel the Walsh columns, preserving
the physical primes and all heights. Thus (4) falls under one of the
two linked arbitrary-weight phase theorems.

Those theorems, including their certificates' avoidance of row zero,
give

```text
m <= B_C [1+(W/log(2+W))^(2/3)]                              (6)
```

for a constant depending only on `C`. Bounded small dimensions are
absorbed in `B_C`. Combining (3) and (6),

```text
M <= 2^(k+2q) m <= 4*2^k*(K+1)^2 m,
```

which proves (2), also covering the small-dimensional case separated
above. No actual point, common factor or constant prime weight has
been assigned a smaller radius in this deduction.

For example, with fixed `k` and `K<=M^delta`, `0<=delta<1/2`, (2)
gives

```text
M=O_(C,k,delta)((log R/loglog R)^(2/[3(1-2delta)]))
```

as `R` tends to infinity. This exponent is below one for `delta<1/6`.
The approximation by a symmetric affine assignment remains a
substantive hypothesis in all these statements. The general uniform
endpoint bound is still unproved.

The [checker](check_walsh_symmetric_affine_restriction_growth.py) tests
the rank loss, exception-free coset selection, nondegenerate restriction,
affine translation and literal physical-column dictionary over `F_2`.
