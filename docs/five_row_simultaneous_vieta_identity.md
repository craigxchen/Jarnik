# A projective identity for the five simultaneous other roots

Let `Q` be a simultaneous `SL_2` invariant, homogeneous of degree two
in each of five projective rows. Let the old rows `v_1,...,v_5` be
distinct and satisfy `Q(v)=0`. Assume every row quadratic is unramified
at its old root, and let `v_i'` be its other root, with the other four
rows held at their **old** values. Then

```text
sum_(i<j) det(v_i,v_j) det(v_i',v_j')
          / (det(v_i,v_i') det(v_j,v_j')) = 6.           (1)
```

All five replacements in (1) come from the same section `Q` and the
same source tuple. They are not five successive replacements. The
formula allows an infinite old or new root and is unchanged by
independent rescaling of any of the ten row representatives.

This gives a joint geometric constraint that the individual
[Vieta content bounds](five_row_vieta_swap_content_and_residue_audit.md)
do not supply. It does not yet give a saving in the radius or in the
primitive-residue growth exponent.

## A polynomial first- and second-derivative identity

Choose an affine chart in which every old row is `v_i=(1,x_i)` and
write `q(x)=Q((1,x_1),...,(1,x_5))`. Set

```text
f(t)=product_j(t-x_j),   d_i=f'(x_i),   S=sum_j x_j,
lambda_i=partial_i q,   a_i=(partial_i^2 q)/2.
```

The following identity holds as a polynomial identity, before imposing
`q=0`:

```text
6 d_i lambda_i
 - sum_j d_j(2x_i+3x_j-S)a_j = K_i q,                 (2)

K_i=8x_i^3-6s_1 x_i^2+4s_2 x_i
    -2s_1^3+8s_1 s_2-18s_3,
```

where `s_r` is the elementary symmetric polynomial of degree `r` in
the four coordinates other than `x_i`.

Here is a finite exact proof of (2). The degree-two-each invariant
space has dimension six and is spanned by the six bracket graphs
constructed in the
[integral graph-basis certificate](five_row_gradient_discriminant_core_count.md).
Both sides of (2) are linear in `q`. The accompanying
[symbolic checker](check_five_row_simultaneous_vieta_identity.py)
expands those six literal graph polynomials and verifies equality of
every rational coefficient in (2), for each of the five indices.
Thus its thirty polynomial identities prove (2) for every section in
the invariant space; these are not tests at selected numerical points.

At `q=0`, (2) becomes the rank-two formula

```text
d_i lambda_i=(x_i/3) U + V/2-SU/6,
U=sum_j d_j a_j,          V=sum_j x_j d_j a_j.          (3)
```

## Derivation of the five-swap formula

The three infinitesimal invariance identities are

```text
sum_i partial_i q=0,
sum_i x_i partial_i q=5q,
sum_i x_i^2 partial_i q=2S q.
```

At a zero of `q`, they imply

```text
lambda_i=(A x_i+B)/d_i                                (4)
```

for two scalars `A,B`. Indeed, the kernel of the three-row Vandermonde
matrix has dimension two. Its basis is given by the vectors
`(1/d_i)_i` and `(x_i/d_i)_i`, using the partial-fraction identities
`sum_i x_i^r/d_i=0` for `r=0,1,2,3`. Unramifiedness gives every
`lambda_i!=0`, so `(A,B)!=(0,0)`.

Put `u_i=a_i/lambda_i` and `U_r=sum_i x_i^r u_i`. Comparing (3) and
(4) gives `U=3A` and `V=2B+SA`, hence

```text
A U_1+B U_0=3A,
A U_2+B U_1=SA+2B.
```

The determinant of these two homogeneous equations in `(A,B)` is zero:

```text
(U_1-3)(U_1-2)=U_0(U_2-S).                            (5)
```

An other-root representative, valid also when `a_i=0`, is

```text
v_i'=(a_i, a_i x_i-lambda_i).
```

It has `det(v_i,v_i')=-lambda_i!=0`. Its `ij` summand in (1) is

```text
(x_i-x_j)^2 u_i u_j + (x_j-x_i)(u_j-u_i).
```

Summing these expressions gives

```text
U_0 U_2-U_1^2+5U_1-SU_0=6,
```

by (5). This proves (1) without division by any `a_i`. A common
projectivity can put any five distinct old real or rational points in
a finite chart. Every summand is invariant under that projectivity
and under independent row rescaling, which proves the full projective
statement.

## The identity is sufficient on the unramified locus

For the same six graph basis sections `q_0,...,q_5`, let `J(x)` be the
matrix whose first row is `(q_r(x))_r` and whose next five rows are
`((partial_i^2 q_r)(x)/2)_r`. There is an exact determinant identity

```text
det J(x)=6 product_(i<j)(x_j-x_i)^2.                    (6)
```

The checker proves (6) by polynomial determinant expansion, including
its sign for the specified basis order. Consequently prescribing the
value of `q` and its five diagonal second derivatives determines a
unique section whenever the old coordinates are distinct.

This gives a converse to (1). Start with five proposed new directions
distinct from their respective old directions. Write each uniquely as
`(u_i,u_i x_i-1)`; an infinite new root has `u_i=0`. Suppose (1)
holds. Equivalently, the moment matrix

```text
M=[ U_1-3     U_0   ]
  [ U_2-S     U_1-2 ]
```

has determinant zero. It has rank exactly one: a zero matrix would
require both `U_1=3` and `U_1=2`. Let `(A,B)` span its kernel. Impose
the additional open conditions `A x_i+B!=0` for all five old points.
Set

```text
lambda_i=(A x_i+B)/d_i,    a_i=u_i lambda_i.
```

By (6), there is a unique section with `q(x)=0` and these five `a_i`.
The moment equations and (3) make its actual first derivatives equal
to the prescribed `lambda_i`. Therefore it is unramified in every
row, and its five other roots are exactly the proposed directions.
Changing the kernel generator scales the section, so it is unique
up to scalar. Over rational source and new directions this construction
gives a rational section, which can be scaled to integer coefficients.

Thus on this open locus, (1) is the complete algebraic compatibility
condition for five other roots of one invariant section. An additional
independent projective equation cannot supply the missing arithmetic
estimate there. The inverse construction gives no small bound for the
integer coefficients or their contents.

## Scope and checks

The checker also constructs exact rational numerical zeros, computes
their five row quadratics, verifies their other roots, and checks (1)
under separate row rescalings and a common projectivity. It includes
cases with an infinite new root. Reducibility of `Q` is permitted in
this identity; the required assumptions are distinct old directions
and nonzero row derivatives.

For example, at `x=(-4,-3,-2,-1,1)`, the first two graph-basis values
are `6,8`. For `Q=8G_0-6G_1`,

```text
lambda=(-24,96,-136,72,-8),
a=(24,96,-8,-36,-4),
other roots=(-3,-4,-19,1,-1).
```

All five self brackets are nonzero, and (1) gives exactly six.
This fixture is an algebraic check, not an endpoint configuration
with a small coefficient section.

The identity alone places no upper bound on the Gaussian denominator
needed to put a new phase on the old circle. That denominator remains
the unresolved cost in using a near replacement for descent.
