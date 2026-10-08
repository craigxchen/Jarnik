# Coherent boundary basepoints in an actual numerical relation space

Common points of all first-smaller-cut images do not, by themselves,
exclude boundary points. There is an explicit model inside the actual
degree-two-each invariant system: it consists of numerical zero relations
at distinct rational directions, can have dimension greater than half the
full invariant space, has polynomial gcd one, and has nonzero restriction
images. Its common points at every cut nevertheless arise from one global
degeneration of the labelled directions. At every cut with nonzero image,
all common points are outside the distinct-direction locus.

This is a geometric compatibility model, not a counterexample to the
arithmetic short-relation criterion. In particular, no required coefficient
height is proved for a basis of the large relation space, and the space
need not be the entire collection of relations below a given height. No
full conductor profile or endpoint localization is asserted.

The invariant spaces and restriction maps are those of
[the all-cut hierarchy](all_cut_invariant_relation_hierarchy.md).
The earlier [shared-factor example](smaller_cut_relation_rank.md) gives a
two-dimensional polynomial space without imposing numerical vanishing.
The construction below adds actual numerical zeros, polynomial gcd one,
large dimension, and a coherent choice of all boundary points.

## 1. A quartet numerical zero and the large space

Fix an even `m=2q>=16` and distinct rational binary directions
`P_1,...,P_m`, using integer representatives. Let `V_m` be the rational
SL2 invariant space of degree two in each row. Let `K_k(m)` vanish on
every collapsed restriction with inside cardinality at least `q-k+1`.
Let `Lambda` be the kernel of numerical evaluation at the fixed directions.

For a quartet `A={1,2,3,4}`, set

```text
u=[12][34](P),       v=[13][24](P),
Q_A = v^2 [12]^2[34]^2-u^2 [13]^2[24]^2.
```

The two matching polynomials are independent; thus `Q_A` is a nonzero
polynomial, of degree two in each quartet row. It evaluates to zero, and
both integers `u,v` are nonzero. Write `B=A^c` and define

```text
U = (Lambda intersect K_2(m)) + Q_A K_1(m-4 on B).       (1)
```

Then `U subset Lambda intersect K_1(m)`. A four-row invariant can have
a nonzero collapsed restriction only when at most two of its rows are
inside. An invariant in `K_1(m-4)` can have a nonzero restriction only
when at most `q-3` of its rows are inside. Their product therefore
vanishes whenever at least `q` rows are inside, proving membership in
`K_1(m)`. Numerical vanishing follows from the factor `Q_A`.

Evaluation is nonzero on `K_2(m)`: a product of four disjoint triangle
invariants and squared matching brackets on the remaining rows lies in
that space and is nonzero at distinct directions. Hence

```text
dim U >= dim K_2(m)-1
      >= r-s-binom(m,q-1)t-1,
r=[x^m](1+x+x^2)^m-[x^(m-1)](1+x+x^2)^m,
s=binom(m,q)/2,       n=q+1,       t=n(n-3)/2.          (2)
```

For example the lower bound is `670484982` at `m=24`, while
`r=834086421`; it is `424540705110` at `m=30`, while
`r=439742222071`. The exact inequality making the lower bound greater
than `r/2` also holds at `m=400`. These are dimensions of the specified
space, not counts of relations below an arithmetic height threshold.

## 2. Every nonzero image has an explicit boundary point

Fix a first-smaller cut `S`, with `|S|=q-1`, and let
`a=|S intersect A|`. The term `Lambda intersect K_2(m)` has zero
restriction at every such cut. For the second term in (1):

* If `a<=1`, then at least `q-2` rows of `B` are inside, so the
  `K_1(m-4)` factor restricts to zero.
* If `a>=3`, the four-row invariant `Q_A` restricts to zero.
* If `a=2`, `Q_A` restricts to a scalar `c_I`, where `I=S intersect A`,
  and the other factor has its first-smaller restriction on the `n-2`
  outside rows from `B`.

The integral surjectivity proved in
[the full-restriction theorem](full_invariant_restriction_image_and_index.md)
therefore gives exactly

```text
T_S(U) = c_I W_(n-2,2),                               (3)
```

embedded as quadrics depending only on those `B` variables, when `a=2`.
Here `W_(r,2)` is the squarefree translation-invariant quadratic space.
All other cut images are zero. For `I={1,2}`, the scalar is `-u^2`,
so (3) is nonzero and has dimension `(n-2)(n-5)/2`.

Let `X_n` be the projective image closure of the matching quadrics on
`P(Q^n/Q1)`. At a cut with `a=2`, write `i,j` for the two outside rows
from `A`. Set all outside `B` variables to zero, and give `z_i,z_j`
nonzero values. This is not a universal basepoint of the quadratic map:
using two different outside `B` rows `b,c`, the matching coordinate
`(z_i-z_b)(z_j-z_c)` is nonzero. Every quadric in (3), however, vanishes.
Thus this point of `X_n` annihilates the whole restriction space.

Its projective image is independent of the two nonzero values. Indeed a
squarefree quadric evaluated at `alpha e_i+beta e_j` is `alpha beta`
times its `z_i z_j` coefficient. Denote this rational image point by
`p_ij`; geometrically it is the image of the line between two frame
basepoints after resolution.

If (3) is nonzero, it has no common point coming from `n` distinct outside
directions. Four distinct outside `B` coordinates give a nonzero matching
quadric in `W_(n-2,2)`. Thus all its common points lie outside that open
image locus. The boundary issue is intrinsic to this model, rather than
an arbitrary bad choice of one available common point.

## 3. The choices are compatible with one labelled degeneration

Choose pairwise distinct nonzero rational numbers `a_i` for `i in A`
and pairwise distinct rational numbers `b_j` for `j in B`. Consider

```text
z_i(t)=a_i       if i in A,
z_j(t)=t b_j     if j in B.                            (4)
```

For generic nonzero `t` all `m` directions are distinct. For every cut,
apply its matching-quadrics map to the outside coordinates from (4),
and take the projective limit as `t` tends to zero. Such a limit exists:
divide the coordinate polynomial vector by its smallest nonzero power
of `t` and then specialize. It lies in `X_n`.

When the restriction image is nonzero, `a=2` and the limit is exactly
`p_ij` from Section 2. At every other cut the restriction image is zero,
so its limit is also a common point. In particular, when only one `A`
row is outside the raw quadratic map at `t=0` is a frame basepoint;
the leading nonzero coefficient supplies its well-defined projective
limit. When no `A` row is outside, all quadratic coordinates scale by
`t^2` and the limit is the image of the distinct `b_j`.

Consequently the tuple of common points lies in the closure of the
simultaneous map from **one** labelled direction configuration to all
the cut images. An overlap condition satisfied by every such limit
cannot exclude this example. This does not assert that arbitrary
independently chosen local common points have that compatibility.

## 4. The global polynomial gcd is one

For every even `r>=12`, the polynomial gcd of the entire space `K_2(r)`
is one. It contains all permutations of the graph monomial consisting
of four disjoint triangles and doubled matching edges on the remaining
rows. All factors of these monomials are irreducible brackets. Any fixed
edge can be avoided by one of the permutations, so no bracket divides
every monomial. There is therefore no common nonconstant factor.

Choose a second quartet `A'` disjoint from `A`, and form its nonzero
numerical zero `Q_(A')` by the same formula as in Section 1. Inside
`Lambda intersect K_2(m)` are the two spaces

```text
Q_A K_2(m-4 on A^c),
Q_(A') K_2(m-4 on (A')^c).                            (5)
```

Their membership follows from the same inside-count argument: the maximum
is now `2+(q-4)=q-2`. Each polynomial gcd in (5) is precisely its quartet
factor. Those factors use disjoint variable sets and are coprime.
Thus `Lambda intersect K_2(m)`, and hence `U`, has gcd one in the full
polynomial ring `Q[X_1,Y_1,...,X_m,Y_m]`.

This rules out a global common-factor explanation for the whole space
`U`. It does not remove the visible quartet tensor in its quotient by
`K_2(m)`. Controlling the height of that quotient can be additional
arithmetic information and is not supplied by the gcd calculation.

## 5. Exact scope of the model

The starting directions can be any distinct rational directions. Standard
Gaussian denominator clearing realizes them as distinct points on an
actual integer circle, and common Gaussian division makes the tuple
primitive. Thus the numerical relation system here is the same actual
invariant system used by the arithmetic argument.

What is not asserted is that a basis of `U` has the coefficient height
required by the conductor-sensitive Chow criterion. Nor is `U` asserted
to equal the span of every numerical zero below that height. Additional
short relations can enlarge a restriction image and eliminate its common
points. The model supplies no full fair source profile or endpoint arc.

The proved limitation is specific: large dimension, actual numerical
vanishing, polynomial gcd one, and even compatibility of all common
points with one global degeneration do not by themselves force an
interior common point or a contradiction. A continuation must retain
quantitative arithmetic information or impose a condition excluding
these coherent collision limits.

The subsequent [quartet height descent](quartet_boundary_height_descent.md)
supplies such quantitative information for this particular construction.
In a full-profile system, any integer element of the larger space
`K_2(m)+Lambda_A K_1(B)` outside `K_2(m)` has coefficient logarithm at
least `2^(m-3)(1-eta)w-16 sigma`. Thus the nonzero extension in (1)
cannot contribute a restriction at the Chow-form short height. This
does not remove short vectors inside `K_2` or classify arbitrary
numerical relations outside that quartet-supported space.

The [exact checker](check_coherent_boundary_basepoints_of_numeric_relations.py)
checks all `11440` first-smaller cuts at `m=16`, the common-point limits,
the nonzero image rank in (3), the graph gcd certificates, and the exact
dimension inequalities in (2). The general claims follow from the
arguments above; no huge invariant basis is numerically approximated.

## 6. Even coherent interior points do not force one quartet space

There is a separate linear dimension existence theorem already at `m=10`.
It does not use the deeper kernel: the hierarchy terminates at `K_2(10)=0`.
The exact dimensions and number of first-smaller cuts are

```text
dim V_10=603,       dim K_1(10)=603-binom(10,5)/2=477,
number of cuts |S|=4: binom(10,4)=210,
n=6,               dim W_(6,2)=9.                    (6)
```

Fix any distinct rational directions `P_1,...,P_10` for numerical
evaluation, and independently fix one globally labelled set of distinct
rational affine coordinates `b_1,...,b_10`. At every cut let `x_S in X_6`
be the matching-quadrics image of its six outside `b` coordinates.
These are interior image points, and they are coherent by construction:
all come from the same single configuration `b`, without a degeneration.

Define the rational vector space

```text
L_b={W in K_1(10): T_S(W)(x_S)=0 for every |S|=4},
U_(P,b)=L_b intersect Lambda_P.                       (7)
```

Each cut imposes at most one rational linear equation. Numerical
evaluation at `P` imposes at most one more. Therefore

```text
dim L_b >= 477-210=267,
dim U_(P,b) >= 266.                                  (8)
```

Every member of `U_(P,b)` is an actual degree-two-each numerical zero.
At every cut its entire restriction space has the specified common
interior point `x_S`. Rational bases can be made integral by clearing
denominators, but no bound for those denominators or coefficient heights
is asserted.

For any fixed quartet `A`, the space from the height-descent theorem
simplifies here to

```text
E_A=Lambda_A K_1(A^c),
dim E_A=2*5=10.                                      (9)
```

Indeed the quartet invariant space has dimension three and its nonzero
numerical evaluation has kernel dimension two; the complementary six-row
balanced kernel has dimension `15-10=5`. Multiplication on disjoint
variable sets is injective on their tensor product. Since `K_2(10)=0`,
there is no quotient ambiguity in (9).

It follows from (8)--(9) that `U_(P,b)` is not contained in any one `E_A`.
It is not contained even in their finite union: a rational vector space
over the infinite field `Q` cannot be a union of finitely many proper
linear subspaces. Thus actual numerical zero relations with all the
coherent interior common-point conditions can lie outside every fixed
quartet-supported space, even when the deeper kernel is zero.

This is an existence theorem from linear dimensions, not an explicit
basis construction or a short-height/full-profile example. It does not
refute a classification that additionally uses the arithmetic coefficient
threshold. It shows precisely that excluding boundary points and imposing
coherence of the common points are insufficient, by themselves, to deduce
membership in a fixed `E_A`; quantitative arithmetic remains necessary.

The checker verifies the dimension data (6)--(9) and all 210 nonzero
interior matching-coordinate vectors for `b_i=i+1`. It does not claim to
construct a basis meeting the height required by the endpoint argument.
