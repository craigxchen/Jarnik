# All admissible maps from one large endpoint tuple have the same normal

Fix a Gaussian-primitive source tuple `Z` of `m` distinct points on squared
radius `N`, in an arc of length at most `2N^(1/4)`. Call a rational
projective map admissible if it is nonconformal and its full image has
actual least squared radius `N'<=N` and endpoint constant at most two.
No point is deleted and no condition on intermediate points of an arc is
imposed.

If

```text
m>=512,       N>2^744,                               (1)
```

then all admissible maps have the same **oriented primitive physical
projection normal**. Moreover any two admissible maps differ only by a
conformal transformation on the output. In particular their target
configurations are conformally equivalent, with exactly the same least
radius.

Two structural consequences follow:

* the projective automorphism group of `Z` has order at most two, without
  any binary-allocation hypothesis;
* if a labelled projective endpoint orbit with `m>=512` contains even
  one realization with squared radius exceeding `2^744`, the whole orbit
  has at most two conformal equivalence classes of endpoint realizations.

These theorems constrain maps and representatives already known to exist.
They do not construct a nonconformal admissible map or prove a uniform
endpoint point count.

## 1. The physical normal is independent of actual anchoring

In actual source and target phase coordinates, write a primitive integral
triangular representative as

```text
M=((a,b),(0,d)),       a>0, d!=0,
T=a^2+b^2+d^2, W=a^2-b^2-d^2+2iab, K=Norm(W)>0.
```

At the actual physical source anchor `z_j`, put

```text
C=W z_j,       g=gcd_Z(Re C,Im C),       n=C/g.        (2)
```

This is the normal in the physical stretch projection formula. The
[reciprocal-stretch conductor theorem](mobius_reciprocal_stretch_grid.md)
gives `Q=Nclip|g`.

Here is the exact covariance, including orientation. If `H` represents
the relative phase `z_k/z_j=H/conjugate(H)`, source reanchoring multiplies
the matrix on the right by complex multiplication by `H`. Its Gram
coefficient becomes `W conjugate(H)^2`. Hence

```text
[W conjugate(H)^2] z_k = Norm(H) W z_j.              (3)
```

Output rotations, reflections, and scalar normalizations multiply the
Gram matrix by a positive scalar. Ordinary primitive matrix reduction
does the same. Thus (2) defines one oriented primitive integer vector
for the original map on this fixed physical source tuple, regardless of
the actual source or target anchors used to compute it. Even half-angle
norms in (3) cause no missing radius or sign factor.

For matrices written in a common source phase frame, equality of these
normals means their Gram matrices have the same ordered eigenlines:
the physical normal records the doubled angle of the expanding
half-angle eigenline. The sign distinguishes the expanding and
contracting lines.

## 2. Small height and concentration near the source arc

Use the constants from
[the finite-sample critical triangular reduction](mobius_finite_sample_critical_form.md):

```text
h_m=g_m+b_m,
a,|d|<=D_m N^h_m,
Q>=2^(-m q_m)N^(1-q_m),
N^z_m/512 <= kappa <= 4(2N)^u_m.
```

All source units and the actual least target norm are included in these
estimates. Since `sqrt(K)<=T` and
`T<=5a|d|(2N)^u_m`, equation (2) gives

```text
|n| <= H_m N^nu_m,
H_m=5D_m^2 2^(u_m+m q_m),
nu_m=2h_m+u_m-1/2+q_m
    =23/(2m)+O(m^(-2)).                             (4)
```

At this selected actual anchor, the normalized squared stretch is
`x_j=a/|d|<=D_mN^h_m`. The physical direction of `C` relative to `z_j`
is the direction of `W`. Exactly,

```text
|C/|C|+z_j/sqrt(N)|^2
 =2+2(a^2-b^2-d^2)/sqrt(K)
 <=4a^2/sqrt(K)
 =4x_j/(kappa-kappa^(-1)).                          (5)
```

For `kappa>=2`, the lower condition-number bound gives

```text
|C/|C|+z_j/sqrt(N)|
 <=64 sqrt(D_m) N^mu_m,
mu_m=(h_m-z_m)/2.                                  (6)
```

This is a finite-sample stretch statement at the selected point. It does
not assume that a pole lies outside either sampled arc.

## 3. Integer determinants force one oriented normal

Take any two admissible maps. Their selected source anchors may differ,
but those anchors have normalized physical distance at most `2N^(-1/4)`.
Equations (4)--(6) therefore imply

```text
|det(n_1,n_2)|
 <=H_m^2(128sqrt(D_m)+2) N^(2nu_m+mu_m),
2nu_m+mu_m=-1/4+209/(8m)+O(m^(-2)).                 (7)
```

In combining the two angular terms, we used `mu_m>=-1/4`, which follows
from `h_m>=0` and `z_m<=1/2`.

Here are coarse explicit certificates for (1). The previously proved
uniform estimates give

```text
D_m<=2^14, h_m<=7/m, q_m<=8/m,
u_m<=1/2+2/m, z_m>=1/2-6/m.
```

For `m>=512`, they imply `H_m<2^39`, `128sqrt(D_m)+2<2^15`, and

```text
nu_m<=24/m,
mu_m<=-1/4+13/(2m),
2nu_m+mu_m<=-1/4+109/(2m)<=-1/8.
```

For the height constant one may use `u_m<=2/3` and the exact certificate
`5^3*2^110<2^117`. Also `z_m>=1/4`, so (1) guarantees `kappa>=2`.
Thus (7) gives

```text
|det(n_1,n_2)|<2^93 N^(-1/8)<1.                    (8)
```

The determinant is an integer, so it vanishes. Since the vectors are
primitive they differ by at most a sign. The unit-direction distance in
the same proof is less than `2^15N^(-1/8)<2`, excluding the opposite
sign. Hence `n_1=n_2`.

## 4. There is at most one nonconformal image class

Represent two admissible maps `A,B` in one common source half-angle
frame. The normal theorem says that `A^T A` and `B^T B` have a common
ordered eigenbasis. Their output orthogonal factors, including possible
reflections, do not affect singular values. Consequently

```text
kappa(B A^(-1))
 =max(kappa(B)/kappa(A),kappa(A)/kappa(B))
 <=2048 * 2^u_m N^(u_m-z_m).                        (9)
```

This is a statement about the original rational composite. Diagonalizing
the two Gram matrices for this calculation does not change the rational
coordinates to which an arithmetic theorem is applied.

Both actual target norms satisfy the arbitrary-unit retention bound
`N_A,N_B>=N^(1-4/m)/16`. Order them by size, and apply the nonconformal
lower condition-number theorem to the composite in the radius-nonincreasing
direction. If the composite were nonconformal, it would give

```text
kappa(B A^(-1)) >=2^(-9-4z_m) N^[(1-4/m)z_m].        (10)
```

For `m>=512`, the same coarse estimates imply

```text
(2-4/m)z_m-u_m >=1/4,
20+u_m+4z_m <23.
```

Combining (9)--(10) would require `N^(1/4)<2^23`, contrary to (1).
The composite is therefore conformal. This proves the uniqueness of
the nonconformal output class.

A conformal correspondence between two primitive physical tuples is a
common complex similarity or conjugate similarity. Its multiplier and
its inverse are Gaussian integral by Bezout for the two primitive tuples;
the multiplier is therefore a Gaussian unit. In particular it preserves
the exact least squared radius.

## 5. Self-maps: no binary hypothesis and no Klein four group

For a self-map, both its forward map and inverse are admissible from the
same fixed physical source tuple. They may each select different actual
anchors, but (3) transports their normals to the same source coordinates.
Thus inverse closure introduces no unrelated physical gauge.

There is also an algebraic description of this special case. In a common
orthonormal eigenbasis, let `A^T A=diag(lambda,mu)` with `lambda>mu`.
If `A` and `A^(-1)` have the same oriented normal, then
`A A^T=diag(mu,lambda)`. The identity
`A(A^T A)=(A A^T)A` forces `A` to be off-diagonal in this basis, and
therefore `A^2` is scalar. This covers either determinant sign and is
unchanged by projective rescaling.

Two such involutions with the same axes have diagonal product. Equal
axes alone do **not** make that product conformal; its finite projective
order as a self-permutation does. A real diagonal finite-order projective
matrix has eigenvalue ratio `+1` or `-1`, so the product is conformal.

Finally, the primitive conformal self-group is trivial for this source:
unit rotations cannot preserve the short arc, and a unit-axis reflection
would give at most two points by the
[reflection radial bound](reflection_union_radial_count_obstruction.md).
Thus two distinct nonidentity self-maps cannot occur. Equivalently, the
output-class theorem directly says their quotient is a conformal self-map
and hence the identity. The automorphism group has order at most two.

This does not imply that an involution acts freely for arbitrary source
contents. The stronger free-action assertion in the
[binary group theorem](binary_endpoint_projective_group_classification.md)
retains its binary hypothesis.

## 6. At most two endpoint representatives in a large labelled orbit

Suppose three labelled projectively equivalent endpoint configurations
are pairwise conformally inequivalent and all have `m>=512` points.
Choose the one with largest least squared radius as source. If that
radius exceeded `2^744`, the other two would be images under two
admissible maps, and Section 4 would make them conformally equivalent.
Therefore all three radii must be at most `2^744`.

In particular, if the endpoint orbit has even one realization above that
threshold, it has at most two conformal classes in total: combine that
realization with representatives of any two other classes and apply the
preceding argument to the largest of the three.

This leaves open whether the second class exists. Radius minimization,
symmetry, and reanchoring alone do not provide it. The original uniform
endpoint point-count problem therefore remains unproved.

The [checker](check_critical_projection_normal_rigidity.py) verifies exact
physical normal covariance, the inverse/common-axis algebra, the
finite-order product implication, and the rational exponent certificates.
