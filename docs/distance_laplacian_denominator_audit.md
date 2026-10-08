# Squared-distance Laplacians: an exact full-denominator obstruction

The positive spanning-tree sum of all squared unit-circle chords is a
natural all-pairs integer certificate. Its arithmetic content can be
computed exactly on a large class of actual rational circle configurations:
the reduced denominator is the largest permitted denominator, `N^(k-1)`.
There is no conductor content to cancel, even on a four-point family with
bounded endpoint normalization. This rules out one specific proposed
Laplacian denominator estimate; it does not prove or disprove the uniform
short-arc conjecture.

The distance matrix and its centered positive semidefinite matrix also
have exact low-rank formulas. Their nonzero determinants reduce to
triangle areas. The positive sum of their order-two centered minors gives
the usual triangle exponent unless further arithmetic content is supplied.

## 1. A positive all-pairs invariant and the precise required estimate

Let `x_1,...,x_k` be distinct rational unit-circle points, put

```text
d_ij=|x_i-x_j|^2,       delta=max_(i,j) |x_i-x_j|,
N=lcm_(i<j) denominator(d_ij),
tau=sum_(spanning trees T on [k]) product_(ij in T) d_ij.
```

Here `N` is the least primitive squared radius, as in
[homogeneous_central_content_radius.md](homogeneous_central_content_radius.md).
Equivalently, `tau` is any principal cofactor of the weighted Laplacian
with off-diagonal entries `-d_ij` and diagonal entries `sum_(j!=i)d_ij`.
The identity follows by expanding the determinant using an oriented
incidence matrix: a set of `k-1` edges contributes squared incidence
determinant one exactly when it is a tree, and zero otherwise.

There are `k^(k-2)` labelled trees, so

```text
0 < tau <= k^(k-2) delta^(2(k-1)).                    (1)
```

Write the reduced fraction as `tau=P/D`, with `P,D` positive coprime
integers. Then `P>=1` gives

```text
delta >= k^(-(k-2)/(2(k-1))) D^(-1/(2(k-1))).        (2)
```

Consequently a bound `D<=B_k N^((k-1)/2)` would imply the desired
`N^(-1/4)` exponent for this fixed number of points. The automatic
denominator estimate is only `D | N^(k-1)`. The next section shows that
it is attained exactly; positivity does not improve it.

This distinction is separate from real cancellation. For fixed distinct
real numbers `t_i`, taking circle angles `epsilon*t_i` gives

```text
tau = epsilon^(2(k-1))
      * sum_T product_(ij in T)(t_i-t_j)^2
      + O_k,t(epsilon^(2k)).                         (3)
```

The leading coefficient is strictly positive. The squared-distance
Laplacian therefore has exactly its ordinary homogeneous vanishing order
on a shrinking arc; the rank-three distance matrix does not give it an
additional curvature factor.

## 2. Full reduced denominator on actual circle configurations

Choose distinct even integers `A_i` such that the odd positive integers

```text
q_i=A_i^2+1
```

are pairwise coprime. Use the rational phases

```text
z_i=(A_i+i)/(A_i-i).
```

Their squared chords are exactly

```text
d_ij=4(A_i-A_j)^2/(q_i q_j).                         (4)
```

If a prime divides `q_i` and `A_i-A_j`, it also divides `q_j`.
Thus (4) is reduced at every prime of its denominator, and

```text
N=product_i q_i.                                    (5)
```

For a tree `T`, its denominator divides
`product_i q_i^(deg_T(i))`; factors from other edges can cancel part of
this denominator. Fix `p|q_i` and write `e=v_p(q_i)>0`. Every edge
incident to `i` has valuation exactly `-e`, and all other edges have
nonnegative valuation. The unique tree of degree `k-1` at vertex `i`
is the star centered at `i`. That summand has valuation `-(k-1)e`, while
all other summands have strictly larger valuation. There can be no
cancellation at this unique minimum. Hence

```text
v_p(tau)=-(k-1)e,       D=N^(k-1).                  (6)
```

This is an all-prime content audit, since every denominator prime divides
one of the pairwise coprime `q_i`. Equivalently the explicit integer

```text
P=4^(k-1) sum_T [product_(ij in T)(A_i-A_j)^2
                 * product_i q_i^(k-1-deg_T(i))]
  =N^(k-1) tau
```

satisfies `gcd(P,N)=1`. Reduction modulo a prime of `q_i` leaves only
the star term, which is a unit. In particular the frequently proposed
division of this positive integer by a positive power of `N` is false.

These are actual Gaussian integer-circle configurations: set
`h_i=A_i+i` and `Z_i=h_i product_(j!=i) conjugate(h_j)`. Then
`|Z_i|^2=N` and their relative phases are the phases above. Formula (5)
shows this realization already has the least squared radius.

For every fixed `k`, such examples can lie in arbitrarily small angular
arcs. For example set `A_1=2T` and recursively

```text
A_i=2T product_(j<i)q_j.
```

Every new `q_i` is one modulo all previous `q_j`, so they are pairwise
coprime; all phases tend to one as `T` tends to infinity. The recursive
heights grow rapidly. This example is not claimed to have bounded
`angular_span*N^(1/4)`.

## 3. Four points at a bounded endpoint scale

There is also a precise fixed-endpoint example for `k=4`. Take

```text
A_i=2(T+i),       i=0,1,2,3,       T=0 mod 5, T>0.
```

If an odd prime `p` divides two norms at index difference `d`, then,
comparing the two roots of `-1` modulo `p`, either `p|d` or
`p|(d^2+1)`. For `d=1,2,3`, the only possible prime congruent to one
modulo four is five. On the chosen progression only `i=1` has norm
divisible by five. The four norms are therefore pairwise coprime.

Writing `Theta` for their angular span gives

```text
Theta=2 arctan(1/(2T))-2 arctan(1/(2(T+3))),
N=product_(i=0)^3 (4(T+i)^2+1),
Theta*N^(1/4) -> 12,       denominator(tau)=N^3.      (7)
```

Thus, for any fixed `C>12`, this family eventually consists of four
integer points on an arc of length at most `C sqrt(R)`, with `R=sqrt(N)`,
and still violates every estimate `D<=B N^(3/2)` with fixed `B`.
The bounded endpoint constant agrees with the four-point family in
[ordered_consecutive_ptolemy_nonlinear_audit.md](ordered_consecutive_ptolemy_nonlinear_audit.md),
but (6) is a different arithmetic calculation.

There is no contradiction to uniform boundedness of the number of points:
this family always has exactly four points. It also says nothing against
an endpoint-restricted content estimate for sufficiently large `k` or a
smaller specified endpoint constant. Such an estimate would need to use
those restrictions, rather than Laplacian positivity alone.

## 4. What distance rank and conditional negativity actually give

Let `X` have rows `x_i`, `1` be the all-ones column, and
`H=I-11^T/k`. The squared-distance matrix and centered Gram matrix are

```text
D_dist=2*11^T-2XX^T,       B=-H D_dist H/2=HXX^TH.   (8)
```

Thus `D_dist` has rank three for at least three distinct circle points,
`B` is positive semidefinite of rank two, and for `sum_i c_i=0`,

```text
c^T D_dist c=-2|sum_i c_i x_i|^2 <= 0.              (9)
```

Let `Delta_ijk=det(x_j-x_i,x_k-x_i)`. If `Y` has rows
`(1,x_i^(1),x_i^(2))`, then
`D_dist=Y diag(2,-2,-2)Y^T`. Every order-three minor is therefore
eight times a product of two signed triangle determinants. In particular

```text
det(D_dist[{i,j,k},{i,j,k}])
 =8 Delta_ijk^2=2 d_ij d_ik d_jk.                   (10)
```

All larger distance minors vanish. The circle identity used here is
`Delta_ijk^2=d_ij d_ik d_jk/4`, obtained from the triangle
circumradius formula, or directly from sine half-angle differences.

The sum `e_2(B)` of all order-two principal minors has the exact positive
formula

```text
e_2(B)=1/k sum_(i<j<l) Delta_ijl^2
      =1/(4k) sum_(i<j<l) d_ij d_il d_jl.           (11)
```

Indeed the two nonzero eigenvalues of `B` are those of the centered
two-by-two covariance matrix `X^T H X`; taking its determinant and
applying Cauchy--Binet to `Y` proves (11).

In a least integral realization of squared radius `N`, write `a_ijl`
for its integer signed triangle determinant. Scaling gives

```text
k N^2 e_2(B)=sum_(i<j<l) a_ijl^2.                   (12)
```

Every summand is a positive integer because three distinct circle points
are not collinear. Combining (11)--(12) and `d_ij<=delta^2` yields only

```text
delta^6 >= 4/N^2,       delta >= 2^(1/3) N^(-1/3).  (13)
```

The point-count factors cancel exactly. Primitivizing the positive sum
in (12) could change this bound only with an additional, justified content
estimate; its rank and positivity do not supply one. These calculations
cover the displayed distance minors and centered-minor sum, not every
possible invariant built from distance data.

The later [primitive distance Smith theorem](primitive_distance_smith_form.md)
computes all Smith factors of the complete integer distance matrix,
including arbitrary row units and the exact parity-dependent content.
Its endpoint consequence still has the triangle exponent above.

## Verification

[check_distance_laplacian_denominator.py](check_distance_laplacian_denominator.py)
checks the tree sum against its Laplacian determinant, the exact full
reduced denominator on pairwise coprime norm examples, the four-point
endpoint progression, (10)--(12), and the conditionally negative identity.
All arithmetic identities and denominator assertions use exact fractions.
