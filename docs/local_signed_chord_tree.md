# The local signed chord tree

This note records the exact good-prime possibilities for the edge defect

```text
lambda_ij = v_p(q_i)+v_p(q_j)-2v_p(r)-2v_p(Delta_ij).
```

At an odd prime with `p` not dividing `r`, this is the negative valuation of
the squared normalized chord.  Put `a_i=v_p(q_i)` and
`d_ij=v_p(Delta_ij)`, so `lambda_ij=a_i+a_j-2d_ij`.

## First residue level

Two primitive vectors have `d_ij>0` exactly when their projective reductions
in `P^1(F_p)` coincide.  If they are in different residue classes then
`d_ij=0` and

```text
lambda_ij=a_i+a_j.                                  (1)
```

For `Q=I`, an isotropic residue class is a class on which `x^2+y^2=0`.
When `p = 3 mod 4` (for example `p=3,7`), `-1` is nonsquare, so there are no such classes and every
primitive vector has `a_i=0`.

When `p = 1 mod 4` (for example `p=5,13`), `-1` is a square and there are exactly two isotropic
classes.  For `p=5` they are `[1:2]` and `[1:-2]`; for `p=13` they are
`[1:5]` and `[1:-5]`.  A nonisotropic class has `a_i=0`, whereas an
isotropic node has `a_i>=1`.

Thus a determinant-unit edge is positive whenever at least one endpoint is
isotropic.  It is positive for mixed isotropic/nonisotropic pairs as well as
for pairs on the two different isotropic branches.  It is zero between two
distinct nonisotropic classes.

## All depths on a split branch

Choose `iota in Z_p` with `iota^2=-1`.  On the branch
`x-iota*y=0 mod p`, put

```text
z(v)=(x-iota*y)/(x+iota*y).
```

The denominator is a unit, and direct factorization gives, for nodes on this
branch,

```text
a_i=v_p(z_i),             d_ij=v_p(z_i-z_j).          (2)
```

The other branch is identical with the two factors interchanged.  Therefore

```text
lambda_ij=v_p(z_i)+v_p(z_j)-2v_p(z_i-z_j).            (3)
```

If `a_i != a_j`, the ultrametric inequality gives
`d_ij=min(a_i,a_j)` and hence

```text
lambda_ij=|a_i-a_j| > 0.                               (4)
```

If `a_i=a_j=a`, write `z_i=p^a u_i` with `u_i` a unit.  Then

```text
lambda_ij=-2v_p(u_i-u_j).                              (5)
```

This is zero for distinct unit children modulo `p`, and is `-2h` when the
two nodes remain together for exactly `h` further levels.  Negative edges
inside an isotropic branch therefore occur precisely in equal-depth child
clusters.  Negative edges inside a nonisotropic residue class are simply
`-2d_ij`.  The residue classes and the deeper child balls are disjoint as
containers, but the edge sets inside one ball need not be matchings.

For inert primes the formula is simpler.  Since `a_i=0` for every primitive
node,

```text
lambda_ij = -2d_ij       if the residue classes coincide,
lambda_ij = 0           otherwise.                      (6)
```

There are no positive edges at an inert prime.

## General positive Gram form

Let `Q` be primitive positive integral with `det(Q)=r^2`, and assume
`p` does not divide `2r`.  The binary form is regular modulo `p`.  If it is
inert, the anisotropic argument above applies: a primitive vector cannot
have `q(v)=0 mod p`, so `a_i=0`.

If it is split, choose a primitive local basis in which the two isotropic
linear factors are defined over `Z_p`.  Their determinant is a unit.  A
projective coordinate on either branch consequently has unit Jacobian, so
the valuations in (2)--(5) are unchanged.  This is just factorization of a
regular binary form over `Z_p`; the required roots can also be obtained by
Hensel lifting a simple root modulo `p`.  No assertion is made when `p|r`.

## Radius versus the actual cotangent denominator

The intrinsic squared-radius lcm only sees

```text
v_p(N)=max_ij max(0,lambda_ij)
```

at odd `p`.  A negative edge is invisible to this maximum because its chord
is `p`-adically small.  It is visible in an all-edge rational cotangent
normalization: the pair cotangent is, up to a unit,
`B_ij/(r Delta_ij)`.  For a negative edge at a good odd prime, the Gram
identity implies

```text
v_p(B_ij/(r Delta_ij)) = lambda_ij/2 < 0,
```

so clearing every actual pair cotangent forces the common cotangent
denominator `L` to contain `p^{-lambda_ij/2}`.  Thus a family with fixed
positive maximum but arbitrarily negative edges has unbounded `v_p(L)`.  This
rules out fixed `L`; it does not by itself rule out `log L=o(w)` when the
negative depth is itself `o(w)`.  More precisely, every negative edge obeys
`(-lambda_ij/2) log(p) <= v_p(L) log(p)`; summing over the fifteen edges gives
an upper bound `15 log L` for the total negative depth.  Thus an allowed
subpower-`L` error remains `o(w)`.

For example, at `p=5`, for every `M>=1` take

```text
v=( (1,0), (1,5^M), (1,2), (1,-2), (2,1), (1,1) ).   (7)
```

The only nonunit determinant valuations are `d_12=M` and `d_45=1`, while
the norm valuations are `(0,0,1,1,1,0)`.  Hence
`lambda_12=-2M`, the positive maximum is `2`, and the central valuations are

```text
v_5(u_i)=(0,0,M+2,M+1,M+1,M).
```

The all-edge cotangent normalization needs `v_5(L)>=M` because of edge 12.
Moreover, its actual all-edge radius lcm satisfies

```text
log N = log(1+5^(2M)) + O(1) = 2M log(5) + O(1).       (8)
```

Indeed edge 12 has norm `(1+5^(2M))/2`; every edge involving node 2 has
norm dividing `(1+5^(2M))` times a fixed norm, and all edges avoiding node 2
are fixed.  Thus `log L>=M log(5)` is a linear charge in this family's own
radius logarithm, rather than a subpower-`L` example.  Replacing `5` by `13`
and `(1,2),(1,-2),(2,1)` by `(1,5),(1,-5),(5,1)` gives the same table at
`13`; in that case the same argument gives
`N=13^2(1+13^(2M))/2`.

Two small same-branch checks at `5` are

```text
(1,2),(1,7):    (a_i,a_j,d,lambda)=(1,2,1, 1),
(1,7),(1,132):  (a_i,a_j,d,lambda)=(2,2,3,-2).
```

The first is positive despite being on one isotropic branch; the second is a
negative equal-depth child collision.
