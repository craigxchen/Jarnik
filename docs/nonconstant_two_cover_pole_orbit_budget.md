# Pole-orbit budgets for the nonconstant `[2]` cover

Let `C/Q` be the genus-five hyperelliptic curve whose twelve branch points
have Galois group `S_12`, and put

```text
V=J[2]={even branch subsets}/<complement>.
```

Assume the triangle boundary is the union of ten rational closed points,
each with a biquadratic residue field of degree four.  The forty geometric
boundary points lift to `40*1024` geometric poles on the unramified
multiplication-by-two cover.

This note determines a much smaller arithmetic orbit budget after adjoining
one target lift.  It does not identify a boundary fiber with the target
fiber.  A Burnside bound handles every boundary affine cocycle at once.
The result closes the ordinary Runge pole surplus for every target-fiber
orbit except degrees `1`, `12`, and `66`; degree `66` is also closed under
the existing ten-distinct-good-primes hypothesis.

## 1. Boundary and target have different Kummer cocycles

Fix a rational base point `A in C(Q)`, a rational target `P`, and a geometric triangle pole
`Q`.  Choose halves

```text
2T=[P-A],       2U_Q=[Q-A].
```

For an automorphism fixing `Q`, the target and boundary affine cocycles are

```text
kappa_P(sigma)=sigma(T)-T,
kappa_Q(sigma)=sigma(U_Q)-U_Q.                       (1)
```

Put `D_Q=U_Q-T`, so `2D_Q=[Q-P]`.  Then

```text
kappa_Q=kappa_P+(sigma(D_Q)-D_Q).                    (2)
```

After adjoining the chosen target lift `T`, the first term vanishes, but
the second need not:

```text
kappa_(Q,L)(sigma)=sigma(D_Q)-D_Q,
2D_Q=[Q-P].                                         (3)
```

Changing `D_Q` changes (3) by a coboundary.  Thus each boundary fiber
depends on the restricted division class of `[Q-P]`; its affine cocycle and
translation kernel are not determined by the target-fiber type.

There is an exact orbit formula.  Let `Gamma_L` act on the forty geometric
base poles `B`.  For each `Gamma_L`-orbit choose `Q`, and let `H_Q` be the
image of its stabilizer on the affine fiber `V` above `Q`.  Then

```text
r_L=sum_(O in Gamma_L\B) number_of_orbits(H_Q on V). (4)
```

This is the orbit decomposition of a group action on a set fibered over
`B`.  The target-fiber orbit list alone cannot replace (4).

## 2. A cocycle-independent Burnside bound

Let an affine group `Gamma<=V semidirect G` project onto a linear group
`G<=GL(V)`, and let `K=Gamma intersect V`.  Every fiber of
`Gamma->G` has `|K|` elements.  For fixed `g in G`, an affine lift

```text
x |-> g x+t
```

has either no fixed point or exactly `|ker(g-I)|` fixed points.  Burnside's
lemma therefore gives

```text
number_of_orbits(Gamma on V)
 <=(1/|G|) sum_(g in G)|ker(g-I)|
 =number_of_orbits(G on V).                          (5)
```

No hypothesis on the affine cocycle or its translation kernel appears in
(5).

For `S_12`, and also for `A_12`, the linear action on `V` has four orbits,
represented by subset weights `0,2,4,6`.  The `A_12` orbits do not split:
the stabilizer of each subset class contains an odd transposition.

More generally let `H_k` be the stabilizer in `S_12` of a subset class of
weight `k<=6`; for `k=6` it may interchange the two halves.  Direct
two-by-two intersection-count enumeration gives

```text
k                    0   1   2   3   4   5   6
|S_12:H_k|           1  12  66 220 495 792 462
number H_k\V         4   6   9  10  12  12  10.      (6)
```

Intersecting `H_k` with `A_12` does not change the last row.  Indeed, for a
target subset and any representative of a class in `V`, the four joint
membership cells contain twelve labels, so one contains at least three.
A transposition inside that cell is odd and fixes both classes.  Hence no
`H_k`-orbit splits under `H_k intersect A_12`.

The checker `check_s12_even_subset_h1.py` enumerates all 1024 classes and
verifies (6), including the half-interchange at `k=6`.  Formula (5) then
bounds an arbitrary joint boundary affine action by the last row of (6).

## 3. The ten biquadratic base poles remain ten

First consider the linear and odd affine target types.  Their affine image
is a graph isomorphic to `S_12`, so its Galois field is the branch splitting
field `K`.  A target lift in the orbit of weight `k` has residue field

```text
L=K^(H_k).
```

Every `H_k` contains an odd transposition.  The only quadratic subfield of
an `S_12` extension is the sign-discriminant field, and it lies in `L` only
if `H_k<=A_12`.  Thus `L` has no quadratic subfield.  A biquadratic field
`E/Q` has a quadratic subfield in every nontrivial subextension, so

```text
E intersect L=Q.
```

Because `E/Q` is Galois, `E` and `L` are linearly disjoint.  Each of the ten
degree-four base poles therefore remains one closed point over `L`.

Its joint stabilizer has linear image `H_k` or
`H_k intersect A_12`.  To see this, `E intersect K` is a normal subfield of
the `S_12` extension of degree dividing four, hence is either `Q` or the
sign-discriminant field.  The two possibilities give precisely those two
images.  By (5)--(6), each base pole contributes at most the corresponding
entry in the last row of (6).

In the full-translation target type, the affine image is
`V semidirect S_12` and a lift has degree 1024.  Its point stabilizer is a
complement `S_12`.  This stabilizer is maximal: any intermediate subgroup
meets `V` in an `S_12`-invariant subspace, which is zero or `V` because `V`
is simple.  Hence the lift field has no proper intermediate field and is
again disjoint from every biquadratic `E`.  Over `LE`, the linear image is a
normal subgroup of `S_12` of index dividing four, hence `S_12` or `A_12`.
Formula (5) bounds each of the ten pole fibers by four orbits.

## 4. Exact pole and archimedean comparison

Combining the ten preserved base orbits with (6) gives:

| Target fiber type | Lift degree | Pole-orbit upper bound | Archimedean lower bound |
| --- | ---: | ---: | ---: |
| linear, weight 0 | 1 | 40 | 1 |
| odd, weight 1 | 12 | 60 | 6 |
| linear, weight 2 | 66 | 90 | 33 |
| odd, weight 3 | 220 | 100 | 110 |
| linear, weight 4 | 495 | 120 | 248 |
| odd, weight 5 | 792 | 120 | 396 |
| linear, weight 6 | 462 | 100 | 231 |
| full translations | 1024 | 40 | 512 |

Thus archimedean places alone exceed the pole-orbit count in degrees
`220`, `495`, `792`, `462`, and `1024`, independently of all boundary
cocycles.  Ordinary Runge's required strict pole surplus is impossible in
those cases.

There is also a conditional closure of degree 66.  An unramified rational
prime with Frobenius cycle lengths `l_1,...,l_c` in `S_12` has

```text
sum_i floor(l_i/2)+sum_(i<j)gcd(l_i,l_j)             (7)
```

cycles on two-element subsets, hence that many places in the degree-66
field.  Formula (7) is at least six.  If `o` cycle lengths are odd, its
first sum is `(12-o)/2`.  If there is one cycle then `o=0`; otherwise the
cross terms are at least `c(c-1)/2>=c/2>=o/2`.  This proves the bound.

Consequently, if the target has ten distinct good denominator primes, one
for each triangle partition and unramified in the controlled fields, they
supply at least `10*6=60` finite places of the degree-66 lift field.  With
at least 33 archimedean places, the required set has size at least 93,
strictly exceeding the pole bound 90.  This is the same substantive
good-prime hypothesis used in the triangle-pole audits; the full fair
profile alone has not been proved to supply it.

Only the degree-one and degree-twelve target lifts remain outside these
comparisons.  Their pole budgets are 40 and 60, while archimedean places
supply only 1 and 6.  Boundary cocycles may reduce those pole counts, but
(3) shows that the target class does not determine such a reduction.
No blanket no-surplus statement is justified for these two cases without
additional joint division or finite-place information.

These statements concern the indicated cover over `Q`. They do not assume
that the rational base point has controlled height. A controlled
Weierstrass base obtained by adjoining a branch root changes the base
field and its linear Galois image; the full `S_12` calculations cannot be
transferred to that base unchanged. The coefficient-height and descent
issues for that construction are kept separately in
`six_point_two_cover_basepoint_audit.md`.
