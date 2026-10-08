# All-anchor deficit products and their exact common normalization

Coupling actual anchors supplies many automatic square products, but
their squareclasses come from one shared row-label space. There is also
an exact normalization obstruction: for a Gaussian-primitive tuple of
at least three points, the common affine divisor of **all** raw deficit
roots is only one or two, regardless of individual anchor contents.

The Eulerian character structure itself is already present in
[graph traces and spin defects](graph_trace_and_spin_defects.md), and
pair-label sums appear in
[pair-squareclass injectivity](endpoint_pair_squareclass_injectivity.md).
The statements below specialize these to integer raw deficits, a common
monic Runge parameter, and their exact joint affine gcd. They do not
prove a new general endpoint count.

## 1. One parameter and exact cycle square roots

Let `z_0,...,z_(m-1)` be distinct Gaussian integers of common norm `N`.
Choose ordinary primitive half-angle rows `h_i=a_i+ib_i` relative to
the actual anchor `z_0`, with `h_0=1` and

```text
z_i/z_0=h_i/conjugate(h_i),       f_i=a_i^2+b_i^2.
```

For an oriented edge `ij` set

```text
Delta_ij=a_i b_j-a_j b_i,       A_ij=a_i a_j+b_i b_j,
c_ij=N-Re(conjugate(z_i)z_j),       Z=2N.
```

Direct multiplication gives, without any reduced-pair parity convention,

```text
c_ij=Z Delta_ij^2/(f_i f_j),
Z-c_ij=Z A_ij^2/(f_i f_j),
c_ij(Z-c_ij)=Im(conjugate(z_i)z_j)^2.               (1)
```

All `c_ij` are positive integers. For a graph edge set `E` with even
cardinality and even degree at every vertex, (1) gives the explicit
rational square roots

```text
product_(ij in E)c_ij
 =[Z^(|E|/2) product_(ij in E)Delta_ij
                 /product_i f_i^(deg_E(i)/2)]^2,
product_(ij in E)(Z-c_ij)
 =[Z^(|E|/2) product_(ij in E)A_ij
                 /product_i f_i^(deg_E(i)/2)]^2.    (2)
```

Each product is an integer, so its rational square root is an integer.
The signs chosen on edge orientations disappear after squaring. In
particular every even cycle gives a monic integer-square product at
`Z=2N`, without a squareclass pigeonhole argument. Distinct edges need
not have distinct deficits; distinct roots must still be checked before
applying the [even-product bound](runge_even_product_square_bound.md).

## 2. The shared label rank, with a sharp primitive realization

Write squareclasses additively and augment each edge label with a final
coordinate one. If

```text
v_i=([f_i],0),       v_0=0,       w=([2N],1),
```

then (1) says exactly

```text
([c_ij],1)=v_i+v_j+w.                              (3)
```

Thus the entire complete graph has augmented rank at most `m`, rather
than one independent row-label space for each anchor. With independent
formal `v_1,...,v_(m-1),w` and `m>=3`, its rank is exactly `m`: a triangle
sums to `w`, and the star edges then recover all `v_i`. The formal kernel
has dimension `binom(m,2)-m`. Its elements are exactly the even-cardinality
Eulerian edge sets. Actual row-label dependencies can enlarge this kernel.

The rank bound is attained on primitive integer circles. For `m-1` rows
choose even integers recursively by

```text
t_1=2,       t_(i+1)=2 product_(j<=i)(1+t_j^2),
h_i=1+i t_i,       f_i=1+t_i^2,
z_0=product_i conjugate(h_i),
z_i=h_i product_(j!=i)conjugate(h_j),
N=product_i f_i.                                  (4)
```

The `f_i` are odd, pairwise coprime, and nonsquares: each lies strictly
between `t_i^2` and `(t_i+1)^2`, and every later one is one modulo every
earlier one. Their squareclasses are independent. Each `h_i` is coprime
to its conjugate; its prime block is flipped only at row `i`, so the
whole physical tuple is Gaussian-primitive. Its points are distinct
because the positive `t_i` are distinct. Formula (3) therefore has rank
exactly `m`. No factorization of the large integers in (4) is needed.
This is a label-rank fixture, not an asserted endpoint family.

## 3. The all-edge affine divisor is exactly one or two

Assume now that the full tuple is Gaussian-primitive and `m>=3`. Then
`N` is odd. Define

```text
G_plus=gcd_(i<j)(N+Re(conjugate(z_i)z_j)).
```

We have the exact classification

```text
G_plus=2 if all points have the same coordinate-parity orientation,
G_plus=1 otherwise.                               (5)
```

The two orientations mean odd real/even imaginary or even real/odd
imaginary coordinates.

For an odd prime `p|N`, Gaussian primitivity makes its allocation range
attain both endpoints. A pair at opposite endpoints has a scalar product
which is a unit modulo `p`: modulo either Gaussian orientation, exactly
one of the two conjugate products is a unit. Thus `p` does not divide
`G_plus`, at any source-prime exponent.

If an odd prime `p` absent from `N` divided `G_plus`, any three actual
coordinate vectors would have a Gram matrix modulo `p` with diagonal
`N` and off-diagonal entries `-N`. Its determinant is `-4N^3`, nonzero
modulo `p`. This contradicts the rank-at-most-two Gram representation.

Finally, primitive odd `N` is one modulo four. Points in different
coordinate-parity orientations give an odd complement. If all points
have one orientation, every complement is even. Among three points,
two odd coordinates agree modulo four, making their scalar product one
modulo four and their complement two modulo four. This proves (5),
including its exact two-adic exponent.

For a general tuple, let `G` be its Gaussian gcd and divide all rows by
`G`. Scalar products and `N` both scale by `Norm(G)`. Consequently the
general value of `G_plus` is exactly `Norm(G)` times the one-or-two
factor from the primitive tuple.

For any reference edge `e_0`, the maximal integer affine normalization
of all roots `c_e` and the parameter `Z` has effective divisor

```text
gcd(Z-c_(e_0), {c_e-c_(e_0):e in E(K_m)})
 =gcd_e(Z-c_e)=G_plus.                             (6)
```

Thus jointly using all anchors removes the possible large content
available at a single anchor. A larger gcd of root differences alone
gives a rational parameter; clearing its denominator restores exactly
the effective divisor (6), as in the
[direct affine-normalization audit](direct_projection_affine_normalization_height_tradeoff.md).

## 4. Cycle compatibility does not force local remainder divisibility

The [local remainder audit](direct_monic_remainder_local_obstruction.md)
also survives the common-circle constraint, not just separate norm
equations. Let an Eulerian simple graph have `2k` edges, and isolate a
nonisolated vertex `v` of even degree `l`. Then `2<=l<=2k-2`.
At a split prime `p^e||N`, use local isotropic coordinates
`alpha_i beta_i=N` and units `t_i`, choosing

```text
alpha_v=N t_v,       beta_v=1/t_v,
alpha_i=t_i,         beta_i=N/t_i       (i!=v).
```

These are points on one integral local circle with a binary allocation
cut. For the crossing edges,

```text
c_vi=N-(N^2 t_v/t_i+t_i/t_v)/2
     =-t_i/(2t_v)  (mod p^e).
```

The `l` unit roots are freely chosen through the distinct neighbors'
units `t_i`. Every other edge has deficit
`N[1-(t_i/t_j+t_j/t_i)/2]`, divisible by `p^e`. Thus the mixed
occupancy is even. For split `p>8k`, the variable-degree construction
in the local audit can make
`F_k(U)=[t^k]product_(u in U)(1-u t)^(1/2)` a unit. Its exact formula
then gives `R(2N)=-[2^(2k-1)F_k(U)]^2 mod p^e`, a unit remainder.
Changing orientation swaps the isotropic coordinates and preserves
this conclusion. This construction is local; global endpoint
compatibility remains a separate question.

A small global check uses the cycle
`(-8,1),(-8,-1),(-4,7),(-7,4)` on `N=65`. Its four successive
deficits are `2,40,9,5`, whose product is `60^2`; their complementary
product at `Z=130` is `13200^2`. Nevertheless the cleared quartic
Runge remainder is `-41433616`, a unit modulo five. Thus exact cycle
square products do not themselves supply the missing divisibility.
This is a finite circle fixture, not a large endpoint family.

## 5. What this gives, and what remains missing

On an endpoint arc all roots satisfy `c_ij<=C^2 sqrt(N)/2`. For a
primitive tuple with `C<=sqrt(2)`, the large-cohort diameter result in
the affine-normalization audit gives an actual stronger statement.
For any `epsilon>0`, if

```text
ceil((m-1)/2)>3/(4epsilon^2)+4,
```

the all-edge root set has diameter at least `N^(1/2-epsilon)`: it
contains the star at any anchor, to which that result applies.
After maximal integer affine normalization, (5)--(6) therefore leave
diameter at least

```text
(1/2)N^(1/2-epsilon).                              (7)
```

This is independent of the individual contents of the anchors. It
prevents the complete all-anchor raw-deficit model from acquiring
`N^(O(1/m))` root height through one common affine normalization.

The cycle square products in (2) are exact and useful identities, but
they neither add independent prime-support labels nor supply the small
root height used in the conditional reciprocal-map theorem. A selected
small graph or cohort can have a different gcd; further identities or
special arithmetic cancellation remain outside this obstruction. The
general endpoint goal remains unproved.

The [checker](check_all_anchor_deficit_cycle_squareclasses.py) verifies
the exact edge and cycle square roots, formal graph ranks and kernels,
primitive rank fixtures, and the one-or-two all-edge gcd on literal
circles, including nonprimitive scaling.

## 6. Connected sparse graphs in the binary case

Assume now that the whole tuple is Gaussian-primitive and its split-prime
allocations are binary: at every `p^e||N`, each point has allocation
`0` or `e`. Thus `N` is odd, every point has ordinary coordinate content
one, and each prime cut is nonconstant. For a connected simple spanning
graph `E` on at least three distinct points, put

```text
G_E = gcd_(ij in E)(N+Re(bar(z_i)z_j)),
L_E = gcd_(ij in E)(Re(z_i+z_j), Im(z_i+z_j)).
```

These gcds are positive: a connected graph on at least three distinct
points cannot have every edge antipodal. Antipodal edges contribute zero
and are ignored in valuation minima.
The effective integer affine divisor of all its raw deficits is exactly
`G_E`, by the identity (6). The following exact formulas account for the
ramified prime as well as both split-prime orientations:

```text
mixed coordinate-parity orientations:    G_E=L_E^2=K^2,   K=L_E;
uniform coordinate-parity orientation:   G_E=L_E^2/2=2K^2, K=L_E/2.       (8)
```

Here mixed means that some points have odd first coordinate and others
have odd second coordinate. To prove (8), first consider an odd `p|N`.
Connectivity forces an edge across its nonconstant binary allocation
cut. At that edge the complementary deficit is a p-unit, so `p` divides
neither `G_E` nor `L_E`. This argument is unchanged on interchanging the
two Gaussian prime orientations.

For odd `p` not dividing `N`, use the integral p-adic orthogonal basis
`z_i, i z_i` and write `z_j=a z_i+b i z_i`, with `a^2+b^2=1`.
If `p|(1+a)`, then `a-1` is a unit and
`v_p(1+a)=2v_p(b)`. If `1+a` is a unit, both relevant minimum valuations
are zero. Consequently, with the usual infinite valuation at zero,

```text
v_p(N+Re(bar(z_i)z_j))
  =2 min(v_p(Re(z_i+z_j)), v_p(Im(z_i+z_j))).            (9)
```

At two, a mixed-orientation edge has complementary deficit odd. On a
same-orientation edge write `z_i+z_j=2^r u`, where `r>=1` and `u` has
at least one odd coordinate. Its norm is odd: otherwise both coordinates
of `u` are odd, and the identity
`Norm(z_i+z_j)=2 Re(bar(z_i)(z_i+z_j))` would have valuations `2r+1`
and `r+1`, respectively, which is impossible. Hence the complementary
deficit has valuation `2r-1`. Connectivity now proves (8).

If `E` is nonbipartite, an odd cycle and its edge-sum congruences imply
`L_E|2z_i` for every vertex. Whole-tuple primitivity implies `L_E|2`.
Thus

```text
connected nonbipartite spanning graph:    K=1, G_E in {1,2}.             (10)
```

If `E` is bipartite, all vertices in either part are congruent modulo
`L_E`, and the two parts have opposite common residues. For two distinct
vertices in one part, let `H=a+i b` be their ordinary primitive
half-angle representative, of norm `f=epsilon n`, where `epsilon` is
one or two and `n|N` is the odd primitive pair norm. Then

```text
K divides b.                                                        (11)
```

Indeed odd primes of `K` are absent from `N`, so are absent from `f`;
the identity `z_j-z_i=2i z_i b/bar(H)` proves the odd-prime divisibility.
If `K` is even, all points have the same parity orientation. Such a pair
has `f` odd: a half-angle with both coordinates odd would interchange
the physical parity orientation. The factor `z_i/bar(H)` is then a
2-adic unit vector. Since `2K=L_E` divides the physical difference,
(11) also holds at two.

Suppose additionally that all points lie in an arc of length
`C N^(1/4)` with `C<=sqrt(2)` and `N>1`. For any part of size `q>=2`,
the pair chord bound gives

```text
log|b_ij| <= -W mu_ij/4,
W=log N,  mu_ij=1-2 log(n_ij)/W.
```

The allocation-sign Gram matrix, formed at the original radius `N`,
is positive semidefinite with diagonal one. Thus
`sum_(i<j)(-mu_ij)<=q/2`. This is the same all-unit residue budget used
in the [affine-normalization audit](direct_projection_affine_normalization_height_tradeoff.md),
now applied to the common divisor (11). Multiplication over all pairs
in the part proves

```text
K <= N^[1/(4(q-1))],     G_E <= 2 N^[1/(2(q-1))].                      (12)
```

One may take the larger part, so `q>=ceil(m/2)`. The selected part is
not primitively rescaled: its Gram retains the original prime weights.
This extends the common-normalization obstruction to every connected
spanning graph in the binary case. It does not assert a root-diameter
bound for an arbitrary sparse graph.

For the literal cycle in geometric arc order, a diameter bound is
available. Suppose `m>=4` is even, and the large-cohort hypothesis in
Section 5 holds. The closing edge has maximum raw deficit `c_max`, and
some consecutive angular gap is at most the total angle divided by
`m-1`. On a minor arc the sine bounds give
`c_min/c_max<=4/(m-1)^2<=4/9`; hence the cycle root diameter is at least
`c_max/2`. Section 5 implies `c_max>=N^(1/2-epsilon)`. Its bipartition
has `q=m/2`, so (12) leaves normalized diameter at least

```text
(1/4) N^[1/2-epsilon-1/(m-2)].                                       (13)
```

In this small-arc regime all unordered pair deficits are distinct.
Indeed equal chord lengths on a minor arc give equal positive angular
gaps and hence a nontrivial equality of two point products, including
possible repeated factors. The all-unit
[multiplicative rectangle theorem](multiplicative_rectangle_separation.md)
excludes this whenever `C<2sqrt(2)`. Thus (13) applies to the actual
distinct-root cycle polynomial; there is no repeated-factor cancellation
in its stated regime.

Outside this regime, removing repeated factors of even multiplicity can
disconnect the retained graph and change its gcd. For example, the geometrically ordered cycle

```text
(-8,-1), (-7,-4), (-4,-7), (-1,-8), (4,-7), (7,-4),   N=65
```

is Gaussian-primitive and binary. Its deficits are
`[5,9,5,13,9,117]` with complementary gcd one. Removing the two repeated
values leaves `[13,117]`, with complementary gcd thirteen. Its angular
span is approximately `2.4981`, giving endpoint constant approximately
`7.0931`, so this finite fixture is not in the
`C<=sqrt(2)` regime. It demonstrates the algebraic cancellation caveat,
not a counterexample to (12) or (13). Outside that small-arc regime, no
obstruction for the subsequently renormalized squarefree kernel is claimed.
