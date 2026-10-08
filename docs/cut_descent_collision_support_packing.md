# Collision graphs after commuting cut descents

Commuting opposite-twist descents retain more information than the total
factor `Q`.  Every collision edge has an exact support in the selected cuts.
Its source chord gives a height inequality for the product of the supported
cut factors.  Fractionally packing several collision supports yields a
radius-independent exclusion as soon as the packing mass reaches two at
arc constant `C<2`.

This is an application of the elementary common-unit support estimate from
[the opposite-phase exchange note](opposite_phase_exchange_graph.md), now
with collision histories supplying the supports.  It permits arbitrary
nested prime powers, shared rational primes among cut factors, and arbitrary
original row units.  It does not exclude every pair of collisions: strongly
overlapping supports with different signed histories can have packing mass
below two.

## Exact collision-edge factorization

Start with distinct Gaussian integers `z_1,...,z_m` of common norm `N` in
an angular interval of width

```text
Delta <= C N^(-1/4).
```

Perform `r` commuting positive-gap descents.  Let their oriented Gaussian
factors be `alpha_1,...,alpha_r`, with

```text
K_l=Norm(alpha_l),       Q=product_l K_l,       Q|N.
```

For row `j`, record `epsilon_jl=1` if it was divided by `alpha_l` and
`epsilon_jl=-1` if it was divided by `conjugate(alpha_l)`.  Then

```text
beta_j=product_l alpha_l^((1+epsilon_jl)/2)
                 conjugate(alpha_l)^((1-epsilon_jl)/2),
z_j=beta_j w_j,
Norm(beta_j)=Q,         Norm(w_j)=N/Q.                 (1)
```

No coprimality between different `alpha_l` is assumed.  The exact
cut-deletion theorem ensures that these are the literal factors removed in
any order.

Let `ij` be a collision edge: `w_i=w_j=w`.  Its unsigned cut support is

```text
T_ij={l:epsilon_il!=epsilon_jl},
K(T_ij)=product_(l in T_ij) K_l.                       (2)
```

The common-history factors give an exact factorization

```text
beta_i=H A,       beta_j=H conjugate(A),
Norm(H)=Q/K(T_ij),       Norm(A)=K(T_ij).              (3)
```

Here `A` chooses from `alpha_l,conjugate(alpha_l)` according to row `i`
on the differing cuts.  Since the original rows are distinct,
`A-conjugate(A)` is nonzero.  It is twice an imaginary integer, so

```text
z_i-z_j=H(A-conjugate(A))w,
abs(z_i-z_j) >= 2 sqrt(N/K(T_ij)).                     (4)
```

The source chord is at most `sqrt(N) Delta`.  Therefore every collision
edge satisfies

```text
K(T_ij) >= 4/Delta^2.                                  (5)
```

This retains the individual cut factors rather than only the cumulative
pair threshold `Q>=4/Delta^2`.

## A joint packing inequality

Give collision edges nonnegative weights `q_ij` such that every selected
cut has load at most one:

```text
sum_(ij:l in T_ij) q_ij <= 1             for every l,
P=sum_ij q_ij.                                           (6)
```

Multiplying (5) with these weights gives the exact joint inequality

```text
(4/Delta^2)^P
 <= product_l K_l^[sum_(ij:l in T_ij)q_ij]
 <= Q <= N.                                             (7)
```

At the endpoint scale this becomes

```text
(4/C^2)^P N^(P/2) <= N.                                (8)
```

Consequently, when `C<2`, the collision-support hypergraph has fractional
packing number strictly less than two.  In particular, two collision edges
with disjoint cut supports are impossible at `C<2`.  This is a
radius-independent restriction beyond applying the pair-distance lower
bound to each edge separately.

The condition is not automatic from having two collisions.  For example,
two supports `{1}` and `{1,2}` have packing optimum one.  Three supports
`{1,2},{1,3},{2,3}` have packing optimum `3/2`.  Thus (7) leaves genuine
overlapping collision graphs.

## Repeated signed histories give a forbidden rectangle

Orient a collision edge and record its signed difference word

```text
d_ij,l=(epsilon_il-epsilon_jl)/2 in {-1,0,1}.
```

Equation (1) gives

```text
z_i/z_j=beta_i/beta_j
       =product_l (alpha_l/conjugate(alpha_l))^d_ij,l. (9)
```

If two distinct collision edges `ij` and `kl` have the same signed word,
then

```text
z_i/z_j=z_k/z_l,             z_i z_l=z_j z_k.          (10)
```

Distinct edges with the same nonzero word have four distinct endpoints,
even when history cells are not singletons. A shared corresponding endpoint
(`i=k` or `j=l`) makes the other endpoints equal by (9), so the edges would
be identical. A shared opposite endpoint, say `i=l`, gives
`2 epsilon_i=epsilon_j+epsilon_k` coordinatewise. Since every coordinate
is `1` or `-1`, all three coordinates must agree, forcing the word to be
zero. The case `j=k` is the same. A zero word cannot support a collision
between distinct source points, because it gives `beta_i=beta_j` in (1).
The
[multiplicative rectangle theorem](multiplicative_rectangle_separation.md)
excludes (10) whenever

```text
Delta N^(1/4) < 2 sqrt(2).                              (11)
```

This condition also guarantees the theorem's minor-arc hypothesis:
`N>=1` implies `Delta<2 sqrt(2)<pi`.

Thus at `C<2 sqrt(2)` the evaluated signed words on collision edges are
all distinct.  There are only `(3^r-1)/2` words up to reversing an edge, so

```text
sum_over_collision_classes binom(class_size,2)
 <= (3^r-1)/2.                                         (12)
```

For a single descent every cross-collision edge has the same signed word.
Division is injective within either side, so two distinct collision edges
would have four distinct endpoints.  Hence one positive-gap descent from
the original target arc can create at most one collision when
`C<2 sqrt(2)`.  This conclusion uses the joint product relation, not just
the two individual chord bounds.

## A collision triple

If three rows collide to the same `w`, their three `beta` values are
distinct points of norm `Q` in an angular interval of width `Delta`.
Their real argument lifts are obtained by subtracting one fixed argument
of `w` from the source lifts, without choosing principal arguments anew.
All `beta_j` are congruent modulo `2`, since every `alpha_l` is congruent
to its conjugate modulo `2`.  The two difference rows in their affine
determinant are therefore divisible by two. Three distinct points on a
circle are noncollinear, so this determinant is nonzero even if the triple
is not Gaussian-primitive. Hence

```text
abs(D_beta)>=4.
```

Multiplication by `w` scales the determinant by `Norm(w)=N/Q`.  The usual
short-arc triangle upper bound now gives

```text
4N/Q <= abs(D_z) <= N Delta^3/8,
Q >= 32/Delta^3.                                       (13)
```

At `Delta<=C N^(-1/4)`, this is
`Q>=32 C^(-3)N^(3/4)`.  It is a stronger lower threshold than a single
collision, but by itself is compatible with `Q<=N` for large `N`.

## Actual six-row boundary fixture

The ordered primitive tuple

```text
N=14365,
z=((-107,-54),(-98,-69),(-91,-78),
   (-78,-91),(-69,-98),(-54,-107))
```

has a cut `S={0,1,3}`, `S^c={2,4,5}` with

```text
alpha=4-i,       K=17.
```

Its one-step descent is

```text
((-22,-19),(-19,-22),(-26,-13),
 (-13,-26),(-22,-19),(-19,-22)).
```

There are two collision pairs, `0~4` and `1~5`.  They have the same signed
one-cut word, and exactly

```text
z_0 z_5=z_1 z_4.                                       (14)
```

This tuple has normalized arc constant about `6.963`, above the rectangle
threshold.  It shows that the repeated-word relation is a real arithmetic
possibility and that the small source arc is essential.  It is not a
counterexample in the `C=1/2` class.

The [checker](check_cut_descent_collision_support_packing.py) verifies the
history factorization with cut factors sharing rational primes, the parity
and determinant identities for a collision triple, the elementary packing
examples, and the exact six-row fixture.
