# Maximal-gap descent for arbitrary bipartitions

The oriented opposite-twist descent extends from the ordered outer-pair cut
to every bipartition of a primitive equal-norm Gaussian tuple.  The product
of all bipartition gap factors is exactly the squared radius.  A maximal
gap therefore gives a strict norm descent, and an original short arc becomes
two shorter clusters.  At the target constant this first descent is also
collision-free.

This supplies a genuine binary decomposition, but norm descent alone does
not close an induction for the number of points: the two child
cardinalities add, and subsequent global descents may compare points from
different clusters.  The precise missing input is stated below.

## All-cut product identity

Let `z_1,...,z_m` be distinct Gaussian integers of common norm `N`, primitive
as a full tuple.  As in the
[five-point descent](ordered_gap_opposite_twist_descent.md), primitivity
implies that `N` is odd and supported on split rational primes.  Fix one
factor `pi` above every such prime `p`.  Write

```text
e_p=v_p(N),             t_j(p)=v_pi(z_j).
```

For every `p`, full primitivity gives

```text
min_j t_j(p)=0,         max_j t_j(p)=e_p.              (1)
```

For an unordered nontrivial bipartition `S|S^c`, let

```text
U_S(p)=[min_(j in S)t_j(p),max_(j in S)t_j(p)],
h_p(S)=distance(U_S(p),U_(S^c)(p)),
K_S=product_p p^h_p(S).                                (2)
```

There are `q_m=2^(m-1)-1` such cuts.  They satisfy the exact identity

```text
product_(unordered cuts S|S^c) K_S = N.                (3)
```

Indeed, for each fixed `p` and each level `a=1,...,e_p`, the upper set

```text
S_a={j:t_j(p)>=a}
```

is nonempty and proper by (1), hence defines one unordered cut.  A level
`a` contributes to `h_p(S)` exactly when `S_a` is `S` or `S^c`: all rows
on one side then lie above that level and all rows on the other lie below
it.  Thus

```text
sum_(cuts S|S^c) h_p(S)=e_p,
```

which proves (3) prime by prime.  In particular some maximal cut obeys

```text
K_max >= N^(1/q_m).                                    (4)
```

For `N>1` it has `K_max>1` and therefore gives a strict norm reduction.

## Opposite-twist descent for a cut

Fix a cut with `K_S>1`.  At every prime with disjoint intervals, include
`pi^h` in `alpha_S` if `U_S` lies above `U_(S^c)`, and include
`conjugate(pi)^h` if it lies below.  Then

```text
Norm(alpha_S)=K_S,
alpha_S | z_j                 (j in S),
conjugate(alpha_S) | z_j      (j in S^c).              (5)
```

Divide accordingly:

```text
w_j=z_j/alpha_S                    (j in S),
w_j=z_j/conjugate(alpha_S)         (j in S^c).          (6)
```

The primewise interval proof from the five-point case applies verbatim:
the `w_j` have common norm

```text
M=N/K_S,                                               (7)
```

are primitive as a full tuple, and the chosen cut has zero gap after the
descent.

Suppose the original tuple lies on an arc of physical length

```text
L <= C N^(1/4).
```

Within either side of the cut, (6) is one common similarity of ratio
`K_S^(-1/2)`.  Each side is therefore contained in an arc of physical
length at most

```text
L/sqrt(K_S) <= C K_S^(-1/4) M^(1/4).                  (8)
```

Thus maximal-gap descent strictly reduces the norm and preserves all
labels in at most two short angular clusters, each with improved normalized
constant `C K_S^(-1/4)`.  The relative rotation between those clusters is
uncontrolled.

There is an equivalent independent-child version.  Let

```text
G_S=gcd_G(z_j:j in S),       A_S=Norm(G_S),
G_T=gcd_G(z_j:j in S^c),     A_T=Norm(G_T).
```

Equation (5) gives `A_S,A_T>=K_S`.  Dividing each child by its own gcd
produces two primitive one-arc configurations of norms at most `N/K_S`
and normalized constants at most `C K_S^(-1/4)`.

## The first descent is collision-free at the target constant

Opposite rotations can in principle identify a row from `S` with a row
from `S^c`.  On one sufficiently short source arc this cannot happen
unless both sides are singletons.

The factors in `alpha_S` use at most one of `pi,conjugate(pi)` over each
rational prime.  Hence

```text
gcd_G(alpha_S,conjugate(alpha_S))=1.
```

Write `alpha_S=a+ib`.  Since `K_S>1`, neither `a` nor `b` is zero; in
particular `|b|>=1`.  If a cross pair collides in (6), say both descended
rows equal `w`, then

```text
z_i=alpha_S w,       z_j=conjugate(alpha_S) w,
|z_i-z_j|=2 |b| sqrt(M) >= 2 sqrt(N/K_S).              (9)
```

Both original points lie on the source arc, so (9) and its length bound
give

```text
K_S >= 4 C^(-2) N^(1/2).                              (10)
```

Consequently every within-side arc after (6) has length

```text
L/sqrt(K_S) <= L^2/(2 sqrt(N)) <= C^2/2.              (11)
```

Two distinct equal-norm Gaussian integers have distance at least `sqrt(2)`.
If

```text
C < 2^(3/4),                                            (12)
```

then (11) is shorter than `sqrt(2)`, so each side contains at most one
row.  Therefore for `m>=3` no cross collision is possible.  In particular,
every positive-gap first descent of a tuple in the `C=1/2` arc class
retains all distinct points.

The conclusion concerns a descent from one common short arc.  It cannot be
reapplied blindly after (6), because the new tuple generally occupies two
separated clusters.

## Exact finite sequences and history cells

A descent has a stronger labeled invariant.  At a fixed prime, assume
`U_S` lies above `U_(S^c)`, with gap `h`.  Subtracting `h` from the
allocations on `S` deletes exactly the `h` threshold layers whose upper set
is `S`.  Every other threshold upper set and its multiplicity is unchanged.
The reverse orientation is identical after complementation.  Therefore a
descent at the cut `S|S^c` preserves the primewise size and orientation of
every other positive gap and changes the all-cut factors by

```text
K_S -> 1,             K_T -> K_T for every T!=S.       (13)
```

No labeled cut factor is created or enlarged by a previous descent, so the
positive cuts may be selected in any order.

Consider `r` such descents, with factors `alpha_1,...,alpha_r`, and put

```text
Q=product_l Norm(alpha_l)=product_l K_l.
```

For each original row `j`, let `beta_j` be the product of `alpha_l` or
`conjugate(alpha_l)` according to the side containing `j` at step `l`.
Then

```text
z_j=beta_j w_j,          Norm(beta_j)=Q.                (14)
```

Moreover `alpha_l` and `conjugate(alpha_l)` are congruent modulo the
rational Gaussian integer `2`.  Hence all `beta_j` are congruent modulo
`2`.  If two distinct original rows collide at the end, so that
`w_i=w_j=w`, then

```text
z_i-z_j=(beta_i-beta_j)w,
abs(beta_i-beta_j)>=2,
Q >= 4 C^(-2) N^(1/2).                                 (15)
```

Rows with one common sequence of side choices form a history cell.  There
are at most `2^r` cells.  Within a cell, `beta_j` is constant and divides
every original pair difference.  If a cell contains two distinct rows,
Gaussian separation and the source arc bound give

```text
2Q <= abs(z_i-z_j)^2 <= L^2 <= C^2 N^(1/2),
Q <= (C^2/2) N^(1/2).                                  (16)
```

The upper bound (16) and collision lower bound (15) contradict each other
exactly when `C^4<8`, equivalently `C<2^(3/4)`.  Therefore a final collision
forces every history cell to be a singleton and hence

```text
m <= 2^r.                                               (17)
```

A collision at an intermediate step persists: equal labeled rows cannot
lie on opposite sides of any later positive-gap cut, because their common
allocation makes the two cut intervals meet at every prime.  Thus if
`m>2^r`, the entire `r`-step sequence is collision-free.  More locally, a
final collision is impossible whenever even one history cell still
contains two rows.

There is a critical-exponent consequence which indicates what an induction
would have to improve.  Fix three labels.  Descend through every positive
cut which does not split that triple.  By (13), the accumulated factor is
the product of the original `K_S` over precisely those cuts.  The triple
stays in one history cell and remains a distinct short-arc triple.  Its
source angular span `Delta=L/sqrt(N)` is unchanged by the common division.
Its nonzero even determinant gives

```text
Q_triple <= N Delta^3/16 <= (C^3/16) N^(1/4).          (18)
```

This accumulated factor has a familiar exact meaning.  At each prime, the
threshold layers which do not split the triple are precisely those below
its minimum allocation or above its maximum allocation.  They form its
constant Gaussian layers.  Therefore

```text
Q_triple=Norm(gcd_G(z_i,z_j,z_k)).                     (18a)
```

The sequential construction is an operational realization of ordinary
subgroup Gaussian normalization here.  Equations (18)--(20) rephrase its
usual cut-counting average and do not create a new exponent.

Average `log Q_triple` over all triples.  A cut of sizes `s,m-s` leaves
exactly `binom(s,3)+binom(m-s,3)` triples unsplit.  Using (3), some triple
therefore has

```text
Q_triple >= N^rho_m,
rho_m=min_s [binom(s,3)+binom(m-s,3)]/binom(m,3).       (19)
```

Writing `m=2k` or `m=2k+1`, respectively,

```text
rho_(2k)   =(k-2)/(2(2k-1)),
rho_(2k+1) =(k-1)/(2(2k+1)).                           (20)
```

These exponents approach `1/4` from below.  Thus descent plus an unweighted
all-cut average reaches the critical exponent of (18), but never exceeds
it.  A viable uniform induction along this route would need a geometric
bias in the gap weights which improves (19), or a stronger joint estimate
than applying (18) to one history triple.

## Why the binary induction does not yet close

Equations (4), (7), and (8) give a rigorous maximal-gap algorithm.  At a
node containing `s` points it finds a cut with

```text
K >= N_node^(1/(2^(s-1)-1))                            (21)
```

and replaces the node by two independently primitive, shorter-constant
children of smaller norm.  This is a valid induction on the norm for any
statement already known to be subadditive under disjoint unions.  A uniform
point bound is not subadditive with the same constant: even if each child
has at most `B` points, their union can have `2B`.

The gain in (21) also deteriorates exponentially with `s`.  Pure numerical
iteration therefore allows arbitrarily large binary trees when the root
norm is chosen sufficiently large.  The improvement of the normalized
constant does not supply a terminal case: for every positive constant and
arbitrarily large norms, a short arc can still contain two points.  An
elementary primitive example is

```text
z_-=t-i,       z_+=t+i,       N=t^2+1                 (t even).
```

Their minor-arc length is
`2 sqrt(t^2+1) arctan(1/t)`, so their normalized constant is

```text
2 (t^2+1)^(1/4) arctan(1/t) -> 0.
```

Thus no positive normalized-constant threshold makes the leaves
singletons.

Keeping the children on one common circle does not fix this.  A second
global cut can mix the two clusters, and a collision may then involve
points in different clusters, for which the source-arc estimate in (9) is
unavailable.  This phenomenon occurs in exact arithmetic.  The ordered
`N=1105,K=5` fixture from the five-point note first descends to

```text
((-14,5),(11,10),(10,11),(5,14),(-14,-5)),            (22)
```

of norm `221`.  For the next cut

```text
S={0,2,3},       S^c={1,4},       alpha=4+i,
```

the gap is `K_S=17`, the new norm is `13`, and the rows at indices one and
three both descend to `2+3i`.  This fixture has a much larger normalized
constant than `1/2`; it demonstrates the exact recombination mechanism,
not a counterexample in the target arc class.

Thus the maximal-gap construction does preserve a bounded number of short
clusters in one step and gives a strict norm descent.  To turn it into a
uniform proof one still needs a new statement controlling recombination
across previously separated clusters, or a potential whose decrease pays
for the addition of the two child cardinalities.  Neither the all-cut
product (3) nor the normalized-length gain (8) supplies that statement by
itself.

The [checker](check_general_bipartition_opposite_twist_descent.py) exhausts
the all-cut product and cut-deletion identities for small allocation tables,
verifies random multi-prime descents and subgroup divisibility, checks the
triple-averaging constants, and verifies the exact two-stage history
collision fixture.
