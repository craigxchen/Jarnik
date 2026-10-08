# Positive squareclasses: a pair bound for integer-radius circles

Assume the radius is an integer `R`. The squareclass decomposition then has a
useful extra feature: the nonsquare factor is a positive squarefree rational
integer. This note proves a pairwise restriction on distinct squareclasses in
one unit class. It is a genuine consequence of the common angular interval,
but it does not bound the number of low squareclasses.

## 1. Factorization with positive squarefree radicand

If `z in Z[i]` satisfies `N(z)=R^2`, unique factorization gives

```text
z = epsilon d w^2,       epsilon in {+1,-1,+i,-i},    (1)
```

where `d` is a positive squarefree rational integer, `d|R`, and

```text
N(w)=R/d.                                                   (2)
```

For a split prime `p=pi bar(pi)` with `v_p(R)=e`, the two Gaussian
valuations of `z` sum to `2e` and hence have the same parity. If that parity
is odd, remove the positive factor `p=pi bar(pi)` from `z`; the remaining two
valuations are even. For an inert prime `p=3 mod 4`, the Gaussian valuation of
`z` is `v_p(R)`, so remove `p` exactly when this is odd. At `1+i`, the
valuation of `z` is `2v_2(R)`, already even. The unit is the remaining factor.
This proves (1)--(2), including ramified and inert primes.

## 2. Pair bound for distinct positive squareclasses

Let `z_1,z_2` be distinct points in an arc of angular width
`Delta<=C/sqrt(R)`, and suppose they have the same unit in (1):

```text
z_j = epsilon d_j w_j^2,       d_1 != d_2.                (3)
```

Then

```text
d_1 d_2 <= C^2 R/4.                                      (4)
```

### Proof

The vectors `w_1,w_2` cannot be real-collinear. Otherwise, after choosing a
primitive Gaussian generator of their common rational line, write
`w_1=a u`, `w_2=b u` with nonzero integers `a,b`. Equality of the two point
norms gives `d_1a^2=d_2b^2`; two squarefree positive integers in a rational
square ratio are equal, contrary to (3).

Put `B=det(w_1,w_2)`, a nonzero integer. The common unit and positivity of
the `d_j` imply that the angular difference of `w_1,w_2`, after changing one
root sign if necessary, is at most `Delta/2`. Therefore

```text
1 <= |B| = |w_1||w_2| |sin(arg(w_2)-arg(w_1))|
       <= (R/sqrt(d_1 d_2)) (Delta/2)
       <= C sqrt(R)/(2 sqrt(d_1 d_2)).                    (5)
```

Rearranging proves (4).

The estimate is sharp in scale. For `R=65`,

```text
z_1=5(3+2i)^2=25+60i,       d_1=5,
z_2=13(2+i)^2=39+52i,       d_2=13,
```

both have norm `R^2`, and `|det((3,2),(2,1))|=1`. Their angular
separation is exactly `2 asin(1/sqrt(65))`, approximately `0.24871`
radians. The normalized arc constant is
`2 sqrt(65) asin(1/sqrt(65))`, approximately `2.00516`, just above
the lower bound `2 sqrt(d_1 d_2/R)=2`. The small difference is the
sine-versus-angle slack.

## 3. Consequences and the exact remaining gap

Partitioning an integer-radius cluster by the four possible units and by
literal positive squarefree radicand gives the per-class root bound

```text
number in class (epsilon,d) <= 1+C/(2 sqrt(2d)),       (6)
```

from the square-root packing argument. The new pair bound says that, within
one unit class, at most one occupied radicand can exceed `C sqrt(R)/2`, since
two such radicands would have product greater than `C^2R/4`.

This does not control the number of small radicands, even for small fixed
`C`. Let `Q` be a product of arbitrarily many distinct odd split primes,
choose a further split prime `P>4Q/C^2`, and put `R=QP`. Every divisor
`d|Q` then satisfies `d<=Q<C sqrt(R)/2`, so every pair among these
`2^{omega(Q)}` radicands satisfies (4). Each class separately can be
represented on the integer-radius circle, since `R/d` is a product of
split primes and hence a Gaussian norm. This supplies no common short arc.
It shows precisely that positivity, individual realizability and the
pair inequality (4) by themselves do not bound the number of classes.
No impossibility claim is made about further determinant relations among
simultaneously chosen roots.

For nonintegral radius the decomposition (1) need not exist with an integer
`R`; this note makes no claim in that case. Even for integer radius, a further
cross-class relation among the `w_j` or among the positive radicands `d_j` is
needed to turn (4)--(6) into a radius-independent total count.
