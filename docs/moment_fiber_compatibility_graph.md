# Finite compatibility graphs for odd-moment fibers

Fix distinct positive numbers

```
0 < x_1 < ... < x_m
```

and consider sign vectors `e in {+1,-1}^m` for which the odd moments

```
sum_j e_j x_j^(2r-1),     1 <= r <= s
```

are equal.  If `e` and `f` are two different vectors in the same fiber,
their difference has signed support

```
T = {j : e_j != f_j},       b_j = e_j x_j  (j in T).
```

The Newton-identity argument gives `|T| >= 2s+1`: for `|T| <= 2s`, all
odd elementary coefficients through the relevant degree vanish.  Odd
support would force the nonzero product of the roots to vanish; even
support would make the root polynomial even and hence force opposite
roots, impossible for distinct absolute values.

At the two smallest possible supports, the
[critical and successor root-order lemmas](critical_odd_moment_root_order.md)
give additional necessary sign patterns. Read the
entries below in increasing order of the coordinates in `T`; global sign
negation is allowed.

```
t = 2s+1:  P = ++--++--...       (truncate at t)

t = 2s+2:  A = +--++--++...      (truncate at t)
            B = +++--++--...      (truncate at t)
```

For a fixed `(m,s)`, define the finite graph `G(m,s)` as follows.
Its vertices are all `2^m` sign vectors.  Two vertices `e,f` are adjacent
exactly when, with `T={j:e_j != f_j}` and `t=|T|`,

* `t <= 2s` is rejected;
* `t=2s+1` is accepted only when `(e_j)_{j in T}` is `P` or `-P`;
* `t=2s+2` is accepted only when it is `A`, `-A`, `B`, or `-B`;
* `t > 2s+2` is accepted without further restrictions.

The last rule is deliberate: this checker records only the two proved
small-support root-order restrictions.  Larger supports can satisfy
additional variation constraints, so `G(m,s)` is a relaxation.  Every
actual odd-moment fiber injects into a clique of `G(m,s)`, for arbitrary
distinct positive nodes.  Thus a maximum-clique computation gives a
uniform necessary upper bound, while a graph clique need not be realized
by any choice of the `x_j`.

The [exact exhaustive checker](check_moment_fiber_compatibility_graph.py)
computes all adjacency masks and uses a greedy-color branch-and-bound
maximum-clique search.  It verifies the following maxima:

```
G(6,1):  5
G(10,2): 7
G(8,2):  3
```

The `(10,2)` value is the relevant `m=4s+2` case in this small range.
The `(8,2)` computation is included as a nearby even-dimensional check.
These are finite necessary-code bounds only; they do not prove that a
fiber of size five or seven exists for positive nodes.

In particular, eight sign rows on ten distinct positive nodes cannot share
their first and third moments. This conclusion requires no orthogonality
or Hadamard completion. It remains a bound at these fixed dimensions;
the computations do not bound fibers uniformly as `m` and `s` grow.
