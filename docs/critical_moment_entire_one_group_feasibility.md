# Entire residual identities allow arbitrarily many rows against one row

The constant-even cross identity and linear-odd within-group identity
have nondegenerate entire-function realizations with arbitrarily many
rows in one group and one row in the other. All source zeros are real,
simple, nonzero, and have pairwise distinct absolute values. This is
an analytic functional construction, not a finite orthogonal sign code,
not a `4+4` configuration, and not a lattice-circle construction.
It does not assert approximation by full finite moment templates.

## A real-rooted pencil

Set

```
B(z)=cos z+sin z,       C(z)=cos z-sin z,
A_t(z)=B(z)+t z C(z),          t>0.
```

Then `B(-z)=C(z)`, `C(-z)=B(z)`, and `B^2+C^2=2`.
Choose any distinct positive parameters `t_1,...,t_m` and put
`A_i=A_(t_i)`. All `A_i(0)=C(0)=1`.

Every zero of `A_t` is real and simple. Away from the real zeros of `C`,

```
A_t(z)/C(z)=tan(z+pi/4)+t z.
```

In the upper half-plane the imaginary part of both terms is positive:

```
Im tan(x+iy)=sinh(2y)/(cos(2x)+cosh(2y))>0,     y>0.
```

Thus the quotient has no nonreal zero, using conjugation for the lower
half-plane. At a zero of `C`, the value of `A_t` is `B!=0`. On each
real interval between successive zeros of `C`, the quotient increases
strictly from minus infinity to plus infinity, since its derivative is
`sec^2(z+pi/4)+t>0`. This proves simplicity and gives exactly one zero
per interval. The zeros of `C` itself are simple and real.

In fact the zeros of all the source functions `C,A_1,...,A_m` have
distinct absolute values. Two different `A` functions cannot share a
zero, since their difference is `(t_i-t_j)zC(z)` and their value at
zero is one. Nor can `A_i(z)` and `A_j(-z)` vanish at the same real
point: the corresponding equations in `B(z),C(z)` have determinant
`1+t_i t_j z^2>0`. This argument includes `i=j`.
Zeros shared with `C(z)` or `C(-z)=B(z)` are likewise impossible.
Finally `B,C` have no common zero. These observations also exclude
opposite zeros within `C`.

## Row factors and exact residual identities

Let `A={1,...,m}` and let `b` denote the single row in the other group.
Define entire row factors normalized to one at zero by

```
E_i(z)=A_i(-z) C(z) product_(j!=i) A_j(z),
E_b(z)=C(-z) product_j A_j(z).
```

Every row satisfies the same reflection-product identity

```
E_i(z)E_i(-z)=E_b(z)E_b(-z)
 =C(z)C(-z) product_j A_j(z)A_j(-z).                  (1)
```

Between row `i` and row `b`, the common factors are the other `A_j`.
The differing residual is `A_i(-z)C(z)`, and exactly

```
A_i(-z)C(z)+A_i(z)C(-z)=2.                           (2)
```

Indeed this residual is
`1-sin(2z)-t_i z cos(2z)`, visibly one plus an odd function.
Between rows `i,j` in `A`, the common factors are `C` and the other
`A` factors. The differing residuals satisfy

```
A_i(-z)A_j(z)-A_i(z)A_j(-z)=2(t_j-t_i)z.             (3)
```

Thus the cross even coefficients vanish and all within-group odd
coefficients above degree one vanish, just as in the normalized finite
residual identities. The differing residuals have only simple real
zeros with distinct absolute values by the preceding argument.

There is no coalescence in the row derivatives. Write
`beta_i=-E_i'(0)`, `beta_b=-E_b'(0)`, and
`T=sum_j(1+t_j)`. Direct differentiation gives

```
beta_i=3+2t_i-T,       beta_b=-1-T,
beta_i-beta_b=4+2t_i>0.
```

All `m+1` derivatives are therefore distinct. These derivative labels
are not presented as finite reciprocal sums; no choice of a canonical
product's real linear exponential factor is being suppressed.
The assertions above concern the displayed entire functions themselves.

## The common partner is unique within this pencil construction

Fix two distinct parameters `t_i,t_j`. If an entire function `D` obeys

```
A_i(-z)D(z)+A_i(z)D(-z)=2,
A_j(-z)D(z)+A_j(z)D(-z)=2,
```

then `D=C`. Indeed, subtraction gives
`C(z)D(-z)=B(z)D(z)`. Substituting this relation into either equation,
and multiplying by `C`, gives
`(C^2+B^2)D=2C`. Since `C^2+B^2=2`, the conclusion follows.
This proves that the displayed pencil cannot acquire a second distinct
partner satisfying these same factor-level equations. It does not
classify other factorizations or exclude a general `4+4` entire system.

## What is missing from a finite orthogonal template

The source functions describe sign choices on infinitely many real
zeros, but supply neither a finite common degree nor the required
punctured orthogonal row counts. One can see a count mismatch directly
for natural symmetric magnitude truncations. Each source factor has
`2R/pi+O(1)` zeros of absolute value at most `R`, since `A_t` has one
zero between successive zeros of `C`. There are `m+1` source factors,
and every row pair differs on exactly two of them. Thus its differing
fraction tends to `2/(m+1)`, which equals the orthogonal proportion
`1/2` only when `m=3`. This observation does not rule out other tail
completions or assert a finite construction even for `m=3`.

In particular, taking `m=7` gives eight noncoalescing entire row factors,
but with group sizes `7+1`, not `4+4`, and no finite orthogonality.
Consequently the displayed entire residual identities alone do not
bound the total number of rows by four when group sizes are unrestricted.
Whether the balanced eight-row structure adds an obstruction remains open.

The [exact checker](check_critical_moment_entire_one_group_feasibility.py)
verifies the symbolic residual identities, the opposite-root determinant,
the normalized derivatives, and reflection products at rational algebraic
fixtures. The real-zero and distinct-absolute-value assertions are proved
analytically above, not certified by numerical root sampling.
