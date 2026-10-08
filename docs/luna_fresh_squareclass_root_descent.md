# Exact squareclass root descent and the moving-class gap

This note isolates the part of central truncation that descends the endpoint
scale when the retained Gaussian exponents are even. It is a geometric lemma
for an exact Gaussian squareclass, followed by the precise obstruction when
the squareclass varies.

## Exact root lemma

Let `z_1,...,z_m` be distinct Gaussian integers of common modulus `R>0` in
an arc of length at most `C sqrt(R)`. Suppose that for one fixed nonzero
Gaussian integer `b` and one fixed unit `epsilon`,

```text
z_j = epsilon b w_j^2,       w_j in Z[i].                 (1)
```

Then

```text
m <= 1 + C/(2 sqrt(2 |b|)).                             (2)
```

In particular, an exact squareclass has a radius-independent count for every
fixed `C`. To apply this with `b=1` to an endpoint allocation, all retained
Gaussian prime exponents must be even, up to the fixed unit.

### Proof

Write `z_j=R exp(i theta_j)` using arguments in the interval defining the
arc, of width `Delta<=C/sqrt(R)`. Choose for each `w_j` the square-root
argument satisfying

```text
arg(w_j) = (theta_j-arg(epsilon b))/2
```

inside one interval of width `Delta/2`; changing a root by its sign gives
the other choice. Every `w_j` has common modulus

```text
rho = sqrt(R/|b|),
```

and the physical length of their root arc is at most

```text
rho Delta/2 <= C/(2 sqrt(|b|)).                         (3)
```

The `w_j` are distinct after this sign choice: if two chosen roots agree,
their squares, and hence the corresponding `z_j`, agree. Distinct Gaussian
integers have squared distance a positive integer. For equal-norm Gaussian
integers it is in fact a positive even integer, because
`|w_j-w_k|^2=2(rho^2-Re(w_j conjugate(w_k)))`. It is therefore at least two.
Packing
successive points along the root arc therefore gives (2).

The same argument applies separately to the two unit classes modulo unit
squares: `{1,-1}` and `{i,-i}`. If `C<2 sqrt(2)`, the bound is less than
two for every nonzero integral representative `b`, so distinct points on
such an arc have distinct Gaussian squareclasses, including the unit class.

## Why this does not yet handle central truncation

When the retained exponents are even, central truncation gives

```text
z_j = B_j w_j^2,       |B_j|=D,                         (4)
```

For example, this parity condition holds when the original squared radius
is an integer square: each split-prime exponent `e_p` is even and truncation
replaces it by `e_p-2r_p`. In (4), `D` has logarithm small relative to the
original conductor weight,
but the correcting factors `B_j` may vary with `j`. The exact lemma applies
after partitioning rows by the literal value of `B_j` (and its unit
squareclass), giving

```text
m <= (1 + C/(2 sqrt(2))) * (# squareclasses).             (5)
```

This is a genuine uniform bound per class. The number of possible Gaussian
integers of modulus `D` can still grow with `D`, and `log D=o(log R)` does
not make that number bounded. Hence no conclusion about the total number of
rows follows from (4) without a new argument controlling the number of
classes or forcing repeated correction factors.

The point of the lemma is therefore diagnostic as well as positive: the root
map resolves the even-exponent part of pure endpoint allocations, while the
moving corrections and odd retained exponents are the remaining arithmetic
gap. Any argument that silently replaces `B_j` by one, or silently takes
roots when a retained exponent is odd, has dropped this gap.

## A finite test of the scale

For `b=1`, the root arc has length at most `C/2`, independent of `R`. The
bound (2) gives at most `1+C/(2 sqrt(2))` points by one-dimensional packing.
The spacing argument uses only that two equal-norm Gaussian integers have
positive even squared distance; it does not require an odd split-prime
circle or a radius threshold.

The unresolved statement is a class-control theorem: show that an endpoint
cluster with factors (4) can occupy only `O_C(1)` squareclasses, or derive a
second descent from the differences between the `B_j`. Small logarithmic
correction height alone is insufficient for either assertion.
