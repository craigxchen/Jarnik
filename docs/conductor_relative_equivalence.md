# Conductor-relative exceptional curves: equivalence audit

## Status

The conductor-relative exceptional-curve statement is a useful formulation, but
its archimedean form is not yet an independent theorem.  Once the endpoint
condition is imposed, bounded exceptional degree and bounded cluster size are
equivalent by an elementary grid argument.

## 1. Endpoint divisor set

Fix a Gaussian integer `G` and a constant `C > 0`.  Define

```text
A_G(C) = { A / conjugate(A) : A | G,
           |A / conjugate(A) - 1| <= C |G|^(-1/2) }.
```

After the usual common-factor normalization, the ratios associated to lattice
points in a `C sqrt R` arc form a subset of `A_G(C)` for the common Gaussian
conductor `G`.

## 2. Bounded degree implies bounded cardinality

Suppose there is a nonzero polynomial `F_G(X,Y)` of total degree at most `D`
which vanishes at every ordered off-diagonal pair from `A_G(C)`.

A nonzero degree-`D` polynomial has at most `D |A_G(C)|` zeros on the grid
`A_G(C) x A_G(C)`.  Since all off-diagonal pairs are zeros,

```text
|A_G(C)| (|A_G(C)| - 1) <= D |A_G(C)|,
```

and hence

```text
|A_G(C)| <= D + 1.
```

This is the grid endgame formalized in
`GaussianChain/UniformExceptionalReduction.lean`.

## 3. Bounded cardinality implies bounded degree

Conversely, if `|A_G(C)| <= B`, then

```text
F_G(X,Y) = product_{u in A_G(C)} (X-u)
```

has degree at most `B` and vanishes on every pair whose first coordinate lies
in `A_G(C)`, in particular on every ordered off-diagonal pair.

Thus the archimedean conductor-relative exceptional statement

```text
all endpoint pairs lie on a curve of uniformly bounded degree
```

is equivalent, up to the harmless `+1`, to the desired uniform bound for the
endpoint divisor set itself.

## 4. Consequence for the Ru--Vojta route

The finite beta calculation remains genuinely useful.  At the very small level

```text
N = 1,  d = 6,  k = 5,
```

the section space has dimension `27`,

```text
beta = 50/27,
6 - 3 beta = 4/9 < 1/2.
```

These identities are formalized in `GaussianChain/FiniteBeta.lean`.

What is not yet proved is the load-bearing Diophantine statement that converts
this strict numerical margin into a single bounded-degree exceptional curve for
all divisor ratios from one conductor.  Aggregating conductor weights controls
the Harder--Narasimhan slope data, but it does not determine the local values
`|L(P)|_v` entering the twisted height.  The prime identities remain visible in
those evaluations.

Therefore the next theorem must be a genuinely conductor-sensitive,
one-scale max-over-bases Subspace Theorem, not merely a compactness theorem for
aggregate slopes.

## 5. Exact remaining target

A sufficient new theorem is:

> For the fixed 27-dimensional section space `H^0(Bl_P P^2, 6H-E)`, there is a
> uniform constant `T` such that, for every Gaussian conductor `G`, all
> conductor-height solutions of the associated strict max-over-adapted-bases
> inequality at the single scale `Q = |G|` lie in at most `T` proper linear
> subspaces.

Pulling those subspaces back gives a curve of degree at most `6T`, and the grid
argument gives the uniform arc bound.

This one-scale theorem is narrower than the false unrestricted uniform-in-`S`
statement, but it is not yet a consequence of the slope aggregation lemmas
alone.
