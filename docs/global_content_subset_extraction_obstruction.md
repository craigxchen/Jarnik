# No bounded row extraction for the global quartet content

For a row set `T` of at least four distinct Gaussian integers, write

```text
G(T) = gcd_Q gamma_Q,
```

where `Q` ranges over four-element subsets of `T`.  A subset has fewer
quartets, so `G(T)` divides `G(S)` whenever `S` is contained in `T`.  The
question is whether a fixed number of rows can always be retained with
`G(S)=G(T)` for primitive equal-norm circle rows.

The answer is no without using the endpoint angular hypothesis.  For every
`m>=4`, there is a primitive `m`-row equal-norm circle tuple for which every
row is needed.  Therefore a bound such as six or eight cannot follow from
the circle equation, primitivity, or the source valuation profiles alone.

## Singleton source blocks

Fix `m` and let

```text
M = lcm { j^2-k^2 : 1 <= k < j <= m },
κ_j = 2 j M q + ι       (1 <= j <= m),
n_j = Norm(κ_j).
```

The positive integers `n_j` are pairwise coprime.  Indeed, a prime dividing
both `n_j` and `n_k` cannot divide `2Mq`; subtracting the two equations then
puts it into `j^2-k^2`, hence into `M`, a contradiction.  Every odd prime
`p<=m` divides `M`: use `(p-1)^2-1=p(p-2)`.  Since `n_j=1 mod p` for
`p|M`, every prime factor of every `n_j` is consequently greater than `m`.

Define

```text
z_i = conjugate(κ_i) * product(κ_j for j != i).
```

All rows have norm `N=product(n_j)`.  The rows are distinct, since equality
of `z_i` and `z_k` would imply `conjugate(κ_i)κ_k=κ_i conjugate(κ_k)`,
and hence `i=k`.  Their Gaussian gcd is one: a prime factor of `κ_j`
occurs on the `m-1` majority rows, while the singleton row carries the
conjugate factor; the two are coprime because `n_j` is odd.  Thus this is a
primitive equal-norm circle tuple.

## Why each singleton is necessary

Let `p^e || n_j`, and choose a Gaussian prime `rho` over `p` dividing
`κ_j`. On every majority row,
`rho^e` divides the row; on row `j`, `rho` does not divide the row.  The
minimum `rho`-valuation of a majority difference is exactly `e`.  To see
the exactness, write for `i != j`

```text
A_i = conjugate(κ_i) * product(κ_l for l not in {i,j}),
z_i = κ_j A_i.
```

Modulo `rho`, the relation `ι = -2jMq` gives

```text
A_i/A_k = ((i+j)(k-j))/((i-j)(k+j)).
```

The denominators are nonzero because the block norms are pairwise coprime.
For distinct majority labels `i,k`, this ratio is not one: its numerator
minus denominator is `2j(k-i)`, which is nonzero modulo p since `p>m`.
Hence every majority difference has valuation exactly `e`.

The full tuple has `v_rho(G)=e`.  A quartet containing row `j` and three
majority rows has one cross edge of valuation zero and one majority edge of
valuation `e`, attaining this value.  Every quartet has valuation at least
`e`, because it either contains the singleton or consists entirely of
majority rows.  If a retained subset omits row `j`, every retained
difference has valuation at least `e`, so every quartet matching product has
valuation at least `2e`.  Its quartet gcd therefore has strictly larger
`rho`-valuation than `G(T)`.

This applies independently to every `j`; preserving the full Gaussian ideal
forces retention of all `m` singleton rows.  The argument is a direct
set-cover obstruction: the source factor attached to `j` has the unique
requirement “retain `j` and any three other labels.”  Since the singleton
labels are all different, their union is the entire row set.

For comparison, at a two-level source prime whose cut has at least two
rows on each side, the full global matching content has valuation zero
at both conjugate orientations. Preserving that zero contribution after
retaining a subset `T` requires

```text
2 <= |T intersect S| <= |T|-2.
```

If the profile contains every `floor(m/2)`-element cut, a retained set `T`
must have size at least
`floor(m/2)+2+(m mod 2)`: with fewer rows, one can choose a balanced cut
containing zero or one selected label (or its complement containing zero or
one).  This lower bound is attained by any set of that size, since every
balanced cut then meets it in at least two labels and leaves at least two
labels in the complementary side. Thus a complete balanced cut profile has
an exact linear abstract cover cost for preserving its source valuations.
The actual primitive realizations in `quartet_matching_gcd_cut_budget.md`
therefore also require at least this many rows to preserve their full G,
even though `gcd(N,Norm(G))=1`. The abstract threshold is not asserted to
suffice for preserving all accidental prime valuations in those realizations.

## Small exact fixture

For a compact five-row instance, use the blocks

```text
(4+i), (6+i), (10+i), (14+i), (16+i),
```

whose pairwise coprime norms are `17,37,101,197,257`.  The resulting rows
are

```text
(56046,8675), (53954,17475), (51210,24371),
(49746,27235), (49254,28115).
```

They are primitive and have common norm
`3,216,409,741`.  Exact Gaussian arithmetic gives

```text
G = (-544464,725920),   Norm(G) = 823400893696,
```

and every one of the five `4`-row deletions has a strictly larger gcd; the
minimum preserving subset has size five.  The reusable checker
`check_global_content_subset_extraction.py` exhaustively verifies this
fixture family for `m=4,...,8` and checks the source-prime valuations.
For the displayed `2jMq+i` construction with fixed m and q tending to
infinity, the angular span is
`2(arctan(1/(2Mq))-arctan(1/(2mMq)))`, of order `1/q`, whereas the
radius has order `q^m`. Its arc length divided by `sqrt(R)` grows as
`q^(m/2-1)`. Thus it is not an endpoint counterexample.

Consequently, the existing balanced-cut product bound cannot be converted
to a bounded quartet witness by a general gcd extraction lemma.  Any
endpoint theorem must use the actual phase accuracy of the `C sqrt(R)` arc,
or some additional structure beyond these source cuts.
