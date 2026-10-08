# No punctured orthogonal eight-row sign template at length sixteen

Let `T` be any 8-by-16 sign matrix with orthogonal rows. Delete one
constant column and one column balanced on its eight rows, leaving an
8-by-14 matrix `S`. There are no nonzero real coefficients `a_j` with
pairwise distinct absolute values satisfying

```text
S(a) = alpha_1 1,
S(a^3) = alpha_2 1,
S(a^5) = alpha_3 1.                                      (1)
```

This covers **all partial orthogonal sign matrices** of this size; no
Hadamard completion is assumed. The result uses the exact critical and
successor [root sign orders](critical_odd_moment_root_order.md), together
with the common [spectral coefficient ordering](critical_moment_spectral_orientation.md).
It is a finite sign-pattern exclusion, not an endpoint-radius theorem.

## Exact pair requirements

The deleted balanced column partitions rows into two four-row groups.
Every row pair disagrees in eight columns of `T`. After the deletions,
cross-group pairs disagree in seven retained columns, while pairs within
one group disagree in eight.

For a cross pair, the signed sequence
`((S_ij-S_kj)/2) sign(a_j)` on its seven supported columns, read in
increasing `|a_j|`, must be `++--++-` or its global sign negation. For an eight-column
pair the successor rule permits type A `+--++--+` or type B
`+++--++-`, again up to global sign negation. The total sign is zero in type A and
`+2` or `-2` in type B.

The seven-column critical patterns make the four row sums in at least one
deleted-column group equal. Call that group `A`; the other is `B`.
The rows of `B` fall into two classes according to whether their row sum
is two above or below the common `A` sum. Let the sizes be `q` and `4-q`.
Within each class, and within `A`, pairs have successor type A. Pairs
between the two `B` classes have type B. All sixteen cross pairs are
critical.

The spectral coefficient lemma makes the first pattern orientations
depend on one strict row ordering. Here `s=3`; label the rows
in descending spectral-coefficient order. Then for every pair `i<j`,
the first sign is `+` for a critical or successor-B pair and `-` for
a successor-A pair. The three classes are consecutive in this order;
row relabeling leaves only `q=0,1,2,3,4` to check.

## Exhaustive exact search

The [checker](check_partial_orthogonal16_critical_patterns.cpp) appends
retained signed columns in increasing order of `|a_j|`. It tries every
eight-bit signed column. A constant column is included as a no-op, and
repeated columns are allowed. For each of the 28 row pairs, it records
the exact number of separating columns already appended. A separating
column must give the next prescribed sign of that pair's critical,
successor-A, or successor-B pattern. Counts may never exceed seven or
eight, and must be able to reach those targets with the remaining
columns.

Two prefixes at the same depth with the same 28 counts have identical
possible suffixes, since every pair's type and orientation is already
fixed. Memoizing these states is therefore exact. At depth fourteen,
the checker accepts only if all sixteen cross counts equal seven and
all twelve within-group counts equal eight. It exhausts the following
dead states without finding an assignment:

```text
q (B class size)          0      1      2      3      4
exact dead states     59684  87000  77636  87000  59684
```

These individual distance requirements matter. A coarser automaton
that remembers counts only modulo four can reach a terminal pattern
state and the correct aggregate pair-defect value at length sixteen,
but those paths need not give every pair its required distance. The
exact search above closes that gap.

Run the checker with

```text
c++ -O3 -std=c++17 docs/check_partial_orthogonal16_critical_patterns.cpp -o /tmp/check_partial16
/tmp/check_partial16
```
