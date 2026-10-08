# The normalized odd-positive tensor family fails the critical sign order

Every normalized odd-positive eight-row subset of the Paley-12 tensor
Sylvester-2 matrix is excluded by the five odd-moment equations after
deleting the constant column and **any** balanced column. This extends
the [single-template exclusion](tensor24_odd_column_critical_order_obstruction.md)
to a finite family. It does not classify other sign matrices or assert a
uniform lattice-circle bound.

The [exact C++ checker](check_tensor24_odd_column_family_critical_order.cpp)
reconstructs the 24-by-24 matrix and enumerates all `binom(24,8)=735471`
normalized row subsets. Exactly 1,320 have an odd column and at least
fourteen balanced columns. Each has sixteen balanced columns, yielding
21,120 retained 8-by-22 matrices after deleting the constant and one
balanced column.

For each retained matrix, the deleted balanced column partitions the
eight rows into two four-row sets. The sixteen cross-pairs differ in
exactly eleven retained columns. Under the equations

```text
S(a^(2k+1)) = alpha_k 1,       0 <= k <= 4,
```

each cross-pair imposes the [critical root-order pattern](critical_odd_moment_root_order.md)
`++--++--++-` or its global sign negation on its eleven supported columns, ordered
by increasing `|a_j|`. The global order and coefficient signs must serve
all sixteen cross-pairs simultaneously.

The checker first groups matrices by their multisets of unoriented
eight-row column cuts. This allows column permutations and column sign
flips and leaves 9,056 distinct multisets. It then relabels rows so the
deleted balanced cut becomes `{0,1,2,3}` and takes the least key under
all `4! 4! 2 = 1152` permutations within or between the two four-row
sets. This is an exact equivalence for the necessary sign-order test;
row permutations merely rename equations, while column permutations
and sign flips merely relabel or negate the coefficients. Four classes
remain. The checker's representatives all use selected tensor rows
`(0,1,2,5,15,19,20,22)` and delete balanced column `1`, `20`, `3`,
or `17`, respectively.

The exhaustive prefix search of the single-template note is then run
on each representative. It tries every next column and coefficient sign,
recording the selected-column mask and the chosen orientation of each
already touched cross-pair pattern. Those data determine every allowed
continuation, so memoization is exact. All four searches fail:

```text
deleted column       1     20      3     17
reachable states   3029   1479   1335   1279
maximum depth        13     12     12     12
```

No distance-twelve certificate or numerical coefficient search is used.
The exclusion concerns nonzero real coefficients with pairwise distinct
absolute values for this normalized tensor family.

Run the checker with

```text
c++ -O3 -std=c++17 docs/check_tensor24_odd_column_family_critical_order.cpp -o /tmp/check_tensor24_family
/tmp/check_tensor24_family
```
