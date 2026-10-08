# Gap-free ordered binary allocation diagnostic

This is a purely combinatorial check. It does not infer real phases, Gaussian
points, or a circle configuration from binary row order.

For a binary column (w_0,dots,w_{M-1}), a five-row endpoint/interior gap is
present exactly when some ordered five-subsequence has pattern `01110` or
`10001`. The diagnostic enumerates all allowed columns and all multisets of
four columns, requiring every pair of ordered rows to have Hamming distance at
least (4/2=2), the unit-weight near-obtuse condition.

The exact output for (5\le M\le11) was:

```text
M=5 words=30 k4_example=1
M=6 words=50 k4_example=1
M=7 words=70 k4_example=0
M=8 words=90 k4_example=0
M=9 words=110 k4_example=0
M=10 words=130 k4_example=0
M=11 words=150 k4_example=0
```

An exact (M=6) example is given by the four columns

```text
110000  101100  101010  110110
```

All four columns have no forbidden five-subsequence. The six row words are

```text
1111, 1001, 0110, 0101, 0011, 0000,
```

with pairwise Hamming-distance matrix having diagonal zero, off-diagonal
entries (2) except the two opposite pairs with distance (4). Thus all
endpoint/interior intervals can be gap-free while the weighted row metric is
near-obtuse at equality.

The (K=4) enumeration is only a small diagnostic. The general elementary
transition count gives each gap-free column at most five transitions, hence
unit-weight adjacent separation implies (M-1\le10), or (M\le11). This
finite experiment does not improve that bound and does not address whether
opposite-twist descents preserve an actual geometric phase ordering.
