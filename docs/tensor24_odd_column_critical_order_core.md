# A ten-certificate core for the tensor sign-order obstruction

The fixed tensor sign matrix has sixteen distance-11 certificates in
`tensor24_odd_column_critical_order_obstruction.md`. I independently
reimplemented the prefix search with a ternary orientation tuple
`(0, +1, -1)` rather than the original orientation bitmask.

The full 16-certificate breadth-first enumeration reproduces 3,029 reachable
states and no complete 22-column assignment. A greedy deletion pass removes
the first six certificates while preserving infeasibility. The remaining ten
row-pair certificates are

```
(1,4), (2,5), (2,6), (2,7), (3,5),
(3,6), (3,7), (4,5), (4,6), (4,7).
```

The independent 10-certificate search has the following numbers of states at
depth 0 through 21:

```
1, 43, 302, 822, 1672, 3202, 5390, 8582, 13098, 18122,
21838, 23620, 24074, 21628, 15884, 9472, 5466, 3186,
1124, 200, 28, 0.
```

It reaches 177,754 distinct prefix states and dies at depth 21, so the ten
certificates already form a checkable infeasible core. This is a certificate
selection result for the fixed matrix only; minimality below ten certificates
was not claimed.

The greedy deletion sequence was

```
removed pair   remaining certificates   states     last nonempty depth
(0,1)                    15              3,551             13
(0,5)                    14              6,057             16
(0,6)                    13             11,277             16
(0,7)                    12             24,922             18
(1,2)                    11             57,466             20
(1,3)                    10            177,754             20
```

Each reduced search remained infeasible in the exploratory deletion run.
The standalone script reproduces the full and final ten-pair searches;
it does not rerun the five intermediate deletion counts.

The [independent checker](check_tensor24_odd_column_critical_order_core.py)
also verifies the original 16-certificate state count and depth profile before
checking the reduced core.
