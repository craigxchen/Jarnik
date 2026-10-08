# First moments on the compressed Euler words

The exact mask words in [the even-`q=2` graph construction](critical_moment_compressed_graph_even_q2.md) admit **strictly increasing positive integer magnitudes satisfying all eight first-moment equations**. This is only a first-moment statement. For the published eleven-column initial path, the later [reciprocal-prefix obstruction](critical_moment_reciprocal_prefix_obstruction.md) proves that every magnitude tuple from the specific recipe below fails the full moments, even after nonnegative order-preserving shifts of later columns. Other positive convex weights remain unresolved. The automaton's polynomial-value order need not be the order produced by these magnitudes.

Write a mask word as columns `S_j in {−1,+1}^8`, with rows numbered `0,...,7`, and let `H_t=sum_{j<=t}S_j`, `H_0=0`. The compressed height is `r_{t,i}=(H_{t,i}−H_{t,0})/2`, with `r_{t,0}=0`. For even `s`, `q=2`, the terminal height is

```
δ=(0,0,0,0,−1,−1,1,1).
```

If `0<a_1<...<a_n`, summation by parts gives

```
sum_j (S_{j,i}−S_{j,0}) a_j
  =2[a_n δ_i−sum_{t=1}^{n−1}(a_{t+1}−a_t)r_{t,i}].      (1)
```

Thus the first moments agree across all rows precisely when `δ` is a strictly positive convex combination of the proper prefix heights `r_0,...,r_{n−1}`. The coefficient at `r_0=0` is `a_1/a_n`, and the others are `(a_{t+1}−a_t)/a_n`.

Every valid compressed state satisfies `|r_i−δ_i|<=1`. The [exact checker](check_critical_moment_first_moment_barycentric.py) gives active paths of length at most three from the terminal state to each of the fourteen axial heights `δ±e_i`, `i=1,...,7`, and active return paths of length at most three. The terminal lies in the graph's closed giant component. The positive active circulation used by the Euler construction puts positive multiplicity on **every** active giant-component edge, so every resulting Euler word visits all fourteen axial heights at proper prefixes. In particular `δ` lies in the relative interior of the convex hull of its proper prefix heights; the affine span is the full seven-dimensional height space.

Here is an exact constructive choice of magnitudes for **any** such word. Put `v_i=sum_{t=0}^{n−1}(δ_i−r_{t,i})`; then `|v_i|<=n`. Give every prefix initial rational weight `1/(8n)`. At one proper prefix with height `δ+e_i`, add `(n+v_i)/(16n)`, and at one with height `δ−e_i`, add `(n−v_i)/(16n)`. All weights are strictly positive. Their sum is `1/8+7/8=1`, and the weighted mean of the heights is exactly `δ`. Multiplying weights by `16n` gives positive integer increments: the base increment is `2`, and the added increments are `n±v_i`. Define `a_j` as the sum of the first `j` increments. Then

```
2<=a_1<a_2<...<a_n=16n,
sum_j S_{j,i}a_j = sum_j S_{j,0}a_j  for every i.      (2)
```

This construction preserves the Euler word's exact pair-distance labels because it changes only the magnitudes, not the signed columns. It does **not** provide a solution of the full odd-moment system or actual circle points.

The checker replays an 87-column finite automaton path containing all fourteen axial heights, computes the integer increments above, and verifies (1)–(2) directly. This short fixture is **not** an exact-distance punctured orthogonal template: its length is not `4s+2`. The positive-circulation argument applies the same algebra to each much longer exact-distance Euler word without explicitly materializing it.
