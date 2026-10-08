# Every nonzero homogeneous binomial obeys the fair cut budget

Let `R_I` be the eight-point phase coordinates in
[the Newton note](gale_torus_newton_cut_orders.md), and use all 127
nontrivial cuts modulo complement with equal weight. If `M,M'` are
products of `d` coordinates and

```text
f=M+c M'
```

is not the zero polynomial, then

```text
sum_S nu_S(f) <= 53d.                              (1)
```

Here `c` is any constant, and the orders are the generic normalized cut
orders. Thus the quadratic census extends to all homogeneous binomials.
This does not bound sums of three or more monomials and proves no uniform
arc theorem.

## 1. Factor exponents and possible cancellations

For a coordinate product, let `c_i` count its index triples containing
label `i`, and let `e_ij` count its complementary Vandermonde factors
containing edge `ij`. Then

```text
M=product_i z_i^(2c_i) product_(i<j)(z_j-z_i)^e_ij,
sum_(j different from i)e_ij=4(d-c_i).              (2)
```

As proved in [the independent binomial audit](gale_binomial_independent_audit.md),
proportional initial forms at even one cut force the same `c_i` for the
two products: inside-variable minimum degrees are `2c_i`, and
outside-variable maximum degrees are `4d-2c_i`.

If the row weights differ, no cut admits a leading cancellation.
The total order of each individual product is `53d`, so summing their
cutwise minimum orders proves (1) immediately.

Suppose therefore that the row weights agree, and set

```text
b_ij=e_ij-e'_ij,       b_ji=b_ij,       b_ii=0.
```

Equation (2) gives the exact signed row sums

```text
sum_j b_ij=0.                                       (3)
```

If all `b_ij` vanish, unique factorization makes `M=M'`; their nonzero
linear combination is a constant multiple of one product and satisfies
(1). This handles coordinate-monomial identities such as the
degree-four Pasch trade. Henceforth `b` is nonzero.

Let `G` be the graph of edges with `b_ij!=0`. At a cut, the two leading
forms are proportional exactly when every edge of `G` crosses the cut.
Indeed retained within-side difference exponents must agree; conversely
the all-cross condition and (3) make all outside monomial powers agree
as well. Call the set of these cuts `C`.

At a cut in `C`, a suitable sign may cancel the leading forms. The
first-order coefficient after cancellation is a nonzero linear
combination of the distinct Laurent monomials `x_i/y_j` with
coefficients `b_ij`. Therefore the extra order is exactly one. The total
extra order for any fixed scalar `c` is at most `|C|`.

## 2. Two disjoint injections pay for every possible cancellation

Write the difference of the product orders as

```text
t(S)=nu_S(M)-nu_S(M')=sum_(i<j, i,j in S)b_ij.        (4)
```

The variable powers and common primitive baseline cancel because the
row weights agree. The zero row sums imply `t(S)=t(S^c)`, so this is a
function on cuts modulo complement.

If `C` is empty there is no gain to account for. Otherwise `G` is
bipartite. Every nonisolated vertex has at least two neighbors: a
single nonzero incident coefficient would contradict (3). Choose two
adjacent nonzero edges `ij` and `jk` in one component, with distinct
vertices `i,j,k`.

For `S in C`, define two cuts by flipping the endpoints of these edges:

```text
phi_1(S)=S symmetric_difference {i,j},
phi_2(S)=S symmetric_difference {j,k}.               (5)
```

Each edge crosses `S`, so a flip exchanges one included and one excluded
vertex. In particular the cut remains nontrivial and has the same
cardinality. The operations commute with taking complements and are
individually injective on the cut classes.

For the first map, suppose `i in S` and `j notin S`. Originally all
supported edges cross. The only nonzero internal coefficients in the
new cut are those joining the newly included vertex `j` to
`S minus {i}`. By (3), their sum is `-b_ij`. Reversing the roles of
`i,j` gives the same conclusion. Therefore

```text
t(phi_1(S))=-b_ij,       t(phi_2(S))=-b_jk.          (6)
```

The two image sets in (5) are disjoint, even modulo complement. Two
good cuts differ by unions of whole connected components of `G`
(isolated vertices count as components). Equality of an image of
`phi_1` and an image of `phi_2` would make the underlying good cuts
differ by `{i,k}`, or its complement. Neither can be a union of whole
components: `i,k` lie in the component containing `j`, while `j` is not
in `{i,k}`. Taking a complement does not change this failure.

Since the nonzero integers `b_ij,b_jk` have absolute value at least one,
the two disjoint image families give

```text
sum_S |t(S)| >= (|b_ij|+|b_jk|)|C| >= 2|C|.         (7)
```

This is the arithmetic step: integer edge multiplicities supply the
unit lower bound. No claim about arbitrary real edge coefficients is
needed.

## 3. Completing the order bound

Since both products have total order `53d`,

```text
sum_S min(nu_S(M),nu_S(M'))
  =53d-(1/2)sum_S |t(S)|.
```

The possible cancellation gain is at most `|C|`, and (7) pays for it.
This proves (1). In the Newton-support notation, every nonzero
homogeneous binomial pullback consequently has `sum_S h_S>=320d`.

The proof uses the full equally weighted cut collection. It does not
assert the same comparison for arbitrary cut weights. As with the
underlying generic-order formulas, coefficient heights, denominators,
and exceptional-prime specializations require their arithmetic error
budget when applying this result to actual integer configurations.

[check_gale_binomial_cut_flip_bound.py](check_gale_binomial_cut_flip_bound.py)
checks the graph lemma on all 19682 nonzero signed four-by-four
bipartite matrices obtained by choosing the first three-by-three block
from `{-1,0,1}` and forcing every row and column sum to zero. It verifies
the complementary-cut identity, both injections, their disjointness,
the exact values (6), and (7). The proof above applies to every integer
matrix satisfying (3), not just this finite collection.
