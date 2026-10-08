# A coefficient floor from Weyl-related binary residuals

This note proves a higher-rank extension of the actual-cut coefficient
floor under the standard highest-root description of the conformal-block
(CB) annihilator and the extracted Gaussian core hypotheses of the
actual-cut theorem: oriented core factors have pairwise disjoint prime
support and `log Norm(H_T) >= (1-eta)w`. The projected-block isolation
argument is proved in Section 7 of
[`coherent_cut_evaluation_differentials.md`](coherent_cut_evaluation_differentials.md).
The finite checker verifies the combinatorics and coefficient
bookkeeping, not the arithmetic core hypotheses.

Let `n>=2`, `d>=2`, and `m=nd`. Let `Q` be an integral multilinear
`SL_n` invariant in `m` vector rows. Its `GL_n` relative weight is
`det^d`. Thus every monomial assigns exactly `d` labels to each of the
`n` coordinate families. Write its coefficient norm `C_Q` in the
standard row-coordinate monomial basis.

## The CB condition makes every binary residual vanish

For a coordinate pair `(a,b)`, define

```text
D_(z,ab)=sum_i z_i v_i^a partial_(v_i^b).
```

The CB-annihilator condition is `D_(z,ab)^d Q=0`. It holds for every
coordinate pair by Weyl conjugation of the highest-root condition:
`Q` transforms by the determinant character, so conjugating the
operator preserves its annihilation. Multilinearity and degree `d` in
coordinate family `b` give

```text
D_(z,ab)^d Q = d! Q(..., v_i^b=z_i*v_i^a, ...).
```

On the dense chart `v_i^a != 0`, factor `v_i^a` from each row. This
says `Q` is zero with pair coordinates `(1,z_i)` and all other
coordinates arbitrary. Polynomial continuation removes the chart
restriction. Therefore every
coefficient residual in

```text
R_z Q = Q((1,z_i,x_i^2,...,x_i^(n-1))_i)
```

has value zero at the actual binary rows. For homogeneous rows
`(P_i,Y_i)` with `z_i=Y_i/P_i`, an outside residual `r_O` satisfies

```text
r_O(P,Y) = (product_(i in O) P_i) * r_O(1,z).
```

The actual core-profile rows have `P_i != 0`, as required for this
homogeneous normalization. This conclusion needs the full
CB-annihilator condition. A single numerical equality
`Q((P_i^(n-1),...,Y_i^(n-1))_i)=0` does not imply that the residual
values vanish; the six-row sign countermodel in
[`terminal_ordered_circuit_sign_countermodel.md`](terminal_ordered_circuit_sign_countermodel.md)
demonstrates that distinction.

## Coefficient extraction preserves height

Under the row Veronese substitution

```text
(P_i,Y_i) |-> (P_i^(n-1), P_i^(n-2)Y_i, ..., Y_i^(n-1)),
```

a source monomial assigning coordinate `j_i` to row `i` becomes
`product_i P_i^(n-1-j_i) Y_i^j_i`. Its exponent vector `(j_i)_i`
recovers the assignment, so distinct source monomials never collide.
The pullback has exactly the same coefficient `l1` norm `C_Q`.

Fix disjoint label sets `S_2,...,S_(n-1)`, each of size `d`, for the
free coordinate families. The coefficient of

```text
product_(j=2)^(n-1) product_(i in S_j) x_i^j
```

is a binary residual on

```text
O = [m] minus union_(j=2)^(n-1) S_j,     |O|=2d.
```

It is homogeneous of degree one in every outside row and is `SL_2`
invariant: the block-diagonal `SL_2` on the chosen coordinate pair
fixes all free-coordinate basis rows. Distinct free-set assignments
partition the source monomials, so the sum of the residual coefficient
`l1` norms is exactly `C_Q`. In particular, a nonzero `Q` has a
nonzero residual polynomial. Since every such residual vanishes at
the actual rows, the projected-block isolation argument applies to at
least one nonzero residual.

For one fixed free-set assignment, the cut has inside size
`s=(n-2)d` and outside size `2d`. The actual-cut conductor argument
therefore gives

```text
log C_Q >= 2^((n-2)d+1) * (1-eta) * w
           - 2 * sum_(i in O) log|K_i|.                  (1)
```

The first `2^s` is the number of full-core blocks with a fixed outside
intersection; the second factor of two comes from the equal-magnitude
complementary binary coefficient of an `SL_2` invariant. The core
profile and row-correction terms are those in the actual-cut theorem;
the coefficient norm here is the residual monomial `l1` norm, so no
matching-basis conversion is used. Formula (1) is not a stand-alone
claim for arbitrary integer rows. The modular isolation argument in
Section 7 cited above uses only an integral binary residual and the
`2^s` projected blocks, so it extends to odd `m=nd`.

## One source coefficient across all coordinate pairs

There is a stronger floor under the same hypotheses. Fix one nonzero
source coefficient `c`, indexed by a partition of the labels into
coordinate blocks `B_1,...,B_n`, each of size `d`. A global coordinate
permutation preserves its absolute value. For each ordered pair `a!=b`,
use the binary residual on `O=B_a union B_b`, fixing the other
coordinate families at their blocks. Orient the pair so the
residual coefficient of `Y_(B_a) P_(B_b)` is `+/- c`. Projected-block
isolation gives the associated core factors for

```text
T contains B_a,       T intersect B_b = empty.
```

The union over ordered pairs consists exactly of subsets `T` that
contain at least one whole coordinate block and omit at least one
whole coordinate block. Its cardinality is

```text
U(n,d)=2^(nd)-2*(2^d-1)^n+(2^d-2)^n.                  (2)
```

The subtracted terms count subsets with no full block and with no
empty block; their intersection has `(2^d-2)^n` choices. Assign
each union block `T` to one witness ordered pair. Pairwise disjoint
core prime supports let the assigned divisors accumulate in the same
integer coefficient `c`. For pair `(a,b)`, only the `P` rows in `B_b`
contribute the row-correction loss, bounded by
`2 sum_(i in B_b) log|K_i|`. Summing over all ordered pairs bounds
the total loss by

```text
2*(n-1)*sum_i log|K_i|.
```

Thus

```text
log C_Q >= log|c| >= U(n,d)*(1-eta)*w
           - 2*(n-1)*sum_i log|K_i|.                   (3)
```

This uses one fixed coefficient throughout; it does not assume
different residuals have independent coefficients. The claim uses the
core-factor disjointness and row-correction hypotheses in Section 7
cited above. The exact count and correction bookkeeping are checked in
[`check_higher_rank_weyl_union_floor.py`](check_higher_rank_weyl_union_floor.py).

For context, the generalized CB annihilator and finite fusion-height
comparison are developed in
[`higher_rank_terminal_fusion_grid.md`](higher_rank_terminal_fusion_grid.md).
