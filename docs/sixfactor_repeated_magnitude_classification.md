# Repeated magnitudes in the six-factor equal-first-moment model

The missing repeated-magnitude case of the signed six-factor model has
only two coefficient rays once five rows are affinely independent.
Both are specializations of the existing six-factor representative, so
its [primitive-content exclusion](sixfactor_moving_content_c_half_exclusion.md)
applies without a new content estimate.

## Exact finite classification

Let `alpha_0,...,alpha_4` be positive and pairwise distinct. Suppose five
rows in

```text
{0,1,2} x {0,1}^4
```

are affinely independent and have a common weighted sum `sum_k alpha_k x_k`.
Up to a positive common coefficient scale and permutation of the four
singleton coordinates, exactly two cases are possible:

| Doubled coefficient first | Five-row fiber sums at the displayed scale |
| --- | --- |
| `(1,2,3,4,5)` | `7` or `9` |
| `(2,1,3,4,5)` | `8` or `9` |

Each listed fiber contains exactly five rows, and the two fibers in each
case are coordinatewise complements using the bounds `(2,1,1,1,1)`.

Here is an exhaustive exact proof, implemented in the
[standard-library checker](check_sixfactor_repeated_magnitude_classification.py).
Two distinct equal-sum rows cannot have distance one or two in the
integer `l1` metric: distance one contradicts positivity; distance two
either contradicts positivity or forces two different coefficients equal.
Enumerate five-cliques of the distance-at-least-three graph on all 48
box points. There are exactly 3360 such cliques.

For each clique subtract a reference row, giving a `4 x 5` integer
matrix `D`. Its signed `4 x 4` minors form a kernel vector. The vector
is nonzero exactly when the rows have affine rank four. In this case
the real kernel is one-dimensional, so it contains a positive,
pairwise-distinct coefficient vector exactly when its entries have one
strict sign and are pairwise distinct. Divide their absolute values by
their greatest common divisor. Exactly 96 cliques pass this test.
Sorting only the last four coordinates of the coefficient vector gives
the two displayed rays, each 48 times. This calculation bounds neither
the input real coefficients nor their denominators: the one-dimensional
kernel determines the entire coefficient ray exactly.

Finally enumerate all 48 row sums for each ray. The only fibers with at
least five rows are those in the table, and each has exactly five rows.
These are the claimed complete alternatives.

## Exact identification with the previous representative

The earlier representative uses

```text
(a_1,...,a_6)=(u,2u,v,u+v,2u+v,3u+v),
S_0={4,6}, S_1={1,3,6}, S_2={1,4,5},
S_3={2,3,5}, S_4={1,2,3,4}.
```

For the first ray take `(u,v)=(2,-5)`. Its signed coefficients are
`(2,4,-5,-3,-1,1)`. Complement the membership bit for each negative
coefficient, then merge the two absolute-one coordinates by adding their
bits. In coordinate order `(1,2,3,4,5)`, the resulting five aggregate
rows are exactly the fiber of sum seven.

For the second ray take `(u,v)=(1,2)`, giving
`(1,2,2,3,4,5)`. Merge the two coefficient-two bits and order coordinates
as `(2,1,3,4,5)`. The five rows are exactly the fiber of sum eight.
The checker verifies these two set equalities literally. The remaining
fiber is the global complement in each case.

For signed products

```text
w_i(T)=(-1)^|S_i| product_(j in S_i)(a_j+iT)
                   product_(j notin S_i)(a_j-iT),
```

replacing a negative coefficient by its absolute value and complementing
its bit preserves the actual product: the factor changes sign, as does
the row prefactor. Interchanging identical factors has no effect.
Globally complementing all six bits conjugates the entire tuple because
the number of factors is even. Thus these identifications preserve the
actual Gaussian tuple, up to conjugation and row permutation. A common
coefficient scale is absorbed into the corresponding `u,v`.

## Consequence on a small arc

Consider six nonzero rational coefficients and a rational parameter `T`,
allowing repeated absolute coefficient values, and five distinct signed
products of the displayed form with exact equal first moments. Clear a
common rational denominator and
divide by their full Gaussian gcd. These primitive points cannot lie on
an arc with normalized length at most `1/2`.

To see why affine independence is available here, first turn every
coefficient positive by the bit complement above, and aggregate equal
coefficients. Write the resulting row vectors as `x_i`. At every chosen
split Gaussian prime, the valuation of the row is an affine function of
these aggregate counts:

```text
v_pi(w_i) = constant + sum_k x_ik
                   [v_pi(alpha_k+iT)-v_pi(alpha_k-iT)].
```

Rational clearing and common Gaussian division only change the constant.
Consequently an affine relation among the `x_i` would give the same
affine relation among actual prime-allocation vectors. The
[allocation rigidity theorem](linear_allocation_affine_rigidity.md)
excludes such a relation for five distinct points at normalized constant
`1/2`, even with arbitrary Gaussian units. The aggregate rows therefore
have affine rank four.

Their common first moment is a nontrivial affine hyperplane equation,
so the number `k` of distinct positive coefficients is at least five.
If `k=6`, the existing distinct-magnitude classification applies. If
`k=5`, exactly one coefficient is doubled and the finite classification
above applies. Both repeated cases identify the tuple with the old
representative, whose content proof requires nonzero coefficients but
**does not require distinct absolute values**. That proof gives every
containing arc the lower bound

```text
normalized length >= 2/sqrt(5 sqrt(2)) > 1/2,
```

a contradiction. This argument covers rational parameters that vary;
it does not assume a fixed limiting coefficient ray in advance. Its
scope is exact equal first moments in this signed six-factor model.
It does not assert a transfer from arbitrary circle tuples or from
approximate moment equalities.

## Fixed constructions of effective width at most nine

Together with the [repeated-factor moment theorem](five_row_nine_factor_moment_exclusion.md),
this closes the remaining small-width case for fixed rational signed
factor constructions. Five fixed distinct aggregate rows of effective
width at most nine cannot lie on arcs of normalized length at most
`1/2` for arbitrarily large positive integral `T`, after primitive
normalization.

Here effective width is the sum of the coordinate ranges after grouping
equal positive magnitudes and removing common polynomial factors. The
remaining specialized Gaussian gcd is bounded by polynomial Bezout,
so the primitive radius has order `T^n_eff`.

If `n_eff<=2`, the number of possible aggregate vectors is at most
`product_j(m_j+1)<=2^n_eff<=4`. For `3<=n_eff<=5`, endpoint scale
forces equal first moments. Distinct aggregate rows then have integer
distance at least three. Threshold embedding gives binary rows of length
`n_eff`. For length at most four, the ten pair distances sum to at least
thirty but at most twenty-four. For length five, append row parity;
the ten new distances are even and at least four, whereas six columns
can contribute at most thirty-six in total. Both cases are impossible.

At width six, endpoint scale again forces exact equal first moments,
and the six-factor theorem above excludes the five points. At widths
seven through nine it forces both first and third moments equal, and
the repeated-factor moment theorem excludes five aggregate rows.
This proves the stated fixed-construction conclusion. Passing to five
rows from a larger fixed family is harmless: any additional common gcd
of that subset only decreases its primitive radius and normalized span.

The large-`T` threshold may depend on the fixed coefficients. This is
not a uniform theorem in the number of factors or in moving coefficients,
and it does not improve the general circle-point growth bound.

Sol independently enumerated the clique graph and the two coefficient
rays, reconstructed the two representative specializations, and audited
the primitive-allocation and content transfer. Root replayed the full
checker and supplied the first ray's signed-parameter identification.

Run the finite classification and exact representative identification with

```sh
python3 docs/check_sixfactor_repeated_magnitude_classification.py
```
