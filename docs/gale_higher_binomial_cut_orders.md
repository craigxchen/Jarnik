# Cubic Gale binomials have no weighted cut-order gain

The exact phase pullback and all-cut filtration from
[`eight_point_weighted_cut_sections.md`](eight_point_weighted_cut_sections.md)
make cubic binomials a finite calculation.  The result is sharp:

```text
Every nonzero degree-three binomial in coordinate products has
sum of generic orders over the 127 clean cuts at most 159.          (1)
```

Since `159=3*53`, no cubic binomial exceeds the fair-profile
archimedean budget.  Some binomials attain equality, but none gives a
strict gain.

## 1. Factor data and cut initials

Recall the actual circle-moment pullback

```text
R_I=(product_(i in I) z_i^2)
    product_(j<k, j,k not in I)(z_k-z_j),       |I|=3.       (2)
```

For a product `M=product_(a=1)^d R_(I_a)`, put

```text
c_i=# {a:i in I_a},
e_ij=# {a:i,j not in I_a}.                              (3)
```

Unique factorization gives

```text
M=z_1^(2c_1)...z_8^(2c_8)
  product_(i<j)(z_j-z_i)^(e_ij).                       (4)
```

Thus `(c_i,e_ij)` is the full factor signature of a coordinate product.
For degree three, direct enumeration verifies that the signatures of all
`binomial(58,3)=30856` unordered products are distinct.  Hence two distinct
cubic coordinate products are distinct polynomial functions.  This
qualification matters in higher degree, where factor trades may first
produce identical products.

At a cut `S`, write `z_i=epsilon*x_i` on `S` and `z_j=y_j` off `S`.
Each cross factor has the form

```text
y_j-epsilon*x_i = y_j(1-epsilon*x_i/y_j),
epsilon*x_j-y_i = -y_i(1-epsilon*x_j/y_i).            (5)
```

The cut initial of (4) is therefore recorded exactly by:

- its orientation sign;
- the monomial powers of the eight cut-unit variables;
- the powers `e_ij` on all same-side differences.

These are sparse integer data.  No polynomial expansion is required.

If two product initials are proportional, their incidence vectors `c`
are equal.  Indeed, for `i in S` the least exponent of `x_i` in the initial
is `2c_i`.  For `i` outside `S`, the greatest exponent of `y_i` is
`4d-2c_i`: each factor containing `i` contributes two, while each factor
excluding it contributes four across its Vandermonde factors.  Proportional
initials have the same extrema and hence recover the same `c_i` on both
sides of the cut.

Suppose now that two distinct products of the same row weight have
proportional initials.  Their same-side edge powers agree.  Their full
factor signatures are distinct, so at least one cross-edge power differs.
After removing the common initial, (5) shows that the difference of their
first logarithmic jets is

```text
-sum_(i in S, j outside S) (e_ij(M)-e_ij(N)) x_i/y_j,   (6)
```

with the indices oriented only to determine the already recorded sign.
The monomials `x_i/y_j` are independent, so (6) is nonzero.  Consequently
the unique choice of coefficient sign which cancels proportional initials
raises the cut order by exactly one, never by two or more.  Because all
initial coefficients are signs, this cancellation requires coefficient
ratio `+1` or `-1`.  Arbitrary other binomial coefficients give no bonus.

## 2. The exact binomial budget

Every degree-three coordinate product has

```text
sum_S nu_S(M)=159.                                    (7)
```

For two such products let their 127 order vectors be `u` and `v`.  Before
initial cancellation, the binomial has total order

```text
sum_S min(u_S,v_S)
 =159-(1/2)sum_S |u_S-v_S|.                           (8)
```

For each coefficient sign, add one for every equal-order cut on which the
factored initials cancel.  Equations (6)-(8) are the complete calculation.
Products of different row weight cannot have proportional initials, so
only products in the same incidence class need to be paired.

Up to relabeling, the possible sorted cubic incidence vectors and the sizes
of their product classes are:

| incidence type | number of labeled weights | products per weight |
|---|---:|---:|
| `33300000` | 56 | 1 |
| `33210000` | 840 | 1 |
| `33111000` | 560 | 1 |
| `32220000` | 280 | 1 |
| `32211000` | 1680 | 3 |
| `32111100` | 840 | 6 |
| `31111110` | 56 | 15 |
| `22221000` | 280 | 6 |
| `22211100` | 560 | 16 |
| `22111110` | 168 | 40 |
| `21111111` | 8 | 105 |

These are 11 types and 5328 labeled row weights.  The largest canonical
calculation has only 105 products.  Scanning both signs for every pair in
one representative of each type gives 12996 signed cases.  The maximum of
(8) plus the exact cancellation count is `159`.  There are 264 attaining
signed cases among the canonical representatives.  Typical equality cases
are a quadratic Plucker exchange multiplied by a third coordinate; after
the exchange the binomial is another coordinate product.

The checker
[`check_gale_higher_binomial_cut_orders.py`](check_gale_higher_binomial_cut_orders.py)
verifies all counts, the injectivity of the degree-three factor signatures,
all monomial order sums, every factored cut initial, and the complete
canonical signed-pair scan.  It uses exact integers and finishes without a
dense coefficient matrix or sampled boundary points.

## 3. Scope

Result (1) rules out a strict weighted-order certificate made from two
degree-three coordinate products.  It does not bound general cubic sums:
three or more terms can cancel leading layers in patterns not represented
by one binomial pair.  It also does not settle degree-four binomials, where
distinct triple multisets must first be quotiented by possible identical
factor products.  Finally, as in degree two, weighted order and absence of
common zeros on the positive actual-circle locus are separate requirements.
