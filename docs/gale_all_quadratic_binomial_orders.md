# Every quadratic binomial stays within the fair cut budget

Let `R_I=(product_(i in I)z_i^2)V_(I^c)` be the eight-point polynomial
coordinates in [the Newton cut formulation](gale_torus_newton_cut_orders.md).
For every nonzero quadratic binomial in these coordinates, with arbitrary
fixed rational coefficients, the sum of its 127 generic primitive
clean-cut orders is at most 106. This includes mixtures of different
row weights. It does not bound arbitrary sums of three or more monomials.

## 1. Factorization determines all possible leading cancellations

For a product of `d` coordinates, let `c_i` be the number of selected
triples containing label `i`, and let `e_ij` count how many complementary
five-subsets contain the pair `i,j`. Its exact factorization is

```text
M(z)=product_i z_i^(2c_i) product_(i<j)(z_j-z_i)^(e_ij).
```

In a cut `z_i=epsilon x_i` on `S`, `z_j=y_j` off `S`, factors within
one side remain differences. A crossing difference becomes, to leading
order, plus or minus its outside variable. Thus its initial form is a
sign times a monomial times products of distinct irreducible differences
within each side. Unique factorization determines exactly when two such
initial forms are proportional; their proportionality constant is `+1`
or `-1`.

The initial form also determines the row weight. For an inside variable
its minimum exponent is exactly `2c_i`. For an outside variable its
maximum exponent is exactly `4d-2c_i`, since its total incident difference
degree is `4(d-c_i)`. Proportional initial forms therefore require all
the row weights `c_i` to agree. Different row weights never have a
two-term leading cancellation on any cut.

If two distinct polynomial products have equal initial forms after
adjusting the sign, factor this common initial form from their expansions.
The coefficient of the next power of `epsilon` is the difference of
their crossing-edge multiplicities, multiplying the distinct rational
monomials `x_i/y_j` for `i in S`, `j` outside. It is nonzero unless all
crossing multiplicities agree. In that exceptional case their retained
factors and monomial powers also agree, so the two full polynomial
products were identical. Hence distinct products gain exactly one order
when their leading forms cancel. This argument applies in any degree;
coordinate monomials which become identical polynomial products must
always be removed first.

## 2. The complete quadratic count

There are 56 coordinates and 1,596 products of two coordinates with
repetition. Their full irreducible-factor signatures are all distinct.
Each product has total primitive cut order 106. For two such products
with valuation lists `u_S` and `v_S`, put

```text
loss=(1/2) sum_S |u_S-v_S|.
```

Their generic minimum-order sum is `106-loss`. If the coefficient ratio
is not `+1` or `-1`, that is the binomial's exact order sum. For either
sign, add the number of cuts on which the initial forms cancel, exactly
one per cut by Section 1.

Only 2,100 pairs have a possible initial-form cancellation. All have one
common row weight. Up to relabeling, their exact data are:

| Weight type | Number of pairs | Loss | Gains for the two signs | Resulting sums |
|---|---:|---:|---:|---:|
| one repeated label, four single labels | 840 | 16 | 0 and 16 | 90 and 106 |
| six single labels | 1,260 | 24 | 0 and 4 | 82 and 86 |

Every omitted pair has no gain and hence total at most 106. A single
monomial has total exactly 106. This proves the claimed bound for all
quadratic binomials, including arbitrary coefficients and all row weights.

The checker
[check_gale_binomial_all_weights.py](check_gale_binomial_all_weights.py)
enumerates every product and every cut using integer factor signatures,
checks the nonzero first-jet criterion for all 18,480 leading-form
coincidences, and verifies the table. It uses neither finite-field
sampling nor numerical evaluation. The finite enumeration is complete
for this specified binomial class. The subsequent
[support-multiplicity theorem](polarized_coefficient_support_multiplicity.md)
now settles arbitrary cross-weight cancellations and every degree, with
the same sharp budget `sum nu_S<=53d`.
