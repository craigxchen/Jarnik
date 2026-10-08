# Independent audit of the quadratic binomial cut census

The factor-based census in
[check_gale_binomial_all_weights.py](check_gale_binomial_all_weights.py)
has an exact generic-order interpretation. The independent checker
below reproduces its counts and also checks five selected binomials by
literal polynomial expansion at every one of the 127 cuts.

For a product of `d` phase coordinates, let `c_i` be the number of its
index triples containing `i`, and let `e_ij` count the factors whose
complementary five-subset contains the edge `ij`. The product is exactly

```text
M(z)=product_i z_i^(2c_i) product_(i<j)(z_j-z_i)^e_ij.
```

In particular its total degree in variable `z_i` is
`2c_i+4(d-c_i)=4d-2c_i`.

## Proportional cut initial forms determine the row weights

At a cut put `z_i=epsilon x_i` inside and `z_j=y_j` outside. After
removing the lowest power of `epsilon`, each cross edge becomes a signed
outside variable, and edges lying within one side remain differences.
The initial form is a product of variables and irreducible differences.
These factors are pairwise nonassociate over the rational or complex
coefficient field. Its factor exponents therefore determine it up to a
constant; the constants occurring here are exactly `+1` and `-1`.

For an inside label, the minimum exponent of its variable in this
initial polynomial is `2c_i`. None of the remaining difference factors
is divisible by that variable. For an outside label, the maximum
exponent is the outside monomial power plus all incident remaining edge
multiplicities, namely `4d-2c_i`. Products preserve these extreme
degrees. Thus two proportional initial forms determine the same `c_i`
for every label. This argument works in every degree, not only two.

## Exactly one additional order for distinct polynomial products

After factoring out its signed initial form, the remaining contribution
to a product is

```text
product_(i inside, j outside)(1-epsilon x_i/y_j)^e_ij.
```

If two initial forms are proportional and their signed leading terms
are cancelled, their first-order difference has coefficient

```text
-sum_(i inside, j outside)(e_ij-e'_ij) x_i/y_j.
```

The ratios `x_i/y_j` are distinct Laurent monomials, hence linearly
independent. The coefficient is nonzero whenever one cross-edge
multiplicity differs. If none differs, equality of initial forms
already gives the same row powers and internal-edge multiplicities,
and the two original polynomials are identical. Therefore cancelling
proportional initial forms of two distinct polynomial products gives
exactly one extra generic order.

The qualification about distinct polynomial products matters in higher
degree. For example the degree-four coordinate products indexed by

```text
012,034,135,245       and       013,024,125,345
```

are identical polynomials. This is the six-label Pasch trade: their
vertex and pair-incidence counts agree, so all variable and edge factor
exponents agree. Their difference is identically zero, not a nonzero
binomial with one extra order. In degree two, the enumeration verifies
that all 1596 coordinate products have distinct factor exponents.

## Counts and comparison with full polynomial supports

Every quadratic coordinate product has total normalized order 106.
For two such products the sum of their cutwise minimum orders is
`106-(1/2)sum_S|nu_S(M)-nu_S(M')|`. Only coefficients `+1` or `-1`
can cancel any proportional initial forms. The extra-order formula
above therefore completes the census without sampling.

There are 2100 product pairs with some proportional initial form and
18480 pair-cut coincidences. Their two signed binomials have order-sum
histogram

| Total order | Number |
|---:|---:|
| 82 | 1260 |
| 86 | 1260 |
| 90 | 840 |
| 106 | 840 |

Pairs with no proportional initial forms have no extra orders and
cannot exceed 106. Other nonzero coefficients cannot improve the
generic cut order. Thus no nonzero quadratic binomial gives a strict
gain over the degree-two fair height budget. This concerns generic
characteristic-zero orders; coefficient denominators and exceptional
primes still require their ordinary arithmetic accounting.

[check_gale_binomial_all_weights_independent.py](check_gale_binomial_all_weights_independent.py)
independently reconstructs factors from coordinate triples and checks
the full census. For one signed binomial from each histogram row, it
expands both coordinate products into their actual monomials, combines
coefficients, and obtains all 127 cut orders directly from the resulting
Newton support. It also checks a coefficient-two binomial with no
coincident initial form, and verifies the higher-degree Pasch example.
