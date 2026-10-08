# Small-width affine circuits and fixed-factor families

This note converts short-product rigidity into a count in one bounded
exponent-width case. It does not settle arbitrary conductor patterns.

## Integer circuit lemma

Let distinct integer vectors lie in a rectangular box whose integer
coordinate widths have sum at most four. If the vectors are affinely
dependent, some of them have a nonzero integer relation

```text
sum c_j = 0,       sum c_j v_j = 0,
sum_(c_j>0) c_j = sum_(c_j<0) (-c_j) <= 10.          (1)
```

To prove this, take a minimally dependent subset of affine rank `r`.
It has `r+2` elements, with `1<=r<=4`. Select `r` coordinate functions
which, together with the constant function, have full rank on this
subset. The selected positive widths `w_1,...,w_r` still sum to at
most four. Signed maximal minors of the resulting `(r+1)`-by-`(r+2)`
integer matrix give a nonzero integer kernel vector `c`.
The selected rows span every coordinate row on this subset, so their
kernel relation annihilates all original coordinates as required.

For any maximal minor, subtract the midpoint of each selected
coordinate interval times the constant row from that coordinate row.
Hadamard's determinant inequality gives

```text
|minor| <= (r+1)^((r+1)/2) product(w_i) / 2^r.
```

Since each minor is integral, the following bounds suffice:

| r | Largest possible product of selected widths | Bound on each minor |
|---|---:|---:|
| 1 | 4 | 4 |
| 2 | 4 | 5 |
| 3 | 2 | 4 |
| 4 | 1 | 3 |

There are `r+2` signed coefficients and their sum is zero. The side
with fewer nonzero coefficients has at most `floor((r+2)/2)` terms.
The common positive and negative mass is consequently at most
`4,10,8,9`, respectively. This proves (1), without a bound on the
locations of the coordinate intervals or an assumption of binary
coordinates.

## A count for persistent fixed-exponent families

Consider any family of nonzero Gaussian integers

```text
z_j(n) = c_j product_s H_s(n)^(a_js)
                       conjugate(H_s(n))^(e_s-a_js),
```

where the integer exponent vectors and nonzero Gaussian prefactors
`c_j` are fixed, all `c_j` have the same modulus, and all evaluated
factors are Gaussian integers. No Pell recurrence is assumed here.
Suppose the total active exponent width is at most four. Divide the
rows by their common Gaussian gcd. Assume their primitive radius
tends to infinity along a sequence, their rows are distinct, and
they lie in arcs of length at most `C sqrt(R)` along that sequence.

If `10 C^2 <= 8`, there are at most **five** rows.

First, identical exponent vectors would give a fixed quotient of
prefactors. Shrinking angular widths force that quotient to be one,
contrary to distinctness. The exponent vectors are therefore distinct.
If they were affinely dependent, write `b_j` for the coefficients in
(1). The exponent relation and
zero coefficient sum give the exact identity

```text
product_j z_j(n)^(b_j) = product_j c_j^(b_j),
```

where the `c_j` are still the fixed Gaussian prefactors.
The right side is fixed and has modulus
one. The left side tends to one because the angular width tends to
zero and `sum b_j=0`. The right side is thus exactly one. Dividing
all rows by their common Gaussian gcd preserves this identity.

There is consequently an equality between two products of at most
ten primitive arc points, with disjoint nonempty multisets of row
indices. The [angular moment theorem](angular_moment_product_rigidity.md)
excludes this when `10 C^2<=8`. Hence the exponent vectors are
affinely independent. At most four coordinates are active, so there
can be at most five vectors.

## Finite configurations with bounded prefactors

Persistence can be removed if the prefactor height is controlled.
Suppose an individual configuration has the form

```text
z_j = g c_j product_s H_s^(a_js) conjugate(H_s)^(e_s-a_js),
```

with nonzero Gaussian integers `g,c_j,H_s`, fixed integer exponents
within this configuration, total active width at most four, and all
`c_j` of one common modulus at most `H`, where `H>=1`. The points
are distinct, of common radius `R`, and occupy an arc of length at
most `C sqrt(R)`. Neither a recurrence nor pairwise coprimality of
the factors is assumed. If

```text
10 C^2 <= 8,       R > 50 C^2 H^20,                 (2)
```

then the configuration has at most **five** points.

First, repeated exponent vectors would give `z_i/z_j=c_i/c_j`.
Since distinct equal-modulus Gaussian integers have distance at
least `sqrt(2)`, this would imply

```text
sqrt(2)/H <= |z_i/z_j-1| <= C/sqrt(R),
R <= C^2 H^2/2,
```

contrary to (2). Thus the exponent vectors are distinct.

If they have an affine dependence, choose (1), with common positive
and negative mass `q<=10`. Let `A` and `B` be the products of the
prefactors on the positive and negative sides. They have the same
modulus, at most `H^q`, and factor cancellation gives

```text
product_j z_j^(b_j) = A/B.
```

If `A!=B`, equal-norm Gaussian separation and the angular interval
bound yield

```text
sqrt(2)/H^q <= |A/B-1| <= q Delta <= q C/sqrt(R),
R <= (q^2 C^2/2) H^(2q) <= 50 C^2 H^20.
```

No angular unwrapping is needed here: choose any lifts in the arc,
use `sum b_j=0` to bound their combination by `q Delta`, and then
use `|exp(i x)-1|<=|x|`. If `A=B`, the exact product identity is
excluded by rigidity through ten factors. Thus there is no affine
dependence. At most four active coordinates again give at most five
points.

### Arbitrary unit prefactors: no radius threshold

If every `c_j` is a Gaussian unit, set `H=1`. Under `10 C^2<=8`,
the same five-point bound holds for **every** radius. Otherwise six
distinct points would require arc length at least `5 sqrt(2)`, by
applying equal-norm separation to successive points. Hence

```text
R >= 50/C^2 > 50 C^2,
```

because `C^2<=4/5<1`. The finite result (2) then rules out those
six points. This gives a finite, nonasymptotic five-point bound for
every such width-four factor configuration, with arbitrary block
values, overlaps, common factors, and Gaussian units.

## Consequence for the quadratic-unit templates

The [fixed quadratic-unit theorem](quadratic_unit_template_bound.md)
proves active width at most four whenever `lambda>sqrt(11)`.
At `C=1/2`, its count therefore improves to **five** for every such
unit, including the canonical unit `9+4 sqrt(5)`.

For arbitrary quadratic units, that theorem permits width up to
twelve; the earlier bound of six at `C=1/2` is not improved by this
argument. No extraction theorem places arbitrary endpoint clusters
in a persistent fixed-exponent family of width four. In particular,
fixed-prefactor cancellation in the persistent-family argument uses
persistence. For isolated configurations the explicit height condition
(2), or the unit-prefactor hypothesis, supplies a different sufficient
condition; unrestricted changing prefactors are still not covered.

The circuit and family arguments are prose proofs. The determinant
table uses Hadamard's inequality and integer rounding. The equal-norm
Gaussian separation and relative prefactor-separation inequalities
are checked in `GaussianChain/DependentEndpoint.lean`; the circuit,
height threshold, and five-point theorems remain prose proofs.
