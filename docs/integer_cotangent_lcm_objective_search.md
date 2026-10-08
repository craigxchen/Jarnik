# A bounded search for low all-edge cotangent height ratios

The all-edge lcm objective

```text
N/(A/L)^4
```

was searched by a divisor graph at one small normalization, `L=6`. The
search used every positive integer `X<=10000`, generated edges from
`d | y^2+L^2`, and enumerated all 193,836 four-cliques and 85,642
five-cliques. Among four finite coordinates it found the exact record

```text
L=6,
X=(1146,1257,1842,4104),
A/L=191,
N=2163838625,
N/(A/L)^4 = 2163838625/1330863361 ~= 1.6258908979.
```

Here `N` is the least squared radius from **all** anchor and finite-finite
pair denominators. The six finite-finite cotangents are

```text
-12978, -3033, -1590, -3958, -1812, -3342,
```

and the prime factorization of the all-edge lcm is

```text
N = 5^3 * 13 * 17 * 29 * 37 * 73.
```

The minimum anchor has

```text
S=A^2+L^2=1313352,
d=X-A=(111,696,2958),
S/d=(11832,1887,444).
```

Thus this record lies on the ordinary divisor-cofactor branch of the
four-offset equations. Its triangle invariant from the existing
four-offset note is `J=164250`. Complement reanchoring produces the equal
objective tuple

```text
L=6,
X=(1146,1590,3033,12978).
```

This is a finite arithmetic record, not an endpoint counterexample or a
global minimization theorem. The search supplies its divisor data but no
Pell parametrization or classification. A single record cannot establish
whether it belongs to a Pell family.

For comparison, the lowest five-finite tuple in the same small search window
is the already structured extension

```text
L=6,
X=(9,12,18,22,48),
N=325,
N/(A/L)^4=5200/81.
```

It contains the known four-finite set `(12,18,22,48)`. No claim is made that
either displayed tuple minimizes over other `L`, larger boxes, or all
clique sizes. The exact checker verifies every pair quotient and every edge
denominator before comparing the objective.

## Verification

[check_integer_cotangent_lcm_objective_search.py](check_integer_cotangent_lcm_objective_search.py)
performs the bounded `L=6`, `X<=10000` four- and five-clique searches. It
checks the pair cotangents, offset invariant and complement formula, and
compares the all-edge lcm of the two fixtures with an independent Gaussian
denominator calculation.
