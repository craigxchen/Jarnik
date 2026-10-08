# Aggregate parity for mixed layers

The even-layer certificate extends to a restricted mixed profile when parity
is imposed prime by prime rather than layer by layer. Singleton and triple
layers may occur, but at every odd split prime the total allocation exponent
across the eight rows is even. This does not cover an arbitrary odd profile.

Write `S`, `P`, `T`, and `B` for the total log-norm weights of minority sizes
`1`, `2`, `3`, and `4`, and `D` for complete common Gaussian content. Thus
`W=log R^2=D+S+P+T+B`. Assume a common source unit,
`Delta<=C/sqrt(R)`, `C<=2`, and eight distinct actual points. As in
[the even-layer proof](even_column_balanced_character_height.md), affine
allocation rigidity then supplies nonvanishing for all 35 balanced
characters. For each odd split prime `p`, let `a_i(p)` be its row allocation
exponent and assume

```
sum_i a_i(p) == 0 (mod 2).                              (1)
```

The same condition is required after grouping nested layers belonging to one
prime; it is not enough that the total number of odd-cardinality layers be
even across different primes. Complete-content removal and a common Gaussian
scaling change the allocation total by a multiple of eight; reversing a
prime orientation replaces it by `8e-sum_i a_i(p)`. None changes this parity.

## Character height

For balanced `lambda`, the signed source exponent at `p` is
`n_lambda(p)=sum_i lambda_i a_i(p)`. Since every `lambda_i` is odd, (1) makes
`n_lambda(p)` even. Therefore the exact character product has a Gaussian
square-root identity

```
product_i z_i^lambda_i
  = u_lambda (Beta_lambda/conjugate(Beta_lambda))^2,
```

after conjugate-primitive reduction, with prime exponent `n_lambda(p)/2`.
Affine allocation nonvanishing excludes a unit `Beta_lambda` for every
nonzero balanced character in the stated endpoint range. Its integer
axis perpendicular coordinate gives

```
sum_p |n_lambda(p)|/2 * log p >= W/2 - 2 log C.           (2)
```

This is the exact prime-height inequality; it does not delete odd layers or
alter their phases.

## Layer decomposition and cancellation loss

Choose the unit threshold decomposition of each allocation exponent:
each integer threshold occurs once, with weight `log p`. Expand grouped
repetitions before using the following formula, so multiplicity is retained.
Define

```
K = sum_p log p * (
      sum_lambda sum_(layers at p) |lambda^t S_j|/4
      - (1/2) sum_lambda |sum_i lambda_i a_i(p)|
    ).
```

The triangle inequality makes `K>=0`. Across all 35 characters, the absolute
quarter-coefficient sums for cut sizes one through four are respectively
`35/2`, `15`, `45/2`, and `18`. These values hold for every individual cut,
irrespective of which row is designated first. Summing (2) therefore gives

```
(35/2)S + 15P + (45/2)T + 18B - K
    >= (35/2)(D+S+P+T+B) - 70 log C.
```

Rearranging proves the sharper exact certificate

```
B + 10T >= 5P + 35D + 2K - 140 log C.
```

In particular,

```
B + 10T >= 5P + 35D - 140 log C.                       (3)
```

The singleton coefficient cancels. Formula (3) is valid for arbitrary nested prime
powers satisfying (1), including overlapping rational support. The sharper
form above retains `K` on the favorable side; dropping `K` gives (3).
For arbitrary row units, assume `C<=sqrt(2)` and even total row-unit parity;
the same argument replaces every `log C` by `log(sqrt(2)C)`.

This extends the even theorem within the aggregate-parity subclass, but
does not control the full intrinsic target `7S+2P-T-B`.
For the singleton-only allocation `(0,x,x,x,x,x,x,x+y)`, direct character
averaging gives `K_p=15 min(x,y) log p`. It detects cancellation between
opposing singleton layers. But `(0,0,0,0,0,0,0,2)` has singleton weight
`2 log p`, satisfies aggregate parity, and has `K_p=0`. A uniform lower
bound for `K` proportional to singleton weight is therefore false for
these allocation profiles. No endpoint realization of this example is
being asserted.

## Why balanced-character ratios do not repair arbitrary parity

It is tempting to pair two balanced characters so their odd exponents become
integral. Every Hamming-distance-two pair differs by swapping two rows, and
its ratio is exactly a square of a two-point ratio `z_i/z_j`. Hence this
construction recovers only the existing pair separation constraint; it cannot
supply the missing eight-row balanced height.

For completeness, with the 35 characters normalized by fixing row zero to be
positive, the distance-two pair totals for a cut of size `r` are

```
                 r=1   r=2   r=3   r=4
cut contains 0     0    60   100   120
cut avoids 0      60   100   120   120
```

The dependence on whether the cut contains the anchored row is why averaging
only the first representative cut gives a false mixed-layer theorem.
Within the 35 chosen representatives there are 210 pairs at distance two
and 70 at distance six. For the former use `lambda-mu`; for the latter use
`lambda+mu`. Each resulting coefficient vector is twice an ordinary
row difference, and every one of the 28 row pairs occurs ten times.
The cut totals are thus `(70,120,150,160)=10r(8-r)`. These give only the
existing pair-average constraint. Equivalently, all 70 balanced vectors
give 560 distance-two edges, twenty copies of each row pair.

Thus (3) is the valid bounded mixed extension in the aggregate-parity
subclass. Outside this subclass, recovering `5P+35D` with only a finite
penalty `k(S+T)` remains open. The single-character phase can be shifted by
the odd Gaussian factor, so exact-height deletion is not legitimate. The
[fixed-twist Pell example](odd_layer_character_quadratic_twist_obstruction.md)
shows why a uniform integer perpendicular-coordinate gap cannot simply be
assumed in that branch.

[The exact checker](check_aggregate_parity_balanced_character_height.py)
checks all 254 nonconstant sign columns, the complete anchored and folded
pair incidence counts, aggregate-parity Gaussian identities on actual
nested products, and the source-prime cancellation formula. These are
arithmetic fixtures, not generated endpoint clusters.
