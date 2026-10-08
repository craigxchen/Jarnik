# Exact bounded search for four vertical cofactors

This note records a finite search for four distinct positive integers

```text
w_j = -d + i t_j,
```

with `d` in `{1,2,5,10}`.  For each tuple, the checker computes the Gaussian
lcm `A` exactly and the four quotients `v_j=A/w_j`.  It retains only tuples
for which every `v_j` is primitive for the ordinary integer content
`gcd(|Re(v_j)|,|Im(v_j)|)=1`.

The script is [`check_four_vertical_lcm_search.py`](check_four_vertical_lcm_search.py).
Its default run is reproducible and bounded:

* exhaustive `1 <= t_1 < ... < t_4 <= 40` for each of the four `d`;
* arithmetic progressions with largest `t <= 100`;
* 1,000 pseudorandom tuples with `t <= 1,000` for each `d`;
* targeted scaled configurations `(d,t)=(2,(6,8,16,42)),
  (5,(15,25,65,155)), (10,(30,40,80,210))`;
* the Pell orbit from [`vertical_lcm_pell_counterfamily.md`](vertical_lcm_pell_counterfamily.md), through 20 steps.

All Gaussian division uses integer nearest-quotient arithmetic.  No floating
point operation is used in gcd or lcm computation, so the large Pell values
are exact.  Ranking uses the exact square of the requested lcm ratio,

```text
( |A| / (min(t_j/d))^2 )^2
    = N(A) d^4 / min(t_j)^4,
```

For a ramified lift, the actual integral circle uses `B=(1+i)A`; its ratio is
reported separately by the checker.

## What counts as an actual factorization

The endpoint identity is `2P=v_j w_j`.  The minimal lcm `A` need not itself
be `2P`.  The checker tests the following exact possibilities:

* If both coordinates of `A` are even, take `B=A`, so `P=B/2` and
  `B/w_j=A/w_j` remains primitive.
* If both coordinates of `A` are odd and every quotient `A/w_j` has odd norm,
  take `B=(1+i)A`.  Its coordinates are even, and `(1+i)(A/w_j)` is still
  ordinary-primitive.  This is the possible ramified-prime lift.
* If `A` has mixed coordinate parity, every multiplier that makes `B` even
  contains `2`, which makes every complement non-primitive.  Such a tuple is
  reported as abstract lcm data only.

Thus replacing an odd lcm by `2A` is not accepted: it doubles every quotient.
The displayed `A` is always the exact lcm; `B` and the `B/w_j` quotients are
shown separately when a ramified lift is used.

## Best small exact configurations

The lowest lcm ratio among actual-factorization tuples in the bounded and
targeted checks is the following `d=2` tuple.

```text
w = (-2+6i, -2+8i, -2+16i, -2+42i)
A = 86-38i,                 N(A)=8840
B = A,                      P = 43-19i,       N(P)=2210
B/w = (-10-11i, -7-9i, -3-5i, -1-2i)
Q = (23-41i, 29-37i, 37-29i, 41-23i)
```

The four `Q` are distinct and all have norm `2210`.  Every complement is
ordinary-primitive.  The exact lcm ratio identity is

```text
( |A| / (min(t)/d)^2 )^2 = 8840/81,
|A| / (min(t)/d)^2       = sqrt(8840/81) = 10.4468082431...
```

The pair gcd norms for this tuple are
`4,20,8,4,68,52`; every triple gcd norm is `4`, and the four-way gcd norm is
`4`.  Its cotangent spread is `H=42-6=36`.

Scaling every coordinate by five gives a useful `d=10` instance:

```text
t = (30,40,80,210),       A=B=430-190i,
P = 215-95i,              N(P)=55250,
B/w = (-10-11i, -7-9i, -3-5i, -1-2i),
Q = (115-205i,145-185i,185-145i,205-115i).
```

Here the exact squared lcm ratio is `221000/81`, so the ratio is
`52.2340412157...`.  The pair gcd norms are
`100,500,200,100,1700,1300`; every triple gcd norm is `100`, and the
four-way gcd norm is `100`; its spread is `H=210-30=180`.  This is exactly
the scaled `d=2` pattern.

The best `d=1` actual tuple in the exhaustive box uses the ramified lift:

```text
t = (3,5,13,31),
A = -179-223i,
B = (1+i)A = 44-402i,
P = 22-201i,
B/w = (-125+27i, -79+7i, -31-i, -13-i).
```

Its exact squared lcm ratio is `81770/81`, ratio `31.7727268713...`; the
actual `B` ratio is `sqrt(163540/81)=44.9334212550...`.  All unlifted
quotients `A/w` have odd norm, and the `1+i` lift preserves primitivity;
declaring `A/2` to be the center would be invalid.

For `d=5`, the first actual tuple in the default box is also a ramified lift,
but it is much larger in ratio.  The arithmetic-progression scan finds
`t=(25,35,45,55)` with lcm ratio `255.0020...` and actual `B` ratio
`360.6272...`; scaling the `d=1` example by five gives `t=(15,25,65,155)`
with lcm ratio `158.8636...` and actual `B` ratio `224.6671...` (the latter
is a targeted check rather than part of the default `t<=40` exhaustive box).

## Pell-family check

The four-factor Pell orbit from the companion note was checked with exact
integer Gaussian arithmetic.  For steps `n=1,...,19`, all four `t_j` are
positive and distinct and all four lcm quotients are ordinary-primitive.  The
lcm ratio takes approximately

```text
81.6506661456...  (n = 1)
81.6311896062...  (n = 6,11,16,...)
182.5328890437... (the other checked indices).
```

The latter values stabilize to the displayed digits as `n` grows.  No trend
towards zero was observed.  The three-factor subfamily in
`vertical_lcm_pell_counterfamily.md` has a different behavior and can drive
an exponent above `3/2` to failure; that does not produce a four-factor
quadratic counterfamily here.

## Finite-search conclusions

The exhaustive `t<=40` run retained the best 12 actual tuples for each of
`d=(1,2,5,10)` (the checker also reports the total number of primitive lcm
tuples encountered), with abstract-only primitive candidates present.  Within
that box the best `d=2` actual tuple is `t=(6,8,10,36)`, with lcm ratio
`sqrt(44200/81)=23.359773...`; the lower `t=(6,8,16,42)` value above is an
explicit targeted check outside the box.  The arithmetic-progression and
random searches produced no smaller lcm ratio than that targeted tuple.

These computations do not prove a uniform bound or exclude larger sporadic
examples.  They provide finite evidence about exact realizability, gcd
sharing, and parity, not an asymptotic argument.
