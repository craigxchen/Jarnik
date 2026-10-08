# A certified real ten-factor moment solution with an outer-pair column

The six-factor ordered cut exclusion does **not** extend to all real
ten-factor critical moment solutions. There are ten positive real numbers
with distinct magnitudes and five distinct subsets whose first and third
power sums agree, whose fifth power sums are distinct, and for which one
source column separates the first and last rows in fifth-sum order from
all three middle rows.

The construction below is an exact existence proof by a rational
contraction certificate. It is not merely a floating-point solution.
It does not establish rational coefficients, Gaussian lattice points, or
an endpoint constant less than `1/2`. Its scale-invariant formal endpoint
statistic lies between `6.93` and `6.94`.

## Signed formulation and the exact box

Use the following five sign rows, with both rows and columns numbered
from zero:

```text
S = [ - - - + - - + - + +
      - + - - + - + + + -
      + - - - + - - + - +
      + + + - - - + - - +
      + + - - - + - - + - ].                         (1)
```

Put `A=(S_1-S_0,...,S_4-S_0)/2`, as a four-by-ten matrix. We require

```text
A a=0,                    A a^3=0,                  (2)
```

where powers are coordinatewise. Fix `a_5=1,a_9=-2/5` and introduce
`x=(a_4,a_6,a_7,a_8)`. All four linear equations in (2) hold identically
under the affine parametrization

```text
a_0=-1+x_0+x_1+x_2,        a_1=-7/5+x_0+x_2-x_3,
a_2=7/5-x_1+x_3,           a_3=-1+2x_0+2x_2-x_3,
a_4=x_0,                  a_5=1,
a_6=x_1,                  a_7=x_2,
a_8=x_3,                  a_9=-2/5.                  (3)
```

Define the rational center and box radius by

```text
x* = (-48557368900,927791797667,805289163286,174670899854)/10^12,
r=10^(-8),                ||x-x*||_infinity<=r.       (4)
```

There is a unique solution of the four remaining cubic equations in
this box. Its coefficient vector, displayed only for orientation, is
approximately

```text
a = ( .684523592054, -.817939105468, .646879102187,
      .338792688919, -.048557368900, 1,
      .927791797667, .805289163286, .174670899854, -.4 ). (5)
```

All subsequent claims concern the exact root in (4), not these rounded
coordinates.

## Rational contraction proof

Write (3) as `a=c+E x` and define the four cubic polynomials
`F(x)=A(c+Ex)^3`. Let `J*=DF(x*)` and `B=(J*)^(-1)`, calculated by
exact rational row reduction. Set `L_j=sum_k |E_jk|`, so
`|a_j(x)-a_j(x*)|<=L_j r` throughout the box. The derivative difference
has the explicit row-sum norm bound

```text
||DF(x)-J*||_infinity
 <= max_i sum_j 3 |A_ij| L_j
                 (2 |a_j(x*)| L_j r+L_j^2 r^2) = V.  (6)
```

The [checker](check_tenfactor_real_outercut_certificate.py) verifies,
using rational arithmetic only,

```text
eta=||BF(x*)||_infinity < 10^(-10),
q=||B||_infinity V < 1/1000,
eta+q r<r.                                             (7)
```

Consequently `T(x)=x-BF(x)` is a contraction of the closed box into
itself. Banach's theorem gives its unique fixed point, and invertibility
of `B` makes that fixed point an exact zero of `F`. This proves the
claimed existence without a numerical solver assumption.

## Strict fifth-sum order and positive coefficients

The intervals

```text
|a_j(x*)|-L_j r <= |a_j| <= |a_j(x*)|+L_j r
```

are pairwise disjoint and all have lower endpoint greater than `0.048`.
Moreover, every two absolute values differ by more than `0.0126`, as
certified by the interval endpoints. Thus nonvanishing and separation
hold with explicit margins.
For each signed fifth sum `f_i=sum_j S_ij a_j^5`, the mean-value theorem
bounds its displacement from the rational center by

```text
epsilon_5=sum_j 5(|a_j(x*)|+L_j r)^4 L_j r.             (8)
```

Exact rational comparisons of the centers with these errors certify

```text
f_2 < f_3 < f_1 < f_0 < f_4.                           (9)
```

Column six is `(+,+,-,+,-)` in the original row labels. It is therefore
`(-,+,+,+,-)` in the order (9), the required extreme-pair cut.

For a positive-coefficient subset formulation, put `b_j=|a_j|` and
`S'_ij=S_ij sign(a_j)`. The signs of `a` are fixed throughout (4).
Then for every odd integer `m`,

```text
sum_j S'_ij b_j^m = sum_j S_ij a_j^m.
```

Use the positive-entry subsets of `S'`. Equal signed first and third
sums are equivalent to equal subset sums, and signed fifth sums differ
from twice the subset fifth sums by one common constant. Column six has
`a_6>0`, so its pattern is unchanged. This proves the positive real
subset assertion stated at the beginning.

## Quantitative scope

Let `P=abs(product_j a_j)` and define the scale-invariant statistic

```text
C_formal=(f_4-f_2)/(5 sqrt(P))
        =(2/5)(largest subset fifth sum-smallest subset fifth sum)/sqrt(P).
```

Using the coefficient intervals and (8), the checker verifies by squaring
positive rational inequalities that

```text
6.93 < C_formal < 6.94.                               (10)
```

This is the usual coefficient-normalized expression for a critical
ten-factor formal family. It is not an asserted primitive lattice arc
constant: no rational or integral realization of these coefficients,
and no Gaussian gcd normalization, has been supplied.

For the literal complex products `z_i(T)=product_j(T+i S'_ij b_j)`, with
all five row-unit multipliers equal to one, the
two exact lower odd-moment equalities yield

```text
arg z_i(T)-arg z_k(T)
  = (f_i-f_k)/(5T^5)+O_b(T^(-7)).                      (11)
```

Thus their eventual real angular order is indeed (9), with a genuine
outer-pair factor column. These products need not be Gaussian integers.
In particular, (11) and the real-root construction cannot be used to
claim a growing arithmetic forced factor `K` on actual lattice circles.

The result precisely rules out extending the six-factor absence theorem
using only distinct positive magnitudes, the exact first and third
moment equations, and eventual angular order. Additional arithmetic or
a quantitative small-constant hypothesis would still be needed.
