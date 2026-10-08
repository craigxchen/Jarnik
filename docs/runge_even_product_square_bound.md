# An elementary bound for square values of an even product

Let `k,D>=1` be integers, and let `c_1,...,c_(2k)` be distinct positive
integers at most `D`. Set `P(X)=product_j (X-c_j)`. If an integer
`Z>D` satisfies `P(Z)=Y^2` for an integer `Y`, then

```text
Z <= (8kD)^(2k+1).                                  (1)
```

The exponent is linear in the number of factors. Both the integer roots
and their distinctness matter. This elementary polynomial square-root
argument applies to moving coefficients, with all constants displayed.

## 1. Coefficients and one common dyadic denominator

Expand at zero

```text
f(t)=product_j (1-c_j t)^(1/2)=sum_(a>=0) s_a t^a,
s_0=1,
S(X)=sum_(a=0)^k s_a X^(k-a),
F=2kD,       q=2^(2k-1).
```

The binomial coefficients of `(1-t)^(1/2)` have absolute value at most
one. There are `binom(2k+a-1,a)` weak compositions of `a` into `2k`
parts. The sorted words corresponding to these compositions form a subset
of all words of length `a` on `2k` symbols. Consequently, for every `a>=0`,

```text
|s_a| <= binom(2k+a-1,a)D^a <= (2kD)^a=F^a.          (2)
```

Write `product_j(1-c_j t)=sum p_a t^a`, with integral `p_a`.
The recurrence

```text
2s_a=p_a-sum_(b=1)^(a-1) s_b s_(a-b)
```

shows by induction that the denominator of `s_a` divides `2^(2a-1)`
for `a>=1`. A product in the sum has denominator dividing `2^(2a-2)`;
division by two proves the induction step. Thus `qS` is an integer
polynomial. The denominator grows exponentially in `k`, independently
of the sizes of the roots.

## 2. Rounding the square root leaves a polynomial equation

For real `Z>2F`, the convergent series is the positive square root of
`P(Z)`. Equation (2) gives

```text
|sqrt(P(Z))-S(Z)|
 <= Z^k sum_(a>=k+1) (F/Z)^a
 <= 2F^(k+1)/Z.                                    (3)
```

If also `Z>2qF^(k+1)` and `P(Z)` is an integer square, (3) places
its nonnegative square root at distance less than `1/q` from `S(Z)`.
Both belong to `q^(-1) Z`, so they are equal. Therefore

```text
R(Z)=0,       R(X)=q^2(P(X)-S(X)^2) in Z[X].         (4)
```

All terms of degrees at least `k` cancel in `P-S^2`, so its degree is
at most `k-1`. It is nonzero, because a polynomial with distinct simple
roots cannot be a square. For `k=1`, (4) is already impossible since
`R` is a nonzero constant.

For the remaining cases, (2) bounds the coefficient at degree `2k-a`
of `S^2` by `(k+1)F^a`. The corresponding coefficient of `P` has
absolute value at most `F^a`. Thus the coefficient height of `R` is
at most

```text
q^2(k+2)F^(2k) <= (k+2)(8kD)^(2k).                 (5)
```

A nonzero integer polynomial with coefficient height `B` has every
root of absolute value at most `1+B`: divide by its nonzero leading
coefficient, whose absolute value is at least one, and bound the other
terms by a geometric series. Equations (4)--(5) imply

```text
Z <= 1+(k+2)(8kD)^(2k) <= (8kD)^(2k+1).            (6)
```

The last bound also exceeds `2F` and `2qF^(k+1)`. Hence assuming that
`Z` exceeds it gives a contradiction in every preliminary range. This
proves (1), including the strict inequalities required for rounding.

## 3. Why complementary factors must first be removed

If roots repeat in pairs, `P` can be a polynomial square and (4) can
vanish identically. For example `(X-5)^2(X-4)^2` is a square at every
integer. There is then no bound on `Z` from this argument.

In the [reciprocal squareclass application](reciprocal_squareclass_runge_growth.md),
equal roots are exactly complementary divisor labels. The application
selects distinct roots before applying (1), excluding that exception.

The [exact checker](check_reciprocal_squareclass_runge.py) tests the
dyadic coefficients, cancellation, nonzero remainder and coefficient
bounds, with a repeated-root control. Root, Astra and Luna independently
derived the linear exponent; the finite checks supplement the proof.
