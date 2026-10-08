# Real endpoint configurations pass all odd Vandermonde lower bounds

The classical odd-cardinality determinant inequalities, even when imposed on
every subset, do not by themselves force a uniform number of points on an
endpoint arc.  Equally spaced real points satisfy all of them while their
cardinality tends to infinity as slowly as desired below
`log(R)/log log(R)`.  These configurations are not lattice points; they test
only the strength of the geometric determinant consequences.

## 1. The family of inequalities being tested

For `n=2k+1` distinct integer points on an origin circle of radius `R`, the
classical Gaussian Laurent-Vandermonde argument gives

```text
product_(i<j) |P_i-P_j| >= 2^(k^2) R^(k^2).          (1)
```

The powers agree with the familiar first cases:

```text
n=3:  product of the three chords >= 2R,
n=5:  product of the ten chords >= 16R^4.
```

Here is a self-contained determinant normalization.  Evaluate the integer
real-basis columns

```text
1;  x,y;  x^2,xy; ...; x^k,x^(k-1)y
```

at the `2k+1` points.  Their determinant `D` is a nonzero integer: after
restriction to the circle these columns span the real trigonometric
polynomials of degree at most `k`, and a nonzero such polynomial has at most
`2k` zeros on the circle.  At degree `j`, change the top harmonic pair to
`z^j,̅z^j`; terms of lower harmonic degree are eliminated using earlier
columns and `x^2+y^2=R^2`.  The absolute determinant of this pairwise column
change is `2^(1-2j)`.  Multiplication over `j=1,...,k` gives `2^(-k^2)`.

For the Laurent columns `z^(-k),...,1,...,z^k`, the Vandermonde determinant
has absolute value equal to the chord product divided by `R^(k(2k+1))`.
Replacing its negative powers by `̅z^j/R^(2j)` contributes
`R^(k(k+1))`; hence

```text
|D|=product_(i<j)|P_i-P_j|/(2^(k^2)R^(k^2)).
```

Since `|D|>=1`, this proves (1), including its power of two.  In the rest of
the note (1) is imposed as an axiom on every odd subset of a real point
configuration.

## 2. Equally spaced points on an endpoint arc

Fix `C>0`.  For an integer `m>=3`, put `h=C sqrt(R)/(m-1)` and take points
`P_0,...,P_(m-1)` at arc-length parameters

```text
0,h,2h,...,(m-1)h
```

on a circle of radius `R`.  Their total arc length is exactly `C sqrt(R)`.
For `R` sufficiently large in terms of `C`, this is less than a semicircle.
If two indices differ by `d`, their chord length satisfies

```text
|P_i-P_j|=2R sin(dh/(2R)) >= (2/pi) d h,             (2)
```

using `sin x>=2x/pi` for `0<=x<=pi/2`.

Consider any subset of `n` indices `j_0<...<j_(n-1)`.  Since
`j_b-j_a>=b-a`, multiplying (2) over its pairs gives

```text
product_(a<b)|P_(j_a)-P_(j_b)|
 >= ((2C)/(pi(m-1)))^p R^(p/2) V_n,                  (3)
```

where

```text
p=n(n-1)/2,
V_n=product_(0<=a<b<n)(b-a)=product_(d=1)^(n-1)d^(n-d)>=1.
```

Thus consecutive indices minimize the elementary lower bound; no separately
chosen sparse odd subset is worse.

## 3. Simultaneous verification for every odd subset

Let `n=2k+1`, so `p=k(2k+1)`.  After discarding the helpful factor `V_n`,
(3) implies (1) provided

```text
(k/2) log R
 >= k^2 log 2 + k(2k+1) log(pi(m-1)/(2C)).            (4)
```

After division by `k`, it is enough that

```text
(1/2)log R
 >= k log 2 +(2k+1) log_+(pi(m-1)/(2C)),             (5)
```

where replacing the logarithm by its positive part only strengthens the
requirement.  The right side increases with `k`, so it suffices to check
`k<= (m-1)/2` at the endpoint.  A convenient uniform sufficient condition is

```text
(1/2)log R >= m[log 2+log_+(pi(m-1)/(2C))].           (6)
```

Consequently, if

```text
m=m(R) -> infinity,
m(R) log m(R)=o(log R),                               (7)
```

then for all sufficiently large `R`, every odd subset of the `m(R)` equally
spaced real points satisfies its corresponding inequality (1).

More quantitatively, for every fixed `0<c<1/2`,

```text
m(R)=floor(c log R/log log R)                         (8)
```

satisfies (6) for all sufficiently large `R`.  Indeed
`m log m=(c+o(1))log R`, while the additional `m log 2` and fixed constants
are `o(log R)`.

## 4. The full shape-dependent four-cluster inequality already follows

Take any four ordered points spanning an arc of length `S`, and a distinct
fifth point at distance `d` from an anchor among the four.  The product of the
six chords among the four base points is at most

```text
S^6/16.                                               (9)
```

Indeed, if the three intervening arc gaps are `a,b,c`, the six chords are at
most their corresponding arc distances, whose product is
`abc(a+b)(b+c)(a+b+c)`.  This is at most `abc S^3<=S^6/27`, which is stronger
than (9).

Apply the `n=5` case of (1) to these five points.  Each of the four remaining
chords, from the fifth point to a base point, is at most `d+S`.  Consequently

```text
16R^4 <= (S^6/16)(d+S)^4,
d+S >= 4R/S^(3/2).                                   (10)
```

Thus every real configuration satisfying all five-point inequalities already
satisfies the full shape-dependent four-cluster/fifth-point inequality,
simultaneously for every `S`; no fixedness assumption on `K=S/R^(1/3)` is
needed.  In particular, substituting `S=K R^(1/3)` in (10) gives

```text
d+K R^(1/3) >= 4K^(-3/2)sqrt(R),                     (11)
```

the explicit local-isolation inequality.

For intuition, in the equally spaced construction even the fixed-`K`
condition is eventually vacuous.  Its consecutive arc spacing satisfies

```text
h/R^(1/3)=C R^(1/6)/(m-1) -> infinity                (12)
```

under (7), so four points cannot lie in a fixed multiple of `R^(1/3)`.

This is the precise compatibility between the two geometric inputs.  The
odd-subset Vandermonde bounds permit the points to spread almost evenly over
the full `sqrt(R)` scale, and their five-point member already contains every
shape-dependent inequality (10).

## 5. Scope

The construction does not give lattice points and is not a counterexample to
a uniform lattice-point theorem.  It proves a narrower logical statement:
all classical odd-subset chord-product lower bounds, together with the full
shape-dependent four-cluster/fifth-point inequality, are consistent
with an unbounded number of real points on an arc of length `C sqrt(R)`.
Therefore a uniform lattice bound cannot follow from those inequalities alone;
some additional arithmetic input must distinguish lattice configurations from
these real ones.
