# Critical odd moments force restricted root sign orders

Let `s>=1`, and let `b_1,...,b_(2s+1)` be nonzero real numbers with
pairwise distinct absolute values. Suppose

```
sum_j b_j^(2k+1)=0,       0<=k<s.                         (1)
```

Order the numbers by increasing absolute value. Their signs are necessarily

```
epsilon * (+,+,-,-,+,+,-,-,...),
```

truncated after `2s+1` entries, where `epsilon` is either sign. Thus there
are exactly `s` sign changes, and the majority sign occurs `s+1` times.
For eleven entries the only possibilities are `++--++--++-` and its
overall sign negation.

This strengthens the general Chebyshev sign-variation condition at the
critical support size. It is a necessary condition on the full nonlinear
moment equations, not a sufficient condition on magnitudes. It does not
give a uniform lattice-circle bound or by itself exclude the degree-22
odd-column code.

## Proof by real polynomial fibers

Newton's identities applied to the monic root polynomial give

```
P(t)=product_j(t-b_j)=F(t)+r,
F(t)=t Q(t^2) monic and odd,        r!=0.                 (2)
```

Indeed, the odd elementary symmetric coefficients of indices
`1,3,...,2s-1` vanish successively; only the constant coefficient of
even degree remains. Nonzero roots imply `r!=0`.

The derivative `F'=P'` has exactly `2s` distinct real roots by Rolle's
theorem: one lies in each interval between successive simple roots of
`P`. Since the derivative is even, its roots come in nonzero opposite
pairs. They divide the real line into `2s+1` monotonicity intervals.
Each interval contains exactly one root of `F+r`.

Reflection shows that `F-r=-P(-t)` also has `2s+1` simple real roots.
Consequently the image of every monotonicity interval contains both
levels `-r` and `r`, strictly away from its finite endpoint values.
Every intermediate fiber

```
P_u(t)=F(t)+u r,              -1<=u<=1,
```

therefore has one simple root on each such interval. In particular

```
F(t)=t product_(k=1)^s (t^2-d_k^2),
0<d_1<...<d_s.                                         (3)
```

For `u>0`, none of the roots is zero, and no two are opposite: the
equations `F(x)+u r=F(-x)+u r=0` would imply `2u r=0`.
Roots move continuously along the monotonicity intervals, so their
signs and their order by absolute value cannot change as `u` moves
from a small positive number to one.

It remains to compute the order for small positive `u`. The root
near zero has sign

```
epsilon = -sign(r) (-1)^s,
```

because `sign(F'(0))=(-1)^s`. The derivatives at the two roots
`+d_k` and `-d_k` are equal and have sign `(-1)^(s-k)`.
Both roots move in the direction `-r/F'(d_k)`. Of this pair, the
one with smaller absolute value therefore has sign
`sign(r)(-1)^(s-k)`, and the other has its opposite sign.

Beginning with the central root, the resulting sequence is

```
epsilon, (epsilon,-epsilon), (-epsilon,epsilon), ...,
```

which is exactly the displayed paired sign order. Continuity transfers
it to the original roots at `u=1`, proving the claim.

## Application to minimal odd-moment certificates

Suppose `S(a^(2k+1))` is constant across rows for `0<=k<s`, and
`c=S^t mu` with `sum mu=0` has exactly `2s+1` nonzero entries, each
equal to `+1` or `-1`. Set `b_j=c_j a_j` on its support. Equation (1)
holds, so after ordering that support by increasing `|a_j|`, the signs
`c_j sign(a_j)` must have the paired pattern above.

Every such certificate must use the same global ordering and the same
column signs of `a`. Separate choices for each certificate are not
allowed. This supplies a finite combinatorial necessary test beyond
[the five-variation test](tensor24_odd_column_variation.md), while the
remaining coefficient magnitudes still have to solve the moment
equations. No assertion about their feasibility is made here.

## One additional root permits only two patterns

There is also a useful successor case. Suppose there are `2s+2` roots,
with the same nonzero, distinct-absolute-value and moment hypotheses.
Their signs in increasing absolute order must, up to overall sign negation,
have one of the following two forms:

```
A: (+,-,-,+,+,-,-,+,+,...),       truncated to length 2s+2;
B: (+,+) followed by pattern A of length 2s.
```

For twelve roots these are `+--++--++--+` and `+++--++--++-`.
Pattern A has equally many positive and negative signs; pattern B has
two more of its initial sign.

Here is a proof. Put `m=s+1>=2`. Newton's identities now give

```
P(t)=Q(t^2)+r t,       deg Q=m,       Q(0)!=0, r!=0.
```

If `r=0`, the roots would occur in opposite pairs, which is excluded.
The polynomial `Q(x)^2-r^2 x` has precisely the `2m` distinct positive
roots `b_j^2`. Write them as `x_1<...<x_(2m)`. Its negative intervals
are `(x_1,x_2),...,(x_(2m-1),x_(2m))`.

On `x>0`, put `H(x)=Q(x)/sqrt(x)`. On every intervening positive
interval, `H` has the same sign at both ends and magnitude greater
than `|r|` inside; hence there is a critical point. This gives `m-1`
critical points. Each negative interval whose two endpoints have the
same `H` sign supplies another critical point. If `v` negative intervals
instead have opposite endpoint signs, this gives at least `2m-1-v`
distinct positive zeros of `H'`.

The numerator `2x Q'(x)-Q(x)` of `H'` has degree `m`, so `v>=m-1`.
Each of those `v` intervals contains a positive root of `Q` of odd
multiplicity. Thus either `Q` has `m` simple positive roots, one in
each negative interval, or it has `m-1` simple positive roots and
one negative root. A remaining positive or repeated root in the latter
case would contradict the endpoint sign counts; a zero root is excluded.

In the first case, the endpoint signs of `H` change on each negative
interval and remain the same across each intervening gap. Since the
original root sign satisfies `sign(b)=-sign(Q(b^2)/r)`, this gives A.

In the second case write the roots of `Q` as `-c,d_1,...,d_(m-1)`,
with `c>0` and `0<d_1<...<d_(m-1)`. For `x>d_1` away from these roots,

```
(log|H|)' = 1/(x+c) + sum_k 1/(x-d_k) - 1/(2x),
(log|H|)'' < -1/(x-d_1)^2 + 1/(2x^2) < 0.
```

Thus each interval between consecutive positive roots of `Q` has
exactly one critical point of `H`. After its last positive root,
the first derivative of `log|H|` is strictly positive, since
`sum_k 1/(x-d_k)>(m-1)/x>=1/x`. The one negative interval with equal
endpoint signs must consequently occur before `d_1`: between two
positive roots it would require three critical points, and after
the last it would require one. Its two equal signs, followed by the
first endpoint of the next interval, give three equal initial signs.
Thereafter signs follow A, proving B.

This successor restriction applies to unit-coefficient certificates of
support twelve in the degree-22 code, alongside the eleven-term rule.
For arbitrary sign matrices, simultaneous ordering and magnitudes still
require analysis. For the odd-positive tensor family, the
[simultaneous eleven-term test alone](tensor24_odd_column_family_critical_order.md)
now excludes the full nonlinear moment system.

## Exact verification

[The checker](check_critical_odd_moment_root_order.py) constructs rational
polynomials `F+r` for both signs of `r` and `1<=s<=8`. Disjoint rational
intervals certify all roots are real and simple. Rational bisection
certifies their absolute ordering and signs; Newton's identities check
the vanishing odd moments without approximate root calculations.
Additional checks cover successor pattern A in degrees four through
eighteen, exact rational examples of pattern B in degrees four and six,
and a degree-twelve pattern-B polynomial whose roots are isolated using
rational intervals. The proofs cover all tuples with the stated hypotheses.
