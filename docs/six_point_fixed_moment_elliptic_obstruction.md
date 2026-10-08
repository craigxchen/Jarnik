# Fixed six-point moment coefficients permit angular coalescence

There are infinitely many primitive Gaussian integer six-point configurations
whose angular spans tend to zero and whose primitive ordinary central moment
vector is always

```text
u = (5,-5,5,-5,-1,1).
```

This is an obstruction to obtaining an angular height bound from the central
coefficient height alone. It is **not** an endpoint counterexample: the family
has exact two-factor product collisions, and its least primitive radius
compensates for the shrinking angle. Its unbalanced conductor cuts are bounded
in terms of the residue normalization, so it does not have the full fair
six-label profile.

## 1. An exact rational family

Suppose rational `p,q` satisfy

```text
p^2 (99-q^2) = q^2 (3-q^2).                         (1)
```

Put

```text
Z1 = (1-pq) + i(p+q),
Z2 = (1+pq) + i(p-q),
Z3 = (1-q^2) + 10ip,
Z  = (Z1, conjugate(Z1), Z2, conjugate(Z2), Z3, conjugate(Z3)).
```

All six rows have squared modulus `D=(1+p^2)(1+q^2)`. For `Z3`, the
difference from `D` is exactly the left side of (1) minus its right side.
The vector `u` satisfies

```text
sum_i u_i Z_i^j = 0,       j=0,1,2.                 (2)
```

Indeed, conjugate pairing reduces the linear identity to
`5[(p+q)+(p-q)]-10p=0`. For the quadratic identity use

```text
(1-pq)(p+q)+(1+pq)(p-q)=2p(1-q^2).
```

For six distinct nonzero points on a circle, these ordinary real moment
identities have a one-dimensional kernel. To see this, conjugation supplies
the negative moments, and the five Laurent powers `-2,-1,0,1,2` have rank five
on six distinct nodes: multiply columns by the square of the node and use
Vandermonde rank. Since `gcd(u_i)=1`, (2) is the actual primitive central
integer vector, up to sign. This does not merely exhibit a nonprimitive
relation among the rows.

## 2. Infinitely many rational parameters approaching zero

Consider the nonsingular elliptic curve

```text
E: Y^2 = X(X-3)(X-99),       P=(1,14).
```

If `X=q^2` is a nonzero rational square and `X!=99`, put
`p=Y/(99-X)`. The elliptic equation immediately gives (1).

The point `P` has infinite order. Exact duplication gives

```text
2P = (5476/49, -135050/343).
```

The integral translation `X=x+34` gives the short integral model
`Y^2=x^3-3171x-68510`. The translated first coordinate of `2P` is
`3810/49`, which is not an integer. Lutz--Nagell therefore excludes torsion;
see [Milne, Elliptic Curves, II, Theorem 5.1](https://www.jmilne.org/math/Books/EC2.pdf).

Every positive multiple of `P` has a square first coordinate. Here is a
direct proof avoiding a descent assumption. If a nonvertical line
`Y=aX+b` meets this cubic at three points, the product of their first
coordinates is `b^2`, by the constant term of the intersection polynomial.
The third point has the same first coordinate as the sum of the first two.
Starting with `X(P)=1`, this proves the assertion inductively by adding `P`.
The initial tangent is allowed in this product argument. Subsequent
secants cannot be vertical, since `mP=-P` would contradict infinite order;
the other coincidence `mP=P` is equally excluded for `m>1`. Nor can a
nonzero multiple have first coordinate zero, since `(0,0)` is torsion.

The real curve has two components: the bounded oval `0<=X<=3`, and the
component containing infinity with `X>=99`. The point `P` is on the oval,
so `2P` lies in the identity component. That component is a compact
connected one-dimensional Lie group, hence a circle. One can see the last
description by integrating its nonvanishing invariant differential around
the component, which identifies it with `R/TZ` for its positive period `T`.
An infinite-order point on a circle generates a dense subgroup: every proper
closed subgroup is finite. Thus the positive multiples of `2P` are dense
in the identity component, and the odd positive multiples of `P` are dense
in the bounded oval.

In particular, there is a sequence of odd multiples approaching `(0,0)`
on its branch with `Y>0`. Take the positive rational square root `q` of
their first coordinate. Then

```text
q -> 0,       p/q -> 1/sqrt(33).
```

All six `Z_i` tend to `1`. Their imaginary coordinates divided by `q`
tend, in conjugate pairs, to

```text
±(1+1/sqrt(33)), ±(1/sqrt(33)-1), ±10/sqrt(33).
```

These six real numbers are distinct. Hence all six points are distinct
for sufficiently small nonzero `q`, and their angular span tends to zero.
Their actual pair cotangent magnitudes all tend to infinity.

## 3. Common primitive normalization and the radius obstruction

For each rational tuple, clear one common ordinary denominator and divide
all rows by their common Gaussian gcd. The resulting tuple `z` is primitive,
has equal positive modulus, and retains both (2) and its angular span.
Any Gaussian integer realization of the same Gaussian rational projective
tuple is an integral Gaussian multiple of this one: a Bezout combination
of the primitive coordinates shows that its common rational multiplier is
Gaussian integral. Thus this normalization has the least radius.

The construction has three identical pair products:

```text
Z1 conjugate(Z1) = Z2 conjugate(Z2) = Z3 conjugate(Z3) = D.
```

Common complex scaling preserves these equalities. Consequently the least
primitive tuple still has a nontrivial multiplicative rectangle. If `R` is
its radius and `Delta<pi` its angular span, the sharp
[rectangle theorem](multiplicative_rectangle_separation.md) gives

```text
Delta sqrt(R) >= 2 sqrt(2),
N=R^2 >= 64/Delta^4.                                (3)
```

Thus bounded moment coefficients and angular coalescence coexist with
precisely the radius growth an endpoint argument must retain. The family
does not disprove a height bound that includes the full least radius.

## 4. Exact failure of the full fair cut hypothesis

Choose an integer cotangent normalization `L` for the six labels. At a split
prime `p` not dividing `2L`, write the sorted physical allocation exponents
as `e1<=...<=e6`, permuting the associated coefficient entries `u_i`
along with them. Formula (8) of
[the eigenvector height note](integer_cotangent_eigenvector_heights.md) gives

```text
v_p(u_i) = [sum_j |e_i-e_j| - H3(p)]/2,              (4)
H3(p) = e4+e5+e6-e1-e2-e3.
```

If also `p!=5`, every left side is zero. Every exponent must therefore
minimize the sum of absolute deviations, so all exponents lie in the median
interval `[e3,e4]`. Equivalently,

```text
e1=e2=e3,       e4=e5=e6.
```

Only a balanced `3|3` cut can occur at such a prime. At `p=5`, every
coefficient valuation is at most one. In a two-level profile, a minority
singleton contributes `2a` to its coefficient and a minority pair contributes
`a` to both coefficients. Good singleton cuts are therefore impossible;
good pair cuts are supported only at `5` and have exponent at most one.

There is also a bound including primes dividing the normalization. Let
`W_unbal` be the total log norm of all allocation layers with minority size
one or two; arbitrary multi-level profiles are decomposed into their ordered
successive layers. The same local content proof as in the height note gives

```text
|v_p(u_i) - [sum_j |e_i-e_j|-H3(p)]/2|
    <= 5 v_p(2L)                                    (5)
```

at every split prime. To verify the error, each of the five differences
in a barycentric denominator has extra valuation between zero and
`v_pi(2L)`. Subtracting the minimum coordinate valuation changes the raw
baseline by at most `5 v_pi(2L)`. Average this statement at the two conjugate
Gaussian primes; the vector is ordinary integral.

Summing the baseline in (5) over all six coordinates counts every unbalanced
layer with multiplicity two and every balanced layer with multiplicity zero.
Since `product_i |u_i|=5^4`, it follows that

```text
W_unbal <= 2 log(5) + 15 log(2L).                    (6)
```

A full fair six-label profile has six singleton cuts and fifteen pair cuts,
each of log norm `w+o(w)`, hence `W_unbal=21w+o(w)`. Equation (6) excludes
such a profile here when `log L=o(w)`. Therefore this obstruction does not
dispose of the outstanding full-fair problem, where the central coefficients
have height `exp(7w+o(w))` rather than being fixed.

## 5. Exact finite verification

Run `python3 docs/check_six_point_fixed_moment_elliptic.py`. The checker uses
only rational arithmetic. On the first fifteen multiples of `P`, it checks
square first coordinates, the parameter equation, distinct equal-radius rows,
the three moments, primitive barycentric central coefficients, common Gaussian
normalization, pair-product collisions, the independent all-edge ordinary
lcm formula for the least squared radius, and the central valuation error
bound at eleven split primes. Infinitude and coalescence follow
from the proof above, not from these finite samples.
