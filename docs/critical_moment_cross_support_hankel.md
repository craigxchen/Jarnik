# A sign-free Hankel condition on every critical cross support

This is a necessary condition for the full punctured orthogonal moment
template, not a sufficiency theorem or an exclusion of all templates.
It concerns each individual cross-group differing support and is distinct
from projectors on the span of several sign rows. Its inputs are only
the positive magnitudes on that support.

Let `D` be a cross support, `d=|D|=2s+1`, and write
`z_h=S_ih a_h`, `x_h=a_h^2`. The full moment equations give

```
sum_D z_h^(2r+1)=0,       0<=r<s,
R(X)=product_D(X-z_h)=F(X)+r,       F odd, r!=0.
```

For `0<=k<=s`, let `H_k` be the span in `R^d` of the vectors
`(x_h^j)_h`, `1<=j<=k`, and put

```
p_k=1-Proj_(H_k) 1.
```

Then the sharp conditions are

```
||p_s||^2=1,
||p_k||^2>1                 for 0<=k<s.               (1)
```

In explicit moment form, set `m_l=sum_D x_h^l`. For `k>=1`, define

```
M_k=(m_(i+j))_(1<=i,j<=k),      h_k=(m_1,...,m_k)^t.
```

Distinct positive nodes make `M_k` positive definite. Equation (1)
becomes the sign-free necessary test

```
h_s^t M_s^-1 h_s=2s,
h_k^t M_k^-1 h_k<2s          for 1<=k<s.              (2)
```

Equivalently,
`det([[d-1,h_s^t],[h_s,M_s]])=0`; the corresponding determinants
with `k<s` are positive. These conditions follow from the full
residual-polynomial identities and are not asserted to be logically
independent of those identities. Their sufficiency is not claimed.

## Proof

Let `v_h=1/z_h`. The forward odd moments give
`<v,x^j>=sum_D z_h^(2j-1)=0` for `1<=j<=s`, so `v` is orthogonal
to `H_s`. At a root of `R`, the equation `F(z_h)=-r` gives

```
v_h=-(F(z_h)/z_h)/r.
```

The right side is a polynomial in `x_h` of degree `s`. Thus `v` lies
in `span(1,H_s)` and is a nonzero scalar multiple of `p_s`.
The reciprocal identity proved in
[the reciprocal Gram note](critical_moment_reciprocal_gram.md) is
`<v,1>^2=||v||^2`. Since `<p_s,1>=||p_s||^2`, substitution gives
`||p_s||^2=1`.

Projection onto nested spaces gives `||p_k||^2>=||p_s||^2`. Equality
would make `p_k=p_s`, so `v` would be the values of a polynomial `q`
of degree at most `k`. Then all the `d` distinct signed roots `z_h`
would solve `X q(X^2)-1=0`, a nonzero polynomial of degree at most
`2k+1<d`, a contradiction. This proves strictness. The normal equations
for the orthogonal projector give
`||p_k||^2=d-h_k^t M_k^-1 h_k`, proving (2).

## What the first moment and reciprocal identity alone imply

Even without the higher odd moments, suppose a cross support satisfies
`sum_D z_h=0` and `<v,1>^2=||v||^2`. Then `v` is orthogonal to
`x=(a_h^2)_h`. Cauchy--Schwarz applied to
`<v,1>=<v,p_1>` gives `||p_1||^2>=1`, hence

```
(sum_D a_h^2)^2 <= (d-1) sum_D a_h^4.                 (3)
```

For `d>=5`, the inequality is strict. Equality would make
`1/z_h=lambda(1-m_1 z_h^2/m_2)` for every `h`, so the `d` distinct
signed roots would solve one nonzero cubic. For `d=3`, the fixture
`z=(-3,1,2)` attains equality.

The exact degree-five fixture `z=(1,5,9,-7,-8)` has zero first and
third moments and distinct absolute values. It gives
`h_2^t M_2^-1 h_2=4`, whereas `h_1^t M_1^-1 h_1=1100/311<4`.
The [exact checker](check_critical_moment_cross_support_hankel.py)
verifies both fixtures, the residual shape, and the reciprocal and
projection identities. These are individual cross-pair fixtures, not
full eight-row solutions.

## Equal-weight quadrature and a bound against excessive concentration

There is a further consequence of the complete residual polynomial,
in addition to the projection equations. Write `F(X)=X H(X^2)`, with
`H` monic of degree `s`, and set `U(Y)=product_D(Y-x_h)`. Then

```
U(Y)=Y H(Y)^2-r^2.                                   (4)
```

Indeed, since `d` is odd,
`R(X)R(-X)=-U(X^2)=r^2-F(X)^2`.

The `s` roots `xi_1,...,xi_s` of `H` are distinct and positive.
To see this without assuming it, both `F+r=R` and
`F-r=-R(-X)` have `d` distinct real roots. Their shared derivative
`F'` has `d-1` distinct real roots by Rolle's theorem. At every local
maximum, both real-rooted level polynomials must be positive; at every
local minimum they must be negative. Thus the local maxima of `F`
are greater than `|r|` and the local minima less than `-|r|`.
Every intermediate level, including zero, has `d` distinct real roots.
The odd polynomial `F` therefore has its simple root zero and `s`
opposite nonzero root pairs, proving the assertion about `H`.

The monic degree-`d` polynomials `U(Y)` and `Y H(Y)^2` differ only
in their constant coefficient. Newton's identities consequently give
the exact equal-weight quadrature equations

```
m_l = 2 sum_(j=1)^s xi_j^l,           1<=l<=2s.       (5)
```

Including degree zero, this says that the `d` unit atoms at the
`x_h` have the same moments through degree `2s` as one unit atom at
zero and weight two at each `xi_j`.

For `s>=2` and positive integers `r_1,r_2` with `r_1+r_2<=2s`,
positivity of the distinct quadrature nodes now gives

```
m_(r_1) m_(r_2) > 2 m_(r_1+r_2).                    (6)
```

The difference is exactly
`4 sum_(i!=j) xi_i^(r_1) xi_j^(r_2)>0`. For `s=1` it is zero.
In particular, combining (6) with (2) at degree one gives

```
2 < m_1^2/m_2 < 2s                 when s>=2.         (7)
```

Thus the full equations impose both sides of this concentration
bound. More generally (6) gives magnitude-only tests up to the even
power `a^(4s)`. These tests are consequences of the full polynomial
identity; no logical independence from every Hankel condition is asserted.

For example, let `L` and `B` be the largest and second-largest squared
magnitudes in a cross support. Taking `r_1=r_2=s` gives
`m_s>sqrt(2) L^s`, whereas `m_s<=L^s+(d-1)B^s`. Hence

```
L/B < [(d-1)/(sqrt(2)-1)]^(1/s)       when s>=2.       (8)
```

The checker also verifies (4)--(6) exactly on the two local fixtures,
using Newton power sums of `H` without numerical root approximations.
All sixteen cross supports of a full template must satisfy these
conditions simultaneously, but no incompatibility of those overlapping
conditions has been proved here.
