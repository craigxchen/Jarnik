# Explicit magnitudes pass the scalar and determinant relaxations

This note concerns only the punctured orthogonal template and the
[multirow Christoffel inequalities](critical_moment_multirow_christoffel.md).
It does not construct an odd-moment solution or imply a bound for arbitrary
lattice circles. The full condition `K+Q<=I` is not asserted below.

There is an exact limitation of the scalar and `s`-subset determinant
relaxations: **for `s>=100`, they admit strictly increasing positive
magnitudes if and only if the centered row space contains no nonzero vector
supported on at most `s` columns.** Whenever this support condition holds,
one explicit choice, independent of the column labels, works:

```
n=4s+2,     theta_j=(2j-1)pi/(2n),
x_j=1-cos(theta_j)/20,       a_j=sqrt(x_j),       1<=j<=n.     (1)
```

The existing admissible compressed words for `s=2+20,805,120k`, `k>=1`,
satisfy the support condition. Thus neither all scalar inequalities nor
all the coupled `s`-subset determinant inequalities exclude those words.

## Every template passes the scalar inequalities for s>=5

Keep `V,C,Q,K,H` as in the Christoffel note. The row-space projector has

```
q_j <= 8/(4s-4)=2/(s-1),                              (2)
```

because every centered column has squared norm at most eight, and the
smallest positive eigenvalue of `CC^t` is `4s-4`.

Let `K_0` be the unweighted regression projector for polynomials of degree
at most `s-1` evaluated at the nodes `cos(theta_j)`. Discrete Chebyshev
orthogonality gives the exact diagonal formula

```
(K_0)_jj = [1+2 sum_(r=1)^(s-1) cos(r theta_j)^2]/n
           <= (2s-1)/(4s+2) < 1/2.                   (3)
```

Affine transformation of nodes does not change this polynomial evaluation
space. The actual matrix `U` adds row weights `sqrt(x_j)`. Since
`19/20<x_j<21/20`, comparison of weighted Gram matrices gives

```
K_jj <= (21/19)(K_0)_jj
     <= (21/19)(2s-1)/(4s+2) < 21/38 < 3/5.          (4)
```

Combining (2) and (4) yields `q_j+K_jj<1` for every `s>=5`.
After clearing the positive denominator `19(s-1)(4s+2)`, the difference
from one is `34s^2-127s-135`; at `s=5+t` this is
`34t^2+213t+80>0`. Consequently all scalar Christoffel inequalities hold
for (1), on every punctured orthogonal template in these degrees.

## The determinant inequalities reduce to a support test

Use the seven integer row differences

```
B_i=V_i-V_0,          1<=i<=7.
```

They span the centered row space. Set `G=BB^t`. Its determinant depends
only on the degree, not on the column labels:

```
D_s := det G = 8(4s+4)^6(4s-4).                      (5)
```

Indeed, write `W_i=e_i-e_0`, `lambda=4s+4`, and use
`G=lambda WW^t-(Wb)(Wb)^t`, `det(WW^t)=8`, and
`(Wb)^t(WW^t)^(-1)(Wb)=8`. The last identity holds because `b` lies
in the row space `1^perp` of `W`.

The matrix determinant lemma gives, for every `s`-subset `J`,

```
det(I_J-Q_J) = det(B_(Jc) B_(Jc)^t)/D_s.             (6)
```

The numerator vanishes precisely when a nonzero vector in the centered
row space is supported in `J`. If it vanishes, the determinant inequality
cannot hold for any positive distinct nodes: `det K_J` is a positive
weighted Vandermonde determinant divided by `det H`.

Otherwise the numerator in (6) is a positive integer, so

```
det(I_J-Q_J) >= 1/D_s.                              (7)
```

For the explicit nodes (1), Hadamard's inequality and (4) give

```
det K_J <= product_(j in J) K_jj < (3/5)^s.           (8)
```

For every integer `s>=100`,

```
(3/5)^s D_s < 1.                                   (9)
```

The inequality at `s=100` is checked by the exact integer comparison
`3^100 D_100 < 5^100`. To propagate it, use

```
D_(s+1)/D_s = [(s+2)/(s+1)]^6 s/(s-1)
             <= (100/99)^7 < 5/3.
```

Equations (7)--(9) prove every `s`-subset determinant inequality at once.
Together with the zero-determinant obstruction, this proves the stated
if-and-only-if result. It characterizes these relaxations, not the full
moment equations.

## The exact admissible word family passes the support test

This uses the established integer words in
[the even-q=2 graph note](critical_moment_compressed_graph_even_q2.md).
Their retained row inner products are zero across the two groups and
`-2` within either group, with diagonal `4s+2`. Adjoining the constant
column and the balanced group column therefore gives
`TT^t=(4s+4)I`, as required for the projector formulas.
In one active circulation, every particular mask of size one or seven
occurs 403,200 times, and every particular mask of size two through six
occurs 322,560 times. Constants do not occur. Counting the fixed base path
and signed cycle correction in the saved JSON gives total count ten,
minimum mask correction -6, and maximum mask correction 5.
Therefore in the degree `s=2+20,805,120k` word, each mask occurs at most

```
m=403,200k+5.                                       (10)
```

The following general fourth-moment bound turns (10) into a robust
support bound. For any `mu in 1^perp`, write `c_j=mu^t V_j`. Summing over
the 256 possible sign masks gives

```
sum_j c_j^4 <= 256m [3||mu||^4-2 sum_i mu_i^4]
             <= 768m ||mu||^4.
```

Thus for every `s`-subset `J`, Cauchy--Schwarz and the full row Gram yield

```
sum_(j in J) c_j^2 / sum_j c_j^2
  <= sqrt(768sm)/(4s-4) =: rho.                      (11)
```

For the values of `s,m` above, `rho<39/40` for every integer `k>=1`.
This is an exact quadratic inequality after squaring and clearing positive
denominators; the accompanying checker verifies its three positive
coefficients after substituting `k=1+t`. Hence no nonzero centered row
combination can be supported on `s` columns. In fact (11) implies the
stronger uniform bound

```
||Q_J|| < 39/40,       det(I_J-Q_J) > 40^(-7).
```

Since these degrees exceed 100, the explicit magnitudes (1) pass every
scalar and `s`-subset determinant relaxation, in the order of any of the
established admissible Euler words. The sign-pattern and discrete prefix
transitivity tests encoded by the graph continue to hold in this word
order. This does not assert that the actual polynomial values for the
assigned magnitudes realize the graph's predicted value orders, or that
their odd moments vanish.

[The exact checker](check_critical_moment_christoffel_relaxation_limit.py)
audits the arithmetic thresholds, reads the existing mask witnesses to
verify (10), and proves the `rho` bound for every `k>=1` by polynomial
coefficients. It also checks (5)--(7) on all 364 three-column subsets of
a partial Walsh fixture, including 144 zero and 220 positive complementary
Gram determinants. It does not enumerate a giant Euler word or claim that
the Chebyshev magnitudes satisfy the original equations.
