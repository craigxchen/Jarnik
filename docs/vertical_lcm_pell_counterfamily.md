# Three vertical cofactors can have only cubic-over-quadratic growth

There is an explicit infinite family of **three cofactors with the same
fixed real part**, with primitive ordinary-integer complements, for which

```text
min_j t_j/d ~ lambda^(2n),       |lcm_G(w_j)| ~ lambda^(3n).
```

Consequently no bound `|lcm_G(w_j)| >= c(d) (min_j t_j/d)^alpha`
with `alpha>3/2` and positive `c(d)` can hold for every three-cofactor
configuration of this kind. In particular the proposed quadratic
strengthening fails. Gaussian normalization gives the same obstruction
on full primitive circle tuples, with the three fixed contents `(1,2,1)`.

This is a four-point counterfamily, not an unbounded cluster. It says
nothing against a simultaneous inequality requiring more than three
cofactors, or against a uniform endpoint point-count bound.

## 1. An integral transformation preserving all three real parts

Identify `x+iy` with its coordinate column. Put

```text
g = 1+2i,               M = [[17,4],[4,1]],
a_0 = i,               b_0 = 1+3i,          c_0 = 4-3i,
a_n = M^n a_0,         b_n = M^n b_0,      c_n = M^n c_0.
```

The symmetric bilinear form `B(x,y)=Re(gxy)` has matrix

```text
S = [[1,-2],[-2,-1]],          M^T S M = S,       det M = 1.
```

Direct evaluation at index zero therefore proves, for every `n`,

```text
Re(g a_n b_n) = Re(g a_n c_n) = Re(g b_n c_n) = -5.       (1)
```

All three coordinate pairs remain primitive because `M` is unimodular.
Also `det(a_n,b_n)=-1`; in particular their Gaussian gcd is a unit.

Define

```text
P_n = g a_n b_n c_n,
(v_1,v_2,v_3) = (c_n,b_n,a_n),
(w_1,w_2,w_3) = (2g a_n b_n, 2g a_n c_n, 2g b_n c_n),
Q_j = P_n + 10 v_j.
```

Then, exactly,

```text
2P_n = v_j w_j,             w_j = -10+i t_j,
|Q_j| = |P_n|,              gcd(Re(Q_j-P_n),Im(Q_j-P_n)) = 10.
                                                               (2)
```

The norm identity follows either by expanding (2), or from
`2 Re(P_n conjugate(v_j))=-10 |v_j|^2`. These are genuine endpoint
factorizations with the stated ordinary contents.

The complementary-factor gcd is one, hence the exact lcm identity is

```text
lcm_G(w_1,w_2,w_3) = 2P_n                         up to a unit. (3)
```

There is no discarded pair-gcd product in (3).

## 2. Growth and one-sidedness

Let `lambda=9+4 sqrt(5)` and `s=sqrt(5)-2`. The expanding eigenvector of
`M` is `(1,s)`, with eigenvalue `lambda`. All three initial vectors have
positive projections on it:

```text
s > 0,                  1+3s > 0,                4-3s > 0.
```

Moreover `g(1+is)^2` is positive imaginary. Thus

```text
t_j ~ lambda^(2n),          |P_n| ~ lambda^(3n),           (4)
```

where `~` here means bounded above and below by fixed positive
constants. Positivity already holds at `n=1`. Explicitly the three
sequences have initial pairs

```text
(t_1(0),t_1(1)) = (-10,550),
(t_2(0),t_2(1)) = (20,1060),
(t_3(0),t_3(1)) = (70,7670),
```

and each satisfies `t(n+2)=322 t(n+1)-t(n)`. This follows by diagonalizing
`M`; the imaginary part of the mixed expanding/contracting term is zero.
It can equally be checked by direct matrix multiplication. The initial
values and recurrence prove positivity and increasing order for `n>=1`.

The `Q_j` lie on one side of `P_n` with angular distances
`2 arctan(10/t_j)`. Their containing arc consequently has length

```text
O(|P_n| / min_j t_j) = O(lambda^n) = O(|P_n|^(1/3)).     (5)
```

In particular its length divided by `sqrt(|P_n|)` tends to zero.
Combining (3) and (4), for every fixed `alpha>3/2`,

```text
|lcm_G(w_j)| / (min_j t_j/10)^alpha
     ~ lambda^((3-2 alpha)n) --> 0.                      (6)
```

All `t_j<=2|P_n|`, as required by their actual endpoint factorization.

The three imaginary parts are mutually comparable. On such a family,
the existing pairwise bound

```text
|lcm_G(w_1,w_2,w_3)|
 >= product_j |w_j| / product_(i<j) sqrt(d |t_i-t_j|)
```

has order at least `T^(3/2)/d^(3/2)` when all `t_j` and their spread are
`O(T)` and `t_j>=T`. This family's exact lcm has order `T^(3/2)` with
`d=10`. Thus that exponent is sharp for three comparable same-content
cofactors, even after imposing primitive complements. Higher gcds
cannot universally improve its exponent in this scope.

## 3. Exact Gaussian normalization

The whole four-point tuple has Gaussian gcd `h=5(1+i)`, up to a unit,
at every index. Here is a proof retaining the ordinary contents.

Since `gcd_G(a_n,b_n)=1`,

```text
gcd_G(P_n,Q_1,Q_2,Q_3) = gcd_G(P_n,10).                   (7)
```

Modulo two, `M` is the identity. Thus `a_n,c_n` have even real and odd
imaginary coordinates, while `b_n` has both coordinates odd. The latter
has exactly one factor `1+i`, and the former two and `g` have none.
Therefore the `1+i`-valuation of `P_n` is exactly one.

For the prime five, write

```text
M=9I+4K,            K=[[2,1],[1,-2]],            K^2=5I.
```

Since `K b_0=(5,-5)`, modulo five one has
`b_n=(-1)^n b_0`. The factor `b_0=1+3i` is divisible by `conjugate(g)`.
Thus `g b_n` is divisible by the ordinary integer five. It follows that
`5|P_n`. Equation (7), whose other argument is ten, now proves the
asserted exact gcd `h`.

Divide all four points by the fixed `h`, and write the resulting points
as `P'_n,Q'_j`. Their chord differences are

```text
Q'_j-P'_n = (1-i) v_j.
```

For a primitive ordinary vector, multiplication by `1-i` creates
ordinary content two precisely when both original coordinates are odd;
otherwise its ordinary content is one. Thus the three contents are
exactly

```text
(d'_1,d'_2,d'_3) = (1,2,1).                              (8)
```

Their genuine endpoint cofactors are

```text
w'_1=w_1/10,             w'_2=w_2/5,             w'_3=w_3/10,
Im(w'_j)/d'_j = t_j/10.                                  (9)
```

The primitive complementary vectors have Gaussian gcd one: the first
and third have gcd associated to `1-i`, while the second is odd in
Gaussian norm. Consequently

```text
lcm_G(w'_1,w'_2,w'_3)=2P'_n                  up to a unit. (10)
```

Equations (6) and (9)--(10) give the same failure of every exponent
`alpha>3/2` on these **full Gaussian-primitive tuples**, with contents
fixed independently of the radius. The normalization does not preserve
equal contents; that distinction is essential.

## 4. Four same-content factors can attain quadratic growth

There is also a four-cofactor variant, which attains rather than
contradicts a proposed quadratic lower bound. Put

```text
e_0=2+i,                  e_n=M^n e_0,
A_n=2g a_n b_n c_n e_n,
(w_1,w_2,w_3,w_4)=2g (a_n b_n,a_n c_n,b_n c_n,a_n e_n).
```

Their real parts are all `-10`, by invariance of `B`. The complements
`A_n/w_j` are respectively

```text
(c_n e_n,b_n e_n,a_n e_n,b_n c_n).                       (11)
```

These are ordinary-primitive on the progression `n=0 mod 5`.
Indeed their invariant `Re(g xy)` values are respectively
`15,-15,-5,-5`, so any ordinary prime dividing both coordinates of
a product must be three or five. An inert prime such as three cannot
divide a product of two ordinary-primitive Gaussian integers. Modulo
five, `M^n=(-1)^n(I+nK)`; on this progression all four vectors equal
their initial vectors up to sign. The initial products in (11) are

```text
11-2i,                 -1+7i,             -1+2i,          13+9i,
```

none of which has ordinary content divisible by five.

Take `P=A_n/2`, `Q_j=P+10(A_n/w_j)` to obtain the corresponding
five-point integer circle tuple. For positive indices in the progression,
all four imaginary parts are positive, distinct, and of order
`lambda^(2n)`. Meanwhile `|A_n|` has order `lambda^(4n)`.
The fourth imaginary-part sequence has initial values `(0,720)` and
the same recurrence with coefficient 322; comparison with the initial
values in Section 2 proves `0<t_1<t_4<t_2<t_3` for every `n>=1`.

The gcd of (11) is `gcd_G(e_n,b_n c_n)` because
`gcd_G(a_n,b_n,c_n)=1`. Its norm is at most 50: the fixed determinants
`det(e_n,b_n)=5` and `det(e_n,c_n)=-10` bound the two pair-gcd norms by
5 and 10, and
`gcd_G(e_n,b_n c_n) | gcd_G(e_n,b_n) gcd_G(e_n,c_n)`.
It follows that

```text
|A_n|/sqrt(50) <= |lcm_G(w_1,w_2,w_3,w_4)| <= |A_n|,
|lcm_G(w_1,w_2,w_3,w_4)| ~ (min_j Im(w_j)/10)^2.          (12)
```

Thus exponent two, if a universal four-cofactor lower bound with
constants depending only on the fixed content can be proved, would be
sharp. This construction supplies no such lower bound. Its arc has
length of order `sqrt(|P|)`, with a fixed limiting constant; it supplies
no unbounded-cardinality example.

A more direct extension of Section 1 fails for a specific algebraic
reason. If four vectors all evolve by one fixed hyperbolic matrix with
eigenvalues `lambda,lambda^(-1)`, the product of any three is a linear
combination of `lambda^(3n),lambda^n,lambda^(-n),lambda^(-3n)`.
Multiplication by a fixed Gaussian number and taking real parts adds
no constant mode. Such a sequence cannot equal a fixed nonzero real
part on unbounded indices: dominance first removes all positive modes,
and the remaining negative modes tend to zero. Consequently the ansatz
`A=g a_1 a_2 a_3 a_4`, `v_j=a_j`, cannot preserve nonzero fixed
cofactor real parts along one common hyperbolic orbit. This only
excludes that ansatz, not arbitrary sharing among four factors.

## 5. Scope relative to the existing corpus

The corpus already has four-point Pell examples on arcs of size
`O(R^(1/3))`, including [the common-unit construction](four_point_bonus_counterexample.md).
The new obstruction recorded here is their relevant mechanism in a
particularly simple vertical-line model: three **equal-real-part**
cofactors, ordinary-primitive complements, and an **exact** lcm rather
than a bound losing all pair gcds. Simultaneous cancellation can be
as large as the example requires even under these constraints.

It would therefore be invalid to infer a quadratic lower bound merely
from three same-content factors or from full-tuple primitivity with
bounded contents. A larger-cardinality statement remains a separate
question; this construction does not settle it.

Run `python3 docs/check_vertical_lcm_pell_counterfamily.py` for exact
integer checks. The identities and growth arguments above give the
infinite-family proof.
