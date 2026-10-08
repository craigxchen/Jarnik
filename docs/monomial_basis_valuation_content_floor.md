# A radius floor for every integral affine monomial basis change

Prime-by-prime valuation content gives an exact obstruction to descent
by simultaneous monomial changes of the rows. If every prime allocation
in a primitive tuple has only two levels, the source squared radius
divides the primitive squared radius after **every nonsingular integer
affine monomial row change**. This includes arbitrary sequences of
one-row reflections and partial reflections of several rows.

There are actual four-point families with normalized arc constant
tending to zero that satisfy this two-level hypothesis. They are
therefore radius minima in the entire monomial basis orbit. This is
not a counterexample for arbitrarily many short-arc points, and it
does not settle uniformity.

## 1. The invariant and its exact normalization

Let `z_1,...,z_m in Z[i]` be a primitive Gaussian tuple, of common
norm N. Its prime factorization has the form

```text
z_i = epsilon_i product_p pi_p^(e_pi) bar(pi_p)^(s_p-e_pi),
N = product_p p^s_p,            0<=e_pi<=s_p,
min_i e_pi=0,                  max_i e_pi=s_p.
```

Only odd split primes remain: an inert prime or the ramified prime
would otherwise divide every coordinate. For each p define

```text
g_p = gcd_(i,j)(e_pi-e_pj) > 0,
B(z) = product_p p^g_p.                                 (1)
```

Let A be an m-by-m integer matrix with `A 1=1` and define

```text
w_i = product_j z_j^(A_ij).                             (2)
```

Negative exponents are allowed. Every w_i has norm N. Clear one
common Gaussian denominator and divide the resulting whole tuple
by its exact Gaussian gcd. The primitive squared radius N_A has
the formula

```text
N_A = product_p p^(max_i(A e_p)_i-min_i(A e_p)_i).        (3)
```

Indeed the two conjugate valuation vectors before normalization are
`A e_p` and `s_p 1-A e_p`. Clearing the denominator shifts their
minima to nonnegative numbers. Removing the common gcd subtracts
both minima, leaving total norm exponent exactly the displayed
width. Units do not affect this calculation.

Write `e_p=e_p1 1+g_p v_p` with v_p integral. Then

```text
A e_p=e_p1 1+g_p A v_p,
g_p divides every difference of entries of A e_p.       (4)
```

If A is nonsingular, `A e_p` cannot be constant: since `A 1=1`, its
inverse fixes 1, so a constant image would force e_p constant. The
width in (3) is therefore a positive multiple of g_p. Consequently

```text
B(z) divides N_A.                                      (5)
```

If A is unimodular, `A^-1` is also integral and fixes 1. Applying
(4) in both directions proves the stronger equality

```text
gcd_(i,j)((A e_p)_i-(A e_p)_j)=g_p.                     (6)
```

Thus B(z) is an exact invariant of the affine-unimodular monomial
orbit. Formula (5) is valid even for nonsingular integer A with
determinant other than plus or minus one. No assertion is made that
the general lower bound B(z) is attained simultaneously at all primes.

## 2. Two-level allocations are already global radius minima

Suppose every e_p takes only its two extreme values 0 and s_p.
Then `g_p=s_p`, so (1) gives `B(z)=N`. Equation (5) becomes

```text
N divides N_A,                  hence N_A>=N.            (7)
```

In particular, every primitive tuple of squarefree norm N satisfies
(7). Squarefreeness is unnecessary: a prime may have an arbitrary
exponent as long as every row takes either all or none of its chosen
Gaussian orientation.

A one-row reflection `z_l -> z_i z_j/z_l`, with anchors distinct
from l and possibly equal to each other, is such a unimodular
matrix: it fixes the other rows, sends the l-th row to
`e_i+e_j-e_l`, fixes 1, and has determinant -1. Its own inverse is
the same matrix. Simultaneously reflecting any subset of non-anchor
rows has determinant `(-1)^(subset size)` and also fixes 1. Products
of these operations remain in the same group.

Primitive normalization between operations only multiplies every
row by one common Gaussian rational. A later affine row change
preserves that common multiplier because its row sums are one.
Consequently normalization at intermediate steps does not evade
(6) or (7). Intermediate radii may grow and shrink, but no member
of this orbit has primitive norm below the initial two-level norm.

## 3. Actual four-point and five-point examples

The family in
[four_point_bonus_counterexample.md](four_point_bonus_counterexample.md)
has `U+V sqrt(5)=(9+4 sqrt(5))^n`, `n=1 mod 10`, and the blocks

```text
P=-1+2i,
A=(U+V)+2Vi,
B=2V+(U-V)i,
C=(2U+8V)+(3U+V)i.
```

Its four rows are

```text
-P A B C,
-bar(P) A bar(B) bar(C),
-bar(P) bar(A) B bar(C),
-bar(P) bar(A) bar(B) C.
```

The cited note proves that the four block norms are pairwise
coprime and each block is coprime to its conjugate. Every prime is
therefore assigned entirely to one of its two orientations in each
row: all valuation columns have precisely the two levels required
in Section 2, even when the block norms have repeated prime factors.
It also proves primitivity, distinctness, and that the normalized
arc constant tends to zero.

Every member of this infinite family therefore satisfies (7) for
all nonsingular integer affine monomial row changes. In particular,
no reflection, no partial reflection, and no arbitrarily long
sequence of these operations can lower its primitive squared radius.
This is an actual short-arc example, not merely a valuation model;
it has exactly four source points.

The five-point family in
[five_point_affine_shape_cubic_family.md](five_point_affine_shape_cubic_family.md)
also has the two-level property for every permitted parameter
`T=60n`, `10|n`. Its six Gaussian blocks are `1+i b_a n`, with
`b_a in {60,30,20,15,12,10}`. A prime dividing two signed blocks
must divide a nonzero `b_a-b_c` or `b_a+b_c`; their prime factors
belong to `{2,3,5,7,11}`. Inert primes cannot divide a Gaussian
integer of real part one, and `10|n` excludes the primes above two
and five. Opposite signs in a single block have gcd dividing two.
Hence the six norms are pairwise coprime and every block is
conjugate-coprime. Each block has both orientations among the five
rows, proving the two-level assertion without any conjecture about
squarefree values.

These five-point configurations therefore also satisfy (7) for all
parameters, while their normalized arc constant tends to `2sqrt(5)`.
That positive limit is essential: they do not supply five-point
examples whose normalized constants tend to zero.

## 4. Scope for a backwards minimal-counterexample argument

The earlier
[symmetry descent note](symmetry_descent_parallel.md) studied a
single generated orbit point on the original circle and the cost of
adjoining it. The present invariant instead permits simultaneous
replacement of all rows and exact common normalization. It gives
a divisibility floor throughout the resulting basis orbit.

Thus radius minimality alone imposes no new condition on a two-level
configuration: it holds automatically within this entire class of
mutations, regardless of the angles. Actual four-point examples show
that even arbitrarily small normalized arc constants do not remove
the obstruction in that cardinality.

A proposed argument at some larger fixed number of points would
still need to use small angles to exclude the two-level case or to
produce a useful operation outside this class. Singular row maps
can erase a valuation column and are outside the nonsingularity
hypothesis, but they require a separate audit of the number of
distinct output points and the new arc span. Neither a larger
short-arc counterfamily nor a uniform bound is proved here.
