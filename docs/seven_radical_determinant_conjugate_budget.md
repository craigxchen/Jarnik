# All conjugates of the normalized seven-row radical determinants

The middle `3 x 3` radical moment determinant uses every sign conjugate
more efficiently than the scalar radical. After the full finite-place
normalization, its conjugate exponents sum to **exactly zero**. This remains
true under arbitrary higher squareclass dependencies when the seven
radicand squareclasses are pairwise distinct. Thus no uncounted
archimedean gain survives in this determinant's leading exponent.

## 1. Determinants and their complete content normalization

Use the seven-point primitive marked normalization, with positive integral
metric values `c_i`, primitive marked vectors `v_i`, and `Q_i=q(v_i)`.
Fix the signs giving the small physical branch, and put

```text
u_(k-1)(x,y)=(x^(k-1),x^(k-2)y,...,y^(k-1)),
H_k=sum_i epsilon_i sqrt(c_i)
              u_(k-1)(v_i) u_(k-1)(v_i)^T / Q_i^(k-1).
```

The matrix is `k x k`. Its Cauchy--Binet expansion is the sum over
`k`-subsets `I` of terms

```text
w_I = product_(i in I)(epsilon_i sqrt(c_i))
      * product_(i<j in I)d_ij^2 / product_(i in I)Q_i^(k-1).
```

Each `w_I^2` is a positive rational number. Let `g_k` be their positive
rational gcd, defined prime by prime by taking the minimum valuation.
Then every `w_I/sqrt(g_k)` has integral square and is an algebraic integer.
Consequently

```text
Z_k=det(H_k)/sqrt(g_k)
```

is an algebraic integer. This includes denominator primes of all the
terms, rather than normalizing just their clean-prime numerators.

At a clean cut of minority size `s` and depth `e`, for
`a=|I intersect minority|` one has exactly

```text
v_p(w_I^2)=2a(a+7/2-s-k)e,
max(0,k-7+s)<=a<=min(k,s).                         (1)
```

Taking the minimum, and adding the `7,21,35` cut masses, gives

```text
log sqrt(g_1)=0*w+o(w),
log sqrt(g_2)=-(35/2)w+o(w),
log sqrt(g_3)=-63w+o(w).                          (2)
```

The error assertion uses the full marked-profile error bound. For `k=2`
the per-depth minimum is `-1/2` only at triple cuts. For `k=3` it is
`-1/2` at pair cuts and `-3/2` at triple cuts. Their sums are `-35/2`
and `-21/2-105/2=-63` respectively.

## 2. Rank perturbations, including the close-node Vandermonde

The physical branch satisfies

```text
each entry of H_k = O(exp(-77w+o(w))),
sqrt(c_i)=exp(19w+o(w)),
|t_i-t_j|=exp(-16w+o(w)).                          (3)
```

The normalized features in the definition of `H_k` are smooth Veronese
features of the slopes, with coefficients and derivatives of subpower
size in the reduced metric. Changing the signs of a set `J` replaces the
matrix by

```text
H_k - 2 sum_(i in J) epsilon_i sqrt(c_i) f_i f_i^T,
f_i=u_(k-1)(v_i)/Q_i^((k-1)/2).                   (4)
```

This is a perturbation of rank at most `|J|`. In a determinant term using
`m` distinct perturbing vectors, the remaining cofactor has size
`O(exp(-77(k-m)w+o(w)))`. The product of its `m` radical coefficients has
size `exp(19mw+o(w))`. Finally, every `m`-vector exterior minor of the
Veronese features contains their Vandermonde. Its square supplies

```text
exp(-16m(m-1)w+o(w)).
```

Thus that determinant term is at most

```text
exp([-77(k-m)+19m-16m(m-1)]w+o(w)).                (5)
```

Dropping the Vandermonde in this step would lose an essential part of the
budget. For `k<=3`, (5) increases strictly as `m` runs from zero to `k`.
If `r` signs are flipped, the largest permissible value is therefore
`m=min(k,r)`.

Global sign reversal negates `H_k` and preserves `Z_k^2`. Represent each
sign class by its minority set of flipped signs, so `r=0,1,2,3`.
The `64` full sign classes have counts

```text
1, 7, 21, 35.                                    (6)
```

## 3. All 64 normalized exponent bounds

Subtracting (2) from (5) gives the following exponents for `|Z_k|`:

| k | r=0 | r=1 | r=2 | r=3 | Sum with counts (6) |
|---|---:|---:|---:|---:|---:|
| 1 | -77 | 19 | 19 | 19 | 1120 |
| 2 | -273/2 | -81/2 | 47/2 | 47/2 | 896 |
| 3 | -168 | -72 | -8 | 24 | 0 |

For the product of `Z_k^2` over the full sign group, these sums double:
the exponent bounds are `2240`, `1792`, and `0`. In particular the
apparently negative unnormalized `k=3` total is `-4032`, exactly canceled
by subtracting `64*(-63)` from the content normalization.

The `k=3` bounds are sharp in the full fair archimedean regime.
After a real metric normalization and centering the common limiting
slope at zero, its physical small matrix has leading coefficient matrix

```text
C_3 = [[1/16, 0, -1/8],
       [0, -1/8, 0],
       [-1/8, 0, 1/2]].
```

These are the coefficients
`(C_3)_(a,b)=[t^(6-a-b)](1+t^2)^(1/2)` for `a,b=0,1,2`.
The leading complementary coefficients for `r=0,1,2,3` flipped vectors
are respectively

```text
det C_3=-1/512,
det C_3[{1,2},{1,2}]=-1/16,
(C_3)_(2,2)=1/2,
1.
```

They are nonzero. Since `r<=3`, the top perturbation term uses exactly
all the flipped vectors; there is no sum over alternative subsets of
the same size that could cancel it. Its coefficient includes their
nonzero Vandermonde square. Terms with fewer perturbing vectors have a
strictly smaller exponential order. Thus each `k=3` exponent in the
table is an equality up to `o(w)`, not merely a generic upper estimate.

## 4. Actual Galois dependencies

The full `64`-sign product should not automatically be identified with
the field norm when squareclass dependencies are present. This issue can
be resolved exactly for `k=3`.

Let `eta_i=+1` or `-1` be a Galois sign on `sqrt(c_i)`, relative to the
physical branch. The normalized exponent in the last row is

```text
f(eta)=-8 sum_(i<j) eta_i eta_j.                   (7)
```

Indeed, when `r` signs are negative,
`sum_(i<j)eta_i eta_j=21-14r+2r^2`, and (7) is
`-168+112r-16r^2`. It is invariant under global reversal, so either
representative of a sign class gives the same exponent.

Average (7) over the actual Galois group of
`K=Q(sqrt(c_1),...,sqrt(c_7))`. A pair character `eta_i eta_j` has mean
zero unless `c_i c_j` is a rational square, in which case its mean is one.
Consequently

```text
mean_sigma f(eta(sigma))
 = -8 * #{i<j : c_i c_j is a rational square}.     (8)
```

The possible sign of `sqrt(g_3)` disappears on squaring, and `Z_3^2`
belongs to the pair-radical subfield of `K`. Averaging its logarithmic
conjugates over `K` gives the normalized logarithm of its positive
integer norm, even if some conjugates are repeated. Equations (7)--(8)
therefore imply a negative norm exponent if any two radicands share a
squareclass. If their classes are pairwise distinct, the exact leading
norm exponent is zero, regardless of dependencies involving four or more
classes.

In the intended small-arc problem, pairwise squareclass distinction is
also supplied by the separate square-root descent argument. Thus this
degree restriction does not itself resolve the remaining case. On that
case, the middle determinant's complete archimedean and finite-place
budgets are exactly equal at leading order.

The checker `check_positive_radical_norm_budget.py` verifies the content
minima, all three weighted exponent sums, and all 64 Fourier coefficients
of (7). The only nonzero Fourier coefficients are `-8` on pair characters,
which also checks the dependence-free averaging identity (8).
