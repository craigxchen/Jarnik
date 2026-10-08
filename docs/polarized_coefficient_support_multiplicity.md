# Multiplicity and average corner distance of coefficient support

This note proves the multiplicity extension of the invariant cut-width
bound. Its input is an ordinary polynomial and its output is a bound on
the Newton support after all cancellations. It uses the polynomial grid
lemma proved in `sl2_invariant_average_cut_width.md`, Section 1.

## 1. Coefficient-indicator polynomial

Let `f(z_1,...,z_n)` be nonzero over a characteristic-zero field, with
individual degrees at most `d_i`. Put `N=sum_i d_i`. Polarize each variable
into a group of `d_i` variables to obtain the separately symmetric,
multilinear polynomial

```text
F(x) = sum_(A subset {1,...,N}) C_A x_A.
```

For group counts `alpha_i=|A intersect group i|`, the coefficient is

```text
C_A = c_alpha / product_i binomial(d_i,alpha_i),
f(z) = sum_alpha c_alpha z^alpha.                  (1)
```

Suppose `f` vanishes at `(1,...,1)` to order at least `m`. Then `F` does
also. Indeed `F(1+y)` is separately symmetric and multilinear, and its
restriction to group diagonals is `f(1+z)`. Diagonal restriction is
injective on the separately symmetric multilinear polynomials, degree by
degree, since the basis `product_i e_(alpha_i)(y_i)` restricts to distinct
monomials times nonzero binomial coefficients. Thus all homogeneous parts
of total degree less than `m` vanish in `F(1+y)`.

Form another multilinear polynomial

```text
Q(x) = sum_A (-1)^|A| C_A x_A product_(j notin A)(1-x_j).  (2)
```

At a Boolean vertex `1_A`, its value is `(-1)^|A| C_A`. In particular its
Boolean nonzero locus is the coefficient support of `F`. The coefficient
of the monomial `x_T` in (2) is

```text
[x_T] Q = (-1)^|T| sum_(A subset T) C_A
        = (-1)^|T| F(1_T).                         (3)
```

If `|T|>N-m`, then the point `1_T` differs from the all-ones point in fewer
than `m` coordinates. Evaluating the multilinear Taylor expansion of `F`
about all-ones at this point gives zero: every nonzero Taylor monomial
uses at least `m` distinct coordinates. Therefore (3) vanishes, proving

```text
deg Q <= N-m.                                     (4)
```

## 2. Compression to the grid

The polynomial `Q` is separately symmetric in the same groups. Express
it in the basis `product_i e_(alpha_i)(x_i)` and replace each elementary
symmetric polynomial by `binomial(t_i,alpha_i)`. This gives an ordinary
polynomial `q(t)` with individual degrees at most `d_i` and total degree
at most `N-m`. Its value at an integer grid point

```text
t in product_i {0,...,d_i}
```

is the value of `Q` at any Boolean vertex with those group counts.
By (1)--(2),

```text
q(t) = (-1)^|t| c_t / product_i binomial(d_i,t_i).  (5)
```

Hence its nonzero grid locus is exactly `Supp(f)`. Applying the polynomial
grid lemma gives

```text
E_B dist_1(B,Supp(f)) <= (N-m)/2,                  (6)
```

where the corner coordinates `B_i` are independently `0` or `d_i` with
equal probabilities. This conclusion does not require homogeneity or
reciprocity.

## 3. Homogeneous reciprocal support

Now suppose `f` is homogeneous of total degree `E`, `N=2E`, and write

```text
h(S) = max_(alpha in Supp(f)) [alpha(S)-d(S)/2].
```

For the corner with coordinates `d_i` on `S` and zero off `S`,

```text
dist_1(B(S),Supp(f)) = E-2h(S).                    (7)
```

Combining (6) and (7) yields

```text
E_S h(S) >= m/4.                                 (8)
```

If the support is invariant under `alpha -> d-alpha`, then
`h(S)=h(S^c)` and `h(empty)=h(all)=0`. Thus

```text
sum_(unoriented nontrivial cuts S) h(S) >= 2^(n-3) m.  (9)
```

No assertion that its Newton polytope contains a bracket-product polytope
is used or needed.

## 4. Application to the eight-point Gale pullback

For a nonzero homogeneous degree-`D` coordinate polynomial pullback
`f=F(R_I(z))`, the established formulas give

```text
deg f=16D,    d_i=4D,    sum_i d_i=32D=2deg(f),
Supp(f) invariant under alpha -> 4D*1-alpha,
mult_(1,...,1) f >= 10D.
```

The last inequality holds because every `R_I` is a monomial times a
five-point Vandermonde, which has ten linear difference factors at the
diagonal. Products of `D` coordinates vanish to order at least `10D`,
and their sums retain that lower bound.

Equation (9), with `n=8`, now proves

```text
W(f) >= 32 mult_1(f) >= 320D,
sum_S nu_S(F) = 373D-W(f) <= 53D.                 (10)
```

This applies to all nonzero pullbacks, including arbitrary cancellation
between different row weights. It rules out a strict gain in this stated
generic clean-cut weighted-order budget. It does not by itself prove a
lattice-circle endpoint bound or exclude a different arithmetic mechanism.
