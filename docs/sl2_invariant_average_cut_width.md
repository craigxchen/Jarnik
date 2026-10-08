# An average cut-width bound for binary-form invariants

This note proves a general inequality for the Newton support of an
`SL(2)`-invariant bracket polynomial. It applies to an individual row-weight
component of the eight-point Gale pullback after its common monomial factor
is removed. It does not assert a bound for sums of different row weights,
or a lattice-circle endpoint theorem.

## 1. A polynomial grid lemma

Let `d_1,...,d_n` be nonnegative integers and let

```text
Q = product_i {0,1,...,d_i}.
```

Suppose `p` is a nonzero polynomial over a characteristic-zero field, its
degree in `t_i` is at most `d_i`, and its total degree is at most `k`.
Write `A={a in Q : p(a) != 0}`. Choose a corner `B` by independently taking
`B_i=0` or `d_i`, each with probability one half. Then

```text
E_B min_(a in A) sum_i |B_i-a_i| <= k/2.             (1)
```

Here repeated corners caused by `d_i=0` retain their natural probabilities.

**Proof.** Induct on the number of coordinates, omitting coordinates of
length zero. The assertion for no coordinates is immediate. Put `d=d_n`.
Let `a` and `b` be the first and last indices whose slices contain a point
where `p` is nonzero. Such slices exist because a polynomial with the stated
individual degree bounds is determined by its values on the grid.

Every slice at an integer `j<a` or `j>b` vanishes identically as a polynomial
in the other variables: it vanishes on their entire interpolation grid.
Consequently

```text
p(t) = [product_(0<=j<a or b<j<=d) (t_n-j)] q(t),
r = a+d-b,
deg(q) <= k-r.                                    (2)
```

The polynomials `q(t',a)` and `q(t',b)` are nonzero, obey the remaining
individual degree bounds, and have total degree at most `k-r`. Their
nonzero grid loci agree with the corresponding slices of `p`.

When `B_n=0`, find a nonzero point using slice `a`; when `B_n=d`, use slice
`b`. The induction hypothesis bounds the expected distance in the other
coordinates by `(k-r)/2` in either case. The mean last-coordinate cost is
`(a+d-b)/2=r/2`. Thus the mean distance to the whole nonzero locus is at
most `k/2`. This proves (1).

## 2. Turning invariant coefficients into grid values

Let `g(z_1,...,z_n)` be nonzero, homogeneous of total degree `E`, with
individual degrees at most `d_i`, where

```text
sum_i d_i = 2E,
g(z_1+c,...,z_n+c) = g(z_1,...,z_n).                (3)
```

In particular, (3) holds for the affine evaluation
`g(z)=G((1,z_1),...,(1,z_n))` of an `SL(2)`-invariant polynomial `G`
multihomogeneous of degrees `d_i`.

Polarize each variable `z_i` into `d_i` variables `x_(i,1),...,x_(i,d_i)`.
Explicitly, if `g=sum_alpha c_alpha z^alpha`, set

```text
g_tilde(x) = sum_alpha c_alpha
            product_i e_(alpha_i)(x_(i,1),...,x_(i,d_i))
                       / binomial(d_i,alpha_i).         (4)
```

The resulting polynomial is multilinear, separately symmetric in each
group, homogeneous of total degree `E`, and has `2E` variables.
Polarization is unique among such separately symmetric multilinear
polynomials with the given value on each group's diagonal. Translating
every `x` by the same `c` preserves separate symmetry and multilinearity;
on the group diagonals it preserves the value by (3). Therefore
`g_tilde(x+c)=g_tilde(x)`. In particular,

```text
g_tilde(1-x) = g_tilde(-x) = (-1)^E g_tilde(x).     (5)
```

At a Boolean point with fewer than `E` ones, homogeneity and multilinearity
force `g_tilde` to vanish. At one with more than `E` ones, its complement
has fewer than `E` ones, so (5) forces the same conclusion.

By separate symmetry the Boolean value depends only on the group counts
`t_i`. It is represented on `Q` by the ordinary polynomial

```text
p(t) = sum_alpha c_alpha
       product_i binomial(t_i,alpha_i)/binomial(d_i,alpha_i),  (6)
```

where `binomial(t,a)=t(t-1)...(t-a+1)/a!`. It has total degree at most `E`
and individual degrees at most `d_i`. It vanishes on grid points with
`sum t_i != E`. At a grid point `t` with `sum t_i=E`, only `alpha=t` can
contribute to (6), and hence

```text
p(t) = c_t / product_i binomial(d_i,t_i).           (7)
```

Thus `p` is nonzero and its nonzero grid locus is exactly `Supp(g)`.
There is no claim here that a general polynomial's coefficient support
is its evaluation support; (3)--(5) provide that special property.

## 3. The cut-width inequality

For `S subset {1,...,n}` write

```text
h(S) = max_(alpha in Supp(g)) [alpha(S)-d(S)/2].    (8)
```

The corner `B(S)` with entries `d_i` on `S` and zero off `S` satisfies,
for every support exponent `alpha`,

```text
sum_i |B(S)_i-alpha_i| = d(S)+E-2alpha(S).
```

Taking the minimum over the support gives

```text
dist_1(B(S),Supp(g)) = E-2h(S).                    (9)
```

Apply (1) to (6), using (7). For uniform random `S`,

```text
E_S h(S) >= E/4.                                  (10)
```

For an `SL(2)` invariant, applying the matrix
`[[0,-1],[1,0]]` to its binary-vector arguments gives

```text
g(z) = (-1)^E product_i z_i^(d_i) g(1/z).          (11)
```

Its support is therefore invariant under `alpha -> d-alpha`. This makes
the centered support symmetric and gives `h(S)=h(S^c)`. Also
`h(empty)=h(all)=0`. Consequently (10) is exactly

```text
sum_(unoriented nontrivial cuts S) h(S) >= 2^(n-3) E.  (12)
```

Products of `E` brackets attain equality. Indeed a factor `(z_j-z_i)`
contributes one half to `h(S)` exactly when its endpoints are separated
by the cut. A uniformly random cut separates them with probability one
half. Newton support functions add under multiplication, so the average
is `E/4`, with multiplicities allowed.

## 4. Scope for the eight-point problem

In a homogeneous degree-`D` expression supported in one fixed row weight,
let `r_i` count occurrences of label `i` among its triple indices.
Every coordinate monomial has the same factor `product_i z_i^(2r_i)`.
After this factor is removed, each remaining term is a product of `D`
five-point Vandermondes. It is an `SL(2)` bracket invariant with

```text
d_i=4(D-r_i),     sum_i r_i=3D,
sum_i d_i=20D,    E=10D.                          (13)
```

Its centered exponent `alpha_i-d_i/2` is exactly the centered exponent
of the original pullback after subtracting `2D` in each coordinate.
Therefore every nonzero such fixed-row-weight pullback satisfies

```text
W >= 2^5 * 10D = 320D,
sum_S nu_S = 373D-W <= 53D.                       (14)
```

This permits arbitrary cancellations and any number of terms within the
fixed row weight, in every degree. It does not apply directly to a sum
of different row weights: the monomial factor and the invariant
multidegrees in (13) then vary between components.
